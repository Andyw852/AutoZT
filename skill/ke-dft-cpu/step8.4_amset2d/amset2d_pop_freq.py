#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""amset2d_pop_freq.py —— 取二维体系 Γ 点**面内极性**声子模的有效频率（THz）。

为什么不用 amset 自带的 phonon-frequency：
    AMSET 的 effective_phonon_frequency 对**三维**球面方向做加权平均，在 slab 里
    面外极性模（ZO/A1）也会进权重；二维 Fröhlich 耦合只用面内极性模（E' 之类）。
    实测 CrS2 单层：AMEET 口径给出 12.45 THz（含 2.4% 的 15.07 THz 面外模），
    面内口径给出 12.38 THz（纯 E' 简并模）。

做法（与 AMSET 的 get_phonon_weight 同公式，只把 q 方向限制在 xy 面内）：
    w_n ∝ |Σ_atoms q̂·Z*_a·e_a(n) / √(M_a ω_n)|，q̂ 取面内若干方向做平均；
    去掉声学模与权重可忽略的模后，按权重加权平均频率。
    pop_frequency（二维）取该值：插件再乘 √(ε̃0(q)/ε̃∞(q))（二维 LST）得到 ω_LO(q)。

用法：python amset2d_pop_freq.py [OUTCAR] [vasprun.xml]
      需要 amset 运行环境（import amset.tools.phonon_frequency 的同款解析链）。
输出：一行 JSON 到 stdout（gen 脚本解析）。
"""
import json
import sys

import numpy as np

# [SKILL_REV] 与 gen_step14_amset2d.py 同批版本
_SKILL_REV = "2026-09-16-2d-inplane-pop"
_N_DIRS = 72          # 面内 q 方向数（0..2π 均布）
_KEEP_RATIO = 1e-3    # 权重小于 max(w)*该比例的模视为非极性，剔除


def inplane_directions(n=_N_DIRS):
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([np.cos(ang), np.sin(ang), np.zeros(n)], axis=1) * 0.01


def spheresph_directions(n=_N_DIRS):
    """三维参考口径：球面上 n 个方向（与 AMSET 的 DEFAULT_QUAD 不同，仅作对照）。"""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1) * 0.01


def weights_for(directions, freqs, eigvecs, born, structure):
    """每个模在给定 q 方向集合上的平均耦合权重（与 AMSET 的公式一致）。"""
    out = []
    for ev, f in zip(eigvecs, freqs):
        if f <= 0:                     # 声学模/虚频：权重置 0
            out.append(0.0)
            continue
        zes = []
        for ab, av, site in zip(born, ev, structure):
            zes.append(np.dot(ab, av) / np.sqrt(site.specie.atomic_mass * f))
        zes = np.array(zes)
        all_w = np.abs(np.einsum("ai,di->d", zes, directions))
        val = float(np.mean(all_w))
        out.append(0.0 if not np.isfinite(val) else val)
    return np.array(out)


def effective(freqs, weights):
    """按权重加权平均；剔除声学模（f<=0）与权重可忽略的模后归一化。"""
    keep = (weights > 0) & (freqs > 0) & (weights >= _KEEP_RATIO * weights.max())
    if not keep.any():
        return None, keep
    w = weights[keep] / weights[keep].sum()
    return float(np.sum(w * freqs[keep])), keep


def main():
    outcar_path = sys.argv[1] if len(sys.argv) > 1 else "OUTCAR"
    vasprun_path = sys.argv[2] if len(sys.argv) > 2 else "vasprun.xml"

    from pymatgen.io.vasp import Outcar, Vasprun

    vasprun = Vasprun(vasprun_path)
    outcar = Outcar(outcar_path)
    outcar.read_lepsilon()
    born = np.array(outcar.born)

    eigenvalues = -vasprun.normalmode_eigenvals[::-1]
    eigenvectors = vasprun.normalmode_eigenvecs[::-1]
    freqs = np.sqrt(np.abs(eigenvalues)) * np.sign(eigenvalues)
    if int(vasprun.vasp_version.split(".")[0]) < 6:
        freqs = freqs * 15.633302
    structure = vasprun.final_structure

    w2 = weights_for(inplane_directions(), freqs, eigenvectors, born, structure)
    w3 = weights_for(spheresph_directions(), freqs, eigenvectors, born, structure)
    f2, keep2 = effective(freqs, w2)
    f3, _ = effective(freqs, w3)

    res = {
        "skill_rev": _SKILL_REV,
        "n_modes": int(len(freqs)),
        "frequencies_THz": [round(float(x), 4) for x in freqs],
        "weights_inplane": [float("%.3e" % x) for x in w2],
        "weights_sphere3d": [float("%.3e" % x) for x in w3],
        "polar_inplane_modes": [int(i) for i in np.where(keep2)[0]],
        "pop_frequency_THz": None if f2 is None else round(f2, 4),
        "pop_frequency_3d_ref_THz": None if f3 is None else round(f3, 4),
        "note": "pop_frequency_THz = Γ 点面内极性模权重平均（二维口径）",
    }
    print(json.dumps(res, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
