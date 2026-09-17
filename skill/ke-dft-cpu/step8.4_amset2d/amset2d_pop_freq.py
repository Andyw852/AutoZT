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


# ---- Δr 交叉核对（2026-09-16 新增）--------------------------------------
# 路线①（介电）：Δr = c(ε0,∥ − ε∞,∥)/2     —— gen_step14 已在用
# 路线②（声子模求和）：Δr_αβ = (2π/A)·Σ_ν X_να X_νβ / ω_ν²     —— 本函数
#     X_να = Σ_κ Z*_κ,αβ' e_νκβ' / √M_κ        （原子单位：e=1、m_e=1、ħ=1）
#     · e 用 **质量加权归一** 的本征矢（vasprun 的 normalmode_eigenvecs，
#       Σ_κ|e_κ|²=1，即 VASP "after division by SQRT(m)" 那一组）；
#       若改用 OUTCAR 里未除质量的那组，必须自己除 √M，别混用。
#     · M 用**电子质量**（amu × 1822.888）—— 路线①③四位有效数字全对的前提；
#       用 amu 会差 M 倍（实测 1125 Å vs 0.617 Å）。
#     · 声学模必须剔除（|ω|→0，1/ω² 把数值噪声放大）。
# 实测 MoS2 单层：路线① 0.6173 Å、路线② 0.617319 Å（5 位有效数字一致）。
_BOHR = 0.529177210903          # Å/bohr
_THZ_TO_HA = 4.135667696e-3 / 27.211386245988   # THz -> Hartree
_AMU_TO_ME = 1822.888486209     # amu -> m_e
_ACOUSTIC_THZ = 1.0             # 频率绝对值小于此值视为声学模，剔除


def _inplane_area_bohr(structure):
    lat = np.array(structure.lattice.matrix)
    return float(np.linalg.norm(np.cross(lat[0], lat[1]))) / _BOHR ** 2


def delta_r_mode_sum(freqs_thz, eigvecs, born, structure):
    """路线②：Γ 点声子模求和给出 Δr 张量（Å，2×2）。声学模剔除。"""
    A = _inplane_area_bohr(structure)
    m_me = np.array([s.specie.atomic_mass * _AMU_TO_ME for s in structure])
    ev = np.asarray(eigvecs, float)
    w = np.asarray(freqs_thz, float)
    x = np.einsum("kab,nkb->na", np.asarray(born, float),
                  ev / np.sqrt(m_me)[None, :, None])          # (nmode, 3)
    keep = np.abs(w) > _ACOUSTIC_THZ
    w_ha = w[keep] * _THZ_TO_HA
    xm = x[keep][:, :2]
    tensor = (2 * np.pi / A) * np.einsum("na,nb,n->ab", xm, xm, 1.0 / w_ha ** 2)
    return tensor * _BOHR, int(keep.sum()), int(len(w) - keep.sum())


def delta_r_dielectric(outcar, structure):
    """路线①：c(ε0,∥ − ε∞,∥)/2 与 r∞ = c(ε∞,∥ − 1)/2（Å，面内对角平均）。"""
    lat = np.array(structure.lattice.matrix)
    c_len = float(abs(np.dot(lat[2], np.cross(lat[0], lat[1])))
                  / np.linalg.norm(np.cross(lat[0], lat[1])))
    eps_inf = np.array(outcar.dielectric_tensor, float)
    eps_ion = np.array(outcar.dielectric_ionic_tensor, float)
    d_inf = (eps_inf[0, 0] + eps_inf[1, 1]) / 2.0
    d_ion = (eps_ion[0, 0] + eps_ion[1, 1]) / 2.0
    return {"c_A": round(c_len, 4),
            "eps_inf_inplane": round(float(d_inf), 5),
            "eps_ionic_inplane": round(float(d_ion), 5),
            "r_inf_A": round(c_len * (d_inf - 1.0) / 2.0, 4),
            "delta_r_A": round(c_len * d_ion / 2.0, 6)}


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
    # ---- Δr 两种算法交叉核对 ----
    try:
        dr2, n_used, n_ac = delta_r_mode_sum(freqs, eigenvectors, born, structure)
        dr2_iso = float((dr2[0, 0] + dr2[1, 1]) / 2.0)
        dr1 = delta_r_dielectric(outcar, structure)
        res["delta_r_check"] = {
            "route1_dielectric_A": dr1["delta_r_A"],
            "route2_mode_sum_A": round(dr2_iso, 6),
            "route2_tensor_A": np.round(dr2, 6).tolist(),
            "deviation_pct": (round(100.0 * (dr2_iso - dr1["delta_r_A"])
                                    / dr1["delta_r_A"], 3)
                              if dr1["delta_r_A"] else None),
            "modes_used": n_used, "acoustic_removed": n_ac,
            "r_inf_A": dr1["r_inf_A"],
            "eps_inf_inplane": dr1["eps_inf_inplane"],
            "eps_ionic_inplane": dr1["eps_ionic_inplane"],
            "anchors_literature": {"MoS2_delta_r_A": 0.53, "MoS2_r_inf_A": 46.5},
            "note": ("路线①= c(ε0-ε∞)/2（DFPT 介电）；路线②= (2π/A)ΣX²/ω²（声子模求和，"
                     "原子单位、M 用电子质量、e 用质量加权归一的本征矢）。两条应一致；"
                     "差几个百分点以上说明本征矢口径或声学模剔除有问题。"),
        }
        if dr1["delta_r_A"] and abs(res["delta_r_check"]["deviation_pct"]) > 5.0:
            print("[amset2d][WARN] Δr 两条算法差 %.2f%%（>5%%）：检查本征矢口径"
                  "（是否混用了除/未除质量的版本）与声学模剔除"
                  % res["delta_r_check"]["deviation_pct"], file=sys.stderr)
    except Exception as _e:                                    # noqa: BLE001
        res["delta_r_check"] = {"error": "%s: %s" % (type(_e).__name__, _e)}
        print("[amset2d][WARN] Δr 交叉核对失败：%s" % _e, file=sys.stderr)

    print(json.dumps(res, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
