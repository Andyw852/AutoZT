#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kz_overlap_probe.py —— 2D：AMSET 沿 k_z 插值波函数系数造成的重叠假象有多大（V126，秒级，不提交作业）。

背景（V125 的 MoS₂ 实测）：KZ_CAP_2D 把插值网格的 k_z 从 21 层截到 3 层后，n 型 ADP 迁移率 −5.7%、
p 型 +0.5%，而合成体系（重叠恒为 1）只有 ≤0.1%。差别只能来自重叠。AMSET 求重叠时，对不在 h5 网格上的
k 点把复系数**线性插值**（实部虚部分开）再归一化；h5 沿 k_z 只有 3 层（0, ±1/3），旧跑法 21 层里
大多数层落在两层之间。VASP 每个 k 点的整体相位是任意的，两层系数又不严格平行（slab 形状因子），
插出来的"中间态"与真实态的重叠会被压低，压低多少取决于随机相位 —— 散射率偏低，迁移率偏高。
k_z 截到 3 层时，插值层正好落在 h5 的网格层上，不需要沿 k_z 插值。

本工具用 AMSET 自己的 WavefunctionOverlapCalculator（与作业同一条路径，含 desym 补丁），在带边附近
取若干面内 k 点，比较两种 k_z 层集合上的平均重叠：
    旧：k_z = j/N_old（默认 21 层）；新：k_z = j/3。
2D 核只用 q∥，所以散射率 ∝ 同一面内位置上各 k_z 层之间重叠的平均；迁移率按 τ = 1/率 平均：
    μ_新/μ_旧（预测）= mean_i[1/mean_j M_new] / mean_i[1/mean_j M_old]
真正单层的极限是 q_z = 0（只有一层，重叠 = 1）；新集合的平均是 (1+2F)/3，F = 网格层之间的真实重叠。

用法（在 S8.4 运行目录）：
    python <skill>/ke-dft-cpu/tools/kz_overlap_probe.py [--dir .] [--nkz-old 21] [--emax 0.15] [--json out.json]
能带窗口从 amset.log 的 "Interpolating spin-up bands L-H" 读（或 --bands L:H，1-based）。
"""
import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np


def band_window(run_dir, explicit=None):
    if explicit:
        lo, hi = (int(x) for x in explicit.replace("-", ":").split(":"))
        return lo, hi
    p = Path(run_dir) / "amset.log"
    if p.is_file():
        m = re.search(r"Interpolating spin-(?:up|down) bands\s+(\d+)\s*[-—–~]\s*(\d+)",
                      p.read_text(errors="ignore"))
        if m:
            return int(m.group(1)), int(m.group(2))
    raise ValueError("读不到能带窗口：%s 里没有 'Interpolating ... bands L-H'，请用 --bands L:H" % p)


def load_calculator(h5, run_dir=None, desym_fix=True):
    """与作业同一路径：挂 desym 补丁（运行目录里有插件时）后 WavefunctionOverlapCalculator.from_file。"""
    if not hasattr(np, "string_"):
        np.string_ = np.bytes_
    if desym_fix and run_dir and (Path(run_dir) / "amset_desym_fix.py").is_file():
        os.environ.setdefault("AZ_DESYM_FIX", "1")
        sys.path.insert(0, str(Path(run_dir).resolve()))
        import amset_desym_fix  # noqa: F401
    from amset.interpolation.wavefunction import WavefunctionOverlapCalculator
    return WavefunctionOverlapCalculator.from_file(str(h5))


def kz_set(n):
    kz = np.arange(n) / float(n)
    return kz - np.rint(kz)


def pair_matrix(calc, spin, band, k_inplane, kzs, axis=2):
    """M[i, j] = |<ψ_band(k_i)|ψ_band(k_j)>|²，k_i 与 k_inplane 面内相同、k_z = kzs[i]（AMSET get_overlap）。"""
    ks = np.repeat(np.asarray(k_inplane, float)[None], len(kzs), axis=0)
    ks[:, axis] = kzs
    n = len(kzs)
    M = np.zeros((n, n))
    for i in range(n):
        M[i] = calc.get_overlap(spin, int(band), ks[i], np.full(n, int(band)), ks)
    return M


def summarize(M):
    row = M.mean(axis=1)
    return {"mean": float(M.mean()), "min": float(M.min()), "inv_tau_weight": float(np.mean(1.0 / row))}


def edge_samples(energies, kpoints, band_lo, band_hi, efermi, emax=0.15, nsample=12, axis=2, degen_tol=0.005):
    """带边附近、k_z = 0 的 (载流子, 自旋, h5 带下标, k, 离带边能量) 样本。
    energies: {spin: (nb_abs, nk)} eV（vasprun 全部带）；band_lo/band_hi：h5 窗口（1-based，含两端）；
    efermi：vasprun 的费米能（半导体在价带顶附近/带隙里）。"""
    kp = np.asarray(kpoints, float)
    on_plane = np.abs(kp[:, axis] - np.rint(kp[:, axis])) < 1e-6
    allE = np.concatenate([np.asarray(e) for e in energies.values()], axis=1)
    nb = allE.shape[0]
    below = [b for b in range(nb) if allE[b].max() <= efermi + 0.05]
    if not below or below[-1] + 1 >= nb or allE[below[-1] + 1].min() <= allE[below[-1]].max():
        raise ValueError("按费米能 %.3f eV 找不到带隙（金属？）" % efermi)
    vb = below[-1]
    cb = vb + 1
    vbm, cbm = allE[vb].max(), allE[cb].min()
    out = []
    for car, bands, ref, sgn in (("electron", range(cb, nb), cbm, 1), ("hole", range(vb, -1, -1), vbm, -1)):
        cand = []
        for spin, E in energies.items():
            E = np.asarray(E)
            for b in bands:
                if not band_lo - 1 <= b <= band_hi - 1:
                    continue
                d = sgn * (E[b] - ref)
                for k in np.where(on_plane & (d >= -1e-6) & (d <= emax))[0]:
                    others = np.delete(E[:, k], b)
                    if np.min(np.abs(others - E[b, k])) < degen_tol:
                        continue
                    cand.append((float(d[k]), spin, b - (band_lo - 1), kp[k].copy()))
        cand.sort(key=lambda x: x[0])
        if len(cand) > nsample:
            idx = np.linspace(0, len(cand) - 1, nsample).round().astype(int)
            cand = [cand[i] for i in idx]
        out += [(car, s, b, k, e) for e, s, b, k in cand]
    return out, float(vbm), float(cbm)


def probe(calc, samples, nkz_old=21, nkz_new=3, axis=2):
    old_k, new_k = kz_set(nkz_old), kz_set(nkz_new)
    rep = {}
    for car in ("electron", "hole"):
        rows = []
        for c, spin, b, k, e in samples:
            if c != car:
                continue
            so = summarize(pair_matrix(calc, spin, b, k, old_k, axis))
            Mn = pair_matrix(calc, spin, b, k, new_k, axis)
            sn = summarize(Mn)
            rows.append({"E_from_edge_eV": round(e, 4), "k": [round(float(x), 5) for x in k], "band_h5": int(b),
                         "old": so, "new": sn, "F_grid": float(Mn[0, 1:].mean()) if nkz_new > 1 else 1.0})
        if not rows:
            continue
        mo = np.mean([r["old"]["mean"] for r in rows])
        mn = np.mean([r["new"]["mean"] for r in rows])
        wo = np.mean([r["old"]["inv_tau_weight"] for r in rows])
        wn = np.mean([r["new"]["inv_tau_weight"] for r in rows])
        rep[car] = {"n": len(rows), "mean_overlap_old": round(float(mo), 5), "mean_overlap_new": round(float(mn), 5),
                    "min_overlap_old": round(float(min(r["old"]["min"] for r in rows)), 5),
                    "F_grid_mean": round(float(np.mean([r["F_grid"] for r in rows])), 5),
                    "mu_new_over_old_pred": round(float(wn / wo), 4), "rows": rows}
    return rep


def fermi_from_occupations(eigenvalues, occ_tol=0.5):
    """占据数定带隙中点：VBM = 占据态（占据数 > occ_tol）最高能，CBM = 空态最低能。
    eigenvalues: {spin: (nk, nb, 2)}（vasprun 的 [能量, 占据数]）。没有带隙 -> None。"""
    occ_e, emp_e = [], []
    for v in eigenvalues.values():
        v = np.asarray(v)
        occ = v[..., 1] > occ_tol
        occ_e.append(v[..., 0][occ])
        emp_e.append(v[..., 0][~occ])
    occ_e, emp_e = np.concatenate(occ_e), np.concatenate(emp_e)
    if not len(occ_e) or not len(emp_e):
        return None
    vbm, cbm = float(occ_e.max()), float(emp_e.min())
    return 0.5 * (vbm + cbm) if cbm > vbm else None


def _energies_from_vasprun(path):
    from pymatgen.io.vasp.outputs import Vasprun
    vr = Vasprun(str(path), parse_dos=False, parse_potcar_file=False)
    E = {s: np.asarray(v)[:, :, 0].T for s, v in vr.eigenvalues.items()}   # (nb, nk)
    # parse_dos=False 时 pymatgen 不读 <dos> 里的 efermi（一直是 None）；有的 vasprun 干脆没有这个标签
    ef = vr.efermi if vr.efermi is not None else fermi_from_occupations(vr.eigenvalues)
    if ef is None:
        raise ValueError("%s：没有 efermi，按占据数也找不到带隙（金属？）" % path)
    return E, np.asarray(vr.actual_kpoints, float), float(ef)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=".")
    ap.add_argument("--h5", default=None)
    ap.add_argument("--vasprun", default=None)
    ap.add_argument("--bands", default=None, help="h5 的能带窗口 L:H（1-based）；默认读 amset.log")
    ap.add_argument("--nkz-old", type=int, default=21)
    ap.add_argument("--nkz-new", type=int, default=3)
    ap.add_argument("--emax", type=float, default=0.15, help="带边以内多少 eV 取样")
    ap.add_argument("--nsample", type=int, default=12)
    ap.add_argument("--no-desym-fix", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    d = Path(a.dir)
    try:
        lo, hi = band_window(d, a.bands)
        E, kp, ef = _energies_from_vasprun(a.vasprun or d / "vasprun.xml")
        samples, vbm, cbm = edge_samples(E, kp, lo, hi, ef, a.emax, a.nsample)
        calc = load_calculator(a.h5 or d / "wavefunction.h5", d, not a.no_desym_fix)
    except (ValueError, OSError, KeyError) as e:
        print("[ERROR] %s" % e, file=sys.stderr)
        return 2
    rep = probe(calc, samples, a.nkz_old, a.nkz_new)
    print("能带窗口 %d-%d，VBM %.3f / CBM %.3f eV；k_z 层：旧 %d，新 %d" % (lo, hi, vbm, cbm, a.nkz_old, a.nkz_new))
    print("%-9s %4s %12s %12s %12s %10s %18s" % ("载流子", "样本", "旧平均重叠", "新平均重叠", "旧最小重叠",
                                              "网格层F", "预测 μ新/μ旧"))
    for car, r in rep.items():
        print("%-9s %4d %12.4f %12.4f %12.4f %10.4f %17.3f"
              % (car, r["n"], r["mean_overlap_old"], r["mean_overlap_new"], r["min_overlap_old"],
                 r["F_grid_mean"], r["mu_new_over_old_pred"]))
    print("解读：网格层 F 接近 1 而旧平均明显更低 -> 旧 %d 层结果含 k_z 插值假象（散射率偏低、迁移率偏高），"
          "截到 %d 层更接近单层物理（q_z=0 时重叠 = 1）。预测比与 compare_transport_json 的实测差对得上，"
          "就说明 KZ_CAP_2D 的差别来自这里。" % (a.nkz_old, a.nkz_new))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            fh.write(json.dumps({"band_window": [lo, hi], "vbm": vbm, "cbm": cbm, "probe": rep},
                                ensure_ascii=False, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
