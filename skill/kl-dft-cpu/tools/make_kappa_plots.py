#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kl-dft-cpu 结果 -> κL 图 + 表（对【已有结果】补生成，不必重算）。

★ step6_kappa 只产出 kappa_summary.json / kappa-m*.hdf5，技能内【没有】κL 绘图的 gen 步骤；
  本脚本补齐这部分，并把 S5.1_plot 已算好的声子谱图一并归集。

用法:
    python make_kappa_plots.py <kl-dft-cpu 目录> [<输出目录>] [<材料名>]

<kl-dft-cpu 目录>  含 step6_kappa/kappa_summary.json 的目录（本地 result/ 或集群工作目录均可）。
                   也接受材料根目录（自动找 ./*/step6_kappa 或 ./step6_kappa）。
<输出目录>          缺省 = <材料根>/kl_plots（存在同名目录时先清空旧图，保留 README 以外的旧产物覆盖写）。
<材料名>            缺省 = 由 kappa_summary 的 path 推断，失败则用目录名。

产物:
    <M>_kl_vs_T.png              κL(T) 面内 + 2D 归一化（双轴）
    <M>_kl_vs_T.csv              温度依赖长表
    <M>_kl_mesh_convergence.png  q 网格收敛 + κ∞ 外推
    <M>_kl_cutoff_selection.png  FC3 截断扫描（含 plateau_ok 标记）
    <M>_kl_anisotropy.png        xx/yy/zz 分量
    <M>_kl_summary.csv/.json     300 K 标量汇总
    phonon_band/                 S5.1_plot 声子谱图（若存在则一并复制）
"""
from __future__ import annotations
import csv
import json
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _find_summary(root: Path) -> Path:
    """在 root 或 root/*/ 下找 step6_kappa/kappa_summary.json。"""
    cands = [root / "step6_kappa" / "kappa_summary.json"]
    cands += sorted(root.glob("*/step6_kappa/kappa_summary.json"))
    for c in cands:
        if c.is_file():
            return c
    sys.exit(f"[ERROR] 在 {root} 下找不到 step6_kappa/kappa_summary.json")


def _find_plot_dir(root: Path) -> Path | None:
    """找 S5.1_plot 产物目录（step5_phonon_plot）。"""
    for c in [root / "step5_phonon_plot"] + sorted(root.glob("*/step5_phonon_plot")):
        if c.is_dir() and any(c.glob("*.png")):
            return c
    return None


def _arr1(d, key, idx=0):
    """取 [[a,b],[c,d]] 的第 idx 列，或直接返回一维。"""
    v = d[key]
    out = []
    for x in v:
        out.append(x[idx] if isinstance(x, (list, tuple)) else x)
    return np.asarray(out, dtype=float)


def fig_kappa_vs_T(T, k_ii, k_n, mat, out: Path):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(T, k_ii, "o-", color="#1f77b4", lw=2, ms=5,
            label=r"$\kappa_{L}$ in-plane (W/mK)")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$\kappa_{L}$ (W m$^{-1}$ K$^{-1}$)", color="#1f77b4")
    ax.tick_params(axis="y", labelcolor="#1f77b4")
    ax.grid(alpha=.3)
    ax2 = ax.twinx()
    ax2.plot(T, k_n, "s--", color="#d62728", lw=1.6, ms=4,
             label=r"2D-normalized $\kappa_{L}\cdot h/d$")
    ax2.set_ylabel(r"2D-normalized $\kappa_{L}$ (W/mK)", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    if len(T):
        i300 = int(np.argmin(np.abs(T - 300)))
        ax.annotate(f"{k_ii[i300]:.2f}\n@300K", (T[i300], k_ii[i300]),
                    textcoords="offset points", xytext=(8, 10),
                    fontsize=8, color="#1f77b4")
        ax2.annotate(f"{k_n[i300]:.2f}\n@300K", (T[i300], k_n[i300]),
                     textcoords="offset points", xytext=(8, -18),
                     fontsize=8, color="#d62728")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], loc="upper right",
              fontsize=9, framealpha=.9)
    ax.set_title(f"{mat}: lattice thermal conductivity (kl-dft-cpu, RTA)", fontsize=11.5)
    fig.tight_layout()
    p = out / f"{mat}_kl_vs_T.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    print("[fig]", p)


def fig_mesh(d, mat, out: Path):
    mc = (d.get("mesh_convergence") or {}).get(d.get("bte_method") or "rta", {})
    meshes = mc.get("meshes") or []
    k300 = mc.get("kappa300_inplane") or []
    if not meshes or not k300:
        print("[skip] 无网格收敛数据")
        return
    xs = [int(str(m).split()[0]) for m in meshes]
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(xs, k300, "o-", lw=2, ms=6, color="#2ca02c")
    kinf = mc.get("kappa_inf_extrap_300K_inplane")
    if kinf:
        ax.axhline(kinf, ls="--", lw=1.2, color="#888",
                   label=rf"extrap $\kappa_\infty$={kinf:.2f}")
    for x, y in zip(xs, k300):
        ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=8)
    ax.set_xlabel("q-mesh (N x N x 1)")
    ax.set_ylabel(r"$\kappa_{L}$(300K) in-plane (W/mK)")
    ax.set_title(f"{mat}: BTE mesh convergence "
                 f"(max rel. change {mc.get('max_rel_change_pct', 0):.2f}%)", fontsize=11)
    ax.grid(alpha=.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    p = out / f"{mat}_kl_mesh_convergence.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    print("[fig]", p)


def fig_cutoff(d, mat, out: Path):
    cs = d.get("cutoff_selection") or {}
    recs = cs.get("records") or []
    cuts, kaps, errs, oks = [], [], [], []
    for r in recs:
        c, k = r.get("cut"), r.get("kappa")
        if c is None or k is None:
            continue
        cuts.append(float(c)); kaps.append(float(k))
        errs.append(float(r.get("kappa_err") or 0.0))
        oks.append(bool(r.get("plateau_ok")))
    if not cuts:
        print("[skip] 无 cut3 记录")
        return
    order = np.argsort(cuts)
    cuts = np.array(cuts)[order]; kaps = np.array(kaps)[order]
    errs = np.array(errs)[order]; oks = [oks[i] for i in order]
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.errorbar(cuts, kaps, yerr=errs, fmt="o-", lw=2, ms=6, capsize=3,
                color="#9467bd", ecolor="#bbbbbb")
    for c, k, ok in zip(cuts, kaps, oks):
        if ok:
            ax.plot([c], [k], "o", ms=13, mfc="none", mec="#2ca02c", mew=1.6)
    chosen = d.get("chosen_cutoff_A") or cs.get("chosen_cut")
    if chosen:
        ax.axvline(chosen, ls="--", color="#d62728", lw=1.2,
                   label=f"chosen cut3 = {chosen} A")
    ax.plot([], [], "o", ms=11, mfc="none", mec="#2ca02c", label="plateau_ok")
    ax.set_xlabel("FC3 pair cutoff (A)")
    ax.set_ylabel(r"$\kappa_{L}$(300K) in-plane (W/mK)")
    ax.set_title(f"{mat}: 3rd-order FC cutoff scan (status: {cs.get('status')})",
                 fontsize=11)
    ax.grid(alpha=.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    p = out / f"{mat}_kl_cutoff_selection.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    print("[fig]", p)


def fig_aniso(d, T, mat, out: Path):
    k_xyz = np.array(d.get("kappa_xx_yy_zz") or [], dtype=float)
    if k_xyz.ndim != 2 or k_xyz.shape[1] < 3 or not len(T):
        print("[skip] 无 xx/yy/zz 分量")
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(T, k_xyz[:, 0], "o-", lw=2, ms=5, label=r"$\kappa_{xx}$")
    ax.plot(T, k_xyz[:, 1], "s--", lw=1.6, ms=4, label=r"$\kappa_{yy}$")
    ax.plot(T, k_xyz[:, 2], "^:", lw=1.4, ms=4,
            label=r"$\kappa_{zz}$ (no physical meaning in 2D)")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$\kappa$ (W/mK)")
    ax.set_title(f"{mat}: directional components", fontsize=11)
    ax.grid(alpha=.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    p = out / f"{mat}_kl_anisotropy.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    print("[fig]", p)


def write_tables(d, T, k_ii, k_n, mat, out: Path):
    with (out / f"{mat}_kl_vs_T.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_K", "kappa_inplane_W_mK", "kappa_2d_normalized_W_mK",
                    "kappa_xx_W_mK", "kappa_yy_W_mK", "kappa_zz_W_mK"])
        k_xyz = d.get("kappa_xx_yy_zz") or []
        for i, t in enumerate(T):
            row = [t, k_ii[i], k_n[i]]
            if i < len(k_xyz):
                row += list(k_xyz[i][:3])
            w.writerow(row)
    print("[table]", out / f"{mat}_kl_vs_T.csv")

    i300 = int(np.argmin(np.abs(T - 300))) if len(T) else 0
    mc = (d.get("mesh_convergence") or {}).get(d.get("bte_method") or "rta", {})
    rows = [
        ["material", mat],
        ["dim", d.get("dim")],
        ["bte_method", d.get("bte_method")],
        ["chosen_cutoff_A", d.get("chosen_cutoff_A")],
        ["mesh", d.get("mesh")],
        ["thickness_d_A", round(float(d.get("thickness_d_A") or 0), 4)],
        ["thickness_convention", d.get("thickness_convention")],
        ["stability_verdict", d.get("stability_verdict")],
        ["kappa_inplane_300K_W_mK", round(float(d.get("kappa_inplane_300K") or 0), 4)],
        ["kappa_2d_normalized_300K_W_mK",
         round(float(d.get("kappa_2d_normalized_inplane_300K") or 0), 4)],
        ["kappa_inf_extrap_300K_W_mK", round(float(mc.get("kappa_inf_extrap_300K_inplane") or 0), 4)],
        ["mesh_converged_5pct", mc.get("converged_5pct")],
        ["kappa_300K_from_T_series", round(float(k_ii[i300]), 4) if len(T) else ""],
    ]
    with (out / f"{mat}_kl_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["property", "value"])
        w.writerows(rows)
    (out / f"{mat}_kl_summary.json").write_text(
        json.dumps(dict(rows), ensure_ascii=False, indent=2), encoding="utf-8")
    print("[table]", out / f"{mat}_kl_summary.csv")


def copy_phonon_plots(root: Path, mat: str, out: Path):
    src = _find_plot_dir(root)
    if src is None:
        print("[skip] 无 S5.1_plot 声子谱产物")
        return
    dst = out / "phonon_band"
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in sorted(src.iterdir()):
        if f.is_file() and f.suffix.lower() in (".png", ".dat", ".yaml", ".json"):
            shutil.copy2(f, dst / f.name)
            n += 1
    print(f"[copy] {n} 个声子谱文件 -> {dst}")


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    root = Path(argv[1]).expanduser().resolve()
    out = Path(argv[2]).expanduser().resolve() if len(argv) > 2 else None
    mat = argv[3] if len(argv) > 3 else None

    sp = _find_summary(root)
    d = json.loads(sp.read_text(encoding="utf-8"))
    if mat is None:
        mat = root.name if root.name not in ("result", "kl-dft-cpu") else root.parent.name
    if out is None:
        base = root if root.name not in ("result", "kl-dft-cpu") else root.parent
        out = base / "kl_plots"
    out.mkdir(parents=True, exist_ok=True)
    print(f"[info] summary = {sp}")
    print(f"[info] out     = {out}")

    T = np.asarray(d.get("temperatures") or [], dtype=float)
    k_ii = _arr1(d, "kappa_inplane_xx_yy") if T.size else np.array([])
    k_n = _arr1(d, "kappa_2d_normalized_inplane_xx_yy") if T.size else np.array([])

    if T.size:
        fig_kappa_vs_T(T, k_ii, k_n, mat, out)
    fig_mesh(d, mat, out)
    fig_cutoff(d, mat, out)
    if T.size:
        fig_aniso(d, T, mat, out)
        write_tables(d, T, k_ii, k_n, mat, out)
    copy_phonon_plots(root, mat, out)
    print(f"[DONE] kappa plots -> {out}")


if __name__ == "__main__":
    main(sys.argv)