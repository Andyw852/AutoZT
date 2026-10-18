#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S6.1_plot: kl-dft-cpu 全流程出图归集（run: gen，登录节点，不交作业）。

★ 设计原则（2026-10-18 用户定）：
  1) **不改动**既有各步自动出的图（S5.1_plot 的声子谱 band_*.png 保持原样、仍在 step5_phonon_plot/）；
  2) 本步只做**归集**：把 S5.1 声子谱 + S6 κL 图/表，统一复制到
     `<材料>/kl-dft-cpu/phonon-transport_plot/`，作为交付给用户/论文的成品目录；
  3) 幂等：重复跑只覆盖同名文件，不删其它内容。

产物目录 `<材料>/kl-dft-cpu/phonon-transport_plot/`（与 stepN_* 同级）:
    phonon_band/            声子色散（S5.1 原图直拷，不改不重画）
    kappa/                  κL 图 + 表（本步绘制）
    README.md               说明（自动生成）

完成标记: phonon-transport_plot/plot_manifest.json
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

OUT_NAME = "phonon-transport_plot"


# ---------------- 定位 ----------------
def find_kappa_summary(cwd: Path) -> Path:
    cands = [cwd / "step6_kappa" / "kappa_summary.json"]
    cands += sorted(cwd.glob("*/step6_kappa/kappa_summary.json"))
    for c in cands:
        if c.is_file():
            return c
    sys.exit("[ERROR] 找不到 step6_kappa/kappa_summary.json（S6 未完成？）")


def find_phonon_plot_dir(cwd: Path) -> Path | None:
    for c in [cwd / "step5_phonon_plot"] + sorted(cwd.glob("*/step5_phonon_plot")):
        if c.is_dir() and any(c.glob("band_*.png")):
            return c
    return None


def _col(d, key, idx=0):
    v = d.get(key) or []
    return np.asarray([(x[idx] if isinstance(x, (list, tuple)) else x) for x in v],
                      dtype=float)


# ---------------- κL 绘图 ----------------
def plot_kappa(d, mat: str, out: Path):
    T = np.asarray(d.get("temperatures") or [], dtype=float)
    k_ii = _col(d, "kappa_inplane_xx_yy")
    k_n = _col(d, "kappa_2d_normalized_inplane_xx_yy")
    made = []

    if T.size and k_ii.size:
        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.plot(T, k_ii, "o-", color="#1f77b4", lw=2, ms=5,
                label=r"$\kappa_{L}$ in-plane (W/mK)")
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"$\kappa_{L}$ (W m$^{-1}$ K$^{-1}$)", color="#1f77b4")
        ax.tick_params(axis="y", labelcolor="#1f77b4")
        ax.grid(alpha=.3)
        ax2 = ax.twinx()
        if k_n.size:
            ax2.plot(T, k_n, "s--", color="#d62728", lw=1.6, ms=4,
                     label=r"2D-normalized $\kappa_{L}\cdot h/d$")
            ax2.set_ylabel(r"2D-normalized $\kappa_{L}$ (W/mK)", color="#d62728")
            ax2.tick_params(axis="y", labelcolor="#d62728")
        i300 = int(np.argmin(np.abs(T - 300)))
        ax.annotate(f"{k_ii[i300]:.2f}\n@300K", (T[i300], k_ii[i300]),
                    textcoords="offset points", xytext=(8, 10),
                    fontsize=8, color="#1f77b4")
        if k_n.size:
            ax2.annotate(f"{k_n[i300]:.2f}\n@300K", (T[i300], k_n[i300]),
                         textcoords="offset points", xytext=(8, -18),
                         fontsize=8, color="#d62728")
        lines = ax.get_lines() + ax2.get_lines()
        ax.legend(lines, [l.get_label() for l in lines], loc="upper right",
                  fontsize=9, framealpha=.9)
        ax.set_title(f"{mat}: lattice thermal conductivity (kl-dft-cpu, RTA)",
                     fontsize=11.5)
        fig.tight_layout()
        p = out / f"{mat}_kl_vs_T.png"
        fig.savefig(p, dpi=160); plt.close(fig); made.append(p)

        with (out / f"{mat}_kl_vs_T.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["T_K", "kappa_inplane_W_mK", "kappa_2d_normalized_W_mK",
                        "kappa_xx_W_mK", "kappa_yy_W_mK", "kappa_zz_W_mK"])
            k_xyz = d.get("kappa_xx_yy_zz") or []
            for i, t in enumerate(T):
                row = [t, k_ii[i], (k_n[i] if i < k_n.size else "")]
                if i < len(k_xyz):
                    row += list(k_xyz[i][:3])
                w.writerow(row)
        made.append(out / f"{mat}_kl_vs_T.csv")

    # 网格收敛
    mc = (d.get("mesh_convergence") or {}).get(d.get("bte_method") or "rta", {})
    meshes, k300 = mc.get("meshes") or [], mc.get("kappa300_inplane") or []
    if meshes and k300:
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
        ax.grid(alpha=.3); ax.legend(fontsize=9)
        fig.tight_layout()
        p = out / f"{mat}_kl_mesh_convergence.png"
        fig.savefig(p, dpi=160); plt.close(fig); made.append(p)

    # FC3 截断扫描
    cs = d.get("cutoff_selection") or {}
    cuts, kaps, errs, oks = [], [], [], []
    for r in (cs.get("records") or []):
        c, k = r.get("cut"), r.get("kappa")
        if c is None or k is None:
            continue
        cuts.append(float(c)); kaps.append(float(k))
        errs.append(float(r.get("kappa_err") or 0.0))
        oks.append(bool(r.get("plateau_ok")))
    if cuts:
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
        ax.set_title(f"{mat}: 3rd-order FC cutoff scan "
                     f"(status: {cs.get('status')})", fontsize=11)
        ax.grid(alpha=.3); ax.legend(fontsize=9)
        fig.tight_layout()
        p = out / f"{mat}_kl_cutoff_selection.png"
        fig.savefig(p, dpi=160); plt.close(fig); made.append(p)

    # 各向异性
    k_xyz = np.asarray(d.get("kappa_xx_yy_zz") or [], dtype=float)
    if T.size and k_xyz.ndim == 2 and k_xyz.shape[1] >= 3:
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        ax.plot(T, k_xyz[:, 0], "o-", lw=2, ms=5, label=r"$\kappa_{xx}$")
        ax.plot(T, k_xyz[:, 1], "s--", lw=1.6, ms=4, label=r"$\kappa_{yy}$")
        ax.plot(T, k_xyz[:, 2], "^:", lw=1.4, ms=4,
                label=r"$\kappa_{zz}$ (no physical meaning in 2D)")
        ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$\kappa$ (W/mK)")
        ax.set_title(f"{mat}: directional components", fontsize=11)
        ax.grid(alpha=.3); ax.legend(fontsize=9)
        fig.tight_layout()
        p = out / f"{mat}_kl_anisotropy.png"
        fig.savefig(p, dpi=160); plt.close(fig); made.append(p)

    # 汇总表
    if T.size:
        i300 = int(np.argmin(np.abs(T - 300)))
        rows = [
            ["material", mat], ["dim", d.get("dim")],
            ["bte_method", d.get("bte_method")],
            ["chosen_cutoff_A", d.get("chosen_cutoff_A")],
            ["mesh", d.get("mesh")],
            ["thickness_d_A", round(float(d.get("thickness_d_A") or 0), 4)],
            ["thickness_convention", d.get("thickness_convention")],
            ["stability_verdict", d.get("stability_verdict")],
            ["kappa_inplane_300K_W_mK", round(float(d.get("kappa_inplane_300K") or 0), 4)],
            ["kappa_2d_normalized_300K_W_mK",
             round(float(d.get("kappa_2d_normalized_inplane_300K") or 0), 4)],
            ["kappa_inf_extrap_300K_W_mK",
             round(float(mc.get("kappa_inf_extrap_300K_inplane") or 0), 4)],
            ["mesh_converged_5pct", mc.get("converged_5pct")],
            ["kappa_300K_from_T_series", round(float(k_ii[i300]), 4) if k_ii.size else ""],
        ]
        with (out / f"{mat}_kl_summary.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["property", "value"]); w.writerows(rows)
        (out / f"{mat}_kl_summary.json").write_text(
            json.dumps(dict(rows), ensure_ascii=False, indent=2), encoding="utf-8")
        made += [out / f"{mat}_kl_summary.csv", out / f"{mat}_kl_summary.json"]
    return made


# ---------------- 主流程 ----------------
def main():
    # autozt 的 run: gen 步骤在【技能目录】内执行：cwd = <材料>/kl-dft-cpu/
    # 所以产物必须写 cwd/phonon-transport_plot/（= <材料>/kl-dft-cpu/phonon-transport_plot/），
    # 材料名取父目录名。兼容被手动从材料根或 result 目录调起的情况。
    cwd = Path.cwd()
    if (cwd / "step6_kappa" / "kappa_summary.json").is_file():
        skill_dir = cwd                      # 正常：cwd 就是 kl-dft-cpu/
    elif (cwd / "kl-dft-cpu" / "step6_kappa" / "kappa_summary.json").is_file():
        skill_dir = cwd / "kl-dft-cpu"       # 从材料根调起
    else:
        skill_dir = cwd
    mat = skill_dir.parent.name if skill_dir.name in ("kl-dft-cpu", "result") else skill_dir.name
    out_root = skill_dir / OUT_NAME
    phonon_out = out_root / "phonon_band"
    kappa_out = out_root / "kappa"
    for d in (phonon_out, kappa_out):
        d.mkdir(parents=True, exist_ok=True)

    manifest = {"material": mat, "generator": "gen_step6.1_plot.py",
                "phonon_band": [], "kappa": [], "note": "S5.1 原图直拷、未重画"}

    # 1) 音子谱：原样复制（不改、不重画）
    src = find_phonon_plot_dir(skill_dir)
    if src:
        n = 0
        for f in sorted(src.iterdir()):
            if f.is_file() and f.suffix.lower() in (".png", ".dat", ".yaml", ".path", ".json"):
                shutil.copy2(f, phonon_out / f.name)
                manifest["phonon_band"].append(f.name)
                n += 1
        print(f"[copy] 声子谱 {n} 个文件 <- {src}（原图直拷，未重画）")
    else:
        print("[warn] 未找到 step5_phonon_plot 声子谱产物")

    # 2) κL 图/表
    sp = find_kappa_summary(skill_dir)
    d = json.loads(sp.read_text(encoding="utf-8"))
    if not str(d.get("KAPPA_DONE", "")).lower() in ("true", "1"):
        sys.exit(f"[ERROR] kappa_summary 未标记 KAPPA_DONE: {sp}")
    made = plot_kappa(d, mat, kappa_out)
    manifest["kappa"] = [p.name for p in made]
    print(f"[plot] κL 图/表 {len(made)} 个 -> {kappa_out}")

    # 3) README
    readme = [
        f"# {mat} —— phonon + transport plots",
        "",
        "本目录由 `kl-dft-cpu` 的 **S6.1_plot** 步骤自动生成（`run: gen`，不交作业）。",
        "",
        "- `phonon_band/`：S5.1_plot 已算好的声子色散图**原样复制**（不改动、不重画）",
        "- `kappa/`：S6_kappa 结果的 κL 图与表（本步绘制）",
        "",
        "## 文件",
        "",
    ]
    for f in manifest["phonon_band"]:
        readme.append(f"- `phonon_band/{f}`")
    for f in manifest["kappa"]:
        readme.append(f"- `kappa/{f}`")
    readme += ["", "## 说明", "",
               "- 声子谱图来自 S5.1_plot，保持技能原有出图不变；",
               "- κL 原始值依赖 Lz；2D 归一化 κ·h/d 与 Lz 无关，跨算例可比；",
               "- 2D 的 κzz 无物理意义（真空方向），仅作参考。", ""]
    (out_root / "README.md").write_text("\n".join(readme), encoding="utf-8")

    manifest["done"] = True
    (out_root / "plot_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[DONE] 出图归集 -> {out_root}")


if __name__ == "__main__":
    main()