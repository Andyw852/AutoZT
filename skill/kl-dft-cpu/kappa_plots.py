"""kl-dft-cpu 绘图核心（ke 风格：面板化 + 统一配色 + dpi150 + grid alpha .3）。

把两类收敛【分开画】，避免把 q 网格收敛和 FC3 截断收敛混为一谈：
  - q 网格收敛：BTE 的 q 点密度（本次 0.16%，收敛）
  - FC3 截断收敛：三阶力常数对距离截断（本次未到平台，status=not_converged）
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ke 同风格配色
C_KAPPA = "#1f77b4"
C_2D = "#d62728"
C_MESH = "#2ca02c"
C_CUT = "#9467bd"
C_OK = "#2ca02c"
C_WARN = "#ff7f0e"
C_REF = "#888888"

RC = {
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.grid": False,
    "figure.autolayout": False,
}


def _col(d, key, idx=0):
    v = d.get(key) or []
    return np.asarray([(x[idx] if isinstance(x, (list, tuple)) else x) for x in v],
                      dtype=float)


def _save(fig, path: Path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_overview(d, mat: str, out: Path):
    """一图总览（ke 风格多面板）：κL(T) / 2D 归一 / q 网格 / FC3 截断。"""
    with plt.rc_context(RC):
        T = np.asarray(d.get("temperatures") or [], dtype=float)
        k_ii = _col(d, "kappa_inplane_xx_yy")
        k_n = _col(d, "kappa_2d_normalized_inplane_xx_yy")
        bte = d.get("bte_method") or "rta"
        mc = (d.get("mesh_convergence") or {}).get(bte, {})
        cs = d.get("cutoff_selection") or {}

        fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
        ax1, ax2, ax3, ax4 = axes.flat

        # (a) kappa_L(T) 面内
        if T.size and k_ii.size:
            ax1.plot(T, k_ii, "o-", color=C_KAPPA, lw=2, ms=5, label=r"$\kappa_{L}$ in-plane")
            i300 = int(np.argmin(np.abs(T - 300)))
            ax1.plot([T[i300]], [k_ii[i300]], "o", ms=11, mfc="none",
                     mec=C_WARN, mew=2, zorder=5)
            ax1.annotate(f"{k_ii[i300]:.2f} W/mK @300K", (T[i300], k_ii[i300]),
                         textcoords="offset points", xytext=(12, 6),
                         fontsize=8.5, color=C_WARN)
        ax1.set_xlabel("Temperature (K)")
        ax1.set_ylabel(r"$\kappa_{L}$ (W m$^{-1}$ K$^{-1}$)")
        ax1.set_title("(a) in-plane thermal conductivity")
        ax1.grid(alpha=.3, which="both")
        ax1.legend(loc="upper right")

        # (b) 2D 归一化
        if T.size and k_n.size:
            ax2.plot(T, k_n, "s-", color=C_2D, lw=2, ms=4,
                     label=r"$\kappa_{L}\cdot h/d$")
            i300 = int(np.argmin(np.abs(T - 300)))
            ax2.annotate(f"{k_n[i300]:.1f} @300K", (T[i300], k_n[i300]),
                         textcoords="offset points", xytext=(12, 6),
                         fontsize=8.5, color=C_2D)
        ax2.set_xlabel("Temperature (K)")
        ax2.set_ylabel(r"2D-normalized $\kappa_{L}$ (W/mK)")
        ax2.set_title(r"(b) 2D-normalized $\kappa_{L}\cdot h/d$  (d=%.2f Å)"
                      % float(d.get("thickness_d_A") or 0))
        ax2.grid(alpha=.3, which="both")
        ax2.legend(loc="upper right")

        # (c) q 网格收敛
        meshes = mc.get("meshes") or []
        k300 = mc.get("kappa300_inplane") or []
        if meshes and k300:
            xs = [int(str(m).split()[0]) for m in meshes]
            ax3.plot(xs, k300, "o-", color=C_MESH, lw=2, ms=6)
            kinf = mc.get("kappa_inf_extrap_300K_inplane")
            if kinf:
                ax3.axhline(kinf, ls="--", lw=1.2, color=C_REF,
                            label=r"$\kappa_\infty$ = %.2f" % kinf)
            for x, y in zip(xs, k300):
                ax3.annotate("%.2f" % y, (x, y), textcoords="offset points",
                             xytext=(0, 8), ha="center", fontsize=7.5)
            ok = mc.get("converged_5pct")
            ax3.set_title("(c) q-mesh convergence  [%s]" % ("converged" if ok else "not converged"),
                          color=C_OK if ok else C_WARN)
            ax3.legend(loc="best")
        else:
            ax3.set_title("(c) q-mesh convergence  [no data]")
        ax3.set_xlabel("q-mesh N (N x N x 1)")
        ax3.set_ylabel(r"$\kappa_{L}$(300 K) (W/mK)")
        ax3.grid(alpha=.3, which="both")

        # (d) FC3 截断收敛（重点是 plateea 判据）
        recs = cs.get("records") or []
        cuts, kaps, errs, oks = [], [], [], []
        for r in recs:
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
            ax4.errorbar(cuts, kaps, yerr=errs, fmt="o-", lw=2, ms=5, capsize=3,
                         color=C_CUT, ecolor="#cccccc")
            for c, k, ok in zip(cuts, kaps, oks):
                if ok:
                    ax4.plot([c], [k], "o", ms=13, mfc="none", mec=C_OK, mew=1.8)
            chosen = d.get("chosen_cutoff_A") or cs.get("chosen_cut")
            if chosen:
                ax4.axvline(float(chosen), ls="--", color=C_2D, lw=1.1,
                            label="chosen = %s Å" % chosen)
            ax4.plot([], [], "o", ms=10, mfc="none", mec=C_OK,
                     label="plateau_ok (converged)")
            st = cs.get("status")
            ax4.set_title("(d) FC3 cutoff convergence  [%s]" % st,
                          color=C_OK if st == "converged" else C_WARN)
            ax4.legend(loc="best", fontsize=7.5)
        else:
            ax4.set_title("(d) FC3 cutoff convergence  [no data]")
        ax4.set_xlabel(r"FC3 pair cutoff (Å)")
        ax4.set_ylabel(r"$\kappa_{L}$(300 K) (W/mK)")
        ax4.grid(alpha=.3, which="both")

        fig.suptitle("%s — kl-dft-cpu: %s, thickness d = %.2f Å"
                     % (mat, bte.upper(), float(d.get("thickness_d_A") or 0)),
                     fontsize=12.5, y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.975])
        return _save(fig, out / f"{mat}_kl_overview.png")


def plot_anisotropy(d, mat: str, out: Path):
    """各向异性（2D 的 zz 明确标注无物理意义）。"""
    with plt.rc_context(RC):
        T = np.asarray(d.get("temperatures") or [], dtype=float)
        k_xyz = np.asarray(d.get("kappa_xx_yy_zz") or [], dtype=float)
        if not T.size or k_xyz.ndim != 2 or k_xyz.shape[1] < 3:
            return None
        dim2 = str(d.get("dim") or "").lower().startswith("2")
        fig, ax = plt.subplots(figsize=(6.6, 4.6))
        ax.plot(T, k_xyz[:, 0], "o-", lw=2, ms=5, color=C_KAPPA, label=r"$\kappa_{xx}$")
        ax.plot(T, k_xyz[:, 1], "s--", lw=1.6, ms=4, color=C_2D, label=r"$\kappa_{yy}$")
        ax.plot(T, k_xyz[:, 2], "^:", lw=1.3, ms=4, color=C_REF,
                label=r"$\kappa_{zz}$" + (" (vacuum, no physical meaning)" if dim2 else ""))
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
        ax.set_title("%s — directional components" % mat)
        ax.grid(alpha=.3, which="both")
        ax.legend(loc="upper right")
        return _save(fig, out / f"{mat}_kl_anisotropy.png")


def write_tables(d, mat: str, out: Path):
    T = np.asarray(d.get("temperatures") or [], dtype=float)
    k_ii = _col(d, "kappa_inplane_xx_yy")
    k_n = _col(d, "kappa_2d_normalized_inplane_xx_yy")
    made = []
    if T.size:
        p = out / f"{mat}_kl_vs_T.csv"
        with p.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["T_K", "kappa_inplane_W_mK", "kappa_2d_normalized_W_mK",
                        "kappa_xx_W_mK", "kappa_yy_W_mK", "kappa_zz_W_mK"])
            k_xyz = d.get("kappa_xx_yy_zz") or []
            for i, t in enumerate(T):
                row = [t, k_ii[i], (k_n[i] if i < k_n.size else "")]
                if i < len(k_xyz):
                    row += list(k_xyz[i][:3])
                w.writerow(row)
        made.append(p)

    bte = d.get("bte_method") or "rta"
    mc = (d.get("mesh_convergence") or {}).get(bte, {})
    cs = d.get("cutoff_selection") or {}
    i300 = int(np.argmin(np.abs(T - 300))) if T.size else 0
    rows = [
        ["material", mat], ["dim", d.get("dim")],
        ["bte_method", bte],
        ["thickness_d_A", round(float(d.get("thickness_d_A") or 0), 4)],
        ["thickness_convention", d.get("thickness_convention")],
        ["stability_verdict", d.get("stability_verdict")],
        ["kappa_inplane_300K_W_mK", round(float(d.get("kappa_inplane_300K") or 0), 4)],
        ["kappa_2d_normalized_300K_W_mK",
         round(float(d.get("kappa_2d_normalized_inplane_300K") or 0), 4)],
        ["mesh", d.get("mesh")],
        ["mesh_converged_5pct", mc.get("converged_5pct")],
        ["mesh_max_rel_change_pct", round(float(mc.get("max_rel_change_pct") or 0), 4)],  # 已是百分数
        ["kappa_inf_extrap_300K_W_mK",
         round(float(mc.get("kappa_inf_extrap_300K_inplane") or 0), 4)],
        ["fc3_cutoff_chosen_A", d.get("chosen_cutoff_A") or cs.get("chosen_cut")],
        ["fc3_cutoff_status", cs.get("status")],
        ["kappa_300K_from_T_series", round(float(k_ii[i300]), 4) if k_ii.size else ""],
    ]
    p = out / f"{mat}_kl_summary.csv"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["property", "value"]); w.writerows(rows)
    (out / f"{mat}_kl_summary.json").write_text(
        json.dumps(dict(rows), ensure_ascii=False, indent=2), encoding="utf-8")
    made += [p, out / f"{mat}_kl_summary.json"]
    return made