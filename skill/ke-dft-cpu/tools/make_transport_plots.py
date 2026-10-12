#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ke-dft-cpu 结果 -> 规定格式的绘图输出（图 + 同目录数据）。

★ 同一套输出已【固化进 S8.3 技能】：gen_step13_output.write_transport_plots() 会在 S8.3 跑完时
  自动写 <材料>/transport_plots/，并由 autozt 的 fetch_extra（skill.yaml S8.3）拉回本地
  result/transport_plots/。本独立脚本用于：
    ① 对【已有结果】补生成（不必重跑 S8.3）；
    ② 离线/独立验证同一套图与数据。
  两者产物格式与文件名一致。

用法:
    python make_transport_plots.py <材料目录> <输出目录> [材料名]

<材料目录>  含各 step 目录的 ke-dft-cpu 根（本地已 fetch 的 result/ 或集群目录）；
            也可直接给以 step* 开头的目录（即 <M>/ke-dft-cpu/ 本身）。
<输出目录>  独立输出目录（约定: <M>/transport_plots/），生成:
              <M>_transport.png / <M>_transport_data_<T>K.csv   9温度x5参数(mu/S/sigma/PF/kappa_e)，
                        数据集 4 条：overall / intrinsic(ADP+POP) / ADP_only / BT2xDPT-tau
                        （intrinsic 图已并入本图，不再单独输出）
              <M>_params_<T>K.csv                          每个温度一张参数表(E1/C2D/m*/tau/mu)

数据来源:
    transport 图  <- step8.1_boltztrap/boltztrap_crta.json + step8.2_dpt/dpt_result.json
                   + step8.4_amset2d/transport.json + step8.4_amset2d/intrinsic_ADP.json
    intrinsic 图  <- step8.3_output/intrinsic/intrinsic_transport_table.csv
    params 表     <- step8.3_output/comparison_vs_T.csv

输出位置约定: 只写 <输出目录>，绝不写进 ke-dft-cpu/step8.3_output/。
"""
import csv, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

E = 1.602176634e-19
L = 2.44e-8
T0 = 300.0
NFLOOR = 5e17
QUANTS_T = [("mu", r"$\mu$  (cm$^2$/V/s)"), ("S", "S  (uV/K)"),
            ("sigma", r"$\sigma$  (S/m)"), ("PF", r"PF  (W/m/K$^2$)"),
            ("kappa_e", r"$\kappa_e$  (W/m/K)")]
QUANTS_I = [("mu_cm2/Vs", "mu  (cm^2/V/s)"), ("S_uV/K", "S  (uV/K)"), ("sigma_S/m", "sigma  (S/m)"),
            ("PF_W/mK2", "PF  (W/m/K$^2$)"), ("kappa_e_W/mK", "kappa_e  (W/m/K)")]  # 与 QUANTS_T 同列
DATASETS = [("overall", "tab:red", "--", "s"), ("intrinsic", "tab:orange", "-.", "D"),
            ("ADP_only", "tab:green", ":", "^"), ("BT2_DPTtau", "tab:blue", "-", "o")]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def xs(n):
    n = np.asarray(n, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sign(n) * np.log10(np.maximum(np.abs(n), 1e17) / 1e17)


def seg(ax, n, y, **kw):
    n = np.asarray(n, float); y = np.asarray(y, float)
    for msk in (n < 0, n > 0):
        mm = msk & np.isfinite(y) & np.isfinite(n)
        if mm.sum() < 2:
            continue
        o = np.argsort(xs(n[mm]))
        ax.plot(xs(n[mm])[o], y[mm][o], **kw)


def axes_grid(nrow, ncol, rowlab, collab):
    fig, axes = plt.subplots(nrow, ncol, figsize=(19.5, 2.35 * nrow), squeeze=False)
    for rr in range(nrow):
        for j in range(ncol):
            ax = axes[rr][j]
            ax.set_xlim(-4.15, 4.15)
            ax.set_xticks([-4, -3, -2, -1, 0, 1, 2, 3, 4])
            ax.set_xticklabels(["1e21", "1e20", "1e19", "1e18", "1e17",
                                "1e18", "1e19", "1e20", "1e21"], fontsize=6.5)
            ax.axvline(0, color="0.55", lw=0.6, ls=":")
            ax.grid(alpha=0.25, which="both", lw=0.4)
            if rr == 0:
                ax.set_title(collab[j], fontsize=10)
            if j == 0:
                ax.set_ylabel(rowlab[rr], fontsize=9, fontweight="bold")
            if rr == nrow - 1:
                ax.set_xlabel(r"$n$ (cm$^{-3}$)", fontsize=7)
    return fig, axes


def make_transport(mdir, out, M):
    bt_p = os.path.join(mdir, "step8.1_boltztrap", "boltztrap_crta.json")
    dp_p = os.path.join(mdir, "step8.2_dpt", "dpt_result.json")
    am_p = os.path.join(mdir, "step8.4_amset2d", "transport.json")
    # [0038] 只 ADP 的子集文件：subset_ADP.json（SUBSETS 生成）优先，退回旧的 intrinsic_ADP.json
    ia_p = os.path.join(mdir, "step8.4_amset2d", "subset_ADP.json")
    if not os.path.isfile(ia_p):
        ia_p = os.path.join(mdir, "step8.4_amset2d", "intrinsic_ADP.json")
    if not (os.path.isfile(bt_p) and os.path.isfile(am_p)):
        print("[skip] %s: 缺 boltztrap_crta.json 或 transport.json" % M); return False
    d = json.load(open(bt_p)); fac = float(d["twoD_c_over_t"])
    btT = [float(x) for x in d["temperatures_K"]]; mu_arr = d["mu_rel_efermi_eV"]
    carr_by_T = d["carrier_conc_cm-3"]; sot = d["sigma_over_tau_xx"]; sbk = d["seebeck_xx_V_per_K"]
    tau300, muedge = {}, {}
    if os.path.isfile(dp_p):
        for r in json.load(open(dp_p)).get("results", []):
            bd = r.get("by_direction") or {}
            if bd.get("status") == "ok" and isinstance(bd.get("x"), dict):
                tau300[r["carrier"]] = bd["x"]["tau_s"]; muedge[r["carrier"]] = bd["x"]["mobility_cm2_Vs"]
    am = json.load(open(am_p)); amT = [float(x) for x in am["temperatures"]]
    dop = [float(x) for x in am["doping"]]
    cond, asbk, kap, mob = am["conductivity"], am["seebeck"], am["electronic_thermal_conductivity"], am["mobility"]
    iad = json.load(open(ia_p))["intrinsic"] if os.path.isfile(ia_p) else None

    def c3(a, i, j):
        return float(a[i][j][0][0])

    rows_out = [["T_K", "source", "carrier", "n_cm-3", "S_uVK", "sigma_Sm", "PF_WmK2", "kappa_e_WmK", "mu_cm2Vs"]]
    data = {}
    for T in btT:
        iT = btT.index(T)
        jT = min(range(len(amT)), key=lambda k: abs(amT[k] - T))
        carr = carr_by_T[iT]
        nd, Sd, gd, pd, kd = [], [], [], [], []
        for i in range(len(mu_arr)):
            n = float(carr[i]); br = "electron" if n < 0 else "hole"
            S = float(sbk[iT][i]) * 1e6
            nd.append(n); Sd.append(S)
            tb = tau300.get(br)
            if tb is not None:
                sg = float(sot[iT][i]) * tb * (T0 / T) * fac
                gd.append(sg); pd.append((S * 1e-6) ** 2 * sg); kd.append(L * sg * T)
            else:
                gd.append(np.nan); pd.append(np.nan); kd.append(np.nan)

        def collect(cc, ss, kk, mm):
            na, Sa, ga, pa, ka, ma = [], [], [], [], [], []
            for i in range(len(dop)):
                s = c3(ss, i, jT) if ss else np.nan
                g = c3(cc, i, jT) * fac if cc else np.nan
                k = c3(kk, i, jT) * fac if kk else np.nan
                na.append(dop[i]); Sa.append(s); ga.append(g); ka.append(k)
                pa.append((s * 1e-6) ** 2 * g if np.isfinite(s) and np.isfinite(g) else np.nan)
                ma.append(float(mm[i][jT][0][0]))
            return dict(na=np.array(na), Sa=np.array(Sa), ga=np.array(ga), pa=np.array(pa), ka=np.array(ka), ma=np.array(ma))
        OV = collect(cond, asbk, kap, mob["overall"])
        AD = collect(iad["conductivity_S_m"], iad["seebeck_uV_K"], iad["electronic_thermal_conductivity_W_mK"],
                     iad["mobility_cm2_Vs_s"]["ADP"]) if iad else None
        data[T] = dict(nd=np.array(nd), Sd=np.array(Sd), gd=np.array(gd), pd=np.array(pd), kd=np.array(kd), OV=OV, AD=AD)
        for i in range(len(nd)):
            br = "electron" if nd[i] < 0 else "hole"
            rows_out.append([T, "DPT", br, nd[i], Sd[i], gd[i], pd[i], kd[i], ""])
        for tag, X in (("AMSET_overall", OV), ("AMSET_ADP", AD)):
            if not X:
                continue
            for i in range(len(X["na"])):
                br = "electron" if X["na"][i] < 0 else "hole"
                rows_out.append([T, tag, br, X["na"][i], X["Sa"][i], X["ga"][i], X["pa"][i], X["ka"][i], X["ma"][i]])
    for T in btT:                                   # 每温度一个 CSV
        rows_T = [rows_out[0]] + [r for r in rows_out[1:] if abs(float(r[0]) - T) < 1e-6]
        with open(os.path.join(out, "%s_transport_data_%dK.csv" % (M, int(round(T)))), "w", newline="") as fh:
            csv.writer(fh).writerows(rows_T)

    fig, axes = axes_grid(len(btT), 5, ["%d K" % t for t in btT], [q[1] for q in QUANTS_T])
    for rr, T in enumerate(btT):
        dd = data[T]
        for j, (q, _) in enumerate(QUANTS_T):
            ax = axes[rr][j]
            n = dd["nd"]; keep = np.abs(n) >= NFLOOR
            if q != "mu":
                seg(ax, n[keep], {"S": dd["Sd"], "sigma": dd["gd"], "PF": dd["pd"], "kappa_e": dd["kd"]}[q][keep],
                    color="tab:blue", lw=1.4)
            key = {"mu": "ma", "S": "Sa", "sigma": "ga", "PF": "pa", "kappa_e": "ka"}[q]
            seg(ax, dd["OV"]["na"], dd["OV"][key], color="tab:red", lw=1.2, ls="--", marker="s", ms=3.5)
            if dd["AD"] is not None:
                seg(ax, dd["AD"]["na"], dd["AD"][key], color="tab:green", lw=1.2, ls=":", marker="^", ms=3.5)
            if q == "mu":
                for c, col in (("electron", "tab:blue"), ("hole", "tab:purple")):
                    if c in muedge:
                        ax.axhline(muedge[c] * (T0 / T), color=col, ls="-.", lw=0.9, alpha=0.85)
                ax.set_yscale("log"); ax.set_ylim(2, 50000)
            elif q == "S":
                ax.set_ylim(-950, 950); ax.axhline(0, color="k", lw=0.5)
            else:
                ax.set_yscale("log")
                ax.set_ylim({"sigma": (1e0, 1e10), "PF": (1e-6, 1e1), "kappa_e": (1e-4, 1e5)}[q])
    hh = [Line2D([], [], color="tab:blue", lw=1.6, label="DPT route (S8.1 CRTA x DPT-tau)"),
          Line2D([], [], color="tab:red", lw=1.2, ls="--", marker="s", ms=4, label="AMSET overall"),
          Line2D([], [], color="tab:green", lw=1.2, ls=":", marker="^", ms=4, label="AMSET ADP only (intrinsic_ADP.json)"),
          Line2D([], [], color="tab:blue", lw=1.0, ls="-.", label="DPT band-edge e-"),
          Line2D([], [], color="tab:purple", lw=1.0, ls="-.", label="DPT band-edge h+")]
    fig.legend(handles=hh, loc="lower center", ncol=5, fontsize=9, frameon=False)
    fig.suptitle("%s  electronic transport  (100-900 K)   --   x = carrier conc. n,  left = n-type, right = p-type" % M, fontsize=12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.985])
    fig.savefig(os.path.join(out, "%s_transport.png" % M), dpi=130); plt.close(fig)
    return True


def make_intrinsic(mdir, out, M):
    src = os.path.join(mdir, "step8.3_output", "intrinsic", "intrinsic_transport_table.csv")
    if not os.path.isfile(src):
        print("[skip] %s: 缺 intrinsic_transport_table.csv（S8.3 跑 0032+ 才有）" % M); return False
    rows = list(csv.DictReader([l for l in open(src) if not l.startswith("#") and l.strip()]))
    temps = sorted(set(float(r["T_K"]) for r in rows))
    cols = ["T_K", "dataset", "type", "n_cell_cm-3", "n_2D_cm-2", "n_lit_cm-3",
            "carriers_per_cell", "carriers_per_atom", "in_interpret_range",
            "S_uV/K", "sigma_S/m", "PF_W/mK2", "kappa_e_W/mK", "mu_cm2/Vs", "flag"]
    for T in temps:                                 # 每温度一个 CSV
        with open(os.path.join(out, "%s_transport_data_%dK.csv" % (M, int(round(T)))), "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(cols)
            for r in rows:
                if abs(float(r["T_K"]) - T) < 1e-6:
                    w.writerow([r.get(c, "") for c in cols])
    fig, axes = axes_grid(len(temps), 5, ["%d K" % t for t in temps], [q[1] for q in QUANTS_I])
    for rr, T in enumerate(temps):
        sub = [r for r in rows if abs(float(r["T_K"]) - T) < 1e-6]
        for j, (col, _) in enumerate(QUANTS_I):
            ax = axes[rr][j]
            for ds, color, ls, mk in DATASETS:
                dd = [r for r in sub if r["dataset"] == ds]
                if not dd:
                    continue
                n = np.array([fnum(r["n_cell_cm-3"]) for r in dd])
                y = np.array([fnum(r[col]) for r in dd])
                keep = np.abs(n) >= NFLOOR
                seg(ax, n[keep], y[keep], color=color, lw=1.3, ls=ls, marker=mk, ms=3.2)
            if col == "S_uV/K":
                ax.set_ylim(-1000, 1000); ax.axhline(0, color="k", lw=0.5)
            else:
                ax.set_yscale("log")
    hh = [Line2D([], [], color=c, lw=1.4, ls=ls, marker=mk, ms=4, label=ds) for ds, c, ls, mk in DATASETS]
    hh.append(Line2D([], [], color="0.35", lw=1.0, ls=(0, (1, 1)), label="DPT (band-edge, mu only)"))
    fig.legend(handles=hh, loc="lower center", ncol=5, fontsize=9, frameon=False)
    fig.suptitle("%s  transport  (overall / ADP+POP / ADP only / BT2xDPT-tau)  100-900 K" % M, fontsize=12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.985])
    fig.savefig(os.path.join(out, "%s_transport.png" % M), dpi=130); plt.close(fig)
    return True


def make_params(mdir, out, M):
    src = os.path.join(mdir, "step8.3_output", "comparison_vs_T.csv")
    if not os.path.isfile(src):
        return False
    R = list(csv.DictReader([l for l in open(src) if not l.startswith("#") and l.strip()]))
    temps = sorted(set(int(float(r["T_K"])) for r in R))
    for T in temps:
        sub = [r for r in R if int(float(r["T_K"])) == T]
        with open(os.path.join(out, "%s_params_%dK.csv" % (M, T)), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Type", "Carrier", "Direction", "E1 (eV)", "C2D (N/m)", "m* (me)", "mu (cm2/V/s)", "tau (1e-15 s)"])
            for c in ("electron", "hole"):
                for r in [x for x in sub if x["carrier"] == c]:
                    w.writerow([M, "N" if c == "electron" else "P", r["direction"], r["E1_eV"],
                                r["C_2D_N_per_m"], r["m_eff_m0"], r["mu_DPT"], r["tau_DPT_fs"]])
    return True


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    mdir = os.path.abspath(sys.argv[1]); out = os.path.abspath(sys.argv[2])
    M = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(os.path.dirname(mdir).rstrip("/")) or os.path.basename(mdir.rstrip("/"))
    os.makedirs(out, exist_ok=True)
    b = make_intrinsic(mdir, out, M)
    c = make_params(mdir, out, M)
    print("[%s] transport=%s params=%s -> %s" % (M, b, c, out))


if __name__ == "__main__":
    main()
