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
              transport_carrier_tuning_1e17-1e21/<M>_transport_carrier_tuning.png + _<T>K.csv
              transport_low_doping_<NLO>-1e16/<M>_transport_low_doping.png + _<T>K.csv
                        9温度x5参数(mu/S/sigma/PF/kappa_e)：DPT 路线 / AMSET phonon + impurity / ADP only；
                        [0056] 与 phonon_limited 一样按浓度拆两段、各一个文件夹
              phonon_limited_carrier_tuning_1e17-1e21/<M>_phonon_limited_carrier_tuning.png + _<T>K.csv
              phonon_limited_low_doping_<NLO>-1e16/<M>_phonon_limited_low_doping.png + _<T>K.csv
                        [0056] 2026-10-17 起按浓度拆两段、各一个文件夹（载流子调控区 / 低浓度段，后者有 DOPING_LOW 才有），
                        不再出整轴 <M>_phonon_limited.png；只含上列 4 条基础曲线
              <M>_params_<T>K.csv                          每个温度一张参数表(E1/C2D/m*/tau/mu)

数据来源:
    transport 图  <- step8.1_boltztrap/boltztrap_crta.json + step8.2_dpt/dpt_result.json
                   + step8.4_amset2d/transport.json + step8.4_amset2d/intrinsic_ADP.json
    phonon_limited 图 <- step8.3_output/intrinsic/intrinsic_transport_table.csv（内部文件名不变）
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
QUANTS_T = [("mu", r"$\mu$  (cm$^2$/V/s)"), ("S", r"S  ($\mu$V/K)"),
            ("sigma", r"$\sigma$  (S/m)"), ("PF", r"PF  (W/m/K$^2$)"),
            ("kappa_e", r"$\kappa_e$  (W/m/K)")]
QUANTS_I = [("mu_cm2/Vs", r"$\mu$  (cm$^2$/V/s)"), ("S_uV/K", r"S  ($\mu$V/K)"), ("sigma_S/m", r"$\sigma$  (S/m)"),
            ("PF_W/mK2", r"PF  (W/m/K$^2$)"), ("kappa_e_W/mK", r"$\kappa_e$  (W/m/K)")]  # 与 QUANTS_T 同列（符号统一为 mathtext）
DATASETS = [("overall", "tab:red", "--", "s"), ("intrinsic", "tab:orange", "-.", "D"),
            ("ADP_only", "tab:green", ":", "^"), ("BT2_DPTtau", "tab:blue", "-", "o")]
DS_LABELS = {"overall": "phonon + impurity", "intrinsic": "phonon-limited",      # [0054] 与 S8.3 统一命名
             "ADP_only": "ADP only (acoustic limit)", "BT2_DPTtau": "BT2×DPT-τ"}


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


AX = {"NLO": 1e17, "XT": [-4, -3, -2, -1, 0, 1, 2, 3, 4],
      "XTL": ["1e21", "1e20", "1e19", "1e18", "1e17", "1e18", "1e19", "1e20", "1e21"]}


def _axis_ticks(nlo, n_top):
    """正负对数浓度轴刻度（与 gen_step13_output._axis_ticks 同算法）。"""
    import math
    ndec = int(round(math.log10(n_top / nlo)))
    step = 1 if ndec <= 4 else 2
    pos = list(range(ndec, -1, -step))[::-1]
    if pos[-1] != ndec:
        pos.append(ndec)
    if pos[0] != 0:
        pos.insert(0, 0)
    xt = [-k for k in pos[::-1] if k] + pos
    lo = int(round(math.log10(nlo)))
    return xt, ["1e%d" % (lo + abs(k)) for k in xt]


# [0056] 与 gen_step13_output._phonon_limited_ranges 同：(标签, 下限, 上限, S 纵轴半宽)，上下限含端点
PL_TUNE_LO, PL_LOW_HI = 1e17, 1e16


def phonon_limited_ranges(nlo, n_top=1e21):
    import math
    out = [("carrier_tuning_1e17-1e%d" % round(math.log10(n_top)), PL_TUNE_LO, n_top, 1000.0)]   # [0056-r3] 上限 = 轴端
    if nlo < PL_LOW_HI:
        out.append(("low_doping_1e%d-1e16" % round(math.log10(nlo)), nlo, PL_LOW_HI, 2500.0))
    return out


# [0056-r3] 图注：每条曲线的散射机制、缩写全称、灰底含义（与 gen_step13_output._mech_caption 同口径）
MECH_DESC = {"ADP": "acoustic deformation potential (LA/TA acoustic phonons)",
             "POP": "polar optical phonon (LO, Fr\u00f6hlich)",
             "ODP": "optical deformation potential (non-polar optical phonons)",
             "IVS": "intervalley scattering", "IMP": "ionized impurity (N_ii = |n - p|)", "PIE": "piezoelectric"}


def _mechs(mdir):
    """(全机制, 去 IMP 机制)：全机制 = S8.4 transport.json 的 mobility 键；去 IMP 取 S8.3 本征表头注释，退回全机制 - IMP。"""
    import re
    allm = []
    for run in ("step8.4_amset2d", "step8_amset"):
        p = os.path.join(mdir, run, "transport.json")
        if os.path.isfile(p):
            try:
                allm = [str(k).upper() for k in (json.load(open(p)).get("mobility") or {}) if str(k).lower() != "overall"]
            except (OSError, ValueError):
                allm = []
            break
    intr = [m for m in allm if m != "IMP"]
    src = os.path.join(mdir, "step8.3_output", "intrinsic", "intrinsic_transport_table.csv")
    if os.path.isfile(src):
        for l in open(src, encoding="utf-8", errors="ignore"):
            if not l.startswith("#"):
                break
            mm = re.search(r"去 IMP：([A-Za-z+]+)", l)
            if mm:
                intr = mm.group(1).upper().split("+")
    return allm, intr


def _bad_ranges(mdir):
    """{n 型?: 最低的 in_interpret_range = NO 浓度}（S8.3 本征表；刚性能带近似不成立，与 S8.4 同判据）。"""
    src = os.path.join(mdir, "step8.3_output", "intrinsic", "intrinsic_transport_table.csv")
    bad = {}
    if os.path.isfile(src):
        for r in csv.DictReader([l for l in open(src) if not l.startswith("#") and l.strip()]):
            n = fnum(r.get("n_cell_cm-3"))
            if r.get("in_interpret_range") == "NO" and n == n and n:
                bad[n < 0] = min(bad.get(n < 0, abs(n)), abs(n))
    return bad


def _shade(axes, bad, _in):
    xe = AX["XT"][-1] + 0.15
    for neg, nb in bad.items():
        if not _in(np.array([nb]))[0]:
            continue
        x0 = float(xs(-nb if neg else nb))
        for row in axes:
            for ax in row:
                ax.axvspan(min(x0, -xe if neg else xe), max(x0, -xe if neg else xe), color="0.85", alpha=0.6, lw=0, zorder=0)


def _finish(fig, hh, cap, title, path, fs=11):
    """图注（下）-> 图例 -> 子图，按实际渲染高度排，不重叠。"""
    ct = fig.text(0.5, 0.004, "\n".join(cap), ha="center", va="bottom", fontsize=8.2, color="0.2", linespacing=1.35)
    rd = fig.canvas.get_renderer()
    y = ct.get_window_extent(rd).height / fig.bbox.height + 0.006
    lg = fig.legend(handles=hh, loc="lower center", bbox_to_anchor=(0.5, y), ncol=3, fontsize=9, frameon=False)
    y += lg.get_window_extent(rd).height / fig.bbox.height + 0.004
    fig.suptitle(title, fontsize=fs)
    fig.tight_layout(rect=[0, y, 1, 0.985 - 0.009 * title.count("\n")])
    fig.savefig(path, dpi=130); plt.close(fig)


def _abbr(used):
    ab = [m for m in dict.fromkeys(used) if m in MECH_DESC]
    return ["; ".join("%s = %s" % (m, MECH_DESC[m]) for m in ab)] if ab else []


X_CAP = ("x: doping n_cell (cm$^{-3}$, AMSET/BoltzTraP2 cell volume, incl. vacuum for 2D); "
         "left = n-type (electrons), right = p-type (holes)")
GRAY_CAP = ("Gray: beyond rigid-band range (in_interpret_range = NO in the S8.3 table, same criterion as S8.4) "
            "- results there are not reliable")


def set_axis(mdir, n_top=1e21):
    """[0055] 浓度轴跟着 AMSET 掺杂网格（与 gen_step13_output._doping_axis 同算法）：最小非零 |n| ≥ 1e17 时与旧版相同；
    有 DOPING_LOW 低浓度段时原点下探到最小点所在十年，十年数 > 4 时每 2 个十年一个刻度。"""
    import math
    global NFLOOR
    nz = []
    for run in ("step8.4_amset2d", "step8_amset"):
        p = os.path.join(mdir, run, "transport.json")
        if os.path.isfile(p):
            try:
                d = json.load(open(p))["doping"]
                nz = [abs(float(x)) for x in d if float(x) != 0.0]
            except (OSError, ValueError, KeyError, TypeError):
                nz = []
            break
    if nz and min(nz) < 1e17:
        nlo = 10.0 ** math.floor(math.log10(min(nz)) + 1e-9)
        NFLOOR = nlo * 0.999
    else:
        nlo, NFLOOR = 1e17, 5e17
    xt, xtl = _axis_ticks(nlo, n_top)
    AX.update(NLO=nlo, XT=xt, XTL=xtl)


def xs(n):
    n = np.asarray(n, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sign(n) * np.log10(np.maximum(np.abs(n), AX["NLO"]) / AX["NLO"])


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
            ax.set_xlim(AX["XT"][0] - 0.15, AX["XT"][-1] + 0.15)
            ax.set_xticks(AX["XT"])
            ax.set_xticklabels(AX["XTL"], fontsize=6.5)
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
    set_axis(mdir)                                          # [0055]
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
        # [0056-r4] 2026-10-18 用户定：BT2 曲线只画到 AMSET 掺杂网格的范围（同号支，|n| 在 AMSET 最小~最大档之间，
        #   与 S8.3 BT2×DPT-τ「不外推」同规则）。BT2 化学势扫描能一路扫到 1e21 以上，那一段 AMSET 没有点可比，
        #   低温（100–200 K）还会因 DOS 直方图分辨率不够出现锯齿（kT 小于 / 接近能量格宽）。
        nd = np.array(nd)
        _ok = np.zeros(len(nd), bool)
        for _neg in (True, False):
            _a = np.abs(np.array([x for x in OV["na"] if x != 0 and (x < 0) == _neg]))
            if _a.size:
                _ok |= ((nd < 0) == _neg) & (nd != 0) & (np.abs(nd) >= _a.min() * 0.999) & (np.abs(nd) <= _a.max() * 1.001)
        nd = np.where(_ok, nd, np.nan)
        data[T] = dict(nd=nd, Sd=np.array(Sd), gd=np.array(gd), pd=np.array(pd), kd=np.array(kd), OV=OV, AD=AD)
        for i in np.flatnonzero(_ok):
            br = "electron" if nd[i] < 0 else "hole"
            rows_out.append([T, "DPT", br, nd[i], Sd[i], gd[i], pd[i], kd[i], ""])
        for tag, X in (("AMSET_overall", OV), ("AMSET_ADP", AD)):
            if not X:
                continue
            for i in range(len(X["na"])):
                br = "electron" if X["na"][i] < 0 else "hole"
                rows_out.append([T, tag, br, X["na"][i], X["Sa"][i], X["ga"][i], X["pa"][i], X["ka"][i], X["ma"][i]])
    nlo0 = AX["NLO"]
    for tag, lo, hi, slim in phonon_limited_ranges(nlo0):   # [0056] 按浓度拆两段，各一个文件夹（图 + 每温度表）
        _in = lambda n: (np.abs(np.asarray(n, float)) >= lo * 0.999) & (np.abs(np.asarray(n, float)) <= hi * 1.001)
        xt, xtl = _axis_ticks(lo, hi)
        AX.update(NLO=lo, XT=xt, XTL=xtl)
        kind = tag.split("_1e")[0]
        pdir = os.path.join(out, "transport_" + tag)
        os.makedirs(pdir, exist_ok=True)
        for T in btT:                               # 每温度一个 CSV
            rows_T = [rows_out[0]] + [r for r in rows_out[1:] if abs(float(r[0]) - T) < 1e-6 and _in([r[3]])[0]]
            with open(os.path.join(pdir, "%s_transport_%s_%dK.csv" % (M, kind, int(round(T)))), "w", newline="") as fh:
                csv.writer(fh).writerows(rows_T)
        _plot_transport(data, btT, muedge, _in, 950 if lo >= 1e17 else slim, tag,
                        os.path.join(pdir, "%s_transport_%s.png" % (M, kind)), M, mdir)
    set_axis(mdir)
    return True


def _plot_transport(data, btT, muedge, _in, slim, tag, path, M, mdir):
    def cut(X, key):
        k = _in(X["na"])
        return X["na"][k], X[key][k]
    fig, axes = axes_grid(len(btT), 5, ["%d K" % t for t in btT], [q[1] for q in QUANTS_T])
    for rr, T in enumerate(btT):
        dd = data[T]
        for j, (q, _) in enumerate(QUANTS_T):
            ax = axes[rr][j]
            n = dd["nd"]; keep = _in(n)
            if q != "mu":
                seg(ax, n[keep], {"S": dd["Sd"], "sigma": dd["gd"], "PF": dd["pd"], "kappa_e": dd["kd"]}[q][keep],
                    color="tab:blue", lw=1.4)
            key = {"mu": "ma", "S": "Sa", "sigma": "ga", "PF": "pa", "kappa_e": "ka"}[q]
            seg(ax, *cut(dd["OV"], key), color="tab:red", lw=1.2, ls="--", marker="s", ms=3.5)
            if dd["AD"] is not None:
                seg(ax, *cut(dd["AD"], key), color="tab:green", lw=1.2, ls=":", marker="^", ms=3.5)
            if q == "mu":
                # [0056-r4] 用户 2026-10-18：横轴是掺杂浓度，左半 = n 型只该有电子、右半 = p 型只该有空穴；线上直接标注
                _xe = AX["XT"][-1] + 0.15
                for c, col, x0, x1, txt in (("electron", "tab:blue", -_xe, 0.0, "DPT $\\mu_{e}$ (band edge)"),
                                            ("hole", "tab:purple", 0.0, _xe, "DPT $\\mu_{h}$ (band edge)")):
                    if c in muedge:
                        yv = muedge[c] * (T0 / T)
                        ax.plot([x0, x1], [yv, yv], color=col, ls="-.", lw=0.9, alpha=0.85)
                        ax.text(x0 + 0.1 if c == "electron" else x1 - 0.1, yv * 1.15, txt, color=col, fontsize=6.5,
                                ha="left" if c == "electron" else "right", va="bottom",
                                bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.4))
                ax.set_yscale("log"); ax.set_ylim(2, 50000)
            elif q == "S":
                ax.set_ylim(-slim, slim); ax.axhline(0, color="k", lw=0.5)
            else:
                ax.set_yscale("log")
                if tag.startswith("carrier_tuning"):        # 低浓度段 σ/PF/κe 低好几个量级，用自动范围
                    ax.set_ylim({"sigma": (1e0, 1e10), "PF": (1e-6, 1e1), "kappa_e": (1e-4, 1e5)}[q])
    bad = _bad_ranges(mdir)
    _shade(axes, bad, _in)
    allm, _ = _mechs(mdir)
    _ov = "+".join(allm) if allm else "mechanisms not recorded"
    # [0056-r3] 用户 2026-10-18 看不懂旧图例（DPT route / DPT band-edge e- h+）：改成物理含义 + 图注逐条解释
    hh = [Line2D([], [], color="tab:red", lw=1.2, ls="--", marker="s", ms=4, label="AMSET phonon + impurity (%s)" % _ov),
          Line2D([], [], color="tab:green", lw=1.2, ls=":", marker="^", ms=4, label="AMSET ADP only (acoustic limit)"),
          Line2D([], [], color="tab:blue", lw=1.6, label="BT2×DPT-τ (S, σ, PF, κe)"),
          Line2D([], [], color="tab:blue", lw=1.0, ls="-.", label="DPT μ, electrons (band edge)"),
          Line2D([], [], color="tab:purple", lw=1.0, ls="-.", label="DPT μ, holes (band edge)")]
    cap = ["AMSET phonon + impurity: AMSET, %s (all mechanisms of this run).  AMSET ADP only: AMSET, ADP alone" % _ov,
           "BT2×DPT-τ: BoltzTraP2 band transport with one constant τ per carrier from DPT (Bardeen-Shockley, ADP only),"
           " τ(T) = τ(300 K)·300/T; x = BoltzTraP2 net carrier density n(μ) of each chemical-potential point",
           "DPT μ (dash-dot horizontal lines, μ panel only; electrons on the n-type half, holes on the p-type half):"
           " Bardeen-Shockley band-edge mobility μ = eτ/m*, μ(T) = μ(300 K)·300/T - one number per carrier and T,"
           " independent of n, hence horizontal",
           "BT2×DPT-τ is drawn only within the AMSET doping range of each side (no extrapolation);"
           " κe of BT2×DPT-τ = L0·σ·T with L0 = 2.44e-8 W Ω/K² (Wiedemann-Franz)"]
    cap += _abbr(allm + ["ADP"])
    if bad:
        cap.append(GRAY_CAP)
    cap.append(X_CAP)
    _trng = "%g-%g K" % (min(btT), max(btT)) if btT else ""     # [0052-r1] 温度范围按实际网格，不再写死 100-900 K
    _what = {"carrier_tuning": "carrier tuning (doping) range", "low_doping": "low-doping range (non-degenerate)"}
    _finish(fig, hh, cap, "%s  electronic transport, %s  n = %s - %s cm$^{-3}$  (%s)"
            % (M, _what[tag.split("_1e")[0]], AX["XTL"][len(AX["XTL"]) // 2], AX["XTL"][0], _trng), path, fs=12)


def make_intrinsic(mdir, out, M):
    set_axis(mdir)                                          # [0055]
    src = os.path.join(mdir, "step8.3_output", "intrinsic", "intrinsic_transport_table.csv")
    if not os.path.isfile(src):
        print("[skip] %s: 缺 intrinsic_transport_table.csv（S8.3 跑 0032+ 才有）" % M); return False
    rows = list(csv.DictReader([l for l in open(src) if not l.startswith("#") and l.strip()]))
    temps = sorted(set(float(r["T_K"]) for r in rows))
    cols = ["T_K", "dataset", "type", "n_cell_cm-3", "n_2D_cm-2", "n_lit_cm-3",
            "carriers_per_cell", "carriers_per_atom", "in_interpret_range",
            "S_uV/K", "sigma_S/m", "PF_W/mK2", "kappa_e_W/mK", "mu_cm2/Vs", "flag"]
    _names = {ds for ds, _c, _ls, _mk in DATASETS}   # [2026-10-15 选 2] intrinsic CSV 只含 4 条基础曲线
    nlo0 = AX["NLO"]
    for tag, lo, hi, slim in phonon_limited_ranges(nlo0):   # [0056] 按浓度拆两段，各一个文件夹
        _in = lambda n: (np.abs(n) >= lo * 0.999) & (np.abs(n) <= hi * 1.001)
        xt, xtl = _axis_ticks(lo, hi)
        AX.update(NLO=lo, XT=xt, XTL=xtl)
        kind = tag.split("_1e")[0]
        pdir = os.path.join(out, "phonon_limited_" + tag)
        os.makedirs(pdir, exist_ok=True)
        for T in temps:                             # 每温度一个 CSV
            with open(os.path.join(pdir, "%s_phonon_limited_%s_%dK.csv" % (M, kind, int(round(T)))), "w",
                      newline="") as fh:
                w = csv.writer(fh); w.writerow(cols)
                for r in rows:
                    if (abs(float(r["T_K"]) - T) > 1e-6 or r["dataset"] not in _names
                            or not _in(np.array([fnum(r["n_cell_cm-3"])]))[0]):
                        continue
                    w.writerow([r.get(c, "") for c in cols])
        _plot_phonon_limited(rows, temps, _in, slim, tag, os.path.join(pdir, "%s_phonon_limited_%s.png" % (M, kind)), M,
                             mdir)
    AX.update(NLO=nlo0)
    set_axis(mdir)
    return True


def _plot_phonon_limited(rows, temps, _in, slim, tag, path, M, mdir):
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
                keep = _in(n)
                seg(ax, n[keep], y[keep], color=color, lw=1.3, ls=ls, marker=mk, ms=3.2)
            if col == "S_uV/K":
                ax.set_ylim(-slim, slim); ax.axhline(0, color="k", lw=0.5)
            else:
                ax.set_yscale("log")
    bad = _bad_ranges(mdir)
    _shade(axes, bad, _in)
    allm, intr = _mechs(mdir)
    lab = dict(DS_LABELS)
    if allm:
        lab["overall"] = "phonon + impurity (%s)" % "+".join(allm)
    if intr:
        lab["intrinsic"] = "phonon-limited (%s)" % "+".join(intr)
    hh = [Line2D([], [], color=c, lw=1.4, ls=ls, marker=mk, ms=4, label=lab.get(ds, ds)) for ds, c, ls, mk in DATASETS]
    cap = ["phonon + impurity: AMSET, %s (all mechanisms of this run)" % ("+".join(allm) or "?"),
           "phonon-limited: AMSET, %s (IMP removed)" % ("+".join(intr) or "?"),
           "ADP only: AMSET, ADP alone (acoustic-phonon limit)",
           "BT2×DPT-τ: BoltzTraP2 band transport, constant relaxation time; τ from DPT (Bardeen-Shockley, ADP only,"
           " one energy-independent τ per carrier); μ masked where bipolar"]
    cap += _abbr(allm + ["ADP"])
    if bad:
        cap.append(GRAY_CAP)
    cap.append(X_CAP)
    # [0052-r1][0054] 文案与权威脚本 gen_step13_output.py 对齐：去 IMP = phonon-limited；
    #   "intrinsic" 只指本征载流子（transport_plots/intrinsic_carrier/）。
    _trng = "%g-%g K" % (min(temps), max(temps)) if temps else ""
    _what = {"carrier_tuning": "carrier tuning (doping) range", "low_doping": "low-doping range (non-degenerate)"}
    _finish(fig, hh, cap, "%s  phonon-limited transport vs doping, %s  n = %s - %s cm$^{-3}$  %s\n"
            "base 4 sets: phonon + impurity / phonon-limited / ADP only / BT2xDPT-tau;"
            " phonon-limited = all phonon mechanisms, no ionized impurity"
            % (M, _what[tag.split("_1e")[0]], AX["XTL"][len(AX["XTL"]) // 2], AX["XTL"][0], _trng), path)


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
    a = make_transport(mdir, out, M)          # [B] 两张图都要：transport + phonon_limited
    b = make_intrinsic(mdir, out, M)
    c = make_params(mdir, out, M)
    print("[%s] transport=%s phonon_limited=%s params=%s -> %s" % (M, a, b, c, out))


if __name__ == "__main__":
    main()
