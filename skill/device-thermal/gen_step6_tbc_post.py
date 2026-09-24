#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step6_tbc_post.py -- (可选, GPU) 从 GPUMD NEMD 产出提取界面热导 G（通用版）.

修复要点：热流 q 直接由源/漏恒温器的累积传能给出
  q = |dE_thermostat/dt| / A     (E 来自 compute.out 最后两列, eV)
不再需要文献 kappa（旧版用 kappa*slope 反推，非稳态时会差几个量级）。
kappa 若提供，仅作交叉核对。

读:
  step5_tbc/tbc_inputs.json   几何/源漏/面积/步长/轴
  step5_tbc/compute.out       恒温器累积传能（测量段的最后两列）
  step5_tbc/compute_chunk.out 逐 bin 温度剖面
产物:
  step6_tbc_post/tbc_result.json  含 G / dT_i / q / quality
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stepconf  # noqa: E402

STEP = "step6_tbc_post"
OUTDIR = STEP
SPEC = {
    "TBC_FIT_MARGIN_A": (2.0, "float"),
    "TBC_MIN_BINS": (3, "int"),
    "TBC_STEADY_TOL": (0.2, "float"),
    "TBC_MIN_WINDOW_FRAC": (0.2, "float"),
    "TBC_DT_TOL": (0.3, "float"),
    "TBC_MARGIN_TOL": (1.5, "float"),
    "TBC_KAPPA_A_W_mK": (None, "float"),
    "TBC_KAPPA_B_W_mK": (None, "float"),
}
EV_PER_FS_TO_W = 1.602176634e-4
A2_TO_M2 = 1e-20


def fit(xs, ys):
    n = len(xs)
    if n < 2:
        return 0.0, (ys[0] if ys else 0.0)
    mx = sum(xs) / n
    my = sum(ys) / n
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return 0.0, my
    a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
    return a, my - a * mx


def fit_stats(xs, ys, a, b, x0, n_blocks=5):
    """最小二乘拟合的斜率标准误与 x0 处预测值标准误；用**块平均**估计。

    相邻 bin 的 T 是空间相关的，把它们当独立点会用残差公式严重低估 ΔT_i 的误差；
    这里把按 x 排序的点切成 n_blocks 块、每块单独拟合，用各块斜率/预测值的离散度
    除以 sqrt(n_blocks) 作为标准误。点数不足返回 (None, None)。"""
    n = len(xs)
    if n <= 2:
        return None, None
    order = sorted(range(n), key=lambda i: xs[i])
    xs_s = [xs[i] for i in order]
    ys_s = [ys[i] for i in order]
    blk_slopes, blk_preds = [], []
    for ib in range(n_blocks):
        lo = ib * n // n_blocks
        hi = (ib + 1) * n // n_blocks
        if hi - lo < 2:
            continue
        ab, bb = fit(xs_s[lo:hi], ys_s[lo:hi])
        blk_slopes.append(ab)
        blk_preds.append(ab * x0 + bb)
    m = len(blk_slopes)
    if m < 2:
        return None, None
    se_slope = (sum((s - a) ** 2 for s in blk_slopes) / (m - 1)) ** 0.5 / (m ** 0.5)
    se_pred = (sum((p - (a * x0 + b)) ** 2 for p in blk_preds) / (m - 1)) ** 0.5 / (m ** 0.5)
    return se_slope, se_pred


def read_compute(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            p = ln.split()
            if len(p) >= 2:
                try:
                    rows.append([float(v) for v in p])
                except ValueError:
                    pass
    return rows


def read_profile(path, start_block=0):
    """读逐 bin 温度剖面；start_block 起只平均稳态窗口内的 block（与 q 的时间窗对齐）。"""
    blocks, cur = [], {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            try:
                cid = int(float(p[0]))
                c = float(p[1])
                cnt = float(p[2])
                t = float(p[3])
            except ValueError:
                continue
            if cid == 0 and cur:
                blocks.append(cur)
                cur = {}
            cur[cid] = (c, t, cnt)
    if cur:
        blocks.append(cur)
    if not blocks:
        sys.exit("[ERROR] compute_chunk.out 无有效数据（NEMD 是否已运行？）")
    sel = blocks[start_block:]
    out = {}
    for cid in sorted(blocks[0]):
        cs = [b[cid][0] for b in sel if cid in b]
        ts = [b[cid][1] for b in sel if cid in b]
        ns = [b[cid][2] for b in sel if cid in b]
        if ts:
            out[cid] = (sum(cs) / len(cs), sum(ts) / len(ts), sum(ns) / len(ns))
    return out, len(sel)


def best_steady_window(Es, Ek, minlen, tol):
    """在 [k,n) 后缀里找源/漏热流自洽的窗口；优先返回**最长**的可接受窗口。
    找不到可接受窗口时返回最自洽的那个（供 quality=poor 报告）。
    k 从 n//3 起（至少丢掉前 1/3 暂态）：源/漏热流自洽在 Langevin 恒温器下太容易满足，
    单凭它不足以证明到了稳态。"""
    n = len(Es)
    best_any = None
    best_ok = None
    k_start = n // 3
    k_hi = max(0, n - minlen)
    if k_start > k_hi:
        k_start = k_hi
    for k in range(k_start, k_hi + 1):
        xs = list(range(k, n))
        aS, _ = fit(xs, Es[k:])
        aK, _ = fit(xs, Ek[k:])
        qS, qK = abs(aS), abs(aK)
        incons = abs(qS - qK) / max(qS, qK, 1e-30)
        if best_any is None or incons < best_any[3]:
            best_any = (k, qS, qK, incons)
        if incons <= tol and (best_ok is None or k < best_ok[0]):
            best_ok = (k, qS, qK, incons)
    return best_ok if best_ok is not None else best_any


def main():
    conf = stepconf.load(SPEC, STEP, strict=False)
    cwd = os.getcwd()
    meta_p = os.path.join(cwd, "step5_tbc", "tbc_inputs.json")
    comp_p = os.path.join(cwd, "step5_tbc", "compute.out")
    prof_p = os.path.join(cwd, "step5_tbc", "compute_chunk.out")
    for p, what in ((meta_p, "step5_tbc/tbc_inputs.json（先跑 S5_tbc）"),
                    (comp_p, "step5_tbc/compute.out（先跑 GPUMD NEMD）"),
                    (prof_p, "step5_tbc/compute_chunk.out（先跑 GPUMD NEMD）")):
        if not os.path.isfile(p):
            sys.exit("[ERROR] 缺少 %s" % what)
    with open(meta_p, encoding="utf-8") as fh:
        meta = json.load(fh)

    axis = meta.get("axis", "z")
    cif = float(meta["interface_coord_A"])
    cmin, cmax = float(meta["coord_min_A"]), float(meta["coord_max_A"])
    src_th = float(meta.get("source_thickness_A") or 0.0)
    snk_th = float(meta.get("sink_thickness_A") or 0.0)
    side = (meta.get("source_side") or "high").lower()
    area = float(meta["area_A2"])
    dt_fs = float(meta.get("time_step_fs", 1.0)) * int(meta.get("sample_interval", 1)) \
        * int(meta.get("output_interval", 1))

    rows = read_compute(comp_p)
    if len(rows) < 4:
        sys.exit("[ERROR] compute.out 数据太少（%d 行）；加长 TBC_RUN_STEPS 或减小 "
                 "TBC_SAMPLE/OUTPUT_INTERVAL" % len(rows))
    ncol = len(rows[0])
    Es = [r[ncol - 2] for r in rows]      # 源恒温器累积传能 (eV)
    Ek = [r[ncol - 1] for r in rows]      # 漏恒温器累积传能 (eV)
    n = len(rows)
    minlen = max(4, int(float(conf["TBC_MIN_WINDOW_FRAC"]) * n))
    k, qS, qK, incons = best_steady_window(Es, Ek, minlen, float(conf["TBC_STEADY_TOL"]))
    qsrc = (qS / dt_fs) * EV_PER_FS_TO_W / (area * A2_TO_M2)   # eV/行 -> eV/fs -> W/m2
    qsnk = (qK / dt_fs) * EV_PER_FS_TO_W / (area * A2_TO_M2)
    q = 0.5 * (qsrc + qsnk)

    prof, nblocks = read_profile(prof_p, start_block=k)   # 与 q 的稳态窗口对齐
    minbins = int(conf["TBC_MIN_BINS"])
    reasons = []

    def evaluate(margin):
        """按给定拟合余量算一次（点数不足时 dTi=0）。"""
        if side == "high":
            hlo, hhi = cif + margin, cmax - src_th - margin
            clo, chi = cmin + snk_th + margin, cif - margin
        else:
            hlo, hhi = cmin + src_th + margin, cif - margin
            clo, chi = cif + margin, cmax - snk_th - margin
        hot = [(v[0], v[1]) for v in prof.values() if hlo < v[0] < hhi and v[2] > 0.5]
        cold = [(v[0], v[1]) for v in prof.values() if clo < v[0] < chi and v[2] > 0.5]
        d = {"margin": margin, "n_hot": len(hot), "n_cold": len(cold),
             "aH": 0.0, "aC": 0.0, "Th": 0.0, "Tc": 0.0, "dTi": 0.0,
             "slH": None, "slC": None, "seH": None, "seC": None,
             "hot": hot, "cold": cold}
        if len(hot) >= 2 and len(cold) >= 2:
            xsH = [c for c, _ in hot]
            ysH = [t for _, t in hot]
            xsC = [c for c, _ in cold]
            ysC = [t for _, t in cold]
            aH, bH = fit(xsH, ysH)
            aC, bC = fit(xsC, ysC)
            slH, seH = fit_stats(xsH, ysH, aH, bH, cif)
            slC, seC = fit_stats(xsC, ysC, aC, bC, cif)
            d.update(aH=aH, aC=aC, Th=aH * cif + bH, Tc=aC * cif + bC,
                     dTi=(aH * cif + bH) - (aC * cif + bC),
                     slH=slH, slC=slC, seH=seH, seC=seC)
        return d

    m0 = float(conf["TBC_FIT_MARGIN_A"])
    ev = evaluate(m0)
    hot, cold = ev["hot"], ev["cold"]
    aH, aC, Th, Tc = ev["aH"], ev["aC"], ev["Th"], ev["Tc"]
    slH, slC, seH, seC = ev["slH"], ev["slC"], ev["seH"], ev["seC"]
    dTi = ev["dTi"]
    G = q / abs(dTi) if dTi else 0.0
    dT_unc = ((seH or 0.0) ** 2 + (seC or 0.0) ** 2) ** 0.5 if (seH is not None and seC is not None) else None
    G_rel = (dT_unc / abs(dTi)) if (dT_unc is not None and dTi) else None

    # 拟合余量敏感性：真正扩散的引线 G 不该随余量变化；短/弹道引线会剧烈变化
    scan = []
    for mg in sorted(set([1.0, 2.0, 4.0, m0])):
        e = evaluate(mg)
        scan.append({"margin_A": mg, "G_W_m2K": (q / abs(e["dTi"]) if e["dTi"] else None),
                     "dT_interface_K": e["dTi"], "n_bins_hot": e["n_hot"],
                     "n_bins_cold": e["n_cold"]})
    gs = [s["G_W_m2K"] for s in scan if s["G_W_m2K"]]
    g_ratio = (max(gs) / min(gs)) if len(gs) >= 2 and min(gs) > 0 else None

    tol = float(conf["TBC_STEADY_TOL"])
    dt_tol = float(conf["TBC_DT_TOL"])
    mg_tol = float(conf["TBC_MARGIN_TOL"])
    if len(hot) < minbins:
        reasons.append("hot bins %d < %d" % (len(hot), minbins))
    if len(cold) < minbins:
        reasons.append("cold bins %d < %d" % (len(cold), minbins))
    if incons > tol:
        reasons.append("q 源漏不自洽 %.2f > %.2f" % (incons, tol))
    if dTi <= 0:
        reasons.append("dT_i<=0（热侧未更热/方向反了）")
    if G_rel is not None and G_rel > dt_tol:
        reasons.append("dT_i 外推相对误差 %.2f > %.2f" % (G_rel, dt_tol))
    if g_ratio is not None and g_ratio > mg_tol:
        reasons.append("G 对拟合余量敏感 max/min=%.2f > %.1f（引线太短/非扩散区，G 不可靠）"
                       % (g_ratio, mg_tol))
    quality = "good" if not reasons else "poor"

    # kappa 交叉核对（可选；短引线/弹道区会差很多，仅信息）
    ka = conf["TBC_KAPPA_A_W_mK"] if conf["TBC_KAPPA_A_W_mK"] is not None else meta.get("kappa_a_W_mK")
    kb = conf["TBC_KAPPA_B_W_mK"] if conf["TBC_KAPPA_B_W_mK"] is not None else meta.get("kappa_b_W_mK")
    xcheck = {}
    if ka and kb and len(hot) >= 2 and len(cold) >= 2:
        qk_h = abs(float(ka) * aH * 1e10)
        qk_c = abs(float(kb) * aC * 1e10)
        xcheck = {"kappa_hot_W_mK": ka, "kappa_cold_W_mK": kb,
                  "q_kappa_hot_W_m2": qk_h, "q_kappa_cold_W_m2": qk_c,
                  "ratio_thermo_over_kappa_hot": (q / qk_h if qk_h else None),
                  "ratio_thermo_over_kappa_cold": (q / qk_c if qk_c else None),
                  "note": "kappa 法仅在扩散极限成立；短引线/弹道区比值会远离 1"}

    out = {"TBC_POST_DONE": True, "axis": axis,
           "G_W_m2K": G, "R_interface_m2K_W": (1.0 / G if G else None),
           "q_W_m2": q, "q_from_source_W_m2": qsrc, "q_from_sink_W_m2": qsnk,
           "q_consistency": incons, "quality": quality, "quality_reasons": reasons,
           "dT_interface_K": dTi, "T_hot_K": Th, "T_cold_K": Tc,
           "slope_hot_K_per_A": aH, "slope_cold_K_per_A": aC,
           "slope_hot_std_K_per_A": slH, "slope_cold_std_K_per_A": slC,
           "dT_interface_uncertainty_K": dT_unc, "G_relative_uncertainty": G_rel,
           "G_margin_scan": scan, "G_margin_ratio": g_ratio,
           "source_side": side, "interface_coord_A": cif, "area_A2": area,
           "n_bins_hot": len(hot), "n_bins_cold": len(cold), "n_blocks_averaged": nblocks,
           "steady_window": {"rows_used": n - k, "rows_total": n,
                             "start_ps": k * dt_fs * 1e-3, "total_ps": n * dt_fs * 1e-3,
                             "discarded_frac": (k / n if n else 0.0)},
           "kappa_crosscheck": xcheck,
           "note": "q 取自源/漏恒温器累积传能（不依赖文献 kappa）；dT_i 为两侧体区线性外推。"
                   "quality=poor 时 G 仅作跑通性参考。"}

    od = os.path.join(cwd, OUTDIR)
    os.makedirs(od, exist_ok=True)
    p = os.path.join(od, "tbc_result.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)
    print("[OK] %s" % p)
    print("    q = %.4g W/m2 (源 %.4g / 漏 %.4g, 自洽度 %.3f, 用后 %d/%d 行)"
          % (q, qsrc, qsnk, incons, n - k, n))
    print("    dT_interface = %.3f K (T_hot=%.2f T_cold=%.2f) ; G = %.4g MW/m2K (R=%.3g)"
          % (dTi, Th, Tc, G / 1e6, out["R_interface_m2K_W"] or 0))
    print("    quality=%s%s" % (quality, ("  原因: " + "; ".join(reasons)) if reasons else ""))


if __name__ == "__main__":
    main()
