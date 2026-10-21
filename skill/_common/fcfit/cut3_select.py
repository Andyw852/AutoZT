#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""cut3_select.py -- pick the fc3 cutoff from per-candidate fit + kappa records.

Three criteria (same as kl-dft-cpu S6): (1) residual 1-SE, (2) per-shell
stability (cut <= stable_upper_cut), (3) kappa plateau between adjacent cuts.
This module is deliberately dependency-free (stdlib math only) so it can be
pushed to a compute node alongside kappa_driver.py.  It mirrors
kl_common.select_cutoff -- keep the two in sync (TODO: move to skill/_common).
"""
import math


def _plateau_pair(a, b, kappa_tol_pct, floor):
    """相邻两截断是否在平台容差内。

    有 kappa_vec（xx/yy/zz）时**逐分量**判：各向异性体系面外 zz 往往收敛更慢，只看
    面内 κ 会把它掩盖（同 S2 网格收敛的 per_component 判据）。分母取
    max(|κ_ii|, floor·max_j|κ_jj|)：2D 真空轴、极小分量不会因自身相对噪声卡住平台。
    没有 kappa_vec（老记录）时退回只看标量 kappa。"""
    frac = float(kappa_tol_pct) / 100.0
    err = math.hypot(a["kappa_err"], b["kappa_err"])
    va, vb = a.get("kappa_vec"), b.get("kappa_vec")
    if va and vb and len(va) == len(vb):
        # kappa_err 是标量 κ 的绝对误差，按相对量折算到各分量
        rel_err = err / abs(a["kappa"]) if a["kappa"] else 0.0
        big = max([abs(x) for x in va if x is not None] or [0.0])
        for x, y in zip(va, vb):
            if x is None or y is None:
                continue
            den = max(abs(x), float(floor) * big)
            if den > 0 and abs(y - x) > max(frac, rel_err) * den + 1e-15:
                return False
        return True
    ka, kb = a["kappa"], b["kappa"]
    if ka is None or kb is None:
        return True
    tol = max(frac * abs(ka), err)
    return abs(kb - ka) <= tol + 1e-15


def _plateau_runs(usable, pair_ok, kappa_tol_pct, floor):
    """连续稳定区间 [(i, j), ...]（usable 下标，j > i）。

    先按相邻判据 pair_ok 切出最长的连续段，再把每段起点右移到"与段内其后每一档
    累计偏差都在容差内"的第一个位置——相邻都 <5% 但一路同向漂移的段不会整段算平。"""
    runs, i, n = [], 0, len(usable)
    while i < n - 1:
        if not pair_ok[i]:
            i += 1
            continue
        j = i + 1
        while j < n - 1 and pair_ok[j]:
            j += 1
        for k in range(i, j):
            if all(_plateau_pair(usable[k], usable[m], kappa_tol_pct, floor)
                   for m in range(k + 1, j + 1)):
                runs.append((k, j))
                break
        i = j
    return runs


# 数据量门槛（方程数/参数数），默认 3（用户定的标准，与 kl S4 的 OVERSAMPLE=3 对齐）。
#   可用 CUT3_MIN_RATIO 调高：合成数据检验（2D 5x5 超胞、36 帧、1 meV/Å 噪声）里
#   ratio≈3.5 的档 κ 仍被噪声压低约 10%，ratio≈2.2 压低近 40%。
MIN_RATIO_DEFAULT = 3.0


def select_cutoff(records, kappa_tol_pct=5.0, se_mult=1.0, stability_thr=0.3,
                  pick="largest", component_floor=0.1, min_ratio=MIN_RATIO_DEFAULT):
    """按三判据自动选最小可行截断。records 是每个候选截断一条 dict：

        cut               截断 (Å)
        ratio             方程数/参数数（< min_ratio，默认 3 = 数据不足 → 不参与）
        train_rel_err     拟合残差（★in-sample，**不是**交叉验证）：pheasy -f 在【全部帧】
                          上的 Relative error；判据①的输入
        train_rel_se      该残差在帧 bootstrap 重抽样下的标准误（1-SE 用；缺省 0）
        err_kind          "in-sample"（默认，判据①不具否决权）| "held-out"（真留出 CV）
        stable_upper_cut  数据能确定的最大截断 (Å)：只用 σ/|mean| < stability_thr 的壳层
                          （判据②；该截断 ≤ 它才算"数据能确定"）
        kappa             300 K 面内 κ
        kappa_err         κ 的 bootstrap 误差棒（判据③；缺省 0）

    判据③（κ 平台）找的是【连续稳定区间】：相邻 Δκ 与区间内累计 Δκ 都在容差内的一段
    连续截断；区间之后的跳变只记进 report["after_plateau"]，不影响选择。

    pick：平台区间里满足判据的截断有多个时，"largest"（2026-10-15 起默认：取区间内
          能确定的最大截断 —— 三阶相互作用最全、平台内 κ 最稳，且与网格收敛的参考档更可能
          重合，不再需要为选中档单独复核网格）或 "smallest"（旧默认：取区间内最小）。

    返回 (chosen_cut, report)：chosen_cut 为 None 表示无可用截断。report 给出
    每个截断是否 data_ok / stable_ok / err_ok / plateau_ok 与最终判定依据。
    report["err_gate_applied"]=False 表示 err_kind 不是 held-out → 判据①只展示、不筛选。
    """
    recs = []
    for r in records:
        d = dict(r)
        d["cut"] = float(d.get("cut"))
        d["ratio"] = float(d.get("ratio") or 0.0)
        d["train_rel_err"] = (None if d.get("train_rel_err") is None
                              else float(d["train_rel_err"]))
        d["train_rel_se"] = float(d.get("train_rel_se") or 0.0)
        d["err_kind"] = str(d.get("err_kind") or "in-sample")
        d["stable_upper_cut"] = (None if d.get("stable_upper_cut") is None
                                 else float(d["stable_upper_cut"]))
        # stability_measured=False：没有 bootstrap 重拟合（CUT3_BOOTSTRAP=0），判据②
        # 无从计算——跳过它、只靠判据③（κ 平台），与 phono3py 逐档扫描同一口径。
        # 以前 None 一律判 stable_ok=False，导致默认能跑的扫描全部"无可用截断"。
        d["stability_measured"] = bool(d.get("stability_measured", True))
        d["kappa"] = (None if d.get("kappa") is None else float(d["kappa"]))
        d["kappa_err"] = float(d.get("kappa_err") or 0.0)
        kv = d.get("kappa_vec")
        d["kappa_vec"] = ([None if x is None else float(x) for x in kv][:3]
                          if kv else None)
        recs.append(d)
    recs.sort(key=lambda r: r["cut"])
    _min_ratio = float(MIN_RATIO_DEFAULT if min_ratio is None else min_ratio)
    for r in recs:
        r["data_ok"] = r["ratio"] >= _min_ratio
        r["stable_ok"] = ((not r["stability_measured"])
                          or (r["stable_upper_cut"] is not None
                              and r["cut"] <= r["stable_upper_cut"] + 1e-9))
    usable = [r for r in recs if r["data_ok"] and r["stable_ok"]]

    # 判据①（残差 1-SE）：**只在误差是真·留出(CV)误差时才具否决权**。
    #   ★ 2026-09-27 修：该字段旧名 cv_err，但值来自 pheasy -f 在【全部帧】上的
    #     Relative error —— 是 in-sample 训练残差，会随参数增多单调下降；拿它做
    #     "≤ min+SE"的否决会系统性偏向大截断（长程噪声项）。故 err_kind != "held-out"
    #     时判据①只展示、不参与筛选（err_gate_applied=False）。
    err_gate = bool(usable) and all(
        str(r.get("err_kind") or "in-sample").lower() == "held-out" for r in usable)
    if not usable:
        err_ok = set()
    elif not err_gate:
        err_ok = {r["cut"] for r in usable}
    else:
        _e = [r["train_rel_err"] for r in usable if r["train_rel_err"] is not None]
        if _e:
            _emin = min(_e)
            _rec = next(r for r in usable if r["train_rel_err"] == _emin)
            _thr = _emin + float(se_mult) * float(_rec["train_rel_se"])
            err_ok = {r["cut"] for r in usable
                      if r["train_rel_err"] is not None
                      and r["train_rel_err"] <= _thr + 1e-15}
        else:
            err_ok = {r["cut"] for r in usable}

    # 判据③：κ 平台 = 一段【连续稳定区间】（2026-10-14 用户定）。
    #   相邻两档 Δκ ≤ max(tol·|κ|, 误差) 的连续档组成平台段；段内再从小往大找起点：
    #   起点与段内其后每一档的【累计】偏差也都在容差内（防小步慢漂）。
    #   平台段之后的跳变（如最外一档突跌近 20%：长程三阶项数据定不住）只写进报告，
    #   不再像旧规则那样要求"一直平到最外侧"而把整条判成 not_converged、退回最外档。
    #   有多段平台时取最长的一段（档数相同取截断小的一段）；单独一档不构成平台。
    # κ 缺失（该档 BTE 没算出来）不算稳定——旧规则把 None 当"平"放行
    pair_ok = [a["kappa"] is not None and b["kappa"] is not None
               and _plateau_pair(a, b, kappa_tol_pct, component_floor)
               for a, b in zip(usable, usable[1:])]
    runs = _plateau_runs(usable, pair_ok, kappa_tol_pct, component_floor)
    runs.sort(key=lambda ij: (-(ij[1] - ij[0]), ij[0]))
    plateau, run_used = set(), None
    for i, j in runs:
        pts = {r["cut"] for r in usable[i:j + 1]}
        if pts & err_ok:
            plateau, run_used = pts, (i, j)
            break

    for r in recs:
        r["err_ok"] = r["cut"] in err_ok
        r["plateau_ok"] = r["cut"] in plateau

    inter = sorted(err_ok & plateau)
    after = []
    if inter:
        if str(pick).lower() == "largest":
            chosen = inter[-1]     # Mg8C120 口径：平台区间内能确定的最大截断
        else:
            chosen = inter[0]      # 总则：取平台区间内满足判据的最小截断
        status = "ok"
        _end = usable[run_used[1]]
        for r in usable[run_used[1] + 1:]:
            if _end["kappa"] and r["kappa"] is not None:
                after.append({"cut": r["cut"], "kappa": r["kappa"],
                              "change_pct_vs_plateau_end":
                                  round(100.0 * (r["kappa"] - _end["kappa"]) / abs(_end["kappa"]), 2)})
    elif usable:
        chosen, status = usable[-1]["cut"], "not_converged"
    else:
        chosen, status = None, "no_usable_cutoff"
    report = {
        "status": status,
        "err_gate_applied": err_gate,
        "stability_gate_applied": any(r["stability_measured"] for r in recs),
        "min_ratio": _min_ratio,
        "criteria": {"data": "方程数/参数数 >= %g" % _min_ratio,
                     "rel_err_1se": ("train_rel_err <= min(train_rel_err) + %g*SE"
                                     "（err_kind=held-out 时才启用）" % float(se_mult)),
                     "stability": "sigma/|mean| < %g" % float(stability_thr),
                     "plateau": ("连续稳定区间：相邻 Δκ 且起点到区间内各档的累计 Δκ "
                                 "<= max(%g%%, sqrt(err_i^2+err_j^2))" % float(kappa_tol_pct))},
        "plateau_runs": [[usable[i]["cut"], usable[j]["cut"]] for i, j in sorted(runs)],
        "plateau_interval": ([usable[run_used[0]]["cut"], usable[run_used[1]]["cut"]]
                             if run_used else None),
        "after_plateau": after,
        "records": [
            {"cut": r["cut"], "ratio": r["ratio"],
             "train_rel_err": r["train_rel_err"], "train_rel_se": r["train_rel_se"],
             "err_kind": r["err_kind"],
             "stable_upper_cut": r["stable_upper_cut"],
             "kappa": r["kappa"], "kappa_err": r["kappa_err"],
             "stability_measured": r["stability_measured"],
             "data_ok": r["data_ok"], "stable_ok": r["stable_ok"],
             "err_ok": r["err_ok"], "plateau_ok": r["plateau_ok"]}
            for r in recs],
        "chosen_cut": chosen,
        "pick": str(pick).lower(),
        "kappa_tol_pct": float(kappa_tol_pct),
        "reason": ("κ 平台区间 %.2f–%.2f Å，取%s截断 = %.2f Å%s" % (
                       usable[run_used[0]]["cut"], usable[run_used[1]]["cut"],
                       "最大" if str(pick).lower() == "largest" else "最小",
                       chosen,
                       "" if err_gate else "（判据①为 in-sample 残差，已跳过否决）")
                   if status == "ok"
                   else ("κ 在数据能确定的范围内（≤ %.2f Å）仍未到平台，"
                         "报能确定范围内最大截断的 κ；更长程三阶项当前数据无法确定"
                         % (usable[-1]["cut"] if usable else float("nan"))
                         if status == "not_converged"
                         else "没有数据量/稳定性达标的截断：需补帧或加大位移")),
    }
    _low = [r for r in recs if not r["data_ok"]]
    if _low:
        report["reason"] += "；%s 方程数/参数数 < %g（数据不足，未参与选择；补帧后可纳入）" % (
            "、".join("%.2f Å（%.2f）" % (r["cut"], r["ratio"]) for r in _low),
            _min_ratio)
    if after:
        report["reason"] += "；平台之后 %s 跳变（%s），长程三阶项数据定不住，未采用" % (
            "、".join("%.2f Å" % a["cut"] for a in after),
            "、".join("%+.1f%%" % a["change_pct_vs_plateau_end"] for a in after))
    if recs and not report["stability_gate_applied"]:
        report["reason"] += ("；未做 bootstrap 重拟合（CUT3_BOOTSTRAP=0），稳定性判据②"
                             "未参与，只按 κ 平台判断")
    return chosen, report
