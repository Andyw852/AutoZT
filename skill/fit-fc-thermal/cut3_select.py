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


def select_cutoff(records, kappa_tol_pct=5.0, se_mult=1.0, stability_thr=0.3,
                  pick="smallest", component_floor=0.1):
    """按三判据自动选最小可行截断。records 是每个候选截断一条 dict：

        cut               截断 (Å)
        ratio             方程数/参数数（<3 = 数据不足 → 不参与）
        train_rel_err     拟合残差（★in-sample，**不是**交叉验证）：pheasy -f 在【全部帧】
                          上的 Relative error；判据①的输入
        train_rel_se      该残差在帧 bootstrap 重抽样下的标准误（1-SE 用；缺省 0）
        err_kind          "in-sample"（默认，判据①不具否决权）| "held-out"（真留出 CV）
        stable_upper_cut  数据能确定的最大截断 (Å)：只用 σ/|mean| < stability_thr 的壳层
                          （判据②；该截断 ≤ 它才算"数据能确定"）
        kappa             300 K 面内 κ
        kappa_err         κ 的 bootstrap 误差棒（判据③；缺省 0）

    pick：三判据都满足的截断有多个时，"smallest"（默认，总则：取最小）或
          "largest"（Mg8C120 口径：平的就取能确定范围里最大的那个）。

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
    for r in recs:
        r["data_ok"] = r["ratio"] >= 3.0
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

    # 判据③：从某截断起到最外侧，相邻 Δκ ≤ max(5%·|κ|, 两误差的平方和根) 视为平台。
    # ★ 最外侧可用截断**没有更大的可比对象**，不算"平台证据" —— 那正是 user 说的
    #   "κ 到能确定的边界还在变"，必须判 not_converged 而不是假称收敛。
    plateau = set()
    for i, r in enumerate(usable):
        pairs = list(zip(usable[i:], usable[i + 1:]))
        if not pairs:
            continue
        ok = True
        for a, b in pairs:
            if not _plateau_pair(a, b, kappa_tol_pct, component_floor):
                ok = False
                break
        if ok:
            plateau.add(r["cut"])

    for r in recs:
        r["err_ok"] = r["cut"] in err_ok
        r["plateau_ok"] = r["cut"] in plateau

    inter = sorted(err_ok & plateau)
    if inter:
        if str(pick).lower() == "largest":
            # Mg8C120 口径：κ 平就取【能确定范围内】最大的截断 —— 条件是"倒数第二个
            # 可用截断也在平台上"（说明趋势一路平到能确定的边界），否则退回 plateau 起点。
            if len(usable) >= 2 and usable[-2]["cut"] in plateau:
                chosen = usable[-1]["cut"]
            else:
                chosen = inter[-1]
        else:
            chosen = inter[0]      # 总则：取满足三判据的最小截断
        status = "ok"
    elif usable:
        chosen, status = usable[-1]["cut"], "not_converged"
    else:
        chosen, status = None, "no_usable_cutoff"
    report = {
        "status": status,
        "err_gate_applied": err_gate,
        "stability_gate_applied": any(r["stability_measured"] for r in recs),
        "criteria": {"rel_err_1se": ("train_rel_err <= min(train_rel_err) + %g*SE"
                                     "（err_kind=held-out 时才启用）" % float(se_mult)),
                     "stability": "sigma/|mean| < %g" % float(stability_thr),
                     "plateau": "相邻 Δκ <= max(%g%%, sqrt(err_i^2+err_j^2))"
                                % float(kappa_tol_pct)},
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
        "reason": ("判据②③同时满足的最小截断 = %.2f Å%s" % (
                       chosen,
                       "" if err_gate else "（判据①为 in-sample 残差，已跳过否决）")
                   if status == "ok"
                   else ("κ 在数据能确定的范围内（≤ %.2f Å）仍未到平台，"
                         "报能确定范围内最大截断的 κ；更长程三阶项当前数据无法确定"
                         % (usable[-1]["cut"] if usable else float("nan"))
                         if status == "not_converged"
                         else "没有数据量/稳定性达标的截断：需补帧或加大位移")),
    }
    if recs and not report["stability_gate_applied"]:
        report["reason"] += ("；未做 bootstrap 重拟合（CUT3_BOOTSTRAP=0），稳定性判据②"
                             "未参与，只按 κ 平台判断")
    return chosen, report
