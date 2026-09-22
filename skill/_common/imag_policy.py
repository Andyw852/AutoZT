#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""imag_policy.py —— 虚频判据的【唯一真源】（2026-09-22）。

## 为什么有它

同一个"结构有没有虚频"此前在四处各判一次、阈值还各不相同：

| 位置 | 旧阈值 | 作用 |
|---|---|---|
| kl-dft-cpu S5  kl_fc_backends._stability_gate | IMAG_THR = 0.10 THz | 挡 S6 的门禁 |
| kl-dft-cpu S5.1 gen_step5.1_plot_phonon | IMAG_TOL = -0.05 THz（硬编码） | 只标图 |
| kl-dft-cpu S6  lattice_kappa | imag_thr = 0.5 THz | 另一套 |
| mlff  phonon_fit_driver | 硬编码 -0.10 THz | 稳定判据 |

结果互相打架（MoS₂ 实测：S5 说 stable、S5.1 说 imaginary）。本模块把它们收成一个函数，
S5 是唯一调用方/唯一裁判，其余位置读 S5 的结论。

## 判据

**近 Γ 区**：`|q|_frac < IMAG_QGAMMA`（**严格小于**；Petretto 窗口）。
★ **声学支的"有效近 Γ 窗口"= `IMAG_QGAMMA × IMAG_QGAMMA_GRACE`**（默认 1.2，用户 2026-09-22 批准）：
  目的是让**同一条近 Γ 软模**跨越 0.05 边界时**不再被升格为 fail**。
  实测背景：未对称化 MoS₂ 的 ZA 软模在 |q|∈(0,0.061) 为负，101 点路径上一个点落在
  |q|=0.0519（ν=−0.0503）—— 只比窗口大 3.8%，却因"区外 + |ν|>STRICT"把整体从 warn 升成 fail。
  `IMAG_QGAMMA_GRACE=1.0` 即退回 §二 的原始严格口径。
`|q|_frac` 的定义 = 分数坐标向量**去掉真空轴分量**后的欧氏范数（3D 用全部三个分量）。
依据：Petretto et al., Sci. Data 5, 180065 (2018) 用 `0 < |q| < 0.05` 的窗口做声子失稳标记。

**声学支**：每个 q 点的**前 3 支**（下标 0/1/2）。
★ 函数**不重排**每个 q 点的频率：负值一旦排序必然落到最前，会把"虚频光学模"误判成声学支
  （例如 Γ 点三支声学 ≈ 0、某个光学模 −0.06：排序后 −0.06 变下标 0，就再也分不出光学）。
  调用方按 phonopy 的支序给出即可（phonopy 通常已升序 ⇒ 前 3 支就是 TA/TA/LA 或 2D 的 ZA/TA/LA）。

对每个负频 `(q, 支)` 取 `|ν|`，按以下顺序归类：

1. `|ν| ≤ IMAG_THR_STRICT`（0.05 THz）                       → `noise`（不计入，verdict 仍是 pass）
2. 近 Γ 区、声学支、`|ν| ≤ IMAG_THR`（0.15 THz）            → `near_gamma_acoustic` → **warn**
3. 近 Γ 区、声学支、`|ν| > IMAG_THR`                        → `near_gamma_large` → **fail**
4. 近 Γ 区外 或 光学支，且 `|ν| > IMAG_THR_STRICT`          → `off_gamma` / `optical` → **fail**

整体 `verdict` / `imag_class` 取所有负频里**最坏**的一档。

## 阈值依据

- **IMAG_THR = 0.15 THz = 5 cm⁻¹**：对齐 Petretto 2018 的失稳标记阈值。
  （语义已从"全局阈值"改为"**近 Γ 声学支上限**"。）
- **IMAG_THR_STRICT = 0.05 THz**：数值噪声底；也是**近 Γ 区外与光学支**的上限。
  近 Γ 区外比 MP 的入库标注更严 —— 本函数是**拦 S6 的门禁**，不是入库标注。
  Lin, Poncé, Marzari, npj Comput. Mater. 8, 236 (2022) 把远离 Γ 的虚频视为真失稳特征。
- **IMAG_QGAMMA = 0.05**（分数坐标）：Petretto 2018 的近 Γ 窗口半径。
  与 `ZA_QMAX` 的关系：`ZA_QMAX` 也是分数坐标，但它是 **ZA 幂律拟合的取样上限**（默认 0.05），
  与本参数**物理含义不同**（一个定"近 Γ 邻域"，一个定"拟合窗口"），故不合并成一个参数。
"""
import math

POLICY_VERSION = "imag_policy/1.0"
THZ_TO_CM1 = 33.35641                      # 1 THz = 33.35641 cm^-1

DEFAULTS = {
    "IMAG_THR": 0.15,          # 近 Γ 声学支上限 (THz) = 5 cm^-1（Petretto 2018）
    "IMAG_THR_STRICT": 0.05,   # 噪声底；近 Γ 区外与光学支的上限 (THz)
    "IMAG_QGAMMA": 0.05,       # 近 Γ 半径（分数坐标，Petretto 2018）
    "IMAG_QGAMMA_GRACE": 1.2,  # 声学支有效近 Γ 窗口 = IMAG_QGAMMA × 它（1.0 = 退回严格口径）
}

POLICY_REFS = [
    "Petretto et al., Scientific Data 5, 180065 (2018)：0 < |q| < 0.05 失稳窗口、5 cm^-1 阈值",
    "Lin, Poncé, Marzari, npj Comput. Mater. 8, 236 (2022)：远离 Γ 的虚频是真失稳特征",
]

# 类别严重度（越大越坏）；verdict 由类别直接映射
_CLASS_RANK = {"none": 0, "noise": 1, "near_gamma_acoustic": 2,
               "near_gamma_large": 3, "off_gamma": 3, "optical": 3}
_VERDICT_OF_CLASS = {"none": "pass", "noise": "pass",
                     "near_gamma_acoustic": "warn",
                     "near_gamma_large": "fail", "off_gamma": "fail", "optical": "fail"}
# 最终 stability_verdict 的排序：pass < warn < needs_review < fail
_VERDICT_ORDER = {"pass": 0, "warn": 1, "needs_review": 2, "fail": 3}


def _get(cfg, key):
    """从 dict / StepConf / None 里取值，缺省回落到 DEFAULTS。"""
    if cfg is not None:
        try:
            v = cfg[key]
            if v is not None and str(v) != "":
                return float(v)
        except (KeyError, TypeError, ValueError):
            pass
    return float(DEFAULTS[key])


def thresholds(cfg=None):
    return {k: _get(cfg, k) for k in DEFAULTS}


def q_norm_frac(q_frac, is2d, vac_axis):
    """|q|_frac：分数坐标去掉真空轴分量后的欧氏范数（3D 用全部三个分量）。"""
    if q_frac is None:
        return None
    try:
        v = [float(x) for x in q_frac]
    except (TypeError, ValueError):
        return None
    if is2d and vac_axis is not None:
        try:
            k = int(vac_axis)
        except (TypeError, ValueError):
            k = None
        if k is not None and 0 <= k < len(v):
            v = [x for i, x in enumerate(v) if i != k]
    return math.sqrt(sum(x * x for x in v))


def classify_imag(freqs, qpoints_frac, is2d, vac_axis, cfg=None):
    """对一批 q 点的频率做统一的虚频判定。**只此一处实现**。

    参数
    ----
    freqs        : 每个 q 点的**全部支**频率 (THz，负值=虚频)；每个 q 点内部按频率升序。
                   越界/空项按跳过处理（容错，不抛）。
    qpoints_frac : 每个 q 点在原胞倒格基下的分数坐标 (长度 3)，与 freqs 一一对应。
    is2d         : 是否 2D（决定 q 范数是否去掉真空轴分量）。
    vac_axis     : 真空轴在**分数坐标/原胞基矢**里的下标（2D 必填；3D 可为 None）。
    cfg          : 含 IMAG_THR / IMAG_THR_STRICT / IMAG_QGAMMA 的配置（dict/StepConf）。

    返回
    ----
    dict，字段见模块 docstring：verdict, imag_class, min_freq_THz, min_freq_cm1,
    q_at_min_frac, q_norm_at_min, branch_at_min, n_neg_qpoints, thresholds, policy_version。
    """
    thr = _get(cfg, "IMAG_THR")
    strict = _get(cfg, "IMAG_THR_STRICT")
    qg = _get(cfg, "IMAG_QGAMMA")
    grace = _get(cfg, "IMAG_QGAMMA_GRACE") or 1.0
    qg_eff = qg * grace             # 声学支的有效近 Γ 半径（grace=1.0 时等于 Petretto 窗口）
    freqs = list(freqs or [])
    qpts = list(qpoints_frac or [])

    worst_cls, worst_abs = "none", -1.0
    n_neg_q = 0
    gmin = None                      # (freq, qi, bi, qnorm)

    for qi in range(len(freqs)):
        try:
            # ★ 不重排：支序即调用方给的顺序（见模块 docstring）。负值排序后必然在最前，
            #   会把虚频光学模误判成声学支；声学支一律取【下标 < 3】。
            fr = [float(x) for x in freqs[qi]]
        except (TypeError, ValueError):
            continue
        if not fr:
            continue
        qf = qpts[qi] if qi < len(qpts) else None
        qnorm = q_norm_frac(qf, is2d, vac_axis)
        # 声学支用带 grace 的有效窗口；光学支不受影响（见上方 docstring）
        near = (qnorm is not None and qnorm < qg_eff)
        near_strict = (qnorm is not None and qnorm < qg)
        if any(x < 0.0 for x in fr):
            n_neg_q += 1
        for bi, nu in enumerate(fr):
            if gmin is None or nu < gmin[0]:
                gmin = (nu, qi, bi, qnorm)
            if nu >= 0.0:
                continue
            a = abs(nu)
            if a <= strict:                       # 规则 1：噪声底
                cls = "noise"
            elif near and bi < 3:                 # 规则 2/3：近 Γ 声学支
                cls = "near_gamma_acoustic" if a <= thr else "near_gamma_large"
            else:                                 # 规则 4：近 Γ 区外 / 光学支
                cls = "optical" if bi >= 3 else "off_gamma"
            if (_CLASS_RANK[cls], a) > (_CLASS_RANK[worst_cls], worst_abs):
                worst_cls, worst_abs = cls, a

    q_at_min = None
    if gmin is not None and gmin[1] < len(qpts) and qpts[gmin[1]] is not None:
        try:
            q_at_min = [float(x) for x in qpts[gmin[1]]]
        except (TypeError, ValueError):
            q_at_min = None
    return {
        "verdict": _VERDICT_OF_CLASS[worst_cls],
        "imag_class": worst_cls,
        "min_freq_THz": (gmin[0] if gmin else None),
        "min_freq_cm1": (gmin[0] * THZ_TO_CM1 if gmin else None),
        "q_at_min_frac": q_at_min,
        "q_norm_at_min": (gmin[3] if gmin else None),
        "branch_at_min": (gmin[2] if gmin else None),
        "n_neg_qpoints": n_neg_q,
        "thresholds": {"IMAG_THR": thr, "IMAG_THR_STRICT": strict,
                       "IMAG_QGAMMA": qg, "IMAG_QGAMMA_GRACE": grace,
                       "IMAG_QGAMMA_EFFECTIVE": qg_eff,
                       "near_gamma": ("声学支：|q|_frac < IMAG_QGAMMA × IMAG_QGAMMA_GRACE（=%.4g）；"
                                      "Petretto 严格窗口 = |q|_frac < IMAG_QGAMMA（=%.4g）"
                                      % (qg_eff, qg)),
                       "q_norm_def": "分数坐标去掉真空轴分量后的欧氏范数（3D 用三分量）",
                       "IMAG_THR_semantics": "近 Γ 声学支上限（不是全局阈值）",
                       "IMAG_QGAMMA_GRACE_semantics": ("用户 2026-09-22 批准：同一条近 Γ 软模跨 0.05 "
                                                       "边界时不被升格为 fail；1.0 = 严格口径"),
                       "refs": list(POLICY_REFS)},
        "policy_version": POLICY_VERSION,
    }


def combine_verdict(*verdicts):
    """取最坏者；排序 pass < warn < needs_review < fail（未知值按 fail 处理）。"""
    best = "pass"
    for v in verdicts:
        if v is None:
            continue
        if _VERDICT_ORDER.get(str(v), 3) > _VERDICT_ORDER[best]:
            best = str(v)
    return best


def za_verdict(za):
    """把 kl 的 ZA 检查 dict 映射成 pass/warn/needs_review/fail。

    历史语义保持不变：没跑成（error）不拦；ok=False 才是 fail；needs_review 放行。
    """
    if not za:
        return "pass"
    if "error" in za:
        return "pass"
    if not za.get("ok"):
        return "fail"
    if za.get("za_needs_review"):
        return "needs_review"
    return "pass"


def is_stable(verdict):
    """只有 fail 挡 S6；warn 与 needs_review 都放行（保持历史行为）。"""
    return str(verdict) != "fail"


__all__ = ["POLICY_VERSION", "THZ_TO_CM1", "DEFAULTS", "POLICY_REFS",
           "classify_imag", "q_norm_frac", "thresholds",
           "combine_verdict", "za_verdict", "is_stable"]
