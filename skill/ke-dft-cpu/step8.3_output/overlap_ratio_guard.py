#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""overlap_ratio_guard.py —— 不依赖根因的「重叠比值」报警（VERIFICATION V23）。

【为什么要有它】重叠满足 |I|^2 <= 1，因此 unity_overlap（全部取 1）与真实重叠
（乘上真实 |I|^2）之间必有两条对**任何材料**都成立的约束：

  ADP（弹性、q 无关，整条等能线上的点对都参与）：
      mu_real / mu_unity >= 1，上界**不是 N_v**，而是
      **N_ch = "unity 计入、而真实重叠不计入的所有通道之比"**（2026-09-18 修正，V23.5）：
          [1, N_v]        绿 —— 正常
          (N_v, 上界]     黄 —— 放行，但要求人工确认
          > 上界          红 —— 拦截，不得进 zT 汇总
      上界的取法：**有实测 N_ch 就用它**（`--valley-ratio`，来自"按谷区分"测试：
      同谷 |I|^2=1、谷间=0，其余输入不变）；没有才退回 2*N_v 作兜底代理。
      为什么 N_v 当上界有误：unity 不只多算谷间那一项，还按 |I|^2=1 计入了所有
      带间 / 跨 kz 通道，而真实重叠里这些都是被压低的 -> 比值可以超过 N_v。
      MoS2 实测：真实重叠/unity = 2.52，按谷区分/unity = **3.13** -> 2.52 < 3.13 判通过；
      **N_v 只是 N_ch 的下界**。下界 1 是因为 |I|^2 <= 1 只会让散射率变小、迁移率变大。

  POP（非弹性，以谷内小 q 为主）：
      **只记录比值，不参与红/黄判定**（2026-09-19 用户修正）：MoS2 实测真实重叠连谷内
      通道也压低（|I|^2 均值约 0.9），期望本就不该贴近 1；用 +-20% 会把正常结果判红。

分级（2026-09-19 用户修正）：
    ADP  ratio < 1           红 —— 与 |I|^2<=1 矛盾，拦截
    ADP  1 <= ratio <= N_v   绿
    ADP  ratio > N_v         黄 —— 要求人工确认；实测 N_ch 仅作解释，**不作放行阈值**
                               （SS 实测 N_ch=23.3 太大，用它放行等于没有红线）
    POP / 其它机制           只记录
**红线只有 ADP 与“比值不得小于 1”这两条。**

用法：
  python overlap_ratio_guard.py <运行目录1> <运行目录2> ...     # 各自含 transport*.json + settings.yaml
  python overlap_ratio_guard.py --material <材料result目录>      # 自动找 step8_amset / step8.4_amset2d
"""
import argparse
import glob
import json
import os
import re
import sys

ADP_LO = 1.0
# ★ 上界口径（2026-09-17 修正，VERIFICATION V23.5）：
#   用户原先给的是 [1, N_v]。但标准答案重叠（全网格 h5）实测 —— **同谷**相邻格点的
#   |I|^2 也不是 1（导带均值 0.905、价带 0.912，最差低到 0.04），所以真实重叠不仅关掉
#   谷间通道，也**压低谷内**通道 -> 比值会超过 N_v。全网格 + 真实重叠的实测比值是 2.52，
#   而 N_v = 2 —— 用 [1, N_v] 会把"已经修好的正确结果"误判成标红。
#   真正的上界是 N_ch（= unity 计入而真实重叠不计入的所有通道之比），只能实测（按谷区分测试）。
#   **N_v 只是它的下界**：unity 还按 |I|^2=1 计入了所有带间/跨 kz 通道。
#   本常数给出**兜底代理**上界 = N_v * ADP_HI_FACTOR（没做按谷区分测试时用）；
#   一旦有实测 N_ch（--valley-ratio），优先用它（见 check_pair）。
#   分级（2026-09-18 用户修正）：[1, N_v] 绿；(N_v, 上界] **黄**（人工确认）；
#   超过上界 **红**（拦截，不得进 zT 汇总）。
ADP_HI_FACTOR = 2.0
YELLOW_FACTOR = 1.0                       # 黄区上界 = N_v * (1 + YELLOW_FACTOR) 的保守写法见 check_pair
ADP_HI_DEFAULT = 2.0 * ADP_HI_FACTOR      # 默认按 N_v = 2 算；调用方可传 --nv 覆盖
POP_TOL = 0.20
DOPING_REF = 1e17


def _read_run(d):
    """读一个运行目录的 transport*.json + settings.yaml，返回摘要 dict 或 None。"""
    trs = sorted(glob.glob(os.path.join(d, "transport*.json")))
    if not trs:
        return None
    st = os.path.join(d, "settings.yaml")
    unity, interp = None, None
    if os.path.isfile(st):
        s = open(st, errors="ignore").read()
        m = re.search(r"^unity_overlap:\s*(\S+)", s, re.M)
        # ★ 缺这一行时按 AMSET 0.4.19 的默认值 False（真实重叠）处理 —— 旧 gen 不写这一行，
        #   若记成 None，8.3 / CLI 会把它当"未知"跳过 -> 三维黄色、二维红色都不会触发。
        #   （settings.yaml 整个缺失时保持 None = 未知，见上面的初值。）
        unity = (m.group(1).lower() == "true") if m else False
        m = re.search(r"^interpolation_factor:\s*(\S+)", s, re.M)
        interp = float(m.group(1)) if m else None
        m = re.search(r"^use_projections:\s*(\S+)", s, re.M)
        use_proj = (m.group(1).lower() == "true") if m else False
    try:
        j = json.load(open(trs[0]))
    except Exception:
        return None
    dop = j.get("doping")
    temps = j.get("temperatures")
    if not dop or not temps:
        return None
    import numpy as np
    dop = np.array(dop)
    i_t = int(np.argmin(np.abs(np.array(temps) - 300.0)))
    i_n = int(np.argmin(np.abs(dop + DOPING_REF)))     # 负掺杂 = n 型（电子）
    out = {"dir": d, "unity_overlap": unity, "interpolation_factor": interp,
           "use_projections": use_proj, "json": trs[0], "mu": {}}
    for mech, arr in (j.get("mobility") or {}).items():
        try:
            a = np.array(arr)
            out["mu"][mech] = float(a[i_n, i_t, 0, 0])
        except Exception:
            pass
    return out


def find_runs(material_dir):
    cands = []
    for pat in ("step8_amset", "step8.4_amset2d", "step8*", "*amset*", "*wave*"):
        cands += glob.glob(os.path.join(material_dir, pat))
    runs = []
    for d in sorted(set(cands)):
        if os.path.isdir(d):
            r = _read_run(d)
            if r:
                runs.append(r)
    return runs


def check_pair(unity, real, nv=ADP_HI_DEFAULT, valley_ratio=None):
    """返回 (verdict, lines)。verdict: 'ok' | 'red' | 'yellow' | 'skip'

    ★ 2026-09-19 用户修正的分级：
      ADP: ratio < 1            -> 红（与 |I|^2<=1 矛盾）
      ADP: 1 <= ratio <= N_v    -> 绿
      ADP: ratio > N_v          -> 黄（**要求人工确认**）；实测 N_ch 只在黄里作为解释打印，
                                   不再当放行阈值（SS 实测 N_ch=23.3 太大，用它放行等于没红线）
      POP / 其它                -> 只记录，不参与红/黄
    红线只有 ADP 与“比值不得小于 1”。
    """
    lines = []
    # use_projections（轨道投影）：|I|^2 <= 1 的先验**不成立**（实测 mu_real/mu_unity = 0.71/0.41 < 1，
    # 见 VERIFICATION V26）—— 它是另一种近似，不是"同一个重叠的粗略版"。比值报警对它不适用。
    if unity.get("use_projections") or real.get("use_projections"):
        return "skip", ["  使用了 use_projections（轨道投影）：|I|^2<=1 的先验不成立，比值报警不适用"
                        "（实测可 <1）—— 跳过。"]
    red = False
    yellow = False
    checked = 0
    for mech in sorted(set(unity["mu"]) & set(real["mu"])):
        u, r = unity["mu"][mech], real["mu"][mech]
        if not u or u <= 0:
            continue
        ratio = r / u
        checked += 1
        if mech == "ADP":
            if ratio < ADP_LO:
                tier = "红"; red = True
                note = "**低于下界 1** —— 与 |I|^2<=1 矛盾，重叠路径必然有问题（拦截）"
            elif ratio <= nv:
                tier, note = "绿", "正常（在带边能谷数 N_v 内）"
            else:
                tier = "黄"; yellow = True
                if valley_ratio:
                    note = ("超出 N_v，**要求人工确认**；按谷区分实测 N_ch=%.2f 可作解释"
                            "（不自动放行）" % float(valley_ratio))
                else:
                    note = "超出 N_v，**要求人工确认**（建议做一次按谷区分测试作解释）"
            lines.append("  %-4s mu_real/mu_unity = %8.2f   [绿<=%.1f | >%.1f 黄(人工确认)]  %s  %s"
                         % (mech, ratio, nv, nv, tier, note))
        elif mech == "POP":
            lines.append("  %-4s mu_real/mu_unity = %8.2f   （只记录，不参与红/黄；真实重叠也压低谷内通道）"
                         % (mech, ratio))
        else:
            lines.append("  %-4s mu_real/mu_unity = %8.2f   （无先验范围，仅记录）" % (mech, ratio))
    if checked == 0:
        return "skip", ["  两个运行没有可比的机制，跳过"]
    return ("red" if red else ("yellow" if yellow else "ok")), lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--material", help="材料 result 目录（自动找 amset 运行）")
    ap.add_argument("--nv", type=float, default=2.0,
                    help="带边等价能谷数 N_v（MoS2 K 谷 = 2；Q 谷六重则取 6）")
    ap.add_argument("--hi-factor", type=float, default=ADP_HI_FACTOR,
                    help="已弃用（2026-09-19 分级修正后不再用乘性上界；仅为兼容保留）")
    ap.add_argument("--valley-ratio", type=float, default=None,
                    help="按谷区分测试实测的 N_ch —— 仅作解释打印，不作放行阈值（见 VERIFICATION V27）")
    a = ap.parse_args()

    if a.material:
        runs = find_runs(a.material)
    else:
        runs = [r for r in (_read_run(d) for d in a.dirs) if r]
    unity = [r for r in runs if r["unity_overlap"] is True]
    real = [r for r in runs if r["unity_overlap"] is False]
    _nch = ("；实测 N_ch=%.2f（仅作解释，不作放行阈值）" % a.valley_ratio) if a.valley_ratio else ""
    print("重叠比值报警（VERIFICATION V23.5，2026-09-19 分级修正）：找到 unity 运行 %d 个、"
          "真实重叠运行 %d 个（N_v=%.0f -> ADP：绿 [1,%.0f]；>%.0f 黄（要求人工确认）；<1 红；"
          "POP 只记录%s）" % (len(unity), len(real), a.nv, a.nv, a.nv, _nch))
    if not unity or not real:
        print("  只有一种模式 —— 本检查需要两种都存在才有意义；当前不触发（生产默认只跑 unity）。")
        return 0
    verdict, lines = check_pair(unity[0], real[0], a.nv, a.valley_ratio)
    print("  对照：%s  vs  %s" % (os.path.basename(unity[0]["dir"]), os.path.basename(real[0]["dir"])))
    for ln in lines:
        print(ln)
    print()
    if verdict == "red":
        print("  结论：**红色** —— ADP 比值 < 1，与 |I|^2<=1 矛盾，**拦截、不得进入 zT 汇总**。")
        return 1
    if verdict == "yellow":
        print("  结论：**黄色** —— ADP 比值超过 N_v，放行但**要求人工确认**后再采纳")
        print("        实测 N_ch 仅作解释（不自动放行）；人工确认后再进入 zT 汇总。")
        return 0
    if verdict == "skip":
        print("  结论：**跳过**（见上；use_projections / 无可比机制时不适用本报警）。")
        return 0
    print("  结论：绿色（ADP 比值在 [1, N_v] 内）。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
