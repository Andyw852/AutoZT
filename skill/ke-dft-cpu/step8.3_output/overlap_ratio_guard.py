#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""overlap_ratio_guard.py —— 不依赖根因的「重叠比值」报警（VERIFICATION V23）。

【为什么要有它】重叠满足 |I|^2 <= 1，因此 unity_overlap（全部取 1）与真实重叠
（乘上真实 |I|^2）之间必有两条对**任何材料**都成立的约束：

  ADP（弹性、q 无关，整条等能线上的点对都参与）：
      mu_real / mu_unity >= 1，上界不钉死 N_v，而是**分级**（用户 2026-09-17 定）：
          [1, N_v]      绿 —— 正常
          (N_v, 2N_v]   黄 —— 放行，但要求人工确认
          > 2N_v        红 —— 拦截，不得进 zT 汇总
      N_v = 带边几个 kT 以内的**等价能谷数**（MoS2 的 K 谷是 2；六重 Q 谷会更大）。
      下界 1 是因为 |I|^2 <= 1 只会让散射率变小、迁移率变大；
      超过 N_v 是因为真实重叠不仅关掉谷间通道（unity 那边按 |I|^2=1 计入），
      连**谷内**相邻格点的 |I|^2 也小于 1（全网格实测导带均值 0.905、价带 0.912），
      所以"已修好的正确结果"也可能超过 N_v —— 这正是要分黄区而不是一刀切的原因。

  POP（非弹性，以谷内小 q 为主）：
      mu_real / mu_unity 应接近 1（本条取 +-20%）。

超出范围就标红：说明重叠这条路径给出的结果不可信，
**该结果不应进入 zT 汇总**，即使根因还没查清。

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
#   故上界取 N_v * ADP_HI_FACTOR。因子 2 的依据：谷内 |I|^2 均值约 0.9 只能解释 ~1.1x，
#   余下由大 |Δk| 通道（|I|^2 低到 0.2~0.4）贡献；取 2 留足余量，同时仍能抓住
#   病态情形（旧 h5 实测 10.6 与 30.1，全网格 2.52 判过）。
#   分级（2026-09-17 用户定）：[1, N_v] 正常；(N_v, 2N_v] **黄色**（放行但要人工确认）；
#   超过 2N_v **红色**（拦截，不得进 zT 汇总）。不再用"事后放宽单一阈值"的做法。
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
        unity = (m.group(1).lower() == "true") if m else None
        m = re.search(r"^interpolation_factor:\s*(\S+)", s, re.M)
        interp = float(m.group(1)) if m else None
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
           "json": trs[0], "mu": {}}
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


def check_pair(unity, real, nv=ADP_HI_DEFAULT):
    """返回 (verdict, lines)。verdict: 'ok' | 'red' | 'warn' | 'skip'"""
    lines = []
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
            # 分级：<=N_v 正常；(N_v, 2N_v] 黄色；> 2N_v 红色
            if ratio < ADP_LO:
                tier = "红"; red = True
                note = "**低于下界 1** —— 与 |I|^2<=1 矛盾，重叠路径必然有问题"
            elif ratio <= nv:
                tier, note = "绿", "正常"
            elif ratio <= 2.0 * nv:
                tier, note = "黄", "放行但**要求人工确认**（超出 N_v，可能是谷内 |I|^2<1 的正常效应，也可能是新问题）"
            else:
                tier, note = "红", "**拦截**（远超 2*N_v，重叠路径不可信）"
                red = True
            yellow |= (tier == "黄")   # 黄：放行但要求人工确认（用户 2026-09-17 定的分级）
            lines.append("  %-4s mu_real/mu_unity = %8.2f   [绿<=%.0f | 黄<=%.0f | 红>%.0f]  %s  %s"
                         % (mech, ratio, nv, 2.0 * nv, 2.0 * nv, tier, note))
        elif mech == "POP":
            ok = abs(ratio - 1.0) <= POP_TOL
            lines.append("  %-4s mu_real/mu_unity = %8.2f   允许 1 +-%.0f%%      %s"
                         % (mech, ratio, POP_TOL * 100, "OK" if ok else "<<< 标红"))
            red |= not ok
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
                    help="上界 = N_v * 该因子（默认 %.1f；见文件头说明）"
                         % ADP_HI_FACTOR)
    a = ap.parse_args()

    if a.material:
        runs = find_runs(a.material)
    else:
        runs = [r for r in (_read_run(d) for d in a.dirs) if r]
    unity = [r for r in runs if r["unity_overlap"] is True]
    real = [r for r in runs if r["unity_overlap"] is False]
    print("重叠比值报警（VERIFICATION V23）：找到 unity 运行 %d 个、真实重叠运行 %d 个"
          "（N_v=%.0f -> ADP 分级：绿<=%.0f | 黄<=%.0f | 红>%.0f）"
          % (len(unity), len(real), a.nv, a.nv, 2 * a.nv, 2 * a.nv))
    if not unity or not real:
        print("  只有一种模式 —— 本检查需要两种都存在才有意义；当前不触发（生产默认只跑 unity）。")
        return 0
    verdict, lines = check_pair(unity[0], real[0], a.nv)
    print("  对照：%s  vs  %s" % (os.path.basename(unity[0]["dir"]), os.path.basename(real[0]["dir"])))
    for ln in lines:
        print(ln)
    print()
    if verdict == "red":
        print("  结论：**红色** —— 真实重叠结果不可信，**拦截、不得进入 zT 汇总**。")
        return 1
    if verdict == "yellow":
        print("  结论：**黄色** —— 放行，但**要求人工确认**后再采纳")
        print("        （超出 N_v 可能是谷内 |I|^2<1 的正常效应；也可能是新的问题，需查）。")
        return 0
    print("  结论：绿色，比值在 [1, N_v] 内。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
