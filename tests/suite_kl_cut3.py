#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""suite_kl_cut3 —— kl 三阶截断自动选择（候选/判据/选择）的单元自测。

覆盖：
  1. neighbor_shells / cut3_candidates：壳层枚举 + 相邻壳层中点（不落在壳层距离上）
  2. resolve_cut3_candidates：auto / off / 显式列表 三种写法
  3. select_cutoff：三判据都满足取最小、平的时候 largest 口径、κ 未收敛、无可用截断
  4. fc3_shell_stability：逐壳层 Φ³ 离散度 + stable_upper_cut
不碰集群、不依赖 phono3py。
"""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "skill" / "kl-dft-cpu"))
sys.path.insert(0, str(ROOT / "skill" / "_common" / "opt"))

import kl_common as kc  # noqa: E402

FAIL = []
N = 0


def check(name, cond, extra=""):
    global N
    N += 1
    if cond:
        print("  PASS %s" % name)
    else:
        print("  FAIL %s %s" % (name, extra))
        FAIL.append(name)


def _write_poscar(d, a=3.0):
    p = Path(d) / "POSCAR"
    p.write_text(
        "Si\n1.0\n%.6f 0 0\n0 %.6f 0\n0 0 %.6f\nSi\n1\nDirect\n0 0 0\n"
        % (a, a, a), encoding="utf-8")
    return p


def test_shells_and_candidates():
    print("[1] 壳层枚举 + 候选中点")
    with tempfile.TemporaryDirectory() as d:
        poscar = _write_poscar(d)
        shells, nat = kc.neighbor_shells(poscar, r_max=7.0, tol=0.05)
        check("neighbor_shells 返回非空且 nat=1", len(shells) >= 3 and nat == 1,
              str((shells, nat)))
        if len(shells) >= 3:
            check("最近邻壳层 = 3.0 Å", abs(shells[0] - 3.0) < 1e-6, str(shells[:3]))
            check("壳层严格升序",
                  all(shells[i] < shells[i + 1] for i in range(len(shells) - 1)))
        cands, sh, note = kc.cut3_candidates(poscar, r_max=7.0, tol=0.05)
        check("cut3_candidates 有候选", len(cands) >= 2, str(cands))
        # 每个候选必是相邻两壳层的中点，且不等于任一壳层距离
        ok_mid = True
        for c in cands:
            below = [s for s in sh if s < c]
            above = [s for s in sh if s > c]
            if not below or not above:
                ok_mid = False
                break
            if abs(c - round(0.5 * (below[-1] + above[0]), 2)) > 1e-6:
                ok_mid = False
                break
            if any(abs(c - s) < 1e-6 for s in sh):
                ok_mid = False
                break
        check("候选都是相邻壳层中点、不落在壳层上", ok_mid, str((cands, sh)))


def test_resolve():
    print("[2] resolve_cut3_candidates 解析")
    with tempfile.TemporaryDirectory() as d:
        poscar = _write_poscar(d)
        c_off, _, n1 = kc.resolve_cut3_candidates(poscar, "off")
        check("off → 空", c_off == [], str(c_off))
        c_none, _, _ = kc.resolve_cut3_candidates(poscar, "")
        check("空串 → 空（不自动）", c_none == [])
        c_auto, _, n3 = kc.resolve_cut3_candidates(poscar, "auto")
        check("auto → 与 cut3_candidates 一致", c_auto == kc.cut3_candidates(poscar)[0])
        c_list, _, _ = kc.resolve_cut3_candidates(poscar, "4.8 5.3 5.9")
        check("显式列表按序", c_list == [4.8, 5.3, 5.9], str(c_list))


def _rec(cut, kappa, cv=0.10, se=0.01, upper=6.5, ratio=5.0):
    return {"cut": cut, "ratio": ratio, "cv_err": cv, "cv_se": se,
            "stable_upper_cut": upper, "kappa": kappa, "kappa_err": 0.05 * abs(kappa)}


def test_select():
    print("[3] select_cutoff 三判据")
    flat = [_rec(4.8, 10.0), _rec(5.3, 10.1), _rec(5.9, 10.2)]
    c1, r1 = kc.select_cutoff(flat, kappa_tol_pct=5.0, pick="smallest")
    check("平坦三档 smallest → 4.8", abs(c1 - 4.8) < 1e-9 and r1["status"] == "ok",
          str((c1, r1["status"])))
    c2, r2 = kc.select_cutoff(flat, kappa_tol_pct=5.0, pick="largest")
    check("平坦三档 largest → 5.9（能确定范围内最大）", abs(c2 - 5.9) < 1e-9,
          str(c2))
    drop = [_rec(4.8, 10.0), _rec(5.3, 12.0), _rec(5.9, 13.0)]
    c3, r3 = kc.select_cutoff(drop, kappa_tol_pct=5.0)
    check("κ 还在涨 → not_converged 且报最大可用", r3["status"] == "not_converged"
          and abs(c3 - 5.9) < 1e-9, str((c3, r3["status"])))
    bad = [_rec(4.8, 10.0, ratio=2.0), _rec(5.3, 10.0, ratio=2.5)]
    c4, r4 = kc.select_cutoff(bad)
    check("数据量全不足 → no_usable_cutoff", c4 is None
          and r4["status"] == "no_usable_cutoff", str(r4["status"]))
    unstable = [_rec(4.8, 10.0, upper=4.2), _rec(5.3, 10.0, upper=4.2)]
    c5, r5 = kc.select_cutoff(unstable)
    check("稳定性上限低于候选 → no_usable_cutoff", c5 is None, str(c5))


def test_fc3_stability():
    print("[4] fc3_shell_stability 逐壳层离散度")
    import numpy as np
    cell = np.diag([3.0, 3.0, 3.0])
    frac = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5]])
    rng = np.random.default_rng(0)
    base = rng.normal(size=(4, 4, 4, 3, 3, 3))
    fits = [base + 0.02 * rng.normal(size=base.shape) for _ in range(6)]
    stats, upper = kc.fc3_shell_stability(cell, frac, fits, tol=0.05,
                                          stability_thr=0.3, r_max=8.0)
    check("返回壳层统计且 shell 升序", len(stats) >= 1
          and all(stats[i]["shell"] <= stats[i + 1]["shell"]
                  for i in range(len(stats) - 1)), str(stats[:3]))
    check("低噪声 → stable_upper_cut 有限且 ≥ 最小壳层",
          upper > 0 and (not stats or upper >= stats[0]["shell"] - 0.05),
          str(upper))
    check("每次拟合均值都算进去了（rel_std 有定义）",
          all(s["rel_std"] is not None for s in stats), str(stats))


def main():
    test_shells_and_candidates()
    test_resolve()
    test_select()
    test_fc3_stability()
    print("\nsuite_kl_cut3: %s（%d 项）" % ("ALL PASS" if not FAIL else "FAIL", N))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
