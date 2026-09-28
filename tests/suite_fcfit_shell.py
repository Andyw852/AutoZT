#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""suite_fcfit_shell -- fc-fit 的「壳层 / 三阶截断确定」单元自测。

覆盖（从 kl-dft-cpu 移植过来、保留同样的判据语义）：
  1. fc_common.neighbor_shells / cut3_candidates：壳层枚举 + 相邻壳层中点
  2. fc_common.resolve_cut3_candidates：auto / off / 显式列表
  3. fc_common.supercell_safe_cutoff：周期镜像安全截断
  4. fc_fit_driver.fc3_shell_index / bins / means：逐壳层 |Phi^3|
  5. fc_fit_driver.fc3_shell_stats_from_means / _recommend_cut：σ/|mean| 与
     stable_upper_cut（数据能确定的最大截断）

不碰集群、不依赖 pheasy / phono3py。
"""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "skill" / "fc-fit"))

import fc_common as fc  # noqa: E402

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
        shells, nat = fc.neighbor_shells(poscar, r_max=7.0, tol=0.05)
        check("neighbor_shells 非空且 nat=1", len(shells) >= 3 and nat == 1,
              str((shells, nat)))
        if len(shells) >= 3:
            check("最近邻壳层 = 3.0 A", abs(shells[0] - 3.0) < 1e-6, str(shells[:3]))
            check("壳层严格升序",
                  all(shells[i] < shells[i + 1] for i in range(len(shells) - 1)))
        cands, sh, note = fc.cut3_candidates(poscar, r_max=7.0, tol=0.05)
        check("cut3_candidates 有候选", len(cands) >= 2, str(cands))
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
        c_off, _, _ = fc.resolve_cut3_candidates(poscar, "off")
        check("off -> 空", c_off == [], str(c_off))
        c_none, _, _ = fc.resolve_cut3_candidates(poscar, "")
        check("空串 -> 空（不自动）", c_none == [])
        c_auto, _, _ = fc.resolve_cut3_candidates(poscar, "auto")
        check("auto -> 与 cut3_candidates 一致",
              c_auto == fc.cut3_candidates(poscar)[0])
        c_list, _, _ = fc.resolve_cut3_candidates(poscar, "4.8 5.3 5.9")
        check("显式列表按序", c_list == [4.8, 5.3, 5.9], str(c_list))


def test_safe_cutoff():
    print("[3] 超胞安全截断")
    with tempfile.TemporaryDirectory() as d:
        poscar = _write_poscar(d, a=10.0)
        lat, frac = fc.read_poscar_cell_frac(poscar)
        safe = fc.supercell_safe_cutoff(lat, margin=0.1)
        check("10 A 立方胞安全截断 = 4.9 A", abs(safe - 4.9) < 1e-6, str(safe))


def test_fc3_shell_stability():
    print("[4/5] 逐壳层 |Phi^3| 与 stable_upper_cut")
    import numpy as np
    import fc_fit_driver as drv
    cell = np.diag([3.0, 3.0, 3.0])
    frac = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0],
                     [0.0, 0.5, 0.5]])
    shell_dist, tri = drv.fc3_shell_index(cell, frac, tol=0.05, r_max=8.0)
    check("fc3_shell_index 返回壳层且升序", len(shell_dist) >= 1
          and all(shell_dist[i] <= shell_dist[i + 1]
                  for i in range(len(shell_dist) - 1)), str(shell_dist[:4]))
    check("tri 形状正确", tri.shape == (4, 4, 4), str(tri.shape))
    bins = drv.fc3_shell_bins(tri, len(shell_dist))
    check("每个壳层都有 bin 索引", len(bins) == len(shell_dist))
    rng = np.random.default_rng(0)
    base = rng.normal(size=(4, 4, 4, 3, 3, 3))
    fits = [base + 0.02 * rng.normal(size=base.shape) for _ in range(6)]
    per_means = [drv.fc3_shell_means(f, bins) for f in fits]
    stats, upper = drv.fc3_shell_stats_from_means(shell_dist, per_means,
                                                  stability_thr=0.3)
    check("低噪声 -> stable_upper_cut 有限", upper > 0 and upper != float("inf"),
          str(upper))
    check("rel_std 有定义", all(s["rel_std"] is not None for s in stats), str(stats[:2]))

    # high scatter: bootstrap fits that do not agree -> the cutoff the data can
    # determine collapses to (at most) the nearest shell
    noisy = [base + (0.5 + i) * rng.normal(size=base.shape) for i in range(6)]
    pn = [drv.fc3_shell_means(f, bins) for f in noisy]
    _, upper_n = drv.fc3_shell_stats_from_means(shell_dist, pn, stability_thr=0.3)
    check("高噪声 -> stable_upper_cut 收到最近邻附近", upper_n <= shell_dist[0] + 0.05,
          str(upper_n))

    recs = [{"cut": 4.0, "stable_upper_cut": 5.0},
            {"cut": 5.0, "stable_upper_cut": 5.0},
            {"cut": 6.0, "stable_upper_cut": 5.0}]
    check("_recommend_cut 取自身稳定范围内的最大候选", drv._recommend_cut(recs) == 5.0)
    check("_recommend_cut 全不可用时 None",
          drv._recommend_cut([{"cut": 6.0, "stable_upper_cut": 4.0}]) is None)


def main():
    test_shells_and_candidates()
    test_resolve()
    test_safe_cutoff()
    test_fc3_shell_stability()
    print("\nsuite_fcfit_shell: %s（%d 项）"
          % ("ALL PASS" if not FAIL else "FAIL", N))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
