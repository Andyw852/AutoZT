#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_dp_valley_probe.py —— V141：tools/dp_valley_probe.py 的谷拆分与判读（纯 numpy，不需要 AMSET / 数据）。

用法：python test_dp_valley_probe.py      退出码 0 = 全部 PASS。

合成六方单层导带：K/K' 谷在带边，D_iso = 7.4 eV；Q 谷（6 个）在其上 ΔE，D_iso = 1.0 eV。
  · ΔE = 30 meV：Q 谷约占一半载流子，而 μ ∝ 1/D² 并联 -> 倍数 E1²·⟨1/D²⟩ ≈ 27（MoSe2 电子实测 28.5 的量级），
    ADP-only 迁移率几乎全来自 Q 谷；只看 rms D（≈ 5.4）会严重低估；
  · ΔE = 400 meV：Q 谷在窗口外 -> 倍数 = 1。
另测 −k 匹配与空窗口。
"""
import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "tools"))
import dp_valley_probe as P  # noqa: E402

N = 24
A = 3.3
LAT = np.array([[A, 0, 0], [-A / 2, A * np.sqrt(3) / 2, 0], [0, 0, 20.0]])
RECIP = 2 * np.pi * np.linalg.inv(LAT).T
K_PTS = [np.array([1 / 3, 1 / 3, 0]), np.array([2 / 3, 2 / 3, 0])]
Q_PTS = [np.array(q + [0]) for q in ([1 / 6, 1 / 6], [5 / 6, 5 / 6], [1 / 6, 2 / 3], [2 / 3, 1 / 6],
                                     [5 / 6, 1 / 3], [1 / 3, 5 / 6])]


def _grid():
    g = np.array([[i / N, j / N, 0.0] for i in range(N) for j in range(N)])
    return g


def _band(dq_ev, m=0.5):
    """一条导带：K 谷与 Q 谷都是抛物线（有效质量 m），返回能量与 D 张量（只填 xx/yy）。"""
    kp = _grid()
    hb2m = 3.81 / m                      # ħ²/2m，eV·Å²
    E = np.full(len(kp), 10.0)
    D = np.zeros((len(kp), 3, 3))
    for centers, e0, d in ((K_PTS, 0.0, 7.4), (Q_PTS, dq_ev, 1.0)):
        for c in centers:
            dk = kp - c
            dk -= np.rint(dk)
            e = e0 + hb2m * np.linalg.norm(dk @ RECIP, axis=1) ** 2
            better = e < E
            E[better] = e[better]
            D[better, 0, 0] = d
            D[better, 1, 1] = d
    return kp, E[None, :], D[None, :]


class ValleyTests(unittest.TestCase):
    def test_floor(self):
        kp, E, D = _band(0.4)
        D = D.copy()
        D[0, :, 0, 0] = D[0, :, 1, 1] = 0.0                          # D 全 0：不发散，计数
        r = P.valley_stats(E, D, kp, 0.0, "electron", recip=RECIP)
        self.assertEqual(r["mu_factor"], 1.0)
        self.assertEqual(r["n_below_floor"], r["n_states"])

    def test_low_q_valley_dominates(self):
        kp, E, D = _band(0.030)
        r = P.valley_stats(E, D, kp, 0.0, "electron", T=300.0, window=0.3, recip=RECIP)
        self.assertAlmostEqual(r["edge_D"], 7.4, places=3)
        q = [v for v in r["valleys"] if abs(v["rms_D"] - 1.0) < 1e-6]
        k = [v for v in r["valleys"] if abs(v["rms_D"] - 7.4) < 1e-6]
        self.assertEqual((len(q), len(k)), (6, 2))
        wq = sum(v["weight"] for v in q)
        self.assertTrue(0.3 < wq < 0.7, wq)
        # 解析值：倍数 = 7.4² · (wK/7.4² + wQ/1²)
        self.assertAlmostEqual(r["mu_factor"], 1 - wq + wq * 7.4 ** 2, delta=0.05)
        self.assertGreater(r["mu_factor"], 15.0)
        self.assertGreater(sum(v["mu_share"] for v in q), 0.95)        # ADP-only 迁移率几乎全来自 Q 谷
        self.assertGreater(r["rms_D"], 4.0)                           # rms 看不出来
        self.assertAlmostEqual(min(v["dE_meV"] for v in q), 30.0, delta=15.0)

    def test_far_q_valley_no_effect(self):
        kp, E, D = _band(0.400)
        r = P.valley_stats(E, D, kp, 0.0, "electron", T=300.0, window=0.3, recip=RECIP)
        self.assertAlmostEqual(r["rms_D"], 7.4, places=3)
        self.assertAlmostEqual(r["mu_factor"], 1.0, places=3)
        self.assertEqual(len(r["valleys"]), 2)                      # 只有 K、K'

    def test_hole_side_and_empty(self):
        kp, E, D = _band(0.4)
        r = P.valley_stats(-E, D, kp, 0.0, "hole", window=0.3, recip=RECIP)   # 价带：能量向下
        self.assertAlmostEqual(r["edge_D"], 7.4, places=3)
        self.assertIsNone(P.valley_stats(E + 5.0, D, kp, 0.0, "electron", window=0.3))

    def test_match_kpoints_time_reversal(self):
        kbs = np.array([[0, 0, 0], [0.25, 0, 0], [1 / 3, 1 / 3, 0]])
        kh5 = np.array([[1.0, 0, 0], [-0.25, 0, 0], [2 / 3, 2 / 3, 0], [0.5, 0.5, 0]])
        self.assertEqual(P.match_kpoints(kh5, kbs).tolist(), [0, 1, 2, -1])

    def test_report_and_missing_files(self):
        kp, E, D = _band(0.03)
        r = P.valley_stats(E, D, kp, 0.0, "electron", recip=RECIP)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            P._report(r)
        self.assertIn("E1²·⟨1/D²⟩", buf.getvalue())
        self.assertIn("迁移率份", buf.getvalue())
        self.assertIn("偏大约 5 倍", buf.getvalue())                    # V142：绝对值不可直接当结论
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(P.main([tempfile.mkdtemp()]), 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
