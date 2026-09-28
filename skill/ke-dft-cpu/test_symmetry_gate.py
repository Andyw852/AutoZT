#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_symmetry_gate.py —— ke_common.symmetry_gate 的自检（**不依赖任何本地数据**）。

用法：python test_symmetry_gate.py     退出码 0 = 全部 PASS。

★ 2026-09-28 换判据：旧规则"有反演或 τ≡0 就可走 IBZ"两个方向都错（见 V115）。
新判据是**逐操作**的：非 TR 操作要 (R+I)τ≡0，TR 操作要 (R−I)τ≡0（mod 1）。
下面的期望值全部来自独立检验 tmp/amset_desym_phase_check.py 的逐操作实测
（模型势直接对角化 + 与 AMSET 去对称化结果比较），12 个构型 100% 吻合。

注意：结构**在测试里构造**，不读 /mnt/d 下的任何文件（那些路径会被推到公开仓库）。
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
import ke_common as kc


def shifted(st, s):
    x = st.copy()
    x.translate_sites(range(len(x)), s, frac_coords=True)
    return x


def gan_std():
    """GaN P6_3mc，6_3 轴过原点（标准设置）—— 原公式 0/24 错，IBZ 精确。"""
    return Structure(Lattice.hexagonal(3.19, 5.19), ["Ga", "Ga", "N", "N"],
                     [[1/3, 2/3, 0], [2/3, 1/3, .5], [1/3, 2/3, .377], [2/3, 1/3, .877]])


def mos2_aligned():
    return Structure(Lattice.hexagonal(3.19, 12.0), ["Mo", "S", "S"],
                     [[0, 0, .5], [1/3, 2/3, .63], [1/3, 2/3, .37]])


def si_on_atom():
    return Structure(Lattice([[0, 2.715, 2.715], [2.715, 0, 2.715], [2.715, 2.715, 0]]),
                     ["Si", "Si"], [[0, 0, 0], [.25, .25, .25]])


def hcp_std():
    return Structure(Lattice.hexagonal(3.21, 5.21), ["Mg", "Mg"],
                     [[1/3, 2/3, .25], [2/3, 1/3, .75]])


def zincblende_std():
    return Structure(Lattice([[0, 2.825, 2.825], [2.825, 0, 2.825], [2.825, 2.825, 0]]),
                     ["Ga", "N"], [[0, 0, 0], [.25, .25, .25]])


# (名字, 结构, bad_ops 期望, needs_full_grid 期望) —— bad_ops 来自独立实测
CASES = [
    ("GaN 标准原点",       gan_std(),                          0,  False),
    ("GaN 原点在 Ga 上",   shifted(gan_std(), [-1/3, -2/3, 0]), 6,  True),
    ("GaN 任意平移",       shifted(gan_std(), [.13, .07, .2]), 15,  True),
    ("MoS2 对齐后",        mos2_aligned(),                     0,  False),
    ("MoS2 对齐前",        shifted(mos2_aligned(), [.1, .2, 0]), 12, True),
    ("Si 原点在原子",      si_on_atom(),                       0,  False),
    # ★ 关键回归：有反演的 Si，原点落在反演中心时**原公式错 24/48**，旧规则却放行 IBZ
    ("Si 原点在反演中心",  shifted(si_on_atom(), [-.125] * 3), 24,  True),
    ("Si 任意平移",        shifted(si_on_atom(), [.11, .05, .19]), 28, True),
    ("hcp 标准原点",       hcp_std(),                          0,  False),
    # ★ 关键回归：有反演但任意平移 -> 错 8/24，旧规则同样会放行
    ("hcp 任意平移(有反演)", shifted(hcp_std(), [.12, .31, .07]), 8, True),
    ("闪锌矿 标准原点",    zincblende_std(),                   0,  False),
    ("闪锌矿 任意平移",    shifted(zincblende_std(), [.11, .05, .19]), 37, True),
]


class SymmetryGateTests(unittest.TestCase):
    def test_ground_truth_table(self):
        """12 个构型的 bad_ops 必须与独立检验逐一对上。"""
        for name, st, bad, need in CASES:
            g = kc.symmetry_gate(st)
            self.assertEqual(g["bad_ops"], bad, "%s: %s" % (name, g))
            self.assertEqual(g["needs_full_grid"], need, "%s: %s" % (name, g))

    def test_none_structure_is_conservative(self):
        g = kc.symmetry_gate(None)
        self.assertTrue(g["needs_full_grid"])
        self.assertIsNone(g["bad_ops"])
    
    def test_op_phase_exact_helper(self):
        I = np.eye(3)
        # 非 TR：R = -I 时 (R+I)=0 -> 对任意 tau 都满足
        self.assertTrue(kc.op_phase_exact(-I, np.array([0.3, 0.0, 0.0]), False))
        # 非 TR：R = I 时 (R+I)=2I -> tau=0.3 给出 0.6 mod 1 -> 不满足
        self.assertFalse(kc.op_phase_exact(I, np.array([0.3, 0.0, 0.0]), False))
        self.assertTrue(kc.op_phase_exact(I, np.array([0.5, 0.0, 0.0]), False))
        # TR：R = I 时 (R-I)=0 -> 恒满足
        self.assertTrue(kc.op_phase_exact(I, np.array([0.137, 0.2, 0.31]), True))
        # TR：R = -I 时 (R-I)=-2I -> tau=0.3 -> -0.6 mod 1 不满足
        self.assertFalse(kc.op_phase_exact(-I, np.array([0.3, 0.0, 0.0]), True))

    def test_write_full_grid_marker(self):
        """S3 用它落 full_grid_needed.json —— 键与结论都要对。"""
        for name, st, bad, need in (("GaN标准", gan_std(), 0, False),
                                    ("GaN平移", shifted(gan_std(), [.13, .07, .2]), 15, True)):
            d = Path(tempfile.mkdtemp())
            from pymatgen.io.vasp.inputs import Poscar
            p = d / "POSCAR"
            Poscar(st).write_file(str(p))
            g = kc.write_full_grid_marker(p, d)
            self.assertIsNotNone(g, name)
            self.assertEqual(g["needs_full_grid"], need, g)
            f = d / "full_grid_needed.json"
            self.assertTrue(f.is_file(), name)
            j = json.loads(f.read_text(encoding="utf-8"))
            self.assertEqual(j["needs_full_grid"], need, j)
            self.assertEqual(j["bad_ops"], bad, j)
            self.assertIn("reason", j)

    def test_marker_returns_none_on_bad_input(self):
        self.assertIsNone(kc.write_full_grid_marker(Path("/nonexistent/POSCAR")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
