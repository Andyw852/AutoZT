#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_preflight_symmetry_verdict.py —— overlap_preflight 的**统一判据**接线自检。

用法：python test_preflight_symmetry_verdict.py     退出码 0 = 全部 PASS。

背景：preflight 的 ② 原本只按"h5 是否完整网格"一刀切，会把原点已对齐（或本来就精确）
的体系也拦下；现在它调 ke_common.symmetry_gate()。本测试固定这条接线，**结构在测试里
构造**（不读 /mnt/d 下的任何文件）。
"""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.inputs import Poscar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT), str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import overlap_preflight as pf


def shifted(st, s):
    x = st.copy()
    x.translate_sites(range(len(x)), s, frac_coords=True)
    return x


def gan_std():
    return Structure(Lattice.hexagonal(3.19, 5.19), ["Ga", "Ga", "N", "N"],
                     [[1/3, 2/3, 0], [2/3, 1/3, .5], [1/3, 2/3, .377], [2/3, 1/3, .877]])


def si_on_atom():
    return Structure(Lattice([[0, 2.715, 2.715], [2.715, 0, 2.715], [2.715, 2.715, 0]]),
                     ["Si", "Si"], [[0, 0, 0], [.25, .25, .25]])


def _material_dir(st):
    """把结构写成 <tmp>/step3_uniform/POSCAR，返回材料目录。"""
    d = Path(tempfile.mkdtemp())
    (d / "step3_uniform").mkdir(parents=True)
    Poscar(st).write_file(str(d / "step3_uniform" / "POSCAR"))
    return d


class PreflightVerdictTests(unittest.TestCase):
    def test_gan_std_released_via_ke_common(self):
        """GaN 标准原点：原公式 0/24 错 -> 放行，且**确实走 ke_common**。"""
        need, why = pf._structure_verdict(_material_dir(gan_std()))
        self.assertFalse(need, why)
        self.assertIn("ke_common", why)

    def test_si_at_inversion_centre_requires_full_grid(self):
        """★ 关键：有反演但原点在反演中心的 Si 原公式错 24/48 —— 必须拦。"""
        st = shifted(si_on_atom(), [-.125, -.125, -.125])
        need, why = pf._structure_verdict(_material_dir(st))
        self.assertTrue(need, why)
        self.assertIn("ke_common", why)

    def test_missing_structure_is_conservative(self):
        with tempfile.TemporaryDirectory() as tmp:
            need, why = pf._structure_verdict(Path(tmp))
        self.assertTrue(need, why)
        self.assertIn("读不到结构", why)

    def test_fallback_matches_ke_common(self):
        cases = [("GaN std", gan_std(), False),
                 ("Si inv centre", shifted(si_on_atom(), [-.125] * 3), True)]
        for name, st, expect in cases:
            need, why = pf._fallback_verdict(st)
            self.assertEqual(need, expect, "%s: %s" % (name, why))
            self.assertIn("fallback", why)


if __name__ == "__main__":
    unittest.main(verbosity=2)
