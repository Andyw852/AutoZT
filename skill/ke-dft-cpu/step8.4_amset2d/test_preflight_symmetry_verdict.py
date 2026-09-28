#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_preflight_symmetry_verdict.py —— overlap_preflight 的**统一判据**接线自检。

用法：python test_preflight_symmetry_verdict.py     退出码 0 = 全部 PASS。

背景（2026-09-28 第二批 (a)）：preflight 的 ② 原本只按"h5 是否完整网格"一刀切，
会把**原点已对齐**（τ≡0）的 2D 体系也拦下 —— 那种情况走 IBZ 是正确的（V109/V111）。
现在它调 ke_common.symmetry_gate()；本测试固定这条接线：
  · GaAs（无反演、τ≡0）-> 放行；
  · GaN （无反演、τ=0.5）-> 拦截；
  · 读不到结构        -> 保守拦截；
  · 退回实现 _fallback_verdict 与 ke_common 结论一致。
"""
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT), str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import overlap_preflight as pf

RES = {
    "GaAs": "/mnt/d/tf_data/work_AutoZT/GaAs/ke-dft-cpu/result",
    "GaN": "/mnt/d/tf_data/work_AutoZT/GaN/ke-dft-cpu/result",
}


class PreflightVerdictTests(unittest.TestCase):
    def _verdict(self, name):
        d = Path(RES[name])
        if not d.is_dir():
            self.skipTest("缺 %s" % d)
        return pf._structure_verdict(d)

    def test_gaas_tau_zero_released(self):
        need, why = self._verdict("GaAs")
        self.assertFalse(need, why)
        self.assertIn("ke_common", why)          # 走的是唯一来源，不是 fallback

    def test_gan_requires_full_grid(self):
        need, why = self._verdict("GaN")
        self.assertTrue(need, why)
        self.assertIn("ke_common", why)

    def test_missing_structure_is_conservative(self):
        with tempfile.TemporaryDirectory() as tmp:
            need, why = pf._structure_verdict(Path(tmp))
        self.assertTrue(need, why)
        self.assertIn("读不到结构", why)

    def test_fallback_matches_ke_common(self):
        from pymatgen.core import Structure
        for name, expect in (("GaAs", False), ("GaN", True)):
            p = Path(RES[name]) / "step3_uniform" / "POSCAR"
            if not p.is_file():
                continue
            st = Structure.from_file(str(p))
            need, why = pf._fallback_verdict(st)
            self.assertEqual(need, expect, "%s: %s" % (name, why))
            self.assertIn("fallback", why)


if __name__ == "__main__":
    unittest.main(verbosity=2)
