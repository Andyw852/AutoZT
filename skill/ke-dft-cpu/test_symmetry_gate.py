#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_symmetry_gate.py —— ke_common.symmetry_gate 的自检（用真实结构，只读）。

用法：python test_symmetry_gate.py     退出码 0 = 全部 PASS。
覆盖第二批 (a) 要统一的那条判据：
  · GaAs F-43m  ：有反演、|tau|=0        -> IBZ 即可
  · Si   Fd-3m  ：有反演、|tau|=0.25     -> IBZ 即可（TR 被点群吸收）
  · GaN  P6_3mc ：无反演、|tau|=0.5      -> **必须全网格**（非简单空间群，平移对不齐）
  · MoS2 P-6m2 原点已对齐：无反演、|tau|≈0 -> IBZ 即可（V109/V111）
  · 判不出来（structure=None）          -> 保守按需要全网格
"""
import sys
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
import ke_common as kc

STRUCTS = {
    "GaAs": "/mnt/d/tf_data/work_AutoZT/GaAs/ke-dft-cpu/result/step3_uniform/POSCAR",
    "Si": "/mnt/d/tf_data/work_AutoZT/Si_ke34/ke-dft-cpu/result/step3_uniform/POSCAR",
    "GaN": "/mnt/d/tf_data/work_AutoZT/GaN/ke-dft-cpu/result/step3_uniform/POSCAR",
    "MoS2": "/mnt/d/tf_data/work_AutoZT/MoS2/ke-dft-cpu/result/step3_uniform/POSCAR",
}


def _load(name):
    from pymatgen.core import Structure
    return Structure.from_file(STRUCTS[name])


class SymmetryGateTests(unittest.TestCase):
    def _gate(self, name):
        p = Path(STRUCTS[name])
        if not p.is_file():
            self.skipTest("缺 %s" % p)
        return kc.symmetry_gate(_load(name))

    def test_none_structure_is_conservative(self):
        g = kc.symmetry_gate(None)
        self.assertTrue(g["needs_full_grid"])
        self.assertIsNone(g["max_tau"])

    def test_gaas_no_inversion_but_tau_zero_ibz_ok(self):
        """闪锌矿 F-43m **没有**反演中心，但全部 |tau|=0 -> 去对称化照样对，IBZ 即可。"""
        g = self._gate("GaAs")
        self.assertFalse(g["has_inversion"], g)
        self.assertFalse(g["needs_full_grid"], g)
        self.assertEqual(g["tau_ops"], 0, g)

    def test_si_inversion_ibz_ok(self):
        g = self._gate("Si")
        self.assertTrue(g["has_inversion"], g)
        self.assertFalse(g["needs_full_grid"], g)
        self.assertGreater(g["max_tau"], 0.1, g)      # Fd-3m 有 τ≠0，但有反演 -> 仍可 IBZ

    def test_gan_requires_full_grid(self):
        g = self._gate("GaN")
        self.assertFalse(g["has_inversion"], g)
        self.assertTrue(g["needs_full_grid"], g)
        self.assertGreater(g["tau_ops"], 0, g)

    def test_mos2_aligned_ibz_ok(self):
        g = self._gate("MoS2")
        if g["max_tau"] > kc.TAU_TOL:
            self.skipTest("该 MoS2 POSCAR 尚未对齐原点（|tau|max=%.4f）" % g["max_tau"])
        self.assertFalse(g["has_inversion"], g)
        self.assertFalse(g["needs_full_grid"], g)

    def test_translation_flips_the_verdict(self):
        """GaAs（无反演、τ≡0）人为平移 0.13 -> τ≠0 -> 判据翻成"必须全网格"。

        这就是 MoS2 的处境：P-6m2 本身是简单空间群，只是所给原点不在高对称原子上；
        S3 的 align_origin 把原点对回去之后，同一个结构又回到 IBZ 可用态（V109/V111）。
        """
        p = Path(STRUCTS["GaAs"])
        if not p.is_file():
            self.skipTest("缺 %s" % p)
        st = _load("GaAs")
        if kc.symmetry_gate(st)["needs_full_grid"]:
            self.skipTest("起点不是 IBZ 可用态")
        st2 = st.copy()
        st2.translate_sites(list(range(len(st2))), [0.13, 0.07, 0.05], frac_coords=True)
        g2 = kc.symmetry_gate(st2)
        self.assertTrue(g2["needs_full_grid"], g2)
        self.assertFalse(g2["has_inversion"], g2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
