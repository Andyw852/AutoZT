#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_hybrid_params.py —— V165：杂化参数不是标准 HSE06/HSE03 要告警，并记进 band_summary（不需要 VASP）。

用法：python test_hybrid_params.py      退出码 0 = 全部 PASS。

起因：Si 的 S2.3 INCAR 是 HFSCREEN = 0.11（HSE06 的 0.2 Å⁻¹ 换成 bohr⁻¹ 才是 0.106），带隙 1.340 eV；
同一材料按 0.2 重跑是 1.091 eV。这个值进了 band_summary、S8 的剪刀差和手稿的表，一路没人报警。
"""
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


def _incar(text):
    p = Path(tempfile.mkdtemp()) / "INCAR"
    p.write_text(text)
    return p


class HybridTests(unittest.TestCase):
    def test_si_bohr_value_warned(self):
        h = kc.hybrid_params(_incar("LHFCALC = .TRUE.\nAEXX = 0.25\nHFSCREEN = 0.11   # HSE06\n"))
        self.assertFalse(h["standard"])
        self.assertIn("bohr", h["warning"])
        self.assertIn("0.2", h["warning"])

    def test_standard(self):
        for text, label in (("LHFCALC=.TRUE.; HFSCREEN=0.2\n", "HSE06/HSEsol"),
                            ("LHFCALC = T\nAEXX = 0.25\nHFSCREEN = 0.3\n", "HSE03")):
            h = kc.hybrid_params(_incar(text))
            self.assertTrue(h["standard"], text)
            self.assertEqual(h["label"], label)
            self.assertIsNone(h["warning"])

    def test_not_hybrid_or_pbe0(self):
        self.assertIsNone(kc.hybrid_params(_incar("ISMEAR = 0\n# LHFCALC = .TRUE.\n")))
        self.assertIsNone(kc.hybrid_params(Path(tempfile.mkdtemp()) / "missing"))
        h = kc.hybrid_params(_incar("LHFCALC = .TRUE.\nAEXX = 0.25\n"))
        self.assertIn("PBE0", h["label"])
        self.assertIsNone(h["warning"])

    def test_wiring(self):
        gen = (ROOT / "step2_bandgap/step2.3_hse/gen_step4_HSE.py").read_text(encoding="utf-8")
        i = gen.index('with open(os.path.join(out_dir, "INCAR"), "w") as f:')
        self.assertLess(i, gen.index("ke_common.hybrid_params(", i))
        plot = (ROOT / "step2_bandgap/step2.3_hse_plot/gen_step4.1_plot_band.py").read_text(encoding="utf-8")
        k = plot.index('_kc.hybrid_params(out / "INCAR")')
        self.assertLess(k, plot.index('(out / "band_summary.json").write_text('))
        self.assertIn('result["hybrid"] = _hy', plot[k:k + 400])


if __name__ == "__main__":
    unittest.main(verbosity=2)
