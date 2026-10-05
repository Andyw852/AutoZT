#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_core_ref_kind.py —— V169：S7.1 逐构型判定芯能级参考来源（平均静电芯势 / 1s 芯能级），混用报错（纯本地）。

起因（CrS₂ S7，2026-10-06）：V168 的回显核对报出旧 INCAR 带 ICORELEVEL=1。这时 VASP 不写平均静电芯势块，
AMSET 的 get_reference_energy 会静默改用 1s 芯能级——部分构型重算之后就是两种量混着减。

用法：python test_core_ref_kind.py      退出码 0 = 全部 PASS。
"""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
import ke_common as kc  # noqa: E402

AVG = """ average (electrostatic) potential at core
  the test charge radii are     0.8726  0.7552
  (the norm of the test charge is              1.0000)
      1 -40.1234      2 -71.5678
 E-fermi :  -1.2345
"""
CORE = """ the core state eigenenergies are
    1-    1s  -5988.1234   2s   -697.3456
 E-fermi :  -1.2345
"""


class CoreRefKindTests(unittest.TestCase):
    def test_kinds(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / "a").write_text("head\n" + AVG + "General timing and accounting informations\n")
            (d / "b").write_text("head\n" + CORE)
            (d / "c").write_text("head\nGeneral timing and accounting informations\n")
            with gzip.open(d / "e.gz", "wt") as f:
                f.write(AVG)
            self.assertEqual(kc.outcar_core_ref_kind(d / "a"), "avg_core")
            self.assertEqual(kc.outcar_core_ref_kind(d / "b"), "1s")
            self.assertIsNone(kc.outcar_core_ref_kind(d / "c"))
            self.assertIsNone(kc.outcar_core_ref_kind(d / "missing"))
            self.assertEqual(kc.outcar_core_ref_kind(d / "e"), "avg_core")       # 只剩 OUTCAR.gz 也认


class S71WiringTests(unittest.TestCase):
    SRC = (_ROOT / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py").read_text(encoding="utf-8")

    def test_checked_before_amset_read_and_recorded(self):
        s = self.SRC
        i = s.index("kc.outcar_core_ref_kind(")
        self.assertLess(s.index("dp_subs = [_dp_folder(p)"), i)                  # 判的是 amset 实际读的目录
        self.assertLess(i, s.index("amset deform read undeformed"))
        self.assertLess(s.index("芯能级参考来源不一致"), s.index("amset deform read undeformed"))
        self.assertIn('_be["core_reference"] = _core_kind', s)


if __name__ == "__main__":
    unittest.main(verbosity=2)
