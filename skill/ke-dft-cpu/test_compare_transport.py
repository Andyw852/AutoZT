#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_compare_transport.py —— tools/compare_transport_json.py 的判定（纯本地、无数据）。

用法：python test_compare_transport.py      退出码 0 = 全部 PASS。
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "tools"))
import compare_transport_json as C  # noqa: E402


def _write(path, mob_scale=1.0, xx_yy=(1.0, 1.0), doping=(-1e19, 1e19), extra_mech=False):
    nd, nt = len(doping), 2
    t = np.zeros((nd, nt, 3, 3))
    t[..., 0, 0], t[..., 1, 1], t[..., 2, 2] = xx_yy[0], xx_yy[1], 1e-3
    base = np.arange(1, nd * nt + 1, dtype=float).reshape(nd, nt)[..., None, None]
    mob = {"overall": (base * t * 100 * mob_scale).tolist(), "ADP": (base * t * 300 * mob_scale).tolist()}
    if extra_mech:
        mob["PIE"] = (base * t * 900).tolist()
    Path(path).write_text(json.dumps({"doping": list(doping), "temperatures": [300, 600],
                                      "conductivity": (base * t * 1e4 * mob_scale).tolist(),
                                      "seebeck": (base * t * -200).tolist(), "mobility": mob}))


class CompareTransportTests(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp())
        _write(self.d / "old.json", xx_yy=(1.02, 0.98))               # 旧：面内 ±2% 数值噪声

    def test_symmetrized_new_passes_on_average(self):
        _write(self.d / "new.json")                                    # 新：xx=yy，平均不变
        rep = C.compare(self.d / "old.json", self.d / "new.json")
        self.assertTrue(rep["pass"], rep)
        q = rep["quantities"]["mobility:overall"]
        self.assertLess(q["max_rel"], 1e-9)
        self.assertAlmostEqual(q["aniso_old"], 0.04, places=6)
        self.assertEqual(q["aniso_new"], 0.0)
        self.assertEqual(C.main([str(self.d / "old.json"), str(self.d / "new.json")]), 0)

    def test_real_change_fails(self):
        _write(self.d / "new.json", mob_scale=1.05)
        rep = C.compare(self.d / "old.json", self.d / "new.json")
        self.assertFalse(rep["pass"])
        self.assertFalse(rep["quantities"]["mobility:ADP"]["pass"])
        self.assertTrue(rep["quantities"]["seebeck"]["pass"])
        self.assertEqual(C.main([str(self.d / "old.json"), str(self.d / "new.json")]), 1)

    def test_grid_mismatch_is_input_error(self):
        _write(self.d / "new.json", doping=(-1e19, 1e20))
        self.assertEqual(C.main([str(self.d / "old.json"), str(self.d / "new.json")]), 2)

    def test_mechanism_sets_reported(self):
        _write(self.d / "new.json", extra_mech=True)
        rep = C.compare(self.d / "old.json", self.d / "new.json")
        self.assertEqual(rep["only_new"], ["mobility:PIE"])
        self.assertTrue(rep["pass"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
