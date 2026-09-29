#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_compare_deformation.py —— tools/compare_deformation_h5.py 的分诊判定（纯本地、无数据）。

用法：python test_compare_deformation.py      退出码 0 = 全部 PASS。
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "tools"))
import compare_deformation_h5 as C  # noqa: E402


def _h5(path, D, kp, sym=0):
    import h5py
    with h5py.File(str(path), "w") as f:
        f.create_dataset("deformation_potentials_up", data=D)
        f["kpoints"] = kp
        f["structure"] = np.bytes_("{}")
        f.attrs["dp_symmetrized"] = sym


class CompareTests(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp())
        rng = np.random.default_rng(0)
        self.kp = rng.random((40, 3))
        D = rng.normal(size=(3, 40, 3, 3)) * 5
        self.D = 0.5 * (D + np.swapaxes(D, -1, -2))
        _h5(self.d / "old.h5", self.D, self.kp)

    def test_identical_is_small(self):
        _h5(self.d / "new.h5", self.D, self.kp, sym=1)
        self.assertEqual(C.main([str(self.d / "old.h5"), str(self.d / "new.h5")]), 0)

    def test_mos2_like_doubling_needs_rerun(self):
        """MoS2：电子 ADP 迁移率 ×2 ≈ ⟨|D|²⟩ 减半 ≈ D × 0.707。"""
        _h5(self.d / "new.h5", self.D * 0.707, self.kp, sym=1)
        rep = C.compare(self.d / "old.h5", self.d / "new.h5")
        s = rep["spins"]["deformation_potentials_up"]
        self.assertAlmostEqual(s["R_mean_D2"], 0.5, places=2)
        self.assertAlmostEqual(s["mu_adp_factor_rough"], 2.0, places=1)
        self.assertTrue(rep["rerun"])
        self.assertEqual(C.main([str(self.d / "old.h5"), str(self.d / "new.h5")]), 1)

    def test_small_symmetrization_noise_is_small(self):
        rng = np.random.default_rng(1)
        noise = 1 + 0.03 * rng.normal(size=self.D.shape)
        _h5(self.d / "new.h5", self.D * noise, self.kp, sym=1)
        rep = C.compare(self.d / "old.h5", self.d / "new.h5")
        self.assertFalse(rep["rerun"], rep)

    def test_different_grid_is_input_error(self):
        _h5(self.d / "new.h5", self.D[:, :30], self.kp[:30], sym=1)
        self.assertEqual(C.main([str(self.d / "old.h5"), str(self.d / "new.h5")]), 2)

    def test_band_edge_hits_reported(self):
        _h5(self.d / "new.h5", self.D * 0.9, self.kp, sym=1)
        be = self.d / "band_edges.json"
        be.write_text(json.dumps({"electron": {"hits": [{"spin": "up", "band_h5": 1, "k_h5": 5}]}}))
        rep = C.compare(self.d / "old.h5", self.d / "new.h5", be)
        self.assertEqual(len(rep["band_edge_hits"]), 1)
        self.assertEqual(rep["new_attrs"]["dp_symmetrized"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
