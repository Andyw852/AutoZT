# -*- coding: utf-8 -*-
"""离线单测：nspin_norm_fix 的纯逻辑 + h5 集成（不开 AMSET、不碰项目数据）。

跑法（worktree 根目录）：
    PYTHONPATH=skill/ke-dft-cpu \
      ~/miniconda3/envs/atomate2_p_a/bin/python \
      -m unittest tests.test_nspin_norm_fix -v

h5 集成部分只造临时合成 h5（deformation_potentials_up/down 全 1.0），
不读任何真实项目数据；无 h5py 时整类 SkipTest。
"""
import os
import sys
import tempfile
import unittest

# 双保险：即便没设 PYTHONPATH，也能从 worktree 的 skill/ke-dft-cpu 导入
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "skill", "ke-dft-cpu"))

import nspin_norm_fix as nnf  # noqa: E402

try:
    import h5py
    import numpy as np
except ImportError:                                       # pragma: no cover
    h5py = None
    np = None


class TestPureLogic(unittest.TestCase):
    """decide() 的四种情形 + 版本解析工具。"""

    def test_amset_0419_nspin2_should_fix(self):
        self.assertEqual(nnf.decide("0.4.19", 2, False),
                         (True, "amset<0.5.1 & nspin=2"))

    def test_amset_0419_nspin1_no_fix(self):
        self.assertEqual(nnf.decide("0.4.19", 1, False), (False, "nspin=1"))

    def test_amset_051_nspin2_no_fix(self):
        self.assertEqual(nnf.decide("0.5.1", 2, False), (False, "amset>=0.5.1"))

    def test_marked_no_fix(self):
        self.assertEqual(nnf.decide("0.4.19", 2, True), (False, "already-marked"))

    def test_version_helpers(self):
        self.assertTrue(nnf.version_lt("0.4.19", "0.5.1"))
        self.assertTrue(nnf.version_lt("0.5", "0.5.1"))
        self.assertTrue(nnf.version_lt("0.5.0", "0.5.1"))
        self.assertFalse(nnf.version_lt("0.5.1", "0.5.1"))
        self.assertFalse(nnf.version_lt("0.5.2", "0.5.1"))
        self.assertEqual(nnf.parse_version("0.4.19"), (0, 4, 19))
        self.assertEqual(nnf.parse_version("0.5.1.dev0+g123"), (0, 5, 1))
        self.assertEqual(nnf.parse_version("0.5.1rc1"), (0, 5, 1))


class TestGate(unittest.TestCase):
    """硬闸门判据（CrSe2_hex 实测数）。"""

    def test_gate_rejects_halved_value(self):
        # h5 vac=2.4465 是减半值，band_edges E1_vac=4.8808 是正确值 -> 差 ~49.9%
        self.assertFalse(nnf.gate_ok(2.4465, 4.8808))

    def test_gate_accepts_corrected_value(self):
        # 乘回 2 后 4.893 对 4.8808，差 ~0.25% < 5%
        self.assertTrue(nnf.gate_ok(4.89, 4.8808))

    def test_relative_diff(self):
        self.assertAlmostEqual(nnf.relative_diff(4.89, 4.8808), 0.001885, places=5)
        self.assertEqual(nnf.relative_diff(0.0, 0.0), 0.0)
        self.assertEqual(nnf.relative_diff(1.0, 0.0), float("inf"))


def _make_h5(path, spins, value=1.0):
    """造一份合成形变势 h5：每个自旋道 dataset 全 value，另加 kpoints 占位。"""
    with h5py.File(path, "w") as fh:
        for s in spins:
            fh.create_dataset("deformation_potentials_" + s,
                              data=np.full((2, 3, 3, 3), value, dtype="f8"))
        fh.create_dataset("kpoints", data=np.zeros((3, 3), dtype="f8"))


@unittest.skipUnless(h5py is not None, "h5py 不可用，跳过 h5 集成测试")
class TestFixH5(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="nspin_fix_")
        self.addCleanup(self.tmp.cleanup)

    def _path(self, name="deform.h5"):
        return os.path.join(self.tmp.name, name)

    def test_double_spin_fix_then_idempotent(self):
        p = self._path()
        _make_h5(p, ["up", "down"])
        fixed, reason = nnf.fix_h5(p, "0.4.19")
        self.assertTrue(fixed)
        self.assertEqual(reason, "amset<0.5.1 & nspin=2")
        with h5py.File(p, "r") as fh:
            self.assertTrue(np.allclose(fh["deformation_potentials_up"][...], 2.0))
            self.assertTrue(np.allclose(fh["deformation_potentials_down"][...], 2.0))
            self.assertTrue(bool(fh.attrs["nspin_norm_fixed"]))
            self.assertEqual(int(fh.attrs["nspin_norm_fix_factor"]), 2)
            self.assertEqual(str(fh.attrs["nspin_norm_fix_version"]), "0.4.19")
            self.assertEqual(str(fh.attrs["nspin_norm_fix_reason"]),
                             "amset<0.5.1 & nspin=2")
        # 第二次调用：already-marked，数值不再翻倍（幂等）
        fixed2, reason2 = nnf.fix_h5(p, "0.4.19")
        self.assertFalse(fixed2)
        self.assertEqual(reason2, "already-marked")
        with h5py.File(p, "r") as fh:
            self.assertTrue(np.allclose(fh["deformation_potentials_up"][...], 2.0))
            self.assertTrue(np.allclose(fh["deformation_potentials_down"][...], 2.0))

    def test_single_spin_unchanged(self):
        p = self._path("single.h5")
        _make_h5(p, ["up"])
        fixed, reason = nnf.fix_h5(p, "0.4.19")
        self.assertFalse(fixed)
        self.assertEqual(reason, "nspin=1")
        with h5py.File(p, "r") as fh:
            self.assertTrue(np.allclose(fh["deformation_potentials_up"][...], 1.0))
            self.assertEqual(str(fh.attrs["nspin_norm_fixed"]), "skipped")

    def test_new_amset_unchanged(self):
        p = self._path("new.h5")
        _make_h5(p, ["up", "down"])
        fixed, reason = nnf.fix_h5(p, "0.5.1")
        self.assertFalse(fixed)
        self.assertEqual(reason, "amset>=0.5.1")
        with h5py.File(p, "r") as fh:
            self.assertTrue(np.allclose(fh["deformation_potentials_up"][...], 1.0))
            self.assertTrue(np.allclose(fh["deformation_potentials_down"][...], 1.0))
            self.assertEqual(str(fh.attrs["nspin_norm_fixed"]), "skipped")


if __name__ == "__main__":
    unittest.main(verbosity=2)
