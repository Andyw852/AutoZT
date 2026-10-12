#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_pie3d.py —— [0037] 三维压电散射（PIE）的单元自检。

    cd skill/ke-dft-cpu/step8_amset && python test_pie3d.py

核（需要 amset；没装就跳过）：
1. ε 各向同性、不屏蔽：amset_pie3d 的核与原版 AMSET PiezoelectricScattering 逐点相同（锁住绝对归一化）；
2. 屏蔽：factor(κ) / factor(0) = q⁴/(q² + κ/ε)²，q→0 趋 0（原版趋常数）；
3. 张量整体旋转不变（e、ε、C、q 一起转）；
4. 纵场只看 q̂·ε·q̂：q̂ 沿 z 时改 ε_xx、ε_xz 不影响结果。
gen（纯 numpy，不需要 amset）：
5. read_piezo 取 including ionic 块并按标准 Voigt 重排；
6. resolve_pie_3d：没 PIE 不动、缺张量/张量为零/2D/缺静态介电时去掉 PIE、正常时给记录；
7. write_settings 写 piezoelectric_constant；_install_pie_3d 写记录 + 复制插件、原版核时删残留。
"""
import importlib.util
import io
import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT / "step8.4_amset2d"), str(ROOT), str(ROOT.parent / "_common" / "opt"),
           str(ROOT / "step2_bandgap" / "step2.15_discriminant")):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


try:
    import spglib as _sg                     # 新 spglib 去掉了 spglib.spglib，老 amset 还在 import 它
    if not hasattr(_sg, "spglib"):
        _sg.spglib = _sg
        sys.modules["spglib.spglib"] = _sg
    from amset.constants import coulomb_to_au, gpa_to_au, m_to_bohr
    from amset.scattering.elastic import PiezoelectricScattering as _STOCK
    from amset.util import cast_elastic_tensor, cast_piezoelectric_tensor
    _old = os.getcwd()
    _d = tempfile.mkdtemp()
    os.chdir(_d)
    Path(_d, "pie_3d.json").write_text('{"kernel": "az", "screening": true}')
    with redirect_stdout(io.StringIO()):
        P = _load("amset_pie3d", HERE / "amset_pie3d.py")
    os.chdir(_old)
    HAVE_AMSET = True
except ImportError:
    HAVE_AMSET = False

# 纤锌矿 GaN 量级（C/m²、GPa，标准 Voigt）
E_VOIGT = [[0, 0, 0, 0, -0.40, 0], [0, 0, 0, -0.40, 0, 0], [-0.49, -0.49, 0.73, 0, 0, 0]]
C_VOIGT = [[390, 145, 106, 0, 0, 0], [145, 390, 106, 0, 0, 0], [106, 106, 398, 0, 0, 0],
           [0, 0, 0, 105, 0, 0], [0, 0, 0, 0, 105, 0], [0, 0, 0, 0, 0, 122.5]]


def _unit(n, seed=0):
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    return u / np.linalg.norm(u, axis=1)[:, None], rng.uniform(1e-4, 0.5, n)


@unittest.skipUnless(HAVE_AMSET, "没装 amset")
class TestKernel(unittest.TestCase):
    def setUp(self):
        self.e = cast_piezoelectric_tensor(E_VOIGT) * coulomb_to_au / m_to_bohr ** 2
        self.c = cast_elastic_tensor(C_VOIGT) * gpa_to_au

    def test_matches_stock_isotropic_unscreened(self):
        eps = 9.5
        uq, qsq = _unit(60)
        st = _STOCK.__new__(_STOCK)
        st.elastic_constant, st._shape = self.c, np.ones((1, 1))
        st.inverse_screening_length_sq = np.full((1, 1), 1e-12)
        st.piezoelectric_constant = np.einsum("mn,mkl->nkl", np.linalg.inv(np.eye(3) * eps / (4 * np.pi)), self.e)
        ref = st.factor(uq, qsq, None, 0, None, None)
        got = P.pie_factor(self.e, np.eye(3) * eps, self.c, np.zeros((1, 1)), uq, qsq)
        np.testing.assert_allclose(got, ref, rtol=1e-8)

    def test_screening_form(self):
        eps, kap = 9.5, 0.3
        uq, qsq = _unit(40, 1)
        f0 = P.pie_factor(self.e, np.eye(3) * eps, self.c, np.zeros((1, 1)), uq, qsq)
        fk = P.pie_factor(self.e, np.eye(3) * eps, self.c, np.full((1, 1), kap), uq, qsq)
        np.testing.assert_allclose(fk / f0, (qsq ** 2 / (qsq + kap / eps) ** 2)[None, None], rtol=1e-12)
        tiny = P.pie_factor(self.e, np.eye(3) * eps, self.c, np.full((1, 1), kap), uq[:3], np.full(3, 1e-10))
        amp = P.pie_factor(self.e, np.eye(3) * eps, self.c, np.zeros((1, 1)), uq[:3], np.ones(3))   # 不屏蔽 q=1
        self.assertTrue(np.all(tiny / amp < 1e-6))     # q→0：q²/(κ/ε)² -> 0（原版趋常数 1/β²）

    def test_rotation_invariant(self):
        from scipy.spatial.transform import Rotation
        R = Rotation.from_euler("zyx", [0.3, -0.7, 1.1]).as_matrix()
        eps = np.diag([8.0, 9.0, 11.0])
        uq, qsq = _unit(30, 2)
        kap = np.array([[0.0, 0.2]])
        a = P.pie_factor(self.e, eps, self.c, kap, uq, qsq)
        e2 = np.einsum("ia,jb,kc,abc->ijk", R, R, R, self.e)
        c2 = np.einsum("ia,jb,kc,ld,abcd->ijkl", R, R, R, R, self.c)
        b = P.pie_factor(e2, R @ eps @ R.T, c2, kap, uq @ R.T, qsq)
        np.testing.assert_allclose(a, b, rtol=1e-9)

    def test_longitudinal_projection(self):
        uq, qsq = np.array([[0.0, 0.0, 1.0]]), np.array([0.01])
        f1 = P.pie_factor(self.e, np.diag([9.0, 9.0, 11.0]), self.c, np.zeros((1, 1)), uq, qsq)
        f2 = P.pie_factor(self.e, np.array([[20.0, 0, 3.0], [0, 20.0, 0], [3.0, 0, 11.0]]),
                          self.c, np.zeros((1, 1)), uq, qsq)
        np.testing.assert_allclose(f1, f2, rtol=1e-12)
        # 解析：q̂ = z 只有 LA（v = z）耦合 e33 -> (4π e33 / ε_zz)² / C33 / q²
        ref = (4 * np.pi * self.e[2, 2, 2] / 11.0) ** 2 / self.c[2, 2, 2, 2] / 0.01
        np.testing.assert_allclose(f1[0, 0, 0], ref, rtol=1e-10)


OUTCAR = """
 PIEZOELECTRIC TENSOR  for field in x, y, z        (C/m^2)
           XX          YY          ZZ          XY          YZ          ZX
  ---------------------------------------------------------------------------------
  x     0.00000     0.00000     0.00000     0.00000     0.00000    -0.20000
  y     0.00000     0.00000     0.00000     0.00000    -0.20000     0.00000
  z    -0.30000    -0.30000     0.50000     0.00000     0.00000     0.00000

 PIEZOELECTRIC TENSOR (including ionic contribution)  for field in x, y, z        (C/m^2)
           XX          YY          ZZ          XY          YZ          ZX
  ---------------------------------------------------------------------------------
  x     0.00000     0.00000     0.00000     0.00000     0.00000    -0.40000
  y     0.00000     0.00000     0.00000     0.00000    -0.40000     0.00000
  z    -0.49000    -0.49000     0.73000     0.00000     0.00000     0.00000
"""


class TestGen(unittest.TestCase):
    def setUp(self):
        self.g = _load("g10_pie3d", HERE / "gen_step10_amset.py")
        self.g.SCATTERING = ["ADP", "IMP", "POP", "PIE"]
        self.g.PIE_KERNEL, self.g.PIE_SCREENING = "az", True
        self.cwd = Path(tempfile.mkdtemp())
        (self.cwd / "step6_elastic").mkdir()

    def _outcar(self, text=OUTCAR):
        (self.cwd / "step6_elastic" / "OUTCAR").write_text(text)

    def _resolve(self, is_2d=False, elastic=C_VOIGT, eps=[[9, 0, 0], [0, 9, 0], [0, 0, 10]]):
        buf = io.StringIO()
        with redirect_stdout(buf):
            out = self.g.resolve_pie_3d(is_2d, self.cwd, elastic, eps)
        return out, buf.getvalue()

    def test_spec_keys(self):
        self.assertIn("PIE_KERNEL", self.g.SPEC)
        self.assertIn("PIE_SCREENING", self.g.SPEC)

    def test_read_piezo_ionic_voigt(self):
        self._outcar()
        with redirect_stdout(io.StringIO()):
            e = self.g.read_piezo(self.cwd)
        self.assertEqual(e, E_VOIGT)          # VASP (XX YY ZZ XY YZ ZX) -> Voigt (XX YY ZZ YZ XZ XY)

    def test_not_asked(self):
        self.g.SCATTERING = ["ADP", "IMP"]
        (piezo, rec), _ = self._resolve()
        self.assertIsNone(piezo)
        self.assertIsNone(rec)
        self.assertEqual(self.g.SCATTERING, ["ADP", "IMP"])

    def test_missing_tensor_drops(self):
        (piezo, rec), log = self._resolve()
        self.assertIsNone(rec)
        self.assertNotIn("PIE", self.g.SCATTERING)
        self.assertIn("LCALCEPS", log)

    def test_zero_tensor_drops(self):
        self._outcar(OUTCAR.replace("-0.40000", " 0.00000").replace("-0.49000", " 0.00000")
                     .replace("0.73000", "0.00000").replace("-0.20000", " 0.00000")
                     .replace("-0.30000", " 0.00000").replace("0.50000", "0.00000"))
        (_, rec), log = self._resolve()
        self.assertIsNone(rec)
        self.assertNotIn("PIE", self.g.SCATTERING)
        self.assertIn("中心对称", log)

    def test_2d_drops(self):
        self._outcar()
        (_, rec), log = self._resolve(is_2d=True)
        self.assertIsNone(rec)
        self.assertNotIn("PIE", self.g.SCATTERING)
        self.assertIn("S8.4", log)

    def test_no_static_dielectric_drops(self):
        self._outcar()
        (_, rec), _ = self._resolve(eps=None)
        self.assertIsNone(rec)
        self.assertNotIn("PIE", self.g.SCATTERING)

    def test_ok_and_settings(self):
        self._outcar()
        (piezo, rec), log = self._resolve()
        self.assertEqual(piezo, E_VOIGT)
        self.assertEqual(rec["kernel"], "az")
        self.assertTrue(rec["screening"])
        self.assertAlmostEqual(rec["e_max_C_m2"], 0.73)
        self.assertIn("PIE", self.g.SCATTERING)
        out = Path(tempfile.mkdtemp())
        buf = io.StringIO()
        with redirect_stdout(buf):
            ok = self.g._install_pie_3d(out, rec)
        self.assertTrue(ok)
        self.assertTrue((out / "amset_pie3d.py").is_file())
        self.assertEqual(json.loads((out / "pie_3d.json").read_text())["kernel"], "az")
        # 原版核：不装插件、删残留
        with redirect_stdout(io.StringIO()):
            self.assertFalse(self.g._install_pie_3d(out, dict(rec, kernel="amset")))
        self.assertFalse((out / "pie_3d.json").exists())

    def test_bad_kernel(self):
        self._outcar()
        self.g.PIE_KERNEL = "foo"
        with self.assertRaises(SystemExit):
            self._resolve()

    def test_settings_wiring(self):
        """write_settings 写 piezoelectric_constant；作业命令 import amset_pie3d。"""
        src = (HERE / "gen_step10_amset.py").read_text(encoding="utf-8")
        self.assertIn('lines.append("piezoelectric_constant:")', src)
        self.assertIn("piezo=_piezo", src)
        self.assertIn('"import amset_pie3d; " if _has_pie3d', src)

if __name__ == "__main__":
    unittest.main(verbosity=1)
