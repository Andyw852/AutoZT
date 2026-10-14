# -*- coding: utf-8 -*-
"""[V187] 形变势构型口径：3D 也做离子弛豫、弛豫力判据收紧 + 核对、z 方向按维度处理（纯本地、无数据）。

  · ke_common.ionrelax_force_check：首/末步最大力、跑满 NSW、回收不足、判据偏松；
  · gen_step9_deform.select_ionrelax_dirs：3D 全部分量、2D 只要面内（面外拉的是真空）、IONRELAX_3D=false 不建；
  · gen_step9_deform._ionrelax_params / _ionrelax_incars：step.conf 覆盖力判据，弛豫段写进 INCAR.relax；
  · 下游：S8.4 / S8.3 不再把 2D 的 "ionrelax(in-plane)+clamped(out-of-plane)" 当成非离子弛豫告警；
    S7.1 的 3D E1_iso 取三个对角。
"""
import ast
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SKILL = REPO / "skill" / "ke-dft-cpu"
for _p in (str(REPO), str(SKILL), str(SKILL / "step7_deform")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import numpy as np  # noqa: F401
    HAVE_NP = True
except ImportError:                                            # pragma: no cover
    HAVE_NP = False


def _kc():
    import ke_common
    return ke_common


def _gen9():
    # 模块级要解析 AMSET 环境（只拼字符串，不执行）；测试里给个假的，免得 SystemExit
    os.environ.setdefault("AUTOZT_CONDA_SH", "/dev/null")
    os.environ.setdefault("AUTOZT_AMSET_ENV", "amset-test")
    import gen_step9_deform
    return gen_step9_deform


def _force_block(forces):
    rows = "\n".join("     0.00000      0.00000      0.00000      %.6f      %.6f      %.6f"
                     % (fx, fy, fz) for fx, fy, fz in forces)
    return (" POSITION                                       TOTAL-FORCE (eV/Angst)\n"
            " -----------------------------------------------------------------------------------\n"
            + rows + "\n"
            " -----------------------------------------------------------------------------------\n"
            "    total drift:                                0.000000      0.000000      0.000000\n")


def _write_relax(d, steps, ediffg="-0.005", nsw=60):
    d.mkdir(parents=True, exist_ok=True)
    (d / "INCAR.relax").write_text("IBRION = 2\nNSW    = %d\nEDIFFG = %s\n" % (nsw, ediffg))
    (d / "OUTCAR.relax").write_text("".join(_force_block(f) for f in steps))


class ForceCheckTests(unittest.TestCase):
    def test_max_force_per_ionic_step(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t) / "ionrelax"
            _write_relax(d, [[(0.03, 0.04, 0.0), (0.0, 0.0, 0.01)], [(0.0, 0.0, 0.002)]])
            fs = _kc().outcar_max_forces(d / "OUTCAR.relax")
            self.assertEqual(len(fs), 2)
            self.assertAlmostEqual(fs[0], 0.05, places=6)
            self.assertAlmostEqual(fs[1], 0.002, places=6)

    def test_ok(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t) / "ionrelax"
            _write_relax(d, [[(0.08, 0, 0)], [(0.02, 0, 0)], [(0.004, 0, 0)]])
            r = _kc().ionrelax_force_check(d)
            self.assertEqual(r["status"], "ok", r)
            self.assertEqual(r["n_ionic"], 3)
            self.assertAlmostEqual(r["ratio_last_first"], 0.05, places=4)

    def test_weak_recovery_with_old_loose_criterion(self):
        # MoS2 量级：应变诱导的力 ~0.04 eV/Å，旧判据 -0.02 收敛时可能只回收了一半
        with tempfile.TemporaryDirectory() as t:
            d = Path(t) / "ionrelax"
            _write_relax(d, [[(0.04, 0, 0)], [(0.018, 0, 0)]], ediffg="-0.02")
            r = _kc().ionrelax_force_check(d)
            self.assertEqual(r["status"], "weak", r)
            self.assertIn("loose", r["flags"])

    def test_nsw_hit(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t) / "ionrelax"
            _write_relax(d, [[(0.1, 0, 0)], [(0.05, 0, 0)], [(0.03, 0, 0)]], nsw=3)
            r = _kc().ionrelax_force_check(d)
            self.assertEqual(r["status"], "nsw_hit", r)

    def test_small_initial_force_is_ok(self):
        # 高对称 3D（岩盐等）：受力本来就小于判据，弛豫 1 步就停——不是"回收不足"
        with tempfile.TemporaryDirectory() as t:
            d = Path(t) / "ionrelax"
            _write_relax(d, [[(0.0005, 0, 0)]])
            self.assertEqual(_kc().ionrelax_force_check(d)["status"], "ok")

    def test_no_relax_outcar(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t) / "ionrelax"
            d.mkdir()
            (d / "INCAR.relax").write_text("EDIFFG = -0.005\n")
            self.assertEqual(_kc().ionrelax_force_check(d)["status"], "no_data")


_CELL = ((5.0, 0.0, 0.0), (0.0, 5.0, 0.0), (0.0, 0.0, 20.0))


def _poscar(path, eps):
    """eps：3×3 对称小应变；晶格行矢量 a_i' = a_i (I + eps)。"""
    import numpy as np
    L = np.array(_CELL) @ (np.eye(3) + np.array(eps))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x\n1.0\n" + "".join("  %.10f %.10f %.10f\n" % tuple(r) for r in L)
                    + "Si\n1\nDirect\n0 0 0.5\n")


def _strain(i, j, v=0.005):
    e = [[0.0] * 3 for _ in range(3)]
    e[i][j] += v
    if i != j:
        e[j][i] += v
    return e


@unittest.skipUnless(HAVE_NP, "需要 numpy")
class SelectDirsTests(unittest.TestCase):
    def _tree(self, t):
        out = Path(t) / "step7_deform"
        _poscar(out / "undeformed" / "POSCAR", [[0] * 3] * 3)
        comps = {"deform-01": (0, 0), "deform-02": (1, 1), "deform-03": (2, 2),
                 "deform-04": (1, 2), "deform-05": (0, 2), "deform-06": (0, 1)}
        for name, (i, j) in comps.items():
            _poscar(out / name / "POSCAR", _strain(i, j))
        return out

    def test_3d_all_components(self):
        g = _gen9()
        with tempfile.TemporaryDirectory() as t:
            out = self._tree(t)
            self.assertEqual(g.select_ionrelax_dirs(out, "3d", True),
                             {"deform-%02d" % i for i in range(1, 7)})

    def test_3d_opt_out_keeps_clamped(self):
        g = _gen9()
        with tempfile.TemporaryDirectory() as t:
            self.assertEqual(g.select_ionrelax_dirs(self._tree(t), "3d", False), set())

    def test_2d_in_plane_only(self):
        g = _gen9()
        with tempfile.TemporaryDirectory() as t:
            self.assertEqual(g.select_ionrelax_dirs(self._tree(t), "2d", True),
                             {"deform-01", "deform-02", "deform-06"})


class IonrelaxParamsTests(unittest.TestCase):
    def test_defaults(self):
        g = _gen9()
        self.assertEqual(g.IONRELAX_EDIFFG, "-0.005")
        self.assertEqual(float(g.IONRELAX_EDIFFG), -_kc().IONRELAX_EDIFFG_REF)
        self.assertTrue(g.IONRELAX_3D)
        with tempfile.TemporaryDirectory() as t:
            self.assertEqual(g._ionrelax_params(t), ("-0.005", True, False))

    def test_step_conf_override(self):
        g = _gen9()
        with tempfile.TemporaryDirectory() as t:
            Path(t, "step.conf").write_text("[params]\nIONRELAX_EDIFFG = -0.003\nIONRELAX_3D = false\n"
                                            "IONRELAX_REBUILD_LOOSE = true\nFUNC = pbesol\n")
            self.assertEqual(g._ionrelax_params(t), ("-0.003", False, True))

    def test_positive_ediffg_rejected(self):
        g = _gen9()
        with tempfile.TemporaryDirectory() as t:
            Path(t, "step.conf").write_text("[params]\nIONRELAX_EDIFFG = 0.003\n")
            with self.assertRaises(SystemExit):
                g._ionrelax_params(t)

    def test_relax_segment_uses_ediffg(self):
        g = _gen9()
        base = ("ISTART = 0\nICHARG = 2\nIBRION = -1\nNSW    = 0\nISIF   = 2\nISYM   = 0\n"
                "EDIFF  = 1E-8\nLWAVE  = .FALSE.\nLCHARG = .FALSE.\nLVHAR  = .TRUE.\n")
        relax, stat = g._ionrelax_incars(base)
        self.assertRegex(relax, r"EDIFFG = -0\.005\b")
        self.assertRegex(relax, r"IBRION = 2\b")
        self.assertRegex(relax, r"ISIF   = 2")               # 固定晶胞：应变不被弛豫掉
        self.assertNotIn("EDIFFG", stat)
        relax2, stat2 = g._ionrelax_incars(base, "-0.002")
        self.assertRegex(relax2, r"EDIFFG = -0\.002\b")
        self.assertEqual(stat, stat2)                         # 静态段与力判据无关（V168 复用判据不受影响）


def _func_from_source(path, name):
    src = Path(path).read_text(encoding="utf-8")
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            ns = {}
            exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), ns)
            return ns[name]
    raise AssertionError("%s 里没有 %s" % (path, name))


class DownstreamTests(unittest.TestCase):
    def test_s84_accepts_2d_inplane_relaxed(self):
        ok = _func_from_source(SKILL / "step8.4_amset2d" / "gen_step14_amset2d.py", "_deform_geom_ok_2d")
        self.assertTrue(ok("ionrelax"))
        self.assertTrue(ok("ionrelax(in-plane)+clamped(out-of-plane)"))
        self.assertFalse(ok("clamped"))
        self.assertFalse(ok("mixed"))
        self.assertFalse(ok(None))

    def test_s71_3d_iso_uses_three_diagonals(self):
        src = (SKILL / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py").read_text(encoding="utf-8")
        self.assertIn("(t[0, 0] + t[1, 1] + t[2, 2]) / 3 if _DIM == \"3d\"", src)
        self.assertIn('"E1_zz_eV"', src)
        self.assertIn("shlex.quote(_dim_now or \"\")", src)      # 维度传进嵌入脚本
        self.assertIn('_be["ionrelax_check"] = _ir_check', src)

    def test_band_edges_rev_bumped(self):
        self.assertNotEqual(_kc().BAND_EDGES_REV, "2026-09-16-deform-ref-vacuum")


if __name__ == "__main__":
    unittest.main()
