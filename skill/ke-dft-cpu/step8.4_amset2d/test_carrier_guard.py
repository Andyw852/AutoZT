#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_carrier_guard.py —— 载流子数防护的自检（离线，2D 与 3D 两个 gen 一起测）。

用法：python test_carrier_guard.py     退出码 0 = 全部 PASS。

覆盖（2026-09-28 用户指正后的口径：按**原子**归一）：
  · C60 网络原胞（约 146 A^2、60 原子）在出厂上限 1e14 cm^-2 下**不报**（按原胞会误拦）；
  · MoS2 那次单位事故（1.003469e19 cm^-2）仍然**报错**；
  · Si 用**原胞**（2 原子）与**惯用胞**（8 原子）算，结论一致；
  · Mg2C60 在 1e21 cm^-3 下不报；量级写错（1e23）报错；
  · 阈值：>0.5 个/原子报错、>0.05 个/原子告警。
"""
import importlib.util
import io
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent          # step8.4_amset2d
ROOT = HERE.parent                              # skill/ke-dft-cpu
for _p in (str(HERE), str(ROOT), str(ROOT / "step8_amset"),
           str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


G14 = _load("gen_step14_amset2d", HERE / "gen_step14_amset2d.py")
G10 = _load("gen_step10_amset", ROOT / "step8_amset" / "gen_step10_amset.py")

AREA_MOS2 = 8.65e-16          # MoS2 原胞面内面积 cm^2
AREA_C60 = 146e-16            # C60 网络原胞面内面积 cm^2（约 146 A^2）


def _call(fn, *a, **kw):
    """返回 (退出错误串 or None, 打印输出)。"""
    buf = io.StringIO()
    err = None
    with redirect_stdout(buf):
        try:
            fn(*a, **kw)
        except SystemExit as e:
            err = str(e)
    return err, buf.getvalue().strip()


class PerAtomGuardTests(unittest.TestCase):
    # ---------------- 2D ----------------
    def test_c60_network_not_blocked(self):
        """C60 网络：60 原子、146 A^2、出厂上限 1e14 cm^-2 -> 0.024 个/原子，不报不警。"""
        err, out = _call(G14.check_carriers_per_atom, "-1e14:-1e11:5, 1e11:1e14:5",
                         AREA_C60, 60, tag="C60net")
        self.assertIsNone(err, err)
        self.assertEqual(out, "", out)

    def test_mos2_accident_still_errors(self):
        err, _ = _call(G14.check_carriers_per_atom, "-1.003469e19:1.003469e19:2",
                       AREA_MOS2, 3, tag="MoS2")
        self.assertIsNotNone(err)
        self.assertIn("> 0.5", err)
        self.assertIn("每原子载流子数", err)

    def test_mos2_normal_doping_silent(self):
        err, out = _call(G14.check_carriers_per_atom, "-1e13:1e13:2", AREA_MOS2, 3)
        self.assertIsNone(err)
        self.assertEqual(out, "")

    def test_warn_band(self):
        """0.05 < 每原子 <= 0.5 -> 只告警。"""
        err, out = _call(G14.check_carriers_per_atom, "-2e14:2e14:2", AREA_MOS2, 3)
        self.assertIsNone(err)
        self.assertIn("[WARN]", out)

    def test_missing_inputs_silent(self):
        for args in ((AREA_MOS2, None), (None, 3)):
            err, out = _call(G14.check_carriers_per_atom, "-1e19:1e19:2", *args)
            self.assertIsNone(err)
            self.assertEqual(out, "")

    # ---------------- 3D ----------------
    def test_si_primitive_and_conventional_agree(self):
        """Si 原胞（2 原子、40.03 A^3）与惯用胞（8 原子、160.12 A^3）结论必须一致。"""
        v_prim = 40.03e-24      # cm^3
        v_conv = 160.12e-24
        for vol, nat in ((v_prim, 2), (v_conv, 8)):
            err, out = _call(G10.check_carriers_per_atom, "-1e22:1e22:2", vol, nat)
            self.assertIsNone(err, err)
            self.assertIn("[WARN]", out)
            # 每原子数必须相同（0.200）
            self.assertIn("0.200", out, out)

    def test_mg2c60_1e21_not_blocked(self):
        """Mg2C60 这类大原胞在 1e21 cm^-3 下每原子很小 -> 不报。"""
        err, out = _call(G10.check_carriers_per_atom, "-1e21:-1e17:5, 1e17:1e21:5",
                         2744e-24, 62)
        self.assertIsNone(err, err)
        self.assertEqual(out, "", out)

    def test_absurd_3d_errors(self):
        err, _ = _call(G10.check_carriers_per_atom, "-1e23:1e23:2", 2744e-24, 62)
        self.assertIsNotNone(err)
        self.assertIn("> 0.5", err)

    # ---------------- _count_atoms ----------------
    def test_count_atoms_vasp5_and_vasp4(self):
        import tempfile
        tmp = Path(tempfile.mkdtemp())
        d = tmp / "step3_uniform"
        d.mkdir(parents=True)
        hdr5 = "C60\n1.0\n10 0 0\n0 10 0\n0 0 20\nC\n60\nDirect\n"
        (d / "POSCAR").write_text(hdr5 + "0 0 0\n" * 60, encoding="utf-8")
        self.assertEqual(G14._count_atoms(tmp), 60)
        d2 = tmp / "step7_deform"
        d2.mkdir(parents=True)
        # 同时存在时 STRUCT_CANDS 先取 step3_uniform；换个目录测 VASP4
        tmp2 = Path(tempfile.mkdtemp())
        e = tmp2 / "step3_uniform"
        e.mkdir(parents=True)
        (e / "POSCAR").write_text("Si\n5.43\n5.43 0 0\n0 5.43 0\n0 0 5.43\n8\nDirect\n"
                                  + "0 0 0\n" * 8, encoding="utf-8")
        self.assertEqual(G10._count_atoms(tmp2), 8)


if __name__ == "__main__":
    unittest.main(verbosity=2)
