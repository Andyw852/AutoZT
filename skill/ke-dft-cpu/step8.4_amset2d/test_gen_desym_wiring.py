#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_gen_desym_wiring.py —— S8 / S8.4 两个 gen 的 V115 接线自检（不跑 VASP/AMSET，结构在测试里造）。

用法：python test_gen_desym_wiring.py      退出码 0 = 全部 PASS。

固定四件事：
  1) 全网格判定用**逐操作判据**（旧规则"有反演或 τ≡0 就走 IBZ"两个方向都错）：
       Si 原点在反演中心 -> 必须全网格（旧规则放行 = 静默算错）
       GaN 标准原点      -> IBZ 即可（旧规则要全网格 = 白算）
     DESYM_FIX 开 -> 都走 IBZ，但判据摘要里照打 bad_ops/total_ops；
  2) AMSET 命令里有且只有一个 @AMSET_PLUGINS@ 占位（3D 也是 python -c 入口，日志写法同 2D）；
  3) DESYM_FIX 开着却找不到插件 -> 拒绝生成；
  4) 2D ε∞ 覆盖的 layer/slab 口径换算与"离子部分不变"；带隙覆盖通道与来源记录。
"""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT), str(ROOT / "step8_amset"), str(ROOT.parent / "_common" / "opt"),
           str(ROOT / "step2_bandgap" / "step2.15_discriminant")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pymatgen.io.vasp.inputs import Poscar  # noqa: E402
import test_symmetry_gate as tsg            # noqa: E402  (结构构造器)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


G10 = _load("g10_under_test", ROOT / "step8_amset" / "gen_step10_amset.py")
G14 = _load("g14_under_test", HERE / "gen_step14_amset2d.py")


def _matdir(st):
    d = Path(tempfile.mkdtemp())
    (d / "step3_uniform").mkdir()
    Poscar(st).write_file(str(d / "step3_uniform" / "POSCAR"))
    return d


def _decide(mod, st, fix):
    mod.WAVEFUNCTION_FULL = None
    mod._WF_EXPLICIT = None
    mod._DESYM_FIX_ON = fix
    if hasattr(mod, "_GATE_CACHE"):
        mod._GATE_CACHE.clear()
    if mod is G10:
        mod.UNITY_OVERLAP = "false"          # 真实重叠
    else:
        mod.UNITY_OVERLAP = False
    mod._apply_inversion_rule(_matdir(st))
    return mod.WAVEFUNCTION_FULL, mod._GATE_SUMMARY


class PerOpRuleTests(unittest.TestCase):
    CASES = [("Si 原点在反演中心", tsg.shifted(tsg.si_on_atom(), [-.125] * 3), True, 24),
             ("GaN 标准原点", tsg.gan_std(), False, 0),
             ("MoS2 对齐前", tsg.shifted(tsg.mos2_aligned(), [.1, .2, 0]), True, 12)]

    def test_rule_fix_off(self):
        for mod in (G10, G14):
            for name, st, need, bad in self.CASES:
                full, summ = _decide(mod, st, False)
                self.assertEqual(full, need, "%s %s: %s" % (mod.__name__, name, summ))
                self.assertIn("bad_ops=%d/" % bad, summ)

    def test_rule_fix_on_releases_ibz_but_logs(self):
        for mod in (G10, G14):
            for name, st, _need, bad in self.CASES:
                full, summ = _decide(mod, st, True)
                self.assertFalse(full, "%s %s: %s" % (mod.__name__, name, summ))
                self.assertIn("bad_ops=%d/" % bad, summ)
                self.assertIn("DESYM_FIX=on", summ)


class CommandTests(unittest.TestCase):
    def test_plugin_placeholder_and_logging(self):
        for mod in (G10, G14):
            self.assertEqual(mod.AMSET_CMD.count("@AMSET_PLUGINS@"), 1, mod.__name__)
            self.assertIn("initialize_amset_logger as L; L();", mod.AMSET_CMD)
            self.assertIn("Runner.from_directory('.').run()", mod.AMSET_CMD)
            self.assertIn(">> amset.log 2>&1", mod.AMSET_CMD)
            self.assertNotIn("amset run", mod.AMSET_CMD)
        # 3D：费米窗口检查仍在 cp 成 transport.json 之前（V63 的顺序）
        c = G10.AMSET_CMD
        self.assertLess(c.index("fermi_window_check"), c.index("cp -f"))

    def test_install_refuses_when_on_and_missing(self):
        for mod in (G10, G14):
            d = Path(tempfile.mkdtemp())
            mod._DESYM_FIX_ON = True
            real = mod.__file__
            try:
                mod.__file__ = str(d / "gen.py")      # 让 here 指向空目录
                import os
                old = os.getcwd()
                os.chdir(d)
                try:
                    with self.assertRaises(SystemExit):
                        mod._install_desym_fix(d)
                    mod._DESYM_FIX_ON = False
                    self.assertFalse(mod._install_desym_fix(d))
                finally:
                    os.chdir(old)
            finally:
                mod.__file__ = real
                mod._DESYM_FIX_ON = False

    def test_install_copies_plugin(self):
        # 3D S8 的插件住在 step8.4_amset2d/，gen 时靠 find_asset 复制进运行目录（cwd），
        # 不在 gen 自己的目录里 —— 这里模拟那份复制，验证 _install_desym_fix 能落盘。
        import os
        import shutil
        plugin = HERE / "amset_desym_fix.py"
        for mod in (G10, G14):
            d = Path(tempfile.mkdtemp())
            shutil.copyfile(plugin, d / "amset_desym_fix.py")
            old = os.getcwd()
            os.chdir(d)
            try:
                self.assertTrue(mod._install_desym_fix(d))
            finally:
                os.chdir(old)
            self.assertTrue((d / "amset_desym_fix.py").is_file())


class OverrideTests(unittest.TestCase):
    def test_eps_override_2d_layer_and_slab(self):
        inf = [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 1.5]]
        stat = [[4.5, 0, 0], [0, 4.5, 0], [0, 0, 1.6]]
        f = 3.0                                          # c/t
        try:
            G14.EPS_INF_OVERRIDE, G14.EPS_INF_OVERRIDE_BASIS = "7.0", "layer"
            ni, ns, rec = G14._apply_eps_inf_override_2d(inf, stat, f)
            self.assertAlmostEqual(ni[0][0], 1 + 6.0 / 3)      # layer 7 -> slab 3
            self.assertAlmostEqual(ni[2][2], 1.5)              # 标量不动 zz
            self.assertAlmostEqual(ns[0][0] - ni[0][0], 0.5)   # 离子部分不变
            self.assertAlmostEqual(G14._EPS_OV_LAYER[0][0][0], 7.0)
            self.assertEqual(rec["basis"], "layer")
            G14.EPS_INF_OVERRIDE_BASIS = "slab"
            ni, _ns, _rec = G14._apply_eps_inf_override_2d(inf, stat, f)
            self.assertAlmostEqual(ni[0][0], 7.0)
            self.assertAlmostEqual(G14._EPS_OV_LAYER[0][0][0], 1 + 3 * 6.0)
            G14.EPS_INF_OVERRIDE_BASIS = ""
            with self.assertRaises(SystemExit):              # 2D 不给口径 -> 停
                G14._apply_eps_inf_override_2d(inf, stat, f)
        finally:
            G14.EPS_INF_OVERRIDE, G14.EPS_INF_OVERRIDE_BASIS = None, ""

    def test_eps_override_3d_path_of_gen10(self):
        inf = [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]]
        stat = [[12.0, 0, 0], [0, 12.0, 0], [0, 0, 12.0]]
        try:
            G10.EPS_INF_OVERRIDE = "8.0"
            ni, ns = G10._apply_eps_inf_override(inf, stat)
            self.assertAlmostEqual(ni[0][0], 8.0)
            self.assertAlmostEqual(ns[0][0], 10.0)
            with self.assertRaises(SystemExit):              # 2D 不给口径 -> 停
                G10._apply_eps_inf_override(inf, stat, is_2d=True, out=Path("."))
        finally:
            G10.EPS_INF_OVERRIDE, G10.EPS_INF_OVERRIDE_BASIS = None, ""

    def test_bandgap_override(self):
        d = Path(tempfile.mkdtemp())
        for mod in (G10, G14):
            try:
                mod.BANDGAP_OVERRIDE, mod.BANDGAP_NOTE = 1.23, "PBE0"
                self.assertAlmostEqual(mod.read_bandgap(d), 1.23)
                self.assertIn("BANDGAP_OVERRIDE", mod._BANDGAP_SRC)
                self.assertIn("PBE0", mod._BANDGAP_SRC)
                mod.BANDGAP_OVERRIDE = 0.0
                with self.assertRaises(SystemExit):
                    mod.read_bandgap(d)
            finally:
                mod.BANDGAP_OVERRIDE, mod.BANDGAP_NOTE, mod._BANDGAP_SRC = None, "", None

    def test_2d_mesh_target_default(self):
        self.assertEqual(G14.MESH_MIN, 263)
        self.assertGreaterEqual(G14.MESH_MIN_FMAX, 120)


if __name__ == "__main__":
    unittest.main(verbosity=2)
