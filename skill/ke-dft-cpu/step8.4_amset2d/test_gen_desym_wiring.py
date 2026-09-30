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

    def test_2d_mesh_density_floor(self):
        """patch_mesh_density（V118）：面内下限 = min(263, ceil(|b|/0.00865))，只降不升。"""
        from pymatgen.core import Lattice, Structure
        def st(a):
            return Structure(Lattice.hexagonal(a, 20.0), ["Mo"], [[0, 0, .5]])
        dk = G14.MESH_MIN_DK
        self.assertAlmostEqual(dk, 0.00865)
        self.assertEqual(G14.mesh_need_by_density(st(3.19), [263, 263, 3], dk), [263, 263, 3])  # V112 MoS2
        self.assertEqual(G14.mesh_need_by_density(st(2.46), [263, 263, 3], dk), [263, 263, 3])  # 小晶胞不升
        self.assertEqual(G14.mesh_need_by_density(st(16.0), [263, 263, 3], dk), [53, 53, 3])    # 大晶胞
        self.assertEqual(G14.mesh_need_by_density(st(16.0), [263, 263, 3], 0), [263, 263, 3])   # 关掉
        self.assertEqual(G14.mesh_need_by_density(st(16.0), [0, 0, 3], dk), [0, 0, 3])          # MESH_MIN=0
        self.assertIn("MESH_MIN_DK", G14.SPEC)


class TensorSymmetryTests(unittest.TestCase):
    """patch_tensor_symmetry（V120）：六方 MoS2 的弹性张量带 9% C11/C22 噪声 -> 对称化后 C11=C22、
    C66=(C11-C12)/2、C16=0；对称输入不变；MANUAL_ELASTIC 不动；两个 gen 都在 read_elastic 之后调用。"""

    def _hex_dir(self):
        from pymatgen.core import Lattice, Structure
        st = Structure(Lattice.hexagonal(3.178, 34.879), ["Mo", "S", "S"],
                       [[0, 0, .5], [1/3, 2/3, .545], [1/3, 2/3, .455]])
        d = Path(tempfile.mkdtemp())
        (d / "step6_elastic").mkdir()
        Poscar(st).write_file(str(d / "step6_elastic" / "POSCAR"))
        return d, st

    def test_noisy_hexagonal_elastic_is_projected(self):
        import numpy as np
        import ke_common as kc
        d, st = self._hex_dir()
        C = np.zeros((6, 6))
        C[0, 0], C[1, 1], C[0, 1] = 45.0, 45.0 * 0.91, 11.0
        C[1, 0] = C[0, 1]
        C[5, 5], C[2, 2], C[3, 3], C[4, 4] = 17.8, 0.3, 0.05, 0.05
        Cs = np.asarray(kc.symmetrize_elastic_for(d, C.tolist(), ("step6_elastic",)))
        self.assertAlmostEqual(Cs[0, 0], Cs[1, 1], places=3)
        self.assertAlmostEqual(Cs[5, 5], (Cs[0, 0] - Cs[0, 1]) / 2, places=3)
        again = np.asarray(kc.symmetrize_elastic_for(d, Cs.tolist(), ("step6_elastic",)))
        self.assertTrue(np.allclose(again, Cs, atol=2e-3))                 # 对称输入不变
        off = kc.symmetrize_elastic_for(d, C.tolist(), ("step6_elastic",), enabled=False)
        self.assertEqual(off, C.tolist())                                   # 开关
        self.assertEqual(kc.symmetrize_elastic_for(d, 120.0, ("step6_elastic",)), 120.0)  # 标量不动

    def test_dielectric_asymmetry_diagnostic(self):
        import ke_common as kc
        _d, st = self._hex_dir()
        self.assertLess(kc.rank2_asymmetry([[4, 0, 0], [0, 4, 0], [0, 0, 1.5]], st), 1e-9)
        self.assertGreater(kc.rank2_asymmetry([[4, 0, 0], [0, 3.6, 0], [0, 0, 1.5]], st), 0.04)

    def test_gens_call_symmetrize_after_read_elastic(self):
        for mod in (G10, G14):
            src = Path(mod.__file__).read_text(encoding="utf-8")
            i = src.index("    elastic = read_elastic(cwd)")
            j = src.index("kc.symmetrize_elastic_for(", i)
            self.assertLess(j - i, 600, mod.__name__)
            self.assertIn("if MANUAL_ELASTIC is None:", src[i:j])
            self.assertIn("SYMMETRIZE_ELASTIC", mod.SPEC)
            self.assertTrue(mod.SYMMETRIZE_ELASTIC)


class FullGridWasteGuardTests(unittest.TestCase):
    """patch_full_grid_waste（2026-09-29）：补丁默认开以后，S3b 全网格对非 SOC 材料是白算 ——
    gen 时用同一判据拦下；SOC + TR 操作、补丁被关且有坏操作时照常放行。"""

    @classmethod
    def setUpClass(cls):
        cls.G5B = _load("g5b_under_test", ROOT / "step3b_uniform_full" / "gen_step5b_uniform_full.py")

    def _case(self, st, desym, soc=False):
        d = _matdir(st)
        if soc:
            (d / "step3_uniform" / "INCAR").write_text("LSORBIT = .TRUE.\n")
        return self.G5B._full_grid_needed(d / "step3_uniform" / "POSCAR", d, desym)

    def test_default_on_blocks_full_grid(self):
        import os
        old = os.environ.pop("AZ_DESYM_FIX", None)
        try:
            for st in (tsg.shifted(tsg.si_on_atom(), [-.125] * 3), tsg.gan_std(),
                       tsg.shifted(tsg.mos2_aligned(), [.1, .2, 0])):
                need, why = self._case(st, "auto")                 # 出厂默认 = 开
                self.assertFalse(need, why)
                self.assertIn("DESYM_FIX=on", why)
        finally:
            if old is not None:
                os.environ["AZ_DESYM_FIX"] = old

    def test_still_needed_cases(self):
        need, why = self._case(tsg.shifted(tsg.si_on_atom(), [-.125] * 3), "off")
        self.assertTrue(need, why)                                   # 补丁关 + 24/48 坏操作
        need, why = self._case(tsg.gan_std(), "off")
        self.assertFalse(need, why)                                  # 补丁关但 0 坏操作：也不需要
        need, why = self._case(tsg.shifted(tsg.mos2_aligned(), [.1, .2, 0]), "on", soc=True)
        self.assertTrue(need, why)                                   # SOC + TR 操作：补丁未验证

    def test_force_key_in_spec(self):
        self.assertIn("FULL_GRID_FORCE", self.G5B.SPEC)
        self.assertIn("DESYM_FIX", self.G5B.SPEC)
        self.assertFalse(self.G5B.FULL_GRID_FORCE)


class FullGridFactorTests(unittest.TestCase):
    """patch_full_grid_factor（V124）：vasprun 是全网格时出厂 factor 按 IBZ 口径换算 ——
    全网格 / IBZ 两条路径落在同一张插值网格上；旧行为（全网格直接 f=10）每方向多出一截。
    实例：WSe2/WS2 单层全网格 6075/6627 点 f=10 -> ~299×299×37（3.3M 点），验证密度只要 ~253×253×31。"""

    @staticmethod
    def _grid(mesh, shift=(0, 0, 0)):
        import numpy as np
        g = (np.array(list(np.ndindex(*mesh)), float) + 0.5 * np.array(shift)) / np.array(mesh)
        return g - np.rint(g)                                  # VASP 的 (-0.5, 0.5] 写法

    @staticmethod
    def _ibz(st, mesh):
        import numpy as np
        import spglib
        cell = (st.lattice.matrix, st.frac_coords, [s.Z for s in st.species])
        mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
        return grid[np.unique(mapping)] / np.array(mesh, float)

    def _run(self, mod, st, kpts, **glob):
        import json
        import pymatgen.io.vasp.outputs as O

        class _VR:
            def __init__(self, *a, **k):
                self.actual_kpoints = [list(map(float, x)) for x in kpts]
                self.final_structure = st
        glob.setdefault("INTERPOLATION_FACTOR", 10)
        glob.setdefault("INTERPOLATION_FACTOR_EXPLICIT", False)
        saved = {k: getattr(mod, k) for k in glob}
        old = O.Vasprun
        O.Vasprun = _VR
        try:
            for k, v in glob.items():
                setattr(mod, k, v)
            out = Path(tempfile.mkdtemp())
            f, mesh = mod.apply_mesh_min(out / "vasprun.xml", out)
            return f, mesh, json.loads((out / "interpolation_info.json").read_text())
        finally:
            O.Vasprun = old
            for k, v in saved.items():
                setattr(mod, k, v)

    def test_detect_full_grid(self):
        import ke_common as kc
        st = tsg.mos2_aligned()
        r = kc.full_grid_ibz_count(st, self._grid([12, 12, 3]))
        self.assertEqual(r, (len(self._ibz(st, [12, 12, 3])), 432, [12, 12, 3]))
        self.assertIsNone(kc.full_grid_ibz_count(st, self._ibz(st, [12, 12, 3])))   # IBZ 列表
        self.assertIsNone(kc.full_grid_ibz_count(st, self._grid([12, 12, 3])[1:]))  # 缺一个点
        r = kc.full_grid_ibz_count(tsg.gan_std(), self._grid([4, 4, 4], (1, 1, 1)))  # 平移 MP
        self.assertEqual(r[1:], (64, [4, 4, 4]))

    def _paths(self, mod, st, grid, **glob):
        f_full, m_full, i_full = self._run(mod, st, self._grid(grid), **glob)
        f_ibz, m_ibz, i_ibz = self._run(mod, st, self._ibz(st, grid), **glob)
        self.assertTrue(i_full["full_grid"])
        self.assertEqual(i_full["nk_ir"], len(self._ibz(st, grid)))
        self.assertFalse(i_ibz["full_grid"])
        self.assertEqual(i_ibz["factor_start"], 10)            # IBZ 路径行为不变
        self.assertTrue(i_full["satisfied"] and i_ibz["satisfied"])
        self.assertLess(f_full, 10)
        old = mod._interp_mesh(st, i_full["nk"], 10)            # 旧行为：全网格直接 f=10
        # 两条路径同一张网格（差一两个奇数格点的步长粒度），旧行为每方向多出一截
        self.assertLessEqual(abs(int(m_full[0]) - int(m_ibz[0])), 4, (m_full, m_ibz))
        self.assertGreater(int(old[0]), int(m_full[0]) + 8, (old, m_full))
        return m_full

    def test_2d_full_grid_matches_ibz_path(self):
        m = self._paths(G14, tsg.mos2_aligned(), [12, 12, 3], MESH_MIN=81, MESH_MIN_KZ=3, MESH_MIN_DK=0,
                        MESH_MAX=400, MESH_MIN_FMAX=150)
        self.assertGreaterEqual(int(m[0]), 81)

    def test_3d_full_grid_matches_ibz_path(self):
        m = self._paths(G10, tsg.gan_std(), [6, 6, 4], MESH_MIN=39, MESH_MIN_KZ=None, MESH_MAX=200,
                        MESH_MIN_FMAX=60)
        self.assertGreaterEqual(int(m[0]), 39)

    def test_explicit_factor_not_scaled(self):
        st = tsg.mos2_aligned()
        f, mesh, info = self._run(G14, st, self._grid([12, 12, 3]), INTERPOLATION_FACTOR_EXPLICIT=True,
                                  MESH_MIN=81, MESH_MIN_KZ=3, MESH_MIN_DK=0, MESH_MAX=400, MESH_MIN_FMAX=150)
        self.assertEqual(f, 10)
        self.assertFalse(info["full_grid"])
        self.assertEqual(list(mesh), [int(x) for x in G14._interp_mesh(st, 432, 10)])



class AmsetEnvLazyTests(unittest.TestCase):
    """V129：gen_step14 在 import 时不取 AMSET_ENV（材料目录外 import 不能被 sys.exit 打断，
    本文件与 test_ir_fix_kzcap / test_kernels 都要 import 它）；生成作业时缺 AMSET_ENV 仍报错退出。"""

    def test_import_is_lazy_and_generation_still_strict(self):
        import os
        d = Path(tempfile.mkdtemp())
        old = os.getcwd()
        try:
            os.chdir(d)
            g = _load("g14_env_lazy", HERE / "gen_step14_amset2d.py")
            self.assertFalse(hasattr(g, "AMSET_ENV_NAME"))
            with self.assertRaises(SystemExit):
                g._amset_env()
            (d / "step.conf").write_text("[params]\nAMSET_ENV = amset051\n", encoding="utf-8")
            self.assertEqual(g._amset_env(), "amset051")
        finally:
            os.chdir(old)


if __name__ == "__main__":
    unittest.main(verbosity=2)
