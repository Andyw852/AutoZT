#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_symmetry_gate.py —— ke_common.symmetry_gate 的自检（**不依赖任何本地数据**）。

用法：python test_symmetry_gate.py     退出码 0 = 全部 PASS。

★ 2026-09-28 换判据：旧规则"有反演或 τ≡0 就可走 IBZ"两个方向都错（见 V115）。
新判据是**逐操作**的：非 TR 操作要 (R+I)τ≡0，TR 操作要 (R−I)τ≡0（mod 1）。
下面的期望值全部来自独立检验 tmp/amset_desym_phase_check.py 的逐操作实测
（模型势直接对角化 + 与 AMSET 去对称化结果比较），12 个构型 100% 吻合。

注意：结构**在测试里构造**，不读 /mnt/d 下的任何文件（那些路径会被推到公开仓库）。
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
import ke_common as kc


def shifted(st, s):
    x = st.copy()
    x.translate_sites(range(len(x)), s, frac_coords=True)
    return x


def gan_std():
    """GaN P6_3mc，6_3 轴过原点（标准设置）—— 原公式 0/24 错，IBZ 精确。"""
    return Structure(Lattice.hexagonal(3.19, 5.19), ["Ga", "Ga", "N", "N"],
                     [[1/3, 2/3, 0], [2/3, 1/3, .5], [1/3, 2/3, .377], [2/3, 1/3, .877]])


def mos2_aligned():
    return Structure(Lattice.hexagonal(3.19, 12.0), ["Mo", "S", "S"],
                     [[0, 0, .5], [1/3, 2/3, .63], [1/3, 2/3, .37]])


def si_on_atom():
    return Structure(Lattice([[0, 2.715, 2.715], [2.715, 0, 2.715], [2.715, 2.715, 0]]),
                     ["Si", "Si"], [[0, 0, 0], [.25, .25, .25]])


def hcp_std():
    return Structure(Lattice.hexagonal(3.21, 5.21), ["Mg", "Mg"],
                     [[1/3, 2/3, .25], [2/3, 1/3, .75]])


def zincblende_std():
    return Structure(Lattice([[0, 2.825, 2.825], [2.825, 0, 2.825], [2.825, 2.825, 0]]),
                     ["Ga", "N"], [[0, 0, 0], [.25, .25, .25]])


# (名字, 结构, bad_ops 期望, needs_full_grid 期望) —— bad_ops 来自独立实测
CASES = [
    ("GaN 标准原点",       gan_std(),                          0,  False),
    ("GaN 原点在 Ga 上",   shifted(gan_std(), [-1/3, -2/3, 0]), 6,  True),
    ("GaN 任意平移",       shifted(gan_std(), [.13, .07, .2]), 15,  True),
    ("MoS2 对齐后",        mos2_aligned(),                     0,  False),
    ("MoS2 对齐前",        shifted(mos2_aligned(), [.1, .2, 0]), 12, True),
    ("Si 原点在原子",      si_on_atom(),                       0,  False),
    # ★ 关键回归：有反演的 Si，原点落在反演中心时**原公式错 24/48**，旧规则却放行 IBZ
    ("Si 原点在反演中心",  shifted(si_on_atom(), [-.125] * 3), 24,  True),
    ("Si 任意平移",        shifted(si_on_atom(), [.11, .05, .19]), 28, True),
    ("hcp 标准原点",       hcp_std(),                          0,  False),
    # ★ 关键回归：有反演但任意平移 -> 错 8/24，旧规则同样会放行
    ("hcp 任意平移(有反演)", shifted(hcp_std(), [.12, .31, .07]), 8, True),
    ("闪锌矿 标准原点",    zincblende_std(),                   0,  False),
    ("闪锌矿 任意平移",    shifted(zincblende_std(), [.11, .05, .19]), 37, True),
]


class SymmetryGateTests(unittest.TestCase):
    def test_ground_truth_table(self):
        """12 个构型的 bad_ops 必须与独立检验逐一对上。"""
        for name, st, bad, need in CASES:
            g = kc.symmetry_gate(st)
            self.assertEqual(g["bad_ops"], bad, "%s: %s" % (name, g))
            self.assertEqual(g["needs_full_grid"], need, "%s: %s" % (name, g))

    def test_none_structure_is_conservative(self):
        g = kc.symmetry_gate(None)
        self.assertTrue(g["needs_full_grid"])
        self.assertIsNone(g["bad_ops"])
    
    def test_op_phase_exact_helper(self):
        I = np.eye(3)
        # 非 TR：R = -I 时 (R+I)=0 -> 对任意 tau 都满足
        self.assertTrue(kc.op_phase_exact(-I, np.array([0.3, 0.0, 0.0]), False))
        # 非 TR：R = I 时 (R+I)=2I -> tau=0.3 给出 0.6 mod 1 -> 不满足
        self.assertFalse(kc.op_phase_exact(I, np.array([0.3, 0.0, 0.0]), False))
        self.assertTrue(kc.op_phase_exact(I, np.array([0.5, 0.0, 0.0]), False))
        # TR：R = I 时 (R-I)=0 -> 恒满足
        self.assertTrue(kc.op_phase_exact(I, np.array([0.137, 0.2, 0.31]), True))
        # TR：R = -I 时 (R-I)=-2I -> tau=0.3 -> -0.6 mod 1 不满足
        self.assertFalse(kc.op_phase_exact(-I, np.array([0.3, 0.0, 0.0]), True))

    def test_write_full_grid_marker(self):
        """S3 用它落 full_grid_needed.json —— 键与结论都要对。"""
        for name, st, bad, need in (("GaN标准", gan_std(), 0, False),
                                    ("GaN平移", shifted(gan_std(), [.13, .07, .2]), 15, True)):
            d = Path(tempfile.mkdtemp())
            from pymatgen.io.vasp.inputs import Poscar
            p = d / "POSCAR"
            Poscar(st).write_file(str(p))
            g = kc.write_full_grid_marker(p, d)
            self.assertIsNotNone(g, name)
            self.assertEqual(g["needs_full_grid"], need, g)
            f = d / "full_grid_needed.json"
            self.assertTrue(f.is_file(), name)
            j = json.loads(f.read_text(encoding="utf-8"))
            self.assertEqual(j["needs_full_grid"], need, j)
            self.assertEqual(j["bad_ops"], bad, j)
            self.assertIn("reason", j)

    def test_marker_returns_none_on_bad_input(self):
        self.assertIsNone(kc.write_full_grid_marker(Path("/nonexistent/POSCAR")))

    # ---- [patch_desym_fix] 相位补丁开关与判据联动 --------------------------------------
    def test_desym_fix_releases_ibz_but_keeps_counts(self):
        """补丁开了：放行 IBZ，但 bad_ops/total_ops 照算照打（日志要能看到）。"""
        for name, st, bad, need in CASES:
            g = kc.symmetry_gate(st, desym_fix=True)
            self.assertFalse(g["needs_full_grid"], "%s: %s" % (name, g))
            self.assertEqual(g["needs_full_grid_original"], need, "%s: %s" % (name, g))
            self.assertEqual(g["bad_ops"], bad, "%s: %s" % (name, g))
            self.assertIn("bad_ops=%d/" % bad, kc.gate_log_line(g))

    def test_soc_with_tr_ops_always_full_grid(self):
        """SOC + 无反演（有 TR 操作）：ncl 分支不处理 TR -> 补丁开了也要全网格。"""
        g = kc.symmetry_gate(gan_std(), desym_fix=True, ncl=True)
        self.assertGreater(g["tr_ops"], 0, g)
        self.assertTrue(g["needs_full_grid"], g)
        # 有反演（TR 副本被去重，tr_ops=0）的 SOC 体系：补丁开了可以走 IBZ，但注明未验证
        g = kc.symmetry_gate(shifted(si_on_atom(), [-.125] * 3), desym_fix=True, ncl=True)
        self.assertEqual(g["tr_ops"], 0, g)
        self.assertFalse(g["needs_full_grid"], g)
        self.assertIn("未经实数据验证", g["reason"])

    def test_desym_fix_setting_priority(self):
        import os
        old = os.environ.pop(kc.DESYM_FIX_ENV, None)
        try:
            self.assertEqual(kc.desym_fix_setting(None)[0], kc.DESYM_FIX_DEFAULT)
            os.environ[kc.DESYM_FIX_ENV] = "1"
            self.assertTrue(kc.desym_fix_setting(None)[0])
            self.assertTrue(kc.desym_fix_setting("auto")[0])
            self.assertFalse(kc.desym_fix_setting("off")[0])      # step.conf 显式优先
            os.environ.pop(kc.DESYM_FIX_ENV)
            self.assertTrue(kc.desym_fix_setting("on")[0])
        finally:
            os.environ.pop(kc.DESYM_FIX_ENV, None)
            if old is not None:
                os.environ[kc.DESYM_FIX_ENV] = old
        self.assertIn("export AZ_DESYM_FIX=1", kc.desym_fix_cmd_prefix(True))
        self.assertIn("unset AZ_DESYM_FIX", kc.desym_fix_cmd_prefix(False))

    def test_incar_ncl_detection(self):
        d = Path(tempfile.mkdtemp())
        (d / "step3_uniform").mkdir()
        (d / "step3_uniform" / "INCAR").write_text("ISYM = 2\nLSORBIT = .TRUE.\n")
        self.assertTrue(kc.detect_ncl(d))
        (d / "step3_uniform" / "INCAR").write_text("ISYM = 2\n# LSORBIT = .TRUE.\n")
        self.assertFalse(kc.detect_ncl(d))

    # ---- [patch_symprec_unify] 判据的操作集 == AMSET 的操作集 ------------------------------
    def test_ops_match_amset(self):
        """复刻路径（登录节点没装 amset）与 AMSET 自己的 get_reciprocal_point_group_operations
        必须逐操作相同；没装 amset 时跳过。"""
        try:
            from amset.electronic_structure.symmetry import \
                get_reciprocal_point_group_operations as gpo
        except ImportError:
            self.skipTest("没装 amset")
        for name, st, _bad, _need in CASES:
            rots, taus, _src = kc.amset_symmetry_dataset(st)
            mine = kc.amset_op_set(rots, taus)
            r, t, tr = gpo(st, symprec=kc.AMSET_SYMPREC, time_reversal=True)
            theirs = [((-a if b3 else a).T, (-b2 if b3 else b2), bool(b3))
                      for a, b2, b3 in zip(np.asarray(r, float), np.asarray(t, float), tr)]
            key = lambda o: (o[2], tuple(np.round(o[0], 6).ravel()))   # noqa: E731
            A = sorted(mine, key=key)
            B = sorted(theirs, key=key)
            self.assertEqual(len(A), len(B), name)
            for (Ra, ta, fa), (Rb, tb, fb) in zip(A, B):
                self.assertEqual(fa, fb, name)
                self.assertTrue(np.allclose(Ra, Rb), name)
                self.assertTrue(np.allclose(ta, tb, atol=1e-8), name)

    def test_strict_ops_diagnostic(self):
        """近似对称结构：AMSET(0.01) 认出的操作比严格容差多 -> reason 里要提示。"""
        st = gan_std()
        st.translate_sites([0], [0.001, 0.0, 0.0], frac_coords=True)   # ~3 mÅ 畸变
        g = kc.symmetry_gate(st)
        self.assertIsNotNone(g["strict_ops"])
        self.assertLess(g["strict_ops"], 12, g)
        self.assertIn("严格容差", g["reason"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
