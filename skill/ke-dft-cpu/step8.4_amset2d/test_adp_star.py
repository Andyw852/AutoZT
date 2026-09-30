#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_adp_star.py —— V127：IR_FIX 下 ADP 的星平均（不需要 VASP）。

用法：python test_adp_star.py      退出码 0 = 全部 PASS。

背景：AMSET 的 ADP 核逐元素取 |D|、q 用分数坐标取整定像、D 在粗网格上线性插值 —— 三者都不随旋转协变，
同一颗星各成员的散射率本来就不等。旧跑法（转置晶格的小群）把成员分别算；IR_FIX 只算代表点，
MoS₂ n 型 ADP 因而 −4.28%。星平均在代表点上对"新群/旧群"的陪集取平均，复现旧跑法。

1) 陪集：六方单层 6 个（D3h×T / {E,σh}×T），纤锌矿 >1，fcc 原胞 1（什么都不做）；AZ_ADP_STAR=0 关。
2) 星平均因子：对星成员 S·k（q 按 AMSET 取整重取像、D 在成员上取值）结果不变；
   等于逐个成员算 AMSET 原因子再平均；标量形变势不动；幂等。
3) 接线：导入 amset_ir_fix 后 scattering_worker / calculate_rate 被包装（子进程会导入本模块）；
   只导入 amset2d_plugin（= spawn 子进程）也会带上 amset_ir_fix。
"""
import json
import os
import subprocess
import sys
import textwrap
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

os.environ.setdefault("AZ_IR_FIX", "0")          # 本进程只用纯函数：不给 AMSET 打补丁
import amset_ir_fix as F  # noqa: E402

try:
    import spglib
    from pymatgen.core import Lattice, Structure
    _HAS = True
except ImportError:                                # noqa: BLE001
    _HAS = False


def mos2(c=20.0):
    z = 1.56 / c
    return Structure(Lattice.hexagonal(3.19, c), ["Mo", "S", "S"],
                     [[0, 0, .5], [1 / 3, 2 / 3, .5 + z], [1 / 3, 2 / 3, .5 - z]])


def gan():
    return Structure(Lattice.hexagonal(3.19, 5.19), ["Ga", "Ga", "N", "N"],
                     [[1 / 3, 2 / 3, 0], [2 / 3, 1 / 3, .5], [1 / 3, 2 / 3, .377], [2 / 3, 1 / 3, .877]])


def si():
    return Structure(Lattice([[0, 2.73, 2.73], [2.73, 0, 2.73], [2.73, 2.73, 0]]), ["Si", "Si"],
                     [[0, 0, 0], [.25, .25, .25]])


def cells(st):
    """(正确 cell, AMSET 原式的转置 cell)"""
    L = np.ascontiguousarray(st.lattice.matrix)
    pos, num = st.frac_coords, [s.Z for s in st.species]
    return (L, pos, num), (np.ascontiguousarray(L.T), pos, num)


def groups(st, symprec=0.01):
    good, bad = cells(st)
    L = good[0]
    G = F.cart_rotations(L, spglib.get_symmetry(good, symprec=symprec)["rotations"])
    H = F.cart_rotations(L, spglib.get_symmetry(bad, symprec=symprec)["rotations"])
    return G, H


@unittest.skipUnless(_HAS, "需要 spglib + pymatgen")
class CosetTests(unittest.TestCase):
    def test_counts(self):
        G, H = groups(mos2())
        self.assertEqual(len(F._uniq(G)), 12)
        self.assertEqual(len(F._uniq(H)), 2)                          # {E, σh}
        self.assertEqual(len(F.coset_rotations(G, H)), 6)
        Gg, Hg = groups(gan())
        self.assertGreater(len(F.coset_rotations(Gg, Hg)), 1)
        Gs, Hs = groups(si())
        self.assertEqual(len(F.coset_rotations(Gs, Hs)), 1)          # 晶格矩阵对称：新旧同群
        for R in G + H:
            self.assertLess(np.abs(R.T @ R - np.eye(3)).max(), 1e-8)

    def test_cosets_cover_group(self):
        G, H = groups(mos2())
        reps = F.coset_rotations(G, H)
        Hp = F._uniq(H + [-h for h in H])
        cover = F._uniq([h @ R for R in reps for h in Hp])
        self.assertEqual(len(cover), len(F._uniq(G + [-g for g in G])))

    def test_publish_env(self):
        env0 = {k: os.environ.get(k) for k in (F.STAR_ENV, F.STAR_SWITCH)}
        try:
            os.environ.pop(F.STAR_SWITCH, None)
            good, bad = cells(mos2())
            F._publish_star_rotations(spglib, good, bad, 0.01, True)
            self.assertEqual(len(F.star_rotations()), 6)
            good, bad = cells(si())
            F._publish_star_rotations(spglib, good, bad, 0.01, True)
            self.assertNotIn(F.STAR_ENV, os.environ)
            os.environ[F.STAR_SWITCH] = "0"
            good, bad = cells(mos2())
            F._publish_star_rotations(spglib, good, bad, 0.01, True)
            self.assertNotIn(F.STAR_ENV, os.environ)
        finally:
            for k, v in env0.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
            F._STAR["loaded"] = False


class _DP:
    """协变的解析形变势场 D(k) = d0·I + d1·g gᵀ，g(k) = 六方第一壳层 cos 和的梯度（笛卡尔）。"""

    def __init__(self, st, d0=-3.0, d1=40.0):
        self.rl = st.lattice.reciprocal_lattice.matrix
        self.R = [np.array([a, b, 0.0]) @ st.lattice.matrix for a, b in ((1, 0), (0, 1), (1, 1))]
        self.d0, self.d1 = d0, d1

    def interpolate(self, spin, bands, kpoints):
        out = []
        for k in np.atleast_2d(kpoints):
            kc = np.asarray(k, float) @ self.rl
            g = -sum(np.sin(kc @ R) * R for R in self.R)
            out.append(self.d0 * np.eye(3) + self.d1 * np.outer(g, g) / 100.0)
        return np.array(out)


class _ADP:
    """AMSET 式 ADP 因子（|D| 逐元素 + v⊗v，纵/横两支、面内各向同性）。"""
    name = "ADP"

    def __init__(self, dp):
        self.deformation_potential = dp
        self.fermi_levels = np.zeros((1, 1))

    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        D = np.abs(self.deformation_potential.interpolate(spin, [band_idx], [kpoint])[0])
        D = D + np.outer(velocity, velocity)
        u = np.asarray(unit_q)
        t = np.cross([0, 0, 1.0], u)
        sl = np.einsum("ni,nj->nij", u, u)
        st = 0.5 * (np.einsum("ni,nj->nij", u, t) + np.einsum("ni,nj->nij", t, u))
        f = np.tensordot(sl, D) ** 2 / 3.0 + np.tensordot(st, D) ** 2 / 1.0
        return f[None, None] * np.ones(self.fermi_levels.shape + np.shape(norm_q_sq))


def _amset_q(k, kprimes, rlat):
    """AMSET 在起点 k 上的 q：分数坐标差取整（pbc_diff），再转笛卡尔。"""
    f = np.asarray(kprimes) - np.asarray(k)
    return (f - np.round(f)) @ rlat


@unittest.skipUnless(_HAS, "需要 spglib + pymatgen")
class StarFactorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.st = mos2()
        cls.rl = cls.st.lattice.reciprocal_lattice.matrix
        G, H = groups(cls.st)
        cls.G, cls.reps = F._uniq(G + [-g for g in G]), F.coset_rotations(G, H)
        rng = np.random.default_rng(7)
        cls.k = np.array([0.29, 0.37, 0.0])                                 # K 谷附近的一般点
        kp = rng.uniform(-0.5, 0.5, (40, 3))
        kp[:, 2] = rng.choice([-1 / 3, 0, 1 / 3], 40)
        cls.kp = kp
        cls.v = np.array([0.03, -0.05, 0.0])

    def _eval(self, adp, k, v):
        q = _amset_q(k, self.kp @ np.eye(3), self.rl)
        n = np.sum(q ** 2, axis=1)
        return adp.factor(q / np.sqrt(n)[:, None], n, "up", 0, k, v)

    def _member(self, S, k):
        kc = (k @ self.rl) @ S.T
        kS = kc @ np.linalg.inv(self.rl)
        return kS - np.round(kS)

    def test_rep_independent(self):
        dp = _DP(self.st)
        a = _ADP(dp)
        self.assertTrue(F.star_average_factor(a, self.reps, self.rl))
        self.assertFalse(F.star_average_factor(a, self.reps, self.rl))          # 幂等
        ref = None
        spreads = []
        for S in self.G:
            kS = self._member(S, self.k)
            kpS = (self.kp @ self.rl) @ S.T @ np.linalg.inv(self.rl)
            q = _amset_q(kS, kpS, self.rl)
            n = np.sum(q ** 2, axis=1)
            out = a.factor(q / np.sqrt(n)[:, None], n, "up", 0, kS, S @ self.v)
            plain = _ADP(dp).factor(q / np.sqrt(n)[:, None], n, "up", 0, kS, S @ self.v)
            spreads.append(plain.sum())
            if ref is None:
                ref = out
            np.testing.assert_allclose(out, ref, rtol=1e-9)
        self.assertGreater((max(spreads) - min(spreads)) / np.mean(spreads), 0.02)   # 原因子确实不协变

    def test_equals_explicit_member_average(self):
        dp = _DP(self.st)
        a = _ADP(dp)
        F.star_average_factor(a, self.reps, self.rl)
        got = self._eval(a, self.k, self.v).sum()
        exp = []
        for S in self.reps:
            kS = self._member(S, self.k)
            kpS = (self.kp @ self.rl) @ S.T @ np.linalg.inv(self.rl)
            q = _amset_q(kS, kpS, self.rl)
            n = np.sum(q ** 2, axis=1)
            exp.append(_ADP(dp).factor(q / np.sqrt(n)[:, None], n, "up", 0, kS, S @ self.v).sum())
        self.assertAlmostEqual(got, np.mean(exp), delta=1e-9 * abs(got))
        self.assertIsNone(dp.__dict__.get("interpolate"))                    # 临时覆盖已撤掉

    def test_scalar_dp_untouched(self):
        class _S:
            name = "ADP"
            deformation_potential = (6.0, 6.0)

            def factor(self, *a):
                return 1.0
        s = _S()
        self.assertFalse(F.star_average_factor(s, self.reps, self.rl))
        self.assertNotIn("factor", s.__dict__)


def _stub_env():
    import importlib.util
    import tempfile
    e = dict(os.environ)
    paths = [str(HERE)]
    if importlib.util.find_spec("sumo") is None:
        d = Path(tempfile.mkdtemp()) / "sumo"
        d.mkdir()
        (d / "__init__.py").write_text("")
        (d / "symmetry.py").write_text("class Kpath: pass\nclass PymatgenKpath(Kpath): pass\n")
        paths.append(str(d.parent))
    e["PYTHONPATH"] = os.pathsep.join(paths + [e.get("PYTHONPATH", "")])
    return e


class WiringTests(unittest.TestCase):
    PRE = textwrap.dedent('''
        import json, os, sys, numpy as np
        for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
            if not hasattr(np, _a):
                setattr(np, _a, getattr(np, _b))
    ''')

    def _run(self, code, env=None, cwd=None):
        e = _stub_env()
        e.update(env or {})
        r = subprocess.run([sys.executable, "-c", self.PRE + textwrap.dedent(code)], capture_output=True,
                           text=True, env=e, cwd=cwd, timeout=600)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        return json.loads(r.stdout.strip().splitlines()[-1]), r.stderr

    def test_3d_worker_wrapped(self):
        out, _ = self._run('''
            import amset_ir_fix
            import amset.scattering.calculate as C
            print(json.dumps({"w": C.scattering_worker.__module__, "cr": getattr(C.calculate_rate, "_az_star", False)}))
        ''', env={"AZ_IR_FIX": "1"})
        self.assertEqual(out, {"w": "amset_ir_fix", "cr": True})

    def test_off_keeps_original(self):
        out, _ = self._run('''
            import amset_ir_fix
            import amset.scattering.calculate as C
            print(json.dumps({"w": C.scattering_worker.__module__, "cr": getattr(C.calculate_rate, "_az_star", False)}))
        ''', env={"AZ_IR_FIX": "0"})
        self.assertEqual(out, {"w": "amset.scattering.calculate", "cr": False})

    def test_plugin_brings_ir_fix_into_child(self):
        import shutil
        import tempfile
        d = Path(tempfile.mkdtemp())
        shutil.copyfile(HERE / "example_2d_correction.json", d / "2d_correction.json")
        out, _ = self._run('''
            import amset2d_plugin                      # spawn 子进程只导入插件
            import amset.scattering.calculate as C
            print(json.dumps({"irfix": "amset_ir_fix" in sys.modules,
                              "cr": getattr(C.calculate_rate, "_az_star", False),
                              "w": C.scattering_worker.__module__}))
        ''', env={"AZ_IR_FIX": "1", "AMSET2D_RECORD": str(d / "2d_correction.json")}, cwd=str(d))
        self.assertEqual(out, {"irfix": True, "cr": True, "w": "amset2d_plugin"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
