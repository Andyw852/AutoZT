#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_desym_fix.py —— amset_desym_fix 的模型检验（**不需要 VASP、不读任何本地数据**）。

用法：python test_desym_fix.py      退出码 0 = 全部 PASS；没装 amset 时整体 skip。

真值怎么来（与被测公式无关）：每个原子放两个球对称高斯 s 轨道，跳跃只依赖
(元素对, 距离) 的紧束缚哈密顿量 —— 它天然具有晶体的全部空间群对称性。在**每个全网格
k 点上直接对角化**，本征矢乘高斯轨道的解析傅里叶变换得到平面波系数
    c_nk(G) = Σ_a v_a^n(k) · φ̂_a(|k+G|) · e^{−2πi (k+G)·x_a}
—— 不经过任何对称变换公式；c_{k+b}(G)=c_k(G+b) 与 VASP 的 k 周期性一致。
然后走与实数据检验同一条路（tools/desym_fix_validate.validate）：spglib IBZ → 整份传给
AMSET 的 expand_kpoints → AMSET 原式 / 本补丁的修正式展开回全网格 → 与解析真值比。

期望：
  · 修正式在所有构型、所有 k 点上 |cos| ≈ 1（>= 1-1e-8）；
  · 原式出错的 k 点**全部**来自逐操作判据判为"坏"的操作（归因必须干净）；
  · AMSET 自己的操作集上数出来的坏操作数 == ke_common.symmetry_gate 的 bad_ops。
"""
import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT), str(ROOT / "tools"), str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import amset  # noqa: F401
    from pymatgen.electronic_structure.core import Spin
    _HAS_AMSET = True
except Exception:                                            # noqa: BLE001
    _HAS_AMSET = False


def tb_coeffs(st, kpoints, gpoints, rcut=4.0, seed=0):
    """s 轨道紧束缚模型的**本征态**平面波系数与本征值。

    每种元素两个球对称高斯轨道（宽度/在位能按元素定），跳跃只依赖 (元素对, 距离) ——
    哈密顿量天然具有晶体的全部空间群对称性（与原点、是否对齐无关）。
    Bloch 和用约定 I（相位只含格矢 R），所以 c_{k+b}(G) = c_k(G+b)，与 VASP 的 k 点周期性一致；
    本征态满足 O_g ψ_k = ψ_{gk}（简并子空间内差一个幺正变换，比较时用子空间投影）。
    返回 coeffs (nb, nk, ng) 与 energies (nb, nk)（eV 量级）。
    """
    rng = np.random.default_rng(seed)
    species = sorted({s.Z for s in st.species})
    par = {Z: [(rng.uniform(-3, 3), rng.uniform(0.45, 0.8)) for _ in range(2)] for Z in species}
    orbs = [(a, o) for a in range(len(st)) for o in range(2)]
    zs = [st.species[a].Z for a in range(len(st))]
    hop = {}                                                 # (Za,oa,Zb,ob) -> (amp, decay)
    for Za in species:
        for Zb in species:
            for oa in range(2):
                for ob in range(2):
                    key = tuple(sorted([(Za, oa), (Zb, ob)]))
                    if key not in hop:
                        hop[key] = (rng.uniform(0.5, 1.5), rng.uniform(0.8, 1.6))
    lat = st.lattice.matrix
    x = st.frac_coords
    imgs = [np.array(v) for v in np.ndindex(5, 5, 5)]
    imgs = [v - 2 for v in imgs]
    terms = []                                               # (i, j, R, t)
    for i, (a, oa) in enumerate(orbs):
        for j, (b, ob) in enumerate(orbs):
            for R in imgs:
                d = np.linalg.norm((x[b] + R - x[a]) @ lat)
                if d < 1e-6 or d > rcut:
                    continue
                amp, lam = hop[tuple(sorted([(zs[a], oa), (zs[b], ob)]))]
                terms.append((i, j, R, -amp * np.exp(-d / lam)))
    rec = st.lattice.reciprocal_lattice.matrix               # 含 2π，行矢量
    no = len(orbs)
    coeffs = np.zeros((no, len(kpoints), len(gpoints)), complex)
    energies = np.zeros((no, len(kpoints)))
    for ik, k in enumerate(kpoints):
        H = np.diag([par[zs[a]][o][0] for a, o in orbs]).astype(complex)
        for i, j, R, t in terms:
            H[i, j] += t * np.exp(2j * np.pi * np.dot(k, R))
        H = 0.5 * (H + H.conj().T)
        w, v = np.linalg.eigh(H)
        kg = gpoints + k
        q = np.linalg.norm(kg @ rec, axis=1)
        basis = np.stack([np.exp(-0.5 * (q * par[zs[a]][o][1]) ** 2)
                          * np.exp(-2j * np.pi * kg @ x[a]) for a, o in orbs])   # (no, ng)
        coeffs[:, ik] = v.T @ basis
        energies[:, ik] = w
    return coeffs, energies


def gsphere(st, qmax=11.0):
    rec = st.lattice.reciprocal_lattice.matrix
    n = [int(np.ceil(qmax / np.linalg.norm(rec[i]))) + 1 for i in range(3)]
    g = np.array([(i, j, k) for i in range(-n[0], n[0] + 1)
                  for j in range(-n[1], n[1] + 1) for k in range(-n[2], n[2] + 1)])
    return g[np.linalg.norm(g @ rec, axis=1) <= qmax]


def full_mesh(mesh):
    m = np.asarray(mesh)
    grid = np.array(list(np.ndindex(*mesh)), float) / m
    return grid - np.rint(grid)


@unittest.skipUnless(_HAS_AMSET, "需要 amset（集群 AMSET 环境里跑）")
class DesymFixModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import test_symmetry_gate as tsg
        cls.tsg = tsg
        import desym_fix_validate as V
        cls.V = V

    def _run(self, name, st, mesh, expect_bad_ops):
        kpts = full_mesh(mesh)
        gp = gsphere(st)
        cf, en = tb_coeffs(st, kpts, gp)
        rep = self.V.validate({Spin.up: cf}, gp, kpts, st, energies={Spin.up: en},
                              mesh=mesh, verbose=False)
        fx, og = rep["formulas"]["fixed"], rep["formulas"]["orig"]
        msg = "%s: %s" % (name, rep)
        self.assertGreater(fx["min_score"], 1 - 1e-6, msg)
        self.assertEqual(og["bad_but_not_predicted"], 0, msg)
        self.assertTrue(og["attribution_ok"], msg)
        self.assertEqual(rep["ops"]["bad_ops"], expect_bad_ops, msg)
        import ke_common as kc
        self.assertEqual(kc.symmetry_gate(st)["bad_ops"], expect_bad_ops, msg)
        if expect_bad_ops:
            self.assertGreater(og["n_bad_k"], 0, msg)       # 原式确实出错
            self.assertLess(og["min_score"], 0.99, msg)
        else:
            self.assertGreater(og["min_score"], 1 - 1e-6, msg)
        return rep

    # 小网格就够：要测的是每个操作的相位，不是收敛
    MESHES = {"GaN": [4, 4, 3], "MoS2": [6, 6, 1], "Si": [4, 4, 4], "hcp": [4, 4, 3],
              "闪锌矿": [4, 4, 4]}

    def test_all_ground_truth_cases(self):
        """test_symmetry_gate.CASES 的 12 个构型全部过一遍（2026-09-28 本地实测：
        修正式 12/12 最小 |cos| = 1.00000000；原式出错 k 点与预测坏点逐点相同）。"""
        for name, st, bad, _need in self.tsg.CASES:
            mesh = next(m for k, m in self.MESHES.items() if name.startswith(k))
            with self.subTest(name=name):
                self._run(name, st, mesh, bad)

    def test_cli_on_h5(self):
        """实数据入口 main()：写一个真 h5（amset 的 write_coefficients），走完整读盘路径。"""
        import json
        import tempfile
        from amset.wavefunction.io import write_coefficients
        if not hasattr(np, "string_"):          # amset 0.5.1 的 io 用 np.string_（numpy<2）
            np.string_ = np.bytes_
        st = self.tsg.shifted(self.tsg.mos2_aligned(), [.1, .2, 0])
        mesh = [6, 6, 1]
        kp = full_mesh(mesh)
        gp = gsphere(st)
        cf, en = tb_coeffs(st, kp, gp)
        d = Path(tempfile.mkdtemp())
        write_coefficients({Spin.up: cf}, gp, kp, st, filename=str(d / "wavefunction.h5"))
        old = self.V._energies_from_vasprun
        self.V._energies_from_vasprun = lambda *a, **k: {Spin.up: en}
        try:
            rc = self.V.main(["--h5", str(d / "wavefunction.h5"), "--vasprun", "stub",
                              "--json", str(d / "out.json")])
        finally:
            self.V._energies_from_vasprun = old
        self.assertEqual(rc, 0)
        j = json.loads((d / "out.json").read_text(encoding="utf-8"))
        self.assertTrue(j["formulas"]["fixed"]["passed"], j)
        self.assertTrue(j["formulas"]["orig"]["attribution_ok"], j)
        self.assertEqual(j["ops"]["bad_ops"], 12, j)
        self.assertEqual(j["ke_common_gate"]["bad_ops"], 12, j)

    # ---- V116：对照组判据 / 原点平移 ------------------------------------------------
    @staticmethod
    def _noisy(cf, eps, seed=1, sigma=0.0):
        """每个 (带, k) 态独立加相对幅度 eps 的复高斯噪声 —— 模拟 VASP 波函数在对称相关
        k 点间"不精确对称"的数值地板（与公式无关、与原点无关）。sigma>0：各态幅度按
        对数正态抽（异质地板，只有一部分点掉到 0.999 以下，像 MoS2 实测的 1.2%）。"""
        rng = np.random.default_rng(seed)
        n = rng.normal(size=cf.shape) + 1j * rng.normal(size=cf.shape)
        n *= np.linalg.norm(cf, axis=-1, keepdims=True) / np.linalg.norm(n, axis=-1, keepdims=True)
        if sigma:
            eps = eps * np.exp(sigma * rng.normal(size=cf.shape[:2]))[..., None]
        return cf + eps * n

    def test_control_criterion_on_noisy_model(self):
        """带数值地板的数据：严格判据必然 FAIL；对照组判据应 PASS（修正式在坏操作点 ≈ 好操作点），
        而**同一判据套在原式上必须 FAIL**。2D（均匀地板 1-|cos|~2e-3，全部点 < 0.999）与
        3D（异质地板，约一半点 < 0.999）各一例。
        （2026-09-29 本地扫描：两种结构 × 4 种地板 × 10 个种子，修正式 80/80 PASS、原式 80/80 FAIL。）"""
        cases = [("MoS2 2D", self.tsg.shifted(self.tsg.mos2_aligned(), [.1, .2, 0]), [12, 12, 1],
                  dict(eps=0.045)),
                 ("Si 3D", self.tsg.shifted(self.tsg.si_on_atom(), [-.125] * 3), [6, 6, 6],
                  dict(eps=0.015, sigma=0.6))]
        for name, st, mesh, nz in cases:
            with self.subTest(name=name):
                kp = full_mesh(mesh)
                gp = gsphere(st)
                cf, en = tb_coeffs(st, kp, gp)
                rep = self.V.validate({Spin.up: self._noisy(cf, **nz)}, gp, kp, st,
                                      energies={Spin.up: en}, mesh=mesh, verbose=False)
                fx, og, ctl = rep["formulas"]["fixed"], rep["formulas"]["orig"], rep["control"]
                msg = "%s: %s" % (name, {k: rep[k] for k in ("control", "verdict_strict")})
                self.assertLess(fx["min_score"], 0.999, msg)        # 地板确实模拟出来了
                self.assertFalse(rep["verdict_strict"], msg)
                self.assertGreaterEqual(fx["by_class"]["good"]["n"], self.V.CONTROL_MIN_N, msg)
                self.assertGreater(fx["by_class"]["bad"]["n"], 0, msg)
                self.assertTrue(ctl["ok"], msg)
                self.assertIs(ctl["orig_would_pass"], False, msg)   # 判据能区分对错
                self.assertEqual(og["gross_bad_but_not_predicted"], 0, msg)
                self.assertTrue(rep["verdict"], msg)

    def test_shift_origin_is_exact(self):
        """shift_origin(系数, t) 必须等于"平移后结构"的解析真值（符号/约定核对）。"""
        for st, mesh in ((self.tsg.gan_std(), [4, 4, 3]), (self.tsg.mos2_aligned(), [6, 6, 1])):
            t = np.array([.1, .2, .05 if mesh[2] > 1 else 0.0])
            kp = full_mesh(mesh)
            gp = gsphere(st)
            cf, en = tb_coeffs(st, kp, gp)
            c = {Spin.up: cf.copy()}
            st_s = self.V.shift_origin(c, gp, kp, st, t)
            truth, en_s = tb_coeffs(st_s, kp, gp)
            np.testing.assert_allclose(en_s, en, atol=1e-9)
            for ik in range(len(kp)):
                s = self.V.subspace_scores(truth[:, ik], c[Spin.up][:, ik],
                                           self.V.degenerate_groups(en[:, ik]))
                self.assertGreater(s.min(), 1 - 1e-9, "%s k=%s" % (st.formula, kp[ik]))

    def test_shift_makes_null_case_discriminating(self):
        """标准原点 GaN（3D）/ 已对齐 MoS2（2D）：坏操作 0 -> 原样检验无区分力；
        平移原点后同一份数据出现坏操作，修正式逐点分数**平移不变**（含噪声地板），原式不然。"""
        for st, mesh, t in ((self.tsg.gan_std(), [4, 4, 3], [.1, .2, .05]),
                            (self.tsg.mos2_aligned(), [6, 6, 1], [.1, .2, 0])):
            kp = full_mesh(mesh)
            gp = gsphere(st)
            cf, en = tb_coeffs(st, kp, gp)
            with self.subTest(st=st.formula):
                r0 = self.V.validate({Spin.up: cf.copy()}, gp, kp, st, energies={Spin.up: en},
                                     mesh=mesh, verbose=False)
                self.assertEqual(r0["ops"]["bad_ops"], 0)
                self.assertIsNone(r0["control"]["ok"])           # 无区分力
                for noisy in (False, True):
                    c = self._noisy(cf, 0.045) if noisy else cf.copy()
                    rep = self.V.validate_with_shift({Spin.up: c}, gp, kp, st, t,
                                                     energies={Spin.up: en}, mesh=mesh,
                                                     verbose=False)
                    sh = rep["shift"]
                    msg = "%s noisy=%s %s" % (st.formula, noisy, sh)
                    self.assertEqual(sh["bad_ops_before"], 0, msg)
                    self.assertGreater(sh["bad_ops_after"], 0, msg)
                    self.assertTrue(sh["invariant"], msg)
                    self.assertGreater(sh["orig_max_abs_diff"], 0.05, msg)
                    self.assertGreater(rep["formulas"]["orig"]["n_bad_k"], 0, msg)
                    if not noisy:
                        self.assertTrue(rep["verdict_strict"], msg)
                        self.assertTrue(rep["verdict"], msg)

    def test_plugin_switch(self):
        """AZ_DESYM_FIX 未设 -> apply() 不动 AMSET；设了 -> 两个模块都被替换。"""
        import os
        import importlib
        import amset.interpolation.wavefunction as iw
        import amset.wavefunction.common as wc
        import amset_desym_fix as fix
        orig_wc, orig_iw = wc.desymmetrize_coefficients, iw.desymmetrize_coefficients
        old = os.environ.pop("AZ_DESYM_FIX", None)
        try:
            self.assertFalse(fix.apply())
            self.assertIs(wc.desymmetrize_coefficients, orig_wc)
            os.environ["AZ_DESYM_FIX"] = "1"
            self.assertTrue(fix.apply())
            self.assertIs(wc.desymmetrize_coefficients, fix.desymmetrize_coefficients)
            self.assertIs(iw.desymmetrize_coefficients, fix.desymmetrize_coefficients)
        finally:
            wc.desymmetrize_coefficients, iw.desymmetrize_coefficients = orig_wc, orig_iw
            fix.STATE["applied"] = False
            if old is None:
                os.environ.pop("AZ_DESYM_FIX", None)
            else:
                os.environ["AZ_DESYM_FIX"] = old
            importlib.invalidate_caches()


if __name__ == "__main__":
    unittest.main(verbosity=2)
