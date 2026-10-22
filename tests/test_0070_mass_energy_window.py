"""[0070] S8.2 的 3D 带边质量：各向异性谷用能量窗口取点（|ΔE| ≤ 0.05→0.15 eV，二次型 + 线性），各向同性谷仍用旧的
球形半径 + 各向同性 q⁴。

合成能带用 Kane 非抛物色散 E(1+αE) = Σ ħ²q_i²/2m_i（主轴任意取向、谷心在一般位置或被钉住的区边界）；真值是带边质量。
旧做法对长椭球谷把纵向质量拉小（Si 型 Δ 谷 34³ 偏 24%、Ge 型 L 谷偏 70%，EPM 硅 34³ m_l 0.405 vs 0.873）；
能量窗口在网格分辨得出的情况下误差 ≤ 10%；分辨不出（轻质量谷、粗网格）时拒绝并在 caveats 里标 m* 不可靠，不给错数。
需要 numpy（EPM 用例另需 scipy）；端到端用例需要 pymatgen + spglib。
"""
import contextlib
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(KE / "step8.2_dpt"), str(KE / "tools"),
           str(Path(__file__).resolve().parent)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

try:
    import spglib  # noqa: F401
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Spin
    HAVE_PMG = True
except ImportError:                                   # pragma: no cover
    HAVE_PMG = False

HB = 3.80998
A0 = 5.43
L_FCC = 0.5 * A0 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
B_FCC = 2 * np.pi * np.linalg.inv(L_FCC).T


def _D():
    import gen_step12_dpt as D
    return D


def axis_rot(v):
    """第三主轴对准 v 的正交阵（列 = 主轴）。"""
    v = np.asarray(v, float) / np.linalg.norm(v)
    t = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(v, t)
    e1 /= np.linalg.norm(e1)
    return np.column_stack([e1, np.cross(v, e1), v])


def rand_rot(seed):
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return q


def kane(kc, kv, m, R, alpha, sign):
    d = (kc - kv) @ R
    Q = HB * np.sum(d ** 2 / np.asarray(m), axis=1)
    E = Q if alpha == 0 else (np.sqrt(1 + 4 * alpha * Q) - 1) / (2 * alpha)
    return sign * E


def grid_near(kv_frac, N, B=B_FCC, rmax=0.45):
    g = np.arange(-N // 2, N // 2 + 1) / N
    kk = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T + np.round(np.asarray(kv_frac) * N) / N
    kc = kk @ B
    sel = np.linalg.norm(kc - np.asarray(kv_frac) @ B, axis=1) < rmax
    return kk[sel], kc[sel]


def fit_masses(D, m, kv_frac, R, alpha, sign, N, pinned, window):
    kf, kc = grid_near(kv_frac, N)
    e = kane(kc, np.asarray(kv_frac) @ B_FCC, m, R, alpha, sign)
    k0 = int(np.argmin(sign * e))
    with contextlib.redirect_stdout(io.StringIO()):
        fit, err = D._quad_fit_3d(kf, B_FCC, e, k0, mesh=np.array([N] * 3), need_lin=not pinned, window=window)
    if err:
        return None, err
    mm, merr = D._mass_from_hessian_3d(fit["H"], "electron" if sign > 0 else "hole")
    if merr:
        return None, merr
    return sorted(np.linalg.eigvalsh(np.array(mm["M"]))), fit


def max_err(got, m):
    return max(abs(g / t - 1) for g, t in zip(got, sorted(m)))


DELTA = np.array([0.4265, 0.4265, 0.0])            # Si 型 Δ 谷（不在格点上），谷轴 z
LPT = np.array([0.5, 0.5, 0.5])                     # L 点（区边界，被对称性钉住），谷轴 [111]
GEN = np.array([0.137, 0.291, 0.073])               # 一般位置
CASES = [  # 名称, 质量, 谷心, 主轴, α, 符号, 是否钉住
    ("Si 型 Δ 谷", (0.19, 0.19, 0.92), DELTA, axis_rot([0, 0, 1]), 0.5, +1, False),
    ("Ge 型 L 谷", (0.08, 0.08, 1.6), LPT, axis_rot([1, 1, 1]), 0.5, +1, True),
    ("扁椭球", (1.0, 1.0, 0.2), GEN, rand_rot(1), 0.5, +1, False),
    ("三轴一般位置", (0.2, 0.5, 1.5), GEN, rand_rot(2), 0.5, +1, False),
    ("空穴一般位置", (0.4, 0.4, 2.0), GEN, rand_rot(5), 0.5, -1, False),
]


class EnergyWindowFitTests(unittest.TestCase):
    def test_anisotropic_valleys_recovered(self):
        D = _D()
        for name, m, kv, R, a, sgn, pin in CASES:
            for N in (34, 48):
                got, fit = fit_masses(D, m, kv, R, a, sgn, N, pin, "energy")
                self.assertIsNotNone(got, (name, N, fit))
                self.assertLess(max_err(got, m), 0.10, (name, N, got))
                self.assertLessEqual(fit["Ew"], 0.15)
                self.assertLessEqual(fit["Ew_dev"], 0.15)                  # 与下一档窗口收敛

    def test_old_radial_fit_biased_for_prolate(self):
        """动机：旧做法对长椭球谷的纵向质量明显偏小（34³）。"""
        D = _D()
        for name, m, kv, R, a, sgn, pin in CASES[:2]:
            got, _ = fit_masses(D, m, kv, R, a, sgn, 34, pin, "radial")
            self.assertGreater(abs(got[2] / max(m) - 1), 0.2, (name, got))

    def test_parabolic_exact(self):
        D = _D()
        m = (0.2, 0.5, 1.5)
        for window in ("energy", "radial"):
            got, _ = fit_masses(D, m, GEN, rand_rot(6), 0.0, +1, 34, False, window)
            self.assertLess(max_err(got, m), 1e-3, window)               # 质量张量输出保留 4 位小数

    def test_unresolvable_light_valley_refused(self):
        """PbTe 型轻谷（m_t = 0.024）：0.15 eV 内只有一两个网格点，能量窗口拒绝而不是给错数。"""
        D = _D()
        got, why = fit_masses(D, (0.024, 0.024, 0.24), LPT, axis_rot([1, 1, 1]), 2.0, +1, 34, True, "energy")
        self.assertIsNone(got)
        self.assertIn("能量窗口", why)


@unittest.skipUnless(HAVE_PMG, "需要 pymatgen + spglib")
class SiteSymmetryTests(unittest.TestCase):
    def test_isotropic_sites(self):
        D = _D()
        si = Structure(Lattice(L_FCC), ["Si", "Si"], [[0, 0, 0], [0.25, 0.25, 0.25]])
        rots = D._cart_rots(si)
        self.assertTrue(D._site_isotropic([0, 0, 0], si, rots))           # 立方 Γ
        self.assertFalse(D._site_isotropic(DELTA, si, rots))              # Δ：C4v
        self.assertFalse(D._site_isotropic(LPT, si, rots))                # L：D3d
        tet = Structure(Lattice.tetragonal(4.0, 6.0), ["Si"], [[0, 0, 0]])
        self.assertFalse(D._site_isotropic([0, 0, 0], tet, D._cart_rots(tet)))   # 四方 Γ：单轴


class EPMSiTests(unittest.TestCase):
    def test_epm_si_conduction_band(self):
        """EPM 硅导带（真值 m_l = 0.873、m_t = 0.184）在 34³ 网格上：能量窗口（限本谷区域）≤ 5%，旧做法 m_l ≈ 0.40。
        Δ 谷离 X 只有 0.15·2π/a：不限本谷区域时 X 另一侧（相邻等价谷）的点进窗口，m_l 拟合成 6–19。"""
        try:
            import _epm_si as E
        except ImportError:                                               # pragma: no cover
            self.skipTest("没有 tests/_epm_si.py")
        D = _D()
        B = 2 * np.pi * np.linalg.inv(E.L0).T
        kv = np.array([0.4265, 0.4265, 0.0])
        g = np.arange(-17, 18) / 34.0
        kk = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T + np.round(kv * 34) / 34
        kc = kk @ B
        sel = np.linalg.norm(kc - kv @ B, axis=1) < 0.36                  # 包含 X 另一侧（相邻等价谷）的点
        kf = kk[sel]
        e = np.array([E.solve(k, E.L0, 80.0)[1][E.NOCC] for k in kf])
        k0 = int(np.argmin(e))
        if not HAVE_PMG:
            self.skipTest("需要 pymatgen（等价谷）")
        st = Structure(Lattice(E.L0), ["Si", "Si"], [E.TAU, 1 - E.TAU])
        star = D._star_images(kf[k0], st, D._cart_rots(st))
        self.assertEqual(len(star), 6)
        out = {}
        for window in ("energy", "radial"):
            with contextlib.redirect_stdout(io.StringIO()):
                fit, err = D._quad_fit_3d(kf, B, e, k0, mesh=np.array([34] * 3), need_lin=True, window=window,
                                          star=star)
            self.assertIsNone(err, window)
            mm, _ = D._mass_from_hessian_3d(fit["H"], "electron")
            out[window] = sorted(np.linalg.eigvalsh(np.array(mm["M"])))
        self.assertLess(abs(out["energy"][2] / 0.873 - 1), 0.05, out)
        self.assertLess(abs(out["energy"][0] / 0.184 - 1), 0.10, out)
        self.assertLess(out["radial"][2], 0.5, out)


def build_x_kane(m, N, alpha=0.5, ml=0.9, mt=0.2):
    """简立方（a = 4 Å）导带底在三个等价 X 点，Kane 非抛物、各谷沿自己的轴 m_l、横向 m_t；价带 Γ 各向同性 0.5。"""
    import spglib
    st = Structure(Lattice.cubic(4.0), ["Si"], [[0, 0, 0]])
    mp, grid = spglib.get_ir_reciprocal_mesh([N] * 3, (st.lattice.matrix, st.frac_coords, [14]), is_shift=[0, 0, 0])
    kf = grid[np.unique(mp)] / float(N)
    kf -= np.rint(kf)
    B = st.lattice.reciprocal_lattice.matrix
    cbs = []
    for v in range(3):
        X = np.eye(3)[v] * 0.5
        dq = kf - X
        dq -= np.rint(dq)
        ms = [ml if i == v else mt for i in range(3)]
        cbs.append(1.5 + kane(dq @ B, np.zeros(3), ms, np.eye(3), alpha, +1))
    cb = np.min(np.stack(cbs), axis=0)
    vb = -HB * np.sum((kf @ B) ** 2, axis=1) / 0.5
    ene = np.stack([vb - 1.0, vb, cb], axis=1)
    occ = np.stack([np.ones(len(kf)), np.ones(len(kf)), np.zeros(len(kf))], axis=1)
    (m / "step3_uniform").mkdir(parents=True)
    (m / "step3_uniform" / "vasprun.xml").write_text("<x/>")
    return types.SimpleNamespace(final_structure=st, actual_kpoints=kf.tolist(), parameters={"ISPIN": 1},
                                 eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})


@unittest.skipUnless(HAVE_PMG, "需要 pymatgen + spglib")
class EndToEndTests(unittest.TestCase):
    def _mass(self, N):
        D = _D()
        m = Path(tempfile.mkdtemp())
        fake = build_x_kane(m, N)
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            return D.get_effective_mass_3d(m, "electron")

    def test_resolved_grid_uses_energy_window(self):
        mres, prov = self._mass(40)
        self.assertIsNotNone(mres, prov)
        self.assertIn("能量窗口", prov)
        self.assertNotIn("caveat", mres)
        mv = sorted(mres["m_valley"])
        self.assertLess(abs(mv[2] / 0.9 - 1), 0.10, mv)
        self.assertLess(abs(mv[0] / 0.2 - 1), 0.10, mv)

    def test_coarse_grid_falls_back_with_caveat(self):
        mres, prov = self._mass(12)
        self.assertIsNotNone(mres, prov)
        self.assertIn("回退", prov)
        self.assertIn("m* 不可靠", mres["caveat"])

    def test_single_window_cannot_check_convergence(self):
        """只剩一个窗口可用时不拿它当结果（没法验收敛），原因写明是网格偏粗，不是"互相不收敛"。"""
        import test_0025_dpt_3d as T25
        D = _D()
        m = Path(tempfile.mkdtemp()) / "m"
        fake = T25.build_cubic_x(m, mesh=(24, 24, 24))                  # 纯抛物 X 谷：基区内只有 0.2 eV 窗口够点
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            mres, prov = D.get_effective_mass_3d(m, "electron")
        self.assertIsNotNone(mres, prov)
        self.assertIn("一个能量窗口可用", mres["caveat"])
        self.assertNotIn("互相不收敛", mres["caveat"])

    def test_isotropic_gamma_valley_keeps_radial_fit(self):
        """价带顶在立方 Γ：小群强制各向同性，仍走球形半径 + q⁴（不加 caveat）。"""
        D = _D()
        m = Path(tempfile.mkdtemp())
        fake = build_x_kane(m, 24)
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            mres, prov = D.get_effective_mass_3d(m, "hole")
        self.assertIsNotNone(mres, prov)
        self.assertIn("各向同性谷", prov)
        self.assertNotIn("caveat", mres)
        self.assertAlmostEqual(mres["m_d"], 0.5, places=2)


if __name__ == "__main__":
    unittest.main()
