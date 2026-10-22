"""[0071] S8.2 的 2D 单谷 DPT：面内各向异性谷的表头 m* 改用面内二次型拟合的 m_d；面内拟合的 cond 按列归一化。

① 表头 m* 以前取方向分解"最轻两支"的几何平均。面内各向异性谷里最轻的两支都在轻轴上（+x 与 −x），取到的是轻轴质量：
  磷烯型 0.17/1.1 给 0.170，m_d 应为 0.432，2D 表头 μ ∝ 1/m*² 虚高约 6 倍。带边小群不强制面内各向同性时改用
  by_direction 的 m_d = √(m_x·m_y)（与 [0025] 的 3D 同一口径）；六方 K 这类各向同性谷表头不变。
② _quad_fit_2d 的 cond 不归一化时随网格加密变大（Δk²、Δk、1 量纲不同），密网格 + 线性项被误判病态、直接报错；
  归一化后只反映取点几何（十字形取点定不了 xy 交叉项，仍拒绝）。需要 pymatgen + spglib。
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


def _D():
    import gen_step12_dpt as D
    return D


def _rect():
    return Structure(Lattice.orthorhombic(3.30, 4.62, 20.0), ["P"], [[0, 0, 0.5]])          # 磷烯型矩形胞


def _hex():
    return Structure(Lattice.hexagonal(3.19, 20.0), ["Mo"], [[0, 0, 0.5]])                  # MoS2 型六方胞


def _kane(q, m, alpha):
    Q = HB * (q[:, 0] ** 2 / m[0] + q[:, 1] ** 2 / m[1])
    return (np.sqrt(1 + 4 * alpha * Q) - 1) / (2 * alpha)


def build_2d(cwd, st, mesh, valleys, m, alpha=0.3):
    """S3 假 vasprun（IBZ 点 + 能量）：导带在 valleys（分数坐标）各有一个谷，主轴沿 x/y，质量 m=(m_x, m_y)，Kane α；
    价带顶在 Γ（各向同性 0.5）。"""
    import spglib
    mp, grid = spglib.get_ir_reciprocal_mesh(list(mesh), (st.lattice.matrix, st.frac_coords, [1] * len(st)),
                                             is_shift=[0, 0, 0])
    kf = grid[np.unique(mp)] / np.array(mesh, float)
    kf -= np.rint(kf)
    B = st.lattice.reciprocal_lattice.matrix
    cbs = []
    for kv in valleys:
        d = kf - np.asarray(kv, float)
        d -= np.rint(d)
        cbs.append(1.0 + _kane(d @ B, m, alpha))
    cb = np.min(np.stack(cbs), axis=0)
    vb = -HB * np.sum((kf @ B) ** 2, axis=1) / 0.5
    ene = np.stack([vb - 1.0, vb, cb], axis=1)
    occ = np.stack([np.ones(len(kf)), np.ones(len(kf)), np.zeros(len(kf))], axis=1)
    (cwd / "step3_uniform").mkdir(parents=True)
    (cwd / "step3_uniform" / "vasprun.xml").write_text("<x/>")
    return types.SimpleNamespace(final_structure=st, actual_kpoints=kf.tolist(), parameters={"ISPIN": 1},
                                 eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})


@unittest.skipUnless(HAVE_PMG, "需要 pymatgen + spglib")
class SiteIsotropy2DTests(unittest.TestCase):
    def test_inplane_isotropy(self):
        D = _D()
        hx, rc = _hex(), _rect()
        sq = Structure(Lattice.tetragonal(3.5, 20.0), ["Si"], [[0, 0, 0.5]])
        cases = [(hx, [1 / 3, 1 / 3, 0], True),           # 六方 K：C₃
                 (hx, [0, 0, 0], True),                    # 六方 Γ
                 (hx, [0.17, 0.17, 0], False),             # TMD Q 谷
                 (rc, [0, 0, 0], False),                   # 矩形胞 Γ：x、y 不等价
                 (sq, [0, 0, 0], True)]                    # 四方 Γ：C₄
        for st, kf, want in cases:
            self.assertEqual(D._site_isotropic(np.array(kf), st, D._cart_rots(st), dim=2), want, (st.lattice.abc, kf))


@unittest.skipUnless(HAVE_PMG, "需要 pymatgen + spglib")
class HeadlineMass2DTests(unittest.TestCase):
    def _run(self, st, mesh, valleys, m, manual_aniso=None):
        D = _D()
        cwd = Path(tempfile.mkdtemp())
        fake = build_2d(cwd, st, mesh, valleys, m)
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                mock.patch.object(D, "get_C", lambda c, is2d: (100.0, "N/m", "test")), \
                mock.patch.object(D, "get_C_aniso", lambda c, is2d: ((100.0, 100.0), "test")), \
                mock.patch.object(D, "get_E1", lambda c, car: (5.0, "test")), \
                mock.patch.object(D, "get_E1_aniso", lambda c, car: ((5.0, 5.0), "test")), \
                mock.patch.object(D, "_multi_valley_block", lambda *a, **k: {"status": "test"}), \
                mock.patch.object(D, "_single_valley_scope", lambda *a, **k: None), \
                mock.patch.dict(D.MANUAL_ANISO, {"m_eff_electron_xy": manual_aniso}), \
                contextlib.redirect_stdout(io.StringIO()):
            return D, D._one_carrier(cwd, True, "electron", 300.0)

    def test_anisotropic_valley_uses_md(self):
        D, rec = self._run(_rect(), (36, 26, 1), [[0, 0, 0]], (0.17, 1.1))
        bd = rec["by_direction"]
        self.assertEqual(bd["status"], "ok")
        md = np.sqrt(0.17 * 1.1)
        self.assertEqual(rec["inputs"]["m_eff_m0"], bd["m_d_m0"])
        self.assertLess(abs(rec["inputs"]["m_eff_m0"] / md - 1), 0.05, rec["inputs"]["m_eff_m0"])
        self.assertIn("面内二次型拟合 m_d", rec["inputs"]["m_provenance"])
        self.assertAlmostEqual(rec["mobility_cm2_Vs"],
                               D.mobility_dpt(True, 100.0, bd["m_d_m0"], 5.0, 300.0), places=2)

    def test_isotropic_k_valley_unchanged(self):
        """六方 K 谷：小群 C₃ 强制面内各向同性，表头仍是方向分解值。"""
        D, rec = self._run(_hex(), (36, 36, 1), [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]], (0.45, 0.45))
        self.assertIn("方向分解", rec["inputs"]["m_provenance"])
        self.assertNotIn("面内二次型拟合 m_d", rec["inputs"]["m_provenance"])
        self.assertAlmostEqual(rec["inputs"]["m_eff_m0"], 0.45, delta=0.01)

    def test_manual_aniso_mass_keeps_headline(self):
        """手填了分方向 m*：没做面内拟合、判不了小群，表头照旧（不拿手填的 m_d 替换自动值）。"""
        D, rec = self._run(_rect(), (36, 26, 1), [[0, 0, 0]], (0.17, 1.1), manual_aniso=(0.2, 0.9))
        self.assertNotIn("面内二次型拟合 m_d", rec["inputs"]["m_provenance"])


@unittest.skipUnless(HAVE_PMG, "需要 pymatgen + spglib")
class Cond2DTests(unittest.TestCase):
    def _grid(self, st, N, kv, m, th=0.3, alpha=0.5):
        B = st.lattice.reciprocal_lattice.matrix
        g = np.arange(N) / N
        kf = np.array([[x, y, 0.0] for x in g for y in g])
        kf -= np.round(kf)
        d = kf - np.asarray(kv, float)
        d -= np.round(d)
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        e = _kane((d @ B)[:, :2] @ R, m, alpha)
        return kf, B, e, int(np.argmin(e))

    def test_fine_grid_with_linear_terms_not_rejected(self):
        """矩形胞一般位置（带线性项）：48、72 网格以前 cond 1.5e5 / 7.7e5 被拒；现在出数，主质量误差 < 6%。"""
        D = _D()
        for N in (48, 72):
            kf, B, e, k0 = self._grid(_rect(), N, [0.137, 0.21, 0], (0.3, 0.9))
            fit, err = D._quad_fit_2d(kf, B, e, k0, need_lin=True, mesh=np.array([N, N, 1]))
            self.assertIsNone(err, (N, err))
            H = np.array([[2 * fit["a"], fit["c"]], [fit["c"], 2 * fit["b"]]])
            got = sorted(2 * HB / np.linalg.eigvalsh(H))
            self.assertLess(max(abs(g / t - 1) for g, t in zip(got, (0.3, 0.9))), 0.06, (N, got))
            self.assertLess(fit["cond"], 1e4)

    def test_cross_shaped_sampling_still_rejected(self):
        """只在两条坐标轴上取到点：xy 交叉项完全无约束，归一化后 cond 仍是无穷大 -> 拒绝。"""
        D = _D()
        N = 36
        kf, B, e, k0 = self._grid(_rect(), N, [0, 0, 0], (0.3, 0.9), th=0.0)
        on_axes = (np.abs(kf[:, 0]) < 1e-9) | (np.abs(kf[:, 1]) < 1e-9)
        kf, e = kf[on_axes], e[on_axes]
        k0 = int(np.argmin(e))
        fit, err = D._quad_fit_2d(kf, B, e, k0, need_lin=False, mesh=np.array([N, N, 1]))
        self.assertIsNone(fit)
        self.assertIn("病态", err)


if __name__ == "__main__":
    unittest.main()
