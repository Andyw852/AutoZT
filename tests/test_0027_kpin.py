"""[0027] 带边是否被对称性"钉住"按晶体小群判断，不再查 ⅓/¼ 坐标表。

SS 实例：矩形胞，导带底落在网格点 (0, ⅓)、真实极小在 ≈(0, 0.36)。坐标表把 (0, ⅓) 当成高对称点 ->
NSTEP=2 只取到 <18 个点 -> 线性项关掉 -> rel=0.23 拒绝。改按小群判断后带线性项，m* 回到真值。
六方 K 等真正被钉住的点行为不变。纯本地，需要 pymatgen + spglib。
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
for _d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(KE / "step8.2_dpt")):
    if _d not in sys.path:
        sys.path.insert(0, _d)

try:
    import spglib  # noqa: F401
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Spin
    HAVE = True
except ImportError:                                   # pragma: no cover
    HAVE = False

HB = 3.80998
MX, MY, YSTAR = 0.97, 1.06, 0.3636                    # SS 诊断里 NSTEP=3 的量级；真实极小不在网格上


def _D():
    import gen_step12_dpt as D
    return D


def _rect():
    return Structure(Lattice.orthorhombic(4.3, 4.0, 20.0), ["Si"], [[0, 0, 0.5]])


def _hex():
    return Structure(Lattice.hexagonal(3.2, 20.0), ["Mo", "S", "S"],
                     [[0, 0, 0.5], [1 / 3, 2 / 3, 0.58], [1 / 3, 2 / 3, 0.42]])


def _grid(N):
    g = np.array([[i / N[0], j / N[1], 0.0] for i in range(N[0]) for j in range(N[1])])
    return g - np.round(g)


def _band(g, recip, center, m, q4=0.0):
    d = g - np.asarray(center, float)
    d -= np.round(d)
    q = d @ recip
    return 1.5 + HB * (q[:, 0] ** 2 / m[0] + q[:, 1] ** 2 / m[1]) + q4 * (q[:, 0] ** 2 + q[:, 1] ** 2) ** 2


def _masses(fit):
    W = np.array([[2 * fit["a"], fit["c"]], [fit["c"], 2 * fit["b"]]]) / (2 * HB)
    return 1 / W[0, 0], 1 / W[1, 1]


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class PinnedTests(unittest.TestCase):
    def test_rectangular(self):
        D = _D()
        kops = D._mv_kops(_rect(), True)
        for kf, want in (((0, 1 / 3, 0), False), ((1 / 3, 0, 0), False), ((0.25, 0, 0), False),
                         ((0, 0, 0), True), ((0, 0.5, 0), True), ((0.5, 0.5, 0), True)):
            self.assertEqual(D._k_pinned(kf, kops, 2), want, kf)
        # 磁性体系不用时间反演时，Y 点仍被 C₂/镜面钉住
        self.assertTrue(D._k_pinned((0, 0.5, 0), D._mv_kops(_rect(), False), 2))

    def test_hexagonal(self):
        D = _D()
        kops = D._mv_kops(_hex(), True)
        for kf, want in (((1 / 3, 1 / 3, 0), True), ((0.5, 0, 0), True), ((0, 0, 0), True),
                         ((1 / 3, 0, 0), False), ((1 / 6, 1 / 6, 0), False)):
            self.assertEqual(D._k_pinned(kf, kops, 2), want, kf)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class Fit2DTests(unittest.TestCase):
    def test_ss_like_off_grid_valley(self):
        """坐标表（kops=None）在 (0, ⅓) 拒绝；按小群判断带上线性项，质量回到真值（能带完全抛物也一样）。"""
        D = _D()
        st = _rect()
        recip = st.lattice.reciprocal_lattice.matrix
        N = (15, 15)
        g = _grid(N)
        for q4 in (0.0, -3.0):
            e = _band(g, recip, (0, YSTAR, 0), (MX, MY), q4)
            k0 = int(np.argmin(e))
            self.assertTrue(np.allclose(g[k0][:2], (0, 1 / 3)))
            with contextlib.redirect_stdout(io.StringIO()):
                old, oerr = D._quad_fit_2d(g, recip, e, k0, mesh=np.array([N[0], N[1], 1]))
                new, nerr = D._quad_fit_2d(g, recip, e, k0, mesh=np.array([N[0], N[1], 1]),
                                           kops=D._mv_kops(st, True))
            self.assertIsNone(old)
            self.assertIn("残差过大", oerr)
            self.assertIsNone(nerr)
            self.assertIn("d", new["names"])                       # 带线性项
            mx, my = _masses(new)
            self.assertAlmostEqual(mx / MX, 1.0, delta=0.02)
            self.assertAlmostEqual(my / MY, 1.0, delta=0.02)

    def test_hex_K_unchanged(self):
        """六方 K（被 C₃ 钉住）：新旧判断一致，拟合结果逐项相同。"""
        D = _D()
        st = _hex()
        recip = st.lattice.reciprocal_lattice.matrix
        N = (15, 15)
        g = _grid(N)
        e = _band(g, recip, (1 / 3, 1 / 3, 0), (0.5, 0.5), -2.0)
        k0 = int(np.argmin(e))
        with contextlib.redirect_stdout(io.StringIO()):
            old, _ = D._quad_fit_2d(g, recip, e, k0, mesh=np.array([N[0], N[1], 1]))
            new, _ = D._quad_fit_2d(g, recip, e, k0, mesh=np.array([N[0], N[1], 1]), kops=D._mv_kops(st, True))
        for k in ("a", "b", "c", "npt", "rel", "names"):
            self.assertEqual(old[k], new[k], k)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class WiringTests(unittest.TestCase):
    """get_effective_mass_aniso 真的把对称操作传给了拟合（假 Vasprun，矩形 2D，导带底在 (0, ⅓) 附近）。"""

    def _build(self, m):
        import spglib
        st = _rect()
        mesh = (15, 15, 1)
        mp, grid = spglib.get_ir_reciprocal_mesh(list(mesh), (st.lattice.matrix, st.frac_coords, [14]),
                                                 is_shift=[0, 0, 0])
        kf = grid[np.unique(mp)] / np.array(mesh, float)
        kf -= np.rint(kf)
        recip = st.lattice.reciprocal_lattice.matrix
        cb = _band(kf, recip, (0, YSTAR, 0), (MX, MY))
        vb = -(_band(kf, recip, (0, 0, 0), (0.6, 0.7)) - 1.5)
        ene = np.stack([vb - 1.0, vb, cb], axis=1)
        occ = np.stack([np.ones(len(kf)), np.ones(len(kf)), np.zeros(len(kf))], axis=1)
        (m / "step1_opt").mkdir(parents=True)
        (m / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\n")
        (m / "step3_uniform").mkdir()
        (m / "step3_uniform" / "vasprun.xml").write_text("<x/>")
        return types.SimpleNamespace(final_structure=st, actual_kpoints=kf.tolist(), parameters={"ISPIN": 1},
                                     eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})

    def test_aniso_mass(self):
        D = _D()
        m = Path(tempfile.mkdtemp())
        fake = self._build(m)
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            got, prov = D.get_effective_mass_aniso(m, "electron", True)
            with mock.patch.object(D, "_mv_kops", side_effect=RuntimeError("无对称操作")):
                old, oprov = D.get_effective_mass_aniso(m, "electron", True)
        self.assertIsNotNone(got, prov)
        self.assertAlmostEqual(got[0] / MX, 1.0, delta=0.02)
        self.assertAlmostEqual(got[1] / MY, 1.0, delta=0.02)
        self.assertIsNone(old)                                     # 拿不到对称操作时退回坐标表 = 旧行为
        self.assertIn("残差过大", oprov)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class Fit3DTests(unittest.TestCase):
    def test_3d_off_grid_valley(self):
        D = _D()
        st = Structure(Lattice.orthorhombic(4.3, 4.0, 4.6), ["Si"], [[0, 0, 0]])
        recip = st.lattice.reciprocal_lattice.matrix
        N = 15
        g = np.array([[i / N, j / N, k / N] for i in range(N) for j in range(N) for k in range(N)])
        g -= np.round(g)
        d = g - np.array([0, YSTAR, 0])
        d -= np.round(d)
        q = d @ recip
        ms = (0.9, 1.1, 1.3)
        e = 1.5 + HB * sum(q[:, i] ** 2 / ms[i] for i in range(3))
        k0 = int(np.argmin(e))
        self.assertTrue(np.allclose(g[k0], (0, 1 / 3, 0)))
        with contextlib.redirect_stdout(io.StringIO()):
            new, nerr = D._quad_fit_3d(g, recip, e, k0, mesh=np.array([N, N, N]), kops=D._kops_3d(st, True))
        self.assertIsNone(nerr)
        self.assertIn("g", new["names"])                           # 带线性项
        W = np.asarray(new["H"]) / (2 * HB)
        for i in range(3):
            self.assertAlmostEqual(1 / W[i, i] / ms[i], 1.0, delta=0.02)
        self.assertFalse(D._k_pinned((0, 1 / 3, 0), D._kops_3d(st, True), 3))
        self.assertTrue(D._k_pinned((0, 0.5, 0), D._kops_3d(st, True), 3))


if __name__ == "__main__":
    unittest.main()
