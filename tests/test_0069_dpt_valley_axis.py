"""[0069] S8.2 的 3D 多谷 DPT：E1 取带边谷沿谷轴的 LA 形变势 |â·D·â|，C 取沿谷轴的纵向模量 C_L(â)。

以前对等价谷先把张量对角平均、再进 1/E1²：Si 导带（Δ 谷，D = Ξ_d·I + Ξ_u·ââ）成了静水值 Ξ_d + Ξ_u/3
（实跑 2.99 vs 谷轴 8.72 eV），μ 高估 8.5 倍。合成数据：简立方 X 谷（沿坐标轴，幅值文件即可）、面心立方 L 谷
（沿 [111] 型轴，非对角有正负号：带符号文件与幅值文件两条路都要给出 Ξ_d + Ξ_u）；全弹性张量的纵向模量；
简并带边与单谷不受影响；表头取分方向平均；空穴 m* 回退值标"不可靠"。需要 pymatgen + spglib + h5py。
"""
import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(KE / "step8.2_dpt"), str(KE / "tools"),
           str(Path(__file__).resolve().parent)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

try:
    import h5py
    import numpy as np
    import spglib  # noqa: F401
    from pymatgen.core import Lattice, Structure
    HAVE = True
except ImportError:                                   # pragma: no cover
    HAVE = False

if HAVE:
    import test_0025_dpt_3d as T25

C11, C12, C44 = 1650.0, 640.0, 790.0                  # kBar（Si 量级）-> 165 / 64 / 79 GPa


def _D():
    import gen_step12_dpt as D
    return D


def _outcar(d, c11=C11, c12=C12, c44=C44, extra=None):
    """立方弹性张量，VASP 顺序 XX YY ZZ XY YZ ZX（kBar）。extra：{(i, j): v} 覆盖（测试 Voigt 重排）。"""
    M = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            M[i, j] = c11 if i == j else c12
    for i in range(3, 6):
        M[i, i] = c44
    for (i, j), v in (extra or {}).items():
        M[i, j] = M[j, i] = v
    rows = "\n".join(" %s %s" % (n, " ".join("%.4f" % x for x in M[i])) for i, n in enumerate(
        ("XX", "YY", "ZZ", "XY", "YZ", "ZX")))
    (d / "step6_elastic").mkdir(parents=True, exist_ok=True)
    (d / "step6_elastic" / "OUTCAR").write_text(" TOTAL ELASTIC MODULI (kBar)\n Direction XX YY ZZ XY YZ ZX\n"
                                               + rows + "\n")


def build_fcc_L(m, xid=-2.0, xiu=12.0, signed=True, amp_tensor=True, hit_k=2):
    """面心立方原胞，导带底在 4 个 L 点（分数坐标 (½,½,½)、(½,0,0)、(0,½,0)、(0,0,½)，笛卡尔轴 [111] 型）。
    每个 L 点 D = Ξ_d·I + Ξ_u·ââ（â = 该谷轴，非对角带正负号）；带边 hit 落在 hit_k（默认 (½,0,0)，轴 (−1,1,1)/√3）。
    signed：写 deformation_signed.h5（AMSET 约定：−dE/dε，非对角 = 2·D_ij）；amp_tensor：幅值 h5 的非对角是否已 ÷2。"""
    a = 5.43
    st = Structure(Lattice([[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]]), ["Si"], [[0, 0, 0]])
    B = st.lattice.reciprocal_lattice.matrix
    kp = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]], float)
    D = np.zeros((2, len(kp), 3, 3))
    D[0, 0] = 3.0 * np.eye(3)
    for j in range(1, len(kp)):
        ax = kp[j] @ B
        ax /= np.linalg.norm(ax)
        D[1, j] = xid * np.eye(3) + xiu * np.outer(ax, ax)
    off = ~np.eye(3, dtype=bool)
    d = m / "step7b_deform_read"
    d.mkdir(parents=True)
    amp = np.abs(D)
    if not amp_tensor:
        amp[..., off] *= 2.0
    with h5py.File(str(d / "deformation.h5"), "w") as f:
        f["deformation_potentials_up"] = amp
        f["kpoints"] = kp
        f["structure"] = st.to_json()
        if amp_tensor:
            f.attrs["dp_offdiag"] = "tensor"
    if signed:
        S = -D.copy()
        S[..., off] *= 2.0
        with h5py.File(str(d / "deformation_signed.h5"), "w") as f:
            f["signed_up"] = S
            f["kpoints"] = kp
    hit = {"spin": "up", "band_global": 2, "band_h5": 1, "k_h5": hit_k, "k_frac": kp[hit_k].tolist()}
    (d / "band_edges.json").write_text(json.dumps({
        "electron": {"n_hits": 1, "hits": [hit]},
        "hole": {"n_hits": 3, "hits": [dict(hit, band_h5=0, k_h5=0, k_frac=[0, 0, 0])] * 1
                 + [dict(hit, band_h5=b, k_h5=0, k_frac=[0, 0, 0]) for b in (0, 1)]},
        "vac_align": {"status": "skipped", "reason": "dim=3d"}}))
    _outcar(m)
    return st


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class LongitudinalModulusTests(unittest.TestCase):
    def test_cubic_directions(self):
        D = _D()
        m = Path(tempfile.mkdtemp())
        _outcar(m)
        Cv = D.read_elastic_voigt_GPa(m)
        c11, c12, c44 = C11 / 10, C12 / 10, C44 / 10
        self.assertAlmostEqual(D.longitudinal_modulus_GPa(Cv, [1, 0, 0]), c11, places=6)
        self.assertAlmostEqual(D.longitudinal_modulus_GPa(Cv, [0, 0, -2]), c11, places=6)
        self.assertAlmostEqual(D.longitudinal_modulus_GPa(Cv, [1, 1, 1]), (c11 + 2 * c12 + 4 * c44) / 3, places=6)
        self.assertAlmostEqual(D.longitudinal_modulus_GPa(Cv, [1, -1, 0]), (c11 + c12 + 2 * c44) / 2, places=6)

    def test_voigt_reorder(self):
        """VASP 的 XY/YZ/ZX 要重排成 Voigt 的 6/4/5：只给 YZ 行（VASP 第 5 行）一个不同的剪切模量。"""
        D = _D()
        m = Path(tempfile.mkdtemp())
        _outcar(m, extra={(4, 4): 500.0})
        Cv = D.read_elastic_voigt_GPa(m)
        self.assertAlmostEqual(Cv[3][3], 50.0)                     # Voigt 4 = YZ
        self.assertAlmostEqual(Cv[4][4], C44 / 10)
        self.assertAlmostEqual(Cv[5][5], C44 / 10)
        self.assertEqual(D.read_elastic_diag_GPa(m), (C11 / 10,) * 3)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class ValleyAxisTests(unittest.TestCase):
    def test_L_valleys_signed(self):
        D = _D()
        m = Path(tempfile.mkdtemp())
        build_fcc_L(m)
        va, why = D._valley_axis_3d(m, "electron")
        self.assertIsNone(why)
        self.assertAlmostEqual(va["E1_eV"], 10.0, places=4)        # Ξ_d + Ξ_u = −2 + 12
        self.assertEqual(va["n_star"], 4)
        self.assertIn("带符号", va["tensor_source"])
        ax = np.abs(va["axis"])
        self.assertTrue(np.allclose(ax, 1 / math.sqrt(3), atol=1e-3), va["axis"])
        self.assertAlmostEqual(va["C_GPa"], (C11 + 2 * C12 + 4 * C44) / 30.0, places=2)
        self.assertEqual(sorted(round(x, 3) for x in va["D_eig_eV"]), [-2.0, -2.0, 10.0])

    def test_L_valleys_amplitude_only(self):
        """没有带符号文件：幅值丢了非对角的正负号，按 |â| 投影仍给 Ξ_d + Ξ_u（Ξ_d + Ξ_u/3 ≥ 0、Ξ_u > 0 时）。
        幅值文件是 AMSET 原生约定（非对角 = 2·D_ij）时也要先 ÷2。"""
        D = _D()
        for amp_tensor in (True, False):
            m = Path(tempfile.mkdtemp())
            build_fcc_L(m, signed=False, amp_tensor=amp_tensor)
            va, why = D._valley_axis_3d(m, "electron")
            self.assertIsNone(why)
            self.assertAlmostEqual(va["E1_eV"], 10.0, places=4, msg=amp_tensor)
            self.assertIn("幅值", va["tensor_source"])

    def test_degenerate_and_off_skip(self):
        D = _D()
        m = Path(tempfile.mkdtemp())
        build_fcc_L(m)
        va, why = D._valley_axis_3d(m, "hole")
        self.assertIsNone(va)
        self.assertIn("简并", why)
        with mock.patch.object(D, "DPT_VALLEY_AXIS", "off"):
            va, why = D._valley_axis_3d(m, "electron")
        self.assertIsNone(va)
        self.assertIn("off", why)

    def test_single_gamma_valley_unchanged(self):
        """[0025] 的 Γ 点正交晶体：单谷，不走谷轴口径（分方向 E1 仍是对角）。"""
        D = _D()
        m = Path(tempfile.mkdtemp())
        fake = T25.build(m)
        va, why = D._valley_axis_3d(m, "electron")
        self.assertIsNone(va)
        self.assertTrue("单谷" in why or "kpoints" in why, why)
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            bd = D._aniso_block(m, False, "electron", 300.0)
        self.assertNotIn("valley_axis", bd)
        self.assertEqual([bd[d]["E1_eV"] for d in "xyz"], list(T25.DC))


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class XValleyDPTTests(unittest.TestCase):
    """简立方 X 谷（[0025] 的 build_cubic_x：Ξ_d = 2、Ξ_u = 8，C11 = 250 GPa，m_l = 0.9、m_t = 0.2）。"""

    def _run(self, fn):
        D = _D()
        m = Path(tempfile.mkdtemp())
        fake = T25.build_cubic_x(m, mesh=(32, 32, 32))                   # [0070] 网格够细，能量窗口能做收敛检查（不触发 m* 回退标注）
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            return D, fn(D, m)

    def test_by_direction_valley_axis(self):
        D, bd = self._run(lambda D, m: D._aniso_block(m, False, "electron", 300.0))
        self.assertEqual(bd["status"], "ok", bd.get("status"))
        mc = 3.0 / (1.0 / T25.ML + 2.0 / T25.MT)
        md = (T25.ML * T25.MT * T25.MT) ** (1 / 3)
        e1 = T25.XID + T25.XIU
        want = T25._mu3(250.0, mc, md, e1)
        self.assertEqual(bd["valley_axis"]["n_star"], 3)
        self.assertIn("多谷谷轴", bd["E1_provenance"])
        for d in "xyz":
            self.assertAlmostEqual(bd[d]["E1_eV"], e1, places=4, msg=d)
            self.assertAlmostEqual(bd[d]["C_3D_GPa"], 250.0, places=3, msg=d)
            self.assertAlmostEqual(bd[d]["m_eff_m0"], mc, places=3, msg=d)
            self.assertAlmostEqual(bd[d]["mobility_cm2_Vs"] / want, 1.0, places=3, msg=d)
        # 旧口径（静水值）：μ 偏大 (E1_axis / E1_hydro)²
        old = T25.XID + T25.XIU / 3.0
        self.assertAlmostEqual(T25._mu3(250.0, mc, md, old) / want, (e1 / old) ** 2, places=6)

    def test_headline_is_direction_mean(self):
        D, rec = self._run(lambda D, m: D._one_carrier(m, False, "electron", 300.0))
        bd = rec["by_direction"]
        self.assertEqual(rec["status"], "ok")
        self.assertAlmostEqual(rec["mobility_cm2_Vs"],
                               sum(bd[d]["mobility_cm2_Vs"] for d in "xyz") / 3.0, places=2)
        self.assertAlmostEqual(rec["inputs"]["E1_eV"], T25.XID + T25.XIU, places=4)
        self.assertAlmostEqual(rec["inputs"]["C_3D_Pa"], 250e9, delta=1e3)
        self.assertAlmostEqual(rec["inputs"]["m_c_m0"], 3.0 / (1.0 / T25.ML + 2.0 / T25.MT), places=3)
        self.assertIn("多谷谷轴", rec["formula"])
        self.assertNotIn("caveats", rec)

    def test_manual_inputs(self):
        """手填 E1：表头用手填值、各向同性公式，C 仍取谷轴纵向模量；只手填 m*：E1 仍是谷轴值（不退回静水值）。"""
        def fn_e1(D, m):
            with mock.patch.dict(D.MANUAL, {"E1_electron_eV": 5.0}):
                return D._one_carrier(m, False, "electron", 300.0)
        D, rec = self._run(fn_e1)
        self.assertEqual(rec["inputs"]["E1_eV"], 5.0)
        self.assertAlmostEqual(rec["inputs"]["C_3D_Pa"], 250e9, delta=1e3)
        self.assertNotIn("formula", rec)
        self.assertAlmostEqual(rec["mobility_cm2_Vs"] / D.mobility_dpt(False, 250e9, rec["inputs"]["m_eff_m0"], 5.0, 300.0),
                               1.0, places=4)

        def fn_m(D, m):
            with mock.patch.dict(D.MANUAL, {"m_eff_electron": 0.3}):
                return D._one_carrier(m, False, "electron", 300.0)
        D, rec = self._run(fn_m)
        self.assertAlmostEqual(rec["inputs"]["E1_eV"], T25.XID + T25.XIU, places=4)
        self.assertAlmostEqual(rec["mobility_cm2_Vs"] / D.mobility_dpt(False, 250e9, 0.3, T25.XID + T25.XIU, 300.0),
                               1.0, places=4)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class CaveatTests(unittest.TestCase):
    def test_mass_fallback_marked_unreliable(self):
        """三维二次型拟合失败（Si 价带顶 rel = 0.44 那种）：表头 m* 是回退值，结果里明确标出来。"""
        D = _D()
        m = Path(tempfile.mkdtemp())
        T25.build(m)
        bd = {"status": "缺输入：m*", "m_provenance": "三维二次型拟合残差过大 rel=0.44（带边附近非二次型）"}

        def fake_mass(cwd, carrier, is_2d):
            D._EDGE_DEGEN[carrier] = 3
            return 0.5, "带边方向分解"
        with mock.patch.object(D, "get_effective_mass", fake_mass), mock.patch.object(D, "_aniso_block",
                                                                                     lambda *a, **k: bd), \
                contextlib.redirect_stdout(io.StringIO()):
            rec = D._one_carrier(m, False, "hole", 300.0)
        cav = " ".join(rec["caveats"])
        self.assertIn("3 重简并", cav)
        self.assertIn("回退值", cav)
        self.assertIn("rel=0.44", cav)


if __name__ == "__main__":
    unittest.main()
