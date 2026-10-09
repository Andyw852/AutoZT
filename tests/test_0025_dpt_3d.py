"""[0025] S8.2 的 3D 分方向 DPT（x/y/z）+ 3D 各向同性值改用三个对角（C33、D_zz 以前漏了）。

合成正交晶体（a=4.0, b=4.6, c=5.2 Å，mmm），Γ 点的抛物导带/价带三个方向质量不同；deformation.h5 在带边处的
3×3 形变势对角、OUTCAR 的 C11/C22/C33 也三个方向不同。期望值按
μ_α = 2√(2π) e ħ⁴ C_α / (3 (k_BT)^{3/2} m_α m_d^{3/2} E1_α²) 独立算（不调被测代码）。纯本地，需要 pymatgen + spglib + h5py。
"""
import contextlib
import importlib.util
import io
import json
import math
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(KE / "step8.2_dpt"), str(KE / "tools")):
    if _d not in sys.path:
        sys.path.insert(0, _d)

try:
    import h5py
    import numpy as np
    import spglib  # noqa: F401
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Spin
    HAVE = True
except ImportError:                                   # pragma: no cover
    HAVE = False

E_C, HBAR, KB, M0 = 1.602176634e-19, 1.054571817e-34, 1.380649e-23, 9.1093837015e-31
HB = 3.80998
MC, MV = (0.25, 0.40, 0.70), (0.50, 0.90, 1.60)       # 导带 / 价带 x, y, z 质量
DC, DV = (5.0, 6.0, 7.0), (2.0, 3.0, 4.0)             # 导带 / 价带 D_xx, D_yy, D_zz（eV）
CKBAR = (2000.0, 2500.0, 3000.0)                      # C11, C22, C33（kBar）= 200 / 250 / 300 GPa


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _D():
    import gen_step12_dpt as D
    return D


def _mu3(C_GPa, m_a, m_d, E1, T=300.0):
    """独立实现：μ_α = 2√(2π) e ħ⁴ C_α / (3 (k_BT)^{3/2} m_α m_d^{3/2} E1_α²)，cm²/Vs。"""
    return (2 * math.sqrt(2 * math.pi) * E_C * HBAR ** 4 * C_GPa * 1e9
            / (3 * (KB * T) ** 1.5 * m_a * M0 * (m_d * M0) ** 1.5 * (E1 * E_C) ** 2)) * 1e4


def build(m, mc=MC, mv=MV, dc=DC, dv=DV, ck=CKBAR, mesh=(14, 12, 10), h5=True):
    """材料目录：S1 的 DIM=3d、S3 vasprun（占位，mock 给内容）、S6 OUTCAR、S7.1 band_edges.json + deformation.h5。"""
    import spglib
    st = Structure(Lattice.orthorhombic(4.0, 4.6, 5.2), ["Si"], [[0, 0, 0]])
    mp, grid = spglib.get_ir_reciprocal_mesh(list(mesh), (st.lattice.matrix, st.frac_coords, [14]),
                                             is_shift=[0, 0, 0])
    kf = grid[np.unique(mp)] / np.array(mesh, float)
    kf -= np.rint(kf)
    q = kf @ st.lattice.reciprocal_lattice.matrix
    cb = 1.5 + HB * sum(q[:, i] ** 2 / mc[i] for i in range(3))
    vb = -HB * sum(q[:, i] ** 2 / mv[i] for i in range(3))
    ene = np.stack([vb - 1.0, vb, cb], axis=1)
    occ = np.stack([np.ones(len(kf)), np.ones(len(kf)), np.zeros(len(kf))], axis=1)
    fake = types.SimpleNamespace(final_structure=st, actual_kpoints=kf.tolist(), parameters={"ISPIN": 1},
                                 eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})
    (m / "step1_opt").mkdir(parents=True)
    (m / "step1_opt" / "workflow_method.txt").write_text("DIM=3d\n")
    (m / "step3_uniform").mkdir()
    (m / "step3_uniform" / "vasprun.xml").write_text("<x/>")
    rows = "\n".join(" XX %s" % " ".join("%.4f" % (ck[i] if i == j and i < 3 else 10.0) for j in range(6))
                     for i in range(6))
    (m / "step6_elastic").mkdir()
    (m / "step6_elastic" / "OUTCAR").write_text(" TOTAL ELASTIC MODULI (kBar)\n Direction XX YY ZZ XY YZ ZX\n"
                                               + rows + "\n")
    d = m / "step7b_deform_read"
    d.mkdir()
    hit = lambda b, bh: {"spin": "up", "band_global": b, "band_h5": bh, "k_h5": 0, "k_frac": [0.0, 0.0, 0.0]}
    (d / "band_edges.json").write_text(json.dumps({
        "electron": {"n_hits": 1, "hits": [hit(2, 1)], "E1_iso_eV": (dc[0] + dc[1]) / 2},
        "hole": {"n_hits": 1, "hits": [hit(1, 0)], "E1_iso_eV": (dv[0] + dv[1]) / 2},
        "vac_align": {"status": "skipped", "reason": "dim=3d"}}))
    if h5:
        t = np.zeros((2, 1, 3, 3))
        t[0, 0] = np.diag(dv) + 0.3 * (1 - np.eye(3))       # 非对角项不该进 E1
        t[1, 0] = np.diag(dc) + 0.5 * (1 - np.eye(3))
        with h5py.File(str(d / "deformation.h5"), "w") as f:
            f["deformation_potentials_up"] = t
    return fake


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class Aniso3DTests(unittest.TestCase):
    def setUp(self):
        self.D = _D()
        self.m = Path(tempfile.mkdtemp())

    def _blk(self, fake, carrier):
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            return self.D._aniso_block(self.m, False, carrier, 300.0)

    def test_xyz_analytic(self):
        fake = build(self.m)
        for car, ms, ds in (("electron", MC, DC), ("hole", MV, DV)):
            bd = self._blk(fake, car)
            self.assertEqual(bd["status"], "ok", bd.get("status"))
            md = (ms[0] * ms[1] * ms[2]) ** (1 / 3)
            self.assertAlmostEqual(bd["m_d_m0"], md, places=3)
            for i, d in enumerate("xyz"):
                self.assertAlmostEqual(bd[d]["m_eff_m0"], ms[i], places=3, msg=(car, d))
                self.assertAlmostEqual(bd[d]["E1_eV"], ds[i], places=4)
                self.assertAlmostEqual(bd[d]["C_3D_GPa"], CKBAR[i] / 10.0, places=3)
                want = _mu3(CKBAR[i] / 10.0, ms[i], md, ds[i])
                self.assertAlmostEqual(bd[d]["mobility_cm2_Vs"] / want, 1.0, places=3, msg=(car, d))
                self.assertAlmostEqual(bd[d]["tau_s"] / (bd[d]["mobility_cm2_Vs"] * 1e-4 * ms[i] * M0 / E_C),
                                       1.0, places=4)

    def test_isotropic_limit_equals_3d_formula(self):
        """三个方向质量、形变势、弹性都相同 -> 与现有 3D 各向同性公式（m*^{5/2}）完全一致。"""
        fake = build(self.m, mc=(0.3, 0.3, 0.3), dc=(5.0, 5.0, 5.0), ck=(2500.0, 2500.0, 2500.0))
        bd = self._blk(fake, "electron")
        iso = self.D.mobility_dpt(False, 250e9, 0.3, 5.0, 300.0)
        for d in "xyz":
            self.assertAlmostEqual(bd[d]["mobility_cm2_Vs"] / iso, 1.0, places=3)

    def test_saddle_refused(self):
        """曲率不全同号（鞍点 / 带边定位错）-> 拒绝；导带正定、价带负定才收。"""
        H = np.diag([2 * HB / 0.3, 2 * HB / 0.4, -2 * HB / 0.5])
        got, why = self.D._mass_from_hessian_3d(H, "electron")
        self.assertIsNone(got)
        self.assertIn("不全是正号", why)
        got, why = self.D._mass_from_hessian_3d(-np.diag([2 * HB / 0.3, 2 * HB / 0.4, 2 * HB / 0.5]), "hole")
        self.assertIsNone(why)
        self.assertEqual(got["m"], (0.3, 0.4, 0.5))

    def test_isotropic_inputs_use_three_diagonals(self):
        """3D 各向同性：C = (C11+C22+C33)/3、E1 = (D_xx+D_yy+D_zz)/3（以前是面内两项）。"""
        build(self.m)
        C, unit, prov = self.D.get_C(self.m, False)
        self.assertAlmostEqual(C, 250e9, delta=1.0)
        self.assertIn("C33", prov)
        e1, prov = self.D.get_E1(self.m, "electron")
        self.assertAlmostEqual(e1, 6.0, places=4)
        self.assertIn("(D_xx+D_yy+D_zz)/3", prov)

    def test_missing_h5(self):
        fake = build(self.m, h5=False)
        e1, prov = self.D.get_E1(self.m, "electron")
        self.assertAlmostEqual(e1, 5.5, places=4)                  # 退回面内平均，并注明
        self.assertIn("非立方体系有偏差", prov)
        bd = self._blk(fake, "electron")
        self.assertIn("E1", bd["status"])

    def test_mobility_vs_T_3d(self):
        fake = build(self.m)
        bd = self._blk(fake, "electron")
        rows = self.D.mobility_vs_T({"mobility_cm2_Vs": 100.0, "by_direction": bd}, False, 300.0, [300.0, 600.0])
        r3, r6 = rows
        self.assertAlmostEqual(r3["mobility_mean_cm2_Vs"],
                               (r3["mobility_x_cm2_Vs"] + r3["mobility_y_cm2_Vs"] + r3["mobility_z_cm2_Vs"]) / 3,
                               places=2)
        self.assertAlmostEqual(r6["mobility_z_cm2_Vs"] / r3["mobility_z_cm2_Vs"], 0.5 ** 1.5, places=3)
        self.assertNotIn("mobility_inplane_cm2_Vs", r3)


ML, MT, XID, XIU = 0.9, 0.2, 2.0, 8.0     # Si 型 X 谷：纵/横质量，Ξ_d、Ξ_u（eV）


def build_cubic_x(m, mesh=(12, 12, 12)):
    """简立方（a=4.0 Å，Pm-3m）导带底在三个等价 X 点（(½,0,0) 等），每个谷沿自己的轴是 m_l、横向 m_t；
    价带顶在 Γ（各向同性 0.5）。h5 带 kpoints/structure，X 谷 D = Ξ_d·I + Ξ_u·e e（e = 谷轴）；
    band_edges.json 只记 (½,0,0) 一个带边点（另外两个要靠对称性在 h5 里找到）。"""
    import spglib
    st = Structure(Lattice.cubic(4.0), ["Si"], [[0, 0, 0]])
    mp, grid = spglib.get_ir_reciprocal_mesh(list(mesh), (st.lattice.matrix, st.frac_coords, [14]),
                                             is_shift=[0, 0, 0])
    kf = grid[np.unique(mp)] / np.array(mesh, float)
    kf -= np.rint(kf)
    B = st.lattice.reciprocal_lattice.matrix
    X = np.eye(3) * 0.5
    cbs = []
    for v in range(3):
        dq = kf - X[v]
        dq -= np.rint(dq)
        q = dq @ B
        ms = [ML if i == v else MT for i in range(3)]
        cbs.append(1.5 + HB * sum(q[:, i] ** 2 / ms[i] for i in range(3)))
    cb = np.min(np.stack(cbs), axis=0)
    q = kf @ B
    vb = -HB * np.sum(q ** 2, axis=1) / 0.5
    ene = np.stack([vb - 1.0, vb, cb], axis=1)
    occ = np.stack([np.ones(len(kf)), np.ones(len(kf)), np.zeros(len(kf))], axis=1)
    fake = types.SimpleNamespace(final_structure=st, actual_kpoints=kf.tolist(), parameters={"ISPIN": 1},
                                 eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})
    (m / "step1_opt").mkdir(parents=True)
    (m / "step1_opt" / "workflow_method.txt").write_text("DIM=3d\n")
    (m / "step3_uniform").mkdir()
    (m / "step3_uniform" / "vasprun.xml").write_text("<x/>")
    rows = "\n".join(" XX %s" % " ".join("%.4f" % (2500.0 if i == j and i < 3 else 10.0) for j in range(6))
                     for i in range(6))
    (m / "step6_elastic").mkdir()
    (m / "step6_elastic" / "OUTCAR").write_text(" TOTAL ELASTIC MODULI (kBar)\n Direction XX YY ZZ XY YZ ZX\n"
                                               + rows + "\n")
    d = m / "step7b_deform_read"
    d.mkdir()
    kp = np.array([[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]], float)
    t = np.zeros((2, 4, 3, 3))
    t[0, 0] = np.eye(3) * 3.0
    for v in range(3):
        t[1, 1 + v] = XID * np.eye(3) + XIU * np.outer(np.eye(3)[v], np.eye(3)[v])
    with h5py.File(str(d / "deformation.h5"), "w") as f:
        f["deformation_potentials_up"] = t
        f["kpoints"] = kp
        f["structure"] = st.to_json()
    hit = lambda b, bh, k, kfr: {"spin": "up", "band_global": b, "band_h5": bh, "k_h5": k, "k_frac": kfr}
    (d / "band_edges.json").write_text(json.dumps({
        "electron": {"n_hits": 1, "hits": [hit(2, 1, 1, [0.5, 0, 0])], "E1_iso_eV": XID + XIU / 2},
        "hole": {"n_hits": 1, "hits": [hit(1, 0, 0, [0, 0, 0])], "E1_iso_eV": 3.0},
        "vac_align": {"status": "skipped", "reason": "dim=3d"}}))
    return fake


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class EquivalentValley3DTests(unittest.TestCase):
    """立方晶体的带边在 X（三个等价谷）：晶体各方向必须相同。只取带边那一个谷会给 x 方向 m_l、y/z 方向 m_t
    （μ 差几倍，立方晶体不可能）；按点群平均后各方向都是电导质量 3/(1/m_l+2/m_t)、E1 = Ξ_d + Ξ_u/3。"""

    def test_cubic_x_valleys_isotropic(self):
        D = _D()
        m = Path(tempfile.mkdtemp())
        fake = build_cubic_x(m)
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            bd = D._aniso_block(m, False, "electron", 300.0)
        self.assertEqual(bd["status"], "ok", bd.get("status"))
        mc = 3.0 / (1.0 / ML + 2.0 / MT)
        md = (ML * MT * MT) ** (1 / 3)
        e1 = XID + XIU / 3.0
        self.assertEqual(bd["n_equiv_valleys"], 3)
        self.assertAlmostEqual(bd["m_d_m0"], md, places=3)
        self.assertEqual(sorted(round(x, 2) for x in bd["m_valley_m0"]), [MT, MT, ML])
        self.assertIn("3 个等价谷", bd["m_provenance"])
        self.assertIn("共 3 个", bd["E1_provenance"])
        want = _mu3(250.0, mc, md, e1)
        for d in "xyz":
            self.assertAlmostEqual(bd[d]["m_eff_m0"], mc, places=3, msg=d)
            self.assertAlmostEqual(bd[d]["E1_eV"], e1, places=4, msg=d)
            self.assertAlmostEqual(bd[d]["mobility_cm2_Vs"] / want, 1.0, places=3, msg=d)
        # 表头各向同性 E1 = 三对角平均（迹/3，与谷取向无关）
        with contextlib.redirect_stdout(io.StringIO()):
            e1h, prov = D.get_E1(m, "electron")
        self.assertAlmostEqual(e1h, (3 * XID + XIU) / 3.0, places=4)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class Main3DTests(unittest.TestCase):
    """main()：3D 的 json/summary 有 x/y/z；表头 m* 用三维拟合的 m_d（旧方向分解在各向异性时偏轻，各向同性时一致）。"""

    def _run(self, **kw):
        import os
        D = _D()
        m = Path(tempfile.mkdtemp())
        fake = build(m, **kw)
        (m / "step3_uniform" / "KPOINTS").write_text("auto\n0\nGamma\n14 12 10\n0 0 0\n")
        cwd = os.getcwd()
        os.chdir(m)
        try:
            with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                    mock.patch.object(D, "TEMPERATURES", "300,600"), contextlib.redirect_stdout(io.StringIO()):
                D.main()
        finally:
            os.chdir(cwd)
        return (json.loads((m / D.OUTDIR_NAME / "dpt_result.json").read_text()),
                (m / D.OUTDIR_NAME / "dpt_summary.txt").read_text())

    def test_main_3d(self):
        j, txt = self._run()
        e = next(r for r in j["results"] if r["carrier"] == "electron")
        md = (MC[0] * MC[1] * MC[2]) ** (1 / 3)
        self.assertAlmostEqual(e["inputs"]["m_eff_m0"], md, places=3)
        self.assertIn("三维二次型拟合 m_d", e["inputs"]["m_provenance"])
        self.assertAlmostEqual(e["mobility_cm2_Vs"] / _D().mobility_dpt(False, 250e9, round(md, 4), 6.0, 300.0),
                               1.0, places=4)
        self.assertEqual(sorted(k for k in ("x", "y", "z") if k in e["by_direction"]), ["x", "y", "z"])
        self.assertIn("mobility_mean_cm2_Vs", e["mobility_vs_T"][0])
        self.assertIn("[z] mu=", txt)
        self.assertIn("(μx+μy+μz)/3", txt)
        self.assertTrue(e["multi_valley"]["status"].startswith("多谷 DPT 只实现 2D"))

    def test_isotropic_headline_unchanged(self):
        """各向同性（立方 Γ 谷）：表头 m* 新旧一致，只是来源变成三维拟合。"""
        j, _ = self._run(mc=(0.3, 0.3, 0.3))
        e = next(r for r in j["results"] if r["carrier"] == "electron")
        self.assertAlmostEqual(e["inputs"]["m_eff_m0"], 0.3, places=3)
        self.assertIn("旧的方向分解值 0.3", e["inputs"]["m_provenance"])


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib + h5py")
class Downstream3DTests(unittest.TestCase):
    BD = {"status": "ok", "m_d_m0": 0.4,
          "x": {"mobility_cm2_Vs": 300.0, "tau_s": 1e-14}, "y": {"mobility_cm2_Vs": 200.0, "tau_s": 1e-14},
          "z": {"mobility_cm2_Vs": 100.0, "tau_s": 1e-14}}

    def test_s83_and_tool_average_xyz(self):
        O = _load("out_0025", KE / "step8.3_output" / "gen_step13_output.py")
        self.assertAlmostEqual(O._dpt_inplane_mu(self.BD, {}), 200.0)
        m = Path(tempfile.mkdtemp())
        (m / "step8.2_dpt").mkdir()
        (m / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({
            "is_2d": False, "temperature_K": 300.0,
            "results": [{"carrier": "electron", "mobility_cm2_Vs": 180.0, "by_direction": self.BD}]}))
        T0, dpt, how = O._dpt_mu_T(m, False)
        self.assertEqual(dpt["electron"], {"x": 300.0, "y": 200.0, "z": 100.0, "mean": 200.0})
        self.assertEqual(how["electron"], "by_direction(3D)")
        import mobility_vs_dpt as T
        self.assertEqual(T.dpt_mu({"by_direction": self.BD}), (200.0, "by_direction"))


if __name__ == "__main__":
    unittest.main()
