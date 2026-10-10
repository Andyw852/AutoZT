"""[0026] 3D 的 zz 方向打通 S8.1 -> S8.3；顺带恢复 0024 的 S8.3 汇总两行。

S8.1 用一个假的 BoltzTraP2（只给出 3×3 输运张量，三个对角各不相同），走真实的 run_boltztrap_crta、
文献口径扫描与 S8.3 的 load_bt2 / build_table / write_table。期望值按"σ_α = (σ/τ)_αα · τ_α、
ZT_α = S_α² σ_α T / (κL_α + L σ_α T)"独立算。2D 的输出（键、列、顺序）必须和以前一样。纯本地。
"""
import contextlib
import csv
import importlib.util
import io
import json
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(ROOT)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

SOT = (1.0e17, 2.0e17, 3.0e17)        # σ/τ 的 xx/yy/zz（S/m/s）
SBK = (2.0e-4, 1.5e-4, 1.0e-4)        # |S| 的 xx/yy/zz（V/K）；电子侧取负
KAP = (1.0e3, 2.0e3, 3.0e3)           # κe/τ
TAU = {"electron": (1.0e-14, 2.0e-14, 4.0e-14), "hole": (3.0e-14, 3.0e-14, 5.0e-14)}
KL = (1.5, 2.0, 0.8)                  # κL xx/yy/zz（W/mK）
L0 = 2.44e-8
NELECT, VOL_BOHR3 = 8.0, 300.0
T300 = 300.0

BT2_2D_KEYS = ["scissor", "temperatures_K", "mu_rel_efermi_eV", "cell_volume_cm3", "nelect",
               "sigma_over_tau_S_per_m_s", "seebeck_V_per_K", "kappa_e_over_tau_W_per_m_K_s",
               "sigma_over_tau_xx", "sigma_over_tau_yy", "seebeck_xx_V_per_K", "seebeck_yy_V_per_K",
               "kappa_e_over_tau_xx", "kappa_e_over_tau_yy", "power_factor_over_tau", "lorenz_WOhm_per_K2",
               "carrier_conc_cm-3", "note"]
SCAN_2D_COLS = ["mu_rel_efermi_eV", "carrier_conc_cm-3", "n_2D_cm-2", "branch",
                "tau_xx_fs", "tau_yy_fs", "S_xx_uV/K", "S_yy_uV/K",
                "sigma_xx_S/m", "sigma_yy_S/m", "sheet_sigma_xx_S", "sheet_sigma_yy_S",
                "PF_xx_W/mK2", "PF_yy_W/mK2", "kappa_e_WF_xx_W/mK", "kappa_e_WF_yy_W/mK"]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _g11():
    return _load("g11_0026", KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py")


def _g13():
    return _load("g13_0026", KE / "step8.3_output" / "gen_step13_output.py")


def _fake_btp2(nT, nmu):
    """假的 BoltzTraP2：只有 run_boltztrap_crta 用到的几个函数。载流子 ±1e17…3e21 cm⁻³（对数间隔）。"""
    pkg = types.ModuleType("BoltzTraP2")
    dft, sphere, fite, bandlib, units = (types.ModuleType("BoltzTraP2." + n)
                                         for n in ("dft", "sphere", "fite", "bandlib", "units"))

    class DFTData:
        def __init__(self, path):
            self.nelect, self.kpoints, self.atoms, self.magmom = NELECT, np.zeros((4, 3)), None, None
            self.fermi, self.dosweight = 0.0, 2.0

        def get_volume(self):
            return VOL_BOHR3

        def get_lattvec(self):
            return np.eye(3)

    dft.DFTData = DFTData
    sphere.get_equivalences = lambda *a, **k: None
    fite.fitde3D = lambda *a, **k: None
    fite.getBTPbands = lambda *a, **k: (None, None, None)
    bandlib.BTPDOS = lambda *a, **k: (None, None, None, None)
    vuc_cm3 = VOL_BOHR3 * (0.529177210903e-8) ** 3
    half = np.logspace(17, 21.5, nmu // 2)
    n_cm3 = np.concatenate([-half[::-1], half])
    N = np.tile(n_cm3 * vuc_cm3 - NELECT, (nT, 1))

    def fermiintegrals(*a, **k):
        return N, None, None, None, None

    def onsager(*a, **k):
        sig = np.zeros((nT, nmu, 3, 3))
        sbk = np.zeros((nT, nmu, 3, 3))
        kap = np.zeros((nT, nmu, 3, 3))
        sgn = np.where(n_cm3 < 0, -1.0, 1.0)
        for i in range(3):
            sig[..., i, i] = SOT[i]
            sbk[..., i, i] = SBK[i] * sgn
            kap[..., i, i] = KAP[i]
        return sig, sbk, kap, None

    bandlib.fermiintegrals = fermiintegrals
    bandlib.calc_Onsager_coefficients = onsager
    for n, m in (("dft", dft), ("sphere", sphere), ("fite", fite), ("bandlib", bandlib), ("units", units)):
        setattr(pkg, n, m)
    return {"BoltzTraP2": pkg, "BoltzTraP2.dft": dft, "BoltzTraP2.sphere": sphere, "BoltzTraP2.fite": fite,
            "BoltzTraP2.bandlib": bandlib, "BoltzTraP2.units": units}


def _run_bt2(B, is_2d, nmu=40):
    from unittest import mock
    with mock.patch.dict(sys.modules, _fake_btp2(1, nmu)), contextlib.redirect_stdout(io.StringIO()):
        return B.run_boltztrap_crta("unused", [T300], 0.5, nmu, 100, 1, nelect=NELECT, is_2d=is_2d)


def _material(tmp, dim):
    """<tmp>/ke（材料目录）+ 同级 kl-dft-cpu 的 kappa_summary.json + S8.2 的 dpt_result.json。"""
    m = Path(tmp) / "ke"
    (m / "step1_opt").mkdir(parents=True)
    (m / "step1_opt" / "workflow_method.txt").write_text("DIM=%s\n" % dim)
    kl = Path(tmp) / "kl-dft-cpu" / "step6_kappa"
    kl.mkdir(parents=True)
    (kl / "kappa_summary.json").write_text(json.dumps(
        {"temperatures": [200.0, 400.0], "kappa_xx_yy_zz": [list(KL), list(KL)]}))
    dirs = "xyz" if dim == "3d" else "xy"
    res = []
    for car in ("electron", "hole"):
        bd = {"status": "ok", "m_d_m0": 0.5}
        for i, d in enumerate(dirs):
            bd[d] = {"mobility_cm2_Vs": 100.0 * (i + 1), "tau_s": TAU[car][i], "m_eff_m0": 0.5}
        res.append({"carrier": car, "mobility_cm2_Vs": 150.0, "inputs": {"m_eff_m0": 0.5}, "by_direction": bd})
    (m / "step8.2_dpt").mkdir()
    (m / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(
        {"temperature_K": T300, "is_2d": dim == "2d", "dim": dim, "results": res}))
    return m


@contextlib.contextmanager
def _cwd(p):
    old = os.getcwd()
    os.chdir(p)
    try:
        yield
    finally:
        os.chdir(old)


def _read_csv(p):
    rows = [ln for ln in Path(p).read_text(encoding="utf-8").splitlines() if not ln.startswith("#")]
    return list(csv.DictReader(rows))


class S81ArraysTests(unittest.TestCase):
    def test_3d_writes_zz_2d_unchanged(self):
        B = _g11()
        r3 = _run_bt2(B, False)
        self.assertEqual(list(r3)[:len(BT2_2D_KEYS)], BT2_2D_KEYS)
        self.assertEqual(list(r3)[len(BT2_2D_KEYS):],
                         ["sigma_over_tau_zz", "seebeck_zz_V_per_K", "kappa_e_over_tau_zz"])
        self.assertAlmostEqual(r3["sigma_over_tau_zz"][0][0], SOT[2])
        self.assertAlmostEqual(r3["kappa_e_over_tau_zz"][0][-1], KAP[2])
        self.assertAlmostEqual(r3["seebeck_zz_V_per_K"][0][0], -SBK[2])          # 电子侧
        self.assertAlmostEqual(r3["sigma_over_tau_S_per_m_s"][0][0], sum(SOT) / 3)  # 3D 平均本来就含 zz
        r2 = _run_bt2(B, True)
        self.assertEqual(list(r2), BT2_2D_KEYS)                                   # 2D：键与顺序都不变
        self.assertAlmostEqual(r2["sigma_over_tau_S_per_m_s"][0][0], (SOT[0] + SOT[1]) / 2)


class S81PaperScanTests(unittest.TestCase):
    def _scan(self, dim):
        B = _g11()
        tmp = tempfile.mkdtemp()
        m = _material(tmp, dim)
        res = _run_bt2(B, dim == "2d")
        B.NORM_2D = "cell"
        B.DOPING_ALIGN = None
        out = m / "step8.1_boltztrap"
        out.mkdir()
        with contextlib.redirect_stdout(io.StringIO()):
            lines = B.write_paper_scan(out, m, res, dim == "2d")
        return B, lines, out / "paper_mu_scan_300K.csv"

    def test_3d_zz_columns_and_values(self):
        B, lines, p = self._scan("3d")
        head = Path(p).read_text(encoding="utf-8").splitlines()[1].split(",")
        for k in ("tau_%s_fs", "S_%s_uV/K", "sigma_%s_S/m", "PF_%s_W/mK2", "kappa_e_WF_%s_W/mK"):
            self.assertEqual(head.index(k % "zz"), head.index(k % "yy") + 1, k)
        self.assertNotIn("sheet_sigma_zz_S", head)
        self.assertEqual(head[-3:], ["ZT_xx", "ZT_yy", "ZT_zz"])
        rows = _read_csv(p)
        e = next(r for r in rows if r["branch"] == "electron")
        h = next(r for r in rows if r["branch"] == "hole")
        for r, car in ((e, "electron"), (h, "hole")):
            sig = SOT[2] * TAU[car][2]
            self.assertAlmostEqual(float(r["tau_zz_fs"]), TAU[car][2] * 1e15, places=6)
            self.assertAlmostEqual(float(r["sigma_zz_S/m"]) / sig, 1.0, places=9)
            ke = L0 * sig * T300
            zt = (SBK[2] ** 2) * sig * T300 / (KL[2] + ke)
            self.assertAlmostEqual(float(r["ZT_zz"]) / zt, 1.0, places=2)        # L 取脚本的 Sommerfeld 值
        txt = "\n".join(lines)
        self.assertIn("z=%.1f fs" % (TAU["electron"][2] * 1e15), txt)
        self.assertIn("ZT_zz", txt)
        self.assertIn("zz=%.4f" % KL[2], txt)

    def test_2d_columns_unchanged(self):
        B, lines, p = self._scan("2d")
        head = Path(p).read_text(encoding="utf-8").splitlines()[1].split(",")
        self.assertEqual(head, SCAN_2D_COLS + ["ZT_xx", "ZT_yy"])
        self.assertNotIn("z=", "\n".join(x for x in lines if x.startswith("# tau_")))
        self.assertNotIn("zz=", "\n".join(lines))


class S83ChainTests(unittest.TestCase):
    """S8.1 写的 boltztrap_crta.json -> S8.3 的 load_bt2 / build_table / write_table。"""

    def _chain(self, dim):
        B, O = _g11(), _g13()
        tmp = tempfile.mkdtemp()
        m = _material(tmp, dim)
        (m / O.BT2_DIR).mkdir()
        (m / O.BT2_DIR / "boltztrap_crta.json").write_text(json.dumps(_run_bt2(B, dim == "2d")))
        out = m / "step8.3_output"
        out.mkdir()
        with _cwd(m), contextlib.redirect_stdout(io.StringIO()):
            bt, dpt = O.load_bt2(m), O.load_dpt(m)
            rows = O.build_table(None, bt, dpt)
            O.write_table(out, rows, None, bt, dpt)
        return rows, _read_csv(out / "comparison_300K.csv")

    def test_3d_zz_columns(self):
        rows, csvrows = self._chain("3d")
        e = next(r for r in rows if r["carrier_conc_cm-3"] == -1e19)
        sig = SOT[2] * TAU["electron"][2]
        self.assertAlmostEqual(e["bt2_sigma_zz_S/m"] / sig, 1.0, places=9)
        self.assertAlmostEqual(e["bt2_S_zz_uV/K"], -SBK[2] * 1e6, places=6)
        zt = (SBK[2] ** 2) * sig * T300 / (KL[2] + e["bt2_kappa_e_WF_zz_W/mK"])
        self.assertAlmostEqual(e["bt2_ZT_zz"], zt, places=12)
        # 3D 的平均 ZT 用三个方向 κL 的平均（S、σ 本来就是三方向平均）
        kl_avg = sum(KL) / 3.0
        zt_avg = (e["bt2_S_uV/K"] * 1e-6) ** 2 * e["bt2_sigma_S/m"] * T300 / (kl_avg + e["bt2_kappa_e_WF_W/mK"])
        self.assertAlmostEqual(e["bt2_ZT"], zt_avg, places=12)
        head = list(csvrows[0])
        for k in ("bt2_S_%s_uV/K", "bt2_sigma_%s_S/m", "bt2_PF_%s_W/mK2", "bt2_kappa_e_WF_%s_W/mK"):
            self.assertEqual(head.index(k % "zz"), head.index(k % "yy") + 1, k)
        self.assertEqual(head.index("bt2_ZT_zz"), head.index("bt2_ZT_yy") + 1)

    def test_2d_no_zz(self):
        rows, csvrows = self._chain("2d")
        self.assertFalse(any("zz" in k for r in rows for k in r))
        self.assertFalse(any("zz" in k for k in csvrows[0]))
        e = next(r for r in rows if r["carrier_conc_cm-3"] == -1e19)
        zt_avg = (e["bt2_S_uV/K"] * 1e-6) ** 2 * e["bt2_sigma_S/m"] * T300 / (
            (KL[0] + KL[1]) / 2.0 + e["bt2_kappa_e_WF_W/mK"])
        self.assertAlmostEqual(e["bt2_ZT"], zt_avg, places=12)                   # 2D 平均仍只用面内 κL


class S83SummaryLinesTests(unittest.TestCase):
    """0024 的 [DPT 单谷] / [AMSET 机制] 两行要真的写进 comparison_summary.txt（以前只测了函数本身，
    本机合并 0030+V181 时调用被冲掉也没测出来）。"""

    def test_write_table_has_0024_lines(self):
        O = _g13()
        m = _material(tempfile.mkdtemp(), "2d")
        p = m / "step8.2_dpt" / "dpt_result.json"
        j = json.loads(p.read_text())
        j["temperatures_K"] = [300.0]
        j["results"][0]["single_valley_scope"] = {
            "edge_star": "K", "T0_K": 300.0, "share_T0": 0.4, "others": ["Q"],
            "shares_vs_T": [{"T_K": 300.0, "share": 0.4}], "note": "x"}
        p.write_text(json.dumps(j))
        run = m / "step8.4_amset2d"
        run.mkdir()
        diag = lambda v: [[v, 0, 0], [0, v, 0], [0, 0, 1]]
        (run / "transport.json").write_text(json.dumps({
            "doping": [-1e17, 1e17], "temperatures": [300.0],
            "mobility": {"overall": [[diag(30.0)]] * 2, "ADP": [[diag(400.0)]] * 2, "POP": [[diag(50.0)]] * 2}}))
        out = m / "step8.3_output"
        out.mkdir()
        with _cwd(m), contextlib.redirect_stdout(io.StringIO()):
            dpt = O.load_dpt(m)
            O.write_table(out, O.build_table(None, None, dpt), None, None, dpt)
        txt = (out / "comparison_summary.txt").read_text(encoding="utf-8")
        self.assertIn("[DPT 单谷] ", txt)
        self.assertIn("只代表 K 谷", txt)
        self.assertIn("[AMSET 机制] electron", txt)
        self.assertIn("限制因素 POP", txt)
        self.assertIn('_SKILL_REV = "2026-10-10-0026"',
                      (KE / "step8.3_output" / "gen_step13_output.py").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
