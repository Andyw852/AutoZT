"""[0023] 多谷 DPT：有贡献的谷之间形变势悬殊即拒绝出数；带符号 dE/dε 诊断；Jin et al. PRB 90, 045422 (2014) 文献对照。

起因：0022 实跑四个 TMD，电子 Q 谷 E1 只有 0.5–0.8 eV、K 谷 7.4–9.4 eV（真空口径），多谷 μ 被放大 7× 到 6×10⁵×
（WSe2 电子 639 357 cm²/Vs）；空穴（K+Γ）多/单 0.87–0.98，正常。
"""
import contextlib
import importlib.util
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


def _load(name, path):
    for d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(KE / "step8.2_dpt"), str(KE / "tools")):
        if d not in sys.path:
            sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T22 = _load("t22_for_0023", ROOT / "tests" / "test_0022_multi_valley_dpt.py")
E_C, HBAR, KB, M0, AMU = 1.602176634e-19, 1.054571817e-34, 1.380649e-23, 9.1093837015e-31, 1.66053906660e-27


def _D():
    import gen_step12_dpt as D
    return D


# ---------------- Jin et al. PRB 90, 045422 (2014)（LDA、无 SOC、QE；DFPT + 全能带蒙特卡洛） ----------------
#   Table I：a(Å)、E^c_QK、E^v_ΓK(meV)、m^c_K、m^c_Q（几何平均）、m^v_K、m^v_Γ（m0）
#   Table II 题注：v_LA（10⁵ cm/s）；Table IV/V：谷内声学 D1（eV，TA+LA 合并拟合）
#   Table III：μc（含全部散射）、μc,K（去掉 K–Q 谷间）、μv（cm²/Vs，300 K）
JIN = {
    "MoS2":  dict(a=3.14, dQK=81, dGK=148, mcK=0.51, mcQ=0.76, mvK=0.58, mvG=4.05, vLA=6.6,
                  DeK=4.5, DeQ=2.8, DhK=2.5, DhG=2.7, mu_c=130, mu_cK=320, mu_v=270, M="Mo", X="S"),
    "MoSe2": dict(a=3.27, dQK=28, dGK=374, mcK=0.64, mcQ=0.80, mvK=0.71, mvG=7.76, vLA=4.1,
                  DeK=3.4, DeQ=3.1, DhK=2.8, DhG=1.5, mu_c=25, mu_cK=180, mu_v=90, M="Mo", X="Se"),
    "WS2":   dict(a=3.10, dQK=67, dGK=173, mcK=0.31, mcQ=0.60, mvK=0.42, mvG=4.07, vLA=4.3,
                  DeK=3.2, DeQ=1.8, DhK=1.7, DhG=1.9, mu_c=320, mu_cK=690, mu_v=540, M="W", X="S"),
    "WSe2":  dict(a=3.25, dQK=16, dGK=427, mcK=0.39, mcQ=0.64, mvK=0.51, mvG=7.77, vLA=3.3,
                  DeK=3.2, DeQ=1.9, DhK=2.1, DhG=1.8, mu_c=30, mu_cK=250, mu_v=270, M="W", X="Se"),
}
AMASS = {"Mo": 95.95, "W": 183.84, "S": 32.06, "Se": 78.971}


def _C2D(p):
    """C = ρ v_LA²，ρ = 原胞质量 / 原胞面积（Jin 没给 ρ，按 Table I 的晶格常数与原子量算）。"""
    rho = (AMASS[p["M"]] + 2 * AMASS[p["X"]]) * AMU / (math.sqrt(3) / 2 * (p["a"] * 1e-10) ** 2)
    return rho * (p["vLA"] * 1e3) ** 2


def _star(label, g, m, D1, dE):
    """各向同性谷（Jin 的 Q 质量是纵横几何平均，文中说与电导质量很接近；D1 是标量）。"""
    W = [[1.0 / m, 0.0], [0.0, 1.0 / m]]
    return {"label": label, "members": [{"m_d": m, "dE": dE, "E1": [D1, D1], "W": W} for _ in range(g)]}


def _indep(vals, C, T=300.0):
    """独立实现（不调被测代码）：μ = (e/m0)·Σ g m e^{-ΔE/kT} τ/m / Σ g m e^{-ΔE/kT}，τ = ħ³C/(kT m D²)。"""
    kT = KB * T / E_C
    sn = sntw = 0.0
    for g, m, D, dE in vals:
        n = g * m * math.exp(-dE / kT)
        sn += n
        sntw += n * HBAR ** 3 * C / (KB * T * m * M0 * (D * E_C) ** 2) / m
    return E_C / M0 * sntw / sn * 1e4


class JinLiteratureTests(unittest.TestCase):
    """0022 的谷内公式 = Jin 式 (3)：1/τ = m D1² k_BT/(ħ³ ρ v_s²)（谷的态密度质量；Q 6 重、K 2 重）。
    用 Jin 自己的逐谷参数代进生产代码 _mv_mobility，核对：① 与独立实现一致；② 只 K ≥ Jin μc,K、多谷 ≥ Jin μc
    （DPT 只有谷内声学，必是上限）；③ Jin 的形变势（Q/K 比 0.6–0.9）不触发 0023 的闸门，实跑的真空口径（≈0.07）触发。"""

    def test_production_weighting_matches_jin_formula_and_bounds(self):
        D = _D()
        expect = {"MoS2": (548, 563, 1358), "MoSe2": (344, 300, 412), "WS2": (1980, 1886, 3799),
                  "WSe2": (924, 960, 1255)}
        for name, p in JIN.items():
            C = _C2D(p)
            eK = [_star("K", 2, p["mcK"], p["DeK"], 0.0)]
            emv = eK + [_star("Q", 6, p["mcQ"], p["DeQ"], p["dQK"] / 1000.0)]
            hmv = [_star("K", 2, p["mvK"], p["DhK"], 0.0), _star("Γ", 1, p["mvG"], p["DhG"], p["dGK"] / 1000.0)]
            mu_eK = D._mv_mobility(eK, (C, C), 300.0, True)["mu"][0]
            mu_e = D._mv_mobility(emv, (C, C), 300.0, True)["mu"][0]
            mu_h = D._mv_mobility(hmv, (C, C), 300.0, True)["mu"][0]
            self.assertAlmostEqual(mu_eK / _indep([(2, p["mcK"], p["DeK"], 0.0)], C), 1.0, places=9)
            self.assertAlmostEqual(mu_e / _indep([(2, p["mcK"], p["DeK"], 0.0),
                                                   (6, p["mcQ"], p["DeQ"], p["dQK"] / 1000.0)], C), 1.0, places=9)
            for got, want in zip((mu_eK, mu_e, mu_h), expect[name]):
                self.assertAlmostEqual(got, want, delta=1.0, msg=name)
            self.assertGreaterEqual(mu_eK, p["mu_cK"], name)        # 只 K 谷内声学 ≥ Jin（+光学、K↔K′）
            self.assertGreaterEqual(mu_e, p["mu_c"], name)          # 多谷无谷间 ≥ Jin 全部散射
            self.assertGreaterEqual(mu_h, p["mu_v"], name)
            temps = [100.0, 300.0, 900.0]
            self.assertIsNone(D._mv_e1_gate(emv, (C, C), temps, True)[0], name)
            self.assertIsNone(D._mv_e1_gate(hmv, (C, C), temps, True)[0], name)

    def test_wse2_intervalley_gap_is_large(self):
        """Jin：WSe2 电子含谷间 30、只 K 250 —— 多谷 DPT（无谷间）只能当上限，这里差 30 倍。"""
        D = _D()
        p = JIN["WSe2"]
        C = _C2D(p)
        mu = D._mv_mobility([_star("K", 2, p["mcK"], p["DeK"], 0.0),
                             _star("Q", 6, p["mcQ"], p["DeQ"], p["dQK"] / 1000.0)], (C, C), 300.0, True)
        self.assertGreater(mu["mu"][0] / p["mu_c"], 20.0)
        self.assertGreater(mu["w"]["Q"], 0.7)                     # ΔE_QK = 16 meV：Q 谷占七成以上

    def test_vacuum_e1_ratio_from_real_run_trips_gate(self):
        """实跑的真空口径 E1（WSe2：Q 0.72 eV、K 8.64 eV）代进 Jin 的质量/能量差 -> 闸门拒绝。"""
        D = _D()
        p = JIN["WSe2"]
        C = _C2D(p)
        stars = [_star("K", 2, p["mcK"], 8.64, 0.0), _star("Q", 6, p["mcQ"], 0.72, p["dQK"] / 1000.0)]
        msg, info = D._mv_e1_gate(stars, (C, C), [100.0, 300.0, 900.0], True)
        self.assertIsNotNone(msg)
        self.assertTrue(msg.startswith("拒绝"))
        self.assertIn("Q E1=0.72", msg)
        self.assertAlmostEqual(info["K"]["E1_eff_eV"], 8.64, places=3)


class GateEndToEndTests(unittest.TestCase):
    """合成六方单层（0022 的模型），把 Q 谷形变势改成实跑那样小：多谷拒绝出数、只留诊断；空穴不受影响。"""

    def setUp(self):
        self.D = _D()
        self.D._MV_CACHE.clear()
        self.m = Path(tempfile.mkdtemp())
        self._p = [mock.patch.dict(self.D.MANUAL_ANISO, {"C_2D_xy_N_per_m": (T22.C2D, T22.C2D)}),
                   mock.patch.object(self.D, "VALLEY_WINDOW_EV", 0.30),
                   mock.patch.object(self.D, "DPT_VALLEYS", "auto"),
                   mock.patch.object(self.D, "E1_SOURCE", "vac")]
        for q in self._p:
            q.start()

    def tearDown(self):
        for q in self._p:
            q.stop()

    @staticmethod
    def _model(eK_cb=0.03, DK=8.6, DQ=(0.7, 0.5)):
        import numpy as np
        model = T22.Model(eK_cb=eK_cb)
        for v in model.cb:
            if v["lab"] == "K":
                v["D"] = DK * np.eye(2)
            else:
                t = math.atan2(*((v["k"] @ model.recip)[:2][::-1]))
                R = T22._rot(t)
                v["D"] = R @ np.diag(DQ) @ R.T
        return model

    def _run(self, model, carrier, **kw):
        vr = T22.build(model, self.m)
        with mock.patch("pymatgen.io.vasp.Vasprun", vr), contextlib.redirect_stdout(io.StringIO()) as buf:
            bd = self.D._aniso_block(self.m, True, carrier, 300.0)
            mv = self.D._multi_valley_block(self.m, True, carrier, 300.0, [100.0, 300.0, 600.0], bd)
        return mv, buf.getvalue()

    def test_wse2_like_q_edge_refused_with_signed_diagnostics(self):
        mv, out = self._run(self._model(), "electron")
        self.assertTrue(mv["status"].startswith("拒绝"), mv["status"])
        for k in ("x", "y", "mobility_vs_T", "mobility_inplane_cm2_Vs"):
            self.assertNotIn(k, mv)                              # 下游拿不到任何 μ/τ
        labs = {v["label"]: v for v in mv["valleys"]}
        self.assertEqual(set(labs), {"Q", "K"})
        self.assertAlmostEqual(labs["K"]["E1_eff_eV"], 8.6, places=2)
        self.assertLess(labs["Q"]["E1_eff_eV"], 0.8)
        for v in mv["valleys"]:                                  # raw − 真空项 = 带符号的真空口径 E1
            for (rx, ry), (sx, sy) in zip(v["dEdeps_raw_xy_eV"], v["E1_signed_xy_eV"]):
                self.assertAlmostEqual(rx - sx, T22.DVAC["xx"], places=3)
                self.assertAlmostEqual(ry - sy, T22.DVAC["yy"], places=3)
        self.assertEqual(mv["dvac_eV_per_unit_strain"], T22.DVAC)
        self.assertIn("拒绝", out)

    def test_mos2_like_small_q_population_still_refused(self):
        """K 带边、Q 高 80 meV（300 K 布居只有几个百分点），Q 的 τ 长两个数量级照样主导 -> 拒绝。"""
        mv, _ = self._run(self._model(eK_cb=-0.08), "electron")
        self.assertTrue(mv["status"].startswith("拒绝"), mv["status"])
        q = next(v for v in mv["valleys"] if v["label"] == "Q")
        self.assertLess(q["weight_T0"], 0.25)                    # 300 K 布居不到四分之一
        self.assertGreaterEqual(q["pop_max"], self.D.VALLEY_W_MIN)

    def test_holes_unaffected(self):
        mv, _ = self._run(self._model(), "hole")
        self.assertEqual(mv["status"], "ok", mv["status"])
        self.assertIn("E1_eff_eV", mv["valleys"][0])

    def test_ratio_override(self):
        with mock.patch.object(self.D, "VALLEY_E1_RATIO_MAX", 100.0):
            mv, _ = self._run(self._model(), "electron")
        self.assertEqual(mv["status"], "ok", mv["status"])

    def test_s81_gets_no_electron_tau_when_refused(self):
        B = _load("bt2_0023", KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py")
        mv, _ = self._run(self._model(), "electron")
        d = self.m / "step8.2_dpt"
        d.mkdir(exist_ok=True)
        (d / "dpt_result.json").write_text(json.dumps({
            "is_2d": True, "temperature_K": 300.0,
            "results": [{"carrier": "electron", "mobility_cm2_Vs": None, "inputs": {"m_eff_m0": None},
                         "by_direction": {"status": "非高对称谷"}, "multi_valley": mv}]}, ensure_ascii=False))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertNotIn("electron", B._tau_from_json(self.m, 300.0, True))


if __name__ == "__main__":
    unittest.main()
