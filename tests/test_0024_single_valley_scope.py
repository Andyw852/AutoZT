"""[0024] 单谷 DPT 只代表带边那组谷：S8.2 写明它在各温度占多少载流子（多谷出没出数都写），S8.1/S8.3/mobility_vs_dpt 显示；
S8.3 加 AMSET 各散射机制的迁移率与限制因素。只报布居、不设阈值、不拦截。

起因：0023 实跑后 MoS2/WS2/MoSe2 的电子单谷值仍标 ok、S8.1 照用、S8.3 当主栏，但它只是 K 谷的值 ——
MoSe2 300 K 时六成电子在 Q 谷（Jin PRB 90, 045422：MoSe2 电子含谷间 25、只 K 180 cm²/Vs）。
"""
import contextlib
import importlib.util
import io
import json
import math
import os
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


T22 = _load("t22_for_0024", ROOT / "tests" / "test_0022_multi_valley_dpt.py")
T23 = _load("t23_for_0024", ROOT / "tests" / "test_0023_mv_e1_gate.py")
KB_EV = 1.380649e-23 / 1.602176634e-19


def _D():
    import gen_step12_dpt as D
    return D


def _share_K(dQ, T, mK=0.40, mQ=math.sqrt(0.8 * 0.5)):
    """独立算：K 星（2 个谷）在 T 时占的载流子比例，Q 星 6 个谷高 dQ。"""
    nK = 2 * mK
    nQ = 6 * mQ * math.exp(-dQ / (KB_EV * T))
    return nK / (nK + nQ)


class S82ScopeTests(unittest.TestCase):
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

    def _mv(self, model, carrier, temps=(100.0, 300.0, 600.0)):
        vr = T22.build(model, self.m)
        with mock.patch("pymatgen.io.vasp.Vasprun", vr), contextlib.redirect_stdout(io.StringIO()):
            bd = self.D._aniso_block(self.m, True, carrier, 300.0)
            return self.D._multi_valley_block(self.m, True, carrier, 300.0, list(temps), bd)

    def test_refused_multi_valley_still_reports_scope(self):
        """MoS2 型：K 带边、Q 高 80 meV 且形变势小 -> 多谷拒绝，但单谷值只代表 K 谷、占多少照样写出。"""
        mv = self._mv(T23.GateEndToEndTests._model(eK_cb=-0.08), "electron")
        self.assertTrue(mv["status"].startswith("拒绝"))
        self.assertEqual(mv["edge_star"], "K")
        sc = self.D._single_valley_scope(mv, 300.0)
        self.assertEqual(sc["edge_star"], "K")
        self.assertAlmostEqual(sc["share"], _share_K(0.08, 300.0), places=3)
        for row in sc["shares_vs_T"]:
            self.assertAlmostEqual(row["share"], _share_K(0.08, row["T_K"]), places=3)
        self.assertIn("单谷值只代表 K 谷", sc["note"])
        self.assertIn("600 K 时降到", sc["note"])                 # 高温时 K 占比最低
        self.assertIn("Q", sc["others"])

    def test_ok_multi_valley_and_q_edge(self):
        mv = self._mv(T22.Model(eK_cb=-0.03), "electron")        # K 带边、Q 高 30 meV、形变势正常 -> 多谷 ok
        self.assertEqual(mv["status"], "ok", mv["status"])
        sc = self.D._single_valley_scope(mv, 300.0)
        self.assertAlmostEqual(sc["share"], _share_K(0.03, 300.0), places=3)
        self.D._MV_CACHE.clear()
        self.m = Path(tempfile.mkdtemp())
        mv = self._mv(T22.Model(), "electron")                    # WSe2 型：带边在 Q
        self.assertEqual(self.D._single_valley_scope(mv, 300.0)["edge_star"], "Q")

    def test_single_star_has_no_scope(self):
        with mock.patch.object(self.D, "VALLEY_WINDOW_EV", 0.05):
            mv = self._mv(T22.Model(), "hole")                    # 价带只剩 K
        self.assertEqual([v["label"] for v in mv["valleys"]], ["K"])
        self.assertIsNone(self.D._single_valley_scope(mv, 300.0))
        self.assertIsNone(self.D._single_valley_scope({"status": "off（step.conf DPT_VALLEYS = off）"}, 300.0))


class S82MainTests(unittest.TestCase):
    def test_main_writes_scope_and_summary_line(self):
        import numpy as np
        D = _D()
        D._MV_CACHE.clear()
        m = Path(tempfile.mkdtemp())
        model = T23.GateEndToEndTests._model(eK_cb=-0.08)
        vr = T22.build(model, m)
        (m / "step1_opt").mkdir()
        (m / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\n")
        model.st.to(filename=str(m / D.UNIFORM_DIR / "POSCAR"), fmt="poscar")
        (m / D.UNIFORM_DIR / "KPOINTS").write_text("auto\n0\nGamma\n48 48 1\n0 0 0\n")
        L = model.st.lattice.matrix
        h = abs(np.dot(L[2], np.cross(L[0], L[1]))) / np.linalg.norm(np.cross(L[0], L[1]))
        ck = T22.C2D / (h * 1e-10) / 1e9 * 10.0
        rows = "\n".join(" XX %s" % " ".join("%.4f" % (ck if i == j and i < 2 else 1.0) for j in range(6))
                         for i in range(6))
        (m / D.ELASTIC_DIR).mkdir()
        (m / D.ELASTIC_DIR / "OUTCAR").write_text(" TOTAL ELASTIC MODULI (kBar)\n Direction XX YY\n" + rows + "\n")
        cwd = os.getcwd()
        os.chdir(m)
        try:
            with mock.patch("pymatgen.io.vasp.Vasprun", vr), mock.patch.object(D, "TEMPERATURES", "300,600"), \
                    mock.patch.object(D, "VALLEY_WINDOW_EV", 0.30), mock.patch.object(D, "DPT_VALLEYS", "auto"), \
                    contextlib.redirect_stdout(io.StringIO()) as buf:
                D.main()
        finally:
            os.chdir(cwd)
        j = json.loads((m / D.OUTDIR_NAME / "dpt_result.json").read_text())
        e = next(r for r in j["results"] if r["carrier"] == "electron")
        self.assertEqual(e["by_direction"]["status"], "ok")          # K 带边：单谷照常出数
        self.assertTrue(e["multi_valley"]["status"].startswith("拒绝"))
        self.assertAlmostEqual(e["single_valley_scope"]["share"], _share_K(0.08, 300.0), places=3)
        txt = (m / D.OUTDIR_NAME / "dpt_summary.txt").read_text()
        self.assertIn("[单谷代表] 单谷值只代表 K 谷", txt)
        self.assertIn("单谷值只代表 K 谷", buf.getvalue())


def _dpt_json(m, share300=0.40, share600=0.30):
    d = m / "step8.2_dpt"
    d.mkdir(parents=True, exist_ok=True)
    bd = {"status": "ok", "m_d_m0": 0.6,
          "x": {"mobility_cm2_Vs": 160.0, "tau_s": 5.0e-14}, "y": {"mobility_cm2_Vs": 160.0, "tau_s": 5.0e-14}}
    scope = {"edge_star": "K", "T_K": 300.0, "share": share300, "others": {"Q": 1 - share300},
             "shares_vs_T": [{"T_K": 300.0, "share": share300}, {"T_K": 600.0, "share": share600}],
             "note": "单谷值只代表 K 谷"}
    res = {"dim": "2d", "is_2d": True, "temperature_K": 300.0, "temperatures_K": [300.0, 600.0],
           "results": [{"carrier": "electron", "mobility_cm2_Vs": 150.0, "inputs": {"m_eff_m0": 0.6},
                        "by_direction": bd, "single_valley_scope": scope,
                        "multi_valley": {"status": "拒绝：……", "edge_star": "K"}}]}
    (d / "dpt_result.json").write_text(json.dumps(res, ensure_ascii=False))


class S81Tests(unittest.TestCase):
    def test_tau_note(self):
        B = _load("bt2_0024", KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py")
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            t = B._tau_from_json(m, 600.0, True)
        self.assertAlmostEqual(t["electron"]["x"], 2.5e-14, delta=1e-21)   # τ 不变（单谷分方向 × T0/T）
        self.assertIn("只代表 K 谷：600 K 时占载流子 30.0%", B._TAU_SCOPE["electron"])
        self.assertIn("[WARN] electron", buf.getvalue())               # 不到一半 -> WARN


class S83Tests(unittest.TestCase):
    def _mod(self):
        return _load("out_0024", KE / "step8.3_output" / "gen_step13_output.py")

    def test_scope_lines_and_mu_T_column(self):
        O = self._mod()
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        with contextlib.redirect_stdout(io.StringIO()):
            d = O.load_dpt(m)
        lines = O._dpt_scope_lines(d)
        self.assertEqual(len(lines), 1)
        self.assertIn("★ electron：dpt_mu 只代表 K 谷（300 K 占载流子 40.0%；Q 未计入", lines[0])
        res = O.build_mu_vs_T(m, True)
        r600 = next(r for r in res["rows"] if r["T_K"] == 600.0 and r["direction"] == "x")  # [0030] 分方向
        self.assertAlmostEqual(r600["DPT_edge_share"], 0.30)
        self.assertEqual(res["meta"]["dpt_scope"], {"electron": "K 谷"})

    def test_amset_mechanism_lines(self):
        O = self._mod()
        m = Path(tempfile.mkdtemp())
        run = m / "step8.4_amset2d"
        run.mkdir()
        diag = lambda v: [[v, 0, 0], [0, v, 0], [0, 0, 1]]
        temps = [300.0, 600.0]
        dop = [-1e19, -1e17, 1e17]
        mob = {"overall": [[diag(30.0)] * 2] * 3,
               "ADP": [[diag(4640.0)] * 2] * 3,
               "POP": [[diag(335.0)] * 2] * 3,
               "IMP": [[diag(86.0)] * 2, [diag(86.4)] * 2, [diag(900.0)] * 2]}
        (run / "transport.json").write_text(json.dumps({"doping": dop, "temperatures": temps, "mobility": mob}))
        lines = O._amset_mech_lines(m, True, 300.0)
        e = next(x for x in lines if "electron" in x)
        self.assertIn("掺杂 -1.00e+17", e)                         # |掺杂| 最低一档
        self.assertIn("IMP 86.4 < POP 335 < ADP 4640 → 限制因素 IMP；本征（不计 IMP）限制因素 POP", e)
        h = next(x for x in lines if "hole" in x)
        self.assertIn("限制因素 POP", h)
        self.assertNotIn("本征", h)                                  # 本征限制因素与总的相同，不重复

    def test_multi_valley_main_column_drops_scope(self):
        """主栏已改用多谷（单谷算不出）时，单谷的适用范围不再显示。"""
        O = self._mod()
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        p = m / "step8.2_dpt" / "dpt_result.json"
        j = json.loads(p.read_text())
        r = j["results"][0]
        r["by_direction"], r["mobility_cm2_Vs"] = {"status": "非高对称谷"}, None
        r["multi_valley"] = {"status": "ok", "x": {"m_c_m0": 0.6}, "y": {"m_c_m0": 0.6}, "valleys": [],
                             "mobility_vs_T": [{"T_K": 300.0, "mobility_x_cm2_Vs": 400.0, "mobility_y_cm2_Vs": 400.0,
                                                "mobility_inplane_cm2_Vs": 400.0, "tau_x_s": 1e-13, "tau_y_s": 1e-13}]}
        p.write_text(json.dumps(j, ensure_ascii=False))
        with contextlib.redirect_stdout(io.StringIO()):
            d = O.load_dpt(m)
        self.assertEqual(d["src_electron"], "multi_valley")
        self.assertNotIn("scope_electron", d)
        self.assertEqual(O._dpt_scope_lines(d), [])


class ToolTests(unittest.TestCase):
    def test_scope_line_printed(self):
        import mobility_vs_dpt as T
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        scopes = T.dpt_scope(m / "step8.2_dpt" / "dpt_result.json")
        self.assertIn("只代表 K 谷：600 K 时占载流子 30.0%", T.scope_line("electron", scopes["electron"], 600.0))
        diag = lambda v: [[v, 0, 0], [0, v, 0], [0, 0, 1]]
        tj = {"doping": [-1e17], "temperatures": [300.0], "seebeck": [[diag(-300.0)]],
              "mobility": {"overall": [[diag(30.0)]], "ADP": [[diag(800.0)]]}}
        rs, mechs = T.rows(tj, 300.0, True, {"electron": 160.0, "hole": None}, None, None, dpt_T0=300.0)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            T._print_block(Path("t.json"), True, rs, mechs, {"electron": 160.0}, {}, 300.0, set(), set(), 300.0,
                           scopes=scopes)
        self.assertIn("★ DPT(电子(n)) 只代表 K 谷：300 K 时占载流子 40.0%", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
