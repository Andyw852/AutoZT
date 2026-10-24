"""[0078] S8.3 mu_compare/：随化学势 u（带隙中点 = 0）的对比图 —— σ、S、PF、κe、μ 各两张
（BT2×DPT-τ + AMSET ADP only + 文献；再加 AMSET 各机制组合），列 = 方向。合成数据，纯本地。"""
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")
np = pytest.importorskip("numpy")

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
KB = 8.617333262e-5
C_A, CT = 20.0, 3.0
GAP, EF0 = 1.0, -1.6                       # AMSET 能量零点下本征能级在 −1.6 eV
TAU = {"electron": {"x": 1e-14, "y": 2e-14}, "hole": {"x": 3e-14, "y": 4e-14}}


def _s83():
    for d in (ROOT, KE, KE / "step8.4_amset2d", ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("s83_0078", str(KE / "step8.3_output" / "gen_step13_output.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _n_of_u(u, T):
    """带符号元胞浓度（负 = 电子）：抛物带 Boltzmann，N = 1e19。"""
    kt = KB * T
    return 1e19 * (math.exp((-GAP / 2 - u) / kt) - math.exp((u - GAP / 2) / kt))


def _material(tmp_path, lit=True):
    cwd = tmp_path / "MoX2" / "ke-dft-cpu"
    run = cwd / "step8.4_amset2d"
    run.mkdir(parents=True)
    (cwd / "step8.1_boltztrap").mkdir()
    temps = [300.0, 600.0]
    # BT2：S8.1 零点 = 带隙中点 + 0.2 eV（检验减中点）
    mu = [round(-2.0 + 0.01 * i, 3) for i in range(401)]
    ns = [[_n_of_u(m - 0.2, T) for m in mu] for T in temps]
    sot = [[abs(n) * 1e-3 + 1.0 for n in row] for row in ns]
    see = [[(-200e-6 if n < 0 else 200e-6) for n in row] for row in ns]
    bj = {"temperatures_K": temps, "mu_rel_efermi_eV": mu, "carrier_conc_cm-3": ns,
          "band_edges_rel_efermi_eV": {"vbm": 0.2 - GAP / 2, "cbm": 0.2 + GAP / 2},
          "sigma_over_tau_S_per_m_s": sot, "seebeck_V_per_K": see, "kappa_e_over_tau_W_per_m_K_s": sot,
          "sigma_over_tau_xx": sot, "sigma_over_tau_yy": [[2 * v for v in r] for r in sot],
          "seebeck_xx_V_per_K": see, "seebeck_yy_V_per_K": see,
          "kappa_e_over_tau_xx": [[v * 2.44e-8 * T for v in r] for r, T in zip(sot, temps)],
          "kappa_e_over_tau_yy": [[v * 2.44e-8 * T for v in r] for r, T in zip(sot, temps)],
          "lorenz_WOhm_per_K2": [[2.44e-8] * len(mu)] * 2}
    (cwd / "step8.1_boltztrap" / "boltztrap_crta.json").write_text(json.dumps(bj))
    # AMSET：掺杂档 = 出厂 [0077] 网格（三维等效）换回元胞口径；E_F 由同一抛物带反解
    dop3 = [s * 10.0 ** p for p in range(11, 22) for s in (-1, 1)]
    dop = [x / CT for x in dop3] + [0.0]

    def _ef(x, T):
        if not x:
            return EF0
        lo, hi = -GAP, GAP
        for _ in range(200):
            m = (lo + hi) / 2
            lo, hi = (m, hi) if _n_of_u(m, T) > x else (lo, m)
        return EF0 + (lo + hi) / 2

    def _t(f):
        return [[[[f(x, T), 0, 0], [0, 2 * f(x, T), 0], [0, 0, 0]] for T in temps] for x in dop]
    tj = {"doping": dop, "temperatures": temps, "fermi_levels": [[_ef(x, T) for T in temps] for x in dop],
          "conductivity": _t(lambda x, T: abs(x) * 1e-3 + 1.0), "seebeck": _t(lambda x, T: -150.0 if x < 0 else 150.0),
          "electronic_thermal_conductivity": _t(lambda x, T: 1e-3),
          "mobility": {"overall": _t(lambda x, T: 50.0), "ADP": _t(lambda x, T: 100.0), "IMP": _t(lambda x, T: 1e3)}}
    (run / "transport.json").write_text(json.dumps(tj))
    (run / "2d_correction.json").write_text(json.dumps({"cell_c_A": C_A, "elastic_rescale_factor_c_over_t": CT,
                                                        "num_atoms": 3}))
    (run / "undoped_ni.json").write_text(json.dumps({"temperatures_K": temps, "E_F_eV": [EF0, EF0],
                                                     "n_i_cm-3": [1.0, 1.0], "E_v_eV": [EF0 - GAP / 2] * 2,
                                                     "E_c_eV": [EF0 + GAP / 2] * 2}))
    if lit:
        d = tmp_path / "MoX2" / "literature"
        d.mkdir()
        (d / "paper.csv").write_text(
            "# digitized, u: mid-gap = 0\n"
            "quantity,axis,u_eV,value,T_K,label,material\n"
            "sigma,xx,-1.0,1e6,600,lit A,\n"
            "sigma,xx,-0.9,2e6,600,lit A,\n"
            "sigma,xx,-0.8,3e6,300,lit A,\n"            # 温度不对：丢
            "PF,yy,0.9,0.01,,lit A,MoX2\n"
            "PF,yy,0.8,0.02,,lit A,OtherMat\n"          # 别的材料：丢
            "seebeck,xx,0.5,-100,,,\n"                  # 别名 + 缺 label -> literature
            "foo,xx,0.5,1,,,\n")                        # 不认识的量：丢
    return cwd


@pytest.fixture
def S(monkeypatch):
    g = _s83()

    def _dpt(cwd, T=None, quiet=False):
        out = {}
        for c, m_ in (("electron", 0.5), ("hole", 0.8)):
            out[c], out["m_" + c] = 100.0, m_
            out["dir_" + c] = {"status": "ok", "x": {"tau_s": TAU[c]["x"]}, "y": {"tau_s": TAU[c]["y"]}}
        return out
    monkeypatch.setattr(g, "load_dpt", _dpt)
    return g


def test_literature_reader(tmp_path, S):
    cwd = _material(tmp_path)
    lit, src = S.load_literature(cwd, "MoX2", 600.0)
    assert sorted(lit) == [("PF", "yy"), ("S", "xx"), ("sigma", "xx")]
    assert lit[("sigma", "xx")] == {"lit A": [(-1.0, 1e6), (-0.9, 2e6)]}
    assert lit[("PF", "yy")] == {"lit A": [(0.9, 0.01)]}
    assert lit[("S", "xx")] == {"literature": [(0.5, -100.0)]}
    assert src == ["literature/paper.csv (4 pts)"]
    assert S.load_literature(cwd, "MoX2", 300.0)[0][("sigma", "xx")] == {"lit A": [(-0.8, 3e6)]}
    assert S.load_literature(tmp_path / "nowhere", "MoX2", 600.0) == ({}, [])


def test_data_frames_and_values(tmp_path, S):
    cwd = _material(tmp_path, lit=False)
    d = S.u_compare_data(cwd, True)
    assert d["T"] == 600.0 and d["dirs"] == ("xx", "yy") and d["fac"] == pytest.approx(CT) and d["nfac"] == pytest.approx(CT)
    bt = d["bt2"]
    assert bt["u"][0] == pytest.approx(-2.2) and bt["u"][-1] == pytest.approx(1.8)    # 减去 S8.1 零点下的带隙中点
    i = min(range(len(bt["u"])), key=lambda k: abs(bt["u"][k] - 0.7))               # n 型：电子 τ，分方向
    n = bt["n"][i]
    assert n < 0
    assert bt["vals"]["xx"]["sigma_S/m"][i] == pytest.approx((abs(n) * 1e-3 + 1) * TAU["electron"]["x"] * CT)
    assert bt["vals"]["yy"]["sigma_S/m"][i] == pytest.approx(2 * (abs(n) * 1e-3 + 1) * TAU["electron"]["y"] * CT)
    j = min(range(len(bt["u"])), key=lambda k: abs(bt["u"][k] + 0.7))               # p 型：空穴 τ
    assert bt["vals"]["xx"]["sigma_S/m"][j] == pytest.approx((abs(bt["n"][j]) * 1e-3 + 1) * TAU["hole"]["x"] * CT)
    k0 = min(range(len(bt["u"])), key=lambda k: abs(bt["u"][k]))
    assert bt["vals"]["xx"]["mu_cm2/Vs"][k0] is None                                # 带隙中央不给 μ
    # AMSET：u = E_F − AMSET 带隙中点；本征点不画；ADP only 只有 μ -> σ = |n| e μ × c/t（钉在本征能级的不给）
    assert set(d["amset"]) == {"overall", "ADP_only"}
    pts = {r["n"]: r for r in d["amset"]["ADP_only"]}
    r = pts[-1e21 / CT]
    assert r["u"] == pytest.approx(next(u for x, u in d["pairs"] if x == -1e21 / CT))
    assert r["u"] > 0.3 and not r["pinned"]
    assert r["vals"]["xx"]["sigma_S/m"] == pytest.approx(1e21 / CT * 1e6 * 1.602176634e-19 * 100.0 * 1e-4 * CT)
    assert r["vals"]["xx"]["S_uV/K"] is None
    low = pts[-1e11 / CT]
    assert low["pinned"] and low["vals"]["xx"]["sigma_S/m"] is None and low["vals"]["xx"]["mu_cm2/Vs"] == 100.0
    ov = {r["n"]: r for r in d["amset"]["overall"]}[1e18 / CT]
    assert ov["vals"]["yy"]["sigma_S/m"] == pytest.approx(2 * (1e18 / CT * 1e-3 + 1) * CT)
    assert ov["vals"]["xx"]["PF_W/mK2"] == pytest.approx((150e-6) ** 2 * (1e18 / CT * 1e-3 + 1) * CT)


def test_top_axis_ticks_sit_on_amset_points(tmp_path, S):
    cwd = _material(tmp_path, lit=False)
    d = S.u_compare_data(cwd, True)
    tk = dict((p, u) for u, p in S._ucmp_ticks(d, (-0.9, 0.9), min_frac=0.0))
    for p in (17, 19, 21):                                   # 掺杂网格内：刻度 = AMSET 该档的 E_F
        u_p = [u for x, u in d["pairs"] if abs(abs(x) * CT / 10.0 ** p - 1) < 1e-9]
        assert len(u_p) == 2 and min(abs(tk[p] - u) for u in u_p) < 1e-9
    full = S._ucmp_ticks(d, (-2.2, 1.8), min_frac=0.0)
    assert any(p == 22 for _, p in full)                     # 网格外由 BT2 的 n(u) 补
    thin = S._ucmp_ticks(d, (-2.2, 1.8))
    us = sorted(u for u, _ in thin)
    assert all(b - a >= 0.055 * 4.0 - 1e-9 for a, b in zip(us, us[1:]))   # 太挤的丢掉


def test_bt2_inversion(S):
    T = 600.0
    u = [round(-1.5 + 0.005 * i, 3) for i in range(601)]
    bt = {"u": u, "n": [_n_of_u(x, T) for x in u]}
    for x in (-1e18, 1e18, -1e20, 1e20):
        assert _n_of_u(S._bt2_u_of_n(bt, x), T) == pytest.approx(x, rel=0.02)
    assert S._bt2_u_of_n(bt, 1e30) is None


def test_write_outputs(tmp_path, S):
    cwd = _material(tmp_path)
    assert S.write_u_compare(cwd, cwd / "transport_plots", "MoX2", True)
    od = cwd / "transport_plots" / "mu_compare"
    for q in ("sigma", "S", "PF", "kappa_e", "mu"):
        assert (od / ("MoX2_%s_vs_u_600K.png" % q)).is_file()
        assert (od / ("MoX2_%s_vs_u_mechanisms_600K.png" % q)).is_file()
    rows = list(csv.DictReader([l for l in (od / "MoX2_vs_u_600K.csv").read_text().splitlines()
                                if not l.startswith("#")]))
    assert {r["source"] for r in rows} == {"BT2", "AMSET", "literature"}
    assert {r["dataset"] for r in rows if r["source"] == "AMSET"} == {"overall", "ADP_only"}
    lit = [r for r in rows if r["source"] == "literature"]
    assert len(lit) == 4 and {r["label"] for r in lit} == {"lit A", "literature"}
    am = [r for r in rows if r["source"] == "AMSET" and r["dataset"] == "overall" and r["axis"] == "xx"]
    assert sorted(float(r["n_axis_cm-3"]) for r in am if float(r["n_axis_cm-3"]) > 0) == pytest.approx(
        [10.0 ** p for p in range(11, 22)])


def test_wired_into_transport_plots():
    src = (KE / "step8.3_output" / "gen_step13_output.py").read_text(encoding="utf-8")
    i = src.index("write_mu_scan(cwd, out, plots, M, is_2d)\n    except")
    assert 0 < src.index("write_u_compare(cwd, plots, M, is_2d)", i) - i < 400
    import yaml
    steps = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    s83 = [s for grp in _walk(steps) for s in grp if isinstance(s, dict) and s.get("name") == "step8.3_output"]
    assert s83 and "literature" in (s83[0].get("push_paths") or [])


def _walk(o):
    if isinstance(o, list):
        yield o
        for x in o:
            yield from _walk(x)
    elif isinstance(o, dict):
        for x in o.values():
            yield from _walk(x)
