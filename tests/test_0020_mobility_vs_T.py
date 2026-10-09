"""[0020] 迁移率多温度：S8.2 的 mobility_vs_T（2D ∝1/T、3D ∝T^-1.5）、S8.3 的 μ(T) 三值表
（DPT / AMSET 只 ADP / AMSET 本征 / AMSET 全机制）、tools/mobility_vs_dpt.py 的 --T all 与 DPT 换温。纯本地。"""
import contextlib
import csv
import importlib.util
import io
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _load(name, path):
    for d in (ROOT, KE, ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _dpt_mod():
    return _load("dpt_0020", KE / "step8.2_dpt" / "gen_step12_dpt.py")


def _diag(v, z=1.0):
    return [[v, 0, 0], [0, v * 1.2, 0], [0, 0, z]]


TEMPS = [300.0, 600.0, 900.0]


def _material(tmp_path, intrinsic=True, T0=300.0):
    run = tmp_path / "step8.4_amset2d"
    run.mkdir()
    dop = [-1e17, -1e19, 1e17, 1e19]

    def grid(base):          # base[掺杂] × 温度 -> μ ∝ 1/T
        return [[_diag(b * 300.0 / T) for T in TEMPS] for b in base]
    tj = {"doping": dop, "temperatures": TEMPS,
          "seebeck": [[_diag(-300.0 if x < 0 else 300.0) for _ in TEMPS] for x in dop],
          "mobility": {"overall": grid([30.0, 20.0, 90.0, 60.0]),
                       "ADP": grid([61.0, 50.0, 201.0, 150.0])}}
    (run / "transport.json").write_text(json.dumps(tj))
    (run / "2d_correction.json").write_text(json.dumps({"cell_c_A": 20.0, "areal_density_factor_cm": 2e-7}))
    if intrinsic:
        for name, base, kept in (("intrinsic_ADP.json", [60.0, 50.0, 200.0, 150.0], ["ADP"]),
                                 ("intrinsic_transport.json", [45.0, 35.0, 120.0, 100.0], ["ADP", "POP"])):
            (run / name).write_text(json.dumps({
                "doping_cm3": dop, "temperatures_K": TEMPS, "intrinsic_labels": kept,
                "intrinsic_dropped": ["IMP"] if kept == ["ADP", "POP"] else ["IMP", "POP", "PIE"],
                "intrinsic": {"mobility_cm2_Vs_s": {"overall": grid(base)}}}))
    res = []
    for car, head, bx, by in (("electron", 150.0, 100.0, 120.0), ("hole", 500.0, 300.0, 340.0)):
        res.append({"carrier": car, "mobility_cm2_Vs": head, "inputs": {"m_eff_m0": 0.9},
                    "by_direction": {"status": "ok", "m_d_m0": 0.95,
                                     "x": {"mobility_cm2_Vs": bx, "tau_s": 1e-14},
                                     "y": {"mobility_cm2_Vs": by, "tau_s": 2e-14}}})
    (tmp_path / "step8.2_dpt").mkdir()
    (tmp_path / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(
        {"temperature_K": T0, "is_2d": True, "results": res}))
    return tmp_path


# ---------------------------------------------------------------- S8.2
def test_parse_temps_range_and_list():
    D = _dpt_mod()
    assert D.parse_temps("100:900:9") == [100.0 * i for i in range(1, 10)]
    assert D.parse_temps("100:700:7") == [100.0 * i for i in range(1, 8)]
    assert D.parse_temps("600, 300 ,300") == [300.0, 600.0]
    for bad in ("100:900", "0:300:3", "abc"):
        with pytest.raises(ValueError):
            D.parse_temps(bad)


def test_mobility_vs_T_scaling_2d_and_3d():
    D = _dpt_mod()
    rec = {"mobility_cm2_Vs": 150.0,
           "by_direction": {"status": "ok", "x": {"mobility_cm2_Vs": 100.0, "tau_s": 1e-14},
                            "y": {"mobility_cm2_Vs": 120.0, "tau_s": 2e-14}}}
    rows = {r["T_K"]: r for r in D.mobility_vs_T(rec, True, 300.0, [300.0, 600.0])}
    assert rows[300.0]["mobility_cm2_Vs"] == 150.0
    assert rows[600.0]["mobility_cm2_Vs"] == 75.0
    assert rows[600.0]["mobility_x_cm2_Vs"] == 50.0 and rows[600.0]["mobility_y_cm2_Vs"] == 60.0
    assert rows[600.0]["mobility_inplane_cm2_Vs"] == 55.0
    assert rows[600.0]["tau_x_s"] == pytest.approx(0.5e-14)
    r3 = D.mobility_vs_T({"mobility_cm2_Vs": 100.0, "by_direction": {"status": "3D"}}, False, 300.0, [1200.0])[0]
    assert r3["mobility_cm2_Vs"] == pytest.approx(100.0 / 8.0)
    assert "mobility_x_cm2_Vs" not in r3


def test_step_conf_temperatures(tmp_path):
    D = _dpt_mod()
    (tmp_path / "step.conf").write_text("[params]\nSTEP = step8.2_dpt\nTEMPERATURES = 100:700:7\n"
                                        "TEMPERATURE_K = 600\n")
    rec = D.apply_conf(tmp_path)
    assert D.TEMPERATURES == "100:700:7" and D.TEMPERATURE_K == 600.0
    assert rec["TEMPERATURES"]["source"].startswith("step.conf")
    (tmp_path / "step.conf").write_text("[params]\nSTEP = step8.2_dpt\nTEMPERATURES = 100:700\n")
    with pytest.raises(SystemExit):
        D.apply_conf(tmp_path)


# ---------------------------------------------------------------- S8.3
def _s83():
    return _load("s83_0020", KE / "step8.3_output" / "gen_step13_output.py")


def test_s83_mu_vs_T_three_values(tmp_path):
    S = _s83()
    mat = _material(tmp_path)
    res = S.build_mu_vs_T(mat, True)
    rows = {(r["T_K"], r["carrier"], r["direction"]): r for r in res["rows"]}
    assert sorted({r["T_K"] for r in res["rows"]}) == TEMPS
    e600 = rows[(600.0, "electron", "x")]                        # [0030] 按方向，不出平均
    assert e600["doping_cm3"] == -1e17                           # 最低掺杂
    assert e600["n2D_cm2"] == pytest.approx(1e17 * 2e-7)
    assert e600["mu_DPT"] == pytest.approx(100.0 / 2)            # 100 × 300/600
    assert e600["mu_AMSET_ADP"] == pytest.approx(60.0 / 2)       # intrinsic_ADP.json，不是 transport 的 61
    assert e600["mu_AMSET_intrinsic"] == pytest.approx(45.0 / 2)
    assert e600["mu_AMSET_all"] == pytest.approx(30.0 / 2)
    assert e600["m_d_m0"] == 0.95 and e600["tau_DPT_fs"] == pytest.approx(10.0)
    assert e600["E1_eV"] is None and e600["C_2D_N_per_m"] is None  # 夹具 by_direction 未带这些字段
    assert "ADP_over_DPT" not in e600                            # [0030] 不再输出比值
    assert rows[(600.0, "electron", "y")]["mu_DPT"] == pytest.approx(120.0 / 2)
    h300x = rows[(300.0, "hole", "x")]
    assert h300x["mu_DPT"] == 300.0 and h300x["mu_AMSET_ADP"] == 200.0 and h300x["doping_cm3"] == 1e17
    assert "intrinsic_ADP.json" in res["meta"]["amset_sources"]["mu_AMSET_ADP"]
    assert res["meta"]["intrinsic_labels"] == ["ADP", "POP"]


def test_s83_falls_back_to_transport_adp(tmp_path):
    S = _s83()
    res = S.build_mu_vs_T(_material(tmp_path, intrinsic=False), True)
    r = [x for x in res["rows"] if x["T_K"] == 300.0 and x["carrier"] == "electron" and x["direction"] == "x"][0]
    assert r["mu_AMSET_ADP"] == 61.0 and r["mu_AMSET_intrinsic"] is None
    assert "transport.json mobility.ADP" in res["meta"]["amset_sources"]["mu_AMSET_ADP"]


def test_s83_write_outputs(tmp_path):
    S = _s83()
    mat = _material(tmp_path)
    out = mat / "step8.3_output"
    out.mkdir()
    (out / "comparison_summary.txt").write_text("# 300 K\n")
    with contextlib.redirect_stdout(io.StringIO()):
        S.write_mu_vs_T(out, mat, True)
    lines = [ln for ln in (out / S.MU_T_CSV).read_text().splitlines() if not ln.startswith("#")]
    rs = list(csv.DictReader(lines))
    assert len(rs) == len(TEMPS) * 2 * 2                         # 温度 × 载流子 × (x, y)，不出平均
    assert "ADP_over_DPT" not in rs[0] and "E1_eV" in rs[0]      # [0030] 分方向参数列在、比值列没了
    assert json.loads((out / S.MU_T_JSON).read_text())["meta"]["dpt_scaling"] == "1/T"
    summ = (out / "comparison_summary.txt").read_text()
    assert summ.startswith("# 300 K") and "迁移率随温度" in summ
    assert (out / S.MU_T_PNG).is_file()


# ---------------------------------------------------------------- tools/mobility_vs_dpt.py
def _tool():
    return _load("mvd_0020", KE / "tools" / "mobility_vs_dpt.py")


def test_tool_scales_dpt_to_requested_T(tmp_path):
    M = _tool()
    mat = _material(tmp_path)
    out = tmp_path / "o.json"
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        rc = M.main([str(mat), "--T", "600", "--json", str(out)])
    assert rc == 0 and "由 300 K 按 1/T 换算" in buf.getvalue()
    rows = {(r["carrier"], r["doping"]): r for r in json.loads(out.read_text())}
    e = rows[("electron", -1e17)]
    assert e["T"] == 600.0 and e["DPT"] == 55.0
    assert e["ADP_over_DPT"] == round(e["ADP"] / 55.0, 3)


def test_tool_T_all(tmp_path):
    M = _tool()
    mat = _material(tmp_path)
    out = tmp_path / "o.json"
    with contextlib.redirect_stdout(io.StringIO()):
        assert M.main([str(mat), "--T", "all", "--json", str(out)]) == 0
    rs = json.loads(out.read_text())
    assert sorted({r["T"] for r in rs}) == TEMPS and len(rs) == 4 * len(TEMPS)
    dpt = {r["T"]: r["DPT"] for r in rs if r["carrier"] == "hole"}
    assert dpt[300.0] == 320.0 and dpt[900.0] == pytest.approx(320.0 / 3, abs=1e-3)
    with contextlib.redirect_stderr(io.StringIO()):
        assert M.main([str(mat), "--T", "hot"]) == 2
