"""[0056] S8.3 transport / phonon-limited 掺杂扫描按浓度拆两段，各一个文件夹（图 + 每温度表）：
<前缀>_carrier_tuning_1e17-1e21/（载流子调控区）与 <前缀>_low_doping_<NLO>-1e16/（低浓度段，有 DOPING_LOW 才有）。
合成 intrinsic_transport_table.csv + transport.json 掺杂网格，纯本地。"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"

COLS = ["T_K", "dataset", "type", "n_cell_cm-3", "n_2D_cm-2", "n_lit_cm-3", "carriers_per_cell",
        "carriers_per_atom", "in_interpret_range", "S_uV/K", "sigma_S/m", "PF_W/mK2", "kappa_e_W/mK",
        "mu_cm2/Vs", "flag"]


def _s83():
    for d in (ROOT, KE, KE / "step8.4_amset2d", ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("s83_0056", str(KE / "step8.3_output" / "gen_step13_output.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

LOW = [1e10, 1e12, 1e14, 1e16]
HIGH = [1e17, 1e18, 1e19, 1e20, 1e21]


def _material(tmp_path, dop):
    cwd = tmp_path / "MoS2" / "ke-dft-cpu"
    out = cwd / "step8.3_output"
    (out / "intrinsic").mkdir(parents=True)
    (cwd / "step8.4_amset2d").mkdir()
    grid = sorted([-x for x in dop] + [0.0] + dop)
    (cwd / "step8.4_amset2d" / "transport.json").write_text(json.dumps({"doping": grid, "temperatures": [300.0]}))
    with open(out / "intrinsic" / "intrinsic_transport_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for T in (300.0, 600.0):
            for ds in ("overall", "intrinsic", "ADP_only", "BT2_DPTtau", "ADP_POP"):
                for n in grid:
                    w.writerow([T, ds, "n" if n < 0 else "p", n, "", "", "", "", "yes",
                                -200 if n < 0 else 200, 1e4 * (abs(n) + 1) / 1e18, 1e-3, 0.1, 100, ""])
    return cwd, out


def _ns(p):
    rows = list(csv.DictReader(p.read_text().splitlines()))
    assert {r["dataset"] for r in rows} <= {"overall", "intrinsic", "ADP_only", "BT2_DPTtau"}
    return sorted({abs(float(r["n_cell_cm-3"])) for r in rows})


def test_two_folders_with_low_doping(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, LOW + HIGH)
    g.write_transport_plots(cwd, out, True)
    tp = cwd.parent / "transport_plots"
    hi = tp / "phonon_limited_carrier_tuning_1e17-1e21"
    lo = tp / "phonon_limited_low_doping_1e10-1e16"
    assert (hi / "MoS2_phonon_limited_carrier_tuning.png").is_file()
    assert (lo / "MoS2_phonon_limited_low_doping.png").is_file()
    for T in (300, 600):
        assert _ns(hi / ("MoS2_phonon_limited_carrier_tuning_%dK.csv" % T)) == HIGH
        assert _ns(lo / ("MoS2_phonon_limited_low_doping_%dK.csv" % T)) == LOW   # 本征点 n = 0 不进任何一段
    assert not list(tp.glob("MoS2_phonon_limited*")) and not list(tp.glob("MoS2_transport*"))   # 整轴旧图/旧表不再出
    thi = tp / "transport_carrier_tuning_1e17-1e21"
    tlo = tp / "transport_low_doping_1e10-1e16"
    assert (thi / "MoS2_transport_carrier_tuning.png").is_file()
    assert (tlo / "MoS2_transport_low_doping.png").is_file()
    rows = list(csv.DictReader((thi / "MoS2_transport_carrier_tuning_300K.csv").read_text().splitlines()))
    assert "ADP_POP" in {r["dataset"] for r in rows}            # transport 表含机制对照数据集
    assert sorted({abs(float(r["n_cell_cm-3"])) for r in rows}) == HIGH
    rows = list(csv.DictReader((tlo / "MoS2_transport_low_doping_300K.csv").read_text().splitlines()))
    assert sorted({abs(float(r["n_cell_cm-3"])) for r in rows}) == LOW


def test_only_carrier_tuning_without_low_doping(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, HIGH)
    g.write_transport_plots(cwd, out, True)
    tp = cwd.parent / "transport_plots"
    assert sorted(p.name for p in tp.glob("phonon_limited_*")) == ["phonon_limited_carrier_tuning_1e17-1e21"]
    assert sorted(p.name for p in tp.glob("transport_*")) == ["transport_carrier_tuning_1e17-1e21"]


def test_range_axes():
    g = _s83()
    r = g._phonon_limited_ranges(1e9)
    assert [x[0] for x in r] == ["carrier_tuning_1e17-1e21", "low_doping_1e9-1e16"]
    assert r[0][4] == ["1e21", "1e20", "1e19", "1e18", "1e17", "1e18", "1e19", "1e20", "1e21"]
    assert r[1][4][0] == "1e16" and r[1][4][len(r[1][4]) // 2] == "1e9" and r[1][4][-1] == "1e16"
    assert len(g._phonon_limited_ranges(1e17)) == 1


def _tool():
    spec = importlib.util.spec_from_file_location("mtp_0056", str(KE / "tools" / "make_transport_plots.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_tool_transport_split(tmp_path):
    """tools/make_transport_plots.py 的 transport 图（DPT 路线 + AMSET）同样按两段出文件夹。"""
    t = _tool()
    cwd, out = _material(tmp_path, LOW + HIGH)
    grid = json.loads((cwd / "step8.4_amset2d" / "transport.json").read_text())["doping"]
    temps = [300.0, 600.0]
    t3 = lambda v: [[v, 0, 0], [0, v, 0], [0, 0, v]]
    arr = lambda f: [[t3(f(n)) for _ in temps] for n in grid]
    (cwd / "step8.4_amset2d" / "transport.json").write_text(json.dumps({
        "doping": grid, "temperatures": temps, "conductivity": arr(lambda n: abs(n) * 1e-14 + 1e-6),
        "seebeck": arr(lambda n: -200.0 if n < 0 else 200.0), "electronic_thermal_conductivity": arr(lambda n: 1e-3),
        "mobility": {"overall": arr(lambda n: 100.0)}}))
    (cwd / "step8.4_amset2d" / "subset_ADP.json").write_text(json.dumps({"intrinsic": {
        "conductivity_S_m": arr(lambda n: abs(n) * 2e-14 + 1e-6), "seebeck_uV_K": arr(lambda n: 210.0),
        "electronic_thermal_conductivity_W_mK": arr(lambda n: 2e-3), "mobility_cm2_Vs_s": {"ADP": arr(lambda n: 300.0)}}}))
    (cwd / "step8.1_boltztrap").mkdir()
    nb = [s * 10 ** k for s in (-1, 1) for k in range(9, 22)]
    (cwd / "step8.1_boltztrap" / "boltztrap_crta.json").write_text(json.dumps({
        "twoD_c_over_t": 1.0, "temperatures_K": temps, "mu_rel_efermi_eV": list(range(len(nb))),
        "carrier_conc_cm-3": [nb for _ in temps], "sigma_over_tau_xx": [[abs(n) * 1e-1 for n in nb] for _ in temps],
        "seebeck_xx_V_per_K": [[-2e-4 if n < 0 else 2e-4 for n in nb] for _ in temps]}))
    (cwd / "step8.2_dpt").mkdir()
    (cwd / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": c, "by_direction": {"status": "ok", "x": {"tau_s": 1e-13, "mobility_cm2_Vs": m}}}
        for c, m in (("electron", 200.0), ("hole", 900.0))]}))
    od = tmp_path / "tp"
    od.mkdir()
    assert t.make_transport(str(cwd), str(od), "MoS2")
    hi = od / "transport_carrier_tuning_1e17-1e21"
    assert (hi / "MoS2_transport_carrier_tuning.png").is_file()
    assert not list(od.glob("MoS2_transport*"))
    rows = list(csv.reader((hi / "MoS2_transport_carrier_tuning_300K.csv").read_text().splitlines()))[1:]
    assert rows and min(abs(float(r[3])) for r in rows) >= 1e17 * 0.999
    dpt = [r for r in rows if r[1] == "DPT"]          # [0056-r4] BT2 只在 AMSET 掺杂范围内（不外推，1e21 那一档不画）
    assert dpt and max(abs(float(r[3])) for r in dpt) <= max(HIGH) * 1.001
    assert min(abs(float(r[3])) for r in dpt) >= min(HIGH) * 0.999
    lo = sorted(od.glob("transport_low_doping_*"))     # 原点 = AMSET 网格最小点所在十年
    assert len(lo) == 1
    rows = list(csv.reader((lo[0] / "MoS2_transport_low_doping_300K.csv").read_text().splitlines()))[1:]
    assert rows and max(abs(float(r[3])) for r in rows) <= 1e16 * 1.001


def test_mech_caption_lists_mechanisms():
    g = _s83()
    meta = {"overall_labels": ["ADP", "IMP", "POP", "ODP"],
            "sets": {"intrinsic": {"labels": ["ADP", "POP", "ODP"]}, "ADP_only": {"labels": ["ADP"]}}}
    lab = {"overall": "phonon + impurity (ADP+IMP+POP+ODP)", "intrinsic": "phonon-limited (ADP+POP+ODP)",
           "ADP_only": "ADP only (acoustic limit)", "BT2_DPTtau": "BT2×DPT-τ"}
    ds = [(k, "", "", "") for k in ("overall", "intrinsic", "ADP_only", "BT2_DPTtau")]
    cap = g._mech_caption(ds, meta, lab)
    assert cap[0] == "phonon + impurity: AMSET, ADP+IMP+POP+ODP (all mechanisms of this run)"
    assert cap[1] == "phonon-limited: AMSET, ADP+POP+ODP (IMP removed)"
    assert "Bardeen-Shockley" in cap[3]
    assert cap[4].startswith("ADP = ") and "IMP = ionized impurity" in cap[4] and "ODP = " in cap[4]


def test_gray_shading_for_rigid_band_limit(tmp_path):
    """in_interpret_range = NO 的档画灰底（只验证不崩、图注里有 Gray 说明）。"""
    g = _s83()
    cwd, out = _material(tmp_path, LOW + HIGH)
    tab = out / "intrinsic" / "intrinsic_transport_table.csv"
    rows = list(csv.DictReader(tab.read_text().splitlines()))
    for r in rows:
        r["in_interpret_range"] = "NO" if abs(float(r["n_cell_cm-3"])) >= 1e20 else "yes"
    with open(tab, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    seen = []
    import matplotlib.figure as mf
    orig = mf.Figure.text
    mf.Figure.text = lambda self, *a, **k: (seen.append(a[2]), orig(self, *a, **k))[1]
    try:
        g.write_transport_plots(cwd, out, True)
    finally:
        mf.Figure.text = orig
    assert any("Gray: beyond rigid-band range" in t for t in seen)
