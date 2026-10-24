"""[0056][0057] S8.3 掺杂扫描按浓度拆两段，各一个文件夹（图 + 每温度表）：
<前缀>_carrier_tuning_1e16-1e21/（载流子调控区）与 <前缀>_low_doping_<NLO>-1e15/（低浓度段，有 DOPING_LOW 才有；
[0077] 前是 1e17-1e21 / <NLO>-1e16）；
[0057] 前缀 transport_models（原 phonon_limited，基础 4 条 x 5 参数）与 mobility_mechanisms（原 transport，只画 μ：
单机制 + 总值 + DPT 带边）。合成 intrinsic_transport_table.csv + transport.json，纯本地。"""
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

LOW = [1e11, 1e12, 1e13, 1e14, 1e15]                  # [0077] 出厂低浓度段（近本征）
HIGH = [1e16, 1e17, 1e18, 1e19, 1e20, 1e21]           # [0077] 出厂掺杂段


def _material(tmp_path, dop):
    cwd = tmp_path / "MoS2" / "ke-dft-cpu"
    out = cwd / "step8.3_output"
    (out / "intrinsic").mkdir(parents=True)
    (cwd / "step8.4_amset2d").mkdir()
    grid = sorted([-x for x in dop] + [0.0] + dop)
    temps = [300.0, 600.0]

    def _t(f):                                       # [nd][nT][3][3]，面内各向同性
        return [[[[f(n, T), 0, 0], [0, f(n, T), 0], [0, 0, 0]] for T in temps] for n in grid]
    mob = {"ADP": _t(lambda n, T: 500.0 * 300 / T), "POP": _t(lambda n, T: 300.0 * 300 / T),
           "IMP": _t(lambda n, T: 1e3 * (1e21 / max(abs(n), 1.0)) ** 0.5),
           "overall": _t(lambda n, T: 1.0 / (T / 150000 + T / 90000 + 1.0 / (1e3 * (1e21 / max(abs(n), 1.0)) ** 0.5)))}
    (cwd / "step8.4_amset2d" / "transport.json").write_text(json.dumps(
        {"doping": grid, "temperatures": temps, "mobility": mob}))
    with open(out / "intrinsic" / "intrinsic_transport_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for T in (300.0, 600.0):
            for ds in ("overall", "intrinsic", "ADP_only", "BT2_DPTtau", "ADP_POP"):
                for n in grid:
                    w.writerow([T, ds, "n" if n < 0 else "p", n, abs(n) * 1e-7, "", "", "", "yes",
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
    tp = cwd / "transport_plots"
    hi = tp / "transport_models_carrier_tuning_1e16-1e21"
    lo = tp / "transport_models_low_doping_1e11-1e15"
    assert (hi / "MoS2_transport_models_carrier_tuning.png").is_file()
    assert (lo / "MoS2_transport_models_low_doping.png").is_file()
    for T in (300, 600):
        assert _ns(hi / ("MoS2_transport_models_carrier_tuning_%dK.csv" % T)) == HIGH
        assert _ns(lo / ("MoS2_transport_models_low_doping_%dK.csv" % T)) == LOW   # 本征点 n = 0 不进任何一段
    assert not list(tp.glob("MoS2_*.png"))                    # 整轴旧图不再出
    assert not list(tp.glob("phonon_limited_*")) and not list(tp.glob("transport_carrier*"))   # [0057] 旧目录名不再出
    mhi = tp / "mobility_mechanisms_carrier_tuning_1e16-1e21"
    mlo = tp / "mobility_mechanisms_low_doping_1e11-1e15"
    assert (mhi / "MoS2_mobility_mechanisms_carrier_tuning.png").is_file()
    assert (mlo / "MoS2_mobility_mechanisms_low_doping.png").is_file()
    rows = list(csv.DictReader((mhi / "MoS2_mobility_mechanisms_carrier_tuning_300K.csv").read_text().splitlines()))
    assert {"mu_ADP_only_cm2/Vs", "mu_POP_only_cm2/Vs", "mu_IMP_only_cm2/Vs", "mu_total_cm2/Vs"} <= set(rows[0])
    assert sorted({abs(float(r["n_cell_cm-3"])) for r in rows}) == HIGH
    assert float(rows[0]["mu_ADP_only_cm2/Vs"]) == pytest.approx(500.0)
    assert all(r["n_2D_cm-2"] for r in rows)                  # 浓度口径列从本征表配上（[0057-r1]）
    rows = list(csv.DictReader((mlo / "MoS2_mobility_mechanisms_low_doping_600K.csv").read_text().splitlines()))
    assert sorted({abs(float(r["n_cell_cm-3"])) for r in rows}) == LOW
    assert float(rows[0]["mu_ADP_only_cm2/Vs"]) == pytest.approx(250.0)


def test_only_carrier_tuning_without_low_doping(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, HIGH)
    g.write_transport_plots(cwd, out, True)
    tp = cwd / "transport_plots"
    assert sorted(p.name for p in tp.glob("transport_models_*")) == ["transport_models_carrier_tuning_1e16-1e21"]
    assert sorted(p.name for p in tp.glob("mobility_mechanisms_*")) == ["mobility_mechanisms_carrier_tuning_1e16-1e21"]


def test_range_axes():
    g = _s83()
    r = g._phonon_limited_ranges(1e9)
    assert [x[0] for x in r] == ["carrier_tuning_1e16-1e21", "low_doping_1e9-1e15"]
    assert r[0][4] == ["1e21", "1e20", "1e19", "1e18", "1e17", "1e16", "1e17", "1e18", "1e19", "1e20", "1e21"]
    assert r[1][4][0] == "1e15" and r[1][4][len(r[1][4]) // 2] == "1e9" and r[1][4][-1] == "1e15"
    assert [x[0] for x in g._phonon_limited_ranges(1e11)] == ["carrier_tuning_1e16-1e21", "low_doping_1e11-1e15"]
    assert len(g._phonon_limited_ranges(1e16)) == 1
    assert g._phonon_limited_ranges(1e11, 1e22)[0][0] == "carrier_tuning_1e16-1e22"   # 网格超出 1e21 时上端跟着放


def _tool():
    spec = importlib.util.spec_from_file_location("mtp_0056", str(KE / "tools" / "make_transport_plots.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _add_dpt_rows(out):
    tab = out / "intrinsic" / "intrinsic_transport_table.csv"
    with open(tab, "a", newline="") as fh:
        w = csv.writer(fh)
        for T in (300.0, 600.0):
            for typ, mu in (("n", 200.0 * 300 / T), ("p", 900.0 * 300 / T)):
                w.writerow([T, "DPT", typ] + [""] * 10 + [mu, ""])


def test_tool_is_s83(tmp_path):
    """[0056-r5] 用户定「统一」：独立脚本 = S8.3 的 write_transport_plots，产物逐字相同。"""
    g, t = _s83(), _tool()
    cwd, out = _material(tmp_path, LOW + HIGH)
    _add_dpt_rows(out)
    (cwd / "step1_opt").mkdir()
    (cwd / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\n")      # 独立脚本按 S1 判维（同 S8.3 main）
    g.write_transport_plots(cwd, out, True)
    od = tmp_path / "tool_out"
    assert t.run(str(cwd), str(od), "MoS2")
    tp = cwd / "transport_plots"
    a = sorted(str(x.relative_to(tp)) for x in tp.rglob("*.csv"))
    b = sorted(str(x.relative_to(od)) for x in od.rglob("*.csv"))
    assert a == b and a
    for rel in a:
        assert (tp / rel).read_text() == (od / rel).read_text(), rel
    assert (od / "mobility_mechanisms_carrier_tuning_1e16-1e21" / "MoS2_mobility_mechanisms_carrier_tuning.png").is_file()


def test_tool_skips_without_s83_table(tmp_path):
    t = _tool()
    cwd = tmp_path / "MoS2" / "ke-dft-cpu"
    cwd.mkdir(parents=True)
    assert t.run(str(cwd)) is False


def test_dpt_lines_half_side_and_label_headroom(tmp_path):
    """DPT 带边 μ：电子线只在左半（x ≤ 0）、空穴只在右半（x ≥ 0）；μ 轴顶部留空白放标注，曲线不进标注区
    （[0057] 在 mobility_mechanisms 图里；低浓度 IMP 远超总值的部分不画，也不能顶进标注区）。"""
    g = _s83()
    cwd, out = _material(tmp_path, LOW + HIGH)
    _add_dpt_rows(out)
    import matplotlib.figure as mf
    got = []
    orig = mf.Figure.savefig

    def spy(self, *a, **k):
        if "mobility_mechanisms_carrier_tuning" in str(a[0]):
            got.append(self.axes[0])                  # 第一个温度的面板
        return orig(self, *a, **k)
    mf.Figure.savefig = spy
    try:
        g.write_transport_plots(cwd, out, True)
    finally:
        mf.Figure.savefig = orig
    ax = got[0]
    dash = [l for l in ax.get_lines() if l.get_color() in ("navy", "tab:purple")]
    assert len(dash) == 2
    for l in dash:
        x = l.get_xdata()
        assert max(x) <= 0 or min(x) >= 0                   # 不跨过中线
    lo, hi = ax.get_ylim()
    import math
    top_data = max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata()))
    frac = (math.log10(top_data) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    assert frac <= 0.80                                     # 标注区（顶部约 20%）里没有曲线
    texts = [tx.get_text() for tx in ax.texts]
    assert any("mu_{e}" in tx for tx in texts) and any("mu_{h}" in tx for tx in texts)


def test_mech_caption_lists_mechanisms():
    g = _s83()
    meta = {"overall_labels": ["ADP", "IMP", "POP", "ODP"],
            "sets": {"intrinsic": {"labels": ["ADP", "POP", "ODP"]}, "ADP_only": {"labels": ["ADP"]}}}
    lab = dict(g.amset_labels(meta), BT2_DPTtau=g.BT2_LABEL)          # [0079] 图例 = 实际散射机制组合
    assert lab["overall"] == "AMSET ADP+POP+ODP+IMP" and lab["BT2_DPTtau"] == "BT2×DPT-τ (ADP)"
    ds = [(k, "", "", "") for k in ("overall", "intrinsic", "ADP_only", "BT2_DPTtau")]
    cap = g._mech_caption(ds, meta, lab)                              # [0080] 精简：图例已是机制组合
    assert cap[0].startswith("AMSET: scattering mechanisms as labeled")
    assert cap[1].startswith("BT2×DPT-τ (ADP): BoltzTraP2") and "Bardeen-Shockley" in cap[1] and len(cap) == 3
    assert not any(w in c for c in cap for w in ("phonon +", "phonon-limited", "partial phonon", "together"))
    assert cap[2] == ("ADP = acoustic deformation potential; POP = polar optical phonon; "
                      "ODP = optical deformation potential; IMP = ionized impurity")   # 按 MECH_ORDER


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
    assert any("Gray: beyond the rigid-band range" in t for t in seen)
