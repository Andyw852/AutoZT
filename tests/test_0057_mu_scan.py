"""[0057] S8.3 transport_plots/mu_scan/：BT2 与 AMSET 放在同一根 μ − E_VBM 轴上（原 bt2_muscan/ + amset_scan/）。
换算：BT2 x = μ_rel − VBM_rel（同为 S8.1 的 E_F 零点）；AMSET x = fermi_levels − E_v（undoped_ni.json，剪刀后）。
合成数据，纯本地。"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _s83():
    for d in (ROOT, KE, KE / "step8.4_amset2d", ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("s83_0057", str(KE / "step8.3_output" / "gen_step13_output.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


EV, EC = -1.20, 0.60          # AMSET 带边（eV，与 fermi_levels 同零点）
DOP = [-1e19, -1e17, 0.0, 1e17, 1e19]
EF = [0.55, 0.40, -0.30, -1.05, -1.15]   # 每个掺杂档的 E_F（600 K）


def _material(tmp_path, bt2_edges=None, scissor=None):
    cwd = tmp_path / "MoS2" / "ke-dft-cpu"
    out = cwd / "step8.3_output"
    (out / "intrinsic").mkdir(parents=True)
    run = cwd / "step8.4_amset2d"
    run.mkdir()
    (run / "2d_correction.json").write_text(json.dumps({"cell_c_A": 20.0, "elastic_rescale_factor_c_over_t": 3.0}))
    (run / "transport.json").write_text(json.dumps({"doping": DOP, "temperatures": [300.0, 600.0],
                                                    "fermi_levels": [[e + 0.01, e] for e in EF]}))
    (run / "undoped_ni.json").write_text(json.dumps({"temperatures_K": [300.0, 600.0], "E_F_eV": [-0.29, -0.30],
                                                     "n_i_cm-3": [1.0, 1e5], "E_v_eV": [EV, EV], "E_c_eV": [EC, EC]}))
    cols = ["T_K", "dataset", "type", "n_cell_cm-3", "in_interpret_range", "S_uV/K", "sigma_S/m", "PF_W/mK2",
            "kappa_e_W/mK", "mu_cm2/Vs"]
    with open(out / "intrinsic" / "intrinsic_transport_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for ds in ("overall", "intrinsic", "BT2_DPTtau"):
            for n in DOP:
                if n:
                    w.writerow([600.0, ds, "n" if n < 0 else "p", n, "yes", -200 if n < 0 else 200, 1e4, 1e-3, 0.1, 100])
    mu = [-2.0 + 0.01 * i for i in range(401)]           # μ_rel，-2 .. +2 eV
    j = {"temperatures_K": [300.0, 600.0], "mu_rel_efermi_eV": mu,
         "carrier_conc_cm-3": [[(-1e18 if m > 0 else 1e18) for m in mu]] * 2,
         "seebeck_V_per_K": [[(-2e-4 if m > 0 else 2e-4) for m in mu]] * 2,
         "lorenz_WOhm_per_K2": [[2e-8] * len(mu)] * 2,
         "sigma_over_tau_S_per_m_s": [[1e18] * len(mu)] * 2,
         "kappa_e_over_tau_W_per_m_K_s": [[1e10] * len(mu)] * 2}
    if bt2_edges:
        j["band_edges_rel_efermi_eV"] = bt2_edges
    if scissor:
        j["scissor"] = scissor
    (cwd / "step8.1_boltztrap").mkdir()
    (cwd / "step8.1_boltztrap" / "boltztrap_crta.json").write_text(json.dumps(j))
    return cwd, out


def test_bt2_edges_sources():
    g = _s83()
    assert g._bt2_edges({"band_edges_rel_efermi_eV": {"vbm": -0.1, "cbm": 1.6}})[:2] == (-0.1, 1.6)
    assert g._bt2_edges({"scissor": {"applied": True, "target_gap_eV": 1.8}})[:2] == (-0.9, 0.9)   # E_F 在带隙正中
    assert g._bt2_edges({"scissor": {"applied": False}})[0] is None


def test_mu_scan_conversion(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, scissor={"applied": True, "target_gap_eV": 1.8})
    d = g.mu_scan_data(cwd, out, True)
    assert d["T"] == 600.0 and d["amset_vbm"] == EV and d["bt2_vbm"] == -0.9
    am = {(r["dataset"], float(r["n_cell_cm-3"])): r["mu_minus_VBM_eV"] for r in d["amset"]}
    for n, ef in zip(DOP, EF):
        if n:
            assert am[("overall", n)] == pytest.approx(ef - EV)       # 600 K 那一列 fermi_levels，减 AMSET VBM
    assert ("BT2_DPTtau", 1e17) not in am                            # BT2 不从表里取点，画自己的连续曲线
    bt = {round(r["mu_rel_efermi_eV"], 4): r["mu_minus_VBM_eV"] for r in d["bt2"]}
    assert bt[0.0] == pytest.approx(0.9) and bt[-0.9] == pytest.approx(0.0)
    assert not any("WARN" in n for n in d["notes"])                  # 本征点 E_F 两份记录一致


def test_mu_scan_files_and_old_dirs_gone(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, bt2_edges={"vbm": -0.9, "cbm": 0.9})
    g.write_transport_plots(cwd, out, True)
    tp = cwd / "transport_plots"
    assert (tp / "mu_scan" / "MoS2_mu_scan_600K.png").is_file()
    rows = [r for r in csv.DictReader([l for l in (tp / "mu_scan" / "MoS2_mu_scan_600K.csv").read_text().splitlines()
                                       if not l.startswith("#")])]
    assert {r["source"] for r in rows} == {"BT2", "AMSET"}
    assert not (tp / "bt2_muscan").exists() and not (tp / "amset_scan").exists() and not (tp / "datasets").exists()


def test_bt2_without_reference_is_skipped(tmp_path):
    """没做剪刀、也没存带边：μ 零点离 VBM 多远不知道 —— 不画 BT2，不猜。"""
    g = _s83()
    cwd, out = _material(tmp_path)
    d = g.mu_scan_data(cwd, out, True)
    assert d["bt2"] == [] and d["amset"]
    assert any(n.startswith("BT2：") for n in d["notes"])
