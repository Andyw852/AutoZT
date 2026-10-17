"""[0056] S8.3 phonon-limited 掺杂扫描按浓度拆两段，各一个文件夹（图 + 每温度表）：
carrier_tuning_1e17-1e21/（载流子调控区）与 low_doping_<NLO>-1e16/（低浓度段，有 DOPING_LOW 才有）。
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
    assert not list(tp.glob("MoS2_phonon_limited*"))          # 整轴旧图/旧表不再出
    assert (tp / "MoS2_transport.png").is_file() and (tp / "MoS2_transport_data_300K.csv").is_file()


def test_only_carrier_tuning_without_low_doping(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, HIGH)
    g.write_transport_plots(cwd, out, True)
    tp = cwd.parent / "transport_plots"
    assert sorted(p.name for p in tp.glob("phonon_limited_*")) == ["phonon_limited_carrier_tuning_1e17-1e21"]


def test_range_axes():
    g = _s83()
    r = g._phonon_limited_ranges(1e9)
    assert [x[0] for x in r] == ["carrier_tuning_1e17-1e21", "low_doping_1e9-1e16"]
    assert r[0][4] == ["1e21", "1e20", "1e19", "1e18", "1e17", "1e18", "1e19", "1e20", "1e21"]
    assert r[1][4][0] == "1e16" and r[1][4][len(r[1][4]) // 2] == "1e9" and r[1][4][-1] == "1e16"
    assert len(g._phonon_limited_ranges(1e17)) == 1
