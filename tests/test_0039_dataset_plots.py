"""[0039] S8.3 transport_plots/datasets/：每个数据集（overall / intrinsic / ADP only / 子集 / BT2）单独一张
σ、S、PF、κe 随浓度图（每温度一条线）+ 同名 csv。合成 intrinsic_transport_table.csv，纯本地。"""
import csv
import importlib.util
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
    spec = importlib.util.spec_from_file_location("s83_0039", str(KE / "step8.3_output" / "gen_step13_output.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _material(tmp_path, datasets):
    cwd = tmp_path / "MoS2" / "ke-dft-cpu"
    out = cwd / "step8.3_output"
    (out / "intrinsic").mkdir(parents=True)
    with open(out / "intrinsic" / "intrinsic_transport_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for T in (300.0, 600.0):
            for ds in datasets:
                for n in (-1e20, -1e19, -1e18, 1e18, 1e19, 1e20):
                    w.writerow([T, ds, "n" if n < 0 else "p", n, "", "", "", "", "yes",
                                -200 if n < 0 else 200, 1e4 * abs(n) / 1e18, 1e-3, 0.1, 100, ""])
    return cwd, out


def test_one_figure_per_dataset(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, ["overall", "intrinsic", "ADP_only", "BT2_DPTtau"])
    g.write_transport_plots(cwd, out, True)
    ddir = cwd.parent / "transport_plots" / "datasets"
    for ds, fn in (("overall", "overall"), ("intrinsic", "noIMP"), ("ADP_only", "ADP_only"),
                   ("BT2_DPTtau", "BT2_DPTtau")):
        assert (ddir / ("MoS2_%s.png" % fn)).is_file(), ds          # [0040] 去 IMP 子集文件名 noIMP
        rows = (ddir / ("MoS2_%s.csv" % fn)).read_text().splitlines()
        assert rows[0].startswith("# dataset=%s" % ds)
        assert rows[1].split(",")[:3] == ["T_K", "type", "n_cell_cm-3"]
        assert len(rows) == 2 + 2 * 6                      # 2 温度 x 6 浓度
    assert (cwd.parent / "transport_plots" / "MoS2_transport.png").is_file()   # 主图照旧


def test_stale_dataset_removed(tmp_path):
    g = _s83()
    cwd, out = _material(tmp_path, ["overall", "intrinsic", "ADP_only"])
    ddir = cwd.parent / "transport_plots" / "datasets"
    ddir.mkdir(parents=True)
    (ddir / "MoS2_ADP_POP.png").write_bytes(b"old")
    g.write_transport_plots(cwd, out, True)
    assert not (ddir / "MoS2_ADP_POP.png").exists()
    assert (ddir / "MoS2_noIMP.png").is_file()
