"""[0017] tools/mobility_vs_dpt.py：2D 的 DPT 迁移率取 by_direction 面内平均（与 S8.3、JAP Table I 同口径），
没有 by_direction 才退回 header 标量；来源写进表头与 --json。纯本地。"""
import contextlib
import importlib.util
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _tool():
    for d in (KE, ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("mvd_0017", str(KE / "tools" / "mobility_vs_dpt.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _diag(v):
    return [[v, 0, 0], [0, v, 0], [0, 0, 1.0]]


def _material(tmp_path, with_bd=True):
    run = tmp_path / "step8.4_amset2d"
    run.mkdir()
    dop, temps = [-1e17, 1e17], [300.0]
    tj = {"doping": dop, "temperatures": temps,
          "seebeck": [[_diag(-300.0)], [_diag(300.0)]],
          "mobility": {"overall": [[_diag(30.0)], [_diag(90.0)]],
                       "ADP": [[_diag(60.0)], [_diag(200.0)]]}}
    (run / "transport.json").write_text(json.dumps(tj))
    (run / "2d_correction.json").write_text(json.dumps({"cell_c_A": 20.0}))
    res = []
    for car, head, bx, by in (("electron", 150.0, 100.0, 120.0), ("hole", 500.0, 300.0, 340.0)):
        r = {"carrier": car, "mobility_cm2_Vs": head, "inputs": {"m_eff_m0": 0.9}}
        if with_bd:
            r["by_direction"] = {"status": "ok", "m_d_m0": 0.95,
                                 "x": {"mobility_cm2_Vs": bx}, "y": {"mobility_cm2_Vs": by}}
        res.append(r)
    (tmp_path / "step8.2_dpt").mkdir()
    (tmp_path / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": res}))
    return tmp_path


def test_dpt_uses_by_direction_inplane_mean(tmp_path):
    M = _tool()
    mat = _material(tmp_path)
    out = tmp_path / "o.json"
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        rc = M.main([str(mat), "--json", str(out)])
    assert rc == 0 and "by_direction" in buf.getvalue()
    rows = {r["carrier"]: r for r in json.loads(out.read_text())}
    assert rows["electron"]["DPT"] == 110.0 and rows["hole"]["DPT"] == 320.0
    assert rows["electron"]["ADP_over_DPT"] == round(60.0 / 110.0, 3)
    assert rows["electron"]["DPT_source"] == "by_direction"
    assert M.dpt_masses(mat / "step8.2_dpt" / "dpt_result.json")["electron"] == 0.95


def test_dpt_falls_back_to_header_without_by_direction(tmp_path):
    M = _tool()
    mat = _material(tmp_path, with_bd=False)
    out = tmp_path / "o.json"
    with contextlib.redirect_stdout(io.StringIO()):
        M.main([str(mat), "--json", str(out)])
    rows = {r["carrier"]: r for r in json.loads(out.read_text())}
    assert rows["electron"]["DPT"] == 150.0 and rows["electron"]["DPT_source"] == "header"
