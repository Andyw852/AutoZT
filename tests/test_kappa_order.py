# -*- coding: utf-8 -*-
"""fit-fc-thermal S2_kappa 收敛顺序：q 网格（参考截断）→ 截断（收敛网格上）→ 温度。"""
import importlib.util
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SK = os.path.join(ROOT, "skill", "_common", "fcfit")   # 共享引擎（原 skill/fit-fc-thermal/）
sys.path.insert(0, os.path.join(ROOT, "skill", "_common"))


def _load():
    sys.path.insert(0, SK)
    spec = importlib.util.spec_from_file_location("kappa_driver_t", os.path.join(SK, "kappa_driver.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _summary(cut, L):
    k = 10.0 / cut                       # 截断越大 κ 越小，模拟平台前的下降
    return {"mesh": "%d %d 1" % (L, L), "temperatures": [100.0, 300.0, 800.0],
            "kappa_inplane_300K": k, "kappa_300K_xx_yy_zz": [k, k, 0.0]}


def _run(tmp_path, monkeypatch, chosen):
    kd = _load()
    calls = []

    def fake_conv(out, yaml, fc2, fc3, cfg, src, first=None):
        cut = kd._cut_of_tag(fc2.split("/")[0][5:])
        calls.append(("conv", cut, None if first is None else first[0]))
        s = _summary(cut, 70)
        s["mesh_converged"] = True
        s["mesh_convergence"] = {"records": [{"length": 40.0}, {"length": 70.0}]}
        return s

    def fake_bte(out, yaml, fc2, fc3, cfg, src, write_kappa, spec=None):
        cut = kd._cut_of_tag(fc2.split("/")[0][5:])
        calls.append(("bte", cut, spec))
        return _summary(cut, int(spec[1]) if spec else 40)

    monkeypatch.setattr(kd, "_converge_mesh", fake_conv)
    monkeypatch.setattr(kd, "_run_bte", fake_bte)
    monkeypatch.setattr(kd, "_write_kappa_vs_cutoff", lambda *a, **k: None)
    monkeypatch.setattr(kd, "_select", lambda out, cfg, per: (
        chosen, {"chosen_cut": chosen}, dict(per)))
    (tmp_path / "kappa_config.json").write_text(json.dumps({
        "disp_yaml": "x.yaml", "scan": True, "mesh": "auto", "mesh_length": 40,
        "cut_scans": ["4p00", "5p00", "6p00"], "nominal_cut3_A": 5.0}))
    monkeypatch.chdir(tmp_path)
    kd.main()
    return calls, json.loads((tmp_path / "kappa_summary.json").read_text())


def test_mesh_first_then_cutoffs_at_converged_mesh(tmp_path, monkeypatch):
    calls, s = _run(tmp_path, monkeypatch, chosen=5.0)
    assert calls[0] == ("conv", 5.0, None)                    # ① 参考=标称截断，先收敛网格
    assert {c[1] for c in calls[1:]} == {4.0, 6.0}            # ② 其余档
    assert all(c[0] == "bte" and c[2] == ("length", 70.0) for c in calls[1:])  # 都在收敛网格上
    assert s["chosen_cutoff_A"] == 5.0 and s["mesh"] == "70 70 1"
    assert s["temperatures"] == [100.0, 300.0, 800.0]          # ③ 完整 κ(T)
    assert s["convergence_order"][0].startswith("q_mesh")


def test_chosen_not_reference_rechecks_mesh(tmp_path, monkeypatch):
    calls, s = _run(tmp_path, monkeypatch, chosen=4.0)
    assert calls[-1] == ("conv", 4.0, 70.0)                   # 从收敛长度起复核
    assert s["chosen_cutoff_A"] == 4.0


def test_reference_falls_back_to_largest():
    kd = _load()
    assert kd._reference_tag(["4p00", "6p50", "5p00"], None) == "6p50"
    assert kd._reference_tag(["4p00", "6p50"], 5.0) == "6p50"
    assert kd._reference_tag(["4p00", "5p00"], 5.0) == "5p00"
