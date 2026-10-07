"""组合技能的跨技能结果依赖（needs_results）+ 结果来源戳：按材料目录 + 结构指纹定位，
同名材料 / 旧结构的结果不会被误读。纯本地。"""
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import autozt  # noqa: E402
from autozt import workflow as W  # noqa: E402

SPEC = {"transport": ["ke-dft-cpu/result/step8_amset/transport.json",
                      "ke-dft-cpu/result/step8.4_amset2d/transport.json"],
        "kappa_L": ["kl-dft-cpu/result/step6_kappa/kappa_summary.json",
                    "kl-mlff-cpu/result/step4_kappa/kappa_summary.json"]}


def _material(root, name="Si", poscar="Si\n1.0\n"):
    lp = root / name
    lp.mkdir(parents=True, exist_ok=True)
    (lp / "POSCAR").write_text(poscar)
    return lp


def _fetched(lp, skill, step, fname, poscar_owner=None):
    """模拟 autozt 把 <skill> 的 <step> 拉回本地（含来源戳）。"""
    m = {"name": lp.name, "qualified_name": "proj/" + lp.name, "tt": skill,
         "lpath": str(poscar_owner or lp), "result_dir": str(lp / skill / "result")}
    d = lp / skill / "result" / step
    d.mkdir(parents=True, exist_ok=True)
    (d / fname).write_text("{}")
    W.write_result_origin(m, {"name": step, "label": step})
    return d


def _m(lp):
    return {"name": lp.name, "lpath": str(lp)}


def test_needs_results_by_material_dir_and_structure(tmp_path):
    lp = _material(tmp_path / "projA")
    chosen, miss = W.resolve_needs_results(_m(lp), SPEC)
    assert chosen == {} and len(miss) == 2 and "等 ke-dft-cpu" in miss[0]

    _fetched(lp, "ke-dft-cpu", "step8_amset", "transport.json")
    _fetched(lp, "kl-dft-cpu", "step6_kappa", "kappa_summary.json")
    chosen, miss = W.resolve_needs_results(_m(lp), SPEC)
    assert miss == [] and chosen["kappa_L"][1] == SPEC["kappa_L"][0]

    # 结构改了（同名材料重新算）：旧结果不再满足依赖
    (lp / "POSCAR").write_text("Si\n1.01\n")
    chosen, miss = W.resolve_needs_results(_m(lp), SPEC)
    assert chosen == {} and all("结构已变" in x for x in miss)


def test_same_name_in_another_project_is_not_used(tmp_path):
    a = _material(tmp_path / "projA")
    b = _material(tmp_path / "projB")                 # 同名、同结构，另一个项目
    _fetched(b, "ke-dft-cpu", "step8_amset", "transport.json")
    # 把 B 的结果拷到 A 的目录里冒充（来源戳仍指向 B）
    src = b / "ke-dft-cpu"
    os.makedirs(a / "ke-dft-cpu", exist_ok=True)
    import shutil
    shutil.copytree(src / "result", a / "ke-dft-cpu" / "result")
    chosen, miss = W.resolve_needs_results(_m(a), {"transport": SPEC["transport"]})
    assert chosen == {} and "另一个材料目录" in miss[0]


def test_legacy_results_without_origin_must_be_refetched(tmp_path):
    lp = _material(tmp_path / "p")
    d = lp / "kl-dft-cpu" / "result" / "step6_kappa"
    d.mkdir(parents=True)
    (d / "kappa_summary.json").write_text("{}")       # 旧版本拉回的，没有来源戳
    _c, miss = W.resolve_needs_results(_m(lp), {"kappa_L": SPEC["kappa_L"]})
    assert "缺来源戳" in miss[0]


def test_alternatives_take_first_valid(tmp_path):
    lp = _material(tmp_path / "p")
    _fetched(lp, "kl-mlff-cpu", "step4_kappa", "kappa_summary.json")
    chosen, miss = W.resolve_needs_results(_m(lp), {"kappa_L": SPEC["kappa_L"]})
    assert miss == [] and chosen["kappa_L"][1].startswith("kl-mlff-cpu/")
    _fetched(lp, "kl-dft-cpu", "step6_kappa", "kappa_summary.json")
    chosen, _ = W.resolve_needs_results(_m(lp), {"kappa_L": SPEC["kappa_L"]})
    assert chosen["kappa_L"][1].startswith("kl-dft-cpu/")   # 列表顺序 = 优先级


def test_dag_blocks_until_results_and_graph_says_why(tmp_path, monkeypatch):
    lp = _material(tmp_path / "p")
    t = {"key": "zt-dft-cpu", "steps": [{"name": "step20_zt", "label": "S20_zt",
                                         "needs": [], "needs_results": SPEC}]}
    m = {"name": "Si", "lpath": str(lp), "tt": "zt-dft-cpu",
         "steps": [{"name": "step20_zt", "label": "S20_zt", "job": None, "done": False,
                    "has_outcar": False, "has_slurm_out": False, "has_incar": False}]}
    W._dag_recompute(t, m)
    s = m["steps"][0]
    assert s["kind"] == "WAIT" and len(s["_missing_results"]) == 2
    from autozt import taskgraph
    g = taskgraph.material_graph({"key": "zt-dft-cpu"}, m)
    assert "等跨技能结果" in g["nodes"][0]["status"]
    _fetched(lp, "ke-dft-cpu", "step8_amset", "transport.json")
    _fetched(lp, "kl-dft-cpu", "step6_kappa", "kappa_summary.json")
    W._dag_recompute(t, m)
    assert m["steps"][0]["kind"] == "PREP" and m["steps"][0]["_missing_results"] == []


def test_origin_written_with_fetch_receipt(tmp_path):
    lp = _material(tmp_path / "p")
    m = {"name": "Si", "qualified_name": "p/Si", "tt": "ke-dft-cpu", "lpath": str(lp),
         "result_dir": str(lp / "ke-dft-cpu" / "result"), "fetch_files": []}
    s = {"name": "step8_amset", "label": "S8_kappa", "done": True, "job": None,
         "dir": "/remote/Si/ke-dft-cpu/step8_amset"}
    W._fetch_receipt_write({}, m, s)
    o = json.loads((lp / "ke-dft-cpu" / "result" / "step8_amset" / W.ORIGIN_NAME).read_text())
    assert o["material_dir"] == os.path.realpath(lp) and o["skill"] == "ke-dft-cpu"
    assert o["poscar_sha256"] == W.material_poscar_sha({"lpath": str(lp)})
