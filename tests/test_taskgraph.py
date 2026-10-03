"""任务图（autozt graph / MCP task_graph / 审批卡片里的图）——纯本地单元测试。"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from autozt import agentgate as G  # noqa: E402
from autozt import taskgraph as TG  # noqa: E402


def _st(name, label, kind, deps, job=None, diag="", lt=None, missing=None):
    s = {"name": name, "label": label, "kind": kind, "_deps": deps, "job": job,
         "diag": diag, "label_txt": lt or kind}
    if missing:
        s["_missing_deps"] = missing
    return s


def _mat():
    return {"name": "Bi2Te3", "qualified_name": "proj/Bi2Te3", "hpc_name": "jzzn",
            "tt": "ke-dft-cpu", "steps": [
                _st("step1_opt", "S1_opt", "OK", []),
                _st("step3_uniform", "S3_uniform", "OK", ["step1_opt"]),
                _st("step4_wave", "S4_wave", "R", ["step3_uniform"],
                    {"id": "5756", "state": "R", "info": "cu12", "time": "2:13:04"}),
                _st("step5_dielect", "S5_dielect", "FAIL", ["step1_opt"],
                    diag="force not converged"),
                _st("step6_elastic", "S6_elastic", "PD", ["step1_opt"],
                    {"id": "5757", "state": "PD", "info": "(Resources)"}),
                _st("step7_deform", "S7_deform", "OK", ["step1_opt"], lt="10/10"),
                _st("step7b_deform_read", "S7.1_read", "SCANCEL", ["step7_deform"]),
                _st("step8_kappa", "S8_kappa", "WAIT",
                    ["step4_wave", "step5_dielect", "step6_elastic", "step7b_deform_read"],
                    missing=["S2.3_hseplot"]),
            ]}


def test_material_graph_nodes_edges_and_status():
    g = TG.material_graph({"key": "ke-dft-cpu"}, _mat())
    assert g["view"] == "material" and g["id"] == "proj/Bi2Te3"
    assert g["done"] == 3 and g["total"] == 8
    by = {n["label"]: n for n in g["nodes"]}
    assert by["S4_wave"]["status"].startswith("运行 #5756 @cu12")
    assert by["S6_elastic"]["status"] == "排队 #5757 (Resources)"
    assert "force not converged" in by["S5_dielect"]["status"]
    assert by["S7_deform"]["status"] == "完成 [10/10]"           # 扇出计数
    # 等上游只列没完成的依赖；主父节点 = 最深的依赖，其余进 extra_deps
    assert by["S8_kappa"]["pending_deps"] == ["S4_wave", "S5_dielect", "S6_elastic",
                                              "S7.1_read"]
    assert by["S8_kappa"]["parent"] == "S7.1_read"
    assert ["S1_opt", "S3_uniform"] in g["edges"] and len(g["edges"]) == 10


def test_render_text_tree_marks_and_ascii_fallback():
    g = TG.material_graph({"key": "ke-dft-cpu"}, _mat())
    txt = TG.render_text(g, marks={"S6_elastic": "◀ 本次：取消"}, ascii_mode=False)
    lines = txt.splitlines()
    assert lines[0].startswith("ke-dft-cpu · proj/Bi2Te3 · jzzn   3/8 完成")
    assert any(ln.startswith("├─◷ S6_elastic") and ln.endswith("◀ 本次：取消") for ln in lines)
    assert any("└─· S8_kappa" in ln and "缺 S2.3_hseplot" in ln for ln in lines)
    assert "另依赖" not in txt                     # WAIT 的状态里已经列出未完成上游
    asc = TG.render_text(g, ascii_mode=True)
    icons = [ln.split()[0] for ln in asc.splitlines()[1:-1]]
    assert all(all(ord(c) < 128 for c in tok) for tok in icons), icons


def test_mermaid_keeps_line_breaks_and_escapes_quotes():
    m = _mat()
    m["steps"][3]["diag"] = 'bad "INCAR" <tag>'
    g = TG.material_graph({"key": "x"}, m)
    mer = TG.render_mermaid(g, marks={"S6_elastic": "◀ 本次：取消"})
    assert mer.startswith("flowchart TD")
    assert "<br/>" in mer and "&lt;br/&gt;" not in mer
    assert "#quot;INCAR#quot;" in mer and "&lt;tag&gt;" in mer
    assert "class n4 target" in mer and "n0 --> n1" in mer


def test_skill_graph_counts_and_descendants():
    m1, m2 = _mat(), _mat()
    m2["steps"] = [dict(s, kind="OK") for s in m2["steps"]]
    sk = TG.skill_graph({"key": "ke-dft-cpu", "materials": [m1, m2]})
    by = {n["label"]: n for n in sk["nodes"]}
    assert sk["materials"] == 2 and sk["materials_done"] == 1
    assert by["S5_dielect"]["counts"] == {"FAIL": 1, "OK": 1}
    assert sk["fails"][0]["step"] == "S5_dielect"
    assert "S5_dielect" in TG.render_text(sk, ascii_mode=False)
    g = TG.material_graph({"key": "x"}, m1)
    assert TG.descendants(g, ["S7_deform"]) == ["S7.1_read", "S8_kappa"]


def test_graph_without_dag_info_falls_back_to_previous_step():
    m = {"name": "A", "steps": [{"name": "s1", "label": "S1", "kind": "OK"},
                                {"name": "s2", "label": "S2", "kind": "R"}]}
    g = TG.material_graph({"key": "x"}, m)
    assert g["edges"] == [["S1", "S2"]]


def test_card_marks_target_and_downstream():
    g = TG.material_graph({"key": "ke"}, _mat())
    mk = G.agent_card_marks(g, "stop", ["-p", "Bi2Te3", "-j", "S6_elastic", "stop"])
    assert mk["S6_elastic"] == "◀ 本次：取消"
    assert mk["S8_kappa"].startswith("↳")
    mk = G.agent_card_marks(g, "stop", ["-p", "Bi2Te3", "-j", "S7_deform", "stop"])
    assert "S7.1_read" not in mk                   # 已取消的下游不再标
    mk = G.agent_card_marks(g, "stop", ["-p", "Bi2Te3", "stop"])   # 无 -j：有作业的步骤
    assert {k for k, v in mk.items() if v.startswith("◀")} == {"S4_wave", "S6_elastic"}
    mk = G.agent_card_marks(g, "rerun", ["-p", "Bi2Te3", "rerun"])  # 整材料重来
    assert all(v.startswith("◀") for v in mk.values()) and len(mk) == 8
    mk = G.agent_card_marks(g, "rerun", ["-p", "Bi2Te3", "-j", "5", "rerun"])
    assert mk["S5_dielect"] == "◀ 本次：删除重算" and mk["S8_kappa"].startswith("↳")
    assert "整材料重来" in G.agent_consequences("rerun", ["-p", "X", "rerun"])
    assert "-f" in G.agent_consequences("start", ["-p", "X", "start", "-f"])
