import json

from autozt.science import preflight, query_results, research_plan


def test_research_plan_requires_grids_for_zt():
    got = research_plan("比较二维材料 zT", [{"name": "zt-dft-cpu"}], dimension="2D")
    assert got["status"] == "needs_input"
    assert got["conversation_state"] == "needs_user_input"
    assert got["next_action"] == "provide_inputs"
    assert {x["code"] for x in got["gaps"]} == {
        "missing_temperature_grid", "missing_carrier_grid", "missing_material"
    }
    ready = research_plan("比较二维材料 zT", [{"name": "zt-dft-cpu"}], dimension="2D",
                          material="MoS2", temperature=[300], carrier=[1e18])
    assert ready["status"] == "ready"
    assert ready["conversation_state"] == "awaiting_confirmation"
    assert ready["requires_user_confirmation"] is True
    assert ready["confirmation_payload"]["will_submit_jobs"] is True
    assert ready["review_card"]["proceed"]["label"] == "Proceed"
    assert ready["review_card"]["proceed"]["enabled"] is True
    assert ready["review_card"]["execution_boundary"]["this_plan_submits_jobs"] is False
    assert ready["action_surface"]["mcp_safe_actions"] == [ready["actions"][-1]]
    assert len(ready["action_surface"]["cli_only_actions"]) == 2
    assert ready["execution"]["submits_jobs"] is False
    assert {a["action"] for a in ready["actions"]} == {
        "init_project", "start_workflow", "run_ready_steps"
    }
    assert ready["actions"][-1]["command"][-1] == "advance"


def test_preflight_reports_cross_workflow_inputs(tmp_path):
    sample = "tmp/_zt_offline/Si/zt-dft-cpu/step20_zt"
    got = preflight(sample, dimension="3D", temperature=[300], carrier=[1e18])
    assert got["status"] == "review"
    assert got["conversation_state"] == "preflight_review"
    assert any(item["code"] == "electronic_file" and not item["ok"] for item in got["checks"])
    bad = preflight(sample, dimension="2D", thickness=6.0, temperature=[12345])
    assert bad["status"] == "review"
    assert bad["next_action"] == "review_preflight"
    assert any(item["code"] == "temperature_grid" and not item["ok"] for item in bad["checks"])


def test_results_keep_source_and_apply_context_filter(tmp_path):
    (tmp_path / "zt_summary.json").write_text(json.dumps({
        "records": [{"temperature": 300, "carrier": 1e18, "zt": 0.2},
                    {"temperature": 600, "carrier": 1e18, "zt": 0.5}]}), encoding="utf-8")
    got = query_results(str(tmp_path), "zt", temperature=600, carrier=1e18)
    assert got["count"] == 1
    assert got["conversation_state"] == "completed"
    assert got["results"][0]["value"] == 0.5
    assert got["results"][0]["source"].endswith("zt_summary.json")
    assert got["results"][0]["unit"] == "1"
    assert "validator" in got["results"][0]
    assert "provenance" in got["results"][0]


def test_research_plan_selects_non_zt_skills():
    skills = [{"name": name} for name in (
        "opt-dft-cpu", "phonon-dft-cpu", "kl-dft-cpu", "elastic-dft-cpu",
        "defect-dft-cpu", "band-dft-cpu")]
    cases = {
        "帮我优化Si结构": "opt-dft-cpu",
        "计算Si声子谱": "phonon-dft-cpu",
        "计算Si晶格热导率": "kl-dft-cpu",
        "计算Si缺陷形成能": "defect-dft-cpu",
        "计算Si弹性常数": "elastic-dft-cpu",
        "计算Si能带": "band-dft-cpu",
    }
    for goal, expected in cases.items():
        got = research_plan(goal, skills, material="Si", dimension="3D")
        assert got["status"] == "ready", (goal, got)
        assert got["selected_skills"] == [expected], (goal, got)
        assert got["actions"][1]["command"] == [
            "autozt", "-tt", expected, "-p", "Si", "start"
        ]


def test_research_plan_does_not_claim_unknown_goal_is_ready():
    got = research_plan("帮我做一个未知的计算", [{"name": "opt-dft-cpu"}], material="Si")
    assert got["status"] == "needs_input"
    assert got["next_action"] == "provide_inputs"
    assert "no_matching_skill" in {item["code"] for item in got["gaps"]}
    assert got["actions"] == []
