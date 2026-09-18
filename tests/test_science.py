import json

from autozt.science import preflight, query_results, research_plan


def test_research_plan_requires_grids_for_zt():
    got = research_plan("比较二维材料 zT", [{"name": "zt-dft-cpu"}], dimension="2D")
    assert got["status"] == "needs_input"
    assert {x["code"] for x in got["gaps"]} == {
        "missing_temperature_grid", "missing_carrier_grid", "missing_material"
    }
    ready = research_plan("比较二维材料 zT", [{"name": "zt-dft-cpu"}], dimension="2D",
                          material="MoS2", temperature=[300], carrier=[1e18])
    assert ready["status"] == "ready"
    assert ready["execution"]["submits_jobs"] is False
    assert {a["action"] for a in ready["actions"]} == {
        "init_project", "start_workflow", "advance_ready"
    }


def test_preflight_reports_cross_workflow_inputs(tmp_path):
    sample = "tmp/_zt_offline/Si/zt-dft-cpu/step20_zt"
    got = preflight(sample, dimension="3D", temperature=[300], carrier=[1e18])
    assert got["status"] == "review"
    assert any(item["code"] == "electronic_file" and not item["ok"] for item in got["checks"])
    bad = preflight(sample, dimension="2D", thickness=6.0, temperature=[12345])
    assert bad["status"] == "review"
    assert any(item["code"] == "temperature_grid" and not item["ok"] for item in bad["checks"])


def test_results_keep_source_and_apply_context_filter(tmp_path):
    (tmp_path / "zt_summary.json").write_text(json.dumps({
        "records": [{"temperature": 300, "carrier": 1e18, "zt": 0.2},
                    {"temperature": 600, "carrier": 1e18, "zt": 0.5}]}), encoding="utf-8")
    got = query_results(str(tmp_path), "zt", temperature=600, carrier=1e18)
    assert got["count"] == 1
    assert got["results"][0]["value"] == 0.5
    assert got["results"][0]["source"].endswith("zt_summary.json")
    assert got["results"][0]["unit"] == "1"
    assert "validator" in got["results"][0]
    assert "provenance" in got["results"][0]
