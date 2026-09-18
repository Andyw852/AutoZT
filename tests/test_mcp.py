"""MCP 接口测试（不需要集群）：工具表、协议握手、风险分级、网关拒绝。

MCP 是 AutoZT 对外给 agent 用的唯一入口，所以这层必须锁死：
  1. 工具是通用动词，数量封顶 21（不许按技能加工具）；
  2. 危险动作必须在协议边界被拒绝并给出批准命令（agent 不能自我批准）；
  3. 只读调用不触发任何提交。
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "autozt")
sys.path.insert(0, ROOT)
from autozt import mcp as M  # noqa: E402


def test_tool_table_is_generic_and_capped():
    assert len(M.TOOLS) >= 24, "当前协议应包含通用工具和科学工作流工具"
    names = [t[0] for t in M.TOOLS]
    assert len(names) == len(set(names))
    # 不许出现"某技能专用"的工具名
    for n in names:
        assert "dft" not in n and "mace" not in n and "phonon" not in n, n


def test_risk_tiers_cover_three_levels():
    risks = {t[3] for t in M.TOOLS}
    assert risks == {"read", "mutate", "destructive"}


def test_every_tool_has_schema():
    for name, desc, schema, risk, verb in M.TOOLS:
        assert schema.get("type") == "object", name
        assert desc and verb, name


def test_initialize_and_tools_list():
    env = dict(os.environ)
    env["AUTOZT_MCP_PROFILE"] = "full"
    p = subprocess.run([sys.executable, PROG, "mcp"], input=
                       json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize"}) + "\n" +
                       json.dumps({"jsonrpc": "2.0", "id": 2, "method": "tools/list"}) + "\n",
                       capture_output=True, text=True, cwd=ROOT, timeout=120, env=env)
    lines = [json.loads(x) for x in p.stdout.splitlines() if x.strip()]
    assert lines[0]["result"]["serverInfo"]["name"] == "autozt"
    assert len(lines[1]["result"]["tools"]) == len(M.TOOLS)


def test_unknown_method_is_an_error():
    resp = M.handle({"jsonrpc": "2.0", "id": 9, "method": "nope"})
    assert resp["error"]["code"] == -32601


def test_destructive_tool_goes_through_gateway():
    # 不连集群：只断言这条路走的是 act 网关（agent 会话下会被拒）
    rc, out, err = M._run(["act", "-tt", "opt-dft-cpu", "-p", "Si_none", "clean", "-y"])
    text = out + err
    assert rc != 0, "破坏性动作在 agent 会话下必须被拒"
    assert "approve" in text or "批准" in text, text[:200]


def test_resources_list_exposes_skill_docs():
    items = M._resources_list()
    uris = [r["uri"] for r in items]
    assert items and any(u.startswith("autozt://skill/") for u in uris)
    assert any(u.endswith("/skill.yaml") for u in uris)


def test_resources_read_returns_file_content():
    uri = "autozt://doc/README.en.md"
    got = M._resource_read(uri)
    assert got and got["contents"][0]["text"].startswith("# AutoZT")


def test_resources_read_refuses_traversal_and_unknown():
    assert M._resource_read("autozt://doc/../../etc/passwd") is None
    assert M._resource_read("autozt://doc/does-not-exist.md") is None
    assert M._resource_read("http://example.com/x") is None


def test_resources_include_configured_skill_paths(tmp_path, monkeypatch):
    skill_root = tmp_path / "skill"
    skill = skill_root / "demo-skill"
    skill.mkdir(parents=True)
    (skill / "skill.yaml").write_text("name: demo-skill\nsteps:\n  - label: S1\n",
                                       encoding="utf-8")
    (skill / "README.md").write_text("# demo\n", encoding="utf-8")
    cfg = tmp_path / "tf.yaml"
    cfg.write_text("skill_paths:\n  - %s\n" % skill_root, encoding="utf-8")
    monkeypatch.setenv("AUTOZT_CONFIG", str(cfg))
    names = [item["name"] for item in M._resources_list()]
    assert "demo-skill/skill.yaml" in names
    assert "demo-skill/README.md" in names
    got = M._resource_read("autozt://skill/demo-skill/README.md")
    assert got["contents"][0]["text"].startswith("# demo")


def test_handle_serves_resource_methods():
    resp = M.handle({"jsonrpc": "2.0", "id": 7, "method": "resources/list"})
    assert resp["result"]["resources"]
    bad = M.handle({"jsonrpc": "2.0", "id": 8, "method": "resources/read",
                    "params": {"uri": "autozt://doc/nope.md"}})
    assert bad["error"]["code"] == -32602


def test_prompts_list_and_get():
    items = M._prompts_list()
    assert {p["name"] for p in items} == {
        "triage-failures", "review-before-submit", "plan-2d-zt",
        "run-validated-workflow", "explain-result-provenance"}
    got = M._prompt_get("review-before-submit", {"material": "Si_demo", "step": "S1_opt"})
    assert "Si_demo" in got["messages"][0]["content"]["text"]
    assert M._prompt_get("nope") is None
    triage = M._prompt_get("triage-failures")
    triage_text = triage["messages"][0]["content"]["text"]
    assert "inspect" in triage_text
    assert "get_summary" not in triage_text and "get_status" not in triage_text
    review_text = got["messages"][0]["content"]["text"]
    assert "inspect" in review_text
    assert "get_status" not in review_text
    resp = M.handle({"jsonrpc": "2.0", "id": 9, "method": "prompts/get",
                     "params": {"name": "nope"}})
    assert resp["error"]["code"] == -32602


def test_readonly_profile_hides_hazardous_tools(monkeypatch=None):
    import os
    os.environ["AUTOZT_MCP_READONLY"] = "1"
    try:
        names = [t["name"] for t in M._tools_list()]
        assert names, "只读档不该为空"
        assert all(n in {x[0] for x in M.TOOLS if x[3] == "read"} for n in names)
        res = M.call_tool("clean_material", {"material": "Si_x"})
        assert res["isError"], "只读档下危险工具必须被拒"
    finally:
        os.environ.pop("AUTOZT_MCP_READONLY", None)
    full = [t["name"] for t in M._tools_list()]
    assert len(full) > len(names), "默认档应比只读档暴露更多工具"


def _fake_status_json(kind="R", diag="", suggested_action="none"):
    return json.dumps({"schema_version": 2, "tf_version": "test", "types": [{
        "key": "demo", "materials": [{"name": "Si_demo", "hpc_name": "jzzn",
        "dim": "3D", "active": {"label": "S1_opt"}, "action": "start S1_opt",
        "steps": [{"label": "S1_opt", "name": "step1", "kind": kind,
                    "diag": diag, "diag_code": "force_not_converged" if diag else "none",
                    "suggested_action": suggested_action,
                    "action_reason": "test reason",
                    "job": {"id": "123", "state": "R", "info": "node01"}}]}]}],
        "queue": {"R": 1, "PD": 0, "total": 1}})


def test_tools_have_output_schema_and_compact_profile():
    old = os.environ.get("AUTOZT_MCP_PROFILE")
    try:
        os.environ["AUTOZT_MCP_PROFILE"] = "compact"
        tools = M._tools_list()
        assert {x["name"] for x in tools} == {
            "schema", "list_skills", "describe_skill", "get_snapshot", "probe_step",
            "propose_actions", "apply_actions"}
        assert all(x["outputSchema"]["type"] == "object" for x in tools)
    finally:
        if old is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old


def test_default_profile_is_workflow_sized():
    old = os.environ.pop("AUTOZT_MCP_PROFILE", None)
    try:
        names = {item["name"] for item in M._tools_list()}
        assert len(names) == 12
        assert "inspect" in names and "cycle" in names
        assert "stop_step" not in names and "rerun_step" not in names
    finally:
        if old is not None:
            os.environ["AUTOZT_MCP_PROFILE"] = old


def test_workflow_profile_exposes_one_shot_control_loop(monkeypatch):
    old = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "workflow"
    raw = _fake_status_json("FAIL", "force not converged", "retry")
    monkeypatch.setattr(M, "_run", lambda argv, timeout=1800: (0, raw, ""))
    try:
        names = [x["name"] for x in M._tools_list()]
        assert "inspect" in names and "cycle" in names
        got = M.call_tool("inspect", {})
        assert got["structuredContent"]["data"]["actions"][0]["action"] == "retry_step"
        policy = got["structuredContent"]["data"]["policy"]
        assert "retry_step" in policy["allowed_actions"]
        assert "retry_step" not in policy["cycle_default_actions"]
        assert policy["retry"]["opt_in"] == "include_retry"
        cycle = M.call_tool("cycle", {"execute": False})
        assert cycle["structuredContent"]["data"]["mode"] == "dry_run"
    finally:
        if old is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old


def test_capabilities_tool_returns_fixed_agent_contract():
    got = M.call_tool("capabilities", {})
    assert not got["isError"]
    data = got["structuredContent"]["data"]
    assert data["schema_version"] == "agent/4"
    assert "request_schema" in M.call_tool("schema", {})["structuredContent"]["data"]
    assert "mcp-stdio" in data["transport"]
    assert data["actions"]["default"] == "dry_run"
    assert "rerun" in data["actions"]["destructive"]


def test_workflow_profile_has_fixed_eight_tools():
    old = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "workflow"
    try:
        names = {item["name"] for item in M._tools_list()}
        assert len(names) == 12
        assert names == {
            "schema", "capabilities", "list_skills", "describe_skill", "get_snapshot",
            "inspect", "probe_step", "cycle", "apply_actions", "research_plan",
            "preflight", "results",
        }
    finally:
        if old is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old


def test_monitor_profile_is_the_smallest_steady_state_surface():
    old = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "monitor"
    try:
        names = {item["name"] for item in M._tools_list()}
        assert names == {"get_snapshot", "cycle"}
    finally:
        if old is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old


def test_high_level_cycle_can_execute_hidden_primitives(monkeypatch):
    old_profile = os.environ.get("AUTOZT_MCP_PROFILE")
    old_run = M._run
    calls = []
    os.environ["AUTOZT_MCP_PROFILE"] = "monitor"
    raw = _fake_status_json("TODO", "", "start")

    def fake_run(argv, timeout=1800):
        calls.append(argv)
        if len(calls) <= 3:
            return 0, raw, ""
        return 0, "started", ""

    M._run = fake_run
    try:
        got = M.call_tool("cycle", {"execute": True})
        assert not got["isError"]
        assert got["structuredContent"]["data"]["mode"] == "execute"
        assert calls[-1][0] == "act"
    finally:
        M._run = old_run
        if old_profile is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old_profile


def test_profile_hides_direct_calls_to_unlisted_tools():
    old = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "workflow"
    try:
        got = M.call_tool("get_status", {"material": "Si_demo"})
        assert got["isError"]
        assert "not exposed" in got["structuredContent"]["error"]
    finally:
        if old is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old


def test_cycle_execute_failure_keeps_structured_execution_data(monkeypatch):
    old_run = M._run
    calls = []
    raw = _fake_status_json("FAIL", "force not converged", "retry")

    def fake_run(argv, timeout=1800):
        calls.append(argv)
        if len(calls) <= 3:
            return 0, raw, ""
        return 1, "", "simulated action failure"

    M._run = fake_run
    try:
        got = M.call_tool("cycle", {"execute": True, "include_retry": True})
        assert got["isError"]
        payload = got["structuredContent"]
        assert payload["returncode"] == 1
        assert payload["data"]["mode"] == "execute"
        assert payload["data"]["execution"]["failed"] is True
        assert "batch stopped" in payload["error"]
        result = payload["data"]["execution"]["results"][0]
        assert result["returncode"] == 1
        assert "simulated action failure" in result["error"]
    finally:
        M._run = old_run


def test_compact_results_do_not_duplicate_large_structured_payload():
    old = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "compact"
    try:
        got = M._result("read", data={"payload": "x" * 5000})
        assert got["structuredContent"]["data"]["payload"] == "x" * 5000
        assert "data_in_structured_content" in got["content"][0]["text"]
        assert len(got["content"][0]["text"]) < 300
    finally:
        if old is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old


def test_partial_skill_contract_gets_explicit_inferred_graph():
    got = M._compact_contract({"skill": "new-skill", "steps": [
        {"label": "S1", "seq": 1, "check": "outcar"}
    ]})
    assert got["contract_status"]["complete"] is False
    assert "io_schema" in got["contract_status"]["missing"]
    assert got["inferred"]["step_graph"][0]["label"] == "S1"


def test_tool_validation_rejects_unknown_and_destructive_batch_actions():
    bad = M.call_tool("describe_skill", {"tt": "demo", "extra": 1})
    assert bad["isError"] and "unknown argument" in bad["structuredContent"]["error"]
    old = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "full"
    try:
        bad = M.call_tool("propose_actions", {"max_actions": True})
        assert bad["isError"] and "must be integer" in bad["structuredContent"]["error"]
    finally:
        if old is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old
    bad = M.call_tool("apply_actions", {"actions": [{"action": "rerun_step"}]})
    assert bad["isError"] and "destructive" in bad["structuredContent"]["error"]


def test_snapshot_returns_cursor_then_only_changes(monkeypatch=None):
    old_run = M._run
    M._SNAPSHOTS.clear()
    state = [_fake_status_json("R"), _fake_status_json("OK")]
    def fake_run(argv, timeout=1800):
        return 0, state.pop(0), ""
    M._run = fake_run
    try:
        first = M.call_tool("get_snapshot", {"view": "attention"})
        assert not first["isError"]
        cursor = first["structuredContent"]["data"]["cursor"]
        second = M.call_tool("get_snapshot", {"cursor": cursor})
        assert not second["isError"]
        payload = second["structuredContent"]["data"]
        assert payload["unchanged"] is False
        step_changes = [x for x in payload["changes"] if x["step"] == "S1_opt"]
        assert step_changes[0]["new"]["kind"] == "OK"
    finally:
        M._run = old_run
        M._SNAPSHOTS.clear()


def test_model_state_collection_uses_read_only_list():
    args = {"tt": "demo", "material": "Si_demo", "status": "error"}
    assert M._state_argv(args) == ["-tt", "demo", "-p", "Si_demo",
                                   "-status", "error", "list", "--json"]


def test_scope_filters_material_and_status_before_model_context():
    raw = {"queue": {}, "types": [
        {"key": "demo", "materials": [
            {"name": "Si_demo", "steps": [{"label": "S1", "kind": "FAIL"}]},
            {"name": "Ge_demo", "steps": [{"label": "S1", "kind": "OK"}]},
        ]},
        {"key": "other", "materials": [
            {"name": "Si_demo", "steps": [{"label": "S1", "kind": "FAIL"}]},
        ]},
    ]}
    got = M._compact_state(raw, {"tt": "demo", "material": "Si_demo",
                                 "status": "error"})
    assert [t["key"] for t in got["types"]] == ["demo"]
    assert [m["name"] for m in got["types"][0]["materials"]] == ["Si_demo"]


def test_probe_routes_non_defect_skills_to_generic_diagnose(monkeypatch=None):
    old_run = M._run
    seen = []
    def fake_run(argv, timeout=1800):
        seen.append(argv)
        return 0, json.dumps({"steps": []}), ""
    M._run = fake_run
    try:
        got = M.call_tool("probe_step", {"tt": "band-dft-cpu", "material": "Si_demo"})
        assert not got["isError"]
        assert seen[-1][-1] == "diagnose"
    finally:
        M._run = old_run


def test_propose_actions_is_advisory(monkeypatch=None):
    old_run = M._run
    old_profile = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "full"
    M._run = lambda argv, timeout=1800: (0, _fake_status_json(
        "FAIL", "force not converged", "retry"), "")
    try:
        got = M.call_tool("propose_actions", {})
        assert not got["isError"]
        proposal = got["structuredContent"]["data"]["proposals"][0]
        assert proposal["tool"] == "retry_step"
        assert proposal["requires_approval"] is False
    finally:
        M._run = old_run
        if old_profile is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old_profile


def test_propose_actions_is_bounded(monkeypatch=None):
    old_run = M._run
    old_profile = os.environ.get("AUTOZT_MCP_PROFILE")
    os.environ["AUTOZT_MCP_PROFILE"] = "full"
    raw = {"queue": {}, "types": [{"key": "demo", "materials": [
        {"name": "Si_%d" % i, "active": {"label": "S1_opt"},
         "steps": [{"label": "S1_opt", "kind": "TODO"}]}
        for i in range(8)
    ]}]}
    M._run = lambda argv, timeout=1800: (0, json.dumps(raw), "")
    try:
        got = M.call_tool("propose_actions", {"max_actions": 3})
        data = got["structuredContent"]["data"]
        assert len(data["proposals"]) == 3
        assert data["proposal_limit"] == 3
    finally:
        M._run = old_run
        if old_profile is None:
            os.environ.pop("AUTOZT_MCP_PROFILE", None)
        else:
            os.environ["AUTOZT_MCP_PROFILE"] = old_profile


def test_apply_actions_dry_run_is_side_effect_free():
    got = M.call_tool("apply_actions", {"dry_run": True, "actions": [
        {"action": "retry_step", "tt": "demo", "material": "Si_demo", "step": "S1_opt"}
    ]})
    assert not got["isError"]
    assert got["structuredContent"]["data"]["dry_run"] is True


def test_apply_actions_rejects_stale_cursor(monkeypatch):
    compact = M._compact_state(json.loads(_fake_status_json("R")))
    monkeypatch.setattr(M, "_mcp_load_compact",
                        lambda _args: (compact, None, 0))
    got = M.call_tool("apply_actions", {
        "cursor": "definitely-old", "tt": "demo",
        "actions": [{"action": "start_step", "material": "Si_demo",
                      "step": "S1_opt"}],
    })
    assert got["isError"]
    assert got["structuredContent"]["data"]["stale_plan"] is True
