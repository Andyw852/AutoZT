"""纯本地 agent CLI 测试：快照持久化、计划校验和稳定输出形状。"""

import argparse
import io
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "autozt")
sys.path.insert(0, ROOT)

from autozt import agent_cli as A  # noqa: E402
from autozt import agentgate as G  # noqa: E402


def _args(tmp_path, **extra):
    base = {"config": str(tmp_path / "tf.yaml"), "tt": None, "material": None,
            "status": None, "view": "attention", "cursor": None}
    base.update(extra)
    return argparse.Namespace(**base)


def _state(kind="R"):
    return {"queue": {"R": 1, "PD": 0}, "types": [{
        "key": "demo", "materials": [{"name": "Si_demo", "hpc_name": "jzzn",
        "active": {"label": "S1_opt"}, "steps": [{"label": "S1_opt", "kind": kind,
        "diag_code": "force_not_converged" if kind == "FAIL" else "none",
        "suggested_action": "retry" if kind == "FAIL" else "start",
        "job": {"id": "1", "state": "R", "info": "node01"}}]}]}]}


def test_plan_rejects_destructive_and_accepts_proposals():
    actions, error = A._normalise_actions({"proposals": [{
        "tool": "prepare_step", "tt": "demo", "material": "Si_demo", "step": "S1_opt"
    }]})
    assert error is None
    assert actions == [{"action": "prepare_step", "tt": "demo", "material": "Si_demo",
                        "step": "S1_opt"}]
    actions, error = A._normalise_actions({"actions": [{"action": "rebuild_step",
                                                         "material": "Si_demo"}]})
    assert actions is None
    assert "destructive" in error


def test_contract_marks_missing_sections_and_infers_graph():
    got = A._compact_contract({"skill": "new", "steps": [{"label": "S1", "seq": 1}]})
    assert got["contract_status"]["complete"] is False
    assert got["inferred"]["step_graph"][0]["label"] == "S1"


def test_science_envelope_promotes_conversation_metadata():
    got = A._envelope("research_plan", {
        "conversation_state": "awaiting_confirmation",
        "message": "计划已生成", "next_action": "confirm_plan",
        "requires_user_confirmation": True,
        "confirmation_payload": {"will_submit_jobs": True},
        "plan_id": "plan-test",
    })
    assert got["conversation"]["conversation_state"] == "awaiting_confirmation"
    assert got["conversation"]["next_action"] == "confirm_plan"
    assert got["data"]["plan_id"] == "plan-test"


def test_snapshot_persists_and_returns_only_changes(tmp_path, monkeypatch):
    state = [_state("R"), _state("OK")]
    monkeypatch.setattr(A, "_status", lambda args: (A._compact_state(state.pop(0)), None, 0))
    args = _args(tmp_path)
    first, rc = A._cmd_snapshot(args)
    assert rc == 0
    assert first["data"]["first_snapshot"] is True
    second, rc = A._cmd_snapshot(args)
    assert rc == 0
    assert second["data"]["unchanged"] is False
    assert any(x["step"] == "S1_opt" for x in second["data"]["changes"])
    assert os.path.isfile(second["data"]["snapshot_file"])


def test_parser_accepts_context_before_agent():
    parser = A.build_parser()
    got = parser.parse_args(["-tt", "demo", "snapshot", "--view", "active"])
    assert got.tt == "demo"
    assert got.agent_command == "snapshot"
    assert got.view == "active"


def test_agent_state_collection_uses_read_only_list(tmp_path):
    args = _args(tmp_path, tt="demo", material="Si_demo", status="error")
    assert A._context_argv(args) + ["list", "--json"] == [
        "-tt", "demo", "-p", "Si_demo", "-status", "error", "list", "--json"]


def test_agent_compact_state_applies_material_scope():
    raw = {"queue": {}, "types": [{"key": "demo", "materials": [
        {"name": "Si_demo", "steps": [{"label": "S1", "kind": "FAIL"}]},
        {"name": "Ge_demo", "steps": [{"label": "S1", "kind": "OK"}]},
    ]}]}
    got = A._compact_state(raw, {"tt": "demo", "material": "Si_demo",
                                 "status": "error"})
    assert [m["name"] for t in got["types"] for m in t["materials"]] == ["Si_demo"]


def test_inspect_and_cycle_are_one_shot_and_dry_by_default(monkeypatch, tmp_path):
    raw = _state("FAIL")
    monkeypatch.setattr(A, "_status", lambda args: (A._compact_state(raw, args), None, 0))
    inspect_args = _args(tmp_path, view="attention", include_monitoring=False, max_actions=20)
    inspected, rc = A._cmd_inspect(inspect_args)
    assert rc == 0
    assert inspected["data"]["actions"][0]["action"] == "prepare_step"
    cycle_args = _args(tmp_path, view="attention", include_monitoring=False,
                       max_actions=20, execute=False, dry_run=True,
                       agent_command="cycle", include_retry=False)
    cycled, rc = A._cmd_cycle(cycle_args)
    assert rc == 0
    assert cycled["data"]["mode"] == "dry_run"
    assert cycled["data"]["execution"]["failed"] is False
    assert cycled["data"]["actions"] == []


def test_propose_uses_shared_bounded_policy(monkeypatch, tmp_path):
    raw = _state("FAIL")
    monkeypatch.setattr(A, "_status", lambda args: (A._compact_state(raw, args), None, 0))
    args = _args(tmp_path, tt="demo", material="Si_demo", max_actions=1,
                 include_monitoring=False)
    result, rc = A._cmd_propose(args)
    assert rc == 0
    assert result["data"]["proposal_limit"] == 1
    assert len(result["data"]["proposals"]) == 1


def test_apply_rejects_plan_when_state_cursor_changed(monkeypatch, tmp_path):
    args = _args(tmp_path)
    current = A._compact_state(_state("R"))
    monkeypatch.setattr(A, "_status", lambda _args: (current, None, 0))
    called = []
    monkeypatch.setattr(A, "_execute_actions",
                        lambda *values: called.append(values) or {"failed": False})
    result, rc = A._apply_payload(args, {
        "cursor": "definitely-old", "scope": {},
        "actions": [{"action": "start_step", "tt": "demo",
                      "material": "Si_demo", "step": "S1_opt"}],
    })
    assert rc == 1
    assert result["data"]["stale_plan"] is True
    assert result["data"]["current_cursor"] == A._protocol.state_cursor(current)
    assert called == []


def test_service_outputs_cursor_for_plan():
    current = A._compact_state(_state("R"))
    data, error, rc = A._service.plan(lambda _scope: (current, None, 0), {})
    assert error is None and rc == 0
    assert data["cursor"] == A._protocol.state_cursor(current)


def test_service_cycle_rechecks_state_before_execute():
    states = [A._compact_state(_state("R")), A._compact_state(_state("OK"))]
    executed = []
    def load(_scope):
        return states.pop(0), None, 0
    data, error, rc = A._service.cycle(
        load, lambda actions, dry_run: executed.append((actions, dry_run)) or {}, {},
        execute=True)
    assert rc == 1
    assert "state changed" in error
    assert data["execution"]["stale_plan"] is True
    assert executed == []


def test_run_ready_steps_uses_one_shot_advance_command(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(A, "_run", lambda argv, config=None: seen.append(argv) or (0, "", ""))
    result = A._execute_actions(
        [{"action": "run_ready_steps", "tt": "demo"}],
        str(tmp_path / "tf.yaml"), False, None, {"tt": "demo"})
    assert result["failed"] is False
    assert seen == [["act", "-tt", "demo", "advance"]]


def test_gateway_policy_lists_advance_and_fails_closed():
    policy = G.render_agent_policy()
    assert "advance" in policy
    assert G.agent_classify("migrate-subdir", ["migrate-subdir"])[0] == "destructive"
    assert "迁移/重排技能目录" in G.agent_classify(
        "migrate-subdir", ["migrate-subdir"])[1][0]
    assert G.agent_classify("push", ["push"])[0] == "destructive"
    assert G.agent_classify("future_command", ["future_command"])[0] == "destructive"


def test_capabilities_and_stdin_request_are_stable(monkeypatch, tmp_path):
    args = _args(tmp_path)
    capabilities, rc = A._cmd_capabilities(args)
    assert rc == 0
    assert capabilities["data"]["schema_version"] == "agent/6"
    assert "propose" in capabilities["data"]["read"]["commands"]
    assert capabilities["data"]["request"]["schema"]["properties"]["op"]["enum"][-1] == "apply"
    assert capabilities["data"]["entrypoints"]["jsonl"] == "autozt agent serve"
    assert "schema" in capabilities["data"]["read"]["commands"]


def test_schema_is_self_describing_and_run_is_exposed(monkeypatch):
    result, rc = A._cmd_schema(argparse.Namespace())
    assert rc == 0
    data = result["data"]
    assert data["schema_version"] == "agent/6"
    assert data["operations"]["run"]["default"] == "dry_run"
    assert "run" in data["request_schema"]["properties"]["op"]["enum"]
    monkeypatch.setattr(sys, "stdin", io.StringIO('{"op":"capabilities"}'))
    request_args = A.build_parser().parse_args(["request", "-"])
    result, rc = A._cmd_request(request_args)
    assert rc == 0
    assert result["command"] == "capabilities"


def test_jsonl_serve_is_one_response_per_input_line(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(sys, "stdin", io.StringIO(
        '{"op":"capabilities"}\nnot-json\n{"op":"capabilities"}\n'))
    args = A.build_parser().parse_args(["serve", "-c", str(tmp_path / "tf.yaml")])
    rc = A._cmd_serve(args)
    lines = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert rc == 2
    assert len(lines) == 3
    assert lines[0]["command"] == "capabilities"
    assert lines[1]["ok"] is False
    assert lines[2]["command"] == "capabilities"


def test_request_rejects_unknown_fields():
    result, rc = A._dispatch_request({"op": "capabilities", "typo": True},
                                     argparse.Namespace(config=None))
    assert rc == 2
    assert result["ok"] is False
    assert "unknown request field" in result["error"]


def test_plain_cli_discovery_works_without_agent_transport():
    for argv in (("--help",), ("skills",), ("schema", "--json")):
        got = subprocess.run([sys.executable, PROG, *argv],
                             capture_output=True, text=True, cwd=ROOT, timeout=60)
        assert got.returncode == 0, (argv, got.stdout, got.stderr)
        assert got.stdout.strip()


def test_evidence_is_read_only_json_wrapper(monkeypatch, tmp_path):
    args = _args(tmp_path, material="Si_demo", tt="demo", step="S1_opt")
    monkeypatch.setattr(A, "_call_json",
                        lambda argv, config: ({"diag_code": "force_not_converged"}, None, 0))
    result, rc = A._cmd_evidence(args)
    assert rc == 0
    assert result["data"]["diag_code"] == "force_not_converged"


def test_run_execute_uses_the_act_gateway(monkeypatch, tmp_path):
    """The no-MCP CLI must retain the same mutation boundary as MCP."""
    raw = _state("TODO")
    calls = []
    monkeypatch.setattr(A, "_status", lambda _args: (A._compact_state(raw), None, 0))
    monkeypatch.setattr(A, "_run",
                        lambda argv, config, timeout=1800: calls.append((argv, config))
                        or (0, "started", ""))
    args = _args(tmp_path, execute=True, dry_run=False, include_retry=False,
                 include_monitoring=False, max_actions=1, view="attention")
    result, rc = A._cmd_run(args)
    assert rc == 0
    assert result["data"]["mode"] == "execute"
    assert result["data"]["execution"]["failed"] is False
    assert calls and calls[0][0][0] == "act"
    assert calls[0][0][-1] == "start"
