"""AI 接入与 2026-09 实测问题的回归测试（纯本地，不连超算）。

覆盖：采集 payload 走 stdin（E2BIG）、调优旋钮进配置、项目级 max_jobs、
被屏蔽项目同名不挡路、稳定材料 id、pip 入口路由 agent/mcp、进度文件、
doctor、默认开启的 agent 网关、结构化失败码、monitor 单轮异常不退出。
"""

import argparse
import contextlib
import io
import json
import os
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import autozt  # noqa: E402
from autozt import agentgate as G  # noqa: E402
from autozt import bootstrap as B  # noqa: E402
from autozt import progress as P  # noqa: E402


@pytest.fixture(autouse=True)
def _clean_conflicts():
    B.reset_config_conflicts()
    yield
    B.reset_config_conflicts()


# ---------------------------------------------------------------- 采集 / E2BIG
def test_collector_payload_goes_through_stdin_not_argv(tmp_path):
    """payload 远超单参数 128 KiB 上限也能采集（以前 --config64 走 argv → E2BIG）。"""
    root = tmp_path / "remote"
    (root / "M0001" / "step1").mkdir(parents=True)
    (root / "M0001" / "step1" / "done.txt").write_text("all DONE here")
    names = ["M%04d_%s" % (i, "x" * 120) for i in range(1500)] + ["M0001"]
    types = [{"key": "demo", "desc": "demo", "root": str(root),
              "steps": [{"name": "step1", "label": "S1", "check": "marker",
                         "marker": "done.txt:DONE", "submit": "submit.sh"}],
              "materials": names}]
    payload_len = len(json.dumps(types))
    assert payload_len * 4 // 3 > 128 * 1024   # base64 后远超单参数上限
    with patch.object(autozt, "skill_checks_for", return_value={}):
        data = autozt.collect({"user": "nobody"}, types, host=None)
    got = {m["name"]: m for m in data["types"][0]["materials"]}
    assert len(got) == len(names)
    assert got["M0001"]["steps"][0]["done"] is True


def test_tuning_knobs_env_over_config_over_default(monkeypatch):
    monkeypatch.delenv("AUTOZT_COLLECT_CHUNK", raising=False)
    assert autozt.tuning_value("collect_chunk", {}, True) == (100, "default")
    assert autozt.tuning_value("collect_chunk", {"collect_chunk": 50}) == 50
    assert autozt.tuning_value("collect_chunk", {"tuning": {"collect_chunk": "30"}}) == 30
    monkeypatch.setenv("AUTOZT_COLLECT_CHUNK", "7")
    assert autozt.tuning_value("collect_chunk", {"collect_chunk": 50}, True) == \
        (7, "env:AUTOZT_COLLECT_CHUNK")
    assert autozt.tuning_value("collect_chunk", {"collect_chunk": "junk"}) == 7


# ---------------------------------------------------------------- max_jobs
def _mat(name, proj_cfg, busy, cap=None):
    steps = [{"kind": "R"} for _ in range(busy)] + [{"kind": "TODO"}]
    return {"name": name, "steps": steps,
            "_seg": {"_from": proj_cfg, "max_jobs": cap}}


def test_project_max_jobs_is_enforced_alongside_skill_cap():
    data = {"types": [
        {"key": "opt", "materials": [_mat("a1", "/p/tf_pilot.yaml", 2, cap=2)]},
        {"key": "opt", "materials": [_mat("b1", "/p/tf_batch.yaml", 3)]},
    ]}
    cfg = {"task_types": {"opt": {"max_jobs": 6}}}
    gate = autozt._SkillGate(cfg, data)
    pilot = data["types"][0]["materials"][0]
    batch = data["types"][1]["materials"][0]
    assert gate.try_acquire("opt", pilot) is False          # 项目级 2 已满
    assert gate.limit_hit("opt", pilot)[0] == "project"
    assert "tf_pilot.yaml" in gate.describe("opt", pilot)
    assert gate.try_acquire("opt", batch) is True           # 别的项目不受影响
    assert gate.busy("opt") == 6
    assert gate.try_acquire("opt", batch) is False          # 技能级 6 已满
    assert gate.limit_hit("opt", batch)[0] == "skill"
    gate.release("opt", batch)
    assert gate.try_acquire("opt", batch) is True


def test_get_types_marks_only_explicit_project_max_jobs():
    cfg = {"task_types": {"opt": {"steps": [{"name": "s1"}], "max_jobs": 100,
                                  "_segments": [
                                      {"local_root": "/a", "_from": "/a/tf_a.yaml", "max_jobs": 16},
                                      {"local_root": "/b", "_from": "/b/tf_b.yaml"}]}}}
    segs = {t["_from"]: t for t in autozt.get_types(cfg) if t.get("_from")}
    assert segs["/a/tf_a.yaml"]["_project_max_jobs"] == 16
    assert segs["/b/tf_b.yaml"]["_project_max_jobs"] is None


# ---------------------------------------------------------------- 屏蔽项目 / 同名
def _project(root, name, skill, mats):
    d = root / name
    (d / "project_setting").mkdir(parents=True)
    (d / "project_setting" / ("tf_%s.yaml" % name)).write_text(
        "task_types:\n  %s:\n    local_root: materials\n" % skill)
    for m in mats:
        (d / "materials" / m).mkdir(parents=True)
        (d / "materials" / m / "POSCAR").write_text("x")
    return d


def _cfg(root):
    cfg = {"project_roots": [str(root)],
           "task_types": {"band": {"steps": [{"name": "s1", "label": "S1"}]}}}
    with contextlib.redirect_stderr(io.StringIO()):
        return B.merge_project_configs(cfg)


def test_stale_blocked_project_does_not_shadow_live_twin(tmp_path):
    _project(tmp_path, "batch", "band", ["Al/X"])
    _project(tmp_path, "pilot_old", "no-such-skill", ["Al/X", "Ca/Q"])
    cfg = _cfg(tmp_path)
    types = autozt.get_types(cfg, quiet=True)
    # basename 与屏蔽项目撞名，但正常项目里也有 → 放行，并解析到正常项目
    B.reject_config_conflict_targets("X", "band", types_fn=lambda: types)
    lp = autozt.resolve_mat_dir(cfg, types, "band", "X")
    assert lp == str((tmp_path / "batch/materials/Al/X").resolve())
    # 只在屏蔽项目里有的材料仍然拒绝，并且说出是哪个项目
    with pytest.raises(SystemExit) as exc:
        B.reject_config_conflict_targets("Q", "band", types_fn=lambda: types)
    assert "pilot_old" in str(exc.value.code) and "已屏蔽" in str(exc.value.code)
    # 明确指向屏蔽项目的写法一律拒绝
    with pytest.raises(SystemExit):
        B.reject_config_conflict_targets("pilot_old/Al/X", "band", types_fn=lambda: types)


def test_block_candidates_are_cached(tmp_path):
    _project(tmp_path, "old", "no-such-skill", ["Al/X"])
    _cfg(tmp_path)
    calls = []
    real = B.discover_local

    def counting(*a, **k):
        calls.append(a)
        return real(*a, **k)
    with patch.object(B, "discover_local", side_effect=counting):
        for _ in range(50):
            B.config_block_hits("X", "band")
    assert len(calls) <= 1


def test_find_material_accepts_qualified_id_and_explains_ambiguity():
    data = {"types": [{"key": "opt", "materials": [
        {"name": "Al/X", "project": "batch", "qualified_name": "batch/Al/X"},
        {"name": "Al/X2", "project": "pilot", "qualified_name": "pilot/Al/X2"},
    ]}]}
    assert autozt.find_material(data, "batch/Al/X")[1]["project"] == "batch"
    data["types"][0]["materials"][1].update(name="Al/X", qualified_name="pilot/Al/X")
    with pytest.raises(SystemExit) as exc:
        autozt.find_material(data, "X")
    msg = str(exc.value.code)
    assert "batch/Al/X" in msg and "pilot/Al/X" in msg and "--project" in msg


def test_skill_local_mats_qualifies_cross_project_duplicates(tmp_path):
    _project(tmp_path, "a", "band", ["Al/X", "Al/Y"])
    _project(tmp_path, "b", "band", ["Al/X"])
    cfg = _cfg(tmp_path)
    names = autozt._skill_local_mats(cfg, autozt.get_types(cfg, quiet=True), "band")
    assert sorted(names) == ["Al/Y", "a/Al/X", "b/Al/X"]


def test_auto_on_dry_run_writes_nothing_and_is_idempotent(tmp_path, capsys):
    _project(tmp_path, "a", "band", ["Al/X"])
    ps = tmp_path / "a/materials/Al/X/project_setting"
    ps.mkdir()
    (ps / "setting.yaml").write_text("auto_advance: false\n")
    cfg = _cfg(tmp_path)
    types = autozt.get_types(cfg, quiet=True)
    assert autozt.cmd_auto_project(cfg, types, "Al/X", "band", "on", dry=True) == 0
    assert (ps / "setting.yaml").read_text() == "auto_advance: false\n"
    assert "[dry-run]" in capsys.readouterr().out
    autozt.cmd_auto_project(cfg, types, "Al/X", "band", "on")
    assert "auto_advance: true" in (ps / "setting.yaml").read_text()
    mtime = (ps / "setting.yaml").stat().st_mtime_ns
    autozt.cmd_auto_project(cfg, types, "Al/X", "band", "on")
    assert (ps / "setting.yaml").stat().st_mtime_ns == mtime   # 已是目标值：不重写


def test_legacy_alias_warning_is_aggregated(capsys):
    for i in range(30):
        B._SKILL_MIGRATIONS.add(("/p%d/tf.yaml" % i, "opt-mace-cpu", "opt-mlff-cpu"))
    B._warn_skill_configs()
    err = capsys.readouterr().err
    assert "opt-mace-cpu→opt-mlff-cpu ×30" in err
    assert len(err.splitlines()) <= 3


# ---------------------------------------------------------------- 入口路由
def test_pip_entrypoint_routes_agent_before_collection(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["autozt", "agent", "capabilities"])
    with patch.object(autozt, "collect_data", side_effect=AssertionError("collected")):
        with pytest.raises(SystemExit) as exc:
            autozt.cli.main()
    assert exc.value.code == 0
    out = json.loads(capsys.readouterr().out)
    assert out["ok"] is True and "progress" in out["data"]["read"]["commands"]


def test_route_ignores_material_named_agent():
    from autozt.cli import route_subcommand
    assert route_subcommand(["-p", "agent", "status"]) is None
    assert route_subcommand(["--project", "mcp", "summary"]) is None


# ---------------------------------------------------------------- 结构化失败码
@pytest.mark.parametrize("diag,code,cls,act", [
    ("S2 段中断（有 .started 无 .done）：VASP 报 ZBRENT: fatal error", "segment_zbrent",
     "NOT_CONVERGED", "retry"),
    ("S1 段中断（有 .started 无 .done）：被 run_relax 看门狗判卡死杀掉", "segment_watchdog",
     "INTERRUPTED", "retry"),
    ("已完成 1/3 段", "segments_incomplete", "INTERRUPTED", "retry"),
    ("3 段全部跑完但未收敛 —— 看 OUTCAR", "segments_exhausted", "NOT_CONVERGED",
     "human_review"),
    ("力已收敛但压强 12.0 kB > 5 kB -> cp", "pressure_not_converged", "NOT_CONVERGED",
     "retry"),
    ("NODE_FAIL on cn12", "node_fail", "NODE_FAIL", "retry"),
    ("imaginary frequency (-0.3)", "imaginary_freq", "PHYSICS", "human_review"),
    ("something odd", "unknown", "UNKNOWN", "human_review"),
])
def test_diag_codes_are_structured(diag, code, cls, act):
    got = autozt._diag_code(diag)
    assert got == code
    assert autozt._diag_class(got) == cls
    assert autozt._suggested_action(got)[0] == act


# ---------------------------------------------------------------- 进度文件
def _data(kind="FAIL", diag="已完成 0/1 段", name="Al/X", proj="batch"):
    return {"queue": {"R": 0}, "types": [{"key": "opt", "materials": [{
        "name": name, "project": proj, "qualified_name": "%s/%s" % (proj, name),
        "hpc_name": "jzzn", "action": "retry S1",
        "active": {"label": "S1", "name": "step1", "kind": kind},
        "steps": [{"label": "S1", "name": "step1", "kind": kind, "diag": diag}]}]}]}


def test_progress_file_roundtrip_merge_and_retry_escalation(tmp_path):
    cfg = {"_config_dir": str(tmp_path), "max_auto_retries": 2}
    P.write_progress(cfg, _data(), full_scope=True)
    P.write_progress(cfg, _data(kind="R", diag="", name="Ca/Z", proj="pilot"))
    doc = P.load_progress(cfg)
    assert {m["id"] for m in doc["materials"]} == {"batch/Al/X", "pilot/Ca/Z"}
    fail = [m for m in doc["materials"] if m["id"] == "batch/Al/X"][0]["fails"][0]
    assert (fail["code"], fail["class"], fail["action"]) == \
        ("segments_incomplete", "INTERRUPTED", "retry")
    # 同一步骤历史上已失败 3 次 → 不再建议自动 retry
    hist = tmp_path / "history.jsonl"
    hist.write_text("".join(json.dumps({"mat": "Al/X", "skill": "opt", "step": "S1",
                                        "ev": "fail", "f": "R", "t": "FAIL"}) + "\n"
                            for _ in range(3)))
    P.write_progress(cfg, _data(), full_scope=False)
    doc = P.load_progress(cfg)
    fail = [m for m in doc["materials"] if m["id"] == "batch/Al/X"][0]["fails"][0]
    assert fail["fail_count"] == 3 and fail["action"] == "human_review"
    view = P.filter_progress(doc, project="pilot")
    assert [m["id"] for m in view["materials"]] == ["pilot/Ca/Z"]
    assert view["projects"] == {"pilot": {"opt": {"materials": 1, "states": {"R": 1}}}}


def test_progress_eta_from_history_medians(tmp_path):
    cfg = {"_config_dir": str(tmp_path)}
    (tmp_path / "history.jsonl").write_text("".join(
        json.dumps({"ev": "finish", "skill": "opt", "step": "S1", "dur": d}) + "\n"
        for d in (1000, 2000, 3000)))
    (tmp_path / ".tf_history_state.json").write_text(json.dumps(
        {"batch/Al/X\topt\tstep1": {"k": "R", "s": "2000-01-01T00:00:00"}}))
    data = _data(kind="R", diag="")
    data["types"][0]["materials"][0]["name"] = "batch/Al/X"
    entry = P.build_entries(cfg, data)["batch/Al/X"]
    assert entry["since"].startswith("2000-01-01") and entry["eta_s"] == 0
    assert "3 finished" in entry["eta_basis"]


def test_agent_progress_command_reads_file_without_subprocess(tmp_path, capsys):
    (tmp_path / "tf.yaml").write_text("project_roots: []\n")
    P.write_progress({"_config_dir": str(tmp_path)}, _data(), full_scope=True)
    from autozt import agent_cli as A
    with patch.object(A, "_run", side_effect=AssertionError("subprocess")):
        rc = A.main(["-c", str(tmp_path / "tf.yaml"), "progress", "--project", "batch"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["ok"] and out["data"]["materials"][0]["id"] == "batch/Al/X"
    assert "liveness" in out["data"]


def test_mcp_get_progress_tool(tmp_path, monkeypatch):
    (tmp_path / "tf.yaml").write_text("project_roots: []\n")
    P.write_progress({"_config_dir": str(tmp_path)}, _data(), full_scope=True)
    monkeypatch.setenv("AUTOZT_CONFIG", str(tmp_path / "tf.yaml"))
    from autozt import mcp as M
    got = M.call_tool("get_progress", {"status": "error"})
    assert not got["isError"]
    assert got["structuredContent"]["data"]["materials"][0]["fails"][0]["class"] == "INTERRUPTED"


# ---------------------------------------------------------------- 网关默认开启
def test_gate_blocks_non_tty_destructive_without_actor(tmp_path, monkeypatch):
    for k in ("AUTOZT_ACTOR", "AUTOZT_AGENT_GATE", "AUTOZT_AGENT_STRICT",
              "AUTOZT_AGENT_GATEWAY") + G.AGENT_ENV_MARKERS:
        monkeypatch.delenv(k, raising=False)
    cfg = {"_config_dir": str(tmp_path)}
    monkeypatch.setattr(G, "_stdin_isatty", lambda: False)
    with contextlib.redirect_stdout(io.StringIO()):
        assert G.agent_direct_gate(cfg, "clean", ["-p", "X", "clean", "-y"]) == 3
    assert G.agent_direct_gate(cfg, "summary", ["summary"]) is None     # 只读不拦
    monkeypatch.setattr(G, "_stdin_isatty", lambda: True)
    assert G.agent_direct_gate(cfg, "clean", ["-p", "X", "clean", "-y"]) is None  # 人在终端
    monkeypatch.setattr(G, "_stdin_isatty", lambda: False)
    monkeypatch.setenv("AUTOZT_AGENT_GATE", "off")
    assert G.agent_direct_gate(cfg, "clean", ["-p", "X", "clean", "-y"]) is None


def test_gate_detects_known_and_configured_agent_markers(monkeypatch):
    for k in ("AUTOZT_ACTOR", "AUTOZT_AGENT_GATE") + G.AGENT_ENV_MARKERS:
        monkeypatch.delenv(k, raising=False)
    assert G.agent_detect({}) == ("", None)
    monkeypatch.setenv("DSH_SESSION", "1")
    assert G.agent_detect({"agent_env_markers": ["DSH_SESSION"]})[0] == "agent:dsh_session"
    monkeypatch.setenv("AUTOZT_ACTOR", "bot")
    assert G.agent_detect({}) == ("bot", "env:AUTOZT_ACTOR")
    assert G.agent_classify("progress", [])[0] == "read"
    assert G.agent_classify("doctor", [])[0] == "read"


# ---------------------------------------------------------------- doctor
def test_doctor_reports_project_caps_knob_sources_and_shadowing(tmp_path, monkeypatch):
    monkeypatch.setenv("AUTOZT_COLLECT_CHUNK", "50")
    d = _project(tmp_path, "pilot", "band", ["Al/X"])
    (d / "project_setting/tf_pilot.yaml").write_text(
        "task_types:\n  band:\n    local_root: materials\n    max_jobs: 16\n")
    _project(tmp_path, "old", "no-such-skill", ["Al/X"])
    cfg = _cfg(tmp_path)
    cfg["task_types"]["band"]["max_jobs"] = 100
    cfg["_config_dir"] = str(tmp_path)
    rep = autozt.run_doctor(cfg, autozt.get_types(cfg, quiet=True))
    row = [r for r in rep["max_jobs"] if r["project"] == "pilot"][0]
    assert (row["skill_cap"], row["project_cap"], row["effective_cap"]) == (100, 16, 16)
    codes = {f["code"] for f in rep["findings"]}
    assert {"knob_env_only", "project_blocked", "blocked_project_shadow"} <= codes


# ---------------------------------------------------------------- monitor
def test_monitor_round_error_does_not_kill_loop(tmp_path, capsys):
    cfg = {"_config_dir": str(tmp_path)}
    calls = {"n": 0}

    def flaky_collect(c, t):
        calls["n"] += 1
        if calls["n"] == 1:
            raise SystemExit("错误：连接超时（host=x）。")
        return _data()

    sleeps = {"n": 0}

    def fake_sleep(_s):
        sleeps["n"] += 1
        if sleeps["n"] >= 2:
            raise KeyboardInterrupt

    with patch.object(autozt, "collect_data", side_effect=flaky_collect), \
            patch.object(autozt, "fill_local_dim"), \
            patch.object(autozt, "auto_fetch"), patch.object(autozt, "auto_advance"), \
            patch.object(autozt.ops, "auto_recover_hung"), \
            patch.object(autozt.ops, "_watch_cfg_sig", return_value=()), \
            patch("time.sleep", side_effect=fake_sleep):
        autozt.cmd_watch(cfg, [], [], None, 1)
    out = capsys.readouterr().out
    assert "第 1 轮异常" in out and calls["n"] == 2
    rounds = [json.loads(l[len("[round] "):]) for l in out.splitlines()
              if l.startswith("[round] ")]
    assert rounds[0]["error"].startswith("错误：连接超时")
    assert rounds[1]["round"] == 2 and rounds[1]["steps"] == {"FAIL": 1}
    assert P.load_progress(cfg)["writer"] == "monitor"


def test_monitor_daemon_check_skips_project_scan_when_running(tmp_path, monkeypatch, capsys):
    cfg_path = tmp_path / "tf.yaml"
    cfg_path.write_text("project_roots: []\n")
    monkeypatch.setattr(sys, "argv", ["autozt", "-c", str(cfg_path), "monitor", "-d"])
    with patch.object(autozt, "_watch_running_pid", return_value=(4242, "x")), \
            patch.object(autozt, "merge_project_configs",
                         side_effect=AssertionError("scanned")), \
            patch.object(autozt, "agent_direct_gate", return_value=None):
        autozt.cli.main()
    assert "PID 4242" in capsys.readouterr().out


# ---------------------------------------------------------------- 技能库命令不扫项目
@pytest.mark.parametrize("argv", [["skills"], ["schema", "--json"],
                                  ["schema", "band-dft-cpu", "--json"],
                                  ["skill"], ["history"]])
def test_skill_library_commands_skip_project_scan(tmp_path, monkeypatch, capsys, argv):
    """skills/schema/skill（及不带目标的 history）只读技能库，不能触发 merge_project_configs
    ——9p 上它每次 1–2 分钟，agent skills/contract 和 MCP list_skills/describe_skill
    都经 `schema --json`，以前每调一次都白扫一遍。"""
    cfg_path = tmp_path / "tf.yaml"
    cfg_path.write_text("project_roots: [%s]\n" % tmp_path)
    monkeypatch.setattr(sys, "argv", ["autozt", "-c", str(cfg_path)] + argv)
    with patch.object(autozt, "merge_project_configs",
                      side_effect=AssertionError("project scan")), \
            patch.object(autozt, "agent_direct_gate", return_value=None):
        try:
            rc = autozt.cli.main()
        except SystemExit as exc:
            rc = exc.code
    assert rc in (None, 0)
    out = capsys.readouterr().out
    if "--json" in argv and argv[0] == "schema":
        skills = json.loads(out)["skills"]
        assert skills and all("skill" in s for s in skills)
        if len(argv) == 3:
            assert [s["skill"] for s in skills] == ["band-dft-cpu"]


def test_history_with_target_still_checks_blocked_projects(tmp_path, monkeypatch):
    cfg_path = tmp_path / "tf.yaml"
    cfg_path.write_text("project_roots: [%s]\n" % tmp_path)
    monkeypatch.setattr(sys, "argv", ["autozt", "-c", str(cfg_path), "-p", "X", "history"])
    with patch.object(autozt, "merge_project_configs",
                      side_effect=AssertionError("project scan")), \
            patch.object(autozt, "agent_direct_gate", return_value=None):
        with pytest.raises(AssertionError, match="project scan"):
            autozt.cli.main()


def test_skill_library_commands_with_explicit_target_still_check_blocks(tmp_path, monkeypatch):
    cfg_path = tmp_path / "tf.yaml"
    cfg_path.write_text("project_roots: [%s]\n" % tmp_path)
    for argv in (["-p", "X", "schema", "--json"], ["-p", "X", "skills"]):
        monkeypatch.setattr(sys, "argv", ["autozt", "-c", str(cfg_path)] + argv)
        with patch.object(autozt, "merge_project_configs",
                          side_effect=AssertionError("project scan")), \
                patch.object(autozt, "agent_direct_gate", return_value=None):
            with pytest.raises(AssertionError, match="project scan"):
                autozt.cli.main()


# ---------------------------------------------------------------- history 状态按范围合并
def _hist_data(mat, tt, kind):
    return {"types": [{"key": tt, "materials": [
        {"name": mat, "steps": [{"name": "step1", "label": "S1", "kind": kind, "diag": ""}]}]}]}


def test_history_state_merges_across_scopes(tmp_path):
    """两个不同范围的采集者（如 -tt 限定的 monitor + 另一个项目的 monitor / 一次 -p 查询）
    共用配置目录时，互相不能抹掉对方的基线，否则状态转移（和 FAIL 次数）会丢。"""
    from autozt.history import history_record, history_load
    cfg = {"_config_dir": str(tmp_path)}
    history_record(cfg, _hist_data("A", "opt", "R"))       # 首轮：只落基线
    history_record(cfg, _hist_data("B", "band", "R"))      # 另一个范围
    history_record(cfg, _hist_data("A", "opt", "FAIL"))    # A 的转移必须被记下
    history_record(cfg, _hist_data("B", "band", "OK"))     # B 的转移也必须被记下
    evs, _ = history_load(cfg)
    got = sorted((e["mat"], e["f"], e["t"]) for e in evs)
    assert got == [("A", "R", "FAIL"), ("B", "R", "OK")]
