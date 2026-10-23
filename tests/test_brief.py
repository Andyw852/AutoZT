# -*- coding: utf-8 -*-
"""材料基础信息卡（autozt brief / agent brief / MCP material_brief）——纯本地，不连超算。"""
import json
import os
import subprocess
import sys
import time

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from autozt import brief as BR  # noqa: E402

_AGENT_ENV = ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT", "GEMINI_CLI", "CODEX_SANDBOX",
              "AUTOZT_ACTOR", "AUTOZT_JSON_ERRORS")
MID = "Si1_fit-fc-thermal/Si1"


def _env(**kw):
    env = {k: v for k, v in os.environ.items() if k not in _AGENT_ENV}
    env.update({k: str(v) for k, v in kw.items()})
    return env


@pytest.fixture()
def sb(tmp_path):
    proj = tmp_path / "proj"
    proj.mkdir()
    cfg = tmp_path / "tf.yaml"
    cfg.write_text('host: ""\nproject_roots:\n  - %s\ntask_types:\n'
                   '  fit-fc-thermal: {hpc: {name: local, ssh_host: ""}}\n' % proj)
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "POSCAR").write_text("Si\n1.0\n0 2.7 2.7\n2.7 0 2.7\n2.7 2.7 0\nSi\n2\n"
                               "Direct\n0 0 0\n0.25 0.25 0.25\n")
    box = {"cfg": str(cfg), "tmp": tmp_path, "ds": ds}
    r = _cli(box, "-tt", "fit-fc-thermal", "-p", "Si1", "register", "--dataset", str(ds),
             "--cluster", "local")
    assert r.returncode == 0, r.stdout + r.stderr
    r = _cli(box, "-tt", "fit-fc-thermal", "-p", "Si1", "-j", "step1_fit", "conf",
             "--set", "params.FIT_ENGINE=phono3py")
    assert r.returncode == 0, r.stdout + r.stderr
    return box


def _cli(sb, *args, **env):
    return subprocess.run([sys.executable, "-m", "autozt", "-c", sb["cfg"]] + list(args),
                          capture_output=True, text=True, cwd=ROOT,
                          env=_env(HOME=sb["tmp"], **env), timeout=300)


def _fake_job(sb, state="PD", jid="28903807", info="(Priority)"):
    now = time.time()
    entry = {"id": MID, "name": "Si1", "tt": "fit-fc-thermal", "hpc": "local",
             "state": state, "active": "S1_fit", "done_steps": 0, "total_steps": 3,
             "obs": now - 60, "steps": {"S1_fit": state, "S2_kappa": "WAIT"}}
    if state in ("PD", "R"):
        entry["job"] = {"id": jid, "state": state, "info": info}
    if state == "FAIL":
        entry["fails"] = [{"step": "S1_fit", "code": "scf_not_converged", "action": "retry"}]
    (sb["tmp"] / ".tf_progress.json").write_text(json.dumps(
        {"schema": "autozt-progress/1", "ts": now - 60, "materials": [entry]}))
    (sb["tmp"] / ".tf_submit_ledger.jsonl").write_text(json.dumps(
        {"ts": now - 600, "tt": "fit-fc-thermal", "id": MID, "step": "S1_fit",
         "jobs": 1, "jid": jid}) + "\n")


def test_brief_json_reports_cluster_workdir_params_and_job(sb):
    _fake_job(sb)
    t0 = time.time()
    r = _cli(sb, "-p", "Si1", "brief", "--json")
    assert r.returncode == 0, r.stdout + r.stderr
    assert time.time() - t0 < 60                       # 只读配置 + 本地文件，不采集
    b = json.loads(r.stdout)
    assert b["id"] == MID and b["tt"] == "fit-fc-thermal"
    assert b["hpc"]["name"] == "local" and b["hpc"]["ssh_host"] == "本机"
    assert b["hpc"]["source"].endswith(os.path.join("project_setting", "hpc.yaml"))
    assert b["work_dir"]["path"] and b["work_dir"]["source"]
    params = {p["key"]: p for p in b["params"]}
    assert params["FIT_ENGINE"]["value"] == "phono3py"
    assert params["FIT_ENGINE"]["default"] == "pheasy"     # 技能出厂默认
    assert params["FIT_ENGINE"]["steps"] == ["S1_fit"]
    assert params["FIT_INPUT_DIR"]["value"] == str(sb["ds"])
    # init 原样拷来的默认值不算"改过"
    assert set(params) == {"FIT_ENGINE", "FIT_INPUT_DIR"}, sorted(params)
    job = b["job"]
    assert job["state"] == "PD" and job["job"]["id"] == "28903807"
    assert job["last_submit"]["jid"] == "28903807"
    assert b["auto_advance"] == {"on": False, "source": os.path.join(
        "fit-fc-thermal", "project_setting", "setting.yaml")}
    md = b["brief_md"].splitlines()
    assert md[0] == "**【AutoZT 基础信息】%s（fit-fc-thermal）**" % MID
    assert md[1].startswith("1. hpc → local ✅")
    assert md[2].startswith("2. work_dir → `")
    assert "FIT_ENGINE=phono3py（S1_fit）" in md[3]
    assert md[4].startswith("4. 作业 → S1_fit 已提交（28903807，") and "排队中（Priority）" in md[4]
    assert md[5].startswith("5. 自动推进 → 关")


def test_brief_text_and_step_scope(sb):
    _fake_job(sb, state="FAIL")
    r = _cli(sb, "-tt", "fit-fc-thermal", "-p", "Si1", "brief")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "❌ 失败：scf_not_converged，建议 retry" in r.stdout
    r = _cli(sb, "-p", "Si1", "-j", "2", "brief")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "3. 参数 → 本步全部为技能默认值" in r.stdout
    assert "4. 作业 → S2_kappa · 等上游" in r.stdout


def test_brief_without_progress_points_to_refresh(sb):
    r = _cli(sb, "-p", "Si1", "brief")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "暂无作业记录" in r.stdout and "list --refresh" in r.stdout


def test_brief_needs_material_and_rejects_unknown(sb):
    assert _cli(sb, "brief").returncode != 0
    assert _cli(sb, "-p", "Nope", "brief").returncode != 0


def test_agent_brief_envelope_and_request(sb):
    _fake_job(sb, state="R", info="cu12")
    r = _cli(sb, "agent", "brief", "-p", "Si1")
    assert r.returncode == 0, r.stdout + r.stderr
    env = json.loads(r.stdout)
    assert env["ok"] and env["command"] == "brief"
    assert "▶ 运行中（@cu12）" in env["data"]["brief_md"]
    req = subprocess.run([sys.executable, "-m", "autozt", "-c", sb["cfg"], "agent",
                          "request", "-"], input=json.dumps({"op": "brief", "material": "Si1",
                                                             "step": "S1_fit"}),
                         capture_output=True, text=True, cwd=ROOT,
                         env=_env(HOME=sb["tmp"]), timeout=300)
    env = json.loads(req.stdout)
    assert env["ok"] and env["data"]["job"]["step_state"] == "R", req.stdout + req.stderr
    bad = json.loads(_cli(sb, "agent", "brief").stdout)
    assert bad["ok"] is False and "-p" in bad["error"]


def test_mcp_material_brief_text_is_the_card(sb, monkeypatch):
    from autozt import mcp as M
    _fake_job(sb)
    monkeypatch.setenv("AUTOZT_CONFIG", sb["cfg"])
    monkeypatch.setenv("HOME", str(sb["tmp"]))
    monkeypatch.delenv("AUTOZT_MCP_PROFILE", raising=False)
    for k in _AGENT_ENV:
        monkeypatch.delenv(k, raising=False)
    assert "material_brief" in M.CORE_TOOLS and M.STAGE_OF["material_brief"] == "report"
    got = M.call_tool("material_brief", {"material": "Si1"})
    assert not got["isError"], got
    card = got["structuredContent"]["data"]["brief_md"]
    assert got["content"][0]["text"] == card and "FIT_ENGINE=phono3py" in card
    assert "material_brief" in M.HARNESS_RULES


def test_render_md_done_and_override_summary():
    b = {"id": "p/X", "tt": "kl-dft-cpu",
         "hpc": {"name": "jzzn", "source": "project_setting/hpc.yaml", "ssh_host": "jzzn",
                 "step_overrides": [{"step": "S7_deform", "hpc": "a800"}]},
         "work_dir": {"path": "/w", "source": "project_setting/setting.yaml"},
         "params": [{"section": "submit", "key": "ntasks_per_node", "value": "24",
                     "steps": ["S1", "S2", "S3"]}],
         "job": {"progress": True, "state": "DONE", "done_steps": 3, "total_steps": 3},
         "auto_advance": {"on": True, "source": "全局 tf.yaml"}, "_nsteps": 3}
    md = BR.render_md(b)
    assert "分步：S7_deform@a800" in md
    assert "submit.ntasks_per_node=24（全部步骤）" in md
    assert "4. 作业 → ✅ 全部 3/3 步完成" in md
    assert "5. 自动推进 → 开（全局 tf.yaml）" in md
