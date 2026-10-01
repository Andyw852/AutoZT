# -*- coding: utf-8 -*-
"""agent 快路径：进度文件读状态、-p 只采目标材料、max_jobs 基线/账本、分组执行 + 前置条件。

端到端用例在 tmp 里搭一个本地模式小集群：假 sbatch/squeue（队列就是一个文本文件），
opt-dft-cpu 技能，项目 big（项目级 max_jobs=3）8 个材料，其中 1 个已在跑。
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from autozt import agent_state as S  # noqa: E402
from autozt import agentgate as G  # noqa: E402
from autozt import progress as P  # noqa: E402
from autozt import workflow as W  # noqa: E402

PROG = str(ROOT / "bin" / "autozt")
N_MATS = 8


@pytest.fixture
def farm(tmp_path):
    root = tmp_path
    (root / "setting").mkdir()
    (root / "setting" / "tf.yaml").write_text(
        "project_roots: [%s]\nhang_check: false\ntask_types:\n  opt-dft-cpu:\n"
        "    max_jobs: 100\n"
        # 内联本机集群：否则回退技能默认 hpc=jzzn，开发机上真实的（gitignore 的）
        # setting/jzzn.yaml 会让测试去 ssh jzzn
        "    hpc: {name: farm, ssh_host: \"\"}\n" % (root / "projects"), encoding="utf-8")
    ps = root / "projects" / "big" / "project_setting"
    ps.mkdir(parents=True)
    (ps / "tf_big.yaml").write_text(
        "task_types:\n  opt-dft-cpu:\n    local_root: materials\n    work_dir: %s\n"
        "    max_jobs: 3\n" % (root / "remote"), encoding="utf-8")
    for i in range(N_MATS):
        m = root / "projects" / "big" / "materials" / ("E%d" % (i % 2)) / ("m%d" % i)
        (m / "project_setting").mkdir(parents=True)
        (m / "opt-dft-cpu" / "result").mkdir(parents=True)
        (m / "POSCAR").write_text("m\n1.0\n1 0 0\n0 1 0\n0 0 1\nSi\n1\nDirect\n0 0 0\n")
        (m / "project_setting" / "setting.yaml").write_text("auto_advance: false\n")
        sd = root / "remote" / ("E%d" % (i % 2)) / ("m%d" % i) / "opt-dft-cpu" / "step1_opt"
        sd.mkdir(parents=True)
        for f in ("INCAR", "POSCAR", "KPOINTS", "POTCAR"):
            (sd / f).write_text("x\n")
        (sd / "submit.sh").write_text("#!/bin/bash\necho run\n")
    queue = root / "queue.txt"
    running = root / "remote" / "E0" / "m0" / "opt-dft-cpu" / "step1_opt"
    queue.write_text("900|j|RUNNING|1:00|1|node1|%s\n" % running)
    bindir = root / "bin"
    bindir.mkdir()
    (bindir / "squeue").write_text("#!/bin/bash\ncat %s\n" % queue)
    (bindir / "sbatch").write_text(
        "#!/bin/bash\nn=$(( $(wc -l < %s) + 1000 ))\n"
        "echo \"$n|job|PENDING|0:00|1|(Priority)|$PWD\" >> %s\n"
        "echo \"Submitted batch job $n\"\n" % (queue, queue))
    for f in ("squeue", "sbatch"):
        os.chmod(bindir / f, 0o755)
    env = dict(os.environ)
    env["PATH"] = "%s:%s" % (bindir, env.get("PATH", ""))
    env["AUTOZT_AGENT_EVENTS"] = "1"
    env.pop("AUTOZT_NARROW_MAX_AGE", None)
    env["AUTOZT_AGENT_SOURCE"] = "auto"
    env["AUTOZT_SKIP_ENV_CHECK"] = "1"   # 夹具没有 VASP；本机依赖闸门不是这里要测的

    def run(*argv, **extra):
        e = dict(env)
        e.update(extra)
        return subprocess.run([sys.executable, PROG, "-c", str(root / "setting" / "tf.yaml")]
                              + list(argv), capture_output=True, text=True, env=e,
                              timeout=300, cwd=str(root))

    def events(out):
        return S.parse_events(out)

    return argparse.Namespace(root=root, run=run, events=events, queue=queue,
                              cfg={"_config_dir": str(root / "setting")})


def _ids(n):
    return ",".join("big/E%d/m%d" % (i % 2, i) for i in n)


def test_full_collect_records_coverage_and_busy(farm):
    r = farm.run("-tt", "opt-dft-cpu", "list", "--refresh", "--json")
    assert r.returncode == 0, r.stderr[-2000:]
    doc = P.load_progress(farm.cfg)
    assert "opt-dft-cpu" in doc["coverage"]
    busy = {e["id"]: e["busy"] for e in doc["materials"]}
    assert busy["big/E0/m0"] == 1 and sum(busy.values()) == 1
    assert all("obs" in e for e in doc["materials"])


def test_narrowed_start_counts_uncollected_busy_and_ledger(farm):
    assert farm.run("-tt", "opt-dft-cpu", "list", "--refresh", "--json").returncode == 0
    # -p 只采目标材料：list 只返回 1 个材料
    r = farm.run("-tt", "opt-dft-cpu", "-p", "big/E1/m5", "list", "--json")
    assert sum(len(t["materials"]) for t in json.loads(r.stdout)["types"]) == 1
    # 项目上限 3，m0 在跑（不在采集范围里，靠进度文件基线）→ 只能再交 2 个
    r = farm.run("-tt", "opt-dft-cpu", "-p", _ids([2, 3, 4]), "-j", "S1_opt", "start")
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    kinds = [ev["event"] for ev in farm.events(r.stdout)]
    assert kinds.count("submitted") == 2 and kinds.count("deferred") == 1
    # 第二次：刚交的 2 个在进度文件里还是 TODO（采集在提交之前），账本把它们补上 → 满
    r = farm.run("-tt", "opt-dft-cpu", "-p", _ids([6]), "-j", "S1_opt", "start")
    assert [ev["event"] for ev in farm.events(r.stdout)] == ["deferred"]
    led = P.ledger_recent(farm.cfg, 0)
    assert len(led) == 2 and all(x["tt"] == "opt-dft-cpu" for x in led)


def test_narrowing_off_without_fresh_full_collection(farm):
    # 没有整技能采集记录（只有 -p 采集）→ 全量采集，JSON 里是全部材料
    r = farm.run("-tt", "opt-dft-cpu", "-p", "big/E1/m5", "list", "--json")
    assert sum(len(t["materials"]) for t in json.loads(r.stdout)["types"]) == N_MATS
    farm.run("-tt", "opt-dft-cpu", "list", "--refresh", "--json")
    r = farm.run("-tt", "opt-dft-cpu", "-p", "big/E1/m5", "list", "--json",
                 AUTOZT_NARROW_MAX_AGE="0")
    assert sum(len(t["materials"]) for t in json.loads(r.stdout)["types"]) == N_MATS


def test_explicit_start_honours_cap_and_shared_gate_without_narrowing(farm):
    # 无进度文件（全量采集）时，显式 -p A,B,C -j 也受项目上限约束（以前完全不过闸门）
    r = farm.run("-tt", "opt-dft-cpu", "-p", _ids([1, 2, 3, 4]), "-j", "S1_opt", "start")
    kinds = [ev["event"] for ev in farm.events(r.stdout)]
    assert kinds.count("submitted") == 2 and kinds.count("deferred") == 2


def test_plan_from_progress_then_apply_skips_changed_material(farm):
    farm.run("-tt", "opt-dft-cpu", "list", "--refresh", "--json")
    r = farm.run("agent", "plan", "-tt", "opt-dft-cpu", "--project", "big",
                 "--max-actions", "2")
    plan = json.loads(r.stdout)
    assert plan["data"]["source"]["kind"] == "progress"
    acts = plan["data"]["actions"]
    assert len(acts) == 2 and all(a["action"] == "start_step" for a in acts)
    # 计划之后、执行之前，第二个材料被别人提交了
    tgt = acts[1]["material"].split("/", 1)[1]
    with open(farm.queue, "a") as f:
        f.write("77|j|RUNNING|1:00|1|n|%s\n"
                % (farm.root / "remote" / tgt / "opt-dft-cpu" / "step1_opt"))
    planfile = farm.root / "plan.json"
    planfile.write_text(r.stdout)
    r = farm.run("agent", "apply", str(planfile))
    out = json.loads(r.stdout)
    assert out["ok"], out
    assert len(out["data"]["groups"]) == 1
    assert "--expect-state=TODO,PREP,SCANCEL" in out["data"]["groups"][0]["argv"]
    assert [x["outcome"] for x in out["data"]["results"]] == ["submitted", "skipped_stale"]
    # 读视图把刚提交的叠加成 PD，下一轮 inspect 不再提议它
    r = farm.run("agent", "inspect", "-tt", "opt-dft-cpu", "--project", "big")
    again = [a["material"] for a in json.loads(r.stdout)["data"]["actions"]]
    assert acts[0]["material"] not in again


def test_expect_state_requires_targets(farm):
    r = farm.run("-tt", "opt-dft-cpu", "-j", "S1_opt", "--expect-state", "TODO", "start")
    assert r.returncode != 0 and "--expect-state" in (r.stdout + r.stderr)


# ---------------------------------------------------------------------------
# 单元
# ---------------------------------------------------------------------------
def _gm(name, pfrom, busy=0):
    return {"name": name, "tt": "k", "_seg": {"_from": pfrom, "max_jobs": 3},
            "steps": [{"kind": "R"}] * busy + [{"kind": "TODO"}]}


def test_gate_adds_busy_baseline():
    data = {"types": [{"key": "k", "materials": [_gm("a", "/p/tf_x.yaml")]}],
            "_busy_baseline": {"skill": {"k": 2}, "project": {"k\t/p/tf_x.yaml": 3}}}
    gate = W._SkillGate({"task_types": {"k": {"max_jobs": 10}}}, data)
    assert gate.busy("k") == 2
    assert gate.limit_hit("k", data["types"][0]["materials"][0])[0] == "project"


def test_busy_baseline_uses_progress_and_newer_ledger(tmp_path):
    cfg = {"_config_dir": str(tmp_path)}
    now = time.time()
    doc = {"coverage": {"k": now - 100}, "materials": [
        {"id": "p/a", "tt": "k", "project": "p", "busy": 1, "obs": now - 100},
        {"id": "p/b", "tt": "k", "project": "p", "busy": 0, "obs": now - 100},
        {"id": "p/c", "tt": "k", "project": "p", "busy": 5, "obs": now - 100}]}
    with open(P.ledger_path(cfg), "w") as f:
        f.write(json.dumps({"ts": now - 50, "tt": "k", "id": "p/b", "jobs": 1,
                            "project": "p", "from": "/f"}) + "\n")
        f.write(json.dumps({"ts": now - 500, "tt": "k", "id": "p/a", "jobs": 1}) + "\n")
    data = {"types": [{"key": "k", "materials": [{"name": "c", "qualified_name": "p/c"}]}]}
    base = P.busy_baseline(cfg, doc, data, {"k"}, {"p": "/f"})
    # a 的 1 + b 的账本 1（旧于观测的那条不算；c 在采集范围里不进基线）
    assert base["skill"] == {"k": 2}
    assert base["project"] == {"k\t/f": 2}


def test_progress_view_overlays_submissions(tmp_path):
    cfg = {"_config_dir": str(tmp_path)}
    data = {"types": [{"key": "k", "materials": [{
        "name": "a", "qualified_name": "p/a", "project": "p",
        "active": {"label": "S1", "kind": "TODO"},
        "steps": [{"label": "S1", "name": "s1", "kind": "TODO"}]}]}]}
    P.write_progress(cfg, data, observed=time.time() - 10)
    P.ledger_append(cfg, {"tt": "k", "qualified_name": "p/a"}, {"label": "S1"}, "42")
    raw = P.load_progress(cfg)["materials"][0]
    view = P.load_progress_view(cfg)["materials"][0]
    assert raw["steps"]["S1"] == "TODO"
    assert view["steps"]["S1"] == "PD" and view["job"]["id"] == "42"
    assert view["submitted_after_collect"] is True


def test_compact_from_progress_feeds_proposals():
    doc = {"queue": {}, "materials": [
        {"id": "p/a", "name": "a", "project": "p", "tt": "k", "active": "S1",
         "steps": {"S1": "FAIL"}, "fails": [{"step": "S1", "code": "x", "class": "C",
                                             "action": "human_review", "reason": "r"}]},
        {"id": "p/b", "name": "b", "project": "p", "tt": "k", "active": "S1",
         "steps": {"S1": "TODO"}}]}
    compact = S.compact_from_progress(doc, {"tt": "k"})
    props = {p["material"]: p for p in S.protocol.proposals(compact)}
    assert props["p/a"]["requires_human_review"] is True
    assert props["p/b"]["tool"] == "start_step"


def test_execute_grouped_batches_and_maps_outcomes():
    calls = []

    def run(argv):
        calls.append(argv)
        if argv[-1] != "start":
            return 0, "fetched", ""
        return 0, "\n".join([
            '[agent-event] {"event": "submitted", "id": "p/a", "name": "a"}',
            '[agent-event] {"event": "deferred", "id": "p/b", "name": "b"}',
            '[agent-event] {"event": "skipped", "material": "p/c", "id": "p/c"}']), ""
    acts = [{"action": "start_step", "tt": "k", "material": m, "step": "S1"}
            for m in ("p/a", "p/b", "p/c")] + [{"action": "sync_results", "tt": "k",
                                                   "material": "p/a"}]
    out = S.execute_grouped(acts, run)
    assert len(calls) == 2
    assert calls[0] == ["act", "-tt", "k", "-p", "p/a,p/b,p/c", "-j", "S1",
                        "--expect-state=TODO,PREP,SCANCEL", "start"]
    assert calls[1] == ["act", "-tt", "k", "-p", "p/a", "fetch"]
    assert [r["outcome"] for r in out["results"]] == [
        "submitted", "deferred", "skipped_stale", "ok"]
    assert out["outcomes"] == {"submitted": 1, "deferred": 1, "skipped_stale": 1, "ok": 1}
    assert "[agent-event]" not in out["groups"][0]["output"]


def test_gate_parses_project_value_flag():
    assert G.agent_command(["--project", "big", "-tt", "k", "advance"]) == "advance"
    assert G.agent_classify("advance", ["--project", "big", "advance"])[0] == "mutate"
    assert G.agent_classify("start", ["-p", "a", "--expect-state=TODO", "start"])[0] == "mutate"


def test_agent_setup_prints_absolute_mcp_entry(tmp_path):
    cfgp = tmp_path / "tf.yaml"
    cfgp.write_text("hang_check: false\n")
    r = subprocess.run([sys.executable, PROG, "-c", str(cfgp), "agent", "setup"],
                       capture_output=True, text=True, timeout=60)
    data = json.loads(r.stdout)["data"]
    srv = data["mcp_server"]["mcpServers"]["autozt"]
    assert os.path.isabs(srv["command"]) and srv["args"][-1] == "mcp"
    assert srv["env"]["AUTOZT_CONFIG"] == str(cfgp)
    assert "skipped_stale" in data["rules_md"]
    r = subprocess.run([sys.executable, PROG, "-c", str(cfgp), "agent", "setup", "--rules"],
                       capture_output=True, text=True, timeout=60)
    assert r.stdout.startswith("# AutoZT 调用规则")
