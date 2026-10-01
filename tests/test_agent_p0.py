# -*- coding: utf-8 -*-
"""P0 agent-usability fixes: local backend (fakeslurm), error mirroring, unknown -p,
register / conf_set, standalone gen helpers."""
import importlib.util
import json
import os
import subprocess
import sys
import time

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHIM = os.path.join(ROOT, "tools", "fakeslurm")
_AGENT_ENV = ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT", "GEMINI_CLI", "CODEX_SANDBOX",
              "AUTOZT_ACTOR", "AUTOZT_JSON_ERRORS")


def _env(**kw):
    env = {k: v for k, v in os.environ.items() if k not in _AGENT_ENV}
    env.update({k: str(v) for k, v in kw.items()})
    return env


def _shim(cmd, args, env, cwd=None):
    return subprocess.run([os.path.join(SHIM, cmd)] + list(args), capture_output=True,
                          text=True, env=env, cwd=cwd, timeout=60)


# ---------------------------------------------------------------- fakeslurm
def test_fakeslurm_lifecycle(tmp_path):
    env = _env(FAKESLURM_HOME=tmp_path / "fs", FAKESLURM_POLL="0.1", FAKESLURM_MAX_CPUS="1")
    wd = tmp_path / "w"
    wd.mkdir()
    (wd / "ok.sh").write_text("#!/bin/bash\n#SBATCH -J okjob\n#SBATCH --cpus-per-task=8\n"
                              "echo hi $SLURM_JOB_ID $SLURM_CPUS_PER_TASK\n")
    (wd / "bad.sh").write_text("#!/bin/bash\nexit 3\n")
    (wd / "long.sh").write_text("#!/bin/bash\nsleep 60\n")
    r = _shim("sbatch", ["ok.sh"], env, cwd=wd)
    assert r.returncode == 0 and r.stdout.startswith("Submitted batch job ")
    j1 = r.stdout.split()[-1]
    j2 = _shim("sbatch", ["bad.sh"], env, cwd=wd).stdout.split()[-1]
    deadline = time.time() + 20
    while time.time() < deadline:
        out = _shim("sacct", ["-j", "%s,%s" % (j1, j2), "-n", "-P", "-X", "-o",
                              "JobID,State"], env).stdout
        if "COMPLETED" in out and "FAILED" in out:
            break
        time.sleep(0.2)
    assert "%s|COMPLETED" % j1 in out and "%s|FAILED" % j2 in out
    # cpus clipped to the budget (8 -> 1) so the job could start at all
    assert (wd / ("slurm-%s.out" % j1)).read_text().split() == ["hi", j1, "1"]
    j3 = _shim("sbatch", ["long.sh"], env, cwd=wd).stdout.split()[-1]
    time.sleep(0.5)
    q = _shim("squeue", ["-u", os.environ.get("USER") or __import__("getpass").getuser(),
                         "-h", "-o", "%i|%Z|%T"], env).stdout.strip().splitlines()
    assert any(ln.startswith(j3 + "|" + str(wd) + "|") for ln in q)
    assert _shim("scancel", [j3], env).returncode == 0
    out = _shim("sacct", ["-j", j3, "-n", "-P", "-X", "-o", "JobID,State"], env).stdout
    assert out.strip() == "%s|CANCELLED" % j3


def test_local_path_prefix_toggle(monkeypatch):
    import autozt  # noqa: F401  (包命名空间里 collect 是同名函数，取模块本身)
    collect = sys.modules["autozt.collect"]
    monkeypatch.setenv("AUTOZT_FAKESLURM", "1")
    assert collect._effective_remote_path_prefix({}, None).endswith("tools/fakeslurm")
    assert collect._effective_remote_path_prefix(
        {"remote_path_prefix": "/opt/x/bin"}, "").startswith("/opt/x/bin:")
    monkeypatch.setenv("AUTOZT_FAKESLURM", "0")
    assert collect._effective_remote_path_prefix({}, None) == ""
    rc, out = collect.run_remote({"remote_path_prefix": "/opt/x/bin"}, "echo $PATH", host="")
    assert rc == 0 and out.startswith("/opt/x/bin:")


# ---------------------------------------------------------------- CLI: register + errors
@pytest.fixture()
def sandbox(tmp_path):
    proj = tmp_path / "proj"
    proj.mkdir()
    cfg = tmp_path / "tf.yaml"
    cfg.write_text('host: ""\nproject_roots:\n  - %s\ntask_types:\n'
                   '  fit-fc-thermal: {hpc: {name: local, ssh_host: ""}}\n' % proj)
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "POSCAR").write_text("Si\n1.0\n0 2.7 2.7\n2.7 0 2.7\n2.7 2.7 0\nSi\n2\n"
                               "Direct\n0 0 0\n0.25 0.25 0.25\n")
    return {"cfg": str(cfg), "proj": proj, "ds": ds, "tmp": tmp_path}


def _cli(sb, *args, **env):
    return subprocess.run([sys.executable, "-m", "autozt", "-c", sb["cfg"]] + list(args),
                          capture_output=True, text=True, cwd=ROOT,
                          env=_env(HOME=sb["tmp"], **env), timeout=300)


def test_register_and_unknown_material(sandbox):
    r = _cli(sandbox, "-tt", "fit-fc-thermal", "-p", "Si1", "register",
             "--dataset", str(sandbox["ds"]), "--cluster", "local", "--json")
    assert r.returncode == 0, r.stdout + r.stderr
    rep = json.loads(r.stdout.strip().splitlines()[-1])
    assert rep["ok"] and rep["hpc"] == "local"
    ps = sandbox["proj"] / "Si1" / "fit-fc-thermal" / "project_setting"
    conf = (ps / "templates" / "step1_fit" / "step.conf").read_text()
    assert "FIT_INPUT_DIR = %s" % sandbox["ds"] in conf
    assert (sandbox["proj"] / "Si1" / "POSCAR").is_file()
    # idempotent
    assert _cli(sandbox, "-tt", "fit-fc-thermal", "-p", "Si1", "register").returncode == 0
    # the new material is visible; a wrong -p is an error (was: silent empty table, rc=0)
    ok = _cli(sandbox, "-tt", "fit-fc-thermal", "-p", "Si1", "summary", "--json")
    assert ok.returncode == 0, ok.stdout + ok.stderr
    bad = _cli(sandbox, "-p", "Si2", "summary", "--json")
    assert bad.returncode == 1
    err = json.loads(bad.stdout.strip().splitlines()[-1])
    assert err["ok"] is False and "Si2" in err["error"] and "Si1" in err["error"]
    # agent env: mirrored to stdout as a plain line; human terminal: stderr only
    a = _cli(sandbox, "-p", "Si2", "summary", CLAUDECODE="1")
    assert a.returncode == 1 and a.stdout.startswith("[autozt error]")
    h = _cli(sandbox, "-p", "Si2", "summary")
    assert h.returncode == 1 and h.stdout == "" and "找不到材料" in h.stderr


def test_register_rejects_bad_input(sandbox):
    r = _cli(sandbox, "-tt", "band-dft-cpu", "-p", "X", "register",
             "--dataset", str(sandbox["ds"]), "--json")
    assert r.returncode == 1 and not (sandbox["proj"] / "X").exists()
    r = _cli(sandbox, "-tt", "fit-fc-thermal", "-p", "../evil", "register", "--json")
    assert r.returncode == 1


# ---------------------------------------------------------------- MCP surface
def test_mcp_workflow_has_register_and_conf_set(monkeypatch):
    monkeypatch.delenv("AUTOZT_MCP_PROFILE", raising=False)
    monkeypatch.delenv("AUTOZT_MCP_READONLY", raising=False)
    from autozt import mcp
    names = {t["name"] for t in mcp._tools_list()}
    assert {"register_material", "conf_set", "conf_get"} <= names
    assert mcp.BY_NAME["conf_set"][3] == "mutate"
    assert mcp.BY_NAME["register_material"][3] == "mutate"
    r = mcp.call_tool("conf_set", {"material": "m", "step": "s", "key": "A;rm -rf",
                                   "value": "1"})["structuredContent"]
    assert r["ok"] is False and "key" in r["error"]


# ---------------------------------------------------------------- standalone gen helpers
def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_stepconf_warn_mode(tmp_path, capsys):
    sc = _load(os.path.join(ROOT, "skill", "_common", "opt", "stepconf.py"), "sc_p0")
    (tmp_path / "step.conf").write_text("[params]\nFUNC = auto\nFIT_METHODS = x\n")
    conf = sc.load({"FIT_METHODS": ("auto", "str")}, "s1", str(tmp_path), strict="warn")
    assert conf["FIT_METHODS"] == "x" and "FUNC" in capsys.readouterr().err
    (tmp_path / "step.conf").write_text("[params]\nFIT_METHOD = x\n")
    with pytest.raises(SystemExit, match="FIT_METHODS"):
        sc.load({"FIT_METHODS": ("auto", "str")}, "s1", str(tmp_path), strict="warn")
    with pytest.raises(SystemExit):     # strict=True keeps the old behaviour
        sc.load({"OTHER": ("auto", "str")}, "s1", str(tmp_path))


def test_resolve_submit_step_subdir():
    sys.path.insert(0, os.path.join(ROOT, "skill", "fit-fc-thermal"))
    try:
        fc = _load(os.path.join(ROOT, "skill", "fit-fc-thermal", "fc_common.py"), "fc_p0")
        p = fc.resolve_submit(os.path.join(ROOT, "skill", "fit-fc-thermal"), "submit_kappa")
        assert str(p).endswith(os.path.join("templates", "step2_kappa", "submit_kappa.tpl"))
    finally:
        sys.path.pop(0)
