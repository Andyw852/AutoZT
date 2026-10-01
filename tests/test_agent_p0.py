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


# ---------------------------------------------------------------- local backend vs global host
@pytest.fixture()
def remote_default(tmp_path):
    """用户真实布局：全局 host 指向集群、技能类型带集群 work_dir；本机没有 ssh 可用。"""
    proj = tmp_path / "proj"
    proj.mkdir()
    cfg = tmp_path / "tf.yaml"
    cfg.write_text("host: no-such-cluster\nproject_roots:\n  - %s\ntask_types:\n"
                   "  fit-fc-thermal:\n    work_dir: /public/home/nobody/work\n" % proj)
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "POSCAR").write_text("Si\n1.0\n0 2.7 2.7\n2.7 0 2.7\n2.7 2.7 0\nSi\n2\n"
                               "Direct\n0 0 0\n0.25 0.25 0.25\n")
    return {"cfg": str(cfg), "proj": proj, "ds": ds, "tmp": tmp_path}


def test_register_local_ignores_global_host(remote_default):
    sb = remote_default
    r = _cli(sb, "-tt", "fit-fc-thermal", "-p", "Si9", "register",
             "--dataset", str(sb["ds"]), "--cluster", "local", "--json")
    assert r.returncode == 0, r.stdout + r.stderr
    # conf / dir only read configuration: no ssh (the host does not exist), fast
    t0 = time.time()
    c = _cli(sb, "-tt", "fit-fc-thermal", "-p", "Si9", "-j", "step1_fit", "conf",
             "--set", "params.FIT_RMSE_FRAMES=5")
    assert c.returncode == 0, c.stdout + c.stderr
    assert time.time() - t0 < 60
    assert "no-such-cluster" not in c.stdout + c.stderr
    # cluster layer (setting/local.yaml) beats the skill's hard-coded conda path
    for ln in c.stdout.splitlines():
        if ln.startswith(("CONDA_SH", "CONDA_ENV")):
            assert "[cluster:local]" in ln, ln
    d = _cli(sb, "-tt", "fit-fc-thermal", "-p", "Si9", "dir")
    assert d.returncode == 0, d.stdout + d.stderr
    wd = d.stdout.strip().splitlines()[-1]
    assert not wd.startswith("/public/") or os.access("/", os.W_OK)
    assert wd.endswith(os.path.join("Si9", "fit-fc-thermal"))


def test_local_host_marker_runs_locally():
    import autozt  # noqa: F401
    collect = sys.modules["autozt.collect"]
    cfg = {"host": "no-such-cluster"}
    assert collect._norm_host(cfg, collect.LOCAL_HOST) == ""
    assert collect._norm_host(cfg, "__default__") == "no-such-cluster"
    assert collect._ssh_cmd(cfg, collect.LOCAL_HOST, ["echo", "hi"])[:2] == ["bash", "-c"]
    rc, out = collect.run_remote(cfg, "echo local-ok", host=collect.LOCAL_HOST)
    assert rc == 0 and out.endswith("local-ok")


def test_local_work_dir_skips_unwritable(tmp_path, monkeypatch, capsys):
    import autozt  # noqa: F401
    bs = sys.modules["autozt.bootstrap"]
    ok_dir = tmp_path / "w"
    monkeypatch.setattr(bs, "_writable_here",
                        lambda p: str(p).startswith(str(tmp_path)))
    m = {"name": "X"}
    bs._local_work_dir(m, [("project_setting/setting.yaml", "/public/a"),
                           ("hpc.yaml", None), ("task_types.k.work_dir", str(ok_dir))])
    assert m["work_dir_eff"] == str(ok_dir)
    assert m["work_dir_src"] == "task_types.k.work_dir"
    assert "/public/a" in capsys.readouterr().err


def test_preflight_conda_regex_does_not_cross_lines():
    from autozt import preflight
    tpl = open(os.path.join(ROOT, "skill", "fit-fc-thermal", "templates", "step1_fit",
                            "submit_fcfit_pheasy.tpl")).read()
    empty = tpl.replace("{{CONDA_SH}}", "").replace("{{CONDA_ENV}}", "")
    assert preflight.parse_conda_activations(empty) == []
    full = tpl.replace("{{CONDA_SH}}", "/x/etc/profile.d/conda.sh").replace(
        "{{CONDA_ENV}}", "env1")
    assert preflight.parse_conda_activations(full) == [
        ("sh", "/x/etc/profile.d/conda.sh"), ("env", "env1")]


# ---------------------------------------------------------------- envcheck / local templates
def _ec():
    import autozt.envcheck  # noqa: F401  (不在包命名空间注入表里，按子模块导入)
    return sys.modules["autozt.envcheck"]


def test_envcheck_reports_missing_and_hints():
    ec = _ec()
    cfg = {"task_types": {
        "fake-py": {"_skill_requires": {"python": ["numpy", "no-such-pkg-zz"],
                                        "optional_python": ["also-missing-zz"],
                                        "exe": ["sh", "no_such_exe_zz"]}},
        "fake-vasp": {"_skill_requires": {"exe": ["vasp_std"], "potcar": "potcar_dir"}},
        "fake-gpu": {"_skill_requires": {"python": ["numpy"], "gpu": True}},
        "fake-none": {}}}
    rep = ec.check(cfg, ["fake-py", "fake-vasp", "fake-gpu", "fake-none"], "local")
    assert rep["ok"], rep
    by = {s["skill"]: s for s in rep["skills"]}
    assert by["fake-py"]["status"] == "missing"
    assert set(by["fake-py"]["missing"]) == {"py:no-such-pkg-zz", "exe:no_such_exe_zz"}
    assert by["fake-py"]["optional_missing"] == ["py:also-missing-zz"]
    assert any("pip install no-such-pkg-zz" in h for h in by["fake-py"]["hints"])
    assert "potcar_dir" in by["fake-vasp"]["missing"]          # local.yaml potcar_dir is ""
    import shutil as _sh
    assert ("gpu" in by["fake-gpu"]["missing"]) == (not _sh.which("nvidia-smi"))
    assert by["fake-none"]["status"] == "unknown"
    assert "合计" in ec.render(rep)


def test_check_cli_json_and_strict(remote_default):
    r = _cli(remote_default, "-tt", "te-screen", "check", "--json", "--strict")
    doc = json.loads(r.stdout)
    assert doc["cluster"] == "local" and doc["skills"][0]["skill"] == "te-screen"
    assert r.returncode == (0 if doc["skills"][0]["status"] == "ready" else 1)
    bad = _cli(remote_default, "check", "--cluster", "no-such-cluster")
    assert bad.returncode == 1


def test_local_start_gate(monkeypatch, capsys):
    import autozt  # noqa: F401
    wf = sys.modules["autozt.workflow"]
    ec = _ec()
    monkeypatch.delenv("AUTOZT_SKIP_ENV_CHECK", raising=False)
    monkeypatch.setattr(ec, "missing_for", lambda cfg, k, c: {
        "status": "missing", "missing": ["exe:vasp_std"], "hints": ["装 VASP"]})
    m = {"host_eff": "@local", "hpc_name": "local", "tt": "band-dft-cpu"}
    assert wf._local_env_ready({}, m, {}, "t") is False
    assert "exe:vasp_std" in capsys.readouterr().out
    monkeypatch.setenv("AUTOZT_SKIP_ENV_CHECK", "1")
    assert wf._local_env_ready({}, m, {}, "t") is True
    monkeypatch.delenv("AUTOZT_SKIP_ENV_CHECK")
    assert wf._local_env_ready({}, dict(m, host_eff="jzzn"), {}, "t") is True  # clusters: no gate


def test_hpc_yaml_writer_roundtrip(tmp_path):
    import autozt  # noqa: F401
    ops = sys.modules["autozt.ops"]
    from autozt import _load_yaml_file
    d = _load_yaml_file(os.path.join(ROOT, "setting", "local.yaml"))
    d["x"] = "a: b"
    p = str(tmp_path / "hpc.yaml")
    ops._write_hpc_yaml(p, d, "t")
    assert _load_yaml_file(p) == d
    assert d["template_map"] == {} and d["ssh_host"] == ""


def test_register_local_copies_no_cluster_templates(remote_default):
    sb = remote_default
    poscar = sb["ds"] / "POSCAR"
    r = _cli(sb, "-tt", "band-dft-cpu", "-p", "SiV", "register", "--poscar", str(poscar),
             "--cluster", "local", "--json")
    assert r.returncode == 0, r.stdout + r.stderr
    ps = sb["proj"] / "SiV" / "band-dft-cpu" / "project_setting"
    assert not list(ps.rglob("submit_*.tpl"))
    from autozt import _load_yaml_file
    assert _load_yaml_file(str(ps / "hpc.yaml")).get("template_map") == {}


def test_local_templates_render():
    import re as _re
    import autozt  # noqa: F401
    wf = sys.modules["autozt.workflow"]
    from autozt import _load_yaml_file
    tdir = os.path.join(ROOT, "setting", "local", "templates")
    keys = {   # placeholders each gen fills (see the gen scripts' write_submit/render calls)
        "submit_mlff_relax.tpl": ["JOBNAME", "CONDA_SH", "CONDA_ENV", "MACE_CMD", "LOG"],
        "submit_mlff_opt.tpl": ["JOBNAME", "CONDA_SH", "CONDA_ENV", "MACE_CMD", "LOG"],
        "submit_mlff.tpl": ["JOBNAME", "CONDA_SH", "CONDA_ENV", "MACE_CMD"],
        "submit_fc.tpl": ["JOBNAME", "CONDA_SH", "CONDA_ENV"],
        "submit_shengbte_klm.tpl": ["JOBNAME", "CONDA_SH", "CONDA_ENV", "SHENGBTE_EXE",
                                    "NTASKS", "CPUS_PER_TASK", "SB_EXTRACT"],
        "submit_hamgnn.tpl": ["JOBNAME", "CONDA_SH", "CONDA_ENV", "NTHREADS", "CMD"],
        "submit_hamgnn_mpi.tpl": ["JOBNAME", "CONDA_SH", "CONDA_ENV", "NPROC", "CMD"],
        "submit_amset.tpl": ["JOBNAME", "AMSET_CMD", "AMSET_ENV"],
    }
    profiles = _load_yaml_file(os.path.join(ROOT, "setting", "local.yaml"))["vasp"]
    names = sorted(os.listdir(tdir))
    assert len(names) == 17
    for n in names:
        text = open(os.path.join(tdir, n), encoding="utf-8").read()
        if _re.fullmatch(r"submit_(std|gam|ncl)_(0d|2d|3d)\.tpl", n):
            text = wf.render_vasp_template(text, n, "step2_static", profiles)
            text = text.replace("{{JOBNAME}}", "j")
        else:
            for k in keys[n]:
                text = text.replace("{{%s}}" % k, "2")
        assert not _re.findall(r"\{\{[A-Z_0-9]+\}\}", text), n
        r = subprocess.run(["bash", "-n"], input=text, capture_output=True, text=True)
        assert r.returncode == 0, (n, r.stderr)


def test_agent_setup_write_merges(tmp_path, monkeypatch):
    import argparse
    from autozt import agent_cli as A
    monkeypatch.setattr(A, "ROOT", str(tmp_path))
    (tmp_path / ".mcp.json").write_text(json.dumps({"mcpServers": {"other": {"command": "x"}}}))
    res, rc = A._cmd_setup(argparse.Namespace(config=None, write="all"))
    assert rc == 0
    doc = json.loads((tmp_path / ".mcp.json").read_text())
    assert set(doc["mcpServers"]) == {"other", "autozt"}
    assert doc["mcpServers"]["autozt"]["args"][-1] == "mcp"
    assert json.loads((tmp_path / ".cursor" / "mcp.json").read_text())["mcpServers"]["autozt"]
    (tmp_path / ".mcp.json").write_text("{broken")
    _res, rc = A._cmd_setup(argparse.Namespace(config=None, write="claude"))
    assert rc == 1 and (tmp_path / ".mcp.json").read_text() == "{broken"


def test_hang_helpers_use_job_host(monkeypatch):
    import autozt  # noqa: F401
    ops = sys.modules["autozt.ops"]
    seen = []
    monkeypatch.setattr(autozt, "run_remote", lambda cfg, line, host="__default__", **k:
                        (seen.append(host) or (0, "")))
    ops._hung_scancel_wait({}, "1", timeout=1, host="a800")
    ops._hung_incar_fix({}, "/w", 1, host="a800")
    assert seen and all(h == "a800" for h in seen)
