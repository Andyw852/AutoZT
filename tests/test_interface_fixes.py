# -*- coding: utf-8 -*-
"""Regression tests for the 2026-10 agent-interface fixes (ported to local master).

Covers the five still-present bugs:
2. conf_get schema declares `step` required
3. doctor judges clusters by what materials actually use, not by unused skill defaults
4. preflight: input check before execution, result validation after recovery
5. probe reads local-cluster materials locally (never ssh); missing ssh is a clean error
6. hpc / adopt confirmations give a -y hint instead of an EOFError traceback without a TTY

All sandboxed: local cluster, temporary config, no ssh, no job submission.
"""
import json
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
_AGENT_ENV = ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT", "GEMINI_CLI", "CODEX_SANDBOX",
              "CODEX_SANDBOX_NETWORK_DISABLED", "AUTOZT_ACTOR", "AUTOZT_JSON_ERRORS",
              "AUTOZT_CONFIG", "AUTOZT_CFG", "AUTOZT_MCP_PROFILE", "AUTOZT_MCP_READONLY")

POSCAR_2D_C = """MoS2 monolayer, vacuum along c
1.0
 3.19 0.0 0.0
-1.595 2.7626 0.0
 0.0 0.0 20.0
Mo S
1 2
Direct
0.0 0.0 0.5
0.33333 0.66667 0.578
0.33333 0.66667 0.422
"""
POSCAR_2D_A = """MoS2 monolayer, vacuum along a
1.0
 20.0 0.0 0.0
 0.0 3.19 0.0
 0.0 -1.595 2.7626
Mo S
1 2
Direct
0.5 0.0 0.0
0.578 0.33333 0.66667
0.422 0.33333 0.66667
"""
POSCAR_3D = """Si diamond
1.0
0 2.715 2.715
2.715 0 2.715
2.715 2.715 0
Si
2
Direct
0 0 0
0.25 0.25 0.25
"""


def _env(home, **kw):
    env = {k: v for k, v in os.environ.items() if k not in _AGENT_ENV}
    env.update({"HOME": str(home), "AUTOZT_APPROVAL_GRAPH": "0", "AUTOZT_MCP_ELICITATION": "0"})
    env.update({k: str(v) for k, v in kw.items()})
    return env


def _cli(sb, *args, stdin=subprocess.DEVNULL, **env):
    return subprocess.run([sys.executable, "-m", "autozt", "-c", sb["cfg"]] + list(args),
                          capture_output=True, text=True, cwd=ROOT, stdin=stdin,
                          env=_env(sb["home"], **env), timeout=300)


def _mcp(sb, name, args, **env):
    p = subprocess.run([sys.executable, "-m", "autozt", "mcp", "--call", name,
                        json.dumps(args, ensure_ascii=False)],
                       capture_output=True, text=True, cwd=ROOT, timeout=300,
                       env=_env(sb["home"], AUTOZT_CONFIG=sb["cfg"], AUTOZT_MCP_PROFILE="full",
                                **env))
    got = json.loads(p.stdout)
    return got.get("structuredContent", got)


@pytest.fixture(scope="module")
def sb(tmp_path_factory):
    """One local-cluster opt-dft-cpu material (2D MoS2) in a throw-away config."""
    tmp = tmp_path_factory.mktemp("iface")
    proj, home = tmp / "proj", tmp / "home"
    proj.mkdir()
    home.mkdir()
    cfg = tmp / "tf.yaml"
    cfg.write_text('''host: ""
project_roots:
  - %s
auto_advance: false
auto_watch: false
''' % proj)
    pos = tmp / "POSCAR_2D"
    pos.write_text(POSCAR_2D_C)
    (tmp / "empty_results").mkdir()
    box = {"cfg": str(cfg), "proj": proj, "home": home, "tmp": tmp, "poscar": str(pos),
           "results": str(tmp / "empty_results")}
    r = _cli(box, "-tt", "opt-dft-cpu", "-p", "MoS2", "register", "--poscar", str(pos),
             "--cluster", "local", "--json")
    assert r.returncode == 0, r.stdout + r.stderr
    return box


# ---------------------------------------------------------------- 2. schema conformance
def test_conf_get_declares_step():
    from autozt import mcp
    assert mcp._validate_args("conf_get", {"material": "X"}) == "missing required argument: step"
    assert mcp._validate_args("conf_get", {"material": "X", "step": "S1_opt"}) is None


# ---------------------------------------------------------------- 3. doctor
def test_doctor_checks_clusters_materials_actually_use(tmp_path):
    proj, home = tmp_path / "proj", tmp_path / "home"
    proj.mkdir()
    home.mkdir()
    cfg = tmp_path / "tf.yaml"
    cfg.write_text('''host: ""
project_roots:
  - %s
auto_advance: false
task_types:
  opt-dft-cpu: {hpc: no-such-cluster-xyz}
''' % proj)
    pos = tmp_path / "POSCAR"
    pos.write_text(POSCAR_3D)
    box = {"cfg": str(cfg), "home": home}
    assert _cli(box, "-tt", "opt-dft-cpu", "-p", "SiLocal", "register", "--poscar", str(pos),
                "--cluster", "local", "--json").returncode == 0
    r = _cli(box, "doctor", "--json")
    rep = json.loads(r.stdout)
    codes = {(f["level"], f["code"]) for f in rep["findings"]}
    # the only material runs locally: the unused skill default is a warning, not an error
    assert ("warn", "cluster_default_unconfigured") in codes
    assert not any(c == "cluster_config_missing" for _lv, c in codes)
    assert r.returncode == 0, r.stdout
    # a material that really uses the unconfigured cluster is an error that names it
    assert _cli(box, "-tt", "opt-dft-cpu", "-p", "SiRemote", "register", "--poscar", str(pos),
                "--json").returncode == 0
    r = _cli(box, "doctor", "--json")
    rep = json.loads(r.stdout)
    miss = [f for f in rep["findings"] if f["code"] == "cluster_config_missing"]
    assert r.returncode == 1 and miss and miss[0]["level"] == "error"
    assert any("SiRemote" in x for x in miss[0]["examples"])
    assert not any("SiLocal" in x for x in miss[0]["examples"])


# ---------------------------------------------------------------- 4. preflight
def _poscar(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_preflight_inputs_before_execution(tmp_path):
    from autozt.science import preflight_inputs
    ok = preflight_inputs(_poscar(tmp_path, "c", POSCAR_2D_C), dimension="2D",
                          temperature=[300, 600], carrier=[1e19, 1e20])
    assert ok["stage"] == "inputs" and ok["status"] == "pass"
    assert ok["conversation_state"] == "ready_to_execute"
    assert ok["structure"]["dimension_detected"] == "2D"
    assert ok["structure"]["vacuum_axis"] == "c"

    def failed(res):
        return {c["code"] for c in res["checks"] if not c["ok"] and c["severity"] == "error"}
    si = _poscar(tmp_path, "si", POSCAR_3D)
    assert failed(preflight_inputs(si, dimension="2D")) == {"dimension"}
    along_a = _poscar(tmp_path, "a", POSCAR_2D_A)
    assert failed(preflight_inputs(along_a, dimension="2D")) == {"vacuum_axis"}
    bad = preflight_inputs(_poscar(tmp_path, "c2", POSCAR_2D_C), temperature=[-5], carrier=[0])
    assert failed(bad) == {"temperature_grid", "carrier_grid"}
    assert bad["conversation_state"] == "preflight_review"
    assert failed(preflight_inputs(str(tmp_path / "missing"))) == {"poscar"}
    # an unusual carrier range is a warning, it does not block
    warn = preflight_inputs(_poscar(tmp_path, "c3", POSCAR_2D_C), carrier=[1e30])
    assert warn["status"] == "pass"
    assert any(c["code"] == "carrier_range" and not c["ok"] for c in warn["checks"])


def test_preflight_result_stage_means_results_are_usable(tmp_path):
    from autozt.science import preflight
    units = {"S": "uV/K", "sigma": "S/m", "kappa_e": "W/m/K", "kappa_L": "W/m/K", "T": "K",
             "PF": "W/m/K^2"}
    (tmp_path / "zt_summary.json").write_text(json.dumps({
        "ZT_DONE": True, "units": units, "dim": "3D", "temperatures": [300],
        "doping_cm-3": [1e18], "sources": {"transport": "transport.json"},
        "formula": "zT = S^2 sigma T / (kappa_e + kappa_L)"}))
    (tmp_path / "transport.json").write_text(json.dumps({"temperatures": [300]}))
    (tmp_path / "kappa_summary.json").write_text(json.dumps({"KAPPA_DONE": True,
                                                              "temperatures": [300]}))
    got = preflight(str(tmp_path), dimension="3D", temperature=[300], carrier=[1e18])
    assert got["stage"] == "results" and got["status"] == "pass", got["checks"]
    # validated results are ready to report, not "ready to execute"
    assert got["conversation_state"] == "completed"
    assert got["next_action"] == "review_result"


def test_preflight_mcp_modes_and_workflow_order(sb):
    got = _mcp(sb, "preflight", {})
    assert not got["ok"] and "result_dir" in got["error"] and "poscar" in got["error"]
    got = _mcp(sb, "preflight", {"material": "MoS2", "tt": "opt-dft-cpu", "dimension": "2D",
                                 "temperature": [300], "carrier": [1e19]})
    assert got["ok"] and got["data"]["stage"] == "inputs"
    assert got["data"]["conversation_state"] == "ready_to_execute"
    got = _mcp(sb, "preflight", {"result_dir": sb["results"], "dimension": "2D"})
    assert got["ok"] and got["data"]["stage"] == "results"
    from autozt import mcp, science
    text = json.dumps(mcp._prompt_get("compute-zt"), ensure_ascii=False)
    assert text.index("preflight") < text.index("cycle(execute=false)") < text.index(
        "preflight(result_dir)") < text.rindex("results")
    plan = science.research_plan("compute zT", [{"name": "zt-dft-cpu"}], dimension="2D",
                                 material="MoS2", temperature=[300], carrier=[1e19])
    modes = [(c["tool"], c.get("mode")) for c in plan["review_card"]["planned_calls"]]
    assert modes.index(("preflight", "inputs")) < modes.index(("cycle", "execute")) < modes.index(
        ("preflight", "results")) < modes.index(("results", None))


def test_agent_cli_preflight_input_mode(sb):
    p = subprocess.run([sys.executable, "-m", "autozt", "agent", "preflight", "-c", sb["cfg"],
                        "-tt", "opt-dft-cpu", "-p", "MoS2", "--dimension", "2D",
                        "--temperature", "300", "--carrier", "1e19"],
                       capture_output=True, text=True, cwd=ROOT, env=_env(sb["home"]), timeout=300)
    got = json.loads(p.stdout)
    assert p.returncode == 0 and got["ok"], p.stdout + p.stderr
    assert got["data"]["stage"] == "inputs" and got["data"]["status"] == "pass"
    p = subprocess.run([sys.executable, "-m", "autozt", "agent", "preflight", "-c", sb["cfg"]],
                       capture_output=True, text=True, cwd=ROOT, env=_env(sb["home"]), timeout=300)
    assert p.returncode == 2 and "--result-dir" in p.stdout


# ---------------------------------------------------------------- 5. probe
def _bin_with(tmp_path, tools, fake_ssh=None):
    """A PATH holding only the given tools (symlinked), plus an optional fake ssh."""
    b = tmp_path / "bin"
    b.mkdir(exist_ok=True)
    for tool in tools:
        src = shutil.which(tool)
        if src and not (b / tool).exists():
            os.symlink(src, b / tool)
    if fake_ssh:
        (b / "ssh").write_text('''#!/bin/sh
touch %s
exit 255
''' % fake_ssh)
        (b / "ssh").chmod(0o755)
    return str(b)


def test_probe_local_material_never_uses_ssh(sb, tmp_path):
    marker = tmp_path / "ssh_was_called"
    path = _bin_with(tmp_path, ("bash", "sh", "env", "python3"), fake_ssh=marker)
    r = _cli(sb, "-tt", "opt-dft-cpu", "-p", "MoS2", "probe", PATH=path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert not marker.exists()
    assert json.loads(r.stdout[r.stdout.index("{"):r.stdout.rindex("}") + 1])["jobs"] == []


def test_probe_missing_ssh_is_a_clean_error(sb, tmp_path):
    path = _bin_with(tmp_path, ("bash", "sh", "env", "python3"))
    r = _cli(sb, "-tt", "opt-dft-cpu", "-p", "MoS2", "--host", "some-cluster", "probe", PATH=path)
    assert r.returncode == 1
    assert "Traceback" not in r.stdout + r.stderr
    assert "ssh" in r.stderr


# ---------------------------------------------------------------- 6. confirmations
def test_confirmation_without_tty_gives_hint_not_traceback(sb):
    r = _cli(sb, "-tt", "opt-dft-cpu", "-p", "MoS2", "hpc", "local")
    out = r.stdout + r.stderr
    assert r.returncode != 0
    assert "Traceback" not in out and "EOFError" not in out
    assert "-y" in out


def test_ask_confirm_names_the_command(monkeypatch):
    import io
    from autozt.workflow import _ask_confirm
    monkeypatch.setattr(sys, "stdin", io.StringIO(""))
    with pytest.raises(SystemExit) as exc:
        _ask_confirm("clean? [y/N] ", example="tf -p M clean -y")
    assert "tf -p M clean -y" in str(exc.value)
