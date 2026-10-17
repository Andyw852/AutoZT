# -*- coding: utf-8 -*-
"""[amset_env_python] S7 / S7.1 在 conda activate 之后用 $CONDA_PREFIX/bin/python（amset），不靠 PATH 里的 python3。

2026-10-17 Mg2C60/Mg4C60 S7.1：deform_read.log 只有 "ModuleNotFoundError: No module named 'numpy'"，
而手动 conda activate amset 后同一条命令正常。PATH 最前面有别的 python3（remote_path_prefix，或登录 shell
已激活同一环境、再 activate 不重排 PATH）时就会这样。这里用假的 conda.sh 复现：activate 只设 CONDA_PREFIX、
不动 PATH，PATH 最前面放一个"没有 numpy"的 python3。
"""
import importlib.util
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
S7 = REPO / "skill" / "ke-dft-cpu" / "step7_deform"
S71 = S7 / "step7b_read" / "gen_step9b_deform_read.py"


def _load(path, name):
    for p in (str(REPO), str(REPO / "skill" / "ke-dft-cpu"), str(S7)):
        if p not in sys.path:
            sys.path.insert(0, p)
    spec = importlib.util.spec_from_file_location(name, str(path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture()
def fake_env(monkeypatch):
    d = Path(tempfile.mkdtemp())
    env = d / "envs" / "amset" / "bin"
    env.mkdir(parents=True)
    os.symlink(sys.executable, env / "python")
    (env / "amset").write_text("#!/bin/sh\necho AMSET_OK \"$@\"\n")
    (env / "amset").chmod(0o755)
    bad = d / "pybin"
    bad.mkdir()
    (bad / "python3").write_text("#!/bin/sh\necho \"ModuleNotFoundError: No module named 'numpy'\" >&2\nexit 1\n")
    (bad / "python3").chmod(0o755)
    (bad / "amset").write_text("#!/bin/sh\necho WRONG_AMSET >&2\nexit 1\n")
    (bad / "amset").chmod(0o755)
    (d / "conda.sh").write_text('conda () { [ "$1" = activate ] && export CONDA_PREFIX="%s"; }\n'
                                % (d / "envs" / "amset"))
    monkeypatch.setenv("AUTOZT_CONDA_SH", str(d / "conda.sh"))
    monkeypatch.setenv("AUTOZT_AMSET_ENV", "amset")
    monkeypatch.chdir(d)                                   # 没有 step.conf -> 走环境变量
    yield d, bad
    shutil.rmtree(d, ignore_errors=True)


@pytest.mark.skipif(not shutil.which("bash"), reason="需要 bash")
def test_s71_uses_conda_prefix_python(fake_env):
    d, bad = fake_env
    m = _load(S71, "gen9b_envtest")
    cmd = "%s && cd %s && { %s; %s -c %s ; }" % (
        m.AMSET_ENV_SRC, shlex.quote(str(d)), m._ENV_ECHO, m.AMSET_PY,
        shlex.quote("import sys; print('PY', sys.executable)"))
    env = dict(os.environ, PATH="%s:%s" % (bad, os.environ["PATH"]))
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    assert "PY " in r.stdout and "[env] CONDA_PREFIX=" in r.stdout
    assert "(PATH 里的 python3=%s)" % (bad / "python3") in r.stdout      # 日志里看得出 PATH 被谁占了


def test_s71_no_bare_python3_after_activate():
    src = S71.read_text(encoding="utf-8")
    assert "&& python3 -c" not in src and "cd %s && python3" not in src
    assert src.count("AMSET_PY, ") + src.count("AMSET_PY,\n") >= 3      # NBANDS 包装 / band_edges / 真空 h5


@pytest.mark.skipif(not shutil.which("bash"), reason="需要 bash")
def test_s7_deform_create_uses_conda_prefix_amset(fake_env):
    d, bad = fake_env
    s = (S7 / "gen_step9_deform.py").read_text(encoding="utf-8")
    assert '"${CONDA_PREFIX:?conda activate 没生效}/bin/amset" deform create' in s
    m = _load(S7 / "gen_step9_deform.py", "gen9_envtest")
    cmd = ('%s && cd %s && { echo "[env] CONDA_PREFIX=$CONDA_PREFIX"; '
           '"${CONDA_PREFIX:?conda activate 没生效}/bin/amset" deform create -f POSCAR ; }'
           % (m.AMSET_ENV_SRC, shlex.quote(str(d))))
    env = dict(os.environ, PATH="%s:%s" % (bad, os.environ["PATH"]))
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, env=env)
    assert r.returncode == 0 and "AMSET_OK deform create" in r.stdout, (r.stdout, r.stderr)
