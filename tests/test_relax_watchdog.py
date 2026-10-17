# -*- coding: utf-8 -*-
"""run_relax.sh 看门狗（relax_common.RUN_RELAX_HELPERS 的 _watchdog）实跑：STALL_MIN / FREEZE_MIN 的开关语义。

  · STALL_MIN=0（关）+ FREEZE_MIN>0：以前 STALL_MIN=0 直接 return，加了 FREEZE_MIN 之后循环照跑，
    stall 分支拿 0 当阈值 -> 第一分钟就杀。现在只按 FREEZE_MIN 杀。
  · step.conf 的 FREEZE_MINUTES 落到 FREEZE_MIN。
sleep 被换成空函数，"分钟"即循环轮数。
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
OPT = REPO / "skill" / "_common" / "opt"          # autozt 推到超算的是公共池这一份（仓库根的同名文件是旧副本）


def _rc():
    import importlib.util
    if str(OPT) not in sys.path:
        sys.path.insert(0, str(OPT))
    spec = importlib.util.spec_from_file_location("relax_common_opt", str(OPT / "relax_common.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _helpers():
    rc = _rc()
    src = rc.RUN_RELAX_HELPERS
    i = src.index("_cpu_seconds () {")
    j = src.index("\n}\n", src.index("_watchdog () {")) + 3
    return src[i:j]


def _run(stall, freeze):
    d = tempfile.mkdtemp()
    try:
        Path(d, "OUTCAR").write_text("x\n")
        Path(d, "OSZICAR").write_text("x\n")
        script = ("set -u\nSTALL_MIN=%d\nFREEZE_MIN=%d\n%s\n"
                  "command sleep 30 & P=$!\n"
                  "n=0; sleep () { n=$((n+1)); [ \"$n\" -gt 50 ] && kill -KILL \"$P\" 2>/dev/null; :; }\n"
                  "_watchdog \"$P\"\necho \"ROUNDS=$n\"\n" % (stall, freeze, _helpers()))
        r = subprocess.run(["bash", "-c", script], cwd=d, capture_output=True, text=True, timeout=60)
        err = Path(d, "queue.err").read_text() if Path(d, "queue.err").is_file() else ""
        rounds = int(r.stdout.strip().split("ROUNDS=")[-1])
        return rounds, err
    finally:
        shutil.rmtree(d, ignore_errors=True)


@pytest.mark.skipif(not shutil.which("bash") or not shutil.which("ps"), reason="需要 bash/ps")
def test_stall_off_freeze_on_kills_at_freeze_min():
    rounds, err = _run(0, 5)
    assert "同时冻结" in err and "判定卡死" not in err, err
    assert rounds >= 5                                   # 第 5 轮才判；以前 STALL_MIN=0 第 1 轮就杀


@pytest.mark.skipif(not shutil.which("bash") or not shutil.which("ps"), reason="需要 bash/ps")
def test_both_off_returns_immediately():
    rounds, err = _run(0, 0)
    assert rounds == 0 and err == ""


def test_freeze_minutes_from_step_conf(monkeypatch):
    rc = _rc()
    monkeypatch.setattr(rc, "STEP_PARAMS", {"FREEZE_MINUTES": 0, "STALL_MINUTES": 240})
    monkeypatch.setattr(rc, "FREEZE_MIN", rc.FREEZE_MIN)
    monkeypatch.setattr(rc, "STALL_MIN", rc.STALL_MIN)
    rc.apply_step_params()
    assert rc.FREEZE_MIN == 0 and rc.STALL_MIN == 240
    assert "FREEZE_MINUTES" in rc.CONF_SPEC


def test_freeze_default_not_below_stall():
    # 冻结判据不看 CPU、比 STALL 宽：默认值小于 STALL_MIN 会把 09-17 实测的长时间文件系统停顿重新误杀
    rc = _rc()
    assert rc.FREEZE_MIN == 240 and rc.FREEZE_MIN >= rc.STALL_MIN
