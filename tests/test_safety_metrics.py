"""安全指标矩阵（CI 可跑，不需要集群）：拦截率 100%、误操作 0、只读不被误拦、输出确定。"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "scripts", "safety_metrics.py")


def _report():
    p = subprocess.run([sys.executable, SCRIPT, "--json"], capture_output=True,
                       text=True, cwd=ROOT, timeout=900)
    assert p.returncode == 0, p.stderr[-500:]
    return json.loads(p.stdout)


def test_every_hazardous_action_is_intercepted():
    rep = _report()
    assert rep["interception"]["interception_rate"] == 100.0
    assert rep["interception"]["executed"] == 0


def test_every_denial_offers_an_approval_path():
    rep = _report()
    assert rep["interception"]["every_denial_offers_approval_path"] is True


def test_readonly_actions_still_work_for_agents():
    rep = _report()
    assert rep["readonly"]["all_allowed"] is True


def test_outputs_are_deterministic():
    rep = _report()
    assert rep["determinism"]["determinism_rate"] == 100.0
