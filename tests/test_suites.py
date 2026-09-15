"""把 tests/ 下的独立自测脚本纳入 pytest。

不带标记的用例只在本机跑（不碰集群）；dropin 需要真集群，默认跳过：
    pytest -q            # 单元 + 本机自测
    pytest -q -m cluster # 端到端（需要已配置好的集群）
"""
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL = ["agentgate", "session", "autodeps", "history", "skillspec",
         "correct_cli", "prov", "uniform", "supercell", "templates"]
CLUSTER = ["dropin"]


def _run(name):
    p = subprocess.run([sys.executable, os.path.join(ROOT, "tests", "test_v1_%s.py" % name)],
                       capture_output=True, text=True, cwd=ROOT, timeout=1800)
    print(((p.stdout or "") + (p.stderr or ""))[-1500:])
    return p.returncode


@pytest.mark.parametrize("name", LOCAL)
def test_local_suite(name):
    assert _run(name) == 0, "suite %s failed" % name


@pytest.mark.cluster
@pytest.mark.parametrize("name", CLUSTER)
def test_cluster_suite(name):
    assert _run(name) == 0, "suite %s failed" % name
