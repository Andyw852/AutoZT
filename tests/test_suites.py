"""把 tests/ 下的独立自测脚本纳入 pytest。

不带标记的用例只在本机跑（不碰集群）；dropin 需要真集群，默认跳过：
    pytest -q            # 单元 + 本机自测
    pytest -q -m cluster # 端到端（需要已配置好的集群）

★ 2026-09-16 修正：本文件原来去找 tests/test_v1_<名字>.py，那些文件**并不存在**
（套件实际叫 tests/suite_<名字>.py），于是"测试全过"是假的——什么都没跑。
现在直接按 suite_*.py 枚举，并显式维护两张清单。
没有 pytest 时也可以直接跑：for t in tests/suite_*.py; do python3 $t || echo FAIL $t; done

★ 2026-09-22 教训：**不要直接 `python3 tests/test_suites.py`** —— 本模块只有
@pytest.mark.parametrize 用例、没有 __main__，直接跑什么都不做却退出 0，会被误当成"全绿"
（本仓库前面几处"test_suites RC=0"就是这么来的）。正确跑法：
    python3 -m pytest tests/test_suites.py -q            # 本机套件
    python3 -m pytest tests/test_suites.py -q -m cluster # 端到端（需集群）
"""
import glob
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 本机跑（不碰集群）的套件；名字 = tests/suite_<名字>.py
LOCAL = ["agentgate", "session", "autodeps", "history", "skillspec",
         "correct_cli", "prov", "uniform", "supercell", "templates",
         "ke_common", "asset_lookup", "poscar_sync", "io_schema", "s6_marker",
         "symmetry_audit", "structure_health", "template_drift", "stepconf"]
CLUSTER = ["dropin"]
# 独立的 pytest 测试模块（不是 suite_*.py 形式），也要纳入回归
LOCAL_PYTEST = ["test_imag_policy"]


# 干净克隆里没有 tmp/（gitignore），依赖本机测试配置的套件自动跳过
NEEDS_TMP = {"autodeps": "tmp/tf_jzzn_smoke.yaml",
             "uniform": "tmp/tf_jzzn_si_all.yaml"}
# 需要本机私有配置（setting/<集群>/templates 等未入库文件）的套件
NEEDS_LOCAL = {"templates": "setting/jzzn/templates"}


def _suite_path(name):
    p = os.path.join(ROOT, "tests", "suite_%s.py" % name)
    if not os.path.isfile(p):
        # 兜底：不在显式清单里的套件也按 suite_*.py 找
        hits = [x for x in glob.glob(os.path.join(ROOT, "tests", "suite_*.py"))
                if os.path.basename(x) == "suite_%s.py" % name]
        p = hits[0] if hits else p
    return p


def _run(name):
    p = subprocess.run([sys.executable, _suite_path(name)],
                       capture_output=True, text=True, cwd=ROOT, timeout=1800)
    print(((p.stdout or "") + (p.stderr or ""))[-1500:])
    return p.returncode


@pytest.mark.parametrize("name", LOCAL)
def test_local_suite(name):
    if not os.path.isfile(_suite_path(name)):
        pytest.skip("suite %s 不存在" % name)
    need = NEEDS_TMP.get(name) or NEEDS_LOCAL.get(name)
    if need and not os.path.exists(os.path.join(ROOT, need)):
        pytest.skip("clean clone: %s is not shipped" % need)
    assert _run(name) == 0, "suite %s failed" % name


@pytest.mark.parametrize("name", LOCAL_PYTEST)
def test_local_pytest_module(name):
    p = os.path.join(ROOT, "tests", "%s.py" % name)
    if not os.path.isfile(p):
        pytest.skip("pytest module %s 不存在" % name)
    r = subprocess.run([sys.executable, "-m", "pytest", p, "-q", "-p", "no:cacheprovider"],
                       capture_output=True, text=True, cwd=ROOT, timeout=1800)
    print(((r.stdout or "") + (r.stderr or ""))[-1200:])
    assert r.returncode == 0, "pytest module %s failed" % name


@pytest.mark.cluster
@pytest.mark.parametrize("name", CLUSTER)
def test_cluster_suite(name):
    assert _run(name) == 0, "suite %s failed" % name
