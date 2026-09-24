#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""device-zt 回归套件入口：真正实现在 tests/device_zt/run_cases.py。

放进 suite_*.py 命名，是为了让 tests/test_suites.py 的 pytest 参数化能发现它
（device-zt 与 device-thermal 共用 skill/_common/thermal/device_common.py 的求解器，
改公共池后必须两个技能一起回归）。
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "zt_run_cases", os.path.join(HERE, "device_zt", "run_cases.py"))
_m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m)

if __name__ == "__main__":
    sys.exit(_m.main(sys.argv[1:]))
