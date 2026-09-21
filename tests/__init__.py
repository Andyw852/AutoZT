# -*- coding: utf-8 -*-
"""让 tests/ 成为常规包。

2026-09-21：nspin_norm_fix 离线测试用
    python -m unittest tests.test_nspin_norm_fix -v
运行。某些环境（如 atomate2_p_a）的 site-packages 里存在同名的 tests 常规包，
而命名空间包（无 __init__.py 的 tests/）优先级低于它，导致 unittest 解析到
site-packages/tests 而报 ModuleNotFoundError。加上本文件后，仓库根目录的
tests 成为常规包并优先于 site-packages。
"""
