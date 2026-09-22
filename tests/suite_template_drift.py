# -*- coding: utf-8 -*-
"""自测：init 阶段模板占位符漂移检查（autozt.ops._template_drift）。

判据：项目级模板若是技能模板的旧拷贝（缺了 {{占位符}}），init 时就该 WARN ——
否则要到 gen 阶段才被 render_tpl/require 拦下，那时前序作业多已提交。
"""
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from autozt.ops import _template_drift  # noqa: E402

FAILS = []


def ok(cond, msg):
    print(("  . " if cond else "  x ") + msg)
    if not cond:
        FAILS.append(msg)


def _w(p, t):
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(t)


def main():
    d = tempfile.mkdtemp(prefix="tmpl_drift_")
    try:
        s = os.path.join(d, "skill.tpl")
        p = os.path.join(d, "proj.tpl")
        _w(s, "A={{A}}\nB={{B}}\nC={{C}}\n# 说明里的 {{DOC}} 也算\n")
        _w(p, "A={{A}}\n")
        ok(_template_drift(s, p) == {"B", "C", "DOC"},
           "旧拷贝缺 3 个占位符 -> 完整报出：%s" % sorted(_template_drift(s, p)))
        _w(p, "A={{A}}\nB={{B}}\nC={{C}}\n# 说明里的 {{DOC}} 也算\n")
        ok(_template_drift(s, p) == set(), "一致 -> 无漂移")
        _w(p, "A={{A}}\nB={{B}}\nC={{C}}\nDOC={{DOC}}\nEXTRA={{EXTRA}}\n")
        ok(_template_drift(s, p) == set(), "项目副本多出的占位符不算漂移")
        ok(_template_drift(os.path.join(d, "nope.tpl"), p) == set(),
           "读不到源文件时返回空集（不拦 init）")
    finally:
        shutil.rmtree(d, ignore_errors=True)
    print()
    if FAILS:
        print("FAILED: %d 项 -> %s" % (len(FAILS), ", ".join(FAILS)))
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
