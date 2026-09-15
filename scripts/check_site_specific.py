#!/usr/bin/env python3
"""列出仓库里写死的站点/个人路径，供发布前自查（CI 里跑，非阻断）。

用法：
    python scripts/check_site_specific.py            # 全仓库摘要
    python scripts/check_site_specific.py --strict   # 有任何一处在【执行代码】里就退出 1
"""
import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAT = re.compile(r"(/home/[A-Za-z0-9_.-]+|/public/home/[A-Za-z0-9_.-]+|/scratch/[A-Za-z0-9_.-]+)")
SKIP_DIRS = {".git", "__pycache__", "tmp", "output", "versions", "hfeshell_tmp"}
CODE_EXT = (".py", ".sh", ".tpl")


def scan():
    hits = {"code": [], "config": [], "docs": []}
    for base, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for fn in files:
            p = os.path.join(base, fn)
            rel = os.path.relpath(p, ROOT)
            if not fn.endswith((".py", ".sh", ".tpl", ".yaml", ".yml", ".md", ".conf")):
                continue
            try:
                txt = open(p, encoding="utf-8", errors="replace").read()
            except OSError:
                continue
            for i, line in enumerate(txt.splitlines(), 1):
                if PAT.search(line):
                    kind = ("code" if fn.endswith(CODE_EXT) and not rel.startswith("tests/")
                            else "docs" if fn.endswith(".md") else "config")
                    if rel.startswith("setting/template"):
                        kind = "docs"          # 模板里的占位示例不算
                    hits[kind].append("%s:%d" % (rel, i))
    return hits


h = scan()
print("站点相关路径统计（%d 处）：" % sum(len(v) for v in h.values()))
print("  执行代码/脚本/模板 : %d" % len(h["code"]))
print("  配置文件           : %d" % len(h["config"]))
print("  文档与模板示例     : %d" % len(h["docs"]))
for k in ("code", "config"):
    for x in h[k][:8]:
        print("   %-8s %s" % (k, x))
if "--strict" in sys.argv:
    sys.exit(1 if h["code"] else 0)
