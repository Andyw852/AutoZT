#!/usr/bin/env python3
"""脱敏：按 conf 的「源词=目标词」映射替换（长词优先防子串误伤），
并把常见的站点/个人绝对路径统一成占位符（与 check_site_specific.py 的 PAT 对应）。
用法: python3 scripts/sanitize.py setting/git-sanitize.conf"""
import re
import subprocess
import sys

conf = sys.argv[1]
repls = []
for line in open(conf, encoding="utf-8"):
    line = line.strip()
    if not line or line.startswith("#") or "=" not in line:
        continue
    a, b = line.split("=", 1)
    repls.append((a, b))
repls.sort(key=lambda x: -len(x[0]))  # 长词优先

# 站点/个人绝对路径 -> 占位符；替换后的 <user> 不再被自身或 check_site_specific.py 命中
SITE_RE = re.compile(r"/(public/home|scratch|home)/[A-Za-z0-9_.-]+")


def grep_files(pattern):
    out = subprocess.run(
        ["git", "grep", "-l", "-E", pattern, "--", "."],
        capture_output=True, text=True,
    )
    return set(f for f in out.stdout.split() if f)


files = set()
if repls:
    files |= grep_files("|".join(r[0] for r in repls))
files |= grep_files(SITE_RE.pattern)          # 只含站点路径的文件也要处理

for f in sorted(files):
    try:
        t = open(f, encoding="utf-8").read()
    except (UnicodeDecodeError, IsADirectoryError, OSError):
        continue
    for a, b in repls:
        t = t.replace(a, b)
    t = SITE_RE.sub(r"/\1/<user>", t)
    open(f, "w", encoding="utf-8").write(t)
print("脱敏文件数:", len(files))
