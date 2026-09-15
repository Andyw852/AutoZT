#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""盘点运行期中文串（i18n 迁移进度）。

用法：
    python scripts/i18n_report.py            # 按文件统计
    python scripts/i18n_report.py --list 20  # 另附前 20 条待翻译字符串
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PAT = re.compile(r'"[^"\n]*[\u4e00-\u9fff][^"\n]*"')
SKIP = {"i18n.py"}


def main():
    rows, samples = [], []
    for p in sorted((ROOT / "phonoagent").glob("*.py")):
        if p.name in SKIP:
            continue
        txt = p.read_text(encoding="utf-8")
        hits = PAT.findall(txt)
        # 只统计真正的 t(...) 调用：前面不是标识符字符（避免把 print("中文... 也算进来）
        translated = len(re.findall(
            r'(?<![A-Za-z0-9_])t\(\s*"[^"\n]*[\u4e00-\u9fff]', txt))
        if hits:
            rows.append((p.name, len(hits), translated))
            for h in hits[:3]:
                samples.append("%s: %s" % (p.name, h.strip('"')))
    total = sum(r[1] for r in rows)
    done = sum(r[2] for r in rows)
    print("运行期中文字面量：%d 条（已接 t() 的调用点：%d）" % (total, done))
    print("迁移进度：%.1f%%" % (100.0 * done / total if total else 0.0))
    print()
    print("%-22s %6s %8s" % ("file", "zh", "t()"))
    for name, n, tr in sorted(rows, key=lambda r: -r[1]):
        print("%-22s %6d %8d" % (name, n, tr))
    if "--list" in sys.argv:
        k = 20
        try:
            k = int(sys.argv[sys.argv.index("--list") + 1])
        except (IndexError, ValueError):
            pass
        print("\n前 %d 条待翻译：" % k)
        for s in samples[:k]:
            print("   ", s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
