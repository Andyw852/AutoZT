#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从技能定义自动导出论文/文档用的参数表（不手抄）。

数据来源（都在 skill/<技能>/ 下）：
  1) skill.yaml             —— 技能名/说明、步骤清单（label），以及 schema 2 的
                              io_schema.params 声明（name/where/default/values/desc）
  2) templates/**/step.conf —— 技能级与步骤级的**出厂默认值**（用户不写就用它）
  3) *.py 的 SPEC 字典      —— step.conf 里没有、但 gen 脚本内建默认值的键

用法（仓库根）：
    python3 scripts/export_params.py                       # Markdown 打到 stdout
    python3 scripts/export_params.py --out docs/PARAMETERS.md
    python3 scripts/export_params.py --csv tmp/params.csv  # 机器可读（论文附录/表格脚本）

本脚本只读仓库文件，不连超算、不跑 gen。改了 step.conf 后重跑即可刷新，
不要手工编辑生成出来的 Markdown。
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys

try:
    import yaml
except ImportError:                                   # pragma: no cover
    sys.exit("需要 pyyaml：pip install pyyaml")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKILL_DIR = os.path.join(ROOT, "skill")

KV = re.compile(r'^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*(?:#\s*(.*))?$')
SPECDEF = re.compile(r'"([A-Z][A-Z0-9_]*)"\s*:\s*\(\s*([^,()]+?)\s*,\s*"([a-z]+)"')
SKIP_DIRS = {"__pycache__"}


def git_rev() -> str:
    """当前提交（工作区脏则加 +dirty 后缀）—— 供论文溯源引用。"""
    try:
        out = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, check=True)
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return out.stdout.strip() + ("+dirty" if dirty else "")
    except Exception:                                  # noqa: BLE001
        return "(no git)"


def step_confs(skill: str):
    """[(相对路径, {键: (值, 备注)})]，技能级 templates/step.conf 排最前。"""
    base = os.path.join(SKILL_DIR, skill)
    found = []
    for dp, dn, fn in os.walk(base):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        for f in sorted(fn):
            if f != "step.conf":
                continue
            p = os.path.join(dp, f)
            vals = {}
            for line in open(p, encoding="utf-8", errors="replace"):
                m = KV.match(line.rstrip())
                if m and m.group(1):
                    vals[m.group(1)] = (m.group(2).strip(), (m.group(3) or "").strip())
            found.append((os.path.relpath(p, base), vals))
    found.sort(key=lambda x: (x[0].count(os.sep), x[0]))
    return found


def spec_defaults(skill: str):
    """{键: (默认值, 类型, 来源文件)}，来自 gen 脚本的 SPEC 字典。"""
    base = os.path.join(SKILL_DIR, skill)
    out = {}
    for dp, dn, fn in os.walk(base):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        for f in sorted(fn):
            if not f.endswith(".py"):
                continue
            p = os.path.join(dp, f)
            txt = open(p, encoding="utf-8", errors="replace").read()
            for m in SPECDEF.finditer(txt):
                out.setdefault(m.group(1), (m.group(2).strip(), m.group(3),
                                            os.path.relpath(p, base)))
    return out


def skills():
    return [s for s in sorted(os.listdir(SKILL_DIR))
            if os.path.isdir(os.path.join(SKILL_DIR, s)) and not s.startswith("_")]


def collect():
    rows = []            # (技能, 来源, 参数, 默认值, 类型, 备注)
    declared = []        # (技能, 参数, 声明默认, 取值, 位置, 说明)
    overview = []        # (技能, 说明, 步骤数, 步骤标签)
    for sk in skills():
        sp = os.path.join(SKILL_DIR, sk, "skill.yaml")
        if not os.path.isfile(sp):
            continue
        spec = yaml.safe_load(open(sp, encoding="utf-8")) or {}
        steps = [s.get("label") or s.get("name") for s in (spec.get("steps") or [])]
        overview.append((sk, str(spec.get("desc") or "").strip(), len(steps),
                         ", ".join(str(x) for x in steps)))
        seen = set()
        for rel, vals in step_confs(sk):
            is_step = os.path.dirname(rel) not in ("", "templates")
            for key, (val, note) in sorted(vals.items()):
                rows.append((sk, rel + ("（步骤级）" if is_step else "（技能级）"),
                             key, val or "(空)", "", note))
                seen.add(key)
        for key, (val, typ, src) in sorted(spec_defaults(sk).items()):
            if key in seen:
                continue
            rows.append((sk, src, key, val.strip('"'), typ, "SPEC 内建默认"))
        for prm in ((spec.get("io_schema") or {}).get("params") or []):
            declared.append((sk, prm.get("name", ""), prm.get("default", ""),
                             ", ".join(str(v) for v in prm["values"]) if prm.get("values") else "",
                             prm.get("where", ""), (prm.get("desc") or "").replace("|", "/")))
    return overview, rows, declared


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        cells = [str(c).replace("|", "/").replace("\n", " ") or "-" for c in r]
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description="导出技能参数表")
    ap.add_argument("--out", help="写到文件（Markdown）")
    ap.add_argument("--csv", help="同时导出扁平 CSV（技能,来源,参数,默认值,类型,备注）")
    args = ap.parse_args()

    overview, rows, declared = collect()
    md = []
    md.append("# 技能参数表（自动生成）")
    md.append("")
    md.append("> 由 `python3 scripts/export_params.py --out docs/PARAMETERS.md` 生成；"
              "来源：`skill/*/skill.yaml`、`skill/*/templates/**/step.conf`、gen 脚本的 SPEC。")
    md.append("> 生成时的仓库版本：`%s`。请勿手工编辑本文件，改 step.conf 后重跑即可。"
              % git_rev())
    md.append("")
    md.append("## 1. 技能与步骤")
    md.append("")
    md.append(md_table(["技能", "说明", "步骤数", "步骤"], overview))
    md.append("")
    md.append("## 2. io_schema 声明的参数（schema 2，机器可读）")
    md.append("")
    md.append(md_table(["技能", "参数", "声明默认值", "取值域", "位置", "说明"], declared))
    md.append("")
    md.append("## 3. 出厂默认值总表（step.conf 与 SPEC 内建）")
    md.append("")
    md.append(md_table(["技能", "来源", "参数", "默认值", "类型", "备注"], rows))
    md.append("")
    md.append("统计：%d 个技能，%d 条 step.conf/SPEC 默认值，%d 条声明参数。"
              % (len(overview), len(rows), len(declared)))
    text = "\n".join(md) + "\n"

    if args.out:
        path = args.out if os.path.isabs(args.out) else os.path.join(ROOT, args.out)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        print("已写出 %s（%d 行）" % (os.path.relpath(path, ROOT), text.count("\n")))
    else:
        sys.stdout.write(text)
    if args.csv:
        path = args.csv if os.path.isabs(args.csv) else os.path.join(ROOT, args.csv)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["skill", "source", "param", "default", "type", "note"])
            w.writerows(rows)
        print("已写出 %s（%d 条）" % (os.path.relpath(path, ROOT), len(rows)))


if __name__ == "__main__":
    main()
