# -*- coding: utf-8 -*-
"""自测：每个技能 skill.yaml 的 io_schema.params 声明默认值，必须与实际默认值一致。

为什么要有这条闸（2026-09-18 定）：
    schema 2 的 io_schema.params 是"技能自带机器可读参数说明"，论文/审稿人会拿它当
    接口文档看。声明一旦和 step.conf / SPEC 的真实默认值漂移（历史上 kl-dft-cpu 的
    KAPPA_MESH 声明 15 15 15 实际 auto、FIT_ENGINE 声明 phono3py 实际 auto、
    kl-mace-gpu 的 MACE_MODEL/DEVICE、unihamgnn 的 SOC、cohp 的 DENSIFY 都漂过），
    说明就是错的，比没有更糟。本套件把它变成 CI 级别的硬闸。

实际默认值的取法（"用户不写 step.conf 时拿到什么"）：
    1) skill/<技能>/templates/*/step.conf 与 templates/step.conf 里该键的**全部**出现处；
       只要**任意一处**等于声明值即算通过（同一键在不同步骤可以各不相同，例如
    kl-mace-gpu 的 DEVICE 在 step2 是 cuda、在技能级是 auto —— 声明写技能级默认值）；
    2) step.conf 里没有该键时，退回该技能 .py 里 SPEC 形式的 "KEY": (default, "type")。
    两处都找不到 → 记为"无法核对"（不算失败，但会打印，便于人工补）。
    数值按 float 值比较（0.1 == 0.10，0.002 == 2e-3）。

跑法（仓库根）：
    python3 tests/suite_io_schema.py
全绿 = exit 0；任一不一致 = exit 1（最后打印 FAIL 清单）。
"""
import os
import re
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

SKILL_DIR = os.path.join(ROOT, "skill")
KV = re.compile(r'^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*(?:#.*)?$')
SPECDEF = re.compile(r'"([A-Z][A-Z0-9_]*)"\s*:\s*\(\s*([^,()]+?)\s*,')
FAILS = []


def ok(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("  x %s" % msg)
    else:
        print("  . %s" % msg)


def same(a, b):
    """声明值 vs 实际值：先按数值比，再按字符串比（忽略大小写/引号）。"""
    sa = str(a).strip().strip('"').strip("'")
    sb = str(b).strip().strip('"').strip("'")
    try:
        return float(sa) == float(sb)
    except (TypeError, ValueError):
        return sa.lower() == sb.lower()


def conf_values(skill, key):
    """该技能所有 step.conf 里 key 的值（列表）。"""
    vals = []
    base = os.path.join(SKILL_DIR, skill)
    for dp, dn, fn in os.walk(base):
        if "__pycache__" in dp:
            continue
        for f in sorted(fn):
            if f != "step.conf":
                continue
            p = os.path.join(dp, f)
            for line in open(p, encoding="utf-8", errors="replace"):
                m = KV.match(line.split("#")[0].rstrip())
                if m and m.group(1).upper() == key.upper():
                    vals.append((m.group(2).strip(), os.path.relpath(p, base)))
    return vals


def spec_value(skill, key):
    """该技能 .py 的 SPEC 里 key 的默认值（找不到返回 (None, None)）。"""
    base = os.path.join(SKILL_DIR, skill)
    for dp, dn, fn in os.walk(base):
        if "__pycache__" in dp:
            continue
        for f in sorted(fn):
            if not f.endswith(".py"):
                continue
            txt = open(os.path.join(dp, f), encoding="utf-8", errors="replace").read()
            for m in SPECDEF.finditer(txt):
                if m.group(1) == key.upper():
                    return m.group(2).strip(), os.path.relpath(os.path.join(dp, f), base)
    return None, None


print("[1] 枚举技能并核对 io_schema.params 的声明默认值")
skills = [s for s in sorted(os.listdir(SKILL_DIR))
          if os.path.isdir(os.path.join(SKILL_DIR, s)) and not s.startswith("_")]

checked = unknown = 0
for skill in skills:
    sp = os.path.join(SKILL_DIR, skill, "skill.yaml")
    if not os.path.isfile(sp):
        continue
    spec = yaml.safe_load(open(sp, encoding="utf-8")) or {}
    for prm in ((spec.get("io_schema") or {}).get("params") or []):
        name = prm.get("name")
        if not name or "default" not in prm:
            continue
        if "step.conf" not in str(prm.get("where") or ""):
            continue                      # 项目/tf.yaml 级开关，不在本检查范围
        want = prm["default"]
        vals = conf_values(skill, name)
        if vals:
            if any(same(want, v) for v, _ in vals):
                checked += 1
                continue
            checked += 1
            src = "; ".join("%s=%r" % (f, v) for v, f in vals)
            ok(False, "%s.%s 声明 default=%r，step.conf 实际 %s" % (skill, name, want, src))
            continue
        got, src = spec_value(skill, name)
        if got is None:
            unknown += 1
            print("  ? %s.%s 无法核对（step.conf 与 SPEC 都没有该键，声明=%r）" % (skill, name, want))
            continue
        checked += 1
        ok(same(want, got), "%s.%s 声明 default=%r，SPEC 实际 %r (%s)"
           % (skill, name, want, got, src))

print("[2] 统计")
print("  核对 %d 项，无法核对 %d 项" % (checked, unknown))
ok(checked >= 20, "至少核对了 20 项（防止枚举逻辑失效导致空转）：实际 %d" % checked)

if FAILS:
    print("\nFAIL %d 项：" % len(FAILS))
    for m in FAILS:
        print("  - %s" % m)
    sys.exit(1)
print("\n全部通过。")
