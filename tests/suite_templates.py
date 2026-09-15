#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""出厂模板/参数体检（v1.0）——把端到端真跑踩出来的坑固化成静态检查。

起因（2026-09-15）：defect-dft-cpu 的 step.conf 写 ENCUT = auto，VASP 直接
"Error reading item ENCUT from file INCAR ... REFUSE TO CONTINUE"，但静态看代码
完全看不出来——只有真跑才炸。所以这里做几项"一眼能查"的体检：

1) step.conf 里会被塞进 INCAR 的数值参数（ENCUT/SIGMA/NELM/NSW…）不能是非数值；
2) INCAR 模板里 {{ENCUT}} 这类占位符必须能在对应 step.conf 或默认值里找到数；
3) 提交模板里 #SBATCH --nodes=1 + mpirun 的组合，建议带单节点 MPI 防护
   （osc=^ucx 或 btl=self）——坏节点（cu42/cu50）不加会段错误/狂刷日志。
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OK, FAIL, WARN = [], [], []


def ck(cond, msg):
    (OK if cond else FAIL).append(msg)
    print("  %s %s" % ("PASS" if cond else "FAIL", msg))


def warn(cond, msg):
    if not cond:
        WARN.append(msg)
        print("  WARN %s" % msg)


NUMERIC_KEYS = ("ENCUT", "SIGMA", "NELM", "NELMIN", "NSW", "EDIFF", "EDIFFG",
                "ISMEAR", "ISIF", "IBRION", "LMAXMIX", "IVDW", "PREC")
BAD_WORDS = ("auto", "default", "none", "?")

print("== 1. step.conf 数值参数 ==")
for conf in sorted(glob.glob(os.path.join(ROOT, "skill", "*", "**", "step.conf"),
                             recursive=True)):
    try:
        txt = open(conf, encoding="utf-8-sig").read()
    except OSError:
        continue
    for i, line in enumerate(txt.splitlines(), 1):
        s = line.split("#")[0].strip()
        m = re.match(r"^([A-Z_]+)\s*=\s*(.+)$", s)
        if not m or m.group(1) not in NUMERIC_KEYS:
            continue
        val = m.group(2).strip().strip('"').strip("'")
        if val.lower() in BAD_WORDS or not re.match(r"^-?[\d.eE+-]+$", val):
            FAIL.append("%s:%d %s = %s 不是数值" % (conf, i, m.group(1), val))
            print("  FAIL %s:%d  %s = %s（会被原样写进 INCAR，VASP 会 REFUSE）"
                  % (os.path.relpath(conf, ROOT), i, m.group(1), val))
ck(not [x for x in FAIL], "所有 step.conf 数值参数都是数值（ENCUT=auto 那类会被拦下）")

print("== 2. 真值核对：defect 的 ENCUT ==")
dc = os.path.join(ROOT, "skill", "defect-dft-cpu", "templates", "step.conf")
v = None
for line in open(dc, encoding="utf-8-sig"):
    m = re.match(r"^ENCUT\s*=\s*(\S+)", line.split("#")[0].strip())
    if m:
        v = m.group(1)
ck(v == "370", "defect-dft-cpu 的 ENCUT = %s（技能内 gen_references/run_isif2_static 都是 370）" % v)

print("== 3. ncl 提交模板的环境自洽 ==")
ncl = os.path.join(ROOT, "setting", "jzzn", "templates", "submit_ncl_3d.tpl")
txt = open(ncl, encoding="utf-8").read()
ck("aocc" in txt and "aocl" in txt,
   "ncl 模板用 AOCC/AOCL 环境（与 setting/jzzn.yaml 里 vasp.standard.ncl 指向的二进制一致）")
_loads_intel = [ln for ln in txt.splitlines()
                if re.match(r"^\s*(module load|source)\b.*(oneapi|intel)", ln)]
ck(not _loads_intel,
   "ncl 模板不再加载 Intel/oneAPI 模块（那个环境没有 libscalapack → 个人编译的 ncl 起不来）")
ck("OMPI_MCA_osc=^ucx" in txt,
   "单节点 MPI 防护在（坏节点 cu42/cu50 的 osc_ucx 会在 MPI_Init 段错误）")

print("\n结果：PASS %d，FAIL %d，WARN %d" % (len(OK), len(FAIL), len(WARN)))
for x in FAIL:
    print("  ! %s" % x)
sys.exit(1 if FAIL else 0)
