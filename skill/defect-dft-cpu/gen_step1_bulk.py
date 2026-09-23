#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""defect-dft-cpu step1：按输入及显式 SUPERCELL 生成体相，SOC/泛函可配，弛豫（LVHAR=开，供静电势对齐）。"""
import sys, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import defects_common as D

STEP = "step1_bulk"

def main():
    conf = D.load_stepconf()
    if not os.path.exists("POSCAR"):
        raise SystemExit("[错误] 材料目录缺 POSCAR（请放入已弛豫原胞 POSCAR_B）")
    prim = D.parse_poscar("POSCAR")
    n = tuple(int(x) for x in conf.get("SUPERCELL", "1 1 1").split())
    if len(n) != 3 or any(x < 1 for x in n):
        raise SystemExit("[错误] SUPERCELL 必须为三个正整数")
    sc = D.supercell(prim, n)
    print("[计划] 输入%d原子 × %dx%dx%d -> %d原子；尺寸收敛未由原子数保证"
          % (len(prim["atoms"]), *n, len(sc["atoms"])))
    D.build_job(STEP, sc, conf, {
        "SYSTEM": "defect-dft-cpu step1_bulk %dx%dx%d" % n,
        "IBRION": "2", "ISIF": "2", "NSW": "60",
        "EDIFFG_LINE": "EDIFFG = -0.02",
        "LVHAR_LINE": "LVHAR = .TRUE.",
    }, "def_%dx%dx%d_bulk" % n)
    print("[OK] 生成 %s/（完美超胞 %d 原子）" % (STEP, len(sc["atoms"])))

if __name__ == "__main__":
    main()
