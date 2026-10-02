#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lineage_check.py —— 上游同源核对（V140）：各步结构是否都来自 S1 的当前 CONTCAR、派生产物是否旧于来源。

背景：S1 重新弛豫后整条 rerun 时，
  · autozt rerun 先 rm -rf 步骤目录再 gen，V122/V136 的"gen 见到旧产物 -> 下游失效"不触发
    （V140 起各 gen 改为无条件失效；在此之前已经 rerun 过的，用 --invalidate-from 补）；
  · S2.3 HSE 的 gen 从 S2.2 拷 POSCAR/WAVECAR，S2.2 没先重跑就 rerun S2.3，HSE 算的是旧结构；
  · S4 的 h5、S7.1 的形变势、S2 画图的 band_summary.json 仍是旧的"完成"。
S8/S8.4 的 gen 做同一核对（STRUCTURE_GUARD），不通过就不生成。本工具只读（除非 --invalidate-from）。

用法：
    python lineage_check.py <材料目录或项目根> [--glob "*/ke-dft-cpu"]
        列出每个材料的问题；退出码 1 = 有问题。
    python lineage_check.py <材料目录> --invalidate-from step3_uniform step7_deform ...
        已经 rerun 过这些上游：把它们下游的完成标记改名归档（*.stale-upstream-*，不删），
        之后 auto-advance 会在上游跑完后按顺序重做 S4 / S7.1 / S8 …。
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE.parent), str(HERE.parent.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


def find_materials(root, pattern=None):
    root = Path(root)
    if any((root / d).is_dir() for d in kc.S1_DIRS):
        return [root]
    cands = root.glob(pattern) if pattern else (p.parent for d in kc.S1_DIRS for p in root.rglob(d))
    return sorted({Path(p) for p in cands if any((Path(p) / d).is_dir() for d in kc.S1_DIRS)})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("path")
    ap.add_argument("--glob", default=None, help="材料目录的 glob（相对 path），默认递归找含 step1_opt 的目录")
    ap.add_argument("--invalidate-from", nargs="+", default=None, metavar="STEP",
                    help="只对单个材料目录：这些上游已经重算 -> 归档它们下游的完成标记")
    a = ap.parse_args(argv)

    if a.invalidate_from:
        mat = Path(a.path)
        if not any((mat / d).is_dir() for d in kc.S1_DIRS):
            sys.exit("[ERROR] --invalidate-from 只接受单个材料目录（里面要有 step1_opt）：%s" % mat)
        unknown = [s for s in a.invalidate_from if s not in kc.DOWNSTREAM]
        if unknown:
            sys.exit("[ERROR] 不认识的步骤 %s；可用：%s" % (unknown, ", ".join(sorted(kc.DOWNSTREAM))))
        done = []
        for step in a.invalidate_from:
            done += kc.invalidate_downstream(mat, step, "lineage_check --invalidate-from")
        print("[OK] 归档了 %d 个完成标记%s" % (len(done), "" if done else "（下游本来就没有完成标记）"))
        for step, m in done:
            print("     %s/%s" % (step, m))

    bad = 0
    mats = find_materials(a.path, a.glob)
    if not mats:
        print("没找到材料目录（含 step1_opt / step1_std_opt）：%s" % a.path)
        return 0
    for mat in mats:
        ref, probs = kc.structure_lineage(mat)
        if not probs:
            print("[OK]  %s%s" % (mat, "" if ref else "（没有 S1 CONTCAR，只核对了派生产物）"))
            continue
        bad += 1
        print("[★]   %s：%d 处" % (mat, len(probs)))
        for step, why in probs:
            print("        - %s：%s" % (step, why))
    if bad:
        print("\n%d/%d 个材料的上游不同源。结构不同的 VASP 步骤 rerun（S2.3 要等 S2.2、S2.2 要等 S2.1）；"
              "派生步骤在来源跑完后 rerun，或用 --invalidate-from 让 auto-advance 接手。" % (bad, len(mats)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
