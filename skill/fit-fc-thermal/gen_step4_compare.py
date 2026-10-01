#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step4_compare.py -- S4_compare (login node, run: gen).

Gathers step3_sweep/m-*/sweep_result.json into
step4_compare/fit_compare.json (+ fit_compare.csv): one row per
engine x method x neighbour shell with the fit gate, the fit residual and
kappa(300 K).  Pure bookkeeping, no physics is recomputed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# 脱离 autozt 直接运行（python skill/fit-fc-thermal/gen_*.py）时，gen_need 里
# 来自 skill/_common 的 stepconf/dim_common/thickness_2d 等不在脚本旁边——
# 追加到 sys.path 末尾兜底（autozt 推送的同目录拷贝仍优先）。
for _d in ("_common/opt", "_common"):
    _p = Path(__file__).resolve().parent.parent / _d
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.append(str(_p))
import sweep_driver

OUTDIR = "step4_compare"
SWEEP_DIR = "step3_sweep"


def main():
    cwd = Path.cwd()
    sweep = cwd / SWEEP_DIR
    if not (sweep / "sweep_plan.json").is_file():
        sys.exit("[ERROR] %s/sweep_plan.json missing -- run S3_sweep first" % SWEEP_DIR)
    doc = sweep_driver.compare(sweep, cwd / OUTDIR)
    if doc["missing"]:
        # the step is re-run by autozt once the missing variants finish
        print("[WARN] 汇总不完整：%d/%d 个变体完成" % (doc["n_done"], doc["n_variants"]))


if __name__ == "__main__":
    main()
