#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step6_kappa.py —— kl-dft-cpu S6_kappa：BTE 晶格热导率（共享 fit-fc-thermal 引擎）。

取 step5_fc 的 fc2/fc3 → 三声子 BTE（RTA/LBTE，自动 q 网格收敛 + 逐档截断三判据选择），
写 kappa_summary.json（marker: KAPPA_DONE）。计算节点跑 kappa_driver.py。

注意：换共享引擎后，旧引擎的 SOLVER=shengbte/fourphonon 尚未并入共享 kappa_driver，
后续单独补（与 kl-mlff 的路线 C 一致）。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# 共享引擎（_common/fcfit）+ stepconf/dim_common（_common/opt）+ _common
_base = Path(__file__).resolve().parent.parent
for _d in (_base / "fcfit", _base / "opt", _base):
    if _d.is_dir() and str(_d) not in sys.path:
        sys.path.append(str(_d))

import gen_step2_kappa as g2  # noqa: E402

OUTDIR = "step6_kappa"
STEP = "step6_kappa"
FIT_DIR = "step5_fc"


def main():
    g2.main(outdir=OUTDIR, fit_dir=FIT_DIR, step=STEP)


if __name__ == "__main__":
    main()
