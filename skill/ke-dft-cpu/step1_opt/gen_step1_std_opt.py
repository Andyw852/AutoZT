#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ke-dft-cpu 技能 step1：与 elastic-dft-cpu 同策略（原来两个文件字节相同）。

放置：skill/ke-dft-cpu/step1_opt/gen_step1_std_opt.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relax_common as R

# [patch_rerun_cascade V153] S1 目录是新建的（rerun）-> 下游 S2/S3/S3b/S3c/S5/S6/S7（及 zt 的 kl 分支）整目录
#   归档，等 S1 跑完后按新结构重新生成。S1 原地重新生成（已有 OUTCAR）不动下游。必须在 R.run 写文件之前。
try:
    for _up in Path(__file__).resolve().parents[:2]:
        if (_up / "ke_common.py").is_file() and str(_up) not in sys.path:
            sys.path.insert(0, str(_up))
    import ke_common as _kc
except ImportError:
    print("[WARN] 没有 ke_common.py：S1 重算不会让下游自动重排（skill.yaml 的 gen_need 漏了它？）")
else:
    _kc.cascade_on_rerun(Path.cwd(), "step1_opt")

R.run(
    OUTDIR_SINGLE="step1_opt",
    SCRIPT_NAME="gen_step1_std_opt.py",
    NEXT_STEP="gen_step2_static.py",
    STAGE_MODE="in_job",
    CELL_POLICY="primitive",        # 默认原胞（四技能统一）；需惯用晶轴取向(C_ij/输运张量)时，在该项目 step.conf 设 CELL_POLICY=standard
    STD_CELL="primitive_standard",
    VACUUM_AXIS_POLICY="rotate",
    MOL_BRANCH=True,                # step1 支持 0D；AMSET/形变势步骤用 require_dim 拦住
)
