#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_kappa.py —— GPUMD MD-κ 作业准备（step3_kappa）。

从 step2_nep 接力 nep.txt / structure_reference.xyz，写 gpumd_params.json 与 submit.sh。
判据：kappa_summary.json 的 "KAPPA_DONE": true。
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc
import stepconf

OUTDIR = "step3_kappa"
STEP = "step3_kappa"

SPEC = {
    "GPUMD_BIN": ("/home/user_3090/gpumd/src/gpumd", "str"),
    "NEP_BIN": ("/home/user_3090/gpumd/src/nep", "str"),
    "CONDA_SH": (gc.DEFAULT_CONDA_SH, "str"),
    "CONDA_ENV": (gc.DEFAULT_CONDA_ENV, "str"),
    # --- MD 协议 ---
    "TIME_STEP": (1.0, "float"),
    "T_COUP": (100.0, "float"),
    "T_LIST": ("300", "str"),
    "T_FINAL": ("", "str"),          # 温度 ramp 终值；留空 = 与 T_LIST 相同（恒温）
    "EQUIL_STEPS": (200000, "int"),
    "PROD_STEPS": (2000000, "int"),
    "OUTPUT_INTERVAL": (1000, "int"),
    "DRIVING_FORCE": (0.00001, "float"),
    "CELL_REPLICATE": ("1 1 1", "str"),
    "EMD_ENABLED": (True, "bool"),
    "EMD_SAMPLE": (10, "int"),
    "EMD_CORR_STEPS": (100000, "int"),
    "EMD_OUTPUT": (100, "int"),
    "EMD_PROD_STEPS": (2000000, "int"),
    # --- 谱量子修正（可选 SHC）：compute_shc 谱热流 -> κ(ω) -> f(ω) 逐频修正 ---
    "SHC_ENABLED": (False, "bool"),
    "SHC_SAMPLE": (10, "int"),
    "SHC_NC": (200, "int"),
    "SHC_NUM_OMEGA": (200, "int"),
    "SHC_MAX_OMEGA": (50.0, "float"),
    "SHC_PROD_STEPS": (2000000, "int"),
    "DEBYE_T": (0.0, "float"),
    "IFC_KAPPA_RTA_DIR": ("", "str"),
    "KAPPA_SG_AVG": (0.0, "float"),
}


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    prev = cwd / "step2_nep"
    for f in ("nep.txt", "structure_reference.xyz"):
        src = prev / f
        if not src.is_file():
            sys.exit("[ERROR] 缺 %s（S2 是否已完成？）" % src)
        shutil.copyfile(str(src), str(out / f))
    params = dict(conf)
    params["T_LIST"] = [float(x) for x in str(conf["T_LIST"]).split()]
    params["CELL_REPLICATE"] = [int(x) for x in str(conf["CELL_REPLICATE"]).split()]
    gc.save_params(out, params)
    here = Path(__file__).resolve().parent
    for f in ("gpumd_kappa.py", "gpumd_common.py"):
        if not (here / f).is_file():
            sys.exit("[ERROR] 缺 %s（gen_need 里漏了？）" % f)
        shutil.copyfile(str(here / f), str(out / f))
    tpl = gc.resolve_submit(here, "submit_gpumd")
    gc.write_submit(tpl, out / "submit.sh",
                    {"JOBNAME": gc.new_jobname(cwd, "S3kl"),
                     "CONDA_SH": conf["CONDA_SH"] or gc.DEFAULT_CONDA_SH,
                     "CONDA_ENV": conf["CONDA_ENV"] or gc.DEFAULT_CONDA_ENV,
                     "GPUMD_CMD": "python gpumd_kappa.py --step-dir .",
                     "LOG": "gpumd.log"})
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    print("[DONE] %s：submit.sh 就绪（T=%s，Fe=%g Å^-1，EMD=%s）"
          % (OUTDIR, conf["T_LIST"], conf["DRIVING_FORCE"], conf["EMD_ENABLED"]))


if __name__ == "__main__":
    main()
