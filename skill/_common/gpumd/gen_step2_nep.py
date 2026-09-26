#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_nep.py —— NEP（微调）作业准备（step2_nep）。

把引擎脚本 + 参数 + submit.sh 准备好；作业里先 OUTCAR→extxyz，再 micro-tune/训练 NEP。
判据：nep_summary.json 的 "NEP_DONE": true（力 RMSE ≤ FORCE_RMSE_MAX）。
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc
import stepconf

OUTDIR = "step2_nep"
STEP = "step2_nep"

SPEC = {
    "GPUMD_BIN": ("/home/<user>/gpumd/src/gpumd", "str"),
    "NEP_BIN": ("/home/<user>/gpumd/src/nep", "str"),
    "CONDA_SH": (gc.DEFAULT_CONDA_SH, "str"),
    "CONDA_ENV": (gc.DEFAULT_CONDA_ENV, "str"),
    "DATA_DIR": ("", "str"),
    "DATA_GLOB": ("POSCAR-*/OUTCAR", "str"),
    "TRAIN_FRAC": (0.9, "float"),
    "SPLIT_SEED": (20260925, "int"),
    # --- 微调基座（GPUMD 自带 nep89_20250409；两个文件都要）---
    "PRETRAINED_NEP": ("", "str"),
    "PRETRAINED_RESTART": ("", "str"),
    # 1 = 只用预训练势、完全不训练：跳过 OUTCAR→xyz 与 nep 训练，
    #     直接把 PRETRAINED_NEP 落成 nep.txt，并用材料 POSCAR 作 MD 参考结构。
    "PRETRAINED_ONLY": (0, "int"),
    # --- 架构：微调时**必须**与基座一致（不可改）---
    "NEP_VERSION": (4, "int"),
    "NEP_ZBL": (2, "int"),
    "CUTOFF": ("6 5", "str"),
    "N_MAX": ("4 4", "str"),
    "BASIS_SIZE": ("8 8", "str"),
    "L_MAX": ("4 2 1", "str"),
    "NEURON": (80, "int"),
    # --- 可调 ---
    "LAMBDA_1": (0.0, "float"),
    "LAMBDA_E": (1.0, "float"),
    "LAMBDA_F": (1.0, "float"),
    "LAMBDA_V": (0.0, "float"),
    "BATCH": (5000, "int"),
    "POPULATION": (50, "int"),
    "GENERATION": (5000, "int"),
    "SAVE_EVERY": (1000, "int"),
    "SAVE_START": (0, "int"),
    "SAVE_END": (1, "int"),
    "FORCE_RMSE_MAX": (0.05, "float"),
    # 相对力 RMSE 门槛（力 RMSE / DFT 力 RMS）；>0 时优先于绝对 FORCE_RMSE_MAX 判合格。
    "FORCE_RMSE_REL_MAX": (0.05, "float"),
    "USE_VIRIAL": (0, "int"),
    "ELEMENTS": ("Mg C", "str"),
}


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    # 上一步的参数作兜底：本步 step.conf 未显式设置的键（取值仍等于 SPEC 默认）才继承。
    # ★ StepConf 是**只读 dict**（没有 setdefault/__setitem__，见 stepconf.py 注释），
    #   所以这里显式合并——原写法 conf.setdefault(...) 会 AttributeError，从未真正跑过。
    prev = cwd / "step1_struct" / gc.PARAMS_JSON
    params = dict(conf)
    if prev.is_file():
        base = gc.load_params(cwd / "step1_struct")
        for key, (default, _typ) in SPEC.items():
            if key in base and default is not None and params.get(key) == default:
                params[key] = base[key]
    params["ELEMENTS"] = [x for x in str(conf["ELEMENTS"]).split() if x]
    for key in ("CUTOFF", "N_MAX", "BASIS_SIZE", "L_MAX"):
        params[key] = [int(x) for x in str(conf[key]).split()]
    gc.save_params(out, params)
    if int(conf.get("PRETRAINED_ONLY") or 0):
        ref = cwd / "step1_struct" / "POSCAR"
        if not ref.is_file():
            sys.exit("[ERROR] PRETRAINED_ONLY=1 需要 step1_struct/POSCAR，缺 %s" % ref)
        shutil.copyfile(str(ref), str(out / "reference_POSCAR"))
    here = Path(__file__).resolve().parent
    for f in ("vasp_to_xyz.py", "nep_train.py", "gpumd_common.py"):
        if not (here / f).is_file():
            sys.exit("[ERROR] 缺 %s（gen_need 里漏了？）" % f)
        shutil.copyfile(str(here / f), str(out / f))
    cmd = "python vasp_to_xyz.py --step-dir . --use-virial %d && python nep_train.py --step-dir ."
    cmd = cmd % int(conf["USE_VIRIAL"] or 0)
    tpl = gc.resolve_submit(here, "submit_nep")
    gc.write_submit(tpl, out / "submit.sh",
                    {"JOBNAME": gc.new_jobname(cwd, "S2nep"),
                     "CONDA_SH": conf["CONDA_SH"] or gc.DEFAULT_CONDA_SH,
                     "CONDA_ENV": conf["CONDA_ENV"] or gc.DEFAULT_CONDA_ENV,
                     "NEP_CMD": cmd, "LOG": "train.log"})
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    mode = ("仅预训练、不训练" if int(conf.get("PRETRAINED_ONLY") or 0)
            else ("微调基座=%s" % conf["PRETRAINED_NEP"] if conf["PRETRAINED_NEP"]
                  else "无（从头训练）"))
    print("[DONE] %s：submit.sh 就绪（%s）" % (OUTDIR, mode))


if __name__ == "__main__":
    main()
