#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_struct.py —— 结构/数据检查（step1_struct，run: gen，不提交作业）。

检查材料 POSCAR 与训练数据目录，落盘 gpumd_params.json，写 struct_summary.json。
真正的 OUTCAR→extxyz 转换在 S2 的作业里做（数据可能要先从别的集群拷过来）。
"""
import glob
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc
import stepconf

OUTDIR = "step1_struct"
STEP = "step1_struct"

SPEC = {
    "GPUMD_BIN": ("/home/<user>/gpumd/src/gpumd", "str"),
    "NEP_BIN": ("/home/<user>/gpumd/src/nep", "str"),
    "CONDA_SH": (gc.DEFAULT_CONDA_SH, "str"),
    "CONDA_ENV": (gc.DEFAULT_CONDA_ENV, "str"),
    "DATA_DIR": ("", "str"),
    "DATA_GLOB": ("POSCAR-*/OUTCAR", "str"),
    "TRAIN_FRAC": (0.9, "float"),
    "SPLIT_SEED": (20260925, "int"),
    "ELEMENTS": ("Mg C", "str"),
    "STRUCT_SRC": ("dft_ref", "str"),
    "IFC_KAPPA_RTA_DIR": ("", "str"),
    "KAPPA_SG_AVG": (0.0, "float"),
}


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    poscar = cwd / "POSCAR"
    if not poscar.is_file():
        sys.exit("[ERROR] 材料目录下没有 POSCAR")
    shutil.copyfile(str(poscar), str(out / "POSCAR"))
    data_dir = conf["DATA_DIR"] or ""
    n_outcar = 0
    if data_dir:
        n_outcar = len(glob.glob(os.path.join(data_dir, conf["DATA_GLOB"])))
        if n_outcar == 0:
            sys.exit("[ERROR] %s 下没找到 %s（数据还没从 jzzn 拷过来？）"
                     % (data_dir, conf["DATA_GLOB"]))
    for k in ("GPUMD_BIN", "NEP_BIN"):
        if not os.access(conf[k], os.X_OK):
            sys.exit("[ERROR] %s=%s 不存在或不可执行" % (k, conf[k]))
    params = dict(conf)
    params["ELEMENTS"] = [x for x in str(conf["ELEMENTS"]).split() if x]
    gc.save_params(out, params)
    gc.write_method(out / gc.METHOD_FILE, FUNC="gpumd", DIM="3d",
                    MODEL="NEP(fine_tune)", DATA_DIR=data_dir)
    gc.write_summary(out, "struct_summary.json",
                     dict(STRUCT_DONE=True, n_outcar=n_outcar,
                          data_dir=data_dir, poscar=str(poscar),
                          gpumd_bin=conf["GPUMD_BIN"], nep_bin=conf["NEP_BIN"]))
    print("[DONE] %s：参数与结构就绪（%d 个 OUTCAR）" % (OUTDIR, n_outcar))


if __name__ == "__main__":
    main()
