#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare an AutoZT 3090 job that samples thermal structures with the seed NEP."""
import shlex
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc
import stepconf

STEP = "step1_thermal_sample"
OUTDIR = STEP
SPEC = {
    "GPUMD_BIN": ("/home/<user>/gpumd/src/gpumd", "str"),
    "NEP_BIN": ("/home/<user>/gpumd/src/nep", "str"),
    "CONDA_SH": (gc.DEFAULT_CONDA_SH, "str"),
    "CONDA_ENV": (gc.DEFAULT_CONDA_ENV, "str"),
    "SEED_NEP": ("nep_Mn2In2Se5_gen1000.txt", "str"),
    "TEMPERATURES": ("300,500,650,800", "str"),
    "COUNTS": ("13,12,12,13", "str"),
    "EQUIL_STEPS": (30000, "int"),
    "PRODUCTION_STEPS": (100000, "int"),
    "SAMPLE_INTERVAL": (5000, "int"),
    "TIME_STEP_FS": (1.0, "float"),
    "TAU_FS": (100.0, "float"),
    "GPU_DEVICE": (1, "int"),
    "SAMPLE_SUBDIR": ("sampled", "str"),
    "REFERENCE_TRAIN": ("", "str"),
    "REFERENCE_TEST": ("", "str"),
    "DEDUP_RMSD": (0.05, "float"),
}


def main():
    cwd = Path.cwd()
    conf = stepconf.load(SPEC, None)
    step_name = str(conf["STEP"] or STEP)
    out = cwd / step_name
    out.mkdir(exist_ok=True)
    poscar = cwd / "POSCAR"
    seed = cwd / str(conf["SEED_NEP"])
    if not poscar.is_file() or not seed.is_file():
        sys.exit("[ERROR] 缺材料 POSCAR 或采样势 %s" % seed.name)
    shutil.copy2(poscar, out / "POSCAR")
    shutil.copy2(seed, out / seed.name)
    helper = Path(__file__).with_name("thermal_md_sampler.py")
    if not helper.is_file():
        sys.exit("[ERROR] gen_need 缺 thermal_md_sampler.py")
    shutil.copy2(helper, out / helper.name)
    reference_files = []
    for key in ("REFERENCE_TRAIN", "REFERENCE_TEST"):
        value = str(conf[key] or "").strip()
        if not value:
            continue
        source = Path(value)
        if not source.is_absolute():
            source = cwd / source
        if not source.is_file():
            sys.exit("[ERROR] 去重参考数据不存在：%s" % source)
        target = out / source.name
        shutil.copy2(source, target)
        reference_files.append(target.name)
    gc.save_params(out, dict(conf))
    args = [
        "python thermal_md_sampler.py",
        "--poscar POSCAR",
        "--potential %s" % shlex.quote(seed.name),
        "--gpumd %s" % shlex.quote(str(conf["GPUMD_BIN"])),
        "--outdir %s" % shlex.quote(str(conf["SAMPLE_SUBDIR"])),
        "--temperatures %s" % shlex.quote(str(conf["TEMPERATURES"])),
        "--counts %s" % shlex.quote(str(conf["COUNTS"])),
        "--equil-steps %s" % int(conf["EQUIL_STEPS"]),
        "--production-steps %s" % int(conf["PRODUCTION_STEPS"]),
        "--sample-interval %s" % int(conf["SAMPLE_INTERVAL"]),
        "--time-step-fs %s" % float(conf["TIME_STEP_FS"]),
        "--tau-fs %s" % float(conf["TAU_FS"]),
        "--dedup-rmsd %s" % float(conf["DEDUP_RMSD"]),
    ]
    for ref_name in reference_files:
        args.append("--reference-xyz %s" % shlex.quote(ref_name))
    sample_subdir = shlex.quote(str(conf["SAMPLE_SUBDIR"]))
    cmd = " ".join(args) + " && cp %s/thermal_manifest.json thermal_manifest.json" % sample_subdir
    tpl = gc.resolve_submit(Path(__file__).resolve().parent, "submit_gpumd")
    gc.write_submit(tpl, out / "submit.sh", {
        "JOBNAME": gc.new_jobname(cwd, "Tsample"),
        "CONDA_SH": conf["CONDA_SH"] or gc.DEFAULT_CONDA_SH,
        "CONDA_ENV": conf["CONDA_ENV"] or gc.DEFAULT_CONDA_ENV,
        "GPUMD_CMD": cmd,
        "LOG": "sampler.log",
    })
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    # fakeslurm pins jobs to logical GPU 0, which maps to the host's faulty
    # physical GPU 4 on this node. Pin this sampling job to healthy GPU 1.
    submit = out / "submit.sh"
    submit_text = submit.read_text(encoding="utf-8")
    omp_line = "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}"
    if omp_line not in submit_text:
        sys.exit("[ERROR] submit_gpumd.tpl 缺 OMP_NUM_THREADS 行，无法设置健康 GPU")
    submit_text = submit_text.replace(
        omp_line, omp_line + "\nexport CUDA_VISIBLE_DEVICES=%d" % int(conf["GPU_DEVICE"]), 1)
    submit.write_text(submit_text, encoding="utf-8", newline="\n")
    print("[DONE] %s：50 帧热态采样作业输入已生成" % step_name)


if __name__ == "__main__":
    main()
