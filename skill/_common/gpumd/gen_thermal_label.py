#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fan out thermal MD proposals into reproducible single-point VASP jobs."""
import json
import hashlib
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stepconf

STEP = "step1_thermal_label"
OUTDIR = STEP
SPEC = {
    "GPUMD_BIN": ("/home/<user>/gpumd/src/gpumd", "str"),
    "NEP_BIN": ("/home/<user>/gpumd/src/nep", "str"),
    "INCAR_SOURCE": ("INCAR_thermal_source", "str"),
    "KPOINTS_SOURCE": ("KPOINTS_thermal_source", "str"),
    "POTCAR_SOURCE": ("POTCAR_thermal_source", "str"),
    "SAMPLE_STEP": ("step1_thermal_sample", "str"),
    "SAMPLE_SUBDIR": ("sampled", "str"),
    "EXPECTED_FRAMES": (50, "int"),
    "EXPECTED_TRAIN": (40, "int"),
    "EXPECTED_VALIDATION": (10, "int"),
    "LABEL_OFFSET": (0, "int"),
}


def set_incar(text, updates):
    lines = []
    seen = set()
    for line in text.splitlines():
        m = re.match(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=", line)
        key = m.group(1).upper() if m else None
        if key in updates:
            if key not in seen:
                lines.append("%-10s = %s" % (key, updates[key]))
                seen.add(key)
            continue
        lines.append(line)
    for key, value in updates.items():
        if key not in seen:
            lines.append("%-10s = %s" % (key, value))
    return "\n".join(lines) + "\n"


def poscar_species(path):
    lines = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    syms = lines[5].split()
    nums = [int(x) for x in lines[6].split()]
    if len(syms) != len(nums):
        raise ValueError("POSCAR species/count lines differ")
    return syms, nums


def main():
    cwd = Path.cwd()
    conf = stepconf.load(SPEC, None)
    step_name = str(conf["STEP"] or STEP)
    out = cwd / step_name
    out.mkdir(exist_ok=True)
    sample = cwd / str(conf["SAMPLE_STEP"]) / str(conf["SAMPLE_SUBDIR"])
    mp = sample / "thermal_manifest.json"
    if not mp.is_file():
        sys.exit("[ERROR] 缺少采样清单 %s（先完成热态 MD）" % mp)
    manifest = json.loads(mp.read_text(encoding="utf-8"))
    entries = manifest.get("frames", [])
    if len(entries) != int(conf["EXPECTED_FRAMES"]):
        sys.exit("[ERROR] 采样帧数 %d ≠ 预期 %d" %
                 (len(entries), int(conf["EXPECTED_FRAMES"])))
    train_n = sum(x.get("split") == "train" for x in entries)
    val_n = sum(x.get("split") == "validation" for x in entries)
    expected_train = int(conf["EXPECTED_TRAIN"])
    expected_validation = int(conf["EXPECTED_VALIDATION"])
    if (train_n, val_n) != (expected_train, expected_validation):
        sys.exit("[ERROR] train/validation=%d/%d，期望 %d/%d"
                 % (train_n, val_n, expected_train, expected_validation))
    sources = [Path(str(conf[k])) for k in
               ("INCAR_SOURCE", "KPOINTS_SOURCE", "POTCAR_SOURCE")]
    if not all(x.is_file() for x in sources):
        sys.exit("[ERROR] 缺少 INCAR/KPOINTS/POTCAR 来源文件：%s" %
                 ", ".join(str(x) for x in sources if not x.is_file()))

    incar_text = sources[0].read_text(encoding="utf-8", errors="replace")
    # Same electronic/magnetic/DFT+U and k mesh as the existing 51-frame dataset.
    # Output-only flags are disabled; labels are fixed-cell static energies and forces.
    incar_text = set_incar(incar_text, {
        "IBRION": "-1", "NSW": "0", "ISYM": "0",
        "LWAVE": ".FALSE.", "LCHARG": ".FALSE.",
    })
    kp_text = sources[1].read_text(encoding="utf-8", errors="replace")
    pc_text = sources[2].read_text(encoding="utf-8", errors="replace")
    template = cwd / "submit_std_3d.tpl"
    if not template.is_file():
        sys.exit("[ERROR] AutoZT 未推送 jzzn submit_std_3d.tpl")
    template_text = template.read_text(encoding="utf-8")

    n_new = n_kept = 0
    offset = int(conf["LABEL_OFFSET"])
    for i, ent in enumerate(entries, 1):
        serial = offset + i
        source_poscar = sample / ent["file"]
        if not source_poscar.is_file():
            sys.exit("[ERROR] 采样结构缺失：%s" % source_poscar)
        cfg = out / ("cfg-%03d" % serial)
        cfg.mkdir(exist_ok=True)
        stamp = cfg / ".thermal_frame.json"
        meta = {"frame_id": ent.get("id"), "temperature_K": ent.get("temperature_K"),
                "split": ent.get("split"), "source_frame_index": ent.get("source_frame_index"),
                "poscar_sha256": hashlib.sha256(source_poscar.read_bytes()).hexdigest()}
        submit = template_text.replace("{{JOBNAME}}", "Mn2In2Se5-t%03d" % serial)
        meta["dft_inputs_sha256"] = hashlib.sha256(
            (incar_text + kp_text + pc_text + submit +
             json.dumps(dict(conf.submit), sort_keys=True)).encode("utf-8")).hexdigest()
        if (cfg / "OUTCAR").is_file() and \
                "General timing and accounting informations" in \
                (cfg / "OUTCAR").read_text(errors="ignore")[-300000:]:
            if stamp.is_file() and json.loads(stamp.read_text()) == meta:
                n_kept += 1
                continue
            sys.exit("[ERROR] %s 有已完成 OUTCAR，但结构指纹标记不同/缺失；为避免把旧力标签"
                     "配给新结构，保留目录并停止。请先人工归档该 cfg 的旧结果再 retry。"
                     % cfg.name)
        shutil.copy2(source_poscar, cfg / "POSCAR")
        syms, nums = poscar_species(cfg / "POSCAR")
        if sum(nums) != int(manifest["natoms"]):
            sys.exit("[ERROR] %s 原子数变化：%s" % (cfg.name, nums))
        potcar_syms = []
        for line in pc_text.splitlines():
            if "TITEL" in line and "PAW_PBE" in line:
                potcar_syms.append(line.split("PAW_PBE", 1)[1].strip().split()[0])
        potcar_base = [x.split("_")[0] for x in potcar_syms]
        if potcar_base != syms:
            sys.exit("[ERROR] POTCAR 元素序 %s 与 POSCAR %s 不一致" %
                     (potcar_base, syms))
        (cfg / "INCAR").write_text(incar_text, encoding="utf-8", newline="\n")
        (cfg / "KPOINTS").write_text(kp_text, encoding="utf-8", newline="\n")
        (cfg / "POTCAR").write_text(pc_text, encoding="utf-8", newline="\n")
        stamp.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n",
                         encoding="utf-8")
        (cfg / "submit.sh").write_text(submit, encoding="utf-8", newline="\n")
        (cfg / "submit.sh").chmod(0o755)
        stepconf.apply_submit(cfg / "submit.sh", conf.submit)
        n_new += 1
    print("[DONE] %s：结构标签检查完成，保留已完成 %d 帧，新建/补全 %d 帧输入"
          % (step_name, n_kept, n_new))


if __name__ == "__main__":
    main()
