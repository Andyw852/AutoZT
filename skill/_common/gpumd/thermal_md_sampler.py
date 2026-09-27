#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sample thermally displaced structures with a GPUMD NEP and write VASP POSCARs.

This is a proposal sampler only: every emitted structure still needs a converged DFT label.
Run from an isolated working directory so GPUMD append-only outputs cannot contaminate reruns.
"""
import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


def parse_csv_ints(s):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def stratified_validation_indices(counts, n_validation=10):
    """Allocate an exact validation budget proportionally, then space it in each stratum."""
    total = sum(counts)
    if total < n_validation or any(n < 1 for n in counts):
        raise ValueError("counts must contain at least n_validation candidate frames")
    quotas = [n_validation * n / total for n in counts]
    allocated = [int(q) for q in quotas]
    left = n_validation - sum(allocated)
    order = sorted(range(len(counts)), key=lambda i: (-(quotas[i] - allocated[i]), i))
    for i in order[:left]:
        allocated[i] += 1
    result = {}
    for i, (n, keep) in enumerate(zip(counts, allocated)):
        if keep == 0:
            result[i] = set()
        elif keep == 1:
            result[i] = {n // 2}
        else:
            result[i] = {round((k + 1) * (n - 1) / (keep + 1)) for k in range(keep)}
        if len(result[i]) != keep:
            raise ValueError("could not place unique validation indices")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poscar", required=True)
    ap.add_argument("--potential", required=True)
    ap.add_argument("--gpumd", required=True)
    ap.add_argument("--outdir", default="thermal_sample")
    ap.add_argument("--temperatures", default="300,500,650,800")
    ap.add_argument("--counts", default="13,12,12,13")
    ap.add_argument("--equil-steps", type=int, default=30000)
    ap.add_argument("--production-steps", type=int, default=100000)
    ap.add_argument("--sample-interval", type=int, default=5000)
    ap.add_argument("--time-step-fs", type=float, default=1.0)
    ap.add_argument("--tau-fs", type=float, default=100.0)
    args = ap.parse_args()

    temps = [float(x.strip()) for x in args.temperatures.split(",") if x.strip()]
    counts = parse_csv_ints(args.counts)
    if len(temps) != len(counts) or sum(counts) != 50:
        sys.exit("[ERROR] temperatures/counts 长度必须相同，且 counts 总和必须为 50")
    if any(n < 1 for n in counts) or args.production_steps // args.sample_interval < max(counts):
        sys.exit("[ERROR] production_steps/sample_interval 不足以提供所需帧数")

    from ase.io import read, write
    root = Path(args.outdir).resolve()
    if root.exists() and any(root.iterdir()):
        if (root / "thermal_manifest.json").is_file():
            sys.exit("[ERROR] 已存在完整 thermal_manifest.json；拒绝覆盖已生成的 50 帧")
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        archive = root.with_name(root.name + ".incomplete." + stamp)
        if archive.exists():
            sys.exit("[ERROR] 不完整输出归档目标已存在：%s" % archive)
        root.rename(archive)
        print("[WARN] 已保留不完整采样目录为 %s" % archive)
    root.mkdir(parents=True, exist_ok=True)
    poscar = Path(args.poscar).resolve()
    atoms = read(str(poscar), format="vasp")
    poscar_md5 = hashlib.md5(poscar.read_bytes()).hexdigest()
    model = root / "model.xyz"
    write(str(model), atoms, format="extxyz")
    potential = Path(args.potential).resolve()
    if not potential.is_file():
        sys.exit("[ERROR] 找不到采样 NEP 势：%s" % potential)

    candidate_records = []
    for T, need in zip(temps, counts):
        case = root / ("T%g" % T)
        case.mkdir()
        shutil.copy2(str(potential), str(case / "nep.txt"))
        shutil.copy2(str(model), str(case / "model.xyz"))
        run_in = "\n".join([
            "potential nep.txt", "velocity %g" % T,
            "time_step %g" % args.time_step_fs,
            "ensemble nvt_ber %g %g %g" % (T, T, args.tau_fs / args.time_step_fs),
            "dump_thermo 1000", "run %d" % args.equil_steps,
            "dump_xyz %d trajectory.xyz precision double" % args.sample_interval,
            "run %d" % args.production_steps, "",
        ])
        (case / "run.in").write_text(run_in, encoding="utf-8")
        log = case / "gpumd.log"
        with log.open("w", encoding="utf-8") as fh:
            rc = subprocess.call(["stdbuf", "-o0", "-e0", str(Path(args.gpumd).resolve())],
                                 cwd=str(case), stdout=fh, stderr=subprocess.STDOUT)
        if rc != 0 or not (case / "trajectory.xyz").is_file():
            sys.exit("[ERROR] T=%g K 采样失败，rc=%d；查看 %s" % (T, rc, log))
        frames = read(str(case / "trajectory.xyz"), index=":", format="extxyz")
        if len(frames) < need:
            sys.exit("[ERROR] T=%g K 只有 %d 帧，目标 %d" % (T, len(frames), need))
        if need == 1:
            selected = [frames[len(frames) // 2]]
        else:
            idx = [round(i * (len(frames) - 1) / (need - 1)) for i in range(need)]
            selected = [frames[i] for i in idx]
        for j, frame in enumerate(selected):
            candidate_records.append((T, j, frame))

    frame_dir = root / "frames"
    frame_dir.mkdir()
    manifest = {"schema": 1, "sampler": "GPUMD+NEP", "proposal_only": True,
                "temperatures_K": temps, "counts": counts, "total": 50,
                "train": 40, "validation": 10, "natoms": len(atoms),
                "reference_poscar_md5": poscar_md5, "frames": []}
    # Hold out exactly 10 frames, proportionally stratified by temperature.
    val_indices = stratified_validation_indices(counts, 10)
    serial = 0
    for candidate_idx, (T, j, frame) in enumerate(candidate_records):
        # candidate_records is ordered by temperature and local frame index.
        # Determine stratum index from the cumulative configured counts.
        stratum = 0
        offset = candidate_idx
        while offset >= counts[stratum]:
            offset -= counts[stratum]
            stratum += 1
        split = "validation" if j in val_indices[stratum] else "train"
        serial += 1
        name = "POSCAR-%03d" % serial
        write(str(frame_dir / name), frame, format="vasp", direct=True, vasp5=True, sort=False)
        manifest["frames"].append({"id": name, "file": str(Path("frames") / name),
                                   "temperature_K": T, "split": split,
                                   "source_frame_index": j})
    if sum(x["split"] == "validation" for x in manifest["frames"]) != 10:
        sys.exit("[ERROR] 50 帧的分层 train/validation 划分不是 40/10")
    (root / "thermal_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("[DONE] 50 帧：40 train + 10 validation；温度 %s K；%d atoms" %
          (", ".join(str(int(x)) for x in temps), len(atoms)))


if __name__ == "__main__":
    main()
