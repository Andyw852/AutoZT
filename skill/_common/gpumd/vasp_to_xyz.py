#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""vasp_to_xyz.py —— VASP OUTCAR 目录 → GPUMD/NEP 用的 extxyz（能量+力，可选 virial）。

用法（作业侧，在含 gpumd_params.json 的步骤目录里跑）：
  python vasp_to_xyz.py            # 读 gpumd_params.json
产物：train.xyz / test.xyz / structure_reference.xyz / xyz_summary.json
"""
import argparse
import glob
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc  # noqa: E402


def read_outcar(path):
    """读 OUTCAR 的最后一个离子步：能量/力/晶格/应力的原始张量。"""
    import numpy as np
    from ase.io import read
    atoms = read(path, index=-1, format="vasp-out")
    energy = atoms.get_potential_energy()
    forces = np.asarray(atoms.get_forces(), dtype=float)
    lat = np.asarray(atoms.cell.array, dtype=float)
    return atoms, lat, energy, forces


def parse_stress_kb(path):
    """OUTCAR 里最后一处 "in kB" 的 6 个应力分量（xx yy zz xy yz zx），失败返回 None。"""
    import numpy as np
    last = None
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if "in kB" in line:
                parts = line.split()
                try:
                    last = [float(x) for x in parts[-6:]]
                except ValueError:
                    pass
    return np.array(last, dtype=float) if last else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step-dir", default=".")
    ap.add_argument("--use-virial", type=int, default=0)
    args = ap.parse_args()
    p = gc.load_params(args.step_dir)
    if int(p.get("PRETRAINED_ONLY") or 0):
        # 只用预训练势：不读任何 OUTCAR，直接把材料 POSCAR 当 MD 参考结构。
        import numpy as np
        from ase.io import read
        ref_poscar = Path(args.step_dir, "reference_POSCAR")
        if not ref_poscar.is_file():
            sys.exit("[ERROR] PRETRAINED_ONLY=1 但缺 reference_POSCAR")
        atoms = read(str(ref_poscar))
        lat = np.asarray(atoms.cell.array, dtype=float)
        zeros = np.zeros((len(atoms), 3))
        Path(args.step_dir, "structure_reference.xyz").write_text(
            gc.extxyz_frame(lat, atoms.get_chemical_symbols(),
                            atoms.get_positions(), zeros, 0.0, None, "reference"))
        gc.write_summary(args.step_dir, "xyz_summary.json",
                         dict(XYZ_DONE=True, pretrained_only=True, n_frames=1,
                              n_train=0, n_test=0, natoms=len(atoms),
                              reference=str(ref_poscar)))
        print("[xyz] PRETRAINED_ONLY：以 %s 作 MD 参考结构（%d 原子）"
              % (ref_poscar, len(atoms)))
        return
    import numpy as np
    data_dir = p["DATA_DIR"]
    pattern = p.get("DATA_GLOB", "POSCAR-*/OUTCAR")
    files = sorted(glob.glob(os.path.join(data_dir, pattern)))
    if not files:
        sys.exit("[ERROR] 没找到 OUTCAR：%s/%s" % (data_dir, pattern))
    print("[xyz] 发现 %d 个 OUTCAR" % len(files))
    frames = []
    for i, f in enumerate(files):
        try:
            atoms, lat, energy, forces = read_outcar(f)
        except Exception as exc:
            print("[WARN] 跳过 %s：%s" % (f, exc))
            continue
        vol = float(abs(np.linalg.det(lat)))
        sigma = parse_stress_kb(f) if args.use_virial else None
        virial = None
        if sigma is not None:
            # OUTCAR stress in kB -> eV/Å^3 = kB/1602.1766208; virial = -sigma * V
            s = sigma / 1602.1766208   # xx yy zz xy yz zx
            voigt = np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]],
                              [s[5], s[4], s[2]]])
            virial = (-voigt * vol).reshape(-1).tolist()
        frames.append(dict(file=os.path.relpath(f, data_dir), atoms=atoms,
                           lat=lat, energy=float(energy), forces=forces,
                           virial=virial, volume=vol))
    if not frames:
        sys.exit("[ERROR] 没有任何可用的 OUTCAR 帧")
    energies = np.array([fr["energy"] for fr in frames])
    i_ref = int(np.argmin(energies))
    ref = frames[i_ref]
    rng = np.random.default_rng(int(p.get("SPLIT_SEED", 0)))
    order = rng.permutation(len(frames))
    n_train = int(round(len(frames) * float(p.get("TRAIN_FRAC", 0.9))))
    train_idx = set(order[:n_train].tolist())
    out = Path(args.step_dir)
    n_tr = n_te = 0
    force_sq = []
    with open(out / "train.xyz", "w") as ftr, open(out / "test.xyz", "w") as fte:
        ftr.write(gc.extxyz_frame(ref["lat"], ref["atoms"].get_chemical_symbols(),
                                  ref["atoms"].get_positions(), ref["forces"],
                                  ref["energy"], ref["virial"], "reference"))
        for k, fr in enumerate(frames):
            if k == i_ref:
                continue
            force_sq.append(float((fr["forces"] ** 2).mean()))   # 每分量的 F^2
            txt = gc.extxyz_frame(fr["lat"], fr["atoms"].get_chemical_symbols(),
                                  fr["atoms"].get_positions(), fr["forces"],
                                  fr["energy"], fr["virial"],
                                  "train" if k in train_idx else "test")
            if k in train_idx:
                ftr.write(txt); n_tr += 1
            else:
                fte.write(txt); n_te += 1
    (out / "structure_reference.xyz").write_text(
        gc.extxyz_frame(ref["lat"], ref["atoms"].get_chemical_symbols(),
                        ref["atoms"].get_positions(), ref["forces"], ref["energy"],
                        ref["virial"], "reference"))
    f_rms = float(np.sqrt(np.mean(np.array(force_sq)))) if force_sq else float("nan")
    natoms = len(ref["atoms"])
    summary = dict(XYZ_DONE=True, n_frames=len(frames), n_train=n_tr, n_test=n_te,
                   natoms=natoms, energy_min=float(energies.min()),
                   energy_max=float(energies.max()),
                   energy_span=float(energies.max() - energies.min()),
                   force_rms=float(f_rms),
                   reference=ref["file"], use_virial=bool(args.use_virial))
    gc.write_summary(args.step_dir, "xyz_summary.json", summary)
    if f_rms > 2.0:
        print("[WARN] 力 RMS = %.3f eV/Å 偏大，检查是否混入了未收敛/错误构型" % f_rms)


if __name__ == "__main__":
    main()
