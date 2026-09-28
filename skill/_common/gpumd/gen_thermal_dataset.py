#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a fixed train/test split from old rattle data plus thermal DFT labels."""
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc
import stepconf

STEP = "step1_thermal_dataset"
OUTDIR = STEP
SPEC = {
    "GPUMD_BIN": ("/home/wangchaoyue852/gpumd/src/gpumd", "str"),
    "NEP_BIN": ("/home/wangchaoyue852/gpumd/src/nep", "str"),
    "SAMPLE_STEP": ("step1_thermal_sample", "str"),
    "LABEL_STEP": ("step1_thermal_label", "str"),
    "OLD_TRAIN": ("train_old.xyz", "str"),
    "OLD_TEST": ("test_old.xyz", "str"),
    "FOUNDATION_NEP": ("nep89_20250409.txt", "str"),
    "FOUNDATION_RESTART": ("nep89_20250409.restart", "str"),
    "EXPECTED_NEW_TRAIN": (40, "int"),
    "EXPECTED_NEW_TEST": (10, "int"),
    "NATOMS": (81, "int"),
}


def xyz_frame_count(path, expected_natoms):
    n = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        while True:
            first = fh.readline()
            if not first:
                break
            try:
                nat = int(first.strip())
            except ValueError:
                raise ValueError("%s: frame %d has invalid atom count" % (path, n + 1))
            header = fh.readline()
            if nat != expected_natoms or not header.startswith("Lattice="):
                raise ValueError("%s: frame %d header/natoms invalid" % (path, n + 1))
            for _ in range(nat):
                if not fh.readline():
                    raise ValueError("%s: frame %d is truncated" % (path, n + 1))
            n += 1
    if n == 0:
        raise ValueError("%s contains no XYZ frames" % path)
    return n


def _poscar_data(path):
    """Read the generated VASP5 Direct POSCAR with stdlib only."""
    lines = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    scale = float(lines[1].split()[0])
    cell = [[scale * float(v) for v in lines[i].split()[:3]] for i in (2, 3, 4)]
    symbols, counts = lines[5].split(), [int(x) for x in lines[6].split()]
    if len(symbols) != len(counts):
        raise ValueError("POSCAR species/count rows differ: %s" % path)
    atom_symbols = [el for el, n in zip(symbols, counts) for _ in range(n)]
    i = 7
    if lines[i].lower().startswith("selective"):
        i += 1
    mode = lines[i].lower()
    if not mode.startswith(("direct", "cartesian", "k")):
        raise ValueError("POSCAR has no Direct/Cartesian coordinate mode: %s" % path)
    positions = []
    for row in lines[i + 1:i + 1 + len(atom_symbols)]:
        xyz = [float(x) for x in row.split()[:3]]
        if mode.startswith("d"):
            positions.append([sum(xyz[k] * cell[k][j] for k in range(3))
                              for j in range(3)])
        else:
            positions.append([x * scale for x in xyz])
    if len(positions) != len(atom_symbols):
        raise ValueError("POSCAR coordinate count differs: %s" % path)
    return atom_symbols, positions, cell


def _read_static_outcar(path, natoms):
    """Read final VASP static energy(sigma->0), cell, positions, and forces.

    VASP OUTCAR uses fixed-width lattice columns; adjacent negative values can
    have no separating whitespace, so scan numeric tokens instead of split().
    """
    import re
    number = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")
    energy_pat = re.compile(r"energy\(sigma->0\)\s*=\s*([-+0-9.EeDd]+)")
    energy, cell, positions, forces = None, None, None, None
    lines = path.read_text(errors="ignore").splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if "energy(sigma->0)" in line:
            vals = energy_pat.findall(line)
            if vals:
                energy = float(vals[-1].replace("D", "E").replace("d", "e"))
        if "direct lattice vectors" in line.lower():
            rows = []
            for row in lines[i + 1:i + 4]:
                vals = number.findall(row)
                if len(vals) < 3:
                    rows = []
                    break
                rows.append([float(x.replace("D", "E").replace("d", "e"))
                             for x in vals[:3]])
            if len(rows) == 3:
                cell = rows
        if "POSITION" in line and "TOTAL-FORCE" in line:
            pos, frc = [], []
            j = i + 2  # header followed by a separator line
            for row in lines[j:j + natoms]:
                vals = number.findall(row)
                if len(vals) < 6:
                    pos, frc = [], []
                    break
                xyzf = [float(x.replace("D", "E").replace("d", "e"))
                        for x in vals[:6]]
                pos.append(xyzf[:3])
                frc.append(xyzf[3:6])
            if len(pos) == natoms:
                positions, forces = pos, frc
        i += 1
    if energy is None or cell is None or positions is None or forces is None:
        raise ValueError("OUTCAR misses final energy/cell/position/force block")
    return energy, cell, positions, forces

def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    sample_dir = cwd / str(conf["SAMPLE_STEP"]) / "sampled"
    label_dir = cwd / str(conf["LABEL_STEP"])
    manifest_path = sample_dir / "thermal_manifest.json"
    if not manifest_path.is_file():
        sys.exit("[ERROR] 缺热态 manifest：%s" % manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries = manifest.get("frames", [])
    tr = [x for x in entries if x.get("split") == "train"]
    te = [x for x in entries if x.get("split") == "validation"]
    if (len(tr), len(te)) != (int(conf["EXPECTED_NEW_TRAIN"]),
                              int(conf["EXPECTED_NEW_TEST"])):
        sys.exit("[ERROR] 热态 train/test=%d/%d，与设置不符" % (len(tr), len(te)))

    old_train, old_test = Path(str(conf["OLD_TRAIN"])), Path(str(conf["OLD_TEST"]))
    for f in (old_train, old_test):
        if not f.is_file():
            sys.exit("[ERROR] 缺原有训练数据 %s（由 AutoZT gen_need 推送）" % f)
    foundation_files = [Path(str(conf["FOUNDATION_NEP"])),
                        Path(str(conf["FOUNDATION_RESTART"]))]
    if not all(f.is_file() for f in foundation_files):
        sys.exit("[ERROR] NEP89 基座 txt/restart 未推送，不能保证一致微调")
    nat = int(conf["NATOMS"])
    old_ntr, old_nte = xyz_frame_count(old_train, nat), xyz_frame_count(old_test, nat)
    train_lines = old_train.read_text(encoding="utf-8").rstrip() + "\n"
    test_lines = old_test.read_text(encoding="utf-8").rstrip() + "\n"
    train_extra, test_extra = [], []
    force_sq, nforce = 0.0, 0
    label_input_hashes = None
    for i, ent in enumerate(entries, 1):
        cfg = label_dir / ("cfg-%03d" % i)
        hashes = {name: hashlib.sha256((cfg / name).read_bytes()).hexdigest()
                  for name in ("INCAR", "KPOINTS", "POTCAR")}
        if label_input_hashes is None:
            label_input_hashes = hashes
        elif hashes != label_input_hashes:
            sys.exit("[ERROR] DFT 标注输入在 cfg-%03d 与其它帧之间不一致" % i)
        outcar = cfg / "OUTCAR"
        if not outcar.is_file() or \
                "General timing and accounting informations" not in \
                outcar.read_text(errors="ignore")[-300000:]:
            sys.exit("[ERROR] DFT 帧未完整：%s" % outcar)
        outcar_text = outcar.read_text(errors="ignore")
        if "aborting loop because EDIFF is reached" not in outcar_text[-800000:]:
            sys.exit("[ERROR] %s 未找到 EDIFF 收敛标记；不把未收敛标签并入 NEP"
                     % outcar)
        # Ensure geometry and species are the exact frame that was labeled.
        symbols, ref_positions, ref_cell = _poscar_data(cfg / "POSCAR")
        energy, cell, positions, forces = _read_static_outcar(outcar, nat)
        if len(symbols) != nat or len(forces) != nat:
            sys.exit("[ERROR] %s 的原子数/力数组不符" % outcar)
        if max(abs(cell[i][j] - ref_cell[i][j]) for i in range(3) for j in range(3)) > 1e-6:
            sys.exit("[ERROR] %s 的 OUTCAR 晶胞与输入 POSCAR 不同" % outcar)
        if max(abs(positions[i][j] - ref_positions[i][j])
               for i in range(nat) for j in range(3)) > 1e-5:
            sys.exit("[ERROR] %s 的最终坐标偏离静态单点输入结构" % outcar)
        if len(positions) != nat or any(len(f) != 3 for f in forces):
            sys.exit("[ERROR] %s 元素顺序与输入 POSCAR 不同" % outcar)
        txt = gc.extxyz_frame(cell, symbols, positions, forces, energy, None,
                              "thermal_%gK" % float(ent["temperature_K"]))
        if ent["split"] == "train":
            train_extra.append(txt)
        elif ent["split"] == "validation":
            test_extra.append(txt)
        else:
            sys.exit("[ERROR] %s split 无效：%r" % (ent.get("id"), ent.get("split")))
        force_sq += sum(x * x for row in forces for x in row)
        nforce += sum(len(row) for row in forces)

    train_path, test_path = out / "train.xyz", out / "test.xyz"
    for f in foundation_files:
        # Ship the matching model/restart with this fetched dataset, so the next
        # AutoZT stage can relay the pair across jzzn -> 3090 without manual scp.
        import shutil
        shutil.copy2(f, out / f.name)
    train_path.write_text(train_lines + "".join(train_extra), encoding="utf-8", newline="\n")
    test_path.write_text(test_lines + "".join(test_extra), encoding="utf-8", newline="\n")
    ntr, nte = xyz_frame_count(train_path, nat), xyz_frame_count(test_path, nat)
    summary = {
        "DATASET_DONE": True,
        "natoms": nat,
        "old_train_frames": old_ntr,
        "old_test_frames": old_nte,
        "thermal_train_frames": len(train_extra),
        "thermal_validation_frames": len(test_extra),
        "train_frames": ntr,
        "test_frames": nte,
        "thermal_force_rms_eV_A": (force_sq / nforce) ** 0.5,
        "train_sha256": hashlib.sha256(train_path.read_bytes()).hexdigest(),
        "test_sha256": hashlib.sha256(test_path.read_bytes()).hexdigest(),
        "foundation_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                              for f in foundation_files},
        "dft_input_sha256": label_input_hashes,
        "temperatures_K": manifest["temperatures_K"],
        "split": "old 46/5 + thermal 40/10; validation held out by temperature",
    }
    gc.write_summary(out, "dataset_summary.json", summary)
    print("[DONE] 数据集：train=%d (旧%d+热态%d), test=%d (旧%d+热态%d)"
          % (ntr, old_ntr, len(train_extra), nte, old_nte, len(test_extra)))


if __name__ == "__main__":
    main()
