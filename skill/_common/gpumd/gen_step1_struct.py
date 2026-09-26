#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_struct.py —— 结构/数据检查（step1_struct，run: gen，不提交作业）。

检查材料 POSCAR 与训练数据目录，落盘 gpumd_params.json，写 struct_summary.json。
真正的 OUTCAR→extxyz 转换在 S2 的作业里做（数据可能要先从别的集群拷过来）。
"""
import glob
import math
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
    # 1 = 必须能找到训练数据（微调 / 从头训练）；0 = 允许没有（PRETRAINED_ONLY=1 路线）。
    "REQUIRE_TRAINING_DATA": (1, "int"),
    "TRAIN_FRAC": (0.9, "float"),
    "SPLIT_SEED": (20260925, "int"),
    "ELEMENTS": ("Mg C", "str"),
    "STRUCT_SRC": ("dft_ref", "str"),
    "IFC_KAPPA_RTA_DIR": ("", "str"),
    "KAPPA_SG_AVG": (0.0, "float"),
    # ★ 2026-09-26 user 要求：训练数据 k 网格体检（S1_struct 是唯一"数据入口"，
    #   此前对 DATA_DIR 里的 k 网格零校验）。
    "KDELTA_CHECK": ("on", "str"),     # on/off：是否抽样检查训练帧的 KPOINTS
    "KDELTA_TOL": (0.20, "float"),     # Δk 相对容差（帧间不一致判据）
    "KDELTA_BLOCK": (0, "int"),        # 1 = 体检不过 → sys.exit（默认 0 只 WARN）
}


def _lattice_of(frame_dir):
    """取帧目录的晶格矩阵（3x3，Å）。优先 POSCAR/CONTCAR，回退 OUTCAR 的
    'direct lattice vectors' 头三行。取不到返回 None。"""
    for name in ("POSCAR", "CONTCAR"):
        p = frame_dir / name
        if not p.is_file():
            continue
        try:
            ln = [x for x in p.read_text(errors="ignore").splitlines() if x.strip()]
            sc = float(ln[1].split()[0])
            return [[float(x) * sc for x in ln[2 + i].split()[:3]] for i in range(3)]
        except Exception:                                # noqa: BLE001
            continue
    p = frame_dir / "OUTCAR"
    if not p.is_file():
        return None
    try:
        txt = p.read_text(errors="ignore")
    except OSError:
        return None
    idx = txt.find("direct lattice vectors")
    if idx < 0:
        return None
    out = []
    for r in txt[idx:].splitlines()[1:4]:
        t = r.split()
        try:
            out.append([float(t[0]), float(t[1]), float(t[2])])
        except (ValueError, IndexError):
            return None
    return out if len(out) == 3 else None


def _mesh_of(frame_dir):
    """取帧目录 KPOINTS 的自动网格行（第 4 行）。没有 / 非自动网格 → None。"""
    p = frame_dir / "KPOINTS"
    if not p.is_file():
        return None
    try:
        ln = p.read_text(errors="ignore").splitlines()
        if len(ln) >= 4 and ln[1].strip().startswith("0"):
            return " ".join(ln[3].split()[:3])
    except OSError:
        return None
    return None


def _kdelta(mesh, lat):
    """胞不变 k 点间距 Δk = 2π·max_i(|b_i|/N_i)（Å⁻¹）；越大 = k 采样越粗。

    ★ 与 mlff/dft_settings.kspacing_density（= 2π/min_i(N_i·|b_i|)，∝ 胞长/N，
    **不是胞不变量**）不同：只有 Δk 能在"原胞 ↔ m 倍超胞"之间直接比较。
    """
    if not mesh or lat is None:
        return 0.0
    try:
        nums = [int(x) for x in mesh.split()[:3]]
    except ValueError:
        return 0.0
    a = lat
    det = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
           - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    if abs(det) < 1e-12:
        return 0.0
    c0 = [a[1][1] * a[2][2] - a[1][2] * a[2][1],
          a[1][2] * a[2][0] - a[1][0] * a[2][2],
          a[1][0] * a[2][1] - a[1][1] * a[2][0]]
    c1 = [a[2][1] * a[0][2] - a[2][2] * a[0][1],
          a[2][2] * a[0][0] - a[2][0] * a[0][2],
          a[2][0] * a[0][1] - a[2][1] * a[0][0]]
    c2 = [a[0][1] * a[1][2] - a[0][2] * a[1][1],
          a[0][2] * a[1][0] - a[0][0] * a[1][2],
          a[0][0] * a[1][1] - a[0][1] * a[1][0]]
    inv = [[c0[0] / det, c1[0] / det, c2[0] / det],
           [c0[1] / det, c1[1] / det, c2[1] / det],
           [c0[2] / det, c1[2] / det, c2[2] / det]]
    rec = [(sum(inv[i][j] ** 2 for j in range(3))) ** 0.5 for i in range(3)]
    return 2.0 * math.pi * max(rec[i] / max(nums[i], 1) for i in range(3))


def check_training_data_kpoints(data_dir, glob_pat, tol=0.20, block=False):
    """★ 训练数据 k 网格体检（2026-09-26，user 要求）。返回 (summary, bad)。

    为什么必须查：S1_struct 此前只要有 OUTCAR 就放行，对 DATA_DIR 里的 k 网格
    **零校验**。而 2026-08 那批 51 帧 / 243 原子数据集正是 Γ-only（超胞 Γ-only
    ≡ 原胞 3×3×3）：实测参考构型 max|F| = 0.103 eV/Å，而 6×6×6 只要 1.0e-4，
    差 1000 倍 —— 这个误差会原样进 NEP/MACE 的训练标签。
    block=False（默认）只 WARN 放行；block=True 时未通过 → sys.exit。
    """
    frames = sorted(glob.glob(os.path.join(data_dir, glob_pat)))
    dirs = sorted({os.path.dirname(f) for f in frames})
    seen, dks, n_nok = {}, [], 0
    for d in dirs:
        fd = Path(d)
        mesh = _mesh_of(fd)
        if mesh is None:
            n_nok += 1
            continue
        seen[mesh] = seen.get(mesh, 0) + 1
        dk = _kdelta(mesh, _lattice_of(fd))
        if dk > 0:
            dks.append(dk)
    summ = {"n_frames": len(frames), "n_dirs": len(dirs), "meshes": seen,
            "n_no_kpoints": n_nok,
            "dk_1A_min": (min(dks) if dks else 0.0),
            "dk_1A_max": (max(dks) if dks else 0.0)}
    if not seen:
        print("[WARN] 训练数据体检：%d 个帧目录里一个 KPOINTS 都没读到，无法判定 k 网格"
              % len(dirs))
        return summ, False
    print("[..] 训练数据体检：%d 帧 / %d 目录；KPOINTS 网格 %s；Δk %.5f–%.5f 1/Å"
          % (len(frames), len(dirs),
             ", ".join("%s x%d" % (k, v) for k, v in
                       sorted(seen.items(), key=lambda kv: -kv[1])),
             summ["dk_1A_min"], summ["dk_1A_max"]))
    bad = False
    if len(seen) > 1:
        bad = True
        print("[WARN] 训练数据的 k 网格不唯一（%d 种）—— 各帧的 k 点系统误差不同，"
              "会原样进训练标签" % len(seen))
    if dks and min(dks) > 0 and (max(dks) / min(dks) - 1.0) > tol:
        bad = True
        print("[WARN] 训练数据的 k 点间距 Δk 跨度 %.0f%% > 容差 %.0f%%"
              % ((max(dks) / min(dks) - 1.0) * 100, tol * 100))
    if list(seen) == ["1 1 1"]:
        bad = True
        print("[WARN] ★★ 训练数据全是 Γ-only —— 对超胞取力这等价于只在原胞 3×3×3 "
              "采样。实测会把参考构型 max|F| 从 1e-4 抬到 0.1 eV/Å 量级。"
              "强烈建议改用与生产原胞网格配套的超胞网格重做数据。")
    if bad and block:
        sys.exit("[ERROR] 训练数据 k 网格体检未通过（KDELTA_BLOCK=1）。\n"
                 "        详见上面 WARN 行；确认要放行就把 KDELTA_BLOCK 设 0。")
    if not bad:
        print("[OK] 训练数据 k 网格体检通过")
    return summ, bad


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
            if int(conf.get("REQUIRE_TRAINING_DATA", 1) or 0):
                sys.exit("[ERROR] %s 下没找到 %s（数据还没从 jzzn 拷过来？）"
                         % (data_dir, conf["DATA_GLOB"]))
            print("[WARN] %s 下没找到 %s；REQUIRE_TRAINING_DATA=0 继续"
                  "（PRETRAINED_ONLY 路线不需要训练数据）"
                  % (data_dir, conf["DATA_GLOB"]))
    # ★ 训练数据 k 网格体检（2026-09-26，user 要求）：S1 是唯一数据入口，先把住。
    ksumm = {}
    if data_dir and n_outcar and str(conf.get("KDELTA_CHECK", "on")).strip().lower() \
            not in ("off", "0", "false", "no"):
        ksumm, _kbad = check_training_data_kpoints(
            data_dir, conf["DATA_GLOB"], tol=float(conf["KDELTA_TOL"]),
            block=int(conf.get("KDELTA_BLOCK", 0) or 0) == 1)
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
                          gpumd_bin=conf["GPUMD_BIN"], nep_bin=conf["NEP_BIN"],
                          kpoints_check=ksumm))
    print("[DONE] %s：参数与结构就绪（%d 个 OUTCAR）" % (OUTDIR, n_outcar))


if __name__ == "__main__":
    main()
