#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""nep_train.py —— 写 nep.in + 跑 nep（GPUMD v5.6）→ nep_summary.json。

微调（推荐）：fine_tune <foundation.txt> <foundation.restart>；架构键（version/zbl/
cutoff/n_max/basis_size/l_max/neuron）必须与基座一致、不可改（GPUMD 官方文档）。
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc  # noqa: E402


def foundation_types(path):
    """从基座 NEP .txt 首行读 (n_types, symbols)：格式 'nep4 89 H He Li ...'。
    微调时 nep.in 的 type 必须与基座一致（否则 nep 直接报 inconsistent）。"""
    with open(path) as fh:
        first = fh.readline().split()
    if len(first) < 3 or not first[1].isdigit():
        raise SystemExit("[ERROR] 无法从基座首行解析类型：%s" % " ".join(first[:4]))
    n = int(first[1])
    syms = first[2:2 + n]
    if len(syms) != n:
        raise SystemExit("[ERROR] 基座类型列表长度不符：%s" % path)
    return n, syms


def build_nep_in(p):
    lines = []
    if p.get("PRETRAINED_NEP"):
        lines.append("fine_tune %s %s" % (p["PRETRAINED_NEP"],
                                          p.get("PRETRAINED_RESTART", "")))
        _n, _syms = foundation_types(p["PRETRAINED_NEP"])
    else:
        _n, _syms = len(p["ELEMENTS"]), list(p["ELEMENTS"])
    lines.append("type %d %s" % (_n, " ".join(_syms)))
    # 训练侧的 type 是「元素种类数」（nep: "number of types should be integer"），
    # 元素符号由 train_xyz 里的 species 决定；预测侧才写符号列表。
    # 训练侧 type 需要 ≥2 个参数：类型数 + 元素符号列表；微调时**必须与基座模型一致**。
    lines.append("version %s" % p.get("NEP_VERSION", 4))
    if p.get("NEP_ZBL"):
        lines.append("zbl %s" % p["NEP_ZBL"])
    lines.append("cutoff %s" % " ".join(str(x) for x in p.get("CUTOFF", [6, 5])))
    lines.append("n_max %s" % " ".join(str(x) for x in p.get("N_MAX", [4, 4])))
    lines.append("basis_size %s" % " ".join(str(x) for x in p.get("BASIS_SIZE", [8, 8])))
    lines.append("l_max %s" % " ".join(str(x) for x in p.get("L_MAX", [4, 2, 1])))
    lines.append("neuron %s" % p.get("NEURON", 80))
    for key, dflt in (("LAMBDA_1", 0), ("LAMBDA_E", 1), ("LAMBDA_F", 1),
                      ("LAMBDA_V", 0), ("BATCH", 5000), ("POPULATION", 50),
                      ("GENERATION", 5000)):
        lines.append("%s %s" % (key.lower(), p.get(key, dflt)))
    # nep 不接受 train/test 文件参数：它在工作目录里固定读 train.xyz（test.xyz 可选）。
    lines.append("save_potential %s %s %s"
                 % (p.get("SAVE_EVERY", 1000), p.get("SAVE_START", 0),
                    p.get("SAVE_END", 1)))
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step-dir", default=".")
    ap.add_argument("--summary-only", action="store_true")
    args = ap.parse_args()
    p = gc.load_params(args.step_dir)
    nep_bin = p["NEP_BIN"]
    if not args.summary_only:
        text = build_nep_in(p)
        Path(args.step_dir, "nep.in").write_text(text, encoding="utf-8")
        print("[nep.in]\n%s" % text)
        if not os.path.isfile(Path(args.step_dir, "train.xyz")):
            sys.exit("[ERROR] 缺 train.xyz —— 先跑 step1_struct/S2 的数据转换")
        log = Path(args.step_dir, "train.log")
        with open(log, "w") as fh:
            rc = subprocess.call(["stdbuf", "-o0", "-e0", nep_bin], cwd=args.step_dir, stdout=fh,
                                 stderr=subprocess.STDOUT)
        print("[nep] rc=%d log=%s" % (rc, log))
    log_path = Path(args.step_dir, "train.log")
    text = log_path.read_text(errors="ignore") if log_path.is_file() else ""
    metrics = gc.parse_nep_log(text)
    loss_path = Path(args.step_dir, "loss.out")
    loss = gc.parse_loss_out(str(loss_path)) if loss_path.is_file() else {}
    fmax = float(p.get("FORCE_RMSE_MAX", 0.05))
    rel_max = float(p.get("FORCE_RMSE_REL_MAX", 0.05))
    # 力 RMSE 以 loss.out 最后一行为准（train.log 不含 RMSE）；优先用 test 力 RMSE 判合格。
    fr_train = loss.get("rmse_force_train")
    fr_test = loss.get("rmse_force_test")
    fr = fr_test if fr_test is not None else fr_train
    # DFT 力 RMS 来自 force_test.out（没有就用 force_train.out），用于相对误差门槛。
    rms_dft = None
    for fname in ("force_test.out", "force_train.out"):
        fp = Path(args.step_dir, fname)
        if fp.is_file():
            rms_dft = gc.force_file_rms(str(fp))
            if rms_dft:
                break
    fr_rel = (fr / rms_dft) if (fr is not None and rms_dft) else None
    # 合格判据：设了 FORCE_RMSE_REL_MAX 就按相对误差（力 RMSE / DFT 力 RMS），否则按绝对门槛。
    if rel_max > 0:
        gate_ok = fr_rel is None or fr_rel <= rel_max
    else:
        gate_ok = fr is None or fr <= fmax
    ok = Path(args.step_dir, "nep.txt").is_file() and gate_ok
    summary = dict(NEP_DONE=bool(ok), fine_tune=bool(p.get("PRETRAINED_NEP")),
                   generation=loss.get("generation", metrics.get("generation")),
                   rmse_energy=loss.get("rmse_energy_train", metrics.get("energy")),
                   rmse_force=fr_train, rmse_force_test=fr_test,
                   force_rms_dft_eV_A=rms_dft, force_rmse_rel=fr_rel,
                   force_rmse_max=fmax,
                   force_rmse_rel_max=(rel_max if rel_max > 0 else None),
                   model="nep.txt",
                   reason=None if ok else ("力 RMSE 不合格或没生成 nep.txt"))
    gc.write_summary(args.step_dir, "nep_summary.json", summary)
    if not ok:
        sys.exit("[ERROR] NEP 训练未通过：%s" % summary["reason"])


if __name__ == "__main__":
    main()
