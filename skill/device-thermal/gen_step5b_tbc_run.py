#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step5b_tbc_run.py -- (可选, GPU) 跑 GPUMD NEMD，并按"块数达标"判定完成.

为什么单独一步：S5_tbc 只生成 deck，S6_tbc_post 只做后处理。中间这段在 GPU 上跑
（3090 无真 SLURM，直接 bash run_gpumd.sh），可能几十分钟。放进 gen 步是让 autozt
的完成标记反映**真实跑完**，而不是"文件存在"——半截的 compute_chunk.out / 崩溃残留
都会被判成未完成。

完成判据（不要只看文件在不在）：
    compute_chunk.out 的块数 >= 测量段理论块数
    理论块数 = TBC_RUN_STEPS // (TBC_SAMPLE_INTERVAL * TBC_OUTPUT_INTERVAL)
      （GPUMD compute_chunk 每 sample_interval 步采一次、每 output_interval 个样本输出一次，
        所以每 si*oi 步出一块；2026-09-26 smoke 实测：20000 步 / (10*100) = 20 块）
    块数 = 剖面文件里第一列（bin 编号）为 0 的行数：每个 block 从 bin 0 开始。

读 step5_tbc/tbc_inputs.json（S5 产物）；产物 step5b_run/run_done.json
（marker: "TBC_RUN_DONE": true）。
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stepconf  # noqa: E402

STEP = "step5b_run"
OUTDIR = STEP
S5DIR = "step5_tbc"

SPEC = {
    "TBC_RUN_TIMEOUT_S": (0, "int"),   # 0 = 不限时
}


def count_blocks(path):
    """数 compute_chunk.out 的 block 数（第一列 == 0 的行数）。"""
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            try:
                cid = int(float(p[0]))
            except ValueError:
                continue
            if cid == 0:
                n += 1
    return n


def write_marker(**kw):
    os.makedirs(OUTDIR, exist_ok=True)
    p = os.path.join(OUTDIR, "run_done.json")
    d = {"TBC_RUN_DONE": True}
    d.update(kw)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(d, fh, ensure_ascii=False, indent=2)
    print("[OK] %s" % p)
    return p


def main():
    conf = stepconf.load(SPEC, STEP, strict=False)
    meta_p = os.path.join(S5DIR, "tbc_inputs.json")
    if not os.path.exists(meta_p):
        sys.exit("[ERROR] 缺 %s：请先跑 S5_tbc 生成 deck" % meta_p)
    with open(meta_p, encoding="utf-8") as fh:
        meta = json.load(fh)
    nrun = int(meta.get("run_steps") or 0)
    si = int(meta.get("sample_interval") or 0)
    oi = int(meta.get("output_interval") or 0)
    if nrun <= 0 or si <= 0 or oi <= 0:
        sys.exit("[ERROR] tbc_inputs.json 缺 run_steps/sample_interval/output_interval"
                 "（无法定理论块数）")
    expected = max(1, nrun // (si * oi))
    prof = os.path.join(S5DIR, "compute_chunk.out")
    nblk = count_blocks(prof) if os.path.exists(prof) else 0
    print("[S5b] 理论块数=%d（run_steps=%d / (sample=%d * output=%d)），现有=%d"
          % (expected, nrun, si, oi, nblk))

    if nblk < expected:
        runner = os.path.join(S5DIR, "run_gpumd.sh")
        if not os.path.exists(runner):
            sys.exit("[ERROR] 缺 %s：请先跑 S5_tbc 生成 run_gpumd.sh" % runner)
        print("[S5b] 块数不足，运行 GPUMD: bash %s" % runner)
        to = int(conf["TBC_RUN_TIMEOUT_S"] or 0)
        try:
            r = subprocess.run(["bash", runner], timeout=(to if to > 0 else None))
        except subprocess.TimeoutExpired:
            sys.exit("[ERROR] GPUMD 超时（TBC_RUN_TIMEOUT_S=%d s）；未完成，不写标记" % to)
        if r.returncode != 0:
            sys.exit("[ERROR] GPUMD 退出码 %d（见 %s/gpumd.log）；未完成，不写标记"
                     % (r.returncode, S5DIR))
        nblk = count_blocks(prof) if os.path.exists(prof) else 0
        print("[S5b] GPUMD 结束，现有块数=%d" % nblk)

    if nblk < expected:
        sys.exit("[ERROR] 只有 %d 块 < 理论 %d 块：NEMD 未跑完或输出被截断，"
                 "不写完成标记（可重跑本步或检查 %s/gpumd.log）" % (nblk, expected, S5DIR))

    write_marker(blocks=nblk, expected_blocks=expected, run_steps=nrun,
                 sample_interval=si, output_interval=oi, seed=meta.get("seed"),
                 structure=meta.get("structure"), note="块数达标")


if __name__ == "__main__":
    main()
