#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ke-dft-cpu 结果 -> 规定格式的绘图输出（图 + 同目录数据），对【已有结果】补生成，不必重跑 S8.3。

★ [0056-r5] 2026-10-18 用户定「统一」：本脚本不再有自己的一套算法，直接调用 S8.3 的
  step8.3_output/gen_step13_output.write_transport_plots() —— 与 S8.3 跑完时自动写的 <材料>/transport_plots/
  逐字同一份代码、同一套产物。以前本脚本自己算 BT2×DPT-τ（BT2 自己的 n(μ) 横轴、只取 xx 分量、
  τ(300 K)×300/T、κe 用 Wiedemann–Franz L0·σ·T、AMSET 也只取 xx），与 S8.3（BT2 插值到 AMSET 掺杂档、
  面内平均、每个温度的 DPT τ、BT2 自己的 κe/τ）两张图的蓝线对不上。

用法:
    python make_transport_plots.py <材料目录> [<输出目录>] [材料名]

<材料目录>  ke-dft-cpu 根（含 step8.3_output/、step8.4_amset2d/ 等；本地已 fetch 的 result/ 或集群目录）。
            需要 S8.3 已经写出 step8.3_output/intrinsic/intrinsic_transport_table.csv（S8.3 [0032] 起都有）。
<输出目录>  缺省 = <材料目录>/transport_plots（ke 目录内，与 S8.3 相同；目录名 transport_plots 时先清空旧结果）。
[材料名]    缺省 = <材料目录> 的上一级目录名。

产物（与 S8.3 相同，详见 gen_step13_output.write_transport_plots 的说明）:
    {transport_models,mobility_mechanisms}_carrier_tuning_1e16-1e21/<M>_<前缀>_carrier_tuning.png + _<T>K.csv
    {transport_models,mobility_mechanisms}_low_doping_<NLO>-1e15/<M>_<前缀>_low_doping.png + _<T>K.csv
    <M>_params_<T>K.csv、temperature_series/、mu_scan/、intrinsic_carrier/、defect_equilibrium/

输出位置约定: 只写 <输出目录>，绝不写进 ke-dft-cpu/step8.3_output/。
"""
import importlib.util
import os
import sys
from pathlib import Path

KE = Path(__file__).resolve().parents[1]


def _s83():
    for d in (KE.parents[1], KE, KE / "step8.4_amset2d", KE.parent / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("gen_step13_output", str(KE / "step8.3_output" / "gen_step13_output.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(mdir, out=None, M=None):
    mdir = Path(os.path.abspath(mdir))
    M = M or mdir.parent.name
    S = _s83()
    src = mdir / S.OUTDIR_NAME
    if not (src / "intrinsic" / "intrinsic_transport_table.csv").is_file():
        print("[skip] %s: 缺 %s/intrinsic/intrinsic_transport_table.csv（先跑 S8.3）" % (M, S.OUTDIR_NAME))
        return False
    out = Path(os.path.abspath(out)) if out else mdir / "transport_plots"
    S.write_transport_plots(mdir, src, S._read_dim(mdir) == "2d", plots=out, M=M)
    print("[%s] -> %s" % (M, out))
    return True


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ok = run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None, sys.argv[3] if len(sys.argv) > 3 else None)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
