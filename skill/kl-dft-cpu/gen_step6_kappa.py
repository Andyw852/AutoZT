#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step6_kappa.py —— kl-dft-cpu S6_kappa：BTE 晶格热导率（共享 fit-fc-thermal 引擎）。

取 step5_fc 的 fc2/fc3 → 三声子 BTE（RTA/LBTE，自动 q 网格收敛 + 逐档截断三判据选择），
写 kappa_summary.json（marker: KAPPA_DONE）。计算节点跑 kappa_driver.py。

注意：换共享引擎后，旧引擎的 SOLVER=shengbte/fourphonon 尚未并入共享 kappa_driver，
后续单独补（与 kl-mlff 的路线 C 一致）。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# 共享引擎（_common/fcfit）+ stepconf/dim_common（_common/opt）+ _common
_base = Path(__file__).resolve().parent.parent
for _d in (_base / "fcfit", _base / "opt", _base):
    if _d.is_dir() and str(_d) not in sys.path:
        sys.path.append(str(_d))

import gen_step2_kappa as g2  # noqa: E402

OUTDIR = "step6_kappa"
STEP = "step6_kappa"
FIT_DIR = "step5_fc"


METHOD_NAME = "kl_method.json"   # 本链几何用的泛函，随 step6_kappa 拉回，S20_zt 拿它和 ke 的 FUNC 对照


def write_kl_method(cwd=None):
    """2026-10-14：把 step1 的 FUNC / FUNC_REASON / DIM 写进 step6_kappa/kl_method.json。

    zt-dft-cpu 只拿到 ke 的 workflow_method.txt 与 kl 的 kappa_summary.json，后者不带泛函，
    κ_e（ke）与 κ_L（kl）是不是同一个泛函的晶格以前无从核对。本文件随 S6 结果目录拉回、推到
    S20 的 inputs/kappa_L/，S20 据此告警。读不到 step1 就不写（S20 按"未知"处理）。"""
    import json
    cwd = Path(cwd or Path.cwd())
    for f in sorted(cwd.glob("step1*/workflow_method.txt")):
        kv = {}
        for ln in f.read_text(encoding="utf-8", errors="ignore").splitlines():
            if "=" in ln:
                k, v = ln.split("=", 1)
                kv[k.strip().upper()] = v.strip()
        if kv.get("FUNC"):
            out = cwd / OUTDIR
            out.mkdir(exist_ok=True)
            rec = {"skill": "kl-dft-cpu", "func": kv["FUNC"].lower(),
                   "gga": kv.get("GGA"), "ivdw": kv.get("IVDW"),
                   "func_reason": kv.get("FUNC_REASON"), "dim": kv.get("DIM"),
                   "source": "%s/%s" % (f.parent.name, f.name)}
            (out / METHOD_NAME).write_text(json.dumps(rec, ensure_ascii=False, indent=1),
                                           encoding="utf-8")
            print("[..] %s/%s：FUNC=%s（S20_zt 据此核对 ke/kl 泛函一致）"
                  % (OUTDIR, METHOD_NAME, rec["func"]))
            return rec
    print("[WARN] 读不到 step1*/workflow_method.txt 的 FUNC —— 不写 %s，S20_zt 无法核对 ke/kl 泛函"
          % METHOD_NAME)
    return None


def main():
    g2.main(outdir=OUTDIR, fit_dir=FIT_DIR, step=STEP)
    write_kl_method()


if __name__ == "__main__":
    main()
