#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_transport_json.py —— 两份 AMSET transport.json 逐 (掺杂, 温度) 对比（V125 真实数据验证用，不跑计算）。

用法：
    python <skill>/ke-dft-cpu/tools/compare_transport_json.py 旧/transport.json 新/transport.json \\
        [--dims 2|3] [--tol 0.02] [--json 报告.json]

比较量：迁移率（overall + 每个散射机制）、电导率、Seebeck。2D（--dims 2，默认）取面内平均 (xx+yy)/2，
3D 取 tr/3 —— 平均量不受"六方面内 xx/yy 数值噪声"影响（V122：~2%，IR_FIX 开了以后这部分被对称化掉）。
另外报告新旧各自的面内各向异性 max|xx−yy|/avg，供核对。
判定：所有平均量的最大相对差 ≤ tol -> PASS（退出码 0）；否则 FAIL（1）；网格/键对不上 -> 2。
|旧值| 小于该量全表最大值 1e-6 的点不计相对差。
"""
import argparse
import json
import sys

import numpy as np


def _load(path):
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    out = {"doping": np.asarray(d["doping"], float), "temperatures": np.asarray(d["temperatures"], float)}
    for k in ("conductivity", "seebeck"):
        out[k] = np.asarray(d[k], float)
    for mech, arr in (d.get("mobility") or {}).items():
        out["mobility:%s" % mech] = np.asarray(arr, float)
    return out


def _avg(t, dims):
    return (t[..., 0, 0] + t[..., 1, 1]) / 2.0 if dims == 2 else np.trace(t, axis1=-2, axis2=-1) / 3.0


def _aniso(t):
    avg = (t[..., 0, 0] + t[..., 1, 1]) / 2.0
    big = np.abs(avg) > 1e-6 * max(np.abs(avg).max(), 1e-300)
    return float((np.abs(t[..., 0, 0] - t[..., 1, 1])[big] / np.abs(avg[big])).max()) if big.any() else 0.0


def compare(old_path, new_path, dims=2, tol=0.02):
    a, b = _load(old_path), _load(new_path)
    for k in ("doping", "temperatures"):
        if a[k].shape != b[k].shape or not np.allclose(a[k], b[k], rtol=1e-6):
            raise ValueError("%s 不同：%s vs %s —— 不是同一套设置，不能逐点比" % (k, a[k], b[k]))
    rep = {"dims": dims, "tol": tol, "quantities": {}, "only_old": [], "only_new": []}
    keys = [k for k in a if k not in ("doping", "temperatures")]
    rep["only_old"] = [k for k in keys if k not in b]
    rep["only_new"] = [k for k in b if k not in a and k not in ("doping", "temperatures")]
    ok = True
    for k in keys:
        if k not in b:
            continue
        if a[k].shape != b[k].shape:
            raise ValueError("%s 形状不同 %s vs %s" % (k, a[k].shape, b[k].shape))
        va, vb = _avg(a[k], dims), _avg(b[k], dims)
        big = np.abs(va) > 1e-6 * max(np.abs(va).max(), 1e-300)
        rel = np.zeros_like(va)
        rel[big] = np.abs(vb[big] - va[big]) / np.abs(va[big])
        i, j = np.unravel_index(int(np.argmax(rel)), rel.shape)
        q = {"max_rel": round(float(rel.max()), 5), "at_doping": float(a["doping"][i]),
             "at_T": float(a["temperatures"][j]), "old": float(va[i, j]), "new": float(vb[i, j]),
             "aniso_old": round(_aniso(a[k]), 5), "aniso_new": round(_aniso(b[k]), 5)}
        q["pass"] = bool(q["max_rel"] <= tol)
        ok &= q["pass"]
        rep["quantities"][k] = q
    rep["pass"] = bool(ok)
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--dims", type=int, choices=(2, 3), default=2)
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    try:
        rep = compare(a.old, a.new, a.dims, a.tol)
    except (ValueError, KeyError, OSError) as e:
        print("[ERROR] %s" % e, file=sys.stderr)
        return 2
    # V143：掺杂旁标出载流子（AMSET 约定：负 = 电子），别再靠符号去猜
    def _car(x):
        return {-1: "电子", 1: "空穴"}.get(int((x > 0) - (x < 0)), "本征")
    print("%-22s %9s %13s %5s %7s %11s %11s  %s" % ("量（%s）" % ("面内平均" if a.dims == 2 else "tr/3"), "最大相对差",
                                                 "掺杂", "载流子", "T", "各向异性旧", "各向异性新", "判定"))
    for k, q in rep["quantities"].items():
        print("%-22s %8.2f%% %13.4g %5s %7.0f %10.2f%% %10.2f%%  %s"
              % (k, 100 * q["max_rel"], q["at_doping"], _car(q["at_doping"]), q["at_T"], 100 * q["aniso_old"],
                 100 * q["aniso_new"], "PASS" if q["pass"] else "FAIL"))
    for k in rep["only_old"]:
        print("  只在旧结果里：%s" % k)
    for k in rep["only_new"]:
        print("  只在新结果里：%s" % k)
    print("结论：%s（阈值 %.1f%%）" % ("PASS" if rep["pass"] else "FAIL", 100 * a.tol))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            fh.write(json.dumps(rep, ensure_ascii=False, indent=2))
    return 0 if rep["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
