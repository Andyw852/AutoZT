#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_deformation_h5.py —— 重新 gen S7.1 前后两份形变势 h5 的整场对比（V123 分诊用，不跑计算）。

为什么需要：MoS₂ 单层重新 gen S7.1 后带边 E1_vac iso 只从 8.52 变到 8.536，电子 ADP 迁移率却 198 -> 406
（V122）—— 只看带边 E1 判断不了要不要重跑 S8/S8.4，得看谷内整个 D(k) 场。

用法（在材料目录；先把旧文件备份再 retry S7.1）：
    cp step7b_deform_read/deformation_vac.h5 step7b_deform_read/deformation_vac.h5.pre-V121   # 2D
    cp step7b_deform_read/deformation.h5     step7b_deform_read/deformation.h5.pre-V121       # 3D
    （retry S7.1）
    python <skill>/ke-dft-cpu/tools/compare_deformation_h5.py \\
        step7b_deform_read/deformation_vac.h5.pre-V121 step7b_deform_read/deformation_vac.h5 \\
        [--band-edges step7b_deform_read/band_edges.json]

输出（每个自旋）：逐 (带, k) 的相对变化 |ΔD|_F / |D_old|_F 的中位 / p90 / 最大；
ADP 粗估因子 R = ⟨|D_new|²⟩ / ⟨|D_old|²⟩（全窗口等权；ADP 迁移率大致按 1/R 变）；给了 band_edges.json 时
另报带边 hit 点的 D_xx/D_yy 新旧值。
判定：|1 − R| > RERUN_R 或 p90 相对变化 > RERUN_P90 -> "需要重跑"；否则 "变化小，可不重跑"。
这是**粗估**（没有按费米窗口加权）；拿不准时以重跑为准。退出码：0 = 变化小，1 = 需要重跑，2 = 输入问题。
"""
import argparse
import json
import sys

import numpy as np

RERUN_R = 0.05        # ⟨|D|²⟩ 变 5% ≈ ADP 迁移率变 5%
RERUN_P90 = 0.10      # 九成 (带, k) 点的 D 变化都在 10% 以内才算"变化小"


def _load(path):
    import h5py
    with h5py.File(str(path), "r") as f:
        d = {k: np.array(f[k]) for k in f.keys() if k.startswith("deformation_potentials_")}
        kp = np.array(f["kpoints"])
        attrs = {k: (int(v) if np.ndim(v) == 0 else v) for k, v in f.attrs.items()}
    return d, kp, attrs


def compare(old_path, new_path, band_edges=None):
    d0, k0, a0 = _load(old_path)
    d1, k1, a1 = _load(new_path)
    if k0.shape != k1.shape or not np.allclose(np.mod(k0 + 1e-6, 1), np.mod(k1 + 1e-6, 1), atol=1e-4):
        raise ValueError("两份 h5 的 k 网格不同（%s vs %s）—— 不是同一批形变单点，无法逐点对比"
                         % (k0.shape, k1.shape))
    rep = {"old_attrs": a0, "new_attrs": a1, "spins": {}}
    need = False
    for key in sorted(set(d0) & set(d1)):
        D0, D1 = d0[key], d1[key]
        if D0.shape != D1.shape:
            raise ValueError("%s 形状不同 %s vs %s（带窗口变了？）" % (key, D0.shape, D1.shape))
        n0 = np.linalg.norm(D0, axis=(-2, -1))
        n1 = np.linalg.norm(D1, axis=(-2, -1))
        rel = np.linalg.norm(D1 - D0, axis=(-2, -1)) / np.maximum(n0, 1e-6 * max(n0.max(), 1e-12))
        R = float((n1 ** 2).mean() / max((n0 ** 2).mean(), 1e-30))
        s = {"R_mean_D2": round(R, 4), "mu_adp_factor_rough": round(1.0 / R, 3) if R else None,
             "rel_median": round(float(np.median(rel)), 4), "rel_p90": round(float(np.percentile(rel, 90)), 4),
             "rel_max": round(float(rel.max()), 4)}
        s["rerun"] = bool(abs(1 - R) > RERUN_R or s["rel_p90"] > RERUN_P90)
        need |= s["rerun"]
        rep["spins"][key] = s
    if band_edges:
        with open(band_edges, encoding="utf-8") as fh:
            be = json.load(fh)
        rows = []
        for car in ("electron", "hole"):
            for h in (be.get(car) or {}).get("hits") or []:
                ds = "deformation_potentials_%s" % h.get("spin")
                if ds in d0 and ds in d1:
                    b, k = int(h["band_h5"]), int(h["k_h5"])
                    rows.append({"carrier": car, "band": b, "k": k,
                                 "old_xx_yy": [round(float(d0[ds][b, k, 0, 0]), 4), round(float(d0[ds][b, k, 1, 1]), 4)],
                                 "new_xx_yy": [round(float(d1[ds][b, k, 0, 0]), 4), round(float(d1[ds][b, k, 1, 1]), 4)]})
        rep["band_edge_hits"] = rows
    rep["rerun"] = need
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--band-edges", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    try:
        rep = compare(a.old, a.new, a.band_edges)
    except (ValueError, KeyError, OSError) as e:
        print("[ERROR] %s" % e, file=sys.stderr)
        return 2
    for key, s in rep["spins"].items():
        print("%s：⟨|D|²⟩ 新/旧 = %.3f（ADP 迁移率粗估 ×%.2f）；逐点相对变化 中位 %.1f%% p90 %.1f%% 最大 %.1f%% -> %s"
              % (key, s["R_mean_D2"], s["mu_adp_factor_rough"] or float("nan"), 100 * s["rel_median"],
                 100 * s["rel_p90"], 100 * s["rel_max"], "需要重跑" if s["rerun"] else "变化小"))
    for r in rep.get("band_edge_hits", []):
        print("  带边 %-8s band=%d k=%d  D_xx/D_yy 旧 %s -> 新 %s" % (r["carrier"], r["band"], r["k"],
                                                                     r["old_xx_yy"], r["new_xx_yy"]))
    print("结论：%s" % ("需要重跑 S8/S8.4（形变势整场变化超阈值）" if rep["rerun"]
                     else "形变势整场变化小，旧的 S8/S8.4 结果可保留（粗估；拿不准时以重跑为准）"))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            fh.write(json.dumps(rep, ensure_ascii=False, indent=2, default=str))
    return 1 if rep["rerun"] else 0


if __name__ == "__main__":
    sys.exit(main())
