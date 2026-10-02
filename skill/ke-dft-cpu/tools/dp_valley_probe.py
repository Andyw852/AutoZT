#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dp_valley_probe.py —— 带边附近形变势 D(k) 的分布，按谷拆开（V141；只读、秒级、不提交作业）。

起因（MoSe2 单层 S8.4，300 K）：电子 ADP-only 迁移率是 DPT 的 28.5 倍（4640 vs 163），
MoS2 同一流程只有 1.9 倍；两者的空穴都约 0.7 倍。S8.4 用的是真空口径 deformation_vac.h5，S7 是
ISYM=0 + LVHAR，带边 E1 也过了 S7.1 的 5% 硬闸门 —— 带边那一点的 D 是对的。
DPT 只用带边一点的 E1；AMSET 的 ADP 用热分布范围内每个态自己的 D(k)。所以有两种可能：
  · 多谷：另一个低能谷（单层 TMD 的 Q/Λ 谷离 K 谷常只有几十 meV）占了相当多载流子、而那里的 D 小
    —— 这是物理，DPT 的单谷公式本来就看不到；
  · D 场本身不对：带边那一点对、它周围的 k 点却很小（形变单点之间带的次序/k 点对应出错等）。
本工具把这两种情况分开。

做法：读 S7.1 的 h5（默认 deformation_vac.h5）与 step7_deform/undeformed 的能带（AMSET 的
parse_calculation；h5 的 k 点比能带多时按 deform read 同样的方式 expand_bandstructure），band 轴按
get_ibands 映射到 h5，k 点按分数坐标对上（允许 −k）。对电子（CBM 以上 window 内）/空穴（VBM 以下）：
  · 权重 w = exp(−ΔE/kT)，D_iso = |D_xx + D_yy| / 2（与 band_edges.json 的 E1_iso 同一口径）；
  · 倍数 mu_factor = E1_带边² · Σ(w / D_iso²) / Σ w：各态（各谷）并联导电，μ_i ∝ τ_i ∝ 1/D_i²，
    所以平均的是 1/D²，不是 D²（只算 D 的分布，不管有效质量差别）≈ ADP-only 迁移率相对 DPT 的倍数。
    例：K 谷 D=7.4、Q 谷 D=1.0 且只高 30 meV（6 个 Q 对 2 个 K，Q 占约一半载流子）-> 约 27 倍；
  · 另报 rms_D = √(Σ w·D_iso² / Σ w) 供参考；
  · 按倒空间距离把态归成谷，逐谷列出：最低能量、权重占比、加权 rms D、谷底的 D_xx / D_yy。
判读（看"迁移率份"列：各谷对 ADP-only 迁移率的贡献占比）：
  · 贡献主要来自一个非带边谷（ΔE 几十 meV、D 小）-> 多谷物理：ADP-only 与单谷 DPT 本来就不可比，
    报告时注明；总迁移率里 ADP 只占一部分（POP/IMP 往往主导），影响要看 overall；
  · 带边谷自己的 rms D 就远小于带边 D -> D 场有问题，把输出发回来。
  · 倍数接近实测（MoSe2 电子 28.5）就说明 D 的分布已经解释了差异；差很远则还有别的原因（有效质量、
    谷间散射核等）。
用法（在材料目录，AMSET 0.5.1 环境）：
    python <skill>/ke-dft-cpu/tools/dp_valley_probe.py [--h5 step7b_deform_read/deformation_vac.h5]
                                                      [--T 300] [--window 0.3] [--json out.json]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

KB_EV = 8.617333262e-5
D_FLOOR = 0.05      # eV：D_iso 低于它按它算（1/D² 不发散），并单独计数


def match_kpoints(k_h5, k_bs, tol=1e-4):
    """h5 的每个 k -> k_bs 里同一点（或 −k，时间反演）的下标；找不到为 −1。"""
    k_bs = np.asarray(k_bs, float)
    out = np.full(len(k_h5), -1, int)
    for i, k in enumerate(np.asarray(k_h5, float)):
        for sgn in (1.0, -1.0):
            d = k_bs - sgn * k
            d -= np.rint(d)
            dd = np.linalg.norm(d, axis=1)
            j = int(np.argmin(dd))
            if dd[j] < tol:
                out[i] = j
                break
    return out


def _kdist(a, b, recip=None):
    d = np.asarray(a, float) - np.asarray(b, float)
    d -= np.rint(d)
    return float(np.linalg.norm(np.dot(d, recip) if recip is not None else d))


def valley_stats(energies, dtens, kfrac, edge, carrier, T=300.0, window=0.3, kclust=None, recip=None):
    """energies (nb, nk) eV 与 dtens (nb, nk, 3, 3) 在同一 band/k 轴上，kfrac (nk, 3)。
    kclust：归谷的倒空间距离（有 recip 时单位 Å⁻¹，含 2π；默认取面内最短倒格矢的 15%）。
    返回 dict 或 None（窗口内没有态）。"""
    E = np.asarray(energies, float)
    D = np.asarray(dtens, float)
    kfrac = np.asarray(kfrac, float)
    de = (E - edge) if carrier == "electron" else (edge - E)
    b, k = np.nonzero((de >= -1e-3) & (de <= window))
    if not len(b):
        return None
    dE = np.clip(de[b, k], 0.0, None)
    diso = np.abs(D[b, k, 0, 0] + D[b, k, 1, 1]) / 2.0
    w = np.exp(-dE / (KB_EV * T))
    rms = float(np.sqrt((w * diso ** 2).sum() / w.sum()))
    inv2 = 1.0 / np.maximum(diso, D_FLOOR) ** 2
    if kclust is None:
        if recip is not None:
            kclust = 0.15 * float(min(np.linalg.norm(recip[0]), np.linalg.norm(recip[1])))
        else:
            kclust = 0.15
    valleys = []
    for i in np.argsort(dE, kind="stable"):
        kf = kfrac[k[i]]
        for v in valleys:
            if _kdist(kf, v["k"], recip) < kclust:
                v["idx"].append(int(i))
                break
        else:
            valleys.append({"k": kf, "idx": [int(i)]})
    rows = []
    for v in valleys:
        ii = np.asarray(v["idx"])
        wv = w[ii]
        j = ii[0]
        rows.append({"k_frac": np.round(v["k"], 4).tolist(), "dE_meV": round(1000.0 * float(dE[j]), 1),
                     "weight": round(float(wv.sum() / w.sum()), 4),
                     "mu_share": round(float((wv * inv2[ii]).sum() / (w * inv2).sum()), 4),
                     "rms_D": round(float(np.sqrt((wv * diso[ii] ** 2).sum() / wv.sum())), 4),
                     "D_xx": round(float(D[b[j], k[j], 0, 0]), 4), "D_yy": round(float(D[b[j], k[j], 1, 1]), 4),
                     "n_states": int(len(ii))})
    rows.sort(key=lambda r: -r["weight"])
    edge_D = float(diso[int(np.argmin(dE))])
    mu_factor = max(edge_D, D_FLOOR) ** 2 * float((w * inv2).sum() / w.sum())
    return {"carrier": carrier, "T": T, "window_eV": window, "n_states": int(len(b)),
            "edge_D": round(edge_D, 4), "rms_D": round(rms, 4), "mu_factor": round(mu_factor, 3),
            "n_below_floor": int((diso < D_FLOOR).sum()), "valleys": rows}


def _report(r, top=8):
    print("== %s（%g K，带边 %s %.2f eV 以内，%d 个态）"
          % (r["carrier"], r["T"], "以上" if r["carrier"] == "electron" else "以下", r["window_eV"], r["n_states"]))
    print("   带边 D_iso = %.3f eV；倍数 E1²·⟨1/D²⟩ = %.2f（≈ 只因 D 的分布，ADP-only 迁移率相对 DPT 的倍数）；"
          "加权 rms D = %.3f eV%s"
          % (r["edge_D"], r["mu_factor"], r["rms_D"],
             "；%d 个态 D < %.2f eV" % (r["n_below_floor"], D_FLOOR) if r["n_below_floor"] else ""))
    print("   %-24s %9s %8s %8s %8s %8s %8s %5s"
          % ("谷（谷底 k_frac）", "ΔE/meV", "载流子", "迁移率份", "rmsD", "D_xx", "D_yy", "态数"))
    for v in r["valleys"][:top]:
        print("   %-24s %9.1f %8.3f %8.3f %8.3f %8.3f %8.3f %5d"
              % (str(v["k_frac"]), v["dE_meV"], v["weight"], v["mu_share"], v["rms_D"], v["D_xx"], v["D_yy"],
                 v["n_states"]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("path", nargs="?", default=".", help="材料目录（默认当前目录）")
    ap.add_argument("--h5", default=None, help="默认 step7b_deform_read/deformation_vac.h5，没有则 deformation.h5")
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--window", type=float, default=0.3, help="带边以内多少 eV 的态（默认 0.3）")
    ap.add_argument("--kclust", type=float, default=None, help="归谷距离（Å⁻¹，含 2π）；默认面内最短倒格矢的 15%%")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    mat = Path(a.path)
    rd = mat / "step7b_deform_read"
    h5 = Path(a.h5) if a.h5 else (rd / "deformation_vac.h5" if (rd / "deformation_vac.h5").is_file()
                                  else rd / "deformation.h5")
    und = mat / "step7_deform" / "undeformed"
    if not h5.is_file() or not und.is_dir():
        print("[ERROR] 找不到 %s 或 %s" % (h5, und), file=sys.stderr)
        return 2
    from amset.constants import defaults
    from amset.deformation.io import load_deformation_potentials, parse_calculation
    from amset.electronic_structure.common import get_ibands
    dp, kpoints, structure = load_deformation_potentials(str(h5))
    bs = parse_calculation(str(und))["bandstructure"]
    kbs = np.array([kk.frac_coords for kk in bs.kpoints])
    kmap = match_kpoints(kpoints, kbs)
    if (kmap < 0).any():
        from amset.electronic_structure.symmetry import expand_bandstructure
        bs = expand_bandstructure(bs, symprec=defaults["symprec"])
        kbs = np.array([kk.frac_coords for kk in bs.kpoints])
        kmap = match_kpoints(kpoints, kbs)
    if (kmap < 0).any():
        print("[ERROR] h5 的 %d 个 k 点在 undeformed 能带里找不到（不是同一批计算？）"
              % int((kmap < 0).sum()), file=sys.stderr)
        return 2
    ib = get_ibands(defaults["energy_cutoff"], bs)
    recip = np.asarray(structure.lattice.reciprocal_lattice.matrix)
    cbm, vbm = bs.get_cbm()["energy"], bs.get_vbm()["energy"]
    if cbm is None or vbm is None:
        print("[ERROR] 能带没有带边（金属？）", file=sys.stderr)
        return 2
    print("h5：%s（%d 个 k 点）；能带：%s；CBM %.4f eV，VBM %.4f eV" % (h5, len(kpoints), und, cbm, vbm))
    be = rd / "band_edges.json"
    if be.is_file():
        bej = json.loads(be.read_text(encoding="utf-8"))
        print("band_edges.json：电子 E1_iso %s eV，空穴 E1_iso %s eV（真空口径原值见文件）"
              % ((bej.get("electron") or {}).get("E1_iso_eV"), (bej.get("hole") or {}).get("E1_iso_eV")))
    out = []
    for spin, D in dp.items():
        idx = np.asarray(ib[spin])
        if D.shape[0] != len(idx):
            print("[ERROR] band 轴与 h5 不一致（h5=%d, ibands=%d）" % (D.shape[0], len(idx)), file=sys.stderr)
            return 2
        E = np.asarray(bs.bands[spin])[idx][:, kmap]
        for carrier, edge in (("electron", cbm), ("hole", vbm)):
            r = valley_stats(E, D, kpoints, edge, carrier, a.T, a.window, a.kclust, recip)
            if r is None:
                continue
            r["spin"] = spin.name
            out.append(r)
            _report(r)
    if a.json:
        Path(a.json).write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
