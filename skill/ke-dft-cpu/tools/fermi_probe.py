#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fermi_probe.py —— AMSET 求不到费米能级的诊断（V131，不提交作业，几分钟）。

背景：WS2（KZ_CAP_2D 开、IR_FIX 关）在 set_doping_and_temperatures 报
    "Could not find fermi within 100.0% of concentration=-2.759e-07" / "Could not calculate Fermi level position"
同配置的 MoS₂、WSe2 没事；合成 DOS（带隙 1.5–2.5 eV、100–1000 K、同一浓度）上 AMSET 的贪心搜索也都能解，
所以问题在 WS2 的 DOS 本身。另一处隐患：AMSET 的 get_fermi 在中间量是 NaN 时不报错（NaN > tol 为假），
直接返回 NaN 当费米能级 —— 下游输运就会在 BoltzTraP2 的 pinv 报 "SVD did not converge"
（WS2 上一次 IR_FIX 开时正是这样崩的）。

本工具在 S8.4 运行目录按作业同一套插值重建能带与 DOS（加载同目录插件；IR_FIX 状态见 --ir-fix），
不算散射，然后：
  · DOS 诊断：本征费米能级、能带的 VBM/CBM、DOS 看到的带边、带隙里的态数（插值假态）、nelect、总态数；
  · 逐个 (掺杂, 温度)：
      AMSET 自己的 get_fermi（同样的容差阶梯 1e-5…1）—— 成/败、E_F、是否 NaN；
      AMSET 贪心搜索停在哪里、相对误差多大；
      二分法解 doping(E_F) = 目标（N(E_F) 单调，只要导带/价带有态就一定有解）—— E_F 与实现浓度。

用法（在 S8.4 运行目录）：
    python <skill>/ke-dft-cpu/tools/fermi_probe.py [--dir .] [--ir-fix auto|on|off] [--json fermi_probe.json]
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "step8.4_amset2d"))

HA_EV = 27.211386245988


def amset_ladder(dos, conc, T):
    """与 AmsetData.set_doping_and_temperatures 相同：tol = 1e-5 … 1 依次试 dos.get_fermi。
    返回 dict(ok, tol, ef_ha, nan)。"""
    for tol in np.logspace(-5, 0, 6):
        try:
            ef = float(dos.get_fermi(conc, T, tol=tol, precision=10))
        except ValueError:
            continue
        return {"ok": bool(np.isfinite(ef)), "tol": float(tol), "ef_ha": ef, "nan": not np.isfinite(ef)}
    return {"ok": False, "tol": None, "ef_ha": None, "nan": False}


def greedy_trace(dos, conc, T, nstep=50, step=0.1, precision=10):
    """照抄 AMSET FermiDos.get_fermi 的贪心网格搜索，返回它停下的 E_F（Ha）与最小相对误差。"""
    fermi = dos.efermi
    err = np.array([np.inf])
    for _ in range(precision):
        fr = np.arange(-nstep, nstep + 1) * step + fermi
        calc = np.array([dos.get_doping(f, T) for f in fr])
        err = np.abs(calc / conc - 1.0)
        fermi = fr[np.nanargmin(err)] if np.any(np.isfinite(err)) else fr[0]
        step /= 10.0
    return float(fermi), float(np.nanmin(err)) if np.any(np.isfinite(err)) else float("nan")


def bisect_fermi(dos, conc, T, pad_ev=1.0, xtol_ha=1e-12):
    """doping(E_F) − conc 的二分解。doping 随 E_F 单调减（N(E_F) 单调增）。无变号返回 None。"""
    from scipy.optimize import brentq
    lo = float(np.min(dos.energies)) - pad_ev / HA_EV
    hi = float(np.max(dos.energies)) + pad_ev / HA_EV
    g = lambda e: dos.get_doping(e, T) - conc  # noqa: E731
    glo, ghi = g(lo), g(hi)
    if not (np.isfinite(glo) and np.isfinite(ghi)) or glo * ghi > 0:
        return None
    return float(brentq(g, lo, hi, xtol=xtol_ha, maxiter=500))


def dos_diagnostics(dos, vbm_ha, cbm_ha, volume_bohr3):
    """DOS 的带边、带隙内的态数（每胞）、nelect、总态数。"""
    E, d, de = np.asarray(dos.energies), np.asarray(dos.tdos), float(dos.de)
    thr = 1e-8 * float(np.max(d)) if np.max(d) > 0 else 0.0
    ef = float(dos.efermi)
    below = (E <= ef) & (d > thr)
    above = (E > ef) & (d > thr)
    gap_mask = (E > vbm_ha + de) & (E < cbm_ha - de)
    return {"efermi_intrinsic_eV": ef * HA_EV, "efermi_finite": bool(np.isfinite(ef)),
            "vbm_bands_eV": vbm_ha * HA_EV, "cbm_bands_eV": cbm_ha * HA_EV,
            "gap_bands_eV": (cbm_ha - vbm_ha) * HA_EV,
            "vbm_dos_eV": float(E[below].max()) * HA_EV if below.any() else None,
            "cbm_dos_eV": float(E[above].min()) * HA_EV if above.any() else None,
            "states_in_gap_per_cell": float(d[gap_mask].sum() * de),
            "nelect": float(dos.nelect), "total_states": float(d.sum() * de),
            "volume_bohr3": float(volume_bohr3)}


def probe_dos(dos, dopings_cm3, temps, vbm_ha, cbm_ha, volume_bohr3):
    from amset.constants import cm_to_bohr
    rows = []
    for c_cm3 in dopings_cm3:
        conc = float(c_cm3) * (1 / cm_to_bohr) ** 3
        for T in temps:
            a = amset_ladder(dos, conc, float(T))
            gf, gerr = greedy_trace(dos, conc, float(T))
            b = bisect_fermi(dos, conc, float(T))
            ach = float(dos.get_doping(b, float(T))) if b is not None else None
            rows.append({"doping_cm3": float(c_cm3), "T": float(T), "conc_bohr3": conc,
                         "electrons_per_cell": conc * volume_bohr3,
                         "amset_ok": a["ok"], "amset_tol": a["tol"], "amset_nan": a["nan"],
                         "amset_ef_eV": a["ef_ha"] * HA_EV if a["ef_ha"] is not None and np.isfinite(a["ef_ha"]) else None,
                         "greedy_stop_eV": gf * HA_EV, "greedy_min_rel_err": gerr,
                         "bisect_ef_eV": b * HA_EV if b is not None else None,
                         "bisect_rel_err": abs(ach / conc - 1) if ach is not None else None})
    return {"dos": dos_diagnostics(dos, vbm_ha, cbm_ha, volume_bohr3), "rows": rows}


def band_edges(ad):
    vbm, cbm = -np.inf, np.inf
    for spin, e in ad.energies.items():
        vb = int(ad.vb_idx[spin]) if isinstance(ad.vb_idx, dict) else int(ad.vb_idx)
        vbm = max(vbm, float(np.max(e[: vb + 1])))
        cbm = min(cbm, float(np.min(e[vb + 1:])))
    return vbm, cbm


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=".")
    ap.add_argument("--ir-fix", default="auto", choices=("auto", "on", "off"),
                    help="重建时的不可约映射：auto = 看环境变量 AZ_IR_FIX（未设 = 关，与 V130 出厂默认一致）")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    import postprocess_intrinsic as P
    ir = {"on": True, "off": False}.get(a.ir_fix,
                                        os.environ.get("AZ_IR_FIX", "").strip().lower() in ("1", "true", "on", "yes"))
    ad, s, ver = P.build_amset_data(str(Path(a.dir).resolve()), ir_fix=ir, set_doping=False)
    vbm, cbm = band_edges(ad)
    rep = probe_dos(ad.dos, np.atleast_1d(s["doping"]), np.atleast_1d(s["temperatures"]), vbm, cbm,
                    ad.structure.volume)
    rep.update({"amset": ver, "ir_fix": ir, "dir": str(Path(a.dir).resolve())})
    d = rep["dos"]
    print("DOS：本征 E_F %.4f eV%s | 能带 VBM %.4f / CBM %.4f eV（带隙 %.3f）| DOS 带边 %s / %s eV"
          % (d["efermi_intrinsic_eV"], "" if d["efermi_finite"] else "（NaN！）", d["vbm_bands_eV"],
             d["cbm_bands_eV"], d["gap_bands_eV"], d["vbm_dos_eV"], d["cbm_dos_eV"]))
    print("     带隙内态数 %.3e /胞 | nelect %.6f | 总态数 %.3f" % (d["states_in_gap_per_cell"], d["nelect"],
                                                                  d["total_states"]))
    print("%12s %6s %11s %6s %12s %12s %11s %12s" % ("掺杂 cm^-3", "T", "载流子/胞", "AMSET", "AMSET E_F",
                                                      "贪心停在", "贪心误差", "二分 E_F"))
    for r in rep["rows"]:
        st = "NaN!" if r["amset_nan"] else ("ok" if r["amset_ok"] else "失败")
        f = lambda x: "—" if x is None else "%.4f" % x  # noqa: E731
        print("%12.3e %6.0f %11.3e %6s %12s %12.4f %11.3e %12s"
              % (r["doping_cm3"], r["T"], r["electrons_per_cell"], st, f(r["amset_ef_eV"]), r["greedy_stop_eV"],
                 r["greedy_min_rel_err"], f(r["bisect_ef_eV"])))
    bad = [r for r in rep["rows"] if not r["amset_ok"]]
    if not d["efermi_finite"]:
        print("判断：本征费米能级是 NaN —— AMSET 的 get_fermi 会静默返回 NaN，下游 pinv 报 SVD 不收敛。")
    elif bad and d["states_in_gap_per_cell"] > 0.1 * min(abs(r["electrons_per_cell"]) for r in bad):
        print("判断：带隙内有插值假态（%.2e /胞），与失败点的载流子数（≥ %.2e /胞）可比 —— 低掺杂点受它支配。"
              % (d["states_in_gap_per_cell"], min(abs(r["electrons_per_cell"]) for r in bad)))
    elif bad:
        print("判断：DOS 正常而 AMSET 的贪心搜索失败 —— 看'贪心停在'与'二分 E_F'差多少（搜索被局部极小困住）。")
    else:
        print("判断：所有 (掺杂, 温度) AMSET 都能解 —— 失败不在费米能级这一步，或与作业配置不同（--ir-fix？）。")
    if a.json:
        Path(a.json).write_text(json.dumps(rep, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
