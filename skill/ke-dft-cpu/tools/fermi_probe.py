#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fermi_probe.py —— AMSET 求不到费米能级的诊断（V131，不提交作业，几分钟）。

背景：WS2（KZ_CAP_2D 开、IR_FIX 关）在 set_doping_and_temperatures 报
    "Could not find fermi within 100.0% of concentration=-2.759e-07" / "Could not calculate Fermi level position"
同配置的 MoS₂、WSe2 没事；合成 DOS（带隙 1.5–2.5 eV、100–1000 K、同一浓度）上 AMSET 的贪心搜索也都能解，
所以问题在 WS2 的 DOS 本身。另一处隐患：AMSET 的 get_fermi 在中间量是 NaN 时不报错（NaN > tol 为假），
直接返回 NaN 当费米能级 —— 下游输运就会在 BoltzTraP2 的 pinv 报 "SVD did not converge"
（WS2 上一次 IR_FIX 开时正是这样崩的）。

V132 结论（WS2 实测）：DOS 正常，是 AMSET 的贪心搜索被双精度舍入困住（低温时第二层 0.272 eV 网格上
离目标最近的点载流子数低于 ~nelect·1e-16 /胞，整个带隙误差都"等于 1"）。修法见 step8.4_amset2d/amset_fermi_fix.py。
"AMSET"一列始终是**原搜索**（不经 amset_fermi_fix），"实际误差"是原搜索被接受时的掺杂偏差 —— 用来核对旧作业
有没有在宽松容差下接受了偏差大的解（> 1% 的点旧结果不可信）。3D 的 S8 运行目录同样可用。

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


def _original_get_fermi(dos):
    """AMSET 原搜索（剥掉 amset_fermi_fix 与旧 amset2d_plugin 的外壳）。"""
    import amset_fermi_fix
    return amset_fermi_fix.unwrap(type(dos).get_fermi)


def amset_ladder(dos, conc, T):
    """与 AmsetData.set_doping_and_temperatures 相同：tol = 1e-5 … 1 依次试 AMSET 原 get_fermi。
    返回 dict(ok, tol, ef_ha, nan, rel_err)；rel_err = 被接受的解的实际掺杂偏差。"""
    fn = _original_get_fermi(dos)
    for tol in np.logspace(-5, 0, 6):
        try:
            with np.errstate(all="ignore"):
                ef = float(fn(dos, conc, T, tol=tol, precision=10))
        except ValueError:
            continue
        ok = bool(np.isfinite(ef))
        with np.errstate(all="ignore"):
            err = float(abs(dos.get_doping(ef, T) / conc - 1.0)) if ok else None
        return {"ok": ok, "tol": float(tol), "ef_ha": ef, "nan": not ok, "rel_err": err}
    return {"ok": False, "tol": None, "ef_ha": None, "nan": False, "rel_err": None}


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


def bisect_fermi(dos, conc, T):
    """doping(E_F) − conc 的二分解（Ha）—— 与 amset_fermi_fix 作业里用的是同一个函数。无变号返回 None。"""
    import amset_fermi_fix
    return amset_fermi_fix.solve_fermi(dos, conc, T)


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
                         "amset_rel_err": a["rel_err"],
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
    print("%12s %6s %11s %6s %12s %10s %12s %11s %12s" % ("掺杂 cm^-3", "T", "载流子/胞", "AMSET", "AMSET E_F",
                                                           "实际误差", "贪心停在", "贪心误差", "二分 E_F"))
    for r in rep["rows"]:
        st = "NaN!" if r["amset_nan"] else ("失败" if not r["amset_ok"] else
                                             ("宽松!" if r["amset_rel_err"] > 0.01 else "ok"))
        f = lambda x: "—" if x is None else "%.4f" % x  # noqa: E731
        print("%12.3e %6.0f %11.3e %6s %12s %10s %12.4f %11.3e %12s"
              % (r["doping_cm3"], r["T"], r["electrons_per_cell"], st, f(r["amset_ef_eV"]),
                 "—" if r["amset_rel_err"] is None else "%.1e" % r["amset_rel_err"], r["greedy_stop_eV"],
                 r["greedy_min_rel_err"], f(r["bisect_ef_eV"])))
    bad = [r for r in rep["rows"] if not r["amset_ok"]]
    loose = [r for r in rep["rows"] if r["amset_ok"] and r["amset_rel_err"] > 0.01]
    if loose:
        print("判断：AMSET 原式在宽松容差下接受了 %d 个掺杂偏差 > 1%% 的解（表中'宽松!'）—— 旧作业这些点的结果不可信，"
              "V132（amset_fermi_fix）起改用二分解，重跑即可。" % len(loose))
    if not d["efermi_finite"]:
        print("判断：本征费米能级是 NaN —— AMSET 的 get_fermi 会静默返回 NaN，下游 pinv 报 SVD 不收敛。")
    elif bad and d["states_in_gap_per_cell"] > 0.1 * min(abs(r["electrons_per_cell"]) for r in bad):
        print("判断：带隙内有插值假态（%.2e /胞），与失败点的载流子数（≥ %.2e /胞）可比 —— 低掺杂点受它支配。"
              % (d["states_in_gap_per_cell"], min(abs(r["electrons_per_cell"]) for r in bad)))
    elif bad and all(r["bisect_ef_eV"] is not None for r in bad):
        print("判断：DOS 正常而 AMSET 的贪心搜索失败（低温时被双精度舍入困住，见 amset_fermi_fix）—— 二分法都能解；"
              "V132 起作业自动改用二分解，不必改掺杂/温度。")
    elif bad:
        print("判断：有 %d 个点连二分法也无解 —— 浓度超出 DOS 能容纳的范围（能量窗口/带数不够？）。"
              % sum(r["bisect_ef_eV"] is None for r in bad))
    elif not loose:
        print("判断：所有 (掺杂, 温度) AMSET 都能解 —— 失败不在费米能级这一步，或与作业配置不同（--ir-fix？）。")
    if a.json:
        Path(a.json).write_text(json.dumps(rep, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
