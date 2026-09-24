#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_report.py -- 器件电-热自洽结果汇总报告

产物 step3_report/device_zt_report.json（判据 marker: "ZT_REPORT_DONE"）。
只读 S2 的 device_zt_solution.json，算几个派生指标并打印人读摘要。
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stepconf  # noqa: E402

OUTDIR = "step3_report"
STEP = "step3_report"
SOL = "step2_selftrap/device_zt_solution.json"
SPEC = {}


def _doping_scan(elec, s):
    """用当前器件的晶格 kappa 给各档掺杂估 ZT（近似：kappa_L 视为与掺杂无关），
    并标出最优档。PF = S^2*sigma 与 kappa_e 都用 transport.json 的胞口径——
    ZT 是比值，与电导口径无关；kappa_e 需换算成物理层值才能和 kappa_L 相加。"""
    scan = elec.get("doping_scan") or []
    if not scan:
        return None
    kl = float(s.get("kappa_lattice_W_mK") or 0.0)
    f = float(s.get("ke_cell_to_layer_factor") or 1.0)
    rows = []
    for r in scan:
        ke2d = float(r["kappa_e_W_mK"]) * f
        ktot = kl + ke2d
        T = float(r["T_K"])
        rows.append({"doping_cm3": r["doping_cm3"], "seebeck_uV_K": r["seebeck_uV_K"],
                     "sigma_S_m": r["sigma_S_m"], "PF_W_mK2": r["PF_W_mK2"],
                     "kappa_e_2D_W_mK": ke2d, "kappa_total_W_mK": ktot,
                     "zt_estimate": (float(r["PF_W_mK2"]) * T / ktot) if ktot > 0 else None})
    best = max(rows, key=lambda x: x["zt_estimate"] or -1)
    return {"T_K": scan[0]["T_K"], "kappa_lattice_W_mK": kl,
            "note": "近似：kappa_L 取当前器件值、视为与掺杂无关；ZT 与电导口径无关",
            "rows": rows, "best_doping_cm3": best["doping_cm3"],
            "best_zt_estimate": best["zt_estimate"]}


def main():
    cwd = os.getcwd()
    stepconf.load(SPEC, STEP, strict=False)
    out = os.path.join(cwd, OUTDIR)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(cwd, SOL), encoding="utf-8") as f:
        s = json.load(f)
    elec = {}
    ep = s.get("elec_props_source")
    if ep and os.path.isfile(ep):
        with open(ep, encoding="utf-8") as f:
            elec = json.load(f)

    kl = float(s["kappa_lattice_W_mK"])
    ke = float(s["kappa_e_W_mK"])
    ktot = kl + ke
    share = (ke / ktot) if ktot > 0 else float("nan")

    # 器件热阻：把 dT 归到 Q 上（有效面热阻），和 device-thermal 的解析串联值不同口径
    q = float(s["Q_joule_W_m2"])
    dT = float(s["dT_peak_K"] or 0.0)
    r_th_eff = (dT / q) if q > 0 else None

    rep = {
        "ZT_REPORT_DONE": True,
        "device_name": s.get("device_name"),
        "geometry": s.get("geometry"),
        "electrical": {
            "V_bias_V": s["V_bias_V"], "I_A": s["I_A"], "P_W": s["P_W"],
            "P_per_width_W_m": s["P_per_width_W_m"],
            "P_max_W": s.get("P_max_W"), "V_at_P_max_V": s.get("V_at_P_max_V"),
            "R_electrical_ohm": s["R_electrical_ohm"],
            "R_channel_ohm": s.get("R_channel_ohm"),
            "R_contact_ohm": s.get("R_contact_ohm"),
            "contact_share": s.get("contact_share"),
            "P_channel_W": s.get("P_channel_W"), "P_contact_W": s.get("P_contact_W"),
            "sheet_conductance_S": s["sheet_conductance_S"],
            "sigma_S_m": s["sigma_S_m"], "seebeck_V_K": s["seebeck_V_K"],
        },
        "thermal": {
            "Q_joule_W_m2": q, "dT_peak_K": dT,
            "T_eff_final_K": s["T_eff_final_K"],
            "T_channel_mean_K": s["T_channel_mean_K"],
            "T_channel_max_K": s["T_channel_max_K"],
            "kappa_lattice_W_mK": kl, "kappa_e_W_mK": ke,
            "kappa_total_W_mK": ktot, "electronic_share": share,
            "R_th_eff_m2K_per_W": r_th_eff,
            "tbc_W_m2K": s.get("tbc_W_m2K"), "tbc_source": s.get("tbc_source"),
        },
        "zt": {"zt_material": s["zt_material"],
               "doping_scan": _doping_scan(elec, s)},
        "limits": {
            "dt_limit_K": s.get("dt_limit_K"), "T_limit_K": s.get("T_limit_K"),
            "within_limit": s.get("within_limit"),
            "V_max_safe_V": s.get("V_max_safe_V"),
            "note": "超限不代表算错：只是工作点超出热/可靠性预算；V_max_safe 是扫描内仍满足判据的最大偏压",
        },
        "consistency": {
            "converged": s["converged"], "iterations": s["iterations"],
            "sor_iterations": s.get("sor_iterations"),
            "sigma_convention": s.get("sigma_convention"),
            "sigma_convention_source": s.get("sigma_convention_source"),
            "kappa_source_mode": s.get("kappa_source_mode"),
        },
        "iv_sweep": s.get("iv_sweep"),
        "upstream": {"kappa": s.get("upstream_kappa"), "elec": s.get("elec_props_source")},
    }
    p = os.path.join(out, "device_zt_report.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(rep, f, ensure_ascii=False, indent=2)

    print("=" * 68)
    print("device-zt 电-热自洽报告  %s" % rep["device_name"])
    print("=" * 68)
    g = rep["geometry"]
    ox_nm = float(s.get("oxide", {}).get("thickness_m", 0.0)) * 1e9
    print("  几何 : L=%.0f nm  W=%.2f um  %s on %s/%.0f nm on %s   T_amb=%.0f K"
          % (g["L_channel_m"] * 1e9, g["W_device_m"] * 1e6, g["channel_material"],
             s.get("oxide", {}).get("material", "oxide"), ox_nm,
             g["substrate_material"], g["T_amb_K"]))
    e, t = rep["electrical"], rep["thermal"]
    print("  电学 : V=%.4g V  I=%.4g A  P=%.4g W (%.4g W/mm)  R=%.4g ohm"
          % (e["V_bias_V"], e["I_A"], e["P_W"], e["P_per_width_W_m"] * 1e-3, e["R_electrical_ohm"]))
    print("         sigma=%.4g S/m  S=%.4g uV/K  G_sheet=%.4g S"
          % (e["sigma_S_m"], e["seebeck_V_K"] * 1e6, e["sheet_conductance_S"]))
    print("  热学 : Q=%.4g W/m2  dT_peak=%.4f K  T_eff=%.3f K  R_th_eff=%.4g m2K/W"
          % (t["Q_joule_W_m2"], t["dT_peak_K"], t["T_eff_final_K"], t["R_th_eff_m2K_per_W"] or 0.0))
    print("         kappa_L=%.4g  kappa_e=%.4g  kappa_tot=%.4g W/mK  （电子占 %.1f%%）"
          % (kl, ke, ktot, share * 100))
    print("  ZT   : S^2 sigma T/kappa = %.4g  @ T_eff=%.1f K" % (s["zt_material"], t["T_eff_final_K"]))
    if e.get("R_contact_ohm"):
        print("  接触 : R_contact=%.4g ohm（占总电阻 %.1f%%）；接触耗散 %.3g W 不进沟道"
              % (e["R_contact_ohm"], (e.get("contact_share") or 0) * 100,
                 e.get("P_contact_W") or 0.0))
    lim = rep["limits"]
    print("  上限 : dT<=%.4g K 且 T_max<=%.4g K -> %s；扫描内最大可用偏压 %s"
          % (lim["dt_limit_K"], lim["T_limit_K"],
             "满足" if lim["within_limit"] else "**超限**",
             ("%.3f V" % lim["V_max_safe_V"]) if lim["V_max_safe_V"] is not None
             else "无（全部超限）"))
    if not lim["within_limit"]:
        print("         [warn] 当前工作点 V=%.4g V 超出温度上限（dT=%.3f K, T_max=%.3f K）"
              % (e["V_bias_V"], t["dT_peak_K"], s.get("T_channel_max_K") or 0.0))
    print("  自洽 : converged=%s  迭代=%d 次  SOR=%s it"
          % (s["converged"], s["iterations"], s.get("sor_iterations")))
    print("  P_max: %.4g W @ V=%.3f V" % (s.get("P_max_W") or 0.0, s.get("V_at_P_max_V") or 0.0))
    ds = rep["zt"].get("doping_scan")
    if ds:
        print()
        print("  掺杂扫描（近似：kappa_L 固定 %.4g W/mK；ZT 与电导口径无关）" % ds["kappa_lattice_W_mK"])
        print("    %-12s %10s %10s %10s %10s" % ("掺杂(cm^-3)", "S(uV/K)", "sigma(S/m)", "PF(mW/mK2)", "ZT"))
        for r in ds["rows"]:
            mark = "  <== 最优" if abs(r["doping_cm3"] - ds["best_doping_cm3"]) < 1e-30 else ""
            print("    %-12.3g %10.3f %10.3g %10.4g %10.4g%s"
                  % (r["doping_cm3"], r["seebeck_uV_K"], r["sigma_S_m"],
                     r["PF_W_mK2"] * 1e3, r["zt_estimate"] or 0.0, mark))
        print("    最优掺杂 %.4g cm^-3 -> ZT ~ %.4g" % (ds["best_doping_cm3"], ds["best_zt_estimate"]))
    if not s["converged"]:
        print("  [warn] 自洽未收敛（达到 SC_MAXIT）——提高 SC_MAXIT 或放小 SC_RELAX 后重跑")
    print("[OK] %s" % p)


if __name__ == "__main__":
    main()
