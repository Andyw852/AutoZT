#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_selftrap.py -- 焦耳自热 <-> 温度依赖物性 的自洽迭代（device-zt 核心步）

autozt 约定：run: gen 步骤的 cwd = 材料/device-zt/ ；
产物 step2_selftrap/device_zt_solution.json（判据 marker: "ZT_SOLVE_DONE"）+ T_field.json。

复用：传热 FVM 用公共池 skill/_common/thermal/device_common.py 的 solve_device
（与 device-thermal 同一份代码）；热物性组装用 device_common.resolve_thermal_props。
本步只新增"电学 -> 焦耳热 -> 温度 -> 电学"这一层自洽循环。
"""
import copy
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import device_common as DC  # noqa: E402
import stepconf            # noqa: E402
import zt_common as ZC     # noqa: E402

OUTDIR = "step2_selftrap"
STEP = "step2_selftrap"
ELEC = "step1_elec/elec_props.json"

SPEC = {
    # ---- 器件几何/热边界（与 device-thermal 同名同义，便于对照）----
    "CHANNEL_MATERIAL": ("MoS2", "str"),
    "CHANNEL_THICKNESS": (None, "float"),
    "OXIDE_MATERIAL": ("SiO2", "str"),
    "OXIDE_THICKNESS": (300e-9, "float"),
    "SUBSTRATE_MATERIAL": ("Si", "str"),
    "CONTACT_BC": ("sink", "str"),
    "L_CHANNEL": (200e-9, "float"),
    "W_DEVICE": (1e-6, "float"),
    "T_AMB": (300.0, "float"),
    "TBC_OVERRIDE": (None, "float"),
    "KAPPA_SOURCE": ("auto", "str"),
    "UPSTREAM_SKILL": ("kl-dft-cpu", "str"),
    "UPSTREAM_STEP": ("step6_kappa", "str"),
    # ---- 电学驱动 ----
    "V_BIAS": (0.1, "float"),           # 自洽求解的工作偏压（V）
    "V_SWEEP_MIN": (0.0, "float"),
    "V_SWEEP_MAX": (0.5, "float"),
    "V_SWEEP_N": (6, "int"),
    "R_CONTACT_OHM": (0.0, "float"),    # 接触电阻（串联，Ω）；0=理想接触。接触耗散不进沟道
    "DT_LIMIT_K": (100.0, "float"),     # 温升上限判据（K）：超过就在扫描里标红
    "T_LIMIT_K": (450.0, "float"),      # 沟道绝对温度上限判据（K）
    # ---- 自洽 ----
    "SC_MAXIT": (40, "int"),
    "SC_TOL_K": (0.05, "float"),
    "SC_RELAX": (0.5, "float"),
    # ---- 网格/SOR ----
    "N_GRID_X": (41, "int"),
    "N_GRID_SUB": (60, "int"),
    "SOR_OMEGA": (1.7, "float"),
    "SOR_TOL": (1e-8, "float"),
    "SOR_MAXIT": (200000, "int"),
}


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    out = os.path.join(cwd, OUTDIR)
    os.makedirs(out, exist_ok=True)

    ep = os.path.join(cwd, ELEC)
    if not os.path.isfile(ep):
        sys.exit("[ERROR] 缺 %s，先跑 S1_elec" % ep)
    with open(ep, encoding="utf-8") as f:
        elec = json.load(f)

    T_amb = float(conf["T_AMB"])
    ch, ox, tbc, tbc_ref, up_info, ksrc = DC.resolve_thermal_props(
        cwd, conf["CHANNEL_MATERIAL"], conf["OXIDE_MATERIAL"], T_amb,
        up_skill=conf["UPSTREAM_SKILL"], up_step=conf["UPSTREAM_STEP"],
        kappa_source=conf["KAPPA_SOURCE"],
        thickness_override=conf["CHANNEL_THICKNESS"],
        tbc_override=conf["TBC_OVERRIDE"],
        oxide_thickness=conf["OXIDE_THICKNESS"])

    spec = {
        "device_name": elec.get("device_name", "device"),
        "channel_material": ch["material"], "oxide_material": ox["material"],
        "substrate_material": conf["SUBSTRATE_MATERIAL"], "contact_bc": conf["CONTACT_BC"],
        "L_channel_m": float(conf["L_CHANNEL"]), "W_device_m": float(conf["W_DEVICE"]),
        "T_amb_K": T_amb, "Q_joule_W_m2": 0.0,
        "R_contact_ohm": float(conf["R_CONTACT_OHM"]),
    }
    props = {"channel": ch, "oxide": ox, "tbc_W_m2K": tbc}

    tr = {"T_K": elec["T_K"], "sigma_S_m": elec["sigma_S_m"],
          "seebeck_V_K": elec["seebeck_V_K"], "kappa_e_W_mK": elec["kappa_e_W_mK"]}
    t_conv = float(elec["t_conv_m"])
    opts = {"maxit": int(conf["SC_MAXIT"]), "tol_K": float(conf["SC_TOL_K"]),
            "relax": float(conf["SC_RELAX"]),
            "n_grid_x": int(conf["N_GRID_X"]), "n_grid_sub": int(conf["N_GRID_SUB"]),
            "sor_omega": float(conf["SOR_OMEGA"]), "sor_tol": float(conf["SOR_TOL"]),
            "sor_maxit": int(conf["SOR_MAXIT"]), "do_transient": False,
            "dt_limit_K": float(conf["DT_LIMIT_K"]), "T_limit_K": float(conf["T_LIMIT_K"])}

    print("[..] 热物性: channel %s kx=%.4g ky=%.4g t=%.4g m | tbc=%.4g W/m2K"
          % (ch["material"], ch["kx_W_mK"], ch["ky_W_mK"], ch["thickness_m"], tbc))
    print("     kx 来源: %s" % ch.get("source"))
    if ch.get("kappa_note"):
        print("     kappa 口径: %s" % ch["kappa_note"])
    if ch.get("ky_source"):
        print("     跨面 kx: %s" % ch["ky_source"])
    print("     TBC 来源: %s" % tbc_ref)
    print("[..] 电学口径: %s  t_conv=%.4e m  doping=%.4g cm^-3"
          % (elec["sigma_convention"], t_conv, elec.get("doping_cm3", float("nan"))))

    # ---- 工作点自洽 ----
    net = {"k_lattice_W_mK": ch["kx_W_mK"], "thickness_m": ch["thickness_m"]}
    res, T, dy, lay = ZC.self_consistent_joule(spec, props, tr, t_conv, net,
                                              dict(opts, V_bias_V=float(conf["V_BIAS"])))
    print("[..] 自洽：%d 次迭代 converged=%s  T_eff=%.3f K  dT_peak=%.4f K"
          % (res["iterations"], res["converged"], res["T_eff_final_K"], res["dT_peak_K"]))
    print("     Q=%.4g W/m2  I=%.4g A  P=%.4g W  G_sheet=%.4g S"
          % (res["Q_joule_W_m2"], res["I_A"], res["P_W"], res["sheet_conductance_S"]))

    # ---- I-V / P-V 扫描（每个偏压各自自洽）----
    vmin, vmax, nv = float(conf["V_SWEEP_MIN"]), float(conf["V_SWEEP_MAX"]), int(conf["V_SWEEP_N"])
    vs = list(np.linspace(vmin, vmax, max(2, nv)))
    sweep = []
    for V in vs:
        r1, _, _, _ = ZC.self_consistent_joule(spec, props, tr, t_conv, net,
                                               dict(opts, V_bias_V=float(V)))
        sweep.append({"V_V": float(V), "I_A": r1["I_A"], "P_W": r1["P_W"],
                      "P_per_width_W_m": r1["P_per_width_W_m"],
                      "Q_joule_W_m2": r1["Q_joule_W_m2"],
                      "dT_peak_K": r1["dT_peak_K"], "T_eff_K": r1["T_eff_final_K"],
                      "T_channel_max_K": r1["T_channel_max_K"],
                      "sigma_S_m": r1["sigma_S_m"], "zt_material": r1["zt_material"],
                      "converged": r1["converged"], "within_limit": r1["within_limit"]})
    imax = int(np.argmax([s["P_W"] for s in sweep]))
    print("[..] 扫描 %d 点：P_max=%.4g W @ V=%.3f V" % (len(sweep), sweep[imax]["P_W"], sweep[imax]["V_V"]))
    safe = [s for s in sweep if s["within_limit"]]
    v_safe = max(s["V_V"] for s in safe) if safe else None
    print("[..] 温度上限判据 dT<=%.4g K 且 T_max<=%.4g K：最大可用偏压 = %s"
          % (float(conf["DT_LIMIT_K"]), float(conf["T_LIMIT_K"]),
             ("%.3f V" % v_safe) if v_safe is not None else "扫描范围内全部超限"))

    res["iv_sweep"] = sweep
    res["P_max_W"] = sweep[imax]["P_W"]
    res["V_at_P_max_V"] = sweep[imax]["V_V"]
    res["V_max_safe_V"] = (max(s["V_V"] for s in sweep if s["within_limit"])
                           if any(s["within_limit"] for s in sweep) else None)
    res["dt_limit_K"] = float(conf["DT_LIMIT_K"])
    res["T_limit_K"] = float(conf["T_LIMIT_K"])
    res["R_contact_ohm"] = float(conf["R_CONTACT_OHM"])
    res["device_name"] = spec["device_name"]
    res["geometry"] = {k: spec[k] for k in ("channel_material", "oxide_material",
                                            "substrate_material", "contact_bc",
                                            "L_channel_m", "W_device_m", "T_amb_K")}
    res["oxide"] = {k: ox[k] for k in ("material", "kx_W_mK", "ky_W_mK", "thickness_m")}
    res["tbc_W_m2K"] = tbc
    res["tbc_source"] = tbc_ref
    res["channel"] = {k: ch[k] for k in ("material", "kx_W_mK", "ky_W_mK", "thickness_m")}
    res["kappa_source_mode"] = ksrc
    res["upstream_kappa"] = up_info
    res["elec_props_source"] = ep
    res["sigma_convention"] = elec["sigma_convention"]
    res["sigma_convention_source"] = elec["sigma_convention_source"]

    p = os.path.join(out, "device_zt_solution.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    tf = {"T_K": np.asarray(T).tolist(), "dy_m": np.asarray(dy).tolist(),
          "layer_of_row": np.asarray(lay).tolist(),
          "Nx": int(T.shape[0]), "Ny": int(T.shape[1]),
          "L_channel_m": spec["L_channel_m"], "T_amb_K": T_amb}
    with open(os.path.join(out, "T_field.json"), "w", encoding="utf-8") as f:
        json.dump(tf, f, ensure_ascii=False, indent=2)
    print("[OK] %s" % p)


if __name__ == "__main__":
    main()
