#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_solve.py -- 器件级稳态/瞬态传热 FVM 求解。

产出 step3_solve/device_thermal_summary.json（marker: "SOLVE_DONE"）+ T_field.json。
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import device_common as DC  # noqa: E402
import stepconf  # noqa: E402

OUTDIR = "step3_solve"
STEP = "step3_solve"
PREV_SPEC = "step1_device_spec"
PREV_PROPS = "step2_props"
SPEC = {
    "N_GRID_X": (41, "int"),
    "N_GRID_SUB": (60, "int"),
    "DO_TRANSIENT": (True, "bool"),
    "SOR_OMEGA": (1.7, "float"),
    "SOR_TOL": (1e-10, "float"),
    "SOR_MAXIT": (200000, "int"),
}


def _pick(v, fallback):
    return fallback if v is None else v


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    with open(os.path.join(cwd, PREV_SPEC, "device_spec.json"), encoding="utf-8") as f:
        spec = json.load(f)
    with open(os.path.join(cwd, PREV_PROPS, "thermal_props.json"), encoding="utf-8") as f:
        props = json.load(f)
    summary, T, dy, lay = DC.solve_device(
        spec, props,
        n_grid_x=int(_pick(conf["N_GRID_X"], spec.get("n_grid_x", 41))),
        n_grid_sub=int(_pick(conf["N_GRID_SUB"], spec.get("n_grid_sub", 60))),
        do_transient=bool(_pick(conf["DO_TRANSIENT"], spec.get("do_transient", True))),
        sor_omega=float(_pick(conf["SOR_OMEGA"], spec.get("sor_omega", 1.7))),
        sor_tol=float(_pick(conf["SOR_TOL"], spec.get("sor_tol", 1e-10))),
        sor_maxit=int(_pick(conf["SOR_MAXIT"], spec.get("sor_maxit", 200000))))
    summary["device_name"] = spec["device_name"]
    summary["channel_material"] = props["channel"]["material"]
    summary["oxide_material"] = props["oxide"]["material"]
    summary["tbc_W_m2K"] = props["tbc_W_m2K"]
    summary["kappa_source"] = props["channel"]["source"]
    out = os.path.join(cwd, OUTDIR)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "device_thermal_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    Nx, Ny = summary["Nx"], summary["Ny"]
    x_nm = [i * (spec["L_channel_m"] / max(1, Nx - 1)) * 1e9 for i in range(Nx)]
    with open(os.path.join(out, "T_field.json"), "w", encoding="utf-8") as f:
        json.dump({"Nx": Nx, "Ny": Ny, "x_nm": x_nm, "y_nm": summary["y_nm"],
                   "T_K": [[float(T[i, j]) for j in range(Ny)] for i in range(Nx)],
                   "layer_of": [int(x) for x in lay]}, f, ensure_ascii=False, indent=2)
    print("[OK] %s/device_thermal_summary.json (+ T_field.json)" % OUTDIR)
    print("    dT_peak=%.4f K (FVM) vs %.4f K (1D net, ratio %.3f) ; SOR %d it ; unknowns=%d"
          % (summary["dT_peak_K"], summary["dT_peak_1D_network_K"],
             summary["dT_ratio_FVM_over_1D"], summary["sor_iterations"], summary["unknowns"]))
    print("    power=%.3f uW @ Q=%.3g W/m2 ; R_th_eff=%.4g m2K/W ; contact=%s"
          % (summary["power_uW"], spec["Q_joule_W_m2"],
             summary["R_th_eff_m2K_per_W"], summary["contact_bc"]))
    if summary.get("tau_63_s") is not None:
        print("    transient tau_fast=%.4g s tau_RC=%.4g s tau_63=%.4g s (dt0=%.3g s, n=%d)"
              % (summary["tau_fast_s"], summary["tau_RC_s"], summary["tau_63_s"],
                 summary["transient_dt0_s"], summary["transient_nsteps"]))


if __name__ == "__main__":
    main()
