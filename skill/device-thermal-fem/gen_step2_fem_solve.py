#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_fem_solve.py -- 用可插拔后端求解 FEM 模型。

产出 step2_fem_solve/fem_summary.json（marker: "FEM_SOLVE_DONE"）+ T_field_fem.json。
后端由 step.conf 的 FEM_BACKEND 决定：skfem（默认）/ comsol（接口已留）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fem_common as FC  # noqa: E402
import stepconf  # noqa: E402

OUTDIR = "step2_fem_solve"
STEP = "step2_fem_solve"
PREV = "step1_fem_model"
SPEC = {
    "FEM_BACKEND": ("skfem", "str"),
    "FEM_DO_1D_CHECK": (True, "bool"),
    "FEM_NX": (60, "int"),
    "FEM_N_CHANNEL": (8, "int"),
    "FEM_N_INTERFACE": (4, "int"),
    "FEM_N_OXIDE": (200, "int"),
}


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    mp = os.path.join(cwd, PREV, "fem_model.json")
    if not os.path.isfile(mp):
        sys.exit("[ERROR] 找不到 %s，先跑 %s" % (mp, PREV))
    model = FC.load_json(mp)
    backend_name = str(conf["FEM_BACKEND"] or model.get("backend") or "skfem").strip().lower()
    backend = FC.get_backend(backend_name)

    mesh = {"nx": int(conf["FEM_NX"]), "n_channel": int(conf["FEM_N_CHANNEL"]),
            "n_interface": int(conf["FEM_N_INTERFACE"]), "n_oxide": int(conf["FEM_N_OXIDE"])}
    main_res = backend.solve(model, workdir=cwd, opts={"mesh": mesh})
    out = {"FEM_SOLVE_DONE": True, "backend": backend_name,
           "model": model, "main": main_res}
    if bool(conf["FEM_DO_1D_CHECK"]):
        out["ref_1d"] = backend.solve(model, workdir=cwd,
                                      opts={"mesh": mesh, "sink": False})
    FC.dump_json(os.path.join(cwd, OUTDIR, "fem_summary.json"), out)
    FC.dump_json(os.path.join(cwd, OUTDIR, "T_field_fem.json"),
                 {"y_nodes_m": main_res["y_nodes_m"], "T_mid_x_K": main_res["T_mid_x_K"],
                  "mesh": main_res["mesh"], "sink": main_res["sink"]})

    print("[OK] %s/fem_summary.json (+ T_field_fem.json)" % OUTDIR)
    print("    backend=%s ; dT_peak=%.6f K ; 源功率=%.6g W/m (名义 %.6g, 比 %.6f)"
          % (backend_name, main_res["dT_peak_K"], main_res["P_source_W_per_m"],
             main_res["P_nominal_W_per_m"], main_res.get("source_ratio") or 0.0))
    print("    mesh: %s" % main_res["mesh"])
    if "ref_1d" in out:
        r = out["ref_1d"]; c = r.get("analytic_check") or {}
        print("    1D 自检(sink=False): dT=%.6f K vs 解析 %.6f K -> 相对差 %s"
              % (r["dT_peak_K"], c.get("dT_analytic_K") or 0.0,
                 ("%.4f%%" % (100.0 * c["rel_err"])) if c.get("rel_err") is not None else "n/a"))


if __name__ == "__main__":
    main()
