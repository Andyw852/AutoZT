#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_fem_model.py -- 复用 device-thermal 的器件定义，建 FEM 模型。

产出 step1_fem_model/fem_model.json（marker: "FEM_MODEL_DONE"）。
cwd 约定 = <材料>/device-thermal-fem/（与 device-thermal 一致）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fem_common as FC  # noqa: E402
import stepconf  # noqa: E402

OUTDIR = "step1_fem_model"
STEP = "step1_fem_model"
SPEC = {
    "FEM_BACKEND": ("skfem", "str"),
    "FEM_FROM_DEVICE_THERMAL_DIR": (None, "str"),
    "FEM_INTERFACE_LAYER_FRAC": (1e-3, "float"),
    "FEM_NX": (60, "int"),
    "FEM_N_CHANNEL": (8, "int"),
    "FEM_N_INTERFACE": (4, "int"),
    "FEM_N_OXIDE": (200, "int"),
}


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    backend = str(conf["FEM_BACKEND"] or "skfem").strip().lower()
    if backend not in FC.BACKENDS:
        sys.exit("[ERROR] FEM_BACKEND 只能是 %s，收到 %r"
                 % (" / ".join(FC.BACKENDS), backend))

    spec_p, prop_p = FC.find_artifacts(
        cwd, "device-thermal", conf["FEM_FROM_DEVICE_THERMAL_DIR"])
    if not spec_p or not prop_p:
        sys.exit(
            "[ERROR] 找不到 device-thermal 的产物：\n"
            "        device_spec.json   : %s\n"
            "        thermal_props.json : %s\n"
            "        先跑 device-thermal 到 S2_props 以上：\n"
            "          autozt -tt device-thermal -p <材料> -j S1_spec start\n"
            "          autozt -tt device-thermal -p <材料> -j S2_props start\n"
            "        或用 params.FEM_FROM_DEVICE_THERMAL_DIR 指定其目录。" % (spec_p, prop_p))

    spec = FC.load_json(spec_p)
    props = FC.load_json(prop_p)
    t_ox = float(props["oxide"]["thickness_m"])
    mesh = {
        "interface_layer_m": t_ox * float(conf["FEM_INTERFACE_LAYER_FRAC"]),
        "nx": int(conf["FEM_NX"]),
        "n_channel": int(conf["FEM_N_CHANNEL"]),
        "n_interface": int(conf["FEM_N_INTERFACE"]),
        "n_oxide": int(conf["FEM_N_OXIDE"]),
    }
    model = FC.build_model(spec, props, backend=backend, mesh=mesh)
    model["FEM_MODEL_DONE"] = True
    model["inputs"] = {"device_spec_json": spec_p, "thermal_props_json": prop_p}
    p = FC.dump_json(os.path.join(cwd, OUTDIR, "fem_model.json"), model)

    g, an = model["geometry"], model["analytic"]
    print("[OK] %s" % p)
    print("    backend=%s ; 几何 L=%.4g m t_ch=%.4g m t_ox=%.4g m"
          % (backend, g["L_m"], g["t_channel_m"], g["t_oxide_m"]))
    print("    界面 G=%.4g W/m2K ; Q=%.4g W/m2 ; T_amb=%.4g K ; 接触=%s"
          % (model["interface_G_W_m2K"], model["Q_W_m2"], model["T_amb_K"],
             model["contact_bc"]))
    print("    解析: R_series=%.4g -> dT_1D=%.4f K ; R_dist=%.4g -> dT_1D=%.4f K"
          % (an["R_series_m2K_W"], an["dT_series_K"],
             an["R_distributed_m2K_W"], an["dT_distributed_K"]))


if __name__ == "__main__":
    main()
