#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_device_spec.py -- 把 step.conf 的器件参数固化成 device_spec.json

autozt 约定：run: gen 步骤的 cwd = 材料/device-thermal/ ；
产出 step1_device_spec/device_spec.json（判据 marker: "DEVICE_SPEC"）。
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stepconf  # noqa: E402

OUTDIR = "step1_device_spec"
STEP = "step1_device_spec"
SPEC = {
    "DEVICE_NAME": ("device", "str"),
    "CHANNEL_MATERIAL": ("MoS2", "str"),
    "CHANNEL_THICKNESS": (None, "float"),   # 可选：显式覆盖沟道厚度(m)；空=用材料库/上游
    "TBC_OVERRIDE": (None, "float"),        # 可选：显式覆盖界面热导 G (W/m^2K)；空=查材料库
    "OXIDE_MATERIAL": ("SiO2", "str"),
    "OXIDE_THICKNESS": (300e-9, "float"),
    "SUBSTRATE_MATERIAL": ("Si", "str"),
    "CONTACT_BC": ("sink", "str"),
    "L_CHANNEL": (200e-9, "float"),
    "W_DEVICE": (1e-6, "float"),
    "Q_JOULE": (1e8, "float"),
    "T_AMB": (300.0, "float"),
    "KAPPA_SOURCE": ("auto", "str"),
    "UPSTREAM_SKILL": ("kl-dft-cpu", "str"),
    "UPSTREAM_STEP": ("step6_kappa", "str"),
    "N_GRID_X": (41, "int"),
    "N_GRID_SUB": (60, "int"),
    "DO_TRANSIENT": (True, "bool"),
    "SOR_OMEGA": (1.7, "float"),
    "SOR_TOL": (1e-10, "float"),
    "SOR_MAXIT": (200000, "int"),
}


def main():
    conf = stepconf.load(SPEC, STEP, strict=False)
    out = os.path.join(os.getcwd(), OUTDIR)
    os.makedirs(out, exist_ok=True)
    spec = {
        "DEVICE_SPEC": True,
        "device_name": conf["DEVICE_NAME"],
        "channel_material": conf["CHANNEL_MATERIAL"],
        "channel_thickness_m": conf["CHANNEL_THICKNESS"],
        "tbc_override_W_m2K": conf["TBC_OVERRIDE"],
        "oxide_material": conf["OXIDE_MATERIAL"],
        "oxide_thickness_m": conf["OXIDE_THICKNESS"],
        "substrate_material": conf["SUBSTRATE_MATERIAL"],
        "contact_bc": conf["CONTACT_BC"],
        "L_channel_m": conf["L_CHANNEL"],
        "W_device_m": conf["W_DEVICE"],
        "Q_joule_W_m2": conf["Q_JOULE"],
        "T_amb_K": conf["T_AMB"],
        "kappa_source": conf["KAPPA_SOURCE"],
        "upstream_skill": conf["UPSTREAM_SKILL"],
        "upstream_step": conf["UPSTREAM_STEP"],
        "n_grid_x": conf["N_GRID_X"],
        "n_grid_sub": conf["N_GRID_SUB"],
        "do_transient": conf["DO_TRANSIENT"],
        "sor_omega": conf["SOR_OMEGA"],
        "sor_tol": conf["SOR_TOL"],
        "sor_maxit": conf["SOR_MAXIT"],
    }
    p = os.path.join(out, "device_spec.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(spec, f, ensure_ascii=False, indent=2)
    print("[OK] %s  channel=%s oxide=%s/%.0f nm L=%.0f nm W=%.2f um Q=%.3g W/m2 contact=%s"
          % (p, spec["channel_material"], spec["oxide_material"],
             spec["oxide_thickness_m"] * 1e9, spec["L_channel_m"] * 1e9,
             spec["W_device_m"] * 1e6, spec["Q_joule_W_m2"], spec["contact_bc"]))


if __name__ == "__main__":
    main()
