#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step4_report.py -- 器件热仿真汇总报告 + 温度场图。

产出 step4_report/device_report.json（marker: "REPORT_DONE"）；有 matplotlib 时另出 PNG。
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stepconf  # noqa: E402

OUTDIR = "step4_report"
STEP = "step4_report"
SPEC = {"MAKE_PLOT": (True, "bool")}


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    with open(os.path.join(cwd, "step1_device_spec", "device_spec.json"), encoding="utf-8") as f:
        spec = json.load(f)
    with open(os.path.join(cwd, "step2_props", "thermal_props.json"), encoding="utf-8") as f:
        props = json.load(f)
    with open(os.path.join(cwd, "step3_solve", "device_thermal_summary.json"), encoding="utf-8") as f:
        summ = json.load(f)
    with open(os.path.join(cwd, "step3_solve", "T_field.json"), encoding="utf-8") as f:
        tf = json.load(f)
    out = os.path.join(cwd, OUTDIR)
    os.makedirs(out, exist_ok=True)

    png, plot_err = None, None
    if conf["MAKE_PLOT"] is None or conf["MAKE_PLOT"]:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
            T = np.array(tf["T_K"])
            x = np.array(tf["x_nm"])
            y = np.array(tf["y_nm"])
            fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
            im = ax[0].pcolormesh(x, y, T.T, shading="auto", cmap="inferno")
            ax[0].set_yscale("symlog", linthresh=1.0)
            ax[0].set_xlabel("x along channel (nm)")
            ax[0].set_ylabel("depth y (nm)")
            ax[0].set_title("%s: T map (K)" % spec["device_name"])
            fig.colorbar(im, ax=ax[0])
            ax[1].plot(summ["T_x_mid_K"], x, "-")
            ax[1].set_xlabel("T (K)")
            ax[1].set_ylabel("x (nm)")
            ax[1].set_title("channel T(x) at surface")
            ax[2].plot(y, summ["T_profile_mid_K"], "o-", ms=3)
            ax[2].set_xlabel("depth y (nm)")
            ax[2].set_ylabel("T (K)")
            ax[2].set_title("vertical T(y) at mid-channel")
            ax[2].set_xscale("symlog", linthresh=1.0)
            fig.tight_layout()
            png = os.path.join(out, "device_thermal_Tmap.png")
            fig.savefig(png, dpi=130)
            plt.close(fig)
        except Exception as e:  # noqa: BLE001
            plot_err = str(e)
            print("[warn] 画图失败（不影响 JSON 报告）: %s" % e)

    report = {"REPORT_DONE": True, "device_name": spec["device_name"],
              "channel": props["channel"]["material"], "oxide": props["oxide"]["material"],
              "kappa_source": props["channel"]["source"], "tbc_W_m2K": props["tbc_W_m2K"],
              "dT_peak_K": summ["dT_peak_K"],
              "dT_peak_1D_network_K": summ["dT_peak_1D_network_K"],
              "R_th_eff_m2K_per_W": summ["R_th_eff_m2K_per_W"],
              "power_uW": summ["power_uW"], "tau_63_s": summ.get("tau_63_s"),
              "png": os.path.basename(png) if png else None, "plot_error": plot_err}
    with open(os.path.join(out, "device_report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print("[OK] %s/device_report.json%s" % (OUTDIR, ("  + " + os.path.basename(png)) if png else ""))
    print("    peak dT=%.4f K, power=%.3f uW, tau=%s" % (report["dT_peak_K"], report["power_uW"], report["tau_63_s"]))


if __name__ == "__main__":
    main()
