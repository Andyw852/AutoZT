#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot_kappa.py —— S4_plot：汇总 MD-κ，画 κ(T)，可选与 IFC 路线 κ / κ_sg 对比。"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc

OUTDIR = "step3_kappa"


def load_ifc_kappa(dirpath):
    """从 IFC 路线 step6_kappa 的 BTE.KappaTensorVsT_RTA 里抽 κ(T)（可选对比）。"""
    f = Path(dirpath) if dirpath else None
    if not f or not (f / "BTE.KappaTensorVsT_RTA").is_file():
        return None
    import numpy as np
    a = np.loadtxt(f / "BTE.KappaTensorVsT_RTA")
    return dict(T=a[:, 0].tolist(),
                kxx=a[:, 1].tolist(), kyy=a[:, 5].tolist(), kzz=a[:, 9].tolist())


def main():
    cwd = Path.cwd()
    step = cwd / OUTDIR
    ks = gc.load_params(step) if False else None
    import json
    summ = json.loads((step / "kappa_summary.json").read_text())
    temps = sorted(float(t) for t in summ["temperatures"])
    kxx = [summ["temperatures"]["%g" % t]["x"].get("kappa") for t in temps]
    kyy = [summ["temperatures"]["%g" % t]["y"].get("kappa") for t in temps]
    kzz = [summ["temperatures"]["%g" % t]["z"].get("kappa") for t in temps]
    payload = dict(KAPPA_PLOT_DONE=True, temperatures=temps,
                   kappa_xx=kxx, kappa_yy=kyy, kappa_zz=kzz)
    cfg_file = cwd / "step1_struct" / gc.PARAMS_JSON
    if cfg_file.is_file():
        cfg = gc.load_params(cwd / "step1_struct")
        ifc = load_ifc_kappa(cfg.get("IFC_KAPPA_RTA_DIR"))
        if ifc:
            payload["ifc_kappa"] = ifc
        if cfg.get("KAPPA_SG_AVG"):
            payload["kappa_sg_avg"] = cfg["KAPPA_SG_AVG"]
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.2, 3.6), dpi=150)
        ax.plot(temps, kxx, "o-", label="MD NEP HNEMD $\\kappa_{xx}$")
        ax.plot(temps, kyy, "s-", label="$\\kappa_{yy}$")
        ax.plot(temps, kzz, "^-", label="$\\kappa_{zz}$")
        if payload.get("ifc_kappa"):
            ax.plot(payload["ifc_kappa"]["T"], payload["ifc_kappa"]["kxx"],
                    "k--", lw=1, label="IFC RTA $\\kappa_{xx}$")
        if payload.get("kappa_sg_avg"):
            ax.axhline(payload["kappa_sg_avg"], color="grey", ls=":",
                       label="$\\kappa_{sg}$(IFC)")
        ax.set_xlabel("T (K)")
        ax.set_ylabel("$\\kappa$ (W/mK)")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(cwd / "kappa_T.png")
        payload["figure"] = "kappa_T.png"
    except Exception as exc:
        payload["figure_error"] = str(exc)
    gc.write_summary(cwd, "compare_summary.json", payload)


if __name__ == "__main__":
    main()
