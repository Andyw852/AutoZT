#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kappa_driver.py -- step2_kappa compute-node driver: fc2/fc3 -> BTE -> kappa.

Nominal mode: one BTE run on the step's fc2/fc3 -> kappa_summary.json.

Scan mode (CUT3_KAPPA_SCAN, needs step1_fit to have written
cutoff_scan/cut3_<tag>/): run the BTE once per candidate fc3 cutoff, then pick
the cutoff with the three criteria - residual 1-SE, per-shell stability and the
kappa plateau - via cut3_select.select_cutoff, and promote the chosen kappa to
the top-level kappa_summary.json (with a cutoff_selection report).  This is the
"S1_fit determines the shells -> S2_kappa picks the cutoff by kappa" loop.
"""
import json
import sys
from pathlib import Path

import numpy as np


def _load_fc(ph3, out, fc2, fc3, enable_fc):
    from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
    ph3.fc2 = read_fc2_from_hdf5(str(out / fc2))
    if int(enable_fc or 3) >= 3:
        ph3.fc3 = read_fc3_from_hdf5(str(out / fc3))


def _flat_kappa(kappa):
    """Normalise the phono3py kappa array to (n_T, 3) [xx, yy, zz]."""
    k = np.asarray(kappa, float)
    if k.ndim == 3 and k.shape[-2:] == (3, 3):
        return np.stack([np.diag(k[i]) for i in range(k.shape[0])])
    if k.ndim == 1:
        k = k.reshape(1, -1)
    else:
        k = k.reshape(-1, k.shape[-1])
    if k.shape[1] >= 9:
        k = k[:, [0, 4, 8]]
    elif k.shape[1] >= 6:
        k = k[:, 0:3]
    elif k.shape[1] != 3:
        sys.exit("[ERROR] unexpected kappa shape %s" % (k.shape,))
    return k


def _summary(tc, cfg, source_fc):
    T = [float(x) for x in np.asarray(tc.temperatures, float).reshape(-1)]
    raw = [[float(v) for v in row] for row in _flat_kappa(tc.kappa)]
    n = min(len(T), len(raw))
    T, raw = T[:n], raw[:n]
    j = int(np.argmin(np.abs(np.asarray(T) - 300.0)))
    return {
        "KAPPA_DONE": True,
        "solver": "phono3py",
        "bte_method": str(cfg.get("bte_method") or "rta"),
        "mesh": cfg["mesh"],
        "nac": bool(cfg.get("nac")),
        "isotope": bool(cfg.get("isotope", True)),
        "temperatures": T,
        "kappa_xx_yy_zz": raw,
        "kappa_inplane_xx_yy": [[v[0], v[1]] for v in raw],
        "kappa_300K_xx_yy_zz": raw[j],
        "kappa_inplane_300K": 0.5 * (raw[j][0] + raw[j][1]),
        "source_fc": source_fc,
    }


def _run_bte(out, yaml_name, fc2, fc3, cfg, source_fc, write_kappa):
    import phono3py
    mesh = [int(x) for x in str(cfg["mesh"]).split()]
    ph3 = phono3py.load(str(out / yaml_name), produce_fc=False,
                        is_nac=bool(cfg.get("nac")), log_level=0)
    _load_fc(ph3, out, fc2, fc3, cfg.get("enable_fc"))
    ph3.mesh_numbers = mesh
    ph3.init_phph_interaction()
    ph3.run_thermal_conductivity(
        is_LBTE=(str(cfg.get("bte_method") or "rta").lower() == "lbte"),
        temperatures=[float(x) for x in cfg["temperatures"]],
        is_isotope=bool(cfg.get("isotope", True)),
        write_kappa=write_kappa, log_level=1)
    return _summary(ph3.thermal_conductivity, cfg, source_fc)


def _cut_of_tag(tag):
    return float(tag.replace("p", "."))


def _select(out, cfg, per_cut):
    """Three-criteria cutoff choice from the per-cut kappa + cutoff_scan.json."""
    import cut3_select
    scan = json.loads((out / "cutoff_scan.json").read_text(encoding="utf-8"))
    rows = {round(float(r["cut"]), 3): r
            for r in (scan.get("records") or []) if r.get("cut") is not None}
    recs = []
    for tag, s in per_cut:
        cut = _cut_of_tag(tag)
        r = rows.get(round(cut, 3), {})
        rel = 0.0
        for st in (r.get("shell_stats") or []):
            if (st.get("rel_std") is not None and st.get("shell") is not None
                    and float(st["shell"]) <= cut):
                rel = max(rel, float(st["rel_std"]))
        k = s.get("kappa_inplane_300K")
        recs.append({
            "cut": cut, "ratio": r.get("ratio"),
            "train_rel_err": r.get("train_rel_err"),
            "train_rel_se": r.get("train_rel_se"),
            "err_kind": r.get("err_kind"),
            "stable_upper_cut": r.get("stable_upper_cut"),
            "kappa": k, "kappa_err": (abs(k) * rel if k is not None else None),
        })
    chosen, rep = cut3_select.select_cutoff(
        recs,
        kappa_tol_pct=float(cfg.get("cut3_kappa_tol_pct") or 5.0),
        se_mult=float(cfg.get("cut3_cv_1se_mult") or 1.0),
        stability_thr=float(cfg.get("cut3_stability_thr") or 0.3),
        pick=str(cfg.get("cut3_pick") or "smallest"))
    rep["kappa_err_source"] = ("fc3 per-shell sigma/|mean| propagated (coarse; "
                               "one kappa run per cutoff)")
    (out / "cutoff_selection.json").write_text(
        json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8",
        newline="\n")
    print("[选截断] " + str(rep.get("reason", "")), flush=True)
    for r in rep.get("records", []):
        print("  c3=%-5.2f kappa=%-10s data=%s stable=%s cv=%s plateau=%s"
              % (r["cut"], r["kappa"], r["data_ok"], r["stable_ok"],
                 r["err_ok"], r["plateau_ok"]), flush=True)
    return chosen, rep, dict(per_cut)


def main():
    out = Path.cwd()
    cfg = json.loads((out / "kappa_config.json").read_text(encoding="utf-8"))
    yaml_name = cfg["disp_yaml"]

    if cfg.get("scan") and cfg.get("cut_scans"):
        per_cut = []
        for tag in cfg["cut_scans"]:
            print("[..] BTE @ cutoff %s" % _cut_of_tag(tag), flush=True)
            s = _run_bte(out, yaml_name, "cut3_%s/fc2.hdf5" % tag,
                         "cut3_%s/fc3.hdf5" % tag, cfg,
                         "cutoff_scan/cut3_%s" % tag, write_kappa=False)
            _cd = out / ("cut3_%s" % tag)
            _cd.mkdir(exist_ok=True)
            (_cd / "kappa_summary.json").write_text(
                json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
                newline="\n")
            per_cut.append((tag, s))
        chosen, rep, by_tag = _select(out, cfg, per_cut)
        if chosen is None:
            (out / "kappa_summary.json").write_text(
                json.dumps({"KAPPA_DONE": False, "cutoff_selection": rep},
                           indent=2, ensure_ascii=False), encoding="utf-8",
                newline="\n")
            sys.exit("[ERROR] 没有可用截断（数据量/稳定性不足）—— 见 cutoff_selection.json")
        tag = ("%.2f" % chosen).replace(".", "p")
        s = by_tag[tag]
        s["cutoff_selection"] = rep
        s["chosen_cutoff_A"] = chosen
        (out / "kappa_summary.json").write_text(
            json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
            newline="\n")
        print("[DONE] 选中 c3=%.2f A | 300K in-plane %.4f W/mK -> kappa_summary.json"
              % (chosen, s["kappa_inplane_300K"]), flush=True)
        return

    # nominal: single cutoff
    s = _run_bte(out, yaml_name, "fc2.hdf5", "fc3.hdf5", cfg,
                 cfg.get("source_fc"), write_kappa=True)
    (out / "kappa_summary.json").write_text(
        json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
        newline="\n")
    print("[DONE] kappa_summary.json | 300K in-plane %.4f W/mK "
          "(xx=%.4f yy=%.4f zz=%.4f)"
          % (s["kappa_inplane_300K"], s["kappa_300K_xx_yy_zz"][0],
             s["kappa_300K_xx_yy_zz"][1], s["kappa_300K_xx_yy_zz"][2]),
          flush=True)


if __name__ == "__main__":
    main()
