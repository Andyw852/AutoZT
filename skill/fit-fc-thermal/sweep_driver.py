#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sweep_driver.py -- compute-node / login-node helper of the method sweep.

    python sweep_driver.py kappa-prep <variant dir>   after the fit (compute node)
    python sweep_driver.py result     <variant dir>   -> sweep_result.json (marker)
    python sweep_driver.py compare    <step3_sweep dir> [<out dir>]
                                                      -> fit_compare.json / .csv

kappa-prep finishes the S2_kappa recipe gen_step3_sweep.py left in
kappa/kappa_recipe.json with what only exists after the fit: fc2/fc3, the
phono3py YAML (copied or synthesised exactly like gen_step2_kappa), BORN/NAC and
the 2D flag.  It writes kappa/SKIP instead when there is nothing sensible to run
(no fc3, imaginary modes and SWEEP_KAPPA_IF_IMAG=false, SWEEP_KAPPA=false).
"""
import csv
import json
import shutil
import sys
from pathlib import Path


def _load(p):
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:
        return {}


def _write(p, obj):
    Path(p).write_text(json.dumps(obj, indent=2, ensure_ascii=False),
                       encoding="utf-8", newline="\n")


def kappa_prep(vdir):
    vdir = Path(vdir).resolve()
    kdir = vdir / "kappa"
    kdir.mkdir(exist_ok=True)
    skip = kdir / "SKIP"
    if skip.exists():
        skip.unlink()
    recipe = _load(kdir / "kappa_recipe.json")
    ps = _load(vdir / "phonon_summary.json")

    def _skip(why):
        skip.write_text(why + "\n", encoding="utf-8")
        print("[..] kappa skipped: %s" % why, flush=True)

    if not recipe:
        return _skip("SWEEP_KAPPA=false")
    if not ((vdir / "fc2.hdf5").is_file() and (vdir / "fc3.hdf5").is_file()):
        return _skip("the fit wrote no fc2.hdf5/fc3.hdf5")
    if ps.get("stable") is False and not recipe.get("kappa_if_imag"):
        return _skip("imaginary modes (min %s THz); SWEEP_KAPPA_IF_IMAG=false"
                     % ps.get("min_frequency_THz"))
    for f in ("fc2.hdf5", "fc3.hdf5"):
        shutil.copyfile(str(vdir / f), str(kdir / f))
    req = str(recipe.pop("nac_request", "auto") or "auto")
    born = vdir / "BORN"
    nac = born.is_file() if req in ("auto", "") else req in ("on", "true", "1", "yes")
    if nac and born.is_file():
        shutil.copyfile(str(born), str(kdir / "BORN"))
    elif nac:
        print("[WARN] NAC=on but the fit directory has no BORN -> nonac", flush=True)
        nac = False
    sys.path.insert(0, str(vdir))
    import gen_step2_kappa as g2      # same YAML lookup / synthesis as S2_kappa
    recipe["disp_yaml"] = g2._find_disp_yaml(vdir, kdir)
    recipe["nac"] = nac
    recipe["is_2d"] = bool(ps.get("is_2d", recipe.get("is_2d")))
    recipe.pop("kappa_if_imag", None)
    _write(kdir / "kappa_config.json", recipe)
    print("[OK] kappa/kappa_config.json (yaml=%s nac=%s 2D=%s mesh=%s)"
          % (recipe["disp_yaml"], nac, recipe["is_2d"], recipe.get("mesh")),
          flush=True)


def result(vdir):
    vdir = Path(vdir).resolve()
    var = _load(vdir / "sweep_variant.json")
    fs = _load(vdir / "fc_fit_summary.json") or _load(vdir / "phonon_summary.json")
    ks = _load(vdir / "kappa" / "kappa_summary.json")
    skip = vdir / "kappa" / "SKIP"
    row = dict(var)
    row.update({
        "SWEEP_DONE": True,
        "stable": fs.get("stable"),
        "min_frequency_THz": fs.get("min_frequency_THz"),
        "fit_rmse_relative": fs.get("fit_rmse_relative"),
        "fit_rmse_eV_per_A": fs.get("fit_rmse_eV_per_A"),
        "pheasy_relative_error": fs.get("pheasy_relative_error"),
        "pheasy_worst_force_correlation": fs.get("pheasy_worst_force_correlation"),
        "fit_quality_gate_failed": fs.get("fit_quality_gate_failed"),
        "kappa_done": bool(ks.get("KAPPA_DONE")),
        "kappa_skipped": (skip.read_text(encoding="utf-8").strip()
                          if skip.is_file() else None),
        "kappa_300K_xx_yy_zz": ks.get("kappa_300K_xx_yy_zz"),
        "kappa_inplane_300K": ks.get("kappa_inplane_300K"),
        "kappa_avg_300K": ks.get("kappa_avg_300K"),
        "kappa_principal_300K": ks.get("kappa_principal_300K"),
        "mesh": ks.get("mesh"),
        "mesh_converged": ks.get("mesh_converged"),
    })
    if not fs:
        # prep/fit/post ran under set -e, so a missing summary is a tool error:
        # leave no marker so the fanout step shows this variant as failed.
        sys.exit("[ERROR] %s: no fc_fit_summary.json -- the fit did not finish" % vdir)
    _write(vdir / "sweep_result.json", row)
    print("[DONE] sweep_result.json  %s: stable=%s  kappa_inplane_300K=%s"
          % (row.get("tag"), row["stable"], row["kappa_inplane_300K"]), flush=True)


COLS = ("tag", "engine", "method", "c3_shell", "c3_cutoff_A", "c2", "stable",
        "min_frequency_THz", "fit_rmse_relative", "pheasy_relative_error",
        "kappa_inplane_300K", "kappa_avg_300K", "kappa_xx_300K", "kappa_yy_300K",
        "kappa_zz_300K", "mesh", "mesh_converged", "kappa_skipped")


def compare(sweep_dir, out_dir=None):
    sweep_dir = Path(sweep_dir).resolve()
    out = Path(out_dir).resolve() if out_dir else sweep_dir
    out.mkdir(parents=True, exist_ok=True)
    plan = _load(sweep_dir / "sweep_plan.json")
    rows, missing = [], []
    tags = [v["tag"] for v in plan.get("variants", [])] or sorted(
        p.name for p in sweep_dir.glob("m-*") if p.is_dir())
    for tag in tags:
        r = _load(sweep_dir / tag / "sweep_result.json")
        if not r.get("SWEEP_DONE"):
            missing.append(tag)
            continue
        k = r.get("kappa_300K_xx_yy_zz") or [None, None, None]
        r["kappa_xx_300K"], r["kappa_yy_300K"], r["kappa_zz_300K"] = (list(k) + [None] * 3)[:3]
        rows.append(r)
    rows.sort(key=lambda r: (r.get("engine") or "", r.get("method") or "",
                             r.get("c3_shell") or 0))
    doc = {"COMPARE_DONE": True, "n_variants": len(tags), "n_done": len(rows),
           "missing": missing,
           "shell_convention": plan.get("convention"),
           "shell_table": plan.get("shell_table"),
           "supercell_safe_cutoff_A": plan.get("supercell_safe_cutoff_A"),
           "rows": rows}
    _write(out / "fit_compare.json", doc)
    with open(out / "fit_compare.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for r in rows:
            w.writerow(["" if r.get(c) is None else r.get(c) for c in COLS])
    print("%-34s %4s %7s %6s %10s %10s" % ("variant", "N", "c3(A)", "stable",
                                          "rmse_rel", "k_in300K"))
    for r in rows:
        def _f(x, fmt):
            return "-" if x is None else fmt % x
        print("%-34s %4s %7s %6s %10s %10s"
              % (r.get("tag"), r.get("c3_shell"), _f(r.get("c3_cutoff_A"), "%.3f"),
                 r.get("stable"), _f(r.get("fit_rmse_relative"), "%.4f"),
                 _f(r.get("kappa_inplane_300K"), "%.3f")))
    if missing:
        print("[WARN] %d variant(s) without sweep_result.json: %s"
              % (len(missing), ", ".join(missing)))
    print("[DONE] fit_compare.json / fit_compare.csv (%d/%d variants)"
          % (len(rows), len(tags)))
    return doc


def main(argv):
    if len(argv) < 3 or argv[1] not in ("kappa-prep", "result", "compare"):
        sys.exit("usage: sweep_driver.py kappa-prep|result <variant dir> | "
                 "compare <step3_sweep dir> [<out dir>]")
    if argv[1] == "kappa-prep":
        kappa_prep(argv[2])
    elif argv[1] == "result":
        result(argv[2])
    else:
        compare(argv[2], argv[3] if len(argv) > 3 else None)


if __name__ == "__main__":
    main(sys.argv)
