#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_sweep.py -- fit-fc-thermal method sweep (step3_sweep, fanout m-*).

Fits the SAME dataset with several engines/methods and third-order cutoffs and
runs the phono3py BTE on every fit, one independent job per variant:

    step3_sweep/m-<engine>-<method>-n<N>/      one variant = one sbatch
        fit_config.json submit.sh ...          = an S1_fit recipe (gen_step1_fit.main)
        kappa/kappa_config.json                = the S2_kappa recipe, finished on
                                                 the compute node after the fit
        sweep_result.json                      marker ("SWEEP_DONE": true)

S4_compare (login node) gathers every sweep_result.json into fit_compare.json /
fit_compare.csv.

Which variants (step3_sweep/step.conf, or `conf --set params.X=...`):
    SWEEP_METHODS    all | "phono3py:symfc pheasy:OLS,RIDGE hiphive:ridge"
                     (engine alone = that engine's default method,
                      engine:all = every method of the engine)
    SWEEP_C3_SHELLS  auto | "5" | "4 5 6" | "3-6"   -- N-th nearest neighbours,
                     the ShengBTE thirdorder.py "-n" convention (per atom, midpoint
                     between the N-th and (N+1)-th neighbour distance, max over
                     atoms), so "fc3 up to the 5th neighbours" reads the same as in
                     the literature.  auto = the largest N inside the supercell
                     safe cutoff.
    SWEEP_C2         all | N   -- fc2 range for pheasy/hiphive (phono3py/symfc
                     always fits the full fc2)
Every S1_fit / S2_kappa parameter (dataset, supercell, conda env, MESH, T range,
HIPHIVE_CELL, pheasy tuning ...) can be set in this step.conf too; the sweep only
overrides engine, method and cutoffs, and switches the per-variant cutoff scan
off (the shell list IS the scan).
"""
import json
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fc_common as fc
import gen_step1_fit as g1
import gen_step2_kappa as g2
import stepconf

OUTDIR = "step3_sweep"
STEP = "step3_sweep"

ENGINE_METHODS = {
    "phono3py": ("symfc", "alm"),
    "pheasy": g1.PHEASY_METHODS,
    "hiphive": g1.HIPHIVE_METHODS,
}
DEFAULT_METHOD = {"phono3py": "symfc", "pheasy": "RFE", "hiphive": "ridge"}

SWEEP_SPEC = {
    "SWEEP_METHODS":      ("all", "str"),
    "SWEEP_C3_SHELLS":    ("auto", "str"),
    "SWEEP_C2":           ("all", "str"),
    # shells closer than this (A) are one shell; thirdorder.py uses np.allclose
    # (~1e-5 relative).  Raise to ~0.05 to merge near-degenerate shells of a
    # slightly distorted cell (then N no longer matches thirdorder's count).
    "SWEEP_SHELL_TOL":    (1e-4, "float"),
    "SWEEP_KAPPA":        (True, "bool"),       # run the BTE on every fit
    "SWEEP_KAPPA_IF_IMAG": (False, "bool"),     # ...even when the fit has imaginary modes
}
SPEC = dict(g1.SPEC)
SPEC.update({k: v for k, v in g2.SPEC.items() if k not in SPEC})
SPEC.update(SWEEP_SPEC)

# per-variant values that must not be inherited from the sweep step.conf
FORCED = {"ENABLE_FC": 3, "CUT3_SCAN": "off", "CUT3_CANDIDATES": "off",
          "EXPORT_KL_BUNDLE": False}

HELPERS = ("sweep_driver.py", "gen_step2_kappa.py", "stepconf.py", "fc_common.py",
           # kappa-prep 在计算节点上算 2D 层厚归一化（gen_step2_kappa.two_d_norm）
           "dim_common.py", "thickness_2d.py", "vdw_radii.py")


class Overlay(object):
    """StepConf with per-variant overrides (read-only, same protocol)."""

    def __init__(self, base, over):
        self._b, self._o = base, dict(over)

    def __getitem__(self, k):
        return self._o[k] if k in self._o else self._b[k]

    def get(self, k, default=None):
        try:
            return self[k]
        except KeyError:
            return default

    def __contains__(self, k):
        return k in self._o or k in self._b

    @property
    def submit(self):
        return self._b.submit


def parse_methods(spec):
    """'all' | 'phono3py:symfc pheasy:OLS,RIDGE hiphive' -> [(engine, method)]."""
    s = str(spec or "all").strip()
    out = []
    toks = re.split(r"[\s;]+", s)
    for tok in toks:
        if not tok:
            continue
        if tok.lower() == "all":
            for e, ms in ENGINE_METHODS.items():
                out += [(e, m) for m in ms]
            continue
        eng, _, meths = tok.partition(":")
        eng = eng.strip().lower()
        if eng not in ENGINE_METHODS:
            sys.exit("[ERROR] SWEEP_METHODS: unknown engine %r (phono3py | pheasy | "
                     "hiphive)" % eng)
        allowed = ENGINE_METHODS[eng]
        if not meths:
            ms = [DEFAULT_METHOD[eng]]
        elif meths.strip().lower() == "all":
            ms = list(allowed)
        else:
            ms = [x.strip() for x in meths.split(",") if x.strip()]
        for m in ms:
            norm = m.upper() if eng == "pheasy" else m.lower()
            if norm not in allowed:
                sys.exit("[ERROR] SWEEP_METHODS: %s has no method %r (%s)"
                         % (eng, m, " | ".join(allowed)))
            out.append((eng, norm))
    seen, uniq = set(), []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    if not uniq:
        sys.exit("[ERROR] SWEEP_METHODS is empty")
    return uniq


def parse_shells(spec):
    """'auto' -> None; '5' / '4 5 6' / '3-6' / '4,6' -> sorted [N...]."""
    s = str(spec or "auto").strip().lower()
    if s in ("", "auto"):
        return None
    bad = ("[ERROR] SWEEP_C3_SHELLS=%r: use auto, N, 'N M ...' or 'N-M' (N >= 1)"
           % spec)
    vals = set()
    for tok in re.split(r"[\s,;]+", s):
        if not tok:
            continue
        m = re.fullmatch(r"(\d+)-(\d+)", tok)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            vals.update(range(min(a, b), max(a, b) + 1))
        elif tok.isdigit():
            vals.add(int(tok))
        else:
            sys.exit(bad)
    if not vals or min(vals) < 1:
        sys.exit(bad)
    return sorted(vals)


def variant_overrides(engine, method, c3, c2):
    """Engine/method/cutoff overrides for one variant.  c2 None = all."""
    o = dict(FORCED)
    o["FIT_ENGINE"] = engine
    if engine == "phono3py":
        o["FC_CALC"] = method
        o["FC3_CUTOFF"] = "%.4f" % c3
    elif engine == "pheasy":
        o["PHEASY_FIT_METHOD"] = method
        o["PHEASY_C3_CUTOFF"] = "%.4f" % c3
        o["PHEASY_C2_CUTOFF"] = "" if c2 is None else "%.4f" % c2
    else:
        o["HIPHIVE_FIT_METHOD"] = method
        o["HIPHIVE_CUTOFF3"] = c3
        o["HIPHIVE_CUTOFF2"] = "" if c2 is None else c2   # "" -> gen clamps to safe
    return o


def kappa_recipe(conf, vdir):
    """The S2_kappa recipe minus what only exists after the fit (fc2/fc3, the
    phono3py YAML, BORN, is_2d) -- sweep_driver.py kappa-prep fills those in."""
    mcfg = g2._mesh_cfg(conf, vdir)
    method = str(conf["BTE_METHOD"] or "rta").lower()
    if method not in ("rta", "lbte"):
        sys.exit("[ERROR] BTE_METHOD must be rta or lbte")
    t0, t1, dt = float(conf["T_MIN"]), float(conf["T_MAX"]), float(conf["T_STEP"])
    if dt <= 0 or t1 < t0:
        sys.exit("[ERROR] T_MIN/T_MAX/T_STEP 非法：%s %s %s" % (t0, t1, dt))
    temps, t = [], t0
    while t <= t1 + 1e-9:
        temps.append(round(t, 6))
        t += dt
    return dict(mcfg, **{
        "disp_yaml": None, "scan": False, "cut_scans": [],
        "cut3_kappa_tol_pct": float(conf["CUT3_KAPPA_TOL_PCT"]),
        "cut3_cv_1se_mult": float(conf["CUT3_CV_1SE_MULT"]),
        "cut3_stability_thr": float(conf["CUT3_STABILITY_THR"]),
        "cut3_pick": str(conf["CUT3_PICK"]),
        "temperatures": temps, "isotope": bool(conf["ISOTOPE"]),
        "nac_request": str(conf["NAC"] or "auto").strip().lower(),
        "bte_method": method, "enable_fc": 3, "source_fc": "..",
        "kappa_2d_thickness": str(conf["KAPPA_2D_THICKNESS"] or "vdw"),
    })


SWEEP_TAIL = """
# ---- fit-fc-thermal method sweep (gen_step3_sweep.py) ----------------------
# fit done above -> finish the kappa recipe on this node, run the BTE, then
# write sweep_result.json (the fanout marker) whatever the gate said.
python sweep_driver.py kappa-prep .
if [ -f kappa/kappa_config.json ] && [ ! -f kappa/SKIP ]; then
    ( cd kappa
      export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-${NCPU:-8}}
      export OMP_PROC_BIND=close OMP_PLACES=cores
      export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMPY_NUM_THREADS=1
      python kappa_driver.py kappa_config.json )
fi
python sweep_driver.py result .
"""


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)

    pairs = parse_methods(conf["SWEEP_METHODS"])
    want_n = parse_shells(conf["SWEEP_C3_SHELLS"])
    tol = float(conf["SWEEP_SHELL_TOL"] or 1e-4)

    src, _sig = g1.resolve_dataset(cwd, out, conf["FIT_INPUT_DIR"])
    pos = src / "POSCAR" if (src / "POSCAR").is_file() else src / "SPOSCAR"
    if not (src / "SPOSCAR").is_file():
        sys.exit("[ERROR] %s has no SPOSCAR -- the sweep needs it for the supercell "
                 "safe cutoff" % src)
    sc_cell, _ = fc.read_poscar_cell_frac(src / "SPOSCAR")
    safe = fc.supercell_safe_cutoff(sc_cell)

    # shells: enumerate far enough to cover the request (auto: up to the safe cutoff)
    nmax = max(want_n) if want_n else 30
    try:
        cuts, per_atom = fc.nth_neighbor_cutoffs(pos, nmax, tol=tol)
    except RuntimeError:
        if want_n:
            raise
        cuts, per_atom = fc.nth_neighbor_cutoffs(pos, 12, tol=tol)
    table = [{"n": i + 1, "cutoff_A": c, "within_safe": c <= safe + 1e-9}
             for i, c in enumerate(cuts)]
    print("[..] 第 N 近邻截断（thirdorder.py -n 约定，容差 %g A）；超胞安全截断 %.2f A"
          % (tol, safe))
    for r in table:
        print("     N=%-2d  c3=%.3f A%s" % (r["n"], r["cutoff_A"],
                                            "" if r["within_safe"] else "  > 安全截断"))
    if want_n is None:
        ok = [r["n"] for r in table if r["within_safe"]]
        if not ok:
            sys.exit("[ERROR] even the 1st-neighbour cutoff %.2f A exceeds the supercell "
                     "safe cutoff %.2f A -- enlarge the supercell" % (cuts[0], safe))
        want_n = [max(ok)]
        print("[..] SWEEP_C3_SHELLS=auto -> N=%d" % want_n[0])
    over = [n for n in want_n if cuts[n - 1] > safe + 1e-9]
    if over:
        print("[WARN] N=%s 的截断超过超胞安全截断 %.2f A（周期镜像重复计数）-> 跳过；"
              "要算更远的近邻请用更大的超胞" % (over, safe))
    shells = [n for n in want_n if n not in over]
    if not shells:
        sys.exit("[ERROR] no requested shell fits the supercell (safe cutoff %.2f A)"
                 % safe)

    c2s = str(conf["SWEEP_C2"] or "all").strip().lower()
    if c2s in ("", "all"):
        c2, c2_label = None, "all"
    else:
        try:
            n2 = int(c2s)
        except ValueError:
            sys.exit("[ERROR] SWEEP_C2 must be all or a neighbour number N")
        c2 = fc.nth_neighbor_cutoffs(pos, n2, tol=tol)[0][n2 - 1]
        if c2 > safe + 1e-9:
            print("[WARN] SWEEP_C2=%d (%.2f A) > 安全截断 -> 用 all" % (n2, c2))
            c2, c2_label = None, "all"
        else:
            c2_label = "n%d (%.3f A)" % (n2, c2)

    here = Path(__file__).resolve().parent
    for f in HELPERS + ("kappa_driver.py", "cut3_select.py"):
        if not (here / f).is_file():
            sys.exit("[ERROR] %s missing -- is it listed in gen_need?" % f)

    plan, kept, done = [], [], []
    for eng, meth in pairs:
        for n in shells:
            c3 = cuts[n - 1]
            tag = "m-%s-%s-n%d" % (eng, meth, n)
            vdir = out / tag
            rec = {"tag": tag, "engine": eng, "method": meth, "c3_shell": n,
                   "c3_cutoff_A": c3, "c2": c2_label,
                   "c2_cutoff_A": c2}
            plan.append(rec)
            res = vdir / "sweep_result.json"
            if res.is_file() and '"SWEEP_DONE": true' in res.read_text(
                    encoding="utf-8", errors="ignore"):
                done.append(tag)      # idempotent: never touch a finished variant
                continue
            print("\n[..] ==== %s  (c3 = %.3f A, N=%d; c2 = %s) ====" % (tag, c3, n, c2_label))
            if eng == "hiphive" and not str(conf["HIPHIVE_CELL"] or "").strip():
                print("[..] 提示：层状/薄原胞 hiphive 会因原胞太薄报错，需要时设 "
                      "HIPHIVE_CELL=supercell")
            vconf = Overlay(conf, variant_overrides(eng, meth, c3, c2))
            g1.main(conf=vconf, out=vdir, job_label="S3sweep")
            (vdir / "sweep_variant.json").write_text(
                json.dumps(rec, indent=2, ensure_ascii=False), encoding="utf-8",
                newline="\n")
            kdir = vdir / "kappa"
            kdir.mkdir(exist_ok=True)
            if bool(conf["SWEEP_KAPPA"]):
                kc = kappa_recipe(conf, vdir)
                kc["kappa_if_imag"] = bool(conf["SWEEP_KAPPA_IF_IMAG"])
                (kdir / "kappa_recipe.json").write_text(
                    json.dumps(kc, indent=2, ensure_ascii=False), encoding="utf-8",
                    newline="\n")
                for f in ("kappa_driver.py", "cut3_select.py"):
                    shutil.copyfile(str(here / f), str(kdir / f))
            for f in HELPERS:
                shutil.copyfile(str(here / f), str(vdir / f))
            with open(vdir / "submit.sh", "a", encoding="utf-8", newline="\n") as fh:
                fh.write(SWEEP_TAIL)
            kept.append(tag)

    (out / "sweep_plan.json").write_text(json.dumps({
        "convention": "N-th nearest neighbours, ShengBTE thirdorder.py -n "
                      "(per atom midpoint of N-th/(N+1)-th distance, max over atoms)",
        "shell_tol_A": tol, "supercell_safe_cutoff_A": round(safe, 4),
        "shell_table": table, "per_atom_shells_A": per_atom,
        "variants": plan}, indent=2, ensure_ascii=False),
        encoding="utf-8", newline="\n")
    print("\n[DONE] %s: %d 个变体（%d 方法 × %d 壳层 N=%s）；本次生成 %d，已完成跳过 %d。"
          "\n       每个 m-*/ 一个作业；全部完成后 S4_compare 汇总 fit_compare.json。"
          % (OUTDIR, len(plan), len(pairs), len(shells), shells, len(kept), len(done)))


if __name__ == "__main__":
    main()
