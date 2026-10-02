#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_fit.py -- prepare the force-constant fitting job (S1_fit).

Runs on the login node in the material's skill directory (the step's own
directory is ./step1_fit).  It does three things and nothing heavy:

  1. Locate and sanity-check the displacement+force dataset (FIT_INPUT_DIR, or
     an automatic search of the usual sibling step directories).
  2. Parse step.conf into step1_fit/fit_config.json -- the single input the
     compute-node driver reads.
  3. Copy fc_fit_driver.py into step1_fit/ and render the engine-specific
     submit.sh, then let stepconf apply the [submit] overrides.

The fitting itself happens in the submitted job (fc_fit_driver.py prep|fit|post).
"""
import glob
import json
import os
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# 脱离 autozt 直接运行（python skill/fit-fc-thermal/gen_*.py）时，gen_need 里
# 来自 skill/_common 的 stepconf/dim_common/thickness_2d 等不在脚本旁边——
# 追加到 sys.path 末尾兜底（autozt 推送的同目录拷贝仍优先）。
for _d in ("_common/opt", "_common"):
    _p = Path(__file__).resolve().parent.parent / _d
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.append(str(_p))
import fc_common as fc
import stepconf

OUTDIR = "step1_fit"
STEP = "step1_fit"

# Dataset directory search order, relative to the skill directory (the gen cwd).
# Sibling skills keep their step directories next to ours under the material
# directory, so the cross-skill entries are the common case.
AUTO_DIRS = (
    "step4_disp", "step2_disp_force", "step2_disp", "step3_disp",
    "../kl-dft-cpu/step4_disp",
    "../kl-mlff-cpu/step2_disp_force",
    "../kl-mlff-gpu/step2_disp_force",
    "../phonon-dft-cpu/step2_disp",
    "../phonon-mlff-cpu/step2_disp_force",
    "../phonon-mlff-gpu/step2_disp_force",
    ".", "..", "../..",
)

SIGNATURES = (
    ("phono3py_params", ("phono3py_params.yaml",)),
    ("phono3py_disp_forces", ("phono3py_disp.yaml", "FORCES_FC3")),
    ("arrays", ("dataset_disps.npy", "dataset_forces.npy")),
    ("pkl", ("disp_matrix.pkl", "force_matrix.pkl")),
)

SPEC = {
    # ---- dataset ----
    "FIT_INPUT_DIR": ("auto", "str"),       # dataset dir: auto | path (skill-dir or absolute)
    "SUBTRACT_EQUILIBRIUM": (True, "bool"),  # subtract the equilibrium residual forces
    "EQ_FORCE_MAX": (0.2, "float"),          # max |F| of the equilibrium frame (eV/A); 0 = off
    "EQUILIBRIUM_FORCES_NPY": ("", "str"),   # optional explicit (natom,3) npy
    "COORDS": ("auto", "str"),               # auto | cartesian | fractional (displacement input)
    "SUPERCELL": ("", "str"),                # 对角 "n n n"；或一般矩阵 9 个数
                                             # "n11 n12 n13 n21 n22 n23 n31 n32 n33"
                                             # （行主序，与 phonopy/phono3py --dim 同义）
    # ---- engine ----
    # 拟合方法：FIT_METHODS 是总开关（auto = 用下面的 FIT_ENGINE + 该引擎的方法键）。
    #   写一个（如 pheasy:RIDGE）= 切换；写多个 / all = 全都算，第一个是主方法
    #   （S1 判据、S2 主结果），其余同一作业内串行拟合，S2 各算 κ 并出方法对比。
    "FIT_METHODS": ("auto", "str"),
    "FIT_ENGINE": ("pheasy", "str"),         # phono3py | pheasy | hiphive
    "ENABLE_FC": (3, "int"),                 # 2 | 3 (highest order to fit)
    "DIM": ("auto", "str"),                  # auto | 2d | 3d (NAC verdict only)
    # ---- shell / third-order cutoff determination (ported from kl-dft-cpu
    #      S4/S5, 2026-09-24 user 流程) ----
    # Candidates are the midpoints between adjacent neighbour shells of the
    # unit cell (never on a shell distance).  With CUT3_SCAN on and >= 2
    # candidates the pheasy engine refits each candidate, does a frame
    # bootstrap and reports, per shell, the scatter of |Phi^3| -- i.e. how far
    # the data can actually determine the fc3 (stable_upper_cut).  See
    # fc_fit_driver._shell_scan_pheasy and cutoff_scan.json.
    "CUT3_CANDIDATES": ("auto", "str"),      # auto | off | "3.6 4.2 4.8 ..."
    "CUT3_MAX": (7.0, "float"),              # largest radius scanned for shells (A)
    "CUT3_MIN_SHELLS": (2, "int"),           # skip the nearest-shell-only midpoint
    "CUT3_GAP_TOL": (0.05, "float"),         # shell clustering tolerance (A)
    "CUT3_MIN_GAP": (0.05, "float"),         # skip a midpoint in a split double shell
    # 默认 auto，与 templates/step1_fit/step.conf 一致（以前代码里是 off：项目 step.conf
    # 没写这个键就悄悄不扫，S2 只剩单点、画不出 κ 随截断的变化图）
    "CUT3_SCAN": ("auto", "str"),            # off | auto/on: refit each candidate
    "CUT3_BOOTSTRAP": (10, "int"),           # frame bootstrap per candidate (0 = off)
    "CUT3_STABILITY_THR": (0.3, "float"),    # sigma/|mean| below this = "determined"
    # ---- phono3py ----
    "FC_CALC": ("symfc", "str"),             # symfc | alm
    "FC3_CUTOFF": ("", "str"),               # fc3 cutoff in A; empty = no cutoff
    # ---- pheasy ----
    "PHEASY_FIT_METHOD": ("ALASSO", "str"),  # OLS|LASSO|ALASSO|RFE-OLS|RFE-OLS-TSQR|RIDGE|ARDR|RVM
    "PHEASY_BIN": ("pheasy-gpu", "str"),  # pheasy-gpu (CPU build removed)
    "PHEASY_C2_CUTOFF": ("", "str"),         # fc2 cutoff in A; empty = none
    "PHEASY_C3_CUTOFF": ("", "str"),         # fc3 cutoff in A; empty = none
    "NULL_SPACE_EPS": (0.001, "float"),
    "PHEASY_RASR": ("BHH", "str"),           # BH | H | BHH | none
    "PHEASY_STD": (False, "bool"),           # standardise the training data
    # '' = auto: on when the resident GPU LASSO is requested (that backend only
    # accepts a TwoLevelSM, and pheasy builds a dense/CSR matrix for LASSO unless
    # this is on).  true/false force it.
    "PHEASY_LASSO_TWOLEVEL": ("", "str"),    # '' | true | false
    # How many GPUs the job gets (--gres=gpu:N) and how many pheasy splits the
    # two-level sparse matvec across (PHEASY_GPU_SM_NGPU + SM_DEVICES=0..N-1).
    # Empty = keep whatever the submit template says.  The two must agree, so the
    # gen writes gres from this value unless [submit] already sets one -- a job
    # that holds fewer cards than pheasy is told to use fails at runtime, and one
    # that holds more just wastes them.
    "PHEASY_NGPU": ("", "str"),
    # pheasy's resident GPU backend for LASSO/ALASSO ('' = whatever the submit
    # template exports).  It requires the two-level sensing matrix, which the
    # driver turns on automatically when this is true.
    "PHEASY_GPU_LASSO_RESIDENT": ("", "str"),   # '' | true | false
    # pheasy's cross-validation iteration budget.  The CV stops each alpha at
    # min(max_iter, 400) iterations (resident LASSO) with tol = max(tol, 1e-3);
    # at that cap the CV curve is unresolved and alpha can end up on the grid
    # boundary.  Empty = pheasy's own defaults.
    "PHEASY_CV_MAX_ITER": ("", "str"),       # e.g. 2000
    "PHEASY_CV_TOL": ("", "str"),            # e.g. 1e-4
    # safe = only memory/threading/CV-grouping knobs; kl = the production solver
    # settings from submit_fit_pheasy.tpl, whose PHEASY_OLS_RIDGE=1e-4 measurably
    # degraded the local test fit (58% vs 1.3% relative error).  See
    # fc_fit_driver._pheasy_env for the measurements.
    "PHEASY_TUNING": ("safe", "str"),        # safe | kl
    "PHEASY_OLS_RIDGE": ("", "str"),         # override the OLS ridge (e.g. 0)
    "PHEASY_OLS_MAXITER": ("", "str"),       # OLS LSMR iteration cap; empty = pheasy default (5000)
    # RFE-OLS-TSQR only: how many features to keep.  pheasy's default since the
    # RFE -> RFE-OLS rename is 'aic' with n = the force-component rows; 'bic' is
    # the same with a heavier penalty, 'cv' reproduces RFE-OLS exactly (CV +
    # 1-SE).  Empty = pheasy's own default.
    "PHEASY_TSQR_CRITERION": ("", "str"),    # '' | aic | bic | cv
    # The LASSO grid and the RFE grid want different ranges; empty means "use
    # the method default" (see PHEASY_GRID_DEFAULTS), which reproduces the
    # values the pheasy author's own runs use.
    "PHEASY_CV": (5, "int"),
    "PHEASY_NMU": ("", "str"),
    # alpha exponents: integers for pheasy's argparse (see _int_str)
    "PHEASY_MU_MIN": (-8, "int"),
    "PHEASY_MU_MAX": ("", "str"),
    "PHEASY_MAX_ITER": ("", "str"),
    "PHEASY_TOL": ("", "str"),
    "PHEASY_SEED": (666666, "int"),
    # ---- hiphive ----
    "HIPHIVE_CUTOFF2": (6.0, "float"),       # fc2 cutoff in A (required)
    "HIPHIVE_CUTOFF3": (6.0, "float"),       # fc3 cutoff in A (ENABLE_FC=3)
    "HIPHIVE_FIT_METHOD": ("ridge", "str"),  # ols | ridge | lasso | ard | bayes
    "HIPHIVE_ALPHA": (1e-10, "float"),
    "HIPHIVE_SYMPREC": (1e-5, "float"),
    "HIPHIVE_ENFORCE_ASR": (True, "bool"),   # project rotational sum rules
    "HIPHIVE_N_CONFIGS": (0, "int"),         # 0 = use every frame
    # hiphive 参考胞："" = 原胞（默认）；"supercell" = 数据集超胞。原胞最小周期
    # 宽度 < 2*cutoff 时（尖锐/层状原胞）原胞建不起 cluster space，可设 supercell。
    "HIPHIVE_CELL": ("", "str"),
    # ---- output / gate ----
    "EXPORT_SHENGBTE": (True, "bool"),
    # 额外产出一份 kl_bundle/（fc2/fc3 + shengbte 力常数 + phonon_summary.json
    # + phono3py_disp.yaml 等），供 kl-dft-cpu 的 S5_fc 跨集群直接导入（见
    # kl skill 的 FC_IMPORT_DIR / push_paths）。体积小，默认开。
    "EXPORT_KL_BUNDLE": (True, "bool"),
    "FC3_LOAD_GB_LIMIT": (8.0, "float"),     # skip materialising fc3 above this (ShengBTE/RMSE)
    "BAND_POINTS": (51, "int"),
    "IMAG_THR": (0.10, "float"),             # imaginary-frequency threshold (THz)
    "FIT_RMSE_FRAMES": (10, "int"),          # frames for the force residual (also per
                                             # cutoff -> S2 RMSE panel); 0 = off
    # ---- environment ----
    "CONDA_SH": ("/public/home/<user>/miniconda3/etc/profile.d/conda.sh", "str"),
    "CONDA_ENV": ("atomate2_p_a", "str"),
}

ENGINES = ("phono3py", "pheasy", "hiphive")
# Canonical pheasy method names.  "RFE-OLS" is the current name of the OLS-based
# recursive-feature-elimination fitter; the pre-rename spelling "RFE" is NOT
# accepted -- write RFE-OLS.
PHEASY_METHODS = ("OLS", "LASSO", "ALASSO", "RFE-OLS", "RFE-OLS-TSQR", "RIDGE",
                   "ARDR", "RVM")


def normalize_pheasy_method(value):
    """Canonical spelling of a pheasy method name (upper case, no spaces)."""
    return str(value or "").strip().upper().replace(" ", "")
HIPHIVE_METHODS = ("ols", "ridge", "lasso", "ard", "bayes")

# Per-method defaults for the regularisation grid: (nmu, tol, max_iter).  A
# sparser penalty grid is enough for RFE because it refits on a shrinking
# active set, while plain LASSO needs to resolve the alpha that minimises the
# CV error.  These reproduce the settings the pheasy author uses.
PHEASY_GRID_DEFAULTS = {
    "LASSO": (40, 1e-5, 100000),
    "ALASSO": (40, 1e-5, 100000),
    "RFE-OLS": (5, 1e-3, 1000),
    "RFE-OLS-TSQR": (5, 1e-3, 1000),
}
PHEASY_MU_MAX_DEFAULT = {2: 0, 3: -5}   # upper bound of the alpha grid

# pheasy's --mu_min/--mu_max/--nmu/--max_iter are argparse ints: "-8.0" is
# rejected outright ("invalid int value"), which is how the first LASSO run
# failed.  Keep them as plain integer strings all the way to the command line.
def _int_str(value, key):
    try:
        return str(int(str(value).strip()))
    except (TypeError, ValueError):
        sys.exit("[ERROR] %s must be an integer (alpha exponent), got %r"
                 % (key, value))


def pheasy_grid(conf, method, enable_fc):
    """Resolve nmu / tol / max_iter / mu_max, honouring explicit settings."""
    nmu_d, tol_d, mi_d = PHEASY_GRID_DEFAULTS.get(
        method, PHEASY_GRID_DEFAULTS["RFE-OLS" if method.startswith("RFE") else "LASSO"])

    def pick(key, default, cast):
        v = conf[key]
        return default if v in (None, "") else cast(v)
    mu_max = pick("PHEASY_MU_MAX", PHEASY_MU_MAX_DEFAULT.get(enable_fc, -5), int)
    mu_min = int(conf["PHEASY_MU_MIN"])
    nmu = pick("PHEASY_NMU", nmu_d, int)
    tol = pick("PHEASY_TOL", tol_d, float)
    max_iter = pick("PHEASY_MAX_ITER", mi_d, int)
    if mu_max <= mu_min:
        sys.exit("[ERROR] PHEASY_MU_MAX=%d must exceed PHEASY_MU_MIN=%d"
                 % (mu_max, mu_min))
    return nmu, tol, max_iter, mu_max


# --------------------------------------------------------------------------
def detect_signature(d):
    d = Path(d)
    if (d / "phono3py_params.yaml").is_file():
        return "phono3py_params"
    if (d / "phono3py_disp.yaml").is_file() and (d / "FORCES_FC3").is_file():
        return "phono3py_disp_forces"
    if list(d.glob("disp-*/vasprun.xml")):
        return "vasprun"
    if (d / "dataset_disps.npy").is_file() and (d / "dataset_forces.npy").is_file():
        return "arrays"
    if (d / "disp_matrix.pkl").is_file() and (d / "force_matrix.pkl").is_file():
        return "pkl"
    if (d / "phono3py_disp.yaml").is_file():
        return "phono3py_disp"
    return None


def _yaml_read(path):
    try:
        import yaml
        return yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None


def supercell_diag(dataset_dir):
    """Diagonal supercell repetitions recorded in the dataset, if any."""
    for name in ("phono3py_params.yaml", "phono3py_disp.yaml"):
        p = Path(dataset_dir) / name
        if not p.is_file():
            continue
        y = _yaml_read(p)
        if not isinstance(y, dict):
            continue
        scm = y.get("supercell_matrix")
        if scm is None:
            continue
        try:
            rows = [[float(x) for x in row] for row in scm]
            if len(rows) == 3 and all(abs(rows[i][j]) < 1e-8
                                      for i in range(3) for j in range(3) if i != j):
                return " ".join(str(int(round(rows[i][i]))) for i in range(3)), rows
        except Exception:
            pass
    return None, None


def primitive_matrix(dataset_dir):
    """primitive_matrix recorded in the dataset YAML ('auto'/identity -> None).

    Phonopy/phono3py need it to rebuild the primitive cell from the unit cell;
    carrying it through keeps the gate and the plots faithful to the dataset."""
    for name in ("phono3py_params.yaml", "phono3py_disp.yaml"):
        p = Path(dataset_dir) / name
        if not p.is_file():
            continue
        y = _yaml_read(p)
        if not isinstance(y, dict) or y.get("primitive_matrix") is None:
            continue
        pm = y["primitive_matrix"]
        if isinstance(pm, str):
            return None if pm.strip().lower() in ("auto", "none", "") else None
        try:
            rows = [[float(x) for x in row] for row in pm]
        except Exception:
            return None
        if len(rows) != 3:
            return None
        if all(abs(rows[i][j] - (1.0 if i == j else 0.0)) < 1e-8
               for i in range(3) for j in range(3)):
            return None
        return rows
    return None


def resolve_dim(param, dataset_dir, material_dir):
    """DIM=auto -> inherit from a sibling stream's bookkeeping, else detect."""
    mode = str(param or "auto").lower()
    if mode in ("2d", "3d"):
        return mode
    for base in (dataset_dir, material_dir):
        for f in ("kl_params.txt", "workflow_method.txt", "klmlff_params.txt"):
            p = Path(base) / f if base else None
            if p and p.is_file():
                for ln in p.read_text(errors="ignore").splitlines():
                    if ln.strip().upper().startswith("DIM"):
                        v = ln.split("=", 1)[-1].strip().lower()
                        if v in ("2d", "3d"):
                            print("[OK] DIM=%s inherited from %s" % (v, p))
                            return v
    poscar = Path(dataset_dir) / "POSCAR" if dataset_dir else None
    if poscar and poscar.is_file():
        try:
            import dim_common
            dim, axis, _ = dim_common.detect_dimension(str(poscar))
            print("[OK] DIM=%s detected from %s (vacuum axis %s)"
                  % (dim, poscar, axis))
            return dim
        except Exception as e:
            print("[..] dimensional detection unavailable (%s)" % e)
    return "3d"


def resolve_dataset(cfg_dir, step_dir, want):
    """Return (abs_path, signature) for the dataset directory."""
    want = str(want or "auto").strip()
    if want and want.lower() != "auto":
        p = Path(want).expanduser()
        if not p.is_absolute():
            p = (Path(cfg_dir) / p)
        if not p.is_dir():
            sys.exit("[ERROR] FIT_INPUT_DIR=%s is not a directory" % p)
        sig = detect_signature(p)
        if not sig:
            sys.exit("[ERROR] FIT_INPUT_DIR=%s holds no recognised displacement "
                     "dataset (looked for %s)"
                     % (p, ", ".join(s[1][0] for s in SIGNATURES)))
        return p.resolve(), sig
    # The fixed list covers the skills that exist today; the glob covers any
    # other sibling skill whose displacement step is named differently, so a
    # new producer does not need this file to be edited.
    cands = [Path(cfg_dir) / rel for rel in AUTO_DIRS]
    matdir = Path(cfg_dir).resolve().parent
    # A dataset dropped at the material root (<material>/step4_disp, the layout
    # a hand-assembled dataset usually gets) comes right after the fixed list.
    for pat in ("step*disp*", "step*force*",
                "*/step*disp*", "*/step*force*", "*/step*fit*"):
        for hit in sorted(glob.glob(str(matdir / pat))):
            cands.append(Path(hit))
    # Never feed this skill its own output back: step1_fit/ holds the
    # normalised dataset + a synthesised phono3py_params.yaml from the previous
    # run, which ranks above the arrays and silently replaced the real dataset.
    own = Path(cfg_dir).resolve()
    seen, tried = set(), []
    for p in cands:
        key = str(p)
        if key in seen:
            continue
        seen.add(key)
        if not p.is_dir():
            continue
        rp = p.resolve()
        if rp == Path(step_dir).resolve() or (
                rp.parent == own and not any(
                    tok in rp.name for tok in ("disp", "force"))):
            continue
        sig = detect_signature(p)
        if sig:
            print("[OK] dataset auto-detected: %s (%s)" % (p.resolve(), sig))
            return p.resolve(), sig
        tried.append(str(p))
    sys.exit("[ERROR] no displacement+force dataset found.  Searched:\n  %s\n"
             "        Point FIT_INPUT_DIR at the dataset directory explicitly:\n"
             "        tf -tt fit-fc-thermal -p <material> -j step1_fit conf --set "
             "params.FIT_INPUT_DIR=../kl-dft-cpu/step4_disp"
             % "\n  ".join(tried))


ENGINE_METHODS = {
    "phono3py": ("symfc", "alm"),
    "pheasy": PHEASY_METHODS,
    "hiphive": HIPHIVE_METHODS,
}
DEFAULT_METHOD = {"phono3py": "symfc", "pheasy": "ALASSO", "hiphive": "ridge"}
_METHOD_KEY = {"phono3py": "FC_CALC", "pheasy": "PHEASY_FIT_METHOD",
               "hiphive": "HIPHIVE_FIT_METHOD"}


def parse_methods(spec):
    """'all' | 'ALASSO' | 'pheasy:OLS,RIDGE hiphive' | 'pheasy:all' -> [(engine, method)].

    A bare method name (ALASSO, ridge, symfc ...) is looked up across engines;
    an engine alone means its default method (phono3py:symfc, pheasy:ALASSO,
    hiphive:ridge)."""
    s = str(spec or "all").strip()
    out = []
    for tok in re.split(r"[\s;]+", s):
        if not tok:
            continue
        if tok.lower() == "all":
            for e, ms in ENGINE_METHODS.items():
                out += [(e, m) for m in ms]
            continue
        eng, _, meths = tok.partition(":")
        eng = eng.strip().lower()
        if eng not in ENGINE_METHODS and not meths:
            # bare method name: ALASSO / RIDGE / ridge / symfc ...  (normalise
            # case/space before matching; RFE-OLS must be written in full)
            _want = normalize_pheasy_method(eng).lower()
            hits = [(e, m) for e, ms in ENGINE_METHODS.items() for m in ms
                    if normalize_pheasy_method(m).lower() == _want]
            if len(hits) != 1:
                sys.exit("[ERROR] FIT_METHODS: %r is %s -- write engine:method "
                         "(e.g. pheasy:RIDGE, hiphive:ridge)"
                         % (tok, "ambiguous" if hits else "unknown"))
            out += hits
            continue
        if eng not in ENGINE_METHODS:
            sys.exit("[ERROR] FIT_METHODS: unknown engine %r (phono3py | pheasy | "
                     "hiphive)" % eng)
        allowed = ENGINE_METHODS[eng]
        if not meths:
            ms = [DEFAULT_METHOD[eng]]
        elif meths.strip().lower() == "all":
            ms = list(allowed)
        else:
            ms = [x.strip() for x in meths.split(",") if x.strip()]
        for m in ms:
            norm = normalize_pheasy_method(m) if eng == "pheasy" else m.lower()
            if norm not in allowed:
                sys.exit("[ERROR] FIT_METHODS: %s has no method %r (%s)"
                         % (eng, m, " | ".join(allowed)))
            out.append((eng, norm))
    seen, uniq = set(), []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    if not uniq:
        sys.exit("[ERROR] FIT_METHODS is empty")
    return uniq


def method_overrides(engine, method):
    return {"FIT_ENGINE": engine, _METHOD_KEY[engine]: method, "FIT_METHODS": "auto"}


def method_tag(engine, method):
    return "m-%s-%s" % (engine, method)


class Overlay(object):
    """StepConf with per-method overrides (read-only, same protocol)."""

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


def _cfg_method(cfg):
    eng = str(cfg.get("engine") or "")
    return eng, {"phono3py": cfg.get("fc_calc"), "pheasy": cfg.get("pheasy_method"),
                 "hiphive": cfg.get("hiphive_fit_method")}.get(eng)


METHODS_TAIL = """
# ---- FIT_METHODS: the other fitting methods, serially, in methods/<tag>/ ------
# Non-fatal: a failing method is logged and skipped; the primary fit above owns
# the S1 gate.  Each method directory is a complete S1 recipe (fit_config.json).
set +e
for d in %s; do
    ( cd "methods/$d" && python fc_fit_driver.py prep fit_config.json \\
        && python fc_fit_driver.py fit fit_config.json \\
        && python fc_fit_driver.py post fit_config.json ) > "methods/$d/run.log" 2>&1 \\
        || echo "[WARN] method $d failed -- see methods/$d/run.log"
done
set -e
"""


def main(conf=None, out=None, job_label="S1fit"):
    """S1_fit gen.

    Called without arguments it is the plain S1_fit gen: FIT_METHODS picks the
    method(s) -- auto = FIT_ENGINE + its method key; one method = switch; a list
    or `all` = the first is the primary (step1_fit/ itself, the S1 gate), the
    others go to step1_fit/methods/<tag>/ and run in the same job.
    conf/out/job_label let gen_step3_sweep.py reuse the recipe for one sweep
    variant (out = step3_sweep/m-<tag>/)."""
    cwd = Path.cwd()
    if conf is None:
        conf = stepconf.load(SPEC, STEP, strict="warn")
    if out is not None:
        return _gen_one(conf, Path(out), job_label)
    out = cwd / OUTDIR
    spec = str(conf.get("FIT_METHODS") or "auto").strip()
    methods = None if spec.lower() in ("", "auto") else parse_methods(spec)
    if not methods:
        cfg = _gen_one(conf, out, job_label)
        entries = [dict(zip(("engine", "method"), _cfg_method(cfg)), dir=".",
                        primary=True)]
    else:
        print("[..] FIT_METHODS=%s -> %d 个方法：%s（主方法 %s）"
              % (spec, len(methods), ", ".join("%s:%s" % m for m in methods),
                 "%s:%s" % methods[0]), flush=True)
        cfg = _gen_one(Overlay(conf, method_overrides(*methods[0])), out, job_label)
        entries = [{"engine": methods[0][0], "method": methods[0][1], "dir": ".",
                    "primary": True}]
        for eng, meth in methods[1:]:
            tag = method_tag(eng, meth)
            print("\n[..] ==== 方法 %s ====" % tag, flush=True)
            _gen_one(Overlay(conf, method_overrides(eng, meth)),
                     out / "methods" / tag, job_label)
            entries.append({"engine": eng, "method": meth,
                            "dir": "methods/" + tag, "primary": False})
        extra = [e["dir"].split("/", 1)[1] for e in entries if not e["primary"]]
        if extra:
            with open(out / "submit.sh", "a", encoding="utf-8", newline="\n") as fh:
                fh.write(METHODS_TAIL % " ".join(extra))
            print("[..] 其余 %d 个方法与主方法同一作业串行运行（多方法耗时累加；"
                  "大批量对比更适合 S3_sweep 并行）" % len(extra), flush=True)
    for e in entries:
        e["tag"] = method_tag(e["engine"], e["method"])
    (out / "methods.json").write_text(
        json.dumps({"methods": entries}, indent=2, ensure_ascii=False),
        encoding="utf-8", newline="\n")
    return cfg


def _gen_one(conf, out, job_label="S1fit"):
    """One S1 recipe (one engine/method) into `out`."""
    cwd = Path.cwd()
    out.mkdir(parents=True, exist_ok=True)

    engine = str(conf["FIT_ENGINE"] or "pheasy").lower()
    if engine not in ENGINES:
        sys.exit("[ERROR] FIT_ENGINE must be one of %s" % " | ".join(ENGINES))
    enable = int(conf["ENABLE_FC"] or 3)
    if enable not in (2, 3):
        sys.exit("[ERROR] ENABLE_FC must be 2 or 3 (this skill fits fc2/fc3)")
    p_bin = str(conf["PHEASY_BIN"] or "pheasy-gpu").lower()
    if p_bin != "pheasy-gpu":
        sys.exit("[ERROR] PHEASY_BIN must be pheasy-gpu (the CPU build was removed)")
    p_ngpu = str(conf["PHEASY_NGPU"] or "").strip()
    if p_ngpu:
        try:
            _n = int(p_ngpu)
        except ValueError:
            sys.exit("[ERROR] PHEASY_NGPU must be an integer (GPU count), got %r"
                     % p_ngpu)
        if _n < 1:
            sys.exit("[ERROR] PHEASY_NGPU must be >= 1, got %d" % _n)
        p_ngpu = str(_n)
    fc_calc = str(conf["FC_CALC"] or "symfc").lower()
    if fc_calc not in ("symfc", "alm"):
        sys.exit("[ERROR] FC_CALC must be symfc or alm")
    p_method = normalize_pheasy_method(conf["PHEASY_FIT_METHOD"] or "ALASSO")
    if engine == "pheasy" and p_method not in PHEASY_METHODS:
        sys.exit("[ERROR] PHEASY_FIT_METHOD must be one of %s"
                 % " | ".join(PHEASY_METHODS))
    h_method = str(conf["HIPHIVE_FIT_METHOD"] or "ridge").lower()
    if engine == "hiphive" and h_method not in HIPHIVE_METHODS:
        sys.exit("[ERROR] HIPHIVE_FIT_METHOD must be one of %s"
                 % " | ".join(HIPHIVE_METHODS))
    coords = str(conf["COORDS"] or "auto").lower()
    if coords not in ("auto", "cartesian", "fractional"):
        sys.exit("[ERROR] COORDS must be auto, cartesian or fractional")

    src, sig = resolve_dataset(cwd, out, conf["FIT_INPUT_DIR"])

    # Fail fast on the one dataset/engine combination that cannot work: a
    # disp-*/vasprun.xml layout has neither forces embedded in a phono3py YAML
    # nor a FORCES_FC3 file, so the phono3py engine would only discover it at
    # fit time ("no forces available for the phono3py engine").  Guard exactly
    # this layout: a phono3py_params.yaml / phono3py_disp.yaml + FORCES_FC3 /
    # npy / pkl dataset still goes through untouched.
    if sig == "vasprun" and engine == "phono3py":
        sys.exit(
            "[ERROR] FIT_ENGINE=phono3py cannot use a disp-*/vasprun.xml "
            "dataset (%s).\n"
            "        phono3py needs the forces either embedded in a phono3py "
            "YAML\n"
            "        (phono3py_params.yaml / phono3py_disp.yaml) or in a "
            "FORCES_FC3 file,\n"
            "        and this layout carries neither -- the fit would abort "
            "later with\n"
            "        'no forces available for the phono3py engine'.\n"
            "        Use a regression engine (it reads the vasprun forces):\n"
            "          tf -tt fit-fc-thermal -p <material> -j step1_fit conf "
            "--set params.FIT_ENGINE=pheasy\n"
            "          tf -tt fit-fc-thermal -p <material> -j step1_fit conf "
            "--set params.FIT_ENGINE=hiphive\n"
            "        or point FIT_INPUT_DIR at a dataset whose YAML embeds the "
            "forces\n"
            "        / that ships FORCES_FC3:\n"
            "          tf -tt fit-fc-thermal -p <material> -j step1_fit conf "
            "--set params.FIT_INPUT_DIR=<yaml-or-FORCES_FC3 dataset>\n"
            "        (dataset: %s)" % (src, src))

    # Random-displacement (type-2) datasets are required by the regression
    # engines; finite-difference datasets only carry enough information for the
    # phono3py rebuild.
    sc_str, sc_rows = supercell_diag(src)
    if not sc_str:
        # A dataset without a phonopy YAML (the npy / pkl layouts) still ships
        # POSCAR + SPOSCAR, so the repetitions are recoverable from the edge
        # ratio -- pheasy needs them explicitly.
        _miss = [n for n in ("POSCAR", "SPOSCAR") if not (src / n).is_file()]
        if _miss:
            sys.exit("[ERROR] dataset %s has no phonopy YAML and lacks %s -- the "
                     "supercell cannot be deduced.\n"
                     "        Add POSCAR (unit cell) + SPOSCAR (supercell) next to the "
                     "arrays, or point FIT_INPUT_DIR at a complete dataset:\n"
                     "          autozt -tt fit-fc-thermal -p <material> -j step1_fit conf "
                     "--set params.FIT_INPUT_DIR=<dir>" % (src, " + ".join(_miss)))
        reps = fc.supercell_reps(src / "POSCAR", src / "SPOSCAR")
        if reps:
            sc_str = " ".join(str(x) for x in reps)
            sc_rows = [[float(reps[0]), 0, 0], [0, float(reps[1]), 0],
                       [0, 0, float(reps[2])]]
            print("[OK] supercell %s deduced from the POSCAR/SPOSCAR edge ratio"
                  % sc_str)
    if str(conf["SUPERCELL"] or "").strip():
        # 3 个数=对角扩胞；9 个数=3×3 矩阵（行主序，与 phono3py --dim 同义，
        # hiphive/pheasy 直接吃这个 3×3；phono3py 那一路按原样透传 --dim）
        _sc_vals = [int(x) for x in str(conf["SUPERCELL"]).split()]
        if len(_sc_vals) == 9:
            sc_str = " ".join(str(x) for x in _sc_vals)
            sc_rows = [[float(x) for x in _sc_vals[0:3]],
                       [float(x) for x in _sc_vals[3:6]],
                       [float(x) for x in _sc_vals[6:9]]]
        elif len(_sc_vals) == 3:
            sc_str = " ".join(str(x) for x in _sc_vals)
            sc_rows = [[float(_sc_vals[0]), 0.0, 0.0],
                       [0.0, float(_sc_vals[1]), 0.0],
                       [0.0, 0.0, float(_sc_vals[2])]]
        else:
            sys.exit("[ERROR] SUPERCELL 要写 3 个（对角，如 \"3 3 3\"）或 9 个"
                     "（3×3 矩阵，行主序，如 \"2 1 0 -1 2 0 0 0 1\"）整数，"
                     "收到 %r" % conf["SUPERCELL"])
    if engine in ("pheasy", "hiphive") and sig in ("vasprun", "phono3py_disp_forces"):
        print("[WARN] the dataset looks like a finite-displacement set; "
              "%s fits random-displacement (rattled) data.  The driver will stop "
              "unless the YAML actually carries a type-2 displacement array."
              % engine)
    if engine == "pheasy" and not sc_str:
        sys.exit("[ERROR] FIT_ENGINE=pheasy needs the supercell repetitions; set "
                 "SUPERCELL=\"n n n\" in step.conf (the dataset records none)")

    p_nmu, p_tol, p_max_iter, p_mu_max = pheasy_grid(conf, p_method, enable)
    dim = resolve_dim(conf["DIM"], src, src.parent)

    # ---- shell / third-order cutoff candidates (ported from kl-dft-cpu S4/S5) ----
    #   Cheap and engine independent: enumerate the neighbour shells of the unit
    #   cell and take the midpoints between adjacent shells as cutoff candidates.
    #   The compute-node driver then decides (pheasy: refit + frame bootstrap)
    #   how far the data can determine the fc3.  A failure here degrades to "no
    #   candidates" -- it must never block a fit that would otherwise run.
    _pos = src / "POSCAR"
    if not _pos.is_file():
        _pos = src / "SPOSCAR"
    _cands, _shells, _cnote = [], [], "no POSCAR/SPOSCAR in the dataset"
    try:
        _cands, _shells, _cnote = fc.resolve_cut3_candidates(
            _pos, conf["CUT3_CANDIDATES"], cut3_max=float(conf["CUT3_MAX"]),
            tol=float(conf["CUT3_GAP_TOL"]), min_shells=int(conf["CUT3_MIN_SHELLS"]),
            min_gap=float(conf["CUT3_MIN_GAP"]))
    except Exception as _e:                              # noqa: BLE001
        _cnote = "shell enumeration failed: %s" % _e
    if not _shells and _pos.is_file():
        # explicit / off candidates: still record the shells, S2's kappa-vs-cutoff
        # figure labels every cutoff with the number of shells inside it
        try:
            _shells = fc.neighbor_shells(_pos, r_max=float(conf["CUT3_MAX"]),
                                         tol=float(conf["CUT3_GAP_TOL"]))[0]
        except Exception:                                # noqa: BLE001
            _shells = []
    print("[..] 壳层/截断候选：%s" % _cnote)
    _safe, _nsc, _eff = None, 0, {}
    try:
        _sc_cell, _sc_frac = fc.read_poscar_cell_frac(src / "SPOSCAR")
        _nsc = len(_sc_frac)
        _safe = fc.supercell_safe_cutoff(_sc_cell)
    except Exception:
        _safe = None
    if _safe is not None:
        _over = [c for c in _cands if c > _safe + 1e-9]
        if _over:
            print("[WARN] 候选截断 %s 超过超胞安全截断 %.2f Å（周期镜像会重复计数），"
                  "已剔除" % (_over, _safe))
        _cands = [c for c in _cands if c <= _safe + 1e-9]
        # 补齐「安全上限」候选：可能存在壳层落在最后一个中点候选之上、但在安全
        # 截断之内（例 In2MnSe4：5.10 Å 壳层，安全上限 5.165 Å，而它到下一壳层
        # 5.41 的中点 5.258 Å 已越界）——只取中点会把这条壳层整条漏掉。
        # 仅在 CUT3_CANDIDATES=auto 时自动补，显式列表按用户给的跑。
        if str(conf["CUT3_CANDIDATES"]).strip().lower() == "auto" \
                and _shells:
            _top = max([c for c in _cands], default=0.0)
            _skipped = [s for s in _shells
                        if _top + 1e-9 < s <= _safe + 1e-9]
            _extra = float(int(_safe * 100)) / 100.0   # 2 位小数，仍 < safe
            if _skipped and _extra > _top + 1e-9:
                _cands = sorted(set(_cands + [_extra]))
                print("[..] 补齐安全上限候选 %.2f Å（含安全范围内被中点方案漏掉的"
                      "壳层 %s）" % (_extra,
                                 ", ".join("%.2f" % s for s in _skipped)))
        print("[..] 超胞安全截断 = %.2f Å → 有效候选 %s" % (_safe, _cands))
        # 截断一律不越界：超过超胞安全截断的值自动压到安全截断（周期镜像会重复
        # 计数，力常数是错的）；fc3 截断留空（= 不截断，同样越界且大超胞内存爆）
        # 也取安全截断。hiphive 的 cutoff2/3 同理。
        _cap = float(int(_safe * 100)) / 100.0
        _keys = {"phono3py": ("FC3_CUTOFF",),
                 "pheasy": ("PHEASY_C3_CUTOFF",),
                 "hiphive": ("HIPHIVE_CUTOFF2", "HIPHIVE_CUTOFF3")}.get(engine, ())
        for _k in _keys:
            if _k in ("FC3_CUTOFF", "PHEASY_C3_CUTOFF", "HIPHIVE_CUTOFF3") \
                    and enable < 3:
                continue
            _raw = str(conf[_k] if conf[_k] is not None else "").strip()
            if _raw == "":
                _eff[_k] = _cap
                print("[..] %s 为空（不截断会越过超胞镜像，%d 原子超胞内存也吃不消）"
                      "→ 取超胞安全截断 %.2f Å" % (_k, _nsc, _cap), flush=True)
                continue
            try:
                _v = float(_raw)
            except ValueError:
                continue
            if _v > _safe + 1e-9:
                _eff[_k] = _cap
                print("[WARN] %s=%.2f > 超胞安全截断 %.2f Å（周期镜像会重复计数）"
                      "→ 已压到 %.2f Å" % (_k, _v, _safe, _cap), flush=True)
    _scan_req = str(conf["CUT3_SCAN"] or "off").strip().lower()
    if _scan_req not in ("off", "none", "false", "no", "0", "auto", "on", "true",
                         "1", "yes"):
        sys.exit("[ERROR] CUT3_SCAN 只允许 off 或 auto/on，收到 %r" % conf["CUT3_SCAN"])
    # The scan measures the |Phi^3| spread across refits, so it only means
    # something when fc3 is actually fitted (ENABLE_FC=3) and there are >= 2
    # candidates.
    _scan = (_scan_req in ("auto", "on", "true", "1", "yes")
             and len(_cands) >= 2 and enable >= 3)
    print("[..] CUT3_SCAN=%s（候选 %d 档）→ %s"
          % (conf["CUT3_SCAN"], len(_cands), "逐档扫描" if _scan else "单截断"))
    if _scan and engine != "pheasy":
        print("[..] 拟合器=%s：%s" % (engine, (
            "逐档 refit 写 cutoff_scan/cut3_<c>/（无 bootstrap，S2_kappa 只按 κ 平台选）"
            if engine == "phono3py" else
            "本步只做单档壳层报告；逐档 refit + bootstrap 仅 pheasy 有")))

    cfg = {
        "engine": engine,
        "dataset_dir": os.path.relpath(str(src), str(out)),
        "dataset_dir_abs": str(src),
        "dataset_signature": sig,
        "supercell": sc_str,
        "supercell_matrix": sc_rows,
        "primitive_matrix": primitive_matrix(src),
        "dim": dim,
        "enable_fc": enable,
        "subtract_equilibrium": bool(conf["SUBTRACT_EQUILIBRIUM"]),
        "eq_force_max": float(conf["EQ_FORCE_MAX"] or 0.0),
        "equilibrium_forces_npy": str(conf["EQUILIBRIUM_FORCES_NPY"] or ""),
        "coords": coords,
        # shell / third-order cutoff determination
        "cut3_candidates": [float(x) for x in _cands],
        "cut3_shells": [float(x) for x in _shells],
        "cut3_scan": bool(_scan),
        "cut3_bootstrap": int(conf["CUT3_BOOTSTRAP"] or 0),
        "cut3_stability_thr": float(conf["CUT3_STABILITY_THR"] or 0.3),
        "cut3_safe_cutoff": (None if _safe is None else float(_safe)),
        # phono3py
        "fc_calc": fc_calc,
        "fc3_cutoff": str(_eff.get("FC3_CUTOFF", conf["FC3_CUTOFF"]) or ""),
        # pheasy
        "pheasy_method": p_method,
        "pheasy_bin": p_bin,
        "pheasy_c2_cutoff": str(conf["PHEASY_C2_CUTOFF"] or ""),
        "pheasy_c3_cutoff": str(_eff.get("PHEASY_C3_CUTOFF",
                                         conf["PHEASY_C3_CUTOFF"]) or ""),
        "null_space_eps": float(conf["NULL_SPACE_EPS"] or 0.001),
        "pheasy_rasr": str(conf["PHEASY_RASR"] or ""),
        "pheasy_std": bool(conf["PHEASY_STD"]),
        "pheasy_lasso_twolevel": str(conf["PHEASY_LASSO_TWOLEVEL"] or ""),
        "pheasy_ngpu": p_ngpu,
        "pheasy_gpu_lasso_resident": str(conf["PHEASY_GPU_LASSO_RESIDENT"] or ""),
        "pheasy_cv_max_iter": str(conf["PHEASY_CV_MAX_ITER"] or ""),
        "pheasy_cv_tol": str(conf["PHEASY_CV_TOL"] or ""),
        "pheasy_tsqr_criterion": str(conf["PHEASY_TSQR_CRITERION"] or ""),
        "pheasy_tuning": str(conf["PHEASY_TUNING"] or "safe"),
        "pheasy_ols_ridge": str(conf["PHEASY_OLS_RIDGE"] or ""),
        "pheasy_ols_maxiter": str(conf["PHEASY_OLS_MAXITER"] or ""),
        "pheasy_cv": int(conf["PHEASY_CV"] or 5),
        "pheasy_nmu": p_nmu,
        "pheasy_mu_min": int(conf["PHEASY_MU_MIN"]),
        "pheasy_mu_max": p_mu_max,
        "pheasy_max_iter": p_max_iter,
        "pheasy_tol": p_tol,
        "pheasy_seed": int(conf["PHEASY_SEED"]),
        # hiphive
        "hiphive_cutoff2": float(_eff.get("HIPHIVE_CUTOFF2", conf["HIPHIVE_CUTOFF2"])),
        "hiphive_cutoff3": float(_eff.get("HIPHIVE_CUTOFF3", conf["HIPHIVE_CUTOFF3"])),
        "hiphive_fit_method": h_method,
        "hiphive_alpha": float(conf["HIPHIVE_ALPHA"]),
        "hiphive_symprec": float(conf["HIPHIVE_SYMPREC"]),
        "hiphive_enforce_asr": bool(conf["HIPHIVE_ENFORCE_ASR"]),
        "hiphive_n_configs": int(conf["HIPHIVE_N_CONFIGS"] or 0),
        "hiphive_cell": str(conf["HIPHIVE_CELL"] or ""),
        # output / gate
        "export_shengbte": bool(conf["EXPORT_SHENGBTE"]),
        "export_kl_bundle": bool(conf["EXPORT_KL_BUNDLE"]),
        "fc3_load_gb_limit": float(conf["FC3_LOAD_GB_LIMIT"] or 8.0),
        "band_points": int(conf["BAND_POINTS"]),
        "imag_thr": float(conf["IMAG_THR"]),
        "fit_rmse_frames": int(conf["FIT_RMSE_FRAMES"] or 0),
    }
    (out / "fit_config.json").write_text(
        json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")

    # The driver runs on the compute node: gen_need only places it next to this
    # gen script, so copy it into the step directory explicitly.
    here = Path(__file__).resolve().parent
    if not (here / "fc_fit_driver.py").is_file():
        sys.exit("[ERROR] fc_fit_driver.py missing -- is it listed in gen_need?")
    shutil.copyfile(str(here / "fc_fit_driver.py"), str(out / "fc_fit_driver.py"))

    kind = "submit_fcfit_%s" % ("p3py" if engine == "phono3py" else engine)
    if engine == "pheasy" and p_bin == "pheasy-gpu":
        kind = "submit_fcfit_pheasy_gpu"
    try:
        tpl = fc.resolve_submit(here, kind)
    except SystemExit:
        if kind == "submit_fcfit_pheasy_gpu":
            sys.exit("[ERROR] PHEASY_BIN=pheasy-gpu needs the GPU template "
                     "submit_fcfit_pheasy_gpu.tpl.\n"
                     "        If this cluster has no GPU nodes, drop a "
                     "cluster-specific copy into setting/<hpc>/templates/.")
        raise
    subs = {
        "JOBNAME": fc.new_jobname(cwd, job_label),
        "CONDA_SH": str(conf["CONDA_SH"] or ""),
        "CONDA_ENV": str(conf["CONDA_ENV"] or ""),
        "ENGINE": engine,
    }
    fc.write_submit(tpl, out / "submit.sh", subs)
    sub = dict(conf.submit or {})
    if p_ngpu and p_bin == "pheasy-gpu":
        # Keep --gres=gpu:N in step with the pheasy device count (an explicit
        # [submit] gres still wins).
        sub.setdefault("gres", "gpu:%d" % int(p_ngpu))
    stepconf.apply_submit(out / "submit.sh", sub)

    print("[..] engine=%s enable_fc=%d dim=%s supercell=%s"
          % (engine, enable, dim.upper(), sc_str or "?"))
    print("[DONE] %s: fit_config.json + submit.sh + fc_fit_driver.py ready.  "
          "tf submits the job; the compute node writes fc2.hdf5 (+ fc3.hdf5), "
          "shengbte/ and phonon_summary.json, which the built-in phonon judge "
          "reads (stable | imaginary frequency | tool error)."
          % (OUTDIR if out == cwd / OUTDIR else os.path.relpath(str(out), str(cwd))))
    return cfg


if __name__ == "__main__":
    main()
