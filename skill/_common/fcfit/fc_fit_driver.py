#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fc_fit_driver.py -- compute-node driver for the fit-fc-thermal skill (S1_fit).

This file is the force-constant fitting branch extracted from
kl-dft-cpu/kl_fc_backends.py, generalised into a standalone skill and extended
with a third fitting engine (hiphive).  submit.sh calls its sub-commands in
order; every sub-command runs inside the step directory (step1_fit/):

  prep    Normalise the input dataset into a canonical set of files:
            POSCAR, SPOSCAR                       unit cell / supercell
            dataset_disps.npy, dataset_forces.npy (n+1, Nsc, 3) Cartesian
                                                  displacements / forces, the
                                                  trailing frame being the zero
                                                  equilibrium reference
            disp_matrix.pkl, force_matrix.pkl     (n, Nsc, 3) equilibrium-
                                                  subtracted, for pheasy
            phono3py_disp.yaml / phono3py_params.yaml  copied through when the
                                                  source has them (phono3py engine)
            fc_dataset.json                       provenance / dataset audit
          Supported input signatures, tried in this order inside the dataset dir:
            1. phono3py_params.yaml          (forces embedded)
            2. phono3py_disp.yaml + FORCES_FC3
            3. phono3py_disp.yaml            (forces embedded)
            4. dataset_disps.npy + dataset_forces.npy (+ POSCAR/SPOSCAR)
            5. disp_matrix.pkl + force_matrix.pkl (+ POSCAR/SPOSCAR)
            6. disp-*/vasprun.xml            (+ phono3py_disp.yaml or SPOSCAR)
          Sources 1-3 and 6 need phono3py; sources 4-5 need only phonopy.

  fit     Dispatch on FIT_ENGINE and write fc2.hdf5 (+ fc3.hdf5):
            phono3py  phono3py + symfc/alm least squares
            pheasy    pheasy CLI, four steps (cluster space -> symmetry
                      constraints -> sensing matrix -> fit)
            hiphive   hiphive cluster space + linear regression, optional
                      rotational-sum-rule projection

  post    Export (optional) + imaginary-frequency gate:
            shengbte/FORCE_CONSTANTS_2ND, _3RD, POSCAR   (EXPORT_SHENGBTE=true)
            band-dft-cpu.yaml                             (best effort)
            fc_fit_summary.json      <- S1 marker: "FIT_DONE": true
          Exit status: a genuine imaginary frequency exits 0 (the marker is
          simply not satisfied, so downstream steps are held back); a tool
          error (mesh could not be computed at all) exits non-zero so tf shows
          the step as error.

Engine-independent fit quality: FIT_RMSE_FRAMES > 0 evaluates the fitted force
constants with hiphive's ForceConstantCalculator on that many training frames
and records RMSE / relative RMSE in fc_fit_summary.json.
"""
import glob
import json
import math
import os
import pickle
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import imag_policy       # 近 Γ 声学支虚频判据唯一真源（_common/imag_policy.py）
import phonon_stability  # 虚频/ZA 编排唯一真源（_common/phonon_stability.py，import za_2d）

OUTDIR = "step1_fit"
SB_SUB = "shengbte"
P3_SUB = "phono3py"

# Recognised dataset signatures, most specific first.
SIGNATURES = (
    ("phono3py_params", ("phono3py_params.yaml",)),
    ("phono3py_disp_forces", ("phono3py_disp.yaml", "FORCES_FC3")),
    ("phono3py_disp", ("phono3py_disp.yaml",)),
    ("arrays", ("dataset_disps.npy", "dataset_forces.npy")),
    ("pkl", ("disp_matrix.pkl", "force_matrix.pkl")),
    ("vasprun", ("phono3py_disp.yaml",)),   # validated further by disp-* below
)


# ==========================================================================
# Small helpers
# ==========================================================================
def load_cfg(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def run(cmd, **kw):
    """Run a command in the foreground; exit on a non-zero return code."""
    print("[cmd] %s" % (cmd if isinstance(cmd, str) else " ".join(cmd)), flush=True)
    r = subprocess.run(cmd, shell=isinstance(cmd, str), **kw)
    if r.returncode != 0:
        sys.exit("[ERROR] command failed (rc=%d): %s" % (r.returncode, cmd))
    return r


def _die_with_parent():
    """preexec_fn: the child gets SIGTERM when this driver dies (Linux only).

    Under SLURM scancel reaps the whole job cgroup, but a driver run by hand or
    by a local orchestrator that times it out (subprocess.run(timeout=...) sends
    SIGKILL to the driver only) used to leave the pheasy grandchild running as an
    orphan at 100 % CPU.  Best effort: silently a no-op where prctl is missing.
    """
    try:
        import ctypes
        import signal
        ctypes.CDLL("libc.so.6", use_errno=True).prctl(1, signal.SIGTERM)  # PR_SET_PDEATHSIG
    except Exception:
        pass


def _ensure_env_bin_on_path(binary):
    """Find a console script next to the running interpreter when PATH lacks it.

    submit.sh activates the conda env, so PATH is right on a compute node; a
    driver started as <env>/bin/python without activation (local runs, agents)
    has the env's python but not its bin/ on PATH, and pheasy was "not found".
    """
    if shutil.which(binary):
        return
    env_bin = str(Path(sys.executable).parent)   # no resolve(): keep venv symlinks
    if (Path(env_bin) / binary).is_file():
        os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")
        print("[..] %s found in %s (not on PATH) -- prepended it to PATH"
              % (binary, env_bin), flush=True)


def _rel_or_abs(p, base):
    p = Path(p).resolve()
    try:
        return str(p.relative_to(Path(base).resolve()))
    except ValueError:
        return str(p)


def _is_none(v):
    return v in (None, "", "None", "none", "null", "False", "false")


def _as_float_or_none(v):
    return None if _is_none(v) else float(v)


def _truthy(v, default=False):
    if v is None:
        return default
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("true", "1", "yes", "y", "on")


# ==========================================================================
# Reading a phonopy/phono3py-style YAML without importing phono3py
# ==========================================================================
def _atoms_from_phonopy_block(block):
    """Build an ase.Atoms from a phonopy YAML 'cell' block (lattice + points)."""
    from ase import Atoms
    lat = np.asarray(block["lattice"], float)
    syms = [p["symbol"] for p in block["points"]]
    frac = np.asarray([p["coordinates"] for p in block["points"]], float)
    return Atoms(symbols=syms, cell=lat, scaled_positions=frac, pbc=True)


def _load_yaml(path):
    """Load a phono3py YAML file. Returns None when neither PyYAML nor the file
    is usable."""
    try:
        import yaml
    except ImportError:
        return None
    try:
        return yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None


def _dataset_from_yaml(y):
    """Extract the displacement dataset from a parsed phono3py YAML.

    Returns (kind, disps, forces) where kind is 'random' (type-2, an
    (n, Nsc, 3) pair of arrays) or 'findiff' (type-1, first_atoms)."""
    if not isinstance(y, dict):
        return None, None, None
    for key in ("dataset", "displacement_dataset"):
        ds = y.get(key)
        if not isinstance(ds, dict):
            continue
        if ds.get("displacements") is not None:
            disps = np.asarray(ds["displacements"], float)
            forces = ds.get("forces")
            forces = None if forces is None else np.asarray(forces, float)
            return "random", disps, forces
        if ds.get("first_atoms"):
            return "findiff", None, None
    return None, None, None


# ==========================================================================
# prep: normalise the dataset
# ==========================================================================
def _resolve_dataset_dir(cfg, cwd):
    """Return (path, signature) for the dataset directory, or exit with help."""
    raw = str(cfg.get("dataset_dir") or "").strip()
    if not raw:
        sys.exit("[ERROR] fit_config.json has no dataset_dir -- rerun the gen step")
    cands = [Path(raw)] if raw else []
    for c in cands:
        if c.is_dir():
            kind = _detect_signature(c)
            if kind:
                return c, kind
    sys.exit("[ERROR] dataset dir %s does not exist or holds no recognised "
             "displacement dataset (looked for %s)"
             % (raw, ", ".join(s[0] for s in SIGNATURES)))


def _detect_signature(d):
    """Classify a directory by the dataset files it contains (None if unknown)."""
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


def _write_poscar(path, atoms, comment="generated by fit-fc-thermal"):
    """Minimal VASP POSCAR writer (keeps us independent of phonopy for this)."""
    lat = np.asarray(atoms.cell, float)
    fr = np.asarray(atoms.get_scaled_positions(wrap=True), float)
    syms = list(atoms.get_chemical_symbols())
    order, counts = [], []
    for s in syms:
        if s not in order:
            order.append(s)
            counts.append(syms.count(s))
    lines = [comment, "   1.0"]
    for v in lat:
        lines.append("  %20.16f %20.16f %20.16f" % tuple(v))
    lines.append("  " + "  ".join(order))
    lines.append("  " + "  ".join(str(c) for c in counts))
    lines.append("Direct")
    for s in order:
        for i, si in enumerate(syms):
            if si != s:
                continue
            lines.append("  %20.16f %20.16f %20.16f" % tuple(fr[i]))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _read_poscar(path):
    from ase.io import read as ase_read
    return ase_read(str(path), format="vasp")


def _pheasy_supercell_order(uc, sc, dim):
    """Permutation from this dataset's supercell order to pheasy's own.

    pheasy does not read SPOSCAR: it rebuilds the supercell from POSCAR and
    SUPERCELL as per-primitive-atom blocks (all periodic images of primitive
    atom 0 first, then atom 1, ...), each block in the ndindex(dim[::-1])
    translation order.  phonopy/phono3py happen to agree with that layout, but
    ASE Atoms.repeat() interleaves the images and does not -- every atom except
    the first is then scrambled, the fit silently lands near 50% residual with
    a force correlation around 0.7, and pheasy still exits 0.

    Returns the index array 'perm' such that row i of the pheasy-ordered array
    is row perm[i] of the dataset array, or None when the two supercells cannot
    be matched (different shape, rotated cell, missing atoms).
    """
    dim = [int(x) for x in dim]
    if len(dim) != 3 or len(uc) * int(np.prod(dim)) != len(sc):
        return None
    # Cartesian positions of pheasy's supercell: (trans + pos_uc) . cell_uc
    cell = np.asarray(uc.cell, float)
    sc_cell = np.asarray(sc.cell, float)
    try:
        inv = np.linalg.inv(sc_cell)
    except np.linalg.LinAlgError:
        return None
    frac_sc = np.asarray(sc.get_scaled_positions(wrap=True), float) % 1.0
    trans = [list(t)[::-1] for t in np.ndindex(*dim[::-1])]
    perm = []
    for pos in np.asarray(uc.get_scaled_positions(wrap=True), float):
        for t in trans:
            cart = (np.asarray(t, float) + pos) @ cell
            f = (cart @ inv) % 1.0
            d = frac_sc - f
            d -= np.round(d)
            j = int(np.argmin(np.abs(d).sum(axis=1)))
            if np.abs(d[j]).max() > 1e-5:
                return None
            perm.append(j)
    if sorted(perm) != list(range(len(sc))):
        return None
    return np.asarray(perm, int)


def _recorded_pheasy_perm(out):
    """The permutation prep froze in fc_dataset.json, or None."""
    try:
        rec = json.loads((Path(out) / "fc_dataset.json").read_text(
            encoding="utf-8")).get("pheasy_atom_order")
    except Exception:
        return None
    if not isinstance(rec, dict) or not rec.get("permutation"):
        return None
    perm = np.asarray(rec["permutation"], int)
    return None if np.array_equal(perm, np.arange(len(perm))) else perm


def apply_pheasy_order(cfg, out):
    """Put disp_matrix.pkl / force_matrix.pkl into pheasy's atom order.

    Only the two pickles are touched -- they are the pheasy-specific inputs --
    so hiphive (which aligns the supercell internally) and phono3py (which uses
    the YAML/SPOSCAR that prep wrote, in the dataset order) are unaffected.
    Idempotent: the decision is recorded in .pheasy_order.json.
    """
    out = Path(out)
    rec = out / ".pheasy_order.json"
    if rec.is_file():
        return json.loads(rec.read_text(encoding="utf-8")).get("applied", False)
    dp, fp = out / "disp_matrix.pkl", out / "force_matrix.pkl"
    if not (dp.is_file() and fp.is_file()):
        return False
    # The permutation was frozen by prep, while SPOSCAR still described the
    # dataset: pheasy overwrites SPOSCAR with its own atom order, so reading
    # it back here would compare pheasy against pheasy and always say "fine".
    rec_order = None
    try:
        rec_order = json.loads((out / "fc_dataset.json").read_text(
            encoding="utf-8")).get("pheasy_atom_order")
    except Exception:
        rec_order = None
    if isinstance(rec_order, dict) and rec_order.get("permutation"):
        perm = np.asarray(rec_order["permutation"], int)
    else:
        # older prep output (or a hand-built step dir): fall back to the files
        scm, _pm = _cells_cfg(cfg, out)
        if scm is None or np.shape(scm) != (3, 3):
            scm = np.asarray(_diag_supercell_matrix(out), float)
        d = np.diag(scm)
        if np.abs(scm - np.diag(d)).max() > 1e-8:
            print("[WARN] the dataset supercell matrix is not diagonal; pheasy's "
                  "--dim cannot express it, so the atom order cannot be checked",
                  flush=True)
            perm = None
        else:
            perm = _pheasy_supercell_order(_read_poscar(out / "POSCAR"),
                                           _read_poscar(out / "SPOSCAR"),
                                           [int(round(x)) for x in d])
    applied = False
    if perm is None:
        print("[WARN] could not verify the supercell atom order against the one "
              "pheasy builds from POSCAR + SUPERCELL -- if the fit residual is "
              "near 50%% with a force correlation around 0.7, the dataset was "
              "written in ASE repeat() order instead", flush=True)
    elif not np.array_equal(perm, np.arange(len(perm))):
        with open(dp, "rb") as fh:
            disps = pickle.load(fh)
        with open(fp, "rb") as fh:
            forces = pickle.load(fh)
        with open(dp, "wb") as fh:
            pickle.dump(np.asarray(disps, float)[:, perm, :], fh)
        with open(fp, "wb") as fh:
            pickle.dump(np.asarray(forces, float)[:, perm, :], fh)
        applied = True
        print("[WARN] the dataset supercell is not in the atom order pheasy "
              "rebuilds from POSCAR + SUPERCELL, which would have scrambled "
              "every atom but the first.  disp_matrix.pkl / force_matrix.pkl "
              "were permuted into pheasy's order (npy files, SPOSCAR and every "
              "other engine are untouched).", flush=True)
    rec.write_text(json.dumps({"applied": applied,
                               "permutation": perm.tolist() if perm is not None else None},
                              indent=2), encoding="utf-8", newline="\n")
    return applied

def apply_pheasy_fc_order(cfg, out):
    """Permute pheasy's fitted force constants back to the dataset atom order.

    pheasy rebuilds the supercell from POSCAR + SUPERCELL in its own atom order
    (recorded as `pheasy_atom_order.permutation`) and writes fc2/fc3/fc4.hdf5 in
    that order.  Every downstream consumer that rebuilds the supercell from
    POSCAR + SUPERCELL (phono3py-loaded YAML, ShengBTE, hiphive, the phonon
    gate) expects the *dataset* order, so leaving pheasy's order here silently
    scrambles the force constants: the training residual stays ~1% and the
    phonon gate can even pass, while kappa is off by tens of percent
    (Mn2In2Se5 3x3x1: kappa_xx 1.503 vs 1.150).  This rewrites the hdf5 files
    in place and restores SPOSCAR to the dataset order saved by cmd_fit_pheasy.
    Idempotent: each file's post-permutation md5 is recorded in
    .pheasy_fc_reordered.json, so a fresh fit is re-permuted while a skipped fit
    is not permuted twice.
    """
    import hashlib
    import h5py
    perm = _recorded_pheasy_perm(out)
    if perm is None:
        return False
    inv = np.argsort(perm)
    sig = hashlib.md5(np.asarray(inv, dtype="<i8").tobytes()).hexdigest()[:12]
    marker = out / ".pheasy_fc_reordered.json"
    state = {}
    if marker.is_file():
        try:
            state = json.loads(marker.read_text(encoding="utf-8"))
        except Exception:
            state = {}

    def _md5(path):
        h = hashlib.md5()
        with h5py.File(str(path), "r") as fh:
            for k in sorted(fh.keys()):
                if k == "version":
                    continue
                h.update(np.ascontiguousarray(fh[k][()]).tobytes())
        return h.hexdigest()[:12]

    changed = False
    for name, keys in (("fc2.hdf5", ("fc2", "force_constants")),
                       ("fc3.hdf5", ("fc3", "force_constants_third")),
                       ("fc4.hdf5", ("fc4", "force_constants_fourth"))):
        p = out / name
        if not p.is_file():
            continue
        cur = _md5(p)
        old = state.get(name) or {}
        if old.get("post_md5") == cur and old.get("sigma_md5") == sig:
            continue
        with h5py.File(str(p), "r+") as fh:
            key = next((k for k in keys if k in fh), None)
            if key is None:
                continue
            a = np.asarray(fh[key])
            if a.ndim == 0 or a.shape[0] != len(inv):
                continue
            nax = a.ndim // 2
            o = np.take(a, inv, axis=0)
            for ax in range(1, nax):
                o = np.take(o, inv, axis=ax)
            fh[key][...] = o
        state[name] = {"sigma_md5": sig, "post_md5": _md5(p)}
        changed = True
        print("[fc-order] %s -> dataset order (%d/%d atoms moved)"
              % (name, int((inv != np.arange(len(inv))).sum()), len(inv)), flush=True)
    marker.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")
    # pheasy overwrote SPOSCAR with its own order; put the dataset's back.
    ds = out / "SPOSCAR.dataset"
    if ds.is_file():
        shutil.copyfile(str(ds), str(out / "SPOSCAR"))
    return changed


def _phonopy_unitcell(path):
    """Unit cell in the form phonopy 2.x accepts.

    phonopy wants its own Atoms (or an explicit (lattice, scaled_positions,
    numbers) tuple); handing it an ASE Atoms fails later with the misleading
    "'Atoms' object has no attribute 'scaled_positions'"."""
    from phonopy.structure.atoms import PhonopyAtoms
    a = _read_poscar(path)
    return PhonopyAtoms(symbols=a.get_chemical_symbols(),
                        cell=np.asarray(a.cell, float),
                        scaled_positions=a.get_scaled_positions(wrap=True))


def _equilibrium_from_vasprun(disp_dir, dirstr):
    """Return (files, missing, equilibrium_file) for disp-*/vasprun.xml."""
    base = Path(disp_dir)
    subs = sorted(glob.glob(str(base / "disp-*")),
                  key=lambda p: int(re.search(r"disp-(\d+)", p).group(1)))
    files, missing, eq = [], [], None
    for d in subs:
        num = re.search(r"disp-(\d+)", d).group(1)
        vr = Path(d) / "vasprun.xml"
        ok = vr.is_file() and vr.stat().st_size
        rel = "%s/disp-%s/vasprun.xml" % (dirstr, num)
        if int(num) == 0:
            eq = rel if ok else None
        elif ok:
            files.append(rel)
        else:
            missing.append(num)
    return files, missing, eq


def _stage_yaml_cells(src, out):
    """Copy POSCAR/SPOSCAR and the phonopy YAMLs from a yaml-based source.

    Returns a dict with the cell matrices when they could be recovered."""
    info = {}
    for f in ("POSCAR", "SPOSCAR", "phono3py_disp.yaml", "phono3py_params.yaml",
              "FORCES_FC3", "BORN"):
        if (src / f).is_file():
            shutil.copyfile(str(src / f), str(out / f))
    y = None
    for f in ("phono3py_params.yaml", "phono3py_disp.yaml"):
        if (out / f).is_file():
            y = _load_yaml(out / f)
            if y:
                break
    if y:
        for key, mat in (("supercell_matrix", "supercell_matrix"),
                         ("primitive_matrix", "primitive_matrix")):
            if y.get(key) is not None:
                info[mat] = [[float(x) for x in row] for row in np.asarray(y[key], float)]
        # Synthesise the cells from the YAML only when the source did not ship
        # them as POSCAR/SPOSCAR (the phonopy YAMLs call them unit_cell /
        # supercell and store fractional coordinates).
        for blk, dst in (("unit_cell", "POSCAR"), ("supercell", "SPOSCAR")):
            if (out / dst).is_file() or not isinstance(y.get(blk), dict):
                continue
            try:
                _write_poscar(out / dst, _atoms_from_phonopy_block(y[blk]))
                print("[..] %s written from the YAML block %s" % (dst, blk), flush=True)
            except Exception as e:
                print("[WARN] could not write %s from %s: %s" % (dst, blk, e))
    return info


def _min_image_cart(frac_delta, cell):
    u = np.asarray(frac_delta, float)
    u = u - np.round(u)
    return u @ np.asarray(cell, float)


def _detect_coords_mode(d, cell, frac_ref):
    # Decide whether d holds absolute fractional coordinates or Cartesian
    # displacements, with no user input.  Absolute fractional coordinates
    # live in [0, 1) and the frame closest to the ideal SPOSCAR sits within
    # ~1 A of it; Cartesian displacements are centred on zero, so subtracting
    # the ideal POSITIONS leaves a several-Angstrom residual and ~half the
    # values are negative.  Returns (mode, closest_to_ideal_A, frac_like).
    d = np.asarray(d, float)
    frac_ref = np.asarray(frac_ref, float)
    devs = []
    for i in range(len(d)):
        cart = _min_image_cart(d[i] - frac_ref, cell)
        devs.append(float(np.sqrt((cart ** 2).sum(axis=1).mean())))
    min_dev = min(devs) if devs else 0.0
    flat = d.reshape(-1)
    frac_like = float(np.mean((flat >= -1e-6) & (flat < 1.0)))
    mode = 'fractional' if (min_dev < 1.0 or frac_like > 0.95) else 'cartesian'
    return mode, min_dev, frac_like


def _normalize_displacement_frames(d, f, cell, frac_ref, mode):
    # Normalise a raw (m, N, 3) frame set into the fit-fc-thermal canonical form.
    # The equilibrium frame is wherever it is (first, middle, last or absent):
    # the frame closest to the ideal supercell is located by its RMS
    # displacement.  If it really is the ideal structure it becomes the
    # reference and is dropped from training; otherwise every frame is kept
    # and the closest frame's forces are the equilibrium approximation when
    # it is close enough.  Returns (disps, forces, eq_forces, info).
    d = np.asarray(d, float)
    f = np.asarray(f, float)
    frac_ref = np.asarray(frac_ref, float)
    if mode == 'fractional':
        disp_all = np.stack([_min_image_cart(d[i] - frac_ref, cell) for i in range(len(d))])
    else:
        disp_all = d.copy()
    rms = np.sqrt((disp_all ** 2).sum(axis=2).mean(axis=1))
    info = {'coords_detected': mode}
    if not len(rms):
        return disp_all, f, None, info
    iref = int(np.argmin(rms))
    info['reference_frame_index'] = iref
    if float(rms[iref]) < 1e-4:
        keep = [i for i in range(len(d)) if i != iref]
        info['equilibrium_source'] = 'ideal frame #%d (%s)' % (iref, mode)
        return disp_all[keep], f[keep], f[iref].copy(), info
    info['equilibrium_source'] = None
    eq_forces = None
    if float(rms[iref]) < 0.25:
        eq_forces = f[iref].copy()
        info['equilibrium_source'] = 'nearest frame #%d to the ideal (rms %.3f A)' % (iref, float(rms[iref]))
    return disp_all, f, eq_forces, info


def cmd_prep(cfg, out):
    out = Path(out)
    src, kind = _resolve_dataset_dir(cfg, out)
    print("[..] dataset source: %s (%s)" % (src, kind), flush=True)

    info = {"dataset_dir": str(src), "signature": kind}
    disps = forces = None
    eq_forces = None

    if kind in ("phono3py_params", "phono3py_disp", "phono3py_disp_forces"):
        info.update(_stage_yaml_cells(src, out))
        y = None
        for f in ("phono3py_params.yaml", "phono3py_disp.yaml"):
            if (out / f).is_file():
                y = _load_yaml(out / f)
                if y:
                    break
        dkind, disps, forces = _dataset_from_yaml(y)
        info["dataset_type"] = dkind
        if disps is None:
            # Either a finite-displacement dataset (phono3py rebuilds it from
            # FORCES_FC3) or forces missing entirely -- the phono3py engine can
            # still consume the YAML directly, so this is not fatal for it.
            print("[WARN] no type-2 (random) displacement array in the YAML; "
                  "FIT_ENGINE=pheasy / hiphive need one", flush=True)
    elif kind == "vasprun":
        info.update(_stage_yaml_cells(src, out))
        dirstr = _rel_or_abs(src, out)
        files, missing, eqf = _equilibrium_from_vasprun(src, dirstr)
        if not files:
            sys.exit("[ERROR] no disp-*/vasprun.xml under %s" % src)
        if eqf is None:
            sys.exit("[ERROR] missing equilibrium frame disp-00000/vasprun.xml "
                     "under %s -- it anchors the residual-force subtraction" % src)
        from phonopy.interface.vasp import parse_set_of_forces
        sc = _read_poscar(out / "SPOSCAR") if (out / "SPOSCAR").is_file() else None
        if sc is None:
            sys.exit("[ERROR] no SPOSCAR next to the disp-* directories")
        nsc = len(sc)
        rd = parse_set_of_forces(nsc, [str(Path(f)) for f in files + [eqf]],
                                 verbose=False)
        allf = np.asarray(rd["forces"], float)
        forces, eq_forces = allf[:-1], allf[-1]
        # Displacements come from the phono3py YAML when present.
        y = _load_yaml(out / "phono3py_disp.yaml") if (out / "phono3py_disp.yaml").is_file() else None
        dkind, d, _ = _dataset_from_yaml(y)
        if d is None:
            sys.exit("[ERROR] the vasprun source needs phono3py_disp.yaml with a "
                     "type-2 displacement array (finite-displacement datasets are "
                     "rebuilt by the phono3py engine, not here)")
        idx = [int(re.search(r"disp-(\d+)", f).group(1)) - 1 for f in files]
        disps = np.asarray(d, float)[idx]
        if missing:
            print("[WARN] %d frame(s) missing (disp-%s); random-displacement "
                  "regression tolerates this" % (len(missing), ",".join(missing[:8])))
        info["missing_frames"] = missing
        info["equilibrium_source"] = "disp-00000/vasprun.xml"
    elif kind == "arrays":
        y = None
        for f in ("phono3py_params.yaml", "phono3py_disp.yaml"):
            if (src / f).is_file():
                shutil.copyfile(str(src / f), str(out / f))
        for f in ("POSCAR", "SPOSCAR", "BORN", "FORCES_FC3"):
            if (src / f).is_file():
                shutil.copyfile(str(src / f), str(out / f))
        d = np.load(str(src / "dataset_disps.npy"))
        f_ = np.load(str(src / "dataset_forces.npy"))
        if d.ndim != 3 or f_.ndim != 3:
            sys.exit("[ERROR] dataset_disps/forces.npy must be (n, natom, 3)")
        if len(d) != len(f_):
            sys.exit("[ERROR] dataset_disps.npy has %d frames but "
                     "dataset_forces.npy has %d" % (len(d), len(f_)))
        coords = str(cfg.get("coords") or "auto").strip().lower()
        if coords not in ("auto", "cartesian", "fractional"):
            sys.exit("[ERROR] COORDS must be auto|cartesian|fractional, got %r"
                     % coords)
        sp_path = out / "SPOSCAR"
        if sp_path.is_file():
            sc0 = _read_poscar(sp_path)
            cell0 = np.asarray(sc0.cell, float)
            frac_ref0 = np.asarray(sc0.get_scaled_positions(wrap=True), float)
        else:
            cell0 = frac_ref0 = None
        if coords == "auto" and cell0 is None:
            coords = "cartesian"
            print("[WARN] COORDS=auto needs SPOSCAR to compare against; falling "
                  "back to Cartesian displacements.", flush=True)
        if coords == "fractional" and cell0 is None:
            sys.exit("[ERROR] COORDS=fractional needs SPOSCAR next to the dataset "
                     "to convert fractional coordinates -- none found in %s." % src)
        if coords == "auto":
            coords, _dev0, _fl0 = _detect_coords_mode(d, cell0, frac_ref0)
            print("[..] COORDS=auto -> %s (closest frame to the ideal supercell "
                  "%.3f A, %.0f%% of values in [0,1))"
                  % (coords, _dev0, 100.0 * _fl0), flush=True)
        info["coords_mode"] = coords
        if coords == "fractional":
            # dataset_disps.npy holds ABSOLUTE fractional coordinates: measure
            # every frame against the ideal SPOSCAR (minimum-image wrap), convert
            # to Cartesian, and use the frame closest to the ideal as the
            # equilibrium reference wherever it sits in the file.
            disps, forces, eq_forces, _ninfo = _normalize_displacement_frames(
                d, f_, cell0, frac_ref0, "fractional")
        else:
            # dataset_disps.npy already holds Cartesian displacements; a frame
            # of ~zero displacement (anywhere) is the equilibrium reference.
            disps, forces, eq_forces, _ninfo = _normalize_displacement_frames(
                d, f_, None, None, "cartesian")
        info.update(_ninfo)
        info["dataset_type"] = "random"
    else:   # pkl
        for f in ("POSCAR", "SPOSCAR", "BORN", "phono3py_params.yaml",
                  "phono3py_disp.yaml"):
            if (src / f).is_file():
                shutil.copyfile(str(src / f), str(out / f))
        disps = np.asarray(pickle.loads((src / "disp_matrix.pkl").read_bytes()), float)
        forces = np.asarray(pickle.loads((src / "force_matrix.pkl").read_bytes()), float)
        info["dataset_type"] = "random"

    if disps is None or forces is None:
        n_frames = 0
    else:
        disps = np.asarray(disps, float)
        forces = np.asarray(forces, float)
        n_frames = len(disps)
        if disps.shape != forces.shape:
            sys.exit("[ERROR] displacement/force shapes differ: %s vs %s"
                     % (disps.shape, forces.shape))

    # ---- structure: require POSCAR + SPOSCAR for every engine ----
    if not (out / "POSCAR").is_file() or not (out / "SPOSCAR").is_file():
        sys.exit("[ERROR] could not obtain POSCAR/SPOSCAR from %s (signature %s). "
                 "The fit needs both the unit cell and the supercell." % (src, kind))
    sc = _read_poscar(out / "SPOSCAR")
    nsc = len(sc)
    if n_frames and disps.shape[1] != nsc:
        sys.exit("[ERROR] dataset has %d atoms per frame but SPOSCAR has %d -- the "
                 "displacement set and the supercell do not belong together"
                 % (disps.shape[1], nsc))

    # ---- freeze the dataset's own supercell atom order -------------------
    # pheasy rebuilds the supercell from POSCAR + SUPERCELL in its own atom
    # order and *overwrites SPOSCAR* while fitting, so the order has to be
    # captured now, while SPOSCAR still describes the dataset.  Everything that
    # needs to know whether the dataset agrees with pheasy's convention reads
    # this record instead of the file.
    try:
        _scm_guess = info.get("supercell_matrix") or _diag_supercell_matrix(out)
        _d = np.diag(np.asarray(_scm_guess, float))
        _uc, _sp = _read_poscar(out / "POSCAR"), _read_poscar(out / "SPOSCAR")
        if np.abs(np.asarray(_scm_guess, float)
                  - np.diag(_d)).max() <= 1e-8:
            _perm = _pheasy_supercell_order(_uc, _sp, [int(round(x)) for x in _d])
        else:
            _perm = None
        info["pheasy_atom_order"] = (
            None if _perm is None else {
                "is_pheasy_order": bool(np.array_equal(_perm, np.arange(len(_perm)))),
                "permutation": [int(x) for x in _perm]})
        if _perm is not None and not np.array_equal(_perm, np.arange(len(_perm))):
            print("[WARN] the dataset supercell is not in the atom order pheasy "
                  "rebuilds from POSCAR + SUPERCELL (ASE repeat() interleaves the "
                  "images, pheasy blocks them per primitive atom).  The pheasy "
                  "engine will permute its input into pheasy's order; hiphive and "
                  "phono3py are unaffected.", flush=True)
    except Exception as e:
        info["pheasy_atom_order"] = None
        print("[WARN] could not determine the supercell atom order: %s" % e,
              flush=True)

    # ---- equilibrium reference -------------------------------------------
    # Priority: EQUILIBRIUM_FORCES_NPY, the reference frame the dataset itself
    # carried (disp-00000 or a trailing zero-displacement frame), then a
    # reference frame shipped next to the dataset.  A zero-displacement frame is
    # never a training sample: it is dropped and only used as the reference.
    subtract = _truthy(cfg.get("subtract_equilibrium"), True)
    eq_applied = False
    if subtract and eq_forces is None and n_frames:
        eq_forces = _read_equilibrium_npy(cfg, out, nsc)
        if eq_forces is not None:
            info["equilibrium_source"] = "EQUILIBRIUM_FORCES_NPY"

    zero_idx = [i for i in range(n_frames) if np.allclose(disps[i], 0.0)] if n_frames else []
    if zero_idx:
        if eq_forces is None:
            eq_forces = np.asarray(forces, float)[zero_idx[0]].copy()
            info["equilibrium_source"] = "zero-displacement frame in the dataset"
        keep = [i for i in range(n_frames) if i not in zero_idx]
        disps = np.asarray(disps, float)[keep]
        forces = np.asarray(forces, float)[keep]
        n_frames = len(keep)
        info["dropped_reference_frames"] = zero_idx
        print("[..] %d zero-displacement frame(s) kept as the equilibrium "
              "reference and excluded from training" % len(zero_idx), flush=True)

    if subtract and eq_forces is None and n_frames:
        eq_forces = _equilibrium_from_npy_pair(src, nsc)
        if eq_forces is not None:
            info["equilibrium_source"] = "reference frame in dataset_disps/forces.npy"

    if subtract and eq_forces is not None and n_frames:
        forces = forces - np.asarray(eq_forces, float)[None]
        eq_applied = True
        info["equilibrium_max_force"] = float(
            np.linalg.norm(np.asarray(eq_forces, float), axis=1).max())
        print("[OK] equilibrium residual subtracted from %s (max|F_eq| = %.4f eV/A)"
              % (info.get("equilibrium_source", "reference"), info["equilibrium_max_force"]),
              flush=True)
        # 残余力闸：扣除只修正一阶项，参考结构离能量极小太远时 fc2 本身就偏了。
        _fmax = float(cfg.get("eq_force_max") or 0.0)
        if _fmax > 0 and info["equilibrium_max_force"] > _fmax:
            sys.exit("[ERROR] 平衡帧残余力 max|F_eq| = %.4f eV/A > EQ_FORCE_MAX = %.3f "
                     "—— 参考结构没弛豫好，先重新弛豫再生成数据集（确要硬跑：调大 "
                     "EQ_FORCE_MAX，0 = 关闭此闸）"
                     % (info["equilibrium_max_force"], _fmax))
    elif n_frames:
        print("[..] no equilibrium reference available -- forces used as given; "
              "set EQUILIBRIUM_FORCES_NPY if the reference cell is not force-free",
              flush=True)

    # ---- write the canonical dataset ----
    if n_frames:
        zeros = np.zeros((1, nsc, 3))
        np.save(out / "dataset_disps.npy", np.concatenate([disps, zeros], axis=0))
        np.save(out / "dataset_forces.npy",
                np.concatenate([forces, np.zeros((1, nsc, 3))], axis=0))
        # Always rewrite the pickles in the step directory: they are what
        # pheasy reads through --disp_file, and they must hold exactly the
        # equilibrium-subtracted frames used for every other engine.
        with open(out / "disp_matrix.pkl", "wb") as fh:
            pickle.dump(disps, fh)
        with open(out / "force_matrix.pkl", "wb") as fh:
            pickle.dump(forces, fh)
        rms = float(np.sqrt((disps ** 2).mean()))
        if not (1e-8 < rms < 1.0):
            sys.exit("[ERROR] displacement RMS %.3e A is not physical -- check the "
                     "dataset units (Cartesian A expected)" % rms)
    else:
        rms = None

    # ---- unit cell / supercell bookkeeping for the gate ----
    scm = info.get("supercell_matrix")
    if scm is None:
        scm = _diag_supercell_matrix(out)
        info["supercell_matrix"] = scm
        info["supercell_matrix_guess"] = True
    info.update({"n_frames": int(n_frames), "natom_super": int(nsc),
                 "natom_unit": int(len(_read_poscar(out / "POSCAR"))),
                 "displacement_rms": rms,
                 "equilibrium_subtracted": bool(eq_applied)})

    # ---- phono3py engine on an array/pkl dataset: synthesise the params YAML ----
    #   Layouts 4-5 carry no YAML, but the phono3py engine reads one.  Build a
    #   type-2 (random-displacement) phono3py_params.yaml from the canonical
    #   arrays so a bare .npy/.pkl dataset also works with FIT_ENGINE=phono3py.
    if (n_frames and str(cfg.get("engine") or "").lower() == "phono3py"
            and not (out / "phono3py_params.yaml").is_file()
            and (kind == "vasprun" or not (out / "phono3py_disp.yaml").is_file())):
        try:
            _build_phono3py_params(out, disps, forces, scm)
        except Exception as _e:  # noqa: BLE001
            print("[WARN] could not build phono3py_params.yaml: %s" % _e,
                  flush=True)

    (out / "fc_dataset.json").write_text(
        json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")
    print("[OK] dataset normalised: %d frames x %d atoms, RMS displacement %s"
          % (n_frames, nsc, ("%.4f A" % rms) if rms else "n/a"), flush=True)
    print("[DONE] prep", flush=True)


def _equilibrium_from_npy_pair(src, nsc):
    """Reference-frame forces from dataset_disps/forces.npy next to the dataset.

    Those files end with an all-zero displacement frame by convention (the kl
    skills write it that way); its forces are the equilibrium residual.  Returns
    None when the pair is absent, mismatched, or has no such frame."""
    try:
        d = np.load(str(Path(src) / "dataset_disps.npy"))
        f = np.load(str(Path(src) / "dataset_forces.npy"))
    except Exception:
        return None
    if d.ndim != 3 or f.ndim != 3 or len(d) != len(f) or len(d) == 0:
        return None
    if d.shape[1] != nsc or not np.allclose(d[-1], 0.0):
        return None
    return np.asarray(f[-1], float).copy()


def _read_equilibrium_npy(cfg, out, nsc):
    p = cfg.get("equilibrium_forces_npy")
    if not p:
        return None
    p = Path(p)
    if not p.is_absolute():
        p = out / p
    if not p.is_file():
        return None
    arr = np.load(str(p))
    if arr.shape != (nsc, 3):
        print("[WARN] equilibrium_forces_npy %s has shape %s, expected %s"
              % (p, arr.shape, (nsc, 3)))
        return None
    return np.asarray(arr, float)


def _diag_supercell_matrix(out):
    """Fall back to the diagonal ratio of SPOSCAR/POSCAR lengths."""
    try:
        uc = _read_poscar(out / "POSCAR")
        sc = _read_poscar(out / "SPOSCAR")
        ru = np.linalg.norm(uc.cell, axis=1)
        rs = np.linalg.norm(sc.cell, axis=1)
        reps = [int(round(rs[i] / ru[i])) for i in range(3)]
        print("[WARN] supercell_matrix not recorded in the source YAML -- guessing "
              "diagonal %s from the POSCAR/SPOSCAR edge lengths" % reps, flush=True)
        return [[reps[0], 0, 0], [0, reps[1], 0], [0, 0, reps[2]]]
    except Exception as e:
        print("[WARN] could not guess the supercell matrix: %s" % e)
        return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def _load_dataset(out):
    """Read the canonical dataset written by prep."""
    d = np.load(out / "dataset_disps.npy")
    f = np.load(out / "dataset_forces.npy")
    return np.asarray(d[:-1], float), np.asarray(f[:-1], float)


# ==========================================================================
# Shell / third-order cutoff determination (ported from kl-dft-cpu S5_fc)
# --------------------------------------------------------------------------
# The candidate cutoffs come from the gen (midpoints between adjacent neighbour
# shells of the unit cell); here the fit decides how far the data can actually
# determine the fc3.  For pheasy every candidate is refitted and a frame
# bootstrap gives, per shell, the scatter of |Phi^3| -- the outermost shell
# whose sigma/|mean| is below the threshold fixes stable_upper_cut.  The other
# engines get a single-fit per-shell report (no statistics).
# ==========================================================================
def _min_image_shifts(cell, r_max):
    """Periodic image shifts (-k..k per axis) needed to cover r_max."""
    import itertools
    cell = np.asarray(cell, float)
    vol = abs(np.linalg.det(cell))
    ks = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        area = float(np.linalg.norm(np.cross(cell[j], cell[k])))
        perp = vol / area if area > 1e-12 else r_max
        ks.append(max(1, int(math.ceil(float(r_max) / max(perp, 1e-6)))))
    out = [(0, 0, 0)]
    for s in itertools.product(*[range(-k, k + 1) for k in ks]):
        if any(s):
            out.append(tuple(s))
    return out


def fc3_shell_index(cell, frac, tol=0.05, r_max=None):
    """Return (shell_dist, tri) for a supercell fc3.

    shell_dist is the sorted list of clustered atom-pair distances; tri[i,j,k]
    is the shell of the triplet (i,j,k) = the largest of its three pair shells
    (a self pair is -1)."""
    cell = np.asarray(cell, float)
    frac = np.asarray(frac, float)
    n = len(frac)
    if r_max is None:
        r_max = 12.0
    shifts = np.array(_min_image_shifts(cell, max(1.0, min(float(r_max), 20.0))),
                      dtype=float)
    D = np.full((n, n), np.inf)
    for shift in shifts:
        df = frac[None, :, :] + shift[None, None, :] - frac[:, None, :]
        cart = df @ cell
        D = np.minimum(D, np.sqrt((cart ** 2).sum(-1)))
    np.fill_diagonal(D, 0.0)
    iu = np.triu_indices(n, k=1)
    ds = np.sort(D[iu])
    ds = ds[ds < float(r_max)]
    shell_edges = []
    for x in ds:
        if not shell_edges or x - shell_edges[-1][-1] > tol:
            shell_edges.append([x])
        else:
            shell_edges[-1].append(x)
    shell_dist = [float(np.mean(s)) for s in shell_edges]
    if not shell_dist:
        return [], np.full((n, n, n), -1, dtype=int)
    idx = np.clip(np.searchsorted(shell_dist, D), 0, len(shell_dist) - 1)
    idx = np.where(D < 1e-9, -1, idx)
    a = idx[:, :, None]
    b = idx[:, None, :]
    c = np.transpose(idx, (1, 0))[None, :, :]
    tri = np.maximum(np.maximum(a, b), c)
    return shell_dist, tri


def fc3_shell_bins(tri, nshell):
    """Per-shell flat indices (streamed so one fc3 stays resident at a time)."""
    flat = np.asarray(tri).reshape(-1)
    return [np.nonzero(flat == s)[0] for s in range(int(nshell))]


def fc3_shell_means(fc3, bins):
    """Per-shell mean |Phi^3| of one force-constant array."""
    F3 = np.asarray(fc3, dtype=float)
    Fn = np.sqrt((F3 ** 2).sum(axis=(3, 4, 5))).reshape(-1)
    return [float(Fn[b].mean()) if len(b) else 0.0 for b in bins]


def fc3_shell_stats_from_means(shell_dist, per_shell_mean, stability_thr=0.3,
                               tol=0.05):
    """Per-shell sigma/|mean| and the largest cutoff the data can determine."""
    if not shell_dist:
        return [], float("inf")
    ns = len(shell_dist)
    arr = np.asarray(per_shell_mean, dtype=float)
    stats = []
    for s in range(ns):
        col = arr[:, s] if arr.size else np.zeros(1)
        mu = float(col.mean())
        sd = float(col.std(ddof=1)) if len(col) > 1 else 0.0
        stats.append({"shell": round(shell_dist[s], 4), "mean_abs": mu,
                      "std_abs": sd, "rel_std": (sd / abs(mu) if mu else None)})
    last_ok = -1
    for s, st in enumerate(stats):
        r = st["rel_std"]
        if r is not None and r < stability_thr:
            last_ok = s
        else:
            break
    if last_ok < 0:
        upper = shell_dist[0] - tol              # even the nearest shell is unstable
    elif last_ok >= ns - 1:
        upper = shell_dist[-1] + tol
    else:
        upper = 0.5 * (shell_dist[last_ok] + shell_dist[last_ok + 1])
    return stats, float(upper)


def _shell_cell_frac(out):
    from ase.io import read as _ase_read
    p = Path(out) / "SPOSCAR"
    if not p.is_file():
        p = Path(out) / "POSCAR"
    a = _ase_read(str(p))
    return np.asarray(a.cell, float), np.asarray(a.get_scaled_positions(), float)


def _recommend_cut(records):
    """Largest candidate inside its own stable_upper_cut (no kappa step here)."""
    ok = [r for r in records
          if r.get("stable_upper_cut") is not None
          and r.get("cut") is not None
          and r["cut"] <= r["stable_upper_cut"] + 1e-9]
    return max(r["cut"] for r in ok) if ok else None


def _write_cutoff_scan(out, payload):
    p = Path(out) / "cutoff_scan.json"
    p.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                 encoding="utf-8", newline="\n")
    print("[OK] cutoff_scan.json（%d 档，mode=%s）" % (len(payload.get("records") or []),
                                                       payload.get("mode")), flush=True)


def _clear_pheasy_cache(out):
    """Delete pheasy's cluster-space / null-space caches in ``out``.

    pheasy -s/-c reuse cs.pkl / ns_*.npz found in the working directory and
    only *warn* when the .meta.json sidecar is missing, so a second run with a
    different --c3 silently inherits the first cutoff's cluster space (the
    CUT3_SCAN bug: every candidate logged "2.73 / 15 clusters").  Always start
    each pheasy -s from a clean directory."""
    out = Path(out)
    for pat in ("cs.pkl", "cs.pkl.meta.json", "ns_*.npz", "ns_*.npz.meta.json"):
        for f in out.glob(pat):
            try:
                f.unlink()
            except OSError as e:
                sys.exit("[ERROR] cannot clear stale pheasy cache %s: %s" % (f, e))


def _check_cutoff_honoured(txt, c3v, what):
    """Abort when pheasy -s did not build the cluster space for the requested c3.

    pheasy prints "Cutoff distance (A) for 3-order IFCs: <c3>"; a mismatch (or a
    missing line while c3 was passed) means a stale cache was used, so any fc3
    from this run would belong to a different model than its label."""
    if c3v is None:
        return
    m = re.search(r"Cutoff distance \(A\) for 3-order IFCs:\s*([0-9.eE+-]+)", txt or "")
    if not m:
        print("[WARN] %s: pheasy -s log has no 3-order cutoff line; cannot verify "
              "c3=%.2f" % (what, c3v), flush=True)
        return
    got = float(m.group(1))
    if abs(got - float(c3v)) > 0.01 + 0.005 * abs(float(c3v)):
        sys.exit("[ERROR] %s: requested c3=%.2f but pheasy built the cluster space "
                 "for c3=%.2f (stale cs.pkl/ns cache?) -- refusing to continue"
                 % (what, float(c3v), got))


def _pheasy_scan(cfg, out, cands, base_for, rasr_flag, fit_flags, make_env):
    """Refit every candidate cutoff; frame-bootstrap the per-shell |Phi^3|.

    Runs in the step directory, in pheasy's own atom order (SPOSCAR as pheasy
    wrote it).  The fitted fc2/fc3 of the last candidate are overwritten by the
    nominal fit that cmd_fit_pheasy runs afterwards."""
    boot = int(cfg.get("cut3_bootstrap") or 0)
    thr = float(cfg.get("cut3_stability_thr") or 0.3)
    scan_dir = Path(out) / "cutoff_scan"
    scan_dir.mkdir(exist_ok=True)
    with open(Path(out) / "disp_matrix.pkl", "rb") as fh:
        d0 = np.asarray(pickle.load(fh))
    with open(Path(out) / "force_matrix.pkl", "rb") as fh:
        f0 = np.asarray(pickle.load(fh))
    nd = len(d0)
    rng = np.random.default_rng(int(cfg.get("pheasy_seed") or 20250924))

    def _run(label, phase, cmd, logfile=None):
        env = make_env(phase)
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env,
                           preexec_fn=_die_with_parent)
        txt = (r.stdout or "") + (r.stderr or "")
        (Path(out) / logfile).write_text(txt, encoding="utf-8") if logfile else None
        sys.stdout.write(txt)
        sys.stdout.flush()
        if r.returncode != 0:
            sys.exit("[ERROR] pheasy %s failed (rc=%d): %s" % (label, r.returncode, cmd))
        if label == "symmetry constraints" and rasr_flag:
            if not re.search(r"Imposing rotational invariance"
                             r"|Imposing equilibrium conditions", txt):
                sys.exit("[ERROR] RASR=%s requested but pheasy -c never logged "
                         "imposing it" % cfg.get("pheasy_rasr"))
        return txt

    def _fit_stages(c3v, tag):
        base = base_for(c3v)
        _clear_pheasy_cache(out)
        _s_txt = _run("cluster space", "setup", base + " -s", "pheasy_s_%s.log" % tag)
        _check_cutoff_honoured(_s_txt, c3v, "scan c3=%.2f" % c3v)
        _run("symmetry constraints", "setup", base + " -c" + rasr_flag,
             "pheasy_c_%s.log" % tag)
        _run("displacement matrix", "displacement",
             "%s -d --ndata %d --disp_file" % (base, nd), "pheasy_d_%s.log" % tag)
        rc, txt = _run_streaming("%s -f --ndata %d %s" % (base, nd, " ".join(fit_flags)),
                                 make_env("fit"))
        (Path(out) / ("pheasy_f_%s.log" % tag)).write_text(txt, encoding="utf-8")
        if rc != 0:
            sys.exit("[ERROR] pheasy -f (c3=%.2f) failed" % c3v)
        return txt

    records, shell_dist, bins = [], None, None
    for c in cands:
        tag = ("%.2f" % c).replace(".", "p")
        cdir = scan_dir / ("cut3_%s" % tag)
        cdir.mkdir(parents=True, exist_ok=True)
        fit_txt = _fit_stages(c, tag)
        # Bring pheasy's fc back to the DATASET atom order (and restore SPOSCAR)
        # before saving/copying: phono3py/ShengBTE/hiphive rebuild the supercell
        # from POSCAR + SUPERCELL in phonopy order, so leaving pheasy's order here
        # would silently give the downstream kappa wrong force constants.
        apply_pheasy_fc_order(cfg, out)
        if bins is None:
            # Shell index on the (now dataset-order) supercell.
            cell, frac = _shell_cell_frac(out)
            shell_dist, tri = fc3_shell_index(cell, frac, tol=0.05,
                                              r_max=max(cands) + 1.0)
            bins = fc3_shell_bins(tri, len(shell_dist))
        rec = {"cut": float(c), "natom_super": int(len(frac)), "ndata": int(nd),
               "shell_stats": [], "stable_upper_cut": None,
               "dir": str(cdir.relative_to(out))}
        m = _pheasy_metrics(fit_txt)
        rec["train_rel_err"] = m.get("pheasy_relative_error")
        rec["free_ifcs"] = m.get("pheasy_free_ifcs")
        rec["err_kind"] = "in-sample"          # pheasy -f residual on all frames
        # equations / parameters -- kl-dft-cpu S6's select_cutoff gates on
        # ratio >= 3 (data_ok); without this field every candidate is rejected.
        _free = m.get("pheasy_free_ifcs")
        rec["ratio"] = (float(nd) * 3.0 * len(frac) / float(_free)) if _free else None
        if (Path(out) / "fc2.hdf5").is_file():
            shutil.copyfile(str(Path(out) / "fc2.hdf5"), str(cdir / "fc2.hdf5"))
        if (Path(out) / "fc3.hdf5").is_file():
            shutil.copyfile(str(Path(out) / "fc3.hdf5"), str(cdir / "fc3.hdf5"))
        rels = [rec["train_rel_err"]] if rec["train_rel_err"] is not None else []
        # Only the bootstrap refits feed the stability statistics: with
        # CUT3_BOOTSTRAP=0 there is one fit per candidate and no spread to
        # measure, so stable_upper_cut must stay None (mirrors kl's S5 scan,
        # where the nominal fit is deliberately not counted as a sample).
        per_means = []
        for b in range(boot):
            idx = rng.integers(0, nd, size=nd)
            with open(Path(out) / "disp_matrix.pkl", "wb") as fh:
                pickle.dump(d0[idx], fh)
            with open(Path(out) / "force_matrix.pkl", "wb") as fh:
                pickle.dump(f0[idx], fh)
            base = base_for(c)
            _run("displacement matrix", "displacement",
                 "%s -d --ndata %d --disp_file" % (base, nd),
                 "pheasy_d_%s_b%02d.log" % (tag, b))
            rc, txt = _run_streaming("%s -f --ndata %d %s"
                                     % (base, nd, " ".join(fit_flags)), make_env("fit"))
            (Path(out) / ("pheasy_f_%s_b%02d.log" % (tag, b))).write_text(
                txt, encoding="utf-8")
            if rc != 0:
                print("[WARN] bootstrap #%d c3=%.2f -f failed; sample skipped"
                      % (b, c), flush=True)
                continue
            mm = _pheasy_metrics(txt)
            if mm.get("pheasy_relative_error") is not None:
                rels.append(mm["pheasy_relative_error"])
            apply_pheasy_fc_order(cfg, out)      # bootstrap fc -> dataset atom order
            try:
                fc3 = _read_fc3_any(out)
                if fc3 is not None:
                    per_means.append(fc3_shell_means(fc3, bins))
            except Exception as e:                       # noqa: BLE001
                print("[WARN] reading bootstrap fc3 failed: %s" % e, flush=True)
        with open(Path(out) / "disp_matrix.pkl", "wb") as fh:
            pickle.dump(d0, fh)
        with open(Path(out) / "force_matrix.pkl", "wb") as fh:
            pickle.dump(f0, fh)
        if per_means:
            stats, upper = fc3_shell_stats_from_means(shell_dist, per_means,
                                                      stability_thr=thr)
            rec["shell_stats"] = stats
            rec["stable_upper_cut"] = upper
        rec["stability_measured"] = boot > 0   # boot>0 但样本全失败 = 测了、不可用
        if len(rels) >= 2:
            rec["train_rel_se"] = float(np.std(rels, ddof=1))
        elif rels:
            rec["train_rel_se"] = 0.0
        records.append(rec)
        print("[..] c3=%.2f A: shells=%d, stable_upper_cut=%s"
              % (c, len(shell_dist or []), rec["stable_upper_cut"]), flush=True)
        _write_cutoff_scan(out, _scan_payload(cands, shell_dist, boot, thr, records))
    return records


def _scan_payload(cands, shell_dist, boot, thr, records):
    return {"mode": "scan", "candidates": [float(x) for x in cands],
            "shells": [float(x) for x in (shell_dist or [])],
            "bootstrap": int(boot), "stability_thr": float(thr),
            "recommended_cut": _recommend_cut(records),
            "records": records,
            "note": ("pheasy per-candidate refit + frame bootstrap: a candidate is "
                     "usable when cut <= its own stable_upper_cut (outermost shell "
                     "with sigma/|mean| < stability_thr).  The final cutoff choice "
                     "belongs to a downstream kappa step.")}


def _shell_single_report(cfg, out):
    """Per-shell report from the nominal fit (no bootstrap) for non-pheasy runs."""
    cands = [float(x) for x in (cfg.get("cut3_candidates") or [])]
    if not cands or (Path(out) / "cutoff_scan.json").is_file():
        return
    try:
        fc3 = _read_fc3_any(out, compact_ok=True)
        if fc3 is None:
            return
        cell, frac = _shell_cell_frac(out)
        shell_dist, tri = fc3_shell_index(cell, frac, tol=0.05,
                                          r_max=max(cands) + 1.0)
        if fc3.shape[0] != fc3.shape[1]:
            # compact fc3: only the rows of the primitive atoms are stored; by
            # translation symmetry their per-shell mean equals the full one.
            tri = tri[_fc3_p2s(out)]
        bins = fc3_shell_bins(tri, len(shell_dist))
        means = fc3_shell_means(fc3, bins)
    except Exception as e:                               # noqa: BLE001
        print("[..] shell report skipped: %s" % e, flush=True)
        return
    stats = [{"shell": round(shell_dist[s], 4), "mean_abs": means[s],
              "std_abs": None, "rel_std": None} for s in range(len(shell_dist))]
    payload = {"mode": "single", "candidates": [float(x) for x in cands],
               "shells": [float(x) for x in shell_dist],
               "bootstrap": 0,
               "stability_thr": float(cfg.get("cut3_stability_thr") or 0.3),
               "recommended_cut": None,
               "records": [{"cut": None, "natom_super": int(len(frac)),
                            "shell_stats": stats, "stable_upper_cut": None}],
               "note": ("nominal-fit shell report: per-shell mean|Phi^3| only.  "
                        "stable_upper_cut needs the pheasy refit + bootstrap scan "
                        "(set CUT3_SCAN=auto).")}
    _write_cutoff_scan(out, payload)


# ==========================================================================
# fit: phono3py + symfc / alm
# ==========================================================================
def _build_phono3py_params(out, disps, forces, scm):
    # phono3py engine needs a YAML; an array/pkl dataset (layouts 4-5) carries
    # only displacements+forces.  Synthesise phono3py_params.yaml (a type-2
    # random-displacement dataset) so a bare .npy/.pkl works with phono3py too.
    from phonopy.interface.vasp import read_vasp
    from phono3py import Phono3py
    out = Path(out)
    uc = read_vasp(str(out / "POSCAR"))
    ph3 = Phono3py(uc, supercell_matrix=np.asarray(scm, float).tolist())
    ph3.dataset = {"displacements": np.asarray(disps, float),
                   "forces": np.asarray(forces, float),
                   "natom": len(ph3.supercell)}
    ph3.save(str(out / "phono3py_params.yaml"))
    print("[OK] phono3py_params.yaml built from the array dataset "
          "(%d frames x %d atoms)" % (len(disps), len(ph3.supercell)), flush=True)
    return ph3

def _load_ph3(out):
    """Build a phono3py object carrying the forces.

    phono3py_params.yaml (written by prep or supplied by the source) embeds the
    forces and is preferred; phono3py_disp.yaml + FORCES_FC3 is the fallback."""
    import phono3py
    pm = out / "phono3py_params.yaml"
    if pm.is_file():
        try:
            ph3 = phono3py.load(str(pm), produce_fc=False, is_nac=False, log_level=0)
            ds = ph3.dataset or {}
            if (ds.get("forces") is not None
                    or (ds.get("first_atoms") and "forces" in ds["first_atoms"][0])):
                print("[..] forces read from %s" % pm.name, flush=True)
                return ph3
        except Exception as e:
            print("[WARN] phono3py.load(%s) failed: %s" % (pm.name, e), flush=True)
    dy = out / "phono3py_disp.yaml"
    if not dy.is_file():
        sys.exit("[ERROR] FIT_ENGINE=phono3py needs phono3py_disp.yaml or "
                 "phono3py_params.yaml in the dataset")
    ph3 = phono3py.load(str(dy), produce_fc=False, is_nac=False, log_level=1)
    ds = ph3.dataset or {}
    has = (ds.get("forces") is not None
           or (ds.get("first_atoms") and "forces" in ds["first_atoms"][0]))
    if not has:
        ph3 = phono3py.load(str(dy), forces_fc3_filename="FORCES_FC3",
                            produce_fc=False, is_nac=False, log_level=1)
        ds = ph3.dataset or {}
        has = ds.get("forces") is not None or bool(ds.get("first_atoms"))
    if not has:
        sys.exit("[ERROR] no forces available for the phono3py engine "
                 "(neither embedded in the YAML nor readable from FORCES_FC3)")
    return ph3


def _phono3py_scan(cfg, out, ph3, calc, enable, fc2_full, p2s):
    # Refit fc3 at every candidate cutoff and save per-cut fc2/fc3 dirs.
    # phono3py's fc3 cutoff is a physical truncation, not a data limit (the
    # symfc/alm fit is determined), so every candidate is usable and S2_kappa
    # picks the cutoff by the kappa plateau.  Same cutoff_scan/cut3_<c>/
    # layout as the pheasy scan.
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
    cands = [float(x) for x in (cfg.get("cut3_candidates") or [])]
    if not cfg.get("cut3_scan") or len(cands) < 2 or enable < 3:
        return
    safe = cfg.get("cut3_safe_cutoff")
    safe = float(safe) if safe else max(cands)
    scan_dir = Path(out) / "cutoff_scan"
    scan_dir.mkdir(exist_ok=True)
    records = []
    print("[..] phono3py 逐档截断扫描（%d 档）：%s" % (len(cands), cands), flush=True)
    for c in cands:
        tag = ("%.2f" % c).replace(".", "p")
        cd = scan_dir / ("cut3_%s" % tag)
        cd.mkdir(exist_ok=True)
        try:
            ph3.produce_fc3(fc_calculator=calc,
                            fc_calculator_options="cutoff = %s" % c,
                            is_compact_fc=True)
        except Exception as e:  # noqa: BLE001
            print("[WARN] c3=%.2f fc3 refit failed: %s" % (c, e), flush=True)
            continue
        write_fc2_to_hdf5(fc2_full, filename=str(cd / "fc2.hdf5"))
        write_fc3_to_hdf5(ph3.fc3, filename=str(cd / "fc3.hdf5"), p2s_map=p2s)
        records.append({"cut": c, "natom_super": len(ph3.supercell),
                        "ratio": 999.0, "train_rel_err": None,
                        "train_rel_se": 0.0, "err_kind": "in-sample",
                        "stable_upper_cut": safe, "shell_stats": [],
                        "note": "phono3py per-cut fc3 refit"})
        print("  c3=%.2f A -> cutoff_scan/cut3_%s/" % (c, tag), flush=True)
    _write_cutoff_scan(out, {
        "mode": "scan", "engine": "phono3py",
        "candidates": [float(x) for x in cands],
        "shells": [float(x) for x in (cfg.get("cut3_shells") or [])],
        "stability_thr": float(cfg.get("cut3_stability_thr") or 0.3),
        "safe_cutoff": safe,
        "note": ("phono3py per-cut fc3 refits: the cutoff is a physical "
                 "truncation, every candidate is usable; S2_kappa picks by "
                 "the kappa plateau."),
        "records": records})

def cmd_fit_phono3py(cfg, out):
    import phono3py                                    # noqa: F401
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5

    calc = str(cfg.get("fc_calc") or "symfc").lower()
    if calc not in ("symfc", "alm"):
        sys.exit("[ERROR] FC_CALC must be symfc or alm")
    enable = int(cfg.get("enable_fc") or 3)
    cutoff = _as_float_or_none(cfg.get("fc3_cutoff"))

    ph3 = _load_ph3(out)
    print("[..] phono3py fit: calculator=%s fc3_cutoff=%s" % (calc, cutoff), flush=True)
    # fc2 stays full (N,N,3,3): small, and every phonopy/plot/ShengBTE consumer
    # reads it as is.  fc3 is written COMPACT (n_prim,N,N,3,3,3) with its
    # p2s_map: the full array is N/n_prim times larger (243-atom supercell:
    # 3.1 GB vs 0.7 MB compressed) and producing it needs >15 GB.  phono3py's
    # BTE reads compact fc3 natively; _read_fc3_any expands it on demand.
    ph3.produce_fc2(fc_calculator=calc, is_compact_fc=False)
    fc2_full = np.array(ph3.fc2, copy=True)
    opts = None if cutoff is None else "cutoff = %s" % cutoff
    ph3.produce_fc3(fc_calculator=calc, fc_calculator_options=opts,
                    is_compact_fc=True)
    p2s = np.asarray(ph3.primitive.p2s_map, dtype="int64")
    write_fc2_to_hdf5(fc2_full, filename=str(out / "fc2.hdf5"))
    write_fc3_to_hdf5(ph3.fc3, filename=str(out / "fc3.hdf5"), p2s_map=p2s)
    for f in ("fc2.hdf5", "fc3.hdf5"):
        if not (out / f).is_file():
            sys.exit("[ERROR] phono3py produced no %s" % f)
    _phono3py_scan(cfg, out, ph3, calc, enable, fc2_full, p2s)
    # NAC from a BORN file, then a self-describing params file for the gate and
    # for whoever picks the force constants up next (matches kl-dft-cpu S5_fc).
    born = out / "BORN"
    if born.is_file():
        try:
            from phonopy.file_IO import parse_BORN
            nac = parse_BORN(ph3.primitive, filename=str(born))
            if isinstance(nac, dict) and not nac.get("factor"):
                nac["factor"] = 14.399652
            ph3.nac_params = nac
            print("[OK] BORN -> nac_params", flush=True)
        except Exception as e:
            print("[WARN] could not read BORN, saving without nac_params: %s" % e,
                  flush=True)
    try:
        ph3.save(str(out / P3_SUB / "phono3py_params.yaml"))
    except Exception:
        (out / P3_SUB).mkdir(exist_ok=True)
        try:
            ph3.save(str(out / P3_SUB / "phono3py_params.yaml"))
        except Exception as e:
            print("[WARN] could not save phono3py_params.yaml: %s" % e, flush=True)
    _write_fit_metrics(out, {"phono3py_fc_calculator": calc,
                             "phono3py_fc3_cutoff": cutoff})
    print("[DONE] fit_phono3py: fc2.hdf5 + fc3.hdf5", flush=True)


# ==========================================================================
# fit: pheasy
# ==========================================================================
# Canonical names; the pre-rename spelling "RFE" is not accepted.  Kept local
# so the compute-node driver does not import the gen module.
PHEASY_METHODS = ("OLS", "LASSO", "ALASSO", "RFE-OLS", "RFE-OLS-TSQR", "RIDGE",
                   "ARDR", "RVM")


def normalize_pheasy_method(value):
    """Canonical spelling of a pheasy method name (upper case, no spaces)."""
    return str(value or "").strip().upper().replace(" ", "")


def _pheasy_env(method, phase, ncpu, natom_super, tuning="safe", ols_ridge=None,
                ols_maxiter=None, cv_max_iter=None, cv_tol=None, tsqr_criterion=None):
    """Environment for one pheasy sub-step.

    'phase' is 'setup' for -s/-c, 'displacement' for -d and 'fit' for -f: the
    phases want opposite thread settings (the displacement matrix is built with
    many independent jobs and one BLAS thread each; the fit wants one job with
    all the BLAS threads).

    tuning='safe' (default) sets only what affects memory layout, threading and
    cross-validation grouping -- nothing that can move the solution.
    tuning='kl' additionally reproduces the solver settings hard-coded in the
    kl-dft-cpu submit_fit_pheasy.tpl production template.

    Measured on the local test dataset (BaS, 250 atoms, 10 frames, fc2+fc3,
    c3=4.0, OLS, --rasr BHH), end to end through this driver:

        safe   relative error 0.0216, worst force correlation 0.9997 -> gate stable
        kl     relative error 0.5760, worst force correlation 0.8910 -> gate imaginary

    Bisected to a single variable: PHEASY_OLS_RIDGE=1e-4.  Everything else in
    the kl block (MAXITER, ATOL, BTOL, ILP64, TWOLEVEL) is harmless here.

    Root cause (pheasy core/optimizer.py:_ols_lsmr): the value is not a
    normalised sklearn ridge, it is Tikhonov damping
    damp = sqrt(ridge * n_samples) appended to the LSMR system, so its effect
    grows with the number of equations.  At ndata=7500 that is damp = 0.87,
    which swamps the design matrix.  Dose-response on one dataset, one SM
    cache, one seed:

        ridge      LSMR iters   relative error   correlation
        0               509         0.0216          0.9997
        1e-10           508         0.0216          0.9997
        1e-8            342         0.0306          0.9995
        1e-6             72         0.0676          0.9977
        1e-4             18         0.5805          0.8991

    ridge > 0 also disables pheasy's resident GPU OLS (it falls back to CPU).
    The four production templates that carried it (jzzn / a800 / 3090 /
    hanhai25 submit_fit_pheasy*.tpl) were corrected to PHEASY_OLS_RIDGE=0 on
    2026-09-11, so both tuning profiles now leave it at pheasy's default 0.
    PHEASY_OLS_RIDGE is still honoured as an explicit override -- do not raise
    it without reading the correlation pheasy reports at the end of the fit.
"""
    env = dict(os.environ)
    nblas = min(int(ncpu), 32)
    env.update({
        "PHEASY_SVD_THRESHOLD": "500",
        "PHEASY_ASR_COMBINED": "1",
        "PHEASY_ASR_LWORK_LIMIT": "1500000000",
        "PHEASY_SM_DTYPE": "float32",
        "PHEASY_SM_THR": "1e-12",
        "PHEASY_ASR_SPARSE": "1",
        "PHEASY_ASR_SPARSE_THR": "1e-10",
        "PHEASY_ASR_COL_BLOCK": "5000",
        "PHEASY_BLAS_THREADS": str(nblas),
        "PHEASY_USE_CELER": "1",
        # group the cross-validation folds by configuration, otherwise the
        # 3*natom rows of one frame leak between training and validation
        "PHEASY_CV_GROUP_SIZE": str(3 * int(natom_super)),
    })
    if method == "RFE-OLS":
        env.update({
            "PHEASY_USE_RFE": "1", "PHEASY_RFE_TWOLEVEL": "1",
            "MKL_INTERFACE_LAYER": "ILP64", "PHEASY_RFE_MKL": "1",
            "PHEASY_RFE_STEP": "0.1", "PHEASY_RFE_RIDGE_ALPHA": "1e-11",
            "PHEASY_RFE_CV": "5", "PHEASY_RFE_LSMR_MAXITER": "60000",
            "PHEASY_RFE_WARM_START": "1", "PHEASY_COLNORM_FRAMES": "24",
            "PHEASY_COLNORM_EXACT": "0", "PHEASY_RFE_ONE_SE": "1",
        })
    elif method == "OLS" and tuning == "kl":
        env.update({
            "MKL_INTERFACE_LAYER": "ILP64",
            "PHEASY_OLS_ATOL": "1e-6", "PHEASY_OLS_BTOL": "1e-6",
        })
    if method == "OLS":
        if ols_ridge is not None:
            env["PHEASY_OLS_RIDGE"] = str(ols_ridge)
        else:
            # Never inherit a stray ridge: it is sqrt(ridge*ndata) damping, not
            # a normalised ridge, and 1e-4 silently wrecked a production fit.
            env.pop("PHEASY_OLS_RIDGE", None)
        if ols_maxiter is not None:
            env["PHEASY_OLS_MAXITER"] = str(ols_maxiter)
    elif method == "RIDGE":
        env["PHEASY_USE_CELER"] = "0"
    # RFE-OLS-TSQR only: the feature-count criterion (aic | bic | cv).  pheasy's
    # default since the rename is aic with n = the force-component rows; 'cv'
    # reproduces RFE-OLS exactly.  Empty = leave pheasy's default.
    if tsqr_criterion:
        env["PHEASY_TSQR_CRITERION"] = str(tsqr_criterion).strip().lower()
    if phase == "displacement":
        env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1", "PHEASY_N_JOBS": str(int(ncpu))})
    if cv_max_iter is not None:
        # The LASSO/ALASSO CV is the most expensive part of a pheasy fit and
        # pheasy caps it at min(max_iter, 800) with tol >= 1e-3; on a large
        # system every CV solve then stops at the cap (relative KKT ~2e-2) and
        # the alpha ranking is a convergence artefact.
        env["PHEASY_CV_MAX_ITER"] = str(int(cv_max_iter))
    if cv_tol is not None:
        env["PHEASY_CV_TOL"] = str(float(cv_tol))
    if phase == "fit":
        env.update({"OPENBLAS_NUM_THREADS": str(nblas), "OMP_NUM_THREADS": str(nblas),
                    "MKL_NUM_THREADS": str(nblas), "PHEASY_N_JOBS": "1",
                    "LOKY_MAX_CPU_COUNT": "1", "PHEASY_DOT_THREADS": str(int(ncpu)),
                    "OMP_NESTED": "FALSE", "MKL_DYNAMIC": "FALSE"})
    else:
        env.update({"OPENBLAS_NUM_THREADS": str(nblas), "OMP_NUM_THREADS": str(nblas),
                    "MKL_NUM_THREADS": str(nblas)})
    return env


def _apply_cv_knobs(env, cfg):
    """PHEASY_CV_MAX_ITER / PHEASY_CV_TOL from step.conf.

    pheasy's cross-validation stops every alpha at min(max_iter, 400) iterations
    (resident LASSO backend) with tol = max(tol, 1e-3).  At that cap the CV curve
    is not resolved, so the selected alpha can sit on the grid boundary and the
    fit looks far worse than the model is -- measured on the 3090 (Mg4C60,
    512 atoms, 128 frames): 400 iterations, relative KKT ~5e-3,
    alpha_opt = grid minimum, 8.2 % relative error.  pheasy's own warning tells
    the user to raise/lower exactly these two; this makes them reachable from
    step.conf instead of the submit template.
    """
    for key, env_key in (("pheasy_cv_max_iter", "PHEASY_CV_MAX_ITER"),
                         ("pheasy_cv_tol", "PHEASY_CV_TOL")):
        v = str(cfg.get(key) or "").strip()
        if v:
            env[env_key] = v


def _run_streaming(cmd, env):
    """Run a command, echoing its output line by line while collecting it.

    The fit is the one step that can run for an hour and print only at the end;
    capture_output=True hid even those lines until the process exited, so a
    running job looked dead (measured on the 3090: a 40-minute dense-LASSO fit
    with an empty queue.out).  Returns (returncode, text).
    """
    proc = subprocess.Popen(cmd, shell=True, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1,
                            preexec_fn=_die_with_parent)
    chunks = []
    try:
        for line in proc.stdout:
            chunks.append(line)
            sys.stdout.write(line)
            sys.stdout.flush()
    finally:
        proc.stdout.close()
        rc = proc.wait()
    return rc, "".join(chunks)


def _reconcile_gpu_env(env, method, cfg):
    """Keep the requested GPU knobs consistent with the sensing matrix pheasy builds.

    PHEASY_GPU_LASSO_RESIDENT (pheasy's initial spelling:
    PHEASY_GPU_TWOLEVEL_LASSO) exists only for a TwoLevelSM input.  With the
    dense/CSR matrix pheasy builds for LASSO by default it raises

        NotImplementedError: Resident GPU LASSO/ALASSO requires TwoLevelSM input

    and it raises it *after* the sensing matrix has been materialised, i.e. tens
    of minutes into the job (measured on the 3090: the run died 2 min into the
    fit, after a 10.8 GB sm_dense.npy).  run_pheasy.py only builds a TwoLevelSM
    when PHEASY_LASSO_TWOLEVEL=1 (default 0 for LASSO, 1 for OLS via
    PHEASY_OLS_TWOLEVEL), so asking for the GPU on a LASSO fit has to select the
    two-level path as well -- otherwise the request is dropped with a warning
    instead of killing the job.

    Returns the keys that belong in fit_metrics.json.
    """
    info = {"pheasy_gpu_requested": False, "pheasy_gpu_resident": False,
            "pheasy_lasso_twolevel": False}
    ngpu = str(cfg.get("pheasy_ngpu") or "").strip()
    if ngpu:
        try:
            n = int(ngpu)
        except ValueError:
            sys.exit("[ERROR] PHEASY_NGPU must be an integer (GPU count), got %r"
                     % ngpu)
        if n < 1:
            sys.exit("[ERROR] PHEASY_NGPU must be >= 1, got %d" % n)
        # The explicit count wins over whatever the submit template exported, and
        # it is the same number the job's --gres was rendered from (the gen writes
        # that from PHEASY_NGPU).  1..6 devices is what pheasy's own acceptance
        # suite covers; the ordinals are CUDA ordinals inside the allocation,
        # which is what PHEASY_GPU_SM_DEVICES expects.
        env["PHEASY_GPU_SM_NGPU"] = str(n)
        env["PHEASY_GPU_SM_DEVICES"] = ",".join(str(i) for i in range(n))
        info["pheasy_ngpu"] = n
        print("[..] pheasy GPU: %d device(s), PHEASY_GPU_SM_DEVICES=%s"
              % (n, env["PHEASY_GPU_SM_DEVICES"]), flush=True)
    info["pheasy_gpu_requested"] = any(
        _truthy(env.get(k)) for k in ("PHEASY_USE_GPU", "PHEASY_GPU_SM"))
    # The *resident* flag is LASSO/ALASSO-only (pheasy raises
    # NotImplementedError for it elsewhere), but PHEASY_LASSO_TWOLEVEL is
    # honoured for the whole LASSO family: run_pheasy.py accepts it for
    # LASSO / ALASSO / RIDGE.  Gating both on ("LASSO", "ALASSO") silently
    # dropped the two-level request for RIDGE, which then materialised the dense
    # sensing matrix instead of wrapping SM_prime/NS in a TwoLevelSM -- measured
    # at c2=6.5/c3=4.5 that is a 89.5 GB dense SM plus a 64 GB dense NS, no
    # two-level matvec and therefore no GPU path at all.
    resident = False
    if method in ("LASSO", "ALASSO"):
        _res = str(cfg.get("pheasy_gpu_lasso_resident") or "").strip().lower()
        if _res in ("1", "true", "yes", "y", "on"):
            env["PHEASY_GPU_LASSO_RESIDENT"] = "1"
        elif _res in ("0", "false", "no", "n", "off"):
            # Explicitly off also clears pheasy's initial spelling of the flag.
            for k in ("PHEASY_GPU_LASSO_RESIDENT", "PHEASY_GPU_TWOLEVEL_LASSO"):
                env.pop(k, None)
        resident = any(_truthy(env.get(k)) for k in
                       ("PHEASY_GPU_LASSO_RESIDENT", "PHEASY_GPU_TWOLEVEL_LASSO"))
    raw = str(cfg.get("pheasy_lasso_twolevel") or "").strip().lower()
    explicit = raw in ("1", "true", "yes", "y", "on", "0", "false", "no", "n", "off")
    want_tl = _truthy(cfg.get("pheasy_lasso_twolevel"), False)
    if resident and explicit and not want_tl:
        for k in ("PHEASY_GPU_LASSO_RESIDENT", "PHEASY_GPU_TWOLEVEL_LASSO"):
            env.pop(k, None)
        print("[WARN] PHEASY_GPU_LASSO_RESIDENT dropped: pheasy's resident LASSO "
              "needs a TwoLevelSM input, and PHEASY_LASSO_TWOLEVEL is explicitly "
              "off.  Drop the flag or set PHEASY_LASSO_TWOLEVEL=true.",
              flush=True)
        return info
    if resident and not want_tl:
        want_tl = True
        print("[..] PHEASY_GPU_LASSO_RESIDENT is set: selecting the two-level "
              "sensing matrix as well (PHEASY_LASSO_TWOLEVEL=1), otherwise the "
              "resident backend refuses the dense/CSR matrix.", flush=True)
    if want_tl:
        env["PHEASY_LASSO_TWOLEVEL"] = "1"
        info["pheasy_lasso_twolevel"] = True
        info["pheasy_gpu_resident"] = bool(resident)
    return info


def resolve_rasr(value, dim):
    """Resolve PHEASY_RASR, including the 'auto' policy.

    RASR (Born-Huang rotational invariance + Huang equilibrium constraints) is
    read by pheasy ONLY in the null-space construction step (-c), and it is
    MANDATORY for a 2D slab: without it the ZA branch is unconstrained, omega
    goes linear near Gamma instead of q^2, and kappa is wrong.

    On a bulk crystal those conditions are not merely unnecessary -- they
    measurably wreck the fit.  MnIn2Se4 (R-3m, 189-atom 3x3x3, pheasy OLS,
    c3 = 5.16 A), one dataset, one run per setting:
        PHEASY_RASR = BHH  -> pheasy relative error 8.9 %, worst force
                              correlation 0.994, min freq -1.57 THz
                              (spurious imaginary modes from the poorer fc2)
        PHEASY_RASR = none -> 0.62 %, correlation 1.000, min freq -0.03 THz
    The cutoff was irrelevant (4.99 A gave the same numbers as 5.16 A), so
    'auto' resolves BHH for DIM=2d and none for bulk.  An explicit BH / H /
    BHH / none is returned unchanged.
    """
    v = str(value or "").strip()
    if v.lower() != "auto":
        return v
    return "BHH" if str(dim or "").strip().lower() == "2d" else "none"


def cmd_fit_pheasy(cfg, out):
    """Four pheasy CLI steps: cluster space / symmetry constraints / sensing
    matrix / fit.  Mirrors templates in kl-dft-cpu and _common/mlff."""
    method = normalize_pheasy_method(cfg.get("pheasy_method") or "RFE-OLS")
    tsqr_criterion = (str(cfg.get("pheasy_tsqr_criterion") or "").strip() or None)
    if method not in PHEASY_METHODS:
        sys.exit("[ERROR] PHEASY_FIT_METHOD must be one of %s"
                 % ", ".join(PHEASY_METHODS))
    binary = str(cfg.get("pheasy_bin") or "pheasy")
    _ensure_env_bin_on_path(binary)
    if not shutil.which(binary):
        sys.exit("[ERROR] %r is not on PATH -- check CONDA_ENV/CONDA_SH in "
                 "step.conf (the submit template activates it) and PHEASY_BIN "
                 "(%s)" % (binary, "pheasy-gpu"))
    dim = str(cfg.get("supercell") or "").strip()
    if not dim:
        scm = np.asarray((cfg.get("supercell_matrix") or []), float)
        if scm.size == 9:
            dim = " ".join(str(int(round(x))) for x in np.diag(scm))
    if not dim:
        sys.exit("[ERROR] unknown supercell for pheasy -- set SUPERCELL in step.conf")
    enable = int(cfg.get("enable_fc") or 3)
    tuning = str(cfg.get("pheasy_tuning") or "safe").strip().lower()
    if tuning not in ("safe", "kl"):
        sys.exit("[ERROR] PHEASY_TUNING must be 'safe' or 'kl'")
    ols_ridge = _as_float_or_none(cfg.get("pheasy_ols_ridge"))
    ols_maxiter = None
    _om = str(cfg.get("pheasy_ols_maxiter") or "").strip()
    if _om:
        try:
            ols_maxiter = int(_om)
        except ValueError:
            sys.exit("[ERROR] PHEASY_OLS_MAXITER must be an integer, got %r" % _om)
    if tuning == "kl":
        print("[WARN] PHEASY_TUNING=kl: reproducing the kl-dft-cpu production "
              "solver settings (MAXITER/ATOL/BTOL/ILP64/TWOLEVEL). The "
              "PHEASY_OLS_RIDGE=1e-4 that used to come with them is NOT set: "
              "it is sqrt(ridge*ndata) LSMR damping that measurably degraded "
              "the fit, and the production templates were corrected on "
              "2026-09-11.", flush=True)
    eps = float(cfg.get("null_space_eps") or 0.001)
    apply_pheasy_order(cfg, out)
    disps, forces = _load_dataset(out)

    c2 = _as_float_or_none(cfg.get("pheasy_c2_cutoff"))
    c3 = _as_float_or_none(cfg.get("pheasy_c3_cutoff"))
    wflag = "-w %d" % enable

    natom_super = len(_read_poscar(out / "SPOSCAR"))
    try:
        ncpu = int(os.environ.get("SLURM_CPUS_PER_TASK")
                   or len(os.sched_getaffinity(0)) or 4)
    except Exception:
        ncpu = int(os.environ.get("SLURM_CPUS_PER_TASK") or 4)

    def base_for(c3v):
        """pheasy command prefix for one fc3 cutoff (the scan refits each)."""
        flags = []
        if c2 is not None:
            flags += ["--c2", str(c2)]
        if c3v is not None and enable >= 3:
            flags += ["--c3", str(c3v)]
        return "%s --dim %s %s %s --eps %s" % (binary, dim, wflag, " ".join(flags), eps)

    base = base_for(c3)

    # fit step: -l LASSO for RFE is deliberate -- the RFE strategy itself is
    # selected through PHEASY_USE_RFE in the environment
    fit_flags = ["--full_ifc", "-l", ("LASSO" if method == "RFE-OLS" else method),
                 "--hdf5"]
    # RASR (Born-Huang rotational invariance + Huang equilibrium) is read by
    # pheasy ONLY in the null-space construction step (-c).  Passing it to -f was
    # a silent no-op: -f loads ns_*.npz and never re-imposes the constraints, so
    # on a 2D slab the ZA branch stayed unconstrained -- omega goes linear near
    # Gamma instead of q^2, the frequencies can all be positive (so the imaginary
    # frequency gate passes) and kappa is still wrong.  Flag goes on -c now, and
    # the step output is checked to confirm the constraint really got imposed.
    rasr = resolve_rasr(cfg.get("pheasy_rasr"), cfg.get("dim"))
    if str(cfg.get("pheasy_rasr") or "").strip().lower() == "auto":
        print("[..] PHEASY_RASR=auto -> %s (DIM=%s)"
              % (rasr, cfg.get("dim") or "unknown"), flush=True)
    rasr_flag = ""
    if rasr and rasr.lower() not in ("none", "false", "off"):
        rasr_flag = " --rasr %s" % rasr
    # ARDR / RVM are sparse Bayesian fits that prune on column scale, so they
    # always get --std (unit variance) regardless of PHEASY_STD.
    if _truthy(cfg.get("pheasy_std"), False) or method in ("ARDR", "RVM"):
        fit_flags.append("--std")
    seed = cfg.get("pheasy_seed")
    if seed is not None:
        fit_flags += ["--seed", str(int(seed))]
    # RIDGE belongs here too.  Without it RIDGE silently fell back to
    # run_pheasy's own defaults (CV=5, NALPHA=50, a 4-decade default grid), so
    # every pheasy_cv / pheasy_nmu / pheasy_mu_min / pheasy_mu_max in
    # fit_config.json was ignored -- measured: a ridge sweep logged
    # "[RIDGE-CV] alpha 27/50" with 5 folds while fit_config.json said
    # nmu=1, mu=4.052069, cv=2.  Everything downstream (which alpha the solver
    # ever sees, and therefore whether the fit is regularised at all) followed
    # from that one missing string.
    if method in ("LASSO", "ALASSO", "RFE-OLS", "RFE-OLS-TSQR", "RIDGE"):
        # pheasy parses these with argparse type=int, so "-8.0" is rejected
        # outright -- coerce here as well so a hand-written fit_config.json
        # cannot get past the gen step with a float and fail on the node.
        def _gi(key, dflt):
            v = cfg.get(key)
            if v in (None, ""):
                return int(dflt)
            try:
                return int(str(v).strip())
            except ValueError:
                sys.exit("[ERROR] %s must be an integer (alpha exponent for "
                         "pheasy's argparse), got %r" % (key, v))
        fit_flags += ["--cv", str(_gi("pheasy_cv", 5)),
                      "--nmu", str(_gi("pheasy_nmu", 40)),
                      "--mu_min", str(_gi("pheasy_mu_min", -8)),
                      "--mu_max", str(_gi("pheasy_mu_max", -5)),
                      "--max_iter", str(_gi("pheasy_max_iter", 100000)),
                      "--tol", str(float(cfg.get("pheasy_tol") or 1e-5))]
    fit_step = "%s -f --ndata %d %s" % (base, len(disps), " ".join(fit_flags))
    disp_step = "%s -d --ndata %d --disp_file" % (base, len(disps))

    # pheasy's cluster-space step rewrites SPOSCAR in its own atom order; keep the
    # dataset's copy so apply_pheasy_fc_order can restore it -- including during
    # the cutoff scan below, which must save pheasy's fc in the DATASET atom order.
    if (out / "SPOSCAR").is_file():
        shutil.copyfile(str(out / "SPOSCAR"), str(out / "SPOSCAR.dataset"))

    # ---- third-order cutoff scan (ported from kl-dft-cpu S5_fc) ----
    #   Refit every candidate cutoff and frame-bootstrap it BEFORE the nominal
    #   fit, so the nominal fit overwrites the last candidate's fc2/fc3 in cwd.
    #   cutoff_scan.json carries the per-shell |Phi^3| scatter and
    #   stable_upper_cut (= the largest cutoff the data can determine).
    _cands = [float(x) for x in (cfg.get("cut3_candidates") or [])]
    _scan_on = bool(cfg.get("cut3_scan")) and len(_cands) >= 2 and enable >= 3
    if _scan_on:
        def _make_env(phase):
            e = _pheasy_env(
                method, phase, ncpu, natom_super, tuning, ols_ridge, ols_maxiter,
                cv_max_iter=(str(cfg.get("pheasy_cv_max_iter") or "").strip() or None),
                cv_tol=(str(cfg.get("pheasy_cv_tol") or "").strip() or None),
                tsqr_criterion=tsqr_criterion)
            if phase == "fit":
                _apply_cv_knobs(e, cfg)
                _reconcile_gpu_env(e, method, cfg)
            return e
        print("[..] 三阶截断扫描：候选 %s，bootstrap=%d"
              % (_cands, int(cfg.get("cut3_bootstrap") or 0)), flush=True)
        _pheasy_scan(cfg, out, _cands, base_for, rasr_flag, fit_flags, _make_env)
    elif cfg.get("cut3_candidates"):
        print("[..] 单截断拟合（CUT3_SCAN=%s，候选 %d 档）"
              % (cfg.get("cut3_scan"), len(_cands)), flush=True)

    run_steps = (("cluster space", "setup", base + " -s"),
                 ("symmetry constraints", "setup", base + " -c" + rasr_flag),
                 ("displacement matrix", "displacement", disp_step))
    for label, phase, cmd in run_steps:
        print("[pheasy] %s: %s" % (label, cmd), flush=True)
        if label == "cluster space":
            _clear_pheasy_cache(out)
        # capture output: the -c step is where RASR has to show up, and we need
        # its log to prove the constraint was imposed (see below)
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                           env=_pheasy_env(method, phase, ncpu, natom_super,
                                           tuning, ols_ridge, ols_maxiter,
                                           tsqr_criterion=tsqr_criterion),
                           preexec_fn=_die_with_parent)
        _out = (r.stdout or "") + (r.stderr or "")
        sys.stdout.write(_out)
        if r.returncode != 0:
            sys.exit("[ERROR] pheasy %s failed (rc=%d)" % (label, r.returncode))
        if label == "cluster space" and enable >= 3:
            _check_cutoff_honoured(_out, c3, "nominal c3=%s" % c3)
        if label == "symmetry constraints" and rasr_flag:
            open("pheasy_c.log", "w").write(_out)
            if not re.search(r"Imposing rotational invariance"
                             r"|Imposing equilibrium conditions", _out):
                sys.exit("[ERROR] RASR=%s requested but pheasy -c never logged "
                         "imposing rotational invariance / equilibrium "
                         "(see pheasy_c.log) -- force constants are not trustworthy"
                         % rasr)
            print("[OK] RASR=%s imposed in the null-space construction step (-c)" % rasr,
                  flush=True)
    print("[pheasy] fit: %s" % fit_step, flush=True)
    env = _pheasy_env(method, "fit", ncpu, natom_super, tuning, ols_ridge,
                      ols_maxiter,
                      cv_max_iter=(str(cfg.get("pheasy_cv_max_iter") or "").strip() or None),
                      cv_tol=(str(cfg.get("pheasy_cv_tol") or "").strip() or None),
                      tsqr_criterion=tsqr_criterion)
    _apply_cv_knobs(env, cfg)
    gpu_info = _reconcile_gpu_env(env, method, cfg)
    rc, log_txt = _run_streaming(fit_step, env)
    if rc != 0:
        sys.exit("[ERROR] pheasy fit failed (rc=%d)" % rc)
    # Bring fc2/fc3 (+ SPOSCAR) back to the dataset atom order before anything
    # downstream reads them.
    apply_pheasy_fc_order(cfg, out)

    # pheasy checks the per-configuration force correlation itself and warns
    # when it drops.  Interpret it carefully:
    #   corr < ~0.5   the mapping between the dataset and pheasy's supercell is
    #                 broken -- usually the atom order (pheasy blocks the images
    #                 per primitive atom, ASE repeat() interleaves them), so the
    #                 fit is worthless even though pheasy exits 0.  Gate failure.
    #   ~0.5 .. 0.98  ordinary model error: an fc2-only fit of anharmonic data,
    #                 or cutoffs too small.  Recorded, not fatal -- otherwise
    #                 every legitimate ENABLE_FC=2 run would be rejected.
    metrics = _pheasy_metrics(log_txt)
    metrics["pheasy_method"] = method
    metrics.update(gpu_info)
    if gpu_info.get("pheasy_gpu_requested") and not metrics.get("pheasy_gpu_used"):
        print("[WARN] PHEASY_USE_GPU/PHEASY_GPU_SM were requested but pheasy's "
              "log shows no GPU code path: this fit ran on the CPU.  The GPU "
              "backend only accelerates the two-level sparse matvec, so for "
              "LASSO it needs PHEASY_LASSO_TWOLEVEL=true (see the skill README).",
              flush=True)
    if method == "OLS" and metrics.get("pheasy_lsmr_istop") == 7:
        print("[..] OLS solver stopped at the iteration limit (istop=7, itn=%s): "
              "the residual is at the cap, not converged to atol.  Raise "
              "PHEASY_OLS_MAXITER to keep iterating (it leaves the fit otherwise "
              "unchanged)." % metrics.get("pheasy_lsmr_itn"), flush=True)
    mort = metrics.get("pheasy_worst_force_correlation")
    if mort is not None and mort < 0.5:
        print("[FAIL] pheasy's per-configuration force correlation collapsed to "
              "%.3f: the fitted constants do not reproduce the training forces, so "
              "the dataset and pheasy's supercell almost certainly disagree on the "
              "atom order (pheasy blocks the periodic images per primitive atom, "
              "while ASE repeat() interleaves them).\n"
              "        Compare POSCAR + SUPERCELL against the ordering the dataset "
              "was produced with before trusting anything else." % mort, flush=True)
        (out / ".fit_gate_fail").write_text(
            "pheasy force correlation %.3f (atom-order mismatch?)\n" % mort)
    elif mort is not None and mort < 0.98:
        print("[..] pheasy force correlation %.3f -- ordinary model error (fc2 "
              "only, or cutoffs too small); widen PHEASY_C2_CUTOFF / "
              "PHEASY_C3_CUTOFF or fit fc3 as well if it matters" % mort, flush=True)
    _write_fit_metrics(out, metrics)

    # Relative-error hard gate: an atom-order mismatch (~0.10) or a too-tight
    # cutoff (~0.05) lands here first.  The old code only parsed the number into
    # metrics and never failed on it, so a scrambled fit could pass the phonon
    # gate (kappa then silently off by tens of percent).
    rel = metrics.get("pheasy_relative_error")
    _rel_thr = float(cfg.get("fit_rel_err_fail") or 0.02)
    if rel is not None and _rel_thr > 0 and rel > _rel_thr:
        print("[FAIL] pheasy relative error %.4f > %.4f -- the fitted constants do "
              "not reproduce the training forces; check the supercell atom order and "
              "the C2/C3 cutoffs before trusting kappa"
              % (rel, _rel_thr), flush=True)
        (out / ".fit_gate_fail").write_text("pheasy relative error %.4f\n" % rel)

    # LASSO/ALASSO selection sanity now comes from pheasy's own manifest
    # (>= 35c0467), not from re-deriving log10(alpha*) against a hard-coded
    # [mu_min, mu_max] grid.  --mu_min/--mu_max are ignored for LASSO/ALASSO when
    # alpha is auto-selected (only RIDGE uses them), so that comparison was both
    # meaningless and unable to tell a genuine grid-minimum pin apart from the
    # data simply supporting no L1 penalty.
    if method in ("LASSO", "ALASSO"):
        _mf = {}
        _mf_path = out / "fit_manifest.json"
        if _mf_path.is_file():
            try:
                _mf = json.loads(_mf_path.read_text(encoding="utf-8"))
            except Exception as _e:
                print("[WARN] fit_manifest.json unreadable (%s); treating it as "
                      "missing" % _e, flush=True)
        _results = _mf.get("results") or {}
        if "cv_selection" not in _results:
            # Old pheasy (< 35c0467) never measured the alpha -> 0 OLS reference,
            # so an alpha* at the grid minimum is indistinguishable from a real
            # selection.  That absence is itself the red flag: require the upgrade.
            print("[FAIL] fit_manifest.json carries no cv_selection: this pheasy "
                  "is too old (before 35c0467, no per-fold alpha->0 OLS reference) "
                  "to certify the LASSO/ALASSO selection.  Upgrade pheasy to "
                  ">= 35c0467 and refit.", flush=True)
            (out / ".fit_gate_fail").write_text(
                "pheasy too old: no cv_selection in fit_manifest.json\n")
        elif _results.get("cv_selected_ols_limit"):
            # The CV resolved the exact alpha -> 0 limit at least as well as every
            # resolved alpha: the delivered fit IS OLS, no L1 penalty.  Certified,
            # not a truncated grid edge.
            print("[OK] CV selected the exact OLS limit (alpha -> 0): the data "
                  "support no L1 penalty; the delivered fit is OLS.", flush=True)
        elif _results.get("alpha_at_grid_edge"):
            # alpha* pinned at the grid bottom while the CV curve was still
            # falling: the selection is not trustworthy.
            print("[FAIL] alpha* sits on the grid boundary and the CV curve was "
                  "still falling -- the LASSO selection is not trustworthy; add "
                  "frames or widen the alpha grid.", flush=True)
            (out / ".fit_gate_fail").write_text("alpha on grid boundary\n")

    for f in ("fc2.hdf5",):
        if not (out / f).is_file():
            sys.exit("[ERROR] pheasy produced no %s" % f)
    if enable >= 3 and not (out / "fc3.hdf5").is_file():
        sys.exit("[ERROR] pheasy produced no fc3.hdf5 (ENABLE_FC=3 requested)")
    print("[DONE] fit_pheasy: method=%s ndata=%d dim=%s"
          % (method, len(disps), dim), flush=True)


def _pheasy_metrics(log_txt):
    """Pull the numbers pheasy reports out of its log."""
    m = {}
    for key, pat, cast in (
            ("pheasy_rmse_eV_per_A", r"\bRMSE:\s*([\d.eE+-]+)", float),
            ("pheasy_relative_error", r"Relative error:\s*([\d.eE+-]+)", float),
            ("pheasy_worst_force_correlation", r"worst corr=([\d.]+)", float),
            ("pheasy_free_ifcs", r"Free IFC terms:\s*(\d+)", int),
            ("pheasy_best_alpha", r"best alpha=\s*([\d.eE+-]+)", float),
            # pheasy prints "- alpha_opt: <x>" (the dense path prints
            # "best alpha="), so without these two the boundary gate below and
            # the recorded metric both saw nothing for FIT_ENGINE=pheasy.
            ("pheasy_alpha_opt", r"alpha_opt:[ ]*([0-9.eE+-]+)", float),
            ("pheasy_alpha_min", r"alpha_min:[ ]*([0-9.eE+-]+)", float),
            ("pheasy_alpha_max", r"alpha_max:[ ]*([0-9.eE+-]+)", float),
            ("pheasy_lsmr_istop", r"\[LSMR\] istop=(\d+)", int),
            ("pheasy_lsmr_itn", r"\[LSMR\] istop=\d+ itn=(\d+)", int)):
        mm = re.search(pat, log_txt)
        if mm:
            try:
                m[key] = cast(mm.group(1))
            except ValueError:
                pass
    if "pheasy_best_alpha" not in m and "pheasy_alpha_opt" in m:
        m["pheasy_best_alpha"] = m["pheasy_alpha_opt"]
    # A LASSO alpha pinned to the grid edge is only meaningful when the CV
    # solver converged; pheasy says so itself when FISTA hits its cap.
    m["pheasy_cv_hit_cap"] = "FISTA did not converge" in log_txt
    m["pheasy_cv_warning"] = "[CV] WARNING" in log_txt
    # Which code path pheasy actually took.  The GPU backend only accelerates the
    # two-level sparse matvec, so a GPU request without one of these markers is a
    # CPU fit -- silent until now (an hour went into discovering that on the 3090,
    # where PHEASY_USE_GPU=1 + 4 allocated cards still produced a CPU fit).
    ev = []
    for name, pat in (
            ("gpu_sm_matvec", r"\[GPU-SM\]\s*SpMV on"),
            ("gpu_resident_lasso", r"\[gpu_resident\]"),
            ("gpu_cv_folds", r"\[GPU\]\s*RIDGE CV"),
            ("twolevel_sm", r"\[SM-twolevel\]")):
        if re.search(pat, log_txt):
            ev.append(name)
    if re.search(r"\[GPU-SM\]\s*SpMV unavailable", log_txt):
        ev.append("gpu_sm_unavailable")
    m["pheasy_gpu_evidence"] = ev
    m["pheasy_gpu_used"] = any(n.startswith("gpu_") and n != "gpu_sm_unavailable"
                               for n in ev)
    return m


def _write_fit_metrics(out, metrics):
    """Merge engine-specific numbers into step1_fit/fit_metrics.json.

    cmd_post folds this file into phonon_summary.json, so the fit quality is
    visible next to the gate verdict instead of buried in queue.out."""
    if not metrics:
        return
    p = out / "fit_metrics.json"
    old = {}
    if p.is_file():
        try:
            old = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            old = {}
    old.update(metrics)
    p.write_text(json.dumps(old, indent=2, ensure_ascii=False), encoding="utf-8",
                 newline="\n")


# ==========================================================================
# fit: hiphive
# ==========================================================================
def _hiphive_fit(A, y, method, alpha):
    """Fit the hiphive design matrix.  Works with either hiphive generation:
    the old 'hiphive.fitting.Optimizer' API or a plain scikit-learn regressor
    (hiphive >= 1.5 dropped the fitting module)."""
    method = method.lower()
    try:
        from hiphive.fitting import Optimizer            # hiphive < 1.5
        return np.asarray(Optimizer((A, y)).train().parameters, float), "hiphive.Optimizer"
    except ImportError:
        pass
    from sklearn.linear_model import LinearRegression, Ridge, Lasso, BayesianRidge
    if method in ("ols", "linear"):
        reg = LinearRegression(fit_intercept=False)
    elif method in ("ridge", "bayes", "bayesian"):
        reg = BayesianRidge(fit_intercept=False) if method in ("bayes", "bayesian") \
            else Ridge(alpha=float(alpha), fit_intercept=False)
    elif method == "lasso":
        reg = Lasso(alpha=float(alpha), fit_intercept=False, max_iter=20000)
    elif method in ("ard", "ardregression"):
        from sklearn.linear_model import ARDRegression
        reg = ARDRegression(fit_intercept=False)
    else:
        sys.exit("[ERROR] HIPHIVE_FIT_METHOD must be ols | ridge | lasso | ard | bayes")
    reg.fit(A, y)
    return np.asarray(reg.coef_, float).ravel(), "sklearn." + type(reg).__name__


def _min_periodic_width(cell):
    # Smallest interplanar spacing of the lattice (Angstrom).  hiphive's
    # ClusterSpace requires 2*cutoff < this width; an acute / layered primitive
    # cell can be too thin for a physical cutoff, which makes hiphive's orbit
    # code fail with a cryptic "(0, N) is not in list".
    lat = np.asarray(cell, float)
    V = abs(float(np.linalg.det(lat)))
    widths = []
    for i, j, k in ((0, 1, 2), (1, 0, 2), (2, 0, 1)):
        n = float(np.linalg.norm(np.cross(lat[j], lat[k])))
        if n > 0:
            widths.append(V / n)
    return min(widths) if widths else 0.0

def cmd_fit_hiphive(cfg, out):
    """Fit fc2/fc3 from the displaced configurations with hiphive.

    hiphive determines its own (cutoff-bounded) atom list; the configurations
    handed to it can be any supercell of the primitive cell, because hiphive
    aligns each structure internally.  We therefore feed the full dataset
    supercell with the documented per-atom 'displacements' and 'forces' arrays
    and lift the fitted parameters back onto the same supercell."""
    from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
    import hiphive

    enable = int(cfg.get("enable_fc") or 3)
    c2 = cfg.get("hiphive_cutoff2")
    c3 = cfg.get("hiphive_cutoff3")
    if c2 in (None, ""):
        sys.exit("[ERROR] HIPHIVE_CUTOFF2 is required for FIT_ENGINE=hiphive")
    cutoffs = [float(c2)] + ([float(c3)] if enable >= 3 else [])
    symprec = float(cfg.get("hiphive_symprec") or 1e-5)
    disps, forces = _load_dataset(out)
    uc = _read_poscar(out / "POSCAR")
    sc = _read_poscar(out / "SPOSCAR")
    nsample = int(cfg.get("hiphive_n_configs") or 0)
    if nsample and nsample < len(disps):
        idx = np.linspace(0, len(disps) - 1, nsample).round().astype(int)
        disps, forces = disps[idx], forces[idx]
        print("[..] hiphive uses %d of the available frames" % len(disps), flush=True)

    # Reference structure for the cluster space: the primitive cell derived from
    # the unit cell (identity when they coincide).
    import spglib
    from ase import Atoms
    prim_matrix = cfg.get("primitive_matrix")
    if prim_matrix and np.allclose(np.asarray(prim_matrix, float), np.eye(3)):
        prim = uc
    else:
        lat = np.asarray(uc.cell, float)
        std = spglib.standardize_cell((lat, uc.get_scaled_positions(),
                                       uc.get_atomic_numbers()), to_primitive=True,
                                      no_idealize=True, symprec=symprec)
        if std is None:
            print("[WARN] spglib could not find a primitive cell -- using POSCAR",
                  flush=True)
            prim = uc
        else:
            prim = Atoms(numbers=std[2], cell=std[0], scaled_positions=std[1], pbc=True)
    print("[..] hiphive: primitive=%d atoms, supercell=%d atoms, cutoffs=%s"
          % (len(prim), len(sc), cutoffs), flush=True)

    # hiphive needs 2*cutoff < the reference cell's smallest periodic width.  An
    # acute/layered primitive (e.g. a 7-atom rhombohedral cell with 3.5 A
    # interplanar spacing) is too thin for any physical cutoff and hiphive fails
    # deep inside orbit construction with a cryptic "(0, N) is not in list".
    # Catch it here with an actionable message; HIPHIVE_CELL=supercell opts into
    # using the dataset supercell as the reference instead (parameter blow-up).
    _w = _min_periodic_width(prim.cell)
    if 2.0 * max(cutoffs) >= _w and str(cfg.get("hiphive_cell") or "").lower() \
            in ("supercell", "sc", "sposcar"):
        print("[WARN] hiphive 参考胞改用数据集超胞（%d 原子）：原胞周期宽度 %.2f A "
              "< 2*cutoff %.2f A" % (len(sc), _w, max(cutoffs)), flush=True)
        prim = sc
        _w = _min_periodic_width(prim.cell)
    if 2.0 * max(cutoffs) >= _w:
        sys.exit("[ERROR] hiphive 参考胞的最小周期宽度 %.2f A < 2*cutoff %.2f A，"
                 "hiphive 的 ClusterSpace 建不起来（内部轨道构造会崩）。\n"
                 "        处理：① 换 FIT_ENGINE=phono3py 或 pheasy（推荐，本结构可用）；\n"
                 "              ② 把 HIPHIVE_CUTOFF2/3 降到 < %.2f A；\n"
                 "              ③ 设 HIPHIVE_CELL=supercell 用数据集超胞当参考胞"
                 "（参数数量暴增，谨慎）。" % (_w, 2.0 * max(cutoffs), _w / 2.0))

    cs = ClusterSpace(prim, cutoffs, symprec=symprec)
    container = StructureContainer(cs)
    for u, f in zip(disps, forces):
        a = sc.copy()                      # ideal supercell: hiphive reads the
        a.new_array("displacements", u)    # displacement/force arrays, positions
        a.new_array("forces", f)           # are only used for alignment
        container.add_structure(a)
    A, y = container.get_fit_data()
    print("[..] hiphive design matrix %s, target %s, %d parameters"
          % (A.shape, y.shape, cs.n_dofs), flush=True)

    params, backend = _hiphive_fit(A, y, str(cfg.get("hiphive_fit_method") or "ridge"),
                                   cfg.get("hiphive_alpha", 1e-10))
    print("[OK] fitted with %s" % backend, flush=True)

    if _truthy(cfg.get("hiphive_enforce_asr"), True):
        try:
            params = hiphive.enforce_rotational_sum_rules(
                cs, params, sum_rules=["Huang", "Born-Huang"])
            print("[OK] rotational sum rules projected onto the parameters", flush=True)
        except Exception as e:
            print("[WARN] rotational-sum-rule projection skipped: %s" % e, flush=True)

    fcp = ForceConstantPotential(cs, params)
    fcs = fcp.get_force_constants(sc)

    # fc2: dense array is small; write phonopy hdf5 (dataset 'force_constants')
    # and the phonopy text file for maximum compatibility.
    fc2 = np.asarray(fcs.get_fc_array(order=2, format="phonopy"), float)
    _write_fc2_hdf5(fc2, out / "fc2.hdf5")
    try:
        fcs.write_to_phonopy(str(out / "FORCE_CONSTANTS"), format="text")
    except Exception as e:
        print("[WARN] could not write FORCE_CONSTANTS: %s" % e, flush=True)

    if enable >= 3:
        # fc3 can be huge (n^3 * 27); use hiphive's streaming writer instead of
        # materialising the dense array through get_fc_array(order=3).
        try:
            fcs.write_to_phono3py(str(out / "fc3.hdf5"))
            print("[OK] fc3.hdf5 written by hiphive", flush=True)
        except Exception as e:
            print("[WARN] hiphive write_to_phono3py failed: %s" % e, flush=True)
        try:
            (out / SB_SUB).mkdir(exist_ok=True)
            fcs.write_to_shengBTE(str(out / SB_SUB / "FORCE_CONSTANTS_3RD"), prim)
            print("[OK] shengbte/FORCE_CONSTANTS_3RD written by hiphive", flush=True)
        except Exception as e:
            print("[WARN] hiphive write_to_shengBTE failed: %s" % e, flush=True)

    np.save(out / "fc2_parameters.npy", params)
    if not (out / "fc2.hdf5").is_file():
        sys.exit("[ERROR] hiphive produced no fc2.hdf5")
    _write_fit_metrics(out, {"hiphive_backend": backend,
                             "hiphive_parameters": int(cs.n_dofs),
                             "hiphive_design_matrix": "%d x %d" % (A.shape[0], A.shape[1]),
                             "hiphive_cutoffs": list(cutoffs),
                             "hiphive_sum_rules_enforced":
                                 bool(_truthy(cfg.get("hiphive_enforce_asr"), True))})
    print("[DONE] fit_hiphive", flush=True)


def _write_fc2_hdf5(fc2, path):
    """Write fc2 in the phono3py/phonopy hdf5 layouts (dataset 'fc2' plus the
    phonopy-compatible 'force_constants')."""
    try:
        from phono3py.file_IO import write_fc2_to_hdf5
        write_fc2_to_hdf5(fc2, filename=str(path))
        return
    except Exception:
        pass
    import h5py
    with h5py.File(str(path), "w") as h:
        h.create_dataset("fc2", data=fc2, compression="gzip")
        h.create_dataset("force_constants", data=fc2, compression="gzip")


def cmd_fit(cfg, out):
    engine = str(cfg.get("engine") or "phono3py").lower()
    print("[..] FIT_ENGINE=%s" % engine, flush=True)
    if engine == "phono3py":
        cmd_fit_phono3py(cfg, out)
    elif engine == "pheasy":
        cmd_fit_pheasy(cfg, out)
    elif engine == "hiphive":
        cmd_fit_hiphive(cfg, out)
    else:
        sys.exit("[ERROR] FIT_ENGINE must be phono3py | pheasy | hiphive")


# ==========================================================================
# post: export + imaginary-frequency gate
# ==========================================================================
def _read_fc2_any(out):
    """Read fc2 as a numpy array from fc2.hdf5 / FORCE_CONSTANTS / fc2.npy."""
    import h5py
    p = out / "fc2.hdf5"
    if p.is_file():
        with h5py.File(str(p), "r") as h:
            for k in ("fc2", "force_constants"):
                if k in h:
                    return np.asarray(h[k][()], float)
    p = out / "fc2.npy"
    if p.is_file():
        return np.asarray(np.load(str(p)), float)
    p = out / "FORCE_CONSTANTS"
    if p.is_file():
        from phonopy.file_IO import parse_FORCE_CONSTANTS
        return np.asarray(parse_FORCE_CONSTANTS(str(p)), float)
    sys.exit("[ERROR] no fc2 artifact found (fc2.hdf5 / FORCE_CONSTANTS / fc2.npy)")


def _fc3_p2s(out):
    """p2s_map of a compact fc3.hdf5 (n_prim,N,N,3,3,3), or None when full."""
    import h5py
    p = Path(out) / "fc3.hdf5"
    if not p.is_file():
        return None
    with h5py.File(str(p), "r") as h:
        if "fc3" not in h:
            return None
        shp = h["fc3"].shape
        if shp[0] == shp[1]:
            return None
        if "p2s_map" in h:
            return np.asarray(h["p2s_map"][()], dtype="int64")
    # Compact without a stored map: take it from the phono3py primitive.
    return np.asarray(_fc3_primitive(out).p2s_map, dtype="int64")


def _fc3_primitive(out):
    """phono3py primitive cell of the fit (the basis compact fc3 refers to)."""
    import phono3py
    out = Path(out)
    for y in (out / P3_SUB / "phono3py_params.yaml", out / "phono3py_params.yaml",
              out / "phono3py_disp.yaml"):
        if y.is_file():
            return phono3py.load(str(y), produce_fc=False, is_nac=False,
                                 log_level=0).primitive
    sys.exit("[ERROR] compact fc3.hdf5 needs a phono3py YAML to expand it")


def _expand_compact_fc3(fc3c, prim):
    """(n_prim,N,N,3,3,3) -> (N,N,N,3,3,3) by lattice translations."""
    try:
        from phono3py.phonon3.fc3 import compact_fc3_to_full_fc3
        return compact_fc3_to_full_fc3(prim, fc3c)
    except ImportError:
        pass
    # Older phono3py: fc3[t(s), t(j), t(k)] = fc3[s, j, k] for every lattice
    # translation t of the supercell (atomic_permutations).
    n = fc3c.shape[1]
    full = np.zeros((n, n, n, 3, 3, 3), dtype="double")
    p2s = np.asarray(prim.p2s_map, dtype="int64")
    for t in np.asarray(prim.atomic_permutations, dtype="int64"):
        for ip, s in enumerate(p2s):
            full[t[s]][np.ix_(t, t)] = fc3c[ip]
    return full


def _read_fc3_any(out, compact_ok=False, prim_dir=None):
    """fc3 array from fc3.hdf5; a compact file is expanded to full unless
    compact_ok (then the caller gets (n_prim,N,N,3,3,3) and uses _fc3_p2s).
    prim_dir: where the phono3py YAML for the expansion lives (default `out`;
    cutoff_scan/cut3_<c>/ holds only the hdf5 files)."""
    import h5py
    p = out / "fc3.hdf5"
    if not p.is_file():
        return None
    with h5py.File(str(p), "r") as h:
        if "fc3" not in h:
            return None
        fc3 = np.asarray(h["fc3"][()], float)
    if compact_ok or fc3.shape[0] == fc3.shape[1]:
        return fc3
    return _expand_compact_fc3(fc3, _fc3_primitive(prim_dir or out))


def _fc3_load_bytes(out):
    """Byte count of the FULL fc3 array once materialised (no load performed).

    A full (N,N,N,3,3,3) array is ~8*N^3*27 bytes -- ~29 GB for a 512-atom
    supercell; a compact file is counted at its expanded size, since the
    consumers that ask (ShengBTE export, residual) need the full array.
    """
    import h5py
    p = out / "fc3.hdf5"
    if not p.is_file():
        return None
    with h5py.File(str(p), "r") as h:
        if "fc3" in h:
            d = h["fc3"]
            n = int(d.size) * int(d.dtype.itemsize)
            if d.shape[0] != d.shape[1]:
                n = n // int(d.shape[0]) * int(d.shape[1])
            return n
    return None


def _wrap_scaled_positions(atoms, eps=1e-9):
    """Wrap atoms into the cell; hiphive's ShengBTE writer wants strictly
    in-cell fractional positions (0 <= s < 1)."""
    fr = atoms.get_scaled_positions(wrap=True)
    fr = np.where(fr >= 1.0 - eps, 0.0, fr)
    fr = np.where(fr < eps, 0.0, fr)
    atoms.set_scaled_positions(fr)
    return atoms


def _export_shengbte(cfg, out):
    """Export fc2/fc3 to ShengBTE, then check the set and write CONTROL.

    The force-constant files alone are not runnable: ShengBTE also needs a
    CONTROL whose natoms/scell agree with them.  _shengbte_finalize derives that
    CONTROL from the very cells the files were written from and records the
    checks in shengbte/shengbte_manifest.json."""
    ok = _export_shengbte_files(cfg, out)
    if ok:
        ok = _shengbte_finalize(cfg, out)
    return ok


# ShengBTE layout (Src/input.f90 read2fc / split_index, Test-VASP example):
#   FORCE_CONSTANTS_2ND  phonopy FORCE_CONSTANTS text for the FULL supercell;
#                        header ntot must equal natoms * scell(1)*scell(2)*scell(3)
#                        (e.g. 9 atoms x 3x3x1 -> "81 81"); supercell index =
#                        ((iatom*nz + iz)*ny + iy)*nx + ix  (x fastest, atom slowest)
#   FORCE_CONSTANTS_3RD  thirdorder format: atom indices of the UNIT cell (1..natoms)
#                        plus Cartesian R vectors of the 2nd/3rd atom's cells
#   CONTROL              &crystal lattvec/positions = that unit cell, scell = the
#                        supercell of FORCE_CONSTANTS_2ND
# A "supercell 2ND + unit-cell 3RD/POSCAR" set is therefore the correct one; folding
# the 2ND file down to the unit cell makes ShengBTE stop with "wrong number of force
# constants for the specified scell".


def _shengbte_scell(uc, sc):
    """Integer supercell matrix sc = M . uc (rows = lattice vectors)."""
    m = np.asarray(sc.cell) @ np.linalg.inv(np.asarray(uc.cell))
    mi = np.rint(m).astype(int)
    if not np.allclose(m, mi, atol=1e-3):
        return None
    return mi


def _read_fc2_text(path):
    lines = Path(path).read_text().split("\n")
    head = lines[0].split()
    n = int(head[0])
    fc = np.zeros((n, n, 3, 3))
    i = 1
    for _ in range(n * n):
        a, b = (int(x) - 1 for x in lines[i].split()[:2])
        fc[a, b] = [[float(x) for x in lines[i + k].split()[:3]] for k in (1, 2, 3)]
        i += 4
    return fc


def _write_fc2_text(fc, path):
    n = fc.shape[0]
    with open(path, "w") as f:
        f.write("%d %d\n" % (n, n))
        for a in range(n):
            for b in range(n):
                f.write("%d %d\n" % (a + 1, b + 1))
                for r in fc[a, b]:
                    f.write("%22.15f%22.15f%22.15f\n" % tuple(r))


def _fc3_text_stats(path):
    """(number of triplet blocks, largest atom index) of FORCE_CONSTANTS_3RD."""
    lines = [ln.split() for ln in Path(path).read_text().splitlines()]
    nblk, i, imax = int(lines[0][0]), 1, 0
    for _ in range(nblk):
        while i < len(lines) and not lines[i]:
            i += 1
        idx = lines[i + 3]                       # block no., R2, R3, "i j k"
        imax = max(imax, *(int(x) for x in idx[:3]))
        i += 4 + 27
    return nblk, imax


def _shengbte_ngrid(cfg, uc):
    spec = str(cfg.get("shengbte_ngrid") or "auto").split()
    if len(spec) == 3:
        return [max(1, int(x)) for x in spec]
    length = float(spec[0]) if spec and spec[0] not in ("auto", "") else 40.0
    rec = np.linalg.inv(np.asarray(uc.cell)).T        # rows = b_i / 2pi (1/A)
    return [max(1, int(round(length * np.linalg.norm(b)))) for b in rec]


def _shengbte_born(uc, dims, out):
    """(epsilon 3x3, born[natoms,3,3]) from BORN, or None."""
    born = out / "BORN"
    if not born.is_file():
        return None
    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
        from phonopy.file_IO import parse_BORN
        pa = PhonopyAtoms(symbols=uc.get_chemical_symbols(), cell=uc.cell,
                          scaled_positions=uc.get_scaled_positions())
        ph = Phonopy(pa, np.diag(dims))
        nac = parse_BORN(ph.primitive, filename=str(born))
        if len(nac["born"]) != len(uc):
            return None
        return np.asarray(nac["dielectric"]), np.asarray(nac["born"])
    except Exception as e:
        print("[WARN] BORN not usable for CONTROL (%s) -- nonanalytic=.false." % e,
              flush=True)
        return None


def _write_shengbte_control(path, cfg, uc, dims, ngrid, nac):
    elems = []
    for s in uc.get_chemical_symbols():
        if s not in elems:
            elems.append(s)
    types = [elems.index(s) + 1 for s in uc.get_chemical_symbols()]
    t = [float(x) for x in str(cfg.get("shengbte_t") or "300").split()]
    L = ["&allocations",
         "        nelements=%d," % len(elems),
         "        natoms=%d," % len(uc),
         "        ngrid(:)=%d %d %d" % tuple(ngrid),
         "&end", "&crystal",
         "        lfactor=0.1,",                       # lattvec in A -> nm
         ]
    for i, v in enumerate(np.asarray(uc.cell)):
        L.append("        lattvec(:,%d)=%.10f %.10f %.10f," % ((i + 1,) + tuple(v)))
    L.append("        elements=%s" % " ".join('"%s"' % e for e in elems))
    L.append("        types=%s," % " ".join(str(x) for x in types))
    for i, f in enumerate(uc.get_scaled_positions(wrap=True)):
        L.append("        positions(:,%d)=%.10f %.10f %.10f," % ((i + 1,) + tuple(f)))
    if nac is not None:
        eps, born = nac
        for j in range(3):
            L.append("        epsilon(:,%d)=%.6f %.6f %.6f," % ((j + 1,) + tuple(eps[:, j])))
        for a in range(len(uc)):
            for j in range(3):
                L.append("        born(:,%d,%d)=%.6f %.6f %.6f,"
                         % ((j + 1, a + 1) + tuple(born[a][:, j])))
    L.append("        scell(:)=%d %d %d" % tuple(dims))
    L += ["&end", "&parameters"]
    if len(t) == 3:
        L.append("        T_min=%g, T_max=%g, T_step=%g" % tuple(t))
    else:
        L.append("        T=%g" % t[0])
    L += ["        scalebroad=%g" % float(cfg.get("shengbte_scalebroad") or 1.0),
          "&end", "&flags",
          # 迭代解（CONV）+ 同批的 RTA 都写；汇总 CONV 优先、缺温度点回退 RTA
          "        convergence=.TRUE.,",
          "        isotopes=%s," % (".TRUE." if cfg.get("isotope", True) else ".FALSE."),
          "        nonanalytic=%s," % (".TRUE." if nac is not None else ".FALSE."),
          "        nanowires=.FALSE.", "&end", ""]
    Path(path).write_text("\n".join(L))


def _shengbte_finalize(cfg, out):
    """Check FORCE_CONSTANTS_2ND/_3RD against the unit cell, put the 2ND file in
    ShengBTE supercell order when needed, and write CONTROL + a manifest."""
    sbdir = out / SB_SUB
    f2, f3 = sbdir / "FORCE_CONSTANTS_2ND", sbdir / "FORCE_CONSTANTS_3RD"
    man = {"format": "ShengBTE",
           "layout": {"FORCE_CONSTANTS_2ND": "full supercell (natoms*prod(scell) atoms; "
                                             "ShengBTE read2fc requires this)",
                      "FORCE_CONSTANTS_3RD": "unit-cell atom indices + R vectors",
                      "POSCAR/CONTROL": "unit cell"},
           "checks": {}, "usable": False}
    try:
        uc = _read_poscar(sbdir / "POSCAR")
        sc = _read_poscar(out / "SPOSCAR")
        M = _shengbte_scell(uc, sc)
        if M is None or np.count_nonzero(M - np.diag(np.diag(M))):
            man["error"] = ("supercell matrix %s is not diagonal -- ShengBTE scell "
                            "cannot express it" % (None if M is None else M.tolist()))
            raise ValueError(man["error"])
        dims = [int(x) for x in np.diag(M)]
        nat = len(uc)
        man.update({"natoms": nat, "scell": dims, "supercell_atoms": nat * int(np.prod(dims))})
        n2 = int(f2.read_text().split(None, 1)[0])
        man["checks"]["fc2_header"] = n2
        if n2 != nat * int(np.prod(dims)):
            man["error"] = ("FORCE_CONSTANTS_2ND has %d atoms, expected natoms*prod(scell)"
                            " = %d*%d = %d" % (n2, nat, int(np.prod(dims)),
                                               nat * int(np.prod(dims))))
            raise ValueError(man["error"])
        # 判定+重排都交给 fc_common（按力常数的平移不变性判断文件现在是哪种顺序，
        # 幂等：已是 ShengBTE 顺序的文件——pheasy 原生输出、新版导出、重复补导——不会被
        # 再次打乱；kappa_driver / kl-mlff S4 用的是同一个函数）
        _here = str(Path(__file__).resolve().parent)       # fc_common.py 与本脚本同目录
        if _here not in sys.path:
            sys.path.insert(0, _here)
        from fc_common import shengbte_fc2_ensure_order_files
        od = shengbte_fc2_ensure_order_files(f2, sbdir / "POSCAR", out / "SPOSCAR")
        man["checks"]["fc2_order"] = od["action"]
        if od["action"] == "reordered":
            man["checks"]["fc2_reordered"] = True
            print("[OK] FORCE_CONSTANTS_2ND re-ordered to ShengBTE supercell order",
                  flush=True)
        if f3.is_file():
            nblk, imax = _fc3_text_stats(f3)
            man["checks"].update({"fc3_blocks": nblk, "fc3_max_index": imax})
            if imax > nat:
                man["error"] = ("FORCE_CONSTANTS_3RD uses atom index %d > natoms %d"
                                % (imax, nat))
                raise ValueError(man["error"])
        ngrid = _shengbte_ngrid(cfg, uc)
        if str(cfg.get("dim") or "").upper().startswith("2") and \
                len(str(cfg.get("shengbte_ngrid") or "auto").split()) != 3:
            # 2D：真空方向（最长晶格矢）只取 1 个 q 点
            ngrid[int(np.argmax(np.linalg.norm(np.asarray(uc.cell), axis=1)))] = 1
        nac = _shengbte_born(uc, dims, out)
        _write_shengbte_control(sbdir / "CONTROL", cfg, uc, dims, ngrid, nac)
        man.update({"ngrid": ngrid, "nonanalytic": nac is not None, "usable": True,
                    "control": "CONTROL"})
        if str(cfg.get("dim") or "").upper().startswith("2"):
            man["note_2d"] = ("ShengBTE normalises kappa by the full (vacuum) cell "
                              "volume: rescale by c / thickness for the in-plane value")
        print("[OK] shengbte/CONTROL written (natoms=%d scell=%s ngrid=%s nonanalytic=%s);"
              " FC_2ND %d atoms = %d x %d"
              % (nat, dims, ngrid, nac is not None, n2, nat, int(np.prod(dims))), flush=True)
    except Exception as e:
        man.setdefault("error", str(e))
        print("[WARN] ShengBTE set not runnable: %s" % man["error"], flush=True)
    (sbdir / "shengbte_manifest.json").write_text(
        json.dumps(man, ensure_ascii=False, indent=1))
    return bool(man["usable"])


def _export_shengbte_files(cfg, out):
    """Export fc2/fc3 to the ShengBTE text formats.

    pheasy already writes FORCE_CONSTANTS_2ND / FORCE_CONSTANTS_3RD from its
    compact cluster representation (no dense fc3 materialisation), so those
    native files are reused directly.  phono3py / hiphive go through the
    hiphive writer, guarded by FC3_LOAD_GB_LIMIT.
    """
    want = _truthy(cfg.get("export_shengbte"), True)
    if not want:
        print("[..] EXPORT_SHENGBTE=false -- skipping", flush=True)
        return None
    sbdir = out / SB_SUB
    sbdir.mkdir(exist_ok=True)
    f2, f3 = sbdir / "FORCE_CONSTANTS_2ND", sbdir / "FORCE_CONSTANTS_3RD"

    # 1) pheasy writes these directly from compact IFCs.  Require the complete
    #    set for the requested order; never reuse a partial/stale file from
    #    another engine or from a failed previous fit.
    engine = str(cfg.get("engine") or "").strip().lower()
    n2, n3 = out / "FORCE_CONSTANTS_2ND", out / "FORCE_CONSTANTS_3RD"
    enable = int(cfg.get("enable_fc") or 3)
    # Native text is written in pheasy's internal atom order; once
    # apply_pheasy_fc_order has put the hdf5 into the dataset order the native
    # text is stale, so force the hiphive re-export from the corrected hdf5.
    _reordered = (out / ".pheasy_fc_reordered.json").is_file()
    native_ok = (engine == "pheasy" and not _reordered and n2.is_file() and
                 (enable < 3 or n3.is_file()))
    if native_ok:
        shutil.copyfile(str(n2), str(f2))
        if enable >= 3:
            shutil.copyfile(str(n3), str(f3))
        try:
            uc = _read_poscar(out / "POSCAR")
            _wrap_scaled_positions(uc)
            from ase.io import write as ase_write
            ase_write(str(sbdir / "POSCAR"), uc, format="vasp", direct=True,
                      sort=False)
        except Exception as e:
            print("[WARN] could not write shengbte/POSCAR: %s" % e, flush=True)
            return False
        ok = f2.is_file() and (enable < 3 or f3.is_file())
        print("[%s] ShengBTE export reused native FORCE_CONSTANTS_2ND/_3RD "
              "(written from compact IFCs)" % ("OK" if ok else "WARN"), flush=True)
        return ok

    if f2.is_file() and f3.is_file():
        print("[OK] shengbte/ already written by the hiphive engine", flush=True)
        return True
    try:
        from hiphive import ForceConstants
    except Exception as e:
        print("[WARN] hiphive unavailable -- skipping the ShengBTE export: %s" % e,
              flush=True)
        return False
    try:
        uc = _read_poscar(out / "POSCAR")
        sc = _read_poscar(out / "SPOSCAR")
        fc2 = _read_fc2_any(out)
        fc3 = None
        _fc3_bytes = _fc3_load_bytes(out)
        _fc3_limit = float(cfg.get("fc3_load_gb_limit") or 8.0) * 1e9
        if _fc3_bytes is not None and _fc3_bytes > _fc3_limit:
            print("[WARN] fc3.hdf5 would need %.1f GB in memory (limit %.1f GB); "
                  "skipping the ShengBTE fc3 export.  Raise FC3_LOAD_GB_LIMIT if "
                  "this node has the memory." % (_fc3_bytes / 1e9, _fc3_limit / 1e9),
                  flush=True)
        else:
            fc3 = _read_fc3_any(out)
        if fc2.shape[0] != len(sc):
            print("[WARN] fc2 is not a full supercell array (%s vs %d atoms) -- "
                  "ShengBTE export needs full force constants; skipping"
                  % (fc2.shape, len(sc)), flush=True)
            return False
        arrays = {"fc2_array": fc2}
        if fc3 is not None:
            arrays["fc3_array"] = fc3
        fcs = ForceConstants.from_arrays(sc, **arrays)

        _wrap_scaled_positions(uc)
        try:
            _wrap_scaled_positions(fcs._supercell)
        except Exception:
            pass
        fcs.write_to_phonopy(str(f2), format="text")
        if fc3 is not None:
            fcs.write_to_shengBTE(str(f3), uc)
        from ase.io import write as ase_write
        ase_write(str(sbdir / "POSCAR"), uc, format="vasp", direct=True, sort=False)
        ok = f2.is_file() and (f3.is_file() if fc3 is not None else True)
        print("[%s] ShengBTE export %s" % ("OK" if ok else "WARN",
                                           "done" if ok else "incomplete"), flush=True)
        return ok
    except Exception as e:
        print("[WARN] ShengBTE export failed (the phono3py artifacts are unaffected): "
              "%s" % e, flush=True)
        return False


def _parse_min_freq(band_yaml):
    p = Path(band_yaml)
    if not p.is_file():
        return None
    fr = []
    for ln in p.read_text(errors="ignore").splitlines():
        m = re.match(r"\s*frequency:\s*([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)\s*$", ln)
        if m:
            try:
                fr.append(float(m.group(1).replace("d", "e").replace("D", "e")))
            except ValueError:
                pass
    return min(fr) if fr else None


def _cells_cfg(cfg, out):
    """(supercell_matrix, primitive_matrix) for the phonon machinery.

    fc_dataset.json wins: prep resolves the supercell there (from the dataset
    YAML when it records one, otherwise from the POSCAR/SPOSCAR edge ratio), so
    fit_config.json may legitimately hold no matrix at all.  Returns
    (None, None) when neither source knows it."""
    scm = pm = None
    p = out / "fc_dataset.json"
    if p.is_file():
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            d = {}
        scm, pm = d.get("supercell_matrix"), d.get("primitive_matrix")
    if scm is None:
        scm = cfg.get("supercell_matrix")
    if pm is None:
        pm = cfg.get("primitive_matrix")
    if not scm:
        return None, None
    return (np.asarray(scm, float),
            None if pm is None else np.asarray(pm, float))


def _build_phonopy(scm, pm, uc, fc2):
    """Phonopy object for the stability gate / band structure."""
    from phonopy import Phonopy
    ph = Phonopy(uc, supercell_matrix=np.asarray(scm, float),
                 primitive_matrix=pm)
    ph.force_constants = fc2
    return ph


def _apply_nac(out, ph):
    born = out / "BORN"
    if not born.is_file():
        return False
    try:
        from phonopy.file_IO import parse_BORN
        nac = parse_BORN(ph.primitive, filename=str(born))
        if isinstance(nac, dict) and not nac.get("factor"):
            nac["factor"] = 14.399652
        ph.nac_params = nac
        return True
    except Exception as e:
        print("[WARN] could not read BORN (%s) -- running without NAC" % e, flush=True)
        return False


def _stability_gate(cfg, out):
    """q-mesh minimum frequency + imag_policy near-Gamma classification (merged
    from kl-dft-cpu S5_fc).  2D uses the no-NAC verdict (the 3D Coulomb kernel
    fakes imaginary frequencies near Gamma) and an explicit integer mesh
    (vacuum axis = 1) so near-Gamma soft modes are not missed."""
    dim = str(cfg.get("dim") or "3d").lower()
    is2d = dim.startswith("2")

    def _bad(note):
        return {"tool_ok": False, "stable": False, "min_freq": None,
                "min_freq_nonac": None, "min_freq_nac": None, "nac_used": False,
                "is_2d": is2d, "status": "tool_error", "note": note}

    scm, pm = _cells_cfg(cfg, out)
    if scm is None:
        return _bad("no supercell matrix in fc_dataset.json or fit_config.json")
    try:
        fc2 = _read_fc2_any(out)
        uc = _phonopy_unitcell(out / "POSCAR")
    except Exception as e:
        return _bad("could not prepare the phonopy input: %s" % e)

    # 2D 真空轴：取最长的原胞格矢（近似 dim_common.detect_dimension；共享引擎不
    #   携带真空轴）。用于显式整数网格的真空轴设 1 + ZA 面内方向投影。
    vax = 2
    if is2d:
        _cell = np.asarray(uc.cell, float)
        vax = int(np.argmax(np.linalg.norm(_cell, axis=1)))

    def _mesh_numbers():
        if not is2d:
            return None
        m = [60, 60, 60]
        m[vax] = 1
        return m

    def _minfreq(with_nac):
        p = _build_phonopy(scm, pm, uc, fc2)
        got = _apply_nac(out, p) if with_nac else False
        _mn = _mesh_numbers()
        if _mn is None:
            p.run_mesh(mesh=60.0, with_eigenvectors=False, is_mesh_symmetry=True)
        else:
            p.run_mesh(mesh=_mn, with_eigenvectors=False, is_mesh_symmetry=True,
                       is_gamma_center=True)
        return float(np.min(p.get_mesh_dict()["frequencies"])), p, got

    try:
        mf_nonac, ph_nonac, _ = _minfreq(False)
    except Exception as e:
        return _bad("phonopy mesh (no NAC) failed: %s" % e)

    mf_nac = None
    if (out / "BORN").is_file():
        try:
            mf_nac, _, _ = _minfreq(True)
        except Exception as e:
            print("[WARN] mesh with NAC failed, using the no-NAC verdict: %s" % e,
                  flush=True)

    try:
        ph_nonac.auto_band_structure(plot=False, write_yaml=True,
                                     npoints=int(cfg.get("band_points") or 51),
                                     filename=str(out / "band-dft-cpu.yaml"))
        bf = _parse_min_freq(out / "band-dft-cpu.yaml")
        if bf is not None:
            mf_nonac = min(mf_nonac, bf)
    except Exception as e:
        print("[..] band structure skipped (does not affect the verdict): %s" % e,
              flush=True)

    # ---- imag_policy 单一裁判（近 Γ 声学支 warn/fail、区外/光学 fail）----
    #   只对 mesh 采样（显式整数网格的真空轴分量恒 0，q 范数不受真空轴影响，故
    #   vac_axis 传 None 即可）。
    try:
        md = ph_nonac.get_mesh_dict()
        freqs = [[float(x) for x in row] for row in md["frequencies"]]
        qpts = [[float(x) for x in q] for q in md["qpoints"]]
        imag = imag_policy.classify_imag(freqs, qpts, is2d, None, cfg)
    except Exception as e:
        return _bad("imag_policy classify failed: %s" % e)

    # ---- 2D ZA 二次性（P2-2）：只看最小频率不够，ZA 线性化时频率全正也能漏判 ----
    za = None
    if is2d:
        try:
            za = phonon_stability.za_check(cfg, ph_nonac, vax, is2d)
            if za is not None and "error" in za:
                print("[WARN] ZA 二次性检查没跑成（不拦）：%s" % za["error"], flush=True)
        except Exception as e:
            print("[WARN] ZA 二次性检查没跑成（不拦）：%s" % e, flush=True)

    stability_verdict, stable = phonon_stability.finalize(imag, za)
    _za_v = imag_policy.za_verdict(za)
    mf_used = imag["min_freq_THz"] if imag.get("min_freq_THz") is not None else mf_nonac
    nac_used = False
    parts = ["min_freq(no-NAC)=%.3f THz" % mf_nonac]
    if mf_nac is not None:
        parts.append("min_freq(NAC)=%.3f THz" % mf_nac)
    parts.append("imag_policy=%s/%s" % (imag["verdict"], imag["imag_class"]))
    if is2d:
        parts.append("ZA verdict=%s" % _za_v)
    note = "; ".join(parts) + " -> " + ("no significant imaginary frequency"
                                        if stable else "imaginary frequency present")
    return {"tool_ok": True, "stable": stable, "min_freq": mf_used,
            "min_freq_nonac": mf_nonac, "min_freq_nac": mf_nac, "nac_used": nac_used,
            "is_2d": is2d, "imag": imag, "za_exponent": za, "za_verdict": _za_v,
            "stability_verdict": stability_verdict,
            "status": ("stable" if stable else
                       ("imaginary" if imag["verdict"] == "fail" else "za_not_quadratic")),
            "note": note}


def _fit_rmse(cfg, out, fc_dir=None, quiet=False):
    """Engine-independent training residual: predict the forces of a few frames
    from the fitted force constants and compare with the dataset.

    fc_dir: where fc2.hdf5/fc3.hdf5 live (default: the step directory); the
    dataset and SPOSCAR always come from the step directory.  Used per
    cutoff_scan/cut3_<c>/ too (those force constants are in dataset order)."""
    n = int(cfg.get("fit_rmse_frames") or 0)
    if n <= 0:
        return {}
    fcd = Path(fc_dir) if fc_dir is not None else Path(out)
    try:
        from hiphive import ForceConstants
        from hiphive.calculators import ForceConstantCalculator
        disps, forces = _load_dataset(out)
        if len(disps) == 0:
            return {}
        # apply_pheasy_fc_order (called in cmd_fit_pheasy) already normalises
        # fc2/fc3.hdf5 and SPOSCAR back to the dataset atom order, so
        # dataset_*.npy, SPOSCAR and the force constants share one order.
        idx = np.linspace(0, len(disps) - 1, min(n, len(disps))).round().astype(int)
        sc = _read_poscar(out / "SPOSCAR")
        fc2 = _read_fc2_any(fcd)
        fc3 = None
        _fc3_bytes = _fc3_load_bytes(fcd)
        _fc3_limit = float(cfg.get("fc3_load_gb_limit") or 8.0) * 1e9
        if _fc3_bytes is not None and _fc3_bytes > _fc3_limit:
            print("[..] fc3.hdf5 too large to materialise (%.1f GB > %.1f GB); "
                  "evaluating the fc2-only residual" % (_fc3_bytes / 1e9, _fc3_limit / 1e9),
                  flush=True)
        else:
            fc3 = _read_fc3_any(fcd, prim_dir=out)
        arrays = {"fc2_array": fc2}
        if fc3 is not None:
            arrays["fc3_array"] = fc3
        fcs = ForceConstants.from_arrays(sc, **arrays)
        calc = ForceConstantCalculator(fcs)
        se = ss = 0.0
        cnt = 0
        for i in idx:
            a = sc.copy()
            a.positions += disps[i]
            a.calc = calc
            p = a.get_forces()
            se += float(((p - forces[i]) ** 2).sum())
            ss += float((forces[i] ** 2).sum())
            cnt += p.size
        rmse = float(np.sqrt(se / max(cnt, 1)))
        rel = float(np.sqrt(se / ss)) if ss > 0 else None
        if not quiet:
            print("[OK] fit residual over %d frame(s): RMSE %.5f eV/A (relative %.4f)"
                  % (len(idx), rmse, rel if rel is not None else float("nan")),
                  flush=True)
        return {"fit_rmse_eV_per_A": rmse, "fit_rmse_relative": rel,
                "fit_rmse_frames_used": int(len(idx))}
    except Exception as e:
        print("[..] fit-residual evaluation skipped: %s" % e, flush=True)
        return {}


def _scan_rmse(cfg, out):
    """Force RMSE of every cutoff_scan/cut3_<c>/ fit, written back into
    cutoff_scan.json (force_rmse_eV_per_A / force_rmse_relative per record) --
    the right-hand panel of S2's kappa_vs_cutoff figure."""
    p = Path(out) / "cutoff_scan.json"
    if not p.is_file() or int(cfg.get("fit_rmse_frames") or 0) <= 0:
        return
    try:
        scan = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return
    changed = False
    for r in scan.get("records") or []:
        if r.get("cut") is None:
            continue
        d = Path(out) / (r.get("dir") or ("cutoff_scan/cut3_%s"
                                          % ("%.2f" % float(r["cut"])).replace(".", "p")))
        if not (d / "fc2.hdf5").is_file():
            continue
        try:
            m = _fit_rmse(cfg, out, fc_dir=d, quiet=True)
        except (Exception, SystemExit) as e:          # never fail post for a panel
            print("[..] c3=%s: force RMSE skipped (%s)" % (r["cut"], e), flush=True)
            m = {}
        if m:
            r["force_rmse_eV_per_A"] = m["fit_rmse_eV_per_A"]
            r["force_rmse_relative"] = m["fit_rmse_relative"]
            changed = True
            print("[..] c3=%.2f A: force RMSE %.2f meV/A (relative %.2f%%)"
                  % (float(r["cut"]), 1e3 * m["fit_rmse_eV_per_A"],
                     100 * (m["fit_rmse_relative"] or float("nan"))), flush=True)
    if changed:
        p.write_text(json.dumps(scan, ensure_ascii=False, indent=2),
                     encoding="utf-8", newline="\n")


def _write_kl_bundle(out, cfg):
    # Assemble the small bundle kl-dft-cpu's S5_fc import mode consumes.
    # Written next to the step outputs as kl_bundle/ so the kl side can pull it
    # across clusters (tf gen push_paths) without dragging the hundreds of MB of
    # fit intermediates (sm_dense.npy, ns_anharm3.npz, cs.pkl ...).  Idempotent.
    if not _truthy(cfg.get("export_kl_bundle"), True):
        print("[..] EXPORT_KL_BUNDLE=false -- skipping kl_bundle", flush=True)
        return None
    out = Path(out)
    b = out / "kl_bundle"
    if b.exists():
        shutil.rmtree(str(b))
    b.mkdir(parents=True, exist_ok=True)
    for f in ("fc2.hdf5", "fc3.hdf5", "POSCAR", "phonon_summary.json",
              "fc_fit_summary.json", "cutoff_scan.json", "cutoff_selection.json",
              "fc_dataset.json", "phono3py_disp.yaml", "phono3py_params.yaml"):
        if (out / f).is_file():
            shutil.copyfile(str(out / f), str(b / f))
    for d in ("shengbte", "cutoff_scan"):
        s = out / d
        if s.is_dir():
            shutil.copytree(str(s), str(b / d))
    for s in sorted(out.glob("cut3_*")):
        if s.is_dir():
            shutil.copytree(str(s), str(b / s.name))
    # phono3py per-engine params live in a phono3py/ subdir for some engines
    p3 = out / "phono3py"
    if p3.is_dir():
        for f in ("phono3py_disp.yaml", "phono3py_params.yaml"):
            if (p3 / f).is_file() and not (b / f).is_file():
                shutil.copyfile(str(p3 / f), str(b / f))
    n = sum(1 for _ in b.rglob("*") if _.is_file())
    print("[OK] kl_bundle/（%d 文件）—— kl-dft-cpu 跨集群复用 fc2/fc3 用"
          % n, flush=True)
    return b

def cmd_post(cfg, out):
    # prep records the dataset audit; fit_config.json only carries the request.
    ds = {}
    p = out / "fc_dataset.json"
    if p.is_file():
        try:
            ds = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            ds = {}
    cfg = dict(cfg)
    for k, src_key in (("n_frames", "n_frames"), ("natom_super", "natom_super")):
        if ds.get(src_key) is not None:
            cfg[k] = ds[src_key]
    sb_ok = _export_shengbte(cfg, out)
    rmse = _fit_rmse(cfg, out)
    _scan_rmse(cfg, out)
    # Shell report: the pheasy scan already wrote the fuller cutoff_scan.json;
    # for phono3py / hiphive this records the per-shell mean |Phi^3| of the
    # nominal fit (no bootstrap -- stable_upper_cut needs the pheasy refits).
    _shell_single_report(cfg, out)
    g = _stability_gate(cfg, out)
    stable, tool_ok = g["stable"], g["tool_ok"]
    print("[%s] %s" % ("OK" if stable else "FAIL", g["note"]), flush=True)

    gate_fail = (out / ".fit_gate_fail").is_file()
    if gate_fail:
        stable = False
        print("[FAIL] a fit-quality gate flagged this run -- see the log above",
              flush=True)

    summary = {
        "FIT_DONE": True,
        "engine": cfg.get("engine"),
        "fit_method": {"phono3py": cfg.get("fc_calc"),
                       "pheasy": cfg.get("pheasy_method"),
                       "hiphive": cfg.get("hiphive_fit_method")}.get(
                           str(cfg.get("engine") or "").lower()),
        "enable_fc": cfg.get("enable_fc"),
        "n_frames": cfg.get("n_frames"),
        "natom_super": cfg.get("natom_super"),
        "stable": bool(stable),
        "status": g["status"],
        "imaginary_frequency": bool(tool_ok and not stable),
        "min_frequency_THz": g["min_freq"],
        "min_freq_nonac_THz": g["min_freq_nonac"],
        "min_freq_nac_THz": g["min_freq_nac"],
        "nac_used_for_verdict": g["nac_used"],
        "is_2d": g["is_2d"],
        "shengbte_export": sb_ok,
        "fit_quality_gate_failed": bool(gate_fail),
        "tool_ok": tool_ok,
        "stability_verdict": g.get("stability_verdict"),
        "imag": g.get("imag"),
        "imag_class": (g.get("imag") or {}).get("imag_class"),
        "min_freq_THz": (g.get("imag") or {}).get("min_freq_THz"),
        "note": g["note"],
    }
    summary.update(rmse)
    try:
        summary.update(json.loads((out / "fit_metrics.json").read_text(encoding="utf-8")))
    except Exception:
        pass
    text = json.dumps(summary, indent=2, ensure_ascii=False)
    # phonon_summary.json is what tf's built-in "phonon" judge reads for this
    # step (tool_ok / stable / min_frequency_THz, giving the three-way
    # stable | imaginary | error verdict); fc_fit_summary.json is the skill's
    # own report and carries the same content.
    for name in ("phonon_summary.json", "fc_fit_summary.json"):
        (out / name).write_text(text, encoding="utf-8", newline="\n")
    print("[DONE] post: phonon_summary.json + fc_fit_summary.json ready "
          "(stable=%s, status=%s)" % (str(stable).lower(), g["status"]), flush=True)
    # kl-dft-cpu 跨集群复用：把 fc2/fc3 + shengbte 力常数打包成一个小 bundle，
    # 供 kl 的 S5_fc 导入模式在目标集群直接取用（见 kl skill 的 push_paths）。
    _write_kl_bundle(out, cfg)

    # A tool error means the gate itself could not run: fail the job so tf marks
    # the step as error.  A genuine imaginary frequency exits 0 on purpose --
    # the marker is simply unsatisfied and downstream steps stay held back.
    if not tool_ok:
        sys.exit("[ERROR] stability gate could not be evaluated -- see the log above")


# ==========================================================================
def cmd_shengbte(cfg, out):
    """Re-export / re-check shengbte/ (+ CONTROL) in a finished S1_fit directory
    without refitting -- e.g. in the fetched result/step1_fit/."""
    cfg = dict(cfg)
    cfg["export_shengbte"] = True
    # 从源头（pheasy 原生文件 / fc2/fc3.hdf5）重新导出，不复用 shengbte/ 里的旧文件：
    # 旧文件可能是老版本按 SPOSCAR 顺序写的，甚至被老版本错误重排过
    try:
        import hiphive  # noqa: F401
        _can_hdf5 = (out / "fc2.hdf5").is_file()
    except Exception:
        _can_hdf5 = False
    for f in ("FORCE_CONSTANTS_2ND", "FORCE_CONSTANTS_3RD"):
        p = out / SB_SUB / f
        if p.is_file() and ((out / f).is_file() or _can_hdf5):
            p.unlink()
    sys.exit(0 if _export_shengbte(cfg, out) else 1)


COMMANDS = {"prep": cmd_prep, "fit": cmd_fit, "post": cmd_post,
            "shengbte": cmd_shengbte}


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in COMMANDS:
        sys.exit("usage: fc_fit_driver.py <%s> fit_config.json"
                 % "|".join(COMMANDS))
    cfg = load_cfg(sys.argv[2])
    out = Path.cwd()
    COMMANDS[sys.argv[1]](cfg, out)


if __name__ == "__main__":
    main()
