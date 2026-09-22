#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fc_plot_phonon.py -- phonon band plot for the fc-fit skill (S2_plot).

Login-node step (run: gen, no job submitted).  Reads the force constants fitted
by S1_fit and draws the dispersion:

  <skill>/step1_fit/fc2.hdf5              fitted second-order force constants
  <skill>/step1_fit/POSCAR                unit cell
  <skill>/step1_fit/fc_dataset.json       supercell / primitive matrices
  <skill>/step1_fit/BORN                  optional; enables the NAC correction

Outputs (in <skill>/phonon_band_plot/):
  phonon_band_full.png        full frequency range
  phonon_band_lowfreq.png     0..10 THz zoom, where soft modes show up
  band-dft-cpu.yaml           the band structure phonopy computed
  phonon_band_summary.json    done marker + min/max frequencies
                              + ZA bending-branch exponent p per in-plane direction
                              + R^2 / residual / boundary margin / review flag

The imaginary-frequency gate in S1_fit already records the q-mesh minimum; this
step exists so a human can see where a soft mode lives.
"""
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SRC = "step1_fit"
OUTDIR = "phonon_band_plot"
LOWF_MAX = 10.0
FULL_PAD = 3.0
ZA_QMAX = 0.05            # reduced q cut for the ZA log-log fit (kl-dft-cpu default)
ZA_NQ = 12                # q points taken in (0, ZA_QMAX]
from za_2d import (ZA_P_RANGE, ZA_ZFRAC_MIN, ZA_MARGIN_MIN, ZA_R2_MIN,
                   ZA_SYMPREC_DEFAULT, ZA_SYMPREC_RELAX,
                   za_power_law, za_power_law_eig,
                   k_sum_sign as _k_sum_sign,
                   inplane_qdirs as _inplane_qdirs,
                   vacuum_axis_in_primitive as _vacuum_axis_in_primitive,
                   eig_out_of_plane as _eig_out_of_plane,
                   cell_symmetry as _cell_symmetry,
                   relaxed_symprec as _relaxed_symprec,
                   clone_phonopy as _clone_phonopy)


def _ensure_phonopy():
    """Re-run this script with an interpreter that actually has phonopy.

    The step runs wherever tf's remote PATH points, and that is not always a
    python that can plot: on the 3090 the gen python is the amset env
    (numpy/pymatgen/ase, *no* phonopy) while the mace-gpu env next to it has
    phonopy.  Without this the plot fails for an environment reason that has
    nothing to do with the fit.  FC_PLOT_PYTHON overrides the search and the
    FC_PLOT_REEXEC marker stops the recursion (a second failure surfaces as the
    plain ImportError below).
    """
    try:
        import phonopy  # noqa: F401
        return
    except Exception:
        pass
    if os.environ.get("FC_PLOT_REEXEC"):
        return
    import glob
    import subprocess
    home = os.path.expanduser("~")
    cands = [os.environ.get("FC_PLOT_PYTHON") or ""]
    for pat in ("miniconda3/envs/*/bin/python", "anaconda3/envs/*/bin/python",
                "miniforge3/envs/*/bin/python", "mambaforge/envs/*/bin/python",
                "/opt/miniconda3/envs/*/bin/python"):
        full = pat if pat.startswith("/") else os.path.join(home, pat)
        cands.extend(sorted(glob.glob(full)))
    seen = set()
    for p in cands:
        if not p or p in seen or not os.path.exists(p):
            continue
        seen.add(p)
        try:
            if os.path.realpath(p) == os.path.realpath(sys.executable):
                continue
            ok = subprocess.run([p, "-c", "import phonopy"],
                                capture_output=True, timeout=180).returncode == 0
        except Exception:
            continue
        if ok:
            print("[..] phonopy missing in %s; re-running the plot with %s"
                  % (sys.executable, p), flush=True)
            os.execve(p, [p, os.path.abspath(__file__)] + sys.argv[1:],
                      dict(os.environ, FC_PLOT_REEXEC="1"))


def _emit(result, code):
    print(json.dumps(result, ensure_ascii=False), flush=True)
    sys.exit(code)


def _read_fc2(path):
    import h5py
    with h5py.File(str(path), "r") as h:
        for k in ("fc2", "force_constants"):
            if k in h:
                return np.asarray(h[k][()], float)
    _emit({"status": "error", "reason": "%s has no fc2/force_constants dataset" % path}, 40)


def _matrices(src):
    """Supercell / primitive matrices recorded by the gen or by prep."""
    scm = pm = None
    for name in ("fc_dataset.json", "fit_config.json"):
        p = src / name
        if not p.is_file():
            continue
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        if scm is None and d.get("supercell_matrix"):
            scm = np.asarray(d["supercell_matrix"], float)
        if pm is None and d.get("primitive_matrix"):
            pm = np.asarray(d["primitive_matrix"], float)
    return scm, pm


def _import_dim_common():
    """Import the skill's existing dimension detector (skill/_common/opt).

    dim_common is the single implementation gen_step1_fit.resolve_dim uses; the
    plot reuses it instead of carrying a second vacuum heuristic.  On the cluster
    it is pushed next to this script by S2_plot's gen_need; in a plain checkout we
    fall back to the shared pool path.
    """
    try:
        import dim_common
        return dim_common
    except Exception:
        pass
    for cand in (Path(__file__).resolve().parent.parent / "_common" / "opt",
                 Path(__file__).resolve().parent.parent / "_common"):
        if (cand / "dim_common.py").is_file():
            sys.path.insert(0, str(cand))
            try:
                import dim_common
                return dim_common
            except Exception:
                return None
    return None


def _dim_and_axis(src):
    """(is_2d, cart_vacuum_axis, dim) from the dimension S1_fit already resolved.

    The fit gen froze DIM in fit_config.json (honouring an explicit 2d/3d
    setting, else dim_common.detect_dimension).  Re-read it and re-run the very
    same detector on the fitted POSCAR only to recover the vacuum axis -- no new
    heuristic is introduced.

    The second value is the CARTESIAN vacuum axis (an index of the POSCAR cell),
    returned even when the fit labelled the cell "3d".  It must be mapped into
    the primitive basis before it is used as a direction index -- see
    _vacuum_axis_in_primitive and the bug note on _inplane_qdirs.
    """
    dim = None
    for name in ("fit_config.json", "fc_dataset.json"):
        p = src / name
        if not p.is_file():
            continue
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        if dim is None and d.get("dim"):
            dim = str(d["dim"]).strip().lower()
    detected, axis = None, None
    dc = _import_dim_common()
    poscar = src / "POSCAR"
    if dc is not None and poscar.is_file():
        try:
            detected, axis, _ = dc.detect_dimension(str(poscar))
        except Exception as e:
            print("[WARN] dim_common.detect_dimension failed (%s)" % e)
    if dim is None:
        dim = detected
    is2d = (dim == "2d")
    # Return the Cartesian axis even for dim=="3d": the physical vacuum is a
    # property of the cell, and the primitive-basis mapping needs it.  Passing
    # this Cartesian index straight to _inplane_qdirs was the old bug.
    return is2d, (int(axis) if axis is not None else None), dim


def _za_summary(ph, src):
    """ZA bending-branch fields merged into phonon_band_summary.json.

    The ZA (flexural) branch is identified by its eigenvector: at every q the
    out-of-plane fraction of each branch is computed (Cartesian vacuum-axis
    component) and the branch used is the LOWEST-frequency branch carrying at
    least ZA_ZFRAC_MIN of its weight out of plane.  At Gamma the optical A1'/A2''
    branches are also ~100 % out of plane, so the global argmax would pick an
    optical mode; taking the lowest such branch is what makes this robust.  The
    old "lowest acoustic branch" result is kept side by side for comparison.

    The cell symmetry is audited at phonopy's default symprec and at
    ZA_SYMPREC_RELAX.  When the relaxed tolerance restores symmetry operations
    the cell is only numerically micro-distorted, and the ZA fit is redone on a
    clone built with that tolerance (the default-construction exponent stays in
    za_exponent_default_qX for comparison).  Without that repair the near-Gamma
    ZA branch comes out spuriously LINEAR, p ~ 1.19 instead of ~2.

    Never raises: a direction with no positive frequency (p is None) or a failed
    fit is recorded as None with a note, so the plot step still produces its
    marker.
    """
    res = {
        "za_qmax": ZA_QMAX,
        "za_n_qpoints": ZA_NQ,
        "za_p_range": [ZA_P_RANGE[0], ZA_P_RANGE[1]],
        "za_exponent_q1": None,
        "za_exponent_q2": None,
        "za_minfreq_q1_THz": None,
        "za_minfreq_q2_THz": None,
        "za_qdir_q1": None,
        "za_qdir_q2": None,
        "za_is_2d": None,
        "za_vacuum_axis": None,
        "za_vacuum_axis_in_primitive": None,
        # fit quality of the log-log power law, per direction
        "za_r2_q1": None,
        "za_r2_q2": None,
        "za_log_resid_rms_q1": None,
        "za_log_resid_rms_q2": None,
        # distance of p to the nearer edge of za_p_range, per direction
        "za_margin_q1": None,
        "za_margin_q2": None,
        # True when a bound-adjacent p or a noisy fit says "look by hand"
        "za_needs_review": None,
        # the review thresholds actually applied (traceability)
        "za_margin_min": ZA_MARGIN_MIN,
        "za_r2_min": ZA_R2_MIN,
        # eigenvector identification of the ZA branch (the new criterion)
        "za_branch_index_q1": None,
        "za_branch_index_q2": None,
        "za_zfrac_q1": None,
        "za_zfrac_q2": None,
        "za_zfrac_min_q1": None,
        "za_zfrac_min_q2": None,
        "za_branch_note_q1": "",
        "za_branch_note_q2": "",
        # the old criterion kept alongside: fit the lowest-frequency branch
        "za_lowest_exponent_q1": None,
        "za_lowest_exponent_q2": None,
        "za_lowest_r2_q1": None,
        "za_lowest_r2_q2": None,
        "za_lowest_log_resid_rms_q1": None,
        "za_lowest_log_resid_rms_q2": None,
        # symmetry audit: a numerically micro-distorted cell makes phonopy's
        # default symprec under-detect the space group and the ZA branch linear
        "za_symprec": ZA_SYMPREC_DEFAULT,
        "za_symprec_relax": ZA_SYMPREC_RELAX,
        "za_spacegroup_default": None,
        "za_spacegroup_relax": None,
        "za_symmetry_relaxed": None,
        "za_exponent_default_q1": None,
        "za_exponent_default_q2": None,
        # coarse machine verdict (backward compatible); see za_note
        "za_ok": None,
        "za_note": "",
    }
    try:
        is2d, cart_axis, dim = _dim_and_axis(src)
    except Exception as e:
        res["za_note"] = "dimension detection failed: %s" % e
        return res
    res["za_is_2d"] = bool(is2d)
    res["za_vacuum_axis"] = cart_axis
    # ★ direction bug fix: the Cartesian vacuum axis must be mapped into the
    # PRIMITIVE basis before it is used as a direction index.  The production
    # primitive_matrix puts the 25 A vacuum on primitive vector 0, so the old
    # ax = 2 selected (1,0,0)=vacuum + (0,1,0)=one in-plane direction.
    try:
        ax = _vacuum_axis_in_primitive(ph, cart_axis)
    except Exception as e:
        res["za_note"] = "vacuum axis -> primitive basis mapping failed: %s" % e
        return res
    res["za_vacuum_axis_in_primitive"] = int(ax)
    try:
        qdirs = _inplane_qdirs(ph, ax)
    except Exception as e:
        res["za_note"] = "in-plane direction selection failed: %s" % e
        return res
    res["za_qdir_q1"] = list(qdirs[0])
    res["za_qdir_q2"] = list(qdirs[1])

    # --- symmetry audit -----------------------------------------------------
    # A DFT-relaxed cell is hexagonal only to numerical precision, so phonopy's
    # default symprec can under-detect it; the nearest images are then no longer
    # symmetry averaged and the flexural ZA branch turns spuriously linear.  When
    # the relaxed tolerance restores operations the ZA fit is redone on a clone.
    notes = []
    res["za_spacegroup_default"] = _cell_symmetry(ph, ZA_SYMPREC_DEFAULT)
    res["za_spacegroup_relax"] = _cell_symmetry(ph, ZA_SYMPREC_RELAX)
    relaxed = _relaxed_symprec(ph, ZA_SYMPREC_RELAX)
    ph_za = ph
    if relaxed is not None:
        try:
            ph_za = _clone_phonopy(ph, relaxed)
            res["za_symprec"] = relaxed
        except Exception as e:
            relaxed = None
            notes.append("could not rebuild phonopy at symprec=%.0e (%s)"
                         % (ZA_SYMPREC_RELAX, e))
    res["za_symmetry_relaxed"] = bool(ph_za is not ph)
    if ph_za is not ph:
        notes.append(
            "phonopy default symprec=%.0e sees %s while symprec=%.0e sees %s: "
            "the cell is hexagonal only to numerical precision, so the ZA fit "
            "used the relaxed tolerance (compare za_exponent_default_q1/q2)"
            % (ZA_SYMPREC_DEFAULT, res["za_spacegroup_default"],
               ZA_SYMPREC_RELAX, res["za_spacegroup_relax"]))
    # ZA is the transverse (out-of-plane) branch: as in kl-dft-cpu's gate, fit
    # with the NAC kernel switched off so the near-Gamma LO-TO term cannot touch
    # it.  The primary result uses the eigenvector criterion (new); the old
    # lowest-branch fit is recorded alongside for comparison.
    had_nac = getattr(ph, "nac_params", None)
    tags = ("q1", "q2")
    ps, wmins, fits = [], [], []
    try:
        ph.nac_params = None
        try:
            ph_za.nac_params = None
        except Exception:
            pass
        for i, d in enumerate(qdirs):
            eig = None
            try:
                eig = za_power_law_eig(ph_za, qdir=d, qmax=ZA_QMAX, n=ZA_NQ,
                                       cart_vac_axis=cart_axis)
            except Exception as e:
                notes.append("qdir %s eigenvector ZA fit failed: %s"
                             % (list(d), e))
            try:
                lp, lw, lf = za_power_law(ph_za, qdir=d, qmax=ZA_QMAX, n=ZA_NQ)
            except Exception as e:
                lp, lw, lf = None, None, None
                notes.append("qdir %s lowest-branch fit failed: %s" % (list(d), e))
            if eig is not None:
                p, wmin, fit, extra = eig
                res["za_branch_index_%s" % tags[i]] = extra["branch_index"]
                res["za_zfrac_%s" % tags[i]] = extra["zfrac_first"]
                res["za_zfrac_min_%s" % tags[i]] = extra["zfrac_min"]
                res["za_branch_note_%s" % tags[i]] = extra["note"]
                if extra["note"]:
                    notes.append("%s %s" % (tags[i], extra["note"]))
            else:
                # no eigenvectors (e.g. the fake phonopy in tests): lowest branch
                p, wmin, fit = lp, lw, lf
            ps.append(p)
            wmins.append(wmin)
            fits.append(fit)
            res["za_lowest_exponent_%s" % tags[i]] = lp
            res["za_lowest_r2_%s" % tags[i]] = lf["r2"] if lf else None
            res["za_lowest_log_resid_rms_%s" % tags[i]] = (
                lf["log_resid_rms"] if lf else None)
            # p with the DEFAULT (un-repaired) construction, for comparison
            if ph_za is not ph:
                dflt = None
                try:
                    dflt = za_power_law_eig(ph, qdir=d, qmax=ZA_QMAX, n=ZA_NQ,
                                            cart_vac_axis=cart_axis)
                except Exception:
                    dflt = None
                res["za_exponent_default_%s" % tags[i]] = (
                    dflt[0] if dflt is not None else None)
            if p is None and wmin is not None:
                notes.append("qdir %s has no positive frequency "
                             "(omega_min=%.4f THz)" % (list(d), wmin))
    finally:
        try:
            ph.nac_params = had_nac
        except Exception:
            pass
    if ph_za is ph:
        res["za_exponent_default_q1"] = ps[0]
        res["za_exponent_default_q2"] = ps[1]
    res["za_exponent_q1"], res["za_exponent_q2"] = ps
    res["za_minfreq_q1_THz"], res["za_minfreq_q2_THz"] = wmins
    res["za_r2_q1"], res["za_r2_q2"] = [
        (f["r2"] if f else None) for f in fits]
    res["za_log_resid_rms_q1"], res["za_log_resid_rms_q2"] = [
        (f["log_resid_rms"] if f else None) for f in fits]
    lo, hi = ZA_P_RANGE
    # margin = how far p sits from the nearer bound of za_p_range; None when the
    # direction has no positive-frequency power law to fit.
    margins = [None if p is None else float(min(p - lo, hi - p)) for p in ps]
    res["za_margin_q1"], res["za_margin_q2"] = margins
    res["za_ok"] = bool(all(p is not None and lo < p < hi for p in ps))
    pstr = ["%.3f" % p if p is not None else "None" for p in ps]

    # Human-review hints.  The margin to a criterion bound and the R^2 of the
    # power law are the two ways the coarse za_ok can lie: p within
    # ZA_MARGIN_MIN of 1.7/2.3 can flip on a tiny change of dataset/q-window,
    # and a noisy log-log fit makes the fitted p uncertain whatever its value.
    # Both only make sense for a 2D layer (for 3D the criterion does not apply).
    review = False
    if is2d:
        for tag, p, fit, m in zip(("q1", "q2"), ps, fits, margins):
            if p is None or m is None:
                continue
            if m < ZA_MARGIN_MIN:
                review = True
                notes.append(
                    "%s p=%.3f is only %.3f from the [%.1f, %.1f] criterion "
                    "boundary (margin < %.2f)"
                    % (tag, p, m, lo, hi, ZA_MARGIN_MIN))
            r2 = fit.get("r2") if fit else None
            if r2 is not None and r2 < ZA_R2_MIN:
                review = True
                notes.append(
                    "%s log-log fit is too noisy (R2=%.4f < %.2f, log residual "
                    "RMS=%.4f), so its p is unreliable"
                    % (tag, r2, ZA_R2_MIN, fit.get("log_resid_rms")))
    res["za_needs_review"] = bool(review)
    if not is2d:
        notes.append("DIM=%s is not 2D; the quadratic-ZA criterion is only "
                     "meaningful for a 2D layer, exponents are informational"
                     % (dim or "unknown"))
    elif res["za_ok"]:
        notes.append("both in-plane directions have %.1f < p < %.1f" % (lo, hi))
    else:
        notes.append("the ZA bending branch is not quadratic in every direction")
    if review:
        notes.append("建议人工复核，不要仅凭 za_ok 下结论 (manual review "
                     "recommended: do not decide from the coarse machine "
                     "verdict za_ok alone)")
    res["za_note"] = ("ZA exponents p(q1,q2)=%s (za_ok = coarse machine "
                      "verdict); %s" % (pstr, "; ".join(notes)))
    return res


def main():
    _ensure_phonopy()
    root = Path.cwd().resolve()
    src = root / SRC
    if not src.is_dir():
        _emit({"status": "error", "reason": "missing %s -- run S1_fit first" % src}, 40)
    fc2p = src / "fc2.hdf5"
    if not fc2p.is_file():
        _emit({"status": "error", "reason": "missing %s (S1_fit did not finish)" % fc2p}, 40)

    scm, pm = _matrices(src)
    if scm is None:
        _emit({"status": "error",
               "reason": "no supercell_matrix in fc_dataset.json / fit_config.json"}, 40)

    from ase.io import read as ase_read
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    a = ase_read(str(src / "POSCAR"), format="vasp")
    # phonopy 2.x needs its own Atoms type, not an ASE Atoms
    uc = PhonopyAtoms(symbols=a.get_chemical_symbols(),
                      cell=np.asarray(a.cell, float),
                      scaled_positions=a.get_scaled_positions(wrap=True))
    ph = Phonopy(uc, supercell_matrix=scm,
                 primitive_matrix=None if pm is None else pm)
    ph.force_constants = _read_fc2(fc2p)
    nac = False
    born = src / "BORN"
    if born.is_file():
        try:
            from phonopy.file_IO import parse_BORN
            params = parse_BORN(ph.primitive, filename=str(born))
            if isinstance(params, dict) and not params.get("factor"):
                params["factor"] = 14.399652
            ph.nac_params = params
            nac = True
        except Exception as e:
            print("[WARN] BORN present but unusable (%s); plotting without NAC" % e)

    try:
        ph.auto_band_structure(plot=False, write_yaml=True, npoints=51,
                               filename=str(src / "band-dft-cpu.yaml"))
        bands = ph.get_band_structure_dict()
    except Exception as e:
        _emit({"status": "error", "reason": "phonopy band structure failed: %s" % e}, 40)

    freqs = np.asarray(bands["frequencies"], dtype=object)
    dists = np.asarray(bands["distances"], dtype=object)
    labels = list(bands.get("labels") or [])
    ticks = []
    for i, lb in enumerate(labels):
        if lb:
            ticks.append((float(dists[i][0]), str(lb).replace("$", "")))
    allf = np.concatenate([np.asarray(f, float).ravel() for f in freqs])
    fmin, fmax = float(allf.min()), float(allf.max())

    outdir = root / OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)

    def _draw(fname, lo, hi):
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        for q, f in zip(dists, freqs):
            # phonopy returns per-segment (n_qpoints, n_branches) arrays
            q = np.asarray(q, float)
            f = np.asarray(f, float)
            for ib in range(f.shape[1]):
                ax.plot(q, f[:, ib], "-", lw=1.0, color="#1f4e79")
        ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
        for x, lb in ticks:
            ax.axvline(x, color="0.7", lw=0.8)
        if len(ticks) > 1:
            ax.set_xticks([t[0] for t in ticks])
            ax.set_xticklabels([t[1] for t in ticks])
        ax.set_xlim(float(dists[0][0]), float(dists[-1][-1]))
        ax.set_ylim(lo, hi)
        ax.set_ylabel("Frequency (THz)")
        ax.set_title("Phonon dispersion%s" % (" (NAC)" if nac else ""), fontsize=11)
        fig.tight_layout()
        fig.savefig(str(outdir / fname), dpi=200)
        plt.close(fig)

    _draw("phonon_band_full.png", fmin - FULL_PAD, fmax + FULL_PAD)
    _draw("phonon_band_lowfreq.png", min(fmin - FULL_PAD, -1.0), LOWF_MAX)

    summary = {
        "status": "ok",
        "min_frequency_THz": fmin,
        "max_frequency_THz": fmax,
        "nac": nac,
        "n_qpoints": int(len(allf) // max(len(freqs[0][0]), 1)),
        "figures": ["phonon_band_full.png", "phonon_band_lowfreq.png"],
        "band_yaml": str((src / "band-dft-cpu.yaml").relative_to(root)),
        "stable_here": bool(fmin >= -0.10),
    }
    # ZA bending-branch exponent(s).  Strictly a 2D criterion, computed from the
    # force constants already in memory; any failure is recorded, never raised,
    # so the plot still writes its marker.
    try:
        summary.update(_za_summary(ph, src))
    except Exception as e:
        print("[WARN] ZA exponent analysis failed (%s); writing null fields" % e)
        summary.update({
            "za_qmax": ZA_QMAX, "za_n_qpoints": ZA_NQ,
            "za_p_range": [ZA_P_RANGE[0], ZA_P_RANGE[1]],
            "za_exponent_q1": None, "za_exponent_q2": None,
            "za_minfreq_q1_THz": None, "za_minfreq_q2_THz": None,
            "za_qdir_q1": None, "za_qdir_q2": None,
            "za_is_2d": None, "za_vacuum_axis": None,
            "za_vacuum_axis_in_primitive": None, "za_ok": None,
            "za_r2_q1": None, "za_r2_q2": None,
            "za_log_resid_rms_q1": None, "za_log_resid_rms_q2": None,
            "za_margin_q1": None, "za_margin_q2": None,
            "za_branch_index_q1": None, "za_branch_index_q2": None,
            "za_zfrac_q1": None, "za_zfrac_q2": None,
            "za_zfrac_min_q1": None, "za_zfrac_min_q2": None,
            "za_branch_note_q1": "", "za_branch_note_q2": "",
            "za_lowest_exponent_q1": None, "za_lowest_exponent_q2": None,
            "za_lowest_r2_q1": None, "za_lowest_r2_q2": None,
            "za_lowest_log_resid_rms_q1": None,
            "za_lowest_log_resid_rms_q2": None,
            "za_symprec": ZA_SYMPREC_DEFAULT,
            "za_symprec_relax": ZA_SYMPREC_RELAX,
            "za_spacegroup_default": None, "za_spacegroup_relax": None,
            "za_symmetry_relaxed": None,
            "za_exponent_default_q1": None, "za_exponent_default_q2": None,
            "za_needs_review": None,
            "za_margin_min": ZA_MARGIN_MIN, "za_r2_min": ZA_R2_MIN,
            "za_note": "ZA analysis failed: %s" % e,
        })
    (outdir / "phonon_band_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")
    def _fp(v):
        return "%.3f" % v if v is not None else "None"
    print("[DONE] phonon band: min %.4f THz, max %.4f THz%s | ZA p(q1,q2)=%s/%s "
          "(lowest-branch %s/%s, symprec=%s)"
          % (fmin, fmax, " (NAC)" if nac else "",
             _fp(summary["za_exponent_q1"]), _fp(summary["za_exponent_q2"]),
             _fp(summary["za_lowest_exponent_q1"]),
             _fp(summary["za_lowest_exponent_q2"]),
             summary.get("za_symprec")),
          flush=True)
    if summary.get("za_needs_review"):
        print("[WARN] ZA verdict flagged for human review: %s"
              % summary.get("za_note"), flush=True)
    _emit(summary, 0)


if __name__ == "__main__":
    main()
