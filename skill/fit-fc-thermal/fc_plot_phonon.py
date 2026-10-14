#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fc_plot_phonon.py -- phonon band plot for the fit-fc-thermal skill (S2_plot).

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
OUTDIR = "phonon_band_plot"      # == skill.yaml S2_plot 的步骤目录名（done_marker 在这里找）
LOWF_MAX = 10.0
FULL_PAD = 3.0
import phonon_stability  # 2D ZA 二次性 + 虚频编排唯一真源（_common/phonon_stability.py）


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

    outdir = root / OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)

    # 2026-10-14：2D 的声子谱走 kz=0 的二维路径（za_2d.band_path_2d，与 kl-dft-cpu S5.1 同一套）；
    #   原来一律 auto_band_structure（seekpath 是 3D 工具），2D 图里多出 Γ-A 这类沿真空方向的线段。
    try:
        _plot2d, _plot_axis, _ = _dim_and_axis(src)
    except Exception as e:  # noqa: BLE001 —— 判不出维度按 3D 画，并说明
        print("[WARN] could not resolve the dimension for the band path (%s); using the 3D path" % e)
        _plot2d, _plot_axis = False, None
    band_path = {"kind": None}

    def _draw_band(fc2_path, dst):
        # band structure + the two figures for one set of fc2; |-> (fmin, fmax)
        dst.mkdir(parents=True, exist_ok=True)
        ph.force_constants = _read_fc2(fc2_path)
        try:
            if _plot2d:
                import za_2d
                _lat = za_2d.write_band_yaml_2d(ph, dst / "band-dft-cpu.yaml", npoints=51,
                                                cart_axis=_plot_axis)
                band_path["kind"] = "2d-%s" % _lat
            else:
                ph.auto_band_structure(plot=False, write_yaml=True, npoints=51,
                                       filename=str(dst / "band-dft-cpu.yaml"))
                band_path["kind"] = "seekpath"
            bd = ph.get_band_structure_dict()
        except Exception as e:
            _emit({"status": "error",
                   "reason": "phonopy band structure failed: %s" % e}, 40)
        fr = np.asarray(bd["frequencies"], dtype=object)
        ds = np.asarray(bd["distances"], dtype=object)
        tk = [(float(ds[i][0]), str(l).replace("$", ""))
              for i, l in enumerate(bd.get("labels") or []) if l]
        af = np.concatenate([np.asarray(f, float).ravel() for f in fr])
        lo0, hi0 = float(af.min()), float(af.max())

        def dr(fname, lo, hi):
            fig, ax = plt.subplots(figsize=(6.4, 4.8))
            for q, f in zip(ds, fr):
                # phonopy returns per-segment (n_qpoints, n_branches) arrays
                q = np.asarray(q, float)
                f = np.asarray(f, float)
                for ib in range(f.shape[1]):
                    ax.plot(q, f[:, ib], "-", lw=1.0, color="#1f4e79")
            ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
            for x, _lb in tk:
                ax.axvline(x, color="0.7", lw=0.8)
            if len(tk) > 1:
                ax.set_xticks([t[0] for t in tk])
                ax.set_xticklabels([t[1] for t in tk])
            ax.set_xlim(float(ds[0][0]), float(ds[-1][-1]))
            ax.set_ylim(lo, hi)
            ax.set_ylabel("Frequency (THz)")
            ax.set_title("Phonon dispersion%s" % (" (NAC)" if nac else ""),
                         fontsize=11)
            fig.tight_layout()
            fig.savefig(str(dst / fname), dpi=200)
            plt.close(fig)

        dr("phonon_band_full.png", lo0 - FULL_PAD, hi0 + FULL_PAD)
        dr("phonon_band_lowfreq.png", min(lo0 - FULL_PAD, -1.0), LOWF_MAX)
        return lo0, hi0, max(len(s) for s in ds)

    fmin, fmax, _nq = _draw_band(fc2p, outdir)

    # Per-cut phonon bands: every fc3 cutoff candidate written by the S1_fit scan
    # (cutoff_scan/cut3_<c>/fc2.hdf5) gets its own dispersion under
    # phonon_band_plot/cut3_<c>/, so a cutoff sweep always ships the spectra.
    cut_bands = []
    for _cd in sorted((src / "cutoff_scan").glob("cut3_*")):
        if not (_cd / "fc2.hdf5").is_file():
            continue
        try:
            _lo, _hi, _nq2 = _draw_band(_cd / "fc2.hdf5", outdir / _cd.name)
            cut_bands.append({"cut": float(_cd.name[5:].replace("p", ".")),
                              "dir": _cd.name, "min_frequency_THz": _lo,
                              "max_frequency_THz": _hi})
        except Exception as _e:  # noqa: BLE001
            print("[WARN] per-cut band %s failed: %s" % (_cd.name, _e))
    if cut_bands:
        print("[OK] per-cut phonon bands: %d" % len(cut_bands), flush=True)

    summary = {
        "status": "ok",
        "min_frequency_THz": fmin,
        "max_frequency_THz": fmax,
        "nac": nac,
        "n_qpoints": int(_nq),
        "figures": ["phonon_band_full.png", "phonon_band_lowfreq.png"],
        "band_yaml": str((src / "band-dft-cpu.yaml").relative_to(root)),
        "stable_here": bool(fmin >= -0.10),
        "cut3_bands": cut_bands,
        "band_path": band_path["kind"],
    }
    # ZA bending-branch exponent(s).  Strictly a 2D criterion; the analysis is the
    # single shared implementation in _common/phonon_stability.py (the same code
    # kl-dft-cpu / kl-mlff / phonon-dft use).  Failure is recorded, never raised,
    # so the plot still writes its marker.  3D (e.g. Si) -> za stays None.
    za = None
    try:
        _is2d, _cart_axis, _dim = _dim_and_axis(src)
        za = phonon_stability.za_check({}, ph, _cart_axis, _is2d)
    except Exception as e:
        print("[WARN] ZA exponent analysis failed (%s)" % e)
    if za is not None:
        summary["za"] = za
    (outdir / "phonon_band_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")
    def _fp(v):
        return "%.3f" % v if v is not None else "None"
    _p = (za or {}).get("p") or []
    _pl = (za or {}).get("p_lowest") or []
    print("[DONE] phonon band: min %.4f THz, max %.4f THz%s | ZA p(q1,q2)=%s/%s "
          "(lowest-branch %s/%s, symprec=%s)"
          % (fmin, fmax, " (NAC)" if nac else "",
             _fp(_p[0] if len(_p) > 0 else None),
             _fp(_p[1] if len(_p) > 1 else None),
             _fp(_pl[0] if len(_pl) > 0 else None),
             _fp(_pl[1] if len(_pl) > 1 else None),
             (za or {}).get("za_symprec")),
          flush=True)
    if za and za.get("za_needs_review"):
        print("[WARN] ZA verdict flagged for human review: %s"
              % za.get("note"), flush=True)
    _emit(summary, 0)


if __name__ == "__main__":
    main()
