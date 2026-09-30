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

q-mesh (MESH in step.conf):
  "n n n"  explicit mesh, used as is (one BTE run);
  auto     phonopy length convention n_i = max(1, round(L*|b_i|)) (|b_i| in 1/A,
           no 2pi) on the phono3py PRIMITIVE cell, symmetry-consistent via
           phono3py's own length->mesh; L = MESH_LENGTH;
  <float>  the same with that length.
  2D: the vacuum axis (the primitive axis with the shortest |b_i|) is pinned to 1.
Mesh convergence (MESH_CONV=auto, only for auto/<float> meshes): the BTE is
rerun on denser meshes (L *= MESH_CONV_FACTOR) until the whole kappa tensor at
MESH_CONV_T changes by < MESH_CONV_TOL_PCT between two consecutive meshes, or
MESH_CONV_MAX_LENGTH / MESH_CONV_MAX_POINTS is hit.  The denser mesh of the
converged pair is reported; mesh_convergence.json keeps every point.  With a
cutoff scan the scan runs at the starting mesh, then the chosen cutoff's fc is
converged in mesh.
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


def _voigt(kappa):
    """(n_T, 6) Voigt kappa [xx, yy, zz, yz, xz, xy] of the first sigma, or None."""
    k = np.asarray(kappa, float)
    if k.ndim >= 2 and k.shape[-1] == 6:
        return k.reshape(-1, k.shape[-2], 6)[0]
    return None


def _tensor(v):
    xx, yy, zz, yz, xz, xy = v
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], float)


def _summary(tc, cfg, source_fc, mesh=None):
    T = [float(x) for x in np.asarray(tc.temperatures, float).reshape(-1)]
    raw = [[float(v) for v in row] for row in _flat_kappa(tc.kappa)]
    n = min(len(T), len(raw))
    T, raw = T[:n], raw[:n]
    j = int(np.argmin(np.abs(np.asarray(T) - 300.0)))
    extra = {}
    voigt = _voigt(tc.kappa)
    if voigt is not None and len(voigt) >= n:
        voigt = voigt[:n]
        # The Cartesian xx/yy/zz follow the POSCAR orientation; for a cell that
        # is not in the standard setting (e.g. a rhombohedral primitive cell)
        # they are not the crystal axes.  The eigenvalues are.
        extra = {
            "kappa_voigt_xx_yy_zz_yz_xz_xy": [[float(x) for x in r] for r in voigt],
            "kappa_principal_300K": sorted(
                float(x) for x in np.linalg.eigvalsh(_tensor(voigt[j]))),
            "kappa_avg_300K": float(np.trace(_tensor(voigt[j])) / 3.0),
        }
    return {
        "KAPPA_DONE": True,
        "solver": "phono3py",
        "bte_method": str(cfg.get("bte_method") or "rta"),
        "mesh": mesh if mesh is not None else cfg["mesh"],
        "nac": bool(cfg.get("nac")),
        "isotope": bool(cfg.get("isotope", True)),
        "temperatures": T,
        "kappa_xx_yy_zz": raw,
        "kappa_inplane_xx_yy": [[v[0], v[1]] for v in raw],
        "kappa_300K_xx_yy_zz": raw[j],
        "kappa_inplane_300K": 0.5 * (raw[j][0] + raw[j][1]),
        "source_fc": source_fc,
        **extra,
    }


def _mesh_spec(cfg):
    """('explicit', [n,n,n]) | ('length', L) from cfg['mesh']."""
    m = str(cfg.get("mesh") or "auto").strip().lower()
    if m == "auto":
        return "length", float(cfg.get("mesh_length") or 50.0)
    parts = m.split()
    if len(parts) == 3:
        return "explicit", [int(x) for x in parts]
    if len(parts) == 1:
        return "length", float(parts[0])
    sys.exit("[ERROR] MESH must be auto, a length or three integers, got %r"
             % cfg.get("mesh"))


def _vacuum_axis(ph3):
    """Primitive axis with the shortest reciprocal vector (largest height)."""
    rec = np.linalg.inv(np.asarray(ph3.primitive.cell, float))  # columns = b_i
    return int(np.argmin(np.linalg.norm(rec, axis=0)))


def _apply_mesh(ph3, spec, cfg):
    """Set ph3.mesh_numbers from spec and return the mesh actually used."""
    kind, val = spec
    if kind == "explicit":
        ph3.mesh_numbers = list(val)
    elif cfg.get("is_2d") and cfg.get("mesh_2d_vacuum", True):
        # 2D: build the diagonal mesh ourselves (vacuum axis = 1).  Squashing the
        # scalar result instead is wrong when phono3py returns a generalized
        # (non-diagonal) grid: its D_diag entries are not per-axis counts.
        rec = np.linalg.norm(np.linalg.inv(np.asarray(ph3.primitive.cell, float)),
                             axis=0)
        m = [max(1, int(round(float(val) * b))) for b in rec]
        m[_vacuum_axis(ph3)] = 1
        try:
            ph3.mesh_numbers = m
        except Exception as e:  # noqa: BLE001
            sys.exit("[ERROR] 2D 网格 %s 不满足晶格对称性（%s）—— POSCAR 的真空轴大概"
                     "不垂直于面内基矢，先标准化晶胞" % (m, e))
    else:
        # 3D: phono3py turns a scalar into a symmetry-consistent grid on the
        # primitive cell (phonopy length convention).  For a basis whose axes
        # are coupled by symmetry (e.g. a tilted c) it returns a generalized
        # regular grid: D_diag like [1, 11, 22] is then NOT "points per axis";
        # the real grid is grid_matrix, the q-point count prod(D_diag).
        ph3.mesh_numbers = float(val)
    return [int(x) for x in np.asarray(ph3.mesh_numbers).reshape(-1)[:3]]


def _grid_info(ph3):
    """grid_matrix (3x3) + total q-point count of the grid in use."""
    try:
        gm = np.asarray(ph3.grid.grid_matrix).astype(int).tolist()
    except Exception:  # noqa: BLE001
        gm = None
    nq = int(np.prod(np.asarray(ph3.mesh_numbers).reshape(-1)[:3]))
    return gm, nq


def _run_bte(out, yaml_name, fc2, fc3, cfg, source_fc, write_kappa, spec=None):
    import phono3py
    ph3 = phono3py.load(str(out / yaml_name), produce_fc=False,
                        is_nac=bool(cfg.get("nac")), log_level=0)
    _load_fc(ph3, out, fc2, fc3, cfg.get("enable_fc"))
    mesh = _apply_mesh(ph3, spec or _mesh_spec(cfg), cfg)
    gm, nq = _grid_info(ph3)
    diag = gm is not None and all(gm[i][j] == 0 for i in range(3)
                                  for j in range(3) if i != j)
    print("[..] BTE mesh %s (%d q-points%s)"
          % (mesh, nq, "" if diag or gm is None
             else "; generalized grid, grid_matrix=%s" % gm), flush=True)
    ph3.init_phph_interaction()
    ph3.run_thermal_conductivity(
        is_LBTE=(str(cfg.get("bte_method") or "rta").lower() == "lbte"),
        temperatures=[float(x) for x in cfg["temperatures"]],
        is_isotope=bool(cfg.get("isotope", True)),
        write_kappa=write_kappa, log_level=1)
    s = _summary(ph3.thermal_conductivity, cfg, source_fc,
                 mesh=" ".join(str(x) for x in mesh))
    # "mesh" is phono3py's D_diag; for a generalized grid it is not points per
    # axis, so keep the full grid matrix and the q-point count alongside.
    s["mesh_grid_matrix"] = gm
    s["n_qpoints"] = nq
    return s


def _kappa_at(s, t):
    """Voigt kappa (or xx/yy/zz) at the temperature closest to t."""
    T = np.asarray(s["temperatures"], float)
    j = int(np.argmin(np.abs(T - float(t))))
    rows = s.get("kappa_voigt_xx_yy_zz_yz_xz_xy") or s["kappa_xx_yy_zz"]
    return np.asarray(rows[j], float), float(T[j])


def _mesh_change(prev, cur, t, mode="max_ii"):
    """Relative change of kappa between two runs at temperature t.

    mode='max_ii' (default): rel = max|d kappa_ij| / max|kappa_ii| (whole
    tensor).  Converges fast, but can mask a slowly-converging out-of-plane zz
    when it is much smaller than the in-plane kappa; the caller warns on that.
    mode='per_component': rel = max_i |d kappa_ii| / |kappa_ii| (each diagonal
    component normalised by *its own* value), so an anisotropic zz cannot be
    hidden by the larger in-plane one.

    Returns (rel, rel_vec, tt); rel_vec holds the per-component relative
    changes (None where the reference component is ~0).
    """
    a, _ = _kappa_at(prev, t)
    b, tt = _kappa_at(cur, t)
    n = min(len(a), len(b))
    nd = min(3, n)
    rel_vec = [abs(b[i] - a[i]) / abs(a[i]) if abs(a[i]) > 1e-12 else None
               for i in range(nd)]
    if str(mode).lower() == "max_ii":
        scale = float(np.max(np.abs(b[:3])))
        rel = None if scale <= 0 else float(np.max(np.abs(b - a)) / scale)
        return rel, rel_vec, tt
    vals = [x for x in rel_vec if x is not None]
    rel = max(vals) if vals else None
    # Off-diagonal components are still checked, but normalised by the largest
    # diagonal (their own values can be ~0, so a per-component ratio is not
    # meaningful there).
    if n > 3:
        scale = float(np.max(np.abs(b[:3])))
        if scale > 0:
            off = float(np.max(np.abs(b[3:n] - a[3:n])) / scale)
            rel = off if rel is None else max(rel, off)
    return rel, rel_vec, tt


def _converge_mesh(out, yaml_name, fc2, fc3, cfg, source_fc, first=None):
    """Densify the mesh until kappa stops changing.  Returns the summary of the
    last (densest converged) mesh, with a mesh_convergence report attached.

    first: optional (length, summary) already computed at the starting length
    (reused so the cutoff scan's run is not repeated)."""
    kind, L = _mesh_spec(cfg)
    conv = str(cfg.get("mesh_conv") or "off").lower() not in (
        "off", "false", "0", "no", "none", "")
    if kind == "explicit" or not conv:
        if first is not None:
            return first[1]
        return _run_bte(out, yaml_name, fc2, fc3, cfg, source_fc, True,
                        spec=(kind, L))
    tol = float(cfg.get("mesh_conv_tol_pct") or 3.0) / 100.0
    mode = str(cfg.get("mesh_conv_mode") or "max_ii").strip().lower()
    fac = float(cfg.get("mesh_conv_factor") or 1.25)
    lmax = float(cfg.get("mesh_conv_max_length") or 150.0)
    pmax = int(cfg.get("mesh_conv_max_points") or 64000)
    t_chk = float(cfg.get("mesh_conv_t") or 300.0)
    if fac <= 1.0:
        sys.exit("[ERROR] MESH_CONV_FACTOR must be > 1")

    import phono3py
    probe = phono3py.load(str(out / yaml_name), produce_fc=False, is_nac=False,
                          log_level=0)

    def _next_length(L, done_mesh):
        """Smallest L' >= L (in factor steps) whose mesh differs from done_mesh."""
        while L <= lmax + 1e-9:
            m = _apply_mesh(probe, ("length", L), cfg)
            if " ".join(map(str, m)) != done_mesh:
                return L, m
            L *= fac
        return None, None

    recs, prev, converged, stop = [], None, False, ""
    if first is not None:
        prev = first[1]
        recs.append({"length": round(float(first[0]), 3), "mesh": prev["mesh"],
                     "kappa_at_T": _kappa_at(prev, t_chk)[0].tolist(),
                     "rel_change": None})
        L = float(first[0]) * fac
    while True:
        L, m = _next_length(L, prev["mesh"] if prev else None)
        if L is None:
            stop = "MESH_CONV_MAX_LENGTH=%g reached" % lmax
            break
        if prev is not None and int(np.prod(m)) > pmax:
            stop = "MESH_CONV_MAX_POINTS=%d exceeded by mesh %s" % (pmax, m)
            break
        print("[..] 网格收敛：L=%.1f A -> mesh %s" % (L, m), flush=True)
        s = _run_bte(out, yaml_name, fc2, fc3, cfg, source_fc, True,
                     spec=("length", L))
        if prev is None:
            rel, rel_vec, tt = None, None, t_chk
        else:
            rel, rel_vec, tt = _mesh_change(prev, s, t_chk, mode)
        recs.append({"length": round(L, 3), "mesh": s["mesh"],
                     "kappa_at_T": _kappa_at(s, t_chk)[0].tolist(),
                     "rel_change": rel, "rel_change_per_component": rel_vec})
        print("    mesh %-10s kappa(%gK)=%s  change=%s%s"
              % (s["mesh"], tt, np.round(_kappa_at(s, t_chk)[0][:3], 4),
                 "-" if rel is None else "%.2f%%" % (100 * rel),
                 "" if rel_vec is None else " (per-component %s)"
                 % ["%.2f%%" % (100 * x) if x is not None else "-"
                    for x in rel_vec]), flush=True)
        prev = s
        if rel is not None and rel < tol:
            if str(mode).lower() == "max_ii" and rel_vec:
                _names = ("xx", "yy", "zz")
                over = [(_names[i], 100 * rel_vec[i])
                        for i in range(min(3, len(rel_vec)))
                        if rel_vec[i] is not None and rel_vec[i] > tol]
                if over:
                    print("[WARN] q 网格按 max_ii 已收敛，但分量 %s 的相对变化仍超阈（%s）；"
                          "该分量可能未收敛，见 mesh_convergence.json 的 rel_change_per_component"
                          % (", ".join("%s %.2f%%" % (n, v) for n, v in over),
                             "%.1f%%" % (100 * tol)), flush=True)
            converged = True
            break
        L *= fac
    if prev is None:
        sys.exit("[ERROR] 起始网格长度 MESH_LENGTH 已超过 MESH_CONV_MAX_LENGTH")
    rep = {"converged": converged, "tol_pct": 100 * tol, "check_T": t_chk,
           "factor": fac, "stop_reason": stop or "converged",
           "norm": "max_ii" if mode == "max_ii" else "per_component",
           "criterion": (
               "max|d kappa_ij| / max|kappa_ii| between consecutive meshes "
               "(whole tensor, first sigma)"
               if mode == "max_ii" else
               "per-component max_i |d kappa_ii| / |kappa_ii| between "
               "consecutive meshes (xx, yy, zz each normalised by itself)"),
           "records": recs}
    (out / "mesh_convergence.json").write_text(
        json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8",
        newline="\n")
    if not converged:
        print("[WARN] q 网格未收敛（%s）—— 结果取最密网格 %s，见 mesh_convergence.json；"
              "可调大 MESH_CONV_MAX_LENGTH/MESH_CONV_MAX_POINTS 或放宽 MESH_CONV_TOL_PCT"
              % (stop, prev["mesh"]), flush=True)
    else:
        print("[OK] q 网格收敛：mesh %s（相邻变化 %.2f%% < %.1f%%）"
              % (prev["mesh"], 100 * recs[-1]["rel_change"], 100 * tol), flush=True)
    prev["mesh_converged"] = converged
    prev["mesh_convergence"] = rep
    return prev


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
        # The scan ran at the starting mesh; now converge the mesh on the
        # chosen cutoff's force constants (the scan point is reused).
        s = _converge_mesh(out, yaml_name, "cut3_%s/fc2.hdf5" % tag,
                           "cut3_%s/fc3.hdf5" % tag, cfg,
                           "cutoff_scan/cut3_%s" % tag,
                           first=(_mesh_spec(cfg)[1], s))
        s["cutoff_selection"] = rep
        s["chosen_cutoff_A"] = chosen
        (out / "kappa_summary.json").write_text(
            json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
            newline="\n")
        print("[DONE] 选中 c3=%.2f A | mesh %s | 300K in-plane %.4f W/mK "
              "-> kappa_summary.json"
              % (chosen, s["mesh"], s["kappa_inplane_300K"]), flush=True)
        return

    # nominal: single cutoff
    s = _converge_mesh(out, yaml_name, "fc2.hdf5", "fc3.hdf5", cfg,
                       cfg.get("source_fc"))
    (out / "kappa_summary.json").write_text(
        json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
        newline="\n")
    print("[DONE] kappa_summary.json | mesh %s | 300K in-plane %.4f W/mK "
          "(xx=%.4f yy=%.4f zz=%.4f; principal %s)"
          % (s["mesh"], s["kappa_inplane_300K"], s["kappa_300K_xx_yy_zz"][0],
             s["kappa_300K_xx_yy_zz"][1], s["kappa_300K_xx_yy_zz"][2],
             s.get("kappa_principal_300K")),
          flush=True)


if __name__ == "__main__":
    main()
