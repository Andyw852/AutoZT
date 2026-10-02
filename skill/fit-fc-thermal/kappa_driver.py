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
cutoff scan the order is: (1) converge the mesh on a reference cutoff (the
nominal one, else the largest), (2) run every cutoff at that converged mesh and
pick the cutoff there, (3) every BTE run carries all temperatures, so the chosen
cutoff's result is the full kappa(T).  If the chosen cutoff differs from the
reference, its mesh is re-checked starting from the converged length.
"""
import json
import math
import shutil
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
        **_two_d_fields(cfg, raw, j),
    }


def _two_d_fields(cfg, raw, j):
    """2D: kappa_2d_normalized_* = raw kappa * h_perp/d (kl-dft-cpu convention).

    phono3py divides by the whole cell volume, vacuum included, so the raw
    in-plane kappa of a slab is diluted by h_perp/d; the zz component of a 2D
    layer has no physical meaning (no dispersion across the vacuum)."""
    g = cfg.get("kappa_2d_norm") or None
    if not g or not g.get("kappa_2d_norm_factor"):
        return {}
    f = float(g["kappa_2d_norm_factor"])
    out = {
        "kappa_2d_norm": g,
        "kappa_2d_normalized_xx_yy_zz": [[v * f for v in r] for r in raw],
        "kappa_2d_normalized_inplane_xx_yy": [[r[0] * f, r[1] * f] for r in raw],
        "kappa_2d_normalized_inplane_300K": 0.5 * (raw[j][0] + raw[j][1]) * f,
    }
    # 面热导 sheet conductance = κ·d = κ_raw·h⊥（W/K）：与层厚取法无关，跨文献可直接比
    # （Wu et al., arXiv:1607.06542）。h⊥ = V/A，即 phono3py 归一用的胞高。
    h = g.get("h_perp_A")
    if h:
        hm = float(h) * 1e-10
        out["sheet_conductance_W_per_K_xx_yy"] = [[r[0] * hm, r[1] * hm] for r in raw]
        out["sheet_conductance_inplane_300K_W_per_K"] = 0.5 * (raw[j][0] + raw[j][1]) * hm
    return out


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


def _mesh_change(prev, cur, t, mode="per_component", floor=0.1):
    """Relative change of kappa between two runs at temperature t.

    mode='per_component' (default): rel = max_i |d kappa_ii| / D_i with
    D_i = max(|kappa_ii|, floor * max_j |kappa_jj|): every diagonal component
    is judged against *its own* value, so a slowly converging out-of-plane zz
    of an anisotropic crystal cannot hide behind the larger in-plane kappa.
    The floor keeps a negligible component (a 2D vacuum axis, zz < 10 % of xx
    with floor=0.1) from blocking convergence on its own relative noise -- the
    reason the default was once switched to max_ii.  floor=0 is the strict
    per-component criterion.
    mode='max_ii': rel = max|d kappa_ij| / max|kappa_ii| (whole tensor, legacy).

    Returns (rel, rel_vec, tt); rel_vec holds the plain per-component relative
    changes |d kappa_ii| / |kappa_ii| (None where the reference is ~0).
    """
    a, _ = _kappa_at(prev, t)
    b, tt = _kappa_at(cur, t)
    n = min(len(a), len(b))
    nd = min(3, n)
    rel_vec = [abs(b[i] - a[i]) / abs(a[i]) if abs(a[i]) > 1e-12 else None
               for i in range(nd)]
    scale = float(np.max(np.abs(b[:nd]))) if nd else 0.0
    if str(mode).lower() == "max_ii":
        rel = None if scale <= 0 else float(np.max(np.abs(b[:n] - a[:n])) / scale)
        return rel, rel_vec, tt
    if scale <= 0:
        return None, rel_vec, tt
    lo = max(float(floor or 0.0), 0.0) * scale
    vals = []
    for i in range(nd):
        den = max(abs(a[i]), lo)
        if den > 1e-12:
            vals.append(abs(b[i] - a[i]) / den)
    rel = max(vals) if vals else None
    # Off-diagonal components are still checked, normalised by the largest
    # diagonal (their own values can be ~0).
    if n > 3:
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
    mode = str(cfg.get("mesh_conv_mode") or "per_component").strip().lower()
    if mode not in ("per_component", "max_ii"):
        sys.exit("[ERROR] MESH_CONV_MODE must be per_component or max_ii, got %r" % mode)
    floor = cfg.get("mesh_conv_floor")
    floor = 0.1 if floor is None else float(floor)
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

    recs, prev, converged, stop, unconv = [], None, False, "", []
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
            rel, rel_vec, tt = _mesh_change(prev, s, t_chk, mode, floor)
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
            if rel_vec:
                _names = ("xx", "yy", "zz")
                over = [(_names[i], 100 * rel_vec[i])
                        for i in range(min(3, len(rel_vec)))
                        if rel_vec[i] is not None and rel_vec[i] > tol]
                if over:
                    # max_ii: any slow component; per_component: only one below
                    # the floor.  Either way it is recorded, not just printed.
                    unconv = [n for n, _v in over]
                    print("[WARN] q 网格按 %s 已收敛，但分量 %s 自身的相对变化仍超阈（%s）；"
                          "该分量可能未收敛（kappa_summary.json 的 "
                          "mesh_convergence.components_over_tol）"
                          % (mode, ", ".join("%s %.2f%%" % (n, v) for n, v in over),
                             "%.1f%%" % (100 * tol)), flush=True)
            converged = True
            break
        L *= fac
    if prev is None:
        sys.exit("[ERROR] 起始网格长度 MESH_LENGTH 已超过 MESH_CONV_MAX_LENGTH")
    rep = {"converged": converged, "tol_pct": 100 * tol, "check_T": t_chk,
           "factor": fac, "stop_reason": stop or "converged",
           "norm": "max_ii" if mode == "max_ii" else "per_component",
           "floor": None if mode == "max_ii" else floor,
           "components_over_tol": unconv,
           "criterion": (
               "max|d kappa_ij| / max|kappa_ii| between consecutive meshes "
               "(whole tensor, first sigma)"
               if mode == "max_ii" else
               "per-component max_i |d kappa_ii| / max(|kappa_ii|, floor*max_j "
               "|kappa_jj|) between consecutive meshes (xx, yy, zz each "
               "normalised by itself; components below the floor by "
               "floor*max)"),
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


def _reference_tag(tags, nominal=None):
    """网格收敛用的参考截断：标称截断在扫描里就用它，否则取最大一档（包含的三阶
    相互作用最全，散射最强、收敛最慢，用它定的网格对其余档偏保守）。"""
    if nominal is not None:
        want = ("%.2f" % float(nominal)).replace(".", "p")
        if want in tags:
            return want
    return max(tags, key=_cut_of_tag)


def _cut_of_tag(tag):
    return float(tag.replace("p", "."))


def _select(out, cfg, per_cut):
    """Three-criteria cutoff choice from the per-cut kappa + cutoff_scan.json."""
    import cut3_select
    scan = json.loads((out / "cutoff_scan.json").read_text(encoding="utf-8"))
    rows = {round(float(r["cut"]), 3): r
            for r in (scan.get("records") or []) if r.get("cut") is not None}
    # pheasy 扫描 CUT3_BOOTSTRAP=0：没有重拟合样本，stable_upper_cut 恒为 None ——
    # 是"没测"，不是"不稳定"。旧的 cutoff_scan.json 没有 stability_measured 字段，
    # 按 bootstrap 次数推断，已有结果不用重拟合也能选截断。
    _boot = scan.get("bootstrap")
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
            "stability_measured": r.get("stability_measured",
                                        not (_boot is not None and int(_boot) == 0
                                             and r.get("stable_upper_cut") is None)),
            "kappa": k, "kappa_err": (abs(k) * rel if k is not None else None),
            "kappa_vec": (s.get("kappa_300K_xx_yy_zz") or None),
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


# 作图（期刊风格）：默认参考配色第一槽（dataviz palette.md），黑色细框、刻度朝内
_LINE = "#2a78d6"
_INK, _INK2 = "#000000", "#4d4d4d"
_SCI_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.8, "axes.edgecolor": _INK, "axes.labelcolor": _INK,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.minor.visible": False, "ytick.minor.visible": True,
    "xtick.color": _INK, "ytick.color": _INK,
    "axes.grid": False, "legend.frameon": False,
    "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42,
    "mathtext.default": "regular",
}


def _write_kappa_vs_cutoff(out, per_cut, rep=None, cfg=None):
    """kappa(300 K) and the fit error vs third-order cutoff.

    kappa_vs_cutoff.json always; kappa_vs_cutoff.png (300 dpi) + .pdf when
    matplotlib is available.  Left column: one panel per component (xx, yy,
    zz), every point labelled with its kappa and every segment with its
    relative change (bold inside the plateau tolerance).  Right column: the
    force RMSE of each cutoff's fit (cutoff_scan.json force_rmse_eV_per_A,
    written by S1 when FIT_RMSE_FRAMES > 0), or pheasy's relative error when
    that is all there is.  A 2D layer gets xx and yy only -- its zz has no
    physical meaning -- thickness-normalised when the run carries the factor.
    Written before the cutoff choice can fail."""
    cfg = cfg or {}
    scan = {}
    try:
        scan = json.loads((out / "cutoff_scan.json").read_text(encoding="utf-8"))
    except Exception:
        scan = {}
    shells = [float(x) for x in (scan.get("shells") or [])]
    srec = {round(float(r["cut"]), 3): r for r in (scan.get("records") or [])
            if r.get("cut") is not None}
    norm = cfg.get("kappa_2d_norm") or None
    f2d = float(norm["kappa_2d_norm_factor"]) if norm and norm.get(
        "kappa_2d_norm_factor") else None
    flags = {round(float(r["cut"]), 3): r for r in ((rep or {}).get("records") or [])}
    rows = []
    for tag, s in sorted(per_cut, key=lambda x: _cut_of_tag(x[0])):
        cut = _cut_of_tag(tag)
        k3 = list(s.get("kappa_300K_xx_yy_zz") or [None, None, None]) + [None] * 3
        f = flags.get(round(cut, 3), {})
        sr = srec.get(round(cut, 3), {})
        rm = sr.get("force_rmse_eV_per_A")
        te = sr.get("train_rel_err")
        row = {"cut_A": cut,
               "n_shells": (sum(1 for d in shells if d <= cut + 1e-9)
                            if shells else None),
               "kappa_xx": k3[0], "kappa_yy": k3[1], "kappa_zz": k3[2],
               "kappa_inplane": s.get("kappa_inplane_300K"),
               "mesh": s.get("mesh"),
               "force_rmse_meV_per_A": None if rm is None else 1e3 * float(rm),
               "force_rmse_relative": sr.get("force_rmse_relative"),
               "fit_rel_err": None if te is None else float(te),
               "stable_ok": f.get("stable_ok"), "plateau_ok": f.get("plateau_ok")}
        if f2d:
            for k in ("kappa_xx", "kappa_yy", "kappa_inplane"):
                row[k + "_2d_norm"] = None if row[k] is None else row[k] * f2d
        rows.append(row)
    if len(rows) == 1 and rows[0]["force_rmse_meV_per_A"] is None:
        # no cutoff scan: the nominal fit's own residual (S2 copies fc_fit_summary.json)
        try:
            fs = json.loads((out / "fc_fit_summary.json").read_text(encoding="utf-8"))
            if fs.get("fit_rmse_eV_per_A") is not None:
                rows[0]["force_rmse_meV_per_A"] = 1e3 * float(fs["fit_rmse_eV_per_A"])
                rows[0]["force_rmse_relative"] = fs.get("fit_rmse_relative")
            elif fs.get("pheasy_relative_error") is not None:
                rows[0]["fit_rel_err"] = float(fs["pheasy_relative_error"])
        except Exception:
            pass
    for a_, b_ in zip(rows, rows[1:]):
        b_["pct_change_from_prev"] = {
            k: (None if (a_[k] in (None, 0) or b_[k] is None
                         or abs(a_[k]) < 1e-12)
                else round(100.0 * (b_[k] - a_[k]) / abs(a_[k]), 3))
            for k in ("kappa_xx", "kappa_yy", "kappa_zz", "kappa_inplane")}
    big = max([abs(r[k] or 0.0) for r in rows for k in ("kappa_xx", "kappa_yy")] or [0])
    is_2d = bool(cfg.get("is_2d")) or (bool(rows) and all(
        abs(r["kappa_zz"] or 0.0) < 1e-3 * max(big, 1e-12) for r in rows))
    doc = {"T_K": 300, "rows": rows, "is_2d": is_2d,
           "method": cfg.get("method_label"),
           "kappa_2d_norm_factor": f2d,
           "chosen_cut_A": (rep or {}).get("chosen_cut"),
           "status": (rep or {}).get("status"),
           "note": "each cutoff at the scan's starting mesh; only the chosen one "
                   "is mesh-converged afterwards"}
    (out / "kappa_vs_cutoff.json").write_text(
        json.dumps(doc, indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")
    print("[..] kappa vs cutoff (300 K):", flush=True)
    for r in rows:
        print("     c3=%.2f A%s  in-plane %s%s  xx/yy/zz %s%s"
              % (r["cut_A"], "" if r["n_shells"] is None else " (%d shells)" % r["n_shells"],
                 r["kappa_inplane"],
                 "" if not f2d else " (2D-normalised %.3f)" % r["kappa_inplane_2d_norm"],
                 [r["kappa_xx"], r["kappa_yy"], r["kappa_zz"]],
                 "" if r["force_rmse_meV_per_A"] is None
                 else "  RMSE %.2f meV/A" % r["force_rmse_meV_per_A"]), flush=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator
    except Exception as e:                                 # noqa: BLE001
        print("[..] matplotlib unavailable (%s) -- kappa_vs_cutoff.json only" % e,
              flush=True)
        return doc
    tol = float((rep or {}).get("kappa_tol_pct") or 5.0)
    suf = "_2d_norm" if (is_2d and f2d) else ""
    xs = [r["cut_A"] for r in rows]
    ch = doc["chosen_cut_A"]
    if any(r["force_rmse_meV_per_A"] is not None for r in rows):
        ekey, elab = "force_rmse_meV_per_A", r"Force RMSE (meV $\AA^{-1}$)"
    elif any(r["fit_rel_err"] is not None for r in rows):
        ekey, elab = "fit_rel_err", "Fit relative error"
    else:
        ekey, elab = None, r"Force RMSE (meV $\AA^{-1}$)"
    xticks = (["%.2f\n(%d)" % (r["cut_A"], r["n_shells"]) for r in rows]
              if rows and rows[0]["n_shells"] is not None else None)
    with plt.rc_context(_SCI_RC):
        # always 2 x 2: (a) kxx | (b) kzz  /  (c) kyy | (d) force RMSE
        fig, grid = plt.subplots(2, 2, figsize=(7.0, 5.4), sharex=True,
                                 gridspec_kw={"wspace": 0.42, "hspace": 0.10})
        ax_xx, ax_zz = grid[0]
        ax_yy, ax_rm = grid[1]

        def _series(ax, key, fmt, label_pct):
            pts = [(x, r[key]) for x, r in zip(xs, rows) if r[key] is not None]
            ys = [p[1] for p in pts]
            if not ys:
                return False
            ax.plot([p[0] for p in pts], ys, color=_LINE, lw=1.2, zorder=3)
            ax.plot([p[0] for p in pts], ys, ls="none", marker="o", ms=4.5,
                    mfc=_LINE, mec=_INK, mew=0.6, zorder=4)
            lo, hi = min(ys), max(ys)
            # never zoom below 10 % of the value: a 0.01 % wiggle must not look
            # like a collapse
            span = max(hi - lo, 0.10 * max(abs(hi), abs(lo), 1e-12))
            mid = 0.5 * (hi + lo)
            ax.set_ylim(mid - 0.75 * span, mid + 0.75 * span)
            if fmt is None:      # enough decimals to tell neighbouring values apart
                diffs = [abs(b - a) for a, b in zip(ys, ys[1:]) if b != a]
                dmin = min(diffs) if diffs else abs(hi) or 1.0
                nd = max(1, min(4, int(math.ceil(-math.log10(dmin))) + 1)) \
                    if dmin > 0 else 2
                fmt = "%%.%df" % nd
            for k, (x, y) in enumerate(pts):
                # value under the point, on the side the curve does not run to
                nxt = pts[k + 1][1] if k + 1 < len(pts) else None
                prv = pts[k - 1][1] if k > 0 else None
                down = (nxt is not None and nxt < y) or (nxt is None and prv is not None
                                                         and prv > y)
                ax.annotate(fmt % y, (x, y), xytext=(-4 if down else 4, -5),
                            textcoords="offset points", ha="right" if down else "left",
                            va="top", fontsize=6.5, color=_INK)
            if label_pct:
                for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
                    if abs(y0) < 1e-12:
                        continue
                    pct = 100.0 * (y1 - y0) / abs(y0)
                    inside = abs(pct) <= tol
                    ax.annotate("%+.1f%%" % pct, ((x0 + x1) / 2, (y0 + y1) / 2),
                                xytext=(0, 7), textcoords="offset points",
                                ha="center", va="bottom", fontsize=6.5,
                                color=_INK if inside else _INK2,
                                fontweight="bold" if inside else "normal")
            return True

        _series(ax_xx, "kappa_xx" + suf, "%.1f", True)
        _series(ax_yy, "kappa_yy" + suf, "%.1f", True)
        if is_2d:
            ax_zz.set_yticks([])        # 2D: kappa_zz is not physical -> frame only
            ax_zz.yaxis.set_minor_locator(plt.NullLocator())
        else:
            _series(ax_zz, "kappa_zz", "%.1f", True)
        if not (ekey and _series(ax_rm, ekey, None, False)):
            ax_rm.set_yticks([])
            ax_rm.yaxis.set_minor_locator(plt.NullLocator())
        for ax, lab in ((ax_xx, "xx"), (ax_yy, "yy"), (ax_zz, "zz")):
            ax.set_ylabel(r"$\kappa_{%s}$ (W m$^{-1}$ K$^{-1}$)" % lab)
        ax_rm.set_ylabel(elab)
        for ax, let in ((ax_xx, "a"), (ax_zz, "b"), (ax_yy, "c"), (ax_rm, "d")):
            ax.text(0.03, 0.95, "(%s)" % let, transform=ax.transAxes, ha="left",
                    va="top", fontsize=9, fontweight="bold")
            if ax.get_yticks().size:
                ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            if ch is not None:
                ax.axvline(float(ch), color=_INK2, lw=0.7, ls=(0, (4, 2)), zorder=1)
        if len(xs) > 1:
            span = xs[-1] - xs[0]
            ax_xx.set_xlim(xs[0] - 0.22 * span, xs[-1] + 0.10 * span)
        else:                           # a single (nominal) cutoff: centre it
            ax_xx.set_xlim(xs[0] - 0.6, xs[0] + 0.6)
        for ax in (ax_yy, ax_rm):
            ax.set_xlabel(r"Third-order cutoff $r_c$ ($\AA$)")
            ax.set_xticks(xs)
            ax.set_xticklabels(xticks or ["%.2f" % x for x in xs])
        fig.subplots_adjust(top=0.93, bottom=0.13, left=0.10, right=0.98)
        if doc.get("method"):
            fig.text(0.10, 0.955, doc["method"], ha="left", va="bottom",
                     fontsize=10, fontweight="bold", transform=fig.transFigure)
        notes = ["T = 300 K; (n) = neighbour shells within $r_c$; "
                 "bold: |change| $\\leq$ %g%%" % tol]
        if len(rows) == 1:
            notes.insert(0, "single cutoff only (no cutoff scan): set CUT3_SCAN = auto "
                            "in step1_fit/step.conf to get kappa vs. cutoff")
        if ch is not None:
            notes.append("dashed: %s $r_c$ = %.2f $\\AA$ (%s)"
                         % ("chosen" if len(rows) > 1 else "nominal", float(ch),
                            doc["status"] or "nominal"))
        if is_2d:
            notes.append(("2D layer, d = %.2f $\\AA$; $\\kappa_{zz}$ not physical "
                          "(frame left empty)" % float(norm["thickness_d_A"])) if f2d else
                         "2D layer: raw cell-volume kappa (not thickness-normalised); "
                         "$\\kappa_{zz}$ not physical (frame left empty)")
        fig.text(0.0, -0.01, "\n".join(notes), ha="left", va="top", fontsize=6.5,
                 color=_INK2, transform=fig.transFigure)
        fig.savefig(str(out / "kappa_vs_cutoff.png"))
        fig.savefig(str(out / "kappa_vs_cutoff.pdf"))
        plt.close(fig)
    print("[OK] kappa_vs_cutoff.png / .pdf", flush=True)
    return doc


def compare_methods(out):
    """FIT_METHODS with several methods: kappa(300 K) and the fit's force RMSE
    per method -> methods_compare.json (data only), and every method's own 2x2
    kappa/RMSE figure copied to figures/kappa_vs_cutoff_<tag>.{png,pdf}."""
    out = Path(out)
    try:
        entries = json.loads((out / "methods.json").read_text(
            encoding="utf-8")).get("methods") or []
    except Exception:
        print("[..] no methods.json -- nothing to compare", flush=True)
        return None
    rows = []
    for e in entries:
        d = out / e.get("dir", ".")
        try:
            ks = json.loads((d / "kappa_summary.json").read_text(encoding="utf-8"))
        except Exception:
            ks = {}
        try:
            fs = json.loads((d / "fc_fit_summary.json").read_text(encoding="utf-8"))
        except Exception:
            fs = {}
        k2d = ks.get("kappa_2d_normalized_inplane_300K")
        rm = fs.get("fit_rmse_eV_per_A")
        rows.append({
            "tag": e.get("tag"), "engine": e.get("engine"), "method": e.get("method"),
            "primary": bool(e.get("primary")), "skipped": e.get("skipped"),
            "kappa_done": bool(ks.get("KAPPA_DONE")),
            "stable": fs.get("stable"),
            "chosen_cutoff_A": ks.get("chosen_cutoff_A"),
            "kappa_inplane_300K": ks.get("kappa_inplane_300K"),
            "kappa_2d_normalized_inplane_300K": k2d,
            "kappa_avg_300K": ks.get("kappa_avg_300K"),
            "kappa_300K_xx_yy_zz": ks.get("kappa_300K_xx_yy_zz"),
            "sheet_conductance_inplane_300K_W_per_K":
                ks.get("sheet_conductance_inplane_300K_W_per_K"),
            "fit_rmse_meV_per_A": None if rm is None else 1e3 * float(rm),
            "fit_rmse_relative": fs.get("fit_rmse_relative"),
            "pheasy_relative_error": fs.get("pheasy_relative_error"),
        })
    is_2d = any(r["kappa_2d_normalized_inplane_300K"] is not None for r in rows)
    if is_2d:
        kkey = "kappa_2d_normalized_inplane_300K"
    elif any(r["kappa_avg_300K"] is not None for r in rows):
        kkey = "kappa_avg_300K"            # trace/3: orientation independent
    else:
        kkey = "kappa_inplane_300K"
    doc = {"T_K": 300, "kappa_key": kkey, "rows": rows}
    (out / "methods_compare.json").write_text(
        json.dumps(doc, indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")
    print("[..] method comparison (300 K, %s):" % kkey, flush=True)
    for r in rows:
        print("     %s%-28s kappa=%-10s RMSE=%s meV/A  stable=%s rc=%s"
              % ("*" if r["primary"] else " ", r["tag"], r[kkey],
                 r["fit_rmse_meV_per_A"], r["stable"], r["chosen_cutoff_A"]),
              flush=True)
    # one 2x2 kappa/RMSE figure per method, collected side by side
    figs = out / "figures"
    n = 0
    for r, e in zip(rows, entries):
        d = out / e.get("dir", ".")
        for ext in ("png", "pdf"):
            f = d / ("kappa_vs_cutoff.%s" % ext)
            if f.is_file():
                figs.mkdir(exist_ok=True)
                shutil.copyfile(str(f), str(figs / ("kappa_vs_cutoff_%s.%s"
                                                    % (r["tag"], ext))))
                n += ext == "png"
    print("[OK] methods_compare.json; %d per-method figure(s) in figures/" % n,
          flush=True)
    return doc



def main():
    out = Path.cwd()
    cfg = json.loads((out / "kappa_config.json").read_text(encoding="utf-8"))
    yaml_name = cfg["disp_yaml"]

    if cfg.get("scan") and cfg.get("cut_scans"):
        # 收敛顺序（2026-10-02 用户定，文献通行做法）：
        #   ① q 网格：在参考截断上做网格收敛，得到收敛网格长度 L*；
        #   ② 截断：每档截断都在 L* 上跑 BTE，κ 随截断的平台/逐壳层稳定性在收敛网格上判
        #      （以前在起始粗网格上选截断，再只对选中档收敛网格——选截断用的 κ 本身没收敛，
        #       kappa_vs_cutoff 图上的数值也和最终报告的 κ 对不上）；
        #   ③ 温度：每次 BTE 都按 T_MIN..T_MAX 全部温度算（RTA 下多温度几乎不加成本），
        #      最终报告即选中截断 + 收敛网格的完整 κ(T)。
        # 选中档若不是参考档，再从 L* 起对它复核一次网格（至少多一档加密）。
        tags = list(cfg["cut_scans"])
        ref = _reference_tag(tags, cfg.get("nominal_cut3_A"))
        kind, _L0 = _mesh_spec(cfg)
        print("[..] ① q 网格收敛（参考截断 %s A）" % _cut_of_tag(ref), flush=True)
        s_ref = _converge_mesh(out, yaml_name, "cut3_%s/fc2.hdf5" % ref,
                               "cut3_%s/fc3.hdf5" % ref, cfg,
                               "cutoff_scan/cut3_%s" % ref)
        mesh_rep = s_ref.get("mesh_convergence")
        if mesh_rep:
            (out / "mesh_convergence_reference.json").write_text(
                json.dumps(dict(mesh_rep, reference_cut3_A=_cut_of_tag(ref)),
                           indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")
        spec = (("length", float(mesh_rep["records"][-1]["length"]))
                if (kind == "length" and mesh_rep and mesh_rep.get("records"))
                else None)
        print("[..] ② 截断扫描（收敛网格 %s）" % s_ref["mesh"], flush=True)
        per_cut = []
        for tag in tags:
            if tag == ref:
                s = s_ref
            else:
                print("[..] BTE @ cutoff %s" % _cut_of_tag(tag), flush=True)
                s = _run_bte(out, yaml_name, "cut3_%s/fc2.hdf5" % tag,
                             "cut3_%s/fc3.hdf5" % tag, cfg,
                             "cutoff_scan/cut3_%s" % tag, write_kappa=False,
                             spec=spec)
            _cd = out / ("cut3_%s" % tag)
            _cd.mkdir(exist_ok=True)
            (_cd / "kappa_summary.json").write_text(
                json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
                newline="\n")
            per_cut.append((tag, s))
        chosen, rep, by_tag = _select(out, cfg, per_cut)
        try:
            _write_kappa_vs_cutoff(out, per_cut, rep, cfg)
        except Exception as e:                             # noqa: BLE001
            print("[WARN] kappa_vs_cutoff not written: %s" % e, flush=True)
        if chosen is None:
            (out / "kappa_summary.json").write_text(
                json.dumps({"KAPPA_DONE": False, "cutoff_selection": rep},
                           indent=2, ensure_ascii=False), encoding="utf-8",
                newline="\n")
            sys.exit("[ERROR] 没有可用截断（数据量/稳定性不足）—— 见 cutoff_selection.json")
        tag = ("%.2f" % chosen).replace(".", "p")
        s = by_tag[tag]
        if tag != ref and spec is not None:
            print("[..] 选中截断 %.2f A ≠ 参考截断：从收敛网格起复核一次" % chosen, flush=True)
            s = _converge_mesh(out, yaml_name, "cut3_%s/fc2.hdf5" % tag,
                               "cut3_%s/fc3.hdf5" % tag, cfg,
                               "cutoff_scan/cut3_%s" % tag, first=(spec[1], s))
        elif mesh_rep:
            s["mesh_converged"] = s_ref.get("mesh_converged")
            s["mesh_convergence"] = mesh_rep
        s["convergence_order"] = ["q_mesh@reference_cut", "cutoff@converged_mesh",
                                  "temperatures"]
        s["mesh_reference_cut3_A"] = _cut_of_tag(ref)
        s["cutoff_selection"] = rep
        s["chosen_cutoff_A"] = chosen
        (out / "kappa_summary.json").write_text(
            json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
            newline="\n")
        print("[DONE] ③ 选中 c3=%.2f A | mesh %s | 300K in-plane %.4f W/mK "
              "| T=%s K -> kappa_summary.json"
              % (chosen, s["mesh"], s["kappa_inplane_300K"],
                 "/".join("%g" % t for t in s.get("temperatures") or [])), flush=True)
        return

    # nominal: single cutoff
    print("[WARN] 没有截断扫描（step1_fit 的 CUT3_SCAN=off，或只有一档候选）：只算标称截断，"
          "kappa_vs_cutoff 图只有一个点。要 κ 随截断/壳层的变化：\n"
          "       autozt -tt fit-fc-thermal -p <材料> -j step1_fit conf --set params.CUT3_SCAN=auto"
          "  然后 retry S1_fit、再跑 S2_kappa", flush=True)
    s = _converge_mesh(out, yaml_name, "fc2.hdf5", "fc3.hdf5", cfg,
                       cfg.get("source_fc"))
    (out / "kappa_summary.json").write_text(
        json.dumps(s, indent=2, ensure_ascii=False), encoding="utf-8",
        newline="\n")
    if cfg.get("nominal_cut3_A") is not None:
        # every fitting method gets its 2x2 figure, a single point without a scan
        c = float(cfg["nominal_cut3_A"])
        try:
            _write_kappa_vs_cutoff(out, [(("%.2f" % c).replace(".", "p"), s)],
                                   {"chosen_cut": c, "status": "nominal"}, cfg)
        except Exception as e:                             # noqa: BLE001
            print("[WARN] kappa_vs_cutoff not written: %s" % e, flush=True)
    print("[DONE] kappa_summary.json | mesh %s | 300K in-plane %.4f W/mK "
          "(xx=%.4f yy=%.4f zz=%.4f; principal %s)"
          % (s["mesh"], s["kappa_inplane_300K"], s["kappa_300K_xx_yy_zz"][0],
             s["kappa_300K_xx_yy_zz"][1], s["kappa_300K_xx_yy_zz"][2],
             s.get("kappa_principal_300K")),
          flush=True)
    if s.get("kappa_2d_normalized_inplane_300K") is not None:
        print("[DONE] 2D 层厚归一化后 300K 面内 %.4f W/mK（d=%.3f Å，factor=%.4f）；"
              "面热导 %.4g W/K（与层厚取法无关）；κzz 对 2D 无物理意义"
              % (s["kappa_2d_normalized_inplane_300K"],
                 s["kappa_2d_norm"]["thickness_d_A"],
                 s["kappa_2d_norm"]["kappa_2d_norm_factor"],
                 s.get("sheet_conductance_inplane_300K_W_per_K") or float("nan")),
              flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "compare":
        compare_methods(Path.cwd())
    else:
        main()
