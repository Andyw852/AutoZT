#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""defect-dft-cpu 核心（优先纯标准库，超算登录节点可跑）。

职责：POSCAR 解析/写出、超胞、通用缺陷枚举（对称性不等价位优先用 spglib，
缺失时尝试 pymatgen；没有可靠后端则拒绝枚举）、
INCAR/KPOINTS/POTCAR/submit.sh 渲染。

通用化：不再写死任何元素/结构类型；元素角色、缺陷类型、电荷态窗口、
SOC/IVDW/ISPIN/NCORE/KPAR 等全部由 step.conf 控制。
对称性枚举必须有 spglib 或 pymatgen；POSCAR IO 和 MIC 只需标准库。
"""
import os
import re
import math
from functools import lru_cache
from pathlib import Path

# ---------------------------------------------------------------- POSCAR IO
def parse_poscar(path):
    """Read VASP5 POSCAR geometry; selective flags are intentionally not retained.

    A single positive scale or negative target volume is supported. Three
    anisotropic scale factors and VASP4 (no species names) are rejected.
    """
    L = Path(path).read_text(encoding="utf-8-sig").splitlines()
    if len(L) < 8:
        raise ValueError("Incomplete POSCAR header")
    comment = L[0].strip()
    scales = L[1].split()
    if len(scales) != 1:
        raise ValueError("POSCAR requires one scalar scale factor")
    scale = float(scales[0])
    lat = [[float(x) for x in L[i].split()[:3]] for i in range(2, 5)]
    if (not math.isfinite(scale) or scale == 0 or
            any(len(row) != 3 or not all(math.isfinite(x) for x in row) for row in lat)):
        raise ValueError("Invalid POSCAR scale or lattice")
    a = lat
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
           - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
           + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    if not math.isfinite(det) or det == 0:
        raise ValueError("Singular POSCAR lattice")
    if scale < 0:
        scale = (-scale / abs(det)) ** (1.0 / 3.0)
    syms = L[5].split()
    cnts = [int(x) for x in L[6].split()]
    if (not syms or len(syms) != len(cnts) or any(n < 0 for n in cnts)
            or any(not re.fullmatch(r"[A-Z][a-z]?", s) for s in syms)):
        raise ValueError("POSCAR requires explicit species names and matching counts")
    mode_line = 7
    if L[mode_line].strip().lower().startswith("s"):
        mode_line += 1
    if mode_line >= len(L):
        raise ValueError("Missing POSCAR coordinate mode")
    mode = L[mode_line].strip().lower()
    if not mode.startswith(("d", "c", "k")):
        raise ValueError("Unknown POSCAR coordinate mode")
    start = mode_line + 1
    coords = [[float(x) for x in line.split()[:3]]
              for line in L[start:start + sum(cnts)]]
    if (len(coords) != sum(cnts) or
            any(len(c) != 3 or not all(math.isfinite(x) for x in c) for c in coords)):
        raise ValueError("Missing or invalid POSCAR atomic coordinates")
    atoms = [s for s, n in zip(syms, cnts) for _ in range(n)]
    if mode.startswith(("c", "k")):
        # Both lattice and Cartesian coordinates scale identically. Convert
        # against the unscaled ROW-vector lattice: fractional = cart @ inv(a).
        inv = [[0.0]*3 for _ in range(3)]
        inv[0][0] = (a[1][1]*a[2][2]-a[1][2]*a[2][1])/det
        inv[0][1] = (a[0][2]*a[2][1]-a[0][1]*a[2][2])/det
        inv[0][2] = (a[0][1]*a[1][2]-a[0][2]*a[1][1])/det
        inv[1][0] = (a[1][2]*a[2][0]-a[1][0]*a[2][2])/det
        inv[1][1] = (a[0][0]*a[2][2]-a[0][2]*a[2][0])/det
        inv[1][2] = (a[0][2]*a[1][0]-a[0][0]*a[1][2])/det
        inv[2][0] = (a[1][0]*a[2][1]-a[1][1]*a[2][0])/det
        inv[2][1] = (a[0][1]*a[2][0]-a[0][0]*a[2][1])/det
        inv[2][2] = (a[0][0]*a[1][1]-a[0][1]*a[1][0])/det
        coords = [[sum(c[j] * inv[j][i] for j in range(3)) for i in range(3)]
                  for c in coords]
    return {"comment": comment, "lat": [[x * scale for x in row] for row in lat],
            "atoms": atoms, "coords": [[x % 1 for x in c] for c in coords]}


def write_poscar(path, struct, comment=None):
    lat = struct["lat"]
    atoms = struct["atoms"]
    order = []
    for a in atoms:
        if a not in order:
            order.append(a)
    counts = [atoms.count(a) for a in order]
    idx = {a: i for i, a in enumerate(order)}
    lines = [comment or struct.get("comment", "defect"),
             "1.0",
             "  %16.10f  %16.10f  %16.10f" % tuple(lat[0]),
             "  %16.10f  %16.10f  %16.10f" % tuple(lat[1]),
             "  %16.10f  %16.10f  %16.10f" % tuple(lat[2]),
             "  " + "  ".join(order),
             "  " + "  ".join(str(counts[i]) for i in range(len(order))),
             "Direct"]
    grouped = [[] for _ in order]
    for a, c in zip(atoms, struct["coords"]):
        grouped[idx[a]].append(c)
    for g in grouped:
        for c in g:
            lines.append("  %16.10f  %16.10f  %16.10f" % tuple(c))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")

# ---------------------------------------------------------------- 超胞 / 缺陷
def supercell(struct, n=(3, 3, 1)):
    lat = struct["lat"]
    new_lat = [[lat[i][j] * n[i] for j in range(3)] for i in range(3)]
    atoms, coords = [], []
    for a, c in zip(struct["atoms"], struct["coords"]):
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    atoms.append(a)
                    coords.append([(c[0] + i) / n[0],
                                   (c[1] + j) / n[1],
                                   (c[2] + k) / n[2]])
    return {"lat": new_lat, "atoms": atoms, "coords": coords}

# ------------------------------------------------- 元素表（对称性 / 缺陷角色）
# 原子序数（spglib 需要）
_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
    "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
    "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
    "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
    "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94,
}

# 非金属/准非金属（作阴离子），其余按阳离子。启发式：Te 在碲化物里是阴离子、
# Sb/Bi 是阳离子，故 Sb/Bi 不在此表。step.conf 显式 CATIONS/ANIONS 永远优先。
_ANION_LIKE = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br",
               "Te", "I", "At", "He", "Ne", "Ar", "Kr", "Xe", "Rn"}


def _warn(msg):
    import sys
    sys.stderr.write("[defect] WARNING: %s\n" % msg)


def element_roles(conf, species):
    """返回 (cations, anions)。step.conf 的 CATIONS/ANIONS（逗号分隔）优先；
    缺省按 _ANION_LIKE 启发式；只给一边时另一边取补集。"""
    conf = conf or {}

    def _split(k):
        v = (conf.get(k, "") or "").strip()
        return [x.strip() for x in v.split(",") if x.strip()] if v else None

    cats, ans = _split("CATIONS"), _split("ANIONS")
    if cats is None and ans is None:
        ans = [a for a in species if a in _ANION_LIKE]
        cats = [a for a in species if a not in set(ans)]
        _warn("CATIONS/ANIONS 未配置，按内置表判定：阳离子=%s 阴离子=%s"
              "（不符请写入 step.conf 覆盖）" % (cats, ans))
    elif cats is None:
        cats = [a for a in species if a not in set(ans)]
    elif ans is None:
        ans = [a for a in species if a not in set(cats)]
    return cats, ans


def _zkey(c, tol=1e-3):
    z = c[2] % 1.0
    z = min(z, 1.0 - z)          # 旧的显示标签，不定义对称性等价关系
    return round(z / tol) * tol


def _legacy_orbits(struct):
    """旧 (元素,z) 启发式，仅保留私有接口；不可用于可靠缺陷枚举。"""
    seen, reps = {}, []
    for i, (a, c) in enumerate(zip(struct["atoms"], struct["coords"])):
        key = (a, _zkey(c))
        if key not in seen:
            seen[key] = i
            reps.append((i, key))
    return reps


def _symmetry_data(struct, symprec=1e-3, mode="auto"):
    """Reliable atom equivalence and fractional space-group operations.

    The legacy (species,z) heuristic is not a symmetry backend. Missing or
    failed backends fail closed instead of silently losing defect sites.
    """
    if mode not in ("auto", "spglib", "pymatgen"):
        raise ValueError("SYMMODE must be auto, spglib or pymatgen; legacy is unsafe")
    if not math.isfinite(symprec) or symprec <= 0:
        raise ValueError("SYMPREC must be finite and positive")
    lat, pos, atoms = struct["lat"], struct["coords"], struct["atoms"]
    failures = []
    for backend in (("spglib", "pymatgen") if mode == "auto" else (mode,)):
        try:
            if backend == "spglib":
                import spglib
                # Distinct positive species IDs suffice; no incomplete Z table.
                ids = {a: i + 1 for i, a in enumerate(sorted(set(atoms)))}
                ds = spglib.get_symmetry_dataset(
                    (lat, pos, [ids[a] for a in atoms]), symprec=symprec)
            else:
                from pymatgen.core import Lattice, Structure
                from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
                st = Structure(Lattice(lat), atoms, pos, coords_are_cartesian=False)
                ds = SpacegroupAnalyzer(st, symprec=symprec).get_symmetry_dataset()
            if ds is None:
                raise ValueError("no symmetry dataset returned")
            def field(key):
                return ds[key] if isinstance(ds, dict) else getattr(ds, key)
            eq = [int(x) for x in field("equivalent_atoms")]
            rotations = [[[int(x) for x in row] for row in rot]
                         for rot in field("rotations")]
            translations = [[float(x) for x in t] for t in field("translations")]
            if (len(eq) != len(atoms) or not rotations or
                    len(rotations) != len(translations) or
                    any(i < 0 or i >= len(atoms) or atoms[i] != atoms[k]
                        for k, i in enumerate(eq))):
                raise ValueError("invalid symmetry dataset")
            return eq, list(zip(rotations, translations)), backend
        except Exception as exc:
            failures.append("%s: %s" % (backend, exc))
    raise RuntimeError("Reliable symmetry backend required (install spglib or pymatgen); "
                       + "; ".join(failures))


def symmetry_orbits(struct, symprec=1e-3, mode="auto"):
    """Return (representative, species, display z, cell multiplicity).

    Only spglib/pymatgen equivalence defines orbits; z is a display label.
    """
    eq, _, _ = _symmetry_data(struct, symprec, mode)
    return [(i, struct["atoms"][i], _zkey(struct["coords"][i]), eq.count(i))
            for i in sorted(set(eq))]


def inequivalent_sites(struct, symprec=1e-3, mode="auto"):
    """兼容旧接口：返回 [(原子下标, (元素, z层)), ...]。"""
    return [(i, (a, z)) for i, a, z, m in symmetry_orbits(struct, symprec, mode)]

def vdw_gap_z(struct):
    """旧实现（已废弃）：找最大分数 z 间隙中点，会把间隙原子压到离宿主 ~1 Å。
    保留仅为兼容；间隙缺陷请用 find_voids。"""
    zs = sorted(_zkey(c) for c in struct["coords"])
    best, best_w = 0.5, 0.0
    for z1, z2 in zip(zs, zs[1:]):
        if z2 - z1 > best_w:
            best, best_w = 0.5 * (z1 + z2), z2 - z1
    return best


def _frac_to_cart(frac, lat):
    return [sum(frac[k] * lat[k][j] for k in range(3)) for j in range(3)]


def _metric(lat):
    """metric 张量 G = lat·lat^T，用于由分数坐标差直接求笛卡尔距离平方。"""
    return [[sum(lat[i][k] * lat[j][k] for k in range(3)) for j in range(3)]
            for i in range(3)]


@lru_cache(maxsize=32)
def _mic_factor(metric):
    """Upper Cholesky factor R: ||R df||² = df^T G df (cached per cell)."""
    if len(metric) != 9 or not all(math.isfinite(x) for x in metric):
        raise ValueError("MIC requires a finite 3x3 positive-definite metric")
    G = [metric[i:i + 3] for i in (0, 3, 6)]
    L = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(i + 1):
            value = G[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if value <= 0:
                    raise ValueError("MIC requires a nonsingular positive-definite cell")
                L[i][j] = math.sqrt(value)
            else:
                L[i][j] = value / L[j][j]
    return tuple(tuple(L[j][i] for j in range(3)) for i in range(3))


def _mic_distance(df, R):
    """Closest lattice vector by bounded 3D sphere enumeration.

    Wrapping components or testing only 27 images is not valid for arbitrary
    skew/unreduced cells. A Babai upper bound gives a finite sphere; triangular
    back substitution enumerates all images that can improve it. Severely
    ill-conditioned cells may be slow or lose floating-point precision.
    """
    if not all(math.isfinite(x) for x in df):
        raise ValueError("MIC requires finite fractional coordinates")
    d = [x - round(x) for x in df]
    residual = [0.0] * 3
    best = 0.0
    for i in (2, 1, 0):
        tail = sum(R[i][j] * residual[j] for j in range(i + 1, 3))
        center = d[i] + tail / R[i][i]
        residual[i] = d[i] - round(center)
        best += (R[i][i] * residual[i] + tail) ** 2

    def search(i, partial):
        nonlocal best
        if i < 0:
            best = min(best, partial)
            return
        tail = sum(R[i][j] * residual[j] for j in range(i + 1, 3))
        center = d[i] + tail / R[i][i]
        radius = math.sqrt(max(0.0, best - partial)) / R[i][i]
        # Outward padding protects images exactly on the current bound.
        lo = math.ceil(center - radius - 1e-12)
        hi = math.floor(center + radius + 1e-12)
        for image in range(lo, hi + 1):
            residual[i] = d[i] - image
            value = partial + (R[i][i] * residual[i] + tail) ** 2
            if value <= best + 1e-14 * max(1.0, best):
                search(i - 1, value)
    search(2, 0.0)
    return math.sqrt(best)


def _min_dist_to_atoms(frac_pos, host_frac, G):
    """Minimum Cartesian periodic distance to host atoms in a skew cell."""
    R = _mic_factor(tuple(x for row in G for x in row))
    return min((_mic_distance([frac_pos[k] - h[k] for k in range(3)], R)
                for h in host_frac), default=float("inf"))


def _frac_sep(f1, f2, G):
    """Minimum Cartesian periodic separation, valid for skew/unreduced cells."""
    R = _mic_factor(tuple(x for row in G for x in row))
    return _mic_distance([f1[k] - f2[k] for k in range(3)], R)


def interstitial_orbits(struct, seeds, symprec=1e-3, mode="auto"):
    """Deduplicate a fixed raw seed budget; never refill discarded seeds.

    site_count is the number of distinct images in the full host-cell orbit,
    NOT the number of sampled seeds. Grid seeds are not optimized void centers;
    their multiplicity need not equal that of a subsequently relaxed site.
    """
    _, operations, backend = _symmetry_data(struct, symprec, mode)
    G = _metric(struct["lat"])
    groups = []
    for seed_index, seed in enumerate(seeds):
        frac = tuple(float(x) % 1.0 for x in seed)
        for group in groups:
            if any(_frac_sep(frac, image, G) <= symprec for image in group["orbit"]):
                group["seed_indices"].append(seed_index)
                break
        else:
            orbit = []
            for rotation, translation in operations:
                image = tuple((sum(rotation[i][j] * frac[j] for j in range(3))
                               + translation[i]) % 1.0 for i in range(3))
                if not any(_frac_sep(image, old, G) <= symprec for old in orbit):
                    orbit.append(image)
            groups.append({"frac": frac, "orbit": orbit, "seed_indices": [seed_index]})
    return [(group["frac"], {
        "site_count": len(group["orbit"]),
        "multiplicity_basis": "host_cell_symmetry_orbit_of_grid_seed",
        "seed_indices": group["seed_indices"],
        "seed_count": len(group["seed_indices"]),
        "raw_seed_count": len(seeds),
        "symmetry_backend": backend,
        "symprec": symprec,
    }) for group in groups]


def find_voids(struct, n=2, grid=20, min_sep=2.0):
    """搜索最多 n 个原始网格种子，随后调用 interstitial_orbits 去重。
    不保证穷尽局域空隙极大值；去重后不可补附近网格点凑足 n 个。"""
    if n <= 0:
        return []
    if grid <= 0 or not math.isfinite(min_sep) or min_sep < 0:
        raise ValueError("VOID_GRID must be positive and VOID_MIN_SEP finite/nonnegative")
    lat = struct["lat"]
    G = _metric(lat)
    host = struct["coords"]  # 分数坐标
    cands = []
    for i in range(grid):
        for j in range(grid):
            for k in range(grid):
                f = ((i + 0.5) / grid, (j + 0.5) / grid, (k + 0.5) / grid)
                cands.append((_min_dist_to_atoms(f, host, G), f))
    cands.sort(key=lambda x: -x[0])
    chosen = []
    for d, f in cands:
        if all(_frac_sep(f, cf, G) > min_sep for _, cf in chosen):
            chosen.append((d, f))
            if len(chosen) >= n:
                break
    return [f for _, f in chosen]

def remove_atom(struct, idx):
    st = {"lat": [r[:] for r in struct["lat"]],
          "atoms": struct["atoms"][:], "coords": [c[:] for c in struct["coords"]]}
    st["atoms"].pop(idx); st["coords"].pop(idx)
    return st

def substitute_atom(struct, idx, el):
    st = {"lat": [r[:] for r in struct["lat"]],
          "atoms": struct["atoms"][:], "coords": [c[:] for c in struct["coords"]]}
    st["atoms"][idx] = el
    return st

def insert_atom(struct, el, frac):
    st = {"lat": [r[:] for r in struct["lat"]],
          "atoms": struct["atoms"][:] + [el],
          "coords": [c[:] for c in struct["coords"]] + [list(frac)]}
    return st

DEFECT_TYPES_ALL = ("vacancy", "antisite", "pair", "interstitial")


def enumerate_defects(sc, conf=None):
    """通用缺陷枚举。返回 [(dir_suffix, 显示名, 结构, meta)]。

    meta = {"type":..., "element":..., "sub":..., "site_count":...}
    step.conf 可配：
      DEFECT_TYPES   = vacancy,antisite,pair,interstitial  # 生成哪几类
      DEFECT_SPECIES = C,Mg          # 只生成"与这些元素相关"的缺陷：
                                     #   空位=宿主在该集合；反位=宿主或替换者任一在该集合；
                                     #   间隙=插入元素在该集合。缺省全部元素。
      ANTISITE_PAIRS = C-Mg,Mg_C     # 只生成这些替换对（缺省全部非同种）
      CATIONS/ANIONS = Mg / C        # 决定阳离子互占位对（pair）
      N_INTERSTITIALS= 2             # 原始种子预算；去重后可少于此数
      SYMPREC        = 1e-3
    不再写死任何元素：Mg4C60、A2B2Te5 等任意结构走同一套逻辑。
    """
    conf = conf or {}

    def _split(k, default=None):
        v = (conf.get(k, "") or "").strip()
        if not v:
            return default
        return [x.strip() for x in v.split(",") if x.strip()]

    sites = symmetry_orbits(sc, symprec=float(conf.get("SYMPREC", "1e-3") or 1e-3),
                            mode=(conf.get("SYMMODE", "auto") or "auto"))
    species = sorted({a for a in sc["atoms"]})
    want = set(_split("DEFECT_TYPES", list(DEFECT_TYPES_ALL)))
    only_el = set(_split("DEFECT_SPECIES", species) or species)
    pair_pairs = None
    ap = _split("ANTISITE_PAIRS")
    if ap:
        pair_pairs = set()
        for x in ap:
            parts = re.split(r"[-_>]+", x)
            if len(parts) == 2:
                pair_pairs.add((parts[0].strip(), parts[1].strip()))
    cats, ans = element_roles(conf, species)
    nvoid = int(conf.get("N_INTERSTITIALS", "2") or 2)
    out = []

    if "vacancy" in want:
        for idx, a, z, m in sites:
            if a not in only_el:
                continue
            out.append(("v_%s" % a, "v_%s(z=%.3f)" % (a, z), remove_atom(sc, idx),
                        {"type": "vacancy", "element": a, "site_count": m}))

    if "antisite" in want:
        for idx, a, z, m in sites:
            for b in species:
                if b == a:
                    continue
                if b not in only_el and a not in only_el:
                    continue
                if pair_pairs is not None and (b, a) not in pair_pairs:
                    continue
                out.append(("%s_%s" % (b, a), "%s_%s(z=%.3f)" % (b, a, z),
                            substitute_atom(sc, idx, b),
                            {"type": "antisite", "element": a, "sub": b,
                             "site_count": m}))

    if "pair" in want:
        pair_cats = [c for c in cats if c in set(sc["atoms"])]
        for ii in range(len(pair_cats)):
            for jj in range(ii + 1, len(pair_cats)):
                A, B = pair_cats[ii], pair_cats[jj]
                if A not in only_el and B not in only_el:
                    continue
                ia = next(i for i, x in enumerate(sc["atoms"]) if x == A)
                ib = next(i for i, x in enumerate(sc["atoms"]) if x == B)
                pair = substitute_atom(substitute_atom(sc, ia, B), ib, A)
                out.append(("%s_%s__%s_%s_pair" % (A, B, B, A),
                            "%s_%s+%s_%s pair" % (A, B, B, A), pair,
                            {"type": "pair", "element": A, "sub": B, "site_count": 2}))

    if "interstitial" in want and nvoid > 0:
        voids = find_voids(sc, n=nvoid, grid=int(conf.get("VOID_GRID", "20") or 20),
                           min_sep=float(conf.get("VOID_MIN_SEP", "2.0") or 2.0))
        void_orbits = interstitial_orbits(
            sc, voids, symprec=float(conf.get("SYMPREC", "1e-3") or 1e-3),
            mode=(conf.get("SYMMODE", "auto") or "auto"))
        tags = "abcdefghij"
        for el in species:
            if el not in only_el:
                continue
            for k, (frac, orbit_meta) in enumerate(void_orbits):
                tag = tags[k] if k < len(tags) else str(k)
                out.append(("%s_i_%s" % (el, tag), "%s_i@void%s" % (el, tag),
                            insert_atom(sc, el, frac),
                            {"type": "interstitial", "element": el, "void": tag,
                              "frac_coords": list(frac), **orbit_meta}))
    return out


def charge_states(name, conf=None, species=None):
    """缺陷电荷态清单（通用，不再写死价电子表）。

    CHARGE_MODE:
      window（默认，最保险）：返回 QMIN..QMAX（默认 -2..2）全部整数；
        带电稳定态交给 S4 形成能分析挑选 -> 不会漏态。
      valence（省算力）：按阴/阳离子角色给单边区间，等价旧 A2B2Te5 口径。
    金属体系可 QMIN=0 QMAX=0（金属里带电态意义有限）。
    """
    conf = conf or {}
    mode = (conf.get("CHARGE_MODE", "window") or "window").strip().lower()
    if mode == "valence":
        return _charge_states_valence(name, conf, species)
    qmin = int(conf.get("QMIN", "-2") or -2)
    qmax = int(conf.get("QMAX", "2") or 2)
    return list(range(qmin, qmax + 1))


def _charge_states_valence(name, conf, species):
    conf = conf or {}
    species = species or sorted(set(re.findall(r"[A-Z][a-z]?", name)))
    cats, ans = element_roles(conf, species)
    Q = abs(int(conf.get("QMAX_MAG", "2") or 2))
    plus = [0] + list(range(1, Q + 1))
    minus = [0] + list(range(-1, -Q - 1, -1))
    if "pair" in name:
        return [0]
    if "_i_" in name:
        el = name.split("_i_")[0]
        return minus if el in ans else plus
    if name.startswith("v_"):
        el = name[2:]
        return plus if el in ans else minus
    parts = name.split("_q")[0].split("_")
    if len(parts) >= 2:
        X, Y = parts[0], parts[1]
        if X in ans and Y in cats:
            return plus
        if X in cats and Y in ans:
            return minus
    return [0]

# ---------------------------------------------------------------- 输入渲染
def render_incar(tpl, mapping):
    text = Path(tpl).read_text(encoding="utf-8")
    for k, v in mapping.items():
        text = text.replace("{{" + k + "}}", str(v))
    return text

def write_kpoints(path, mesh=(2, 2, 2)):
    text = ("Auto mesh (supercell)\n0\nGamma\n%d %d %d\n0 0 0\n" % mesh)
    Path(path).write_text(text, encoding="utf-8")

def potcar_variant(el):
    # 本体系推荐的 PAW：主族重元素用 _d（含 d 半芯态），Te/Sb/Bi 默认
    return {"Pb": "Pb_d", "Sn": "Sn_d", "Sb": "Sb", "Bi": "Bi_d", "Te": "Te"}.get(el, el)

def assemble_potcar(symbols, potcar_dir, variants=None, out_path="POTCAR"):
    """symbols 按 POSCAR 元素顺序；potcar_dir 下含 potpaw_PBE/POTCAR_<variant>。
    直接写到 out_path（避免先写 cwd 再移动导致的竞态）。"""
    potcar_dir = os.path.expanduser(str(potcar_dir))
    variants = variants or {}
    blocks = []
    for el in symbols:
        v = variants.get(el, potcar_variant(el))
        cands = [
            os.path.join(potcar_dir, "potpaw_PBE_64", v, "POTCAR"),   # jzzn 子目录布局
            os.path.join(potcar_dir, "potpaw_PBE_54", v, "POTCAR"),
            os.path.join(potcar_dir, "potpaw_PBE", v, "POTCAR"),      # 通用子目录
            os.path.join(potcar_dir, "potpaw_PBE", "POTCAR_" + v),    # 平铺文件
            os.path.join(potcar_dir, "potpaw_PBE.54", "POTCAR_" + v),
        ]
        p = next((c for c in cands if os.path.exists(c)), None)
        if p is None:
            raise SystemExit("[错误] 找不到 POTCAR(%s)：试过 %s（请核对 step.conf 的 POTCAR_DIR / variant）" % (v, potcar_dir))
        blocks.append(Path(p).read_text(encoding="utf-8"))
    # 完整拼接：每个 POTCAR 都以 TITEL 行开头，VASP 按块顺序读取，无需删行
    content = "".join(b if b.endswith("\n") else b + "\n" for b in blocks)
    Path(out_path).write_text(content, encoding="utf-8")

def render_submit(tpl_path, out_path, jobname):
    text = Path(tpl_path).read_text(encoding="utf-8")
    text = text.replace("{{JOBNAME}}", jobname)
    left = set(re.findall(r"\{\{(\w+)\}\}", text))
    if left:
        raise SystemExit("[错误] submit 模板仍有未填充占位符：%s" % left)
    Path(out_path).write_text(text, encoding="utf-8")

# ---------------------------------------------------------------- step.conf + 作业组装
def load_stepconf(path="step.conf"):
    import configparser
    c = configparser.ConfigParser(inline_comment_prefixes=("#",))
    c.optionxform = str          # 保留键名大小写（SUPERCELL 等）
    c.read(path)
    return {k: v.strip() for k, v in c.items("params")}

def parse_variant(conf):
    out = {}
    for item in conf.get("POTCAR_VARIANT", "").split(","):
        item = item.strip()
        if ":" in item:
            k, v = item.split(":", 1)
            out[k.strip()] = v.strip()
    return out

def find_incar_tpl():
    """优先 cwd（tf 从技能 templates/ 推到材料目录），其次技能目录 templates/。"""
    for d in (Path.cwd(), Path(__file__).resolve().parent / "templates"):
        p = d / "incar_defect.tpl"
        if p.exists():
            return str(p)
    raise SystemExit("[错误] 找不到 incar_defect.tpl")


def find_submit_tpl(use_ncl=True, dim="3d"):
    """优先 cwd（tf 从 setting/<hpc>/templates 推来），其次技能目录 templates/。"""
    base = "submit_ncl" if use_ncl else "submit_std"
    names = ["%s_%s.tpl" % (base, dim), "%s.tpl" % base]
    cwd = Path.cwd()
    skill = Path(__file__).resolve().parent / "templates"
    for d in (cwd, skill):
        for n in names:
            p = d / n
            if p.exists():
                return str(p)
    raise SystemExit("[错误] 找不到提交模板 %s（请在 setting/<hpc>/templates 或技能 templates/ 放置）" % names[0])

def sanitize_label(text):
    import re
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip()).strip("_.-") or "defect"

def build_job(outdir, struct, conf, incar_map, jobname, use_ncl=None):
    """在 outdir 下写 POSCAR/INCAR/KPOINTS/POTCAR/submit.sh。struct 需已排好元素序。

    use_ncl=None（默认）时按 conf 的 SOC 自动选：SOC=1 -> submit_ncl_3d（vasp_ncl），
    SOC=0 -> submit_std_3d（vasp_std）。显式传 True/False 可覆盖。"""
    os.makedirs(outdir, exist_ok=True)
    # 元素顺序（保持 POSCAR 一致）
    order = []
    for a in struct["atoms"]:
        if a not in order:
            order.append(a)
    write_poscar(os.path.join(outdir, "POSCAR"), struct, comment=incar_map.get("SYSTEM", "defect"))
    natoms = len(struct["atoms"])
    gga = "PE"
    vdw = "IVDW = 12" if conf.get("FUNC", "") == "pbe-d3" else ""
    soc = str(conf.get("SOC", "1")).strip().lower() not in ("0", "false", "no", "off")
    if use_ncl is None:
        use_ncl = soc          # SOC 关闭 -> vasp_std；否则 vasp_ncl
    m = dict(incar_map)
    m.setdefault("ENCUT", conf.get("ENCUT", "370"))
    m.setdefault("GGA", gga)
    m.setdefault("VDW_LINE", vdw)
    m.setdefault("IBRION", "-1"); m.setdefault("ISIF", "0"); m.setdefault("NSW", "0")
    m.setdefault("EDIFFG_LINE", "")
    m.setdefault("LVHAR_LINE", "")
    m.setdefault("ICHARG_LINE", "ICHARG = 2")
    m.setdefault("SIGMA_LINE", "SIGMA = 0.05")
    ncomp = 3 if soc else 1          # 非共线 3 分量 / 共线 1 分量
    m.setdefault("MAGMOM", "%d*0" % (ncomp * natoms))
    m.setdefault("LSORBIT_LINE", "LSORBIT = .TRUE." if soc else "")
    m.setdefault("GGA_COMPAT_LINE", "GGA_COMPAT = .FALSE." if soc else "")
    # 共线 + 非零初始磁矩 -> 必须 ISPIN=2 才生效（非共线 SOC 时 ISPIN 恒为 1）
    mag = str(m.get("MAGMOM", "")).strip()
    is_zero_seed = bool(re.match(r"^\d+\*0(\.0)?$", mag))
    if soc or is_zero_seed:
        m.setdefault("ISPIN", conf.get("ISPIN", "1"))
    else:
        # SOC 关闭但给了非零初始磁矩：必须 ISPIN=2，否则 MAGMOM 被 VASP 忽略
        m.setdefault("ISPIN", "2")
    m.setdefault("NCORE", conf.get("NCORE", "6"))
    m.setdefault("KPAR", conf.get("KPAR", "4"))
    m.setdefault("NELECT_LINE", "")
    incar = render_incar(find_incar_tpl(), m)
    Path(os.path.join(outdir, "INCAR")).write_text(incar, encoding="utf-8")
    write_kpoints(os.path.join(outdir, "KPOINTS"), tuple(int(x) for x in conf.get("KMESH", "2 2 1").split()))
    assemble_potcar(order, conf["POTCAR_DIR"], parse_variant(conf),
                    os.path.join(outdir, "POTCAR"))
    tpl = find_submit_tpl(use_ncl)
    render_submit(tpl, os.path.join(outdir, "submit.sh"), sanitize_label(jobname))
    return order

def spin_seed_magmom(natoms, seed_atom=0, moment=1.0, ncomp=3):
    """SOC(非共线)下给 seed_atom 一个初始磁矩，打破自旋对称。

    奇数电子数的带电态（q=±1,±3…）有未配对自旋（doublet），MAGMOM=0 会把
    它强行算成非磁（闭壳配对）→ 能量偏高。这里 seed 一个 1 μB 的初始矩让 SCF
    探索自旋极化分支（矩会弛豫到缺陷局域）。偶数电子态不受影响（矩自动归零）。
    返回 MAGMOM 字符串（每原子 mx my mz）。"""
    parts = []
    if ncomp >= 3:
        for i in range(natoms):
            parts.append("%.1f 0 0" % moment if i == seed_atom else "0 0 0")
    else:
        for i in range(natoms):
            parts.append("%.1f" % moment if i == seed_atom else "0.0")
    return " ".join(parts)

def read_nelect_from_potcar(path):
    """从 POTCAR 汇总价电子数（每个 POTCAR 块的 ZVAL 之和）。"""
    zvals = [float(x) for x in re.findall(r"ZVAL\s*=\s*([\d.]+)", Path(path).read_text(encoding="utf-8"))]
    # ZVAL 每个元素一个；POTCAR 串联后每个块一个 ZVAL
    # 直接读 'number of ions per type' 不可靠，这里返回每个 ZVAL 列表供调用方按 counts 加权
    return zvals


def guarded_sbatch(workdir, jobname, submit_file="submit.sh", user=None,
                   lock_timeout=120, vis_window=60):
    """带锁防重复 sbatch（纯标准库，登录节点可跑）。

    与 autozt 远端提交守卫等价：flock 目录锁串行化，锁内实时 squeue 按 jobname
    匹配所有活动态，查询失败 fail closed；再用持久回执 + sacct 兜住 sbatch
    可见性延迟。返回 (ok, msg)：ok=True 表示已提交或已存在（幂等，调用方继续），
    ok=False 表示查询/提交失败，调用方应告警（不重试不改状态）。
    """
    import fcntl, time, getpass, subprocess
    user = user or os.environ.get("USER") or getpass.getuser()
    workdir = os.path.abspath(workdir)
    rcpt_path = os.path.join(workdir, ".tf_job_receipt")
    lock_path = os.path.join(workdir, ".tf_submit.lock")
    TERMINAL = {"COMPLETED", "CANCELLED", "FAILED", "TIMEOUT", "NODE_FAIL",
                "BOOT_FAIL", "OUT_OF_MEMORY", "DEADLINE", "PREEMPTED",
                "CD", "CA", "F", "TO", "NF", "OOM"}

    def _read_receipt():
        try:
            with open(rcpt_path, encoding="utf-8") as f:
                parts = f.read().strip().split()
            return (parts[0], float(parts[1])) if len(parts) >= 2 else None
        except (OSError, ValueError):
            return None

    def _write_receipt(jid):
        try:
            tmp = rcpt_path + ".tmp"
            with open(tmp, "w", encoding="utf-8") as f:
                f.write("%s %.6f\n" % (jid, time.time()))
            os.replace(tmp, rcpt_path)
        except OSError:
            pass

    def _sacct_state(jid):
        try:
            p = subprocess.run(["sacct", "-j", jid, "-n", "-P", "-X", "-o", "JobID,State"],
                               capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.TimeoutExpired):
            return None
        if p.returncode != 0:
            return None
        for ln in (p.stdout or "").splitlines():
            f = ln.split("|")
            if f and f[0].strip() == jid:
                return (f[1].strip() if len(f) > 1 else "")
        return ""

    fd = open(lock_path, "a+")
    try:
        got = False
        if lock_timeout <= 0:
            fcntl.flock(fd, fcntl.LOCK_EX)
            got = True
        else:
            deadline = time.time() + lock_timeout
            while True:
                try:
                    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    got = True
                    break
                except OSError:
                    if time.time() >= deadline:
                        break
                    time.sleep(0.2)
        if not got:
            return False, "提交锁获取超时（%ss）" % lock_timeout
        try:
            p = subprocess.run(["squeue", "-u", user, "-h", "-o", "%j"],
                               capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.TimeoutExpired):
            return False, "squeue 查询失败，拒绝提交（fail closed）"
        if p.returncode != 0:
            return False, "squeue 查询失败，拒绝提交（fail closed）"
        if jobname and jobname in set((p.stdout or "").split()):
            return True, "已在队列（jobname=%s）" % jobname
        rcpt = _read_receipt()
        if rcpt:
            jid, ts = rcpt
            if time.time() - ts < vis_window:
                st = _sacct_state(jid)
                if st is None:
                    return False, "存在近期回执 %s 但 sacct 查询失败，拒绝重复提交" % jid
                if st == "" or st.upper() not in TERMINAL:
                    return True, "已有近期提交回执 %s（%s），跳过" % (jid, st or "尚不可见")
        try:
            p = subprocess.run(["sbatch", submit_file], cwd=workdir,
                               capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.TimeoutExpired) as e:
            return False, "sbatch 失败：" + str(e)
        if p.returncode != 0:
            return False, "sbatch 失败：" + (p.stdout + p.stderr).strip()
        m = re.search(r"Submitted batch job\s+(\d+)", p.stdout or "")
        if m:
            _write_receipt(m.group(1))
        return True, (p.stdout or "").strip() or "已提交"
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
            fd.close()
        except Exception:
            pass