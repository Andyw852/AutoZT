# -*- coding: utf-8 -*-
"""fc_common.py -- shared helpers for the fit-fc-thermal skill.

Deliberately tiny and dependency-free (standard library only) so the login-node
gen script runs under a plain system python.  The heavy lifting (dataset
parsing, fitting, gates) lives in fc_fit_driver.py, which runs on the compute
node inside the fitting job.

Placed next to the gen script by gen_need (the skill directory wins over the
common pool in tf's asset lookup chain).
"""
import math
import re
import sys
from pathlib import Path

# Cross-step bookkeeping files written by the phonon / thermal-conductivity
# skills; read for inheritance only (dimension, supercell, method).
METHOD_FILE = "workflow_method.txt"
KL_PARAMS = "kl_params.txt"


def engine_dir():
    """Directory of the shared fitting/kappa engine (fc_fit_driver.py,
    kappa_driver.py, cut3_select.py).  fc_common.py lives in that same
    directory, so its own parent is correct both in autozt's flat step dir and
    in the source tree (_common/fcfit/)."""
    return Path(__file__).resolve().parent


# ==========================================================================
# Simple KEY = VALUE files
# ==========================================================================
def read_keyval(path):
    """Read a KEY=VALUE file into an upper-cased dict (missing file -> {})."""
    d = {}
    p = Path(path)
    if not p.is_file():
        return d
    for ln in p.read_text(errors="ignore").splitlines():
        s = ln.strip()
        if s and not s.startswith("#") and "=" in s:
            k, v = s.split("=", 1)
            d[k.strip().upper()] = v.strip()
    return d


# ==========================================================================
# Structures (no numpy / ase on the login node)
# ==========================================================================
def poscar_lattice_lengths(path):
    """(|a1|, |a2|, |a3|) of a VASP POSCAR, scale factor applied.

    Stdlib-only: the gen script may run under a bare system python, and the
    only thing it needs from a structure is the edge-length ratio between
    POSCAR and SPOSCAR (to recover the supercell repetitions when the dataset
    ships no phonopy YAML)."""
    lines = Path(path).read_text(errors="ignore").splitlines()
    if len(lines) < 5:
        return None
    try:
        scale = [float(x) for x in lines[1].split()]
        lat = [[float(x) for x in lines[2 + i].split()[:3]] for i in range(3)]
    except (ValueError, IndexError):
        return None
    if len(scale) == 3:                      # VASP allows a per-axis scale
        lat = [[lat[i][k] * scale[k] for k in range(3)] for i in range(3)]
    elif scale:
        lat = [[x * scale[0] for x in v] for v in lat]
    return [math.sqrt(sum(c * c for c in v)) for v in lat]


def supercell_reps(uc_path, sc_path, tol=1e-3):
    """Diagonal supercell repetitions from the POSCAR/SPOSCAR edge ratio."""
    ru = poscar_lattice_lengths(uc_path)
    rs = poscar_lattice_lengths(sc_path)
    if not ru or not rs:
        return None
    reps = []
    for i in range(3):
        if ru[i] <= 0:
            return None
        ratio = rs[i] / ru[i]
        n = int(round(ratio))
        if n < 1 or abs(ratio - n) > tol:
            return None
        reps.append(n)
    return reps


# ==========================================================================
# ShengBTE FORCE_CONSTANTS_2ND atom order
# --------------------------------------------------------------------------
# ShengBTE (Src/input.f90 read2fc/split_index) reads the 2ND file as the full
# supercell and *assumes* phonopy's layout: index-1 = ((ia*nz + iz)*ny + iy)*nx + ix
# (x fastest, unit-cell atom slowest).  hiphive/phonopy writers emit the rows in
# whatever order the supercell they were given has -- a dataset SPOSCAR built by
# ASE repeat() interleaves the images, and ShengBTE then silently reads scrambled
# force constants (Mn2In2Se5 / MnIn2Se4, 2026-10).  shengbte_fc2_ensure_order()
# decides which layout the file is actually in from the force constants
# themselves (a correct labelling is invariant under lattice translations, a
# scrambled one is not), so it is idempotent: run it after every export, on
# re-used files, on files written by older versions -- it only ever reorders a
# file that is in SPOSCAR order, and refuses one that matches neither layout.
# Stdlib only (runs on the login node too).
# ==========================================================================
def _frac_key(f, n_digits=3):
    return tuple(int(round(((x % 1.0) % 1.0) * 10 ** n_digits)) % 10 ** n_digits for x in f)


class _FracIndex(object):
    """fractional position -> row, periodic, tolerant (dict first, scan fallback)."""

    def __init__(self, frac, tol=1e-3):
        self.frac, self.tol = frac, tol
        self.key = {_frac_key(f): k for k, f in enumerate(frac)}

    def get(self, f):
        k = self.key.get(_frac_key(f))
        if k is not None:
            return k
        for j, g in enumerate(self.frac):
            if all(abs(((f[d] - g[d]) + 0.5) % 1.0 - 0.5) < self.tol for d in range(3)):
                return j
        return None


def _matmul3(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def shengbte_dims(uc_lat, sc_lat, tol=1e-3):
    """Diagonal supercell repetitions (ShengBTE scell) or raise ValueError."""
    m = _matmul3(sc_lat, _inv3(uc_lat))
    mi = [[int(round(x)) for x in r] for r in m]
    if any(abs(m[i][j] - mi[i][j]) > tol for i in range(3) for j in range(3)):
        raise ValueError("supercell is not an integer multiple of the unit cell")
    if any(mi[i][j] for i in range(3) for j in range(3) if i != j):
        raise ValueError("supercell matrix %s is not diagonal -- ShengBTE scell cannot "
                         "express it" % mi)
    return [mi[0][0], mi[1][1], mi[2][2]]


def _read_fc2_blocks(path):
    lines = Path(path).read_text(errors="ignore").splitlines()
    n = int(lines[0].split()[0])
    if len(lines[0].split()) > 1 and int(lines[0].split()[1]) != n:
        raise ValueError("compact FORCE_CONSTANTS_2ND (%s) -- ShengBTE needs the full "
                         "supercell array" % lines[0].strip())
    blocks, i = {}, 1
    while i + 3 < len(lines) + 1 and len(blocks) < n * n:
        ab = lines[i].split()
        if len(ab) < 2:
            i += 1
            continue
        a, b = int(ab[0]) - 1, int(ab[1]) - 1
        blocks[(a, b)] = lines[i + 1:i + 4]
        i += 4
    if len(blocks) != n * n:
        raise ValueError("FORCE_CONSTANTS_2ND has %d blocks, expected %d" % (len(blocks), n * n))
    return n, blocks


def _block_vals(lines):
    return [float(x) for ln in lines for x in ln.split()[:3]]


def _translation_mismatch(blocks, frac, n, dims):
    """Relative mismatch of fc2 under one lattice translation, for a labelling
    index -> fractional supercell position.  ~0 for the true labelling."""
    where = _FracIndex(frac)
    axis = max(range(3), key=lambda i: dims[i])
    step = [0.0, 0.0, 0.0]
    step[axis] = 1.0 / dims[axis]
    # 只比非对角块（on-site 自作用块在同种原子间几乎一样，区分不了顺序）
    num = den = 0.0
    picks = list(range(0, n, max(1, n // 40)))
    for i in picks:
        ti = where.get([frac[i][d] + step[d] for d in range(3)])
        if ti is None:
            return float("inf")
        for j in range(n):
            if j == i:
                continue
            tj = where.get([frac[j][d] + step[d] for d in range(3)])
            if tj is None:
                return float("inf")
            a, b = _block_vals(blocks[(i, j)]), _block_vals(blocks[(ti, tj)])
            num += sum(abs(x - y) for x, y in zip(a, b))
            den += sum(abs(x) for x in a)
    return num / den if den > 0 else 0.0


def shengbte_fc2_ensure_order(fc2_path, uc_lat, uc_frac, sc_lat, sc_frac, tol=1e-3):
    """Make FORCE_CONSTANTS_2ND follow ShengBTE's supercell order (idempotent).

    uc_* = the unit cell ShengBTE's CONTROL describes, sc_* = the supercell the
    file was written from (e.g. SPOSCAR).  Returns a dict (dims, natoms, n,
    action = kept | reordered, mismatch_*); raises ValueError when the set
    cannot be run by ShengBTE."""
    dims = shengbte_dims(uc_lat, sc_lat)
    nat = len(uc_frac)
    n, blocks = _read_fc2_blocks(fc2_path)
    ncell = dims[0] * dims[1] * dims[2]
    if n != nat * ncell:
        raise ValueError("FORCE_CONSTANTS_2ND has %d atoms, expected natoms*prod(scell) "
                         "= %d*%d = %d" % (n, nat, ncell, nat * ncell))
    if len(sc_frac) != n:
        raise ValueError("supercell has %d atoms but FORCE_CONSTANTS_2ND %d"
                         % (len(sc_frac), n))
    sheng = []                                    # ShengBTE label k -> sc fractional
    for ia in range(nat):
        for iz in range(dims[2]):
            for iy in range(dims[1]):
                for ix in range(dims[0]):
                    sheng.append([(uc_frac[ia][0] + ix) / dims[0],
                                  (uc_frac[ia][1] + iy) / dims[1],
                                  (uc_frac[ia][2] + iz) / dims[2]])
    where = _FracIndex(sc_frac)
    perm = [where.get(f) for f in sheng]               # ShengBTE k -> sc row
    if None in perm or len(set(perm)) != n:
        raise ValueError("supercell atoms do not map onto unit cell + scell translations")
    info = {"dims": dims, "natoms": nat, "n": n}
    if perm == list(range(n)):
        info["action"] = "kept"                   # both layouts coincide
        return info
    m_sheng = _translation_mismatch(blocks, sheng, n, dims)
    m_sc = _translation_mismatch(blocks, sc_frac, n, dims)
    info.update({"mismatch_shengbte_order": m_sheng, "mismatch_supercell_order": m_sc})
    # 正确的标注下平移前后一致（拟合力常数 ~1e-10），错误的标注明显不一致（实测 ~1e-2）；
    # 两种都不一致时宁可报错也不猜
    if m_sheng <= 1e-4 and m_sheng <= m_sc:
        info["action"] = "kept"                   # already in ShengBTE order
        return info
    if m_sc > 1e-4:
        raise ValueError("FORCE_CONSTANTS_2ND matches neither the ShengBTE nor the "
                         "supercell atom order (translation mismatch %.2g / %.2g)"
                         % (m_sheng, m_sc))
    with open(str(fc2_path), "w") as f:
        f.write("%d %d\n" % (n, n))
        for a in range(n):
            for b in range(n):
                f.write("%d %d\n" % (a + 1, b + 1))
                f.write("\n".join(blocks[(perm[a], perm[b])]) + "\n")
    info["action"] = "reordered"
    return info


def shengbte_fc2_ensure_order_files(fc2_path, unit_poscar, super_poscar):
    """shengbte_fc2_ensure_order with the two cells read from POSCAR files."""
    ul, uf = read_poscar_cell_frac(unit_poscar)
    sl, sf = read_poscar_cell_frac(super_poscar)
    return shengbte_fc2_ensure_order(fc2_path, ul, uf, sl, sf)


# ==========================================================================
# Neighbour shells -> third-order cutoff candidates
# --------------------------------------------------------------------------
# Ported from kl-dft-cpu/kl_common.py (2026-09-24 user 流程): the third-order
# cutoff is taken at the MIDPOINT between two adjacent neighbour shells (never
# on a shell distance), and the fit then decides how far the data can actually
# determine it.  Stdlib-only here so the login-node gen can run under a bare
# system python (ASE is used for the enumeration when it is importable).
# ==========================================================================
def _det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


def _inv3(m):
    """3x3 matrix inverse (Cramer).  Rows are the lattice vectors."""
    d = _det3(m)
    if abs(d) < 1e-12:
        raise ValueError("singular lattice (zero volume)")
    c = [[0.0] * 3 for _ in range(3)]
    c[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d
    c[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / d
    c[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d
    c[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / d
    c[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d
    c[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / d
    c[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d
    c[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / d
    c[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d
    return c


def read_poscar_cell_frac(path):
    """(lattice 3x3, fractional positions Nx3) of a VASP POSCAR, scale applied.

    Handles the VASP-4 (counts only) and VASP-5 (symbols + counts) headers, the
    optional Selective-dynamics line and both Direct / Cartesian coordinates.
    """
    lines = Path(path).read_text(errors="ignore").splitlines()
    if len(lines) < 8:
        raise ValueError("POSCAR too short: %s" % path)
    scal = [float(x) for x in lines[1].split()]
    lat = [[float(x) for x in lines[2 + i].split()[:3]] for i in range(3)]
    if len(scal) == 3:
        lat = [[lat[i][k] * scal[k] for k in range(3)] for i in range(3)]
    elif scal:
        lat = [[x * scal[0] for x in v] for v in lat]

    toks5 = lines[5].split()
    if toks5 and all(re.fullmatch(r"[+-]?\d+", t) for t in toks5):
        counts = [int(t) for t in toks5]          # VASP-4: line 5 is the counts
        cur = 6
    else:
        counts = [int(t) for t in lines[6].split()]
        cur = 7
    s = lines[cur].strip().lower() if cur < len(lines) else ""
    if s.startswith("s"):                          # Selective dynamics
        cur += 1
        s = lines[cur].strip().lower() if cur < len(lines) else ""
    direct = not s.startswith(("c", "k"))          # Direct / Cartesian
    cur += 1
    n = sum(counts)
    frac = []
    for i in range(n):
        p = [float(x) for x in lines[cur + i].split()[:3]]
        if not direct:
            cart = p
            inv = _inv3(lat)
            p = [sum(cart[k] * inv[k][m] for k in range(3)) for m in range(3)]
        frac.append(p)
    return lat, frac


def _min_image_shifts(cell, r_max):
    """Periodic image shifts (-k..k per axis) needed to cover r_max."""
    import itertools
    vol = abs(_det3(cell))
    ks = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        area = _norm(_cross(cell[j], cell[k]))
        perp = vol / area if area > 1e-12 else r_max
        ks.append(max(1, int(math.ceil(r_max / max(perp, 1e-6)))))
    out = [(0, 0, 0)]
    for s in itertools.product(*[range(-k, k + 1) for k in ks]):
        if any(s):
            out.append(tuple(s))
    return out


# Quasi-continuous distance spectrum (low symmetry / large cells, e.g. a 128-atom
# fullerene network: pair distances 1.37 -> 7 A follow each other with no shell
# gap): the tol chaining merges several A of distances into one "shell", and the
# midpoints of shell means then land among many pair distances.  A shell wider
# than QUASI_SPAN flags it; auto then uses the explicit radius grid QUASI_GRID
# and lets the kappa plateau pick (2026-10-15).  Crystals have spans ~0 (Si
# 0.00, Mn2In2Se5 0.09 A), so 0.5 A leaves a wide margin.
QUASI_SPAN = 0.5
QUASI_GRID = (3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
QUASI_TAG = "[quasi-continuous]"
# Ordinary crystals: the auto candidates are the ShengBTE thirdorder.py -n
# N-th neighbour cutoffs (2026-10-15), see nth_neighbor_cutoffs.  NN_PAD is the
# extra enumeration radius beyond CUT3_MAX so the last u_{N+1} is listed.
NN_TOL = 1e-4
NN_PAD = 3.0


def _per_atom_shells(lat, frac, r_max, tol):
    """Per atom, the distinct neighbour distances u_1 < u_2 < ... below r_max
    (distances within tol are one shell)."""
    natom = len(frac)
    cart_of = [[sum(f[k] * lat[k][m] for k in range(3)) for m in range(3)]
               for f in frac]
    shifts = _min_image_shifts(lat, r_max)
    tvec = [[sum(s[k] * lat[k][m] for k in range(3)) for m in range(3)]
            for s in shifts]
    per_atom = []
    for i in range(natom):
        ds = []
        for j in range(natom):
            for t in tvec:
                d = math.sqrt(sum((cart_of[j][m] + t[m] - cart_of[i][m]) ** 2
                                  for m in range(3)))
                if 1e-6 < d < r_max:
                    ds.append(d)
        ds.sort()
        u = []
        for d in ds:
            if not u or d - u[-1] > tol:
                u.append(d)
        per_atom.append(u)
    return per_atom


def _snap_cut(c, pair_ds):
    """Round a -n cutoff to 2 decimals (dir name cut3_<%.2f>) without crossing
    a pair distance: round / floor / ceil are tried in turn and the first that
    keeps the same set of pairs as the exact value wins; None when the
    neighbouring distances are < 0.01 A apart."""
    import bisect
    lo = bisect.bisect_left(pair_ds, c)
    for x in (round(c, 2), math.floor(c * 100) / 100.0, math.ceil(c * 100) / 100.0):
        if bisect.bisect_left(pair_ds, x) == lo and (lo == len(pair_ds) or pair_ds[lo] != x):
            return round(x, 2)
    return None


def _shell_groups(poscar, r_max=7.0, tol=0.05, max_atoms=600):
    """Cluster the interatomic distances below r_max into shells.

    Returns (list of distances per shell, n_atoms).  Distances within tol merge into one
    shell: WS2's 8th/9th shells differ by only 0.02 A and must not split, or a
    candidate cutoff would fall in the 6.34-6.36 gap.
    """
    ds, natom = [], 0
    try:
        from ase.io import read as _ase_read
        from ase.neighborlist import neighbor_list as _nl
        atoms = _ase_read(str(poscar))
        if 0 < len(atoms) <= max_atoms:
            _i, _j, _d = _nl("ijd", atoms, r_max)
            ds = [float(x) for x in _d]
            natom = len(atoms)
    except Exception:                              # noqa: BLE001
        ds = []
    if not ds:
        lat, frac = read_poscar_cell_frac(poscar)
        natom = len(frac)
        if natom > max_atoms:
            raise RuntimeError("cell has %d atoms > %d, shell enumeration too big"
                               % (natom, max_atoms))
        shifts = _min_image_shifts(lat, r_max)
        for i in range(natom):
            ri = frac[i]
            for j in range(natom):
                for s in shifts:
                    p = [frac[j][k] + s[k] - ri[k] for k in range(3)]
                    cart = [sum(p[k] * lat[k][m] for k in range(3)) for m in range(3)]
                    d = math.sqrt(sum(x * x for x in cart))
                    if 1e-6 < d < r_max:
                        ds.append(d)
    ds.sort()
    shells = []
    for x in ds:
        if not shells or x - shells[-1][-1] > tol:
            shells.append([x])
        else:
            shells[-1].append(x)
    return shells, natom


def neighbor_shells(poscar, r_max=7.0, tol=0.05, max_atoms=600):
    """Shell mean distances (shell_distances, n_atoms); see _shell_groups."""
    groups, natom = _shell_groups(poscar, r_max=r_max, tol=tol, max_atoms=max_atoms)
    return [round(sum(s) / len(s), 4) for s in groups], natom


def cut3_candidates(poscar, r_max=7.0, tol=0.05, min_shells=2, min_gap=0.05,
                    max_candidates=12):
    """Candidate third-order cutoffs = the thirdorder.py -n N-th neighbour
    cutoffs for N = min_shells ... up to r_max.

    min_shells=2 starts at the 2nd neighbours.  Cutoffs are kept to 2 decimals
    with the same pair set as the exact -n value (_snap_cut); one within
    min_gap/2 of a pair distance is kept and only flagged in the note (nearly
    degenerate, a small distortion crosses it).  A quasi-continuous spectrum (a tol-clustered shell wider than
    QUASI_SPAN) has no meaningful N either: the candidates become QUASI_GRID
    and the note carries QUASI_TAG.  The caller drops candidates beyond the
    supercell safe cutoff.
    Returns (candidates, shells, note); shells are still the tol-clustered
    shell means (used to label the kappa-vs-cutoff figure).
    """
    groups, natom = _shell_groups(poscar, r_max=r_max, tol=tol)
    shells = [round(sum(s) / len(s), 4) for s in groups]
    wide = max(groups, key=lambda s: s[-1] - s[0]) if groups else None
    if wide is not None and wide[-1] - wide[0] > QUASI_SPAN:
        cands = [c for c in QUASI_GRID if c < r_max]
        note = ("%s distance spectrum: %.2f-%.2f A merged into one shell (span "
                "%.2f A > %.1f, %d pairs), shell midpoints are meaningless -> "
                "explicit radius grid %s, picked by the kappa plateau; set "
                "CUT3_CANDIDATES explicitly for another grid" % (
                    QUASI_TAG, wide[0], wide[-1], wide[-1] - wide[0], QUASI_SPAN,
                    len(wide), cands))
        print("[WARN] " + note)
        return cands, shells, note
    lat, frac = read_poscar_cell_frac(poscar)
    per_atom = _per_atom_shells(lat, frac, float(r_max) + NN_PAD, NN_TOL)
    nmax = min(len(u) for u in per_atom) - 1 if per_atom else 0
    pair_ds = sorted(d for u in per_atom for d in u)
    cands, tags, near, lost = [], [], [], []
    for n in range(max(1, int(min_shells)), nmax + 1):
        c_exact = max(0.5 * (u[n - 1] + u[n]) for u in per_atom)
        if c_exact > float(r_max):
            break
        c = _snap_cut(c_exact, pair_ds)
        if c is None:
            lost.append("N=%d:%.4f" % (n, c_exact))
            continue
        if any(abs(c - d) < 0.5 * float(min_gap) for d in pair_ds):
            near.append("N=%d:%.2f" % (n, c))
        if not cands or c > cands[-1] + 1e-9:
            cands.append(c)
            tags.append("N=%d:%.2f" % (n, c))
    note = "thirdorder -n N-th neighbour cutoffs (N>=%d, <= %.2f A): %s" % (
        int(min_shells), float(r_max), " ".join(tags) or "none")
    if near:
        note += "; note: %s within %.3f A of a pair distance (nearly degenerate)" % (
            " ".join(near), 0.5 * float(min_gap))
    if lost:
        note += "; %s: pair distances < 0.01 A apart, not representable with 2 " \
                "decimals, skipped" % " ".join(lost)
    if len(cands) > int(max_candidates):
        note += "; %d candidates > %d, keeping the %d smallest" % (
            len(cands), int(max_candidates), int(max_candidates))
        cands = cands[:int(max_candidates)]
    return cands, shells, note


def resolve_cut3_candidates(poscar, spec, cut3_max=7.0, tol=0.05, min_shells=2,
                            min_gap=0.05, max_candidates=12):
    """Parse CUT3_CANDIDATES: auto | off | "3.6 4.2 4.8".

    Returns (candidates, shells, note).  off/none/empty -> ([], [], note);
    auto -> enumerate from the structure; an explicit list is kept as given.
    """
    s = "auto" if spec is None else str(spec).strip()   # empty = off, not auto
    low = s.lower()
    if low in ("", "off", "none", "false", "no", "0"):
        return [], [], "CUT3_CANDIDATES=%s: no cutoff scan (single cutoff)" % (s or "off")
    if low == "auto":
        cands, shells, note = cut3_candidates(
            poscar, r_max=cut3_max, tol=tol, min_shells=min_shells,
            min_gap=min_gap, max_candidates=max_candidates)
        return cands, shells, "auto: " + note
    vals = []
    for tok in s.replace(",", " ").replace(";", " ").split():
        try:
            vals.append(round(float(tok), 2))
        except ValueError:
            sys.exit("[ERROR] CUT3_CANDIDATES=%r is neither auto/off nor a number list" % s)
    vals = sorted(set(vals))
    return vals, [], "explicit list: candidates %s" % vals


def nth_neighbor_cutoffs(poscar, nmax, tol=1e-4, max_atoms=600):
    """Cutoff (A) that keeps the 1st..N-th nearest neighbours, N = 1..nmax.

    Same definition as ShengBTE thirdorder.py ``-n`` (thirdorder_common.
    calc_frange), so "fc3 up to the N-th neighbours" is comparable with the
    literature: for every atom of the cell the distinct neighbour distances
    u_1 < u_2 < ... are listed (distances closer than ``tol`` are one shell;
    thirdorder uses np.allclose, i.e. ~1e-5 relative), and

        cutoff(N) = max over atoms of (u_N + u_{N+1}) / 2

    i.e. the midpoint between the N-th and (N+1)-th shell of the atom whose
    N-th shell reaches furthest.  Distances are taken in the infinite crystal
    (thirdorder takes them inside the supercell - identical as long as the
    cutoff stays below the supercell safe cutoff, which the caller enforces).

    Returns (cuts, per_atom_shells): cuts[N-1] is the cutoff for N; the lists
    in per_atom_shells hold each atom's first nmax+1 shell distances.
    """
    lat, frac = read_poscar_cell_frac(poscar)
    natom = len(frac)
    if natom > max_atoms:
        raise RuntimeError("cell has %d atoms > %d, shell enumeration too big"
                           % (natom, max_atoms))
    nmax = int(nmax)
    if nmax < 1:
        raise ValueError("nmax must be >= 1")
    r_max = 6.0
    while True:
        per_atom = [u[:nmax + 1] for u in _per_atom_shells(lat, frac, r_max, tol)]
        if all(len(u) >= nmax + 1 for u in per_atom):
            break
        if r_max > 40.0:
            raise RuntimeError("could not find %d neighbour shells within %.0f A"
                               % (nmax + 1, r_max))
        r_max *= 1.5
    cuts = [round(max(0.5 * (u[n - 1] + u[n]) for u in per_atom), 4)
            for n in range(1, nmax + 1)]
    return cuts, [[round(x, 4) for x in u] for u in per_atom]


def supercell_safe_cutoff(cell, margin=0.1):
    """Largest cutoff that avoids periodic-image double counting (A).

    0.5 * the shortest perpendicular cell height - margin.  Pure stdlib so the
    gen can warn before the job is submitted.
    """
    vol = abs(_det3(cell))
    perp = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        area = _norm(_cross(cell[j], cell[k]))
        perp.append(vol / area if area > 1e-12 else 0.0)
    return 0.5 * min(perp) - float(margin)


# ==========================================================================
# Job naming / template rendering
# ==========================================================================
def new_jobname(cwd, step_label):
    return "%s-fit-fc-thermal-%s" % (Path(cwd).name, step_label)


def _strip_doc_placeholders(text):
    """Neutralise {{X}} tokens that appear in explanatory comments.

    A multi-line substitution landing inside a comment would be emitted before
    the first #SBATCH line, and SLURM stops parsing directives at the first
    executable statement -- partition / cpus / time / output would all be
    silently ignored.
    """
    out = []
    for ln in text.split("\n"):
        s = ln.lstrip()
        if s.startswith("#") and not s.startswith("#SBATCH") and "{{" in ln:
            ln = re.sub(r"\{\{(\w+)\}\}", r"\1", ln)
        out.append(ln)
    return "\n".join(out)


def _check_sbatch_order(path):
    """Every #SBATCH directive must precede the first executable statement."""
    lines = Path(path).read_text(encoding="utf-8").split("\n")
    first = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s and not s.startswith("#"):
            first = i
            break
    if first is None:
        return
    bad = [i + 1 for i, ln in enumerate(lines[first + 1:], start=first + 1)
           if ln.strip().startswith("#SBATCH")]
    if bad:
        sys.exit("[ERROR] %s: #SBATCH on line(s) %s appear after the first "
                 "executable statement (line %d); SLURM ignores them silently. "
                 "This usually means a multi-line placeholder leaked into a "
                 "header comment -- check the template."
                 % (Path(path).name, ",".join(str(b) for b in bad), first + 1))


def write_submit(tpl_path, out_path, subs):
    """Render a submit template ({{KEY}} substitution) and validate the result."""
    text = _strip_doc_placeholders(Path(tpl_path).read_text(encoding="utf-8"))
    for k, v in subs.items():
        text = text.replace("{{%s}}" % k, str(v))
    left = sorted(set(re.findall(r"\{\{([A-Z_]+)\}\}", text)))
    if left:
        sys.exit("[ERROR] template %s still has placeholders: %s"
                 % (Path(tpl_path).name, ", ".join(left)))
    Path(out_path).write_text(text, encoding="utf-8", newline="\n")
    _check_sbatch_order(out_path)
    print("[OK] %s" % Path(out_path).name)


def resolve_submit(base_dir, kind):
    """Locate <kind>.tpl either in <base_dir>/templates/ or next to the gen
    script (gen_need may have flattened the skill's template directory)."""
    base = Path(base_dir)
    for cand in (base / "templates" / ("%s.tpl" % kind),
                 base / ("%s.tpl" % kind),
                 base / "templates" / ("%s_step1_fit.tpl" % kind)):
        if cand.is_file():
            return cand
    # 技能源码目录里模板按步骤分子目录（templates/step1_fit/submit_*.tpl）；
    # 脱离 autozt 直接运行 gen 时 gen_need 没把它们摊平到脚本旁边。
    for cand in sorted((base / "templates").glob("*/%s.tpl" % kind)):
        if cand.is_file():
            return cand
    sys.exit("[ERROR] cannot find the submit template %s.tpl (looked in %s and "
             "%s) -- is it listed in gen_need?" % (kind, base / "templates", base))
