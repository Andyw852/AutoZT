#!/usr/bin/env python3
"""Shared generators for the QE/Wannier90/Perturbo workflow.

The generator is deliberately parameter driven.  The skill defaults remain a
small smoke test, while the official Perturbo epr9 settings and convergence
cases are selected through the merged ``step.conf`` that AutoZT pushes into
each step directory.
"""
from __future__ import annotations

import json
import math
import re
import shlex
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Standard atomic weights used only when a project does not provide an
# ``ATOMIC_MASSES`` override.  Isotopically enriched calculations must set the
# masses explicitly in the project step.conf.
PERIODIC_MASSES = {
    "H": 1.008, "He": 4.002602, "Li": 6.94, "Be": 9.0121831,
    "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998403163, "Ne": 20.1797, "Na": 22.98976928,
    "Mg": 24.305, "Al": 26.9815385, "Si": 28.085, "P": 30.973761998,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.0983,
    "Ca": 40.078, "Sc": 44.955908, "Ti": 47.867, "V": 50.9415,
    "Cr": 51.9961, "Mn": 54.938044, "Fe": 55.845, "Co": 58.933194,
    "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723,
    "Ge": 72.630, "As": 74.921595, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90584,
    "Zr": 91.224, "Nb": 92.90637, "Mo": 95.95, "Tc": 98.0,
    "Ru": 101.07, "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682,
    "Cd": 112.414, "In": 114.818, "Sn": 118.710, "Sb": 121.760,
    "Te": 127.60, "I": 126.90447, "Xe": 131.293, "Cs": 132.90545196,
    "Ba": 137.327, "La": 138.90547, "Ce": 140.116, "Pr": 140.90766,
    "Nd": 144.242, "Pm": 145.0, "Sm": 150.36, "Eu": 151.964,
    "Gd": 157.25, "Tb": 158.92535, "Dy": 162.500, "Ho": 164.93033,
    "Er": 167.259, "Tm": 168.93422, "Yb": 173.045, "Lu": 174.9668,
    "Hf": 178.49, "Ta": 180.94788, "W": 183.84, "Re": 186.207,
    "Os": 190.23, "Ir": 192.217, "Pt": 195.084, "Au": 196.966569,
    "Hg": 200.592, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98040,
    "Po": 209.0, "At": 210.0, "Rn": 222.0, "Fr": 223.0,
    "Ra": 226.0, "Ac": 227.0, "Th": 232.0377, "Pa": 231.03588,
    "U": 238.02891, "Np": 237.0, "Pu": 244.0, "Am": 243.0,
    "Cm": 247.0, "Bk": 247.0, "Cf": 251.0, "Es": 252.0,
    "Fm": 257.0, "Md": 258.0, "No": 259.0, "Lr": 266.0,
    "Rf": 267.0, "Db": 268.0, "Sg": 269.0, "Bh": 270.0,
    "Hs": 269.0, "Mt": 278.0, "Ds": 281.0, "Rg": 282.0,
    "Cn": 285.0, "Nh": 286.0, "Fl": 289.0, "Mc": 290.0,
    "Lv": 293.0, "Ts": 294.0, "Og": 294.0,
}
JZZN_ENVIRONMENT = r'''module purge
module load gcc/10.2.0 openmpi/4.0.5 hdf5/1.14.6 aocl/4.1.0 fftw/3.3.8 wannier/3.1.0
export OPENMPI_BIN=/public/software/openmpi/4.0.5/bin
export INTEL_MPI_BIN=/public/software/intel/2022.3/mpi/2021.7.0/bin
export MKLROOT=/public/software/intel/2022.3/mkl/2022.2.0
# qe2pert's AOCL BLIS zdotc path is unstable on jzzn.  Preload MKL's LP64
# BLAS symbols for this smoke test while retaining AOCL for unresolved deps.
export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:/public/software/aocc/4.1.0/aocc-compiler-4.1.0/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_rt.so${LD_PRELOAD:+:$LD_PRELOAD}
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export QE_BIN=/public/home/.../software/AutoZT/qe-7.3-perturbo/bin
export WANNIER_BIN=/public/software/wannier/3.1.0/bin
export PERTURBO_BIN=/public/home/.../software/AutoZT/qe-7.3-perturbo/perturbo/bin
export PSEUDO_DIR=/public/home/.../software/AutoZT/pseudo
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
# jzzn compute nodes may expose an unusable OpenFabrics device.  These
# single-node tests use shared-memory Open MPI transport deliberately.
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader
# AOCL's Zen dot kernel is unstable with this qe2pert build on jzzn.
export BLIS_ARCH_TYPE=generic
export BLIS_NUM_THREADS=1
export PATH=$OPENMPI_BIN:$PATH
for d in "$QE_BIN" "$WANNIER_BIN" "$PERTURBO_BIN"; do
  [ -z "$d" ] || export PATH="$d:$PATH"
done
run_mpi() {
  # QE/ph.x can enter the same AOCL/MKL complex-dot path while reducing the
  # dynamical matrix (not only qe2pert).  Use the scoped shim for all MPI
  # programs so a dense q mesh does not fail after several valid q points.
  run_perturbo_omp "$OPENMPI_BIN/mpirun" -np "${SLURM_NTASKS:?SLURM_NTASKS is required}" "$@"
}
run_perturbo_omp() {
  # The tested jzzn qe2pert binary can segfault in its OpenMP zdotc path.
  # Keep the four-core allocation for QE, but run Perturbo serially for the
  # small smoke test; larger production runs need a separately validated
  # threaded Perturbo build.
  # Work around a reproducible BLAS zdotc crash in the site AOCL/MKL stack.
  # The shim implements only the GNU Fortran complex-return ABI and is scoped
  # to this job's temporary directory; QE/Perturbo binaries remain untouched.
  local shim_dir="${TMPDIR:-/tmp}/eph-qe-zdotc-${USER:-user}"
  mkdir -p "$shim_dir"
  if [ ! -s "$shim_dir/zdotc_shim_v2.so" ]; then
    cat > "$shim_dir/zdotc_shim.c" <<'EOF'
#include <complex.h>
typedef double _Complex zc;
zc zdotc_(const int *n, const zc *x, const int *incx,
          const zc *y, const int *incy) {
  zc sum = 0.0 + 0.0 * I;
  int ix = 0, iy = 0;
  for (int i = 0; i < *n; ++i) {
    sum += conj(x[ix]) * y[iy];
    ix += *incx; iy += *incy;
  }
  return sum;
}
zc zdotu_(const int *n, const zc *x, const int *incx,
          const zc *y, const int *incy) {
  zc sum = 0.0 + 0.0 * I;
  int ix = 0, iy = 0;
  for (int i = 0; i < *n; ++i) {
    sum += x[ix] * y[iy];
    ix += *incx; iy += *incy;
  }
  return sum;
}
EOF
    gcc -shared -fPIC -O2 "$shim_dir/zdotc_shim.c" -o "$shim_dir/zdotc_shim_v2.so"
  fi
  LD_PRELOAD="$shim_dir/zdotc_shim_v2.so${LD_PRELOAD:+:$LD_PRELOAD}" OMP_NUM_THREADS=1 "$@"
}
run_wannier() {
  "$INTEL_MPI_BIN/mpirun" -np 1 "$@"
}
'''
STEP = {
    "preflight": "step0_preflight", "scf": "step1_scf", "wannier": "step2_wannier",
    "phonon": "step3_phonon", "qe2pert": "step4_qe2pert", "perturbo": "step5_perturbo",
}

# Keep these defaults compatible with the original four-core smoke test.  The
# official epr9 reference is documented in README.md and is selected by a
# project-level step.conf rather than silently changing the smoke test.
CONFIG_DEFAULTS = {
    "QE_PREFIX": ("auto", "str"),
    "PSEUDO_FILES": ("auto", "str"),
    "ATOMIC_MASSES": ("auto", "str"),
    "WANNIER_PROJECTIONS": ("auto", "str"),
    "ECUTWFC": (40.0, "float"),
    "ECUTRHO": (320.0, "float"),
    "SCF_CONV_THR": (1.0e-10, "float"),
    "NSCF_CONV_THR": (1.0e-10, "float"),
    "TR2_PH": (1.0e-14, "float"),
    "QE_IBRAV": (0, "int"),
    "QE_CELLDm1": (10.264, "float"),
    "SCF_KGRID": ("4 4 4", "str"),
    "NSCF_KGRID": ("4 4 4", "str"),
    "PH_QGRID": ("2 2 2", "str"),
    "PH_KGRID": ("4 4 4", "str"),
    "Q2PERT_KGRID": ("4 4 4", "str"),
    "PH_SEARCH_SYM": ("false", "str"),
    "PH_EPSIL": ("true", "str"),
    "PH_LQDIR": ("true", "str"),
    "NUM_BANDS": (16, "int"),
    "NUM_WANN": (8, "int"),
    "DFT_BAND_MIN": (1, "int"),
    "DFT_BAND_MAX": (0, "int"),
    "DIS_NUM_ITER": (500, "int"),
    "DIS_WIN_MIN": (-100.0, "float"),
    "DIS_WIN_MAX": (17.2, "float"),
    "DIS_FROZ_MIN": (-100.0, "float"),
    "DIS_FROZ_MAX": (9.0, "float"),
    "WANNIER_NUM_ITER": (10000, "int"),
    "OCCUPATIONS": ("fixed", "str"),
    "SMEARING": ("gauss", "str"),
    "DEGAUSS": (0.02, "float"),
    "PERT_LIST_MODE": ("gamma", "str"),
    "PERT_KGRID": ("1 1 1", "str"),
    "PERT_QGRID": ("1 1 1", "str"),
    "PERT_BAND_MIN": (1, "int"),
    "PERT_BAND_MAX": (1, "int"),
    "PHFREQ_CUTOFF": (1.0, "float"),
    "PERT_DELTA_SMEAR": (10.0, "float"),
}


def _load_params() -> dict[str, object]:
    """Read the merged step.conf, with a small local fallback for unit tests."""
    values = {key: default for key, (default, _typ) in CONFIG_DEFAULTS.items()}
    conf = Path("step.conf")
    if not conf.is_file():
        return values
    try:
        import stepconf  # supplied through gen_need by AutoZT

        loaded = stepconf.load(CONFIG_DEFAULTS, cwd=".", strict=False)
        values.update(loaded.params)
        return values
    except Exception:
        # Keep direct ``python gen_step*.py`` diagnostics useful even when the
        # shared stepconf helper was not copied into a scratch directory.
        text = conf.read_text(encoding="utf-8-sig", errors="replace")
        in_params = False
        raw = {}
        for line in text.splitlines():
            line = line.split("#", 1)[0].split("!", 1)[0].strip()
            if line.startswith("["):
                in_params = line.lower() == "[params]"
            elif in_params and "=" in line:
                key, val = line.split("=", 1)
                raw[key.strip().upper()] = val.strip()
        for key, (_default, typ) in CONFIG_DEFAULTS.items():
            if key not in raw or raw[key] == "":
                continue
            try:
                if typ == "int":
                    values[key] = int(raw[key], 0)
                elif typ == "float":
                    values[key] = float(raw[key])
                else:
                    values[key] = raw[key]
            except ValueError as exc:
                raise SystemExit(f"[错误] step.conf 的 {key}={raw[key]!r} 无法解析: {exc}")
        return values


def grid(value: object, name: str) -> tuple[int, int, int]:
    """Parse ``n1 n2 n3`` (also accepting ``n1x n2x n3``) safely."""
    if isinstance(value, (tuple, list)):
        tokens = [str(item).strip() for item in value]
    else:
        tokens = str(value).lower().replace("x", " ").replace(",", " ").split()
    if len(tokens) != 3:
        raise SystemExit(f"[错误] {name} 必须是三个正整数，收到 {value!r}")
    try:
        out = tuple(int(token, 0) for token in tokens)
    except ValueError as exc:
        raise SystemExit(f"[错误] {name} 不是整数网格: {value!r}") from exc
    if any(n < 1 for n in out):
        raise SystemExit(f"[错误] {name} 必须全为正整数，收到 {value!r}")
    return out  # type: ignore[return-value]


def _fmt_grid(g: tuple[int, int, int]) -> str:
    return " ".join(str(x) for x in g)


def _qe_bool(value: object, name: str) -> str:
    """Normalize a step.conf boolean to the Fortran spelling QE accepts."""
    token = str(value).strip().lower()
    if token in {"1", "true", ".true.", "yes", "on"}:
        return ".true."
    if token in {"0", "false", ".false.", "no", "off"}:
        return ".false."
    raise SystemExit(f"[错误] {name} 必须是 true/false，收到 {value!r}")


def _canonical_symbol(value: object) -> str:
    """Return the canonical element symbol used by QE and the config maps."""
    token = str(value).strip()
    if not re.fullmatch(r"[A-Za-z]{1,2}", token):
        raise SystemExit(f"[错误] 元素符号无效：{value!r}")
    symbol = token[0].upper() + token[1:].lower()
    if symbol not in PERIODIC_MASSES:
        raise SystemExit(f"[错误] 不认识的元素符号：{value!r}")
    return symbol


def _parse_mapping(value: object, name: str, cast) -> dict[str, object]:
    """Parse ``Element:value`` pairs used by project-level step.conf."""
    text = str(value).strip()
    if not text or text.lower() == "auto":
        return {}
    # Permit the readable ``C : C.upf`` spelling as well as ``C:C.upf``.
    text = re.sub(r"([A-Za-z]{1,2})\s*([:=])\s*", r"\1\2", text)
    result: dict[str, object] = {}
    for token in re.split(r"[;,\s]+", text):
        token = token.strip().strip("'\"")
        if not token:
            continue
        if ":" in token:
            key, raw = token.split(":", 1)
        elif "=" in token:
            key, raw = token.split("=", 1)
        else:
            raise SystemExit(f"[错误] {name} 必须写成 Element:value，收到 {value!r}")
        key = _canonical_symbol(key)
        raw = raw.strip().strip("'\"")
        if not raw:
            raise SystemExit(f"[错误] {name} 中 {key} 的值为空。")
        if key in result:
            raise SystemExit(f"[错误] {name} 中重复定义了 {key}。")
        try:
            result[key] = cast(raw)
        except (TypeError, ValueError) as exc:
            raise SystemExit(f"[错误] {name} 中 {key}:{raw!r} 无法解析：{exc}") from exc
    return result


def _validate_prefix(value: object) -> str:
    prefix = str(value).strip()
    if not prefix or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", prefix):
        raise SystemExit(
            "[错误] QE_PREFIX 必须以字母或数字开头，只能包含字母、数字、下划线和连字符，"
            f"收到 {value!r}"
        )
    return prefix


def _validate_pseudo_name(value: object, element: str) -> str:
    name = str(value).strip().strip("'\"")
    if not name or Path(name).name != name or not re.fullmatch(r"[A-Za-z0-9_.+@-]+", name):
        raise SystemExit(
            f"[错误] PSEUDO_FILES 中 {element} 的文件名无效：{value!r}；"
            "只允许 PSEUDO_DIR 下的单个文件名。"
        )
    return name


def _config_manifest(params: dict[str, object], context: dict[str, object]) -> dict[str, object]:
    """Return resolved parameters and structure data for downstream steps."""
    out: dict[str, object] = {}
    for key, value in params.items():
        if key.endswith("GRID"):
            out[key] = list(grid(value, key))
        else:
            out[key] = value
    out["workflow_defaults"] = "smoke-compatible"
    out.update({
        "prefix": context["prefix"],
        "species": context["species"],
        "counts": context["counts"],
        "nat": context["nat"],
        "ntyp": context["ntyp"],
        "atom_symbols": context["atom_symbols"],
        "pseudo_files": context["pseudo_files"],
        "atomic_masses": context["atomic_masses"],
        "wannier_projections": context["wannier_projections"],
        "dft_band_min": context["dft_band_min"],
        "dft_band_max": context["dft_band_max"],
    })
    return out


def _write_manifest(dest: Path, params: dict[str, object], context: dict[str, object]) -> None:
    write(dest / "eph_parameters.json", json.dumps(
        _config_manifest(params, context), indent=2
    ) + "\n")


# Fixed official epr9 random lists.  Keeping the coordinates in the generator
# makes the reference case reproducible without requiring a remote data path.
OFFICIAL_K_POINTS = """0.41590 0.23200 0.44382 1.0
0.01231 0.94492 0.92699 1.0
0.74245 0.82150 0.78308 1.0
0.41229 0.73562 0.67200 1.0
0.40727 0.12882 0.26821 1.0
0.58833 0.11904 0.96421 1.0
0.06981 0.49834 0.94893 1.0
0.90108 0.51877 0.67738 1.0
0.41915 0.83242 0.75735 1.0
0.82294 0.38240 0.24800 1.0
0.30474 0.70112 0.73272 1.0
0.21087 0.59605 0.52871 1.0
0.53058 0.49604 0.11970 1.0
0.39451 0.69645 0.63207 1.0
0.19735 0.66564 0.58967 1.0
0.50792 0.25265 0.39838 1.0
0.56073 0.24573 0.57818 1.0
0.15488 0.35203 0.97552 1.0
0.40710 0.03491 0.19385 1.0
0.01082 0.28881 0.77249 1.0
0.03942 0.54512 0.68440 1.0
0.14551 0.56228 0.08296 1.0
0.37236 0.45945 0.81430 1.0
0.01384 0.82308 0.38737 1.0
0.70167 0.86882 0.51744 1.0"""

OFFICIAL_Q_POINTS = """0.13843 0.58899 0.50753 1.0
0.78749 0.90377 0.94826 1.0
0.89376 0.13029 0.30233 1.0
0.10014 0.18902 0.07285 1.0
0.97397 0.00532 0.13379 1.0
0.76655 0.58433 0.89558 1.0
0.64521 0.44463 0.93646 1.0
0.12299 0.85290 0.39485 1.0
0.36861 0.89833 0.68771 1.0
0.50179 0.08177 0.64884 1.0
0.18061 0.01049 0.14836 1.0
0.84423 0.33316 0.22903 1.0
0.37995 0.01986 0.24198 1.0
0.31187 0.79472 0.77745 1.0
0.54355 0.05365 0.08732 1.0
0.88533 0.42545 0.41231 1.0
0.91184 0.53854 0.10149 1.0
0.44342 0.56113 0.64574 1.0
0.59596 0.73704 0.29634 1.0
0.49715 0.37671 0.37870 1.0
0.23038 0.84571 0.12669 1.0
0.48633 0.20884 0.54150 1.0
0.85383 0.51393 0.87155 1.0
0.08836 0.85234 0.43671 1.0
0.00544 0.82280 0.35852 1.0"""

def _poscar() -> tuple[list[list[float]], list[str], list[int], list[list[float]]]:
    p = Path("POSCAR")
    if not p.is_file():
        raise SystemExit("[错误] 材料根目录缺 POSCAR。")
    ls = [x.strip() for x in p.read_text(encoding="utf-8", errors="replace").splitlines() if x.strip()]
    if len(ls) < 8:
        raise SystemExit("[错误] POSCAR 不完整。")
    try:
        scale = float(ls[1])
        raw_cell = [[float(x) for x in ls[i].split()[:3]] for i in range(2, 5)]
        species = ls[5].split()
        counts = [int(x) for x in ls[6].split()]
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"[错误] 仅支持 VASP 5 格式 POSCAR：{exc}")
    if len(species) != len(counts) or not species or any(n < 1 for n in counts):
        raise SystemExit("[错误] POSCAR 元素数和原子数不匹配，且每种元素至少要有一个原子。")
    if scale == 0:
        raise SystemExit("[错误] POSCAR 缩放因子不能为 0。")
    raw_det = (
        raw_cell[0][0] * (raw_cell[1][1] * raw_cell[2][2] - raw_cell[1][2] * raw_cell[2][1])
        - raw_cell[0][1] * (raw_cell[1][0] * raw_cell[2][2] - raw_cell[1][2] * raw_cell[2][0])
        + raw_cell[0][2] * (raw_cell[1][0] * raw_cell[2][1] - raw_cell[1][1] * raw_cell[2][0])
    )
    if abs(raw_det) < 1.0e-12:
        raise SystemExit("[错误] POSCAR 晶胞奇异。")
    lattice_scale = scale if scale > 0 else (abs(scale) / abs(raw_det)) ** (1.0 / 3.0)
    cell = [[component * lattice_scale for component in vector] for vector in raw_cell]
    pos = 7
    if ls[pos].lower().startswith("s"):
        pos += 1
    if pos >= len(ls) or not ls[pos].lower().startswith(("d", "c")):
        raise SystemExit("[错误] POSCAR 缺 Direct/Cartesian 坐标标记。")
    direct = ls[pos].lower().startswith("d"); pos += 1
    n = sum(counts)
    if len(ls) < pos + n:
        raise SystemExit(f"[错误] POSCAR 坐标行不足：需要 {n} 行，实际只有 {len(ls) - pos} 行。")
    try:
        coords = [[float(x) for x in ls[i].split()[:3]] for i in range(pos, pos + n)]
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"[错误] POSCAR 原子坐标解析失败：{exc}")
    if not direct:
        # POSCAR 晶格矢量按行保存，Cartesian = fractional @ cell。
        coords = [[component * lattice_scale for component in vector] for vector in coords]
        a, b, c = cell
        det = sum(a[i] * (b[(i+1) % 3] * c[(i+2) % 3] - b[(i+2) % 3] * c[(i+1) % 3]) for i in range(3))
        if abs(det) < 1.0e-12:
            raise SystemExit("[错误] POSCAR 晶胞奇异。")
        inv = [
            [(b[1]*c[2]-b[2]*c[1])/det, (a[2]*c[1]-a[1]*c[2])/det, (a[1]*b[2]-a[2]*b[1])/det],
            [(b[2]*c[0]-b[0]*c[2])/det, (a[0]*c[2]-a[2]*c[0])/det, (a[2]*b[0]-a[0]*b[2])/det],
            [(b[0]*c[1]-b[1]*c[0])/det, (a[1]*c[0]-a[0]*c[1])/det, (a[0]*b[1]-a[1]*b[0])/det],
        ]
        coords = [[sum(v[j] * inv[j][i] for j in range(3)) for i in range(3)] for v in coords]
    return cell, species, counts, coords


def _material_context(params: dict[str, object] | None = None) -> dict[str, object]:
    """Resolve POSCAR and project-level mappings once for every generator."""
    params = _load_params() if params is None else params
    cell, raw_species, counts, coords = _poscar()
    species = [_canonical_symbol(symbol) for symbol in raw_species]
    if len(set(species)) != len(species):
        raise SystemExit("[错误] POSCAR 元素列表含重复元素；请把同一元素合并为一个种类。")
    atom_symbols = [symbol for symbol, count in zip(species, counts) for _ in range(count)]

    raw_prefix = str(params.get("QE_PREFIX", "auto")).strip()
    if not raw_prefix or raw_prefix.lower() == "auto":
        if species == ["Si"] and counts == [2]:
            prefix = "si"
        else:
            prefix = "".join(f"{symbol.lower()}{count}" for symbol, count in zip(species, counts))
    else:
        prefix = raw_prefix
    prefix = _validate_prefix(prefix)

    pseudo_map = _parse_mapping(params.get("PSEUDO_FILES", "auto"), "PSEUDO_FILES", str)
    if not pseudo_map:
        pseudo_map = {symbol: f"{symbol}.upf" for symbol in species}
    missing = [symbol for symbol in species if symbol not in pseudo_map]
    if missing:
        raise SystemExit(
            "[错误] PSEUDO_FILES 缺少元素映射：" + ", ".join(missing)
            + "；例如 PSEUDO_FILES = C:C.upf Si:Si.upf"
        )
    pseudo_files = {symbol: _validate_pseudo_name(pseudo_map[symbol], symbol) for symbol in species}

    mass_map = _parse_mapping(params.get("ATOMIC_MASSES", "auto"), "ATOMIC_MASSES", float)
    if not mass_map:
        mass_map = {symbol: PERIODIC_MASSES[symbol] for symbol in species}
    missing = [symbol for symbol in species if symbol not in mass_map]
    if missing:
        raise SystemExit("[错误] ATOMIC_MASSES 缺少元素映射：" + ", ".join(missing))
    atomic_masses = {symbol: float(mass_map[symbol]) for symbol in species}
    if any(not math.isfinite(mass) or mass <= 0 for mass in atomic_masses.values()):
        raise SystemExit("[错误] ATOMIC_MASSES 必须全部为有限正数。")

    projection_map = _parse_mapping(
        params.get("WANNIER_PROJECTIONS", "auto"), "WANNIER_PROJECTIONS", str
    )
    if not projection_map:
        if species == ["Si"] and counts == [2]:
            projection_map = {"Si": "sp3"}
        else:
            raise SystemExit(
                "[错误] 无法从该结构安全推断 WANNIER_PROJECTIONS；"
                "请按元素显式设置，例如 WANNIER_PROJECTIONS = C:sp2 Si:sp3。"
            )
    missing = [symbol for symbol in species if symbol not in projection_map]
    if missing:
        raise SystemExit("[错误] WANNIER_PROJECTIONS 缺少元素映射：" + ", ".join(missing))
    wannier_projections = {symbol: str(projection_map[symbol]) for symbol in species}

    try:
        nbands = int(params.get("NUM_BANDS", 16))
        num_wann = int(params.get("NUM_WANN", 8))
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"[错误] NUM_BANDS/NUM_WANN 必须为整数：{exc}") from exc
    if nbands < 1 or num_wann < 1:
        raise SystemExit("[错误] NUM_BANDS 和 NUM_WANN 必须为正整数。")
    if num_wann > nbands:
        raise SystemExit(f"[错误] NUM_WANN={num_wann} 不能大于 NUM_BANDS={nbands}。")
    try:
        dft_band_min = int(params.get("DFT_BAND_MIN", 1))
        dft_band_max = int(params.get("DFT_BAND_MAX", 0)) or nbands
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"[错误] DFT_BAND_MIN/DFT_BAND_MAX 必须为整数：{exc}") from exc
    if dft_band_min < 1 or dft_band_max < dft_band_min or dft_band_max > nbands:
        raise SystemExit(
            f"[错误] DFT 能带窗无效：{dft_band_min}..{dft_band_max}，NUM_BANDS={nbands}。"
        )
    return {
        "cell": cell, "species": species, "counts": counts, "coords": coords,
        "nat": len(atom_symbols), "ntyp": len(species), "atom_symbols": atom_symbols,
        "prefix": prefix, "pseudo_files": pseudo_files, "atomic_masses": atomic_masses,
        "wannier_projections": wannier_projections, "num_bands": nbands,
        "num_wann": num_wann, "dft_band_min": dft_band_min, "dft_band_max": dft_band_max,
    }


def structure_qe(include_cell: bool = True, params: dict[str, object] | None = None,
                 context: dict[str, object] | None = None) -> str:
    context = _material_context(params) if context is None else context
    text = ""
    if include_cell:
        text += "CELL_PARAMETERS angstrom\n" + "\n".join(
            "  %.10f %.10f %.10f" % tuple(v) for v in context["cell"]
        ) + "\n"
    text += "\nATOMIC_SPECIES\n"
    for symbol in context["species"]:
        text += f"{symbol} {context['atomic_masses'][symbol]:.8g} {context['pseudo_files'][symbol]}\n"
    text += "ATOMIC_POSITIONS crystal\n"
    return text + "\n".join(
        f"{symbol} %.10f %.10f %.10f" % tuple(coord)
        for symbol, coord in zip(context["atom_symbols"], context["coords"])
    ) + "\n"

def write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8", newline="\n")

def render_submit(dest: Path, job: str, command: str) -> None:
    tpl = HERE / "submit_qe_test.tpl"
    if not tpl.is_file(): raise SystemExit("[错误] 缺 submit_qe_test.tpl。")
    write(dest / "submit.sh", tpl.read_text(encoding="utf-8").replace("{{JOBNAME}}", job).replace("{{COMMAND}}", command))
    (dest / "submit.sh").chmod(0o755)

def make_step(kind: str, command: str, files: dict[str, str],
              params: dict[str, object] | None = None) -> Path:
    params = _load_params() if params is None else params
    context = _material_context(params)
    dest = Path(STEP[kind]); dest.mkdir(parents=True, exist_ok=True)
    for name, text in files.items(): write(dest / name, text)
    _write_manifest(dest, params, context)
    render_submit(dest, f"{context['prefix']}-eph-{kind}", JZZN_ENVIRONMENT + command)
    print(f"[DONE] {dest} 输入已生成；待 AutoZT 提交。")
    return dest

def uniform_kpoints(grid_value: object) -> str:
    """Return a full uniform grid with QE/Perturbo weights."""
    n1, n2, n3 = grid(grid_value, "KGRID")
    total = n1 * n2 * n3
    points = [
        (i / n1, j / n2, k / n3, 1.0 / total)
        for i in range(n1) for j in range(n2) for k in range(n3)
    ]
    return "K_POINTS crystal\n%d\n" % len(points) + "\n".join(
        "  %.10f %.10f %.10f %.12e" % point for point in points
    ) + "\n"

def qe_input(calculation: str, kgrid: object, params: dict[str, object], nbnd: int | None = None) -> str:
    context = _material_context(params)
    # Keep the NSCF mesh as an automatic grid so qe2pert can read nk1/nk2/nk3
    # from QE's XML metadata; Wannier90 still receives the explicit list in
    # its .win file below.
    kgrid_text = _fmt_grid(grid(kgrid, f"{calculation.upper()}_KGRID"))
    kpoints = (uniform_kpoints(kgrid) if calculation == "nscf"
               else "K_POINTS automatic\n  " + kgrid_text + " 0 0 0\n")
    # Wannier90 requires an explicit full-BZ NSCF mesh.  Keep the SCF save
    # symmetry-consistent for ph.x: qe2pert reconstructs the complete q mesh
    # from its irreducible dynamical-matrix files and their q-stars.
    nscf_flags = "  nosym = .true.\n  noinv = .true.\n" if calculation == "nscf" else ""
    occupations = str(params.get("OCCUPATIONS", "fixed")).strip().lower()
    smearing = ""
    if occupations not in {"fixed", "none"}:
        smearing = (f"  occupations = '{occupations}'\n"
                    f"  smearing = '{str(params.get('SMEARING', 'gauss')).strip()}'\n"
                    f"  degauss = {float(params.get('DEGAUSS', 0.02)):.8g}\n")
    nbnd_value = context["num_bands"] if nbnd is None else int(nbnd)
    if nbnd_value < 1:
        raise SystemExit("[错误] QE 的 nbnd 必须为正整数。")
    conv_key = "NSCF_CONV_THR" if calculation == "nscf" else "SCF_CONV_THR"
    conv_thr = float(params.get(conv_key, 1.0e-10))
    qe_ibrav = int(params.get("QE_IBRAV", 0))
    if qe_ibrav not in (0, 2):
        raise SystemExit(f"[错误] QE_IBRAV 目前仅支持 0 或 2，收到 {qe_ibrav}")
    if qe_ibrav == 2:
        lattice = (
            "  ibrav = 2\n"
            f"  celldm(1) = {float(params.get('QE_CELLDm1', 10.264)):.8g}\n"
        )
        cell_block = structure_qe(include_cell=False, context=context)
    else:
        lattice = "  ibrav = 0\n"
        cell_block = structure_qe(include_cell=True, context=context)
    return f"""&CONTROL
  calculation = '{calculation}'
  prefix = '{context['prefix']}'
  outdir = './tmp'
  pseudo_dir = './pseudo'
  wf_collect = .true.
/
&SYSTEM
{lattice}  nat = {context['nat']}
  ntyp = {context['ntyp']}
  ecutwfc = {float(params.get('ECUTWFC', 40.0)):.8g}
  ecutrho = {float(params.get('ECUTRHO', 320.0)):.8g}
  nbnd = {nbnd_value}
{smearing}{nscf_flags}/
&ELECTRONS
  conv_thr = {conv_thr:.8e}
  mixing_beta = 0.7
/
{cell_block}{kpoints}"""

def wannier_cell(params: dict[str, object], context: dict[str, object] | None = None) -> tuple[str, list[tuple[float, float, float]]]:
    """Return the cell in the same convention that QE writes to its XML.

    QE ``ibrav=2`` derives its primitive vectors from ``celldm(1)`` rather
    than from the POSCAR vectors.  Wannier90 passes ``unit_cell_cart`` back to
    ``pw2wannier90.x``, which checks the lattice exactly; using the POSCAR
    representation here therefore fails after switching to the official epr9
    input even when both cells describe the same physical crystal.
    """
    if int(params.get("QE_IBRAV", 0)) == 2:
        alat = float(params.get("QE_CELLDm1", 10.264))
        return "bohr", [
            (-0.5 * alat, 0.0, 0.5 * alat),
            (0.0, 0.5 * alat, 0.5 * alat),
            (-0.5 * alat, 0.5 * alat, 0.0),
        ]
    context = _material_context(params) if context is None else context
    return "angstrom", [tuple(v) for v in context["cell"]]


def _pseudo_links(context: dict[str, object]) -> str:
    lines = ["mkdir -p pseudo"]
    for symbol in context["species"]:
        filename = shlex.quote(str(context["pseudo_files"][symbol]))
        lines.append(f'ln -sfn "$PSEUDO_DIR/{filename}" "pseudo/{filename}"')
    return "\n".join(lines)


def _replace_prefix(text: str, prefix: str) -> str:
    return text.replace("__PREFIX__", prefix)

def main(kind: str) -> int:
    params = _load_params()
    context = _material_context(params)
    prefix = str(context["prefix"])
    if kind == "preflight":
        probe = HERE / "run_preflight.py"
        if not probe.is_file():
            raise SystemExit("[错误] 缺 run_preflight.py。")
        make_step(kind, "python3 run_preflight.py",
                  {"run_preflight.py": probe.read_text(encoding="utf-8")}, params)
    elif kind == "scf":
        command = _pseudo_links(context) + "\nrun_mpi pw.x -in scf.in > scf.out\ngrep -q 'JOB DONE' scf.out\nprintf 'OK\\n' > scf.ok"
        make_step(kind, command, {"scf.in": qe_input("scf", params["SCF_KGRID"], params)}, params)
    elif kind == "wannier":
        nscf_grid = grid(params["NSCF_KGRID"], "NSCF_KGRID")
        kpoints = uniform_kpoints(nscf_grid).splitlines()[2:]
        cell_unit, cell = wannier_cell(params, context)
        projections = "\n".join(f"{symbol}:{context['wannier_projections'][symbol]}" for symbol in context["species"])
        atoms = "\n".join(f"{symbol} %.10f %.10f %.10f" % tuple(coord) for symbol, coord in zip(context["atom_symbols"], context["coords"]))
        win = f"""num_bands = {context['num_bands']}
num_wann = {context['num_wann']}
write_u_matrices = true
write_hr = true
write_xyz = true
dis_num_iter = {int(params['DIS_NUM_ITER'])}
dis_win_min = {float(params['DIS_WIN_MIN']):.8g}
dis_win_max = {float(params['DIS_WIN_MAX']):.8g}
dis_froz_min = {float(params['DIS_FROZ_MIN']):.8g}
dis_froz_max = {float(params['DIS_FROZ_MAX']):.8g}
num_iter = {int(params['WANNIER_NUM_ITER'])}
begin projections
{projections}
end projections
mp_grid = {nscf_grid[0]} {nscf_grid[1]} {nscf_grid[2]}
begin unit_cell_cart
{cell_unit}
""" + "\n".join("%.10f %.10f %.10f" % tuple(x) for x in cell) + f"\nend unit_cell_cart\nbegin atoms_frac\n{atoms}\nend atoms_frac\nbegin kpoints\n" + "\n".join(" ".join(line.split()[:3]) for line in kpoints) + "\nend kpoints\n"
        pw2 = f"&inputpp\n  outdir='./tmp'\n  prefix='{prefix}'\n  seedname='{prefix}'\n  write_mmn=.true.\n  write_amn=.true.\n  write_unk=.false.\n/\n"
        command = _replace_prefix("ln -sfn ../step1_scf/tmp tmp\n" + _pseudo_links(context) + "\nrun_mpi pw.x -in nscf.in > nscf.out\n# The site Wannier90 binary is linked against Intel MPI; use its launcher.\nrun_wannier wannier90.x -pp __PREFIX__ > wannier_pp.out\nrun_mpi pw2wannier90.x -in pw2wannier.in > pw2wannier.out\nrun_wannier wannier90.x __PREFIX__ > wannier.out\nprintf 'OK\\n' > wannier.ok", prefix)
        make_step(kind, command, {"nscf.in": qe_input("nscf", nscf_grid, params), f"{prefix}.win": win, "pw2wannier.in": pw2}, params)
    elif kind == "phonon":
        qgrid = grid(params["PH_QGRID"], "PH_QGRID")
        phgrid = grid(params["PH_KGRID"], "PH_KGRID")
        search_sym = _qe_bool(params.get("PH_SEARCH_SYM", "false"), "PH_SEARCH_SYM")
        epsil = _qe_bool(params.get("PH_EPSIL", "true"), "PH_EPSIL")
        lqdir = _qe_bool(params.get("PH_LQDIR", "true"), "PH_LQDIR")
        ph = f"""&inputph
  verbosity='debug'
  prefix='{prefix}'
  outdir='../step1_scf/tmp'
  fildyn='{prefix}.dyn.xml'
  fildvscf='dvscf'
  epsil={epsil}
  ldisp=.true.
  lqdir={lqdir}
  search_sym={search_sym}
  tr2_ph={float(params.get("TR2_PH", 1.0e-14)):.8e}
  nq1={qgrid[0]}, nq2={qgrid[1]}, nq3={qgrid[2]}
  nk1={phgrid[0]}, nk2={phgrid[1]}, nk3={phgrid[2]}
 /
"""
        # qe2pert reads a flattened phonon directory, not QE's _ph0 tree.
        # Preserve the QE-generated numbering so each irreducible q point has
        # its matching dynmat, displacement pattern, and dvscf perturbation.
        expected_qgrid = "%d %d %d" % qgrid
        expected_nq = qgrid[0] * qgrid[1] * qgrid[2]
        collect = r'''# Reuse a completed DFPT payload when a previous run only failed
# during the local collection/validation tail.  This avoids repeating the
# expensive ph.x calculation after a harmless checker fix.
rm -f phonon.ok
reuse_saved=0
expected_qgrid="__EXPECTED_QGRID__"
saved_qgrid=""
if [ -s save/__PREFIX__.dyn0 ]; then
  saved_qgrid=$(head -n 1 save/__PREFIX__.dyn0 | xargs)
fi
saved_nq=""
if [ -s save/__PREFIX__.dyn0 ]; then
  saved_nq=$(sed -n '2p' save/__PREFIX__.dyn0 | awk '{print $1}')
fi
expected_nq="__EXPECTED_NQ__"
if [ "$saved_qgrid" = "$expected_qgrid" ] && [ "$saved_nq" = "$expected_nq" ] \
   && [ -d save/__PREFIX__.phsave ] && [ -s save/__PREFIX__.dvscf_q1 ]; then
  reuse_saved=1
fi
mkdir -p save
if [ "$reuse_saved" -eq 0 ]; then
# Archive only QE's phonon scratch tree for this step.  A stale _ph0
# can contain dvscf files from a different FFT grid; keep it recoverable while
# ensuring ph.x starts with a clean tree.
if [ -d ../step1_scf/tmp/_ph0 ]; then
  mv ../step1_scf/tmp/_ph0 ../step1_scf/tmp/_ph0.stale.${SLURM_JOB_ID:-manual}
fi
run_mpi ph.x -in ph.in > ph.out
phroot='../step1_scf/tmp/_ph0'
cp -f __PREFIX__.dyn* save/
test -d "$phroot/__PREFIX__.phsave"
cp -a "$phroot/__PREFIX__.phsave" save/

shopt -s nullglob
for qdir in "$phroot"/__PREFIX__.q_*; do
  iq=${qdir##*.q_}
  case "$iq" in (*[!0-9]*|'') exit 1;; esac
  dvscf=("$qdir"/__PREFIX__.dvscf*)
  if [ ${#dvscf[@]} -gt 0 ]; then
    cp -f "${dvscf[0]}" "save/__PREFIX__.dvscf_q${iq}"
  elif [ "$iq" = 1 ]; then
    # QE stores Gamma's response in _ph0/ rather than a q_1 directory on this build.
    gamma_dvscf=("$phroot"/__PREFIX__.dvscf*)
    test ${#gamma_dvscf[@]} -gt 0
    cp -f "${gamma_dvscf[0]}" "save/__PREFIX__.dvscf_q${iq}"
  else
    exit 1
  fi
done
# QE may store the Gamma-point response directly in _ph0 rather than a q_1 directory.
if [ ! -s save/__PREFIX__.dvscf_q1 ]; then
  gamma_dvscf=("$phroot"/__PREFIX__.dvscf*)
  test ${#gamma_dvscf[@]} -gt 0
  cp -f "${gamma_dvscf[0]}" save/__PREFIX__.dvscf_q1
fi
fi

test -s save/__PREFIX__.dyn0
nq_files=$(sed -n '2p' save/__PREFIX__.dyn0 | awk '{print $1}')
case "$nq_files" in (*[!0-9]*|'') exit 1;; esac
# QE stores only irreducible q points here; qe2pert reads each XML q-star and
# verifies that their union fills the configured uniform mesh.  Requiring the
# number of files to equal nq1*nq2*nq3 would incorrectly force no-symmetry
# phonons and duplicates q/-q entries for dense grids.
test "$nq_files" -gt 0
test "$nq_files" -le "__EXPECTED_NQ__"
for iq in $(seq 1 "$nq_files"); do
  test -s "save/__PREFIX__.dyn${iq}.xml"
  test -s "save/__PREFIX__.phsave/patterns.${iq}.xml"
  test -s "save/__PREFIX__.dvscf_q${iq}"
done
        printf 'OK\n' > phonon.ok'''
        collect = collect.replace("__EXPECTED_QGRID__", expected_qgrid)
        collect = collect.replace("__EXPECTED_NQ__", str(expected_nq))
        collect = _replace_prefix(collect, prefix)
        make_step(kind, collect, {"ph.in": ph}, params)
    elif kind == "qe2pert":
        nscf_grid = grid(params["NSCF_KGRID"], "NSCF_KGRID")
        ph_grid = grid(params["PH_QGRID"], "PH_QGRID")
        q2p = grid(params["Q2PERT_KGRID"], "Q2PERT_KGRID")
        if nscf_grid != q2p:
            raise SystemExit(
                "[错误] qe2pert 的 nk1/nk2/nk3 必须与 NSCF_KGRID 一致，"
                f"收到 NSCF_KGRID={_fmt_grid(nscf_grid)}、Q2PERT_KGRID={_fmt_grid(q2p)}"
            )
        if any(k % q != 0 for k, q in zip(nscf_grid, ph_grid)):
            raise SystemExit(
                "[错误] NSCF_KGRID 必须是 PH_QGRID 的整数倍，"
                "以保证 qe2pert 的 k/q 网格可相容；"
                f"收到 NSCF_KGRID={_fmt_grid(nscf_grid)}、PH_QGRID={_fmt_grid(ph_grid)}"
            )
        inp = f"""&qe2pert
  prefix='{prefix}'
  outdir='../step1_scf/tmp'
  phdir='../step3_phonon/save'
  nk1={nscf_grid[0]}, nk2={nscf_grid[1]}, nk3={nscf_grid[2]}
  dft_band_min={context['dft_band_min']}
  dft_band_max={context['dft_band_max']}
  num_wann={context['num_wann']}
  lwannier=.true.
/
 """
        # qe2pert writes prefix_epr.h5 to its QE outdir, which is this step.
        expected_qgrid = " ".join(str(x) for x in ph_grid)
        qe2pert_cmd = """ln -sfn ../step2_wannier/__PREFIX__.win __PREFIX__.win
ln -sfn ../step2_wannier/__PREFIX___u.mat __PREFIX___u.mat
ln -sfn ../step2_wannier/__PREFIX___u_dis.mat __PREFIX___u_dis.mat
ln -sfn ../step2_wannier/__PREFIX___centres.xyz __PREFIX___centres.xyz
test -s ../step3_phonon/save/__PREFIX__.dyn0
saved_qgrid=$(head -n 1 ../step3_phonon/save/__PREFIX__.dyn0 | xargs)
test \"$saved_qgrid\" = \"__EXPECTED_QGRID__\" || { echo \"wrong phonon q grid: got '$saved_qgrid', expected '__EXPECTED_QGRID__'\" >&2; exit 1; }
run_perturbo_omp qe2pert.x -npools 1 -in qe2pert.in > qe2pert.out
test -s __PREFIX___epr.h5
cp -f __PREFIX___epr.h5 eph.h5
printf 'OK\\n' > qe2pert.ok"""
        qe2pert_cmd = qe2pert_cmd.replace("__EXPECTED_QGRID__", expected_qgrid)
        qe2pert_cmd = _replace_prefix(qe2pert_cmd, prefix)
        make_step(kind, qe2pert_cmd, {"qe2pert.in": inp}, params)
    elif kind == "perturbo":
        mode = str(params.get("PERT_LIST_MODE", "gamma")).strip().lower()
        if mode == "official":
            kfile = "25\n" + OFFICIAL_K_POINTS + "\n"
            qfile = "25\n" + OFFICIAL_Q_POINTS + "\n"
        elif mode == "grid":
            kfile = uniform_kpoints(params["PERT_KGRID"]).replace("K_POINTS crystal\n", "")
            qfile = uniform_kpoints(params["PERT_QGRID"]).replace("K_POINTS crystal\n", "")
        elif mode in {"gamma", "single"}:
            kfile = "1\n0.0 0.0 0.0 1.0\n"
            qfile = "1\n0.0 0.0 0.0 1.0\n"
        else:
            raise SystemExit("[错误] PERT_LIST_MODE 只能是 gamma、grid 或 official。")
        pert_min = int(params["PERT_BAND_MIN"])
        pert_max = int(params["PERT_BAND_MAX"])
        if pert_min < 1 or pert_max < pert_min or pert_max > context["num_bands"]:
            raise SystemExit(
                f"[错误] Perturbo 能带窗无效：{pert_min}..{pert_max}，"
                f"NUM_BANDS={context['num_bands']}。"
            )
        inp = f"""&perturbo
  prefix='{prefix}'
  calc_mode='ephmat'
  fklist='eph.kpt'
  fqlist='eph.qpt'
  band_min={pert_min}
  band_max={pert_max}
  phfreq_cutoff={float(params['PHFREQ_CUTOFF']):.8g}
  delta_smear={float(params['PERT_DELTA_SMEAR']):.8g}
/
 """
        command = _replace_prefix(
            "ln -sfn ../step4_qe2pert/eph.h5 __PREFIX___epr.h5\n"
            "run_perturbo_omp perturbo.x -npools 1 -in pert.in > pert.out\n"
            "test -s __PREFIX__.ephmat\n"
            "cp -f __PREFIX__.ephmat ephmat\n"
            "test -s __PREFIX___ephmat.yml\n"
            "cp -f __PREFIX___ephmat.yml ephmat.yml\n"
            "printf 'OK\\n' > perturbo.ok",
            prefix,
        )
        make_step(kind, command, {"pert.in": inp, "eph.kpt": kfile, "eph.qpt": qfile}, params)
    else: raise SystemExit(f"unknown generator kind: {kind}")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
