#!/usr/bin/env python3
"""Shared generators for the reproducible QE/Perturbo Si workflow.

The generator is deliberately parameter driven.  The skill defaults remain a
small smoke test, while the official Perturbo epr9 settings and convergence
cases are selected through the merged ``step.conf`` that AutoZT pushes into
each step directory.
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREFIX = "si"
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
export QE_BIN=/public/home/wangchao/software/AutoZT/qe-7.3-perturbo/bin
export WANNIER_BIN=/public/software/wannier/3.1.0/bin
export PERTURBO_BIN=/public/home/wangchao/software/AutoZT/qe-7.3-perturbo/perturbo/bin
export PSEUDO_DIR=/public/home/wangchao/software/AutoZT/pseudo
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
  # small Si smoke test; larger production runs need a separately validated
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
    "NUM_BANDS": (16, "int"),
    "NUM_WANN": (8, "int"),
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


def _config_manifest(params: dict[str, object]) -> dict[str, object]:
    """Return JSON-safe values, including normalized grids, for provenance."""
    out: dict[str, object] = {}
    for key, value in params.items():
        if key.endswith("GRID"):
            out[key] = list(grid(value, key))
        else:
            out[key] = value
    out["workflow_defaults"] = "smoke-compatible"
    return out


def _write_manifest(dest: Path, params: dict[str, object]) -> None:
    write(dest / "eph_parameters.json", json.dumps(_config_manifest(params), indent=2) + "\n")


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
    if len(ls) < 9:
        raise SystemExit("[错误] POSCAR 不完整。")
    try:
        scale = float(ls[1]); cell = [[float(x) * scale for x in ls[i].split()[:3]] for i in range(2, 5)]
        species = ls[5].split(); counts = [int(x) for x in ls[6].split()]
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"[错误] 仅支持 VASP 5 格式 POSCAR：{exc}")
    pos = 7
    if ls[pos].lower().startswith("s"):
        pos += 1
    if not ls[pos].lower().startswith(("d", "c")):
        raise SystemExit("[错误] POSCAR 缺 Direct/Cartesian 坐标标记。")
    direct = ls[pos].lower().startswith("d"); pos += 1
    n = sum(counts)
    try:
        coords = [[float(x) for x in ls[i].split()[:3]] for i in range(pos, pos + n)]
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"[错误] POSCAR 原子坐标解析失败：{exc}")
    if not direct:
        # QE 仍接收 crystal；将 Cartesian 转 fractional，避免写一个只对 Si 有效的晶胞。
        a, b, c = cell
        det = sum(a[i] * (b[(i+1)%3] * c[(i+2)%3] - b[(i+2)%3] * c[(i+1)%3]) for i in range(3))
        if abs(det) < 1e-12: raise SystemExit("[错误] POSCAR 晶胞奇异。")
        # POSCAR 晶格矢量按行保存，Cartesian = fractional @ cell。
        inv = [
            [(b[1]*c[2]-b[2]*c[1])/det, (a[2]*c[1]-a[1]*c[2])/det, (a[1]*b[2]-a[2]*b[1])/det],
            [(b[2]*c[0]-b[0]*c[2])/det, (a[0]*c[2]-a[2]*c[0])/det, (a[2]*b[0]-a[0]*b[2])/det],
            [(b[0]*c[1]-b[1]*c[0])/det, (a[1]*c[0]-a[0]*c[1])/det, (a[0]*b[1]-a[1]*b[0])/det],
        ]
        coords = [[sum(v[j] * inv[j][i] for j in range(3)) for i in range(3)] for v in coords]
    return cell, species, counts, coords

def structure_qe(include_cell: bool = True) -> str:
    cell, species, counts, coords = _poscar()
    if species != ["Si"] or counts != [2]:
        raise SystemExit("[错误] v0.1 的端到端验证仅支持两原子金刚石 Si POSCAR；通用元素映射将在后续版本实现。")
    text = ""
    if include_cell:
        text += "CELL_PARAMETERS angstrom\n" + "\n".join("  %.10f %.10f %.10f" % tuple(v) for v in cell)
        text += "\n"
    text += "\nATOMIC_SPECIES\nSi 28.085 Si.upf\nATOMIC_POSITIONS crystal\n"
    return text + "\n".join("Si %.10f %.10f %.10f" % tuple(v) for v in coords) + "\n"

def write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8", newline="\n")

def render_submit(dest: Path, job: str, command: str) -> None:
    tpl = HERE / "submit_qe_test.tpl"
    if not tpl.is_file(): raise SystemExit("[错误] 缺 submit_qe_test.tpl。")
    write(dest / "submit.sh", tpl.read_text(encoding="utf-8").replace("{{JOBNAME}}", job).replace("{{COMMAND}}", command))
    (dest / "submit.sh").chmod(0o755)

def make_step(kind: str, command: str, files: dict[str, str]) -> Path:
    dest = Path(STEP[kind]); dest.mkdir(parents=True, exist_ok=True)
    for name, text in files.items(): write(dest / name, text)
    _write_manifest(dest, _load_params())
    render_submit(dest, f"si-eph-{kind}", JZZN_ENVIRONMENT + command)
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
    # Keep the NSCF mesh as an automatic grid so qe2pert can read nk1/nk2/nk3
    # from QE's XML metadata; Wannier90 still receives the explicit list in
    # its .win file below.
    kgrid_text = _fmt_grid(grid(kgrid, f"{calculation.upper()}_KGRID"))
    kpoints = (uniform_kpoints(kgrid) if calculation == "nscf"
               else "K_POINTS automatic\n  " + kgrid_text + " 0 0 0\n")
    nscf_flags = "  nosym = .true.\n  noinv = .true.\n" if calculation == "nscf" else ""
    occupations = str(params.get("OCCUPATIONS", "fixed")).strip().lower()
    smearing = ""
    if occupations not in {"fixed", "none"}:
        smearing = (f"  occupations = '{occupations}'\n"
                    f"  smearing = '{str(params.get('SMEARING', 'gauss')).strip()}'\n"
                    f"  degauss = {float(params.get('DEGAUSS', 0.02)):.8g}\n")
    nbnd_value = int(params.get("NUM_BANDS", 16) if nbnd is None else nbnd)
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
        cell_block = structure_qe(include_cell=False)
    else:
        lattice = "  ibrav = 0\n"
        cell_block = structure_qe(include_cell=True)
    return f"""&CONTROL
  calculation = '{calculation}'
  prefix = '{PREFIX}'
  outdir = './tmp'
  pseudo_dir = './pseudo'
  wf_collect = .true.
/
&SYSTEM
{lattice}  nat = 2
  ntyp = 1
  ecutwfc = {float(params.get('ECUTWFC', 40.0)):.8g}
  ecutrho = {float(params.get('ECUTRHO', 320.0)):.8g}
  nbnd = {nbnd_value}
{smearing}{nscf_flags}/
&ELECTRONS
  conv_thr = {conv_thr:.8e}
  mixing_beta = 0.7
/
{cell_block}{kpoints}"""

def wannier_cell(params: dict[str, object]) -> tuple[str, list[tuple[float, float, float]]]:
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
    return "angstrom", [tuple(v) for v in _poscar()[0]]

def main(kind: str) -> int:
    if kind == "preflight":
        probe = HERE / "run_preflight.py"
        if not probe.is_file():
            raise SystemExit("[错误] 缺 run_preflight.py。")
        make_step(kind, "python3 run_preflight.py",
                  {"run_preflight.py": probe.read_text(encoding="utf-8")})
    elif kind == "scf":
        params = _load_params()
        make_step(kind, "mkdir -p pseudo && ln -sfn \"$PSEUDO_DIR/Si.upf\" pseudo/Si.upf\nrun_mpi pw.x -in scf.in > scf.out\ngrep -q 'JOB DONE' scf.out\nprintf 'OK\\n' > scf.ok", {"scf.in": qe_input("scf", params["SCF_KGRID"], params)})
    elif kind == "wannier":
        params = _load_params()
        nscf_grid = grid(params["NSCF_KGRID"], "NSCF_KGRID")
        kpoints = uniform_kpoints(nscf_grid).splitlines()[2:]
        cell_unit, cell = wannier_cell(params)
        win = f"""num_bands = {int(params['NUM_BANDS'])}
num_wann = {int(params['NUM_WANN'])}
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
Si:sp3
end projections
mp_grid = {nscf_grid[0]} {nscf_grid[1]} {nscf_grid[2]}
begin unit_cell_cart
{cell_unit}
""" + "\n".join("%.10f %.10f %.10f" % tuple(x) for x in cell) + "\nend unit_cell_cart\nbegin atoms_frac\n" + "\n".join("Si %.10f %.10f %.10f" % tuple(x) for x in _poscar()[3]) + "\nend atoms_frac\nbegin kpoints\n" + "\n".join(" ".join(line.split()[:3]) for line in kpoints) + "\nend kpoints\n"
        pw2 = "&inputpp\n  outdir='./tmp'\n  prefix='si'\n  seedname='si'\n  write_mmn=.true.\n  write_amn=.true.\n  write_unk=.false.\n/\n"
        make_step(kind, "ln -sfn ../step1_scf/tmp tmp\nmkdir -p pseudo && ln -sfn \"$PSEUDO_DIR/Si.upf\" pseudo/Si.upf\nrun_mpi pw.x -in nscf.in > nscf.out\n# The site Wannier90 binary is linked against Intel MPI; use its launcher.\nrun_wannier wannier90.x -pp si > wannier_pp.out\nrun_mpi pw2wannier90.x -in pw2wannier.in > pw2wannier.out\nrun_wannier wannier90.x si > wannier.out\nprintf 'OK\\n' > wannier.ok", {"nscf.in": qe_input("nscf", nscf_grid, params), "si.win": win, "pw2wannier.in": pw2})
    elif kind == "phonon":
        params = _load_params()
        qgrid = grid(params["PH_QGRID"], "PH_QGRID")
        phgrid = grid(params["PH_KGRID"], "PH_KGRID")
        ph = """&inputph
  verbosity='debug'
  prefix='si'
  outdir='../step1_scf/tmp'
  fildyn='si.dyn.xml'
  fildvscf='dvscf'
  epsil=.true.
  ldisp=.true.
  lqdir=.true.
  tr2_ph=%.8e
  nq1=%d, nq2=%d, nq3=%d
  nk1=%d, nk2=%d, nk3=%d
 /
""" % (float(params.get("TR2_PH", 1.0e-14)), *qgrid, *phgrid)
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
if [ -s save/si.dyn0 ]; then
  saved_qgrid=$(head -n 1 save/si.dyn0 | xargs)
fi
saved_nq=""
if [ -s save/si.dyn0 ]; then
  saved_nq=$(sed -n '2p' save/si.dyn0 | awk '{print $1}')
fi
expected_nq="__EXPECTED_NQ__"
if [ "$saved_qgrid" = "$expected_qgrid" ] && [ "$saved_nq" = "$expected_nq" ] \
   && [ -d save/si.phsave ] && [ -s save/si.dvscf_q1 ]; then
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
cp -f si.dyn* save/
test -d "$phroot/si.phsave"
cp -a "$phroot/si.phsave" save/

shopt -s nullglob
for qdir in "$phroot"/si.q_*; do
  iq=${qdir##*.q_}
  case "$iq" in (*[!0-9]*|'') exit 1;; esac
  dvscf=("$qdir"/si.dvscf*)
  if [ ${#dvscf[@]} -gt 0 ]; then
    cp -f "${dvscf[0]}" "save/si.dvscf_q${iq}"
  elif [ "$iq" = 1 ]; then
    # QE stores Gamma's response in _ph0/ rather than si.q_1 on this build.
    gamma_dvscf=("$phroot"/si.dvscf*)
    test ${#gamma_dvscf[@]} -gt 0
    cp -f "${gamma_dvscf[0]}" "save/si.dvscf_q${iq}"
  else
    exit 1
  fi
done
# QE may store the Gamma-point response directly in _ph0/ rather than si.q_1.
if [ ! -s save/si.dvscf_q1 ]; then
  gamma_dvscf=("$phroot"/si.dvscf*)
  test ${#gamma_dvscf[@]} -gt 0
  cp -f "${gamma_dvscf[0]}" save/si.dvscf_q1
fi
fi

test -s save/si.dyn0
nq_irr=$(sed -n '2p' save/si.dyn0 | awk '{print $1}')
case "$nq_irr" in (*[!0-9]*|'') exit 1;; esac
test "$nq_irr" -gt 0
for iq in $(seq 1 "$nq_irr"); do
  test -s "save/si.dyn${iq}.xml"
  test -s "save/si.phsave/patterns.${iq}.xml"
  test -s "save/si.dvscf_q${iq}"
done
        printf 'OK\\n' > phonon.ok'''
        collect = collect.replace("__EXPECTED_QGRID__", expected_qgrid)
        collect = collect.replace("__EXPECTED_NQ__", str(expected_nq))
        make_step(kind, collect, {"ph.in": ph})
    elif kind == "qe2pert":
        params = _load_params()
        phgrid = grid(params["PH_QGRID"], "PH_QGRID")
        q2p = grid(params["Q2PERT_KGRID"], "Q2PERT_KGRID")
        if phgrid != q2p:
            raise SystemExit(
                "[错误] qe2pert 要求均匀 PH_QGRID 与 Q2PERT_KGRID 一致，"
                f"收到 PH_QGRID={_fmt_grid(phgrid)}、Q2PERT_KGRID={_fmt_grid(q2p)}"
            )
        inp = """&qe2pert
  prefix='si'
  outdir='../step1_scf/tmp'
  phdir='../step3_phonon/save'
  nk1=%d, nk2=%d, nk3=%d
  dft_band_min=1
  dft_band_max=%d
  num_wann=%d
  lwannier=.true.
/
 """ % (*q2p, int(params["NUM_BANDS"]), int(params["NUM_WANN"]))
        # qe2pert writes prefix_epr.h5 to its QE outdir, which is this step.
        expected_qgrid = " ".join(str(x) for x in phgrid)
        qe2pert_cmd = """ln -sfn ../step2_wannier/si.win si.win
ln -sfn ../step2_wannier/si_u.mat si_u.mat
ln -sfn ../step2_wannier/si_u_dis.mat si_u_dis.mat
ln -sfn ../step2_wannier/si_centres.xyz si_centres.xyz
test -s ../step3_phonon/save/si.dyn0
saved_qgrid=$(head -n 1 ../step3_phonon/save/si.dyn0 | xargs)
test \"$saved_qgrid\" = \"%s\" || { echo \"wrong phonon q grid: got '$saved_qgrid', expected '%s'\" >&2; exit 1; }
run_perturbo_omp qe2pert.x -npools 1 -in qe2pert.in > qe2pert.out
test -s si_epr.h5
printf 'OK\\n' > qe2pert.ok""" % (expected_qgrid, expected_qgrid)
        make_step(kind, qe2pert_cmd, {"qe2pert.in": inp})
    elif kind == "perturbo":
        params = _load_params()
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
        inp = """&perturbo
  prefix='si'
  calc_mode='ephmat'
  fklist='eph.kpt'
  fqlist='eph.qpt'
  band_min=%d
  band_max=%d
  phfreq_cutoff=%.8g
  delta_smear=%.8g
/
 """ % (int(params["PERT_BAND_MIN"]), int(params["PERT_BAND_MAX"]), float(params["PHFREQ_CUTOFF"]), float(params["PERT_DELTA_SMEAR"]))
        make_step(kind, "ln -sfn ../step4_qe2pert/si_epr.h5 si_epr.h5\nrun_perturbo_omp perturbo.x -npools 1 -in pert.in > pert.out\nprintf 'OK\\n' > perturbo.ok", {"pert.in": inp, "eph.kpt": kfile, "eph.qpt": qfile})
    else: raise SystemExit(f"unknown generator kind: {kind}")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
