#!/usr/bin/env python3
"""Shared generators for the minimal, reproducible QE/Perturbo Si workflow."""
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
export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:/public/software/aocc/4.1.0/aocc-compiler-4.1.0/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export QE_BIN=/public/home/wangchao/software/AutoZT/qe-7.3-perturbo/bin
export WANNIER_BIN=/public/software/wannier/3.1.0/bin
export PERTURBO_BIN=/public/home/wangchao/software/AutoZT/qe-7.3-perturbo/perturbo/bin
export PSEUDO_DIR=/public/home/wangchao/software/AutoZT/pseudo
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
# AOCL's Zen dot kernel is unstable with this qe2pert build on jzzn.
export BLIS_ARCH_TYPE=generic
export BLIS_NUM_THREADS=1
export PATH=$OPENMPI_BIN:$PATH
for d in "$QE_BIN" "$WANNIER_BIN" "$PERTURBO_BIN"; do
  [ -z "$d" ] || export PATH="$d:$PATH"
done
run_mpi() {
  "$OPENMPI_BIN/mpirun" -np "${SLURM_NTASKS:?SLURM_NTASKS is required}" "$@"
}
run_perturbo_omp() {
  # The tested jzzn qe2pert binary can segfault in its OpenMP zdotc path.
  # Keep the four-core allocation for QE, but run Perturbo serially for the
  # small Si smoke test; larger production runs need a separately validated
  # threaded Perturbo build.
  OMP_NUM_THREADS=1 "$@"
}
run_wannier() {
  "$INTEL_MPI_BIN/mpirun" -np 1 "$@"
}
'''
STEP = {
    "preflight": "step0_preflight", "scf": "step1_scf", "wannier": "step2_wannier",
    "phonon": "step3_phonon", "qe2pert": "step4_qe2pert", "perturbo": "step5_perturbo",
}

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

def structure_qe() -> str:
    cell, species, counts, coords = _poscar()
    if species != ["Si"] or counts != [2]:
        raise SystemExit("[错误] v0.1 的端到端验证仅支持两原子金刚石 Si POSCAR；通用元素映射将在后续版本实现。")
    text = "CELL_PARAMETERS angstrom\n" + "\n".join("  %.10f %.10f %.10f" % tuple(v) for v in cell)
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
    render_submit(dest, f"si-eph-{kind}", JZZN_ENVIRONMENT + command)
    print(f"[DONE] {dest} 输入已生成；待 AutoZT 提交。")
    return dest

def uniform_kpoints(n: int) -> str:
    """Return the full uniform grid required by the Wannier90 interface."""
    points = [
        (i / n, j / n, k / n, 1.0 / (n ** 3))
        for i in range(n) for j in range(n) for k in range(n)
    ]
    return "K_POINTS crystal\n%d\n" % len(points) + "\n".join(
        "  %.10f %.10f %.10f %.12e" % point for point in points
    ) + "\n"

def qe_input(calculation: str, kgrid: str, nbnd: int = 16) -> str:
    # Keep the NSCF mesh as an automatic grid so qe2pert can read nk1/nk2/nk3
    # from QE's XML metadata; Wannier90 still receives the explicit list in
    # its .win file below.
    kpoints = (uniform_kpoints(int(kgrid.split()[0])) if calculation == "nscf"
               else "K_POINTS automatic\n  " + kgrid + " 0 0 0\n")
    nscf_flags = "  nosym = .true.\n  noinv = .true.\n" if calculation == "nscf" else ""
    return f"""&CONTROL
  calculation = '{calculation}'
  prefix = '{PREFIX}'
  outdir = './tmp'
  pseudo_dir = './pseudo'
  wf_collect = .true.
/
&SYSTEM
  ibrav = 0
  nat = 2
  ntyp = 1
  ecutwfc = 40.0
  ecutrho = 320.0
  nbnd = {nbnd}
{nscf_flags}/
&ELECTRONS
  conv_thr = 1.0d-10
  mixing_beta = 0.7
/
{structure_qe()}{kpoints}"""

def main(kind: str) -> int:
    if kind == "preflight":
        probe = HERE / "run_preflight.py"
        if not probe.is_file():
            raise SystemExit("[错误] 缺 run_preflight.py。")
        make_step(kind, "python3 run_preflight.py",
                  {"run_preflight.py": probe.read_text(encoding="utf-8")})
    elif kind == "scf":
        make_step(kind, "mkdir -p pseudo && ln -sfn \"$PSEUDO_DIR/Si.upf\" pseudo/Si.upf\nrun_mpi pw.x -in scf.in > scf.out\ngrep -q 'JOB DONE' scf.out\nprintf 'OK\\n' > scf.ok", {"scf.in": qe_input("scf", "4 4 4")})
    elif kind == "wannier":
        kpoints = uniform_kpoints(4).splitlines()[2:]
        win = """num_bands = 16
num_wann = 8
write_u_matrices = true
write_hr = true
write_xyz = true
begin projections
Si:sp3
end projections
mp_grid = 4 4 4
begin unit_cell_cart
angstrom
""" + "\n".join("%.10f %.10f %.10f" % tuple(x) for x in _poscar()[0]) + "\nend unit_cell_cart\nbegin atoms_frac\n" + "\n".join("Si %.10f %.10f %.10f" % tuple(x) for x in _poscar()[3]) + "\nend atoms_frac\nbegin kpoints\n" + "\n".join(" ".join(line.split()[:3]) for line in kpoints) + "\nend kpoints\n"
        pw2 = "&inputpp\n  outdir='./tmp'\n  prefix='si'\n  seedname='si'\n  write_mmn=.true.\n  write_amn=.true.\n  write_unk=.false.\n/\n"
        make_step(kind, "ln -sfn ../step1_scf/tmp tmp\nmkdir -p pseudo && ln -sfn \"$PSEUDO_DIR/Si.upf\" pseudo/Si.upf\nrun_mpi pw.x -in nscf.in > nscf.out\n# The site Wannier90 binary is linked against Intel MPI; use its launcher.\nrun_wannier wannier90.x -pp si > wannier_pp.out\nrun_mpi pw2wannier90.x -in pw2wannier.in > pw2wannier.out\nrun_wannier wannier90.x si > wannier.out\nprintf 'OK\\n' > wannier.ok", {"nscf.in": qe_input("nscf", "4 4 4"), "si.win": win, "pw2wannier.in": pw2})
    elif kind == "phonon":
        ph = """&inputph
  prefix='si'
  outdir='../step1_scf/tmp'
  fildyn='si.dyn.xml'
  fildvscf='dvscf'
  ldisp=.true.
  lqdir=.true.
  search_sym=.false.
  tr2_ph=1.0d-14
  nq1=2, nq2=2, nq3=2
  nk1=4, nk2=4, nk3=4
/
"""
        # qe2pert reads a flattened phonon directory, not QE's _ph0 tree.
        # Preserve the QE-generated numbering so each irreducible q point has
        # its matching dynmat, displacement pattern, and dvscf perturbation.
        collect = r'''# Archive only QE's phonon scratch tree for this step.  A stale _ph0
# can contain dvscf files from a different FFT grid; keep it recoverable while
# ensuring ph.x starts with a clean tree.
if [ -d ../step1_scf/tmp/_ph0 ]; then
  mv ../step1_scf/tmp/_ph0 ../step1_scf/tmp/_ph0.stale.${SLURM_JOB_ID:-manual}
fi
run_mpi ph.x -in ph.in > ph.out
phroot='../step1_scf/tmp/_ph0'
mkdir -p save
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

test -s save/si.dyn0
nq_irr=$(sed -n '2p' save/si.dyn0 | awk '{print $1}')
case "$nq_irr" in (*[!0-9]*|'') exit 1;; esac
test "$nq_irr" -gt 0
for iq in $(seq 1 "$nq_irr"); do
  test -s "save/si.dyn${iq}"
  test -s "save/si.phsave/patterns.${iq}.xml"
  test -s "save/si.dvscf_q${iq}"
done
printf 'OK\\n' > phonon.ok'''
        make_step(kind, collect, {"ph.in": ph})
    elif kind == "qe2pert":
        inp = """&qe2pert
  prefix='si'
  outdir='../step1_scf/tmp'
  phdir='../step3_phonon/save'
  nk1=4, nk2=4, nk3=4
  dft_band_min=1
  dft_band_max=16
  num_wann=8
  lwannier=.true.
/
"""
        # qe2pert writes prefix_epr.h5 to its QE outdir, not its CWD.  Keep a
        # step-local link so AutoZT can fetch and validate the declared output.
        qe2pert_cmd = "ln -sfn ../step2_wannier/si.win si.win\nln -sfn ../step2_wannier/si_u.mat si_u.mat\nln -sfn ../step2_wannier/si_u_dis.mat si_u_dis.mat\nln -sfn ../step2_wannier/si_centres.xyz si_centres.xyz\nrun_perturbo_omp qe2pert.x -npools 1 -in qe2pert.in > qe2pert.out\ntest -s ../step1_scf/tmp/si_epr.h5\nln -sfn ../step1_scf/tmp/si_epr.h5 si_epr.h5\ntest -s si_epr.h5\nprintf 'OK\\n' > qe2pert.ok"
        make_step(kind, qe2pert_cmd, {"qe2pert.in": inp})
    elif kind == "perturbo":
        inp = """&perturbo
  prefix='si'
  calc_mode='ephmat'
  fklist='eph.kpt'
  fqlist='eph.qpt'
  band_min=1
  band_max=1
  phfreq_cutoff=1
/
"""
        make_step(kind, "ln -sfn ../step4_qe2pert/si_epr.h5 si_epr.h5\nrun_perturbo_omp perturbo.x -npools 1 -in pert.in > pert.out\nprintf 'OK\\n' > perturbo.ok", {"pert.in": inp, "eph.kpt": "1\n0.0 0.0 0.0\n", "eph.qpt": "1\n0.0 0.0 0.0\n"})
    else: raise SystemExit(f"unknown generator kind: {kind}")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
