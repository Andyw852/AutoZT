#!/bin/bash
# Four-core, isolated smoke-test template for eph-qe-cpu.  It never changes
# the site-wide VASP templates.
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=regular

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# This template is self-contained for the private QE 7.3/Perturbo install.
# The generator repeats the environment so generated submit.sh files remain
# valid if a project has an older copied template.
module purge
module load gcc/10.2.0 openmpi/4.0.5 hdf5/1.14.6 aocl/4.1.0 fftw/3.3.8 wannier/3.1.0
export OPENMPI_BIN=/public/software/openmpi/4.0.5/bin
export LD_LIBRARY_PATH=/public/software/aocc/4.1.0/aocc-compiler-4.1.0/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export QE_BIN=/public/home/wangchao/software/AutoZT/qe-7.3-perturbo/bin
export WANNIER_BIN=/public/software/wannier/3.1.0/bin
export PERTURBO_BIN=/public/home/wangchao/software/AutoZT/qe-7.3-perturbo/perturbo/bin
export PSEUDO_DIR=/public/home/wangchao/software/AutoZT/pseudo
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
export BLIS_ARCH_TYPE=generic
export BLIS_NUM_THREADS=1
export PATH="$OPENMPI_BIN:$QE_BIN:$WANNIER_BIN:$PERTURBO_BIN:$PATH"
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

{{COMMAND}}
