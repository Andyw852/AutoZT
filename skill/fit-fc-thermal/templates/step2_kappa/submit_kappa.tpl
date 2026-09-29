#!/bin/bash
# fit-fc-thermal S2_kappa job template -- phono3py three-phonon BTE (RTA/LBTE).
# Placeholders: JOBNAME, CONDA_SH, CONDA_ENV, CPUS_PER_TASK, QOS.
# Resource overrides: [submit] in step2_kappa/step.conf, e.g.
#   [submit]
#   cpus_per_task = 96
#   mem           = 128G
#   qos           = premium
#   time          = 48:00:00
#
# phono3py is OpenMP-only (C extension links libgomp, no mpi4py/MPI): the only
# useful single-node parallel is 1 process x N OpenMP threads -> --ntasks=1
# --cpus-per-task=N.  Do NOT use mpirun (it would start N full calculations that
# overwrite each other's output).
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={{CPUS_PER_TASK}}
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos={{QOS}}
cd $SLURM_SUBMIT_DIR

if [ -n "{{CONDA_ENV}}" ] && [ -x "{{CONDA_ENV}}/bin/python" ]; then
    source "{{CONDA_ENV}}/bin/activate"
else
    source {{CONDA_SH}}
    conda activate {{CONDA_ENV}}
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_NESTED=FALSE
export OMP_STACKSIZE=1G GOMP_STACKSIZE=1G
# single process + big OMP: BLAS threads would fight OMP for cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMPY_NUM_THREADS=1

echo "[env] host=$(hostname) env=${CONDA_DEFAULT_ENV:-${VIRTUAL_ENV:-}} threads=$OMP_NUM_THREADS"
python -c "import phono3py; print('phono3py', phono3py.__version__)" || exit 1

set -e
python kappa_driver.py kappa_config.json
