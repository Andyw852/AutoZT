#!/bin/bash
# fit-fc-thermal S1_fit job template -- FIT_ENGINE=pheasy (CPU build).
# Placeholders: JOBNAME, CONDA_SH, CONDA_ENV, ENGINE.
# pheasy is core- and memory-hungry (cluster space + sensing matrix + fit);
# override via step1_fit/step.conf:
#   [submit]
#   cpus_per_task = 96
#   mem           = 256G
#   time          = 48:00:00
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=regular
cd $SLURM_SUBMIT_DIR

# 环境：CONDA_ENV 是 venv 目录 → 直接激活；否则 source CONDA_SH 再 conda activate。
# 两者都留空（如 setting/local.yaml）= 用 PATH 里现成的 python，不激活任何环境。
if [ -n "{{CONDA_ENV}}" ] && [ -x "{{CONDA_ENV}}/bin/python" ]; then
    source "{{CONDA_ENV}}/bin/activate"
elif [ -n "{{CONDA_SH}}" ]; then
    if [ -f {{CONDA_SH}} ]; then
        source {{CONDA_SH}}
        if [ -n "{{CONDA_ENV}}" ]; then
            conda activate {{CONDA_ENV}}
        fi
    else
        echo "[WARN] conda.sh not found: {{CONDA_SH}} -- using python from PATH" >&2
    fi
fi

NCPU=${SLURM_CPUS_PER_TASK:-48}
export OMP_NUM_THREADS=${NCPU}
export OPENBLAS_NUM_THREADS=${NCPU}
export MKL_NUM_THREADS=${NCPU}

echo "[env] host=$(hostname) env=${CONDA_DEFAULT_ENV:-${VIRTUAL_ENV:-}} threads=${NCPU}"
for _m in numpy scipy phonopy spglib h5py; do
    python -c "import ${_m}" 2>/dev/null && echo "  ok  ${_m}" \
        || echo "  MISSING ${_m} (required by fit-fc-thermal)"
done
python -c "import pheasy" 2>/dev/null && echo "  ok  pheasy" \
    || echo "  MISSING pheasy in this environment"
python -c "from celer import Lasso" 2>/dev/null && echo "  ok  celer (LASSO/ALASSO/RFE)" \
    || echo "  note: celer absent -- PHEASY_FIT_METHOD must be OLS/RIDGE"

set -e
# prep also writes disp_matrix.pkl / force_matrix.pkl, which is what pheasy
# reads through --disp_file.
python fc_fit_driver.py prep fit_config.json
python fc_fit_driver.py fit  fit_config.json
python fc_fit_driver.py post fit_config.json
