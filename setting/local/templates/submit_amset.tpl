#!/bin/bash
# =====================================================================
# submit_amset.tpl —— AutoZT 本机后端（setting/local.yaml）通用提交模板
# 由 tools/fakeslurm（或本机真 SLURM）执行；分区/qos 对 fakeslurm 无意义，不写。
# 核数：step.conf 的 [submit] cpus_per_task 覆盖；fakeslurm 按 FAKESLURM_MAX_CPUS 限流。
# 本文件不含站点路径，可随仓库分发；某台机器要特殊处理就复制到项目
# project_setting/templates/ 下改（项目级优先于这里）。
# =====================================================================
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=queue.out
#SBATCH --error=queue.err
cd "${SLURM_SUBMIT_DIR:-$PWD}"
NCPU=${SLURM_CPUS_PER_TASK:-8}
# 没装 MPI 的机器：mpirun -np N cmd → 直接跑 cmd（单进程）
command -v mpirun >/dev/null 2>&1 || mpirun() { while [ $# -gt 0 ] && [[ "$1" == -* ]]; do case "$1" in -np|-n) shift 2;; *) shift;; esac; done; "$@"; }
# 只有 python3 的系统（Debian/Ubuntu）：gen 写的命令里用的是 python
command -v python >/dev/null 2>&1 || python() { python3 "$@"; }
# AMSET_ENV：venv 目录或 conda 环境名（setting/local.yaml 的 amset_env）
AMSET_ENV="{{AMSET_ENV}}"
if [ -x "${AMSET_ENV}/bin/python" ]; then
    source "${AMSET_ENV}/bin/activate"
else
    for _sh in "${CONDA_SH_LOCAL:-}" "$HOME/miniconda3/etc/profile.d/conda.sh" \
               "$HOME/anaconda3/etc/profile.d/conda.sh" "$HOME/miniforge3/etc/profile.d/conda.sh"; do
        [ -n "$_sh" ] && [ -f "$_sh" ] && { source "$_sh"; break; }
    done
    if command -v conda >/dev/null 2>&1; then
        conda activate "${AMSET_ENV}" || echo "[WARN] conda env ${AMSET_ENV} not found -- using PATH" >&2
    else
        echo "[WARN] no conda -- using amset from PATH" >&2
    fi
fi
export OMP_NUM_THREADS=${NCPU} OPENBLAS_NUM_THREADS=${NCPU} MKL_NUM_THREADS=${NCPU}
{{AMSET_CMD}}
