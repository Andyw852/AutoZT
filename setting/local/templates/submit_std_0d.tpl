#!/bin/bash
# =====================================================================
# submit_std_0d.tpl —— AutoZT 本机后端（setting/local.yaml）通用提交模板
# 由 tools/fakeslurm（或本机真 SLURM）执行；分区/qos 对 fakeslurm 无意义，不写。
# 核数：step.conf 的 [submit] cpus_per_task 覆盖；fakeslurm 按 FAKESLURM_MAX_CPUS 限流。
# 本文件不含站点路径，可随仓库分发；某台机器要特殊处理就复制到项目
# project_setting/templates/ 下改（项目级优先于这里）。
# =====================================================================
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=queue.out
#SBATCH --error=queue.err
cd "${SLURM_SUBMIT_DIR:-$PWD}"
NCPU=${SLURM_CPUS_PER_TASK:-4}
# 没装 MPI 的机器：mpirun -np N cmd → 直接跑 cmd（单进程）
command -v mpirun >/dev/null 2>&1 || mpirun() { while [ $# -gt 0 ] && [[ "$1" == -* ]]; do case "$1" in -np|-n) shift 2;; *) shift;; esac; done; "$@"; }
# 只有 python3 的系统（Debian/Ubuntu）：gen 写的命令里用的是 python
command -v python >/dev/null 2>&1 || python() { python3 "$@"; }
ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=1
# 并行：本机按 MPI 进程数 = cpus_per_task（无 MPI 时单进程）
mpirun -np ${NCPU} vasp_std
