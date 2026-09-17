#!/bin/bash
# 声子谱绘图作业（step5_phonon_plot）。占位符（由 gen 脚本替换）：
#   {{JOBNAME}} {{QOS}} {{NTASKS}} {{CONDA_SH}} {{CONDA_ENV}} {{PLOT_CMD}}
#
# 为什么这一步是【作业】而不是登录节点 gen（2026-09-17 从 taskflow-v2.0 同步）：
#   128 原子胞（512 原子超胞）的一条完整 band structure 在登录节点要几分钟到十几分钟，
#   超过 autozt 给 gen 的预算，会被中途杀掉（实测连 band_<tag>.yaml 都写不完 → 永远失败）。
#   放计算节点没有这个限制，且脚本支持从已算好的 band_<tag>.yaml 断点续画（重跑几秒）。
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos={{QOS}}
cd $SLURM_SUBMIT_DIR
source {{CONDA_SH}}
conda activate {{CONDA_ENV}}
export OMP_NUM_THREADS=1
echo "[plot] start $(date) host=$(hostname) cores=$SLURM_NTASKS_PER_NODE"
{{PLOT_CMD}}
_rc=$?
echo "[plot] done $(date) rc=$_rc"
# 必须把工作进程的退出码传出去：否则 sbatch 会把失败当成功（曾经这样"18 秒完成"过）
exit $_rc
