#!/bin/bash
# =====================================================================
# submit_nep.tpl —— NEP 数据转换 + 训练提交模板（kl-gpumd-3090 / step2_nep）
# 占位符：JOBNAME / CONDA_SH / CONDA_ENV / NEP_CMD / LOG
#
# GPUMD v5.6 的 nep 是 CPU/OpenMP 实现（CUDA 版要单独编译），默认走 cpu192；
# 换队列/核数用 step.conf 的 [submit] 段覆盖，不要改模板。
# =====================================================================
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=regular

cd ${SLURM_SUBMIT_DIR}

if [ -f "{{CONDA_ENV}}/bin/activate" ]; then
  source "{{CONDA_ENV}}/bin/activate"
else
  source {{CONDA_SH}}
  conda activate {{CONDA_ENV}}
fi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "[env] host=$(hostname) threads=$OMP_NUM_THREADS"
{{NEP_CMD}} 2>&1 | tee -a {{LOG}}
echo "NEP_EXIT=${PIPESTATUS[0]}"
tail -20 {{LOG}} || true
