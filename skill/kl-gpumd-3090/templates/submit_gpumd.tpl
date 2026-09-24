#!/bin/bash
# =====================================================================
# submit_gpumd.tpl —— GPUMD MD-κ 提交模板（kl-gpumd-3090 / step3_kappa）
# 占位符：JOBNAME / CONDA_SH / CONDA_ENV / GPUMD_CMD / LOG
#
# 注意：3090 机上 ~/gpumd/src/gpumd 是 **CUDA 版**（启动日志会列 GPU），
# 所以这里默认投 gpu 分区、单卡；它仍然需要几个 CPU 核做邻居表和输入输出。
# 若换成 CPU 版（make -f makefile / makefile.gcc），用 step.conf 的 [submit] 段
# 把 partition 改回 cpu192 并去掉 gres 即可，模板不用动。
# =====================================================================
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=queue.out
#SBATCH --error=queue.err

cd ${SLURM_SUBMIT_DIR}

if [ -f "{{CONDA_ENV}}/bin/activate" ]; then
  source "{{CONDA_ENV}}/bin/activate"
else
  source {{CONDA_SH}}
  conda activate {{CONDA_ENV}}
fi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "[env] host=$(hostname) threads=$OMP_NUM_THREADS cuda=${CUDA_VISIBLE_DEVICES:-none}"
nvidia-smi --query-gpu=index,memory.free --format=csv,noheader || true
{{GPUMD_CMD}} 2>&1 | tee {{LOG}}
echo "GPUMD_EXIT=${PIPESTATUS[0]}"

