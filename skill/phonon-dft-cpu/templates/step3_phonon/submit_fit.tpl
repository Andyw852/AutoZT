#!/bin/bash
# phonon-dft-cpu S3 拟合作业模板（phonopy 收力 + fc2 + 声子谱）。占位符 {{JOBNAME}}。
# 集群 conda 走 step.conf 注入的 {{CONDA_SH}}/{{CONDA_ENV}}，不写死个人路径/环境名。
# 资源可用 step.conf 的 [submit] 段覆盖：cpus_per_task / qos / partition / time。
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=regular
cd $SLURM_SUBMIT_DIR
source {{CONDA_SH}}
conda activate {{CONDA_ENV}}

NCPU=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=${NCPU}
export OPENBLAS_NUM_THREADS=${NCPU}
export MKL_NUM_THREADS=${NCPU}
echo "[env] $CONDA_DEFAULT_ENV | phonopy=$(python -c 'import phonopy;print(phonopy.__version__)' 2>&1)"

# ★ 不用管道吞退出码：先落日志再显式 exit $_rc，作业失败必须让 SLURM 看到非 0。
python -u phonon_fit_driver.py > phonon_fit.log 2>&1
_rc=$?
tail -40 phonon_fit.log
exit $_rc
