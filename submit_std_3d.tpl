#!/bin/bash
# =====================================================================
# submit_std_3d.tpl —— vasp_std 提交模板（3D，step1/step2；官方 VASP）
# 集群参数在此模板中直接写死；换机器/换队列/换 VASP 时只改本文件。
# 脚本只填充 JOBNAME 占位符（按体系与步骤自动生成，可用 --jobname 覆盖）。
# =====================================================================
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=queue.out
#SBATCH --error=queue-%j.err
#SBATCH --qos=regular

cd $SLURM_SUBMIT_DIR

# 加载与编译时一致的 AOCC 环境（不要用 Intel）
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited

# 标准共线版本：step1/step2 使用 vasp_std
export AUTOZT_CELL_CONSTRAINT=none
mpirun -np $SLURM_NTASKS /public/home/<user>/software/vasp.6.4.3/bin/vasp_std
