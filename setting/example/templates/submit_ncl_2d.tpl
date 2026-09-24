#!/bin/bash
# =====================================================================
# submit_ncl_2d.tpl —— 【示例/脱敏】非共线 2D VASP 提交模板
# 用法：复制为 setting/<集群名>/templates/submit_ncl_2d.tpl 后，把下方占位符换成真实值。
# 本文件不含任何站点路径/账号/分区信息，可随仓库分发。
# =====================================================================
#SBATCH --partition=<your-partition>
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=<your-cores>
#SBATCH --output=queue.out
#SBATCH --error=queue-%j.err
#SBATCH --qos=<your-qos>
cd $SLURM_SUBMIT_DIR
# 加载编译器/MPI 环境（示例，按站点改）：module load <compiler> <mpi>
ulimit -s unlimited
export OMP_NUM_THREADS=1
mpirun -np <your-cores> /path/to/vasp
