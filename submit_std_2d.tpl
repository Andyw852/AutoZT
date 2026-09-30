#!/bin/bash
# =====================================================================
# submit_std_2d.tpl —— vasp_std 提交模板（2D，step1/step2；optcell 补丁版）
# 集群参数在此模板中直接写死；换机器/换队列/换 VASP 时只改本文件。
# 脚本只填充 JOBNAME 占位符（按体系与步骤自动生成，可用 --jobname 覆盖）。
# 使用 optcell 补丁版：弛豫步若放置 OPTCELL 文件可做约束变胞，
# 无 OPTCELL 文件时行为与官方版完全一致。
# =====================================================================
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=queue.out
#SBATCH --error=queue-%j.err
#SBATCH --qos=regular
# --time 不能超过 QOS 的 MaxWall，超了作业会永远 PD 在
# QOSMaxWallDurationPerJobLimit。先查上限再改：
#     sacctmgr show qos format=Name,MaxWall,MaxTRESPerJob
# 不写这一行 = 直接吃 QOS 上限（jzzn 的 regular 是 24h）。
#SBATCH --time=24:00:00
cd $SLURM_SUBMIT_DIR
# 加载与编译时一致的 AOCC 环境（不要用 Intel）
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1
# 走 UCX，禁用已废弃的 openib BTL
export OMPI_MCA_pml=ucx
export OMPI_MCA_btl=^openib,uct
export UCX_TLS=rc,sm,self
export UCX_NET_DEVICES=mlx5_0:1
# 标准共线版本：step1/step2 使用 vasp_std
export AUTOZT_CELL_CONSTRAINT=ioptcell_tag
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3-optcell/bin/vasp_std
