#!/bin/bash
# =====================================================================
# submit_std_0d.tpl —— vasp_std 提交模板（0D 孤立分子/团簇，step1/step2）
# 逻辑名 submit_std_0d.tpl，实际文件名 submit_jzzn_vaspstd_0d.tpl，
# 映射写在 hpc.yaml 的 template_map 里；换机器只改本文件。
#
# 与 3D 版的区别：分子只有 Γ 点、原子数少，96 核往往跑不满，
# 反而通信开销占比上升。这里默认半个节点；71 个体系批量跑时
# 用更少的核换更高的吞吐通常更划算。按实际 scaling 调。
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

# 只有一个不可约 k 点：INCAR 里 KPAR=1，别开 k 点并行。
# 想再省一半内存/时间可以换 vasp_gam（Γ-only 专用二进制），
# 但要确认 KPOINTS 确实是 1x1x1，且 relax_common 的 extract_vasp_cmd
# 认得这一行（正则匹配 vasp_(std|ncl|gam)，vasp_gam 也认）。
export AUTOZT_CELL_CONSTRAINT=none
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_std
