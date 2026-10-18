#!/bin/bash
# eph-qe-cpu [0059] ivs_dfpt（S11_ivsscf / S12_ivsph）的通用提交模板：不写集群专用的分区、模块和路径。
#   核数：step.conf 的 IVS_NP（按 --ntasks 申请；无 SLURM 的 fakeslurm 垫片同样按它限核）。
#   额外 #SBATCH 行：step.conf 的 IVS_SBATCH（分号分隔，如 jzzn：--partition=cpu192; --qos=regular）。
#   QE 路径与环境：setting/<集群>.yaml 的 qe_bin / qe_pseudo_dir / qe_setup（autozt 注入 step.conf 的
#   QE_BIN / QE_PSEUDO_DIR / QE_SETUP，由 gen 写进下面的命令段）。
#   集群确要专用版本：放 setting/<集群>/templates/submit_ivs_qe.tpl（同名覆盖，保留 JOBNAME/NP/SBATCH_EXTRA/COMMAND
#   四个双花括号占位符）。
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks={{NP}}
#SBATCH --output=queue.out
#SBATCH --error=queue.err
{{SBATCH_EXTRA}}

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

{{COMMAND}}
