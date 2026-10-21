#!/usr/bin/env bash
# =============================================================================
# probe_cluster.sh —— 只读侦察一个超算节点/账号，给 AutoZT 接入（setting/<hpc>.yaml）
# 收集事实。不提交作业、不写远端文件（--test-only 除外：仅向调度器做"空提交预检"）。
#
# 用法：
#   tools/probe_cluster.sh                       # 本机
#   tools/probe_cluster.sh <ssh别名>             # 远端（如 jzzn / eshell111）
#   tools/probe_cluster.sh <ssh别名> --test-only # 额外逐分区做 sbatch --test-only 可提交性探测
#   tools/probe_cluster.sh <ssh别名> --test-only --max-partitions 8
#
# 说明：
#   * 对"有没有 SLURM / 有没有 module 系统 / 有没有 GPU"都不做假设，逐项探测；
#     没有的命令直接跳过并在输出里写明，不会报错退出。
#   * 输出给人看，也便于直接抄进 setting/<hpc>.yaml；末尾给填写提示。
# =============================================================================
set -u
set -o pipefail

TARGET="${1:-local}"
[ $# -gt 0 ] && shift
TEST_ONLY=0
MAX_PART=24
while [ $# -gt 0 ]; do
  case "$1" in
    --test-only) TEST_ONLY=1 ;;
    --max-partitions) MAX_PART="${2:-24}"; shift ;;
    -h|--help) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "未知参数: $1" >&2; exit 2 ;;
  esac
  shift
done

if [ "$TARGET" = "local" ]; then
  SSH=""
else
  SSH="ssh -o BatchMode=yes -o ConnectTimeout=15"
  if ! $SSH "$TARGET" true 2>/dev/null; then
    echo "[ERROR] 连不上 ssh 目标 '$TARGET'（BatchMode）。先配好 ~/.ssh/config 与密钥。" >&2
    exit 1
  fi
fi

RUN() {
  # 用 bash -s 经 stdin 执行：不依赖登录 shell 类型、也不重复执行（否则同一段输出会打两遍）
  if [ -z "$SSH" ]; then printf '%s\n' "$1" | bash -s; else printf '%s\n' "$1" | $SSH "$TARGET" bash -s; fi
}

echo "======================================================================"
echo " AutoZT 集群侦察报告 —— target=$TARGET   $(date '+%F %T')"
echo "======================================================================"

echo
echo "## 0. 基本信息"
RUN '
echo "hostname      : $(hostname 2>/dev/null)"
echo "user          : $(whoami 2>/dev/null)"
echo "login shell   : $(getent passwd $(whoami) 2>/dev/null | cut -d: -f7)"
echo "os            : $(grep PRETTY_NAME /etc/os-release 2>/dev/null | cut -d= -f2- | tr -d "\"")"
echo "kernel        : $(uname -r 2>/dev/null)  arch: $(uname -m 2>/dev/null)"
echo "cpus(login)   : $(nproc 2>/dev/null)"
echo "mem(login)    : $(free -g 2>/dev/null | awk "/^Mem:/{print \$2\" GB\"}")"
echo "home          : $HOME"
'

echo
echo "## 1. 批处理调度器"
RUN '
for c in sbatch squeue sinfo sacct scontrol scancel; do command -v $c >/dev/null 2>&1 && echo "slurm  : $c -> $(command -v $c)"; done
command -v qsub   >/dev/null 2>&1 && echo "PBS/SGE: qsub -> $(command -v qsub)"
command -v bsub   >/dev/null 2>&1 && echo "LSF    : bsub -> $(command -v bsub)"
command -v sinfo >/dev/null 2>&1 && echo "slurm 版本: $(sinfo -V 2>/dev/null | head -1)"
if command -v sbatch >/dev/null 2>&1; then
  echo "SLURM 可用（AutoZT 直接用 start 提交）"
elif command -v qsub >/dev/null 2>&1; then
  echo "非 SLURM 批处理（AutoZT 需自备 submit 模板或适配层，见文档第 9 节）"
else
  echo "无批处理调度器：AutoZT 可用 local 模式 / tools/fakeslurm（AGENTS 铁律 1）"
fi
'

echo
echo "## 2. SLURM 分区（AV|分区|时限|节点数|核/节点|内存MB|GRES）"
RUN '
if command -v sinfo >/dev/null 2>&1; then
  sinfo -a -h -o "%P|%a|%l|%D|%c|%m|%G" 2>/dev/null | awk -F"|" "
    { key=\$1; cnt[key]+=\$4; cpu[key]=\$5; mem[key]=\$6; gres[key]=\$7; avail[key]=\$2; tl[key]=\$3 }
    END { for (k in cnt) printf \"%-22s %-4s %-12s nodes=%-5d %s-core %sMB %s\n\", k, avail[k], tl[k], cnt[k], cpu[k], mem[k], gres[k] }" | sort
  echo
  echo "分区状态汇总:"
  sinfo -a -h -o "%t" 2>/dev/null | sed "s/[^a-z].*//" | sort | uniq -c | head -8
  echo
  echo "节点状态样例（前 5）:"
  sinfo -a -h -o "%N %P %t %c %G" 2>/dev/null | head -5
else
  echo "（无 sinfo：非 SLURM 集群）"
fi
'

echo
echo "## 3. 账户关联 / 配额（能读到多少算多少）"
RUN '
if command -v sacctmgr >/dev/null 2>&1; then
  echo "--- assoc（account/user/partition/QOS/GrpTRES/MaxTRES/MaxJobs）---"
  sacctmgr -n -P show assoc where user=$USER format=Account,User,Partition,QOS,GrpTRES,MaxTRES,MaxJobs 2>/dev/null | head -5
  echo "--- QOS ---"
  sacctmgr -n -P show qos where name="user_$USER" 2>/dev/null | head -3
fi
if command -v sshare >/dev/null 2>&1; then
  echo "--- sshare ---"; sshare -P -u $USER 2>/dev/null | head -3
fi
echo "提示：GrpTRES 里的 cpu= 就是账户并发上限；单作业上限可用 --test-only 逐档试（第 10 节）"
'

echo
echo "## 4. 存储"
RUN '
echo "--- df -h ---"
df -h $HOME 2>/dev/null | tail -2
for d in /work /public /scratch /gpfs /lustre; do [ -d "$d" ] && df -h "$d" 2>/dev/null | tail -1; done
echo "--- quota（若有）---"
( quota -s 2>/dev/null || lfs quota -h $HOME 2>/dev/null ) | head -6
echo "--- home 写权限 ---"
[ -w "$HOME" ] && echo "能写 $HOME" || echo "不能写 $HOME"
'

echo
echo "## 5. 模块系统 / 编译器 / MPI"
RUN '
if command -v module >/dev/null 2>&1; then
  echo "module 系统: $(module --version 2>&1 | head -2 | tr "\n" " ")"
  echo "可加载模块数: $(module avail 2>&1 | wc -l)"
  echo "--- 编译器/MPI/数学库相关（节选）---"
  module avail 2>&1 | grep -oE "(compiler|mathlib|mpi)/[^ ]+" | sort -u | head -25
else
  echo "（无 module 系统：编译器/库来自 PATH 或 /opt）"
fi
for c in gcc g++ gfortran icc icpc ifort ifx mpiicc mpicc mpif90 mpirun mpiexec cmake make; do
  p=$(command -v $c 2>/dev/null) && echo "$c -> $p ($($c --version 2>&1 | head -1 | cut -c1-60))"
done
'

echo
echo "## 6. 数学库 / FFT（VASP 等源码构建要用）"
RUN '
echo "MKLROOT=${MKLROOT:-（未设）}"
for d in /opt/intel /opt/hpc/software/compiler/intel /public/software/mathlib; do
  [ -d "$d" ] && echo "目录 $d: $(ls "$d" 2>/dev/null | head -6 | tr "\n" " ")"
done
echo "--- 动态库可见性 ---"
( ldconfig -p 2>/dev/null | grep -E "libmkl_core|libopenblas|liblapack|libblas|libfftw3|libscalapack|libsymspg" | head -12 ) || true
'

echo
echo "## 7. 加速卡"
RUN '
if command -v nvidia-smi >/dev/null 2>&1; then echo "NVIDIA: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -4 | tr "\n" "; ")"; else echo "NVIDIA: 无 nvidia-smi"; fi
if [ -x /opt/hyhal/bin/hy-smi ]; then echo "Hygon DCU: /opt/hyhal/bin/hy-smi 可用"; /opt/hyhal/bin/hy-smi -L 2>/dev/null | head -4; fi
command -v rocm-smi >/dev/null 2>&1 && echo "AMD ROCm: $(command -v rocm-smi)"
[ -e /dev/kfd ] && echo "/dev/kfd 存在（ROCm/DCU 类）"
'

echo
echo "## 8. Python 环境（AutoZT 采集要求【非交互】PATH 里有 python3）"
RUN '
echo "python3 -> $(command -v python3 2>/dev/null || echo 未找到)"
python3 -V 2>/dev/null
python3 -c "import numpy;print(\"numpy\",numpy.__version__)" 2>/dev/null || echo "numpy: 无"
python3 -c "import scipy;print(\"scipy\",scipy.__version__)" 2>/dev/null || true
python3 -c "import pymatgen;print(\"pymatgen ok\")" 2>/dev/null || true
python3 -c "import spglib;print(\"spglib ok\")" 2>/dev/null || true
python3 -c "import h5py;print(\"h5py ok\")" 2>/dev/null || true
for c in conda mamba micromamba; do command -v $c >/dev/null 2>&1 && echo "$c -> $(command -v $c)"; done
echo "--- PATH 里第一个可写目录（可放 python3 软链）---"
IFS=":"; for d in $PATH; do [ -d "$d" ] && [ -w "$d" ] && { echo "$d"; break; }; done; unset IFS
'

echo
echo "## 9. AutoZT 现状（是否已部署过）"
RUN '
if [ -d "$HOME/software/AutoZT" ]; then
  echo "已存在 ~/software/AutoZT："
  du -sh "$HOME"/software/AutoZT/* 2>/dev/null | tail -12
else
  echo "未发现 ~/software/AutoZT（按 AGENTS 铁律 11 新建）"
fi
echo "--- ~/.vaspkit ---"; [ -f "$HOME/.vaspkit" ] && head -4 "$HOME/.vaspkit" || echo "无"
echo "--- POTCAR 候选目录（前 3）---"
find "$HOME/software" -maxdepth 5 -type d -name "potpaw*" 2>/dev/null | head -3
'

if [ "$TEST_ONLY" = "1" ]; then
  echo
  echo "## 10. 可提交性预检（sbatch --test-only，不占资源；每分区 1 核/10 分钟）"
  RUN "
if ! command -v sbatch >/dev/null 2>&1; then echo '（非 SLURM，跳过）'; exit 0; fi
tmp=\$(mktemp -d 2>/dev/null) || tmp=\$HOME
i=0
for p in \$(sinfo -a -h -o '%P' 2>/dev/null | tr -d '*' | sort -u); do
  i=\$((i+1)); [ \$i -gt $MAX_PART ] && break
  f=\$tmp/probe.slurm
  printf '#!/bin/bash\n#SBATCH -p %s\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -t 00:10:00\nhostname\n' \"\$p\" > \"\$f\"
  out=\$(sbatch --test-only \"\$f\" 2>&1 | head -1)
  echo \"P=\$p => \$out\"
  sleep 1
done
rm -f \"\$tmp/probe.slurm\" 2>/dev/null; rmdir \"\$tmp\" 2>/dev/null
"
fi

echo
echo "======================================================================"
echo " 填写提示：把上面结果抄进 setting/<hpc>.yaml（字段说明见 docs/DEPLOY-NEW-CLUSTER.md 第 8 节）"
echo " 必填：name / ssh_host / work_dir / partition / max_cpus / vasp.* / potcar_dir / conda_env"
echo " 建议：先跑 autozt check --cluster <hpc> 当验收，再提交冒烟作业"
echo "======================================================================"
