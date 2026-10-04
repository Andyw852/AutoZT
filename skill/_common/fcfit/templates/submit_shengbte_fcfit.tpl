#!/bin/bash
# fit-fc-thermal S2_kappa 的 ShengBTE 提交模板（SOLVER=shengbte，共享引擎 _common/fcfit）。
# 占位符：JOBNAME / CONDA_SH / CONDA_ENV / SHENGBTE_EXE / NTASKS / CPUS_PER_TASK / QOS
#
# 流程（与 phono3py 的 submit_kappa.tpl 不同：ShengBTE 必须由 mpirun 直接启动才能吃到
#   MPI 按 q 点并行，不能套在 python 子进程里）：
#   ① kappa_driver.py --prepare-shengbte  fc2/fc3.hdf5 -> FORCE_CONSTANTS_2ND/3RD + CONTROL
#   ② mpirun（全套 MPI/OMP/BLIS/AOCL/UCX 加速参数）跑 ShengBTE
#   ③ kappa_driver.py --collect-shengbte  解析 BTE.KappaTensorVsT_RTA -> kappa_summary.json
#
# 本文件是 jzzn（AOCL 版）的默认模板：module gcc/14.1 + openmpi/4.0.1 + aocl-gcc BLAS。
# 换集群时在 setting/<hpc>/templates/ 放同名文件覆盖（模块名 / AOCL 路径 / UCX 网卡各不同）。
#
# =====================================================================
# 并行布局：这一步决定「算不算得出来」，不是性能调优
#
# 坑 1（最致命）：ShengBTE 靠 MPI 按 q 点并行。mpirun -n 1 = 单进程 = 串行，日志自报
#   "running on 1 MPI process(es)" 后长时间停在 "about to obtain the spectrum" 零产物
#   （实测 11.8 h）。正确布局 = 一个 MPI rank 占一个 NUMA 节点，且 rank 数必须 <= NUMA 节点数：
#       ntasks-per-node = 总核数 / 每 NUMA 节点核数
#       cpus-per-task   = 每 NUMA 节点核数
#   jzzn 计算节点实测：192 核 = 2 socket x 96 = 8 个 NUMA 节点 x 24 核。
#       192 核 -> 8 rank x 24 线程；96 核 -> 4 rank x 24 线程。
#   NTASKS / CPUS_PER_TASK 由 gen_step2_kappa.py 按 SHENGBTE_TOTAL_CORES /
#   SHENGBTE_CORES_PER_NUMA 算好填入 —— 写死 ntasks=1 会退化成串行。
#
# 坑 2：OpenMPI 4.x 下 --map-by numa:PE=N 与 --bind-to ... 是冲突写法（--bind-to core
#   报 "binding more processes than cpus"，--bind-to numa 报 conflicting binding policy），
#   PE 也不能超过 NUMA 核数（jzzn=24）。所以用裸 --map-by numa：让 OpenMPI 按 Slurm 给
#   每个 task 的核集把 rank 铺到各 NUMA 节点上，不再额外 --bind-to 或 PE=N。
# =====================================================================
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={{NTASKS}}
#SBATCH --cpus-per-task={{CPUS_PER_TASK}}
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos={{QOS}}

cd $SLURM_SUBMIT_DIR

module purge
module load gcc/14.1
module load openmpi/4.0.1

# conda：CONDA_ENV 是 venv 目录就直接激活；否则 source CONDA_SH 再 conda activate。
# 两者都留空（如 setting/local.yaml）= 用 PATH 里现成的 python，不激活任何环境。
if [ -n "{{CONDA_ENV}}" ] && [ -x "{{CONDA_ENV}}/bin/python" ]; then
    source "{{CONDA_ENV}}/bin/activate"
elif [ -n "{{CONDA_SH}}" ]; then
    if [ -f "{{CONDA_SH}}" ]; then
        source "{{CONDA_SH}}"
        if [ -n "{{CONDA_ENV}}" ]; then
            conda activate "{{CONDA_ENV}}"
        fi
    else
        echo "[WARN] conda.sh not found: {{CONDA_SH}} -- using python from PATH" >&2
    fi
fi

# ① 准备：fc2/fc3.hdf5 -> FORCE_CONSTANTS_2ND/3RD + CONTROL（需要 phono3py/hiphive/ase）
python -c "import phono3py, hiphive; print('phono3py', phono3py.__version__)" || exit 1
python kappa_driver.py kappa_config.json --prepare-shengbte || exit 1

for f in CONTROL FORCE_CONSTANTS_2ND FORCE_CONSTANTS_3RD; do
    [ ! -f "$f" ] && echo "MISSING $f" >&2 && exit 1
done

# ② ShengBTE 加速环境（jzzn AOCL 版；换集群用 setting/<hpc>/templates/ 覆盖）
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_NESTED=FALSE
export BLIS_NUM_THREADS=$OMP_NUM_THREADS
export AOCL_ENABLE_INSTRUCTIONS=AVX512
export LD_LIBRARY_PATH=/public/home/wangchao/software/aocl-gcc/5.0.0/gcc/lib_LP64:$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
unset MKL_NUM_THREADS MKL_DEBUG_CPU_TYPE I_MPI_PMI_LIBRARY
ulimit -s unlimited
# [FIX P49/OMP] 第一个温度 RTA 写完、rank 0 进 cumulative-κ 段就 SIGSEGV（rc=139）的根因 =
#   OpenMP 工作线程栈溢出：libgomp 线程栈由 OMP_STACKSIZE 决定（ulimit -s 只管主线程）。
#   results(3,3,nbands,nticks) 在 nbands 大时每线程一份可达数 MB，设 1G 兜底。
export OMP_STACKSIZE=1G
export GOMP_STACKSIZE=1G

# 退化保护：rank=1 就是串行，宁可当场失败也不要白跑一夜
if [ "$SLURM_NTASKS" -le 1 ]; then
    echo "SLURM_NTASKS=$SLURM_NTASKS —— ShengBTE 靠 MPI 按 q 点并行，单进程等于串行" >&2
    echo "（实测 11.8 h 停在 about to obtain the spectrum 无产物）。请检查 step.conf" >&2
    echo "的 SHENGBTE_TOTAL_CORES / SHENGBTE_CORES_PER_NUMA。" >&2
    exit 2
fi
echo "[shengbte] exe={{SHENGBTE_EXE}}"
echo "[shengbte] MPI ranks=$SLURM_NTASKS  每 rank OMP 线程=$SLURM_CPUS_PER_TASK  合计 $((SLURM_NTASKS*SLURM_CPUS_PER_TASK)) 核  开始 $(date)"

# --map-by numa：一个 rank 一个 NUMA 域。不要再加 --bind-to 或 PE=N（见坑 2）
mpirun -n $SLURM_NTASKS \
    --map-by numa \
    --report-bindings \
    --mca pml ucx --mca osc ucx --mca btl ^openib,tcp \
    -x UCX_TLS=rc,sm,self -x UCX_NET_DEVICES=mlx5_0:1 -x UCX_LOG_LEVEL=error \
    -x OMP_NUM_THREADS -x OMP_PROC_BIND -x OMP_PLACES \
    -x OMP_STACKSIZE -x GOMP_STACKSIZE \
    -x BLIS_NUM_THREADS -x AOCL_ENABLE_INSTRUCTIONS -x LD_LIBRARY_PATH \
    {{SHENGBTE_EXE}} > shengbte.log 2>&1
_rc=$?
echo "[shengbte] 结束 $(date)  rc=$_rc"
[ "$_rc" -eq 0 ] || { echo "ShengBTE 失败（rc=$_rc），看 shengbte.log" >&2; exit $_rc; }

# ③ 汇总：解析 RTA κ -> kappa_summary.json
python kappa_driver.py kappa_config.json --collect-shengbte || exit 1
