#!/bin/bash
# phono3py BTE 提交模板（step6_kappa, solver=phono3py）。
# 占位符：{{JOBNAME}} {{CPUS_PER_TASK}} {{QOS}} {{CONDA_SH}} {{CONDA_ENV}} {{P3PY_CMD}}
#
# 并行布局（2026-09-18 实测后修正）——本环境的 phono3py 是 **OpenMP-only**：
#   C 扩展只链 libgomp，包内没有任何 mpi4py/MPI 引用（实测 phono3py 3.24.0：
#   mpi4py 未安装、_phono3py*.so 的 ldd 里无 MPI 库，包内 0 处 MPI 代码）。
#   因此：
#     * **绝对不要用 mpirun**：那会启动 N 个各自算全量的进程，互相覆盖输出文件；
#     * 单节点唯一有效的并行 = 「1 个进程 × N 个 OpenMP 线程」。
#   所以这里申请 --ntasks=1 --cpus-per-task=N（N = step.conf 的 P3PY_OMP_THREADS，
#   默认 96 = jzzn 192 逻辑核节点的物理核数），并把 OMP_NUM_THREADS 设成同一个数。
#   旧模板写 --ntasks-per-node=48 且不写 --cpus-per-task：Slurm 分 48 个 task、每 task
#   只有 1 核，而实际只跑一个进程 → 48 个 OMP 线程全挤在 1 个核上，整节点被白占。
#   以后若换成 mpi4py/MPI 版 phono3py，再改成 mpirun -n <ranks> + OMP=1。
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={{CPUS_PER_TASK}}
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos={{QOS}}
cd $SLURM_SUBMIT_DIR
source {{CONDA_SH}}
conda activate {{CONDA_ENV}}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-96}
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_NESTED=FALSE
export OMP_STACKSIZE=1G GOMP_STACKSIZE=1G
# 单进程 + 大 OMP：BLAS 再开线程会与 OMP 抢核（嵌套超订），统一压成 1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMPY_NUM_THREADS=1
echo "[phono3py] 单进程 OMP 线程=$OMP_NUM_THREADS（OpenMP-only：本模板不能用 mpirun）env=$CONDA_DEFAULT_ENV"
phono3py --version || true
{{P3PY_CMD}}
