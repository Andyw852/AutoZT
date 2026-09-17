#!/bin/bash
# fourphonon(ShengBTE 同源, 3ph RTA) 提交模板 —— step6_kappa SOLVER=fourphonon
# 占位符: {{JOBNAME}} {{FOURPHONON_EXE}} {{FOURPHONON_NGPU}} {{FOURPHONON_CPUS_PER_GPU}}
#         {{KAPPA_2D_FACTOR}} {{KAPPA_2D_META}} {{KAPPA_2D_THICK2D}}（2D 归一化）
# 输入 CONTROL + FORCE_CONSTANTS_2ND/3RD 由 gen_step6 备好(力常数拷自 S5_fc/shengbte/)。
#
# fourphonon GPU 版要点(2026-09 实战验证, 见 kl-dft-cpu/README.md):
#   1. 必须 multi-GPU 修复版(acc_set_device_num 自动选卡 + MPI_IN_PLACE 归约);
#      旧单 rank 版枚举段极慢(GPU 永不点火)。
#   2. rank 数 = GPU 数, 不要超过机器卡数。
#   3. Intel 2019 mpiexec.hydra 部分机器损坏(连 hostname 都 FPE) -> 用 NVHPC OpenMPI。
#   4. AOCC/flang 编译版大胞枚举会卡死(12h 零进展), 勿用于 >20 原子体系。
#SBATCH --partition=gpu
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={{FOURPHONON_NGPU}}
#SBATCH --cpus-per-task={{FOURPHONON_CPUS_PER_GPU}}
#SBATCH --gres=gpu:{{FOURPHONON_NGPU}}
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=normal
#SBATCH --time=48:00:00
cd $SLURM_SUBMIT_DIR

for f in CONTROL FORCE_CONSTANTS_2ND FORCE_CONSTANTS_3RD; do
    [ ! -f "$f" ] && echo "missing $f" >&2 && exit 1
done

# ---- NVHPC OpenMPI + CUDA + spglib/openblas(按机器改路径) ----
NVHPC=\${NVHPC_ROOT:-/opt/nvhpc/24.11/Linux_x86_64/24.11}
OMPI=$NVHPC/comm_libs/12.6/openmpi4/openmpi-4.1.5
CUDA=\${CUDA_ROOT:-/usr/local/cuda}
SPG=\${SPGLIB_DIR:-/path/to/spglib/lib}
BLAS=\${BLAS_DIR:-/path/to/openblas/lib}
export PATH=$OMPI/bin:$CUDA/bin:$PATH
export LD_LIBRARY_PATH=$SPG:$BLAS:$NVHPC/compilers/lib:$OMPI/lib:$NVHPC/math_libs/12.6/lib64:$CUDA/lib64:$LD_LIBRARY_PATH

export OMP_NUM_THREADS={{FOURPHONON_CPUS_PER_GPU}}
export OMP_PROC_BIND=spread OMP_PLACES=threads
export OMP_STACKSIZE=1G
ulimit -s unlimited

echo "fourphonon {{FOURPHONON_NGPU}}-GPU start: $(date)"

# rank 数 = GPU 数; 程序内 acc_set_device_num(myid) 自动选卡
# 每 rank 占 CPUS_PER_GPU 个核(map-by numa:PE)供枚举段 OMP 使用。
#
# ★ 不要再加 --bind-to —— OpenMPI 4.x 下它与 map-by PE=N 是**冲突写法**：
#     配 --bind-to core -> "would result in binding more processes than cpus"
#     配 --bind-to numa -> "a conflicting binding policy was specified"
#     不给             -> 它自己报 "no binding need be specified"（默认就是对的）
#   而且 PE 不能超过该节点 NUMA 域里的核数（jzzn 是 8 核/域，写 PE=16 会 8 秒退出）。
mpirun -n {{FOURPHONON_NGPU}} --map-by numa:PE={{FOURPHONON_CPUS_PER_GPU}} \
    {{FOURPHONON_EXE}} > fourphonon.log 2>&1
echo "EXIT=$? end: $(date)"

# 汇总 RTA 结果 -> kappa_summary.json
python - <<'PY'
import glob, json
def _read_kt(fname):
    rows = []
    for l in open(fname):
        l = l.strip()
        if not l or l.startswith('#'):
            continue
        parts = l.split()
        if len(parts) < 10:
            continue
        try:
            vals = [float(x) for x in parts[:10]]
        except ValueError:
            continue
        if abs(vals[1]) < 1e6 and abs(vals[5]) < 1e6 and abs(vals[9]) < 1e6:
            rows.append(vals)
    return rows

rows = _read_kt('BTE.KappaTensorVsT_RTA') if glob.glob('BTE.KappaTensorVsT_RTA') else []
d = {'KAPPA_DONE': bool(rows), 'source': 'BTE.KappaTensorVsT_RTA' if rows else None}
# 2D 厚度归一化：fourphonon 与 ShengBTE 同源，分母同样是【含真空的胞体积】
#   (h_perp=V/A)，面内 κ 被胞高稀释，乘 FACTOR = h_perp/d 才是层本身的量。
FACTOR = {{KAPPA_2D_FACTOR}}
META = json.loads(r'''{{KAPPA_2D_META}}''')
THICK2D = json.loads(r'''{{KAPPA_2D_THICK2D}}''')
d.update(META)
if THICK2D:
    json.dump(THICK2D, open('thickness_2d.json', 'w'), ensure_ascii=False, indent=2)
    d['thickness_2d_source'] = 'kl-dft-cpu:step6_kappa'
    # 材料级 thickness_2d.json（文档 P0-2(b) 的 kl/ke 共用契约）：谁先跑谁写，后跑方读它做对照。
    #   作业 cwd 是 <材料>/<步骤>，材料根 = ".."；先确认 ".." 里真有别的步骤目录才写，
    #   避免路径猜错写到别处。best-effort，失败只打一行，不拦作业。
    try:
        import glob as _g, os as _os
        if _g.glob(_os.path.join("..", "step*")):
            from thickness_2d import write_material_level
            print(write_material_level("..", THICK2D))
        else:
            print("[..] 父目录里没识别出步骤目录，跳过材料级 thickness_2d.json")
    except Exception as _e:
        print("[WARN] 材料级 thickness_2d.json 写入跳过: %s" % _e)
if rows:
    d['temperatures'] = [r[0] for r in rows]
    d['kappa_xx_yy_zz'] = [[r[1], r[5], r[9]] for r in rows]
    # 2D 的 zz 分量无物理意义，面内单独出一份
    d['kappa_inplane_xx_yy'] = [[r[1], r[5]] for r in rows]
    if abs(FACTOR - 1.0) > 1e-9:
        d['kappa_2d_normalized_xx_yy_zz'] = [[float(v) * FACTOR for v in row]
                                            for row in d['kappa_xx_yy_zz']]
        d['kappa_2d_normalized_inplane_xx_yy'] = [[r[1] * FACTOR, r[5] * FACTOR] for r in rows]
json.dump(d, open('kappa_summary.json', 'w'), ensure_ascii=False, indent=2)
print('KAPPA_DONE' if d['KAPPA_DONE'] else 'NO_KAPPA')
PY
