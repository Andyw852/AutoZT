#!/bin/bash
# ShengBTE BTE 提交模板（step6_kappa, SOLVER=shengbte）。占位符 {{JOBNAME}} {{SHENGBTE_EXE}}
#   {{KAPPA_2D_FACTOR}} {{KAPPA_2D_META}} {{KAPPA_2D_THICK2D}}（2D 厚度归一化，3D 时 FACTOR=1.0）
#   {{NTASKS}} {{CPUS_PER_TASK}} {{QOS}} —— 由 gen_step6_kappa.py 按集群 NUMA 拓扑与
#   step.conf 的 SBATCH_QOS 算好填入（2026-09-16 从 taskflow-v2.0 同步；见下"并行布局"）
# 输入 CONTROL + FORCE_CONSTANTS_2ND/3RD 已由 gen_step6 备好（力常数拷自 S5_fc/shengbte/）。
#
# =====================================================================
# ★★ 并行布局：这一步决定「算不算得出来」，不是性能调优 ★★
#
# 坑 1（最致命）：ShengBTE 靠 **MPI 按 q 点并行**。
#     写 mpirun -n 1 = 单进程 = **串行**。日志会自报
#     "running on 1 MPI process(es)"，然后长时间停在
#     "about to obtain the spectrum" 一个产物都不出（实测 11.8 h 零产物）。
#     正确做法 = **一个 MPI rank 占一个 NUMA 节点**，且 **rank 数必须 ≤ NUMA 节点数**：
#         ntasks-per-node = 总核数 / 每个 NUMA 节点的核数
#         cpus-per-task   = 每个 NUMA 节点的核数
#     jzzn 计算节点实测：192 核 = 2 socket x 96 = **8 个 NUMA 节点 x 24 核**。
#         192 核 -> 8 rank x 24 线程；96 核 -> 4 rank x 24 线程。
#     ★ rank 数 > NUMA 节点数时 OpenMPI 会**循环复用**节点并把多个 rank 绑到同一组核上：
#       实测 24 rank x 8 线程时每个 NUMA 节点塞进 3 个 rank、24 个线程挤在 8 个核上，
#       每线程只拿到 33% 的核，整作业只用 64/192 核（用 srun --overlap 进作业读
#       /proc/<pid>/status 的 Cpus_allowed_list 可复现：并集只有 64 个 CPU）。
#       所以 CORES_PER_NUMA 必须填**真实每 NUMA 核数 24**（不是 8）。
#
# 坑 2：OpenMPI 4.x 下 --map-by numa:PE=N 与任何 --bind-to ... 是**冲突写法**：
#       · 配 --bind-to core -> "would result in binding more processes than cpus"
#       · 配 --bind-to numa -> "a conflicting binding policy was specified"
#       · 干脆不给        -> 报 "no binding need be specified"（它自带的默认就是对的）
#     而且 PE 不能超过 NUMA 节点内的核数（jzzn = 24），写大了会报
#         "a directive was also given to map to an object level that has less cpus"
#     并 8 秒退出。所以这里用**裸 --map-by numa**：不带 PE、不带 --bind-to，
#     让 OpenMPI 按 Slurm 给每个 task 的核集把 rank 铺到各 NUMA 节点上。
#
# 换集群时只需改 step.conf 的 SHENGBTE_TOTAL_CORES / SHENGBTE_CORES_PER_NUMA，
# 不用动本模板。
# =====================================================================
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={{NTASKS}}
#SBATCH --cpus-per-task={{CPUS_PER_TASK}}
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos={{QOS}}
module purge
module load gcc/14.1
module load openmpi/4.0.1
source /public/home/<user>/miniconda3/etc/profile.d/conda.sh
conda activate atomate2_p_a
cd $SLURM_SUBMIT_DIR

for f in CONTROL FORCE_CONSTANTS_2ND FORCE_CONSTANTS_3RD; do
    [ ! -f "$f" ] && echo "缺 $f" >&2 && exit 1
done

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_NESTED=FALSE
export BLIS_NUM_THREADS=$OMP_NUM_THREADS
export AOCL_ENABLE_INSTRUCTIONS=AVX512
export LD_LIBRARY_PATH=/public/home/<user>/software/aocl-gcc/5.0.0/gcc/lib_LP64:$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
unset MKL_NUM_THREADS MKL_DEBUG_CPU_TYPE I_MPI_PMI_LIBRARY
ulimit -s unlimited
# [FIX P49/OMP 2026-09-18，移植自 taskflow-v2.0] "第一个温度（100 K）的 RTA 写完、
# rank 0 进入 cumulative-κ 段就 SIGSEGV（rc=139）"的根因 = **OpenMP 工作线程栈溢出**：
# conductivity.f90 的 CumulativeTConduct 用 reduction(+:results)，而
# results(3,3,nbands,nticks) 在 nbands=384、nticks=100（默认）时约 2.8 MB **每线程一份**；
# ulimit -s unlimited 只管主线程，libgomp 的线程栈由 OMP_STACKSIZE 决定。
# 实测（jzzn，1x1x1 网格 + Mg4C60 128 原子胞，每轮数十秒）：
#   OMP=1 未设 STACKSIZE → normal exit；OMP=2/4 未设 → SIGSEGV（固定在 100 K 之后）；
#   OMP=2/4/8 且 OMP_STACKSIZE=8M/64M/256M → 全部 normal exit；
#   与 MPI rank 数（1/2/4/12）、q 网格（1³/2³/7³）无关；自带的 InAs 例子（nbands=6）
#   私有副本仅约 43 KB，OMP=2 下正常，所以此前一直没暴露。
export OMP_STACKSIZE=1G
export GOMP_STACKSIZE=1G

# ---- 退化保护：rank 数 = 1 就是串行，宁可当场失败也不要白跑一夜 ----
#   来源（taskflow-v2.0，2026-09-16 同步）：ShengBTE 靠 MPI 按 q 点并行，写 mpirun -n 1
#   = 单进程 = 串行，日志自报 "running on 1 MPI process(es)" 后长时间停在
#   "about to obtain the spectrum"，实测 11.8 h 零产物。
if [ "$SLURM_NTASKS" -le 1 ]; then
    echo "SLURM_NTASKS=$SLURM_NTASKS —— ShengBTE 靠 MPI 按 q 点并行，单进程等于串行" >&2
    echo "（实测 11.8 h 停在 about to obtain the spectrum 无产物）。请检查 step.conf 的" >&2
    echo "SHENGBTE_TOTAL_CORES / SHENGBTE_CORES_PER_NUMA。" >&2
    exit 2
fi
echo "[shengbte] MPI ranks=$SLURM_NTASKS  每 rank OMP 线程=$SLURM_CPUS_PER_TASK  合计 $((SLURM_NTASKS*SLURM_CPUS_PER_TASK)) 核"
echo "[shengbte] 节点=$SLURM_JOB_NODELIST  开始 $(date)"

# --map-by numa：一个 rank 一个 NUMA 域。不要再加 --bind-to 或 PE=N
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

# 汇总：优先 _CONV（迭代解），退回 _RTA
python - <<'PY'
import glob, json
def _read_kt(fname):
    rows = [l.split() for l in open(fname) if l.strip() and not l.startswith("#")]
    out = []
    for r in rows:
        try:
            vals = [float(x) for x in r[:10]]
        except ValueError:
            continue
        if abs(vals[1]) < 1e6 and abs(vals[5]) < 1e6 and abs(vals[9]) < 1e6:
            out.append(vals)
    return out

cand = ["BTE.KappaTensorVsT_CONV", "BTE.KappaTensorVsT_RTA"]
best = None
for c in cand:
    rows = _read_kt(c) if glob.glob(c) else []
    if not rows:
        continue
    if len(rows) >= 3:  # CONV 发散（1e147 量级）被过滤后可能 <3 行 -> 回退 RTA
        best = (c, rows)
        if c.endswith("_CONV"):
            # CONV 正常温度点数（CONTROL T 扫描计数）不足时：发散点被过滤 =
            # 曲线不完整，回退 RTA（更稳定、全温度）
            # T_min 原来写死成 100.0：step.conf 把 T_MIN 改成别的值时温度点数预期就
            # 算错，正常的 CONV 结果会被误判成"发散被过滤"而退回 RTA。三个都从 CONTROL 读。
            _tmin = _tmax = _ts = None
            try:
                for ln in open("CONTROL"):
                    s = ln.strip()
                    for _key in ("T_min", "T_max", "T_step"):
                        if s.startswith(_key + "="):
                            _val = float(s.split("=")[1].split(",")[0].strip())
                            if _key == "T_min":
                                _tmin = _val
                            elif _key == "T_max":
                                _tmax = _val
                            else:
                                _ts = _val
            except Exception:
                pass
            _nt_expected = 8
            if None not in (_tmin, _tmax, _ts) and _ts:
                _nt_expected = int(round((_tmax - _tmin) / _ts)) + 1
            else:
                print("[shengbte] CONTROL 的 T_min/T_max/T_step 没读全，温度点数预期退回 %d"
                      % _nt_expected)
            if len(rows) < _nt_expected:
                print("[shengbte] CONV 温度点 %d < 预期 %d（发散被过滤），回退 RTA" % (len(rows), _nt_expected))
                continue
        break
f, rows = best if best else (None, [])
d = {"KAPPA_DONE": bool(f)}
# 2D 厚度归一化：ShengBTE 的分母和 phono3py 一样是【含真空的胞体积】(h_perp=V/A)，
#   所以面内 κ 被胞高稀释，必须乘 FACTOR = h_perp/d（d=层有效厚度）才是层本身的量。
#   FACTOR 由 gen_step6_kappa.py 用 _common/thickness_2d.py 算好传进来；3D 时 =1.0。
FACTOR = {{KAPPA_2D_FACTOR}}
META = json.loads(r'''{{KAPPA_2D_META}}''')
THICK2D = json.loads(r'''{{KAPPA_2D_THICK2D}}''')
d.update(META)
if THICK2D:
    json.dump(THICK2D, open("thickness_2d.json", "w"), ensure_ascii=False, indent=2)
    d["thickness_2d_source"] = "kl-dft-cpu:step6_kappa"
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
if f:
    d["source"] = f
    d["temperatures"] = [r[0] for r in rows]
    # ShengBTE KappaTensorVsT：col0=T，col1..9=kappa 张量 xx xy xz yx yy yz zx zy zz
    d["kappa_xx_yy_zz"] = [[r[1], r[5], r[9]] for r in rows]
    # 2D 的 zz 分量无物理意义，面内单独出一份方便下游直接用
    d["kappa_inplane_xx_yy"] = [[r[1], r[5]] for r in rows]
    if abs(FACTOR - 1.0) > 1e-9:
        d["kappa_2d_normalized_xx_yy_zz"] = [[float(v) * FACTOR for v in row]
                                            for row in d["kappa_xx_yy_zz"]]
        d["kappa_2d_normalized_inplane_xx_yy"] = [[r[1] * FACTOR, r[5] * FACTOR] for r in rows]
if not d["KAPPA_DONE"]:
    # [FIX P49 2026-09-18，移植自 taskflow-v2.0] 只有 small-grain limit 时不要静默失败：
    # 记下"三声子段未完成"，让 KAPPA_DONE=false 变成**有信息**的失败，而不是看起来像
    # "什么都没算"的空白。BTE.KappaTensorVsT_sg 是纯谐性、完全不含三声子散射的最小 κ
    # （说明见 gen_step6_kappa.py 的 P49 注释），单独存在时**不可**当作热导率交付或引用。
    try:
        _sg = _read_kt("BTE.KappaTensorVsT_sg")
    except Exception:
        _sg = []
    if _sg:
        d["small_grain_limit_present"] = True
        d["small_grain_limit_rows"] = len(_sg)
        d["warning"] = ("只找到 BTE.KappaTensorVsT_sg = small-grain limit（纯谐性、"
                        "无三声子散射），不是热导率；RTA(_RTA)/迭代(_CONV) 段未完成，"
                        "不可当作 κ 引用或对外交付")
json.dump(d, open("kappa_summary.json", "w"), ensure_ascii=False, indent=2)
print("KAPPA_DONE" if d["KAPPA_DONE"] else "NO_KAPPA")
PY
