#!/bin/bash
# S5_fc 拟合作业模板（FIT_ENGINE=pheasy）。移植自用户 pheasy 拟合脚本，去掉 VCA/内存监控/
# 尾部自动提交 κ；输入准备(prep)与收尾(collect/post)交给 kl_fc_backends.py。
# 占位符：{{JOBNAME}} {{DIM}} {{FIT_METHOD}} {{ENABLE_FC}} {{C3_CUTOFF}} {{CUT3_CANDIDATES}} {{CUT3_BOOTSTRAP}} {{NULL_SPACE_EPS}} {{RASR}}
# 资源（cpus_per_task/qos）建议按体系用 step.conf 的 [submit] 段覆盖；pheasy 较吃核与内存。
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=premium
cd $SLURM_SUBMIT_DIR
source /public/home/<user>/miniconda3/etc/profile.d/conda.sh
conda activate atomate2_p_a
set -e

# ===== 环境自检 =====
for _m in numpy scipy phonopy spglib phono3py pheasy; do
    python -c "import ${_m}" 2>/dev/null && echo "  ✅ ${_m}" \
        || { echo "  ❌ ${_m} 导入失败，检查 atomate2_p_a 环境"; exit 1; }
done

# ===== 用户参数（gen 从 step.conf 注入）=====
DIM="{{DIM}}"                 # 超胞对角三整数
FIT_METHOD="{{FIT_METHOD}}"   # LASSO | RFE | OLS
ENABLE_FC={{ENABLE_FC}}       # 2|3|4
C3_CUTOFF="{{C3_CUTOFF}}"     # fc3 截断 Å，None=不截断
CUT3_CANDIDATES="{{CUT3_CANDIDATES}}"   # S4 记的三阶截断候选（空格分隔；空=单截断老路）
CUT3_BOOTSTRAP={{CUT3_BOOTSTRAP}}       # 每档帧 bootstrap 次数（只重跑 -d/-f；-c 复用）
NULL_SPACE_EPS={{NULL_SPACE_EPS}}
# RASR = 旋转不变性(Born-Huang) + 平衡条件(Huang, 零应力)。gen 按 DIM 注入：
#   2D→BHH（必加，否则 ZA 近 Γ 线性化/出虚频），3D→none（文献结论：对体材料可忽略）。
#   pheasy 只在【零空间构造步 -c】读 RASR；-f 读的是 ns_*.npz，把 --rasr 放 -f 上无效。
#   pheasy 的 --rasr 只接受 BH/H/BHH，none = 不传这个开关（不是传 "none"）。
RASR="{{RASR}}"
RASR_FLAGS=""
if [ -n "${RASR}" ] && [ "${RASR}" != "none" ]; then
    [[ "${RASR}" =~ ^(BHH|BH|H)$ ]] || { echo "❌ RASR 非法: ${RASR}（BHH|BH|H|none）"; exit 1; }
    RASR_FLAGS="--rasr ${RASR}"
fi
FIT_ORDER=${ENABLE_FC}

# ===== 并行 =====
NCPU=${SLURM_CPUS_PER_TASK:-48}
NCPU_BLAS=$(( NCPU>32 ? 32 : NCPU ))
NCPU_DISP=${NCPU}
if [ "${FIT_METHOD}" = "RFE" ] || [ "${FIT_METHOD}" = "OLS" ]; then
    NCPU_LOKY=1; NCPU_FIT_BLAS=${NCPU}
else
    NCPU_LOKY=1; NCPU_FIT_BLAS=${NCPU}
fi
[[ "${FIT_METHOD}" =~ ^(LASSO|RFE|OLS)$ ]] || { echo "❌ FIT_METHOD 非法: ${FIT_METHOD}"; exit 1; }
[[ "${ENABLE_FC}" =~ ^[234]$ ]] || { echo "❌ ENABLE_FC 非法: ${ENABLE_FC}"; exit 1; }

# ===== 输入准备：vasprun → FORCES_FC3 → POSCAR/SPOSCAR/dataset_*.npy =====
python kl_fc_backends.py prep fit_config.json
for f in POSCAR SPOSCAR dataset_disps.npy dataset_forces.npy; do
    [ ! -f "$f" ] && echo "❌ prep 后仍缺 $f" && exit 1
done

# ===== 加速 / BLAS 线程 =====
export PHEASY_SVD_THRESHOLD=500
export PHEASY_ASR_COMBINED=1
export PHEASY_ASR_LWORK_LIMIT=1500000000
export PHEASY_SM_DTYPE=float32
export PHEASY_SM_THR=1e-12
export PHEASY_ASR_SPARSE=1
export PHEASY_ASR_SPARSE_THR=1e-10
export PHEASY_ASR_COL_BLOCK=5000
# 零空间秩判据：施加 RASR(BHH/BH/H) 时约束直接写进零空间，秩容差决定"哪些方向被当成约束"
#   —— 放宽会漏约束、收紧会吃掉真实自由度。1e-6 与 MACE 链
#   (_common/mlff/gen_step3_fc.py) 取同一个值，两条链的 2D 行为才对得上。
export PHEASY_NS_RANK_TOL=1e-6
export OPENBLAS_NUM_THREADS=${NCPU_BLAS}
export OMP_NUM_THREADS=${NCPU_BLAS}
export MKL_NUM_THREADS=${NCPU_BLAS}
export PHEASY_BLAS_THREADS=${NCPU_BLAS}
export PHEASY_USE_CELER=1

# RFE 触发（复用用户脚本的 PHEASY_USE_RFE 重定向）
if [ "${FIT_METHOD}" = "RFE" ]; then
    export PHEASY_USE_RFE=1 PHEASY_RFE_TWOLEVEL=1 MKL_INTERFACE_LAYER=ILP64 PHEASY_RFE_MKL=1
    export PHEASY_RFE_STEP=0.1 PHEASY_RFE_RIDGE_ALPHA=1e-11 PHEASY_RFE_CV=5
    export PHEASY_RFE_LSMR_MAXITER=60000 PHEASY_RFE_WARM_START=1
    export PHEASY_COLNORM_FRAMES=24 PHEASY_COLNORM_EXACT=0 PHEASY_RFE_ONE_SE=1
    echo "RFE 启用：step=0.1 ridge=1e-11 cv=5"
elif [ "${FIT_METHOD}" = "OLS" ]; then
    export MKL_INTERFACE_LAYER=ILP64 PHEASY_OLS_TWOLEVEL=1 PHEASY_RFE_MKL=1
    # PHEASY_OLS_RIDGE 必须保持 0。pheasy 的 OLS 把它当 Tikhonov 阻尼，并且按
    # damp = sqrt(ridge*ndata) 放大（core/optimizer.py:_ols_lsmr）：1e-4 在
    # ndata=7500 时 damp≈0.87，直接把系数按尺度压偏，而 pheasy 仍然退 0。
    # 同一数据集实测（BaS 250 原子 / 10 帧 / fc2+fc3 / OLS）：
    #   ridge=0    相对误差 0.022  correlation 0.9997   (LSMR 509 次迭代)
    #   ridge=1e-4 相对误差 0.580  correlation 0.8991   (LSMR  18 次迭代)
    # 另外 ridge>0 会让 GPU 常驻 OLS 直接回退 CPU（optimizer.py:2643）。
    # 对比：RFE 的 PHEASY_RFE_RIDGE_ALPHA=1e-11 走 sklearn Ridge 的归一化路径，
    # 不放大，可以保留。
    export PHEASY_OLS_MAXITER=500 PHEASY_OLS_RIDGE=0 PHEASY_OLS_ATOL=1e-6 PHEASY_OLS_BTOL=1e-6
    echo "OLS 启用：两级 matvec, MKL ILP64（内存较高）"
else
    python -c "from celer import Lasso" 2>/dev/null || { echo "❌ LASSO 需 celer，禁止提交"; exit 1; }
fi

# ===== 数据处理：dataset_*.npy → disp_matrix.pkl/force_matrix.pkl + ndata =====
python << 'PYEOF'
import numpy as np, pickle
from phonopy.interface.vasp import read_vasp
sup = read_vasp('SPOSCAR'); natom = len(sup.numbers)
d = np.load('dataset_disps.npy'); f = np.load('dataset_forces.npy')
eq = d[-1]
# prep 写的是笛卡尔位移 + 末帧零平衡帧；这里统一扣末帧还原（幂等）
is_frac = (d.min() >= -0.15 and d.max() <= 1.15 and 0.2 < d.mean() < 0.8)
if is_frac:
    dd = d - eq; dd = np.where(dd > 0.5, dd-1, np.where(dd < -0.5, dd+1, dd))
    dcart = np.einsum('ij,njk->nik', sup.cell.T, np.transpose(dd,(0,2,1)))
    dcart = np.transpose(dcart,(0,2,1))
else:
    dcart = d
disps = (dcart - dcart[-1])[:-1]
forces = (f - f[-1])[:-1]
rms = np.sqrt((disps**2).mean())
assert 1e-6 < rms < 1.0, "RMS位移异常 %.3e" % rms
print("参与拟合 %d 帧, RMS位移=%.4f Å, natom_super=%d" % (len(disps), rms, natom))
# 照 doc2 落全 4 个文件：pheasy --disp_file 依赖 pkl；npy 供其它读取路径/对拍，一并写全以防口径不一致
np.save('dataset_disps_cartesian.npy', disps)
np.save('dataset_forces_corrected.npy', forces)
pickle.dump(disps,  open('disp_matrix.pkl','wb'))
pickle.dump(forces, open('force_matrix.pkl','wb'))
open('ndata_total.txt','w').write(str(len(disps)))
open('natom_super.txt','w').write(str(natom))
PYEOF
NDATA=$(cat ndata_total.txt)
export PHEASY_CV_GROUP_SIZE=$(( 3 * $(cat natom_super.txt) ))
echo "ndata=${NDATA}  CV_GROUP_SIZE=${PHEASY_CV_GROUP_SIZE}"

# ===== 三阶截断扫描（2026-09-24 user 定）：候选≥2 时逐档拟合 + 逐壳层稳定性 =====
#   nominal 截断（PHEASY_C3_CUTOFF）的 fc2/fc3 + pheasy_{c,f}.log 由 scan_pheasy
#   留在 cwd，下面 collect/metrics/post 原样复用。
if [ -n "${CUT3_CANDIDATES}" ]; then
    echo "【截断扫描】候选=${CUT3_CANDIDATES} bootstrap=${CUT3_BOOTSTRAP}"
    python kl_fc_backends.py scan_pheasy fit_config.json
else
# ===== 参数拼装 =====
C_FLAG=""
[ "${C3_CUTOFF}" != "None" ] && [ "${C3_CUTOFF}" != "none" ] && [ -n "${C3_CUTOFF}" ] \
    && [ "${FIT_ORDER}" -ge 3 ] && C_FLAG="--c3 ${C3_CUTOFF}"
W_FLAG="-w ${FIT_ORDER}"
CLI_METHOD="${FIT_METHOD}"; [ "${FIT_METHOD}" = "RFE" ] && CLI_METHOD="LASSO"   # RFE 走 env 重定向
# 注意：这里【不再】带 --rasr —— RASR 只在 -c（零空间构造）步生效，放 -f 上是死参数。
FIT_FLAGS="--full_ifc -l ${CLI_METHOD} --hdf5"
if [ "${FIT_METHOD}" = "LASSO" ]; then
    # P0-3：列标准化 --std 必须开。不同阶力常数量级差很大，ℓ1 惩罚对系数尺度敏感 ——
    #   不标准化时 fc3 幅值系统性偏低（Mg2C60 实测 fc3max 3.64 vs 参考 37.04，低 90%），
    #   κ 相应偏高；开 --std + 去偏后 37.00（−0.13%）。--std 也决定 alpha 网格的尺度
    #   （derive_alpha_grid 收到 standardize=True，按标准化空间的 alpha_max 锚定）。
    #   原来写的 --mu_min/--mu_max 在 --alpha_auto（pheasy 默认）打开时是死参数，
    #   只在自动锚定失败时兜底，这里保留作兜底并显式写 --alpha_auto 表明意图。
    export PHEASY_LASSO_DEBIAS=1   # LS 选支撑集 + 支撑集上 OLS 去偏（pheasy 默认已开，显式钉住）
    if [ "${FIT_ORDER}" -eq 2 ]; then FIT_FLAGS="${FIT_FLAGS} --mu_min -8 --mu_max 0"
    else FIT_FLAGS="${FIT_FLAGS} --mu_min -8 --mu_max -5"; fi
    FIT_FLAGS="${FIT_FLAGS} --std --alpha_auto --alpha_decades 4.0 --max_iter 100000 --cv 5 --nmu 40 --tol 0.00001"
    echo "LASSO：--std 列标准化 + 去偏(PHEASY_LASSO_DEBIAS=${PHEASY_LASSO_DEBIAS}) + alpha_auto 锚定网格"
elif [ "${FIT_METHOD}" = "RFE" ]; then
    FIT_FLAGS="${FIT_FLAGS} --mu_min -8 --mu_max -5 --max_iter 1000 --cv 5 --nmu 5 --tol 0.001"
fi

# ===== pheasy 四步：cluster space → 对称约束 → 位移矩阵 → 拟合 =====
echo "【pheasy】阶次=${FIT_ORDER} 方法=${FIT_METHOD} ndata=${NDATA} C_FLAG='${C_FLAG}' RASR=${RASR}"
rm -f fc2.hdf5 fc3.hdf5 fc4.hdf5

export OPENBLAS_NUM_THREADS=${NCPU_BLAS} OMP_NUM_THREADS=${NCPU_BLAS} MKL_NUM_THREADS=${NCPU_BLAS}
pheasy --dim ${DIM} ${W_FLAG} -s ${C_FLAG} --eps ${NULL_SPACE_EPS}
# ★ RASR 必须挂在这一步（-c 构造零空间）。挂到 -f 上是死参数（-f 直接读 ns_*.npz）。
pheasy --dim ${DIM} ${W_FLAG} -c ${C_FLAG} --eps ${NULL_SPACE_EPS} ${RASR_FLAGS} 2>&1 | tee pheasy_c.log
_c_rc=${PIPESTATUS[0]}
if [ "${_c_rc}" -ne 0 ]; then echo "❌ pheasy -c 失败 rc=${_c_rc}" >&2; exit 1; fi

# RASR 守卫：要求施加旋转不变性/平衡条件却没在日志里看到 → 直接失败。
#   理由：RASR 生效与否决定 ZA 是 ω∝q² 还是线性（后者频率可能全为正、能过虚频闸，
#   但 κ 是错的），这种错误下游看不出来，只能在源头拦住。
if [ -n "${RASR_FLAGS}" ]; then
    if grep -q "Imposing rotational invariance\|Imposing equilibrium conditions" pheasy_c.log; then
        echo "✅ RASR=${RASR} 已在零空间构造步施加（pheasy_c.log）"
        grep -n "Imposing rotational invariance\|Imposing equilibrium conditions" pheasy_c.log | tail -3
    else
        echo "❌ RASR=${RASR} 但 pheasy -c 日志里没有施加记录 → 力常数不可信，作业失败" >&2
        echo "   常见原因：环境里的 pheasy 太旧（symmetry_constraints.py 没打 RASR 补丁）；" >&2
        echo "   或 --rasr 没传到 -c 步。请核对 pheasy_c.log 与 pheasy --help 的 choices。" >&2
        tail -30 pheasy_c.log >&2 || true
        exit 1
    fi
fi

# ===== 数据量闸（review #1）：方程数 / 参数数 ≥ 3 =====
# pheasy -c 已构造完参数空间，pheasy_c.log 里有实际 free IFC 数；-f 之前先判，
# 欠定就停下并给出至少需要的帧数（WS₂ 实测 12 帧只有 1.72，需 ≥21 帧）。
_NATOM=$(cat natom_super.txt)
_FREE_IFC=$(sed -n 's/.*Total number of free IFCs:[[:space:]]*\([0-9]*\).*/\1/p' pheasy_c.log | tail -1)
if [ -z "${_FREE_IFC}" ] || [ "${_FREE_IFC}" = "0" ]; then
    echo "⚠ 未能从 pheasy_c.log 读到 free IFC 数，跳过量纲闸" >&2
else
    _EQ=$(( ${NDATA} * 3 * ${_NATOM} ))
    _NEED=$(( ( ${_FREE_IFC} + ${_NATOM} - 1 ) / ${_NATOM} ))
    _RATIO_PCT=$(( ${_EQ} * 100 / ${_FREE_IFC} ))
    if [ "${_RATIO_PCT}" -lt 300 ]; then
        echo "❌ 数据量不足：方程数/参数数 = ${_EQ}/${_FREE_IFC} ≈ ${_RATIO_PCT}/100 < 3。" >&2
        echo "   当前 ${NDATA} 帧，至少需要 ${_NEED} 帧。请加大 S4 帧数（OVERSAMPLE/截断），或缩小 S5 截断。" >&2
        exit 1
    fi
    echo "✅ 数据量：方程数/参数数 = ${_EQ}/${_FREE_IFC} ≈ ${_RATIO_PCT}/100（≥3 通过；至少需 ${_NEED} 帧）"
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PHEASY_N_JOBS=${NCPU_DISP}
pheasy --dim ${DIM} ${W_FLAG} -d ${C_FLAG} --ndata ${NDATA} --disp_file --eps ${NULL_SPACE_EPS}

export LOKY_MAX_CPU_COUNT=${NCPU_LOKY} OPENBLAS_NUM_THREADS=${NCPU_FIT_BLAS}
export OMP_NUM_THREADS=${NCPU_FIT_BLAS} MKL_NUM_THREADS=${NCPU_FIT_BLAS}
export PHEASY_N_JOBS=${NCPU_LOKY} PHEASY_DOT_THREADS=${NCPU} OMP_NESTED=FALSE MKL_DYNAMIC=FALSE
pheasy --dim ${DIM} ${W_FLAG} -f ${C_FLAG} --ndata ${NDATA} --eps ${NULL_SPACE_EPS} ${FIT_FLAGS} 2>&1 | tee pheasy_f.log
_f_rc=${PIPESTATUS[0]}
if [ "${_f_rc}" -ne 0 ]; then echo "❌ pheasy -f 失败 rc=${_f_rc}" >&2; exit 1; fi

# ===== 输出检查 =====
sync
[ ! -f fc2.hdf5 ] && echo "❌ 缺 fc2.hdf5" && exit 1
[ "${FIT_ORDER}" -ge 3 ] && [ ! -f fc3.hdf5 ] && echo "❌ 缺 fc3.hdf5" && exit 1
ls -lh fc2.hdf5 fc3.hdf5 2>/dev/null
fi  # CUT3_CANDIDATES 分支

# ===== 收尾：搬产物 + shengbte 导出 + 虚频闸 =====
python kl_fc_backends.py collect_pheasy fit_config.json

# ===== 拟合质量指标留档（A/B 对照用）=====
#   pheasy 自报的 RMSE / Relative error / alpha，加上力常数量级、零空间自由度，
#   统一落 fit_metrics.json；post 步会并进 phonon_summary.json。
#   动力：2026-09-16 做 RASR=none vs BHH 的 A/B 时发现这些量一个都没记，
#   只能靠翻日志，无法两臂对照。
python - <<'PYEOF'
import json, re, os, glob
import numpy as np
m = {}
log = ""
for _f in ('pheasy_f.log', 'pheasy_c.log'):
    if os.path.isfile(_f):
        log += open(_f, encoding='utf-8', errors='ignore').read()
for _k, _p, _c in (
        ('pheasy_rmse_eV_per_A', r'\bRMSE:\s*([\d.eE+-]+)', float),
        ('pheasy_relative_error', r'Relative error:\s*([\d.eE+-]+)', float),
        ('pheasy_worst_force_correlation', r'worst corr=([\d.]+)', float),
        ('pheasy_free_ifcs', r'Free IFC terms:\s*(\d+)', int),
        ('pheasy_best_alpha', r'best alpha=\s*([\d.eE+-]+)', float),
        ('pheasy_alpha_opt', r'alpha_opt:[ ]*([0-9.eE+-]+)', float),
        ('pheasy_alpha_min', r'alpha_min:[ ]*([0-9.eE+-]+)', float),
        ('pheasy_alpha_max', r'alpha_max:[ ]*([0-9.eE+-]+)', float)):
    _mm = re.search(_p, log)
    if _mm:
        try:
            m[_k] = _c(_mm.group(1))
        except ValueError:
            pass
# RASR 是否真在 -c 步施加（与模板里的守卫同一判据，这里留证）
m['rasr_applied'] = bool(re.search(r'Imposing rotational invariance|Imposing equilibrium conditions', log))
m['rasr_requested'] = os.environ.get('RASR', '') 
m['ns_rank_tol'] = os.environ.get('PHEASY_NS_RANK_TOL', '')
def _h5max(p):
    if not os.path.isfile(p):
        return None
    try:
        import h5py
        with h5py.File(p) as h:
            v = [float(np.abs(np.array(d)).max()) for d in h.values() if hasattr(d, 'shape')]
        return max(v) if v else None
    except Exception:
        try:
            return float(np.abs(np.load(p)).max())
        except Exception:
            return None
for _tag, _p in (('fc2', 'phono3py/fc2.hdf5'), ('fc2', 'fc2.hdf5'),
                 ('fc3', 'phono3py/fc3.hdf5'), ('fc3', 'fc3.hdf5')):
    _v = _h5max(_p)
    if _v is not None and ('%smax' % _tag) not in m:
        m['%smax' % _tag] = _v
# 零空间自由度：-c 步产出的 ns_*.npz
for _p in sorted(glob.glob('ns_*.npz') + glob.glob('phono3py/ns_*.npz')):
    try:
        _z = np.load(_p)
        m['%s_shape' % os.path.basename(_p)[:-4]] = [list(np.shape(_z[k])) for k in _z.files]
        m['ns_harm_free_params' if 'harm' in _p else os.path.basename(_p)[:-4] + '_ndof'] = \
            int(sum(int(np.prod(np.shape(_z[k])[1:])) if np.ndim(_z[k]) > 1 else 1 for k in _z.files))
    except Exception as _e:
        m['%s_error' % os.path.basename(_p)] = str(_e)
open('fit_metrics.json', 'w').write(json.dumps(m, ensure_ascii=False, indent=2))
print('[OK] fit_metrics.json：%s' % json.dumps(
      {k: m[k] for k in ('pheasy_relative_error', 'pheasy_rmse_eV_per_A',
                         'fc2max', 'fc3max', 'rasr_applied') if k in m}, ensure_ascii=False))
PYEOF

python kl_fc_backends.py post           fit_config.json
