"""
amset_ir_fix.py —— AMSET 密网格不可约 k 点改用正确的晶格（运行时插件，不改 AMSET 源码）

用法（S8 / S8.4 的 gen 已把本文件复制进运行目录，并在 python -c 里最先 import）：
    AZ_IR_FIX=1 python -c "import amset_ir_fix; ...; Runner.from_directory('.').run()"   （V130 起必须显式设 1）

开关：环境变量 AZ_IR_FIX（1/true/on/yes = 修；未设 / 0/false/off/no = AMSET 原行为；V130 起未设 = 关）。
    gen 按 step.conf 的 IR_FIX（auto/on/off，auto = on）在命令前写 export AZ_IR_FIX=…。

------------------------------------------------------------------------------
问题（V125）：amset/electronic_structure/kpoints.py::get_kpoints_tetrahedral 给 spglib 的是
    cell = (np.array(atoms.get_cell().T, ...),   # lattice
            atoms.get_scaled_positions(), ...)
ASE 的 get_cell() 按行存晶格矢量，spglib 的 cell[0] 也要求按行 —— 多出来的 .T 等于把晶格换成了
另一个。spglib 于是只找得到两种晶格都容许、又把原子映回原子的操作：
  · 晶格矩阵对称（立方、fcc/bcc 原胞、正交/四方常规胞）-> 没影响；
  · 六方/三方（2D TMD、纤锌矿 GaN/AlN…）-> 只剩 {E, σh}×时间反演。
找到的都是真对称操作，所以**结果是对的**，只是不可约 k 点多了 5.5–5.9 倍（实测：WSe2 单层
263×263×33 587945 vs 100232；GaN 61×61×37 35359 vs 6479）。散射按不可约点逐个算，于是白算 5–6 倍。
合成六方单层 AMSET 端到端（ADP+IMP，85×85×11）：不可约点 21678 -> 3870，散射 121 s -> 25 s，
迁移率一致到 0.05%。

修法：只替换 kpoints 模块里 get_ir_reciprocal_mesh 的调用 —— 把收到的 lattice 再转置一次（= 行矢量）。
该模块只有 get_kpoints_tetrahedral 这一处调用 spglib（护栏里核对）。
护栏：源码里必须能找到 "atoms.get_cell().T" 且 spglib 只被调用这一处；对不上（AMSET 已修或改了
实现）-> 不打补丁、照原样运行（原行为结果是对的，只是慢），并在 stderr 说明。

V127 更正：上面"迁移率一致"只对协变的散射核成立。AMSET 的 ADP 核（|D| 逐元素、分数取整的 q、
D 的线性插值）不协变，张量形变势下代表点≠成员平均（MoS₂ n 型 ADP −4.28%）。所以 IR_FIX 开时同时挂
"ADP 星平均"（见下文 V127 段；AZ_ADP_STAR=0 关）。
------------------------------------------------------------------------------
"""

import inspect
import json
import os
import sys

import numpy as np

ENV_NAME = "AZ_IR_FIX"
_OFF = ("0", "false", "off", "no")
BUG_SNIPPET = "atoms.get_cell().T"
STATE = {"applied": False, "calls": 0, "last": None}


_ON = ("1", "true", "on", "yes")


def enabled():
    """V130：未设 = 关（与 gen 的出厂默认一致）；只有显式 1/true/on/yes 才打补丁。"""
    return os.environ.get(ENV_NAME, "").strip().lower() in _ON


def _squash(s):
    return "".join(str(s).split())


class _SpglibProxy:
    """只改 get_ir_reciprocal_mesh：lattice 转回行矢量；其余属性原样转发。"""

    def __init__(self, real):
        self._real = real

    def __getattr__(self, name):
        return getattr(self._real, name)

    def get_ir_reciprocal_mesh(self, mesh, cell, *args, **kwargs):
        lat = np.ascontiguousarray(np.asarray(cell[0], dtype="double").T)
        fixed = (lat,) + tuple(cell[1:])
        mapping, grid = self._real.get_ir_reciprocal_mesh(mesh, fixed, *args, **kwargs)
        STATE["calls"] += 1
        try:
            # spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift, is_time_reversal, symprec, ...)
            pos = dict(zip(("is_shift", "is_time_reversal", "symprec"), args))
            pos.update(kwargs)
            _publish_star_rotations(self._real, fixed, cell, pos.get("symprec", 1e-5),
                                    pos.get("is_time_reversal", True))
        except Exception as e:                           # noqa: BLE001  星平均失败 -> 退回代表点（照 V125）
            os.environ.pop(STAR_ENV, None)
            print("[ir_fix][WARN] ADP 星平均准备失败（%s）-> 只用代表点" % e, file=sys.stderr, flush=True)
        try:
            m0, _ = self._real.get_ir_reciprocal_mesh(mesh, cell, *args, **kwargs)
            n_old = int(len(np.unique(m0)))
        except Exception:                                  # noqa: BLE001  只用于日志
            n_old = None
        n_new = int(len(np.unique(mapping)))
        STATE["last"] = {"mesh": [int(x) for x in mesh], "n_full": int(np.prod(mesh)),
                         "ir_fixed": n_new, "ir_amset_original": n_old}
        print("[ir_fix] mesh %s：不可约 k 点 %d（AMSET 原式 %s，少算 x%s）"
              % ("x".join(str(int(x)) for x in mesh), n_new, n_old,
                 ("%.1f" % (n_old / n_new)) if n_old and n_new else "?"),
              file=sys.stderr, flush=True)
        return mapping, grid


# ------------------------------------------------------------------------------
# V127：ADP 星平均 —— IR_FIX 只该改变"算哪些 k 点"，不该改变结果
#
# MoS₂ 实测（IR_FIX 开、KZ_CAP 关）：n 型 ADP 面内平均 −4.28%、p 型 −0.04%，xx/yy 各向异性 1.98% -> 0。
# 原因不在 IR_FIX：AMSET 的 ADP 核有三处不随旋转协变
#     deform = np.abs(D(k)) + v⊗v        逐元素取绝对值，|R D Rᵀ| ≠ R |D| Rᵀ
#     q = pbc_diff(...)                  分数坐标取整定像（六方的分数胞不是 Wigner–Seitz 胞）
#     D(k)                               deformation.h5 粗网格上线性插值
# （amset/scattering/elastic.py、calculate.py；amset2d_plugin 照抄前两处），同一颗星各成员的率本来就不等。
# 旧跑法（转置晶格，只剩 {E, σh}×时间反演）把星里的成员分别算 —— 等于在星上做了平均，于是 xx≠yy；
# IR_FIX 只算代表点，结果就取决于 spglib 挑的是哪个成员。合成六方单层（各向异性形变势场，
# proto3）：n 型 ±19% 量级的差都能造出来。
#
# 修法：代表点上把 ADP 被积因子对"新群比旧群多出来的那些操作"（右陪集 H\G 的代表 S，含时间反演）
# 取平均：  f̄(q) = mean_S f(q̂_S, D(S·k), S v(k))，q_S = S q 按 AMSET 的分数取整重新取像，
# D 在星成员 S·k 上重新插值 —— 就是旧跑法在成员 S·k 上算的那个被积函数（v 由 FFT 得到、本来协变；
# 弹性张量在点群下不变）。残差只剩"对率平均"与"对 τ 平均"之差（二阶，Jensen）。
# 立方等晶格矩阵对称的体系 G = H，只有一个陪集 -> 什么都不做（与原 AMSET 逐位相同）。
# 父进程在建密网格时把陪集旋转写进环境变量 AZ_IR_FIX_STAR（JSON），spawn 子进程读；
# AZ_ADP_STAR=0 关掉（A/B 用）。只动张量形变势（deformation.h5）；标量形变势本来协变。
# ------------------------------------------------------------------------------
STAR_ENV = "AZ_IR_FIX_STAR"
STAR_SWITCH = "AZ_ADP_STAR"
_STAR = {"rots": None, "loaded": False}


def _quiet():
    try:
        import multiprocessing
        # spawn 出来的散射子进程不刷屏。子进程在反序列化 worker 入口时导入本模块，那时 parent_process()
        # 还没设好（_bootstrap 之后才有），但 current_process()._inheriting 已经是 True。
        return (multiprocessing.parent_process() is not None
                or bool(getattr(multiprocessing.current_process(), "_inheriting", False)))
    except Exception:                                              # noqa: BLE001
        return False


def cart_rotations(lattice, rotations):
    """spglib 的分数坐标旋转 W（作用在原子分数坐标）-> 笛卡尔 R = Lᵀ W L⁻ᵀ（L 按行存晶格矢量）。"""
    Lt = np.asarray(lattice, float).T
    Li = np.linalg.inv(Lt)
    return [Lt @ np.asarray(W, float) @ Li for W in rotations]


def _uniq(mats, tol=1e-6):
    out = []
    for m in mats:
        if not any(np.abs(m - o).max() < tol for o in out):
            out.append(np.asarray(m, float))
    return out


def coset_rotations(G, H, time_reversal=True, tol=1e-6):
    """k 空间点群 G±、H±（时间反演 = -E）的右陪集 H·S 的代表 S。
    旧 AMSET 每个 H 轨道算一个率，所以旧结果 = 对这些代表的平均。"""
    def pm(X):
        X = _uniq(X, tol)
        return _uniq(X + [-x for x in X], tol) if time_reversal else X
    Gp, Hp = pm(G), pm(H)
    reps = []
    for S in Gp:
        if not any(any(np.abs(S - h @ R).max() < tol for h in Hp) for R in reps):
            reps.append(S)
    return reps


def _publish_star_rotations(real_spglib, fixed, orig, symprec, time_reversal):
    if os.environ.get(STAR_SWITCH, "").strip().lower() in _OFF:
        os.environ.pop(STAR_ENV, None)
        print("[ir_fix] %s=off -> ADP 只用代表点（不做星平均）" % STAR_SWITCH, file=sys.stderr, flush=True)
        return
    L = np.asarray(fixed[0], float)
    G = cart_rotations(L, real_spglib.get_symmetry(fixed, symprec=symprec)["rotations"])
    # 旧式（转置晶格）找到的分数操作，作用在**真实**晶体上的笛卡尔形式
    H = cart_rotations(L, real_spglib.get_symmetry(orig, symprec=symprec)["rotations"])
    bad = [R for R in G + H if np.abs(R.T @ R - np.eye(3)).max() > 1e-3]
    if bad:
        raise ValueError("有 %d 个操作在真实晶格上不正交" % len(bad))
    reps = coset_rotations(G, H, bool(time_reversal))
    if len(reps) <= 1:
        os.environ.pop(STAR_ENV, None)
        return
    os.environ[STAR_ENV] = json.dumps([np.round(R, 12).tolist() for R in reps])
    _STAR["loaded"] = False
    print("[ir_fix] ADP 星平均：%d 个陪集（新点群 %d 个操作 / AMSET 原式 %d，不计时间反演）"
          % (len(reps), len(_uniq(G)), len(_uniq(H))), file=sys.stderr, flush=True)


def star_rotations():
    if not _STAR["loaded"]:
        _STAR["loaded"] = True
        s = os.environ.get(STAR_ENV)
        _STAR["rots"] = [np.array(r, float) for r in json.loads(s)] if s else None
    return _STAR["rots"]


def reimage_q(q, rlat):
    """AMSET 在每个起点 k 上用分数坐标取整（pymatgen pbc_diff：f - round(f)）定 q 的像，不是笛卡尔最短像。
    旋转到星成员后重新取整，才是 AMSET 在那个成员上会用的 q。q (n,3) 笛卡尔 Å⁻¹，rlat 按行（含 2π）。"""
    f = np.asarray(q) @ np.linalg.inv(rlat)
    return (f - np.round(f)) @ rlat


def star_average_factor(scatterer, rots, rlat=None):
    """把 ADP 实例的 factor 换成星平均：
        f̄ = mean_S f(q̂_S, D(S·k), S v)，q_S = S q 按 AMSET 的分数取整重新取像，D(S·k) 在星成员上重新插值
        （rlat 为 None 时退回 S D(k) Sᵀ、不重取像）。
    AMSET 核里的 |D|（逐元素）、分数取整的 q、粗网格上线性插值的 D 都不随旋转协变 —— 所以要逐个成员算。
    幂等；标量形变势（本来协变）不动。"""
    if getattr(scatterer, "_az_star", False):
        return False
    scatterer._az_star = True
    dp = getattr(scatterer, "deformation_potential", None)
    if not hasattr(dp, "interpolate") or not rots or len(rots) < 2:
        return False
    orig = scatterer.factor
    rlat = None if rlat is None else np.asarray(rlat, float)

    interp = dp.interpolate
    rinv = None if rlat is None else np.linalg.inv(rlat)

    def factor(unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        uq, nq, v = np.asarray(unit_q), np.asarray(norm_q_sq), np.asarray(velocity)
        q = uq * np.sqrt(nq)[:, None]
        if rlat is None:
            D0 = np.asarray(interp(spin, [band_idx], [kpoint]))
        else:
            kc = np.asarray(kpoint, float) @ rlat
        acc = None
        try:
            for R in rots:
                if rlat is None:
                    DR = np.einsum("ij,njk,lk->nil", R, D0, R)
                else:
                    # 在星成员 S·k 上重新插值（与 AMSET 在那个成员上做的一样；粗网格线性插值也不协变）
                    kS = (kc @ R.T) @ rinv
                    DR = np.asarray(interp(spin, [band_idx], [kS - np.round(kS)]))
                dp.interpolate = (lambda *_a, _D=DR, **_k: _D)
                if rlat is None or np.allclose(R, np.eye(3)):
                    uR, nR = uq @ R.T, nq
                else:
                    qR = reimage_q(q @ R.T, rlat)
                    nR = np.sum(qR ** 2, axis=-1)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        uR = qR / np.sqrt(nR)[:, None]
                out = orig(uR, nR, spin, band_idx, kpoint, R @ v)
                acc = out if acc is None else acc + out
        finally:
            dp.__dict__.pop("interpolate", None)
        return acc / len(rots)

    scatterer.factor = factor
    return True


def worker(*args, **kwargs):
    """散射子进程入口：spawn 按引用反序列化本函数 -> 子进程会导入本模块（星平均在子进程里生效）。"""
    return _ORIG["worker"](*args, **kwargs)


_ORIG = {"worker": None, "calculate_rate": None}


def _wrap_scattering():
    import amset.scattering.calculate as C
    if not getattr(C.calculate_rate, "_az_star", False):
        _ORIG["calculate_rate"] = C.calculate_rate

        def calculate_rate(tbs, overlap_calculator, mrta_calculator, elastic_scatterers, *a, **k):
            rots = star_rotations()
            if rots:
                try:
                    rlat = a[1].structure.lattice.reciprocal_lattice.matrix       # amset_data_min
                except (AttributeError, IndexError):
                    rlat = None
                for s in elastic_scatterers:
                    if getattr(s, "name", "") == "ADP":
                        star_average_factor(s, rots, rlat)
            return _ORIG["calculate_rate"](tbs, overlap_calculator, mrta_calculator, elastic_scatterers, *a, **k)

        calculate_rate._az_star = True
        C.calculate_rate = calculate_rate
    if C.scattering_worker is not worker and _ORIG["worker"] is None:
        _ORIG["worker"] = C.scattering_worker
        C.scattering_worker = worker


def apply(force=False):
    """给 amset.electronic_structure.kpoints 换上 spglib 代理。返回是否打了补丁。"""
    if not (force or enabled()):
        if not _quiet():
            print("[ir_fix] %s=off -> 不打补丁（AMSET 原不可约 k 点）" % ENV_NAME, file=sys.stderr, flush=True)
        return False
    import amset
    import amset.electronic_structure.kpoints as _kp

    if isinstance(_kp.spglib, _SpglibProxy):
        _wrap_scattering()
        return True
    try:
        fsrc = inspect.getsource(_kp.get_kpoints_tetrahedral)
        msrc = inspect.getsource(_kp)
    except (OSError, TypeError) as e:
        print("[ir_fix] 读不到 AMSET 源码（%s）-> 不打补丁（原行为：结果对，只是慢）" % e,
              file=sys.stderr, flush=True)
        return False
    n_calls = _squash(msrc).count("spglib.")
    if _squash(BUG_SNIPPET) not in _squash(fsrc) or n_calls != 1:
        print("[ir_fix] amset %s 的 get_kpoints_tetrahedral 与已核对的实现不同（%s 出现=%s，spglib 调用 %d 处）"
              " -> 不打补丁（原行为：结果对，只是可能慢）"
              % (getattr(amset, "__version__", "?"), BUG_SNIPPET, _squash(BUG_SNIPPET) in _squash(fsrc), n_calls),
              file=sys.stderr, flush=True)
        return False
    _kp.spglib = _SpglibProxy(_kp.spglib)
    STATE["applied"] = True
    _wrap_scattering()
    if not _quiet():
        print("[ir_fix] amset %s：密网格不可约 k 点改用正确晶格（行矢量）" % getattr(amset, "__version__", "?"),
              file=sys.stderr, flush=True)
    return True


apply()
