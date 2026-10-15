"""
amset_desym_fix.py —— AMSET 去对称化相位修正（运行时插件，不改 AMSET 源码）

用法（在 AMSET 运行目录内；S8 / S8.4 的 gen 已把本文件复制进运行目录）：
    AZ_DESYM_FIX=1 python -c "import amset_desym_fix; \
        from amset.log import initialize_amset_logger as L; L(); \
        from amset.core.run import Runner; Runner.from_directory('.').run()"

开关：环境变量 AZ_DESYM_FIX（1/true/on/yes = 打补丁；未设或其它值 = **不打补丁、AMSET 原行为**）。
    验证期默认关。gen 按 step.conf 的 DESYM_FIX 决定是否在命令前加 `export AZ_DESYM_FIX=1`。

------------------------------------------------------------------------------
问题（V115）：amset/wavefunction/common.py::desymmetrize_coefficients 第 126 行
    factor = np.exp(-1j * 2 * np.pi * np.dot(rot_gpoints + rot_kpoint, tau))
AMSET 在倒空间用的旋转是 M = Rᵀ（spglib 实空间分数坐标 {R|τ} 的转置），把 k 映到 Mk 的
实空间操作是 {R|τ}⁻¹ = {R⁻¹ | −R⁻¹τ}。平面波系数的正确变换为
    c_{Mk}(MG) = c_k(G) · exp(+2πi (k+G)·τ)                      （非 TR 操作）
    c_{−Rᵀk}(−RᵀG) = conj[ c_k(G) · exp(+2πi (k+G)·τ) ]           （TR 操作；代码里 tau 已翻号）
即相位里应是**未旋转的 k+G**、符号为 +。原式等价于 exp(−2πi (k+G)·Rτ)，与正确值差
    exp(2πi G·(R+I)τ)（非 TR）   /   exp(2πi G·(R−I)τ)（TR），
因此原式精确 ⇔ 非 TR：(R+I)τ ≡ 0，TR：(R−I)τ ≡ 0（mod 1）—— 就是 ke_common.op_phase_exact
的逐操作判据。与 k 有关的那部分只是每个 k 的整体相位，不影响 |<ψ|ψ'>|。

修法（一行）：
    factor = np.exp(+1j * 2 * np.pi * np.dot(gpoints + kpoint, tau))
其中 tau 取 `if tr: tau = -tau` 之后的值（与原代码位置相同）。

挂法：desymmetrize_coefficients 是 amset.wavefunction.common 的模块级函数，被
amset.interpolation.wavefunction 按名字导入（from ... import desymmetrize_coefficients），
在 WavefunctionOverlapCalculator.from_coefficients 里调用。所以两个模块的名字**都要替换**，
只换一个不生效。0.4.19 与 0.5.1 的这个函数逐字节相同（已 diff 核对），同一份补丁两边通用。

SOC（非共线，系数最后一维是 2）分支另有两处错，V200（2026-10-15）一并修掉：
  ② 自旋转动：AMSET 用 rotation_matrix_to_su2 —— get_rotation_angle 只给无符号转角、
     get_rotation_axis 把轴强制翻到 z>=0 半球，两者组合约一半操作的 SU(2) 转反了方向
     （Si 48 个操作里 40 个与正确 SU(2) 不等价）。改用四元数构造 su2_from_rotation。
  ③ 时间反演：AMSET 的 ncl 分支对 TR 操作既不取复共轭也不乘 iσ_y（第 145 行只对非 ncl 取共轭）。
     改为先共轭再 (c_up, c_dn) -> (c_dn, −c_up)，即 Θ = iσ_y·K。
  检验：SOC 旋量紧束缚模型（Kane–Mele 型 iλ(d1×d2)·σ，天然满足空间群 + TR），12 个构型
  （test_symmetry_gate.CASES）修正式最小 |cos| >= 0.999999，原 AMSET 式在 14–43 个 k 点上出错（多数 |cos|≈0）；
  有反演、无 TR 操作的 Si/hcp 原式同样出错（②的后果）。**SOC 实数据尚未验证**。

版本护栏：打补丁前核对原函数源码里确实有上面那行原式；对不上（AMSET 改了实现）就**抛异常
  停跑**，绝不静默 —— 因为开了 DESYM_FIX 之后判据会放行 IBZ，补丁没挂上就是静默算错。
------------------------------------------------------------------------------
"""

import inspect
import logging
import os
import sys
import time

import numpy as np

ENV_NAME = "AZ_DESYM_FIX"
_TRUE = ("1", "true", "on", "yes")
# 原函数里必须能找到的那一行（去掉空白后比较）
ORIG_FACTOR_LINE = "factor = np.exp(-1j * 2 * np.pi * np.dot(rot_gpoints + rot_kpoint, tau))"
FIXED_FACTOR_LINE = "factor = np.exp(1j * 2 * np.pi * np.dot(gpoints + kpoint, tau))"
OP_TOL = 1e-4          # 与 ke_common.OP_TOL 同值（只用于日志统计，不影响计算）

logger = logging.getLogger("amset.desym_fix")   # amset 的子 logger -> 进 amset.log
STATE = {"applied": False, "calls": 0, "last": None}


def enabled():
    return os.environ.get(ENV_NAME, "").strip().lower() in _TRUE


def _squash(s):
    return "".join(str(s).split())


def _op_exact(rot_recip, tau_amset, tr, atol=OP_TOL):
    """AMSET 原式在这个操作上是否精确（换回实空间 R、τ 后套逐操作判据）。

    AMSET 数组：非 TR 为 (Rᵀ, τ)，TR 为 (−Rᵀ, −τ)。
    """
    R = (-rot_recip if tr else rot_recip).T
    tau = -tau_amset if tr else tau_amset
    v = (R - np.eye(3)) @ tau if tr else (R + np.eye(3)) @ tau
    return bool(np.allclose(v - np.rint(v), 0.0, atol=atol))


def su2_from_rotation(rot_cart):
    """笛卡尔 3×3（真/非真）转动 -> SU(2)：U = q0·I − i(q1σx + q2σy + q3σz)。

    非真转动先乘反演（反演不作用于自旋）。用四元数（Shepperd 法）直接构造，带符号的转角与
    转轴一次给出 —— 替代 AMSET 的 rotation_matrix_to_su2：后者的 get_rotation_angle 只给
    无符号转角、get_rotation_axis 又把轴强制翻到 z≥0 半球，两者组合在约一半的操作上给出
    反向转动（模型检验：Si 48 个操作里 40 个不对）。整体 ±1 不影响任何 |重叠|。
    """
    R = np.asarray(rot_cart, float)
    if np.linalg.det(R) < 0:
        R = -R
    t = np.trace(R)
    if t > 0:
        s = 2.0 * np.sqrt(1.0 + t)
        q = np.array([0.25 * s, (R[2, 1] - R[1, 2]) / s,
                      (R[0, 2] - R[2, 0]) / s, (R[1, 0] - R[0, 1]) / s])
    else:
        i = int(np.argmax(np.diag(R)))
        j, k = (i + 1) % 3, (i + 2) % 3
        s = 2.0 * np.sqrt(max(1.0 + R[i, i] - R[j, j] - R[k, k], 0.0))
        q = np.zeros(4)
        q[0] = (R[k, j] - R[j, k]) / s
        q[1 + i] = 0.25 * s
        q[1 + j] = (R[j, i] + R[i, j]) / s
        q[1 + k] = (R[k, i] + R[i, k]) / s
    q = q / np.linalg.norm(q)
    return np.array([[q[0] - 1j * q[3], -q[2] - 1j * q[1]],
                     [q[2] - 1j * q[1], q[0] + 1j * q[3]]], dtype=np.complex128)


def op_audit(rotations, translations, is_tr, op_mapping):
    """统计本次去对称化**实际用到**的操作里，原式会算错的有几个、影响多少个 k 点。"""
    rotations = np.asarray(rotations)
    translations = np.asarray(translations, float)
    is_tr = np.asarray(is_tr, bool)
    op_mapping = np.asarray(op_mapping, int)
    exact = np.array([_op_exact(rotations[i], translations[i], bool(is_tr[i]))
                      for i in range(len(rotations))], bool)
    used = np.unique(op_mapping)
    bad_used = [int(i) for i in used if not exact[i]]
    return {
        "total_ops": int(len(rotations)),
        "bad_ops": int((~exact).sum()),
        "used_ops": int(len(used)),
        "bad_used_ops": len(bad_used),
        "tr_ops_used": int(is_tr[used].sum()),
        "nk_full": int(len(op_mapping)),
        "nk_bad": int((~exact[op_mapping]).sum()),
    }


def desymmetrize_coefficients(
    coeffs,
    gpoints,
    kpoints,
    structure,
    rotations,
    translations,
    is_tr,
    op_mapping,
    kp_mapping,
    pbar=True,
):
    """amset.wavefunction.common.desymmetrize_coefficients 的逐行拷贝：改 factor 一行（相位），
    SOC 分支另改 SU(2) 构造与旋量 TR（见模块说明 ②③）。"""
    from amset.electronic_structure.symmetry import similarity_transformation
    from amset.log import log_time_taken
    from amset.util import get_progress_bar
    from amset.wavefunction.common import is_ncl

    logger.info("Desymmetrizing wavefunction coefficients "
                "[AutoZT phase fix: exp(+2*pi*i*(k+G).tau)]")
    t0 = time.perf_counter()

    ncl = is_ncl(coeffs)
    _audit = op_audit(rotations, translations, is_tr, op_mapping)
    _audit["ncl"] = bool(ncl)
    STATE["calls"] += 1
    STATE["last"] = _audit
    logger.info(
        "  desym_fix: ops %d, original-formula-inexact ops %d/%d "
        "(used %d/%d, affecting %d/%d full-grid k-points), TR ops used %d"
        % (_audit["total_ops"], _audit["bad_ops"], _audit["total_ops"],
           _audit["bad_used_ops"], _audit["used_ops"], _audit["nk_bad"],
           _audit["nk_full"], _audit["tr_ops_used"]))
    if ncl:
        logger.warning(
            "  ★ desym_fix: SOC/non-collinear coefficients: phase fix + quaternion SU(2) + "
            "spinor time reversal (i*sigma_y*K, %d TR ops used). Validated on a spinor "
            "tight-binding model (12/12 structures); NOT yet validated on real SOC data."
            % _audit["tr_ops_used"])

    rots = rotations[op_mapping]
    taus = translations[op_mapping]
    trs = is_tr[op_mapping]

    su2s = None
    if ncl:
        # get cartesian rotation matrix
        r_cart = [
            similarity_transformation(structure.lattice.matrix.T, r.T)
            for r in rotations
        ]

        # calculate SU(2)
        # ---- AutoZT 修正 ②（SOC）：不用 AMSET 的 rotation_matrix_to_su2（转角无符号、转轴强制翻
        #      到 z>=0 半球，约一半操作给出反向转动），改用四元数构造 -------------------------
        su2_no_dagger = np.array([su2_from_rotation(r) for r in r_cart])

        # calculate SU(2)^{dagger}
        su2 = np.conjugate(su2_no_dagger).transpose((0, 2, 1))
        su2s = su2[op_mapping]

    g_mesh = (np.abs(gpoints).max(axis=0) + 3) * 2
    g1, g2, g3 = (gpoints + g_mesh / 2).astype(int).T  # indices of g-points to keep

    all_rot_coeffs = {}
    for spin, spin_coeffs in coeffs.items():
        coeff_shape = (len(spin_coeffs), len(rots)) + tuple(g_mesh)
        if ncl:
            coeff_shape += (2,)
        rot_coeffs = np.zeros(coeff_shape, dtype=complex)

        state_idxs = list(range(len(rots)))
        if pbar:
            state_idxs = get_progress_bar(state_idxs, desc="progress")

        for k_idx in state_idxs:
            map_idx = kp_mapping[k_idx]

            rot = rots[k_idx]
            tau = taus[k_idx]
            tr = trs[k_idx]
            kpoint = kpoints[map_idx]

            rot_kpoint = np.dot(rot, kpoint)
            kdiff = np.around(rot_kpoint)
            rot_kpoint -= kdiff

            edges = np.around(rot_kpoint, 5) == -0.5
            rot_kpoint += edges
            kdiff -= edges

            rot_gpoints = np.dot(rot, gpoints.T).T
            rot_gpoints = np.around(rot_gpoints).astype(int)
            rot_gpoints += kdiff.astype(int)

            if tr:
                tau = -tau

            # ---- AutoZT 修正（唯一改动）：未旋转的 k+G，符号 + -------------------
            # 原式：factor = np.exp(-1j * 2 * np.pi * np.dot(rot_gpoints + rot_kpoint, tau))
            factor = np.exp(1j * 2 * np.pi * np.dot(gpoints + kpoint, tau))
            rg1, rg2, rg3 = (rot_gpoints + g_mesh / 2).astype(int).T

            if ncl:
                # perform rotation in spin space
                su2 = su2s[k_idx]
                rc = np.zeros_like(spin_coeffs[:, map_idx])
                rc[:, :, 0] = (
                    su2[0, 0] * spin_coeffs[:, map_idx, :, 0]
                    + su2[0, 1] * spin_coeffs[:, map_idx, :, 1]
                )
                rc[:, :, 1] = (
                    su2[1, 0] * spin_coeffs[:, map_idx, :, 0]
                    + su2[1, 1] * spin_coeffs[:, map_idx, :, 1]
                )
                rot_coeffs[:, k_idx, rg1, rg2, rg3] = factor[None, :, None] * rc
            else:
                rot_coeffs[:, k_idx, rg1, rg2, rg3] = spin_coeffs[:, map_idx] * factor

            if tr and not ncl:
                rot_coeffs[:, k_idx] = np.conjugate(rot_coeffs[:, k_idx])
            elif tr and ncl:
                # ---- AutoZT 修正 ③（SOC）：旋量的时间反演 Θ = iσ_y·K —— 先取复共轭，
                #      再 (c_up, c_dn) -> (c_dn, −c_up)。AMSET 原 ncl 分支对 TR 操作两样都没做。
                rc_tr = np.conjugate(rot_coeffs[:, k_idx])
                rot_coeffs[:, k_idx, ..., 0] = rc_tr[..., 1]
                rot_coeffs[:, k_idx, ..., 1] = -rc_tr[..., 0]

        all_rot_coeffs[spin] = rot_coeffs[:, :, g1, g2, g3]

    log_time_taken(t0)
    return all_rot_coeffs


def apply(force=False):
    """把两个模块里的 desymmetrize_coefficients 换成修正版。返回是否打了补丁。

    未开启（且非 force）时什么都不做，返回 False —— AMSET 原行为。
    开启了但原函数对不上（AMSET 改了实现）时抛 RuntimeError（绝不静默）。
    """
    if not (force or enabled()):
        print("[desym_fix] %s 未设 -> 不打补丁（AMSET 原去对称化公式）" % ENV_NAME,
              file=sys.stderr, flush=True)
        return False
    import amset
    import amset.interpolation.wavefunction as _iw
    import amset.wavefunction.common as _wc

    if STATE["applied"] and _wc.desymmetrize_coefficients is desymmetrize_coefficients:
        return True
    orig = _wc.desymmetrize_coefficients
    try:
        src = inspect.getsource(orig)
    except (OSError, TypeError) as e:
        raise RuntimeError("[desym_fix] 读不到 AMSET 原函数源码（%s），无法核对版本 -> 拒绝打补丁"
                           % e)
    if _squash(ORIG_FACTOR_LINE) not in _squash(src):
        raise RuntimeError(
            "[desym_fix] amset %s 的 desymmetrize_coefficients 里找不到原式\n    %s\n"
            "AMSET 实现已变（或已被修过），本补丁不适用 -> 停跑。请重新核对后再开 %s。"
            % (getattr(amset, "__version__", "?"), ORIG_FACTOR_LINE, ENV_NAME))
    if _iw.desymmetrize_coefficients is not orig:
        raise RuntimeError("[desym_fix] amset.interpolation.wavefunction 里的 "
                           "desymmetrize_coefficients 已被别的东西替换过 -> 拒绝叠加补丁")
    _wc.desymmetrize_coefficients = desymmetrize_coefficients
    _iw.desymmetrize_coefficients = desymmetrize_coefficients
    STATE["applied"] = True
    print("[desym_fix] amset %s：去对称化相位已修正（%s）；SOC 分支另修 SU(2) 与 iσ_y·K（模型检验通过，实数据未验证）"
          % (getattr(amset, "__version__", "?"), FIXED_FACTOR_LINE),
          file=sys.stderr, flush=True)
    return True


apply()
