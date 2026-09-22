#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""za_2d.py —— 2D ZA 弯曲支判定的【单一真源】（三副本合并，2026-09-21）。

背景：同一套「原胞基矢真空轴映射 + 60° 面内基的第二方向 + 本征矢量识别 ZA 支 +
symprec 审计」逻辑此前在三个文件里各写了一份：
  · skill/fc-fit/fc_plot_phonon.py
  · skill/kl-dft-cpu/kl_fc_backends.py
  · skill/_common/mlff/klmlff_common.py   ← 仍是坏的（未迁移，见下）
修一个 bug 要改三处，故收敛到这里。两个已修的文件改为 import 本模块；
mlff 那份由拥有它的会话迁移（TODO）。

关键判据（为什么这些函数长这样）：
  · ZA 弯曲支在无应力 2D 必须 ω ∝ q²（Born-Huang 旋转不变 + Huang 零应力）；
    只看"最小频率没有负值"不够 —— ZA 被线性化（p≈1）时频率全正、照样过虚频闸，κ 却整体错。
  · 生产 primitive_matrix 会把 25 Å 真空映射到【原胞第 0 个基矢】，而 q 方向 / 能带路径
    索引的都是"原胞基矢"下标；把 POSCAR 的 Cartesian 下标直接传进去就是那个方向 bug。
  · 生产原胞的【面内基是 60°】（不是 120° identity），第二面内方向必须取 e_i + s·e_j
    且 s 使它是 Γ-K（否则得到第二个 Γ-M）。

用法（消费者只做别名，不再自带实现）：
    import za_2d as _za
    ZA_P_RANGE = _za.ZA_P_RANGE           # 判据常量也在这里
    za_power_law = _za.za_power_law
    _inplane_qdirs = _za.inplane_qdirs
    ...
"""
import numpy as np

# ---------------------------------------------------------------------------
# 判据常量（唯一真源）
# ---------------------------------------------------------------------------
ZA_P_RANGE = (1.7, 2.3)       # 无应力 2D 的 ZA 弯曲支 p 应落在此区间
ZA_ZFRAC_MIN = 0.5            # 本征矢量沿 Cartesian 真空轴的占比阈值（ZA 支识别）
ZA_MARGIN_MIN = 0.20          # p 距判据边界多近就要人工复核
ZA_R2_MIN = 0.98              # log-log 拟合 R^2 低于此值认为 p 不可信
ZA_SYMPREC_DEFAULT = 1e-5     # phonopy 默认容差
ZA_SYMPREC_RELAX = 1e-4       # DFT 弛豫后数值微畸变原胞的放宽容差


def za_power_law(ph, qdir=(1, 0, 0), qmax=0.05, n=12):
    """沿面内方向取 q ∈ (0, qmax]（倒格子约化单位），拟合最低支 ω ∝ q^p。

    返回 (p, ω_min, fit)，fit = {"r2","log_resid_rms"}；有非正频率时返回
    (None, ω_min, None)。这只是"旧判据"（全局最低支）；生产判据走 za_power_law_eig。
    """
    qs = np.linspace(float(qmax) / int(n), float(qmax), int(n))
    pts = [np.array(qdir, float) * q for q in qs]
    ph.run_qpoints(pts)
    w = np.asarray(ph.get_qpoints_dict()["frequencies"])[:, 0]
    if np.any(w <= 0):
        return None, float(np.min(w)), None
    x = np.log(qs)
    y = np.log(w)
    coeff = np.polyfit(x, y, 1)
    p = float(coeff[0])
    yfit = np.polyval(coeff, x)
    ss_res = float(np.sum((y - yfit) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    fit = {"r2": float(r2), "log_resid_rms": float(np.sqrt(ss_res / max(len(y), 1)))}
    return p, float(np.min(w)), fit


def k_sum_sign(cell, i, j):
    """Sign s such that e_i + s*e_j is the Gamma-K in-plane direction.

    2D 六方倒格子的两个高对称方向是 e_i（Γ-M）与 e_i ± e_j；离 Γ 更远（|倒格矢| 更大）
    的那个才是 Γ-K：两倒格基矢夹角 ≤90°（实空间 120° 原胞）取 +e_j，>90°（生产
    primitive_matrix 给的 60° 面内基）取 -e_j。旧的 e_i+e_j 在 60° 原胞上会得到第二个 Γ-M。
    """
    B = 2.0 * np.pi * np.linalg.inv(np.asarray(cell, float)).T
    return 1.0 if np.dot(B[i], B[j]) >= 0.0 else -1.0


def inplane_qdirs(ph, vac_axis=2):
    """两个"不等价"的面内 q 方向（约化坐标）。

    正交晶格 → (1,0,0)/(0,1,0)；六方等非正交 → Γ-M 与 Γ-K（对称性不同，只查一个会漏
    各向异性）。★ vac_axis 是【原胞基矢】下标（ph.primitive.cell 的行），不是 Cartesian 下标；
    用 vacuum_axis_in_primitive 先映射。
    """
    idx = [i for i in range(3) if i != int(vac_axis)]
    try:
        cell = np.asarray(ph.primitive.cell, float)
        v1, v2 = cell[idx[0]], cell[idx[1]]
        cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    except Exception:
        cell, cos = np.eye(3), 0.0
    d1 = [0.0, 0.0, 0.0]
    d1[idx[0]] = 1.0
    if abs(cos) > 0.2:
        s = k_sum_sign(cell, idx[0], idx[1])
        d2 = list(d1)
        d2[idx[1]] = s
    else:
        d2 = [0.0, 0.0, 0.0]
        d2[idx[1]] = 1.0
    return [tuple(d1), tuple(d2)]


def vacuum_axis_in_primitive(ph, cart_axis):
    """把 Cartesian 真空轴映射成【原胞基矢】下标（面外分量占比最大者）。

    primitive_matrix 可以任意排列原胞基矢；生产 MoS₂ 把真空放在第 0 个基矢。
    cart_axis=None 时退化为"最长基矢即真空"，仍返回合法下标。
    """
    cell = np.asarray(ph.primitive.cell, float)
    norms = np.linalg.norm(cell, axis=1)
    if cart_axis is None:
        return int(np.argmax(norms))
    vhat = np.zeros(3, float)
    vhat[int(cart_axis)] = 1.0
    frac = np.abs(cell @ vhat) / np.maximum(norms, 1e-300)
    return int(np.argmax(frac))


def eig_out_of_plane(ev, natom, cart_vac_axis):
    """每个声子支的本征矢量里，沿 Cartesian 真空轴的权重占比。

    phonopy 的 eigenvectors[q] 形状 (3*natom, 3*natom)：行 = atom*3+cart，列 = 支号。
    """
    ev = np.asarray(ev)
    rows = [a * 3 + int(cart_vac_axis) for a in range(int(natom))]
    power = np.abs(ev) ** 2
    num = np.sum(power[rows, :], axis=0)
    den = np.sum(power, axis=0)
    return np.asarray(num / np.maximum(den, 1e-300), float)


def za_power_law_eig(ph, qdir=(1, 0, 0), qmax=0.05, n=12, cart_vac_axis=2):
    """用本征矢量识别 ZA 支后拟合 ω ~ q^p。

    逐 q 算每支的面外占比，取"面外占比 ≥ ZA_ZFRAC_MIN 的最低频支"（Γ 点光学支也几乎
    100% 面外，全局 argmax 会挑到光学支）。返回 (p, ω_min, fit, extra)；对象无本征矢量
    （假 phonopy）时返回 None，由调用方退回最低支。
    """
    if cart_vac_axis is None:
        return None
    qs = np.linspace(float(qmax) / int(n), float(qmax), int(n))
    pts = [np.array(qdir, float) * q for q in qs]
    try:
        ph.run_qpoints(pts, with_eigenvectors=True)
        d = ph.get_qpoints_dict()
        evs = d.get("eigenvectors")
    except TypeError:
        ph.run_qpoints(pts)
        return None
    if evs is None:
        return None
    w = np.asarray(d["frequencies"], float)
    natom = len(ph.primitive)
    idxs, zsel = [], []
    for i in range(len(qs)):
        z = eig_out_of_plane(evs[i], natom, cart_vac_axis)
        cand = [b for b in range(len(z)) if z[b] >= ZA_ZFRAC_MIN]
        b = int(cand[0]) if cand else int(np.argmax(z))
        idxs.append(b)
        zsel.append(float(z[b]))
    ws = np.asarray([w[i, idxs[i]] for i in range(len(qs))], float)
    branch_note = ""
    if len(set(idxs)) > 1:
        branch_note = ("ZA branch index changes across q (%s): possible branch "
                       "crossing, p may mix two branches" % idxs)
    elif min(zsel) < ZA_ZFRAC_MIN:
        branch_note = ("selected branch is only %.3f out of plane at its worst q "
                       "(< %.2f): in-plane mixing" % (min(zsel), ZA_ZFRAC_MIN))
    p, fit = None, None
    if np.all(ws > 0):
        x = np.log(qs)
        y = np.log(ws)
        c = np.polyfit(x, y, 1)
        p = float(c[0])
        yf = np.polyval(c, x)
        ssr = float(np.sum((y - yf) ** 2))
        sst = float(np.sum((y - float(np.mean(y))) ** 2))
        fit = {"r2": float(1.0 - ssr / sst) if sst > 0 else 1.0,
               "log_resid_rms": float(np.sqrt(ssr / max(len(y), 1)))}
    extra = {"branch_index": idxs, "zfrac_first": zsel[0],
             "zfrac_min": float(min(zsel)), "zfrac_all": zsel, "note": branch_note}
    return p, float(np.min(ws)), fit, extra


def cell_symmetry(ph, symprec):
    """'<intl> (#num) N ops' of the phonopy unit cell at symprec, or None."""
    try:
        from phonopy.structure.symmetry import Symmetry
        ds = Symmetry(ph.unitcell, symprec=symprec).dataset
        return "%s (#%d) %d ops" % (ds.international, ds.number, len(ds.rotations))
    except Exception:
        return None


def relaxed_symprec(ph, cand=ZA_SYMPREC_RELAX):
    """cand when phonopy's default symprec under-detects the cell symmetry, else None.

    DFT 弛豫后的胞只在数值精度内是六方，默认 symprec=1e-5 会丢掉放宽容差能恢复的操作；
    此时 phonopy 在高容差下做的最近像平均才是物理正确的，ZA 拟合用 cand 重做。
    """
    try:
        from phonopy.structure.symmetry import Symmetry
        n0 = len(Symmetry(ph.unitcell, symprec=ZA_SYMPREC_DEFAULT).dataset.rotations)
        n1 = len(Symmetry(ph.unitcell, symprec=cand).dataset.rotations)
        return cand if n1 > n0 else None
    except Exception:
        return None


def clone_phonopy(ph, symprec):
    """Rebuild a Phonopy with the same cell/FC but a different symprec."""
    from phonopy import Phonopy
    ph2 = Phonopy(ph.unitcell, supercell_matrix=ph.supercell_matrix,
                  primitive_matrix=ph.primitive_matrix, symprec=symprec)
    ph2.force_constants = ph.force_constants
    return ph2


__all__ = ["ZA_P_RANGE", "ZA_ZFRAC_MIN", "ZA_MARGIN_MIN", "ZA_R2_MIN",
           "ZA_SYMPREC_DEFAULT", "ZA_SYMPREC_RELAX", "za_power_law", "k_sum_sign",
           "inplane_qdirs", "vacuum_axis_in_primitive", "eig_out_of_plane",
           "za_power_law_eig", "cell_symmetry", "relaxed_symprec", "clone_phonopy"]
