#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""adp_degen.py —— [0061] 简并价带顶（立方晶体 Γ 点三重简并，无 SOC）的声学形变势（ADP）带间修正。

问题：AMSET 的 ADP 矩阵元是 |ê:D_n(k)|² · |⟨m,k′|n,k⟩|²：只用初态那条带的形变势 D_n（S7 形变构型里按带号排序后的
能移），再乘波函数重叠。对非简并的带这没问题；但 Si/Ge/金刚石/闪锌矿的价带顶是 p 型三重简并（Γ25′/Γ15），剪切应变
（Bir–Pikus 的 b、d）的主要作用是把重、轻空穴**互相耦合**——这是应变算符在简并子空间里的非对角矩阵元，AMSET 的写法
里没有：重、轻空穴在相近 k 上几乎正交，重叠一乘就没了。结果是空穴的 ADP 散射偏弱、迁移率偏高（Si：k·p 模型估计只算
ADP 时高估约 2 倍，见 METHODOLOGY [0061]）。

修正（本模块）：在同一个 k·p + Bir–Pikus 模型里，对每条带 n、每个方向 k̂ 算两种写法的总散射率之比
    r_n(k̂) = Σ_m∫ |⟨m,k′| H_ε(q̂,û_λ) |n,k⟩|² / C_λ   ÷   Σ_m∫ |ê:|D_n(k)||² |⟨m,k′|n,k⟩|² / C_λ      （弹性、均分）
（完整矩阵元 ÷ AMSET 写法；ADP 的 |M|² 与 |k| 无关，比值只依赖方向与带），运行期乘到 AMSET 自己算的 ADP 率上
（amset_dpo3d 的 ADP 类，按初态 (n, k) 乘一个标量；AMSET 用总散射率当 1/τ，所以按初态乘与改矩阵元等价）。
AMSET 的 DFT 形变势、能带、态密度、重叠都不动，只补上它缺的那部分结构。

模型：
  · k·p（Dresselhaus–Kip–Kittel，无 SOC 三带）：H_ii = L k_i² + M(k² − k_i²)，H_ij = N k_i k_j（单位 ħ²/2m0，电子能量）。
    L、M、N 用 AMSET 自己的（插值后）能带在 Γ 附近拟合（fit_kp）。
  · Bir–Pikus（同一个 p 型基）：H_ii = l ε_ii + m(Tr ε − ε_ii)，H_ij = n ε_ij（ε 为张量分量）；a = (l+2m)/3、b = (l−m)/3。
    l、m、|n| 由 S7.1 的形变势在 Γ 点读出（gamma_bir_pikus）：排序后的三条带在 ±δ 应变下的中心差分（带符号文件，
    deformation_signed.h5）给 m = −S_中、l = −2S_外 + S_中；剪切分量的幅值（deformation.h5）给 |n|。
    n 的符号 Γ 点定不出来（劈裂只有 ±|n|），用 Γ 附近非简并 k 点的带符号剪切形变势与 k·p 本征矢的预测比对（sign_of_n）。
  · 弹性：Christoffel 方程（AMSET 同一个 4 阶弹性张量），三支声学模。

约定（与 AMSET / S7.1 一致）：形变势文件存 −ΔE/δ（电子能量；ΔE 含参考能级对齐）；AMSET 原生非对角 = 2·D_ij，
V189 起 S7.1 把幅值文件换成张量分量（attrs dp_offdiag = "tensor"）；带符号文件始终是原生约定。
固定分数坐标的 k 在应变下移动：k′ = (1 − ε)k，所以带符号形变势还含一项 −Σ ε_ab k_b ∂_aE（forward_signed 里一并算）。
"""
import math

import numpy as np

HBAR2_2M = 3.80998212          # eV·Å²，ħ²/2m0
HARTREE_EV = 27.211386245988
DEG_TOL_EV = 0.005             # Γ 点三重简并的判据（eV）
EFIT_EV = 0.15                 # k·p 拟合用的空穴能量上限（eV，低于价带顶）
EMIN_EV = 0.002                # 下限（贴着 Γ 的点能量太小，插值噪声占比大）
FIT_REL_TOL = 0.10             # k·p 拟合的相对残差上限：超过说明不是这个模型（非立方、有 SOC、带太非抛物）
FIT_COND_MAX = 1.0e4           # 拟合雅可比的条件数上限：超过 = 数据定不住 L、M、N（Γ 附近只有高对称方向的点、缺轻空穴带）
EFIT_STEPS = (1.0, 2.0, 4.0)   # 定不住时把拟合窗口放宽到这些倍数再试
FIT_MIN_PER_BAND = 6           # 三条带每条至少这么多个能级：缺一条（通常是轻空穴）时 L、M、N 有一个组合定不住，
#                                而且解落在带序交换的折点上，数值条件数看不出来（合成网格实测）
KCUT_FRAC = 0.25               # 只修正 |k| < KCUT_FRAC × 最短倒格矢 的态（Γ 附近，k·p 适用）
SIGN_KCUT_FRAC = 0.30          # 定 n 的符号时用的 k 点范围（S7 网格粗，放宽一点）
SIGN_MIN_CONF = 0.2            # 两种符号残差的相对差小于它 = 定不出来
NANG = 240                     # 方向网格（Fibonacci 球面）
CARRIERS = ("electron", "hole")


# ---------------------------------------------------------------------------
# k·p 与应变算符
# ---------------------------------------------------------------------------
def kp_matrix(k, L, M, N):
    """DKK 三带 k·p（电子能量，eV）。k (..., 3) 笛卡尔 Å⁻¹。"""
    k = np.asarray(k, float)
    kx, ky, kz = k[..., 0], k[..., 1], k[..., 2]
    k2 = kx * kx + ky * ky + kz * kz
    H = np.empty(k.shape[:-1] + (3, 3))
    for i, ki in enumerate((kx, ky, kz)):
        H[..., i, i] = L * ki * ki + M * (k2 - ki * ki)
    H[..., 0, 1] = H[..., 1, 0] = N * kx * ky
    H[..., 0, 2] = H[..., 2, 0] = N * kx * kz
    H[..., 1, 2] = H[..., 2, 1] = N * ky * kz
    return HBAR2_2M * H


def kp_dmatrix(k, L, M, N):
    """∂H/∂k_a（eV·Å），(..., 3[a], 3, 3)。"""
    k = np.asarray(k, float)
    out = np.zeros(k.shape[:-1] + (3, 3, 3))
    for a in range(3):
        for i in range(3):
            out[..., a, i, i] = 2 * (L if i == a else M) * k[..., a]
        for i in range(3):
            for j in range(3):
                if i != j:
                    if a == i:
                        out[..., a, i, j] += N * k[..., j]
                    elif a == j:
                        out[..., a, i, j] += N * k[..., i]
    return HBAR2_2M * out


def kp_bands(k, L, M, N):
    """-> (能量 (..., 3) 升序，本征矢 (..., 3 分量, 3 带))。"""
    return np.linalg.eigh(kp_matrix(k, L, M, N))


def strain_operator(l, m, n):
    """Bir–Pikus 应变算符 Θ[i, j]（3×3 矩阵，作用在 p 型基上）：H_ε = Σ_ij Θ[i,j] ε_ij（ε 对称、两个非对角都计）。
    -> (3, 3, 3, 3)：前两个指标 = 应变分量，后两个 = 矩阵元。H_xx = l ε_xx + m(ε_yy + ε_zz)，H_xy = n ε_xy。"""
    Th = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for p in range(3):
            Th[i, i, p, p] = l if p == i else m
    for i in range(3):
        for j in range(3):
            if i != j:
                Th[i, j, i, j] = Th[i, j, j, i] = n / 2.0
    return Th


def diag_potentials(c, Th):
    """对角形变势张量 D_ij = c†Θ_ij c（张量约定：ΔE = Σ_ij D_ij ε_ij）。c (..., 3) -> (..., 3, 3)。"""
    return np.einsum("...a,ijab,...b->...ij", c, Th, c)


def forward_signed(kcart, L, M, N, l, m, n):
    """S7.1 带符号形变势（原生约定）的模型值：(nk, 3 带升序, 3, 3)。
    分量 (i,i)：应变 ε_ii = δ；分量 (i,j)：ε_ij = ε_ji = δ（AMSET 的剪切应变，原生 = 2·D_ij）。
    存的是 −ΔE/δ，ΔE = c†Θ(ε)c − Σ_ab ε_ab k_b ∂_aE（固定分数坐标的 k 随应变移动）。"""
    k = np.atleast_2d(np.asarray(kcart, float))
    e, V = kp_bands(k, L, M, N)
    dH = kp_dmatrix(k, L, M, N)
    grad = np.einsum("kin,kaij,kjn->kna", V, dH, V)                         # (nk, 带, a) ∂_aE
    Th = strain_operator(l, m, n)
    out = np.zeros((len(k), 3, 3, 3))
    for i in range(3):
        for j in range(3):
            P = np.zeros((3, 3))
            P[i, j] = P[j, i] = 1.0
            op = np.einsum("ab,abxy->xy", P, Th)
            th = np.einsum("kxn,xy,kyn->kn", V, op, V)
            shift = np.einsum("ab,kb,kna->kn", P, k, grad)
            out[:, :, i, j] = -(th - shift)
    return out


# ---------------------------------------------------------------------------
# 从数据定参数
# ---------------------------------------------------------------------------
def fit_kp(kcart, energy, rank, starts=None):
    """Γ 附近的能带 -> (L, M, N, 信息)。kcart (n,3) Å⁻¹；energy (n,) 电子能量（eV，相对价带顶，≤ 0）；
    rank (n,) 0/1/2 = 三重态里的第几条（能量升序）。最小二乘（多起点，取残差最小）。"""
    k = np.asarray(kcart, float)
    E = np.asarray(energy, float)
    r = np.asarray(rank, int)
    if len(E) < 6:
        raise ValueError("Γ 附近只有 %d 个能级可拟合（至少 6 个）" % len(E))
    k2 = np.sum(k * k, axis=1)
    A0 = float(np.mean(E / (HBAR2_2M * k2)))        # ~ (L+2M)/3 的量级（负）
    s = abs(A0) if A0 else 3.0
    starts = starts or [(-1.3, -0.85, -2.0), (-1.3, -0.85, 2.0), (-0.85, -1.3, -2.0), (-0.85, -1.3, 2.0),
                        (-2.5, -0.5, -2.5), (-0.5, -2.5, -2.5), (-1.0, -1.0, -0.5), (-1.0, -1.0, 0.5)]

    def resid(p):
        ev = np.linalg.eigvalsh(kp_matrix(k, *p))
        return (ev[np.arange(len(r)), r] - E) / np.maximum(np.abs(E), 1e-3)

    best = None
    try:
        from scipy.optimize import least_squares
        for st in starts:
            sol = least_squares(resid, np.asarray(st) * s, method="lm" if len(E) > 3 else "trf")
            if best is None or sol.cost < best[1]:
                best = (sol.x, sol.cost)
    except ImportError:                       # 没有 scipy：Gauss–Newton（数值雅可比）
        for st in starts:
            p = np.asarray(st, float) * s
            for _ in range(200):
                f = resid(p)
                J = np.empty((len(f), 3))
                for a in range(3):
                    dp = np.zeros(3)
                    dp[a] = 1e-6 * max(1.0, abs(p[a]))
                    J[:, a] = (resid(p + dp) - f) / dp[a]
                step = np.linalg.lstsq(J, -f, rcond=None)[0]
                p = p + step
                if np.max(np.abs(step)) < 1e-10:
                    break
            cost = 0.5 * float(np.sum(resid(p) ** 2))
            if best is None or cost < best[1]:
                best = (p, cost)
    p = best[0]
    ev = np.linalg.eigvalsh(kp_matrix(k, *p))[np.arange(len(r)), r]
    rel = float(np.sqrt(np.mean(((ev - E) / np.maximum(np.abs(E), 1e-3)) ** 2)))
    # 可辨识性：数据只沿高对称方向、又缺某条带时，L、M、N 有一个组合定不住（残差照样是 0）
    f0 = resid(p)
    J = np.empty((len(f0), 3))
    for a in range(3):
        dp = np.zeros(3)
        dp[a] = 1e-5 * max(1.0, abs(p[a]))
        J[:, a] = (resid(p + dp) - f0) / dp[a]
    sv = np.linalg.svd(J, compute_uv=False)
    cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf")
    return float(p[0]), float(p[1]), float(p[2]), {"n_points": int(len(E)), "rel_rms": rel, "cond": cond,
                                                    "points_per_band": [int(np.sum(r == i)) for i in range(3)],
                                                    "rms_meV": float(1e3 * np.sqrt(np.mean((ev - E) ** 2)))}


def gamma_bir_pikus(S_gamma, A_gamma, offdiag="tensor"):
    """Γ 点三重态（带按能量升序）的形变势 -> (l, m, |n|, 信息)。
    S_gamma (3,3,3)：带符号中心差分（原生约定）；排序后的三条带在 ε_ii = ±δ 下：中间那条 = −m，两头 = −(l+m)/2。
    A_gamma (3,3,3)：幅值（offdiag = "tensor"：非对角为张量分量，剪切下两头 = |n|/2；"native"：= |n|）。"""
    S = np.asarray(S_gamma, float)
    A = np.asarray(A_gamma, float)
    mid = np.array([S[1, i, i] for i in range(3)])
    outer = np.array([0.5 * (S[0, i, i] + S[2, i, i]) for i in range(3)])
    m = -float(np.mean(mid))
    l = float(-2.0 * np.mean(outer) + np.mean(mid))
    off = [(i, j) for i in range(3) for j in range(3) if i != j]
    n_out = np.array([0.5 * (A[0, i, j] + A[2, i, j]) for i, j in off])
    n_abs = float(np.mean(n_out)) * (2.0 if offdiag == "tensor" else 1.0)
    info = {"diag_spread_eV": float(max(np.ptp(mid), np.ptp(outer))),
            "outer_asym_eV": float(max(abs(S[0, i, i] - S[2, i, i]) for i in range(3))),
            "shear_mid_eV": float(max(abs(A[1, i, j]) for i, j in off)) * (2.0 if offdiag == "tensor" else 1.0),
            "shear_spread_eV": float(np.ptp(n_out)) * (2.0 if offdiag == "tensor" else 1.0)}
    return l, m, n_abs, info


def sign_of_n(kcart, S, usable, L, M, N, l, m, n_abs):
    """Γ 附近非简并 (带, k) 的带符号剪切形变势 vs k·p 预测，定 n 的符号。
    kcart (nk,3)；S (nk, 3 带, 3, 3) 原生约定；usable (nk, 3) 布尔。-> (±1 或 0, 置信度, 用了几个值)。"""
    k = np.atleast_2d(np.asarray(kcart, float))
    use = np.asarray(usable, bool)
    if not use.any():
        return 0, 0.0, 0
    off = [(0, 1), (0, 2), (1, 2)]
    meas = np.asarray(S, float)
    res = {}
    for sg in (+1.0, -1.0):
        pred = forward_signed(k, L, M, N, l, m, sg * n_abs)
        d = [(meas[:, :, i, j] - pred[:, :, i, j])[use] for i, j in off]
        res[sg] = float(np.sum(np.concatenate(d) ** 2))
    best = min(res, key=res.get)
    other = res[-best]
    conf = (other - res[best]) / other if other > 0 else 0.0
    return (int(best) if conf >= SIGN_MIN_CONF else 0), float(conf), int(use.sum() * len(off))


# ---------------------------------------------------------------------------
# 散射率之比
# ---------------------------------------------------------------------------
def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def _cases(L, M, N, dirs):
    """方向 -> (空穴能量系数 h (nd,3)，本征矢 (nd,3,3))；带按电子能量升序（第 0 条最深）。"""
    e, V = kp_bands(dirs, L, M, N)
    h = -e / HBAR2_2M
    if np.any(h <= 0):
        raise ValueError("k·p 参数给出的价带不是向下弯的（有方向 h ≤ 0）：L、M、N 不对")
    return h, V


def rate_ratio(L, M, N, l, m, n, cijkl, nang=NANG, offdiag="tensor", use_abs=True):
    """-> {"dirs": (nang,3), "ratio": (nang,3), "R_full", "R_diag"}：每个初态（方向、带）完整矩阵元与 AMSET 写法的总率之比。
    offdiag：AMSET 运行期用的形变势约定（"tensor" = V189 起的张量分量；"native" = AMSET 原生，非对角 ×2）。
    use_abs=True：AMSET 运行期对形变势逐分量取 |·|（np.abs），这里同样取。"""
    C = np.asarray(cijkl, float)
    dirs = fib_sphere(nang)
    h, V = _cases(L, M, N, dirs)
    Th = strain_operator(l, m, n)
    fac = np.ones((3, 3))
    if offdiag == "native":
        fac += np.ones((3, 3)) - 2 * np.eye(3)            # 非对角 ×2
    w = 4 * np.pi / nang
    kf = dirs[:, None, :] / np.sqrt(h)[:, :, None]       # (nd, 带, 3)
    dos = 1.0 / (2.0 * h ** 1.5)                          # 等能面上的态密度因子（E = 1）
    Rf = np.zeros((nang, 3))
    Rd = np.zeros((nang, 3))
    for i0 in range(nang):
        for b in range(3):
            c0 = V[i0, :, b]
            k0 = kf[i0, b]
            q = kf - k0                                   # (nd, 带, 3)
            qn = np.linalg.norm(q, axis=-1)
            ok = qn > 1e-9
            qh = np.where(ok[..., None], q / np.where(ok, qn, 1.0)[..., None], 0.0)
            G = np.einsum("ijkl,fmj,fml->fmik", C, qh, qh)
            cl, U = np.linalg.eigh(G)                     # (nd, 带, λ)、(nd, 带, i, λ)
            cl = np.where(cl > 0, cl, np.inf)
            Esym = 0.5 * (np.einsum("fmi,fmjl->fmlij", qh, U) + np.einsum("fmj,fmil->fmlij", qh, U))
            Hm = np.einsum("ijab,fmlij->fmlab", Th, Esym)
            amp = np.einsum("fam,fmlab,b->fml", V, Hm, c0)
            D = diag_potentials(c0, Th)
            D = (np.abs(D) if use_abs else D) * fac
            proj = np.einsum("fmi,fmjl,ij->fml", qh, U, D)
            ovl = np.einsum("fam,a->fm", V, c0)
            wt = (dos * ok)[:, :, None] / cl
            Rf[i0, b] = w * float(np.sum(np.abs(amp) ** 2 * wt))
            Rd[i0, b] = w * float(np.sum(proj ** 2 * (ovl ** 2)[:, :, None] * wt))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(Rd > 0, Rf / Rd, 1.0)
    return {"dirs": dirs, "ratio": ratio, "R_full": Rf, "R_diag": Rd}


def model_mobility(L, M, N, l, m, n, cijkl, T=300.0, nang=NANG, which="full"):
    """模型自身的 ADP 限制空穴迁移率（cm²/Vs，非简并统计、弛豫时间 = 散射率）：给测试与文档用（绝对值核对）。
    cijkl 单位 Pa（kT/(ρv²) = kT/C_λ，不需要密度）。"""
    hbar, m0, e, kB = 1.054571817e-34, 9.1093837015e-31, 1.602176634e-19, 1.380649e-23
    c = hbar ** 2 / (2 * m0)
    kT = kB * T
    out = rate_ratio(L, M, N, l, m, n, cijkl, nang=nang)
    R = out["R_full" if which == "full" else "R_diag"]
    dirs = out["dirs"]
    h, V = _cases(L, M, N, dirs)
    # R 的单位：∫dΩ |X|²/C · 1/(2h^1.5)，X 用 eV、C 用 Pa、能量取 E = 1（以 c 为单位）。换成 SI：
    pref = (2 * np.pi / hbar) * kT / (2 * np.pi) ** 3 * e ** 2 / c ** 1.5
    num = den = 0.0
    w = 4 * np.pi / nang
    dH = kp_dmatrix(dirs, L, M, N) / HBAR2_2M            # 以 c 为单位、|k| = 1 时的 ∂H/∂k
    for i in range(nang):
        for b in range(3):
            gx = float(V[i, :, b] @ dH[i, 0] @ V[i, :, b]) * c
            rate = pref * R[i, b]                         # 1/τ = rate·E^{1/2}（E 以 J 计）
            num += w * gx * gx / (2 * (h[i, b] * c) ** 2.5 * rate * hbar ** 2)
            den += w / (2 * (h[i, b] * c) ** 1.5)
    return float(e / kT * kT ** 0.5 / (math.sqrt(math.pi) / 2) * num / den * 1e4)


def lookup(table, kcart, rank):
    """方向表里最近的方向 -> 比值（k = Γ 时取该带各方向的平均）。"""
    k = np.asarray(kcart, float)
    kn = float(np.linalg.norm(k))
    r = np.asarray(table["ratio"])[:, rank]
    if kn < 1e-9:
        return float(np.mean(r))
    return float(r[int(np.argmax(np.asarray(table["dirs"]) @ (k / kn)))])


# ---------------------------------------------------------------------------
# 运行期：AMSET 数据 + S7.1 形变势 -> 修正表
# ---------------------------------------------------------------------------
def _wrap(kf):
    kf = np.asarray(kf, float)
    return kf - np.round(kf)


def _is_cubic(structure):
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        return SpacegroupAnalyzer(structure, symprec=0.01).get_crystal_system() == "cubic"
    except Exception:                                     # noqa: BLE001 —— 判不了就交给 k·p 拟合质量
        return None


def _read_h5(path, prefix, spin_name):
    import h5py
    with h5py.File(str(path), "r") as f:
        key = next((k for k in f.keys() if k.startswith(prefix) and k.endswith(spin_name)), None)
        if key is None:
            raise KeyError("%s 里没有 %s*%s" % (path, prefix, spin_name))
        arr = np.array(f[key], float)
        kp = np.array(f["kpoints"], float)
        attrs = dict(f.attrs)
        gap = np.array(f["gap_" + spin_name], float) if ("gap_" + spin_name) in f else None
    return arr, kp, attrs, gap


def setup(energies, kpoints, vb_idx, structure, cijkl, deform_h5, signed_h5, spin_name="up",
          efit=EFIT_EV, nang=NANG, is_soc=False):
    """-> (修正表 dict 或 None, 报告 dict)。不适用 / 数据不够时返回 None 并在报告 "status" 写原因（调用方按 auto/on 处理）。
    energies (nb, nk) Hartree（AMSET 内部单位）；kpoints (nk,3) 分数坐标；vb_idx 价带顶带号；structure pymatgen；
    cijkl (3,3,3,3) 任意单位（只用比例）；deform_h5 / signed_h5：S7.1 的幅值 / 带符号形变势文件。"""
    rep = {"status": None}
    if is_soc:
        rep["status"] = "SOC 能带：本修正只做无 SOC 的 p 型三重态"
        return None, rep
    E = np.asarray(energies, float) * HARTREE_EV
    kf = _wrap(kpoints)
    vb = int(vb_idx)
    if vb < 2:
        rep["status"] = "价带不足三条"
        return None, rep
    iG = np.where(np.linalg.norm(kf, axis=1) < 1e-6)[0]
    if not len(iG):
        rep["status"] = "AMSET 网格不含 Γ 点"
        return None, rep
    iG = int(iG[0])
    tri = E[vb - 2:vb + 1, iG]
    vbm = float(np.max(E[vb]))
    rep["gamma_triplet_eV"] = [float(x - vbm) for x in tri]
    if np.ptp(tri) > DEG_TOL_EV or abs(tri[-1] - vbm) > 1e-3:
        rep["status"] = "价带顶不是 Γ 点的三重简并（Γ 三条带 %s eV，价带顶 %.4f eV）" % (
            ", ".join("%.4f" % x for x in tri), vbm)
        return None, rep
    if vb >= 3 and E[vb - 3, iG] > tri[0] - DEG_TOL_EV:
        rep["status"] = "Γ 点简并度多于 3（第 4 条价带也在简并组里）"
        return None, rep
    cubic = _is_cubic(structure)
    rep["cubic"] = cubic
    if cubic is False:
        rep["status"] = "不是立方晶系：三带 k·p / Bir–Pikus 的立方形式不适用"
        return None, rep
    rlat = np.asarray(structure.lattice.reciprocal_lattice.matrix, float)
    bmin = float(np.min(np.linalg.norm(rlat, axis=1)))
    kc = kf @ rlat
    kn = np.linalg.norm(kc, axis=1)
    # 1) k·p（定不住就放宽窗口重拟合）
    fit = None
    for mult in EFIT_STEPS:
        ks, es, rs = [], [], []
        for r in range(3):
            hole = vbm - E[vb - 2 + r]
            sel = (hole > EMIN_EV) & (hole < efit * mult) & (kn < KCUT_FRAC * bmin)
            ks.append(kc[sel])
            es.append(E[vb - 2 + r][sel] - vbm)
            rs.append(np.full(int(sel.sum()), r))
        try:
            fit = fit_kp(np.concatenate(ks), np.concatenate(es), np.concatenate(rs)) + (efit * mult,)
        except ValueError as e:
            fit = (None, None, None, {"error": str(e), "cond": float("inf")}, efit * mult)
            continue
        if fit[3]["cond"] <= FIT_COND_MAX and min(fit[3]["points_per_band"]) >= FIT_MIN_PER_BAND:
            break
    L, M, N, finfo, used = fit
    if L is None:
        rep["status"] = "k·p 拟合失败：%s" % finfo["error"]
        return None, rep
    rep["kp"] = {"L": L, "M": M, "N": N, **finfo, "efit_eV": used}
    if finfo["cond"] > FIT_COND_MAX or min(finfo["points_per_band"]) < FIT_MIN_PER_BAND:
        rep["status"] = ("Γ 附近的能级定不住 k·p 参数（条件数 %.1e，各带点数 %s，每条至少 %d；窗口已放宽到 %.2f eV）："
                         "加密 AMSET 插值网格" % (finfo["cond"], finfo["points_per_band"], FIT_MIN_PER_BAND, used))
        return None, rep
    if finfo["rel_rms"] > FIT_REL_TOL:
        rep["status"] = "k·p 拟合残差 %.1f%% > %.0f%%：价带顶不是立方 p 型三重态的形状" % (
            100 * finfo["rel_rms"], 100 * FIT_REL_TOL)
        return None, rep
    # 2) Bir–Pikus
    try:
        A, kpA, attrsA, _g = _read_h5(deform_h5, "deformation_potentials", spin_name)
        S, kpS, attrsS, gap = _read_h5(signed_h5, "signed_", spin_name)
    except (OSError, KeyError) as e:
        rep["status"] = "读不了形变势文件：%s" % e
        return None, rep
    if A.shape[0] != E.shape[0] or S.shape[0] != E.shape[0]:
        rep["status"] = "形变势文件的带数（%d / %d）与 AMSET 的 %d 不同" % (A.shape[0], S.shape[0], E.shape[0])
        return None, rep
    offdiag = "tensor" if str(attrsA.get("dp_offdiag", "")) == "tensor" else "native"
    kS = _wrap(kpS)
    jG = np.where(np.linalg.norm(kS, axis=1) < 1e-6)[0]
    jGA = np.where(np.linalg.norm(_wrap(kpA), axis=1) < 1e-6)[0]
    if not len(jG) or not len(jGA):
        rep["status"] = "S7 的 k 网格不含 Γ 点：读不出 Γ 点的形变势"
        return None, rep
    l, m, n_abs, binfo = gamma_bir_pikus(S[vb - 2:vb + 1, int(jG[0])], A[vb - 2:vb + 1, int(jGA[0])], offdiag)
    rep["bir_pikus_gamma"] = {"l": l, "m": m, "n_abs": n_abs, "a": (l + 2 * m) / 3, "b": (l - m) / 3,
                              "d_abs": n_abs / math.sqrt(3), **binfo, "offdiag": offdiag}
    kSc = kS @ rlat
    near = np.linalg.norm(kSc, axis=1) < SIGN_KCUT_FRAC * bmin
    near &= np.linalg.norm(kSc, axis=1) > 1e-6
    delta = float(attrsS.get("strain_amount") or 0.005)
    St = S[vb - 2:vb + 1].transpose(1, 0, 2, 3)          # (nk, 3, 3, 3)
    if gap is not None:
        g = gap[vb - 2:vb + 1].T                          # (nk, 3)
        thr = 2 * delta * np.max(np.abs(St), axis=(2, 3)) + 1e-3
        usable = near[:, None] & (g > thr)
    else:
        usable = np.repeat(near[:, None], 3, axis=1)
    sg, conf, nval = sign_of_n(kSc, St, usable, L, M, N, l, m, n_abs)
    rep["n_sign"] = {"sign": sg, "confidence": conf, "n_values": nval}
    if sg == 0:
        sg = 1 if N > 0 else -1                           # 定不出来：取与 N 同号（Si/Ge/GaAs 的 d、N 都为负）
        rep["n_sign"]["fallback"] = "与 N 同号"
    n = sg * n_abs
    rep["bir_pikus"] = {"l": l, "m": m, "n": n, "a": (l + 2 * m) / 3, "b": (l - m) / 3, "d": n / math.sqrt(3)}
    # 3) 比值表
    tab = rate_ratio(L, M, N, l, m, n, cijkl, nang=nang, offdiag=offdiag)
    rep["ratio_mean"] = [float(np.mean(tab["ratio"][:, r])) for r in range(3)]
    rep["ratio_range"] = [[float(np.min(tab["ratio"][:, r])), float(np.max(tab["ratio"][:, r]))] for r in range(3)]
    rep["status"] = "ok"
    table = {"bands": [vb - 2, vb - 1, vb], "rlat": rlat, "kcut": KCUT_FRAC * bmin,
             "dirs": tab["dirs"], "ratio": tab["ratio"]}
    return table, rep


def scale(table, band_idx, kfrac):
    """初态 (带, k) 的 ADP 率要乘的系数：不在三重态里或离 Γ 太远 -> 1。"""
    if not table or band_idx not in table["bands"]:
        return 1.0
    k = _wrap(kfrac) @ table["rlat"]
    if np.linalg.norm(k) > table["kcut"]:
        return 1.0
    return lookup(table, k, table["bands"].index(band_idx))
