#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""adp_degen.py —— [0061] 简并 / 近简并带边的声学形变势（ADP）带间修正。

通用：任意晶系、任意谷位置（含非 Γ 的谷星）、价带与导带、无 SOC 与 SOC 旋量、三维（S8）与二维（S8.4）。

问题：AMSET 的 ADP 矩阵元是 |ê:|D_n(k)||² · |⟨m,k′|n,k⟩|²：只用初态那条带的对角形变势（S7 形变构型里按带号排序后的
能移，逐分量取 |·|），再乘波函数重叠。带边多重简并（或几条带只差几个 kT）时，应变的主要作用常常是把这几条带**互相耦合**
—— 应变算符在这个子空间里的非对角矩阵元（立方 p 型价带的 Bir–Pikus b、d；六方的 D3–D6；二维 Γ 点 p_x/p_y 双重态的
剪切项……）。AMSET 的写法里没有这部分：初态与末态是子空间里不同的组合，重叠一乘就没了，ADP 偏弱、迁移率偏高
（Si 空穴只算 ADP 时约 1.6–2 倍，见 METHODOLOGY [0061]）。

做法（不预设对称性，全部从本材料的 DFT 数据里读出来）：
  1. 谷：DFT 全网格上带边带的局部极值（离带边 efit 以内），按点群 + 时间反演归成谷星；谷中心 k0 处与带边带逐条相隔
     < WIN 的带构成子空间（d 条；d = 1 不需要修正）。
  2. k·p 矩阵：k0 附近的 DFT 网格点上，用 AMSET 的波函数系数把各带投影到 k0 处的子空间基（Löwdin 正交化），
     H_DFT(δk) = Σ_n (E_n − E_edge)|p_n⟩⟨p_n| —— 与各网格点上波函数的相位、简并子空间里的组合方式都无关。
  3. 对称性：子空间在所有 k 上若可分块（例如二维的 p_z 与 p_x/p_y 不耦合），先按公共不变子空间分块；每块由 H_DFT 解出
     k0 小群里每个操作的表示矩阵 U（H(Mδk) = U H(δk) U†，反幺正操作带复共轭），闭合成群。
  4. k·p 模型 H(δk) = D0 + Σ δk_a B_a + Σ δk_a δk_b Q_ab 限制在对称允许的子空间里：先用网格点上的 H_DFT 线性拟合，
     再用 AMSET 自己的（插值后更密的）能带做本征值最小二乘精修。
  5. 应变算符 H_ε = Σ_ij ε_ij Y_ij（d×d 厄米，同样只取对称允许的部分 —— Bir–Pikus 原理：ε_ij 与 δk_i δk_j 同样变换）：
     (a) k0 点：严格模拟 S7.1 的做法（±δ 应变 -> 按能量排序 -> 能移 / δ），对上带符号中心差分 S 与 AMSET 幅值 A；
     (b) k0 附近的 S7 网格点：非简并组 ⟨c_n|X^(P)|c_n⟩ − Σ P_ab δk_b ∂_aE_n = −S_n,P（固定分数坐标的 k 随应变移动那一项），
         简并组取组内之和（与组内排序、组合无关），离 k0 越远权重越小（k·p 越不准）；
     多起点最小二乘。数据定不住的离散多解（如立方 Bir–Pikus d 的符号）按 (b) 判，判不了取"与 k·p 二次项同号"并写进报告。
  6. 比值表：同一个模型里，对每个能量、方向、带算"完整矩阵元 Σ_λ|⟨m,k′|H_ε(q̂,û_λ)|n,k⟩|²/C_λ"与"AMSET 写法
     Σ_λ(ε_λ:|D_n|)²|⟨m,k′|n,k⟩|²/C_λ"在等能面上的总和之比（弹性、均分；3D 三支声学模，2D 面内两支，与 amset2d_plugin 同口径）。
  7. 运行期：AMSET 的 ADP 率按初态 (带, k) 乘这个比值（amset_dpo3d / amset2d_plugin 的 ADP 类里乘，谷星成员按对称操作
     映回代表谷）。AMSET 用总散射率当 1/τ，按初态乘与改矩阵元等价；DFT 形变势、能带、态密度、重叠都不动。

约定（与 AMSET / S7.1 一致）：形变势存 −ΔE/δ（电子能量，含参考能级对齐）；带符号文件是原生约定（剪切 ε_ij = ε_ji = δ，
存 −⟨X^(ij)⟩，X^(ij) = 2 Y_ij）；AMSET 运行期的幅值文件 V189 起是张量分量（D_ij，ΔE = Σ_ij D_ij ε_ij），运行期逐分量取 |·|。
"""
import itertools
import math

import numpy as np

HARTREE_EV = 27.211386245988
SQ2 = math.sqrt(2.0)
HBAR2_2M_EV = 3.80998212       # ħ²/2m_e（eV·Å²）：m_e v⊗v = ∇E∇E / (2·ħ²/2m_e)
WIN_EV = 0.05                  # 子空间：k0 处与带边带逐条相隔 < WIN（还要大于应变造成的能移，见 _win_eff）
EFIT_EV = 0.15                 # 谷的能量窗口、k·p 拟合与比值表的载流子能量上限（eV）
EFIT_STEPS = (1.0, 2.0, 4.0)   # 点不够时放宽到这些倍数
W_MIN = 0.8                    # 网格点在 k0 子空间里的投影权重下限（低了说明有子空间外的带混进来）
KCUT_FRAC = 0.25               # 只在 |δk| < KCUT_FRAC × 最短（面内）倒格矢 内拟合 / 修正（还不超过到别的谷一半的距离）
TAPER = 0.2                    # 修正在 rcut 前 TAPER·rcut 内线性过渡到 1
FIT_REL_TOL = 0.15             # k·p 拟合的相对残差上限
STRAIN_REL_TOL = 0.2           # 应变算符在 k0 点（严格模拟）的相对残差上限
NB_REL_TOL = 0.35              # 没有 k0 点数据时，邻近 S7 点的相对残差上限
NULL_TOL = 2e-3                # 表示矩阵 / 公共不变子空间：零空间判据（相对最大奇异值）
COMM_TOL = 1e-3
LAMBDA_NODE = 0.1              # 本征值精修时网格点 H_DFT 的锚定权重
AMBIG_TOL = 0.05               # 两个极小的代价相对差小于它 = 数据分不开
AMBIG_SPREAD = 0.10            # 分不开的几个解给出的平均放大系数相差不超过它才取平均，否则不修正
N_START = 24                   # 应变拟合的随机起点数
NANG = {3: 240, 2: 96}         # 方向网格
NE = 8                         # 能量网格（对数均匀）
E_MIN_EV = 0.005
NT = 240                       # 射线上的采样点（找等能点）
KEY_RES = 1e-4                 # k 点匹配的分辨率（分数坐标）


# ---------------------------------------------------------------------------
# 厄米矩阵基、单项式与应变坐标
# ---------------------------------------------------------------------------
def herm_basis(d):
    """d×d 厄米矩阵的正交归一基（Tr(B_p B_q) = δ_pq），d² 个。"""
    out = []
    for i in range(d):
        M = np.zeros((d, d), complex)
        M[i, i] = 1.0
        out.append(M)
    for i in range(d):
        for j in range(i + 1, d):
            M = np.zeros((d, d), complex)
            M[i, j] = M[j, i] = 1 / SQ2
            out.append(M)
            M = np.zeros((d, d), complex)
            M[i, j], M[j, i] = -1j / SQ2, 1j / SQ2
            out.append(M)
    return np.array(out)


def herm_coords(H, basis):
    """厄米矩阵 (..., d, d) -> 实坐标 (..., d²)：x_p = Tr(B_p H)。"""
    return np.real(np.einsum("pji,...ij->...p", basis, H))


def pairs(dim):
    """对称指标对 (a ≤ b)：3D 六个，2D 面内三个。"""
    return [(a, b) for a in range(dim) for b in range(a, dim)]


def n_terms(dim):
    return dim + len(pairs(dim))


def monomials(k, dim):
    """δk (..., 3) -> 模型单项式 (..., nt)：δk_a（a < dim），再是 δk δkᵀ 的正交归一坐标（δk_a²；a < b 时 √2 δk_a δk_b）。"""
    k = np.asarray(k, float)
    out = [k[..., a] for a in range(dim)]
    out += [k[..., a] ** 2 if a == b else SQ2 * k[..., a] * k[..., b] for a, b in pairs(dim)]
    return np.stack(out, -1)


def dmonomials(k, dim):
    """∂(单项式)/∂δk_a，(..., 3, nt)。"""
    k = np.asarray(k, float)
    out = np.zeros(k.shape[:-1] + (3, n_terms(dim)))
    for a in range(dim):
        out[..., a, a] = 1.0
    for t, (a, b) in enumerate(pairs(dim)):
        c = dim + t
        if a == b:
            out[..., a, c] = 2 * k[..., a]
        else:
            out[..., a, c] = SQ2 * k[..., b]
            out[..., b, c] = SQ2 * k[..., a]
    return out


def strain_coords(eps, dim):
    """对称应变 (..., 3, 3) -> 正交归一坐标 (..., np)（与 δk δkᵀ 同一套：ε_aa；a < b 时 √2 ε_ab）。"""
    eps = np.asarray(eps)
    return np.stack([eps[..., a, a] if a == b else SQ2 * eps[..., a, b] for a, b in pairs(dim)], -1)


# ---------------------------------------------------------------------------
# k·p 模型与应变算符
# ---------------------------------------------------------------------------
def model_H(model, dk):
    """模型 {"dim", "D0" (d,d), "C" (nt,d,d)}：H(δk) = D0 + Σ_t f_t(δk) C_t（eV，相对带边）。δk 笛卡尔 (..., 3)。"""
    f = monomials(dk, model["dim"])
    return model["D0"] + np.einsum("...t,tij->...ij", f, model["C"])


def model_dH(model, dk):
    """∂H/∂δk_a，(..., 3, d, d)。"""
    return np.einsum("...at,tij->...aij", dmonomials(dk, model["dim"]), model["C"])


def strain_H(Z, eps, dim):
    """H_ε = Σ_t e_t(ε) Z_t = Σ_ij ε_ij Y_ij。Z (np, d, d)。"""
    return np.einsum("...t,tij->...ij", strain_coords(eps, dim), Z)


def pattern_op(Z, dim, P):
    """应变图样 P = (a, b) 的算符 X^(P)（原生约定：剪切 ε_ab = ε_ba = 1）。"""
    t = pairs(dim).index(tuple(P))
    return Z[..., t, :, :] * (1.0 if P[0] == P[1] else SQ2)


def tensor_Y(Z, dim):
    """Z -> Y (3,3,d,d)：Y_aa = Z_aa，Y_ab = Y_ba = Z_ab/√2（张量约定，H_ε = Σ_ij ε_ij Y_ij）。"""
    d = Z.shape[-1]
    Y = np.zeros((3, 3, d, d), complex)
    for t, (a, b) in enumerate(pairs(dim)):
        if a == b:
            Y[a, a] = Z[t]
        else:
            Y[a, b] = Y[b, a] = Z[t] / SQ2
    return Y


def Z_from_Y(Y, dim):
    return np.array([Y[a, a] if a == b else SQ2 * Y[a, b] for a, b in pairs(dim)])


# ---------------------------------------------------------------------------
# 对称性：表示矩阵、闭合成群、对称允许的子空间、公共不变子空间
# ---------------------------------------------------------------------------
def _ad(U, C, anti):
    """群元在系数矩阵上的作用：U† C U（反幺正再取复共轭）。C (..., d, d)。"""
    out = np.einsum("ji,...jk,kl->...il", U.conj(), C, U)
    return out.conj() if anti else out


def term_T(M, dim):
    """f(Mδk) = T f(δk) 的 T（nt × nt）。前 dim 个是线性项，后面是二次项（也是应变坐标的变换）。"""
    rng = np.random.default_rng(1)
    nt = n_terms(dim)
    K = rng.normal(size=(nt + 8, 3))
    if dim == 2:
        K[:, 2] = 0.0
    F = monomials(K, dim)
    F2 = monomials(K @ np.asarray(M, float).T, dim)
    return np.linalg.lstsq(F, F2, rcond=None)[0].T


def layer_op(M):
    """二维：只保留不把面内与 z 混合的操作。"""
    M = np.asarray(M, float)
    return bool(np.all(np.abs(M[2, :2]) < 1e-6) and np.all(np.abs(M[:2, 2]) < 1e-6))


def solve_rep(pairs_H, d, anti):
    """H(Mδk) = U H(δk)[*] U† 的 U。pairs_H = [(H(δk), H(Mδk)), ...]。零空间一维（U 唯一到相位）-> (U, 残差)，否则 (None, 原因)。"""
    I = np.eye(d)
    A = np.concatenate([np.kron(H2, I) - np.kron(I, (H.conj() if anti else H).T) for H, H2 in pairs_H])
    _u, s, vh = np.linalg.svd(A, full_matrices=False)
    smax = s[0] if s[0] > 0 else 1.0
    if s[-1] > NULL_TOL * smax:
        return None, "无解"
    if d > 1 and s[-2] < 5 * NULL_TOL * smax:
        return None, "不唯一"
    D = vh[-1].conj().reshape(d, d)
    W, sv, Vh = np.linalg.svd(D)
    if sv[-1] < 0.9 * sv[0]:
        return None, "非幺正"
    U = W @ Vh
    err = max(np.linalg.norm(H2 - U @ (H.conj() if anti else H) @ U.conj().T) / (np.linalg.norm(H2) + 1e-12)
              for H, H2 in pairs_H)
    if err > 5e-2:
        return None, "残差 %.2g" % err
    return U, float(err)


def _mkey(M, anti):
    return (tuple((np.round(np.asarray(M, float), 4) + 0.0).ravel().tolist()), bool(anti))


def close_group(gens, d, max_size=512):
    """生成元 [(M, U, anti)] -> 闭合的群（含单位元）；同一个 (M, anti) 得到不同的 U（相位之外）-> None。
    复合 g∘h：M = M_g M_h，U = U_g · (U_h 的复共轭，若 g 反幺正)，anti = anti_g ⊕ anti_h。"""
    e = (np.eye(3), np.eye(d, dtype=complex), False)
    G = {_mkey(e[0], False): e}
    frontier = [e]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                M = g[0] @ h[0]
                U = g[1] @ (h[1].conj() if g[2] else h[1])
                a = bool(g[2]) ^ bool(h[2])
                key = _mkey(M, a)
                if key in G:
                    if abs(abs(np.trace(G[key][1].conj().T @ U)) / d - 1.0) > 1e-2:
                        return None
                    continue
                G[key] = (M, U, a)
                new.append(G[key])
                if len(G) > max_size:
                    return None
        frontier = new
    return list(G.values())


def equivariant_basis(group, d, Ts):
    """对称允许的系数子空间。参数 = len(T) 个 d×d 厄米矩阵的实坐标；Ts[g] 是该群元的项变换矩阵。
    对称化 C′_s = (1/|G|) Σ_g Σ_t T_g[t,s] Ad_g(C_t) 是正交投影；返回它的值域的正交归一基 (npar, neq)。"""
    basis = herm_basis(d)
    nt = Ts[0].shape[0]
    npar = nt * d * d
    Pm = np.zeros((npar, npar))
    for j in range(npar):
        x = np.zeros(npar)
        x[j] = 1.0
        C = np.einsum("tp,pij->tij", x.reshape(nt, d * d), basis)
        acc = np.zeros_like(C)
        for (M, U, anti), T in zip(group, Ts):
            acc += np.einsum("ts,tij->sij", T, _ad(U, C, anti))
        Pm[:, j] = herm_coords(acc / len(group), basis).ravel()
    w, V = np.linalg.eigh(0.5 * (Pm + Pm.T))
    return V[:, w > 0.5]


def commutant_blocks(Hs, d):
    """一组厄米矩阵的公共不变子空间 -> (W 幺正（列 = 新基）, [块的下标列表])。只有单位矩阵与之对易 -> 一块。"""
    if d == 1:
        return np.eye(1, dtype=complex), [[0]]
    I = np.eye(d)
    tl = sum(np.linalg.norm(H - np.trace(H) / d * I) ** 2 for H in Hs)
    if tl <= 1e-18 * (sum(np.linalg.norm(H) ** 2 for H in Hs) + 1e-30):
        return np.eye(d, dtype=complex), [[i] for i in range(d)]           # 处处 ∝ 单位阵（Kramers 对）：任意分解都行
    A = np.concatenate([np.kron(H, I) - np.kron(I, H.T) for H in Hs])
    _u, s, vh = np.linalg.svd(A, full_matrices=False)
    smax = s[0] if s[0] > 0 else 1.0
    null = [v.conj().reshape(d, d) for v, x in zip(vh, s) if x < COMM_TOL * smax]
    if len(null) <= 1:
        return np.eye(d, dtype=complex), [list(range(d))]
    rng = np.random.default_rng(7)
    Ch = sum(rng.normal() * (N + N.conj().T) + rng.normal() * 1j * (N - N.conj().T) for N in null)
    w, W = np.linalg.eigh(Ch)
    tol = 1e-3 * (np.ptp(w) + 1e-12)
    blocks, cur = [], [0]
    for i, g in enumerate(np.diff(w)):
        if g > tol:
            blocks.append(cur)
            cur = [i + 1]
        else:
            cur.append(i + 1)
    blocks.append(cur)
    mask = _block_mask(blocks, d)
    off = tot = 0.0
    for H in Hs:
        Hr = W.conj().T @ H @ W
        off += np.linalg.norm(Hr[~mask]) ** 2
        tot += np.linalg.norm(Hr) ** 2
    if off > 1e-4 * (tot + 1e-30):
        return np.eye(d, dtype=complex), [list(range(d))]
    return W, blocks


def _block_mask(blocks, d):
    m = np.zeros((d, d), bool)
    for b in blocks:
        m[np.ix_(b, b)] = True
    return m


def lsq(fun, x0, jac=None, max_iter=200, tol=1e-12):
    """Levenberg–Marquardt（纯 numpy，不依赖 scipy）。-> (x, 残差, 雅可比)。jac 缺省用前向差分。"""
    def J_of(x, r):
        if jac is not None:
            return np.asarray(jac(x), float)
        J = np.empty((len(r), len(x)))
        for i in range(len(x)):
            h = 1e-7 * max(1.0, abs(x[i]))
            xp = x.copy()
            xp[i] += h
            J[:, i] = (fun(xp) - r) / h
        return J
    x = np.array(x0, float)
    r = fun(x)
    cost = float(r @ r)
    J = J_of(x, r)
    lam = None
    for _ in range(max_iter):
        A = J.T @ J
        g = J.T @ r
        dA = np.diag(A).copy()
        dA[dA <= 0] = 1.0
        if lam is None:
            lam = 1e-3 * float(np.max(dA))
        improved = False
        for _try in range(30):
            try:
                step = -np.linalg.solve(A + lam * np.diag(dA), g)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            xn = x + step
            rn = fun(xn)
            cn = float(rn @ rn)
            if cn < cost:
                improved = True
                break
            lam *= 4
        if not improved:
            break
        rel = (cost - cn) / max(cost, 1e-300)
        x, r, cost = xn, rn, cn
        J = J_of(x, r)
        lam = max(lam / 3, 1e-15)
        if rel < tol or np.linalg.norm(step) < 1e-12 * (np.linalg.norm(x) + 1e-12):
            break
    return x, r, J


# ---------------------------------------------------------------------------
# 声学模、方向、等能面上的比值表
# ---------------------------------------------------------------------------
def directions(n, dim):
    """3D：Fibonacci 球面；2D：xy 平面上的圆。(n, 3)。"""
    if dim == 2:
        th = 2 * np.pi * (np.arange(n) + 0.5) / n
        return np.stack([np.cos(th), np.sin(th), np.zeros(n)], 1)
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def acoustic_modes(C, qh, dim):
    """Christoffel -> (刚度 (n, nb), 极化 (n, 3, nb))。3D 三支；2D 只取面内 2×2 子块（面内传播、面内极化的 LA/TA，
    与 amset2d_plugin._inplane_modes 同口径）。C (3,3,3,3) 任意单位（只用比例）。"""
    qh = np.asarray(qh, float)
    G = np.einsum("ijkl,...j,...l->...ik", np.asarray(C, float), qh, qh)
    if dim == 2:
        w, U2 = np.linalg.eigh(G[..., :2, :2])
        U = np.zeros(U2.shape[:-2] + (3, 2))
        U[..., :2, :] = U2
        return w, U
    return np.linalg.eigh(G)


def _crossings(eps, ts, E):
    """eps (nd, nt, d)（载流子能量，按载流子序）沿射线 -> 第一次从 < E 穿到 ≥ E 的 t（线性插值）与是否存在。"""
    f = eps - E
    cr = (f[:, :-1, :] < 0) & (f[:, 1:, :] >= 0)
    has = cr.any(axis=1)
    j = np.argmax(cr, axis=1)
    f0 = np.take_along_axis(f, j[:, None, :], 1)[:, 0]
    f1 = np.take_along_axis(f, (j + 1)[:, None, :], 1)[:, 0]
    t = ts[j] + (ts[j + 1] - ts[j]) * (-f0) / np.where(has, f1 - f0, 1.0)
    return np.where(has, t, np.nan), has


def _ray_states(model, dirs, t, ok, E, sign, rmax):
    """牛顿法把等能点精化到 ε = E；-> (k (nd,d,3), 本征矢 (nd,d,d), 等能面权重 (nd,d), ok)。"""
    nd, d = t.shape
    dim = model["dim"]
    ar = np.arange(d)
    band_e = ar if sign > 0 else ar[::-1]                     # 载流子序 n -> 电子能量序
    tt = np.where(ok, t, 0.0)
    for it in range(4):
        if it == 3:
            # 简并带（Kramers 对、高对称方向）落在同一点：必须用同一次对角化的本征矢，否则简并子空间里的基各取各的、互不正交
            for n in range(1, d):
                for m in range(n):
                    same = np.abs(tt[:, n] - tt[:, m]) < 1e-9 * rmax
                    tt[:, n] = np.where(same, tt[:, m], tt[:, n])
        k = dirs[:, None, :] * tt[..., None]
        e, V = np.linalg.eigh(model_H(model, k))
        en = sign * e[:, ar, band_e]
        v = np.swapaxes(V, 2, 3)[:, ar, band_e, :]            # (nd, d, comp)
        proj = np.einsum("ia,inajk->injk", dirs, model_dH(model, k))
        g = sign * np.real(np.einsum("inj,injk,ink->in", v.conj(), proj, v))
        if it == 3:
            break
        step = np.where(ok & (np.abs(g) > 1e-12), (E - en) / np.where(np.abs(g) > 1e-12, g, 1.0), 0.0)
        tt = np.clip(tt + step, 0.0, rmax)
    ok = ok & (g > 1e-9) & (np.abs(en - E) < 1e-3 * max(E, 1e-3)) & (tt > 0)
    w = np.where(ok, tt ** (dim - 1) / np.where(ok, g, 1.0), 0.0)
    return k, v, w, ok


def rate_ratio_table(model, Y, C, sign, energies, nang=None, offdiag="tensor", amset_terms=True):
    """-> {"dirs", "E", "ratio" (nE, nd, d)}：初态（能量、方向、带）完整矩阵元的率 ÷ AMSET 写法的率。
    带按载流子能量升序（空穴：电子能量降序）。没有这个初态（等能面不到这个方向 / 这条带）记 1。
    Y (3,3,d,d) 张量约定；offdiag = "native" 时 AMSET 写法的非对角 |D| 取两倍（旧的 S7.1 h5）。
    amset_terms（出厂开）：两边都照 AMSET 弹性散射的实际算法 ——
      · 末态乘动量弛豫因子 1 − v·v′/|v|²（AMSET calculate_rate 对弹性机制总是乘）；
      · AMSET 写法的形变势是 S7 文件里固定分数坐标的值再取 |·|、加 m_e v⊗v：|Y_n − sym(δk⊗∇E_n)| + m_e v⊗v；
        完整矩阵元 = ⟨m,k′|H_ε|n,k⟩ + ε:(m_e v⊗v − sym(δk⊗∇E))_n ⟨m,k′|n,k⟩（运动学项与 v⊗v 按初态的带内项，
        与 AMSET 同口径；带间部分只来自应变算符）。
    这样比值只差在"带间应变耦合"这一处，乘到 AMSET 自己的率上不引入别的口径差。amset_terms=False：两边都只用 H_ε
    与总散射率（[0061] 首版的口径，测试用）。"""
    dim = model["dim"]
    d = np.shape(model["D0"])[0]
    nang = int(nang or NANG[dim])
    dirs = directions(nang, dim)
    rmax = float(model["rmax"])
    ts = np.linspace(0.0, rmax, NT + 1)
    eps = sign * np.linalg.eigvalsh(model_H(model, dirs[:, None, :] * ts[None, :, None]))
    if sign < 0:
        eps = eps[..., ::-1]
    fac = np.ones((3, 3)) if offdiag == "tensor" else (2.0 - np.eye(3))
    Yf = np.asarray(Y).reshape(9, d, d)
    ratio = np.ones((len(energies), nang, d))
    for ie, E in enumerate(energies):
        t, ok = _crossings(eps, ts, float(E))
        k, c, w, ok = _ray_states(model, dirs, t, ok, float(E), sign, rmax)
        if not ok.any():
            continue
        idx = np.argwhere(ok)
        kf, cf, wf = k[ok], c[ok], w[ok]
        Zf = np.einsum("fa,pab->fpb", cf.conj(), Yf)                       # c_f† Y_p
        Dall = np.real(np.einsum("fpb,fb->fp", Zf, cf)).reshape(-1, 3, 3)  # 各态的对角形变势张量（本征）
        if amset_terms:
            g = np.real(np.einsum("fa,fxab,fb->fx", cf.conj(), model_dH(model, kf), cf))   # ∇E（eV·Å）
            Tk = 0.5 * (np.einsum("fb,fa->fab", kf, g) + np.einsum("fa,fb->fab", kf, g))   # 固定分数坐标的运动学项
            VV = np.einsum("fa,fb->fab", g, g) / (2.0 * HBAR2_2M_EV)       # m_e v⊗v = ∇E∇E/(ħ²/m_e)
            g2 = np.einsum("fa,fa->f", g, g)
        for n0, (i0, b0) in enumerate(idx):
            q = kf - kf[n0]
            qn = np.linalg.norm(q, axis=1)
            sel = qn > 1e-9 * rmax
            qh = np.where(sel[:, None], q / np.where(sel, qn, 1.0)[:, None], np.array([1.0, 0.0, 0.0]))
            if dim == 2:
                qh[:, 2] = 0.0
            Cl, U = acoustic_modes(C, qh, dim)
            Es = 0.5 * (np.einsum("fi,fjl->flij", qh, U) + np.einsum("fj,fil->flij", qh, U))
            nb = U.shape[-1]
            ov = cf.conj() @ cf[n0]
            amp = np.einsum("flp,fp->fl", Es.reshape(len(qh), nb, 9), Zf @ cf[n0])
            wt = (wf * sel)[:, None] / np.where(Cl > 0, Cl, np.inf)
            if amset_terms:
                amp = amp + np.einsum("flij,ij->fl", Es, VV[n0] - Tk[n0]) * ov[:, None]
                Dam = np.abs(Dall[n0] - Tk[n0]) * fac + VV[n0]
                if g2[n0] > 0:
                    wt = wt * (1.0 - (g @ g[n0]) / g2[n0])[:, None]
            else:
                Dam = np.abs(Dall[n0]) * fac
            proj = np.einsum("flij,ij->fl", Es, Dam)
            num = float(np.sum(np.abs(amp) ** 2 * wt))
            den = float(np.sum(proj ** 2 * (np.abs(ov) ** 2)[:, None] * wt))
            ratio[ie, i0, b0] = num / den if den > 0 else 1.0
    return {"dirs": dirs, "E": np.asarray(energies, float), "ratio": ratio}


# ---------------------------------------------------------------------------
# 网格、谷与谷星
# ---------------------------------------------------------------------------
def _wrap(kf):
    kf = np.asarray(kf, float)
    return kf - np.round(kf)


def _key(kf):
    n = int(round(1 / KEY_RES))
    a = np.mod(np.round(np.asarray(kf, float), 6), 1.0)
    return tuple((np.round(a / KEY_RES).astype(int) % n).tolist())


def local_extrema(Eb, kfrac, sign, axes, reach=2):
    """带 Eb (nk,) 在完整规则网格上的局部极值（sign < 0 极大，> 0 极小）下标；与 axes 方向 ±reach 步以内的点比
    （只比最近邻时，很平的方向上会出现网格假极值：重空穴沿 [110] 走一步能量不变、绕到 [111] 才降）。"""
    kf = np.mod(np.round(np.asarray(kfrac, float), 6), 1.0)
    n = np.ones(3, int)
    addr = np.zeros((len(kf), 3), int)
    for a in range(3):
        u = np.unique(kf[:, a])
        n[a] = len(u)
        if n[a] > 1:
            s = (kf[0, a] * n[a]) % 1.0
            addr[:, a] = np.round(kf[:, a] * n[a] - s).astype(int) % n[a]
    if np.prod(n) != len(kf) or len({tuple(x) for x in addr.tolist()}) != len(kf):
        raise ValueError("DFT 的 k 点不是完整的规则网格")
    grid = np.full(tuple(n), np.nan)
    grid[tuple(addr.T)] = Eb
    ext = np.ones(tuple(n), bool)
    for off in itertools.product(range(-reach, reach + 1), repeat=len(axes)):
        if not any(off):
            continue
        sh = [0, 0, 0]
        for a, o in zip(axes, off):
            sh[a] = o
        ext &= sign * grid <= sign * np.roll(grid, shift=sh, axis=(0, 1, 2)) + 1e-9
    return np.where(ext[tuple(addr.T)])[0]


def group_stars(idx, kfrac, rlat, ops, Eb=None, tol=1e-3):
    """极值点下标 -> [(代表, [(成员, M)])]：M k_代表 ≡ k_成员（笛卡尔，模倒格矢），ops 已含时间反演（−R）。
    给了带能量 Eb 时成员还要与代表同能量（差 < tol eV；磁性体系里 −k 不一定等价）。"""
    rinv = np.linalg.inv(rlat)
    left = list(idx)
    stars = []
    while left:
        r = left.pop(0)
        kc = np.asarray(kfrac[r], float) @ rlat
        members, rest = [(r, np.eye(3))], []
        for m in left:
            hit = None
            if Eb is not None and abs(float(Eb[m]) - float(Eb[r])) > tol:
                rest.append(m)
                continue
            for M in ops:
                f = (M @ kc) @ rinv - kfrac[m]
                if np.max(np.abs(f - np.round(f))) < 1e-4:
                    hit = np.asarray(M, float)
                    break
            if hit is None:
                rest.append(m)
            else:
                members.append((m, hit))
        left = rest
        stars.append((r, members))
    return stars


def manifold(E, k0i, edge_band, sign, win):
    """k0 处从带边带往里逐条相隔 < win 的带（升序带号）。"""
    bands, b, step = [edge_band], edge_band, (-1 if sign < 0 else 1)
    while 0 <= b + step < E.shape[0] and abs(E[b + step, k0i] - E[b, k0i]) < win:
        b += step
        bands.append(b)
    return sorted(bands)


def _dk_cart(kf, k0f, rlat, dim, vac_axis):
    """-> (δk = k − k0 笛卡尔（2D 去掉 z）, 是否与 k0 同一个面（2D 真空方向分数坐标相同）)。"""
    df = _wrap(np.asarray(kf, float) - np.asarray(k0f, float))
    dk = df @ rlat
    if dim == 2:
        dk = dk.copy()
        dk[..., 2] = 0.0
        return dk, np.abs(df[..., vac_axis]) < 1e-6
    return dk, np.ones(dk.shape[:-1], bool)


def subspace_H(C0, Ck, Ek):
    """k0 处子空间基系数 C0 (d, …)、网格点上同一组带的系数 Ck (d, …)、能量 Ek (d,)（相对带边）->
    (H_DFT (d,d), 投影权重最小值)。投影 P_αn = ⟨α|n⟩，Löwdin 正交化后 H = P diag(E) P†。"""
    C0, Ck = np.asarray(C0), np.asarray(Ck)
    P = C0.reshape(len(C0), -1).conj() @ Ck.reshape(len(Ck), -1).T
    w, U = np.linalg.eigh(P.conj().T @ P)
    if np.min(w) <= 1e-8:
        return None, float(np.min(w))
    Po = P @ (U @ np.diag(w ** -0.5) @ U.conj().T)
    return Po @ np.diag(np.asarray(Ek, float)) @ Po.conj().T, float(np.min(w))


def align_coefficients(c, g, g_maps):
    """与 AMSET get_overlap 同一处理：网格点与参考点隔着倒格矢 g（整数 3 元组）时，把系数按 g_maps 平移到参考点的 G 基上。
    c (ncoeff[, 2])；g_maps (3,3,3,ncoeff) 的 -1 指向末尾补的零。"""
    g = np.asarray(g, int)
    if not np.any(g):
        return c
    if np.any(np.abs(g) > 1):
        raise ValueError("k 与参考点相差的倒格矢超过一个（%s）" % g)
    pad = np.concatenate([c, np.zeros((1,) + c.shape[1:], c.dtype)])
    return pad[np.asarray(g_maps)[1 - g[0], 1 - g[1], 1 - g[2], :]]


# ---------------------------------------------------------------------------
# 一个谷：投影、分块、对称性、k·p、应变算符、比值表
# ---------------------------------------------------------------------------
def _node_data(k0i, bands, E, kfrac, rlat, coeff_fn, rcut, dim, vac_axis, E_edge, chunk=64):
    k0f = np.asarray(kfrac[k0i], float)
    dk, plane = _dk_cart(kfrac, k0f, rlat, dim, vac_axis)
    kn = np.linalg.norm(dk, axis=1)
    idx = np.where(plane & (kn < rcut))[0]
    d = len(bands)
    C0 = np.asarray(coeff_fn(list(bands), np.repeat(k0f[None], d, 0), k0f))
    keep, Hs, ws = [], [], []
    for s in range(0, len(idx), chunk):
        J = idx[s:s + chunk]
        Ck = np.asarray(coeff_fn(np.tile(bands, len(J)), np.repeat(np.asarray(kfrac)[J], d, 0), k0f))
        Ck = Ck.reshape((len(J), d) + Ck.shape[1:])
        for j, i in enumerate(J):
            H, w = subspace_H(C0, Ck[j], E[bands, i] - E_edge)
            if H is not None and w >= W_MIN:
                keep.append(i)
                Hs.append(0.5 * (H + H.conj().T))
                ws.append(w)
    keep = np.array(keep, int)
    return {"idx": keep, "dk": dk[keep], "H": np.array(Hs).reshape(-1, d, d), "w": np.array(ws),
            "keys": [_key(kfrac[i]) for i in keep]}


def _little_group(k0f, rlat, rots, dim):
    rinv = np.linalg.inv(rlat)
    k0c = np.asarray(k0f, float) @ rlat
    out, seen = [], set()
    for R in rots:
        for M, anti in ((np.asarray(R, float), False), (-np.asarray(R, float), True)):
            if dim == 2 and not layer_op(M):
                continue
            f = (M @ k0c) @ rinv - k0f
            if np.max(np.abs(f - np.round(f))) < 1e-4 and _mkey(M, anti) not in seen:
                seen.add(_mkey(M, anti))
                out.append((M, anti))
    return out


def _block_group(nodes, Hb, k0f, rlat, ops, d):
    """一块：每个小群操作由网格点数据解出表示矩阵，闭合成群。-> (群, 报告)。"""
    rinv = np.linalg.inv(rlat)
    lut = {k: j for j, k in enumerate(nodes["keys"])}
    gens, skipped = [], 0
    for M, anti in ops:
        if np.allclose(M, np.eye(3)) and not anti:
            continue
        prs = []
        for j, dk in enumerate(nodes["dk"]):
            j2 = lut.get(_key(k0f + (M @ dk) @ rinv))
            if j2 is not None:
                prs.append((Hb[j], Hb[j2]))
        if len(prs) < 2:
            skipped += 1
            continue
        U, _info = solve_rep(prs, d, anti)
        if U is None:
            skipped += 1
            continue
        gens.append((M, U, anti))
    G = close_group(gens, d) if gens else None
    if G is None:
        return [(np.eye(3), np.eye(d, dtype=complex), False)], {"n_ops": 1, "skipped": skipped,
                                                                "closure": "failed" if gens else "none"}
    return G, {"n_ops": len(G), "skipped": skipped}


def _embed(blocks, d, mats):
    """[(nt, dβ, dβ)] -> (nt, d, d) 块对角。"""
    nt = mats[0].shape[0]
    out = np.zeros((nt, d, d), complex)
    for b, m in zip(blocks, mats):
        out[:, np.array(b)[:, None], np.array(b)[None, :]] = m
    return out


def _param_mats(blocks, d, bases, nt):
    """每块的对称允许基 -> 全局参数矩阵 (nq, nt, d, d)（块对角）。"""
    out = []
    for bi, (b, V) in enumerate(zip(blocks, bases)):
        db = len(b)
        hb = herm_basis(db)
        for q in range(V.shape[1]):
            m = np.einsum("tp,pij->tij", V[:, q].reshape(nt, db * db), hb)
            mats = [np.zeros((nt, len(bb), len(bb)), complex) for bb in blocks]
            mats[bi] = m
            out.append(_embed(blocks, d, mats))
    return np.array(out).reshape(-1, nt, d, d)


def _fit_model(nodes, Hr, D0r, Cq, dim, sign, efit, dense, rcut):
    """k·p 模型（对称允许的参数 γ）。网格点按载流子能量加权（w = 1/(1 + (e/efit)²)，k·p 只在带边附近准），
    有 AMSET 能带时再做本征值精修（窗口不够点时放宽）。-> (model, 信息)。"""
    d = D0r.shape[0]
    basis = herm_basis(d)
    nq = Cq.shape[0]
    F = monomials(nodes["dk"], dim)
    A = herm_coords(np.einsum("jt,qtab->jqab", F, Cq), basis)              # (nn, nq, d²)
    y = herm_coords(Hr - D0r[None], basis)                                 # (nn, d²)
    emax = np.max(sign * np.linalg.eigvalsh(Hr), axis=1)
    wn = 1.0 / (1.0 + (np.maximum(emax, 0) / efit) ** 2)
    Aw = (A * np.sqrt(wn)[:, None, None]).transpose(0, 2, 1).reshape(-1, nq)
    g0 = np.linalg.lstsq(Aw, (y * np.sqrt(wn)[:, None]).reshape(-1), rcond=None)[0]
    info = {"n_nodes": int(len(y)), "n_params": int(nq), "rank_nodes": int(np.linalg.matrix_rank(Aw))}

    def mk(g):
        return {"dim": dim, "D0": D0r, "C": np.einsum("q,qtab->tab", g, Cq), "rmax": rcut}

    def node_rel(g):
        r = np.einsum("q,jqp->jp", g, A) - y
        return float(np.sqrt(np.sum(wn[:, None] * r ** 2)) / (np.sqrt(np.sum(wn[:, None] * y ** 2)) or 1.0))

    info["node_rel_res"] = node_rel(g0)
    sel = None
    if dense is not None:
        dkd, Ed = dense                                                    # (Nd,3)，(Nd, d) 升序、相对带边
        em = np.max(sign * Ed, axis=1)
        for mult in EFIT_STEPS:
            s = (np.linalg.norm(dkd, axis=1) < rcut) & (em < efit * mult)
            if s.sum() >= max(4 * nq, 20):
                sel, info["efit_fit_eV"] = s, efit * mult
                break
    if sel is None:
        info.update(rel_res=info["node_rel_res"], source="DFT 网格点 H_DFT")
        if dense is not None:
            info["dense_note"] = "AMSET 能带在窗口内的点太少"
        return mk(g0), info
    dkd, Ed = dkd[sel], Ed[sel]
    dHq = np.einsum("kt,qtab->kqab", monomials(dkd, dim), Cq)
    nd_rows = Ed.size
    lam = math.sqrt(LAMBDA_NODE / max(float(np.sum(wn)) * d * d, 1.0))
    sw = np.sqrt(wn)

    def res(g):
        e = np.linalg.eigvalsh(D0r + np.einsum("q,kqab->kab", g, dHq))
        r1 = (e - Ed).ravel() / math.sqrt(nd_rows)
        r2 = ((np.einsum("q,jqp->jp", g, A) - y) * sw[:, None]).ravel() * lam
        return np.concatenate([r1, r2])

    def jac(g):
        _e, V = np.linalg.eigh(D0r + np.einsum("q,kqab->kab", g, dHq))
        J1 = np.real(np.einsum("kan,kqab,kbn->knq", V.conj(), dHq, V)).reshape(-1, nq) / math.sqrt(nd_rows)
        J2 = (A * sw[:, None, None]).transpose(0, 2, 1).reshape(-1, nq) * lam
        return np.concatenate([J1, J2])

    g = lsq(res, g0, jac=jac)[0]
    e = np.linalg.eigvalsh(D0r + np.einsum("q,kqab->kab", g, dHq))
    info.update(n_dense=int(len(dkd)), source="AMSET 能带 + DFT 网格点 H_DFT",
                rel_res=float(np.linalg.norm(e - Ed) / (np.linalg.norm(Ed) or 1.0)), node_rel_res=node_rel(g))
    return mk(g), info


def split_groups(e, Sp, delta):
    """S7 点上哪些带能单独用：相邻两带（能量 e 升序）在 ±δ 应变下不交换次序 <=> 能隙 > δ·|S_n − S_{n+1}|（各图样取最大）。
    交换了的话排序后的中心差分是两带混在一起的值，只能用组内之和（与排序无关）；留 1.5 倍余量 + 1 meV。
    -> [(组, 权重)]：不交换但挨得近时，中心差分还带 (δ|X|/能隙)² 量级的三阶误差（带间应变耦合），
    按 x = δ|ΔS|/能隙 降权 1/(1 + (3x)⁴)（x = 1/3 时一半，临界 x = 2/3 时约 6%）。"""
    d = len(e)
    x = np.zeros(max(d - 1, 0))
    cut = np.zeros(max(d - 1, 0), bool)
    for n in range(1, d):
        dS = float(np.max(np.abs(Sp[n] - Sp[n - 1])))
        g = float(e[n] - e[n - 1])
        cut[n - 1] = g > 1.5 * delta * dS + 1e-3
        x[n - 1] = delta * dS / g if g > 0 else np.inf
    out, cur = [], [0]
    for n in range(1, d):
        if cut[n - 1]:
            out.append(cur)
            cur = [n]
        else:
            cur.append(n)
    out.append(cur)
    res = []
    for G in out:
        xs = [x[G[0] - 1]] if G[0] > 0 else []
        xs += [x[G[-1]]] if G[-1] < d - 1 else []
        xm = max(xs) if xs else 0.0
        res.append((G, 1.0 / (1.0 + (3.0 * xm) ** 4)))
    return res


def forward_k0(D0, X, delta):
    """k0 点：S7.1 / AMSET 的做法严格重演。±δ 应变下按能量排序 -> v_s = (E0 − E_s)/(sδ)；
    -> (S = ⟨v_s⟩ 带符号中心差分, A = ⟨|v_s|⟩ AMSET 幅值)，各 (d,)。"""
    e0 = np.linalg.eigvalsh(D0)
    vp = (e0 - np.linalg.eigvalsh(D0 + delta * X)) / delta
    vm = (e0 - np.linalg.eigvalsh(D0 - delta * X)) / (-delta)
    return 0.5 * (vp + vm), 0.5 * (np.abs(vp) + np.abs(vm))


def _fit_strain(model, Zq, dim, sign, bands, k0f, rlat, sign_h5, amp_h5, vac_axis, rcut, efit, rng_seed=0):
    """应变算符（对称允许的参数 γ）。-> (Z (np,d,d), 信息) 或 (None, 信息)。"""
    nx = Zq.shape[0]
    P_list = pairs(dim)
    delta = float(sign_h5.get("delta") or 0.005)
    S = np.asarray(sign_h5["S"])
    gap = sign_h5.get("gap")
    kS = np.asarray(sign_h5["kpoints"], float)
    bands = np.asarray(bands, int)
    info = {"n_params": int(nx), "delta": delta}
    Xq = {P: np.array([pattern_op(Zq[q], dim, P) for q in range(nx)]) for P in P_list}
    # (a) k0 点
    j0 = {_key(k): j for j, k in enumerate(kS)}.get(_key(k0f))
    S0 = A0 = None
    if j0 is not None:
        S0 = S[bands][:, j0]
        if amp_h5 is not None:
            jA = {_key(k): j for j, k in enumerate(np.asarray(amp_h5["kpoints"], float))}.get(_key(k0f))
            if jA is not None:
                A0 = np.array(np.asarray(amp_h5["A"])[bands][:, jA], float)
                if amp_h5.get("offdiag") == "tensor":
                    A0 = A0 * (2.0 - np.eye(3))                  # 张量约定 -> 原生（非对角 ×2）
    info["k0_in_s7"] = j0 is not None
    info["k0_amplitude"] = A0 is not None
    # (b) 邻近 S7 点：线性行
    dS, plane = _dk_cart(kS, k0f, rlat, dim, vac_axis)
    nS = np.linalg.norm(dS, axis=1)
    near = np.where(plane & (nS < rcut) & (nS > 1e-9))[0]
    Pidx = [(a, b) for a, b in P_list]
    rows, rhs, wts = [], [], []
    if len(near):
        e, V = np.linalg.eigh(model_H(model, dS[near]))
        grad = np.real(np.einsum("kin,kaij,kjn->kna", V.conj(), model_dH(model, dS[near]), V))
        for kk, j in enumerate(near):
            Sp = np.array([[S[b, j, a_, b_] for a_, b_ in Pidx] for b in bands])
            if not np.all(np.isfinite(Sp)):
                continue
            bad = set()
            if gap is not None:                                            # 子空间外的带挨得太近：排序会混进来
                for pos, out in ((0, bands[0] - 1), (len(bands) - 1, bands[-1] + 1)):
                    if not 0 <= out < S.shape[0]:
                        continue
                    g_in = (e[kk][1] - e[kk][0]) if pos == 0 else (e[kk][-1] - e[kk][-2]) if len(bands) > 1 else np.inf
                    g_f = float(np.asarray(gap)[bands[pos], j])
                    if g_f < g_in - 1e-4:
                        dSo = np.max(np.abs(Sp[pos] - np.array([S[out, j, a_, b_] for a_, b_ in Pidx])))
                        if g_f < 1.5 * delta * dSo + 1e-3:
                            bad.add(pos)
            for G, wx in split_groups(e[kk], Sp, delta):
                if bad & set(G):
                    continue
                Pg = V[kk][:, G] @ V[kk][:, G].conj().T
                eG = float(np.max(sign * e[kk][G]))
                wG = wx / (1.0 + (max(eG, 0.0) / efit) ** 2)
                for iP, P in enumerate(P_list):
                    Pm = np.zeros((3, 3))
                    Pm[P[0], P[1]] = Pm[P[1], P[0]] = 1.0
                    shift = float(np.einsum("ab,b,na->", Pm, dS[j], grad[kk][G]))
                    rows.append(np.real(np.einsum("ij,qji->q", Pg, Xq[P])))
                    rhs.append(-float(np.sum(Sp[G, iP])) + shift)
                    wts.append(wG)
    rows, rhs, wts = np.array(rows).reshape(-1, nx), np.array(rhs), np.array(wts)
    info["n_neighbor_rows"] = int(len(rhs))
    # 线性初值：k0 点的简并组之和（严格）+ 邻近行
    lr, lb, lw = [rows], [rhs], [wts]
    if S0 is not None:
        e0, V0 = np.linalg.eigh(model["D0"])
        for G, _w in split_groups(e0, np.array([[S0[n, a_, b_] for a_, b_ in Pidx] for n in range(len(bands))]), delta):
            Pg = V0[:, G] @ V0[:, G].conj().T
            for P in P_list:
                lr.append(np.real(np.einsum("ij,qji->q", Pg, Xq[P]))[None])
                lb.append(np.array([-float(np.sum(S0[G, P[0], P[1]]))]))
                lw.append(np.array([1.0]))
    Lr, Lb, Lw = np.concatenate(lr), np.concatenate(lb), np.concatenate(lw)
    if not len(Lb):
        info["status"] = "S7 网格在谷附近没有可用的点"
        return None, info
    sw = np.sqrt(Lw)
    g_lin = np.linalg.lstsq(Lr * sw[:, None], Lb * sw, rcond=None)[0]

    def res(g, parts=False):
        ra = []
        if S0 is not None:
            for P in P_list:
                Sm, Am = forward_k0(model["D0"], np.einsum("q,qab->ab", g, Xq[P]), delta)
                ra.append(Sm - S0[:, P[0], P[1]])
                if A0 is not None:
                    ra.append(Am - A0[:, P[0], P[1]])
        ra = np.concatenate(ra) if ra else np.zeros(0)
        rb = (rows @ g - rhs) * np.sqrt(wts) if len(rhs) else np.zeros(0)
        return (ra, rb) if parts else np.concatenate([ra, rb])

    data_a = []
    if S0 is not None:
        for P in P_list:
            data_a.append(S0[:, P[0], P[1]])
            if A0 is not None:
                data_a.append(A0[:, P[0], P[1]])
    na = float(np.linalg.norm(np.concatenate(data_a))) if data_a else 0.0
    nb_ = float(np.linalg.norm(rhs * np.sqrt(wts))) if len(rhs) else 0.0
    rng = np.random.default_rng(rng_seed)
    scale = (na + nb_) / math.sqrt(max(nx, 1)) + 1e-6
    sols = []
    for x0 in [g_lin] + [rng.normal(size=nx) * scale for _ in range(N_START)]:
        x, r, J = lsq(res, x0)
        sols.append((float(r @ r), x, J))
    sols.sort(key=lambda s: s[0])
    best_c, best_g, best_J = sols[0]
    sv = np.linalg.svd(best_J, compute_uv=False) if best_J.size else np.zeros(1)
    info["jac_cond"] = float(sv[0] / sv[-1]) if len(sv) >= nx and sv[-1] > 0 else float("inf")
    distinct, alt = [best_g], []
    for c, g, _J in sols[1:]:
        if all(np.linalg.norm(g - x) > 1e-2 * (np.linalg.norm(x) + 1e-9) for x in distinct):
            distinct.append(g)
            alt.append((c, g))
    tie = [(c, g) for c, g in alt if c <= best_c * (1 + AMBIG_TOL) + 1e-10 * (na + nb_ + 1) ** 2]
    info["n_minima"] = len(distinct)
    if alt:
        info["confidence"] = float((alt[0][0] - best_c) / alt[0][0]) if alt[0][0] > 0 else 0.0
    chosen = best_g
    if tie:                                    # 数据分不开的几个解：交给 build_valley 比较它们给出的放大系数
        info["candidates"] = [np.einsum("q,qtab->tab", g, Zq) for g in [best_g] + [g for _c, g in tie]]
    ra, rb = res(chosen, parts=True)
    info["k0_rel_res"] = float(np.linalg.norm(ra) / na) if na > 0 else None
    info["neighbor_rel_res"] = float(np.linalg.norm(rb) / nb_) if nb_ > 0 else None
    if info["jac_cond"] > 1e8:
        info["status"] = "S7 形变势定不住应变算符（雅可比条件数 %.1e）：k0 不在 S7 网格上或邻近点太少" % info["jac_cond"]
        return None, info
    if info["k0_rel_res"] is not None and info["k0_rel_res"] > STRAIN_REL_TOL:
        info["status"] = "k0 点的形变势与子空间模型对不上（相对残差 %.0f%% > %.0f%%）" % (
            100 * info["k0_rel_res"], 100 * STRAIN_REL_TOL)
        return None, info
    if S0 is None and info["neighbor_rel_res"] is not None and info["neighbor_rel_res"] > NB_REL_TOL:
        info["status"] = "k0 不在 S7 网格上，邻近点的形变势与模型对不上（相对残差 %.0f%%）" % (100 * info["neighbor_rel_res"])
        return None, info
    info["status"] = "ok"
    return np.einsum("q,qtab->tab", chosen, Zq), info


def _win_eff(S, b, j0, delta):
    """子空间的带间判据：至少 WIN_EV，且大于 k0 处应变造成的能移（排序才不会和子空间外的带混）。"""
    if j0 is None:
        return WIN_EV
    lo, hi = max(0, b - 4), min(S.shape[0], b + 5)
    return max(WIN_EV, 2 * delta * float(np.nanmax(np.abs(S[lo:hi, j0]))) + 0.002)


def build_valley(k0i, members, edge_band, sign, E, kfrac, coeff_fn, rlat, rots, cijkl, sign_h5, amp_h5,
                 dense, dim, vac_axis, efit, nang, rcut):
    """一个谷星 -> (谷表 dict 或 None, 报告)。"""
    k0f = np.asarray(kfrac[k0i], float)
    E_edge = float(E[edge_band, k0i])
    rep = {"k0": [float(x) for x in _wrap(k0f)], "n_members": len(members), "rcut": float(rcut)}
    delta = float(sign_h5.get("delta") or 0.005)
    j0 = {_key(k): j for j, k in enumerate(np.asarray(sign_h5["kpoints"], float))}.get(_key(k0f))
    win = _win_eff(np.asarray(sign_h5["S"]), edge_band, j0, delta)
    bands = manifold(E, k0i, edge_band, sign, win)
    d = len(bands)
    rep.update(bands=bands, d=d, win_eV=float(win), split_eV=[float(E[b, k0i] - E_edge) for b in bands])
    if d < 2:
        rep["status"] = "谷中心不简并（子空间只有 1 条带）：AMSET 的写法就是对的，不需要修正"
        return None, rep
    nodes = _node_data(k0i, bands, E, kfrac, rlat, coeff_fn, rcut, dim, vac_axis, E_edge)
    rep["n_nodes"] = int(len(nodes["idx"]))
    if len(nodes["idx"]) < 3:
        rep["status"] = "谷附近可用的 DFT 网格点太少（%d 个；投影权重 ≥ %.1f）" % (len(nodes["idx"]), W_MIN)
        return None, rep
    W, blocks = commutant_blocks(list(nodes["H"]), d)
    Hr = np.einsum("ai,kab,bj->kij", W.conj(), nodes["H"], W)
    Hr[:, ~_block_mask(blocks, d)] = 0.0
    D0r = W.conj().T @ np.diag(E[bands, k0i] - E_edge) @ W
    D0r[~_block_mask(blocks, d)] = 0.0
    D0r = 0.5 * (D0r + D0r.conj().T)
    rep["blocks"] = [len(b) for b in blocks]
    if max(len(b) for b in blocks) < 2:
        rep["status"] = "子空间在所有 k 上都分解成一维块（Kramers 简并或对称性禁止耦合）：不需要修正"
        return None, rep
    ops = _little_group(k0f, rlat, rots, dim)
    nt, np_ = n_terms(dim), len(pairs(dim))
    mbases, sbases, grep = [], [], []
    for b in blocks:
        ix = np.array(b)
        G, gi = _block_group(nodes, Hr[:, ix[:, None], ix[None, :]], k0f, rlat, ops, len(b))
        Ts = [term_T(M, dim) for M, _U, _a in G]
        mbases.append(equivariant_basis(G, len(b), Ts))
        sbases.append(equivariant_basis(G, len(b), [T[dim:, dim:] for T in Ts]))
        grep.append(dict(gi, n_kp_params=int(mbases[-1].shape[1]), n_strain_params=int(sbases[-1].shape[1])))
    rep["symmetry"] = grep
    Cq = _param_mats(blocks, d, mbases, nt)
    Zq = _param_mats(blocks, d, sbases, np_)
    dense_v = None
    if dense is not None:
        dkd, pl = _dk_cart(dense[0], k0f, rlat, dim, vac_axis)
        s = pl & (np.linalg.norm(dkd, axis=1) < rcut)
        dense_v = (dkd[s], np.asarray(dense[1])[bands][:, s].T - E_edge)
    model, kinfo = _fit_model(nodes, Hr, D0r, Cq, dim, sign, efit, dense_v, rcut)
    rep["kp"] = kinfo
    if kinfo["rel_res"] > FIT_REL_TOL:
        rep["status"] = "k·p 拟合残差 %.1f%% > %.0f%%（%s）" % (100 * kinfo["rel_res"], 100 * FIT_REL_TOL, kinfo["source"])
        return None, rep
    Z, sinfo = _fit_strain(model, Zq, dim, sign, bands, k0f, rlat, sign_h5, amp_h5, vac_axis, rcut, efit)
    rep["strain"] = sinfo
    if Z is None:
        rep["status"] = sinfo["status"]
        return None, rep
    energies = np.geomspace(E_MIN_EV, max(efit, 2 * E_MIN_EV), NE)
    offd = (amp_h5 or {}).get("offdiag", "tensor")
    bms = []
    for b in blocks:
        ix = np.array(b)
        bms.append({"dim": dim, "D0": D0r[np.ix_(ix, ix)], "C": model["C"][:, ix[:, None], ix[None, :]], "rmax": rcut})

    def tables_for(Zc):
        out = []
        for b, mb in zip(blocks, bms):
            ix = np.array(b)
            out.append(None if len(b) < 2 else
                       rate_ratio_table(mb, tensor_Y(Zc[:, ix[:, None], ix[None, :]], dim), cijkl, sign, energies, nang,
                                        offdiag=offd))
        return out

    cands = sinfo.pop("candidates", None) or [Z]
    alltabs = [tables_for(Zc) for Zc in cands]
    btabs = alltabs[0]
    if len(cands) > 1:
        # 数据分不开的几个应变算符：放大系数差得不多就取平均，差得多就不修正（不靠猜）
        means = np.array([np.concatenate([np.mean(t["ratio"], axis=(0, 1)) for t in tl if t is not None])
                          for tl in alltabs])                              # 各解、各块各带的平均放大系数
        spread = float(np.max((means.max(axis=0) - means.min(axis=0)) / means.mean(axis=0)))
        sinfo["ambiguous"] = {"n": len(cands), "ratio_means": means.tolist(), "spread": spread}
        if spread > AMBIG_SPREAD:
            rep["status"] = ("S7 形变势分不开 %d 个应变算符解（放大系数差 %.0f%%）：加密 S7 的 k 网格，让谷附近有非简并的点"
                             % (len(cands), 100 * spread))
            return None, rep
        btabs = [None if t is None else dict(t, ratio=np.mean([tl[i]["ratio"] for tl in alltabs], axis=0))
                 for i, t in enumerate(btabs)]
        sinfo["pick"] = "数据分不开 %d 个解（放大系数差 %.0f%%），取平均" % (len(cands), 100 * spread)
    rep["ratio_mean"] = [None if t is None else [float(np.mean(t["ratio"][:, :, n])) for n in range(t["ratio"].shape[2])]
                         for t in btabs]
    rep["ratio_range"] = [None if t is None else [float(np.min(t["ratio"])), float(np.max(t["ratio"]))] for t in btabs]
    Y = tensor_Y(Z, dim)
    rep["strain"]["Y_eig_eV"] = {"%d%d" % (a, b): np.round(np.linalg.eigvalsh(Y[a, b]), 4).tolist()
                                 for a, b in pairs(dim)}                       # 各分量在子空间里的本征值
    rep["status"] = "ok"
    valley = {"k0": k0f, "members": [(np.asarray(kfrac[m], float), np.asarray(M, float)) for m, M in members],
              "bands": [int(x) for x in bands], "sign": sign, "rcut": float(rcut), "dim": dim, "model": model,
              "Z": Z, "blocks": [list(b) for b in blocks], "block_models": bms, "tables": btabs}
    return valley, rep


# ---------------------------------------------------------------------------
# 总入口：DFT 全网格 + 波函数 + S7.1 形变势 -> 修正表
# ---------------------------------------------------------------------------
def setup(kfrac, E, coeff_fn, vb, rlat, rots, cijkl, sign_h5, amp_h5=None, dense=None, dim=3, vac_axis=2,
          efit=EFIT_EV, nang=None):
    """-> (修正表 dict 或 None, 报告 dict)。
    kfrac (N,3)：DFT 全网格（= AMSET 波函数的网格）；E (nb, N) eV，带号与 AMSET 相同（能量窗口里的带）；vb 价带顶带号；
    coeff_fn(带号列表, kfrac (n,3), 参考 kfrac) -> 逐对 (带, k) 的归一化波函数系数（已按参考点对齐倒格矢）；
    rlat (3,3) 倒格矢（行，任意单位）；rots 点群的笛卡尔旋转（含非真旋转，不含时间反演）；cijkl (3,3,3,3) 弹性张量（只用比例）；
    sign_h5 {"kpoints", "S" (nb, M, 3, 3) 原生约定, "delta", "gap"?}；amp_h5 {"kpoints", "A", "offdiag"} 运行期幅值文件（可缺）；
    dense (kfrac (Nd,3), E (nb, Nd) eV)：AMSET 插值后的能带（可缺）；dim 3 / 2（2D 的真空方向 = 分数坐标 vac_axis，笛卡尔 z）。"""
    E = np.asarray(E, float)
    kfrac = np.asarray(kfrac, float)
    rlat = np.asarray(rlat, float)
    rots = [np.asarray(R, float) for R in rots]
    if dim == 2:
        rots = [R for R in rots if layer_op(R)]
    axes = [a for a in range(3) if dim == 3 or a != vac_axis]
    bmin = float(np.min(np.linalg.norm(rlat[axes], axis=1)))
    ops_star = [s * R for R in rots for s in (1.0, -1.0)]
    report = {"dim": dim, "valleys": []}
    valleys = []
    if np.asarray(sign_h5["S"]).shape[0] != E.shape[0]:
        report["status"] = "带符号形变势的带数（%d）与 AMSET 的 %d 不同" % (np.asarray(sign_h5["S"]).shape[0], E.shape[0])
        return None, report
    for carrier, b, sign in (("hole", vb, -1), ("electron", vb + 1, +1)):
        if b < 0 or b >= E.shape[0]:
            continue
        ext = float(np.max(E[b]) if sign < 0 else np.min(E[b]))
        try:
            cand = local_extrema(E[b], kfrac, sign, axes)
        except ValueError as e:
            report["status"] = str(e)
            return None, report
        cand = sorted((int(i) for i in cand if sign * (E[b, i] - ext) < efit), key=lambda i: sign * E[b, i])
        stars = group_stars(cand, kfrac, rlat, ops_star, Eb=E[b])
        centers = np.array([kfrac[m] for _r, ms in stars for m, _M in ms])
        for r, ms in stars:
            dk, _pl = _dk_cart(centers, kfrac[r], rlat, dim, vac_axis)
            dist = np.linalg.norm(dk, axis=1)
            dist = dist[dist > 1e-9]
            rcut = min(KCUT_FRAC * bmin, 0.5 * float(np.min(dist)) if len(dist) else np.inf)
            v, vrep = build_valley(r, ms, b, sign, E, kfrac, coeff_fn, rlat, rots, cijkl, sign_h5, amp_h5,
                                   dense, dim, vac_axis, efit, nang, rcut)
            vrep.update(carrier=carrier, E_rel_edge_eV=float(E[b, r] - ext))
            report["valleys"].append(vrep)
            if v is not None:
                valleys.append(v)
    if not valleys:
        why = [v["status"] for v in report["valleys"]]
        report["status"] = "；".join(sorted(set(why))) if why else "没有找到带边谷"
        return None, report
    report["status"] = "ok"
    owner, cen = [], []
    for vi, v in enumerate(valleys):
        for mi, (kf, _M) in enumerate(v["members"]):
            owner.append((vi, mi))
            cen.append(kf)
    table = {"rlat": rlat, "dim": dim, "vac_axis": vac_axis, "valleys": valleys,
             "centers": np.array(cen), "owner": owner}
    return table, report


def scale(table, band_idx, kfrac):
    """初态 (带, k) 的 ADP 率要乘的系数：不在修正过的谷 / 子空间块里、或离谷中心超过 rcut -> 1。"""
    if not table:
        return 1.0
    dim = table["dim"]
    dk, plane = _dk_cart(np.asarray(kfrac, float)[None], table["centers"], table["rlat"], dim, table["vac_axis"])
    dist = np.where(plane, np.linalg.norm(dk, axis=1), np.inf)
    j = int(np.argmin(dist))
    vi, mi = table["owner"][j]
    v = table["valleys"][vi]
    r = float(dist[j])
    if band_idx not in v["bands"] or r >= v["rcut"]:
        return 1.0
    dkr = v["members"][mi][1].T @ dk[j]
    if dim == 2:
        dkr[2] = 0.0
    e, V = np.linalg.eigh(model_H(v["model"], dkr))
    i = v["bands"].index(band_idx)
    bi = int(np.argmax([float(np.sum(np.abs(V[b, i]) ** 2)) for b in v["blocks"]]))
    tab = v["tables"][bi]
    if tab is None:
        return 1.0
    eb = np.linalg.eigvalsh(model_H(v["block_models"][bi], dkr))
    re_ = int(np.argmin(np.abs(eb - e[i])))
    rank = re_ if v["sign"] > 0 else len(eb) - 1 - re_
    eps = max(v["sign"] * float(e[i]), 0.0)
    R = np.asarray(tab["ratio"])[:, :, rank]
    if r < 1e-12:
        col = R.mean(axis=1)
    else:
        col = R[:, int(np.argmax(np.asarray(tab["dirs"]) @ (dkr / np.linalg.norm(dkr))))]
    Eg = np.asarray(tab["E"])
    s = float(np.interp(np.log(max(eps, Eg[0])), np.log(Eg), col))
    w = min(max((v["rcut"] - r) / (TAPER * v["rcut"]), 0.0), 1.0)
    return 1.0 + w * (s - 1.0)


# ---------------------------------------------------------------------------
# 运行期（插件调用）：AMSET 数据 + vasprun + S7.1 文件
# ---------------------------------------------------------------------------
def read_signed(path, spin_name):
    import h5py
    with h5py.File(str(path), "r") as f:
        key = "signed_" + spin_name
        if key not in f:
            raise KeyError("%s 里没有 %s" % (path, key))
        out = {"S": np.array(f[key], float), "kpoints": np.array(f["kpoints"], float),
               "delta": float(f.attrs.get("strain_amount") or 0.005)}
        if ("gap_" + spin_name) in f:
            out["gap"] = np.array(f["gap_" + spin_name], float)
    return out


def read_amp(path, spin_name):
    import h5py
    with h5py.File(str(path), "r") as f:
        key = next((k for k in f.keys() if k.startswith("deformation_potentials") and k.endswith(spin_name)), None)
        if key is None:
            raise KeyError("%s 里没有 deformation_potentials*%s" % (path, spin_name))
        return {"A": np.array(f[key], float), "kpoints": np.array(f["kpoints"], float),
                "offdiag": "tensor" if str(f.attrs.get("dp_offdiag", "")) == "tensor" else "native"}


def point_rotations(structure, symprec=0.01):
    """点群操作的笛卡尔矩阵（含非真旋转）。"""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    out = []
    for op in SpacegroupAnalyzer(structure, symprec=symprec).get_point_group_operations(cartesian=True):
        R = np.asarray(op.rotation_matrix, float)
        if not any(np.allclose(R, x, atol=1e-6) for x in out):
            out.append(R)
    return out


def load_dft_nodes(path="vasprun.xml", symprec=0.01):
    """vasprun.xml（AMSET 的输入，与 wavefunction.h5 同源）-> DFT 全网格能带：{"kpoints", "E": {spin: (nb, N) eV}, "vb": {spin: 带号}}。"""
    from pymatgen.io.vasp import Vasprun
    from amset.electronic_structure.common import get_band_structure
    from amset.electronic_structure.symmetry import expand_bandstructure
    vr = Vasprun(str(path), parse_dos=False, parse_projected_eigen=False, parse_potcar_file=False)
    bs = expand_bandstructure(get_band_structure(vr), symprec=symprec)
    vbm = bs.get_vbm()
    out = {"kpoints": np.array([k.frac_coords for k in bs.kpoints], float), "E": {}, "vb": {}}
    for spin, b in bs.bands.items():
        out["E"][spin.name] = np.asarray(b, float)
        idx = (vbm.get("band_index") or {}).get(spin) or []
        out["vb"][spin.name] = (int(max(idx)) if len(idx)
                                else int(np.sum(np.max(b, axis=1) <= vbm["energy"] + 1e-6)) - 1)
    return out


def make_coeff_fn(oc, spin):
    """AMSET WavefunctionOverlapCalculator -> coeff_fn（在网格点上线性插值 = 节点值；按参考点对齐倒格矢）。"""
    g_maps = np.asarray(oc.g_maps)

    def fn(bands, kf, kref):
        kf = np.asarray(kf, float)
        c = np.array(oc.get_coefficients(spin, np.asarray(bands, int), kf))
        g = np.round(kf - np.asarray(kref, float)).astype(int)
        for i in np.where(np.any(g != 0, axis=1))[0]:
            c[i] = align_coefficients(c[i], g[i], g_maps)
        return c
    return fn


def runtime_tables(amset_data, rec, cijkl, dim=3, tag="[adp_degen]", report_path="adp_degen_fit.json"):
    """插件入口（主进程 from_amset_data 里调一次）：-> {spin.name: 修正表} 或 None。
    rec = gen 写的 "adp_degen" 记录（mode / signed_h5 / deform_h5 / efit_eV / nang / vasprun）。
    不适用时：mode = on 报错，auto 打印原因、照原版。报告写 report_path。"""
    import json
    mode = str(rec.get("mode") or "auto").lower()
    efit = float(rec.get("efit_eV") or EFIT_EV)
    nang = rec.get("nang")
    reports, tables = {}, {}
    oc = getattr(amset_data, "overlap_calculator", None)
    why, dft, rlat, rots = None, None, None, None
    if oc is None or not hasattr(oc, "get_coefficients") or not hasattr(oc, "g_maps"):
        why = "AMSET 没用波函数重叠（unity_overlap / use_projections）：读不到子空间结构"
    elif getattr(amset_data, "is_metal", False) or not getattr(amset_data, "vb_idx", None):
        why = "金属：没有带边"
    if why is None:
        try:
            dft = load_dft_nodes(rec.get("vasprun") or "vasprun.xml")
        except Exception as e:                                         # noqa: BLE001
            why = "读不了 vasprun.xml（%s: %s）" % (type(e).__name__, e)
    if why is None:
        rlat = np.asarray(amset_data.structure.lattice.reciprocal_lattice.matrix, float)
        try:
            rots = point_rotations(amset_data.structure)
        except Exception as e:                                         # noqa: BLE001
            rots = [np.eye(3)]
            print("%s[WARN] 取不到点群（%s），按无对称性处理" % (tag, e), flush=True)
    for spin in amset_data.spins:
        if why is not None:
            reports[spin.name] = {"status": why}
            continue
        try:
            vb = int(amset_data.vb_idx[spin])
            Ea = np.asarray(amset_data.energies[spin], float) * HARTREE_EV
            nb = Ea.shape[0]
            off = dft["vb"][spin.name] - vb
            Ed = dft["E"][spin.name]
            if off < 0 or off + nb > Ed.shape[0]:
                raise ValueError("DFT 带号对不上 AMSET 的能量窗口（偏移 %d，AMSET %d 条，DFT %d 条）" % (off, nb, Ed.shape[0]))
            En = Ed[off:off + nb]
            if abs(np.max(En[vb]) - np.max(Ea[vb])) > 0.05:
                raise ValueError("DFT 与 AMSET 的价带顶差 %.3f eV：带号对不上" % (np.max(En[vb]) - np.max(Ea[vb])))
            sgn = read_signed(rec.get("signed_h5") or "deformation_signed.h5", spin.name)
            try:
                amp = read_amp(rec.get("deform_h5") or "deformation.h5", spin.name)
            except (OSError, KeyError):
                amp = None
            t, rp = setup(dft["kpoints"], En, make_coeff_fn(oc, spin), vb, rlat, rots, np.asarray(cijkl, float), sgn,
                          amp, dense=(np.asarray(amset_data.kpoints, float), Ea), dim=dim, efit=efit, nang=nang)
        except Exception as e:                                         # noqa: BLE001
            t, rp = None, {"status": "%s: %s" % (type(e).__name__, e)}
        reports[spin.name] = rp
        if t is not None:
            tables[spin.name] = t
    try:
        with open(report_path, "w") as fh:
            json.dump({"mode": mode, "dim": dim, "spins": reports}, fh, indent=1, ensure_ascii=False, default=_jdef)
    except OSError:
        pass
    for sn, rp in reports.items():
        if rp.get("status") == "ok":
            for v in rp["valleys"]:
                if v.get("status") == "ok":
                    print("%s %s" % (tag, summary_line(sn, v)), flush=True)
        else:
            msg = "%s ADP 简并修正（%s）不适用：%s" % (tag, sn, rp.get("status"))
            if mode == "on":
                raise RuntimeError(msg + "（step.conf ADP_DEGEN = on；确认不需要就写 auto 或 off）")
            print(msg + "（ADP_DEGEN = auto：照原版 ADP 算）", flush=True)
    return tables or None


def summary_line(spin_name, v):
    """一个修正过的谷的一行摘要（插件日志与 S8.3 汇总共用）。"""
    kp, st = v["kp"], v["strain"]
    rm = [x for r in (v.get("ratio_mean") or []) if r for x in r]
    return ("ADP 简并修正（%s，%s谷 k0=(%s)，%d 条带，分块 %s）：k·p 残差 %.1f%%，k0 形变势残差 %s，散射率放大 %s%s" % (
        spin_name, "空穴" if v.get("carrier") == "hole" else "电子", ",".join("%.3f" % x for x in v["k0"]), v["d"],
        v["blocks"], 100 * kp["rel_res"],
        "—" if st.get("k0_rel_res") is None else "%.1f%%" % (100 * st["k0_rel_res"]),
        " / ".join("%.2f" % x for x in rm), "；" + st["pick"] if st.get("pick") else ""))


def _jdef(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, complex):
        return [o.real, o.imag]
    return str(o)
