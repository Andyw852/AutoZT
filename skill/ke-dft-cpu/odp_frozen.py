# -*- coding: utf-8 -*-
"""odp_frozen.py —— 冻声子法自动提取谷内光学形变势（ODP）的纯 numpy 部分（[0035]）。

S7.2（gen_step9c_odp.py）用它选模、生成位移；S8.21（gen_step12c_odp_read.py）用它从各构型的本征值
算 D；S8.4（gen_step14_amset2d.py）用它把 odp_params.json 变成 dp_optical 的谷内通道。
不依赖 pymatgen/amset，单元测试直接测它（tests/test_0035_odp_frozen.py）。

物理（与 dp_optical.py 的常数形变势模型同一个归一化）：
    dp_optical 的谷内率  1/τ = π D²/(M ω)·[N·G(x+ħω) + (N+1)·G(x−ħω)]
    来自 |g|² = ħ D² / (2 N M ω)，M = 原胞质量。Γ 点声子 ν 的位移（e_κ 为质量加权、归一的本征矢，Σ_κ|e_κ|² = 1）：
        Δ_κ = u · √(M/m_κ) · e_κ                                    （u 是长度）
    此时质量加权的简正坐标 Q = Σ_κ √m_κ e_κ·Δ_κ = u√M，故 dE/dQ = (dE/du)/√M，代入 |g|² = ħ|dE/dQ|²/(2Nω)
    正好得到上式 —— 所以 D = dE/du（eV/Å）直接能用，不用再换算。纯平移（e_κ ∝ √m_κ）时 Δ_κ = u（全部原子平移 u）。
    D 用真空能级对齐的带边能量（E − V_vac），与 S7.1/S8.2 的 E1 同一口径。

    非简并带：中心差分 D = [E(+u) − E(−u)] / 2u，二次项相消。对称性禁戒的模（该谷的小群下不是全对称表示）自动给 D ≈ 0。
    简并带（d 重，u = 0 时相差 < DEG_TOL_EV）：微扰把它们劈裂成 E₀ + λ_i u。按能量排序后
        中心差分会把 +λ 和 −λ 配错对（得到 0），所以改用迹恒等式（与排序无关）：
            Σλ² = Tr(V²) = [Σ_i(E_i(+u)−ē)² + Σ_i(E_i(−u)−ē)² − 2Σ_i(E_i(0)−ē)²] / 2u²
        dp_optical 的目标谷态密度 G 统计该谷的全部能带（d 条都算），每个态的散射率是 Σ_j|M_ij|²·G/d 的平均，
        即 D_eff² = Tr(M²)/d² = Σλ²/d²（d = 1 时退回 λ²）。
    组外邻带近到位移下可能交叉时做核对（排序配对与迹恒等式一致才算数），不一致就报错要求减小位移：
        非简并带的 D 是对角元 ⟨n|∂H/∂u|n⟩（小位移极限），交叉/强混合时它与位移大小有关，不能硬给一个数。
    同一组简并声子（ħω 相同）的各个伙伴模的 D² 相加 —— 和与简并子空间里本征矢的取法无关（正交变换不变）。
    星的各成员谷：对称等价，Σ_伙伴 D² 应相同；取平均并记录离散度。

极性模（红外活性）默认不算：它们的长程（Fröhlich）耦合由 POP 管；短程部分要算就设 ODP_POLAR = include。

[0036] 三维：模型与归一化不变（上式与维数无关）。区别只有两处：
  · 没有真空：本征值不对齐、直接差分（dvac = 0）。Γ 点光学模不改变晶胞，VASP 的平均静电势（G = 0 分量）约定在 ±u
    构型里相同；±u 互为对称像的模（如 Si 的 Γ25′）参考点的变化是 u 的偶函数，中心差分里本来就相消。
  · 没有 S8.2 多谷谷表（只实现了 2D）：find_valleys_3d 直接在 S7 undeformed 的网格（ISYM=0，全网格、仅按 ±k 约化）上找
    窗口内的带边极值点，按（能量、带、自旋）归星；或者用 step.conf 手填的 DP_OPTICAL_VALLEYS_E/H。
"""
import hashlib
import json
import math

import numpy as np

THZ_TO_MEV = 4.135667696          # h·1 THz = 4.1357 meV
DEG_TOL_EV = 0.005                # 带简并判据：u = 0 时同一 k 点、同一自旋两条带相差小于它
GAP_MARGIN_EV = 0.002             # 简并组外最近的带离得太近（< 2·最大位移 + 它）就做交叉核对
CROSS_REL = 0.05                  # 交叉核对：排序配对平方和与迹恒等式的相对差上限
CROSS_ABS = 0.05                  # 交叉核对：绝对差上限（(eV/Å)²）
FREQ_DEG_THZ = 0.05               # 声子简并判据（ħω 相差 < 0.2 meV）
ACOUSTIC_WEIGHT = 0.5             # 本征矢在平移子空间上的投影权重 > 它 = 声学模
POLAR_TOL_E = 0.1                 # 模有效电荷 |Σ Z*·U|/|U|（e）> 它 = 极性模
# [0035-r3] 面外分量单独判：面外偶极弱（MoS₂ 的 A2″ 总模有效电荷 < 0.1 e）也要认出来。它不在二维 POP 里（POP 只用面内
#   极性模），有 σh 的材料里它对谷内一阶耦合是对称性禁戒的；位移还会产生面外偶极、真空区出现斜坡。默认不算。
#   对称的非极性模（A1′/E″）的面外偶极按对称性精确为零（数值上 ~1e-3 e），0.02 e 留了一个数量级的余量。
POLAR_Z_TOL_E = 0.02
FMIN_THZ = 0.5                    # 低于它的光学模（含虚频）不算：ħω -> 0 时 1/ω 发散，参数也不可信
DMIN_EV_A = 0.05                  # S8.4 用：D 小于它的通道不进 SCATTERING（对称性禁戒的模，贡献 < 0.1%）
SPREAD_WARN = 0.10                # 星内各成员谷 Σ D² 的相对离散度超过它 -> 告警（对称等价的谷应相同）
SCHEMA = "odp_params/1"


# ---------------------------------------------------------------------------
# 声子模：分类、简并分组、选模
# ---------------------------------------------------------------------------
def check_eigvecs(eigvecs, tol_norm=1e-3, tol_orth=1e-2):
    """质量加权本征矢 (nmode, nat, 3)：每个归一、两两正交。不满足 -> ValueError（多半拿成了"除以 √m"的那组）。"""
    e = np.asarray(eigvecs, float)
    if e.ndim != 3 or e.shape[2] != 3 or e.shape[0] != 3 * e.shape[1]:
        raise ValueError("本征矢形状 %s 不是 (3N, N, 3)" % (e.shape,))
    flat = e.reshape(e.shape[0], -1)
    gram = flat @ flat.T
    dn = float(np.max(np.abs(np.diag(gram) - 1.0)))
    off = gram - np.diag(np.diag(gram))
    do = float(np.max(np.abs(off))) if off.size else 0.0
    if dn > tol_norm or do > tol_orth:
        raise ValueError("本征矢不是正交归一的（|e|²−1 最大 %.2e，交叠最大 %.2e）：应当用 vasprun.xml 的"
                         "质量加权本征矢（dynmat eigenvectors），不是 OUTCAR 里除过 √m 的那组" % (dn, do))
    return dn, do


def mode_table(freqs_thz, eigvecs, masses_amu, born=None, polar_tol=POLAR_TOL_E, fmin_thz=FMIN_THZ,
               deg_thz=FREQ_DEG_THZ, include_polar=False, polar_z_tol=POLAR_Z_TOL_E):
    """Γ 点声子模 -> (modes, groups, notes)。按频率升序编号（index 从 1 起）。

    modes[i]：index、freq_THz、hw_meV、trans_weight（平移子空间投影）、acoustic、Z_mode_e（模有效电荷，无 Born 电荷为 None）、
              Z_mode_z_e（面外分量）、polar、group（声学模为 None）。
    groups[j]：group、modes（index 列表）、freq_THz、hw_meV、polar、use（是否要算）、why（不算的原因）。"""
    f = np.asarray(freqs_thz, float)
    e = np.asarray(eigvecs, float)
    m = np.asarray(masses_amu, float)
    nmode, nat = e.shape[0], e.shape[1]
    if len(f) != nmode or len(m) != nat:
        raise ValueError("频率 %d 个、本征矢 %d 个、原子 %d 个对不上" % (len(f), nmode, len(m)))
    check_eigvecs(e)
    order = np.argsort(f, kind="stable")
    # 平移子空间：T_α = √m_κ ê_α / √M（三个正交归一矢量）
    sq = np.sqrt(m / m.sum())
    trans = np.zeros((3, nat, 3))
    for a in range(3):
        trans[a, :, a] = sq
    tw = np.einsum("nka,bka->nb", e, trans) ** 2
    tw = tw.sum(axis=1)
    notes = []
    acoustic = np.zeros(nmode, bool)
    cand = [i for i in range(nmode) if tw[i] > ACOUSTIC_WEIGHT]
    if len(cand) != 3:
        cand = list(np.argsort(-tw)[:3])
        notes.append("平移投影 > %.1f 的模不是正好 3 个，按投影最大的 3 个当声学模（%s）"
                     % (ACOUSTIC_WEIGHT, ", ".join("%.3f" % tw[i] for i in cand)))
    acoustic[cand] = True
    zmode = zz = None
    if born is not None:
        z = np.asarray(born, float)
        if z.shape != (nat, 3, 3):
            raise ValueError("Born 有效电荷形状 %s 不是 (N, 3, 3)" % (z.shape,))
        U = e / np.sqrt(m)[None, :, None]
        dip = np.einsum("kab,nkb->na", z, U)
        un = np.sqrt(np.einsum("nka,nka->n", U, U))
        zmode = np.linalg.norm(dip, axis=1) / un
        zz = np.abs(dip[:, 2]) / un
    modes = []
    for rank, i in enumerate(order):
        modes.append({"index": rank + 1, "freq_THz": float(f[i]), "hw_meV": float(f[i] * THZ_TO_MEV),
                      "trans_weight": float(tw[i]), "acoustic": bool(acoustic[i]),
                      "Z_mode_e": None if zmode is None else float(zmode[i]),
                      "Z_mode_z_e": None if zz is None else float(zz[i]),
                      "polar": bool(zmode is not None and (zmode[i] > polar_tol or zz[i] > polar_z_tol)),
                      "polar_kind": (None if zmode is None or not (zmode[i] > polar_tol or zz[i] > polar_z_tol)
                                     else ("面外" if zz[i] > polar_z_tol and zz[i] >= 0.5 * zmode[i]
                                           else ("体相" if math.isinf(polar_z_tol) else "面内"))),   # [0036] 3D 不分面内外
                      "group": None, "_src": int(i)})
    # 简并分组：频率相差 < deg_thz 且极性相同（真简并的伙伴属于同一个不可约表示，极性必相同；
    #   偶然近简并的极性/非极性两个模不能并组，否则非极性的那个会跟着被当成极性剔掉）。
    #   分组只影响汇报与选模的整组补全：通道里各模的 D² 本来就按 ħω 相加，拆开或并起来率都一样。
    groups = []
    for md in modes:
        if md["acoustic"]:
            continue
        if (groups and abs(md["freq_THz"] - groups[-1]["_flast"]) < deg_thz
                and md["polar"] == groups[-1]["_polar"]):
            g = groups[-1]
            g["modes"].append(md["index"])
            g["_f"].append(md["freq_THz"])
            g["_flast"] = md["freq_THz"]
        else:
            groups.append({"group": len(groups) + 1, "modes": [md["index"]], "_f": [md["freq_THz"]],
                           "_flast": md["freq_THz"], "_polar": md["polar"]})
        md["group"] = groups[-1]["group"]
    by_index = {md["index"]: md for md in modes}
    for g in groups:
        fm = float(np.mean(g.pop("_f")))
        g.pop("_flast")
        g.pop("_polar")
        g["freq_THz"] = fm
        g["hw_meV"] = fm * THZ_TO_MEV
        g["polar"] = any(by_index[i]["polar"] for i in g["modes"])
        if born is None:
            g["polar"] = None
        why = None
        if min(by_index[i]["freq_THz"] for i in g["modes"]) <= fmin_thz:
            why = "频率 %.3f THz ≤ %.2f（虚频或太软）" % (fm, fmin_thz)
        elif g["polar"] is None and not include_polar:
            why = "没有 Born 有效电荷，判不了极性（ODP_POLAR=include 则不判、全算）"
        elif g["polar"] and not include_polar:
            kinds = {by_index[i]["polar_kind"] for i in g["modes"]}
            zt = max(by_index[i]["Z_mode_e"] for i in g["modes"])
            zz = max(by_index[i]["Z_mode_z_e"] for i in g["modes"])
            if kinds == {"面外"}:
                why = ("面外极性模（面外模有效电荷 %.3f e > %.2f）：二维 POP 不含它；有 σh 的材料里它对谷内一阶耦合是"
                       "对称性禁戒的（D≈0），位移还会产生面外偶极（ODP_POLAR=include 才算）" % (zz, polar_z_tol))
            else:
                why = "极性模（模有效电荷 %.2f e > %.2f）：长程耦合由 POP 管（ODP_POLAR=include 才算）" % (zt, polar_tol)
        g["use"] = why is None
        g["why"] = why
    return modes, groups, notes


def displacement(e_mode, masses_amu, disp_max_a):
    """一个模的冻声子位移：Δ_κ = u·√(M/m_κ)·e_κ，u 取到最大的原子位移 = disp_max_a。返回 (Δ (nat,3) Å, u Å)。"""
    e = np.asarray(e_mode, float)
    m = np.asarray(masses_amu, float)
    s = np.sqrt(m.sum() / m)[:, None] * e
    amp = float(np.max(np.linalg.norm(s, axis=1)))
    if not amp > 0:
        raise ValueError("本征矢为零")
    u = float(disp_max_a) / amp
    return u * s, u


# ---------------------------------------------------------------------------
# 带边 D：一个 k 点、一个自旋
# ---------------------------------------------------------------------------
class BandCrossing(ValueError):
    """简并组外的邻带在 ±u 构型里与本组交叉/强混合：这个位移下单带的 D 解不出来（与 u 有关），要减小位移重算。"""


def trace_sq(e0, ep, em, u):
    """一组能带的 Σλ²（迹恒等式，与排序、交叉无关）：Σ_i(E_i(+u)−ē)² + Σ_i(E_i(−u)−ē)² − 2Σ_i(E_i(0)−ē)² 除以 2u²，
    ē = 该组 u = 0 的平均能量。H(u) = H₀ + uV 时它等于 Tr(V²)（含组内的非对角耦合），误差 O(u²)。"""
    e0, ep, em = (np.asarray(x, float) for x in (e0, ep, em))
    c = float(np.mean(e0))
    return float((np.sum((ep - c) ** 2) + np.sum((em - c) ** 2) - 2.0 * np.sum((e0 - c) ** 2)) / (2.0 * u * u))


def _sorted_sq(e0, ep, em, u):
    """按能量排序逐条配对的斜率平方和（两侧平均）。不交叉、组内没有强非对角耦合时 = trace_sq（差 O(u²)）。"""
    sp = (np.asarray(ep, float) - np.asarray(e0, float)) / u
    sm = (np.asarray(em, float) - np.asarray(e0, float)) / (-u)
    return 0.5 * (float(np.sum(sp * sp)) + float(np.sum(sm * sm)))


def edge_coupling(e0, ep, em, band, u, dvac_p=0.0, dvac_m=0.0, deg_tol=DEG_TOL_EV, margin=GAP_MARGIN_EV):
    """e0/ep/em：该 k 点、该自旋全部能带的本征值（eV，升序）——u = 0、+u、−u 三个构型；
    dvac_p/dvac_m：+u/−u 构型的真空能级相对 u = 0 的变化（eV），从能量里扣掉（真空对齐）。

    简并组 P = u = 0 时与 band 相差 ≤ deg_tol 的带（物理判据，与位移大小无关）：
      d = 1：D² = λ²，λ = 中心差分；d > 1：D² = Σλ²/d²，Σλ² 用迹恒等式（与组内排序无关）。
    组外邻带离得比 2×最大带边移动 + margin 还近（位移下可能交叉）时，把它们一起算一遍：排序配对的平方和与迹恒等式
    一致（差 ≤ CROSS_REL·值 + CROSS_ABS）= 没交叉、没有强组间混合，P 的结果可信；不一致 -> 抛 BandCrossing
    （这个位移下 P 的 D 与 u 有关，要减小 ODP_DISP_A 重算 S7.2）。
    返回 dict：D2、d、bands [lo, hi]（P）、near [lo, hi]（参与交叉核对的范围）、method、lam、shift_max_eV、gap_out_eV。"""
    e0 = np.asarray(e0, float)
    ep = np.asarray(ep, float) - float(dvac_p)
    em = np.asarray(em, float) - float(dvac_m)
    nb = len(e0)
    b = int(band)
    if not (0 <= b < nb) or len(ep) != nb or len(em) != nb:
        raise ValueError("带号 %d 越界或三个构型的带数不同（%d/%d/%d）" % (b, nb, len(ep), len(em)))
    lo = hi = b
    while lo > 0 and abs(e0[lo - 1] - e0[b]) <= deg_tol:
        lo -= 1
    while hi < nb - 1 and abs(e0[hi + 1] - e0[b]) <= deg_tol:
        hi += 1
    idx = slice(lo, hi + 1)
    d = hi - lo + 1
    shift = float(max(np.max(np.abs(ep[idx] - e0[idx])), np.max(np.abs(em[idx] - e0[idx]))))
    wlo, whi = lo, hi
    for _ in range(nb):
        w = slice(wlo, whi + 1)
        sw = float(max(np.max(np.abs(ep[w] - e0[w])), np.max(np.abs(em[w] - e0[w]))))
        need = 2.0 * sw + margin
        grow = False
        if wlo > 0 and e0[wlo] - e0[wlo - 1] < need:
            wlo -= 1
            grow = True
        if whi < nb - 1 and e0[whi + 1] - e0[whi] < need:
            whi += 1
            grow = True
        if not grow:
            break
    gaps = []
    if lo > 0:
        gaps.append(float(e0[lo] - e0[lo - 1]))
    if hi < nb - 1:
        gaps.append(float(e0[hi + 1] - e0[hi]))
    out = {"d": d, "bands": [lo, hi], "near": [wlo, whi], "shift_max_eV": shift,
           "gap_out_eV": min(gaps) if gaps else None, "widened": (wlo, whi) != (lo, hi)}
    if out["widened"]:
        w = slice(wlo, whi + 1)
        s_w, t_w = _sorted_sq(e0[w], ep[w], em[w], u), trace_sq(e0[w], ep[w], em[w], u)
        out["near_check"] = {"sorted": s_w, "trace": t_w}
        if abs(s_w - t_w) > CROSS_REL * max(abs(s_w), abs(t_w)) + CROSS_ABS:
            raise BandCrossing(
                "带 %d 附近的带 %d–%d 在 ±u 构型里交叉或强混合（排序配对 Σλ² = %.3g，迹恒等式 %.3g (eV/Å)²；"
                "带边移动 %.0f meV，邻带间隔 %.0f meV）：这个位移下单带的 D 解不出来"
                % (b, wlo, whi, s_w, t_w, 1000 * shift, 1000 * min(gaps) if gaps else 0.0))
    if d == 1:
        lam = (ep[b] - em[b]) / (2.0 * u)
        out.update(D2=float(lam * lam), method="central", lam=[float(lam)])
        return out
    s2 = trace_sq(e0[idx], ep[idx], em[idx], u)
    sp = (ep[idx] - e0[idx]) / u
    out.update(D2=s2 / (d * d), method="degenerate(%d)" % d, lam=[float(x) for x in sp],
               trace_pm=[float(np.sum(sp)), float(np.sum((em[idx] - e0[idx]) / (-u)))])
    return out


# ---------------------------------------------------------------------------
# 真空能级（与 S7.1 的 _vac_window/_vac_level 同一算法：最大真空间隙中心 ± 宽度×half_ratio，平坦度检查）
# ---------------------------------------------------------------------------
def vacuum_window(zfrac, half_ratio=0.25):
    f = np.sort(np.asarray(zfrac, float) % 1.0)
    gaps = np.diff(f).tolist() + [f[0] + 1.0 - f[-1]]
    i = int(np.argmax(gaps))
    lo = f[i]
    hi = f[i + 1] if i < len(f) - 1 else f[0] + 1.0
    c, w = (lo + hi) / 2.0, (hi - lo) * half_ratio
    return (c - w) % 1.0, (c + w) % 1.0


def vacuum_level(planar, zfrac, half_ratio=0.25, flat_tol=1e-3, return_slope=False):
    """planar：沿真空轴的平面平均静电势（eV，n 个格点，第 j 个在分数坐标 j/n）。返回 (真空区中点的值, 去线性后的残差 std)，
    return_slope=True 时再加斜率（eV/分数坐标）。窗口以真空间隙中心为对称中心，可以跨 z = 1 边界。

    [0035-r3] 平坦度按"扣掉线性斜坡后的残差"判：周期性 slab 有面外偶极时（面外极性模的位移、Janus），真空区电势是一段
    直线（补偿电场），中点值 = 两侧真空能级的平均。对有 σh 的材料，±u 两个构型的中点值相同，面外模的 D 正确地为 0。
    残差 > flat_tol（取到 slab 尾巴、窗口里不是直线）-> ValueError。"""
    v = np.asarray(planar, float)
    n = len(v)
    a, b = vacuum_window(zfrac, half_ratio)
    z = np.arange(n) / float(n)
    sel = (z >= a) & (z <= b) if a <= b else (z >= a) | (z <= b)
    if sel.sum() < 3:
        raise ValueError("真空窗口里只有 %d 个格点（LOCPOT 网格太粗？）" % int(sel.sum()))
    c = (a + ((b - a) % 1.0) / 2.0) % 1.0                         # 窗口（真空间隙）中心
    dz = z[sel] - c
    dz -= np.round(dz)                                             # 以中心为原点、跨边界展开
    seg = v[sel]
    slope, mid = np.polyfit(dz, seg, 1)
    sd = float(np.std(seg - (mid + slope * dz)))
    if sd > flat_tol:
        raise ValueError("真空区平面平均势扣掉线性斜坡后仍不平坦 std=%.3g eV（> %.0e）：可能取到了 slab 尾巴"
                         % (sd, flat_tol))
    return (float(mid), sd, float(slope)) if return_slope else (float(mid), sd)


# ---------------------------------------------------------------------------
# 汇总：一个谷（星）、一个声子简并组
# ---------------------------------------------------------------------------
def combine_members(d2_by_mode):
    """d2_by_mode：{模 index: [各成员谷的 D²]}（同一简并组的各伙伴）。返回 (D_eV_A, 各成员 ΣD², 相对离散度)。"""
    arr = np.array([v for v in d2_by_mode.values()], float)          # (n_mode, n_member)
    per_member = arr.sum(axis=0)
    mean = float(per_member.mean())
    spread = float((per_member.max() - per_member.min()) / mean) if mean > 1e-12 else 0.0
    return math.sqrt(max(mean, 0.0)), [float(x) for x in per_member], spread


# ---------------------------------------------------------------------------
# S8.4：odp_params.json -> dp_optical 谷内通道
# ---------------------------------------------------------------------------
def _carrier(params, carrier):
    c = ((params or {}).get("carriers") or {}).get(carrier) or {}
    return c if c.get("status") == "ok" else None


def frozen_channels(params, carrier, keep=None, norm=None, dmin=DMIN_EV_A, alias=None):
    """odp_params.json 里某载流子的谷内参数 -> ([通道], [丢掉的], [提示])。
    keep：S8.4 谷表合并（同中心星）后保留的标签集合（None = 全要）；norm：标签归一函数（dp_optical.norm_label）；
    alias：被合并掉的标签 -> 保留的标签（dp_optical.merge_coincident 的第二个返回值）。
    同一个谷（含被并进来的星：同一组 k 上的另一条带或另一个自旋）的各星先把 ħω 相同（取到 0.001 meV）的组 D² 相加，
    再在星之间取平均（插件对该谷的所有带、所有自旋用同一组通道）；星之间相差 > 10% 写进提示。
    通道：{"from": L, "to": L, "kind": "intra", "D_eV_A", "hw_meV", "source": "frozen"}；D < dmin 的丢掉（对称性禁戒）。"""
    c = _carrier(params, carrier)
    if not c:
        return [], [], []
    norm = norm or (lambda s: s)
    alias = alias or {}
    per, hwv = {}, {}
    for v in c.get("valleys") or []:
        lab = norm(v.get("label"))
        tgt = alias.get(lab, lab)
        if keep is not None and tgt not in keep:
            continue
        src = per.setdefault(tgt, {}).setdefault(lab, {})
        for g in v.get("groups") or []:
            D, hw = float(g["D_eV_A"]), float(g["hw_meV"])
            if not (hw > 0 and math.isfinite(D)):
                continue
            key = round(hw, 3)
            src[key] = src.get(key, 0.0) + D * D
            hwv.setdefault(key, hw)
    chs, dropped, notes = [], [], []
    for tgt in sorted(per):
        srcs = per[tgt]
        for key in sorted(set().union(*[set(x) for x in srcs.values()])):
            vals = [srcs[s_].get(key, 0.0) for s_ in sorted(srcs)]
            d2 = float(np.mean(vals))
            D = math.sqrt(d2)
            if D < dmin:
                dropped.append({"valley": tgt, "hw_meV": hwv[key], "D_eV_A": D})
                continue
            if len(vals) > 1 and (max(vals) - min(vals)) / d2 > SPREAD_WARN:
                notes.append("谷 %s（并进来的星 %s）ħω=%.1f meV：各星 D 相差大（%s eV/Å），取 D² 平均"
                             % (tgt, ",".join(sorted(srcs)), hwv[key], "/".join("%.2f" % math.sqrt(x) for x in vals)))
            chs.append({"from": tgt, "to": tgt, "kind": "intra", "D_eV_A": D, "hw_meV": hwv[key], "source": "frozen"})
    return chs, dropped, notes


def params_digest(params, carrier):
    """某载流子参数的指纹（谷标签 + 各组 D、ħω，四舍五入）。S8.4 记下它；S8.21 重算后比对，变了才让 S8.4 重排。"""
    c = _carrier(params, carrier)
    if not c:
        return None
    rows = []
    for v in c.get("valleys") or []:
        for g in v.get("groups") or []:
            rows.append([str(v.get("label")), round(float(g["D_eV_A"]), 4), round(float(g["hw_meV"]), 3)])
    rows.sort()
    return hashlib.sha256(json.dumps(rows, ensure_ascii=False).encode("utf-8")).hexdigest()[:16]


def describe(params, carrier, norm=None, keep=None, dmin=DMIN_EV_A, alias=None):
    chs, _d, _n = frozen_channels(params, carrier, keep=keep, norm=norm, dmin=dmin, alias=alias)
    return ", ".join("%s %.2f eV/Å@%.1f meV" % (c["from"], c["D_eV_A"], c["hw_meV"]) for c in chs) or "无（全部 < %.2f eV/Å）" % dmin


# ---------------------------------------------------------------------------
# [0036] 三维谷表：S7 undeformed 网格上的带边极值点
# ---------------------------------------------------------------------------
VALLEY_WINDOW_EV = 0.30           # 只要离带边这么近的谷（900 K 时 e^(−0.3/kT) ≈ 2%）
STAR_TOL_EV = 0.002               # 同一个星：同带、同自旋、极值能量相差小于它（对称等价的谷能量相同）


def grid_from_kpoints(kpts):
    """ISYM=0 的 k 列表（分数坐标，只按 ±k 约化）-> (N (3,), 偏移 (3,), {整数网格点: (k 下标, ±1)})。
    网格由各方向的最小非零间隔推出；偏移 = k₀·N 的小数部分（Γ 心为 0，MP 偶网格为 ½）。不是规则网格 -> ValueError。"""
    k = np.asarray(kpts, float)
    N = []
    for a in range(3):
        x = np.mod(k[:, a], 1.0)
        x = np.unique(np.round(np.concatenate([x, np.mod(-x, 1.0)]), 6))
        d = np.diff(np.concatenate([x, [x[0] + 1.0]]))
        d = d[d > 1e-5]
        N.append(int(round(1.0 / float(d.min()))) if len(d) else 1)
    N = np.asarray(N, int)
    off = np.mod(k[0] * N, 1.0)
    off = np.where(np.abs(off - 1.0) < 1e-4, 0.0, off)
    table = {}
    for i, kk in enumerate(k):
        for sg in (1, -1):
            g = sg * kk * N - off
            gi = np.round(g)
            if float(np.max(np.abs(g - gi))) > 1e-3:
                raise ValueError("k 点 %s 不在推出的 %s 网格（偏移 %s）上：S7 的 k 列表不是规则网格"
                                 % (np.round(kk, 5).tolist(), N.tolist(), np.round(off, 3).tolist()))
            key = tuple(int(x) for x in np.mod(gi.astype(int), N))
            table.setdefault(key, (i, sg))
    if len(table) != int(np.prod(N)):
        raise ValueError("S7 的 k 点（含 −k）只覆盖 %d / %d 个网格点：不是 ISYM=0 的全网格（ISYM≥1 约化过？）"
                         % (len(table), int(np.prod(N))))
    return N, off, table


def _k_of(key, N, off):
    f = (np.asarray(key, float) + off) / N
    return f - np.round(f)                         # 取到 (−½, ½]


def find_valleys_3d(kpts, eig, nocc, sign, recip, window=VALLEY_WINDOW_EV, star_tol=STAR_TOL_EV):
    """三维带边谷表（S8.21 用；S8 的 ODP 通道用同一份谷中心）。
    kpts：S7 undeformed 的 k（ISYM=0）；eig：{自旋: (nk, nb) eV}；nocc：{自旋: 占据带数}；sign：+1 电子、−1 空穴。
    在每个自旋、本侧离带边 window 以内的每条带上找网格局部极值（26 个近邻里最小，载流子能量 x = sign·E），
    相邻的极值点（平台）只留一个；同一 k 上多条带/两个自旋（非磁 ISPIN=2）都有极值时只留最靠近带隙的那条（同一个谷：
    dp_optical 按谷中心归属 k 点、态密度统计全部带）；剩下的按（带、自旋、能量 ±star_tol）归星。
    返回 [{"label","g","k_frac","k_grid","k_S7","spin","band_index","dE_eV"}]，按 dE 升序；标签 Γ 或 V1、V2…"""
    N, off, table = grid_from_kpoints(kpts)
    keys = list(table)
    nbr = [(i, j, l) for i in (-1, 0, 1) for j in (-1, 0, 1) for l in (-1, 0, 1) if (i, j, l) != (0, 0, 0)]
    xs = {}
    edge = None
    for sp, e in eig.items():
        e = np.asarray(e, float)
        side = range(nocc[sp], e.shape[1]) if sign > 0 else range(nocc[sp] - 1, -1, -1)
        for b in side:
            col = sign * e[:, b]
            x = {key: float(col[table[key][0]]) for key in keys}
            xs[(sp, b)] = x
            m = min(x.values())
            edge = m if edge is None else min(edge, m)
    if edge is None:
        raise ValueError("没有%s带" % ("导" if sign > 0 else "价"))
    found = []
    for (sp, b), x in xs.items():
        if min(x.values()) > edge + window:
            continue
        seen = set()
        for key in sorted(keys, key=lambda q: x[q]):
            v = x[key]
            if v > edge + window:
                break
            ok = True
            for d in nbr:
                q = tuple(int(t) for t in np.mod(np.asarray(key) + d, N))
                if q != key and x[q] < v - 1e-6:
                    ok = False
                    break
            if not ok or key in seen:
                continue
            # 平台：与之相连、能量相同的网格点都算这一个极值
            stack = [key]
            while stack:
                c = stack.pop()
                if c in seen:
                    continue
                seen.add(c)
                for d in nbr:
                    q = tuple(int(t) for t in np.mod(np.asarray(c) + d, N))
                    if q not in seen and abs(x[q] - v) < 1e-6:
                        stack.append(q)
            found.append({"spin": sp, "band": b, "key": key, "x": v})
    # 同一 k：只留最靠近带隙的带（电子最低、空穴最高的带号）、第一个自旋
    best = {}
    for f in found:
        k0 = f["key"]
        cur = best.get(k0)
        rank = (f["band"] * sign, list(eig).index(f["spin"]))
        if cur is None or rank < (cur["band"] * sign, list(eig).index(cur["spin"])):
            best[k0] = f
    pts = sorted(best.values(), key=lambda f: (f["x"], f["spin"], f["band"]))
    stars = []
    for f in pts:
        for s_ in stars:
            if s_["spin"] == f["spin"] and s_["band"] == f["band"] and abs(s_["x"] - f["x"]) < star_tol:
                s_["keys"].append(f["key"])
                break
        else:
            stars.append({"spin": f["spin"], "band": f["band"], "x": f["x"], "keys": [f["key"]]})
    out, n = [], 0
    for s_ in stars:
        ks = [_k_of(q, N, off) for q in s_["keys"]]
        if len(ks) == 1 and np.all(np.abs(ks[0]) < 1e-9):
            lab = "Γ"
        else:
            n += 1
            lab = "V%d" % n
        kl = [np.round(k, 6).tolist() for k in ks]
        out.append({"label": lab, "g": len(ks), "k_frac": kl, "k_grid": kl, "k_S7": kl,
                    "spin": str(s_["spin"]), "band_index": int(s_["band"]), "dE_eV": round(s_["x"] - edge, 6)})
    return out


def valley_sig_3d(valleys):
    """三维谷表签名：每个星的中心集合（三个分量、模 1、4 位），与标签、顺序无关。S8 用它核对谷表同源。"""
    stars = set()
    for v in valleys or []:
        ks = []
        for k in v.get("k_frac") or []:
            f = []
            for x in list(k)[:3]:
                y = round(round(float(x), 4) % 1.0, 4)
                f.append(0.0 if y >= 1.0 else y)
            ks.append(tuple(f))
        stars.add(tuple(sorted(ks)))
    return [list(map(list, s_)) for s_ in sorted(stars)]
