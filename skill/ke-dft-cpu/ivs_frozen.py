# -*- coding: utf-8 -*-
"""ivs_frozen.py —— 冻声子法自动提取谷间零阶形变势（IVS，K↔K′）的纯 numpy 部分（[0050]）。

S7.3（gen_step9d_ivs_ph.py）建超胞、跑超胞 Γ̃ 点 DFPT；S7.4（gen_step9e_ivs.py）用它挑区边界模、生成 ±u 位移；
S8.22（gen_step12d_ivs_read.py）用它给超胞 Γ̃ 的能级标来源、从 K/K′ 二重态的劈裂算 D；S8.4 用它把 ivs_params.json
变成 dp_optical 的谷间通道（kind = inter_star，X>X'）。不依赖 pymatgen/amset，单元测试直接测它（tests/test_0050_ivs_frozen.py）。

思路（劈裂法，Piscanec 等 PRL 93, 185503 (2004) 求石墨烯 K 点耦合的同一类做法）：
  一个谷星的各成员（K、K′）在超胞里都折叠到 Γ̃（超胞矩阵 M 满足 M·k ∈ Z²），原胞的区边界声子（q = K）也折叠成超胞
  Γ̃ 点的模。冻结这样一个模（实的驻波），它没有 q = 0 分量，所以 K/K′ 二重态的对角元为 0，只有谷间非对角元 V_KK′：
  二重态劈裂成 ±|V|·u —— 一阶、线性、与参考能级无关（同一构型、同一 k 点两条能级的差）。

归一化（与 odp_frozen / dp_optical 同一口径，D = dE/du 直接能用）：
  超胞当原胞用：Δ_j = u·√(M_s/m_j)·ẽ_j（ẽ 为超胞质量加权、归一的本征矢，M_s = 超胞质量 = det·M），
  |g|² = ħ D_s²/(2 N_s M_s ω) = ħ D_s²/(2 N M ω)（N_s M_s = N M）—— 超胞口径的 D 就是原胞口径的 D。
  实驻波 C/S 与复声子 A（q 连 K→K′）的关系：⟨K′|V_C|K⟩ = ⟨K′|V_A|K⟩/√2，⟨K′|V_S|K⟩ = ⟨K′|V_A|K⟩/(i√2)，
  所以一对简并驻波的 |V|² 之和 = |V_A|²。每个模的 Σλ²（迹，组内去掉共同平移）= Tr(V²) = Σ_{i≠j}|V_ij|²；
  dp_optical 的 inter_star 率 = πD²/(Mω)·G_星·(g−1)/g，等价于每个成员对其它 g−1 个成员各 |V|²，
  故  D² = Σ_{简并伙伴模} Σλ² / (g(g−1))，g = 星的成员数（K 星 g = 2）。
  三角晶格紧束缚玩具模型（tests）：给定在位耦合 α，按这套流程算出的 D 等于 α。

[0050] 只做 A 方案：非 SOC、共线能带（ISPIN = 1 或 2 都行，按自旋道分别算）。SOC 能带下 K/K′ 二重态被 Δso 劈开、
  自旋守恒的谷间散射只能落到同自旋子带（Xiao 等 PRL 108, 196802 (2012)），劈裂法与 dp_optical 的星平均态密度都不适用：
  gen 与 S8.4 遇到 SOC 直接拒绝。只算零阶（q 正好在 K），一阶谷间耦合（∝ |q − K|）与 K↔Q、K↔Γ 都不在内。
"""
import hashlib
import json
import math

import numpy as np

from odp_frozen import (FREQ_DEG_THZ, THZ_TO_MEV, BandCrossing, _sorted_sq, check_eigvecs,  # noqa: F401
                        displacement)

SCHEMA = "ivs_params/1"
FOLD_TOL = 1e-4                   # M·k 离整数的容差（分数坐标）
MAX_DET = 9                       # 超胞最大原胞数（K 只要 3；再大单点太贵，报错）
DEG_TOL_EV = 0.005                # u = 0 时谷星各成员的同一条带必须简并到这个程度（非 SOC、时间反演）
GAP_MARGIN_EV = 0.002             # 组外最近的能级离得太近（< 2·最大移动 + 它）就做交叉核对
CROSS_REL = 0.05
CROSS_ABS = 0.05
ZB_WEIGHT = 0.5                   # 本征矢在"原胞平移不变"子空间上的投影权重 < 它 = 区边界（折叠来的）模
LABEL_TOL_EV = 0.05               # 超胞 Γ̃ 能级与原胞折叠能级（对齐后）的最大偏差
LABEL_WINDOW_EV = 1.0             # 标来源只看带边上下这么宽（更远的带 NBANDS 截断、收敛差，不用）
AMBIG_EV = 0.01                   # 目标二重态附近这么近有别的来源的能级 -> 标不清，报错
FMIN_THZ = 0.5                    # 区边界模频率低于它不算（虚频/太软）
DMIN_EV_A = 0.05                  # S8.4 用：D 小于它的通道不进 SCATTERING（对称性禁戒）


# ---------------------------------------------------------------------------
# 超胞：让谷星的所有成员折叠到 Γ̃
# ---------------------------------------------------------------------------
def fold_matrix(members_frac, lat2=None, max_det=MAX_DET, tol=FOLD_TOL):
    """谷星各成员 k（原胞倒格分数坐标，取前两个分量）-> 面内超胞矩阵 M（2×2 整数，行 = 超胞面内基矢在原胞基矢下的系数）。
    条件：M·k ∈ Z²（k 是超胞的倒格矢 = 折叠到 Γ̃）。取 |det| 最小的；同 det 里按 lat2（原胞面内 2×2 笛卡尔基矢）
    取最紧凑的（两基矢长度平方和最小，再看夹角接近 60°/120°）。det 取正（右手）。找不到 -> ValueError。"""
    ks = np.asarray(members_frac, float)
    if ks.ndim != 2 or ks.shape[1] < 2:
        raise ValueError("谷星成员 k 的形状 %s 不对" % (ks.shape,))
    if ks.shape[1] >= 3 and np.any(np.abs(ks[:, 2] - np.rint(ks[:, 2])) > tol):
        raise ValueError("谷星成员 k 的面外分量不为 0（%s）：本方法只处理面内折叠" % ks[:, 2].tolist())
    k2 = ks[:, :2]
    L = np.eye(2) if lat2 is None else np.asarray(lat2, float)[:2, :2]
    rng = range(-max_det, max_det + 1)
    rows = [np.array([a, b]) for a in rng for b in rng if (a, b) != (0, 0)]
    ok = [r for r in rows if np.all(np.abs((k2 @ r) - np.rint(k2 @ r)) < tol)]
    best = None
    for i, r1 in enumerate(ok):
        for r2 in ok[i + 1:]:
            d = int(round(r1[0] * r2[1] - r1[1] * r2[0]))
            if d == 0 or abs(d) > max_det:
                continue
            m = np.array([r1, r2]) if d > 0 else np.array([r2, r1])
            v = m @ L
            l2 = float(np.sum(v * v))
            cosang = abs(float(v[0] @ v[1]) / math.sqrt(float(v[0] @ v[0]) * float(v[1] @ v[1])))
            key = (abs(d), round(l2, 6), round(abs(cosang - 0.5), 6), -int(m.sum()))   # 同样紧凑时取系数为正的
            if best is None or key < best[0]:
                best = (key, m)
    if best is None:
        raise ValueError("找不到 |det| ≤ %d 的超胞让谷星 %s 折叠到 Γ̃（谷不在可公度的点上，如 Q 谷）"
                         % (max_det, np.round(k2, 4).tolist()))
    return best[1].astype(int)


def folded_kpoints(M2):
    """M（2×2）下折叠到 Γ̃ 的全部原胞 k（分数坐标，模 1，第三分量 0），共 det 个；第一个是 Γ。"""
    M2 = np.asarray(M2, int)
    det = int(round(np.linalg.det(M2)))
    inv = np.linalg.inv(M2.astype(float))
    pts = []
    for a in range(det):
        for b in range(det):
            k = inv @ np.array([a, b], float)          # M·k = (a, b) ∈ Z²
            k = k - np.floor(k + 1e-9)
            if not any(np.allclose(k, p[:2], atol=1e-6) for p in pts):
                pts.append(np.array([k[0], k[1], 0.0]))
    pts.sort(key=lambda p: (round(p[0], 6), round(p[1], 6)))
    if len(pts) != det:
        raise ValueError("折叠 k 点数 %d ≠ det %d" % (len(pts), det))
    return pts


def build_supercell(lattice, frac, species, M2):
    """原胞（lattice 3×3 行矢量、frac (n,3)、species 列表）-> 超胞（面内 M，面外不变）。
    返回 dict：lattice、frac (n·det, 3)、species、origin [(原胞原子号, 平移 t (2,) 整数)]、translations [t…]、M3。
    原子按"原胞原子号优先、再平移"排序：原胞同元素连续，超胞也连续（POSCAR 要求）。"""
    M2 = np.asarray(M2, int)
    det = int(round(np.linalg.det(M2)))
    if det <= 0:
        raise ValueError("超胞矩阵 det = %d（要正）" % det)
    M3 = np.eye(3)
    M3[:2, :2] = M2
    lat = np.asarray(lattice, float)
    slat = M3 @ lat
    inv3 = np.linalg.inv(M3)
    f0 = np.asarray(frac, float)
    span = range(-det - 1, det + 2)
    trans = []
    for a in span:                               # 不等价的原胞平移（模超胞）：t·M⁻¹ 在 [0,1)² 里的那一个代表
        for b in span:
            t = np.array([a, b, 0.0])
            ft = t @ inv3
            if np.all(ft[:2] > -1e-9) and np.all(ft[:2] < 1 - 1e-9):
                trans.append(np.array([a, b], int))
    if len(trans) != det:
        raise ValueError("超胞里找到 %d 个原胞平移（应为 %d）" % (len(trans), det))
    trans.sort(key=lambda t: (abs(int(t[0])) + abs(int(t[1])), int(t[0]), int(t[1])))   # (0,0) 排第一
    sfrac, sspec, origin = [], [], []
    for i, f in enumerate(f0):
        for t in trans:
            fs = (f + np.array([t[0], t[1], 0.0])) @ inv3
            fs[:2] -= np.floor(fs[:2] + 1e-12)
            sfrac.append(fs)
            sspec.append(species[i])
            origin.append((i, t.copy()))
    return {"lattice": slat, "frac": np.array(sfrac), "species": list(sspec), "origin": origin,
            "translations": trans, "M3": M3, "det": det}


def translation_perms(sc):
    """超胞里每个原胞平移 T 对应的原子置换 perm_T：原子 j 平移 T 后落在 perm_T[j]。返回 [perm…]（第一个是恒等）。"""
    inv3 = np.linalg.inv(sc["M3"])
    f = np.asarray(sc["frac"], float)
    perms = []
    for T in sc["translations"]:
        shift = np.array([T[0], T[1], 0.0]) @ inv3
        g = f + shift
        p = []
        for j in range(len(f)):
            d = f - g[j]
            d[:, :2] -= np.rint(d[:, :2])
            idx = np.where(np.all(np.abs(d) < 1e-6, axis=1))[0]
            if len(idx) != 1 or sc["species"][idx[0]] != sc["species"][j]:
                raise ValueError("原胞平移 %s 下原子 %d 找不到唯一对应" % (T.tolist(), j))
            p.append(int(idx[0]))
        perms.append(np.array(p))
    return perms


def gamma_weight(eigvecs, perms):
    """超胞 Γ̃ 模在"原胞平移不变"子空间上的投影权重：原胞 Γ 来的模 ≈ 1，区边界（K 等）折叠来的模 ≈ 0。
    eigvecs：(nmode, nat, 3) 质量加权、归一（同元素原子质量相同，平移不改质量权重）。"""
    e = np.asarray(eigvecs, float)
    proj = np.zeros_like(e)
    for p in perms:
        moved = np.zeros_like(e)
        moved[:, p, :] = e                           # 原子 j 的位移挪到 perm[j]
        proj += moved
    proj /= float(len(perms))
    return np.einsum("nka,nka->n", proj, proj)


def zone_boundary_groups(freqs_thz, eigvecs, perms, fmin_thz=FMIN_THZ, deg_thz=FREQ_DEG_THZ):
    """超胞 Γ̃ 模 -> (modes, groups, notes)。只留区边界模（gamma_weight < ZB_WEIGHT），按频率升序编号、简并分组。
    modes[i]：index（全体模按频率升序从 1 数）、freq_THz、hw_meV、gamma_weight、zone_boundary、group、_src。
    groups[j]：group、modes、freq_THz、hw_meV、use、why。"""
    f = np.asarray(freqs_thz, float)
    e = np.asarray(eigvecs, float)
    check_eigvecs(e)
    w = gamma_weight(e, perms)
    order = np.argsort(f, kind="stable")
    det = len(perms)
    nat = e.shape[1]
    notes = []
    modes = []
    for rank, i in enumerate(order):
        modes.append({"index": rank + 1, "freq_THz": float(f[i]), "hw_meV": float(f[i] * THZ_TO_MEV),
                      "gamma_weight": float(w[i]), "zone_boundary": bool(w[i] < ZB_WEIGHT), "group": None,
                      "_src": int(i)})
    nzb = sum(m["zone_boundary"] for m in modes)
    want = 3 * nat - 3 * nat // det
    if nzb != want:
        notes.append("区边界模数 %d ≠ 期望 %d（3·N_sc − 3·N_prim）：平移投影 %s"
                     % (nzb, want, ", ".join("%.2f" % m["gamma_weight"] for m in modes)))
    groups = []
    for md in modes:
        if not md["zone_boundary"]:
            continue
        if groups and abs(md["freq_THz"] - groups[-1]["_flast"]) < deg_thz:
            g = groups[-1]
            g["modes"].append(md["index"])
            g["_f"].append(md["freq_THz"])
            g["_flast"] = md["freq_THz"]
        else:
            groups.append({"group": len(groups) + 1, "modes": [md["index"]], "_f": [md["freq_THz"]],
                           "_flast": md["freq_THz"]})
        md["group"] = groups[-1]["group"]
    for g in groups:
        fm = float(np.mean(g.pop("_f")))
        g.pop("_flast")
        g["freq_THz"] = fm
        g["hw_meV"] = fm * THZ_TO_MEV
        why = None
        if fm <= fmin_thz:
            why = "频率 %.3f THz ≤ %.2f（虚频或太软）" % (fm, fmin_thz)
        g["use"] = why is None
        g["why"] = why
    return modes, groups, notes


# ---------------------------------------------------------------------------
# 超胞 Γ̃ 能级标来源（原胞折叠来的哪个 k、哪条带）
# ---------------------------------------------------------------------------
def label_folded(e_sc, prim, edge_eV, window=LABEL_WINDOW_EV, tol=LABEL_TOL_EV):
    """e_sc：超胞参考构型 Γ̃ 点某自旋的本征值（eV，升序）；prim：{k 名: 原胞同一自旋在该折叠 k 的本征值（eV，升序）}，
    覆盖全部折叠 k；edge_eV：原胞口径的带边参考能量（只标它上下 window 以内的能级）。
    原胞各折叠 k 的能级并起来排序 = 超胞 Γ̃ 的能级（差一个常数平移 + 收敛误差）：按序配对、取中位数平移对齐，
    残差 > tol -> ValueError（两边输入不一致，或 NBANDS 截断进了窗口）。
    返回 (labels, shift)：labels[i] = (k 名, 原胞带号) 或 None（窗口外），shift = e_sc − e_prim 的中位数。"""
    e_sc = np.sort(np.asarray(e_sc, float))
    pred = sorted((float(x), kn, b) for kn, arr in prim.items() for b, x in enumerate(np.asarray(arr, float)))
    n = min(len(e_sc), len(pred))
    pe = np.array([p[0] for p in pred[:n]])
    lo, hi = edge_eV - window, edge_eV + window
    top = min(float(np.max(np.asarray(arr, float))) for arr in prim.values())
    if hi >= top:                                   # 某个折叠 k 的能带在窗口里就截断了：并起来的能级不全，按序配对会错位
        raise ValueError("窗口顶端 %.3f eV 高于原胞某折叠 k 的最高能带 %.3f eV（NBANDS 截断）：加大 NBANDS"
                         % (hi, top))
    sel = (pe >= lo) & (pe <= hi)
    if sel.sum() < 2:
        raise ValueError("带边 %.3f eV 上下 %.1f eV 内的折叠能级不到 2 个（NBANDS 太少？）" % (edge_eV, window))
    if sel[n - 1]:                                  # 最高的那条（NBANDS 截断处，Davidson 收敛最差）进了窗口
        raise ValueError("窗口顶端碰到了参与比对的最高能带（NBANDS 截断）：加大 NBANDS 或减小窗口")
    shift = float(np.median(e_sc[:n][sel] - pe[sel]))
    res = np.abs(e_sc[:n][sel] - pe[sel] - shift)
    if float(res.max()) > tol:
        raise ValueError("超胞 Γ̃ 能级与原胞折叠能级对不上（窗口内最大偏差 %.1f meV > %.0f）：两边输入/结构不一致，"
                         "或折叠 k 点不全" % (1000 * float(res.max()), 1000 * tol))
    labels = [None] * len(e_sc)
    for i in range(n):
        if sel[i]:
            labels[i] = (pred[i][1], pred[i][2])
    return labels, shift


def edge_group(labels, e_sc, star_names, band, ambig=AMBIG_EV, deg_tol=DEG_TOL_EV):
    """谷星成员（star_names）的原胞带 band 在超胞 Γ̃ 的能级下标（升序）。核对：每个成员正好一个、彼此简并（≤ deg_tol）、
    附近 ambig 以内没有别的来源的能级（否则按序配对可能把来源配错）。不满足 -> ValueError。"""
    e_sc = np.sort(np.asarray(e_sc, float))
    idx = [i for i, lab in enumerate(labels) if lab is not None and lab[0] in star_names and lab[1] == band]
    if len(idx) != len(star_names):
        raise ValueError("谷星 %s 的带 %d 在超胞 Γ̃ 找到 %d 个能级（应为 %d）" % (list(star_names), band, len(idx),
                                                                  len(star_names)))
    ee = e_sc[idx]
    if float(ee.max() - ee.min()) > deg_tol:
        raise ValueError("谷星 %s 的带 %d 在超胞里不简并（相差 %.1f meV > %.0f）：SOC/磁性破坏了 K/K′ 简并？"
                         "本方法只支持非 SOC 共线能带" % (list(star_names), band, 1000 * float(ee.max() - ee.min()),
                                                  1000 * deg_tol))
    c = float(ee.mean())
    for i, lab in enumerate(labels):
        if i in idx or lab is None:
            continue
        if abs(e_sc[i] - c) < ambig:
            raise ValueError("谷星 %s 的带 %d 附近 %.0f meV 内有别的来源的能级（%s 带 %d）：标不清，劈裂法不可靠"
                             % (list(star_names), band, 1000 * ambig, lab[0], lab[1]))
    return idx


# ---------------------------------------------------------------------------
# 劈裂 -> Σλ²
# ---------------------------------------------------------------------------
def _centered_sq(e):
    e = np.asarray(e, float)
    return float(np.sum((e - e.mean()) ** 2))


def split_sq(e0, ep, em, grp, u, deg_tol=DEG_TOL_EV, margin=GAP_MARGIN_EV):
    """e0/ep/em：u = 0、+u、−u 三个构型同一自旋 Γ̃ 点全部能级（eV，升序）；grp：谷星二重（多重）态的下标（u = 0 时）。
    区边界模没有 q = 0 分量：组内对角元为 0，劈裂只来自谷间非对角元。每个构型先减去组内平均（去掉共同平移，
    与参考能级、真空无关），再用迹：Σλ² = [Σ(ep−⟨ep⟩)² + Σ(em−⟨em⟩)² − 2Σ(e0−⟨e0⟩)²] / 2u²。
    组外能级离得比 2×最大移动 + margin 还近 -> 连同它们排序配对再算一遍，与迹不一致 = 交叉/强混合 -> BandCrossing。
    返回 dict：S2、g、shift_max_eV、gap_out_eV、widened、near。"""
    e0 = np.asarray(e0, float)
    ep = np.asarray(ep, float)
    em = np.asarray(em, float)
    nb = len(e0)
    if len(ep) != nb or len(em) != nb:
        raise ValueError("三个构型的能级数不同（%d/%d/%d）" % (nb, len(ep), len(em)))
    grp = sorted(int(i) for i in grp)
    lo, hi = grp[0], grp[-1]
    if hi - lo + 1 != len(grp):
        raise ValueError("谷星能级 %s 在 u = 0 时不连续：中间夹着别的能级" % grp)
    if float(e0[hi] - e0[lo]) > deg_tol:
        raise ValueError("谷星能级 %s 在 u = 0 时不简并（%.1f meV）" % (grp, 1000 * float(e0[hi] - e0[lo])))
    g = len(grp)
    s = slice(lo, hi + 1)
    S2 = (_centered_sq(ep[s]) + _centered_sq(em[s]) - 2.0 * _centered_sq(e0[s])) / (2.0 * u * u)
    shift = float(max(np.max(np.abs(ep[s] - e0[s])), np.max(np.abs(em[s] - e0[s]))))
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
    out = {"S2": float(S2), "g": g, "shift_max_eV": shift, "gap_out_eV": min(gaps) if gaps else None,
           "widened": (wlo, whi) != (lo, hi), "near": [wlo, whi]}
    if out["widened"]:
        w = slice(wlo, whi + 1)
        c0 = e0[w] - e0[s].mean()
        cp = ep[w] - ep[s].mean()
        cm = em[w] - em[s].mean()
        s_w = _sorted_sq(c0, cp, cm, u)
        t_w = (float(np.sum(cp ** 2)) + float(np.sum(cm ** 2)) - 2.0 * float(np.sum(c0 ** 2))) / (2.0 * u * u)
        out["near_check"] = {"sorted": s_w, "trace": t_w}
        if abs(s_w - t_w) > CROSS_REL * max(abs(s_w), abs(t_w)) + CROSS_ABS:
            raise BandCrossing(
                "谷星能级 %s 附近的 %d–%d 在 ±u 构型里交叉或强混合（排序配对 %.3g，迹 %.3g (eV/Å)²；移动 %.0f meV，"
                "邻近间隔 %.0f meV）：这个位移下劈裂解不出来" % (grp, wlo, whi, s_w, t_w, 1000 * shift,
                                                     1000 * min(gaps) if gaps else 0.0))
    return out


def group_D(s2_by_mode, g):
    """同一简并组各伙伴模的 Σλ² -> D（eV/Å）：D² = Σ_模 Σλ² / (g(g−1))。g = 谷星成员数（≥ 2）。"""
    if g < 2:
        raise ValueError("谷星只有 %d 个成员：没有谷间" % g)
    tot = float(sum(s2_by_mode.values()))
    return math.sqrt(max(tot, 0.0) / (g * (g - 1)))


# ---------------------------------------------------------------------------
# NBANDS：超胞与原胞要有足够的空带给"标来源"用（两边按同一规则取，折叠后数目对得上）
# ---------------------------------------------------------------------------
def nelect_from_potcar(potcar_text, counts):
    """POTCAR 文本（各元素依次拼接）+ 各元素原子数（POSCAR 顺序）-> 总价电子数 NELECT。"""
    import re
    zv = [float(x) for x in re.findall(r"ZVAL\s*=\s*([-+]?\d+(?:\.\d*)?)", potcar_text)]
    if len(zv) != len(counts):
        raise ValueError("POTCAR 里 ZVAL %d 个，POSCAR 元素 %d 种：对不上" % (len(zv), len(counts)))
    return float(sum(z * n for z, n in zip(zv, counts)))


def nbands_for(nelect, ncell, extra_per_cell):
    """NBANDS = ⌈0.6·NELECT⌉ + extra·ncell（ncell = 原胞数）。原胞与超胞用同一规则：折叠后空带数相当。"""
    return int(math.ceil(0.6 * float(nelect))) + int(extra_per_cell) * int(ncell)


# ---------------------------------------------------------------------------
# S8.4：ivs_params.json -> dp_optical 谷间通道
# ---------------------------------------------------------------------------
def _carrier(params, carrier):
    c = ((params or {}).get("carriers") or {}).get(carrier) or {}
    return c if c.get("status") == "ok" else None


def inter_channels(params, carrier, keep=None, norm=None, dmin=DMIN_EV_A, alias=None):
    """ivs_params.json 里某载流子的谷间参数 -> ([通道], [丢掉的], [提示])。同一个谷的各组按 ħω（0.001 meV）把 D² 相加
    （同中心并进来的星取平均，同 odp_frozen.frozen_channels）。通道：{"from": L, "to": L, "kind": "inter_star",
    "D_eV_A", "hw_meV", "source": "frozen_ivs"}（dp_optical 写成 L>L'）；D < dmin 的丢掉（对称性禁戒）。"""
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
        for gr in v.get("groups") or []:
            D, hw = float(gr["D_eV_A"]), float(gr["hw_meV"])
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
            D = math.sqrt(float(np.mean(vals)))
            if D < dmin:
                dropped.append({"valley": tgt, "hw_meV": hwv[key], "D_eV_A": D})
                continue
            chs.append({"from": tgt, "to": tgt, "kind": "inter_star", "D_eV_A": D, "hw_meV": hwv[key],
                        "source": "frozen_ivs"})
    return chs, dropped, notes


def params_digest(params, carrier):
    c = _carrier(params, carrier)
    if not c:
        return None
    rows = []
    for v in c.get("valleys") or []:
        for gr in v.get("groups") or []:
            rows.append([str(v.get("label")), round(float(gr["D_eV_A"]), 4), round(float(gr["hw_meV"]), 3)])
    rows.sort()
    return hashlib.sha256(json.dumps(rows, ensure_ascii=False).encode("utf-8")).hexdigest()[:16]


def describe(params, carrier, norm=None, keep=None, dmin=DMIN_EV_A, alias=None):
    chs, _d, _n = inter_channels(params, carrier, keep=keep, norm=norm, dmin=dmin, alias=alias)
    return ", ".join("%s↔%s′ %.2f eV/Å@%.1f meV" % (c["from"], c["from"], c["D_eV_A"], c["hw_meV"]) for c in chs) \
        or "无（全部 < %.2f eV/Å）" % dmin
