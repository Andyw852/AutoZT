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
        中心差分会把 +λ 和 −λ 配错对（得到 0），所以改用单侧斜率的平方和（两侧平均，线性误差项相消）：
            Σλ² ≈ ½[Σ_i ((E_i(+u)−E₀)/u)² + Σ_i ((E_i(−u)−E₀)/(−u))²]
        dp_optical 的目标谷态密度 G 统计该谷的全部能带（d 条都算），每个态的散射率是 Σ_j|M_ij|²·G/d 的平均，
        即 D_eff² = Tr(M²)/d² = Σλ²/d²（d = 1 时退回 λ²）。
    同一组简并声子（ħω 相同）的各个伙伴模的 D² 相加 —— 和与简并子空间里本征矢的取法无关（正交变换不变）。
    星的各成员谷：对称等价，Σ_伙伴 D² 应相同；取平均并记录离散度。

极性模（红外活性）默认不算：它们的长程（Fröhlich）耦合由 POP 管；短程部分要算就设 ODP_POLAR = include。
"""
import hashlib
import json
import math

import numpy as np

THZ_TO_MEV = 4.135667696          # h·1 THz = 4.1357 meV
DEG_TOL_EV = 0.005                # 带简并判据：u = 0 时同一 k 点、同一自旋两条带相差小于它
GAP_MARGIN_EV = 0.002             # 简并组外最近的带离得太近（< 2·最大位移 + 它）就并进来，免得排序换位
FREQ_DEG_THZ = 0.05               # 声子简并判据（ħω 相差 < 0.2 meV）
ACOUSTIC_WEIGHT = 0.5             # 本征矢在平移子空间上的投影权重 > 它 = 声学模
POLAR_TOL_E = 0.1                 # 模有效电荷 |Σ Z*·U|/|U|（e）> 它 = 极性模
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
               deg_thz=FREQ_DEG_THZ, include_polar=False):
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
                      "polar": bool(zmode is not None and zmode[i] > polar_tol),
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
        elif g["polar"] is None:
            why = "没有 Born 有效电荷，判不了极性"
        elif g["polar"] and not include_polar:
            why = "极性模（模有效电荷 %.2f e > %.2f）：长程耦合由 POP 管（ODP_POLAR=include 才算）" % (
                max(by_index[i]["Z_mode_e"] for i in g["modes"]), polar_tol)
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
def edge_coupling(e0, ep, em, band, u, dvac_p=0.0, dvac_m=0.0, deg_tol=DEG_TOL_EV, margin=GAP_MARGIN_EV):
    """e0/ep/em：该 k 点、该自旋全部能带的本征值（eV，升序）——u = 0、+u、−u 三个构型；
    dvac_p/dvac_m：+u/−u 构型的真空能级相对 u = 0 的变化（eV），从能量里扣掉（真空对齐）。
    返回 dict：D2（(eV/Å)²，按模块说明的 D_eff²）、d（简并组带数）、bands [lo, hi]、method、
              lam（一阶斜率：非简并为中心差分值；简并为 +u 侧排序斜率）、shift_max_eV、gap_out_eV、widened。"""
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
    widened = False
    for _ in range(nb):
        idx = slice(lo, hi + 1)
        shift = float(max(np.max(np.abs(ep[idx] - e0[idx])), np.max(np.abs(em[idx] - e0[idx]))))
        need = 2.0 * shift + margin
        grow = False
        if lo > 0 and e0[lo] - e0[lo - 1] < need:
            lo -= 1
            grow = True
        if hi < nb - 1 and e0[hi + 1] - e0[hi] < need:
            hi += 1
            grow = True
        if not grow:
            break
        widened = True
    idx = slice(lo, hi + 1)
    d = hi - lo + 1
    gaps = []
    if lo > 0:
        gaps.append(float(e0[lo] - e0[lo - 1]))
    if hi < nb - 1:
        gaps.append(float(e0[hi + 1] - e0[hi]))
    shift = float(max(np.max(np.abs(ep[idx] - e0[idx])), np.max(np.abs(em[idx] - e0[idx]))))
    out = {"d": d, "bands": [lo, hi], "shift_max_eV": shift, "gap_out_eV": min(gaps) if gaps else None,
           "widened": widened}
    if d == 1:
        lam = (ep[b] - em[b]) / (2.0 * u)
        out.update(D2=float(lam * lam), method="central", lam=[float(lam)])
        return out
    sp = (np.sort(ep[idx]) - np.sort(e0[idx])) / u
    sm = (np.sort(em[idx]) - np.sort(e0[idx])) / (-u)
    s2 = 0.5 * (float(np.sum(sp * sp)) + float(np.sum(sm * sm)))
    out.update(D2=s2 / (d * d), method="degenerate(%d)" % d, lam=[float(x) for x in sp],
               trace_pm=[float(np.sum(sp)), float(np.sum(sm))])
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


def vacuum_level(planar, zfrac, half_ratio=0.25, flat_tol=1e-3):
    """planar：沿真空轴的平面平均静电势（eV，n 个格点，第 j 个在分数坐标 j/n）。返回 (平均, 标准差)。
    窗口可以跨 z = 1 边界。不平坦（std > flat_tol）-> ValueError（偶极层/取到 slab 尾巴，对齐不可信）。"""
    v = np.asarray(planar, float)
    n = len(v)
    a, b = vacuum_window(zfrac, half_ratio)
    z = np.arange(n) / float(n)
    sel = (z >= a) & (z <= b) if a <= b else (z >= a) | (z <= b)
    if sel.sum() < 2:
        raise ValueError("真空窗口里只有 %d 个格点（LOCPOT 网格太粗？）" % int(sel.sum()))
    seg = v[sel]
    sd = float(np.std(seg))
    if sd > flat_tol:
        raise ValueError("真空区平面平均势不平坦 std=%.3g eV（> %.0e）：可能有面外偶极（极性面外模）或取到了 slab 尾巴"
                         % (sd, flat_tol))
    return float(np.mean(seg)), sd


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


def frozen_channels(params, carrier, keep=None, norm=None, dmin=DMIN_EV_A):
    """odp_params.json 里某载流子的谷内参数 -> ([通道], [丢掉的])。
    keep：S8.4 谷表合并（同中心星）后保留的标签集合（None = 全要）；norm：标签归一函数（dp_optical.norm_label）。
    通道：{"from": L, "to": L, "kind": "intra", "D_eV_A", "hw_meV", "source": "frozen"}；同一谷、ħω 相同（四舍五入到
    0.001 meV）的几组 D² 相加（dp_optical 不许重复通道）。D < dmin 的丢掉（对称性禁戒，记进第二个返回值）。"""
    c = _carrier(params, carrier)
    if not c:
        return [], []
    norm = norm or (lambda s: s)
    acc, dropped = {}, []
    seen = set()
    for v in c.get("valleys") or []:
        lab = norm(v.get("label"))
        if keep is not None and lab not in keep:
            continue
        if lab in seen:            # 同一个标签出现两次（S8.2 同中心的两个星归一后同名）：参数相同，只取第一份
            continue
        seen.add(lab)
        for g in v.get("groups") or []:
            D, hw = float(g["D_eV_A"]), float(g["hw_meV"])
            if not (D >= dmin and hw > 0 and math.isfinite(D)):
                dropped.append({"valley": lab, "group": g.get("group"), "D_eV_A": D, "hw_meV": hw})
                continue
            key = (lab, round(hw, 3))
            acc[key] = acc.get(key, 0.0) + D * D
    chs = [{"from": lab, "to": lab, "kind": "intra", "D_eV_A": math.sqrt(d2), "hw_meV": hw, "source": "frozen"}
           for (lab, hw), d2 in sorted(acc.items(), key=lambda t: (t[0][0], t[0][1]))]
    return chs, dropped


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


def describe(params, carrier, norm=None, keep=None, dmin=DMIN_EV_A):
    chs, _ = frozen_channels(params, carrier, keep=keep, norm=norm, dmin=dmin)
    return ", ".join("%s %.2f eV/Å@%.1f meV" % (c["from"], c["D_eV_A"], c["hw_meV"]) for c in chs) or "无（全部 < %.2f eV/Å）" % dmin
