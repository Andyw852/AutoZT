# -*- coding: utf-8 -*-
"""dp_optical.py —— 光学形变势（ODP，谷内）与谷间散射（IVS）的常数形变势模型（[0033]）。

纯 numpy，不依赖 amset：amset2d_plugin.py 在 AMSET 运行里调用它给每个 k 点算散射率
（AMSET 的 AbstractBasicScattering 接口，单位 1/s）；gen_step14 用它解析 step.conf；单元测试直接测它。

模型（Jin et al. PRB 90, 045422 (2014)；Kaasbjerg et al. PRB 85, 115317 (2012) 的简化形式）：
    矩阵元与声子波矢 q 无关：|M|² = ħ D² / (2 ρ V ω) · (N + ½ ∓ ½)
    费米黄金定则对目标谷求和后，谷 i 里能量 x 的态的散射率（每个通道 c = i→j）：
        1/τ_c(x) = π D_c² / (M ω_c) · [ N_c · G_c(x + ħω_c) + (N_c + 1) · G_c(x − ħω_c) ]
    M = 原胞质量（kg）；G_c = 目标谷集合的态密度（每原胞、每自旋、每焦耳）；D_c 单位 J/m；
    N_c = 1/(e^{ħω/kT} − 1)。体积（含真空）在 D²·G/M 里抵消，与真空长度无关。
    x 是"载流子能量"：电子 x = E，空穴 x = −E —— 吸收声子 x 升高，发射声子 x 降低。
    发射/吸收的能量阈值：目标谷在 x ± ħω 处没有态，G = 0（G 做了 6 meV 高斯展宽，阈值附近 ±2σ 内被抹圆）。
    矩阵元与 q 无关 => 各向同性散射，动量弛豫率 = 散射率。
    [0042] 简并统计（泡利阻塞）：传入 fermi_levels（每个掺杂、温度的电子费米能级）时，按线性化出射率
        1/τ_c(x) = π D_c²/(M ω_c) · [ (N_c + f(x + ħω_c))·G_c(x + ħω_c) + (N_c + 1 − f(x − ħω_c))·G_c(x − ħω_c) ]
    f 是载流子占据（电子 f_e；空穴 1 − f_e）。它与 Kaasbjerg Eq.30 的 Σ P_kk′(1 − f_k′)/(1 − f_k) 恒等，
    也与 AMSET 自己的 POP（N + 1 − f / N + f）同一口径；不传 fermi_levels 时退回非简并极限（Kaasbjerg Eq.31）。
    D 取冻声子 Γ 点（q = 0）的裸耦合时是速率的上限：真实 |M_q| 离开 Γ 点会下降（Kaasbjerg MoS₂ 同极模：
    Γ 点 ≈4.8 eV/Å，拟合散射率的等效常数 4.1 eV/Å，速率差 ≈1.37 倍）。

谷的约定（与 S8.2 multi_valley 的谷表一致）：
    星（star）= 一组对称等价的谷，如 K 星 2 个谷（K、K′）、Q 星 6 个、Γ 星 1 个。每个成员谷有自己的 k 中心。
    通道写法：X>X   谷内（同一个成员谷，ODP）
              X>X'  同一星的其它成员（如 K→K′，谷间）
              X>Y   到星 Y 的全部成员（如 K→Q，谷间）
    态密度按星平均（等价谷的 DOS 相同，平均能压掉直方图噪声）：谷内用 G_星/g，X>X' 用 G_星·(g−1)/g，X>Y 用 G_Y。

step.conf 写法（D 的单位 eV/Å = 1e8 eV/cm；ħω 单位 meV）：
    DP_OPTICAL_E = K>K:4.1@50; K>K':2.6@46; K>Q:2.0@48; Q>Q:3.0@50; Q>Q':2.0@45; Q>K:2.0@48
    DP_OPTICAL_H = K>K:...; K>Γ:...; Γ>K:...
可选手填谷中心（不写就读 S8.2 的 dpt_result.json multi_valley 谷表）：
    DP_OPTICAL_VALLEYS_E = K:0.333333,0.333333,0|-0.333333,-0.333333,0; Q:...

[0036] 三维（S8 原版 AMSET 的 3D 路径，插件 step8_amset/amset_dpo3d.py）：模型与公式不变 —— 上式本来就与维数无关
（|M|² 里的 1/V 与对末态求和得到的 V 抵消，G 用每原胞态密度）。三维只改"谷"的几何：谷中心不限 kz = 0，
k 点归属与"同一个谷"的判断用完整三维距离（27 个周期像），函数统一加 full3d=True。full3d=False（出厂）与以前逐字节相同。
    DP_OPTICAL_VALLEYS_E = X:0.425,0,0.425|... ; L:0.5,0.5,0.5|...     # 3D 必须写 3 个分量（S3 晶胞分数坐标）
"""
import math
import re

import numpy as np

E_C = 1.602176634e-19        # J/eV
HA_J = 4.3597447222071e-18   # J/Hartree
HA_EV = HA_J / E_C
HBAR = 1.054571817e-34       # J·s
KB = 1.380649e-23            # J/K
AMU = 1.66053906660e-27      # kg
# 没有通道/在阈值下的态给一个可忽略的下限（1/s）。必须 > 1e6：AMSET 的 _interpolate_zero_rates 把
# ≤ 1e6 s⁻¹ 的率当"没算到"、用最近邻 k 点的率顶替 —— 对解析率来说阈值下的 0 是真 0，被顶替就错了。
# 2e6 s⁻¹（τ = 0.5 μs）比任何真实散射率小 5 个量级以上，对总率没有影响。
RATE_FLOOR = 2.0e6
DOS_BIN_EV = 0.002           # 态密度直方图步长
DOS_SIGMA_EV = 0.006         # 高斯展宽（压直方图噪声；台阶边缘相应变圆）
DOS_WINDOW_EV = 1.5          # 从带边往上统计多宽（再加最大 ħω）

_ALIAS = {"G": "Γ", "GAMMA": "Γ", "Γ": "Γ", "GM": "Γ–M", "Γ-M": "Γ–M", "Γ–M": "Γ–M",
          "G-M": "Γ–M", "G–M": "Γ–M", "GAMMA-M": "Γ–M"}          # S8.2 的标签：Γ/K/M/Q/Γ–M


def norm_label(s):
    """谷标签归一：Γ 的各种写法统一成 Γ；其余去空白、ASCII 字母转大写（K、Q、M…）。"""
    t = str(s).strip()
    u = t.upper()
    if u in _ALIAS:
        return _ALIAS[u]
    return u if re.fullmatch(r"[A-Za-z0-9_()'.,+-]+", t) else t


def parse_channels(spec):
    """'K>K:4.1@50; K>K':2.6@46; K>Q:2.0@48' -> [{"from","to","D_eV_A","hw_meV","kind"}]。
    kind：intra（X>X）、inter_star（X>X'）、inter（X>Y）。同一对谷的不同声子支各写一条（ħω 不同），率相加。
    格式错/数值非正/完全重复 -> ValueError（不静默）。"""
    out = []
    for item in [p.strip() for p in str(spec or "").replace("\n", ";").split(";")]:
        if not item:
            continue
        m = re.fullmatch(r"(.+?)\s*>\s*(.+?)\s*:\s*([0-9.eE+-]+)\s*@\s*([0-9.eE+-]+)", item)
        if not m:
            raise ValueError("通道 %r 写法不对，应为 起>止:D@ħω（D 单位 eV/Å，ħω 单位 meV），如 K>Q:2.0@48" % item)
        a, b = norm_label(m.group(1)), m.group(2).strip()
        D, hw = float(m.group(3)), float(m.group(4))
        if not (D > 0 and hw > 0 and math.isfinite(D) and math.isfinite(hw)):
            raise ValueError("通道 %r：D 和 ħω 必须为正数" % item)
        if b.endswith("'") or b.endswith("′"):
            bb = norm_label(b[:-1])
            if bb != a:
                raise ValueError("通道 %r：X' 只能指同一星的其它谷（写成 %s>%s'）" % (item, a, a))
            kind, b = "inter_star", bb
        else:
            b = norm_label(b)
            kind = "intra" if b == a else "inter"
        out.append({"from": a, "to": b, "D_eV_A": D, "hw_meV": hw, "kind": kind})
    # 同一对谷可以写多条（文献按声子支分列：TA/LA/LO… 各自的 D 与 ħω），各支的率相加；
    # 连 ħω 都相同的两条多半是笔误 -> 报错
    keys = [(c["from"], c["to"], c["kind"], c["hw_meV"]) for c in out]
    if len(set(keys)) != len(keys):
        raise ValueError("通道有重复（同一对谷、同一 ħω 写了两次）：%s" % spec)
    return out


def parse_valleys(spec, full3d=False):
    """'K:0.3333,0.3333,0|-0.3333,-0.3333,0; Q:...' -> [{"label","g","k_frac"}]。
    full3d=True（[0036] 三维）：每个 k 点必须写 3 个分量（2 个分量会被当成 kz=0，三维里多半是漏写）。"""
    out = []
    for item in [p.strip() for p in str(spec or "").replace("\n", ";").split(";")]:
        if not item:
            continue
        if ":" not in item:
            raise ValueError("谷 %r 写法不对，应为 标签:kx,ky,kz|kx,ky,kz（分数坐标）" % item)
        lab, rest = item.split(":", 1)
        ks = []
        for kk in rest.split("|"):
            v = [float(x) for x in kk.replace(" ", "").split(",") if x != ""]
            if len(v) == 2 and full3d:
                raise ValueError("谷 %r 的 k 点 %r 只有 2 个分量：三维谷中心必须写 kx,ky,kz" % (item, kk))
            if len(v) == 2:
                v.append(0.0)
            if len(v) != 3:
                raise ValueError("谷 %r 的 k 点 %r 不是 2 或 3 个分量" % (item, kk))
            ks.append(v)
        out.append({"label": norm_label(lab), "g": len(ks), "k_frac": ks})
    return out


def valleys_from_dpt(dpt, carrier):
    """从 S8.2 的 dpt_result.json 取某载流子的谷表（multi_valley.valleys；拒绝出数时谷表照样在）。"""
    for r in (dpt or {}).get("results", []):
        if r.get("carrier") != carrier:
            continue
        vs = ((r.get("multi_valley") or {}).get("valleys")) or []
        out = []
        for v in vs:
            ks = v.get("k_frac") or []
            if ks:
                row = {"label": norm_label(v.get("label")), "g": len(ks),
                       "k_frac": [[float(x) for x in k] for k in ks], "dE_eV": v.get("dE_eV")}
                kg = v.get("k_grid") or []
                if len(kg) == len(ks):           # [0034-r2] 网格上的极值点：判"同一个谷"用它（拟合中心按带各自浮动）
                    row["k_grid"] = [[float(x) for x in k] for k in kg]
                out.append(row)
        return out
    return []


def _same_k(a, b, tol=1e-4, full3d=False):
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d - np.rint(d)
    if not full3d:
        d[2:] = 0.0
    return bool(np.all(np.abs(d) < tol))


def grid_key(k, full3d=False):
    """网格点的比较键：面内分数坐标取 5 位、模 1（S8.2 的 k_grid 是网格上的精确点，按带/自旋不浮动）。
    full3d=True：三个分量都参与。"""
    out = []
    for x in list(k)[:3 if full3d else 2]:
        y = round(round(float(x), 5) % 1.0, 5)
        out.append(0.0 if y >= 1.0 else y)
    return tuple(out)


def _member_match(v, u, full3d=False):
    """v 的每个成员是否落在 u 的某个成员上：两边都有 k_grid 就按网格点精确比，否则按拟合中心（容差 1e-4）。"""
    if v.get("k_grid") and u.get("k_grid"):
        uk = {grid_key(w, full3d) for w in u["k_grid"]}
        return [grid_key(k, full3d) in uk for k in v["k_grid"]]
    return [any(_same_k(k, w, full3d=full3d) for w in u["k_frac"]) for k in (v.get("k_frac") or [])]


def merge_coincident(valleys, full3d=False):
    """S8.2 按（自旋, 能带）分星：两条带（或两个自旋）在同一组 k 上都有极值时，会出现同一组谷的两个星，
    第二个标成 "K(up带27)" 这样。本模型按"离谷中心最近"归属 k 点、DOS 统计所有带，所以它们是同一个谷：
    后出现的那个并进前一个（返回 (新谷表, {被并掉的标签: 保留的标签})）。只部分重合 -> ValueError（谷表自相矛盾）。
    [0034-r2] 判同优先用 S8.2 的 k_grid（网格极值点）：非高对称谷（Q 等）的拟合中心按带各自浮动、差一个网格步长的
    零头，按拟合中心判会漏合并（谷被中垂线劈成两半、一半拿不到通道）或误报部分重合。"""
    out, alias = [], {}
    for v in valleys:
        ks = v.get("k_frac") or []
        hit = None
        for u in out:
            same = _member_match(v, u, full3d)
            if same and all(same) and len(ks) == len(u["k_frac"]):
                hit = u
                break
            if any(same):
                raise ValueError("谷表里 %s 与 %s 只有部分中心重合（%d/%d）——谷表自相矛盾，请手填 DP_OPTICAL_VALLEYS_*"
                                 % (v.get("label"), u["label"], sum(same), len(ks)))
        if hit is None:
            out.append(dict(v))
        else:
            alias[v["label"]] = hit["label"]
    return out, alias


def check_channels(channels, valleys):
    """通道两端的谷都要在谷表里；X>X' 要求该星至少 2 个谷。不通过 -> ValueError。"""
    have = {v["label"]: v for v in valleys}
    for c in channels:
        for lab in (c["from"], c["to"]):
            if lab not in have:
                raise ValueError("通道 %s>%s：谷表里没有 %s（有：%s）" % (c["from"], c["to"], lab, ", ".join(have) or "无"))
        if c["kind"] == "inter_star" and have[c["from"]]["g"] < 2:
            raise ValueError("通道 %s>%s'：%s 星只有 1 个谷，没有同星的其它谷" % (c["from"], c["from"], c["from"]))


def member_table(valleys):
    """谷表 -> (中心 (n,3), 每个成员所属星的下标, 星标签列表, 每星成员数)。"""
    centers, star_of, labels, g = [], [], [], []
    for si, v in enumerate(valleys):
        labels.append(v["label"])
        g.append(len(v["k_frac"]))
        for k in v["k_frac"]:
            centers.append([float(x) for x in k])
            star_of.append(si)
    return np.asarray(centers, float), np.asarray(star_of, int), labels, np.asarray(g, int)


def _image_shifts(full3d):
    """最近像搜索用的分数坐标平移：2D 面内 9 个，3D 27 个。"""
    if full3d:
        return [(i, j, l) for i in (-1, 0, 1) for j in (-1, 0, 1) for l in (-1, 0, 1)]
    return [(i, j, 0) for i in (-1, 0, 1) for j in (-1, 0, 1)]


def assign_members(kfrac, centers, rlat, full3d=False):
    """每个 k 点（分数坐标 (n,3)）分给最近的谷中心（面内笛卡尔距离，含 ±1 周期像）。返回成员下标 (n,)。
    rlat：倒格矢矩阵（行 = b1,b2,b3，任意长度单位）。只用面内（z 分量置零：中心都在 kz=0 层）。
    full3d=True（[0036] 三维）：完整三维笛卡尔距离、27 个周期像。"""
    k = np.atleast_2d(np.asarray(kfrac, float))
    c = np.atleast_2d(np.asarray(centers, float))
    r = np.asarray(rlat, float)
    best = np.full(len(k), np.inf)
    idx = np.zeros(len(k), int)
    shifts = _image_shifts(full3d)
    for ci, cc in enumerate(c):
        d = k - cc
        d = d - np.rint(d)
        if not full3d:
            d[:, 2] = 0.0
        dmin = np.full(len(k), np.inf)
        for i, j, l in shifts:
            dd = d + np.array([i, j, float(l)])
            dmin = np.minimum(dmin, np.einsum("ni,ni->n", dd @ r, dd @ r))
        better = dmin < best
        best[better] = dmin[better]
        idx[better] = ci
    return idx


def same_member(k0, kps, centers, rlat, full3d=False):
    """k0 (3,) 与每个 kps (n,3)（分数坐标）是否分到同一个成员谷 -> bool (n,)。ADP 谷内掩码用。"""
    m0 = assign_members(np.asarray(k0, float)[None], centers, rlat, full3d)[0]
    return assign_members(kps, centers, rlat, full3d) == m0


def star_dos(x, star_of_k, n_star, n_k, grid):
    """按星统计态密度：x (nb, nk) 载流子能量（Hartree），star_of_k (nk,)。
    返回 (n_star, ngrid)：每原胞、每自旋、每 Hartree 的态数（直方图 + 高斯展宽）。"""
    x = np.asarray(x, float)
    dx = grid[1] - grid[0]
    edges = np.concatenate([grid - dx / 2, [grid[-1] + dx / 2]])
    sig = DOS_SIGMA_EV / HA_EV
    kern_x = np.arange(-4 * sig, 4 * sig + dx / 2, dx)
    kern = np.exp(-0.5 * (kern_x / sig) ** 2)
    kern /= kern.sum()
    out = np.zeros((n_star, len(grid)))
    for s in range(n_star):
        sel = x[:, star_of_k == s].ravel()
        h, _ = np.histogram(sel, bins=edges)
        out[s] = np.convolve(h, kern, mode="same") / (n_k * dx)
    return out


def occupation(x, x_f, kt):
    """载流子占据 f_c(x) = 1/(exp((x − x_F)/kT) + 1)（x、x_F、kT 同单位）；tanh 形式不溢出。
    电子 x = E、x_F = E_F；空穴 x = −E、x_F = −E_F（即 1 − f_e）。"""
    return 0.5 * (1.0 - np.tanh((np.asarray(x, float) - x_f) / (2.0 * kt)))


def bose(hw_J, T):
    T = np.asarray(T, float)
    with np.errstate(over="ignore", divide="ignore"):
        return 1.0 / np.expm1(hw_J / (KB * np.maximum(T, 1e-9)))


def carrier_rates(energies, kfrac, rlat, valleys, channels, temps, cell_mass_kg, kinds, sign=+1, full3d=False,
                  fermi_levels=None):
    """一种载流子（一组能带）的散射率。
    energies (nb, nk)：Hartree（AMSET 的 energies[spin][这些带]）；sign=+1 电子、−1 空穴。
    kinds：要算的通道种类集合（{"intra"} = ODP，{"inter_star","inter"} = IVS）。
    fermi_levels：None -> 非简并极限，返回 (nT, nb, nk)；
                  (ndop, nT) 的电子费米能级（Hartree，AMSET 的 amset_data.fermi_levels）-> [0042] 带泡利阻塞，
                  返回 (ndop, nT, nb, nk)。
    单位 1/s（没有通道的态为 0）。"""
    x = sign * np.asarray(energies, float)
    nb, nk = x.shape
    temps = np.atleast_1d(np.asarray(temps, float))
    fl = None if fermi_levels is None else np.atleast_2d(np.asarray(fermi_levels, float))
    if fl is not None and fl.shape[1] != len(temps):
        raise ValueError("fermi_levels 形状 %s 与温度数 %d 不符（应为 (ndop, nT)）" % (fl.shape, len(temps)))
    rates = np.zeros((len(temps), nb, nk) if fl is None else (fl.shape[0], len(temps), nb, nk))
    chs = [c for c in channels if c["kind"] in kinds]
    if not chs or not valleys:
        return rates
    centers, star_of_member, labels, g = member_table(valleys)
    member = assign_members(kfrac, centers, rlat, full3d)
    star_k = star_of_member[member]
    lab_idx = {lab: i for i, lab in enumerate(labels)}
    hw_max = max(c["hw_meV"] for c in chs) * 1e-3 / HA_EV
    x0 = float(np.nanmin(x))
    dx = DOS_BIN_EV / HA_EV
    grid = np.arange(x0 - 5 * DOS_SIGMA_EV / HA_EV, x0 + DOS_WINDOW_EV / HA_EV + 2 * hw_max, dx)
    G = star_dos(x, star_k, len(labels), nk, grid) / HA_J                 # 每原胞、每自旋、每焦耳
    for c in chs:
        si, ti = lab_idx[c["from"]], lab_idx[c["to"]]
        if c["kind"] == "intra":
            Gt = G[ti] / g[ti]
        elif c["kind"] == "inter_star":
            Gt = G[ti] * (g[ti] - 1) / g[ti]
        else:
            Gt = G[ti]
        D = c["D_eV_A"] * E_C / 1e-10                                     # J/m
        hw_J = c["hw_meV"] * 1e-3 * E_C
        om = hw_J / HBAR
        hw_H = hw_J / HA_J
        pref = math.pi * D * D / (cell_mass_kg * om)
        sel = (star_k == si)
        if not sel.any():
            continue
        xs = x[:, sel]
        g_abs = np.interp(xs + hw_H, grid, Gt, left=0.0, right=Gt[-1])
        g_em = np.interp(xs - hw_H, grid, Gt, left=0.0, right=Gt[-1])
        N = bose(hw_J, temps)
        if fl is None:
            r = pref * (N[:, None, None] * g_abs[None] + (N[:, None, None] + 1.0) * g_em[None])
            rates[:, :, sel] += r
            continue
        for it, T in enumerate(temps):                    # [0042] 末态占据：吸收 N + f(x+ħω)，发射 N + 1 − f(x−ħω)
            kt = KB * max(float(T), 1e-9) / HA_J
            for idop in range(fl.shape[0]):
                xf = sign * fl[idop, it]
                f_up = occupation(xs + hw_H, xf, kt)
                f_dn = occupation(xs - hw_H, xf, kt)
                rates[idop, it][:, sel] += pref * ((N[it] + f_up) * g_abs + (N[it] + 1.0 - f_dn) * g_em)
    return rates


def valley_report(energies, kfrac, rlat, valleys, sign=+1, full3d=False):
    """谷表与能带对得上吗：每个星的区域里载流子能量最低的 k 点，离本星最近成员中心多远（|Δk|/|b1|），
    比整体带边高多少（eV）。谷中心与 AMSET 网格错位（晶胞设置不同、手填错）时距离会很大。
    energies (nb, nk) Hartree；返回 [{"label", "dk_rel", "dE_eV"}]。"""
    x = (sign * np.asarray(energies, float)).min(axis=0)
    k = np.atleast_2d(np.asarray(kfrac, float))
    r = np.asarray(rlat, float)
    centers, star_of_member, labels, _ = member_table(valleys)
    member = assign_members(k, centers, r, full3d)
    star = star_of_member[member]
    b1 = float(np.linalg.norm(r[0]))
    gmin = float(np.nanmin(x))
    out = []
    for si, lab in enumerate(labels):
        idx = np.where(star == si)[0]
        if not len(idx):
            out.append({"label": lab, "dk_rel": None, "dE_eV": None})
            continue
        j = idx[int(np.argmin(x[idx]))]
        d = k[j] - centers[member[j]]
        d = d - np.rint(d)
        if not full3d:
            d[2] = 0.0
        dk = min(float(np.linalg.norm((d + np.array([i, jj, float(ll)])) @ r))
                 for i, jj, ll in _image_shifts(full3d))
        out.append({"label": lab, "dk_rel": dk / b1, "dE_eV": (float(x[j]) - gmin) * HA_EV})
    return out


def analytic_rate(x_J, valley_min_J, m_d, g_target, area_m2, D_eV_A, hw_meV, T, cell_mass_kg):
    """解析参考（抛物带、二维常数态密度）：给单元测试和量级估算用。
    x_J：初态载流子能量（J）；valley_min_J：目标谷底（J）；g_target：目标谷集合的谷数。"""
    G = g_target * m_d * 9.1093837015e-31 * area_m2 / (2 * math.pi * HBAR ** 2)   # 每原胞每自旋每焦耳
    D = D_eV_A * E_C / 1e-10
    hw = hw_meV * 1e-3 * E_C
    om = hw / HBAR
    N = 1.0 / math.expm1(hw / (KB * T))
    ga = G if x_J + hw >= valley_min_J else 0.0
    ge = G if x_J - hw >= valley_min_J else 0.0
    return math.pi * D * D / (cell_mass_kg * om) * (N * ga + (N + 1) * ge)


def analytic_rate_3d(x_J, valley_min_J, m_d, g_target, volume_m3, D_eV_A, hw_meV, T, cell_mass_kg):
    """[0036] 解析参考（抛物带、三维 √E 态密度）：给单元测试和量级估算用（Conwell 的非极性光学声子公式）。
    volume_m3：原胞体积；m_d：态密度有效质量（m_e）；其余同 analytic_rate。"""
    pref_g = g_target * volume_m3 * (2 * m_d * 9.1093837015e-31) ** 1.5 / (4 * math.pi ** 2 * HBAR ** 3)

    def G(e):                                      # 每原胞每自旋每焦耳
        return pref_g * math.sqrt(e) if e > 0 else 0.0
    D = D_eV_A * E_C / 1e-10
    hw = hw_meV * 1e-3 * E_C
    om = hw / HBAR
    N = 1.0 / math.expm1(hw / (KB * T))
    return math.pi * D * D / (cell_mass_kg * om) * (N * G(x_J + hw - valley_min_J)
                                                    + (N + 1) * G(x_J - hw - valley_min_J))


def describe(rec):
    """2d_correction.json 的 dp_optical 记录 -> 一行横幅文字。"""
    parts = []
    for car in ("electron", "hole"):
        c = (rec or {}).get(car) or {}
        chs = c.get("channels") or []
        if chs:
            parts.append("%s[%s] %s" % (car, ",".join("%s×%d" % (v["label"], v["g"]) for v in c.get("valleys") or []),
                                        "; ".join("%s>%s%s:%.3g@%.3g" % (h["from"], h["to"],
                                                                         "'" if h["kind"] == "inter_star" else "",
                                                                         h["D_eV_A"], h["hw_meV"]) for h in chs)))
    return " | ".join(parts) or "无"
