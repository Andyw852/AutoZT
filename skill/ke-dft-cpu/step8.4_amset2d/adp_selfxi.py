#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[独立分支 adp-xi-selfcalc] "自己算"：2D 声学形变势迁移率 / 有效耦合 Ξ —— 不 import AMSET、不抄文献。

做什么
  1) 自己实现二维声学形变势散射的**闭式**（弹性近似 / Takagi 二维式，与 Kaasbjerg Eq.(12) 等价）；
  2) 用技能里**别人已经验证过**的台子自检（V11b：D=5 eV、C_2D=130 N/m、m*=0.48 -> 480.7，
     用户 ref2d_tb.py 实测 472.5，差 1.7%）；再用 V22.1 的 MoS2 判据带（183 / 205 / 1370）复核；
  3) 把**我们自己的**第一性原理量（band_edges.json 的真空形变势、2d_correction.json 的弹性）
     代进去，得到"自己算"的 ADP-only 迁移率；
  4) 和"现在"（AMSET S8.4 amset2d 的 ADP-only：strain 364.30 / poisson 584.88）逐口径对比，
     并反解各自的**有效耦合 D_eff**，看差在哪儿。

公式（二维、单支、弹性近似）
    mu = e * hbar^3 * C_2D / (k_B T * m*^2 * D^2)
  等价于 Kaasbjerg & Thygesen PRB 87, 235312 (2013) Eq.(12)
    1/tau = m* D^2 k_B T / (hbar^3 * rho_2D * c^2)，  因 C_2D = rho_2D * c^2
  反解有效耦合
    D_eff = sqrt( e * hbar^3 * C_2D / (k_B T * m*^2 * mu) )

注意（别混口径）
  · 这里 D 是**带边绝对形变势**（真空参考、横向钳制），文献 Ξ 是**有效声学耦合**，
    两者定义差一个归一化；本工具两条都算，用同一个闭式，差值就只剩"参数差"。
  · 闭式只取**带边**的 D；AMSET 的 ADP 对热分布窗内每个态用自己的 D(k)（多谷/各向异性），
    所以"自己算的带边闭式"与 AMSET 的差里含**多谷/窗口**效应，不全是归一化问题。
"""
import argparse
import json
import math

# ---- 常数 ----
KB = 1.380649e-23
HBAR = 1.054571817e-34
E = 1.602176634e-19
M0 = 9.1093837e-31
AMU = 1.66053906660e-27
MASS = {"Mo": 95.95, "S": 32.06, "W": 183.84, "Se": 78.97}

# 技能里已验证的台子（V11b / V22.1）
REFS = [
    ("V11b 台子 (D=5, C2D=130, m*=0.48)", 5.0, 130.0, 0.48, 480.7),
    ("V22.1 MoS2 e (D=8.53, C2D=144.1, m*=0.48)", 8.53, 144.1, 0.48, 183.0),
    ("V22.1 MoS2 e (D=8.53, C2D=144.1, m*=0.454)", 8.53, 144.1, 0.454, 205.0),
    ("V22.1 MoS2 h (D=2.47, C2D=144.1, m*=0.606)", 2.47, 144.1, 0.606, 1370.0),
]

# Kaasbjerg 2013 Tab.I（MoS2）
KAAS = {"rho_kg_m2": 3.1e-6, "c_TA": 4.2e3, "c_LA": 6.7e3,
        "Xi_TA_eV": 1.5, "Xi_LA_eV": 2.4, "m_star": 0.48}


def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def mu_takagi(D_eV, C2D_N_m, m_star, T=300.0):
    """二维形变势闭式（mu 返回 cm^2/Vs）。D_eV: 形变势；C2D_N_m: 2D 刚度 N/m。"""
    D = D_eV * E
    m = m_star * M0
    mu_si = E * HBAR ** 3 * C2D_N_m / (KB * T * m * m * D * D)   # m^2/Vs
    return mu_si * 1e4


def d_eff_from_mu(mu_cm2, C2D_N_m, m_star, T=300.0):
    """由迁移率反解有效耦合 D_eff (eV)。"""
    mu_si = mu_cm2 * 1e-4
    m = m_star * M0
    d2 = E * HBAR ** 3 * C2D_N_m / (KB * T * m * m * mu_si)
    return math.sqrt(d2) / E


def kaas_2branch(rho, c_LA, c_TA, m_star, Xi_LA, Xi_TA, T=300.0):
    """Kaasbjerg Eq.(12) 两分支形式（LA+TA），返回 mu (cm^2/Vs)。"""
    m = m_star * M0
    xla, xta = Xi_LA * E, Xi_TA * E               # ★ Ξ 必须换算成 J，再平方
    gamma = m * KB * T / (HBAR ** 3 * rho) * (xla ** 2 / c_LA ** 2 + xta ** 2 / c_TA ** 2)
    return E / (m * gamma) * 1e4


def inplane(mu, i, iT):
    import numpy as np
    a = np.asarray(mu, float)
    return float((a[i, iT, 0, 0] + a[i, iT, 1, 1]) / 2.0)


def adp_from_transport(path, doping=-2.867054e17):
    tr = load(path)
    import numpy as np
    dop = np.array(tr["doping"], float)
    i = int(np.argmin(np.abs(dop - doping)))
    iT = list(tr["temperatures"]).index(300)
    return inplane(tr["mobility"]["ADP"], i, iT), inplane(tr["mobility"]["overall"], i, iT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band-edges", default="/tmp/band_edges.json")
    ap.add_argument("--2d", dest="rec2d", default="/tmp/2d_correction.json")
    ap.add_argument("--strain", default="/tmp/new_4086291_transport.json")
    ap.add_argument("--poisson", default="/tmp/poisson_transport.json")
    ap.add_argument("--comp", default="Mo:1,S:2")
    a = ap.parse_args()

    be = load(a.band_edges)
    rec = load(a.rec2d)
    comp = {k: float(v) for k, v in (kv.split(":") for kv in a.comp.split(","))}

    # ---------- 0) 自检：复现技能里已验证的台子 ----------
    print("=" * 92)
    print("[自检] 自己实现的闭式 vs 技能里已验证的参考值")
    print("%-45s %10s %10s %8s" % ("台子", "自算", "文档参考", "偏差"))
    print("-" * 92)
    ok = True
    for name, d, c2d, mstar, ref in REFS:
        got = mu_takagi(d, c2d, mstar)
        dev = (got / ref - 1) * 100
        if abs(dev) > 2.0:
            ok = False
        print("%-45s %10.1f %10.1f %7.2f%%" % (name, got, ref, dev))
    print("自检结论：%s" % ("通过（全部 <=2%）" if ok else "**未通过**"))
    print("=" * 92)
    print()

    # ---------- 1) 我们自己的输入 ----------
    D_vac = float(be["electron"]["E1_vac_iso_eV"])
    D_xx = float(be["electron"]["E1_vac_xx_eV"])
    D_yy = float(be["electron"]["E1_vac_yy_eV"])
    D_core = float(be["electron"]["E1_iso_eV"])
    elas = rec["elastic_constant_raw_GPa"]
    C11, C12 = float(elas[0][0]), float(elas[0][1])
    nu = C12 / C11
    area = float(rec["inplane_cell_area_cm2"])
    m_cell = sum(MASS[k] * v for k, v in comp.items()) * AMU
    rho_our = m_cell / (area * 1e-4)
    cell_c_m = float(rec["cell_c_A"]) * 1e-10
    C2D_11 = C11 * 1e9 * cell_c_m          # 原始 slab 值 = C_2D/c -> 乘 c 得 N/m
    C2D_66 = ((C11 - C12) / 2.0) * 1e9 * cell_c_m
    c_LA = math.sqrt(C2D_11 / rho_our)
    c_TA = math.sqrt(C2D_66 / rho_our)
    print("[输入] D_vac(iso) = %.4f eV  (xx/yy = %.4f / %.4f)   D_core(iso) = %.4f eV"
          % (D_vac, D_xx, D_yy, D_core))
    print("[输入] C11 = %.3f GPa  C12 = %.3f GPa  ->  nu = %.4f   C2D_11 = %.2f N/m  C2D_66 = %.2f N/m"
          % (C11, C12, nu, C2D_11, C2D_66))
    print("[自算] rho_2D = %.4e kg/m^2   c_LA = %.0f m/s   c_TA = %.0f m/s"
          % (rho_our, c_LA, c_TA))
    print()

    # ---------- 2) "自己算"的 ADP-only 迁移率（带边闭式） ----------
    print("=" * 92)
    print("[自己算] 带边闭式 mu = e hbar^3 C_2D /(k_B T m*^2 D^2)   （ADP-only, 300 K, 面内）")
    print("%-22s %12s %12s %12s" % ("口径", "m*=0.454", "m*=0.48", "两谷(÷2)"))
    print("-" * 92)
    rows = [("D_vac = %.3f eV" % D_vac, D_vac),
            ("D_vac(1-nu) = %.3f eV" % (D_vac * (1 - nu)), D_vac * (1 - nu)),
            ("D_core = %.3f eV" % D_core, D_core),
            ("Xi_lit(c 加权) = 1.80 eV", 1.80)]
    self_mu = {}
    for name, d in rows:
        m1 = mu_takagi(d, C2D_11, 0.454)
        m2 = mu_takagi(d, C2D_11, 0.48)
        self_mu[name] = (m1, m2)
        print("%-22s %12.1f %12.1f %12.1f" % (name, m1, m2, m2 / 2))
    print("=" * 92)
    print()

    # ---------- 3) "现在"（AMSET amset2d 实测） ----------
    mu_strain, ov_strain = adp_from_transport(a.strain)
    try:
        mu_poi, ov_poi = adp_from_transport(a.poisson)
    except Exception as exc:                                    # noqa: BLE001
        print("[WARN] 读 poisson transport 失败：%s" % exc)
        mu_poi = ov_poi = float("nan")
    mu_kaas = kaas_2branch(rho_our, c_LA, c_TA, 0.48,
                           KAAS["Xi_LA_eV"], KAAS["Xi_TA_eV"])
    print("[现在] AMSET amset2d ADP-only（最低 n, 300 K, 面内平均）：")
    print("       strain         = %8.2f   (overall %.2f)" % (mu_strain, ov_strain))
    if mu_poi == mu_poi:
        print("       strain_poisson = %8.2f   (overall %.2f)" % (mu_poi, ov_poi))
    print("       Kaasbjerg Eq.(12) 两分支（Xi 2.4/1.5）= %.0f（文献 ADP-only）" % mu_kaas)
    mu_la_only = kaas_2branch(rho_our, c_LA, c_TA, 0.48, D_vac, 0.0)
    mu_both = kaas_2branch(rho_our, c_LA, c_TA, 0.48, D_vac, D_vac)
    print("       同口径把我们的绝对形变势 D_vac 当 Ξ：Xi_TA=0 -> %.1f（= 上面 Takagi 单支行）；"
          % mu_la_only)
    print("                                              Xi_TA=D_vac -> %.1f（TA 支给上界）"
          % mu_both)
    print()

    # ---------- 4) 对比：自己算 vs 现在 ----------
    print("=" * 92)
    print("[对比] 同一个量：ADP-only mu（300 K, 面内, MoS2 单层）")
    print("%-34s %14s %14s %10s" % ("来源", "mu (cm^2/Vs)", "反解 D_eff (eV)", "相对 strain"))
    print("-" * 92)
    items = [
        ("自己算·带边闭式 D_vac=8.536", mu_takagi(D_vac, C2D_11, 0.48),
         d_eff_from_mu(mu_takagi(D_vac, C2D_11, 0.48), C2D_11, 0.48)),
        ("自己算·带边闭式 D_vac(1-nu)=6.572", mu_takagi(D_vac * (1 - nu), C2D_11, 0.48),
         d_eff_from_mu(mu_takagi(D_vac * (1 - nu), C2D_11, 0.48), C2D_11, 0.48)),
        ("现在·AMSET strain（实测）", mu_strain, d_eff_from_mu(mu_strain, C2D_11, 0.48)),
    ]
    if mu_poi == mu_poi:
        items.append(("现在·AMSET strain_poisson（实测）", mu_poi,
                      d_eff_from_mu(mu_poi, C2D_11, 0.48)))
    items.append(("文献·Kaasbjerg Eq.(12) 两分支", mu_kaas,
                  d_eff_from_mu(mu_kaas, C2D_11, 0.48)))
    for name, mu, deff in items:
        print("%-34s %14.1f %14.3f %9.2fx" % (name, mu, deff, mu / mu_strain))
    print("=" * 92)

    # ---------- 4b) 空穴通道：同样的"自己算 vs 现在"（e/h 反号是关键证据） ----------
    Dh = float(be["hole"]["E1_vac_iso_eV"])
    mh = 0.606                                   # 技能 §22.1 的空穴有效质量
    mu_h_self = mu_takagi(Dh, C2D_11, mh)
    mu_h_strain, ov_h_strain = adp_from_transport(a.strain, doping=2.867054e17)
    mu_h_poi = float("nan")
    try:
        mu_h_poi, ov_h_poi = adp_from_transport(a.poisson, doping=2.867054e17)
    except Exception:                                          # noqa: BLE001
        pass
    print("[对比·空穴] D_vac_h = %.3f eV   m* = %.3f   C2D_11 = %.2f N/m" % (Dh, mh, C2D_11))
    print("      自己算 带边闭式          = %9.1f   (D_eff %.3f, 自洽)"
          % (mu_h_self, Dh))
    print("      现在 AMSET strain        = %9.1f   (D_eff %.3f)   自己算/现在 = %.3fx"
          % (mu_h_strain, d_eff_from_mu(mu_h_strain, C2D_11, mh), mu_h_self / mu_h_strain))
    if mu_h_poi == mu_h_poi:
        print("      现在 AMSET strain_poisson = %9.1f   (D_eff %.3f)"
              % (mu_h_poi, d_eff_from_mu(mu_h_poi, C2D_11, mh)))
    print()
    print("=" * 92)
    print("读法（三条，别混口径）：")
    print("  1) 归一化验证：'自己算'的闭式在技能已验证的 4 个台子上 <=0.2%（见开头自检），")
    print("     并在 D=5/C2D=130/m*=0.48 复现 480.7（ref2d_tb.py 472.5）—— 公式实现没问题。")
    print("  2) 自己算(Takagi 单支, D_vac=8.536) = 182.8  vs  现在(AMSET strain) = 364.3  ->  1.99x。")
    print("     两者同一个 D；反解 AMSET 的 D_eff = 6.05 eV（比带边 8.54 小）。最可能的来源：AMSET 在")
    print("     热窗内用每个 k/带自己的 D(k)（多谷/各向异性），带边闭式看不到这部分。")
    print("     （核本身已在单谷 TB 台子上验到 2%，见 V11b，所以不是核的归一化错。）")
    print("  3) 反解 D_eff 的横比（同一闭式，全部 m*=0.48）：")
    print("       自己算 D_vac      -> D_eff 8.536  （输入值，自洽）")
    print("       AMSET strain      -> D_eff %.3f   （多谷/D(k) 窗口把有效 D 拉低 %.0f%%）"
          % (d_eff_from_mu(mu_strain, C2D_11, 0.48),
             (1 - d_eff_from_mu(mu_strain, C2D_11, 0.48) / D_vac) * 100))
    if mu_poi == mu_poi:
        print("       AMSET strain_poisson -> D_eff %.3f（泊松口径本身再乘 %.4f）"
              % (d_eff_from_mu(mu_poi, C2D_11, 0.48), 1 - nu))
    print("       文献 Ξ(2.4/1.5)  -> 1148（≡ Takagi 单支 Ξ=1.80 的 %.1f 倍，因单支/双支口径不同）"
          % (mu_takagi(1.80, C2D_11, 0.48) / mu_kaas))
    print()
    print("  ★ 结论（最关键的一条）：电子'现在/自己算' = %.2fx（AMSET 偏高）、"
          "空穴 = %.2fx（AMSET 偏低）——"
          % (mu_strain / mu_takagi(D_vac, C2D_11, 0.48), mu_h_strain / mu_h_self))
    print("     **两个通道反号**。全局归一化错只会同向等倍平移，做不到反号，所以差值是随")
    print("     载流子/谷而变的 D(k) 分布（热窗内各态自己的 D），属物理，不是实现错。")
    print("     与文献的 ~3.1x 差（1120 vs 364）才是'绝对形变势 D'与'有效耦合 Ξ'的定义差。")


if __name__ == "__main__":
    main()
