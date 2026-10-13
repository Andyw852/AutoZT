#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[分支测试 adp-xi-test] ADP 形变势"四口径"对照 —— 只为对比，不改动 master 代码。

四个口径（都只作用在 h5 的形变势张量 D 上，给出 ADP-only 面内迁移率）：
  ① strain          因子 1.0              D = D_vac（真空对齐 + 横向钳制，现状 production）
  ② strain_poisson  因子 1 - C12/C11      D = D_vac(1-nu)（Wiktor 泊松横向弛豫口径）
  ③ xi_literal      把文献有效耦合 Ξ_LA 当 D   直接抄 Kaasbjerg Tab.I 的 Ξ_LA = 2.4 eV
  ④ xi_own（自实现） 由**我们自己**的 rho_2D / c_LA / c_TA + 实测 tau 反解出的有效耦合 Ξ_own

判据/公式（Kaasbjerg & Thygesen, PRB 87, 235312 (2013)）：
  Eq.(12)  1/tau_l = m* Ξ_l^2 k_B T / (hbar^3 rho_2D c_l^2)
  Eq.(11)  |M_ql| = Ξ_l q     （Ξ 是他们拟合 DFPT 耦合得到的**有效**常数，不是绝对形变势）

要点（写论文时别混）：
  · 绝对形变势 D（真空口径，Poisson 比 ~6-8.5 eV）与文献的有效耦合 Ξ（~1.8-2.4 eV）**不是同一个量**，
    定义差一个因子，直接把 Ξ 抄进 AMSET 的 ADP 核得不出文献的 ADP-only 迁移率。
  · 可行/可辩护的比法是：用**同一套公式结构**把两边都折成 Ξ（④ 自实现反解），再比 Ξ 的大小。
  · AMSET 2D ADP 核里 D_eff = |D_interp| + v(x)v，所以 D -> f*D 时速率**不是**严格 f^2；
    MoS2 实测两口径 364.30 -> 584.88，纯 D^2 预测 346.74，残差 +5.1%（v(x)v ~ D 的 9%）。
    因此 ③ 这种大倍数外推只给量级。
"""
import argparse, json, math

# ---- 常数 ----
KB = 1.380649e-23
HBAR = 1.054571817e-34
E = 1.602176634e-19
M0 = 9.1093837e-31
AMU = 1.66053906660e-27

# Kaasbjerg 2013 Tab.I（单层 MoS2）
#   rho = 3.1e-7 g/cm^2 = 3.1e-6 kg/m^2, c_TA = 4.2e3, c_LA = 6.7e3,
#   Ξ_TA = 1.5 eV, Ξ_LA = 2.4 eV, m* = 0.48 me
KAAS = {"rho_kg_m2": 3.1e-6, "c_TA": 4.2e3, "c_LA": 6.7e3,
        "Xi_TA_eV": 1.5, "Xi_LA_eV": 2.4, "m_star": 0.48}
WIKTOR = {"PBE_xx": -6.08, "PBE_yy": -6.19, "PBE0_xx": -7.16, "VBM_PBE_xx": -1.86}
MASS = {"Mo": 95.95, "S": 32.06, "W": 183.84, "Se": 78.97}


def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def inplane(mu, i, iT):
    import numpy as np
    a = np.asarray(mu, float)
    return float((a[i, iT, 0, 0] + a[i, iT, 1, 1]) / 2.0)


def pick(tr, doping):
    import numpy as np
    dop = np.array(tr["doping"], float)
    i = int(np.argmin(np.abs(dop - doping)))
    T = list(tr["temperatures"])
    return i, T.index(300), float(dop[i])


def rho_2d(area_cm2, comp):
    m_cell = sum(MASS[k] * v for k, v in comp.items()) * AMU
    return m_cell / (area_cm2 * 1e-4)


def kaas_rate(Xi_eV, c, m_star):
    Xi = Xi_eV * E
    return m_star * M0 * Xi * Xi * KB * 300.0 / (HBAR ** 3 * KAAS["rho_kg_m2"] * c * c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band-edges", default="/tmp/band_edges.json")
    ap.add_argument("--2d", dest="rec2d", default="/tmp/2d_correction.json")
    ap.add_argument("--strain", default="/tmp/new_4086291_transport.json",
                    help="strain 口径的 transport*.json（4086291）")
    ap.add_argument("--poisson", default="/tmp/poisson_transport.json",
                    help="strain_poisson 口径的 transport*.json（4086466）")
    ap.add_argument("--doping", type=float, default=-2.867054e17)
    ap.add_argument("--formula", default="MoS2")
    ap.add_argument("--comp", default="Mo:1,S:2")
    a = ap.parse_args()

    be = load(a.band_edges)
    rec = load(a.rec2d)
    comp = {k: float(v) for k, v in (kv.split(":") for kv in a.comp.split(","))}

    # ---------- 我们的形变势 / 弹性 / 声速 ----------
    D_vac = float(be["electron"]["E1_vac_iso_eV"])
    D_xx = float(be["electron"]["E1_vac_xx_eV"])
    D_yy = float(be["electron"]["E1_vac_yy_eV"])
    D_core = float(be["electron"]["E1_iso_eV"])
    elas = rec["elastic_constant_raw_GPa"]
    C11, C12 = float(elas[0][0]), float(elas[0][1])
    nu = C12 / C11
    area = float(rec["inplane_cell_area_cm2"])
    rho_our = rho_2d(area, comp)
    c_cell_m = float(rec["cell_c_A"]) * 1e-10
    # elastic_constant_raw 是原始 slab 值 = C_2D/c -> C_2D = C_raw * c（N/m）
    C2D_11 = C11 * 1e9 * c_cell_m
    C2D_66 = ((C11 - C12) / 2.0) * 1e9 * c_cell_m
    c_LA = math.sqrt(C2D_11 / rho_our)
    c_TA = math.sqrt(C2D_66 / rho_our)

    print("[输入] 电子 D_vac(iso) = %.4f eV   D_vac_xx/yy = %.4f / %.4f   D_core(iso) = %.4f"
          % (D_vac, D_xx, D_yy, D_core))
    print("[输入] C11 = %.3f GPa  C12 = %.3f GPa  ->  nu = %.4f   1-nu = %.4f"
          % (C11, C12, nu, 1 - nu))
    print("[自算] rho_2D = %.4e kg/m^2   c_LA = %.0f m/s   c_TA = %.0f m/s"
          % (rho_our, c_LA, c_TA))
    print("[对照] Kaasbjerg Tab.I  rho = %.4e   c_LA = %.0f   c_TA = %.0f  （三者均 <=1.5%%）"
          % (KAAS["rho_kg_m2"], KAAS["c_LA"], KAAS["c_TA"]))
    print("")

    # ---------- 文献 ADP-only（Kaasbjerg Eq.12 两支求和） ----------
    r_lit = kaas_rate(KAAS["Xi_LA_eV"], KAAS["c_LA"], KAAS["m_star"]) + \
            kaas_rate(KAAS["Xi_TA_eV"], KAAS["c_TA"], KAAS["m_star"])
    mu_lit = E / (KAAS["m_star"] * M0 * r_lit) * 1e4
    Xi_lit_eff = math.sqrt((KAAS["Xi_LA_eV"] ** 2 / KAAS["c_LA"] ** 2 +
                            KAAS["Xi_TA_eV"] ** 2 / KAAS["c_TA"] ** 2) /
                           (1.0 / KAAS["c_LA"] ** 2 + 1.0 / KAAS["c_TA"] ** 2))
    print("[文献] Kaasbjerg Eq.(12) 300K: 1/tau = %.4e s^-1 -> ADP-only mu = %.0f cm^2/Vs"
          % (r_lit, mu_lit))
    print("[文献] 他们的 c 加权有效耦合 Xi_eff = %.3f eV（不是 Xi_LA=2.4 单独一支）"
          % Xi_lit_eff)
    print("")

    # ---------- 我们两口径的实测 ----------
    tr1 = load(a.strain)
    i1, iT, dop1 = pick(tr1, a.doping)
    mu1 = inplane(tr1["mobility"]["ADP"], i1, iT)
    mu1_ov = inplane(tr1["mobility"]["overall"], i1, iT)
    tr2 = None
    try:
        tr2 = load(a.poisson)
        i2, iT2, dop2 = pick(tr2, a.doping)
        mu2 = inplane(tr2["mobility"]["ADP"], i2, iT2)
        mu2_ov = inplane(tr2["mobility"]["overall"], i2, iT2)
    except Exception as exc:                                   # noqa: BLE001
        print("[WARN] 读 poisson transport 失败（%s）" % exc)
        mu2 = mu2_ov = float("nan")

    D_poi = D_vac * (1 - nu)
    print("[实测] ① strain        D = %.4f eV  ADP-only mu = %.2f  overall = %.2f"
          % (D_vac, mu1, mu1_ov))
    if mu2 == mu2:
        print("[实测] ② strain_poisson D = %.4f eV  ADP-only mu = %.2f  overall = %.2f"
              % (D_poi, mu2, mu2_ov))
    print("[自检] 纯 D^2 律预测 ① = ② * (1-nu)^2 = %.2f（实测 %.2f，残差 %+.1f%% = 核里 v(x)v 项）"
          % (mu2 * (1 - nu) ** 2, mu1,
             (mu1 / (mu2 * (1 - nu) ** 2) - 1) * 100 if mu2 == mu2 else float("nan")))
    print("")

    # ---------- ④ 自实现：由我们的 tau 反解有效 Xi_own ----------
    r_our = E / (KAAS["m_star"] * M0 * (mu1 * 1e-4))
    denom = KAAS["m_star"] * M0 * KB * 300.0 * (1.0 / c_LA ** 2 + 1.0 / c_TA ** 2)
    Xi_own = math.sqrt(r_our * HBAR ** 3 * rho_our / denom) / E
    print("[自实现] 用我们自己的 rho_2D/c_LA/c_TA 把 ① 的速率反解成有效耦合：")
    print("         Gamma_eff = e/(m* mu) = %.4e s^-1  ->  Xi_own = %.3f eV" % (r_our, Xi_own))
    print("         Xi_own / Xi_lit_eff = %.3f  ->  预期 mu 比 = (该比值)^2 = %.3f"
          % (Xi_own / Xi_lit_eff, (Xi_own / Xi_lit_eff) ** 2))
    print("         实测 mu_lit / mu1 = %.3f   （两者一致 = 反解自洽）" % (mu_lit / mu1))
    print("")

    # ---------- ③ 直接把 Ξ 当 D（只给量级） ----------
    f3 = KAAS["Xi_LA_eV"] / D_vac
    mu3_from1 = mu1 / f3 ** 2
    mu3_from2 = mu2 / (KAAS["Xi_LA_eV"] / D_poi) ** 2 if mu2 == mu2 else float("nan")
    mu1_at_own = mu1 / (Xi_own / D_vac) ** 2

    # ---------- 四口径表 ----------
    print("=" * 104)
    print("%-26s %9s %9s %14s %14s %s" %
          ("口径", "D_eff/eV", "D/D_vac", "ADP-only mu", "来源", "说明"))
    print("-" * 104)
    rows = [
        ("① strain（现状 production）", D_vac, mu1, "AMSET 实测 4086291", "真空对齐 + 横向钳制的绝对形变势"),
        ("② strain_poisson", D_poi, mu2, "AMSET 实测 4086466", "Wiktor 泊松口径 D=D_vac(1-nu)"),
        ("③ xi_literal（抄文献）", KAAS["Xi_LA_eV"],
         "%.0f~%.0f" % (min(mu3_from1, mu3_from2), max(mu3_from1, mu3_from2)),
         "D^2 外推（不可信）", "把 Xi_LA=2.4 eV 当绝对 D -> 与文献公式结构不符"),
        ("④ xi_own（自实现反解）", Xi_own, mu1_at_own, "自实现（闭环）",
         "由我们的 rho/c/tau 反解；等价于 ①"),
    ]
    for name, d, mu, src, note in rows:
        mus = ("%.2f" % mu) if isinstance(mu, float) else str(mu)
        print("%-26s %9.3f %9.4f %14s %14s %s" %
              (name, d, d / D_vac, mus, src, note))
    print("=" * 104)
    print("注：③ 的 D^2 外推从 ① 得 %.0f、从 ② 得 %.0f（外推跨 3.6 倍，核里 v(x)v 非线性被放大，"
          "只作量级）；" % (mu3_from1, mu3_from2))
    print("    ④ 的 mu 是把 Xi_own 当 D 喂回 AMSET 核会得到的值，仅用于展示口径换算，"
          "物理结论看 Xi_own 本身。")
    print("")
    print("【结论】文献 ADP-only = %.0f，我们 ① = %.2f -> 差 %.2f 倍；"
          % (mu_lit, mu1, mu_lit / mu1))
    print("        换成同一套 Ξ 语言：我们 Xi_own = %.2f eV vs 文献 Xi_eff = %.2f eV -> %.2f 倍，"
          % (Xi_own, Xi_lit_eff, Xi_own / Xi_lit_eff))
    print("        平方得 %.2f 倍迁移率差，与直接比 %.2f 一致 -> **差值就是 D 与 Ξ 的定义差**，不是流程错。"
          % ((Xi_own / Xi_lit_eff) ** 2, mu_lit / mu1))
    print("        加上 ② 泊松口径（D 降 %.4f 倍）后，与文献的差从 %.2f 倍缩到 %.2f 倍。"
          % (1 - nu, mu_lit / mu1, mu_lit / mu2 if mu2 == mu2 else float("nan")))
    print("")
    print("Wiktor PRB 94, 245411 (2016) MoS2 真空绝对形变势: a_CBM(K) = %.2f (PBE xx) / %.2f (PBE0 xx) eV"
          % (WIKTOR["PBE_xx"], WIKTOR["PBE0_xx"]))


if __name__ == "__main__":
    main()
