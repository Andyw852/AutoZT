#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[分支测试 adp-xi-test] ADP 形变势四口径对照 —— 只为对比，不改动 master 代码。

四个口径（都作用在 h5 的 D 上，无量纲因子）：
  ① strain          因子 1.0                    D = D_vac（横向钳制，现状）
  ② strain_poisson  因子 1 - C12/C11            D = Wiktor 那种"横向泊弛豫"口径
  ③ xi_literature   因子 Ξ_lit / D_vac          D = Kaasbjerg 的有效声学耦合 Ξ
  ④ xi_own          由我们自己的 Γ(E) 反解 Kaasbjerg Eq.(12) 得到的 Ξ_own

判据/公式（Kaasbjerg & Thygesen, PRB 87, 235312 (2013)）：
  Eq.(12)  1/tau_l = m* Ξ_l^2 k_B T / (hbar^3 rho_2D c_l^2)
  Eq.(11)  |M_ql| = Ξ_l q     （Ξ 是他们拟合 DFPT 耦合得到的**有效**常数）

附录里另用 Wiktor & Pasquarello, PRB 94, 245411 (2016) 的绝对形变势（真空口径）做对照。
"""
import argparse, json, math

# ---- 常数 ----
KB = 1.380649e-23
HBAR = 1.054571817e-34
E = 1.602176634e-19
M0 = 9.1093837e-31
AMU = 1.66053906660e-27

# Kaasbjerg 2013 Tab.I（单层 MoS2）：rho=3.1e-7 g/cm^2, c_TA=4.2e3, c_LA=6.7e3,
#   Xi_TA=1.5 eV, Xi_LA=2.4 eV, m*=0.48 me, sigma=5.41 A, e11=3.0e-11 C/m
KAAS = {"rho_kg_m2": 3.1e-7 * 10.0, "c_TA": 4.2e3, "c_LA": 6.7e3,
        "Xi_TA_eV": 1.5, "Xi_LA_eV": 2.4, "m_star": 0.48}
WIKTOR = {"PBE_xx": -6.08, "PBE_yy": -6.19, "PBE0_xx": -7.16, "VBM_PBE_xx": -1.86}


def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def rho_2d(the, area_cm2, formula=None):
    """由元胞原子数与面积算 2D 质量密度 [kg/m^2]。formula 缺省 MoS2。"""
    mm = {"Mo": 95.95, "S": 32.06, "W": 183.84, "Se": 78.97}
    comp = formula or {"Mo": 1, "S": 2}
    m_cell = sum(mm[k] * v for k, v in comp.items()) * AMU
    return m_cell / (area_cm2 * 1e-4)


def kaas_rate(Xi_eV, c, m_star):
    Xi = Xi_eV * E
    m = m_star * M0
    return m * Xi * Xi * KB * 300.0 / (HBAR ** 3 * KAAS["rho_kg_m2"] * c * c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band-edges", default="/tmp/band_edges.json")
    ap.add_argument("--2d", dest="rec2d", default="/tmp/2d_correction.json")
    ap.add_argument("--transport", default="/tmp/new_4086291_transport.json")
    ap.add_argument("--tau", default="/tmp/tau_energy.json")
    ap.add_argument("--formula", default="MoS2")
    a = ap.parse_args()

    be = load(a.band_edges)
    rec = load(a.rec2d)

    # ---- 我们的输入量 ----
    D_vac = float(be["electron"]["E1_vac_iso_eV"])
    D_vac_xx = float(be["electron"]["E1_vac_xx_eV"])
    D_vac_yy = float(be["electron"]["E1_vac_yy_eV"])
    D_core = float(be["electron"]["E1_iso_eV"])
    elas = rec["elastic_constant_raw_GPa"]
    C11, C12 = float(elas[0][0]), float(elas[0][1])
    area = float(rec["inplane_cell_area_cm2"])
    rho_our = rho_2d(rec, area, {"Mo": 1, "S": 2})
    # ★ 2D 刚度：elastic_constant_raw 是"原始 slab 值"= C_2D/c，故 C_2D = C_raw × c（N/m）
    c_cell_m = float(rec["cell_c_A"]) * 1e-10
    C2D_11 = C11 * 1e9 * c_cell_m
    C2D_66 = ((C11 - C12) / 2.0) * 1e9 * c_cell_m
    c_LA = math.sqrt(C2D_11 / rho_our)      # [N/m]/[kg/m^2] = m^2/s^2
    c_TA = math.sqrt(C2D_66 / rho_our)
    nu = C12 / C11

    print("[输入] D_vac(iso) = %.4f eV   D_vac_xx/yy = %.4f / %.4f   D_core(iso) = %.4f"
          % (D_vac, D_vac_xx, D_vac_yy, D_core))
    print("[输入] C11 = %.3f GPa  C12 = %.3f GPa  ->  nu = C12/C11 = %.4f" % (C11, C12, nu))
    print("[自算] 2D 质量密度 rho = %.3e kg/m^2   c_LA = %.0f m/s  c_TA = %.0f m/s"
          % (rho_our, c_LA, c_TA))
    print("[对照] Kaasbjerg Tab.I rho = %.3e kg/m^2  c_LA = %.0f  c_TA = %.0f"
          % (KAAS["rho_kg_m2"], KAAS["c_LA"], KAAS["c_TA"]))
    print("")

    # ---- 他们的 ADP-only 迁移率（Eq.12 实算）----
    r_lit = kaas_rate(KAAS["Xi_LA_eV"], KAAS["c_LA"], KAAS["m_star"]) + \
            kaas_rate(KAAS["Xi_TA_eV"], KAAS["c_TA"], KAAS["m_star"])
    mu_lit = E / (KAAS["m_star"] * M0 * r_lit) * 1e4
    print("[对照] Kaasbjerg Eq.(12) 300K: 1/tau = %.3e s^-1  -> ADP-only mu = %.0f cm^2/Vs"
          % (r_lit, mu_lit))
    print("")

    # ---- 我们的 ADP 实测（transport.json 的 mobility.ADP，最低 n 掺杂）----
    mu_our = None
    try:
        tr = load(a.transport)
        import numpy as np
        dop = np.array(tr["doping"], float)
        i = int(np.argmin(np.abs(dop - (-2.867054e17))))
        m = np.asarray(tr["mobility"]["ADP"], float)
        mu_our = float((m[i, 2, 0, 0] + m[i, 2, 1, 1]) / 2.0)
    except Exception as e:                                    # noqa: BLE001
        print("[WARN] 读 transport.json 失败：%s" % e)
    print("[实测] 我们的 ADP-only（最低 n, 300K, strain 口径）= %s cm^2/Vs"
          % ("%.1f" % mu_our if mu_our else "n/a"))

    # ---- 由我们的速率反解 Ξ_own（用 Kaasbjerg Eq.12 + **我们自己的** rho/c）----
    Xi_own = None
    if mu_our:
        # 有效速率由迁移率反推：mu = e/(m* Γ) -> Γ = e/(m* mu)
        r_our = E / (KAAS["m_star"] * M0 * (mu_our * 1e-4))
        # 与 Kaasbjerg Eq.(12) 同结构反解 Ξ（用**我们自己的** rho 与 c，LA/TA 按 1/c^2 加权）
        denom = KAAS["m_star"] * M0 * KB * 300.0 * (1.0 / c_LA ** 2 + 1.0 / c_TA ** 2)
        Xi_own = math.sqrt(r_our * HBAR ** 3 * rho_our / denom) / E
        r_lit_eff = math.sqrt((KAAS["Xi_LA_eV"] ** 2 / KAAS["c_LA"] ** 2 +
                               KAAS["Xi_TA_eV"] ** 2 / KAAS["c_TA"] ** 2) /
                              (1.0 / KAAS["c_LA"] ** 2 + 1.0 / KAAS["c_TA"] ** 2))
        print("[自算] 我们 ADP 的有效速率 Γ_eff = e/(m* mu) = %.3e s^-1" % r_our)
        print("[自算] 由 Γ_eff 反解 Ξ_own（我们自己的 rho/c）= %.2f eV" % Xi_own)
        print("[对照] Kaasbjerg 同口径的 c 加权有效 Ξ = %.2f eV  ->  我们/他们 = %.2f"
              % (r_lit_eff, Xi_own / r_lit_eff))
    print("")

    # ---- 四口径表 ----
    f_match = math.sqrt(mu_our / mu_lit) if mu_our else float("nan")
    fac = [("① strain        ", 1.0, D_vac),
           ("② strain_poisson", 1.0 - nu, D_vac * (1 - nu)),
           ("③ xi_literal    ", KAAS["Xi_LA_eV"] / D_vac, KAAS["Xi_LA_eV"]),
           ("④ xi_matched    ", f_match, D_vac * f_match)]
    print("=" * 96)
    print("%-18s %10s %10s %12s %12s %s" % ("口径", "因子", "D(eV)", "D/D_vac", "ADP-only μ*", "说明"))
    print("-" * 96)
    for name, f, d in fac:
        if mu_our and name.startswith("①"):
            mu = mu_our
        elif mu_our:
            mu = mu_our / (f ** 2)          # 率 ∝ D^2 → mu ∝ 1/f^2
        else:
            mu = float("nan")
        note = {"①": "现状（横向钳制，真空参考）",
                "②": "Wiktor 泊松绝对形变势口径（1-nu=%.4f）" % (1 - nu),
                "③": "直接抄 Ξ_LA=2.4 eV（**注意归一化不同**）",
                "④": "标定到 Kaasbjerg 的 ADP-only=%d（因子=sqrt(μ_our/μ_lit)）" % round(mu_lit)}[name[0]]
        print("%-18s %10.4f %10.3f %12.4f %12.1f %s" % (name, f, d, d / D_vac, mu, note))
    print("=" * 96)
    print("* ADP-only μ 由 D^2 反比缩放估计（仅量级）；②的实际值以跑完的作业为准。")
    print("")
    print("Wiktor PRB 94, 245411 (2016) MoS2 绝对形变势（真空）: "
          "a_CBM(K) = %.2f (PBE xx) / %.2f (PBE0 xx) eV" % (WIKTOR["PBE_xx"], WIKTOR["PBE0_xx"]))


if __name__ == "__main__":
    main()
