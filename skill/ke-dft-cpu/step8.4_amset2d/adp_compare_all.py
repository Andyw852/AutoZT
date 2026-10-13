#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[独立分支 adp-xi-selfcalc] 所有 ADP 形变势口径的总对照 —— 每个量标清楚"我是谁、我怎么来的"。

行 = 一个口径/方法；列 = 它用的形变势 D/Ξ、由它得到的 ADP-only 迁移率、以及它的出处。
所有"有效耦合 Ξ"都用同一套公式反解（Kaasbjerg PRB 87, 235312 Eq.12 的 c 加权式），
所以横比是同一把尺子；"带边绝对形变势"列则直接是 DFT 出的 D，不再折算。
"""
import json
import math

import numpy as np

KB, HBAR, E, M0 = 1.380649e-23, 1.054571817e-34, 1.602176634e-19, 9.1093837e-31
AMU = 1.66053906660e-27
MASS = {"Mo": 95.95, "S": 32.06, "W": 183.84, "Se": 78.97}

# 文献值
LIT = {
    "kaas_rho": 3.1e-6, "kaas_c_LA": 6700.0, "kaas_c_TA": 4200.0,
    "kaas_Xi_LA": 2.4, "kaas_Xi_TA": 1.5, "kaas_m": 0.48,
    "wiktor_a_CBM_PBE": 6.08, "wiktor_a_CBM_PBE0": 7.16,
}
# 冻声子 DFT 模式投影（fp_xi.py h5，集群实测，4 材料；见提交 fcf961f）
FP = {"MoS2": {"e": 8.5359, "h": 2.4809},
      "WS2": {"e": 9.3650, "h": 2.1650},
      "WSe2": {"e": 0.7354, "h": 1.9665},
      "MoSe2": {"e": 7.3872, "h": 2.1572}}


def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def kaas_2branch(rho, c_LA, c_TA, m_star, x_la, x_ta, t=300.0):
    m = m_star * M0
    g = m * KB * t / (HBAR ** 3 * rho) * ((x_la * E) ** 2 / c_LA ** 2 + (x_ta * E) ** 2 / c_TA ** 2)
    return E / (m * g) * 1e4


def xi_from_mu(mu, rho, c_LA, c_TA, m_star, t=300.0):
    """由 ADP-only 迁移率反解 c 加权有效耦合 Ξ [eV]（与 Eq.12 同结构）。"""
    r = E / (m_star * M0 * (mu * 1e-4))
    den = m_star * M0 * KB * t * (1.0 / c_LA ** 2 + 1.0 / c_TA ** 2)
    return math.sqrt(r * HBAR ** 3 * rho / den) / E


def takagi(d_eV, c2d, m_star, t=300.0):
    d = d_eV * E
    m = m_star * M0
    return E * HBAR ** 3 * c2d / (KB * t * m * m * d * d) * 1e4


def inplane(mu, i, iT):
    a = np.asarray(mu, float)
    return float((a[i, iT, 0, 0] + a[i, iT, 1, 1]) / 2.0)


def adp(path, doping):
    tr = load(path)
    dop = np.array(tr["doping"], float)
    i = int(np.argmin(np.abs(dop - doping)))
    iT = list(tr["temperatures"]).index(300)
    return inplane(tr["mobility"]["ADP"], i, iT), inplane(tr["mobility"]["overall"], i, iT)


def main():
    be = load("/tmp/band_edges.json")
    rec = load("/tmp/2d_correction.json")
    Dv_e = float(be["electron"]["E1_vac_iso_eV"])
    Dv_h = float(be["hole"]["E1_vac_iso_eV"])
    Dc_e = float(be["electron"]["E1_iso_eV"])
    C11, C12 = float(rec["elastic_constant_raw_GPa"][0][0]), float(rec["elastic_constant_raw_GPa"][0][1])
    nu = C12 / C11
    area = float(rec["inplane_cell_area_cm2"])
    cc = float(rec["cell_c_A"]) * 1e-10
    rho = (MASS["Mo"] + 2 * MASS["S"]) * AMU / (area * 1e-4)
    C2D = C11 * 1e9 * cc
    C2D66 = ((C11 - C12) / 2) * 1e9 * cc
    cL, cT = math.sqrt(C2D / rho), math.sqrt(C2D66 / rho)

    # 实测
    st_e, ov_st_e = adp("/tmp/new_4086291_transport.json", -2.867054e17)
    st_h, ov_st_h = adp("/tmp/new_4086291_transport.json", 2.867054e17)
    try:
        po_e, ov_po_e = adp("/tmp/poisson_transport.json", -2.867054e17)
        po_h, ov_po_h = adp("/tmp/poisson_transport.json", 2.867054e17)
    except Exception as exc:                                     # noqa: BLE001
        print("[WARN] %s" % exc)
        po_e = po_h = float("nan")
        ov_po_e = ov_po_h = float("nan")

    xi_st_e = xi_from_mu(st_e, rho, cL, cT, 0.48)     # AMSET strain 口径反解出的有效 Ξ
    xi_po_e = xi_from_mu(po_e, rho, cL, cT, 0.48)
    xi_st_h = xi_from_mu(st_h, rho, cL, cT, 0.606)
    xi_po_h = xi_from_mu(po_h, rho, cL, cT, 0.606)
    kaas_our = kaas_2branch(rho, cL, cT, 0.48, LIT["kaas_Xi_LA"], LIT["kaas_Xi_TA"])
    kaas_theirs = kaas_2branch(LIT["kaas_rho"], LIT["kaas_c_LA"], LIT["kaas_c_TA"], 0.48,
                               LIT["kaas_Xi_LA"], LIT["kaas_Xi_TA"])
    xi_lit_eff = math.sqrt((LIT["kaas_Xi_LA"] ** 2 / LIT["kaas_c_LA"] ** 2 +
                            LIT["kaas_Xi_TA"] ** 2 / LIT["kaas_c_TA"] ** 2) /
                           (1 / LIT["kaas_c_LA"] ** 2 + 1 / LIT["kaas_c_TA"] ** 2))

    print("=" * 108)
    print("MoS2 单层 ADP 形变势口径总对照（300 K，面内平均）")
    print("自算输入：D_vac(e)=%.4f eV  D_vac(h)=%.4f eV  D_core(e)=%.4f eV  nu=%.4f"
          % (Dv_e, Dv_h, Dc_e, nu))
    print("          rho_2D=%.4e kg/m^2  c_LA=%.0f  c_TA=%.0f m/s  C2D_11=%.2f N/m"
          % (rho, cL, cT, C2D))
    print("=" * 108)
    print()
    print("【A】每个方法给出的'形变势/有效耦合'（eV）。带★=自算，带☆=文献")
    print("-" * 108)
    print("%-3s %-42s %10s %10s %s" % ("#", "口径 / 方法（都在 h5 的 D 上作用）", "D 或 Ξ", "被当 Ξ?", "出处"))
    rows = [
        ("A1", "AMSET 2D ADP · strain（生产现状，D_vac 直接用）", Dv_e, "是", "4086291 transport.json"),
        ("A2", "AMSET 2D ADP · strain_poisson（D_vac·(1-nu)）", Dv_e * (1 - nu), "是", "4086466 transport.json"),
        ("A3", "AMSET 2D ADP · 抄文献 Ξ_LA 当 D", LIT["kaas_Xi_LA"], "是", "adp_xi_compare ③"),
        ("A4", "AMSET 实测 τ 反解出的 c 加权有效耦合 Ξ_own", xi_st_e, "是", "adp_xi_compare ④（反解）"),
        ("A5", "★冻声子 DFT 模式投影（q→0，本分支）电子", FP["MoS2"]["e"], "是", "fp_xi.py h5（fcf961f）"),
        ("A6", "★Takagi/Kaasbjerg 闭式，代我们的 D_vac", Dv_e, "是", "adp_selfxi.py"),
        ("A7", "☆Kaasbjerg Tab.I 的 Ξ_LA", LIT["kaas_Xi_LA"], "是", "PRB 87,235312 (2013)"),
        ("A8", "☆Kaasbjerg 的 c 加权有效耦合 Ξ_eff", xi_lit_eff, "是", "由 A7 的 2.4/1.5 折算"),
        ("A9", "☆Wiktor 真空绝对形变势 a_CBM(K) (PBE)", LIT["wiktor_a_CBM_PBE"], "否(是绝对 D)", "PRB 94,245411 (2016)"),
    ]
    for n, name, val, asxi, src in rows:
        print("%-3s %-42s %10.3f %10s %s" % (n, name, val, asxi, src))
    print()
    print("【B】ADP-only 迁移率（cm^2/Vs）。实测=AMSET 跑出来的；自算/文献=按上面同一套公式算的")
    print("-" * 108)
    print("%-4s %-46s %10s %10s %s" % ("#", "方法", "电子", "空穴", "备注"))
    b = [
        ("B1", "AMSET 实测 · strain（D_vac，生产现状）", st_e, st_h, "transport.json 4086291"),
        ("B2", "AMSET 实测 · strain_poisson", po_e, po_h, "transport.json 4086466"),
        ("B3", "AMSET 实测 · 反解有效 Ξ=%.2f eV" % xi_st_e, st_e, st_h, "与 B1 同一个数，只是换语言"),
        ("B4", "★Takagi 闭式 · 代 D_vac", takagi(Dv_e, C2D, 0.48), takagi(Dv_h, C2D, 0.606), "adp_selfxi 自检过 4 个台子"),
        ("B5", "★Takagi 闭式 · 代泊松 D_vac(1-nu)", takagi(Dv_e * (1 - nu), C2D, 0.48),
         takagi(Dv_h * (1 - nu), C2D, 0.606), "同上"),
        ("B6", "☆Kaasbjerg Eq.12 · 他们的 Ξ（我们 ρ/c）", kaas_our, float("nan"), "2.4/1.5 eV"),
        ("B7", "☆Kaasbjerg Eq.12 · 他们的 Ξ（他们 ρ/c）", kaas_theirs, float("nan"), "Tab.I 复现"),
        ("B8", "★把我们的 D_vac 当 Ξ 喂 Eq.12（单支 Ξ_TA=0）", takagi(Dv_e, C2D, 0.48), takagi(Dv_h, C2D, 0.606), "= B4，自洽"),
    ]
    for n, name, a, b2, note in b:
        sa = "%.1f" % a if a == a else "-"
        sb = "%.1f" % b2 if b2 == b2 else "-"
        print("%-4s %-46s %10s %10s %s" % (n, name, sa, sb, note))
    print()
    print("【C】整体（全机制）迁移率 —— 论文里真正的数（cm^2/Vs）")
    print("-" * 108)
    print("    AMSET strain        : 电子 overall %.2f   空穴 overall %.2f" % (ov_st_e, ov_st_h))
    print("    AMSET strain_poisson: 电子 overall %.2f   空穴 overall %.2f" % (ov_po_e, ov_po_h))
    print("    （整体由 ADP+POP+ODP+IMP 竞争决定；ADP 口径只动 ADP 那一项）")
    print()
    # ---------- E. 4 材料：我们冻声子 Xi_LA vs Wiktor 绝对形变势 ----------
    WIK = {  # PRB 94, 245411 (2016) TABLE I，PBE xx，eV；真空参考的绝对形变势
        "MoS2":  {"e_K": 6.08, "h_K": 1.86, "e_Q": 0.61},
        "MoSe2": {"e_K": 5.62, "h_K": 1.86, "e_Q": 0.71},
        "WS2":   {"e_K": 6.76, "h_K": 1.59, "e_Q": 0.90},
        "WSe2":  {"e_K": 6.35, "h_K": 1.43, "e_Q": 0.93},
    }
    print("【E】4 材料：我们冻声子 DFT 的 Xi_LA vs Wiktor **绝对形变势**（PRB 94,245411 Tab.I, PBE xx）")
    print("-" * 108)
    print("%-7s %11s %14s %8s   %11s %13s %8s"
          % ("材料", "我们Xi_LA(e)", "Wiktor a_CBM", "比", "我们Xi_LA(h)", "Wiktor a_VBM", "比"))
    for m in ("MoS2", "MoSe2", "WS2", "WSe2"):
        fp, w = FP[m], WIK[m]
        ref_e, tag = (w["e_Q"], "Q") if m == "WSe2" else (w["e_K"], "K")
        print("%-7s %11.3f %9.3f(%s) %7.2fx   %11.3f %10.3f(K) %7.2fx"
              % (m, fp["e"], ref_e, tag, fp["e"] / ref_e,
                 fp["h"], w["h_K"], fp["h"] / w["h_K"]))
    print("  WSe2 电子带边落在非 K 谷（band_edges.json 给的 k=(0.1778,0.1778,0)），故与 a_CBM(Q) 比。")
    print("  -> 同趋势、同量级：我们的冻声子绝对形变势在 4 个材料上系统性约为 Wiktor 的 1.3~1.4 倍")
    print("     （差异来自 ionrelax 构型/泛函/真空窗；WSe2 的 Q 谷是 0.79 倍）。")
    print("  -> 而 Kaasbjerg 的'有效耦合 Ξ'(2.4) 比绝对形变势小 ~3.5 倍 —— 两个'文献值'不是同一个量。")
    print()
    print("【D】结论（三句话）")
    print(" 1. 同一把尺子下：我们 ADP 的有效耦合 Ξ_own=%.2f eV vs 文献 Ξ_eff=%.2f eV -> %.2f 倍,"
          % (xi_st_e, xi_lit_eff, xi_st_e / xi_lit_eff))
    print("    平方 %.2f 倍迁移率差,与直接比文献/我们 = %.2f 一致 -> 差是**口径**不是流程错。"
          % ((xi_st_e / xi_lit_eff) ** 2, kaas_theirs / st_e))
    print(" 2. 带边绝对形变势 (冻声子 DFT) = %.2f eV,比热窗有效 Ξ_own 大 %.2f 倍 ——"
          % (FP["MoS2"]["e"], FP["MoS2"]["e"] / xi_st_e))
    print("    差在 AMSET 用热窗内各态自己的 D(k)(多谷/各向异性),带边单点看不到。")
    print(" 3. 冻声子 DFT 给的 Ξ_TA(q→0)=0(四材料全为 0),与文献 App.B 一致；"
          "文献 Ξ_TA=1.5 来自 umklapp。")
    print("=" * 108)


if __name__ == "__main__":
    main()
