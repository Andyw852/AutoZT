#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zt_common.py -- 器件级电-热自洽 ZT 公共库。

与 device-thermal 的分工：
  - 传热求解（FVM 红黑 SOR）**不在这里**，直接 import 公共池的 device_common
    （skill/_common/thermal/device_common.py），避免重复实现、防止两处漂移。
  - 本模块只管【电学】：读 ke 的 transport.json、按温度插值、面电导/焦耳热、
    以及"焦耳自热 -> T 升高 -> sigma(T)/S(T)/k_e(T) 变 -> 焦耳热变"的自洽迭代。

★ 电导口径（重要）：
  phono3py 之外，AMSET/BoltzTraP 的 sigma 也是【含真空的胞值】（AMSET 用胞体积；
  BoltzTraP 路径 gen_step11 已显式乘 c/t，两条路径口径不同）。
  器件只关心【面电导】 G_s = sigma * t_conv（S，per square）——只要 sigma 与
  t_conv 取自同一口径，算出的电阻、焦耳热、ZT 都与该口径无关：
      sigma_cell * h_perp  ==  sigma_2D * d
  本模块因此统一按 G_s 工作，并把实际用的 t_conv 与口径来源写进输出。
  ★ 但【kappa_e 并入热导】时必须先换成物理层值：kappa_e_2D = kappa_e_cell * t_conv/d，
  否则会把胞值（MoS2 实测 4.1 W/mK）加到物理口径的晶格 kappa（63.7 W/mK）上——混口径。
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import device_common as DC  # noqa: E402  公共池：传热 FVM + 材料库

# transport.json 里 seebeck 的单位（AMSET/BoltzTraP 均为 uV/K）
SEEBECK_UNIT_V_PER_K = 1e-6


def blank(v):
    """把 step.conf 里写成 ""/''/none/null/- 的值当成空。

    ★ 坑：stepconf 不剥引号，KE_UPSTREAM_STEP = "" 解析出来是【两个字面量引号】
    的字符串（真值），会被当成"指定了一个叫空串的步骤"，于是找不到上游文件。
    """
    if v is None:
        return None
    s = str(v).strip().strip('"').strip("'").strip()
    return None if s.lower() in ("", "none", "null", "nil", "-") else s

# transport.json 可能出现的步骤目录（AMSET 标准步 / AMSET-2D 步 / BoltzTraP 旁路）
TRANSPORT_STEPS = ("step8_amset", "step8.4_amset2d", "step8.1_boltztrap")


# ------------------------------------------------------------------ 定位/读取
def find_transport(cwd, skill="ke-dft-cpu", steps=None):
    """在 <材料>/<skill>/<step>/ 与 <材料>/<skill>/result/<step>/ 下找 transport.json。"""
    for st in (steps or TRANSPORT_STEPS):
        p = DC.find_upstream(cwd, skill, st, "transport.json")
        if p:
            return p, st
    return None, None


def _nearest(xs, x):
    return int(np.argmin([abs(float(v) - x) for v in xs]))


def _inplane(tensor, mode="avg"):
    """3x3 张量的面内代表值：avg=(xx+yy)/2，xx=只取 xx。"""
    m = np.asarray(tensor, float)
    if m.ndim != 2 or m.shape[0] < 2:
        return float(m)
    return 0.5 * (m[0, 0] + m[1, 1]) if mode == "avg" else float(m[0, 0])


def read_transport(path, doping_cm3, inplane="avg"):
    """读 transport.json，取【最接近 doping_cm3 的那条掺杂曲线】的全部温度点。

    返回 dict：T_K[] 与同长度的 sigma_S_m / seebeck_uV_K / kappa_e_W_mK（面内代表值），
    外加该掺杂下标与实际的 doping 值（便于核对没选错曲线）。
    """
    with open(path, encoding="utf-8") as f:
        j = json.load(f)
    dop = [float(x) for x in j["doping"]]
    Ts = [float(x) for x in j["temperatures"]]
    # 掺杂按 signed 值就近选（n/p 不能混）
    iD = int(np.argmin([abs(float(d) - float(doping_cm3)) for d in dop]))
    cond = np.asarray(j["conductivity"], float)          # (nD, nT, 3, 3)  S/m
    seeb = np.asarray(j["seebeck"], float)               # (nD, nT, 3, 3)  uV/K
    kel = np.asarray(j["electronic_thermal_conductivity"], float)  # (nD,nT,3,3) W/mK
    sig = [_inplane(cond[iD, i], inplane) for i in range(len(Ts))]
    se = [_inplane(seeb[iD, i], inplane) for i in range(len(Ts))]
    ke = [_inplane(kel[iD, i], inplane) for i in range(len(Ts))]
    return {
        "source": path, "doping_all_cm3": dop, "doping_index": iD,
        "doping_cm3": dop[iD], "doping_requested_cm3": float(doping_cm3),
        "T_K": Ts, "sigma_S_m": sig, "seebeck_uV_K": se, "seebeck_V_K": [v * SEEBECK_UNIT_V_PER_K for v in se],
        "kappa_e_W_mK": ke, "n_T": len(Ts),
        "seebeck_unit": "uV/K（原值）；seebeck_V_K 为 SI",
    }



def read_all_dopings(path, T_ref, inplane="avg"):
    """读 transport.json 的**全部掺杂档**，每档在 T_ref（就近温度）取面内代表值。

    用于"S 与 σ 的折中在哪一档掺杂"——器件 ZT 对掺杂极敏感（重掺 S 塌成几 uV/K，
    ZT 归零），而 transport.json 里本来就带 10 档，不该只看一档。
    返回 [{"doping_cm3","T_K","sigma_S_m","seebeck_uV_K","kappa_e_W_mK","PF_W_mK2"}, ...]
    """
    with open(path, encoding="utf-8") as f:
        j = json.load(f)
    dop = [float(x) for x in j["doping"]]
    Ts = [float(x) for x in j["temperatures"]]
    iT = _nearest(Ts, float(T_ref))
    cond = np.asarray(j["conductivity"], float)
    seeb = np.asarray(j["seebeck"], float)
    kel = np.asarray(j["electronic_thermal_conductivity"], float)
    out = []
    for iD in range(len(dop)):
        s = _inplane(cond[iD, iT], inplane)
        se = _inplane(seeb[iD, iT], inplane)
        ke = _inplane(kel[iD, iT], inplane)
        sv = se * SEEBECK_UNIT_V_PER_K
        out.append({"doping_cm3": dop[iD], "T_K": Ts[iT],
                    "sigma_S_m": s, "seebeck_uV_K": se, "kappa_e_W_mK": ke,
                    "PF_W_mK2": sv * sv * s})       # S^2 * sigma (W/m/K^2)
    return out

def interp_T(xs, ys, T):
    """线性插值；超出网格【不外推】——夹到端点并标注（避免外推给出假的高 ZT）。"""
    if not xs:
        return None, "empty"
    if T <= xs[0]:
        return float(ys[0]), "clamped_low"
    if T >= xs[-1]:
        return float(ys[-1]), "clamped_high"
    return float(np.interp(T, xs, ys)), "ok"


# ------------------------------------------------------------------ 电学量
def sheet_conductance(sigma_S_m, t_conv_m):
    """面电导 G_s [S]（per square）。sigma 与 t_conv 必须同口径。"""
    return float(sigma_S_m) * float(t_conv_m)


def sheet_resistance(sigma_S_m, t_conv_m):
    g = sheet_conductance(sigma_S_m, t_conv_m)
    return (1.0 / g) if g > 0 else float("inf")


def channel_resistance(sigma_S_m, t_conv_m, L_m, W_m):
    """沟道体电阻 R_ch = (L/W) * R_sheet = L / (G_s * W)。"""
    R = sheet_resistance(sigma_S_m, t_conv_m) * (float(L_m) / float(W_m))
    return R if np.isfinite(R) else float("inf")


def circuit(V_bias, sigma_S_m, t_conv_m, L_m, W_m, R_contact_ohm=0.0):
    """串联电路：沟道体电阻 + 接触电阻（简单串联，不改传热模型结构）。

    只有【沟道体内】耗散的功率进温度场：
        Q_channel = I^2 * R_ch / (L*W)     [W/m^2]
    接触耗散 I^2*R_contact 留在接触处——本技能的 CONTACT_BC 默认 sink（接触即热沉），
    所以那部分热直接进热沉、不加热沟道。这是保守且物理上讲得通的简化；
    若把接触热也算进沟道，dT 会被高估。两项都在输出里单列，便于核对。

    返回 dict：I / R_channel / R_total / P_channel / P_contact / P_total / Q_channel / Q_contact。
    """
    R_ch = channel_resistance(sigma_S_m, t_conv_m, L_m, W_m)
    R_ct = max(0.0, float(R_contact_ohm or 0.0))
    R_tot = R_ch + R_ct
    V = float(V_bias)
    I = (V / R_tot) if (R_tot > 0 and np.isfinite(R_tot)) else 0.0
    L, W = float(L_m), float(W_m)
    area = L * W
    P_ch, P_ct = I * I * R_ch, I * I * R_ct
    return {"I_A": I, "R_channel_ohm": R_ch, "R_contact_ohm": R_ct, "R_total_ohm": R_tot,
            "P_channel_W": P_ch, "P_contact_W": P_ct, "P_total_W": P_ch + P_ct,
            "Q_channel_W_m2": (P_ch / area) if area > 0 else 0.0,
            "Q_contact_W_m2": (P_ct / area) if area > 0 else 0.0,
            "contact_share": (P_ct / (P_ch + P_ct)) if (P_ch + P_ct) > 0 else 0.0}


def zt_material(S_V_K, sigma_S_m, kappa_total_W_mK, T_K):
    """材料 ZT = S^2 * sigma * T / kappa（比值，与电导口径无关）。"""
    k = float(kappa_total_W_mK)
    if k <= 0:
        return float("nan")
    return float(S_V_K) ** 2 * float(sigma_S_m) * float(T_K) / k


# ------------------------------------------------------------------ 自洽迭代
def self_consistent_joule(spec, props, tr, t_conv_m, network, opts=None):
    """焦耳自热 <-> 温度依赖物性 的自洽迭代。

    每一步：
      1) 用当前有效温度 T_eff 从 transport 表查 sigma/S/k_e（不外推）；
      2) Q = V^2 * G_s(T_eff) / L^2；G_s = sigma * t_conv；
      3) 把 k_e 并入沟道面内热导（k_total = k_lattice + k_e），跑 device_common FVM；
      4) T_eff_new = 沟道层单元温度均值（集总模型的标准选择）；阻尼更新直到收敛。

    network: {"k_lattice_W_mK": float, "t_conv_m": float}
    opts: {"V_bias_V","maxit","tol_K","relax","sor_*","n_grid_x","n_grid_sub","do_transient"}
    返回 (result dict, T_field ndarray, dy, lay)
    """
    import copy
    o = dict(opts or {})
    V = float(o.get("V_bias_V", 0.1))
    maxit = int(o.get("maxit", 40))
    tol = float(o.get("tol_K", 0.05))
    relax = float(o.get("relax", 0.5))
    T_amb = float(spec["T_amb_K"])
    L = float(spec["L_channel_m"])
    kl = float(network["k_lattice_W_mK"])
    # kappa_e 与 sigma 同为【含真空胞口径】；加到物理口径的晶格 kappa 上必须先换算：
    #   kappa_e_2D = kappa_e_cell * t_conv / d      （d = 沟道物理层厚）
    d_layer = float(network.get("thickness_m") or 0.0)
    ke_expand = (float(t_conv_m) / d_layer) if d_layer > 0 else 1.0
    dt_lim = float(o.get("dt_limit_K", 100.0))
    T_lim = float(o.get("T_limit_K", 450.0))

    Ts_ = tr["T_K"]
    T_eff = T_amb
    hist = []
    converged = False
    T = None
    summary = None
    dy = lay = None
    for it in range(1, maxit + 1):
        sig, f1 = interp_T(Ts_, tr["sigma_S_m"], T_eff)
        se, f2 = interp_T(Ts_, tr["seebeck_V_K"], T_eff)
        ke, f3 = interp_T(Ts_, tr["kappa_e_W_mK"], T_eff)
        cir = circuit(V, sig, t_conv_m, L, float(spec["W_device_m"]),
                      float(spec.get("R_contact_ohm") or 0.0))
        I, R, Q = cir["I_A"], cir["R_total_ohm"], cir["Q_channel_W_m2"]

        ke_2d = ke * ke_expand                    # 胞值 -> 物理层值
        p2 = copy.deepcopy(props)
        p2["channel"]["kx_W_mK"] = kl + ke_2d     # 电子热导并入面内（同为物理口径）
        p2["channel"]["ky_W_mK"] = float(props["channel"]["ky_W_mK"])
        sp = copy.deepcopy(spec)
        sp["Q_joule_W_m2"] = Q
        summary, T, dy, lay = DC.solve_device(
            sp, p2, n_grid_x=int(o.get("n_grid_x", 41)),
            n_grid_sub=int(o.get("n_grid_sub", 60)),
            do_transient=bool(o.get("do_transient", False)),
            sor_omega=float(o.get("sor_omega", 1.7)),
            sor_tol=float(o.get("sor_tol", 1e-8)),
            sor_maxit=int(o.get("sor_maxit", 200000)))
        T_ch = float(np.mean(T[:, lay == 0]))
        dT = T_ch - T_eff
        hist.append({"iter": it, "T_eff_K": T_eff, "sigma_S_m": sig,
                     "Q_joule_W_m2": Q, "I_A": I, "dT_peak_K": summary["dT_peak_K"],
                     "T_channel_mean_K": T_ch,
                     "kappa_e_cell_W_mK": ke, "kappa_e_2D_W_mK": ke * ke_expand,
                     "clamp": [f1, f2, f3]})
        T_eff = T_eff + relax * dT
        if abs(dT) < tol:
            converged = True
            break

    # 收敛后再算一次一致性量（用最终 T_eff）
    sig, _ = interp_T(Ts_, tr["sigma_S_m"], T_eff)
    se, _ = interp_T(Ts_, tr["seebeck_V_K"], T_eff)
    ke, _ = interp_T(Ts_, tr["kappa_e_W_mK"], T_eff)
    cir = circuit(V, sig, t_conv_m, L, float(spec["W_device_m"]),
                  float(spec.get("R_contact_ohm") or 0.0))
    Q, I, R = cir["Q_channel_W_m2"], cir["I_A"], cir["R_total_ohm"]
    out = {
        "ZT_SOLVE_DONE": True,
        "converged": bool(converged),
        "iterations": len(hist),
        "history": hist,
        "V_bias_V": V,
        "T_eff_final_K": T_eff,
        "T_channel_mean_K": float(np.mean(T[:, lay == 0])) if T is not None else None,
        "T_channel_max_K": float(np.max(T[:, lay == 0])) if T is not None else None,
        "dT_peak_K": summary["dT_peak_K"] if summary else None,
        "Q_joule_W_m2": Q,
        "sigma_S_m": sig, "seebeck_V_K": se,
        "kappa_e_cell_W_mK": ke, "kappa_e_W_mK": ke * ke_expand,
        "ke_cell_to_layer_factor": ke_expand,
        "kappa_lattice_W_mK": kl, "kappa_total_W_mK": kl + ke * ke_expand,
        "sheet_conductance_S": sheet_conductance(sig, t_conv_m),
        "R_electrical_ohm": R, "I_A": I,
        "R_channel_ohm": cir["R_channel_ohm"],
        "R_contact_ohm": cir["R_contact_ohm"],
        "P_W": cir["P_total_W"], "P_channel_W": cir["P_channel_W"],
        "P_contact_W": cir["P_contact_W"], "contact_share": cir["contact_share"],
        "Q_contact_W_m2": cir["Q_contact_W_m2"],
        "P_per_width_W_m": (cir["P_total_W"] / float(spec["W_device_m"])) if I else 0.0,
        # ---- 温度/击穿上限判据（超限即标记，便于选安全偏压）----
        "dt_limit_K": dt_lim, "T_limit_K": T_lim,
        "within_limit": bool((summary["dT_peak_K"] <= dt_lim)
                             and (float(np.max(T[:, lay == 0])) <= T_lim)) if T is not None else None,
        "zt_material": zt_material(se, sig, kl + ke * ke_expand, T_eff),
        "t_conv_m": float(t_conv_m),
        "sor_iterations": summary["sor_iterations"] if summary else None,
    }
    return out, T, dy, lay
