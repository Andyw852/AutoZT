#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zt_common.py —— ZT 汇总的纯逻辑（无 IO、无 matplotlib），供 gen_step20_zt.py
与离线单测共用。

物理口径（与上游技能保持一致，逐条对齐 ke-dft-cpu/step8.3_output）：
  ZT = S²σT / (κ_e + κ_L)          S 取 μV/K → 乘 1e-6 转 V/K
  · S / σ / κ_e  来自 step8_amset/transport.json（AMSET 口径：S=μV/K、σ=S/m、κ=W/m/K）
  · κ_L          来自 step6_kappa/kappa_summary.json（phono3py BTE，W/m/K，元胞口径）
  · 张量各向同性化：3D → (xx+yy+zz)/3；2D → (xx+yy)/2（剔除真空 zz，理由同 step8.3）
  · κ_L 与 σ/κ_e 必须【同口径】：本模块只用元胞口径的 kappa_xx_yy_zz，
    不用 kl 的 kappa_2d_normalized_*（那是按厚度重标度的另一套口径）。
  · κ_L(T) 线性插值到 AMSET 的温度网格；**不外推**（超出 kl 温区返回 None）。

本文件不 import numpy，纯 Python 列表运算，任何环境都能跑。
"""

# AMSET 的 Seebeck 单位是 μV/K（amset/core/transport.py 里 seebeck *= 1e6）
UV_PER_V = 1e-6


# ===== 基础 =====
def unwrap(v):
    """monty 的 {@class:array,data:[...]} → 原始 list；{n:[...],p:[...]} → n+p 拼接。"""
    if isinstance(v, dict):
        if v.get("@class") == "array":
            return unwrap(v.get("data"))
        if "n" in v or "p" in v:
            return list(unwrap(v.get("n")) or []) + list(unwrap(v.get("p")) or [])
        return None
    return v


def diag_avg(t33, is_2d=False):
    """张量各向同性化。3D: (xx+yy+zz)/3；2D: (xx+yy)/2。与 step8.3 逐字一致。"""
    if is_2d:
        return (t33[0][0] + t33[1][1]) / 2.0
    return (t33[0][0] + t33[1][1] + t33[2][2]) / 3.0


def diag3(t33):
    return (t33[0][0], t33[1][1], t33[2][2])


def interp_T(temps, vals, T):
    """按温度线性插值；超出范围返回 None（不外推）。temps/vals 等长。"""
    if not temps or not vals or len(temps) != len(vals):
        return None
    pairs = sorted((float(t), float(v)) for t, v in zip(temps, vals)
                   if v is not None and _finite(v))
    if not pairs:
        return None
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    if T < xs[0] - 1e-9 or T > xs[-1] + 1e-9:
        return None
    for i in range(len(xs) - 1):
        if xs[i] - 1e-9 <= T <= xs[i + 1] + 1e-9:
            if abs(xs[i + 1] - xs[i]) < 1e-12:
                return ys[i]
            w = (T - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] * (1 - w) + ys[i + 1] * w
    return ys[-1]


def _finite(x):
    try:
        f = float(x)
    except (TypeError, ValueError):
        return False
    return f == f and abs(f) != float("inf")


def zt_value(S_uVK, sigma_Sm, kappa_e, kappa_L, T):
    """ZT = S²σT/(κ_e+κ_L)。任一缺项/非有限/分母<=0 → None。"""
    if None in (S_uVK, sigma_Sm, kappa_e, kappa_L, T):
        return None
    if not all(_finite(x) for x in (S_uVK, sigma_Sm, kappa_e, kappa_L, T)):
        return None
    denom = float(kappa_e) + float(kappa_L)
    if denom <= 0:
        return None
    return (float(S_uVK) * UV_PER_V) ** 2 * float(sigma_Sm) * float(T) / denom


def power_factor(S_uVK, sigma_Sm):
    """PF = S²σ  [W/m/K²]；缺项返回 None。"""
    if S_uVK is None or sigma_Sm is None:
        return None
    if not all(_finite(x) for x in (S_uVK, sigma_Sm)):
        return None
    return (float(S_uVK) * UV_PER_V) ** 2 * float(sigma_Sm)


def carrier_type(doping_cm3):
    """AMSET 约定：负 doping = n 型，正 = p 型。"""
    return "n" if float(doping_cm3) < 0 else "p"


# ===== 读入 =====
def load_transport(path):
    """读 AMSET transport.json → {doping, temperatures, sigma, seebeck, kappa_e, ...}。

    sigma/seebeck/kappa_e 的索引是 [i_doping][i_T][3][3]（解析失败抛 ValueError）。
    """
    import json
    with open(path, encoding="utf-8") as fh:
        j = json.load(fh)
    doping = unwrap(j.get("doping"))
    temps = unwrap(j.get("temperatures"))
    if not doping or not temps:
        raise ValueError("transport.json 缺 doping/temperatures")
    # 掺杂值只做显示层归一（AMSET 写出的是 -9.999999999999998e19 这类浮点尾巴），
    # 10 位有效数字远高于物理精度，不影响任何计算。
    out = {"doping": [float("%.10g" % x) for x in doping],
           "temperatures": [float(x) for x in temps]}
    for key, name in (("conductivity", "sigma"),
                      ("seebeck", "seebeck"),
                      ("electronic_thermal_conductivity", "kappa_e")):
        q = unwrap(j.get(key))
        if q is None:
            raise ValueError("transport.json 缺 %s" % key)
        if len(q) != len(out["doping"]):
            raise ValueError("%s 的掺杂档数 %d != doping 档数 %d"
                             % (name, len(q), len(out["doping"])))
        out[name] = q
    return out


def load_kappa(path):
    """读 kl 的 kappa_summary.json → {temperatures, k_xx, k_yy, k_zz, k_scalar(3D/2D 都存原始)}。

    只用元胞口径 kappa_xx_yy_zz；缺该字段抛 ValueError。
    """
    import json
    with open(path, encoding="utf-8") as fh:
        j = json.load(fh)
    temps = j.get("temperatures")
    raw = j.get("kappa_xx_yy_zz")
    single_T = False
    if not temps or not raw or len(temps) != len(raw):
        # kl-mlff-* 的老格式只写 300 K 单值（kappa_300K_xx_yy_zz），没有惰温数组。
        # 这里降级成"单温度点"网格并打标：interp 模式下就只剩 300 K 可算，
        # 想要全 T 曲线请用 KTEMP_MODE=const300（本技能文档里写明的那条路）。
        k300 = j.get("kappa_300K_xx_yy_zz")
        if k300 and len(k300) >= 3:
            temps, raw, single_T = [300.0], [list(k300[:3])], True
        else:
            raise ValueError("kappa_summary.json 缺 kappa_xx_yy_zz/temperatures，"
                             "也没有 kappa_300K_xx_yy_zz（单值兜底）")
    out = {"temperatures": [float(t) for t in temps],
           "k_xx": [float(r[0]) for r in raw],
           "k_yy": [float(r[1]) for r in raw],
           "k_zz": [float(r[2]) for r in raw],
           "kappa_300K_xx_yy_zz": j.get("kappa_300K_xx_yy_zz"),
           "kappa_inplane_xx_yy": j.get("kappa_inplane_xx_yy"),
           "Lz_ang": j.get("Lz_ang"), "thickness_d_ang": j.get("thickness_d_ang"),
           "thickness_convention": j.get("thickness_convention"),
           "kappa_2d_normalized_300K": j.get("kappa_2d_normalized_300K"),
           "single_T_fallback": single_T,
           "dim": j.get("dim"), "raw_keys": sorted(j.keys())}
    return out


def kappa_scalar(kl, is_2d):
    """κ_L(T) 各向同性化：2D 取面内 (xx+yy)/2，3D 取 (xx+yy+zz)/3。"""
    if is_2d:
        return [(a + b) / 2.0 for a, b in zip(kl["k_xx"], kl["k_yy"])]
    return [(a + b + c) / 3.0
            for a, b, c in zip(kl["k_xx"], kl["k_yy"], kl["k_zz"])]


# ===== 组装 ZT 栅格 =====
def build_grid(tr, kl, is_2d=False, ktemp_mode="interp", const_T=300.0):
    """把电子输运栅格与 κ_L 合成 ZT 栅格。

    返回 dict：
      temperatures, doping(带符号, cm^-3), kinds(n/p),
      rows: 每个掺杂档一行 dict（S/sigma/kappa_e/kappa_L/ZT 的列表，按温度对齐）
      grid_zt: [[ZT(掺杂 i, 温度 j)]]
      notes: 口径说明行
    """
    temps = [float(t) for t in tr["temperatures"]]
    notes = []
    ksc = kappa_scalar(kl, is_2d)
    if ktemp_mode == "const300":
        k300 = interp_T(kl["temperatures"], ksc, const_T)
        kl_used = [k300] * len(temps)
        notes.append("κ_L 取 %g K 单值 %.4g W/m/K（KTEMP_MODE=const300）"
                     % (const_T, k300 if k300 is not None else float("nan")))
    else:
        kl_used = [interp_T(kl["temperatures"], ksc, T) for T in temps]
        miss = [T for T, v in zip(temps, kl_used) if v is None]
        notes.append("κ_L 按温度线性插值到电子段温度网格（KTEMP_MODE=interp，不外推）")
        if miss:
            notes.append("κ_L 温区 %.0f~%.0f K 覆盖不到的温度点 ZT 记为 null：%s"
                         % (min(kl["temperatures"]), max(kl["temperatures"]),
                            ", ".join("%g" % t for t in miss)))
    notes.append("κ_L 用元胞口径 kappa_xx_yy_zz（与 AMSET 的 σ/κ_e 同口径；"
                 "不用 kappa_2d_normalized_*）")
    if kl.get("single_T_fallback"):
        notes.append("★ κ_L 来源只有 300 K 单值（kl-mlff-* 老格式，无惰温数组）："
                     "interp 模式下只有 300 K 能算 ZT；要全温曲线请用 KTEMP_MODE=const300。")
    notes.append("横坐标 doping 是 AMSET 的【输入掺杂网格】（名义 cm^-3，负=n 型）——"
                 "transport.json 没有自洽载流子浓度字段；要真实 n 请用 "
                 "step8.1_boltztrap/boltztrap_crta.json 的 carrier_conc_cm-3 另行核对。")
    notes.append("κ_e 为 AMSET electronic_thermal_conductivity（对全部能带积分并扣零电流项，"
                 "近本征/窄带隙时天然含双极项；transport.json 未把双极项单独拆出）。")
    notes.append("张量约化：%s" % ("2D 面内 (xx+yy)/2（剔除真空 zz）" if is_2d
                                  else "3D 对角 (xx+yy+zz)/3"))

    rows, grid = [], []
    for i, dop in enumerate(tr["doping"]):
        r = {"doping_cm-3": dop, "type": carrier_type(dop),
             "seebeck_uV/K": [], "sigma_S/m": [], "kappa_e_W/mK": [],
             "kappa_L_W/mK": [], "kappa_tot_W/mK": [],
             "PF_W/mK2": [], "ZT": [],
             "seebeck_xx": [], "seebeck_yy": [], "sigma_xx": [], "sigma_yy": [],
             "kappa_e_xx": [], "kappa_e_yy": []}
        zrow = []
        for jT, T in enumerate(temps):
            se = diag_avg(tr["seebeck"][i][jT], is_2d)
            sg = diag_avg(tr["sigma"][i][jT], is_2d)
            ke = diag_avg(tr["kappa_e"][i][jT], is_2d)
            klv = kl_used[jT]
            z = zt_value(se, sg, ke, klv, T)
            r["seebeck_uV/K"].append(se)
            r["sigma_S/m"].append(sg)
            r["kappa_e_W/mK"].append(ke)
            r["kappa_L_W/mK"].append(klv)
            r["kappa_tot_W/mK"].append(None if klv is None else ke + klv)
            r["PF_W/mK2"].append(power_factor(se, sg))
            r["ZT"].append(z)
            zrow.append(z)
            r["seebeck_xx"].append(tr["seebeck"][i][jT][0][0])
            r["seebeck_yy"].append(tr["seebeck"][i][jT][1][1])
            r["sigma_xx"].append(tr["sigma"][i][jT][0][0])
            r["sigma_yy"].append(tr["sigma"][i][jT][1][1])
            r["kappa_e_xx"].append(tr["kappa_e"][i][jT][0][0])
            r["kappa_e_yy"].append(tr["kappa_e"][i][jT][1][1])
        rows.append(r)
        grid.append(zrow)
    return {"temperatures": temps, "doping": list(tr["doping"]),
            "kappa_L_temperatures": list(kl["temperatures"]),
            "kappa_L_scalar": ksc, "kappa_L_used": kl_used,
            "rows": rows, "grid_zt": grid, "notes": notes}


def peak_zt(grid):
    """每种载流子类型的峰值 ZT：{'n': {...}, 'p': {...}}，找不到给 None。"""
    out = {}
    for kind in ("n", "p"):
        best = None
        for i, dop in enumerate(grid["doping"]):
            if carrier_type(dop) != kind:
                continue
            for jT, T in enumerate(grid["temperatures"]):
                z = grid["grid_zt"][i][jT]
                if z is None or not _finite(z):
                    continue
                if best is None or z > best["ZT"]:
                    best = {"ZT": z, "T_K": T, "doping_cm-3": dop,
                            "kappa_L_W/mK": grid["rows"][i]["kappa_L_W/mK"][jT]}
            # 同时记该类型的“最优掺杂档”（按其峰值）
        out[kind] = best
    return out


def zt_at_doping(grid, kind, target_cm3):
    """取某类型里 |doping| 最接近 target 的档（用于出 ZT-T 曲线）。"""
    cand = [i for i, d in enumerate(grid["doping"])
            if carrier_type(d) == kind]
    if not cand:
        return None
    i = min(cand, key=lambda k: abs(abs(grid["doping"][k]) - abs(target_cm3)))
    return {"index": i, "doping_cm-3": grid["doping"][i],
            "T": grid["temperatures"], "ZT": grid["grid_zt"][i]}


def _fmt(x, nd=4):
    if x is None:
        return "-"
    try:
        f = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not _finite(f):
        return "nan"
    if f != 0 and (abs(f) >= 1e5 or abs(f) < 1e-3):
        return "%.3e" % f
    return ("%%.%df" % nd) % f


def text_report(grid, peaks, meta):
    """人读的汇总（写 zt_summary.txt 用），纯字符串。"""
    L = []
    L.append("=" * 78)
    L.append("ZT 汇总  material=%s  skill=zt-dft-cpu  step=step20_zt" % meta.get("material"))
    L.append("ZT = S^2 * sigma * T / (kappa_e + kappa_L)     "
             "[S: uV/K -> V/K, sigma: S/m, kappa: W/m/K]")
    L.append("电子段 %s    晶格段 %s" % (meta.get("transport_src"), meta.get("kappa_src")))
    L.append("DIM=%s（%s）" % (meta.get("dim"), meta.get("dim_note")))
    L.append("=" * 78)
    for n in grid["notes"]:
        L.append("  · " + n)
    L.append("")
    L.append("-- 峰值 ZT ------------------------------------------------------------")
    for kind in ("n", "p"):
        p = peaks.get(kind)
        if not p:
            L.append("  %s 型：无可用点（κ_L 温区与电子段不重叠？）" % kind)
            continue
        L.append("  %s 型：ZT_max = %s  @ T = %s K, n = %s cm^-3 (κ_L = %s W/m/K)"
                 % (kind, _fmt(p["ZT"]), _fmt(p["T_K"], 1), _fmt(p["doping_cm-3"], 3),
                    _fmt(p["kappa_L_W/mK"])))
    L.append("")
    for i, r in enumerate(grid["rows"]):
        L.append("-- %s 型  %s cm^-3 -------------------------------------------------"
                 % (r["type"], _fmt(r["doping_cm-3"], 3)))
        if grid["rows"][i]["ZT"] and all(v is None for v in r["ZT"]):
            L.append("   （该档所有温度都缺 κ_L，ZT 不计算）")
            continue
        hdr = "   %-8s %-11s %-12s %-12s %-12s %-11s %-9s" % (
            "T(K)", "S(uV/K)", "sigma(S/m)", "k_e(W/mK)", "k_L(W/mK)",
            "PF(W/mK2)", "ZT")
        L.append(hdr)
        for jT, T in enumerate(grid["temperatures"]):
            L.append("   %-8s %-11s %-12s %-12s %-12s %-11s %-9s" % (
                _fmt(T, 0), _fmt(r["seebeck_uV/K"][jT]), _fmt(r["sigma_S/m"][jT]),
                _fmt(r["kappa_e_W/mK"][jT]), _fmt(r["kappa_L_W/mK"][jT]),
                _fmt(r["PF_W/mK2"][jT]), _fmt(r["ZT"][jT], 3)))
        L.append("")
    return "\n".join(L) + "\n"
