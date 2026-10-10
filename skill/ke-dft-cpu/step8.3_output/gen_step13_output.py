#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step13_output.py —— 三方电子输运对比（step8.3_output）。

run:gen 步骤：登录节点直接跑、秒级。读取
  · step8_amset/transport.json        (amset：σ/S/κ_e/迁移率，多温度多掺杂，绝对值)
  · step8.1_boltztrap/boltztrap_crta.json (BoltzTraP2 CRTA：σ/τ、S、κ_e/τ、Lorenz)
  · step8.2_dpt/dpt_result.json       (DPT：ADP 声学支迁移率，300K)
在 **300 K 附近** 并排出 σ / S / κ_e / Lorenz / 迁移率 的对比（表 + 图）。

可比性说明（写进产物里）：
  · S 和 Lorenz 与弛豫时间无关 → amset 与 BoltzTraP2 可**直接对比**（核心交叉验证）。
  · σ、κ_e：amset 是绝对值，BoltzTraP2 是 per-τ（除以弛豫时间）→ 只比**形状/趋势**。
  · 迁移率：amset(overall) 随掺杂变；DPT 是单值(仅 ADP)、随 T 按 1/T 标度 → 画成水平线。
  · 张量各向同性化：3D 取对角平均 (xx+yy+zz)/3；2D 取面内 (xx+yy)/2（剔除真空 zz），
    与 step8.1_boltztrap 一致。剔除 zz 的理由：真空方向群速度≈0，σ_zz/κ_zz≈0 是稀释项，
    而 S_zz=α_zz/σ_zz 是 0/0 病态比值——把它塞进平均会不可控地污染 S。
    2D 面内各向异性(xx≠yy)时，另出 S_xx/S_yy、σ_xx/σ_yy、κ_xx/κ_yy 分方向列，ZT 分方向算。

产物（done_marker）：comparison_300K.png + comparison_300K.csv + comparison_summary.txt
[0020] 另出迁移率随温度表（尽力而为，不影响完成标记）：comparison_vs_T.csv / .json + mobility_vs_T.png ——
  每个温度（AMSET 的温度网格，默认 100–900 K）× 载流子（取最低掺杂）× 方向（2D x / y / 面内平均）：
  DPT（由 dpt_result.json 的 temperature_K 换算，2D ∝1/T、3D ∝T^-1.5）、AMSET 只 ADP（intrinsic_ADP.json，
  没有就用 transport.json 的 ADP 单机制）、AMSET 本征（intrinsic_transport.json，去 IMP）、AMSET 全机制（overall）。
缺任一输入都不崩，能画多少画多少。
"""
import json
import math
import os
import sys
from pathlib import Path

# [SKILL_REV] 版本戳：写进 comparison_summary.txt。每次改本脚本逻辑后更新。
#   陈旧副本已咬人三次（step12/step9b 盖戳后，step13 是最后一个没盖的），
#   这里盖戳便于从产物反查到底跑的是哪份 skill 副本。
_SKILL_REV = "2026-10-10-0026"
#   bump 记录：2026-10-10-0026 —— 3D 的 BT2 分方向加 zz（σ/S/PF/κe_WF/ZT），3D 的平均 ZT 用三个方向的 κL；
#   恢复 0024 的 [DPT 单谷] / [AMSET 机制] 汇总行（本机合并 0030+V181 时调用被冲掉，函数还在）。
#   （本机 5911b60 的 "2026-10-09-0025" 是 0030+0025 的合并版，与纯 0025 不是同一份内容。）
#   bump 记录：2026-10-07-v181 —— 汇总里写明泛函（csv 的 functional 列 + summary 的 [口径] 行）。
#   bump 记录：2026-10-09-0021 —— DPT 的 μ/τ 换到 TARGET_T（S8.2 的 TEMPERATURE_K 可改了）；
#   2D 只有真实重叠、没有 unity 对照不再判红（2026-09-26 起真实重叠是出厂默认）。
#   bump 记录：2026-10-09-mu-vs-T —— [0020] 增加 μ(T) 三值表 comparison_vs_T.csv/json + mobility_vs_T.png。
#   bump 记录：2026-09-16-amset2d —— 增加 step8.4_amset2d（二维散射核）对比列；
#   仅在目录存在时读取，**不写进 needs**，以免 3D 项目因缺这个步骤而卡住。

OUTDIR_NAME = "step8.3_output"
AMSET_DIR = "step8_amset"
# patch_amset2d：二维散射核（插件版）的结果目录。存在就一起对比，不存在就跳过。
AMSET2D_DIR = "step8.4_amset2d"
# ---- 掺杂可解读范围（2026-09-16 二次修订）----------------------------------
# 判据必须**与原胞大小无关**，所以不用面密度、用【每原胞载流子数】：
#     n_cell = n_2D × A_cell（2D）或 n_3D × V_cell（3D）
# 理由：同样 1e14 cm^-2，在 MoS2（A = 8.7 Å²）是每胞 0.087 个载流子，合理；
# 在富勒烯网络 / Mg2C60 这类大胞（A 大一个量级以上）就是每胞好几个 ——
# 早超出刚性能带近似（也超出实验可达），只按面密度判会**漏报**。
# 上限取每原胞约 0.1 个载流子；超过的档位只在表里留数、标 NO，不作定量解读、
# 不进论文图。物理依据见 VERIFICATION §15.10。
NCELL_INTERPRET_MAX = 0.1         # 每原胞载流子数（2D 与 3D 同一判据）
N2D_INTERPRET_MAX_REF = 1.0e14    # cm^-2：仅作 2D 输出的参考对照（MoS2 量级），不参与判据


def _cell_geometry(cwd):
    """(A_cell [cm²], V_cell [cm³])：从 step3_uniform 的 CONTCAR/POSCAR 读晶胞。

    2D 用 A = |a×b|（面内），3D 用 V = |a·(b×c)|。读不到返回 (None, None)——
    那种情况下不做每胞判据，只保留面密度列。
    """
    import glob as _glob
    cands = []
    for _d in ("step3_uniform", "step1_opt", "step1_std_opt", "."):
        cands += sorted(_glob.glob(str(Path(cwd) / _d / "CONTCAR")))
        cands += sorted(_glob.glob(str(Path(cwd) / _d / "POSCAR")))
    for p in cands:
        try:
            ln = Path(p).read_text(errors="ignore").splitlines()
            s = float(ln[1].split()[0])
            L = [[float(x) * s for x in ln[2 + i].split()[:3]] for i in range(3)]
        except (OSError, IndexError, ValueError):
            continue
        a, b, c = L
        def _cross(u, v):
            return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0]]
        def _norm(u):
            return sum(x * x for x in u) ** 0.5
        cr = _cross(a, b)
        A = _norm(cr) * 1e-16                                # Å² -> cm²
        # ★ 体积必须是三重积 c·(a×b)。初版误写成 a·(a×b) 恒为 0（自检当场抓到：
        #   V_cell 打印成 0.0000 Å³），那会让 3D 的每胞判据永远"通过"。
        V = abs(sum(c[i] * cr[i] for i in range(3))) * 1e-24  # Å³ -> cm³
        return A, V
    return None, None


def _deform_geometry(cwd):
    """形变势构型口径（ionrelax/clamped/mixed）：step7b 的 band_edges.json 里记着实际用
    哪一套。汇总表要输出这一列 —— 回退到 clamped 时不同项目口径不一致，必须看得见。"""
    for rel in ("step7b_deform_read", "step7_deform"):
        p = Path(cwd) / rel / "band_edges.json"
        if p.is_file():
            try:
                v = (_load_json(p) or {}).get("deform_geometry")
            except Exception:                                  # noqa: BLE001
                v = None
            if v:
                return str(v)
    return None


# [V181] 泛函标签。本步 gen_need 为空、不依赖 ke_common 模块：先取 S8.2/S8.4 已写进结果的，
#   都没有再直接读 step1 的 workflow_method.txt（FUNC=），标签与 ke_common.FUNC_LABEL 保持一致。
_FUNC_LABEL = {"pbe": "PBE（GGA=PE，无色散修正）",
               "pbesol": "PBEsol（GGA=PS，无色散修正）",
               "pbe-d3": "PBE+D3(BJ)（GGA=PE，IVDW=12）"}


def _functional(cwd):
    """返回 (标签, 来源)；读不到返回 (None, None)。"""
    for rel, key in (("step8.2_dpt/dpt_result.json", "functional"),
                     ("step8.4_amset2d/2d_correction.json", "functional")):
        f = (_load_json(Path(cwd) / rel) or {}).get(key) if (Path(cwd) / rel).is_file() else None
        if isinstance(f, dict) and f.get("func"):
            return f.get("label") or f["func"], rel
    for d in ("step1_opt", "step1"):
        mf = Path(cwd) / d / "workflow_method.txt"
        if mf.is_file():
            for ln in mf.read_text(errors="ignore").splitlines():
                if ln.strip().upper().startswith("FUNC="):
                    v = ln.split("=", 1)[1].strip().lower()
                    return _FUNC_LABEL.get(v, v), "%s/workflow_method.txt" % d
    return None, None


def _areal_factor_cm(cwd):
    """2D 面密度换算因子 c（cm）：优先 amset2d 的 2d_correction.json，再 step8_amset。"""
    for rel in (AMSET2D_DIR, AMSET_DIR):
        p = Path(cwd) / rel / "2d_correction.json"
        if p.is_file():
            try:
                v = (_load_json(p) or {}).get("areal_density_factor_cm")
            except Exception:                                  # noqa: BLE001
                v = None
            if v:
                return float(v)
    return None
BT2_DIR   = "step8.1_boltztrap"
DPT_DIR   = "step8.2_dpt"
TARGET_T  = 300.0
# 表里报的代表性载流子浓度（cm^-3，负=n型电子，正=p型空穴）
TARGET_N = [-1e20, -1e19, -1e18, 1e18, 1e19, 1e20]

SOMMERFELD_L = 2.44e-8   # WΩ/K²
# patch_kl_auto：晶格热导率来源。
#   "auto" = 自动读 kl-dft-cpu 技能的 step6_kappa/kappa_summary.json，取 TARGET_T 处的
#            **原始** kappa_xx_yy_zz（元胞口径）。本步 amset 的 σ/κ_e 也是元胞
#            口径，两边天然同口径，所以这里不做 c/t 重标度。
#            （要文献口径的绝对值请看 step8.1 的 NORM_2D。）
#   None   = 只用下面手填的值。
KAPPA_L_SOURCE = "auto"
# 多链探测：kl-dft-cpu → kl-mlff-cpu → kl-mlff-gpu（sibling 子目录，取第一个有结果的）。
KL_CHAIN_DIRS  = {
    "kl-dft-cpu":  ("step6_kappa", "step6_kappa_shengbte"),
    "kl-mlff-cpu": ("step4_kappa",),
    "kl-mlff-gpu": ("step4_kappa",),
}
KL_CHAIN_ORDER = ("kl-dft-cpu", "kl-mlff-cpu", "kl-mlff-gpu")
KL_ROOT        = None   # None=按 KL_CHAIN_ORDER 探测 sibling；填路径=只找该目录（配 KL_KAPPA_DIRS）
KL_KAPPA_DIRS  = ("step6_kappa", "step6_kappa_shengbte")   # 仅 KL_ROOT 手填时用
# 手填值（W/mK，元胞口径）。非 None 时优先于 auto。
KAPPA_L_XX_W_MK = None
KAPPA_L_YY_W_MK = None
E_C = 1.602176634e-19    # C
M0  = 9.1093837015e-31   # kg


# ---------- 通用 ----------
def _load_json(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return None


def _num(v):
    """数值才返回 float，其它（None/nan/字符串）返回 None。"""
    return float(v) if isinstance(v, (int, float)) and v == v else None


def _unwrap(v):
    """monty 的 {@class:array,data:[...]} → 原始 list；否则原样返回。"""
    if isinstance(v, dict) and v.get("@class") == "array":
        return v.get("data")
    return v


# ---------- 维度判定 & 张量约化 ----------
_STEP1_DIM_CANDS = ("step1_opt", "step1_std_opt",
                    "step1c_PBE_opt", "step1b_PBE_opt", "step1a_PBE_opt")
ANISO_TOL = 0.05   # 面内 xx/yy 相对差 > 5% 视为各向异性 → 分方向报值


def _read_dim(cwd):
    """从 step1 的 workflow_method.txt 读 DIM=（2d/3d/0d）；读不到返回 None。
    与 step8.1_boltztrap 的判维口径一致。"""
    for name in _STEP1_DIM_CANDS:
        mf = Path(cwd) / name / "workflow_method.txt"
        if not mf.is_file():
            continue
        for ln in mf.read_text(errors="ignore").splitlines():
            if ln.strip().upper().startswith("DIM="):
                return ln.split("=", 1)[1].strip().lower()
    return None


def overlap_grade_only_real(dim):
    """只有"真实重叠"运行、没有 unity 对照时的判级：返回 "2d" 或 "3d"。

    2d：[0021] 不再判红。V23（去对称化把重叠算坏）由 S8.4 的 overlap_preflight 按反演判据把关
        （全网格 h5 / IBZ + 相位补丁，不安全就拦截），2026-09-26 起真实重叠是 2D 出厂默认。
    3d / 维度读不到（本模块其它地方一律"按 3D 处理"）：黄 —— 受影响面已确认，
        是否真算错待 Si 全网格对照（V25）。
    """
    return "2d" if str(dim).lower() == "2d" else "3d"


def _diag_avg(t33, is_2d=False):
    """张量各向同性化。
    3D: (xx+yy+zz)/3。
    2D: (xx+yy)/2 —— 剔除真空方向 zz。zz 沿真空群速度 v_z≈0：对 σ/κ 而言
        分量≈0，是把正确值稀释成 2/3 的项；对 S=α/σ 而言是 0/0 的病态比值，
        污染方向和幅度都不可控，必须剔除。与 step8.1_boltztrap 保持一致。
    """
    if is_2d:
        return (t33[0][0] + t33[1][1]) / 2.0
    return (t33[0][0] + t33[1][1] + t33[2][2]) / 3.0


def _inplane_anisotropy(out):
    """S/σ/κ_e 面内 xx-yy 的最大相对差；判定是否需要分方向报值。"""
    def _reldiff(xx, yy):
        m = 0.0
        for a, b in zip(xx, yy):
            if a == a and b == b:            # 跳过 NaN
                d = max(abs(a), abs(b))
                if d > 0:
                    m = max(m, abs(a - b) / d)
        return m
    r = {"seebeck": _reldiff(out["seebeck_xx"], out["seebeck_yy"]),
         "sigma":   _reldiff(out["sigma_xx"],   out["sigma_yy"]),
         "kappa_e": _reldiff(out["kappa_e_xx"], out["kappa_e_yy"])}
    r["max"] = max(r["seebeck"], r["sigma"], r["kappa_e"])
    r["anisotropic"] = r["max"] > ANISO_TOL
    return r


def _nearest_T_index(temps, target):
    return min(range(len(temps)), key=lambda i: abs(temps[i] - target))


def _interp_loglog(x_signed, y, target_signed):
    """在同号子集里，对 |x| 做 log 插值取 target。x_signed/target_signed 同号才有效。"""
    import numpy as np
    if target_signed < 0:
        idx = [i for i in range(len(x_signed)) if x_signed[i] < 0]
    else:
        idx = [i for i in range(len(x_signed)) if x_signed[i] > 0]
    if len(idx) < 2:
        return None
    xs = np.array([math.log10(abs(x_signed[i])) for i in idx])
    ys = np.array([y[i] for i in idx])
    order = np.argsort(xs)
    xs, ys = xs[order], ys[order]
    xt = math.log10(abs(target_signed))
    if xt < xs.min() or xt > xs.max():
        return None      # 不外推
    # patch_logy：σ/κ_e/μ 在相邻掺杂点间能差一个数量级以上，y 线性插值会错得很远。
    # 正定且跨度 >10 倍的量改在 log10(y) 上插。S 可正可负、Lorenz 跨度小，走原路。
    if np.all(ys > 0) and (ys.max() / ys.min()) > 10.0:
        return float(10.0 ** np.interp(xt, xs, np.log10(ys)))
    return float(np.interp(xt, xs, ys))


# ---------- 读三方 ----------
def load_amset(cwd, is_2d=False, dirname=None):
    j = _load_json(Path(cwd) / (dirname or AMSET_DIR) / "transport.json")
    if not j:
        return None
    try:
        doping = _unwrap(j["doping"])
        temps = _unwrap(j["temperatures"])
        iT = _nearest_T_index(temps, TARGET_T)
        cond = _unwrap(j["conductivity"])
        seeb = _unwrap(j["seebeck"])
        kel = _unwrap(j["electronic_thermal_conductivity"])
        mob = j.get("mobility")
        mob_over = mob_adp = None
        if isinstance(mob, dict):
            mob_over = _unwrap(mob.get("overall") or mob.get("Overall"))
            mob_adp = _unwrap(mob.get("ADP") or mob.get("ACD"))  # 纯声学形变势
        out = {"T_used": temps[iT], "is_2d": bool(is_2d),
               "doping_signed": [float(x) for x in doping],
               "sigma": [], "seebeck": [], "kappa_e": [], "lorenz": [],
               "mobility": [], "mobility_adp": [],
               # 面内原始对角分量：各向异性诊断 + 分方向报值（3D 时也存着，不用而已）
               "sigma_xx": [], "sigma_yy": [],
               "seebeck_xx": [], "seebeck_yy": [],
               "kappa_e_xx": [], "kappa_e_yy": []}
        for n in range(len(doping)):
            s = _diag_avg(cond[n][iT], is_2d); se = _diag_avg(seeb[n][iT], is_2d)
            k = _diag_avg(kel[n][iT], is_2d)
            out["sigma"].append(s)
            out["seebeck"].append(se)      # μV/K（amset 单位）
            out["kappa_e"].append(k)
            out["lorenz"].append(k / (s * temps[iT]) if s > 0 else float("nan"))
            out["mobility"].append(_diag_avg(mob_over[n][iT], is_2d) if mob_over else float("nan"))
            out["mobility_adp"].append(_diag_avg(mob_adp[n][iT], is_2d) if mob_adp else float("nan"))
            out["sigma_xx"].append(cond[n][iT][0][0]);  out["sigma_yy"].append(cond[n][iT][1][1])
            out["seebeck_xx"].append(seeb[n][iT][0][0]); out["seebeck_yy"].append(seeb[n][iT][1][1])
            out["kappa_e_xx"].append(kel[n][iT][0][0]);  out["kappa_e_yy"].append(kel[n][iT][1][1])
        out["aniso"] = _inplane_anisotropy(out) if is_2d else None
        return out
    except Exception as e:
        print("[WARN] amset transport.json 解析失败：%s: %s" % (type(e).__name__, e))
        return None


def load_bt2(cwd):
    j = _load_json(Path(cwd) / BT2_DIR / "boltztrap_crta.json")
    if not j:
        return None
    try:
        temps = j["temperatures_K"]
        iT = _nearest_T_index(temps, TARGET_T)
        carr = j.get("carrier_conc_cm-3")
        if carr is None:
            return None
        # step8.1 修正后 carrier_conc 已是 amset 约定（负=n型电子、正=p型空穴），不再取负
        n_signed = list(carr[iT])
        res = {"T_used": temps[iT], "n_signed": n_signed,
               "seebeck": [v * 1e6 for v in j["seebeck_V_per_K"][iT]],  # V/K→µV/K
               "lorenz": j["lorenz_WOhm_per_K2"][iT],
               "sigma_over_tau": j["sigma_over_tau_S_per_m_s"][iT],
               "kappa_e_over_tau": j["kappa_e_over_tau_W_per_m_K_s"][iT]}

        def _dir(key, scale=1.0):        # patch_bt2_dir：8.1 打过 aniso 补丁才有
            a = j.get(key)
            try:
                return [v * scale for v in a[iT]] if a else None
            except (IndexError, TypeError):
                return None
        res["sigma_over_tau_xx"] = _dir("sigma_over_tau_xx")
        res["sigma_over_tau_yy"] = _dir("sigma_over_tau_yy")
        res["seebeck_xx"] = _dir("seebeck_xx_V_per_K", 1e6)
        res["seebeck_yy"] = _dir("seebeck_yy_V_per_K", 1e6)
        res["kappa_e_over_tau_xx"] = _dir("kappa_e_over_tau_xx")
        res["kappa_e_over_tau_yy"] = _dir("kappa_e_over_tau_yy")
        res["sigma_over_tau_zz"] = _dir("sigma_over_tau_zz")         # [0026] 3D 才有
        res["seebeck_zz"] = _dir("seebeck_zz_V_per_K", 1e6)
        res["kappa_e_over_tau_zz"] = _dir("kappa_e_over_tau_zz")
        return res
    except Exception as e:
        print("[WARN] boltztrap_crta.json 解析失败：%s: %s" % (type(e).__name__, e))
        return None


def _dpt_inplane_mu(bd, r):
    """by_direction 的面内平均 μ (cm²/Vs)；缺/非 2D 回退 header mobility_cm2_Vs。

    by_direction 用 full-BZ 二次型拟合的 m_d（m_d=√(m_x·m_y)），header 的
    mobility_cm2_Vs 是 3 点抛物 m* 的旧口径，两者 μ 能差近一倍（CrS2 45.6 vs 87.5）。
    主口径一律取 by_direction。3D 时 by_direction.status 非 "ok"，自然回退。"""
    if isinstance(bd, dict) and bd.get("status") == "ok":
        vals = [bd[d]["mobility_cm2_Vs"] for d in ("x", "y", "z")          # [0025] 3D 有 z
                if isinstance(bd.get(d), dict)
                and isinstance(bd[d].get("mobility_cm2_Vs"), (int, float))]
        if vals:
            return sum(vals) / len(vals)
    return r.get("mobility_cm2_Vs")


def _dpt_m_d(bd, r):
    """by_direction 的态密度质量 m_d（full-BZ 二次型）；缺/非 2D 回退 header m_eff。"""
    if isinstance(bd, dict) and isinstance(bd.get("m_d_m0"), (int, float)):
        return bd["m_d_m0"]
    return r.get("inputs", {}).get("m_eff_m0")


def _dpt_T_factor(j, T):
    """[0021] dpt_result.json 的 μ/τ 是它 temperature_K 下的值，换到 T：2D ×T0/T，3D ×(T0/T)^1.5。"""
    T0 = float(j.get("temperature_K") or 300.0)
    is2d = j.get("is_2d") if j.get("is_2d") is not None else (str(j.get("dim")).lower() == "2d")
    return ((T0 / T) if is2d else (T0 / T) ** 1.5), T0


def _scale_dir(bd, f):
    """by_direction 的 x/y μ 与 τ 同乘 f（返回新 dict，不改原数据）。"""
    if not (isinstance(bd, dict) and bd.get("status") == "ok") or f == 1.0:
        return bd
    out = dict(bd)
    for d in ("x", "y", "z"):                                   # [0025] 3D 有 z
        if isinstance(bd.get(d), dict):
            q = dict(bd[d])
            for k in ("mobility_cm2_Vs", "tau_s"):
                if isinstance(q.get(k), (int, float)):
                    q[k] = q[k] * f
            out[d] = q
    return out


def _mv_at(mv, T):
    """[0022] dpt_result.json 多谷块（multi_valley）在温度 T 的值 {x, y, mean, tau_x, tau_y}。多谷不是 ∝1/T
    （谷间布居随 T 变）：正好有这一档就取，夹在两档之间按 log–log 线性插值，网格外 -> None（不外推）。"""
    import math
    if not (isinstance(mv, dict) and mv.get("status") == "ok"):
        return None
    rows = sorted(mv.get("mobility_vs_T") or [], key=lambda r: r["T_K"])
    keys = (("x", "mobility_x_cm2_Vs"), ("y", "mobility_y_cm2_Vs"), ("tau_x", "tau_x_s"), ("tau_y", "tau_y_s"))
    hit = None
    for r in rows:
        if abs(float(r["T_K"]) - float(T)) < 1e-6:
            hit = {k: float(r[c]) for k, c in keys}
    for a, b in zip(rows, rows[1:]):
        if hit is None and float(a["T_K"]) < float(T) < float(b["T_K"]):
            t = math.log(float(T) / a["T_K"]) / math.log(b["T_K"] / a["T_K"])
            hit = {k: math.exp(math.log(a[c]) + t * math.log(b[c] / a[c])) for k, c in keys}
    if hit:
        hit["mean"] = (hit["x"] + hit["y"]) / 2.0
    return hit


def _scope_share(scope, T):
    """[0024] single_valley_scope 在温度 T（S8.2 网格上最近的一档）的 (档温度, 带边谷占比)；没有 -> (None, None)。"""
    rows = (scope or {}).get("shares_vs_T") or []
    if not rows:
        return None, None
    r = min(rows, key=lambda q: abs(float(q["T_K"]) - float(T)))
    return float(r["T_K"]), float(r["share"])


def _dpt_scope_lines(dpt):
    """[0024] 汇总里的"单谷 DPT 只代表哪个谷"。"""
    out = []
    for c in ("electron", "hole"):
        sc = (dpt or {}).get("scope_" + c)
        if not sc:
            continue
        T_, sh = _scope_share(sc, (dpt or {}).get("T_used") or TARGET_T)
        others = ", ".join("%s" % k for k in sc.get("others") or {})
        out.append("# [DPT 单谷] %s%s：dpt_mu 只代表 %s 谷（%.0f K 占载流子 %.1f%%；%s 未计入，谷间散射也不含）"
                   "—— 材料的迁移率看 AMSET 栏" % ("★ " if sh is not None and sh < 0.5 else "", c,
                                               sc["edge_star"], T_ or 0.0, 100 * (sh or 0.0), others or "其余谷"))
    return out


def _amset_mech_lines(cwd, is_2d, T=None):
    """[0024] AMSET 各散射机制单独的迁移率（每种载流子取 |掺杂| 最低一档、最接近 T 的温度）与限制因素。
    Matthiessen：1/μ ≈ Σ 1/μ_i，μ_i 最小的机制限制迁移率；"本征"不计 IMP（随掺杂变）。"""
    T = TARGET_T if T is None else T
    for run in (AMSET2D_DIR, AMSET_DIR):
        j = _load_json(Path(cwd) / run / "transport.json")
        if j:
            break
    else:
        return []
    try:
        dop = [float(x) for x in _unwrap(j["doping"])]
        temps = [float(t) for t in _unwrap(j["temperatures"])]
        mob = j.get("mobility") or {}
    except (KeyError, TypeError, ValueError):
        return []
    it = min(range(len(temps)), key=lambda i: abs(temps[i] - T))
    mechs = [m for m in mob if m.lower() != "overall" and mob[m] is not None]
    if not mechs:
        return []
    out = []
    for car, sgn in (("electron", -1), ("hole", 1)):
        idx = [i for i, x in enumerate(dop) if x * sgn > 0]
        if not idx:
            continue
        i = min(idx, key=lambda k: abs(dop[k]))
        vals = {}
        for m in mechs:
            try:
                v = _diag_avg(_unwrap(mob[m])[i][it], is_2d)
            except (IndexError, TypeError, KeyError):
                continue
            if isinstance(v, (int, float)) and v == v and v > 0:
                vals[m] = float(v)
        if not vals:
            continue
        order = sorted(vals.items(), key=lambda kv: kv[1])
        intr = [kv for kv in order if kv[0].upper() != "IMP"]
        out.append("# [AMSET 机制] %s %.0f K（掺杂 %.2e cm⁻³，%s）：%s → 限制因素 %s%s" % (
            car, temps[it], dop[i], run, " < ".join("%s %.4g" % kv for kv in order), order[0][0],
            ("；本征（不计 IMP）限制因素 %s" % intr[0][0]) if intr and intr[0][0] != order[0][0] else ""))
    return out


def _mv_desc(mv):
    """'Q×6 ΔE=0.000 / K×2 ΔE=0.030' —— 多谷块里有哪些谷。"""
    return " / ".join("%s×%d ΔE=%.3f" % (v["label"], v["g"], v["dE_eV"]) for v in (mv or {}).get("valleys") or [])


def load_dpt(cwd):
    j = _load_json(Path(cwd) / DPT_DIR / "dpt_result.json")
    if not j:
        return None
    # [0021] 对比在 TARGET_T：S8.2 的 TEMPERATURE_K 不是它时，μ 与 τ 一并换温（以前直接拿来比）
    f, T0 = _dpt_T_factor(j, TARGET_T)
    out = {"T_used": TARGET_T if abs(f - 1.0) > 1e-12 else T0, "T0": T0}
    if abs(f - 1.0) > 1e-12:
        print("[..] DPT 结果在 %.0f K，换算到 %.0f K（×%.4f）再对比" % (T0, TARGET_T, f))
    for r in j.get("results", []):
        carrier = r["carrier"]
        bd = r.get("by_direction")
        # [fix m*口径] 主 μ/m* 取 by_direction 的 m_d + 面内平均 μ，
        #   不再用 inputs 的 3 点抛物 m*（后者 μ 系统性偏高一倍）。
        mu = _dpt_inplane_mu(bd, r)
        out[carrier] = (mu * f) if isinstance(mu, (int, float)) else mu
        out["m_" + carrier] = _dpt_m_d(bd, r)
        out["dir_" + carrier] = _scale_dir(bd, f)   # patch_bt2_dir
        if r.get("single_valley_scope"):            # [0024] 单谷值只代表带边那组谷
            out["scope_" + carrier] = r["single_valley_scope"]
        # [0022] 多谷 DPT：另列一栏对照；单谷算不出（带边在 Q 这类非高对称谷）时主栏也用它
        mv = r.get("multi_valley") or {}
        mvT = _mv_at(mv, out["T_used"])
        if mvT:
            out["mv_" + carrier] = mvT
            out["mv_desc_" + carrier] = _mv_desc(mv)
            if not isinstance(out[carrier], (int, float)):
                out[carrier] = mvT["mean"]
                out["m_" + carrier] = (mv["x"]["m_c_m0"] + mv["y"]["m_c_m0"]) / 2.0
                out["dir_" + carrier] = {"status": "ok", "source": "multi_valley",
                                         "x": {"mobility_cm2_Vs": mvT["x"], "tau_s": mvT["tau_x"]},
                                         "y": {"mobility_cm2_Vs": mvT["y"], "tau_s": mvT["tau_y"]}}
                out["src_" + carrier] = "multi_valley"
                out.pop("scope_" + carrier, None)       # [0024] 主栏已是多谷，单谷的适用范围不适用
        elif mv.get("status") == "ok":
            print("[WARN] %s 多谷 DPT 没有 %.0f K 这一档（S8.2 的 TEMPERATURES 不含它，不外推）"
                  % (carrier, out["T_used"]))
        # [C6/C7] 真空对齐的 edge_flip / vac_align / provenance
        out["E1_prov_" + r["carrier"]] = r.get("inputs", {}).get("E1_provenance", "")
    # step7b_deform_read/band_edges.json 的 edge_flip / vac_align / window_scan
    be_path = Path(cwd) / "step7b_deform_read" / "band_edges.json"
    if be_path.is_file():
        be = _load_json(be_path)
        if be:
            out["vac_align"] = be.get("vac_align")
            for c in ("electron", "hole"):
                if c in be:
                    out.setdefault("edge_flip_" + c, be[c].get("edge_flip"))
                    out.setdefault("vac_scan_" + c, be[c].get("vac_window_scan"))
    return out


def _dpt_tau_s(dpt, carrier):
    """由 DPT 反推弛豫时间 τ = μ·m*/e（μ cm²/Vs→m²/Vs, m*=m_eff·m0）。返回秒或 None。"""
    if not dpt:
        return None
    mu = dpt.get(carrier)
    meff = dpt.get("m_" + carrier)
    if not (isinstance(mu, (int, float)) and isinstance(meff, (int, float))
            and mu > 0 and meff > 0):
        return None
    return (mu * 1e-4) * (meff * M0) / E_C


# ---------- 表 ----------
# === patch_zt：ZT = S²σT/(κ_L+κ_e)，分方向优先 ===
def _zt(S_uV, sigma, kappa_e, kappa_L, T):
    if None in (S_uV, sigma, kappa_e, kappa_L):
        return None
    if any(v != v for v in (S_uV, sigma, kappa_e)):      # NaN
        return None
    denom = kappa_L + kappa_e
    if denom <= 0:
        return None
    return (S_uV * 1e-6) ** 2 * sigma * T / denom


# === patch_kl_auto：从 kl-dft-cpu 自动取 κ_L ===
_KL_CACHE = {}


def _interp_T(temps, vals, T):
    """按温度线性插值；超范围返回 None（不外推）。"""
    if not temps or len(temps) != len(vals):
        return None
    pairs = sorted(zip(temps, vals))
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


def _find_kl_kappa(cwd):
    """多链探测 kl 的 kappa_summary.json。返回 (src, root) 或 (None, root)。"""
    if KL_ROOT:
        root = Path(KL_ROOT)
        for d in KL_KAPPA_DIRS:
            p = root / d / "kappa_summary.json"
            if p.is_file():
                return p, root
        return None, root
    base = Path(cwd).parent
    for chain in KL_CHAIN_ORDER:
        root = base / chain
        for d in KL_CHAIN_DIRS.get(chain, ()):
            p = root / d / "kappa_summary.json"
            if p.is_file():
                return p, root
    return None, base


def resolve_kappa_L(cwd):
    """返回 (kxx, kyy, lines)。手填优先；否则读 kl 链的 kappa_summary.json。
    只取原始 kappa_xx_yy_zz（元胞口径），与 amset 的 σ/κ_e 同口径。"""
    if "v" in _KL_CACHE:
        return _KL_CACHE["v"]
    if KAPPA_L_XX_W_MK is not None or KAPPA_L_YY_W_MK is not None:
        r = (KAPPA_L_XX_W_MK, KAPPA_L_YY_W_MK,
             ["kappa_L 来源：手填（元胞口径）"])
        _KL_CACHE["v"] = r
        return r
    if not KAPPA_L_SOURCE or str(KAPPA_L_SOURCE).lower() != "auto":
        r = (None, None, ["kappa_L 未提供 —— 不出 ZT 列"])
        _KL_CACHE["v"] = r
        return r
    src, root = _find_kl_kappa(cwd)
    if src is None:
        r = (None, None, ["kappa_L 未找到：sibling kl-dft-cpu/kl-mlff-cpu/kl-mlff-gpu 下都没有 kappa_summary.json",
                          "  kl 链跑完了吗？或用 KL_ROOT 显式指定 kl 的材料目录"])
        _KL_CACHE["v"] = r
        return r
    j = _load_json(src) or {}
    raw, temps = j.get("kappa_xx_yy_zz"), j.get("temperatures")
    if not raw or not temps:
        r = (None, None, ["kappa_L 读取失败：%s 缺 kappa_xx_yy_zz/temperatures"
                          % src.name])
        _KL_CACHE["v"] = r
        return r
    kxx = _interp_T(temps, [x[0] for x in raw], TARGET_T)
    kyy = _interp_T(temps, [x[1] for x in raw], TARGET_T)
    if _read_dim(cwd) == "3d" and all(len(x) > 2 for x in raw):            # [0026] 3D 才用 zz
        _KL_CACHE["zz"] = _interp_T(temps, [x[2] for x in raw], TARGET_T)
    lines = ["kappa_L 来源：%s（原始元胞口径，与 amset 同口径，不重标度）"
             % os.path.relpath(str(src), str(root))]
    if kxx is None or kyy is None:
        lines.append("  T=%.0f K 超出 kl-dft-cpu 的温度网格 [%.0f, %.0f] —— 不外推，不出 ZT"
                     % (TARGET_T, min(temps), max(temps)))
        r = (None, None, lines)
        _KL_CACHE["v"] = r
        return r
    lines.append("  kappa_L(%.0f K)  xx=%.4f  yy=%.4f W/mK" % (TARGET_T, kxx, kyy))
    if _KL_CACHE.get("zz") is not None:
        lines[-1] += "  zz=%.4f" % _KL_CACHE["zz"]
    # 元胞闸门：kl-dft-cpu 的 Lz 与 ke-dft-cpu 的 c
    # [0016] 2D 默认不跑 S8：胞高改读 S8.4 的记录
    rec = (_load_json(Path(cwd) / AMSET_DIR / "2d_correction.json")
           or _load_json(Path(cwd) / AMSET2D_DIR / "2d_correction.json") or {})
    Lz, cA = j.get("Lz_ang"), rec.get("cell_c_A")
    if Lz and cA:
        dc = abs(float(Lz) - float(cA)) / float(cA) * 100
        tag = "一致" if dc <= 1.0 else "★不一致★"
        lines.append("  [闸门] 元胞%s：kl-dft-cpu Lz=%.4f Å vs ke-dft-cpu c=%.4f Å（差 %.3f%%）"
                     % (tag, Lz, cA, dc))
        if dc > 1.0:
            lines.append("         两边 step1 都冻结 c，同一份输入 POSCAR 不该出现"
                         "这种差异——查 step.conf 的 CELL_POLICY/VACUUM_AXIS_POLICY")
    r = (kxx, kyy, lines)
    _KL_CACHE["v"] = r
    return r


def _add_zt_columns(row, has_am2=False):
    """按 KAPPA_L_* 补 ZT 列。amset 分方向 + 面内平均；bt2 走文献口径
    （σ=(σ/τ)·τ_DPT，κ_e=LσT）；has_am2=True 时同时补 amset2d 的 ZT。"""
    kxx, kyy, _ = resolve_kappa_L(Path.cwd())        # patch_kl_auto
    if kxx is None and kyy is None:
        return
    for d, kl in (("xx", kxx), ("yy", kyy)):
        if kl is None:
            continue
        row["amset_ZT_%s" % d] = _zt(row.get("amset_S_%s_uV/K" % d),
                                     row.get("amset_sigma_%s_S/m" % d),
                                     row.get("amset_kappa_e_%s_W/mK" % d),
                                     kl, TARGET_T)
        row["bt2_ZT_%s" % d] = _zt(row.get("bt2_S_%s_uV/K" % d),      # patch_bt2_dir
                                   row.get("bt2_sigma_%s_S/m" % d),
                                   row.get("bt2_kappa_e_WF_%s_W/mK" % d),
                                   kl, TARGET_T)
        if has_am2:                                                   # patch_amset2d
            row["amset2d_ZT_%s" % d] = _zt(row.get("amset2d_S_%s_uV/K" % d),
                                           row.get("amset2d_sigma_%s_S/m" % d),
                                           row.get("amset2d_kappa_e_%s_W/mK" % d),
                                           kl, TARGET_T)
    kzz = _KL_CACHE.get("zz")                                         # [0026] 只有 3D + auto 来源才有
    if kzz is not None and row.get("bt2_sigma_zz_S/m") is not None:
        row["bt2_ZT_zz"] = _zt(row.get("bt2_S_zz_uV/K"), row.get("bt2_sigma_zz_S/m"),
                               row.get("bt2_kappa_e_WF_zz_W/mK"), kzz, TARGET_T)
    kl_avg = ([k for k in (kxx, kyy, kzz) if k is not None])        # 3D：S、σ 是三方向平均，κL 也取三方向
    kl_avg = sum(kl_avg) / len(kl_avg)
    row["amset_ZT"] = _zt(row.get("amset_S_uV/K"), row.get("amset_sigma_S/m"),
                          row.get("amset_kappa_e_W/mK"), kl_avg, TARGET_T)
    row["bt2_ZT"] = _zt(row.get("bt2_S_uV/K"), row.get("bt2_sigma_S/m"),
                        row.get("bt2_kappa_e_WF_W/mK"), kl_avg, TARGET_T)
    if has_am2:                                                       # patch_amset2d
        row["amset2d_ZT"] = _zt(row.get("amset2d_S_uV/K"), row.get("amset2d_sigma_S/m"),
                                row.get("amset2d_kappa_e_W/mK"), kl_avg, TARGET_T)


def _dpt_tau_dir(dpt, carrier, d):
    """分方向 τ（秒）。step8.2 没打各向异性补丁时回落到各向同性 τ。"""
    bd = (dpt or {}).get("dir_" + carrier) or {}
    if bd.get("status") == "ok" and isinstance(bd.get(d), dict):
        return bd[d].get("tau_s")
    return _dpt_tau_s(dpt, carrier)


def _amset_cols(row, tag, a, nt):
    """把一路 amset 结果（原版 / amset2d）插值进 row，列名前缀 = tag。"""
    row["%s_sigma_S/m" % tag] = _interp_loglog(a["doping_signed"], a["sigma"], nt)
    row["%s_S_uV/K" % tag] = _interp_loglog(a["doping_signed"], a["seebeck"], nt)
    row["%s_kappa_e_W/mK" % tag] = _interp_loglog(a["doping_signed"], a["kappa_e"], nt)
    row["%s_Lorenz" % tag] = _interp_loglog(a["doping_signed"], a["lorenz"], nt)
    row["%s_mu_cm2/Vs" % tag] = _interp_loglog(a["doping_signed"], a["mobility"], nt)
    if any(v == v for v in a.get("mobility_adp", [])):
        row["%s_mu_ADP_cm2/Vs" % tag] = _interp_loglog(
            a["doping_signed"], a["mobility_adp"], nt)
    if a.get("aniso") and a["aniso"]["anisotropic"]:
        for d in ("xx", "yy"):
            row["%s_S_%s_uV/K" % (tag, d)] = _interp_loglog(
                a["doping_signed"], a["seebeck_%s" % d], nt)
            row["%s_sigma_%s_S/m" % (tag, d)] = _interp_loglog(
                a["doping_signed"], a["sigma_%s" % d], nt)
            row["%s_kappa_e_%s_W/mK" % (tag, d)] = _interp_loglog(
                a["doping_signed"], a["kappa_e_%s" % d], nt)


def build_table(am, bt, dpt, am2=None):
    rows = []
    for nt in TARGET_N:
        typ = "n(电子)" if nt < 0 else "p(空穴)"
        row = {"carrier_conc_cm-3": nt, "type": typ}
        if am:
            row["amset_sigma_S/m"] = _interp_loglog(am["doping_signed"], am["sigma"], nt)
            row["amset_S_uV/K"] = _interp_loglog(am["doping_signed"], am["seebeck"], nt)
            row["amset_kappa_e_W/mK"] = _interp_loglog(am["doping_signed"], am["kappa_e"], nt)
            row["amset_Lorenz"] = _interp_loglog(am["doping_signed"], am["lorenz"], nt)
            row["amset_mu_cm2/Vs"] = _interp_loglog(am["doping_signed"], am["mobility"], nt)
            if any(v == v for v in am.get("mobility_adp", [])):
                row["amset_mu_ADP_cm2/Vs"] = _interp_loglog(
                    am["doping_signed"], am["mobility_adp"], nt)
            # 2D 面内各向异性：附上分方向值（面内平均会抹掉它们）
            if am.get("aniso") and am["aniso"]["anisotropic"]:
                row["amset_S_xx_uV/K"] = _interp_loglog(am["doping_signed"], am["seebeck_xx"], nt)
                row["amset_S_yy_uV/K"] = _interp_loglog(am["doping_signed"], am["seebeck_yy"], nt)
                row["amset_sigma_xx_S/m"] = _interp_loglog(am["doping_signed"], am["sigma_xx"], nt)
                row["amset_sigma_yy_S/m"] = _interp_loglog(am["doping_signed"], am["sigma_yy"], nt)
                row["amset_kappa_e_xx_W/mK"] = _interp_loglog(am["doping_signed"], am["kappa_e_xx"], nt)
                row["amset_kappa_e_yy_W/mK"] = _interp_loglog(am["doping_signed"], am["kappa_e_yy"], nt)
        if am2:                                   # patch_amset2d：二维散射核
            _amset_cols(row, "amset2d", am2, nt)
        if bt:
            row["bt2_S_uV/K"] = _interp_loglog(bt["n_signed"], bt["seebeck"], nt)
            row["bt2_Lorenz"] = _interp_loglog(bt["n_signed"], bt["lorenz"], nt)
            row["bt2_sigma/tau"] = _interp_loglog(bt["n_signed"], bt["sigma_over_tau"], nt)
        if dpt:
            row["dpt_mu_cm2/Vs"] = dpt.get("electron") if nt < 0 else dpt.get("hole")
            _mv = dpt.get("mv_electron" if nt < 0 else "mv_hole")      # [0022]
            if _mv:
                row["dpt_mv_mu_cm2/Vs"] = _mv["mean"]
        # κ_e 三方（BT2 用 DPT-τ 反推、DPT 用 WF 估）
        carrier = "electron" if nt < 0 else "hole"
        if bt and dpt:
            tau = _dpt_tau_s(dpt, carrier)
            ket = _interp_loglog(bt["n_signed"], bt["kappa_e_over_tau"], nt)
            row["bt2_kappa_e_W/mK"] = (ket * tau) if (ket is not None and tau) else None
            # patch_bt2_abs：CRTA+DPT（= 文献做法）的绝对 σ、PF 与 WF-κ_e。
            #   注意 bt2_kappa_e_W/mK 是零电流下的 κ_e（BoltzTraP2 已扣 S²σT），
            #   bt2_kappa_e_WF_W/mK 才是文献 Eq.(5) 的 κ_e = LσT，两者不是同一个量。
            sot = _interp_loglog(bt["n_signed"], bt["sigma_over_tau"], nt)
            sig = (sot * tau) if (sot is not None and tau) else None
            row["bt2_sigma_S/m"] = sig
            _s = row.get("bt2_S_uV/K")
            if sig is not None and _s is not None and _s == _s:
                row["bt2_PF_W/mK2"] = (_s * 1e-6) ** 2 * sig
                row["bt2_kappa_e_WF_W/mK"] = SOMMERFELD_L * sig * TARGET_T
            # patch_bt2_dir：每个方向用它自己的 τ_α —— 这才是文献的做法
            for _d, _dd in (("xx", "x"), ("yy", "y"), ("zz", "z")):     # [0026] zz 只有 3D 的 S8.1 才写
                _sot = bt.get("sigma_over_tau_%s" % _d)
                _sbk = bt.get("seebeck_%s" % _d)
                if not _sot:
                    continue
                _td = _dpt_tau_dir(dpt, carrier, _dd)
                _sd = _interp_loglog(bt["n_signed"], _sot, nt)
                _sgd = (_sd * _td) if (_sd is not None and _td) else None
                row["bt2_sigma_%s_S/m" % _d] = _sgd
                _sv = _interp_loglog(bt["n_signed"], _sbk, nt) if _sbk else None
                row["bt2_S_%s_uV/K" % _d] = _sv
                if _sgd is not None and _sv is not None and _sv == _sv:
                    row["bt2_PF_%s_W/mK2" % _d] = (_sv * 1e-6) ** 2 * _sgd
                    row["bt2_kappa_e_WF_%s_W/mK" % _d] = SOMMERFELD_L * _sgd * TARGET_T
        if dpt:
            mu = dpt.get(carrier)
            if isinstance(mu, (int, float)) and mu > 0:
                sig = abs(nt) * 1e6 * E_C * (mu * 1e-4)      # n cm^-3→m^-3; S/m
                row["dpt_kappa_e_WF_W/mK"] = SOMMERFELD_L * sig * TARGET_T
        _add_zt_columns(row, has_am2=bool(am2))        # patch_zt / patch_amset2d
        rows.append(row)
    return rows


def write_table(out, rows, am, bt, dpt, am2=None):
    import csv
    cols = ["carrier_conc_cm-3", "type", "amset_sigma_S/m", "amset_S_uV/K",
            "amset_kappa_e_W/mK", "amset_Lorenz", "amset_mu_cm2/Vs",
            "amset_mu_ADP_cm2/Vs"]
    if am2:                                   # patch_amset2d：二维散射核一路
        cols += ["amset2d_sigma_S/m", "amset2d_S_uV/K", "amset2d_kappa_e_W/mK",
                 "amset2d_Lorenz", "amset2d_mu_cm2/Vs", "amset2d_mu_ADP_cm2/Vs",
                 "amset2d_S_xx_uV/K", "amset2d_S_yy_uV/K",
                 "amset2d_sigma_xx_S/m", "amset2d_sigma_yy_S/m",
                 "amset2d_kappa_e_xx_W/mK", "amset2d_kappa_e_yy_W/mK"]
    # 2D 面内各向异性时，把分方向列插在平均列后面
    if am and am.get("aniso") and am["aniso"]["anisotropic"]:
        cols += ["amset_S_xx_uV/K", "amset_S_yy_uV/K",
                 "amset_sigma_xx_S/m", "amset_sigma_yy_S/m",
                 "amset_kappa_e_xx_W/mK", "amset_kappa_e_yy_W/mK"]
    cols += ["bt2_S_uV/K", "bt2_Lorenz", "bt2_sigma/tau", "bt2_sigma_S/m",
             "bt2_PF_W/mK2", "bt2_kappa_e_W/mK", "bt2_kappa_e_WF_W/mK",
             "dpt_mu_cm2/Vs", "dpt_kappa_e_WF_W/mK"]
    cols += [c for c in ("dpt_mv_mu_cm2/Vs",) if any(c in r for r in rows)]       # [0022]
    cols += [c for c in ("bt2_S_xx_uV/K", "bt2_S_yy_uV/K", "bt2_S_zz_uV/K",            # patch_bt2_dir
                         "bt2_sigma_xx_S/m", "bt2_sigma_yy_S/m", "bt2_sigma_zz_S/m",
                         "bt2_PF_xx_W/mK2", "bt2_PF_yy_W/mK2", "bt2_PF_zz_W/mK2",
                         "bt2_kappa_e_WF_xx_W/mK", "bt2_kappa_e_WF_yy_W/mK", "bt2_kappa_e_WF_zz_W/mK")
             if any(c in r for r in rows)]
    if any(resolve_kappa_L(Path.cwd())[:2]):                          # patch_kl_auto
        cols += [c for c in ("amset_ZT_xx", "amset_ZT_yy", "amset_ZT",
                             "amset2d_ZT_xx", "amset2d_ZT_yy", "amset2d_ZT",
                             "bt2_ZT_xx", "bt2_ZT_yy", "bt2_ZT_zz", "bt2_ZT")
                 if any(c in r for r in rows)]
    # [patch_dose_range] 每原胞载流子数 + 可解读范围标记（2D/3D 同一判据）
    _c_cm = _areal_factor_cm(Path.cwd())
    _A_cm2, _V_cm3 = _cell_geometry(Path.cwd())
    _add_dose = bool(am or am2 or bt or dpt) and bool(_A_cm2 or _V_cm3)
    if _add_dose:
        _dcols = (["n_2D_cm-2"] if (_c_cm and (am or am2)) else [])
        cols = _dcols + ["carriers_per_cell", "in_interpret_range"] + cols
    # [patch_deform_geom] 形变势构型口径整列输出（每行同值；回退到 clamped 时一眼可见）
    _dgeom = _deform_geometry(Path.cwd())
    if _dgeom:
        cols = cols + ["deform_geometry"]
    _func, _func_src = _functional(Path.cwd())             # [V181]
    if _func:
        cols = cols + ["functional"]
    with open(out / "comparison_300K.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            _row = {k: ("" if r.get(k) is None else r.get(k)) for k in cols}
            if _add_dose and r.get("carrier_conc_cm-3") is not None:
                try:
                    _n3a = abs(float(r["carrier_conc_cm-3"]))
                    # [2026-09-16 修订] 2D/3D **同一个公式**：n_cell = n_3D × V_cell。
                    # 因为 V_cell = A_cell·h⊥、n_2D = n_3D·h⊥，所以 n_2D×A_cell 恒等于
                    # n_3D×V_cell —— 少一个分支就少一处出错的可能（原先 2D 走 A_cell、
                    # 3D 走 V_cell 是多余的）。n_2D 那一列只作展示。
                    _ncell = _n3a * _V_cm3 if _V_cm3 else None
                    if _c_cm and (am or am2):
                        _row["n_2D_cm-2"] = "%.3g" % (_n3a * _c_cm)
                    if _ncell is not None:
                        _row["carriers_per_cell"] = "%.3g" % _ncell
                        _row["in_interpret_range"] = (
                            "yes" if _ncell <= NCELL_INTERPRET_MAX else
                            "NO(>%.2g/cell：超出刚性能带近似/实验范围)" % NCELL_INTERPRET_MAX)
                except (TypeError, ValueError):
                    pass
            if _dgeom:
                _row["deform_geometry"] = _dgeom
            if _func:
                _row["functional"] = _func
            w.writerow(_row)

    def fmt(v, f="%.4g"):
        return (f % v) if isinstance(v, (int, float)) and v == v else "—"
    _red = "面内(xx+yy)/2 [2D]" if (am and am.get("is_2d")) else "对角(xx+yy+zz)/3 [3D]"
    lines = ["# 300 K 三方电子输运对比（详见 comparison_300K.csv）",
             "# skill_rev: %s" % _SKILL_REV,
             "# 有 amset=%s  amset2d=%s  BoltzTraP2=%s  DPT=%s"
             % (bool(am), bool(am2), bool(bt), bool(dpt)),
             "# 张量约化：%s" % _red,
             "# S/Lorenz 可直接比；σ/κ_e amset是绝对值、BT2是per-τ只比趋势；DPT迁移率仅ADP",
             "# [DPT μ] m* 取 full-BZ 二次型 m_d（面内平均），非 3 点抛物拟合"]
    for _c in ("electron", "hole"):                                   # [0022] 多谷 DPT
        if dpt and dpt.get("mv_" + _c):
            lines.append("# [DPT 多谷] %s：μ=%.1f cm²/Vs（%s；无谷间散射，上限）%s" % (
                _c, dpt["mv_" + _c]["mean"], dpt.get("mv_desc_" + _c, ""),
                "  ← 单谷算不出（带边在非高对称谷），dpt_mu 主栏与 BT2 的 τ 都用它"
                if dpt.get("src_" + _c) == "multi_valley" else ""))
    lines += _dpt_scope_lines(dpt)                                     # [0024]
    lines += _amset_mech_lines(Path.cwd(), _read_dim(Path.cwd()) == "2d")
    if not _func:
        lines.append("# [口径] 泛函 = 未知（step1 的 workflow_method.txt 读不到）—— 先查 step1 INCAR 的 GGA/IVDW，"
                     "再和文献比")
    else:
        _fnote = ("；GGA=PS 是 PBEsol 不是 PBE，晶格通常比 PBE 小约 1%，C_2D 偏硬，m* 和谷间能差也会变"
                  if _func.startswith("PBEsol") else "")
        lines.append("# [口径] 泛函 = %s（来源 %s）—— 和文献比较前先核对文献用的泛函%s"
                     % (_func, _func_src, _fnote))
    if _dgeom:
        _warn = ("" if _dgeom == "ionrelax" else
                 "  ⚠ 非离子弛豫口径：与 TOTAL ELASTIC MODULI 不同源，ADP 绝对值跨项目不可直接比")
        lines.append("# [口径] 形变势构型 deform_geometry = %s（ionrelax=离子弛豫；clamped=离子固定）%s"
                     % (_dgeom, _warn))
    if _add_dose:
        lines.append("# [范围] 可解读上限 = 每原胞 %.2g 个载流子"
                     "（n_cell = n_2D×A_cell 或 n_3D×V_cell；本胞 A = %s Å², V = %s Å³）%s。"
                     "超出档位在 csv 里标 in_interpret_range=NO：刚性能带近似已不成立、"
                     "实验也达不到，只看趋势，不进论文图。"
                     % (NCELL_INTERPRET_MAX,
                        ("%.3g" % (_A_cm2 * 1e16)) if _A_cm2 else "?",
                        ("%.3g" % (_V_cm3 * 1e24)) if _V_cm3 else "?",
                        ("；面密度换算因子 c = %.4f cm" % _c_cm) if _c_cm else ""))
    # [C6/C7] DPT 真空对齐结局 / edge_flip / 来源
    if dpt:
        va = dpt.get("vac_align") or {}
        if va.get("status") == "skipped":
            lines.append("# [DPT 真空对齐] 跳过：%s" % va.get("reason", "未知原因"))
        elif va.get("status") == "ok":
            lines.append("# [DPT 真空对齐] 成功 (window_half=%.2f, flat_tol=%.4g eV)" % (
                va.get("window_half", 0.25), va.get("flat_tol_eV", 1e-3)))
        for c in ("electron", "hole"):
            flip = dpt.get("edge_flip_" + c)
            if flip:
                lines.append("# [WARN] %s 带边占据态翻转（形变后带序交换）：%s" % (c, "; ".join(flip)))
            prov = dpt.get("E1_prov_" + c, "")
            if "真空对齐" in prov:
                lines.append("# [DPT] %s E1 来源: %s" % (c, prov))
    if am and am.get("aniso") and am["aniso"]["anisotropic"]:
        a = am["aniso"]
        lines.append("# [各向异性] 面内 xx-yy 最大相对差 %.1f%%"
                     "（S %.1f%% / σ %.1f%% / κ_e %.1f%%）："
                     "见 CSV 的 *_xx/*_yy 列，ZT 建议分方向算"
                     % (a["max"] * 100, a["seebeck"] * 100,
                        a["sigma"] * 100, a["kappa_e"] * 100))
    if am2:                                   # patch_amset2d：逐档给 amset2d vs amset 的比值
        lines.append("# [amset2d] 二维散射核（插件）vs 原版 amset 的三维核：")
        lines.append("#   弹性常数 = 原始 slab 值（标准 Voigt），散射核为二维形式；"
                     "差异来源 = 核的 q 依赖 + 原始 slab 弹性 + 面内 pop_frequency。")
        _m2 = [r for r in rows if _num(r.get("amset2d_mu_cm2/Vs"))]
        if _m2 and am:
            _rt = []
            for r in _m2:
                a1, a2 = _num(r.get("amset_mu_cm2/Vs")), _num(r.get("amset2d_mu_cm2/Vs"))
                if a1 and a2:
                    _rt.append(a2 / a1)
            if _rt:
                lines.append("#   μ(amset2d)/μ(amset)：%.2f ~ %.2f（中位 %.2f）"
                             % (min(_rt), max(_rt), sorted(_rt)[len(_rt) // 2]))
        if am2.get("aniso") and am2["aniso"]["anisotropic"]:
            _a = am2["aniso"]
            lines.append("#   [amset2d 各向异性] 面内 xx-yy 最大相对差 %.1f%%"
                         "（S %.1f%% / σ %.1f%% / κ_e %.1f%%）"
                         % (_a["max"] * 100, _a["seebeck"] * 100,
                            _a["sigma"] * 100, _a["kappa_e"] * 100))
    lines.append("")
    hdr = "%-13s %-8s | %-10s %-9s %-9s |" % (
        "n(cm^-3)", "type", "amsetS", "amsetL", "amsetμ")
    if am2:
        hdr += " %-10s %-9s |" % ("2dS", "2dμ")
    hdr += " %-9s %-9s | %-9s" % ("bt2 S", "bt2 L", "dptμ")
    lines.append(hdr); lines.append("-" * len(hdr))
    for r in rows:
        ln = "%-13.2e %-8s | %-10s %-9s %-9s |" % (
            r["carrier_conc_cm-3"], r["type"],
            fmt(r.get("amset_S_uV/K")), fmt(r.get("amset_Lorenz")),
            fmt(r.get("amset_mu_cm2/Vs")))
        if am2:
            ln += " %-10s %-9s |" % (fmt(r.get("amset2d_S_uV/K")),
                                     fmt(r.get("amset2d_mu_cm2/Vs")))
        ln += " %-9s %-9s | %-9s" % (fmt(r.get("bt2_S_uV/K")),
                                     fmt(r.get("bt2_Lorenz")),
                                     fmt(r.get("dpt_mu_cm2/Vs")))
        lines.append(ln)
    (out / "comparison_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------- 图 ----------
# ---------- [0020] 迁移率随温度：DPT / AMSET 只 ADP / AMSET 本征 / AMSET 全机制 ----------
MU_T_CSV = "comparison_vs_T.csv"
MU_T_JSON = "comparison_vs_T.json"
MU_T_PNG = "mobility_vs_T.png"
MU_T_COLS = ("mu_DPT", "mu_DPT_mv", "mu_AMSET_ADP", "mu_AMSET_intrinsic", "mu_AMSET_all")   # [0022] +mu_DPT_mv
# [0030] 分方向参数列（来自 S8.2 by_direction，逐方向）。2D=x/y，3D=x/y/z；不出比值/平均。
def _c_col(is_2d):
    """C 列名（单位随维度变）：2D 用 C_2D_N_per_m（N/m），3D 用 C_3D_GPa（GPa）。"""
    return "C_2D_N_per_m" if is_2d else "C_3D_GPa"


def _param_cols(is_2d):
    """[0030] 分方向参数列：E1 / C / m* / m_d / tau。"""
    return ("E1_eV", _c_col(is_2d), "m_eff_m0", "m_d_m0", "tau_DPT_fs")


def _dirs(is_2d):
    """分方向列表（**不出平均**）：2D 面内 x/y；3D x/y/z。"""
    return ("x", "y") if is_2d else ("x", "y", "z")


def _t_dirs(t33, is_2d):
    """3×3 张量 -> {"x","y",("z")}（按方向，不出平均）。"""
    d = {"x": float(t33[0][0]), "y": float(t33[1][1])}
    if not is_2d:
        d["z"] = float(t33[2][2])
    return d


def _amset_mu_T_sources(cwd):
    """选 AMSET 运行目录（2D 有 step8.4_amset2d 就用它）并读三套迁移率。
    返回 {"run", "doping", "temps", "cols": {列名: arr[ndop][nT][3][3]}, "src": {列名: 说明}, "intrinsic_labels"}；没有 -> None。"""
    for run in (AMSET2D_DIR, AMSET_DIR):
        j = _load_json(Path(cwd) / run / "transport.json")
        if j:
            break
    else:
        return None
    doping = [float(x) for x in _unwrap(j["doping"])]
    temps = [float(t) for t in _unwrap(j["temperatures"])]
    mob = j.get("mobility") or {}
    cols, src = {}, {}
    over = _unwrap(mob.get("overall") or mob.get("Overall"))
    if over is not None:
        cols["mu_AMSET_all"], src["mu_AMSET_all"] = over, "%s/transport.json mobility.overall" % run

    def _intr(name):
        q = _load_json(Path(cwd) / run / name)
        if not q:
            return None, None
        try:
            ok = (len(q["doping_cm3"]) == len(doping) and len(q["temperatures_K"]) == len(temps)
                  and all(abs(float(a) - b) <= 1e-6 * max(1.0, abs(b)) for a, b in zip(q["doping_cm3"], doping))
                  and all(abs(float(a) - b) < 1e-6 for a, b in zip(q["temperatures_K"], temps)))
        except (KeyError, TypeError, ValueError):
            ok = False
        if not ok:
            print("[WARN] %s/%s 的掺杂/温度网格与 transport.json 不同（不是同一次运行？）——不用" % (run, name))
            return None, None
        return q["intrinsic"]["mobility_cm2_Vs_s"].get("overall"), q

    adp, q_adp = _intr("intrinsic_ADP.json")
    if adp is not None:
        cols["mu_AMSET_ADP"] = adp
        src["mu_AMSET_ADP"] = "%s/intrinsic_ADP.json（只留 %s 重积分）" % (
            run, "+".join(q_adp.get("intrinsic_labels") or ["ADP"]))
    else:
        a = _unwrap(mob.get("ADP") or mob.get("ACD"))
        if a is not None:
            cols["mu_AMSET_ADP"] = a
            src["mu_AMSET_ADP"] = "%s/transport.json mobility.ADP（单机制）" % run
    intr, q_in = _intr("intrinsic_transport.json")
    labels = None
    if intr is not None:
        labels = list(q_in.get("intrinsic_labels") or [])
        cols["mu_AMSET_intrinsic"] = intr
        src["mu_AMSET_intrinsic"] = "%s/intrinsic_transport.json（%s，去 %s）" % (
            run, "+".join(labels) or "?", "+".join(q_in.get("intrinsic_dropped") or []) or "?")
    return {"run": run, "doping": doping, "temps": temps, "cols": cols, "src": src,
            "intrinsic_labels": labels}


def _dpt_mu_T(cwd, is_2d):
    """dpt_result.json -> (T0, {载流子: {"x","y",("z"),"mean"} @T0}, 来源)；没有 -> (None, {}, None)。"""
    j = _load_json(Path(cwd) / DPT_DIR / "dpt_result.json")
    if not j:
        return None, {}, None
    T0 = float(j.get("temperature_K") or 300.0)
    out, how = {}, {}
    for r in j.get("results", []):
        bd = r.get("by_direction")
        if is_2d and isinstance(bd, dict) and bd.get("status") == "ok":
            d = {"x": bd["x"]["mobility_cm2_Vs"], "y": bd["y"]["mobility_cm2_Vs"]}
            d["mean"] = (d["x"] + d["y"]) / 2.0
            how[r["carrier"]] = "by_direction"
        elif (not is_2d) and isinstance(bd, dict) and bd.get("status") == "ok" and isinstance(bd.get("z"), dict):
            d = {k: bd[k]["mobility_cm2_Vs"] for k in ("x", "y", "z")}        # [0025] 3D 分方向
            d["mean"] = (d["x"] + d["y"] + d["z"]) / 3.0
            how[r["carrier"]] = "by_direction(3D)"
        elif isinstance(r.get("mobility_cm2_Vs"), (int, float)):
            v = float(r["mobility_cm2_Vs"])
            d = {k: v for k in _dirs(is_2d)}
            how[r["carrier"]] = "header（各向同性）"
        else:
            continue
        out[r["carrier"]] = d
    return T0, out, how


def _dpt_params_dir(cwd, carrier, is_2d):
    """[0030] 从 dpt_result.json 的 by_direction 取**分方向**参数（E1/C2D/m*/m_d/tau）。
    返回 {direction: {...}}；没有 by_direction(ok) 返回 {}。2D=x/y，3D=x/y/z；不出平均、不出比值。"""
    j = _load_json(Path(cwd) / DPT_DIR / "dpt_result.json") or {}
    for r in j.get("results", []):
        if r.get("carrier") != carrier:
            continue
        bd = r.get("by_direction") or {}
        if bd.get("status") != "ok":
            return {}
        out, md = {}, bd.get("m_d_m0")
        for d in _dirs(is_2d):
            q = bd.get(d)
            if not isinstance(q, dict):
                continue
            t = q.get("tau_s")
            ck = _c_col(is_2d)
            out[d] = {"E1_eV": q.get("E1_eV"), ck: q.get(ck),
                      "m_eff_m0": q.get("m_eff_m0"), "m_d_m0": md,
                      "tau_DPT_fs": (t * 1e15) if t is not None else None}
        return out
    return {}


def _dpt_scopes(cwd):
    """[0024] {载流子: single_valley_scope}。"""
    j = _load_json(Path(cwd) / DPT_DIR / "dpt_result.json") or {}
    return {r["carrier"]: r["single_valley_scope"] for r in j.get("results", []) if r.get("single_valley_scope")}


def _dpt_mv_blocks(cwd, is_2d):
    """[0022] {载流子: multi_valley 块}（只 2D、status ok 的）。"""
    j = _load_json(Path(cwd) / DPT_DIR / "dpt_result.json") or {}
    if not is_2d:
        return {}
    return {r["carrier"]: r["multi_valley"] for r in j.get("results", [])
            if (r.get("multi_valley") or {}).get("status") == "ok"}


def build_mu_vs_T(cwd, is_2d):
    """[0020] μ(T) 三值表的行。温度网格 = AMSET 的；没有 AMSET 就用 dpt_result.json 的 temperatures_K。"""
    am = _amset_mu_T_sources(cwd)
    T0, dpt, dpt_how = _dpt_mu_T(cwd, is_2d)
    mvb = _dpt_mv_blocks(cwd, is_2d)
    scopes = _dpt_scopes(cwd)                                          # [0024]
    if not am and not dpt and not mvb:
        return None
    cfac = _areal_factor_cm(cwd) if is_2d else None
    if am:
        temps = am["temps"]
        picks = {}
        for i, x in enumerate(am["doping"]):
            car = "electron" if x < 0 else ("hole" if x > 0 else None)
            if car and (car not in picks or abs(x) < abs(am["doping"][picks[car]])):
                picks[car] = i
    else:
        jd = _load_json(Path(cwd) / DPT_DIR / "dpt_result.json") or {}
        temps = [float(t) for t in (jd.get("temperatures_K") or [T0])]
        picks = {c: None for c in list(dpt) + list(mvb)}
    rows = []
    for car in ("electron", "hole"):
        if car not in picks and car not in dpt and car not in mvb:
            continue
        i = picks.get(car)
        params = _dpt_params_dir(cwd, car, is_2d)              # [0030] 分方向 E1/C2D/m*/m_d/tau
        for it, T in enumerate(temps):
            f = ((T0 / T) if is_2d else (T0 / T) ** 1.5) if T0 else None
            vals = {}
            if car in dpt and f:
                vals["mu_DPT"] = {k: v * f for k, v in dpt[car].items()}
            mvT = _mv_at(mvb.get(car), T)
            if mvT:
                vals["mu_DPT_mv"] = {"x": mvT["x"], "y": mvT["y"], "mean": mvT["mean"]}
            if am and i is not None:
                for c, arr in am["cols"].items():
                    try:
                        vals[c] = _t_dirs(arr[i][it], is_2d)
                    except (IndexError, TypeError):
                        pass
            x = am["doping"][i] if (am and i is not None) else None
            for d in _dirs(is_2d):
                row = {"T_K": T, "carrier": car, "direction": d, "doping_cm3": x,
                       "n2D_cm2": (abs(x) * cfac if (x is not None and cfac) else None)}
                for c in MU_T_COLS:
                    v = (vals.get(c) or {}).get(d)
                    row[c] = v if (v is not None and v == v) else None
                for c in _param_cols(is_2d):                            # [0030] 分方向参数
                    row[c] = (params.get(d) or {}).get(c)
                _ts, _sh = _scope_share(scopes.get(car), T)             # [0024] mu_DPT 代表的谷占多少载流子
                row["DPT_edge_share"] = (_sh if (row["mu_DPT"] is not None and _ts is not None
                                                 and abs(_ts - T) < 1e-6) else None)
                rows.append(row)
    meta = {"is_2d": is_2d, "dpt_T0_K": T0, "dpt_source": dpt_how,
            "dpt_scaling": ("1/T" if is_2d else "T^-1.5"),
            "dpt_mv": ({c: _mv_desc(b) for c, b in mvb.items()} or None),
            "dpt_scope": ({c: "%s 谷" % sc["edge_star"] for c, sc in scopes.items()} or None),
            "amset_run": am["run"] if am else None, "amset_sources": am["src"] if am else {},
            "intrinsic_labels": am["intrinsic_labels"] if am else None,
            "doping_choice": "每种载流子取 |掺杂| 最低的一档（最接近 DPT 的非简并口径）",
            "temperatures_K": temps, "skill_rev": _SKILL_REV}
    return {"meta": meta, "rows": rows}


def write_mu_vs_T(out, cwd, is_2d):
    """[0020] 写 comparison_vs_T.csv / .json + mobility_vs_T.png，并在 comparison_summary.txt 末尾追加 μ(T) 表。"""
    import csv
    res = build_mu_vs_T(cwd, is_2d)
    if not res or not res["rows"]:
        print("[..] μ(T) 表：没有 AMSET 也没有 DPT 结果，跳过")
        return None
    head = (["T_K", "carrier", "direction", "doping_cm3", "n2D_cm2"]
            + list(_param_cols(is_2d)) + list(MU_T_COLS) + ["DPT_edge_share"])

    def _fmt(v):
        return "" if v is None else (("%.6g" % v) if isinstance(v, float) else v)
    with open(out / MU_T_CSV, "w", newline="") as fh:
        m = res["meta"]
        fh.write("# 迁移率随温度（cm²/Vs）；DPT 由 %s K 按 %s 换算（%s）；AMSET=%s；%s\n" % (
            m["dpt_T0_K"], m["dpt_scaling"], m["dpt_source"], m["amset_run"], m["doping_choice"]))
        if m.get("dpt_scope"):                                 # [0024]
            fh.write("# DPT_edge_share = mu_DPT（单谷）所代表的谷在该温度占载流子的比例：%s；其余谷与谷间散射不含\n"
                     % "；".join("%s %s" % kv for kv in sorted(m["dpt_scope"].items())))
        if m.get("dpt_mv"):                                    # [0022]
            fh.write("# mu_DPT_mv = S8.2 多谷 DPT（逐温度算，不是 1/T 换算；无谷间散射，上限）：%s\n"
                     % "；".join("%s %s" % kv for kv in sorted(m["dpt_mv"].items())))
        for c, sdesc in sorted(m["amset_sources"].items()):
            fh.write("# %s <- %s\n" % (c, sdesc))
        w = csv.writer(fh)
        w.writerow(head)
        for r in res["rows"]:
            w.writerow([_fmt(r[k]) for k in head])
    (out / MU_T_JSON).write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    lines = ["", "# [0030] 迁移率随温度 + 分方向参数（%s；E1 eV / C2D N/m / m* mₑ / tau fs / mu cm²/Vs；全表见 %s）"
             % ("2D：面内 x/y" if is_2d else "3D：x/y/z", MU_T_CSV),
             "# %-8s %4s %5s %8s %8s %7s %7s %8s %9s %9s %9s %9s %9s" % (
                 "carrier", "dir", "T(K)", "E1", "C2D", "m*", "m_d", "tau",
                 "DPT", "DPT_mv", "AMSET_ADP", "intrinsic", "all")]
    for r in res["rows"]:
        lines.append("  %-8s %4s %5.0f %8s %8s %7s %7s %8s %9s %9s %9s %9s %9s" % (
            r["carrier"], r["direction"], r["T_K"],
            *[("%.3f" % r[c]) if r.get(c) is not None else "—" for c in _param_cols(is_2d)],
            *[("%.1f" % r[c]) if r[c] is not None else "—" for c in MU_T_COLS]))
    try:
        with open(out / "comparison_summary.txt", "a", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
    except OSError:
        pass
    try:
        _plot_mu_vs_T(out, res)
    except Exception as e:                                     # noqa: BLE001
        print("[WARN] μ(T) 图失败（%s: %s）——表已生成" % (type(e).__name__, e))
    print("[OK] μ(T) 表：%s / %s（%d 个温度）" % (MU_T_CSV, MU_T_JSON, len(res["meta"]["temperatures_K"])))
    return res


def _plot_mu_vs_T(out, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style = {"mu_DPT": ("k", "-", "DPT"), "mu_DPT_mv": ("0.5", "-", "DPT multi-valley"),
             "mu_AMSET_ADP": ("C0", "--", "AMSET ADP"),
             "mu_AMSET_intrinsic": ("C2", "-.", "AMSET intrinsic"), "mu_AMSET_all": ("C3", ":", "AMSET all")}
    cars = [c for c in ("electron", "hole") if any(r["carrier"] == c for r in res["rows"])]
    fig, axs = plt.subplots(1, len(cars), figsize=(5.5 * len(cars), 4.2), squeeze=False)
    for ax, car in zip(axs[0], cars):
        rs = [r for r in res["rows"] if r["carrier"] == car]
        dirs = [d for d in ("x", "y", "z") if any(r["direction"] == d for r in rs)]
        for c, (col, ls, lab) in style.items():
            for d in dirs:                                   # [0030] 每个方法按 x/y(/z) 各画一条
                pts = [(r["T_K"], r[c]) for r in rs if r["direction"] == d and r[c]]
                if pts:
                    ax.plot(*zip(*pts), color=col, ls=ls, ms=3,
                            marker=("o" if d == "x" else "s"),
                            alpha=(1.0 if d == "x" else 0.55),
                            label="%s %s" % (lab, d))
        ax.set_yscale("log")
        ax.set_xlabel("T (K)")
        ax.set_ylabel("mobility (cm$^2$/Vs)")
        ax.set_title(car)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / MU_T_PNG, dpi=130)
    plt.close(fig)


def make_figure(out, am, bt, dpt, am2=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def split(sig, y):
        n = [(abs(sig[i]), y[i]) for i in range(len(sig))
             if sig[i] < 0 and y[i] == y[i]]
        p = [(abs(sig[i]), y[i]) for i in range(len(sig))
             if sig[i] > 0 and y[i] == y[i]]
        return sorted(n), sorted(p)

    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9))
    # (0,0) Seebeck
    a = ax[0, 0]
    if am:
        n, p = split(am["doping_signed"], am["seebeck"])
        if n: a.plot(*zip(*[(x, abs(y)) for x, y in n]), "o-", c="C0", label="amset n")
        if p: a.plot(*zip(*[(x, abs(y)) for x, y in p]), "s-", c="C1", label="amset p")
    if am2:                                   # patch_amset2d
        n, p = split(am2["doping_signed"], am2["seebeck"])
        if n: a.plot(*zip(*[(x, abs(y)) for x, y in n]), "-.", c="C2", label="amset2d n")
        if p: a.plot(*zip(*[(x, abs(y)) for x, y in p]), "-.", c="C3", label="amset2d p")
    if bt:
        n, p = split(bt["n_signed"], bt["seebeck"])
        if n: a.plot(*zip(*[(x, abs(y)) for x, y in n]), "--", c="C0", alpha=.7, label="BT2 n")
        if p: a.plot(*zip(*[(x, abs(y)) for x, y in p]), "--", c="C1", alpha=.7, label="BT2 p")
    a.set_xscale("log"); a.set_xlabel("|carrier conc| (cm$^{-3}$)")
    a.set_ylabel("|Seebeck| (uV/K)"); a.set_title("Seebeck @300K (directly comparable)")
    a.legend(fontsize=8); a.grid(alpha=.3)
    # (0,1) Lorenz
    a = ax[0, 1]
    if am:
        n, p = split(am["doping_signed"], am["lorenz"])
        if n: a.plot(*zip(*n), "o-", c="C0", label="amset n")
        if p: a.plot(*zip(*p), "s-", c="C1", label="amset p")
    if bt:
        n, p = split(bt["n_signed"], bt["lorenz"])
        if n: a.plot(*zip(*n), "--", c="C0", alpha=.7, label="BT2 n")
        if p: a.plot(*zip(*p), "--", c="C1", alpha=.7, label="BT2 p")
    a.axhline(SOMMERFELD_L, ls=":", c="k", label="Sommerfeld 2.44e-8")
    a.set_xscale("log"); a.set_xlabel("|carrier conc| (cm$^{-3}$)")
    a.set_ylabel("Lorenz (W.Ohm/K^2)"); a.set_title("Lorenz @300K (directly comparable)")
    a.legend(fontsize=8); a.grid(alpha=.3)
    # (0,2) sigma: amset absolute + BT2 per-tau (right axis)
    a = ax[0, 2]
    if am:
        n, p = split(am["doping_signed"], am["sigma"])
        if n: a.plot(*zip(*n), "o-", c="C0", label="amset n (abs)")
        if p: a.plot(*zip(*p), "s-", c="C1", label="amset p (abs)")
        a.set_ylabel("amset sigma (S/m)")
    if am2:                                   # patch_amset2d
        n, p = split(am2["doping_signed"], am2["sigma"])
        if n: a.plot(*zip(*n), "-.", c="C2", label="amset2d n (abs)")
        if p: a.plot(*zip(*p), "-.", c="C3", label="amset2d p (abs)")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("|carrier conc| (cm$^{-3}$)")
    a.set_title("sigma @300K (amset abs; BT2 dashed=sigma/tau, shape only)")
    if bt:
        a2 = a.twinx()
        n, p = split(bt["n_signed"], bt["sigma_over_tau"])
        if n: a2.plot(*zip(*n), "--", c="C0", alpha=.6)
        if p: a2.plot(*zip(*p), "--", c="C1", alpha=.6)
        a2.set_yscale("log"); a2.set_ylabel("BT2 sigma/tau (S/m/s)")
    a.legend(fontsize=8); a.grid(alpha=.3)
    # (1,0) kappa_e 三方：amset真值 / BT2用DPT-τ反推 / DPT用WF估
    a = ax[1, 0]
    if am:
        n, p = split(am["doping_signed"], am["kappa_e"])
        if n: a.plot(*zip(*n), "o-", c="C0", label="amset n (true)")
        if p: a.plot(*zip(*p), "s-", c="C1", label="amset p (true)")
    if am2:                                   # patch_amset2d
        n, p = split(am2["doping_signed"], am2["kappa_e"])
        if n: a.plot(*zip(*n), "-.", c="C2", label="amset2d n")
        if p: a.plot(*zip(*p), "-.", c="C3", label="amset2d p")
    if bt and dpt:                       # BT2 κ_e = (κ_e/τ)_BT × τ_DPT
        for carr, col, nm, neg in (("electron", "C0", "n", True),
                                   ("hole", "C1", "p", False)):
            tau = _dpt_tau_s(dpt, carr)
            if tau is None:
                continue
            pts = sorted((abs(bt["n_signed"][i]), bt["kappa_e_over_tau"][i] * tau)
                         for i in range(len(bt["n_signed"]))
                         if (bt["n_signed"][i] < 0) == neg
                         and bt["kappa_e_over_tau"][i] == bt["kappa_e_over_tau"][i])
            if pts:
                a.plot(*zip(*pts), "--", c=col, alpha=.7, label="BT2 %s (DPTτ)" % nm)
    if dpt:                              # DPT κ_e = L·σ·T (WF), σ=n e μ_DPT
        import numpy as np
        ntar = np.logspace(17, 21, 25)
        for carr, col, nm in (("electron", "C0", "n"), ("hole", "C1", "p")):
            mu = dpt.get(carr)
            if isinstance(mu, (int, float)) and mu > 0:
                sig = ntar * 1e6 * E_C * (mu * 1e-4)
                a.plot(ntar, SOMMERFELD_L * sig * 300.0, ":", c=col, alpha=.6,
                       label="DPT %s (WF)" % nm)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("|carrier conc| (cm$^{-3}$)")
    a.set_ylabel("kappa_e (W/m/K)"); a.set_title("kappa_e @300K (3-way, see notes)")
    a.legend(fontsize=7); a.grid(alpha=.3)
    # (1,1) mobility: amset total(实线) + amset-ADP(点划,与DPT同为纯声学) + DPT(虚线水平)
    a = ax[1, 1]
    if am2:                                   # patch_amset2d：二维核一路
        n, p = split(am2["doping_signed"], am2["mobility"])
        if n: a.plot(*zip(*n), "-.", c="C2", marker="^", label="amset2d n (total)")
        if p: a.plot(*zip(*p), "-.", c="C3", marker="v", label="amset2d p (total)")
        na, pa = split(am2["doping_signed"], am2["mobility_adp"])
        if na: a.plot(*zip(*na), ":", c="C2", alpha=.9, label="amset2d n (ADP-only)")
        if pa: a.plot(*zip(*pa), ":", c="C3", alpha=.9, label="amset2d p (ADP-only)")
    if am:
        n, p = split(am["doping_signed"], am["mobility"])
        if n: a.plot(*zip(*n), "o-", c="C0", label="amset n (total)")
        if p: a.plot(*zip(*p), "s-", c="C1", label="amset p (total)")
        if any(v == v for v in am.get("mobility_adp", [])):   # 有 ADP 分量才画
            na, pa = split(am["doping_signed"], am["mobility_adp"])
            if na: a.plot(*zip(*na), "-.", c="C0", alpha=.8, label="amset n (ADP-only)")
            if pa: a.plot(*zip(*pa), "-.", c="C1", alpha=.8, label="amset p (ADP-only)")
    if dpt:
        if isinstance(dpt.get("electron"), (int, float)):
            a.axhline(dpt["electron"], ls="--", c="C0", label="DPT n (ADP,300K)")
        if isinstance(dpt.get("hole"), (int, float)):
            a.axhline(dpt["hole"], ls="--", c="C1", label="DPT p (ADP,300K)")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("|carrier conc| (cm$^{-3}$)")
    a.set_ylabel("mobility (cm^2/Vs)"); a.set_title("Mobility @300K (amset vs DPT-ADP)")
    a.legend(fontsize=8); a.grid(alpha=.3)
    # (1,2) notes / provenance
    a = ax[1, 2]; a.axis("off")
    notes = (
        "Mobility comparison key:\n"
        "  DPT is ADP-only (acoustic).\n"
        "  Compare DPT vs amset ADP-only\n"
        "  (dash-dot), NOT vs amset total.\n"
        "  amset total = ADP+IMP+POP\n"
        "  (+PIE/ODP if enabled).\n"
        "  If total << ADP-only, IMP/POP\n"
        "  suppress mobility (esp. polar).\n"
        "\n"
        "kappa_e provenance:\n"
        "  amset : true (full scattering).\n"
        "  BT2   : (kappa_e/tau)*tau_DPT.\n"
        "  DPT   : WF, L*sigma*T (rough).\n"
        "\n"
        "Caveat: the 'amset' curves (S8)\n"
        "use the 3D Frohlich kernel ->\n"
        "its total may be over-suppressed\n"
        "for polar 2D; 'amset2d' (S8.4)\n"
        "uses the 2D kernel.\n"
        "S & Lorenz (tau-free) are the\n"
        "clean amset-vs-BT2 check."
    )
    a.text(0.0, 0.98, notes, va="top", ha="left", fontsize=8.5,
           family="monospace", transform=a.transAxes)

    fig.suptitle("Electronic transport comparison @ ~300 K "
                 "(amset / BoltzTraP2-CRTA / DPT)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out / "comparison_300K.png", dpi=130)
    plt.close(fig)


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    dim = _read_dim(cwd)
    is_2d = (dim == "2d")
    print("[..] 维度 DIM=%s → amset 张量约化 %s"
          % (dim or "未知(按3D处理)",
             "面内(xx+yy)/2" if is_2d else "对角(xx+yy+zz)/3"))
    am, bt, dpt = load_amset(cwd, is_2d=is_2d), load_bt2(cwd), load_dpt(cwd)
    # patch_amset2d：step8.4_amset2d 存在才读（不进 needs，3D 项目缺它不会卡）
    am2 = None
    if (cwd / AMSET2D_DIR / "transport.json").is_file():
        am2 = load_amset(cwd, is_2d=is_2d, dirname=AMSET2D_DIR)
        if am2:
            print("[..] 读到 amset2d（step8.4_amset2d，二维散射核）——一并对比")
        else:
            print("[WARN] %s/transport.json 存在但解析失败" % AMSET2D_DIR)
    if not any([am, bt, dpt]):
        sys.exit("[ERROR] step8/8.1/8.2 三个结果都读不到——先把它们跑出来再做对比。")
    if am and am.get("aniso") and am["aniso"]["anisotropic"]:
        a = am["aniso"]
        print("[WARN] 2D 面内各向异性：S/σ/κ_e 的 xx-yy 最大相对差 %.1f%% (>%.0f%%)。"
              % (a["max"] * 100, ANISO_TOL * 100))
        print("       面内平均会抹掉方向信息——CSV 已附 *_xx/*_yy 分方向列，ZT 建议分方向算。")
    got = [n for n, v in (("amset", am), ("amset2d", am2), ("BoltzTraP2", bt),
                          ("DPT", dpt)) if v]
    print("[..] 读到：%s；在 %.0f K 附近对比" % ("、".join(got), TARGET_T))
    rows = build_table(am, bt, dpt, am2=am2)
    write_table(out, rows, am, bt, dpt, am2=am2)
    try:
        make_figure(out, am, bt, dpt, am2=am2)
    except Exception as e:
        print("[WARN] 画图失败（%s: %s）——表已生成，图跳过" % (type(e).__name__, e))
    try:                                              # [0020] μ(T)，尽力而为，不影响完成标记
        write_mu_vs_T(out, cwd, is_2d)
    except Exception as e:                            # noqa: BLE001
        print("[WARN] μ(T) 表失败（%s: %s）——300 K 对比不受影响" % (type(e).__name__, e))
    for _l in resolve_kappa_L(cwd)[2]:                # patch_kl_auto
        print("[..] " + _l)

    # ---- patch_overlap_guard：重叠比值报警（VERIFICATION V23 / V27）----------------
    # 不依赖根因：重叠满足 |I|^2 <= 1，因此 unity 与真实重叠之间必有先验约束 ——
    #   ADP：绿区 [1, N_v]（N_v = 带边几个 kT 内的等价能谷数）
    #        ratio < 1     -> 红（与 |I|^2<=1 矛盾；拦截，不得进 zT）
    #        ratio > N_v   -> 黄（放行但要求人工确认；实测 N_ch 仅作解释，不作放行阈值）
    #   POP：只记录比值，不参与红/黄（真实重叠也压低谷内通道，期望本就不该贴近 1）
    # 生产默认只跑 unity_overlap，所以正常情况下这里会打印"跳过"。
    try:
        import overlap_ratio_guard as _org
        _runs = _org.find_runs(cwd)
        _u = [r for r in _runs if r["unity_overlap"] is True]
        _r = [r for r in _runs if r["unity_overlap"] is False]
        if _u and _r:
            # V27：valley_ratio.json（{"n_ch": ...}）里的 N_ch 现在只作**解释打印**，
            # 不再当放行阈值（SS 实测 N_ch=23.3 太大，用它放行等于没有红线）。
            _nch = None
            for _cand in (Path(cwd) / "valley_ratio.json",
                          Path(_u[0]["dir"]) / "valley_ratio.json",
                          Path(_r[0]["dir"]) / "valley_ratio.json"):
                try:
                    if _cand.is_file():
                        _nch = float(json.loads(_cand.read_text())["n_ch"])
                        break
                except Exception:
                    _nch = None
            if _nch:
                print("[重叠报警] 读到按谷区分实测 N_ch=%.2f（valley_ratio.json，仅作解释）" % _nch)
            # nv 显式传 2（文档默认；不传会落到 check_pair 的兜底 4.0）
            _v, _lines = _org.check_pair(_u[0], _r[0], nv=2.0, valley_ratio=_nch)
            for _l in _lines:
                print("[重叠报警] " + _l)
            if _v == "red":
                print("[重叠报警] ★ 标红：ADP 比值 < 1，与 |I|^2<=1 矛盾，**不得进入 zT 汇总**。")
            elif _v == "yellow":
                print("[重叠报警] **黄色**：ADP 比值超出 N_v，**要求人工确认**后再采纳"
                      "（实测 N_ch 仅作解释，不自动放行）。")
            elif _v == "ok":
                print("[重叠报警] ADP 比值在 [1, N_v] 内。")
        elif _r and not _u:
            # 只跑了真实重叠、没有 unity 对照 -> 按维数分档（V25）：
            #   2D：[0021] 只提示（以前判红，要求 unity 重跑 —— 与 09-26 起的出厂默认矛盾）
            #   3D：黄色（Si 全网格对照已出：ADP/overall 差 3~9%、IMP 逐机制高 1.3~1.7 倍）
            if overlap_grade_only_real(dim) == "3d":
                print("[重叠报警] **黄色**：三维 + 真实重叠（无 unity 对照）——"
                      "Si 全网格对照（V26）：ADP/overall 只差 3~9%，但 IMP 逐机制偏高 1.3~1.7 倍。"
                      "**采纳前需人工确认逐机制值**。")
            else:
                # [0021] 2026-09-26 起 2D 出厂就是真实重叠：V23 的根因（无反演体系去对称化算坏重叠）由 S8.4 的
                #   overlap_preflight 按反演判据处理 —— 全网格 h5（WAVEFUNCTION_FULL）或 IBZ + 相位补丁，不安全就拦截。
                #   以前这里一律判红、要求 unity 重跑，与出厂默认矛盾（unity 的 μ 系统性偏低，只作数量级筛选）。
                print("[重叠报警] 2D + 真实重叠（出厂默认，无 unity 对照）：不做比值报警。去对称化风险由 S8.4 的"
                      "overlap_preflight 按反演判据把关（全网格 h5 或 IBZ + 相位补丁，不安全就拦截）；"
                      "要核对就看 S8.4 生成日志里 WAVEFUNCTION_FULL / 真实重叠那几行。")
        else:
            print("[..] 重叠报警：只找到 %d 个 unity / %d 个真实重叠运行，跳过。"
                  % (len(_u), len(_r)))
    except Exception as _e:
        print("[WARN] 重叠报警检查失败（%s: %s）——不影响主表" % (type(_e).__name__, _e))

    print("[DONE] %s：comparison_300K.png / .csv / summary.txt 已生成（μ(T)：%s）" % (OUTDIR_NAME, MU_T_CSV))


if __name__ == "__main__":
    main()