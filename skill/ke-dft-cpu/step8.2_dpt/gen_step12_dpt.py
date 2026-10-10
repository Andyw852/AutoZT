#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step12_dpt.py —— 形变势理论 DPT（Bardeen-Shockley）迁移率（step8.2_dpt）。

run:gen 步骤：登录节点直接跑、不提交 SLURM，秒级。
用经典形变势理论算声学支限制的载流子迁移率——这是 2D 输运文献最主流的做法，
用来和 step8 的 amset(ADP) 做对标/交叉验证。

公式（各向同性）：
  2D:  μ = e ℏ³ C_2D / (k_B T m* m_d E1²)          （Bardeen-Shockley 2D）
  3D:  μ = 2√(2π) e ℏ⁴ C_3D / (3 (k_B T)^{3/2} m*^{5/2} E1²)
其中 C 为弹性模量、m* 有效质量、E1 形变势。

输入来源（都可手填覆盖：写**本材料的 step.conf**，见下方 SPEC；自动值仅尽力而为、务必核对）：
  · C   ← step6_elastic/OUTCAR（可靠自动；2D 再乘层厚得 C_2D[N/m]）
  · m*  ← step3_uniform 能带曲率（BoltzTraP2 尽力自动）
  · E1  ← step7b_deform_read/deformation.h5 的带边形变势（尽力自动）

【重要诚实说明】m* 与 E1 的自动提取精度有限。amset 的 ADP 已经用 deformation.h5 做了
更严格的形变势散射；本步是"经典闭式对标"。要发表级 DPT，建议手填 m*/E1/C
（多数 DPT 论文就是这么做的），或核对自动值后再用。手填写本材料的 step.conf：
    tf -p <材料> -j step8.2_dpt conf --set M_EFF_ELECTRON=0.030
不要改本脚本的 MANUAL：本脚本每次 gen 都原样推到每个材料，改了就作用到所有材料（V148）。

产出（done_marker）：本步目录下的 dpt_result.json（+ summary.txt）。缺输入时写出
带指引的部分结果，绝不崩溃、绝不编造。
"""
import json
import math
import os
import sys
from pathlib import Path

# =========================== 可改参数区 ===========================
# [SKILL_REV] 版本戳：写进 dpt_result.json，铺开时验证跑的是哪份 skill 副本。
# 每次改本脚本逻辑后更新（如 "2026-08-29-nstep-linear"）。
_SKILL_REV = "2026-10-10-0028-spin"
# [R-SCAN] 强制选点壳层 NSTEP（None=自动 2/3/4/5；填 2/3/4/5=只试该值，做 R-scan 用）。
# 用法：tf -p <材料> -j step8.2_dpt conf --set FORCE_NSTEP=3 → retry + start，对比 dpt_result.json 的 m_d
#   与 m_provenance 里的 R（V148 起写 step.conf，不改本脚本）。
FORCE_NSTEP = None
OUTDIR_NAME = "step8.2_dpt"
UNIFORM_DIR = "step3_uniform"
ELASTIC_DIR = "step6_elastic"
DEFORM_READ_DIR = "step7b_deform_read"
DEFORM_DIR = "step7_deform"          # [0022] 逐谷 E1 读 S7 形变构型的 vasprun（S7.1 读的同一份）
AMSET_DIR   = "step8_amset"          # 读 2d_correction.json 拿层厚 t、c/t

TEMPERATURE_K = 300.0                # DPT 迁移率报此温度（μ∝1/T，可换算）
# [0020] 多温度：mobility_vs_T 按 μ∝1/T（2D）/ T^-1.5（3D）从 TEMPERATURE_K 的值精确换算（m*/E1/C 与 T 无关）。
#   默认与 S8/S8.4 的 AMSET 网格（100:900:9）、S8.1 的 TEMPERATURES 一致，S8.3 才能逐温度对齐；
#   只要 100–700 K 就在本步 step.conf 写 TEMPERATURES = 100:700:7。格式 "起:止:个数" 或逗号列表。
TEMPERATURES = "100:900:9"
CARRIER = "both"                     # electron / hole / both
# [0022] 多谷 DPT（dpt_result.json 的 multi_valley，2D）：自动找离带边 VALLEY_WINDOW_EV 以内的所有谷（K/Q/Γ/M…），
#   逐谷质量张量与 E1、按玻尔兹曼布居加权，每个温度单独算。带边在 Q 这类非高对称谷时单谷值拒绝出数，这里照样能算。
#   关掉：step.conf 写 DPT_VALLEYS = off。窗口：VALLEY_WINDOW_EV（eV，默认 0.30；900 K 时 e^(−0.3/kT)≈2%）。
DPT_VALLEYS = "auto"
VALLEY_WINDOW_EV = 0.30
# [0023] 有贡献的谷（任一温度布居 ≥ VALLEY_W_MIN）之间等效形变势相差 > VALLEY_E1_RATIO_MAX 倍 -> 多谷拒绝出数。
#   实跑四个 TMD：电子 Q 谷 E1 只有 0.5–0.8 eV、K 谷 7.4–9.4 eV（真空口径），τ∝1/E1² 让 Q 谷 τ 长两个数量级，
#   多谷 μ 被放大 7× 到 6×10⁵×。E1 小说明该谷的 LA 膨胀形变势耦合几乎抵消，它的输运由 DPT 不含的
#   TA/谷间/光学散射决定（Jin et al. PRB 90, 045422：DFPT 拟合的 Q/K 声学形变势比 0.6–0.9；WSe2 电子含谷间 30、
#   只 K 250 cm²/Vs）——不是 E1 取法能修的，DPT 在这里失效。3 倍 = 两谷 τ 差一个数量级。
VALLEY_E1_RATIO_MAX = 3.0
VALLEY_W_MIN = 0.01
# E1 自动提取口径开关（手填 MANUAL/MANUAL_ANISO 时此开关失效）：
#   "vac"   = 真空对齐 E1_vac（文献 Eq.9 dE_edge/dγ 口径，仅 2D 且 LOCPOT 可用）
#   "amset" = amset h5 平均芯势口径（⟨|D|⟩）
# [C5] 符号约定：两套口径写出的字段【都恒非负】——vac 侧生成端已做
#   vals.append(abs(raw - dvac))，amset 侧 h5 本身是 ⟨|D|⟩。消费端仍补 abs()
#   兜底，但不要指望 E1_vac_* 带符号（μ 只依赖 E1²，符号无物理用途）。
# [C2] E1_SOURCE 是【硬选择】不是偏好：要 "vac" 却拿不到 E1_vac_* 时直接报错，
#   绝不静默回退 amset 口径（两者数值能差 3 倍，μ 差一个量级）。
E1_SOURCE = "vac"
# [V174] m* 直接拟合 step3_uniform 的网格点：S3 面内网格不合现行规则（2D 0.05 / 3D 0.06 Å⁻¹）就拒绝生成。
#   CrS₂_hex 15×15×1 -> 48×48×3：m* 0.967/1.009 -> 0.866/0.883，DPT +25%/+31%（拟合窗口 4×最近邻 = 0.64 Å⁻¹，
#   吃进了非抛物区；V147 的最近一壳 3.85 kT 刚好没过 4 kT 的线）。只为复现旧结果才写 step.conf ALLOW_COARSE_S3 = true。
ALLOW_COARSE_S3 = False

# —— 手动覆盖（填了就用手填值，最可靠；None=尝试自动）——
# ★ V148：**只给一个材料用就写它的 step.conf**（下方 SPEC 的键），不要改这里。这里的值作用到所有材料
#   （本脚本每次 gen 都从 skill 原样推到每个材料目录），不是全 None 时每次 gen 都 ★ 告警。
MANUAL = {
    "m_eff_electron": None,   # m*/m0（各向同性；也可只填电子或空穴）
    "m_eff_hole":     None,
    "E1_electron_eV": None,   # 形变势 eV
    "E1_hole_eV":     None,
    "C_2D_N_per_m":   None,   # 2D 弹性模量 N/m（填了就不从 step6 推）
    "C_3D_GPa":       None,   # 3D 弹性模量 GPa
    "thickness_A":    None,   # 2D 层厚 Å（None=读 2d_correction.json）
}

# [V148] 材料级覆盖：本材料的 step.conf（tf -p <材料> -j step8.2_dpt conf --set 键=值）。
#   起因：GaAs 要拿 amset eff-mass 的 m*=0.030 对照 ADP/DPT，只能改共享脚本的 MANUAL，改着的那段时间里
#   任何材料的 S8.2 都会用上 GaAs 的 0.030（gen 脚本总是从 skill 覆盖推送）。step.conf 的值优先于脚本。
STEP = OUTDIR_NAME
SPEC = {
    "M_EFF_ELECTRON": (None, "float"),    # m*/m0
    "M_EFF_HOLE":     (None, "float"),
    "E1_ELECTRON_EV": (None, "float"),    # 形变势 eV
    "E1_HOLE_EV":     (None, "float"),
    "C_2D_N_PER_M":   (None, "float"),    # 2D 弹性模量 N/m
    "C_3D_GPA":       (None, "float"),    # 3D 弹性模量 GPa
    "THICKNESS_A":    (None, "float"),    # 2D 层厚 Å
    "E1_SOURCE":      (None, "str"),      # vac / amset
    "FORCE_NSTEP":    (None, "int"),      # 2/3/4/5
    "ALLOW_COARSE_S3": (False, "bool"),   # [V174] S3 面内网格不合现行规则也照算
    "TEMPERATURE_K":  (None, "float"),    # [0020] 主报告温度（K），默认 300
    "TEMPERATURES":   (None, "str"),      # [0020] μ(T) 网格："100:900:9" 或 "100,300,600"
    "DPT_VALLEYS":    (None, "str"),      # [0022] auto / off
    "VALLEY_WINDOW_EV": (None, "float"),  # [0022] 多谷窗口（eV）
    "VALLEY_E1_RATIO_MAX": (None, "float"),  # [0023] 贡献谷之间 E1 最大/最小 超过它就拒绝
}
_CONF_TO_MANUAL = {"M_EFF_ELECTRON": "m_eff_electron", "M_EFF_HOLE": "m_eff_hole",
                   "E1_ELECTRON_EV": "E1_electron_eV", "E1_HOLE_EV": "E1_hole_eV",
                   "C_2D_N_PER_M": "C_2D_N_per_m", "C_3D_GPA": "C_3D_GPa", "THICKNESS_A": "thickness_A"}
_MANUAL_SRC = {}                          # MANUAL 键 -> "step.conf" / "脚本"
# =================================================================

# 物理常数（SI）
E_C   = 1.602176634e-19
HBAR  = 1.054571817e-34
KB    = 1.380649e-23
M0    = 9.1093837015e-31

_STEP1_CANDS = ("step1_opt", "step1_std_opt",
                "step1c_PBE_opt", "step1b_PBE_opt", "step1a_PBE_opt")


def _exc_brief(e):
    """[V150] 兜底 except 的说明：类型 + 信息 + 本脚本里出错的那一行（以前只写类型名，WS2 的
    "二次型拟合异常：IndexError" 无从定位）。"""
    import traceback
    msg = "%s: %s" % (type(e).__name__, str(e)[:200]) if str(e) else type(e).__name__
    here = [f for f in traceback.extract_tb(e.__traceback__)
            if Path(f.filename).name == Path(__file__).name]
    if here:
        f = here[-1]
        msg += "（%s:%d `%s`）" % (Path(f.filename).name, f.lineno, (f.line or "").strip()[:120])
    return msg


def _manual_prov(key):
    """[V148] 手填值的来源：step.conf（本材料）/ 脚本 MANUAL（所有材料）。"""
    return {"step.conf": "manual(step.conf)",
            "脚本": "manual(脚本 MANUAL，作用到所有材料)"}.get(_MANUAL_SRC.get(key), "manual")


def apply_conf(cwd):
    """[V148] 读本材料 step.conf 的覆盖写进 MANUAL / E1_SOURCE / FORCE_NSTEP；脚本里 MANUAL 被改过就 ★ 告警。
    返回 {键: {"value", "source"}}（写进 dpt_result.json 的 overrides）。"""
    global E1_SOURCE, FORCE_NSTEP, ALLOW_COARSE_S3, TEMPERATURE_K, TEMPERATURES, DPT_VALLEYS, VALLEY_WINDOW_EV
    global VALLEY_E1_RATIO_MAX
    rec = {}
    for k, v in MANUAL.items():
        if v is not None:
            _MANUAL_SRC[k] = "脚本"
            rec[k] = {"value": v, "source": "gen_step12_dpt.py 的 MANUAL（作用到所有材料）"}
    for k, v in MANUAL_ANISO.items():
        if v is not None:
            rec[k] = {"value": list(v), "source": "gen_step12_dpt.py 的 MANUAL_ANISO（作用到所有材料）"}
    if rec:
        print("[WARN] ★ 共享脚本 gen_step12_dpt.py 的 MANUAL/MANUAL_ANISO 被改过（%s）：本脚本每次 gen 都原样推到每个材料，"
              "这些值会用到**所有**材料的 S8.2。只给一个材料用：还原脚本，改写本材料的 step.conf"
              "（tf -p <材料> -j step8.2_dpt conf --set M_EFF_ELECTRON=…）"
              % ", ".join("%s=%s" % (k, r["value"]) for k, r in sorted(rec.items())))
    if not (Path(cwd) / "step.conf").is_file():
        return rec
    try:
        import stepconf
    except ImportError:
        print("[WARN] 有 step.conf 但没有 stepconf.py：step.conf 的覆盖未生效（skill.yaml 的 gen_need 漏了它？）")
        return rec
    p = stepconf.load(SPEC, STEP, str(cwd), strict=False)
    for ck, mk in _CONF_TO_MANUAL.items():
        if p[ck] is None:
            continue
        if mk in rec:
            print("[WARN] %s：step.conf 的 %s=%s 盖过脚本 MANUAL 的 %s" % (mk, ck, p[ck], rec[mk]["value"]))
        MANUAL[mk] = float(p[ck])
        _MANUAL_SRC[mk] = "step.conf"
        rec[mk] = {"value": float(p[ck]), "source": "step.conf %s" % ck}
        print("[OK] %s = %s（step.conf %s）" % (mk, p[ck], ck))
    if p["E1_SOURCE"]:
        v = str(p["E1_SOURCE"]).strip().lower()
        if v not in ("vac", "amset"):
            sys.exit("[ERROR] step.conf 的 E1_SOURCE=%r：只能是 vac / amset" % p["E1_SOURCE"])
        E1_SOURCE = v
        rec["E1_SOURCE"] = {"value": v, "source": "step.conf E1_SOURCE"}
        print("[OK] E1_SOURCE = %s（step.conf）" % v)
    if p["FORCE_NSTEP"] is not None:
        if int(p["FORCE_NSTEP"]) not in (2, 3, 4, 5):
            sys.exit("[ERROR] step.conf 的 FORCE_NSTEP=%r：只能是 2/3/4/5" % p["FORCE_NSTEP"])
        FORCE_NSTEP = int(p["FORCE_NSTEP"])
        rec["FORCE_NSTEP"] = {"value": FORCE_NSTEP, "source": "step.conf FORCE_NSTEP"}
        print("[OK] FORCE_NSTEP = %d（step.conf）" % FORCE_NSTEP)
    if p["ALLOW_COARSE_S3"]:
        ALLOW_COARSE_S3 = True
        rec["ALLOW_COARSE_S3"] = {"value": True, "source": "step.conf ALLOW_COARSE_S3"}
        print("[OK] ALLOW_COARSE_S3 = true（step.conf）")
    if p["TEMPERATURE_K"] is not None:                    # [0020]
        if not float(p["TEMPERATURE_K"]) > 0:
            sys.exit("[ERROR] step.conf 的 TEMPERATURE_K=%r：必须 > 0" % p["TEMPERATURE_K"])
        TEMPERATURE_K = float(p["TEMPERATURE_K"])
        rec["TEMPERATURE_K"] = {"value": TEMPERATURE_K, "source": "step.conf TEMPERATURE_K"}
        print("[OK] TEMPERATURE_K = %g（step.conf）" % TEMPERATURE_K)
    if p["TEMPERATURES"]:
        try:
            parse_temps(p["TEMPERATURES"])
        except ValueError as e:
            sys.exit("[ERROR] step.conf 的 TEMPERATURES=%r：%s" % (p["TEMPERATURES"], e))
        TEMPERATURES = str(p["TEMPERATURES"]).strip()
        rec["TEMPERATURES"] = {"value": TEMPERATURES, "source": "step.conf TEMPERATURES"}
        print("[OK] TEMPERATURES = %s（step.conf）" % TEMPERATURES)
    if p["DPT_VALLEYS"]:                                  # [0022]
        v = str(p["DPT_VALLEYS"]).strip().lower()
        if v not in ("auto", "off"):
            sys.exit("[ERROR] step.conf 的 DPT_VALLEYS=%r：只能是 auto / off" % p["DPT_VALLEYS"])
        DPT_VALLEYS = v
        rec["DPT_VALLEYS"] = {"value": v, "source": "step.conf DPT_VALLEYS"}
        print("[OK] DPT_VALLEYS = %s（step.conf）" % v)
    if p["VALLEY_WINDOW_EV"] is not None:
        if not 0.0 < float(p["VALLEY_WINDOW_EV"]) <= 2.0:
            sys.exit("[ERROR] step.conf 的 VALLEY_WINDOW_EV=%r：必须在 (0, 2] eV" % p["VALLEY_WINDOW_EV"])
        VALLEY_WINDOW_EV = float(p["VALLEY_WINDOW_EV"])
        rec["VALLEY_WINDOW_EV"] = {"value": VALLEY_WINDOW_EV, "source": "step.conf VALLEY_WINDOW_EV"}
        print("[OK] VALLEY_WINDOW_EV = %g（step.conf）" % VALLEY_WINDOW_EV)
    if p["VALLEY_E1_RATIO_MAX"] is not None:                # [0023]
        if not float(p["VALLEY_E1_RATIO_MAX"]) >= 1.0:
            sys.exit("[ERROR] step.conf 的 VALLEY_E1_RATIO_MAX=%r：必须 ≥ 1" % p["VALLEY_E1_RATIO_MAX"])
        VALLEY_E1_RATIO_MAX = float(p["VALLEY_E1_RATIO_MAX"])
        rec["VALLEY_E1_RATIO_MAX"] = {"value": VALLEY_E1_RATIO_MAX, "source": "step.conf VALLEY_E1_RATIO_MAX"}
        print("[OK] VALLEY_E1_RATIO_MAX = %g（step.conf）" % VALLEY_E1_RATIO_MAX)
    return rec


def parse_temps(spec):
    """[0020] "100:900:9"（起:止:个数，与 AMSET 的 TEMPERATURES 同写法）或 "100,300,600" -> 升序 K 列表。"""
    if isinstance(spec, (list, tuple)):
        vals = [float(x) for x in spec]
    else:
        txt = str(spec).strip()
        if ":" in txt:
            parts = txt.split(":")
            if len(parts) != 3:
                raise ValueError("区间写法是 起:止:个数，例如 100:900:9")
            a, b, n = float(parts[0]), float(parts[1]), int(float(parts[2]))
            if n < 1:
                raise ValueError("个数必须 >= 1")
            vals = [a] if n == 1 else [a + (b - a) * i / (n - 1) for i in range(n)]
        else:
            vals = [float(x) for x in txt.replace(";", ",").replace(" ", ",").split(",") if x]
    vals = sorted({round(v, 6) for v in vals})
    if not vals or vals[0] <= 0:
        raise ValueError("温度必须 > 0 K")
    return vals


def _scale_T(is_2d, T0, T):
    """[0020] μ(T)/μ(T0)：2D Bardeen-Shockley ∝ 1/T，3D ∝ T^-1.5（τ = μm*/e 同比例）。"""
    return (T0 / T) if is_2d else (T0 / T) ** 1.5


def mobility_vs_T(rec, is_2d, T0, temps):
    """[0020] 由 T0 的结果换算各温度：μ（主值）、2D 分方向 μ/τ 与面内平均 (μx+μy)/2。"""
    bd = rec.get("by_direction") or {}
    aniso = bd.get("status") == "ok"
    rows = []
    for T in temps:
        f = _scale_T(is_2d, T0, T)
        mu = rec.get("mobility_cm2_Vs")
        row = {"T_K": T, "mobility_cm2_Vs": (round(mu * f, 3) if mu is not None else None)}
        if aniso:
            dirs = [d for d in ("x", "y", "z") if isinstance(bd.get(d), dict)]     # [0025] 3D 有 z
            for d in dirs:
                row["mobility_%s_cm2_Vs" % d] = round(bd[d]["mobility_cm2_Vs"] * f, 3)
                row["tau_%s_s" % d] = bd[d]["tau_s"] * f
            if is_2d:
                row["mobility_inplane_cm2_Vs"] = round(
                    (row["mobility_x_cm2_Vs"] + row["mobility_y_cm2_Vs"]) / 2.0, 3)
            else:
                row["mobility_mean_cm2_Vs"] = round(sum(row["mobility_%s_cm2_Vs" % d] for d in dirs) / len(dirs), 3)
        rows.append(row)
    return rows


def invalidate_consumers(cwd):
    """[patch_post_invalidate V148] 本步重新生成 -> S8.1（用本步的 τ）与 S8.3（对比图）的完成标记归档、重新排队。
    以前 DOWNSTREAM 表里有这两条，但没有 gen 调：GaAs 的 S8.2 重跑 4 次，S8.1 一直用 V145 之前的空穴 τ。"""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        import ke_common as _kc
    except ImportError:
        print("[WARN] 没有 ke_common.py，S8.1/S8.3 不会自动重排（skill.yaml 的 gen_need 漏了它？）")
        return []
    return _kc.invalidate_downstream(cwd, OUTDIR_NAME, "%s 重新生成（DPT 重算）" % OUTDIR_NAME)


def _read_dim(cwd):
    for name in _STEP1_CANDS:
        mf = Path(cwd) / name / "workflow_method.txt"
        if not mf.is_file():
            continue
        for ln in mf.read_text(errors="ignore").splitlines():
            if ln.strip().upper().startswith("DIM="):
                return ln.split("=", 1)[1].strip().lower()
    return None


def _guard_not_0d(cwd):
    if _read_dim(cwd) == "0d":
        sys.exit("[ERROR] step8.2_dpt 不支持 0D 体系（无能带色散，形变势无意义）。")


def _s3_grid_gate(cwd, dim):
    """[V174] S3 面内网格不合现行规则 -> 拒绝（ALLOW_COARSE_S3 时告警照算）。真空轴 kz 不影响面内 m*，不查。
    返回写进 dpt_result.json 的 {"mesh", "issues", "allow_coarse"}；没有 ke_common 返回 None。"""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        import ke_common as _kc
    except ImportError:
        print("[WARN] 没有 ke_common.py，S3 网格没查（skill.yaml 的 gen_need 漏了它？）")
        return None
    iss = _kc.s3_grid_gate(cwd, "S8.2", allow=ALLOW_COARSE_S3, dim=dim if dim in ("2d", "3d") else None,
                           check_kz=False, uses="m* 直接拟合这张网格上的点（拟合窗口随最近邻间距放大）",
                           redo="本步（E1 跟着 S7 -> S7.1 的新网格再更新一次）")
    mesh = _kc.read_kpoints_mesh(Path(cwd) / UNIFORM_DIR / "KPOINTS")
    return {"mesh": [int(x) for x in mesh] if mesh else None,
            "issues": [m for _, m in iss], "allow_coarse": bool(ALLOW_COARSE_S3)}


def _functional(cwd):
    """[V181] 本材料用的泛函（step1 定下、下游继承），写进 dpt_result.json 和摘要；没有 ke_common 返回 None。"""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        import ke_common as _kc
    except ImportError:
        return None
    f = _kc.material_functional(cwd)
    print("[..] 泛函：%s（来源：%s）" % (f["label"], f["source"] or "—"))
    return f


def _hex_cell(recip):
    """面内两个倒格矢等长、夹角 60°/120° -> 六方胞。"""
    import numpy as np
    _b1 = float(np.linalg.norm(recip[0]))
    _b2 = float(np.linalg.norm(recip[1]))
    _ang = float(np.degrees(np.arccos(max(-1.0, min(1.0, float(np.dot(recip[0], recip[1])) / (_b1 * _b2))))))
    return (abs(_b1 - _b2) < 0.05 * max(_b1, _b2)
            and (abs(_ang - 60.0) < 5.0 or abs(_ang - 120.0) < 5.0))


def _hex_offgrid(recip, kf, kgrid=None):
    """[GUARD-2026-09-23 / V152] 六方胞的带边锚在非高对称点 -> (说明, False)；否则 (None, 是否高对称点)。
    各向同性与分方向两条 m* 路径共用（V152 之前分方向路径没有这道闸：WSe2 的电子分方向给出
    x=0.48 / y=0.68 且 status ok，而 C₃ 晶体的迁移率张量面内各向同性）。
    kgrid（本路径用的 k 点）里有 K 点 -> 不是网格问题：带边在 Q/Λ 这类非高对称谷（WSe2 45×45：K=15/45 在网上，
    带边在 8/45 = 0.178，Γ–K 中点附近）；重做 S3 没用，原来的"K 离网"提示把人引到了错的方向。"""
    import numpy as np
    _hex = _hex_cell(recip)
    _kf = np.asarray(kf, float) % 1.0
    _at_hs = bool(np.all(np.abs(_kf[:2] * 6.0 - np.round(_kf[:2] * 6.0)) < 0.01))
    if not _hex or _at_hs:
        return None, _at_hs
    k_on_grid = False
    if kgrid is not None:
        g = np.asarray(kgrid, float)[:, :2] * 3.0
        on3 = np.all(np.abs(g - np.round(g)) < 1e-3, axis=1)
        k_on_grid = bool(np.any(on3 & np.any(np.abs(np.round(g) % 3) > 0.5, axis=1)))
    if k_on_grid:
        return ("六方胞带边在非高对称点 k_frac=%s，而 K 点在网格上：不是网格问题，带边多半在 Q/Λ 这类"
                "非高对称谷（6 个等价谷）。单谷 DPT 的 m*（方向分解与分方向 x/y）都不适用，重做 step3_uniform "
                "没用 —— 看 dpt_result.json 的 multi_valley（[0022] 逐谷质量张量与 E1、按布居加权；[0023] Q 与 K 的形变势"
                "悬殊时它也拒绝出数：那是 DPT 对这类带边不适用，不是缺参数）；不要手填 M_EFF_*（单个质量混了 K/Q "
                "两类谷，也不随温度变）" % np.round(_kf, 4)), False
    return ("六方胞带边锚在非高对称点 k_frac=%s（K 离网：step3_uniform 的 N 不是"
            " 3 的倍数）——请按 DK_MAX 规则用 N 为 6 的倍数重新生成 step3_uniform，"
            "或在本材料 step.conf 手填 M_EFF_*" % np.round(_kf, 4)), False


# ---------- [PATCH-KE-2026] 自旋道感知的带边定位 ----------
def _spin_arrays(v):
    """按 (up, down) 固定顺序返回各自旋道本征值 [(ispin, arr(nk,nb,2)), ...]。
    ISPIN=1 时只有一道，行为与打补丁前完全一致。"""
    import numpy as np
    eigmap = v.eigenvalues
    try:
        from pymatgen.electronic_structure.core import Spin
        order = [s for s in (Spin.up, Spin.down) if s in eigmap]
    except Exception:
        order = []
    if not order:
        order = list(eigmap.keys())
    return [(i, np.asarray(eigmap[s])) for i, s in enumerate(order)]


# [V145] 带边简并时取哪一支：价带取简并里**编号最大**的一支、导带取**编号最小**的一支 ——
#   VASP 按能量给带编号，离开带边后它们正是"最平 = 最重 = 态密度最大"的那支（重空穴 / 重电子）。
#   原来 argmax 遇到并列取第一个，即价带简并里编号最小 = 下降最快的轻空穴：GaAs（无 SOC，Γ 点三重简并）
#   的 DPT 空穴迁移率因此是 276 万 cm²/Vs。简并重数记在 _EDGE_DEGEN 里，写进 m_provenance。
DEGEN_TOL_EV = 0.005
_EDGE_DEGEN = {}
# [V147] 网格分辨不出带边曲率：最近网格点高出带边 > NN_KT_MAX·kT 时，拟合量到的是能量平均的质量，
#   带边非抛物时偏重。GaAs（PBEsol 带隙 0.418 eV）：S3 最近点高出 CBM 约 0.25 eV（~10 kT），拟合 0.0523，
#   AMSET 插值密网格（amset eff-mass）给 0.030；ADP ∝ m*^(−5/2) -> ADP/DPT 正好差 4 倍。
NN_KT_MAX = 4.0
_EDGE_NN_DE = {}


def _find_band_edge(v, carrier, occ_tol=0.5):
    """跨【所有】自旋道找全局 CBM/VBM。

    打补丁前用的是 list(v.eigenvalues.values())[0]，ISPIN=2 时永远只看
    spin-up，spin-down 的带边被丢掉 —— 对 FM 半导体/半金属会直接取错带。
    返回 (ispin, k0, b0, ene(nk,nb))；全占据/全空(金属)时返回 None。"""
    import numpy as np
    best = None
    for isp, eig in _spin_arrays(v):
        ene, occ = eig[:, :, 0], eig[:, :, 1]
        mask = (occ < occ_tol) if carrier == "electron" else (occ > occ_tol)
        if not mask.any():
            continue
        if carrier == "electron":
            cand = np.where(mask, ene, np.inf)
            k0, b0 = np.unravel_index(np.argmin(cand), cand.shape)
            val, better = cand[k0, b0], (best is None or cand[k0, b0] < best[0])
        else:
            cand = np.where(mask, ene, -np.inf)
            k0, b0 = np.unravel_index(np.argmax(cand), cand.shape)
            val, better = cand[k0, b0], (best is None or cand[k0, b0] > best[0])
        if better:
            best = (val, isp, int(k0), int(b0), ene)
    if best is None:
        return None
    val, isp, k0, b0, ene = best
    deg = [b for b in range(ene.shape[1]) if abs(float(ene[k0, b]) - float(val)) < DEGEN_TOL_EV]
    _EDGE_DEGEN[carrier] = len(deg) if deg else 1
    if len(deg) > 1:
        b0 = max(deg) if carrier == "hole" else min(deg)
    return isp, k0, int(b0), ene


def _sorted_dp_keys(f):
    """deformation.h5 里的 deformation_potentials 数据集，up 在前 down 在后。"""
    keys = [k for k in f.keys() if "deformation_potentials" in k]
    return sorted(keys, key=lambda s: ("down" in s.lower(), s))


# ---------- C：从 step6_elastic/OUTCAR 读面内弹性 ----------
def read_elastic_inplane_GPa(cwd):
    """返回 (C11, C22) GPa 面内对角；读不到返回 (None, None)。"""
    oc = Path(cwd) / ELASTIC_DIR / "OUTCAR"
    if not oc.is_file():
        return None, None
    import re
    txt = oc.read_text(errors="ignore")
    i = txt.rfind("TOTAL ELASTIC MODULI")
    if i < 0:
        return None, None
    rows = []
    for ln in txt[i:].splitlines()[1:]:
        nums = re.findall(r"-?\d+\.\d+", ln)
        if len(nums) >= 6:
            rows.append([float(x) for x in nums[:6]])
        if len(rows) == 6:
            break
    if len(rows) < 6:
        return None, None
    return rows[0][0] / 10.0, rows[1][1] / 10.0   # kBar->GPa, C11 C22


def read_elastic_diag_GPa(cwd):
    """[0025] step6_elastic/OUTCAR 的 TOTAL ELASTIC MODULI 对角 (C11, C22, C33)，GPa；读不到返回 None。
    VASP 的顺序是 XX YY ZZ XY YZ ZX，所以第三行第三列就是 C33。"""
    oc = Path(cwd) / ELASTIC_DIR / "OUTCAR"
    if not oc.is_file():
        return None
    import re
    txt = oc.read_text(errors="ignore")
    i = txt.rfind("TOTAL ELASTIC MODULI")
    if i < 0:
        return None
    rows = []
    for ln in txt[i:].splitlines()[1:]:
        nums = re.findall(r"-?\d+\.\d+", ln)
        if len(nums) >= 6:
            rows.append([float(x) for x in nums[:6]])
        if len(rows) == 6:
            break
    if len(rows) < 6:
        return None
    return rows[0][0] / 10.0, rows[1][1] / 10.0, rows[2][2] / 10.0   # kBar -> GPa


def _thickness_A(cwd):
    if MANUAL["thickness_A"]:
        return float(MANUAL["thickness_A"])
    p = Path(cwd) / AMSET_DIR / "2d_correction.json"
    try:
        rec = json.loads(p.read_text())
        return float((rec.get("layer_thickness") or {}).get("thickness_used_A"))
    except (OSError, KeyError, ValueError, TypeError):
        return None


def _cell_height_A(cwd):
    """2D 元胞垂直层面的高度 h = V/A（体积÷面内面积），对非正交/c 倾斜也正确。
    从 POSCAR/CONTCAR 的三个格矢直接算；读不到再退回 2d_correction.json 的 cell_c_A。"""
    for name in ("step3_uniform", ELASTIC_DIR, "step8_amset"):
        for fn in ("POSCAR", "CONTCAR"):
            f = Path(cwd) / name / fn
            if not f.is_file():
                continue
            try:
                import numpy as np
                ln = f.read_text().splitlines()
                scale = float(ln[1].split()[0])
                a = np.array([float(x) for x in ln[2].split()[:3]]) * scale
                b = np.array([float(x) for x in ln[3].split()[:3]]) * scale
                c = np.array([float(x) for x in ln[4].split()[:3]]) * scale
                V = abs(np.dot(a, np.cross(b, c)))     # 体积
                A = np.linalg.norm(np.cross(a, b))     # 面内面积 |a×b|
                if A > 0 and V > 0:
                    return float(V / A)                # 垂直高度 h
            except (IndexError, ValueError):
                pass
    # 退回：json 里的 cell_c_A（注意它是 |c|，仅当 c 垂直层面时才等于 h）
    p = Path(cwd) / AMSET_DIR / "2d_correction.json"
    try:
        return float(json.loads(p.read_text())["cell_c_A"])
    except (OSError, KeyError, ValueError, TypeError):
        return None


def get_C(cwd, is_2d):
    """返回 (C_value, C_unit, provenance)。2D 给 C_2D[N/m]，3D 给 C_3D[Pa]。
    2D 刚度 C_2D[N/m] = C_3D[Pa] × 元胞垂直高度 h=V/A —— VASP 的 C_3D 按整胞体积
    A·h 归一化，乘 h 得单位面积的 2D 模量（与厚度歧义无关，2D-DPT 标准定义，
    如 Qiao 等黑磷工作）。用 V/A 而非 |c|，对非正交/c 倾斜的胞也正确。"""
    if is_2d:
        if MANUAL["C_2D_N_per_m"]:
            return float(MANUAL["C_2D_N_per_m"]), "N/m", _manual_prov("C_2D_N_per_m")
        c11, c22 = read_elastic_inplane_GPa(cwd)
        h = _cell_height_A(cwd)
        if c11 is None or c22 is None or not h:
            return None, "N/m", "缺 step6 弹性或元胞高度 h=V/A"
        C_inplane_Pa = (c11 + c22) / 2.0 * 1e9
        C_2D = C_inplane_Pa * (h * 1e-10)          # Pa·m = N/m，用 h=V/A
        return C_2D, "N/m", "step6 面内(C11+C22)/2 × h(V/A)=%.3fÅ" % h
    else:
        if MANUAL["C_3D_GPa"]:
            return float(MANUAL["C_3D_GPa"]) * 1e9, "Pa", _manual_prov("C_3D_GPa")
        # [0025] 3D 取三个对角的平均：原来是 (C11+C22)/2（照搬 2D），非立方体系（六方、正交）漏了 C33
        cd = read_elastic_diag_GPa(cwd)
        if cd is None:
            return None, "Pa", "缺 step6 弹性"
        return sum(cd) / 3.0 * 1e9, "Pa", "step6 (C11+C22+C33)/3"


# ---------- m*：能带边抛物拟合（2D/3D 通用，全自动） ----------
def get_effective_mass(cwd, carrier, is_2d):
    """能带边抛物拟合 m*：E=C|Δk|² → m*/m0=3.80998/|C|（k 用倒格矢含2π，rad/Å）。
    2D 只取面内。返回 (m*/m0, provenance)；失败 (None, 原因)。"""
    key = "m_eff_electron" if carrier == "electron" else "m_eff_hole"
    if MANUAL[key]:
        return float(MANUAL[key]), _manual_prov(key)
    vr = None
    for n in ("vasprun.xml", "vasprun.xml.gz"):
        p = Path(cwd) / UNIFORM_DIR / n
        if p.is_file():
            vr = p
            break
    if vr is None:
        return None, "缺 step3_uniform/vasprun.xml"
    try:
        import numpy as np
        from pymatgen.io.vasp import Vasprun
        v = Vasprun(str(vr), parse_dos=False, parse_potcar_file=False)
        recip = v.final_structure.lattice.reciprocal_lattice.matrix  # rad/Å, 含2π
        kfrac = np.array(v.actual_kpoints)                            # (nk,3)
        kcart = kfrac @ recip                                         # (nk,3) rad/Å
        hit = _find_band_edge(v, carrier)        # [PATCH-KE-2026] 跨自旋道
        if hit is None:
            return None, "无未占据/占据态 —— 体系是金属？（带边不存在，m* 无意义）"
        isp, k0, b0, ene = hit
        e0 = ene[k0, b0]
        eband = ene[:, b0]                        # 该带所有 k 的能量
        dk = kcart - kcart[k0]                    # 相对带边 Δk（rad/Å）
        de = eband - e0                           # ΔE（eV）
        # 2D：只取面内（|Δkz| 小），用 x,y 分量拟合
        if is_2d:
            zmask = np.abs(dk[:, 2]) < 0.02
            dk_use = dk[zmask][:, :2]
            de_use = de[zmask]
        else:
            dk_use = dk
            de_use = de
        q2 = np.sum(dk_use**2, axis=1)
        nz = q2 > 1e-9
        q2n, den = q2[nz], de_use[nz]
        if len(q2n) < 3:
            return None, "带边面内点太少（网格太粗）"
        # [V147] 最近一壳网格点离带边多高（中位数）
        _shell = q2n <= float(q2n.min()) * 1.05
        _EDGE_NN_DE[carrier] = float(np.median(np.abs(den[_shell])))

        # [GUARD-2026-09-23] 方向分解 fit 只在带边锚在【真高对称点】上才可靠。
        # 径向（偶次）fit 吸收不了线性梯度：六方胞带边在 K=(1/3,1/3)，若网格 N 不是
        # 3 的倍数则 K 离网，锚落在离网非极值点，「最轻两支」把梯度方向当曲率，
        # 静默给出错误轻质量（CrSe2 46×46 实测 0.52 vs 真值 1.03）。
        _off, _at_hs = _hex_offgrid(recip, kfrac[k0], kfrac)
        if _off:
            return None, _off
        _use_dird = _at_hs

        # [PATCH-2026-09-23] 方向分解 + q^4 外推（依据 tmp/odp/TASK3b_mass_rootcause.md）。
        # 旧做法（下面保留作回退）按半径 0.10->0.50 逐步放宽到"够 3 个点"的 R，再拟合
        # 单一各向同性曲率 C。粗网格下（step3_uniform KSPACING 0.030 -> 15x15x1，
        # 最近邻 0.159 rad/A）会一直吃到 0.276 rad/A 的点（dE 约 0.25 eV，带隙的 27%），
        # 那里能带已明显变平；且把不同方向混进同一个 C，会被重的离轴方向拉高 m*。
        # CrS2/CrSe2 实测偏重 28-37%（论文 0.94/0.96 对 技能 1.285/1.3015）。
        # 新做法：按方向分组，每组用 E = c1 q^2 + c2 q^4 外推到 q->0 取 c1，
        # 最轻的两个方向的几何均值即态密度质量 m_d（DPT 各向同性口径下 m*=m_d）。
        kxy = dk_use[nz]
        _R = np.sqrt(q2n)
        _ang = np.degrees(np.arctan2(kxy[:, 1], kxy[:, 0]))
        # 只保留带边最近几壳层：R > 4*r_nn 的点已在非抛物区，会把 q^4 外推截距拉低。
        _rnz = _R[_R > 1e-6]
        _keep = (_R <= 4.0 * float(_rnz.min())) if len(_rnz) else np.ones(len(_R), bool)
        _Rf, _angf, _denf = _R[_keep], _ang[_keep], den[_keep]
        _grp = {}
        for _i in range(len(_Rf)):
            _grp.setdefault(int(round(_angf[_i] / 5.0)), []).append(_i)
        _dirs, _dirs2 = [], []
        for _idx in _grp.values():
            _idx = np.array(_idx)
            if len(_idx) >= 2:
                _c1 = float(np.polyfit(_Rf[_idx] ** 2, _denf[_idx] / _Rf[_idx] ** 2, 1)[1])
            else:
                _c1 = float(_denf[_idx[0]] / _Rf[_idx[0]] ** 2)
            if abs(_c1) < 1e-9:
                continue
            _rec = (float(3.80998 / abs(_c1)), len(_idx),
                    float(_angf[_idx[0]]), float(_Rf[_idx].max()))
            _dirs.append(_rec)
            if len(_idx) >= 2:
                _dirs2.append(_rec)
        _pool = sorted(_dirs2 if len(_dirs2) >= 2 else _dirs)
        if _use_dird and len(_pool) >= 2:
            _m1, _m2 = _pool[0][0], _pool[1][0]
            _mres = float(np.sqrt(_m1 * _m2))
            return round(_mres, 4), (
                "带边方向分解+q^4外推(%s, %d方向, 最轻两支 %.4f/%.4f, 全组 %.4f..%.4f, spin=%d)" % (
                    "面内" if is_2d else "3D", len(_pool), _m1, _m2,
                    _pool[0][0], _pool[-1][0], isp))
        # 先按半径逐步放宽找邻域；仍不够就退回"最近的若干点"，保证能拟合
        sel, Rused = None, None
        for R in (0.10, 0.15, 0.20, 0.30, 0.40, 0.50):
            m = q2n < R * R
            if np.count_nonzero(m) >= 3:
                sel, Rused = m, R
                break
        if sel is None:
            order = np.argsort(q2n)
            k = min(6, len(q2n))
            sel = np.zeros(len(q2n), bool)
            sel[order[:k]] = True
            Rused = float(np.sqrt(q2n[order[k - 1]]))
        C = np.sum(q2n[sel] * den[sel]) / np.sum(q2n[sel] ** 2)
        if abs(C) < 1e-6:
            return None, "曲率过小/拟合失败"
        m = 3.80998 / abs(C)
        return round(float(m), 4), "能带边抛物拟合(%s, %d点, R<%.2f, spin=%d)" % (
            "面内" if is_2d else "3D", int(np.count_nonzero(sel)), Rused, isp)
    except Exception as e:
        return None, "抛物拟合异常：%s" % _exc_brief(e)


# ---------- E1：从 deformation.h5 带边取（尽力） ----------
def _band_edges_json(cwd):
    """读 step7b_read 生成的 band_edges.json（amset 权威带边定位）。"""
    p = Path(cwd) / DEFORM_READ_DIR / "band_edges.json"
    if not p.is_file():
        return None
    try:
        be = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None
    # [SKILL_REV 交叉校验] band_edges.json 里的 rev 必须与本脚本一致，不一致说明
    # step7b 跑的是另一份 skill 副本（旧副本无 ionrelax 支持，会静默读刚性构型，
    # E1 偏低 20% 且无报错）。沿用 C2 哲学：不一致硬失败，不静默降级。
    be_rev = be.get("skill_rev")
    want = _band_edges_rev()
    if be_rev and want and be_rev != want:
        print("[WARN] band_edges.json 的版本 %s ≠ 当前 S7.1 的 %s —— 它是旧版 S7.1 生成的"
              "（旧副本可能无 ionrelax 支持，E1 偏低且不报错）：重跑 step7b_deform_read" % (be_rev, want))
    return be


def _band_edges_rev():
    """[V153] 当前 S7.1 写 band_edges.json 用的版本（ke_common.BAND_EDGES_REV）；读不到 ke_common -> None（不核对）。"""
    try:
        for _up in Path(__file__).resolve().parents[:2]:
            if (_up / "ke_common.py").is_file() and str(_up) not in sys.path:
                sys.path.insert(0, str(_up))
        import ke_common as _kc
        return getattr(_kc, "BAND_EDGES_REV", None)
    except ImportError:
        return None


def _vac_missing_msg(be, carrier, field):
    """[C2] E1_SOURCE='vac' 但字段缺失时的报错文案：带上生成端的跳过原因。

    step7b_read 会把真空对齐的结局写进 band_edges.json 的 vac_align：
      status=ok      → 字段应该在，缺了说明 json 被改过/口径不一致
      status=skipped → reason 就是根因（3D 无真空 / LOCPOT 缺失 / LVHAR 未开）
    没有 vac_align 说明 json 是旧版本，重跑 step7b_read。"""
    va = (be or {}).get("vac_align") or {}
    st = va.get("status")
    if st == "skipped":
        why = "step7b_read 跳过了真空对齐：%s" % va.get("reason", "未记录原因")
    elif st == "ok":
        why = "step7b_read 声称真空对齐成功，但 %s 缺失（json 不一致，重跑本步上游）" % field
    else:
        why = ("band_edges.json 无 vac_align 段（旧版本 step7b_read）——"
               "重跑 step7b_deform_read 生成新版")
    return ("E1_SOURCE='vac' 但 band_edges.json 无 %s：%s。"
            "看 step7_deform/band_edges.log；确认要用 amset 芯势口径就把 "
            "本材料 step.conf 写 E1_SOURCE = amset（数值口径不同，勿混report）"
            % (field, why))


def get_E1(cwd, carrier):
    """返回 (E1_eV, provenance)。尽力自动；失败返回 (None, 原因)。"""
    key = "E1_electron_eV" if carrier == "electron" else "E1_hole_eV"
    if MANUAL[key]:
        return float(MANUAL[key]), _manual_prov(key)
    # [PATCH-DPT-2026] 优先 amset 权威带边（step7b_read 已定位好，k/band 全部对齐）
    be = _band_edges_json(cwd)
    if _read_dim(cwd) not in (None, "2d"):
        # [0025] 3D：band_edges.json 的 E1_iso_eV 是 |(D_xx+D_yy)/2|（照搬 2D），漏了 D_zz。直接读 deformation.h5
        #   在带边处的 3×3 张量，取对角三项平均（形变势张量的膨胀部分）。h5 读不到才退回旧的面内平均并注明。
        e3, why = _e1_diag_3d(cwd, carrier, be)
        if e3 is not None:
            return round(sum(e3) / 3.0, 4), "deformation.h5 带边 3×3 张量 (D_xx+D_yy+D_zz)/3（芯势口径，%s）" % why
        if be and carrier in be and float(be[carrier].get("E1_iso_eV") or 0) > 0:
            return (round(abs(float(be[carrier]["E1_iso_eV"])), 4),
                    "band_edges.json 的 (D_xx+D_yy)/2（3D 读不到 h5 的 D_zz：%s；非立方体系有偏差）" % why)
    if be and carrier in be and be[carrier].get("n_hits", 0) > 0:
        n = int(be[carrier]["n_hits"])
        # [PATCH-DPT-R6] 按 E1_SOURCE 选口径（vac=真空对齐/amset=平均芯势）。
        # [C2] 显式要了 vac 却没有 → 硬失败，不静默降级到 amset 口径。
        if E1_SOURCE == "vac":
            ev = be[carrier].get("E1_vac_iso_eV")
            if ev is None or abs(float(ev)) == 0:
                # [3D 兜底] 非 2D 体系无真空对齐，vac 口径本就不适用，
                # 自动用 amset 口径（E1_iso_eV，⟨|D|⟩），不属于静默降级。
                if _read_dim(cwd) != "2d":
                    _e1a = float(be[carrier].get("E1_iso_eV") or 0)
                    if _e1a > 0:
                        return round(abs(_e1a), 4), (
                            "band_edges.json 带边(amset 口径,%s 无真空)" % _read_dim(cwd))
                return None, _vac_missing_msg(be, carrier, "E1_vac_iso_eV")
            return round(abs(float(ev)), 4), (
                "band_edges.json 真空对齐(n=%d)面内对角" % n)
        e1 = float(be[carrier]["E1_iso_eV"])
        if e1 > 0:
            return round(abs(e1), 4), (      # [C5] 兜底 abs（amset h5 本身已非负）
                "band_edges.json 带边(简并平均,n=%d)面内对角" % n)
    h5 = Path(cwd) / DEFORM_READ_DIR / "deformation.h5"
    if not h5.is_file():
        return None, "缺 %s/deformation.h5 与 band_edges.json" % DEFORM_READ_DIR
    try:
        import numpy as np
        import h5py
        with h5py.File(str(h5), "r") as f:
            keys = _sorted_dp_keys(f)          # [PATCH-KE-2026] up 在前 down 在后
            if not keys:
                return None, "deformation.h5 无 deformation_potentials"
            probe = f[keys[0]]                 # (nband, nkpt, 3, 3) eV
            bk = _band_edge_index(cwd, carrier, probe.shape[0], probe.shape[1])
            isp = bk[2] if bk is not None else 0
            dset = keys[isp] if isp < len(keys) else keys[0]
            dp = np.array(f[dset])
        if bk is None:
            # [PATCH-DPT-R5] 不再静默走「全带中位数」兜底（那是本次修复的 bug）：
            # 带边定位失败就显式失败，提示修 step7b 或手填 MANUAL_ANISO。
            return None, "带边定位失败：请检查 step7b_read 的 band_edges.json，" \
                         "或在本材料 step.conf 手填 E1_*_EV"
        b, k = bk[0], bk[1]
        e1 = abs((dp[b, k, 0, 0] + dp[b, k, 1, 1]) / 2.0)   # 面内对角均值 |eV|
        return round(float(e1), 4), (
            "deformation.h5 带边(b=%d,k=%d,%s)面内对角(旧逻辑)" % (b, k, dset))
    except Exception as e:
        return None, "读 deformation.h5 异常：%s" % _exc_brief(e)


def _band_edge_index(cwd, carrier, nband, nkpt):
    """用 step3 vasprun 定位 VBM/CBM，返回 (band_idx, kpt_idx, spin_idx)。
    [PATCH-KE-2026] 增加第三个返回值 spin_idx，供选 deformation.h5 的 up/down 数据集。
    失败返回 None。"""
    try:
        import numpy as np
        from pymatgen.io.vasp import Vasprun
        vr = None
        for n in ("vasprun.xml", "vasprun.xml.gz"):
            p = Path(cwd) / UNIFORM_DIR / n
            if p.is_file():
                vr = Vasprun(str(p), parse_dos=False, parse_potcar_file=False)
                break
        if vr is None:
            return None
        hit = _find_band_edge(vr, carrier)          # [PATCH-KE-2026] 跨自旋道
        if hit is None:
            return None
        isp, k, b, _ene = hit
        if b < nband and k < nkpt:
            return int(b), int(k), int(isp)
        return None
    except Exception:
        return None


# ---------- DPT 迁移率公式 ----------
def mobility_dpt(is_2d, C, m_star, E1_eV, T):
    """返回迁移率 cm²/(V·s)。C：2D=N/m，3D=Pa。m_star：m0 倍数。"""
    m = m_star * M0
    E1 = E1_eV * E_C
    if is_2d:
        mu = (E_C * HBAR**3 * C) / (KB * T * m * m * E1**2)     # m²/(V·s)
    else:
        mu = (2.0 * math.sqrt(2 * math.pi) * E_C * HBAR**4 * C) \
             / (3.0 * (KB * T)**1.5 * m**2.5 * E1**2)
    return mu * 1e4                                            # -> cm²/(V·s)


# === patch_dpt_aniso：各向异性 DPT（m*、E1、C 全部分方向）===
# ★ V148：MANUAL_ANISO 没有 step.conf 键（分方向手填少见）；改了它同样作用到所有材料，gen 时 ★ 告警。
MANUAL_ANISO = {          # 分方向手填覆盖，None=自动。x/y 为笛卡尔面内方向
    "m_eff_electron_xy": None,    # (m_x, m_y)，单位 m0，例 (0.98, 0.90)
    "m_eff_hole_xy": None,
    "E1_electron_xy": None,       # (E1_x, E1_y)，eV
    "E1_hole_xy": None,
    "C_2D_xy_N_per_m": None,      # (C_x, C_y)，N/m
}


def get_C_aniso(cwd, is_2d):
    """分方向 2D 刚度 (C_x, C_y) = (C11, C22) × h(V/A)，N/m。"""
    if MANUAL_ANISO["C_2D_xy_N_per_m"]:
        cx, cy = MANUAL_ANISO["C_2D_xy_N_per_m"]
        return (float(cx), float(cy)), "manual"
    c11, c22 = read_elastic_inplane_GPa(cwd)
    if c11 is None or c22 is None:
        return None, "缺 step6 弹性"
    h = _cell_height_A(cwd)
    if not h:
        return None, "缺元胞高度 h=V/A"
    f = 1e9 * (h * 1e-10)
    return (c11 * f, c22 * f), "step6 C11/C22 × h(V/A)=%.3fA" % h


def _time_reversal_ok(v, cwd, carrier):
    """[0022] 时间反演能否当能带对称操作（原在 _expand_full_bz 里，多谷 DPT 归星时共用；逻辑不变）。"""
    use_tr = True
    try:
        ispin = int(v.parameters.get("ISPIN", 1))
        if ispin == 1:
            use_tr = True
        else:
            oc = Path(cwd) / UNIFORM_DIR / "OUTCAR"
            tot_mag = None
            if oc.is_file():
                import re as _re
                m = _re.findall(r"number of electron\s+\S+\s+magnetization\s+(\S+)",
                                oc.read_text(errors="ignore"))
                if m:
                    tot_mag = abs(float(m[-1]))
            if tot_mag is None:
                use_tr = False
                print("[WARN] %s：ISPIN=2 但读不到总磁矩（OUTCAR 缺失/解析失败），"
                      "保守关闭 time_reversal" % carrier)
            else:
                use_tr = (tot_mag <= 0.05)
    except Exception:
        use_tr = True
    return use_tr


def _expand_full_bz(v, cwd, carrier):
    """[0022] IBZ -> 全 BZ（原在 get_effective_mass_aniso 里，多谷 DPT 共用；逻辑不变，只修了 pymatgen 层的
    倒空间旋转）。返回 (kfrac(n_full,3), kp_mapping(n_full,), 方法, None)；失败 (None, None, None, 原因)。"""
    import numpy as np
    kfrac_ibz = np.array(v.actual_kpoints)
    # [P1-2] 展开 IBZ → 全 BZ：带边 K 点坐在 IBZ 楔形边界，近邻点只存在于
    # 楔形内一侧，单侧取点会让 a/b 两个曲率被不对称地吃掉（hex 实测 2.7 倍差）。
    # 展开后 k0 周围点对称，二次型拟合才无偏。失败退回 IBZ（靠 cond/残差验收兜底）。
    # [TR] 时间反演只在净磁矩≈0 时是能带的对称操作。CrS2/CrSe2 塌到零磁矩
    # 无碍，但本 skill 还要复用到 MnIn2Se4/Mn2In2Se5 等磁性体系——自旋分辨
    # 能带在有净磁矩时 TR 会把 up/down 错配。有磁矩就关掉 TR。
    # [TR] 时间反演只在净磁矩≈0 时是能带的对称操作。ISPIN=1 → TR 安全；
    # ISPIN=2 → 读 OUTCAR 总磁矩判断，读不到 → 保守关掉（磁性体系不安全）。
    use_tr = _time_reversal_ok(v, cwd, carrier)

    # [P1-2] 展开 IBZ → 全 BZ，两层：amset（环境）→ pymatgen 点群（不依赖环境）。
    # 【彻底删掉 IBZ 兜底】IBZ 上的 m* 不是"精度差一点"，是系统性错误（单侧
    # 取点，SS 正交胞对称性禁止交叉项却拟合出 off=0.28 的假数）。展开不成功
    # 就硬失败、不出 m*，绝不静默退回 IBZ。
    expanded = False
    _expand_method = None
    try:
        from amset.electronic_structure.symmetry import expand_kpoints
        full_kfrac, _, _, _, _, kp_mapping = expand_kpoints(
            v.final_structure, kfrac_ibz, symprec=0.01,
            time_reversal=use_tr, return_mapping=True)
        kp_mapping = np.asarray(kp_mapping)
        kfrac = np.array(full_kfrac)            # (n_full, 3)
        _expand_method = "amset"
        expanded = True
    except Exception:
        # 层2：subprocess 到 step.conf 的 amset 环境跑 expand_kpoints
        # （jzzn 默认 python 是 atomate2_p_a 无 amset，但 conda 里有 amset_clean）
        _sub_ok = False
        try:
            import subprocess as _sp, shlex as _shlex, json as _json
            _env = None
            try:
                import stepconf as _sc
                _txt = open(_sc.CONF_NAME, encoding="utf-8-sig").read()
                _p = {k.upper(): v for k, v, _ in _sc.parse(_txt, _sc.CONF_NAME).get("params", [])}
                _sh, _e = _p.get("CONDA_SH"), _p.get("AMSET_ENV")
                if _sh and _e:
                    _env = "source %s && conda activate %s" % (_sh, _e)
            except Exception as _ee:                    # noqa: BLE001
                print("[..] %s：读 step.conf 的 amset 环境失败（%s），跳过 subprocess 展开" % (carrier, _exc_brief(_ee)))
            if _env:
                _code = ("import sys,json,numpy as np\n"
                         "from pymatgen.core import Structure\n"
                         "from amset.electronic_structure.symmetry import expand_kpoints\n"
                         "st=Structure.from_dict(json.loads(sys.argv[1]))\n"
                         "kf=np.array(json.loads(sys.argv[2]))\n"
                         "ut=sys.argv[3]=='True'\n"
                         "f,_,_,_,_,kp=expand_kpoints(st,kf,symprec=0.01,time_reversal=ut,return_mapping=True)\n"
                         "print(json.dumps({'full':np.array(f).tolist(),'kp':np.asarray(kp).tolist()}))")
                _cmd = "%s && python3 -c %s %s %s %s" % (
                    _env, _shlex.quote(_code),
                    _shlex.quote(_json.dumps(v.final_structure.as_dict())),
                    _shlex.quote(_json.dumps(kfrac_ibz.tolist())),
                    str(use_tr))
                _out = _sp.run(["bash", "-lc", _cmd], capture_output=True,
                               text=True, timeout=180)
                _lines = [l for l in (_out.stdout or "").strip().splitlines()
                          if l.strip().startswith("{")]
                if _lines:
                    _data = _json.loads(_lines[-1])
                    kp_mapping = np.asarray(_data["kp"])
                    kfrac = np.array(_data["full"])
                    _expand_method = "amset-subproc"
                    expanded = True
                    _sub_ok = True
        except Exception as _es:                        # noqa: BLE001
            print("[..] %s：amset subprocess 展开失败（%s），改用 pymatgen 点群展开" % (carrier, _exc_brief(_es)))
        if not _sub_ok:
            # 层3：pymatgen 点群展开（+ 时间反演，只依赖 pymatgen）
            try:
                from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
                _ops = SpacegroupAnalyzer(
                    v.final_structure, symprec=0.01
                ).get_point_group_operations(cartesian=False)
                # [0022] k 是倒空间分数坐标：实空间分数旋转 R 作用到 k 上是 (R⁻¹)ᵀ，对整个点群取集合
                #   等价于 Rᵀ。以前直接用 op.operate(k)（= R·k），六方胞上展开不满（12×12 网格只得到 97/144 点），
                #   amset 不可用时分方向 m* 一律被拒（AnisoKzTests 在无 amset 的环境里因此失败）。
                _rots = [np.asarray(op.rotation_matrix, float).T for op in _ops]
                _seen, _kf_list, _idx_list = {}, [], []
                for _i, _k in enumerate(kfrac_ibz):
                    _targets = [r @ _k for r in _rots]
                    if use_tr:
                        _targets += [r @ (-_k) for r in _rots]
                    for _kk in _targets:
                        # [0022] 区界点 ±0.5 是同一点：键取 [0,1) 里的值（0.999999 -> 0），存的值落在 [-0.5, 0.5)
                        _w = _kk - np.floor(_kk)
                        _key = tuple(np.round(_w, 6) % 1.0)
                        _kk = _w - (_w >= 0.5 - 1e-9)
                        if _key not in _seen:
                            _seen[_key] = True
                            _kf_list.append(_kk)
                            _idx_list.append(_i)
                if not _kf_list:
                    raise RuntimeError("点群展开得到 0 个点")
                # [网格一致性校验] 与 amset 同等严格：铺满才接受，否则单侧取点
                _df = kfrac_ibz - kfrac_ibz[0]; _df -= np.round(_df)
                _nz = np.abs(_df) > 1e-6
                _mesh = np.array([int(round(1.0 / float(np.min(np.abs(_df[_nz[:, i], i])))))
                                  if _nz[:, i].any() else 1 for i in range(3)])
                _n_expect = int(np.prod(_mesh))
                if len(_kf_list) != _n_expect:
                    raise RuntimeError(
                        "pymatgen 展开只得到 %d/%d 点，网格未铺满"
                        % (len(_kf_list), _n_expect))
                kfrac = np.array(_kf_list)
                kp_mapping = np.asarray(_idx_list)
                _expand_method = "pymatgen"
                expanded = True
            except Exception as _e2:
                return None, None, None, ("全 BZ 展开失败（amset/subprocess/pymatgen 都不可用：%s）"
                                      "——拒绝在 IBZ 上拟合 m*（单侧取点是系统性错误）"
                                      % _exc_brief(_e2))
    return kfrac, kp_mapping, _expand_method, None


def _quad_fit_2d(kfrac, recip, eband, k0, need_lin=None, mesh=None, kops=None):
    """[0022] 带边二次型拟合的取点 + 最小二乘 + 验收（原在 get_effective_mass_aniso 里，多谷 DPT 逐谷共用；逻辑不变）。
    kfrac：全 BZ 分数坐标；eband：同一组点上该带的能量；k0：拟合中心的下标。
    need_lin=None 时按 k0 是否高对称点自动决定；mesh=None 时由 kfrac 反解。
    [0027] kops（k 空间对称操作，_mv_kops）给出时按 k0 的小群判断（_k_pinned），不再查 ⅓/¼ 坐标表。
    返回 ({a, b, c, den=4|a||b|-c², npt, n_shell, Rused, names, cond, rel}, None)；失败 (None, 原因)。"""
    import numpy as np
    kcart = kfrac @ recip
    dk = kcart - kcart[k0]
    de = eband - eband[k0]
    zm = np.abs(dk[:, 2]) < 0.02                 # 只取面内
    dkx, dky, dee = dk[zm][:, 0], dk[zm][:, 1], de[zm]
    q2 = dkx ** 2 + dky ** 2
    qnz = np.sqrt(q2[q2 > 1e-9])
    if len(qnz) == 0:
        return None, "带边面内无有效 k 点（q > 0）"
    q_min = float(qnz.min())
    # [SHELL+AXIS] 网格分割数：从全 BZ k 点反解。半径由「最差采样轴」的笛卡尔
    # 步长定，细长胞（LS a=21.64, N_x≈2）不会圈出单方向长条。壳层用相对距离
    # q/q_min 聚类（容差 1e-3），避免浮点噪声把单壳层劈成两个。
    if mesh is None:
        try:
            from amset.electronic_structure.symmetry import get_mesh_from_kpoint_diff
            _mesh, _ = get_mesh_from_kpoint_diff(kfrac)
            mesh = np.rint(np.asarray(_mesh)).astype(int)
        except Exception:
            _df = kfrac - kfrac[k0]; _df -= np.round(_df)
            _nz = np.abs(_df) > 1e-6
            mesh = np.array([int(round(1.0 / float(np.min(np.abs(_df[_nz[:, i], i])))))
                             if _nz[:, i].any() else 1 for i in range(3)])
    step_cart = np.array([np.linalg.norm(recip[i]) / max(int(mesh[i]), 1)
                          for i in (0, 1)])
    # [SS/LS 线性项] 带边非高对称点（Y-Γ 路径，如 SS k=(0,0.3636)）时，真实
    # 极值点一般不在网格点上，存在真实线性项——继续升 NSTEP 直到点数够 3 倍
    # 冗余（18 点），否则线性项被关掉，k0 偏移被 a/b 吸收（重演 CrS2 曲率偏置，
    # 且这次没有对称性兜底）。
    def _is_high_sym(kf, tol=1e-3):
        _bases = (0.0, 0.5, 1.0 / 3.0, 2.0 / 3.0, 0.25, 0.75)
        for _c in kf[:2]:
            _c = _c % 1.0
            if not any(abs(_c - _b) < tol or abs(_c - _b - 1.0) < tol
                       for _b in _bases):
                return False
        return True
    if need_lin is None:
        # [0027] 坐标表把矩形胞 Y–Γ 线上的 (0, ⅓) 也当成高对称点（⅓ 是给六方 K 用的）：SS 导带底落在网格点 ⅓、
        #   真实极小在 ≈0.36，NSTEP=2 只取到 <18 个点，线性项被关掉，rel=0.23 拒绝。有对称操作时改按小群判断。
        need_lin = (not _k_pinned(kfrac[k0], kops, 2)) if kops else (not _is_high_sym(kfrac[k0]))

    sel = Rused = n_shell = NSTEP_used = None
    _nstep_iter = (FORCE_NSTEP,) if FORCE_NSTEP else (2, 3, 4, 5)
    for NSTEP in _nstep_iter:
        # -0.1 留余量：整数步恰好卡在下一壳层边界，浮点噪声会把下一壳层的
        # 个别点纳入。1.9 步≈0.302 落在第二/三壳层之间。
        Rmax = (NSTEP - 0.1) * float(step_cart.min())
        m = (q2 < Rmax * Rmax) & (q2 > 1e-9)
        n = int(np.count_nonzero(m))
        if n < 8:
            continue
        q_sel = np.sqrt(q2[m])
        shells = np.unique(np.round(q_sel / q_min, 3))
        n_sh = int(len(shells))
        ok = (n_sh >= 2)                        # 至少 2 壳层，才有径向信息
        if need_lin:
            ok = ok and (n >= 18)               # 非高对称点：线性项需 3 倍冗余
        if ok:
            sel, Rused, n_shell, NSTEP_used = m, Rmax, n_sh, NSTEP
            break
    if sel is None:
        return None, ("带边面内点太少或全在单壳层（需 ≥8 点 + ≥2 壳层）——"
                      "加密 step3_uniform 网格或手填 MANUAL_ANISO")
    npt = int(np.count_nonzero(sel))

    # [AXIS 覆盖] 每个面内倒格矢方向都要有非零步的点，否则该方向曲率无约束。
    # LS 的 kx 轴只有 N≈2 分割，会在这里干脆报错而不是返回假 m*_x。
    # [V151] sel 是在面内子集（zm）上算的掩码，必须先取 kfrac[zm]：2D 的 S3 常用 kz≥3（WS2 48×48×3），
    #   全 BZ 6912 点 vs 面内 2304 点 -> IndexError，分方向 m* 整块失败（kz=1 时两者等长，测不出来）。
    _dfrac = kfrac[zm][sel] - kfrac[k0]
    _dfrac -= np.round(_dfrac)
    _steps = np.rint(_dfrac * mesh).astype(int)
    for _i in (0, 1):
        if int(np.abs(_steps[:, _i]).max()) < 1:
            return None, ("倒格矢方向 %d 无非零步取点（该轴分割数 N=%d 过少），"
                          "该方向曲率完全无约束——把 step3_uniform 的 KSPACING "
                          "调小或设 KMIN_DIV≥12" % (_i, int(mesh[_i])))
    # [P1] 设计矩阵：高对称点（K）在网格上时线性项被 C₃ 强制为 0，但 SS/LS
    # 的带边在 Y–Γ 路径上（非高对称点），真实极值点一般不在网格点上，存在
    # 真实线性项——点数够就拟合它，吸收 k0 偏移；否则只拟合二次项。
    # [Q4] ≥2 壳层时再加各向同性四次项 q⁴，把次抛物偏置直接拟掉。
    use_quartic = (n_shell >= 2) and (npt >= 8)
    # 线性项（d,e,f）需 3 倍冗余才稳定：2 倍冗余（14 点）在对称点上会把曲率
    # 拟合进线性项（hex 实测 m* 1.124 vs 0.924）。≥18 点才启用。
    use_linear = (npt >= 18)
    cols = [dkx[sel] ** 2, dky[sel] ** 2, dkx[sel] * dky[sel]]
    names = ["a", "b", "c"]
    if use_quartic:
        cols.append(q2[sel] ** 2)
        names.append("q4")
    if use_linear:
        cols += [dkx[sel], dky[sel], np.ones(npt)]
        names += ["d", "e", "f"]
    A = np.column_stack(cols)
    if A.shape[0] < 2 * A.shape[1]:
        return None, ("拟合点数不足：%d 点 / %d 参数（需 ≥2 倍冗余）——"
                      "加密 step3_uniform 网格" % (A.shape[0], A.shape[1]))
    cond = float(np.linalg.cond(A))
    if cond > 1e5:
        return None, ("设计矩阵病态 cond=%.1e（k 点单侧/共线，"
                      "IBZ 楔形取点不对称）" % cond)
    sol, *_ = np.linalg.lstsq(A, dee[sel], rcond=None)
    a, b, c = float(sol[0]), float(sol[1]), float(sol[2])
    pred = A @ sol
    rel = float(np.linalg.norm(dee[sel] - pred)
                / (np.linalg.norm(dee[sel]) + 1e-12))
    if rel > 0.2:
        return None, "二次型拟合残差过大 rel=%.2f（带边附近非二次型）" % rel
    # [P2] 逆质量张量求逆的对角元（H=[[2a,c],[c,2b]]）。
    # 质量恒正：electron(CBM) 曲率 a,b>0，hole(VBM) 曲率 a,b<0——统一取 |a|,|b|，
    # 否则 hole 会得到负质量（den=4ab-c² 因 ab>0 为正，但分子 4b/4a 带负号）。
    if a * b <= 0:
        return None, ("带边曲率异号（a=%.3f b=%.3f，鞍点？"
                      "带边定位可疑或取点跨了多个带边）" % (a, b))
    aa, bb = abs(a), abs(b)
    den = 4.0 * aa * bb - c * c
    if den <= 1e-12:
        return None, "逆质量张量奇异（4|a||b|-c²=%.2g，能量椭圆退化）" % den
    return {"a": a, "b": b, "c": c, "den": den, "npt": npt, "n_shell": n_shell,
            "Rused": Rused, "names": names, "cond": cond, "rel": rel,
            "coef": {n: float(x) for n, x in zip(names, sol)},
            "step_min": float(step_cart.min())}, None


def get_effective_mass_aniso(cwd, carrier, is_2d):
    """带边二次型拟合 ΔE = a·Δkx² + b·Δky² + c·Δkx·Δky + d·Δkx + e·Δky + f。

    [P2] 分方向质量 = 逆质量张量 H⁻¹ 的对角元（H=∂²E/∂k∂k=[[2a,c],[c,2b]]）：
         m_x/m0 = 3.80998·4b/(4ab-c²)，m_y/m0 = 3.80998·4a/(4ab-c²)。
         c=0 退化回 m_x=3.80998/|a|、m_y=3.80998/|b|（注意 m_x 由 b 决定）。
    [P1] 线性项 d/e/f 吸收 k0 不在真实极值点的偏移；加 cond(A)+相对残差验收，
         病态/残差大就报错，绝不静默返回一个错的 m*。
    返回 ((m_x, m_y), provenance)；失败 (None, 原因)。"""
    key = "m_eff_electron_xy" if carrier == "electron" else "m_eff_hole_xy"
    if MANUAL_ANISO[key]:
        mx, my = MANUAL_ANISO[key]
        return (float(mx), float(my)), "manual"
    if not is_2d:
        return None, "各向异性 m* 目前只实现 2D 面内"
    vr = None
    for n in ("vasprun.xml", "vasprun.xml.gz"):
        p = Path(cwd) / UNIFORM_DIR / n
        if p.is_file():
            vr = p
            break
    if vr is None:
        return None, "缺 step3_uniform/vasprun.xml"
    try:
        import numpy as np
        from pymatgen.io.vasp import Vasprun
        v = Vasprun(str(vr), parse_dos=False, parse_potcar_file=False)
        recip = v.final_structure.lattice.reciprocal_lattice.matrix
        # [CELL 取向] 正交胞（90°）且长轴在 b 时，by_direction 的 x/y 相对文献
        # 转置（数值全对结论全反）。SS/LS 需 a=长轴(调制方向)在 x。
        _lat = v.final_structure.lattice
        if (abs(_lat.alpha - 90) < 1 and abs(_lat.beta - 90) < 1
                and abs(_lat.gamma - 90) < 1 and _lat.b > _lat.a * 1.2):
            print("[WARN] 胞 a=%.2f < b=%.2f：面内长轴不在 x 上——若本体系有"
                  "面内各向异性，by_direction 的 x/y 可能相对文献转置"
                  "（查 step.conf 的 CELL_POLICY）" % (_lat.a, _lat.b))
        kfrac_ibz = np.array(v.actual_kpoints)      # IBZ (nk,3)
        hit = _find_band_edge(v, carrier)        # [PATCH-KE-2026] 跨自旋道
        if hit is None:
            return None, "无未占据/占据态 —— 体系是金属？（带边不存在，m* 无意义）"
        isp, k0_ibz, b0, ene_ibz = hit           # ene_ibz (nk, nb)

        # [交叉校验] m* 链（本步 step3_uniform）与 E1 链（step7 undeformed，
        # band_edges.json 里的 k_frac）必须定位到同一个带边 k 点；不一致说明
        # 两链各自独立定位时在某个体系上分了叉，m* 与 E1 内部不自洽（照样不报错）。
        _be = _band_edges_json(cwd)
        if _be and carrier in _be:
            _hits = _be[carrier].get("hits") or []
            if _hits and _hits[0].get("k_frac") is not None:
                _kf_be = np.array(_hits[0]["k_frac"], dtype=float)
                _kf_m = kfrac_ibz[k0_ibz].astype(float)
                _d = _kf_be - _kf_m
                _d -= np.round(_d)                   # 周期回卷
                # [V167] 本函数只走 2D：c 是真空方向，能带沿 kz 是平的 —— [⅓,⅓,0] 与 [⅓,⅓,⅓] 是同一个 K 谷
                #   （CrSe2：E1 链在 kz=⅓ 层定位、m* 链在 kz=0 层，以前误报"两链带边分叉"）。只比面内分量。
                _d[2] = 0.0
                if float(np.linalg.norm(_d)) > 1e-4:
                    print("[WARN] %s：m* 链带边 k=%s 与 E1 链带边 k=%s 不一致"
                          "（回卷距离 %.3g）——两链带边分叉，m* 与 E1 内部不自洽"
                          % (carrier, np.round(_kf_m, 4).tolist(),
                             np.round(_kf_be, 4).tolist(), float(np.linalg.norm(_d))))

        kfrac, kp_mapping, _expand_method, _xerr = _expand_full_bz(v, cwd, carrier)
        if _xerr:
            return None, _xerr
        # 展开后（amset/pymatgen 共用）：能量映射 + 精确匹配 IBZ 原 k
        ene = ene_ibz[kp_mapping]                   # (n_full, nb)
        d = kfrac - kfrac_ibz[k0_ibz]; d -= np.round(d)
        dd = np.linalg.norm(d, axis=1)
        k0 = int(np.argmin(dd))
        if dd[k0] > 1e-4:
            return None, "展开后找不到原带边 k 点（回卷距离 %.3g）" % dd[k0]
        _off, _ = _hex_offgrid(recip, kfrac[k0], kfrac)                 # [V152] 与各向同性路径同一道闸
        if _off:
            return None, _off

        try:                                                             # [0027]
            _kops = _mv_kops(v.final_structure, _time_reversal_ok(v, cwd, carrier))
        except Exception:                                                # noqa: BLE001
            _kops = None
        fit, _ferr = _quad_fit_2d(kfrac, recip, ene[:, b0], k0, kops=_kops)
        if _ferr:
            return None, _ferr
        a, b, c, den = fit["a"], fit["b"], fit["c"], fit["den"]
        aa, bb = abs(a), abs(b)
        npt, n_shell, Rused, names, cond, rel = (fit[k] for k in ("npt", "n_shell", "Rused", "names", "cond", "rel"))
        mx = 3.80998 * (4.0 * bb) / den
        my = 3.80998 * (4.0 * aa) / den
        if mx <= 0 or my <= 0:
            return None, "质量为负（mx=%.3f my=%.3f）" % (mx, my)
        off = abs(c) / max(abs(a), abs(b))
        prov = ("带边二次型拟合(%s, %d点/%d壳层, R<%.2f, 模型=%s, "
                "交叉项/主项=%.2f, cond=%.0e, rel=%.3f, spin=%d)"
                % (_expand_method, npt, n_shell, Rused,
                   "+".join(names), off, cond, rel, isp))
        # [P2 后] 交叉项判据：正交胞（SS/LS/*_ortho）的点群禁止面内交叉项，
        # off 明显非零 = 拟合本身有问题（取点不对称/带边非极值点/跨带），
        # 不再是"手填 MANUAL_ANISO 就能绕过"的表象问题。
        if off > 0.3:
            print("[WARN] %s：m* 拟合交叉项/主项=%.2f 偏大 —— 六方胞里 C₃ 强制"
                  "各向同性、正交胞里对称性禁止面内交叉项，两种情况下 off 都应"
                  "接近 0；非零说明拟合有问题（取点不对称/带边不是真极值点/"
                  "取点跨了多条带），此时分方向 m* 与由它导出的 μ、τ 都不可信。"
                  % (carrier, off))
        return (round(mx, 4), round(my, 4)), prov
    except Exception as e:
        return None, "二次型拟合异常：%s" % _exc_brief(e)


def get_E1_aniso(cwd, carrier):
    """带边形变势 (|D_xx|, |D_yy|)，eV —— 不再取面内平均。"""
    key = "E1_electron_xy" if carrier == "electron" else "E1_hole_xy"
    if MANUAL_ANISO[key]:
        ex, ey = MANUAL_ANISO[key]
        return (float(ex), float(ey)), "manual"
    # [PATCH-DPT-2026] 优先 amset 权威带边（step7b_read 已定位好，k/band 全部对齐）
    be = _band_edges_json(cwd)
    if be and carrier in be and be[carrier].get("n_hits", 0) > 0:
        n = int(be[carrier]["n_hits"])
        # [PATCH-DPT-R6] 按 E1_SOURCE 选口径。[C5] 两套字段生成端都已取 |·|，
        # 这里的 abs() 只是兜底。[C2] 要了 vac 却没有 → 硬失败，不降级。
        if E1_SOURCE == "vac":
            vx = be[carrier].get("E1_vac_xx_eV")
            vy = be[carrier].get("E1_vac_yy_eV")
            if (vx is None or vy is None
                    or abs(float(vx)) == 0 or abs(float(vy)) == 0):
                return None, _vac_missing_msg(be, carrier,
                                              "E1_vac_xx_eV/E1_vac_yy_eV")
            return ((round(abs(float(vx)), 4), round(abs(float(vy)), 4)),
                    "band_edges.json 真空对齐(n=%d) D_xx/D_yy" % n)
        ex = float(be[carrier]["E1_xx_eV"])
        ey = float(be[carrier]["E1_yy_eV"])
        if ex > 0 and ey > 0:
            return ((round(abs(ex), 4), round(abs(ey), 4)),   # [C5] 兜底 abs
                    "band_edges.json 带边(简并平均,n=%d) D_xx/D_yy" % n)
    h5 = Path(cwd) / DEFORM_READ_DIR / "deformation.h5"
    if not h5.is_file():
        return None, "缺 %s/deformation.h5 与 band_edges.json" % DEFORM_READ_DIR
    try:
        import numpy as np
        import h5py
        with h5py.File(str(h5), "r") as f:
            keys = _sorted_dp_keys(f)          # [PATCH-KE-2026] up 在前 down 在后
            if not keys:
                return None, "deformation.h5 无 deformation_potentials"
            probe = f[keys[0]]
            bk = _band_edge_index(cwd, carrier, probe.shape[0], probe.shape[1])
            isp = bk[2] if bk is not None else 0
            dset = keys[isp] if isp < len(keys) else keys[0]
            dp = np.array(f[dset])
        if bk is None:
            return None, "带边定位失败：请检查 step7b_read 的 band_edges.json，" \
                         "或手填 MANUAL_ANISO 的 E1"
        b, k = bk[0], bk[1]
        return ((round(float(abs(dp[b, k, 0, 0])), 4),
                 round(float(abs(dp[b, k, 1, 1])), 4)),
                "deformation.h5 带边(b=%d,k=%d,%s) D_xx/D_yy(旧逻辑)" % (b, k, dset))
    except Exception as e:
        return None, "读 deformation.h5 异常：%s" % _exc_brief(e)


def mobility_dpt_2d_aniso(C_alpha, m_alpha, m_d, E1_eV, T):
    """μ_α = eℏ³C_2D,α/(k_B T m*_α m_d E1_α²)，cm²/(V·s)。m_d=√(m_x m_y)。"""
    E1 = E1_eV * E_C
    mu = (E_C * HBAR ** 3 * C_alpha) / (KB * T * (m_alpha * M0) * (m_d * M0) * E1 ** 2)
    return mu * 1e4


def _aniso_block(cwd, is_2d, carrier, T):
    """分方向 DPT 结果块。[0025] 3D 走 _aniso_block_3d（x/y/z）。"""
    if not is_2d:
        return _aniso_block_3d(cwd, carrier, T)
    Cxy, Cprov = get_C_aniso(cwd, is_2d)
    mxy, mprov = get_effective_mass_aniso(cwd, carrier, is_2d)
    exy, eprov = get_E1_aniso(cwd, carrier)
    blk = {"C_2D_N_per_m": Cxy, "C_provenance": Cprov,
           "m_eff_m0": mxy, "m_provenance": mprov,
           "E1_eV": exy, "E1_provenance": eprov, "T_K": T,
           "formula": "mu_a = e*hbar^3*C_2D,a/(kB*T*m_a*m_d*E1_a^2); tau_a = mu_a*m_a/e"}
    miss = [n for n, v in (("C", Cxy), ("m*", mxy), ("E1", exy)) if v is None]
    if miss:
        blk["status"] = "缺输入：%s（可在 MANUAL_ANISO 手填）" % "、".join(miss)
        return blk
    m_d = math.sqrt(mxy[0] * mxy[1])
    blk["m_d_m0"] = round(m_d, 4)
    for i, d in enumerate(("x", "y")):
        mu = mobility_dpt_2d_aniso(Cxy[i], mxy[i], m_d, exy[i], T)
        blk[d] = {"C_2D_N_per_m": round(Cxy[i], 3), "m_eff_m0": mxy[i],
                  "E1_eV": exy[i], "mobility_cm2_Vs": round(mu, 3),
                  "tau_s": (mu * 1e-4) * (mxy[i] * M0) / E_C}
    blk["status"] = "ok"
    return blk


# ============================ [0022] 多谷 DPT ============================
# 单谷 DPT（上面的各向同性值与 by_direction）只看全局带边那一个谷：带边在 Q/Λ 这类非高对称点（WSe2 导带）时
# 单谷 m* 拒绝出数；带边附近还有别的谷（MoS2 价带 K 与 Γ、导带 K 与 Q）时单谷值也漏了它们。这里全自动：
#   ① 在 S3 全 BZ 网格的 kz=0 层上找离带边 VALLEY_WINDOW_EV 以内的所有谷（各带、各自旋道的局部极值），
#      按点群（+时间反演）归成星：K 星 2 个谷、Q 星 6 个、Γ 1 个、M 3 个；
#   ② 每个谷单独做二次型拟合（与 by_direction 同一个 _quad_fit_2d）-> 逆质量张量 W=M⁻¹（1/m0）、
#      态密度质量 m_d=√det M、拟合出的极值能量（离网的 Q 谷不取网格点的能量）；
#   ③ 每个谷单独取形变势，与 S7.1 的 E1_vac 同一定义：S7 形变构型（S7.1 实际读的那份，ionrelax 优先）里
#      该谷 k 点、该带本征值的中心差分减去真空势的应变导数，E1_α = |dE/dε_αα − dE_vac/dε_αα|；
#      不用 deformation_vac.h5：AMSET 把对称化转出来的斜向应变也平均进 D_xx/D_yy，各向异性的 Q 谷会被抹平；
#   ④ τ_α,v = ħ³C_α/(k_BT·m_d,v·E1_α,v²)；非简并统计下谷 v 的载流子数 n_v ∝ m_d,v·exp(−ΔE_v/k_BT)：
#        μ_αα(T) = (e/m0)·Σ_v n_v·τ_α,v·W_αα,v / Σ_v n_v       每个温度单独算（谷间布居随 T 变，不再是 μ∝1/T）
#      只有一个星且交叉项为零时（K 谷材料）与 by_direction 完全相同。不含谷间散射 -> 偏上限。
HB2_2M0 = 3.80998          # ħ²/(2m0)，eV·Å²
KINK_TOL_EV = 0.02         # 按能量排序的两条带相交处，上面那条带是 V 形"极小"（不是谷）：靠带边一侧的相邻带
                           #   在这点差 < 20 meV 且它自己在这点不是极值 -> 当交叉点跳过（记进 skipped）
_MV_CACHE = {}             # 本进程内两种载流子共用：S3 vasprun 与 S7 形变构型只解析一次


class _MVStop(Exception):
    """多谷块的"拒绝出数 + 原因"（不是程序错误）。"""


def _vasprun_s3(cwd):
    for n in ("vasprun.xml", "vasprun.xml.gz"):
        p = Path(cwd) / UNIFORM_DIR / n
        if p.is_file():
            return p
    return None


def _cached_vasprun(path):
    st = Path(path).stat()
    key = ("vr", str(Path(path).resolve()), st.st_mtime_ns, st.st_size)
    if key not in _MV_CACHE:
        from pymatgen.io.vasp import Vasprun
        _MV_CACHE[key] = Vasprun(str(path), parse_dos=False, parse_potcar_file=False)
    return _MV_CACHE[key]


def _cached_eig(path):
    """(k 列表, {自旋名: 本征值数组})，按文件缓存。"""
    import numpy as np
    st = Path(path).stat()
    key = ("eig", str(Path(path).resolve()), st.st_mtime_ns, st.st_size)
    if key not in _MV_CACHE:
        v = _cached_vasprun(path)
        _MV_CACHE[key] = (np.asarray(v.actual_kpoints, float), dict(_spin_items(v)))
    return _MV_CACHE[key]


def _spin_items(v):
    """[(自旋名 up/down, arr(nk,nb,2)), ...]，顺序与 _spin_arrays 相同。"""
    import numpy as np
    eigmap = v.eigenvalues
    keys = list(eigmap.keys())
    try:
        from pymatgen.electronic_structure.core import Spin
        keys = [s for s in (Spin.up, Spin.down) if s in eigmap] or keys
    except Exception:
        pass
    return [(str(getattr(s, "name", s)).lower(), np.asarray(eigmap[s])) for s in keys]


def _kidx_tr(kp, kf, tol=1e-4):
    """kf（或 −kf：VASP 的 ISYM=0 仍按时间反演只留半个网格）在 k 列表里的下标；没有 -> None。"""
    import numpy as np
    for sg in (1.0, -1.0):
        d = kp - sg * np.asarray(kf, float)
        d -= np.round(d)
        dd = np.linalg.norm(d, axis=1)
        i = int(np.argmin(dd))
        if dd[i] <= tol:
            return i
    return None


def _mv_plane(kfrac):
    """全 BZ 的 kz=0 层 -> {"N": (N1, N2), "idx": (N1,N2) 网格坐标 -> kfrac 下标}；失败 (None, 原因)。"""
    import numpy as np
    w = np.asarray(kfrac, float)
    w = w - np.round(w)
    sel = np.where(np.abs(w[:, 2]) < 1e-6)[0]
    if len(sel) < 9:
        return None, "全 BZ 的 kz=0 层只有 %d 个点" % len(sel)
    p = w[sel, :2]
    N = []
    for a in (0, 1):
        nz = np.abs(p[:, a]) > 1e-6
        if not nz.any():
            return None, "面内方向 %d 只有 1 个分割" % a
        N.append(int(round(1.0 / float(np.min(np.abs(p[nz, a]))))))
    N = np.array(N)
    g = p * N
    if float(np.max(np.abs(g - np.round(g)))) > 1e-3:
        return None, "step3_uniform 的面内网格不是 Γ 中心的均匀网格"
    ij = np.mod(np.round(g).astype(int), N)
    key = ij[:, 0] * N[1] + ij[:, 1]
    if len(np.unique(key)) != len(key) or len(key) != int(N[0] * N[1]):
        return None, "kz=0 层 %d 点 ≠ %d×%d（全 BZ 展开不完整）" % (len(key), N[0], N[1])
    idx = np.full((int(N[0]), int(N[1])), -1, int)
    idx[ij[:, 0], ij[:, 1]] = sel
    return {"N": (int(N[0]), int(N[1])), "idx": idx}, None


def _mv_extrema(E, sgn):
    """周期网格 E(N1,N2) 上的局部极值（sgn=+1 极小 / −1 极大，8 邻域）。极值正好落在两个网格点中间时
    （如奇数网格的 M）两点等高、都算极值：8 邻接且等高的并成一个谷。返回 [[(i,j), ...], ...]。"""
    import numpy as np
    s = sgn * np.asarray(E, float)
    ok = np.ones(s.shape, bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                ok &= s <= np.roll(s, (-di, -dj), axis=(0, 1)) + 1e-6
    pts = [tuple(int(x) for x in q) for q in np.argwhere(ok)]
    n1, n2 = s.shape
    par = {q: q for q in pts}

    def find(q):
        while par[q] != q:
            par[q] = par[par[q]]
            q = par[q]
        return q
    for (i, j) in pts:
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                q = ((i + di) % n1, (j + dj) % n2)
                if q != (i, j) and q in par and abs(s[i, j] - s[q]) < 1e-3:
                    par[find((i, j))] = find(q)
    groups = {}
    for q in pts:
        groups.setdefault(find(q), []).append(q)
    return [sorted(g) for g in groups.values()]


def _mv_collect(ene, idx, N, nocc, sgn, edge, win):
    """各自旋道、各带（从带边往外）在 kz=0 网格上的极值。ene：{自旋: (n_full, nb)}；idx：网格坐标 -> 下标；
    edge：sgn·E 的带边值。返回 (valleys=[(自旋, 带, [(i,j)...], 网格能量)], seeds=窗口内的下标, skipped=交叉点说明)。
    窗口外的极值也收（归星时找对称像用），带的最低点（sgn 意义下）出窗口就不再往外扫。"""
    import numpy as np
    valleys, seeds, skipped = [], set(), []
    for sp in ene:
        nbands = ene[sp].shape[1]
        for b in (range(nocc[sp], nbands) if sgn > 0 else range(nocc[sp] - 1, -1, -1)):
            E2 = ene[sp][idx, b]
            if float((sgn * E2).min()) > edge + win:
                break
            nb = b - 1 if sgn > 0 else b + 1                  # 靠带边一侧的相邻带（也得是同类带）
            has_nb = (nb >= nocc[sp]) if sgn > 0 else (nb < nocc[sp])
            nb_ext = {q for g in _mv_extrema(ene[sp][idx, nb], sgn) for q in g} if has_nb else set()
            for pts in _mv_extrema(E2, sgn):
                e_here = float(E2[pts[0]])
                if has_nb and abs(e_here - float(ene[sp][idx, nb][pts[0]])) < KINK_TOL_EV \
                        and not (set(pts) & nb_ext):
                    if sgn * e_here <= edge + win:
                        skipped.append("%s 带 %d k=%s：与带 %d 交叉（差 %.0f meV），不是谷" % (
                            sp, b, (np.array(pts[0]) / np.array(N)).round(4).tolist(), nb,
                            1000 * abs(e_here - float(ene[sp][idx, nb][pts[0]]))))
                    continue
                if sgn * e_here <= edge + win:
                    seeds.add(len(valleys))
                valleys.append((sp, b, pts, e_here))
    return valleys, seeds, skipped


def _mv_kops(structure, use_tr):
    """点群在倒空间分数坐标上的作用（Rᵀ，整数矩阵），有时间反演再加 −Rᵀ；只留保 z 轴的（面内）。"""
    import numpy as np
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    ops = SpacegroupAnalyzer(structure, symprec=0.01).get_point_group_operations(cartesian=False)
    out = []
    for op in ops:
        r = np.rint(np.asarray(op.rotation_matrix, float).T).astype(int)
        if np.any(r[2, :2]) or np.any(r[:2, 2]):
            continue
        for rr in ((r, -r) if use_tr else (r,)):
            if not any(np.array_equal(rr, u) for u in out):
                out.append(rr)
    return out


def _kops_3d(structure, use_tr):
    """[0027] 点群在倒空间分数坐标上的作用（Rᵀ，整数矩阵），有时间反演再加 −Rᵀ；3D 用，不筛面内。"""
    import numpy as np
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    out = []
    for op in SpacegroupAnalyzer(structure, symprec=0.01).get_point_group_operations(cartesian=False):
        r = np.rint(np.asarray(op.rotation_matrix, float).T).astype(int)
        for rr in ((r, -r) if use_tr else (r,)):
            if not any(np.array_equal(rr, u) for u in out):
                out.append(rr)
    return out


def _k_pinned(kf, kops, nd):
    """[0027] 对称性是否把 k 处的能量梯度钉成零：取 k 的小群（kops 里把 k 映回自身、模倒格矢的操作），
    梯度只能落在小群矩阵的公共不动子空间里，而小群矩阵的平均就是到这个子空间的投影 —— 平均为零矩阵才算钉住。
    六方 K（C₃）、Γ / M / X / Y（含 −1）是；矩形胞 Y–Γ 线上的 (0, ⅓)（小群只有 E 与 m_x）不是。nd：2 只看面内。"""
    import numpy as np
    kf = np.asarray(kf, float)
    little = []
    for r in kops:
        d = (np.asarray(r) @ kf - kf)[:nd]
        if float(np.max(np.abs(d - np.round(d)))) < 1e-6:
            little.append(np.asarray(r, float)[:nd, :nd])
    return bool(little) and float(np.max(np.abs(sum(little) / len(little)))) < 1e-6


def _mv_rot3(r):
    """面内部分是 ≥3 次转轴（C3/C4/C6）。"""
    import numpy as np
    r2 = np.asarray(r)[:2, :2]
    return int(round(float(np.linalg.det(r2)))) == 1 and abs(int(np.trace(r2))) <= 1


def _mv_site_iso(kf, kops):
    """k 点的小群（含时间反演）有 ≥3 次转轴（K、Γ）-> 该谷的面内形变势张量各向同性，E1_x = E1_y。"""
    import numpy as np
    kf = np.asarray(kf, float)
    for r in kops:
        d = r @ kf - kf
        d -= np.round(d)
        if float(np.max(np.abs(d[:2]))) < 1e-6 and _mv_rot3(r):
            return True
    return False


def _mv_stars(valleys, N, kops, seeds=None):
    """valleys：[(自旋, 带, [(i,j),...], E)]（扫过的带上的全部极值）。从 seeds（窗口内的谷）出发按点群轨道归星
    -> ([[谷下标...], ...], None)；失败 (None, 原因)。"""
    import numpy as np
    N = np.asarray(N)
    where = {}
    for vi, (sp, b, pts, _e) in enumerate(valleys):
        for q in pts:
            where[(sp, b, q)] = vi
    star_of = [-1] * len(valleys)
    stars = []
    for vi, (sp, b, pts, _e) in enumerate(valleys):
        if star_of[vi] >= 0 or (seeds is not None and vi not in seeds):
            continue
        k = np.array([pts[0][0] / N[0], pts[0][1] / N[1], 0.0])
        mem = set()
        for r in kops:
            g = (r @ k)[:2] * N
            if float(np.max(np.abs(g - np.round(g)))) > 1e-6:
                return None, "step3_uniform 的面内网格 %d×%d 在晶体对称操作下不封闭" % (N[0], N[1])
            q = tuple(int(x) for x in np.mod(np.round(g).astype(int), N))
            vj = where.get((sp, b, q))
            if vj is None:
                return None, ("谷 k=%s（%s 带 %d）的对称像 %s 不是网格上的极值点：网格与晶体对称性不一致"
                              % (np.round(k[:2], 4).tolist(), sp, b, (np.array(q) / N).round(4).tolist()))
            mem.add(vj)
        for vj in mem:
            if star_of[vj] >= 0:
                return None, "谷的对称轨道互相重叠（对称操作不成群？）"
            star_of[vj] = len(stars)
        stars.append(sorted(mem))
    return stars, None


def _mv_label(kf, recip, hexcell):
    """星的名字（只作标注）：六方 Γ/K/M/Q(Γ–K 线)/Γ–M；矩形 Γ/X/Y/S；其余写分数坐标。"""
    import numpy as np
    k = np.asarray(kf, float)[:2]
    k = k - np.round(k)
    if np.all(np.abs(k) < 1e-6):
        return "Γ"
    if hexcell:
        kc = np.array([k[0], k[1], 0.0]) @ np.asarray(recip)
        b = float(np.linalg.norm(recip[0]))
        rel = (math.degrees(math.atan2(kc[1], kc[0]))
               - math.degrees(math.atan2(recip[0][1], recip[0][0]))) % 60.0
        r = float(np.linalg.norm(kc))
        if min(rel, 60.0 - rel) < 2.0:                       # Γ–M 方向
            return "M" if abs(r - b / 2.0) < 1e-2 * b else "Γ–M"
        if abs(rel - 30.0) < 2.0:                            # Γ–K 方向
            rk = b / math.sqrt(3.0)
            if abs(r - rk) < 1e-2 * rk:
                return "K"
            if r < rk:
                return "Q"
    else:
        two = 2.0 * k
        if np.all(np.abs(two - np.round(two)) < 1e-6):
            return {(1, 0): "X", (0, 1): "Y", (1, 1): "S"}.get(
                (int(abs(round(two[0]))), int(abs(round(two[1])))), "k")
    return "k(%.3f,%.3f)" % (k[0], k[1])


class _DeformE1:
    """[0022] 逐谷形变势：与 S7.1 的 E1_vac 同一定义（gen_step9b_deform_read.py 的 _edge_raw / dvac 段）：
        E1_α(k, 带, 自旋) = | [(E₊−E₀)/γ + (E₋−E₀)/(−γ)]/2 − dE_vac/dε_α |
    E 取 step7_deform/<构型>/(ionrelax/)vasprun.xml —— band_edges.json 的 read_paths 记着 S7.1 实际读的哪份；
    γ、dE_vac/dε、±应变配对目录都取 band_edges.json 的 vac_align（S7.1 落盘的同一组数）。"""
    AX = ("xx", "yy")

    def __init__(self, cwd, be):
        va = (be or {}).get("vac_align") or {}
        if va.get("status") != "ok":
            raise _MVStop("band_edges.json 的真空对齐不是 ok（%s）：逐谷 E1 用真空口径，先让 S7.1 的真空对齐跑通"
                          % (va.get("reason") or va.get("status") or "无 vac_align 段，旧版 S7.1"))
        try:
            self.pairs = {a: [str(x) for x in va["pairs"][a]] for a in self.AX}
            self.g = {a: float(va["strain_mag"][a]) for a in self.AX}
            self.dvac = {a: float(va["dvac_eV_per_unit_strain"][a]) for a in self.AX}
        except (KeyError, TypeError, ValueError):
            raise _MVStop("band_edges.json 的 vac_align 缺 pairs/strain_mag/dvac_eV_per_unit_strain"
                          "（旧版 S7.1）：重跑 step7b_deform_read")
        self.read_paths = be.get("read_paths") or {}
        self.root = Path(cwd) / DEFORM_DIR
        self.used = {}

    def _path(self, folder):
        rp = (self.read_paths.get(folder) or {}).get("vasprun.xml")
        dirs = {"ionrelax": [self.root / folder / "ionrelax"], "rigid": [self.root / folder]}.get(
            rp, [self.root / folder / "ionrelax", self.root / folder])
        for d in dirs:
            for n in ("vasprun.xml", "vasprun.xml.gz"):
                if (d / n).is_file():
                    return d / n
        raise _MVStop("缺 %s/%s/%svasprun.xml（S7.1 读的那份，逐谷 E1 要它；S7 目录被清理过？）"
                      % (DEFORM_DIR, folder, "ionrelax/" if rp == "ionrelax" else ""))

    def data(self, folder):
        path = self._path(folder)
        self.used[folder] = str(path)
        return _cached_eig(path)

    def energy(self, folder, kf, band, spin):
        kp, eig = self.data(folder)
        i = _kidx_tr(kp, kf)
        if i is None:
            raise _MVStop("k=%s 不在 %s 的 k 网格（含 −k）上" % ([round(float(x), 4) for x in kf], folder))
        if spin not in eig or band >= eig[spin].shape[1]:
            raise _MVStop("%s 没有自旋 %s 的第 %d 条带" % (folder, spin, band))
        return float(eig[spin][i, band, 0])

    def e1_parts(self, ax, kf, band, spin):
        """[0023] (raw, signed)：raw = 本征值对应变的中心差分（VASP 本征值以胞内平均静电势为零点），
        signed = raw − dE_vac/dε（真空口径，带符号）。E1 = |signed|。两者都写进 json：真空口径下 K、Q 谷
        的形变势若一正一负被同一个 dE_vac/dε 平移，小的那个可能被移到 0 附近，看 raw 就知道是不是这种抵消。"""
        dp_, dm_ = self.pairs[ax]
        g = self.g[ax] if self.g[ax] > 1e-6 else 0.005
        e0 = self.energy("undeformed", kf, band, spin)
        raw = ((self.energy(dp_, kf, band, spin) - e0) / g
               + (self.energy(dm_, kf, band, spin) - e0) / (-g)) / 2.0
        return raw, raw - self.dvac[ax]

    def e1(self, ax, kf, band, spin):
        return abs(self.e1_parts(ax, kf, band, spin)[1])

    def snap(self, kf, recip):
        """S3 网格上的谷 k -> S7 网格上最近的点（ISYM=0 的列表 ∪ −列表 = 全网格）。返回 (k, 偏差 Å⁻¹)。"""
        import numpy as np
        kp, _ = self.data("undeformed")
        cand = np.vstack([kp, -kp])
        d = cand - np.asarray(kf, float)
        d -= np.round(d)
        dc = np.linalg.norm(d @ np.asarray(recip), axis=1)
        i = int(np.argmin(dc))
        return np.asarray(kf, float) + d[i], float(dc[i])

    def nocc(self, spin):
        import numpy as np
        _kp, eig = self.data("undeformed")
        if spin not in eig:
            raise _MVStop("S7 undeformed 没有自旋道 %s（S3 与 S7 的 ISPIN 不同）" % spin)
        return int(np.count_nonzero(eig[spin][0, :, 1] > 0.5))

    def selfcheck(self, be, carrier):
        """带边 E1 必须复现 S7.1 写的 E1_vac_*（同一 k、同一组简并带）：不一致说明构型/带号/自旋对错了。"""
        e = be.get(carrier) or {}
        hits = e.get("hits") or []
        if not hits or hits[0].get("k_frac") is None:
            raise _MVStop("band_edges.json 没有 %s 的带边 hits" % carrier)
        kf = hits[0]["k_frac"]
        sb = sorted({(str(h["spin"]).lower(), int(h["band_global"])) for h in hits})
        out = {}
        for ax in self.AX:
            ref = e.get("E1_vac_%s_raw_eV" % ax, e.get("E1_vac_%s_eV" % ax))
            if ref is None:
                raise _MVStop("band_edges.json 没有 %s 的 E1_vac_%s：S7.1 真空对齐没出数" % (carrier, ax))
            val = sum(self.e1(ax, kf, b, sp) for sp, b in sb) / len(sb)
            out[ax] = {"recomputed_eV": round(val, 4), "S7.1_eV": float(ref)}
            if abs(val - float(ref)) > max(0.02, 0.02 * abs(float(ref))):
                raise _MVStop("带边 E1_%s 重算 %.4f eV ≠ S7.1 的 %.4f eV：读的构型/带号/自旋与 S7.1 不一致，"
                              "逐谷 E1 不可信（S7 目录在 S7.1 之后被改过？重跑 step7b_deform_read）"
                              % (ax, val, float(ref)))
        e_edge = self.energy("undeformed", kf, sb[0][1], sb[0][0])
        return out, e_edge


def _mv_mobility(srecs, Cxy, T, iso_crystal, only=None):
    """[0023] 多谷加权（原是 _multi_valley 里的闭包，挪出来给文献对照测试直接调）。
    srecs：[{"label", "members": [{"m_d", "dE", "E1": [x, y], "W": 2×2 逆质量张量(1/m0)}]}]；only：只算这个星。
    μ_αα = (e/m0)·Σ n τ_α W_αα / Σ n，n = m_d·exp(−ΔE/k_BT)，τ_α = ħ³C_α/(k_BT m_d E1_α²)。
    返回 {"mu": [x, y] cm²/Vs, "tau": [x, y] s, "m_c": [x, y] m0, "w": {星: 布居}}。"""
    kT = KB * T / E_C
    sn, snw, sntw, per = 0.0, [0.0, 0.0], [0.0, 0.0], {}
    for s in srecs:
        if only is not None and s is not only:
            continue
        ns = 0.0
        for m in s["members"]:
            n = m["m_d"] * math.exp(-m["dE"] / kT)
            ns += n
            for a in (0, 1):
                tau = HBAR ** 3 * Cxy[a] / (KB * T * (m["m_d"] * M0) * (m["E1"][a] * E_C) ** 2)
                snw[a] += n * m["W"][a][a]
                sntw[a] += n * tau * m["W"][a][a]
        per[s["label"]] = ns
        sn += ns
    if iso_crystal:                                  # 六方/四方：面内输运张量各向同性
        snw = [sum(snw) / 2.0] * 2
        sntw = [sum(sntw) / 2.0] * 2
    mu = [E_C / M0 * sntw[a] / sn * 1e4 for a in (0, 1)]
    return {"mu": mu, "tau": [sntw[a] / snw[a] for a in (0, 1)],
            "m_c": [sn / snw[a] for a in (0, 1)], "w": {k: x / sn for k, x in per.items()}}


def _mv_e1_eff(star):
    """[0023] 星的等效形变势：τ ∝ ⟨1/E1²⟩（成员 × 方向平均）-> E1_eff = ⟨1/E1²⟩^(-1/2)。"""
    vals = [1.0 / (e * e) for m in star["members"] for e in m["E1"]]
    return (sum(vals) / len(vals)) ** -0.5


def _mv_e1_gate(srecs, Cxy, temps, iso_crystal):
    """[0023] 有贡献的星（任一温度布居 ≥ VALLEY_W_MIN）之间 E1_eff 相差 > VALLEY_E1_RATIO_MAX 倍 -> (说明, 明细)；
    否则 (None, 明细)。明细 = {星: {"E1_eff_eV", "pop_max"}}。"""
    pmax = {}
    for T in temps:
        for k, x in _mv_mobility(srecs, Cxy, T, iso_crystal)["w"].items():
            pmax[k] = max(pmax.get(k, 0.0), x)
    info = {s["label"]: {"E1_eff_eV": round(_mv_e1_eff(s), 4), "pop_max": round(pmax[s["label"]], 4)}
            for s in srecs}
    contrib = [k for k in info if info[k]["pop_max"] >= VALLEY_W_MIN]
    if len(contrib) < 2:
        return None, info
    hi = max(contrib, key=lambda k: info[k]["E1_eff_eV"])
    lo = min(contrib, key=lambda k: info[k]["E1_eff_eV"])
    ratio = info[hi]["E1_eff_eV"] / info[lo]["E1_eff_eV"]
    if ratio <= VALLEY_E1_RATIO_MAX:
        return None, info
    return ("拒绝：有贡献的谷之间形变势相差 %.1f 倍（%s E1=%.2f eV、最大布居 %.1f%%；%s E1=%.2f eV；阈值 "
            "VALLEY_E1_RATIO_MAX=%g）。τ∝1/E1²，E1 小的谷 τ 长约 %.0f 倍，布居再小也会主导电导；它的 LA 膨胀形变势"
            "耦合几乎抵消，真实输运由 DPT 不含的 TA/谷间/光学散射决定（Jin et al. PRB 90, 045422：DFPT 拟合的 Q/K "
            "声学形变势比 0.6–0.9，WSe2 电子含谷间 30、只 K 250 cm²/Vs），多谷 DPT 给不出可信数。"
            "单谷结果（若有）照旧可用；各谷带符号的 dE/dε 见 valleys[].E1_signed_xy_eV / dEdeps_raw_xy_eV"
            % (ratio, lo, info[lo]["E1_eff_eV"], 100 * info[lo]["pop_max"], hi, info[hi]["E1_eff_eV"],
               VALLEY_E1_RATIO_MAX, ratio ** 2)), info


def _multi_valley_block(cwd, is_2d, carrier, T0, temps, by_direction=None):
    """[0022] 多谷 DPT 结果块（2D）。拒绝出数时只写 status（原因），绝不崩。"""
    if not is_2d:
        return {"status": "多谷 DPT 只实现 2D（3D 用上面的各向同性值）"}
    if str(DPT_VALLEYS).strip().lower() == "off":
        return {"status": "off（本材料 step.conf DPT_VALLEYS = off）"}
    try:
        return _multi_valley(cwd, carrier, T0, temps, by_direction)
    except _MVStop as e:
        return {"status": str(e)}
    except Exception as e:                       # noqa: BLE001
        return {"status": "多谷 DPT 异常：%s" % _exc_brief(e)}


def _multi_valley(cwd, carrier, T0, temps, by_direction):
    import numpy as np
    if E1_SOURCE != "vac":
        raise _MVStop("多谷 DPT 只实现真空口径 E1（E1_SOURCE=vac）；当前 E1_SOURCE=%s" % E1_SOURCE)
    Cxy, Cprov = get_C_aniso(cwd, True)
    if Cxy is None:
        raise _MVStop("缺 C：%s" % Cprov)
    be = _band_edges_json(cwd)
    if not be:
        raise _MVStop("缺 %s/band_edges.json" % DEFORM_READ_DIR)
    vr = _vasprun_s3(cwd)
    if vr is None:
        raise _MVStop("缺 step3_uniform/vasprun.xml")
    v = _cached_vasprun(vr)
    st = v.final_structure
    recip = np.asarray(st.lattice.reciprocal_lattice.matrix)
    hexcell = _hex_cell(recip)
    sgn = 1.0 if carrier == "electron" else -1.0
    win = float(VALLEY_WINDOW_EV)

    spins = _spin_items(v)
    nocc = {}
    for sp, arr in spins:
        n_k = np.count_nonzero(arr[:, :, 1] > 0.5, axis=1)
        if int(n_k.min()) != int(n_k.max()):
            raise _MVStop("自旋 %s 各 k 点的占据带数不同（%d..%d）：能带交叠（金属/半金属），多谷 DPT 不适用"
                          % (sp, int(n_k.min()), int(n_k.max())))
        nocc[sp] = int(n_k[0])
    # [0028] 非磁体系用 ISPIN=2 算时两个自旋道的能带逐点相同：同一个谷会按 up/down 各算一次，
    #   布居对半分（CrS2/CrSe2 hex 实测"单谷只代表 K 谷 50%"，另一半其实是同一个 K 谷的另一自旋）。
    #   两道占据带数相同、能带处处差 < 1 meV 就只留一道（μ 不变，布居与占比回到正确值）。
    if len(spins) == 2 and nocc[spins[0][0]] == nocc[spins[1][0]]:
        _e0, _e1 = spins[0][1][:, :, 0], spins[1][1][:, :, 0]
        if _e0.shape == _e1.shape and float(np.max(np.abs(_e0 - _e1))) < 1e-3:
            print("[..] %s：两个自旋道能带逐点相同（非磁），多谷只按一道算" % carrier)
            spins = spins[:1]
            nocc = {spins[0][0]: nocc[spins[0][0]]}

    kfrac, kp_map, xmethod, xerr = _expand_full_bz(v, cwd, carrier)
    if xerr:
        raise _MVStop(xerr)
    plane, perr = _mv_plane(kfrac)
    if perr:
        raise _MVStop(perr)
    N, idx = plane["N"], plane["idx"]
    ene = {sp: arr[np.asarray(kp_map), :, 0] for sp, arr in spins}          # (n_full, nb)

    def bands_of(sp):
        nb = ene[sp].shape[1]
        return range(nocc[sp], nb) if sgn > 0 else range(nocc[sp] - 1, -1, -1)
    edge = None
    for sp, _arr in spins:
        for b in bands_of(sp):
            e = float((sgn * ene[sp][idx, b]).min())
            edge = e if edge is None else min(edge, e)
            break
    if edge is None:
        raise _MVStop("没有%s带" % ("导" if sgn > 0 else "价"))

    valleys, seeds, skipped = _mv_collect(ene, idx, N, nocc, sgn, edge, win)
    kops = _mv_kops(st, _time_reversal_ok(v, cwd, carrier))
    stars, serr = _mv_stars(valleys, N, kops, seeds)
    if serr:
        raise _MVStop(serr)
    iso_crystal = any(_mv_rot3(r) for r in kops)

    rd = _DeformE1(cwd, be)
    for sp in nocc:
        if rd.nocc(sp) != nocc[sp]:
            raise _MVStop("S3 与 S7 undeformed 的占据带数不同（自旋 %s：%d vs %d）：NELECT/SOC 不一致，带号对不上"
                          % (sp, nocc[sp], rd.nocc(sp)))
    selfcheck, e7_edge = rd.selfcheck(be, carrier)

    srecs, warns = [], []
    for star in stars:
        members = []
        for vi in star:
            sp, b, pts, e_grid = valleys[vi]
            k0 = int(idx[pts[0]])
            tag = "%s 带 %d k=%s" % (sp, b, np.round(kfrac[k0][:2], 4).tolist())
            fit, ferr = _quad_fit_2d(kfrac, recip, ene[sp][:, b], k0, kops=kops)
            if ferr:
                raise _MVStop("谷 %s（离带边 %.3f eV）拟合失败：%s —— 它在窗口 VALLEY_WINDOW_EV=%.2f 里不能丢；"
                              "确认它不重要就在本材料 step.conf 把 VALLEY_WINDOW_EV 调到它以下"
                              % (tag, float(sgn * e_grid) - edge, ferr, win))
            a, bq, c = fit["a"], fit["b"], fit["c"]
            if sgn * a <= 0:
                raise _MVStop("谷 %s 的拟合曲率符号与%s不符（a=%.3f）" % (tag, "极小" if sgn > 0 else "极大", a))
            W = sgn * np.array([[a, c / 2.0], [c / 2.0, bq]]) / HB2_2M0              # 1/m0
            M = np.linalg.inv(W)
            co = fit["coef"]
            e_fit, xs = e_grid, np.zeros(2)
            if "d" in co:                                        # 极值不在网格点上：取拟合出的极值（能量与位置）
                H = np.array([[2 * a, c], [c, 2 * bq]])
                gl = np.array([co["d"], co["e"]])
                xs = -np.linalg.solve(H, gl)
                e_fit = e_grid + co["f"] + 0.5 * float(gl @ xs)
            shift = float(np.linalg.norm(xs))
            if shift > 1.5 * fit["step_min"]:                   # 网格极值离真极值不会超过 ~1 步：拟合与网格矛盾
                raise _MVStop("谷 %s 的拟合极值离网格点 %.3f Å⁻¹（> 1.5 步 = %.3f）：与网格上的极值矛盾，拟合不可信"
                              % (tag, shift, 1.5 * fit["step_min"]))
            kf3 = np.array(kfrac[k0], float)
            kc_fit = kf3 @ recip + np.array([xs[0], xs[1], 0.0])
            kf_fit = np.linalg.solve(recip.T, kc_fit)            # 拟合极值的分数坐标（标注、取 E1 都用它）
            ks7, koff = rd.snap(kf_fit, recip)
            parts = [rd.e1_parts(ax, ks7, b, sp) for ax in ("xx", "yy")]
            e1 = [abs(q[1]) for q in parts]
            e1_raw = list(e1)
            if _mv_site_iso(kf3, kops):                         # K/Γ（网格点就是该点）：C₃ 使面内形变势张量各向同性
                e1 = [(e1[0] + e1[1]) / 2.0] * 2
            if min(e1) < 1e-3:
                raise _MVStop("谷 %s 的形变势 E1=(%.4f, %.4f) eV ≈ 0：DPT 发散，该谷不能用 DPT" % (tag, e1[0], e1[1]))
            d7 = sgn * (rd.energy("undeformed", ks7, b, sp) - e7_edge)
            members.append({"spin": sp, "band": b, "k": kf_fit, "k_grid": kf3, "e_grid": e_grid, "e_fit": e_fit,
                            "W": W, "M": M,
                            "m_d": 1.0 / math.sqrt(float(np.linalg.det(W))), "E1": e1, "E1_raw": e1_raw,
                            "E1_signed": [q[1] for q in parts], "dEdeps_raw": [q[0] for q in parts],
                            "k_S7": ks7, "k_offset_S7": koff, "dE_S7": d7, "shift": shift,
                            "fit": "%d点/%d壳层, 模型=%s, rel=%.3f" % (fit["npt"], fit["n_shell"],
                                                                     "+".join(fit["names"]), fit["rel"])})
        srecs.append({"members": members})

    ref = min(sgn * m["e_fit"] for s in srecs for m in s["members"])
    used_labels = set()
    for s in srecs:
        m0 = s["members"][0]
        for m in s["members"]:
            m["dE"] = sgn * m["e_fit"] - ref
        lab = _mv_label(m0["k"], recip, hexcell)
        if lab in used_labels:
            lab = "%s(%s带%d)" % (lab, m0["spin"], m0["band"])
        used_labels.add(lab)
        s["label"] = lab
        dEs = [m["dE"] for m in s["members"]]
        if max(dEs) - min(dEs) > 0.005:
            warns.append("星 %s 的等价谷能量相差 %.1f meV（网格/对称性不一致？）" % (lab, 1000 * (max(dEs) - min(dEs))))
        for m in s["members"]:
            dd = abs(m["dE_S7"] - (sgn * m["e_grid"] - edge))
            if dd > 0.05:
                warns.append("星 %s：S7 网格上的谷能量与 S3 差 %.0f meV（S7 网格太粗、谷 k 偏 %.3f Å⁻¹）"
                             % (lab, 1000 * dd, m["k_offset_S7"]))
                break

    def at_T(T, only=None):
        return _mv_mobility(srecs, Cxy, T, iso_crystal, only)

    r0 = at_T(T0)
    gate, e1info = _mv_e1_gate(srecs, Cxy, sorted(set(temps) | {T0}), iso_crystal)
    # [0024] 单谷 DPT 拟合的是 S3 网格上的带边点 -> 它所在的星；各温度布居与 E1 无关（拒绝时也照样可信）
    grid_edge = min(srecs, key=lambda s: min(sgn * m["e_grid"] for m in s["members"]))
    pops = [{"T_K": T, "weights": {k: round(float(x), 4) for k, x in at_T(T)["w"].items()}}
            for T in sorted(set(temps) | {T0})]
    rows = []
    for T in temps:
        r = at_T(T)
        rows.append({"T_K": T,
                     "mobility_x_cm2_Vs": round(float(r["mu"][0]), 3),
                     "mobility_y_cm2_Vs": round(float(r["mu"][1]), 3),
                     "mobility_inplane_cm2_Vs": round(float(r["mu"][0] + r["mu"][1]) / 2.0, 3),
                     "tau_x_s": float(r["tau"][0]), "tau_y_s": float(r["tau"][1]),
                     "weights": {k: round(float(x), 4) for k, x in r["w"].items()}})

    # 自洽核对：只留带边所在的那个星时应等于 by_direction（同一拟合、同一 E1 定义）
    cons = None
    bd = by_direction or {}
    if bd.get("status") == "ok":
        rs = at_T(T0, only=grid_edge)
        cons = {"star": grid_edge["label"],
                "x": round(float(rs["mu"][0] / bd["x"]["mobility_cm2_Vs"]), 4),
                "y": round(float(rs["mu"][1] / bd["y"]["mobility_cm2_Vs"]), 4)}
        if max(abs(cons["x"] - 1), abs(cons["y"] - 1)) > 0.10:
            warns.append("只留 %s 星的多谷 μ 与 by_direction 差 >10%%（x %.3f, y %.3f）：同一个谷两条路径不自洽"
                         % (grid_edge["label"], cons["x"], cons["y"]))

    vout = []
    for s in sorted(srecs, key=lambda s: min(m["dE"] for m in s["members"])):
        m0 = s["members"][0]
        vout.append({
            "label": s["label"], "g": len(s["members"]), "spin": m0["spin"], "band_index": int(m0["band"]),
            "dE_eV": round(float(min(m["dE"] for m in s["members"])), 4),
            "k_frac": [np.round(m["k"], 6).tolist() for m in s["members"]],
            "k_grid": [np.round(m["k_grid"], 6).tolist() for m in s["members"]],
            "m_d_m0": round(float(np.mean([m["m_d"] for m in s["members"]])), 4),
            "M_m0": np.round(m0["M"], 4).tolist(),
            "E1_xy_eV": [[round(float(x), 4) for x in m["E1"]] for m in s["members"]],
            "E1_xy_raw_eV": [[round(float(x), 4) for x in m["E1_raw"]] for m in s["members"]],
            "E1_signed_xy_eV": [[round(float(x), 4) for x in m["E1_signed"]] for m in s["members"]],
            "dEdeps_raw_xy_eV": [[round(float(x), 4) for x in m["dEdeps_raw"]] for m in s["members"]],
            "E1_eff_eV": e1info[s["label"]]["E1_eff_eV"], "pop_max": e1info[s["label"]]["pop_max"],
            "k_S7": [np.round(m["k_S7"], 6).tolist() for m in s["members"]],
            "k_offset_S7_invA": round(float(max(m["k_offset_S7"] for m in s["members"])), 4),
            "extremum_shift_invA": round(float(max(m["shift"] for m in s["members"])), 4),
            "fit": m0["fit"],
            "weight_T0": round(float(r0["w"][s["label"]]), 4)})
    trunc = [T for T in temps if math.exp(-win / (KB * T / E_C)) > 0.01]
    for w_ in warns:
        print("[WARN] %s 多谷：%s" % (carrier, w_))
    if gate:                                             # [0023] 只留诊断，不给任何迁移率/τ（下游只认 status ok）
        print("[WARN] %s 多谷：%s" % (carrier, gate))
        return {"status": gate, "window_eV": win, "valleys": vout,
                "edge_star": grid_edge["label"], "populations_vs_T": pops,
                "dvac_eV_per_unit_strain": dict(rd.dvac), "E1_selfcheck": selfcheck,
                "S7_files": dict(rd.used), "warnings": warns}
    return {
        "status": "ok",
        "method": ("多谷 DPT：S3 全 BZ（kz=0 层）找窗口内各带各自旋的谷、按点群(+TR)归星；逐谷二次型拟合 W=M⁻¹、"
                   "m_d=√det M；逐谷 E1_α=|dE/dε_αα − dE_vac/dε_αα|（S7 形变构型本征值中心差分，与 S7.1 E1_vac 同一定义）；"
                   "μ_αα(T)=(e/m0)Σ n_v τ_α,v W_αα,v/Σ n_v，n_v∝m_d,v·exp(−ΔE_v/k_BT)，τ_α,v=ħ³C_α/(k_BT m_d,v E1_α,v²)"),
        "window_eV": win, "expand": xmethod, "grid_plane": list(N),
        "C_2D_N_per_m": [round(Cxy[0], 3), round(Cxy[1], 3)], "C_provenance": Cprov,
        "valleys": vout, "T_K": T0,
        "x": {"mobility_cm2_Vs": round(float(r0["mu"][0]), 3), "tau_s": float(r0["tau"][0]),
              "m_c_m0": round(float(r0["m_c"][0]), 4)},
        "y": {"mobility_cm2_Vs": round(float(r0["mu"][1]), 3), "tau_s": float(r0["tau"][1]),
              "m_c_m0": round(float(r0["m_c"][1]), 4)},
        "mobility_inplane_cm2_Vs": round(float(r0["mu"][0] + r0["mu"][1]) / 2.0, 3),
        "inplane_isotropic_by_symmetry": bool(iso_crystal),
        "mobility_vs_T": rows,
        "window_truncation_T_K": trunc,
        "E1_selfcheck": selfcheck, "S7_files": dict(rd.used),
        "dvac_eV_per_unit_strain": dict(rd.dvac),
        "edge_star": grid_edge["label"], "populations_vs_T": pops,
        "edge_star_vs_by_direction": cons,
        "skipped_band_crossings": skipped,
        "warnings": warns,
        "note": ("不含谷间散射（K↔Q、K↔K' 的大波矢声子），只有谷内 LA 形变势散射，是上限；"
                 "窗口外的谷不计（window_truncation_T_K 列出窗口边缘玻尔兹曼因子 >1% 的温度）。"),
    }


def _single_valley_scope(mv, T0):
    """[0024] 单谷 DPT（各向同性值与 by_direction）只算带边所在的那组谷。多谷块找到别的非等价谷时（不论多谷出没出数），
    写明单谷值代表的谷在各温度占多少载流子 —— 只报布居、不设阈值、不拦截：K 谷单谷值作为与文献同口径的对照仍然有用，
    但 MoSe2 这类 300 K 时六成电子在 Q 谷的材料，它不能当成材料的电子迁移率（Jin PRB 90, 045422：MoSe2 电子
    含谷间 25、只 K 180 cm²/Vs）。返回 None = 只有一组谷 / 多谷块没找谷。"""
    edge, pv = (mv or {}).get("edge_star"), (mv or {}).get("populations_vs_T") or []
    if not edge or not pv or len(pv[0]["weights"]) < 2:
        return None
    shares = [{"T_K": r["T_K"], "share": r["weights"].get(edge, 0.0)} for r in pv]
    r0 = min(pv, key=lambda r: abs(r["T_K"] - T0))
    lo = min(shares, key=lambda r: r["share"])
    others = {k: v for k, v in r0["weights"].items() if k != edge}
    note = ("单谷值只代表 %s 谷：%.0f K 时占载流子 %.1f%%（%s 未计入）%s；这些谷和谷间散射单谷 DPT 都不含，"
            "材料的迁移率以 S8.4（AMSET）为准"
            % (edge, r0["T_K"], 100 * r0["weights"][edge],
               "、".join("%s %.1f%%" % (k, 100 * v) for k, v in sorted(others.items(), key=lambda kv: -kv[1])),
               ("，%.0f K 时降到 %.1f%%" % (lo["T_K"], 100 * lo["share"])) if lo["T_K"] != r0["T_K"] else ""))
    return {"edge_star": edge, "T_K": r0["T_K"], "share": r0["weights"][edge], "others": others,
            "shares_vs_T": shares, "note": note}


# ============================ [0025] 3D 分方向 DPT ============================
# 2D 的 by_direction 只做面内 x/y；3D 以前只出各向同性值。这里补 x/y/z：
#   · C_α  = step6 的 C11 / C22 / C33（VASP 顺序 XX YY ZZ）；
#   · m*   = 带边三维二次型拟合 E = Σ a_ij Δk_i Δk_j（+q⁴、+线性项）-> 逆质量张量 W，m_α = 1/W_αα（电导质量），
#            m_d = (det M)^(1/3)；带边谷有对称等价谷时 W 先按点群平均（晶体方向的电导质量）；
#   · E1_α = deformation.h5（S7.1 产出，AMSET 原生输出）在带边处 3×3 张量的对角 D_αα，按 band_edges.json 的
#            每个带边点（band_h5、k_h5、spin）及其对称等价 k 点取，平均。3D 没有真空，用芯势口径；S7 的 amset deform create
#            三维本来就生成全部应变分量（含 zz），不用加 VASP 作业，S7.1 也不用改；
#   · μ_α = 2√(2π) e ħ⁴ C_α / (3 (k_BT)^{3/2} m_α m_d^{3/2} E1_α²)：散射率 ∝ 态密度 ∝ m_d^{3/2}，迁移率 = eτ/m_α；
#           m_α = m_d = m* 时就是各向同性的 m*^{5/2}。τ_α = μ_α m_α / e。
def _cart_rots(structure):
    """[0025] 点群的笛卡尔旋转，并上 −R（时间反演），去重。k 的笛卡尔分量与实空间同样变换：k' = R k。"""
    import numpy as np
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    rots = []
    for op in SpacegroupAnalyzer(structure, symprec=0.01).get_point_group_operations(cartesian=True):
        for r in (np.asarray(op.rotation_matrix, float), -np.asarray(op.rotation_matrix, float)):
            if not any(np.allclose(r, u, atol=1e-6) for u in rots):
                rots.append(r)
    return rots


def _kkey(kf, n=10000):
    import numpy as np
    return tuple(int(x) for x in np.mod(np.rint(np.asarray(kf, float) * n).astype(np.int64), n))


def _star_images(kf, structure, rots):
    """[0025] k（分数坐标）在点群 + 时间反演下的等价点（模倒格矢去重），分数坐标列表。"""
    import numpy as np
    B = np.asarray(structure.lattice.reciprocal_lattice.matrix, float)      # 行 = b_i
    Binv = np.linalg.inv(B)
    kc = np.asarray(kf, float) @ B
    out, seen = [], set()
    for r in rots:
        g = (r @ kc) @ Binv
        key = _kkey(g)
        if key not in seen:
            seen.add(key)
            out.append(g - np.floor(g + 1e-9))
    return out


def _e1_diag_3d(cwd, carrier, be=None):
    """带边处形变势张量的对角 (D_xx, D_yy, D_zz)，eV（|·|：AMSET 的 potentials.py 逐分量取 np.abs(diff/strain)）。
    带边谷有对称等价谷（Si 的 X、Bi₂Te₃ 的 6 个谷）时，晶体某方向的 E1 要对全部等价谷平均：h5 存的是逐分量绝对值，
    张量不能整体旋转，所以直接读 h5 网格上每个等价 k 点的原生值（h5 带 kpoints/structure 时）。
    返回 ((dx, dy, dz), 说明)；失败 (None, 原因)。"""
    be = be if be is not None else _band_edges_json(cwd)
    hits = ((be or {}).get(carrier) or {}).get("hits") or []
    if not hits:
        return None, "band_edges.json 没有 %s 带边 hits" % carrier
    h5 = Path(cwd) / DEFORM_READ_DIR / "deformation.h5"
    if not h5.is_file():
        return None, "缺 %s/deformation.h5" % DEFORM_READ_DIR
    try:
        import numpy as np
        import h5py
        vals, n_hit, n_miss, star_how = [], len(hits), 0, "h5 无 kpoints/structure，只取带边点本身"
        with h5py.File(str(h5), "r") as f:
            kp = st = rots = None
            if "kpoints" in f and "structure" in f:
                try:
                    from pymatgen.core import Structure
                    raw = np.array(f["structure"])
                    sj = raw.item() if raw.shape == () else raw
                    sj = sj.decode() if isinstance(sj, (bytes, np.bytes_)) else str(sj)
                    st = Structure.from_str(sj, fmt="json")
                    rots = _cart_rots(st)
                    kp = np.asarray(f["kpoints"], float)
                    index = {}
                    for i, k in enumerate(kp):
                        index.setdefault(_kkey(k), i)
                except Exception as e:           # noqa: BLE001
                    kp = None
                    star_how = "h5 的 structure/kpoints 读不了（%s），只取带边点本身" % _exc_brief(e)
            done = set()
            for h in hits:
                ds = "deformation_potentials_%s" % str(h.get("spin", "up")).lower()
                if ds not in f or h.get("band_h5") is None or h.get("k_h5") is None:
                    return None, "hit 缺 band_h5/k_h5 或 h5 没有 %s" % ds
                b, k = int(h["band_h5"]), int(h["k_h5"])
                ks = [k]
                if kp is not None:
                    ks = []
                    for g in _star_images(kp[k], st, rots):
                        j = index.get(_kkey(g))
                        if j is None:
                            n_miss += 1
                        else:
                            ks.append(j)
                    ks = ks or [k]
                for j in ks:
                    if (ds, b, j) in done:
                        continue
                    done.add((ds, b, j))
                    t = np.asarray(f[ds][b, j], float)
                    if t.shape != (3, 3):
                        return None, "h5 的形变势不是 3×3（shape=%s）" % (t.shape,)
                    vals.append([abs(float(t[i, i])) for i in range(3)])
            if kp is not None:
                star_how = "带边点 %d 个，连同对称等价点共 %d 个%s" % (
                    n_hit, len(vals), ("；%d 个对称像不在 h5 网格上" % n_miss) if n_miss else "")
        d = np.mean(np.asarray(vals), axis=0)
        return (round(float(d[0]), 4), round(float(d[1]), 4), round(float(d[2]), 4)), "n=%d，%s" % (
            len(vals), star_how)
    except Exception as e:                       # noqa: BLE001
        return None, "读 deformation.h5 异常：%s" % _exc_brief(e)


def _quad_fit_3d(kfrac, recip, eband, k0, need_lin=None, mesh=None, kops=None):
    """[0025] 带边三维二次型拟合（_quad_fit_2d 的三维版，判据同构）：
    ΔE = a x² + b y² + c z² + d xy + e yz + f zx（+ q⁴）（+ 线性 gx + hy + iz + j）。
    返回 ({"H": 3×3 Hessian（eV·Å²）, npt, n_shell, Rused, names, cond, rel, step_min}, None)；失败 (None, 原因)。"""
    import numpy as np
    kcart = kfrac @ recip
    dk = kcart - kcart[k0]
    de = eband - eband[k0]
    q2 = np.sum(dk ** 2, axis=1)
    qnz = np.sqrt(q2[q2 > 1e-9])
    if len(qnz) == 0:
        return None, "带边附近无有效 k 点"
    q_min = float(qnz.min())
    if mesh is None:
        try:
            from amset.electronic_structure.symmetry import get_mesh_from_kpoint_diff
            _mesh, _ = get_mesh_from_kpoint_diff(kfrac)
            mesh = np.rint(np.asarray(_mesh)).astype(int)
        except Exception:
            _df = kfrac - kfrac[k0]
            _df -= np.round(_df)
            _nz = np.abs(_df) > 1e-6
            mesh = np.array([int(round(1.0 / float(np.min(np.abs(_df[_nz[:, i], i])))))
                             if _nz[:, i].any() else 1 for i in range(3)])
    step_cart = np.array([np.linalg.norm(recip[i]) / max(int(mesh[i]), 1) for i in range(3)])
    if need_lin is None and kops:                                        # [0027] 按小群判断
        need_lin = not _k_pinned(kfrac[k0], kops, 3)
    if need_lin is None:
        bases = (0.0, 0.5, 1.0 / 3.0, 2.0 / 3.0, 0.25, 0.75)
        need_lin = not all(any(abs((c % 1.0) - b) < 1e-3 or abs((c % 1.0) - b - 1.0) < 1e-3 for b in bases)
                           for c in kfrac[k0])
    sel = Rused = n_shell = None
    for NSTEP in ((FORCE_NSTEP,) if FORCE_NSTEP else (2, 3, 4, 5)):
        Rmax = (NSTEP - 0.1) * float(step_cart.min())
        m = (q2 < Rmax * Rmax) & (q2 > 1e-9)
        n = int(np.count_nonzero(m))
        if n < 14:                                   # 6 个二次项 + q⁴ 的 2 倍冗余
            continue
        n_sh = int(len(np.unique(np.round(np.sqrt(q2[m]) / q_min, 3))))
        ok = n_sh >= 2 and (n >= 30 or not need_lin)   # 非高对称点：线性项要 3 倍冗余（10 参数）
        if ok:
            sel, Rused, n_shell = m, Rmax, n_sh
            break
    if sel is None:
        return None, ("带边附近点太少或全在单壳层（需 ≥14 点 + ≥2 壳层；非高对称点 ≥30）——加密 step3_uniform 网格")
    npt = int(np.count_nonzero(sel))
    steps = np.rint(((kfrac[sel] - kfrac[k0]) - np.round(kfrac[sel] - kfrac[k0])) * mesh).astype(int)
    for i in range(3):
        if int(np.abs(steps[:, i]).max()) < 1:
            return None, "倒格矢方向 %d 无非零步取点（N=%d 过少），该方向曲率无约束" % (i, int(mesh[i]))
    x, y, z = dk[sel, 0], dk[sel, 1], dk[sel, 2]
    cols = [x * x, y * y, z * z, x * y, y * z, z * x]
    names = ["a", "b", "c", "d", "e", "f"]
    if n_shell >= 2:
        cols.append(q2[sel] ** 2)
        names.append("q4")
    if npt >= 30:
        cols += [x, y, z, np.ones(npt)]
        names += ["g", "h", "i", "j"]
    A = np.column_stack(cols)
    if A.shape[0] < 2 * A.shape[1]:
        return None, "拟合点数不足：%d 点 / %d 参数（需 ≥2 倍冗余）" % (A.shape[0], A.shape[1])
    cond = float(np.linalg.cond(A))
    if cond > 1e5:
        return None, "设计矩阵病态 cond=%.1e" % cond
    sol, *_ = np.linalg.lstsq(A, de[sel], rcond=None)
    rel = float(np.linalg.norm(de[sel] - A @ sol) / (np.linalg.norm(de[sel]) + 1e-12))
    if rel > 0.2:
        return None, "三维二次型拟合残差过大 rel=%.2f（带边附近非二次型）" % rel
    a, b, c, d, e, f = (float(v) for v in sol[:6])
    H = np.array([[2 * a, d, f], [d, 2 * b, e], [f, e, 2 * c]])
    return {"H": H, "npt": npt, "n_shell": n_shell, "Rused": Rused, "names": names, "cond": cond, "rel": rel,
            "step_min": float(step_cart.min())}, None


def _mass_from_hessian_3d(H, carrier, rots=None):
    """[0025] Hessian（eV·Å²）-> 逆质量张量 W（1/m0）-> {"m": 电导质量 1/W_αα, "m_d": (det M)^(1/3), "M"}。
    W 必须正定（导带向上弯、价带向下弯）；否则 (None, 原因)。
    rots（点群笛卡尔旋转 + 时间反演）给出时，"m" 取对全部等价谷平均的电导质量 1/⟨R W Rᵀ⟩_αα（对全群平均 =
    对星里每个谷等权平均）：Si 的 X 谷单谷是 m_l/m_t/m_t，晶体各方向都是 3/(1/m_l + 2/m_t)。"m_valley" 是带边那一个谷
    自己的值；m_d 是行列式，旋转不变。"""
    import numpy as np
    sgn = 1.0 if carrier == "electron" else -1.0
    W = sgn * np.asarray(H, float) / (2 * 3.80998)            # ħ²/m0 = 7.61996 eV·Å²
    ev = np.linalg.eigvalsh((W + W.T) / 2.0)
    if ev.min() <= 0:
        return None, "带边曲率不全是%s号（特征值 %s eV·Å²：鞍点？带边定位可疑）" % (
            "正" if sgn > 0 else "负", np.round(sgn * ev * 2 * 3.80998, 3).tolist())
    M = np.linalg.inv(W)
    Ws = (sum(r @ W @ r.T for r in rots) / len(rots)) if rots else W
    return {"m": tuple(round(1.0 / float(Ws[i, i]), 4) for i in range(3)),
            "m_valley": tuple(round(1.0 / float(W[i, i]), 4) for i in range(3)),
            "m_d": round(float(np.linalg.det(M)) ** (1.0 / 3.0), 4),
            "M": np.round(M, 4).tolist()}, None


def get_effective_mass_3d(cwd, carrier):
    """[0025] 3D 带边质量张量。返回 ({"m": (m_x, m_y, m_z) 电导质量 1/W_αα, "m_d": (det M)^(1/3),
    "M": 3×3 质量张量}, 说明)；失败 (None, 原因)。"""
    vr = _vasprun_s3(cwd)
    if vr is None:
        return None, "缺 step3_uniform/vasprun.xml"
    try:
        import numpy as np
        from pymatgen.io.vasp import Vasprun
        v = Vasprun(str(vr), parse_dos=False, parse_potcar_file=False)
        recip = np.asarray(v.final_structure.lattice.reciprocal_lattice.matrix)
        kfrac_ibz = np.array(v.actual_kpoints)
        hit = _find_band_edge(v, carrier)
        if hit is None:
            return None, "无未占据/占据态 —— 体系是金属？"
        isp, k0_ibz, b0, ene_ibz = hit
        kfrac, kp_mapping, how, xerr = _expand_full_bz(v, cwd, carrier)
        if xerr:
            return None, xerr
        ene = ene_ibz[np.asarray(kp_mapping)]
        dd = kfrac - kfrac_ibz[k0_ibz]
        dd -= np.round(dd)
        k0 = int(np.argmin(np.linalg.norm(dd, axis=1)))
        try:                                                             # [0027]
            kops3 = _kops_3d(v.final_structure, _time_reversal_ok(v, cwd, carrier))
        except Exception:                                                # noqa: BLE001
            kops3 = None
        fit, ferr = _quad_fit_3d(kfrac, recip, ene[:, b0], k0, kops=kops3)
        if ferr:
            return None, ferr
        rots = _cart_rots(v.final_structure)
        mres, merr = _mass_from_hessian_3d(fit["H"], carrier, rots)
        if merr:
            return None, merr
        mres["n_star"] = len(_star_images(kfrac[k0], v.final_structure, rots))
        mres["k0_frac"] = np.round(np.asarray(kfrac[k0], float) - np.round(kfrac[k0]), 4).tolist()
        star = ("；带边 k=%s 有 %d 个等价谷，m_α 取对它们平均的电导质量（单谷 %s）"
                % (mres["k0_frac"], mres["n_star"], list(mres["m_valley"]))) if mres["n_star"] > 1 else ""
        return mres, (
            "带边三维二次型拟合(%s, %d点/%d壳层, R<%.2f, 模型=%s, cond=%.0e, rel=%.3f, spin=%d)%s"
            % (how, fit["npt"], fit["n_shell"], fit["Rused"], "+".join(fit["names"]), fit["cond"], fit["rel"], isp,
               star))
    except Exception as e:                       # noqa: BLE001
        return None, "三维二次型拟合异常：%s" % _exc_brief(e)


def mobility_dpt_3d_aniso(C_alpha_Pa, m_alpha, m_d, E1_eV, T):
    """μ_α = 2√(2π) e ħ⁴ C_α / (3 (k_BT)^{3/2} m_α m_d^{3/2} E1_α²)，cm²/(V·s)。"""
    E1 = E1_eV * E_C
    mu = (2.0 * math.sqrt(2 * math.pi) * E_C * HBAR ** 4 * C_alpha_Pa) / (
        3.0 * (KB * T) ** 1.5 * (m_alpha * M0) * (m_d * M0) ** 1.5 * E1 ** 2)
    return mu * 1e4


def _aniso_block_3d(cwd, carrier, T):
    """[0025] 3D 分方向 DPT（x/y/z）。字段与 2D 的 by_direction 同构，C 用 C_3D_GPa。"""
    cd = read_elastic_diag_GPa(cwd)
    mres, mprov = get_effective_mass_3d(cwd, carrier)
    e3, eprov = _e1_diag_3d(cwd, carrier)
    blk = {"dim": "3d", "C_3D_GPa": list(cd) if cd else None, "C_provenance": "step6 C11/C22/C33",
           "m_eff_m0": list(mres["m"]) if mres else None, "m_provenance": mprov,
           "E1_eV": list(e3) if e3 else None,
           "E1_provenance": ("deformation.h5 带边 3×3 张量对角（芯势口径，%s）" % eprov) if e3 else eprov,
           "T_K": T,
           "formula": ("mu_a = 2*sqrt(2pi)*e*hbar^4*C_a/(3*(kB*T)^1.5*m_a*m_d^1.5*E1_a^2); "
                       "m_a = 1/W_aa（电导质量）, m_d = det(M)^(1/3); tau_a = mu_a*m_a/e")}
    miss = [n for n, v in (("C", cd), ("m*", mres), ("E1", e3)) if v is None]
    if miss:
        blk["status"] = "缺输入：%s" % "、".join(miss)
        return blk
    blk["m_d_m0"] = mres["m_d"]
    blk["M_m0"] = mres["M"]                      # 带边那一个谷的质量张量
    blk["m_valley_m0"] = list(mres["m_valley"])
    blk["n_equiv_valleys"] = mres.get("n_star", 1)
    for i, d in enumerate(("x", "y", "z")):
        mu = mobility_dpt_3d_aniso(cd[i] * 1e9, mres["m"][i], mres["m_d"], e3[i], T)
        blk[d] = {"C_3D_GPa": round(cd[i], 3), "m_eff_m0": mres["m"][i], "E1_eV": e3[i],
                  "mobility_cm2_Vs": round(mu, 3), "tau_s": (mu * 1e-4) * (mres["m"][i] * M0) / E_C}
    blk["status"] = "ok"
    return blk


def _grid_resolution_note(carrier, m, mprov, T):
    """[V147] 最近网格点高出带边 > NN_KT_MAX·kT -> m_provenance 注明并告警。返回新的 mprov。"""
    nn = _EDGE_NN_DE.get(carrier)
    kt = KB * T / E_C
    if m is None or str(mprov).startswith("manual") or nn is None or nn <= NN_KT_MAX * kt:
        return mprov
    print("[WARN] %s 的 m* 网格分辨不足：最近网格点高出带边 %.0f meV（%.1f kT），拟合量到的是能量平均的质量，"
          "带边非抛物时偏重（GaAs 实测 0.0523 vs 带边 0.030）。可靠值：在 S8 目录跑 "
          "amset eff-mass vasprun.xml -i <interpolation_factor> -d <掺杂> -t %g --bandgap <带隙> --average，"
          "写进本材料的 step.conf：tf -p <材料> -j step8.2_dpt conf --set M_EFF_%s=<m*>（不要改脚本的 MANUAL）"
          % (carrier, 1000 * nn, nn / kt, T, carrier.upper()))
    return ("%s；网格分辨不出带边曲率（最近网格点高出带边 %.0f meV = %.1f kT）：m* 是能量平均的值，带边非抛物时偏重"
            % (mprov, 1000 * nn, nn / kt))


def _one_carrier(cwd, is_2d, carrier, T, temps=None):
    C, Cunit, Cprov = get_C(cwd, is_2d)
    _EDGE_DEGEN.pop(carrier, None)
    _EDGE_NN_DE.pop(carrier, None)
    m, mprov = get_effective_mass(cwd, carrier, is_2d)
    _nd = _EDGE_DEGEN.get(carrier, 1)
    if m is not None and _nd > 1:                                   # [V145]
        mprov = "%s；带边 %d 重简并，取重带（单带 DPT 对简并带边只是近似）" % (mprov, _nd)
        print("[WARN] %s 带边 %d 重简并（如立方半导体 Γ 点的价带顶）：m* 取重带；单带 DPT 公式只是近似，"
              "与 AMSET 的差别不能当作问题" % (carrier, _nd))
    mprov = _grid_resolution_note(carrier, m, mprov, T)
    bd = _aniso_block(cwd, is_2d, carrier, T)                      # patch_dpt_aniso
    key = "m_eff_electron" if carrier == "electron" else "m_eff_hole"
    if (not is_2d) and bd.get("status") == "ok" and not MANUAL[key]:
        # [0025] 3D 表头 m*：旧的"方向分解 + 最轻两支几何平均"是 2D 的做法，三维各向异性时偏轻（质量 0.25/0.40/0.70
        #   时给 0.18，比最轻的主质量还轻，表头 μ 虚高约 7.5 倍；各向同性时两者一致）。三维二次型拟合成功就用它的
        #   态密度质量 m_d = (det M)^(1/3)；拟合失败才保留旧值。
        print("[..] %s：3D 表头 m* 改用三维二次型拟合的 m_d=%.4f（旧方向分解值 %s）" % (carrier, bd["m_d_m0"], m))
        mprov = "三维二次型拟合 m_d=(det M)^(1/3)（%s）；旧的方向分解值 %s 不再用（3D 各向异性时偏轻）；%s" % (
            bd["m_provenance"], m, mprov)
        m = bd["m_d_m0"]
    e1, e1prov = get_E1(cwd, carrier)
    rec = {"carrier": carrier,
           "inputs": {"C_%s" % ("2D_N_per_m" if is_2d else "3D_Pa"): C,
                      "C_provenance": Cprov,
                      "m_eff_m0": m, "m_provenance": mprov,
                      "E1_eV": e1, "E1_provenance": e1prov,
                      "T_K": T}}
    if None in (C, m, e1):
        miss = [n for n, v in (("C", C), ("m*", m), ("E1", e1)) if v is None]
        rec["mobility_cm2_Vs"] = None
        rec["status"] = ("缺输入：%s —— 在本材料 step.conf 手填（M_EFF_*/E1_*_EV/C_*）后重跑（μ∝1/T，报 %.0fK）"
                         % ("、".join(miss), T))
    else:
        rec["mobility_cm2_Vs"] = round(mobility_dpt(is_2d, C, m, e1, T), 3)
        rec["status"] = "ok"
    rec["by_direction"] = bd
    rec["multi_valley"] = _multi_valley_block(cwd, is_2d, carrier, T, temps or [T],     # [0022]
                                              rec["by_direction"])
    scope = _single_valley_scope(rec["multi_valley"], T)                            # [0024]
    if scope:
        rec["single_valley_scope"] = scope
        print("%s %s：%s" % ("[WARN]" if scope["share"] < 0.5 else "[..]", carrier, scope["note"]))
    return rec


def _za_check(cwd):
    """[0016] σh / ZA 判断（ke_common.za_coupling_check），写进 dpt_result.json；没有 σh 打 WARN。
    DPT 只算 LA 支：有 σh 时 ZA 一阶耦合被对称性禁止，没有 σh 时结果缺这一项（文献 DPT 同样不含）。"""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        import ke_common as _kc
    except ImportError:
        print("[WARN] 没有 ke_common.py，σh/ZA 没判（skill.yaml 的 gen_need 漏了它？）")
        return None
    r = _kc.za_coupling_check(cwd)
    print(("[WARN] " if r.get("sigma_h") is False else "[..] ") + "ZA：" + r["note"])
    return r


def main():
    cwd = Path.cwd()
    _guard_not_0d(cwd)
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    invalidate_consumers(cwd)
    overrides = apply_conf(cwd)
    dim = _read_dim(cwd)
    is_2d = (dim == "2d")
    s3_grid = _s3_grid_gate(cwd, dim)                     # [V174] 在算 m* 之前
    temps = sorted(set(parse_temps(TEMPERATURES)) | {round(TEMPERATURE_K, 6)})   # [0020]
    functional = _functional(cwd)                         # [V181]

    carriers = (["electron", "hole"] if CARRIER == "both" else [CARRIER])
    # [ENV] 记录执行环境，铺开时"哪个环境跑的"从推断变成可读事实
    try:
        import socket as _socket
        _amset_ver = "N/A"
        try:
            import amset as _amset
            _amset_ver = getattr(_amset, "__version__", "N/A")
        except Exception:
            pass
        _env = {"python": sys.executable, "amset": _amset_ver,
                "host": _socket.gethostname()}
    except Exception:
        _env = {"python": sys.executable}
    res = {"dim": dim, "is_2d": is_2d, "temperature_K": TEMPERATURE_K,
           "skill_rev": _SKILL_REV, "env": _env, "overrides": overrides, "s3_grid": s3_grid,
           "functional": functional,
           "formula": ("2D Bardeen-Shockley: μ=eℏ³C_2D/(k_BT m* m_d E1²)" if is_2d
                       else "3D: μ=2√(2π)eℏ⁴C_3D/(3(k_BT)^{3/2}m*^{5/2}E1²)"),
           "results": [_one_carrier(cwd, is_2d, c, TEMPERATURE_K, temps) for c in carriers],
           "note": ("经典 DPT 声学支迁移率，用于和 amset(ADP) 对标。m*/E1 自动值精度"
                    "有限——务必核对，或在本材料 step.conf 手填。μ∝1/T。")}
    res["temperatures_K"] = temps                          # [0020] μ(T) 网格（含 temperature_K）
    for r in res["results"]:
        r["mobility_vs_T"] = mobility_vs_T(r, is_2d, TEMPERATURE_K, temps)
    if is_2d:
        res["za_coupling"] = _za_check(cwd)                # [0016] σh / ZA
    (out / "dpt_result.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = ["# DPT 迁移率摘要（%s, %.0f K）" % (dim, TEMPERATURE_K)]
    if res.get("za_coupling"):                            # [0016]
        lines.append("# ZA：%s" % res["za_coupling"]["note"])
    if functional:                                        # [V181]
        lines.append("# 泛函：%s —— %s" % (functional["label"], functional["note"]))
    if s3_grid and s3_grid["issues"]:                     # [V174] 只在 ALLOW_COARSE_S3 时走到这里
        lines.append("# ★ S3 网格 %s 不合现行规则（ALLOW_COARSE_S3）：m* 只作复现/对照 —— %s"
                     % ("×".join(map(str, s3_grid["mesh"] or [])), "；".join(s3_grid["issues"])))
    for r in res["results"]:
        mu = r["mobility_cm2_Vs"]
        lines.append("%-9s μ = %s cm²/(V·s)   [%s]" % (
            r["carrier"], (("%.1f" % mu) if mu is not None else "—"), r["status"]))
        if r.get("single_valley_scope"):                        # [0024]
            lines.append("           %s [单谷代表] %s" % ("★" if r["single_valley_scope"]["share"] < 0.5 else " ",
                                                         r["single_valley_scope"]["note"]))
        i = r["inputs"]
        lines.append("           C=%s(%s)  m*=%s(%s)  E1=%s eV(%s)" % (
            i.get("C_2D_N_per_m", i.get("C_3D_Pa")), i["C_provenance"],
            i["m_eff_m0"], i["m_provenance"], i["E1_eV"], i["E1_provenance"]))
        bd = r.get("by_direction") or {}                      # patch_dpt_aniso
        if bd.get("status") == "ok":
            for d in [d for d in ("x", "y", "z") if isinstance(bd.get(d), dict)]:   # [0025] 3D 有 z
                q = bd[d]
                lines.append("           [%s] mu=%.1f cm2/Vs  tau=%.1f fs  "
                             "m*=%.3f  E1=%.3f eV  %s" % (
                                 d, q["mobility_cm2_Vs"], q["tau_s"] * 1e15, q["m_eff_m0"], q["E1_eV"],
                                 ("C_2D=%.1f N/m" % q["C_2D_N_per_m"]) if "C_2D_N_per_m" in q
                                 else ("C=%.1f GPa" % q["C_3D_GPa"])))
            lines.append("           m_d=%.3f  [%s | %s]" % (
                bd["m_d_m0"], bd["m_provenance"], bd["E1_provenance"]))
        elif bd.get("status"):
            lines.append("           [各向异性] %s" % bd["status"])
        mv = r.get("multi_valley") or {}                      # [0022]
        if mv.get("status") == "ok":
            lines.append("           [多谷] mu_x=%.1f mu_y=%.1f cm2/Vs  tau_x=%.1f tau_y=%.1f fs  m_c=%.3f/%.3f  "
                         "（%.0f K，窗口 %.2f eV，无谷间散射 -> 上限）" % (
                             mv["x"]["mobility_cm2_Vs"], mv["y"]["mobility_cm2_Vs"],
                             mv["x"]["tau_s"] * 1e15, mv["y"]["tau_s"] * 1e15,
                             mv["x"]["m_c_m0"], mv["y"]["m_c_m0"], mv["T_K"], mv["window_eV"]))
            for vv in mv["valleys"]:
                lines.append("             %-8s ×%d  ΔE=%.3f eV  m_d=%.3f  E1(x,y)=%s eV  布居 %.1f%%" % (
                    vv["label"], vv["g"], vv["dE_eV"], vv["m_d_m0"],
                    "/".join("%.2f" % x for x in vv["E1_xy_eV"][0]), 100 * vv["weight_T0"]))
            cons = mv.get("edge_star_vs_by_direction")
            if cons:
                lines.append("             只留 %s 星 / by_direction = %.3f (x), %.3f (y)（应≈1）"
                             % (cons["star"], cons["x"], cons["y"]))
            for w_ in mv.get("warnings") or []:
                lines.append("             ★ %s" % w_)
        elif mv.get("status") and not str(mv["status"]).startswith(("off", "多谷 DPT 只实现 2D")):
            lines.append("           [多谷] %s" % mv["status"])
            for vv in mv.get("valleys") or []:                 # [0023] 拒绝时也列出各谷，看得出是哪个谷、为什么
                lines.append("             %-8s ×%d  ΔE=%.3f eV  m_d=%.3f  E1_eff=%.2f eV  最大布居 %.1f%%  "
                             "带符号 dE/dε(x,y)=%s eV" % (
                                 vv["label"], vv["g"], vv["dE_eV"], vv["m_d_m0"], vv.get("E1_eff_eV", float("nan")),
                                 100 * vv.get("pop_max", 0.0),
                                 "/".join("%+.2f" % x for x in (vv.get("E1_signed_xy_eV") or [[0, 0]])[0])))
    # [0020] μ(T) 表：2D 给 x / y / 面内平均，3D 给主值
    lines.append("")
    lines.append("# μ(T) cm²/(V·s)，由 %.0f K 换算（%s）" % (TEMPERATURE_K, "∝1/T" if is_2d else "∝T^-1.5"))
    for r in res["results"]:
        tab = r.get("mobility_vs_T") or []
        if not tab or all(t["mobility_cm2_Vs"] is None and "mobility_x_cm2_Vs" not in t for t in tab):
            continue
        has_xy = "mobility_x_cm2_Vs" in tab[0]
        has_z = "mobility_z_cm2_Vs" in tab[0]                        # [0025] 3D
        lines.append("%-9s %7s %10s%s" % (r["carrier"], "T(K)", "μ",
                                         (("%10s %10s %10s %10s" % ("μx", "μy", "μz", "(μx+μy+μz)/3")) if has_z
                                          else ("%10s %10s %10s" % ("μx", "μy", "(μx+μy)/2")) if has_xy else "")))
        for t in tab:
            mu = t["mobility_cm2_Vs"]
            lines.append("%-9s %7.0f %10s%s" % (
                "", t["T_K"], ("%.1f" % mu) if mu is not None else "—",
                ("%10.1f %10.1f %10.1f %10.1f" % (t["mobility_x_cm2_Vs"], t["mobility_y_cm2_Vs"],
                                                  t["mobility_z_cm2_Vs"], t["mobility_mean_cm2_Vs"])) if has_z
                else ("%10.1f %10.1f %10.1f" % (t["mobility_x_cm2_Vs"], t["mobility_y_cm2_Vs"],
                                                t["mobility_inplane_cm2_Vs"])) if has_xy else ""))
    mv_ok = [r for r in res["results"] if (r.get("multi_valley") or {}).get("status") == "ok"]
    if mv_ok:                                              # [0022] 多谷 μ(T)：每个温度单独算，不是 ∝1/T
        lines.append("")
        lines.append("# 多谷 μ(T) cm²/(V·s)：每个温度单独算（谷间布居随 T 变）；不含谷间散射，是上限")
        for r in mv_ok:
            tab = r["multi_valley"]["mobility_vs_T"]
            labs = list(tab[0]["weights"]) if tab else []
            lines.append("%-9s %7s %10s %10s %10s  %s" % (r["carrier"], "T(K)", "μx", "μy", "(μx+μy)/2",
                                                         "  ".join("%s%%" % l for l in labs)))
            for t in tab:
                lines.append("%-9s %7.0f %10.1f %10.1f %10.1f  %s" % (
                    "", t["T_K"], t["mobility_x_cm2_Vs"], t["mobility_y_cm2_Vs"], t["mobility_inplane_cm2_Vs"],
                    "  ".join("%.1f" % (100 * t["weights"][l]) for l in labs)))
    (out / "dpt_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 软失败 → 报错：关键产物(迁移率)一个都算不出，就让 tf 判 error 而非 completed，
    # 并指明缺什么、该跑哪一步（不手填、全自动的前提下）。
    if all(r["mobility_cm2_Vs"] is None and (r.get("multi_valley") or {}).get("status") != "ok"
           for r in res["results"]):
        miss = set()
        for r in res["results"]:
            i = r["inputs"]
            # [V143] 带上具体原因：只写"m*←S3_uniform(vasprun/网格)"时，WS2（S3 网格 47，K 离网）被误读成
            #   "vasprun 缺失"，白查了一轮。原因就在 *_provenance 里，直接打出来。
            if i.get("C_2D_N_per_m") is None and i.get("C_3D_Pa") is None:
                miss.add("C←S6_elastic（%s）" % i.get("C_provenance"))
            if i["m_eff_m0"] is None:
                miss.add("m*←S3_uniform（%s）" % i.get("m_provenance"))
            if i["E1_eV"] is None:
                miss.add("E1←S7.1_read（%s）" % i.get("E1_provenance"))
        hint = "；".join(sorted(miss)) or "见 json 的 provenance"
        sys.exit("[ERROR] DPT 迁移率全部未算出。缺：%s。请先把对应步骤跑好/修好再重跑本步"
                 "（dpt_result.json 已写出，含各输入来源可排查）。" % hint)

    for r in res["results"]:
        if r["mobility_cm2_Vs"] is None and (r.get("multi_valley") or {}).get("status") == "ok":
            print("[OK] %s：单谷 DPT 未出数（%s），多谷 DPT 已算出（dpt_result.json 的 multi_valley）"
                  % (r["carrier"], (r["inputs"].get("m_provenance") or "")[:60]))
    print("[DONE] %s：dpt_result.json 已生成" % OUTDIR_NAME)


if __name__ == "__main__":
    main()