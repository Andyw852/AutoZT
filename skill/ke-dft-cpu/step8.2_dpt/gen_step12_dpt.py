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
_SKILL_REV = "2026-10-09-v0020"
# [R-SCAN] 强制选点壳层 NSTEP（None=自动 2/3/4/5；填 2/3/4/5=只试该值，做 R-scan 用）。
# 用法：tf -p <材料> -j step8.2_dpt conf --set FORCE_NSTEP=3 → retry + start，对比 dpt_result.json 的 m_d
#   与 m_provenance 里的 R（V148 起写 step.conf，不改本脚本）。
FORCE_NSTEP = None
OUTDIR_NAME = "step8.2_dpt"
UNIFORM_DIR = "step3_uniform"
ELASTIC_DIR = "step6_elastic"
DEFORM_READ_DIR = "step7b_deform_read"
AMSET_DIR   = "step8_amset"          # 读 2d_correction.json 拿层厚 t、c/t

TEMPERATURE_K = 300.0                # DPT 迁移率报此温度（μ∝1/T，可换算）
# [0020] 多温度：mobility_vs_T 按 μ∝1/T（2D）/ T^-1.5（3D）从 TEMPERATURE_K 的值精确换算（m*/E1/C 与 T 无关）。
#   默认与 S8/S8.4 的 AMSET 网格（100:900:9）、S8.1 的 TEMPERATURES 一致，S8.3 才能逐温度对齐；
#   只要 100–700 K 就在本步 step.conf 写 TEMPERATURES = 100:700:7。格式 "起:止:个数" 或逗号列表。
TEMPERATURES = "100:900:9"
CARRIER = "both"                     # electron / hole / both
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
    global E1_SOURCE, FORCE_NSTEP, ALLOW_COARSE_S3, TEMPERATURE_K, TEMPERATURES
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
            for d in ("x", "y"):
                row["mobility_%s_cm2_Vs" % d] = round(bd[d]["mobility_cm2_Vs"] * f, 3)
                row["tau_%s_s" % d] = bd[d]["tau_s"] * f
            row["mobility_inplane_cm2_Vs"] = round(
                (row["mobility_x_cm2_Vs"] + row["mobility_y_cm2_Vs"]) / 2.0, 3)
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
                "没用 —— 用 amset eff-mass 的质量写本材料 step.conf（M_EFF_*）" % np.round(_kf, 4)), False
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
        c11, c22 = read_elastic_inplane_GPa(cwd)
        if c11 is None:
            return None, "Pa", "缺 step6 弹性"
        return (c11 + c22) / 2.0 * 1e9, "Pa", "step6 (C11+C22)/2"


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

        # [P1-2] 展开 IBZ → 全 BZ：带边 K 点坐在 IBZ 楔形边界，近邻点只存在于
        # 楔形内一侧，单侧取点会让 a/b 两个曲率被不对称地吃掉（hex 实测 2.7 倍差）。
        # 展开后 k0 周围点对称，二次型拟合才无偏。失败退回 IBZ（靠 cond/残差验收兜底）。
        # [TR] 时间反演只在净磁矩≈0 时是能带的对称操作。CrS2/CrSe2 塌到零磁矩
        # 无碍，但本 skill 还要复用到 MnIn2Se4/Mn2In2Se5 等磁性体系——自旋分辨
        # 能带在有净磁矩时 TR 会把 up/down 错配。有磁矩就关掉 TR。
        # [TR] 时间反演只在净磁矩≈0 时是能带的对称操作。ISPIN=1 → TR 安全；
        # ISPIN=2 → 读 OUTCAR 总磁矩判断，读不到 → 保守关掉（磁性体系不安全）。
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
                    _seen, _kf_list, _idx_list = {}, [], []
                    for _i, _k in enumerate(kfrac_ibz):
                        _targets = [op.operate(_k) for op in _ops]
                        if use_tr:
                            _targets += [op.operate(-_k) for op in _ops]
                        for _kk in _targets:
                            _kk -= np.round(_kk)
                            _key = tuple(np.round(_kk, 6))
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
                    return None, ("全 BZ 展开失败（amset/subprocess/pymatgen 都不可用：%s）"
                                  "——拒绝在 IBZ 上拟合 m*（单侧取点是系统性错误）"
                                  % _exc_brief(_e2))
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

        kcart = kfrac @ recip
        dk = kcart - kcart[k0]
        de = ene[:, b0] - ene[k0, b0]
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
        need_lin = not _is_high_sym(kfrac[k0])

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
    """分方向 DPT 结果块。3D 不走这条分支。"""
    if not is_2d:
        return {"status": "3D 体系不做各向异性分支，用上面的各向同性值"}
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


def _one_carrier(cwd, is_2d, carrier, T):
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
    rec["by_direction"] = _aniso_block(cwd, is_2d, carrier, T)   # patch_dpt_aniso
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
           "results": [_one_carrier(cwd, is_2d, c, TEMPERATURE_K) for c in carriers],
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
        i = r["inputs"]
        lines.append("           C=%s(%s)  m*=%s(%s)  E1=%s eV(%s)" % (
            i.get("C_2D_N_per_m", i.get("C_3D_Pa")), i["C_provenance"],
            i["m_eff_m0"], i["m_provenance"], i["E1_eV"], i["E1_provenance"]))
        bd = r.get("by_direction") or {}                      # patch_dpt_aniso
        if bd.get("status") == "ok":
            for d in ("x", "y"):
                q = bd[d]
                lines.append("           [%s] mu=%.1f cm2/Vs  tau=%.1f fs  "
                             "m*=%.3f  E1=%.3f eV  C_2D=%.1f N/m" % (
                                 d, q["mobility_cm2_Vs"], q["tau_s"] * 1e15,
                                 q["m_eff_m0"], q["E1_eV"], q["C_2D_N_per_m"]))
            lines.append("           m_d=%.3f  [%s | %s]" % (
                bd["m_d_m0"], bd["m_provenance"], bd["E1_provenance"]))
        elif bd.get("status"):
            lines.append("           [各向异性] %s" % bd["status"])
    # [0020] μ(T) 表：2D 给 x / y / 面内平均，3D 给主值
    lines.append("")
    lines.append("# μ(T) cm²/(V·s)，由 %.0f K 换算（%s）" % (TEMPERATURE_K, "∝1/T" if is_2d else "∝T^-1.5"))
    for r in res["results"]:
        tab = r.get("mobility_vs_T") or []
        if not tab or all(t["mobility_cm2_Vs"] is None and "mobility_x_cm2_Vs" not in t for t in tab):
            continue
        has_xy = "mobility_x_cm2_Vs" in tab[0]
        lines.append("%-9s %7s %10s%s" % (r["carrier"], "T(K)", "μ",
                                         ("%10s %10s %10s" % ("μx", "μy", "(μx+μy)/2")) if has_xy else ""))
        for t in tab:
            mu = t["mobility_cm2_Vs"]
            lines.append("%-9s %7.0f %10s%s" % (
                "", t["T_K"], ("%.1f" % mu) if mu is not None else "—",
                ("%10.1f %10.1f %10.1f" % (t["mobility_x_cm2_Vs"], t["mobility_y_cm2_Vs"],
                                           t["mobility_inplane_cm2_Vs"])) if has_xy else ""))
    (out / "dpt_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 软失败 → 报错：关键产物(迁移率)一个都算不出，就让 tf 判 error 而非 completed，
    # 并指明缺什么、该跑哪一步（不手填、全自动的前提下）。
    if all(r["mobility_cm2_Vs"] is None for r in res["results"]):
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

    print("[DONE] %s：dpt_result.json 已生成" % OUTDIR_NAME)


if __name__ == "__main__":
    main()