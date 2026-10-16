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
  没有就用 transport.json 的 ADP 单机制）、AMSET 去 IMP（intrinsic_transport.json）、AMSET 全机制（overall）。
[0032] 另出本征五参数图 + 输运参数表（2D，尽力而为，不影响完成标记）：intrinsic/ 下
  intrinsic_transport_<T>K.png（每个 AMSET 温度一张：n/p × S、σ、κe、Lorenz、μ；包含可用的 ADP、ADP+POP、
  ADP+POP+ODP 及去 IMP/overall 结果，并含 BT2×DPT τ / DPT 参照）、intrinsic_transport_table.csv（全部温度 × 原生掺杂档 × 数据集）、intrinsic_summary.txt。
  2D 图双横轴：下 n_2D（cm^-2），上文献口径 n_2D/t（cm^-3）；comparison_300K.png 同样处理。
[0053] transport_plots/intrinsic_carrier/<M>_mechanism_resolved_mobility.png/.csv（[0054] 前叫 dilute_mobility）：
  只取稀薄极限这一个浓度，phonon-limited μ(T) 按散射机制拆开（单机制 μ、总值、Matthiessen 倒数和、限制机制占比）；
  空穴取 n = p 本征点，电子取最低 n 型点。
[0054] 命名与电声文献统一：去 IMP = phonon-limited（文献 phonon-limited / intrinsic mobility），含 IMP = phonon + impurity，
  "intrinsic" 只指本征载流子（n = p = n_i）。<M>_intrinsic.png -> <M>_phonon_limited.png。内部数据文件名
  （intrinsic_transport.json、step8.3_output/intrinsic/、数据集键 intrinsic/overall）不变，保证旧结果与 zt-dft-cpu 可读。
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
_SKILL_REV = "2026-10-16-0054"   # [0054] 与文献统一命名：去 IMP = phonon-limited（<M>_phonon_limited.png，原 <M>_intrinsic.png）、含 IMP = phonon + impurity、intrinsic 只指本征载流子；<M>_mechanism_resolved_mobility（原 dilute_mobility）；本征点按去 IMP 后机制去重；η=(E_c−E_F)/kT 判 f≈0   # [0053] intrinsic_carrier/<M>_dilute_mobility.png：稀薄极限（f≈0）μ(T) 按散射机制拆开 + Matthiessen 对照（空穴取 n = p 本征点，电子取最低 n 型去 IMP）   # [0052] transport_plots/ 每次整目录重建（先删旧结果）；intrinsic_carrier/ 标注本征浓度 n_i(T)（undoped_ni.json / fermi_fix 日志）；<M>_intrinsic.png 标题写明 intrinsic = 本征散射、横轴 = 掺杂浓度   # [0047] <M>_intrinsic.png 只画基础 4 条本征曲线（此前与 <M>_transport.png 同图仅换标题）+ intrinsic CSV 同步只含 4 条   # [0046] 本征五参数套件支持 3D（读 S8 的 step8_amset + subset_*.json；_geom_3d 体胞几何、c/t=1、对角 (xx+yy+zz)/3）；3D 也出 transport_plots/<M>_transport|<M>_intrinsic.png 9x5 面板   # [0045] bt2_muscan 同时画 thickness(我们)/cell(文献) 两种口径；[0041] transport_plots/defect_equilibrium/（缺陷平衡净浓度 σ/S/PF/κe vs T）   # 2026-10-13-0040   # [0040] transport_plots/intrinsic_carrier/（未掺杂本征载流子 σ/S/PF/κe vs T）；图例/文件名一律用实际机制组合 amset_<机制>；2026-10-13-0039   # [0039] transport_plots/datasets/：每个数据集单独一张 σ/S/PF/κe 图 + csv； [0038] 本征/子集图读 subset_*.json（SUBSETS）；[0036] transport_plots 加 bt2_muscan/ + amset_scan/ 两个扫描图子文件夹；[0035] transport_plots 合并 + τ 对比（tau/）
#   bump 记录：2026-10-11-0035 —— transport_plots 合并：intrinsic 图并入 transport 图（4 条数据集
#   overall / intrinsic / ADP only / BT2xDPT-tau），不再单独出 <M>_intrinsic*.{png,csv}，
#   只留 <M>_transport.png + 每温度 CSV + <M>_params_<T>K.csv。
#   bump 记录：2026-10-11-0035（Claude）—— 新增 τ 对比（tau/）：τ(E) / τ_eff vs DPT τ。
#   bump 记录：2026-10-11-0034-r2 —— μ(T) 用第一份通过网格/同轮核对的 overall_no*.json（以前只看排序后的第一份）。
#   bump 记录：2026-10-11-0034 —— comparison_300K 的可解读判据改每原子（与 S8.4、本征表同一判据，缺原子数才退回每原胞），
#   csv 加 carriers_per_atom；本征表/图、μ(T) 加 S8.4 的 overall_no*.json（overall 去掉 OVERALL_WITHOUT，出厂去 ODP）。
#   bump 记录：2026-10-10-0032-r2 —— 复审修正：本征文件按 transport 指纹判同一轮；可解读判据改每原子（与 S8.4 同）；
#   BT2 μ 双极区屏蔽；DPT 多谷标注；图例用实际机制；0 掺杂跳过；intrinsic/ 先清旧产物；3D 对比图与改动前逐字节相同。
#   bump 记录：2026-10-10-0032 —— 本征五参数图 + 输运参数表（intrinsic/，全部温度，2D）；2D 对比图双横轴
#   （下 n_2D cm^-2、上文献口径 n_2D/t cm^-3）+ 浓度口径图注；Lorenz 面板以 1e-8 为单位；load_bt2/load_dpt 可传 T。
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
NCELL_INTERPRET_MAX = 0.1         # 每原胞载流子数（只在数不出原子数时兜底）
# ★ [0034] 2026-10-11 用户决定：comparison_300K 也按【每原子】载流子数判可解读（文献通行口径，与 S8.4 /
#   本征表同一判据）。每原胞数取决于晶胞大小（原胞 vs 惯用胞差几倍），按原子才与晶胞选择无关。
#   阈值：2D 取 S8.4 2d_correction.json 的 carrier_per_atom_threshold，3D 与缺键时用这个默认值。
CARRIER_PER_ATOM_MAX = 0.02
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


def _cell_atoms(cwd):
    """[0034] 晶胞原子数：与 _cell_geometry 同一批候选结构文件（VASP5 / VASP4 都认）；读不到 -> None。"""
    import glob as _glob
    cands = []
    for _d in ("step3_uniform", "step1_opt", "step1_std_opt", "."):
        cands += sorted(_glob.glob(str(Path(cwd) / _d / "CONTCAR")))
        cands += sorted(_glob.glob(str(Path(cwd) / _d / "POSCAR")))
    for p in cands:
        try:
            ln = Path(p).read_text(errors="ignore").splitlines()
            for row in (ln[5], ln[6]):
                toks = row.split()
                if toks and all(t.isdigit() for t in toks):
                    n = sum(int(t) for t in toks)
                    if n > 0:
                        return n
        except (OSError, IndexError, ValueError):
            continue
    return None


def _atom_rule(cwd):
    """[0034] 可解读判据：(原子数, 每原子阈值, 来源说明)。2D 优先读 S8.4 的 2d_correction.json（与本征表、S8.4 同源），
    否则数结构文件的原子、阈值用 CARRIER_PER_ATOM_MAX；原子数也读不到 -> (None, None, 原因)。"""
    for run in (AMSET2D_DIR, AMSET_DIR):
        rec = _load_json(Path(cwd) / run / "2d_correction.json") or {}
        n, thr = _num(rec.get("num_atoms")), _num(rec.get("carrier_per_atom_threshold"))
        if n and n > 0:
            return int(n), (thr if thr and thr > 0 else CARRIER_PER_ATOM_MAX), "%s/2d_correction.json" % run
    n = _cell_atoms(cwd)
    if n:
        return n, CARRIER_PER_ATOM_MAX, "结构文件（step3_uniform 等）"
    return None, None, "数不出原子数（step3_uniform/step1 的 POSCAR/CONTCAR 都读不到）"


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
    if not target_signed or target_signed != target_signed:   # [0032] 0 / NaN：没有同号分支
        return None
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


def load_bt2(cwd, T=None):
    """[0032] T=None 即 TARGET_T（原行为）；本征五参数图按 AMSET 的每个温度调用。"""
    j = _load_json(Path(cwd) / BT2_DIR / "boltztrap_crta.json")
    if not j:
        return None
    try:
        temps = j["temperatures_K"]
        iT = _nearest_T_index(temps, TARGET_T if T is None else T)
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


def load_dpt(cwd, T=None, quiet=False):
    """[0032] T=None 即 TARGET_T（原行为）；quiet 只关换温提示（逐温度调用时不刷屏）。"""
    j = _load_json(Path(cwd) / DPT_DIR / "dpt_result.json")
    if not j:
        return None
    T = TARGET_T if T is None else float(T)
    # [0021] 对比在 TARGET_T：S8.2 的 TEMPERATURE_K 不是它时，μ 与 τ 一并换温（以前直接拿来比）
    f, T0 = _dpt_T_factor(j, T)
    out = {"T_used": T if abs(f - 1.0) > 1e-12 else T0, "T0": T0}
    if abs(f - 1.0) > 1e-12 and not quiet:
        print("[..] DPT 结果在 %.0f K，换算到 %.0f K（×%.4f）再对比" % (T0, T, f))
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
        elif mv.get("status") == "ok" and not quiet:
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
    _natom, _thr_atom, _atom_src = _atom_rule(Path.cwd())          # [0034] 每原子判据
    if _add_dose:
        _dcols = (["n_2D_cm-2"] if (_c_cm and (am or am2)) else [])
        cols = _dcols + ["carriers_per_cell", "carriers_per_atom", "in_interpret_range"] + cols
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
                        if _natom:                                     # [0034] 每原子（与 S8.4 同判据）
                            _nat = _ncell / _natom
                            _row["carriers_per_atom"] = "%.3g" % _nat
                            _row["in_interpret_range"] = (
                                "yes" if _nat <= _thr_atom else
                                "NO(>%.3g/atom：超出刚性能带近似/实验范围)" % _thr_atom)
                        else:                                          # 数不出原子：退回每原胞
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
    # [V187] 维度按 step1 的 DIM 判：2D 默认不跑 S8（amset3d 关），am 为空时旧写法会把 2D 标成 "[3D]"
    _red = ("面内(xx+yy)/2 [2D]" if ((am and am.get("is_2d")) or _read_dim(Path.cwd()) == "2d")
            else "对角(xx+yy+zz)/3 [3D]")
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
        # [V187] 2D 的 "ionrelax(in-plane)+clamped(out-of-plane)" 是正常口径：面外应变拉的是真空，二维核（S8.4）
        #   只用面内分量；只有 S8（三维核）会用到面外那几项，那一栏另注。
        _ok2d = _read_dim(Path.cwd()) == "2d" and _dgeom.startswith("ionrelax(in-plane)")
        _warn = ("" if _dgeom == "ionrelax" else
                 ("  （面外分量离子固定：二维核不用；只影响 S8 三维核的对照值）" if _ok2d else
                  "  ⚠ 非离子弛豫口径：与 TOTAL ELASTIC MODULI 不同源，ADP 绝对值跨项目不可直接比"))
        lines.append("# [口径] 形变势构型 deform_geometry = %s（ionrelax=离子弛豫；clamped=离子固定）%s"
                     % (_dgeom, _warn))
    if _add_dose:
        _rule = (("每原子 %.3g 个载流子（%d 原子/胞，来源 %s；与 S8.4、本征表同一判据）" % (_thr_atom, _natom, _atom_src))
                 if _natom else ("每原胞 %.2g 个载流子（%s，退回旧判据）" % (NCELL_INTERPRET_MAX, _atom_src)))
        lines.append("# [范围] 可解读上限 = %s"
                     "（n_cell = n_2D×A_cell 或 n_3D×V_cell；本胞 A = %s Å², V = %s Å³）%s。"
                     "超出档位在 csv 里标 in_interpret_range=NO：刚性能带近似已不成立、"
                     "实验也达不到，只看趋势，不进论文图。"
                     % (_rule,
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
MU_T_COLS = ("mu_DPT", "mu_DPT_mv", "mu_AMSET_ADP", "mu_AMSET_intrinsic", "mu_AMSET_all",   # [0022] +mu_DPT_mv
             "mu_AMSET_all_without")   # [0034] overall 去掉 OVERALL_WITHOUT 的机制（S8.4 的 overall_no*.json，出厂去 ODP）
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
        if not _intr_same_run(Path(cwd) / run, q)[0]:                  # [0032] 上一轮残留
            print("[WARN] %s/%s 是上一轮 S8.4 留下的（与当前 transport.json 内容不同）——不用" % (run, name))
            return None, None
        return q["intrinsic"]["mobility_cm2_Vs_s"].get("overall"), q

    _adp_file = _intrinsic_file(Path(cwd) / run, "ADP_only")
    adp, q_adp = _intr(_adp_file) if _adp_file else (None, None)
    if adp is not None:
        cols["mu_AMSET_ADP"] = adp
        src["mu_AMSET_ADP"] = "%s/%s（只留 %s 重积分）" % (
            run, _adp_file, "+".join(q_adp.get("intrinsic_labels") or ["ADP"]))
    else:
        a = _unwrap(mob.get("ADP") or mob.get("ACD"))
        if a is not None:
            cols["mu_AMSET_ADP"] = a
            src["mu_AMSET_ADP"] = "%s/transport.json mobility.ADP（单机制）" % run
    _intr_file = _intrinsic_file(Path(cwd) / run, "intrinsic")
    intr, q_in = _intr(_intr_file) if _intr_file else (None, None)
    labels = None
    if intr is not None:
        labels = list(q_in.get("intrinsic_labels") or [])
        cols["mu_AMSET_intrinsic"] = intr
        src["mu_AMSET_intrinsic"] = "%s/%s（%s，去 %s）" % (
            run, _intr_file, "+".join(labels) or "?", "+".join(q_in.get("intrinsic_dropped") or []) or "?")
    without = None
    for name in _overall_without_files(Path(cwd) / run):              # [0034] 第一份通过网格/同轮核对的
        wo, q_wo = _intr(name)
        if wo is not None:
            without = list(q_wo.get("intrinsic_dropped") or [])
            cols["mu_AMSET_all_without"] = wo
            src["mu_AMSET_all_without"] = "%s/%s（overall 去 %s：%s）" % (
                run, name, "+".join(without) or "?", "+".join(q_wo.get("intrinsic_labels") or []) or "?")
            break
    return {"run": run, "doping": doping, "temps": temps, "cols": cols, "src": src,
            "intrinsic_labels": labels, "without": without}


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
            "all_without": am.get("without") if am else None,             # [0034]
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
             "# %-8s %4s %5s %8s %8s %7s %7s %8s %9s %9s %9s %9s %9s %9s" % (
                 "carrier", "dir", "T(K)", "E1", "C2D", "m*", "m_d", "tau",
                 "DPT", "DPT_mv", "AMSET_ADP", "ph_limit", "ph+imp",
                 ("w/o" + "+".join(res["meta"].get("all_without") or []))[:9] if res["meta"].get("all_without")
                 else "all_w/o")]
    for r in res["rows"]:
        lines.append("  %-8s %4s %5.0f %8s %8s %7s %7s %8s %9s %9s %9s %9s %9s %9s" % (
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
             "mu_AMSET_intrinsic": ("C2", "-.", "AMSET phonon-limited"), "mu_AMSET_all": ("C3", ":", "AMSET phonon + impurity"),
             "mu_AMSET_all_without": ("C4", "--", "AMSET phonon + impurity w/o %s"         # [0034]
                                      % "+".join(res["meta"].get("all_without") or ["?"]))}
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


# =====================================================================
# [0032] 本征五参数图 + 输运参数表（2026-10-10 用户要求）
#   数据：S8.4 同一次运行的结果集 —— intrinsic_ADP.json、可用时的 intrinsic_ADP_POP*.json、intrinsic_transport.json
#   （去 IMP，保留其余已启用机制）和 transport.json（全部已启用机制）；
#   外加 BT2 验证线（CRTA 的 σ/τ、κe/τ × DPT 的 τ）与 DPT 单谷迁移率（水平线）。
#   输出到 step8.3_output/intrinsic/：每个 AMSET 温度一张 2×5 图（行 n 型 / p 型；列 S、σ、κe、Lorenz、μ），
#   一张长表（全部温度 × 原生掺杂档 × 数据集）和一份文字汇总。只做 2D（S8.4 只跑 2D）。
#   2D 横轴：下 = 面密度 n_2D（cm^-2），上 = 文献口径 n_2D/t（cm^-3，t = 层厚）；AutoZT/AMSET 文件里的
#   cm^-3 用的是含真空的胞高 c（n_cell = n_2D/c），两者差 c/t 倍。σ、κe、PF 按层厚 t 归一（= 元胞值 × c/t），
#   μ、S、Lorenz 与口径无关。尽力而为，不影响完成标记（comparison_300K.png）。
# =====================================================================
INTR_SUBDIR = "intrinsic"
INTR_CSV = "intrinsic_transport_table.csv"
INTR_TXT = "intrinsic_summary.txt"
INTR_SETS = (   # (键, 文件, 图例, 线型, 颜色, 标记)
    ("ADP_only", "intrinsic_ADP.json", "ADP only (acoustic limit)", ":", "C2", "^"),
    ("ADP_POP", "intrinsic_ADP_POP.json", "ADP+POP", "-", "C1", "v"),
    ("ADP_POP_ODP", "intrinsic_ADP_POP_ODP.json", "ADP+POP+ODP", "-", "C6", "P"),
    ("intrinsic", "intrinsic_transport.json", "phonon-limited", "-", "C0", "o"),
    ("overall", "transport.json", "phonon + impurity", "-.", "C3", "s"),
)
# [0034] S8.4 的 OVERALL_WITHOUT 那一套（overall_no*.json，出厂 overall_noODP.json）：有就自动加进图、表、汇总。
INTR_WITHOUT_STYLE = (("--", "C4", "D"), ("--", "C5", "v"))
BT2_KEY, BT2_LABEL = "BT2_DPTtau", "BT2×DPT-τ"


def _overall_without_files(run_dir):
    """[0034] S8.4 目录里的 overall_no*.json（按名字排序）。"""
    import glob as _glob
    return [Path(p).name for p in sorted(_glob.glob(str(Path(run_dir) / "overall_no*.json")))]


# [0038] subset_*.json（SUBSETS 生成）的键分类：按 dropped 集合；全机制（dropped 空）= overall（跳过，transport.json 已覆盖）。
_SUBSET_STYLE = {
    "intrinsic": ("phonon-limited", "-", "C0", "o"),
    "ADP_only": ("ADP only (acoustic limit)", ":", "C2", "^"),
    "ADP_POP": ("ADP+POP", "-", "C1", "v"),
    "ADP_POP_ODP": ("ADP+POP+ODP", "-", "C6", "P"),
}
_SUBSET_GENERIC = (("-", "C1", "v"), ("-", "C6", "P"), ("--", "C4", "D"), ("--", "C5", "s"))


def _subset_key(q):
    """subset_*.json -> 数据集键（按 dropped 集合分类）。"""
    kept = [str(x).upper() for x in (q.get("intrinsic_labels") or [])]
    dropped = [str(x).upper() for x in (q.get("intrinsic_dropped") or [])]
    dset = set(dropped)
    if not dset:
        return "overall"
    if dset == {"IMP"}:
        return "intrinsic"
    if kept == ["ADP"]:
        return "ADP_only"
    if set(kept) == {"ADP", "POP"}:
        return "ADP_POP"
    if set(kept) == {"ADP", "POP", "ODP"}:
        return "ADP_POP_ODP"
    return "+".join(kept)


def _subset_file_defs(run_dir):
    """[0038] subset_*.json -> [(key, fname, label, ls, colr, mk), ...]（全机制子集 = overall 跳过）。"""
    import glob as _glob
    out = []
    for i, fp in enumerate(sorted(_glob.glob(str(Path(run_dir) / "subset_*.json")))):
        q = _load_json(fp)
        if not q:
            continue
        key = _subset_key(q)
        if key == "overall":
            continue
        fname = Path(fp).name
        if key in _SUBSET_STYLE:
            label, ls, colr, mk = _SUBSET_STYLE[key]
        else:
            ls, colr, mk = _SUBSET_GENERIC[i % len(_SUBSET_GENERIC)]
            label = "+".join([str(x).upper() for x in (q.get("intrinsic_labels") or [])])
        out.append((key, fname, label, ls, colr, mk))
    return out


def _intrinsic_file(run_dir, key):
    """[0038] key（ADP_only / intrinsic / ADP_POP / ADP_POP_ODP）对应的结果文件名：
    subset_*.json（SUBSETS 生成）优先，退回旧的 intrinsic_*.json。找不到 -> None。"""
    import glob as _glob
    for fp in sorted(_glob.glob(str(Path(run_dir) / "subset_*.json"))):
        q = _load_json(fp)
        if q and _subset_key(q) == key:
            return Path(fp).name
    return dict((a, b) for a, b, *_ in INTR_SETS).get(key)
INTR_PANELS = (   # (行字典里的键, y 轴标签, 对数 y, 取绝对值)
    ("S_uV/K", "|S| (uV/K)", False, True),
    ("sigma_S/m", "sigma (S/m, per layer thickness t)", True, False),
    ("kappa_e_W/mK", "kappa_e (W/mK, per layer thickness t)", True, False),
    ("Lorenz_WOhm/K2", "Lorenz (1e-8 W Ohm/K^2)", False, False),   # 画图时 ×1e8（避免偏移文字压住上轴）
    ("mu_cm2/Vs", "mobility (cm^2/Vs)", True, False),
)
INTR_COLS = ("T_K", "type", "dataset", "n_cell_cm-3", "n_2D_cm-2", "n_lit_cm-3",
             "carriers_per_cell", "carriers_per_atom", "in_interpret_range", "S84_marks",
             "S_uV/K", "S_xx_uV/K", "S_yy_uV/K",
             "sigma_S/m", "sigma_xx_S/m", "sigma_yy_S/m", "sheet_sigma_S",
             "kappa_e_W/mK", "kappa_e_xx_W/mK", "kappa_e_yy_W/mK",
             "Lorenz_WOhm/K2", "PF_W/mK2", "mu_cm2/Vs", "mu_xx_cm2/Vs", "mu_yy_cm2/Vs", "flag")


def _geom_2d(cwd):
    """[0032] 2D 浓度/归一化换算的几何（按材料从 S8.4 的 2d_correction.json 读）：
    {"c_A", "t_A", "c_over_t", "t_method", "src"}；读不到 -> None。
    c = 含真空的胞高（AMSET/BT2 的 cm^-3 按它折算），t = 层厚（文献单层 cm^-3 与 σ、κe 归一按它）。"""
    for run in (AMSET2D_DIR, AMSET_DIR):
        rec = _load_json(Path(cwd) / run / "2d_correction.json")
        if not rec:
            continue
        c, f = _num(rec.get("cell_c_A")), _num(rec.get("elastic_rescale_factor_c_over_t"))
        if not c or not f or f <= 0:
            continue
        lt = rec.get("layer_thickness")
        # 可解读范围：与 S8.4（postprocess_intrinsic._doping_marks，2026-09-28 用户批准）同一判据 ——
        #   每【原子】载流子数 ≤ carrier_per_atom_threshold（默认 0.02）；缺键时退回每原胞 0.1。
        # [0034] 记录里缺 num_atoms / 阈值（老记录）时数结构文件、用默认 0.02，不再退回每原胞
        _na = _num(rec.get("num_atoms")) or _cell_atoms(cwd)
        return {"c_A": c, "t_A": c / f, "c_over_t": f, "src": run,
                "t_method": (lt.get("method") if isinstance(lt, dict) else None) or "?",
                "area_cm2": _num(rec.get("inplane_cell_area_cm2")),
                "num_atoms": _na,
                "thr_atom": (_num(rec.get("carrier_per_atom_threshold")) or CARRIER_PER_ATOM_MAX) if _na else None,
                "is_3d": False}
    return None


def _geom_3d(cwd):
    """[0046] 3D（体材料）浓度/归一化换算的几何：晶体胞没有真空层，
    c = 胞高（= 体积/底面 = 三重积/|a×b|，正交胞即 c 轴长，Å）；t = c（不做层厚归一，c_over_t = 1）。
    σ、κe、PF 用 AMSET 的体材料值；面密度列退化成「每胞截面」口径：
    n_2D = n_cell × c，n_lit = n_2D / t = n_cell（即体浓度本身）。读不到结构文件 -> None。"""
    A_cm2, V_cm3 = _cell_geometry(cwd)
    if not A_cm2 or not V_cm3:
        return None
    c_A = V_cm3 / A_cm2 * 1e8                      # cm -> Å
    if not c_A or c_A <= 0:
        return None
    _na = _cell_atoms(cwd)
    return {"c_A": c_A, "t_A": c_A, "c_over_t": 1.0, "src": AMSET_DIR, "is_3d": True,
            "t_method": "3D bulk cell (no layer thickness)",
            "area_cm2": A_cm2, "num_atoms": _na,
            "thr_atom": CARRIER_PER_ATOM_MAX if _na else None}


def _sha256_file(p):
    import hashlib
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _intr_same_run(run_dir, q):
    """[0032] 本征文件是否与当前 transport.json 同一次运行：postprocess 记下的本轮 transport_*.json 的 sha256
    必须等于当前 transport.json 的（作业末尾由它 cp 而来，内容相同）。不靠 mtime：drvfs/9p 回拉保不住时间戳。
    没有指纹的旧文件（[0032] 之前）只能靠网格判，返回 (True, None)。"""
    fp = (q or {}).get("transport_fingerprint") or {}
    if not fp.get("sha256"):
        return True, None
    tj = Path(run_dir) / "transport.json"
    try:
        return (_sha256_file(tj) == fp["sha256"]), fp.get("file")
    except OSError:
        return True, None


def _same_grid(q_dop, q_T, doping, temps):
    try:
        return (len(q_dop) == len(doping) and len(q_T) == len(temps)
                and all(abs(float(a) - b) <= 1e-6 * max(1.0, abs(b)) for a, b in zip(q_dop, doping))
                and all(abs(float(a) - b) < 1e-6 for a, b in zip(q_T, temps)))
    except (TypeError, ValueError):
        return False


def load_intrinsic_sets(cwd, run_name=AMSET2D_DIR):
    """[0032][0038][0046] 读可用结果集（2D 默认 S8.4 的 step8.4_amset2d；3D 传 run_name=step8_amset）。
    返回 {"doping", "temps", "sets": {键: {...}}, "notes": [...]}；该目录没有可用 transport.json -> None。
    每套：{"S","sigma","kappa","mu": arr[nd][nT][3][3] 或 None, "partial": 只有 μ, "labels", "reproduce"}。
    本征文件的掺杂/温度网格必须与 transport.json 相同（同一次运行），否则不用。
    [0038] 优先读 subset_*.json（SUBSETS 生成）；没 subset 文件的键退回旧的 intrinsic_*.json。"""
    j = _load_json(Path(cwd) / run_name / "transport.json")
    if not j:
        return None
    doping = [float(x) for x in _unwrap(j["doping"])]
    temps = [float(t) for t in _unwrap(j["temperatures"])]
    mob = j.get("mobility") or {}
    sets, notes = {}, []
    sets["overall"] = {"S": _unwrap(j.get("seebeck")), "sigma": _unwrap(j.get("conductivity")),
                       "kappa": _unwrap(j.get("electronic_thermal_conductivity")),
                       "mu": _unwrap(mob.get("overall") or mob.get("Overall")), "partial": False}
    run_dir = Path(cwd) / run_name
    # [0038] 新格式 subset_*.json（按 dropped 集合分类，全机制 = overall 跳过）；没覆盖到的键退回旧的 intrinsic_*.json。
    sub_defs = _subset_file_defs(run_dir)
    covered = {d[0] for d in sub_defs}
    file_defs = list(sub_defs)
    for key, fname, label, ls, colr, mk in INTR_SETS:
        if key != "overall" and key not in covered:
            file_defs.append((key, fname, label, ls, colr, mk))
    _file_of = {d[0]: d[1] for d in file_defs}
    extra = []
    for k_, fname in enumerate(_overall_without_files(run_dir)):   # [0034]
        ls, colr, mk = INTR_WITHOUT_STYLE[k_ % len(INTR_WITHOUT_STYLE)]
        extra.append((Path(fname).stem, fname, "overall w/o " + Path(fname).stem[len("overall_no"):], ls, colr, mk))
    for key, fname, *_ in file_defs + extra:
        q = _load_json(run_dir / fname)
        if not q:
            continue
        if not _same_grid(q.get("doping_cm3") or [], q.get("temperatures_K") or [], doping, temps):
            notes.append("%s：%s 的掺杂/温度网格与 transport.json 不同（不是同一次运行），不用" % (key, fname))
            continue
        same, src_name = _intr_same_run(run_dir, q)
        if not same:
            notes.append("%s：%s 是上一轮 S8.4 留下的（记录的 %s 与当前 transport.json 内容不同），不用"
                         % (key, fname, src_name or "transport"))
            continue
        it = q.get("intrinsic") or {}
        sets[key] = {"S": it.get("seebeck_uV_K"), "sigma": it.get("conductivity_S_m"),
                     "kappa": it.get("electronic_thermal_conductivity_W_mK"),
                     "mu": (it.get("mobility_cm2_Vs_s") or {}).get("overall"), "partial": False,
                     "labels": q.get("intrinsic_labels"), "reproduce": q.get("reproduction_check"),
                     "dropped": q.get("intrinsic_dropped")}
        if sets[key]["kappa"] is None:
            notes.append("%s：%s 是 [0018] 之前的旧文件，没有 κe —— κe/Lorenz 留空" % (key, fname))
    if "ADP_only" not in sets:
        a = _unwrap(mob.get("ADP") or mob.get("ACD"))
        if a is not None:          # 用户已有的"纯 ADP"：只有迁移率（AMSET 只按机制拆 μ，不拆 S/σ/κe）
            sets["ADP_only"] = {"S": None, "sigma": None, "kappa": None, "mu": a, "partial": True,
                                "labels": ["ADP"]}
            notes.append("ADP_only：没有 subset_ADP.json / intrinsic_ADP.json，只有 transport.json 的单机制 μ_ADP；"
                         "S/σ/κe 要 S8.4 开 WRITE_MESH 重跑后才有")
    if "intrinsic" not in sets:
        notes.append("intrinsic：没有可用的 subset（去 IMP）/ intrinsic_transport.json —— S8.4 没开 WRITE_MESH（[0032] 起 2D 默认开），"
                     "或是上一轮残留；重跑 S8.4 后才有")
    marks = {}
    for key in ("intrinsic", "ADP_only"):                   # S8.4 写的 doping_marks（每原子判据 + sign_anomaly）
        if key in sets and not sets[key].get("partial"):
            q = _load_json(run_dir / _file_of[key]) if key in _file_of else None
            marks = (q or {}).get("doping_marks") or {}
            if marks:
                break
    over_labels = [k for k in mob if str(k).lower() != "overall"]
    # 基础顺序按 INTR_SETS；被 subset 覆盖的键改用 subset 的元信息（真实文件名/图例/样式）。
    set_defs = []
    for key, fname, label, ls, colr, mk in INTR_SETS:
        _hit = next((d for d in sub_defs if d[0] == key), None)
        set_defs.append(_hit if _hit else (key, fname, label, ls, colr, mk))
    set_defs += [e for e in extra if e[0] in sets]                  # [0034] 只列真读到的
    for d in sub_defs:                          # [0038] subset 里出现、不在基础键里的组合追加
        if d[0] in sets and d[0] not in [x[0] for x in set_defs]:
            set_defs.append(d)
    # [0035] 动态加 ODP / IVS 单机制迁移率（transport.json 的 mobility 有这俩键才加；只有 μ，无 S/σ/κe）
    for _k, _ls, _colr, _mk in (("ODP", ":", "C6", "P"), ("IVS", "--", "C7", "X")):
        _m = _unwrap(mob.get(_k))
        if _m is not None:
            sets[_k] = {"S": None, "sigma": None, "kappa": None, "mu": _m, "partial": True, "labels": [_k]}
            set_defs.append((_k, None, _k + " (single-mech, mu only)", _ls, _colr, _mk))
    return {"doping": doping, "temps": temps, "sets": sets, "notes": notes, "marks": marks,
            "overall_labels": over_labels, "set_defs": set_defs}


def _tv(arr, i, k, ndim=2):
    """arr[nd][nT][3][3] -> (对角平均, xx, yy)；缺/NaN -> None。
    [0046] 2D 取面内 (xx+yy)/2；3D 取对角 (xx+yy+zz)/3（立方/各向同性时两者相同）。"""
    try:
        t = arr[i][k]
        xx, yy = float(t[0][0]), float(t[1][1])
        zz = float(t[2][2]) if ndim >= 3 else None
    except (TypeError, IndexError, ValueError):
        return None, None, None
    avg = (xx + yy) / 2.0 if zz is None else (xx + yy + zz) / 3.0
    v = [x if x == x else None for x in (avg, xx, yy)]
    return v[0], v[1], v[2]


def _interp_rule(geom):
    """可解读判据的文字（图注/表头用）。"""
    if geom.get("num_atoms") and geom.get("thr_atom"):
        return "> %.3g carriers per atom (S8.4 criterion)" % geom["thr_atom"], "每原子 ≤ %.3g（与 S8.4 同口径）" % geom["thr_atom"]
    return "> %.2g carriers per cell" % NCELL_INTERPRET_MAX, "每原胞 ≤ %.2g（数不出原子数，退回旧判据）" % NCELL_INTERPRET_MAX


def _dose_cols(n3_signed, geom, V_cm3, A_cm2=None):
    """浓度三口径 + 每原胞 / 每原子载流子数；可解读判据与 S8.4 相同（每原子，缺键退回每原胞 0.1）。"""
    n3 = abs(float(n3_signed))
    n2d = n3 * geom["c_A"] * 1e-8
    A = geom.get("area_cm2") or A_cm2
    ncell = (n2d * A) if A else (n3 * V_cm3 if V_cm3 else None)
    natom = (ncell / geom["num_atoms"]) if (ncell is not None and geom.get("num_atoms")) else None
    if natom is not None and geom.get("thr_atom"):
        ok = natom <= geom["thr_atom"]
    elif ncell is not None:
        ok = ncell <= NCELL_INTERPRET_MAX
    else:
        ok = None
    return {"n_cell_cm-3": float(n3_signed), "n_2D_cm-2": n2d, "n_lit_cm-3": n2d / (geom["t_A"] * 1e-8),
            "carriers_per_cell": ncell, "carriers_per_atom": natom,
            "in_interpret_range": (None if ok is None else ("yes" if ok else "NO"))}


def _row_from(vals, T, fac, geom):
    """vals = {"S":(avg,xx,yy), "sigma":(...), "kappa":(...), "mu":(...)}（元胞口径）-> 表行的输运列。"""
    r = {}
    S, sg, kp, mu = (vals.get(k) or (None, None, None) for k in ("S", "sigma", "kappa", "mu"))
    r["S_uV/K"], r["S_xx_uV/K"], r["S_yy_uV/K"] = S
    r["sigma_S/m"], r["sigma_xx_S/m"], r["sigma_yy_S/m"] = (None if v is None else v * fac for v in sg)
    r["kappa_e_W/mK"], r["kappa_e_xx_W/mK"], r["kappa_e_yy_W/mK"] = (None if v is None else v * fac for v in kp)
    r["mu_cm2/Vs"], r["mu_xx_cm2/Vs"], r["mu_yy_cm2/Vs"] = mu
    if sg[0] is not None:
        r["sheet_sigma_S"] = sg[0] * geom["c_A"] * 1e-10             # 片电导 = σ_元胞 × c，与口径无关
    if sg[0] and kp[0] is not None and sg[0] > 0:
        r["Lorenz_WOhm/K2"] = kp[0] / (sg[0] * T)                     # τ、口径都抵消
    if S[0] is not None and r.get("sigma_S/m") is not None:
        r["PF_W/mK2"] = (S[0] * 1e-6) ** 2 * r["sigma_S/m"]
    return r


def _bt2_tau(dpt_T, carrier):
    return _dpt_tau_s(dpt_T, carrier) if dpt_T else None


def build_intrinsic_rows(cwd, data, geom):
    """[0032] 长表：每个温度 × 原生掺杂档 × 数据集一行；BT2 行在同一掺杂档上插值（同号分支、不外推）；
    DPT 每个温度、载流子一行（与掺杂无关，只有 μ）。返回 (rows, {T: (bt, dpt_T)})。"""
    A_cm2, V_cm3 = _cell_geometry(cwd)
    fac = geom["c_over_t"]
    _nd = 3 if geom.get("is_3d") else 2            # [0046] 3D 的对角平均含 zz
    rows, extra, miss = [], {}, {"BT2": [], "DPT": []}
    marks = data.get("marks") or {}
    for k, T in enumerate(data["temps"]):
        dpt_T = load_dpt(cwd, T=T, quiet=True)
        bt = load_bt2(cwd, T=T)
        if bt and abs(float(bt["T_used"]) - T) > 1.0:
            bt = None                                    # BT2 没有这个温度：不跨温度凑
        if not bt:
            miss["BT2"].append(T)
        extra[T] = (bt, dpt_T)
        for i, n3 in enumerate(data["doping"]):
            if not n3 or n3 != n3:                       # 0 / NaN 掺杂：分不出 n/p，log 插值也没定义
                continue
            typ = "n" if n3 < 0 else "p"
            base = dict(_dose_cols(n3, geom, V_cm3, A_cm2), **{"T_K": T, "type": typ})
            if marks.get(str(i)):
                base["S84_marks"] = ";".join(marks[str(i)])
            for key, *_ in data.get("set_defs", INTR_SETS):
                st = data["sets"].get(key)
                if not st:
                    continue
                vals = {q: _tv(st.get(q), i, k, _nd) for q in ("S", "sigma", "kappa", "mu")}
                rows.append(dict(base, dataset=key, **_row_from(vals, T, fac, geom)))
            if bt:
                tau = _bt2_tau(dpt_T, "electron" if n3 < 0 else "hole")
                S = _interp_loglog(bt["n_signed"], bt["seebeck"], n3)
                L = _interp_loglog(bt["n_signed"], bt["lorenz"], n3)
                sot = _interp_loglog(bt["n_signed"], bt["sigma_over_tau"], n3)
                kot = _interp_loglog(bt["n_signed"], bt["kappa_e_over_tau"], n3)
                sg = (sot * tau) if (sot is not None and tau) else None
                kp = (kot * tau) if (kot is not None and tau) else None
                mu = (sg / (E_C * abs(n3) * 1e6) * 1e4) if sg else None
                # μ_BT2 = σ_总/(e|n_净|)：双极区（S 符号与载流子类型相反）分子含少子、分母是净浓度，
                #   与 AMSET 的多子迁移率不是同一个量 —— 不给 μ，打标记。
                bip = S is not None and ((S > 0) if typ == "n" else (S < 0))
                r = dict(base, dataset=BT2_KEY,
                         **_row_from({"S": (S, None, None), "sigma": (sg, None, None),
                                      "kappa": (kp, None, None), "mu": (None if bip else mu, None, None)},
                                     T, fac, geom))
                if bip:
                    r["flag"] = "bipolar(S sign): mu masked"
                if L is not None:
                    r["Lorenz_WOhm/K2"] = L                      # BT2 自己的 L（τ 无关），不依赖 DPT
                rows.append(r)
        for carrier, typ in (("electron", "n"), ("hole", "p")):
            mu = (dpt_T or {}).get(carrier)
            if isinstance(mu, (int, float)):
                mv = (dpt_T or {}).get("src_" + carrier) == "multi_valley"
                rows.append({"T_K": T, "type": typ, "dataset": "DPT_mv" if mv else "DPT", "mu_cm2/Vs": float(mu)})
            else:
                miss["DPT"].append("%.0f K %s" % (T, carrier))
    for key, lst in miss.items():
        if lst:
            data.setdefault("notes", []).append(
                "%s 缺这些温度（不外推）：%s" % (key, ", ".join("%.0f K" % x if isinstance(x, float) else x for x in lst)))
    return rows, extra


def write_intrinsic_table(out_dir, rows, geom, data):
    import csv
    with open(Path(out_dir) / INTR_CSV, "w", newline="", encoding="utf-8") as fh:
        fh.write("# [0032][0046] 本征输运参数表（%s 同一次运行；全部温度 × 原生掺杂档 × 数据集）\n"
                 % (AMSET_DIR if geom.get("is_3d") else AMSET2D_DIR))
        fh.write("# 浓度：n_cell_cm-3 = AMSET/AutoZT 口径（%s c=%.4f Å，负=n 型）；n_2D_cm-2 = n_cell×c；"
                 "n_lit_cm-3 = n_2D/t（%s，t=%.4f Å，%s）；n_lit = n_cell × c/t（c/t=%.4f）\n"
                 % ("晶体胞高" if geom.get("is_3d") else "含真空胞高", geom["c_A"],
                    "体材料浓度本身" if geom.get("is_3d") else "文献单层口径",
                    geom["t_A"], geom["t_method"], geom["c_over_t"]))
        fh.write("# %s；sheet_sigma_S = 片电导（与口径无关）；"
                 "μ、S、Lorenz 与口径无关；分方向列为 xx/yy\n"
                 % ("σ、κe、PF 用 AMSET 体材料值（3D，不做层厚归一，c/t = 1）" if geom.get("is_3d")
                    else "σ、κe、PF 按层厚 t 归一（= 元胞值 × c/t）"))
        fh.write("# in_interpret_range：%s；S84_marks = S8.4 本征文件的 doping_marks；"
                 "BT2 的 μ = σ_总/(e|n_净|)，双极区（S 符号与类型相反）不给（flag 列）\n" % _interp_rule(geom)[1])
        fh.write("# dataset（内部键 -> 图例）：ADP_only = ADP only（声学极限）；intrinsic = phonon-limited（去 IMP：%s）；"
                 "overall = phonon + impurity（全机制）；"
                 "%s = BoltzTraP2 CRTA × DPT τ（同号插值到本档，L 为 BT2 自身）；DPT / DPT_mv = 单谷 / 多谷 Bardeen–Shockley μ\n"
                 % ("+".join((data["sets"].get("intrinsic") or {}).get("labels") or ["ADP", "POP"]), BT2_KEY))
        w = csv.DictWriter(fh, fieldnames=list(INTR_COLS), extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else
                            ("%.6g" % r[c] if isinstance(r[c], float) else r[c])) for c in INTR_COLS})


def _intr_mechs(key, data):
    """数据集实际包含的散射机制（列表）；不知道 -> None。"""
    if key == "overall":
        return data.get("overall_labels") or None
    return ((data.get("sets") or {}).get(key) or {}).get("labels") or None


def _phys_label(mechs, all_mechs=None):
    """[0054] 统一命名（与电声文献一致）：物理含义 + 实际机制组合。
    含 IMP -> "phonon + impurity (ADP+IMP+POP)"；只 ADP -> "ADP only (acoustic limit)"；
    = 本轮全部声子机制（全机制去 IMP）-> "phonon-limited (ADP+POP+ODP)"（文献的 phonon-limited / intrinsic mobility）；
    其余声子子集 -> "partial phonon (ADP+POP)"；不知道本轮全机制 -> "phonon (…)"。mechs 为空 -> None。
    "intrinsic" 一词只留给本征载流子（n = p = n_i，intrinsic_carrier/）。"""
    ms = [str(m).upper() for m in (mechs or [])]
    if not ms:
        return None
    tag = "+".join(ms)
    if "IMP" in ms:
        return "phonon + impurity (%s)" % tag
    if ms == ["ADP"]:
        return "ADP only (acoustic limit)"
    full = {str(m).upper() for m in (all_mechs or [])} - {"IMP"}
    if not full:
        return "phonon (%s)" % tag
    return ("phonon-limited (%s)" if set(ms) == full else "partial phonon (%s)") % tag


def _intr_legend(key, data, default):
    """[0040][0054] 图例 = 物理含义 + 实际散射机制组合（SCATTERING 可被 step.conf 改；非极性材料没有 POP）。"""
    return _phys_label(_intr_mechs(key, data), (data or {}).get("overall_labels")) or default


def _intr_caption(geom, data):
    sets = data["sets"]
    lab = "+".join((sets.get("intrinsic") or {}).get("labels") or ["ADP", "POP"])
    part = " (ADP only: mobility only, no mesh)" if (sets.get("ADP_only") or {}).get("partial") else ""
    return ("Bottom axis: areal density n_2D (cm^-2).  Top axis: literature 3D density n_2D/t (cm^-3), "
            "t = %.3f A (%s).  AutoZT/AMSET files use n_cell = n_2D/c with c = %.3f A incl. vacuum, "
            "so n_lit = n_cell x c/t (c/t = %.3f).\n"
            "sigma, kappa_e normalized by layer thickness t (cell value x c/t); mobility, S, Lorenz do not depend on it.  "
            "Shaded: %s (beyond rigid-band range).\n"
            "ADP only: acoustic deformation potential only (upper limit)%s.  phonon-limited (%s): all phonon mechanisms, "
            "no ionized impurities (= 'intrinsic mobility' in e-ph literature).  "
            "phonon + impurity: all mechanisms incl. IMP%s.  BT2: BoltzTraP2 CRTA with DPT tau (S, Lorenz tau-free); "
            "BT2 mobility = sigma_total/(e|n_net|), masked where bipolar.  DPT: Bardeen-Shockley (single or multi-valley)."
            % (geom["t_A"], geom["t_method"], geom["c_A"], geom["c_over_t"], _interp_rule(geom)[0], part, lab,
               "".join("; w/o %s: same mesh with %s removed" % ("+".join(sets[k].get("dropped") or []),
                                                               "+".join(sets[k].get("dropped") or []))
                       for k in sets if k.startswith("overall_no"))))


def _bt2_curve(bt, dpt_T, typ, geom, xr):
    """BT2 密集曲线（只取 AMSET 掺杂范围附近）：返回 {列键: [(x=n_2D, y)]}。"""
    if not bt:
        return {}
    tau = _bt2_tau(dpt_T, "electron" if typ == "n" else "hole")
    fac, c_cm = geom["c_over_t"], geom["c_A"] * 1e-8
    out = {k: [] for k, *_ in INTR_PANELS}
    for i, n3 in enumerate(bt["n_signed"]):
        if n3 != n3 or n3 == 0 or (n3 < 0) != (typ == "n"):
            continue
        x = abs(n3) * c_cm
        if not (xr[0] / 3.0 <= x <= xr[1] * 3.0):
            continue
        S, L = bt["seebeck"][i], bt["lorenz"][i]
        if S == S:
            out["S_uV/K"].append((x, abs(S)))
        if L == L:
            out["Lorenz_WOhm/K2"].append((x, L))
        if tau:
            sg = bt["sigma_over_tau"][i] * tau
            kp = bt["kappa_e_over_tau"][i] * tau
            if sg == sg and sg > 0:
                out["sigma_S/m"].append((x, sg * fac))
                if not (S == S and ((S > 0) if typ == "n" else (S < 0))):     # 双极区不画 μ
                    out["mu_cm2/Vs"].append((x, sg / (E_C * abs(n3) * 1e6) * 1e4))
            if kp == kp and kp > 0:
                out["kappa_e_W/mK"].append((x, kp * fac))
    return {k: sorted(v) for k, v in out.items() if v}


def plot_intrinsic_T(path, T, rows, geom, data, bt, dpt_T, material=""):
    """[0032] 一个温度一张 2×5 图：行 n 型 / p 型，列 S、σ、κe、Lorenz、μ。2D 双横轴。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rT = [r for r in rows if abs(r["T_K"] - T) < 1e-6]
    c_t_cm = geom["t_A"] * 1e-8
    fig, ax = plt.subplots(2, len(INTR_PANELS), figsize=(23, 9.6))
    for ri, (typ, tname) in enumerate((("n", "n-type (electrons)"), ("p", "p-type (holes)"))):
        xs = [r["n_2D_cm-2"] for r in rT if r["type"] == typ and r.get("n_2D_cm-2")]
        if not xs:
            for a in ax[ri]:
                a.axis("off")
            continue
        xr = (min(xs), max(xs))
        lim = [r["n_2D_cm-2"] for r in rT if r["type"] == typ and r.get("in_interpret_range") == "NO"]
        btc = _bt2_curve(bt, dpt_T, typ, geom, xr)
        for ci, (col, ylab, logy, absval) in enumerate(INTR_PANELS):
            a = ax[ri][ci]
            ysc = 1e8 if col == "Lorenz_WOhm/K2" else 1.0
            for key, _f, leg, ls, colr, mk in data.get("set_defs", INTR_SETS):
                leg = _intr_legend(key, data, leg)
                pts = sorted((r["n_2D_cm-2"], (abs(r[col]) if absval else r[col]) * ysc) for r in rT
                             if r["dataset"] == key and r["type"] == typ and r.get(col) is not None
                             and r[col] == r[col] and (not logy or (abs(r[col]) if absval else r[col]) > 0))
                if pts:
                    a.plot(*zip(*pts), ls, color=colr, marker=mk, ms=5, lw=1.6, label=leg)
            if col in btc:
                a.plot(*zip(*[(x, y * ysc) for x, y in btc[col]]), "--", color="C1", lw=1.3, alpha=.85,
                       label=BT2_LABEL)
            if col == "mu_cm2/Vs":
                _car = "electron" if typ == "n" else "hole"
                mu_d = (dpt_T or {}).get(_car)
                if isinstance(mu_d, (int, float)) and mu_d > 0:
                    a.axhline(mu_d, color="k", lw=1.0, ls=(0, (1, 1.5)),
                              label="DPT (%s)" % ("multi-valley" if (dpt_T or {}).get("src_" + _car) == "multi_valley"
                                                   else "single valley"))
            if col == "Lorenz_WOhm/K2":
                a.axhline(SOMMERFELD_L * 1e8, color="grey", lw=.8, ls=":", label="Sommerfeld 2.44e-8")
            if lim:
                a.axvspan(min(lim), xr[1] * 3.0, color="grey", alpha=.12, lw=0)
            a.set_xscale("log")
            a.set_xlim(xr[0] / 2.0, xr[1] * 2.0)
            if logy:
                a.set_yscale("log")
            a.set_ylabel(("%s:  " % tname.split(" ")[0] if ci == 0 else "") + ylab, fontsize=9)
            a.grid(alpha=.3)
            try:      # matplotlib < 3.1 没有 secondary_xaxis
                top = a.secondary_xaxis("top", functions=(lambda x, _t=c_t_cm: x / _t,
                                                          lambda x, _t=c_t_cm: x * _t))
            except Exception:                          # noqa: BLE001
                top = None
            if ri == 0:
                if top is not None:
                    top.set_xlabel("n_2D / t (cm$^{-3}$, literature)", fontsize=8.5)
                a.set_title(ylab.split(" (")[0], fontsize=11, pad=24)
            a.set_xlabel("n_2D (cm$^{-2}$)", fontsize=9)
            if ci == 0 or col in ("mu_cm2/Vs", "Lorenz_WOhm/K2"):
                a.legend(fontsize=7, loc="best")
    fig.suptitle("Phonon-limited transport @ %.0f K%s  (S8.4 ADP only / phonon-limited / phonon + impurity; BT2 x DPT tau)"
                 % (T, (" - " + material) if material else ""), fontsize=13)
    fig.text(0.006, 0.006, _intr_caption(geom, data), fontsize=8.2, va="bottom", ha="left", wrap=True)
    fig.tight_layout(rect=[0, 0.075, 1, 0.96])
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _intr_summary(out_dir, rows, geom, data, files):
    sets = data["sets"]
    L = ["# [0032] 本征输运汇总（详见 %s；图 %s）" % (INTR_CSV, ", ".join(files[:2]) + (" …" if len(files) > 2 else "")),
         "# 几何：c = %.4f Å（含真空，AMSET 的 cm^-3 用它）；t = %.4f Å（%s，文献 cm^-3 与 σ/κe 归一用它）；c/t = %.4f"
         % (geom["c_A"], geom["t_A"], geom["t_method"], geom["c_over_t"]),
         "# 换算：n_2D = n_cell × c；n_lit = n_2D / t = n_cell × c/t",
         "# 可解读范围：%s" % _interp_rule(geom)[1]]
    for key, fname, leg, *_ in data.get("set_defs", INTR_SETS):
        st = sets.get(key)
        if not st:
            L.append("# %-9s 缺：%s" % (key, fname))
            continue
        rep = st.get("reproduce") or {}
        md = rep.get("max_rel_diff") or rep.get("max_abs_diff") or {}     # [0044] 判据是相对差（旧文件只有绝对差）
        L.append("# %-9s %s%s%s" % (key, fname if not st.get("partial") else "transport.json mobility.ADP（只有 μ）",
                                    (" 机制=" + "+".join(st["labels"])) if st.get("labels") else "",
                                    ("；复现检查 %s（max|diff| %.1e）" % ("PASS" if rep.get("passed") else "FAIL",
                                                                     max(md.values()) if md else float("nan")))
                                    if rep else ""))
    for n_ in data.get("notes") or []:
        L.append("# [注] " + n_)
    wo = [d[0] for d in data.get("set_defs", INTR_SETS) if d[0].startswith("overall_no")]   # [0034]
    hdr = "%5s %4s %10s %10s | %10s %10s %10s%s %10s %8s | %9s %10s %10s %10s" % (
        "T(K)", "type", "n2D(cm-2)", "nlit(cm-3)", "mu_ADP", "mu_intr", "mu_all",
        "".join(" %10s" % ("mu_" + k[len("overall_"):])[:10] for k in wo), "mu_BT2", "mu_DPT",
        "S_intr", "sig_intr", "ke_intr", "PF_intr")
    L.append(hdr)

    def g(T, typ, n2d, ds, col):
        for r in rows:
            if (r["dataset"] == ds and r["type"] == typ and abs(r["T_K"] - T) < 1e-6
                    and (n2d is None or abs(r.get("n_2D_cm-2", 0) - n2d) <= 1e-6 * n2d)):
                return r.get(col)
        return None

    def f(v, fmt="%10.4g"):
        return (fmt % v) if isinstance(v, (int, float)) and v == v else (" " * (int(fmt[1:].split(".")[0]) - 1) + "-")

    for T in data["temps"]:
        for typ in ("n", "p"):
            dpt_mu = g(T, typ, None, "DPT", "mu_cm2/Vs")
            if dpt_mu is None:
                dpt_mu = g(T, typ, None, "DPT_mv", "mu_cm2/Vs")
            for n2d in sorted({r["n_2D_cm-2"] for r in rows if r.get("n_2D_cm-2") and r["type"] == typ
                               and abs(r["T_K"] - T) < 1e-6}):
                nlit = n2d / (geom["t_A"] * 1e-8)
                L.append("%5.0f %4s %10.3g %10.3g | %s %s %s%s %s %s | %s %s %s %s" % (
                    T, typ, n2d, nlit, f(g(T, typ, n2d, "ADP_only", "mu_cm2/Vs")),
                    f(g(T, typ, n2d, "intrinsic", "mu_cm2/Vs")), f(g(T, typ, n2d, "overall", "mu_cm2/Vs")),
                    "".join(" " + f(g(T, typ, n2d, k, "mu_cm2/Vs")) for k in wo),
                    f(g(T, typ, n2d, BT2_KEY, "mu_cm2/Vs")), f(dpt_mu, "%8.4g"),
                    f(g(T, typ, n2d, "intrinsic", "S_uV/K"), "%9.4g"), f(g(T, typ, n2d, "intrinsic", "sigma_S/m")),
                    f(g(T, typ, n2d, "intrinsic", "kappa_e_W/mK")), f(g(T, typ, n2d, "intrinsic", "PF_W/mK2"))))
    L.append("# 单位：μ cm²/Vs；S µV/K；σ S/m、κe W/mK、PF W/mK²（按层厚 t 归一）")
    (Path(out_dir) / INTR_TXT).write_text("\n".join(L) + "\n", encoding="utf-8")
    return L


def write_intrinsic_suite(out, cwd, is_2d):
    """[0032][0046] 入口：2D 读 S8.4、3D 读 S8；缺结果 / 几何时打印原因后跳过。返回生成的文件名列表。
    先清掉 intrinsic/ 里上一轮的图表（温度网格变了、本轮跳过或失败时，旧文件不能冒充本轮结果）。"""
    odir = Path(out) / INTR_SUBDIR
    if odir.is_dir():
        _old = list(odir.glob("intrinsic_transport_*K.png")) + [odir / INTR_CSV, odir / INTR_TXT]
        _n = 0
        for f in _old:
            if f.is_file():
                f.unlink()
                _n += 1
        if _n:
            print("[..] 本征五参数：清掉上一轮的 %d 个文件（%s/）" % (_n, INTR_SUBDIR))
    run_name = AMSET2D_DIR if is_2d else AMSET_DIR      # [0046] 3D 读 S8（step8_amset）
    data = load_intrinsic_sets(cwd, run_name)
    if not data:
        print("[..] 本征五参数图：没有 %s/transport.json，跳过" % run_name)
        return []
    geom = _geom_2d(cwd) if is_2d else _geom_3d(cwd)
    if not geom:
        print("[WARN] 本征五参数图：%s，无法换算浓度口径，跳过"
              % ("读不到 2d_correction.json 的 cell_c_A / c/t" if is_2d else "读不到结构文件（_cell_geometry）"))
        return []
    odir.mkdir(exist_ok=True)
    rows, extra = build_intrinsic_rows(cwd, data, geom)
    write_intrinsic_table(odir, rows, geom, data)
    material = Path(cwd).resolve().parent.name
    files = []
    for T in data["temps"]:
        name = "intrinsic_transport_%dK.png" % round(T)
        try:
            bt, dpt_T = extra.get(T, (None, None))
            plot_intrinsic_T(odir / name, T, rows, geom, data, bt, dpt_T, material=material)
            files.append(name)
        except Exception as e:                        # noqa: BLE001
            print("[WARN] 本征五参数图 %.0f K 失败（%s: %s）" % (T, type(e).__name__, e))
    for l_ in _intr_summary(odir, rows, geom, data, files)[:6]:
        print("[..] " + l_)
    print("[OK] 本征五参数：%s/%s + %s + %d 张图（%s）"
          % (INTR_SUBDIR, INTR_CSV, INTR_TXT, len(files), ", ".join(sorted(data["sets"]))))
    return [INTR_CSV, INTR_TXT] + files


# ---------- [0035] 电子弛豫时间 τ：S8.4 的 tau_energy.json vs S8.2 DPT 的 τ ----------
TAU_SUBDIR = "tau"
TAU_JSON = "tau_energy.json"
TAU_CSV = "tau_comparison.csv"
TAU_TXT = "tau_summary.txt"
_TAU_COLORS = {"ADP": "#1f77b4", "POP": "#d62728", "IMP": "#2ca02c", "PIE": "#9467bd", "ODP": "#ff7f0e", "IVS": "#8c564b"}


def load_tau(cwd):
    """tau_energy.json（status ok、且与当前 transport.json 同一次运行）-> dict；否则打印原因返回 None。"""
    run_dir = Path(cwd) / AMSET2D_DIR
    j = _load_json(run_dir / TAU_JSON)
    if not j:
        return None
    if j.get("status") != "ok":
        print("[WARN] %s/%s：τ 后处理失败（%s）—— 不出 τ 对比" % (AMSET2D_DIR, TAU_JSON, j.get("error", j.get("status"))))
        return None
    same, src = _intr_same_run(run_dir, j)
    if not same:
        print("[WARN] %s/%s 不是当前 transport.json 那次运行的（来自 %s）—— 不出 τ 对比" % (AMSET2D_DIR, TAU_JSON, src))
        return None
    return j


def _dpt_tau_fs(cwd, T):
    """S8.2 的分方向 τ（fs）在温度 T：{载流子: (τx, τy)}；读不到的方向为 None。"""
    d = load_dpt(cwd, T=T, quiet=True)
    out = {}
    for car in ("electron", "hole"):
        tx, ty = _dpt_tau_dir(d, car, "x"), _dpt_tau_dir(d, car, "y")
        out[car] = tuple(None if v is None else v * 1e15 for v in (tx, ty))
    return out


def _mean_xy(t):
    vals = [v for v in (t or ()) if isinstance(v, (int, float))]
    return sum(vals) / len(vals) if vals else None


def _lowest_doping(tj, car):
    idx = [n for n, c in enumerate(tj["carrier"]) if c == car]
    return min(idx, key=lambda n: abs(tj["doping_cm3"][n])) if idx else None


def build_tau_rows(cwd, tj, geom):
    rows = []
    mechs = [m for m in tj["mechanisms"] if m in tj["tau_eff_fs"]]
    for t, T in enumerate(tj["temperatures_K"]):
        dpt = _dpt_tau_fs(cwd, T)
        for n, car in enumerate(tj["carrier"]):
            ov = tj["tau_eff_fs"]["overall"]
            r = {"T_K": T, "carrier": car, "doping_cm-3": tj["doping_cm3"][n]}
            if geom:
                r["n_2D_cm-2"] = abs(tj["doping_cm3"][n]) * geom["c_A"] * 1e-8
            for d in ("xx", "yy", "inplane"):
                r["tau_eff_overall_%s_fs" % d] = ov[d][n][t]
            for m in mechs:
                r["tau_eff_%s_inplane_fs" % m] = tj["tau_eff_fs"][m]["inplane"][n][t]
            r["tau_sigma_fs"] = tj["tau_sigma_fs"][n][t]
            tx, ty = dpt.get(car, (None, None))
            r["tau_DPT_x_fs"], r["tau_DPT_y_fs"] = tx, ty
            r["tau_DPT_inplane_fs"] = _mean_xy((tx, ty))
            a = r.get("tau_eff_ADP_inplane_fs")
            r["ADP_over_DPT"] = (a / r["tau_DPT_inplane_fs"]) if (a and r["tau_DPT_inplane_fs"]) else None
            rows.append(r)
    return rows, mechs


def write_tau_table(odir, rows):
    import csv
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(odir / TAU_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else (("%.6g" % r[k]) if isinstance(r[k], float) else r[k]))
                        for k in keys})


def plot_tau_E(path, tj, T_idx, cwd):
    """T 下每种载流子（最低掺杂档）的逐机制 τ(E)；横线：AMSET τ_eff（overall、ADP）与 S8.2 DPT 的 τ。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    T = tj["temperatures_K"][T_idx]
    x = tj["energy_bins_eV"]["centers"]
    dpt = _dpt_tau_fs(cwd, T)
    cars = [c for c in ("electron", "hole") if _lowest_doping(tj, c) is not None]
    fig, axes = plt.subplots(1, len(cars), figsize=(6.2 * len(cars), 4.6), squeeze=False)
    for ax, car in zip(axes[0], cars):
        n = _lowest_doping(tj, car)
        for m in tj["mechanisms"]:
            y = tj["tau_E_fs"][m][n][T_idx]
            pts = [(a, b) for a, b in zip(x, y) if b is not None and b > 0]
            if pts:
                ax.plot(*zip(*pts), color=_TAU_COLORS.get(m), label="τ_%s(E)" % m)
        pts = [(a, b) for a, b in zip(x, tj["tau_E_total_fs"][n][T_idx]) if b is not None and b > 0]
        if pts:
            ax.plot(*zip(*pts), color="k", lw=2, label="τ_total(E)")
        te = tj["tau_eff_fs"]["overall"]["inplane"][n][T_idx]
        if te:
            ax.axhline(te, color="k", ls="--", lw=1, label="AMSET τ_eff overall = %.3g fs" % te)
        ta = (tj["tau_eff_fs"].get("ADP") or {}).get("inplane")
        if ta and ta[n][T_idx]:
            ax.axhline(ta[n][T_idx], color=_TAU_COLORS["ADP"], ls="--", lw=1,
                       label="AMSET τ_eff ADP = %.3g fs" % ta[n][T_idx])
        td = _mean_xy(dpt.get(car))
        if td:
            ax.axhline(td, color="0.4", ls=":", lw=1.5, label="S8.2 DPT τ (used by BT2) = %.3g fs" % td)
        ax.set_yscale("log")
        ax.set_xlabel("E - CBM (eV)" if car == "electron" else "VBM - E (eV)")
        ax.set_ylabel("τ (fs)")
        ax.set_title("%s  %.0f K  doping %.2e cm$^{-3}$" % (car, T, tj["doping_cm3"][n]))
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7)
    fig.suptitle("τ(E) = 1/<Γ>_E (S8.4 mesh);  τ_eff = constant τ that reproduces the AMSET mobility", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_tau_T(path, rows, tj):
    """最低掺杂档：τ_eff（overall、ADP）与 DPT τ 随温度。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cars = [c for c in ("electron", "hole") if _lowest_doping(tj, c) is not None]
    fig, axes = plt.subplots(1, len(cars), figsize=(5.4 * len(cars), 4.2), squeeze=False)
    for ax, car in zip(axes[0], cars):
        dop = tj["doping_cm3"][_lowest_doping(tj, car)]
        rr = [r for r in rows if r["carrier"] == car and r["doping_cm-3"] == dop]
        for key, lab, sty in (("tau_eff_overall_inplane_fs", "AMSET τ_eff overall", dict(color="k")),
                              ("tau_eff_ADP_inplane_fs", "AMSET τ_eff ADP", dict(color=_TAU_COLORS["ADP"])),
                              ("tau_DPT_inplane_fs", "S8.2 DPT τ", dict(color="0.4", ls=":"))):
            pts = [(r["T_K"], r[key]) for r in rr if r.get(key)]
            if pts:
                ax.plot(*zip(*pts), marker="o", ms=3, label=lab, **sty)
        ax.set_yscale("log")
        ax.set_xlabel("T (K)")
        ax.set_ylabel("τ (fs)")
        ax.set_title("%s  doping %.2e cm$^{-3}$" % (car, dop))
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_tau_suite(out, cwd, is_2d):
    """[0035] 入口：S8.4 有 tau_energy.json 才做。先清 tau/ 里上一轮的产物。返回生成的文件名列表。"""
    odir = Path(out) / TAU_SUBDIR
    if odir.is_dir():
        for f in list(odir.glob("tau_*.png")) + [odir / TAU_CSV, odir / TAU_TXT]:
            if f.is_file():
                f.unlink()
    if not is_2d:
        return []
    tj = load_tau(cwd)
    if not tj:
        print("[..] τ 对比：没有可用的 %s/%s（S8.4 没开 WRITE_MESH/TAU_ENERGY，或还是旧版作业），跳过" % (AMSET2D_DIR, TAU_JSON))
        return []
    odir.mkdir(exist_ok=True)
    rows, mechs = build_tau_rows(cwd, tj, _geom_2d(cwd))
    write_tau_table(odir, rows)
    files = [TAU_CSV]
    t_idx = min(range(len(tj["temperatures_K"])), key=lambda i: abs(tj["temperatures_K"][i] - TARGET_T))
    for name, fn in (("tau_vs_E_%dK.png" % round(tj["temperatures_K"][t_idx]), lambda p: plot_tau_E(p, tj, t_idx, cwd)),
                     ("tau_eff_vs_T.png", lambda p: plot_tau_T(p, rows, tj))):
        try:
            fn(odir / name)
            files.append(name)
        except Exception as e:                        # noqa: BLE001
            print("[WARN] τ 图 %s 失败（%s: %s）" % (name, type(e).__name__, e))
    T = tj["temperatures_K"][t_idx]
    lines = ["# 电子弛豫时间（S8.4 tau_energy.json vs S8.2 DPT），%.0f K、各载流子最低掺杂档" % T,
             "# τ_eff = 常数 τ 的玻尔兹曼输运复现 AMSET 迁移率所需的 τ（面内）；τ_eff,ADP 与 DPT 的 τ 同为声学形变势口径",
             "# τ_σ/τ_eff 最大相对差 %.1f%%（核对项，见 tau_energy.json 的 notes）"
             % (100 * tj["check"].get("tau_sigma_vs_eff_max_rel", float("nan")))]
    for car in ("electron", "hole"):
        n = _lowest_doping(tj, car)
        if n is None:
            continue
        r = next(x for x in rows if x["carrier"] == car and x["T_K"] == T and x["doping_cm-3"] == tj["doping_cm3"][n])
        f = lambda v: "—" if v is None else "%.3g" % v                                   # noqa: E731
        lines.append("%-8s doping %.2e cm^-3：τ_eff overall %s fs；%s；DPT τ %s fs（x %s / y %s）；ADP/DPT = %s"
                     % (car, r["doping_cm-3"], f(r["tau_eff_overall_inplane_fs"]),
                        "；".join("%s %s" % (m, f(r.get("tau_eff_%s_inplane_fs" % m))) for m in mechs),
                        f(r["tau_DPT_inplane_fs"]), f(r["tau_DPT_x_fs"]), f(r["tau_DPT_y_fs"]), f(r["ADP_over_DPT"])))
    lines.append("# ADP/DPT 偏离 1：AMSET 的 ADP 逐 k 点算（含能量依赖、真实能带、波函数重叠），DPT 是抛物带解析式；"
                 "相差 2 倍以内属正常口径差，差太多先查 S8.2 的 m*/E1")
    (odir / TAU_TXT).write_text("\n".join(lines) + "\n", encoding="utf-8")
    files.append(TAU_TXT)
    for l_ in lines[3:]:
        print("[..] " + l_)
    print("[OK] τ 对比：%s/（%s）" % (TAU_SUBDIR, ", ".join(files)))
    return files


def make_figure(out, am, bt, dpt, am2=None, geom=None):
    """geom（[0032]，2D 才给）：横轴改成面密度 n_2D（cm^-2），上轴加文献口径 n_2D/t（cm^-3）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xf = geom["c_A"] * 1e-8 if geom else 1.0          # n_cell(cm^-3) -> n_2D(cm^-2)
    XL = "n_2D (cm$^{-2}$)" if geom else "|carrier conc| (cm$^{-3}$)"

    def split(sig, y):
        n = [(abs(sig[i]) * xf, y[i]) for i in range(len(sig))
             if sig[i] < 0 and y[i] == y[i]]
        p = [(abs(sig[i]) * xf, y[i]) for i in range(len(sig))
             if sig[i] > 0 and y[i] == y[i]]
        return sorted(n), sorted(p)

    def top_axis(a):                                   # [0032] 2D：上轴 = 文献口径 n_2D/t
        if geom:
            _t = geom["t_A"] * 1e-8
            try:      # matplotlib < 3.1 没有 secondary_xaxis：只丢上轴，不能让完成标记这张图整张失败
                a.secondary_xaxis("top", functions=(lambda x: x / _t, lambda x: x * _t)).set_xlabel(
                    "n_2D / t (cm$^{-3}$, literature)", fontsize=8)
            except Exception as _e:                    # noqa: BLE001
                print("[WARN] 上横轴（文献口径 cm^-3）画不出（%s: %s），只留 n_2D 横轴" % (type(_e).__name__, _e))

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
    a.set_xscale("log"); a.set_xlabel(XL); top_axis(a)
    a.set_ylabel("|Seebeck| (uV/K)"); a.set_title("Seebeck @300K (directly comparable)")
    a.legend(fontsize=8); a.grid(alpha=.3)
    # (0,1) Lorenz（[0032] 以 1e-8 为单位画：科学计数偏移文字会压住 2D 的上横轴）
    a = ax[0, 1]
    _L8 = (lambda v: [x * 1e8 for x in v]) if geom else (lambda v: v)      # 3D 不变
    if am:
        n, p = split(am["doping_signed"], _L8(am["lorenz"]))
        if n: a.plot(*zip(*n), "o-", c="C0", label="amset n")
        if p: a.plot(*zip(*p), "s-", c="C1", label="amset p")
    if bt:
        n, p = split(bt["n_signed"], _L8(bt["lorenz"]))
        if n: a.plot(*zip(*n), "--", c="C0", alpha=.7, label="BT2 n")
        if p: a.plot(*zip(*p), "--", c="C1", alpha=.7, label="BT2 p")
    a.axhline(SOMMERFELD_L * (1e8 if geom else 1.0), ls=":", c="k", label="Sommerfeld 2.44e-8")
    a.set_xscale("log"); a.set_xlabel(XL); top_axis(a)
    a.set_ylabel("Lorenz (1e-8 W.Ohm/K^2)" if geom else "Lorenz (W.Ohm/K^2)")
    a.set_title("Lorenz @300K (directly comparable)")
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
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(XL); top_axis(a)
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
            pts = sorted((abs(bt["n_signed"][i]) * xf, bt["kappa_e_over_tau"][i] * tau)
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
                a.plot(ntar * xf, SOMMERFELD_L * sig * 300.0, ":", c=col, alpha=.6,
                       label="DPT %s (WF)" % nm)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(XL); top_axis(a)
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
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(XL); top_axis(a)
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
    if geom:                                           # [0032] 2D 浓度口径图注
        notes += ("\n\n2D: bottom n_2D (cm^-2); top n_2D/t\n"
                  " (lit. cm^-3, t=%.2f A %s). AMSET\n"
                  " cm^-3 = n_2D/c (c=%.2f A, vacuum):\n"
                  " n_lit = n_cell x c/t (x%.3f).\n"
                  " sigma, kappa_e: per cell (x c/t)."
                  % (geom["t_A"], geom["t_method"], geom["c_A"], geom["c_over_t"]))
    a.text(0.0, (1.0 if geom else 0.98), notes, va="top", ha="left", fontsize=(7.2 if geom else 8.5),
           family="monospace", transform=a.transAxes)

    fig.suptitle("Electronic transport comparison @ ~300 K "
                 "(amset / BoltzTraP2-CRTA / DPT)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out / "comparison_300K.png", dpi=130)
    plt.close(fig)




def _ds_file(ds, data):
    """[0040] datasets/ 的文件名后缀 = 实际机制组合（如 ADP+POP+ODP）；不知道机制时用数据集键。"""
    lab = _intr_mechs(ds, data or {})
    return "+".join(lab) if lab else ds


UNDOPED_DIR = "intrinsic_carrier"     # [0040] transport_plots/ 下本征载流子（doping = 0）图的子目录


def _tavg(t, is_2d):
    """3x3 -> 面内 (xx+yy)/2（2D）或 (xx+yy+zz)/3（3D）；坏值 -> None。"""
    try:
        d = [float(t[i][i]) for i in ((0, 1) if is_2d else (0, 1, 2))]
    except (TypeError, IndexError, ValueError):
        return None
    v = sum(d) / len(d)
    return v if v == v else None


UNDOPED_NI_JSON = "undoped_ni.json"   # [0052] amset_fermi_fix 在 AMSET 运行目录写的本征点 n_i(T) 记录
_NI_LOG_RE = None


def undoped_ni(cwd, is_2d):
    """[0052] 本征点的本征载流子浓度 n_i(T) 与 E_F(T)：{T: {"n_i_cm-3", "E_F_eV"}}, 来源说明。
    transport.json 没有载流子浓度字段：先读 amset_fermi_fix 写的 undoped_ni.json；旧结果（[0052] 之前）
    退回解析运行目录日志里的「[fermi_fix] 本征点 <T> K：E_F = … eV，n = p = … cm^-3」行（取最新的日志）。
    n_i 是每胞体积口径（2D = 含真空的元胞）；都读不到 -> ({}, 原因)。"""
    import re
    global _NI_LOG_RE
    run = Path(cwd) / (AMSET2D_DIR if is_2d else AMSET_DIR)
    j = _load_json(run / UNDOPED_NI_JSON)
    if j and j.get("temperatures_K"):
        out = {}
        _ev, _ec = j.get("E_v_eV") or [], j.get("E_c_eV") or []       # [0054] 带边（新记录才有）
        for q, (T, ef, ni) in enumerate(zip(j["temperatures_K"], j.get("E_F_eV") or [], j.get("n_i_cm-3") or [])):
            out[float(T)] = {"n_i_cm-3": _num(ni), "E_F_eV": _num(ef),
                             "E_v_eV": _num(_ev[q]) if q < len(_ev) else None,
                             "E_c_eV": _num(_ec[q]) if q < len(_ec) else None}
        if out:
            return out, "%s/%s" % (run.name, UNDOPED_NI_JSON)
    if _NI_LOG_RE is None:
        # [0052-r1] 放宽：允许各种空白与全/半角冒号逗号（以前只认打印文案里那一种写法，
        #   改一个空格就整条回退失效、n_i 全标不出来）
        _NI_LOG_RE = re.compile(r"\[fermi_fix\]\s*本征点\s*([0-9.eE+-]+)\s*K\s*[：:]"
                                r"\s*E_F\s*=\s*([0-9.eE+-]+)\s*eV\s*[，,;；]?\s*"
                                r"n\s*=\s*p\s*=\s*([0-9.eE+-]+)\s*cm\^-3")
    logs = []
    if run.is_dir():
        for f in run.iterdir():
            if f.is_file() and (f.suffix in (".out", ".log") or f.name.startswith("queue")):
                try:
                    logs.append((f.stat().st_mtime, f))
                except OSError:
                    pass
    for _mt, f in sorted(logs, reverse=True):
        try:
            txt = f.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        out = {}
        for m in _NI_LOG_RE.finditer(txt):
            out[float(m.group(1))] = {"n_i_cm-3": float(m.group(3)), "E_F_eV": float(m.group(2))}
        if out:
            return out, "%s/%s（日志行）" % (run.name, f.name)
    return {}, ("运行目录里没有 %s，日志里也没有 [fermi_fix] 本征点 行（本征点关了，或不是 amset_fermi_fix 跑的）"
                % UNDOPED_NI_JSON)


KB_EV = 8.617333262e-5      # eV/K
DILUTE_ETA_MIN = 5.0        # [0054] η = (E_c − E_F)/k_BT 与 (E_F − E_v)/k_BT 的下限：低于它 f(带边) > ~0.7%，不再是 f ≈ 0


def _eta_of(v, T):
    """[0054] 本征点记录 v（_ni_lookup 的结果）-> (η_c, η_v, f_CBM)；缺带边（旧记录）-> (None, None, None)。"""
    ec, ev, ef = v.get("E_c_eV"), v.get("E_v_eV"), v.get("E_F_eV")
    if ec is None or ev is None or ef is None or not T or T <= 0:
        return None, None, None
    kt = KB_EV * float(T)
    eta_c, eta_v = (ec - ef) / kt, (ef - ev) / kt
    return eta_c, eta_v, 1.0 / (1.0 + math.exp(min(eta_c, 700.0)))


def _shade_not_dilute(ax, bad_T, all_T):
    """[0054] 在 η < DILUTE_ETA_MIN 的温度段打灰底；横轴固定为全部温度（空面板也不被灰底撑歪）。"""
    if not bad_T or not all_T:
        return
    ts = sorted(set(float(t) for t in all_T))         # 行按数据集重复，先去重再算半格宽
    step = (ts[-1] - ts[0]) / max(len(ts) - 1, 1) / 2.0 if len(ts) > 1 else 10.0
    ax.axvspan(min(bad_T) - step, max(bad_T) + step, color="0.88", zorder=0, lw=0)
    ax.set_xlim(ts[0] - step, ts[-1] + step)


def _ni_lookup(ni, T):
    """n_i 记录里与 T 最接近（差 < 0.5 K）的那条；没有 -> {}。"""
    for k, v in ni.items():
        if abs(k - float(T)) < 0.5:
            return v
    return {}


def undoped_rows(cwd, is_2d):
    """[0044] 本征点（AMSET doping = 0，E_F 由 n = p 定）：从主 transport.json 取 0 掺杂行
    （[0040] 起掺杂网格含 0；[0044] 取消独立拆分脚本）。2D 会带全部可用数据集（overall / 去 IMP / 子集）。
    兼容旧结果：主网格没有 0 点时回退读 transport_undoped.json（[0040]~[0043] 的产物）。"""
    try:
        rows, _why = _rows_at_doping(cwd, is_2d, lambda d: abs(d) <= 1e-12,
                                     "transport.json 的掺杂网格里没有 0（UNDOPED_POINT 关了）")
    except Exception:
        rows = []
    if rows:
        return rows, None
    run = Path(cwd) / (AMSET2D_DIR if is_2d else AMSET_DIR)
    j = _load_json(run / "transport_undoped.json")
    if not j:
        return [], "transport.json 的掺杂网格里没有 0（UNDOPED_POINT 关了），也没有旧的 transport_undoped.json"
    temps = [float(t) for t in _unwrap(j["temperatures"])]
    fac = 1.0
    if is_2d:
        geom = _geom_2d(cwd)
        fac = geom["c_over_t"] if geom else 1.0
    S_arr = _unwrap(j.get("seebeck"))
    sg_arr = _unwrap(j.get("conductivity"))
    kp_arr = _unwrap(j.get("electronic_thermal_conductivity"))
    rows = []
    for k, T in enumerate(temps):
        try:
            S = _tavg(S_arr[0][k], is_2d)
            sg = _tavg(sg_arr[0][k], is_2d)
            kp = _tavg(kp_arr[0][k], is_2d)
        except (TypeError, IndexError, KeyError):
            continue
        sg = None if sg is None else sg * fac
        kp = None if kp is None else kp * fac
        pf = (S * 1e-6) ** 2 * sg if (S is not None and sg is not None) else None
        rows.append({"T_K": T, "dataset": "overall", "S_uV/K": S, "sigma_S/m": sg,
                     "PF_W/mK2": pf, "kappa_e_W/mK": kp})
    return rows, None


def _rows_at_doping(cwd, is_2d, match, why_missing):
    """[0040] 本征载流子（AMSET doping = 0，E_F 由 n = p 定）的 σ/S/PF/κe 随温度：
    返回 (rows, note)；transport.json 里没有 0 掺杂点 -> ([], 原因)。
    2D：S8.4 的全部可用数据集（overall / 去 IMP / ADP only / 子集），σ、κe 按层厚 t 归一（× c/t）；
    3D：S8 的 transport.json（overall）。本征点 |n − p| ≈ 0，IMP 不起作用，overall ≈ 去 IMP。"""
    if is_2d:
        data = load_intrinsic_sets(cwd)
        geom = _geom_2d(cwd)
        if not data or not geom:
            return [], "没有 S8.4 transport.json / 2d_correction.json"
        doping, temps, fac = data["doping"], data["temps"], geom["c_over_t"]
        sets = [(k, data["sets"][k]) for k, *_ in data.get("set_defs", INTR_SETS) if k in data["sets"]]
    else:
        j = _load_json(Path(cwd) / AMSET_DIR / "transport.json")
        if not j:
            return [], "没有 S8 transport.json"
        doping = [float(x) for x in _unwrap(j["doping"])]
        temps = [float(t) for t in _unwrap(j["temperatures"])]
        fac = 1.0
        sets = [("overall", {"S": _unwrap(j.get("seebeck")), "sigma": _unwrap(j.get("conductivity")),
                             "kappa": _unwrap(j.get("electronic_thermal_conductivity"))})]
    i0 = next((i for i, x in enumerate(doping) if match(x)), None)
    if i0 is None:
        return [], why_missing
    rows = []
    for key, st in sets:
        if st.get("partial") or st.get("S") is None or st.get("sigma") is None:
            continue
        for k, T in enumerate(temps):
            try:
                S = _tavg(st["S"][i0][k], is_2d)
                sg = _tavg(st["sigma"][i0][k], is_2d)
                kp = _tavg(st["kappa"][i0][k], is_2d) if st.get("kappa") is not None else None
            except (TypeError, IndexError):
                continue
            sg = None if sg is None else sg * fac
            kp = None if kp is None else kp * fac
            pf = (S * 1e-6) ** 2 * sg if (S is not None and sg is not None) else None
            rows.append({"T_K": T, "dataset": key, "S_uV/K": S, "sigma_S/m": sg, "PF_W/mK2": pf,
                         "kappa_e_W/mK": kp})
    return rows, None


DEFECT_DIR_OUT = "defect_equilibrium"   # [0041] transport_plots/ 下缺陷平衡载流子图的子目录


def defect_rows(cwd, is_2d):
    """[0041] 缺陷平衡净浓度点（S8/S8.4 运行目录的 defect_doping.json 记的各化学势顶点）的 σ/S/PF/κe 随温度。
    返回 (rows, note)；每行多 "vertex"、"carrier_type"、"net"、"unit"。"""
    run = Path(cwd) / (AMSET2D_DIR if is_2d else AMSET_DIR)
    rec = _load_json(run / "defect_doping.json")
    if not rec or not rec.get("points"):
        return [], "运行目录里没有 defect_doping.json（DEFECT_DOPING 关了、缺陷技能没结果，或是 [0041] 之前的结果）"
    rows, miss = [], []
    for p in rec["points"]:
        v = float(p["amset_doping_cm3"])
        rr, why = _rows_at_doping(cwd, is_2d, lambda d, v=v: abs(d - v) <= 1e-6 * max(abs(v), 1.0),
                                  "transport.json 的掺杂网格里没有 %.3e（与 defect_doping.json 不是同一轮）" % v)
        if not rr:
            miss.append(why)
        for r in rr:
            r.update({"vertex": p.get("vertex"), "carrier_type": p.get("carrier_type"),
                      "net": p.get("net"), "unit": p.get("unit"), "E_F_eq_eV": p.get("E_F_eq_eV"),
                      "T_defect_K": p.get("T_K")})
            rows.append(r)
    return rows, ("; ".join(miss) if miss and not rows else None)


def write_defect_plots(cwd, plots, M, is_2d):
    """[0041] transport_plots/defect_equilibrium/<M>_defect_equilibrium.png + .csv：缺陷平衡净浓度下
    σ、S、PF、κe 随温度；每个化学势顶点一种颜色，实线 = 全部机制、虚线 = 去掉 IMP 的那套（若有）。"""
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ddir = Path(plots) / DEFECT_DIR_OUT
    rows, why = defect_rows(cwd, is_2d)
    if ddir.is_dir():
        for _old in ddir.glob("%s_*" % M):
            _old.unlink()
    if not rows:
        print("[..] 缺陷平衡载流子图跳过：%s" % why)
        return False
    ddir.mkdir(parents=True, exist_ok=True)
    if is_2d:
        data = load_intrinsic_sets(cwd) or {}
    else:
        _j = _load_json(Path(cwd) / AMSET_DIR / "transport.json") or {}
        data = {"overall_labels": [k for k in (_j.get("mobility") or {}) if str(k).lower() != "overall"]}
    verts = list(dict.fromkeys(r["vertex"] for r in rows))
    panels = [("sigma_S/m", "sigma (S/m)" + (", per layer t" if is_2d else ""), True),
              ("S_uV/K", "S (uV/K)", False), ("PF_W/mK2", "PF (W/m/K$^2$)", True),
              ("kappa_e_W/mK", "kappa_e (W/m/K)" + (", per layer t" if is_2d else ""), True)]
    fig, axs = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.8), squeeze=False)
    for j, (col, ttl, logy) in enumerate(panels):
        a = axs[0][j]
        for vi, vx in enumerate(verts):
            for ds, ls in (("overall", "-"), ("intrinsic", "--")):
                pts = sorted((r["T_K"], r[col]) for r in rows if r["vertex"] == vx and r["dataset"] == ds
                             and r[col] is not None and (not logy or r[col] > 0))
                if not pts:
                    continue
                r0 = next(r for r in rows if r["vertex"] == vx)
                lab = "vertex %s %s %.2e %s, %s" % (vx, r0.get("carrier_type") or "", float(r0["net"]),
                                                     r0.get("unit") or "", _intr_legend(ds, data, ds))
                a.plot([p[0] for p in pts], [p[1] for p in pts], ls=ls, marker="o", ms=3, lw=1.3,
                       color="C%d" % vi, label=lab)
        a.set_title(ttl, fontsize=10)
        a.set_xlabel("T (K)", fontsize=8)
        a.grid(alpha=0.25, lw=0.4)
        if logy:
            a.set_yscale("log")
        else:
            a.axhline(0, color="k", lw=0.5)
    h, l = axs[0][0].get_legend_handles_labels()
    if h:
        fig.legend(h, l, loc="lower center", ncol=2, fontsize=7.5, frameon=False)
    fig.suptitle("%s  defect-equilibrium carriers (defect-dft-cpu net p0-n0, fixed over T)  —  AMSET" % M, fontsize=11)
    fig.tight_layout(rect=[0, 0.08 + 0.035 * len(h) / 2, 1, 0.92])
    fig.savefig(ddir / ("%s_defect_equilibrium.png" % M), dpi=130)
    plt.close(fig)
    cols = ["vertex", "carrier_type", "net", "unit", "E_F_eq_eV", "T_defect_K", "T_K", "dataset",
            "sigma_S/m", "S_uV/K", "PF_W/mK2", "kappa_e_W/mK"]
    with open(ddir / ("%s_defect_equilibrium.csv" % M), "w", newline="") as fh:
        fh.write("# [0041] 缺陷平衡净浓度（defect-dft-cpu 各化学势顶点 p0-n0，在 T_defect_K 下求得）作为 AMSET 掺杂，各温度不变；"
                 "dataset 键：overall = phonon + impurity（全部机制），intrinsic = phonon-limited（去掉 IMP）；"
                 "IMP 按 |n-p| 计，未计补偿缺陷，偏低\n")
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([("%.6g" % r[c]) if isinstance(r.get(c), float) else ("" if r.get(c) is None else r[c])
                        for c in cols])
    print("[plots] %s/%s_defect_equilibrium.png" % (DEFECT_DIR_OUT, M))
    return True


def write_undoped_plots(cwd, plots, M, is_2d):
    """[0040] transport_plots/intrinsic_carrier/<M>_intrinsic_carrier.png + .csv：本征载流子（未掺杂）
    σ、S、PF、κe 随温度，每个数据集一条线。"""
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    udir = Path(plots) / UNDOPED_DIR
    rows, why = undoped_rows(cwd, is_2d)
    if udir.is_dir():
        for _old in udir.glob("%s_*" % M):
            _old.unlink()
    if not rows:
        print("[..] 本征载流子图跳过：%s" % why)
        return False
    udir.mkdir(parents=True, exist_ok=True)
    names = {"overall": "phonon + impurity", "intrinsic": "phonon-limited", "ADP_only": "ADP only (acoustic limit)"}
    if is_2d:
        data = load_intrinsic_sets(cwd)
    else:
        _j = _load_json(Path(cwd) / AMSET_DIR / "transport.json") or {}
        data = {"overall_labels": [k for k in (_j.get("mobility") or {}) if str(k).lower() != "overall"]}
    if data:
        for k in list(names) + [d[0] for d in data.get("set_defs", [])]:
            names[k] = _intr_legend(k, data, names.get(k, k))
    # [0054] 本征点 N_ii = |n − p| = 0：IMP 率为 0，含 IMP 的数据集与去掉 IMP 的同机制集逐点重合 ->
    #   按"去掉 IMP 后的机制组合"去重（优先留本来就不含 IMP 的那套），图例按去 IMP 后的机制命名（_phys_label）。
    _all_m = (data or {}).get("overall_labels")
    keep, eff_seen = {}, {}
    for k in dict.fromkeys(r["dataset"] for r in rows):
        ms = _intr_mechs(k, data) if data else None
        if not ms:
            keep[k] = names.get(k, k)                     # 机制未知（旧回退文件）：原样
            continue
        eff = tuple(m for m in ms if str(m).upper() != "IMP")
        if not eff:
            continue
        has_imp = len(eff) < len(ms)
        prev = eff_seen.get(eff)
        if prev is None or (prev[1] and not has_imp):
            if prev is not None:
                keep.pop(prev[0], None)
            eff_seen[eff] = (k, has_imp)
            keep[k] = (_phys_label(list(eff), _all_m) or "+".join(eff)) + (" [IMP = 0 at n = p]" if has_imp else "")
    rows = [r for r in rows if r["dataset"] in keep]
    names.update(keep)
    keys = list(keep)
    # [0052] 本征浓度 n_i(T)：transport.json 不带载流子浓度 -> 读 amset_fermi_fix 的记录；
    #   第一栏单独画 n_i(T) 并在每个点标数值，CSV 每行带 n_i / E_F。2D 另给面密度与文献口径（× c/t）。
    ni, ni_src = undoped_ni(cwd, is_2d)
    geom = _geom_2d(cwd) if is_2d else None
    for r in rows:
        v = _ni_lookup(ni, r["T_K"])
        r["n_i_cm-3"] = v.get("n_i_cm-3")
        r["E_F_eV"] = v.get("E_F_eV")
        if is_2d:
            ok = geom and r["n_i_cm-3"] is not None
            r["n_i_2D_cm-2"] = r["n_i_cm-3"] * geom["c_A"] * 1e-8 if ok else None
            r["n_i_lit_cm-3"] = r["n_i_cm-3"] * geom["c_over_t"] if ok else None
        r["E_v_eV"], r["E_c_eV"] = v.get("E_v_eV"), v.get("E_c_eV")
        r["eta_c"], r["eta_v"], r["f_CBM"] = _eta_of(v, r["T_K"])      # [0054] f ≈ 0 判据
        r["label"] = names.get(r["dataset"], r["dataset"])
    _eta_T = sorted({r["T_K"] for r in rows if r.get("eta_c") is not None})
    bad_T = sorted({r["T_K"] for r in rows if r.get("eta_c") is not None
                    and min(r["eta_c"], r["eta_v"]) < DILUTE_ETA_MIN})
    if bad_T:
        print("[WARN] 本征点在 %s 的 η = min(E_c−E_F, E_F−E_v)/k_BT < %g：带边占据不再 ≈ 0（窄带隙/高温），"
              "这些温度不是稀薄极限（图上灰底）" % (", ".join("%gK" % T for T in bad_T), DILUTE_ETA_MIN))
    elif not _eta_T:
        print("[..] 本征点没有带边记录（[0054] 之前的 undoped_ni.json）—— η / f(CBM) 留空；重跑 S8/S8.4 后可得")
    if ni:
        print("[..] 本征浓度 n_i(T) 读自 %s" % ni_src)
        _t_all = sorted({r["T_K"] for r in rows})
        _t_ni = sorted({r["T_K"] for r in rows if r.get("n_i_cm-3") is not None})
        if len(_t_ni) < len(_t_all):
            # [0052-r1] 别静默只给半套 n_i：日志回退只认到部分温度（本轮没写完 / 日志被截断）时说清楚
            print("[WARN] n_i 只覆盖 %d/%d 个温度点（缺 %s）—— 缺的点图上不标数值；"
                  "重跑 ke-dft-cpu S8/S8.4 会写全 undoped_ni.json"
                  % (len(_t_ni), len(_t_all),
                     ", ".join("%gK" % T for T in _t_all if T not in _t_ni)))
    else:
        print("[..] 本征浓度 n_i 无法标注：%s" % ni_src)
    panels = [("sigma_S/m", "sigma (S/m)" + (", per layer t" if is_2d else ""), True),
              ("S_uV/K", "S (uV/K)", False), ("PF_W/mK2", "PF (W/m/K$^2$)", True),
              ("kappa_e_W/mK", "kappa_e (W/m/K)" + (", per layer t" if is_2d else ""), True)]
    ni_col = "n_i_lit_cm-3" if is_2d else "n_i_cm-3"
    ni_pts = sorted({(r["T_K"], r[ni_col]) for r in rows if r.get(ni_col)})
    if ni_pts:
        panels.insert(0, (ni_col, ("n_i (cm$^{-3}$, per layer t)" if is_2d else "n_i = n = p (cm$^{-3}$)"), True))
    fig, axs = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.6), squeeze=False)
    for j, (col, ttl, logy) in enumerate(panels):
        a = axs[0][j]
        if col == ni_col:
            a.plot([p[0] for p in ni_pts], [p[1] for p in ni_pts], marker="o", ms=3, lw=1.3, color="k")
            for x, y in ni_pts:
                a.annotate("%.2g" % y, (x, y), textcoords="offset points", xytext=(0, 5),
                           ha="center", fontsize=6.5)
            a.set_title(ttl, fontsize=10)
            a.set_xlabel("T (K)", fontsize=8)
            a.grid(alpha=0.25, lw=0.4)
            a.set_yscale("log")
            a.margins(y=0.25)
            continue
        for ci, key in enumerate(keys):
            pts = sorted((r["T_K"], r[col]) for r in rows if r["dataset"] == key and r[col] is not None
                         and (not logy or r[col] > 0))
            if len(pts) >= 1:
                a.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", ms=3, lw=1.3,
                       color="C%d" % ci, label=names.get(key, key))
        a.set_title(ttl, fontsize=10)
        a.set_xlabel("T (K)", fontsize=8)
        a.grid(alpha=0.25, lw=0.4)
        if logy:
            a.set_yscale("log")
        else:
            a.axhline(0, color="k", lw=0.5)
    for a in axs[0]:
        _shade_not_dilute(a, bad_T, [r["T_K"] for r in rows])
    h, l = max((a.get_legend_handles_labels() for a in axs[0]), key=lambda hl: len(hl[0]))
    if h:
        fig.legend(h, l, loc="lower center", ncol=min(3, len(h)), fontsize=8, frameon=False)
    fig.suptitle("%s  intrinsic carriers (undoped, n = p = n_i; E_F near mid-gap, IMP = 0; S, sigma bipolar)  —  AMSET%s%s"
                 % (M, "" if ni_pts else "   [n_i not recorded: rerun S8 to get it]",
                    ("\nshaded: eta = (E_c - E_F)/kT < %g, band-edge occupation not ~ 0 (not the dilute limit)"
                     % DILUTE_ETA_MIN) if bad_T else ""), fontsize=11)
    fig.tight_layout(rect=[0, 0.1, 1, 0.92])
    fig.savefig(udir / ("%s_intrinsic_carrier.png" % M), dpi=130)
    plt.close(fig)
    cols = ["T_K", "n_i_cm-3"] + (["n_i_2D_cm-2", "n_i_lit_cm-3"] if is_2d else []) + [
        "E_F_eV", "E_v_eV", "E_c_eV", "eta_c", "eta_v", "f_CBM",
        "dataset", "label", "sigma_S/m", "S_uV/K", "PF_W/mK2", "kappa_e_W/mK"]
    with open(udir / ("%s_intrinsic_carrier.csv" % M), "w", newline="") as fh:
        fh.write("# [0040][0054] 本征载流子（AMSET doping = 0，E_F 由 n = p 定）；%s；N_ii = |n-p| = 0，IMP 为零，"
                 "含 IMP 的数据集与去 IMP 的同机制集重合、已去重（label 列按去 IMP 后的机制命名）；"
                 "σ、S 是双极量（电子+空穴）；eta_c/eta_v = (E_c-E_F)/kT、(E_F-E_v)/kT，< %g 即不是稀薄极限 f≈0\n"
                 % ("2D：sigma、kappa_e、PF 按层厚 t 归一（元胞值 × c/t）" if is_2d else "3D：张量对角平均", DILUTE_ETA_MIN))
        fh.write("# [0052] n_i = n = p（%s）：n_i_cm-3 为每胞体积口径%s；来源 %s\n"
                 % ("本征浓度", "（含真空元胞）；n_i_2D_cm-2 = n_i × c，n_i_lit_cm-3 = n_i × c/t" if is_2d else "",
                    ni_src if ni else "无（%s）" % ni_src))
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([("%.6g" % r[c]) if isinstance(r.get(c), float) else ("" if r.get(c) is None else r[c])
                        for c in cols])
    print("[plots] %s/%s_intrinsic_carrier.png" % (UNDOPED_DIR, M))
    return True


# =====================================================================
# [0053][0054] phonon-limited 迁移率按散射机制拆开（mechanism-resolved，稀薄极限 f ≈ 0；2026-10-16 用户要求；对照 Cheng & Liu, JACS 2018, 140, 17895
#   的分模式迁移率 + Matthiessen）。上面的 <M>_intrinsic_carrier 画的是 n = p 时的 σ/S/PF/κe —— 双极（电子+空穴
#   一起）的量，不是单载流子迁移率；本节只取稀薄极限这一个浓度，画 μ(T) 的机制分解。
#   · 空穴：本征点 doping = 0。AMSET 的 _calculate_mobility（_get_band_idx）对 doping ≥ 0 只积价带、除以 hole_conc，
#     所以 0 点给的就是空穴迁移率；本征点 IMP 的 N_ii = |n − p| = 0、屏蔽来自 n_i（可忽略）、E_F 在带隙中部（f ≈ 0）。
#   · 电子：AMSET 在 0 点拿不到 -> 取 |掺杂| 最低的 n 型点，总值用去 IMP 的那套（IMP 随掺杂变，不属于稀薄极限）。
#   · 核对：空穴另取最低 p 型点（去 IMP）、电子取次低 n 型点（去 IMP），与主值差 > 5% 说明还没到非简并平台
#     （或自由载流子屏蔽已起作用），打 WARN。
#   · 口径差异（写进图注）：AMSET 是 RTA 对 −∂f/∂E 的热平均，单机制 μ 的倒数和 ≥ 总值的倒数（Matthiessen 只给上界）；
#     Cheng & Liu 取带边单态 μ = eτ/m*，Matthiessen 严格成立。两者差在热平均，不在占据数。
# =====================================================================
DILUTE_MU_HUGE = 1e10      # cm²/Vs：单机制 μ 超过它 = 该机制前因子为 0（本征点 IMP、非极性材料 POP），不画、不进 Matthiessen
DILUTE_CHECK_TOL = 0.05    # 主值与核对点相对差超过它打 WARN
_DILUTE_COLORS = {"ADP": "C2", "POP": "C1", "ODP": "C6", "IVS": "C7", "PIE": "C9", "IMP": "C3"}
TOTAL_NOIMP, MATTHIESSEN = "phonon-limited total", "Matthiessen sum"   # [0054] 统一命名


def mechanism_resolved_rows(cwd, is_2d):
    """[0053] 稀薄极限迁移率行：[{T_K, carrier, role, source, doping_cm3, n_2D_cm-2, mechanism, mu_cm2/Vs,
    share_1/mu, eta}]；role = main（画实线）/ check（核对点）。没有 AMSET transport.json -> ([], 原因)。
    [0054] n = p 点带 eta = min(η_c, η_v)（undoped_ni.json 有带边时），< DILUTE_ETA_MIN 即不是 f ≈ 0。"""
    run = AMSET2D_DIR if is_2d else AMSET_DIR
    j = _load_json(Path(cwd) / run / "transport.json")
    if not j:
        return [], "没有 %s/transport.json" % run
    data = load_intrinsic_sets(cwd, run_name=run) or {}
    doping = [float(x) for x in _unwrap(j["doping"])]
    temps = [float(t) for t in _unwrap(j["temperatures"])]
    mob = j.get("mobility") or {}
    singles = {str(k).upper(): _unwrap(v) for k, v in mob.items()
               if str(k).lower() != "overall" and v is not None}
    labels = set(singles)
    sets = data.get("sets") or {}
    intr = sets.get("intrinsic") or {}
    tot_arr = intr.get("mu") if not intr.get("partial") else None
    tot_src = "去 IMP（%s）" % "+".join(intr.get("labels") or []) if tot_arr is not None else None
    if tot_arr is None and "IMP" not in labels:
        tot_arr, tot_src = _unwrap(mob.get("overall") or mob.get("Overall")), "overall（本轮没开 IMP）"
    c_cm = (_geom_2d(cwd) or {}).get("c_A", 0) * 1e-8 if is_2d else None
    ni, _ni_src = undoped_ni(cwd, is_2d)

    def _lowest(sign, nth=0):
        idx = sorted((i for i, x in enumerate(doping) if (x < 0 if sign < 0 else x > 0)), key=lambda i: abs(doping[i]))
        return idx[nth] if len(idx) > nth else None

    i0 = next((i for i, x in enumerate(doping) if abs(x) <= 1e-12), None)
    plan = []                                   # (carrier, role, 掺杂索引, 来源说明, 总值数组, 总值说明)
    e_main, e_chk = _lowest(-1), _lowest(-1, 1)
    if e_main is not None:
        plan.append(("electron", "main", e_main, "lowest n-type", tot_arr, tot_src))
    if e_chk is not None and tot_arr is not None:
        plan.append(("electron", "check", e_chk, "2nd-lowest n-type", tot_arr, tot_src))
    h_low = _lowest(+1)
    if i0 is not None:
        # 本征点 |n − p| = 0：IMP 率为 0，overall 本身就是去 IMP 的（有子集用子集，没有用 overall）
        h_tot = tot_arr if tot_arr is not None else _unwrap(mob.get("overall") or mob.get("Overall"))
        plan.append(("hole", "main", i0, "n = p (undoped)", h_tot, tot_src or "overall（N_ii = 0）"))
        if h_low is not None and tot_arr is not None:
            plan.append(("hole", "check", h_low, "lowest p-type", tot_arr, tot_src))
    elif h_low is not None:
        plan.append(("hole", "main", h_low, "lowest p-type (no n = p point)", tot_arr, tot_src))
    if not plan:
        return [], "transport.json 的掺杂网格里没有可用的 n 型 / p 型 / 0 点"
    rows = []
    for car, role, i, src, tarr, tsrc in plan:
        x = doping[i]
        n2d = abs(x) * c_cm if c_cm else None
        for k, T in enumerate(temps):
            base = {"T_K": T, "carrier": car, "role": role, "source": src, "doping_cm3": x, "n_2D_cm-2": n2d,
                    "eta": None}
            if src.startswith("n = p"):
                _ec, _ev, _f = _eta_of(_ni_lookup(ni, T), T)
                base["eta"] = min(_ec, _ev) if _ec is not None else None
            mech = {}
            if role == "main":
                for name in sorted(singles):
                    if name == "IMP":
                        continue
                    try:
                        v = _tavg(singles[name][i][k], is_2d)
                    except (TypeError, IndexError):
                        v = None
                    if v is not None and 0 < v < DILUTE_MU_HUGE:
                        mech[name] = v
                inv = sum(1.0 / v for v in mech.values())
                for name, v in mech.items():
                    rows.append(dict(base, mechanism=name, **{"mu_cm2/Vs": v, "share_1/mu": (1.0 / v) / inv}))
                if mech:
                    rows.append(dict(base, mechanism=MATTHIESSEN, **{"mu_cm2/Vs": 1.0 / inv, "share_1/mu": None}))
            try:
                tv = _tavg(tarr[i][k], is_2d) if tarr is not None else None
            except (TypeError, IndexError):
                tv = None
            if tv is not None and 0 < tv < DILUTE_MU_HUGE:
                rows.append(dict(base, mechanism=TOTAL_NOIMP, total_source=tsrc,
                                 **{"mu_cm2/Vs": tv, "share_1/mu": None}))
    return rows, None


MECH_RES_NAME = "mechanism_resolved_mobility"   # [0054] 原 dilute_mobility（文献：mode-/mechanism-resolved mobility）


def write_mechanism_resolved_plot(cwd, plots, M, is_2d):
    """[0053][0054] transport_plots/intrinsic_carrier/<M>_mechanism_resolved_mobility.png + .csv：
    稀薄极限 phonon-limited μ(T) 的机制分解。"""
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    udir = Path(plots) / UNDOPED_DIR
    for _pat in ("%s_%s.*" % (M, MECH_RES_NAME), "%s_dilute_mobility.*" % M):   # 连 [0053] 旧名一起清
        for _old in (udir.glob(_pat) if udir.is_dir() else []):
            _old.unlink()
    rows, why = mechanism_resolved_rows(cwd, is_2d)
    if not rows:
        print("[..] 稀薄极限迁移率图跳过：%s" % why)
        return False
    udir.mkdir(parents=True, exist_ok=True)

    def _mu(car, role, mech, T):
        return next((r["mu_cm2/Vs"] for r in rows if r["carrier"] == car and r["role"] == role
                     and r["mechanism"] == mech and r["T_K"] == T), None)

    temps = sorted({r["T_K"] for r in rows})
    T_rep = min(temps, key=lambda t: abs(t - 300.0))
    cars = [c for c in ("electron", "hole") if any(r["carrier"] == c and r["role"] == "main" for r in rows)]
    for car in cars:                                       # 核对：主值 vs 核对点（去 IMP 总值）
        for T in temps:
            a, b = _mu(car, "main", TOTAL_NOIMP, T), _mu(car, "check", TOTAL_NOIMP, T)
            if a and b and abs(a / b - 1.0) > DILUTE_CHECK_TOL:
                src = next(r["source"] for r in rows if r["carrier"] == car and r["role"] == "check")
                print("[WARN] %s 稀薄极限 μ 在 %g K 与 %s 差 %.0f%%（%.4g vs %.4g cm²/Vs）—— 还没到非简并平台，"
                      "或自由载流子屏蔽已起作用" % (car, T, src, 100 * (a / b - 1.0), a, b))
            m_ = _mu(car, "main", MATTHIESSEN, T)
            if a and m_ and a > m_ * 1.01:
                # RTA 下 <1/Σ(1/τ_i)> ≤ 1/Σ(1/<τ_i>)（调和平均凹性），总值高于倒数和说明单机制与去 IMP 子集口径不一致
                print("[WARN] %s 稀薄极限 %g K：去 IMP 总值 %.4g > Matthiessen 倒数和 %.4g cm²/Vs —— RTA 下不应发生，"
                      "核对单机制 μ 与子集是否同一轮、同一屏蔽设置" % (car, T, a, m_))
        shares = sorted(((r["share_1/mu"], r["mechanism"]) for r in rows if r["carrier"] == car
                         and r["role"] == "main" and r["T_K"] == T_rep and r["share_1/mu"] is not None), reverse=True)
        if shares:
            print("[..] 稀薄极限 %s @%gK 限制机制：%s" % (
                car, T_rep, "，".join("%s %.0f%%" % (m, 100 * s) for s, m in shares)))
    fig, axs = plt.subplots(1, len(cars), figsize=(5.6 * len(cars), 4.4), squeeze=False)
    for ax, car in zip(axs[0], cars):
        rc = [r for r in rows if r["carrier"] == car]
        mechs = list(dict.fromkeys(r["mechanism"] for r in rc if r["role"] == "main"
                                   and r["mechanism"] not in (TOTAL_NOIMP, MATTHIESSEN)))
        for ci, m in enumerate(mechs):
            pts = sorted((r["T_K"], r["mu_cm2/Vs"]) for r in rc if r["role"] == "main" and r["mechanism"] == m)
            ax.plot(*zip(*pts), ls="--", marker="o", ms=3, lw=1.2,
                    color=_DILUTE_COLORS.get(m, "C%d" % (ci + 3)), label="%s only" % m)
        _bad = sorted({r["T_K"] for r in rc if r.get("eta") is not None and r["eta"] < DILUTE_ETA_MIN})
        if _bad:
            _shade_not_dilute(ax, _bad, [r["T_K"] for r in rc])
            print("[WARN] %s（n = p 点）在 %s 的 η < %g：不是稀薄极限，与文献 f ≈ 0 的值不可直接比（图上灰底）"
                  % (car, ", ".join("%gK" % T for T in _bad), DILUTE_ETA_MIN))
        for m, kw in ((TOTAL_NOIMP, dict(color="k", ls="-", lw=2.0, marker="s", ms=3.5)),
                      (MATTHIESSEN, dict(color="0.45", ls=":", lw=1.6))):
            pts = sorted((r["T_K"], r["mu_cm2/Vs"]) for r in rc if r["role"] == "main" and r["mechanism"] == m)
            if pts:
                ax.plot(*zip(*pts), label=("Matthiessen $(\\sum 1/\\mu_i)^{-1}$" if m == MATTHIESSEN
                                           else "phonon-limited total (no IMP)"), **kw)
        chk = sorted((r["T_K"], r["mu_cm2/Vs"]) for r in rc if r["role"] == "check")
        if chk:
            src = next(r["source"] for r in rc if r["role"] == "check")
            ax.plot(*zip(*chk), ls="none", marker="o", ms=7, mfc="none", mec="k", mew=0.8,
                    label="check: total @ %s" % src)
        r0 = next(r for r in rc if r["role"] == "main")
        if r0["source"].startswith("n = p"):
            where = "n = p (undoped, E_F mid-gap)"
        elif r0.get("n_2D_cm-2"):
            where = "%s, n$_{2D}$ = %.1e cm$^{-2}$" % (r0["source"], r0["n_2D_cm-2"])
        else:
            where = "%s, |doping| = %.1e cm$^{-3}$" % (r0["source"], abs(r0["doping_cm3"]))
        ax.set_title("%s  (%s)" % (car, where), fontsize=9.5)
        ax.set_yscale("log")
        ax.set_xlabel("T (K)")
        ax.set_ylabel("mobility (cm$^2$/Vs)" + (", in-plane" if is_2d else ""))
        ax.grid(alpha=0.25, lw=0.4)
        ax.legend(fontsize=7, loc="best", framealpha=0.85)
    fig.suptitle("%s  mechanism-resolved phonon-limited mobility, dilute limit (f ≈ 0) — AMSET RTA "
                 "(thermal average; Matthiessen sum is an upper bound; shaded: eta < %g, not dilute)"
                 % (M, DILUTE_ETA_MIN), fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(udir / ("%s_%s.png" % (M, MECH_RES_NAME)), dpi=130)
    plt.close(fig)
    cols = ["T_K", "carrier", "role", "source", "doping_cm3"] + (["n_2D_cm-2"] if is_2d else []) + [
        "mechanism", "mu_cm2/Vs", "share_1/mu", "eta", "total_source"]
    with open(udir / ("%s_%s.csv" % (M, MECH_RES_NAME)), "w", newline="") as fh:
        fh.write("# [0053][0054] mechanism-resolved phonon-limited 迁移率（稀薄极限 f ≈ 0）；eta = n = p 点的 "
                 "min(E_c-E_F, E_F-E_v)/kT（< %g 不是稀薄极限）；" % DILUTE_ETA_MIN)
        fh.write("空穴取本征点 doping = 0（AMSET 对 doping ≥ 0 只积价带、"
                 "除以 hole_conc），电子取最低 n 型点（AMSET 在 0 点不给电子 μ）；总值去 IMP\n")
        fh.write("# 单机制 = AMSET transport.json 的 mobility[机制]（只开该机制）；share_1/mu = (1/μ_i)/Σ(1/μ_i)；"
                 "Matthiessen sum = (Σ 1/μ_i)^-1 ≥ 总值（RTA 热平均）；μ > %.0e 的机制（前因子 0）不列；"
                 "role = check 是非简并平台核对点\n" % DILUTE_MU_HUGE)
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([("%.6g" % r[c]) if isinstance(r.get(c), float) else ("" if r.get(c) is None else r[c])
                        for c in cols])
    print("[plots] %s/%s_%s.png" % (UNDOPED_DIR, M, MECH_RES_NAME))
    return True


def _clear_dir(d):
    """[0052] 清空目录内容（保留目录本身）；删不掉的只告警。
    [0052-r1] 目录本身是符号链接时不动手：否则 rmtree/unlink 会顺着链接删到别处的产物。"""
    import shutil
    p = Path(d)
    if p.is_symlink():
        print("[WARN] %s 是符号链接，跳过清空（避免顺着链接删到别处）" % d)
        return 0
    n = 0
    for f in list(p.iterdir()) if p.is_dir() else []:
        try:
            if f.is_dir() and not f.is_symlink():
                shutil.rmtree(f)
            else:
                f.unlink()
            n += 1
        except OSError as e:
            print("[WARN] 删旧产物 %s 失败：%s" % (f, e))
    if n:
        print("[plots] 已删除 %s/ 下上一轮的 %d 项旧结果" % (Path(d).name, n))
    return n


def write_transport_plots(cwd, out, is_2d):
    """[0033]-[0035] 规定格式绘图输出（图 + 同目录数据）到【材料根目录】的 transport_plots/：
      <M>_transport.png                        9温度(行)x5参数(列 mu/S/sigma/PF/kappa_e)；
                                               基础数据集：overall / 去 IMP / ADP only / BT2xDPT-tau；按可用结果加入机制对照
      <M>_transport_data_<T>K.csv               对应数据，每温度一个
      <M>_params_<T>K.csv                       每温度一张(E1/C2D/m*/tau/mu)
      datasets/<M>_<数据集>.png/.csv             [0039] 每个数据集（overall / intrinsic / ADP_only / 子集 / BT2）
                                               单独一张 σ/S/PF/κe 随浓度图，每温度一条线
    [0035] 原 <M>_intrinsic.png 与 <M>_transport.png 数据集重叠（overall/ADP/DPT 三套重复），
      已合并为一张 <M>_transport.png；额外机制对照按 S8.4 可用结果动态显示。
    [0047] 2026-10-15 用户选 2：<M>_intrinsic.png 恢复为【只画基础 4 条】(overall / 去 IMP /
      ADP only / BT2×DPT-τ)，与带机制对照的 <M>_transport.png 内容真正区分；
      _intrinsic_data_<T>K.csv 同步只含这 4 条（_transport_data_<T>K.csv 仍是全数据集）。
    [0034] 输出位置 = cwd.parent/transport_plots/（材料根，不在 step8.3_output 内）；
      skill.yaml 的 S8.3 用 fetch_extra 让 autozt fetch 一并把它拉回本地。
      绝不写 out 根目录的 comparison_300K.*。失败不阻断 S8.3 完成标记。"""
    import csv
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plots = cwd.parent / "transport_plots"          # [0034] 材料根目录（cwd = <材料>/ke-dft-cpu）
    try:
        plots.mkdir(parents=True, exist_ok=True)
    except OSError:
        plots = out / "plots"                       # 回退：材料根不可写时仍写步骤目录
        plots.mkdir(exist_ok=True)
    M = cwd.parent.name
    NFLOOR = 5e17
    # [0052] 先删旧结果：transport_plots/ 是 S8.3 独占的输出目录，每次整目录重建。以前只覆盖同名文件，
    #   温度网格 / 数据集 / 子文件夹变了以后，上一轮的 <T>K.csv、旧格式长表、旧子目录会一直留着冒充当前结果。
    _clear_dir(plots)

    def fnum(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    def xs(n):
        n = np.asarray(n, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.sign(n) * np.log10(np.maximum(np.abs(n), 1e17) / 1e17)

    def seg(ax, n, y, **kw):
        n = np.asarray(n, float); y = np.asarray(y, float)
        for msk in (n < 0, n > 0):
            mm = msk & np.isfinite(y) & np.isfinite(n)
            if mm.sum() < 2:
                continue
            o = np.argsort(xs(n[mm]))
            ax.plot(xs(n[mm])[o], y[mm][o], **kw)

    def grid(nrow, ncol, rowlab, collab):
        fig, axes = plt.subplots(nrow, ncol, figsize=(19.5, 2.35 * nrow), squeeze=False)
        for rr in range(nrow):
            for j in range(ncol):
                ax = axes[rr][j]
                ax.set_xlim(-4.15, 4.15)
                ax.set_xticks([-4, -3, -2, -1, 0, 1, 2, 3, 4])
                ax.set_xticklabels(["1e21", "1e20", "1e19", "1e18", "1e17",
                                    "1e18", "1e19", "1e20", "1e21"], fontsize=6.5)
                ax.axvline(0, color="0.55", lw=0.6, ls=":")
                ax.grid(alpha=0.25, which="both", lw=0.4)
                if rr == 0:
                    ax.set_title(collab[j], fontsize=10)
                if j == 0:
                    ax.set_ylabel(rowlab[rr], fontsize=9, fontweight="bold")
                if rr == nrow - 1:
                    ax.set_xlabel(r"doping $n$ (cm$^{-3}$):  n-type $\leftarrow$ | $\rightarrow$ p-type",
                                  fontsize=7)
        return fig, axes

    # ---- 规定格式主图：从 S8.3 本征长表读（数据集最全）----
    #   [0035] 原 <M>_intrinsic.png 并入 <M>_transport.png（两者数据集重叠）；
    #   基础曲线为 overall / 去 IMP / ADP only / BT2xDPT-tau；另外显示 S8.4 可用的机制组合。
    itab = out / "intrinsic" / "intrinsic_transport_table.csv"
    if itab.is_file():
        rows = list(csv.DictReader([l for l in itab.read_text(errors="ignore").splitlines()
                                    if l.strip() and not l.startswith("#")]))
        temps = sorted(set(float(r["T_K"]) for r in rows))
        cols = ["T_K", "dataset", "type", "n_cell_cm-3", "n_2D_cm-2", "n_lit_cm-3",
                "carriers_per_cell", "carriers_per_atom", "in_interpret_range",
                "S_uV/K", "sigma_S/m", "PF_W/mK2", "kappa_e_W/mK", "mu_cm2/Vs", "flag"]
        _base_names = {"overall", "intrinsic", "ADP_only", "BT2_DPTtau"}   # [选 2] phonon_limited 产物只含这 4 条
        for T in temps:
            # [用户 2026-10-14 选 B / 2026-10-15 选 2] 独立 intrinsic 产物名：
            #   transport CSV = 全数据集；intrinsic CSV = 仅基础 4 条（与 intrinsic 图一致）
            for _sfx in ("transport", "phonon_limited"):            # [0054] 原 intrinsic -> phonon_limited
                with open(plots / ("%s_%s_data_%dK.csv" % (M, _sfx, int(round(T)))), "w", newline="") as fh:
                    w = csv.writer(fh); w.writerow(cols)
                    for r in rows:
                        if abs(float(r["T_K"]) - T) > 1e-6:
                            continue
                        if _sfx == "phonon_limited" and r["dataset"] not in _base_names:
                            continue
                        w.writerow([r.get(c, "") for c in cols])
        quants = [("mu_cm2/Vs", "mu  (cm^2/V/s)"), ("S_uV/K", "S  (uV/K)"), ("sigma_S/m", "sigma  (S/m)"),
                  ("PF_W/mK2", "PF  (W/m/K$^2$)"), ("kappa_e_W/mK", "kappa_e  (W/m/K)")]
        dsets = [("overall", "tab:red", "--", "s"), ("intrinsic", "tab:orange", "-.", "D"),
                 ("ADP_only", "tab:green", ":", "^"), ("BT2_DPTtau", "tab:blue", "-", "o")]
        _base_dsets = list(dsets)   # [用户 2026-10-15 选 2] 本征专用图只画这 4 条基础曲线
        _LBL = {"overall": "phonon + impurity",                  # [0054] 统一命名（_phys_label）
                "intrinsic": "phonon-limited",
                "ADP_only": "ADP only (acoustic limit)",
                "ADP_POP": "partial phonon (ADP+POP)",
                "ADP_POP_ODP": "partial phonon (ADP+POP+ODP)",
                "ODP": "ODP only (single-mech, mu only)", "IVS": "IVS only (single-mech, mu only)",
                "BT2_DPTtau": "BT2×DPT-τ"}
        _set_meta = load_intrinsic_sets(cwd, AMSET2D_DIR if is_2d else AMSET_DIR) or {}
        _all_mechs = _set_meta.get("overall_labels") or []
        _im = ((_set_meta.get("sets") or {}).get("intrinsic") or {}).get("labels") or []
        if _all_mechs:
            _LBL["overall"] = _phys_label(_all_mechs, _all_mechs)
        if _im:
            _LBL["intrinsic"] = _phys_label(_im, _all_mechs)
        for _k in ("ADP_POP", "ADP_POP_ODP"):
            _LBL[_k] = _intr_legend(_k, _set_meta, _LBL[_k])
        _have = {r["dataset"] for r in rows}                       # [0035] 动态加 ODP/IVS：长表里有才画
        for _k, _c, _ls, _mk in (("ADP_POP", "tab:brown", "-", "v"), ("ADP_POP_ODP", "tab:olive", "-", "P"),
                                 ("ODP", "tab:pink", ":", "P"), ("IVS", "tab:gray", "--", "X")):
            if _k in _have:
                dsets.append((_k, _c, _ls, _mk))
        def _render_9x5(dsets_used, title, outpath):
            """[用户 2026-10-15 选 2] 同一套 9x5 网格渲染器：transport 图带本轮机制对照、
            intrinsic 图只画基础 4 条 —— 两张图内容真正不同（此前是同一张 fig 换标题存两次）。"""
            fig, axes = grid(len(temps), 5, ["%d K" % t for t in temps], [q[1] for q in quants])
            for rr, T in enumerate(temps):
                sub = [r for r in rows if abs(float(r["T_K"]) - T) < 1e-6]
                for jc, (col, _) in enumerate(quants):
                    ax = axes[rr][jc]
                    for ds, color, ls, mk in dsets_used:
                        dd = [r for r in sub if r["dataset"] == ds]
                        if not dd:
                            continue
                        n = np.array([fnum(r["n_cell_cm-3"]) for r in dd])
                        y = np.array([fnum(r[col]) for r in dd])
                        keep = np.abs(n) >= NFLOOR
                        seg(ax, n[keep], y[keep], color=color, lw=1.3, ls=ls, marker=mk, ms=3.2)
                    if col == "S_uV/K":
                        ax.set_ylim(-1000, 1000); ax.axhline(0, color="k", lw=0.5)
                    else:
                        ax.set_yscale("log")
            hh = [Line2D([], [], color=c, lw=1.4, ls=ls, marker=mk, ms=4, label=_LBL.get(ds, ds))
                  for ds, c, ls, mk in dsets_used]
            fig.legend(handles=hh, loc="lower center", ncol=3, fontsize=9, frameon=False)
            fig.suptitle(title, fontsize=12)
            fig.tight_layout(rect=[0, 0.03, 1, 0.985 - 0.009 * title.count("\n")])
            fig.savefig(outpath, dpi=130)
            plt.close(fig)

        _trng = "%g-%g K" % (temps[0], temps[-1]) if temps else ""
        _render_9x5(dsets,
                    "%s  transport vs doping  (AMSET mechanisms from this run; BT2 x DPT tau)  %s" % (M, _trng),
                    plots / ("%s_transport.png" % M))
        # [0052] 用户 2026-10-15 问「本征浓度横坐标不是定值吗」；[0054] 2026-10-16 用户定：与文献统一命名，
        #   原 <M>_intrinsic.png（去 IMP 的掺杂扫描）改名 <M>_phonon_limited.png；"intrinsic" 只指本征载流子
        #   （未掺杂，n = p = n_i，intrinsic_carrier/）。
        _render_9x5(_base_dsets,
                    "%s  phonon-limited transport vs doping, base 4 sets  (phonon + impurity / phonon-limited"
                    " / ADP only / BT2xDPT-tau)  %s\n"
                    "x = doping concentration; phonon-limited = all phonon mechanisms, no ionized impurity"
                    " (intrinsic carriers n = p = n_i -> see intrinsic_carrier/)" % (M, _trng),
                    plots / ("%s_phonon_limited.png" % M))
        print("[plots] %s_transport.png (基础 4 条 + 机制对照) + %s_phonon_limited.png (仅基础 4 条) / "
              "_transport|_phonon_limited_data_<T>K.csv x%d" % (M, M, len(temps)))

        # ---- [0039] 每个数据集单独一张 σ/S/PF/κe（随浓度，每温度一条线）：datasets/<M>_<数据集>.png + .csv ----
        #   主图把 overall / intrinsic / ADP only / BT2 挤在一张 9x5 里，intrinsic（去 IMP）等子集没有单独的图。
        try:
            _ddir = plots / "datasets"
            _ddir.mkdir(exist_ok=True)
            for _old in _ddir.glob("%s_*" % M):                 # 清旧产物（数据集可能变了）
                _old.unlink()
            _dq = [("sigma_S/m", "sigma  (S/m)"), ("S_uV/K", "S  (uV/K)"),
                   ("PF_W/mK2", "PF  (W/m/K$^2$)"), ("kappa_e_W/mK", "kappa_e  (W/m/K)")]
            _cm = plt.get_cmap("viridis")
            _dcols = ["T_K", "type", "n_cell_cm-3", "n_2D_cm-2", "carriers_per_atom", "in_interpret_range",
                      "sigma_S/m", "S_uV/K", "PF_W/mK2", "kappa_e_W/mK", "mu_cm2/Vs", "flag"]
            for ds, _c, _ls, _mk in dsets:
                dd_all = [r for r in rows if r["dataset"] == ds]
                if not dd_all:
                    continue
                _f, _axs = plt.subplots(1, len(_dq), figsize=(4.2 * len(_dq), 3.6), squeeze=False)
                for jc, (col, ttl) in enumerate(_dq):
                    a = _axs[0][jc]
                    for it, T in enumerate(temps):
                        dd = [r for r in dd_all if abs(float(r["T_K"]) - T) < 1e-6]
                        n = np.array([fnum(r["n_cell_cm-3"]) for r in dd])
                        y = np.array([fnum(r[col]) for r in dd])
                        keep = np.abs(n) >= NFLOOR
                        seg(a, n[keep], y[keep], color=_cm(it / max(1, len(temps) - 1)), lw=1.2,
                            marker="o", ms=2.6, label=("%d K" % T))
                    a.set_title(ttl, fontsize=10)
                    a.set_xlim(-4.15, 4.15); a.set_xticks([-4, -3, -2, -1, 0, 1, 2, 3, 4])
                    a.set_xticklabels(["1e21", "1e20", "1e19", "1e18", "1e17",
                                       "1e18", "1e19", "1e20", "1e21"], fontsize=6.5)
                    a.axvline(0, color="0.55", lw=0.6, ls=":")
                    a.grid(alpha=0.25, lw=0.4)
                    a.set_xlabel(r"$n$ (cm$^{-3}$)   n-type | p-type", fontsize=8)
                    if col == "S_uV/K":
                        a.set_ylim(-1000, 1000); a.axhline(0, color="k", lw=0.5)
                    else:
                        a.set_yscale("log")
                _hh = [Line2D([], [], color=_cm(it / max(1, len(temps) - 1)), lw=1.4, label="%d K" % T)
                       for it, T in enumerate(temps)]
                _f.legend(handles=_hh, loc="lower center", ncol=min(9, len(_hh)), fontsize=8, frameon=False)
                _f.suptitle("%s  %s" % (M, _LBL.get(ds, ds)), fontsize=11)
                _f.tight_layout(rect=[0, 0.08, 1, 0.93])
                _fn = _ds_file(ds, _set_meta)       # [0040] 文件名 = 实际机制组合（如 ADP+POP+ODP）
                _f.savefig(_ddir / ("%s_%s.png" % (M, _fn)), dpi=130); plt.close(_f)
                with open(_ddir / ("%s_%s.csv" % (M, _fn)), "w", newline="") as fh:
                    w = csv.writer(fh)
                    w.writerow(["# dataset=%s  label=%s" % (ds, _LBL.get(ds, ds))])
                    w.writerow(_dcols)
                    for r in dd_all:
                        w.writerow([r.get(c, "") for c in _dcols])
                print("[plots] datasets/%s_%s.png" % (M, _fn))
        except Exception as _e:                                   # noqa: BLE001
            print("[WARN] 分数据集图失败（%s: %s）" % (type(_e).__name__, _e))

    # ---- ③ 参数表（每温度一张，文献格式）----
    vp = out / MU_T_CSV
    if vp.is_file():
        R = list(csv.DictReader([l for l in vp.read_text(errors="ignore").splitlines()
                                 if l.strip() and not l.startswith("#")]))
        temps = sorted(set(int(float(r["T_K"])) for r in R))
        for T in temps:
            sub = [r for r in R if int(float(r["T_K"])) == T]
            with open(plots / ("%s_params_%dK.csv" % (M, T)), "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Type", "Carrier", "Direction", "E1 (eV)", ("C2D (N/m)" if is_2d else "C3D (GPa)"), "m* (me)", "mu (cm2/V/s)", "tau (1e-15 s)"])
                for c in ("electron", "hole"):
                    for r in [x for x in sub if x["carrier"] == c]:
                        w.writerow([M, "N" if c == "electron" else "P", r["direction"], r["E1_eV"],
                                    r[_c_col(is_2d)], r["m_eff_m0"], r["mu_DPT"], r["tau_DPT_fs"]])
        print("[plots] %s_params_<T>K.csv x%d" % (M, len(temps)))

    # ---- ⑥ bt2 μ 扫描图 + amset 浓度扫描图（各放一个子文件夹）----
    #   [0036] bt2 的 paper_mu_scan（随化学势 μ）与 amset 的输运（随载流子浓度 n）
    #   各出一张图，分别放 transport_plots/{bt2_muscan, amset_scan}/，便于和文献 Fig.5/7/8 对点。
    def _scan_grid(ncol, titles):
        _f, _ax = plt.subplots(1, ncol, figsize=(3.15 * ncol, 2.9), squeeze=False)
        for j in range(ncol):
            _ax[0][j].set_title(titles[j], fontsize=9.5)
            _ax[0][j].grid(alpha=0.25, lw=0.4)
        return _f, _ax[0]

    # ⑥a bt2 μ 扫描（CRTA×DPT-τ，随 μ，x/y 两方向）
    try:
        _mscan = sorted((cwd / BT2_DIR).glob("paper_mu_scan_*.csv"))
        if _mscan:
            _mp = _mscan[-1]
            _bdir = plots / "bt2_muscan"; _bdir.mkdir(exist_ok=True)
            _br = list(csv.DictReader([l for l in _mp.read_text(errors="ignore").splitlines()
                                        if l.strip() and not l.startswith("#")]))
            _T = int(_mp.stem.rsplit("_", 1)[1].rstrip("K"))
            # [0045] 同时画两种口径：thickness（我们，S8.1 主口径）与 cell（文献 Fig.5/7 裸元胞口径）。
            #   σ/PF/κe 差 c/t 倍（factor 记在 CSV 头部注释）；S/ZT 对口径不变。
            _fac = None
            for _hl in _mp.read_text(errors="ignore").splitlines():
                if _hl.startswith("#") and "factor=" in _hl:
                    _seg = _hl.split("factor=")[1].split(",")[0].split()[0]
                    try:
                        _fac = float(_seg)
                    except ValueError:
                        _fac = None
                    break
            _muv = np.array([fnum(r["mu_rel_efermi_eV"]) for r in _br])
            # (x 列, y 列, 标题, 是否依赖口径)
            _bp = [("S_xx_uV/K", "S_yy_uV/K", "S (uV/K)", False),
                   ("sigma_xx_S/m", "sigma_yy_S/m", "sigma (S/m)", True),
                   ("PF_xx_W/mK2", "PF_yy_W/mK2", "PF (W/mK^2)", True),
                   ("kappa_e_WF_xx_W/mK", "kappa_e_WF_yy_W/mK", "kappa_e (W/mK)", True)]
            if "ZT_xx" in (_br[0] if _br else {}):
                _bp.append(("ZT_xx", "ZT_yy", "ZT", False))
            _f, _ax = _scan_grid(len(_bp), [p[2] for p in _bp])
            for j, (cx, cy, _t, _dep) in enumerate(_bp):
                a = _ax[j]
                for col, ls, lab in ((cx, "-", "x"), (cy, "--", "y")):
                    v = np.array([fnum(r[col]) for r in _br])
                    msk = np.isfinite(_muv) & np.isfinite(v)
                    a.plot(_muv[msk], v[msk], lw=1.15, ls=ls, label=lab + " (thickness)")
                    if _dep and _fac and abs(_fac - 1.0) > 1e-9:
                        a.plot(_muv[msk], v[msk] / _fac, lw=0.9, ls=":", label=lab + " (cell)")
                a.set_xlabel("mu (eV)", fontsize=8)
                if j == 0:
                    a.legend(fontsize=6.5, frameon=False, loc="upper right")
            _f.suptitle("%s  bt2 CRTA x DPT-tau  (chemical-potential scan, %d K)" % (M, _T), fontsize=11)
            _f.tight_layout(rect=[0, 0, 1, 0.93])
            _f.savefig(_bdir / ("bt2_muscan_%dK.png" % _T), dpi=130); plt.close(_f)
            (_bdir / ("bt2_muscan_%dK.csv" % _T)).write_text(_mp.read_text(errors="ignore"), encoding="utf-8")
            print("[plots] bt2_muscan/bt2_muscan_%dK.png" % _T)
    except Exception as _e:
        print("[WARN] bt2 μ 扫描图失败（%s: %s）" % (type(_e).__name__, _e))

    # ⑥b amset 浓度扫描（可用的分级机制结果，随 n，600 K）
    try:
        _adir = plots / "amset_scan"; _adir.mkdir(exist_ok=True)
        _tab = out / "intrinsic" / "intrinsic_transport_table.csv"
        if _tab.is_file():
            _r2 = list(csv.DictReader([l for l in _tab.read_text(errors="ignore").splitlines()
                                        if l.strip() and not l.startswith("#")]))
            _sub = [r for r in _r2 if abs(fnum(r["T_K"]) - 600.0) < 1e-6]
            if _sub:
                _aq = [("mu_cm2/Vs", "mu (cm^2/V/s)"), ("S_uV/K", "S (uV/K)"),
                       ("sigma_S/m", "sigma (S/m)"), ("PF_W/mK2", "PF (W/mK^2)"),
                       ("kappa_e_W/mK", "kappa_e (W/mK)")]
                _ads = [("overall", "tab:red", "-"), ("intrinsic", "tab:orange", "-."), ("ADP_only", "tab:green", ":"),
        ("ADP_POP", "tab:brown", "-"), ("ADP_POP_ODP", "tab:olive", "-")]
                _f2, _ax2 = _scan_grid(len(_aq), [q[1] for q in _aq])
                for j, (col, _t) in enumerate(_aq):
                    a = _ax2[j]
                    for ds, cc, ls in _ads:
                        dd = [r for r in _sub if r["dataset"] == ds]
                        if not dd:
                            continue
                        n = np.array([fnum(r["n_cell_cm-3"]) for r in dd])
                        y = np.array([fnum(r[col]) for r in dd])
                        keep = np.abs(n) >= NFLOOR
                        seg(a, n[keep], y[keep], color=cc, lw=1.3, ls=ls, marker="s", ms=3.0)
                    if col == "S_uV/K":
                        a.set_ylim(-1000, 1000); a.axhline(0, color="k", lw=0.5)
                    else:
                        a.set_yscale("log")
                    a.set_xlim(-4.15, 4.15); a.set_xticks([-4, -3, -2, -1, 0, 1, 2, 3, 4])
                    a.set_xticklabels(["1e21", "1e20", "1e19", "1e18", "1e17", "1e18", "1e19", "1e20", "1e21"], fontsize=6)
                    a.axvline(0, color="0.55", lw=0.6, ls=":")
                    a.set_xlabel("n (cm^-3)", fontsize=8)
                _present = {r["dataset"] for r in _sub}
                _handles = [Line2D([], [], color=cc, lw=1.3, ls=ls, marker="s", ms=3.0,
                                   label=_LBL.get(ds, ds))
                            for ds, cc, ls in _ads if ds in _present]
                if _handles:
                    _f2.legend(handles=_handles, loc="lower center", bbox_to_anchor=(0.5, 0.01),
                               ncol=min(5, len(_handles)), fontsize=7, frameon=False)
                _f2.suptitle("%s  AMSET concentration scan, 600 K" % M, fontsize=11)
                _f2.tight_layout(rect=[0, 0.08, 1, 0.93])
                _f2.savefig(_adir / "amset_scan_600K.png", dpi=130); plt.close(_f2)
                with open(_adir / "amset_scan_data_600K.csv", "w", newline="") as _fh:
                    _w = csv.writer(_fh); _w.writerow(list(_r2[0].keys()))
                    for r in _sub:
                        _w.writerow([r.get(c, "") for c in _r2[0].keys()])
                print("[plots] amset_scan/amset_scan_600K.png")
    except Exception as _e:
        print("[WARN] amset 浓度扫描图失败（%s: %s）" % (type(_e).__name__, _e))

    try:                                                          # [0040] 本征载流子（未掺杂）单独出图
        write_undoped_plots(cwd, plots, M, is_2d)
    except Exception as _e:                                       # noqa: BLE001
        print("[WARN] 本征载流子图失败（%s: %s）" % (type(_e).__name__, _e))
    try:                                                          # [0053][0054] 稀薄极限 μ(T) 机制分解（只这一个浓度）
        write_mechanism_resolved_plot(cwd, plots, M, is_2d)
    except Exception as _e:                                       # noqa: BLE001
        print("[WARN] mechanism-resolved 迁移率图失败（%s: %s）" % (type(_e).__name__, _e))
    try:                                                          # [0041] 缺陷平衡载流子单独出图
        write_defect_plots(cwd, plots, M, is_2d)
    except Exception as _e:                                       # noqa: BLE001
        print("[WARN] 缺陷平衡载流子图失败（%s: %s）" % (type(_e).__name__, _e))
    return True



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
        make_figure(out, am, bt, dpt, am2=am2, geom=(_geom_2d(cwd) if is_2d else None))   # [0032] 2D 双横轴
    except Exception as e:
        print("[WARN] 画图失败（%s: %s）——表已生成，图跳过" % (type(e).__name__, e))
    try:                                              # [0020] μ(T)，尽力而为，不影响完成标记
        write_mu_vs_T(out, cwd, is_2d)
    except Exception as e:                            # noqa: BLE001
        print("[WARN] μ(T) 表失败（%s: %s）——300 K 对比不受影响" % (type(e).__name__, e))
    try:                                              # [0032] 本征五参数图 + 表，尽力而为，不影响完成标记
        write_intrinsic_suite(out, cwd, is_2d)
    except Exception as e:                            # noqa: BLE001
        print("[WARN] 本征五参数图/表失败（%s: %s）——300 K 对比不受影响" % (type(e).__name__, e))
    try:                                              # [0033] 规定格式绘图输出（plots/ 子目录），尽力而为，不影响完成标记
        write_transport_plots(cwd, out, is_2d)
    except Exception as e:                            # noqa: BLE001
        print("[WARN] 规定格式绘图输出失败（%s: %s）——comparison_300K 不受影响" % (type(e).__name__, e))
    try:                                              # [0035] τ(E) / τ_eff vs DPT τ，尽力而为，不影响完成标记
        write_tau_suite(out, cwd, is_2d)
    except Exception as e:                            # noqa: BLE001
        print("[WARN] τ 对比失败（%s: %s）——300 K 对比不受影响" % (type(e).__name__, e))
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
            #   3D：黄色（Si 全网格对照已出：ADP 0.92~1.01、overall 0.92~1.06，IMP 逐机制 0.76~1.07 双向噪声）
            if overlap_grade_only_real(dim) == "3d":
                print("[重叠报警] **黄色**：三维 + 真实重叠（无 unity 对照）——"
                      "Si 全网格对照（V26/0.5.1 复跑）：ADP 0.92~1.01、overall 0.92~1.06（重掺档 ±6% 内），"
                      "IMP 逐机制 0.76~1.07 双向噪声、不可靠。**采纳前需人工确认逐机制值**。")
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
