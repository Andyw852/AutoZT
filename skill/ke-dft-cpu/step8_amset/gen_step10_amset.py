#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step8_amset.py —— 写 settings.yaml → amset run → σ/S/κ_e（step8_amset）。

汇集前面所有产物，写 amset 的 settings.yaml，提交到计算节点跑 amset run。
输入软链：
  wavefunction.h5   ← step4_wave
  deformation.h5    ← step7b_deform_read
介电常数从 step5_dielect/OUTCAR 解析，带隙从带隙段或配置读，弹性从 step6_elastic。
产出目录：step8_amset/，产物 transport.json（判据看 thermal_conductivity）。
"""
import glob
import os
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

# [SKILL_REV] 版本戳：写进 settings.yaml 头，便于从结果反查跑的是哪份 skill 副本。
_SKILL_REV = "2026-09-16-voigt-reorder"
import stepconf  # noqa: E402
try:
    import ke_common as kc
    _HAS_KC = True
except Exception:
    _HAS_KC = False

# =========================== 可改参数区 ===========================
# ---------- 体系判别阻断（step2.15_discriminant）----------
# False = 默认拦截 SEMIMETAL/METAL@PBE（这套半导体框架不适用）。
# True  = 强制继续（金属体系也有人要算输运，或 PBE 误判而杂化还没来得及重判）。
# ★ 不是硬停：读不到 discriminant.json 时一律放行，行为与加这道闸门前完全一致。
FORCE_TRANSPORT = False


def _disc_gate():
    """体系判别闸门：默认拦截 SEMIMETAL/METAL@PBE，FORCE_TRANSPORT 可覆盖。"""
    import sys as _s
    from pathlib import Path as _P
    try:
        _s.path.insert(0, str(_P(__file__).resolve().parent))
        import discriminant_common as _dc
        if not _dc.gate(_P.cwd(), "step8_amset", force=FORCE_TRANSPORT):
            _s.exit("[BLOCKED] 体系判别为 SEMIMETAL/METAL@PBE —— step8_amset 已阻断。"
                    "完整提示见 step2_bandgap/step2.15_discriminant/discriminant.json "
                    "的 block.hint；强制继续请把本脚本顶部 FORCE_TRANSPORT 设为 True。")
    except SystemExit:
        raise
    except Exception as _e:
        print("[WARN] 体系判别闸门异常，放行：%s" % _e, file=_s.stderr)


OUTDIR_NAME = "step8_amset"
WAVE_DIR    = "step4_wave"
# patch_wavefunction_full（2026-09-18，V26）：可选改读全网格 h5（step4b_wave_full，
#   ISYM=-1 写出全部 k 点）。三维"单变量"对照用：h5 与 vasprun 必须同源（都取
#   step3b/step4b），否则 h5 点数与 vasprun 推出的网格对不上，AMSET 会退回去对称化，
#   对照就不成立。默认 False = 仍读 step4_wave，行为与以前完全一致。
WAVEFUNCTION_FULL = False
# patch_unity_overlap_override（2026-09-18，V26）：step.conf 可显式覆盖 unity_overlap，
#   用于**受控对照**（同一份 settings 只翻这一个开关）：
#     auto  -> 现行行为（2D 恒 true；3D 不写 = 真实重叠）
#     true  -> 强制 unity（如给 3D / SOC 补一次 unity 比值，检查比值报警）
#     false -> 强制真实重叠（如给 LS 这种 2D 补"同 settings 的 real 对照"）
UNITY_OVERLAP = "auto"
READ_DIR    = "step7b_deform_read"
DIELECT_DIR = "step5_dielect"
# patch_ke_dag：跟上 v1.9 的目录重命名。HSE 优先，其次 PBE，最后兼容老目录。
BANDGAP_PLOT_CANDS = [
    "step2_bandgap/step2.3_hse_plot",
    "step2_bandgap/step2.2_pbe_plot",
    "step4_band_plot",
]
STEP_LABEL  = "S8_kappa"
# patch_overlap_preflight（V24/V25.10）：作业内先跑运行前检查 —— 运行目录里
#   wavefunction.h5 / vasprun.xml / ../step3_uniform 都在，②③④ 三项才真正生效；
#   不过就直接 exit 1，宁可不跑也不出不可信的数。脚本由 gen 复制进运行目录。
AMSET_CMD   = ('python overlap_preflight.py --in-job || exit 1; '
               'amset run >> amset.log 2>&1 && cp -f "$(ls -t transport_*.json 2>/dev/null | head -1)" transport.json && ls -l transport.json')
# --- 输运设置（可改）---
DOPING      = "-1e21:-1e17:5, 1e17:1e21:5"   # n 型 + p 型各 5 点（对数均布）cm^-3
TEMPERATURES = "100:900:9"            # 100,200,...,900 K，每 100 K 一个点
SCATTERING  = ["ADP", "IMP", "POP"]   # 形变势声学 + 电离杂质 + 极性光学 "ADP", "IMP", "POP"
MANUAL_BANDGAP = None                 # None=自动读；或写数值(eV) 覆盖 scissor
# patch_dielec_assert：介电产物必须通过物理性硬断言，否则直接停步（不静默往下传）。
#   检查项：eps_inf 存在 / 对角项 >= 1 / 不是单位矩阵（DFPT 初值签名）/
#   有 IONIC CONTRIBUTION 块（否则 Frohlich 耦合恒为 0、POP 静默丢失）。
#   已知要强行继续时改 False —— 但结果里必须注明 POP 不生效。
REQUIRE_PHYSICAL_DIELECTRIC = True
# patch_require_bandgap：半导体输运没有 scissor 的结果没有意义——
#   读不到带隙直接停步，不再静默用 PBE 带隙跑完。金属体系才设 False。
REQUIRE_BANDGAP = True
# patch_interp_factor：AMSET 的收敛判据在**插值后**的网格上，别吃默认值 5。
#   设 None 则不写这一行（回到 AMSET 默认）。加大前先做收敛测试。
INTERPOLATION_FACTOR = 10
# patch_unity_overlap_3d（2026-09-17，VERIFICATION V25）：
#   显式把 unity_overlap 写进 settings.yaml —— 以前这一行是缺的，AMSET 默认 False
#   （真实重叠），于是所有走本步的项目都默认用"真实重叠 + 去对称化 h5"这个组合。
#   ★ 分维处理（用户 2026-09-17 定）：
#     - 2D：**直接写 unity_overlap: true**。二维已实测去对称化会把重叠算坏
#       （MoS2：ADP 迁移率被抬高 10~16 倍，V23），补丁已在 8.4 落地。
#     - 3D：**不强制**，只打 WARNING。三维受影响面已确认（Si / 225 / Al5N 都走了这条路），
#       但"是否真算错"要等 Si 的全网格对照（V25）。Si 有 6 个等价能谷，
#       若改用 unity 最多会把迁移率低估约 6 倍，代价比二维大得多 ——
#       所以结论出来之前不擅自改。
UNITY_OVERLAP_2D = True
# NWORKERS：AMSET 的并行进程数，会写进 settings.yaml 的 nworkers。
#   ★ 2026-09-16 新增：此前**不写 nworkers**，AMSET 取默认 -1，其源码
#     (amset/interpolation/bandstructure.py:182) 把 -1 解释成
#         nworkers = multiprocessing.cpu_count() if nworkers == -1 else nworkers
#     即**用满节点全部核**，而提交模板只申请 24 核 —— 于是超订节点 8 倍
#     （jzzn 节点 192 核），抢占同节点其他作业的 CPU、违反 SLURM 分配，
#     且耗时随节点负载剧烈波动。
#   实测（Pb2Sb2Te5，其余设置相同）：nworkers=-1 时 4 分钟；nworkers=24 时 7.5 小时。
#   ★ 默认 None = **自动**：由提交模板（叠加 step.conf 的 [submit] 覆盖）推出本次
#     作业真正分配到的核数 ntasks_per_node × cpus_per_task，拿它当 nworkers。
#     为什么要自动：三台机器的分配口径本来就不同 ——
#         jzzn      24 × 1 = 24
#         hanhai25 128 × 1 = 128     （写死 24 = 白占 104 核）
#         3090       1 × 8 = 8       （写死 24 = 超订 3 倍）
#     写死一个常数在三台上必然错两台：少要 = 算力闲置，多要 = 抢别人的 CPU。
#   ★ 要固定值就在 step.conf 写 NWORKERS = <int>（技能侧
#     skill/ke-dft-cpu/step8_amset/step.conf，或某材料的
#     project_setting/templates/step8_amset/step.conf）；与本次分配不一致时打印 WARN。
#     AMSET 在掺杂/温度维度上接近线性加速，想更快就把提交模板的核数提上去 ——
#     自动模式会跟着走，不用再手动同步两个数。
NWORKERS = None
NWORKERS_FALLBACK = 24      # 提交模板里读不到 SLURM 分配（--ntasks/--cpus）时的兜底
# --- 弹性常数来源（amset run 的 ACD 散射需要）---
#   MANUAL_ELASTIC 填了就用它，否则从 ELASTIC_DIR/OUTCAR 自动解析（kBar→GPa）。
#   直接填：单个数（各向同性近似，GPa），或 6x6 列表（完整 Cij，GPa）。
MANUAL_ELASTIC = None
ELASTIC_DIR = "step6_elastic"
# patch_deform_ref（2026-09-16）：形变势参考口径。step7b 若已写出
# deformation_vac.h5（真空静电势零点，二维文献通行做法），2D 默认用它；
# 否则退回 deformation.h5（AMSET 芯态对齐，官方为体材料设计）。
# ★ 实测影响很大：CrS2 单层 E1_core=3.53 eV vs E1_vac=5.86 eV -> μ 差 2.76 倍。
DEFORM_REF = "auto"      # auto = 2D 优先 vacuum；或强制 "vacuum" / "core"
# --- patch_2d_amset：2D 修正 ---------------------------------------------
# TWO_D_MODE: "auto" = 读 step1 的 workflow_method.txt 的 DIM=；"on"/"off" 强制。
TWO_D_MODE = "auto"
# LAYER_THICKNESS: 层的有效厚度（Å）。三种取值：
#   "vdw"  （默认）自动 = 原子核 z 跨度 + 两侧最外层元素各一个范德华半径。
#          这是 2D 文献做 h/d0 重标度时的通行取法，确定、可复现、可引用。
#   "span" 只用原子核 z 跨度。偏小 -> C 偏大 -> 迁移率偏高，一般别用。
#   数值    手填（Å）。已知对应体相层间距时优先用这个。
# 无论哪种，最终用的 c/t/取法/倍数都会落盘到 step8_amset/2d_correction.json。
LAYER_THICKNESS = "vdw"
# patch_layer_thickness_conf：LAYER_THICKNESS 可从 step.conf 覆盖（出厂默认 "vdw"）。
#   SS/LS 超晶格在 project_setting/templates/step8_amset/step.conf 写 LAYER_THICKNESS = 6.73
#   （照文献取 CrS2 6.53 / CrSe2 6.94 的平均，非本结构 vdW 自动值 6.97）。
STEP = "step8_amset"
SPEC = {
    "LAYER_THICKNESS": (LAYER_THICKNESS, "str"),
    # NWORKERS 必须在 SPEC 里声明：step.conf 的严格模式会把没声明的键判成
    # "本脚本不认识的键" 而 SystemExit（下面对 SystemExit 的兜底会把它吞掉 →
    # 静默回默认值，用户在 step.conf 里写了等于没写）。
    # 技能侧的落点是 skill/ke-dft-cpu/step8_amset/step.conf，不是技能全局 step.conf。
    "NWORKERS": (NWORKERS, "int"),
    # 全网格 h5 对照（三维裁定）；默认 False，既有三维项目行为不变。
    "WAVEFUNCTION_FULL": (WAVEFUNCTION_FULL, "bool"),
    # unity_overlap 覆盖（受控对照用）；默认 auto = 现行行为不变。
    "UNITY_OVERLAP": (UNITY_OVERLAP, "str"),
}
# 2D 时给 settings.yaml 写 free_carrier_screening: true
FREE_CARRIER_SCREENING_2D = True
# patch_2d_dielec：2D 介电取面内 (εx+εy)/2 并扣真空 ε^m=1+(L/t)(ε_sup-1)
# （文献标准做法；L/t 复用弹性同款 c/t）。设 False 则用原始三对角平均。
TWO_D_DIELECTRIC_VACUUM = True
# 找结构的候选目录（按序取第一个有 POSCAR/CONTCAR 的）
STRUCT_CANDS = ["step3_uniform", "step7_deform", "step1_opt", "step1_std_opt"]
_STEP1_METHOD_CANDS = ["step1_opt", "step1_std_opt",
                       "step1c_PBE_opt", "step1b_PBE_opt", "step1a_PBE_opt"]
# =================================================================


def read_dielectric(dielect_dir: Path):
    """从 DFPT OUTCAR 读完整 3x3 ε∞ 与 ε₀ 张量。"""
    oc = dielect_dir / "OUTCAR"
    if not oc.is_file():
        return None, None
    txt = oc.read_text(errors="ignore")
    def grab(tag):
        # 取 tag 之后第一块 3x3，对角平均
        i = txt.rfind(tag)
        if i < 0:
            return None
        rows = []
        for ln in txt[i:].splitlines()[1:]:
            nums = re.findall(r"-?\d+\.\d+", ln)
            if len(nums) >= 3:
                rows.append([float(x) for x in nums[:3]])
            if len(rows) == 3:
                break
        if len(rows) < 3:
            return None
        return rows
    eps_inf = grab("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)")
    eps_0 = grab("MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION")

    # ================= patch_dielec_assert：硬断言，不许静默错 =================
    # 这条链上已经出现两次"静默错"伤人的案例：
    #   1) OUTCAR 的 fundamental gap 行在带重叠时给正值 -> 半金属被报成有隙；
    #   2) 下面的 eps_static = eps_inf 兜底 -> Frohlich 耦合 (1/eps_inf - 1/eps_static)
    #      **恒等于 0** -> POP 散射被整个丢掉，只在日志里留一行 WARN。
    #      实测签名：max|mobility/POP| = 2.3e48（分母趋零），而 sigma 被系统性高估。
    # 所以这里**抛错而不往下传**。确属已知情形要强行继续，显式改本开关。
    # ★ NaN 必须**单独**查（2026-09-16 从 taskflow-v2.0 同步）：Python 里
    #   float('nan') < 0 / < 1.0 全是 False，所以下面所有大小比较对 NaN 一律为假 ——
    #   全 NaN 的张量会完美绕过所有守卫。实测 2026-09-14：Pb2Sb2Te5 / Sn2Sb2Te5 的
    #   step5_dielect 张量全是 NaN，正是从这条缝里静默流进 settings.yaml 的。
    def _has_nan(M):
        if M is None:
            return False
        import math as _m
        return any(_m.isnan(float(x)) for r in M for x in r)

    def _diele_fail(msg):
        if REQUIRE_PHYSICAL_DIELECTRIC:
            raise ValueError(
                "[step5_dielect 产物不可用] %s\n"
                "  覆盖方式：本脚本顶部 REQUIRE_PHYSICAL_DIELECTRIC = False\n"
                "  最常见根因：step5_dielect 的 DFPT 崩过/没跑完。先查那里的 OUTCAR 有没有\n"
                "  完整写出 'MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field\n"
                "  effects in DFT)' 与 '... IONIC CONTRIBUTION' 两块，以及 queue.err 里\n"
                "  有没有 signal 11。DFPT 的 INCAR 走 skill 模板（KPAR=1/NCORE=1/LPEAD=.FALSE.）。"
                % msg)
        print("[WARN][已按开关放行] %s" % msg)

    if eps_inf is None:
        _diele_fail("没读到 eps_inf（%s/OUTCAR 里没有 "
                    "'MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects "
                    "in DFT)' 块）—— DFPT 没跑完或崩了。" % dielect_dir)
        return None, None
    if _has_nan(eps_inf):
        _diele_fail("eps_inf 含 NaN（%s/OUTCAR）—— DFPT 发散/未收敛。"
                    "本链实测失效体系：Pb2Sb2Te5、Sn2Sb2Te5（重元素窄隙）；"
                    "根因是 LPEAD=.TRUE.，DFPT 的 INCAR 走 skill 模板（LPEAD=.FALSE.）。"
                    % dielect_dir)
    if _has_nan(eps_0):
        _diele_fail("IONIC 块含 NaN（%s/OUTCAR）—— 离子介电段发散。"
                    "注意：DFPT 只要离子段出 NaN，整张张量都不可用。" % dielect_dir)
    # 电子静态响应必须 >= 1；< 1 在正常绝缘体里不可能出现
    _diag_inf = [eps_inf[i][i] for i in range(3)]
    if min(_diag_inf) < 1.0:
        _diele_fail("eps_inf 对角项 %.4f %.4f %.4f 中有 < 1 的分量 —— "
                    "电子静态介电不可能小于 1，这是垃圾值。" % tuple(_diag_inf))
    # 单位矩阵 = DFPT 的**迭代初值**签名（零响应）。跑挂的作业会把初值留在 OUTCAR 里。
    if max(abs(eps_inf[i][j] - (1.0 if i == j else 0.0))
           for i in range(3) for j in range(3)) < 1e-6:
        _diele_fail("eps_inf 恰为单位矩阵 —— 这是 DFPT 迭代初值的签名（零响应），"
                    "不是结果。DFPT 在写出 MACROSCOPIC 块之前就结束了。")

    # ε₀ = 电子 + 离子
    if eps_0 is None:
        _diele_fail(
            "有 eps_inf 但**没有** IONIC CONTRIBUTION 块（离子介电没算出来）。"
            "若按旧逻辑兜底成 eps_static = eps_inf，则 Frohlich 耦合 "
            "(1/eps_inf - 1/eps_static) 恒为 0，**POP 散射会被静默丢弃**，"
            "sigma 与 mobility 被系统性高估。必须先把离子贡献算出来。")
        # 显式放行时仍退回 eps_inf，但把"POP 已失效"写进返回值之外可查的地方
        print("[WARN] POP 散射本次不生效（eps_static == eps_inf），"
              "结果里应注明该项缺失。")
        return eps_inf, eps_inf
    eps_static = [[eps_inf[i][j] + eps_0[i][j] for j in range(3)] for i in range(3)]
    # patch_dielec_guard：ε_ionic<0 / ε_static<=0 非物理（多为 Γ 近零声学模
    # 污染 DFPT 离子介电）——弃离子项、退回 ε∞ 兜底并告警
    if _has_nan(eps_0):
        _diele_fail("离子介电对角含 NaN —— 下面的 '< 0' 判断对 NaN 恒为假，会漏过；"
                    "这里显式拦截。")
    if any(eps_0[i][i] < 0 for i in range(3)):
        raise ValueError(
            "离子介电对角项为负（eps_ionic=%.2f:%.2f:%.2f，eps_static=%.2f），非物理，"
            "多因 Gamma 近零声学模污染 DFPT 离子介电；必须检查 DFPT 收敛与声子稳定性。"
            % (eps_0[0][0], eps_0[1][1], eps_0[2][2], eps_static[0][0]))
    _coup = [1.0 / eps_inf[i][i] - 1.0 / eps_static[i][i] for i in range(3)]
    if max(abs(c) for c in _coup) < 1e-9:
        if _is_nonpolar():
            # 单元素非极性体系（Si / Ge / 金刚石…）没有 IR 活性声子，离子介电
            # 恒为 0 -> eps_static == eps_inf 是**物理正确**结果（不是 DFPT 失效）。
            # POP 散射会在下游 _amset_scatterers() 里按物理口径被剔除，
            # 所以这里不能当成"step5_dielect 产物不可用"拦下。
            print("[WARN] 非极性体系 eps_static == eps_inf（物理正确）——"
                  "POP 散射将在散射列表里被剔除，结果里应注明该项缺失。")
        else:
            _diele_fail("eps_static 与 eps_inf 对角项相同，Frohlich 耦合 %.3e 恒为 0 -> "
                        "POP 散射不生效。" % max(abs(c) for c in _coup))
    return eps_inf, eps_static


# --------------------------------------------------------------------------
# patch_voigt_order：VASP → 标准 Voigt 顺序重排
# --------------------------------------------------------------------------
# 【已确认 bug】VASP OUTCAR 'TOTAL ELASTIC MODULI (kBar)' 的行/列标签顺序是
#   XX YY ZZ XY YZ ZX，
# 而 AMSET 走 pymatgen 的 Tensor.from_voigt()，按**标准 Voigt** 顺序解析：
#   XX YY ZZ YZ XZ XY（即 1,2,3,4,5,6 = xx,yy,zz,yz,xz,xy）。
# 直接把 VASP 的行序当标准 Voigt 交给 AMSET，等于把第 4 位(XY)与第 6 位(ZX)互换：
#   · 2D slab：面内剪切 C_xyxy 被换成接近 0 的 C_zxzx（实测 CrS2 单层：
#     C66 = 22.29 GPa 被塞进 YZ 槽位，而 YZ 槽位 −0.10 GPa 被当成面内剪切），
#     面内 TA 支刚度近零 → ADP 散射率与迁移率全错。
#   · 3D 低对称（225 相 tetradymite 这类 C44≠C66、C14≠0 的晶系）同样读错；
#     立方晶系三个剪切分量相等，不受影响。
# 行列必须**同时**重排（张量指标重命名），所以用 C'[i][j] = C[p[i]][p[j]]。
#   实测依据：step6_elastic/OUTCAR 的 XY 行对角线 = 222.8673 kBar = 22.29 GPa。
_VASP2VOIGT = (0, 1, 2, 4, 5, 3)      # new[i] = old[_VASP2VOIGT[i]]，纯 python，不引 numpy
_VOIGT_LABELS = ("XX", "YY", "ZZ", "YZ", "XZ", "XY")        # 标准 Voigt（AMSET/pymatgen）
_VASP_DEFAULT_LABELS = ("XX", "YY", "ZZ", "XY", "YZ", "ZX")  # VASP OUTCAR 缺省


def _norm_voigt_label(s):
    """XZ 与 ZX 是同一个剪切分量，统一成 XZ。"""
    s = str(s).upper()
    return "XZ" if s == "ZX" else s


def _voigt_perm(src_labels):
    """src_labels[j] = 来源第 j 位的标签；返回 perm 使 out[i] = src[perm[i]]（标准 Voigt）。
    标签不全时返回 None（调用方走缺省顺序兜底）。"""
    try:
        src = [_norm_voigt_label(s) for s in src_labels]
        return tuple(src.index(_norm_voigt_label(l)) for l in _VOIGT_LABELS)
    except ValueError:
        return None


def _reorder_voigt(mat, col_labels, row_labels=None):
    """按来源列/行标签把 6×6 张量重排成标准 Voigt 顺序。"""
    pc = _voigt_perm(col_labels)
    pr = _voigt_perm(row_labels if row_labels is not None else col_labels)
    return [[mat[pr[i]][pc[j]] for j in range(6)] for i in range(6)]


def vasp_to_voigt_6x6(m):
    """6×6 弹性/刚度张量：VASP 顺序 (XX YY ZZ XY YZ ZX) → 标准 Voigt (XX YY ZZ YZ XZ XY)。

    只在 OUTCAR 表头认不出来时的兜底路径用；正常走 _reorder_voigt（按表头实测顺序）。"""
    a = [[float(x) for x in row] for row in m]
    return [[a[i][j] for j in _VASP2VOIGT] for i in _VASP2VOIGT]


def vasp_to_voigt_piezo(e3x6):
    """3×6 压电张量（VASP 列序 XX YY ZZ XY YZ ZX）→ 标准 Voigt 列序（只重排列）。"""
    a = [[float(x) for x in row] for row in e3x6]
    return [[row[j] for j in _VASP2VOIGT] for row in a]


def _pick_deformation_h5(cwd, read_dir):
    """挑形变势文件：2D 且存在 deformation_vac.h5 时用真空口径，否则 core 口径。

    DEFORM_REF: "auto"（2D 优先 vacuum）/ "vacuum" / "core"。
    缺失时只告警不回退维度闸门 —— 真空口径是**推荐**而非硬要求，
    但要在日志里说清楚本次用的是哪个口径（绝对值差 ~2.8 倍）。
    """
    read_dir = Path(read_dir)
    ref = str(DEFORM_REF).lower()
    vacuum_ok = (read_dir / "deformation_vac.h5").is_file()
    mode = read_dim(cwd)
    if ref == "auto":
        want = (mode == "2d")
    elif ref in ("vacuum", "core"):
        want = (ref == "vacuum")
    else:
        sys.exit("[ERROR] DEFORM_REF=%r 无效，只允许 auto / vacuum / core" % DEFORM_REF)
    if want and vacuum_ok:
        print("[OK] 形变势参考口径：**真空**（deformation_vac.h5，二维文献口径）")
        return "deformation_vac.h5"
    if want and not vacuum_ok:
        print("[WARN] 想要真空口径但 %s/deformation_vac.h5 不存在 -> 退回 **core** 口径"
              "（AMSET 芯态对齐）。要真空口径请先跑 step7b（需要形变目录里的 LOCPOT）。"
              % read_dir)
    else:
        print("[..] 形变势参考口径：core（AMSET 芯态对齐）")
    return "deformation.h5"


def read_elastic(cwd: Path):
    """弹性常数来源：MANUAL_ELASTIC 优先，否则从 step6_elastic/OUTCAR 解析。
    返回 amset settings.yaml 用的值：单标量(GPa) 或 6x6 列表(GPa，**标准 Voigt 顺序**)，
    读不到返回 None。

    MANUAL_ELASTIC 的约定同样是**标准 Voigt 顺序**（XX YY ZZ YZ XZ XY）。
    注意 ELASTIC_DIR/OUTCAR 对应的是 step6_elastic/OUTCAR（最后一块
    TOTAL ELASTIC MODULI，VASP 顺序）。"""
    if MANUAL_ELASTIC is not None:
        return MANUAL_ELASTIC
    oc = cwd / ELASTIC_DIR / "OUTCAR"
    if not oc.is_file():
        return None
    txt = oc.read_text(errors="ignore")
    i = txt.rfind("TOTAL ELASTIC MODULI")
    if i < 0:
        return None
    lines = txt[i:].splitlines()
    # 表头第 2 行给列顺序（" Direction  XX  YY  ZZ  XY  YZ  ZX"），每行行首给行标签。
    # 按表头**实测**顺序建置换，比写死"VASP 一定是 XX YY ZZ XY YZ ZX"稳。
    col_labels = None
    for ln in lines[:4]:
        p = ln.split()
        if len(p) >= 7 and p[0].upper().startswith("DIRECTION"):
            col_labels = [x.upper() for x in p[1:7]]
            break
    rows, row_labels = [], []
    labels = ("XX", "YY", "ZZ", "XY", "YZ", "ZX")
    for ln in lines:
        p = ln.split()
        if p and p[0].upper() in labels and len(p) >= 7:
            try:
                rows.append([float(x) for x in p[1:7]])
                row_labels.append(p[0].upper())
            except ValueError:
                pass
        if len(rows) == 6:
            break
    if len(rows) != 6:
        return None
    # VASP 输出 kBar，amset 要 GPa
    mat = [[v / 10.0 for v in r] for r in rows]
    # patch_voigt_order：来源顺序 -> 标准 Voigt 顺序（行列同时重排，见文件上部说明）
    src = col_labels or row_labels
    if _voigt_perm(src) is None or _voigt_perm(row_labels) is None:
        print("[WARN] 弹性张量的行列标签认不全（表头 %s / 行标签 %s）——"
              "按 VASP 缺省顺序 XX YY ZZ XY YZ ZX 重排" % (col_labels, row_labels))
        mat = vasp_to_voigt_6x6(mat)
        src = list(_VASP_DEFAULT_LABELS)
    else:
        mat = _reorder_voigt(mat, src, row_labels)
    print("[OK] 弹性常数：来源顺序 %s 已重排为标准 Voigt (XX YY ZZ YZ XZ XY)"
          "——不重排会把面内剪切 C66(XY) 与 C_zxzx 互换" % (" ".join(src),))
    return [[round(v, 3) for v in r] for r in mat]


def read_bandgap(cwd: Path):
    if MANUAL_BANDGAP is not None:
        return float(MANUAL_BANDGAP)
    import json
    bs = None
    for _c in BANDGAP_PLOT_CANDS:
        _p = cwd / _c / "band_summary.json"
        if _p.is_file():
            bs = _p
            print("[..] 带隙来源：%s" % _c)
            break
    if bs is not None:
        try:
            d = json.loads(bs.read_text())
            # patch_bandgap_key：band_summary.json 写的键是 gap_eV，原列表全对不上
            for k in ("gap_eV", "band_gap", "bandgap", "gap", "Egap"):
                if k in d:
                    _g = float(d[k])
                    if _g <= 0:
                        print("[WARN] band_summary 判为金属（gap=%.4f eV）——不写 bandgap"
                              % _g)
                        return None
                    print("[OK] 带隙 %.4f eV（键 %s）" % (_g, k))
                    return _g
        except Exception as _e:
            # ★ 不静默：这里失败会**静默退回项目配置里的 bandgap**，即用了另一个
            #   来源的带隙（可能差几十甚至几百 meV）而不自知。必须喊出来。
            print("[WARN] 解析 %s 失败（%s: %s）—— 将回退到项目 setting.yaml 的"
                  " bandgap。两者可能不是同一个数，请核对。"
                  % (bs, type(_e).__name__, _e), file=sys.stderr)
    # 退回项目配置
    for cand in (cwd / "project_setting" / "setting.yaml",):
        if cand.is_file():
            for ln in cand.read_text(errors="ignore").splitlines():
                m = re.match(r"\s*bandgap\s*:\s*([\d.]+)", ln)
                if m:
                    return float(m.group(1))
    return None


def link(out: Path, src: Path, name: str):
    if not src.is_file():
        sys.exit("[ERROR] 缺 %s：%s（前置步骤没完成？）" % (name, src))
    dst = out / name
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(os.path.relpath(src, out))
    print("[OK] 软链 %s" % name)


# --------------------------------------------------------------------------
# patch_2d_amset：2D 几何与弹性常数重标度
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# patch_2d_thickness：范德华半径表与 vdW 层厚
# --------------------------------------------------------------------------
# Bondi(1964) + Alvarez(2013) 常用值，单位 Å。故意内置而不用 pymatgen：
# gen 脚本在登录节点用系统 python 跑，ke_common 的设计就是不引入 pymatgen。
# patch_vdw_shared：单一真源在 skill/_common/vdw_radii.py，kl-dft-cpu 读同一份。
#   登录节点若拿不到公共池，回落最小内置表（值与公共池一致，只保证不崩）。
try:
    from vdw_radii import VDW_FALLBACK, VDW_RADII        # noqa: F401
except ImportError:
    VDW_RADII = {
        "H": 1.20, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55,
        "O": 1.52, "F": 1.47, "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10,
        "P": 1.80, "S": 1.80, "Cl": 1.75, "K": 2.75, "Ca": 2.31, "Ti": 2.11,
        "V": 2.07, "Cr": 2.06, "Mn": 2.05, "Fe": 2.04, "Co": 2.00, "Ni": 1.97,
        "Cu": 1.96, "Zn": 2.01, "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90,
        "Br": 1.85, "Mo": 2.17, "Ru": 2.13, "Rh": 2.10, "Pd": 2.10, "Ag": 2.11,
        "Cd": 2.18, "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98,
        "W": 2.18, "Pt": 2.13, "Au": 2.14, "Hg": 2.23, "Tl": 1.96, "Pb": 2.02,
        "Bi": 2.07,
    }
    VDW_FALLBACK = 2.00
    print("[WARN] 没找到 _common/vdw_radii.py，用内置最小表")


def _read_poscar_species_z(path: Path):
    """返回 [(元素符号, z 坐标 Å), ...] 和 c 轴长度。解析失败返回 (None, None)。"""
    try:
        lines = path.read_text(encoding="utf-8-sig",
                               errors="ignore").splitlines()
        scale = float(lines[1].split()[0])
        vecs = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
        _axb = (vecs[0][1]*vecs[1][2]-vecs[0][2]*vecs[1][1],  # patch_2d_height
                vecs[0][2]*vecs[1][0]-vecs[0][0]*vecs[1][2],
                vecs[0][0]*vecs[1][1]-vecs[0][1]*vecs[1][0])
        _A = (_axb[0]**2 + _axb[1]**2 + _axb[2]**2) ** 0.5
        c_len = (abs(vecs[2][0]*_axb[0]+vecs[2][1]*_axb[1]+vecs[2][2]*_axb[2]) / _A
                 if _A > 1e-9 else
                 (vecs[2][0]**2 + vecs[2][1]**2 + vecs[2][2]**2) ** 0.5)
        row6 = lines[5].split()
        if row6 and row6[0].lstrip("-").isdigit():
            return None, c_len              # VASP4 无元素符号行，认不出元素
        symbols = row6
        counts = [int(x) for x in lines[6].split()]
        start = 8
        mode = lines[7].strip().lower()
        if mode[:1] == "s":                 # Selective dynamics
            mode = lines[8].strip().lower()
            start = 9
        direct = mode[:1] in ("d", "")
        species = []
        for sym, n in zip(symbols, counts):
            species += [sym] * n
        out = []
        for sym, ln in zip(species, lines[start:start + len(species)]):
            p = ln.split()
            if len(p) < 3:
                continue
            z = float(p[2])
            out.append((sym, z * c_len if direct else z))
        return (out or None), c_len
    except Exception as e:                  # noqa: BLE001
        print("[WARN] 解析 %s 的元素/坐标失败：%s" % (path, e))
        return None, None


def vdw_thickness(cwd: Path):
    """t_vdW = (z_max + r_vdW[最上层元素]) - (z_min - r_vdW[最下层元素])。

    返回 (t, 说明字典)；认不出元素时返回 (None, 说明)。
    """
    for d in STRUCT_CANDS:
        for fn in ("CONTCAR", "POSCAR"):
            p = cwd / d / fn
            if not p.is_file():
                continue
            az, c_len = _read_poscar_species_z(p)
            if not az:
                continue
            top = max(az, key=lambda x: x[1])
            bot = min(az, key=lambda x: x[1])
            miss = [s for s in (top[0], bot[0]) if s not in VDW_RADII]
            if miss:
                print("[WARN] 范德华半径表里没有 %s，用 %.2f Å 兜底"
                      % ("/".join(miss), VDW_FALLBACK))
            r_top = VDW_RADII.get(top[0], VDW_FALLBACK)
            r_bot = VDW_RADII.get(bot[0], VDW_FALLBACK)
            t = (top[1] + r_top) - (bot[1] - r_bot)
            info = {
                "source_file": "%s/%s" % (d, fn),
                "z_span_A": round(top[1] - bot[1], 4),
                "top_species": top[0], "top_vdw_radius_A": r_top,
                "bottom_species": bot[0], "bottom_vdw_radius_A": r_bot,
                "vdw_radii_source": "Bondi 1964 / Alvarez 2013",
                "thickness_A": round(t, 4),
            }
            return round(t, 4), info
    return None, {"error": "找不到可解析元素符号的 POSCAR/CONTCAR"}


def write_2d_record(out: Path, record):
    import json
    (out / "2d_correction.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n")
    print("[OK] 2d_correction.json：修正依据已落盘（写论文直接引用）")


def read_dim(cwd: Path):
    """从 step1 的 workflow_method.txt 读 DIM=（'2d'/'3d'），读不到返回 None。"""
    for name in _STEP1_METHOD_CANDS:
        mf = cwd / name / "workflow_method.txt"
        if not mf.is_file():
            continue
        for ln in mf.read_text(errors="ignore").splitlines():
            if ln.strip().upper().startswith("DIM="):
                v = ln.split("=", 1)[1].strip().lower()
                return v if v in ("2d", "3d", "0d") else None
        return None
    return None


def _read_poscar_cz(path: Path):
    """从 POSCAR 读 (h⊥ Å, 原子 z 向跨度 Å)。解析不了返回 (None, None)。

    **h⊥ = V/|a×b| = |c·(a×b)|/|a×b|**，即垂直于 ab 面的胞高，不是 |c|。
    AMSET/phono3py 用的体积是 V = A·h⊥，所以 2D 归一化因子必须是 h⊥/t，
    二维映射里 kz 积分长度 |b₃| = 2π/h⊥ 也应该用 h⊥ —— 正交胞时 h⊥ = |c|，
    c 相对层法向倾斜时两者才不同（本流程强制真空沿第 3 个矢量、层法向 ∥ z，
    dim_common 与插件都会拦）。
    """
    try:
        lines = path.read_text(encoding="utf-8-sig", errors="ignore").splitlines()
        scale = float(lines[1].split()[0])
        vecs = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
        _axb = (vecs[0][1]*vecs[1][2]-vecs[0][2]*vecs[1][1],  # patch_2d_height
                vecs[0][2]*vecs[1][0]-vecs[0][0]*vecs[1][2],
                vecs[0][0]*vecs[1][1]-vecs[0][1]*vecs[1][0])
        _A = (_axb[0]**2 + _axb[1]**2 + _axb[2]**2) ** 0.5
        c_len = (abs(vecs[2][0]*_axb[0]+vecs[2][1]*_axb[1]+vecs[2][2]*_axb[2]) / _A
                 if _A > 1e-9 else
                 (vecs[2][0]**2 + vecs[2][1]**2 + vecs[2][2]**2) ** 0.5)
        row6 = lines[5].split()
        if row6 and row6[0].lstrip("-").isdigit():
            counts = [int(x) for x in row6]
            start = 7
        else:
            counts = [int(x) for x in lines[6].split()]
            start = 8
        n = sum(counts)
        mode = lines[start - 1].strip().lower()
        while mode[:1] in ("s",):                      # Selective dynamics
            mode = lines[start].strip().lower()
            start += 1
        direct = mode[:1] in ("d", "")
        zs = []
        for ln in lines[start:start + n]:
            p = ln.split()
            if len(p) < 3:
                continue
            z = float(p[2])
            zs.append(z * c_len if direct else z)
        if not zs:
            return c_len, None
        return c_len, max(zs) - min(zs)
    except Exception as e:                              # noqa: BLE001
        print("[WARN] 解析 %s 失败（%s），无法自动定 2D 几何" % (path, e))
        return None, None


def get_slab_geometry(cwd: Path):
    """返回 (c 轴长度 Å, 原子 z 跨度 Å)；都找不到返回 (None, None)。"""
    for d in STRUCT_CANDS:
        for fn in ("CONTCAR", "POSCAR"):
            p = cwd / d / fn
            if p.is_file():
                c_len, span = _read_poscar_cz(p)
                if c_len:
                    print("[..] 2D 几何来源：%s/%s" % (d, fn))
                    return c_len, span
    return None, None


def rescale_elastic_2d(elastic, factor):
    """把弹性常数整体乘 factor（= c/t）。标量和 6x6 都支持。"""
    if elastic is None:
        return None
    if isinstance(elastic, (int, float)):
        return round(elastic * factor, 4)
    return [[round(v * factor, 4) for v in row] for row in elastic]


def _doping_endpoints(spec):
    """从 DOPING 串里抠出出现过的浓度数值（cm^-3），用于打印面浓度换算。"""
    vals = []
    for tok in re.split(r"[,\s]+", str(spec)):
        if not tok:
            continue
        for part in tok.split(":"):
            try:
                v = float(part)
            except ValueError:
                continue
            if abs(v) >= 1e10:          # 过滤掉 ":5" 这种点数
                vals.append(v)
    return vals


def apply_2d_corrections(cwd: Path, elastic):
    """返回 (is_2d, elastic_修正后, c_len)。3D 或判定不出时原样返回。"""
    mode = str(TWO_D_MODE).lower()
    if mode == "off":
        return False, elastic, None
    if mode == "auto":
        dim = read_dim(cwd)
        if dim != "2d":
            return False, elastic, None
    elif mode != "on":
        sys.exit("[ERROR] TWO_D_MODE=%r 无效，只允许 auto / on / off" % TWO_D_MODE)

    c_len, span = get_slab_geometry(cwd)
    if c_len is None:
        print("[WARN] 2D 体系但读不到结构，弹性常数未做 c/t 修正——"
              "结果的绝对值不可用。请手填 MANUAL_ELASTIC。")
        return True, elastic, None

    # patch_2d_thickness：层厚三种取法，结果连同依据一起落盘
    mode_t = LAYER_THICKNESS
    t_info = {}
    if isinstance(mode_t, str) and mode_t.lower() == "vdw":
        t, t_info = vdw_thickness(cwd)
        t_info["method"] = "vdw"
        if t is None:
            print("[WARN] 认不出元素（VASP4 格式 POSCAR？），vdW 层厚失败，"
                  "回退到原子 z 跨度 %s Å" % span)
            t, t_info = span, {"method": "span_fallback",
                               "z_span_A": span}
        else:
            print("[OK] vdW 层厚：z 跨度 %.3f Å + %s(%.2f) + %s(%.2f) = %.3f Å"
                  % (t_info["z_span_A"], t_info["top_species"],
                     t_info["top_vdw_radius_A"], t_info["bottom_species"],
                     t_info["bottom_vdw_radius_A"], t))
    elif isinstance(mode_t, str) and mode_t.lower() == "span":
        t, t_info = span, {"method": "span", "z_span_A": span}
        print("[WARN] LAYER_THICKNESS='span'：只用核-核跨度 %s Å，"
              "偏小 -> C 偏大 -> 迁移率偏高" % t)
    elif mode_t is None:
        t, t_info = span, {"method": "span_legacy", "z_span_A": span}
        print("[WARN] LAYER_THICKNESS=None：用核-核跨度 %s Å。"
              "建议改成 \"vdw\" 或手填。" % t)
    else:
        t, t_info = mode_t, {"method": "manual"}
        print("[OK] LAYER_THICKNESS 手填 %s Å" % t)

    if t is None or float(t) <= 0.1:
        print("[WARN] 2D 体系但定不出层厚，弹性常数未做 c/t 修正。"
              "请手填 LAYER_THICKNESS。")
        return True, elastic, c_len
    t = float(t)
    if t <= 0 or t >= c_len:
        sys.exit("[ERROR] LAYER_THICKNESS=%s Å 不合理（c=%.3f Å）" % (t, c_len))

    factor = c_len / t
    new_elastic = rescale_elastic_2d(elastic, factor)
    if elastic is None:
        print("[WARN] 2D：没读到弹性常数，c/t=%.3f 的修正无从施加" % factor)
    else:
        print("[OK] 2D 弹性常数重标度：c=%.3f Å / t=%.3f Å -> ×%.3f"
              % (c_len, t, factor))
    # patch_2d_zero_outofplane_shear（2026-09-18，用户批准）：
    #   二维 slab 的**面外剪切**是真空伪影（CrS2 C44=-0.52、SS 未重排时 C66 槽 0.56、
    #   Mo2S3 C44=-40.2）→ AMSET 的 Christoffel 特征值可为负/极小，形变势因子 deform²/c
    #   被污染（实测 ADP 空穴 2e4 或 0）。二维弹性文献的标准做法是只报面内 C11/C12/C66，
    #   ★ 只置零**负值**（正的即使很小也保留：面内 q 上 C44/C55=0 会让该分支特征值为 0，
    #   AMSET 的 factor 可能出现 0/0 -> NaN）。面内块完全不动；每个被置零的原值都打印。
    zeroed = []
    if new_elastic is not None:
        _M = [list(map(float, r)) for r in new_elastic]
        _seen = set()
        for _i in (3, 4):                      # 标准 Voigt 索引：3=YZ, 4=XZ
            for _j in range(6):
                for _a, _b in ((_i, _j), (_j, _i)):
                    if (_a, _b) in _seen:
                        continue
                    _seen.add((_a, _b))
                    _v = _M[_a][_b]
                    if _v < 0.0:                     # 只清负的（正的小剪切保留，避免 0/0）
                        zeroed.append("C%d%d=%.4f" % (_a + 1, _b + 1, _v))
                        _M[_a][_b] = 0.0
        new_elastic = _M
        if zeroed:
            print("[OK] 2D 面外**负**剪切置零（真空伪影；面内 C11/C12/C22/C66 不动；"
                  "正的剪切保留以免 0/0）：%s" % "、".join(zeroed))
        else:
            print("[..] 2D 面外剪切没有负值，无需置零（正的保留）")
    # patch_2d_thickness：把全部依据落盘，写论文直接引用
    t_info["thickness_used_A"] = round(t, 4)
    write_2d_record(cwd / OUTDIR_NAME, {
        "cell_c_A": round(c_len, 4),
        "layer_thickness": t_info,
        "elastic_rescale_factor_c_over_t": round(factor, 4),
        "elastic_constant_raw_GPa": elastic,
        "elastic_constant_rescaled_GPa": new_elastic,
        "elastic_constant_convention": "standard Voigt (XX YY ZZ YZ XZ XY)，"
                                       "已由 read_elastic 从 VASP 的 (XX YY ZZ XY YZ ZX) 重排",
        "elastic_outofplane_shear_zeroed": zeroed,
        "elastic_outofplane_note": ("二维只用面内 2x2；面外剪切（YZ/XZ 行列）已置零，"
                                    "原值见 elastic_outofplane_shear_zeroed。"
                                    "★ 标准 step8_amset 路径对二维是近似；准确的二维处理请走 "
                                    "step8.4_amset2d 插件（面内 2x2 Christoffel）。"),
        "areal_density_factor_cm": c_len * 1e-8,
        "areal_density_note": "n_2D [cm^-2] = n_3D [cm^-3] x cell_c [cm]",
        "free_carrier_screening": bool(FREE_CARRIER_SCREENING_2D),
        "limitation": ("POP 用三维 Fröhlich 形式；二维极性材料的 Fröhlich "
                       "耦合 q 依赖不同，属定性偏差。结果只可横向比较，"
                       "绝对值不可信。可发表的迁移率请用 Perturbo / EPW。"),
    })
    return True, new_elastic, c_len


# === patch_amset_doping：settings.yaml 要数值列表，不能是 "a:b:n" 冒号语法 ===
def _seg_expand(a, b, n, log):
    import math
    n = int(round(float(n)))
    if n <= 1:
        return [float(b)]
    if log:  # 对数均布，保号（负=n型电子，正=p型空穴）
        s = -1.0 if a < 0 else 1.0
        la, lb = math.log10(abs(a)), math.log10(abs(b))
        return [s * 10 ** (la + (lb - la) * i / (n - 1)) for i in range(n)]
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def expand_spec(spec, log):
    """把 "a:b:n, c:d:m" 范围串或裸数值展开成 float 列表。已是数值则原样收下。"""
    out = []
    for seg in str(spec).split(","):
        seg = seg.strip()
        if not seg:
            continue
        p = seg.split(":")
        if len(p) == 3:
            out += _seg_expand(float(p[0]), float(p[1]), p[2], log)
        elif len(p) == 1:
            out.append(float(p[0]))
        else:
            sys.exit("[ERROR] 无法解析范围段 %r（应为 start:stop:count 或单个数）" % seg)
    return out



# === patch_amset_pop：DFPT OUTCAR 的 Γ 声子频率 -> pop_frequency（THz）===
def read_pop_frequency(dielect_dir):
    """用 AMSET 官方 phonon-frequency 命令计算有效 POP 频率。"""
    import subprocess
    commands = [["amset", "phonon-frequency", "-o", "OUTCAR", "-v", "vasprun.xml"]]
    # 远端 gen 由 taskflow 的 Python 直接执行，非交互 shell 未必加载 conda；
    # 使用配置约定的 amset_clean 作为无 shell 的兜底，避免误报“缺少 POP”。
    commands.append(["/opt/miniconda3/bin/conda", "run", "--no-capture-output", "-n", "amset_clean",
                     "amset", "phonon-frequency", "-o", "OUTCAR", "-v", "vasprun.xml"])
    commands.append(["conda", "run", "--no-capture-output", "-n", "amset_clean",
                     "amset", "phonon-frequency", "-o", "OUTCAR", "-v", "vasprun.xml"])
    for command in commands:
        try:
            p = subprocess.run(command, cwd=str(dielect_dir), capture_output=True,
                               text=True, timeout=300)
        except (OSError, subprocess.TimeoutExpired):
            continue
        text = (p.stdout or "") + "\n" + (p.stderr or "")
        vals = re.findall(r"pop_frequency\s*[:=]\s*([0-9.]+)\s*THz", text, re.I)
        if p.returncode == 0 and vals:
            return round(float(vals[-1]), 4)
    return None


def _is_nonpolar():
    """单元素体系=非极性：无极性光学声子、离子介电≈0。从 POSCAR 元素行(第6行)判断。"""
    try:
        syms = (Path.cwd() / "POSCAR").read_text(errors="ignore").splitlines()[5].split()
        if syms and not any(ch.isdigit() for ch in syms[0]):
            return len(set(syms)) == 1
    except (OSError, IndexError):
        pass
    return False


def _amset_scatterers():
    """校验 POP 所需输入；非极性体系按物理口径移除 POP。"""
    if "POP" in SCATTERING and _is_nonpolar():
        print("[WARN] 单元素非极性体系无极性光学声子——跳过 POP 散射")
        return [s for s in SCATTERING if s != "POP"]
    if "POP" in SCATTERING and read_pop_frequency(Path.cwd() / DIELECT_DIR) is None:
        raise ValueError("AMSET 有效 POP 频率缺失；检查环境和 DFPT 文件，禁止静默删除 POP")
    return list(SCATTERING)



# === patch_2d_dielec：2D 面内介电 + 扣真空 ===
def _grab_diag3(txt, tag):
    """取 tag 后第一块 3x3 的对角 [xx, yy, zz]；取不到返回 None。"""
    i = txt.rfind(tag)
    if i < 0:
        return None
    rows = []
    for ln in txt[i:].splitlines()[1:]:
        nums = re.findall(r"-?\d+\.\d+", ln)
        if len(nums) >= 3:
            rows.append([float(x) for x in nums[:3]])
        if len(rows) == 3:
            break
    return [rows[0][0], rows[1][1], rows[2][2]] if len(rows) == 3 else None


def _dielectric_2d_inplane(dielect_dir, out_dir, eps_inf_fb, eps_static_fb):
    """2D 面内介电 (εx+εy)/2 + 扣真空 ε^m = 1 + (L/t)(ε_sup − 1)。
    L/t 读 2d_correction.json 的 elastic_rescale_factor_c_over_t（与弹性同款）。
    任一步失败/离子介电为负/非极性 → 安全退回原值或仅电子项。"""
    import json as _json
    try:
        rec = _json.loads((Path(out_dir) / "2d_correction.json").read_text())
        factor = float(rec["elastic_rescale_factor_c_over_t"])
    except (OSError, KeyError, ValueError, TypeError):
        print("[WARN] 2D 介电扣真空：读不到 c/t（2d_correction.json）——保持原始三对角平均")
        return eps_inf_fb, eps_static_fb
    try:
        txt = (Path(dielect_dir) / "OUTCAR").read_text(errors="ignore")
    except OSError:
        return eps_inf_fb, eps_static_fb
    inf = _grab_matrix3(txt, "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)")
    ion = _grab_matrix3(txt, "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION")
    if inf is None:
        return eps_inf_fb, eps_static_fb
    stat = [[inf[i][j] + ion[i][j] for j in range(3)] for i in range(3)] if ion else inf
    corr = lambda mat: [[((1 if i == j else 0) + factor * (mat[i][j] - (1 if i == j else 0)))
                         for j in range(3)] for i in range(3)]
    return corr(inf), corr(stat)


def _grab_matrix3(txt, tag):
    """取 tag 后完整 3x3 张量。"""
    i = txt.rfind(tag)
    if i < 0:
        return None
    rows = []
    for ln in txt[i:].splitlines()[1:]:
        nums = re.findall(r"-?\d+(?:\.\d+)?", ln)
        if len(nums) >= 3:
            rows.append([float(x) for x in nums[:3]])
        if len(rows) == 3:
            return rows
    return None
def write_settings(out: Path, eps_inf, eps_static, gap, elastic,
                   is_2d=False, c_len=None):
    lines = ["# amset settings.yaml（gen_step10 自动生成，可手改后重跑本步）",
             "# skill_rev: %s" % _SKILL_REV]
    if is_2d:
        # patch_2d_amset：把 2D 的三件事和一条局限写进文件头，
        # 免得几个月后拿着 transport.json 忘了这是 slab 模型跑出来的。
        lines += [
            "# ===== 2D slab 模型 =====",
            "# 1) elastic_constant 已按 c/t 重标度（还原成层材料的等效三维值）",
            "#    顺序为标准 Voigt (XX YY ZZ YZ XZ XY)，已从 VASP 顺序重排",
            "# 2) free_carrier_screening 已打开（否则 POP 主导时迁移率对浓度不响应）",
            "# 3) 下面的 doping 是体浓度 cm^-3；面浓度 n_2D = n_3D × c",
        ]
        if c_len:
            lines.append("#    本体系 c = %.4f Å = %.4e cm" % (c_len, c_len * 1e-8))
            for v in _doping_endpoints(DOPING):
                lines.append("#      %.3e cm^-3  ->  %.3e cm^-2"
                             % (v, v * c_len * 1e-8))
        lines += [
            "# 【局限】POP 用的是三维 Fröhlich 形式。二维极性材料的 Fröhlich",
            "#   耦合 q 依赖不同，这是定性偏差而非标定误差：结果只可横向比较，",
            "#   绝对值不可信。需要可发表的迁移率请用 Perturbo / EPW。",
            "# ========================",
        ]
    _dop = expand_spec(DOPING, log=True)        # doping 对数均布
    _tmp = expand_spec(TEMPERATURES, log=False)  # 温度线性
    lines += ["doping: [%s]" % ", ".join("%.6e" % v for v in _dop),
              "temperatures: [%s]" % ", ".join("%g" % v for v in _tmp),
              "scattering_type: [%s]" % ", ".join(_amset_scatterers()),
              "deformation_potential: deformation.h5"]
    if INTERPOLATION_FACTOR:     # patch_interp_factor
        lines.append("interpolation_factor: %d" % int(INTERPOLATION_FACTOR))
    # patch_unity_overlap_3d（V25）+ patch_unity_overlap_override（V26）：2D 默认 unity、3D 默认真实重叠，
    #   但 step.conf 的 UNITY_OVERLAP=true/false 可以显式覆盖（用于受控对照）。
    _uo_ov = None if UNITY_OVERLAP == "auto" else (UNITY_OVERLAP == "true")
    if is_2d or _uo_ov is not None:
        _val = UNITY_OVERLAP_2D if _uo_ov is None else _uo_ov
        lines.append("unity_overlap: %s" % ("true" if _val else "false"))
        if _uo_ov is None:
            lines.append("# ^ 2D：去对称化会把真实重叠算坏（V23），一律用 unity_overlap。")
        else:
            lines.append("# ^ step.conf 显式覆盖 UNITY_OVERLAP=%s（受控对照：同一份 settings 只翻这一项）。"
                         % UNITY_OVERLAP)
            # 让 overlap_preflight 知道这是"显式放行的受控对照"（否则 2D + real 会被 ⓪ 拦下）。
            lines.append("# AZ_OVERLAP_CONTROLLED=1")
    else:
        # ★ 必须**显式写 false**：AMSET 默认本就是 False（真实重叠），但只有写出来，
        #   8.3 的 overlap_ratio_guard 才能把这次运行识别为"真实重叠"并给三维标黄；
        #   不写 -> unity_overlap=None -> 8.3 报"0 unity / 0 真实重叠，跳过"（三维黄色告警不会触发）。
        lines.append("unity_overlap: false")
        lines.append("# ^ 三维显式写 false = 真实重叠（与 AMSET 默认一致；写出来是为了让 8.3 能识别）。")
        lines.append("# [WARN] 三维受影响面已确认；Si 全网格对照（V26）：ADP/overall 只差 3~9%、")
        lines.append("#   IMP 逐机制高 1.3~1.7 倍 -> 结论值可用，逐机制值需注明。")
        lines.append("#   若要给某项目改 unity：在本步 step.conf 写 UNITY_OVERLAP = true")
        lines.append("#   （注意 Si 有 6 个谷，unity 最多低估约 6 倍）。")
    # ★ 必须显式写 nworkers：不写则 AMSET 默认 -1 = 用满节点全部核，
    #   会超订 submit 模板申请的核数（见文件头 NWORKERS 处的实测记录）。
    lines.append("nworkers: %d" % int(NWORKERS))
    _popf = read_pop_frequency(Path.cwd() / DIELECT_DIR)
    if _popf is not None and "POP" in SCATTERING:
        lines.append("pop_frequency: %s" % _popf)
    if eps_inf is not None:
        lines.append("high_frequency_dielectric: %s" % eps_inf)
    # patch_nonpolar_dielec：单元素=非极性，离子介电物理上≈0 → 静态介电=ε∞
    if _is_nonpolar() and eps_inf is not None:
        eps_static = eps_inf
    if eps_static is not None:
        lines.append("static_dielectric: %s" % eps_static)
    if gap is not None:
        lines.append("bandgap: %s" % gap)
    if elastic is not None:
        if isinstance(elastic, (int, float)):
            lines.append("elastic_constant: %s" % elastic)   # 各向同性标量
        else:                                                # 6x6 Cij
            lines.append("elastic_constant:")
            for row in elastic:
                lines.append("  - [%s]" % ", ".join("%g" % v for v in row))
    else:
        lines.append("# elastic_constant: 未读到——step6_elastic 没算完，或手填 "
                     "MANUAL_ELASTIC。ACD 散射需要它。")
    if is_2d and FREE_CARRIER_SCREENING_2D:
        lines.append("free_carrier_screening: true")
    (out / "settings.yaml").write_text("\n".join(lines) + "\n",
                                       encoding="utf-8", newline="\n")
    ela = ("标量%s" % elastic if isinstance(elastic, (int, float))
           else ("6x6" if elastic else "无"))
    print("[OK] settings.yaml：ε∞=%s ε₀=%s gap=%s 弹性=%s"
          % (eps_inf, eps_static, gap, ela))


# patch_dim_guard：本步不跑 VASP、也不解析结构，所以没有 dim 变量可用。
# 直接从 step1 的 workflow_method.txt 读 DIM=，0D 就带原因退出，
# 免得 -f 强推时抛一句看不懂的"缺 xxx.h5"。
_STEP1_CANDS = ("step1_opt", "step1_std_opt",
                "step1c_PBE_opt", "step1b_PBE_opt", "step1a_PBE_opt")


def _guard_not_0d(cwd, step_name, why):
    from pathlib import Path as _P
    for name in _STEP1_CANDS:
        mf = _P(cwd) / name / "workflow_method.txt"
        if not mf.is_file():
            continue
        for ln in mf.read_text(errors="ignore").splitlines():
            if ln.strip().upper().startswith("DIM="):
                dim = ln.split("=", 1)[1].strip().lower()
                if dim == "0d":
                    sys.exit("[ERROR] %s 不支持 0D 体系。\n"
                             "        原因：%s\n"
                             "        支持的维度：2D, 3D\n"
                             "        若判定有误，检查 %s 的 DIM=。"
                             % (step_name, why, mf))
                return
        return


# ---------------------------------------------------------------------------
# 提交分配 -> AMSET 进程数
# ---------------------------------------------------------------------------
_SBATCH_NTASKS = re.compile(r"^\s*#SBATCH\s+--ntasks(?:-per-node)?=(\d+)", re.M)
_SBATCH_CPUS = re.compile(r"^\s*#SBATCH\s+--cpus-per-task=(\d+)", re.M)


def _alloc_cores(cwd: Path):
    """本次作业真正拿到的核数 = ntasks_per_node × cpus_per_task，返回 (cores|None, 来源)。

    优先级：step.conf 的 [submit] 覆盖（apply_submit 会按它改写 submit.sh）> 提交模板。
    模板里 OMP/OPENBLAS/MKL_NUM_THREADS 都被压到 1，所以可用核数就是"任务数 × 每任务核数"。
    """
    here = Path(__file__).resolve().parent
    tpl = next((p for p in (here / "submit_amset.tpl", Path(cwd) / "submit_amset.tpl")
                if p.is_file()), None)
    if tpl is None:
        return None, "找不到 submit_amset.tpl"
    text = tpl.read_text(encoding="utf-8", errors="ignore")
    mn, mc = _SBATCH_NTASKS.search(text), _SBATCH_CPUS.search(text)
    ntasks = int(mn.group(1)) if mn else 1
    cpus = int(mc.group(1)) if mc else 1
    src = tpl.name
    try:
        sub = stepconf.read_submit(stepconf.CONF_NAME, used_incar=True)
    except Exception:                                       # noqa: BLE001
        sub = {}
    for key in ("ntasks_per_node", "cpus_per_task"):
        v = str(sub.get(key) or "").strip()
        if v.isdigit():
            if key == "ntasks_per_node":
                ntasks = int(v)
            else:
                cpus = int(v)
            src += " + step.conf[submit].%s" % key
    return max(1, ntasks * cpus), "%s（%d × %d）" % (src, ntasks, cpus)


# ---------------------------------------------------------------------------
# patch_overlap_preflight（V24 / V25.10）：重叠路径的运行前检查
#   ① gen 时（本机）：只有生成的输入，能查的是维度/重叠模式/版本 -> 挡住 2D+真实重叠；
#   ② 作业内：AMSET_CMD 里先跑 --in-job，h5 完整性 / S3-S3b 能带一致性 / 带窗口全查。
# ---------------------------------------------------------------------------
PREFLIGHT_SRC_NAME = "overlap_preflight.py"


def _install_preflight(out):
    """把 overlap_preflight.py 复制进运行目录（作业内 --in-job 要用）。"""
    here = Path(__file__).resolve().parent
    src = next((p for p in (here / PREFLIGHT_SRC_NAME, Path.cwd() / PREFLIGHT_SRC_NAME)
                if p.is_file()), None)
    if src is None:
        print("[WARN] 找不到 %s —— 作业内的运行前检查会跳过（gen 时仍检查）"
              % PREFLIGHT_SRC_NAME)
        return False
    import shutil
    dst = Path(out) / PREFLIGHT_SRC_NAME
    if src.resolve() != dst.resolve():
        shutil.copyfile(src, dst)
    return True


def _preflight_gate(cwd, out, unity):
    """gen 时的运行前检查：能读到的项当场查，读不到的（无 h5/vasprun）自动跳过。"""
    try:
        import overlap_preflight as _pf
    except Exception as e:                                   # noqa: BLE001
        print("[WARN] 运行前检查脚本不可用（%s: %s）——跳过" % (type(e).__name__, e))
        return
    try:
        verdict, lines = _pf.run(cwd, out, unity)
    except Exception as e:                                   # noqa: BLE001
        print("[WARN] 运行前检查失败（%s: %s）——不拦截，请人工确认"
              % (type(e).__name__, e))
        return
    print("[..] 重叠路径运行前检查：")
    for _l in lines:
        print("[..] " + _l)
    if verdict == "error":
        sys.exit("[ERROR] 重叠路径运行前检查未通过 —— 见上面标 ★ 的行。"
                 "\n        这一步会产出不可信的结果，必须先修（V24/V25）。")
    if verdict == "warn":
        print("[WARN] 重叠路径运行前检查有告警（见上）——不拦截，请人工确认。")


def main():
    _disc_gate()
    cwd = Path.cwd()
    global LAYER_THICKNESS, NWORKERS, WAVEFUNCTION_FULL, UNITY_OVERLAP
    _conf_nworkers = None
    if (cwd / "step.conf").is_file():
        # strict=False：材料级 step.conf 是【全技能共用】的一份，含别的步骤的键
        # （FUNC / CELL_POLICY 等）。严格模式会被这些键打死，而它抛的是
        # SystemExit=BaseException —— 原来这里 except SystemExit: pass 把整段
        # 吞掉，于是本函数对 step.conf 的所有覆盖（LAYER_THICKNESS / NWORKERS）
        # **永远静默失效**。与 2026-09-16 S3_uniform 被 FUNC 打死是同一个坑。
        # 另外 StepConf 没有 .get()，原来写的 _p.get("NWORKERS") 一旦真跑到
        # 就是 AttributeError；这里改用 _p[...]（SPEC 里已给默认值）。
        try:
            _p = stepconf.load(SPEC, STEP, str(cwd), strict=False)
            if _p["LAYER_THICKNESS"]:
                LAYER_THICKNESS = _p["LAYER_THICKNESS"]
            # NWORKERS 允许按材料/步骤覆盖；不写（None）= 自动，见下面。
            if _p["NWORKERS"]:
                _conf_nworkers = int(_p["NWORKERS"])
            # 全网格 h5（三维单变量对照）：显式开关，默认 False。
            if _p["WAVEFUNCTION_FULL"]:
                WAVEFUNCTION_FULL = True
                print("[..] WAVEFUNCTION_FULL=true：wavefunction.h5 ← step4b_wave_full"
                      "（全网格，ISYM=-1），vasprun.xml ← step3b_uniform_full")
            # unity_overlap 显式覆盖（受控对照）：auto/true/false
            _uo = str(_p["UNITY_OVERLAP"]).strip().lower()
            if _uo in ("true", "1", "yes", "on"):
                UNITY_OVERLAP = "true"
                print("[..] UNITY_OVERLAP=true（step.conf 覆盖）：强制 unity_overlap: true")
            elif _uo in ("false", "0", "no", "off"):
                UNITY_OVERLAP = "false"
                print("[..] UNITY_OVERLAP=false（step.conf 覆盖）：强制真实重叠")
            elif _uo != "auto":
                print("[WARN] UNITY_OVERLAP=%r 不认识（只认 auto/true/false），按 auto 处理" % _uo)
        except (KeyError, ValueError, TypeError):
            pass  # step.conf 读不成时保持出厂默认（LAYER_THICKNESS="vdw" / NWORKERS 自动）

    # ---- NWORKERS：默认自动 = 本次提交实际分配到的核数 ----
    # 写死一个常数在 jzzn/hanhai25/3090 上必然错两台（少要 = 算力闲置，
    # 多要 = 抢同节点别人的 CPU），所以默认跟着提交分配走。
    _alloc, _src = _alloc_cores(cwd)
    if _conf_nworkers:
        NWORKERS = _conf_nworkers
        if _alloc and NWORKERS != _alloc:
            print("[WARN] step.conf 写死 NWORKERS=%d，但本次提交分配 %d 核（%s）——"
                  "少要 = 算力闲置，多要 = 抢同节点别人的 CPU，请核对这两个数。"
                  % (NWORKERS, _alloc, _src), file=sys.stderr)
        else:
            print("[..] NWORKERS = %d（step.conf 指定）" % NWORKERS)
    elif _alloc:
        NWORKERS = _alloc
        print("[OK] NWORKERS = %d（自动按提交分配：%s）" % (NWORKERS, _src))
    else:
        NWORKERS = NWORKERS_FALLBACK
        print("[WARN] 推不出提交分配（%s）—— NWORKERS 兜底 %d" % (_src, NWORKERS))
    _guard_not_0d(cwd, "step8_amset",
                  "载流子输运建立在能带色散和布里渊区积分上，孤立分子两者都没有")
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    _wdir = "step4b_wave_full" if WAVEFUNCTION_FULL else WAVE_DIR
    link(out, cwd / _wdir / "wavefunction.h5", "wavefunction.h5")
    link(out, cwd / READ_DIR / _pick_deformation_h5(cwd, READ_DIR), "deformation.h5")
    # patch_amset_vasprun：amset run 还需要密网格 vasprun.xml 拿能带色散
    _vdirs = (("step3b_uniform_full",) if WAVEFUNCTION_FULL
              else ("step3_uniform", "step4_wave"))
    _vr = next((cwd / _d / "vasprun.xml"
                for _d in _vdirs
                if (cwd / _d / "vasprun.xml").is_file()), None)
    if _vr is None:
        sys.exit("[ERROR] 找不到 vasprun.xml（%s 没跑完？）" % " / ".join(_vdirs))
    link(out, _vr, "vasprun.xml")
    eps_inf, eps_static = read_dielectric(cwd / DIELECT_DIR)
    gap = read_bandgap(cwd)
    elastic = read_elastic(cwd)
    # patch_2d_amset：2D 时重标度弹性常数并记下 c，供面浓度换算
    is_2d, elastic, c_len = apply_2d_corrections(cwd, elastic)
    if is_2d and c_len and TWO_D_DIELECTRIC_VACUUM:  # patch_2d_dielec
        eps_inf, eps_static = _dielectric_2d_inplane(
            cwd / DIELECT_DIR, out, eps_inf, eps_static)
    if gap is None:
        _msg = ("没读到带隙——settings.yaml 将不写 bandgap，AMSET 会退回 "
                "step3_uniform 的裸 DFT(PBE) 带隙。PBE 带隙偏小会让 300 K 双极导通"
                "占主导，表现为：p 型 Seebeck 变负、σ 随掺杂非单调、Lorenz 数是 "
                "Sommerfeld 的数倍。这种结果不可用。\n"
                "        排查顺序：\n"
                "          1) step2_bandgap/step2.3_hse_plot/band_summary.json 在不在？"
                "里面 gap_eV 是多少？\n"
                "          2) 带隙段关掉了的话，在 project_setting/setting.yaml 写 "
                "bandgap: <eV>\n"
                "          3) 或直接在本脚本顶部写死 MANUAL_BANDGAP = <eV>\n"
                "          4) 确属金属/半金属，把 REQUIRE_BANDGAP 设为 False")
        if REQUIRE_BANDGAP:      # patch_require_bandgap
            sys.exit("[ERROR] %s" % _msg)
        print("[WARN] %s" % _msg)
    if elastic is None:
        print("[WARN] 没读到弹性常数——step6_elastic 没算完，或手填 MANUAL_ELASTIC。"
              "ACD 声学散射需要它。")
    write_settings(out, eps_inf, eps_static, gap, elastic,
                   is_2d=is_2d, c_len=c_len)
    # patch_overlap_preflight（V24/V25.10）：gen 时的配置闸 + 把检查脚本带进运行目录
    _install_preflight(out)
    _preflight_gate(cwd, out, UNITY_OVERLAP_2D if is_2d else False)

    # tf 把 submit_amset.tpl 与本脚本一起推到 gen 运行目录，但按原名推、不会改成
    # submit.sh；本步自己把它渲染成 out/submit.sh（维度步靠各自 gen 的 render，这里同理）。
    here = Path(__file__).resolve().parent
    tpl = next((p for p in (here / "submit_amset.tpl", cwd / "submit_amset.tpl")
                if p.is_file()), None)
    if tpl is None:
        sys.exit("[ERROR] 找不到 submit_amset.tpl（gen_need 里要有它，且应随 gen 脚本一起推送）")
    submit = out / "submit.sh"
    jobname = ("%s-ke-dft-cpu-%s" % (cwd.name, STEP_LABEL)) if not _HAS_KC \
        else kc.new_jobname(cwd, STEP_LABEL)
    text = tpl.read_text(encoding="utf-8")
    text = text.replace("{{JOBNAME}}", jobname).replace("{{AMSET_CMD}}", AMSET_CMD)
    submit.write_text(text, encoding="utf-8", newline="\n")
    stepconf.apply_submit(submit, stepconf.read_submit(stepconf.CONF_NAME))
    print("[DONE] %s：settings.yaml + 软链就绪，提交后 amset run 产出 transport.json"
          % OUTDIR_NAME)


if __name__ == "__main__":
    main()
