#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step14_amset2d.py —— 二维散射核版 AMSET 输入（step8.4_amset2d）。

【来源】本文件由 skill/ke-dft-cpu/step8_amset/gen_step10_amset.py 复制改写。
两份文件的公共部分需保持同步；二维专属逻辑只在本文件里。

相对 gen_step10_amset.py 的差别（其余逻辑逐行一致）：
  1. 维度闸门：DIM != 2d 直接报错退出（三维请用 step8_amset）。
  2. 弹性常数：写**原始 slab 值**（标准 Voigt 顺序，已从 VASP 顺序重排），**不乘 c/t**。
     理由：AMSET 在超胞里正确的输入就是 C_2D/c；乘了 c/t 会把 ADP 散射率压低 t/c 倍。
     c/t 仍然记录在 2d_correction.json 里，供 σ/κ 的厚度归一化用。
  3. 2d_correction.json 扩充：layer_normal / 完整 3×3 eps_inf_slab 与 eps_static_slab /
     eps_env / imp_distance_A / pop_lo_dispersion / mechanisms_2d（插件读这些）。
  4. pop_frequency：取 Γ 点**面内极性模**的有效频率（amset2d_pop_freq.py），
     不再用 AMSET 自带的三维球面平均（那会把面外 ZO 模算进去）。
  5. 启动命令：先把 amset2d_plugin.py 复制到运行目录，再 import amset2d_plugin 后跑
     Runner —— 插件只在本步生效，不改 AMSET 安装，三维项目完全不受影响。

输入软链（与 step8_amset 相同）：
  wavefunction.h5   ← step4_wave
  deformation.h5    ← step7b_deform_read
  vasprun.xml       ← step3_uniform
产出目录：step8.4_amset2d/，产物 transport.json（判据看 thermal_conductivity）。
"""
import glob
import os
import re
import shutil
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

# [SKILL_REV] 版本戳：写进 settings.yaml 头，便于从结果反查跑的是哪份 skill 副本。
_SKILL_REV = "2026-09-16-amset2d"
# 复制来源（gen_step10_amset.py 的 _SKILL_REV），用于人工比对两份文件是否同步
_SOURCE_GEN_REV = "2026-09-16-voigt-reorder"
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
        if not _dc.gate(_P.cwd(), "step8.4_amset2d", force=FORCE_TRANSPORT):
            _s.exit("[BLOCKED] 体系判别为 SEMIMETAL/METAL@PBE —— step8.4_amset2d 已阻断。"
                    "完整提示见 step2_bandgap/step2.15_discriminant/discriminant.json "
                    "的 block.hint；强制继续请把本脚本顶部 FORCE_TRANSPORT 设为 True。")
    except SystemExit:
        raise
    except Exception as _e:
        print("[WARN] 体系判别闸门异常，放行：%s" % _e, file=_s.stderr)


OUTDIR_NAME = "step8.4_amset2d"
WAVE_DIR    = "step4_wave"
READ_DIR    = "step7b_deform_read"
DIELECT_DIR = "step5_dielect"
# patch_ke_dag：跟上 v1.9 的目录重命名。HSE 优先，其次 PBE，最后兼容老目录。
BANDGAP_PLOT_CANDS = [
    "step2_bandgap/step2.3_hse_plot",
    "step2_bandgap/step2.2_pbe_plot",
    "step4_band_plot",
]
STEP_LABEL  = "S8.4_amset2d"
# 启动命令：先 import 插件（改写 amset 的散射注册表），再按原版入口跑 Runner。
# 用 python -c 而不是 amset run，是因为插件必须在 Runner 读取 settings 之前生效。
# patch_overlap_preflight（V24/V25.10）：作业内先跑运行前检查（--in-job）——
#   本目录里 wavefunction.h5 / vasprun.xml / ../step3_uniform 都在，四项检查全生效；
#   不过就 exit 1。脚本由 _install_preflight 复制进运行目录。
AMSET_CMD   = ('python overlap_preflight.py --in-job || exit 1; '
               'python -c "import amset2d_plugin; '
               'from amset.log import initialize_amset_logger as L; L(); '
               'from amset.core.run import Runner; Runner.from_directory(\'.\').run()" '
               '>> amset.log 2>&1 && cp -f "$(ls -t transport_*.json 2>/dev/null | head -1)" '
               'transport.json && ls -l transport.json')
# --- 输运设置（可改）---
DOPING      = "-1e21:-1e17:5, 1e17:1e21:5"   # n 型 + p 型各 5 点（对数均布）cm^-3
TEMPERATURES = "100:900:9"            # 100,200,...,900 K，每 100 K 一个点
SCATTERING  = ["ADP", "IMP", "POP"]   # 形变势声学 + 电离杂质 + 极性光学 "ADP", "IMP", "POP"
_ALLOWED_SCATTERING = {"ADP", "IMP", "POP", "PIE"}   # step.conf SCATTERING 允许的机制
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
# patch_unity_overlap（2026-09-17，VERIFICATION V22/V23）：
#   显式把 unity_overlap 写进 settings.yaml。**以前这一行是缺的** —— AMSET 的默认值是
#   False（真实重叠），于是所有 2D 项目默认踩在"去对称化把重叠算坏"这个坑上。
#   裁决结论：真实重叠必须让 h5 走 from_data（见 wavefunction_full 分支），
#   否则 ADP 迁移率会被抬高 10~30 倍。所以在拿到全网格波函数之前，这里固定 True。
#   要真实重叠：打开 optional_steps.wavefunction_full，并把这个常量改成 False。
#   也可由 step.conf 覆盖：UNITY_OVERLAP = true/false（见 SPEC / main 的 conf 块）。
UNITY_OVERLAP = True
# patch_wavefunction_full（2026-09-20，二维版；照搬 gen_step10_amset.py 的三维逻辑）：
#   可选改读**全网格** h5（step4b_wave_full，ISYM=-1 写出全部 k 点），让 AMSET 走
#   from_data、跳过去对称化（TR×R bug 见 VERIFICATION V33/V37）。打开方式：项目
#   project_setting 里写 wavefunction_full: true 或本步 step.conf 写 WAVEFUNCTION_FULL = true。
#   h5 与 vasprun 必须同源（都取 step3b/step4b），否则 h5 点数与 vasprun 网格对不上。
#   默认 False = 仍读 step4_wave，二维既有行为完全不变。
WAVEFUNCTION_FULL = False
# --- 弹性常数来源（amset run 的 ACD 散射需要）---
#   MANUAL_ELASTIC 填了就用它，否则从 ELASTIC_DIR/OUTCAR 自动解析（kBar→GPa）。
#   直接填：单个数（各向同性近似，GPa），或 6x6 列表（完整 Cij，GPa）。
MANUAL_ELASTIC = None
ELASTIC_DIR = "step6_elastic"
# patch_deform_ref：形变势参考口径（2D 默认真空）。step7b 写出 deformation_vac.h5 就用它，
# 否则退回 deformation.h5（AMSET 芯态对齐）。实测 CrS2：3.53 -> 5.86 eV，μ 差 2.76 倍。
DEFORM_REF = "auto"      # auto = vacuum 优先（本步只服务 2D）；或强制 core 作对比
# --- patch_2d_amset：2D 修正 ---------------------------------------------
# TWO_D_MODE: "auto" = 读 step1 的 workflow_method.txt 的 DIM=；"on"/"off" 强制。
TWO_D_MODE = "on"   # 本步只服务 2D：main() 里另有 DIM 硬闸门（非 2d 直接退出）
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
# NWORKERS：AMSET 的并行进程数，写进 settings.yaml 的 nworkers。
#   ★ 2026-09-16：本步此前**完全不写 nworkers**（gen_step10 已写、这里漏了），
#     AMSET 取默认 -1，其源码 (amset/interpolation/bandstructure.py:182) 把 -1
#     解释成 multiprocessing.cpu_count() = **用满节点全部核** ——
#     jzzn 提交只申请 24 核，于是超订 8 倍，抢同节点其他作业的 CPU。
#   ★ 默认 None = **自动**：由提交模板（叠加 step.conf 的 [submit] 覆盖）推出
#     本次真正分配到的核数 ntasks_per_node × cpus_per_task。
#     三台机器口径不同：jzzn 24×1=24、hanhai25 128×1=128、3090 1×8=8 ——
#     写死常数必然错两台（少要 = 算力闲置，多要 = 抢别人的 CPU）。
#   ★ 要固定值就在 step.conf 写 NWORKERS = <int>，与分配不一致时打印 WARN。
#   与 gen_step10_amset.py 的同名逻辑保持同步（改一处请改两处）。
NWORKERS = None
NWORKERS_FALLBACK = 24      # 提交模板里读不到 SLURM 分配（--ntasks/--cpus）时的兜底
STEP = "step8_amset"
SPEC = {
    "LAYER_THICKNESS": (LAYER_THICKNESS, "str"),
    # NWORKERS 必须在 SPEC 里声明，否则 step.conf 里写了会被判成"不认识的键"。
    "NWORKERS": (NWORKERS, "int"),
    # 全网格 h5 分支（二维真实重叠的唯一干净路径）；默认 False，既有二维项目行为不变。
    "WAVEFUNCTION_FULL": (WAVEFUNCTION_FULL, "bool"),
    # unity_overlap 覆盖（受控对照/真实重叠用）：auto/true/false；默认 auto = 现行行为。
    "UNITY_OVERLAP": ("auto", "str"),
    # 插值因子覆盖（2026-09-20）：AMSET 收敛判据在插值后的网格上。
    #   稠密网格 ~ (factor*nk)^2*nkz，内存/耗时按此平方级增长 ——
    #   实测 factor=10 + nworkers=24 在共享节点上 MaxRSS 292 GB 被 OOM 杀。
    #   项目里按需降到 4（已校准口径：与 factor 10 差 6-11%，见 V23）。
    "INTERPOLATION_FACTOR": (INTERPOLATION_FACTOR, "int"),
    # 散射类型覆盖（2026-09-22 用户批准）：逗号/空格分隔，如 `SCATTERING = ADP,POP` 用于严格本征对照；
    # 空串 = 出厂 ["ADP","IMP","POP"] 行为不变。
    "SCATTERING": ("", "str"),
}
# --- amset2d 插件相关（可改）---
PLUGIN_SRC_NAME = "amset2d_plugin.py"   # 与本脚本同目录，运行时复制到 OUTDIR_NAME
POP_FREQ_HELPER = "amset2d_pop_freq.py"  # Γ 点面内极性模频率的取法（见该文件头）
# 该 helper 的完整 JSON（含 Δr 两种算法交叉核对），供 2d_correction.json 落盘
_LAST_POP_JSON = {}
# 手填 pop_frequency（THz）。None = 自动（面内极性模加权平均）。已知文献值时优先手填。
POP_FREQUENCY_MANUAL = None
# 二维散射机制就在上面的 SCATTERING 里选（同时写进 2d_correction.json 的 mechanisms_2d，
# 插件按它替换注册表）。PIE 需要 settings 里有 piezoelectric_constant，即 step6_elastic 的
# DFPT 必须输出了 PIEZOELECTRIC TENSOR；没有就在 _mechanisms_2d() 里自动去掉并告警。
# 环境介电（上下包覆介质；悬空 = 1.0）。hBN 包覆填 ~4.5，见文档 T6。
EPS_ENV = 1.0
# 杂质到层中心的距离（Å）；层内杂质 = 0。
IMP_DISTANCE_A = 0.0
# 用二维 LST 关系把 ω_TO 修正成 ω_LO(q)：True 时 pop_frequency 取面内极性 TO 模。
POP_LO_DISPERSION = True
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
    if any(eps_0[i][i] < 0 for i in range(3)):
        raise ValueError(
            "离子介电对角项为负（eps_ionic=%.2f:%.2f:%.2f，eps_static=%.2f），非物理，"
            "多因 Gamma 近零声学模污染 DFPT 离子介电；必须检查 DFPT 收敛与声子稳定性。"
            % (eps_0[0][0], eps_0[1][1], eps_0[2][2], eps_static[0][0]))
    _coup = [1.0 / eps_inf[i][i] - 1.0 / eps_static[i][i] for i in range(3)]
    if max(abs(c) for c in _coup) < 1e-9:
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


def _pick_deformation_h5(cwd, read_dir):
    """挑形变势文件：默认用真空口径 deformation_vac.h5（二维文献口径），缺失则 core。

    DEFORM_REF: "auto"（vacuum 优先）/ "vacuum" / "core"。口径写进 settings 注释，
    因为两种口径的绝对值差 ~2.8 倍，必须能从产物反查。
    """
    read_dir = Path(read_dir)
    ref = str(DEFORM_REF).lower()
    if ref not in ("auto", "vacuum", "core"):
        sys.exit("[ERROR] DEFORM_REF=%r 无效，只允许 auto / vacuum / core" % DEFORM_REF)
    want_vacuum = ref in ("auto", "vacuum")
    if want_vacuum and (read_dir / "deformation_vac.h5").is_file():
        print("[OK] 形变势参考口径：**真空**（deformation_vac.h5，二维文献口径）")
        return "deformation_vac.h5"
    if want_vacuum:
        print("[WARN] 想要真空口径但 %s/deformation_vac.h5 不存在 -> 退回 **core** 口径"
              "（AMSET 芯态对齐，为体材料设计）。要真空口径请重跑 step7b"
              "（需要 step7_deform 形变目录里的 LOCPOT）。" % read_dir)
    else:
        print("[..] 形变势参考口径：core（按 DEFORM_REF=%s）" % ref)
    return "deformation.h5"


def vasp_to_voigt_6x6(m):
    """6×6 弹性/刚度张量：VASP 顺序 (XX YY ZZ XY YZ ZX) → 标准 Voigt (XX YY ZZ YZ XZ XY)。

    只在 OUTCAR 表头认不出来时的兜底路径用；正常走 _reorder_voigt（按表头实测顺序）。"""
    a = [[float(x) for x in row] for row in m]
    return [[a[i][j] for j in _VASP2VOIGT] for i in _VASP2VOIGT]


def vasp_to_voigt_piezo(e3x6):
    """3×6 压电张量（VASP 列序 XX YY ZZ XY YZ ZX）→ 标准 Voigt 列序（只重排列）。"""
    a = [[float(x) for x in row] for row in e3x6]
    return [[row[j] for j in _VASP2VOIGT] for row in a]


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
        except Exception:
            pass
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


# --------------------------------------------------------------------------
# amset2d 插件需要的输入：层法向 / 原始 slab 介电 / 面内极性模频率 / 压电张量
# --------------------------------------------------------------------------
def _layer_normal(cwd: Path):
    """层法向 = a×b 的单位向量（取 +z 那一支）。

    插件要求它与笛卡尔 z 平行（不平行直接报错）。这里如实落盘，便于插件核对与排查。
    """
    for d in STRUCT_CANDS:
        for fn in ("CONTCAR", "POSCAR"):
            p = cwd / d / fn
            if not p.is_file():
                continue
            try:
                lines = p.read_text(encoding="utf-8-sig",
                                    errors="ignore").splitlines()
                scale = float(lines[1].split()[0])
                a = [float(x) * scale for x in lines[2].split()[:3]]
                b = [float(x) * scale for x in lines[3].split()[:3]]
                n = [a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]]
                norm = sum(x * x for x in n) ** 0.5
                if norm < 1e-9:
                    continue
                v = [round(x / norm, 6) for x in n]
                if v[2] < 0:
                    v = [-x for x in v]
                return v
            except Exception as e:                          # noqa: BLE001
                print("[WARN] 解析 %s 的晶格失败：%s" % (p, e))
    print("[WARN] 读不到结构，layer_normal 默认 [0, 0, 1]（插件会按实际值核对）")
    return [0.0, 0.0, 1.0]


def _slab_tensors(dielect_dir, eps_inf_fb, eps_static_fb):
    """DFPT OUTCAR 里**原始 slab**（含真空稀释）的完整 3×3 ε∞ 与 ε₀。

    插件要的是没扣真空、没做面内平均的原值：二维极化长度张量 r = c(ε_∥ − 1)/2
    只有用原始张量才对。读不到就退回 read_dielectric 已解析的值。
    """
    try:
        txt = (Path(dielect_dir) / "OUTCAR").read_text(errors="ignore")
    except OSError:
        return eps_inf_fb, eps_static_fb
    inf = _grab_matrix3(txt, "MACROSCOPIC STATIC DIELECTRIC TENSOR "
                             "(including local field effects in DFT)")
    ion = _grab_matrix3(txt, "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION")
    if inf is None:
        print("[WARN] 读不到完整 3×3 eps_inf，插件改用 read_dielectric 的值")
        return eps_inf_fb, eps_static_fb
    if ion is None or any(ion[i][i] < 0 for i in range(3)):
        print("[WARN] 离子介电缺失或对角为负 -> eps_static_slab 退回 eps_inf_slab"
              "（POP 的 1/ε∞ − 1/ε₀ 会变小甚至为 0，注意核对）")
        return inf, inf
    stat = [[inf[i][j] + ion[i][j] for j in range(3)] for i in range(3)]
    return inf, stat


# === patch_delta_r_check：r∞ / Δr 自检（2026-09-16 新增）===
# 定义（与插件、METHODOLOGY §4 一致）：r = c(ε_slab,∥ − 1)/2（真空收敛后与 c 无关），
# 插件 POP 核用的是 G_2D ∝ 2πc·Δr/[(1 + r∞q)(1 + r0q)]，其中
#     Δr = r0 − r∞ = c(ε0,slab − ε∞,slab)/2 。
# 把这两个数落盘并打印：曾经把**体相** ε0/ε∞ 外推出来的 Δr（43.5 Å）当成 Sohier 的
# 0.53 Å 用，单机制迁移率因此差了约两个量级，而只看 ε 张量根本看不出来。
# 锚点（文献）：MoS2 r∞≈46.5 Å / Δr≈0.53 Å（Δr/r∞≈1.1%）；h-BN r∞≈7.6 Å / Δr≈2.8 Å（≈37%）。
_R_ANCHORS = {
    "MoS2": {"r_inf_A": 46.5, "delta_r_A": 0.53},
    "h-BN": {"r_inf_A": 7.64, "delta_r_A": 2.84},
}


def _inplane_mean(mat):
    """3×3 张量的面内对角平均；读不到/非有限返回 None。"""
    try:
        a = [[float(mat[i][j]) for j in range(2)] for i in range(2)]
    except (TypeError, ValueError, IndexError):
        return None
    vals = [a[0][0], a[1][1]]
    if any(v != v or v in (float("inf"), float("-inf")) for v in vals):
        return None
    return (vals[0] + vals[1]) / 2.0


def delta_r_check(c_len, inf33, stat33):
    """由**原始 slab** 介电张量算 r∞ / r0 / Δr（Å），打印自检并返回落盘的 dict。"""
    ei, es = _inplane_mean(inf33), _inplane_mean(stat33)
    if ei is None or es is None:
        print("[WARN] 2D 介电张量不可用 —— r∞/Δr 自检跳过"
              "（POP 的绝对量级没有依据，结果只能看趋势）")
        return {"available": False}
    r_inf = c_len * (ei - 1.0) / 2.0
    r_0 = c_len * (es - 1.0) / 2.0
    dr = r_0 - r_inf
    ratio = (dr / r_inf) if r_inf else float("inf")
    out = {
        "available": True,
        "eps_inf_inplane_mean": round(ei, 4),
        "eps_static_inplane_mean": round(es, 4),
        "r_inf_A": round(r_inf, 4),
        "r_0_A": round(r_0, 4),
        "delta_r_A": round(dr, 4),
        "delta_r_over_r_inf": round(ratio, 5),
        "definition": "r = c*(eps_slab,inplane - 1)/2 (Angstrom); "
                      "Delta r = r0 - r_inf = c*(eps0-eps_inf)/2",
        "anchors": _R_ANCHORS,
    }
    print("[OK] 2D 极化长度：r∞ = %.3f Å，Δr = %.4f Å（Δr/r∞ = %.1f%%）"
          % (r_inf, dr, 100 * ratio))
    print("     锚点参考：MoS2 r∞≈46.5 Å / Δr≈0.53 Å（1.1%）；"
          "h-BN r∞≈7.6 Å / Δr≈2.8 Å（37%）")
    if dr < -1e-9:
        print("[WARN] Δr < 0（ε0 < ε∞）：离子介电缺失或读错，POP 强度会变成 0 或负值，"
              "请核对 %s/OUTCAR 的 IONIC CONTRIBUTION 块" % DIELECT_DIR)
    elif ratio > 0.6:
        print("[WARN] Δr/r∞ = %.0f%% 远大于文献锚点量级 —— 最常见的原因是把**体相** "
              "ε0/ε∞ 外推当作 slab 值（曾出现 43.5 Å vs 0.53 Å）。"
              "请确认 ε0/ε∞ 取自同一次 DFPT 的原始 slab 张量。" % (100 * ratio))
    if r_inf <= 0 or r_inf > 2000:
        print("[WARN] r∞ = %.1f Å 超出合理范围（0-2000 Å），请核对 ε∞ 张量" % r_inf)
    return out


def _mechanisms_2d(cwd: Path):
    """插件要替换的散射机制（落盘给插件读）与 settings 的 scattering_type 一致。"""
    mechs = list(SCATTERING)
    if "PIE" in mechs and read_piezo(cwd) is None:
        print("[WARN] settings 里要 PIE，但 step6_elastic/OUTCAR 没有 PIEZOELECTRIC TENSOR"
              "——本次去掉 PIE（DFPT 需 LEPSILON + IBRION=6）")
        mechs = [m for m in mechs if m != "PIE"]
    return mechs


def read_piezo(cwd: Path):
    """从 ELASTIC_DIR/OUTCAR 读**弛豫离子**压电张量 e（3×6，C/m²），重排成标准 Voigt 列序。

    VASP 的列序与弹性一致（XX YY ZZ XY YZ ZX），AMSET 按标准 Voigt 解析，所以要重排。
    优先取含 'including ionic contribution' 的总张量块；否则取最后一个块。
    """
    oc = cwd / ELASTIC_DIR / "OUTCAR"
    if not oc.is_file():
        return None
    txt = oc.read_text(errors="ignore")
    idxs = [m.start() for m in re.finditer("PIEZOELECTRIC TENSOR", txt)]
    if not idxs:
        return None
    chosen, tag = None, None
    for i in reversed(idxs):
        head = txt[i:i + 200].splitlines()[0]
        rows = []
        for ln in txt[i:].splitlines()[1:]:
            nums = re.findall(r"-?\d+\.\d+", ln)
            if len(nums) >= 6:
                rows.append([float(x) for x in nums[:6]])
            if len(rows) == 3:
                break
        if len(rows) == 3:
            chosen, tag = rows, head.strip()
            if "including ionic" in head.lower():
                break
    if chosen is None:
        return None
    print("[OK] 压电张量来源：%s" % tag)
    return vasp_to_voigt_piezo(chosen)


def read_pop_frequency_2d(dielect_dir):
    """Γ 点**面内极性模**有效频率（THz）：走 amset2d_pop_freq.py（需要 amset 环境）。

    AMSET 自带的 phonon-frequency 对三维球面方向加权，slab 里会把面外极性模算进去
    （实测 CrS2 单层：3D 球面 12.45 THz / 面内 12.38 THz；旧流程误取最高模 15.07 THz）。
    返回 None 表示取不到——调用方会退回三维口径并告警。
    """
    import json as _json
    import subprocess
    # 本进程只需要跑一次 helper（它的 JSON 里还带 Δr 两路线交叉核对，要给 2d_correction.json 用）。
    # ★ 顺序坑（2026-09-16 MoS2 实测）：_amset_scatterers() 里第一次调用发生在
    #   apply_2d_corrections() **之后**，那时 2d_correction.json 已经写完了，
    #   交叉核对字段会静默缺失。这里缓存住，并在 apply_2d_corrections 里先调一次。
    if _LAST_POP_JSON.get("pop_frequency_THz"):
        return round(float(_LAST_POP_JSON["pop_frequency_THz"]), 4)
    here = Path(__file__).resolve().parent
    helper = here / POP_FREQ_HELPER
    if not helper.is_file():
        print("[WARN] 找不到 %s（gen_need 里要有它）" % POP_FREQ_HELPER)
        return None
    # 用**绝对路径**调 helper：cwd 必须是 step5_dielect（OUTCAR/vasprun.xml 在那儿），
    # 而 helper 自己躺在 gen 脚本目录里。
    script = ("%s '%s'" % (helper.name, str(Path(dielect_dir))))
    commands = [["python", str(helper), "OUTCAR", "vasprun.xml"],
                ["/opt/miniconda3/bin/conda", "run", "--no-capture-output",
                 "-n", "amset_clean", "python", str(helper), "OUTCAR", "vasprun.xml"],
                ["conda", "run", "--no-capture-output", "-n", "amset_clean",
                 "python", str(helper), "OUTCAR", "vasprun.xml"]]
    err_tail = ""
    for command in commands:
        try:
            p = subprocess.run(command, cwd=str(dielect_dir), capture_output=True,
                               text=True, timeout=600)
        except (OSError, subprocess.TimeoutExpired) as _e:
            err_tail = str(_e)
            continue
        err_tail = (p.stderr or "").strip().splitlines()[-1:] or [""]
        err_tail = err_tail[0]
        out = (p.stdout or "").strip().splitlines()
        for ln in reversed(out):
            ln = ln.strip()
            if not ln.startswith("{"):
                continue
            try:
                d = _json.loads(ln)
            except ValueError:
                continue
            if d.get("pop_frequency_THz"):
                _LAST_POP_JSON.clear()
                _LAST_POP_JSON.update(d)
                print("[OK] 二维 pop_frequency=%.4f THz（面内极性模 %s；"
                      "三维对照 %.4f）" % (d["pop_frequency_THz"],
                                           d.get("polar_inplane_modes"),
                                           d.get("pop_frequency_3d_ref_THz") or float("nan")))
                _drc = d.get("delta_r_check") or {}
                if _drc.get("route1_dielectric_A") is not None:
                    print("[OK] Δr 交叉核对：介电路线 %.6f Å vs 声子模求和 %.6f Å"
                          "（偏差 %s%%；剔除声学模 %s 个）；r∞ = %.3f Å（文献 MoS2 "
                          "Δr≈0.53 / r∞≈46.5 Å 作量级参考）"
                          % (_drc["route1_dielectric_A"], _drc["route2_mode_sum_A"],
                             _drc.get("deviation_pct"), _drc.get("acoustic_removed"),
                             _drc.get("r_inf_A") or float("nan")))
                elif _drc.get("error"):
                    print("[WARN] Δr 交叉核对失败：%s" % _drc["error"])
                return round(float(d["pop_frequency_THz"]), 4)
    print("[WARN] %s 没给出结果（命令：%s；末行错误：%s）"
          % (POP_FREQ_HELPER, script, err_tail))
    return None


def apply_2d_corrections(cwd: Path, elastic, eps_inf=None, eps_static=None):
    """二维几何/弹性/介电参数落盘。

    返回 (is_2d, elastic_写进 settings 的值, c_len)。
    ★ 与 step8_amset 的关键差别：**elastic 原样返回**（原始 slab 值，不乘 c/t）。
      AMSET 的散射积分在超胞里做，二维核的正确输入就是 C_2D/c = 原始 slab 值；
      乘 c/t 会让 ADP 散射率偏小 t/c 倍（本步插件按 G_s = c·G_2D 自己换算）。
      c/t 仍写进 2d_correction.json，供 σ/κ 的厚度归一化用。
    """
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

    factor = c_len / t                     # 只用于记录与 σ/κ 的厚度归一化
    t_info["thickness_used_A"] = round(t, 4)
    if elastic is None:
        print("[WARN] 2D：没读到弹性常数——ADP 散射无法计算（settings 里不写 "
              "elastic_constant）")
    else:
        print("[OK] 2D 弹性常数：用**原始 slab 值**（标准 Voigt 顺序），"
              "不乘 c/t=%.3f——二维核按 G_s = c·G_2D 自己换算" % factor)

    # 插件需要的**原始 slab** 介电张量：完整 3×3（含非对角，斜方/单斜要用），不扣真空
    inf33, stat33 = _slab_tensors(cwd / DIELECT_DIR, eps_inf, eps_static)
    layer_normal = _layer_normal(cwd)
    mechanisms = _mechanisms_2d(cwd)
    # 先把 Γ 点声子信息拿到：pop_frequency 与 Δr 交叉核对都在同一个 helper 的 JSON 里，
    # 而 settings 那边的调用发生在写 2d_correction.json 之后（顺序坑，见 read_pop_frequency_2d）。
    if "POP" in mechanisms:
        try:
            read_pop_frequency_2d(cwd / DIELECT_DIR)
        except Exception as _e:                                # noqa: BLE001
            print("[WARN] 预取 pop_frequency/Δr 交叉核对失败（不影响本步）：%s" % _e)
    # [2026-09-16] 形变势构型口径：step7b 把实际用的那一套写进 band_edges.json，
    # 这里读回来落盘/落注释 —— 回退到 clamped 时不能静默（用户要求）。
    global _DEFORM_GEOM
    try:
        import json as _json
        _be = cwd / "step7b_deform_read" / "band_edges.json"
        if _be.is_file():
            _DEFORM_GEOM = (_json.loads(_be.read_text()) or {}).get("deform_geometry")
    except Exception as _e:                                    # noqa: BLE001
        print("[WARN] 读 deform_geometry 失败（不影响计算）：%s" % _e)
    if _DEFORM_GEOM:
        print("[OK] 形变势构型口径 deform_geometry=%s（来自 step7b/band_edges.json）"
              % _DEFORM_GEOM)
        if _DEFORM_GEOM != "ionrelax":
            print("[WARN] 形变势不是离子弛豫口径（%s）—— 与 TOTAL ELASTIC MODULI"
                  "（含离子弛豫）不同源，ADP 绝对值跨项目不可直接比。" % _DEFORM_GEOM)
    r_check = delta_r_check(c_len, inf33, stat33)
    # 把 amset2d_pop_freq.py 的两路线 Δr 交叉核对挂进同一条记录（同一次运行拿到）
    if isinstance(r_check, dict) and _LAST_POP_JSON.get("delta_r_check"):
        r_check["delta_r_crosscheck"] = _LAST_POP_JSON["delta_r_check"]

    rec = {
        "cell_c_A": round(c_len, 4),
        "layer_normal": layer_normal,
        "layer_thickness": t_info,
        "elastic_rescale_factor_c_over_t": round(factor, 4),
        "elastic_constant_raw_GPa": elastic,
        "elastic_constant_rescaled_GPa": None,
        "elastic_constant_convention": "standard Voigt (XX YY ZZ YZ XZ XY)，"
                                       "已由 read_elastic 从 VASP 的 (XX YY ZZ XY YZ ZX) 重排",
        "elastic_constant_note": "插件用原始 slab 值（= C_2D/c），不乘 c/t；"
                                 "c/t 只用于 σ/κ 的厚度归一化（见 elastic_rescale_factor_c_over_t）",
        "eps_inf_slab": inf33,
        "eps_static_slab": stat33,
        "r_inf_delta_r_check": r_check,
        "deform_geometry": _DEFORM_GEOM,
        "eps_env": float(EPS_ENV),
        "imp_distance_A": float(IMP_DISTANCE_A),
        "pop_lo_dispersion": bool(POP_LO_DISPERSION),
        "mechanisms_2d": mechanisms,
        "areal_density_factor_cm": c_len * 1e-8,
        "areal_density_note": "n_2D [cm^-2] = n_3D [cm^-3] x cell_c [cm]",
        "free_carrier_screening": bool(FREE_CARRIER_SCREENING_2D),
        "skill_rev": _SKILL_REV,
        "plugin": PLUGIN_SRC_NAME,
        "limitation": ("单一有效极性模、偶极近似（无四极矩、无面外偶极/Janus 耦合）；"
                       "屏蔽用局域场 + 二维 Thomas–Fermi。模型级近似，"
                       "适合筛选与趋势；要发表的绝对值建议用 EPW/Perturbo 对照。"),
    }
    write_2d_record(cwd / OUTDIR_NAME, rec)
    return True, elastic, c_len


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
def read_pop_frequency_3d(dielect_dir):
    """AMSET 官方 phonon-frequency 命令算出的有效 POP 频率（三维球面口径）。

    二维步**只在**面内口径取不到时退回这里（并在日志里告警），默认用
    read_pop_frequency_2d() 的 Γ 点面内极性模频率。
    """
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


def read_pop_frequency(dielect_dir):
    """二维 pop_frequency：优先 Γ 点**面内极性模**频率，取不到退回三维口径并告警。

    POP_FREQUENCY_MANUAL 非 None 时直接用它（已知文献值时优先手填）。
    """
    if POP_FREQUENCY_MANUAL is not None:
        print("[OK] pop_frequency 手填 %.4f THz" % float(POP_FREQUENCY_MANUAL))
        return round(float(POP_FREQUENCY_MANUAL), 4)
    f2 = read_pop_frequency_2d(dielect_dir)
    if f2 is not None:
        return f2
    f3 = read_pop_frequency_3d(dielect_dir)
    if f3 is None:
        return None
    print("[WARN] 面内极性模频率取不到，退回 AMSET 三维口径 pop_frequency=%.4f THz"
          "（slab 里会把面外极性模算进权重，谨慎使用）" % f3)
    return f3


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
    """settings.yaml 的 scattering_type；校验 POP 所需输入，非极性体系移除 POP。"""
    mechs = list(SCATTERING)
    if "POP" in mechs and _is_nonpolar():
        print("[WARN] 单元素非极性体系无极性光学声子——跳过 POP 散射")
        mechs = [s for s in mechs if s != "POP"]
    if "POP" in mechs and read_pop_frequency(Path.cwd() / DIELECT_DIR) is None:
        raise ValueError("AMSET 有效 POP 频率缺失；检查环境和 DFPT 文件，禁止静默删除 POP")
    return mechs



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
                   is_2d=False, c_len=None, piezo=None):
    lines = ["# amset settings.yaml（gen_step14_amset2d 自动生成，可手改后重跑本步）",
             "# skill_rev: %s（复制自 gen_step10_amset.py @ %s）"
             % (_SKILL_REV, _SOURCE_GEN_REV)]
    if is_2d:
        # patch_2d_amset：把 2D 的口径与局限写进文件头，
        # 免得几个月后拿着 transport.json 忘了这是 slab 模型 + 二维散射核跑出来的。
        lines += [
            "# ===== 2D slab 模型 + 二维散射核（amset2d 插件）=====",
            "# 1) elastic_constant = **原始 slab 值**（标准 Voigt 顺序，已从 VASP 顺序重排），",
            "#    不乘 c/t：AMSET 的散射积分在超胞里做，二维核的正确输入就是 C_2D/c。",
            "#    插件按 G_s = c·G_2D 换算，见 skill/ke-dft-cpu/METHODOLOGY.md 第 4 节。",
            "# 2) 运行必须带 amset2d_plugin（submit.sh 里已 import），散射核为二维形式；",
            "#    二维核参数在 2d_correction.json（原始 slab 介电张量等）。",
            "# 3) free_carrier_screening 已打开（否则 POP/IMP 主导时迁移率对浓度不响应）",
            "#   形变势参考口径：%s（vacuum = 真空零点/二维文献口径，core = AMSET 芯态对齐，差约 2.8 倍）"
            % ("vacuum" if (str(DEFORM_REF).lower() in ("auto", "vacuum")
                           and (Path.cwd() / READ_DIR / "deformation_vac.h5").is_file())
               else "core"),
            "# 4) 下面的 doping 是体浓度 cm^-3；面浓度 n_2D = n_3D × c",
        ]
        if c_len:
            lines.append("#    本体系 c = %.4f Å = %.4e cm" % (c_len, c_len * 1e-8))
            for v in _doping_endpoints(DOPING):
                lines.append("#      %.3e cm^-3  ->  %.3e cm^-2"
                             % (v, v * c_len * 1e-8))
        lines += [
            "# 【局限】单一有效极性模；偶极近似（无四极矩；Poncé PRL 2023 给 MoS2 室温",
            "#   Hall 迁移率误差：电子 23%、空穴 76%）；",
            "#   忽略 Janus 结构的面外偶极耦合；屏蔽为局域场 + 二维 Thomas–Fermi",
            "#   （文献对比偏低 10–20%）。模型级近似，适合筛选与趋势，",
            "#   要发表的绝对值建议用 EPW / Perturbo（开二维库仑截断）对照。",
            "# ========================",
        ]
    _dop = expand_spec(DOPING, log=True)        # doping 对数均布
    _tmp = expand_spec(TEMPERATURES, log=False)  # 温度线性
    lines += ["# deform_geometry: %s（形变势构型口径；ionrelax=离子弛豫、与 TOTAL ELASTIC"
              " MODULI 同源；clamped=离子固定，跨项目 ADP 不可直接比；mixed=两者混用）"
              % (_DEFORM_GEOM or "未知（step7b/band_edges.json 里没有这个字段）")]
    if _DEFORM_GEOM and _DEFORM_GEOM != "ionrelax":
        lines.append("# [WARN] 本次形变势不是离子弛豫口径 —— 与弹性常数不同源，"
                     "与其它项目比较 ADP 绝对值时先核对这一行。")
    lines += ["doping: [%s]" % ", ".join("%.6e" % v for v in _dop),
              "temperatures: [%s]" % ", ".join("%g" % v for v in _tmp),
              "scattering_type: [%s]" % ", ".join(_amset_scatterers()),
              "deformation_potential: deformation.h5"]
    if INTERPOLATION_FACTOR:     # patch_interp_factor
        lines.append("interpolation_factor: %d" % int(INTERPOLATION_FACTOR))
    # patch_unity_overlap：显式写，别用 AMSET 的默认值（默认是 False = 真实重叠，
    # 而真实重叠要求 h5 是完整网格，否则去对称化会把重叠算坏，见 V22/V23）。
    lines.append("unity_overlap: %s" % ("true" if UNITY_OVERLAP else "false"))
    if not UNITY_OVERLAP:
        lines.append("# ^ 真实重叠：**必须**让 h5 走 from_data（wavefunction_full 分支），"
                     "否则结果不可信（见 overlap_preflight.py 的拦截）")
        # patch_overlap_controlled（2026-09-20，照搬 gen_step10_amset.py L982）：
        #   step.conf 显式写 UNITY_OVERLAP=false 时，让 overlap_preflight 知道这是
        #   "显式放行的受控对照"，把 ⓪ 的拦截降为告警——否则全网格真实重叠也无法跑。
        lines.append("# AZ_OVERLAP_CONTROLLED=1")
    # ★ 必须显式写 nworkers：不写则 AMSET 默认 -1 = 用满节点全部核，
    #   会超订提交模板申请的核数（见文件头 NWORKERS 处的说明）。
    lines.append("nworkers: %d" % int(NWORKERS))
    _popf = read_pop_frequency(Path.cwd() / DIELECT_DIR)
    if _popf is not None and "POP" in SCATTERING:
        lines.append("# pop_frequency = Γ 点面内极性模有效频率（二维口径）；插件再乘")
        lines.append("# √(ε̃0(q)/ε̃∞(q)) 得到 ω_LO(q)（二维 LST）。")
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
    if piezo is not None and "PIE" in SCATTERING:
        lines.append("# piezoelectric_constant = 原始 slab 弛豫离子张量（C/m²），"
                     "标准 Voigt 列序")
        lines.append("piezoelectric_constant:")
        for row in piezo:
            lines.append("  - [%s]" % ", ".join("%g" % v for v in row))
    if is_2d and FREE_CARRIER_SCREENING_2D:
        lines.append("free_carrier_screening: true")
    (out / "settings.yaml").write_text("\n".join(lines) + "\n",
                                       encoding="utf-8", newline="\n")
    ela = ("标量%s" % elastic if isinstance(elastic, (int, float))
           else ("6x6" if elastic else "无"))
    print("[OK] settings.yaml：ε∞=%s ε₀=%s gap=%s 弹性=%s"
          % (eps_inf, eps_static, gap, ela))


# patch_dim_guard：本步不跑 VASP、也不解析结构，所以没有 dim 变量可用。
# 直接从 step1 的 workflow_method.txt 读 DIM=：本步**只服务 2D**，不是 2d 就退出。
_STEP1_CANDS = ("step1_opt", "step1_std_opt",
                "step1c_PBE_opt", "step1b_PBE_opt", "step1a_PBE_opt")


def _guard_dim_2d(cwd):
    """维度硬闸门：DIM != 2d 直接退出（三维请走原版 step8_amset）。

    读不到 DIM=（老项目/手工目录）不拦，只告警——与其它闸门保持"读不到就放行"的口径。
    """
    dim = read_dim(cwd)
    if dim is None:
        print("[WARN] 读不到 DIM=（step1 的 workflow_method.txt）——按 2D 继续，"
              "但层法向/层厚都按 slab 处理，请自行核对结构")
        return
    if dim != "2d":
        sys.exit("[ERROR] step8.4_amset2d 只支持 2D 体系（本次 DIM=%s）。\n"
                 "        三维请用原版 step8_amset（amset run，不受本插件影响）；\n"
                 "        本步的二维核要求层法向 ∥ z，并需要 2d_correction.json。" % dim)


def _install_plugin(out):
    """把 amset2d_plugin.py 复制到运行目录（AMSET 在 cwd 里 import 它）。"""
    here = Path(__file__).resolve().parent
    src = next((p for p in (here / PLUGIN_SRC_NAME, Path.cwd() / PLUGIN_SRC_NAME)
                if p.is_file()), None)
    if src is None:
        sys.exit("[ERROR] 找不到 %s（skill.yaml 的 gen_need 里要有它，"
                 "且应随 gen 脚本一起推送）" % PLUGIN_SRC_NAME)
    dst = Path(out) / PLUGIN_SRC_NAME
    if src.resolve() != dst.resolve():
        shutil.copyfile(src, dst)
    print("[OK] 二维插件就位：%s（只在本步的 amset 运行里生效）" % dst)


def _install_preflight(out):
    """把 overlap_preflight.py 复制进运行目录（作业内 python overlap_preflight.py --in-job）。"""
    here = Path(__file__).resolve().parent
    name = "overlap_preflight.py"
    src = next((p for p in (here / name, Path.cwd() / name) if p.is_file()), None)
    if src is None:
        print("[WARN] 找不到 %s —— 作业内的运行前检查会跳过（gen 时仍检查）" % name)
        return
    dst = Path(out) / name
    if src.resolve() != dst.resolve():
        shutil.copyfile(src, dst)
    print("[OK] 运行前检查脚本就位：%s" % dst)


# ---------------------------------------------------------------------------
# 提交分配 -> AMSET 进程数（与 gen_step10_amset.py 的同名实现保持一致）
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


def main():
    _disc_gate()
    cwd = Path.cwd()
    global LAYER_THICKNESS, NWORKERS, UNITY_OVERLAP, WAVEFUNCTION_FULL
    global INTERPOLATION_FACTOR, SCATTERING
    _conf_nworkers = None
    if (cwd / "step.conf").is_file():
        # strict=False：材料级 step.conf 是【全技能共用】的一份，含别的步骤的键
        # （FUNC / CELL_POLICY 等）；严格模式会被它们打死，而抛的是
        # SystemExit=BaseException —— 原来 except SystemExit: pass 把整段吞掉，
        # 于是 LAYER_THICKNESS 的 step.conf 覆盖**永远静默失效**（同 S3_uniform）。
        try:
            _p = stepconf.load(SPEC, STEP, str(cwd), strict=False)
            if _p["LAYER_THICKNESS"]:
                LAYER_THICKNESS = _p["LAYER_THICKNESS"]
            if _p["NWORKERS"]:
                _conf_nworkers = int(_p["NWORKERS"])
            # patch_interp_factor 覆盖（2026-09-20）：与 NWORKERS 同款，缺省沿用常量。
            if _p["INTERPOLATION_FACTOR"]:
                _conf_interp = int(_p["INTERPOLATION_FACTOR"])
                if _conf_interp != INTERPOLATION_FACTOR:
                    print("[OK] INTERPOLATION_FACTOR = %d（step.conf 覆盖，出厂 %d）"
                          % (_conf_interp, INTERPOLATION_FACTOR))
                INTERPOLATION_FACTOR = _conf_interp
            # 全网格 h5 分支：显式开关，默认 False。
            if _p["WAVEFUNCTION_FULL"]:
                WAVEFUNCTION_FULL = True
                print("[..] WAVEFUNCTION_FULL=true：wavefunction.h5 <- step4b_wave_full"
                      "（全网格，ISYM=-1），vasprun.xml <- step3b_uniform_full")
            # unity_overlap 显式覆盖：auto/true/false
            _uo = str(_p["UNITY_OVERLAP"]).strip().lower()
            if _uo in ("true", "1", "yes", "on"):
                UNITY_OVERLAP = True
                print("[..] UNITY_OVERLAP=true（step.conf 覆盖）：强制 unity_overlap: true")
            elif _uo in ("false", "0", "no", "off"):
                UNITY_OVERLAP = False
                print("[..] UNITY_OVERLAP=false（step.conf 覆盖）：真实重叠（需全网格 h5）")
            elif _uo != "auto":
                print("[WARN] UNITY_OVERLAP=%r 不认识（只认 auto/true/false），按 auto 处理" % _uo)
            # patch_scattering（2026-09-22 用户批准）：step.conf 覆盖散射类型。
            #   逗号/空格分隔，如 SCATTERING = ADP,POP 用于"严格本征"对照；
            #   空串 = 出厂 [ADP, IMP, POP] 行为不变；含不认识机制则整项忽略并告警。
            _sc_raw = _p["SCATTERING"]
            if isinstance(_sc_raw, (list, tuple)):
                _sc_list = [str(x).strip().upper() for x in _sc_raw if str(x).strip()]
            else:
                _sc = str(_sc_raw or "").strip()
                _sc_list = ([x.strip().upper() for x in _sc.replace(",", " ").split() if x.strip()]
                            if _sc else [])
            if _sc_list:
                _bad = [x for x in _sc_list if x not in _ALLOWED_SCATTERING]
                if _bad:
                    print("[WARN] SCATTERING 含不认识的机制 %s（允许 %s），本项忽略，保持出厂 %s"
                          % (_bad, sorted(_ALLOWED_SCATTERING), SCATTERING), file=sys.stderr)
                elif _sc_list != SCATTERING:
                    print("[OK] SCATTERING = %s（step.conf 覆盖，出厂 %s）" % (_sc_list, SCATTERING))
                    SCATTERING = _sc_list
        except (KeyError, ValueError, TypeError):
            pass  # step.conf 读不成时保持出厂默认（LAYER_THICKNESS="vdw" / NWORKERS 自动）

    # ---- NWORKERS：默认自动 = 本次提交实际分配到的核数 ----
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
    _guard_dim_2d(cwd)
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    # 插件随本步复制到运行目录（只在那一次 amset 运行里生效，不改 AMSET 安装）
    _install_plugin(out)
    _install_preflight(out)
    _wdir = "step4b_wave_full" if WAVEFUNCTION_FULL else WAVE_DIR
    link(out, cwd / _wdir / "wavefunction.h5", "wavefunction.h5")
    link(out, cwd / READ_DIR / _pick_deformation_h5(cwd, READ_DIR), "deformation.h5")
    # patch_amset_vasprun：amset run 还需要密网格 vasprun.xml 拿能带色散
    # 全网格分支时 vasprun 必须同源（step3b_uniform_full）。
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
    # patch_2d_amset：落盘二维修正依据（原始 slab 介电/层法向/机制），弹性保持原始 slab 值
    is_2d, elastic, c_len = apply_2d_corrections(cwd, elastic, eps_inf, eps_static)
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
    piezo = read_piezo(cwd) if "PIE" in SCATTERING else None
    if "PIE" in SCATTERING and piezo is None:
        print("[WARN] settings 要 PIE 但没有压电张量——本次退化为不含 PIE 的散射集"
              "（step6_elastic 的 DFPT 需 LEPSILON + IBRION=6）")
    write_settings(out, eps_inf, eps_static, gap, elastic,
                   is_2d=is_2d, c_len=c_len, piezo=piezo)

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
    # ---- patch_overlap_preflight：重叠路径的运行前检查（VERIFICATION V24）----
    # ① 版本+重叠模式 ② h5 完整性（真实重叠 + 非完整网格 -> 拦截）
    # ③ S3/S3b 能带一致性 ④ amset wave 带窗口是否显式
    try:
        import overlap_preflight as _pf
        _verdict, _lines = _pf.run(cwd, out, UNITY_OVERLAP)
        print("[..] 重叠路径运行前检查：")
        for _l in _lines:
            print("[..] " + _l)
        if _verdict == "error":
            sys.exit("[ERROR] 重叠路径运行前检查未通过 —— 见上面标 ★ 的行。"
                     "\n        这一步会产出不可信的结果，必须先修再提交。")
    except SystemExit:
        raise
    except Exception as _e:
        print("[WARN] 重叠运行前检查失败（%s）——不拦截，但请人工确认" % type(_e).__name__)

    print("[DONE] %s：settings.yaml + 软链 + %s 就绪；"
          "submit.sh 里先 import 插件再跑 Runner，产出 transport.json"
          % (OUTDIR_NAME, PLUGIN_SRC_NAME))


if __name__ == "__main__":
    main()
