#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step7_deform.py —— 形变势单点，扇出（step7_deform）。

流程：
  1. POSCAR ← step1_std_opt/CONTCAR
  2. amset deform create → 在本目录生成 undeformed/ + deform-01..NN/（各含形变 POSCAR）
  3. 给每个子目录补 INCAR（渲染 incar_deform_*.tpl）+ KPOINTS + POTCAR + submit.sh
tf 的扇出机制随后对每个 *deform* 子目录各自 sbatch。
产出目录：step7_deform/{undeformed, deform-01, ...}
"""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ke_common as kc
import stepconf  # noqa: E402
from dim_common import require_dim, resolve_tpl  # noqa: E402

# =========================== 可改参数区 ===========================
OUTDIR_NAME  = "step7_deform"
PREV_CANDS   = ["step1_opt", "step1_std_opt"]
DIMENSION    = "auto"
VASPKIT_EXE  = "vaspkit"
KSCHEME      = "2"
KSPACING     = "0.03"
DK_MAX       = 0.05                   # 与 step3_uniform 同口径：面内笛卡尔 k 间距上限，
                                      # 保证 step7 与 step3 的网格一致（否则 band 截断差 1，
                                      # amset run 报 band 数不匹配——LS 2×11 vs 6×41 实锤）
# [patch_vacuum_kz] 与 step3_uniform 同口径：把真空方向 kz 提到该层数。
#   kc.vaspkit_kpoints(..., dim="2d") 内部会把真空轴压回 1；只有 1 层 kz 时
#   AMSET 沿 kz 外推，插值网格一变结果差 30-45%（V8）。step3_uniform 有这个
#   patch，step7 此前漏了 —— 于是同一材料的 S3(47×47×3) 与 S7(15×15×1) 网格
#   不一致，形变势 h5 比波函数网格粗 3 倍/轴。默认 1 = 不动（与旧行为一致）。
VACUUM_KZ_MIN = 1
FUNC         = "inherit"      # patch_ke_dag: inherit=继承 step1
# 注意：本步必须和 step3_uniform 同泛函，否则形变势里会掺进泛函差异
MANUAL_ENCUT = None
ENCUT_FACTOR = 1.5
STEP_LABEL   = "S7_deform"
# amset 环境：优先读 step.conf 里 tf 注入的集群 CONDA_SH/AMSET_ENV（setting/<集群>.yaml
# 里配 conda_sh/amset_env，每人按自己的机器配置，脚本不硬编码路径）；
# step.conf 没有（旧项目/直跑脚本）才回退主机探测。
import os as _os
def _amset_env_src():
    try:
        import stepconf as _sc
        _txt = open(_sc.CONF_NAME, encoding="utf-8-sig").read()
        _p = {k.upper(): v for k, v, _ in _sc.parse(_txt, _sc.CONF_NAME).get("params", [])}
        _sh, _env = _p.get("CONDA_SH"), _p.get("AMSET_ENV")
        if _sh and _env:
            return "source %s && conda activate %s" % (_sh, _env)
    except Exception:
        pass
    if _os.path.isdir("/home/wangchaoyue852/miniconda3"):
        return "source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh && conda activate amset"
    return "source /public/home/wangchao/miniconda3/etc/profile.d/conda.sh && conda activate amset051"  # 2026-09-22 全局切 0.5.1
AMSET_ENV_SRC = _amset_env_src()
DEFORM_GLOB  = "*deform*"      # 与 skill.yaml 的 fanout 一致
                               # （同时匹配 deform-NNN 和 undeformed，两者都要交）
# patch_deform_fix：形变幅度，amset 默认 0.005 = ±0.5%，形变势对 ±0.5% 取平均
STRAIN_DISTANCE = 0.005
# 对称精度：None = amset 默认 0.01；"N" = 完全关对称（12 个形变全算）；
# 也可写数值。放宽它能减少不等价形变数，但要确认放宽后的空间群是对的。
SYMPREC      = None
# =================================================================
GGA_MAP = {"pbe": "PE", "pbesol": "PS", "pbe-d3": "PE"}


def run_amset_deform_create(out: Path):
    """在 out 目录跑 amset deform create，再把产出的 POSCAR-NNN 摊成子目录。

    patch_deform_fix：amset 的 create 子命令
      - 不吃位置参数，结构文件走 -f/--filename；
      - 产出的是【文件】POSCAR-001 / POSCAR-002 …（位数随总数变），
        写在当前目录，【不建任何子目录，也不建 undeformed/】。
    所以目录结构得我们自己摊：
        undeformed/POSCAR   <- 原始（未形变）结构
        deform-001/POSCAR   <- POSCAR-001
        deform-002/POSCAR   <- POSCAR-002
        ...
    skill.yaml 的 fanout "*deform*" 同时匹配 deform-NNN 和 undeformed，
    两者都会被提交——undeformed 是形变势的参考态，不能漏。
    """
    import glob
    import shutil as _sh

    done = [p for p in glob.glob(str(out / "deform-*")) if os.path.isdir(p)]
    if done and (out / "undeformed").is_dir():
        print("[..] 已有形变子目录 %d 个 + undeformed，跳过 create（幂等）"
              % len(done))
        return

    opts = "-f POSCAR -d %s" % STRAIN_DISTANCE
    if SYMPREC is not None:
        opts += " -s %s" % SYMPREC
    cmd = ("%s && cd %s && amset deform create %s >> deform_create.log 2>&1"
           % (AMSET_ENV_SRC, str(out), opts))
    print("[..] amset deform create %s ..." % opts)
    rc = subprocess.run(["bash", "-lc", cmd]).returncode
    if rc != 0:
        sys.exit("[ERROR] amset deform create 失败，看 %s/deform_create.log" % out)

    poscars = sorted(Path(p) for p in glob.glob(str(out / "POSCAR-*")))
    if not poscars:
        sys.exit("[ERROR] amset deform create 没有产出 POSCAR-NNN，"
                 "看 %s/deform_create.log" % out)

    src_poscar = out / "POSCAR"
    if not src_poscar.is_file():
        sys.exit("[ERROR] %s 缺原始 POSCAR，无法建 undeformed/" % out)
    und = out / "undeformed"
    und.mkdir(exist_ok=True)
    _sh.copy2(str(src_poscar), str(und / "POSCAR"))
    print("[OK] undeformed/POSCAR（参考态）")

    for p in poscars:
        tag = p.name.split("-", 1)[1]          # POSCAR-003 -> 003
        d = out / ("deform-%s" % tag)
        d.mkdir(exist_ok=True)
        _sh.copy2(str(p), str(d / "POSCAR"))
    print("[OK] 摊出 %d 个 deform-NNN 子目录（+ undeformed，共 %d 个单点）"
          % (len(poscars), len(poscars) + 1))


# patch_stale_grid（2026-09-21）：ck_deform 只查 OUTCAR 收尾标志 + ionrelax/OUTCAR，
#   **不比对 KPOINTS**。于是改了网格之后，旧网格算完的子目录仍被判"已完成"，
#   不进 fan_todo；而 remote_sbatch_fanout 只提交 fan_todo -> retry 只补没完成的，
#   结果是一个形变势里混着两套 k 网格。这里在覆盖 KPOINTS 之前，把旧产物改名归档
#   （不删除），使 ck_deform 判不通过、子目录自然进入 fan_todo 被重算。
#   目录里的输入（POSCAR/POTCAR/INCAR）不在归档集合内，绝不触碰。
STALE_OUTPUTS = (
    "OUTCAR", "OSZICAR", "vasprun.xml", "CONTCAR", "CHGCAR", "CHG", "EIGENVAL",
    "DOSCAR", "PROCAR", "WAVECAR", "IBZKPT", "REPORT", "PCDAT", "XDATCAR",
    "LOCPOT", "ELFCAR", "OUTCAR.stale-grid",
)


def _grid_tag(kpts_text):
    """从 KPOINTS 文本里取网格标签（vaspkit 写法的第四行 3 个整数），用于归档后缀。"""
    for _ln in (kpts_text or "").splitlines():
        _p = _ln.split()
        if len(_p) == 3 and all(x.isdigit() for x in _p):
            return "x".join(_p)
    return "old"


# === patch_kpts_align：形变势要求所有构型共用同一套 k 点 ===
def _reference_kpoints(out, dim, vac_axis, n_sub):
    """只在 undeformed/ 生成一次 KPOINTS 作为基准，返回其路径（失败返回 None）。

    形变势是逐 k 点、逐能带做差。vaspkit 按各自晶格常数换算 KSPACING，
    ±0.5% 应变足以让某个方向的细分数跳一档，于是 deform-NN 与 undeformed
    的 k 点数对不上，差分失效。统一拷贝基准 KPOINTS 可根除这个问题
    （应变只有 0.5%，共用网格在收敛性上无损）。"""
    und = out / "undeformed"
    if not (und / "POSCAR").is_file():
        print("[WARN] 没有 undeformed/POSCAR —— 退回逐目录生成 KPOINTS，"
              "各形变的 k 网格可能不一致，形变势会失真")
        return None
    # patch_stale_grid_undeformed（2026-09-21）：undeformed/ 的 KPOINTS 就是基准本体，
    #   紧接着的 vaspkit_kpoints 会【直接覆盖】它；而调用点的循环走到 d=undeformed 时
    #   _same=True 会跳过归档 —— 旧网格的 OUTCAR 于是留存、ck_deform 仍判"已完成"，
    #   fan_todo 少一个（9 而非 10）。所以必须在覆盖【之前】抓下旧内容，事后比对归档。
    _und_kp = und / "KPOINTS"
    try:
        _und_kp_old = (_und_kp.read_text(encoding="utf-8", errors="ignore")
                       if _und_kp.is_file() else None)
    except OSError:
        _und_kp_old = None
    kc.vaspkit_kpoints(und, KSCHEME, KSPACING, VASPKIT_EXE, dim, vac_axis)
    # [patch_vacuum_kz] 与 step3_uniform 同口径：把真空方向 kz 提到 VACUUM_KZ_MIN。
    #   必须在 vaspkit_kpoints 之后（它会把 2D 真空轴强制压回 1）。
    _kzmin = int(VACUUM_KZ_MIN)
    #   ★ _reference_kpoints 是模块级函数，拿不到 main() 里的 cwd；
    #     项目目录 = out 的父目录（out = cwd/OUTDIR_NAME）。
    _proj = out.parent
    _kzc = (_proj / "step.conf")
    if _kzc.is_file():
        try:
            _conf = stepconf.load({"VACUUM_KZ_MIN": (VACUUM_KZ_MIN, "int")},
                                  OUTDIR_NAME, str(_proj), strict=False)
            _kzmin = int(_conf["VACUUM_KZ_MIN"])
        except (KeyError, ValueError, TypeError):
            pass
    # [patch_vacuum_kz_warn] 2026-09-20：覆盖是【按步骤名】合并的，S3 写了
    #   VACUUM_KZ_MIN 不等于 S7 也有（MoS2 实测：S3/S3b=3，S7 缺 step.conf -> 静默退回 1）。
    #   静默退回正是本问题的成因，所以这里显式告警：若 S3/S3b 的网格 kz 已 >1 而本步
    #   仍解析到 1，说明缺 templates/step7_deform/step.conf —— 提示去补，而不是闷头算。
    if dim == "2d" and _kzmin <= 1:
        for _cand in ("step3_uniform", "step3b_uniform_full"):
            _ck = _proj / _cand / "KPOINTS"
            try:
                _cz = int(_ck.read_text().splitlines()[3].split()[vac_axis])
            except (OSError, IndexError, ValueError):
                continue
            if _cz > 1:
                print("[WARN] patch_vacuum_kz：%s 的 %s 已是 kz=%d，但本步解析到 kz=1 —— "
                      "多半缺 project_setting/templates/%s/step.conf 里的 VACUUM_KZ_MIN=%d。"
                      "两条命令补：autozt -p <MAT> -j S3_uniform conf --set params.VACUUM_KZ_MIN=%d "
                      "（S3b 同理）与 -j S7_deform。否则 S7 与 S3 网格不一致。"
                      % (OUTDIR_NAME, _cand, _cz, OUTDIR_NAME, _cz, _cz))
                break
    if dim == "2d" and _kzmin > 1:
        _kp = und / "KPOINTS"
        _ln = _kp.read_text().splitlines()
        _n = [int(x) for x in _ln[3].split()[:3]]
        if _n[vac_axis] < _kzmin:
            _old = _n[vac_axis]
            _n[vac_axis] = _kzmin
            _ln[3] = " %d %d %d" % (_n[0], _n[1], _n[2])
            _kp.write_text("\n".join(_ln) + "\n")
            print("[OK] patch_vacuum_kz：真空轴 kz %d -> %d（原来只有一层时 AMSET 沿 kz 外推，"
                  "插值网格一变结果差 30-45%%，且 step7 会与 step3 网格不一致）"
                  % (_old, _kzmin))
    # [DK_MAX] 与 step3_uniform 同口径加密：面内笛卡尔间距上限，逐轴 max。
    # 否则 step7(undeformed 网格) 与 step3 网格不同 → band 截断差 → amset 报错。
    #
    # 2026-09-20 修正：去掉原来的「六方胞豁免」（_ratio >= 1.5 才加密）。
    #   那条豁免与 step3_uniform 其实【不一致】—— gen_step5_uniform.py 对所有
    #   周期性方向都套 DK_MAX（无 hex 豁免），六方胞照样加密。豁免的后果是
    #   同一材料的 S3 与 S7 网格差 3 倍/轴（MoS2: S3=47×47×3 vs S7=15×15×1；
    #   CrSe2_hex: 46 vs 15），形变势 h5 比波函数网格粗，ADP 插值误差可达
    #   几十个百分点（谷内 |dD|/D0 实测 p90≈38%，见 tmp/amset2d/paper/）。
    #   现在与 step3_uniform 完全同口径：逐轴 max(vaspkit, ceil(|b_i|/DK_MAX))，
    #   对【面内两轴】生效；2D 的真空轴已在上面 patch_vacuum_kz 定好，不动 kz。
    #
    #   ★★ 2026-09-20 二次修正：范围必须保持 dim == "2d" ★★
    #   原代码是 `if DK_MAX and dim == "2d"`，3D 分支【故意不加密】——
    #   gen_step5_uniform.py(S3) 对 3D 用 DK_MAX_3D=0.06 且带 2x 静态下限，
    #   逐轴收紧到 ~24x24x24；S7 若也照做，会把 3D 的形变势单点成本抬 27 倍，
    #   且与既有 3D 项目（225 合金等）已跑的 S7 口径不一致 —— 用户明确要求
    #   「改动绝不能破坏 3D 的 S7」。A/B 夹具实测：去掉这个 guard 会让
    #   cubic 3D 从 8x8x8 变 24x24x24（DIFF）。所以这里恢复 dim == "2d"。
    #   3D 的 S7 网格维持原行为（vaspkit 原样），不随本次二维修复改变。
    if DK_MAX and dim == "2d":
        import numpy as np
        _ln = (und / "POSCAR").read_text().splitlines()
        _s = float(_ln[1].split()[0])
        _a = np.array([float(x) for x in _ln[2].split()[:3]]) * _s
        _b = np.array([float(x) for x in _ln[3].split()[:3]]) * _s
        _c = np.array([float(x) for x in _ln[4].split()[:3]]) * _s
        _vol = abs(float(np.dot(_a, np.cross(_b, _c))))
        _rec = [2.0 * np.pi * np.cross(_b, _c) / _vol,
                2.0 * np.pi * np.cross(_c, _a) / _vol,
                2.0 * np.pi * np.cross(_a, _b) / _vol]
        _len = [float(np.linalg.norm(v)) for v in _rec]
        _axes = (0, 1) if dim == "2d" else (0, 1, 2)   # 2D 真空轴已定，不动
        _kpt = (und / "KPOINTS").read_text().splitlines()
        try:
            _n = [int(x) for x in _kpt[3].split()]
        except (IndexError, ValueError):
            _n = [1, 1, 1]
        while len(_n) < 3:
            _n.append(1)
        _need = list(_n[:3])
        for i in _axes:
            _need[i] = max(_n[i], int(np.ceil(_len[i] / float(DK_MAX))))
        if _need != _n[:3]:
            print("[WARN] undeformed 网格 %dx%dx%d（笛卡尔间距 %.3f/%.3f Å⁻¹）不满足 "
                  "DK_MAX=%.3f，按轴提到 %dx%dx%d（与 step3_uniform 同口径）"
                  % (_n[0], _n[1], _n[2], _len[0] / _n[0], _len[1] / _n[1],
                     float(DK_MAX), _need[0], _need[1], _need[2]))
            _kpt[3] = "  %d  %d  %d" % (_need[0], _need[1], _need[2])
            (und / "KPOINTS").write_text(
                "\n".join(_kpt) + "\n", encoding="utf-8", newline="\n")
    # patch_stale_grid_undeformed：基准网格变了 -> undeformed/ 的 OUTCAR 是旧网格产物。
    #   归档旧产物（改名，不删除），使 ck_deform 判不通过、该子目录进入 fan_todo。
    try:
        _und_kp_new = (und / "KPOINTS").read_text(encoding="utf-8", errors="ignore")
    except OSError:
        _und_kp_new = None
    if _und_kp_old is not None and _und_kp_new is not None and _und_kp_old != _und_kp_new:
        _tag = _grid_tag(_und_kp_old)
        _moved = []
        for _nm in STALE_OUTPUTS:
            _f = und / _nm
            if _f.is_file():
                _f.rename(und / (_nm + ".stale-grid-" + _tag))
                _moved.append(_nm)
        print("[..] undeformed: KPOINTS %s -> %s，归档旧产物 %d 个（*.stale-grid-%s）"
              % (_tag, _grid_tag(_und_kp_new), len(_moved), _tag))
    kpts = und / "KPOINTS"
    if not kpts.is_file():
        print("[WARN] undeformed/KPOINTS 生成失败 —— 退回逐目录生成")
        return None
    print("[OK] 基准 k 网格取自 undeformed/，拷贝给全部 %d 个单点目录" % n_sub)
    return kpts


# ===========================================================================
#  [PATCH-IONRELAX] 参数
# ===========================================================================
# 弛豫段的力判据。★ 不要收紧 ★
#   2026-09-13 实测（Mg4C60 单层形变胞，ISYM=0 → 58 个 k 点）：
#     第 5 个离子步起能量已完全平（dE ~ 5E-4 eV/全胞 = 4E-6 eV/atom），
#     残余力 max 0.03–0.05 eV/Å、RMS 0.016–0.027 eV/Å —— 早已优于常规判据。
#   而旧代码写死 -1E-3，比它紧 30–50 倍，逼 VASP 跑满 NSW=60 个离子步
#   （每步 ~48 min），11 个形变胞白烧约 5E4 core·h。
#   力 0.03 eV/Å 对应残余位移 ~0.003 Å，对带边能量影响仅 ~1E-4 eV，
#   而形变势信号在 0.1–1 eV —— 完全淹没。
#   要更紧就改这一行（不要改下面那几处 re.sub）。
IONRELAX_EDIFFG = "-0.02"
IONRELAX_NSW = "60"

# ---- [PATCH-IONRELAX] 对 xx±/yy± 4 个形变目录生成 ionrelax/ 子目录 ----
def _build_ionrelax(d: Path, encut, subs, submit_body):
    """在 deform-NN 下建 ionrelax/：两段式（弛豫段 IBRION=2 + 静态段 LVHAR）。

    E1 必须用「离子弛豫后的内坐标」取带边能量（刚性形变 S 原子不动，E1 系统性
    偏小 24%），amset h5 仍用本级刚性单点（clamped-ion 口径不变）。弛豫段
    EDIFF=1E-6 省时间，静态段读 CONTCAR+WAVECAR 取 EDIFF=1E-8 精确带边 + LOCPOT。
    两段都保持 ISYM=0。
    力判据用 IONRELAX_EDIFFG（默认 -0.02）—— 详见文件头该常量的注释。"""
    ir = d / "ionrelax"
    ir.mkdir(exist_ok=True)
    # 清理旧 SCF 产物（重跑/幂等时残留的 vasprun/LOCPOT/OUTCAR 等），避免
    # step7b 的 _resolved 读到过期结果
    for _f in list(ir.glob("*")):
        if _f.name not in ("POSCAR", "KPOINTS", "POTCAR"):
            if _f.is_file() or _f.is_symlink():
                _f.unlink()
    shutil.copy2(str(d / "POSCAR"), str(ir / "POSCAR"))
    shutil.copy2(str(d / "KPOINTS"), str(ir / "KPOINTS"))
    shutil.copy2(str(d / "POTCAR"), str(ir / "POTCAR"))

    base = (d / "INCAR").read_text(encoding="utf-8")
    # 弛豫段
    relax = re.sub(r"IBRION\s*=.*", "IBRION = 2            # 离子弛豫（固定晶格）", base)
    relax = re.sub(r"NSW\s*=.*", "NSW    = " + IONRELAX_NSW, relax)
    relax = re.sub(r"EDIFF\s*=.*", "EDIFF  = 1E-6", relax)
    relax = re.sub(r"LVHAR\s*=.*", "LVHAR  = .FALSE.", relax)
    relax = re.sub(r"LWAVE\s*=.*", "LWAVE  = .TRUE.", relax)
    relax = re.sub(r"LCHARG\s*=.*", "LCHARG = .FALSE.", relax)
    # 力判据：用 IONRELAX_EDIFFG（默认 -0.02）。两种情形都处理：
    #   段 INCAR 里已有 EDIFFG（继承自模板）-> 覆盖；没有 -> 追加。
    # 旧代码只在"没有"时追加 -1E-3，等于把一个过严值写死、且没有任何配置出口。
    if re.search(r"EDIFFG\s*=", relax):
        relax = re.sub(r"EDIFFG\s*=.*", "EDIFFG = " + IONRELAX_EDIFFG, relax)
    else:
        relax = relax.rstrip() + "\nEDIFFG = " + IONRELAX_EDIFFG + "\n"
    (ir / "INCAR.relax").write_text(relax, encoding="utf-8", newline="\n")
    # 静态段
    stat = re.sub(r"IBRION\s*=.*", "IBRION = -1           # 静态段（弛豫后精确带边）", base)
    stat = re.sub(r"NSW\s*=.*", "NSW    = 0", stat)
    stat = re.sub(r"ISTART\s*=.*", "ISTART = 1            # 接弛豫段 WAVECAR", stat)
    stat = re.sub(r"ICHARG\s*=.*", "ICHARG = 0", stat)
    if "ISTART" not in stat:
        stat = stat.rstrip() + "\nISTART = 1\n"
    (ir / "INCAR.static").write_text(stat, encoding="utf-8", newline="\n")

    # deform-NN 的 submit.sh 追加 ionrelax 两段（跑完刚性单点后 cd ionrelax 续跑）
    sub = (d / "submit.sh").read_text(encoding="utf-8")
    _mpi = None
    for _ln in sub.splitlines():
        _s = _ln.strip()
        # 只认真正的启动命令：模板注释里也含 "mpirun" 字样
        # （"…会让 mpirun 段错误…"）。不过滤注释会把注释当命令行，
        # 导致 ionrelax 段丢失 VASP 调用、只剩 cp CONTCAR 而必然失败。
        if not _s or _s.startswith("#"):
            continue
        if _s.startswith(("mpirun", "srun", "mpiexec")):
            _mpi = _s
            break
    if not _mpi:
        _mpi = "mpirun -np $SLURM_NTASKS vasp_std"
    extra = (
        "\n# ---- ionrelax 两段（E1 用离子弛豫内坐标）----\n"
        "cd ionrelax || exit 1\n"
        "cp INCAR.relax INCAR\n"
        "set -o pipefail\n"
        + _mpi + " 2>&1 | tee vasp.relax.log | tail -20\n"
        "_relax_rc=$?\n"
        "for _f in OUTCAR OSZICAR CONTCAR vasprun.xml; do\n"
        "  if [ -f \"$_f\" ]; then cp -p \"$_f\" \"$_f.relax\" || exit 1; fi\n"
        "done\n"
        "[ \"$_relax_rc\" -eq 0 ] || exit \"$_relax_rc\"\n"
        "set +o pipefail\n"
        "cp CONTCAR POSCAR\n"
        "cp INCAR.static INCAR\n"
        + _mpi + " 2>&1 | tail -20\n"
        "cp CONTCAR POSCAR\n"
    )
    sub = sub.rstrip() + extra
    (d / "submit.sh").write_text(sub, encoding="utf-8", newline="\n")


def main():
    import glob
    cwd = Path.cwd(); out = cwd / OUTDIR_NAME; out.mkdir(exist_ok=True)
    prev = kc.find_prev_dir(cwd, PREV_CANDS)
    if prev is None:
        sys.exit("[ERROR] 找不到含 CONTCAR 的上一步：%s" % PREV_CANDS)
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR", "step1_opt")
    _func, _subs = kc.resolve_func(prev, FUNC, OUTDIR_NAME)
    dim = kc.read_method_dim(prev / kc.METHOD_FILE) \
        or kc.resolve_dim_for(out / "POSCAR", DIMENSION)[0]
    _, vac_axis = kc.resolve_dim_for(out / "POSCAR", dim)
    require_dim(dim, ('2d', '3d'), "step7_deform",
                why="载流子输运/形变势建立在能带色散上，孤立分子没有色散")
    print("[..] 维度：%s" % dim.upper())
    kc.write_method(out / kc.METHOD_FILE, dim, "形变势单点（扇出）",
                    func=_func)

    run_amset_deform_create(out)

    subs = sorted([Path(p) for p in glob.glob(str(out / DEFORM_GLOB)) if os.path.isdir(p)])
    if not subs:
        sys.exit("[ERROR] amset deform create 没有产出形变子目录")
    print("[..] 形变子目录 %d 个，逐个补输入" % len(subs))

    tpl = Path(__file__).resolve().parent / ("incar_deform_%s.tpl" % dim)
    submit_tpl = resolve_tpl(Path(__file__).resolve().parent, "submit_std", dim)
    submit_body = submit_tpl.read_text(encoding="utf-8")

    ref_kpts = _reference_kpoints(out, dim, vac_axis, len(subs))   # patch_kpts_align

    encut = None
    for d in subs:
        pos = d / "POSCAR"
        if not pos.is_file():
            sys.exit("[ERROR] %s 缺 POSCAR" % d)
        if ref_kpts is not None:
            # patch_kpts_samefile：subs 里包含 undeformed 自己，而基准 KPOINTS
            #   就生成在那儿——拷到自己身上会抛 SameFileError，跳过即可。
            _dst = d / "KPOINTS"
            _same = False
            if _dst.is_file():
                try:
                    _same = os.path.samefile(str(ref_kpts), str(_dst))
                except OSError:
                    _same = (_dst.resolve() == Path(ref_kpts).resolve())
            if not _same:
                # patch_stale_grid：网格变了 -> 该子目录的 OUTCAR 是旧网格的产物。
                # 归档旧产物（改名，不删除），让 ck_deform 判不通过而进入 fan_todo。
                if _dst.is_file():
                    try:
                        _old = _dst.read_text(encoding="utf-8", errors="ignore")
                        _new = Path(ref_kpts).read_text(encoding="utf-8", errors="ignore")
                    except OSError:
                        _old, _new = None, None
                    if _old is not None and _old != _new:
                        _tag = _grid_tag(_old)
                        _moved = []
                        for _nm in STALE_OUTPUTS:
                            _f = d / _nm
                            if _f.is_file():
                                _f.rename(d / (_nm + ".stale-grid-" + _tag))
                                _moved.append(_nm)
                        print("[..] %s: KPOINTS %s -> %s，归档旧产物 %d 个（*.stale-grid-%s）"
                              % (d.name, _tag, _grid_tag(_new), len(_moved), _tag))
                shutil.copy2(str(ref_kpts), str(_dst))
        else:
            kc.vaspkit_kpoints(d, KSCHEME, KSPACING, VASPKIT_EXE, dim, vac_axis)
        kc.vaspkit_potcar(d, VASPKIT_EXE)
        if encut is None:
            encut = MANUAL_ENCUT or kc.encut_from_potcar(d / "POTCAR", ENCUT_FACTOR)
        _sub = {"SYSTEM": "%s deform %s" % (cwd.name, d.name),
                "ENCUT": encut}
        _sub.update(_subs)
        kc.render_tpl(tpl, _sub, d / "INCAR")
        kc.inherit_scf_tags(d / "INCAR", cwd, with_u=True, label=d.name)
        # 并行参数按宿主机自适应（GPU 版强制 NCORE=1/KPAR=1，CPU 保持模板默认）
        kc.apply_parallel_tags(d / "INCAR")
        # step.conf 的 [incar]/[incar.final]/[incar.delete] 覆盖。
        # 此前本技能所有 gen 脚本只调 read_submit()，这三节**写了没人读**，
        # 用户在 step.conf 里改 INCAR 会被静默忽略。2026-09-15 在 step2.3_hse 上
        # 暴露并修复；这里改用 stepconf 里的唯一实现，避免各脚本各写一套。
        _ic_log = []
        stepconf.apply_incar_file(d / "INCAR", log=_ic_log)
        for _m in _ic_log:
            print("[..] %s" % _m)

        sub = d / "submit.sh"
        sub.write_text(submit_body.replace(
            "{{JOBNAME}}", "%s-ke-dft-cpu-%s-%s" % (cwd.name, STEP_LABEL, d.name)),
            encoding="utf-8", newline="\n")
        stepconf.apply_submit(sub, stepconf.read_submit(stepconf.CONF_NAME, used_incar=True))

    # [PATCH-IONRELAX] 给**所有面内分量**的形变目录生成 ionrelax/ 子目录。
    # 2026-09-16 二次修订：原来只建 xx±/yy± 四个。但二维 ADP 核用的是张量面内 2×2
    # 子块（D_xx/D_yy/**D_xy**），面内剪切分量（xy）若仍是离子固定口径，同一张量里
    # 就混了两种来源 —— 按用户定的规则：**面内只要有一个分量是 clamped 就必须报错**。
    # 所以这里按"主导应变分量是否面内"来选目录（而不是只看 xx/yy 配对），
    # 判据与 step7b 的 _inplane_dirs 完全一致。
    if dim == "2d":
        _ir_dirs = set()
        try:
            import numpy as _np
            _und = _np.array(kc.read_lattice_matrix(out / "undeformed" / "POSCAR"))
            for _d in sorted(p for p in out.glob("deform-*") if p.is_dir()):
                _L = _np.array(kc.read_lattice_matrix(_d / "POSCAR"))
                _F = _np.transpose(_np.dot(_np.linalg.inv(_und), _L))
                _E = (_np.dot(_F.T, _F) - _np.eye(3)) / 2.0
                _i, _j = _np.unravel_index(_np.argmax(_np.abs(_E)), _E.shape)
                if (_i, _j) in ((0, 0), (1, 1), (0, 1), (1, 0)):
                    _ir_dirs.add(_d.name)
        except Exception as _e:                                # noqa: BLE001
            print("[WARN] 按应变量判面内分量失败（%s），退回 xx±/yy± 配对" % _e)
            _pairs, _sm = kc.resolve_strain_pairs(out)
            if _pairs:
                _ir_dirs = {_pairs["xx"][0], _pairs["xx"][1],
                            _pairs["yy"][0], _pairs["yy"][1]}
        if _ir_dirs:
            # patch_ionrelax_idempotent（2026-09-20）：ionrelax 已经算完的目录**不要重建**。
            # _build_ionrelax 会清空 ionrelax/ 下除 POSCAR/KPOINTS/POTCAR 之外的全部文件
            # （含 OUTCAR / LOCPOT / vasprun.xml / CONTCAR），所以 retry 只为"补新增的
            # 面内分量"而重建时，会顺带毁掉已完成的旧分量（实测 CrSe2_hex：补 deform-09
            # 的 xy 弛豫会把 deform-01..04 的 ionrelax 结果全删掉）。
            _built, _skipped = [], []
            for _name in sorted(_ir_dirs):
                _irp = out / _name / "ionrelax"
                # 判据与 autozt 的 ck_deform 一致：**必须**在 OUTCAR 末尾看到 VASP 收尾
                # 横幅 "General timing and accounting informations" 才算完成。
                # 只判"文件存在"会把中途挂掉/被杀的 ionrelax 误当完成而跳过（用户
                # 2026-09-20 指出）——挂掉的 OUTCAR/vasprun 也是存在的。
                _done = False
                _oc = _irp / "OUTCAR"
                if _oc.is_file() and (_irp / "vasprun.xml").is_file():
                    try:
                        _done = ("General timing and accounting informations"
                                 in _oc.read_text(errors="ignore")[-400000:])
                    except OSError:
                        _done = False
                # patch_ionrelax_grid（2026-09-20）：**网格也必须一致**。
                #   只判 VASP 横幅会在「升级网格后重跑」时静默保留旧网格的 ionrelax ——
                #   实测 CrSe2_hex：deform-01..04 的 ionrelax 是 15x15x1，升到 46x46x3 后
                #   主单点被刷新、ionrelax 却被跳过，形变势 E1 与其余量混口径。
                #   判据：把 ionrelax/KPOINTS 第 4 行与基准网格（undeformed/KPOINTS 第 4 行）
                #   逐轴比；不一致 -> 视为需重建（旧结果本就该作废）。
                if _done:
                    try:
                        # 修正（2026-09-21）：此处原写 und，但 main() 里没有这个局部名
                        #   （und 只存在于 _reference_kpoints / run_amset_deform_create）。
                        #   只要有 ionrelax 已跑完，这条 NameError 会让 gen 在归档之后、
                        #   建 ionrelax/ 之前整个崩掉，retry 报 exit=1。
                        _ref_ln = ((out / "undeformed" / "KPOINTS").read_text()
                                   .splitlines()[3].split()[:3])
                        _ir_ln = (_irp / "KPOINTS").read_text().splitlines()[3].split()[:3]
                        if [int(x) for x in _ir_ln] != [int(x) for x in _ref_ln]:
                            print("[IONRELAX] %s：ionrelax 网格 %s 与基准 %s 不一致，重建。"
                                  % (_name, "x".join(_ir_ln), "x".join(_ref_ln)))
                            _done = False
                    except (OSError, IndexError, ValueError):
                        _done = False   # 读不出就保守重建，不赌
                if _done:
                    _skipped.append(_name)
                    continue
                _build_ionrelax(out / _name, encut, _subs, submit_body)
                _built.append(_name)
            print("[IONRELAX] 生成 %d 个 ionrelax/（%s）；已存在结果、跳过重建的 %d 个（%s）"
                  % (len(_built), ", ".join(_built) or "无",
                     len(_skipped), ", ".join(_skipped) or "无"))
        else:
            print("[WARN] 认不出面内分量，跳过 ionrelax/（形变势退回刚性单点口径）")

    print("[DONE] %s：%d 个形变子目录输入就绪，tf 会各自提交" % (OUTDIR_NAME, len(subs)))


if __name__ == "__main__":
    main()
