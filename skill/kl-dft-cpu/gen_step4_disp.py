#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step4_disp.py —— 超胞 + 位移生成，扇出（step4_disp）。

从 step1 弛豫结构接力，扩超胞并生成位移超胞，每个位移一个 disp-NNNNN 子目录单点取力。
两种方法（step.conf 的 METHOD）：
  alm     : 随机位移（默认）。位移数 N = max(10, ceil(Σnfree/(3·N_sc)) × OVERSAMPLE)，
            nfree 由 ALM suggest 给出；位移用 hiPhive MC-rattle 生成（高斯幅度 +
            最近邻 d_min 保护 + rattle_std 标定），不是固定模长的 phono3py --rd。
  findiff : phono3py 对称有限位移，位移数由空间群对称约化决定。三阶全对称集在
            大超胞/低对称下轻易上万帧 —— 本步先在内存里数帧数，超 MAX_DISP 直接报错，
            绝不落盘。

关键闸门：MAX_DISP（默认 500）。生成前先算帧数，超限 → sys.exit，不写任何 POSCAR-*/disp-*。
把 DIM/SUPERCELL/MESH/METHOD 写进 kl_params.txt，供 step5/step6 严格继承。
产出目录：step4_disp/{phono3py_disp.yaml, SPOSCAR, disp-00001, ...}
"""
import glob
import json
import logging
import os
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kl_common as kc
from dim_common import require_dim  # noqa: E402
import stepconf
import incar_bench  # noqa: E402 取力 INCAR 基准对照（隔离 side-car，见 incar_bench.py）

OUTDIR = "step4_disp"
STEP   = "step4_disp"
PREV   = ["step1_std_opt"]

SPEC = {
    "FUNC":         ("pbesol", "str"),
    "KSPACING":     ("0.04",   "str"),   # 超胞取力，比原胞略稀
    "KSCHEME":      ("2",      "str"),
    "ENCUT":        (None,     "int"),
    "ENCUT_FACTOR": (1.5,      "float"),
    # 取力精度（P1-3）：力直接进 fc2/fc3，phonopy 官方超胞取力示例用
    #   PREC=Accurate / EDIFF=1E-8 / LREAL=.FALSE.；LREAL=Auto 的实空间投影带格点噪声，
    #   是 fc 假软模的常见来源。超胞几百原子时可用 FORCE_LREAL=Auto 换速度，
    #   但要先在 3 帧上与 .FALSE. 比力（最大偏差 < 1 meV/Å）并显式设 ROPT。
    "FORCE_PREC":   ("Accurate", "str"),
    "FORCE_EDIFF":  ("1E-8",     "str"),
    "FORCE_LREAL":  (".FALSE.",  "str"),
    "VASPKIT_EXE":  ("vaspkit", "str"),
    "METHOD":       ("alm",    "str"),   # alm | findiff
    "SUPERCELL":    (None,     "words"), # 显式 "3 3 3"；空=按 MIN_SC_LEN 自动
    "MIN_SC_LEN":   (12.0,     "float"),
    "MAX_MULTIPLE": (6,        "int"),
    # q 网格（写进 kl_params 供 step6）。auto = 2D 按倒空间长度估：
    #   N_i = max(MESH_MIN, ceil(Q_LEN_2D/|a_i|))，真空轴=1（P1-1）。
    #   写死 15 15 15 对 a≈3 Å 的 2D 材料只覆盖 ~47 Å 倒空间，κ 远未收敛。
    "KAPPA_MESH":   ("auto",   "str"),   # auto | auto3d | "N N N"
    "Q_LEN_2D":     (275.0,    "float"), # auto 时目标倒空间长度(Å)；文献 2D 用 250~300
    "MESH_MIN":     (20,       "int"),   # auto 时每个方向的下限（BZ 采样下限）
    # P1-4：2D 残余面内应力门禁（层内口径 σ_VASP×h⊥/d，kbar）。0 = 关掉门禁。
    #   Huang 条件是零应力条件；有外应力时 ZA 出现线性项、压制二次项，κ 不可信。
    "STRESS_2D_THR": (0.5,     "float"),
    "MAX_DISP":     (500,      "int"),   # ★ 位移帧数硬闸：超过就报错停步
    "FC3_CUTOFF_PAIR": (None,  "float"), # findiff 三阶对距离上限(Å)；空=全对称集
    "FD_DISTANCE":  (0.03,     "float"), # findiff 位移幅度(Å)
    "OVERSAMPLE":   (3,        "int"),   # alm 过采样系数（超定 3 倍）
    "ALM_CUT2":     (None,     "float"), # alm 二阶截断(Å)；空=不截断
    "ALM_CUT3":     (6.0,      "float"), # alm 三阶截断(Å)
    "DISP_RMS":     (0.03,     "float"), # MC-rattle 目标位移模长 RMS(Å)
    "MC_DMIN_SCALE": (0.75,    "float"), # MC d_min = 最近邻 × 此系数
    "MC_NITER":     (10,       "int"),   # MC 迭代数
    "MC_SEED":      (2025,     "int"),
    # INCAR_BENCH（2026-09-18，wangchao 要求）：取力 INCAR 基准对照。
    #   on 时除正常生成生产 disp-* 外，额外在 step4_disp/bench/<variant>/disp-XXXXX
    #   生成少数几帧（同一批位移，直接 copy 生产 POSCAR），比较 LREAL×ADDGRID 与
    #   ENCUT 变体的力一致性，定下生产取力 INCAR。默认 off，不改变现有材料行为。
    #   隔离：bench/ 不在生产 fanout（"disp-*" 非递归）与 S5 帧校验范围内。
    #   分析：python incar_bench.py analyze（或 submit/analyze 见 incar_bench.py）。
    "INCAR_BENCH":  ("off",    "str"),   # off | on
    "INCAR_BENCH_FRAMES": (2,  "int"),   # 每变体取生产前几帧（默认 2）
}
# 结构体检第五条（B）：SYMMETRY_AUDIT/SYMMETRIZE 等键（见 kl_common.SYMMETRY_SPEC）
SPEC.update(kc.SYMMETRY_SPEC)


# ==========================================================================
# phono3py 对象 / 帧数统计
# ==========================================================================
def make_ph3(out, reps):
    """用 phono3py Python API 建对象（等价 CLI 的 --dim=... --pa auto -c POSCAR）。"""
    try:
        from phono3py import Phono3py
        from phonopy.interface.calculator import read_crystal_structure
    except Exception as e:
        sys.exit("[ERROR] 无法 import phono3py/phonopy（%s）—— step4 需要在装了 "
                 "phono3py 的 conda 环境里跑 gen" % e)
    cell, _ = read_crystal_structure(str(out / "POSCAR"), interface_mode="vasp")
    _reps = np.array(reps, dtype=int)
    # 3 个数=对角扩胞；9 个数=3×3 矩阵（行主序，与 phonopy/phono3py --dim 同义）
    scm = _reps.reshape(3, 3) if _reps.size == 9 else np.diag(_reps)
    return Phono3py(cell, supercell_matrix=scm, primitive_matrix="auto")


def count_findiff_frames(ph3, distance, cutoff_pair):
    """只数帧数、不建超胞：一阶位移 + 被 included 的二阶对位移。"""
    ph3.generate_displacements(distance=distance, cutoff_pair_distance=cutoff_pair)
    ds = ph3.dataset
    n = 0
    for fa in ds.get("first_atoms", []):
        n += 1
        for sa in fa.get("second_atoms", []):
            if sa.get("included", True):
                n += 1
    return n


def gate(n, conf, method, extra=""):
    """位移帧数硬闸。超过 MAX_DISP 直接停步，绝不落盘。"""
    cap = int(conf["MAX_DISP"])
    print("[..] 计划位移帧数 = %d（上限 MAX_DISP=%d）%s" % (n, cap, extra))
    if n <= cap:
        return
    if method == "findiff":
        hint = ("findiff 三阶全对称集在这个超胞下爆了。请二选一：\n"
                "  ① 换随机位移：tf -tt kl-dft-cpu -p <材料> -j step4_disp conf --set params.METHOD=alm\n"
                "  ② 收三阶对距离：... conf --set params.FC3_CUTOFF_PAIR=4.0（或 5.0）\n"
                "  ③ 缩超胞：... conf --set params.MIN_SC_LEN=12.0 / params.SUPERCELL=\"2 2 2\"")
    else:
        hint = ("ALM 自由参数太多。请二选一：\n"
                "  ① 收三阶截断：... conf --set params.ALM_CUT3=4.0\n"
                "  ② 缩超胞：... conf --set params.MIN_SC_LEN=12.0 / params.SUPERCELL=\"2 2 2\"\n"
                "  ③ 降过采样：... conf --set params.OVERSAMPLE=3（不建议低于 3）")
    sys.exit("[ERROR] 位移帧数 %d > MAX_DISP=%d，step4 已停止（未写任何 disp-*）。\n%s"
             % (n, cap, hint))


# ==========================================================================
# alm 分支：ALM 定帧数 + hiPhive MC-rattle 生成位移
# ==========================================================================
def plan_alm(out, ph3, conf):
    """ALM suggest → nfree → N = max(10, ceil(Σnfree/(3·N_sc))×OVERSAMPLE)。返回 (N, nfree, atoms)。"""
    import lattice_kappa as lk
    from ase import Atoms
    from phonopy.interface.vasp import write_vasp

    sc = ph3.supercell
    write_vasp(str(out / "SPOSCAR"), sc, direct=True)
    atoms = Atoms(numbers=sc.numbers, positions=sc.positions,
                  cell=sc.cell, pbc=True)
    n_sc = len(atoms)

    orders, cuts = [2, 3], [conf["ALM_CUT2"], conf["ALM_CUT3"]]
    # kl10: 优先用 ALM Python API 取 nfree —— ALM 的命令行可执行文件要单独
    # cmake 编译，pip 装的 Python 包并不产出它（PATH 里的 alm 只是坏入口脚本）。
    # API 拿不到（没装 alm 包等）才回落到"写 alm.in + 跑 alm + 解析日志"的老路。
    lk.write_alm_suggest_input(atoms, orders, cuts, out / "alm.in")  # 留档，便于复查
    try:
        nfree = lk.alm_nfree_via_api(atoms, orders, cuts)
    except Exception as _e:
        print("[..] ALM Python API 不可用（%s），回落到命令行 alm" % _e)
        try:
            lk.run_alm(out)
            nfree = lk.parse_alm_nfree(out / "alm.log", orders)
        except Exception as _e2:
            sys.exit("[ERROR] ALM 不可用（Python API: %s；命令行: %s）。\n"
                     "        ALM 是确定随机位移帧数的必需依赖（算各阶自由参数 nfree），\n"
                     "        未安装就禁止推进 step4。请先安装 ALM（cmake 编译 + Python 绑定\n"
                     "        from alm import ALM），再重跑 S4_disp。"
                     % (_e, _e2))
    n = int(lk.estimate_n_struct(nfree, n_sc, int(conf["OVERSAMPLE"])))
    detail = "  [自由参数 %s / DOF=3×%d=%d × OVERSAMPLE=%d]" % (
        nfree, n_sc, 3 * n_sc, int(conf["OVERSAMPLE"]))
    return n, nfree, atoms, detail


def make_mc_rattle_dataset(ph3, atoms, n, conf):
    """hiPhive MC-rattle（标定 rattle_std + d_min 保护）→ 塞进 phono3py 的 type-2 dataset。"""
    import lattice_kappa as lk
    structs, rattle_std, rms = lk.generate_calibrated_mc_rattle(
        atoms, n, float(conf["DISP_RMS"]), float(conf["MC_DMIN_SCALE"]),
        int(conf["MC_NITER"]), int(conf["MC_SEED"]))
    ref = atoms.get_positions()
    disps = np.array([s.get_positions() - ref for s in structs], dtype=float)
    if disps.shape != (n, len(atoms), 3):
        sys.exit("[ERROR] MC-rattle 位移形状异常 %s" % (disps.shape,))
    if float(np.linalg.norm(disps, axis=2).std()) < 1e-4:
        sys.exit("[ERROR] 位移模长 std≈0 —— 不是高斯分布，MC-rattle 没生效")
    try:
        ph3.dataset = {"displacements": disps}
    except Exception:
        ph3.displacements = disps
    return rattle_std, rms


# ==========================================================================
# 落盘
# ==========================================================================
def write_supercells(out, ph3):
    """写 POSCAR-XXXXX + phono3py_disp.yaml。返回写出的编号列表。"""
    from phonopy.interface.vasp import write_vasp
    nums = []
    for i, s in enumerate(ph3.supercells_with_displacements, 1):
        if s is None:
            continue
        write_vasp(str(out / ("POSCAR-%05d" % i)), s, direct=True)
        nums.append("%05d" % i)
    ph3.save(filename=str(out / "phono3py_disp.yaml"))
    if not (out / "phono3py_disp.yaml").is_file():
        sys.exit("[ERROR] 未产出 phono3py_disp.yaml")
    return nums


def _poscar_cell(path):
    """读 POSCAR 的 3x3 晶格矩阵（纯文本解析，不引依赖）；失败返回 None。"""
    try:
        ln = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
        scale = float(ln[1].split()[0])
        return [[float(x) * scale for x in ln[i].split()[:3]] for i in (2, 3, 4)]
    except Exception:                                   # noqa: BLE001
        return None


def _disp_yaml_cell(path):
    """读 phono3py_disp.yaml 里的 unit_cell.lattice；失败返回 None。"""
    try:
        import yaml
        d = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        return [[float(x) for x in row] for row in d["unit_cell"]["lattice"]]
    except Exception:                                   # noqa: BLE001
        return None


def check_existing_matches_input(out):
    """结构指纹校验（2026-09-17，wangchao 要求）：已有的位移/力数据必须属于【当前输入结构】。

    为什么必须拦：retry 对 fanout 步骤**只补缺失的子目录**，已有 disp-* 一律不动；而
    build_displacements 见到 phono3py_disp.yaml + POSCAR-* 就"幂等跳过"。于是 S1 换了
    结构之后重投 S4，旧位移、旧帧、旧力会被整套沿用 —— S5 拿旧结构的力配旧位移拟合，
    **内部自洽但与当前结构无关**，而且一声不响。实测 2026-09-17 这一批 6 个材料全部
    中招（残留 1902 个 disp-*，含半成品 vasprun.xml）。
    """
    y = out / "phono3py_disp.yaml"
    if not y.is_file():
        return
    old, new = _disp_yaml_cell(y), _poscar_cell(out / "POSCAR")
    if old is None or new is None:
        print("[WARN] 结构指纹校验跳过（读不到 %s 或 %s）"
              % ("phono3py_disp.yaml" if old is None else "POSCAR",
                 "POSCAR" if new is None else "phono3py_disp.yaml"))
        return
    dmax = max(abs(old[i][j] - new[i][j]) for i in range(3) for j in range(3))
    if dmax > 1e-3:
        n_disp = len(glob.glob(str(out / "disp-*")))
        sys.exit(
            "[ERROR] 已有位移/力数据属于【旧结构】—— phono3py_disp.yaml 的单胞与当前输入 "
            "差异最大 %.4f Å（>1e-3 即判定不一致）。\n"
            "        当前有 %d 个 disp-* 目录。retry 对 fanout 只补缺失帧，这些旧帧会被"
            "整套沿用，\n"
            "        于是 S5 会把旧结构的力配上（同样旧的）位移静默混拟合。\n"
            "        处置：确认已归档后清空本步产物，再 retry：\n"
            "          rm -rf %s/disp-* %s/POSCAR-* %s/phono3py_disp.yaml %s/SPOSCAR\n"
            "        或直接 autozt -tt kl-dft-cpu -p <材料> -j S4_disp rerun（推倒重来）。"
            % (dmax, n_disp, out, out, out, out))


def build_displacements(out, reps, method, conf):
    """幂等：已有位移就跳过；否则先算帧数过闸，再落盘。"""
    check_existing_matches_input(out)
    if (out / "phono3py_disp.yaml").is_file() and glob.glob(str(out / "POSCAR-*")):
        print("[..] 已有位移超胞，跳过生成（幂等）")
        return
    ph3 = make_ph3(out, reps)
    info = {"method": method, "n_atoms_sc": len(ph3.supercell)}

    if method == "findiff":
        n = count_findiff_frames(ph3, float(conf["FD_DISTANCE"]),
                                 conf["FC3_CUTOFF_PAIR"])
        gate(n, conf, method, "  [对称有限位移，cutoff_pair=%s]" % conf["FC3_CUTOFF_PAIR"])
    else:
        n, nfree, atoms, detail = plan_alm(out, ph3, conf)
        gate(n, conf, method, detail)
        std, rms = make_mc_rattle_dataset(ph3, atoms, n, conf)
        info.update(nfree=nfree, oversample=int(conf["OVERSAMPLE"]),
                    rattle_std=round(std, 6), disp_rms=round(rms, 5))

    nums = write_supercells(out, ph3)
    if not nums:
        sys.exit("[ERROR] 没有产出 POSCAR-* 位移超胞")
    info["n_disp"] = len(nums)
    (out / "disp_plan.json").write_text(
        json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    print("[OK] 位移生成完毕：%d 帧（%s）" % (len(nums), method))


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)

    prev = kc.find_prev_dir(cwd, PREV)
    if prev is None:
        sys.exit("[ERROR] 找不到 step1_std_opt 的结构")
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR", "step1_std_opt")
    # 结构体检第五条（B，2026-09-20）：S1 弛豫后的数值微畸变审计 + 可选对称化
    # （默认 SYMMETRY_AUDIT=warn 只告警、SYMMETRY_SYMMETRIZE=off 不改结构）
    kc.symmetry_gate(out / "POSCAR", conf)

    meth = kc.read_method(prev / kc.METHOD_FILE)
    dim = (meth.get("DIM", "").lower() or kc.resolve_dim(out / "POSCAR")[0])
    _, vac_axis = kc.resolve_dim(out / "POSCAR", dim)
    require_dim(dim, ('2d', '3d'), "step4_disp",
                why="晶格热导需要声子群速度和布里渊区积分，孤立分子只有分立振动模式")
    func = conf["FUNC"] if conf["FUNC"] not in (None, "", "auto") \
        else meth.get("FUNC", "pbe-d3").lower()
    method = str(conf["METHOD"]).lower()
    if method not in ("findiff", "alm"):
        sys.exit("[ERROR] METHOD 只允许 findiff / alm")

    # 三阶截断半径（alm 的 ALM_CUT3 / findiff 的 FC3_CUTOFF_PAIR）用于超胞内切球判据：
    # 截断半径必须 ≤ 超胞安全截断（0.5×内切球直径 − margin），否则周期镜像污染力常数。
    cut3 = conf["ALM_CUT3"] if method == "alm" else conf["FC3_CUTOFF_PAIR"]
    cut3_label = "ALM_CUT3" if method == "alm" else "FC3_CUTOFF_PAIR"
    if dim == "2d":
        # P1-5：真空隙 ≤ 层厚 → 最近镜像跨真空（Born-Huang 的 r_ij 错）；截断 ≥ 真空隙
        #   → 跨真空三体簇。两件都是硬错误，放在生成位移之前拦。
        _geo = kc.check_2d_vacuum(out / "POSCAR", vac_axis, cut3=cut3)
        # P1-4：残余面内应力门禁（Huang 条件）。应力从 S1 弛豫的 OUTCAR 末态读，
        #   换算成层内口径 σ×h⊥/d 再比阈值。
        _s1 = sorted(cwd.glob("step1*/OUTCAR"))
        if _s1:
            kc.check_2d_stress(_s1[-1], _geo["h_perp_A"], _geo["thickness_d_A"],
                               conf["STRESS_2D_THR"])
        else:
            print("[WARN] 找不到 S1 的 OUTCAR（%s/step1*/OUTCAR），2D 残余应力未核验"
                  % cwd)
    if conf["SUPERCELL"]:
        reps = [int(x) for x in conf["SUPERCELL"]]
        if dim == "2d":
            reps[vac_axis if vac_axis is not None else 2] = 1
        # 显式超胞不自动扩胞，只事后 WARN（对齐 validate_user_supercell 只 warning 语义）
        kc.warn_cutoff_vs_supercell(out / "POSCAR", reps, cut3, label=cut3_label,
                                    dim=dim, vac_axis=vac_axis)
    else:
        reps = kc.supercell_matrix(out / "POSCAR", dim, conf["MIN_SC_LEN"],
                                   conf["MAX_MULTIPLE"], vac_axis if vac_axis is not None else 2,
                                   cutoff=cut3)
    vac_ax = vac_axis if vac_axis is not None else 2
    mesh, mesh_note = kc.auto_mesh(conf["KAPPA_MESH"], dim, vac_ax, out / "POSCAR",
                                   conf["Q_LEN_2D"], conf["MESH_MIN"])
    print("[..] 维度=%s 方法=%s 超胞=%s mesh=%s %s"
          % (dim.upper(), method, kc.dim_str(reps), mesh, mesh_note))
    kc.write_kl_params(out / kc.KL_PARAMS, DIM=dim.upper(), SUPERCELL=kc.dim_str(reps),
                       MESH=mesh, METHOD=method, FUNC=func)

    build_displacements(out, reps, method, conf)

    poscars = sorted(glob.glob(str(out / "POSCAR-*")))
    if not poscars:
        sys.exit("[ERROR] 没有产出 POSCAR-* 位移超胞")
    # 二道闸：目录里本来就躺着一堆旧帧（上一版无闸门跑出来的）也要拦住，别扇出。
    if len(poscars) > int(conf["MAX_DISP"]):
        sys.exit("[ERROR] %s 下已有 %d 个 POSCAR-*，超过 MAX_DISP=%d —— 拒绝扇出。\n"
                 "        这多半是旧版本（无闸门）留下的，先清干净再重跑：\n"
                 "        rm -rf %s && tf -tt kl-dft-cpu -p <材料> -j step4_disp rerun"
                 % (OUTDIR, len(poscars), int(conf["MAX_DISP"]), OUTDIR))
    print("[..] 位移超胞 %d 个，逐个补输入" % len(poscars))

    # kleq：平衡帧 disp-00000 —— 未位移的完美超胞单点，INCAR/K 点/泛函与位移帧完全一致。
    #   ① 拟合时逐帧扣它的残余力（F_disp - F_eq），去掉未完全弛豫 / egg-box 的净力；
    #   ② 它是"缺帧容错"的锚：随机位移可以少几帧，平衡帧不能缺。
    #   phono3py 的位移帧从 POSCAR-00001 起编号，00000 正好空出来给它。
    sposcar = out / "SPOSCAR"
    if not sposcar.is_file():
        from phonopy.interface.vasp import write_vasp as _write_vasp
        _write_vasp(str(sposcar), make_ph3(out, reps).supercell, direct=True)
    frames = [("00000", str(sposcar))]
    frames += [(os.path.basename(p).split("-", 1)[1], p) for p in poscars]

    here = Path(__file__).resolve().parent
    incar_tpl = kc.resolve_submit(here, dim, "incar_force")
    submit_tpl = kc.resolve_submit(here, dim, "submit_std")
    # P2-3：偶极修正必须与 S1(弛豫)/S2(静态) 一致 —— 直接继承它们的 INCAR，
    #   不复算（三处各自判一遍，判据一漂就没意义了）。取力与弛豫处在不同静电
    #   边界条件下时，位移帧的力根本不对应同一势能面。
    _dip_lines, _dip_src = kc.dipole_line_from_incar(cwd / "step2_static" / "INCAR",
                                                     cwd / "step1_std_opt" / "INCAR")
    print("[..] 偶极修正：%s" % (_dip_src or "上游 INCAR 没找到，未施加"))
    encut = None
    for num, pos in frames:
        d = out / ("disp-%s" % num)
        d.mkdir(exist_ok=True)
        if (d / "INCAR").is_file() and (d / "POSCAR").is_file():
            continue
        shutil.copyfile(pos, d / "POSCAR")
        kc.vaspkit_kpoints(d, conf["KSCHEME"], conf["KSPACING"], conf["VASPKIT_EXE"], dim, vac_axis)
        kc.vaspkit_potcar(d, conf["VASPKIT_EXE"])
        if encut is None:
            encut = conf["ENCUT"] or kc.encut_from_potcar(d / "POTCAR", conf["ENCUT_FACTOR"])
        # patch_kl_vdw：step2/step3 都渲染了 VDW_LINE，唯独取力这步漏了 ——
        # FUNC=pbe-d3 时会变成"带 D3 弛豫、无 D3 取力"，平衡位置残留净力，
        # fc2 里混进伪软模，fc3 更不可用。
        kc.render_tpl(incar_tpl, {"SYSTEM": "%s force %s" % (cwd.name, num),
                                  "ENCUT": encut,
                                  "GGA": kc.GGA_MAP.get(func, "PS"),
                                  "FORCE_PREC": conf["FORCE_PREC"],
                                  "FORCE_EDIFF": conf["FORCE_EDIFF"],
                                  "FORCE_LREAL": conf["FORCE_LREAL"],
                                  "DIPOLE_LINE": "\n".join(_dip_lines),
                                  "VDW_LINE": ("IVDW = %s" % kc.VDW_MAP[func])
                                  if kc.VDW_MAP.get(func) else "# no vdW"},
                      d / "INCAR")
        # P1-3：项目级 incar_force_*.tpl 会遮蔽技能模板（find_asset 优先项目根），
        #   老副本没有 {{FORCE_PREC}} 等占位符 → 渲染时被静默跳过。渲染后强制落一遍。
        kc.enforce_incar_tags(d / "INCAR",
                              {"PREC": conf["FORCE_PREC"],
                               "EDIFF": conf["FORCE_EDIFF"],
                               "LREAL": conf["FORCE_LREAL"],
                               "ADDGRID": ".TRUE."},
                              label="step4 取力 %s：" % num)
        if _dip_lines and any("LDIPOL = .TRUE." in x for x in _dip_lines):
            kc.enforce_incar_tags(d / "INCAR", {"LDIPOL": ".TRUE.", "IDIPOL": "3",
                                                "ISYM": "0"},
                                  label="step4 偶极 %s：" % num)
        kc.write_submit(submit_tpl, d / "submit.sh",
                        {"JOBNAME": "%s-kl-dft-cpu-S4-%s" % (cwd.name, num)})
        stepconf.apply_submit(d / "submit.sh", conf.submit)

    # ---- INCAR_BENCH：取力 INCAR 基准对照（默认 off；隔离 side-car）----
    #   bench/ 不在生产 fanout（"disp-*" 非递归）与 S5 帧校验范围内，故它不会被
    #   autozt 当生产帧提交/统计；用独立命令 python incar_bench.py submit/analyze。
    if incar_bench.enabled(conf["INCAR_BENCH"]):
        try:
            _bp = incar_bench.generate(out, conf=conf, dim=dim,
                                       n_frames=int(conf["INCAR_BENCH_FRAMES"] or 2))
        except Exception as _e:                         # noqa: BLE001
            _bp = None
            print("[WARN] INCAR_BENCH 生成失败（%s）—— 生产位移帧不受影响" % _e)
        if _bp is not None:
            print("[..] INCAR_BENCH：bench 帧不随生产 fanout 提交；登录节点跑 "
                  "python incar_bench.py submit 提交、python incar_bench.py analyze "
                  "出对照表（判据 |dF|_rms/|F|_rms < 1%%；绝对 "
                  "max|dF| < 0.5 meV/A 可用 --threshold 选）。")

    print("[DONE] %s：%d 个位移子目录 + 平衡帧 disp-00000 就绪，"
          "tf 各自提交（fanout disp-*）" % (OUTDIR, len(poscars)))


if __name__ == "__main__":
    main()
