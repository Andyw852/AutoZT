#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step4_kappa.py —— BTE 晶格热导率（step4_kappa），提交计算节点。

求解器（step.conf 的 SOLVER）：
  phono3py（默认）：从 step3_fc 取 fc2/fc3.hdf5（+BORN），按 klmlff_params 的 MESH 组
      phono3py-load --br/--lbte 命令，渲染提交模板。
  shengbte       ：从 step3_fc/shengbte/ 取 FORCE_CONSTANTS_2ND/3RD（step3 拟合后统一导出，
      见 fc_fit_driver.py），写 ShengBTE CONTROL，渲染提交模板跑 ShengBTE 三声子 BTE。

作业跑完就地抽 κ 张量写 kappa_summary.json（marker: KAPPA_DONE）。

MESH_SCAN：分号分隔的多套 q 网格，例如 "16 16 16; 20 20 20; 24 24 24"。
DFT 路线不敢做的收敛测试，在这条链上几乎白送——**κ 对 q 网格的收敛性是这类结果最
常被审稿人问的一条**，顺手扫出来比事后补便宜得多。扫描时每套网格的 κ 都进 summary。
（MESH_SCAN 仅 phono3py 支持；shengbte 单套网格。）

MESH_SCAN：分号分隔的多套 q 网格，例如 "16 16 16; 20 20 20; 24 24 24"。
DFT 路线不敢做的收敛测试，在这条链上几乎白送——**κ 对 q 网格的收敛性是这类结果最
常被审稿人问的一条**，顺手扫出来比事后补便宜得多。扫描时每套网格的 κ 都进 summary。
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import klmlff_common as kc
import stepconf

OUTDIR = "step4_kappa"
STEP = "step4_kappa"
SRC = "step3_fc"

SPEC = {
    # ---- 全局 ----
    "MACE_MODEL": ("mace-mp:medium", "str"),
    "MACE_MODEL_DIR": ("", "str"),
    "DEVICE": ("auto", "str"),
    "DTYPE": ("float64", "str"),
    "CONDA_SH": (kc.DEFAULT_CONDA_SH, "str"),
    "CONDA_ENV": (kc.DEFAULT_CONDA_ENV, "str"),
    # ---- 本步 ----
    "MESH_OVERRIDE": (None, "str"),     # 空=用 step2 写进 klmlff_params 的 MESH
    "MESH_SCAN": ("", "str"),           # "16 16 16; 20 20 20"，空=只跑一套
    "T_MIN": (100, "int"),
    "T_MAX": (800, "int"),
    "T_STEP": (100, "int"),
    "ISOTOPE": (True, "bool"),
    "BTE": ("rta", "str"),              # rta（--br）| lbte（--lbte，直接解，贵得多）
    "EXTRA_ARGS": ("", "str"),          # 原样附加，如 "--boundary-mfp 1e6" / "--write-gamma"
    # ---- 求解器（phono3py 默认 | shengbte）----
    "SOLVER": ("phono3py", "str"),      # phono3py（--br/--lbte，默认）| shengbte（三声子 BTE）
    # shengbte 专属（phono3py 忽略）
    "SCALEBROAD": (1.0, "float"),       # shengbte 高斯展宽系数（ShengBTE 默认 1.0；0.1 会漏过程）
    "SHENGBTE_EXE": ("ShengBTE", "str"),# shengbte 可执行（按集群填绝对路径，见 kl-dft README）
    # ---- ShengBTE 的 MPI/OMP 布局（见 setting/<集群>/templates/submit_shengbte_klm.tpl）----
    # ShengBTE 靠 MPI 按 q 点并行：--ntasks-per-node=1 -> $SLURM_NTASKS=1 ->
    # mpirun -n 1 = 串行，会长时间停在 "about to obtain the spectrum" 零产物。
    # 正确布局 = 一个 MPI rank 占一个 NUMA **节点**，且 **rank 数必须 ≤ NUMA 节点数**：
    #     NTASKS = TOTAL_CORES / CORES_PER_NUMA ，CPUS_PER_TASK = CORES_PER_NUMA
    # jzzn 计算节点实测（lscpu）：192 核 = 2 socket x 96 = **8 个 NUMA 节点 x 24 核**
    #   （早先误记为 24 域 x 8 核）。rank 数 > 节点数时 OpenMPI 会循环复用节点、把多个
    #   rank 绑到同一组核上（实测 24 rank x 8 线程时整作业只用 64/192 核）。故：
    #   192 核 -> 8 rank x 24 线程；96 核 -> 4 rank x 24 线程。
    "SHENGBTE_TOTAL_CORES":    (96,     "int"),
    "SHENGBTE_CORES_PER_NUMA": (24,     "int"),
    "SHENGBTE_NTASKS":         ("auto", "str"),
    # phono3py 网格点并行：>=2 时先 --wgp 写出 ir_grid_points.yaml，把其中**不可约格点
    #   编号**均分成 N 份，每份用 --gp="编号列表" + --write-gamma 并行算散射率，最后用
    #   --read-gamma（不带 --gp）汇总出 κ。单机多核白拿的加速，BTE 网格加密后尤其值得。
    #   ★ --gp 接的是「网格点编号列表」，不是「第 i/n 份」：写成 --gp i N 只会算编号 i、N
    #     两个格点（编号越界还会直接报错），κ/γ 静默不完整。=0/1 关闭，走单作业。
    "GP_SPLIT": (0, "int"),
    # MFP 累积 κ 诊断：加 --mfp 让 phono3py 写 kappa-mfp.hdf5，extract 把「累积 κ vs 声子
    #   自由程」与 50%/90% 累积处的 MFP 收进 kappa_summary.json —— 审稿人常问的热输运尺度。
    "MFP_CUMULATIVE": (False, "bool"),
    # 2D κ 厚度归一化：phono3py 用含真空的原胞体积做分母，2D 面内 κ 被 Lz 稀释；
    #   需乘 Lz/d。取法：vdw=原子z跨度+两侧vdW半径 | cell=用Lz(即不归一) | 数值=固定Å
    "KAPPA_2D_THICKNESS": ("vdw",  "str"),
    # 2D NAC：auto=2D默认不用3D-NAC(LO-TO在2D应趋零)/3D随BORN；on/off 强制
    "KAPPA_NAC":          ("auto", "str"),
}

# patch_vdw_shared：单一真源在 skill/_common/vdw_radii.py，ke-dft-cpu/kl-dft-cpu 读同一份。
#   两边算层厚用同一个公式 d = zspan + vdW(top) + vdW(bot)，表必须同源。
try:
    from vdw_radii import VDW_RADII as _VDW        # noqa: F401
except ImportError:
    _VDW = {"H":1.20,"B":1.92,"C":1.70,"N":1.55,"O":1.52,"F":1.47,"Na":2.27,"Mg":1.73,
            "Al":1.84,"Si":2.10,"P":1.80,"S":1.80,"Cl":1.75,"K":2.75,"Ca":2.31,
            "Ti":2.11,"V":2.07,"Cr":2.06,"Mn":2.05,"Fe":2.04,"Co":2.00,"Ni":1.97,
            "Cu":1.96,"Zn":2.01,"Ga":1.87,"Ge":2.11,"As":1.85,"Se":1.90,"Br":1.85,
            "Mo":2.17,"Ru":2.13,"Rh":2.10,"Pd":2.10,"Ag":2.11,"Cd":2.18,"In":1.93,
            "Sn":2.17,"Sb":2.06,"Te":2.06,"I":1.98,"W":2.18,"Pt":2.13,"Au":2.14,
            "Hg":2.23,"Tl":1.96,"Pb":2.02,"Bi":2.07}
    print("[WARN] 没找到 _common/vdw_radii.py，用内置最小表")

# 2D 层厚的唯一真源（kl-dft-cpu S6 与 ke 侧同一份）：h⊥=V/A、以真空隙切口展开、
#   d = span + vdW(top) + vdW(bot)。本步不再自己算，只做字段名映射。
from thickness_2d import slab_geometry


def _poscar_species(poscar):
    """从 POSCAR 读每个原子的元素符号（VASP5）；VASP4 无元素行则返回 None。"""
    import re as _re
    L = Path(poscar).read_text(encoding="utf-8-sig").splitlines()
    line6 = L[5].split()
    if not line6 or _re.fullmatch(r"[+-]?\d+", line6[0]):
        return None
    counts = [int(x) for x in L[6].split()]
    out = []
    for sym, n in zip(line6, counts):
        out += [sym] * n
    return out


def two_d_norm_factor(poscar, vac_axis, mode):
    """2D κ 厚度归一化因子 factor=h⊥/d 与元数据。mode: vdw | cell | 数值(Å)。

    真源是 skill/_common/thickness_2d.slab_geometry（kl-dft-cpu 的 S6 与 ke 侧读同一份），
    本函数只负责从 POSCAR 取晶格+元素表再转回 kappa_summary 的字段名。

    相对旧实现修了三处（与 kl-dft-cpu/gen_step6_kappa.py 的同名函数一致）：
      ① 基准从 Lz=|c| 改成 h⊥=V/A —— phono3py/ShengBTE 的分母是 V 和面内面积 A，
         c 轴倾斜时 |c|≠V/A，因子会算错；
      ② 展开方式改成"最大间隙=真空"切口 —— 旧版直接取 frac 投影 max−min，层跨越
         z=0/1 时 span 被算成接近整个胞高 → d>h⊥ → 因子<1（把 κ 又缩小一遍）；
      ③ 元数据补 h_perp_A / vacuum_gap_A / thickness_d_A。

    ★ 字段名保持 Lz_ang / thickness_d_ang / atomic_zspan_ang 不变：ke 的
      step8.1_boltztrap/gen_step11_boltztrap.py 与 step8.3_output/gen_step13_output.py
      直接读这些键做两条链的元胞/厚度一致性闸门。Lz_ang 仍是 |c|（要跟 ke 的 c 比），
      归一化用的胞高是 h_perp_A，两者分开了。
    """
    from dim_common import read_poscar_cell_frac
    lat, frac = read_poscar_cell_frac(poscar)
    ax = vac_axis if vac_axis is not None else 2
    sp = _poscar_species(poscar)
    if sp is not None and len(sp) != len(frac):
        sp = None
    g = slab_geometry(lat, frac, sp, vac_axis=ax, mode=mode)
    factor = float(g["kappa_2d_norm_factor"])
    meta = {"kappa_2d_norm_factor": g["kappa_2d_norm_factor"],
            "h_perp_A": g["h_perp_A"], "vacuum_gap_A": g["vacuum_gap_A"],
            "thickness_d_A": g["thickness_d_A"],
            # —— 兼容旧字段名（ke 侧在读，勿删）——
            "Lz_ang": g["Lz_A"],
            "thickness_d_ang": g["thickness_d_A"],
            "thickness_convention": g["thickness_convention"],
            "atomic_zspan_ang": g["atomic_span_A"],
            "note": "kappa_2d_normalized = kappa_raw * h_perp/d"
                    "（h_perp=V/A，面内分量才有物理意义）"}
    return factor, meta

# 作业里抽 κ 的小脚本：把所有 kappa-m*.hdf5 都收进 summary（配合 MESH_SCAN）。
# factor!=1（2D）时每套网格额外给 kappa_2d_normalized（原始 κ × Lz/d）。
# 张量审计：读完整 3x3 κ 张量；按 spglib 空间群定晶系，对要求面内各向同性的晶系
# （hex/cubic/trigonal/tetragonal）用 kappa_validation 的敏感判据审 300K 面内块。
# MFP 累积 κ：--mfp 时 phono3py 写 kappa-mfp.hdf5，把累积分数与 50%/90% MFP 收进 summary。
def build_extract(factor, meta):
    import json as _json
    return (
        "python - <<'PY'\n"
        "import glob, json, os\n"
        "import numpy as np, h5py\n"
        "FACTOR=%r\n" % float(factor) +
        "META=json.loads(%r)\n" % _json.dumps(meta, ensure_ascii=False) +
        "out = {\"KAPPA_DONE\": False, \"runs\": []}\n"
        "out.update(META)\n"
        "# ---- 晶系判定 + 面内投影轴（2D 用真空方向，3D 用 None=原始 xy 块）----\n"
        "CRYSTAL, SG, PLANE = None, None, None\n"
        "try:\n"
        "    import spglib\n"
        "    from phonopy.interface.calculator import read_crystal_structure\n"
        "    from kappa_validation import crystal_system_from_spacegroup, audit_kappa_voigt\n"
        "    _cell, _ = read_crystal_structure(\"POSCAR\", interface_mode=\"vasp\")\n"
        "    _sg = spglib.get_spacegroup((_cell.cell, _cell.scaled_positions, _cell.numbers), symprec=1e-3)\n"
        "    CRYSTAL = crystal_system_from_spacegroup(int(_sg.split(\"(\")[-1].rstrip(\")\")))\n"
        "    SG = _sg\n"
        "    out[\"crystal_system\"] = CRYSTAL; out[\"spacegroup\"] = SG\n"
        "    if META.get(\"dim\") == \"2d\":\n"
        "        PLANE = np.array(_cell.cell, float)[int(META.get(\"vac_axis\", 2))]\n"
        "except Exception as _e:\n"
        "    out.setdefault(\"warnings\", []).append(\"kappa symmetry audit unavailable: %s\" % _e)\n"
        "def _audit(K3):\n"
        "    if not CRYSTAL:\n"
        "        return None\n"
        "    voigt = [K3[0,0], K3[1,1], K3[2,2], K3[1,2], K3[0,2], K3[0,1]]\n"
        "    aud = audit_kappa_voigt(voigt, crystal_system=CRYSTAL, plane_normal=PLANE)\n"
        "    return {k: aud.get(k) for k in (\"eigenvalue_ratio\", \"diagonal_relative_mismatch\",\n"
        "                                 \"offdiag_relative_magnitude\", \"determinant_ratio\",\n"
        "                                 \"threshold\", \"gate\", \"crystal_system\")}\n"
        "expected = [\"kappa-m\" + \"\".join(m.split()) + \".hdf5\" for m in META.get(\"meshes\", [])]\n"
        "files = expected or sorted(glob.glob(\"kappa-m*.hdf5\"))\n"
        "for f in files:\n"
        "    if \"kappa-mfp\" in f:\n"
        "        continue\n"
        "    try:\n"
        "        with h5py.File(f, \"r\") as h:\n"
        "            T = np.array(h[\"temperature\"]); Kv = np.array(h[\"kappa\"])\n"
        "        # phono3py kappa-m*.hdf5 的 kappa 是 (n_T,6) Voigt [xx,yy,zz,yz,xz,xy]，不是 (n_T,3,3)\n"
        "        raw3 = [[[float(Kv[i,0]),float(Kv[i,5]),float(Kv[i,4])],\n"
        "                 [float(Kv[i,5]),float(Kv[i,1]),float(Kv[i,3])],\n"
        "                 [float(Kv[i,4]),float(Kv[i,3]),float(Kv[i,2])]] for i in range(len(T))]\n"
        "        raw = [[float(Kv[i,0]),float(Kv[i,1]),float(Kv[i,2])] for i in range(len(T))]\n"
        "        rec = {\"file\": f, \"mesh\": os.path.basename(f)[7:].replace(\".hdf5\",\"\"),\n"
        "               \"temperatures\": T.tolist(), \"kappa_xx_yy_zz\": raw,\n"
        "               \"kappa_tensor_3x3\": raw3}\n"
        "        j = int(np.argmin(np.abs(T-300.0)))\n"
        "        rec[\"kappa_300K_xx_yy_zz\"] = raw[j]\n"
        "        rec[\"kappa_300K_tensor_3x3\"] = raw3[j]\n"
        "        aud = _audit(np.array(raw3[j]))\n"
        "        if aud:\n"
        "            rec[\"kappa_symmetry_audit_300K\"] = aud\n"
        "        if abs(FACTOR-1.0)>1e-9:\n"
        "            nrm=[[v*FACTOR for v in r] for r in raw]\n"
        "            rec[\"kappa_2d_normalized_xx_yy_zz\"]=nrm\n"
        "            rec[\"kappa_2d_normalized_300K_xx_yy_zz\"]=nrm[j]\n"
        "        out[\"runs\"].append(rec)\n"
        "    except Exception as e:\n"
        "        out.setdefault(\"errors\", []).append(\"%s: %s\" % (f, e))\n"
        "out[\"KAPPA_DONE\"] = bool(out[\"runs\"]) and not out.get(\"errors\") and (not expected or len(out[\"runs\"]) == len(expected))\n"
        "if out[\"runs\"]:\n"
        "    r = out[\"runs\"][-1]\n"
        "    out[\"kappa_300K_xx_yy_zz\"] = r[\"kappa_300K_xx_yy_zz\"]\n"
        "    out[\"mesh_last\"] = r[\"mesh\"]\n"
        "    # ---- 下游（device-thermal 等）读的顶层字段：与 kl-dft-cpu S6 口径一致 ----\n"
        "    # 此前这些只在 runs[i] 里，顶层缺失会让 device-thermal 的 read_kappa() 返回\n"
        "    # None 并静默回退内置材料库（文献值）——看着接上了其实没接。\n"
        "    out[\"file\"] = r[\"file\"]\n"
        "    out[\"temperatures\"] = r[\"temperatures\"]\n"
        "    out[\"kappa_xx_yy_zz\"] = r[\"kappa_xx_yy_zz\"]\n"
        "    if \"kappa_2d_normalized_300K_xx_yy_zz\" in r:\n"
        "        out[\"kappa_2d_normalized_300K_xx_yy_zz\"]=r[\"kappa_2d_normalized_300K_xx_yy_zz\"]\n"
        "    if \"kappa_symmetry_audit_300K\" in r:\n"
        "        out[\"kappa_symmetry_audit_300K\"] = r[\"kappa_symmetry_audit_300K\"]\n"
        "# ---- MFP 累积 κ 诊断（MFP_CUMULATIVE/--mfp 时才有 kappa-mfp.hdf5）----\n"
        "for f in sorted(glob.glob(\"kappa-mfp*.hdf5\")):\n"
        "    try:\n"
        "        with h5py.File(f, \"r\") as h:\n"
        "            mfp = np.array(h[\"mfp\"]); km = np.array(h[\"kappa_mfp\"])\n"
        "            T = np.array(h[\"temperature\"])\n"
        "        j = int(np.argmin(np.abs(T-300.0)))\n"
        "        cum = np.array(km[j], float)\n"
        "        tot = np.trace(cum[-1]) / 3.0\n"
        "        frac = np.array([np.trace(c)/3.0/tot for c in cum]) if tot > 0 else np.zeros(len(cum))\n"
        "        def _mfp_at(pct):\n"
        "            idx = np.argmax(frac >= pct)\n"
        "            return float(mfp[idx]) if frac[idx] >= pct else None\n"
        "        out[\"mfp_cumulative_300K\"] = {\"file\": f,\n"
        "            \"mfp_ang\": [float(x) for x in mfp],\n"
        "            \"cumulative_fraction\": [float(x) for x in frac],\n"
        "            \"mfp_at_50pct_ang\": _mfp_at(0.5), \"mfp_at_90pct_ang\": _mfp_at(0.9)}\n"
        "    except Exception as e:\n"
        "        out.setdefault(\"errors\", []).append(\"mfp: %s\" % e)\n"
        "with open(\"kappa_summary.json\", \"w\") as stream:\n"
        "    json.dump(out, stream, ensure_ascii=False, indent=2)\n"
        "print(\"KAPPA_DONE\" if out[\"KAPPA_DONE\"] else \"NO_KAPPA\")\n"
        "PY")



def meshes(conf, params, dim, vac_axis=2):
    # 每套网格都过 mesh_str：3D 材料被误写成 "N N 1"（2D 残留）时自动纠正成 "N N N"，
    # 否则 phono3py 报 "Grid symmetry is broken"；2D 材料则把真空轴压成 1。
    scan = str(conf["MESH_SCAN"] or "").strip()
    if scan:
        return [kc.mesh_str(m.strip().split(), dim, vac_axis)
                for m in scan.split(";") if m.strip()]
    return [kc.mesh_str(
        (conf["MESH_OVERRIDE"] or params.get("MESH") or "24 24 24").split(),
        dim, vac_axis)]




def _shengbte_control_from_poscar(poscar, SUPERCELL, ngrid, tmin, tmax, tstep,
                                  scalebroad, isotope, use_nac, out_path):
    """从 POSCAR(原胞) 写 ShengBTE CONTROL。格式对齐 kl-dft lattice_kappa._write_shengbte_control。"""
    from ase.io import read as ase_read
    from ase.data import chemical_symbols as _cs
    atoms = ase_read(str(poscar), format="vasp")
    cell = atoms.cell
    numbers = atoms.numbers; unique = sorted(set(numbers))
    syms = [_cs[z] for z in unique]; kd = {z: i + 1 for i, z in enumerate(unique)}
    types = [kd[z] for z in numbers]; spos = atoms.get_scaled_positions(wrap=True)
    scell = [int(x) for x in SUPERCELL]
    if len(scell) == 9:
        # 一般矩阵超胞：ShengBTE 的 control 只写对角 scell。对角化的矩阵照常提取，
        # 真·非对角就直接报错——绝不静默近似成对角。
        _m = np.array(scell, dtype=int).reshape(3, 3)
        if np.count_nonzero(_m - np.diag(np.diag(_m))):
            sys.exit("[ERROR] ShengBTE 的 CONTROL 文件里 scell 是三个整数，物理上"
                     "只表达对角超胞；当前 SUPERCELL=%r 是一般矩阵，没法写进去。\n"
                     "        两条路：(1) 该步 SUPERCELL/FC2_SUPERCELL 用 \"n n n\" 对角扩胞；"
                     "(2) 换 kappa 后端（phono3py/phonopy 能直接吃 3×3 矩阵）。"
                     % (SUPERCELL,))
        scell = [int(x) for x in np.diag(_m)]
    ng = [int(x) for x in ngrid.split()]
    L = ["&allocations", "  nelements=%d," % len(unique), "  natoms=%d," % len(atoms),
         "  ngrid(:)=%d %d %d" % (ng[0], ng[1], ng[2]), "&end", "&crystal", "  lfactor=0.1,"]
    for r in range(3):
        L.append("  lattvec(:,%d)=" % (r+1) + " ".join("%.10f" % cell[r, i] for i in range(3)) + ",")
    L.append("  elements=" + " ".join('"%s"' % s for s in syms))
    L.append("  types=" + " ".join(str(t) for t in types) + ",")
    for idx, p in enumerate(spos):
        L.append("  positions(:,%d)=" % (idx+1) + " ".join("%.10f" % x for x in p) + ",")
    if use_nac:
        from phonopy.file_IO import parse_BORN
        from phonopy.structure.atoms import PhonopyAtoms
        primitive = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell,
                                 scaled_positions=spos)
        nac = parse_BORN(primitive, filename=str(Path(out_path).parent / "BORN"))
        if not nac:
            raise ValueError("无法读取 ShengBTE NAC 参数")
        for j in range(3):
            L.append("  epsilon(:,%d)=" % (j+1) + " ".join(str(x) for x in nac["dielectric"][:, j]) + ",")
        for atom, born in enumerate(nac["born"]):
            for j in range(3):
                L.append("  born(:,%d,%d)=" % (j+1, atom+1) + " ".join(str(x) for x in born[:, j]) + ",")
    L.append("  scell(:)=%d %d %d" % (scell[0], scell[1], scell[2]))
    L += ["&end", "&parameters", "  T_min=%.1f" % float(tmin),
          "  T_max=%.1f" % float(tmax), "  T_step=%.1f" % float(tstep),
          "  scalebroad=%s" % scalebroad, "&end", "&flags",
          "  autoisotopes=%s," % ("T" if isotope else "F"),
          "  convergence=F,",          # kl-mlff 链先走 RTA（迭代 CONV 内存/耗时更大）
          "  nonanalytic=%s," % ("T" if use_nac else "F"),
          "  nanowires=F,", "&end"]
    Path(out_path).write_text("\n".join(L) + "\n", encoding="utf-8")


    print("[OK] CONTROL <- %s" % out_path.name)


def _ensure_shengbte_fc(src, sb_dir):
    """从 step3_fc 取 ShengBTE 力常数到 sb_dir。
    优先 step3_fc/shengbte/（fc_fit_driver 统一导出）；缺则现场从 fc2/3.hdf5 用 hiphive 转换。"""
    sb_dir = Path(sb_dir)
    sb_dir.mkdir(exist_ok=True)
    need = ("FORCE_CONSTANTS_2ND", "FORCE_CONSTANTS_3RD")
    ready = [sb_dir / f for f in need]
    if all(p.is_file() for p in ready):
        return ready
    # 现场转换：读 fc2.hdf5/fc3.hdf5 → hiphive ForceConstants → 导出
    print("[..] %s/shengbte/ 缺 ShengBTE 力常数，从 fc2/fc3.hdf5 现场转换（hiphive）" % src)
    import ase, h5py, numpy as np
    import phono3py
    from hiphive import ForceConstants
    yaml = ("phono3py_params.yaml" if (src / "phono3py_params.yaml").is_file()
            else "phono3py_disp.yaml")
    ph3 = phono3py.load(str(src / yaml), produce_fc=False, log_level=0)
    prim, sc = ph3.phonon_primitive, ph3.supercell
    with h5py.File(str(src / "fc2.hdf5"), "r") as h:
        fc2 = np.asarray(h["fc2" if "fc2" in h else "force_constants"][()])
    with h5py.File(str(src / "fc3.hdf5"), "r") as h:
        fc3 = np.asarray(h["fc3"][()])
    prim_ase = ase.Atoms(symbols=prim.symbols, cell=prim.cell,
                         scaled_positions=prim.scaled_positions, pbc=True)
    sc_ase = ase.Atoms(symbols=sc.symbols, cell=sc.cell,
                       scaled_positions=sc.scaled_positions, pbc=True)
    sc2 = ph3.phonon_supercell
    if sc2 is None:
        sc2 = sc
    sc2_ase = ase.Atoms(symbols=sc2.symbols, cell=sc2.cell,
                        scaled_positions=sc2.scaled_positions, pbc=True)
    fcs2 = ForceConstants.from_arrays(sc2_ase, fc2_array=fc2)
    fcs = ForceConstants.from_arrays(sc_ase, fc3_array=fc3)
    fcs2.write_to_phonopy(str(sb_dir / "FORCE_CONSTANTS_2ND"), format="text")
    fcs.write_to_shengBTE(str(sb_dir / "FORCE_CONSTANTS_3RD"), prim_ase)
    print("[OK] FORCE_CONSTANTS_2ND/3RD <- fc2/fc3.hdf5（hiphive 导出，格式已验证）")
    return [sb_dir / f for f in need]


def build_shengbte_submit(cwd, out, src, conf, params, dim, vac_axis=2,
                          use_nac=False, factor=1.0, meta=None):
    """SOLVER=shengbte：备好 ShengBTE 输入（力常数+CONTROL），渲染 submit_shengbte 提交。"""
    if len(meshes(conf, params, dim, vac_axis)) != 1:
        sys.exit("[ERROR] shengbte 求解器暂只支持单套网格（MESH_SCAN 留空）")
    sb_dir = src / "shengbte"
    sb_dir.mkdir(exist_ok=True)
    _ensure_shengbte_fc(src, sb_dir)
    for f in ("POSCAR", kc.KL_PARAMS, kc.METHOD_FILE):
        if (src / f).is_file() and not (out / f).exists():
            shutil.copyfile(str(src / f), str(out / f))
    # CONTROL（依赖本步运行参数 → 在 step4 生成）
    fc2, fc3 = _ensure_shengbte_fc(src, sb_dir)
    shutil.copyfile(str(fc2), str(out / "FORCE_CONSTANTS_2ND"))
    shutil.copyfile(str(fc3), str(out / "FORCE_CONSTANTS_3RD"))
    poscar = out / "POSCAR" if (out / "POSCAR").is_file() else src / "POSCAR"
    supercell = (params.get("FC2_SUPERCELL") or params.get("SUPERCELL") or "4 4 4").split()
    _shengbte_control_from_poscar(
        poscar, supercell, meshes(conf, params, dim, vac_axis)[0],
        conf["T_MIN"], conf["T_MAX"], conf["T_STEP"],
        conf["SCALEBROAD"], conf["ISOTOPE"], use_nac,
        out / "CONTROL")
    # 抽 κ 小脚本：ShengBTE 输出 BTE.KappaTensorVsT_RTA（CONV 预留）
    import json
    extract = "python - <<'PY'\n"
    extract += "import json\nMETA=json.loads(%r)\nFACTOR=%r\n" % (json.dumps(meta or {}), factor)
    extract += (
        "import json, math\n"
        "from pathlib import Path\n"
        "f = Path(\"BTE.KappaTensorVsT_RTA\")\n"
        "d = {\"KAPPA_DONE\": False, \"solver\": \"shengbte\"}\n"
        "d.update(META)\n"
        "if f.is_file():\n"
        "    rows = [l.split() for l in f.read_text().splitlines() if l.strip() and not l.lstrip().startswith(\"#\")]\n"
        "    if rows:\n"
        "        temperatures = [float(r[0]) for r in rows]\n"
        "        raw = [[float(r[1]), float(r[5]), float(r[9])] for r in rows]\n"
        "        if not all(math.isfinite(v) for r in raw for v in r):\n"
        "            raise ValueError(\"non-finite ShengBTE kappa\")\n"
        "        j = min(range(len(temperatures)), key=lambda i: abs(temperatures[i]-300.0))\n"
        "        d.update(KAPPA_DONE=True, source=str(f), temperatures=temperatures, kappa_xx_yy_zz=raw, kappa_300K_xx_yy_zz=raw[j])\n"
        "        if META.get(\"dim\") == \"2d\":\n"
        "            normalized = [[v*FACTOR for v in r] for r in raw]\n"
        "            d[\"kappa_2d_normalized_xx_yy_zz\"] = normalized\n"
        "            d[\"kappa_2d_normalized_300K_xx_yy_zz\"] = normalized[j]\n"
        "Path(\"kappa_summary.json\").write_text(json.dumps(d, ensure_ascii=False, indent=2))\n"
        "print(\"KAPPA_DONE\" if d[\"KAPPA_DONE\"] else \"NO_KAPPA\")\n"
        "PY"
    )

    exe = str(conf["SHENGBTE_EXE"] or "ShengBTE").strip()
    # MPI/OMP 布局：一个 rank 一个 NUMA **节点**（rank 数必须 ≤ NUMA 节点数，否则
    # OpenMPI 循环复用节点、多个 rank 挤同一组核）。rank=1 会退化成串行（ShengBTE
    # 靠 MPI 按 q 点并行），必须显式算出来，不能沿用模板里的常量。
    _cpus_per_numa = int(conf["SHENGBTE_CORES_PER_NUMA"] or 24)
    _total = int(conf["SHENGBTE_TOTAL_CORES"] or 96)
    _nt_raw = str(conf["SHENGBTE_NTASKS"] or "auto").strip()
    if _nt_raw and _nt_raw.lower() != "auto":
        try:
            _ntasks = int(_nt_raw)
        except ValueError:
            sys.exit("[ERROR] SHENGBTE_NTASKS=%r 不是整数也不是 auto" % _nt_raw)
    else:
        _ntasks = max(1, _total // _cpus_per_numa)
    if _ntasks <= 1:
        sys.exit("[ERROR] SHENGBTE_NTASKS 算出来是 %d —— ShengBTE 靠 MPI 按 q 点并行，"
                 "单进程等于串行（实测 11.8 h 零产物）。请检查 step.conf 的 "
                 "SHENGBTE_TOTAL_CORES=%d / SHENGBTE_CORES_PER_NUMA=%d。"
                 % (_ntasks, _total, _cpus_per_numa))
    print("[..] ShengBTE 布局：%d MPI ranks x %d OMP threads = %d 核"
          % (_ntasks, _cpus_per_numa, _ntasks * _cpus_per_numa))
    here = Path(__file__).resolve().parent
    tpl = kc.resolve_submit(here, "submit_shengbte_klm")
    kc.write_submit(tpl, out / "submit.sh",
                    {"JOBNAME": kc.new_jobname(cwd, "S4sheng"),
                     "CONDA_SH": conf["CONDA_SH"] or kc.DEFAULT_CONDA_SH,
                     "CONDA_ENV": conf["CONDA_ENV"] or kc.DEFAULT_CONDA_ENV,
                     "SHENGBTE_EXE": exe,
                     "NTASKS": str(_ntasks),
                     "CPUS_PER_TASK": str(_cpus_per_numa),
                     "SB_EXTRACT": extract})
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    print("[DONE] %s：submit.sh 就绪（ShengBTE RTA），跑完写 kappa_summary.json" % OUTDIR)

# 计算节点上用：--wgp 写出 ir_grid_points.yaml（不可约格点编号 + 权重）。phono3py 的
#   --gp 接的是「网格点编号列表」而不是分块参数，所以先把编号均分成 <=GP_SPLIT 份、
#   每份一行写进 gp_chunks.txt，shell 再逐行 --gp="<该行>" 起并行 --write-gamma；
#   汇总步 --read-gamma 不带 --gp（不带 --gp 时读全部不可约格点的 gamma）。
_GP_SPLIT_PY = r'''python - <<'PY'
import re
_text = open("ir_grid_points.yaml", encoding="utf-8").read()
_gps = [int(x) for x in re.findall(r"-\s*grid_point:\s*(\d+)", _text)]
_n = __GP_SPLIT__
_groups = [[] for _ in range(max(1, _n))]
for _k, _g in enumerate(_gps):
    _groups[_k % _n].append(_g)
_groups = [_g for _g in _groups if _g]
if not _groups:
    raise SystemExit("[gp] ir_grid_points.yaml 里没解析到不可约格点，--gp 无点可算")
with open("gp_chunks.txt", "w", encoding="utf-8") as _f:
    for _g in _groups:
        _f.write(" ".join(str(x) for x in _g) + "\n")
print("[gp] %d ir grid points -> %d chunk(s), sizes=%s"
      % (len(_gps), len(_groups), [len(_g) for _g in _groups]))
PY'''


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    src = cwd / SRC
    if not src.is_dir():
        sys.exit("[ERROR] 找不到 %s" % SRC)

    ps = src / "phonon_summary.json"
    if ps.is_file():
        import json
        try:
            if not json.loads(ps.read_text()).get("stable", True):
                sys.exit("[ERROR] step3 判定声子谱有虚频，κ 无物理意义，已中止。"
                         "先按 step3 的排查顺序处理。")
        except ValueError:
            pass

    for f in ("fc2.hdf5", "fc3.hdf5"):
        if not (src / f).is_file():
            sys.exit("[ERROR] %s 缺 %s（step3 没拟合成）" % (SRC, f))
        kc.link_or_copy(src / f, out / f)          # fc3.hdf5 可以上 GB，软链不复制
    yaml = ("phono3py_params.yaml" if (src / "phono3py_params.yaml").is_file()
            else "phono3py_disp.yaml")
    for f in (yaml, "BORN", "POSCAR", kc.KL_PARAMS, kc.METHOD_FILE):
        if (src / f).is_file():
            shutil.copyfile(str(src / f), str(out / f))
    use_nac = (out / "BORN").is_file()

    params = kc.read_kl_params(out / kc.KL_PARAMS)
    # 维度 + 2D NAC 门槛（phono3py 只有 3D 方案；2D 极性材料 LO-TO 在 q->0 应趋零，
    #   3D-NAC 是随真空变化的伪劈裂）。auto：2D 不用 NAC、3D 随 BORN；on/off 强制。
    dim = (params.get("DIM") or "").lower()
    try:
        _d, vac_axis = kc.resolve_dim(out / "POSCAR", dim or "auto")
        dim = dim or _d
    except Exception:
        vac_axis = None
    nac_mode = str(conf["KAPPA_NAC"] or "auto").strip().lower()
    if nac_mode in ("on", "true", "1", "yes"):
        pass
    elif nac_mode in ("off", "false", "0", "no"):
        use_nac = False
    elif dim == "2d" and use_nac:
        use_nac = False
        print("[WARN] 2D + KAPPA_NAC=auto：默认不对 2D 施加 phono3py 的 3D-NAC")
        print("       （LO-TO 在 2D 应趋零，3D 方案是随真空变化的伪劈裂）。")
        print("       要强制用设 KAPPA_NAC=on；正确的 2D-NAC 需用 QE 的 2D-DFPT。")

    solver = str(conf["SOLVER"] or "phono3py").strip().lower()
    if solver not in ("phono3py", "shengbte"):
        sys.exit("[ERROR] SOLVER 只允许 phono3py / shengbte")
    ms = meshes(conf, params, dim, vac_axis)
    ts = " ".join(str(t) for t in range(conf["T_MIN"], conf["T_MAX"] + 1, conf["T_STEP"]))
    bte = str(conf["BTE"] or "rta").lower()
    if bte not in ("rta", "lbte"):
        sys.exit("[ERROR] BTE 只允许 rta / lbte")
    print("[..] SOLVER=%s  网格=%s  温度=%d~%dK  NAC=%s  DIM=%s"
          % (solver, " | ".join(ms), conf["T_MIN"], conf["T_MAX"], use_nac, dim or "?"))
    ph3_solver = "--br" if bte == "rta" else "--lbte"

    # 2D κ 厚度归一化因子（3D 时 factor=1、不归一）
    factor, meta = 1.0, {"dim": dim or "?", "meshes": ms, "nac": use_nac}
    if dim == "2d":
        try:
            factor, m2 = two_d_norm_factor(out / "POSCAR", vac_axis, conf["KAPPA_2D_THICKNESS"])
            meta.update(m2)
            meta["vac_axis"] = vac_axis          # extract 用它投影面内块做 κ 对称审计
            print("[..] 2D κ 归一化：Lz=%.3f d=%.3f factor=Lz/d=%.4f (%s)"
                  % (meta["Lz_ang"], meta["thickness_d_ang"], factor, meta["thickness_convention"]))
        except Exception as e:
            print("[WARN] 2D 归一化因子算失败，只出原始 κ：%s" % e)
    # ---- 下游契约：把 2D 有效厚度落成 thickness_2d.json ----------------------
    # device-thermal（器件传热）用 read_thickness() 取沟道厚度；kl/ke 共用契约还要在
    # 材料根再留一份。纯增量：不碰 kappa_summary.json 的任何既有字段，也不改归一化逻辑。
    if dim == "2d" and meta.get("thickness_d_A"):
        try:
            import json as _json
            _t2d = {"thickness_d_A": meta.get("thickness_d_A"),
                    "h_perp_A": meta.get("h_perp_A"),
                    "vacuum_gap_A": meta.get("vacuum_gap_A"),
                    "thickness_convention": meta.get("thickness_convention"),
                    "kappa_2d_norm_factor": meta.get("kappa_2d_norm_factor"),
                    "Lz_ang": meta.get("Lz_ang"),
                    "thickness_d_ang": meta.get("thickness_d_ang"),
                    "atomic_zspan_ang": meta.get("atomic_zspan_ang"),
                    "dim": "2d", "source": "kl-mlff:step4_kappa"}
            (out / "thickness_2d.json").write_text(
                _json.dumps(_t2d, ensure_ascii=False, indent=2), encoding="utf-8")
            print("[OK] thickness_2d.json d=%.3f A -> %s"
                  % (_t2d["thickness_d_A"], out / "thickness_2d.json"))
            import thickness_2d as _t2dmod
            print(_t2dmod.write_material_level(cwd.parent, _t2d))
        except Exception as e:
            print("[WARN] thickness_2d.json 写出失败（不影响 κ）：%s" % e)
    if solver == "shengbte":
        build_shengbte_submit(cwd, out, src, conf, params, dim, vac_axis,
                              use_nac=use_nac, factor=factor, meta=meta)
        return
    if bte == "lbte" and len(ms) > 1:
        print("[WARN] LBTE 要存整个碰撞矩阵，内存 ~O(N_mode²)。多网格扫描请用 rta。")

    gp = int(conf["GP_SPLIT"] or 0)
    if gp >= 2:
        print("[..] GP_SPLIT=%d：每套网格先 --wgp 列不可约格点，均分 %d 份并行 --write-gamma，"
              "再 --read-gamma 汇总" % (gp, gp))
    mfp_on = bool(conf["MFP_CUMULATIVE"])
    if mfp_on:
        print("[..] MFP_CUMULATIVE：加 --mfp，extract 收累积 κ vs 自由程诊断")

    def _base(m):
        """一条 phono3py-load 的公共前缀（网格/温度/同位素/NAC/MFP）。"""
        return '%s %s --mesh %s --ts="%s"%s%s%s' % (
            yaml, ph3_solver, m, ts,
            " --isotope" if conf["ISOTOPE"] else "",
            " --nac" if use_nac else " --nonac",
            " --mfp" if mfp_on else "")

    def _wgp(m):
        """--wgp 前缀：只要网格 + NAC 开关，写出 ir_grid_points.yaml（不可约格点编号）。"""
        return '%s --mesh %s%s' % (yaml, m, " --nac" if use_nac else " --nonac")

    lines = []
    for mi, m in enumerate(ms):
        extra = str(conf["EXTRA_ARGS"] or "").strip()
        extra_s = (" " + extra) if extra else ""
        if gp >= 2:
            tag = "kappa_gp_%d" % mi      # 用网格序号命名日志，多套网格不互相覆盖
            # ① --wgp 拿到不可约格点编号；失败即停（否则后面 --gp 无点可算）
            lines.append('phono3py-load %s --wgp > %s_wgp.log 2>&1 '
                         '|| { echo "[gp] mesh %s: --wgp 失败"; exit 1; }'
                         % (_wgp(m), tag, m))
            # ② 把编号均分成 <=GP_SPLIT 份，每份一行（--gp 接编号列表，不是第 i/n 份）
            lines.append(_GP_SPLIT_PY.replace("__GP_SPLIT__", str(gp)))
            lines.append('[ -s gp_chunks.txt ] || { echo "[gp] mesh %s: 格点分块文件缺失/为空，'
                         '中止"; exit 1; }' % m)
            # ③ 每份一个后台 --write-gamma；任一份失败就停，不用缺 γ 的数据汇总
            lines.append('_pids=""; _n=0')
            lines.append('while IFS= read -r _gps; do')
            lines.append('  _n=$((_n+1))')
            lines.append('  phono3py-load %s --gp="$_gps" --write-gamma%s '
                         '> %s_$_n.log 2>&1 &' % (_base(m), extra_s, tag))
            lines.append('  _pids="$_pids $!"')
            lines.append('done < gp_chunks.txt')
            lines.append('_fail=0; for _p in $_pids; do wait "$_p" || _fail=1; done')
            lines.append('if [ $_fail -ne 0 ]; then echo "[gp] mesh %s: 有 --write-gamma '
                         '分块失败，停止汇总（不用缺 γ 的结果出 κ）"; exit 1; fi' % m)
            lines.append('rm -f gp_chunks.txt')
            lines.append('echo "[gp] mesh %s: gamma chunks done, collecting..."' % m)
            # ④ 汇总：--read-gamma 不带 --gp → 读全部不可约格点的 γ
            lines.append('phono3py-load %s --read-gamma%s 2>&1 | tee -a phono3py_kappa.log'
                         % (_base(m), extra_s))
        else:
            lines.append('phono3py-load %s %s 2>&1 | tee -a phono3py_kappa.log'
                         % (_base(m), extra_s))
    cmd = "\n".join(lines) + "\n" + build_extract(factor, meta)

    here = Path(__file__).resolve().parent
    # κ 对称审计模块在 extract 的计算节点进程里 import，必须随 submit.sh 一起推进
    # step 目录（gen_need 只服务本地 gen 阶段、不推送计算节点辅助脚本——
    # 对照 gen_step2 显式复制 mc_rattle.py/mc_rattle_disp.py 的做法）。
    for _aux in ("kappa_validation.py",):
        if (here / _aux).is_file():
            shutil.copyfile(str(here / _aux), str(out / _aux))
    tpl = kc.resolve_submit(here, "submit_p3py")
    kc.write_submit(tpl, out / "submit.sh",
                    {"JOBNAME": kc.new_jobname(cwd, "S4kappa"),
                     "CONDA_SH": conf["CONDA_SH"] or kc.DEFAULT_CONDA_SH,
                     "CONDA_ENV": conf["CONDA_ENV"] or kc.DEFAULT_CONDA_ENV,
                     "P3PY_CMD": cmd})
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    print("[DONE] %s：submit.sh 就绪（%d 套网格），跑完写 kappa_summary.json"
          % (OUTDIR, len(ms)))


if __name__ == "__main__":
    main()
