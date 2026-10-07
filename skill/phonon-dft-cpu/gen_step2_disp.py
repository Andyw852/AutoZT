#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_disp.py —— phonopy 2 阶对称有限位移 + 超胞单点取力（扇出 disp-*）。

从 step1 弛豫结构接力，扩超胞，phonopy 生成 2 阶位移超胞，每个位移一个
disp-NNNNN 子目录单点取力。产出 phonopy_disp.yaml + disp-00001..N。

2026-09-27 与 kl-dft-cpu S4_disp 对齐（只保留 2 阶部分，去掉 fc3/ALM/phono3py）：
  · 取力 INCAR 精度走 FORCE_PREC / FORCE_EDIFF / FORCE_LREAL + ADDGRID（原模板
    注释写 Accurate/1E-8/.FALSE.，正文却是 Normal/Auto/1E-7，自相矛盾）；
  · relay 后结构对称性审计（symmetry_gate），检出数值微畸变即对称化；
  · 已有数据集做结构指纹 + 配置四重比对，不符则**非破坏**移走重生成；
  · 扇出前帧间一致性硬闸（超胞同胞、KPOINTS 网格唯一）；
  · 偶极修正从 S1 INCAR 继承、泛函与 step1 一致性硬检查；
  · 2D 真空隙/残余面内应力门禁。
"""
import glob
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kl_common as kc
from dim_common import require_dim
import stepconf

OUTDIR = "step2_disp"
STEP   = "step2_disp"
PREV   = ["step1_std_opt"]

SPEC = {
    "FUNC":         ("pbesol", "str"),
    "KSPACING":     ("0.04",  "str"),
    "KSCHEME":      ("2",     "str"),
    "ENCUT":        (None,    "int"),
    "ENCUT_FACTOR": (1.5,     "float"),
    "VASPKIT_EXE":  ("vaspkit", "str"),
    "SUPERCELL":    (None,    "words"),
    "MIN_SC_LEN":   (12.0,    "float"),
    "MAX_MULTIPLE": (6,       "int"),
    "SUPERCELL_SYMMETRY": ("strict", "str"),  # 超胞须保晶体全部点群操作：strict 拦 | warn 告警 | off
    "FD_DISTANCE":  (0.01,    "float"),
    "PHONON_MESH":  ("20 20 20", "str"),
    "MAX_DISP":     (200,     "int"),
    # --- 取力精度（对齐 kl-dft-cpu S4；力直接进 fc2，别降）---
    #   phonopy 官方超胞取力示例口径：PREC=Accurate / EDIFF=1E-8 / LREAL=.FALSE.。
    #   LREAL=Auto 的实空间投影带格点噪声，是 fc 假软模的常见来源。超胞几百原子时
    #   可换 Auto，但必须先在 3 帧上与 .FALSE. 比力（最大偏差 < 1 meV/Å）并设 ROPT。
    "FORCE_PREC":   ("Accurate", "str"),
    "FORCE_EDIFF":  ("1E-8",     "str"),
    "FORCE_LREAL":  (".FALSE.",  "str"),
    # --- 2D 残余面内应力门禁（层内口径 σ_VASP×h⊥/d，kbar）。0 = 关掉门禁。---
    "STRESS_2D_THR": (0.5, "float"),
    # --- 帧间一致性硬闸（已有 disp-* 必须与当前计划同胞同网格）。off = 跳过。---
    "FRAME_GATE":   ("on", "str"),
}
# 结构体检第五条：SYMMETRY_AUDIT / SYMMETRY_SYMMETRIZE 等键（见 kl_common.SYMMETRY_SPEC）
SPEC.update(kc.SYMMETRY_SPEC)


# ==========================================================================
# 结构指纹 / 数据集幂等 / 非破坏隔离（对齐 kl-dft-cpu gen_step4_disp）
# ==========================================================================
def _poscar_cell(path):
    """读 POSCAR 的 3x3 晶格矩阵（纯文本解析）；失败返回 None。"""
    try:
        ln = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
        scale = float(ln[1].split()[0])
        return [[float(x) * scale for x in ln[i].split()[:3]] for i in (2, 3, 4)]
    except Exception:                                   # noqa: BLE001
        return None


def _disp_yaml_cell(path):
    """读 phonopy_disp.yaml 里的 unit_cell.lattice；失败返回 None。"""
    try:
        import yaml
        d = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        return [[float(x) for x in row] for row in d["unit_cell"]["lattice"]]
    except Exception:                                   # noqa: BLE001
        return None


def check_existing_matches_input(out):
    """结构指纹校验：已有位移/力数据必须属于【当前输入结构】。

    retry 对 fanout 只补缺失子目录，已有 disp-* 一律不动；而生成逻辑见到
    phonopy_disp.yaml + POSCAR-* 就"幂等跳过"。S1 换了结构后重投 S2，旧位移、
    旧帧会被整套沿用 —— S3 拿旧结构的力配旧位移拟合，内部自洽却与当前结构无关。
    """
    y = out / "phonopy_disp.yaml"
    if not y.is_file():
        return
    old, new = _disp_yaml_cell(y), _poscar_cell(out / "POSCAR")
    if old is None or new is None:
        print("[WARN] 结构指纹校验跳过（读不到 phonopy_disp.yaml 或 POSCAR）")
        return
    dmax = max(abs(old[i][j] - new[i][j]) for i in range(3) for j in range(3))
    if dmax > 1e-3:
        n_disp = len(glob.glob(str(out / "disp-*")))
        sys.exit(
            "[ERROR] 已有位移/力数据属于【旧结构】—— phonopy_disp.yaml 的单胞与当前输入 "
            "差异最大 %.4f Å（>1e-3 即判定不一致）。\n"
            "        当前有 %d 个 disp-* 目录。retry 对 fanout 只补缺失帧，这些旧帧会被"
            "整套沿用，\n        S3 会把旧结构的力配上（同样旧的）位移静默混拟合。\n"
            "        处置：确认已归档后清空本步产物，再 retry：\n"
            "          rm -rf %s/disp-* %s/POSCAR-* %s/phonopy_disp.yaml %s/SPOSCAR\n"
            "        或直接 autozt -tt phonon-dft-cpu -p <材料> -j S2_disp rerun（推倒重来）。"
            % (dmax, n_disp, out, out, out, out))


def _dataset_matches(out, reps, conf):
    """现有数据集是否【已经】对应当前输入 —— 四重比对：
      ① 单胞 == POSCAR（tol 1e-8）  ② 超胞矩阵 == 配置的 reps
      ③ 位移幅度 FD_DISTANCE == 当前 conf  ④ disp_plan.json 留档一致。
    """
    y = out / "phonopy_disp.yaml"
    if not y.is_file():
        return False
    old, new = _disp_yaml_cell(y), _poscar_cell(out / "POSCAR")
    if old is None or new is None:
        return False
    if max(abs(old[i][j] - new[i][j]) for i in range(3) for j in range(3)) > 1e-8:
        return False
    try:
        import yaml
        _d = yaml.safe_load(y.read_text(encoding="utf-8"))
        _sm = _d["supercell_matrix"]
        # phonopy 的 phonopy_disp.yaml 一般存 3×3 嵌套；兼容扁平 9 值写法
        if len(_sm) == 3 and all(hasattr(r, "__len__") for r in _sm):
            sm = [[int(x) for x in row] for row in _sm]
        else:
            _v = [int(x) for x in _sm]
            sm = [_v[0:3], _v[3:6], _v[6:9]]
    except Exception:                                # noqa: BLE001
        return False
    _r = np.array(reps, dtype=int)
    want = (_r.reshape(3, 3) if _r.size == 9 else np.diag(_r)).tolist()
    if sm != [[int(x) for x in row] for row in want]:
        return False
    pf = out / "disp_plan.json"
    if pf.is_file():
        try:
            pl = json.loads(pf.read_text(encoding="utf-8"))
        except Exception:                            # noqa: BLE001
            return False
        if str(pl.get("method")) != "phonopy":
            return False
        if abs(float(pl.get("fd_distance", -1)) - float(conf["FD_DISTANCE"])) > 1e-12:
            return False
    return True


def _quarantine_stale_dataset(out):
    """把【与当前输入/配置不一致】的旧数据集非破坏地移到时间戳子目录，迫使重生成。

    目的目录每次用唯一时间戳：旧版一直往同一目录 move，与上次同名目录撞名会让
    shutil.move 把源搬进目标里面（嵌套）或留下半搬状态。
    """
    _base = out / "symmetry_audit" / "pre_symmetry_dataset"
    dst = _base / (time.strftime("%Y%m%d-%H%M%S")
                   + "-%03d" % (int(time.time() * 1000) % 1000))
    try:
        dst.mkdir(parents=True, exist_ok=False)
    except OSError as e:
        print("[WARN] 无法建 %s（%s），跳过隔离" % (dst, e))
        return False
    moved = False
    names = ["phonopy_disp.yaml", "phono3py_disp.yaml", "SPOSCAR", "disp_plan.json"]
    names += sorted(p.name for p in out.glob("POSCAR-*"))
    names += sorted(p.name for p in out.glob("disp-*"))
    for n in names:
        src = out / n
        if not src.exists():
            continue
        try:
            shutil.move(str(src), str(dst / n))
            moved = True
        except OSError as e:
            print("[WARN] 移动旧数据集 %s 失败：%s" % (n, e))
    if moved:
        print("[..] 旧位移数据集已非破坏移到 %s" % dst)
    return moved


def _kpoints_mesh_line(path):
    """读 KPOINTS 的自动网格行（第 4 行）；非自动网格 / 读失败返回 None。"""
    try:
        ln = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
        if len(ln) >= 4 and ln[1].strip().startswith("0"):
            return " ".join(ln[3].split()[:3])
    except OSError:
        return None
    return None


def check_frame_consistency(out, tol=1e-3):
    """帧间一致性硬闸：已有 disp-* 必须【同胞同网格】。

    retry 对 fanout 只补缺失子目录；计划中途变过（SUPERCELL/MIN_SC_LEN）或同一
    KSPACING 落在 VASPKIT 取整边界时，同一 step2_disp 里会留下两种超胞/网格，
    S3 把两批当同一批数据拟合 → fc2 报废且一声不响。扇出前逐个已有帧比对：
      ① POSCAR 超胞晶格（基准 SPOSCAR；读不到退化为多数派晶格）
      ② KPOINTS 网格行（全体必须唯一）
    """
    entries = [(d.name, _poscar_cell(d / "POSCAR"))
               for d in sorted(out.glob("disp-*")) if d.is_dir()]
    ref, ref_src = _poscar_cell(out / "SPOSCAR"), "SPOSCAR"
    if ref is None and entries:
        _votes = {}
        for _n, _c in entries:
            if _c is None:
                continue
            _key = tuple(round(x, 4) for _row in _c for x in _row)
            _votes.setdefault(_key, [0, _c])
            _votes[_key][0] += 1
        if _votes:
            _best = max(_votes.values(), key=lambda kv: kv[0])
            ref, ref_src = _best[1], "多数派帧（%d/%d）" % (_best[0], len(entries))
    if ref is None:
        print("[WARN] 帧间一致性闸跳过：既读不到 %s/SPOSCAR，也读不到任何 disp-*/POSCAR"
              % out.name)
        return
    bad_cell, mesh_groups, mesh_ref, n_seen = [], {}, None, 0
    for d in sorted(out.glob("disp-*")):
        if not d.is_dir():
            continue
        n_seen += 1
        cell = _poscar_cell(d / "POSCAR")
        if cell is None:
            bad_cell.append((d.name, "读不到 POSCAR 晶格"))
        else:
            dmax = max(abs(cell[i][j] - ref[i][j]) for i in range(3) for j in range(3))
            if dmax > tol:
                bad_cell.append((d.name, "与基准(%s)差 %.4f A（> %.0e）"
                                 % (ref_src, dmax, tol)))
        mesh = _kpoints_mesh_line(d / "KPOINTS")
        if mesh is None:
            continue
        mesh_groups.setdefault(mesh, []).append(d.name)
        if mesh_ref is None:
            mesh_ref = mesh
    if bad_cell:
        print("[ERROR] 帧间一致性闸：以下已有帧的超胞与当前计划不一致（共 %d 帧）："
              % len(bad_cell))
        for _n, _w in bad_cell[:8]:
            print("          %s：%s" % (_n, _w))
        if len(bad_cell) > 8:
            print("          ... 另有 %d 帧" % (len(bad_cell) - 8))
        sys.exit(
            "[ERROR] 混合超胞 = 混合数据集：S3 会把不同胞的力当同一批数据拟合，fc2 报废。\n"
            "        成因：计划中途变过（SUPERCELL / MIN_SC_LEN），而 retry 只补缺失帧、"
            "不动已有帧。\n        处置：确认已归档后清空本步产物再重跑：\n"
            "          rm -rf %s/disp-* %s/POSCAR-* %s/phonopy_disp.yaml %s/SPOSCAR\n"
            "        或直接 autozt -tt phonon-dft-cpu -p <材料> -j S2_disp rerun（推倒重来）。\n"
            "        想先跳过本闸排查：conf --set params.FRAME_GATE=off（不推荐）。"
            % (out, out, out, out))
    if len(mesh_groups) > 1:
        print("[ERROR] 帧间一致性闸：已有帧的 KPOINTS 网格不唯一（共 %d 种）："
              % len(mesh_groups))
        for _m, _names in sorted(mesh_groups.items(), key=lambda kv: -len(kv[1])):
            print("          [%s] %d 帧，如 %s" % (_m, len(_names), ", ".join(_names[:4])))
        sys.exit(
            "[ERROR] 同一个 step2_disp 里出现多种 k 网格 = 各帧的 k 点系统误差不同，\n"
            "        这个差会原样进 fc2（2026-09 那批虚频的同一类成因）。\n"
            "        成因：同一个 KSPACING 落在 VASPKIT 取整边界上，不同帧落到不同的整数网格。\n"
            "        处置：清空本步产物 → retry / rerun，并把 KSPACING 调离取整边界。")
    print("[OK] 帧间一致性闸：%d 个已有帧同胞同网格（mesh=[%s]）" % (n_seen, mesh_ref))


def build_displacements(out, reps, conf):
    """幂等：已有位移且与当前输入/配置一致就跳过；否则先过闸再落盘。"""
    check_existing_matches_input(out)
    if (out / "phonopy_disp.yaml").is_file() and glob.glob(str(out / "POSCAR-*")):
        if _dataset_matches(out, reps, conf):
            print("[..] 已有位移超胞，跳过生成（幂等）")
            return
        if _quarantine_stale_dataset(out):
            print("[WARN] 已有位移数据集与当前输入/配置不一致（单胞/超胞/位移幅度变了）："
                  "已非破坏移走，本步按当前配置重新生成")
    try:
        from phonopy import Phonopy
        from phonopy.interface.vasp import read_vasp, write_vasp
    except Exception as e:
        sys.exit("[ERROR] 无法 import phonopy（%s）—— S2 需要 phonopy 环境" % e)
    cell = read_vasp(str(out / "POSCAR"))
    ph = Phonopy(cell, supercell_matrix=np.array(kc.sc_matrix(reps), dtype=int),
                 primitive_matrix="auto")
    ph.generate_displacements(distance=float(conf["FD_DISTANCE"]))
    n = len(ph.supercells_with_displacements)
    gate(n, conf)
    write_vasp(str(out / "SPOSCAR"), ph.supercell, direct=True)
    for i, s in enumerate(ph.supercells_with_displacements, 1):
        write_vasp(str(out / ("POSCAR-%05d" % i)), s, direct=True)
    ph.save(filename=str(out / "phonopy_disp.yaml"))
    (out / "disp_plan.json").write_text(json.dumps(
        {"method": "phonopy", "fd_distance": float(conf["FD_DISTANCE"]),
         "supercell": kc.dim_str(reps), "n_disp": n},
        ensure_ascii=False, indent=2), encoding="utf-8")
    print("[OK] phonopy 2 阶位移：%d 帧（FD_DISTANCE=%s）" % (n, conf["FD_DISTANCE"]))


def gate(n, conf):
    """位移帧数硬闸。超过 MAX_DISP 直接停步，绝不落盘。"""
    cap = int(conf["MAX_DISP"])
    print("[..] 计划位移帧数 = %d（上限 MAX_DISP=%d）" % (n, cap))
    if n <= cap:
        return
    sys.exit("[ERROR] 位移帧数 %d > MAX_DISP=%d，step2 已停止（未写任何 disp-*）。\n"
             "        2 阶有限位移帧数由空间群约化决定；请扩/缩超胞或调 FD_DISTANCE：\n"
             "          conf --set params.SUPERCELL=\"2 2 2\"（或 params.MIN_SC_LEN=10.0）"
             % (n, cap))


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)

    prev = kc.find_prev_dir(cwd, PREV)
    if prev is None:
        sys.exit("[ERROR] 找不到 step1_std_opt 的结构")
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR", "step1_std_opt")
    # 结构体检第五条：S1 弛豫后的数值微畸变审计 + 可选对称化（默认检出即对称化）
    kc.symmetry_gate(out / "POSCAR", conf)

    meth = kc.read_method(prev / kc.METHOD_FILE)
    dim = (meth.get("DIM", "").lower() or kc.resolve_dim(out / "POSCAR")[0])
    _, vac_axis = kc.resolve_dim(out / "POSCAR", dim)
    require_dim(dim, ("2d", "3d"), "step2_disp", why="声子谱需要周期性边界")
    func = (conf["FUNC"] if conf["FUNC"] not in (None, "", "auto")
            else meth.get("FUNC", "pbe-d3").lower())
    # 泛函一致性硬检查：取力泛函必须与 step1 workflow_method.txt 记录的一致
    _rec_func = str(meth.get("FUNC") or "").strip().lower()
    if _rec_func and func != _rec_func:
        sys.exit("[ERROR] S2 的泛函与 step1 不一致：本步 FUNC=%s，workflow_method.txt 记 "
                 "FUNC=%s。\n        两者必须相同（同一势能面）；拒绝生成。\n"
                 "        处置：本材料 S2 的 FUNC 置 auto（推荐，自动继承），或先核对 S1 的泛函。"
                 % (func, _rec_func))

    if conf["SUPERCELL"]:
        reps = kc.parse_reps(conf["SUPERCELL"], dim,
                             vac_axis if vac_axis is not None else 2)
    else:
        reps = kc.supercell_matrix(out / "POSCAR", dim, conf["MIN_SC_LEN"],
                                   conf["MAX_MULTIPLE"],
                                   vac_axis if vac_axis is not None else 2)
    # 超胞保对称闸（2026-10-07 Mn2In2Se5）：超胞不保晶体全部点群操作时，拟合出的
    #   力常数必破对称（假虚频、kappa 随网格乱跳），在生成位移、花 DFT 机时之前拦下
    try:
        import symmetry_audit as _SA
    except ImportError as _e:
        print("[WARN] symmetry_audit.py 不可用，超胞保对称检查跳过：%s" % _e)
    else:
        _SA.supercell_symmetry_gate(out / "POSCAR", reps, conf["SUPERCELL_SYMMETRY"],
                                    dim, vac_axis if vac_axis is not None else 2,
                                    conf["MAX_MULTIPLE"])
    mesh = kc.mesh_str(conf["PHONON_MESH"].split(), dim,
                       vac_axis if vac_axis is not None else 2)
    print("[..] 维度=%s 超胞=%s mesh=%s FUNC=%s" % (dim.upper(), kc.dim_str(reps), mesh, func))
    kc.write_kl_params(out / kc.KL_PARAMS, DIM=dim.upper(), SUPERCELL=kc.dim_str(reps),
                       MESH=mesh, FUNC=func, FD_DISTANCE=conf["FD_DISTANCE"])

    # 2D 几何硬闸 + 残余面内应力（Huang 零应力条件）：生成位移之前拦
    if dim == "2d":
        _geo = kc.check_2d_vacuum(out / "POSCAR", vac_axis if vac_axis is not None else 2)
        _s1 = sorted(cwd.glob("step1*/OUTCAR"))
        if _s1:
            kc.check_2d_stress(_s1[-1], _geo["h_perp_A"], _geo["thickness_d_A"],
                               conf["STRESS_2D_THR"])
        else:
            print("[WARN] 找不到 S1 的 OUTCAR（%s/step1*/OUTCAR），2D 残余应力未核验" % cwd)

    build_displacements(out, reps, conf)

    poscars = sorted(glob.glob(str(out / "POSCAR-*")))
    if not poscars:
        sys.exit("[ERROR] 没有产出 POSCAR-* 位移超胞")
    # 二道闸：目录里本来就躺着旧版本留下的多余帧也要拦住，别扇出
    if len(poscars) > int(conf["MAX_DISP"]):
        sys.exit("[ERROR] %s 下已有 %d 个 POSCAR-*，超过 MAX_DISP=%d —— 拒绝扇出。\n"
                 "        多半是旧版本（无闸门）留下的，先清干净再重跑：\n"
                 "        rm -rf %s && autozt -tt phonon-dft-cpu -p <材料> -j S2_disp rerun"
                 % (OUTDIR, len(poscars), int(conf["MAX_DISP"]), OUTDIR))

    # 帧间一致性硬闸：扇出之前，已有 disp-* 必须同胞同网格（SPOSCAR 此时已存在）
    if str(conf.get("FRAME_GATE", "on")).strip().lower() not in ("off", "0", "false", "no"):
        check_frame_consistency(out)

    here = Path(__file__).resolve().parent
    incar_tpl = kc.resolve_submit(here, dim, "incar_force")
    submit_tpl = kc.resolve_submit(here, dim, "submit_std")
    # 偶极修正必须与 S1 一致 —— 直接继承其 INCAR，不复算
    _dip_lines, _dip_src = kc.dipole_line_from_incar(cwd / "step1_std_opt" / "INCAR")
    print("[..] 偶极修正：%s" % (_dip_src or "上游 INCAR 没找到，未施加"))
    encut = None
    for p in poscars:
        num = os.path.basename(p).split("-", 1)[1]
        d = out / ("disp-%s" % num)
        d.mkdir(exist_ok=True)
        if (d / "INCAR").is_file() and (d / "POSCAR").is_file():
            _outcar = d / "OUTCAR"
            _done = False
            if _outcar.is_file():
                try:
                    _done = "General timing and accounting informations" in \
                        _outcar.read_text(errors="ignore")[-400000:]
                except OSError:
                    _done = False
            if not _done:
                _ic = stepconf.apply_incar_file(d / "INCAR")
                if _ic:
                    print("[..] %s：待重算帧，按 step.conf [incar] 更新 %d 项"
                          % (d.name, len(_ic)))
            continue
        shutil.copyfile(p, d / "POSCAR")
        kc.vaspkit_kpoints(d, conf["KSCHEME"], conf["KSPACING"],
                           conf["VASPKIT_EXE"], dim, vac_axis)
        kc.vaspkit_potcar(d, conf["VASPKIT_EXE"])
        if encut is None:
            encut = conf["ENCUT"] or kc.encut_from_potcar(d / "POTCAR", conf["ENCUT_FACTOR"])
        kc.render_tpl(incar_tpl,
                      {"SYSTEM": "%s force %s" % (cwd.name, num),
                       "ENCUT": encut,
                       "GGA": kc.GGA_MAP.get(func, "PS"),
                       "FORCE_PREC": conf["FORCE_PREC"],
                       "FORCE_EDIFF": conf["FORCE_EDIFF"],
                       "FORCE_LREAL": conf["FORCE_LREAL"],
                       "DIPOLE_LINE": "\n".join(_dip_lines),
                       "VDW_LINE": ("IVDW = %s" % kc.VDW_MAP[func])
                                   if kc.VDW_MAP.get(func) else "# no vdW"},
                      d / "INCAR")
        # 项目级 incar_force_*.tpl 若遮蔽技能模板、缺新占位符 → 渲染会静默跳过，强制落一遍
        kc.enforce_incar_tags(d / "INCAR",
                              {"PREC": conf["FORCE_PREC"],
                               "EDIFF": conf["FORCE_EDIFF"],
                               "LREAL": conf["FORCE_LREAL"],
                               "ADDGRID": ".TRUE."},
                              label="step2 取力 %s：" % num)
        _ic = stepconf.apply_incar_file(d / "INCAR")
        if _ic:
            print("[..] %s：step.conf [incar] 覆盖 %d 项" % (d.name, len(_ic)))
        kc.write_submit(submit_tpl, d / "submit.sh",
                        {"JOBNAME": "%s-phonon-dft-cpu-S2-%s" % (cwd.name, num)})
        stepconf.apply_submit(d / "submit.sh", conf.submit)
    print("[DONE] %s：%d 个位移子目录就绪（fanout disp-*）" % (OUTDIR, len(poscars)))


if __name__ == "__main__":
    main()
