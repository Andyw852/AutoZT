#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kl_fc_backends.py —— S5_fc 计算节点作业驱动：数据桥接 → 拟合 → 双格式导出 → 虚频闸。

设计：S5_fc 现在是【提交到计算节点】的作业（不再登录节点 gen 里裸跑）。submit.sh 依次调用
本文件的子命令，全部在 step5_fc/ 工作目录内执行：

  prep          收 step4_disp/disp-*/vasprun.xml → FORCES_FC3 → 从 phono3py_disp.yaml 建
                phono3py 对象；写 POSCAR(原胞)/SPOSCAR(超胞)；alm 分支再落 dataset_disps.npy
                / dataset_forces.npy（末帧=零平衡帧，供 pheasy 用）；拷 phono3py_disp.yaml
                到 phono3py/ 子目录，并（若有 step3_nac）生成 BORN。
  fit_phono3py  phono3py + symfc/alm 拟合 fc2/fc3（full）→ phono3py/{fc2.hdf5,fc3.hdf5,
                phono3py_params.yaml}。移植自 run_phono3py_fit.sh（去掉 npy 读入与
                VASP↔phonopy 原子重排：taskflow 的位移本就是 phonopy 超胞原子序）。
  collect_pheasy 把 pheasy CLI 产出的 fc2.hdf5/fc3.hdf5（在 cwd）搬进 phono3py/。
  post          可选 shengbte 导出（fc2/fc3 → hiphive → shengbte/FORCE_CONSTANTS_2ND/3RD）
                + 声子谱虚频闸（phonopy mesh 最小频率）→ step5_fc/phonon_summary.json。

拟合器/求解器由 step.conf 决定，gen_step5_fc.py 把解析结果写进 fit_config.json，本文件只读它。
产出布局（满足"任一拟合器都产出两套格式，放不同文件夹"）：
  step5_fc/phono3py/  fc2.hdf5 fc3.hdf5 phono3py_disp.yaml [phono3py_params.yaml] [BORN]
  step5_fc/shengbte/  FORCE_CONSTANTS_2ND FORCE_CONSTANTS_3RD POSCAR   (CONTROL 由 S6 按运行参数写)
  step5_fc/phonon_summary.json   ← S5 marker: '"stable": true'
"""
import glob
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

DISP_DIR   = "../step4_disp"
NAC_DIR    = "../step3_nac"
P3PY_SUB   = "phono3py"
SB_SUB     = "shengbte"

# --- ZA 弯曲支二次性判据参数 -------------------------------------------------
# ★ TODO（2026-09-20）：ZA 判据现在至少有【三份】各自独立的副本：
#     1) skill/fc-fit/fc_plot_phonon.py       (_inplane_qdirs/_k_sum_sign/
#        _vacuum_axis_in_primitive/_eig_out_of_plane/za_power_law_eig/_za_summary)
#     2) 本文件 skill/kl-dft-cpu/kl_fc_backends.py（已与 1 对齐）
#     3) skill/_common/mlff/klmlff_common.py  (inplane_qdirs/za_check_2d，
#        仍是无条件 e_i+e_j、无 vac 轴映射、无本征矢量判据；被
#        phonon-mlff-cpu/gpu 的 phonon_fit_driver 调用)
#   这份重复已经直接造成过一次"只修了一侧"的测量 bug（本次只修了 2）。
#   长期应把这段公共逻辑挪到 skill/_common/（例如 _common/za_common.py），
#   三处 import 同一实现，否则下次修 ZA 判据仍要修三遍、且很容易漏边。
from za_2d import (ZA_P_RANGE, ZA_ZFRAC_MIN, ZA_MARGIN_MIN, ZA_R2_MIN,
                   ZA_SYMPREC_DEFAULT, ZA_SYMPREC_RELAX,
                   za_power_law, za_power_law_eig,
                   k_sum_sign as _k_sum_sign,
                   inplane_qdirs as _inplane_qdirs,
                   vacuum_axis_in_primitive as _vacuum_axis_in_primitive,
                   eig_out_of_plane as _eig_out_of_plane,
                   cell_symmetry as _cell_symmetry,
                   relaxed_symprec as _relaxed_symprec,
                   clone_phonopy as _clone_phonopy)
# 虚频判据唯一真源（2026-09-22）：S5 是唯一裁判，其余位置只读它的结论。
import imag_policy


# ==========================================================================
# 小工具
# ==========================================================================
def load_cfg(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def run(cmd, **kw):
    """前台跑命令，非零即抛。cmd 为字符串走 shell。"""
    print("[cmd] %s" % (cmd if isinstance(cmd, str) else " ".join(cmd)), flush=True)
    r = subprocess.run(cmd, shell=isinstance(cmd, str), **kw)
    if r.returncode != 0:
        sys.exit("[ERROR] 命令失败(rc=%d): %s" % (r.returncode, cmd))
    return r


def collect_vaspruns(disp_dir):
    """按编号排序收 disp-*/vasprun.xml。

    返回 (位移帧路径列表, 缺帧编号列表, 平衡帧路径或 None)。
    kleq：disp-00000 是平衡帧（未位移的完美超胞），单独摘出 —— 它不属于位移
    dataset，只用来扣残余力；混进位移列表会让帧序整体错位一格。
    """
    subs = sorted(glob.glob(str(Path(disp_dir) / "disp-*")),
                  key=lambda p: int(re.search(r"disp-(\d+)", p).group(1)))
    files, missing, eq = [], [], None
    for d in subs:
        num = re.search(r"disp-(\d+)", d).group(1)
        vr = Path(d) / "vasprun.xml"
        ok = vr.is_file() and vr.stat().st_size
        path = "%s/disp-%s/vasprun.xml" % (disp_dir, num)
        if int(num) == 0:
            eq = path if ok else None
        elif ok:
            files.append(path)
        else:
            missing.append(num)
    return files, missing, eq


def make_born(out):
    """有 step3_nac 的 vasprun 就用 phonopy-vasp-born 生成 out/BORN。返回是否成功。"""
    vr = Path(NAC_DIR) / "vasprun.xml"
    if not vr.is_file():
        print("[..] 无 step3_nac 结果 → 不加 NAC")
        return False
    run("phonopy-vasp-born %s > %s 2>born.log || true" % (str(vr), str(out / "BORN")))
    ok = (out / "BORN").is_file() and (out / "BORN").stat().st_size > 0
    print("[%s] BORN（NAC）%s" % ("OK" if ok else "WARN",
                                  "已生成" if ok else "生成失败，退回无 NAC"))
    if not ok:
        try:
            (out / "BORN").unlink()
        except OSError:
            pass
    return ok


# ==========================================================================
# prep：数据桥接（vasprun + phono3py_disp.yaml → ph 对象 / POSCAR / SPOSCAR / npy）
# ==========================================================================
def _read_forces(files, natom):
    """读一批 vasprun.xml 的力，返回 (n, natom, 3) 的 ndarray。"""
    import numpy as np
    from phonopy.interface.vasp import parse_set_of_forces
    d = parse_set_of_forces(natom, [str(f) for f in files], verbose=False)
    return np.asarray(d["forces"], float)


def _load_ph3_with_forces(disp_yaml):
    """load phono3py_disp.yaml（含位移 dataset）；phono3py.load 见 cwd 里有 FORCES_FC3
    会自动把力读进 dataset（日志里 'Displacement dataset for fc3 was read from FORCES_FC3'）。
    注意：不能用 parse_FORCES_FC3()——那是 type1(有限位移) 专用，要 dataset['natom']/
    'first_atoms'；alm 是 type2(随机位移) dataset，会 KeyError。"""
    import phono3py
    # kleq：prep 现在把（可能缺帧、已扣残余力的）dataset 存成 phono3py_params.yaml，
    #   它自带 forces，优先读它；没有才回落到 disp yaml + FORCES_FC3 的老路。
    _pm = Path(disp_yaml).parent / "phono3py_params.yaml"
    if _pm.is_file():
        _p = phono3py.load(str(_pm), produce_fc=False, is_nac=False, log_level=0)
        _d = _p.dataset or {}
        if (_d.get("forces") is not None
                or (_d.get("first_atoms") and "forces" in _d["first_atoms"][0])):
            print("[..] 力来自 %s（prep 已扣平衡帧残余力）" % _pm.name)
            return _p
    ph3 = phono3py.load(str(disp_yaml), produce_fc=False, is_nac=False, log_level=1)
    ds = ph3.dataset
    has_forces = (("forces" in ds and ds["forces"] is not None)
                  or ("first_atoms" in ds and ds.get("first_atoms")
                      and "forces" in ds["first_atoms"][0]))
    if not has_forces:
        # 兜底：个别版本 load 不自动读力时，显式指定 FORCES_FC3 重载一次（比手撕文本稳）。
        ph3 = phono3py.load(str(disp_yaml), forces_fc3_filename="FORCES_FC3",
                            produce_fc=False, is_nac=False, log_level=1)
        ds = ph3.dataset
        has_forces = "forces" in ds and ds["forces"] is not None
        print("[..] load 未自动带力，已用 forces_fc3_filename 重载 FORCES_FC3")
    if not has_forces:
        sys.exit("[ERROR] dataset 里没有 forces —— FORCES_FC3 未被读入。"
                 "确认 step5_fc/ 下已生成 FORCES_FC3（prep 的 --cf3 步）。")
    return ph3


def cmd_prep(cfg):
    from phonopy.interface.vasp import write_vasp
    import numpy as np

    out = Path.cwd()
    method = cfg["METHOD"]
    engine = cfg["FIT_ENGINE"]
    p3dir = out / P3PY_SUB
    p3dir.mkdir(exist_ok=True)

    disp_yaml = Path(DISP_DIR) / "phono3py_disp.yaml"
    if not disp_yaml.is_file():
        sys.exit("[ERROR] %s 不存在（step4 未生成位移）" % disp_yaml)
    shutil.copyfile(disp_yaml, out / "phono3py_disp.yaml")
    shutil.copyfile(disp_yaml, p3dir / "phono3py_disp.yaml")

    import numpy as np
    import phono3py

    files, missing, eq_file = collect_vaspruns(DISP_DIR)
    if not files:
        sys.exit("[ERROR] %s 下没有 disp-*/vasprun.xml，位移单点还没算完" % DISP_DIR)
    if not eq_file:
        sys.exit("[ERROR] 缺平衡帧 %s/disp-00000/vasprun.xml。\n"
                 "        它是未位移的完美超胞单点：拟合要用它扣残余力，也是缺帧容错的锚。\n"
                 "        补法：tf -tt kl-dft-cpu -p <材料> -j S4_disp retry（gen 幂等，只补出缺的\n"
                 "        disp-00000，已算完的位移帧不动），再 tf ... -j S4_disp start -f。"
                 % DISP_DIR)

    ph3 = phono3py.load(str(out / "phono3py_disp.yaml"),
                        produce_fc=False, is_nac=False, log_level=0)
    nsc = len(ph3.supercell)
    ds = ph3.dataset or {}
    # type-2（随机位移）dataset 有 displacements 数组；type-1（有限位移）是 first_atoms。
    is_rand = ds.get("displacements") is not None

    # kleq：平衡帧的残余力（未完全弛豫 / egg-box 都体现在这里），逐帧扣掉。
    f_eq = _read_forces([eq_file], nsc)[0]
    print("[..] 平衡帧 max|F_eq| = %.4f eV/Å（逐帧扣除）"
          % float(np.linalg.norm(f_eq, axis=1).max()))

    n_all = len(files) + len(missing)
    if missing:
        head = ",".join(missing[:8]) + ("…" if len(missing) > 8 else "")
        if not is_rand:
            sys.exit("[ERROR] 有限位移法（METHOD=findiff）不能缺帧：每个对称不等价位移都要\n"
                     "        参与有限差分重建 fc2/fc3，少一帧就重建不出来。缺 %d 帧（disp-%s）。\n"
                     "        请补算后再来 S5；要容忍失败帧就把 step4 改成 METHOD=alm\n"
                     "        （随机位移 + 回归拟合，天然可以少几帧）。" % (len(missing), head))
        ratio = len(files) / float(max(n_all, 1))
        min_ratio = float(cfg.get("MIN_SUCCESS_RATIO", 0.9) or 0)
        min_frames = int(cfg.get("MIN_SUCCESS_FRAMES", 0) or 0)
        if ratio < min_ratio or len(files) < min_frames:
            sys.exit("[ERROR] 成功帧 %d/%d（%.0f%%）低于下限"
                     "（MIN_SUCCESS_RATIO=%s，MIN_SUCCESS_FRAMES=%s）。\n"
                     "        补算失败帧，或在 step5_fc 的 step.conf 里调低下限。"
                     % (len(files), n_all, ratio * 100, min_ratio, min_frames))
        print("[WARN] 缺 %d 帧（disp-%s）—— 随机位移法按 %d/%d（%.0f%%）继续拟合"
              % (len(missing), head, len(files), n_all, ratio * 100))
    print("[..] 收 %d 个 vasprun.xml（方法=%s，拟合器=%s）" % (len(files), method, engine))

    forces = _read_forces(files, nsc) - f_eq[None]
    if is_rand:
        disps = np.asarray(ds["displacements"], float)
        idx = [int(re.search(r"disp-(\d+)", f).group(1)) - 1 for f in files]
        if max(idx) >= len(disps):
            sys.exit("[ERROR] disp-%05d 超出 phono3py_disp.yaml 里的 %d 帧位移 —— "
                     "step4 的位移集与 disp-* 目录对不上，请 rerun S4_disp。"
                     % (max(idx) + 1, len(disps)))
        ph3.dataset = {"displacements": disps[idx], "forces": forces}
    else:
        # phono3py 4.x：findiff dataset 保留 included=False 的 second_atoms
        # （被 FC3_CUTOFF_PAIR 过滤掉的对位移，未生成 disp 目录、未算力），而
        # forces setter 遍历所有 second_atoms（不跳过 included=False）——需提供
        # 与 dataset 位移数一致的完整数组，included=False 填 0（fc3 拟合会跳过）。
        # 映射：disp 目录编号 == dataset 位移的 id（first_atoms[0].id=1，
        # second_atoms 的 id=2..N），故 forces_full[id-1] = 对应帧的力。
        n_total = 1 + sum(len(fa.get("second_atoms", []))
                          for fa in ph3.dataset.get("first_atoms", []))
        forces_full = np.zeros((n_total, nsc, 3), dtype=float)
        for f, fv in zip(files, forces):
            num = int(re.search(r"disp-(\d+)", f).group(1))
            forces_full[num - 1] = fv
        ph3.forces = forces_full
    # 存成自带力的 phono3py_params.yaml：后面的拟合直接读它，不再走
    #   `phono3py --cf3` + FORCES_FC3（那条路要求文件数==位移数，缺一帧就错位）。
    ph3.save(str(out / "phono3py_params.yaml"))
    if is_rand:
        print("[OK] phono3py_params.yaml：%d 帧力已写入（已扣平衡帧）" % len(forces))
    else:
        print("[OK] phono3py_params.yaml：%d 帧力已写入（已扣平衡帧；另 %d 个 cutoff 外对填 0）"
              % (len(forces), n_total - len(forces)))

    ph3 = _load_ph3_with_forces(out / "phono3py_disp.yaml")
    write_vasp(str(out / "POSCAR"), ph3.unitcell, direct=True)     # 原胞
    write_vasp(str(out / "SPOSCAR"), ph3.supercell, direct=True)   # 理想超胞

    # pheasy 走随机位移路：落 dataset_disps.npy/dataset_forces.npy（末帧=零平衡帧）。
    #   taskflow 的位移已是"相对完美超胞的笛卡尔位移"，直接当 cartesian 喂 pheasy；
    #   追加一帧零位移/零力当"平衡参考"，pheasy 脚本会 (x - x[-1])[:-1] 还原成原样。
    if engine == "pheasy":
        if method != "alm":
            sys.exit("[ERROR] FIT_ENGINE=pheasy 需要随机位移（METHOD=alm）。"
                     "findiff 请用 FIT_ENGINE=phono3py。")
        ds = ph3.dataset
        if "displacements" not in ds or "forces" not in ds:
            sys.exit("[ERROR] 随机位移 dataset 缺 displacements/forces（type-2 期望）")
        disps = np.asarray(ds["displacements"], float)   # (n, Nsc, 3) 笛卡尔 Å
        forces = np.asarray(ds["forces"], float)         # (n, Nsc, 3) eV/Å
        nsc = disps.shape[1]
        zeros = np.zeros((1, nsc, 3))
        np.save(out / "dataset_disps.npy", np.concatenate([disps, zeros], axis=0))
        np.save(out / "dataset_forces.npy", np.concatenate([forces, zeros], axis=0))
        print("[OK] dataset_disps/forces.npy 就绪：%d 帧(+1 平衡) × %d 原子" % (len(disps), nsc))

    make_born(p3dir)
    print("[DONE] prep 完成")


# ==========================================================================
# fit_phono3py：phono3py + symfc/alm → full fc2/fc3 → phono3py/
# ==========================================================================
def cmd_fit_phono3py(cfg):
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5

    out = Path.cwd()
    p3dir = out / P3PY_SUB
    calc = (cfg.get("FC_CALC") or "symfc").lower()
    cutoff = cfg.get("FC3_CUTOFF")
    ph3 = _load_ph3_with_forces(out / "phono3py_disp.yaml")

    print("[..] phono3py 拟合器=%s  fc3_cutoff=%s" % (calc, cutoff))
    # full fc（is_compact_fc=False）：既给 phono3py-load κ，也便于就地喂 hiphive 导 ShengBTE。
    ph3.produce_fc2(fc_calculator=calc, is_compact_fc=False)
    opts = None if cutoff in (None, "None", "none", "") else "cutoff = %s" % cutoff
    ph3.produce_fc3(fc_calculator=calc, fc_calculator_options=opts, is_compact_fc=False)

    write_fc2_to_hdf5(ph3.fc2, filename=str(p3dir / "fc2.hdf5"))   # full，无 p2s_map
    write_fc3_to_hdf5(ph3.fc3, filename=str(p3dir / "fc3.hdf5"))

    # NAC → phono3py_params.yaml（有 BORN 才写；供溯源/kappa 备选）
    born = p3dir / "BORN"
    if born.is_file():
        try:
            from phonopy.file_IO import parse_BORN
            nac = parse_BORN(ph3.primitive, filename=str(born))
            if isinstance(nac, dict) and not nac.get("factor"):
                nac["factor"] = 14.399652
            ph3.nac_params = nac
            print("[OK] BORN → nac_params")
        except Exception as e:
            print("[WARN] 读 BORN 失败，不写 nac_params：%s" % e)
    try:
        ph3.save(str(p3dir / "phono3py_params.yaml"))
    except Exception as e:
        print("[WARN] 写 phono3py_params.yaml 失败（不影响 kappa，用 phono3py_disp.yaml）：%s" % e)

    for f in ("fc2.hdf5", "fc3.hdf5"):
        if not (p3dir / f).is_file():
            sys.exit("[ERROR] phono3py 拟合未产出 %s" % f)

    # ShengBTE 导出：就地用内存里的 full 数组（和 run_phono3py_fit.sh 一致，最稳）
    if str(cfg.get("EXPORT_SHENGBTE", "true")).lower() in ("true", "1", "yes"):
        import numpy as np
        _export_shengbte(out, np.asarray(ph3.fc2), np.asarray(ph3.fc3))
    print("[DONE] fit_phono3py：phono3py/fc2.hdf5 + fc3.hdf5 就绪")


# ==========================================================================
# collect_pheasy：把 pheasy 的产物搬进 phono3py/
# ==========================================================================
def cmd_collect_pheasy(cfg):
    out = Path.cwd()
    p3dir = out / P3PY_SUB
    for f in ("fc2.hdf5", "fc3.hdf5"):
        src = out / f
        if not src.is_file():
            sys.exit("[ERROR] pheasy 没产出 %s（看 queue.err/out）" % f)
        shutil.move(str(src), str(p3dir / f))
    print("[..] collect_pheasy：fc2.hdf5 + fc3.hdf5 → phono3py/")

    # ShengBTE 导出：pheasy --full_ifc 产出 full fc，复读成数组喂 hiphive
    if str(cfg.get("EXPORT_SHENGBTE", "true")).lower() in ("true", "1", "yes"):
        try:
            fc2, fc3 = _read_full_fc(p3dir)
            _export_shengbte(out, fc2, fc3)
        except Exception as e:
            print("[WARN] pheasy fc 复读/导出 ShengBTE 失败（phono3py 路不受影响）：%s" % e)
    print("[DONE] collect_pheasy 完成")


# ==========================================================================
# post：shengbte 导出（可选） + 虚频闸 → phonon_summary.json
# ==========================================================================
def _read_full_fc(p3dir):
    """从 phono3py/{fc2.hdf5,fc3.hdf5} 复读 full fc 数组（供 hiphive）。"""
    import numpy as np
    from phono3py.file_IO import read_fc2_from_hdf5, read_fc3_from_hdf5
    fc2 = np.asarray(read_fc2_from_hdf5(filename=str(p3dir / "fc2.hdf5")))
    fc3 = np.asarray(read_fc3_from_hdf5(filename=str(p3dir / "fc3.hdf5")))
    return fc2, fc3


def _export_shengbte(out, fc2, fc3):
    """给定 full fc2/fc3 数组 → hiphive → shengbte/FORCE_CONSTANTS_2ND + _3RD + POSCAR。
    CONTROL 不在此写：它依赖 S6 的运行参数（T/mesh/scalebroad/NAC），由 S6 生成。"""
    import numpy as np
    from ase import Atoms
    from phonopy.interface.vasp import read_vasp
    try:
        from hiphive import ForceConstants
    except Exception as e:
        print("[WARN] 无 hiphive，跳过 ShengBTE 力常数导出（phono3py 侧不受影响）：%s" % e)
        return False

    sbdir = out / SB_SUB
    sbdir.mkdir(exist_ok=True)

    uc = read_vasp(str(out / "POSCAR"))
    sc = read_vasp(str(out / "SPOSCAR"))
    to_ase = lambda c: Atoms(numbers=c.numbers, scaled_positions=c.scaled_positions,
                             cell=c.cell, pbc=True)
    uc_ase, sc_ase = to_ase(uc), to_ase(sc)
    fc2 = np.asarray(fc2); fc3 = np.asarray(fc3)
    natom = len(sc_ase)

    try:
        if fc2.shape[0] != natom or fc3.shape[0] != natom:
            print("[WARN] fc 非 full（fc2%s fc3%s，超胞 %d 原子），ShengBTE 导出仅支持 full fc；"
                  "跳过。phono3py 路不受影响。" % (fc2.shape, fc3.shape, natom))
            return False
        fcs = ForceConstants.from_arrays(sc_ase, fc2_array=fc2, fc3_array=fc3)
    except Exception as e:
        print("[WARN] 构造 hiphive ForceConstants 失败，跳过 ShengBTE 导出：%s" % e)
        return False

    # hiphive write_to_shengBTE 内部要求原胞分数坐标严格 <1，抹掉 0.9999.../1e-10 噪声
    def _wrap(a, eps=1e-9):
        fr = a.get_scaled_positions(wrap=True)
        fr = np.where(fr >= 1.0 - eps, 0.0, fr)
        fr = np.where(fr < eps, 0.0, fr)
        a.set_scaled_positions(fr)
    _wrap(uc_ase)
    try:
        _wrap(fcs._supercell)
    except Exception:
        pass

    fcs.write_to_phonopy(str(sbdir / "FORCE_CONSTANTS_2ND"), format="text")
    fcs.write_to_shengBTE(str(sbdir / "FORCE_CONSTANTS_3RD"), uc_ase)
    from ase.io import write as ase_write
    ase_write(str(sbdir / "POSCAR"), uc_ase, format="vasp", direct=True, sort=False)
    ok = (sbdir / "FORCE_CONSTANTS_2ND").is_file() and (sbdir / "FORCE_CONSTANTS_3RD").is_file()
    print("[%s] ShengBTE 力常数导出%s（CONTROL 由 S6 按运行参数写）"
          % ("OK" if ok else "WARN", "完成" if ok else "失败"))
    return ok


def _parse_min_freq(band_yaml):
    if not Path(band_yaml).is_file():
        return None
    fr = []
    for ln in Path(band_yaml).read_text(errors="ignore").splitlines():
        m = re.match(r"\s*frequency:\s*([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)\s*$", ln)
        if not m:
            continue
        try:
            fr.append(float(m.group(1).replace("d", "e").replace("D", "e")))
        except ValueError:
            continue
    return min(fr) if fr else None


def _za_check(cfg, ph, cart_vac_axis, is2d):
    """P2-2：2D 的 ZA 二次性闸门。返回写进 phonon_summary.json 的 dict（3D → None）。

    ★ 与 skill/fc-fit/fc_plot_phonon.py 的 _za_summary 同一套判据（同一 bug 的两份
    副本，见文件顶部 TODO）：
      * vac 轴映射：Cartesian 真空轴先映射成 原胞基矢 下标（_vacuum_axis_in_primitive），
        否则生产原胞（真空在第 0 基矢）会量到"真空方向 + Γ-M"；
      * 面内第二方向：非正交原胞用 _k_sum_sign 取 e_i + s*e_j，否则 60° 生产原胞给出
        两个 Γ-M；
      * ZA 支识别：本征矢量面外占比 >= ZA_ZFRAC_MIN 的"最低频支"（旧代码取全局最低支）；
      * 质量字段：za_r2_* / za_log_resid_rms_* / za_margin_* / za_needs_review / 支号 /
        面外占比 / 混合说明；旧最低支结果并列保留但不作为判据；
      * 对称性审计：symprec 1e-5 与 1e-4 空间群不一致时，用 1e-4 的克隆做 ZA 拟合，
        默认构造的 p 记录在 za_exponent_default_qX 供对照。
    """
    mode = str(cfg.get("ZA_CHECK", "auto")).strip().lower()
    if mode in ("off", "false", "0", "no") or (not is2d and mode not in ("on", "true", "1", "yes")):
        return None
    p_lo, p_hi = ZA_P_RANGE
    qmax = float(cfg.get("ZA_QMAX", 0.05))
    res = {
        "mode": mode, "qmax": qmax, "p_range": [p_lo, p_hi], "is_2d": bool(is2d),
        "vacuum_axis_cartesian": (None if cart_vac_axis is None else int(cart_vac_axis)),
        "vacuum_axis_in_primitive": None,
        "dirs": [], "p": [], "min_freq": [], "p_lowest": [],
        "za_r2_q1": None, "za_r2_q2": None,
        "za_log_resid_rms_q1": None, "za_log_resid_rms_q2": None,
        "za_margin_q1": None, "za_margin_q2": None,
        "za_margin_min": ZA_MARGIN_MIN, "za_r2_min": ZA_R2_MIN,
        "za_needs_review": None,
        "za_branch_index_q1": None, "za_branch_index_q2": None,
        "za_zfrac_q1": None, "za_zfrac_q2": None,
        "za_zfrac_min_q1": None, "za_zfrac_min_q2": None,
        "za_branch_note_q1": "", "za_branch_note_q2": "",
        "za_lowest_exponent_q1": None, "za_lowest_exponent_q2": None,
        "za_lowest_r2_q1": None, "za_lowest_r2_q2": None,
        "za_lowest_log_resid_rms_q1": None, "za_lowest_log_resid_rms_q2": None,
        "za_symprec": ZA_SYMPREC_DEFAULT, "za_symprec_relax": ZA_SYMPREC_RELAX,
        "za_spacegroup_default": None, "za_spacegroup_relax": None,
        "za_symmetry_relaxed": None,
        "za_exponent_default_q1": None, "za_exponent_default_q2": None,
        "ok": False, "note": "",
    }
    # ★ 方向 bug 修复：Cartesian 真空轴 → 原胞基矢下标
    try:
        ax = _vacuum_axis_in_primitive(ph, cart_vac_axis)
    except Exception as e:
        res["error"] = "真空轴 → 原胞基矢映射失败：%s" % e
        return res
    res["vacuum_axis_in_primitive"] = int(ax)
    try:
        qdirs = _inplane_qdirs(ph, ax)
    except Exception as e:
        res["error"] = "面内方向选择失败：%s" % e
        return res
    res["dirs"] = [list(d) for d in qdirs]

    # --- 对称性审计（数值微畸变使默认 symprec 低估空间群 → ZA 假线性）---
    notes = []
    res["za_spacegroup_default"] = _cell_symmetry(ph, ZA_SYMPREC_DEFAULT)
    res["za_spacegroup_relax"] = _cell_symmetry(ph, ZA_SYMPREC_RELAX)
    relaxed = _relaxed_symprec(ph, ZA_SYMPREC_RELAX)
    ph_za = ph
    if relaxed is not None:
        try:
            ph_za = _clone_phonopy(ph, relaxed)
            res["za_symprec"] = relaxed
        except Exception as e:
            relaxed = None
            notes.append("无法按 symprec=%.0e 重建 phonopy（%s）" % (ZA_SYMPREC_RELAX, e))
    res["za_symmetry_relaxed"] = bool(ph_za is not ph)
    if ph_za is not ph:
        notes.append(
            "phonopy 默认 symprec=%.0e 看到 %s，而 symprec=%.0e 看到 %s：原胞只是数值上"
            "微畸变，ZA 拟合改用放宽容差（默认构造的 p 见 za_exponent_default_q1/q2）"
            % (ZA_SYMPREC_DEFAULT, res["za_spacegroup_default"],
               ZA_SYMPREC_RELAX, res["za_spacegroup_relax"]))

    had_nac = getattr(ph, "nac_params", None)
    tags = ("q1", "q2")
    ps, wmins, fits = [], [], []
    try:
        ph.nac_params = None
        try:
            ph_za.nac_params = None
        except Exception:
            pass
        for i, d in enumerate(qdirs):
            eig = None
            try:
                eig = za_power_law_eig(ph_za, qdir=d, qmax=qmax, n=12,
                                       cart_vac_axis=cart_vac_axis)
            except Exception as e:
                notes.append("qdir %s 本征矢量 ZA 拟合失败：%s" % (list(d), e))
            try:
                lp, lw, lf = za_power_law(ph_za, qdir=d, qmax=qmax, n=12)
            except Exception as e:
                lp, lw, lf = None, None, None
                notes.append("qdir %s 最低支拟合失败：%s" % (list(d), e))
            if eig is not None:
                p, wmin, fit, extra = eig
                res["za_branch_index_%s" % tags[i]] = extra["branch_index"]
                res["za_zfrac_%s" % tags[i]] = extra["zfrac_first"]
                res["za_zfrac_min_%s" % tags[i]] = extra["zfrac_min"]
                res["za_branch_note_%s" % tags[i]] = extra["note"]
                if extra["note"]:
                    notes.append("%s %s" % (tags[i], extra["note"]))
            else:
                # 无本征矢量（假 phonopy）→ 退回最低支
                p, wmin, fit = lp, lw, lf
            ps.append(p)
            wmins.append(wmin)
            fits.append(fit)
            res["p_lowest"].append(lp)
            res["za_lowest_exponent_%s" % tags[i]] = lp
            res["za_lowest_r2_%s" % tags[i]] = lf["r2"] if lf else None
            res["za_lowest_log_resid_rms_%s" % tags[i]] = (
                lf["log_resid_rms"] if lf else None)
            if ph_za is not ph:
                dflt = None
                try:
                    dflt = za_power_law_eig(ph, qdir=d, qmax=qmax, n=12,
                                            cart_vac_axis=cart_vac_axis)
                except Exception:
                    dflt = None
                res["za_exponent_default_%s" % tags[i]] = (
                    dflt[0] if dflt is not None else None)
            if p is None and wmin is not None:
                notes.append("qdir %s 无正频（omega_min=%.4f THz）" % (list(d), wmin))
    finally:
        try:
            ph.nac_params = had_nac
        except Exception:
            pass
    if ph_za is ph:
        res["za_exponent_default_q1"] = ps[0]
        res["za_exponent_default_q2"] = ps[1]
    res["p"] = list(ps)
    res["min_freq"] = list(wmins)
    res["za_r2_q1"], res["za_r2_q2"] = [(f["r2"] if f else None) for f in fits]
    res["za_log_resid_rms_q1"], res["za_log_resid_rms_q2"] = [
        (f["log_resid_rms"] if f else None) for f in fits]
    margins = [None if p is None else float(min(p - p_lo, p_hi - p)) for p in ps]
    res["za_margin_q1"], res["za_margin_q2"] = margins
    res["ok"] = bool(all(p is not None and p_lo < p < p_hi for p in ps))
    pstr = ["%.3f" % p if p is not None else "None" for p in ps]
    lowstr = ["%.3f" % p if p is not None else "None" for p in res["p_lowest"]]

    # 人工复核提示：p 贴近判据边界、或 log-log 拟合太脏时，粗判 ok 可能骗人
    review = False
    if is2d:
        for tag, p, fit, m in zip(tags, ps, fits, margins):
            if p is None or m is None:
                continue
            if m < ZA_MARGIN_MIN:
                review = True
                notes.append("%s p=%.3f 距 [%.1f, %.1f] 判据边界仅 %.3f（< %.2f）"
                             % (tag, p, p_lo, p_hi, m, ZA_MARGIN_MIN))
            r2 = fit.get("r2") if fit else None
            if r2 is not None and r2 < ZA_R2_MIN:
                review = True
                notes.append("%s log-log 拟合噪声过大（R2=%.4f < %.2f，log 残差 RMS=%.4f），"
                             "p 不可信" % (tag, r2, ZA_R2_MIN, fit.get("log_resid_rms")))
    res["za_needs_review"] = bool(review)
    if res["ok"]:
        notes.append("两个面内方向都有 %.1f < p < %.1f" % (p_lo, p_hi))
    else:
        notes.append("弯曲支在每个面内方向都不是二次色散")
    if review:
        notes.append("建议人工复核，不要仅凭 ok 下结论")
    if not res["ok"]:
        notes.append("弯曲支被线性化/有虚频时 2D 的 κ 会整体失真，S6 不启动。先查 pheasy 的 "
                     "RASR 是否真在 -c 步施加了 BHH"
                     "（pheasy_c.log 里的 Imposing rotational invariance and equilibrium conditions）、"
                     "结构是否还有残余面内应力。")
    res["note"] = ("ZA 指数 p(q1,q2)=%s（新判据=本征矢量面外占比≥%.1f 的最低频支，"
                   "vac 轴原胞下标=%s/笛卡尔=%s）；旧最低支 p=%s；%s"
                   % (pstr, ZA_ZFRAC_MIN, res["vacuum_axis_in_primitive"],
                      res["vacuum_axis_cartesian"], lowstr, "；".join(notes)))
    return res


def _stability_gate(cfg, out):
    # phonopy 判虚频（q-mesh 最小频率）。2D 材料默认用【无 NAC】判据：3D 库仑核的 NAC 在
    # 严格 2D 体系近 Γ 会产生【虚假虚频】(LO-TO 在真 2D 应趋零)，与 step6 KAPPA_NAC=auto 一致。
    # NAC 与无 NAC 两个最小频率都算出来记进 summary 供对照。返回诊断 dict。
    import numpy as np
    p3dir = out / P3PY_SUB
    dim = str(cfg.get("DIM", "")).lower()
    is2d = dim.startswith("2")

    # 真空轴：vax 是 Cartesian 轴（resolve_dim 返回 POSCAR 的 Cartesian 下标）。
    #   凡是要按"原胞基矢下标"索引的地方（ZA 方向、band 路径、虚频网格）都必须先用
    #   _vacuum_axis_in_primitive 映射一次——生产原胞把 25 Å 真空放在第 0 个基矢。
    vax = 2
    for _cand in (out / "POSCAR", p3dir / "POSCAR"):
        if _cand.is_file():
            try:
                import kl_common as _kc
                _, _v = _kc.resolve_dim(_cand, "2d" if is2d else "auto")
                vax = int(_v if _v is not None else 2)
                break
            except Exception:
                pass

    # 虚频门禁的面内网格数（显式整数网格；见 _min_freq 为什么不能用 float mesh）
    IMAG_MESH_N = 60

    def _err(msg):
        return {"tool_ok": False, "stable": False, "min_freq": None,
                "min_freq_nonac": None, "min_freq_nac": None, "nac_used": False,
                "is_2d": is2d, "status": "tool_error", "note": msg}

    try:
        import phono3py
        from phono3py.file_IO import read_fc2_from_hdf5
        ph3 = phono3py.load(str(p3dir / "phono3py_disp.yaml"),
                            produce_fc=False, is_nac=False, log_level=0)
        uc, scm, pm = ph3.unitcell, ph3.supercell_matrix, ph3.primitive_matrix
        fc2 = np.asarray(read_fc2_from_hdf5(filename=str(p3dir / "fc2.hdf5")))
    except Exception as e:
        return _err("建 Phonopy / 读 fc2 失败：%s" % e)

    def _min_freq(with_nac):
        from phonopy import Phonopy
        ph = Phonopy(uc, supercell_matrix=scm, primitive_matrix=pm)
        ph.force_constants = fc2
        if with_nac:
            from phonopy.file_IO import parse_BORN
            nac = parse_BORN(ph.primitive, filename=str(p3dir / "BORN"))
            if isinstance(nac, dict) and not nac.get("factor"):
                nac["factor"] = 14.399652
            ph.nac_params = nac
        else:
            ph.nac_params = None          # 显式建"无 NAC"的动力学矩阵
        # ★ 不要用 float mesh：phonopy 把 float 当【面间距长度】（N_i=nint(l/d_i)）并强制 Γ 中心。
        #   本体系实测 mesh=60.0 → [2,22,22]，面内最近 q 只有 0.105 Å⁻¹，会漏掉近 Γ 的软模
        #   （MoS₂ 的 ZA 软模在 0.133 Å⁻¹：门禁因此报 -6.5e-8，而 S5.1 的路径报 -0.055；
        #   2026-09-22 用 fc2 直接复算查实）。改显式整数网格：真空轴 1、面内 IMAG_MESH_N。
        mesh = None
        vax_prim = None
        if is2d:                      # ★ 3D 绝不能把某轴设成 1（会只采一个 k 点）
            try:
                vax_prim = _vacuum_axis_in_primitive(ph, vax)
                mesh = [IMAG_MESH_N, IMAG_MESH_N, IMAG_MESH_N]
                mesh[int(vax_prim)] = 1
            except Exception as _e:
                print("[WARN] 虚频门禁：真空轴 → 原胞基矢映射失败（%s）" % _e)
                mesh = None
        if mesh is None:
            if is2d:
                print("[WARN] 虚频门禁：退回面间距网格 mesh=60.0（面内仅 ~22，可能漏掉近 Γ 软模）")
            ph.run_mesh(mesh=60.0, with_eigenvectors=False, is_mesh_symmetry=True)
        else:
            print("[..] 虚频门禁网格 mesh=%s（真空轴原胞下标=%d，面内 %d，Γ 中心）"
                  % (" ".join(str(x) for x in mesh), int(vax_prim), IMAG_MESH_N))
            ph.run_mesh(mesh=mesh, with_eigenvectors=False, is_mesh_symmetry=True,
                        is_gamma_center=True)
        return float(np.min(ph.get_mesh_dict()["frequencies"])), ph

    try:
        mf_nonac, ph_nonac = _min_freq(False)
    except Exception as e:
        return _err("phonopy mesh(no-NAC) 计算失败：%s" % e)

    mf_nac = None
    if (p3dir / "BORN").is_file():
        try:
            mf_nac, _ = _min_freq(True)
        except Exception as e:
            print("[WARN] 带 NAC 的 mesh 失败，仅用无 NAC 判据：%s" % e)

    # （vax 已在本函数开头解析，_min_freq 的显式网格与下面的 band_path_2d 都用它）

    # best-effort 出图（用无 NAC 版，避免 2D 的 NAC 假象污染谱图）；其最小值并入无 NAC 判据
    try:
        if is2d:
            # P2-1：seekpath 是 3D 工具，2D 上会给 Γ-A 这类 kz 线段；改用二维路径表。
            # ★ 与 ZA 同源的方向修复：路径表要的是原胞基矢下标，不是 Cartesian 下标。
            import kl_common as _kc
            from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
            _prim_vax = _vacuum_axis_in_primitive(ph_nonac, vax)
            _paths, _labels, _lat = _kc.band_path_2d(ph_nonac.primitive.cell, 101, _prim_vax)
            _bands, _conn = get_band_qpoints_and_path_connections(_paths, npoints=101)
            ph_nonac.run_band_structure(_bands, path_connections=_conn, labels=_labels)
            ph_nonac.write_yaml_band_structure(filename=str(p3dir / "band-dft-cpu.yaml"))
            print("[..] band-dft-cpu.yaml 用 2D 高对称路径（%s 格子，kz=0）" % _lat)
        else:
            ph_nonac.auto_band_structure(plot=False, write_yaml=True,
                                         filename=str(p3dir / "band-dft-cpu.yaml"))
        bf = _parse_min_freq(p3dir / "band-dft-cpu.yaml")
        if bf is not None:
            mf_nonac = min(mf_nonac, bf)
    except Exception as e:
        print("[WARN] band-dft-cpu.yaml 生成失败：%s —— ZA 判据退回【上面那套 q-mesh】的最小值；"
              "该网格若被降级（见'虚频门禁网格'行）可能漏掉近 Γ 软模" % e)

    # ---- 单一裁判（2026-09-22）：把 [IMAG_MESH_N]^3 网格 + 2D 路径上的频率合并，
    #      交给 imag_policy.classify_imag。阈值语义见 skill/_common/imag_policy.py。
    #      2D 默认用无 NAC 的频率（3D 库仑核的 NAC 在严格 2D 近 Γ 会产生虚假虚频，
    #      与 step6 KAPPA_NAC=auto 一致）；NAC 的 mesh 最小值只留作对照。
    all_freqs, all_qpts = [], []
    try:
        _md = ph_nonac.get_mesh_dict()
        all_freqs += [[float(x) for x in row] for row in _md["frequencies"]]
        all_qpts += [[float(x) for x in q] for q in _md["qpoints"]]
    except Exception as e:                          # noqa: BLE001
        print("[WARN] imag_policy：mesh 频率取用失败：%s" % e)
    try:
        _bs = ph_nonac.get_band_structure_dict()
        all_freqs += [[float(x) for x in f]
                      for seg in (_bs.get("frequencies") or []) for f in seg]
        all_qpts += [[float(x) for x in q]
                     for seg in (_bs.get("qpoints") or []) for q in seg]
    except Exception as e:                          # noqa: BLE001
        print("[WARN] imag_policy：band 频率取用失败：%s" % e)
    try:
        _vaxp = _vacuum_axis_in_primitive(ph_nonac, vax) if is2d else None
    except Exception:                               # noqa: BLE001
        _vaxp = None
    if not all_freqs:
        return _err("imag_policy：mesh/band 都没有取到频率（无法判定虚频）")
    imag = imag_policy.classify_imag(all_freqs, all_qpts, is2d, _vaxp, cfg)
    mf_used = imag["min_freq_THz"] if imag.get("min_freq_THz") is not None else mf_nonac
    nac_used = False

    # P2-2：2D 还要查 ZA 弯曲支是不是二次色散（只看最小频率不够 —— 线性化时频率全正）
    za = None
    if is2d:
        za = _za_check(cfg, ph_nonac, vax, is2d)
        if za is not None and "error" in za:
            print("[WARN] ZA 二次性检查没跑成（不拦，但 κ 的 ZA 部分未经验证）：%s" % za["error"])

    # 最终 stability_verdict = imag verdict 与 ZA 结论中更坏者（pass<warn<needs_review<fail）
    _za_v = imag_policy.za_verdict(za)
    stability_verdict = imag_policy.combine_verdict(imag["verdict"], _za_v)
    stable = imag_policy.is_stable(stability_verdict)
    _thr = imag["thresholds"]["IMAG_THR"]
    parts = ["min_freq(no-NAC mesh)=%.4f THz (%.2f cm-1)"
             % (mf_nonac, mf_nonac * imag_policy.THZ_TO_CM1)]
    if mf_nac is not None:
        parts.append("min_freq(NAC mesh)=%.4f THz" % mf_nac)
    parts.append("imag_policy=%s/%s，合并集 min=%.4f THz (%.2f cm-1)，策略 %s"
                 % (imag["verdict"], imag["imag_class"], mf_used,
                    mf_used * imag_policy.THZ_TO_CM1, imag["policy_version"]))
    if is2d:
        parts.append("ZA verdict=%s" % _za_v)
    if za is not None and "error" not in za:
        parts.append(za["note"])
    if stability_verdict == "pass":
        note = "；".join(parts) + " → 无实质虚频，稳定"
    elif stability_verdict == "warn":
        note = "；".join(parts) + (" → 近 Γ 声学支软化（|ν|=%.4f THz，≤ IMAG_THR=%.2f）；"
                                   "warn 放行 S6，但 κ 可能低估近 Γ ZA 贡献" % (abs(mf_used), _thr))
    elif stability_verdict == "needs_review":
        note = "；".join(parts) + " → ZA 弯曲支需人工复核（放行，行为不变）"
    else:
        note = "；".join(parts) + (" → 存在虚频（%s，|ν|=%.4f THz）；动力学不稳定，S6 不启动"
                                   % (imag["imag_class"], abs(mf_used)))
    return {"tool_ok": True, "stable": stable, "min_freq": mf_used,
            "min_freq_nonac": mf_nonac, "min_freq_nac": mf_nac, "nac_used": nac_used,
            "is_2d": is2d, "za_exponent": za, "imag": imag, "za_verdict": _za_v,
            "stability_verdict": stability_verdict,
            "status": ("stable" if stable else
                       ("imaginary" if imag["verdict"] == "fail" else "za_not_quadratic")),
            "note": note}


def cmd_post(cfg):
    out = Path.cwd()
    # ShengBTE 导出已在 fit_phono3py / collect_pheasy 就地完成；这里只查状态。
    want_sb = str(cfg.get("EXPORT_SHENGBTE", "true")).lower() in ("true", "1", "yes")
    sb_ok = (out / SB_SUB / "FORCE_CONSTANTS_3RD").is_file() if want_sb else None
    if want_sb and not sb_ok:
        print("[WARN] EXPORT_SHENGBTE=true 但未生成 shengbte/FORCE_CONSTANTS_3RD"
              "（多半是 hiphive 缺失或 fc 非 full）。SOLVER=shengbte 时 S6 会报缺文件。")

    g = _stability_gate(cfg, out)
    stable, tool_ok = g["stable"], g["tool_ok"]
    print("[%s] %s" % ("OK" if stable else "FAIL", g["note"]))

    # 拟合质量指标（pheasy 引擎的模板会写 fit_metrics.json；phono3py/symfc 没有 →
    #   缺这个文件是正常的，不报错）。A/B 对照（RASR=none vs BHH）靠它读 rel_err、
    #   fc2max/fc3max、零空间自由度。
    _fm = {}
    _fmp = out / "fit_metrics.json"
    if _fmp.is_file():
        try:
            _fm = json.loads(_fmp.read_text(encoding="utf-8"))
            print("[..] fit_metrics.json：rel_err=%s RMSE=%s fc2max=%s fc3max=%s rasr=%s"
                  % (_fm.get("pheasy_relative_error"), _fm.get("pheasy_rmse_eV_per_A"),
                     _fm.get("fc2max"), _fm.get("fc3max"), _fm.get("rasr_applied")))
        except Exception as _e:
            print("[WARN] fit_metrics.json 解析失败：%s" % _e)
    _im = g.get("imag") or {}
    (out / "phonon_summary.json").write_text(json.dumps(
        {"stable": bool(stable), "status": g["status"], **_fm,
         "imaginary_frequency": bool(tool_ok and g.get("stability_verdict") == "fail"),
         # ---- imag_policy（虚频判据唯一真源）：第一条的全部字段 ----
         "stability_verdict": g.get("stability_verdict"),
         "verdict": _im.get("verdict"),
         "imag_class": _im.get("imag_class"),
         "min_freq_THz": _im.get("min_freq_THz"),
         "min_freq_cm1": _im.get("min_freq_cm1"),
         "q_at_min_frac": _im.get("q_at_min_frac"),
         "q_norm_at_min": _im.get("q_norm_at_min"),
         "branch_at_min": _im.get("branch_at_min"),
         "n_neg_qpoints": _im.get("n_neg_qpoints"),
         "thresholds": _im.get("thresholds"),
         "policy_version": _im.get("policy_version"),
         "imag_policy_refs": list(imag_policy.POLICY_REFS),
         "za_verdict": g.get("za_verdict"),
         # ---- 历史字段（保留：下游/对照仍在读，勿删）----
         "min_frequency_THz": g["min_freq"],
         "min_freq_nonac_THz": g["min_freq_nonac"],
         "min_freq_nac_THz": g["min_freq_nac"],
         "nac_used_for_verdict": g["nac_used"], "is_2d": g["is_2d"],
         "fit_engine": cfg.get("FIT_ENGINE"), "fc_calc": cfg.get("FC_CALC"),
         "pheasy_method": cfg.get("PHEASY_FIT_METHOD"),
         "shengbte_export": sb_ok, "tool_ok": tool_ok, "note": g["note"]},
        ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    print("[DONE] post：phonon_summary.json 就绪，status=%s stable=%s stability_verdict=%s"
          % (g["status"], str(stable).lower(), g.get("stability_verdict")))
    # 区分两种"不通过"：
    #   工具错误(tool_ok=False，mesh 没算成) → 作业失败退出，tf 标 error，提醒去查日志；
    #   真有虚频(status=imaginary) → 正常退出，marker 不满足，S6 被合理挡住(非 error)。
    if not tool_ok:
        sys.exit("[ERROR] 虚频闸工具错误，作业按失败结束（见上 note / phono3py/band-dft-cpu.log）")


# ==========================================================================
CMDS = {"prep": cmd_prep, "fit_phono3py": cmd_fit_phono3py,
        "collect_pheasy": cmd_collect_pheasy, "post": cmd_post}


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in CMDS:
        sys.exit("用法：kl_fc_backends.py <%s> fit_config.json" % "|".join(CMDS))
    CMDS[sys.argv[1]](load_cfg(sys.argv[2]))


if __name__ == "__main__":
    main()