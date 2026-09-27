#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""phonon_fit_driver.py —— phonopy 收力 + fc2 拟合 + 声子谱（phonon-dft-cpu S3）。

计算节点作业里跑，cwd = step3_phonon/：读 POSCAR + kl_params（SUPERCELL/FD_DISTANCE/
MESH/IMAG_*）→ 收 ../step2_disp/disp-*/vasprun.xml 的力 → produce_force_constants →
显式整数 Γ 中心网格（2D 真空轴 1、面内 60）+ 2D 高对称路径 band → ZA 二次性（2D）+
imag_policy.classify_imag（mesh+band 合并的【唯一裁判】）→ phonon_summary.json。

2026-09-27 与 kl-dft-cpu S5_fc / phonon-mlff-cpu S3 的虚频判据对齐：
  · 旧版用 20³ mesh 的全局最小频率 + 本地 IMAG_THR 阈值（且 step.conf 的 IMAG_THR 因
    kl_params 没写而永远被忽略）；现在显式整数网格 + 2D 路径并集，交 imag_policy 判定；
  · 2D 增加 ZA 弯曲支二次性检查（ω∝q^p，1.7<p<2.3）：ZA 被线性化时频率全正、照样过
    虚频闸，但 2D 的稳定性/κ 结论是错的。
  · 只有 fail 挡下一步（warn / needs_review 放行），与 kl-dft-cpu 行为一致。
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

import phonon_stability as ps


# ==========================================================================
# 参数 / 读力
# ==========================================================================
def read_params(path):
    d = {}
    p = Path(path)
    if p.is_file():
        for ln in p.read_text(errors="ignore").splitlines():
            s = ln.strip()
            if s and not s.startswith("#") and "=" in s:
                k, v = s.split("=", 1)
                d[k.strip().upper()] = v.strip()
    return d


def read_vasprun_force(path):
    """从 vasprun.xml 读最后一个结构块的力（N 原子 × 3）。"""
    import xml.etree.ElementTree as ET
    tree = ET.parse(str(path))
    forces = None
    for calc in tree.getroot().iter("calculation"):
        for va in calc.iter("varray"):
            if va.get("name") == "forces":
                rows = [[float(x) for x in (v.text or "").split()]
                        for v in va.findall("v")]
                if rows:
                    forces = rows
    if forces is None:
        raise RuntimeError("vasprun.xml 里没有 forces")
    return np.array(forces, dtype=float)


def _export_force_constants(ph, cwd):
    """Best-effort 标准 fc2 导出；导出失败不中止声子计算。"""
    try:
        from phonopy.file_IO import write_force_constants_to_hdf5
        write_force_constants_to_hdf5(ph.force_constants, filename=str(cwd / "fc2.hdf5"))
        print("[OK] 已写出 fc2.hdf5")
    except Exception as e:                                  # noqa: BLE001
        print("[WARN] fc2.hdf5 导出失败：%s" % e)
    try:
        from phonopy.file_IO import write_FORCE_CONSTANTS
        write_FORCE_CONSTANTS(ph.force_constants, filename=str(cwd / "FORCE_CONSTANTS"))
        print("[OK] 已写出 FORCE_CONSTANTS")
    except Exception as e:                                  # noqa: BLE001
        print("[WARN] FORCE_CONSTANTS 导出失败：%s" % e)


def _cleanup(path):
    """删掉失败时留下的空/半成品文件（best-effort）。"""
    try:
        Path(path).unlink()
    except OSError:
        pass


def _load_nac(ph, cwd):
    """有 step2_nac/vasprun.xml → 用 phonopy-vasp-born 生成 BORN → ph.nac_params。

    返回 (ok, note)。任何一步失败都退回"无 NAC"（只告警，不中止声子计算）；
    产物不物理（介电张量非有限/对角非正 —— 金属或 DFPT 发散的常见表现）同样退回。
    """
    import subprocess
    nac_dir = cwd.parent / "step2_nac"
    vr = nac_dir / "vasprun.xml"
    born = cwd / "BORN"
    if not nac_dir.is_dir():
        return False, "未检测到 step2_nac（NAC 未启用或未跑）—— 本谱无 LO-TO"
    if not vr.is_file():
        return False, "step2_nac 存在但缺 vasprun.xml（未算完/未拉回）—— 本谱无 LO-TO"
    if not born.is_file():
        if shutil.which("phonopy-vasp-born") is None:
            return False, ("找不到 phonopy-vasp-born（作业里 conda activate 了吗？"
                           "它是 phonopy 自带的命令）")
        try:
            with open(str(born), "w", encoding="utf-8") as _fh:
                _r = subprocess.run(["phonopy-vasp-born", str(vr)], stdout=_fh,
                                    stderr=subprocess.PIPE, text=True)
        except Exception as e:                              # noqa: BLE001
            _cleanup(born)
            return False, "phonopy-vasp-born 调用失败（%s）" % e
        if _r.returncode != 0 or born.stat().st_size == 0:
            _cleanup(born)
            return False, ("phonopy-vasp-born 失败（rc=%d）：%s"
                           % (_r.returncode, (_r.stderr or "")[-300:]))
    if not born.is_file() or born.stat().st_size == 0:
        _cleanup(born)
        return False, "BORN 生成失败或为空"
    try:
        from phonopy.file_IO import parse_BORN
        nac = parse_BORN(ph.primitive, filename=str(born))
        if isinstance(nac, dict) and not nac.get("factor"):
            nac["factor"] = 14.399652          # e^2/(4π ε0) 单位换算，phonopy 默认
        eps = np.asarray(nac.get("dielectric"), dtype=float)
        if eps.shape != (3, 3) or not np.all(np.isfinite(eps)):
            return False, "BORN 介电张量非法（形状 %s）" % (eps.shape,)
        if np.any(np.diag(eps) <= 0):
            return False, "BORN ε∞ 对角含非正值（金属/DFPT 发散？）—— 退回无 NAC"
        ph.nac_params = nac
        return True, "BORN 已加载并写入 ph.nac_params"
    except Exception as e:                                  # noqa: BLE001
        return False, "读 BORN 失败：%s" % e


def main():
    cwd = Path.cwd()
    sys.path.insert(0, str(cwd))
    import kl_common as kc
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp

    if not (cwd / "POSCAR").is_file():
        sys.exit("[ERROR] 缺 POSCAR")
    params = read_params(cwd / "kl_params.txt")
    _sc = [int(x) for x in (params.get("SUPERCELL") or "1 1 1").split()]
    # 3 个数=对角；9 个数=3×3 矩阵（行主序，Phonopy 的 supercell_matrix 直接用）
    dim = (np.array(_sc, dtype=int).reshape(3, 3) if len(_sc) == 9
           else np.diag(np.array(_sc, dtype=int)))
    fd = float(params.get("FD_DISTANCE") or 0.01)
    bp = int(params.get("BAND_POINTS") or 51)

    uc = read_vasp("POSCAR")
    ph = Phonopy(uc, supercell_matrix=np.array(dim, dtype=int), primitive_matrix="auto")
    ph.generate_displacements(distance=fd)
    n = len(ph.supercells_with_displacements)
    print("[..] 超胞=%s 位移=%d 帧 FD_DISTANCE=%s"
          % (" ".join(str(x) for x in dim), n, fd))

    forces = []
    for i in range(1, n + 1):
        vasprun = cwd.parent / "step2_disp" / ("disp-%05d" % i) / "vasprun.xml"
        if not vasprun.is_file():
            sys.exit("[ERROR] 缺 %s" % vasprun)
        forces.append(read_vasprun_force(vasprun))
    ph.forces = np.array(forces, dtype=float)
    print("[OK] 收力 %d 帧" % len(forces))

    ph.produce_force_constants()
    print("[OK] fc2 拟合完成")
    _export_force_constants(ph, cwd)

    # ---- 稳定性判据一律用【无 NAC】对象（与 kl-dft-cpu 一致）----
    #   3D 库仑核的 NAC 在严格 2D 近 Γ 会产生【虚假虚频】，会污染虚频闸；NAC 只用于成图
    #   （LO-TO 劈裂）。判据对象与成图对象共享同一份 fc2，只是 nac_params=None。
    ph_gate = Phonopy(uc, supercell_matrix=np.array(dim, dtype=int), primitive_matrix="auto")
    ph_gate.force_constants = ph.force_constants
    ph_gate.nac_params = None

    # ---- 维度（2D 的网格/路径/ZA 都要用）----
    dimtag = (params.get("DIM") or "").lower()
    if not dimtag:
        try:
            dimtag, _v = kc.resolve_dim(cwd / "POSCAR", "auto")
        except Exception:                                   # noqa: BLE001
            dimtag = "3d"
    is2d = str(dimtag).startswith("2")
    vac_axis = 2
    if is2d:
        try:
            _, _v = kc.resolve_dim(cwd / "POSCAR", "2d")
            vac_axis = int(_v if _v is not None else 2)
        except Exception:                                   # noqa: BLE001
            vac_axis = 2

    # ---- 虚频门禁 q 网格：单个真源 phonon_stability.mesh_numbers（2D 真空轴=1、面内 60）----
    try:
        _mesh, _vaxp = ps.mesh_numbers(ph_gate, is2d, vac_axis, ps.IMAG_MESH_N)
    except Exception as e:                                  # noqa: BLE001
        print("[WARN] 真空轴 → 原胞基矢映射失败（%s），退回面间距网格" % e)
        _mesh, _vaxp = None, 2
    _mesh_used = None
    if _mesh is not None:                       # 2D：真空轴 1、面内 60，Γ 中心
        _mesh_used = _mesh
        print("[..] 虚频门禁网格 mesh=%s（真空轴原胞下标=%d，面内 %d，Γ 中心）"
              % (" ".join(str(x) for x in _mesh), int(_vaxp), ps.IMAG_MESH_N))
        ph_gate.run_mesh(mesh=_mesh, with_eigenvectors=False, is_mesh_symmetry=True,
                         is_gamma_center=True)
    elif is2d:                                  # 2D 真空轴映射失败：退回面间距网格
        _mesh_used = "60.0(float)"
        print("[WARN] 2D 真空轴映射失败：退回面间距网格 mesh=60.0（面内仅 ~22，可能漏掉近 Γ 软模）")
        ph_gate.run_mesh(mesh=60.0, with_eigenvectors=False, is_mesh_symmetry=True)
    else:                                       # 3D：用 step.conf 的 PHONON_MESH（显式整数，Γ 中心）
        _mesh3 = [int(x) for x in (params.get("MESH") or "20 20 20").split()]
        if len(_mesh3) != 3 or min(_mesh3) < 1:
            _mesh3 = [20, 20, 20]
        _mesh_used = _mesh3
        print("[..] 虚频门禁网格 mesh=%s（3D，Γ 中心；PHONON_MESH）"
              % " ".join(str(x) for x in _mesh3))
        ph_gate.run_mesh(mesh=_mesh3, with_eigenvectors=False, is_mesh_symmetry=True,
                         is_gamma_center=True)
    mesh_freqs = np.asarray(ph_gate.get_mesh_dict()["frequencies"], dtype=float)
    if mesh_freqs.size == 0 or not np.all(np.isfinite(mesh_freqs)):
        sys.exit("[ERROR] mesh frequencies 为空或含非有限值")
    mf = float(np.min(mesh_freqs))

    # ---- NAC：有 step2_nac/vasprun.xml → BORN → ph.nac_params（决定成图是否带 LO-TO）----
    nac_ok, nac_note = _load_nac(ph, cwd)
    if nac_ok and is2d:
        nac_note += ("；★ 2D：phonopy 只有 3D-NAC 方案（LO-TO 在 q->0 应趋零），"
                     "NAC 谱是近似，请对照 band_noNAC.yaml")
    print("[%s] NAC：%s" % ("OK" if nac_ok else "..", nac_note))

    # ---- band 路径（2D 用 kz=0 的二维高对称路径；seekpath 是 3D 工具）----
    def _run_band(_ph, fname=None):
        if is2d:
            from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
            _paths, _labels, _lat = kc.band_path_2d(_ph.primitive.cell, bp, int(_vaxp))
            _bands, _conn = get_band_qpoints_and_path_connections(_paths, npoints=bp)
            _ph.run_band_structure(_bands, path_connections=_conn, labels=_labels)
            if fname:
                _ph.write_yaml_band_structure(filename=fname)
            return _lat
        _ph.auto_band_structure(plot=False, write_yaml=bool(fname),
                                filename=(fname or "band-dft-cpu.yaml"))
        return None

    try:
        # 门禁用无 NAC 路径；无 NAC 时它同时写成 band-dft-cpu.yaml
        _lat = _run_band(ph_gate, None if nac_ok else "band-dft-cpu.yaml")
        if is2d:
            print("[..] band 用 2D 高对称路径（%s 格子，kz=0，%d 点/段）" % (_lat, bp))
        if nac_ok:
            _run_band(ph, "band-dft-cpu.yaml")               # 带 NAC（LO-TO）
            _run_band(ph_gate, "band_noNAC.yaml")            # 无 NAC 对照
            print("[..] band-dft-cpu.yaml = 含 NAC（LO-TO）；band_noNAC.yaml = 无 NAC 对照")
    except Exception as e:                                  # noqa: BLE001
        sys.exit("[ERROR] band-dft-cpu.yaml 生成失败：%s" % e)

    # ---- ZA 二次性（仅 2D；单个真源 phonon_stability.za_check）----
    za = ps.za_check(params, ph_gate, vac_axis, is2d)
    if za is not None:
        print("[..] ZA 检查：%s" % za.get("note", za.get("error", "")))

    # ---- 虚频判据：编排在 phonon_stability（唯一真源 imag_policy 分类 mesh+band 合并集）----
    import imag_policy
    _cfg = {"IMAG_THR": params.get("IMAG_THR"),
            "IMAG_THR_STRICT": params.get("IMAG_THR_STRICT"),
            "IMAG_QGAMMA": params.get("IMAG_QGAMMA"),
            "IMAG_QGAMMA_GRACE": params.get("IMAG_QGAMMA_GRACE")}
    _imag, _vaxp2 = ps.classify(ph_gate, is2d, vac_axis, _cfg)
    _za_v = imag_policy.za_verdict(za)
    stability_verdict, stable = ps.finalize(_imag, za)
    stable_mesh = (_imag["verdict"] != "fail")
    print("[..] imag_policy=%s/%s  min=%.4f THz (%.2f cm-1)  ZA=%s  最终=%s"
          % (_imag["verdict"], _imag["imag_class"], _imag["min_freq_THz"] or 0.0,
             _imag["min_freq_cm1"] or 0.0, _za_v, stability_verdict))

    if stability_verdict == "pass":
        note = "imag_policy=%s/%s；无实质虚频，稳定" % (_imag["verdict"], _imag["imag_class"])
    elif stability_verdict == "warn":
        note = ("imag_policy=%s/%s（近 Γ 声学支软化）：warn 放行，但下游对近 Γ ZA "
                "贡献可能低估" % (_imag["verdict"], _imag["imag_class"]))
    elif stability_verdict == "needs_review":
        note = "imag_policy=%s/%s；ZA 弯曲支需人工复核（放行，行为不变）" % (
            _imag["verdict"], _imag["imag_class"])
    else:
        note = ("imag_policy=%s/%s：存在虚频（%.4f THz），动力学不稳定，不可用于后续计算"
                % (_imag["verdict"], _imag["imag_class"], _imag["min_freq_THz"] or 0.0))
    if za is not None and "error" not in za:
        note += "；" + za["note"]
    note += "；NAC=%s（%s）" % ("on" if nac_ok else "off", nac_note)

    summary = {
        "PHONON_DONE": True,
        "tool_ok": True,
        "dim": str(dimtag).upper(),
        "is_2d": bool(is2d),
        # ---- imag_policy（唯一真源，对齐 kl-dft-cpu S5_fc）----
        "stable": bool(stable),
        "stable_mesh": bool(stable_mesh),
        "stability_verdict": stability_verdict,
        "verdict": _imag["verdict"],
        "imag_class": _imag["imag_class"],
        "min_freq_THz": _imag["min_freq_THz"],
        "min_freq_cm1": _imag["min_freq_cm1"],
        "q_at_min_frac": _imag["q_at_min_frac"],
        "q_norm_at_min": _imag["q_norm_at_min"],
        "branch_at_min": _imag["branch_at_min"],
        "n_neg_qpoints": _imag["n_neg_qpoints"],
        "thresholds": _imag["thresholds"],
        "policy_version": _imag["policy_version"],
        "imag_policy_refs": list(imag_policy.POLICY_REFS),
        # ---- ZA / 采样 ----
        "za_exponent": za,
        "za_verdict": _za_v,
        # ---- NAC（LO-TO 劈裂）：band-dft-cpu.yaml 是否含 NAC + 无 NAC 对照文件 ----
        "nac_applied": bool(nac_ok),
        "nac_note": nac_note,
        "band_file": "band-dft-cpu.yaml",
        "band_noNAC_file": ("band_noNAC.yaml" if nac_ok else None),
        "mesh_min_frequency_THz": mf,
        "mesh_n": ps.IMAG_MESH_N,       # 2D 门禁面内网格常数（兼容旧字段）
        "mesh_used": (" ".join(str(x) for x in _mesh_used)
                      if isinstance(_mesh_used, (list, tuple)) else _mesh_used),
        "band_points": bp,
        # 历史字段（下游 marker/兼容）：合并集最小值
        "min_frequency_THz": (_imag["min_freq_THz"]
                              if _imag["min_freq_THz"] is not None else mf),
        "n_disp": int(n),
        "supercell": " ".join(str(x) for x in dim),
        "n_atoms_sc": int(len(ph.supercell)),
        "fd_distance_A": fd,
        "assessment_policy": ("stable 由 imag_policy.classify_imag（mesh + 2D band 合并；"
                              "近 Γ 声学支 warn、其余 fail）与 ZA 二次性结论中更坏者决定；"
                              "只有 fail 挡下一步。阈值见 skill/_common/imag_policy.py。"),
        "note": note,
    }
    (cwd / "phonon_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + chr(10), encoding="utf-8")
    print("[DONE] stable=%s verdict=%s min_freq=%.4f THz"
          % (str(stable).lower(), stability_verdict, _imag["min_freq_THz"] or 0.0))


if __name__ == "__main__":
    main()
