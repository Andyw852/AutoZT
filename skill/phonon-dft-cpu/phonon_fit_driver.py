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
import sys
from pathlib import Path

import numpy as np


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


# ==========================================================================
# 2D ZA 弯曲支二次性检查（za_2d 是唯一实现，本函数只做编排）
# ==========================================================================
def _za_check_2d(ph, cart_vac_axis, qmax=0.05):
    import za_2d as Z
    p_lo, p_hi = Z.ZA_P_RANGE
    res = {"qmax": float(qmax), "p_range": [p_lo, p_hi], "dirs": [], "p": [],
           "min_freq": [], "r2": [], "za_symprec": None, "za_symmetry_relaxed": False}
    try:
        vax_prim = Z.vacuum_axis_in_primitive(ph, cart_vac_axis)
        dirs = Z.inplane_qdirs(ph, vax_prim)
    except Exception as e:                                  # noqa: BLE001
        res["error"] = "ZA 面内方向解析失败：%s" % e
        return res
    res["dirs"] = [[float(x) for x in d] for d in dirs]
    ph_za = ph
    try:
        _relaxed = Z.relaxed_symprec(ph)
        if _relaxed is not None:
            ph_za = Z.clone_phonopy(ph, _relaxed)
            res["za_symprec"] = _relaxed
    except Exception:                                       # noqa: BLE001
        ph_za = ph
    res["za_symmetry_relaxed"] = bool(ph_za is not ph)
    had_nac = getattr(ph, "nac_params", None)
    try:
        ph.nac_params = None
        try:
            ph_za.nac_params = None
        except Exception:                                   # noqa: BLE001
            pass
        for d in dirs:
            got = None
            try:
                got = Z.za_power_law_eig(ph_za, qdir=d, qmax=qmax, n=12,
                                         cart_vac_axis=cart_vac_axis)
            except Exception as e:                          # noqa: BLE001
                res.setdefault("za_notes", []).append(
                    "qdir %s 本征矢量 ZA 拟合失败：%s" % (list(d), e))
            if got is None:
                res["p"].append(None)
                res["min_freq"].append(None)
                res["r2"].append(None)
                continue
            p, wmin, fit, _extra = got
            res["p"].append(p)
            res["min_freq"].append(wmin)
            r2 = (fit or {}).get("r2") if isinstance(fit, dict) else None
            res["r2"].append(r2)
            if p is None or (r2 is not None and r2 < Z.ZA_R2_MIN):
                res["za_needs_review"] = True
            elif (abs(p - p_lo) < Z.ZA_MARGIN_MIN) or (abs(p - p_hi) < Z.ZA_MARGIN_MIN):
                res["za_needs_review"] = True
    finally:
        try:
            ph.nac_params = had_nac
        except Exception:                                   # noqa: BLE001
            pass
    res["ok"] = all(p is not None and p_lo < p < p_hi for p in res["p"])
    ps = ["%.2f" % p if p is not None else "None" for p in res["p"]]
    res["note"] = ("ZA 二次性满足（p=%s；本征矢量判据）" % ps) if res["ok"] else (
        "ZA 不是二次色散（p=%s，要求 %.1f~%.1f）：弯曲支被线性化或有虚频，2D 的稳定性结论不可信。"
        " DFT 的 fc2 若未施加 Born-Huang 旋转不变约束，近 Γ ZA 会被线性化。"
        % (ps, p_lo, p_hi))
    return res


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

    # ---- 虚频门禁的 q 网格：★ 显式整数 Γ 中心网格，不能用 float mesh ----
    #   phonopy 把 float 当【面间距长度】(N_i=nint(l/d_i)) 并强制 Γ 中心：实测 mesh=60.0
    #   → [2,22,22]，面内最近 q 只有 0.105 Å⁻¹，会漏掉近 Γ 软模而假 pass（kl-dft-cpu
    #   2026-09-22 用 fc2 复算查实）。2D：真空轴=1、面内 60。
    IMAG_MESH_N = 60
    if is2d:
        try:
            _vaxp = kc.vacuum_axis_in_primitive(ph.primitive.cell, vac_axis)
        except Exception as e:                              # noqa: BLE001
            print("[WARN] 真空轴 → 原胞基矢映射失败（%s），退回 c 轴" % e)
            _vaxp = 2
        _mesh = [IMAG_MESH_N, IMAG_MESH_N, IMAG_MESH_N]
        _mesh[int(_vaxp)] = 1
        print("[..] 虚频门禁网格 mesh=%s（真空轴原胞下标=%d，面内 %d，Γ 中心）"
              % (" ".join(str(x) for x in _mesh), int(_vaxp), IMAG_MESH_N))
        ph.run_mesh(mesh=_mesh, with_eigenvectors=False, is_mesh_symmetry=True,
                    is_gamma_center=True)
    else:
        ph.run_mesh(mesh=60.0, with_eigenvectors=False, is_mesh_symmetry=True)
    mesh_freqs = np.asarray(ph.get_mesh_dict()["frequencies"], dtype=float)
    if mesh_freqs.size == 0 or not np.all(np.isfinite(mesh_freqs)):
        sys.exit("[ERROR] mesh frequencies 为空或含非有限值")
    mf = float(np.min(mesh_freqs))

    # ---- band（2D 用 kz=0 的二维高对称路径；seekpath 是 3D 工具）----
    try:
        if is2d:
            from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
            _paths, _labels, _lat = kc.band_path_2d(ph.primitive.cell, bp, int(_vaxp))
            _bands, _conn = get_band_qpoints_and_path_connections(_paths, npoints=bp)
            ph.run_band_structure(_bands, path_connections=_conn, labels=_labels)
            ph.write_yaml_band_structure(filename="band-dft-cpu.yaml")
            print("[..] band-dft-cpu.yaml 用 2D 高对称路径（%s 格子，kz=0，%d 点/段）"
                  % (_lat, bp))
        else:
            ph.auto_band_structure(plot=False, write_yaml=True, filename="band-dft-cpu.yaml")
    except Exception as e:                                  # noqa: BLE001
        sys.exit("[ERROR] band-dft-cpu.yaml 生成失败：%s" % e)

    # ---- ZA 二次性（仅 2D）----
    za = None
    _za_on = str(params.get("ZA_CHECK") or "auto").strip().lower() not in \
        ("off", "false", "0", "no")
    if is2d and _za_on:
        za = _za_check_2d(ph, vac_axis, float(params.get("ZA_QMAX") or 0.05))
        print("[..] ZA 检查：%s" % za.get("note", za.get("error", "")))

    # ---- 虚频判据统一走 imag_policy（唯一真源）：mesh + band 合并 ----
    import imag_policy

    def _imag_samples(_ph):
        _f, _q = [], []
        try:
            _md = _ph.get_mesh_dict()
            _f += [[float(x) for x in row] for row in _md["frequencies"]]
            _q += [[float(x) for x in q] for q in _md["qpoints"]]
        except Exception:                                   # noqa: BLE001
            pass
        try:
            _bd = _ph.get_band_structure_dict()
            _f += [[float(x) for x in x2]
                   for seg in (_bd.get("frequencies") or []) for x2 in seg]
            _q += [[float(x) for x in x2]
                   for seg in (_bd.get("qpoints") or []) for x2 in seg]
        except Exception:                                   # noqa: BLE001
            pass
        return _f, _q

    _fr, _qp = _imag_samples(ph)
    try:
        _vaxp2 = kc.vacuum_axis_in_primitive(ph.primitive.cell, vac_axis) if is2d else None
    except Exception:                                       # noqa: BLE001
        _vaxp2 = None
    _cfg = {"IMAG_THR": params.get("IMAG_THR"),
            "IMAG_THR_STRICT": params.get("IMAG_THR_STRICT"),
            "IMAG_QGAMMA": params.get("IMAG_QGAMMA"),
            "IMAG_QGAMMA_GRACE": params.get("IMAG_QGAMMA_GRACE")}
    _imag = imag_policy.classify_imag(_fr, _qp, is2d, _vaxp2, _cfg)
    _za_v = imag_policy.za_verdict(za)
    stability_verdict = imag_policy.combine_verdict(_imag["verdict"], _za_v)
    stable = imag_policy.is_stable(stability_verdict)
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
        "mesh_min_frequency_THz": mf,
        "mesh_n": IMAG_MESH_N,
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
