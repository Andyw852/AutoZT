#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""phonon_fit_driver.py —— 2 阶力常数拟合 + 声子谱（phonon-mlff-cpu S3）。

计算节点作业里跑，cwd = step3_phonon/：读 POSCAR + disps.npy + forces.npy →
symfc 拟合 fc2 → q-mesh 最小频率（虚频闸）+ band-dft-cpu.yaml → phonon_summary.json。
"""
import json
import re
import sys
from pathlib import Path

import numpy as np


def _read_force_constants(path):
    try:
        from phonopy.file_IO import read_force_constants_hdf5
        return read_force_constants_hdf5(path)
    except ImportError:
        from phonopy.file_IO import parse_FORCE_CONSTANTS
        return parse_FORCE_CONSTANTS(path)


def _export_force_constants(ph, cwd):
    """Best-effort standard fc2 exports; export failures do not abort phonons."""
    try:
        from phonopy.file_IO import write_force_constants_to_hdf5
        write_force_constants_to_hdf5(ph.force_constants, filename=str(cwd / "fc2.hdf5"))
        print("[OK] 已写出 fc2.hdf5")
    except Exception as e:
        print("[WARN] fc2.hdf5 导出失败：%s" % e)
    try:
        from phonopy.file_IO import write_FORCE_CONSTANTS
        write_FORCE_CONSTANTS(ph.force_constants, filename=str(cwd / "FORCE_CONSTANTS"))
        print("[OK] 已写出 FORCE_CONSTANTS")
    except Exception as e:
        print("[WARN] FORCE_CONSTANTS 导出失败：%s" % e)


def _nfree_fc2_alm(ph, cutoff_a):
    """估算给定 cutoff 下不可约 2 阶力常数元数（与 gen_step2 的 ALM 反推同源）。

    位移帧数 n_disp 是否足以确定 symfc 拟合：方程数 = n_disp × 3N_sc（每个原子
    3 个力分量），未知数 = 不可约 2 阶力常数元数 nfree_fc2。方程数 < 未知数即拟合
    欠定，此时虚频不可信。缺 alm/结构异常时返回 None（不阻塞拟合，只把
    fit_underdetermined 记为 null）。
    """
    try:
        from ase import Atoms
        from alm import ALM
    except Exception:
        return None
    try:
        sc = ph.supercell
        atoms = Atoms(numbers=sc.numbers, positions=sc.positions, cell=sc.cell, pbc=True)
        nkd = len(set(atoms.numbers))
        cut = np.full((1, nkd, nkd), float(cutoff_a), dtype=float)
        with ALM(np.array(atoms.cell), atoms.get_scaled_positions(),
                 atoms.get_atomic_numbers(), verbosity=0) as a:
            a.define(1, cutoff_radii=cut)
            a.suggest()
            return int(a._get_number_of_irred_fc_elements(1))
    except Exception as e:
        print("[WARN] ALM 数 2 阶力常数元数失败：%s" % e)
        return None


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


def main():
    cwd = Path.cwd()
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp

    for f in ("POSCAR", "disps.npy", "forces.npy"):
        if not (cwd / f).is_file():
            sys.exit("[ERROR] 缺 %s（step2 取力没跑完？）" % f)

    uc = read_vasp("POSCAR")
    params = read_params(cwd / "klmlff_params.txt")
    _dim = np.array([int(x) for x in (params.get("SUPERCELL") or "1 1 1").split()],
                    dtype=int)
    # 3 个数=对角扩胞；9 个数=3×3 矩阵（行主序，与 phonopy/phono3py --dim 同义）
    scm = _dim.reshape(3, 3) if _dim.size == 9 else np.diag(_dim)
    ph = Phonopy(uc, supercell_matrix=scm, primitive_matrix="P")

    disps = np.ascontiguousarray(np.load("disps.npy"), dtype="double")
    forces = np.ascontiguousarray(np.load("forces.npy"), dtype="double")
    prefit = cwd / "fc2.hdf5"
    if not prefit.is_file():
        ph.dataset = {"displacements": disps, "forces": forces}

    # symfc dense 求解器的正规方程 X^T X 内存 ~ O(basis^2)，大 supercell 会爆内存
    # （qHPC60 976 原子：无 cutoff 需 ~585 GiB）。cutoff 截断 fc2 非零范围，与
    # FC_CUTOFF 是计算/建模截断；MACE receptive fields 与多体混合导数意味着
    # r_max 不能证明 fc2 在其外严格为零，必须通过 cutoff 收敛性检查验证。
    fc_cutoff = params.get("FC_CUTOFF") or "6.0"
    print("[..] %d 帧 × %d 原子，symfc 拟合 fc2 (cutoff=%s A)"
          % (len(disps), len(uc.numbers) * int(np.prod(dim)), fc_cutoff))

    if prefit.is_file():
        print("[..] 使用 pheasy-gpu 预拟合 fc2.hdf5")
        ph.force_constants = _read_force_constants(prefit)
    else:
        ph.produce_force_constants(
            fc_calculator="symfc",
            fc_calculator_options="cutoff = %s" % fc_cutoff,
        )
    print("[OK] fc2 拟合完成")
    _export_force_constants(ph, cwd)

    # ---- 维度先判：虚频门禁的网格与 2D 路径都要用它（2026-09-22）----
    import klmlff_common as kmc
    dimtag, vac_axis = "3d", 2
    try:
        from dim_common import detect_dimension
        dimtag, vac_axis, _dn = detect_dimension("POSCAR")
    except Exception as e:
        print("[WARN] 维度判定失败（按 3D 处理）：%s" % e)
    is2d = str(dimtag).lower().startswith("2")
    # ★ 2D 用显式整数网格（真空轴 1、面内 60）：float mesh 是"面间距长度"语义，
    #   实测得 [2,22,22]、面内最近 q 只 0.105 Å⁻¹，会漏掉近 Γ 软模（与 kl-dft-cpu 同一 bug）。
    if is2d:
        _imag_mesh, _vaxp = kmc.imag_mesh_numbers(ph, vac_axis or 2)
        print("[..] 虚频门禁网格 mesh=%s（真空轴原胞下标=%d）"
              % (" ".join(str(x) for x in _imag_mesh), _vaxp))
        ph.run_mesh(mesh=_imag_mesh, with_eigenvectors=False, is_mesh_symmetry=True,
                    is_gamma_center=True)
    else:
        ph.run_mesh(mesh=60.0, with_eigenvectors=False, is_mesh_symmetry=True)
    mesh_freqs = np.asarray(ph.get_mesh_dict()["frequencies"], dtype=float)
    if mesh_freqs.size == 0 or not np.all(np.isfinite(mesh_freqs)):
        sys.exit("[ERROR] mesh frequencies 为空或含非有限值")
    mf = float(np.min(mesh_freqs))
    mx = float(np.max(mesh_freqs))
    # ---- 2D：用 kz=0 的二维高对称路径（P2-1）；seekpath 是 3D 工具，会给 Γ-A 这类
    #   真空方向的平凡线段，写死的 Γ-M-K-Γ 又只对六方成立。----
    try:
        if is2d:
            from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
            _paths, _labels, _lat = kmc.band_path_2d_prim(ph, 101, vac_axis or 2)
            _bands, _conn = get_band_qpoints_and_path_connections(_paths, npoints=101)
            ph.run_band_structure(_bands, path_connections=_conn, labels=_labels)
            ph.write_yaml_band_structure(filename="band-dft-cpu.yaml")
            print("[..] band-dft-cpu.yaml 用 2D 高对称路径（%s 格子，kz=0）" % _lat)
        else:
            ph.auto_band_structure(plot=False, write_yaml=True, filename="band-dft-cpu.yaml")
    except Exception as e:
        sys.exit("[ERROR] band-dft-cpu.yaml 生成失败：%s" % e)

    band_text = (cwd / "band-dft-cpu.yaml").read_text(encoding="utf-8", errors="replace")
    band_values = [float(x) for x in re.findall(r"^\s*frequency:\s*([^\s#]+)", band_text, re.MULTILINE)]
    band_freqs = np.asarray(band_values, dtype=float)
    if band_freqs.size == 0 or not np.all(np.isfinite(band_freqs)):
        sys.exit("[ERROR] band frequencies empty or nonfinite")
    bmf = float(np.min(band_freqs))
    bmx = float(np.max(band_freqs))
    screen_min = min(mf, bmf)
    # ---- 2D 的 ZA 二次性（P2-2/P2-4）：symfc 的 fc2 没有施加 Born-Huang 旋转不变
    #   约束，ZA 弯曲支可能被线性化成 ω∝q —— 那时频率全为正、虚频闸门照样过，但 2D 的
    #   κ / 稳定性结论是错的。这里把 p 指数记进 summary，不合格直接判不稳定。----
    za = None
    if is2d:
        print("[WARN] 2D + symfc fc2：没有施加 Born-Huang 旋转不变约束，ZA 可能被线性化；"
              "2D 的稳定性/κ 请用 kl-dft-cpu 的 pheasy+BHH 链（PHEASY_RASR=auto）。")
        za = kmc.za_check_2d(ph, uc.cell, vac_axis or 2)
        print("[..] ZA 检查：%s" % za.get("note", za.get("error", "")))

    # ---- 虚频判据统一走 imag_policy（唯一真源，2026-09-22）；阈值读 mlff conf 否则用公共默认 ----
    import imag_policy

    def _imag_samples(_ph):
        _f, _q = [], []
        try:
            _md = _ph.get_mesh_dict()
            _f += [[float(x) for x in row] for row in _md["frequencies"]]
            _q += [[float(x) for x in q] for q in _md["qpoints"]]
        except Exception:                          # noqa: BLE001
            pass
        try:
            _bd = _ph.get_band_structure_dict()
            _f += [[float(x) for x in x2]
                   for seg in (_bd.get("frequencies") or []) for x2 in seg]
            _q += [[float(x) for x in x2]
                   for seg in (_bd.get("qpoints") or []) for x2 in seg]
        except Exception:                          # noqa: BLE001
            pass
        return _f, _q

    _fr, _qp = _imag_samples(ph)
    try:
        _vaxp = kmc.vacuum_axis_in_primitive(ph, vac_axis or 2) if is2d else None
    except Exception:                              # noqa: BLE001
        _vaxp = None
    _imag = imag_policy.classify_imag(
        _fr, _qp, is2d, _vaxp,
        {"IMAG_THR": params.get("IMAG_THR"),
         "IMAG_THR_STRICT": params.get("IMAG_THR_STRICT"),
         "IMAG_QGAMMA": params.get("IMAG_QGAMMA")})
    _za_v = imag_policy.za_verdict(za)
    stability_verdict = imag_policy.combine_verdict(_imag["verdict"], _za_v)
    stable = imag_policy.is_stable(stability_verdict)
    stable_mesh = (_imag["verdict"] != "fail")
    print("[..] imag_policy=%s/%s  min=%.4f THz (%.2f cm-1)  ZA=%s  最终=%s"
          % (_imag["verdict"], _imag["imag_class"], _imag["min_freq_THz"] or 0.0,
             _imag["min_freq_cm1"] or 0.0, _za_v, stability_verdict))

    # Chemistry-specific quality flag only; not proof of a failed fit.
    bound = None
    raw_bound = params.get("MAX_FREQ_MIN_THZ")
    if raw_bound not in (None, ""):
        try:
            bound = float(raw_bound)
        except (TypeError, ValueError):
            sys.exit("[ERROR] MAX_FREQ_MIN_THZ must be finite and positive")
        if not np.isfinite(bound) or bound <= 0:
            sys.exit("[ERROR] MAX_FREQ_MIN_THZ must be finite and positive")
    fit_suspect = bound is not None and max(mx, bmx) < bound
    note = "imag_policy=%s/%s，合并采样 min=%.4f THz (%.2f cm-1)，最终 stability_verdict=%s" % (
        _imag["verdict"], _imag["imag_class"], _imag["min_freq_THz"] or 0.0,
        _imag["min_freq_cm1"] or 0.0, stability_verdict)
    if fit_suspect:
        note += "；最高频率低于项目复核阈值，仅标可疑，不证明拟合失败"
        print("[WARN] " + note)

    try:
        fc_cutoff_A = float(fc_cutoff)
    except (TypeError, ValueError):
        fc_cutoff_A = None

    # 位移帧数是否足以确定 fc2：方程数 = n_disp × 3N_sc，未知数 = nfree_fc2(ALM)
    n_sc = int(len(uc.numbers) * np.prod(dim))
    n_disp = int(len(disps))
    dof = 3 * n_sc
    nfree_fc2 = _nfree_fc2_alm(ph, fc_cutoff_A) if fc_cutoff_A else None
    if nfree_fc2 is None or nfree_fc2 <= 0:
        frame_margin = None
        fit_underdetermined = None
    else:
        frame_margin = round((n_disp * dof) / float(nfree_fc2), 3)
        fit_underdetermined = bool(n_disp * dof < nfree_fc2)
    if fit_underdetermined:
        note += "；位移帧数不足（n_disp×3N < 2阶力常数元数），拟合欠定，虚频不可信"
        print("[WARN] " + note)

    if za is not None and za.get("ok") is False:
        note += "；ZA 弯曲支非二次色散（symfc 没有旋转不变约束），2D 结论不可信"
    summary = {
        "PHONON_DONE": True,
        "dim": str(dimtag).upper(),
        "za_exponent": za,
        "stable": bool(stable),
        "stable_mesh": bool(stable_mesh),
        "mesh_min_frequency_THz": mf,
        "mesh_max_frequency_THz": mx,
        "band_min_frequency_THz": bmf,
        "band_max_frequency_THz": bmx,
        "screening_min_frequency_THz": screen_min,
        "screening_valid": True,
        "screening_status": "valid",
        "stable_new": bool(stable),
        # ---- imag_policy（唯一真源）----
        "stability_verdict": stability_verdict,
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
        "za_verdict": _za_v,
        "criterion_version": "phonon-screen-v3-imag_policy",
        "criterion": "finite sampled mesh and band paths; not a whole-BZ proof",
        "min_frequency_THz": screen_min,
        "max_frequency_THz": max(mx, bmx),
        "fit_suspect": bool(fit_suspect),
        "max_freq_lower_bound_THz": bound,
        "assessment_policy": "stable 由 imag_policy.classify_imag（mesh+band 合并；近 Γ 声学支 warn、其余 fail）与 ZA 结论中更坏者决定；可选 max-frequency bound 只标 fit_suspect；n_disp×3N_sc < nfree_fc2(ALM) 标 fit_underdetermined。阈值见 skill/_common/imag_policy.py。",
        "n_disp": n_disp,
        "supercell": " ".join(str(x) for x in dim),
        "supercell_atoms": n_sc,
        "nfree_fc2_est": nfree_fc2,
        "fit_frame_margin": frame_margin,
        "fit_underdetermined": fit_underdetermined,
        "fc_cutoff_A": fc_cutoff_A,
        "cutoff_convergence_limitation": "FC_CUTOFF 是计算/建模截断；r_max 不证明 fc2 严格为零，需做 cutoff 收敛性检查。",
        "note": note,
    }
    (cwd / "phonon_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    print("[DONE] stable=%s min_freq=%.4f THz max_freq=%.3f THz fit_suspect=%s fit_underdetermined=%s"
          % (str(stable).lower(), mf, mx, str(fit_suspect).lower(),
             str(fit_underdetermined).lower()))


if __name__ == "__main__":
    main()