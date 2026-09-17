#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""phonon_fit_driver.py —— 2 阶力常数拟合 + 声子谱（phonon-mace-cpu S3）。

计算节点作业里跑，cwd = step3_phonon/：读 POSCAR + disps.npy + forces.npy →
symfc 拟合 fc2 → q-mesh 最小频率（虚频闸）+ band-dft-cpu.yaml → phonon_summary.json。
"""
import json
import sys
from pathlib import Path

import numpy as np


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
    params = read_params(cwd / "klmace_params.txt")
    _dim = np.array([int(x) for x in (params.get("SUPERCELL") or "1 1 1").split()],
                    dtype=int)
    # 3 个数=对角扩胞；9 个数=3×3 矩阵（行主序，与 phonopy/phono3py --dim 同义）
    dim = _dim.reshape(3, 3) if _dim.size == 9 else np.diag(_dim)
    ph = Phonopy(uc, supercell_matrix=np.array(dim, dtype=int))

    disps = np.load("disps.npy")
    forces = np.load("forces.npy")
    ph.displacements = disps
    ph.forces = forces
    print("[..] %d 帧 × %d 原子，symfc 拟合 fc2" % (len(disps), len(uc.numbers) * int(np.prod(dim))))

    ph.produce_force_constants(fc_calculator="symfc")
    print("[OK] fc2 拟合完成")

    ph.run_mesh(mesh=60.0, with_eigenvectors=False, is_mesh_symmetry=True)
    mf = float(np.min(ph.get_mesh_dict()["frequencies"]))
    # ---- 2D：走 kz=0 的二维高对称路径（P2-1），并做 ZA 二次性检查（P2-2/P2-4）----
    import klmace_common as kmc
    dimtag, vac_axis = "3d", 2
    try:
        from dim_common import detect_dimension
        dimtag, vac_axis, _ = detect_dimension("POSCAR")
    except Exception as e:
        print("[WARN] 维度判定失败（按 3D 处理）：%s" % e)
    is2d = str(dimtag).lower().startswith("2")
    za = None
    try:
        if is2d:
            from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
            paths, labels, lat = kmc.band_path_2d(uc.cell, 101, vac_axis or 2)
            bands, conn = get_band_qpoints_and_path_connections(paths, npoints=101)
            ph.run_band_structure(bands, path_connections=conn, labels=labels)
            ph.write_yaml_band_structure(filename="band-dft-cpu.yaml")
            print("[..] band-dft-cpu.yaml 用 2D 高对称路径（%s 格子，kz=0）" % lat)
        else:
            ph.auto_band_structure(plot=False, write_yaml=True, filename="band-dft-cpu.yaml")
    except Exception as e:
        print("[WARN] band-dft-cpu.yaml 出图跳过：%s" % e)
    if is2d:
        # symfc 的 fc2 没有施加 Born-Huang 旋转不变/零应力约束：2D 的 ZA 会被线性化。
        # 这里把 p 指数记下来（判据 1.7<p<2.3），并明确提示 2D 该用 DFT 链。
        print("[WARN] 2D + symfc fc2：没有施加 Born-Huang 旋转不变约束，ZA 弯曲支可能被"
              "线性化（ω∝q 而非 q²）。要做 2D 的稳定性/κ，请用 kl-dft-cpu 的 pheasy+BHH 链"
              "（PHEASY_RASR=auto）。")
        za = kmc.za_check_2d(ph, uc.cell, vac_axis or 2)
        print("[..] ZA 检查：%s" % za.get("note", za.get("error", "")))

    stable = mf >= -0.10
    if za is not None and za.get("ok") is False:
        stable = False
    summary = {
        "PHONON_DONE": True,
        "stable": bool(stable),
        "min_frequency_THz": mf,
        "n_disp": int(len(disps)),
        "supercell": " ".join(str(x) for x in dim),
        "dim": str(dimtag).upper(),
        "za_exponent": za,
        "note": "最小声子频率 %.4f THz（阈值 -0.10）：%s"
                % (mf, "无明显虚频" if stable else
                   ("ZA 弯曲支非二次色散（symfc 无旋转约束），2D 结果不可信"
                    if (za and za.get("ok") is False) else "存在虚频，动力学不稳定")),
    }
    (cwd / "phonon_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    print("[DONE] stable=%s min_freq=%.4f THz" % (str(stable).lower(), mf))


if __name__ == "__main__":
    main()
