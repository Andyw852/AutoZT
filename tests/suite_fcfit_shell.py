#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""suite_fcfit_shell -- fit-fc-thermal 的「壳层 / 三阶截断确定」单元自测。

覆盖（从 kl-dft-cpu 移植过来、保留同样的判据语义）：
  1. fc_common.neighbor_shells / cut3_candidates：壳层枚举 + 相邻壳层中点
  2. fc_common.resolve_cut3_candidates：auto / off / 显式列表
  3. fc_common.supercell_safe_cutoff：周期镜像安全截断
  4. fc_fit_driver.fc3_shell_index / bins / means：逐壳层 |Phi^3|
  5. fc_fit_driver.fc3_shell_stats_from_means / _recommend_cut：σ/|mean| 与
     stable_upper_cut（数据能确定的最大截断）
  6. fc_fit_driver._detect_coords_mode / _normalize_displacement_frames：
     COORDS=auto 自动分辨分数坐标 vs 笛卡尔位移，平衡帧可在文件任意位置
     （首/中/末），无精确理想帧时不丢帧并用最近帧力作平衡近似

不碰集群、不依赖 pheasy / phono3py。
"""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "skill" / "fit-fc-thermal"))

import fc_common as fc  # noqa: E402

FAIL = []
N = 0


def check(name, cond, extra=""):
    global N
    N += 1
    if cond:
        print("  PASS %s" % name)
    else:
        print("  FAIL %s %s" % (name, extra))
        FAIL.append(name)


def _write_poscar(d, a=3.0):
    p = Path(d) / "POSCAR"
    p.write_text(
        "Si\n1.0\n%.6f 0 0\n0 %.6f 0\n0 0 %.6f\nSi\n1\nDirect\n0 0 0\n"
        % (a, a, a), encoding="utf-8")
    return p


def test_shells_and_candidates():
    print("[1] 壳层枚举 + 候选中点")
    with tempfile.TemporaryDirectory() as d:
        poscar = _write_poscar(d)
        shells, nat = fc.neighbor_shells(poscar, r_max=7.0, tol=0.05)
        check("neighbor_shells 非空且 nat=1", len(shells) >= 3 and nat == 1,
              str((shells, nat)))
        if len(shells) >= 3:
            check("最近邻壳层 = 3.0 A", abs(shells[0] - 3.0) < 1e-6, str(shells[:3]))
            check("壳层严格升序",
                  all(shells[i] < shells[i + 1] for i in range(len(shells) - 1)))
        cands, sh, note = fc.cut3_candidates(poscar, r_max=7.0, tol=0.05)
        check("cut3_candidates 有候选", len(cands) >= 2, str(cands))
        ok_mid = True
        for c in cands:
            below = [s for s in sh if s < c]
            above = [s for s in sh if s > c]
            if not below or not above:
                ok_mid = False
                break
            if abs(c - round(0.5 * (below[-1] + above[0]), 2)) > 1e-6:
                ok_mid = False
                break
            if any(abs(c - s) < 1e-6 for s in sh):
                ok_mid = False
                break
        check("候选都是相邻壳层中点、不落在壳层上", ok_mid, str((cands, sh)))


def test_resolve():
    print("[2] resolve_cut3_candidates 解析")
    with tempfile.TemporaryDirectory() as d:
        poscar = _write_poscar(d)
        c_off, _, _ = fc.resolve_cut3_candidates(poscar, "off")
        check("off -> 空", c_off == [], str(c_off))
        c_none, _, _ = fc.resolve_cut3_candidates(poscar, "")
        check("空串 -> 空（不自动）", c_none == [])
        c_auto, _, _ = fc.resolve_cut3_candidates(poscar, "auto")
        check("auto -> 与 cut3_candidates 一致",
              c_auto == fc.cut3_candidates(poscar)[0])
        c_list, _, _ = fc.resolve_cut3_candidates(poscar, "4.8 5.3 5.9")
        check("显式列表按序", c_list == [4.8, 5.3, 5.9], str(c_list))


def test_safe_cutoff():
    print("[3] 超胞安全截断")
    with tempfile.TemporaryDirectory() as d:
        poscar = _write_poscar(d, a=10.0)
        lat, frac = fc.read_poscar_cell_frac(poscar)
        safe = fc.supercell_safe_cutoff(lat, margin=0.1)
        check("10 A 立方胞安全截断 = 4.9 A", abs(safe - 4.9) < 1e-6, str(safe))


def test_fc3_shell_stability():
    print("[4/5] 逐壳层 |Phi^3| 与 stable_upper_cut")
    import numpy as np
    import fc_fit_driver as drv
    cell = np.diag([3.0, 3.0, 3.0])
    frac = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0],
                     [0.0, 0.5, 0.5]])
    shell_dist, tri = drv.fc3_shell_index(cell, frac, tol=0.05, r_max=8.0)
    check("fc3_shell_index 返回壳层且升序", len(shell_dist) >= 1
          and all(shell_dist[i] <= shell_dist[i + 1]
                  for i in range(len(shell_dist) - 1)), str(shell_dist[:4]))
    check("tri 形状正确", tri.shape == (4, 4, 4), str(tri.shape))
    bins = drv.fc3_shell_bins(tri, len(shell_dist))
    check("每个壳层都有 bin 索引", len(bins) == len(shell_dist))
    rng = np.random.default_rng(0)
    base = rng.normal(size=(4, 4, 4, 3, 3, 3))
    fits = [base + 0.02 * rng.normal(size=base.shape) for _ in range(6)]
    per_means = [drv.fc3_shell_means(f, bins) for f in fits]
    stats, upper = drv.fc3_shell_stats_from_means(shell_dist, per_means,
                                                  stability_thr=0.3)
    check("低噪声 -> stable_upper_cut 有限", upper > 0 and upper != float("inf"),
          str(upper))
    check("rel_std 有定义", all(s["rel_std"] is not None for s in stats), str(stats[:2]))

    # high scatter: bootstrap fits that do not agree -> the cutoff the data can
    # determine collapses to (at most) the nearest shell
    noisy = [base + (0.5 + i) * rng.normal(size=base.shape) for i in range(6)]
    pn = [drv.fc3_shell_means(f, bins) for f in noisy]
    _, upper_n = drv.fc3_shell_stats_from_means(shell_dist, pn, stability_thr=0.3)
    check("高噪声 -> stable_upper_cut 收到最近邻附近", upper_n <= shell_dist[0] + 0.05,
          str(upper_n))

    recs = [{"cut": 4.0, "stable_upper_cut": 5.0},
            {"cut": 5.0, "stable_upper_cut": 5.0},
            {"cut": 6.0, "stable_upper_cut": 5.0}]
    check("_recommend_cut 取自身稳定范围内的最大候选", drv._recommend_cut(recs) == 5.0)
    check("_recommend_cut 全不可用时 None",
          drv._recommend_cut([{"cut": 6.0, "stable_upper_cut": 4.0}]) is None)


def test_prep_coords_equilibrium():
    print("[6] prep 坐标类型自动判定 + 平衡帧任意位置识别")
    import numpy as np
    import fc_fit_driver as drv
    rng = np.random.default_rng(3)
    cell = np.diag([9.0, 9.0, 15.0])
    N = 12

    # fractional: the ideal frame sits in the MIDDLE of the file
    frac_ref = rng.random((N, 3))
    fr = np.stack([frac_ref + rng.normal(0, 0.001, (N, 3)) for _ in range(9)])
    fr[4] = frac_ref.copy()
    ff = rng.normal(0, 1.0, (9, N, 3))
    mode, dev, _ = drv._detect_coords_mode(fr, cell, frac_ref)
    check("分数坐标 -> auto 判定 fractional", mode == "fractional", mode)
    check("理想帧在中间 -> 最近帧距离 ~0", dev < 1e-4, str(dev))
    d, _f, eq, info = drv._normalize_displacement_frames(fr, ff, cell, frac_ref, mode)
    check("中间平衡帧被识别且剔除", info["reference_frame_index"] == 4
          and len(d) == 8, str((info["reference_frame_index"], len(d))))
    check("平衡力取自该中间帧", eq is not None and np.allclose(eq, ff[4]))
    check("平衡帧已从训练集剔除（无零位移帧）",
          float(np.sqrt((d ** 2).sum(axis=2).mean(axis=1)).min()) > 1e-4)

    # fractional with no exact ideal frame: keep every frame, eq from nearest
    fr2 = np.stack([frac_ref + rng.normal(0, 0.003, (N, 3)) for _ in range(9)])
    ff2 = rng.normal(0, 1.0, (9, N, 3))
    d2, _f2, eq2, i2 = drv._normalize_displacement_frames(fr2, ff2, cell,
                                                          frac_ref, "fractional")
    check("无精确理想帧 -> 不丢帧", len(d2) == 9, str(len(d2)))
    check("无精确理想帧 -> 用最近帧的力当平衡近似", eq2 is not None,
          str(i2.get("equilibrium_source")))

    # Cartesian displacements: the zero frame sits in the MIDDLE
    cad = rng.normal(0, 0.05, (11, N, 3))
    cad[6] = 0.0
    fc3 = rng.normal(0, 1.0, (11, N, 3))
    mode3, dev3, _ = drv._detect_coords_mode(cad, cell, frac_ref)
    check("笛卡尔位移 -> auto 判定 cartesian", mode3 == "cartesian", mode3)
    check("笛卡尔位移与理想位置相差数埃", dev3 > 1.0, str(dev3))
    d3, _f3, eq3, i3 = drv._normalize_displacement_frames(cad, fc3, None, None,
                                                        "cartesian")
    check("中间零位移帧被识别且剔除", i3["reference_frame_index"] == 6
          and len(d3) == 10, str(i3.get("reference_frame_index")))
    check("零位移帧的力作为平衡参考", eq3 is not None and np.allclose(eq3, fc3[6]))

    # Cartesian with no near-zero frame: keep every frame, no equilibrium force
    cad2 = rng.normal(0, 0.5, (11, N, 3))
    d4, _f4, eq4, _i4 = drv._normalize_displacement_frames(cad2, fc3, None, None,
                                                          "cartesian")
    check("无近零位移帧 -> 不丢帧且无平衡力", len(d4) == 11 and eq4 is None,
          str((len(d4), eq4 is not None)))

    # ...but a frame within 0.25 A of zero is still used as the equilibrium proxy
    cad3 = rng.normal(0, 0.02, (11, N, 3))
    _d5, _f5, eq5, i5 = drv._normalize_displacement_frames(cad3, fc3, None, None,
                                                          "cartesian")
    check("最近帧 <0.25 A -> 用作平衡近似", eq5 is not None,
          str(i5.get("equilibrium_source")))

def test_kl_bundle_reuse():
    print("[7] fit-fc-thermal kl_bundle + kl reuse_fcfit 布局")
    import importlib.util
    import fc_fit_driver as drv
    with tempfile.TemporaryDirectory() as d:
        src = Path(d) / "step1_fit"
        src.mkdir()
        for f in ("fc2.hdf5", "fc3.hdf5", "phonon_summary.json",
                  "cutoff_scan.json", "fc_dataset.json", "phono3py_disp.yaml"):
            (src / f).write_bytes(b"x")
        (src / "shengbte").mkdir()
        (src / "shengbte" / "FORCE_CONSTANTS_2ND").write_bytes(b"y")
        drv._write_kl_bundle(src, {})
        b = src / "kl_bundle"
        check("kl_bundle 生成 fc2/fc3",
              (b / "fc2.hdf5").is_file() and (b / "fc3.hdf5").is_file())
        check("kl_bundle 带 shengbte 力常数",
              (b / "shengbte" / "FORCE_CONSTANTS_2ND").is_file())
        try:
            spec = importlib.util.spec_from_file_location(
                "reuse_fcfit", str(ROOT / "skill" / "kl-dft-cpu" / "reuse_fcfit.py"))
            rf = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(rf)
        except Exception as e:  # noqa: BLE001
            print("  SKIP reuse_fcfit 导入失败：%s" % e)
            return
        out = Path(d) / "step5_fc"
        rf.assemble(b, out, supercell="3 3 1", dim="3d", mesh="24 24 24")
        check("reuse_fcfit 产出 phono3py/fc2.hdf5",
              (out / "phono3py" / "fc2.hdf5").is_file())
        check("reuse_fcfit 写 kl_params.txt（SUPERCELL 折叠成 3 3 1）",
              (out / "kl_params.txt").is_file()
              and "SUPERCELL=3 3 1" in (out / "kl_params.txt").read_text())
        check("reuse_fcfit 复制 shengbte",
              (out / "shengbte" / "FORCE_CONSTANTS_2ND").is_file())

def test_kappa_step():
    print("[8] step2_kappa 步骤与 κ 抽取")
    import importlib.util
    import numpy as np
    base = ROOT / "skill" / "fit-fc-thermal"
    text = (base / "skill.yaml").read_text(encoding="utf-8")
    check("skill.yaml 有 step2_kappa 且 needs step1_fit",
          "step2_kappa" in text and "needs: [step1_fit]" in text)
    check("skill.yaml 声明 KAPPA_DONE marker", "KAPPA_DONE" in text)
    for f in ("gen_step2_kappa.py", "kappa_driver.py",
              "templates/step2_kappa/submit_kappa.tpl",
              "templates/step2_kappa/step.conf"):
        check("文件存在 %s" % f, (base / f).is_file())
    spec = importlib.util.spec_from_file_location(
        "kd", str(base / "kappa_driver.py"))
    kd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kd)
    check("κ (T,6) Voigt -> xx/yy/zz",
          kd._flat_kappa(np.arange(6.0).reshape(1, 6)).tolist() == [[0.0, 1.0, 2.0]])
    check("κ (T,3,3) -> 对角",
          kd._flat_kappa(np.arange(9.0).reshape(1, 3, 3)).tolist()
          == [[0.0, 4.0, 8.0]])
    check("phono3py 逐档扫描函数存在（_phono3py_scan）",
          "_phono3py_scan" in (base / "fc_fit_driver.py").read_text(encoding="utf-8"))
    check("每档拟合都出声子谱（fc_plot_phonon 有 cut3_bands）",
          "cut3_bands" in (base / "fc_plot_phonon.py").read_text(encoding="utf-8"))

def test_cut3_select():
    print("[9] cut3_select.select_cutoff 三判据")
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "cut3select", str(ROOT / "skill" / "fit-fc-thermal" / "cut3_select.py"))
    cs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cs)
    recs = [
        {"cut": 3.9, "ratio": 8, "train_rel_err": 0.0062, "train_rel_se": 0.0004,
         "err_kind": "in-sample", "stable_upper_cut": 4.5, "kappa": 3.68, "kappa_err": 0.01},
        {"cut": 4.3, "ratio": 8, "train_rel_err": 0.0063, "train_rel_se": 0.0004,
         "err_kind": "in-sample", "stable_upper_cut": 4.5, "kappa": 3.68, "kappa_err": 0.01},
        {"cut": 4.9, "ratio": 8, "train_rel_err": 0.0064, "train_rel_se": 0.0004,
         "err_kind": "in-sample", "stable_upper_cut": 4.5, "kappa": 2.0, "kappa_err": 0.01},
    ]
    chosen, rep = cs.select_cutoff(recs)
    check("选截断：稳定的最小平台档 (3.9)", chosen == 3.9, str((chosen, rep.get("status"))))
    by = {x["cut"]: x for x in rep["records"]}
    check("超 stable_upper_cut 的档 stable_ok=False", by[4.9]["stable_ok"] is False)
    _c2, rep2 = cs.select_cutoff([dict(recs[0], kappa=3.0), dict(recs[1], kappa=2.0)])
    check("κ 未平台 -> not_converged", rep2["status"] == "not_converged", rep2["status"])
    c3, rep3 = cs.select_cutoff([dict(recs[0], ratio=1.0)])
    check("数据不足 -> no_usable_cutoff",
          c3 is None and rep3["status"] == "no_usable_cutoff", str(rep3["status"]))

def main():
    test_shells_and_candidates()
    test_resolve()
    test_safe_cutoff()
    test_fc3_shell_stability()
    test_prep_coords_equilibrium()
    test_kl_bundle_reuse()
    test_kappa_step()
    test_cut3_select()
    print("\nsuite_fcfit_shell: %s（%d 项）"
          % ("ALL PASS" if not FAIL else "FAIL", N))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
