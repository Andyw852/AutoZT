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

def test_kappa_mesh_auto():
    print("[10] step2_kappa q 网格 auto + 收敛测试")
    import importlib.util
    import json
    import types
    import numpy as np
    base = ROOT / "skill" / "fit-fc-thermal"
    spec = importlib.util.spec_from_file_location(
        "kd_mesh", str(base / "kappa_driver.py"))
    kd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kd)
    check("MESH=auto -> 长度 MESH_LENGTH",
          kd._mesh_spec({"mesh": "auto", "mesh_length": 40}) == ("length", 40.0))
    check("MESH=60 -> 长度 60", kd._mesh_spec({"mesh": "60"}) == ("length", 60.0))
    check("MESH=n n n -> 显式",
          kd._mesh_spec({"mesh": "8 8 4"}) == ("explicit", [8, 8, 4]))

    # 真空轴 = 倒格矢最短的原胞轴（不是看实空间最长边）
    class _Prim:
        cell = np.array([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 20.0]])
    ph = types.SimpleNamespace(primitive=_Prim())
    check("2D 真空轴取倒格矢最短轴", kd._vacuum_axis(ph) == 2)

    # 逐分量判据：面内 xx/yy 与面外 zz 各自用自己的 κ 归一
    # （旧判据 max|Δκ|/max|κ_ii| 会把慢收敛的 zz 掩盖掉）
    a = {"temperatures": [300.0], "kappa_voigt_xx_yy_zz_yz_xz_xy": [[2, 2, 1, 0, 0, 0]]}
    b = {"temperatures": [300.0], "kappa_voigt_xx_yy_zz_yz_xz_xy": [[2, 2, 1, 0, 0, 0.1]]}
    rel, rel_vec, _t = kd._mesh_change(a, b, 300)
    check("收敛判据含非对角元 (0.1/2=5%)", abs(rel - 0.05) < 1e-12, str(rel))
    check("逐分量 rel_vec=[xx,yy,zz]", len(rel_vec) == 3, str(rel_vec))
    a = {"temperatures": [300.0], "kappa_voigt_xx_yy_zz_yz_xz_xy": [[4, 4, 2, 0, 0, 0]]}
    b = {"temperatures": [300.0], "kappa_voigt_xx_yy_zz_yz_xz_xy": [[4, 4, 1.8, 0, 0, 0]]}
    rel_pc, rel_vec_pc, _t = kd._mesh_change(a, b, 300)
    rel_old, _rv, _t = kd._mesh_change(a, b, 300, mode="max_ii")
    check("zz 慢收敛不再被掩盖（逐分量 10% > 旧 5%）",
          abs(rel_pc - 0.10) < 1e-9 and rel_pc > rel_old + 0.04,
          "pc=%s old=%s" % (rel_pc, rel_old))
    # 分母下限：zz 只有 xx 的 2.5%（2D 真空轴 / 极弱层间）时，它自身 20% 的相对
    # 噪声不再卡住收敛——按 0.1×max κ_ii 归一：0.02/0.4 = 5%；floor=0 仍是严格 20%
    a = {"temperatures": [300.0], "kappa_voigt_xx_yy_zz_yz_xz_xy": [[4, 4, 0.1, 0, 0, 0]]}
    b = {"temperatures": [300.0], "kappa_voigt_xx_yy_zz_yz_xz_xy": [[4, 4, 0.08, 0, 0, 0]]}
    rel_f, _rv, _t = kd._mesh_change(a, b, 300)
    rel_s, _rv, _t = kd._mesh_change(a, b, 300, floor=0.0)
    check("小分量按下限归一（0.02/0.4=5%，严格逐分量 20%）",
          abs(rel_f - 0.05) < 1e-9 and abs(rel_s - 0.20) < 1e-9,
          "floor=%s strict=%s" % (rel_f, rel_s))

    # 模拟 BTE：κ = 1 + 1/n，逐档加密直到相邻变化 < tol，取较密那档
    calls = []

    def fake_apply(_ph3, sp, _cfg):
        if sp[0] == "explicit":
            return list(sp[1])
        n = max(1, int(round(sp[1] * 0.3)))
        return [n, n, n]

    def fake_run(_out, _y, _f2, _f3, cfg, _src, _w, spec=None):
        n = fake_apply(None, spec, cfg)[0]
        calls.append(n)
        k = 1.0 + 1.0 / n
        return {"temperatures": [300.0], "mesh": "%d %d %d" % (n, n, n),
                "kappa_xx_yy_zz": [[k, k, k]],
                "kappa_voigt_xx_yy_zz_yz_xz_xy": [[k, k, k, 0, 0, 0]]}

    fake_p3 = types.ModuleType("phono3py")
    fake_p3.load = lambda *a, **k: object()
    old_mod = sys.modules.get("phono3py")
    sys.modules["phono3py"] = fake_p3
    kd._apply_mesh, kd._run_bte = fake_apply, fake_run
    try:
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            cfg = {"mesh": "auto", "mesh_length": 20, "mesh_conv": "auto",
                   "mesh_conv_tol_pct": 3.0, "mesh_conv_factor": 1.25,
                   "mesh_conv_max_length": 300, "mesh_conv_max_points": 10 ** 9,
                   "mesh_conv_t": 300}
            s = kd._converge_mesh(out, "y", "fc2", "fc3", cfg, "src")
            rep = json.loads((out / "mesh_convergence.json").read_text())
            check("网格收敛：converged 且取最后（较密）一档",
                  s["mesh_converged"] and s["mesh"] == rep["records"][-1]["mesh"],
                  str(rep))
            check("网格收敛：相邻变化 < 3%",
                  rep["records"][-1]["rel_change"] < 0.03, str(rep["records"][-1]))
            check("网格收敛：不重复跑同一网格", len(calls) == len(set(calls)), str(calls))
            calls.clear()
            cfg2 = dict(cfg, mesh_conv_max_length=30)
            s2 = kd._converge_mesh(out, "y", "fc2", "fc3", cfg2, "src")
            check("触顶未收敛 -> mesh_converged=False",
                  s2["mesh_converged"] is False, str(s2.get("mesh_convergence")))
            calls.clear()
            s3 = kd._converge_mesh(out, "y", "fc2", "fc3",
                                   dict(cfg, mesh="8 8 8"), "src")
            check("显式网格不做收敛（只跑一次）", len(calls) == 1 and
                  "mesh_converged" not in s3, str(calls))
    finally:
        if old_mod is None:
            sys.modules.pop("phono3py", None)
        else:
            sys.modules["phono3py"] = old_mod

def test_compact_fc3_expand():
    print("[11] 紧凑 fc3 展开（旧版 phono3py 回退路径 == phono3py 官方）")
    try:
        import numpy as np
        from phono3py import Phono3py
        from phono3py.phonon3.fc3 import compact_fc3_to_full_fc3
        from phonopy.structure.atoms import PhonopyAtoms
    except Exception as e:  # noqa: BLE001
        print("  SKIP 需要 phono3py>=4：%s" % e)
        return
    import fc_fit_driver as drv
    uc = PhonopyAtoms(symbols=["Na", "Cl"], cell=np.eye(3) * 4.0,
                      scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]])
    ph3 = Phono3py(uc, supercell_matrix=[2, 2, 2], primitive_matrix="P")
    n = len(ph3.supercell)
    rng = np.random.default_rng(0)
    c = rng.normal(size=(len(ph3.primitive), n, n, 3, 3, 3))
    ref = compact_fc3_to_full_fc3(ph3.primitive, c)
    import builtins
    real_import = builtins.__import__

    def _no_helper(name, *a, **k):
        if name == "phono3py.phonon3.fc3":
            raise ImportError("simulated old phono3py")
        return real_import(name, *a, **k)

    builtins.__import__ = _no_helper
    try:
        got = drv._expand_compact_fc3(c, ph3.primitive)
    finally:
        builtins.__import__ = real_import
    check("回退展开与 compact_fc3_to_full_fc3 一致",
          got.shape == ref.shape and np.allclose(got, ref))


def test_method_sweep():
    print("[12] 方法扫描：第 N 近邻截断（thirdorder.py -n）+ SWEEP_* 解析")
    import importlib
    sys.path.insert(0, str(ROOT / "skill" / "_common" / "opt"))
    with tempfile.TemporaryDirectory() as td:
        # 简单立方 a=2：近邻 2, 2√2, 2√3, 4 -> N=1 截断 (2+2.828)/2
        p = Path(td) / "POSCAR"
        p.write_text("sc\n1.0\n2 0 0\n0 2 0\n0 0 2\nX\n1\nDirect\n0 0 0\n")
        cuts, per_atom = fc.nth_neighbor_cutoffs(p, 3)
        import math
        check("简单立方 N=1..3 截断 = 相邻近邻距离中点",
              all(abs(c - e) < 1e-3 for c, e in zip(
                  cuts, [(2 + 2 * math.sqrt(2)) / 2,
                         (2 * math.sqrt(2) + 2 * math.sqrt(3)) / 2,
                         (2 * math.sqrt(3) + 4) / 2])), str(cuts))
        # 两种原子、近邻距离不同：取原子间最大（thirdorder calc_frange 的 max(tonth)）
        p.write_text("cscl\n1.0\n3 0 0\n0 3 0\n0 0 3.3\nA B\n1 1\nDirect\n"
                     "0 0 0\n0.5 0.5 0.5\n")
        cuts, per_atom = fc.nth_neighbor_cutoffs(p, 2)
        exp = max(0.5 * (u[1] + u[2]) for u in per_atom)
        check("N=2 截断取各原子 (u2+u3)/2 的最大值", abs(cuts[1] - exp) < 1e-4,
              "%s vs %s" % (cuts, exp))
    g3 = importlib.import_module("gen_step3_sweep")
    got = g3.parse_methods("phono3py pheasy:ols,Ridge hiphive:all")
    check("SWEEP_METHODS 解析（默认方法 / 大小写 / engine:all）",
          got == [("phono3py", "symfc"), ("pheasy", "OLS"), ("pheasy", "RIDGE")]
          + [("hiphive", m) for m in g3.ENGINE_METHODS["hiphive"]], str(got))
    check("SWEEP_METHODS=all = 全部 13 种",
          len(g3.parse_methods("all")) == 2 + 6 + 5)
    check("SWEEP_C3_SHELLS 解析", g3.parse_shells("3-5 7,8") == [3, 4, 5, 7, 8]
          and g3.parse_shells("auto") is None)
    try:
        g3.parse_shells("0")
        ok = False
    except SystemExit:
        ok = True
    check("SWEEP_C3_SHELLS=0 被拒绝", ok)
    o = g3.variant_overrides("pheasy", "RIDGE", 5.2761, None)
    check("变体覆盖：c3 写入、单变体截断扫描关闭",
          o["PHEASY_C3_CUTOFF"] == "5.2761" and o["PHEASY_C2_CUTOFF"] == ""
          and o["CUT3_SCAN"] == "off" and o["FIT_ENGINE"] == "pheasy", str(o))


def test_cut3_without_bootstrap_and_curve():
    print("[13] 无 bootstrap 时只按 κ 平台选截断 + kappa_vs_cutoff 输出")
    import importlib.util
    import json
    base = ROOT / "skill" / "fit-fc-thermal"
    sys.path.insert(0, str(base))
    import cut3_select as cs
    recs = [{"cut": c, "ratio": 10, "stable_upper_cut": None,
             "stability_measured": False, "kappa": k}
            for c, k in ((4.2, 31.9), (4.7, 28.9), (5.2, 27.0), (5.9, 24.6), (6.5, 23.6))]
    ch, rep = cs.select_cutoff(recs)
    check("bootstrap=0：跳过判据②，按 κ 平台选到 5.9", ch == 5.9 and rep["status"] == "ok"
          and rep["stability_gate_applied"] is False, "%s %s" % (ch, rep["status"]))
    for r in recs:
        r["stability_measured"] = True
    ch2, rep2 = cs.select_cutoff(recs)
    check("有 bootstrap 但 stable_upper_cut 为空：仍无可用截断",
          ch2 is None and rep2["status"] == "no_usable_cutoff")
    # 逐分量平台：面内已平、zz 还在降 10% → 不能判收敛；2D 的 zz≈0 不卡平台
    vrec = [{"cut": c, "ratio": 10, "stability_measured": False, "kappa": k,
             "kappa_vec": [k, k, z]}
            for c, k, z in ((5.2, 27.0, 6.0), (5.9, 24.6, 5.0), (6.5, 23.6, 4.5))]
    ch4, rep4 = cs.select_cutoff(vrec)
    check("逐分量平台：zz 未平（-10%）→ not_converged",
          rep4["status"] == "not_converged" and ch4 == 6.5, "%s %s" % (ch4, rep4["status"]))
    for r in vrec:
        r["kappa_vec"][2] = 0.0
    ch5, rep5 = cs.select_cutoff(vrec)
    check("2D 的 κzz≈0 不卡平台（按下限归一）", ch5 == 5.9 and rep5["status"] == "ok",
          "%s %s" % (ch5, rep5["status"]))
    spec = importlib.util.spec_from_file_location("kd_curve", str(base / "kappa_driver.py"))
    kd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kd)
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        (out / "cutoff_scan.json").write_text(json.dumps({
            "bootstrap": 0, "shells": [3.0, 4.0, 5.0],
            "records": [{"cut": c, "ratio": 10, "stable_upper_cut": None}
                        for c in (3.5, 4.5, 5.5)]}))
        per = [(("%.2f" % c).replace(".", "p"),
                {"kappa_inplane_300K": k, "kappa_300K_xx_yy_zz": [k, k, 0.0]})
               for c, k in ((3.5, 30.0), (4.5, 25.0), (5.5, 24.5))]
        cwd = os.getcwd()
        os.chdir(td)
        try:
            ch3, rep3, _ = kd._select(out, {}, per)
        finally:
            os.chdir(cwd)
        doc = kd._write_kappa_vs_cutoff(out, per, rep3)
        check("旧 cutoff_scan.json（bootstrap=0）也能选截断", ch3 == 4.5, str(ch3))
        check("kappa_vs_cutoff.json 含每档 κ 与壳层数",
              (out / "kappa_vs_cutoff.json").is_file()
              and [r["n_shells"] for r in doc["rows"]] == [1, 2, 3])
        check("kappa_vs_cutoff.json 含相邻两档的逐分量变化百分比",
              doc["rows"][1]["pct_change_from_prev"]["kappa_xx"] == round(-500 / 30, 3)
              and doc["rows"][1]["pct_change_from_prev"]["kappa_zz"] is None,
              str(doc["rows"][1].get("pct_change_from_prev")))


def test_kappa_2d_thickness_norm():
    print("[14] 2D 面内 κ 层厚归一化（与 kl-dft-cpu 同口径 h⊥/d）")
    import importlib.util
    base = ROOT / "skill" / "fit-fc-thermal"
    for d in (ROOT / "skill" / "_common", ROOT / "skill" / "_common" / "opt"):
        sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("g2_norm", str(base / "gen_step2_kappa.py"))
    g2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g2)
    with tempfile.TemporaryDirectory() as td:
        # MoS2 单层：a=3.14，真空胞 c=25，S 层 ±1.565 Å → d = 3.13 + 2×1.80(S vdW)
        Path(td, "POSCAR").write_text(
            "MoS2\n1.0\n3.14 0 0\n-1.57 2.719330 0\n0 0 25\nMo S\n1 2\nDirect\n"
            "0.333333 0.666667 0.5\n0.666667 0.333333 0.5626\n0.666667 0.333333 0.4374\n")
        g = g2.two_d_norm(td, "vdw")
        check("MoS2 层厚 d≈6.73 Å，factor=h⊥/d≈3.71",
              g and abs(g["thickness_d_A"] - 6.73) < 0.02
              and abs(g["kappa_2d_norm_factor"] - 25 / g["thickness_d_A"]) < 1e-3, str(g))
        check("KAPPA_2D_THICKNESS=cell 不归一", g2.two_d_norm(td, "cell") is None)
    spec = importlib.util.spec_from_file_location("kd_norm", str(base / "kappa_driver.py"))
    kd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kd)
    f = kd._two_d_fields({"kappa_2d_norm": {"kappa_2d_norm_factor": 4.0,
                                            "h_perp_A": 25.0}},
                         [[10.0, 12.0, 0.0], [20.0, 22.0, 0.0]], 1)
    check("kappa_2d_normalized_*：原始 κ × factor，原始字段不动",
          f["kappa_2d_normalized_inplane_300K"] == 84.0
          and f["kappa_2d_normalized_xx_yy_zz"][0] == [40.0, 48.0, 0.0])
    check("面热导 = 原始 κ × h⊥（W/K，与层厚取法无关）",
          abs(f["sheet_conductance_inplane_300K_W_per_K"] - 21.0 * 25e-10) < 1e-18,
          str(f.get("sheet_conductance_inplane_300K_W_per_K")))
    check("3D 不加归一化字段", kd._two_d_fields({}, [[1, 1, 1]], 0) == {})


def main():
    test_shells_and_candidates()
    test_resolve()
    test_safe_cutoff()
    test_fc3_shell_stability()
    test_prep_coords_equilibrium()
    test_kl_bundle_reuse()
    test_kappa_step()
    test_cut3_select()
    test_kappa_mesh_auto()
    test_compact_fc3_expand()
    test_method_sweep()
    test_cut3_without_bootstrap_and_curve()
    test_kappa_2d_thickness_norm()
    print("\nsuite_fcfit_shell: %s（%d 项）"
          % ("ALL PASS" if not FAIL else "FAIL", N))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
