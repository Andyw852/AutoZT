#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tests/test_supercell_symmetry.py —— 超胞保对称 / 力常数对称 / q 网格逐方向收敛（2026-10-07）。

起因：Mn2In2Se5（R-3m）在菱方原胞上取 3x3x1 超胞，超胞晶格只剩 12 个点群操作里的
4 个。拟合出的 fc2 破对称：布里渊区内部出现假虚频，kappa 随 q 网格来回跳；S1 虚频
门禁（粗网格 + 对称约化）没看出来；S2 的 _converge_mesh 拿相邻两档比，[1,14,42] 与
[1,18,54] 的 k_z 分格一样，2.4% 就判了收敛。本文件用合成的 R-3m 结构（不用任何
真实材料数据）覆盖四处修复：
  1. symmetry_audit.supercell_symmetry_gate（生成位移 / S1 拟合前的超胞闸）
  2. fc_fit_driver：门禁 fc2 原子顺序对齐 + 对称等价 q 点频差检查
  3. kappa_driver._supercell_keeps_symmetry（S2 BTE 前的超胞闸）
  4. kappa_driver._converge_mesh：只和"每个方向都更粗"的那档比
"""
import json
import math
import os
import sys
import types

import pytest

np = pytest.importorskip("numpy")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _d in ("skill/_common/fcfit", "skill/_common", "skill/_common/opt", ""):
    sys.path.insert(0, os.path.join(ROOT, _d))
sys.dont_write_bytecode = True

import symmetry_audit as SA      # noqa: E402
import kappa_driver as kd        # noqa: E402

A_LAT = 4.0     # 面内晶格常数 (A)
D_LAYER = 10.0  # 层周期 (A)，六方 c = 3 * D_LAYER
T_Y = 2.5       # Y 原子到 X 原子沿三重轴的距离 (A)


def _lattice():
    """菱方原胞，基矢取法与 Mn2In2Se5 那份一致：a1、a2 在面内成 60 度，
    a3 = (a1 + a2)/3 + 层周期 * z。这种取法下 n n 1 超胞丢三重轴。"""
    a1 = np.array([A_LAT, 0.0, 0.0])
    a2 = np.array([A_LAT / 2, A_LAT * math.sqrt(3) / 2, 0.0])
    a3 = (a1 + a2) / 3 + np.array([0.0, 0.0, D_LAYER])
    return np.array([a1, a2, a3])


def _write_poscar(path, x_frac_one=False):
    """X 在原点、两个 Y 在三重轴上 ±T_Y：空间群 R-3m (166)。
    x_frac_one=True 时 X 写成 (1,0,0)——复现 POSCAR 里坐标正好是 1.0 的情形。"""
    lat = _lattice()
    inv = np.linalg.inv(lat)
    fy1 = np.array([0.0, 0.0, T_Y]) @ inv
    fy2 = np.array([0.0, 0.0, -T_Y]) @ inv
    fx = [1.0, 0.0, 0.0] if x_frac_one else [0.0, 0.0, 0.0]
    lines = ["synthetic R-3m", "1.0"]
    lines += ["  %.10f %.10f %.10f" % tuple(v) for v in lat]
    lines += ["Bi Se", "1 2", "Direct"]
    for f in (fx, fy1 % 1.0, fy2 % 1.0):
        lines.append("  %.10f %.10f %.10f" % tuple(f))
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------
# 1. 生成阶段的超胞闸（symmetry_audit）
# --------------------------------------------------------------------------
def test_supercell_report_and_suggestion(tmp_path):
    pytest.importorskip("spglib")
    st = SA.read_poscar(_write_poscar(tmp_path / "POSCAR"))
    rots, sg = SA.point_group_rotations(st["lattice"], st["frac"], st["symbols"])
    assert sg == "R-3m (166)" and len(rots) == 12

    bad = SA.supercell_symmetry_report(st["lattice"], st["frac"], st["symbols"], [3, 3, 1])
    assert not bad["ok"] and (bad["n_kept"], bad["n_ops"]) == (4, 12)
    assert bad["natoms"] == 27
    # 9 个数的写法与对角写法等价
    bad9 = SA.supercell_symmetry_report(st["lattice"], st["frac"], st["symbols"],
                                        [3, 0, 0, 0, 3, 0, 0, 0, 1])
    assert bad9["n_kept"] == 4
    for reps in ([3, 3, 3], [2, 2, 2], [1, 1, 2]):
        assert SA.supercell_symmetry_report(st["lattice"], st["frac"], st["symbols"],
                                            reps)["ok"], reps

    sug = SA.suggest_symmetric_diag(st["lattice"], st["frac"], st["symbols"], [3, 3, 1])
    assert sug["reps"] == [3, 3, 3] and sug["natoms"] == 81


def test_supercell_gate_modes(tmp_path, capsys):
    pytest.importorskip("spglib")
    pos = _write_poscar(tmp_path / "POSCAR")
    with pytest.raises(SystemExit) as e:
        SA.supercell_symmetry_gate(pos, [3, 3, 1], mode="strict")
    assert "4/12" in str(e.value) and 'SUPERCELL="3 3 3"' in str(e.value)

    rep = SA.supercell_symmetry_gate(pos, [3, 3, 1], mode="warn")
    assert rep is not None and not rep["ok"] and rep["suggestion"]["reps"] == [3, 3, 3]
    assert "[WARN]" in capsys.readouterr().out

    assert SA.supercell_symmetry_gate(pos, [3, 3, 1], mode="off") is None
    assert SA.supercell_symmetry_gate(pos, [3, 3, 3], mode="strict")["ok"]
    with pytest.raises(SystemExit):
        SA.supercell_symmetry_gate(pos, [3, 3, 3], mode="sometimes")


def test_supercell_gate_without_structure_only_warns(tmp_path, capsys):
    # 读不了结构 / 没有 spglib 时只告警跳过，不阻断生成
    assert SA.supercell_symmetry_gate(tmp_path / "missing", [3, 3, 1]) is None
    assert "跳过" in capsys.readouterr().out


# --------------------------------------------------------------------------
# 2. S2（kappa_driver）BTE 前的超胞闸
# --------------------------------------------------------------------------
def test_kappa_driver_supercell_check(tmp_path):
    spglib = pytest.importorskip("spglib")
    st = SA.read_poscar(_write_poscar(tmp_path / "POSCAR"))
    lat = np.asarray(st["lattice"], float)
    cell = (lat, np.asarray(st["frac"], float), SA.symbol_to_z(st["symbols"]))
    rots = list({np.asarray(r).tobytes(): np.asarray(r)
                 for r in spglib.get_symmetry(cell, symprec=1e-5)["rotations"]}.values())
    for reps, want in (([3, 3, 1], (4, 12)), ([3, 3, 3], (12, 12)), ([2, 2, 2], (12, 12))):
        sc = np.diag(reps) @ lat
        assert kd._supercell_keeps_symmetry(lat, rots, sc) == want, reps


def test_kappa_driver_check_modes(monkeypatch):
    """strict 停、warn 照算、off 不查（用假的 phono3py 对象）。"""
    pytest.importorskip("spglib")
    lat = _lattice()
    inv = np.linalg.inv(lat)
    frac = np.array([[0, 0, 0], np.array([0, 0, T_Y]) @ inv % 1.0,
                     np.array([0, 0, -T_Y]) @ inv % 1.0])
    prim = types.SimpleNamespace(cell=lat, scaled_positions=frac, numbers=[83, 34, 34])
    ph3 = types.SimpleNamespace(primitive=prim,
                                supercell=types.SimpleNamespace(cell=np.diag([3, 3, 1]) @ lat),
                                phonon_supercell=None)
    kd._SC_SYM_CHECKED.clear()
    with pytest.raises(SystemExit) as e:
        kd._check_supercell_symmetry(ph3, {})
    assert "4/12" in str(e.value)
    kd._SC_SYM_CHECKED.clear()
    kd._check_supercell_symmetry(ph3, {"supercell_symmetry": "warn"})
    kd._SC_SYM_CHECKED.clear()
    kd._check_supercell_symmetry(ph3, {"supercell_symmetry": "off"})
    assert not kd._SC_SYM_CHECKED          # off 根本不查


# --------------------------------------------------------------------------
# 3. q 网格收敛：每个方向都加密过才判
# --------------------------------------------------------------------------
def test_mesh_axis_counts():
    # 广义网格：grid_matrix 每行 gcd = 惯用胞该轴分格数（Mn2In2Se5 的 [1,14,42]）
    s = {"mesh": "1 14 42", "mesh_grid_matrix": [[0, -14, 0], [-14, 14, 0], [1, 1, -3]]}
    assert kd._mesh_axis_counts(s) == [14, 14, 1]
    s = {"mesh": "2 22 66", "mesh_grid_matrix": [[0, -22, 0], [-22, 22, 0], [2, 2, -6]]}
    assert kd._mesh_axis_counts(s) == [22, 22, 2]
    # 对角网格：就是 D_diag
    assert kd._mesh_axis_counts({"mesh": "8 8 3", "mesh_grid_matrix": None}) == [8, 8, 3]
    assert kd._mesh_axis_counts({"mesh": "8 8 3",
                                 "mesh_grid_matrix": [[8, 0, 0], [0, 8, 0], [0, 0, 3]]}) \
        == [8, 8, 3]


def test_reference_mesh_is_coarser_everywhere():
    h = [{"mesh_axis_counts": [14, 14, 1]}, {"mesh_axis_counts": [18, 18, 1]}]
    # k_z 没变：没有各方向都更粗的参照档 -> 不判
    assert kd._ref_coarser_everywhere(h[:1], {"mesh_axis_counts": [18, 18, 1]}) is None
    # k_z 进位：参照最近一档各方向都更粗的 [18,18,1]
    assert kd._ref_coarser_everywhere(h, {"mesh_axis_counts": [22, 22, 2]}) is h[1]
    h2 = h + [{"mesh_axis_counts": [22, 22, 2]}]
    assert kd._ref_coarser_everywhere(h2, {"mesh_axis_counts": [28, 28, 2]}) is h[1]
    # 冻结轴（2D 真空轴）不参与比较
    h3 = [{"mesh_axis_counts": [20, 20, 1]}]
    assert kd._ref_coarser_everywhere(h3, {"mesh_axis_counts": [25, 25, 1]},
                                      frozen=(2,)) is h3[0]
    # 分格数拿不到：退回旧行为（上一档）
    assert kd._ref_coarser_everywhere([{"mesh": "x"}], {"mesh": "y"}) == {"mesh": "x"}


def test_converge_mesh_waits_for_every_axis(tmp_path, monkeypatch):
    """面内收敛快、k_z 收敛慢的假 BTE：旧判据在 [18,18,1] 就会判收敛（k_z 没测）。"""
    astar, cstar = 0.2864, 0.02054         # Mn2In2Se5 惯用胞倒格矢长度 (1/A)

    def counts(L):
        return [max(1, round(L * astar)), max(1, round(L * astar)),
                max(1, round(L * cstar))]

    def fake_apply(_ph3, spec, _cfg):
        return counts(spec[1])

    def kappa(c):
        n, _, m = c
        return 1.0 + 0.3 / n + 0.15 / m, 0.13 + 0.03 / m

    def fake_run(_o, _y, _f2, _f3, cfg, _src, _w, spec=None):
        c = counts(spec[1])
        kx, kz = kappa(c)
        return {"temperatures": [300.0], "mesh": " ".join(map(str, c)),
                "mesh_grid_matrix": None, "mesh_axis_counts": c,
                "kappa_xx_yy_zz": [[kx, kx, kz]],
                "kappa_voigt_xx_yy_zz_yz_xz_xy": [[kx, kx, kz, 0, 0, 0]]}

    # 旧判据的陷阱确实存在：[14,14,1] -> [18,18,1] 变化不到 3%
    a, b = kappa([14, 14, 1]), kappa([18, 18, 1])
    assert abs(b[0] - a[0]) / a[0] < 0.03 and b[1] == a[1]

    fake_p3 = types.ModuleType("phono3py")
    fake_p3.load = lambda *a, **k: object()
    monkeypatch.setitem(sys.modules, "phono3py", fake_p3)
    monkeypatch.setattr(kd, "_apply_mesh", fake_apply)
    monkeypatch.setattr(kd, "_run_bte", fake_run)
    monkeypatch.setattr(kd, "_write_kappa_vs_mesh", lambda *a, **k: None)
    cfg = {"mesh": "auto", "mesh_length": 50, "mesh_conv": "auto",
           "mesh_conv_tol_pct": 3.0, "mesh_conv_factor": 1.25,
           "mesh_conv_max_length": 300, "mesh_conv_max_points": 10 ** 9,
           "mesh_conv_t": 300}
    s = kd._converge_mesh(tmp_path, "y", "fc2", "fc3", cfg, "src")
    rep = json.loads((tmp_path / "mesh_convergence.json").read_text())
    assert s["mesh_converged"] and rep["frozen_axes"] == []
    assert s["mesh_axis_counts"][2] >= 3, s["mesh"]       # k_z 必须加密过
    for r in rep["records"]:
        if r.get("compared_with"):
            ref = [int(x) for x in r["compared_with"].split()]
            assert all(x < y for x, y in zip(ref, r["axis_counts"])), r
        else:
            assert r["rel_change"] is None


# --------------------------------------------------------------------------
# 4. S1 门禁：fc2 原子顺序对齐 + 对称等价 q 点频差
# --------------------------------------------------------------------------
def _spring_fc2(sc, rc=4.2, k=2.0):
    """超胞上的中心力弹簧模型 fc2（最小镜像，rc 小于超胞 WS 内切半径）。"""
    pos, lat = np.asarray(sc.positions, float), np.asarray(sc.cell, float)
    shifts = np.array([[i, j, l] for i in (-1, 0, 1) for j in (-1, 0, 1)
                       for l in (-1, 0, 1)]) @ lat
    n = len(pos)
    fc = np.zeros((n, n, 3, 3))
    for i in range(n):
        d = pos[None, :] - pos[i] + shifts[:, None, :]          # (27, n, 3)
        r = np.linalg.norm(d, axis=-1)
        r[:, i] = np.inf
        best = np.argmin(r, axis=0)
        for j in range(n):
            rr = r[best[j], j]
            if rr < rc:
                e = d[best[j], j] / rr
                fc[i, j] = -k * np.outer(e, e)
        fc[i, i] = -fc[i].sum(axis=0)
    return fc


def test_gate_fc2_order_and_symmetry(tmp_path):
    pytest.importorskip("spglib")
    pytest.importorskip("ase")
    phonopy = pytest.importorskip("phonopy")
    from phonopy.interface.vasp import read_vasp, write_vasp
    import fc_fit_driver as fd

    pos = _write_poscar(tmp_path / "POSCAR", x_frac_one=True)   # X 写成 (1,0,0)
    # 3x3x1：和 Mn2In2Se5 同款超胞。弹簧截断 4.2 A 小于其 WS 内切半径（5.1 A），
    # 所以弹簧 fc2 本身满足三重对称（超胞保不保对称与此无关）
    scm = np.diag([3, 3, 1])
    # 数据集按"未卷回"的 POSCAR 生成超胞（真实数据集就是这样）
    ds = phonopy.Phonopy(read_vasp(str(pos)), supercell_matrix=scm).supercell
    write_vasp(str(tmp_path / "SPOSCAR"), ds)
    fc2 = _spring_fc2(ds)

    uc = fd._phonopy_unitcell(pos)                # 门禁读法：坐标卷回 [0,1)
    gate_sc = phonopy.Phonopy(uc, supercell_matrix=scm).supercell
    aligned = fd._fc2_in_phonopy_order(fc2, uc, scm, None, tmp_path / "SPOSCAR")
    # 弹簧 fc2 只由位置决定：在门禁的超胞顺序上直接算一遍，必须和重排结果一致；
    # 不重排的原 fc2 则对不上（X 的 9 个像换了位）
    want = _spring_fc2(gate_sc)
    assert np.allclose(aligned, want) and not np.allclose(fc2, want)

    def spread(f):
        ph = phonopy.Phonopy(uc, supercell_matrix=scm)
        ph.force_constants = f
        return fd._fc2_symmetry_spread(ph, {})

    good = spread(aligned)
    assert good["n_ops"] == 12 and good["spread_THz"] < 1e-6, good

    # 破三重对称的 fc2：频差远超 FC2_SYM_TOL
    bad = aligned.copy()
    sc_pos = np.asarray(gate_sc.positions)
    for j in range(len(sc_pos)):
        d = sc_pos[j] - sc_pos[0]
        if np.linalg.norm(d) < 4.2 and d[0] > 1.0:
            bad[0, j] *= 1.3
            bad[j, 0] *= 1.3
    for i in range(len(bad)):
        bad[i, i] = 0.0
        bad[i, i] = -bad[i].sum(axis=0)
    assert spread(bad)["spread_THz"] > fd.FC2_SYM_TOL
