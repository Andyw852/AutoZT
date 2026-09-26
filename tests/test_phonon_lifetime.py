# -*- coding: utf-8 -*-
"""phonon-lifetime 技能回归。

约定：tau[ps] = 1/(2*2*pi*gamma_total[THz])，与 phono3py 官方后处理一致
（scripts/phono3py_kdeplot.py 的 tau = 1/g/(2*2*pi)、other/kaccum.py 的 mfp）。
覆盖三类：换算约定、加权统计/选文件/对称方向等纯函数、gen 端到端。
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

ROOT = Path(__file__).resolve().parent.parent
SKILL = ROOT / "skill" / "phonon-lifetime"
OPT = ROOT / "skill" / "_common" / "opt"
sys.path.insert(0, str(SKILL))
sys.path.insert(0, str(OPT))

LC = pytest.importorskip("lifetime_common")

# 真实 Si 参考数据（本机有才跑，集群 CI 没有就 skip）
# 设 AUTOZT_SI_REF=<Si 的 kl-dft-cpu/result/step6_kappa 目录> 才跑；不设则 skip（不写死本机路径）
SI_REF = Path(os.environ.get("AUTOZT_SI_REF", "/nonexistent/AUTOZT_SI_REF"))


# ---------------------------------------------------------------------------
# 约定 / 纯函数
# ---------------------------------------------------------------------------
def test_tau_convention_matches_phono3py():
    gamma = np.array([1.0, 0.05, 0.0, -0.2])
    tau = LC.gamma_to_tau_ps(gamma)
    assert tau[0] == pytest.approx(1.0 / (2 * 2 * np.pi))
    assert tau[1] == pytest.approx(1.0 / (0.05 * 2 * 2 * np.pi))
    assert np.isinf(tau[2]) and np.isinf(tau[3])
    ref = 1.0 / np.where(gamma > 0, gamma, -1) / (2 * 2 * np.pi)
    assert tau[0] == pytest.approx(ref[0], rel=1e-15)


def test_tau_gamma_roundtrip():
    g = np.array([0.0, 1e-4, 0.02, 0.1])
    tau = LC.gamma_to_tau_ps(g)
    back = LC.tau_to_gamma_thz(tau)
    mask = g > 0
    assert np.allclose(back[mask], g[mask], rtol=1e-12)


def test_mean_free_path_matches_kaccum():
    gv = np.zeros((1, 1, 3))
    gv[0, 0] = [0.0, 0.0, 2.0]
    g = np.array([[1.0]])
    mfp = LC.mean_free_path_ang(gv, g)
    assert mfp[0, 0] == pytest.approx(2.0 / (2 * 2 * np.pi))


def test_effective_gamma_adds_isotope_only_when_positive():
    g = np.array([[1.0, 2.0]])
    assert np.allclose(LC.effective_gamma(g, None), g)
    assert np.allclose(LC.effective_gamma(g, np.zeros((1, 2))), g)
    assert np.allclose(LC.effective_gamma(g, np.full((1, 2), 0.5)), g + 0.5)


def test_boundary_gamma_formula_micrometer():
    """C1：boundary_mfp 单位是微米（phono3py help 原文 "in micrometer"）。"""
    gv = np.zeros((1, 1, 3))
    gv[0, 0] = [0.0, 0.0, 4.0]
    gbd = LC.boundary_gamma(gv, 1.0)              # L = 1 µm = 1e4 Å
    assert gbd[0, 0] == pytest.approx(4.0 / (4 * np.pi * 1e4))
    assert gbd[0, 0] == pytest.approx(4.0 / (4 * np.pi * 1.0 * LC.MICRON_TO_ANGSTROM))
    # 0.1 µm（100 nm）比 1 µm 强 10 倍
    assert LC.boundary_gamma(gv, 0.1)[0, 0] == pytest.approx(10.0 * gbd[0, 0])
    assert LC.boundary_gamma(gv, 0) is None
    assert LC.boundary_gamma(None, 100.0) is None


def test_boundary_gamma_decision_with_micrometer_mfp():
    """C1 回归：κ 由 γ_anh+γ_iso+γ_bd（微米口径）重建时必须判定“计入”且误差 ~0。

    旧代码漏了 1e4 换算，γ_bd 偏大 1e4 倍 -> 计入后 κ 误差更大 -> 静默判定“不计入”，
    本用例会红。
    """
    ngp, nb = 2, 2
    gv2 = np.zeros((ngp, nb, 6))
    gv2[..., 0] = 1.0
    cv = np.ones((1, ngp, nb))
    uc, mesh = np.array(2.0), np.array([2, 2, 2])
    g_total = np.full((ngp, nb), 0.02) + np.full((ngp, nb), 0.01)
    gv = np.zeros((ngp, nb, 3))
    gv[..., 0] = 4.0                              # |v| = 4 THz·Å
    L_um = 0.1
    g_bd_true = gv[..., 0] / (4 * np.pi * L_um * 1e4)   # 显式微米公式
    k = LC.reconstruct_kappa(gv2, cv[0], g_total + g_bd_true, uc, mesh=mesh)
    data = {"gv_by_gv": gv2, "heat_capacity": cv, "kappa": k[None, :],
            "kappa_unit_conversion": uc, "mesh": mesh,
            "temperatures": np.array([300.0]), "boundary_mfp": np.array(L_um)}
    use, info = LC.decide_boundary_gamma(data, g_total, LC.boundary_gamma(gv, L_um),
                                         t_indices=[0])
    assert use is True
    assert info["boundary_mfp_um"] == pytest.approx(L_um)
    assert info["err_with_boundary"] < 1e-12
    assert info["err_no_boundary"] > 1e-3
    # 默认 1e6 µm 是 phono3py 的“关”，直接短路，不去比 κ
    data2 = dict(data)
    data2["boundary_mfp"] = np.array(1e6)
    use2, info2 = LC.decide_boundary_gamma(data2, g_total, LC.boundary_gamma(gv, 1e6),
                                           t_indices=[0])
    assert use2 is False and info2["method"] == "mfp_sentinel"


def test_pick_temperature_indices_warns_on_missing():
    Ts = np.array([100.0, 200.0, 300.0, 400.0])
    assert set(LC.pick_temperature_indices(Ts).values()) == {0, 1, 2, 3}
    got = LC.pick_temperature_indices(Ts, [290.0, 300.0])
    assert sorted(got.values()) == [2]
    assert 300.0 in got
    msgs = []
    got2 = LC.pick_temperature_indices(Ts, [1000.0], log=msgs.append)
    assert got2 == {400.0: 3}
    assert msgs and "1000" in msgs[0] and "400" in msgs[0]   # 不再静默取最近温度


def test_acoustic_mask_default_ignores_frequency_window():
    # 默认只按支编号取最低 3 支：X 点 ~12 THz 的 LA 必须仍算声学（旧版会被踢去光学）
    f = np.array([[1.0, 12.0, 13.0, 14.0]])
    g = np.full_like(f, 0.01)
    mask, fmax = LC.acoustic_mask(f, g)
    assert fmax is None
    assert list(mask[0]) == [True, True, True, False]
    # 显式给上限才启用频率窗
    mask2, fmax2 = LC.acoustic_mask(f, g, acoustic_fmax=5.0)
    assert fmax2 == 5.0
    assert list(mask2[0]) == [True, False, False, False]


def test_acoustic_mask_lowest_three_bands():
    f = np.array([[1.0, 2.0, 3.0, 10.0, 11.0, 12.0]])
    g = np.full_like(f, 0.01)
    mask, fmax = LC.acoustic_mask(f, g, acoustic_fmax=5.0)
    assert list(mask[0]) == [True, True, True, False, False, False]
    assert fmax == 5.0


# ---------------------------------------------------------------------------
# 加权统计（B4）
# ---------------------------------------------------------------------------
def test_summarize_temperature():
    f = np.array([[1.0, 2.0, 3.0, 10.0]])
    g = np.full((1, 1, 4), 0.01)
    s = LC.summarize_temperature(f, g, np.array([300.0]), 0, acoustic_fmax=5.0)
    assert s["T_K"] == 300.0
    assert s["acoustic"]["n_modes"] == 3
    assert s["optical"]["n_modes"] == 1
    assert s["acoustic"]["tau_median_ps"] == pytest.approx(1.0 / (0.01 * 2 * 2 * np.pi))
    assert s["weights_applied"] is False


def test_weighted_median_matches_full_grid_expansion():
    # 3 个不可约 q 点、1 条支、权重 1:1:47 —— 权重必须作用在 q 点维度上
    f = np.array([[1.0], [2.0], [3.0]])
    taus = np.array([[[1.0], [2.0], [100.0]]])
    g = 1.0 / (LC.TAU_FACTOR * taus)
    w = np.array([1.0, 1.0, 47.0])
    s = LC.summarize_temperature(f, g, np.array([300.0]), 0, weights=w)
    rep = np.repeat([1.0, 2.0, 100.0], [1, 1, 47])   # 展开到全网格
    assert s["weights_applied"] is True
    assert s["all"]["tau_median_ps"] == pytest.approx(float(np.median(rep)))
    assert s["all"]["tau_p90_ps"] == pytest.approx(float(np.percentile(rep, 90)))
    assert s["all"]["n_modes"] == 3
    assert s["all"]["n_modes_weighted"] == pytest.approx(49.0)
    # 不给权重 = 按不可约点等权（旧行为）
    s0 = LC.summarize_temperature(f, g, np.array([300.0]), 0)
    assert s0["all"]["tau_median_ps"] == pytest.approx(2.0)


def test_weighted_overdamped_fraction():
    f = np.array([[0.5, 2.0, 8.0]])
    g = np.array([[[5.0, 0.01, 0.01]]])      # band0 线宽巨大 -> omega*tau << 1
    s = LC.summarize_temperature(f, g, np.array([300.0]), 0, weights=np.array([1.0]))
    assert s["all"]["overdamped_frac"] == pytest.approx(1.0 / 3.0)


def test_isotope_shortens_lifetime():
    f = np.array([[1.0, 2.0, 3.0, 10.0]])
    g = np.full((1, 1, 4), 0.01)
    gi = np.full((1, 4), 0.01)
    s = LC.summarize_temperature(f, g, np.array([300.0]), 0,
                                 acoustic_fmax=5.0, gamma_isotope=gi)
    assert s["isotope_included"] is True
    assert s["all"]["tau_median_ps"] == pytest.approx(1.0 / (0.02 * 2 * 2 * np.pi))
    assert s["all"]["tau_median_anharmonic_ps"] == pytest.approx(
        1.0 / (0.01 * 2 * 2 * np.pi))


def test_gamma_extra_is_added_like_boundary():
    f = np.array([[1.0, 2.0]])
    g = np.full((1, 1, 2), 0.01)
    extra = np.array([[0.01, 0.01]])
    s = LC.summarize_temperature(f, g, np.array([300.0]), 0, gamma_extra=extra)
    assert s["boundary_included"] is True
    assert s["all"]["tau_median_ps"] == pytest.approx(1.0 / (0.02 * 2 * 2 * np.pi))


# ---------------------------------------------------------------------------
# 沿 q 方向（B6）
# ---------------------------------------------------------------------------
def test_select_direction_qpoints():
    q = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.25], [0.0, 0.0, 0.5],
                  [0.25, 0.0, 0.0], [0.5, 0.0, 0.0]])
    picked = LC.select_direction_qpoints(q, [0, 0, 1], n_points=2, tol=0.05)
    gp = [i for _, i in picked]
    assert 0 in gp and 2 in gp
    assert 3 not in gp and 4 not in gp


def _cubic_permutations():
    perms = []
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            k = 3 - i - j
            for a, b in ((j, k), (k, j)):
                r = np.zeros((3, 3))
                r[0, i] = 1
                r[1, a] = 1
                r[2, b] = 1
                perms.append(r)
    return np.array(perms)


def test_select_direction_qpoints_uses_symmetry():
    # 不可约楔子只存 (x,0,0)：不给对称操作时 [001] 只剩 Γ（旧 bug）
    q = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.5, 0.0, 0.0],
                  [0.25, 0.25, 0.0]])
    assert [i for _, i in LC.select_direction_qpoints(q, [0, 0, 1], n_points=2, tol=0.05)] == [0]
    rots = _cubic_permutations()
    picked = LC.select_direction_qpoints(q, [0, 0, 1], n_points=2, tol=0.05, rotations=rots)
    gp = [i for _, i in picked]
    assert 0 in gp and 1 in gp and 2 in gp
    assert 3 not in gp


def _full_grid_qpoints(n):
    """n^3 全网格分数坐标，故意把离线点排在前面（|qx|+|qy| 降序）。"""
    ii = np.arange(n) / float(n)
    q = np.stack(np.meshgrid(ii, ii, ii, indexing="ij"), axis=-1).reshape(-1, 3)
    q = q - np.round(q)
    order = np.argsort(-(np.abs(q[:, 0]) + np.abs(q[:, 1])))
    return q[order]


def test_select_direction_qpoints_tight_tolerance_on_dense_grid():
    """C2：31^3 格点间距 1/31≈0.032 < 旧容差 0.06，离线点会被当成线上点。

    离线点排在前面时旧代码会整列取错；新容差 1e-5 只认精确共线的点。
    """
    n = 31
    q = _full_grid_qpoints(n)
    assert abs(q[0, 0]) + abs(q[0, 1]) > 0.9        # 确实离线点在前
    picked = LC.select_direction_qpoints(q, [0, 0, 1], n_points=200)
    assert len(picked) == 16                        # t = 0, 1/31, ..., 15/31
    for t, gp in picked:
        assert abs(q[gp, 0]) < 1e-9 and abs(q[gp, 1]) < 1e-9
        # 线上点的 |z| = t（gp 可能是 -q 那一侧的等价代表）
        assert abs(q[gp, 2]) == pytest.approx(t)
    ts = [t for t, _ in picked]
    assert ts[0] == pytest.approx(0.0)
    assert ts[-1] == pytest.approx(15.0 / 31.0)


def test_select_direction_qpoints_uses_time_reversal():
    """C3：不可约楔里只存 -q 时，靠 τ(q)=τ(-q) 也要能找到（与有无反演中心无关）。"""
    q = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -0.25], [0.0, 0.0, -0.5],
                  [0.25, 0.25, 0.0]])
    # 不给任何对称操作：-q 仍应因时间反演被找到（旧代码这里只剩 Γ）
    picked = LC.select_direction_qpoints(q, [0, 0, 1], n_points=2, tol=1e-6)
    assert sorted(t for t, _ in picked) == pytest.approx([0.0, 0.25, 0.5])
    # 只给恒等操作（相当于没有反演中心）时同样
    ident = np.array([np.eye(3)])
    picked2 = LC.select_direction_qpoints(q, [0, 0, 1], n_points=2, tol=1e-6,
                                          rotations=ident)
    assert sorted(t for t, _ in picked2) == pytest.approx([0.0, 0.25, 0.5])
    # 离线点不能被带进来
    assert 3 not in [i for _, i in picked]
    assert 3 not in [i for _, i in picked2]


# ---------------------------------------------------------------------------
# 选文件（B1/B2）
# ---------------------------------------------------------------------------
def _write_fake_hdf5(path, mesh=(2, 2, 2), n_gp=4, n_band=3, weights=None,
                     with_isotope=False, gamma_is_2d=False):
    h5py = pytest.importorskip("h5py")
    T = np.array([100.0, 300.0])
    freq = np.linspace(0.5, 8.0, n_gp * n_band).reshape(n_gp, n_band)
    freq[0, 0] = 0.0
    if gamma_is_2d:
        gamma = np.full((n_gp, n_band), 0.01)
    else:
        gamma = np.full((2, n_gp, n_band), 0.01)
        gamma[:, 0, :] = 0.0
    with h5py.File(path, "w") as h:
        h["temperature"] = T
        h["frequency"] = freq
        h["gamma"] = gamma
        h["qpoint"] = np.zeros((n_gp, 3))
        h["weight"] = np.ones(n_gp) if weights is None else np.asarray(weights, float)
        h["mesh"] = np.array(mesh)
        h["group_velocity"] = np.ones((n_gp, n_band, 3))
        h["heat_capacity"] = np.ones((2, n_gp, n_band))
        h["kappa"] = np.tile(np.array([10.0, 10.0, 10.0, 0, 0, 0]), (2, 1))
        h["kappa_unit_conversion"] = np.array(1.0)
        if with_isotope:
            h["gamma_isotope"] = np.full((n_gp, n_band), 0.005)
    return freq


def test_select_kappa_hdf5_prefers_larger_mesh(tmp_path):
    pytest.importorskip("h5py")
    a = tmp_path / "kappa-m999.hdf5"
    _write_fake_hdf5(a, mesh=(9, 9, 9), n_gp=2)
    b = tmp_path / "kappa-m151515.hdf5"
    _write_fake_hdf5(b, mesh=(15, 15, 15), n_gp=4)
    c = tmp_path / "kappa-m151515-g0.hdf5"       # --write-gamma --gp 单点文件
    _write_fake_hdf5(c, mesh=(15, 15, 15), n_gp=1, gamma_is_2d=True)
    logs = []
    best, ignored = LC.select_kappa_hdf5([a, b, c], log=logs.append)
    assert best == b
    assert set(ignored) == {a, c}
    assert logs and "kappa-m151515.hdf5" in logs[0] and "忽略" in logs[0]
    # 只有单点文件时也要能用（不静默丢弃）
    best_only, _ = LC.select_kappa_hdf5([c])
    assert best_only == c


def test_find_kappa_hdf5_ignores_non_phono3py(tmp_path):
    pytest.importorskip("h5py")
    d = tmp_path / "mat" / "kl-dft-cpu" / "result" / "step6_kappa"
    d.mkdir(parents=True)
    _write_fake_hdf5(d / "kappa-m222.hdf5")
    (d / "fc2.hdf5").write_bytes(b"not hdf5")     # 垃圾文件必须被跳过
    skill = tmp_path / "mat" / "phonon-lifetime"
    skill.mkdir(parents=True)
    found = LC.find_kappa_hdf5(skill)
    assert found is not None and found.name == "kappa-m222.hdf5"


def test_find_kappa_hdf5_maps_kl_mlff_step(tmp_path):
    """C4：kl-mlff-* 的 κ 在 step4_kappa，不能写死 step6_kappa。"""
    pytest.importorskip("h5py")
    d = tmp_path / "mat" / "kl-mlff-cpu" / "step4_kappa"
    d.mkdir(parents=True)
    _write_fake_hdf5(d / "kappa-m444.hdf5", mesh=(4, 4, 4), n_gp=5)
    skill = tmp_path / "mat" / "phonon-lifetime"
    skill.mkdir(parents=True)
    found = LC.find_kappa_hdf5(skill, "kl-mlff-cpu", "step6_kappa")
    assert found is not None and found.name == "kappa-m444.hdf5"
    assert LC.upstream_steps("kl-mlff-cpu", "step6_kappa")[:2] == ["step6_kappa", "step4_kappa"]
    assert LC.upstream_steps("kl-dft-cpu", "step6_kappa")[:1] == ["step6_kappa"]
    # 步骤名完全对不上时，rglob 兜底也能找到
    d2 = tmp_path / "mat2" / "kl-mlff-cpu" / "step9_weird"
    d2.mkdir(parents=True)
    _write_fake_hdf5(d2 / "kappa-m333.hdf5", mesh=(3, 3, 3), n_gp=3)
    s2 = tmp_path / "mat2" / "phonon-lifetime"
    s2.mkdir(parents=True)
    assert LC.find_kappa_hdf5(s2, "kl-mlff-cpu", "step4_kappa").name == "kappa-m333.hdf5"


def test_fc_dataset_maps_kl_mlff_step(tmp_path):
    """C4：kl-mlff-* 的 fc2/fc3 在 step3_fc。"""
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "mat" / "kl-mlff-cpu" / "step3_fc" / "phono3py"
    d.mkdir(parents=True)
    for name in ("fc2.hdf5", "fc3.hdf5"):
        with h5py.File(str(d / name), "w") as h:
            h["fc"] = np.zeros((1,))
    (d / "phono3py_disp.yaml").write_text("dummy: 1\n", encoding="utf-8")
    skill = tmp_path / "mat" / "phonon-lifetime"
    skill.mkdir(parents=True)
    fc = LC.find_fc_dataset(skill, "kl-mlff-cpu", "step6_kappa")
    assert fc is not None and fc[1].name == "fc2.hdf5"


def test_select_kappa_hdf5_prefers_tetrahedron_over_smearing(tmp_path):
    """C5a：同 mesh 下展宽法 kappa-m*-s0.1.hdf5 让位给四面体法主文件，并告警。"""
    pytest.importorskip("h5py")
    tet = tmp_path / "kappa-m151515.hdf5"
    _write_fake_hdf5(tet, mesh=(15, 15, 15), n_gp=4)
    sme = tmp_path / "kappa-m151515-s0.1.hdf5"
    _write_fake_hdf5(sme, mesh=(15, 15, 15), n_gp=4)
    logs = []
    best, ignored = LC.select_kappa_hdf5([sme, tet], log=logs.append)
    assert best == tet
    assert ignored == [sme]
    assert any("kappa-m151515.hdf5" in m for m in logs)
    # 只有展宽文件时仍可用，但必须告警
    logs2 = []
    best2, _ = LC.select_kappa_hdf5([sme], log=logs2.append)
    assert best2 == sme
    assert logs2 and "展宽" in logs2[0]


# ---------------------------------------------------------------------------
# 边界散射 + κ 重建自检（B8）
# ---------------------------------------------------------------------------
def test_decide_boundary_gamma_uses_kappa_reconstruction():
    ngp, nb = 2, 2
    gv2 = np.zeros((ngp, nb, 6))
    gv2[..., 0] = 1.0
    cv = np.ones((1, ngp, nb))
    uc = np.array(2.0)
    mesh = np.array([2, 2, 2])
    g_anh = np.full((ngp, nb), 0.02)
    gi = np.full((ngp, nb), 0.01)
    g_total = g_anh + gi
    gv = np.zeros((ngp, nb, 3))
    gv[..., 0] = 4.0
    g_bd = LC.boundary_gamma(gv, 1.0)
    base = {"gv_by_gv": gv2, "heat_capacity": cv, "kappa_unit_conversion": uc,
            "mesh": mesh, "temperatures": np.array([300.0])}
    # 上游 kappa 含边界散射 -> 应当计入
    k_with = LC.reconstruct_kappa(gv2, cv[0], g_total + g_bd, uc, mesh=mesh)
    d1 = dict(base); d1["kappa"] = k_with[None, :]; d1["boundary_mfp"] = np.array(1.0)
    use1, info1 = LC.decide_boundary_gamma(d1, g_total, g_bd, t_indices=[0])
    assert use1 is True
    assert info1["err_with_boundary"] < info1["err_no_boundary"]
    # 上游 kappa 不含边界散射（boundary_mfp 只是 1e6 占位）-> 不能计入
    k_no = LC.reconstruct_kappa(gv2, cv[0], g_total, uc, mesh=mesh)
    d2 = dict(base); d2["kappa"] = k_no[None, :]; d2["boundary_mfp"] = np.array(1e6)
    use2, _ = LC.decide_boundary_gamma(d2, g_total, g_bd, t_indices=[0])
    assert use2 is False


def test_kappa_reconstruction_error_reports_small_residual():
    ngp, nb = 3, 2
    rng = np.random.default_rng(0)
    gv2 = rng.random((ngp, nb, 6))
    cv = np.ones((1, ngp, nb))
    gt = rng.random((ngp, nb)) + 0.1
    uc = 3.0
    mesh = np.array([3, 3, 3])
    k = LC.reconstruct_kappa(gv2, cv[0], gt, uc, mesh=mesh)
    data = {"gv_by_gv": gv2, "heat_capacity": cv, "kappa": k[None, :],
            "kappa_unit_conversion": uc, "mesh": mesh,
            "temperatures": np.array([300.0])}
    chk = LC.kappa_reconstruction_error(data, gt[None, ...], t_indices=[0])
    assert chk is not None and chk["max_rel_err"] < 1e-12
    # 口径不对（少了一部分线宽）就该报出大误差
    bad = LC.kappa_reconstruction_error(data, (gt * 2)[None, ...], t_indices=[0])
    assert bad["max_rel_err"] > 0.1


# ---------------------------------------------------------------------------
# fc 路线命令（B7）
# ---------------------------------------------------------------------------
def _run_gen_py(code, cwd=None):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(OPT) + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("LIFETIME_REEXEC", None)
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          env=env, timeout=300, cwd=cwd)


def test_phono3py_cmd_has_isotope_but_no_nac_flag():
    code = (
        "import sys; sys.path.insert(0, %r); sys.path.insert(0, %r);"
        "import gen_step1_lifetime as G;"
        "print(G._phono3py_cmd({'PHONO3PY_CMD': 'phono3py-load', 'INCLUDE_ISOTOPE': True},"
        " '3 3 3', '300'));"
        "print(G._phono3py_cmd({'PHONO3PY_CMD': 'phono3py-load', 'INCLUDE_ISOTOPE': False},"
        " '3 3 3', '300'))" % (str(SKILL), str(OPT))
    )
    r = _run_gen_py(code)
    assert r.returncode == 0, r.stdout + r.stderr
    lines = [ln for ln in r.stdout.splitlines() if "phono3py-load" in ln]
    assert len(lines) == 2, r.stdout
    assert "--isotope" in lines[0]
    assert '--ts="300"' in lines[0]
    assert "--mesh 3 3 3" in lines[0]
    assert "--isotope" not in lines[1]
    assert "--nac " not in lines[0]      # phono3py-load 没有 --nac，靠 BORN 自动开


def test_phono3py_run_copies_born(tmp_path):
    h5py = pytest.importorskip("h5py")
    src = tmp_path / "src"
    src.mkdir()
    for n in ("fc2.hdf5", "fc3.hdf5"):
        with h5py.File(str(src / n), "w") as h:
            h["fc"] = np.zeros((1,))
    (src / "phono3py_disp.yaml").write_text("dummy: 1", encoding="utf-8")
    (src / "BORN").write_text("14.0\n", encoding="utf-8")
    out = tmp_path / "out"
    out.mkdir()
    conf = {"PHONO3PY_MESH": "3 3 3", "PHONO3PY_TS": "300", "PHONO3PY_CMD": "true",
            "PHONO3PY_TIMEOUT": 60, "INCLUDE_ISOTOPE": True}
    code = ("import sys; sys.path.insert(0, %r); sys.path.insert(0, %r);"
            "import gen_step1_lifetime as G;"
            "G._run_phono3py(%r, %r, %r, (%r, %r, %r, %r))"
            % (str(SKILL), str(OPT), conf, str(tmp_path), str(out),
               str(src), str(src / "fc2.hdf5"), str(src / "fc3.hdf5"),
               str(src / "phono3py_disp.yaml")))
    # PHONO3PY_CMD=true 会“成功”但不产 hdf5 -> 预期 sys.exit(报找不到 hdf5)，但 BORN 必须已拷
    r = _run_gen_py(code)
    assert (out / "BORN").is_file()
    assert (out / "phono3py_disp.yaml").is_file()
    assert "没找到含 gamma 的 hdf5" in (r.stdout + r.stderr)


# ---------------------------------------------------------------------------
# gen 端到端（合成 kappa hdf5，不碰集群）
# ---------------------------------------------------------------------------
def _make_project(tmp_path, with_hdf5=True, with_fc=False, extra_conf="",
                  hdf5_name="kappa-m222.hdf5", weights=None,
                  upstream_skill="kl-dft-cpu", upstream_step="step6_kappa",
                  upstream_dir=None):
    skill = tmp_path / "Si" / "phonon-lifetime"
    skill.mkdir(parents=True)
    conf = (SKILL / "templates" / "step.conf").read_text(encoding="utf-8")
    conf = re.sub(r"(?m)^UPSTREAM_SKILL\s*=.*$", "UPSTREAM_SKILL = " + upstream_skill, conf)
    conf = re.sub(r"(?m)^UPSTREAM_STEP\s*=.*$", "UPSTREAM_STEP = " + upstream_step, conf)
    if extra_conf:
        conf = conf + "\n" + extra_conf
    (skill / "step.conf").write_text(conf, encoding="utf-8")
    if upstream_dir is not None:
        d = upstream_dir
    elif upstream_skill.startswith("kl-mlff"):
        d = tmp_path / "Si" / upstream_skill / "step4_kappa"
    else:
        d = tmp_path / "Si" / upstream_skill / "result" / upstream_step
    d.mkdir(parents=True)
    if with_hdf5:
        _write_fake_hdf5(d / hdf5_name, weights=weights)
    if with_fc:
        h5py = pytest.importorskip("h5py")
        for name in ("fc2.hdf5", "fc3.hdf5"):
            with h5py.File(str(d / name), "w") as h:
                h["fc"] = np.zeros((1,))
        (d / "phono3py_disp.yaml").write_text("dummy: 1\n", encoding="utf-8")
    return skill


def _run_gen(skill):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(OPT) + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("LIFETIME_REEXEC", None)
    return subprocess.run([sys.executable, str(SKILL / "gen_step1_lifetime.py")],
                          cwd=str(skill), capture_output=True, text=True, env=env, timeout=300)


def test_gen_end_to_end(tmp_path):
    skill = _make_project(tmp_path)
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    out = skill / "step1_lifetime"
    js = json.loads((out / "lifetime_summary.json").read_text(encoding="utf-8"))
    assert js["LIFETIME_DONE"] is True
    assert js["tau_factor"] == pytest.approx(2 * 2 * np.pi)
    assert js["source"]["kind"] == "upstream_hdf5"
    assert js["source"]["mesh"] == [2, 2, 2]
    assert len(js["temperatures_analyzed_K"]) == 2
    assert js["weights_applied"] is True
    assert js["per_temperature"][0]["acoustic_fmax_THz"] is None
    assert js["boundary_gamma_included"] is False
    assert (out / "lifetime_data.csv").is_file()
    assert (out / "lifetime_T.csv").is_file()
    hdr = (out / "lifetime_data.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "gamma_iso_THz" in hdr and "gamma_total_THz" in hdr and "tau_anh_ps" in hdr
    ac = [d for d in js["per_temperature"] if d["T_K"] == 300.0][0]["acoustic"]
    assert ac["tau_median_ps"] == pytest.approx(1.0 / (0.01 * 2 * 2 * np.pi))
    # 逐模 CSV 默认只写 FOCUS_T，避免大网格 × 全温度写出上千万行
    assert js["csv_temperatures_K"] == [300.0]
    rows = (out / "lifetime_data.csv").read_text(encoding="utf-8").splitlines()[1:]
    assert {row.split(",")[0] for row in rows} == {"300.0000"}


def test_gen_weighted_median_end_to_end(tmp_path):
    # 权重 1:1:47，中位数必须按权重算（旧版按不可约点等权）
    skill = _make_project(tmp_path, weights=[1.0, 1.0, 47.0, 1.0])
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    js = json.loads((skill / "step1_lifetime" / "lifetime_summary.json").read_text(encoding="utf-8"))
    ac = [d for d in js["per_temperature"] if d["T_K"] == 300.0][0]["acoustic"]
    assert ac["n_modes_weighted"] == pytest.approx(49.0 * 3)
    assert ac["n_modes"] == 12 - 3            # 4 gp x 3 band，gp0 的 gamma=0 被排除


def test_gen_empty_acoustic_group_does_not_crash(tmp_path):
    # ACOUSTIC_FMAX 极小 -> 声学组为空；旧版在写 lifetime_T.csv 时 TypeError
    skill = _make_project(tmp_path, extra_conf="ACOUSTIC_FMAX = 0.05")
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    js = json.loads((skill / "step1_lifetime" / "lifetime_summary.json").read_text(encoding="utf-8"))
    assert js["per_temperature"][0]["acoustic"]["n_modes"] == 0
    line = (skill / "step1_lifetime" / "lifetime_T.csv").read_text(encoding="utf-8").splitlines()[1]
    assert ",," in line                       # 空分组写成空串，不是 "nan"/崩


def test_gen_single_point_warning_when_no_isotope(tmp_path):
    skill = _make_project(tmp_path)           # 合成 hdf5 没有 gamma_isotope
    r = _run_gen(skill)
    assert r.returncode == 0
    assert "没有 gamma_isotope" in (r.stdout + r.stderr)


def test_gen_errors_clearly_without_source(tmp_path):
    skill = _make_project(tmp_path, with_hdf5=False)
    r = _run_gen(skill)
    assert r.returncode != 0
    assert "找不到声子寿命数据源" in (r.stdout + r.stderr)


def test_gen_fc_requires_opt_in(tmp_path):
    skill = _make_project(tmp_path, with_hdf5=False, with_fc=True)
    r = _run_gen(skill)
    assert r.returncode != 0
    assert "RUN_PHONO3PY" in (r.stdout + r.stderr)


def test_gen_local_main_file_beats_single_point(tmp_path):
    # 本地目录同时有主文件与 --gp 单点文件：必须读主文件（旧版取排序第一个）
    skill = _make_project(tmp_path)
    _write_fake_hdf5(skill / "kappa-m333.hdf5", mesh=(3, 3, 3), n_gp=4)
    _write_fake_hdf5(skill / "kappa-m222-g0.hdf5", mesh=(2, 2, 2), n_gp=1, gamma_is_2d=True)
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    js = json.loads((skill / "step1_lifetime" / "lifetime_summary.json").read_text(encoding="utf-8"))
    assert js["source"]["kind"] == "local_hdf5"
    assert js["source"]["n_grid_points"] == 4
    assert "kappa-m333.hdf5" in js["source"]["file"]


def test_gen_local_single_point_does_not_shadow_upstream(tmp_path):
    # 本地只有 -g0 单点文件时，应优先用上游完整主文件，而不是拿单点文件去读
    skill = _make_project(tmp_path)
    _write_fake_hdf5(skill / "kappa-m222-g0.hdf5", mesh=(2, 2, 2), n_gp=1, gamma_is_2d=True)
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    js = json.loads((skill / "step1_lifetime" / "lifetime_summary.json").read_text(encoding="utf-8"))
    assert js["source"]["n_grid_points"] == 4
    assert js["source"]["file"].endswith("kappa-m222.hdf5")


def test_gen_kl_mlff_upstream_end_to_end(tmp_path):
    """C4 e2e：上游是 kl-mlff-cpu（κ 在 step4_kappa）时 gen 也能找到数据源。"""
    skill = _make_project(tmp_path, upstream_skill="kl-mlff-cpu")
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    js = json.loads((skill / "step1_lifetime" / "lifetime_summary.json").read_text(encoding="utf-8"))
    assert js["source"]["kind"] == "upstream_hdf5"
    assert "kl-mlff-cpu" in js["source"]["file"]
    assert js["source"]["file"].endswith("kappa-m222.hdf5")


# ---------------------------------------------------------------------------
# 真实 Si 数据（本机有才跑）
# ---------------------------------------------------------------------------
_SI_H5 = SI_REF / "kappa-m151515.hdf5"


@pytest.mark.skipif(not (SI_REF / "fc3.hdf5").is_file(), reason="本机无 Si 的 fc2/fc3")
def test_boundary_gamma_real_phono3py_bmfp(tmp_path):
    """C1 真跑：--bmfp 0.1（微米）后，计入 gamma_bd 的 κ 逐模重建误差 < 1e-6。"""
    pytest.importorskip("h5py")
    exe = shutil.which("phono3py-load")
    if exe is None:
        pytest.skip("PATH 里没有 phono3py-load")
    for n in ("fc2.hdf5", "fc3.hdf5", "phono3py_disp.yaml"):
        shutil.copyfile(str(SI_REF / n), str(tmp_path / n))
    r = subprocess.run([exe, "phono3py_disp.yaml", "--mesh", "3", "3", "3", "--br",
                        "--ts=300", "--isotope", "--bmfp", "0.1"],
                       cwd=str(tmp_path), capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, (r.stdout[-1500:] + r.stderr[-1500:])
    h5 = tmp_path / "kappa-m333.hdf5"
    assert h5.is_file()
    d = LC.read_kappa_hdf5(str(h5))
    L_um = float(np.ravel(d["boundary_mfp"])[0])
    assert L_um == pytest.approx(0.1)                # hdf5 里存的就是微米
    gt = LC.effective_gamma(d["gamma"], d.get("gamma_isotope"))
    g_bd = LC.boundary_gamma(d["group_velocity"], L_um)
    use, info = LC.decide_boundary_gamma(d, gt, g_bd, t_indices=[0])
    assert use is True
    assert info["err_with_boundary"] < 1e-6
    assert info["err_no_boundary"] > 1e-3


@pytest.mark.skipif(not _SI_H5.is_file(), reason="本机无 Si 参考数据")
def test_real_si_kappa_reconstruction_and_linewidth():
    d = LC.read_kappa_hdf5(str(_SI_H5))
    gi = d.get("gamma_isotope")
    gt = LC.effective_gamma(d["gamma"], gi)
    chk = LC.kappa_reconstruction_error(d, gt, t_indices=list(range(len(d["temperatures"]))))
    assert chk is not None and chk["max_rel_err"] < 1e-4      # 逐模重建与上游 kappa 一致
    # 默认 boundary_mfp=1e6 Å 是占位，计入反而变差 -> 不应计入
    g_bd = LC.boundary_gamma(d.get("group_velocity"), float(np.ravel(d["boundary_mfp"])[0]))
    use, _ = LC.decide_boundary_gamma(d, gt, g_bd, t_indices=[2])
    assert use is False
    # Γ 点最高支（光学）线宽应落在 Si 室温 Raman 自然线宽量级（~3 cm^-1）
    q = d["qpoints"]
    ig = np.where(np.all(np.abs(q) < 1e-8, axis=1))[0]
    assert ig.size == 1
    b = int(np.argmax(d["frequency"][ig[0]]))
    fwhm_cm = 2.0 * float(gt[2, ig[0], b]) * 33.35641
    assert 1.5 < fwhm_cm < 6.0, fwhm_cm
    # 权重生效：展开到全网格的中位数应与加权中位数一致
    s = LC.summarize_temperature(d["frequency"], d["gamma"], d["temperatures"], 2,
                                 gamma_isotope=gi, weights=d["weights"])
    gt2 = LC.effective_gamma(d["gamma"][2], gi)
    tau = LC.gamma_to_tau_ps(gt2)
    f2 = d["frequency"]
    w2 = np.broadcast_to(d["weights"][:, None], f2.shape)
    sel = (f2 > 0.1) & (gt2 > 0) & np.isfinite(tau)
    rep = np.repeat(tau[sel], w2[sel].astype(int))
    assert s["all"]["tau_median_ps"] == pytest.approx(float(np.median(rep)))


# ---------------------------------------------------------------------------
# 2D κ 归一化（κ_raw 含真空 -> κ_2d = κ_raw × h⊥/d；因子优先读上游 summary）
# ---------------------------------------------------------------------------
_MOS2_YAML = """primitive_cell:
  lattice:
  - [ 3.1600000000000,  0.0000000000000,  0.0000000000000 ]
  - [ -1.5800000000000, 2.7366403500000,  0.0000000000000 ]
  - [ 0.0000000000000,  0.0000000000000, 25.0000000000000 ]
  points:
  - symbol: Mo
    coordinates: [ 0.3333333333333, 0.6666666666667, 0.5000000000000 ]
  - symbol: S
    coordinates: [ 0.6666666666667, 0.3333333333333, 0.5624000000000 ]
  - symbol: S
    coordinates: [ 0.6666666666667, 0.3333333333333, 0.4376000000000 ]
"""


def _write_2d_upstream(d, with_summary=True, with_yaml=True, factor=3.7215, d_A=6.7177):
    d.mkdir(parents=True, exist_ok=True)
    _write_fake_hdf5(d / "kappa-m881.hdf5", mesh=(8, 8, 1))
    if with_summary:
        (d / "kappa_summary.json").write_text(json.dumps(
            {"KAPPA_DONE": True, "dim": "2d", "kappa_2d_norm_factor": factor,
             "thickness_d_ang": d_A, "h_perp_A": 25.0,
             "thickness_convention": "span 3.12 + vdW(S 1.80 + S 1.80)"}), encoding="utf-8")
    if with_yaml:
        (d / "phono3py.yaml").write_text(_MOS2_YAML, encoding="utf-8")


def test_vacuum_axis_from_mesh():
    assert LC.vacuum_axis_from_mesh([138, 138, 1]) == 2
    assert LC.vacuum_axis_from_mesh([1, 40, 40]) == 0
    assert LC.vacuum_axis_from_mesh([15, 15, 15]) is None
    assert LC.vacuum_axis_from_mesh([1, 1, 40]) is None          # 1D，不当 2D


def test_two_d_uses_upstream_factor(tmp_path):
    d = tmp_path / "step6_kappa"
    _write_2d_upstream(d)
    n = LC.two_d_normalization(d / "kappa-m881.hdf5", [8, 8, 1],
                               yaml_path=d / "phono3py.yaml")
    assert n["dim"] == "2d" and n["source"] == "upstream_kappa_summary"
    assert n["factor"] == pytest.approx(3.7215)
    assert n["thickness_d_A"] == pytest.approx(6.7177)
    if n["check"] is not None:                                   # 有 pyyaml 时两路对照
        assert n["check"]["rel_diff"] < 0.05


def test_two_d_falls_back_to_slab_geometry(tmp_path):
    pytest.importorskip("yaml")
    d = tmp_path / "step6_kappa"
    _write_2d_upstream(d, with_summary=False)
    n = LC.two_d_normalization(d / "kappa-m881.hdf5", [8, 8, 1],
                               yaml_path=d / "phono3py.yaml")
    assert n["dim"] == "2d" and n["source"].startswith("slab_geometry")
    # 与 kl 同一真源：span(3.12) + 2×r_vdW(S)
    assert n["h_perp_A"] == pytest.approx(25.0, abs=1e-3)
    assert n["factor"] == pytest.approx(25.0 / n["thickness_d_A"], rel=1e-4)
    assert 5.5 < n["thickness_d_A"] < 7.5


def test_two_d_fixed_thickness_mode(tmp_path):
    pytest.importorskip("yaml")
    d = tmp_path / "step6_kappa"
    _write_2d_upstream(d, with_summary=False)
    n = LC.two_d_normalization(d / "kappa-m881.hdf5", [8, 8, 1], thickness_mode="6.15",
                               yaml_path=d / "phono3py.yaml")
    assert n["factor"] == pytest.approx(25.0 / 6.15, rel=1e-4)


def test_two_d_respects_explicit_3d_and_upstream_3d(tmp_path):
    d = tmp_path / "step6_kappa"
    _write_2d_upstream(d)
    assert LC.two_d_normalization(d / "kappa-m881.hdf5", [8, 8, 1], dim="3d")["dim"] == "3d"
    (d / "kappa_summary.json").write_text(json.dumps({"dim": "3d"}), encoding="utf-8")
    assert LC.two_d_normalization(d / "kappa-m881.hdf5", [8, 8, 1])["dim"] == "3d"


def test_two_d_without_any_factor_warns(tmp_path):
    d = tmp_path / "step6_kappa"
    _write_2d_upstream(d, with_summary=False, with_yaml=False)
    msgs = []
    n = LC.two_d_normalization(d / "kappa-m881.hdf5", [8, 8, 1], log=msgs.append)
    assert n["dim"] == "2d" and n["factor"] is None
    assert any("不能直接对文献" in m for m in msgs)


def test_kappa_block_2d_inplane():
    data = {"temperatures": np.array([100.0, 300.0]),
            "kappa": np.array([[30.0, 30.0, 0.0, 0, 0, 0], [27.93, 27.93, 1e-9, 0, 0, 0]])}
    norm = {"dim": "2d", "factor": 3.7215, "vac_axis": 2, "thickness_d_A": 6.7177,
            "source": "upstream_kappa_summary"}
    k = LC.kappa_block(data, 300, norm)
    assert k["kappa_raw_inplane_avg"] == pytest.approx(27.93)
    assert k["kappa_2d_inplane_avg"] == pytest.approx(27.93 * 3.7215)   # = 103.94
    assert set(k["kappa_2d_inplane"]) == {"xx", "yy"}                   # zz 不报
    assert k["kappa_2d_inplane_avg_vs_T"]["100.0000"] == pytest.approx(30.0 * 3.7215)
    k3 = LC.kappa_block(data, 300, {"dim": "3d"})
    assert "kappa_2d_inplane_avg" not in k3 and k3["dim"] == "3d"


def test_gen_2d_end_to_end_reports_normalized_kappa(tmp_path):
    up = tmp_path / "Si" / "kl-dft-cpu" / "result" / "step6_kappa"
    skill = _make_project(tmp_path, with_hdf5=False, upstream_dir=up)
    _write_2d_upstream(up)
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    js = json.loads((skill / "step1_lifetime" / "lifetime_summary.json").read_text(encoding="utf-8"))
    assert js["dim"] == "2d"
    kb = js["kappa_from_hdf5"]
    assert kb["kappa_2d_norm_factor"] == pytest.approx(3.7215)
    assert kb["kappa_2d_inplane_avg"] == pytest.approx(10.0 * 3.7215)
    assert "κ_2d" in r.stdout


def test_gen_3d_has_no_2d_block(tmp_path):
    skill = _make_project(tmp_path)
    r = _run_gen(skill)
    assert r.returncode == 0, r.stdout + r.stderr
    js = json.loads((skill / "step1_lifetime" / "lifetime_summary.json").read_text(encoding="utf-8"))
    assert js["dim"] == "3d"
    assert "kappa_2d_inplane_avg" not in js["kappa_from_hdf5"]
