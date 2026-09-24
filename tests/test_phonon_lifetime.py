# -*- coding: utf-8 -*-
"""phonon-lifetime 技能回归。

约定：tau[ps] = 1/(2*2*pi*gamma_total[THz])，与 phono3py 官方后处理一致
（scripts/phono3py_kdeplot.py 的 tau = 1/g/(2*2*pi)、other/kaccum.py 的 mfp）。
覆盖三类：换算约定、加权统计/选文件/对称方向等纯函数、gen 端到端。
"""
import json
import os
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
SI_REF = Path("/mnt/d/tf_data/work_AutoZT/Si/kl-dft-cpu/result/step6_kappa")


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


def test_boundary_gamma_formula():
    gv = np.zeros((1, 1, 3))
    gv[0, 0] = [0.0, 0.0, 4.0]
    gbd = LC.boundary_gamma(gv, 1.0)
    assert gbd[0, 0] == pytest.approx(4.0 / (4 * np.pi))
    assert LC.boundary_gamma(gv, 0) is None
    assert LC.boundary_gamma(None, 100.0) is None


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
                  hdf5_name="kappa-m222.hdf5", weights=None):
    skill = tmp_path / "Si" / "phonon-lifetime"
    skill.mkdir(parents=True)
    conf = (SKILL / "templates" / "step.conf").read_text(encoding="utf-8")
    if extra_conf:
        conf = conf + "\n" + extra_conf
    (skill / "step.conf").write_text(conf, encoding="utf-8")
    d = tmp_path / "Si" / "kl-dft-cpu" / "result" / "step6_kappa"
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


# ---------------------------------------------------------------------------
# 真实 Si 数据（本机有才跑）
# ---------------------------------------------------------------------------
_SI_H5 = SI_REF / "kappa-m151515.hdf5"


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
