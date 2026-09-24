# -*- coding: utf-8 -*-
"""phonon-lifetime 技能回归：寿命约定 tau=1/(2*2*pi*gamma) + gen 端到端。

约定必须与 phono3py 官方后处理一致（否则寿命会差 2*pi 倍）：
    phono3py/scripts/phono3py_kdeplot.py: tau = 1/g/(2*2*pi)
    phono3py/other/kaccum.py:            mfp = v/(2*2*pi*g)
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
sys.path.insert(0, str(SKILL))
sys.path.insert(0, str(ROOT / "skill" / "_common" / "opt"))

LC = pytest.importorskip("lifetime_common")


# ---------------------------------------------------------------------------
# 约定 / 纯函数
# ---------------------------------------------------------------------------
def test_tau_convention_matches_phono3py():
    gamma = np.array([1.0, 0.05, 0.0, -0.2])
    tau = LC.gamma_to_tau_ps(gamma)
    assert tau[0] == pytest.approx(1.0 / (2 * 2 * np.pi))
    assert tau[1] == pytest.approx(1.0 / (0.05 * 2 * 2 * np.pi))
    assert np.isinf(tau[2]) and np.isinf(tau[3])
    # 与 phono3py kdeplot 的写法逐位一致
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


def test_pick_temperature_indices():
    Ts = np.array([100.0, 200.0, 300.0, 400.0])
    assert set(LC.pick_temperature_indices(Ts).values()) == {0, 1, 2, 3}
    got = LC.pick_temperature_indices(Ts, [290.0, 300.0])
    assert sorted(got.values()) == [2]           # 都映射到 300K 去重
    assert 300.0 in got


def test_acoustic_mask_lowest_three_bands():
    f = np.array([[1.0, 2.0, 3.0, 10.0, 11.0, 12.0]])
    g = np.full_like(f, 0.01)
    mask, fmax = LC.acoustic_mask(f, g, acoustic_fmax=5.0)
    assert list(mask[0]) == [True, True, True, False, False, False]
    assert fmax == 5.0


def test_summarize_temperature():
    f = np.array([[1.0, 2.0, 3.0, 10.0]])
    g = np.full((1, 1, 4), 0.01)
    s = LC.summarize_temperature(f, g, np.array([300.0]), 0, acoustic_fmax=5.0)
    assert s["T_K"] == 300.0
    assert s["acoustic"]["n_modes"] == 3
    assert s["optical"]["n_modes"] == 1
    assert s["acoustic"]["tau_median_ps"] == pytest.approx(1.0 / (0.01 * 2 * 2 * np.pi))


def test_effective_gamma_adds_isotope_only_when_positive():
    g = np.array([[1.0, 2.0]])
    assert np.allclose(LC.effective_gamma(g, None), g)
    assert np.allclose(LC.effective_gamma(g, np.zeros((1, 2))), g)
    assert np.allclose(LC.effective_gamma(g, np.full((1, 2), 0.5)), g + 0.5)


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


def test_select_direction_qpoints():
    q = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.25], [0.0, 0.0, 0.5],
                  [0.25, 0.0, 0.0], [0.5, 0.0, 0.0]])
    picked = LC.select_direction_qpoints(q, [0, 0, 1], n_points=2, tol=0.05)
    gp = [i for _, i in picked]
    assert 0 in gp and 2 in gp          # Gamma 与 0.5 c* 都在线上
    assert 3 not in gp and 4 not in gp  # 垂直方向的点被排除


# ---------------------------------------------------------------------------
# gen 端到端（合成 kappa hdf5，不碰集群）
# ---------------------------------------------------------------------------
def _write_fake_kappa(path):
    h5py = pytest.importorskip("h5py")
    T = np.array([100.0, 300.0])
    ngp, nb = 4, 3
    freq = np.array([[0.0, 1.0, 5.0],
                     [0.5, 2.0, 6.0],
                     [1.0, 3.0, 7.0],
                     [1.5, 4.0, 8.0]], dtype=float)
    gamma = np.full((2, ngp, nb), 0.01, dtype=float)
    gamma[:, 0, :] = 0.0
    with h5py.File(path, "w") as h:
        h["temperature"] = T
        h["frequency"] = freq
        h["gamma"] = gamma
        h["qpoint"] = np.array([[0, 0, 0], [0, 0, .25], [0, 0, .5], [0, 0, .75]], dtype=float)
        h["weight"] = np.ones(ngp)
        h["mesh"] = np.array([2, 2, 2])
        h["group_velocity"] = np.ones((ngp, nb, 3))
        h["heat_capacity"] = np.ones((2, ngp, nb))
        h["kappa"] = np.tile(np.array([10.0, 10.0, 10.0, 0, 0, 0]), (2, 1))
        h["kappa_unit_conversion"] = np.array(1.0)


def _make_project(tmp_path, with_hdf5=True, with_fc=False):
    skill = tmp_path / "Si" / "phonon-lifetime"
    skill.mkdir(parents=True)
    shutil.copy(str(SKILL / "templates" / "step.conf"), str(skill / "step.conf"))
    d = tmp_path / "Si" / "kl-dft-cpu" / "result" / "step6_kappa"
    d.mkdir(parents=True)
    if with_hdf5:
        _write_fake_kappa(str(d / "kappa-m222.hdf5"))
    if with_fc:
        h5py = pytest.importorskip("h5py")
        for name in ("fc2.hdf5", "fc3.hdf5"):
            with h5py.File(str(d / name), "w") as h:
                h["fc"] = np.zeros((1,))
        (d / "phono3py_disp.yaml").write_text("dummy: 1\n", encoding="utf-8")
    return skill


def _run_gen(skill):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "skill" / "_common" / "opt") + os.pathsep + \
        env.get("PYTHONPATH", "")
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
    assert (out / "lifetime_data.csv").is_file()
    assert (out / "lifetime_T.csv").is_file()
    hdr = (out / "lifetime_data.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "gamma_iso_THz" in hdr and "gamma_total_THz" in hdr and "tau_anh_ps" in hdr
    # 300K 声学中位寿命应等于 1/(4*pi*0.01)
    ac = [d for d in js["per_temperature"] if d["T_K"] == 300.0][0]["acoustic"]
    assert ac["tau_median_ps"] == pytest.approx(1.0 / (0.01 * 2 * 2 * np.pi))


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
