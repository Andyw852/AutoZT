"""[0072] 谷内 ADP 形变势张量用 DFPT 声学支斜率（eph-qe-cpu IVS_ADP -> ke-dft-cpu S8 的 ADP_DFPT）。

不需要 QE：用合成的声学模本征矢（含内应变混入、简并 TA 任意转动、每个模任意相位）与按模型 G_ν = Σ p_ν,i w_i 造出来的
ahc_gkk，核对 adp_vector / adp_real / adp_fit 把对称张量 D 原样取回（只差整体符号）；星展开（Si 导带 Δ 谷）；
S10 的作业规划（反演中心 / 二维 / IVS_ADP 开关）；S8 的 ADP_DFPT 读取。
"""
import importlib.util
import io
import json
import math
import shutil
import sys
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
IVSDIR = ROOT / "skill" / "_common" / "ivs"
EPH = ROOT / "skill" / "eph-qe-cpu"
KE = ROOT / "skill" / "ke-dft-cpu"
for _p in (str(IVSDIR), str(ROOT / "skill" / "_common" / "opt"), str(EPH), str(KE), str(KE / "step8.4_amset2d")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ivs_dfpt as iv  # noqa: E402

A_SI = 5.479
SI_LAT = np.array([[0, A_SI / 2, A_SI / 2], [A_SI / 2, 0, A_SI / 2], [A_SI / 2, A_SI / 2, 0]])
SI_FRAC = np.array([[0, 0, 0], [0.25, 0.25, 0.25]])
X_VALLEY = [0.0, 0.425, 0.425]          # 0.85 ΓX，笛卡尔 +x


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rand_sym(rng, scale=5.0):
    a = rng.normal(size=(3, 3)) * scale
    return 0.5 * (a + a.T)


def _synthetic_run(rng, D, n_hat, dq, masses=(28.0855, 28.0855)):
    """一个方向的合成 (gkk, e)：2 原子，声学模 = 平移 + 内应变混入（∝ q）+ 声学子空间任意幺正转动 + 每个模任意相位。"""
    m = np.asarray(masses, float)
    M = m.sum()
    nat = len(m)
    T = np.zeros((3 * nat, 3))
    for a in range(3):
        T[a::3, a] = np.sqrt(m / M)
    opt = np.zeros((3 * nat, 3))                         # 相对位移（质心不动），质量加权
    for a in range(3):
        opt[a::3, a] = np.array([math.sqrt(m[1] / M), -math.sqrt(m[0] / M)])
    mix = rng.normal(size=(3, 3)) * 0.3 * dq             # 内应变：∝ q 的光学支混入
    ac = T + opt @ mix
    ac, _ = np.linalg.qr(ac)
    U = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]
    ac = ac @ U * np.exp(1j * rng.uniform(0, 2 * np.pi, 3))[None, :]
    full = np.concatenate([ac, opt.astype(complex)], axis=1)
    e, _ = np.linalg.qr(full)                             # 正交化，前 3 列张成声学子空间
    e[:, :3] = ac / np.linalg.norm(ac, axis=0)
    # 模型：w = i|q|·D·n̂ × 公共相位（电子规范）
    w = 1j * dq * (D @ n_hat) * np.exp(1j * rng.uniform(0, 2 * np.pi))
    P = iv.adp_polarizations(e, m * iv.AMU_RY, [0, 1, 2])
    G = P @ w                                             # 每个声学模的耦合（eV/Å）
    u = e / np.repeat(np.sqrt(m * iv.AMU_RY), 3)[:, None]
    target = G / (math.sqrt(float((m * iv.AMU_RY).sum())) * iv.GKK_TO_EV_A)
    g = np.linalg.pinv(u[:, :3].T) @ target               # 最小范数解：u_acᵀ g = target
    gkk = g.reshape(3 * nat, 1, 1)
    return gkk, e, m * iv.AMU_RY


def test_q_fracs_cartesian_length_and_direction():
    for n, qf in iv.adp_q_fracs(SI_LAT, 0.03):
        qc = np.asarray(qf) @ iv.recip_lattice(SI_LAT)
        assert np.linalg.norm(qc) == pytest.approx(0.03, rel=1e-10)
        assert np.allclose(qc / np.linalg.norm(qc), n)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_vector_and_fit_recover_tensor(seed):
    rng = np.random.default_rng(seed)
    D = _rand_sym(rng)
    dq = 0.03
    vs, dirs, ims = [], [], []
    for n, _qf in iv.adp_q_fracs(SI_LAT, dq):
        gkk, e, m = _synthetic_run(rng, D, n, dq)
        ac, _tw = iv.acoustic_modes(e, m)
        assert ac == [0, 1, 2]
        w = iv.adp_vector(gkk, e, m, 0, 0, ac)
        v, im = iv.adp_real(w)
        assert im < 1e-8
        vs.append(v)
        dirs.append(n)
    Dfit, resid, signs = iv.adp_fit(dirs, vs, dq)
    assert resid < 1e-8
    assert min(np.max(np.abs(Dfit - D)), np.max(np.abs(Dfit + D))) < 1e-6


def test_fit_flags_inconsistent_data():
    rng = np.random.default_rng(5)
    D = _rand_sym(rng)
    dirs = [n for n, _ in iv.adp_q_fracs(SI_LAT, 0.03)]
    vs = [0.03 * (D @ n) for n in dirs]
    vs[3] = vs[3] + rng.normal(size=3) * 2.0                # 一个方向的数据不服从对称张量
    _D, resid, _s = iv.adp_fit(dirs, vs, 0.03)
    assert resid > iv.ADP_FIT_TOL


def test_member_tensors_si_delta_star():
    rots = iv.point_group(SI_LAT, SI_FRAC, ["Si", "Si"])
    members, norb, _ = iv.expand_star([X_VALLEY], iv.recip_ops(rots), full3d=True)
    assert norb == 1 and len(members) == 6
    xid, xiu = 1.01, 8.84                                    # Li 2021：Ξd、Ξu
    kc = np.asarray(X_VALLEY) @ iv.recip_lattice(SI_LAT)
    a = kc / np.linalg.norm(kc)
    D = xid * np.eye(3) + xiu * np.outer(a, a)
    sm = iv.adp_summary(D, X_VALLEY, SI_LAT)
    assert sm["axis_eV"] == pytest.approx(xid + xiu) and sm["hydrostatic_eV"] == pytest.approx(xid + xiu / 3)
    tens = iv.adp_member_tensors(D, X_VALLEY, members, SI_LAT, rots)
    for km, Dm in zip(members, tens):
        c = np.asarray(km) @ iv.recip_lattice(SI_LAT)
        c = c / np.linalg.norm(c)
        assert float(c @ Dm @ c) == pytest.approx(xid + xiu, abs=1e-9)        # 每个成员沿自己的谷轴
        assert np.allclose(np.sort(np.linalg.eigvalsh(Dm)), [xid, xid, xid + xiu])
    with pytest.raises(ValueError):
        iv.adp_member_tensors(D, X_VALLEY, [[0.1, 0.2, 0.3]], SI_LAT, rots)


def _params_with_adp():
    rots = iv.point_group(SI_LAT, SI_FRAC, ["Si", "Si"])
    members, _n, _a = iv.expand_star([X_VALLEY], iv.recip_ops(rots), full3d=True)
    kc = np.asarray(X_VALLEY) @ iv.recip_lattice(SI_LAT)
    a = kc / np.linalg.norm(kc)
    D = 1.01 * np.eye(3) + 8.84 * np.outer(a, a)
    tens = iv.adp_member_tensors(D, X_VALLEY, members, SI_LAT, rots)
    row = {"label": "X", "k0": X_VALLEY, "status": "ok", "D_ref_eV": D, "signs": [1] * 6,
           "checks": {"fit_rel_resid": 0.001, "im_ratio_max": 0.001, "acoustic_hw_meV": [[1.0, 1.0, 2.0]]},
           "warnings": [], "summary": iv.adp_summary(D, X_VALLEY, SI_LAT),
           "members": [{"k_frac": list(map(float, k)), "D_eV": M} for k, M in zip(members, tens)]}
    adp = {"status": "ok", "dq_A": 0.03, "lattice_A": SI_LAT, "carriers": {
        "electron": {"status": "ok", "valleys": [row], "notes": []},
        "hole": {"status": "degenerate", "valleys": [{"label": "Γ", "k0": [0, 0, 0], "status": "degenerate",
                                                      "why": "带边简并"}], "notes": []}}}
    return iv.jsonable({"schema": iv.SCHEMA, "dim": "3d", "carriers": {}, "raw": [], "adp": adp})


def test_consumer_api_and_digest():
    p = _params_with_adp()
    assert iv.adp_params(p, "electron")["valleys"][0]["label"] == "X"
    assert iv.adp_params(p, "hole") is None
    d1 = iv.adp_digest(p, "electron")
    p2 = json.loads(json.dumps(p))
    p2["adp"]["carriers"]["electron"]["valleys"][0]["members"][0]["D_eV"][0][0] += 0.1
    assert d1 and iv.adp_digest(p2, "electron") != d1 and iv.adp_digest(p, "hole") is None


# ---------------------------------------------------------------------------
# S10 规划
# ---------------------------------------------------------------------------
@pytest.fixture()
def eph(tmp_path, monkeypatch):
    for f in ("ivs_qe_common.py", "ivs_qe_job.py", "eph_common.py", "submit_ivs_qe.tpl"):
        shutil.copy(EPH / f, tmp_path)
    shutil.copy(IVSDIR / "ivs_dfpt.py", tmp_path)
    shutil.copy(ROOT / "skill" / "_common" / "opt" / "stepconf.py", tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.syspath_prepend(str(tmp_path))
    for m in ("ivs_qe_common", "stepconf", "eph_common"):
        sys.modules.pop(m, None)
    import ivs_qe_common as c
    return tmp_path, c


def _carriers():
    rots = iv.point_group(SI_LAT, SI_FRAC, ["Si", "Si"])
    members, _n, _a = iv.expand_star([X_VALLEY], iv.recip_ops(rots), full3d=True)
    return rots, {"electron": [{"label": "X", "k": members}], "hole": [{"label": "G", "k": [[0.0, 0.0, 0.0]]}]}


def test_plan_runs_centrosymmetric(eph, capsys):
    _d, c = eph
    rots, car = _carriers()
    runs = [{"id": 0, "q": [0, 0, 0], "k": [], "name": "q-odp"}]
    rec = c.add_adp_runs({"IVS_ADP": "auto", "IVS_ADP_DQ": 0.03}, car, runs, {"lattice": SI_LAT}, rots, "3d")
    assert rec["status"] == "planned" and len(rec["dirs"]) == 6 and len(runs) == 7
    assert [r["name"] for r in runs[1:]] == ["q-adp%d" % i for i in range(1, 7)]
    assert all(len(r["k"]) == 2 for r in runs[1:])           # 电子 X 中心 + 空穴 Γ
    assert [e["carrier"] for e in rec["entries"]] == ["electron", "hole"]


def test_plan_skips_2d_noncentro_and_off(eph, capsys):
    _d, c = eph
    rots, car = _carriers()
    runs = []
    assert c.add_adp_runs({"IVS_ADP": "off"}, car, runs, {"lattice": SI_LAT}, rots, "3d")["status"] == "off"
    assert c.add_adp_runs({"IVS_ADP": "auto"}, car, runs, {"lattice": SI_LAT}, rots, "2d")["status"] == "skipped"
    zb = iv.point_group(SI_LAT, SI_FRAC, ["Ga", "As"])        # 闪锌矿：无反演中心
    assert not any(np.array_equal(np.asarray(W), -np.eye(3, dtype=int)) for W in zb)
    assert c.add_adp_runs({"IVS_ADP": "auto"}, car, runs, {"lattice": SI_LAT}, zb, "3d")["status"] == "skipped"
    assert runs == []
    rec = c.add_adp_runs({"IVS_ADP": "on", "IVS_ADP_DQ": 0.03}, car, runs, {"lattice": SI_LAT}, zb, "3d")
    assert rec["status"] == "planned" and len(runs) == 6
    with pytest.raises(SystemExit):
        c.add_adp_runs({"IVS_ADP": "on"}, car, [], {"lattice": SI_LAT}, rots, "2d")


# ---------------------------------------------------------------------------
# ke-dft-cpu S8
# ---------------------------------------------------------------------------
@pytest.fixture()
def g10(tmp_path):
    g = _load("g10_adpdfpt", KE / "step8_amset" / "gen_step10_amset.py")
    cwd = tmp_path / "ke"
    (cwd / "ivs_dfpt").mkdir(parents=True)
    return g, cwd


def _resolve(g, cwd):
    buf = io.StringIO()
    with redirect_stdout(buf):
        rec = g.resolve_adp_dfpt_3d(cwd)
    return rec, buf.getvalue()


def test_s8_reads_adp_tensor(g10):
    g, cwd = g10
    (cwd / "ivs_dfpt" / "ivs_dfpt_params.json").write_text(json.dumps(_params_with_adp()))
    g.ADP_DFPT, g.ADP_DFPT_EWIN = "auto", 0.3
    rec, log = _resolve(g, cwd)
    assert rec["use"] == {"electron": "used", "hole": "degenerate"}
    e = rec["carriers"]["electron"]
    assert len(e["valleys"][0]["members"]) == 6 and e["digest"]
    assert e["summary"][0]["axis_eV"] == pytest.approx(9.85)
    assert np.allclose(rec["lattice_A"], SI_LAT) and "谷轴 9.85 eV" in log
    assert "ADP_DFPT" in g.SPEC


def test_s8_adp_modes(g10):
    g, cwd = g10
    g.ADP_DFPT = "auto"
    rec, _ = _resolve(g, cwd)                                  # 没有文件
    assert rec["carriers"] == {} and rec["use"]["electron"] == "unavailable"
    g.ADP_DFPT = "on"
    with pytest.raises(SystemExit):
        _resolve(g, cwd)
    p = _params_with_adp()
    p["adp"] = {"status": "skipped", "why": "二维"}
    (cwd / "ivs_dfpt" / "ivs_dfpt_params.json").write_text(json.dumps(p))
    with pytest.raises(SystemExit):
        _resolve(g, cwd)
    g.ADP_DFPT = "off"
    rec, _ = _resolve(g, cwd)
    assert rec["carriers"] == {} and rec["use"]["electron"] == "off"
