"""[0035] 冻声子 ODP 参数（odp_frozen.py、S7.2/S8.21 接线、S8.4 的 ODP_FROZEN）与 τ(E)（tau_energy.py 纯函数、
S8.4 作业链、S8.3 τ 对比）。纯 numpy，不需要 amset / pymatgen / VASP。
真 AMSET 端到端在 skill/ke-dft-cpu/step8.4_amset2d/test_tau_energy_e2e.py（需 amset 环境）。"""
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
S84 = KE / "step8.4_amset2d"
S83 = KE / "step8.3_output"
S72 = KE / "step7c_odp_frozen"


def _load(name, path):
    for d in (ROOT, KE, S84, ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def of():
    return _load("odp_frozen_0035", KE / "odp_frozen.py")


@pytest.fixture(scope="module")
def te():
    return _load("tau_energy_0035", S84 / "tau_energy.py")


@pytest.fixture(scope="module")
def kc():
    return _load("ke_common_0035", KE / "ke_common.py")


# ================================================================ 位移的归一化
def _orthonormal(n, seed):
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(n, n)))
    return q.T


def test_displacement_is_normal_coordinate_u_sqrtM(of):
    """Δ_κ = u√(M/m_κ)e_κ <=> 质量加权简正坐标 Q = Σ√m_κ e_κ·Δ_κ = u√M（D = dE/du 与 |g|² = ħD²/(2NMω) 一致的前提）。"""
    m = np.array([95.95, 32.06, 32.06, 12.0])
    E = _orthonormal(12, 1).reshape(12, 4, 3)
    for e in E[:5]:
        d, u = of.displacement(e, m, 0.01)
        Q = float(np.sum(np.sqrt(m)[:, None] * e * d))
        assert Q == pytest.approx(u * math.sqrt(m.sum()), rel=1e-12)
        assert np.max(np.linalg.norm(d, axis=1)) == pytest.approx(0.01, rel=1e-12)
    # 纯平移模（e_κ ∝ √m_κ）：所有原子平移同一个 u
    t = np.zeros((4, 3))
    t[:, 2] = np.sqrt(m / m.sum())
    d, u = of.displacement(t, m, 0.01)
    assert np.allclose(d[:, 2], u) and np.allclose(d[:, :2], 0)


# ================================================================ 紧束缚玩具：冻声子 D = √M·⟨∂H/∂Q⟩（Hellmann–Feynman）
A = 3.2
MASS = np.array([95.95, 32.06, 32.06])
POS0 = np.array([[0.0, 0.0, 0.0], [1.6, 1.6, 1.5], [1.6, 1.6, -1.5]])
IMG = [np.array([i * A, j * A, 0.0]) for i in (-1, 0, 1) for j in (-1, 0, 1)]
ONSITE = np.array([0.0, -1.0, -1.0])


def _pairs(pos, cut=3.6):
    for i in range(3):
        for j in range(3):
            for R in IMG:
                if i == j and not R.any():
                    continue
                d = pos[j] + R - pos[i]
                if np.linalg.norm(d) < cut:
                    yield i, j, R, d


def _phonons():
    """中心弹簧力常数 -> Γ 点动力学矩阵 -> (频率 THz, 质量加权本征矢 (9,3,3))。"""
    Phi = np.zeros((9, 9))
    for i, j, R, d in _pairs(POS0):
        r = d / np.linalg.norm(d)
        k = 5.0 if {i, j} == {0, 1} or {i, j} == {0, 2} else 2.0
        blk = k * np.outer(r, r)
        Phi[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= blk
        Phi[3 * i:3 * i + 3, 3 * i:3 * i + 3] += blk
    msq = np.repeat(np.sqrt(MASS), 3)
    Dm = Phi / np.outer(msq, msq)
    w2, v = np.linalg.eigh((Dm + Dm.T) / 2)
    f = np.sqrt(np.abs(w2)) * np.sign(w2) * 15.633302
    return f, v.T.reshape(9, 3, 3)


def _H(kf, pos):
    H = np.diag(ONSITE).astype(complex)
    for i, j, R, d in _pairs(pos):
        t = -1.2 * math.exp(-(np.linalg.norm(d) - 2.7) / 0.6)
        H[i, j] += t * np.exp(2j * np.pi * (kf[0] * R[0] + kf[1] * R[1]) / A)
    return (H + H.conj().T) / 2


def _mode(E, pattern):
    p = pattern / np.linalg.norm(pattern)
    return int(np.argmax([abs(np.sum(e * p)) for e in E]))


def test_tight_binding_frozen_phonon_equals_hellmann_feynman(of):
    f, E = _phonons()
    born = np.array([-2 * np.eye(3), np.eye(3), np.eye(3)])
    modes, groups, notes = of.mode_table(f, E, MASS, born)
    assert sum(m["acoustic"] for m in modes) == 3 and not notes
    a1 = np.zeros((3, 3))
    a1[1, 2], a1[2, 2] = 1, -1                              # X1 +z、X2 −z（同极模）
    e2 = np.zeros((3, 3))
    e2[1, 0], e2[2, 0] = 1, -1                              # X 面内反相（σh 下为奇）
    ia, ie = _mode(E, a1), _mode(E, e2)
    order = list(np.argsort(f, kind="stable"))
    by_src = {m["_src"]: m for m in modes}
    assert not by_src[ia]["polar"] and not by_src[ie]["polar"]          # Raman 模：模有效电荷 0
    assert any(m["polar"] for m in modes if not m["acoustic"])           # M 对 X 的模是极性的
    used = {i for g in groups if g["use"] for i in g["modes"]}
    assert by_src[ia]["index"] in used and by_src[ie]["index"] in used
    assert all(not m["polar"] for m in modes if m["index"] in used)
    assert order.index(ia) + 1 == by_src[ia]["index"]
    msum = MASS.sum()
    for kf in ([0.0, 0.0], [0.23, 0.11]):
        w0, V0 = np.linalg.eigh(_H(kf, POS0))
        for src in (ia, ie):
            d, u = of.displacement(E[src], MASS, 0.01)
            ep = np.linalg.eigvalsh(_H(kf, POS0 + d))
            em = np.linalg.eigvalsh(_H(kf, POS0 - d))
            dq = 1e-6
            pat = E[src] / np.sqrt(MASS)[:, None]
            dH = (_H(kf, POS0 + dq * pat) - _H(kf, POS0 - dq * pat)) / (2 * dq)
            for b in range(3):
                r = of.edge_coupling(w0, ep, em, b, u)
                Dref = math.sqrt(msum) * float(np.real(V0[:, b].conj() @ dH @ V0[:, b]))
                assert r["d"] == 1 and r["method"] == "central"
                assert math.sqrt(r["D2"]) == pytest.approx(abs(Dref), rel=1e-3, abs=1e-6)
                if abs(Dref) > 0.1:                         # 误差是 O(u²)：位移缩到 1/5，误差应缩 ~25 倍
                    d2_, u2 = of.displacement(E[src], MASS, 0.002)
                    r2 = of.edge_coupling(w0, np.linalg.eigvalsh(_H(kf, POS0 + d2_)),
                                          np.linalg.eigvalsh(_H(kf, POS0 - d2_)), b, u2)
                    e1, e2_ = abs(math.sqrt(r["D2"]) - abs(Dref)), abs(math.sqrt(r2["D2"]) - abs(Dref))
                    assert e2_ < e1 / 10 or e2_ < 1e-7
            if src == ie and kf == [0.0, 0.0]:                # σh 奇的模在 Γ 点一阶禁戒
                assert all(of.edge_coupling(w0, ep, em, b, u)["D2"] < 1e-10 for b in range(3))
            if src == ia:
                assert max(of.edge_coupling(w0, ep, em, b, u)["D2"] for b in range(3)) > 0.01


# ================================================================ 简并带、近简并、真空对齐
def test_degenerate_manifold_trace(of):
    rng = np.random.default_rng(3)
    d = 3
    V = rng.normal(size=(d, d))
    V = (V + V.T) / 2
    W = rng.normal(size=(d, d)) * 0.3
    W = (W + W.T) / 2
    base = np.array([-2.0, 0.5, 0.5, 0.5, 3.0])
    u = 1e-3
    def ev(x):
        H = np.diag(base)
        H[1:4, 1:4] += x * V + x * x * W
        return np.linalg.eigvalsh(H)
    r = of.edge_coupling(ev(0.0), ev(u), ev(-u), 2, u)
    assert r["d"] == 3 and r["bands"] == [1, 3] and r["method"] == "degenerate(3)"
    assert r["D2"] == pytest.approx(np.trace(V @ V) / d ** 2, rel=1e-4)
    # 真空对齐：两侧真空能级平移 a·u 等价于每条带的一阶项都减 a
    a = 0.7
    r2 = of.edge_coupling(ev(0.0), ev(u) + a * u, ev(-u) - a * u, 2, u, dvac_p=a * u, dvac_m=-a * u)
    assert r2["D2"] == pytest.approx(r["D2"], rel=1e-9)
    # 非简并：真空项直接从 λ 里扣
    e0, ep, em = np.array([0.0, 1.0]), np.array([0.0, 1.0 + 3e-3]), np.array([0.0, 1.0 - 3e-3])
    r3 = of.edge_coupling(e0, ep, em, 1, 1e-3, dvac_p=1e-3, dvac_m=-1e-3)
    assert r3["lam"][0] == pytest.approx(2.0) and r3["D2"] == pytest.approx(4.0)


def _lin(lams, e0, u):
    e0, lams = np.array(e0, float), np.array(lams, float)
    return e0, np.sort(e0 + lams * u), np.sort(e0 - lams * u)


def test_near_band_crossing_is_refused_not_guessed(of):
    """[0035-r1] 邻带近到位移下会交叉：排序配对配错带，给出的 D 随位移变（复审实测 9 / 2.25 / 1.75）。
    现在：没交叉时 D 与位移无关；交叉时报 BandCrossing（要减小 ODP_DISP_A），不硬给一个数。"""
    for u in (0.002, 0.005):                                      # λ=3 的带、20 meV 外 λ=0 的邻带
        r = of.edge_coupling(*_lin([0, 3, 0, 0], [-1, 0, 0.02, 2], u), 1, u)
        assert r["D2"] == pytest.approx(9.0, rel=1e-9) and r["d"] == 1
    assert of.edge_coupling(*_lin([0, 3, 0, 0], [-1, 0, 0.02, 2], 0.005), 1, 0.005)["widened"]
    with pytest.raises(of.BandCrossing):
        of.edge_coupling(*_lin([0, 3, 0, 0], [-1, 0, 0.02, 2], 0.01), 1, 0.01)
    with pytest.raises(of.BandCrossing):                          # 复审的第二个例子：8 meV、λ=±2
        of.edge_coupling(*_lin([0, 2, -2, 0], [0, 1.0, 1.008, 3], 0.01), 1, 0.01)
    r = of.edge_coupling(*_lin([0, 3, 3, 0], [-1, 0, 0.02, 2], 0.01), 1, 0.01)   # SOC 自旋劈裂型：平行移动、不交叉
    assert r["D2"] == pytest.approx(9.0) and r["widened"]
    # 非对角强混合（λ=0、组间耦合 a）：大位移解不出对角元 -> 拒绝；小位移给对角元 0
    def offd(u, a=3.0, dl=0.02):
        H = lambda x: np.array([[0.0, a * x], [a * x, dl]])           # noqa: E731
        ev = [np.r_[-1.0, np.linalg.eigvalsh(H(x)), 2.0] for x in (0.0, u, -u)]
        return of.edge_coupling(ev[0], ev[1], ev[2], 1, u)
    with pytest.raises(of.BandCrossing):
        offd(0.006)
    assert offd(0.002)["D2"] == pytest.approx(0.0, abs=1e-9)


def test_trace_sq_is_ordering_invariant(of):
    rng = np.random.default_rng(5)
    V = rng.normal(size=(3, 3))
    V = (V + V.T) / 2
    u = 1e-3
    e0 = np.array([0.1, 0.13, 0.2])
    ev = [np.linalg.eigvalsh(np.diag(e0) + x * V) for x in (0.0, u, -u)]
    t = of.trace_sq(ev[0], ev[1], ev[2], u)
    assert t == pytest.approx(np.trace(V @ V), rel=1e-6)
    assert of.trace_sq(ev[0], ev[1][::-1], ev[2], u) == pytest.approx(t)


def test_vacuum_level(of):
    n = 200
    z = np.arange(n) / n
    atoms = np.array([0.05, 0.10, 0.95])                     # slab 跨 z=0，真空在中间
    v = np.where((z > 0.2) & (z < 0.85), 4.321, -10.0 * np.exp(-((z - 0.03) % 1.0) ** 2 / 0.002))
    mean, sd = of.vacuum_level(v, atoms)
    assert mean == pytest.approx(4.321) and sd < 1e-12
    ramp = v + np.where((z > 0.2) & (z < 0.85), 0.5 * z, 0.0)  # 面外偶极 -> 真空区斜坡
    with pytest.raises(ValueError):
        of.vacuum_level(ramp, atoms)


def test_mode_table_rejects_non_orthonormal(of):
    f, E = _phonons()
    with pytest.raises(ValueError):
        of.mode_table(f, E / np.sqrt(MASS)[None, :, None], MASS)     # 除过 √m 的那组


def test_mode_table_polar_and_soft_selection(of):
    f, E = _phonons()
    born = np.array([-2 * np.eye(3), np.eye(3), np.eye(3)])
    _m, g_ex, _n = of.mode_table(f, E, MASS, born)
    _m, g_in, _n = of.mode_table(f, E, MASS, born, include_polar=True)
    assert sum(g["use"] for g in g_in) > sum(g["use"] for g in g_ex)
    assert all(("POP" in g["why"]) for g in g_ex if g["polar"])
    _m, g_nb, _n = of.mode_table(f, E, MASS, None)                   # 没有 Born 电荷：一律不算（判不了极性）
    assert not any(g["use"] for g in g_nb)
    _m, g_nbi, _n = of.mode_table(f, E, MASS, None, include_polar=True)   # include：不判极性、全算
    assert all(g["use"] for g in g_nbi)
    _m, g_soft, _n = of.mode_table(f, E, MASS, born, fmin_thz=1e6)
    assert not any(g["use"] for g in g_soft)


def test_combine_members_and_degenerate_pairs(of):
    D, per, spread = of.combine_members({7: [3.0, 1.0], 8: [1.0, 3.0]})   # 简并伙伴：ΣD² 两谷都是 4
    assert D == pytest.approx(2.0) and per == [4.0, 4.0] and spread == 0.0
    D, per, spread = of.combine_members({7: [4.0, 2.0]})
    assert D == pytest.approx(math.sqrt(3.0)) and spread == pytest.approx(2.0 / 3.0)


# ================================================================ odp_params -> S8.4 通道
def _params(D_k=4.0, hw=50.0):
    return {"schema": "odp_params/1", "carriers": {
        "electron": {"status": "ok", "valley_sig": {"use_grid": False, "sig": []}, "valleys": [
            {"label": "K", "groups": [{"group": 3, "D_eV_A": D_k, "hw_meV": hw},
                                      {"group": 1, "D_eV_A": 0.001, "hw_meV": 35.0},
                                      {"group": 4, "D_eV_A": 1.0, "hw_meV": hw + 1e-4}]},
            {"label": "K(up带27)", "groups": [{"group": 3, "D_eV_A": D_k, "hw_meV": hw}]},
            {"label": "G", "groups": [{"group": 3, "D_eV_A": 2.0, "hw_meV": hw}]}]},
        "hole": {"status": "no_valleys"}}}


def test_frozen_channels(of):
    p = _params()
    norm = lambda s: "Γ" if s == "G" else s                             # noqa: E731
    chs, dropped, notes = of.frozen_channels(p, "electron", keep={"K", "Γ"}, norm=norm)
    assert [(c["from"], c["to"], c["kind"]) for c in chs] == [("K", "K", "intra"), ("Γ", "Γ", "intra")]
    k = chs[0]
    assert k["D_eV_A"] == pytest.approx(math.sqrt(16.0 + 1.0)) and k["source"] == "frozen"   # ħω 相同的两组 D² 相加
    assert [d["hw_meV"] for d in dropped] == [35.0] and notes == []      # D < 0.05 丢掉
    assert of.frozen_channels(p, "hole") == ([], [], [])
    # [0035-r1] 被合并的同中心星（另一条带/另一自旋）不再被静默丢掉：与保留的星取 D² 平均
    chs, _d, notes = of.frozen_channels(p, "electron", keep={"K", "Γ"}, norm=norm, alias={"K(up带27)": "K"})
    assert chs[0]["D_eV_A"] == pytest.approx(math.sqrt((17.0 + 16.0) / 2)) and notes == []
    p2 = _params()
    p2["carriers"]["electron"]["valleys"][1]["groups"][0]["D_eV_A"] = 2.0
    chs, _d, notes = of.frozen_channels(p2, "electron", keep={"K", "Γ"}, norm=norm, alias={"K(up带27)": "K"})
    assert chs[0]["D_eV_A"] == pytest.approx(math.sqrt((17.0 + 4.0) / 2)) and len(notes) == 1
    assert of.params_digest(p, "hole") is None
    d0 = of.params_digest(p, "electron")
    assert of.params_digest(_params(D_k=4.0 + 1e-7), "electron") == d0    # 四舍五入以下的抖动不算变
    assert of.params_digest(_params(D_k=4.1), "electron") != d0


def test_odp_consumer_stale(kc, of):
    p = _params()
    dg = of.params_digest(p, "electron")
    st = kc.odp_consumer_stale
    assert st({"frozen_use": {"electron": "used", "hole": "unavailable"}, "digest": {"electron": dg}}, p) == []
    assert st({"frozen_use": {"electron": "used"}, "digest": {"electron": "old"}}, p) == ["electron"]
    assert st({"frozen_use": {"electron": "unavailable"}}, p) == ["electron"]     # 参数刚可用 -> 用上
    assert st({}, p) == ["electron"]                                               # 旧版 S8.4 没记录
    for how in ("manual", "off", "manual_valleys", "not_in_scattering"):
        assert st({"frozen_use": {"electron": how}}, p) == []
    gone = {"schema": "odp_params/1", "carriers": {"electron": {"status": "no_valleys"}}}
    assert st({"frozen_use": {"electron": "used"}, "digest": {"electron": dg}}, gone) == ["electron"]


def test_invalidate_odp_consumer_archives_markers(kc, of, tmp_path, monkeypatch):
    monkeypatch.setattr(kc, "_active_job_dirs", lambda: set())
    d = tmp_path / "step8.4_amset2d"
    d.mkdir()
    (tmp_path / "step8.3_output").mkdir()
    (tmp_path / "step8.3_output" / "comparison_300K.png").write_text("x")
    p = _params()
    rec = {"dp_optical_frozen": {"frozen_use": {"electron": "used"}, "digest": {"electron": of.params_digest(p, "electron")}}}
    (d / "2d_correction.json").write_text(json.dumps(rec))
    (d / "transport.json").write_text("{}")
    (d / "tau_energy.json").write_text("{}")
    assert kc.invalidate_odp_consumer(tmp_path, p) == []                 # 指纹没变：不动
    done = kc.invalidate_odp_consumer(tmp_path, _params(D_k=5.0))
    names = {m for _s, m in done}
    assert {"transport.json", "tau_energy.json", "comparison_300K.png"} <= names
    assert not (d / "transport.json").exists()
    assert list(d.glob("transport.json.stale-upstream-step8.21_odp_read-*"))


# ================================================================ S8.4 gen：ODP_FROZEN
def _gen(monkeypatch, **conf):
    g = _load("s84_0035", S84 / "gen_step14_amset2d.py")
    for k, v in conf.items():
        monkeypatch.setattr(g, k, v)
    monkeypatch.setattr(g, "SCATTERING", list(conf.get("SCATTERING", ["ADP", "IMP", "POP"])))
    return g


def _dpt_and_params(tmp_path, kc, shift=False):
    K = [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]
    vs = [{"label": "K", "g": 2, "dE_eV": 0.0, "k_frac": K, "k_grid": K, "spin": "up", "band_index": 9}]
    (tmp_path / "step8.2_dpt").mkdir(parents=True, exist_ok=True)
    (tmp_path / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": "electron", "multi_valley": {"status": "ok", "valleys": vs}},
        {"carrier": "hole", "multi_valley": {"status": "ok", "valleys": vs}}]}))
    vs_p = [dict(vs[0], k_grid=[[0.3, 0.3, 0], [-0.3, -0.3, 0]])] if shift else vs
    sig = {"use_grid": True, "sig": json.loads(json.dumps(kc._valley_sig(vs_p, True)))}
    p = {"schema": "odp_params/1", "carriers": {
        "electron": {"status": "ok", "valley_sig": sig, "valleys": [
            {"label": "K", "groups": [{"group": 2, "D_eV_A": 4.2, "hw_meV": 49.5},
                                      {"group": 1, "D_eV_A": 0.0, "hw_meV": 35.0}]}]},
        "hole": {"status": "ok", "valley_sig": sig, "valleys": [
            {"label": "K", "groups": [{"group": 2, "D_eV_A": 2.5, "hw_meV": 49.5}]}]}}}
    (tmp_path / "step8.21_odp_read").mkdir(exist_ok=True)
    (tmp_path / "step8.21_odp_read" / "odp_params.json").write_text(json.dumps(p))
    return p


def test_gen_frozen_used(tmp_path, monkeypatch, kc, capsys):
    _dpt_and_params(tmp_path, kc)
    g = _gen(monkeypatch)
    rec = g.resolve_dp_optical(tmp_path)
    assert g.SCATTERING == ["ADP", "IMP", "POP", "ODP"]
    chs = rec["electron"]["channels"]
    assert len(chs) == 1 and chs[0]["kind"] == "intra" and chs[0]["D_eV_A"] == pytest.approx(4.2)
    assert rec["electron"]["frozen"]["dropped_below_dmin"][0]["hw_meV"] == 35.0
    assert rec["adp_intravalley"] == {"electron": False, "hole": False}       # 只有谷内 ODP：ADP 不掩码
    fr = g._DP_OPTICAL_FROZEN_REC
    assert fr["frozen_use"] == {"electron": "used", "hole": "used"} and set(fr["digest"]) == {"electron", "hole"}
    json.dumps(rec)
    json.dumps(fr)
    assert "冻声子 ODP" in capsys.readouterr().out


def test_gen_frozen_manual_wins_and_ivs_combines(tmp_path, monkeypatch, kc):
    _dpt_and_params(tmp_path, kc)
    g = _gen(monkeypatch, DP_OPTICAL_E="K>K:3.0@50", DP_OPTICAL_H="K>K':1.5@45")
    rec = g.resolve_dp_optical(tmp_path)
    fr = g._DP_OPTICAL_FROZEN_REC["frozen_use"]
    assert fr == {"electron": "manual", "hole": "used"}
    assert [c["D_eV_A"] for c in rec["electron"]["channels"]] == [3.0]                       # 手填谷内优先
    kinds = sorted(c["kind"] for c in rec["hole"]["channels"])
    assert kinds == ["inter_star", "intra"] and "IVS" in g.SCATTERING                        # 手填谷间 + 冻声子谷内
    assert rec["adp_intravalley"]["hole"] is True


def test_gen_frozen_off_explicit_and_stale(tmp_path, monkeypatch, kc):
    _dpt_and_params(tmp_path, kc)
    g = _gen(monkeypatch, ODP_FROZEN="off")
    assert g.resolve_dp_optical(tmp_path) is None and g.SCATTERING == ["ADP", "IMP", "POP"]
    assert g._DP_OPTICAL_FROZEN_REC["frozen_use"] == {"electron": "off", "hole": "off"}
    g = _gen(monkeypatch, SCATTERING=["ADP", "IMP"], _SCATTERING_EXPLICIT=True)
    assert g.resolve_dp_optical(tmp_path) is None and g.SCATTERING == ["ADP", "IMP"]
    assert set(g._DP_OPTICAL_FROZEN_REC["frozen_use"].values()) == {"not_in_scattering"}
    assert g._DP_OPTICAL_FROZEN_REC["digest"] == {}
    g = _gen(monkeypatch, DP_OPTICAL_VALLEYS_E="K:0.3333,0.3333,0|-0.3333,-0.3333,0")
    g.resolve_dp_optical(tmp_path)
    assert g._DP_OPTICAL_FROZEN_REC["frozen_use"]["electron"] == "manual_valleys"
    _dpt_and_params(tmp_path, kc, shift=True)                         # 参数按旧谷表算的
    g = _gen(monkeypatch)
    with pytest.raises(SystemExit):
        g.resolve_dp_optical(tmp_path)


def test_gen_frozen_on_requires_file(tmp_path, monkeypatch):
    g = _gen(monkeypatch, ODP_FROZEN="on")
    with pytest.raises(SystemExit):
        g.resolve_dp_optical(tmp_path)
    g = _gen(monkeypatch)                                              # auto 没有文件：与以前完全一样
    assert g.resolve_dp_optical(tmp_path) is None
    assert g._DP_OPTICAL_FROZEN_REC["frozen_use"] == {"electron": "unavailable", "hole": "unavailable"}
    g = _gen(monkeypatch, ODP_FROZEN="maybe")
    with pytest.raises(SystemExit):
        g.resolve_dp_optical(tmp_path)


def test_gen_chain_and_record_source():
    src = (S84 / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert "tau_energy.json mesh_*.h5; " in src                         # 链首清旧 τ
    i = src.index("python tau_energy.py ||")
    # [0035-r1] τ 在 transport.json 落成之后跑：它超时/被杀不会丢掉 AMSET 结果
    assert src.index("_acmd = _acmd + AMSET_TAIL") < i
    assert src.index("python postprocess_intrinsic.py --check-reproduce") < i
    assert '"dp_optical_frozen": _DP_OPTICAL_FROZEN_REC' in src
    assert "_install_tau_energy(out)" in src


# ================================================================ τ(E) 纯函数
def test_bin_tau_and_total(te):
    x = np.array([0.001, 0.004, 0.006, 0.012, 0.6])
    r = np.array([[1e13, 3e13, 2e13, 5e13, 1e13], [1e13, 1e13, 1e13, 1e13, 1e13]])
    tau, tot, cnt = te.bin_tau(x, r, np.array([0.0, 0.005, 0.01, 0.02]))
    assert cnt.tolist() == [2, 1, 1]
    assert tau[0, 0] == pytest.approx(1 / 2e13) and tau[0, 1] == pytest.approx(1 / 2e13)
    assert tot[0] == pytest.approx(1 / 3e13) and tot[2] == pytest.approx(1 / 6e13)
    tau, tot, cnt = te.bin_tau(x[:1], r[:, :1], np.array([0.0, 0.005, 0.01]))
    assert np.isnan(tau[0, 1]) and np.isnan(tot[1])


def test_tau_eff_and_sigma(te):
    mu = np.zeros((2, 3, 3))
    mu[:, 0, 0], mu[:, 1, 1] = [100.0, 50.0], [200.0, 50.0]
    mc = np.zeros((2, 3, 3))
    mc[:, 0, 0], mc[:, 1, 1] = [10.0, 10.0], [10.0, 10.0]
    t = te.tau_eff(mu, mc, tau0=1e-14)
    assert t["xx"].tolist() == pytest.approx([1e-13, 5e-14]) and t["inplane"][0] == pytest.approx(1.5e-13)
    assert te.tau_sigma([2.0, 2.0], [1.0, 3.0], [0.5, 0.1]) == pytest.approx(2.0)
    assert math.isnan(te.tau_sigma([1.0], [0.0], [1.0]))


def test_carrier_states_matches_amset_band_choice(te):
    e = {"up": np.array([[-1.0, -0.9], [-0.2, -0.1], [0.5, 0.6], [1.5, 1.4]])}
    sel, edge = te.carrier_states(e, {"up": 1}, True)
    assert sel["up"].tolist() == [2, 3] and edge == 0.5
    sel, edge = te.carrier_states(e, {"up": 1}, False)
    assert sel["up"].tolist() == [0, 1] and edge == -0.1


def test_fermi_weight_is_minus_df_de_shape(te):
    T = 300.0
    kT = te.KB_HA * T
    e = np.array([-5, 0, 5]) * kT
    w = te.fermi_weight(e, 0.0, T)
    f = 1 / (np.exp(e / kT) + 1)
    assert w == pytest.approx(f * (1 - f), rel=1e-12)


# ================================================================ S8.3 τ 对比（纯 json）
def test_s83_tau_suite(tmp_path, monkeypatch):
    g = _load("s83_0035", S83 / "gen_step13_output.py")
    run = tmp_path / "step8.4_amset2d"
    run.mkdir()
    (run / "transport.json").write_text("{}")
    tj = {"status": "ok", "temperatures_K": [200.0, 300.0], "doping_cm3": [-1e17, 1e17], "carrier": ["electron", "hole"],
          "mechanisms": ["ADP", "ODP"],
          "energy_bins_eV": {"edges": [0, 0.01, 0.02], "centers": [0.005, 0.015]},
          "tau_E_fs": {"ADP": [[[30, 31], [20, 21]]] * 2, "ODP": [[[300, None], [200, 100]]] * 2},
          "tau_E_total_fs": [[[27, 31], [18, 17]]] * 2,
          "tau_eff_fs": {k: {d: [[30.0, 20.0], [40.0, 25.0]] for d in ("xx", "yy", "inplane")}
                         for k in ("overall", "ADP", "ODP")},
          "tau_sigma_fs": [[29.0, 19.5], [39.0, 24.0]], "check": {"tau_sigma_vs_eff_max_rel": 0.03},
          "transport_fingerprint": {"sha256": g._sha256_file(run / "transport.json"), "file": "transport.json"}}
    (run / "tau_energy.json").write_text(json.dumps(tj))
    monkeypatch.setattr(g, "_dpt_tau_fs", lambda cwd, T: {"electron": (20.0, 22.0), "hole": (None, None)})
    monkeypatch.setattr(g, "_geom_2d", lambda cwd: {"c_A": 20.0})
    out = tmp_path / "out"
    out.mkdir()
    files = g.write_tau_suite(out, tmp_path, True)
    assert "tau_comparison.csv" in files and "tau_summary.txt" in files
    rows = (out / "tau" / "tau_comparison.csv").read_text().splitlines()
    assert len(rows) == 1 + 4
    head = rows[0].split(",")
    r300e = dict(zip(head, rows[3].split(",")))                          # 行序：温度 × 掺杂档（200 K 两行在前）
    assert r300e["carrier"] == "electron" and float(r300e["tau_DPT_inplane_fs"]) == pytest.approx(21.0)
    assert float(r300e["ADP_over_DPT"]) == pytest.approx(20.0 / 21.0, rel=1e-4)
    assert float(r300e["n_2D_cm-2"]) == pytest.approx(1e17 * 20e-8)
    # 不是同一次运行的 τ：不出（并清掉上一轮的表）
    (run / "transport.json").write_text('{"x": 1}')
    assert g.write_tau_suite(out, tmp_path, True) == []
    assert not (out / "tau" / "tau_comparison.csv").exists()


# ================================================================ skill.yaml / 接线
def _steps():
    sk = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    return sk, {s["name"]: s for grp in sk["optional_steps"].values() for s in (grp or {}).get("steps") or []}


def test_skill_yaml_group_and_needs():
    sk, by = _steps()
    grp = sk["optional_steps"]["odp_frozen"]
    assert grp["default"] is False and "auto_dim" not in grp            # 默认关：每材料多一批 VASP 单点
    s72, s821, s84 = by["step7c_odp_frozen"], by["step8.21_odp_read"], by["step8.4_amset2d"]
    assert s72["fanout"] == "odp-*" and s72["check"] == "outcar" and "odp_frozen.py" in s72["gen_need"]
    assert s72["needs"] == ["step1_opt", "step5_dielect", "step7_deform"]
    assert s821["run"] == "gen" and s821["done_marker"] == "odp_params.json"
    assert s821["needs"] == ["step7c_odp_frozen", "step8.2_dpt"] and "odp_frozen.py" in s821["gen_need"]
    assert "step8.21_odp_read?" in s84["needs"]
    assert {"odp_frozen.py", "tau_energy.py"} <= set(s84["gen_need"])
    for s in (s72, s821):
        assert (KE / s["src"] / s["gen"]).is_file()
    gens = []
    for grp in sk["optional_steps"].values():
        gens += [str(s["gen"]).split()[0] for s in (grp or {}).get("steps") or [] if s.get("gen")]
    gens += [str(s["gen"]).split()[0] for s in sk["steps"] if s.get("gen")]
    assert len(gens) == len(set(gens))                                   # gen 都推到材料根目录：不许重名


def test_ke_common_tables(kc):
    assert ("step7c_odp_frozen", "dir") in kc.DOWNSTREAM["step1_opt"]
    assert kc.DOWNSTREAM["step7c_odp_frozen"] == [("step8.21_odp_read", "always")]
    assert ("step8.21_odp_read", "always") in kc.DOWNSTREAM["step8.2_dpt"]
    assert kc.DONE_MARKERS["step8.21_odp_read"] == ("odp_params.json",)
    assert "tau_energy.json" in kc.DONE_MARKERS["step8.4_amset2d"]
    assert ("step8.21_odp_read", "step8.4_amset2d") in kc.DOWNSTREAM_EXEMPT
    assert "step8.21_odp_read" not in kc.DOWNSTREAM                         # S8.4 只按指纹有条件失效
    assert ("step7c_odp_frozen", ("odp-ref/POSCAR",)) in kc.LINEAGE_POSCARS
    assert ("step8.21_odp_read", "always") in kc.DOWNSTREAM["step5_dielect"]          # [0035-r1] S5 重算 -> 重读核对


def test_gen_sources_wiring():
    g72 = (S72 / "gen_step9c_odp.py").read_text(encoding="utf-8")
    assert "kc.invalidate_downstream(cwd, OUTDIR_NAME" in g72
    assert g72.index("kc.reexec_with_pymatgen(__file__)") < g72.index("out.mkdir(exist_ok=True)")
    for n in ("INCAR", "KPOINTS", "POTCAR"):
        assert '"%s"' % n in g72
    assert '"_stale_%s_%s"' in g72                                       # 不再算的子目录移出扇出范围
    g821 = (S72 / "odp_read" / "gen_step12c_odp_read.py").read_text(encoding="utf-8")
    assert "kc.invalidate_odp_consumer(cwd, params)" in g821
    assert g821.index("check_inputs(cwd, man)") < g821.index("carrier_params(")
