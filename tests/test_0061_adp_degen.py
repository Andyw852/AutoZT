"""[0061] 简并价带顶（立方 Γ 三重态）的 ADP 带间修正：adp_degen.py + amset_dpo3d 接线 + S8 gen。

物理核对：
  · 约定：H_xy = n ε_xy（张量分量）；Γ 点 ε_ii 劈裂 {l, m, m}、ε_xy 劈裂 {±n, 0}（本地 QE 7.4 PBE Si 离子弛豫实测
    |n|/√3 = 5.46 eV，对上文献 |d| ≈ 4.85–5.3 eV）。
  · forward_signed（k 点随应变移动那一项）对照"应变后 k·p 直接对角化"的有限差分。
  · 只有静水形变势时完整矩阵元 = AMSET 写法（比值 ≡ 1）；Si 参数下比值 > 1。
  · 模型迁移率的绝对值对 Bardeen–Shockley 解析式。
  · 运行期 setup：合成的 AMSET 数据 + S7.1 形变势文件 -> 拟合参数回到真值，放大系数只作用在 Γ 附近的三重态上。
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
KE = ROOT / "skill" / "ke-dft-cpu"
for _p in (str(KE), str(KE / "step8.4_amset2d"), str(ROOT / "skill" / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import adp_degen as ad  # noqa: E402

A_SI = 5.431
SI_LAT = 0.5 * A_SI * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
SI_B = 2 * np.pi * np.linalg.inv(SI_LAT).T
KP_SI = (-5.568, -3.596, -8.457)                 # 本地 QE PBE（无 SOC）拟合
BP_SI = (1.74 + 2 * -2.31, 1.74 + 2.31, -9.46)   # l = a + 2b、m = a − b、n（离子弛豫）


def _cubic(C11=165.7, C12=63.9, C44=79.6):
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if i == j and k == l:
                        C[i, j, k, l] = C11 if i == k else C12
                    elif ((i == k and j == l) or (i == l and j == k)) and i != j:
                        C[i, j, k, l] = C44
    return C * 1e9


def _gamma_arrays(l, m, n, delta=0.005, offdiag="tensor"):
    """AMSET/S7.1 在 Γ 点会得到的值：每个应变 ±δ 下按能量排序的三条带 -> 带符号中心差分 S、幅值 A（原生约定）。"""
    Th = ad.strain_operator(l, m, n)
    S, A = np.zeros((3, 3, 3)), np.zeros((3, 3, 3))
    for i in range(3):
        for j in range(3):
            P = np.zeros((3, 3))
            P[i, j] = P[j, i] = 1.0
            op = np.einsum("ab,abxy->xy", P, Th)
            sp = -np.linalg.eigvalsh(op * delta) / delta
            sm = np.linalg.eigvalsh(-op * delta) / delta
            S[:, i, j] = 0.5 * (sp + sm)
            A[:, i, j] = 0.5 * (np.abs(sp) + np.abs(sm))
    if offdiag == "tensor":
        A = A * (np.eye(3) + 0.5 * (1 - np.eye(3)))
    return S, A


# ---------------------------------------------------------------------------
# 约定与参数提取
# ---------------------------------------------------------------------------
def test_strain_operator_convention():
    Th = ad.strain_operator(-2.0, 4.0, -8.0)
    eps = np.zeros((3, 3))
    eps[0, 1] = eps[1, 0] = 0.01
    H = np.einsum("ijab,ij->ab", Th, eps)
    assert H[0, 1] == pytest.approx(-0.08) and H[1, 0] == pytest.approx(-0.08)          # H_xy = n ε_xy
    assert np.allclose(np.linalg.eigvalsh(H), [-0.08, 0, 0.08])
    eps = np.diag([0.01, 0, 0])
    assert np.allclose(np.diag(np.einsum("ijab,ij->ab", Th, eps)), [-0.02, 0.04, 0.04])  # {l, m, m}·ε_xx
    c = np.array([1.0, 1.0, 1.0]) / math.sqrt(3)                                        # [111] 单态
    D = ad.diag_potentials(c, Th)
    assert D[0, 1] == pytest.approx(-8.0 / 3) and np.trace(D) == pytest.approx(-2.0 + 2 * 4.0)


@pytest.mark.parametrize("l,m,n", [(-2.88, 4.05, -9.46), (6.0, 1.0, 5.0), (-1.0, -3.0, 2.0)])
@pytest.mark.parametrize("offdiag", ["tensor", "native"])
def test_gamma_bir_pikus(l, m, n, offdiag):
    S, A = _gamma_arrays(l, m, n, offdiag=offdiag)
    gl, gm, gn, info = ad.gamma_bir_pikus(S, A, offdiag)
    assert gl == pytest.approx(l, abs=1e-9) and gm == pytest.approx(m, abs=1e-9) and gn == pytest.approx(abs(n), abs=1e-9)
    assert info["outer_asym_eV"] < 1e-9 and info["shear_mid_eV"] < 1e-9


@pytest.mark.parametrize("N", [-8.457, 6.0])
def test_fit_kp_recovers(N):
    rng = np.random.default_rng(1)
    L, M = -5.568, -3.596
    d = rng.normal(size=(40, 3))
    k = d / np.linalg.norm(d, axis=1)[:, None] * rng.uniform(0.03, 0.12, size=(40, 1))
    e = np.linalg.eigvalsh(ad.kp_matrix(k, L, M, N))
    kk = np.repeat(k, 3, axis=0)
    ee = e.ravel()
    rr = np.tile([0, 1, 2], 40)
    fL, fM, fN, info = ad.fit_kp(kk, ee, rr)
    assert (fL, fM, fN) == pytest.approx((L, M, N), rel=1e-4) and info["rel_rms"] < 1e-6


def test_forward_signed_matches_strained_kp():
    """forward_signed（一阶，含 k 点随应变移动项）vs 应变后的 k·p（k′ = (1+ε)^{-T} k，加 Θ(ε)）直接对角化的中心差分。"""
    L, M, N = KP_SI
    l, m, n = BP_SI
    Th = ad.strain_operator(l, m, n)
    k = np.array([[0.05, 0.03, 0.08], [0.07, -0.02, 0.01]])
    delta = 1e-4
    pred = ad.forward_signed(k, L, M, N, l, m, n)
    for i in range(3):
        for j in range(3):
            P = np.zeros((3, 3))
            P[i, j] = P[j, i] = 1.0
            fd = []
            for sg in (1, -1):
                eps = sg * delta * P
                kp = k @ np.linalg.inv(np.eye(3) + eps)          # 行矢量：k′ = (1+ε)^{-T} k（ε 对称）
                E = np.linalg.eigvalsh(ad.kp_matrix(kp, L, M, N) + np.einsum("ab,abxy->xy", eps, Th))
                fd.append(-(E - np.linalg.eigvalsh(ad.kp_matrix(k, L, M, N))) / (sg * delta))
            assert np.allclose(0.5 * (fd[0] + fd[1]), pred[:, :, i, j], rtol=2e-3, atol=2e-3)   # 近简并点有 O(δ²)


@pytest.mark.parametrize("sign", [1, -1])
def test_sign_of_n(sign):
    L, M, N = KP_SI
    l, m, n = BP_SI
    rng = np.random.default_rng(7)
    d = rng.normal(size=(12, 3))
    k = d / np.linalg.norm(d, axis=1)[:, None] * 0.1
    S = ad.forward_signed(k, L, M, N, l, m, sign * abs(n)) + rng.normal(scale=0.2, size=(12, 3, 3, 3))
    sg, conf, nv = ad.sign_of_n(k, S, np.ones((12, 3), bool), L, M, N, l, m, abs(n))
    assert sg == sign and conf > 0.5 and nv == 12 * 3 * 3


# ---------------------------------------------------------------------------
# 散射率之比与模型迁移率
# ---------------------------------------------------------------------------
def test_ratio_hydrostatic_is_one_and_si_above_one():
    C = _cubic()
    t = ad.rate_ratio(*KP_SI, 1.74, 1.74, 0.0, C, nang=48)
    assert np.allclose(t["ratio"], 1.0, atol=1e-9)
    t = ad.rate_ratio(*KP_SI, *BP_SI, C, nang=48)
    assert np.all(t["ratio"] > 0.8) and float(np.mean(t["ratio"])) > 1.3


def test_model_mobility_bardeen_shockley():
    """各向同性单抛物带（L=M、N=0 -> m* = m0/2，三重简并）+ 只有静水 Ξ + 各向同性弹性：对 Bardeen–Shockley。"""
    C11, C12 = 165.7e9, 63.9e9
    C = _cubic(165.7, 63.9, (165.7 - 63.9) / 2)
    mu = ad.model_mobility(-2.0, -2.0, 0.0, 3.0, 3.0, 0.0, C, nang=80)
    hbar, m0, e, kB = 1.054571817e-34, 9.1093837015e-31, 1.602176634e-19, 1.380649e-23
    ms, kT, Xi = 0.5 * m0, kB * 300.0, 3.0 * e
    ref = 2 * math.sqrt(2 * math.pi) * e * hbar ** 4 * C11 / (3 * ms ** 2.5 * kT ** 1.5 * Xi ** 2) * 1e4
    assert mu == pytest.approx(ref, rel=0.03)


def test_model_mobility_si_full_below_amset():
    C = _cubic()
    full = ad.model_mobility(*KP_SI, *BP_SI, C, nang=60, which="full")
    diag = ad.model_mobility(*KP_SI, *BP_SI, C, nang=60, which="diag")
    assert 1.5 < diag / full < 3.0 and 600 < full < 2500


# ---------------------------------------------------------------------------
# 运行期 setup：合成的 AMSET 数据 + S7.1 形变势文件
# ---------------------------------------------------------------------------
class _Lat:
    def __init__(self, B):
        self.reciprocal_lattice = type("R", (), {"matrix": B})()


class _Struct:
    def __init__(self, B):
        self.lattice = _Lat(B)


def _mesh(n):
    g = np.arange(n) / n
    return np.array([[a, b, c] for a in g for b in g for c in g])


def _bands(kfrac, kp=KP_SI):
    """6 条带：深价带、Γ 三重态（k·p）、两条导带。电子能量（eV，价带顶 = 0）。"""
    kf = kfrac - np.round(kfrac)
    kc = kf @ SI_B
    tri = np.linalg.eigvalsh(ad.kp_matrix(kc, *kp))
    k2 = np.sum(kc * kc, axis=1)
    E = np.vstack([np.full(len(kf), -12.0), tri.T, 1.0 + 3.8 * k2, 1.5 + 3.8 * k2])
    return E, kc


def _write_h5(tmp, kfrac7, bp=BP_SI, kp=KP_SI, delta=0.005):
    import h5py
    E, kc = _bands(kfrac7, kp)
    nb, nk = E.shape
    S = np.zeros((nb, nk, 3, 3))
    iG = int(np.argmin(np.linalg.norm(kc, axis=1)))
    nz = np.arange(nk) != iG
    S[1:4, nz] = ad.forward_signed(kc[nz], *kp, *bp).transpose(1, 0, 2, 3)
    Sg, Ag = _gamma_arrays(*bp, delta=delta, offdiag="native")
    S[1:4, iG] = Sg
    A = np.abs(S)
    A[1:4, iG] = Ag
    A = A * (np.eye(3) + 0.5 * (1 - np.eye(3)))           # 幅值文件：V189 张量约定
    gap = np.full((nb, nk), np.inf)
    d = np.abs(np.diff(E, axis=0))
    gap[:-1] = np.minimum(gap[:-1], d)
    gap[1:] = np.minimum(gap[1:], d)
    with h5py.File(tmp / "deformation.h5", "w") as f:
        f.create_dataset("deformation_potentials_up", data=A)
        f["kpoints"] = kfrac7
        f.attrs["dp_offdiag"] = "tensor"
    with h5py.File(tmp / "deformation_signed.h5", "w") as f:
        f.create_dataset("signed_up", data=S)
        f.create_dataset("gap_up", data=gap)
        f["kpoints"] = kfrac7
        f.attrs["strain_amount"] = delta


@pytest.fixture()
def synth(tmp_path):
    _write_h5(tmp_path, _mesh(8))
    kf = _mesh(24)
    E, _kc = _bands(kf)
    return {"E": E / ad.HARTREE_EV, "kf": kf, "dir": tmp_path}


def test_setup_recovers_parameters(synth):
    t, rep = ad.setup(synth["E"], synth["kf"], 3, _Struct(SI_B), _cubic(), synth["dir"] / "deformation.h5",
                      synth["dir"] / "deformation_signed.h5", nang=48)
    assert rep["status"] == "ok", rep
    assert (rep["kp"]["L"], rep["kp"]["M"], rep["kp"]["N"]) == pytest.approx(KP_SI, rel=1e-3)
    bp = rep["bir_pikus"]
    assert (bp["l"], bp["m"], bp["n"]) == pytest.approx(BP_SI, rel=1e-6)
    assert rep["n_sign"]["sign"] == -1 and rep["n_sign"]["confidence"] > 0.5
    assert min(rep["ratio_mean"]) > 1.0 and t["bands"] == [1, 2, 3]
    k_near = np.array([0.02, 0.01, 0.0])                   # 分数坐标，靠近 Γ
    assert ad.scale(t, 3, k_near) > 1.0 and ad.scale(t, 4, k_near) == 1.0 and ad.scale(t, 0, k_near) == 1.0
    assert ad.scale(t, 3, np.array([0.5, 0.0, 0.0])) == 1.0          # 离 Γ 远：不修正
    assert ad.scale(t, 2, np.zeros(3)) == pytest.approx(float(np.mean(t["ratio"][:, 1])))
    json.dumps(rep, default=float)


def test_setup_not_applicable(synth, tmp_path):
    args = (_Struct(SI_B), _cubic(), synth["dir"] / "deformation.h5", synth["dir"] / "deformation_signed.h5")
    t, rep = ad.setup(synth["E"], synth["kf"], 3, *args, is_soc=True)
    assert t is None and "SOC" in rep["status"]
    E, _ = _bands(synth["kf"])
    E[3, 1000] = 0.3                                      # 价带顶不在 Γ
    t, rep = ad.setup(E / ad.HARTREE_EV, synth["kf"], 3, *args)
    assert t is None and "Γ" in rep["status"]
    E, _ = _bands(synth["kf"])
    E[1:4] += np.array([[0.0], [0.03], [0.06]])           # Γ 不简并
    t, rep = ad.setup(E / ad.HARTREE_EV, synth["kf"], 3, *args)
    assert t is None and "三重简并" in rep["status"]
    t, rep = ad.setup(synth["E"], synth["kf"], 3, *args[:2], tmp_path / "nope.h5", args[3])
    assert t is None and "读不了" in rep["status"]


# ---------------------------------------------------------------------------
# 插件接线（不装 AMSET：最小桩模块）
# ---------------------------------------------------------------------------
_STUBS = {
    "amset/__init__.py": "__version__ = '0.5.1'\n",
    "amset/scattering/__init__.py": "",
    "amset/scattering/basic.py": (
        "class AbstractBasicScattering:\n"
        "    def __init__(self, properties, doping, temperatures, nbands):\n"
        "        self.properties = properties\n"),
    "amset/scattering/elastic.py": (
        "import numpy as np\n"
        "class AcousticDeformationPotentialScattering:\n"
        "    name = 'ADP'\n"
        "    def __init__(self, vb_idx, elastic_constant):\n"
        "        self.vb_idx, self.elastic_constant = vb_idx, elastic_constant\n"
        "    @classmethod\n"
        "    def from_amset_data(cls, mp, amset_data):\n"
        "        return cls(amset_data.vb_idx, mp['elastic_constant'])\n"
        "    def to_reference(self):\n"
        "        return (self.vb_idx, self.elastic_constant)\n"
        "    @classmethod\n"
        "    def from_reference(cls, *ref):\n"
        "        return cls(*ref)\n"
        "    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):\n"
        "        return np.ones(np.shape(norm_q_sq))\n"),
    "amset/scattering/calculate.py": (
        "from amset.scattering.elastic import AcousticDeformationPotentialScattering\n"
        "_scattering_mechanisms = {'ADP': AcousticDeformationPotentialScattering}\n"
        "def scattering_worker(*a, **k):\n"
        "    return 'orig'\n"),
}


class _Spin:
    name = "up"


def test_plugin_scales_adp(synth, tmp_path, monkeypatch):
    run = tmp_path / "run"
    run.mkdir()
    for rel, txt in _STUBS.items():
        (run / rel).parent.mkdir(parents=True, exist_ok=True)
        (run / rel).write_text(txt)
    for f in (KE / "step8_amset" / "amset_dpo3d.py", KE / "step8.4_amset2d" / "dp_optical.py", KE / "adp_degen.py"):
        shutil.copy(f, run)
    for f in ("deformation.h5", "deformation_signed.h5"):
        shutil.copy(synth["dir"] / f, run / f)
    (run / "dp_optical_3d.json").write_text(json.dumps(
        {"mechanisms": [], "dp_optical": None,
         "adp_degen": {"mode": "on", "status": "pending_runtime", "efit_eV": 0.15, "nang": 32,
                       "deform_h5": "deformation.h5", "signed_h5": "deformation_signed.h5"}}))
    monkeypatch.chdir(run)
    monkeypatch.syspath_prepend(str(run))
    for name in [m for m in sys.modules if m == "amset" or m.startswith("amset.") or m in
                 ("amset_dpo3d", "adp_degen", "dp_optical", "amset_ir_fix")]:
        monkeypatch.delitem(sys.modules, name)
    with redirect_stdout(io.StringIO()) as buf:
        import amset_dpo3d as P
        import amset.scattering.calculate as C
    assert C._scattering_mechanisms["ADP"] is P.AcousticDeformation3DIntravalley
    assert "ADP 简并价带修正：开" in buf.getvalue()
    spin = _Spin()
    data = type("AD", (), {"spins": [spin], "energies": {spin: synth["E"]}, "kpoints": synth["kf"],
                           "vb_idx": {spin: 3}, "structure": _Struct(SI_B)})()
    with redirect_stdout(io.StringIO()) as buf:
        adp = P.AcousticDeformation3DIntravalley.from_amset_data({"elastic_constant": _cubic()}, data)
    assert "ADP 简并价带修正（up）" in buf.getvalue() and (run / "adp_degen_fit.json").is_file()
    uq, nq = np.tile([[1.0, 0, 0]], (5, 1)), np.full(5, 0.01)
    k_near = np.array([0.02, 0.01, 0.0])
    f3 = adp.factor(uq, nq, spin, 3, k_near, np.zeros(3))
    f4 = adp.factor(uq, nq, spin, 4, k_near, np.zeros(3))
    assert np.all(f3 > 1.0) and np.allclose(f4, 1.0)
    re = P.AcousticDeformation3DIntravalley.from_reference(*adp.to_reference())
    assert np.allclose(re.factor(uq, nq, spin, 3, k_near, np.zeros(3)), f3)
    # 不适用 + mode on -> 报错
    rec = json.loads((run / "dp_optical_3d.json").read_text())
    data.vb_idx = {spin: 4}
    with pytest.raises(RuntimeError, match="不适用"):
        with redirect_stdout(io.StringIO()):
            P.AcousticDeformation3DIntravalley.from_amset_data({"elastic_constant": _cubic()}, data)
    assert rec["adp_degen"]["mode"] == "on"


# ---------------------------------------------------------------------------
# S8 gen：开关、前置条件、装插件
# ---------------------------------------------------------------------------
def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def g10(tmp_path):
    g = _load("g10_adpdegen", KE / "step8_amset" / "gen_step10_amset.py")
    g.SCATTERING = ["ADP", "IMP", "POP"]
    g.ADP_DEGEN = "auto"
    cwd = tmp_path / "ke"
    (cwd / g.READ_DIR).mkdir(parents=True)
    return g, cwd


def _resolve(g, cwd, is_2d=False):
    buf = io.StringIO()
    with redirect_stdout(buf):
        rec = g.resolve_adp_degen_3d(is_2d, cwd)
    return rec, buf.getvalue()


def test_gen_resolve_modes(g10):
    g, cwd = g10
    assert _resolve(g, cwd, is_2d=True)[0] is None
    rec, log = _resolve(g, cwd)
    assert rec["status"] == "no_signed_h5" and "S7.1" in log
    g.ADP_DEGEN = "on"
    with pytest.raises(SystemExit, match="ADP_DEGEN = on"):
        _resolve(g, cwd)
    g.ADP_DEGEN = "off"
    assert _resolve(g, cwd)[0]["status"] == "off"
    g.ADP_DEGEN = "auto"
    (cwd / g.READ_DIR / "deformation_signed.h5").write_bytes(b"x")
    (cwd / g.READ_DIR / "band_edges.json").write_text(json.dumps({"deform_geometry": "clamped"}))
    rec, log = _resolve(g, cwd)
    assert rec["status"] == "pending_runtime" and rec["deform_geometry"] == "clamped" and "钳位离子" in log
    g.SCATTERING = ["IMP", "POP"]
    assert _resolve(g, cwd)[0]["status"] == "no_adp"
    g.SCATTERING = ["ADP", "IMP", "POP"]
    (cwd / "step3_uniform").mkdir()
    (cwd / "step3_uniform" / "INCAR").write_text("LSORBIT = .TRUE.\n")
    assert _resolve(g, cwd)[0]["status"] == "soc"
    with pytest.raises(SystemExit):
        g.ADP_DEGEN = "maybe"
        _resolve(g, cwd)


def test_gen_installs_plugin_for_degen_only(g10, tmp_path):
    g, cwd = g10
    out = tmp_path / "out"
    out.mkdir()
    rec = {"mechanisms": [], "dp_optical": None, "adp_degen": {"mode": "auto", "status": "pending_runtime"}}
    with redirect_stdout(io.StringIO()):
        assert g._install_dp_optical_3d(out, rec) is True
    for f in ("dp_optical_3d.json", "amset_dpo3d.py", "dp_optical.py", "adp_degen.py"):
        assert (out / f).is_file(), f
    out2 = tmp_path / "out2"
    out2.mkdir()
    with redirect_stdout(io.StringIO()):
        assert g._install_dp_optical_3d(out2, dict(rec, adp_degen={"mode": "auto", "status": "no_signed_h5"})) is False
    assert "adp_degen.py" in (KE / "skill.yaml").read_text() and "ADP_DEGEN" in g.SPEC
