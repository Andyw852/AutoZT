"""[0061-m] ADP_DEGEN 的矩阵模式：应变算符直接取 DFT 形变势矩阵（deform_wave），不再从本征值反推。

  · 合成模型（test_0061_adp_degen.Synth）+ 同接口的"网格"（给出合成哈密顿量在固定分数坐标下的精确 ∂H/∂ε）：
    setup 用它得到的应变算符与真模型一致、比值表对上；S7 网格太粗、拟合分不开符号的情形，矩阵模式直接定下来；
    逐对替换用的 Zk 在 k0 处等于 Z、离开 k0 时对角部分 = Z 的对角 − 运动学项 T（与拟合模式同一口径）。
  · EPM 端到端：S7 写 WAVECAR -> extract -> build -> StrainGrid，Γ 点价带顶三重态的 Bᴴ D̂ B 与精确 ∂H/∂ε 对上；
    AMSET 系数用了另一个倒格矢表示时 _grid_shift 认得出来；matrix_grid 的对齐与带号错位时拒绝。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _p in (str(KE), str(Path(__file__).resolve().parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import adp_degen as ad  # noqa: E402
import test_0061_adp_degen as T  # noqa: E402

pytest.importorskip("h5py")


class FakeGrid:
    """Synth 的系数空间里的精确应变算符：X^(P)(k) = Q_d · ∂_δ H_model(k; δ·P) · Q_dᴴ（固定分数坐标，含运动学项）。"""

    def __init__(self, sy, delta=1e-5):
        self.sy, self.delta = sy, delta
        self.pairs = ad.pairs(3)
        self.Qd = sy.Q[:, :sy.d]

    def _X(self, kf, P):
        Pm = np.zeros((3, 3))
        Pm[P[0], P[1]] = Pm[P[1], P[0]] = 1.0
        h = (self.sy._mH(kf, self.delta * Pm) - self.sy._mH(kf, -self.delta * Pm)) / (2 * self.delta)
        return self.Qd @ h @ self.Qd.conj().T

    def apply(self, kf, vecs):
        return np.array([np.asarray(vecs) @ self._X(kf, P).T for P in self.pairs])


@pytest.fixture(scope="module")
def si():
    return T.Synth(T.SI_B, T._oh(), T.dkk(*T.KP_SI), T.bir_pikus(*T.BP_SI), 3, n=12, m=8, nd=24)


def test_synth_matrix_mode_matches_true_model(si):
    tab, rep = ad.setup(si.kf, si.E, si.coeff_fn, si.d, si.B, si.ops, T._cubic_C(), si.sign_h5, si.amp_h5,
                        dense=(si.dense_kf, si.dense_E), dim=3, nang=48, grid=FakeGrid(si))
    v = rep["valleys"][0]
    st = v["strain"]
    assert v["status"] == "ok" and st["source"] == "deformation_matrix.h5"
    assert st["symmetrized"] and st["sym_rel_res"] < 1e-6 and st["cross_block_rel"] < 1e-9
    assert st["k0_rel_res"] < 1e-3
    assert st["n_zk"] > 0 and "candidates" not in st
    T._check(tab, T.dkk(*T.KP_SI), T.bir_pikus(*T.BP_SI), 3, 3, T._cubic_C(), 48, tol=1e-3)
    # 子空间里的本征值与真 Bir–Pikus 算符一致（基无关）
    Y = ad.tensor_Y(tab["valleys"][0]["Z"], 3)
    Yt = T.bir_pikus(*T.BP_SI)
    for a, b in ad.pairs(3):
        assert np.allclose(np.linalg.eigvalsh(Y[a, b]), np.linalg.eigvalsh(Yt[a, b]), atol=1e-4)


def test_synth_matrix_mode_resolves_ambiguous_sign():
    """S7 网格太粗（4³）：拟合分不开 d 的符号（test_0061_adp_degen.test_setup_ambiguous_sign）；矩阵模式直接定，不取平均。"""
    sy = T.Synth(T.SI_B, T._oh(), T.dkk(*T.KP_SI), T.bir_pikus(*T.BP_SI), 3, n=12, m=4, nd=24)
    tab, rep = ad.setup(sy.kf, sy.E, sy.coeff_fn, sy.d, sy.B, sy.ops, T._cubic_C(), sy.sign_h5, sy.amp_h5,
                        dense=(sy.dense_kf, sy.dense_E), dim=3, nang=24, grid=FakeGrid(sy))
    v = rep["valleys"][0]
    assert v["status"] == "ok" and v["strain"]["source"] == "deformation_matrix.h5" and "pick" not in v["strain"]
    T._check(tab, T.dkk(*T.KP_SI), T.bir_pikus(*T.BP_SI), 3, 3, T._cubic_C(), 24, tol=1e-3)


def test_synth_zk_kinematics(si):
    """Zk（逐对替换用）：k0 处 = Z；离开 k0 时在模型本征态上的对角元 = Z 的对角元 − ε:T_n（T = sym(δk⊗∇E)），
    即矩阵模式不再减 T 时与拟合模式的带内项一致。"""
    tab, _rep = ad.setup(si.kf, si.E, si.coeff_fn, si.d, si.B, si.ops, T._cubic_C(), si.sign_h5, si.amp_h5,
                         dense=(si.dense_kf, si.dense_E), dim=3, nang=24, grid=FakeGrid(si))
    v = tab["valleys"][0]
    j0 = v["zk_keys"][ad._key(np.zeros(3))]
    assert np.allclose(v["zk"][j0], v["Z"], atol=1e-6)
    hits = 0
    for kk, j in v["zk_keys"].items():
        kf = np.array(kk, float) * ad.KEY_RES
        dk = ad._wrap(kf) @ si.B
        r = np.linalg.norm(dk)
        if not 0.3 * v["rcut"] < r < 0.6 * v["rcut"]:
            continue
        e, V = np.linalg.eigh(ad.model_H(v["model"], dk))
        if np.min(np.diff(e)) < 1e-3:
            continue
        grad = np.real(np.einsum("in,aij,jn->na", V.conj(), ad.model_dH(v["model"], dk), V))
        for t, (a, b) in enumerate(ad.pairs(3)):
            Pm = np.zeros((3, 3))
            Pm[a, b] = Pm[b, a] = 1.0
            Xk = ad.pattern_op(v["zk"][j], 3, (a, b))
            X0 = ad.pattern_op(v["Z"], 3, (a, b))
            Tn = np.einsum("ab,b,na->n", Pm, dk, grad)
            lhs = np.real(np.einsum("in,ij,jn->n", V.conj(), Xk, V))
            rhs = np.real(np.einsum("in,ij,jn->n", V.conj(), X0, V)) - Tn
            assert np.allclose(lhs, rhs, atol=2e-3 * np.max(np.abs(X0)))
        hits += 1
        if hits >= 6:
            break
    assert hits >= 3


def test_synth_grid_failure_falls_back_to_fit(si):
    class Broken(FakeGrid):
        def apply(self, kf, vecs):
            return None
    tab, rep = ad.setup(si.kf, si.E, si.coeff_fn, si.d, si.B, si.ops, T._cubic_C(), si.sign_h5, si.amp_h5,
                        dense=(si.dense_kf, si.dense_E), dim=3, nang=24, grid=Broken(si))
    st = rep["valleys"][0]["strain"]
    assert rep["valleys"][0]["status"] == "ok" and "matrix_fallback" in st and st.get("source") != "deformation_matrix.h5"
    assert "zk" not in tab["valleys"][0]
    assert "本征值拟合（矩阵不可用" in ad.summary_line("up", dict(rep["valleys"][0], carrier="hole"))


# ---------------------------------------------------------------------------
# EPM 端到端
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def epm(tmp_path_factory):
    import _epm_si as E
    import deform_wave as DW
    import test_0061m_deform_wave as TW
    root, dirs, (kpts, dat) = TW._s7(tmp_path_factory.mktemp("mm"), n=4)
    out = root.parent / "deformation_matrix.h5"
    DW.build(root, dirs, out)
    gk = np.concatenate([E.mesh(4, tr_reduce=False), [[0.5, 0, 0], [-0.5, 0, 0]]])
    gpts = np.unique(np.concatenate([E.gvecs(k, E.L0, TW.ENCUT) for k in gk]
                                    + [E.gvecs(k, E.L0, TW.ENCUT) + s for k in gk for s in ([1, 0, 0], [-1, 0, 0])]),
                     axis=0)
    return E, DW, TW, out, kpts, gpts


def _dense(gpts, G, v):
    gl = {tuple(g): i for i, g in enumerate(gpts)}
    out = np.zeros(len(gpts), complex)
    out[[gl[tuple(g)] for g in G]] = v
    return out


def test_epm_gamma_triplet_operator(epm):
    E, DW, TW, out, _kpts, gpts = epm
    grid = DW.StrainGrid(DW.load_matrix(out, "up"), gpts)
    G, _e, V = E.solve(np.zeros(3), E.L0, TW.ENCUT)
    vb = [E.NOCC - 3, E.NOCC - 2, E.NOCC - 1]
    B = np.array([_dense(gpts, G, V[b]) for b in vb])
    g0, cov = ad._grid_shift(grid, B, np.zeros(3))
    assert cov > 0.999 and not np.any(g0)
    Z = ad._grid_Z(grid, B, np.zeros(3), 3)
    Vb = V[vb]
    for t, (a, b) in enumerate(ad.pairs(3)):
        Xe = Vb.conj() @ E.dH(np.zeros(3), G, (a, b)) @ Vb.T / (1.0 if a == b else ad.SQ2)
        assert np.allclose(Z[t], Xe, atol=1e-2 * np.max(np.abs(Xe)))
    # 剪切把三重态耦合起来（Bir–Pikus d ≠ 0），不是只剩对角
    t = ad.pairs(3).index((0, 1))
    assert np.max(np.abs(Z[t] - np.diag(np.diag(Z[t])))) > 0.2


def test_epm_g_shift_detected(epm):
    """X 点 (½,0,0) 的态写成 (−½,0,0) 的表示（c′(G) = c(G − e1)）：按 (½,0,0) 投影低于 COVER_MIN（会被拒），_grid_shift 认出 g0 = (−1,0,0)。"""
    E, DW, TW, out, _kpts, gpts = epm
    grid = DW.StrainGrid(DW.load_matrix(out, "up"), gpts)
    k = np.array([0.5, 0.0, 0.0])
    G, _e, V = E.solve(k, E.L0, TW.ENCUT)
    m = DW.load_matrix(out, "up")
    nb = [int(b) for b in m["bands"]][:3]
    Bshift = np.array([_dense(gpts, G + np.array([1, 0, 0]), V[b]) for b in nb])
    assert np.min(grid.coverage(k, Bshift)) < ad.COVER_MIN
    g0, cov = ad._grid_shift(grid, Bshift, k)
    assert cov > 0.99 and list(g0) == [-1, 0, 0]


def test_epm_matrix_grid_alignment(epm, tmp_path):
    E, DW, TW, out, kpts, gpts = epm
    NB = TW.NB
    S = np.zeros((NB, len(kpts), 3, 3))
    for ip, P in enumerate(E.PATS):
        ep = np.array([E.solve(k, E.lattice(E.pattern(P, TW.DELTA)), TW.ENCUT)[1][:NB] for k in kpts])
        em = np.array([E.solve(k, E.lattice(E.pattern(P, -TW.DELTA)), TW.ENCUT)[1][:NB] for k in kpts])
        S[:, :, P[0], P[1]] = S[:, :, P[1], P[0]] = (-(ep - em) / (2 * TW.DELTA) + TW.REF_SHIFT[ip]).T
    sgn = {"S": S, "kpoints": kpts, "delta": TW.DELTA, "ibands": np.arange(NB)}
    grid, info = ad.matrix_grid(out, "up", sgn, 0, gpts)
    assert grid is not None and info["status"] == "ok" and info["n_align_states"] > 50
    ref = {"%d%d" % P: TW.REF_SHIFT[i] for i, P in enumerate(E.PATS)}
    assert all(abs(info["align_eV"][k] - ref[k]) < 2e-2 for k in ref)
    # ibands 说带轴从 DFT 带 1 开始（错位一条）：以 ibands 为准、离散度变大 -> 拒绝
    bad, info2 = ad.matrix_grid(out, "up", dict(sgn, ibands=np.arange(1, NB + 1)), 0, gpts)
    assert bad is None and info2["band_offset"] == 1 and "band_offset_note" in info2
