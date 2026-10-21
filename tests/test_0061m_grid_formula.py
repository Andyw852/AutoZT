"""[0061-m] 网格公式：AMSET 自己的系数 × DFT 应变算符，在四面体积分里换全部末态的 ADP 矩阵元。

  · grid_amplitudes（EPM）：⟨m,k′|D̂(k)|n,k⟩ 与精确 ⟨m,k′|∂H(k)/∂ε|n,k⟩ 对上，离格点误差随 S7 网格加密按二阶下降；
    最后一列就是 ⟨f|i⟩。
  · grid_pair_values：按 AMSET 0.5.1 calculate_rate 的中间量（masks、截面插值、g_diff 平移的末态）算 F，
    与"逐顶点 Σ_λ|e_λ·A|²/C_λ 再按截面插值"逐位相等；|⟨f|i⟩|² 那一项与 AMSET 的截面重叠相等（自检），
    自检不过整体关掉；末态在另一个倒格矢表示里也对。
  · matrix_apply：网格可用时全部 q 点换成 F/重叠、乘谷内掩码；不可用时与 [0061-p] 的 pair_apply 逐位相同。
  · tables_to_reference / tables_from_reference：StrainGrid 经共享内存往返后作用相同。
  · install_pair_hooks：0.5.1 形状的 calculate_rate 上装得上，记下初态 / 系数 / masks / 截面插值参数，认出 g_diff。
"""
import sys
import tempfile
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _p in (str(KE), str(Path(__file__).resolve().parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import adp_degen as ad  # noqa: E402
import _epm_si as E  # noqa: E402
import deform_wave as DW  # noqa: E402
import test_0061m_deform_wave as TW  # noqa: E402

pytest.importorskip("h5py")


def _matrix(n):
    root, dirs, _d = TW._s7(Path(tempfile.mkdtemp()) / "s7", n=n)
    out = root.parent / "m.h5"
    DW.build(root, dirs, out, log=lambda m: None)
    return DW.load_matrix(out, "up")


@pytest.fixture(scope="module")
def m6():
    return _matrix(6)


def _gl(gpts):
    return {tuple(g): i for i, g in enumerate(gpts)}


def _dense(gl, n, G, v):
    out = np.zeros(n, complex)
    out[[gl[tuple(g)] for g in G]] = v
    return out


def _amp_errors(m, ks, rng):
    gpts = np.unique(np.concatenate([E.gvecs(k, E.L0, TW.ENCUT + 60) for k in ks]), axis=0)
    gl = _gl(gpts)
    grid = DW.StrainGrid(m, gpts)
    errs = []
    for k in ks:
        G, _e, V = E.solve(k, E.L0, TW.ENCUT)
        kp = k + rng.normal(size=3) * 0.03
        Gp, _ep, Vp = E.solve(kp, E.L0, TW.ENCUT)
        idx = [gl[tuple(g)] for g in G]
        for nb in (1, 3, 4):
            cn = _dense(gl, len(gpts), G, V[nb])
            F = np.array([_dense(gl, len(gpts), Gp, Vp[mb]) for mb in range(6)])
            A = ad.grid_amplitudes(grid, k, cn, F, 3)
            assert np.allclose(A[:, 6], F.conj() @ cn)
            ex = np.array([[F[mb][idx].conj() @ E.dH(k, G, P) @ V[nb] for P in ad.pairs(3)] for mb in range(6)])
            errs.append(np.max(np.abs(A[:, :6] - ex)) / np.max(np.abs(ex)))
    return np.array(errs)


def test_grid_amplitudes_converge(m6):
    rng = np.random.default_rng(5)
    ks = (rng.random((4, 3)) - 0.5) * 0.6
    e4 = _amp_errors(_matrix(4), ks, np.random.default_rng(7))
    e6 = _amp_errors(m6, ks, np.random.default_rng(7))
    assert np.median(e6) < 0.05 and np.median(e4) / np.median(e6) > 1.8


# ---------------------------------------------------------------------------
# AMSET 0.5.1 形状的中间量
# ---------------------------------------------------------------------------
def _g_maps(gpoints):
    """amset.util.get_g_maps 的同一实现：g_maps[i,j,l,idx] = gpoints[idx] + (i−1, j−1, l−1) 的下标（没有 -> −1）。"""
    lo = gpoints.min(0)
    shp = gpoints.max(0) - lo + 3
    grid = -np.ones(shp, int)
    o = -lo + 1
    for i, g in enumerate(gpoints):
        grid[tuple(g + o)] = i
    gm = np.zeros((3, 3, 3, len(gpoints)), int)
    for i, x in enumerate((-1, 0, 1)):
        for j, y in enumerate((-1, 0, 1)):
            for l_, z in enumerate((-1, 0, 1)):
                gm[i, j, l_] = grid[gpoints[:, 0] + x + o[0], gpoints[:, 1] + y + o[1], gpoints[:, 2] + z + o[2]]
    return gm


class _Spin:
    name = "up"


def _scene(m, rng, shift_some=True, nq=30):
    """初态 k（带 3，价带顶三重态里），末态顶点在 k 附近 6 条带 × 8 个 k′（一部分以 k′+e1 的表示存，系数也按那个表示）；
    四面体 = 随机 4 个顶点；细 q 点随机映到四面体。-> (_CTX 内容, 用到的对象)。"""
    k = np.array([0.07, -0.04, 0.03])
    kps = k + rng.normal(size=(8, 3)) * 0.04
    rep = np.zeros((8, 3), int)
    if shift_some:
        rep[::3] = [1, 0, 0]
    allk = np.vstack([k, kps + rep])
    G0 = [E.solve(x, E.L0, TW.ENCUT) for x in np.vstack([k, kps])]
    gpts = np.unique(np.concatenate([g for g, _e, _v in G0] + [g - r for (g, _e, _v), r in zip(G0[1:], rep)]
                                    + [g + s for g, _e, _v in G0 for s in ([1, 0, 0], [-1, 0, 0])]), axis=0)
    gl = _gl(gpts)
    nb, nk = 6, len(allk)
    rows, cmap = [], np.full((nb, nk), 10 ** 6)
    for j in range(nk):
        G, _e, V = G0[j]
        r = rep[j - 1] if j else np.zeros(3, int)
        for b in range(nb):
            cmap[b, j] = len(rows)
            rows.append(_dense(gl, len(gpts), G - r, V[b]))           # 存的表示：c_s(G′) = c(G′ + r)
    coeffs = np.array(rows)
    spin = _Spin()
    bkm = np.zeros((nb, nk), bool)
    bkm[:, 1:] = True
    bm, km = np.where(bkm)
    ntet = 12
    tb = rng.integers(0, nb, size=(ntet, 4))
    tk = rng.integers(1, nk, size=(ntet, 4))
    pm = (tb, tk)
    gm = _g_maps(gpts)
    # AMSET 0.5.1 的 _get_overlap：末态按 −g_diff 平移后与初态点乘
    ci = coeffs[cmap[3, 0]]
    vo = np.zeros((nb, nk))
    for b, kk in zip(bm, km):
        f = coeffs[cmap[b, kk]]
        gd = ((allk[kk] - k) - ad._pbc_diff(allk[kk], k)).astype(int)
        pad = np.concatenate([f, [0]])
        f = pad[gm[1 - gd[0], 1 - gd[1], 1 - gd[2]]]
        vo[b, kk] = abs(f @ ci.conj()) ** 2
    tet_overlap = vo[pm].mean(axis=1)
    mapping = rng.integers(0, ntet, size=nq)
    uq = rng.normal(size=(nq, 3))
    uq /= np.linalg.norm(uq, axis=1)[:, None]
    nq2 = rng.uniform(1e-4, 1e-2, size=nq)
    tbs = types.SimpleNamespace(kpoints=allk, energies={spin: np.zeros((nb, nk))})
    oc = types.SimpleNamespace(g_maps=gm)
    ctx = {"rate": (tbs, oc, {spin: coeffs}, {spin: cmap}, spin, 3, 0), "masks": (pm, bkm, bm, km),
           "tet_contrib": ((), {}), "mapping": mapping, "tet_overlap": tet_overlap, "energy_diff": None}
    return ctx, dict(k=k, gpts=gpts, uq=uq, nq2=nq2, spin=spin, pm=pm, bkm=bkm, bm=bm, km=km, allk=allk,
                     coeffs=coeffs, cmap=cmap, gm=gm, mapping=mapping, tet_overlap=tet_overlap)


def _mean_cs(values, *a, average=True):
    values = np.asarray(values)
    return values.mean(axis=1) if average else values


@pytest.fixture()
def pair_env(monkeypatch):
    monkeypatch.setitem(ad._PAIR, "ok", True)
    monkeypatch.setitem(ad._PAIR, "grid_ok", True)
    monkeypatch.setitem(ad._PAIR, "gdiff", True)
    monkeypatch.setitem(ad._PAIR, "cs", _mean_cs)
    monkeypatch.setitem(ad._PAIR, "grid_off", False)
    return monkeypatch


C_SI = np.einsum("ij,kl->ijkl", np.eye(3), np.eye(3)) * 60.0 + (                # 各向同性弹性张量（只用比例）
    np.einsum("ik,jl->ijkl", np.eye(3), np.eye(3)) + np.einsum("il,jk->ijkl", np.eye(3), np.eye(3))) * 50.0


def test_grid_pair_values_matches_vertex_sum(m6, pair_env):
    rng = np.random.default_rng(3)
    ctx, o = _scene(m6, rng)
    pair_env.setattr(ad, "_CTX", ctx)
    grid = DW.StrainGrid(m6, o["gpts"])
    table = {"grid": grid, "dim": 3, "vb": 3, "grid_signs": [-1], "valleys": []}
    vel = np.array([0.01, -0.02, 0.005])
    res = ad.grid_pair_values(table, 3, o["k"], o["uq"], o["nq2"], C_SI, vel)
    assert res is not None
    Fq, ovl, elig = res
    assert np.all(elig) and np.allclose(ovl, o["tet_overlap"][o["mapping"]])
    # 逐顶点：A（末态按 g_diff 平移）-> Σ_λ |e_λ·A|²/C_λ，四个顶点取平均（= 假的截面插值），再按 mapping 取
    ci = o["coeffs"][o["cmap"][3, 0]]
    Cl, U = ad.acoustic_modes(C_SI, o["uq"], 3)
    Es = 0.5 * (np.einsum("fi,fjl->flij", o["uq"], U) + np.einsum("fj,fil->flij", o["uq"], U))
    e = np.concatenate([np.stack([Es[..., a, b] / ad.HARTREE_EV for a, b in ad.pairs(3)], -1),
                        np.einsum("flij,ij->fl", Es, np.outer(vel, vel))[..., None]], -1)
    tb, tk = o["pm"]
    expF = np.zeros(len(o["uq"]))
    for qi, t in enumerate(o["mapping"]):
        acc = 0.0
        for vtx in range(4):
            b, kk = tb[t, vtx], tk[t, vtx]
            f = o["coeffs"][o["cmap"][b, kk]]
            gd = ((o["allk"][kk] - o["k"]) - ad._pbc_diff(o["allk"][kk], o["k"])).astype(int)
            f = ad.align_coefficients(f, gd, o["gm"])
            A = ad.grid_amplitudes(grid, o["k"], ci, f[None], 3)[0]
            acc += np.sum(np.abs(e[qi] @ A) ** 2 / Cl[qi])
        expF[qi] = acc / 4
    assert np.allclose(Fq, expF, rtol=1e-9)
    assert not ad._PAIR.get("grid_off")
    # 不在网格公式范围：导带（只开了空穴）、换了初态、能量差（非弹性）
    assert ad.grid_pair_values(table, 4, o["k"], o["uq"], o["nq2"], C_SI, vel) is None
    pair_env.setitem(ctx, "energy_diff", 0.01)
    assert ad.grid_pair_values(table, 3, o["k"], o["uq"], o["nq2"], C_SI, vel) is None


def test_grid_self_check_turns_off(m6, pair_env):
    ctx, o = _scene(m6, np.random.default_rng(4))
    ctx["tet_overlap"] = o["tet_overlap"] * 1.01
    pair_env.setattr(ad, "_CTX", ctx)
    table = {"grid": DW.StrainGrid(m6, o["gpts"]), "dim": 3, "vb": 3, "grid_signs": [-1], "valleys": []}
    assert ad.grid_pair_values(table, 3, o["k"], o["uq"], o["nq2"], C_SI, np.zeros(3)) is None
    assert ad._PAIR["grid_off"] and "自检不过" in ad._PAIR["why_grid"]


def test_grid_without_gdiff_shift_fails_self_check(m6, pair_env):
    """0.4.19 不平移末态：若数据其实需要平移（这里故意让一部分末态在另一个表示里），自检会发现并关掉 —— 不会静默算错。"""
    ctx, o = _scene(m6, np.random.default_rng(6))
    pair_env.setitem(ad._PAIR, "gdiff", False)
    pair_env.setattr(ad, "_CTX", ctx)
    table = {"grid": DW.StrainGrid(m6, o["gpts"]), "dim": 3, "vb": 3, "grid_signs": [-1], "valleys": []}
    assert ad.grid_pair_values(table, 3, o["k"], o["uq"], o["nq2"], C_SI, np.zeros(3)) is None
    assert ad._PAIR["grid_off"]


def test_matrix_apply(m6, pair_env):
    ctx, o = _scene(m6, np.random.default_rng(8))
    pair_env.setattr(ad, "_CTX", ctx)
    table = {"grid": DW.StrainGrid(m6, o["gpts"]), "dim": 3, "vb": 3, "grid_signs": [-1],
             "valleys": [], "centers": np.zeros((0, 3))}
    pair_env.setattr(ad, "pair_factor", lambda *a, **k: None)          # 没有 k·p 谷（核心区不混）
    nq = len(o["uq"])
    f0 = np.full((2, 3, nq), 7.0)
    keep = (np.arange(nq) % 2).astype(float)
    vel = np.array([0.01, 0.0, 0.0])
    Fq, ovl, _el = ad.grid_pair_values(table, 3, o["k"], o["uq"], o["nq2"], C_SI, vel)
    pair_env.setitem(ad._PAIR, "grid_off", False)
    got = ad.matrix_apply(f0 * keep, keep, table, 3, o["k"], o["uq"], o["nq2"], C_SI, vel)
    assert np.allclose(got, np.broadcast_to(Fq / ovl * keep, got.shape))
    # 网格不可用：退回 pair_apply(pair_factor) —— 这里 pair_factor 是 None，原样返回
    pair_env.setitem(ad._PAIR, "grid_off", True)
    assert np.array_equal(ad.matrix_apply(f0, keep, table, 3, o["k"], o["uq"], o["nq2"], C_SI, vel), f0)


def test_matrix_apply_core_blend(m6, pair_env):
    """核心区（rcore 以内）子空间里的末态：w·F_kp/重叠 + (1 − w)·F_grid/重叠。"""
    ctx, o = _scene(m6, np.random.default_rng(9))
    pair_env.setattr(ad, "_CTX", ctx)
    table = {"grid": DW.StrainGrid(m6, o["gpts"]), "dim": 3, "vb": 3, "grid_signs": [-1], "valleys": []}
    vel = np.zeros(3)
    Fq, ovl, _el = ad.grid_pair_values(table, 3, o["k"], o["uq"], o["nq2"], C_SI, vel)
    idx = np.array([0, 3, 5])
    w = np.array([1.0, 0.5, 0.0])
    Fk = np.array([2.0, 3.0, 4.0]) * ovl[idx]
    calls = []

    def fake_pf(*a, **k):
        calls.append(k.get("core"))
        return idx, Fk, ovl[idx], w
    pair_env.setattr(ad, "pair_factor", fake_pf)
    got = ad.matrix_apply(np.zeros(len(ovl)), None, table, 3, o["k"], o["uq"], o["nq2"], C_SI, vel)
    exp = Fq / ovl
    exp[idx] = w * np.array([2.0, 3.0, 4.0]) + (1 - w) * exp[idx]
    assert np.allclose(got, exp) and calls == [True]


def test_tables_reference_roundtrip(m6, monkeypatch, tmp_path):
    """共享内存往返（amset.util 的同一实现，multiprocessing.RawArray）：作用不变。"""
    util = tmp_path / "amset"
    util.mkdir()
    (util / "__init__.py").write_text("")
    src = (Path(__file__).resolve().parent / "_amset_util_shared.py")
    (util / "util.py").write_text(src.read_text())
    monkeypatch.syspath_prepend(str(tmp_path))
    for n in [x for x in sys.modules if x == "amset" or x.startswith("amset.")]:
        monkeypatch.delitem(sys.modules, n)
    gpts = np.unique(np.concatenate([E.gvecs(k, E.L0, TW.ENCUT) for k in E.mesh(6, tr_reduce=False)]), axis=0)
    grid = DW.StrainGrid(m6, gpts, offsets=np.linspace(-0.3, 0.3, 6))
    ref = ad.tables_to_reference({"up": {"grid": grid, "dim": 3}, "down": None})
    assert ref["up"]["grid"][0] == "sgrid" and ref["down"] is None
    back = ad.tables_from_reference(ref)["up"]["grid"]
    rng = np.random.default_rng(1)
    v = rng.normal(size=(2, len(gpts))) + 1j * rng.normal(size=(2, len(gpts)))
    k = np.array([0.11, -0.07, 0.2])
    assert np.allclose(back.apply(k, v), grid.apply(k, v)) and np.allclose(back.coverage(k, v), grid.coverage(k, v))


_FAKE_CALC_051 = '''import numpy as np
from amset.electronic_structure.tetrahedron import get_cross_section_values


def get_fine_mesh_qpoints(*a, **k):
    return a[0].qpoints, a[0].weights, a[0].mapping


def pbc_diff(a, b):
    d = np.asarray(a) - np.asarray(b)
    return d - np.round(d)


def calculate_rate(tbs, overlap_calculator, mrta_calculator, elastic_scatterers, inelastic_scatterers,
                   amset_data_min, coeffs, coeffs_mapping, spin, b_idx, k_idx, energy_diff=None):
    velocity = np.zeros(3)
    tet_mask = (tbs.tet_bands, None)
    property_mask, band_kpoint_mask, band_mask, kpoint_mask = tbs.get_masks(spin, tet_mask)
    k = tbs.kpoints[k_idx]
    k_primes = tbs.kpoints[kpoint_mask]
    k_abs_diff = k_primes - k
    k_pbc_diff = pbc_diff(k_primes, k)
    g_diff = (k_abs_diff - k_pbc_diff).astype("int")
    vert_overlap, tet_contributions = tbs.vert_overlap, ()
    tet_overlap = get_cross_section_values(vert_overlap, *tet_contributions)
    qpoints, weights, mapping = get_fine_mesh_qpoints(tbs)
    qpoint_norm_sq = np.sum(qpoints ** 2, axis=-1)
    unit_q = qpoints / np.sqrt(qpoint_norm_sq)[:, None]
    mrta_factor = 1.0
    rates = [s.factor(unit_q, qpoint_norm_sq, spin, b_idx, k, velocity) for s in elastic_scatterers]
    rates = np.array(rates)
    rates *= tet_overlap[mapping] * weights * mrta_factor
    return np.sum(rates, axis=-1)
'''
_FAKE_TET = '''import numpy as np
def get_cross_section_values(values, *a, average=True):
    values = np.asarray(values)
    return values.mean(axis=1) if average else values
class TetrahedralBandStructure:
    def __init__(self, **kw):
        self.__dict__.update(kw)
    def get_masks(self, spin, tetrahedra_mask):
        return self.masks
'''


def test_install_hooks_on_051_shape(tmp_path, monkeypatch):
    pkg = tmp_path / "amset"
    for rel, txt in {"__init__.py": "__version__ = '0.5.1-fake'\n", "scattering/__init__.py": "",
                     "scattering/calculate.py": _FAKE_CALC_051, "electronic_structure/__init__.py": "",
                     "electronic_structure/tetrahedron.py": _FAKE_TET}.items():
        (pkg / rel).parent.mkdir(parents=True, exist_ok=True)
        (pkg / rel).write_text(txt)
    monkeypatch.syspath_prepend(str(tmp_path))
    for n in [x for x in sys.modules if x == "amset" or x.startswith("amset.")]:
        monkeypatch.delitem(sys.modules, n)
    for key in list(ad._PAIR):
        monkeypatch.delitem(ad._PAIR, key)
    monkeypatch.setitem(ad._PAIR, "ok", False)
    monkeypatch.setattr(ad, "_CTX", {})
    ok, why = ad.install_pair_hooks()
    assert ok, why
    assert ad._PAIR["grid_ok"] and ad._PAIR["gdiff"]
    import amset.scattering.calculate as C
    import amset.electronic_structure.tetrahedron as TT
    seen = {}

    class S:
        def factor(self, unit_q, nq2, spin, b, k, v):
            seen.update(ad._CTX)
            return np.ones(len(nq2))
    masks = ("pm", "bkm", "bm", np.array([0]))
    tbs = TT.TetrahedralBandStructure(tet_bands=np.array([0]), kpoints=np.zeros((1, 3)), masks=masks,
                                      vert_overlap=np.ones((1, 4)), qpoints=np.array([[0.1, 0, 0]]),
                                      weights=np.ones(1), mapping=np.array([0]))
    C.calculate_rate(tbs, "OC", None, [S()], [], None, "COEFFS", "CMAP", "up", 2, 0)
    assert seen["rate"] == (tbs, "OC", "COEFFS", "CMAP", "up", 2, 0) and seen["masks"] is masks
    assert seen["tet_contrib"] == ((), {}) and seen["energy_diff"] is None
    assert ad._PAIR["cs"] is not C.get_cross_section_values               # 插值用原函数，不碰钩子
