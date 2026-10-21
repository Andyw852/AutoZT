"""[0061-m] DFT 形变势矩阵（deform_wave.py）：WAVECAR 读写、窗口抽取、投影出应变一阶哈密顿量、对齐、网格上插值算符。

基准是 _epm_si 的经验赝势模型：∂H/∂ε 直接对平面波矩阵求导，任意 k 的矩阵元都有精确值。
  · WAVECAR：写出 / 读回逐位一致（真实 VASP 文件的读法另外对照过 pymatgen 的读取器：标准、自旋极化、非共线）；
  · 投影：X^(P) 与精确矩阵元（同一组本征矢）对上，含 Γ 点三重态内的非对角元（Bir–Pikus 的 b、d）与正确符号；
  · 对齐：带符号形变势里混进的参考能级常数逐图样读回来，离散度 ~0；
  · 网格算符：格点上精确（含时间反演补的点），离格点的误差随网格加密收敛。
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
import _epm_si as E  # noqa: E402
import deform_wave as DW  # noqa: E402

pytest.importorskip("h5py")

ENCUT = 200.0
NB = 12
DELTA = 0.005
REF_SHIFT = np.array([1.3, -0.7, 0.4, 0.0, 0.25, -0.15])          # 每个图样混进去的参考能级常数（eV）


def _s7(tmp, n=4, ewin_k=30.0, ewin_b=1.5):
    """一套 S7：undeformed + 6 个图样 × ±δ，各写 WAVECAR、抽 deform_wave.h5。-> (根目录, 目录清单, 未形变数据)。"""
    kpts = E.mesh(n)
    root = Path(tmp) / "step7_deform"
    dirs, data0 = [], None
    for name, Em in [("undeformed", None)] + [("deform-%02d%s" % (i, "p" if s > 0 else "m"), E.pattern(P, s * DELTA))
                                              for i, P in enumerate(E.PATS) for s in (1, -1)]:
        d = root / name
        d.mkdir(parents=True)
        dat = E.write_wavecar(d / "WAVECAR", kpts, E.lattice(Em), ENCUT, NB)
        assert DW.extract(d / "WAVECAR", d / "deform_wave.h5", ewin_k=ewin_k, ewin_b=ewin_b) is not None
        if name == "undeformed":
            data0 = (kpts, dat)
        else:
            dirs.append(name)
    return root, dirs, data0


def _exact_X(kf, G, V, bands, pats):
    """精确 X^(P)_mn = ⟨m|∂H/∂ε_P|n⟩（同一组本征矢），图样按 pats 的顺序（h5 里是排好序的）。-> (np, nw, nw)。"""
    Vb = V[bands]
    return np.array([Vb.conj() @ E.dH(kf, G, P) @ Vb.T for P in pats])


@pytest.fixture(scope="module")
def s7(tmp_path_factory):
    root, dirs, data0 = _s7(tmp_path_factory.mktemp("dw"))
    out = root.parent / "deformation_matrix.h5"
    summ = DW.build(root, dirs, out)
    return root, dirs, data0, out, summ


def test_wavecar_roundtrip(tmp_path):
    kpts = E.mesh(3)
    dat = E.write_wavecar(tmp_path / "WAVECAR", kpts, E.L0, ENCUT, NB)
    w = DW.Wavecar(tmp_path / "WAVECAR")
    assert (w.nk, w.nb, w.nspin, w.ncl) == (len(kpts), NB, 1, False)
    for k, (G, e, V) in enumerate(dat):
        assert np.array_equal(w.gvectors(k), G)
        assert np.allclose(w.eig[0, k], e[:NB])
        assert np.allclose(w.coeffs(0, k, range(NB)), V[:NB].astype(np.complex64))
    vbm, cbm, _o = DW._band_edges(w.eig, w.occ)
    assert vbm == pytest.approx(max(d[1][3] for d in dat)) and cbm == pytest.approx(min(d[1][4] for d in dat))


def test_strain_pattern_and_extract(s7):
    root, dirs, _d0, _out, summ = s7
    L0 = E.lattice()
    for i, P in enumerate(E.PATS):
        for s in (1, -1):
            sp = DW.strain_pattern(DW.green_lagrange(L0, E.lattice(E.pattern(P, s * DELTA))))
            assert sp[0] == P and sp[1] == pytest.approx(s * DELTA, rel=1e-9)
    assert sorted(summ["patterns"]) == sorted(DW.PAT_NAMES[P] for P in E.PATS)
    assert summ["spins"]["up"]["q_min"] > 0.99
    assert summ["one_sided"] == []


def test_projection_matches_exact(s7):
    """X 与 ⟨m|∂H/∂ε|n⟩ 对上（相对最大值 < 1%：截断球随应变变的 Pulay 效应 + 中心差分 O(δ²)），
    Γ 三重态里的非对角元也对（含符号）。"""
    _root, _dirs, (kpts, dat), out, _s = s7
    m = DW.load_matrix(out, "up")
    kf = m["kpoints"]
    lut = {tuple(np.round(k * 4).astype(int) % 4): j for j, k in enumerate(kpts)}
    worst, gam = 0.0, None
    for jj, k in enumerate(kf):
        j = lut[tuple(np.round(k * 4).astype(int) % 4)]
        G, _e, V = dat[j]
        Xe = _exact_X(k, G, V, m["bands"], m["pairs"])
        Xp = m["X"][:, jj]
        worst = max(worst, float(np.max(np.abs(Xp - Xe)) / np.max(np.abs(Xe))))
        if np.allclose(k, 0):
            gam = (Xp, Xe, m["eig"][jj])
    assert worst < 1e-2
    Xp, Xe, _e = gam
    vb = np.where(np.isin(m["bands"], [E.NOCC - 3, E.NOCC - 2, E.NOCC - 1]))[0]           # Γ 点价带顶三重态
    assert len(vb) == 3
    for ip in range(6):
        assert np.allclose(Xp[ip][np.ix_(vb, vb)], Xe[ip][np.ix_(vb, vb)], atol=1e-2 * np.max(np.abs(Xe)))
    ixy = m["pairs"].index((0, 1))
    off = Xe[ixy][np.ix_(vb, vb)] - np.diag(np.diag(Xe[ixy][np.ix_(vb, vb)]))
    assert np.max(np.abs(off)) > 0.3                                                    # 剪切把三重态耦合起来


def test_alignment_recovers_reference_shift(s7):
    root, _dirs, (kpts, dat), out, _s = s7
    m = DW.load_matrix(out, "up")
    # 带符号文件：S = −dE/dε（排序后本征值中心差分）+ 参考能级常数
    S = np.zeros((NB, len(kpts), 3, 3))
    shift = {P: REF_SHIFT[i] for i, P in enumerate(E.PATS)}
    for ip, P in enumerate(E.PATS):
        ep = np.array([E.solve(k, E.lattice(E.pattern(P, DELTA)), ENCUT)[1][:NB] for k in kpts])
        em = np.array([E.solve(k, E.lattice(E.pattern(P, -DELTA)), ENCUT)[1][:NB] for k in kpts])
        S[:, :, P[0], P[1]] = S[:, :, P[1], P[0]] = (-(ep - em) / (2 * DELTA) + REF_SHIFT[ip]).T
    c, spread, n = DW.alignment(m, S, kpts, band_offset=0, delta=DELTA)
    assert n > 50
    assert np.allclose(c, [shift[P] for P in m["pairs"]], atol=2e-2)
    assert np.max(spread) < 2e-2
    # 带号错位 1：离散度立刻变大（运行期据此拒绝）
    _c2, spread2, _n2 = DW.alignment(m, S, kpts, band_offset=1, delta=DELTA)
    assert spread2 is None or np.max(spread2) > 0.1


def test_strain_grid_on_grid_exact_incl_time_reversal(s7):
    """格点上：⟨m,k|D̂(k)|n,k⟩ = X(k)；时间反演补出的 −k 点：= X(k)*（态取 Θ|n,k⟩）。"""
    _root, _dirs, (kpts, dat), out, _s = s7
    m = DW.load_matrix(out, "up")
    km = E.mesh(4, tr_reduce=False)
    gpts = np.unique(np.concatenate([E.gvecs(k, E.L0, ENCUT) for k in np.concatenate([km, -km])]), axis=0)
    grid = DW.StrainGrid(m, gpts)
    gl = {tuple(g): i for i, g in enumerate(gpts)}

    def dense(G, v):
        out_ = np.zeros(len(gpts), complex)
        out_[[gl[tuple(g)] for g in G]] = v
        return out_

    for jj in range(len(m["kpoints"])):
        k = m["kpoints"][jj]
        G, _e, V = E.solve(k, E.L0, ENCUT)
        vec = np.array([dense(G, V[b]) for b in m["bands"]])
        Dv = grid.apply(k, vec)
        Xg = np.einsum("an,pbn->pab", vec.conj(), Dv)
        assert np.allclose(Xg, m["X"][:, jj], atol=1e-4 * np.max(np.abs(m["X"])))
        # −k：态取时间反演（系数 c(−G)*），X -> X*
        if np.allclose(k, -k - np.round(-k)) is False:
            vtr = np.array([dense(-G, np.conj(V[b])) for b in m["bands"]])
            Xt = np.einsum("an,pbn->pab", vtr.conj(), grid.apply(-k, vtr))
            assert np.allclose(Xt, np.conj(m["X"][:, jj]), atol=1e-4 * np.max(np.abs(m["X"])))


def test_strain_grid_off_grid_converges(tmp_path):
    """离格点：网格算符在价带上的矩阵元 vs 精确 ⟨m,k|∂H/∂ε|n,k⟩，网格从 3³ 加密到 6³ 误差至少降到 1/3（插值误差 O(Δk²)；
    Si 实际的 S7 网格 Δk ≈ 0.06 Å⁻¹，比这里的 6³（≈ 0.33 Å⁻¹）细五倍多）。
    窗口只由未形变构型定、不连锁缩（EPM 的导带到处有交叉，应变把简并劈得各构型不一样）：6³ 时窗口要盖住价带 + 最低两条导带，
    覆盖度 ~1，且不打"缩窗口"的告警。"""
    errs = []
    rng = np.random.default_rng(3)
    ks = (rng.random((4, 3)) - 0.5) * 0.5
    for n in (3, 6):
        root, dirs, _d0 = _s7(tmp_path / ("n%d" % n), n=n, ewin_b=0.8)
        out = root.parent / "m.h5"
        logs = []
        summ = DW.build(root, dirs, out, log=logs.append)
        assert not [x for x in logs if "缩到" in x]
        assert summ["spins"]["up"]["q_min"] > 0.99
        if n == 6:
            lo, hi = summ["spins"]["up"]["bands"]
            assert lo == 0 and hi >= E.NOCC + 1
        m = DW.load_matrix(out, "up")
        gpts = np.unique(np.concatenate([E.gvecs(k, E.L0, ENCUT + 40) for k in ks]), axis=0)
        grid = DW.StrainGrid(m, gpts)
        gl = {tuple(g): i for i, g in enumerate(gpts)}
        e_rel = []
        for k in ks:
            G, e, V = E.solve(k, E.L0, ENCUT)
            vb = [1, 2, 3]
            vec = np.zeros((len(vb), len(gpts)), complex)
            for i, b in enumerate(vb):
                vec[i, [gl[tuple(g)] for g in G]] = V[b]
            Dv = grid.apply(k, vec)
            Xg = np.einsum("an,pbn->pab", vec.conj(), Dv)
            Xe = _exact_X(k, G, V, vb, m["pairs"])
            e_rel.append(np.max(np.abs(Xg - Xe)) / np.max(np.abs(Xe)))
        errs.append(max(e_rel))
    assert errs[1] < errs[0] / 3
    assert errs[1] < 0.05


def test_cli_extract_and_delete(tmp_path, capsys):
    d = tmp_path / "x"
    d.mkdir()
    E.write_wavecar(d / "WAVECAR", E.mesh(2), E.L0, ENCUT, NB)
    assert DW._cli(["extract", "--wavecar", str(d / "WAVECAR"), "--out", str(d / "deform_wave.h5"), "--delete"]) == 0
    assert (d / "deform_wave.h5").is_file() and not (d / "WAVECAR").exists()
    assert DW._cli(["extract", "--wavecar", str(d / "WAVECAR"), "--out", str(d / "o.h5")]) == 0
    assert "没有" in capsys.readouterr().out

def _drop_k_rows(path, drop=(0,)):
    """去掉 deform_wave.h5 里指定 k 行（模拟各构型 IBZ k 点数不同：形变降对称性后 k 集变化）。"""
    import h5py
    with h5py.File(str(path), "r") as f:
        attrs = dict(f.attrs)
        data = {k: np.array(f[k]) for k in f.keys()}
    nk = data["kpoints"].shape[0]
    keep = np.array([j for j in range(nk) if j not in drop], int)
    for k in list(data):
        if k == "kpoints" or k.startswith(("coeff_", "eig_")):
            data[k] = data[k][keep]
    with h5py.File(str(path), "w") as f:
        f.attrs.update(attrs)
        for k, v in data.items():
            f[k] = v


def _worst_vs_exact(m, dat, kpts):
    lut = {tuple(np.round(k * 4).astype(int) % 4): j for j, k in enumerate(kpts)}
    worst = 0.0
    for jj, k in enumerate(m["kpoints"]):
        G, _e, V = dat[lut[tuple(np.round(k * 4).astype(int) % 4)]]
        for ip, P in enumerate(m["pairs"]):
            Xe = _exact_X(k, G, V, m["bands"], [P])[0]
            worst = max(worst, float(np.max(np.abs(m["X"][ip, jj] - Xe)) / np.max(np.abs(Xe))))
    return worst


def test_build_tolerates_different_ibz_k_counts(tmp_path):
    """回归（[0061-m] 真数据 A1 崩溃）：各构型 IBZ k 点数不同（形变降对称性后 k 集变化）时，build 必须按 k
    坐标对齐；旧代码拿"绝对 k 号"去索引 jk 过滤后的未形变数组（E0all[j]），在 Si 上表现为
    IndexError: index 1224 is out of bounds for axis 0 with size 1223。"""
    import shutil

    root, dirs, data0 = _s7(tmp_path, n=4)
    kpts, dat = data0
    p0 = [d for d in dirs if d.startswith("deform-00")]
    assert len(p0) == 2

    # A) 中心差分（± 都在）：只让 pattern-00 的两个构型各少第 0 个 k -> 公共 k 减一
    rootA = tmp_path / "nk_central"
    shutil.copytree(str(root), str(rootA))
    for d in p0:
        _drop_k_rows(rootA / d / "deform_wave.h5", drop=(0,))
    outA = tmp_path / "dm_central.h5"
    summA = DW.build(rootA, dirs, outA, log=lambda *a: None)
    assert summA["n_k"] == len(kpts) - 1
    assert summA["one_sided"] == []
    assert _worst_vs_exact(DW.load_matrix(outA, "up"), dat, kpts) < 1e-2

    # B) 单侧差分 + k 数不同（真数据 Si 的现场）：pattern-00 只剩一侧，且它少了第 0 个 k
    drop_dir = "deform-00p"
    gone = [d for d in p0 if d != drop_dir][0]
    rootB = tmp_path / "nk_onesided"
    shutil.copytree(str(root), str(rootB))
    shutil.rmtree(str(rootB / drop_dir))
    _drop_k_rows(rootB / gone / "deform_wave.h5", drop=(0,))
    dirsB = [d for d in dirs if d != drop_dir]
    outB = tmp_path / "dm_onesided.h5"
    summB = DW.build(rootB, dirsB, outB, log=lambda *a: None)
    assert summB["n_k"] == len(kpts) - 1
    assert summB["one_sided"] == [DW.PAT_NAMES[E.PATS[0]]]
    mB = DW.load_matrix(outB, "up")
    assert np.all(np.isfinite(mB["X"]))
    assert _worst_vs_exact(mB, dat, kpts) < 0.2      # 单侧差分 O(δ) + 截断球 Pulay，比中心差分松一档
