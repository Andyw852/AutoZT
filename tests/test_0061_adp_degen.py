"""[0061] 简并 / 近简并带边的 ADP 带间修正（通用版）：adp_degen.py + amset_dpo3d / amset2d_plugin 接线 + S8 / S8.4 gen。

物理核对（纯 numpy，不需要 amset）：合成的"DFT"——已知子空间 k·p 模型 + 应变算符 -> 全网格能量与（每个网格点随机相位、
简并子空间里随机组合的）波函数系数、AMSET 密网格能带、S7.1 的带符号 / 幅值形变势（±δ 应变、按能量排序、含 k 随应变
移动那一项）。setup 只看这些数据，读出对称性、k·p、应变算符，比值表对上真模型直接算的比值表：
  · 立方 Γ 三重态（Si，无 SOC）：48×2 个对称操作、k·p 3 个参数（L/M/N）、应变 3 个（l/m/n），与旧版立方实现逐位一致；
  · Si + SOC（Γ8 + Γ7，6 条带，Kramers 简并）；
  · 二维六方 Γ 点双重态、K/K′ 谷星（K′ 由时间反演映回 K）；
  · 正交晶系近简并（劈裂 20 meV）的两条带；
  · 二维 p_z + p_x/p_y（面内 k 不耦合 -> 自动分块，只修正双重态）；
  · Kramers 对（处处简并）与不简并的谷：不修正；
  · S7 网格太粗、符号定不住：几个解的放大系数相近取平均，差得多不修正；
  · 远处有一条耦合的带（4×4 k·p）：三带子空间拟合只是近似，放大系数仍在 4% 以内。
"""
import importlib.util
import io
import itertools
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
S84 = KE / "step8.4_amset2d"
for _p in (str(KE), str(S84), str(ROOT / "skill" / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import adp_degen as ad  # noqa: E402

H2M = 3.80998212
A_SI = 5.431
SI_LAT = 0.5 * A_SI * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
SI_B = 2 * np.pi * np.linalg.inv(SI_LAT).T
KP_SI = (-5.568, -3.596, -8.457)                 # 本地 QE PBE（无 SOC）拟合
BP_SI = (1.74 + 2 * -2.31, 1.74 + 2.31, -9.46)   # l = a + 2b、m = a − b、n（离子弛豫）


# ---------------------------------------------------------------------------
# 合成数据
# ---------------------------------------------------------------------------
def _mesh(n, dim=3):
    ax = [np.arange(n) / n for _ in range(dim)] + [np.zeros(1)] * (3 - dim)
    kf = np.array(list(itertools.product(*ax)), float)
    return kf - np.round(kf)


def _oh():
    out = []
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1, -1), repeat=3):
            M = np.zeros((3, 3))
            for i in range(3):
                M[i, p[i]] = s[i]
            out.append(M)
    return out


def _d6h():
    out = []
    mx = np.diag([1.0, -1, 1])
    for i in range(6):
        c, s = math.cos(i * math.pi / 3), math.sin(i * math.pi / 3)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
        for X in (R, R @ mx):
            for z in (1.0, -1.0):
                out.append(X @ np.diag([1, 1, z]))
    return out


def _cubic_C(C11=165.7, C12=63.9, C44=79.6):
    C = np.zeros((3, 3, 3, 3))
    for i, j, k, l in itertools.product(range(3), repeat=4):
        if i == j and k == l:
            C[i, j, k, l] = C11 if i == k else C12
        elif ((i == k and j == l) or (i == l and j == k)) and i != j:
            C[i, j, k, l] = C44
    return C


def _hex_C(C11=130.0, C12=30.0):
    C = np.zeros((3, 3, 3, 3))
    for i, j, k, l in itertools.product(range(2), repeat=4):
        if i == j and k == l:
            C[i, j, k, l] = C11 if i == k else C12
        elif ((i == k and j == l) or (i == l and j == k)) and i != j:
            C[i, j, k, l] = (C11 - C12) / 2
    C[2, 2, 2, 2] = 1.0
    return C


def dkk(L, M, N, h2m=H2M):
    def H(k):
        k = np.asarray(k, float)
        out = np.zeros((3, 3))
        for i in range(3):
            out[i, i] = L * k[i] ** 2 + M * (k @ k - k[i] ** 2)
            for j in range(3):
                if i != j:
                    out[i, j] = N * k[i] * k[j]
        return h2m * out
    return H


def bir_pikus(l, m, n):
    Y = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for p in range(3):
            Y[i, i, p, p] = l if p == i else m
        for j in range(3):
            if i != j:
                Y[i, j, i, j] = Y[i, j, j, i] = n / 2
    return Y


def e_doublet(al, be, h2m=H2M):
    def H(k):
        kx, ky = k[0], k[1]
        return h2m * (al * (kx * kx + ky * ky) * np.eye(2)
                      + be * np.array([[kx * kx - ky * ky, 2 * kx * ky], [2 * kx * ky, ky * ky - kx * kx]]))
    return H


def y_doublet(a, b, c=0.0):
    Y = np.zeros((3, 3, 2, 2))
    sz, sx = np.diag([1.0, -1.0]), np.array([[0, 1.0], [1.0, 0]])
    Y[0, 0], Y[1, 1] = a * np.eye(2) + b * sz, a * np.eye(2) - b * sz
    Y[0, 1] = Y[1, 0] = b * sx
    Y[2, 2] = c * np.eye(2)
    return Y


def hex_B(a=3.2, c=20.0):
    return 2 * np.pi * np.linalg.inv(np.array([[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, c]])).T


def true_model(Hfun, dim, d, rmax):
    rng = np.random.default_rng(3)
    K = rng.normal(size=(60, 3)) * 0.3
    if dim == 2:
        K[:, 2] = 0
    D0 = Hfun(np.zeros(3))
    F = ad.monomials(K, dim)
    Cf = np.linalg.lstsq(F, np.array([Hfun(k) - D0 for k in K]).reshape(len(K), -1), rcond=None)[0]
    return {"dim": dim, "D0": D0, "C": Cf.reshape(-1, d, d), "rmax": rmax}


def _degen_groups(e, tol=1e-9):
    out = [[0]]
    for n in range(1, len(e)):
        if e[n] - e[n - 1] < tol:
            out[-1].append(n)
        else:
            out.append([n])
    return out


class Synth:
    """子空间模型 Hm(δk)（d×d，电子能量，带边 0）+ 应变算符 Ym (3,3,d,d)（张量约定），谷星由 ops（含时间反演）生成；
    另加一条深带、一条远带。能量、系数、AMSET 密网格、S7.1 文件都从它算。"""

    def __init__(self, B, ops, Hm, Ym, d, sign=-1, n=12, m=8, nd=24, dim=3, k0f=(0, 0, 0), seed=0):
        self.B, self.ops, self.Hm, self.Ym, self.d, self.sign, self.dim = B, ops, Hm, Ym, d, sign, dim
        self.rng = np.random.default_rng(seed)
        rinv = np.linalg.inv(B)
        mem = []
        for R in ops:
            for s in (1, -1):
                f = ad._wrap((s * R @ (np.asarray(k0f, float) @ B)) @ rinv)
                if not any(np.max(np.abs(ad._wrap(f - x[0]))) < 1e-6 for x in mem):
                    mem.append((f, s * R))
        self.members = mem
        self.nb, self.ncoeff = d + 2, d + 4
        Q, _ = np.linalg.qr(self.rng.normal(size=(self.ncoeff,) * 2) + 1j * self.rng.normal(size=(self.ncoeff,) * 2))
        self.Q = Q
        self.kf = _mesh(n, dim)
        self.E = self._bands(self.kf, coeffs=True)
        self.dense_kf = _mesh(nd, dim)
        self.dense_E = self._bands(self.dense_kf)
        self.kS = _mesh(m, dim)
        self._s7(self.kS)

    def _local(self, kf):
        best = None
        for f, M in self.members:
            dk = ad._wrap(kf - f) @ self.B
            if self.dim == 2:
                dk[2] = 0
            r = np.linalg.norm(dk)
            if best is None or r < best[0]:
                best = (r, M.T @ dk, M)
        return best[1], best[2]

    def _mH(self, kf, eps=None):
        dkr, M = self._local(kf)
        if eps is None:
            return self.Hm(dkr)
        er = M.T @ eps @ M
        return self.Hm(dkr - er @ dkr) + np.einsum("ij,ijab->ab", er, self.Ym)

    def _others(self, kf, eps=None):
        dk = ad._wrap(kf) @ self.B
        if self.dim == 2:
            dk[2] = 0
        tr = 0.0 if eps is None else np.trace(eps)
        if eps is not None:
            dk = dk - eps @ dk
        if self.sign < 0:
            return -30.0 - 50.0 * H2M * dk @ dk + tr, 1.0 + 0.5 * H2M * dk @ dk - 3.0 * tr
        return -1.0 - 0.5 * H2M * dk @ dk + tr, 30.0 + 50.0 * H2M * dk @ dk - 3.0 * tr

    def _bands(self, kfs, coeffs=False):
        E = np.zeros((self.nb, len(kfs)))
        if coeffs:
            self.cdict = {}
        for i, kf in enumerate(kfs):
            e, V = np.linalg.eigh(self._mH(kf))
            lo, hi = self._others(kf)
            E[:, i] = np.concatenate([[lo], e, [hi]])
            if coeffs:
                Vg = V * np.exp(1j * self.rng.uniform(0, 2 * np.pi, self.d))
                for G in _degen_groups(e):
                    if len(G) > 1:
                        U, _ = np.linalg.qr(self.rng.normal(size=(len(G),) * 2) + 1j * self.rng.normal(size=(len(G),) * 2))
                        Vg[:, G] = Vg[:, G] @ U
                c = np.zeros((self.nb, self.ncoeff), complex)
                c[0], c[-1] = self.Q[:, self.d], self.Q[:, self.d + 1]
                c[1:self.d + 1] = (self.Q[:, :self.d] @ Vg).T
                self.cdict[ad._key(kf)] = c
        return E

    def coeff_fn(self, bands, kfs, kref):
        return np.array([self.cdict[ad._key(kf)][b] for b, kf in zip(bands, np.asarray(kfs))])

    def _s7(self, kS, delta=0.005):
        S, A, gap = (np.zeros((self.nb, len(kS), 3, 3)), np.zeros((self.nb, len(kS), 3, 3)),
                     np.zeros((self.nb, len(kS))))
        for j, kf in enumerate(kS):
            lo, hi = self._others(kf)
            e0 = np.sort(np.concatenate([[lo], np.linalg.eigvalsh(self._mH(kf)), [hi]]))
            gg = np.diff(e0)
            gap[:, j] = np.minimum(np.concatenate([[np.inf], gg]), np.concatenate([gg, [np.inf]]))
            for a, b in ad.pairs(3):
                Pm = np.zeros((3, 3))
                Pm[a, b] = Pm[b, a] = 1.0
                v = []
                for s in (1, -1):
                    lo2, hi2 = self._others(kf, s * delta * Pm)
                    ed = np.sort(np.concatenate([[lo2], np.linalg.eigvalsh(self._mH(kf, s * delta * Pm)), [hi2]]))
                    v.append((e0 - ed) / (s * delta))
                S[:, j, a, b] = S[:, j, b, a] = 0.5 * (v[0] + v[1])
                A[:, j, a, b] = A[:, j, b, a] = 0.5 * (np.abs(v[0]) + np.abs(v[1]))
        self.sign_h5 = {"kpoints": kS, "S": S, "delta": delta, "gap": gap}
        self.amp_h5 = {"kpoints": kS, "A": A / (2.0 - np.eye(3)), "offdiag": "tensor"}   # V189 张量约定（非对角 ÷2）

    def run(self, C, nang=24, efit=0.15, amp=True, dense=True):
        vb = self.d if self.sign < 0 else 0
        return ad.setup(self.kf, self.E, self.coeff_fn, vb, self.B, self.ops, C, self.sign_h5,
                        self.amp_h5 if amp else None, dense=(self.dense_kf, self.dense_E) if dense else None,
                        dim=self.dim, efit=efit, nang=nang)


def _check(tab, Hm, Ym, d, dim, C, nang, sign=-1, tol=5e-3):
    v = tab["valleys"][0]
    t = next(x for x in v["tables"] if x is not None)
    tt = ad.rate_ratio_table(true_model(Hm, dim, d, v["rcut"]), Ym, C, sign, t["E"], nang)
    assert np.max(np.abs(t["ratio"] / tt["ratio"] - 1)) < tol
    return t


@pytest.fixture(scope="module")
def si():
    return Synth(SI_B, _oh(), dkk(*KP_SI), bir_pikus(*BP_SI), 3, n=12, m=8, nd=24)


# ---------------------------------------------------------------------------
# 约定
# ---------------------------------------------------------------------------
def test_strain_conventions_and_forward_k0():
    """张量约定 Y_ab、原生图样 X^(ab) = 2Y_ab（剪切 ε_ab = ε_ba = δ）；k0 点严格重演 S7.1：
    立方 Γ 三重态 ε_xx 劈裂 {l, m, m} -> S = −(l+m)/2, −m, −(l+m)/2；ε_xy 劈裂 {±n, 0} -> S = 0、A = |n|, 0, |n|。"""
    l, m, n = BP_SI
    Z = ad.Z_from_Y(bir_pikus(l, m, n), 3)
    assert np.allclose(ad.tensor_Y(Z, 3), bir_pikus(l, m, n))
    eps = np.random.default_rng(0).normal(size=(3, 3))
    eps = eps + eps.T
    assert np.allclose(ad.strain_H(Z, eps, 3), np.einsum("ij,ijab->ab", eps, bir_pikus(l, m, n)))
    Sxx, _Axx = ad.forward_k0(np.zeros((3, 3)), ad.pattern_op(Z, 3, (0, 0)), 0.005)
    assert np.allclose(Sxx, [-(l + m) / 2, -m, -(l + m) / 2])
    Sxy, Axy = ad.forward_k0(np.zeros((3, 3)), ad.pattern_op(Z, 3, (0, 1)), 0.005)
    assert np.allclose(Sxy, 0, atol=1e-9) and np.allclose(Axy, [abs(n), 0, abs(n)])


def test_ratio_table_matches_cubic_reference():
    """只有静水形变势：完整矩阵元 = AMSET 写法（比值 ≡ 1）；Si：与旧版立方实现（QE 核对过约定）同一数值、与能量无关
    （amset_terms=False = 首版口径：只用 H_ε、总散射率）。"""
    C = _cubic_C()
    m = true_model(dkk(*KP_SI), 3, 3, 0.5)
    t = ad.rate_ratio_table(m, bir_pikus(1.74, 1.74, 0.0), C, -1, [0.02, 0.1], 48, amset_terms=False)
    assert np.allclose(t["ratio"], 1.0, atol=1e-9)
    t = ad.rate_ratio_table(m, bir_pikus(*BP_SI), C, -1, [0.01, 0.05, 0.15], 48, amset_terms=False)
    mean = t["ratio"].mean(axis=1)
    assert np.allclose(mean, mean[0], rtol=1e-9)                     # 纯二次型：比值与能量无关
    assert np.allclose(mean[0], [1.93037739, 3.20945402, 2.26361053], rtol=1e-6)   # 旧版 rate_ratio（nang=48）


def test_ratio_table_amset_terms():
    """出厂口径照 AMSET 弹性散射的实际算法（动量弛豫因子、|固定分数坐标形变势| + m_e v⊗v）：只有静水时带边附近
    仍 ≈ 1（剩下的是 AMSET 对"形变势 − 运动学项"逐分量取 |·| 的差）；Si 的放大系数与首版差几个百分点、随能量略降。"""
    C = _cubic_C()
    m = true_model(dkk(*KP_SI), 3, 3, 0.5)
    t = ad.rate_ratio_table(m, bir_pikus(1.74, 1.74, 0.0), C, -1, [0.01], 48)
    assert np.allclose(t["ratio"].mean(axis=1), 1.0, atol=0.01)
    t = ad.rate_ratio_table(m, bir_pikus(*BP_SI), C, -1, [0.01, 0.1], 48)
    mean = t["ratio"].mean(axis=1)
    assert np.all(mean[1] < mean[0]) and np.allclose(mean[0], [1.93037739, 3.20945402, 2.26361053], rtol=0.08)


def test_symmetry_helpers():
    """表示矩阵唯一到相位；公共不变子空间；群闭合（含反幺正操作）。"""
    rng = np.random.default_rng(5)
    U0, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    Hk = dkk(*KP_SI)
    R = _oh()[7]
    prs = [(U0.conj().T @ Hk(k) @ U0, U0.conj().T @ Hk(R @ k) @ U0) for k in rng.normal(size=(6, 3)) * 0.1]
    U, err = ad.solve_rep(prs, 3, False)
    assert U is not None and err < 1e-9
    assert all(np.allclose(H2, U @ H @ U.conj().T) for H, H2 in prs)
    _W, blocks = ad.commutant_blocks([np.diag([1.0, 2.0]), np.diag([0.5, 3.0])], 2)
    assert sorted(len(b) for b in blocks) == [1, 1]
    _W, blocks = ad.commutant_blocks([np.diag([1.0, 2.0]), np.array([[0.5, 0.3], [0.3, 3.0]])], 2)
    assert blocks == [[0, 1]]
    G = ad.close_group([(np.diag([-1.0, -1, -1]), np.eye(2), True)], 2)
    assert len(G) == 2 and sum(g[2] for g in G) == 1


# ---------------------------------------------------------------------------
# setup：合成数据 -> 修正表
# ---------------------------------------------------------------------------
def test_setup_cubic_gamma_triplet(si):
    tab, rep = si.run(_cubic_C(), nang=48)
    v = rep["valleys"][0]
    assert rep["status"] == "ok" and v["status"] == "ok" and v["d"] == 3 and v["blocks"] == [3]
    assert v["symmetry"][0]["n_ops"] == 96 and v["symmetry"][0]["n_kp_params"] == 3
    assert v["symmetry"][0]["n_strain_params"] == 3                  # l、m、n
    assert v["strain"]["k0_rel_res"] < 1e-3 and v["kp"]["rel_res"] < 1e-6
    _check(tab, dkk(*KP_SI), bir_pikus(*BP_SI), 3, 3, _cubic_C(), 48, tol=1e-3)
    assert any(x.get("carrier") == "electron" and "不简并" in x["status"] for x in rep["valleys"])
    # 运行期查表：子空间里的带 > 1；子空间外的带、离 Γ 太远 -> 1；Γ 点取方向平均
    near = np.array([0.02, 0.01, 0.0])
    assert all(ad.scale(tab, b, near) > 1.2 for b in (1, 2, 3))
    assert ad.scale(tab, 4, near) == 1.0 and ad.scale(tab, 0, near) == 1.0
    assert ad.scale(tab, 3, np.array([0.3, 0.3, 0.1])) == 1.0
    assert ad.scale(tab, 3, np.zeros(3)) > 1.2
    # 修正在 rcut 前线性过渡到 1（连续）
    rc = tab["valleys"][0]["rcut"]
    kd = np.array([1.0, 0.3, 0.2]) / np.linalg.norm([1.0, 0.3, 0.2])
    f = [ad.scale(tab, 3, (x * rc * kd) @ np.linalg.inv(SI_B)) for x in (0.95, 0.999)]
    assert f[0] > 1.0 and abs(f[1] - 1.0) < 0.02


def test_setup_without_amplitude_file(si):
    tab, rep = si.run(_cubic_C(), nang=24, amp=False)
    assert rep["status"] == "ok" and rep["valleys"][0]["strain"]["k0_amplitude"] is False
    _check(tab, dkk(*KP_SI), bir_pikus(*BP_SI), 3, 3, _cubic_C(), 24, tol=1e-3)


def test_setup_soc_gamma8_gamma7():
    """Si + SOC（Δ = 44 meV）：Γ8 + Γ7 六条带一起进子空间（Kramers 对处处简并），旋量表示由数据解出。"""
    Lm = [np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]), np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]]),
          np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])]
    sg = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])]
    Hso = 0.044 / 3 * (sum(np.kron(L, s) for L, s in zip(Lm, sg)) - np.eye(6))
    Hk = dkk(*KP_SI)

    def Hm(k):
        return np.kron(Hk(k), np.eye(2)) + Hso
    Ym = np.einsum("ijab,cd->ijacbd", bir_pikus(*BP_SI), np.eye(2)).reshape(3, 3, 6, 6)
    sy = Synth(SI_B, _oh(), Hm, Ym, 6, n=12, m=8, nd=16)
    tab, rep = sy.run(_cubic_C(), nang=24)
    v = rep["valleys"][0]
    assert v["status"] == "ok" and v["d"] == 6 and v["blocks"] == [6] and v["symmetry"][0]["n_ops"] == 96
    t = _check(tab, Hm, Ym, 6, 3, _cubic_C(), 24, tol=2e-2)
    m = t["ratio"].mean(axis=(0, 1))
    assert np.allclose(m[0::2], m[1::2], rtol=1e-6)                   # Kramers 对的两个态放大系数相同


@pytest.mark.parametrize("k0f, n_members", [((0, 0, 0), 1), ((1 / 3, 1 / 3, 0), 2)])
def test_setup_2d_hexagonal_doublet(k0f, n_members):
    """二维（面内 k、面内两支声学模）：Γ 点 p_x/p_y 型双重态；K 谷（K′ 由时间反演映回，两边放大系数相同）。"""
    Hm, Ym = e_doublet(-2.0, 0.8), y_doublet(-3.0, 2.5, 0.5)
    sy = Synth(hex_B(), _d6h(), Hm, Ym, 2, n=18, m=12, nd=36, dim=2, k0f=k0f)
    tab, rep = sy.run(_hex_C(), nang=48)
    v = rep["valleys"][0]
    assert v["status"] == "ok" and v["n_members"] == n_members and rep["dim"] == 2
    _check(tab, Hm, Ym, 2, 2, _hex_C(), 48, tol=5e-3)
    if n_members == 2:
        k1, k2 = np.array([1 / 3 + 0.01, 1 / 3, 0]), np.array([-1 / 3 - 0.01, -1 / 3, 0])
        assert ad.scale(tab, 1, k1) == pytest.approx(ad.scale(tab, 1, k2), rel=1e-9) and ad.scale(tab, 1, k1) > 1.1


def test_setup_2d_block_decomposition():
    """二维 p_z（低 20 meV）+ p_x/p_y：面内 k、面内应变都不耦合 -> 自动分成 [1] + [2] 两块，只修正双重态。"""
    HE = e_doublet(-2.0, 0.8)

    def Hm(k):
        out = np.zeros((3, 3))
        out[:2, :2] = HE(k)
        out[2, 2] = H2M * -1.5 * (k[0] ** 2 + k[1] ** 2) - 0.02
        return out
    Ym = np.zeros((3, 3, 3, 3))
    Ym[:, :, :2, :2] = y_doublet(-3.0, 2.5, 0.5)
    Ym[0, 0, 2, 2] = Ym[1, 1, 2, 2] = -1.0
    sy = Synth(hex_B(), _d6h(), Hm, Ym, 3, n=18, m=12, nd=36, dim=2)
    tab, rep = sy.run(_hex_C(), nang=48)
    v, tv = rep["valleys"][0], tab["valleys"][0]
    assert v["status"] == "ok" and sorted(v["blocks"]) == [1, 2]
    bi = next(i for i, b in enumerate(tv["blocks"]) if len(b) == 2)
    t = tv["tables"][bi]
    tt = ad.rate_ratio_table(true_model(HE, 2, 2, tv["rcut"]), y_doublet(-3.0, 2.5, 0.5), _hex_C(), -1, t["E"], 48)
    assert np.max(np.abs(t["ratio"] / tt["ratio"] - 1)) < 5e-3
    k = np.array([0.03, 0.01, 0.0])                                  # 带序：双重态下支、p_z、双重态上支
    f = [ad.scale(tab, b, k) for b in (1, 2, 3)]
    assert f[0] > 1.1 and f[1] == 1.0 and f[2] > 0 and abs(f[2] - 1.0) > 0.05   # 上支也修正（动量弛豫后可 < 1）


def test_setup_near_degenerate_orthorhombic():
    """正交晶系（mmm）两条带劈裂 20 meV：k·p 的 k_x k_y 耦合与 ε_xy 耦合的相对符号由邻近 S7 点定出来。"""
    B = 2 * np.pi * np.linalg.inv(np.diag([4.1, 5.3, 6.2])).T
    ops = [np.diag(s) for s in itertools.product((1.0, -1.0), repeat=3)]

    def Hm(k):
        kx, ky, kz = k
        return H2M * np.array([[-1.2 * kx * kx - 0.5 * ky * ky - 0.8 * kz * kz, -0.9 * kx * ky],
                               [-0.9 * kx * ky, -0.6 * kx * kx - 1.4 * ky * ky - 0.7 * kz * kz]]) + np.diag([0, -0.02])
    Ym = np.zeros((3, 3, 2, 2))
    Ym[0, 0], Ym[1, 1], Ym[2, 2] = np.diag([2.0, 4.0]), np.diag([3.5, 1.5]), np.diag([1.0, 1.2])
    Ym[0, 1] = Ym[1, 0] = 3.0 * np.array([[0, 1.0], [1.0, 0]])
    C = np.zeros((3, 3, 3, 3))
    for (i, j), x in {(0, 0): 150, (1, 1): 120, (2, 2): 100, (0, 1): 40, (0, 2): 35, (1, 2): 30}.items():
        C[i, i, j, j] = C[j, j, i, i] = x
    for (i, j), x in {(0, 1): 50, (0, 2): 40, (1, 2): 45}.items():
        C[i, j, i, j] = C[j, i, j, i] = C[i, j, j, i] = C[j, i, i, j] = x
    sy = Synth(B, ops, Hm, Ym, 2, n=12, m=10, nd=30)
    tab, rep = sy.run(C, nang=48)
    v = rep["valleys"][0]
    assert v["status"] == "ok" and v["symmetry"][0]["n_ops"] == 16 and not v["strain"].get("pick")
    assert v["split_eV"] == pytest.approx([-0.02, 0.0], abs=1e-9)
    _check(tab, Hm, Ym, 2, 3, C, 48, tol=5e-3)


def test_setup_not_needed_cases(si):
    """Kramers 对（处处 ∝ 单位阵）与不简并的谷不修正；带数对不上直接报。"""
    Ym = np.zeros((3, 3, 2, 2))
    for i in range(3):
        Ym[i, i] = 2.0 * np.eye(2)
    sy = Synth(SI_B, _oh(), lambda k: H2M * -1.0 * (k @ k) * np.eye(2), Ym, 2, n=10, m=8, nd=12)
    tab, rep = sy.run(_cubic_C())
    assert tab is None and "一维块" in rep["status"]
    sy = Synth(SI_B, _oh(), lambda k: H2M * -1.0 * (k @ k) * np.eye(1), np.zeros((3, 3, 1, 1)), 1, n=10, m=8, nd=12)
    tab, rep = sy.run(_cubic_C())
    assert tab is None and "不简并" in rep["status"]
    bad = dict(si.sign_h5, S=si.sign_h5["S"][1:])
    tab, rep = ad.setup(si.kf, si.E, si.coeff_fn, 3, SI_B, si.ops, _cubic_C(), bad)
    assert tab is None and "带数" in rep["status"]


def test_setup_ambiguous_sign(monkeypatch):
    """S7 网格只有 2×2×2（Γ 附近没有别的点）：Bir–Pikus n 的符号定不住。两个解的放大系数差 15%（> AMBIG_SPREAD）
    -> 不修正并说明要加密 S7；放宽阈值时取两解的平均并写明。"""
    sy = Synth(SI_B, _oh(), dkk(*KP_SI), bir_pikus(*BP_SI), 3, n=12, m=2, nd=24)
    tab, rep = sy.run(_cubic_C(), nang=24)
    assert tab is None and "分不开 2 个" in rep["valleys"][0]["status"] and "加密 S7" in rep["valleys"][0]["status"]
    monkeypatch.setattr(ad, "AMBIG_SPREAD", 0.3)
    tab, rep = sy.run(_cubic_C(), nang=24)
    st = rep["valleys"][0]["strain"]
    assert rep["status"] == "ok" and st["ambiguous"]["n"] == 2 and "取平均" in st["pick"]
    t = tab["valleys"][0]["tables"][0]
    tt = ad.rate_ratio_table(true_model(dkk(*KP_SI), 3, 3, 0.5), bir_pikus(*BP_SI), _cubic_C(), -1, t["E"], 24)
    assert np.allclose(t["ratio"].mean(axis=(0, 1)), tt["ratio"].mean(axis=(0, 1)), rtol=0.1)


def test_setup_remote_band_robust():
    """真 Hamiltonian 是 4×4（价带 + 3.4 eV 外的 s 带，k·p 耦合 6 eV·Å）：三带子空间拟合是近似（k·p 残差约 10%），
    各带平均放大系数仍在 4% 以内（轻空穴受远带影响最大）。"""
    Hv = dkk(*KP_SI)

    def H4(k):
        out = np.zeros((4, 4))
        out[:3, :3] = Hv(k)
        out[:3, 3] = out[3, :3] = 6.0 * np.asarray(k)
        out[3, 3] = 3.4 + H2M * (k @ k)
        return out
    Y4 = np.zeros((3, 3, 4, 4))
    Y4[:, :, :3, :3] = bir_pikus(*BP_SI)
    for i in range(3):
        Y4[i, i, 3, 3] = -8.0
    sy = Synth(SI_B, _oh(), H4, Y4, 4, n=16, m=8, nd=32)
    tab, rep = ad.setup(sy.kf, sy.E, sy.coeff_fn, 3, SI_B, sy.ops, _cubic_C(), sy.sign_h5, sy.amp_h5,
                        dense=(sy.dense_kf, sy.dense_E), nang=24)
    v = rep["valleys"][0]
    assert v["status"] == "ok" and v["d"] == 3 and 0.01 < v["kp"]["rel_res"] < ad.FIT_REL_TOL
    t = tab["valleys"][0]["tables"][0]
    tt = ad.rate_ratio_table(true_model(H4, 3, 4, tab["valleys"][0]["rcut"]), Y4, _cubic_C(), -1, t["E"], 24)
    # 真模型的载流子序 0 是远带（空穴能量为负，没有态 -> 1），价带从序 1 开始
    assert np.allclose(t["ratio"].mean(axis=(0, 1)), tt["ratio"][..., 1:].mean(axis=(0, 1)), rtol=0.04)


# ---------------------------------------------------------------------------
# 插件接线（不装 AMSET：最小桩模块；vasprun / 点群由测试替换）
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
        "        return np.ones(np.shape(norm_q_sq))\n"
        "class IonizedImpurityScattering:\n"
        "    name = 'IMP'\n"
        "class PiezoelectricScattering:\n"
        "    name = 'PIE'\n"
        "def get_christoffel_tensors(c, q):\n"
        "    return np.einsum('ijkl,ni,nl->njk', c, q, q)\n"),
    "amset/scattering/inelastic.py": "class PolarOpticalScattering:\n    name = 'POP'\n",
    "amset/scattering/calculate.py": (
        "from amset.scattering.elastic import AcousticDeformationPotentialScattering\n"
        "_scattering_mechanisms = {'ADP': AcousticDeformationPotentialScattering}\n"
        "def scattering_worker(*a, **k):\n"
        "    return 'orig'\n"),
    "amset/constants.py": ("angstrom_to_bohr = 1.8897259886\nboltzmann_au = 3.166811563e-06\n"
                           "coulomb_to_au = 6.241509074e18\nm_to_bohr = 1.8897259886e10\n"
                           "hartree_to_ev = 27.211386245988\n"),
    "amset/electronic_structure/__init__.py": "",
    "amset/electronic_structure/fd.py": "def fd(*a, **k):\n    return 0.0\n",
    "amset/core/__init__.py": "",
    "amset/core/transport.py": "",
}


class _Spin:
    name = "up"


class _Lat:
    def __init__(self, B):
        self.reciprocal_lattice = type("R", (), {"matrix": B})()


class _OC:
    """WavefunctionOverlapCalculator 的替身（带号是 AMSET 的，与合成数据的带号相同）。"""

    def __init__(self, sy):
        self.sy = sy
        self.g_maps = np.tile(np.arange(sy.ncoeff), (3, 3, 3, 1))

    def get_coefficients(self, spin, bands, kpoints):
        return self.sy.coeff_fn(np.asarray(bands), kpoints, None)


def _amset_data(sy):
    spin = _Spin()
    return spin, type("AD", (), {
        "spins": [spin], "energies": {spin: sy.dense_E / ad.HARTREE_EV}, "kpoints": sy.dense_kf,
        "vb_idx": {spin: sy.d if sy.sign < 0 else 0}, "is_metal": False,
        "structure": type("St", (), {"lattice": _Lat(sy.B)})(), "overlap_calculator": _OC(sy)})()


def _write_h5(sy, run):
    import h5py
    with h5py.File(str(run / "deformation_signed.h5"), "w") as f:
        f["signed_up"] = sy.sign_h5["S"]
        f["gap_up"] = sy.sign_h5["gap"]
        f["kpoints"] = sy.kS
        f.attrs["strain_amount"] = 0.005
    with h5py.File(str(run / "deformation.h5"), "w") as f:
        f["deformation_potentials_up"] = sy.amp_h5["A"]
        f["kpoints"] = sy.kS
        f.attrs["dp_offdiag"] = "tensor"


def _prep_run(tmp_path, monkeypatch, files, record_name, record, stubs=None):
    run = tmp_path / "run"
    run.mkdir()
    for rel, txt in dict(_STUBS, **(stubs or {})).items():
        (run / rel).parent.mkdir(parents=True, exist_ok=True)
        (run / rel).write_text(txt)
    for f in files:
        shutil.copy(f, run)
    (run / record_name).write_text(json.dumps(record))
    monkeypatch.chdir(run)
    monkeypatch.syspath_prepend(str(run))
    for name in [m for m in sys.modules if m == "amset" or m.startswith("amset.") or m in
                 ("amset_dpo3d", "amset2d_plugin", "adp_degen", "dp_optical", "amset_ir_fix", "amset_fermi_fix")]:
        monkeypatch.delitem(sys.modules, name)
    return run


def _patch_runtime(monkeypatch, sy, off):
    """vasprun 的全网格能带（DFT 比 AMSET 多 off 条深带：按价带顶对齐带号）与点群，换成合成数据。"""
    import adp_degen as rt                                          # 运行目录里那份（插件 import 的同一个）
    E = np.concatenate([np.full((off, sy.E.shape[1]), -60.0), sy.E])
    monkeypatch.setattr(rt, "load_dft_nodes", lambda path="vasprun.xml", symprec=0.01: {
        "kpoints": sy.kf, "E": {"up": E}, "vb": {"up": (sy.d if sy.sign < 0 else 0) + off}})
    monkeypatch.setattr(rt, "point_rotations", lambda structure, symprec=0.01: sy.ops)
    return rt


def test_plugin_3d_scales_adp(si, tmp_path, monkeypatch):
    rec = {"mechanisms": [], "dp_optical": None,
           "adp_degen": {"mode": "on", "status": "pending_runtime", "efit_eV": 0.15, "nang": 24,
                         "deform_h5": "deformation.h5", "signed_h5": "deformation_signed.h5"}}
    run = _prep_run(tmp_path, monkeypatch, [KE / "step8_amset" / "amset_dpo3d.py", S84 / "dp_optical.py",
                                            KE / "adp_degen.py"], "dp_optical_3d.json", rec)
    _write_h5(si, run)
    with redirect_stdout(io.StringIO()) as buf:
        import amset_dpo3d as P
        import amset.scattering.calculate as C
    assert C._scattering_mechanisms["ADP"] is P.AcousticDeformation3DIntravalley and "ADP 简并修正：开" in buf.getvalue()
    _patch_runtime(monkeypatch, si, off=2)
    spin, data = _amset_data(si)
    with redirect_stdout(io.StringIO()) as buf:
        adp = P.AcousticDeformation3DIntravalley.from_amset_data({"elastic_constant": _cubic_C()}, data)
    assert "ADP 简并修正（up，空穴谷" in buf.getvalue()
    report = json.loads((run / "adp_degen_fit.json").read_text())
    assert report["spins"]["up"]["status"] == "ok" and report["dim"] == 3
    uq, nq = np.tile([[1.0, 0, 0]], (5, 1)), np.full(5, 0.01)
    k_near = np.array([0.02, 0.01, 0.0])
    f3 = adp.factor(uq, nq, spin, 3, k_near, np.zeros(3))
    f4 = adp.factor(uq, nq, spin, 4, k_near, np.zeros(3))
    assert np.all(f3 > 1.2) and np.allclose(f4, 1.0)
    re = P.AcousticDeformation3DIntravalley.from_reference(*adp.to_reference())
    assert np.allclose(re.factor(uq, nq, spin, 3, k_near, np.zeros(3)), f3)
    # 不适用（没有波函数重叠）：mode on -> 报错；auto -> 照原版
    data.overlap_calculator = object()
    with pytest.raises(RuntimeError, match="不适用"):
        with redirect_stdout(io.StringIO()):
            P.AcousticDeformation3DIntravalley.from_amset_data({"elastic_constant": _cubic_C()}, data)
    P._DEGEN_REC["mode"] = "auto"
    with redirect_stdout(io.StringIO()) as buf:
        adp2 = P.AcousticDeformation3DIntravalley.from_amset_data({"elastic_constant": _cubic_C()}, data)
    assert adp2._degen is None and "照原版" in buf.getvalue()


def test_plugin_2d_scales_adp(tmp_path, monkeypatch):
    Hm, Ym = e_doublet(-2.0, 0.8), y_doublet(-3.0, 2.5, 0.5)
    sy = Synth(hex_B(), _d6h(), Hm, Ym, 2, n=18, m=12, nd=36, dim=2)
    rec = json.loads((S84 / "example_2d_correction.json").read_text())
    rec.update(mechanisms_2d=["ADP"], adp_degen={"mode": "on", "status": "pending_runtime", "efit_eV": 0.15, "nang": 48,
                                                 "deform_h5": "deformation.h5", "signed_h5": "deformation_signed.h5"})
    run = _prep_run(tmp_path, monkeypatch, [S84 / "amset2d_plugin.py", S84 / "dp_optical.py", KE / "adp_degen.py"],
                    "2d_correction.json", rec)
    _write_h5(sy, run)
    with redirect_stdout(io.StringIO()) as buf:
        import amset2d_plugin as P
        import amset.scattering.calculate as C
    assert C.AcousticDeformationPotentialScattering is P.AcousticDeformation2D
    assert "ADP 简并修正（[0061]）：开" in buf.getvalue()
    _patch_runtime(monkeypatch, sy, off=0)
    spin, data = _amset_data(sy)
    with redirect_stdout(io.StringIO()) as buf:
        adp = P.AcousticDeformation2D.from_amset_data({"elastic_constant": _hex_C()}, data)
    assert "[amset2d] ADP 简并修正（up，空穴谷" in buf.getvalue()
    assert json.loads((run / "adp_degen_fit.json").read_text())["dim"] == 2
    uq, nq = np.tile([[1.0, 0, 0]], (4, 1)), np.full(4, 0.01)
    k = np.array([0.02, 0.01, 0.0])
    keep = adp._iv_keep(uq, nq, spin, 1, k)
    assert keep.shape == (4,) and np.all(keep > 1.1)
    assert adp._iv_keep(uq, nq, spin, 3, k) is None                 # 导带：不修正、不掩
    re = P.AcousticDeformation2D.from_reference(*adp.to_reference())
    assert np.allclose(re._iv_keep(uq, nq, spin, 1, k), keep)


# ---------------------------------------------------------------------------
# S8 / S8.4 gen：开关、前置条件、装模块、链接带符号文件；S8.3 汇总
# ---------------------------------------------------------------------------
def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _resolve(fn, *a):
    buf = io.StringIO()
    with redirect_stdout(buf):
        rec = fn(*a)
    return rec, buf.getvalue()


@pytest.fixture()
def g10(tmp_path):
    g = _load("g10_adpdegen", KE / "step8_amset" / "gen_step10_amset.py")
    g.SCATTERING = ["ADP", "IMP", "POP"]
    g.ADP_DEGEN = "auto"
    cwd = tmp_path / "ke"
    (cwd / g.READ_DIR).mkdir(parents=True)
    return g, cwd


def test_gen_s8_resolve_modes(g10):
    g, cwd = g10
    assert _resolve(g.resolve_adp_degen_3d, True, cwd)[0] is None          # 2D 由 S8.4 做
    rec, log = _resolve(g.resolve_adp_degen_3d, False, cwd)
    assert rec["status"] == "no_signed_h5" and "S7.1" in log
    g.ADP_DEGEN = "on"
    with pytest.raises(SystemExit, match="ADP_DEGEN = on"):
        _resolve(g.resolve_adp_degen_3d, False, cwd)
    g.ADP_DEGEN = "off"
    assert _resolve(g.resolve_adp_degen_3d, False, cwd)[0]["status"] == "off"
    g.ADP_DEGEN = "auto"
    (cwd / g.READ_DIR / "deformation_signed.h5").write_bytes(b"x")
    (cwd / g.READ_DIR / "band_edges.json").write_text(json.dumps({"deform_geometry": "clamped"}))
    rec, log = _resolve(g.resolve_adp_degen_3d, False, cwd)
    assert rec["status"] == "pending_runtime" and rec["deform_geometry"] == "clamped" and "钳位离子" in log
    assert rec["vasprun"] == "vasprun.xml" and rec["nang"] is None
    (cwd / "step3_uniform").mkdir()
    (cwd / "step3_uniform" / "INCAR").write_text("LSORBIT = .TRUE.\n")
    assert _resolve(g.resolve_adp_degen_3d, False, cwd)[0]["status"] == "pending_runtime"   # SOC 也做
    g.SCATTERING = ["IMP", "POP"]
    assert _resolve(g.resolve_adp_degen_3d, False, cwd)[0]["status"] == "no_adp"
    with pytest.raises(SystemExit):
        g.ADP_DEGEN = "maybe"
        _resolve(g.resolve_adp_degen_3d, False, cwd)


def test_gen_s8_installs_module(g10, tmp_path):
    g, _cwd = g10
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
    assert "ADP_DEGEN" in g.SPEC


@pytest.fixture()
def g14(tmp_path, monkeypatch):
    for d in (ROOT, KE, S84):
        monkeypatch.syspath_prepend(str(d))
    g = _load("g14_adpdegen", S84 / "gen_step14_amset2d.py")
    monkeypatch.setattr(g, "SCATTERING", ["ADP", "IMP", "POP"])
    cwd = tmp_path / "ke2d"
    (cwd / g.READ_DIR).mkdir(parents=True)
    return g, cwd


def test_gen_s84_resolve_and_install(g14, tmp_path, monkeypatch):
    g, cwd = g14
    rd = cwd / g.READ_DIR
    for f in ("deformation.h5", "deformation_vac.h5", "deformation_signed.h5"):
        (rd / f).write_bytes(b"x")                                    # 带符号文件只有 core 口径：与真空 h5 不配
    rec, log = _resolve(g.resolve_adp_degen_2d, cwd)
    assert rec["status"] == "no_signed_h5" and "deformation_vac_signed.h5" in log
    (rd / "deformation_vac_signed.h5").write_bytes(b"x")
    rec, log = _resolve(g.resolve_adp_degen_2d, cwd)
    assert rec["status"] == "pending_runtime" and rec["signed_src"] == "deformation_vac_signed.h5"
    monkeypatch.setattr(g, "DEFORM_REF", "core")
    assert _resolve(g.resolve_adp_degen_2d, cwd)[0]["signed_src"] == "deformation_signed.h5"
    monkeypatch.setattr(g, "UNITY_OVERLAP", True)
    assert _resolve(g.resolve_adp_degen_2d, cwd)[0]["status"] == "unity_overlap"
    monkeypatch.setattr(g, "ADP_DEGEN", "on")
    with pytest.raises(SystemExit, match="UNITY_OVERLAP"):
        _resolve(g.resolve_adp_degen_2d, cwd)
    monkeypatch.setattr(g, "UNITY_OVERLAP", False)
    monkeypatch.setattr(g, "SCATTERING", ["POP"])
    with pytest.raises(SystemExit, match="ADP"):
        _resolve(g.resolve_adp_degen_2d, cwd)
    out = tmp_path / "out84"
    out.mkdir()
    monkeypatch.setattr(g, "_ADP_DEGEN_REC", {"status": "no_signed_h5"})
    assert g._install_adp_degen(out) is False and not (out / "adp_degen.py").exists()
    monkeypatch.setattr(g, "_ADP_DEGEN_REC", {"status": "pending_runtime"})
    assert g._install_adp_degen(out) is True and (out / "adp_degen.py").is_file()
    src = (S84 / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert '"adp_degen": _ADP_DEGEN_REC' in src
    assert 'link(out, cwd / READ_DIR / _ADP_DEGEN_REC["signed_src"], "deformation_signed.h5")' in src
    assert "ADP_DEGEN" in g.SPEC


def test_skill_yaml_gen_need():
    yml = (KE / "skill.yaml").read_text(encoding="utf-8")
    for gen in ("gen: gen_step10_amset.py", "gen: gen_step14_amset2d.py"):
        blk = yml[yml.index(gen):]
        assert "adp_degen.py" in blk[:blk.index("fetch_all")], gen


def test_s83_summary_lines(tmp_path):
    g13 = _load("g13_adpdegen", KE / "step8.3_output" / "gen_step13_output.py")
    v = {"status": "ok", "carrier": "hole", "k0": [0, 0, 0], "d": 3, "blocks": [3], "ratio_mean": [[1.9, 3.3, 2.2]],
         "kp": {"rel_res": 0.01}, "strain": {"k0_rel_res": 0.002}}
    (tmp_path / g13.AMSET_DIR).mkdir()
    (tmp_path / g13.AMSET_DIR / "adp_degen_fit.json").write_text(json.dumps(
        {"spins": {"up": {"status": "ok", "valleys": [v, {"status": "谷中心不简并"}]}}}))
    (tmp_path / g13.AMSET2D_DIR).mkdir()
    (tmp_path / g13.AMSET2D_DIR / "adp_degen_fit.json").write_text(json.dumps(
        {"spins": {"up": {"status": "没有找到带边谷"}}}))
    lines = g13._adp_degen_lines(tmp_path)
    assert any("S8（自旋 up）：空穴谷" in x and "1.90/3.30/2.20" in x for x in lines)
    assert any("另有 1 个谷不修正" in x for x in lines) and any("S8.4（自旋 up）：不适用" in x for x in lines)


# ---------------------------------------------------------------------------
# [0061-p] 逐对替换矩阵元：在 AMSET 的四面体积分里换子空间内 (初态, 末态) 的矩阵元
# ---------------------------------------------------------------------------
BOHR = 0.529177210903
H2M_BOHR = 13.605693122994                      # ħ²/2m_e（eV·bohr²）：AMSET 运行期 k 用 bohr⁻¹


@pytest.mark.parametrize("case", ["si3d", "hex2d"])
def test_pairwise_matches_ratio_table(case, monkeypatch):
    """同一组等能面末态上：Σ_末态 F_full·w·(1 − v·v′/|v|²) ÷ Σ AMSET 写法·重叠·w·(…) 必须等于比值表（两条路独立实现，
    同一口径：|固定分数坐标形变势| + m_e v⊗v、动量弛豫因子）。运行期单位：k 用 bohr⁻¹、速度原子单位、能量 Hartree。"""
    if case == "si3d":
        B = SI_B * BOHR
        Hm, Ym, d, dim, ops, C = dkk(*KP_SI, h2m=H2M_BOHR), bir_pikus(*BP_SI), 3, 3, _oh(), _cubic_C()
        sy = Synth(B, ops, Hm, Ym, d, n=12, m=8, nd=16)
    else:
        B = hex_B() * BOHR
        Hm, Ym, d, dim, ops, C = e_doublet(-2.0, 0.8, h2m=H2M_BOHR), y_doublet(-3.0, 2.5, 0.5), 2, 2, _d6h(), _hex_C()
        sy = Synth(B, ops, Hm, Ym, d, n=18, m=12, nd=24, dim=2)
    tab, rep = ad.setup(sy.kf, sy.E, sy.coeff_fn, sy.d, B, sy.ops, C, sy.sign_h5, sy.amp_h5,
                        dense=(sy.dense_kf, sy.dense_E), dim=dim, nang=24, kunit="bohr")
    assert rep["status"] == "ok"
    v = tab["valleys"][0]
    bi = next(i for i, b in enumerate(v["blocks"]) if len(b) == d)
    mb = v["block_models"][bi]
    Y = ad.tensor_Y(v["Z"], dim)
    E = 0.02
    nang = 24
    tt = ad.rate_ratio_table(mb, Y, C, -1, [E], nang, kunit="bohr")
    dirs = ad.directions(nang, dim)
    ts = np.linspace(0.0, mb["rmax"], ad.NT + 1)
    eps = (-np.linalg.eigvalsh(ad.model_H(mb, dirs[:, None, :] * ts[None, :, None])))[..., ::-1]
    t, ok = ad._crossings(eps, ts, E)
    kk, cc, ww, ok = ad._ray_states(mb, dirs, t, ok, E, -1, mb["rmax"])
    states = np.argwhere(ok)
    kf, cf, wf = kk[ok], cc[ok], ww[ok]
    g = np.real(np.einsum("fa,fxab,fb->fx", cf.conj(), ad.model_dH(mb, kf), cf))      # ∇E（eV·bohr）
    bands_e = [v["bands"][d - 1 - r] for _i, r in states]                             # 载流子序 -> AMSET 带号
    monkeypatch.setitem(ad._PAIR, "ok", True)
    Binv = np.linalg.inv(B)
    for n0 in (0, len(kf) // 3, len(kf) - 1):
        i0, r0 = states[n0]
        sel = np.arange(len(kf)) != n0
        sel &= np.linalg.norm(kf - kf[n0], axis=1) > 1e-9
        q = (kf - kf[n0])[sel]
        ov = (cf.conj() @ cf[n0])[sel]
        monkeypatch.setattr(ad, "_CTX", {"qpoints": q, "mapping": np.arange(len(q)),
                                         "tet_bands": np.array(bands_e)[sel], "tet_overlap": np.abs(ov) ** 2})
        vel = g[n0] / ad.HARTREE_EV
        nq2 = np.sum(q ** 2, axis=1)
        res = ad.pair_factor(tab, bands_e[n0], kf[n0] @ Binv, q / np.sqrt(nq2)[:, None], nq2, C, vel)
        assert res is not None and len(res[0]) == len(q) and np.allclose(res[3], 1.0)  # 全在子空间、都在过渡区以内
        idx, F, ovl, _w = res
        qh = q / np.sqrt(nq2)[:, None]
        Cl, U = ad.acoustic_modes(C, qh, dim)
        Es = 0.5 * (np.einsum("fi,fjl->flij", qh, U) + np.einsum("fj,fil->flij", qh, U))
        D0 = np.real(np.einsum("a,xyab,b->xy", cf[n0].conj(), Y, cf[n0])) / ad.HARTREE_EV
        T0 = 0.5 * (np.outer(kf[n0], vel) + np.outer(vel, kf[n0]))
        f_am = np.sum(np.einsum("flij,ij->fl", Es, np.abs(D0 - T0) + np.outer(vel, vel)) ** 2 / Cl, axis=1)
        mrta = 1.0 - (g[sel] @ g[n0]) / (g[n0] @ g[n0])
        r_pair = np.sum(F * wf[sel] * mrta) / np.sum(f_am * ovl * wf[sel] * mrta)
        assert r_pair == pytest.approx(tt["ratio"][0, i0, r0], rel=1e-8)


_FAKE_TET = (
    "import numpy as np\n"
    "def get_cross_section_values(values, *a, average=True):\n"
    "    values = np.asarray(values)\n"
    "    return values.mean(axis=1) if average else values\n"
    "class TetrahedralBandStructure:\n"
    "    def __init__(self, **kw):\n"
    "        self.__dict__.update(kw)\n"
    "    def get_masks(self, spin, tetrahedra_mask):\n"
    "        return None, None, None, None\n")
_FAKE_CALC = (
    "import numpy as np\n"
    "from amset.electronic_structure.tetrahedron import get_cross_section_values\n"
    "from amset.scattering.elastic import AcousticDeformationPotentialScattering\n"
    "_scattering_mechanisms = {'ADP': AcousticDeformationPotentialScattering}\n"
    "def scattering_worker(*a, **k):\n"
    "    return 'orig'\n"
    "def get_fine_mesh_qpoints(tbs):\n"
    "    return tbs.qpoints, tbs.weights, tbs.mapping\n"
    "def calculate_rate(tbs, elastic_scatterers, spin, b_idx, k_idx, velocity):\n"
    "    # 与 AMSET 0.4.19 / 0.5.1 的 calculate_rate 同样的调用顺序（只留插件钩子用到的几步）\n"
    "    tet_mask = (tbs.tet_bands, None)\n"
    "    property_mask, band_kpoint_mask, band_mask, kpoint_mask = tbs.get_masks(spin, tet_mask)\n"
    "    k = tbs.kpoints[k_idx]\n"
    "    vert_overlap, tet_contributions = tbs.vert_overlap, ()\n"
    "    tet_overlap = get_cross_section_values(vert_overlap, *tet_contributions)\n"
    "    intersections = get_cross_section_values(np.zeros(vert_overlap.shape + (3,)), *tet_contributions, average=False)\n"
    "    qpoints, weights, mapping = get_fine_mesh_qpoints(tbs)\n"
    "    qpoint_norm_sq = np.sum(qpoints ** 2, axis=-1)\n"
    "    unit_q = qpoints / np.sqrt(qpoint_norm_sq)[:, None]\n"
    "    mrta_factor = tbs.mrta\n"
    "    rates = [s.factor(unit_q, qpoint_norm_sq, spin, b_idx, k, velocity) for s in elastic_scatterers]\n"
    "    rates = np.array(rates)\n"
    "    rates *= tet_overlap[mapping] * weights * mrta_factor\n"
    "    return np.sum(rates, axis=-1)\n")


def test_plugin_3d_pairwise_in_amset_integration(si, tmp_path, monkeypatch):
    """AMSET 形状的 calculate_rate（同样的调用顺序）：钩子记下末态带号 / 截面重叠 / 细网格映射，ADP factor 对子空间内的
    末态返回 F_full/重叠，AMSET 再乘重叠 -> 正好是 F_full（过渡区按 w 混合）；子空间外的末态、子空间外的初态原样。"""
    rec = {"mechanisms": [], "dp_optical": None,
           "adp_degen": {"mode": "on", "status": "pending_runtime", "efit_eV": 0.15, "nang": 24, "pairwise": "auto",
                         "deform_h5": "deformation.h5", "signed_h5": "deformation_signed.h5"}}
    run = _prep_run(tmp_path, monkeypatch, [KE / "step8_amset" / "amset_dpo3d.py", S84 / "dp_optical.py",
                                            KE / "adp_degen.py"], "dp_optical_3d.json", rec,
                    stubs={"amset/scattering/calculate.py": _FAKE_CALC,
                           "amset/electronic_structure/tetrahedron.py": _FAKE_TET})
    _write_h5(si, run)
    with redirect_stdout(io.StringIO()) as buf:
        import amset_dpo3d as P
        import amset.scattering.calculate as C
        import amset.electronic_structure.tetrahedron as TT
    assert P._PAIR_ON[0] and "逐对替换矩阵元" in buf.getvalue()
    rt = _patch_runtime(monkeypatch, si, off=0)
    spin, data = _amset_data(si)
    with redirect_stdout(io.StringIO()) as buf:
        adp = P.AcousticDeformation3DIntravalley.from_amset_data({"elastic_constant": _cubic_C()}, data)
    assert "方式：逐对替换矩阵元" in buf.getvalue()
    assert "逐对替换矩阵元" in json.loads((run / "adp_degen_fit.json").read_text())["method"]
    rng = np.random.default_rng(11)
    k0 = np.array([0.02, 0.01, 0.0])
    near = k0 + rng.normal(size=(40, 3)) * 0.02                      # Γ 附近（rcut 以内）
    far = rng.uniform(-0.5, 0.5, size=(10, 3))                       # 离谷远
    kp = np.vstack([near, far])
    bands = np.concatenate([rng.choice([1, 2, 3, 4], size=40), rng.choice([1, 2, 3], size=10)])
    ntet = len(kp)
    tbs = TT.TetrahedralBandStructure(
        tet_bands=bands, vert_overlap=rng.uniform(0.05, 1.0, size=(ntet, 4)), kpoints=np.array([k0]),
        qpoints=(kp - k0) @ SI_B, weights=rng.uniform(0.5, 1.5, ntet), mapping=np.arange(ntet),
        mrta=rng.uniform(0.2, 1.0, ntet))
    vel = np.array([0.01, -0.02, 0.005])
    for b0, inside in ((3, True), (4, False)):
        got = C.calculate_rate(tbs, [adp], spin, b0, 0, vel)
        ovl = tbs.vert_overlap.mean(axis=1)
        base = ovl * tbs.weights * tbs.mrta
        res = rt.pair_factor(adp._degen["up"], b0, k0, tbs.qpoints / np.linalg.norm(tbs.qpoints, axis=1)[:, None],
                             np.sum(tbs.qpoints ** 2, axis=1), np.asarray(adp.elastic_constant, float), vel)
        if not inside:
            assert res is None and got[0] == pytest.approx(np.sum(base))
            continue
        idx, F, ov, w = res
        assert 0 < len(idx) < 40 and np.all(idx < 40) and np.allclose(ov, ovl[idx])   # 带 4、远处的末态不换
        exp = base.copy()
        exp[idx] = (w * F + (1 - w) * ovl[idx]) * tbs.weights[idx] * tbs.mrta[idx]
        assert got[0] == pytest.approx(np.sum(exp), rel=1e-10)
        assert np.sum(exp) != pytest.approx(np.sum(base), rel=1e-3)


def test_plugin_2d_pairwise_source_wiring():
    src = (S84 / "amset2d_plugin.py").read_text(encoding="utf-8")
    i = src.index("class AcousticDeformation2D(_ADP)")
    adp = src[i:src.index("# ---------------- ODP / IVS", i)]
    assert "_AD.matrix_apply(f, keep, " in adp and "d_scale=_ADP_D_SCALE" in adp        # [0061-m] 内含 pair_apply 的退路
    assert "or _PAIR_ON[0]:" in adp                                    # 逐对模式不再乘比值
    assert "_PAIR_ON[0], _why_pair = _AD_pair.install_pair_hooks()" in src
