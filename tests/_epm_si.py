"""[0061-m] 测试台：金刚石结构的经验赝势（EPM）模型 + VASP 格式的 WAVECAR 写出。

为什么用它：平面波基里的哈密顿量完全已知，应变一阶哈密顿量 ∂H/∂ε（固定分数坐标）可以直接对矩阵求导，
任意 k、k′ 之间的矩阵元都有精确值 —— 用来检验 deform_wave 的投影 / 插值，以及运行期的逐对矩阵元。
势：Cohen–Bergstresser 的 Si 三个壳层（x = G²/(2π/a)² = 3、8、11 -> −0.21、0.04、0.08 Ry），每个壳层换成宽 1.2 的
光滑紧支撑鼓包（C∞；应变只把 |G| 挪一点点，导数处处有定义），x > 12.2 处严格为零 —— 平面波基收敛快，
截断球随应变变化（Pulay 效应）造成的误差小；体积按 Ω0/Ω 缩放。
物理参数不求像真 Si，只要价带顶是 Γ 点三重态（金刚石对称性保证）、各处光滑。
"""
import itertools

import numpy as np

H2M = 3.80998212                       # ħ²/2m（eV Å²）
RY = 13.605693122994
C_VASP = 0.262465831
A0 = 5.431
L0 = 0.5 * A0 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
TAU = np.array([0.125, 0.125, 0.125])  # 两个原子在 ±TAU（分数坐标）
_X = np.array([3.0, 8.0, 11.0])
_V = np.array([-0.21, 0.04, 0.08])
_W = 1.2
PATS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
NOCC = 4


def vff(q2):
    """形状因子（eV），q2 = |G|²（Å⁻²）。"""
    x = np.asarray(q2, float) / (2 * np.pi / A0) ** 2
    out = np.zeros_like(x)
    for x0, v in zip(_X, _V):
        t = (x - x0) / _W
        m = np.abs(t) < 1
        out[m] += v * np.exp(1.0 - 1.0 / (1.0 - t[m] ** 2))
    return RY * out


def lattice(E=None):
    """Green–Lagrange 应变 E -> 形变晶格（行）：F = (I + 2E)^½，a_i' = F a_i。"""
    if E is None:
        return L0.copy()
    w, U = np.linalg.eigh(np.eye(3) + 2 * np.asarray(E, float))
    F = (U * np.sqrt(w)) @ U.T
    return L0 @ F.T


def pattern(P, d):
    E = np.zeros((3, 3))
    E[P[0], P[1]] = E[P[1], P[0]] = d
    return E


def recip(L):
    return 2 * np.pi * np.linalg.inv(L).T


def _axis(m):
    return np.concatenate([np.arange(0, m + 1), np.arange(-m, 0)])


def gvecs(kf, L, encut, nbmax=12):
    """VASP 顺序的 Miller 指数（与 deform_wave.Wavecar.gvectors 同一规则）。"""
    b = recip(L)
    Z, Y, X = np.meshgrid(_axis(nbmax), _axis(nbmax), _axis(nbmax), indexing="ij")
    G = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    e = np.sum(((np.asarray(kf) + G) @ b) ** 2, axis=1) / C_VASP
    return G[e < encut]


def hamiltonian(kf, G, L):
    """H_GG'(k)（eV），基 = Miller 列表 G，晶格 L（固定分数坐标的应变就是换 L）。"""
    b = recip(L)
    kc = (np.asarray(kf) + G) @ b
    vol = abs(np.linalg.det(L))
    vol0 = abs(np.linalg.det(L0))
    dG = G[:, None, :] - G[None, :, :]
    q2 = np.sum((dG @ b) ** 2, axis=-1)
    V = vff(q2) * np.cos(2 * np.pi * (dG @ TAU)) * (vol0 / vol)
    V[np.arange(len(G)), np.arange(len(G))] = 0.0                     # V(G=0) 只是常数平移，去掉
    return np.diag(H2M * np.sum(kc ** 2, axis=1)) + V


def dH(kf, G, P, d=1e-4):
    """∂H(k)/∂ε_P（固定分数坐标、固定 Miller 基；原生约定：剪切 ε_ab = ε_ba）。"""
    return (hamiltonian(kf, G, lattice(pattern(P, d))) - hamiltonian(kf, G, lattice(pattern(P, -d)))) / (2 * d)


def solve(kf, L, encut):
    G = gvecs(kf, L, encut)
    e, V = np.linalg.eigh(hamiltonian(kf, G, L))
    return G, e, V.T                                                       # 行 = 态


def mesh(n, tr_reduce=True):
    """Γ 心 n³ 网格（分数坐标 (−½, ½]）；tr_reduce：±k 只留一个（VASP ISYM = 0 的做法）。"""
    k = np.array(list(itertools.product(range(n), repeat=3)), float) / n
    k = k - np.round(k + 1e-9)
    k[np.isclose(k, -0.5)] = 0.5
    if not tr_reduce:
        return k
    out, seen = [], set()
    for x in k:
        key = tuple(np.round(x * n).astype(int) % n)
        mkey = tuple(np.round(-x * n).astype(int) % n)
        if mkey in seen:
            continue
        seen.add(key)
        out.append(x)
    return np.array(out)


def write_wavecar(path, kpts, L, encut, nb, efermi=0.0):
    """VASP 标准版格式（rtag 45200，complex64）。-> {k 下标: (G, e, V)}（写进去的东西）。"""
    data = [solve(k, L, encut) for k in kpts]
    npmax = max(len(G) for G, _e, _V in data)
    recl8 = max(npmax, 4 + 3 * nb, 13)
    recl = recl8 * 8
    with open(path, "wb") as f:
        def rec(vals, dtype=np.float64):
            buf = np.zeros(recl, np.uint8)
            raw = np.asarray(vals, dtype).tobytes()
            buf[:len(raw)] = np.frombuffer(raw, np.uint8)
            f.write(buf.tobytes())
        rec([recl, 1, 45200])
        rec(np.concatenate([[len(kpts), nb, encut], np.asarray(L).ravel(), [efermi]]))
        for k, (G, e, V) in zip(kpts, data):
            occ = (np.arange(nb) < NOCC).astype(float)
            rec(np.concatenate([[len(G)], k, np.stack([e[:nb], np.zeros(nb), occ], 1).ravel()]))
            for n in range(nb):
                rec(V[n].astype(np.complex64), np.complex64)
    return data
