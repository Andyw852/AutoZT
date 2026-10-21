#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""deform_wave.py —— [0061-m] DFT 形变势矩阵：S7 形变单点的波函数 -> 应变一阶哈密顿量在能带基上的完整矩阵。

为什么要波函数：AMSET（和 S7.1）的形变势只有本征值的差分 D_n(k) —— 简并 / 近简并的带之间应变的耦合（Bir–Pikus 的
b、d，六方的 D3–D6，二维 Γ 点双重态的剪切项……）在本征值里只剩"劈裂"，非对角元的大小与符号读不出来。把形变构型的态
投影回未形变的能带基，应变的一阶哈密顿量就整个出来了：

    X^(P)_mn(k) = ⟨m,k| ∂H/∂ε_P |n,k⟩       （固定分数坐标，应变图样 P 的原生约定：剪切 ε_ab = ε_ba = δ）

即 DFPT 的应变一阶哈密顿量（约化坐标；Hamann, Wu, Rabe, Vanderbilt, PRB 71, 035117 (2005)）在能带基上的矩阵元，这里用
±δ 有限差分算：每个 k 点把形变构型窗口里的态投影到未形变窗口上，Löwdin 正交化后 H^s = S^-½ O diag(E^s) O† S^-½
（O_ml = ⟨m,0|l,s⟩，S = O O†），X^(P) = (H^+ − H^−)/(δ+ − δ−)。一阶以内与能带的规范、简并子空间里怎么组合都无关；
非简并态的对角元与本征值的中心差分严格相等（S7.1 的 deformation_signed.h5），用来定参考能级（芯势 / 真空）对齐的常数并核对。

三段：
  extract（S7 每个单点作业里，VASP 之后）：WAVECAR -> deform_wave.h5（带边 ± EWIN_K 以内有态的 k 点、能带窗口、截到
          ENCUT_STORE 的平面波系数），随后删 WAVECAR（--delete）。
  build（S7.1，amset deform read 之后）：undeformed + 各形变构型的 deform_wave.h5 -> deformation_matrix.h5（X 未对齐、
          未形变的正交化系数、诊断）。
  StrainGrid（S8 / S8.4 运行期，adp_degen 调）：在 S7 网格上对**算符**做三线性插值 D̂(k) = Σ_i α_i Σ_ab |a,k_i⟩X_ab⟨b,k_i|
          （与规范无关），作用到 AMSET 的波函数系数上。

近似：PAW 赝波函数系数（不含补偿项，与 AMSET 的波函数重叠同口径），截到 ENCUT_STORE；X 是中心差分（O(δ²)）。
不依赖 AMSET（extract / build 只要 numpy + h5py），S7 作业在 AMSET 的 conda 环境里跑。
"""
import argparse
import os
import sys

import numpy as np

VERSION = 1
C_VASP = 0.262465831           # 2m/ħ²（eV⁻¹ Å⁻²），与 VASP / pymatgen 相同
EWIN_K = 0.8                   # eV：只抽带边 ± EWIN_K 以内有态的 k 点（输运窗口；800 K 的 10 kT ≈ 0.7 eV）
EWIN_B = 1.2                   # eV：能带窗口（这个范围里有态的带都要；带号连续）
MARGIN = 6                     # 窗口两侧各多抽几条（偶数：SOC 的 Kramers 对不拆）：build 里给窗口扩简并组的余地，
                               #   形变构型那边投影时也用（窗口边上的带被应变推出窗口也不丢）
PAD = 2                        # build：形变构型抽到的带至少比窗口两侧各多 PAD 条（最低的带 0 下面没有带，不要求）
ENCUT_STORE = 300.0            # eV：系数截断（文件 ENCUT 更低时用文件的）
OCC_FRAC = 0.5                 # 占据 > 0.5 × 最大占据 = 占据态
DEG_TOL = 0.05                 # eV：窗口边界两侧的带在某个 k 上相差小于它 = （近）简并，build 把窗口往外扩（只看未形变构型，
                               #   有上限：各形变构型都还得剩 PAD 条余量）。extract 不扩 —— 各构型各自按能量扩，应变劈开的
                               #   简并组两边扩得不一样，导带里到处是交叉，会一路连锁（[0061-m] 离格点测试失败的原因）
Q_MIN = 0.9                    # 形变构型窗口对未形变窗口的覆盖（投影 S 的最小本征值）低于它 -> 报告里点名
GAP_ND = 0.02                  # eV：对齐 / 核对只用离相邻带 > max(GAP_ND, 4δ|X|) 的非简并态
SPIN_NAMES = ("up", "down")
PAT_NAMES = {(0, 0): "xx", (1, 1): "yy", (2, 2): "zz", (1, 2): "yz", (0, 2): "xz", (0, 1): "xy"}


# ---------------------------------------------------------------------------
# WAVECAR（标准版 / 自旋极化 / 非共线；不支持 vasp_gam）
# ---------------------------------------------------------------------------
def _nbmax(b, encut):
    """每个方向的最大 Miller 指数（与 WaveTrans / pymatgen 同一估计，只用来圈出候选，宁大勿小）。"""
    bmag = np.linalg.norm(b, axis=1)
    out = []
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        phi = np.arccos(np.dot(b[i], b[j]) / (bmag[i] * bmag[j]))
        cr = np.cross(b[i], b[j])
        s = np.dot(b[k], cr) / (bmag[k] * np.linalg.norm(cr))
        n = np.sqrt(encut * C_VASP) / bmag
        n[i] /= abs(np.sin(phi))
        n[j] /= abs(np.sin(phi))
        n[k] /= abs(s)
        out.append(n + 1)
    return np.max(out, axis=0).astype(int)


def _axis(m):
    return np.concatenate([np.arange(0, m + 1), np.arange(-m, 0)])


class Wavecar:
    """流式读 WAVECAR：头信息一次读完，系数按 (自旋, k, 带) 定位记录读，不整个载入。"""

    def __init__(self, path="WAVECAR"):
        self.path = str(path)
        with open(self.path, "rb") as f:
            recl, nspin, rtag = np.fromfile(f, np.float64, 3)
            self.recl, self.nspin, self.rtag = int(recl), int(nspin), int(rtag)
            if self.rtag not in (45200, 45210, 53300, 53310):
                raise ValueError("%s：不认识的 rtag %d" % (self.path, self.rtag))
            f.seek(self.recl)
            nk, nb, encut = np.fromfile(f, np.float64, 3)
            self.a = np.fromfile(f, np.float64, 9).reshape(3, 3)
            self.efermi = float(np.fromfile(f, np.float64, 1)[0])
        self.cdtype = np.complex64 if self.rtag in (45200, 53300) else np.complex128
        self.nk, self.nb, self.encut = int(nk), int(nb), float(encut)
        self.b = 2.0 * np.pi * np.linalg.inv(self.a).T
        self.nbmax = _nbmax(self.b, self.encut)
        recl8 = self.recl // 8
        self._nhdr = -(-(4 + 3 * self.nb) // recl8)
        self._blk = self._nhdr + self.nb
        self.kpoints = np.zeros((self.nk, 3))
        self.nplane = np.zeros(self.nk, int)
        self.eig = np.zeros((self.nspin, self.nk, self.nb))
        self.occ = np.zeros((self.nspin, self.nk, self.nb))
        with open(self.path, "rb") as f:
            for s in range(self.nspin):
                for k in range(self.nk):
                    f.seek(self._rec(s, k) * self.recl)
                    h = np.fromfile(f, np.float64, 4 + 3 * self.nb)
                    if s == 0:
                        self.nplane[k], self.kpoints[k] = int(round(h[0])), h[1:4]
                    eo = h[4:].reshape(self.nb, 3)
                    self.eig[s, k], self.occ[s, k] = eo[:, 0], eo[:, 2]
        n0 = len(self.gvectors(0, check=False))
        if n0 == self.nplane[0]:
            self.ncl = False
        elif 2 * n0 == self.nplane[0]:
            self.ncl = True
        else:
            raise ValueError("%s：第一个 k 点的平面波数 %d 与按 ENCUT 生成的 %d 对不上（vasp_gam 的 WAVECAR 不支持）"
                             % (self.path, self.nplane[0], n0))

    def _rec(self, s, k):
        return 2 + (s * self.nk + k) * self._blk

    def gvectors(self, k, check=True):
        """k 点的 Miller 指数 (npw, 3)，VASP 的顺序（x 最快；每个方向 0..n, −n..−1）。"""
        n = self.nbmax
        Z, Y, X = np.meshgrid(_axis(n[2]), _axis(n[1]), _axis(n[0]), indexing="ij")
        G = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
        e = np.sum(((self.kpoints[k] + G) @ self.b) ** 2, axis=1) / C_VASP
        out = G[e < self.encut]
        if check:
            npw = self.nplane[k] // (2 if self.ncl else 1)
            if len(out) != npw:
                raise ValueError("%s：k 点 %d 的平面波数 %d 与文件里的 %d 对不上" % (self.path, k, len(out), npw))
        return out

    def coeffs(self, s, k, bands):
        """(len(bands), npw) 或非共线 (len(bands), npw, 2)。"""
        n = int(self.nplane[k])
        out = np.zeros((len(bands), n), self.cdtype)
        with open(self.path, "rb") as f:
            for i, b in enumerate(bands):
                f.seek((self._rec(s, k) + self._nhdr + int(b)) * self.recl)
                out[i] = np.fromfile(f, self.cdtype, n)
        if self.ncl:
            return out.reshape(len(bands), 2, n // 2).transpose(0, 2, 1)
        return out


def _band_edges(eig, occ):
    occd = occ > OCC_FRAC * float(np.max(occ))
    if occd.all() or not occd.any():
        return None, None, occd
    return float(np.max(eig[occd])), float(np.min(eig[~occd])), occd


def _gkey(M):
    M = np.asarray(M, np.int64) + 512
    return (M[..., 0] * 1024 + M[..., 1]) * 1024 + M[..., 2]


def extract(wavecar="WAVECAR", out="deform_wave.h5", ewin_k=EWIN_K, ewin_b=EWIN_B, margin=MARGIN,
            encut_store=ENCUT_STORE):
    """WAVECAR -> deform_wave.h5。返回摘要字符串；金属 / 没有带隙时不写文件、返回 None。"""
    import h5py
    w = Wavecar(wavecar)
    vbm, cbm, occd = _band_edges(w.eig, w.occ)
    if vbm is None or cbm <= vbm:
        return None
    top = np.where(occd, w.eig, -np.inf).max(axis=2)
    bot = np.where(~occd, w.eig, np.inf).min(axis=2)
    ksel = np.where(np.any((top >= vbm - ewin_k) | (bot <= cbm + ewin_k), axis=0))[0]
    ecut = min(float(encut_store), w.encut) if encut_store else w.encut
    gk, keep = [], []
    for k in ksel:
        G = w.gvectors(k)
        e = np.sum(((w.kpoints[k] + G) @ w.b) ** 2, axis=1) / C_VASP
        gk.append(G)
        keep.append(e < ecut)
    miller = np.unique(np.concatenate([G[m] for G, m in zip(gk, keep)]), axis=0)
    mkey = _gkey(miller)
    order = np.argsort(mkey)
    bands, inner = {}, {}
    for s in range(w.nspin):
        e = w.eig[s][ksel]
        idx = np.where((e.max(axis=0) >= vbm - ewin_b) & (e.min(axis=0) <= cbm + ewin_b))[0]
        lo, hi = int(idx.min()), int(idx.max())                                # 只按能量；简并组在 build 里处理
        inner[s] = (lo, hi)
        bands[s] = np.arange(max(0, lo - margin), min(w.nb - 1, hi + margin) + 1)
    with h5py.File(str(out), "w") as f:
        f.attrs.update(version=VERSION, encut_file=w.encut, encut_store=ecut, nspin=w.nspin, ncl=int(w.ncl),
                       efermi=w.efermi, vbm=vbm, cbm=cbm, nb=w.nb, ewin_k=ewin_k, ewin_b=ewin_b, margin=margin)
        f["lattice"] = w.a
        f["kpoints"] = w.kpoints[ksel]
        f["kindex"] = ksel
        f["miller"] = miller.astype(np.int32)
        shp = (len(ksel), 0, len(miller)) + ((2,) if w.ncl else ())
        dsets = {}
        for s in range(w.nspin):
            sn = SPIN_NAMES[s]
            f["eig_" + sn] = w.eig[s][ksel]
            f["bands_" + sn] = bands[s]
            f["inner_" + sn] = np.array(inner[s])
            dsets[s] = f.create_dataset("coeff_" + sn, (shp[0], len(bands[s])) + shp[2:], np.complex64,
                                        chunks=(1, len(bands[s])) + shp[2:])
        for j, k in enumerate(ksel):
            pos = order[np.searchsorted(mkey[order], _gkey(gk[j][keep[j]]))]
            for s in range(w.nspin):
                c = w.coeffs(s, k, bands[s])[:, keep[j]]
                row = np.zeros((len(bands[s]),) + shp[2:], np.complex64)
                row[:, pos] = c
                dsets[s][j] = row
    return ("%d/%d 个 k 点，带 %s，%d 个平面波（ENCUT_STORE %.0f eV），VBM %.3f / CBM %.3f eV%s"
            % (len(ksel), w.nk, ", ".join("%s %d-%d" % (SPIN_NAMES[s], bands[s][0], bands[s][-1]) for s in bands),
               len(miller), ecut, vbm, cbm, "，非共线" if w.ncl else ""))


# ---------------------------------------------------------------------------
# build：投影 -> deformation_matrix.h5
# ---------------------------------------------------------------------------
def green_lagrange(L0, L):
    """未形变 / 形变晶格（行向量）-> Green–Lagrange 应变（与 AMSET calculate_deformation + pymatgen 同一口径）。"""
    F = (np.linalg.inv(np.asarray(L0, float)) @ np.asarray(L, float)).T
    return 0.5 * (F.T @ F - np.eye(3))


def strain_pattern(E):
    """应变 -> (图样 (a, b) a ≤ b, 带符号幅度)；不是单一分量（AMSET deform create 的做法）-> None。"""
    E = np.asarray(E, float)
    m = float(np.max(np.abs(E)))
    if m < 1e-8:
        return None
    loc = np.argwhere(np.abs(E) > 0.25 * m)
    ab = {(int(min(i, j)), int(max(i, j))) for i, j in loc}
    if len(ab) != 1:
        return None
    a, b = ab.pop()
    return (a, b), float(E[a, b])


def _isqrt_h(M):
    w, U = np.linalg.eigh(0.5 * (M + M.conj().T))
    return (U * np.maximum(w, 1e-300) ** -0.5) @ U.conj().T, w


def _orth(C):
    """行向量组 Löwdin 正交化（赝波函数系数的格拉姆矩阵不是严格的 I）。-> (正交化后, 格拉姆矩阵)。"""
    C2 = C.reshape(len(C), -1)
    G = C2 @ C2.conj().T
    A, _w = _isqrt_h(G)
    return (A @ C2).reshape(C.shape), G


def project_H(C0o, Cs, Es):
    """未形变窗口（已正交化）C0o (n0, …) 与形变构型的态 Cs (ns, …)、能量 Es (ns,) -> (H (n0, n0), 覆盖度 S 的最小本征值)。"""
    Cso, _G = _orth(Cs)
    O = C0o.reshape(len(C0o), -1).conj() @ Cso.reshape(len(Cso), -1).T
    A, w = _isqrt_h(O @ O.conj().T)
    M = A @ O
    return (M * np.asarray(Es, float)) @ M.conj().T, float(np.min(w))


def _keys(kf, n=10000):
    kf = np.asarray(kf, float)
    return [tuple(x) for x in np.round((kf - np.round(kf)) * n).astype(int) % n]


def build(root, dirs, out="deformation_matrix.h5", und="undeformed", fname="deform_wave.h5", log=print):
    """root：step7_deform；dirs：各形变构型目录（相对 root，ionrelax 优先，与 S7.1 同一份清单）。-> 摘要 dict。"""
    import h5py
    root = os.fspath(root)

    def _open(d):
        return h5py.File(os.path.join(root, d, fname), "r")

    f0 = _open(und)
    L0 = np.array(f0["lattice"])
    ncl = bool(f0.attrs["ncl"])
    nspin = int(f0.attrs["nspin"])
    pats = {}
    for d in dirs:
        p = os.path.join(root, d, fname)
        if not os.path.isfile(p):
            raise FileNotFoundError("%s 缺 %s（S7 没开 DEFORM_WAVE，或作业里抽取失败）" % (d, fname))
        with h5py.File(p, "r") as f:
            sp = strain_pattern(green_lagrange(L0, np.array(f["lattice"])))
        if sp is None:
            log("[deform_wave][WARN] %s 不是单一应变分量，跳过" % d)
            continue
        pats.setdefault(sp[0], []).append((sp[1], d))
    if not pats:
        raise ValueError("没有可用的形变构型")
    plist = sorted(pats)
    k0 = _keys(np.array(f0["kpoints"]))
    common = set(k0)
    fs = {d: _open(d) for P in plist for _a, d in pats[P]}
    for d, f in fs.items():
        common &= set(_keys(np.array(f["kpoints"])))
    jk = [j for j, kk in enumerate(k0) if kk in common]
    if not jk:
        raise ValueError("各构型的 k 点没有交集")
    M0 = np.array(f0["miller"])
    gmap = {}
    for d, f in fs.items():
        Ms = np.array(f["miller"])
        ks = _gkey(Ms)
        o = np.argsort(ks)
        pos = np.searchsorted(ks[o], _gkey(M0))
        pos = np.clip(pos, 0, len(ks) - 1)
        hit = ks[o][pos] == _gkey(M0)
        gmap[d] = (np.where(hit)[0], o[pos[hit]])
    kmaps = {d: {kk: j for j, kk in enumerate(_keys(np.array(f["kpoints"])))} for d, f in fs.items()}
    summ = {"patterns": [PAT_NAMES[P] for P in plist], "n_k": len(jk), "ncl": ncl, "spins": {}}
    with h5py.File(os.fspath(out), "w") as g:
        g.attrs.update(version=VERSION, ncl=int(ncl), nspin=nspin, encut_store=float(f0.attrs["encut_store"]),
                       vbm=float(f0.attrs["vbm"]), cbm=float(f0.attrs["cbm"]))
        g["lattice"] = L0
        g["kpoints"] = np.array(f0["kpoints"])[jk]
        g["miller"] = M0
        g["pairs"] = np.array(plist, np.int32)
        g["amounts"] = np.array([[min(a for a, _d in pats[P]), max(a for a, _d in pats[P])] for P in plist], float)
        for s in range(nspin):
            sn = SPIN_NAMES[s]
            bands = np.array(f0["bands_" + sn])
            lo, hi = (int(x) for x in np.array(f0["inner_" + sn]))
            E0all = np.array(f0["eig_" + sn])[jk]
            slo = max(int(np.min(fs[d]["bands_" + sn])) for d in fs)          # 各形变构型都抽到的带
            shi = min(int(np.max(fs[d]["bands_" + sn])) for d in fs)
            blo = max(int(bands[0]), 0 if slo == 0 else slo + PAD)            # 窗口能到的范围：未形变抽到的带以内，
            bhi = min(int(bands[-1]), shi - PAD)                               #   形变构型两侧各剩 PAD 条
            if bhi < blo:
                raise ValueError("自旋 %s：形变构型抽到的带（%d-%d）盖不住未形变的窗口（%d-%d）" % (sn, slo, shi, lo, hi))
            lo2, hi2 = max(lo, blo), min(hi, bhi)
            if hi2 < lo2:
                raise ValueError("自旋 %s：形变构型抽到的带（%d-%d）盖不住未形变的窗口（%d-%d）" % (sn, slo, shi, lo, hi))
            if (lo2, hi2) != (lo, hi):
                log("[deform_wave][WARN] 自旋 %s：窗口 %d-%d 缩到 %d-%d（形变构型只抽到 %d-%d；加大 S7 的 DEFORM_WAVE_MARGIN）"
                    % (sn, lo, hi, lo2, hi2, slo, shi))

            def _deg(i):                                                       # 带 i 与 i+1 在某个 k 上（近）简并
                return np.min(np.abs(E0all[:, i + 1] - E0all[:, i])) < DEG_TOL

            while lo2 > blo and _deg(lo2 - 1):                                 # 窗口边界不切开简并组（只看未形变构型，
                lo2 -= 1                                                       #   到上限为止；不连锁缩）
            while hi2 < bhi and _deg(hi2):
                hi2 += 1
            cut = [("%d|%d" % (i, i + 1)) for i in (lo2 - 1, hi2) if 0 <= i < E0all.shape[1] - 1 and _deg(i)]
            if cut:
                log("[deform_wave] 自旋 %s：窗口 %d-%d 的边界 %s 在个别 k 上近简并（已到可扩的上限）—— 边界带的矩阵元带规范依赖，"
                    "量级同窗口截断误差" % (sn, lo2, hi2, "、".join(cut)))
            inn = np.where((bands >= lo2) & (bands <= hi2))[0]
            nw = len(inn)
            X = g.create_dataset("X_" + sn, (len(plist), len(jk), nw, nw), np.complex64)
            Cd = g.create_dataset("coeff_" + sn, (len(jk), nw) + f0["coeff_" + sn].shape[2:], np.complex64,
                                  chunks=(1, nw) + f0["coeff_" + sn].shape[2:])
            g["bands_" + sn] = bands[inn]
            g["eig_" + sn] = E0all[:, bands[inn]]
            qmin, augmax = np.inf, 0.0
            for jj, j in enumerate(jk):
                C0o, G0 = _orth(np.array(f0["coeff_" + sn][j])[inn])
                augmax = max(augmax, float(np.max(np.abs(G0 - np.diag(np.diag(G0))))))
                Cd[jj] = C0o
                for ip, P in enumerate(plist):
                    Hs, am = [], []
                    for amt, d in sorted(pats[P]):
                        f = fs[d]
                        js = kmaps[d][k0[j]]
                        i0, i1 = gmap[d]
                        Cs_full = np.array(f["coeff_" + sn][js])
                        Cs = np.zeros(Cs_full.shape[:1] + C0o.shape[1:], np.complex64)
                        Cs[:, i0] = Cs_full[:, i1]
                        Es = np.array(f["eig_" + sn][js])[np.array(f["bands_" + sn])]
                        H, q = project_H(C0o, Cs, Es)
                        qmin = min(qmin, q)
                        Hs.append(H)
                        am.append(amt)
                    if len(Hs) >= 2:
                        Xp = (Hs[-1] - Hs[0]) / (am[-1] - am[0])
                    else:
                        Xp = (Hs[0] - np.diag(E0all[j, bands[inn]])) / am[0]
                    X[ip, jj] = 0.5 * (Xp + Xp.conj().T)
            summ["spins"][sn] = {"bands": [int(bands[inn][0]), int(bands[inn][-1])], "q_min": float(qmin),
                                 "aug_offdiag_max": augmax}
            g.attrs["q_min_" + sn] = float(qmin)
            g.attrs["aug_offdiag_max_" + sn] = augmax
            if qmin < Q_MIN:
                log("[deform_wave][WARN] 自旋 %s：形变构型的窗口对未形变窗口覆盖度最低 %.3f < %.2f —— 窗口边上有带被应变推出，"
                    "加大 S7 step.conf 的 DEFORM_WAVE_MARGIN" % (sn, qmin, Q_MIN))
        summ["one_sided"] = [PAT_NAMES[P] for P in plist if len(pats[P]) < 2]
    f0.close()
    for f in fs.values():
        f.close()
    return summ


# ---------------------------------------------------------------------------
# 运行期：读 deformation_matrix.h5、对齐、在 S7 网格上插值应变算符
# ---------------------------------------------------------------------------
def load_matrix(path, spin_name):
    import h5py
    with h5py.File(os.fspath(path), "r") as f:
        if ("X_" + spin_name) not in f:
            raise KeyError("%s 里没有自旋 %s" % (path, spin_name))
        return {"kpoints": np.array(f["kpoints"], float), "miller": np.array(f["miller"], np.int64),
                "pairs": [tuple(int(x) for x in p) for p in np.array(f["pairs"])],
                "bands": np.array(f["bands_" + spin_name], int), "eig": np.array(f["eig_" + spin_name], float),
                "X": np.array(f["X_" + spin_name]), "coeff": np.array(f["coeff_" + spin_name]),
                "ncl": bool(f.attrs["ncl"]), "q_min": float(f.attrs.get("q_min_" + spin_name, np.nan))}


def alignment(data, S_signed, kS, band_offset, delta):
    """参考能级对齐：逐图样一个常数 c_P（X_aligned = X − c_P·I）。非简并态的对角元 = 本征值中心差分，
    带符号文件存 −dE/dε（已对齐）-> c_P = X_nn + S_n。band_offset：矩阵的带号 − 带符号文件的带下标。
    -> (c (np,), 离散度 (np,) eV, 用到的态数)。"""
    S = np.asarray(S_signed, float)
    lut = {kk: j for j, kk in enumerate(_keys(kS))}
    jS = [lut.get(kk) for kk in _keys(data["kpoints"])]
    E = data["eig"]
    nb = len(data["bands"])
    gap = np.full(E.shape, np.inf)
    if nb > 1:
        dE = np.abs(np.diff(E, axis=1))
        gap[:, :-1] = np.minimum(gap[:, :-1], dE)
        gap[:, 1:] = np.minimum(gap[:, 1:], dE)
    gap[:, 0] = gap[:, -1] = 0.0                                           # 窗口两端的带不知道窗口外的邻居
    c, spread, nused = [], [], 0
    for ip, (a, b) in enumerate(data["pairs"]):
        Xd = np.real(np.einsum("knn->kn", data["X"][ip]))
        xm = np.max(np.abs(Xd)) if Xd.size else 0.0
        vals = []
        for jk, js in enumerate(jS):
            if js is None:
                continue
            for n in range(nb):
                sb = int(data["bands"][n]) - band_offset
                if 0 <= sb < S.shape[0] and gap[jk, n] > max(GAP_ND, 4 * delta * xm):
                    vals.append(Xd[jk, n] + S[sb, js, a, b])
        if not vals:
            return None, None, 0
        vals = np.array(vals)
        med = float(np.median(vals))
        c.append(med)
        spread.append(float(1.4826 * np.median(np.abs(vals - med))))
        nused += len(vals)
    return np.array(c), np.array(spread), nused


def infer_mesh(kf):
    """S7 网格（Γ 心或整体平移）-> (N (3,), shift (3,))：k = (g + shift)/N。"""
    kf = np.asarray(kf, float)
    N = np.ones(3, int)
    for a in range(3):
        x = np.unique(np.round(kf[:, a] - np.floor(kf[:, a]), 8))
        d = np.diff(np.concatenate([x, [x[0] + 1.0]]))
        d = d[d > 1e-6]
        N[a] = int(round(1.0 / float(np.min(d)))) if len(d) else 1
    sh = kf[0] * N
    return N, np.round(sh - np.floor(sh + 1e-6), 6)


class StrainGrid:
    """S7 网格上的应变算符 D̂^(P)(k) = Σ_i α_i Σ_ab |a,k_i⟩ X^(P)_ab(k_i) ⟨b,k_i|（三线性插值算符本身，与能带规范无关）。
    缺的格点用时间反演补（ISYM = 0 只给半个网格）：|a,−k⟩ = Θ|a,k⟩，X(−k) = X(k)*。
    gpoints：AMSET 的 Miller 列表（系数就在它上面）；offsets：参考能级对齐常数 c_P（eV）。X 存 eV，apply 返回 eV。"""

    def __init__(self, data, gpoints, offsets=None):
        self.pairs = list(data["pairs"])
        self.ncl = bool(data["ncl"])
        nw = len(data["bands"])
        X = np.asarray(data["X"])
        if offsets is not None:
            X = X - np.asarray(offsets, float)[:, None, None, None] * np.eye(nw)[None, None]
        self.X = X                                                         # (np, nk, nw, nw)
        self.C = np.asarray(data["coeff"])                                 # (nk, nw, NG[, 2])
        self.miller = np.asarray(data["miller"], np.int64)
        self.kf = np.asarray(data["kpoints"], float)
        self.N, self.shift = infer_mesh(self.kf)
        gp = np.asarray(gpoints, np.int64)
        self._gk = _gkey(gp)
        self._go = np.argsort(self._gk)
        self.ng = len(gp)
        self.lut = {}
        g = np.round(self.kf * self.N - self.shift).astype(int)
        for j, gi in enumerate(g):
            self.lut[tuple(gi % self.N)] = (j, False)
        for j, gi in enumerate(g):                                         # −k = (−g − 2·shift)/N
            self.lut.setdefault(tuple(np.round(-gi - 2 * self.shift).astype(int) % self.N), (j, True))
        self._maps = {}

    def _gmap(self, g, tr):
        key = (tuple(int(x) for x in g), tr)
        if key not in self._maps:
            M = (-self.miller if tr else self.miller) - np.asarray(g, np.int64)
            kk = _gkey(M)
            pos = np.clip(np.searchsorted(self._gk[self._go], kk), 0, self.ng - 1)
            hit = self._gk[self._go][pos] == kk
            self._maps[key] = (np.where(hit)[0], self._go[pos[hit]])
        return self._maps[key]

    def corners(self, kf):
        """k（分数坐标，AMSET 的表示）-> [(权重, 存储下标, 是否时间反演, 倒格矢平移 g)]；缺格点 -> None。"""
        x = np.asarray(kf, float) * self.N - self.shift
        i0 = np.floor(x + 1e-9).astype(int)
        t = x - i0
        out = []
        for bits in np.ndindex(2, 2, 2):
            wt = float(np.prod([t[a] if bits[a] else 1.0 - t[a] for a in range(3)]))
            if wt <= 1e-12:
                continue
            c = i0 + np.array(bits)
            hit = self.lut.get(tuple(c % self.N))
            if hit is None:
                return None
            j, tr = hit
            ks = -self.kf[j] if tr else self.kf[j]
            g = np.round((c + self.shift) / self.N - ks).astype(int)
            out.append((wt, j, tr, g))
        return out

    def _basis(self, j, tr, g):
        """格点的正交化基在 AMSET Miller 列表上的稠密表示的 (非零位置, 系数 (nw, nnz[, 2]))。"""
        i0, i1 = self._gmap(g, tr)
        C = self.C[j][:, i0]
        if tr:
            C = np.conj(C) if not self.ncl else np.stack([np.conj(C[..., 1]), -np.conj(C[..., 0])], -1)
        return i1, C

    def coverage(self, kf, vecs):
        """vecs (n, NG[, 2]) 在 k 处网格基（窗口里的带）上的投影范数² 按插值权重平均 -> (n,)；缺格点 -> None。
        ≈1 = 这组态落在窗口里、且 Miller 约定与网格一致（运行期据此认 AMSET 系数的倒格矢表示）。"""
        cs = self.corners(kf)
        if cs is None:
            return None
        vecs = np.asarray(vecs)
        out = np.zeros(len(vecs))
        for wt, j, tr, g in cs:
            pos, C = self._basis(j, tr, g)
            proj = C.reshape(len(C), -1).conj() @ vecs[:, pos].reshape(len(vecs), -1).T
            out += wt * np.sum(np.abs(proj) ** 2, axis=0)
        return out

    def apply(self, kf, vecs):
        """D̂^(P)(k) 作用到行向量组 vecs (n, NG[, 2]) 上 -> (np, n, NG[, 2])（eV）；缺格点 -> None。"""
        cs = self.corners(kf)
        if cs is None:
            return None
        vecs = np.asarray(vecs)
        out = np.zeros((len(self.pairs),) + vecs.shape, np.complex128)
        for wt, j, tr, g in cs:
            pos, C = self._basis(j, tr, g)
            Cf = C.reshape(len(C), -1)
            vv = vecs[:, pos].reshape(len(vecs), -1)
            proj = Cf.conj() @ vv.T                                        # ⟨b|v⟩ (nw, n)
            Xs = np.conj(self.X[:, j]) if tr else self.X[:, j]
            add = np.einsum("pab,bn->pna", Xs, proj) @ Cf                  # (np, n, nnz·s)
            out[:, :, pos] += wt * add.reshape((len(self.pairs), len(vecs), len(pos)) + vecs.shape[2:])
        return out


# ---------------------------------------------------------------------------
# 命令行（S7 作业 / S7.1 调）
# ---------------------------------------------------------------------------
def _cli(argv=None):
    ap = argparse.ArgumentParser(prog="deform_wave.py")
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract", help="WAVECAR -> deform_wave.h5（S7 作业里，VASP 之后）")
    e.add_argument("--wavecar", default="WAVECAR")
    e.add_argument("--out", default="deform_wave.h5")
    e.add_argument("--ewin-k", type=float, default=EWIN_K)
    e.add_argument("--ewin-b", type=float, default=EWIN_B)
    e.add_argument("--margin", type=int, default=MARGIN)
    e.add_argument("--encut", type=float, default=ENCUT_STORE)
    e.add_argument("--delete", action="store_true", help="抽取成功后删 WAVECAR")
    b = sub.add_parser("build", help="各构型的 deform_wave.h5 -> deformation_matrix.h5（S7.1）")
    b.add_argument("--root", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("dirs", nargs="+")
    a = ap.parse_args(argv)
    if a.cmd == "extract":
        if not os.path.isfile(a.wavecar):
            print("[deform_wave][WARN] 没有 %s：不抽取" % a.wavecar, flush=True)
            return 0
        msg = extract(a.wavecar, a.out, a.ewin_k, a.ewin_b, a.margin, a.encut)
        if msg is None:
            print("[deform_wave] 没有带隙（金属）：不写 %s" % a.out, flush=True)
            return 0
        print("[deform_wave] %s：%s" % (a.out, msg), flush=True)
        if a.delete:
            os.remove(a.wavecar)
            print("[deform_wave] 已删 %s" % a.wavecar, flush=True)
        return 0
    s = build(a.root, a.dirs, a.out, log=lambda m: print(m, flush=True))
    print("[deform_wave] %s：图样 %s，%d 个 k 点%s；%s" % (
        a.out, ",".join(s["patterns"]), s["n_k"], "，非共线" if s["ncl"] else "",
        "；".join("%s 带 %d-%d 覆盖度 ≥ %.3f 赝波函数非正交 ≤ %.3f" % (sn, v["bands"][0], v["bands"][1], v["q_min"],
                                                                v["aug_offdiag_max"]) for sn, v in s["spins"].items())),
        flush=True)
    if s["one_sided"]:
        print("[deform_wave][WARN] 图样 %s 只有单侧应变：前向差分（O(δ)）" % ",".join(s["one_sided"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
