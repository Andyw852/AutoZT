#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""adp_k0_probe.py —— [0061-m] 不跑 AMSET，几秒钟重演 ADP_DEGEN 矩阵路径在谷中心 k0 的两道核对，看哪一步对不上。

S8 的 adp_degen_fit.json 报"形变势矩阵不用 / k0 点对不上"时用。读 S8 运行目录里的三个文件（S7.1 产物的软链）：
  deformation_signed.h5（带符号中心差分 S，原生约定，带轴 = AMSET；ibands = DFT 带号）
  deformation.h5（AMSET 运行期的幅值 A；dp_offdiag = tensor 时非对角 ×2 换回原生）
  deformation_matrix.h5（deform_wave build 的 X，未对齐）
打印：
  1. 对齐（adp_degen.matrix_grid 同一判据）：每个图样的参考能级常数 c_P、离散度、用到的态数；带号偏移来源；
  1b. 标签核对（label_scan）：带符号文件的每个图样对上 X 的哪个 ±X_Q / ½(±X_Q1 ±X_Q2)；对上别的图样 = AMSET 对称补应变
      标错（S7.1 早于 2026-10-22 的标签改正，立方 Si 的 xy、xz）；
  2. k0 处子空间（离带边 --win 以内、在 k0 简并的带）：每个图样
       · 带符号文件的 S0、幅值文件的 A0；
       · 矩阵块 X（对齐后）的本征值；
       · 由 X 按 S7.1 的做法重演（±δ、按能量排序、能移/δ，adp_degen.forward_k0）得到的 S、A；
     以及与 adp_degen._strain_from_grid 同一口径的相对残差（> 20% 那条路就会被拒）。
用法（S8 运行目录里，AMSET 的 conda 环境或任何有 numpy + h5py 的 python）：
  python adp_k0_probe.py [--dir .] [--spin up] [--carrier hole|electron] [--k0 0,0,0] [--win 0.05]
退出码 0 = 对齐和 k0 残差都在门槛内，1 = 有不过的，2 = 读不了。
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in (str(HERE.parent), str(Path.cwd())):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _find_k(kf, k0, tol=1e-4):
    d = np.asarray(kf, float) - np.asarray(k0, float)
    d = d - np.round(d)
    i = np.where(np.linalg.norm(d, axis=1) < tol)[0]
    return int(i[0]) if len(i) else None


def label_scan(m, S, kS, off, delta):
    """带符号文件的每个图样 q，拿矩阵的 ½(σ1·X_Q1 + σ2·X_Q2)（含单个 ±X_Q）逐一按对齐判据配：返回
    {q: (对自己 X_q 的离散度, 最好的组合文字, 它的离散度, 门槛)}。AMSET 对称补应变标错（S7.1 早于 2026-10-22 的
    标签改正）时，±两侧都标错的图样对上别的 ±X_Q，只错一侧的对上 ½(X_q ± X_Q)。态取所有图样都能用的非简并态。"""
    import adp_degen as ad
    import deform_wave as DW
    S = np.asarray(S, float)
    lut = {kk: j for j, kk in enumerate(DW._keys(kS))}
    jS = np.array([lut.get(kk, -1) for kk in DW._keys(m["kpoints"])])
    Xd = np.real(np.einsum("pknn->pkn", m["X"]))
    E = m["eig"]
    gap = np.full(E.shape, np.inf)
    if E.shape[1] > 1:
        dE = np.abs(np.diff(E, axis=1))
        gap[:, :-1] = np.minimum(gap[:, :-1], dE)
        gap[:, 1:] = np.minimum(gap[:, 1:], dE)
    gap[:, 0] = gap[:, -1] = 0.0
    sb = np.asarray(m["bands"], int) - off
    ok = (jS[:, None] >= 0) & ((sb >= 0) & (sb < S.shape[0]))[None, :] & np.all(np.isfinite(Xd), axis=0)
    ok &= gap > max(DW.GAP_ND, 4 * delta * float(np.nanmax(np.abs(Xd))))
    jk, n = np.where(ok)
    if not len(jk):
        return {}
    xv = Xd[:, jk, n]
    pats = [tuple(P) for P in m["pairs"]]
    cand = []
    for i, Q1 in enumerate(pats):
        for s1 in (1, -1):
            cand.append(("%sX_%s" % ("+" if s1 > 0 else "−", DW.PAT_NAMES[Q1]), s1 * xv[i]))
            for j in range(i + 1, len(pats)):
                for s2 in (1, -1):
                    cand.append(("½(%sX_%s %sX_%s)" % ("+" if s1 > 0 else "−", DW.PAT_NAMES[Q1], "+" if s2 > 0 else "−",
                                                     DW.PAT_NAMES[pats[j]]), 0.5 * (s1 * xv[i] + s2 * xv[j])))

    def _sp(v):
        return float(1.4826 * np.median(np.abs(v - np.median(v))))

    res = {}
    for q in DW.ALL_PATS:
        a, b = q
        sv = S[sb[n], jS[jk], a, b]
        own = _sp(xv[pats.index(q)] + sv) if q in pats else np.inf
        sps = [_sp(c + sv) for _, c in cand]
        ib = int(np.argmin(sps))
        tol = max(ad.ALIGN_TOL_EV, ad.ALIGN_TOL_REL * float(np.median(np.abs(sv - np.median(sv)))))
        res[q] = (own, cand[ib][0], sps[ib], tol)
    return res


def probe(rundir, spin="up", carrier="hole", k0=(0.0, 0.0, 0.0), win=0.05, out=print):
    import adp_degen as ad
    import deform_wave as DW
    rd = Path(rundir)
    sgn = ad.read_signed(rd / "deformation_signed.h5", spin)
    try:
        amp = ad.read_amp(rd / "deformation.h5", spin)
    except (OSError, KeyError):
        amp = None
    m = DW.load_matrix(rd / "deformation_matrix.h5", spin)
    delta = float(sgn.get("delta") or 0.005)
    ib = sgn.get("ibands")
    if ib is not None and len(ib) and np.all(np.diff(ib) == 1):
        off, src = int(ib[0]), "带符号文件的 ibands"
    else:
        off, src = None, None
    ok = True
    out("== 0. 文件：S %s（带 %d 条）、A %s、X 图样 %s、窗口带 %d-%d（DFT 带号）、q_min %.3f"
        % (tuple(np.shape(sgn["S"])), np.shape(sgn["S"])[0], "有（%s）" % amp["offdiag"] if amp else "无",
           ",".join(DW.PAT_NAMES[P] for P in m["pairs"]), m["bands"][0], m["bands"][-1], m["q_min"]))
    if off is None:
        out("   ★ 带符号文件没有 ibands（S7.1 早于 [0061-m] 的 retry）：带号偏移要靠运行期 vasprun 推，这里按 --offset 或 0")
        off, src = 0, "默认 0"
    out("   带号偏移（DFT − 带符号文件）= %d（%s）" % (off, src))
    # 1. 对齐
    c, spread, n = DW.alignment(m, sgn["S"], sgn["kpoints"], off, delta)
    out("== 1. 参考能级对齐（X_nn + S_n 在非简并态上的中位数）")
    if c is None:
        out("   ★ 一个非简并态都对不上（带号偏移错？）")
        return 1
    Xd = np.abs(np.real(np.einsum("pknn->pkn", m["X"])) - np.asarray(c)[:, None, None])
    tol = np.maximum(ad.ALIGN_TOL_EV, ad.ALIGN_TOL_REL * np.nanmedian(Xd.reshape(len(c), -1), axis=1))
    for P, cc, sp, tt in zip(m["pairs"], c, spread, tol):
        flag = "" if sp <= tt else "  ★ 超过门槛 %.3f（矩阵会被拒）" % tt
        ok &= sp <= tt
        out("   %s：c = %+8.4f eV，离散度 %.4f eV%s" % (DW.PAT_NAMES[P], cc, sp, flag))
    out("   用到 %d 个非简并态" % n)
    for d_ in (-1, 1):
        c2, sp2, _n2 = DW.alignment(m, sgn["S"], sgn["kpoints"], off + d_, delta)
        if sp2 is not None:
            out("   （对照：带号偏移 %+d 时离散度最大 %.3f eV —— 应明显大于上面）" % (d_, float(np.max(sp2))))
    # 1b. 标签核对
    scan = label_scan(m, sgn["S"], sgn["kpoints"], off, delta)
    if scan:
        out("== 1b. 标签核对：带符号文件的每个图样实际对上 X 的哪个组合（AMSET 对称补应变标错时剪切会对到别的图样上）")
        bad = []
        for q, (own, best, sb_, tt) in scan.items():
            nm = DW.PAT_NAMES[q]
            wrong = own > tt and sb_ <= tt and best != "+X_" + nm
            bad += [nm] if wrong else []
            out("   S_%s：对 X_%s 离散度 %s；最好 %s %.4f eV%s" % (
                nm, nm, "—" if not np.isfinite(own) else "%.4f" % own, best, sb_, "  ★ 标签错" if wrong else ""))
        if bad:
            ok = False
            out("   ★ %s 是 AMSET 对称补应变标错的（S7.1 早于 2026-10-22 的标签改正）：deformation.h5 / *_signed.h5 的"
                "这些分量都错，S7.1 换成新版后重跑（只是重读，不重算 S7）" % "、".join(bad))
    # 2. k0
    jm = _find_k(m["kpoints"], k0)
    jS = _find_k(sgn["kpoints"], k0)
    if jm is None or jS is None:
        out("   ★ k0 = %s 不在%s里" % (k0, "矩阵文件" if jm is None else "带符号文件"))
        return 2
    E = m["eig"][jm]
    vbm, cbm = None, None
    try:
        import h5py
        with h5py.File(str(rd / "deformation_matrix.h5"), "r") as f:
            vbm, cbm = float(f.attrs["vbm"]), float(f.attrs["cbm"])
    except (OSError, KeyError):
        pass
    if carrier == "hole":
        edge = float(np.max(E[E <= vbm + 1e-3])) if vbm is not None else float(np.max(E))
    else:
        edge = float(np.min(E[E >= cbm - 1e-3])) if cbm is not None else float(np.min(E))
    blk = np.where(np.abs(E - edge) < win)[0]
    dft_b = m["bands"][blk]
    sb = dft_b - off
    out("== 2. k0 = %s：%s带边 %.4f eV，子空间 %d 条（DFT 带 %s -> 带符号文件带 %s），劈裂 %s meV"
        % (tuple(k0), "价带顶" if carrier == "hole" else "导带底", edge, len(blk), list(map(int, dft_b)),
           list(map(int, sb)), np.round(1000 * (E[blk] - edge), 2).tolist()))
    if np.any(sb < 0) or np.any(sb >= np.shape(sgn["S"])[0]):
        out("   ★ 子空间的带不在带符号文件的带轴里（带号偏移错？）")
        return 1
    D0 = np.diag(E[blk] - edge)
    ra, rb, da = [], [], []
    jA = _find_k(amp["kpoints"], k0) if amp else None
    for ip, P in enumerate(m["pairs"]):
        a_, b_ = P
        if not np.all(np.isfinite(m["X"][ip, jm])):
            out("   %s：★ k0 这一格缺（两侧对称像都没抽到）" % DW.PAT_NAMES[P])
            ok = False
            continue
        Xb = (m["X"][ip, jm] - c[ip] * np.eye(len(E)))[np.ix_(blk, blk)]
        Xb = 0.5 * (Xb + Xb.conj().T)
        Sp, Ap = ad.forward_k0(D0, Xb, delta)
        S0 = np.asarray(sgn["S"])[sb, jS, a_, b_]
        line = "   %s：eig(X) %s | S0 %s  由 X 重演 %s" % (
            DW.PAT_NAMES[P], np.round(np.linalg.eigvalsh(Xb), 3).tolist(), np.round(S0, 3).tolist(),
            np.round(Sp, 3).tolist())
        if jA is not None:
            A0 = np.asarray(amp["A"])[sb, jA, a_, b_] * (2.0 if (amp["offdiag"] == "tensor" and a_ != b_) else 1.0)
            line += " | A0 %s  重演 %s" % (np.round(A0, 3).tolist(), np.round(Ap, 3).tolist())
        out(line)
        ra.append([np.sum(Sp) - np.sum(S0)])
        rb.append(Sp - S0)
        da.append(S0)
        out("        子空间求和：S0 %+.3f  由 X %+.3f" % (np.sum(S0), np.sum(Sp)))
    na = float(np.linalg.norm(np.concatenate(da)))
    rel = float(np.linalg.norm(np.concatenate(ra)) / na) if na > 0 else None
    relb = float(np.linalg.norm(np.concatenate(rb)) / na) if na > 0 else None
    out("   k0 相对残差（子空间内求和，与 _strain_from_grid 同口径）：%s（门槛 %.0f%%）；逐带 %s（只作诊断：简并点上按能量排序的"
        "中心差分带二阶污染）" % ("—" if rel is None else "%.1f%%" % (100 * rel), 100 * ad.STRAIN_REL_TOL,
                          "—" if relb is None else "%.1f%%" % (100 * relb)))
    if rel is not None and rel > ad.STRAIN_REL_TOL:
        ok = False
        out("   ★ 矩阵在 k0 与带符号文件对不上：逐图样看上面哪一行差得多（某图样整体差一个常数 = 对齐；正负号反 = 约定；"
            "只有剪切差 = 单侧应变 / 对称操作）")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=".")
    ap.add_argument("--spin", default="up")
    ap.add_argument("--carrier", default="hole", choices=("hole", "electron"))
    ap.add_argument("--k0", default="0,0,0")
    ap.add_argument("--win", type=float, default=0.05)
    a = ap.parse_args(argv)
    try:
        return probe(a.dir, a.spin, a.carrier, tuple(float(x) for x in a.k0.split(",")), a.win)
    except (OSError, KeyError) as e:
        print("[ERROR] 读不了：%s: %s" % (type(e).__name__, e))
        return 2


if __name__ == "__main__":
    sys.exit(main())
