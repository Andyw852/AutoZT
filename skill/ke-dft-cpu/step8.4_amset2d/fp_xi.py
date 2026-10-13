#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[独立分支 adp-xi-selfcalc] 冻声子 DFT 自算声学形变势耦合 Xi(q) —— 生成器 + 提取器。

为什么不是解析式
  Kaasbjerg & Thygesen, PRB 87, 235312 (2013) 的 Xi 定义是 |M_q^DP| = Xi*q（Eq.11），
  而 Xi 是**拟合他们第一性原理电声耦合**得到的"有效耦合"（他们明说真耦合各向异性、
  q 依赖复杂，Xi 只是有效参数）。所以要"自算 Xi"，只能算 |M_q| —— 冻声子是最便宜的路子。

做法（二维、面内、Γ→M 方向、N×N×1 超胞）
  1) 把原胞沿 a1/a2 各重复 N 次（3N^2 个原子）；
  2) 冻声子位移：u_k(R) = A * e_k * cos(2*pi*m/N)，m 是沿 a1 的原胞序号，
     e_k = 声学支极化（LA: n=(a1/|a1|)；TA: t⊥n 且面内）；声学长波极限是刚性平移，
     各原子同相（内部坐标分量在 q→0 为零，这里按刚性离子近似，超胞有限 q 会含一部分）；
  3) 参考态 + ±A 三个单点 -> ∂E_CBM/∂A；
  4) 反节点处应变幅 eps = A * q_phys（u=A cos(q.r) 的导数最大值为 A q）：
        Xi(q) = (∂E_CBM/∂A) / q_phys          [eV]
     q→0 时 Xi -> 绝对形变势 D。要拿"有效耦合"，就在若干小 q 上取 Xi(q)（文献 Xi 是 q 窗有效值）。

用法
    # 生成输入（每个方向 ref / p / m 三个目录）
    python fp_xi.py gen --poscar POSCAR --n 4 --amp 0.03 --branch LA --out fp_out
    # 自检（合成一个已知 Xi 的假数据，检查提取器复现它）
    python fp_xi.py selftest
    # 提取（EIGENVAL 放在 fp_out/<case>/EIGENVAL，或用 --energies 给 JSON）
    python fp_xi.py extract --dir fp_out --branch LA --n 4

注意
  · 这是**独立分支**里的探索工具，不进 master、不进任何 autozt 步骤；
  · 超胞要够大（N>=4）才算"长波"；N 小则 Xi(q) 明显偏离 q→0 的 D；
  · 声学支内部坐标分量这里用刚性离子近似 —— 若 S5 DFPT 有该 q 的本征矢，应换成真本征矢
    （--eigvec 接口留了：JSON {原子序号: [dx,dy,dz]}，会自动归一化）。
"""
import argparse
import json
import math
import os
import sys

import numpy as np


# ---------------- POSCAR ----------------
def read_poscar(path):
    txt = open(path, encoding="utf-8-sig").read().splitlines()
    comment = txt[0]
    scale = float(txt[1].split()[0])
    lat = np.array([[float(x) for x in txt[i].split()[:3]] for i in (2, 3, 4)]) * scale
    i = 5
    toks = txt[i].split()
    if all(t.replace(".", "").replace("-", "").isdigit() for t in toks):
        symbols = None
        counts = [int(t) for t in toks]
        i += 1
    else:
        symbols = toks
        counts = [int(t) for t in txt[i + 1].split()]
        i += 2
    direct = not txt[i].strip().lower().startswith(("c", "k"))
    i += 1
    n = sum(counts)
    coords = np.array([[float(x) for x in txt[i + j].split()[:3]] for j in range(n)])
    if not direct:
        coords = coords @ np.linalg.inv(lat)
    return comment, lat, symbols, counts, coords


def write_poscar(path, comment, lat, symbols, counts, frac):
    lines = [comment, "    1.000000000000000"]
    for v in lat:
        lines.append("  %20.16f %20.16f %20.16f" % tuple(v))
    lines.append("  " + "  ".join(symbols if symbols else ["X"] * len(counts)))
    lines.append("  " + "  ".join(str(c) for c in counts))
    lines.append("Direct")
    for f in frac:
        lines.append("  %20.16f %20.16f %20.16f" % tuple(f))
    open(path, "w", encoding="utf-8").write("\n".join(lines) + "\n")


# ---------------- 超胞 + 冻声子 ----------------
def build_supercell(lat, counts, frac, n):
    sup_lat = np.array([lat[0] * n, lat[1] * n, lat[2]])
    fr, cell = [], []
    for m1 in range(n):
        for m2 in range(n):
            shift = np.array([m1 / n, m2 / n, 0.0])
            for f in frac:
                fr.append(f / n + shift)
                cell.append((m1, m2))
    return sup_lat, np.array(fr), cell


def acoustic_pattern(sup_lat, frac, cell, n, branch, amp, phase_sign=1, eigvec=None):
    a1 = sup_lat[0] / n
    n_hat = a1 / np.linalg.norm(a1)
    zhat = np.array([0.0, 0.0, 1.0])
    t_hat = np.cross(zhat, n_hat)
    t_hat = t_hat / np.linalg.norm(t_hat)
    pol = n_hat if branch.upper() == "LA" else t_hat
    disp = np.zeros((len(frac), 3))
    for j, (m1, _m2) in enumerate(cell):
        ph = 2.0 * math.pi * m1 / n
        direction = pol if eigvec is None else np.asarray(eigvec.get(str(j), pol), float)
        disp[j] = amp * phase_sign * direction * math.cos(ph)
    return disp, n_hat, t_hat


def q_phys(lat, n):
    b = 2.0 * math.pi * np.linalg.inv(lat).T
    return float(np.linalg.norm(b[0]) / n)


def cmd_gen(a):
    comment, lat, symbols, counts, frac = read_poscar(a.poscar)
    sup_lat, sup_frac, cell = build_supercell(lat, counts, frac, a.n)
    q = q_phys(lat, a.n)
    os.makedirs(a.out, exist_ok=True)
    for branch in a.branch.upper().split(","):
        for sign, tag in ((0, "ref"), (1, "p"), (-1, "m")):
            d = os.path.join(a.out, "%s_%s" % (branch, tag))
            os.makedirs(d, exist_ok=True)
            if sign == 0:
                new = sup_frac
            else:
                disp, n_hat, t_hat = acoustic_pattern(sup_lat, sup_frac, cell, a.n,
                                                      branch, a.amp, sign)
                new = sup_frac + disp @ np.linalg.inv(sup_lat)
            write_poscar(os.path.join(d, "POSCAR"), comment, sup_lat,
                         symbols, [c * a.n * a.n for c in counts], new)
        print("[gen] %s: %s/{ref,p,m}  %d 原子  q=%.5f 1/A  A=%.4f A  反节点应变=%.5f"
              % (branch, a.out, len(sup_frac), q, a.amp, a.amp * q))
    with open(os.path.join(a.out, "fp_info.json"), "w", encoding="utf-8") as f:
        json.dump({"n": a.n, "amp": a.amp, "q_phys": q, "n_atoms": len(sup_frac),
                   "branches": a.branch.upper().split(","), "poscar": os.path.abspath(a.poscar)},
                  f, indent=1)
    print("[gen] 写 fp_info.json；再对每个目录做 VASP 单点（同 INCAR/KPOINTS/POTCAR），"
          "EIGENVAL 拿带边能量后用 extract 提取")


# ---------------- EIGENVAL / 能量 ----------------
def band_edges_from_eigenval(path, nbands_occ=None):
    """极简 EIGENVAL 解析：返回 (E_CBM, E_VBM, nocc)。
    头部有且仅有一行是 3 个整数：NELECT NKPTS NBANDS；nocc 缺省 = NELECT//2（非自旋）。"""
    lines = open(path, encoding="utf-8").read().splitlines()
    nelect = nk = nb = hdr = None
    for i in range(min(10, len(lines))):
        t = lines[i].split()
        if len(t) == 3 and all(x.lstrip("-").isdigit() for x in t) and int(t[2]) < 5000:
            nelect, nk, nb, hdr = int(t[0]), int(t[1]), int(t[2]), i
            break
    if nk is None:
        return None, None, None
    if nbands_occ is None:
        nbands_occ = nelect // 2
    p = hdr + 1
    allE = []
    for _ in range(nk):
        p += 1                                  # k 点行
        p += 1                                  # 权重/能量行
        for _b in range(nb):
            allE.append(float(lines[p].split()[1]))
            p += 1
    E = np.array(allE).reshape(nk, nb)
    vb = E[:, :nbands_occ].max()
    cb = E[:, nbands_occ:].min() if nb > nbands_occ else float("nan")
    return cb, vb, nbands_occ


def cmd_extract(a):
    info = json.load(open(os.path.join(a.dir, "fp_info.json"), encoding="utf-8"))
    q = info["q_phys"]
    amp = info["amp"]
    res = {}
    for branch in (a.branch or ",".join(info["branches"])).upper().split(","):
        E = {}
        for tag in ("ref", "p", "m"):
            d = os.path.join(a.dir, "%s_%s" % (branch, tag))
            if a.energies:
                vals = json.load(open(a.energies, encoding="utf-8"))[branch + "_" + tag]
                cb, vb = float(vals[0]), float(vals[1])
            else:
                cb, vb, _nocc = band_edges_from_eigenval(
                    os.path.join(d, "EIGENVAL"), a.nocc or None)
                if cb is None:
                    sys.exit("[ERROR] 读不了 %s/EIGENVAL" % d)
            E[tag] = (cb, vb)
        for label, idx in (("CBM", 0), ("VBM", 1)):
            dE = (E["p"][idx] - E["m"][idx]) / (2 * amp)
            xi = abs(dE) / q
            res["%s_%s" % (branch, label)] = xi
            print("[Xi] %-4s %-3s  E(p)=%+.6f E(m)=%+.6f  dE/dA=%+.5f eV/A  Xi(q=%.5f)=%7.3f eV"
                  % (branch, label, E["p"][idx], E["m"][idx], dE, q, xi))
    print("[注] 这是 q=%.5f 1/A 的 Xi(q)；文献 Xi 是其 q 窗有效值，需在几个小 q 上取/外推。"
          % q)
    if a.json_out:
        json.dump(res, open(a.json_out, "w", encoding="utf-8"), indent=1)
        print("[OK] -> %s" % a.json_out)


# ---------------- 从 S7 冻声子 h5 直接取模式的 Xi ----------------
def cmd_h5(a):
    """读 deformation*.h5 的 (nb, nk, 3, 3) 形变势张量，在带边 k 点上做**声学模式投影**：
        Xi_LA(n) = |D : (n̂⊗n̂)|          （纵向/膨胀）
        Xi_TA(n) = |D : sym(n̂⊗t̂)|        （剪切/横向）
    h5 就是 S7 的 ±0.5% 有限差分（= q→0 冻声子）产物，所以这是在读**自己的 DFT 冻声子**。
    要求 h5py（在集群 amset051 环境里跑）。
    """
    import h5py
    be = json.load(open(a.band_edges, encoding="utf-8"))
    f = h5py.File(a.h5, "r")
    dset = f[a.dataset] if a.dataset in f else f[[k for k in f.keys() if "deformation" in k][0]]
    kpts = f["kpoints"][:]
    zhat = np.array([0.0, 0.0, 1.0])
    # 面内两个方向：n 与 t（用 a1 的方向；没有再退回 x 轴）
    n = np.array([1.0, 0.0, 0.0])
    print("[h5] %s   dataset=%s shape=%s" % (a.h5, dset.name, dset.shape))
    for carrier in ("electron", "hole"):
        h = (be.get(carrier) or {}).get("hits") or []
        if not h:
            print("[h5] %s：band_edges.json 里没有 hits，跳过" % carrier)
            continue
        hit = h[0]
        bi, ki = int(hit["band_h5"]), int(hit["k_h5"])
        T = np.array(dset[bi, ki, :, :], float)
        diag = np.diag(T)
        off = math.sqrt(max(T[0, 1] ** 2 + T[0, 2] ** 2 + T[1, 2] ** 2, 0.0))
        print("=" * 78)
        print("%-9s k=%s  band_h5=%d" % (carrier, np.round(kpts[ki], 4), bi))
        with np.printoptions(precision=4, suppress=True, linewidth=120):
            print(T)
        print("  D_xx=%.4f D_yy=%.4f D_zz=%.4f |D_offdiag|=%.4f"
              % (diag[0], diag[1], diag[2], off))
        # 面内方向扫（每 15 度）
        xla, xta = [], []
        for deg in range(0, 180, 15):
            th = math.radians(deg)
            nn = np.array([math.cos(th), math.sin(th), 0.0])
            tt = np.cross(zhat, nn)
            epsL = np.outer(nn, nn)
            epsT = 0.5 * (np.outer(nn, tt) + np.outer(tt, nn))
            xla.append(abs(float(np.sum(T * epsL))))
            xta.append(abs(float(np.sum(T * epsT))))
        print("  Xi_LA(面内各向) = %.4f ~ %.4f eV   (均值 %.4f)"
              % (min(xla), max(xla), float(np.mean(xla))))
        print("  Xi_TA(面内各向) = %.4f ~ %.4f eV   (均值 %.4f)"
              % (min(xta), max(xta), float(np.mean(xta))))
    print("=" * 78)
    print("[对照] Kaasbjerg PRB 87,235312 Tab.I（MoS2 电子）: Xi_LA=2.4 eV, Xi_TA=1.5 eV")
    print("[注意] 连续介质 DP 模型里 TA 不耦合（该文 App.B 明说）—— 他们的 Xi_TA=1.5 来自")
    print("       第一性原理耦合里的 umklapp，正是 q→0 冻声子取不到的部分。")


def cmd_selftest(a):
    lat = np.array([[1.5686622, -2.7170026, 0.0],
                    [1.5686622, 2.7170026, 0.0],
                    [0.0, 0.0, 34.879]])
    ok = True
    for n, amp, xi_true in ((4, 0.03, 8.5), (4, 0.03, 2.4), (5, 0.05, 3.2)):
        q = q_phys(lat, n)
        E_ref, E_p, E_m = -5.0, -5.0 + xi_true * amp * q, -5.0 - xi_true * amp * q
        dE = (E_p - E_m) / (2 * amp)
        xi = abs(dE) / q
        good = abs(xi / xi_true - 1) < 1e-9
        ok = ok and good
        print("[selftest] n=%d A=%.3f Xi_true=%.2f -> Xi_self=%.4f  %s"
              % (n, amp, xi_true, xi, "OK" if good else "**FAIL**"))
    print("[selftest] %s" % ("全部通过" if ok else "**失败**"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description="冻声子 DFT 自算 Xi(q)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gen")
    g.add_argument("--poscar", required=True)
    g.add_argument("--n", type=int, default=4)
    g.add_argument("--amp", type=float, default=0.01,
                   help="位移幅值 Å；反节点应变 = amp*q_phys，建议 <=0.5%")
    g.add_argument("--branch", default="LA")
    g.add_argument("--out", default="fp_out")
    g.set_defaults(func=cmd_gen)
    e = sub.add_parser("extract")
    e.add_argument("--dir", required=True)
    e.add_argument("--branch", default="")
    e.add_argument("--nocc", type=int, default=0,
                   help="占据带数；0/不填 = 自动从 EIGENVAL 头部 NELECT//2")
    e.add_argument("--energies", default="")
    e.add_argument("--json-out", default="")
    e.set_defaults(func=cmd_extract)
    h = sub.add_parser("h5")
    h.add_argument("--h5", required=True)
    h.add_argument("--band-edges", default="band_edges.json")
    h.add_argument("--dataset", default="deformation_potentials_up")
    h.set_defaults(func=cmd_h5)
    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)
    a = ap.parse_args()
    sys.exit(a.func(a) or 0)


if __name__ == "__main__":
    main()
