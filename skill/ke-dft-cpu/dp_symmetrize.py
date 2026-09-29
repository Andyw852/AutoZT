#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dp_symmetrize.py —— 形变势场 D(k) 按晶体 Laue 群对称化（V121，2026-09-29）。

为什么要做：`amset deform create` 不带 -s，12 个应变（6 分量 × ±）各自独立做有限差分（离子弛豫构型），
数值噪声让对称等价的分量不再相等。实测 MoS₂ 单层（D3h）带边：电子 E1_vac xx/yy = 8.69/8.35（差 4%），
空穴 2.41/2.55（差 5.9%）；而 K 点有 C3，面内必须 D_xx = D_yy、D_xy = 0。弹性、介电张量都严格对称
（VASP 已对称化），D 是 ADP 唯一不对称的输入 -> 六方晶体上出现 xx/yy 迁移率差 9–10% 的非物理结果。
AMSET 的 get_symmetrized_strain_mapping 只补"没算过"的旋转等价应变，差异不会被平均掉。

做法（投影到对称子空间，对精确对称的数据不改值）：
    应变 ε 下能移 ΔE_n(k) = D_n(k) : ε，在对称操作 s（笛卡尔旋转）下不变
    => D(s k) = s D(k) sᵀ。对 Laue 群（点群 × {±1}，时间反演给 D(-k) = D(k)）取平均：
        D_sym(k) = (1/|G|) Σ_r  rᵀ D(r k) r
    D_sym 满足 D_sym(s k) = s D_sym(k) sᵀ（群平均是投影，幂等）。
AMSET 在 h5 里存的是**带符号**张量（calculate_deformation_potentials），abs 在算散射率时才取，
所以对 h5 里的场做旋转平均是合法的。

h5 原地修改（保留其它 attrs，例如 nspin_norm_fix 的标记），打 attrs["dp_symmetrized"] = 1，重复调用跳过。
"""
import numpy as np

SYMPREC = 0.01            # 与 AMSET 取对称操作的口径一致
KEY_TOL = 1e-4            # k 点匹配容差（分数坐标）
ATTR = "dp_symmetrized"


def laue_rotations(structure, symprec=SYMPREC):
    """点群的笛卡尔旋转矩阵，并上 -r（时间反演），去重。"""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    ops = SpacegroupAnalyzer(structure, symprec=symprec).get_point_group_operations(cartesian=True)
    rots = []
    for op in ops:
        for r in (np.asarray(op.rotation_matrix, float), -np.asarray(op.rotation_matrix, float)):
            if not any(np.allclose(r, x, atol=1e-6) for x in rots):
                rots.append(r)
    return rots


def _key(kf, n=int(round(1 / KEY_TOL))):
    w = np.mod(np.rint(np.asarray(kf, float) * n).astype(np.int64), n)
    return tuple(int(x) for x in w)


def kpoint_maps(kpoints, structure, rots):
    """maps[o, i] = 网格里 r_o k_i 的下标（模倒格矢），找不到记 -1。"""
    B = np.asarray(structure.lattice.reciprocal_lattice.matrix, float)     # 行 = b_i（含 2π）
    Binv = np.linalg.inv(B)
    kp = np.asarray(kpoints, float)
    index = {}
    for i, k in enumerate(kp):
        index.setdefault(_key(k), i)
    kc = kp @ B
    maps = np.full((len(rots), len(kp)), -1, int)
    for o, r in enumerate(rots):
        kf = (kc @ r.T) @ Binv
        for i in range(len(kp)):
            maps[o, i] = index.get(_key(kf[i]), -1)
    return maps


def symmetrize_field(D, maps, rots):
    """D: (nb, nk, 3, 3) -> (D_sym, 每个 k 参与平均的操作数)。"""
    D = np.asarray(D, float)
    acc = np.zeros_like(D)
    cnt = np.zeros(D.shape[1])
    for o, r in enumerate(rots):
        m = maps[o]
        ok = m >= 0
        if not ok.any():
            continue
        acc[:, ok] += np.einsum("ji,bnjk,kl->bnil", r, D[:, m[ok]], r)     # rᵀ D(r k) r
        cnt[ok] += 1
    cnt_safe = np.where(cnt > 0, cnt, 1)
    out = acc / cnt_safe[None, :, None, None]
    out[:, cnt == 0] = D[:, cnt == 0]
    return out, cnt


def symmetrize_inplane_pair(xx, yy, structure, symprec=SYMPREC):
    """带边 E1 的 (xx, yy) 按保 z 轴的点群子群投影（band_edges.json 用；那里没有 zz）。
    六方/四方/立方 -> 两者取平均；正交等 -> 不变。"""
    rots = [r for r in laue_rotations(structure, symprec)
            if np.allclose(r[2, :2], 0, atol=1e-6) and np.allclose(r[:2, 2], 0, atol=1e-6)]
    T = np.diag([float(xx), float(yy), 0.0])
    Ts = sum(r @ T @ r.T for r in rots) / len(rots)
    return float(Ts[0, 0]), float(Ts[1, 1])


def symmetrize_h5(path, symprec=SYMPREC):
    """原地对称化 AMSET 形变势 h5。返回报告 dict；已对称化过则 {"skipped": ...}。"""
    import h5py
    from pymatgen.core import Structure
    with h5py.File(str(path), "r+") as f:
        if int(f.attrs.get(ATTR, 0)):
            return {"skipped": "已对称化过"}
        raw = np.array(f["structure"])
        s = raw.item() if raw.shape == () else raw
        s = s.decode() if isinstance(s, (bytes, np.bytes_)) else str(s)
        st = Structure.from_str(s, fmt="json")
        kp = np.array(f["kpoints"])
        rots = laue_rotations(st, symprec)
        maps = kpoint_maps(kp, st, rots)
        rep = {"n_ops": len(rots), "n_k": int(len(kp)),
               "k_missing_images": int((maps < 0).sum()), "datasets": {}}
        for key in [k for k in f.keys() if k.startswith("deformation_potentials_")]:
            D = np.array(f[key])
            Ds, _cnt = symmetrize_field(D, maps, rots)
            scale = float(np.abs(D).max()) or 1.0
            rep["datasets"][key] = round(float(np.abs(Ds - D).max() / scale), 6)
            f[key][...] = Ds
        f.attrs[ATTR] = 1
        f.attrs["dp_sym_n_ops"] = len(rots)
    return rep


def is_symmetrized(path):
    try:
        import h5py
        with h5py.File(str(path), "r") as f:
            return bool(int(f.attrs.get(ATTR, 0)))
    except Exception:                                          # noqa: BLE001
        return None


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print(p, symmetrize_h5(p))
