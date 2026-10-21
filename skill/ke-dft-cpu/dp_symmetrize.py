#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dp_symmetrize.py —— 形变势场 D(k) 按晶体 Laue 群对称化（V121，2026-09-29；V188 改正，2026-10-13）。

为什么要做：`amset deform create` 不带 -s，12 个应变（6 分量 × ±）各自独立做有限差分（离子弛豫构型），
数值噪声让对称等价的分量不再相等。实测 MoS₂ 单层（D3h）带边：电子 E1_vac xx/yy = 8.69/8.35（差 4%），
空穴 2.41/2.55（差 5.9%）；而 K 点有 C3，面内必须 D_xx = D_yy、D_xy = 0。弹性、介电张量都严格对称
（VASP 已对称化），D 是 ADP 唯一不对称的输入 -> 六方晶体上出现 xx/yy 迁移率差 9–10% 的非物理结果。
AMSET 的 get_symmetrized_strain_mapping 只补"没算过"的旋转等价应变，差异不会被平均掉。

物理：应变 ε 下能移 ΔE_n(k) = D_n(k) : ε，在对称操作 r（笛卡尔旋转）下不变
    => D(r k) = r D(k) rᵀ。对 Laue 群（点群 × {±1}，时间反演给 D(-k) = D(k)）取平均：
        D_sym(k) = (1/|G|) Σ_r  rᵀ D(r k) r            （symmetrize_field）
    D_sym 满足 D_sym(s k) = s D_sym(k) sᵀ（群平均是投影，幂等）。这只对**带符号**张量成立。

[V188] V121 的错：它以为 AMSET 的 h5 存的是带符号张量。实际上 AMSET 0.4.19 / 0.5.1 的
`get_strain_deformation_potential` 是 `np.abs(diff / strain_amount)` —— **逐分量取了绝对值**再对各应变平均，
h5 里是幅值 |D_ij|。对幅值做 rᵀ|D|r：只要 Laue 群里有垂直 x、y、z 的镜面（mmm、4/mmm、m-3m、6/mmm，
MoS₂ 的 -6m2 也是 6/mmm），rᵀ|D|r 的非对角元带着 ± 号被平均，**所有非对角 |D_xy|、|D_xz|、|D_yz| 被平均成 0**
（剪切形变势 -> 横声学支的 ADP 耦合被抹掉；L 谷的 PbTe/SnTe/Ge 的 Ξu/3 整个没了）；
3 次/6 次轴还把 |D_xy| 混进对角元。

改正（本模块现在的做法）：
  ① 幅值只用"带符号置换"类操作（矩阵每行每列恰一个 ±1：mmm 子群、立方的 3 次轴 C3[111] 等），并用 |r| 变换：
     (r D rᵀ)_ij = ±D_π(i)π(j)，所以 |r D rᵀ| = |r| |D| |r|ᵀ —— 对幅值**精确**（symmetrize_field_abs）。
  ② 6 次/3 次轴（六方、三方）把分量线性混合，幅值没了符号就无法精确变换。为此 S7.1 另算一份**带符号中心差分**
     （signed_deformation_potentials：与 AMSET 同一套 k 映射、分量归属，只是不取 |·|；±应变平均 = 中心差分），
     写进 *_signed.h5；对它做全群对称化（精确），再取 |·| 写回幅值 h5（AMSET 运行期本来就取 np.abs）。
  ③ 带符号什么时候不可信：±应变给出相反符号（A = ⟨|d_s|⟩ 远大于 |S| = |⟨d_s⟩|）有两种来源 ——
     (a) 非简并带、D 本身小：d± = −D ∓ ½Kδ（二阶曲率 K）或 −D ± s/δ（未形变参考的残差 s），AMSET 的 |·| 再平均
         给出 max(|D|, ½|K|δ, |s|/δ) —— 一个虚假的底；中心差分 S 把它精确抵消，**S 才对**；
     (b) 简并能级被应变劈裂（Si 的 Γ 价带顶）：排序后的带 E = E0 ± |a|δ，S = 0、A = |a|；导数在这里没有定义，
         AMSET 的 |a| 是合理的幅值。
     所以只在"近简并（与相邻带的间隙 < 应变造成的能移 2δ·A_max + 1 meV）**且**符号相反"的 (带, k) 及其对称像上
     退回 ①（AMSET 的幅值口径），其余一律用 |全群对称化(S)|。Kramers 简并（SOC + 反演）应变不劈裂、符号一致，照用 S。
     （sign_consistent / combine）
  没有 *_signed.h5（旧的 S7.1 产物）时只做 ①，并在报告里写明。
  ④ 约定：`amset deform create` 的剪切应变是张量分量 ε_xy = δ，ΔE = 2 D_xy δ，h5 存的是 −ΔE/δ ——
     **非对角是 2·D_xy**（AMSET 约定，不是张量分量）。置换类操作只在非对角之间换位置，与约定无关；
     ② 的全群平均会把对角与非对角线性混合，必须先换成真张量（非对角 ÷2）再旋转，转回来再 ×2（_to_tensor/_from_tensor）。

h5 原地修改（保留其它 attrs，例如 nspin_norm_fix 的标记），打 attrs["dp_symmetrized"] = 2，重复调用跳过；
attrs = 1 是 V121 旧法对称化过的，无法复原，报告里提示重跑 S7.1。

[V189] 对称化之后 S7.1 再调 to_tensor_h5，把 h5 的非对角原地 ÷2 换成**真张量分量**（attrs["dp_offdiag"] = "tensor"）。
原因：AMSET 运行期（elastic.py 的 ADP，以及 S8.4 的 amset2d_plugin）用 tensordot(q̂⊗û, D) 做缩并；物理上
ΔE = Σ D_ij ε_ij，ε = sym(q̂⊗û)·A，D 对称 => ΔE ∝ Σ D_ij q_i u_j，D 必须是张量分量。AMSET 原生 h5 的非对角是 2·D_xy，
剪切部分的耦合被算大一倍（紧束缚 GaN 沿随机 q̂、横向 û：原生 h5 与真实带边移动差到 53%，张量约定 0.1%）。
amset_ir_fix 的星平均对 D 做 R D Rᵀ，也要求张量。换成张量之后，原版 AMSET、2D 插件、ir_fix 都不用改。
"""
import numpy as np

SYMPREC = 0.01            # 与 AMSET 取对称操作的口径一致
KEY_TOL = 1e-4            # k 点匹配容差（分数坐标）
ATTR = "dp_symmetrized"
ATTR_VERSION = 2          # [V188]；1 = V121 旧法（把幅值当带符号张量旋转）
OFFDIAG_ATTR = "dp_offdiag"   # [V189] "tensor" = 非对角是张量分量 D_xy；没有这个标记 = AMSET 原生（2·D_xy）
# 带符号数据文件：deformation.h5 -> deformation_signed.h5。数据集名**不能**含 "deformation_potentials"：
#   amset.deformation.io.load_deformation_potentials 把名字里含它的数据集都当成自旋道读。
SIGNED_KEY = "signed_"
# ③ 的判据：(A - |S|) <= max(GAP_ATOL, GAP_RTOL * 该 (带,k) 最大分量) 才算符号一致。
#   噪声底（|D| 小于未形变参考残差 s/γ 的分量）落在 GAP_ATOL 里不触发；简并劈裂给出 A ≈ 主值、|S| ≈ 0，必然触发。
GAP_ATOL = 0.2            # eV
GAP_RTOL = 0.2
DEG_ATOL = 1e-3           # eV：近简并判据的绝对项（再加应变造成的能移 2δ·A_max）


def signed_name(path):
    """deformation.h5 -> deformation_signed.h5（同目录）。"""
    s = str(path)
    return (s[:-3] if s.endswith(".h5") else s) + "_signed.h5"


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


def is_signed_permutation(r, tol=1e-6):
    """笛卡尔矩阵是否"带符号置换"：每个元素 ∈ {0, ±1}，每行每列恰一个非零。"""
    a = np.abs(np.asarray(r, float))
    if not np.all((np.abs(a) < tol) | (np.abs(a - 1.0) < tol)):
        return False
    nz = a > 0.5
    return bool(np.all(nz.sum(axis=0) == 1) and np.all(nz.sum(axis=1) == 1))


def _key(kf, n=int(round(1 / KEY_TOL))):
    w = np.mod(np.rint(np.asarray(kf, float) * n).astype(np.int64), n)
    return tuple(int(x) for x in w)


def kpoint_maps(kpoints, structure, rots):
    """maps[o, i] = 网格里 r_o k_i 的下标（模倒格矢），找不到记 -1。"""
    B = np.asarray(structure.lattice.reciprocal_lattice.matrix, float)     # 行 = b_i（含 2π）
    return kpoint_maps_B(kpoints, B, rots)


def kpoint_maps_B(kpoints, B, rots):
    """同 kpoint_maps，直接给倒格矢矩阵 B（行 = b_i），纯 numpy。"""
    B = np.asarray(B, float)
    Binv = np.linalg.inv(B)
    kp = np.asarray(kpoints, float)
    index = {}
    for i, k in enumerate(kp):
        index.setdefault(_key(k), i)
    kc = kp @ B
    maps = np.full((len(rots), len(kp)), -1, int)
    for o, r in enumerate(rots):
        kf = (kc @ np.asarray(r, float).T) @ Binv
        for i in range(len(kp)):
            maps[o, i] = index.get(_key(kf[i]), -1)
    return maps


def _group_average(D, maps, mats, ops):
    D = np.asarray(D, float)
    acc = np.zeros_like(D)
    cnt = np.zeros(D.shape[1])
    for o, r in zip(ops, mats):
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


def symmetrize_field(D, maps, rots):
    """**带符号张量**场 D: (nb, nk, 3, 3) 的全群平均 -> (D_sym, 每个 k 参与平均的操作数)。
    D 必须是真张量（非对角 = D_xy）；AMSET 约定的数据（非对角 = 2·D_xy）先过 _to_tensor。"""
    return _group_average(D, maps, [np.asarray(r, float) for r in rots], range(len(rots)))


_OFF = ~np.eye(3, dtype=bool)


def _to_tensor(X):
    """AMSET 约定（非对角 = 2·D_xy）-> 张量（非对角 = D_xy）。"""
    T = np.array(X, float, copy=True)
    T[..., _OFF] *= 0.5
    return T


def _from_tensor(T):
    X = np.array(T, float, copy=True)
    X[..., _OFF] *= 2.0
    return X


def symmetrize_field_abs(A, maps, rots):
    """**幅值**场 A = |D|（AMSET 的 h5）的对称化：只用带符号置换操作，以 |r| 变换（对幅值精确）。
    返回 (A_sym, 每个 k 参与平均的操作数)。其余操作（3/6 次轴、斜镜面）把分量线性混合，幅值做不了。"""
    ops = [o for o, r in enumerate(rots) if is_signed_permutation(r)]
    mats = [np.abs(np.asarray(rots[o], float)) for o in ops]
    return _group_average(A, maps, mats, ops)


def sign_consistent(A, S, atol=GAP_ATOL, rtol=GAP_RTOL):
    """(nb, nk) 布尔：幅值 A = ⟨|d_s|⟩ 与带符号 S = ⟨d_s⟩ 吻合（±应变给出同号）。
    A ≥ |S| 恒成立；差值超过 max(atol, rtol·该点最大分量) 的分量说明符号不可信（简并劈裂、排序交换）。"""
    A = np.asarray(A, float)
    S = np.asarray(S, float)
    gap = A - np.abs(S)
    tol = np.maximum(atol, rtol * A.max(axis=(-1, -2)))
    return np.all(gap <= tol[..., None, None] + 1e-12, axis=(-1, -2))


def combine(A, S, maps, rots, gap=None, delta=None):
    """幅值 A + 带符号 S（两者都是 AMSET 约定）-> 对称化后的幅值。返回 (out, info)。
    gap：(nb, nk) 与相邻带的最小能隙（未形变构型）；delta：应变幅度。给了它们：只在"近简并且符号相反"处回退
    （见模块说明 ③）；没给（旧文件）：保守地在所有符号相反处回退。回退点连同全部对称像一起回退，保证输出对称。"""
    A = np.asarray(A, float)
    S = np.asarray(S, float)
    Ap, _ = symmetrize_field_abs(A, maps, rots)
    Ss = _from_tensor(symmetrize_field(_to_tensor(S), maps, rots)[0])
    incons = ~sign_consistent(A, S)
    if gap is not None and delta:
        deg = np.asarray(gap, float) < 2.0 * float(delta) * A.max(axis=(-1, -2)) + DEG_ATOL
        bad = incons & deg
    else:
        deg = None
        bad = incons
    bad_all = bad.copy()
    for o in range(len(rots)):
        m = maps[o]
        has = m >= 0
        bad_all[:, has] |= bad[:, m[has]]
    out = np.where(bad_all[..., None, None], Ap, np.abs(Ss))
    info = {"n_points": int(bad_all.size), "n_signed": int((~bad_all).sum()),
            "n_fallback_abs": int(bad_all.sum()), "n_sign_opposite": int(incons.sum()),
            "fallback_rule": "near-degenerate & sign-opposite" if deg is not None else "sign-opposite (no gap data)"}
    if deg is not None:
        info["n_near_degenerate"] = int(deg.sum())
    return out, info


def symmetrize_inplane_pair(xx, yy, structure, symprec=SYMPREC):
    """带边 E1 的 (xx, yy) 按保 z 轴的点群子群投影（band_edges.json 用；那里没有 zz）。
    六方/四方/立方 -> 两者取平均；正交等 -> 不变。只动对角（T = diag(xx, yy, 0)，默认带边 D_xy = 0），
    对角元在置换、镜面下只换位置不变号，3 次轴下 Σ 2cs·D_xy = 0，所以对幅值也成立。"""
    rots = [r for r in laue_rotations(structure, symprec)
            if np.allclose(r[2, :2], 0, atol=1e-6) and np.allclose(r[:2, 2], 0, atol=1e-6)]
    T = np.diag([float(xx), float(yy), 0.0])
    Ts = sum(r @ T @ r.T for r in rots) / len(rots)
    return float(Ts[0, 0]), float(Ts[1, 1])


# ---------------------------------------------------------------------------
# 带符号中心差分（S7.1 在 AMSET 环境里调用；与 AMSET 0.4.19 / 0.5.1 的 potentials.py 同一套映射）
# ---------------------------------------------------------------------------
def _signed_strain_dp(strain, bulk_bs, def_bs, bulk_ref, def_ref):
    """AMSET get_strain_deformation_potential 去掉 np.abs 的版本：(E_bulk - E_def - Δref) / ε。"""
    from amset.deformation.potentials import get_mesh_from_band_structure
    from amset.electronic_structure.kpoints import get_kpoint_indices, get_kpoints_from_bandstructure
    strain = strain.round(5)
    flat = np.asarray(strain).ravel()
    amount = flat[np.abs(flat).argmax()]
    ref_diff = bulk_ref - def_ref
    kpoints = get_kpoints_from_bandstructure(bulk_bs)
    mesh, is_shifted = get_mesh_from_band_structure(bulk_bs)
    keep = get_kpoint_indices(kpoints, mesh, is_shifted=is_shifted)
    dk = get_kpoints_from_bandstructure(def_bs)
    di = get_kpoint_indices(dk, mesh, is_shifted=is_shifted)
    if not set(keep).issubset(set(di)):
        raise RuntimeError("形变构型的 k 点不含 bulk 的全部 k 点（与 AMSET 同一判据）")
    dmap = np.full(np.max(di) + 1, -1)
    dmap[di] = np.arange(len(di))
    sel = dmap[keep]
    return {spin: (e0 - def_bs.bands[spin][:, sel] - ref_diff) / amount
            for spin, e0 in bulk_bs.bands.items()}


def signed_deformation_potentials(bulk_calculation, strain_mapping):
    """与 AMSET calculate_deformation_potentials 相同的分量归属（> 0.25·max 的应变分量）与平均，但不取 |·|。
    ±应变平均即中心差分（未形变参考的残差在其中抵消）。norm 放在自旋循环外：没有 0.4.19 的 ISPIN=2 减半缺陷，
    不需要 nspin_norm_fix。"""
    bulk_bs = bulk_calculation["bandstructure"]
    out = {s: np.zeros(b.shape + (3, 3)) for s, b in bulk_bs.bands.items()}
    norm = np.zeros((3, 3))
    for strain, calc in strain_mapping.items():
        d = _signed_strain_dp(strain, bulk_bs, calc["bandstructure"],
                              bulk_calculation["reference"], calc["reference"])
        st = np.asarray(strain, float)
        loc = np.abs(st) > 0.25 * np.abs(st).max()
        lx, ly = np.where(loc)
        for spin, v in d.items():
            out[spin][:, :, lx, ly] += v[..., None]
        norm += loc
    for spin in out:
        out[spin] /= norm[None, None]
    return out


def band_gaps(bandstructure):
    """{spin: (nb, nk)} 每条带与上下相邻带的最小能隙（eV），近简并判据用。"""
    out = {}
    for spin, e in bandstructure.bands.items():
        e = np.asarray(e, float)
        g = np.full(e.shape, np.inf)
        if e.shape[0] > 1:
            d = np.abs(np.diff(e, axis=0))
            g[:-1] = np.minimum(g[:-1], d)
            g[1:] = np.minimum(g[1:], d)
        out[spin] = g
    return out


def strain_amount(strain_mapping):
    """应变幅度 δ（各应变最大分量的最大值，AMSET 默认 0.005）。"""
    return float(max(np.abs(np.asarray(s, float)).max() for s in strain_mapping.keys()))


def write_signed(path, deformation_potentials, kpoints, ibands=None, gaps=None, delta=None):
    """带符号场写成独立文件（数据集 signed_<spin> / gap_<spin>，不含 "deformation_potentials"，AMSET 不会误读）。
    ibands：与幅值 h5 同一次 extract_bands 的带下标（{spin: idx}），保证带轴一致。"""
    import h5py

    def _ext(d):
        return {s: (np.asarray(v)[ibands[s]] if ibands is not None else np.asarray(v)) for s, v in d.items()}
    with h5py.File(str(path), "w") as f:
        for spin, arr in _ext(deformation_potentials).items():
            f.create_dataset(SIGNED_KEY + spin.name, data=np.asarray(arr, float), compression="gzip")
        if gaps is not None:
            for spin, arr in _ext(gaps).items():
                f.create_dataset("gap_" + spin.name, data=np.asarray(arr, float), compression="gzip")
        f["kpoints"] = np.asarray(kpoints, float)
        if ibands is not None:                     # [0061-m] 带轴 -> DFT 带号（运行期对齐形变势矩阵用）
            for spin, idx in ibands.items():
                f["ibands_" + spin.name] = np.asarray(idx, int)
        if delta:
            f.attrs["strain_amount"] = float(delta)
        f.attrs["note"] = ("signed central-difference deformation potentials (eV, AMSET convention: "
                           "off-diagonal = 2*D_ij); used by dp_symmetrize only, not read by AMSET")
    return str(path)


def _load_signed(signed_path, kp):
    """-> ({spin: (S, gap 或 None)}, delta 或 None, 原因)。"""
    import h5py
    from pathlib import Path
    if not signed_path or not Path(signed_path).is_file():
        return None, None, "没有带符号文件 %s" % signed_path
    with h5py.File(str(signed_path), "r") as g:
        kps = np.array(g["kpoints"], float)
        if kps.shape != np.shape(kp) or not np.allclose(kps, kp, atol=1e-6):
            return None, None, "带符号文件的 k 点与 h5 不一致"
        out = {}
        for k in g.keys():
            if k.startswith(SIGNED_KEY):
                spin = k[len(SIGNED_KEY):]
                out[spin] = (np.array(g[k], float),
                             np.array(g["gap_" + spin], float) if ("gap_" + spin) in g else None)
        delta = g.attrs.get("strain_amount")
        return out, (float(delta) if delta is not None else None), None


def symmetrize_h5(path, signed_path=None, symprec=SYMPREC):
    """原地对称化 AMSET 形变势 h5（幅值）。signed_path（缺省 = signed_name(path)）存在时走"带符号 + 回退"，
    否则只用置换子群。返回报告 dict；已对称化过则 {"skipped": ...}。"""
    import h5py
    signed_path = signed_name(path) if signed_path is None else signed_path
    with h5py.File(str(path), "r+") as f:
        v = int(f.attrs.get(ATTR, 0))
        if v >= ATTR_VERSION:
            return {"skipped": "已对称化过"}
        if str(f.attrs.get(OFFDIAG_ATTR, "")) == "tensor":
            return {"skipped": "已换成张量约定（[V189]）却没对称化：对称化必须在 to_tensor_h5 之前 —— 重跑 S7.1"}
        if v:
            return {"skipped": "这份 h5 被 V121 旧法对称化过（把幅值当带符号张量旋转，非对角被抹成 0），"
                               "无法复原 —— 重跑 S7.1 重新读形变势"}
        raw = np.array(f["structure"])
        s = raw.item() if raw.shape == () else raw
        s = s.decode() if isinstance(s, (bytes, np.bytes_)) else str(s)
        from pymatgen.core import Structure     # 用到才导入：上面几种"跳过"不需要 pymatgen
        st = Structure.from_str(s, fmt="json")
        kp = np.array(f["kpoints"])
        rots = laue_rotations(st, symprec)
        maps = kpoint_maps(kp, st, rots)
        signed, delta, why = _load_signed(signed_path, kp)
        rep = {"n_ops": len(rots), "n_perm_ops": sum(is_signed_permutation(r) for r in rots),
               "n_k": int(len(kp)), "k_missing_images": int((maps < 0).sum()),
               "signed": "used" if signed else why, "datasets": {}}
        for key in [k for k in f.keys() if k.startswith("deformation_potentials_")]:
            A = np.array(f[key])
            spin = key.split("_")[-1]
            S, gap = (signed or {}).get(spin, (None, None))
            if S is not None and S.shape == A.shape:
                out, info = combine(A, S, maps, rots,
                                    gap if (gap is not None and gap.shape == A.shape[:2]) else None, delta)
                info["mode"] = "signed"
            else:
                out, _ = symmetrize_field_abs(A, maps, rots)
                info = {"mode": "abs-perm"}
            scale = float(np.abs(A).max()) or 1.0
            info["max_rel_change"] = round(float(np.abs(out - A).max() / scale), 6)
            rep["datasets"][key] = info
            f[key][...] = out
        f.attrs[ATTR] = ATTR_VERSION
        f.attrs["dp_sym_n_ops"] = len(rots)
        f.attrs["dp_sym_signed"] = 1 if signed else 0
    return rep


def to_tensor_h5(path):
    """[V189] 原地把 h5 的非对角 ÷2（AMSET 原生 2·D_xy -> 张量分量 D_xy），打 attrs["dp_offdiag"] = "tensor"。
    幂等；对角不动（band_edges / S8.2 只读对角，不受影响）。返回报告 dict。"""
    import h5py
    with h5py.File(str(path), "r+") as f:
        if str(f.attrs.get(OFFDIAG_ATTR, "")) == "tensor":
            return {"skipped": "已是张量约定"}
        keys = [k for k in f.keys() if k.startswith("deformation_potentials_")]
        for key in keys:
            f[key][...] = _to_tensor(np.array(f[key]))
        f.attrs[OFFDIAG_ATTR] = "tensor"
    return {"converted": keys}


def is_tensor_convention(path):
    """h5 的非对角是不是张量分量（[V189] 标记）；读不了返回 None。"""
    try:
        import h5py
        with h5py.File(str(path), "r") as f:
            return str(f.attrs.get(OFFDIAG_ATTR, "")) == "tensor"
    except Exception:                                          # noqa: BLE001
        return None


def is_symmetrized(path):
    try:
        import h5py
        with h5py.File(str(path), "r") as f:
            return int(f.attrs.get(ATTR, 0)) >= ATTR_VERSION
    except Exception:                                          # noqa: BLE001
        return None


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print(p, symmetrize_h5(p))
