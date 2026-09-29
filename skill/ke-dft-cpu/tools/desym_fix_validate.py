#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""desym_fix_validate.py —— 用**同一份全网格** h5 检验 AMSET 去对称化相位公式（V115 第三种设计）。

思路（只用一份 ISYM=-1 的全网格 wavefunction.h5，例如 step4b_wave_full/）：
  1. spglib.get_ir_reciprocal_mesh 取 IBZ（同一结构、Γ 心、开时间反演、symprec 与 AMSET 相同 0.01）；
  2. 从全网格 h5 里**取出 IBZ 代表点**的系数；
  3. 把**整份** IBZ 列表一次性传给 amset.expand_kpoints，拿全局 op_mapping / kp_mapping
     （逐点单独传只会得到平凡轨道 —— 这是第二次实现的错误）；
  4. 分块调 desymmetrize_coefficients（原式 / 修正式各一遍）展开回全网格；
  5. 与**同一份**全网格系数逐点比较：简并带（能量差 < 1 meV）用子空间投影范数，其余逐带 |cos|。
数据来自同一次 VASP、同一个结构 —— 被检验的只有相位公式本身；TR / 非 TR 两类操作都会测到。

判据（两套，任一通过即退出码 0；JSON 里两套都写）：
  · 严格（verdict_strict，模型数据用）：修正式全部 k 点、全部带最小值 >= 0.999（--pass-cos），
    且原式出错点（< 0.999）**全部**来自"坏操作"（逐操作判据 (R±I)τ≢0，AMSET 自己的操作集）。
  · 对照组（verdict_control，V116，实数据用）：VASP 波函数在对称相关 k 点间本来就不是
    精确对称（近简并混合、结构只在 symprec 内对称等），AMSET **原生路径**在它本来就算对的
    操作上也达不到 1。所以不拍绝对门槛，而是拿同一份数据里的**好操作点**当对照组 ——
    在好操作上原式与修正式逐 G 相位差 exp(-2πi G·(R+I)τ) ≡ 1，二者都等于 AMSET 原生结果：
    记"公式级错误"亏损线 L = max(1-0.99, 3×对照组最大亏损)（随数据地板浮动）：
      a) 修正式在坏操作点上无公式级错误：最大亏损 1-|cos| <= L；
      b) 修正式在坏操作点上的 p90 亏损 <= 2×对照组 p90 亏损 + 1e-4；
      c) 修正式在坏操作点上 < 0.999 的比例 <= 2×对照组比例 + 0.2%；
      d) 原式在好操作点上也无公式级错误（最大亏损 <= L）—— 即公式级出错点全部来自坏操作，
         数值地板不参与归因；
      e) 检验有区分力：原式在坏操作点上亏损 > L 的比例 >= 5%（否则判"无区分力"，不算通过）。
    同一判据也套在原式上跑一遍（orig_would_pass），它必须不通过 —— 证明判据能区分对错。
    同时报出错点的"最近邻带能距"分布（近简并诊断），用来核对地板是否集中在近简并点。
  · --shift t1 t2 t3（V116）：把原点整体平移 t（原子 x->x+t，系数 c_k(G)->c_k(G)·e^{-2πi(k+G)·t}，
    精确变换）再检验。用途：标准原点的 GaN 坏操作 = 0，原样检验只能证明"两式相同"（无区分力）；
    平移后同一份数据出现坏操作 -> 有区分力。另报**平移不变性**：修正式逐点分数应与平移前
    相同（|Δ| < 1e-5），原式则不然 —— 地板若是数据本身的性质，就与原点无关。

用法（在材料目录，AMSET 环境里；不需要提交作业，~5–10 分钟；--shift 约 1.5 倍）：
    python <skill>/ke-dft-cpu/tools/desym_fix_validate.py \
        --h5 step4b_wave_full/wavefunction.h5 --vasprun step4b_wave_full/vasprun.xml
可选：--mesh 48 48 3（默认从 h5 的 k 点推）、--bands 0:8（只比部分带，省时间）、
      --energy-cutoff（与 amset wave 取带时一致；默认 AMSET 默认值）、--json out.json、
      --shift 0.1 0.2 0.05（3D 用三个分量；2D 第三个分量写 0）。
退出码：0 = 严格或对照组判据通过（--shift 时还要求平移不变）；1 = 未通过；2 = 输入问题。

注：h5 恒为 complex128（amset 的 write_coefficients 把 dtype 写死），不能据此排除 WAVECAR
单精度 —— 但单精度舍入带来的 |cos| 亏损 ~1e-14，与 ~1e-3 的地板无关。

本文件同时被 step8.4_amset2d/test_desym_fix.py 当库用（模型数据、无 VASP）。
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for _p in (_HERE.parent / "step8.4_amset2d", _HERE.parent, Path.cwd()):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

AMSET_SYMPREC = 0.01      # amset/defaults.yaml: symprec；from_file 路径不传 symprec -> 恒用它
DEGEN_TOL_EV = 1e-3       # 简并判据：能量差 < 1 meV 归为同一子空间
PASS_COS = 0.999
# ---- V116 对照组判据（见文件头）----
GROSS_COS = 0.99          # 公式级错误的量级：原式在坏操作上典型 |cos| ~ 0.1–0.9，远在地板之下
NEAR_DEGEN_MEV = 25.0     # 近简并诊断：出错点的最差带离最近非简并带 < 25 meV 记为"近简并"
CONTROL_MIN_N = 20        # 对照组（好操作、非平凡映射）至少这么多点才有统计意义
DISCRIM_MIN_FRAC = 0.05   # 原式在坏操作点上 < GROSS_COS 的比例 >= 5% 才算"有区分力"
SHIFT_TOL = 1e-5          # 平移不变性：修正式逐点分数平移前后差 < 1e-5


# --------------------------------------------------------------------------
# k 点工具
# --------------------------------------------------------------------------
def _wrap(k):
    k = np.asarray(k, float)
    return k - np.rint(k)


def match_kpoints(query, ref, tol=1e-5):
    """query 中每个 k 在 ref 中的下标（模 1 比较）；找不到抛 ValueError。"""
    q = _wrap(query)
    r = _wrap(ref)
    out = np.empty(len(q), int)
    for i, k in enumerate(q):
        d = r - k
        d -= np.rint(d)
        j = np.where(np.all(np.abs(d) < tol, axis=1))[0]
        if len(j) == 0:
            raise ValueError("k 点 %s 不在参考网格里" % np.round(k, 6).tolist())
        out[i] = j[0]
    return out


def mesh_from_kpoints(kpoints):
    """Γ 心整网格的 (n1,n2,n3)：每个方向取 1/最小非零间隔。"""
    k = _wrap(kpoints)
    mesh = []
    for ax in range(3):
        u = np.unique(np.round(k[:, ax], 6))
        if len(u) == 1:
            mesh.append(1)
            continue
        d = np.diff(u)
        d = d[d > 1e-5]
        mesh.append(int(round(1.0 / d.min())))
    return mesh


def ibz_from_mesh(structure, mesh, symprec=AMSET_SYMPREC, time_reversal=True):
    """spglib get_ir_reciprocal_mesh（Γ 心）-> IBZ 代表点的分数坐标（落在 [-0.5,0.5)）。"""
    import spglib
    cell = (structure.lattice.matrix, structure.frac_coords,
            [s.Z for s in structure.species])
    mapping, grid = spglib.get_ir_reciprocal_mesh(
        np.asarray(mesh, int), cell, is_shift=[0, 0, 0],
        is_time_reversal=bool(time_reversal), symprec=symprec)
    reps = np.unique(mapping)
    return _wrap(np.asarray(grid[reps], float) / np.asarray(mesh, float))


# --------------------------------------------------------------------------
# 比较
# --------------------------------------------------------------------------
def degenerate_groups(energies, tol=DEGEN_TOL_EV):
    """energies: (nb,) 某 k 点的能量 -> [[band idx, ...], ...]（按能量链式聚类，< tol 相连）。"""
    if energies is None:
        return None
    e = np.asarray(energies, float)
    order = np.argsort(e)
    groups, cur = [], [int(order[0])]
    for a, b in zip(order[:-1], order[1:]):
        if e[b] - e[a] < tol:
            cur.append(int(b))
        else:
            groups.append(cur)
            cur = [int(b)]
    groups.append(cur)
    return groups


def band_gaps_mev(e_all, sel, tol=DEGEN_TOL_EV):
    """每条被比较带到"最近的非简并带"（|ΔE| >= tol）的能距（meV）。e_all: 该 k 点 h5 全部带能量。"""
    e_all = np.asarray(e_all, float)
    out = np.full(len(sel), np.inf)
    for j, b in enumerate(sel):
        d = np.abs(e_all - e_all[b])
        d = d[d >= tol]
        if len(d):
            out[j] = d.min() * 1000.0
    return out


def _stats(scores, pass_cos=PASS_COS, gross=GROSS_COS):
    """一组分数的摘要（亏损 = 1-|cos|）。"""
    x = np.asarray(scores, float)
    if len(x) == 0:
        return {"n": 0}
    d = 1.0 - x
    return {"n": int(len(x)), "min": float(x.min()),
            "max_deficit": float(d.max()), "p99_deficit": float(np.percentile(d, 99)),
            "p90_deficit": float(np.percentile(d, 90)),
            "median_deficit": float(np.median(d)),
            "frac_below_pass": float((x < pass_cos).mean()),
            "frac_below_gross": float((x < gross).mean())}


def control_checks(worst, worst_orig, classes, pass_cos=PASS_COS):
    """V116 对照组判据（见文件头 a–e）。

    worst / worst_orig：被检公式 / 原式在每个目标 k 上的最差分数；classes：validate 里的类别掩码
    （good = 好操作非平凡映射 = 对照组，bad = 坏操作非平凡映射）。
    "公式级错误"的亏损线 L = max(1-GROSS_COS, 3×对照组最大亏损) —— 随数据地板浮动，
    这样地板本身偏高（数据差）时不会把地板误当成公式错误；地板偏高另在 note 里提示。
    返回 dict(ok=True/False/None, checks, gross_deficit_line, note)；ok=None = 无从判定。
    """
    dg = 1.0 - np.asarray(worst)[classes["good"]]
    db = 1.0 - np.asarray(worst)[classes["bad"]]
    do_g = 1.0 - np.asarray(worst_orig)[classes["good"]]
    do_b = 1.0 - np.asarray(worst_orig)[classes["bad"]]
    if len(db) == 0:
        return {"ok": None, "note": "没有坏操作点 —— 两式相同，检验无区分力（用 --shift 造坏操作）"}
    if len(dg) < CONTROL_MIN_N:
        return {"ok": None, "note": "对照组（好操作、非平凡映射）只有 %d 点 < %d，无从比较"
                                    % (len(dg), CONTROL_MIN_N)}
    L = max(1.0 - GROSS_COS, 3.0 * float(dg.max()))
    thr = 1.0 - pass_cos
    checks = {
        "a_no_gross_error": float(db.max()) <= L,
        # p90 而不是 p99：对照组 / 坏操作组只有几十点时 p99 就是最大值，重尾下会随机翻；
        #   尾部由 a（公式级亏损线）与 c（< pass_cos 的比例）把关。
        "b_p90_like_control": float(np.percentile(db, 90)) <= 2.0 * float(np.percentile(dg, 90)) + 1e-4,
        "c_rate_like_control": float((db > thr).mean()) <= 2.0 * float((dg > thr).mean()) + 0.002,
        # 原式在"好"操作上也不许有公式级错误 —— 否则逐操作判据漏判了坏操作
        "d_orig_errors_attributed": float(do_g.max()) <= L,
        "e_discriminating": float((do_b > L).mean()) >= DISCRIM_MIN_FRAC,
    }
    ok = all(checks.values())
    note = ("修正式在坏操作点上与对照组（AMSET 原生路径在好操作点）同一水平" if ok else
            "未通过：%s" % ", ".join(k for k, v in checks.items() if not v))
    if not checks["e_discriminating"]:
        ok, note = None, "原式在坏操作点上几乎不出错 —— 无区分力（换结构或用 --shift）"
    if float(dg.max()) > 1.0 - GROSS_COS:
        note += "；注意：对照组地板偏高（最大亏损 %.1e），数据本身的对称一致性差" % float(dg.max())
    return {"ok": ok, "checks": {k: bool(v) for k, v in checks.items()},
            "gross_deficit_line": L, "control_max_deficit": float(dg.max()), "note": note}


def shift_origin(coeffs, gpoints, kpoints, structure, t, chunk=256):
    """原点整体平移 t（分数坐标）：原子 x -> x + t；系数 c_k(G) -> c_k(G)·exp(-2πi (k+G)·t)。

    ψ'(r) = ψ(r - t) 的精确变换，不引入任何误差。**就地**改 coeffs（省一份内存），返回新结构。
    """
    t = np.asarray(t, float)
    kpoints = np.asarray(kpoints, float)
    gpoints = np.asarray(gpoints, float)
    for spin, c in coeffs.items():
        for s in range(0, len(kpoints), chunk):
            sl = slice(s, min(len(kpoints), s + chunk))
            ph = np.exp(-2j * np.pi * ((kpoints[sl, None, :] + gpoints[None, :, :]) @ t))
            if c.ndim == 4:                                  # SOC 旋量 (nb, nk, ng, 2)
                c[:, sl] *= ph[None, :, :, None]
            else:
                c[:, sl] *= ph[None]
    st = structure.copy()
    st.translate_sites(list(range(len(st))), t, frac_coords=True, to_unit_cell=True)
    return st


def _flat(v):
    """(ng,) 或 (ng,2)（SOC 旋量）-> 一维。"""
    return np.asarray(v).reshape(-1)


def subspace_scores(ref, test, groups=None):
    """ref/test: (nb, ng[,2]) 同一 k 点。返回每带分数：
    非简并 = |<ref|test>|/(|ref||test|)；简并组内 = ||P_span(ref 组) test|| / ||test||。"""
    nb = ref.shape[0]
    if groups is None:
        groups = [[i] for i in range(nb)]
    out = np.zeros(nb)
    for g in groups:
        A = np.stack([_flat(ref[i]) for i in g], axis=1)          # (N, d)
        Q, _ = np.linalg.qr(A)
        for i in g:
            b = _flat(test[i])
            nb_ = np.linalg.norm(b)
            out[i] = 0.0 if nb_ == 0 else np.linalg.norm(Q.conj().T @ b) / nb_
    return out


# --------------------------------------------------------------------------
# 去对称化（分块，控制内存）
# --------------------------------------------------------------------------
def _get_formulas():
    import amset.wavefunction.common as _wc
    import amset_desym_fix as _fix
    orig = _wc.desymmetrize_coefficients
    if orig is _fix.desymmetrize_coefficients:
        # 进程里已经被打过补丁（AZ_DESYM_FIX=1）—— 从 AMSET 源码重新拿原函数
        import importlib
        orig = importlib.reload(_wc).desymmetrize_coefficients
    import inspect
    if _fix._squash(_fix.ORIG_FACTOR_LINE) not in _fix._squash(inspect.getsource(orig)):
        raise RuntimeError("拿到的'原式'函数里没有原始 factor 行 —— 不是 AMSET 原函数")
    return {"orig": orig, "fixed": _fix.desymmetrize_coefficients}


def desym_chunks(fn, coeffs_ibz, gpoints, k_ibz, structure, rots, taus, is_tr,
                 op_mapping, kp_mapping, chunk):
    """按目标 k 分块调用 AMSET 的 desymmetrize_coefficients（它的中间数组是
    (nb, nk, g_box) 的稠密盒子，一次展开全部 k 会占满内存）。"""
    n = len(op_mapping)
    for s in range(0, n, chunk):
        sl = slice(s, min(n, s + chunk))
        yield sl, fn(coeffs_ibz, gpoints, k_ibz, structure, rots, taus, is_tr,
                     op_mapping[sl], kp_mapping[sl], pbar=False)


def _auto_chunk(coeffs, gpoints, budget_gb=2.0):
    g_mesh = (np.abs(gpoints).max(axis=0) + 3) * 2
    nb = max(v.shape[0] for v in coeffs.values())
    ncl = 2 if len(next(iter(coeffs.values())).shape) == 4 else 1
    per_k = nb * float(np.prod(g_mesh)) * ncl * 16.0
    return int(max(1, min(512, budget_gb * 1024 ** 3 // per_k)))


def validate(coeffs_full, gpoints, kpoints_full, structure, energies=None,
             mesh=None, symprec=AMSET_SYMPREC, bands=None, pass_cos=PASS_COS,
             chunk=None, formulas=None, verbose=True):
    """核心检验。coeffs_full: {spin: (nb, nk_full, ng[,2])}；energies: {spin: (nb, nk_full)} 或 None。

    返回报告 dict（见 main 的打印）。
    """
    from amset.electronic_structure.symmetry import expand_kpoints
    import amset_desym_fix as _fix

    t0 = time.time()
    kpoints_full = np.asarray(kpoints_full, float)
    mesh = list(mesh) if mesh else mesh_from_kpoints(kpoints_full)
    if int(np.prod(mesh)) != len(kpoints_full):
        raise ValueError("h5 的 k 点数 %d != 网格 %s 的点数 %d —— 不是完整网格（需要 ISYM=-1 的 S3b/S4b）"
                         % (len(kpoints_full), mesh, int(np.prod(mesh))))
    formulas = formulas or _get_formulas()

    k_ibz = ibz_from_mesh(structure, mesh, symprec=symprec, time_reversal=True)
    ibz_idx = match_kpoints(k_ibz, kpoints_full)
    k_ibz = kpoints_full[ibz_idx]                     # 用 h5 自己的坐标（避免 ±0.5 边界歧义）
    full_k, rots, taus, is_tr, op_map, kp_map = expand_kpoints(
        structure, k_ibz, symprec=symprec, return_mapping=True, time_reversal=True,
        verbose=False)
    tgt_idx = match_kpoints(full_k, kpoints_full)
    if len(np.unique(tgt_idx)) != len(kpoints_full):
        raise ValueError("expand_kpoints 展开后覆盖 %d/%d 个全网格点"
                         % (len(np.unique(tgt_idx)), len(kpoints_full)))

    audit = _fix.op_audit(rots, taus, is_tr, op_map)
    exact = np.array([_fix._op_exact(rots[i], taus[i], bool(is_tr[i])) for i in range(len(rots))])
    pred_bad_k = ~exact[op_map]                        # 按逐操作判据，原式会错的目标 k
    # V116 对照组：目标点 = 源点本身（IBZ 代表点映到自己）的是平凡映射，不进统计；
    #   其余按操作分"好"（原式精确 = AMSET 原生路径）/"坏"，再细分 TR / 非 TR。
    trivial = (tgt_idx == ibz_idx[np.asarray(kp_map, int)])
    tr_t = np.asarray(is_tr, bool)[op_map]
    classes = {"good": ~trivial & ~pred_bad_k, "bad": ~trivial & pred_bad_k}
    classes.update({"good_tr": classes["good"] & tr_t, "good_nontr": classes["good"] & ~tr_t,
                    "bad_tr": classes["bad"] & tr_t, "bad_nontr": classes["bad"] & ~tr_t})

    spins = list(coeffs_full.keys())
    report = {"mesh": [int(x) for x in mesh], "nk_full": int(len(kpoints_full)),
              "nk_ibz": int(len(k_ibz)), "symprec": symprec, "degen_tol_eV": DEGEN_TOL_EV,
              "pass_cos": pass_cos, "ops": audit,
              "bad_op_ids": [int(i) for i in np.where(~exact)[0]],
              "bad_ops_tr": int(np.sum(~exact & np.asarray(is_tr, bool))),
              "formulas": {}}
    if verbose:
        print("[..] 网格 %s：全网格 %d 点，IBZ %d 点；AMSET 操作 %d 个（其中原式会错 %d 个，TR %d 个）"
              % ("x".join(map(str, mesh)), len(kpoints_full), len(k_ibz),
                 audit["total_ops"], audit["bad_ops"], report["bad_ops_tr"]))

    worst_by_grid = {}
    for name, fn in formulas.items():
        worst_k = np.ones(len(full_k))                 # 每个目标 k 上各带最小分数
        worst_gap = np.full(len(full_k), np.nan)       # 该点最差带到最近非简并带的能距（meV）
        per_band_min = None
        for spin in spins:
            cf = coeffs_full[spin]
            sel = np.arange(cf.shape[0]) if bands is None else np.asarray(bands, int)
            e_spin = None if energies is None else np.asarray(energies[spin])
            cf = cf[sel]
            c_ibz = {spin: cf[:, ibz_idx]}
            ch = chunk or _auto_chunk(c_ibz, gpoints)
            pbm = np.ones(len(sel))
            for sl, rc in desym_chunks(fn, c_ibz, gpoints, k_ibz, structure, rots, taus,
                                       is_tr, op_map, kp_map, ch):
                rc = rc[spin]
                for j, t in enumerate(range(sl.start, sl.stop)):
                    kf = tgt_idx[t]
                    groups = None
                    if e_spin is not None:
                        groups = degenerate_groups(e_spin[sel, kf])
                    s = subspace_scores(cf[:, kf], rc[:, j], groups)
                    jb = int(np.argmin(s))
                    if s[jb] < worst_k[t]:
                        worst_k[t] = float(s[jb])
                        if e_spin is not None:
                            worst_gap[t] = band_gaps_mev(e_spin[:, kf], [sel[jb]])[0]
                    pbm = np.minimum(pbm, s)
            per_band_min = pbm if per_band_min is None else np.minimum(per_band_min, pbm)
        flagged = worst_k < pass_cos
        r = {"min_score": float(worst_k.min()),
             "per_band_min": [round(float(x), 6) for x in per_band_min],
             "n_bad_k": int(flagged.sum()),
             "frac_bad_k": float(flagged.mean()),
             "passed": bool(worst_k.min() >= pass_cos),
             "by_class": {c: _stats(worst_k[m], pass_cos) for c, m in classes.items()}}
        if energies is not None and flagged.any():
            gf = worst_gap[flagged]
            gf = gf[np.isfinite(gf)]
            ga = worst_gap[~trivial]
            ga = ga[np.isfinite(ga)]
            r["flagged_gap_meV"] = {
                "median": float(np.median(gf)) if len(gf) else None,
                "frac_near_degenerate": float((gf < NEAR_DEGEN_MEV).mean()) if len(gf) else None,
                # 基线：全部非平凡目标点里最差带是近简并的比例 —— 出错点若远高于它，地板来自近简并混合
                "baseline_frac_near_degenerate": float((ga < NEAR_DEGEN_MEV).mean()) if len(ga) else None,
                "threshold_meV": NEAR_DEGEN_MEV}
        if name == "orig":
            gross = worst_k < GROSS_COS
            r.update({
                "pred_bad_k": int(pred_bad_k.sum()),
                "bad_and_predicted": int((flagged & pred_bad_k).sum()),
                "bad_but_not_predicted": int((flagged & ~pred_bad_k).sum()),   # 严格判据要求 0
                "predicted_but_ok": int((~flagged & pred_bad_k).sum()),        # 允许 >0（误差恰不可见）
                "bad_k_ops": sorted({int(op_map[t]) for t in np.where(flagged)[0]}),
                # V116：公式级（< GROSS_COS）出错点的归因 —— 数值地板不参与
                "gross_bad_k": int(gross.sum()),
                "gross_bad_but_not_predicted": int((gross & ~pred_bad_k).sum()),
            })
            r["attribution_ok"] = (r["bad_but_not_predicted"] == 0
                                   and set(r["bad_k_ops"]) <= set(report["bad_op_ids"]))
        report["formulas"][name] = r
        wg = np.ones(len(kpoints_full))
        wg[tgt_idx] = worst_k
        worst_by_grid[name] = wg
        if verbose:
            print("[..] %-5s 最小分数 %.6f，出错 k 点 %d/%d（%.1f%%）%s"
                  % (name, r["min_score"], r["n_bad_k"], len(full_k), 100 * r["frac_bad_k"],
                     ("；归因：预测坏 %d，出错且预测坏 %d，出错但未预测 %d，预测坏但未出错 %d；"
                      "公式级(<%.2f) 出错 %d，其中未预测 %d"
                      % (r["pred_bad_k"], r["bad_and_predicted"], r["bad_but_not_predicted"],
                         r["predicted_but_ok"], GROSS_COS, r["gross_bad_k"],
                         r["gross_bad_but_not_predicted"])) if name == "orig" else ""))
            for c in ("good", "bad"):
                st_ = r["by_class"][c]
                if st_["n"]:
                    print("        %-4s 操作点 %5d：最小 %.6f  p99 亏损 %.2e  <%.3f 比例 %.2f%%"
                          % (c, st_["n"], st_["min"], st_["p99_deficit"], pass_cos,
                             100 * st_["frac_below_pass"]))
    fx = report["formulas"].get("fixed", {})
    og = report["formulas"].get("orig", {})
    report["verdict_strict"] = bool(fx.get("passed", True) and og.get("attribution_ok", True))
    if fx and og:
        wf, wo = worst_by_grid["fixed"][tgt_idx], worst_by_grid["orig"][tgt_idx]
        report["control"] = control_checks(wf, wo, classes, pass_cos)
        # 同一套判据套在原式上：必须**不**通过，才说明判据能区分对错
        report["control"]["orig_would_pass"] = control_checks(wo, wo, classes, pass_cos).get("ok")
    else:
        report["control"] = {"ok": None, "note": "只跑了一个公式"}
    report["verdict_control"] = bool(report["control"].get("ok"))
    report["verdict"] = bool(report["verdict_strict"] or report["verdict_control"])
    report["seconds"] = round(time.time() - t0, 1)
    report["_worst_by_grid"] = worst_by_grid          # 非 JSON：平移不变性比较用，main 写盘前去掉
    return report


# --------------------------------------------------------------------------
# 读真实数据
# --------------------------------------------------------------------------
def _energies_from_vasprun(vasprun, kpoints_h5, nb_h5, energy_cutoff=None, bands=None):
    """按 amset wave 的同一套规则（get_band_structure + get_ibands）从 vasprun 取出 h5 里那些带的能量。"""
    from pymatgen.io.vasp import BSVasprun
    from amset.constants import defaults
    from amset.electronic_structure.common import get_band_structure, get_ibands
    vr = BSVasprun(str(vasprun))
    bs = get_band_structure(vr, zero_weighted=defaults["zero_weighted_kpoints"])
    if bands:
        from amset.util import parse_ibands
        ib = parse_ibands(bands)
    else:
        ib = get_ibands(energy_cutoff or defaults["energy_cutoff"], bs)
    kbs = np.array([k.frac_coords for k in bs.kpoints])
    idx = match_kpoints(kpoints_h5, kbs)
    out = {}
    for spin, b in bs.bands.items():
        sel = np.asarray(ib[spin], int)
        if len(sel) != nb_h5[spin]:
            raise ValueError("vasprun 按 energy_cutoff 取到 %d 条带，h5 里是 %d 条 —— "
                             "请用 --energy-cutoff / --amset-bands 对齐 amset wave 当时的取法"
                             % (len(sel), nb_h5[spin]))
        out[spin] = np.asarray(b)[sel][:, idx]
    return out


def validate_with_shift(coeffs, gpoints, kpoints, structure, t, **kw):
    """V116：先在原原点上跑修正式（参照），再把原点平移 t 后两式都跑，并比较修正式逐点分数。

    **就地**改 coeffs（平移相位）。返回平移后的报告，附 report["shift"]。
    """
    kw = dict(kw)
    kw.pop("formulas", None)
    verbose = kw.pop("verbose", True)
    formulas = _get_formulas()
    rep0 = validate(coeffs, gpoints, kpoints, structure, formulas={"fixed": formulas["fixed"]},
                    verbose=False, **kw)
    st_s = shift_origin(coeffs, gpoints, kpoints, structure, t)
    if verbose:
        print("[..] 原点平移 t = %s：平移前坏操作 %d/%d" % (list(map(float, t)),
              rep0["ops"]["bad_ops"], rep0["ops"]["total_ops"]))
    rep = validate(coeffs, gpoints, kpoints, st_s, formulas=formulas, verbose=verbose, **kw)
    w0, w1 = rep0["_worst_by_grid"]["fixed"], rep["_worst_by_grid"]
    d_fix = float(np.abs(w1["fixed"] - w0).max())
    rep["shift"] = {"t": [float(x) for x in t],
                    "bad_ops_before": int(rep0["ops"]["bad_ops"]),
                    "bad_ops_after": int(rep["ops"]["bad_ops"]),
                    "fixed_max_abs_diff": d_fix,
                    "orig_max_abs_diff": float(np.abs(w1["orig"] - w0).max()),
                    "invariant": bool(d_fix < SHIFT_TOL),
                    "shifted_structure_first_site": [round(float(x), 4) for x in st_s.frac_coords[0]]}
    rep["verdict"] = bool(rep["verdict"] and rep["shift"]["invariant"])
    rep["_structure"] = st_s
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--h5", default="step4b_wave_full/wavefunction.h5")
    ap.add_argument("--vasprun", default=None,
                    help="同一全网格计算的 vasprun.xml（用于简并判断；缺省则逐带比较并告警）")
    ap.add_argument("--mesh", nargs=3, type=int, default=None)
    ap.add_argument("--bands", default=None, help="只比 h5 里的这些带（h5 内下标，如 0:8）")
    ap.add_argument("--energy-cutoff", type=float, default=None)
    ap.add_argument("--amset-bands", default=None, help="amset wave 当时用的 -b（覆盖 energy-cutoff）")
    ap.add_argument("--symprec", type=float, default=AMSET_SYMPREC)
    ap.add_argument("--pass-cos", type=float, default=PASS_COS)
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--shift", nargs=3, type=float, default=None,
                    help="原点平移 t（分数坐标），见文件头；2D 第三个分量写 0")
    ap.add_argument("--json", default="desym_fix_validation.json")
    a = ap.parse_args(argv)

    from amset.wavefunction.io import load_coefficients
    h5 = Path(a.h5)
    if not h5.is_file():
        print("[ERROR] 找不到 %s" % h5, file=sys.stderr)
        return 2
    coeffs, gpoints, kpoints, st = load_coefficients(str(h5))
    if gpoints is None:
        print("[ERROR] h5 里没有 gpoints（太老的 amset wave？）", file=sys.stderr)
        return 2
    energies = None
    if a.vasprun:
        energies = _energies_from_vasprun(a.vasprun, kpoints,
                                          {s: v.shape[0] for s, v in coeffs.items()},
                                          a.energy_cutoff, a.amset_bands)
    else:
        print("[WARN] 没给 --vasprun：不做简并子空间比较，简并带会被误报为出错", file=sys.stderr)
    bands = None
    if a.bands:
        lo, _, hi = a.bands.partition(":")
        bands = list(range(int(lo), int(hi))) if hi else [int(lo)]
    kw = dict(energies=energies, mesh=a.mesh, symprec=a.symprec, bands=bands,
              pass_cos=a.pass_cos, chunk=a.chunk)
    if a.shift:
        rep = validate_with_shift(coeffs, gpoints, kpoints, st, a.shift, **kw)
        st = rep.pop("_structure")
    else:
        rep = validate(coeffs, gpoints, kpoints, st, **kw)
    rep.pop("_worst_by_grid", None)
    try:
        import ke_common as _kc
        g = _kc.symmetry_gate(st)
        print("[..] ke_common.symmetry_gate(%s结构)：bad_ops %s/%s（%s）"
              % ("平移后" if a.shift else "h5 ", g.get("bad_ops"), g.get("total_ops"), g.get("reason")))
    except Exception as e:                                   # noqa: BLE001
        g = None
        print("[..] ke_common 不可用（%s），只用 AMSET 自己的操作集" % e)
    rep["h5"] = str(h5)
    rep["first_site"] = [round(float(x), 4) for x in st.frac_coords[0]]
    if g is not None:
        rep["ke_common_gate"] = {k: g.get(k) for k in ("bad_ops", "total_ops", "needs_full_grid")}
    Path(a.json).write_text(json.dumps(rep, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fx, og = rep["formulas"]["fixed"], rep["formulas"]["orig"]
    ctl = rep["control"]
    print("\n==== 结论 ====")
    print("修正式：最小分数 %.6f  -> 严格判据 %s（要求 >= %.3f）"
          % (fx["min_score"], "PASS" if fx["passed"] else "FAIL", a.pass_cos))
    print("原  式：出错 k 点 %d/%d = %.1f%%；出错点全部来自坏操作：%s（坏操作 %d/%d）；"
          "公式级(<%.2f) 出错 %d 点，其中非坏操作 %d 点"
          % (og["n_bad_k"], rep["nk_full"], 100 * og["frac_bad_k"],
             "是" if og["attribution_ok"] else "否", rep["ops"]["bad_ops"], rep["ops"]["total_ops"],
             GROSS_COS, og["gross_bad_k"], og["gross_bad_but_not_predicted"]))
    print("对照组判据（V116）：%s —— %s"
          % ({True: "PASS", False: "FAIL", None: "无从判定"}[ctl.get("ok")], ctl.get("note")))
    if "checks" in ctl:
        print("    %s；同一判据套原式：%s"
              % ("  ".join("%s=%s" % (k, "✓" if v else "✗") for k, v in ctl["checks"].items()),
                 {True: "通过（判据无区分力！）", False: "不通过（应当如此）", None: "无从判定"}
                 [ctl.get("orig_would_pass")]))
    if "flagged_gap_meV" in fx:
        fg = fx["flagged_gap_meV"]
        if fg.get("frac_near_degenerate") is not None:
            print("修正式出错点的近简并诊断：最差带到最近非简并带能距中位 %.1f meV；"
                  "< %.0f meV 的占 %.0f%%（全部非平凡点基线 %.0f%%）"
                  % (fg["median"], fg["threshold_meV"], 100 * fg["frac_near_degenerate"],
                     100 * (fg["baseline_frac_near_degenerate"] or 0)))
    if a.shift:
        sh = rep["shift"]
        print("平移不变性：修正式逐点 |Δ| 最大 %.2e -> %s（要求 < %.0e）；原式 |Δ| 最大 %.2e；"
              "坏操作 %d -> %d"
              % (sh["fixed_max_abs_diff"], "PASS" if sh["invariant"] else "FAIL", SHIFT_TOL,
                 sh["orig_max_abs_diff"], sh["bad_ops_before"], sh["bad_ops_after"]))
    print("总判定：%s（严格 %s / 对照组 %s）"
          % ("PASS" if rep["verdict"] else "FAIL", rep["verdict_strict"], rep["verdict_control"]))
    print("写出 %s（%.1f s）" % (a.json, rep["seconds"]))
    return 0 if rep["verdict"] else 1


if __name__ == "__main__":
    sys.exit(main())
