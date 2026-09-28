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

判据：
  · 修正式：全部 k 点、全部带的最小值 >= 0.999（--pass-cos）才算通过；
  · 原式：报出错 k 点（最小值 < 0.999）比例，并逐点核对是否**正好**是被"坏操作"
    （逐操作判据 (R±I)τ≢0，用 AMSET 自己取到的操作集算）映射过去的点。

用法（在材料目录，AMSET 环境里；不需要提交作业，~5–10 分钟）：
    python <skill>/ke-dft-cpu/tools/desym_fix_validate.py \
        --h5 step4b_wave_full/wavefunction.h5 --vasprun step4b_wave_full/vasprun.xml
可选：--mesh 48 48 3（默认从 h5 的 k 点推）、--bands 0:8（只比部分带，省时间）、
      --energy-cutoff（与 amset wave 取带时一致；默认 AMSET 默认值）、--json out.json。
退出码：0 = 修正式通过且原式的出错点全部可归因；1 = 未通过；2 = 输入问题。

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

    for name, fn in formulas.items():
        worst_k = np.ones(len(full_k))                 # 每个目标 k 上各带最小分数
        worst_band = np.ones(1)
        per_band_min = None
        for spin in spins:
            cf = coeffs_full[spin]
            sel = np.arange(cf.shape[0]) if bands is None else np.asarray(bands, int)
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
                    if energies is not None:
                        groups = degenerate_groups(np.asarray(energies[spin])[sel, kf])
                    s = subspace_scores(cf[:, kf], rc[:, j], groups)
                    worst_k[t] = min(worst_k[t], float(s.min()))
                    pbm = np.minimum(pbm, s)
            per_band_min = pbm if per_band_min is None else np.minimum(per_band_min, pbm)
        flagged = worst_k < pass_cos
        r = {"min_score": float(worst_k.min()),
             "per_band_min": [round(float(x), 6) for x in per_band_min],
             "n_bad_k": int(flagged.sum()),
             "frac_bad_k": float(flagged.mean()),
             "passed": bool(worst_k.min() >= pass_cos)}
        if name == "orig":
            r.update({
                "pred_bad_k": int(pred_bad_k.sum()),
                "bad_and_predicted": int((flagged & pred_bad_k).sum()),
                "bad_but_not_predicted": int((flagged & ~pred_bad_k).sum()),   # 必须 0
                "predicted_but_ok": int((~flagged & pred_bad_k).sum()),        # 允许 >0（误差恰不可见）
                "bad_k_ops": sorted({int(op_map[t]) for t in np.where(flagged)[0]}),
            })
            r["attribution_ok"] = (r["bad_but_not_predicted"] == 0
                                   and set(r["bad_k_ops"]) <= set(report["bad_op_ids"]))
        report["formulas"][name] = r
        if verbose:
            print("[..] %-5s 最小分数 %.6f，出错 k 点 %d/%d（%.1f%%）%s"
                  % (name, r["min_score"], r["n_bad_k"], len(full_k), 100 * r["frac_bad_k"],
                     ("；归因：预测坏 %d，出错且预测坏 %d，出错但未预测 %d，预测坏但未出错 %d"
                      % (r["pred_bad_k"], r["bad_and_predicted"], r["bad_but_not_predicted"],
                         r["predicted_but_ok"])) if name == "orig" else ""))
    fx = report["formulas"].get("fixed", {})
    og = report["formulas"].get("orig", {})
    report["verdict"] = bool(fx.get("passed", True) and og.get("attribution_ok", True))
    report["seconds"] = round(time.time() - t0, 1)
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
    try:
        import ke_common as _kc
        g = _kc.symmetry_gate(st)
        print("[..] ke_common.symmetry_gate(h5 结构)：bad_ops %s/%s（%s）"
              % (g.get("bad_ops"), g.get("total_ops"), g.get("reason")))
    except Exception as e:                                   # noqa: BLE001
        g = None
        print("[..] ke_common 不可用（%s），只用 AMSET 自己的操作集" % e)
    rep = validate(coeffs, gpoints, kpoints, st, energies=energies, mesh=a.mesh,
                   symprec=a.symprec, bands=bands, pass_cos=a.pass_cos, chunk=a.chunk)
    rep["h5"] = str(h5)
    rep["first_site"] = [round(float(x), 4) for x in st.frac_coords[0]]
    if g is not None:
        rep["ke_common_gate"] = {k: g.get(k) for k in ("bad_ops", "total_ops", "needs_full_grid")}
    Path(a.json).write_text(json.dumps(rep, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fx, og = rep["formulas"]["fixed"], rep["formulas"]["orig"]
    print("\n==== 结论 ====")
    print("修正式：最小分数 %.6f  -> %s（要求 >= %.3f）"
          % (fx["min_score"], "PASS" if fx["passed"] else "FAIL", a.pass_cos))
    print("原  式：出错 k 点 %d/%d = %.1f%%；出错点全部来自坏操作：%s（坏操作 %d/%d）"
          % (og["n_bad_k"], rep["nk_full"], 100 * og["frac_bad_k"],
             "是" if og["attribution_ok"] else "否", rep["ops"]["bad_ops"], rep["ops"]["total_ops"]))
    print("写出 %s（%.1f s）" % (a.json, rep["seconds"]))
    return 0 if rep["verdict"] else 1


if __name__ == "__main__":
    sys.exit(main())
