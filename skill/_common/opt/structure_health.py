#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""structure_health.py —— 结构体检四判据（S1 gen 调用；默认只告警）。

wangchao 2026-09-21 列的四条（题干来自手工拦下 AlN / Zn5O3 / BeO 的那一轮）：

  ① 最近邻距离 < 0.7 × 共价半径和 —— 抓"两个原子坐在同一个位"的种子错误
     （MoS2 seed：Mo、S 同面内位 -> d(Mo-S)=1.57 Å，DFT 直接把胞缩了 14%）。
  ② CN=1 悬挂 —— 抓欠配位/断键（Al5N 那种 CN 1~7 的怪结构）。
  ③ 命名 vs 化学计量比 —— 抓"目录名说 AlN、结构其实是别的"。
  ④ z 跨度 vs 层数 —— 抓"声称单层、其实是块体/过厚"（AlN：span 4.34 Å 却当单层）。

第五条（symprec 1e-5 vs 1e-4 双容差空间群）在 symmetry_audit.py，本模块不管。

设计与 symmetry_audit.py 一致：顶层只标准库；numpy 懒加载；**任何异常只跳过、绝不阻断 gen**。
半径取 Cordero et al. (2008) 共价半径（Å），故意用纯字典、不引 pymatgen —— gen 在登录节点
用系统 python 跑，公共池不引入重依赖。
"""
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

__all__ = [
    "COVALENT_RADII", "check_structure", "check_structure_file",
    "format_report", "reduced_formula",
]

# ---------------------------------------------------------------------------
# 判据阈值（可被 check_structure 的关键字覆盖）
# ---------------------------------------------------------------------------
OVERLAP_RATIO   = 0.7     # ① d < ratio×(rcov_i+rcov_j) 判重叠
ABS_MIN_DIST    = 0.50    # ① 绝对下限：任何原子对 d < 0.5 Å 都荒谬
CN_CUTOFF_SCALE = 1.0     # ② CN 判据：d < 1.0×(rcov_i+rcov_j)（=共价半径和）记为成键
                          #   标定：AlN[0-3,d8] / Zn5O3[0-1,d32] / BeO[0-3,d4] 都见悬挂，
                          #   而四个候选 + 两套 MoS₂ 全 d0（见 tests/suite_structure_health.py）
LAYER_GAP_A     = 0.8     # ④ z 方向相邻原子 > 此值算层间断口
SPAN_MAX_2D_A   = 6.0     # ④ 2D 层起伏超过此值 -> "跨度过大"
MONO_THICK_MAX  = 3.0     # ④ 名义单层但起伏超过此值 -> "单层过厚"
VACUUM_MIN_A    = 8.0     # 真空 >= 此值按 2D/slab 处理 span 判据

# Cordero et al., Dalton Trans. (2008) 共价半径（Å）。缺元素退化为 1.4。
COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01,
    "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92,
    "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.40,
    "At": 1.50, "Rn": 1.50, "Fr": 2.60, "Ra": 2.21, "Ac": 2.15, "Th": 2.06,
    "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87, "Am": 1.80,
}
COV_FALLBACK = 1.40


def cov_radius(symbol):
    return COVALENT_RADII.get(symbol, COV_FALLBACK)


# ---------------------------------------------------------------------------
# 结构读取（复用 symmetry_audit 的解析，少一份 POSCAR 解析器）
# ---------------------------------------------------------------------------
def _read_structure(path):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import symmetry_audit as SA
    st = SA.read_poscar(path)
    return st


def _import_numpy():
    try:
        import numpy as np
        return np
    except Exception:                              # noqa: BLE001
        return None


def _reduce_counts(counts):
    g = 0
    for v in counts.values():
        g = math.gcd(g, int(v))
    g = g or 1
    return {k: int(v) // g for k, v in counts.items()}


def reduced_formula(symbols):
    """按元素首次出现顺序给最简化学式（Mo4S6 -> Mo2S3）。"""
    order = []
    for s in symbols:
        if s not in order:
            order.append(s)
    counts = _reduce_counts(Counter(symbols))
    return "".join("%s%s" % (s, "" if counts[s] == 1 else counts[s]) for s in order)


_FORMULA_RE = re.compile(r"^([A-Z][a-z]?\d*)+$")


def _looks_like_formula(s):
    return bool(s) and bool(_FORMULA_RE.match(s))


def _formula_elements_and_ratio(s):
    """把 'Mo2S3' 拆成 ({Mo,S}, {Mo:2,S:3})；非化学式返回 (None, None)。"""
    if not _looks_like_formula(s):
        return None, None
    elems, counts = [], {}
    for m in re.finditer(r"([A-Z][a-z]?)(\d*)", s):
        sym, num = m.group(1), m.group(2)
        if not sym:
            continue
        elems.append(sym)
        counts[sym] = counts.get(sym, 0) + (int(num) if num else 1)
    if not elems:
        return None, None
    return set(elems), counts


def _norm(s):
    return re.sub(r"[^A-Za-z]", "", s or "")


# ---------------------------------------------------------------------------
# ① 最近邻 / ② 配位数
# ---------------------------------------------------------------------------
def min_image_distances(lattice, frac):
    """(N,N) 最小镜像距离矩阵。lattice=3x3(行=晶格矢量)，frac=Nx3。"""
    np = _import_numpy()
    L = np.asarray(lattice, dtype=float)
    F = np.asarray(frac, dtype=float)
    d = F[:, None, :] - F[None, :, :]
    d -= np.round(d)                      # 最小镜像（对一般胞是标准近似）
    cart = d @ L
    return np.sqrt((cart ** 2).sum(axis=-1))


def nearest_pair(dist, symbols):
    """全局最近的一对（忽略自配对）。返回 dict。"""
    np = _import_numpy()
    n = len(symbols)
    best = None
    for i in range(n):
        for j in range(i + 1, n):
            dij = float(dist[i, j])
            thr = OVERLAP_RATIO * (cov_radius(symbols[i]) + cov_radius(symbols[j]))
            ratio = dij / max(1e-9, cov_radius(symbols[i]) + cov_radius(symbols[j]))
            if best is None or dij < best["dist_A"]:
                best = {"i": i, "j": j, "si": symbols[i], "sj": symbols[j],
                        "dist_A": dij, "overlap_thr_A": thr, "dist_over_sum": ratio}
    if best is None:
        best = {"i": None, "j": None, "si": None, "sj": None,
                "dist_A": None, "overlap_thr_A": None, "dist_over_sum": None}
    ok = (best["dist_A"] is not None
          and best["dist_A"] >= ABS_MIN_DIST
          and best["dist_A"] >= best["overlap_thr_A"])
    best["ok"] = bool(ok)
    return best


def coordination_numbers(lattice, frac, symbols, scale=None):
    """周期镜像下的配位数：CN_i = #{(j,T): d(i, j+T) < scale×(rcov_i+rcov_j)}。

    ★ 必须算周期镜像：原胞只有少数原子时（如 MoS₂ 3 原子胞），一个 S 的 3 个 Mo 邻居
    全是同一个 Mo 原子的 3 个周期像；只数胞内会把它错算成 CN=1、误报悬挂。
    scale 默认 1.0（=共价半径和）。返回 (cn list, dangling list)。
    """
    np = _import_numpy()
    scale = CN_CUTOFF_SCALE if scale is None else float(scale)
    L = np.asarray(lattice, dtype=float)
    F = np.asarray(frac, dtype=float)
    C = F @ L
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)],
                      dtype=float)
    images = shifts @ L
    P = C[:, None, :] + images[None, :, :]                       # (N,27,3)
    D = np.linalg.norm(C[:, None, None, :] - P[None, :, :, :], axis=-1)   # (N,N,27)
    rs = np.array([cov_radius(s) for s in symbols], dtype=float)
    cut = scale * (rs[:, None] + rs[None, :])
    mask = D < cut[:, :, None]
    idx0 = int(np.argmin(np.abs(shifts).sum(axis=1)))            # (0,0,0) 的索引
    for i in range(len(symbols)):
        mask[i, i, idx0] = False
    cn = [int(x) for x in mask.sum(axis=(1, 2))]
    dangling = [i for i, c in enumerate(cn) if c <= 1]
    return cn, dangling


# ---------------------------------------------------------------------------
# ③ 命名 vs 计量比
# ---------------------------------------------------------------------------
def stoichiometry_check(symbols, material_name):
    res = {"material_name": material_name, "claim": None, "formula": reduced_formula(symbols),
           "elements_match": None, "ok": True, "note": ""}
    if not material_name:
        res["note"] = "未提供材料名，跳过命名核对"
        return res
    claim = material_name.rstrip("/").split("_")[-1]
    res["claim"] = claim
    if not _looks_like_formula(claim):
        res["note"] = "材料名末段 %r 不是化学式，跳过" % claim
        return res
    claim_elems, claim_counts = _formula_elements_and_ratio(claim)
    actual_elems = set(symbols)
    actual_counts = _reduce_counts(Counter(symbols))
    res["elements_match"] = (claim_elems == actual_elems)
    if claim_elems != actual_elems:
        res["ok"] = False
        res["note"] = ("命名元素 %s != 结构元素 %s（实际最简式 %s）"
                       % (sorted(claim_elems), sorted(actual_elems), res["formula"]))
        return res
    # 元素一致：再比最简比（元素顺序无关）
    claim_red = _reduce_counts(claim_counts)
    if claim_red != actual_counts:
        res["ok"] = False
        res["note"] = ("命名计量比 %s != 结构 %s" % (claim_red, actual_counts))
    else:
        res["note"] = "命名与结构一致（%s）" % res["formula"]
    return res


# ---------------------------------------------------------------------------
# ④ z 跨度 / 层数
# ---------------------------------------------------------------------------
def layer_analysis(lattice, frac, symbols, vac_axis, vacuum_hint=None):
    np = _import_numpy()
    L = np.asarray(lattice, dtype=float)
    F = np.asarray(frac, dtype=float)
    cart = F @ L
    z = cart[:, int(vac_axis)]
    zs = np.sort(z)
    span = float(zs.max() - zs.min())
    gaps = np.diff(zs)
    n_layers = 1 + int((gaps > LAYER_GAP_A).sum())
    axlen = float(np.linalg.norm(L[int(vac_axis)]))
    # 真空优先用 vacuum_info 的严格口径（最大真空隙 × 垂直胞高），否则退化为 |a| - span
    vacuum = float(vacuum_hint) if vacuum_hint is not None else (axlen - span)
    flags = []
    is_2d = vacuum >= VACUUM_MIN_A
    if is_2d:
        if span > SPAN_MAX_2D_A:
            flags.append("2D 层起伏 %.2f Å > %.1f（跨度过大，可能不是单层）" % (span, SPAN_MAX_2D_A))
        if n_layers <= 1 and span > MONO_THICK_MAX:
            flags.append("名义单层但起伏 %.2f Å > %.1f（单层过厚）" % (span, MONO_THICK_MAX))
    return {"axis": int(vac_axis), "span_A": span, "n_layers": n_layers,
            "vacuum_A": vacuum, "is_2d": bool(is_2d), "flags": flags}


# ---------------------------------------------------------------------------
# 顶层
# ---------------------------------------------------------------------------
def check_structure(lattice, frac, symbols, material_name=None, vacuum_min=VACUUM_MIN_A):
    """四判据体检。返回 dict（available/ok/flags/各判据明细）。"""
    np = _import_numpy()
    if np is None:
        return {"available": False, "reason": "numpy 不可用", "ok": True, "flags": []}
    dist = min_image_distances(lattice, frac)
    near = nearest_pair(dist, symbols)
    cn, dangling = coordination_numbers(lattice, frac, symbols)
    stoich = stoichiometry_check(symbols, material_name)
    # 真空轴：与 symmetry_audit.vacuum_info 同一口径（最大真空隙所在轴），不再用"最长基矢"
    ax, vac = _vacuum_axis(lattice, frac)
    layers = layer_analysis(lattice, frac, symbols, ax, vacuum_hint=vac)

    flags = []
    if not near["ok"]:
        flags.append("① 最近邻过短：%s%d-%s%d d=%.3f Å（0.7×(rcov和)=%.3f Å，比值 %.2f）"
                     % (near["si"], near["i"], near["sj"], near["j"],
                        near["dist_A"], near["overlap_thr_A"], near["dist_over_sum"]))
    if dangling:
        flags.append("② CN≤1 悬挂 %d 个（原子 %s；CN 范围 %d-%d）"
                     % (len(dangling), dangling[:8], min(cn), max(cn)))
    if not stoich["ok"]:
        flags.append("③ %s" % stoich["note"])
    flags.extend("④ %s" % f for f in layers["flags"])
    return {
        "available": True,
        "material_name": material_name,
        "n_atoms": len(symbols),
        "composition": dict(Counter(symbols)),
        "formula": reduced_formula(symbols),
        "nearest": near,
        "coordination": {"cn": cn, "cn_min": int(min(cn)) if cn else None,
                         "cn_max": int(max(cn)) if cn else None,
                         "n_dangling": len(dangling), "dangling": dangling,
                         "cutoff_scale": CN_CUTOFF_SCALE},
        "stoichiometry": stoich,
        "layers": layers,
        "flags": flags,
        "ok": not flags,
    }


def _vacuum_axis(lattice, frac):
    """面外轴 + 真空厚度。口径同 symmetry_audit.vacuum_info（最大真空隙所在轴）。"""
    try:
        import symmetry_audit as SA
        vi = SA.vacuum_info(lattice, frac)
        return int(vi["axis"]), float(vi["vacuum_A"])
    except Exception:                                          # noqa: BLE001
        np = _import_numpy()
        L = np.asarray(lattice, dtype=float)
        comps = [float(np.linalg.norm(L[i])) for i in range(3)]
        return (int(max(range(3), key=lambda i: comps[i])) if comps else 2), None


def check_structure_file(path, material_name=None, vacuum_min=VACUUM_MIN_A):
    """从 POSCAR 读结构并体检；spglib 无关。失败返回 {available: False}（不抛）。"""
    try:
        st = _read_structure(path)
    except Exception as exc:                       # noqa: BLE001
        return {"available": False, "reason": "POSCAR 解析失败：%s" % exc,
                "path": str(path), "ok": True, "flags": []}
    return check_structure(st["lattice"], st["frac"], st["symbols"],
                           material_name=material_name, vacuum_min=vacuum_min)


def format_report(res, tag="structure_health"):
    if not res.get("available"):
        return "[%s] 跳过：%s" % (tag, res.get("reason"))
    L = []
    L.append("[%s] n_atoms=%d composition=%s formula=%s"
             % (tag, res["n_atoms"], res["composition"], res["formula"]))
    n = res["nearest"]
    L.append("  ① 最近邻 %s%d-%s%d  d=%.4f Å  阈值 %.4f Å  比值 %.3f  %s"
             % (n["si"], n["i"], n["sj"], n["j"], n["dist_A"] or 0.0,
                n["overlap_thr_A"] or 0.0, n["dist_over_sum"] or 0.0,
                "OK" if n["ok"] else "FAIL"))
    c = res["coordination"]
    L.append("  ② CN 范围 %s-%s  悬挂(≤1) %d 个  %s"
             % (c["cn_min"], c["cn_max"], c["n_dangling"],
                "OK" if c["n_dangling"] == 0 else "FAIL"))
    s = res["stoichiometry"]
    L.append("  ③ 命名 %s vs 结构 %s  %s（%s）"
             % (s["claim"], s["formula"], "OK" if s["ok"] else "FAIL", s["note"]))
    y = res["layers"]
    L.append("  ④ z 跨度 %.3f Å  层数 %d  真空 %.2f Å  2D=%s  %s"
             % (y["span_A"], y["n_layers"], y["vacuum_A"], y["is_2d"],
                "OK" if not y["flags"] else "FAIL"))
    if res["flags"]:
        L.append("  [WARN] 危险信号：")
        L.extend("    - " + f for f in res["flags"])
    else:
        L.append("  四判据全部通过")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_parser():
    import argparse
    p = argparse.ArgumentParser(description="结构体检四判据（重叠/悬挂/计量比/层厚）")
    p.add_argument("--poscar", required=True, help="输入 POSCAR/CONTCAR")
    p.add_argument("--name", default=None, help="材料目录名（用于命名 vs 计量比）")
    p.add_argument("--json", action="store_true", help="stdout 输出 JSON")
    p.add_argument("--json-out", default=None, help="把结果 JSON 写到该路径")
    p.add_argument("--strict-exit", action="store_true",
                   help="有危险信号时退出码 1（默认 0，只告警）")
    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)
    res = check_structure_file(args.poscar, material_name=args.name)
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(res, ensure_ascii=False, indent=2),
                                       encoding="utf-8", newline="\n")
    if args.json:
        print(json.dumps(res, ensure_ascii=False, indent=2))
    else:
        print(format_report(res))
    if args.strict_exit and res.get("available") and not res.get("ok"):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
