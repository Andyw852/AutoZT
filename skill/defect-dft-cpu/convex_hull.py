#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""convex_hull.py —— 第0步：目标化合物的化学势稳定性窗口（凸包）。

纯标准库，登录节点可跑。读各相总能（元素相 + 竞争二元相 + 目标化合物），
算出目标化合物能稳定存在的化学势窗口，产出 chemical_potential_window.json，
并可写 energies.json（形成能脚本的输入）。

输入 phases.json：
  { "target": {"formula": {"Sn":2,"Sb":1,"Te":5}, "E": -xx.xx},   # 目标化合物(每胞)总能
    "elements": {"Sn": -a, "Sb": -b, "Te": -c},                   # 元素相总能(每原子)
    "phases": [ {"formula": {"Sn":1,"Te":1}, "E": -yy}, ... ] }   # 竞争相(每式量)总能

方法（Δμ_i = μ_i − E_i(元素相)，Δμ_i ≤ 0）：
  平衡：      Σ_i n_i(目标)·Δμ_i = ΔH_f(目标)
  竞争相稳定：Σ_i n_i(相)·Δμ_i ≤ ΔH_f(相)
三元体系消去一个变量 → 2D 凸多边形顶点枚举。
"""
import sys, os, json, math
from pathlib import Path

TOL = 1e-7

def load_phases(path="phases.json"):
    if not os.path.exists(path):
        raise SystemExit("[错误] 找不到 %s —— 请先算好各相总能（README §0）" % path)
    return json.load(open(path, encoding="utf-8"))

def formation_energy(formula, E, elements):
    """ΔH_f = E(相) − Σ_i n_i·E_i(元素相每原子)。"""
    if not formula or any(not isinstance(n, (int, float)) or not math.isfinite(n) or n < 0 for n in formula.values()) or sum(formula.values()) <= 0:
        raise SystemExit("[错误] 相成分必须包含有限非负计数且总数为正")
    missing = set(formula) - set(elements)
    if missing:
        raise SystemExit("[错误] 缺少元素参考能: %s" % sorted(missing))
    if not math.isfinite(E) or any(not math.isfinite(elements[el]) for el in formula):
        raise SystemExit("[错误] 相能量/元素参考能必须有限")
    return E - sum(n * elements[el] for el, n in formula.items())

def solve_2x2(a1, b1, c1, a2, b2, c2):
    """解 a1*x+b1*y=c1, a2*x+b2*y=c2；无解返回 None。"""
    det = a1*b2 - a2*b1
    if abs(det) < 1e-12:
        return None
    x = (c1*b2 - c2*b1) / det
    y = (a1*c2 - a2*c1) / det
    return x, y

def _solve_linear(A, b):
    """高斯消元解 A x = b（A 为 m x m）。无唯一解返回 None。"""
    m = len(b)
    M = [list(A[i]) + [b[i]] for i in range(m)]
    for col in range(m):
        piv = None
        for r in range(col, m):
            if abs(M[r][col]) > 1e-12:
                piv = r
                break
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        for j in range(col, m + 1):
            M[col][j] /= pv
        for r in range(m):
            if r != col and abs(M[r][col]) > 1e-15:
                fac = M[r][col]
                for j in range(col, m + 1):
                    M[r][j] -= fac * M[col][j]
    return [M[i][m] for i in range(m)]


def energy_above_hull(target, elements, phases):
    """目标相的标准"凸包以上能量"（每原子，eV/atom，恒 ≥0）：=0 稳定，>0 亚稳。

    在成分空间构造凸包（顶点 = 元素相 ΔH=0 与各竞争相），用单纯形枚举求目标成分处的
    凸包能量。数据不足/退化时返回 None。用于在窗口为空时判断"数据有误"还是"目标亚稳"。
    """
    els = list(elements.keys())
    N = len(els)
    if N == 0:
        return None
    pts = []
    for i in range(N):
        x = [0.0] * N
        x[i] = 1.0
        pts.append((x, 0.0))
    for ph in phases:
        cnt = [float(ph["formula"].get(e, 0)) for e in els]
        tot = sum(cnt)
        if tot <= 0:
            continue
        dH = formation_energy(ph["formula"], ph["E"], elements) / tot
        pts.append(([c / tot for c in cnt], dH))
    tcnt = [float(target["formula"].get(e, 0)) for e in els]
    ttot = sum(tcnt)
    if ttot <= 0 or len(pts) < N:
        return None
    xt = [c / ttot for c in tcnt]
    dHt = formation_energy(target["formula"], target["E"], elements) / ttot
    import itertools
    best = None
    for combo in itertools.combinations(range(len(pts)), N):
        # 解 Σ_p λ_p x_p = x_t（N 个分量、N 个未知量；Σλ=1 自动成立）
        A = [[pts[idx][0][i] for idx in combo] for i in range(N)]
        lam = _solve_linear(A, xt)
        if lam is None:
            continue
        if any(l < -TOL for l in lam):
            continue
        val = sum(lam[k] * pts[combo[k]][1] for k in range(N))
        if best is None or val < best:
            best = val
    if best is None:
        return None
    # 标准约定：目标在其它相凸包之上时为正，落在/低于凸包时为 0
    return max(0.0, dHt - best)


def convex_hull_window(target, elements, phases):
    """通用 N 元（N>=1）化学势稳定性窗口。返回 (顶点列表, dH_target)。

    约束：dmu_i <= 0；sum_i n_i(target)*dmu_i = dH_target（主元元素消去）；
    每个竞争相 sum_i n_i(phase)*dmu_i <= dH(phase)。在 N-1 维自由变量空间枚举顶点。
    二元（Mg-C）、三元（A2B2Te5）、四元都走同一套。
    """
    els = list(elements.keys())
    N = len(els)
    if N == 0:
        raise SystemExit("[错误] elements 为空")
    n_t = [target["formula"].get(e, 0) for e in els]
    dH_target = formation_energy(target["formula"], target["E"], elements)
    piv = None
    for i in range(N - 1, -1, -1):
        if n_t[i] > 0:
            piv = i
            break
    if piv is None:
        raise SystemExit("[错误] 目标化合物成分全为 0")
    free = [i for i in range(N) if i != piv]
    npv = n_t[piv]
    raw = []
    for i in range(N):
        a = [0.0] * N
        a[i] = 1.0
        raw.append((els[i] + "<=0", a, 0.0))
    for ph in phases:
        a = [float(ph["formula"].get(e, 0)) for e in els]
        raw.append((ph.get("name", "phase"), a,
                    formation_energy(ph["formula"], ph["E"], elements)))
    red = []
    for nm, a, g in raw:
        A = [a[i] - a[piv] * n_t[i] / npv for i in free]
        G = g - a[piv] * dH_target / npv
        red.append((nm, A, G))
    d = len(free)
    verts, seen = [], set()
    if d == 0:
        # 无自由变量仍必须满足元素上界和所有竞争相约束。
        verts = [[]] if all(0.0 <= G + TOL for _, _, G in red) else []
    else:
        from itertools import combinations
        for combo in combinations(range(len(red)), d):
            sol = _solve_linear([red[c][1] for c in combo], [red[c][2] for c in combo])
            if sol is None:
                continue
            if all(sum(A[i] * sol[i] for i in range(d)) <= G + TOL for _, A, G in red):
                key = tuple(round(x, 6) for x in sol)
                if key not in seen:
                    seen.add(key)
                    verts.append(sol)
    if not verts:
        eah = energy_above_hull(target, elements, phases)
        if eah is not None and eah > 1e-3:
            raise SystemExit(
                "[错误] 凸包窗口为空：目标相比元素相/竞争相的凸包高 %.3f eV/atom（亚稳相），"
                "不存在同时满足平衡与竞争相稳定的化学势窗口。\n"
                "       处理：① 先核对相总能与成分；② 若确为亚稳相，把 energies.json 的 "
                "\"mu\" 手动设为一组物理化学势（例如 C 取石墨、Mg 取石墨+Mg2C3 共存点），"
                "不要写 mu_vertices；形成能脚本会直接采用该 mu。" % eah)
        raise SystemExit("[错误] 凸包窗口为空（%d 个顶点）—— 检查相总能/成分是否有误" % len(verts))
    if d == 2:
        import math
        cx = sum(v[0] for v in verts) / len(verts)
        cy = sum(v[1] for v in verts) / len(verts)
        verts = sorted(verts, key=lambda v: math.atan2(v[1] - cy, v[0] - cx))
    elif d == 1:
        verts = sorted(verts, key=lambda v: v[0])
    out = []
    for sol in verts:
        dmu = {}
        for k, i in enumerate(free):
            dmu[els[i]] = sol[k]
        dmu[els[piv]] = (dH_target - sum(n_t[i] * dmu[els[i]] for i in free)) / npv
        mu = {e: elements[e] + dmu[e] for e in els}
        out.append({"dmu": {e: round(v, 6) for e, v in dmu.items()},
                    "mu": {e: round(v, 6) for e, v in mu.items()}})
    return out, dH_target

def main():
    ph = load_phases()
    verts, dH_target = convex_hull_window(ph["target"], ph["elements"], ph["phases"])
    els = list(ph["elements"].keys())
    eah = energy_above_hull(ph["target"], ph["elements"], ph["phases"])
    print("[OK] 目标化合物 ΔH_f = %.4f eV/式量" % dH_target)
    if eah is not None:
        print("     目标相 energy above hull = %+.4f eV/atom%s"
              % (eah, "（亚稳）" if eah > 0.001 else "（稳定）"))
    print("     化学势窗口（%d 个顶点，μ 为绝对值 eV/原子）:" % len(verts))
    for k, v in enumerate(verts):
        mu = "  ".join("%s=%.4f" % (el, v["mu"][el]) for el in els)
        print("       顶点%d: %s" % (k, mu))
    json.dump({"target": ph["target"]["formula"], "dH_f": round(dH_target, 6),
               "energy_above_hull": (round(eah, 6) if eah is not None else None),
               "window_vertices": verts},
              open("chemical_potential_window.json", "w"), indent=2, ensure_ascii=False)
    print("     窗口已写 chemical_potential_window.json")
    # 默认取第一个顶点（可用 --vertex N 覆盖）写入 energies.json 供形成能脚本用
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--vertex", type=int, default=0, help="用哪个顶点作为默认化学势")
    args, _ = ap.parse_known_args()
    if 0 <= args.vertex < len(verts):
        ref = {"mu": verts[args.vertex]["mu"]}
        if os.path.exists("energies.json"):
            ref.update(json.load(open("energies.json", encoding="utf-8")))
        json.dump(ref, open("energies.json", "w"), indent=2, ensure_ascii=False)
        print("     已写 energies.json（mu 取顶点%d；E_gap/epsilon 请手动补）" % args.vertex)


def read_energy(outcar):
    """取 OUTCAR 最后 without entropy（E0）。"""
    if not os.path.exists(outcar):
        return None
    e = None
    for line in open(outcar, errors="ignore"):
        if "without entropy" in line:
            try:
                e = float(line.split("=")[1].split()[0])
            except (ValueError, IndexError):
                pass
    return e


def build_phases_from_references(ref_json="references_energy.json",
                                 bulk_poscar="step1_bulk/POSCAR",
                                 bulk_outcar="step1_bulk/OUTCAR"):
    """由 step0 的 references_energy.json + step1_bulk 组装凸包输入。

    目标化合物式量 = 超胞成分 / gcd（如 Sn18 Sb18 Te45 -> Sn2 Sb2 Te5）。
    目标相每式量总能优先取 references_energy.json 的 target_prim（独立 ISIF=3 密网格
    原胞，与参考相同口径）；没有时退回 step1_bulk 超胞/9 —— 那是 ISIF=2 + 粗 k 网格
    的口径，与参考相 ISIF=3 不一致，ΔH_f 不可信，只能作占位。"""
    ref = json.load(open(ref_json, encoding="utf-8"))
    import defects_common as D2
    st = D2.parse_poscar(bulk_poscar)
    counts = {}
    for a in st["atoms"]:
        counts[a] = counts.get(a, 0) + 1
    from math import gcd
    from functools import reduce
    g = reduce(gcd, counts.values())
    formula = {el: n // g for el, n in counts.items()}
    n_super = sum(counts.values())
    n_prim = sum(formula.values())
    E_super = read_energy(bulk_outcar)
    tp = ref.get("target_prim")
    if tp and tp.get("E_per_fu") is not None:
        E_per_fu = float(tp["E_per_fu"])
        src = "target_prim（独立 ISIF=3 密网格原胞，与参考相同口径）"
    else:
        if E_super is None:
            raise SystemExit("[错误] step1_bulk/OUTCAR 无能量（先跑完 step1）")
        E_per_fu = E_super * (n_prim / n_super)
        src = "step1_bulk 超胞/9（⚠ ISIF=2+粗网格，与参考相 ISIF=3 不同口径，ΔH_f 不可信）"
    # δ = E_super/9 - E_target_prim：口径残差。修完 target_prim 后若 δ 不小，
    # 它会通过 μ 直接进每个 E_f（系数 Δn_i），必须量出来并设闸门（阈值 10 meV/fu）。
    delta = None
    if E_super is not None and tp and tp.get("E_per_fu") is not None:
        delta = E_super * (n_prim / n_super) - E_per_fu
    els = set(formula)
    elements = {el: v for el, v in ref["elements"].items() if el in els}
    src_phases = dict(ref.get("binaries") or {})
    src_phases.update(ref.get("phases") or {})
    phases = [{"name": n, "formula": i["formula"], "E": i["E_per_fu"]}
              for n, i in src_phases.items() if set(i["formula"]).issubset(els)]
    return {"target": {"formula": formula, "E": E_per_fu, "energy_source": src,
                       "delta_eV_fu": delta},
            "elements": elements, "phases": phases}


def run_from_references():
    """S4 调用：references_energy.json + step1_bulk -> 化学势窗口 -> energies.json。"""
    ref_json = None
    for cand in ("step0_references/references_energy.json", "references_energy.json"):
        if os.path.exists(cand):
            ref_json = cand
            break
    if ref_json is None:
        raise SystemExit("[错误] 找不到 references_energy.json（step0 参考相未完成）")
    ph = build_phases_from_references(ref_json=ref_json)
    verts, dH_target = convex_hull_window(ph["target"], ph["elements"], ph["phases"])
    els = list(ph["elements"].keys())
    eah = energy_above_hull(ph["target"], ph["elements"], ph["phases"])
    degenerate_line = (len(verts) < len(els))
    print("[OK] %s ΔH_f = %.4f eV/式量，化学势窗口 %d 个顶点"
          % ("".join("%s%d" % (e, ph["target"]["formula"].get(e, 0)) for e in els),
             dH_target, len(verts)))
    print("     目标相能量来源: %s" % ph["target"].get("energy_source", "(未知)"))
    delta = ph["target"].get("delta_eV_fu")
    if delta is not None:
        flag = "⚠ 超 10 meV 阈值，E_f 会被 μ 端污染" if abs(delta) > 0.010 else "OK"
        print("     δ = E_super/9 - E_target_prim = %+.4f eV/fu   [%s]" % (delta, flag))
    for k, v in enumerate(verts):
        print("   顶点%d: %s" % (k, "  ".join("%s=%.4f" % (el, v["mu"][el]) for el in els)))
    if eah is not None:
        print("   目标相 energy above hull = %+.4f eV/atom%s"
              % (eah, "（亚稳）" if eah > 0.001 else "（稳定）"))
    json.dump({"target": ph["target"]["formula"], "dH_f": round(dH_target, 6),
               "energy_above_hull": (round(eah, 6) if eah is not None else None),
               "window_vertices": verts},
              open("chemical_potential_window.json", "w"), indent=2, ensure_ascii=False)
    # 默认化学势取顶点 0（某元素最富的一端）；形成能脚本会遍历全部顶点，
    # 因此默认取哪个不影响最终报告，只是 results['mu'] 的入口值。
    ref = json.load(open("energies.json", encoding="utf-8")) if os.path.exists("energies.json") else {}
    ref["mu"] = verts[0]["mu"]
    ref["mu_vertices"] = verts   # 全部顶点，供形成能脚本遍历各极端化学势
    ref["degenerate_line"] = degenerate_line
    json.dump(ref, open("energies.json", "w"), indent=2, ensure_ascii=False)
    print("     已写 energies.json（mu 取顶点0 + 全部 %d 个顶点 mu_vertices%s）"
          % (len(verts), "，窗口退化为线段" if degenerate_line else ""))


if __name__ == "__main__":
    main()
