#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""update_target_prim.py —— 把独立 ISIF=3 密网格原胞总能写入 references_energy.json 的 target_prim。

凸包的目标相 dH_f 必须用与参考相同口径（ISIF=3 + 密 k 网格）的原胞能量，
不能用 step1_bulk 超胞/9（ISIF=2 + 粗网格）。本脚本负责这个口径修复。

成分与式量数自动从原胞 POSCAR/CONTCAR 求（gcd 约化），不再写死 Sn2Sb2Te5。

用法：python3 update_target_prim.py [原胞目录]
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathlib import Path


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


def main():
    prim_dir = sys.argv[1] if len(sys.argv) > 1 else "convex_hull_references/target_p8"
    outcar = os.path.join(prim_dir, "OUTCAR")
    if not os.path.exists(outcar):
        raise SystemExit("[错误] 找不到 %s（先跑完原胞 ISIF=3 计算）" % outcar)
    if "reached required accuracy" not in open(outcar, errors="ignore").read():
        print("[警告] %s 尚未显示收敛标记，仍取最后能量（请确认收敛）" % outcar)
    E0 = read_energy(outcar)
    if E0 is None:
        raise SystemExit("[错误] %s 无能量" % outcar)
    poscar = os.path.join(prim_dir, "POSCAR")
    if not os.path.exists(poscar):
        poscar = os.path.join(prim_dir, "CONTCAR")
    if not os.path.exists(poscar):
        raise SystemExit("[错误] %s 缺 POSCAR/CONTCAR（无法确定成分）" % prim_dir)
    import defects_common as D
    st = D.parse_poscar(poscar)
    counts = {}
    for a in st["atoms"]:
        counts[a] = counts.get(a, 0) + 1
    g = 0
    for v in counts.values():
        g = math.gcd(g, int(v))
    g = g or 1
    formula = {el: n // g for el, n in counts.items()}
    ref_path = "references_energy.json"
    for cand in ("step0_references/references_energy.json", "references_energy.json"):
        if os.path.exists(cand):
            ref_path = cand
            break
    ref = json.load(open(ref_path, encoding="utf-8"))
    ref["target_prim"] = {"formula": formula, "E_per_fu": round(E0 / g, 8),
                          "dir": prim_dir, "counts": counts, "E_total": round(E0, 8)}
    json.dump(ref, open(ref_path, "w"), indent=2, ensure_ascii=False)
    print("[OK] target_prim: E_per_fu = %.8f eV/式量  formula=%s (%s) -> %s"
          % (E0 / g, formula, prim_dir, ref_path))


if __name__ == "__main__":
    main()
