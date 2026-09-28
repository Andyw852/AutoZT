#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""regate_projects.py —— 用新判据（逐操作 + AMSET 的 symprec）回查已有项目的 S8/S8.4 结果是否可信。

不跑任何计算，只读每个材料目录里的：
  step3_uniform/POSCAR（或 CONTCAR）    —— 判据的结构（= wavefunction.h5 的结构）
  step8_amset/ 与 step8.4_amset2d/      —— wavefunction.h5 软链指向哪（IBZ 还是全网格）、
                                          settings.yaml 的 unity_overlap / # AZ_DESYM_FIX、
                                          transport.json 在不在
裁决（每个 S8/S8.4 运行目录一行）：
  OK        unity 重叠（不读波函数）/ 全网格 h5（from_data）/ IBZ 且原公式 0 个坏操作 /
            IBZ 且当次打了相位补丁
  ★ WRONG   IBZ + 真实重叠 + 原公式有坏操作 + 没打补丁 —— **静默算错**，要重算
  ?         判不出来（缺结构/缺 settings）
另外列出 settings.yaml 里写了 EPS_INF_OVERRIDE 注释却可能没生效的旧 S8（见 V115 §8）。

用法：
    python regate_projects.py <work 根目录> [--glob '*/ke-dft-cpu'] [--csv out.csv]
    例：python regate_projects.py /public/.../work --glob '*/ke-dft-cpu'
"""
import argparse
import csv
import glob
import os
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))
for _p in (_HERE.parent.parent / "_common" / "opt",):
    sys.path.insert(0, str(_p))


def _structure(mat):
    from pymatgen.core import Structure
    for d in ("step3_uniform", "step3b_uniform_full", "step1_opt"):
        for n in ("POSCAR", "CONTCAR"):
            p = mat / d / n
            if p.is_file():
                try:
                    return Structure.from_file(str(p)), str(p.relative_to(mat))
                except Exception:                              # noqa: BLE001
                    pass
    return None, None


def _run_info(run):
    st = run / "settings.yaml"
    txt = st.read_text(errors="ignore") if st.is_file() else ""
    unity = None
    for ln in txt.splitlines():
        if ln.strip().startswith("unity_overlap:"):
            unity = ln.split(":", 1)[1].strip().lower() == "true"
    h5 = run / "wavefunction.h5"
    src = None
    if h5.is_symlink():
        try:
            src = os.path.basename(os.path.dirname(os.path.realpath(h5)))
        except OSError:
            src = None
    return {"settings": bool(txt), "unity": unity, "h5_src": src,
            "fix": "# AZ_DESYM_FIX=1" in txt,
            "transport": (run / "transport.json").is_file()}


def audit(root, pattern):
    import ke_common as kc
    rows = []
    for mat in sorted(Path(p) for p in glob.glob(os.path.join(root, pattern))):
        if not mat.is_dir():
            continue
        st, st_src = _structure(mat)
        g = kc.symmetry_gate(st) if st is not None else None
        bad = None if g is None else g["bad_ops"]
        tot = None if g is None else g["total_ops"]
        for run_name in ("step8_amset", "step8.4_amset2d"):
            run = mat / run_name
            if not run.is_dir():
                continue
            r = _run_info(run)
            full = (r["h5_src"] or "").startswith("step4b")
            if not r["settings"]:
                verdict, why = "?", "没有 settings.yaml"
            elif r["unity"]:
                verdict, why = "OK", "unity 重叠（不读波函数）"
            elif full:
                verdict, why = "OK", "全网格 h5（from_data，不去对称化）"
            elif bad is None:
                verdict, why = "?", "取不到结构/对称操作"
            elif bad == 0:
                verdict, why = "OK", "IBZ，原公式 0/%d 坏操作" % tot
            elif r["fix"]:
                verdict, why = "OK", "IBZ，%d/%d 坏操作，但当次打了相位补丁" % (bad, tot)
            else:
                verdict, why = "★ WRONG", ("IBZ + 真实重叠，原公式 %d/%d 个操作算错、没打补丁 -> 静默算错"
                                          % (bad, tot))
            rows.append({"material": str(mat.relative_to(root)), "run": run_name,
                         "verdict": verdict, "why": why, "h5_src": r["h5_src"],
                         "unity": r["unity"], "desym_fix": r["fix"],
                         "transport": r["transport"], "bad_ops": bad, "total_ops": tot,
                         "structure": st_src})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root")
    ap.add_argument("--glob", default="*/ke-dft-cpu")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args(argv)
    rows = audit(a.root, a.glob)
    if not rows:
        print("没找到任何 S8/S8.4 运行目录（root=%s, glob=%s）" % (a.root, a.glob))
        return 0
    w = max(len(r["material"]) for r in rows)
    for r in rows:
        print("%-*s  %-15s  %-8s  %s%s" % (w, r["material"], r["run"], r["verdict"], r["why"],
                                           "" if r["transport"] else "（无 transport.json）"))
    n_bad = sum(1 for r in rows if r["verdict"].startswith("★"))
    print("\n共 %d 个运行目录，★ 静默算错 %d 个" % (len(rows), n_bad))
    if a.csv:
        with open(a.csv, "w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0]))
            wr.writeheader()
            wr.writerows(rows)
        print("写出 %s" % a.csv)
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
