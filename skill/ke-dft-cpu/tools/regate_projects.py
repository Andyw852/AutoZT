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
  ★ WRONG   IBZ + 真实重叠 + 原公式有坏操作 + 没打补丁 —— **静默算错**，要重算；
            全网格 h5 + 真实重叠 + 有 k 点标签不在 (-0.5, 0.5]（V117：from_data 贴错系数标签）
  ⚠ STALE-DP 形变势 h5 没有 dp_symmetrized 标记 = S7.1 早于 V121 生成（V122：MoS2 同一批形变单点
            新版 S7.1 重读后电子 ADP 迁移率 198 -> 406）-> 重新 gen S7.1（秒级）再重跑本步
  ★ WRONG   （V134）settings.yaml 的弹性张量不正定（2D 看面内块、3D 看 6×6）且散射含 ADP —— ADP 无意义
  ⚠ ELASTIC （V134）2D 面内剪切 C66 不到 max(C11, C22) 的 2% —— 疑似 VASP 顺序（XX YY ZZ XY YZ ZX）没重排，
            AMSET 拿到的面内剪切近零；核对 step6_elastic/OUTCAR 的表头后重新 gen 本步
  ★ WRONG   （V135）settings.yaml 的弹性张量是 VASP 顺序（与 step6_elastic/OUTCAR 逐位对照：剪切对角三位依次是
            XY、YZ、ZX），没重排成标准 Voigt —— 3D 非立方是 C44/C66 互换，2D 是面内剪切近零；重新 gen 本步。
            立方等三个剪切相等的体系分不出（也不受影响），不判。
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
sys.path.insert(0, str(_HERE.parent / "step8.4_amset2d"))       # overlap_preflight（V117 标签检查）
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


def _outcar_elastic_std(mat):
    """用 S8 gen 的 read_elastic 从 step6_elastic/OUTCAR 读弹性张量（按表头重排成标准 Voigt，GPa）。读不到返回 None。"""
    import contextlib
    import io
    try:
        sys.path.insert(0, str(_HERE.parent / "step8_amset"))
        with contextlib.redirect_stdout(io.StringIO()):
            import gen_step10_amset as g10
            el = g10.read_elastic(Path(mat))
    except Exception:                                            # noqa: BLE001
        return None
    return el if isinstance(el, list) and len(el) == 6 else None


def voigt_order_state(settings_el, outcar_std, match=0.10, apart=0.25):
    """settings.yaml 的弹性张量（可能已按 c/t 重标度、对称化、2D 面外取绝对值）是标准 Voigt 还是 VASP 顺序。
    先用 C11/C22 定标度 k，再把剪切对角三位与 OUTCAR 的两种排法比：
        标准 Voigt：(YZ, XZ, XY)；VASP 顺序：(XY, YZ, ZX)。
    返回 "voigt" / "vasp" / None（判不出来：两种排法几乎一样（立方等），或都对不上）。"""
    import numpy as np
    try:
        S = np.abs(np.asarray(settings_el, float))
        R = np.abs(np.asarray(outcar_std, float))
        if S.shape != (6, 6) or R.shape != (6, 6) or min(R[0, 0], R[1, 1]) <= 0:
            return None
    except (TypeError, ValueError):
        return None
    k = 0.5 * (S[0, 0] / R[0, 0] + S[1, 1] / R[1, 1])
    s = np.array([S[3, 3], S[4, 4], S[5, 5]])
    std = k * np.array([R[3, 3], R[4, 4], R[5, 5]])
    vasp = k * np.array([R[5, 5], R[3, 3], R[4, 4]])
    scale = max(float(np.max(std)), 1e-9)
    if float(np.max(np.abs(std - vasp))) / scale < apart:
        return None                                    # 两种排法本来就差不多（立方/近立方）：分不出，也无所谓
    e_std = float(np.max(np.abs(s - std))) / scale
    e_vasp = float(np.max(np.abs(s - vasp))) / scale
    if e_std < match and e_vasp > apart:
        return "voigt"
    if e_vasp < match and e_std > apart:
        return "vasp"
    return None


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
    offconv = None
    if src and src.startswith("step4b") and h5.exists():
        try:                                                     # 只读 kpoints 数据集，很快
            import h5py
            import overlap_preflight as pf
            with h5py.File(str(h5), "r") as f:
                offconv = pf.n_offconvention_kpoints(f["kpoints"][()])
        except Exception:                                        # noqa: BLE001
            offconv = None
    dp_sym = None
    dh5 = run / "deformation.h5"
    if dh5.exists():
        try:
            import h5py
            with h5py.File(str(dh5), "r") as f:
                dp_sym = bool(int(f.attrs.get("dp_symmetrized", 0)))
        except Exception:                                        # noqa: BLE001
            dp_sym = None
    elastic, adp = None, True
    try:                                                         # V134：弹性张量（标准 Voigt，GPa）
        import yaml
        cfg = (yaml.safe_load(txt) or {}) if txt else {}
        elastic = cfg.get("elastic_constant")
        sc = cfg.get("scattering_type")
        adp = sc is None or sc == "auto" or "ADP" in [str(x).upper() for x in (sc or [])]
    except Exception:                                            # noqa: BLE001
        elastic = None
    two_d = (run / "2d_correction.json").is_file()
    return {"settings": bool(txt), "unity": unity, "h5_src": src, "offconv": offconv,
            "dp_sym": dp_sym, "elastic": elastic, "adp": adp, "two_d": two_d,
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
            elif full and r["offconv"]:
                verdict, why = "★ WRONG", ("全网格 h5 有 %d 个 k 点标签不在 (-0.5, 0.5] -> AMSET from_data "
                                          "贴错这些点的系数标签（V117）" % r["offconv"])
            elif full:
                verdict, why = "OK", "全网格 h5（from_data，不去对称化%s）" % (
                    "；标签未能读取" if r["offconv"] is None else "；标签符合 (-0.5, 0.5]")
            elif bad is None:
                verdict, why = "?", "取不到结构/对称操作"
            elif bad == 0:
                verdict, why = "OK", "IBZ，原公式 0/%d 坏操作" % tot
            elif r["fix"]:
                verdict, why = "OK", "IBZ，%d/%d 坏操作，但当次打了相位补丁" % (bad, tot)
            else:
                verdict, why = "★ WRONG", ("IBZ + 真实重叠，原公式 %d/%d 个操作算错、没打补丁 -> 静默算错"
                                          % (bad, tot))
            # V122：判据 OK 但形变势是旧版 S7.1 生成的 -> 结果要重跑（ADP 可差一倍）
            if verdict == "OK" and r["dp_sym"] is False:
                verdict, why = "⚠ STALE-DP", why + "；但形变势 h5 未对称化（S7.1 早于 V121）-> 重新 gen S7.1 再重跑"
            # V134：弹性张量不正定 -> ADP 无意义；2D 面内剪切近零 -> 疑似 Voigt 顺序没重排
            el = None
            if r["elastic"] is not None and r["adp"]:
                try:
                    el = kc.elastic_stability(r["elastic"], is_2d=r["two_d"])
                except Exception:                                # noqa: BLE001
                    el = None
            order = (voigt_order_state(r["elastic"], _outcar_elastic_std(mat))
                     if r["elastic"] is not None and r["adp"] else None)
            if order == "vasp":                          # V135：与 OUTCAR 逐位对照，确定没重排
                verdict = "★ WRONG"
                why += ("；弹性张量是 VASP 顺序（剪切对角三位依次是 XY、YZ、ZX，与 step6_elastic/OUTCAR 对照），"
                        "没重排成标准 Voigt -> 重新 gen 本步（V135）")
                _P = (0, 1, 2, 4, 5, 3)                  # VASP -> 标准 Voigt（与 gen 的 _VASP2VOIGT 相同）
                el2 = kc.elastic_stability([[r["elastic"][i][j] for j in _P] for i in _P], is_2d=r["two_d"])
                if el2 is not None and not el2["ok"]:
                    why += ("；而且重排后%s也不正定（本征值 %s GPa）-> 重新 gen 会被 V134 的闸门拦下，先重做 S6"
                            % ("面内块" if r["two_d"] else "", ", ".join("%.2f" % x for x in el2["eigs"])))
            if el is not None and not el["ok"]:
                verdict = "★ WRONG"
                why += "；弹性张量%s不正定（本征值 %s GPa）-> ADP 无意义（V134）" % (
                    "面内块" if r["two_d"] else "", ", ".join("%.2f" % x for x in el["eigs"]))
            elif el is not None and el["soft_c66"]:
                if verdict == "OK":
                    verdict = "⚠ ELASTIC"
                why += "；面内剪切 C66 只有 max(C11, C22) 的 %.1f%%（疑似 VASP 顺序没重排，V134）" % (
                    100 * el["block"][2][2] / max(el["block"][0][0], el["block"][1][1]))
            rows.append({"material": str(mat.relative_to(root)), "run": run_name,
                         "verdict": verdict, "why": why, "h5_src": r["h5_src"],
                         "unity": r["unity"], "desym_fix": r["fix"],
                         "transport": r["transport"], "bad_ops": bad, "total_ops": tot,
                         "structure": st_src, "elastic_order": order})
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
    n_stale = sum(1 for r in rows if r["verdict"].startswith("⚠"))
    print("\n共 %d 个运行目录，★ 静默算错 %d 个，⚠ 待处理（形变势过期 / 弹性待核对）%d 个"
          % (len(rows), n_bad, n_stale))
    if a.csv:
        with open(a.csv, "w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0]))
            wr.writeheader()
            wr.writerows(rows)
        print("写出 %s" % a.csv)
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
