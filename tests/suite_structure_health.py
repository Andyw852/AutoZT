# -*- coding: utf-8 -*-
"""自测：结构体检四判据（structure_health.py）。

判据（user 2026-09-21 列）：
  ① 最近邻 < 0.7×共价半径和  ② CN≤1 悬挂  ③ 命名 vs 计量比  ④ z 跨度 vs 层数。
第五条（双容差空间群）在 suite_symmetry_audit.py。

本套件用**构造的已知坏结构** + 真实材料（存在才查）做断言；缺 numpy 自动 SKIP（exit 0）。

跑法：
    python3 tests/suite_structure_health.py
    ~/miniconda3/envs/atomate2_p_a/bin/python tests/suite_structure_health.py
"""
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "skill", "_common", "opt"))

import structure_health as SH  # noqa: E402

FAILS = []
P1 = "/mnt/d/tf_data/test_TE/P1"
MOS2 = "/mnt/d/tf_data/test_kl/MoS2_kltest"


def ok(cond, msg):
    print(("  . " if cond else "  x ") + msg)
    if not cond:
        FAILS.append(msg)


def hexlat(a, b, gamma_deg, c=25.0):
    g = math.radians(gamma_deg)
    return [[a, 0.0, 0.0],
            [b * math.cos(g), b * math.sin(g), 0.0],
            [0.0, 0.0, c]]


def mos2(healthy=True):
    lat = hexlat(3.16, 3.16, 120.0, 25.0)
    x = 1.0 / 3 if healthy else 0.0
    y = 2.0 / 3 if healthy else 0.0
    frac = [[0.0, 0.0, 0.5], [x, y, 0.4372], [x, y, 0.5628]]
    return lat, frac, ["Mo", "S", "S"]


def main():
    if SH._import_numpy() is None:
        print("[SKIP] numpy 不可用 —— 跳过 structure_health 自测（不算失败）")
        return 0

    print("[1] ① 重叠判据：MoS₂ 坏种子（Mo、S 同面内位 -> d=1.57 Å）必须 FAIL")
    lat, frac, syms = mos2(healthy=False)
    r = SH.check_structure(lat, frac, syms, material_name="P1_Mo-MoS2_x_MoS2")
    ok(r["ok"] is False, "总体 ok=False")
    ok(r["nearest"]["ok"] is False, "① nearest.ok=False")
    ok(abs(r["nearest"]["dist_A"] - 0.0628 * 25.0) < 0.02,
       "最近邻 d=%.4f Å（0.0628×25）" % r["nearest"]["dist_A"])
    ok(r["nearest"]["dist_over_sum"] < 0.7,
       "d/(rcov和)=%.3f < 0.7" % r["nearest"]["dist_over_sum"])
    ok(any(f.startswith("①") for f in r["flags"]), "flags 含 ①")

    print("[2] ① 健康 MoS₂ 通过；② 周期镜像 CN 必须给出 S=3 / Mo=6（原胞只有 3 原子）")
    lat, frac, syms = mos2(healthy=True)
    r = SH.check_structure(lat, frac, syms, material_name="P1_Mo-MoS2_x_MoS2")
    ok(r["nearest"]["ok"] is True, "① 通过（d=%.3f Å）" % r["nearest"]["dist_A"])
    cn = r["coordination"]["cn"]
    ok(cn == [6, 3, 3], "CN = %s（Mo=6, S=3/3）—— 若为 [1,2] 说明没算周期镜像" % cn)
    ok(r["coordination"]["n_dangling"] == 0, "无悬挂")
    ok(r["ok"] is True and not r["flags"], "四判据全通过")

    print("[3] ③ 命名 vs 计量比：名字元素与结构不符 -> FAIL；一致 -> OK")
    lat, frac, syms = mos2(healthy=True)
    r_bad = SH.check_structure(lat, frac, syms, material_name="P1_Mo-MoS2_x_MoSe2")
    ok(r_bad["stoichiometry"]["ok"] is False, "MoSe2 名字 vs MoS2 结构 -> ③ FAIL")
    r_ok = SH.check_structure(lat, frac, syms, material_name="P1_Mo-MoS2_x_MoS2")
    ok(r_ok["stoichiometry"]["ok"] is True, "MoS2 名字 -> ③ OK")
    r_skip = SH.check_structure(lat, frac, syms, material_name="MoS2_kltest")
    ok(r_skip["stoichiometry"]["ok"] is True and "不是化学式" in r_skip["stoichiometry"]["note"],
       "末段非化学式（kltest）-> 跳过，不误报")

    print("[4] ② 悬挂判据：孤立原子（无近邻）必须 FAIL")
    lat = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 20.0]]
    frac = [[0.0, 0.0, 0.5], [0.5, 0.5, 0.5]]
    r = SH.check_structure(lat, frac, ["He", "He"], material_name=None)
    ok(r["coordination"]["n_dangling"] == 2, "2 个 He 相距 3.54 Å > 2.8 -> 都悬挂")
    ok(any(f.startswith("②") for f in r["flags"]), "flags 含 ②")

    print("[5] ④ z 跨度：2D 层起伏 > 6 Å -> FAIL")
    lat = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 30.0]]
    frac = [[0.0, 0.0, 0.05], [0.0, 0.0, 0.35]]     # 跨度 9 Å，同一列
    r = SH.check_structure(lat, frac, ["C", "C"], material_name=None)
    ok(any(f.startswith("④") for f in r["flags"]), "flags 含 ④（span=%.2f Å）" % r["layers"]["span_A"])

    print("[6] 真实材料（存在才查）：AlN/Zn5O3/BeO 必须被拦，四个候选必须过")
    bad = {"P1_Al-AlN_A1Z4_A-1-4_Al5N": False,
           "P1_Zn-ZnO_Z3A2_Z-3-2_Zn5O3": False,
           "P1_Be-BeO_A2Z2_A-2-2_Be2O": False}
    good = ["P1_Mo-MoS2_A4-3-1_A4-3-1_Mo2S3", "P1_Mo-MoS2_Z3-2-1_Z3-2-1_Mo3S4",
            "P1_Mo-MoS2_Z4-3-1_Z4-3-1_Mo2S3", "P1_Ti-TiS2_Z4-3-1_Z4-3-1_Ti2S3"]
    real_bad = real_good = 0
    for name, _ in bad.items():
        p = os.path.join(P1, name, "POSCAR")
        if not os.path.isfile(p):
            continue
        rr = SH.check_structure_file(p, material_name=name)
        ok(rr["ok"] is False, "%s 被拦：%s" % (name.split("_", 1)[1], "；".join(rr["flags"])))
        real_bad += 1
    for name in good:
        p = os.path.join(P1, name, "POSCAR")
        if not os.path.isfile(p):
            continue
        rr = SH.check_structure_file(p, material_name=name)
        ok(rr["ok"] is True, "%s 通过" % name.split("_", 1)[1])
        real_good += 1
    mp = os.path.join(MOS2, "POSCAR")
    if os.path.isfile(mp):
        rr = SH.check_structure_file(mp, material_name="MoS2_kltest")
        ok(rr["ok"] is True, "MoS2_kltest 通过（CN=%s）" % rr["coordination"]["cn"])
    if real_bad == 0 and real_good == 0:
        print("  (真实材料不在本机，跳过)")

    print()
    if FAILS:
        print("FAIL %d 项：" % len(FAILS))
        for m in FAILS:
            print("  - %s" % m)
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
