# -*- coding: utf-8 -*-
"""自测：结构体检第五条 —— 数值微畸变审计 + spglib 对称化（symmetry_audit.py）。

判据（wangchao 指定）：symprec=1e-5 与 1e-4 给出的空间群不一致 => 数值微畸变。
本套件用**构造的已知微畸变结构**与 MoS2 实测数字做断言；需要 spglib+numpy，
缺了自动跳过（exit 0，不误报失败）。

跑法（建议用装了 spglib 的 python；system python3 会 SKIP）：
    python3 tests/suite_symmetry_audit.py
    /home/wangchao/miniconda3/envs/atomate2_p_a/bin/python tests/suite_symmetry_audit.py
"""
import math
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "skill", "_common", "opt"))

import symmetry_audit as SA  # noqa: E402

FAILS = []


def ok(cond, msg):
    print(("  . " if cond else "  x ") + msg)
    if not cond:
        FAILS.append(msg)


def hexlat(a, b, gamma_deg, c=20.0):
    g = math.radians(gamma_deg)
    return [[a, 0.0, 0.0],
            [b * math.cos(g), b * math.sin(g), 0.0],
            [0.0, 0.0, c]]


def close(x, y, tol=1e-9):
    return abs(x - y) <= tol


def main():
    if not SA.spglib_available():
        print("[SKIP] spglib 不可用 —— 跳过 symmetry_audit 自测（不算失败）")
        return 0

    print("[1] 构造已知微畸变结构：1e-5 与 1e-4 空间群不一致 -> 触发")
    lat = hexlat(3.0, 3.0 + 1e-5, 120.0001)     # a != b、gamma 偏 120 度
    frac = [[0.0, 0.0, 0.5], [1.0 / 3, 2.0 / 3, 0.5]]
    syms = ["B", "N"]
    ent, micro = SA.spacegroup_audit(lat, frac, syms)
    ok(micro is True, "micro_distortion=True（%s）"
       % [(e["symprec"], e["international"], e["n_ops"]) for e in ent])
    ok(ent[0]["number"] == 38 and ent[1]["number"] == 187,
       "symprec 1e-5 -> Amm2(#38)，1e-4 -> P-6m2(#187)")

    print("[2] 对称化后两个容差一致")
    rl, rp, rs, used = SA.symmetrize(lat, frac, syms, symprec=1e-4)
    rl = rl.tolist() if hasattr(rl, "tolist") else rl
    rp = rp.tolist() if hasattr(rp, "tolist") else rp
    ent2, micro2 = SA.spacegroup_audit(rl, rp, rs)
    ok(micro2 is False and all(e["number"] == 187 for e in ent2),
       "对称化后 1e-5 与 1e-4 都给 P-6m2：%s"
       % [(e["symprec"], e["international"]) for e in ent2])
    m = SA.lattice_metrics(rl)
    ok(close(m["a"], m["b"], 1e-12), "a == b (%.12f)" % m["a"])
    ok(close(m["gamma"], 120.0, 1e-8), "gamma == 120 (%.12f)" % m["gamma"])

    print("[3] MoS2 真实 S1 数字：审计与修复")
    a = 3.1401970575143987
    b = 3.140189175214
    g = 120.000083034652
    latm = hexlat(a, b, g, 25.0)
    fracm = [[0.0, 0.0, 0.5],
             [0.3333333329999988, 0.6666666670000012, 0.5623543263992943],
             [0.3333333329999988, 0.6666666670000012, 0.4376456736007057]]
    symm = ["Mo", "S", "S"]
    em, mm = SA.spacegroup_audit(latm, fracm, symm)
    ok(mm is True, "检出微畸变（|a-b|=%.3e A, gamma-120=%.3e deg）"
       % (abs(a - b), g - 120.0))
    rlm, rpm, rsm, _ = SA.symmetrize(latm, fracm, symm, symprec=1e-4)
    rlm = rlm.tolist() if hasattr(rlm, "tolist") else rlm
    rpm = rpm.tolist() if hasattr(rpm, "tolist") else rpm
    conf = SA.confirmations({"lattice": latm, "frac": fracm, "symbols": symm},
                            {"lattice": rlm, "frac": rpm, "symbols": rsm})
    ok(conf["c1_stoichiometry"]["ok"] is True, "① 原子数/计量比不变 （3 原子 MoS2）")
    ok(conf["c3_vacuum_axis"]["ok"] is True, "③ 真空轴未移动（c 轴）")
    ok(conf["c3_vacuum_axis"]["before"]["axis"] == 2
       and conf["c3_vacuum_axis"]["after"]["axis"] == 2, "真空轴下标 2 -> 2")
    ok(conf["c2_energy"]["status"] == "pending", "② 能量确认 pending（待跑 2 个单点）")
    ok(conf["all_ok"] is True and conf["fully_confirmed"] is False,
       "几何两项通过；未跑能量前 fully_confirmed=False（不静默通过）")
    vi = SA.vacuum_info(rlm, rpm)
    ok(vi["is_2d"] is True and close(vi["vacuum_A"], 21.882283680035286, 1e-6),
       "真空厚度 %.4f A（MoS2 25 A 胞）" % vi["vacuum_A"])

    print("[4] 确认①护栏：refine 若改变原胞选取（原子数变）必须不过")
    rl5, rp5, rs5, _ = SA.symmetrize(latm, fracm, symm, symprec=1e-5)
    ok(len(rs5) != len(symm),
       "symprec=1e-5 的 refine_cell 把 3 原子变成 %d 原子（原胞选取变了）" % len(rs5))
    conf_bad = SA.confirmations({"lattice": latm, "frac": fracm, "symbols": symm},
                                {"lattice": rl5.tolist() if hasattr(rl5, "tolist") else rl5,
                                 "frac": rp5.tolist() if hasattr(rp5, "tolist") else rp5,
                                 "symbols": rs5})
    ok(conf_bad["c1_stoichiometry"]["ok"] is False
       and conf_bad["all_ok"] is False, "原子数变化时 ① 报 False、all_ok=False")

    print("[5] 能量确认：OUTCAR 解析 + 噪声内/超阈值")
    d = tempfile.mkdtemp(prefix="symtest_")
    try:
        f1 = os.path.join(d, "OUTCAR.a")
        f2 = os.path.join(d, "OUTCAR.b")
        open(f1, "w").write("  energy(sigma->0) =     -100.00000000\n")
        open(f2, "w").write("  energy(sigma->0) =     -100.00009000\n")
        ok(close(SA.parse_outcar_energy(f1), -100.0, 1e-12), "解析 energy(sigma->0)")
        r_pass = SA.energy_check(f1, f2, natoms=3, tol_mev_per_atom=1.0)
        ok(r_pass["status"] == "pass", "3 原子差 0.09 meV/atom < 1.0 -> pass")
        open(f2, "w").write("  energy(sigma->0) =      -99.0\n")
        r_fail = SA.energy_check(f1, f2, natoms=3, tol_mev_per_atom=1.0)
        ok(r_fail["status"] == "fail", "差 333 meV/atom > 1.0 -> fail")

        print("[6] 面内应力复检：解析 in kB 行")
        o1 = os.path.join(d, "OUTCAR.s1")
        o2 = os.path.join(d, "OUTCAR.s2")
        open(o1, "w").write("  in kB       -12.345   -10.000     1.000    -0.500     0.0     0.0\n"
                            "  external pressure =        7.50 kB\n")
        open(o2, "w").write("  in kB        -0.050    -0.040     0.100     0.010     0.0     0.0\n"
                            "  external pressure =        0.03 kB\n")
        ok(close(SA.parse_inplane_stress_kb(o1), 12.345, 1e-9),
           "面内 max|sigma| = 12.345 kB（取 xx/yy/xy）")
        st = SA.stress_recheck(o1, o2, tol_kb=0.2)
        ok(st["status"] == "pass" and close(st["inplane_after_kB"], 0.05, 1e-9),
           "对称化后面内 0.05 kB < 0.2 -> pass")

        print("[7] 对称化留档：before/after 结构 + JSON + log + inplace 备份")
        pos = os.path.join(d, "POSCAR")
        SA.write_poscar(pos, latm, fracm, symm, title="MoS2 test")
        arch = os.path.join(d, "symmetry_audit")
        res = SA.symmetrize_structure_file(
            pos, out=os.path.join(d, "POSCAR_sym"), archive_dir=arch,
            symprec=1e-4, allow_energy_pending=True, write_prov=False)
        ok(res["triggered"] is True and res["fit_for_use"] is False,
           "triggered=True，fit_for_use=False（能量 pending）")
        for name in ("POSCAR.before", "POSCAR.symmetrized", "symmetry_audit.json",
                     "symmetry_audit.log"):
            ok(os.path.isfile(os.path.join(arch, name)), "留档 %s" % name)
        ok(os.path.isfile(os.path.join(d, "POSCAR_sym")), "对称化结构另存 POSCAR_sym")
        # 往返读回：a==b、gamma==120
        back = SA.read_poscar(os.path.join(d, "POSCAR_sym"))
        mb = SA.lattice_metrics(back["lattice"])
        ok(close(mb["a"], mb["b"], 1e-12) and close(mb["gamma"], 120.0, 1e-8),
           "读回 POSCAR_sym：a==b、gamma==120")
        # inplace + 能量 pending 未放行 -> SystemExit
        raised = False
        try:
            SA.symmetrize_structure_file(pos, inplace=True, archive_dir=arch,
                                         symprec=1e-4, allow_energy_pending=False,
                                         write_prov=False)
        except SystemExit:
            raised = True
        ok(raised, "能量 pending 且未 --allow-energy-pending -> SystemExit（不静默通过）")
        ok(os.path.isfile(os.path.join(d, "POSCAR.pre_symmetry")),
           "inplace 前备份 POSCAR.pre_symmetry")
    finally:
        shutil.rmtree(d, ignore_errors=True)

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
