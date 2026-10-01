#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_isif2_static.py —— 原胞弛豫晶格 -> 超胞 -> ISIF=2 静态单点。

量 delta = E_super/N - E_prim：隔离超胞 vs 原胞的 k 收敛/尺寸残差。
通用版：元素顺序/POTCAR/ENCUT/SOC/NCORE/KPAR 全部从原胞 POSCAR + step.conf 得。

用法：python3 run_isif2_static.py [原胞CONTCAR] [输出目录] [nx ny nz]
"""
import sys, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import defects_common as D


def main():
    prim = sys.argv[1] if len(sys.argv) > 1 else "convex_hull_references/target_p8/CONTCAR"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "convex_hull_references/target_super_static"
    dims = tuple(int(x) for x in sys.argv[3].split()) if len(sys.argv) > 3 else (3, 3, 1)
    if not os.path.exists(prim):
        raise SystemExit("[错误] 找不到 %s（先等原胞弛豫收敛）" % prim)
    conf = D.load_stepconf()
    st = D.parse_poscar(prim)
    order = list(dict.fromkeys(st["atoms"]))
    sc = D.supercell(st, dims)
    os.makedirs(outdir, exist_ok=True)
    D.write_poscar(os.path.join(outdir, "POSCAR"), sc,
                   comment="static supercell %s from relaxed primitive" % "x".join(map(str, dims)))
    soc = str(conf.get("SOC", "1")).lower() not in ("0", "false", ".false.", "no", "off")
    ncomp = 3 if soc else 1
    natoms = len(sc["atoms"])
    incar = "\n".join([
        "SYSTEM = super static (ISIF=2)",
        "ISTART = 0", "ICHARG = 2", "GGA = PE", "IVDW = 12", "PREC = Accurate",
        "ENCUT = %s" % conf.get("ENCUT", 370),
        "LREAL = .FALSE.", "LASPH = .TRUE.", "ALGO = All", "AMIX = 0.1", "BMIX = 0.0001",
        "EDIFF = 1E-6", "NELM = 200", "NELMIN = 6", "ISMEAR = 0", "SIGMA = 0.05",
        "ISPIN = %s" % conf.get("ISPIN", "1"), "ISYM = 0", "IBRION = -1", "ISIF = 2", "NSW = 0",
        "LWAVE = .FALSE.", "LCHARG = .FALSE.",
        "NCORE = %s" % conf.get("NCORE", "8"), "KPAR = %s" % conf.get("KPAR", "2"),
        "LSORBIT = .TRUE." if soc else "# LSORBIT off",
        "GGA_COMPAT = .FALSE." if soc else "# GGA_COMPAT off",
        "LMAXMIX = 4", "MAGMOM = %d*0" % (ncomp * natoms), ""])
    open(os.path.join(outdir, "INCAR"), "w").write(incar)
    open(os.path.join(outdir, "KPOINTS"), "w").write("Auto mesh (supercell)\n0\nGamma\n2 2 2\n0 0 0\n")
    D.assemble_potcar(order, conf.get("POTCAR_DIR", "/public/home/<user>/software/vasp_pseudopotentials"),
                      out_path=os.path.join(outdir, "POTCAR"))
    D.render_submit(D.find_submit_tpl(soc), os.path.join(outdir, "submit.sh"), "static_super")
    print("[OK] ISIF=2 静态超胞就绪: %s" % outdir)
    print("     %d 原子, KMESH 2 2 2, 元素=%s" % (natoms, ",".join(order)))
    print("     跑完用: E_super/%d - E_prim = delta" % (dims[0] * dims[1] * dims[2]))


if __name__ == "__main__":
    main()
