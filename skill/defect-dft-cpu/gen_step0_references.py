#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""defect-dft-cpu step0（run:gen）：凸包参考相生成/提交/收集能量（通用版）。

幂等 + 多材料共享：参考相只在 REFERENCES_DIR（step.conf）算一份。扫描
<REFERENCES_DIR>/convex_hull_references/*/POSCAR 自动发现相，不再写死清单。

输出 references_energy.json：
  { "elements": {"Mg": 每原子 eV, "C": 每原子 eV},
    "phases":   {"Mg2C3": {"formula": {...}, "E_per_fu": eV, "counts": {...}}},
    "binaries": {<同 phases，向后兼容>} }
"""
import sys, os, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import defects_common as D
from gen_references import (discover_phases, resolve_spec, load_conf_for,
                            formula_name)

STEP = "step0_references"


def ref_energy(outcar):
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
    conf = D.load_stepconf()
    refdir = Path(conf.get("REFERENCES_DIR", "/public/home/.../convex_hull_refs"))
    potdir = conf.get("POTCAR_DIR", "/public/home/.../software/vasp_pseudopotentials")
    user = conf.get("USER", os.environ.get("USER", "user"))
    posdir = refdir / "convex_hull_references"
    os.makedirs(STEP, exist_ok=True)

    refconf = load_conf_for(posdir, conf)
    phases = discover_phases(posdir)
    if not phases:
        raise SystemExit("[错误] %s 下没有 <相>/POSCAR" % posdir)

    # 1) 生成 + 提交（幂等：已有 INCAR 不重生成；已在队列/跑过不重提交）
    for name, info in phases.items():
        d = posdir / name
        d.mkdir(parents=True, exist_ok=True)
        spec = resolve_spec(name, info, refconf)
        if not (d / "INCAR").exists():
            soc = str(conf.get("SOC", "1")).lower() not in ("0", "false", ".false.", "no", "off")
            ncomp = 3 if soc else 1
            from gen_references import INCAR
            incar = (INCAR.replace("{{name}}", name)
                          .replace("{{ENCUT}}", str(conf.get("ENCUT", 370)))
                          .replace("{{ismear}}", str(spec["ismear"]))
                          .replace("{{sigma}}", str(spec["sigma"]))
                          .replace("{{ISPIN}}", str(conf.get("ISPIN", "1")))
                          .replace("{{NCORE}}", str(conf.get("NCORE", "8")))
                          .replace("{{KPAR}}", str(conf.get("KPAR", "2")))
                          .replace("{{LSORBIT_LINE}}", "LSORBIT = .TRUE." if soc else "# LSORBIT off")
                          .replace("{{GGA_COMPAT_LINE}}", "GGA_COMPAT = .FALSE." if soc else "# GGA_COMPAT off")
                          .replace("{{magmom}}", "%d*0" % (ncomp * info["natoms"])))
            (d / "INCAR").write_text(incar, encoding="utf-8")
            D.write_kpoints(str(d / "KPOINTS"), spec["kmesh"])
            D.assemble_potcar(info["order"], potdir, out_path=str(d / "POTCAR"))
            D.render_submit(D.find_submit_tpl(soc), str(d / "submit.sh"), "ref_" + name)
        if not (d / "OUTCAR").exists():
            ok, msg = D.guarded_sbatch(str(d), "ref_" + name, user=user)
            print("  [提交] %s %s" % (name, msg))
            if not ok:
                print("  [警告] %s 提交未完成：%s" % (name, msg))

    # 2) 检查收敛 + 收集能量（发现的每个相都进输出，全部收敛才算完成）
    energies = {}
    all_done = True
    for name in phases:
        d = posdir / name
        e = ref_energy(str(d / "OUTCAR"))
        if e is None:
            print("  %-12s 尚未出结果" % name)
            all_done = False
            continue
        conv = "reached required accuracy" in open(str(d / "OUTCAR"), errors="ignore").read()
        if not conv:
            all_done = False
        energies[name] = e
        print("  %-12s E0=%.6f eV  收敛=%s" % (name, e, conv))

    if all_done:
        out = {"elements": {}, "phases": {}, "binaries": {}}
        for name, info in phases.items():
            e = energies.get(name)
            if e is None:
                continue
            els = list(info["counts"])
            if len(els) == 1 and info["natoms"] == info["nfu"]:
                out["elements"][els[0]] = e / info["natoms"]
            else:
                rec = {"formula": info["formula"], "E_per_fu": e / info["nfu"],
                       "counts": info["counts"], "E_total": e}
                key = formula_name(info["formula"])
                if key in out["phases"]:
                    key = name
                out["phases"][key] = rec
                out["binaries"][key] = {"formula": info["formula"], "E_per_fu": e / info["nfu"]}
        json.dump(out, open(refdir / "references_energy.json", "w"), indent=2, ensure_ascii=False)
        json.dump(out, open(os.path.join(STEP, "references_energy.json"), "w"),
                  indent=2, ensure_ascii=False)
        print("[OK] 参考相全部收敛，能量已写 %s" % (refdir / "references_energy.json"))
    else:
        print("[进行中] 参考相未全收敛，下次 watch 自动再查")


if __name__ == "__main__":
    main()
