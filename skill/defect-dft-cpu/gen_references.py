#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_references.py —— 第0步：给凸包参考相生成 VASP 输入（通用版）。

不再写死相清单：扫描 <refdir>/*/POSCAR（或 CONTCAR）自动发现所有相，按成分判类型：
  单元素相 -> 元素相（能量归一化到每原子）
  多元素相 -> 竞争相（能量归一化到每式量 = 原胞原子数 / gcd）

每个相的 ISMEAR/SIGMA/K 网格默认按内置规则给，可用 references.conf（INI）逐相覆盖：
  [default]        # 全局默认
  ismear = 0
  sigma  = 0.05
  [Mg_hcp]         # 按 <相目录名> 覆盖
  ismear = 1
  sigma  = 0.10
  kmesh  = 16 16 12        # 或 kmesh_length = 60

查找顺序：step.conf 的 REFERENCES_CONF -> <refdir>/references.conf -> 本技能
templates/references.conf（内置了 A2B2Te5 历史参数，保证向后兼容）。

INCAR 的 ENCUT/SOC/ISPIN/NCORE/KPAR 与 step.conf 一致（化学势必须同口径）。

用法：
    python3 gen_references.py [--potcar-dir ...] [--refdir convex_hull_references] [--only a,b]
    （--submit 已禁用：请走 autozt 公共提交入口）
"""
import sys, os, argparse, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import defects_common as D

# 仅用于「自动默认 ISMEAR/SIGMA」的金属元素表（可用 references.conf 覆盖）
_METALS = set(("Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo "
               "Tc Ru Rh Pd Ag Cd In Sn Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
               "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Th Pa U Np Pu").split())

INCAR = """SYSTEM = ref {{name}}
ISTART = 0
ICHARG = 2
GGA    = PE
IVDW   = 12
PREC   = Accurate
ENCUT  = {{ENCUT}}
LREAL  = .FALSE.
LASPH  = .TRUE.
ALGO   = All
AMIX   = 0.1
BMIX   = 0.0001
EDIFF  = 1E-6
NELM   = 200
NELMIN = 6
ISMEAR = {{ismear}}
SIGMA  = {{sigma}}
ISPIN  = {{ISPIN}}
ISYM   = 0
IBRION = 2
ISIF   = 3
NSW    = 60
EDIFFG = -0.02
LWAVE  = .FALSE.
LCHARG = .FALSE.
NCORE  = {{NCORE}}
KPAR   = {{KPAR}}
{{LSORBIT_LINE}}
{{GGA_COMPAT_LINE}}
LMAXMIX = 4
MAGMOM = {{magmom}}
"""

# submit 模板不再写死：按 SOC 从 setting/<hpc>/templates（或技能 templates/）取
#   SOC=1 -> submit_ncl_3d.tpl (vasp_ncl)；SOC=0 -> submit_std_3d.tpl (vasp_std)


def lattice_lengths(lat):
    return [math.sqrt(sum(x * x for x in v)) for v in lat]


def auto_kmesh(lat, target=60.0):
    """粗默认：每个方向 n = max(1, round(target / |a_i|))。请用 references.conf 收敛。"""
    return tuple(max(1, int(round(float(target) / L))) for L in lattice_lengths(lat))


def reduced_formula(counts):
    """counts {el:n} -> (约化式 {el:n/g}, g)。"""
    g = 0
    for v in counts.values():
        g = math.gcd(g, int(v))
    g = g or 1
    return {el: int(n) // g for el, n in counts.items()}, g


def formula_name(formula):
    return "".join(("%s%d" % (e, n)) if n > 1 else e for e, n in formula.items())


def parse_ref_conf(*paths):
    import configparser
    cp = configparser.ConfigParser()
    cp.read([str(p) for p in paths if p and os.path.exists(str(p))])
    return {sec: dict(cp.items(sec)) for sec in cp.sections()}


def load_conf_for(refdir, conf=None):
    """按优先级收集 references.conf（后者覆盖前者）：
    技能内置 templates/references.conf < cwd（tf 推送副本）< <refdir>/references.conf
    < step.conf 的 REFERENCES_CONF。"""
    paths = [os.path.join(str(Path(__file__).resolve().parent), "templates", "references.conf"),
             os.path.join(os.getcwd(), "references.conf"),
             os.path.join(str(refdir), "references.conf")]
    if conf and conf.get("REFERENCES_CONF"):
        paths.append(conf["REFERENCES_CONF"])
    return parse_ref_conf(*paths)


def discover_phases(refdir):
    """扫描 <refdir>/*/POSCAR|CONTCAR，返回有序 {name: info}。"""
    import collections
    out = collections.OrderedDict()
    root = Path(refdir)
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        p = d / "POSCAR"
        if not p.exists():
            p = d / "CONTCAR"
        if not p.exists():
            continue
        st = D.parse_poscar(str(p))
        counts = {}
        for a in st["atoms"]:
            counts[a] = counts.get(a, 0) + 1
        formula, nfu = reduced_formula(counts)
        out[d.name] = {"counts": counts, "formula": formula, "nfu": nfu,
                       "natoms": len(st["atoms"]),
                       "order": list(dict.fromkeys(st["atoms"])), "lat": st["lat"]}
    return out


def resolve_spec(name, info, conf):
    """INI 优先级：相 section > [default] > 自动。"""
    els = list(info["counts"])
    is_metal = all(e in _METALS for e in els)
    spec = {"ismear": 1 if is_metal else 0,
            "sigma": 0.10 if is_metal else 0.05,
            "kmesh": auto_kmesh(info["lat"]),
            "kind": "金属" if is_metal else "半导体/绝缘体"}
    for sec in ("default", name):
        s = conf.get(sec)
        if not s:
            continue
        if "ismear" in s:
            spec["ismear"] = int(float(s["ismear"]))
        if "sigma" in s:
            spec["sigma"] = float(s["sigma"])
        if "kmesh_length" in s:
            spec["kmesh"] = auto_kmesh(info["lat"], float(s["kmesh_length"]))
        if "kmesh" in s:
            spec["kmesh"] = tuple(int(x) for x in s["kmesh"].replace(",", " ").split())
        if "kind" in s:
            spec["kind"] = s["kind"]
    return spec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--potcar-dir", default="/public/home/<user>/software/vasp_pseudopotentials")
    ap.add_argument("--submit", action="store_true", help="生成后立即 sbatch")
    ap.add_argument("--refdir", default="convex_hull_references")
    ap.add_argument("--only", default="", help="只处理这些相（逗号分隔），默认全部")
    ap.add_argument("--force", action="store_true", help="强制重写并重提交已完成（已收敛）的相")
    args = ap.parse_args()
    if args.submit:
        ap.error("直接提交已禁用：请用 autozt 公共提交入口，避免参考相重复作业。")

    only = {x.strip() for x in args.only.split(",") if x.strip()} if args.only else None
    refdir = Path(args.refdir)
    if not refdir.is_dir():
        raise SystemExit("[错误] 找不到 %s（请先把 POSCAR 放进 <refdir>/<相>/）" % refdir)

    conf = D.load_stepconf()
    refconf = load_conf_for(refdir, conf)
    phases = discover_phases(refdir)
    if not phases:
        raise SystemExit("[错误] %s 下没有 <相>/POSCAR" % refdir)

    soc = str(conf.get("SOC", "1")).lower() not in ("0", "false", ".false.", "no", "off")
    ncomp = 3 if soc else 1
    lsorbit = "LSORBIT = .TRUE." if soc else "# LSORBIT off (SOC=0)"
    ggacompat = "GGA_COMPAT = .FALSE." if soc else "# GGA_COMPAT off (SOC=0)"

    n = 0
    for name, info in phases.items():
        if only is not None and name not in only:
            continue
        d = refdir / name
        spec = resolve_spec(name, info, refconf)
        outcar = d / "OUTCAR"
        if outcar.exists() and not args.force:
            try:
                if "reached required accuracy" in outcar.read_text(errors="ignore"):
                    print("[跳过] %s 已完成且收敛（--force 重算）" % name)
                    continue
            except OSError:
                pass
        incar = (INCAR.replace("{{name}}", name)
                      .replace("{{ENCUT}}", str(conf.get("ENCUT", 370)))
                      .replace("{{ismear}}", str(spec["ismear"]))
                      .replace("{{sigma}}", str(spec["sigma"]))
                      .replace("{{ISPIN}}", str(conf.get("ISPIN", "1")))
                      .replace("{{NCORE}}", str(conf.get("NCORE", "8")))
                      .replace("{{KPAR}}", str(conf.get("KPAR", "2")))
                      .replace("{{LSORBIT_LINE}}", lsorbit)
                      .replace("{{GGA_COMPAT_LINE}}", ggacompat)
                      .replace("{{magmom}}", "%d*0" % (ncomp * info["natoms"])))
        (d / "INCAR").write_text(incar, encoding="utf-8")
        D.write_kpoints(str(d / "KPOINTS"), spec["kmesh"])
        D.assemble_potcar(info["order"], args.potcar_dir, out_path=str(d / "POTCAR"))
        D.render_submit(D.find_submit_tpl(soc), str(d / "submit.sh"), "ref_" + name)
        print("[OK] %-12s %2d 原子 %-10s %-12s k=%s" %
              (name, info["natoms"], ",".join(info["order"]), spec["kind"], spec["kmesh"]))
        n += 1
    print("共生成 %d 个参考相输入" % n)


if __name__ == "__main__":
    main()
