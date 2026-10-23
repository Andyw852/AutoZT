#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2.15_discriminant.py —— 体系判别分支：独立密网格静态

为什么需要这一步（三条职责，缺一不可）
------------------------------------------------------------------
1) 判定体系是 SEMICONDUCTOR / SEMIMETAL / METAL。判据必须是**带指标法**
   min_k E[nocc](k) - max_k E[nocc-1](k)（nocc = NELECT/2），
   **不能用 OUTCAR 的 "fundamental gap" 行**。实测同一份 OUTCAR
   （Mg4C60_relax3 step2.1_static，V=1228.1221，4x4x2，NBANDS=324）：
       带指标法   = -0.0765 eV   <- 真值，负值 = 带重叠
       占据数法(VASP 原生) = +0.0512 eV  <- 有带重叠时该估计量失效
   机理：Gamma 点上只有 247 个带满占据、1 个带部分占据，E[248]@Gamma 被占据
   阈值划成"空态"，而 E[249]@(0.5,0,0.5) 更低却是占据的。凡读那一行的下游
   脚本都会把半金属误报成有隙。

2) **给 step2.3 杂化步选 ALGO。** 所以本步必须排在杂化**之前**（seq 2.15，
   after: step2.1_static）：Damped 对有部分占据的体系更稳，All（共轭梯度
   全能带）在有隙体系上收敛更快但在金属上容易出问题。
       SEMICONDUCTOR      -> ALGO = All
       SEMIMETAL / METAL  -> ALGO = Damped
   gen_step4_HSE.py 读不到 discriminant.json 时退回现状 Damped，不改默认行为。

3) 顺带给出**密网格上的正确 E_F**。step2.1/2.2 的 E_F 是在粗网格（如 4x4x2）
   上收敛的偏值，画能带时对齐位置是错的；本步的 E_F 写进 discriminant.json
   的 efermi_eV，供 band plot 对齐用。本步同时留下 EIGENVAL 与 CHGCAR
   （step2.1 的 INCAR 本来就 LCHARG=.TRUE.）。

为什么必须换网格
------------------------------------------------------------------
step2.2_pbe 复用的是 step2.1_static 的网格，用它判"网格够不够"等于自证。
本步按 KSPACING 重新生成，并**强制严格比 step2.1 的网格更密**（不密就报错），
即规格里的"never reusing step2.2's mesh"。

复用而非重写
------------------------------------------------------------------
INCAR / POSCAR / POTCAR 直接从 step2.1_static 复制 —— 那样泛函（pbe-d3 /
pbesol / pbe）、磁性（ISPIN/MAGMOM）、维度、ENCUT 等与 step2.1 逐项一致，
本步唯一的差别就是 k 网格。不重复实现 gen_step2_static.py 里那套继承逻辑。
submit.sh 按**本步所在集群**的 submit_std 模板渲染，再套本步 step.conf 的 [submit]（[0075]；
以前从 step2.1_static 复制，本步用 step_hpc 换集群时会带着 S2.1 集群的分区）。
"""

import argparse
import math
import re
import shutil
import sys
from pathlib import Path

STEP21_DIR = "step2_bandgap/step2.1_static"
OUTDIR_NAME = "step2_bandgap/step2.15_discriminant"
STEP_LABEL = "S2.15_discr"

# 目标 k 间距（A^-1，2pi 约定）。0.02 明显密于 step2 默认的 KSPACING=0.03。
KSPACING_MAX = 0.02
# 各维下限：3D 至少 8x8x6；2D 面内至少 8x8、kz=1
FLOOR_3D = (8, 8, 6)
FLOOR_2D = (8, 8, 1)
# 相对 step2.1 网格的最小加密倍数
MIN_FACTOR = 2


def log(m):
    print(m, file=sys.stderr, flush=True)


def read_poscar(p):
    L = [x for x in Path(p).read_text(errors="ignore").replace("\r", "").split("\n")]
    scale = float(L[1].split()[0])
    lat = [[float(x) * scale for x in L[2 + i].split()[:3]] for i in range(3)]
    return lat


def recip_lengths(lat):
    a1, a2, a3 = lat
    cr = [a1[1]*a2[2]-a1[2]*a2[1], a1[2]*a2[0]-a1[0]*a2[2], a1[0]*a2[1]-a1[1]*a2[0]]
    V = abs(a3[0]*cr[0] + a3[1]*cr[1] + a3[2]*cr[2])
    out = []
    for (x, y) in ((a2, a3), (a3, a1), (a1, a2)):
        b = [2*math.pi*(x[1]*y[2]-x[2]*y[1])/V,
             2*math.pi*(x[2]*y[0]-x[0]*y[2])/V,
             2*math.pi*(x[0]*y[1]-x[1]*y[0])/V]
        out.append(math.sqrt(sum(t*t for t in b)))
    return out


def parse_kpoints(p):
    """读 VASP 格式 KPOINTS，返回 (mesh[3], shift[3], centered_str)。"""
    txt = Path(p).read_text(errors="ignore").replace("\r", "").split("\n")
    # 第 2 行是 0 -> 自动网格；第 3 行 Gamma/Monkhorst；第 4 行 mesh；第 5 行 shift
    scheme = None
    for ln in txt[:6]:
        if ln.strip().lower().startswith(("g", "m")):
            scheme = ln.strip()
            break
    mesh = None
    shift = [0.0, 0.0, 0.0]
    for i, ln in enumerate(txt):
        t = ln.split()
        if len(t) == 3 and all(re.fullmatch(r"-?\d+(\.0+)?", x) for x in t) and mesh is None and i >= 3:
            mesh = [int(float(x)) for x in t]
            if i + 1 < len(txt) and len(txt[i+1].split()) == 3:
                shift = [float(x) for x in txt[i+1].split()]
            break
    if mesh is None:
        mesh = [1, 1, 1]
    return mesh, shift, (scheme or "Gamma")


def detect_dim(cwd, mesh):
    """2D 判定：workflow_method.txt 的 DIM= 优先，其次 KPOINTS 的 kz==1。"""
    wm = cwd / "workflow_method.txt"
    if wm.is_file():
        m = re.search(r"^\s*DIM\s*=\s*(\S+)", wm.read_text(errors="ignore"), re.M)
        if m:
            v = m.group(1).upper()
            if "2D" in v:
                return 2
            if "3D" in v:
                return 3
    return 2 if mesh[2] == 1 else 3


def build_dense_mesh(nb, base, dim):
    floor = FLOOR_2D if dim == 2 else FLOOR_3D
    new = []
    for i in range(3):
        if dim == 2 and i == 2:
            new.append(1)
            continue
        by_spacing = int(math.ceil(nb[i] / KSPACING_MAX - 1e-9))
        by_factor = MIN_FACTOR * base[i]
        new.append(max(by_spacing, by_factor, floor[i]))
    return new


def _apply_submit_conf(sub):
    """[0075] 本步 step.conf 的 [submit]（partition/qos/time/ntasks_per_node…）套到 submit.sh；没有 stepconf.py -> 不改。"""
    for d in (Path(__file__).resolve().parent, Path.cwd()):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    try:
        import stepconf
    except ImportError:
        return []
    # used_incar=True：本步的 INCAR 原样抄 S2.1（S2.1 已套过 [incar*]），不该报"写了没人读"
    return stepconf.apply_submit(str(sub), stepconf.read_submit(stepconf.CONF_NAME, used_incar=True)) or []


def _cascade_on_rerun(step):
    """[patch_rerun_cascade V153] 本步目录是新建的（rerun）-> 下游步骤整目录归档、等本步跑完后重新生成。"""
    try:
        for _up in Path(__file__).resolve().parents[:3]:
            if (_up / "ke_common.py").is_file() and str(_up) not in sys.path:
                sys.path.insert(0, str(_up))
        import ke_common as _kc
    except ImportError:
        print("[WARN] 没有 ke_common.py：%s 重算不会让下游自动重排（skill.yaml 的 gen_need 漏了它？）" % step)
        return []
    return _kc.cascade_on_rerun(Path.cwd(), step)


def main():
    ap = argparse.ArgumentParser(description="生成体系判别步的独立密网格静态输入")
    ap.add_argument("--force", action="store_true", help="即使已存在也重建")
    ap.add_argument("--jobname", default=None)
    args = ap.parse_args()

    cwd = Path.cwd()
    s21 = cwd / STEP21_DIR
    out = cwd / OUTDIR_NAME
    _cascade_on_rerun(OUTDIR_NAME)
    if not s21.is_dir():
        sys.exit("[ERROR] 缺 %s —— step2.1_static 没跑完？" % STEP21_DIR)
    for fn in ("INCAR", "POSCAR", "POTCAR"):
        if not (s21 / fn).is_file():
            sys.exit("[ERROR] 缺 %s/%s" % (STEP21_DIR, fn))
    if out.is_dir() and any(out.iterdir()) and not args.force:
        log("[..] %s 已存在，原样保留（要重建加 --force）" % OUTDIR_NAME)
    out.mkdir(parents=True, exist_ok=True)

    # --- 网格：独立生成，强制严格更密 ---
    if (s21 / "KPOINTS").is_file():
        base, shift, scheme = parse_kpoints(s21 / "KPOINTS")
    else:
        base, shift, scheme = [1, 1, 1], [0.0, 0.0, 0.0], "Gamma"
    dim = detect_dim(cwd, base)
    lat = read_poscar(s21 / "POSCAR")
    nb = recip_lengths(lat)
    new = build_dense_mesh(nb, base, dim)
    if new == base or all(new[i] <= base[i] for i in range(3)):
        sys.exit("[ERROR] 生成的网格 %s 未能严格密于 step2.1 的 %s —— "
                 "判别步必须用独立且更密的网格，否则等于自证。" % (new, base))
    log("[OK] |b| = %.4f %.4f %.4f A^-1，dim=%dD" % (nb[0], nb[1], nb[2], dim))
    log("[OK] step2.1 网格 %s -> 判别网格 %s（KSPACING<=%s，下限 %s，加密>=%dx）"
        % (base, new, KSPACING_MAX, FLOOR_2D if dim == 2 else FLOOR_3D, MIN_FACTOR))

    # --- 复制 step2.1 的输入（保证除网格外逐项一致）---
    for fn in ("INCAR", "POSCAR", "POTCAR"):
        shutil.copyfile(s21 / fn, out / fn)
    (out / "KPOINTS").write_text(
        "system-discriminant: independent denser static, KSPACING<=%s "
        "(|b|=%.4f,%.4f,%.4f) base=%s -> %s\n"
        "0\n%s\n %3d %3d %3d\n%.1f  %.1f  %.1f\n"
        % (KSPACING_MAX, nb[0], nb[1], nb[2], base, new, scheme,
           new[0], new[1], new[2], shift[0], shift[1], shift[2]),
        encoding="utf-8", newline="\n")

    # --- submit.sh：按本步集群的 submit_std 模板渲染 + 本步 step.conf [submit]（[0075]）---
    #   以前从 step2.1 复制：项目 step_hpc 把本步换到别的集群（S2.1 在 hfeshell、本步在 jzzn）时，
    #   复制来的 submit.sh 带着 S2.1 集群的分区 -> "sbatch: invalid partition specified: hfacexclu08"，
    #   step.conf 的 [submit] 也从来没人读。autozt 按本步的集群推 submit_std_*.tpl（gen_need），
    #   同集群时就是 S2.1 用的那份；模板没推过来（旧的 gen_need）才退回复制 step2.1 的。
    jobname = args.jobname or ("%s-ke-dft-cpu-%s" % (cwd.name, STEP_LABEL))
    tpl = next((cwd / n for n in ("submit_std_%dd.tpl" % dim, "submit_std.tpl") if (cwd / n).is_file()), None)
    src_sub = s21 / "submit.sh"
    if tpl is not None:
        txt = tpl.read_text(encoding="utf-8", errors="ignore").replace("{{JOBNAME}}", jobname)
        left = sorted(set(re.findall(r"\{\{(\w+)\}\}", txt)))
        if left:
            sys.exit("[ERROR] %s 里还有没填的占位符 %s：本脚本只填 {{JOBNAME}}，其余请直接写进模板" % (tpl.name, left))
        how = "由本集群模板 %s 渲染" % tpl.name
    elif src_sub.is_file():
        txt = src_sub.read_text(errors="ignore")
        txt = re.sub(r"(#SBATCH\s+--job-name=).*", r"\g<1>" + jobname, txt)
        how = "由 step2.1_static 复制（没有 submit_std 模板；本步与 S2.1 不在同一集群时分区会错）"
    else:
        txt = None
        log("[WARN] 没有 submit_std 模板，也没有 step2.1_static/submit.sh，请自行提供 submit.sh")
    if txt is not None:
        (out / "submit.sh").write_text(txt, encoding="utf-8", newline="\n")
        (out / "submit.sh").chmod(0o755)
        changed = _apply_submit_conf(out / "submit.sh")
        log("[OK] submit.sh %s，job-name=%s%s" % (how, jobname, ("；step.conf [submit]：%s" % ", ".join(changed))
                                                if changed else ""))

    log("[DONE] %s 就绪：独立密网格 %s（step2.1 是 %s）" % (OUTDIR_NAME, new, base))
    log("       跑完后 check_discriminant.py 会写 discriminant.json")


if __name__ == "__main__":
    main()
