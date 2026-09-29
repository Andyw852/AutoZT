# -*- coding: utf-8 -*-
"""ke_common.py —— ke-dft-cpu 技能新步骤（uniform / dfpt / deform / amset）的公共工具。

只依赖标准库 + dim_common（同目录，setup 已放好）。故意不碰 pymatgen，
让 gen 脚本在登录节点用系统 python 就能跑。VASPKIT 负责 KPOINTS/POTCAR。

放置：由 skill.yaml 的 gen_need 列出，随每个用它的步骤推到材料目录。
      ——因此本文件要复制进每个用到它的步骤源目录（step3_uniform、
        step5_dielect、step7_deform）。见 setup_ke.sh。
"""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dim_common import (adaptive_parallel_tags, detect_dimension,  # noqa: E402
                        force_kz1, validate_poscar, VACUUM_MIN)

METHOD_FILE = "workflow_method.txt"


# --------------------------------------------------------------------------
# 维度
# --------------------------------------------------------------------------
def resolve_dim_for(poscar: Path, dimension="auto", vacuum_min=VACUUM_MIN):
    """返回 (dim, vac_axis)。dim ∈ {'2d','3d'}；vac_axis 仅 2D 有意义。"""
    mode = str(dimension).lower()
    if mode in ("2d", "3d"):
        return mode, (2 if mode == "2d" else None)
    dim, axis, vacs = detect_dimension(poscar, vacuum_min)
    if dim == "2d" and axis != 2:
        sys.exit("[ERROR] 检测到 2D 但真空不在 c 轴（在 %d 轴）。请把结构旋转成"
                 "真空沿第 3 个晶格矢量再重跑。" % axis)
    return dim, (axis if dim == "2d" else None)


def read_method_dim(method_file: Path):
    """从上一步的 workflow_method.txt 读 DIM=2D/3D（有就返回 '2d'/'3d'，无返回 None）。"""
    if not method_file.is_file():
        return None
    for ln in method_file.read_text(errors="ignore").splitlines():
        if ln.strip().upper().startswith("DIM="):
            v = ln.split("=", 1)[1].strip().lower()
            if v in ("2d", "3d"):
                return v
    return None


def write_method(path: Path, dim: str, note: str, func: str = None):
    # patch_ke_dag：多记一行 FUNC=，让本步的泛函能被更下游的步骤继续继承
    lines = ["DIM=%s" % dim.upper()]
    if func:
        lines.append("FUNC=%s" % func)
    lines.append("# %s" % note)
    path.write_text("\n".join(lines) + "\n",
                    encoding="utf-8", newline="\n")


# --------------------------------------------------------------------------
# VASPKIT：KPOINTS / POTCAR / ENCUT
# --------------------------------------------------------------------------
def vaspkit_kpoints(outdir: Path, kscheme="2", kspacing="0.03",
                    exe="vaspkit", dim="3d", vac_axis=2, force_kz1_2d=True):
    """VASPKIT 1→102→scheme→kspacing 生成 KPOINTS；2D 把真空方向细分压回 1。"""
    print("[..] VASPKIT KPOINTS：1 -> 102 -> %s -> %s" % (kscheme, kspacing))
    subprocess.run([exe], input="1\n102\n%s\n%s\n" % (kscheme, kspacing),
                   text=True, cwd=outdir, check=True)
    if dim == "2d" and force_kz1_2d:
        changed, note = force_kz1(outdir / "KPOINTS", axis=vac_axis if vac_axis is not None else 2)
        print("[%s] 2D KPOINTS 真空方向细分：%s" % ("OK" if changed else "..", note))


# --------------------------------------------------------------------------
# [KALIGN-2026-09-23] Γ 心网格的高对称点对齐 —— S3 / S3b / S7 的唯一实现
#   Γ 心网格点 = n_i/N_i。带边高对称点要落在网格上：
#     六方 K=(1/3,1/3)            → 面内 N % 3 == 0
#     正方/矩形 X、M=(1/2,…)      → N % 2 == 0
#     菱面体 L/F/Z、立方 X/L 等   → N % 2 == 0（3D 的 K/H 类 1/3 点另需 3）
#   2D：面内取 6 的倍数一并覆盖（与 9-23 S3 修复同口径）。
#   3D：默认 off（不动存量 3D 项目：225 合金等已跑的 S3/S7 网格）；
#       even = 偶数（菱面体/立方/四方够用），6 = 六方 3D（a、b 取 6 的倍数，c 取偶数）。
#   只向上取整，不放松密度。前提：KSCHEME=2（Γ 心）；MP 偶数网格不含 Γ/K，本函数无意义。
#   局限：只保证高对称点。Q 谷等一般 k 点带边不在此列（见 step3c_uniform_offgrid）。
# --------------------------------------------------------------------------
def align_kgrid(need, dim, axes, align3d="off", label="", quiet=False):
    """返回对齐后的新网格 list；有改动时打印 [WARN]。need 不被原地修改。"""
    out = list(need)
    if dim == "2d":
        mult = 6
    else:
        mode = str(align3d or "off").strip().lower()
        mult = {"off": 1, "none": 1, "even": 2, "2": 2, "6": 6}.get(mode)
        if mult is None:
            raise SystemExit("[ERROR] KALIGN_3D=%r 不认识（off | even | 6）" % align3d)
    if mult <= 1:
        return out
    for i in axes:
        # 六方 3D（KALIGN_3D=6）：K/H 的 1/3 只在 a、b 面内；c 轴 A=(0,0,1/2) 只需偶数。
        #   约定 c 为第 3 个晶格矢量（与 pymatgen/vaspkit 六方标准胞一致）。
        m = 2 if (dim != "2d" and mult == 6 and i == 2) else mult
        out[i] = int(-(-int(out[i]) // m) * m)
    if out != list(need) and not quiet:
        print("[WARN] %s高对称点对齐（Γ 心，%s 取 %d 的倍数）：%s -> %s"
              % (label + " " if label else "", dim.upper(), mult,
                 "x".join(str(x) for x in need), "x".join(str(x) for x in out)))
    return out


def read_kpoints_mesh(kpoints: Path):
    """读 Γ/MP 自动网格 KPOINTS 第 4 行的 3 个分割数；读不到返回 None。"""
    try:
        v = [int(x) for x in Path(kpoints).read_text(errors="ignore").splitlines()[3].split()[:3]]
        return v if len(v) == 3 else None
    except (OSError, IndexError, ValueError):
        return None


# --------------------------------------------------------------------------
# [patch_stale_grid-2026-09-23] 网格变了 -> 归档旧产物（S3 / S3b 用）
#   ck_wavecar 只查 WAVECAR 存不存在、够不够大，**不比对网格**。所以改了网格之后，
#   旧网格算完的目录仍被判「已完成」，auto/start 不会重算 —— 输入是新网格、产物是旧网格，
#   下游静默用错。这里在【新 KPOINTS 已就绪且全部护栏通过之后】把旧产物改名归档
#   （不删除），使 ck_wavecar 判不通过、步骤重新进队。
#   输入（KPOINTS/INCAR/POSCAR/POTCAR/submit.sh）不在归档集合内，绝不触碰。
#   （S7 的 gen 里有一份等价实现 patch_stale_grid，因 fanout 要逐子目录归档而未合并。）
# --------------------------------------------------------------------------
STALE_GRID_OUTPUTS = (
    "OUTCAR", "OSZICAR", "vasprun.xml", "CONTCAR", "CHGCAR", "CHG", "EIGENVAL",
    "DOSCAR", "PROCAR", "WAVECAR", "IBZKPT", "REPORT", "PCDAT", "XDATCAR",
    "LOCPOT", "ELFCAR",
)


def archive_stale_grid(outdir, old_mesh, new_mesh):
    """旧网格 != 新网格 时把旧产物改名 *.stale-grid-<旧>；返回 (归档数, 旧标签)。

    old_mesh 为 None（首次生成 / 读不到）或与 new_mesh 相同 -> 不动，返回 (0, None)。
    """
    if not old_mesh or list(old_mesh[:3]) == list(new_mesh[:3]):
        return 0, None
    _tag = "x".join(str(x) for x in old_mesh[:3])
    _out = Path(outdir)
    _n = 0
    for _nm in STALE_GRID_OUTPUTS:
        _f = _out / _nm
        if _f.is_file():
            _f.rename(_out / (_nm + ".stale-grid-" + _tag))
            _n += 1
    return _n, _tag


# --------------------------------------------------------------------------
# [patch_stale_input-2026-09-28] 输入变了 -> 旧产物归档；上游重算 -> 下游完成标记失效
#   patch_stale_grid 只管"网格变了"。但 S3 的 POSCAR 也会变（align_origin 平移原点、
#   step1 重弛豫），INCAR 也会变（step.conf 的 [incar]）—— 旧 WAVECAR 还在，ck_wavecar
#   照样判完成，下游静默用旧波函数（MoS₂ 09-27 的 IBZ h5 就是这样"半新半旧"的）。
#   另外 autozt 的完成判据只看各步**自己**的 marker：S3 重算了，S4 的 wavefunction.h5
#   还在 -> S4 仍判完成，S8 继续用旧 h5。这里给两件事各一个工具：
#     snapshot_inputs / inputs_changed   —— gen 开头拍快照、结尾比对（按内容语义比，
#                                          注释/空白/键顺序不算变化）；
#     invalidate_downstream              —— 上游确实要重算时，把下游的完成标记改名
#                                          *.stale-upstream-<上游>-<时间>（不删除）。
# --------------------------------------------------------------------------
STALE_INPUT_FILES = ("POSCAR", "INCAR", "KPOINTS")


def _poscar_key(text):
    """POSCAR -> (晶格 9 数, 元素行, 个数行, 坐标)，数值四舍五入到 1e-6（忽略注释行）。"""
    ln = [x for x in text.splitlines()]
    try:
        s = float(ln[1].split()[0])
        lat = tuple(round(float(v) * s, 6) for r in ln[2:5] for v in r.split()[:3])
        i = 5
        species = ()
        if not ln[i].split()[0].lstrip("-").replace(".", "").isdigit():
            species = tuple(ln[i].split())
            i += 1
        counts = tuple(int(x) for x in ln[i].split())
        i += 1
        if ln[i].strip()[:1] in ("S", "s"):
            i += 1
        mode = ln[i].strip()[:1].upper()
        n = sum(counts)
        xyz = tuple(round(float(v), 6) for r in ln[i + 1:i + 1 + n] for v in r.split()[:3])
        return (lat, species, counts, mode, xyz)
    except (IndexError, ValueError):
        return text.strip()


# 只影响并行/性能/标签、不影响物理结果的 INCAR 键：换宿主机（GPU/CPU）时它们会变，不能算"输入变了"
_INCAR_NONPHYS = {"SYSTEM", "NCORE", "NPAR", "KPAR", "NSIM", "LPLANE", "LSCALU", "NBLOCK",
                  "KBLOCK", "NWRITE"}


def _input_key(name, text):
    if name == "INCAR":
        return tuple(sorted((k, " ".join(v.split())) for k, v in parse_incar(text).items()
                            if k not in _INCAR_NONPHYS))
    if name == "POSCAR":
        return _poscar_key(text)
    if name == "KPOINTS":
        return tuple(" ".join(x.split()) for x in text.splitlines()[1:4])
    return text


def snapshot_inputs(outdir, names=STALE_INPUT_FILES):
    """gen 开头调：记下已有输入的语义快照（没有的文件记 None）。"""
    snap = {}
    for n in names:
        p = Path(outdir) / n
        snap[n] = _input_key(n, p.read_text(errors="ignore")) if p.is_file() else None
    return snap


def inputs_changed(outdir, snap):
    """gen 结尾调：与快照相比语义变了的输入名列表（原来就没有的文件不算变化）。"""
    out = []
    for n, old in (snap or {}).items():
        if old is None:
            continue
        p = Path(outdir) / n
        new = _input_key(n, p.read_text(errors="ignore")) if p.is_file() else None
        if new != old:
            out.append(n)
    return out


def archive_stale_outputs(outdir, tag):
    """把 STALE_GRID_OUTPUTS 里存在的旧产物改名 *.<tag>；返回归档个数。"""
    _out = Path(outdir)
    n = 0
    for nm in STALE_GRID_OUTPUTS:
        f = _out / nm
        if f.is_file():
            f.rename(_out / (nm + "." + tag))
            n += 1
    return n


def finalize_stale_inputs(cwd, outdir, step, snap, n_grid_archived=0):
    """S3/S3b gen 末尾调：输入语义变了 -> 归档旧产物；归档过（网格或输入）-> 下游失效。

    返回 (变了的输入, 归档个数, 失效的下游)。
    """
    import time as _t
    changed = inputs_changed(outdir, snap)
    n = 0
    if changed:
        n = archive_stale_outputs(outdir, "stale-input-" + _t.strftime("%Y%m%d%H%M%S"))
        if n:
            print("[..] patch_stale_input：%s 变了（%s）-> 归档旧产物 %d 个，本步会重算"
                  % (step, "/".join(changed), n))
    down = []
    if n or n_grid_archived:
        down = invalidate_downstream(cwd, step, "输入变了：%s" % ("/".join(changed) or "网格"))
    return changed, n, down


# 下游关系（与 skill.yaml 的 needs 对应）。mode="link"：下游靠软链引用上游产物
# （S4 链 S3 的 WAVECAR/vasprun；S8/S8.4 链 S4/S4b 的 h5 与 S3/S3b 的 vasprun）——
# 只有软链确实指向这个上游时才失效；mode="always"：下游直接读上游目录，一律失效。
DOWNSTREAM = {
    "step3_uniform": [("step4_wave", "link"), ("step8_amset", "link"),
                      ("step8.4_amset2d", "link"), ("step8.1_boltztrap", "always"),
                      ("step8.2_dpt", "always")],
    "step3b_uniform_full": [("step4b_wave_full", "link"), ("step8_amset", "link"),
                            ("step8.4_amset2d", "link")],
    "step4_wave": [("step8_amset", "link"), ("step8.4_amset2d", "link")],
    "step4b_wave_full": [("step8_amset", "link"), ("step8.4_amset2d", "link")],
    "step8_amset": [("step8.3_output", "always")],
    "step8.1_boltztrap": [("step8.3_output", "always")],
    "step8.2_dpt": [("step8.1_boltztrap", "always"), ("step8.3_output", "always")],
}
DONE_MARKERS = {
    "step4_wave": ("wavefunction.h5",),
    "step4b_wave_full": ("wavefunction.h5",),
    "step8_amset": ("transport.json",),
    "step8.4_amset2d": ("transport.json", "intrinsic_transport.json"),
    "step8.1_boltztrap": ("boltztrap_crta.json",),
    "step8.2_dpt": ("dpt_result.json",),
    "step8.3_output": ("comparison_300K.png",),
}


def _links_into(d, upstream_dir):
    """d 里有没有软链指向 upstream_dir（解析后比较）。"""
    try:
        up = Path(upstream_dir).resolve()
        for f in Path(d).iterdir():
            if f.is_symlink():
                try:
                    tgt = f.resolve()
                except OSError:
                    continue
                if up == tgt.parent or up in tgt.parents:
                    return True
    except OSError:
        pass
    return False


def invalidate_downstream(cwd, step, reason, _seen=None):
    """上游 step 要重算：把下游完成标记改名 *.stale-upstream-<step>-<时间>，递归传递。

    返回 [(下游步骤, 改名的文件), ...]。只改名不删除；下游目录不存在/没有标记时什么都不做。
    """
    import time as _t
    cwd = Path(cwd)
    seen = _seen if _seen is not None else set()
    done = []
    tag = "stale-upstream-%s-%s" % (step.replace("/", "_"), _t.strftime("%Y%m%d%H%M%S"))
    for down, mode in DOWNSTREAM.get(step, ()):
        if down in seen:
            continue
        d = cwd / down
        if not d.is_dir():
            continue
        if mode == "link" and not _links_into(d, cwd / step):
            continue
        seen.add(down)
        hit = False
        for m in DONE_MARKERS.get(down, ()):
            f = d / m
            if f.is_file() or f.is_symlink():
                f.rename(d / (m + "." + tag))
                done.append((down, m))
                hit = True
        if hit:
            print("[WARN] 上游 %s 重算（%s）-> 下游 %s 的完成标记已归档（*.%s），会重新排队"
                  % (step, reason, down, tag))
        done += invalidate_downstream(cwd, down, "上游 %s 失效" % step, seen)
    return done


# --------------------------------------------------------------------------
# [patch_mem_guard-2026-09-28] AMSET 作业级内存粗估（第 6 项）
#   AMSET 自报的 max memory 只统计主进程（memory_profiler include_children=False），作业级
#   （sacct MaxRSS，含全部 worker）是它的 2–4 倍；jzzn 的 Slurm 又不按内存调度（V113 §1-2）。
#   标定只有两个点（都是 nworkers=24，作业级 MaxRSS / 最终插值网格点数）：
#     GaAs S8 3D IBZ 109³ = 1.30M 点 -> 150.9 GiB  => 116 GiB / 百万点（上界）
#     MoS₂ S8.4 unity f=61 263×263×21 = 1.45M 点 -> 91.1 GiB  => 62.8 GiB / 百万点（下界）
#   带数、谷结构不同会让系数在两者之间移动，所以给**区间**，不给单值；nworkers 越少越省。
# --------------------------------------------------------------------------
MEM_GIB_PER_MPT = (62.8, 116.1)
MEM_CAL_NOTE = ("标定：MoS₂ S8.4 1.45M 点 91.1 GiB / GaAs S8 1.30M 点 150.9 GiB"
                "（作业级 MaxRSS，nworkers=24，V113）")


def estimate_amset_mem_gib(mesh):
    """最终插值网格 -> (下界, 上界) GiB。"""
    pts = 1.0
    for x in mesh:
        pts *= int(x)
    return (MEM_GIB_PER_MPT[0] * pts / 1e6, MEM_GIB_PER_MPT[1] * pts / 1e6)


def mem_guard(mesh, warn_gib=300, limit_gib=None, label=""):
    """打印内存粗估；返回 (lo, hi, level)，level ∈ ok/warn/error。

    warn：上界 > warn_gib（建议独占节点 / 降 nworkers / 降 MESH_MAX）；
    error：给了 limit_gib 且**下界**都超过它（肯定放不下）。
    """
    lo, hi = estimate_amset_mem_gib(mesh)
    msg = ("[..] %s内存粗估：最终网格 %s -> 作业级约 %.0f–%.0f GiB（%s）"
           % (label + " " if label else "", "x".join(str(int(x)) for x in mesh), lo, hi,
              MEM_CAL_NOTE))
    print(msg)
    if limit_gib and lo > float(limit_gib):
        print("[ERROR] 内存粗估下界 %.0f GiB 已超过 MEM_LIMIT_GIB=%s —— 放不下。"
              "处理：降 MESH_MAX / 显式降 INTERPOLATION_FACTOR / 降 nworkers。" % (lo, limit_gib))
        return lo, hi, "error"
    if warn_gib and hi > float(warn_gib):
        print("[WARN] 内存粗估上界 %.0f GiB > MEM_WARN_GIB=%s。jzzn 的 Slurm 不按内存调度，"
              "与别的作业同节点会 OOM：建议独占节点（step.conf 写 MEM_AUTO_EXCLUSIVE = true "
              "自动加 #SBATCH --exclusive）、或降 nworkers / MESH_MAX。" % (hi, warn_gib))
        return lo, hi, "warn"
    return lo, hi, "ok"


def add_exclusive(submit_path):
    """在 submit.sh 的最后一个 #SBATCH 行后补 #SBATCH --exclusive（已有则不动）。返回是否改动。"""
    p = Path(submit_path)
    if not p.is_file():
        return False
    lines = p.read_text(encoding="utf-8").splitlines()
    if any(re.match(r"^#SBATCH\s+--exclusive\b", x) for x in lines):
        return False
    idx = [i for i, x in enumerate(lines) if x.startswith("#SBATCH")]
    if not idx:
        return False
    lines.insert(idx[-1] + 1, "#SBATCH --exclusive")
    p.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return True


def vaspkit_potcar(outdir: Path, exe="vaspkit"):
    if (outdir / "POTCAR").exists():
        print("[OK] POTCAR 已存在，跳过")
        return
    print("[..] VASPKIT POTCAR：1 -> 103")
    subprocess.run([exe], input="1\n103\n", text=True, cwd=outdir, check=True)


def encut_from_potcar(potcar: Path, factor=1.5, fallback=300):
    vals = []
    for ln in potcar.read_text(errors="ignore").splitlines():
        m = re.search(r"ENMAX\s*=\s*([\d.]+)", ln)
        if m:
            vals.append(float(m.group(1)))
    if not vals:
        return int(fallback)
    import math
    return int(math.ceil(max(vals) * factor / 10.0) * 10)


# --------------------------------------------------------------------------
# INCAR 模板渲染 + 键改写
# --------------------------------------------------------------------------
def render_tpl(tpl_path: Path, subs: dict, out_path: Path):
    """把模板里的 {{KEY}} 占位符替换成 subs[KEY]，写出。"""
    text = tpl_path.read_text(encoding="utf-8")
    for k, v in subs.items():
        text = text.replace("{{%s}}" % k, str(v))
    left = re.findall(r"\{\{([A-Z_]+)\}\}", text)
    if left:
        sys.exit("[ERROR] 模板 %s 还有未填占位符：%s" % (tpl_path.name, ", ".join(set(left))))
    out_path.write_text(text, encoding="utf-8", newline="\n")
    print("[OK] %s" % out_path.name)


# ---- [PATCH-IONRELAX] 应变配对反解（gen 与 step7b 共用，单一真源）----
def read_lattice_matrix(poscar):
    """读 POSCAR 晶格矩阵（3x3，含 scale）。纯 numpy，不碰 pymatgen。"""
    import numpy as np
    ln = Path(poscar).read_text().splitlines()
    scale = float(ln[1].split()[0])
    rows = [[float(x) for x in ln[i].split()[:3]] for i in (2, 3, 4)]
    return np.array(rows) * scale


def resolve_strain_pairs(out):
    """反解 xx±/yy± 形变配对（纯 numpy，与 step7b 的 _be_code 同源同公式）。

    形变梯度 F = (und_lat⁻¹ · def_lat)ᵀ（列矢量约定，与 amset calculate_deformation
    一致）；Green-Lagrange 应变 E = (FᵀF - I)/2；对角元最大者即主形变轴。
    返回 (pairs, strain_mag)：pairs={"xx":[plus,minus],"yy":[plus,minus]}（目录名），
    strain_mag={"xx":γxx,"yy":γyy}（Green–Lagrange 对角元）。反解不完整返回 (None, None)。
    gen_step9_deform.py 建 ionrelax/ 与 step7b 找 ionrelax/ 必须走这一个函数，
    否则两边各自反解会错位（gen 建 01/02、step7b 翻 03/04，静默降级）。"""
    import glob
    import numpy as np
    und_lat = read_lattice_matrix(Path(out) / "undeformed" / "POSCAR")
    pairs = {"xx": [None, None], "yy": [None, None]}
    strain_mag = {"xx": 0.0, "yy": 0.0}
    for d in sorted(glob.glob(os.path.join(str(out), "deform-*"))):
        p = os.path.join(d, "POSCAR")
        if not os.path.isfile(p):
            continue
        F = np.transpose(np.dot(np.linalg.inv(und_lat), read_lattice_matrix(p)))
        E = (np.dot(F.T, F) - np.eye(3)) / 2.0
        diag = [float(E[i, i]) for i in range(3)]
        idx = int(np.argmax(np.abs(diag)))
        val = diag[idx]
        if abs(val) < 1e-6:
            continue
        sign = 0 if val > 0 else 1
        if idx == 0:
            pairs["xx"][sign] = os.path.basename(d)
            strain_mag["xx"] = abs(val)
        elif idx == 1:
            pairs["yy"][sign] = os.path.basename(d)
            strain_mag["yy"] = abs(val)
    if None in pairs["xx"] or None in pairs["yy"]:
        return None, None
    return pairs, strain_mag


def apply_parallel_tags(incar_path):
    """按宿主机行级改写 INCAR 的 NCORE/KPAR（渲染后调用，模板里写死的值也能覆盖）。

    判定来自 dim_common.adaptive_parallel_tags()（单一真源）：
      GPU 版 → NCORE=1/KPAR=1；CPU → 不改（返回 False）。
    返回 True 表示已改写，False 表示未动（CPU 或文件缺失）。
    """
    tags = adaptive_parallel_tags()
    if not tags:
        return False
    p = Path(incar_path)
    if not p.is_file():
        return False
    text = p.read_text(encoding="utf-8", errors="ignore")
    out, seen = [], set()
    for ln in text.splitlines():
        m = re.match(r"\s*([A-Za-z_]+)\s*=", ln)
        if m and m.group(1).upper() in tags:
            k = m.group(1).upper()
            out.append("%-8s = %s" % (k, tags[k]))
            seen.add(k)
            continue
        out.append(ln)
    for k, v in tags.items():
        if k not in seen:
            out.append("%-8s = %s" % (k, v))
    p.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    print("[..] 并行参数已按宿主机覆盖：%s" % ", ".join("%s=%s" % kv for kv in tags.items()))
    return True


# ---- [PATCH-UCONS] 自旋 / U 标签继承 ----------------------------------
# step3_uniform / step5_dielect / step7_deform 的 INCAR 是纯模板渲染，模板里
# 没有 ISPIN/MAGMOM/LMAXMIX/LDAU* —— 结构按 +U + 自旋极化弛豫，输运/介电/形变势
# 却按裸 GGA 非自旋极化算，两者不是同一个哈密顿量。这里从上游步骤把这些标签接过来。
#
# ★ 上游 ISPIN!=2 且没有 LDAU 时【完全不动 INCAR】，非磁无 U 体系零改动。
# ★ 关掉：export AUTOZT_KE_NO_SCF_INHERIT=1
SCF_SPIN_TAGS = ("ISPIN", "MAGMOM", "NUPDOWN", "LMAXMIX")
SCF_U_TAGS = ("LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ", "LDAUPRINT")
SCF_SRC_DIRS = ("step2_bandgap/step2.1_static", "step1_opt",
                "step1_std_opt", "step1_PBE_opt")


def find_scf_source(cwd):
    """上游 SCF 的 INCAR：优先 step2.1_static（那里的 MAGMOM 是收敛值），
    回退 step1（初始高自旋猜测）。都没有返回 None。"""
    for rel in SCF_SRC_DIRS:
        p = Path(cwd) / rel / "INCAR"
        if p.is_file():
            return p
    return None


def inherit_scf_tags(incar_path, cwd, with_u=True, label=""):
    """把上游的自旋/U 标签注入已渲染好的 INCAR。返回注入的键列表（空 = 没动）。"""
    if os.environ.get("AUTOZT_KE_NO_SCF_INHERIT"):
        return []
    src = find_scf_source(cwd)
    if src is None:
        return []
    up = parse_incar(src.read_text(encoding="utf-8", errors="ignore"))
    ispin2 = str(up.get("ISPIN", "1")).split()[0] == "2"
    u_on = str(up.get("LDAU", "")).upper().lstrip(".").startswith("T")
    if not ispin2 and not u_on:
        return []                      # 非磁 + 无 U：什么都不用做
    want = list(SCF_SPIN_TAGS) + (list(SCF_U_TAGS) if with_u else [])
    take = {k: up[k] for k in want if k in up}
    if not ispin2:
        take.pop("MAGMOM", None)
        take.pop("NUPDOWN", None)
    if u_on and take.get("LMAXMIX") is None:
        take["LMAXMIX"] = "4"          # 加 U 的 d/f 混合需要 LMAXMIX>=4
    if not take:
        return []
    p = Path(incar_path)
    keep = [ln for ln in p.read_text(encoding="utf-8").splitlines()
            if not re.match(r"\s*(%s)\s*=" % "|".join(want), ln, re.IGNORECASE)]
    keep.append("")
    keep.append("# ---- 以下由 gen 脚本从 %s 继承（自旋/U 必须与弛豫时一致）----"
                % src.parent.as_posix())
    for k in want:
        if k in take:
            keep.append("%-10s = %s" % (k, take[k]))
    if u_on and not with_u:
        keep.append("# 注意：上游带 LDAU，但本步是 DFPT(IBRION=8)，多数 VASP 版本")
        keep.append("#       不支持 LDA+U 的 DFPT —— 这里【故意没有】继承 LDAU*。")
        keep.append("#       要带 U 的介电常数请改走 IBRION=6 有限差分。")
    p.write_text("\n".join(keep) + "\n", encoding="utf-8", newline="\n")
    tag = (" [%s]" % label) if label else ""
    print("[..] 继承上游自旋/U 标签%s：%s" % (tag, ", ".join(sorted(take))))
    if u_on and not with_u:
        print("     （本步是 DFPT，未继承 LDAU*；见 INCAR 末尾注释）")
    return sorted(take)


def parse_incar(text: str):
    d = {}
    for ln in text.splitlines():
        ln = ln.split("#", 1)[0].split("!", 1)[0].strip()
        if "=" in ln:
            k, v = ln.split("=", 1)
            d[k.strip().upper()] = v.strip()
    return d


def incar_text(d: dict, system="calc"):
    out = ["SYSTEM = %s" % system, ""]
    for k, v in d.items():
        if k == "SYSTEM":
            continue
        out.append("%-10s = %s" % (k, v))
    return "\n".join(out) + "\n"


def merge_incar(base: dict, overrides: dict):
    """overrides 里 value 为 None 表示删除该键。"""
    d = dict(base)
    for k, v in overrides.items():
        k = k.upper()
        if v is None:
            d.pop(k, None)
        else:
            d[k] = str(v)
    return d


# --------------------------------------------------------------------------
# 结构接力
# --------------------------------------------------------------------------
def relay_poscar(prev_contcar: Path, dst_poscar: Path, label="上一步"):
    """把上一步 CONTCAR 拷成本步 POSCAR；缺失就报错退出（绝不静默用旧结构）。"""
    if not prev_contcar.is_file():
        sys.exit("[ERROR] %s 的 CONTCAR 不存在：%s\n"
                 "        请确认上一步已完成再生成本步。" % (label, prev_contcar))
    validate_poscar(prev_contcar)
    shutil.copyfile(prev_contcar, dst_poscar)
    print("[OK] POSCAR ← %s" % prev_contcar)


def _sym_ops_max_tau(structure, atol=1e-4):
    """spglib 取 (rotations, translations)，返回 (max|tau|, rots, taus)；取不到返回 (None,None,None)。"""
    try:
        import numpy as np
        import spglib
        ds = spglib.get_symmetry_dataset(
            (structure.lattice.matrix, structure.frac_coords, [x.Z for x in structure.species]),
            symprec=atol)
        if ds is None:
            return None, None, None
        g = (lambda k: ds[k]) if isinstance(ds, dict) else (lambda k: getattr(ds, k))
        rots = np.asarray(g("rotations"), float)
        taus = np.asarray(g("translations"), float)
        tn = taus - np.rint(taus)
        return (float(np.abs(tn).max()) if tn.size else 0.0), rots, tn
    except Exception:
        return None, None, None


# ---- 统一判据（2026-09-28 用户推导 + 12 构型独立实测；第二批 (a) 的**修正版**）------
#   ★ 旧规则（"有反演，或 τ≡0，就可走 IBZ"）**两个方向都错**：
#     · 太保守：GaN P6_3mc 标准原点（6_3 轴过原点）原公式 0/24 错，本可走 IBZ，
#       旧规则却要求全网格 —— 正在排队的 S3b/S4b 其实不必要；
#     · 太危险：原点落在**反演中心**的 Si（Fd-3m）有反演，旧规则放行 IBZ，
#       而原公式 24/48 个操作算错 —— **静默算错**。
#   正确判据是**逐操作**的（R、τ 为 spglib 实空间分数坐标）：
#       非时间反演操作：(R + I)·τ ≡ 0 (mod 1)
#       时间反演操作  ：(R − I)·τ ≡ 0 (mod 1)
#   推导与逐操作实测见 tmp/amset_desym_phase_check.py（12 构型判据与实测 100% 吻合）。
TAU_TOL = 1e-4
# 逐操作判据的容差。与 TAU_TOL 同量级：好操作（经 POSCAR 往返）偏差 <= 4e-6，
# 坏操作 >= 0.02，中间 3 个数量级空档，见 op_phase_exact 的说明。
OP_TOL = 1e-4
# ---- [patch_symprec_unify] 判据与 AMSET 用**同一套**对称操作（V115 §4）----------------
#   AMSET 去对称化走 WavefunctionOverlapCalculator.from_file -> from_coefficients，
#   **不传 symprec**，所以恒用 defaults["symprec"] = 0.01 Å（与 settings.yaml 里的
#   symprec 无关），操作来自 pymatgen SpacegroupAnalyzer(结构, symprec=0.01)
#   （dataset 为 None 时退 angle_tolerance=-1），结构取自 wavefunction.h5（= S3 的结构）。
#   旧判据用 spglib 直接取、容差 1e-4 —— 对略有畸变的结构两边可能取到不同的操作集。
#   现在判据按 AMSET 的取法取操作；另外记下 1e-4 严格容差下的操作数（strict_ops）作诊断。
AMSET_SYMPREC = 0.01
# ---- [patch_desym_fix] AMSET 去对称化相位补丁的开关（V115 §5）----------------------------
#   补丁本体：step8.4_amset2d/amset_desym_fix.py（运行时替换 desymmetrize_coefficients 的
#   相位因子）。开了以后，任何结构、任何原点 IBZ 都精确（模型检验 12/12，见 test_desym_fix.py）。
#   ★ 已验证：MoS₂（V116）+ GaN（V117）对照组判据 + 平移不变性均 PASS，默认开启。
#   取值优先级：step.conf 的 DESYM_FIX（on/off）> gen 时的环境变量 AZ_DESYM_FIX > 本常量。
DESYM_FIX_DEFAULT = True
DESYM_FIX_ENV = "AZ_DESYM_FIX"
DESYM_FIX_PLUGIN = "amset_desym_fix.py"
_ON = ("1", "true", "on", "yes")
_OFF = ("0", "false", "off", "no")


def desym_fix_setting(conf_value=None):
    """返回 (是否开启去对称化相位补丁, 来源说明)。conf_value = step.conf 的 DESYM_FIX 原值。"""
    v = "" if conf_value is None else str(conf_value).strip().lower()
    if v in _ON:
        return True, "step.conf DESYM_FIX=%s" % conf_value
    if v in _OFF:
        return False, "step.conf DESYM_FIX=%s" % conf_value
    if v not in ("", "auto", "none", "default"):
        print("[WARN] DESYM_FIX=%r 不认识（只认 auto/on/off），按 auto 处理" % conf_value)
    e = os.environ.get(DESYM_FIX_ENV, "").strip().lower()
    if e in _ON:
        return True, "环境变量 %s=%s" % (DESYM_FIX_ENV, e)
    if e in _OFF:
        return False, "环境变量 %s=%s" % (DESYM_FIX_ENV, e)
    return bool(DESYM_FIX_DEFAULT), ("出厂默认 DESYM_FIX_DEFAULT=%s（默认开）"
                                     % DESYM_FIX_DEFAULT)


def desym_fix_cmd_prefix(on):
    """拼进 AMSET 作业命令最前面的环境变量（插件运行时只认环境变量）。"""
    return ("export %s=1; " % DESYM_FIX_ENV) if on else ("unset %s; " % DESYM_FIX_ENV)


def incar_is_ncl(incar):
    """INCAR 里 LSORBIT 或 LNONCOLLINEAR 为真 -> 非共线（SOC）。读不到返回 False。"""
    p = Path(incar)
    if not p.is_file():
        return False
    d = parse_incar(p.read_text(encoding="utf-8", errors="ignore"))
    for k in ("LSORBIT", "LNONCOLLINEAR"):
        if str(d.get(k, "")).upper().lstrip(".").startswith("T"):
            return True
    return False


def detect_ncl(cwd, dirs=("step3_uniform", "step3b_uniform_full")):
    """波函数来源步骤（S3/S3b）是不是 SOC 计算。"""
    return any(incar_is_ncl(Path(cwd) / d / "INCAR") for d in dirs)


def _dataset_getter(ds):
    return (lambda k: ds[k]) if isinstance(ds, dict) else (lambda k: getattr(ds, k))


def amset_symmetry_dataset(structure, symprec=AMSET_SYMPREC):
    """按 AMSET get_reciprocal_point_group_operations 的取法拿 spglib dataset 的
    (rotations, translations)（实空间分数坐标）。取不到返回 (None, None, 说明)。"""
    import numpy as np
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        sga = SpacegroupAnalyzer(structure, symprec=symprec)
        ds = sga.get_symmetry_dataset()
        src = "pymatgen SGA(symprec=%g)" % symprec
        if ds is None:
            sga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=-1)
            ds = sga.get_symmetry_dataset()
            src += " angle_tolerance=-1"
    except ImportError:
        import spglib
        cell = (structure.lattice.matrix, structure.frac_coords,
                [x.Z for x in structure.species])
        ds = spglib.get_symmetry_dataset(cell, symprec=symprec, angle_tolerance=5)
        src = "spglib(symprec=%g)" % symprec
        if ds is None:
            ds = spglib.get_symmetry_dataset(cell, symprec=symprec, angle_tolerance=-1)
            src += " angle_tolerance=-1"
    if ds is None:
        return None, None, "dataset 为 None"
    g = _dataset_getter(ds)
    return np.asarray(g("rotations"), float), np.asarray(g("translations"), float), src


def amset_symmetry_ops(structure, symprec=AMSET_SYMPREC):
    """AMSET 去对称化**实际使用**的操作集 [(R_real, tau_real, is_tr), ...] 与来源说明。

    装了 amset（作业内）就直接调它的 get_reciprocal_point_group_operations；
    否则（登录节点）按同样的取法复刻（amset_op_set）。两条路的结果在 12 个构型上逐一相同
    （test_symmetry_gate.test_ops_match_amset）。取不到返回 (None, 说明)。
    """
    try:
        from amset.electronic_structure.symmetry import \
            get_reciprocal_point_group_operations as _gpo
        import numpy as np
        rots, taus, trs = _gpo(structure, symprec=symprec, time_reversal=True)
        ops = []
        for r, t, tr in zip(np.asarray(rots, float), np.asarray(taus, float), trs):
            tr = bool(tr)
            ops.append(((-r if tr else r).T, (-t if tr else t), tr))
        return ops, "amset.get_reciprocal_point_group_operations(symprec=%g)" % symprec
    except ImportError:
        pass
    rots, taus, src = amset_symmetry_dataset(structure, symprec)
    if rots is None:
        return None, src
    return amset_op_set(rots, taus), src + " + amset_op_set"


def amset_op_set(rots, taus):
    """按 AMSET 的**实际用法**拼操作集：R^T 与 −R^T（平移 ±τ）拼接，按旋转去重、保留先出现者。

    （AMSET symmetry.py 就是这么做的；有反演时 TR 副本在去重时被丢掉。）
    返回 [(R_real, tau_real, is_tr), ...]：R_real/tau_real 是**该操作真实的实空间** (R, τ)。
    """
    import numpy as np
    RT = np.transpose(np.asarray(rots, float), (0, 2, 1))
    rots_a = np.concatenate([RT, -RT])
    taus_a = np.concatenate([np.asarray(taus, float), -np.asarray(taus, float)])
    _half = len(rots_a) // 2
    is_tr = np.array([False] * _half + [True] * _half)
    _, first = np.unique(rots_a.reshape(len(rots_a), -1), axis=0, return_index=True)
    first = np.sort(first)
    ops = []
    for i in first:
        rot, tau, tr = rots_a[i], taus_a[i], bool(is_tr[i])
        ops.append(((-rot if tr else rot).T, (-tau if tr else tau), tr))
    return ops


def op_phase_exact(R, tau, is_tr, atol=OP_TOL):
    """AMSET 原公式在该操作上是否**精确**的充要条件（用户 2026-09-28 推导）。

    容差 `OP_TOL = 1e-4` 是实测选出来的：
      · 「好」操作经过 POSCAR 往返（pymatgen 写 8 位有效数字）后偏差 <= 4e-6；
      · 「坏」操作的最小偏差是 **0.02**（GaN 平移构型）到 0.5（Si 反演中心构型）。
    两者之间有 3 个数量级的空档，1e-4 既能吸收坐标精度噪声，又远低于任何真实违规。
    （一开始用 1e-6 会把往返后的正常结构误判成 13/24 坏 —— 测试里抓到的。）
    """
    import numpy as np
    I = np.eye(3)
    v = (np.asarray(R, float) - I) @ np.asarray(tau, float) if is_tr \
        else (np.asarray(R, float) + I) @ np.asarray(tau, float)
    return bool(np.allclose(v - np.rint(v), 0.0, atol=atol))


def symmetry_gate(structure, atol=TAU_TOL, symprec=AMSET_SYMPREC, desym_fix=False,
                  ncl=False):
    """去对称化路径的统一裁决（逐操作精确判据）。返回 dict：

        has_inversion   : bool|None   —— 是否含反演操作
        max_tau         : float|None  —— 非整数分数平移的最大绝对值
        tau_ops         : int|None    —— |tau|>atol 的对称操作数（供打印，不参与裁决）
        total_ops       : int|None    —— AMSET 实际操作集大小（按 AMSET 的 symprec 取）
        bad_ops         : int|None    —— 其中**原公式会算错**的操作数（与补丁开没开无关，照打）
        tr_ops          : int|None    —— 操作集里的时间反演操作数（有反演时为 0）
        strict_ops      : int|None    —— symprec=1e-4 时 spglib 的操作数（诊断：与 AMSET 不同就提示）
        needs_full_grid : bool        —— 是否需要 S3b/S4b 全网格（WAVEFUNCTION_FULL）
        needs_full_grid_original : bool —— 不开相位补丁时的裁决（对照用）
        desym_fix / ncl / symprec / op_source
        reason          : str

    裁决：
      · ncl（SOC）且操作集里有 TR 操作 -> 必须全网格（AMSET 的 ncl 分支对 TR 操作不取共轭、
        不乘 iσ_y，与相位无关、补丁也不修它 —— 从源码读出，未经数值验证）；
      · bad_ops == 0 -> IBZ 精确（任何空间群、任何原点）；
      · bad_ops > 0 且 desym_fix（amset_desym_fix 补丁生效）-> IBZ 精确（模型检验 12/12）；
      · bad_ops > 0 且未开补丁 -> 必须全网格。
    结构 / 对称性判不出来时保守按"需要全网格"。
    """
    import numpy as np
    out = {"has_inversion": None, "max_tau": None, "tau_ops": None, "total_ops": None,
           "bad_ops": None, "tr_ops": None, "strict_ops": None, "needs_full_grid": True,
           "needs_full_grid_original": True, "desym_fix": bool(desym_fix),
           "ncl": bool(ncl), "symprec": symprec, "op_source": None}
    if structure is None:
        out["reason"] = "取不到结构 -> 保守按需要全网格"
        return out
    try:
        rots_raw, taus_raw, _src = amset_symmetry_dataset(structure, symprec)
        ops, op_src = amset_symmetry_ops(structure, symprec)
    except Exception as e:                                   # noqa: BLE001
        rots_raw, ops, op_src = None, None, "%s: %s" % (type(e).__name__, e)
    if rots_raw is None or ops is None:
        out["reason"] = "取对称操作失败（%s）-> 保守按需要全网格" % op_src
        return out
    tn = taus_raw - np.rint(taus_raw)
    out["op_source"] = op_src
    out["has_inversion"] = bool(any(np.allclose(R, -np.eye(3), atol=1e-5) for R in rots_raw))
    out["max_tau"] = float(np.abs(tn).max()) if tn.size else 0.0
    out["tau_ops"] = int(np.sum(np.any(np.abs(tn) > atol, axis=1))) if len(tn) else 0
    out["total_ops"] = len(ops)
    bad = [i for i, (R, tau, tr) in enumerate(ops) if not op_phase_exact(R, tau, tr, OP_TOL)]
    out["bad_ops"] = len(bad)
    out["tr_ops"] = int(sum(1 for _R, _t, tr in ops if tr))
    _strict = _sym_ops_max_tau(structure, 1e-4)[1]
    out["strict_ops"] = None if _strict is None else int(len(_strict))
    _note = ""
    if out["strict_ops"] is not None and out["strict_ops"] != len(rots_raw):
        _note = ("；注意：AMSET 的 symprec=%g 认出 %d 个空间群操作，严格容差 1e-4 只有 %d 个"
                 "（结构只是近似对称，去对称化用的是近似操作）"
                 % (symprec, len(rots_raw), out["strict_ops"]))
    _tr_bad = sum(1 for i in bad if ops[i][2])
    out["needs_full_grid_original"] = bool(bad) or (bool(ncl) and out["tr_ops"] > 0)
    if ncl and out["tr_ops"]:
        out.update(needs_full_grid=True, reason=(
            "SOC（非共线）且操作集含 %d 个时间反演操作：AMSET 的 ncl 去对称化对 TR 操作"
            "不取共轭、不乘 iσ_y（与相位补丁无关）-> 必须全网格；原公式另有 %d/%d 个操作相位错"
            % (out["tr_ops"], len(bad), len(ops)) + _note))
    elif not bad:
        out.update(needs_full_grid=False, reason=(
            "AMSET 原公式在全部 %d 个操作上精确（逐操作条件 (R±I)τ≡0 全满足）-> IBZ 即可"
            % len(ops) + ("；SOC 分支未经实数据验证" if ncl else "") + _note))
    elif desym_fix:
        out.update(needs_full_grid=False, reason=(
            "AMSET 原公式会在 %d/%d 个操作上算错（其中 TR 分支 %d 个），但 DESYM_FIX 已开"
            "（amset_desym_fix 相位补丁）-> IBZ 精确"
            % (len(bad), len(ops), _tr_bad)
            + ("；★ SOC 分支同一相位因子、未经实数据验证" if ncl else "") + _note))
    else:
        out.update(needs_full_grid=True, reason=(
            "AMSET 原公式会在 %d/%d 个操作上算错（其中 TR 分支 %d 个）"
            "-> 必须全网格；或先把原点平移后复核，或开 DESYM_FIX（相位补丁）"
            % (len(bad), len(ops), _tr_bad) + _note))
    return out


def gate_log_line(g):
    """判据的一行摘要（gen / preflight 的日志统一用它，补丁开没开都打 bad_ops/total_ops）。"""
    if not g or g.get("bad_ops") is None:
        return "对称性判据：%s" % (g or {}).get("reason", "不可用")
    return ("对称性判据：bad_ops=%s/%s，TR 操作 %s，|tau|max=%.4f，DESYM_FIX=%s%s -> %s（%s）"
            % (g["bad_ops"], g["total_ops"], g.get("tr_ops"), g.get("max_tau") or 0.0,
               "on" if g.get("desym_fix") else "off", "，SOC" if g.get("ncl") else "",
               "需要全网格" if g["needs_full_grid"] else "IBZ 即可", g["reason"]))


def write_full_grid_marker(poscar, outdir=None, desym_fix=False, ncl=False):
    """算对称性判据并把结论落盘到 <outdir>/full_grid_needed.json；成功返回判据 dict，否则 None。

    S3（step3_uniform）用它把结论**提前**写下来，下游（S8/S8.4 的 gen、overlap_preflight）
    与事后追查都读同一份，不必各自重判。落的是**判据结论**，不是"分支开没开"：
    gen 在集群上跑，读不到项目里的 optional_steps.wavefunction_full（那个键不上集群）。
    needs_full_grid_original 是不开相位补丁时的结论（S8 gen 时会按当时的 DESYM_FIX 重判）。
    """
    try:
        import json
        from pymatgen.core import Structure
        p = Path(poscar)
        if not p.is_file():
            return None
        g = symmetry_gate(Structure.from_file(str(p)), desym_fix=desym_fix, ncl=ncl)
        out = Path(outdir) if outdir else p.parent
        (out / "full_grid_needed.json").write_text(
            json.dumps({"needs_full_grid": bool(g["needs_full_grid"]),
                        "needs_full_grid_original": bool(g["needs_full_grid_original"]),
                        "bad_ops": g["bad_ops"],
                        "total_ops": g["total_ops"],
                        "tr_ops": g["tr_ops"],
                        "has_inversion": g["has_inversion"],
                        "max_tau": g["max_tau"],
                        "tau_ops": g["tau_ops"],
                        "strict_ops": g["strict_ops"],
                        "desym_fix": g["desym_fix"],
                        "ncl": g["ncl"],
                        "symprec": g["symprec"],
                        "op_source": g["op_source"],
                        "reason": g["reason"],
                        "poscar": str(p)},
                       ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8")
        return g
    except Exception:
        return None


def _slab_center_span(fz, h_perp=None, vac_min=5.0):
    """沿 c 的分数坐标 -> (层心, 跨度, 是否 slab)。按最大空隙切开，跨周期边界也对。

    是否 slab：最大空隙 × h⊥ >= vac_min Å（h⊥ 未给时只看分数 >= 0.25）。
    """
    import numpy as np
    s = np.sort(np.mod(np.asarray(fz, float), 1.0))
    if s.size == 0:
        return 0.5, 0.0, False
    gaps = np.diff(np.concatenate([s, [s[0] + 1.0]]))
    i = int(np.argmax(gaps))
    lo = s[(i + 1) % s.size]                 # 最大空隙上沿 = 层底
    span = 1.0 - float(gaps[i])
    is_slab = (float(gaps[i]) * h_perp >= vac_min) if h_perp else (float(gaps[i]) >= 0.25)
    return float((lo + span / 2.0) % 1.0), span, bool(is_slab)


def align_origin(poscar, atol=1e-4):
    """把晶体原点平移到"让全部对称操作的分数平移 tau 都变 0"的位置。

    ★ 2026-09-27 用户指示，依据 VERIFICATION V109：
    AMSET 的去对称化 bug 只出在**分数平移 tau** 的处理上。简单空间群若所给原点不在
    高对称位置，symmetry ops 就带上分数平移，bug 触发：
        MoS2 P-6m2  |tau|max=0.3036 -> desym 逐带 cos 中位 0.128 / 0.16（坏）
        原点平移到 Mo  |tau|max=0.0000 -> desym 逐带 cos 中位 1.0000（精确）
    于是**只有真正的非简单空间群**（如 Pnma 的 SnSe）才必须走全网格。

    做法：spglib 取 (R, tau)，解最小二乘 (I-R)s = ±tau 得候选平移量，再**逐候选复核**
    （平移后重跑 spglib，要求 |tau|max < atol）；只有复核通过才改写 POSCAR。
    任何异常或复核失败都**原样不动**（安全侧）。

    ★ [patch_align_origin_slab] 2D slab 额外约束：层必须**连续且居中**。
    有 σh（z->-z）时 s_z 只定到 mod 1/2：最小范数解可能把镜面放到 z=0，层被
    to_unit_cell 劈成两半（上半在 z≈0、下半在 z≈0.9）。VASP 不在乎，但 S8/S8.4 的
    层厚（vdw_thickness / _read_poscar_cz 用 max(z)-min(z)）会读成 ≈c -> 2D 归一化
    因子 h⊥/t 错数倍，或 t>=c 直接 exit。所以对 slab：每个候选再试 s_z+1/2 和
    "层心移到 0.5"（点群没有翻 z 的操作时 s_z 任意），复核 tau=0 后选
    (不跨界, |层心-0.5| 最小) 的那个。3D 行为不变。

    ⚠ 物理量（能量/能带/形变势）不受影响；变的是平面波系数的 G 相位，因此
    **整条链必须同口径** —— 放在 S3 改 POSCAR，S4/S8 都从 S3 接力，天然一致，
    且 S8 的 STRUCT_CANDS 第一个就是 step3_uniform。

    Returns:
        (max_tau_before, max_tau_after, shift) —— 未平移返回 None。
    """
    try:
        import numpy as np
        from pymatgen.core import Structure
        from pymatgen.io.vasp.inputs import Poscar
        if not Path(poscar).is_file():
            return None
        st = Structure.from_file(str(poscar))
        before, rots, taus = _sym_ops_max_tau(st, atol)
        if before is None or before < atol or rots is None:
            return None                      # 不需要平移（或判不出来 -> 不动）
        # 候选：最小二乘解 (I-R)s = ±tau，以及"最重原子搬到原点"
        cands = []
        A = np.vstack([np.eye(3) - R for R in rots])
        for sgn in (+1.0, -1.0):
            b = np.concatenate([sgn * t for t in taus])
            try:
                s, *_ = np.linalg.lstsq(A, b, rcond=None)
                cands.append(-s)
            except Exception:
                pass
        _zs = [x.Z for x in st.species]
        cands.append(np.asarray(st.frac_coords[int(np.argmax(_zs))], float))

        # 2D slab：真空沿 c（本流程强制）。给每个候选补 z 方向的等价/自由平移。
        _M = st.lattice.matrix
        _h = abs(np.linalg.det(_M)) / max(np.linalg.norm(np.cross(_M[0], _M[1])), 1e-12)
        zc0, span0, is_slab = _slab_center_span(st.frac_coords[:, 2], _h)
        if is_slab:
            ext = []
            for c0 in cands:
                for dz in (0.0, 0.5):
                    c1 = np.array(c0, float); c1[2] += dz; ext.append(c1)
                c2 = np.array(c0, float); c2[2] = zc0 - 0.5; ext.append(c2)
            cands = ext

        best = None
        for cand in cands:
            st2 = st.copy()
            st2.translate_sites(range(len(st2)), -cand, frac_coords=True, to_unit_cell=True)
            after, _, _ = _sym_ops_max_tau(st2, atol)
            if after is None or after >= atol:
                continue
            if not is_slab:                  # 3D：第一个复核通过的就用（行为同旧版）
                best = (0, 0.0, cand, st2, after)
                break
            fz = np.mod(st2.frac_coords[:, 2], 1.0)
            wraps = int((fz.max() - fz.min()) > span0 + 1e-6)   # 层被周期边界劈开
            zc, _sp, _ = _slab_center_span(fz, _h)
            score = (wraps, abs(zc - 0.5))
            if best is None or score < best[:2]:
                best = (wraps, abs(zc - 0.5), cand, st2, after)
        if best is None:
            return (before, None, None)      # 复核失败
        if is_slab and best[0]:
            return (before, None, None)      # 所有 tau=0 的原点都会劈开层 -> 不动（安全侧）
        Poscar(best[3], sort_structure=False).write_file(str(poscar))
        return (before, best[4], [round(float(x), 6) for x in best[2]])
    except Exception:
        return None


def find_prev_dir(cwd: Path, candidates):
    """按顺序找第一个存在且有 CONTCAR 的目录名。"""
    for name in candidates:
        d = cwd / name
        if (d / "CONTCAR").is_file():
            return d
    return None


def new_jobname(cwd: Path, step_label: str):
    return "%s-ke-dft-cpu-%s" % (cwd.name, step_label)


def amset_env_name(cwd=None, fallback="amset_clean"):
    """AMSET 运行环境名：读 step.conf 的 [params] AMSET_ENV。

    AMSET_ENV 由 autozt 从 setting/<集群>.yaml 的 amset_env 注入 step.conf
    （autozt/report.py，标签 [cluster:<hpc>]），项目级 step.conf 可覆盖。
    gen 的 cwd 是材料目录、step.conf 就在那儿；缺失或解析失败时回退
    fallback（本地直跑 / 旧项目）。用于渲染 submit_amset.tpl 的 {{AMSET_ENV}}。
    """
    try:
        import stepconf as _sc
        base = Path(cwd) if cwd else Path(".")
        txt = (base / _sc.CONF_NAME).read_text(encoding="utf-8-sig")
        for k, v, _t in _sc.parse(txt, _sc.CONF_NAME).get("params", []):
            if k.upper() == "AMSET_ENV" and v:
                return v
    except Exception:
        pass
    # 2026-09-22（全局切 0.5.1）：绝不能再"静默"回退——旧默认 amset_clean 是 0.4.19，
    # 静默回退会让形变势又被减半且毫无提示。这里强制打印醒目警告，并回报实际用的环境名。
    import sys as _sys
    print("[WARN] step.conf 里没读到 AMSET_ENV，回退到 %r —— 若这不是本集群的 AMSET 环境，"
          "本步结果会口径错误（0.4.19 的 ISPIN=2 形变势会减半）。请在 step.conf 写 "
          "AMSET_ENV=<集群环境名>。" % fallback, file=_sys.stderr)
    return fallback


def patch_submit_jobname(submit: Path, jobname: str):
    text = submit.read_text(encoding="utf-8")
    text = text.replace("{{JOBNAME}}", jobname)
    submit.write_text(text, encoding="utf-8", newline="\n")


# --------------------------------------------------------------------------
# patch_ke_dag：泛函继承（从 step1 的 workflow_method.txt 读 FUNC=）
# --------------------------------------------------------------------------
# 为什么要继承：几何是 step1 用某个泛函优化出来的，下游单点用同一泛函才自洽。
#   pbesol -> pbesol      同一套，论文里一句话说得清
#   pbe    -> pbe
#   pbe-d3 -> pbe-d3      D3 只加在总能/力/应力上，不进 KS 哈密顿量，
#                         所以对 uniform / deform 的本征值是恒等操作（保留只为
#                         INCAR 一致）；对 elastic-dft-cpu 是必须保留的（应力有贡献）。
# 例外：DFPT（IBRION=8）。VASP 手册明确写了 vdW 修正不进 DFPT 声子响应，
#   写上 IVDW 只会污染总能而不改 ε₀ 的离子部分，所以 step5 默认剥掉，
#   由 gen 脚本的 KEEP_D3_IN_DFPT 控制。
FUNC_MAP = {
    "pbe": {"GGA": "PE", "IVDW": None,
            "VDW_LINE": "# IVDW disabled: plain PBE"},
    "pbesol": {"GGA": "PS", "IVDW": None,
               "VDW_LINE": "# IVDW disabled: PBEsol"},
    "pbe-d3": {"GGA": "PE", "IVDW": "12",
               "VDW_LINE": "IVDW   = 12            # PBE + DFT-D3(BJ)"},
}
SUPPORTED_FUNCS = tuple(FUNC_MAP)


def read_method_func(method_file: Path):
    """从 workflow_method.txt 读 FUNC=；读不到或不认识返回 None。"""
    if not Path(method_file).is_file():
        return None
    for ln in Path(method_file).read_text(errors="ignore").splitlines():
        if ln.strip().upper().startswith("FUNC="):
            v = ln.split("=", 1)[1].strip().lower()
            return v if v in FUNC_MAP else None
    return None


def sniff_func_from_incar(incar: Path):
    """兜底：从 step1 的 INCAR 反推 GGA/IVDW。推不出返回 None。"""
    if not Path(incar).is_file():
        return None
    gga, ivdw = "", None
    for ln in Path(incar).read_text(errors="ignore").splitlines():
        ln = ln.split("#", 1)[0].split("!", 1)[0].strip()
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        k, v = k.strip().upper(), v.strip()
        if k == "GGA":
            gga = v.upper()
        elif k == "IVDW":
            ivdw = v.split()[0] if v.split() else None
    if gga == "PS" and ivdw is None:
        return "pbesol"
    if gga == "PE" and ivdw == "12":
        return "pbe-d3"
    if gga == "PE" and ivdw is None:
        return "pbe"
    return None


def resolve_func(prev_dir: Path, setting, step_name, drop_d3=False,
                 fallback="pbe"):
    """定出本步用的泛函，返回 (func, subs)。

    setting = "inherit" 时从 prev_dir/workflow_method.txt 读，读不到就嗅探
    prev_dir/INCAR，再读不到用 fallback 并告警。
    setting 直接写死泛函名时原样采用。
    drop_d3=True 会把 pbe-d3 降级成 pbe（DFPT 专用）。
    subs 是给 render_tpl 的占位符字典：{"GGA": ..., "VDW_LINE": ...}
    """
    s = str(setting).lower()
    if s == "inherit":
        func = read_method_func(Path(prev_dir) / METHOD_FILE)
        src = "workflow_method.txt"
        if func is None:
            func = sniff_func_from_incar(Path(prev_dir) / "INCAR")
            src = "嗅探 step1/INCAR"
        if func is None:
            func, src = fallback, "都读不到，回退默认值"
            print("[WARN] %s：无法从 %s 判定泛函，回退 FUNC=%s。"
                  "若 step1 是 pbesol/pbe-d3，结果会不自洽！"
                  % (step_name, prev_dir, func))
    elif s in FUNC_MAP:
        func, src = s, "脚本内写死"
    else:
        sys.exit("[ERROR] %s：FUNC=%r 无效，只允许 inherit / %s"
                 % (step_name, setting, " / ".join(FUNC_MAP)))

    eff = func
    if drop_d3 and func == "pbe-d3":
        eff = "pbe"
        print("[..] %s：DFPT 不支持 vdW 修正（VASP 手册），"
              "pbe-d3 -> pbe（几何仍是 D3 优化的）" % step_name)
    m = FUNC_MAP[eff]
    print("[..] %s：泛函 %s（来源：%s）-> GGA=%s IVDW=%s"
          % (step_name, eff, src, m["GGA"], m["IVDW"] or "off"))
    return eff, {"GGA": m["GGA"], "VDW_LINE": m["VDW_LINE"]}
