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
# [V173] S3（AMSET 能带 / 波函数 / 形变势网格的唯一来源）网格是否满足现行规则。
#   规则在 step3_uniform 的 gen 里（DK_MAX_2D 0.05 / DK_MAX_3D 0.06 Å⁻¹，2D 真空轴 kz ≥ 3），可只在 gen 运行时生效：
#   旧规则下算完的 S3 一直判完成、不会重新生成，S7 照抄它的网格，S8/S8.4 照用 —— 谁都不查。
#   CrS₂_hex：S3 还是 08-26 的 15×15×1（面内间距 0.16 Å⁻¹ = 规则的 3.2 倍、kz 只有一层），CrSe₂ 是 48×48×3；
#   两者 S8.4 的 ADP/DPT 一个 1.7、一个 0.5。kz 只有一层时 AMSET 沿 kz 外推（V8：插值网格一变 ADP 差 30–45%）。
# [V174] S8.2（DPT 的 m* 直接拟合 S3 网格点）、S8.1（BoltzTraP2 插值 S3）也读这张网格，V173 没管：
#   CrS₂_hex 重跑 S3（15×15×1 -> 48×48×3）后 m* 0.967/1.009 -> 0.866/0.883，DPT +25%/+31%：
#   粗网格的拟合窗口吃进了非抛物区，m* 拟重。面内 m* 与真空轴 kz 无关，这两步 check_kz=False。
# [V175] 还有两处漏洞：
#   ① S7 只在自己 gen 时照抄 S3 网格。S3 重算后 S7 没跟着重跑，S3 闸门照样通过（S3 是新的）、lineage 也通过
#      （S7.1 比 S7 新），S8/S8.4 就拿旧粗网格上的形变势去算 —— s7_grid_issues 查 S7 和 S3 是否同一张网格；
#   ② 闸门只在重新生成时起作用，已经算完的旧结果没人再查，CrS₂ 当初就是这样"收口"的 ——
#      grid_lineage 给 tools/lineage_check.py 用。
# --------------------------------------------------------------------------
S3_DK_MAX = {"2d": 0.05, "3d": 0.06}
S3_KZ_MIN_2D = 3
S3_DK_ERR_FACTOR = 2.0          # 面内间距超过规则这么多倍 -> 报错（规则以内的小差别只告警）


def _mesh_rule_issues(mesh, s3, dim, check_kz, where):
    """一张网格（按 S3 的晶胞和维度）是否合现行规则 -> [(level, msg)]；读不到晶胞返回 []。"""
    import numpy as np
    pos = s3 / "POSCAR"
    if not pos.is_file():
        return []
    try:
        dim = dim or read_method_dim(s3 / METHOD_FILE) or resolve_dim_for(pos, "auto")[0]
        vac = resolve_dim_for(pos, dim)[1] if dim == "2d" else None
        rec = 2.0 * np.pi * np.linalg.inv(read_lattice_matrix(pos)).T
    except (OSError, ValueError, IndexError, SystemExit):
        return []
    dk = S3_DK_MAX.get(dim, S3_DK_MAX["3d"])
    out = []
    for i in range(3):
        if i == vac:
            continue
        sp = float(np.linalg.norm(rec[i])) / max(int(mesh[i]), 1)
        if sp > S3_DK_ERR_FACTOR * dk:
            out.append(("error", "%s 第 %d 轴 %d 分：笛卡尔间距 %.3f Å⁻¹，是现行规则 %.2f 的 %.1f 倍"
                        % (where, i + 1, mesh[i], sp, dk, sp / dk)))
        elif sp > 1.1 * dk:
            out.append(("warn", "%s 第 %d 轴 %d 分：笛卡尔间距 %.3f Å⁻¹，比现行规则 %.2f 粗"
                        % (where, i + 1, mesh[i], sp, dk)))
    if check_kz and vac is not None and int(mesh[vac]) < S3_KZ_MIN_2D:
        out.append(("error", "%s 真空轴 kz = %d < %d：AMSET 沿 kz 外推而不是内插"
                    "（V8：插值网格一变 ADP 差 30–45%%）" % (where, mesh[vac], S3_KZ_MIN_2D)))
    return out


def s3_grid_issues(mat_dir, dim=None, check_kz=True):
    """材料目录 -> [(level, msg)]，level 为 "error" / "warn"；读不到 S3 的 KPOINTS / POSCAR 返回 []。
    check_kz=False：只查面内（S8.1 / S8.2 只用面内能带，2D 真空轴 kz 一层不影响它们）。"""
    s3 = Path(mat_dir) / "step3_uniform"
    mesh = read_kpoints_mesh(s3 / "KPOINTS")
    if not mesh:
        return []
    return _mesh_rule_issues(mesh, s3, dim, check_kz, "step3_uniform")


def s7_grid_issues(mat_dir, dim=None):
    """[V175] S7 形变子目录（含 ionrelax/）的网格和 S3 不同 -> [(level, msg)]，每种网格一条：
    S7 这张网格本身不合现行规则 -> error（S3 重算过、S7 没跟着重跑，形变势还在旧粗网格上）；合规只是不同 -> warn。
    没有 S3 KPOINTS / S7 目录返回 []。"""
    mat = Path(mat_dir)
    s3 = mat / "step3_uniform"
    ref = read_kpoints_mesh(s3 / "KPOINTS")
    s7 = mat / "step7_deform"
    if not ref or not s7.is_dir():
        return []
    ref = [int(x) for x in ref]
    diff = {}
    for kp in sorted(set(s7.glob("*/KPOINTS")) | set(s7.glob("*/ionrelax/KPOINTS"))):
        m = read_kpoints_mesh(kp)
        if m and [int(x) for x in m] != ref:
            diff.setdefault(tuple(int(x) for x in m), []).append(str(kp.parent.relative_to(s7)))
    out = []
    for m, dirs in sorted(diff.items()):
        # [V176] 写事实（哪条规则不合、比 S3 粗还是细），原因只在 S7 确实比 S3 粗时才推断：
        #   CrS2_ortho 的 S7 是 48×48×1、S3 是 15×15×1 —— S7 面内反而更细，只是 kz = 1；V175 的报错一律写成
        #   "S3 重算过、S7 没跟着重跑"，agent 照字面理解成了反方向。
        errs = [msg[len("step7_deform "):].split("：AMSET")[0]
                for lv, msg in _mesh_rule_issues(m, s3, dim, True, "step7_deform") if lv == "error"]
        where = ", ".join(dirs[:4]) + (" 等 %d 个" % len(dirs) if len(dirs) > 4 else "")
        head = "step7_deform 的 %s 是 %s 网格，S3 是 %s" % (where, "×".join(map(str, m)), "×".join(map(str, ref)))
        if not errs:
            out.append(("warn", head + "（两张都合现行规则，AMSET 会插值 D(k)；要完全一致就 retry S7 -> S7.1）"))
            continue
        hint = ("S7 有轴比 S3 粗：多半是 S3 重算过、S7 没跟着重跑" if any(m[i] < ref[i] for i in range(3))
                else "S7 每个轴都不比 S3 粗：是 S7 按旧规则自己定的网格（09-24 之前 S7 不照抄 S3）")
        out.append(("error", "%s；S7 这张网格不合现行规则（%s），形变势 D(k) 是在它上面算的 —— %s"
                    % (head, "；".join(errs), hint)))
    return out


def grid_lineage(mat_dir):
    """[V175] tools/lineage_check.py 用 -> [(步骤, 说明)]：只收 error（略粗于规则的告警不算）。
    gen 的闸门只在重新生成时起作用，已经算完的旧结果只能靠这里查出来。"""
    inplane = {m for lv, m in s3_grid_issues(mat_dir, check_kz=False) if lv == "error"}
    probs = [("step3_uniform", m + ("（S4/S7 照抄这张网格，S8/S8.4/S8.1/S8.2 的结果都在它上面）" if m in inplane
                                    else "（S8/S8.4 的 AMSET 结果受影响）"))
             for lv, m in s3_grid_issues(mat_dir) if lv == "error"]
    probs += [("step7_deform", m + "（S7.1 / S8 / S8.4 受影响）") for lv, m in s7_grid_issues(mat_dir) if lv == "error"]
    return probs


def s3_grid_gate(mat_dir, label, allow=False, dim=None, check_kz=True, uses=None, redo=None, with_s7=False):
    """S8 / S8.4 / S8.2 gen 调：S3 网格不满足现行规则 -> 打印；有 error 且 allow=False -> sys.exit。返回 issues。
    uses：本步拿这张网格做什么（报错用）；redo：重跑顺序。None = S8/S8.4（AMSET）的说法。
    with_s7（S8/S8.4）：S7 的网格也要和 S3 一致（s7_grid_issues），不合规同样拦。"""
    issues = s3_grid_issues(mat_dir, dim, check_kz=check_kz)
    s7 = s7_grid_issues(mat_dir, dim) if with_s7 else []
    for lv, msg in issues + s7:
        print("[%s] %s：%s" % ("ERROR" if lv == "error" and not allow else "WARN", label, msg))
    if any(lv == "error" for lv, _ in s7) and not any(lv == "error" for lv, _ in issues) and not allow:
        sys.exit("[ERROR] %s：S3 已是现行网格，但 step7_deform 的网格不合现行规则（见上面几行），"
                 "形变势 D(k) 是在那张网格上算的，结果不可信。\n"
                 "        处理：retry step7_deform（会照抄 S3 的网格）-> S7.1 -> 本步。\n"
                 "        确实要用这套网格（复现旧结果）：本步 step.conf 写 ALLOW_COARSE_S3 = true。" % label)
    if any(lv == "error" for lv, _ in issues + s7):
        if not allow:
            sys.exit("[ERROR] %s：step3_uniform 的网格是旧规则下生成的，%s，结果不可信。\n"
                     "        处理：retry step3_uniform（现行规则：2D 面内 ≤ %.2f Å⁻¹、kz ≥ %d）-> %s。\n"
                     "        确实要用这张网格（复现旧结果）：本步 step.conf 写 ALLOW_COARSE_S3 = true。"
                     % (label, uses or "AMSET 的能带 / 波函数 / 形变势都在这张网格上", S3_DK_MAX["2d"],
                        S3_KZ_MIN_2D, redo or "S4 -> S7 -> S7.1 -> 本步；S7 会照抄新的 S3 网格"))
        print("[WARN] %s：ALLOW_COARSE_S3 = true —— 按旧网格继续，结果只作复现/对照" % label)
    return issues + s7


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
    # [patch_rerun_invalidate V140] autozt rerun 先 rm -rf 本步目录再 gen：快照全空，上面比不出变化，
    #   可下游（S4 的 h5、S8 …）是用被删掉的那一版算的 -> 同样失效。第一次 gen 时下游不存在，什么都不做。
    fresh = bool(snap) and all(v is None for v in snap.values())
    n = 0
    if changed:
        n = archive_stale_outputs(outdir, "stale-input-" + _t.strftime("%Y%m%d%H%M%S"))
        if n:
            print("[..] patch_stale_input：%s 变了（%s）-> 归档旧产物 %d 个，本步会重算"
                  % (step, "/".join(changed), n))
    down = []
    if n or n_grid_archived:
        down = invalidate_downstream(cwd, step, "输入变了：%s" % ("/".join(changed) or "网格"))
    elif fresh:
        down = invalidate_downstream(cwd, step, "%s 目录是新建的（rerun）" % step)
    return changed, n, down


# --------------------------------------------------------------------------
# [V168] 产物是不是用【现在的】INCAR 跑出来的：拿 vasprun.xml 的 <incar>（VASP 读到的 INCAR 原文键）
#   和 <parameters>（含默认值的生效参数）跟现在的 INCAR 比。
#   snapshot_inputs 只能在 gen 覆盖输入之前拍快照；输入已经被上一次 gen 换掉、旧产物还在时就比不出来了
#   （CrS₂ S7：retry 换上带 LVHAR 的新 INCAR，undeformed/deform-05..08 的旧 OUTCAR 照样判完成、
#   不进 fan_todo，永远不重算，S7.1 读不到 LOCPOT）——这时只有运行时的回显能说明问题。
#   比较按语义：并行/标签键（_INCAR_NONPHYS）不算；n*x 展开；T/.TRUE. 等价；数值按相对 1e-6
#   且容忍 vasprun 定点 8 位小数的舍入；字符串大小写不敏感、前缀相同即算一致（VASP 只认前几个字母，
#   <parameters> 里的 PREC 写作 "accura"）。新 INCAR 里有、VASP 两处都没回显的键（VASP 不认识）跳过。
# --------------------------------------------------------------------------
_VASPRUN_HEAD_MAX = 64 * 1024 * 1024
# 只决定"写不写某个输出文件"、不影响 SCF 结果的键：只有新 INCAR 要了旧运行没给的输出才算变化
#   （补 LVHAR 要 LOCPOT -> 要重算；关掉 LWAVE -> 旧产物照用）。
_INCAR_OUTPUT_ONLY = {"LWAVE", "LCHARG", "LVHAR", "LVTOT", "LELF", "LORBIT", "LAECHG"}
_XML_ITEM = re.compile(r"<(i|v)\b([^>]*)>(.*?)</\1>", re.S)
_XML_NAME = re.compile(r'name="([^"]+)"')


def _vasprun_incar(vasprun):
    """vasprun.xml -> (<incar> 字典, <parameters> 字典)；没有 <incar> 返回 (None, None)。只读到 </parameters>。"""
    buf = ""
    try:
        with open(vasprun, encoding="utf-8", errors="ignore") as f:
            while len(buf) < _VASPRUN_HEAD_MAX:
                chunk = f.read(1 << 20)
                if not chunk:
                    break
                buf += chunk
                if "</parameters>" in buf:
                    break
    except OSError:
        return None, None

    def block(tag):
        i = buf.find("<%s>" % tag)
        j = buf.find("</%s>" % tag, i + 1) if i >= 0 else -1
        return buf[i:j] if j > i else None

    def items(txt):
        d = {}
        for m in _XML_ITEM.finditer(txt or ""):
            n = _XML_NAME.search(m.group(2))
            if n:
                d.setdefault(n.group(1).upper(), m.group(3).strip())
        return d

    inc = block("incar")
    if inc is None:
        return None, None
    return items(inc), items(block("parameters"))


def _incar_tokens(v):
    out = []
    for t in str(v).replace(",", " ").split():
        n, star, x = t.partition("*")
        if star and n.isdigit() and x:
            out.extend([x] * int(n))
        else:
            out.append(t)
    return out


def _incar_atom(t):
    s = t.strip().strip("'\"").lower()
    if s in (".true.", "true", "t", ".t."):
        return True
    if s in (".false.", "false", "f", ".f."):
        return False
    try:
        return float(s.replace("d", "e"))
    except ValueError:
        return s


def _incar_same(a, b):
    ta, tb = _incar_tokens(a), _incar_tokens(b)
    if len(ta) != len(tb):
        return False
    for x, y in zip(ta, tb):
        x, y = _incar_atom(x), _incar_atom(y)
        if isinstance(x, bool) or isinstance(y, bool):
            if x is not y:
                return False
        elif isinstance(x, float) and isinstance(y, float):
            if abs(x - y) > max(1e-6 * max(abs(x), abs(y)), 5e-9):
                return False
        elif isinstance(x, str) and isinstance(y, str):
            if not (x.startswith(y) or y.startswith(x)):
                return False
        else:
            return False
    return True


def run_incar_mismatch(incar_text, vasprun):
    """现在的 INCAR 文本 vs 那次运行的 vasprun.xml 回显 -> 变了的键列表 ["LVHAR: F -> .TRUE.", ...]。

    vasprun.xml 读不到 / 没有 <incar>（没跑、跑到一半就挂）-> None（判不了，调用方按原逻辑走）。"""
    old, eff = _vasprun_incar(vasprun)
    if old is None:
        return None
    new = {k: v for k, v in parse_incar(incar_text).items() if k not in _INCAR_NONPHYS}
    diff = []
    for k, v in sorted(new.items()):
        ov = old.get(k, (eff or {}).get(k))
        if ov is None:
            continue
        if _incar_same(v, ov):
            continue
        if k in _INCAR_OUTPUT_ONLY and all(_incar_atom(x) in (False, 0.0) for x in _incar_tokens(v)):
            continue                                       # 新 INCAR 不要这个输出了：旧产物照用
        diff.append("%s: %s -> %s" % (k, ov, v))
    for k in sorted(set(old) - set(new) - _INCAR_NONPHYS - _INCAR_OUTPUT_ONLY):
        diff.append("%s: %s -> (删除)" % (k, old[k]))
    return diff


def archive_if_run_mismatch(d, incar_text=None, vasprun_name="vasprun.xml", label=None):
    """目录 d 里已有的产物不是现在这份 INCAR 跑的 -> 归档 *.stale-input-<时间>。

    返回变了的键（一致 -> []；没有 vasprun.xml 回显、判不了 -> None）。"""
    import time as _t
    d = Path(d)
    if incar_text is None:
        p = d / "INCAR"
        if not p.is_file():
            return None
        incar_text = p.read_text(errors="ignore")
    diff = run_incar_mismatch(incar_text, d / vasprun_name)
    if not diff:
        return diff
    n = archive_stale_outputs(d, "stale-input-" + _t.strftime("%Y%m%d%H%M%S"))
    print("[..] %s：已有产物是旧 INCAR 跑的（%s%s）-> 归档 %d 个，本目录会重算"
          % (label or d.name, "; ".join(diff[:4]), "…" if len(diff) > 4 else "", n))
    return diff


# [V169] AMSET 形变势的"芯能级参考"（deformation.h5）有两种来源（amset/deformation/io.py get_reference_energy）：
#   首选 OUTCAR 的各原子核处平均静电势块（"the norm of the test charge is"）；没有这块（ICORELEVEL=1 时 VASP 不写）
#   就**静默**改用 1s 芯能级本征值（元素依赖、非刚性，模板注释：会把 E1 压到 ~0.3 eV 量级）。
#   两种在构型之间混用时，参考差是两种量之差，形变势没有物理意义；全是 1s 时也不是 AMSET 的标准口径。
def outcar_core_ref_kind(outcar):
    """OUTCAR -> "avg_core"（平均静电芯势）/ "1s"（只有芯能级本征值）/ None（两者都没有或读不了）。"""
    import gzip
    p = Path(outcar)
    if not p.is_file() and Path(str(p) + ".gz").is_file():
        p = Path(str(p) + ".gz")
    has_1s = False
    try:
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if "the norm of the test charge is" in ln:
                    return "avg_core"
                if "the core state eigen" in ln:
                    has_1s = True
    except OSError:
        return None
    return "1s" if has_1s else None


# 下游关系（与 skill.yaml 的 needs 对应）。mode="link"：下游靠软链引用上游产物
# （S4 链 S3 的 WAVECAR/vasprun；S8/S8.4 链 S4/S4b 的 h5 与 S3/S3b 的 vasprun）——
# 只有软链确实指向这个上游时才失效；mode="always"：下游直接读上游目录，一律失效。
DOWNSTREAM = {
    # [patch_rerun_cascade V153] VASP 步骤（S1 / S2 链）重算 -> 下游 VASP 步骤**整目录**归档（mode="dir"）。
    #   autozt 的 rerun 只删本步；下游 VASP 步骤按 OUTCAR 判"完成"，归档某个文件没用（原地重新生成又会
    #   读到旧 WAVECAR/CHGCAR），所以整目录改名 <dir>.stale-upstream-<步骤>-<时间>（数据保留），autozt 看到
    #   目录不在 -> 等上游跑完后重新生成。Mo2S3 / Si_diamond / P1_Mo-MoS2 三次都是 S2.2/S2.3 带着旧结构
    #   判"完成"，要手工 rerun。只在本步目录是新建的（rerun）时触发，见 cascade_on_rerun；下游有作业在
    #   排队/运行时不动它（只告警）。S3b/S3c 的 needs 是 S3，但结构来自 S1，直接挂在 S1 下。
    "step1_opt": [("step2_bandgap/step2.1_static", "dir"), ("step3_uniform", "dir"),
                  ("step3b_uniform_full", "dir"), ("step3c_uniform_offgrid", "dir"),
                  ("step5_dielect", "dir"), ("step6_elastic", "dir"), ("step7_deform", "dir"),
                  ("step1_std_opt", "dir")],                    # 最后一个：zt 的 kl 分支（标准胞重新弛豫）
    "step2_bandgap/step2.1_static": [("step2_bandgap/step2.15_discriminant", "dir"),
                                     ("step2_bandgap/step2.2_pbe", "dir")],
    "step2_bandgap/step2.15_discriminant": [("step2_bandgap/step2.155_discriminant_decide", "dir")],
    "step2_bandgap/step2.2_pbe": [("step2_bandgap/step2.2_pbe_plot", "dir"), ("step2_bandgap/step2.3_hse", "dir")],
    "step2_bandgap/step2.3_hse": [("step2_bandgap/step2.3_hse_plot", "dir")],
    # zt 的 kl 分支：这些步骤的 gen 在 kl-dft-cpu、不调本函数，只靠 S1 的级联递归走到。
    "step1_std_opt": [("step2_static", "dir")],
    "step2_static": [("step3_nac", "dir")],
    "step3_uniform": [("step4_wave", "link"), ("step8_amset", "link"),
                      ("step8.4_amset2d", "link"), ("step8.1_boltztrap", "always"),
                      ("step8.2_dpt", "always")],
    "step3b_uniform_full": [("step4b_wave_full", "link"), ("step8_amset", "link"),
                            ("step8.4_amset2d", "link")],
    "step4_wave": [("step8_amset", "link"), ("step8.4_amset2d", "link")],
    "step4b_wave_full": [("step8_amset", "link"), ("step8.4_amset2d", "link")],
    # [patch_stale_dp V122] S7.1 重新生成 = 形变势 h5 与 band_edges.json 变了：软链它的 S8/S8.4 与读
    #   band_edges.json 的 S8.2（DPT）一并失效。MoS2 实测：同一批形变单点、新版 S7.1 重读后电子 ADP
    #   迁移率 198 -> 406，而旧 transport.json 仍被判"完成"。
    "step7b_deform_read": [("step8_amset", "link"), ("step8.4_amset2d", "link"),
                           ("step8.2_dpt", "always")],
    # [patch_rerun_invalidate V140] S7（形变单点）重新生成 -> S7.1 读出的形变势作废（递归到 S8/S8.4/S8.2）。
    #   S7.1 的产物在材料目录的 step7b_deform_read/，rerun S7 删不到它，此前一直是旧的"完成"。
    "step7_deform": [("step7b_deform_read", "always")],
    # [patch_stale_upstream V136] S5（介电）/ S6（弹性）重新生成：S8/S8.4 在 gen 时把它们的数值写进
    #   settings.yaml（不是软链）-> mode="always"；S8.2（DPT）读 step6_elastic/OUTCAR 的面内弹性；
    #   S5.1 的校验结论随 S5 一起作废。此前没有这两项：S5/S6 重算后旧 transport.json 仍被判"完成"。
    "step5_dielect": [("step5_dielect_validate", "always"), ("step8_amset", "always"),
                      ("step8.4_amset2d", "always")],
    "step6_elastic": [("step8_amset", "always"), ("step8.4_amset2d", "always"),
                      ("step8.2_dpt", "always")],
    # [patch_gap_invalidate V148] S2 画图重新生成 = band_summary.json 的带隙可能变了。S8/S8.4 在 gen 时把带隙写进
    #   settings.yaml（"bandgap:" + "# bandgap_source: <画图目录>/band_summary.json"）——不是软链，以前 HSE 重算后
    #   S8 照样带着旧带隙判"完成"。mode="gap"：下游确实用的是这份 band_summary、且带隙变了 > GAP_TOL_EV 才失效
    #   （S8 一跑几小时，只重画一张图不该让它重排；BANDGAP_OVERRIDE 的下游来源不同，不受影响）。
    #   S8.1（BoltzTraP2 的剪刀差也读这份 band_summary）只要几分钟 -> "always"。
    "step2_bandgap/step2.2_pbe_plot": [("step8_amset", "gap"), ("step8.4_amset2d", "gap"),
                                       ("step8.1_boltztrap", "always")],
    "step2_bandgap/step2.3_hse_plot": [("step8_amset", "gap"), ("step8.4_amset2d", "gap"),
                                       ("step8.1_boltztrap", "always")],
    # step20_zt：zt v0.1（抄 ke 步骤进 zt 目录）的遗留条目。v0.2 起 S20 在 <材料>/zt-dft-cpu/，这里按 ke 目录
    #   找不到它；S20 的过期改由 autozt 判断（[0015] workflow._mark_stale_composites：比对当初的跨技能输入）。
    #   S8.1 的 c/t 由 S8 gen 里的 invalidate_ct_consumers 管（只在 c/t 真变了时重排）。
    "step8_amset": [("step8.3_output", "always"), ("step20_zt", "dir")],
    # [patch_post_invalidate V148] S8.3 的对比图也读 S8.4 的 transport.json；S8/S8.4/S8.1/S8.2 的 gen 现在都会调
    #   invalidate_downstream（以前只在表里、没有一个 gen 调：GaAs 的 S8.2 重跑了 4 次，S8.1 一直是旧 τ）。
    "step8.4_amset2d": [("step8.3_output", "always")],
    "step8.1_boltztrap": [("step8.3_output", "always")],
    "step8.2_dpt": [("step8.1_boltztrap", "always"), ("step8.3_output", "always")],
}
# skill.yaml 里不进 DOWNSTREAM 的 needs 边与理由（test_downstream_wiring.py：每条 needs 边要么在 DOWNSTREAM，
#   要么在这里；"*" = 该上游的全部下游）。
DOWNSTREAM_EXEMPT = {
    ("step3_uniform", "step3b_uniform_full"): "S3b 从 S3 的 CHGCAR 起步后自洽到收敛（ICHARG=1），终点不依赖 S3；结构由 LINEAGE_POSCARS 管",
    ("step3b_uniform_full", "step3c_uniform_offgrid"): "S3c 段 1 自己自洽，不读 S3b 的产物；结构由 LINEAGE_POSCARS 管",
    ("step6_elastic", "step8.1_boltztrap"): "S8.1 的 τ 来自 S8.2：S6 -> S8.2 -> S8.1 递归失效",
    ("step7b_deform_read", "step8.1_boltztrap"): "同上：S7.1 -> S8.2 -> S8.1 递归失效",
    # zt 的 kl 声子链：gen 在 kl-dft-cpu，不在本技能维护；位移数据集不来自 S1（step4_disp 没有 needs）
    ("step4_disp", "*"): "kl 技能的声子链（gen 在 kl-dft-cpu）",
    ("step5_fc", "*"): "kl 技能的声子链（gen 在 kl-dft-cpu）",
    ("step6_kappa", "*"): "kl 技能的声子链（gen 在 kl-dft-cpu）",
}
# DOWNSTREAM 里 gen 不在本技能的上游（zt 的 kl 分支）：它们自己不调失效，只被 S1 的级联递归走到。
DOWNSTREAM_FOREIGN = frozenset({"step1_std_opt", "step2_static"})
RUNNING_GRACE_S = 900       # [V153] 下游目录的 OUTCAR / slurm 输出这么多秒内还在更新 -> 当作在跑，不整目录归档
GAP_TOL_EV = 1e-3
# [V153] band_edges.json 的格式/口径版本：S7.1（gen_step9b）写、S8.2（gen_step12_dpt）核对，两边都读这一个常量。
#   以前两边各自拿脚本自己的 _SKILL_REV 比（S8.2 "2026-08-31-rscan" vs S7.1 "2026-09-16-deform-ref-vacuum"），
#   从来没对上过，每次 S8.2 都告警"step7b 可能跑的是旧副本"——一条永远在响的告警等于没有告警。
#   改 band_edges.json 的字段或口径时改这里，S8.2 会对旧文件告警。
BAND_EDGES_REV = "2026-09-16-deform-ref-vacuum"


def consumer_bandgap(step_dir):
    """S8/S8.4 的 settings.yaml -> (bandgap 或 None, bandgap_source 或 None)。"""
    gap = src = None
    try:
        for ln in (Path(step_dir) / "settings.yaml").read_text(errors="ignore").splitlines():
            s = ln.strip()
            if s.startswith("# bandgap_source:"):
                src = s.split(":", 1)[1].strip()
            elif s.startswith("bandgap:"):
                try:
                    gap = float(s.split(":", 1)[1].split("#", 1)[0])
                except ValueError:
                    pass
    except OSError:
        pass
    return gap, src


def _gap_changed(step_dir, step, gap):
    """下游用的是 step（画图目录）的 band_summary.json、且带隙与新值差 > GAP_TOL_EV -> 旧带隙值；否则 None。"""
    if gap is None:
        return None
    old, src = consumer_bandgap(step_dir)
    if old is None or not src or (Path(step).name + "/") not in src:
        return None
    return old if abs(float(gap) - old) > GAP_TOL_EV else None
DONE_MARKERS = {
    # S7.1：deformation_vac.h5 也要归档 —— 只归档 deformation.h5 的话，新 S7.1 没产出 vac 时
    #   _pick_deformation_h5 会捡到旧的 vac h5。
    "step7b_deform_read": ("deformation.h5", "deformation_vac.h5", "band_edges.json"),
    "step4_wave": ("wavefunction.h5",),
    "step4b_wave_full": ("wavefunction.h5",),
    "step5_dielect_validate": ("dielectric_check.json",),
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


def _active_job_dirs():
    """squeue 里本用户作业的工作目录（realpath 集合）；没有 squeue / 读不到 -> None（退回看文件时间）。"""
    try:
        r = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""), "-o", "%Z"],
                           capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    if r.returncode != 0:
        return None
    return {os.path.realpath(x.strip()) for x in r.stdout.splitlines() if x.strip()}


def _dir_busy(d, jobs):
    """[V153] 下游目录（含扇出子目录）有作业在排队/运行，或 OUTCAR/slurm 输出刚更新过 -> 说明；否则 None。"""
    import time as _t
    rd = os.path.realpath(str(d))
    if jobs:
        hit = sorted(j for j in jobs if j == rd or j.startswith(rd + os.sep))
        if hit:
            return "squeue 里有作业的工作目录在这里（%s）" % (os.path.relpath(hit[0], rd))
    now = _t.time()
    for pat in ("OUTCAR", "slurm-*.out", "*/OUTCAR", "*/slurm-*.out"):
        for f in Path(d).glob(pat):
            try:
                age = now - f.stat().st_mtime
            except OSError:
                continue
            if age < RUNNING_GRACE_S:
                return "%s 在 %.0f 秒前还在更新" % (os.path.relpath(str(f), rd), age)
    return None


def cascade_on_rerun(cwd, step, reason=None):
    """[V153] VASP 步骤的 gen 开头调（写任何文件之前）：本步目录里还没有 OUTCAR（autozt rerun 先 rm -rf 了，
    或第一次生成）-> invalidate_downstream（下游 VASP 步骤整目录归档、run:gen 产物归档）。
    已经有 OUTCAR（init -f、retry 原地重新生成）-> 不动下游，只提示。返回归档清单。"""
    cwd = Path(cwd)
    out = cwd / step
    had = [q for pat in ("OUTCAR", "*/OUTCAR") for q in out.glob(pat)] if out.is_dir() else []
    if had:
        if any((cwd / d).exists() for d, _m in DOWNSTREAM.get(step, ())):
            print("[..] %s 原地重新生成（已有 OUTCAR）：不自动归档下游。结构或输入真的变了就用 rerun，"
                  "或 tools/lineage_check.py <材料目录> --invalidate-from %s" % (step, step))
        return []
    return invalidate_downstream(cwd, step, reason or ("%s 目录是新建的（rerun）" % step))


def invalidate_downstream(cwd, step, reason, _seen=None, gap=None, _jobs=None):
    """上游 step 要重算：把下游完成标记改名 *.stale-upstream-<step>-<时间>，递归传递。

    返回 [(下游步骤, 改名的文件), ...]。只改名不删除；下游目录不存在/没有标记时什么都不做。
    gap：S2 画图的新带隙（mode="gap" 的下游只在带隙变了时失效；不给 = 这类下游一律不动）。
    mode="dir"（V153）：下游整个目录改名 <dir>.stale-upstream-<step>-<时间>，记作 (下游, "<目录>")；
    下游有作业在排队/运行时不动（★ 告警）。
    """
    import time as _t
    cwd = Path(cwd)
    seen = _seen if _seen is not None else set()
    jobs = _jobs if _jobs is not None else []          # [已查?, 结果]：整个递归只查一次 squeue
    done = []
    tag = "stale-upstream-%s-%s" % (step.replace("/", "_"), _t.strftime("%Y%m%d%H%M%S"))
    for down, mode in DOWNSTREAM.get(step, ()):
        if down in seen:
            continue
        d = cwd / down
        if not d.is_dir():
            continue
        if mode == "dir":
            seen.add(down)
            if not jobs:
                jobs.extend([True, _active_job_dirs()])
            busy = _dir_busy(d, jobs[1])
            if busy:
                print("[WARN] ★ 上游 %s 重算，但下游 %s 看起来还在跑（%s）—— 没有归档。先 scancel 它，"
                      "再 tools/lineage_check.py <材料目录> --invalidate-from %s" % (step, down, busy, step))
            else:
                d.rename(d.with_name(d.name + "." + tag))
                done.append((down, "<目录>"))
                print("[WARN] 上游 %s 重算（%s）-> 下游 %s 整个目录归档为 %s.%s，上游跑完后会重新生成"
                      % (step, reason, down, d.name, tag))
            done += invalidate_downstream(cwd, down, "上游 %s 失效" % step, seen, _jobs=jobs)
            continue
        if mode == "link" and not _links_into(d, cwd / step):
            continue
        old_gap = None
        if mode == "gap":
            old_gap = _gap_changed(d, step, gap)
            if old_gap is None:
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
                  % (step, reason if old_gap is None else "带隙 %.4f -> %.4f eV" % (old_gap, gap), down, tag))
        elif old_gap is not None:
            print("[WARN] ★ %s 的 settings.yaml 用的是旧带隙 %.4f eV（%s 现在是 %.4f eV），它还没跑完 —— "
                  "取消那个作业、rerun %s" % (down, old_gap, step, gap, down))
        done += invalidate_downstream(cwd, down, "上游 %s 失效" % step, seen, _jobs=jobs)
    return done


CT_CONSUMERS = ("step8.1_boltztrap",)
# [0016] S8.1 取 c/t 的优先级（与 gen_step11 的 _ct_info 一致，靠前的优先）。下游用的来源比写记录的这一步
#   优先级更高时，这一步的层厚根本不会被用到，不该让下游重排。
_CT_PRIORITY = ("THICKNESS_A", "LAYER_THICKNESS", "step8_amset/", "step8.4_amset2d/", "vdw")


def _ct_rank(src):
    s = str(src or "")
    for i, pfx in enumerate(_CT_PRIORITY):
        if s.startswith(pfx):
            return i
    return len(_CT_PRIORITY)          # 旧版本没记来源：按最低优先级（应当核对）


def invalidate_ct_consumers(cwd, ct, step="step8_amset", rel_tol=1e-3):
    """[0015] S8 写出 2d_correction.json 之后：读 c/t 的下游（S8.1）当时用的 c/t 与新值不同 -> 让它重排。

    S8.1 的 needs 里没有 S8（S8 要等 S4/S5/HSE，常常更晚生成），S8.1 先跑时 c/t 按本步结构的 vdW
    层厚现算；S8 的 LAYER_THICKNESS 若被 step.conf 覆盖（如 SS/LS 的 6.73），两边就不是同一个 t。
    这里只在 c/t 真的不同（或下游没记录 c/t）时归档下游完成标记，递归到 S8.3；下游用的来源比
    本步优先（手填 THICKNESS_A、共用 step.conf 的 LAYER_THICKNESS、S8.4 调用时的 S8 记录）的不动。
    返回 [(下游步骤, 改名的文件), ...]。"""
    import json as _json
    import time as _t
    if not ct:
        return []
    cwd = Path(cwd)
    done = []
    tag = "stale-upstream-%s-%s" % (step.replace("/", "_"), _t.strftime("%Y%m%d%H%M%S"))
    for down in CT_CONSUMERS:
        d = cwd / down
        hit = False
        for m in DONE_MARKERS.get(down, ()):
            f = d / m
            if not f.is_file():
                continue
            try:
                rec = _json.loads(f.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                rec = {}
            if rec.get("dim") not in (None, "2d") \
                    or _ct_rank(rec.get("twoD_c_over_t_source")) < _ct_rank(step + "/"):
                continue                    # 3D，或下游用的是更优先的来源（手填 / 共用 step.conf / S8）
            old = rec.get("twoD_c_over_t")
            if old and abs(float(old) - float(ct)) <= rel_tol * abs(float(ct)):
                continue
            f.rename(d / (m + "." + tag))
            done.append((down, m))
            hit = True
            print("[WARN] %s 的 c/t=%s（%s）与 %s 新写的 c/t=%.4f 不同 -> 完成标记已归档（*.%s），会重新排队"
                  % (down, old, rec.get("twoD_c_over_t_source") or "旧版本未记录来源", step, float(ct), tag))
        if hit:
            done += invalidate_downstream(cwd, down, "上游 %s 的 c/t 变了" % step)
    return done


# --------------------------------------------------------------------------
# [patch_walltime V144] 大 VASP 作业的墙时预估（S6 gen 用；tools/vasp_eta.py 对在跑的作业外推）
#   起因：Mo2S3 S6（IBRION=6，20 原子，NFREE=4）跑满 regular 的 24 h 被 SLURM 杀掉（3918831，
#   TIMEOUT 24:00:21），白跑一天；8 月 Zn5O3 的 S2.3/S6、GaAs S5（HSE）也都撞过墙时。
#   IBRION=6 的 SCF 次数在 gen 时就能算出来；每步耗时用 S1 的离子步粗估，跑起来以后用实测外推。
# --------------------------------------------------------------------------
# 没写 --time 时按 QoS 的上限（jzzn 实测：regular 24 h —— 3918831 在 24:00:21 被杀；premium 48 h）。
WALLTIME_BY_QOS_H = {"regular": 24.0, "premium": 48.0}
WALLTIME_DEFAULT_H = 24.0
WALLTIME_WARN_FRAC = 0.8


def _hms_to_h(v):
    """SLURM 时间格式 -> 小时。认不出返回 None。
    无天数：MM / MM:SS / HH:MM:SS；有天数：D-HH / D-HH:MM / D-HH:MM:SS。"""
    v = str(v).strip()
    try:
        if "-" in v:
            d, rest = v.split("-", 1)
            p = [float(x) for x in rest.split(":")] + [0.0, 0.0]
            return int(d) * 24 + p[0] + p[1] / 60.0 + p[2] / 3600.0
        p = [float(x) for x in v.split(":")]
        if len(p) == 1:
            return p[0] / 60.0
        if len(p) == 2:
            return p[0] / 60.0 + p[1] / 3600.0
        if len(p) == 3:
            return p[0] + p[1] / 60.0 + p[2] / 3600.0
    except (ValueError, IndexError):
        pass
    return None


def walltime_limit_h(submit_sh, by_qos=None, default_h=WALLTIME_DEFAULT_H):
    """submit.sh 的墙时上限（小时）与来源说明：--time/-t 优先，否则按 --qos 查表，否则默认。"""
    by_qos = WALLTIME_BY_QOS_H if by_qos is None else by_qos
    try:
        text = Path(submit_sh).read_text(errors="ignore")
    except OSError:
        return default_h, "读不到 submit.sh，按默认 %g h" % default_h
    m = re.search(r"^#SBATCH\s+(?:--time[=\s]|-t\s+)(\S+)", text, re.M)
    if m:
        h = _hms_to_h(m.group(1))
        if h:
            return h, "submit.sh --time=%s" % m.group(1)
    m = re.search(r"^#SBATCH\s+--qos[=\s](\S+)", text, re.M)
    if m and m.group(1) in by_qos:
        return by_qos[m.group(1)], "QoS %s 的上限" % m.group(1)
    return default_h, "submit.sh 没写 --time%s，按默认 %g h" % (
        "（QoS %s 不在表里）" % m.group(1) if m else "", default_h)


def outcar_loop_times(outcar):
    """OUTCAR 里每个 LOOP+（一个离子步/一个有限差分构型）的 real time（秒）。"""
    out = []
    try:
        with open(outcar, errors="ignore") as fh:
            for ln in fh:
                if "LOOP+" in ln and "real time" in ln:
                    try:
                        out.append(float(ln.rsplit("real time", 1)[1].split()[0]))
                    except (IndexError, ValueError):
                        pass
    except OSError:
        pass
    return out


def ibrion6_steps(n_atoms, nfree, isif=3, n_ineq=None):
    """IBRION=6 的 SCF 次数（上界）：初始 1 次 + NFREE × 3 × 不等价原子数 + （ISIF≥3）NFREE × 6 个应变。
    有对称性时 VASP 只位移不等价原子、而且未必三个方向都要，所以这是上界。"""
    n = int(n_ineq if n_ineq else n_atoms)
    nf = int(nfree)
    return 1 + nf * 3 * n + (nf * 6 if int(isif) >= 3 else 0)


def n_inequivalent_atoms(poscar, symprec=1e-3):
    """spglib 的不等价原子数；取不到返回 None（调用方按全部原子算）。"""
    try:
        import numpy as np
        import spglib
        lat, counts, frac, _ = read_poscar_cell(poscar)
        nums = [i for i, c in enumerate(counts) for _ in range(c)]
        ds = spglib.get_symmetry_dataset((lat, frac, nums), symprec=symprec)
        if ds is None:
            return None
        eq = ds["equivalent_atoms"] if isinstance(ds, dict) else ds.equivalent_atoms
        return int(len(np.unique(eq)))
    except Exception:                                                # noqa: BLE001
        return None


def walltime_verdict(n_steps, t_step_s, limit_h, done_steps=0, elapsed_s=0.0, frac=WALLTIME_WARN_FRAC):
    """(预计总小时, 是否超过 frac × 上限)。done/elapsed 给了就按"已用 + 剩余步数 × 每步"算。"""
    remaining = max(int(n_steps) - int(done_steps), 0)
    total_h = (float(elapsed_s) + remaining * float(t_step_s)) / 3600.0
    return total_h, total_h > frac * float(limit_h)


# --------------------------------------------------------------------------
# [patch_carrier_sign V143] AMSET 掺杂符号 -> 载流子：**唯一真源**，读 transport.json 的脚本一律用它。
#   AMSET 0.5.1 FermiDos 文档："A negative doping concentration indicates the majority carriers are
#   electrons (n-type doping); a positive doping concentration indicates holes are the majority
#   carriers (p-type doping)." 本技能 DOPING 也是负值在前、标 n 型。
#   起因：2026-10-02 用户侧手写脚本把 +/− 对调，MoSe2 的电子/空穴 ADP/DPT 表整张标反，来回核对了两轮。
# --------------------------------------------------------------------------
def carrier_of_doping(doping):
    """AMSET 掺杂值 -> "electron"（负，n 型）/ "hole"（正，p 型）/ None（0）。"""
    x = float(doping)
    if x < 0:
        return "electron"
    if x > 0:
        return "hole"
    return None


# --------------------------------------------------------------------------
# [patch_lineage V140] 上游结构同源 + 派生产物新鲜度（S8/S8.4 gen 的闸门；tools/lineage_check.py 也用）
#   起因（Mo2S3，S1 重新弛豫后整条 rerun）：
#   · autozt rerun 先 rm -rf 步骤目录再 gen，V122/V136 那几处"gen 时见到旧产物 -> 下游失效"一个都不触发，
#     S4 的 h5、S7.1 的形变势、S2 画图的 band_summary.json、S8 仍是旧的"完成"；
#   · S2.3 HSE 的 gen 从 S2.2 拷 POSCAR 与 WAVECAR —— S2.2 还没重跑就先 rerun S2.3，HSE 算的是旧结构。
#   各步 POSCAR 都是 S1 CONTCAR 的原样拷贝（relay_poscar；S2.2 读 S2.1 的 CONTCAR；S2.3 拷 S2.2 的
#   POSCAR；S7 的 undeformed/POSCAR 拷本步 POSCAR），所以逐个按数值比对不会误判。
# --------------------------------------------------------------------------
S1_DIRS = ("step1_opt", "step1_std_opt")
LINEAGE_POSCARS = (
    ("step2_bandgap/step2.1_static", ("POSCAR",)),
    ("step2_bandgap/step2.2_pbe", ("POSCAR",)),
    ("step2_bandgap/step2.3_hse", ("POSCAR", "*/POSCAR")),
    ("step3_uniform", ("POSCAR",)),
    ("step3b_uniform_full", ("POSCAR",)),
    ("step5_dielect", ("POSCAR",)),
    ("step6_elastic", ("POSCAR",)),
    ("step7_deform", ("undeformed/POSCAR",)),
)
# (派生步骤, 它的产物, 来源步骤, 来源里的计算输出)：来源输出比产物新 = 产物是用旧来源做的。
#   只看计算输出（VASP 重算才会变），不看 INCAR/POSCAR（gen 原地重写不代表重算）。
LINEAGE_DERIVED = (
    ("step4_wave", "wavefunction.h5", "step3_uniform", ("vasprun.xml", "WAVECAR")),
    ("step4b_wave_full", "wavefunction.h5", "step3b_uniform_full", ("vasprun.xml", "WAVECAR")),
    ("step7b_deform_read", "deformation.h5", "step7_deform", ("*/vasprun.xml", "*/*/vasprun.xml")),
    ("step2_bandgap/step2.2_pbe_plot", "band_summary.json", "step2_bandgap/step2.2_pbe", ("vasprun.xml",)),
    ("step2_bandgap/step2.3_hse_plot", "band_summary.json", "step2_bandgap/step2.3_hse",
     ("vasprun.xml", "*/vasprun.xml")),
    ("step5_dielect_validate", "dielectric_check.json", "step5_dielect", ("OUTCAR",)),
)
LINEAGE_TOL_A = 1e-4        # Å：CONTCAR 与 POSCAR 的文本精度差远小于它；重新弛豫的差别远大于它
LINEAGE_SLACK_S = 60.0


def read_poscar_cell(path):
    """POSCAR/CONTCAR -> (晶格 3x3 Å, 各元素原子数, 分数坐标 n×3, 元素符号或 None)。纯 numpy。
    支持 VASP4/5、负 scale（= 体积）、Selective dynamics、Direct/Cartesian。"""
    import numpy as np
    ln = Path(path).read_text(errors="ignore").splitlines()
    scale = float(ln[1].split()[0])
    lat0 = np.array([[float(x) for x in ln[i].split()[:3]] for i in (2, 3, 4)])
    f = scale if scale > 0 else (abs(scale) / abs(np.linalg.det(lat0))) ** (1.0 / 3.0)
    lat = lat0 * f
    i = 5
    toks = ln[i].split()
    syms = None
    if not all(t.isdigit() for t in toks):
        syms = tuple(toks)
        i += 1
    counts = tuple(int(x) for x in ln[i].split())
    i += 1
    if ln[i].strip()[:1] in ("s", "S"):
        i += 1
    cart = ln[i].strip()[:1] in ("c", "C", "k", "K")
    i += 1
    n = sum(counts)
    xyz = np.array([[float(x) for x in ln[i + j].split()[:3]] for j in range(n)]).reshape(n, 3)
    frac = np.dot(xyz * f, np.linalg.inv(lat)) if cart else xyz
    return lat, counts, frac, syms


ZA_STRUCT_CANDS = ("step3_uniform", "step6_elastic", "step1_opt", "step1_std_opt")


def za_coupling_check(cwd, symprec=1e-2, cands=ZA_STRUCT_CANDS):
    """[0016] 2D：ZA（弯曲声子）的一阶形变耦合是否被对称性禁止 —— 结构有没有水平镜面 σh（z→−z，含滑移）。

    有 σh：电子态对 z→−z 有确定宇称，ZA 的一阶耦合为零 -> 只算面内声学支（DPT 的 LA、S8.4 的面内
      2×2 Christoffel）是对的；
    没有 σh（翘曲结构、Janus、上下不对称堆叠的超晶格）：ZA 一阶耦合可以不为零，DPT 与 AMSET 都没计入，
      结果缺这一项（模型局限，不是能补的数值误差）。
    判据：spglib 对称操作里存在笛卡尔部分 R = I − 2nnᵀ 的操作（n = 层法向，平移任意）。
    返回 {"sigma_h": True/False/None, "structure": 相对路径, "note": 说明}；判不出 sigma_h=None。"""
    import numpy as np
    pos = next((Path(cwd) / d / fn for d in cands for fn in ("POSCAR", "CONTCAR")
                if (Path(cwd) / d / fn).is_file()), None)
    if pos is None:
        return {"sigma_h": None, "structure": None, "note": "找不到结构，判不出 σh"}
    rel = "%s/%s" % (pos.parent.name, pos.name)
    try:
        import spglib
        lat, counts, frac, _ = read_poscar_cell(pos)
        dim, axis, _vac = detect_dimension(str(pos))
        if dim != "2d":
            return {"sigma_h": None, "structure": rel, "note": "不是 2D（%s），不适用" % dim}
        ip = [i for i in range(3) if i != int(axis)]
        n = np.cross(lat[ip[0]], lat[ip[1]])
        n = n / np.linalg.norm(n)
        mirror = np.eye(3) - 2.0 * np.outer(n, n)
        nums = [i for i, c in enumerate(counts) for _ in range(c)]
        ds = spglib.get_symmetry_dataset((lat, frac, nums), symprec=symprec)
        if ds is None:
            return {"sigma_h": None, "structure": rel, "note": "spglib 取不到对称操作"}
        rots = ds["rotations"] if isinstance(ds, dict) else ds.rotations
        inv = np.linalg.inv(lat)
        has = any(np.allclose(inv @ np.asarray(R).T @ lat, mirror, atol=1e-3) for R in rots)
    except (SystemExit, Exception) as e:                             # noqa: BLE001
        return {"sigma_h": None, "structure": rel, "note": "判不出 σh：%s: %s" % (type(e).__name__, e)}
    note = ("有水平镜面 σh：ZA 一阶耦合被对称性禁止，只算面内声学支成立" if has else
            "没有水平镜面 σh：ZA 一阶耦合可以不为零，DPT 与 AMSET 都没计入（结果缺这一项，报数时注明）")
    return {"sigma_h": bool(has), "structure": rel, "symprec": symprec, "note": note}


def cell_deviation(a, b):
    """两个 read_poscar_cell 结果的最大偏差（Å）：晶格矢量逐分量、原子位置（最小像，**扣除整体平移**）。
    原子数或元素对不上 -> inf。

    [V149] 扣整体平移：S3/S3b 的 gen 调 align_origin 把原点刚性平移到 τ=0 的位置（V109，绕开 AMSET
    去对称化 bug），是同一个结构。V140 不扣平移，WS2（S1 原点不在高对称位置）的 S3 被报"差 4.45 Å、旧结构"
    —— 三个原子同移 [0, 0, 4.4494] Å，S8/S8.4 的闸门会据此拦住正常的 rerun。真正的旧结构（重新弛豫前的
    晶格/相对位置）扣完平移照样报。"""
    import numpy as np
    la, ca, fa, sa = a
    lb, cb, fb, sb = b
    if ca != cb or (sa and sb and sa != sb):
        return float("inf")
    dl = float(np.abs(la - lb).max())
    d = fb - fa
    d -= np.rint(d)
    if len(d):
        d -= d[0]                       # 以第一个原子为参照扣平移，再按最小像折回
        d -= np.rint(d)
        d -= d.mean(axis=0)             # 剩下的是小量，取平均再扣一次（对称分摊噪声）
    dp = float(np.linalg.norm(np.dot(d, la), axis=1).max()) if len(d) else 0.0
    return max(dl, dp)



def s1_contcar(cwd):
    cwd = Path(cwd)
    return next((cwd / d / "CONTCAR" for d in S1_DIRS if (cwd / d / "CONTCAR").is_file()), None)


def structure_lineage(cwd, tol=LINEAGE_TOL_A, slack=LINEAGE_SLACK_S):
    """返回 (S1 CONTCAR 路径或 None, 问题列表 [(步骤, 说明), ...])。
    ① 各步 POSCAR 与 S1 CONTCAR 的数值偏差 > tol；② LINEAGE_DERIVED 里产物比来源的计算输出旧。
    读不了的文件、不存在的步骤不判。"""
    cwd = Path(cwd)
    probs = []
    ref_path = s1_contcar(cwd)
    ref = None
    if ref_path is not None:
        try:
            ref = read_poscar_cell(ref_path)
        except Exception as e:                                       # noqa: BLE001
            probs.append((ref_path.parent.name, "读不了 CONTCAR（%s: %s）" % (type(e).__name__, e)))
    if ref is not None:
        for step, pats in LINEAGE_POSCARS:
            for p in sorted({q for pat in pats for q in cwd.glob(step + "/" + pat) if q.is_file()}):
                try:
                    dev = cell_deviation(ref, read_poscar_cell(p))
                except Exception:                                    # noqa: BLE001
                    continue
                if dev > tol:
                    probs.append((step, "%s 与 %s 的结构不同（%s）—— 用的是旧结构"
                                  % (os.path.relpath(str(p), str(cwd)),
                                     os.path.relpath(str(ref_path), str(cwd)),
                                     "原子数/元素不同" if dev == float("inf") else "最大偏差 %.4f Å" % dev)))
                    break
    probs += _derived_probs(cwd, LINEAGE_DERIVED, slack)
    return ref_path, probs


def _derived_probs(cwd, table, slack=LINEAGE_SLACK_S):
    probs = []
    for step, marker, src, pats in table:
        f = cwd / step / marker
        if not f.is_file():
            continue
        outs = [q for pat in pats for q in (cwd / src).glob(pat) if q.is_file()]
        if not outs:
            continue
        newest = max(outs, key=lambda q: q.stat().st_mtime)
        if newest.stat().st_mtime > f.stat().st_mtime + slack:
            probs.append((step, "%s 比来源 %s 旧 —— 来源后来重算过，本步要在来源跑完后重新生成"
                          % (marker, os.path.relpath(str(newest), str(cwd)))))
    return probs


# [V148] S8 之后的派生产物与带隙一致性：只给 tools/lineage_check.py 用 —— S8/S8.4 的 gen 正要重写自己的
#   settings.yaml，拿旧的去拦自己就再也 rerun 不了；S8.1/S8.2/S8.3 都是秒级到分钟级，查出来直接重跑即可。
POST_DERIVED = (
    ("step8.1_boltztrap", "boltztrap_crta.json", "step8.2_dpt", ("dpt_result.json",)),
    ("step8.3_output", "comparison_300K.png", "step8_amset", ("transport.json",)),
    ("step8.3_output", "comparison_300K.png", "step8.4_amset2d", ("transport.json",)),
    ("step8.3_output", "comparison_300K.png", "step8.1_boltztrap", ("boltztrap_crta.json",)),
    ("step8.3_output", "comparison_300K.png", "step8.2_dpt", ("dpt_result.json",)),
    # [V163] DOWNSTREAM 里 always 边的事后复核（以前 lineage_check 只看 S8 之后的派生，S8/S8.4/S8.2 本身新不新
    #   没人查：P1_Mo-MoS2 "4 项全 OK" 其实没验证 transport.json 是用最新的 S5/S6 算的）。link 边见 _linked_probs。
    ("step8_amset", "transport.json", "step5_dielect", ("OUTCAR",)),
    ("step8_amset", "transport.json", "step6_elastic", ("OUTCAR",)),
    ("step8.4_amset2d", "transport.json", "step5_dielect", ("OUTCAR",)),
    ("step8.4_amset2d", "transport.json", "step6_elastic", ("OUTCAR",)),
    ("step8.2_dpt", "dpt_result.json", "step7b_deform_read", ("band_edges.json", "deformation.h5")),
    ("step8.2_dpt", "dpt_result.json", "step6_elastic", ("OUTCAR",)),
    ("step8.2_dpt", "dpt_result.json", "step3_uniform", ("vasprun.xml",)),
    ("step8.1_boltztrap", "boltztrap_crta.json", "step3_uniform", ("vasprun.xml",)),
    ("step8.1_boltztrap", "boltztrap_crta.json", "step2_bandgap/step2.2_pbe_plot", ("band_summary.json",)),
    ("step8.1_boltztrap", "boltztrap_crta.json", "step2_bandgap/step2.3_hse_plot", ("band_summary.json",)),
)
# [V163] link 边（S3/S3b/S4/S4b/S7.1 -> S8/S8.4）：S8 用的是哪一份（S4 还是 S4b）看它目录里的软链，按链接目标比新旧。
LINKED_PRODUCTS = (("step8_amset", "transport.json"), ("step8.4_amset2d", "transport.json"))


def _linked_probs(cwd, slack=LINEAGE_SLACK_S):
    """S8/S8.4 目录里软链进来的输入（wavefunction.h5、deformation.h5、vasprun.xml …）比 transport.json 新 -> 旧结果。
    断链（本地回拉镜像里大文件没拉）不判。"""
    probs = []
    for step, marker in LINKED_PRODUCTS:
        d = Path(cwd) / step
        f = d / marker
        if not f.is_file():
            continue
        t0 = f.stat().st_mtime
        try:
            links = sorted(p for p in d.iterdir() if p.is_symlink())
        except OSError:
            continue
        for p in links:
            try:
                tgt = p.resolve(strict=True)
            except (OSError, RuntimeError):
                continue
            if tgt.is_file() and tgt.stat().st_mtime > t0 + slack:
                probs.append((step, "%s 比链接的输入 %s 旧 —— 来源后来重算过，本步要重新生成"
                              % (marker, os.path.relpath(str(tgt), str(cwd)))))
                break
    return probs


def band_summary_gap(plot_dir):
    """S2 画图目录 -> band_summary.json 的 gap_eV（读不到为 None）。"""
    import json
    try:
        return float(json.loads((Path(plot_dir) / "band_summary.json").read_text(encoding="utf-8"))["gap_eV"])
    except (OSError, KeyError, ValueError, TypeError):
        return None


def post_lineage(cwd, slack=LINEAGE_SLACK_S):
    """[V148] -> [(步骤, 说明)]：① POST_DERIVED 里产物旧于来源；② S8/S8.4 settings.yaml 的带隙与它来源的
    band_summary.json 差 > GAP_TOL_EV（S2 画图后来重新生成过，S8 还是旧带隙）。"""
    cwd = Path(cwd)
    probs = _derived_probs(cwd, POST_DERIVED, slack) + _linked_probs(cwd, slack)
    for down in ("step8_amset", "step8.4_amset2d"):
        old, src = consumer_bandgap(cwd / down)
        if old is None or not src:
            continue
        for plot in ("step2_bandgap/step2.2_pbe_plot", "step2_bandgap/step2.3_hse_plot"):
            if (Path(plot).name + "/") not in src:
                continue
            new = band_summary_gap(cwd / plot)
            if new is not None and abs(new - old) > GAP_TOL_EV:
                probs.append((down, "settings.yaml 的带隙 %.4f eV ≠ %s/band_summary.json 的 %.4f eV —— S2 画图后来"
                                    "重新生成过；--invalidate-from %s 让它重排" % (old, plot, new, plot)))
    return probs


def check_lineage(cwd, enabled=True, label="S8"):
    """gen 用：有问题就列出来；enabled -> 退出（不提交），否则只告警。返回问题列表。"""
    ref, probs = structure_lineage(cwd)
    if ref is None:
        print("[..] 上游同源核对（V140）：没有 S1 的 CONTCAR，只核对派生产物的新旧")
    if not probs:
        print("[OK] 上游同源核对（V140）：各步结构与 S1 一致、派生产物不旧于来源")
        return probs
    print("[%s] ★ 上游同源核对（V140）：%d 处用的是旧结构或旧来源 —— %s 的输入混了两版计算："
          % ("ERROR" if enabled else "WARN", len(probs), label))
    for step, why in probs:
        print("        - %s：%s" % (step, why))
    print("      处理：rerun **最上游**那个结构不同的步骤即可（V153 起它的 gen 会把下游 VASP 步骤整目录归档，"
          "上游跑完后自动重新生成）；派生步骤（S4 / S4b / S7.1 / S2 画图 / S5.1）在来源跑完后 rerun；"
          "都 OK 后再 gen 本步。")
    print("      已经 rerun 过的上游可用 tools/lineage_check.py <材料目录> --invalidate-from <步骤> 归档它的下游完成标记。")
    if enabled:
        sys.exit("[ERROR] %s：上游不同源，不生成（确需照跑：step.conf 写 STRUCTURE_GUARD = false）" % label)
    return probs


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


# ---- [patch_ir_fix] AMSET 密网格不可约 k 点用正确晶格（V125）----------------------------------
#   插件本体：step8.4_amset2d/amset_ir_fix.py。AMSET 给 spglib 的晶格多转置了一次，六方/三方只用上
#   {E, σh}×TR，不可约 k 点多 5.5–5.9 倍（散射白算 5–6 倍）；立方、fcc/bcc 原胞、正交常规胞不受影响。
#   纯簿记修正（原来找到的操作也都是真对称操作），合成体系端到端迁移率一致到 0.05%（标量形变势）。
#   V130 改为默认关：张量形变势下 AMSET 的 ADP 核不协变，代表点≠成员平均；V127 的星平均在 MoS₂ 上
#   仍差 n −1.81% / p −3.12%（ADP，对 V121）。2D 的提速改由 KZ_CAP_2D 提供（结果不依赖不可约映射）；
#   需要 IR_FIX 的提速（3D 六方）时在 step.conf 写 IR_FIX = on，并按 VERIFICATION V130 先对照。
#   取值优先级：step.conf 的 IR_FIX（on/off）> gen 时的环境变量 AZ_IR_FIX > 本常量。
IR_FIX_DEFAULT = False
IR_FIX_ENV = "AZ_IR_FIX"
IR_FIX_PLUGIN = "amset_ir_fix.py"


# ---- [patch_fermi_fix] 求费米能级的稳健化（V132）----------------------------------------------
#   插件本体：step8.4_amset2d/amset_fermi_fix.py，S8 与 S8.4 都装进运行目录并在 python -c 里 import。
#   AMSET 的贪心搜索在低温低掺杂时被双精度舍入困住（WS2 n 型 100 K），容差阶梯还会接受偏差至多 100%
#   的解；插件在原搜索失败时改用二分法、宽松解偏差 > 1% 报错、NaN 报错。原搜索能解时结果逐位不变。
#   开关：作业环境变量 AZ_FERMI_FIX（未设 = 开；0 = 关，仅供对照）。
FERMI_FIX_PLUGIN = "amset_fermi_fix.py"


def ir_fix_setting(conf_value=None):
    """返回 (是否修正不可约 k 点, 来源说明)。conf_value = step.conf 的 IR_FIX 原值。"""
    v = "" if conf_value is None else str(conf_value).strip().lower()
    if v in _ON:
        return True, "step.conf IR_FIX=%s" % conf_value
    if v in _OFF:
        return False, "step.conf IR_FIX=%s" % conf_value
    if v not in ("", "auto", "none", "default"):
        print("[WARN] IR_FIX=%r 不认识（只认 auto/on/off），按 auto 处理" % conf_value)
    e = os.environ.get(IR_FIX_ENV, "").strip().lower()
    if e in _ON:
        return True, "环境变量 %s=%s" % (IR_FIX_ENV, e)
    if e in _OFF:
        return False, "环境变量 %s=%s" % (IR_FIX_ENV, e)
    return bool(IR_FIX_DEFAULT), "出厂默认 IR_FIX_DEFAULT=%s（V130 起默认关）" % IR_FIX_DEFAULT


def ir_fix_cmd_prefix(on):
    """插件在环境变量未设时也会打补丁，所以关的时候必须显式 export 0。"""
    return "export %s=%d; " % (IR_FIX_ENV, 1 if on else 0)


def install_run_plugin(name, out, search_dirs):
    """把运行时插件复制进 AMSET 运行目录。找不到返回 False（调用方决定是否致命）。"""
    import shutil as _sh
    src = next((Path(d) / name for d in search_dirs if (Path(d) / name).is_file()), None)
    dst = Path(out) / name
    if src is None:
        if dst.is_file():
            dst.unlink()
        return False
    if src.resolve() != dst.resolve():
        _sh.copyfile(src, dst)
    return True


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


# ---- [patch_tensor_symmetry-2026-09-29] 物理张量按晶体点群对称化（Neumann 原理，V120）----------
#   2026-09-29 MoS₂ 单层（D3h）生产结果：ADP 迁移率 xx/yy 差 9–10%（n 206.9/189.2、p 1095/995），
#   而六方晶体面内任何二阶张量都必须各向同性。n、p 同比例偏 -> 更像两者共用的输入（弹性张量）
#   带着数值噪声；IMP 主导的总迁移率 xx/yy 只差 0.6% -> 网格/重叠不是主因。
#   VASP IBRION=6 的 TOTAL ELASTIC MODULI（含离子弛豫贡献）在六方胞上常见几个百分点的 C11≠C22。
#   ⇒ 读进来的弹性张量先对点群所有操作取平均（pymatgen Tensor.fit_to_structure）：对称的输入
#     不变，带噪声的输入投影到物理允许的子空间；改动 > TENSOR_SYM_WARN 时告警（上游步欠收敛）。
TENSOR_SYM_WARN = 0.02


def tensor_structure(cwd, dirs):
    """张量所在笛卡尔坐标系对应的结构（优先产出该张量的步骤目录里的 POSCAR）。"""
    from pymatgen.core import Structure
    for d in dirs:
        for n in ("POSCAR", "CONTCAR"):
            p = Path(cwd) / d / n
            if p.is_file():
                try:
                    return Structure.from_file(str(p)), "%s/%s" % (d, n)
                except Exception:                              # noqa: BLE001
                    pass
    return None, None


def symmetrize_elastic_voigt(C, structure, symprec=AMSET_SYMPREC):
    """6x6 标准 Voigt 弹性张量对晶体点群取平均。返回 (C_sym, info)；失败返回 (C, {"error": ...})。"""
    import numpy as np
    try:
        from pymatgen.analysis.elasticity.elastic import ElasticTensor
        C0 = np.asarray(C, float)
        Cs = np.asarray(ElasticTensor.from_voigt(C0).fit_to_structure(structure, symprec=symprec).voigt,
                        float)
    except Exception as e:                                     # noqa: BLE001
        return C, {"error": "%s: %s" % (type(e).__name__, e)}
    scale = float(np.abs(C0).max()) or 1.0
    info = {"max_rel_change": float(np.abs(Cs - C0).max() / scale),
            "C11": float(C0[0, 0]), "C22": float(C0[1, 1]), "C66": float(C0[5, 5]),
            "C11_sym": float(Cs[0, 0]), "C22_sym": float(Cs[1, 1]), "C66_sym": float(Cs[5, 5])}
    return [[round(float(v), 3) for v in r] for r in Cs], info


def rank2_asymmetry(T, structure, symprec=AMSET_SYMPREC):
    """3x3 张量（介电等）偏离点群对称的程度：对称化前后最大相对差。只诊断，不改值。"""
    import numpy as np
    try:
        from pymatgen.core.tensors import Tensor
        T0 = np.asarray(T, float)
        Ts = np.asarray(Tensor(T0).fit_to_structure(structure, symprec=symprec), float)
    except Exception:                                          # noqa: BLE001
        return None
    return float(np.abs(Ts - T0).max() / (float(np.abs(T0).max()) or 1.0))


def symmetrize_elastic_for(cwd, elastic, dirs, enabled=True):
    """gen 用：读进来的弹性张量（6x6）按点群对称化并打印改动；标量/None/关掉时原样返回。"""
    import numpy as np
    if not enabled or elastic is None or np.ndim(elastic) != 2:
        return elastic
    st, src = tensor_structure(cwd, dirs)
    if st is None:
        print("[WARN] patch_tensor_symmetry：取不到结构，弹性张量未对称化")
        return elastic
    Cs, info = symmetrize_elastic_voigt(elastic, st)
    if "error" in info:
        print("[WARN] patch_tensor_symmetry：对称化失败（%s），弹性张量原样使用" % info["error"])
        return elastic
    msg = ("弹性张量按点群对称化（结构 %s）：最大相对改动 %.1f%%；C11/C22 %.2f/%.2f -> %.2f/%.2f，"
           "C66 %.2f -> %.2f GPa" % (src, 100 * info["max_rel_change"], info["C11"], info["C22"],
                                    info["C11_sym"], info["C22_sym"], info["C66"], info["C66_sym"]))
    if info["max_rel_change"] > TENSOR_SYM_WARN:
        print("[WARN] " + msg + " —— 原始张量明显破坏晶体对称（step6_elastic 欠收敛？），"
              "未对称化时会造成本应各向同性方向上的迁移率差异")
    else:
        print("[OK] " + msg)
    return Cs


def deformation_h5_symmetrized(path):
    """V121：形变势 h5 是否已按点群对称化（dp_symmetrize 打的 attrs）。读不了返回 None（不告警）。"""
    try:
        import h5py
        with h5py.File(str(path), "r") as f:
            return bool(int(f.attrs.get("dp_symmetrized", 0)))
    except Exception:                                          # noqa: BLE001
        return None


def full_grid_ibz_count(structure, kpoints, symprec=AMSET_SYMPREC):
    """patch_full_grid_factor（V124）：vasprun 的 k 点是不是一整张 Γ/MP 网格（ISYM=-1 全网格）。

    是 -> 返回 (n_ir, n_full, mesh)：n_ir = 同一网格按晶体对称性（含时间反演）约化后的不可约点数；
    IBZ 列表、不完整网格或判不出来 -> None。
    用途：AMSET 的 interpolation_factor 乘的是 vasprun 的 k 点数（equivalence 数 = nk × factor）；
    全网格 vasprun 的点数是 IBZ 的十几倍，但对称副本不带新信息 —— 同一个出厂 factor 在全网格上
    等于 IBZ 口径的 factor × n_full/n_ir（MoS₂：6912 点 f=4 与 IBZ 434 点 f=61 给同一张 263×263×21）。
    """
    import numpy as np
    k = np.asarray(kpoints, float).reshape(-1, 3)
    n = len(k)
    if n < 2:
        return None
    fr = np.mod(np.round(k, 6), 1.0)
    fr[fr > 1 - 1e-6] = 0.0
    mesh, shift = [], []
    for i in range(3):
        u = np.unique(np.round(fr[:, i], 5))
        m = len(u)
        for s in (0, 1):
            if np.allclose(u, (np.arange(m) + 0.5 * s) / m, atol=2e-5):
                break
        else:
            return None
        mesh.append(m)
        shift.append(s)
    if int(np.prod(mesh)) != n or len(np.unique(np.round(fr, 5), axis=0)) != n:
        return None
    import spglib
    cell = (structure.lattice.matrix, structure.frac_coords, [s.Z for s in structure.species])
    mapping, _ = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=shift, is_time_reversal=True,
                                               symprec=symprec)
    return int(len(np.unique(mapping))), int(n), [int(x) for x in mesh]


# ---- [patch_elastic_guard] 弹性张量力学稳定性（V134）--------------------------------------
#   AMSET 的 ADP 用 Christoffel 方程的本征值（声速²）；弹性张量不正定时有负/零本征值，ADP 散射率与
#   迁移率没有意义，而 AMSET 照跑不报错（实测 Mo2S3：settings.yaml 里 C44 = −40.16 GPa）。
#   判据（标准 Voigt 顺序 XX YY ZZ YZ XZ XY，取对称部分）：
#     3D —— 整个 6×6 正定；
#     2D —— 只看面内块 (XX, YY, XY) 的 3×3 正定（面外分量是真空伪影，另有处理）。
#   另：2D 的面内剪切 C66 比 C11/C22 小两个量级以上、而 YZ 槽位却很大 —— 典型的"VASP 顺序没重排"
#   （XY 落在第 4 位），只告警（各向异性很强的材料理论上也可能）。
ELASTIC_SOFT_C66 = 0.02       # 2D：C66 < 0.02·max(C11, C22) 告警
ELASTIC_NEAR_SINGULAR = 1e-3  # 最小本征值 < 1e-3·最大对角元 告警（接近奇异，声速≈0）


def elastic_stability(elastic, is_2d=False):
    """返回 dict(ok, kind, eigs, block, labels, near_singular, soft_c66, misorder_hint)；elastic 为 None 返回 None。"""
    import numpy as _np
    if elastic is None:
        return None
    a = _np.asarray(elastic, dtype=float)
    if a.size == 1:
        v = float(a.reshape(-1)[0])
        return {"ok": v > 0, "kind": "scalar", "eigs": [v], "block": [[v]], "labels": ["C"],
                "near_singular": False, "soft_c66": False, "misorder_hint": False}
    if a.shape != (6, 6):
        return {"ok": False, "kind": "shape", "eigs": [], "block": a.tolist(), "labels": [],
                "near_singular": False, "soft_c66": False, "misorder_hint": False}
    sym = 0.5 * (a + a.T)
    if is_2d:
        idx, labels = [0, 1, 5], ["XX", "YY", "XY"]
    else:
        idx, labels = list(range(6)), ["XX", "YY", "ZZ", "YZ", "XZ", "XY"]
    blk = sym[_np.ix_(idx, idx)]
    eigs = _np.linalg.eigvalsh(blk)
    dmax = float(_np.max(_np.abs(_np.diag(blk)))) or 1.0
    soft = misorder = False
    if is_2d:
        cmax = max(sym[0, 0], sym[1, 1])
        soft = bool(cmax > 0 and 0 < sym[5, 5] < ELASTIC_SOFT_C66 * cmax)     # ≤ 0 由正定性判据处理
        misorder = bool(soft and abs(sym[3, 3]) > 0.1 * cmax and abs(sym[3, 3]) > 5 * sym[5, 5])
    return {"ok": bool(eigs.min() > 0), "kind": "2d" if is_2d else "3d", "eigs": [float(x) for x in eigs],
            "block": blk.tolist(), "labels": labels,
            "near_singular": bool(eigs.min() > 0 and eigs.min() < ELASTIC_NEAR_SINGULAR * dmax),
            "soft_c66": soft, "misorder_hint": misorder}


def check_elastic_stability(elastic, is_2d=False, enabled=True, used=True, label=""):
    """gen 里调：不正定 -> 退出（enabled 且 ADP 要用它）；否则只告警。返回 elastic_stability 的结果。"""
    r = elastic_stability(elastic, is_2d)
    if r is None:
        return None
    tag = "[%s] " % label if label else ""
    what = {"2d": "面内块 (XX, YY, XY)", "3d": "6×6", "scalar": "标量"}.get(r["kind"], "形状不对")
    eig_s = ", ".join("%.3f" % x for x in r["eigs"])
    if r["misorder_hint"] or r["soft_c66"]:
        a = [list(map(float, row)) for row in elastic]
        print("[WARN] %s弹性：面内剪切 C66 = %.3f GPa 只有 max(C11, C22) 的 %.1f%%%s —— 先核对 step6_elastic/OUTCAR "
              "的 TOTAL ELASTIC MODULI 表头（XX YY ZZ XY YZ ZX）是否已重排成标准 Voigt（gen 日志应有"
              "\"[OK] 弹性常数：来源顺序 … 已重排\"）。"
              % (tag, a[5][5], 100 * a[5][5] / max(a[0][0], a[1][1]),
                 "，而 YZ 槽位 |C44| = %.3f GPa 很大（像是 XY 落在了第 4 位）" % abs(a[3][3]) if r["misorder_hint"] else ""))
    if r["ok"]:
        print("[OK] %s弹性张量%s正定（本征值 %s GPa）" % (tag, what, eig_s))
        if r["near_singular"]:
            print("[WARN] %s弹性张量%s最小本征值很小（%s）—— 对应方向声速接近 0，ADP 会被放大" % (tag, what, eig_s))
        return r
    rows = "\n".join("          %s" % "  ".join("%9.3f" % x for x in row) for row in r["block"])
    msg = ("%s弹性张量%s不正定（本征值 %s GPa）：\n%s\n"
           "        AMSET 的 ADP 用 Christoffel 本征值（声速²），这里有负/零本征值，ADP 散射率与迁移率没有意义。\n"
           "        常见原因：① step6_elastic 的结构没弛豫到极小 / Γ 点有软模 -> 离子弛豫贡献把剪切模量压成负值"
           "（看 OUTCAR 的 ELASTIC MODULI CONTR FROM IONIC RELAXATION，收紧 EDIFFG 重新弛豫后重算 S6）；\n"
           "                  ② 弹性表的行列顺序没按表头重排（XX YY ZZ XY YZ ZX -> 标准 Voigt）；\n"
           "                  ③ MANUAL_ELASTIC 填错。\n"
           "        确需照跑（结果的 ADP 不可信）：本步 step.conf 写 ELASTIC_GUARD = false。"
           % (tag, what, eig_s, rows))
    if enabled and used:
        sys.exit("[ERROR] " + msg)
    print("[WARN] " + msg + ("" if used else "\n        （本次散射机制不含 ADP，弹性张量用不上 —— 只告警）"))
    return r


def report_dielectric_symmetry(cwd, eps_inf, eps_static, dirs):
    """介电张量偏离点群对称的诊断（只报告：2D 路径会重新读 OUTCAR，这里改值传不过去）。"""
    st, src = tensor_structure(cwd, dirs)
    if st is None:
        return
    for name, T in (("ε∞", eps_inf), ("ε₀", eps_static)):
        if T is None:
            continue
        a = rank2_asymmetry(T, st)
        if a is not None and a > TENSOR_SYM_WARN:
            print("[WARN] patch_tensor_symmetry：%s 偏离点群对称 %.1f%%（结构 %s）—— step5_dielect 欠收敛？"
                  % (name, 100 * a, src))


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


def amset_env_name(cwd=None, fallback=None):
    """AMSET 运行环境名：读 step.conf 的 [params] AMSET_ENV。

    AMSET_ENV 由 autozt 从 setting/<集群>.yaml 的 amset_env 注入 step.conf
    （autozt/report.py，标签 [cluster:<hpc>]），项目级 step.conf 可覆盖。
    gen 的 cwd 是材料目录、step.conf 就在那儿；缺失或解析失败时回退
    fallback（本地直跑 / 旧项目）。用于渲染 submit_amset.tpl 的 {{AMSET_ENV}}。

    2026-09-30（0.4.19 老环境已删）：fallback 默认改为 None——读不到 AMSET_ENV
    即报错退出，绝不静默落一个可能错误的环境名。各集群 0.5.1 环境名不同
    （jzzn/hanhai25=amset051、a800=amset_env、3090/hfeshell=amset），写死任何一个
    都会在别的集群误触。
    """
    import sys as _sys
    env = None
    try:
        import stepconf as _sc
        base = Path(cwd) if cwd else Path(".")
        txt = (base / _sc.CONF_NAME).read_text(encoding="utf-8-sig")
        for k, v, _t in _sc.parse(txt, _sc.CONF_NAME).get("params", []):
            if k.upper() == "AMSET_ENV" and v:
                env = v
                break
    except Exception:
        pass
    if env:
        if env in OLD_AMSET_ENVS:
            _sys.exit("[ERROR] step.conf 的 AMSET_ENV=%s 是 AMSET 0.4.19 的旧环境（%s）。查材料/项目的 "
                      "templates/step.conf 和 project_setting/hpc.yaml 的 amset_env，改成本集群的 0.5.1 环境名"
                      "（jzzn/hanhai25=amset051、a800=amset_env、3090/hfeshell=amset）。" % (env, OLD_AMSET_ENVS[env]))
        return env
    if fallback:
        print("[WARN] step.conf 里没读到 AMSET_ENV，回退到显式指定的 %r。" % fallback,
              file=_sys.stderr)
        return fallback
    _sys.exit("[ERROR] step.conf 里没读到 AMSET_ENV，无法确定本集群的 amset 环境名。"
              "请在 step.conf 写 AMSET_ENV=<本集群 0.5.1 环境名>（jzzn/hanhai25=amset051、"
              "a800=amset_env、3090/hfeshell=amset）。")


# [V159] AMSET 0.4.19 时代的环境名（已删）：0.4.19 的形变势只有 0.5.1 的一半（V75），绝不能悄悄用上。
OLD_AMSET_ENVS = {"amset_clean": "jzzn 上 0.4.19 的环境，2026-09-30 已删"}
_ACTIVATE_RE = re.compile(r"\b(?:conda|mamba|micromamba)\s+(?:activate|run\s+(?:-n|--name))\s+[\"']?([^\s\"';&|)]+)"
                          r"|\bsource\s+activate\s+[\"']?([^\s\"';&|)]+)")


def amset_submit_envs(text):
    """提交脚本里（非注释行）激活/调用的 conda 环境名；${VAR} 这类 shell 变量不算。"""
    out = []
    for ln in text.splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        for m in _ACTIVATE_RE.finditer(s):
            name = m.group(1) or m.group(2)
            if name and not name.startswith("$") and name not in out:
                out.append(name)
    return out


def check_amset_submit(raw, rendered, env, tpl="submit_amset.tpl"):
    """[V159] 渲染后的 AMSET 提交脚本必须激活 step.conf 的 AMSET_ENV。

    起因：P1_Mo-MoS2 的 S4_wave 作业报 EnvironmentNameNotFound: amset_clean。8 月 init 时拷进
    project_setting/templates/ 的 submit_amset.tpl 是 09-22 之前的版本，写死 conda activate amset_clean
    （0.4.19），没有 {{AMSET_ENV}} 占位符；项目级模板优先于集群/技能模板，gen 只查"占位符有没有残留"，
    于是照常生成、提交。环境还在的集群上会悄悄用 0.4.19 算。
    返回告警列表；激活的环境不对时直接退出。
    """
    found = amset_submit_envs(rendered)
    bad = [n for n in found if n != env]
    stale_hint = ("这份 %s 多半是旧的项目级副本（<材料>/<技能>/templates/ 或 project_setting/templates/ 里的优先于"
                  "集群 setting/<集群>/templates/ 和技能模板）。把它改名为 submit_amset.tpl.stale-<日期>，或者把"
                  "写死的环境名改成 {{AMSET_ENV}}；然后 retry 本步。" % tpl)
    if bad:
        old = [n for n in bad if n in OLD_AMSET_ENVS]
        sys.exit("[ERROR] %s 激活的 AMSET 环境是 %s，而 step.conf 的 AMSET_ENV=%s%s。%s"
                 % (tpl, "/".join(bad), env,
                    "（%s 是 AMSET 0.4.19 的旧环境）" % "/".join(old) if old else "", stale_hint))
    warns = []
    if "{{AMSET_ENV}}" not in raw:
        warns.append("[WARN] %s 没有 {{AMSET_ENV}} 占位符%s —— 换集群或换环境时不会跟着变。%s"
                     % (tpl, "（写死的 %s 与 AMSET_ENV 一致，这次照常生成）" % "/".join(found) if found
                        else "，也没有 conda activate 行，作业用的是 PATH 里的 amset", stale_hint))
    for w in warns:
        print(w, file=sys.stderr)
    return warns


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


# [V181] 结果文件里写明泛函。几何、能带、弹性、形变势全用 step1 定下的同一个泛函（下游继承 FUNC=），
#   但 dpt_result.json / 2d_correction.json / 8.3 汇总以前都不记它：agent 把 GGA=PS 读成"PBE"，
#   和 PBE 文献逐项比时把泛函差异（PBEsol 晶格更小 -> C_2D 更硬、m* 和谷间能差跟着变）当成了别的原因。
FUNC_LABEL = {"pbe": "PBE（GGA=PE，无色散修正）",
              "pbesol": "PBEsol（GGA=PS，无色散修正）",
              "pbe-d3": "PBE+D3(BJ)（GGA=PE，IVDW=12）"}
FUNC_COMPARE_NOTE = {
    "pbesol": "和文献比较前先核对文献用的泛函：GGA=PS 是 PBEsol，不是 PBE。PBEsol 的晶格通常比 PBE "
              "小约 1%，C_2D 随之偏硬，m*、形变势和谷间能差也会变，不能直接当作 PBE 结果比。",
    "pbe": "和文献比较前先核对文献用的泛函。",
    "pbe-d3": "和文献比较前先核对文献用的泛函：D3 只改总能、力和应力（几何），不改给定几何下的能带。",
}


def material_functional(mat_dir):
    """[V181] 本材料用的泛函：{"func", "label", "source", "note"}；读不到时 func=None。

    先读 step1 的 workflow_method.txt（FUNC=，下游各步继承的就是它），没有再从 step1 的 INCAR 反推。
    """
    base = Path(mat_dir)
    for how, fn in (("workflow_method.txt", lambda d: read_method_func(base / d / METHOD_FILE)),
                    ("INCAR（从 GGA/IVDW 反推）", lambda d: sniff_func_from_incar(base / d / "INCAR"))):
        for d in ("step1_opt", "step1"):
            f = fn(d)
            if f:
                return {"func": f, "label": FUNC_LABEL[f], "source": "%s/%s" % (d, how),
                        "note": FUNC_COMPARE_NOTE[f]}
    return {"func": None, "label": "未知（step1 的 workflow_method.txt 和 INCAR 都读不到）",
            "source": None, "note": "泛函未知：先查 step1 的 INCAR（GGA/IVDW），再和文献比。"}


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


# --------------------------------------------------------------------------
# [V165] 杂化泛函参数核对：Si 的 S2.3 是 HFSCREEN = 0.11（HSE06 的 0.2 Å⁻¹ 换成 bohr⁻¹ 才是 0.106），
#   屏蔽变弱 -> 带隙 1.340 eV（标准 HSE06 是 1.091），进了 band_summary、S8 的剪刀差和手稿的表，一路没人报警。
#   VASP 的 HFSCREEN 单位是 Å⁻¹：HSE06/HSEsol = 0.2，HSE03 = 0.3。
# --------------------------------------------------------------------------
HYBRID_STD = {(0.25, 0.2): "HSE06/HSEsol", (0.25, 0.3): "HSE03"}


def hybrid_params(incar):
    """INCAR -> {"lhfcalc", "aexx", "hfscreen", "label", "standard", "warning"}；读不了或没开杂化 -> None。"""
    try:
        text = Path(incar).read_text(errors="ignore")
    except OSError:
        return None
    kv = {}
    for ln in text.splitlines():
        ln = ln.split("#", 1)[0].split("!", 1)[0]
        for part in ln.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                kv[k.strip().upper()] = v.strip()
    if not kv.get("LHFCALC", "").upper().lstrip(".").startswith("T"):
        return None

    def _f(key, default):
        try:
            return float(kv[key].split()[0])
        except (KeyError, ValueError, IndexError):
            return default
    aexx, hfs = _f("AEXX", 0.25), _f("HFSCREEN", 0.0)
    std = HYBRID_STD.get((round(aexx, 4), round(hfs, 4)))
    out = {"lhfcalc": True, "aexx": aexx, "hfscreen": hfs, "standard": bool(std),
           "label": std or ("PBE0 型（不屏蔽）" if hfs == 0 else "非标准杂化（AEXX=%g, HFSCREEN=%g Å⁻¹）" % (aexx, hfs)),
           "warning": None}
    if not std and hfs > 0:
        hint = ""
        if 0.09 <= hfs <= 0.12:
            hint = "；%g 像是 HSE06 的 bohr⁻¹ 值（0.106），VASP 的 HFSCREEN 单位是 Å⁻¹，HSE06 应写 0.2" % hfs
        out["warning"] = ("★ 杂化参数不是标准 HSE06（AEXX=0.25, HFSCREEN=0.2 Å⁻¹）：AEXX=%g, HFSCREEN=%g%s。"
                          "带隙会随之改变；若不是有意为之，查 step.conf 的 [incar] 覆盖" % (aexx, hfs, hint))
    return out
