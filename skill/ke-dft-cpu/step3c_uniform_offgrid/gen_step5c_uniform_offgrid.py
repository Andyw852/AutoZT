#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step5c_uniform_offgrid.py —— 离网格 k 点的两段式计算（step3c_uniform_offgrid）。

【目的】直接检验 AMSET 的**系数插值**是否可靠（VERIFICATION V22.9）：
  比 "AMSET 用粗网格插值出的系数" 与 "DFT 在这些离网格 k 点上真算的系数"。

【两段式，一个作业，**不动 S3b**】
  stage1_scf/  ：自洽计算。KPOINTS/INCAR 直接**复制 step3_uniform**（同 47x47x3、
                 同 ENCUT/泛函/NBANDS/FFT 网格），只强制 LCHARG = .TRUE. 写出 CHGCAR。
                 ★ 用不可约 k 点做自洽即可 —— 收敛后的电荷密度与全网格相同，
                   成本只有全网格的约 1/16；因此**不需要重跑 S3b**。
                 ★ 跑完请核对总能量与 step3_uniform 的差 < 1e-4 eV。
  stage2_nscf/ ：非自洽。ICHARG = 11 读 stage1 的 CHGCAR（电势固定），ISYM = -1，
                 只算 14 个离网格 k 点，LWAVE = .TRUE. 供 amset wave 抽系数。
  跑完把 stage2 的 WAVECAR/vasprun.xml/OUTCAR 复制到步骤根目录，供下游与检查使用。

【口径一致（用户 2026-09-17 要求）】INCAR 整份复制自 step3_uniform，只打补丁；
  NBANDS 从 step3_uniform/OUTCAR 读出来显式写死（不靠 VASP 自动决定，否则
  k 点数一变 VASP 会选不同的 NBANDS，平面波集合就不是同一套）。
"""
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ke_common as kc
import stepconf  # noqa: E402
from dim_common import resolve_tpl  # noqa: E402

# =========================== 可改参数区 ===========================
OUTDIR_NAME = "step3c_uniform_offgrid"
SRC_DIR     = "step3_uniform"        # KPOINTS/INCAR/POSCAR/POTCAR 与 NBANDS 来源
MESH        = 47                      # 面内网格数
STEP_LABEL  = "S3c_offgrid"
K_CELL      = (16, 15)                # K 点所在网格单元（MoS2: K = 16/47, 15/47）
KZ_EXTRA    = 1.0 / 6.0               # 粗网格 kz 只有 {0, 1/3}，取中点
Q_FRAC      = 0.17                    # Q 谷（Gamma-K 线约 0.51 处 -> 0.51/3 ≈ 0.17）
CELL_OFFSETS = [(0, 0), (1, 0), (0, 1), (1, 1),
                (-1, 0), (0, -1), (-1, -1), (-1, 1)]
# =================================================================


def _kpoints_list(mesh):
    """14 个离网格 k 点：K 点本身 + K 谷 8 个单元中点 + 3 个 kz=1/6 + 2 个 Q 谷。"""
    kx0, ky0 = K_CELL
    pts = []
    # ① K 点本身（47 点网格不经过 K：47/3 ≈ 15.67）
    pts.append((1.0 / 3.0, 1.0 / 3.0, 0.0, "K 点本身"))
    # ② K 谷 8 个网格单元中点（kz = 0）
    for dx, dy in CELL_OFFSETS:
        pts.append(((kx0 + dx + 0.5) / mesh, (ky0 + dy + 0.5) / mesh, 0.0, "K 谷单元中点"))
    # ③ kz = 1/6（kz 方向粗网格只有 {0,1/3}，这是最可能出问题的方向）
    for fx, fy in ((1.0 / 3.0, 1.0 / 3.0),
                   ((kx0 + 0.5) / mesh, (ky0 + 0.5) / mesh),
                   ((kx0 + 1.5) / mesh, (ky0 + 1.5) / mesh)):
        pts.append((fx, fy, KZ_EXTRA[0] if isinstance(KZ_EXTRA, tuple) else KZ_EXTRA, "kz=1/6"))
    # ④ Q 谷（比 K 谷高约 0.12 eV，300 K 以上参与输运）
    pts.append((Q_FRAC, Q_FRAC, 0.0, "Q 谷"))
    pts.append((Q_FRAC, Q_FRAC, KZ_EXTRA, "Q 谷 kz=1/6"))
    return pts


def read_nbands(outcar: Path):
    if not outcar.is_file():
        return None
    t = outcar.read_text(errors="ignore")
    m = re.findall(r"NBANDS=\s*(\d+)", t)
    return int(m[-1]) if m else None


def patch_incar(text, sets, removes=()):
    out = []
    seen = set()
    for ln in text.splitlines():
        key = ln.split("=")[0].strip().upper() if "=" in ln and not ln.strip().startswith("#") else None
        if key in removes:
            continue
        if key in sets:
            out.append("%-7s= %s" % (key, sets[key]))
            seen.add(key)
            continue
        out.append(ln)
    for k, v in sets.items():
        if k not in seen:
            out.append("%-7s= %s" % (k, v))
    return "\n".join(out) + "\n"


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    src = cwd / SRC_DIR
    for f in ("INCAR", "KPOINTS", "POSCAR", "POTCAR"):
        if not (src / f).is_file():
            sys.exit("[ERROR] %s 缺失：%s（需要 step3_uniform 已完成）" % (f, src / f))

    nbands = read_nbands(src / "OUTCAR")
    if nbands is None:
        sys.exit("[ERROR] 从 %s 读不到 NBANDS —— 口径一致要求显式写死 NBANDS，不能靠自动。" % (src / "OUTCAR"))
    print("[OK] 从 %s/OUTCAR 读到 NBANDS = %d（两段都显式写死）" % (SRC_DIR, nbands))

    base = (src / "INCAR").read_text(encoding="utf-8")
    s1 = out / "stage1_scf"
    s2 = out / "stage2_nscf"
    for d in (s1, s2):
        d.mkdir(exist_ok=True)
        for f in ("POSCAR", "POTCAR", "KPOINTS"):
            dp = d / f
            if dp.is_symlink() or dp.exists():
                dp.unlink()
            dp.symlink_to(os.path.relpath(src / f, d))

    # ---- 段 1：自洽（不可约 k 点，同 S3），只加 LCHARG ----
    (s1 / "INCAR").write_text(
        patch_incar(base, {"LCHARG": ".TRUE.", "LWAVE": ".FALSE.", "NBANDS": nbands}),
        encoding="utf-8", newline="\n")
    # ---- 段 2：非自洽（离网格点），ICHARG=11 + ISYM=-1 ----
    (s2 / "INCAR").write_text(
        patch_incar(base, {"ICHARG": "11", "ISYM": "-1", "LCHARG": ".FALSE.",
                           "LWAVE": ".TRUE.", "NBANDS": nbands, "NELM": 200, "EDIFF": "1E-7"}),
        encoding="utf-8", newline="\n")
    chg = s2 / "CHGCAR"
    if chg.is_symlink() or chg.exists():
        chg.unlink()
    chg.symlink_to(os.path.relpath(s1 / "CHGCAR", s2))
    print("[OK] 段1 INCAR（LCHARG=.TRUE.）/ 段2 INCAR（ICHARG=11, ISYM=-1）已按 step3_uniform 打补丁")

    # ---- 段 2 的离网格 KPOINTS ----
    pts = _kpoints_list(MESH)
    lines = ["off-grid k-points near K / kz=1/6 / Q  (NOT on the %dx%dx3 mesh)" % (MESH, MESH),
             str(len(pts)), "Reciprocal"]
    lines += ["  %.10f  %.10f  %.10f   1.0" % (p[0], p[1], p[2]) for p in pts]
    # 注意：显式模式覆盖 POSCAR 链接来的 KPOINTS（这里直接写文件，不软链）
    if (s2 / "KPOINTS").is_symlink():
        (s2 / "KPOINTS").unlink()
    (s2 / "KPOINTS").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print("[OK] 段2 KPOINTS：%d 个离网格点 —— %s" % (len(pts), "、".join(sorted({p[3] for p in pts}))))

    # ---- 一个作业跑两段 ----
    tpl = resolve_tpl(Path(__file__).resolve().parent, "submit_std", "2d")
    text = tpl.read_text(encoding="utf-8")
    m = re.search(r"^mpirun .*vasp.*$", text, re.M)
    if not m:
        sys.exit("[ERROR] 提交模板里找不到 mpirun 行：%s" % tpl)
    vasp = m.group(0).split()[-1]
    two = ("VASP=%s\n"
           "echo \"===== 段 1/2：自洽（不可约 k 点，写 CHGCAR）=====\"\n"
           "cd stage1_scf && mpirun -np $SLURM_NTASKS $VASP ; cd ..\n"
           "echo \"===== 段 2/2：非自洽（%d 个离网格 k 点，ICHARG=11）=====\"\n"
           "cd stage2_nscf && mpirun -np $SLURM_NTASKS $VASP ; cd ..\n"
           "cp -f stage2_nscf/WAVECAR . 2>/dev/null\n"
           "cp -f stage2_nscf/vasprun.xml . 2>/dev/null\n"
           "cp -f stage2_nscf/OUTCAR OUTCAR 2>/dev/null\n"
           "cp -f stage2_nscf/OSZICAR OSZICAR 2>/dev/null\n"
           "cp -f stage1_scf/OUTCAR OUTCAR.scf 2>/dev/null\n"
           "echo \"[step3c] 两段完成：OUTCAR.scf 是自洽段，OUTCAR 是非自洽段\"\n"
           % (vasp, len(pts)))
    text = text[:m.start()] + two + text[m.end():]
    submit = out / "submit.sh"
    submit.write_text(text, encoding="utf-8", newline="\n")
    kc.patch_submit_jobname(submit, kc.new_jobname(cwd, STEP_LABEL))
    stepconf.apply_submit(submit, stepconf.read_submit(stepconf.CONF_NAME, used_incar=True))
    print("[DONE] %s：两段式就绪（stage1_scf + stage2_nscf，一个作业），不动 S3b" % OUTDIR_NAME)


if __name__ == "__main__":
    main()
