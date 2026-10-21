#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""nscf_dp_check.py —— 形变势能不能用"粗网格自洽 + 密网格非自洽"算（不改 S7 流程的独立验证）。

背景：3D 的 S7 用 vaspkit 的粗网格（Si 11³），带边谷不在格点上时 D(k) 只能插值（Si 电子 E1 曾只得 0.08 eV）。
  直接把 S7 加密到 S3 的网格（DK_MAX_3D = 0.06，Si 34³）k 点多约 30 倍。候选做法：自洽仍用粗网格（只要电荷密度收敛），
  再用这份 CHGCAR 在密网格上做非自洽（ICHARG = 11）取本征值。要先核对四件事，本工具量其中三件：
    1. 芯能级参考：AMSET 用 OUTCAR 的 "average (electrostatic) potential at core"（各原子平均，取第一次出现，
       与 amset.deformation.io.get_reference_energy 同一读法）。非自洽的 OUTCAR 里有没有、与自洽差多少；
    2. 形变势：同一密网格上，非自洽 vs 直接自洽的 D_n(k) = (E0 − E_s − (ref0 − ref_s))/δ_s、±δ 平均
       （AMSET 的写法，带符号、不取 |·|），带边附近差多少；顺带给粗网格自洽在离谷最近格点上的值作对照；
    3. 墙钟时间。
  （第四件 —— 自洽段多写 CHGCAR 怎样不触发 S7 整批重算 —— 是改流程时的事，这里不涉及。杂化泛函不能用 ICHARG = 11，
   make 时 INCAR 有 LHFCALC = .TRUE. 直接拒绝。）

每个构型三段（同一个作业里顺序跑）：
  scf_coarse  S7 的 INCAR 原样（只改 LCHARG = .TRUE.、LWAVE = .FALSE.）+ S7 的粗 KPOINTS -> CHGCAR
  nscf_dense  ICHARG = 11 读 scf_coarse 的 CHGCAR + 密网格（Γ 心，与 S7 同为 ISYM = 0）
  scf_dense   与 scf_coarse 相同的 INCAR（LCHARG = .FALSE.）+ 密网格：参照
构型：undeformed + 一个应变图样的 +δ、−δ（默认 xx；从各目录 POSCAR 的晶格算 Green–Lagrange 应变认图样）。
用的是 S7 顶层（刚性）构型的输入 —— 比较的是算法，不是离子弛豫口径；S7 只要 gen 过（有输入）即可，不需要已经算完。

材料目录只读：输入一律**拷贝**到 <out_dir>（不得在 step7_deform 里面）。本工具不提交作业：
  python nscf_dp_check.py make <step7_deform 目录> <out_dir> [--pattern xx] [--dk 0.06 | --mesh N N N] [--time 24:00:00]
      -> <out_dir>/<构型>/{scf_coarse,nscf_dense,scf_dense}/ 与 run.sh（bash run.sh <VASP 命令>）；
         S7 构型里有 submit.sh 时另写 <out_dir>/<构型>/submit.sh（同一个 SBATCH 头与环境、作业名改成 nscfdp-<构型>，
         --time 可覆盖墙钟上限），在 <out_dir>/<构型> 里 sbatch submit.sh 即可。
  python nscf_dp_check.py analyze <out_dir> [--ewin 0.1] [--tol-ref 0.001] [--tol-d 0.02]
      -> 打印对比、写 <out_dir>/nscf_dp_report.json。退出码 0 = 三项都过（非自洽有芯能级参考且与自洽差 ≤ tol-ref eV；
         带边 D 的最大差 ≤ tol-d × 带边 |D| 的最大值），1 = 有不过的，2 = 算不了（缺文件 / 没跑完）。
只用 numpy（AMSET 的 conda 环境或任何有 numpy 的 python 都行）。
"""
import argparse
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np

META = "nscf_dp_meta.json"
REPORT = "nscf_dp_report.json"
STAGES = ("scf_coarse", "nscf_dense", "scf_dense")
PAT_NAMES = {(0, 0): "xx", (1, 1): "yy", (2, 2): "zz", (1, 2): "yz", (0, 2): "xz", (0, 1): "xy"}
DONE_TAG = "General timing and accounting informations"
RUN_SH = r'''#!/bin/bash
# nscf_dp_check：scf_coarse -> nscf_dense（读它的 CHGCAR）-> scf_dense。用法：bash run.sh <VASP 命令>
set -uo pipefail
[ $# -ge 1 ] || { echo "用法：bash run.sh <VASP 命令>，例如 bash run.sh mpirun -np 24 vasp_std" >&2; exit 2; }
here="$(cd "$(dirname "$0")" && pwd)"
run() {
  cd "$here/$1" || exit 1
  t0=$(date +%s)
  "${@:2}" > vasp.log 2>&1
  rc=$?
  echo "$1 $(( $(date +%s) - t0 )) $rc" >> "$here/timing.txt"
  grep -q "General timing and accounting informations" OUTCAR 2>/dev/null || { echo "[ERROR] $1 没跑完（rc=$rc）" >&2; return 1; }
}
: > "$here/timing.txt"
run scf_coarse "$@" || exit 1
[ -s "$here/scf_coarse/CHGCAR" ] || { echo "[ERROR] scf_coarse 没写出 CHGCAR" >&2; exit 1; }
cp "$here/scf_coarse/CHGCAR" "$here/nscf_dense/CHGCAR"
run nscf_dense "$@"
run scf_dense "$@"
'''


# ---------------------------------------------------------------------------
# 读写 VASP 文件（只依赖 numpy）
# ---------------------------------------------------------------------------
def read_lattice(poscar):
    """POSCAR -> 晶格（行向量，Å）。缩放因子为负 = 目标体积。"""
    ln = Path(poscar).read_text().splitlines()
    s = float(ln[1].split()[0])
    L = np.array([[float(x) for x in ln[i].split()[:3]] for i in (2, 3, 4)])
    if s < 0:
        return L * (abs(s) / abs(np.linalg.det(L))) ** (1.0 / 3.0)
    return L * s


def green_lagrange(L0, L):
    """与 AMSET calculate_deformation + pymatgen Deformation.green_lagrange_strain 同一口径。"""
    F = (np.linalg.inv(L0) @ L).T
    return 0.5 * (F.T @ F - np.eye(3))


def strain_pattern(E):
    """-> ((a, b), 带符号幅度) 或 None（不是单一分量）。"""
    m = float(np.max(np.abs(E)))
    if m < 1e-8:
        return None
    ab = {(int(min(i, j)), int(max(i, j))) for i, j in np.argwhere(np.abs(E) > 0.25 * m)}
    if len(ab) != 1:
        return None
    a, b = ab.pop()
    return (a, b), float(E[a, b])


def set_tag(text, key, value):
    pat = re.compile(r"^[ \t]*%s[ \t]*=.*$" % key, re.MULTILINE)
    line = "%s = %s" % (key, value)
    return pat.sub(line, text) if pat.search(text) else text.rstrip() + "\n" + line + "\n"


def incar_value(text, key):
    m = re.search(r"^[ \t]*%s[ \t]*=[ \t]*([^#!\n]*)" % key, text, re.MULTILINE | re.IGNORECASE)
    return m.group(1).strip() if m else None


def is_true(v):
    return v is not None and v.strip().strip(".").upper().startswith("T")


def read_mesh(kpoints):
    """KPOINTS（自动网格）-> [N1, N2, N3] 或 None。"""
    try:
        ln = Path(kpoints).read_text().splitlines()
        return [int(x) for x in ln[3].split()[:3]]
    except (OSError, IndexError, ValueError):
        return None


def dense_mesh(L0, coarse, dk):
    """逐轴 max(粗网格, ceil(|b_i|/dk))（|b| 含 2π，Å⁻¹；与 S3 / S7 的 DK_MAX 同一口径）。"""
    b = 2 * math.pi * np.linalg.inv(L0).T
    need = [int(math.ceil(np.linalg.norm(b[i]) / dk - 1e-9)) for i in range(3)]
    return [max(int(c), n) for c, n in zip(coarse or [1, 1, 1], need)]


def kpoints_text(mesh, tag):
    return "%s\n0\nGamma\n %d %d %d\n 0 0 0\n" % (tag, mesh[0], mesh[1], mesh[2])


def launcher_line(sub_text):
    """submit.sh 里第一条真正的 VASP 启动命令（跳过注释）。"""
    for ln in sub_text.splitlines():
        s = ln.strip()
        if s and not s.startswith("#") and s.startswith(("mpirun", "srun", "mpiexec")):
            return s
    return None


def read_core_potentials(outcar):
    """OUTCAR 的 "average (electrostatic) potential at core" -> 第一次出现的各原子值（list）；没有 -> []。
    与 pymatgen Outcar.read_avg_core_poten 同一读法：每 17 个字符一项（8 位序号 + 9 位数值），读到 E-fermi 为止。"""
    try:
        lines = Path(outcar).read_text(errors="ignore").splitlines()
    except OSError:
        return []
    for i, ln in enumerate(lines):
        if "the norm of the test charge is" not in ln:
            continue
        out = []
        for l2 in lines[i + 1:]:
            if "E-fermi" in l2:
                break
            if not l2.strip():
                if out:
                    break
                continue
            for j in range(len(l2) // 17):
                try:
                    out.append(float(l2[j * 17 + 8:j * 17 + 17]))
                except ValueError:
                    return out
        return out
    return []


def read_eigenval(path):
    """EIGENVAL -> (kpoints (nk,3), E (nspin, nk, nb), occ (nspin, nk, nb) 或 None, nelect)。"""
    ln = Path(path).read_text().splitlines()
    ispin = int(ln[0].split()[3])
    nelect, nk, nb = (int(float(x)) for x in ln[5].split()[:3])
    E = np.zeros((ispin, nk, nb))
    occ = np.zeros((ispin, nk, nb))
    has_occ = True
    kp = np.zeros((nk, 3))
    i = 6
    for k in range(nk):
        while not ln[i].strip():
            i += 1
        kp[k] = [float(x) for x in ln[i].split()[:3]]
        i += 1
        for b in range(nb):
            t = ln[i].split()
            i += 1
            for s in range(ispin):
                E[s, k, b] = float(t[1 + s])
            if len(t) >= 1 + 2 * ispin:
                for s in range(ispin):
                    occ[s, k, b] = float(t[1 + ispin + s])
            else:
                has_occ = False
    return kp, E, (occ if has_occ else None), nelect


def elapsed(outcar):
    try:
        m = re.search(r"Elapsed time \(sec\):\s*([\d.]+)", Path(outcar).read_text(errors="ignore")[-20000:])
        return float(m.group(1)) if m else None
    except OSError:
        return None


def _key(k, n=10000):
    k = np.asarray(k, float)
    return tuple(np.round((k - np.round(k)) * n).astype(int) % n)


# ---------------------------------------------------------------------------
# make
# ---------------------------------------------------------------------------
def cmd_make(a):
    s7 = Path(a.s7).resolve()
    out = Path(a.out).resolve()
    if s7 == out or s7 in out.parents:
        sys.exit("[ERROR] out_dir 不能在 %s 里面（S7 的扇出会把它当构型、gen 会清理它）" % s7)
    und = s7 / "undeformed"
    if not (und / "POSCAR").is_file():
        sys.exit("[ERROR] %s 没有 undeformed/POSCAR —— S7 还没 gen（或不是 step7_deform 目录）" % s7)
    L0 = read_lattice(und / "POSCAR")
    want = next((P for P, n in PAT_NAMES.items() if n == a.pattern), None)
    if want is None:
        sys.exit("[ERROR] --pattern 只能是 %s" % "/".join(PAT_NAMES.values()))
    found = {}
    for d in sorted(p for p in s7.glob("deform-*") if (p / "POSCAR").is_file()):
        sp = strain_pattern(green_lagrange(L0, read_lattice(d / "POSCAR")))
        if sp:
            found.setdefault(sp[0], []).append((sp[1], d))
    pair = sorted(found.get(want, []))
    if not (len(pair) >= 2 and pair[0][0] < 0 < pair[-1][0]):
        sys.exit("[ERROR] 图样 %s 没有 +δ 与 −δ 两个构型（有的图样：%s）"
                 % (a.pattern, ", ".join("%s%s" % (PAT_NAMES[P], [round(x, 5) for x, _d in v]) for P, v in found.items())
                    or "无"))
    configs = [("undeformed", und, 0.0), (pair[-1][1].name, pair[-1][1], pair[-1][0]),
               (pair[0][1].name, pair[0][1], pair[0][0])]
    for _n, d, _x in configs:
        for f in ("INCAR", "POSCAR", "POTCAR", "KPOINTS"):
            if not (d / f).is_file():
                sys.exit("[ERROR] %s 缺 %s（S7 的 gen 没跑完？）" % (d, f))
    inc0 = (und / "INCAR").read_text()
    if is_true(incar_value(inc0, "LHFCALC")):
        sys.exit("[ERROR] S7 的 INCAR 是杂化泛函（LHFCALC = .TRUE.）：ICHARG = 11 不适用，这个方案不成立")
    if (incar_value(inc0, "ISYM") or "").split()[:1] != ["0"]:
        print("[WARN] S7 的 INCAR 不是 ISYM = 0（%s）：各构型的不可约 k 点可能不一样" % incar_value(inc0, "ISYM"))
    if incar_value(inc0, "ICORELEVEL"):
        print("[WARN] S7 的 INCAR 写了 ICORELEVEL = %s：ICORELEVEL = 1 时 VASP 不写平均静电芯势，AMSET 会改用 1s 芯能级"
              % incar_value(inc0, "ICORELEVEL"))
    coarse = read_mesh(und / "KPOINTS")
    mesh = [int(x) for x in a.mesh] if a.mesh else dense_mesh(L0, coarse, float(a.dk))
    out.mkdir(parents=True, exist_ok=True)
    meta = {"s7": str(s7), "pattern": a.pattern, "coarse_mesh": coarse, "dense_mesh": mesh, "configs": []}
    for name, d, amt in configs:
        cd = out / name
        for st in STAGES:
            sd = cd / st
            if sd.exists():
                shutil.rmtree(sd)
            sd.mkdir(parents=True)
            for f in ("POSCAR", "POTCAR"):
                shutil.copyfile(d / f, sd / f)
        inc = (d / "INCAR").read_text()
        base = inc
        for k, v in (("NSW", "0"), ("IBRION", "-1"), ("LWAVE", ".FALSE.")):
            base = set_tag(base, k, v)
        (cd / "scf_coarse" / "INCAR").write_text(set_tag(base, "LCHARG", ".TRUE."))
        shutil.copyfile(d / "KPOINTS", cd / "scf_coarse" / "KPOINTS")
        nscf = base
        for k, v in (("ICHARG", "11"), ("ISTART", "0"), ("LCHARG", ".FALSE.")):
            nscf = set_tag(nscf, k, v)
        (cd / "nscf_dense" / "INCAR").write_text(nscf)
        (cd / "scf_dense" / "INCAR").write_text(set_tag(base, "LCHARG", ".FALSE."))
        for st in ("nscf_dense", "scf_dense"):
            (cd / st / "KPOINTS").write_text(kpoints_text(mesh, "nscf_dp_check dense %s" % st))
        (cd / "run.sh").write_text(RUN_SH)
        (cd / "run.sh").chmod(0o755)
        sub = d / "submit.sh"
        wrote_sub = False
        if sub.is_file():
            st_ = sub.read_text()
            mpi = launcher_line(st_)
            if mpi:
                head = []
                for ln in st_.splitlines():
                    if ln.strip() == mpi:
                        break
                    ln = re.sub(r"(--job-name[= ])\S+", r"\g<1>nscfdp-" + name, ln)
                    ln = re.sub(r"^(#SBATCH\s+-J\s+)\S+", r"\g<1>nscfdp-" + name, ln)
                    if a.time:
                        ln = re.sub(r"(--time[= ])\S+", r"\g<1>" + a.time, ln)
                        ln = re.sub(r"^(#SBATCH\s+-t\s+)\S+", r"\g<1>" + a.time, ln)
                    head.append(ln)
                body = 'cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1      # 在 <out_dir>/<构型> 里 sbatch\n'
                (cd / "submit.sh").write_text("\n".join(head).rstrip() + "\n" + body + "bash run.sh " + mpi + "\n")
                wrote_sub = True
        meta["configs"].append({"name": name, "src": str(d), "amount": amt, "submit": wrote_sub})
    (out / META).write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    b = 2 * math.pi * np.linalg.inv(L0).T
    print("[OK] %s：图样 %s，构型 %s" % (out, a.pattern, ", ".join("%s(%+.4f)" % (c["name"], c["amount"])
                                                              for c in meta["configs"])))
    print("     粗网格 %s（S7 原样）-> 密网格 %s（%s）；k 点约 %.0f 倍"
          % ("×".join(map(str, coarse or ["?"])), "×".join(map(str, mesh)),
             "--mesh" if a.mesh else "|b_i|/%.3f Å⁻¹，|b| = %s" % (float(a.dk), "/".join("%.3f" % np.linalg.norm(x) for x in b)),
             (np.prod(mesh) / np.prod(coarse)) if coarse else float("nan")))
    for c in meta["configs"]:
        print("     %s：%s" % (c["name"], ("cd %s && sbatch submit.sh" % (out / c["name"])) if c["submit"]
                                 else "提交脚本里写 bash %s/run.sh <VASP 命令>" % (out / c["name"])))
    print("     跑完：python %s analyze %s" % (Path(__file__).name, out))
    return 0


# ---------------------------------------------------------------------------
# analyze
# ---------------------------------------------------------------------------
def _load_stage(sd):
    oc, ev = sd / "OUTCAR", sd / "EIGENVAL"
    done = oc.is_file() and DONE_TAG in oc.read_text(errors="ignore")[-400000:]
    if not done or not ev.is_file():
        return None
    kp, E, occ, nel = read_eigenval(ev)
    core = read_core_potentials(oc)
    return {"k": kp, "E": E, "occ": occ, "nelect": nel, "core": core,
            "ref": (float(np.mean(core)) if core else None), "elapsed": elapsed(oc)}


def _edges(d0, ewin):
    """未形变的密网格 -> 价带顶 / 导带底附近的 (自旋, k, 带) 与带边能量。"""
    E = d0["E"]
    if d0["occ"] is not None:
        occd = d0["occ"] > 0.5 * float(np.max(d0["occ"]))
    else:
        nocc = int(round(d0["nelect"] / 2))          # 没有占据数列（VASP < 5.4）：按非自旋极化、无 SOC 估
        occd = np.zeros(E.shape, bool)
        occd[..., :nocc] = True
    vbm, cbm = float(np.max(E[occd])), float(np.min(E[~occd]))
    vb = [tuple(x) for x in np.argwhere(occd & (E >= vbm - ewin))]
    cb = [tuple(x) for x in np.argwhere(~occd & (E <= cbm + ewin))]
    return vbm, cbm, vb, cb


def _dp(stages, states):
    """{构型: 段数据} -> 每个态的 D（AMSET 写法、±δ 平均）。stages: [(幅度, 数据)]，第一个是未形变。"""
    d0, rest = stages[0][1], stages[1:]
    lut = [{_key(k): j for j, k in enumerate(d["k"])} for _a, d in rest]
    out = []
    for s, k, b in states:
        vals = []
        for (amt, d), lk in zip(rest, lut):
            j = lk.get(_key(d0["k"][k]))
            if j is None:
                return None
            vals.append((d0["E"][s, k, b] - d["E"][s, j, b] - (d0["ref"] - d["ref"])) / amt)
        out.append(float(np.mean(vals)))
    return np.array(out)


def cmd_analyze(a):
    out = Path(a.out).resolve()
    try:
        meta = json.loads((out / META).read_text())
    except (OSError, ValueError):
        print("[ERROR] %s 没有 %s（先 make）" % (out, META))
        return 2
    data, miss = {}, []
    for c in meta["configs"]:
        for st in STAGES:
            r = _load_stage(out / c["name"] / st)
            if r is None:
                miss.append("%s/%s" % (c["name"], st))
            data[(c["name"], st)] = r
    rep = {"pattern": meta["pattern"], "coarse_mesh": meta["coarse_mesh"], "dense_mesh": meta["dense_mesh"]}
    if miss:
        print("[ERROR] 没跑完 / 缺 OUTCAR、EIGENVAL：%s" % ", ".join(miss))
        rep["missing"] = miss
        (out / REPORT).write_text(json.dumps(rep, indent=1, ensure_ascii=False))
        return 2
    names = [c["name"] for c in meta["configs"]]
    amts = {c["name"]: float(c["amount"]) for c in meta["configs"]}
    ok = True
    # 1. 芯能级参考
    print("== 1. 芯能级参考（各原子平均静电势的平均，eV；AMSET 对齐用它）")
    rep["core_ref"] = {}
    for n in names:
        refs = {st: data[(n, st)]["ref"] for st in STAGES}
        nat = {st: len(data[(n, st)]["core"]) for st in STAGES}
        d_ns = (refs["nscf_dense"] - refs["scf_coarse"]) if None not in (refs["nscf_dense"], refs["scf_coarse"]) else None
        d_sd = (refs["scf_dense"] - refs["scf_coarse"]) if None not in (refs["scf_dense"], refs["scf_coarse"]) else None
        rep["core_ref"][n] = {"ref": refs, "n_atoms": nat, "nscf_minus_scf_coarse": d_ns, "scf_dense_minus_scf_coarse": d_sd}
        print("   %-12s scf_coarse %s  nscf_dense %s  scf_dense %s  | 非自洽−粗自洽 %s  密自洽−粗自洽 %s"
              % (n, *("%.5f" % refs[s] if refs[s] is not None else "缺失" for s in STAGES),
                 "%+.2e" % d_ns if d_ns is not None else "—", "%+.2e" % d_sd if d_sd is not None else "—"))
        if refs["nscf_dense"] is None or d_ns is None or abs(d_ns) > a.tol_ref:
            ok = False
    if any(data[(n, "nscf_dense")]["ref"] is None for n in names):
        print("   ★ 非自洽的 OUTCAR 里没有平均静电芯势：AMSET 读不到参考，这个方案要另取参考（如 scf_coarse 的 OUTCAR）")
    # 2. 形变势（密网格：非自洽 vs 自洽；粗网格对照）
    print("== 2. 带边形变势 D（eV，图样 %s，±δ 平均，带符号；AMSET 运行期再取 |·|）" % meta["pattern"])
    dn = {st: [(amts[n], data[(n, st)]) for n in names] for st in STAGES}
    if any(d["ref"] is None for st in ("nscf_dense", "scf_dense") for _a, d in dn[st]):
        print("   有段缺芯能级参考，算不了形变势")
        (out / REPORT).write_text(json.dumps(rep, indent=1, ensure_ascii=False))
        return 2
    vbm, cbm, vb, cb = _edges(dn["scf_dense"][0][1], a.ewin)
    rep["edges"] = {"vbm": vbm, "cbm": cbm, "n_vb_states": len(vb), "n_cb_states": len(cb), "ewin": a.ewin}
    rep["dp"] = {}
    scale = 0.0
    worst = 0.0
    for lab, states in (("价带顶", vb), ("导带底", cb)):
        Dn, Ds = _dp(dn["nscf_dense"], states), _dp(dn["scf_dense"], states)
        if Dn is None or Ds is None:
            print("   %s：各构型的 k 点对不上（ISYM 不是 0？）" % lab)
            return 2
        E0 = dn["scf_dense"][0][1]["E"]
        edge = vbm if lab == "价带顶" else cbm
        at_edge = np.array([abs(E0[s, k, b] - edge) < 1e-3 for s, k, b in states])
        diff = Dn - Ds
        scale = max(scale, float(np.max(np.abs(Ds))))
        worst = max(worst, float(np.max(np.abs(diff))))
        rec = {"n": len(states), "nscf_edge_mean": float(np.mean(Dn[at_edge])), "scf_edge_mean": float(np.mean(Ds[at_edge])),
               "n_edge": int(at_edge.sum()), "max_abs_diff": float(np.max(np.abs(diff))),
               "rms_diff": float(np.sqrt(np.mean(diff ** 2)))}
        # 粗网格对照：带边态最近的粗格点（同一带号）—— S7 现在的做法在谷附近能给到的值
        dc = dn["scf_coarse"]
        kc = dc[0][1]["k"]
        kd0 = dn["scf_dense"][0][1]["k"]
        near = []
        for s, k, b in [x for x, e in zip(states, at_edge) if e]:
            dist = np.linalg.norm((kc - kd0[k]) - np.round(kc - kd0[k]), axis=1)
            near.append((s, int(np.argmin(dist)), b))
        Dc = _dp(dc, near)
        if Dc is not None:
            rec["coarse_nearest_mean"] = float(np.mean(Dc))
        rep["dp"][lab] = rec
        print("   %s：%d 个态（带边 ±%.2f eV），其中正在带边 %d 个" % (lab, len(states), a.ewin, int(at_edge.sum())))
        print("     带边平均 D：密网格自洽 %.4f  密网格非自洽 %.4f%s" % (
            rec["scf_edge_mean"], rec["nscf_edge_mean"],
            ("  粗网格自洽（最近格点）%.4f" % rec["coarse_nearest_mean"]) if Dc is not None else ""))
        print("     非自洽 − 自洽：最大 %.2e，rms %.2e" % (rec["max_abs_diff"], rec["rms_diff"]))
    rep["dp_rel_worst"] = worst / scale if scale > 0 else None
    if scale > 0 and worst > a.tol_d * scale:
        ok = False
    if meta["pattern"] in ("xx", "yy", "zz"):
        print("   注：立方晶系里正应变图样在一整个谷星上的平均 = 静水形变势（Si 电子 Ξ_d + Ξ_u/3、价带顶 a_v）；"
              "上面的「带边平均」就是它（谷星的态都在带边上时）")
    # 3. 耗时
    print("== 3. 墙钟（s，OUTCAR 的 Elapsed time；timing.txt 是作业里量的）")
    rep["elapsed"] = {}
    for n in names:
        t = {st: data[(n, st)]["elapsed"] for st in STAGES}
        rep["elapsed"][n] = t
        tt = {}
        try:
            for ln in (out / n / "timing.txt").read_text().splitlines():
                p = ln.split()
                if len(p) >= 2:
                    tt[p[0]] = float(p[1])
        except OSError:
            pass
        rep["elapsed"][n + "_timing_txt"] = tt
        r1 = (t["scf_coarse"] + t["nscf_dense"]) / t["scf_dense"] if all(t.values()) else None
        print("   %-12s 粗自洽 %s  密非自洽 %s  密自洽 %s  | (粗自洽+密非自洽)/密自洽 = %s"
              % (n, *("%.0f" % t[s] if t[s] else "?" for s in STAGES), "%.2f" % r1 if r1 else "?"))
    rep["pass"] = ok
    (out / REPORT).write_text(json.dumps(rep, indent=1, ensure_ascii=False))
    print("== 结论：%s（芯能级参考差 ≤ %.0e eV；带边 D 的非自洽−自洽 ≤ %.0f%% × max|D|，实测 %s）-> %s"
          % ("通过" if ok else "不通过", a.tol_ref, 100 * a.tol_d,
             "%.2f%%" % (100 * rep["dp_rel_worst"]) if rep["dp_rel_worst"] is not None else "—", out / REPORT))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("s7", help="材料的 step7_deform 目录（只读）")
    m.add_argument("out", help="测试目录（不能在 step7_deform 里）")
    m.add_argument("--pattern", default="xx", help="应变图样 xx/yy/zz/yz/xz/xy（默认 xx）")
    g = m.add_mutually_exclusive_group()
    g.add_argument("--dk", default=0.06, type=float, help="密网格：逐轴 ceil(|b_i|/dk)，Å⁻¹（默认 0.06 = S3 的 DK_MAX_3D）")
    g.add_argument("--mesh", nargs=3, type=int, help="密网格直接给 N1 N2 N3")
    m.add_argument("--time", help="覆盖 submit.sh 的 --time（密网格自洽可能比 S7 原来的墙钟长）")
    z = sub.add_parser("analyze")
    z.add_argument("out")
    z.add_argument("--ewin", type=float, default=0.1, help="带边窗口（eV）")
    z.add_argument("--tol-ref", dest="tol_ref", type=float, default=1e-3, help="芯能级参考 非自洽−粗自洽 的上限（eV）")
    z.add_argument("--tol-d", dest="tol_d", type=float, default=0.02, help="带边 D 的差相对 max|D| 的上限")
    a = ap.parse_args(argv)
    return cmd_make(a) if a.cmd == "make" else cmd_analyze(a)


if __name__ == "__main__":
    sys.exit(main())
