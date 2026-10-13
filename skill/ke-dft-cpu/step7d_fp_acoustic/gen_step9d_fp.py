#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[独立分支 adp-xi-selfcalc] S7.4_fp —— 有限 q 声学冻声子单点（生成 fp-* 扇出）。

实验性步骤，skill.yaml 里 optional_steps.fp_acoustic **默认关**（default: false）：
只有材料 tf_*.yaml 写了 fp_acoustic: true 才会跑。不改 master。

做什么
  1) 取 step7_deform/undeformed 的原胞 + INCAR/KPOINTS/POTCAR；
  2) 沿 a1/a2 各重复 N 次 -> N×N×1 超胞；
  3) 声学冻声子图案 u_k(R)=A*e_k*cos(2*pi*m/N)（LA 沿 a1、TA 面内正交），
     每方向 ±A + 参考，共 5 个 fp-* 子目录；
  4) INCAR 照抄 S7 但**强制 IBRION=-1 / NSW=0**（冻声子 = 原子固定在位移位置，不弛豫）；
     KPOINTS 换成超胞网格（原胞面内点数 / N，kz=1）；
  5) 写 fp_manifest.json（q、振幅、图案、来源 sha256）。

本脚本自足（只依赖 numpy），不 import 技能里的其它模块。算完用
step8.4_amset2d/fp_xi.py extract 提 Xi(q)。
"""
import hashlib
import json
import math
import shutil
import sys
import time
from pathlib import Path

import numpy as np

OUTDIR_NAME = "step7d_fp_acoustic"
DEFORM_DIR = "step7_deform"
STEP_LABEL = "S7.4_fp"
SUB_PREFIX = "fp-"
N_CELL = 4                 # 超胞倍数（N=4 -> q = Gamma->M 的一半；要更小 q 就得更大 N）
DISP_MAX_A = 0.01          # 位移幅值（Å）；反节点应变 = DISP_MAX_A * q_phys
KZ = 1                     # 真空方向 k 点数


# ---------------- 结构读写 / 超胞 / 冻声子图案（与 fp_xi.py 同源） ----------------
def read_poscar(path):
    txt = Path(path).read_text(encoding="utf-8-sig").splitlines()
    comment = txt[0]
    lat = np.array([[float(x) for x in txt[i].split()[:3]] for i in (2, 3, 4)]) * float(txt[1].split()[0])
    i = 5
    toks = txt[i].split()
    if all(t.replace(".", "").replace("-", "").isdigit() for t in toks):
        symbols, counts, i = None, [int(t) for t in toks], i + 1
    else:
        symbols, counts, i = toks, [int(t) for t in txt[i + 1].split()], i + 2
    direct = not txt[i].strip().lower().startswith(("c", "k"))
    i += 1
    n = sum(counts)
    coords = np.array([[float(x) for x in txt[i + j].split()[:3]] for j in range(n)])
    if not direct:
        coords = coords @ np.linalg.inv(lat)
    return comment, lat, symbols, counts, coords


def write_poscar(path, comment, lat, symbols, counts, frac):
    lines = [comment, "    1.000000000000000"]
    for v in lat:
        lines.append("  %20.16f %20.16f %20.16f" % tuple(v))
    lines.append("  " + "  ".join(symbols if symbols else ["X"] * len(counts)))
    lines.append("  " + "  ".join(str(c) for c in counts))
    lines.append("Direct")
    for f in frac:
        lines.append("  %20.16f %20.16f %20.16f" % tuple(f))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def build_supercell(lat, frac, n):
    sup_lat = np.array([lat[0] * n, lat[1] * n, lat[2]])
    fr, cell = [], []
    for m1 in range(n):
        for m2 in range(n):
            for f in frac:
                # ★ 只把 a1/a2 两个方向的分数坐标除以 n；z 方向没扩胞，
                #   f3 必须原样保留（早先写成 f/n 会把 z 压成 1/n，原子挤到 0.78 Å）。
                fr.append(np.array([f[0] / n + m1 / n, f[1] / n + m2 / n, f[2]]))
                cell.append((m1, m2))
    return sup_lat, np.array(fr), cell


def acoustic_pattern(sup_lat, cell, n, branch, amp, sign):
    a1 = sup_lat[0] / n
    n_hat = a1 / np.linalg.norm(a1)
    t_hat = np.cross([0.0, 0.0, 1.0], n_hat)
    t_hat = t_hat / np.linalg.norm(t_hat)
    pol = n_hat if branch.upper() == "LA" else t_hat
    disp = np.zeros((len(cell), 3))
    for j, (m1, _m2) in enumerate(cell):
        disp[j] = amp * sign * pol * math.cos(2.0 * math.pi * m1 / n)
    return disp


def q_phys(lat, n):
    return float(np.linalg.norm(2.0 * math.pi * np.linalg.inv(lat).T[0]) / n)


# ---------------- 生成 ----------------
def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


def _resolve_submit_tpl(step_dir, dim):
    for cand in ("submit_std_%s.tpl" % dim, "submit_std.tpl"):
        f = Path(step_dir) / cand
        if f.is_file():
            return f
    return None


def _read_mesh(kpoints_path):
    try:
        lines = Path(kpoints_path).read_text(encoding="utf-8", errors="ignore").splitlines()
        for ln in lines[3:6]:
            t = ln.split()
            if len(t) >= 3 and all(x.lstrip("-").isdigit() for x in t[:3]):
                return tuple(int(x) for x in t[:3])
    except Exception:                                          # noqa: BLE001
        pass
    return (1, 1, 1)


def _write_kpoints(path, mesh):
    n1, n2 = max(1, int(round(mesh[0] / N_CELL))), max(1, int(round(mesh[1] / N_CELL)))
    Path(path).write_text("Auto\n0\nGamma\n  %d  %d  %d\n0 0 0\n" % (n1, n2, KZ),
                          encoding="utf-8", newline="\n")
    return (n1, n2, KZ)


def _force_frozen(incar_text):
    out = []
    for ln in incar_text.splitlines():
        key = ln.split("=")[0].strip().upper() if "=" in ln else ""
        if key in ("IBRION", "NSW", "ISIF"):
            continue
        out.append(ln)
    out += ["IBRION = -1              # S7.4_fp: 冻声子，原子固定在位移位置",
            "NSW    = 0", "ISIF   = 2"]
    return "\n".join(out) + "\n"


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    und = cwd / DEFORM_DIR / "undeformed"
    for n in ("INCAR", "KPOINTS", "POTCAR", "POSCAR"):
        if not (und / n).is_file():
            sys.exit("[ERROR] 缺 %s/undeformed/%s：先让 S7_deform 生成输入" % (DEFORM_DIR, n))

    comment, lat, symbols, counts, frac = read_poscar(und / "POSCAR")
    sup_lat, sup_frac, cell = build_supercell(lat, frac, N_CELL)
    q = q_phys(lat, N_CELL)
    mesh = _read_mesh(und / "KPOINTS")
    incar_s7 = (und / "INCAR").read_text(encoding="utf-8", errors="ignore")
    mesh_out = _write_kpoints(out / "KPOINTS", mesh)

    configs = [("ref", None, 0), ("lap", "LA", +1), ("lam", "LA", -1),
               ("tap", "TA", +1), ("tam", "TA", -1)]
    stamp = time.strftime("%Y%m%d%H%M%S")
    keep = {SUB_PREFIX + c[0] for c in configs}
    for d in sorted(p for p in out.glob(SUB_PREFIX + "*") if p.is_dir()):
        if d.name not in keep:
            d.rename(out / ("_stale_%s_%s" % (d.name, stamp)))
            print("[..] %s 这一轮不算，改名 _stale_%s_%s" % (d.name, d.name, stamp))

    tp = _resolve_submit_tpl(Path(__file__).resolve().parent, "2d")
    submit_body = tp.read_text(encoding="utf-8") if tp is not None else None
    if submit_body is None:
        print("[WARN] 找不到 submit_std_*.tpl（应随 gen_need 推来）——不写 submit.sh")

    disp_max = 0.0
    for name, branch, sign in configs:
        d = out / (SUB_PREFIX + name)
        d.mkdir(exist_ok=True)
        if branch is None:
            new = sup_frac
        else:
            disp = acoustic_pattern(sup_lat, cell, N_CELL, branch, DISP_MAX_A, sign)
            new = sup_frac + disp @ np.linalg.inv(sup_lat)
            disp_max = max(disp_max, float(np.max(np.linalg.norm(disp, axis=1))))
        write_poscar(d / "POSCAR", "%s %s %s" % (comment.strip(), STEP_LABEL, name),
                     sup_lat, symbols, [c * N_CELL * N_CELL for c in counts], new)
        shutil.copyfile(und / "POTCAR", d / "POTCAR")
        _write_kpoints(d / "KPOINTS", mesh)
        (d / "INCAR").write_text(_force_frozen(incar_s7), encoding="utf-8", newline="\n")
        if submit_body is not None:
            (d / "submit.sh").write_text(
                submit_body.replace("{{JOBNAME}}",
                                    "%s-ke-dft-cpu-%s-%s" % (cwd.name, STEP_LABEL, name)),
                encoding="utf-8", newline="\n")

    manifest = {
        "schema": "fp_acoustic/1", "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_cell": N_CELL, "disp_max_A": DISP_MAX_A, "q_phys_1_A": q,
        "antinode_strain": DISP_MAX_A * q, "n_atoms": len(sup_frac),
        "kpoints_supercell": mesh_out, "kpoints_from_primitive": mesh,
        "configs": [c[0] for c in configs], "branches": ["LA", "TA"],
        "inputs_from": {"dir": "%s/undeformed" % DEFORM_DIR,
                        "sha256": {n: _sha(und / n) for n in ("INCAR", "KPOINTS", "POTCAR", "POSCAR")}},
        "note": "冻声子单点（IBRION=-1,NSW=0）；算完用 fp_xi.py extract 提 Xi(q)",
    }
    (out / "fp_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1),
                                          encoding="utf-8")
    print("[DONE] %s：%d 个单点（ref + LA± + TA±），%d 原子，q=%.5f 1/A，反节点应变 %.5f，"
          "KPOINTS 超胞网格 %s（原胞 %s）"
          % (OUTDIR_NAME, len(configs), len(sup_frac), q, DISP_MAX_A * q, mesh_out, mesh))


if __name__ == "__main__":
    main()
