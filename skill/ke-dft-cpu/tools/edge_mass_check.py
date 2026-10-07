#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""edge_mass_check.py —— 用带边附近的细 k 点（非自洽）直接量 m*，核对 S8.2 在 S3 网格上拟合的 m* 收敛了没有（V176）。

背景：S8.2 的 m* 直接拟合 step3_uniform 的网格点（窗口 4×最近邻）。CrS₂ 从 15×15×1 换成 48×48×3 后
  m* 变了 10–12%；48×48×3 本身够不够，一直没有对照。
做法：在带边 k0 周围，按 ndir 个面内方向、±q（q ≤ qmax）布点做一次非自洽（ICHARG = 11 读 S3 的 CHGCAR）。
  每个方向取 [E(+q) + E(−q)]/2 − E(k0)（消掉三角翘曲这类奇次项），拟合 c₁q² + c₂q⁴，m = ħ²/(2c₁)。
  和 S8.2 的值差 < 3%、且窗口减半后拟合值变化 < 2%，即可认为 S3 网格对 DPT 的 m* 已经收敛。

只支持 2D（面内）。材料目录只读：输入文件**拷贝**到 <out_dir>（不用硬链接；CHGCAR 也不软链，防止 VASP 写回 S3）。

用法：
  python edge_mass_check.py make <材料目录> <out_dir> [--qmax 0.05] [--nq 10] [--ndir 6]
      -> S3 有可用的 CHGCAR：out_dir 里直接是非自洽的 INCAR / KPOINTS / POSCAR / POTCAR / CHGCAR，
         用 step3_uniform 的同一套 VASP 命令提交（k 点约 1 + 2·ndir·nq 个，几分钟）。
      -> S3 的 CHGCAR 缺失或是 0 字节（出厂 LCHARG = .FALSE.）：两段式（V177）。stage1_scf 在 S3 的副本里从 S3 的
         WAVECAR 起跑自洽写出 CHGCAR（几步收敛），stage2_nscf 再非自洽；提交脚本里只写
         bash <out_dir>/run_two_stage.sh <VASP 命令>。
  python edge_mass_check.py fit <out_dir> [--tol 0.03]
      -> 读 EIGENVAL（6 位小数；没有才退回 vasprun.xml 的 4 位），打印各方向的 m* 和与 S8.2 的对比，
         写 edge_mass_result.json。退出码 0 = 通过，1 = 不通过，2 = 算不了。
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE.parent), str(HERE.parent.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402

HB2M = 3.80998            # ħ²/(2 m0)，eV·Å²（与 S8.2 同一常数）
DEGEN_EV = 0.005
TOTEN_TOL_EV = 1e-4       # stage1 自洽必须复现 S3 的总能量（同网格、同 INCAR）
META = "edge_mass_meta.json"
RESULT = "edge_mass_result.json"
# [V177] 两段式（S3 没有可用的 CHGCAR 时）。VASP 命令由调用者给，本工具不碰集群的提交方式。
RUN_TWO_STAGE = r'''#!/bin/bash
# edge_mass_check 两段式（V177）。用法：bash run_two_stage.sh <VASP 命令>，例如 bash run_two_stage.sh mpirun -np 24 vasp_std
set -euo pipefail
[ $# -ge 1 ] || { echo "用法：bash run_two_stage.sh <VASP 命令>" >&2; exit 2; }
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here/stage1_scf"
"$@"
grep -q "General timing" OUTCAR || { echo "[ERROR] stage1 自洽没跑完" >&2; exit 1; }
[ -s CHGCAR ] || { echo "[ERROR] stage1 没写出 CHGCAR" >&2; exit 1; }
cp CHGCAR "$here/stage2_nscf/CHGCAR"
cd "$here/stage2_nscf"
"$@"
grep -q "General timing" OUTCAR || { echo "[ERROR] stage2 非自洽没跑完" >&2; exit 1; }
echo "[OK] 两段都跑完了"
'''


def _nonempty(p):
    try:
        return Path(p).is_file() and Path(p).stat().st_size > 0
    except OSError:
        return False


def outcar_toten(outcar):
    """OUTCAR（或 OUTCAR.gz）里最后一个 free energy TOTEN（eV）；读不到返回 None。"""
    import gzip
    import re
    p = Path(outcar)
    try:
        if p.is_file():
            txt = p.read_text(errors="ignore")
        elif Path(str(p) + ".gz").is_file():
            with gzip.open(str(p) + ".gz", "rt", errors="ignore") as f:
                txt = f.read()
        else:
            return None
    except OSError:
        return None
    m = re.findall(r"free\s+energy\s+TOTEN\s*=\s*(-?\d+\.\d+)", txt)
    return float(m[-1]) if m else None


# ---------------------------------------------------------------- vasprun（只用标准库）
def read_vasprun_eigen(path):
    """-> (kfrac (nk,3), eig (nspin,nk,nb), occ (nspin,nk,nb))。
    k 点取 <kpoints> 下的 kpointlist；本征值取最后一个 <calculation> 直接下面的 <eigenvalues>（不是 <projected> 里那份）。"""
    import xml.etree.ElementTree as ET
    import numpy as np
    kpts, block, stack = None, None, []
    for ev, el in ET.iterparse(str(path), events=("start", "end")):
        if ev == "start":
            stack.append(el.tag)
            continue
        tag = stack.pop()
        parent = stack[-1] if stack else None
        if tag == "varray" and el.get("name") == "kpointlist" and parent == "kpoints" and kpts is None:
            kpts = np.array([[float(x) for x in v.text.split()] for v in el.findall("v")])
        elif tag == "eigenvalues" and parent == "calculation":
            spins = el.find("array").find("set").findall("set")
            block = np.array([[[[float(x) for x in r.text.split()[:2]] for r in ks.findall("r")]
                               for ks in sp.findall("set")] for sp in spins])
        if tag in ("projected", "dos", "eigenvalues", "calculation"):
            el.clear()
    if kpts is None or block is None:
        raise ValueError("%s 里没有 kpointlist 或 eigenvalues（作业没跑完？）" % path)
    if block.shape[1] != len(kpts):
        raise ValueError("%s：k 点 %d 个，本征值 %d 组" % (path, len(kpts), block.shape[1]))
    return kpts, block[..., 0], block[..., 1]


def read_eigenval(path):
    """EIGENVAL -> (kfrac (nk,3), eig (nspin,nk,nb))。EIGENVAL 的本征值有 6 位小数；vasprun.xml 只有 4 位（0.1 meV），
    小 q 点的能量本身只有零点几 meV，拟合要用这一份。VASP 5（无占据数列）/ VASP 6（有）都认。"""
    import numpy as np
    ln = Path(path).read_text(errors="ignore").splitlines()
    ispin = int(ln[0].split()[3])
    nk, nb = (int(x) for x in ln[5].split()[1:3])
    kf, eig, i = [], np.zeros((ispin, nk, nb)), 6
    for k in range(nk):
        while not ln[i].strip():
            i += 1
        kf.append([float(x) for x in ln[i].split()[:3]])
        i += 1
        for b in range(nb):
            t = ln[i].split()
            eig[:, k, b] = [float(x) for x in t[1:1 + ispin]]
            i += 1
    return np.array(kf), eig


def band_edges(eig, occ, occ_tol=0.5):
    """-> {"vbm": (spin, k, band, E), "cbm": (...)}；没有带隙返回 None。"""
    import numpy as np
    occd = occ > occ_tol
    ev = np.where(occd, eig, -np.inf)
    ec = np.where(~occd, eig, np.inf)
    if not np.isfinite(ev).any() or not np.isfinite(ec).any():
        return None
    v = np.unravel_index(int(np.argmax(ev)), ev.shape)
    c = np.unravel_index(int(np.argmin(ec)), ec.shape)
    if eig[c] <= eig[v]:
        return None
    return {"vbm": tuple(int(x) for x in v) + (float(eig[v]),), "cbm": tuple(int(x) for x in c) + (float(eig[c]),)}


# ---------------------------------------------------------------- 布点与拟合
def star_points(rec, k0, qs, ndir):
    """面内 ndir 个方向（0..180° 均分）× ±q -> [(kfrac, j, sign, q)]，第一个是 k0 本身。rec 的行是 b_i（含 2π）。
    θ 从笛卡尔 x（投影到面内）量起，和 S8.2 by_direction 的 x/y 同一套约定。"""
    import numpy as np
    b1, b2 = np.asarray(rec[0], float), np.asarray(rec[1], float)
    n = np.cross(b1, b2)
    n /= np.linalg.norm(n)
    x = np.array([1.0, 0.0, 0.0])
    u1 = x - (x @ n) * n
    if np.linalg.norm(u1) < 1e-6:                                   # 面法向就是 x：退回 b1
        u1 = b1
    u1 = u1 / np.linalg.norm(u1)
    u2 = np.cross(n, u1)
    inv = np.linalg.inv(rec)
    k0 = np.asarray(k0, float)
    k0c = k0 @ rec
    pts = [(k0, -1, 0, 0.0)]
    for j in range(ndir):
        th = np.pi * j / ndir
        d = np.cos(th) * u1 + np.sin(th) * u2
        for s in (1, -1):
            for q in qs:
                pts.append(((k0c + s * q * d) @ inv, j, s, float(q)))
    return pts


def fit_directions(e, pts, ndir, carrier, qfit):
    """e：各点的能量（与 pts 同序，pts[0] = k0）。-> [{"theta_deg", "m", "c1", "c2", "n"} 或带 "error"]。"""
    import numpy as np
    e0 = e[0]
    out = []
    for j in range(ndir):
        plus = {q: e[i] for i, (_, jj, s, q) in enumerate(pts) if jj == j and s == 1}
        minus = {q: e[i] for i, (_, jj, s, q) in enumerate(pts) if jj == j and s == -1}
        qs = sorted(q for q in plus if q in minus and q <= qfit + 1e-12)
        rec = {"theta_deg": round(180.0 * j / ndir, 3), "n": len(qs)}
        if len(qs) < 3:
            rec["error"] = "窗口内点太少"
            out.append(rec)
            continue
        q = np.array(qs)
        eb = np.array([(plus[x] + minus[x]) / 2.0 - e0 for x in qs])
        c1, c2 = np.linalg.lstsq(np.c_[q ** 2, q ** 4], eb, rcond=None)[0]
        rec.update(c1=float(c1), c2=float(c2))
        if (carrier == "electron" and c1 <= 0) or (carrier == "hole" and c1 >= 0):
            rec["error"] = "曲率符号不对：k0 不是这条带的极值点"
        else:
            rec["m"] = float(HB2M / abs(c1))
        out.append(rec)
    return out


def summarize(dirs):
    import numpy as np
    ms = sorted(d["m"] for d in dirs if "m" in d)
    if len(ms) < 2:
        return None
    return {"m_light2": float(np.sqrt(ms[0] * ms[1])),               # S8.2 的口径：最轻两支的几何均值
            "m_geo": float(np.exp(np.mean(np.log(ms)))), "m_min": ms[0], "m_max": ms[-1],
            "anisotropy": ms[-1] / ms[0]}


# ---------------------------------------------------------------- make
def _dpt_masses(mat):
    p = mat / "step8.2_dpt" / "dpt_result.json"
    try:
        res = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {r["carrier"]: r["inputs"]["m_eff_m0"] for r in res.get("results", [])
            if r.get("inputs", {}).get("m_eff_m0")}


def cmd_make(a):
    import numpy as np
    mat = Path(a.material).resolve()
    s3 = mat / "step3_uniform"
    out = Path(a.out_dir).resolve()
    if out == mat or mat in out.parents:
        sys.exit("[ERROR] out_dir 不能放在材料目录里（%s）：这是一次性核对，放到项目外的临时目录" % out)
    if out.exists() and any(out.iterdir()):
        sys.exit("[ERROR] %s 已存在且不空" % out)
    for f in ("INCAR", "POSCAR", "POTCAR", "KPOINTS", "vasprun.xml"):
        if not _nonempty(s3 / f):
            sys.exit("[ERROR] 缺 %s 或是空文件（S3 没跑完？）" % (s3 / f))
    # [V177] S3 的出厂 INCAR 是 LCHARG = .FALSE.，CHGCAR 只是个 0 字节的占位文件（CrS₂ / CrSe₂ 实测）。
    #   V176 只查 is_file()，0 字节也放行，VASP 要到非自洽那步才会挂。没有可用的 CHGCAR 就走两段式：
    #   stage1 在 S3 的副本里从 S3 的 WAVECAR 起跑自洽（ISTART = 1、ICHARG = 0，几步就收敛）写出 CHGCAR，stage2 再非自洽。
    two_stage = not _nonempty(s3 / "CHGCAR")
    has_wave = _nonempty(s3 / "WAVECAR")
    dim = kc.read_method_dim(s3 / kc.METHOD_FILE) or kc.resolve_dim_for(s3 / "POSCAR", "auto")[0]
    if dim != "2d":
        sys.exit("[ERROR] 只支持 2D（本材料 DIM=%s）" % dim)
    kfrac, eig, occ = read_vasprun_eigen(s3 / "vasprun.xml")
    edges = band_edges(eig, occ)
    if edges is None:
        sys.exit("[ERROR] S3 的 vasprun 里找不到带隙（金属？）")
    rec = 2.0 * np.pi * np.linalg.inv(kc.read_lattice_matrix(s3 / "POSCAR")).T
    qs = [a.qmax * (i + 1) / a.nq for i in range(a.nq)]
    stars, kp_lines, pts_meta = [], [], []
    for name in ("cbm", "vbm"):
        sp, k, b, e = edges[name]
        k0 = [float(x) for x in kfrac[k]]
        for st in stars:
            if np.allclose(st["k0"], k0, atol=1e-6):
                st["edges"].append(name)
                break
        else:
            pts = star_points(rec, k0, qs, a.ndir)
            stars.append({"k0": k0, "edges": [name], "start": len(pts_meta), "count": len(pts)})
            for kf, j, s, q in pts:
                kp_lines.append("  %.10f  %.10f  %.10f  1" % tuple(kf))
                pts_meta.append([len(stars) - 1, j, s, q])
    out.mkdir(parents=True, exist_ok=True)
    base = kc.parse_incar((s3 / "INCAR").read_text(errors="ignore"))
    nb = str(eig.shape[2])
    _quiet = {"LVHAR": None, "LVTOT": None, "LELF": None, "LAECHG": None, "LORBIT": None, "ICORELEVEL": None}
    nscf = out / "stage2_nscf" if two_stage else out
    nscf.mkdir(exist_ok=True)
    if two_stage:
        scf = out / "stage1_scf"
        scf.mkdir()
        for f in ("POSCAR", "POTCAR", "KPOINTS") + (("WAVECAR",) if has_wave else ()):
            shutil.copyfile(s3 / f, scf / f)                         # 拷贝，不链接（VASP 不许写回 S3）
        ov1 = dict(_quiet, ISTART="1" if has_wave else "0", ICHARG="0" if has_wave else "2", LCHARG=".TRUE.",
                   LWAVE=".FALSE.", NBANDS=nb, NSW="0", IBRION="-1")
        (scf / "INCAR").write_text(kc.incar_text(kc.merge_incar(base, ov1), system="edge_mass_check stage1"))
        (out / "run_two_stage.sh").write_text(RUN_TWO_STAGE)
        (out / "run_two_stage.sh").chmod(0o755)
        if not has_wave:
            print("[WARN] S3 也没有 WAVECAR：stage1 只能从头自洽，成本约等于重跑一次 S3")
    for f in ("POSCAR", "POTCAR") + (() if two_stage else ("CHGCAR",)):
        shutil.copyfile(s3 / f, nscf / f)                           # 拷贝，不链接
    ov = dict(_quiet, ICHARG="11", ISTART="0", ISYM="0", NSW="0", IBRION="-1", NBANDS=nb,
              LWAVE=".FALSE.", LCHARG=".FALSE.", KSPACING=None, KGAMMA=None)
    try:
        if float(base.get("EDIFF", "1E-4").lower().replace("d", "e")) > 1e-8:
            ov["EDIFF"] = "1E-8"                                     # 最小的 ΔE 只有 0.1 meV 量级
    except ValueError:
        ov["EDIFF"] = "1E-8"
    if base.get("ALGO", "N").strip()[:1].upper() in ("F", "V"):
        ov["ALGO"] = "Normal"
    (nscf / "INCAR").write_text(kc.incar_text(kc.merge_incar(base, ov), system="edge_mass_check"))
    (nscf / "KPOINTS").write_text("edge_mass_check star (V176)\n%d\nReciprocal\n%s\n" % (len(kp_lines), "\n".join(kp_lines)))
    meta = {"version": 2, "material": str(mat), "s3_mesh": kc.read_kpoints_mesh(s3 / "KPOINTS"),
            "layout": "two_stage" if two_stage else "single", "s3_toten": outcar_toten(s3 / "OUTCAR"),
            "qmax": a.qmax, "nq": a.nq, "ndir": a.ndir,
            "edges": {n: {"spin": edges[n][0], "band": edges[n][2], "E_s3": edges[n][3]} for n in edges},
            "stars": stars, "points": pts_meta, "dpt_m": _dpt_masses(mat)}
    (out / META).write_text(json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
    print("[OK] %s：%d 个 k 点（%d 个星，%s），NBANDS = %s" % (
        out, len(kp_lines), len(stars), "；".join("k0=%s 用于 %s" % (st["k0"], "/".join(st["edges"])) for st in stars), nb))
    if two_stage:
        print("     S3 没有可用的 CHGCAR（LCHARG = .FALSE.）-> 两段式：stage1_scf 从 S3 的 %s 起跑自洽写出 CHGCAR，"
              "stage2_nscf 再非自洽。" % ("WAVECAR" if has_wave else "初始猜测"))
        print("     提交脚本里只写一行：bash %s/run_two_stage.sh <S3 用的 VASP 命令，如 mpirun -np 24 vasp_std>" % out)
    else:
        print("     用 step3_uniform 的同一套 VASP 命令在 %s 里提交" % out)
    print("     跑完：python %s fit %s" % (Path(__file__).name, out))
    return 0


# ---------------------------------------------------------------- fit
def cmd_fit(a):
    import numpy as np
    out = Path(a.out_dir).resolve()
    try:
        meta = json.loads((out / META).read_text(encoding="utf-8"))
        nscf = out / "stage2_nscf" if meta.get("layout") == "two_stage" else out
        if (nscf / "EIGENVAL").is_file():
            kfrac, eig = read_eigenval(nscf / "EIGENVAL")
        else:
            print("[WARN] 没有 EIGENVAL，退回 vasprun.xml：本征值只有 4 位小数（0.1 meV），窗口减半的检查会偏严")
            kfrac, eig, _ = read_vasprun_eigen(nscf / "vasprun.xml")
    except (OSError, ValueError, IndexError) as e:
        print("[ERROR] %s" % e)
        return 2
    if len(kfrac) != len(meta["points"]):
        print("[ERROR] 结果有 %d 个 k 点，meta 记了 %d 个 —— 不是这次 make 的输入" % (len(kfrac), len(meta["points"])))
        return 2
    ok, res = True, {"material": meta["material"], "s3_mesh": meta["s3_mesh"], "qmax": meta["qmax"], "carriers": {}}
    if meta.get("layout") == "two_stage":                            # [V177] stage1 的电荷密度得是 S3 那一份
        e1, e3 = outcar_toten(out / "stage1_scf" / "OUTCAR"), meta.get("s3_toten")
        res["stage1_toten_minus_s3_eV"] = None if None in (e1, e3) else e1 - e3
        if None in (e1, e3):
            print("[WARN] 读不到 stage1 或 S3 的 TOTEN，没法确认 stage1 复现了 S3 的电荷密度")
        elif abs(e1 - e3) > TOTEN_TOL_EV:
            print("[★] stage1 的总能量和 S3 差 %.2e eV（> %.0e）：电荷密度不是 S3 那一份，m* 不可信"
                  % (e1 - e3, TOTEN_TOL_EV))
            ok = False
        else:
            print("[OK] stage1 复现了 S3 的总能量（差 %.1e eV）" % (e1 - e3))
    for name, carrier in (("cbm", "electron"), ("vbm", "hole")):
        ed = meta["edges"][name]
        si = next(i for i, st in enumerate(meta["stars"]) if name in st["edges"])
        st = meta["stars"][si]
        sl = slice(st["start"], st["start"] + st["count"])
        pts = [(None, j, s, q) for _, j, s, q in meta["points"][sl]]
        band = eig[ed["spin"], sl, ed["band"]]
        e0 = band[0]
        nb = eig.shape[2]
        degen = [b for b in (ed["band"] - 1, ed["band"] + 1)
                 if 0 <= b < nb and abs(eig[ed["spin"], st["start"], b] - e0) < DEGEN_EV]
        full = fit_directions(band, pts, meta["ndir"], carrier, meta["qmax"])
        half = fit_directions(band, pts, meta["ndir"], carrier, meta["qmax"] / 2.0)
        sf, sh = summarize(full), summarize(half)
        r = {"directions": full, "summary": sf, "summary_half_window": sh, "degenerate_with": degen,
             "dpt_m": meta["dpt_m"].get(carrier)}
        print("== %s（band %d, spin %d, k0 = %s）" % (carrier, ed["band"], ed["spin"], st["k0"]))
        for d in full:
            print("   θ = %6.1f°  %s" % (d["theta_deg"], ("m = %.4f" % d["m"]) if "m" in d else d["error"]))
        if degen:
            print("   [WARN] 带边和 band %s 简并（< %.0f meV）：只拟合了 band %d" % (degen, DEGEN_EV * 1000, ed["band"]))
        if sf is None or sh is None:
            print("   [ERROR] 有效方向不足 2 个，拟合不了")
            ok = False
            res["carriers"][carrier] = r
            continue
        r["half_window_change"] = sh["m_light2"] / sf["m_light2"] - 1.0
        print("   m*（最轻两支几何均值，S8.2 口径）= %.4f；全部方向几何均值 %.4f；各向异性 %.3f；窗口减半 %+.1f%%"
              % (sf["m_light2"], sf["m_geo"], sf["anisotropy"], 100 * r["half_window_change"]))
        if abs(r["half_window_change"]) > 0.02:
            print("   [WARN] 窗口减半后变化 > 2%%：小 q 点的能量精度不够或非抛物太强，结论要打折扣")
            ok = False
        if r["dpt_m"]:
            r["dpt_over_fine"] = r["dpt_m"] / sf["m_light2"]
            good = abs(r["dpt_over_fine"] - 1.0) < a.tol
            print("   S8.2（S3 网格 %s）m* = %.4f，是细网格值的 %.3f 倍 -> %s"
                  % ("×".join(map(str, meta["s3_mesh"] or [])), r["dpt_m"], r["dpt_over_fine"],
                     "[OK] 已收敛" if good else "[★] 差 ≥ %.0f%%：S3 网格对 m* 还不够" % (100 * a.tol)))
            ok = ok and good
        else:
            print("   （没有 step8.2_dpt/dpt_result.json，只给细网格值）")
        res["carriers"][carrier] = r
    res["pass"] = ok
    (out / RESULT).write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    print("[%s] 结果写到 %s" % ("OK" if ok else "★", out / RESULT))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("material")
    m.add_argument("out_dir")
    m.add_argument("--qmax", type=float, default=0.05, help="Å⁻¹（含 2π），默认 0.05 = DK_MAX_2D")
    m.add_argument("--nq", type=int, default=10)
    m.add_argument("--ndir", type=int, default=6)
    f = sub.add_parser("fit")
    f.add_argument("out_dir")
    f.add_argument("--tol", type=float, default=0.03)
    a = ap.parse_args(argv)
    return cmd_make(a) if a.cmd == "make" else cmd_fit(a)


if __name__ == "__main__":
    sys.exit(main())
