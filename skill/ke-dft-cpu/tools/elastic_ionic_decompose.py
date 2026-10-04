#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""elastic_ionic_decompose.py —— 弹性张量的离子弛豫贡献按本征模拆开（V154；只读，秒级，不需要 VASP）。

**S6 的 TOTAL ELASTIC MODULI 不正定时先跑本工具，再决定要不要重算。**
起因（Mo2S3，S6 = IBRION=6 + ISIF=3，20 原子）：面内 C66 的刚性离子项 +267.80 kBar，离子弛豫项 ≈ −430.7，
TOTAL = −162.94（不正定，V134 拦住 S8）；Γ 点只有 3 个近零虚频（0.006–0.015 THz，即平移模）。
VASP 的离子弛豫项 = −(1/Ω)·Λᵀ K⁻¹ Λ（Λ = 内应变张量，K = 力常数矩阵）。均匀平移和应变不耦合（Λ 在平移方向上
应为 0），但数值上 Λ 有残差、K 在平移方向上的本征值≈0，两个小量相除可以给出任意大的贡献。

做法：读 S6 目录的 OUTCAR —— SECOND DERIVATIVES、INTERNAL STRAIN TENSOR、SYMMETRIZED（刚性离子）/
CONTR FROM IONIC RELAXATION / TOTAL 三块弹性模量、体积：
  1) 把 K 的每个本征模的贡献 −(1/Ω)(uᵀΛ)ᵀ(uᵀΛ)/λ 加起来，**先复现 VASP 自己的离子弛豫项**（自检：对不上 = 本工具
     对 OUTCAR 的约定理解有误，结论不可用，退出码 2）；剪切列的因子（1 / ½ / 2）按最佳复现自动定；
  2) 按与均匀平移的重叠（> 0.9）认出 3 个平移模，单列它们的贡献；
  3) 扣掉平移模 = "投影后的弛豫张量"，判正定（2D 看面内 XX/YY/XY，3D 看 6×6）；
  4) 结论：投影后正定且翻负主要来自平移模 -> 数值假象（不用重算 S6）；否则列出贡献最大的非平移模 -> 真实的内应变
     耦合（Born 不稳定 / 软模），不可用，重算也不会变。
产出：<S6 目录>/elastic_ionic_decompose.json（含投影后的 TOTAL，VASP 顺序 XX YY ZZ XY YZ ZX，kBar）。

用法：
    python <skill>/ke-dft-cpu/tools/elastic_ionic_decompose.py [<材料目录或 step6_elastic 目录>] [--dim 2d|3d]
退出码：0 结论可用；1 投影后仍不正定（真实不稳定）；2 输入问题或复现不了 VASP 的离子项。
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

EV_A3_TO_KBAR = 1602.1766208          # eV/Å³ -> kBar
VASP_ORDER = ("XX", "YY", "ZZ", "XY", "YZ", "ZX")
INPLANE = (0, 1, 3)                   # VASP 顺序里的 XX YY XY
MATCH_TOL = 0.05                      # 复现 VASP 离子项的相对误差上限
TRANS_OVERLAP = 0.9
_LABEL = re.compile(r"^(\d+)([XYZxyz])$")
_FLOAT = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$")


def _is_float(t):
    return bool(_FLOAT.match(t))


def _modulus_block(txt, title):
    """最后一块 <title> 弹性模量 -> 6×6（VASP 行列顺序 XX YY ZZ XY YZ ZX，kBar）；没有 -> None。"""
    i = txt.rfind(title)
    if i < 0:
        return None
    rows = {}
    for ln in txt[i:].splitlines()[1:16]:
        p = ln.split()
        if p and p[0].upper() in VASP_ORDER and len(p) >= 7 and all(_is_float(x) for x in p[1:7]):
            rows[p[0].upper()] = [float(x) for x in p[1:7]]
        if len(rows) == 6:
            break
    if len(rows) != 6:
        return None
    return np.array([rows[k] for k in VASP_ORDER])


def _second_derivatives(txt, n):
    i = txt.rfind("SECOND DERIVATIVES")
    if i < 0:
        return None, None
    labels, vals = [], []
    cur = None
    for ln in txt[i:].splitlines()[1:]:
        toks = ln.split()
        if not toks or set(ln.strip()) <= {"-"}:
            if cur is not None and len(labels) == 3 * n and len(vals[-1]) == 3 * n:
                break
            continue
        if all(_LABEL.match(t) for t in toks):                 # 表头
            continue
        if _LABEL.match(toks[0]):
            cur = toks[0].upper()
            labels.append(cur)
            vals.append([float(x) for x in toks[1:] if _is_float(x)])
        elif cur is not None and all(_is_float(t) for t in toks):   # 折行
            vals[-1].extend(float(x) for x in toks)
        else:
            if len(labels) == 3 * n:
                break
        if len(labels) == 3 * n and len(vals[-1]) == 3 * n:
            break
    if len(labels) != 3 * n or any(len(v) != 3 * n for v in vals):
        return None, None
    return np.array(vals), labels


def _internal_strain(txt, n):
    lam = np.zeros((3 * n, 6))
    seen = set()
    for m in re.finditer(r"INTERNAL STRAIN TENSOR FOR ION\s+(\d+)", txt):
        ion = int(m.group(1)) - 1
        if ion >= n:
            continue
        got = {}
        for ln in txt[m.end():].splitlines()[1:8]:
            p = ln.split()
            if p and p[0].lower() in ("x", "y", "z") and len(p) >= 7 and all(_is_float(x) for x in p[1:7]):
                got[p[0].lower()] = [float(x) for x in p[1:7]]
            if len(got) == 3:
                break
        if len(got) == 3:
            for k, ax in enumerate("xyz"):
                lam[3 * ion + k] = got[ax]
            seen.add(ion)
    return lam if len(seen) == n else None


def parse_outcar(path):
    txt = Path(path).read_text(errors="ignore")
    mn = re.search(r"NIONS\s*=\s*(\d+)", txt)
    vols = re.findall(r"volume of cell\s*:\s*([\d.]+)", txt)
    if not mn or not vols:
        raise ValueError("OUTCAR 里没有 NIONS / volume of cell")
    n = int(mn.group(1))
    s, labels = _second_derivatives(txt, n)
    if s is None:
        raise ValueError("读不到 SECOND DERIVATIVES（%d×%d）—— 不是 IBRION=5/6 的 OUTCAR？" % (3 * n, 3 * n))
    want = ["%d%s" % (a + 1, ax) for a in range(n) for ax in "XYZ"]
    if labels != want:
        raise ValueError("SECOND DERIVATIVES 的自由度顺序不是 1X 1Y 1Z 2X …：%s…" % labels[:6])
    lam = _internal_strain(txt, n)
    if lam is None:
        raise ValueError("读不到每个离子的 INTERNAL STRAIN TENSOR —— 不是 ISIF≥3 的 IBRION=6？")
    blocks = {k: _modulus_block(txt, t) for k, t in (
        ("clamped", "SYMMETRIZED ELASTIC MODULI"), ("ionic", "ELASTIC MODULI CONTR FROM IONIC RELAXATION"),
        ("total", "TOTAL ELASTIC MODULI"))}
    miss = [k for k, v in blocks.items() if v is None]
    if miss:
        raise ValueError("读不到弹性模量块：%s" % ", ".join(miss))
    return {"n": n, "volume": float(vols[-1]), "second": s, "lam": lam, **blocks}


def _is_pd(m):
    m = 0.5 * (np.asarray(m) + np.asarray(m).T)
    return bool(np.all(np.linalg.eigvalsh(m) > 0))


def decompose(d, dim="2d"):
    n, vol = d["n"], d["volume"]
    s = 0.5 * (d["second"] + d["second"].T)
    k = s if np.trace(s) > 0 else -s                        # 力常数矩阵对角为正；VASP 打印的符号两种都接受
    w, u = np.linalg.eigh(k)
    trans = np.zeros((3, 3 * n))
    for ax in range(3):
        trans[ax, ax::3] = 1.0 / np.sqrt(n)
    overlap = np.sum((trans @ u) ** 2, axis=0)              # 每个模与均匀平移的重叠
    is_t = overlap > TRANS_OVERLAP
    scale = max(np.max(np.abs(d["ionic"])), 1e-9)
    best = None
    for f in (1.0, 0.5, 2.0):                               # 剪切列（XY YZ ZX）的约定因子
        lam = d["lam"].copy()
        lam[:, 3:] *= f
        p = u.T @ lam                                       # (3n, 6)
        contrib = np.array([-np.outer(p[m], p[m]) / w[m] for m in range(3 * n)]) * EV_A3_TO_KBAR / vol
        # VASP 求逆时可能含平移模（数值上没去掉），也可能已经去掉：两种都试，哪种复现得上就是哪种
        for incl in (True, False):
            full = contrib.sum(axis=0) if incl else contrib[~is_t].sum(axis=0)
            err = float(np.max(np.abs(full - d["ionic"])) / scale)
            if best is None or err < best[0]:
                best = (err, f, contrib, incl)
    err, f, contrib, incl = best
    proj_ionic = contrib[~is_t].sum(axis=0)
    trans_ionic = contrib[is_t].sum(axis=0)
    proj_total = d["clamped"] + proj_ionic
    idx = list(INPLANE) if dim == "2d" else list(range(6))
    sub = lambda m: np.asarray(m)[np.ix_(idx, idx)]         # noqa: E731
    total_pd, proj_pd = _is_pd(sub(d["total"])), _is_pd(sub(proj_total))
    # 翻负的那个分量：TOTAL 在所查块里最负的对角元
    diag = [(i, d["total"][i, i]) for i in idx]
    worst = min(diag, key=lambda x: x[1])[0]
    nt = [m for m in range(3 * n) if not is_t[m]]
    top = sorted(nt, key=lambda m: contrib[m][worst, worst])[:3]
    res = {
        "dim": dim, "n_ions": n, "volume_A3": vol, "shear_factor": f,
        "reproduce_rel_err": round(err, 4), "reproduced": err <= MATCH_TOL,
        "vasp_includes_translations": bool(incl),
        "n_translation_modes": int(is_t.sum()),
        "translation_eigenvalues_eV_A2": [round(float(w[m]), 6) for m in range(3 * n) if is_t[m]],
        "worst_component": VASP_ORDER[worst],
        "clamped": round(float(d["clamped"][worst, worst]), 2),
        "vasp_ionic": round(float(d["ionic"][worst, worst]), 2),
        "translation_part": round(float(trans_ionic[worst, worst]), 2),
        "other_modes_part": round(float(proj_ionic[worst, worst]), 2),
        "top_other_modes": [{"mode": int(m), "eigenvalue_eV_A2": round(float(w[m]), 5),
                             "contrib_kbar": round(float(contrib[m][worst, worst]), 2)} for m in top],
        "total_pd": total_pd, "projected_pd": proj_pd,
        "projected_total_vasp_order_kbar": np.round(proj_total, 3).tolist(),
    }
    vi = d["ionic"][worst, worst]
    t_share = float(trans_ionic[worst, worst] / vi) if abs(vi) > 1e-9 else 0.0
    res["translation_share"] = round(t_share, 3)
    if not res["reproduced"]:
        res["verdict"] = "unreliable"
    elif incl and proj_pd and not total_pd and t_share > 0.5:
        res["verdict"] = "artifact"
    elif proj_pd:
        res["verdict"] = "ok"
    else:
        res["verdict"] = "unstable"
    return res


def _dim_of(mat):
    for d in ("step1_opt", "step1_std_opt"):
        f = mat / d / "workflow_method.txt"
        if f.is_file():
            for ln in f.read_text(errors="ignore").splitlines():
                if ln.upper().startswith("DIM="):
                    return ln.split("=", 1)[1].strip().lower()
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("path", nargs="?", default=".", help="材料目录或 step6_elastic 目录（默认当前目录）")
    ap.add_argument("--dim", choices=("2d", "3d"), default=None)
    a = ap.parse_args(argv)
    p = Path(a.path)
    s6 = p if (p / "OUTCAR").is_file() and not (p / "step6_elastic").is_dir() else p / "step6_elastic"
    oc = s6 / "OUTCAR"
    if not oc.is_file():
        print("[ERROR] 找不到 %s" % oc, file=sys.stderr)
        return 2
    dim = a.dim or _dim_of(s6.parent) or "3d"
    try:
        r = decompose(parse_outcar(oc), dim)
    except (ValueError, np.linalg.LinAlgError) as e:
        print("[ERROR] %s: %s" % (type(e).__name__, e), file=sys.stderr)
        return 2
    (s6 / "elastic_ionic_decompose.json").write_text(json.dumps(r, ensure_ascii=False, indent=1),
                                                    encoding="utf-8")
    c = r["worst_component"]
    print("== %s（%s，%d 原子；%s 块）" % (oc, dim.upper(), r["n_ions"], "面内 XX/YY/XY" if dim == "2d" else "6×6"))
    print("   复现 VASP 离子弛豫项：相对误差 %.1f%%（剪切列因子 %g；VASP 求逆时%s平移模）%s"
          % (100 * r["reproduce_rel_err"], r["shear_factor"],
             "含" if r["vasp_includes_translations"] else "已去掉", "" if r["reproduced"] else "  ★ 复现不了"))
    print("   %s：刚性离子 %+.2f，离子弛豫 %+.2f = 平移模 %+.2f（%d 个，本征值 %s eV/Å²）+ 其它模 %+.2f"
          % (c, r["clamped"], r["vasp_ionic"], r["translation_part"], r["n_translation_modes"],
             r["translation_eigenvalues_eV_A2"], r["other_modes_part"]))
    for t in r["top_other_modes"]:
        print("     其它模 #%d：本征值 %.4f eV/Å²，贡献 %+.2f kBar" % (t["mode"], t["eigenvalue_eV_A2"], t["contrib_kbar"]))
    print("   扣掉平移模后：%s = %+.2f kBar；%s %s"
          % (c, r["clamped"] + r["other_modes_part"], "面内 3×3" if dim == "2d" else "6×6",
             "正定" if r["projected_pd"] else "仍不正定"))
    v = r["verdict"]
    if v == "unreliable":
        print("[ERROR] ★ 复现不了 VASP 自己的离子弛豫项：本工具对 OUTCAR 约定的理解有误，以上结论不可用")
        return 2
    if v == "artifact":
        print("[OK] 结论：数值假象 —— 翻负的 %.0f%% 来自平移模（均匀平移和应变不耦合，这部分本该是 0）。"
              "扣掉后的弛豫张量正定，不用重算 S6；投影后的 TOTAL 写在 elastic_ionic_decompose.json"
              % (100 * r["translation_share"]))
        return 0
    if v == "unstable":
        print("[WARN] ★ 结论：扣掉平移模仍不正定 —— 真实的内应变耦合（贡献最大的非平移模见上，Born 不稳定或软模），"
              "不可用；重算 S6 不会改变")
        return 1
    print("[OK] 结论：TOTAL 本来就正定（或扣平移模前后都正定）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
