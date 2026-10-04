#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""elastic_ionic_decompose.py —— 弹性张量的离子弛豫贡献按本征模拆开（V154；只读，秒级，不需要 VASP）。

**S6 的 TOTAL ELASTIC MODULI 不正定时先跑本工具，再决定要不要重算。**
起因（Mo2S3，S6 = IBRION=6 + ISIF=3，20 原子）：面内 C66 的刚性离子项 +267.80 kBar，离子弛豫项 ≈ −430.7，
TOTAL = −162.94（不正定，V134 拦住 S8）；Γ 点只有 3 个近零虚频（0.006–0.015 THz，即平移模）。
VASP 的离子弛豫项 = −(1/Ω)·Λᵀ K⁻¹ Λ（Λ = 内应变张量，K = 力常数矩阵）。均匀平移和应变不耦合（Λ 在平移方向上
应为 0），但数值上 Λ 有残差、K 在平移方向上的本征值≈0，两个小量相除可以给出任意大的贡献。

做法：读 S6 目录的 OUTCAR 与 vasprun.xml —— 力常数（OUTCAR 的 SECOND DERIVATIVES；VASP 6 的 IBRION=6 不打印它，
就用 vasprun.xml 的 <dynmat> hessian，去质量加权与否都试）、INTERNAL STRAIN TENSOR（VASP 6 分 FROM STRAINED CELLS /
FROM DISPLACED ATOMS 两套，另试平均）、SYMMETRIZED（刚性离子）与 TOTAL 弹性模量（离子项没单独打印就用两者之差）、体积：
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
    """-> [(口径, Λ 3n×6), ...]。VASP 6 分两段打印：FROM STRAINED CELLS（应变后的力）与 FROM DISPLACED ATOMS
    （位移后的应力）；两段都在时另给平均。没有分段标题（旧版）-> 只有一套。"""
    heads = [(m.start(), lab) for lab, pat in (("strained cells", r"INTERNAL STRAIN TENSORS FROM STRAINED CELLS"),
                                              ("displaced atoms", r"INTERNAL STRAIN TENSORS FROM DISPLACED ATOMS"))
             for m in re.finditer(pat, txt)]
    if not heads:
        lam = _strain_blocks(txt, n, 0, len(txt))
        return [("single", lam)] if lam is not None else []
    heads.sort()
    out = []
    for i, (pos, lab) in enumerate(heads):
        end = heads[i + 1][0] if i + 1 < len(heads) else len(txt)
        lam = _strain_blocks(txt, n, pos, end)
        if lam is not None:
            out.append((lab, lam))
    got = dict(out)
    if "strained cells" in got and "displaced atoms" in got:
        out.append(("average", 0.5 * (got["strained cells"] + got["displaced atoms"])))
    return out


def _strain_blocks(txt, n, start, end):
    lam = np.zeros((3 * n, 6))
    seen = set()
    for m in re.finditer(r"INTERNAL STRAIN TENSOR FOR ION\s+(\d+)", txt[start:end]):
        ion = int(m.group(1)) - 1
        if ion >= n:
            continue
        m_end = start + m.end()
        got = {}
        for ln in txt[m_end:].splitlines()[1:8]:
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


def _vasprun_hessian(path, n):
    """vasprun.xml 的 <dynmat><varray name="hessian">（3n×3n）与每个原子的质量；没有 -> (None, None)。
    VASP 的 IBRION=5–8 都写这一块（phonopy/pymatgen 也从这里读力常数）；OUTCAR 不一定打印 SECOND DERIVATIVES
    （Mo2S3：IBRION=6 + ISIF=3 的 OUTCAR 里就没有）。"""
    import gzip
    import xml.etree.ElementTree as ET
    hess, types, kinds = None, [], []
    try:
        fh = gzip.open(str(path), "rb") if str(path).endswith(".gz") else open(str(path), "rb")
        with fh:
            for el in _iter_end(ET, fh):                    # 逐个处理：文件截断时已读到的部分照样可用
                if el.tag == "array" and el.get("name") == "atomtypes":
                    for rc in el.iter("rc"):
                        c = [x.text.strip() if x.text else "" for x in rc.findall("c")]
                        if len(c) >= 3:
                            types.append(float(c[2]))
                elif el.tag == "array" and el.get("name") == "atoms":
                    for rc in el.iter("rc"):
                        c = [x.text.strip() if x.text else "" for x in rc.findall("c")]
                        if len(c) >= 2:
                            kinds.append(int(c[1]) - 1)
                elif el.tag == "varray" and el.get("name") == "hessian":
                    hess = np.array([[float(x) for x in v.text.split()] for v in el.findall("v")])
    except (ET.ParseError, OSError, ValueError):
        if hess is None:
            return None, None
    if hess is None or hess.shape != (3 * n, 3 * n):
        return None, None
    mass = None
    if len(kinds) == n and types and max(kinds) < len(types):
        mass = np.repeat([types[k] for k in kinds], 3)
    return hess, mass


def _iter_end(ET, fh):
    """只留下要用的元素（atomtypes / atoms / hessian），其余 <calculation> 用完即清，几百 MB 的 vasprun 也不爆内存。"""
    for _ev, el in ET.iterparse(fh, events=("end",)):
        if el.tag == "array" and el.get("name") in ("atomtypes", "atoms"):
            yield el
        elif el.tag == "varray" and el.get("name") == "hessian":
            yield el
        elif el.tag == "calculation":
            el.clear()


def parse_outcar(path):
    """S6 目录的 OUTCAR（+ 同目录 vasprun.xml）-> 拆解要的全部输入；力常数与内应变给出所有候选，交给 decompose
    用"复现 VASP 离子项"来选。"""
    path = Path(path)
    txt = path.read_text(errors="ignore")
    mn = re.search(r"NIONS\s*=\s*(\d+)", txt)
    vols = re.findall(r"volume of cell\s*:\s*([\d.]+)", txt)
    if not mn or not vols:
        raise ValueError("OUTCAR 里没有 NIONS / volume of cell")
    n = int(mn.group(1))
    k_cands = []
    s, labels = _second_derivatives(txt, n)
    if s is not None:
        want = ["%d%s" % (a + 1, ax) for a in range(n) for ax in "XYZ"]
        if labels != want:
            raise ValueError("SECOND DERIVATIVES 的自由度顺序不是 1X 1Y 1Z 2X …：%s…" % labels[:6])
        k_cands.append(("OUTCAR SECOND DERIVATIVES", s))
    vr = next((path.with_name(x) for x in ("vasprun.xml", "vasprun.xml.gz") if path.with_name(x).is_file()),
              path.with_name("vasprun.xml"))
    if vr.is_file():
        h, mass = _vasprun_hessian(vr, n)
        if h is not None:
            if mass is not None:
                k_cands.append(("vasprun hessian × √(m_i m_j)", h * np.sqrt(np.outer(mass, mass))))
            k_cands.append(("vasprun hessian（不去质量加权）", h))
    if not k_cands:
        raise ValueError("拿不到力常数矩阵：OUTCAR 里没有 SECOND DERIVATIVES，%s 里也没有 <dynmat> hessian（%d×%d）"
                         % (vr.name, 3 * n, 3 * n))
    lam_cands = _internal_strain(txt, n)
    if not lam_cands:
        raise ValueError("读不到每个离子的 INTERNAL STRAIN TENSOR —— 不是 ISIF≥3 的 IBRION=6？")
    clamped = _modulus_block(txt, "SYMMETRIZED ELASTIC MODULI")
    if clamped is None:
        clamped = _modulus_block(txt, "ELASTIC MODULI (kBar)")
    total = _modulus_block(txt, "TOTAL ELASTIC MODULI")
    if clamped is None or total is None:
        raise ValueError("读不到弹性模量块：%s" % ", ".join(k for k, v in (("刚性离子", clamped), ("TOTAL", total))
                                                       if v is None))
    ionic = _modulus_block(txt, "ELASTIC MODULI CONTR FROM IONIC RELAXATION")
    src = "ELASTIC MODULI CONTR FROM IONIC RELAXATION"
    if ionic is None:                                       # VASP 6 不单独打印：TOTAL − 刚性离子
        ionic, src = total - clamped, "TOTAL − SYMMETRIZED"
    return {"n": n, "volume": float(vols[-1]), "k_cands": k_cands, "lam_cands": lam_cands,
            "clamped": clamped, "ionic": ionic, "ionic_src": src, "total": total}


def _is_pd(m):
    m = 0.5 * (np.asarray(m) + np.asarray(m).T)
    return bool(np.all(np.linalg.eigvalsh(m) > 0))


def decompose(d, dim="2d"):
    n, vol = d["n"], d["volume"]
    trans = np.zeros((3, 3 * n))
    for ax in range(3):
        trans[ax, ax::3] = 1.0 / np.sqrt(n)
    scale = max(np.max(np.abs(d["ionic"])), 1e-9)
    best = None
    for klab, kraw in d["k_cands"]:
        s = 0.5 * (kraw + kraw.T)
        k = s if np.trace(s) > 0 else -s                    # 力常数矩阵对角为正；打印的符号两种都接受
        w, u = np.linalg.eigh(k)
        is_t_k = np.sum((trans @ u) ** 2, axis=0) > TRANS_OVERLAP     # 与均匀平移的重叠
        for llab, lam0 in d["lam_cands"]:
            for f in (1.0, 0.5, 2.0):                       # 剪切列（XY YZ ZX）的约定因子
                lam = lam0.copy()
                lam[:, 3:] *= f
                p = u.T @ lam                               # (3n, 6)
                contrib = np.array([-np.outer(p[m], p[m]) / w[m] for m in range(3 * n)]) * EV_A3_TO_KBAR / vol
                # VASP 求逆时可能含平移模（数值上没去掉），也可能已经去掉：两种都试，哪种复现得上就是哪种
                for incl in (True, False):
                    full = contrib.sum(axis=0) if incl else contrib[~is_t_k].sum(axis=0)
                    err = float(np.max(np.abs(full - d["ionic"])) / scale)
                    if best is None or err < best[0]:
                        best = (err, f, contrib, incl, w, is_t_k, klab, llab)
    err, f, contrib, incl, w, is_t, klab, llab = best
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
        "force_constants_from": klab, "internal_strain_from": llab, "ionic_reference": d["ionic_src"],
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
    print("   复现 VASP 离子弛豫项（%s）：相对误差 %.1f%%%s"
          % (r["ionic_reference"], 100 * r["reproduce_rel_err"], "" if r["reproduced"] else "  ★ 复现不了"))
    print("     力常数：%s；内应变：%s；剪切列因子 %g；VASP 求逆时%s平移模"
          % (r["force_constants_from"], r["internal_strain_from"], r["shear_factor"],
             "含" if r["vasp_includes_translations"] else "已去掉"))
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
