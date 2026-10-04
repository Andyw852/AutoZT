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
MATCH_TOL = 0.05                      # 复现 VASP 离子项的相对误差上限（Frobenius）
TRANS_OVERLAP = 0.9
VASP_TO_THZ = 15.633302               # sqrt(eV/Å²/amu) -> THz
FREQ_TOL = 0.03                       # 力常数换算出的频率 vs OUTCAR 频率（|f|>1 THz 的模，中位相对误差）
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


def _outcar_freqs(txt, n):
    """OUTCAR 的 Γ 点振动频率（THz，虚频取负），取第一套 3n 个；没有 -> None。"""
    out = {}
    for m in re.finditer(r"^\s*(\d+)\s+f(/i)?\s*=\s*([\d.]+)\s+THz", txt, re.M):
        i = int(m.group(1))
        if i in out:
            if len(out) >= 3 * n:
                break
            continue
        out[i] = -float(m.group(3)) if m.group(2) else float(m.group(3))
    return np.array([out[i] for i in sorted(out)]) if len(out) == 3 * n else None


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
    mass = None
    s, labels = _second_derivatives(txt, n)
    if s is not None:
        want = ["%d%s" % (a + 1, ax) for a in range(n) for ax in "XYZ"]
        if labels != want:
            raise ValueError("SECOND DERIVATIVES 的自由度顺序不是 1X 1Y 1Z 2X …：%s…" % labels[:6])
        k_cands.append(("OUTCAR SECOND DERIVATIVES", s))
    vr = next((path.with_name(x) for x in ("vasprun.xml", "vasprun.xml.gz") if path.with_name(x).is_file()),
              path.with_name("vasprun.xml"))
    if vr.is_file():
        h, mass = _vasprun_hessian(vr, n)                   # mass 也给频率核对用
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
            "clamped": clamped, "ionic": ionic, "ionic_src": src, "total": total,
            "mass": mass, "freqs": _outcar_freqs(txt, n)}


def _is_pd(m):
    m = 0.5 * (np.asarray(m) + np.asarray(m).T)
    return bool(np.all(np.linalg.eigvalsh(m) > 0))


def _freq_calibrate(k, mass, freqs):
    """[V157] 力常数候选 k + 质量 -> 频率谱，与 OUTCAR 的频率谱比（只比 |f| > 1 THz 的模），并定出单位标量。
    -> (scale, err)：k×scale 才是 eV/Å² 的力常数；err = 定标后逐模的中位相对误差。拿不到频率/质量 -> (1.0, None)。
    VASP 写进 vasprun.xml 的 hessian 单位随版本不同：phonopy 按 eV/Å²/amu（−Φ/√(m_i m_j)）读，Mo2S3 这份是
    THz²（−hessian 的本征值直接是 f²）—— 差 15.633² ≈ 244 倍，V156 的固定换算因此两个候选都被拒。
    只允许一个标量：质量加权方式错了，谱的**形状**对不上，定标救不回来。"""
    if mass is None or freqs is None:
        return 1.0, None
    w = np.linalg.eigvalsh(k / np.sqrt(np.outer(mass, mass)))
    f = np.sort(np.sign(w) * np.sqrt(np.abs(w)) * VASP_TO_THZ)
    g = np.sort(freqs)
    sel = (np.abs(g) > 1.0) & (np.abs(f) > 1e-9)
    if not sel.any():
        return 1.0, None
    scale = float(np.median((g[sel] / f[sel]) ** 2))
    f2 = f * np.sqrt(scale)
    return scale, float(np.median(np.abs(f2[sel] - g[sel]) / np.abs(g[sel])))


def _unit_note(scale):
    for val, lab in ((1.0, "eV/Å²"), (1.0 / VASP_TO_THZ ** 2, "THz²（f²）"),
                     (1.0 / (2 * np.pi * VASP_TO_THZ) ** 2, "(2πTHz)²（ω²）")):
        if abs(scale / val - 1.0) < 0.05:
            return lab
    return "未知比例 ×%.4g" % scale


def _translations(n):
    t = np.zeros((3 * n, 3))
    for ax in range(3):
        t[ax::3, ax] = 1.0 / np.sqrt(n)
    return t


def _projected(k, lam, vol):
    """严格去掉均匀平移后的离子项：K 投影到平移的正交补，按本征模求和。-> (离子项, 每模贡献, 本征值, 非平移模下标)"""
    n3 = k.shape[0]
    t = _translations(n3 // 3)
    q = np.eye(n3) - t @ t.T
    kq = q @ k @ q
    kq = 0.5 * (kq + kq.T)
    w, u = np.linalg.eigh(kq)
    null = np.argsort(-np.sum((t.T @ u) ** 2, axis=0))[:3]
    keep = np.array([m for m in range(n3) if m not in set(null)])
    p = u.T @ lam
    contrib = {int(m): -np.outer(p[m], p[m]) / w[m] * EV_A3_TO_KBAR / vol for m in keep}
    return sum(contrib.values()), contrib, w, keep


def _fit_subspace(dmat, b):
    """dmat ≈ bᵀ M b（M 对称 3×3）的最小二乘 -> (拟合矩阵, 残差的 Frobenius 范数)。b 为 3×6。"""
    idx = [(i, j) for i in range(6) for j in range(i, 6)]
    pairs = [(a, c) for a in range(3) for c in range(a, 3)]
    a_mat = np.array([[b[a, i] * b[c, j] + (b[c, i] * b[a, j] if a != c else 0.0) for a, c in pairs]
                      for i, j in idx])
    y = np.array([dmat[i, j] for i, j in idx])
    x, *_ = np.linalg.lstsq(a_mat, y, rcond=None)
    m = np.zeros((3, 3))
    for (a, c), v in zip(pairs, x):
        m[a, c] = m[c, a] = v
    fit = b.T @ m @ b
    return fit, float(np.linalg.norm(dmat - fit))


def decompose(d, dim="2d"):
    """自检 = 能否复现 VASP 的离子项：
    ① 精确复现（VASP 求逆时已去掉平移模）；或 ② 差值 VASP − 本工具 完全落在平移子空间里（bᵀMb，M 对称 3×3、6 个自由参数
    对 21 个独立分量；b = 平移方向上的内应变）—— 即 VASP 多出来的只是平移模的贡献（数值假象，大小取决于 VASP 内部怎么处理
    近零模，无从也无须精确复现）。两者都不成立 -> 结论不可用（并看差值是不是落在某个软模方向上）。"""
    n, vol = d["n"], d["volume"]
    t = _translations(n)
    vi = d["ionic"]
    scale = max(float(np.linalg.norm(vi)), 1e-9)
    tried, best = [], None
    for klab, kraw in d["k_cands"]:
        s_ = 0.5 * (kraw + kraw.T)
        k = s_ if np.trace(s_) > 0 else -s_                 # 力常数矩阵对角为正；打印的符号两种都接受
        fscale, ferr = _freq_calibrate(k, d.get("mass"), d.get("freqs"))
        tried.append((klab, ferr, fscale))
        if ferr is not None and ferr > FREQ_TOL:
            continue
        k = k * fscale                                      # 换成 eV/Å²（频谱定出的单位）
        for llab, lam0 in d["lam_cands"]:
            for f in (1.0, 0.5, 2.0):                       # 剪切列（XY YZ ZX）的约定因子
                lam = lam0.copy()
                lam[:, 3:] *= f
                ours, contrib, w, keep = _projected(k, lam, vol)
                diff = vi - ours
                mode, fit, res = "exact", np.zeros((6, 6)), float(np.linalg.norm(diff))
                # 只用**严格的均匀平移**作基底（b = 平移方向上的内应变 = 它违反平移不变性的那部分）。
                #   不用 K 的近零本征模：它们和软光学模有混合，投影出来指向软模方向，会把软模的差值也"解释"成平移。
                basis = t.T @ lam
                if res / scale > MATCH_TOL and np.linalg.norm(basis) > 1e-6 * max(np.linalg.norm(lam), 1e-12):
                    fit_t, res_t = _fit_subspace(diff, basis)
                    if res_t < res:
                        mode, fit, res = "translation(uniform)", fit_t, res_t
                rel = res / scale
                if best is None or rel < best["rel"]:
                    best = {"rel": rel, "mode": mode, "fit": fit, "ours": ours, "contrib": contrib, "w": w,
                            "keep": keep, "klab": klab, "llab": llab, "f": f, "ferr": ferr, "diff": diff}
    if best is None:
        raise ValueError("所有力常数候选都和 OUTCAR 的振动频率对不上：%s"
                         % "；".join("%s 定标后频率误差 %s（单位按 %s）" % (a, "%.0f%%" % (100 * b) if b is not None else "-",
                                                                     _unit_note(c)) for a, b, c in tried))
    idx = list(INPLANE) if dim == "2d" else list(range(6))
    sub = lambda m: np.asarray(m)[np.ix_(idx, idx)]         # noqa: E731
    proj_total = d["clamped"] + best["ours"]
    total_pd, proj_pd = _is_pd(sub(d["total"])), _is_pd(sub(proj_total))
    worst = min(idx, key=lambda i: d["total"][i, i])        # 翻负的那个分量：所查块里最负的对角元
    contrib, w = best["contrib"], best["w"]
    top = sorted(contrib, key=lambda m: contrib[m][worst, worst])[:3]
    # 差值是不是落在某个软模方向上（软模刚度对数值敏感时的另一种解释）
    soft = None
    if best["rel"] > MATCH_TOL:
        for m in sorted(contrib, key=lambda m: w[m])[:5]:
            pm = contrib[m] / max(abs(contrib[m]).max(), 1e-12)
            c = float(np.sum(best["diff"] * pm) / max(np.sum(pm * pm), 1e-12))
            r_ = float(np.linalg.norm(best["diff"] - c * pm)) / scale
            if soft is None or r_ < soft["rel"]:
                soft = {"mode": int(m), "eigenvalue_eV_A2": round(float(w[m]), 4), "rel": round(r_, 3)}
    res = {
        "dim": dim, "n_ions": n, "volume_A3": vol, "ionic_reference": d["ionic_src"],
        "force_constants_from": best["klab"], "internal_strain_from": best["llab"], "shear_factor": best["f"],
        "freq_check": {a: {"err": (round(b, 4) if b is not None else None), "scale": c, "unit": _unit_note(c)}
                       for a, b, c in tried},
        "match_mode": best["mode"], "reproduce_rel_err": round(best["rel"], 4),
        "reproduced": best["rel"] <= MATCH_TOL,
        "worst_component": VASP_ORDER[worst],
        "clamped": round(float(d["clamped"][worst, worst]), 2),
        "vasp_ionic": round(float(vi[worst, worst]), 2),
        "other_modes_part": round(float(best["ours"][worst, worst]), 2),
        "translation_part": round(float(vi[worst, worst] - best["ours"][worst, worst]), 2),
        "top_other_modes": [{"mode": int(m), "eigenvalue_eV_A2": round(float(w[m]), 5),
                             "contrib_kbar": round(float(contrib[m][worst, worst]), 2)} for m in top],
        "soft_mode_fit": soft,
        "total_pd": total_pd, "projected_pd": proj_pd,
        "projected_total_vasp_order_kbar": np.round(proj_total, 3).tolist(),
    }
    if not res["reproduced"]:
        res["verdict"] = "unreliable"
    elif proj_pd and not total_pd:
        res["verdict"] = "artifact" if best["mode"] != "exact" else "unstable"
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
    print("   力常数候选与 OUTCAR 振动频率谱核对（定标后逐模中位误差）：%s"
          % "；".join("%s %s" % (a, ("%.1f%%（单位按 %s）" % (100 * b["err"], b["unit"])) if b["err"] is not None
                                 else "无法核对") for a, b in r["freq_check"].items()))
    print("   复现 VASP 离子弛豫项（%s）：%s，残差 %.1f%%%s"
          % (r["ionic_reference"], {"exact": "直接复现（VASP 已去掉平移模）",
                                    "translation(uniform)": "差值完全落在平移子空间（均匀平移）"}[r["match_mode"]],
             100 * r["reproduce_rel_err"], "" if r["reproduced"] else "  ★ 复现不了"))
    print("     选中：力常数 %s；内应变 %s；剪切列因子 %g" % (r["force_constants_from"], r["internal_strain_from"],
                                                        r["shear_factor"]))
    print("   %s：刚性离子 %+.2f，VASP 离子弛豫 %+.2f = 非平移模 %+.2f + 平移方向 %+.2f"
          % (c, r["clamped"], r["vasp_ionic"], r["other_modes_part"], r["translation_part"]))
    for t in r["top_other_modes"]:
        print("     非平移模 #%d：本征值 %.4f eV/Å²，贡献 %+.2f kBar" % (t["mode"], t["eigenvalue_eV_A2"], t["contrib_kbar"]))
    print("   扣掉平移方向后：%s = %+.2f kBar；%s %s"
          % (c, r["clamped"] + r["other_modes_part"], "面内 3×3" if dim == "2d" else "6×6",
             "正定" if r["projected_pd"] else "仍不正定"))
    v = r["verdict"]
    if v == "unreliable":
        sm = r.get("soft_mode_fit")
        print("[ERROR] ★ 复现不了 VASP 的离子弛豫项：VASP − 本工具 的差值既不是 0，也不落在平移子空间里，以上结论不可用。")
        if sm and sm["rel"] <= 2 * MATCH_TOL:
            print("        差值主要落在软模 #%d（本征值 %.4f eV/Å²）方向：软模刚度对数值很敏感，单靠后处理分不清，"
                  "要更严的 S6（EDIFF 更小 / POTIM 更合适）或 DFPT 复核" % (sm["mode"], sm["eigenvalue_eV_A2"]))
        return 2
    if v == "artifact":
        print("[OK] 结论：数值假象 —— VASP 离子项里多出的 %+.2f kBar 完全落在平移方向（均匀平移和应变不耦合，这部分本该是 0）。"
              "扣掉后的弛豫张量正定，不用重算 S6；投影后的 TOTAL 写在 elastic_ionic_decompose.json" % r["translation_part"])
        return 0
    if v == "unstable":
        print("[WARN] ★ 结论：扣掉平移方向仍不正定 —— 真实的内应变耦合（贡献最大的非平移模见上，Born 不稳定或软模），"
              "不可用；重算 S6 不会改变")
        return 1
    print("[OK] 结论：TOTAL 本来就正定（或扣平移方向前后都正定）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
