#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""symmetry_audit.py —— 结构体检第五条：数值微畸变审计 + spglib 对称化（公共池）

背景（2026-09-20，MoS2 线 A 实证）
----------------------------------
S1 弛豫后的原胞常有"数值微畸变"：a=3.140197057514、b=3.140189175214 Å
（差 7.88e-6 Å）、gamma=120.000083034652 度。后果：
  * spglib 在 symprec=1e-5 只给 Amm2 (#38)/4 ops（pheasy 沿用这个判断）；
  * symprec>=1e-4 才给真正的 P-6m2 (#187)/12 ops。
对称操作少 2/3 -> 独立力常数大增 -> 允许六方破缺 -> ZA 被算成线性 (p=1.186)；
把晶格拉平 a=b、锁 gamma=120 后，连 1e-5 都立即 P-6m2、ZA 恢复二次 (p=1.993)。
=> 畸变在几何上；进 S4 之前必须对称化。

判据（user 指定，不要固定某个容差当判据）
---------------------------------------------
对结构做 spglib 空间群审计，**symprec=1e-5 与 1e-4 给出的空间群不一致**
=> 判定存在数值微畸变 => 触发对称化建议/动作。

对称化的三项确认（任一不过就报错/告警，绝不静默通过）
----------------------------------------------------
  ① 原子数与化学计量比不变（refine_cell 可能改变原胞选取，原子数会变）；
  ② 对称化前后的单点能量差在数值噪声以内（本模块只做解析对比，
     真正的两个单点由流水线提供；未提供时状态= pending，显式报出，成本 2 个单点）；
  ③ 真空轴是否被移动（refine_cell 可能改变晶格矢量取向/原胞选取；
     若真空轴换了矢量或真空中心移了，必须明确报出，下游认 vac_axis 的代码要跟着走）。

另外：对称化后重新检查面内残余应力（拉平 a/b 会引入一点应力，量级应很小）。

用法（CLI，可被 autozt 步骤/gen/检查脚本调用）
-----------------------------------------------
    # 1) 只审计（默认，退出码 0；micro_distortion=True 时打印 WARN）
    python3 symmetry_audit.py --poscar POSCAR

    # 2) 审计 + 对称化 + 三项确认，两种结构都留档
    python3 symmetry_audit.py --poscar POSCAR --mode symmetrize \
        --out POSCAR_sym --archive-dir symmetry_audit \
        --energy-before OUTCAR.before --energy-after OUTCAR.after

    # 3) 对称化后就地覆盖 POSCAR（原结构备份为 POSCAR.pre_symmetry）
    python3 symmetry_audit.py --poscar POSCAR --mode symmetrize --inplace ...

    # 4) 单独做两个单点的能量对比（也可对上一步 S2/S4 已有的 OUTCAR 复用）
    python3 symmetry_audit.py --mode energy-check \
        --energy-before OUTCAR.a --energy-after OUTCAR.b --natoms 3

退出码：audit=0（除非 --strict-exit）；symmetrize=0 成功；
       3=三项确认未过 / 能量确认缺失且未 --allow-energy-pending；4=能量对比未过；
       5=spglib 缺失 / 解析失败。

开关（step.conf [params]，登记进 relax_common.CONF_SPEC）
--------------------------------------------------------
    SYMMETRY_AUDIT            = warn   # off | warn | error   （S1 gen 体检第五条）
    SYMMETRY_SYMMETRIZE       = off    # off | on            （对称化默认关，显式触发）
    SYMMETRY_AUDIT_SYMPREC    = 1e-4   # 对称化用的容差（必须取两容差中识别出高对称的那个）
    SYMMETRY_ENERGY_TOL_MEV   = 1.0    # 能量确认阈值 meV/atom
读法：CLI --switch 显式指定 > 环境变量 AUTOZT_SYMMETRY_* > step.conf [params] > 默认。

放置：skill/_common/opt/symmetry_audit.py。顶层只依赖标准库；numpy/spglib
在函数内**懒加载**，缺了只让审计"不可用"，绝不让调用方 import 失败。
被 relax_common.py 静态 import（触发公共池同目录依赖自动补推）。
"""
import argparse
import datetime
import json
import math
import os
import re
import sys
from collections import Counter
from pathlib import Path

__all__ = [
    "read_poscar", "write_poscar", "lattice_metrics", "vacuum_info",
    "spacegroup_audit", "symmetrize", "confirmations", "audit_structure",
    "audit_structure_file", "parse_outcar_energy", "parse_inplane_stress_kb",
    "parse_external_pressure_kb", "energy_check", "stress_recheck",
    "read_switches", "symmetrize_enabled", "spglib_available", "main",
]

AXIS_NAMES = ("a", "b", "c")

# 判据用的两个容差：1e-5 与 1e-4（user 指定；不要固定单一容差）
SYMPREC_PRIMARY = 1e-5
SYMPREC_LOOSE = 1e-4
AUDIT_SYMPRECS = (SYMPREC_PRIMARY, SYMPREC_LOOSE)
DEFAULT_SYMPREC = SYMPREC_LOOSE
DEFAULT_ENERGY_TOL_MEV_PER_ATOM = 1.0   # meV/atom
VACUUM_MIN = 8.0                        # Å，与 dim_common 一致
DEFAULT_VAC_CENTER_TOL = 0.02           # 真空中心分数位移报警阈值

# ---------------------------------------------------------------------------
# 元素符号 <-> 原子序数（自包含，不 import pymatgen）
# ---------------------------------------------------------------------------
_ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb "
    "Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
    "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es "
    "Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()
_SYMBOL_TO_Z = {s: i + 1 for i, s in enumerate(_ELEMENTS)}
_Z_TO_SYMBOL = {i + 1: s for i, s in enumerate(_ELEMENTS)}


def symbol_to_z(symbols):
    out = []
    for s in symbols:
        key = str(s).strip().capitalize()
        if key not in _SYMBOL_TO_Z:
            raise ValueError("未知元素符号 %r（POSCAR 元素行？）" % (s,))
        out.append(_SYMBOL_TO_Z[key])
    return out


def z_to_symbol(z):
    return _Z_TO_SYMBOL.get(int(z), "X%d" % int(z))


# ---------------------------------------------------------------------------
# 3x3 线代（纯标准库）
# ---------------------------------------------------------------------------
def _cross(u, v):
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


def _dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _norm(u):
    return _dot(u, u) ** 0.5


def _det3(m):
    return _dot(m[0], _cross(m[1], m[2]))


def _inv3(m):
    d = _det3(m)
    if abs(d) < 1e-12:
        raise ValueError("晶格矩阵奇异（体积为 0）")
    c0, c1, c2 = _cross(m[1], m[2]), _cross(m[2], m[0]), _cross(m[0], m[1])
    return [[c0[0] / d, c1[0] / d, c2[0] / d],
            [c0[1] / d, c1[1] / d, c2[1] / d],
            [c0[2] / d, c1[2] / d, c2[2] / d]]


# ---------------------------------------------------------------------------
# POSCAR 读写（VASP4/5，Direct/Cartesian，负缩放）
# ---------------------------------------------------------------------------
def _read_potcar_symbols(potcar):
    """从 POTCAR 里按 "TITEL = ... PAW_<El>" 或 VRHFIN 提取元素序列（best-effort）。"""
    syms, txt = [], ""
    try:
        txt = Path(potcar).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    for line in txt.splitlines():
        m = re.search(r"VRHFIN\s*=\s*([A-Za-z]{1,2})", line)
        if not m:
            m = re.search(r"PAW[_\s]?([A-Z][a-z]?)", line)
        if m:
            syms.append(m.group(1).capitalize())
    return syms or None


def read_poscar(path, potcar=None):
    """读 POSCAR/CONTCAR，返回 dict：
       {lattice(3x3 已含缩放), frac, symbols(逐原子), title, counts, species}。
    VASP4（无元素行）时用 potcar=/同目录 POTCAR 兜底；再不行报错。"""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError("找不到结构文件：%s" % p)
    lines = p.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 8:
        raise ValueError("POSCAR 行数不足：%s" % p)
    title = lines[0].rstrip()
    scale = float(lines[1].split()[0])
    lattice = [[float(x) for x in lines[2 + i].split()[:3]] for i in range(3)]
    if scale < 0:                       # 负数 = 目标体积
        vol = abs(_det3(lattice))
        scale = (abs(scale) / vol) ** (1.0 / 3.0)
    lattice = [[scale * x for x in row] for row in lattice]

    idx = 5
    tokens = lines[idx].split()
    if not tokens:
        raise ValueError("POSCAR 第 6 行为空")
    species = None
    if not re.fullmatch(r"[+-]?\d+", tokens[0]):
        species = tokens
        idx += 1
    counts = [int(x) for x in lines[idx].split()]
    natoms = sum(counts)
    idx += 1
    if lines[idx].strip()[:1].upper() == "S":       # Selective dynamics
        idx += 1
    mode = lines[idx].strip()[:1].upper()
    idx += 1
    coords = []
    for ln in lines[idx:idx + natoms]:
        t = ln.split()
        if len(t) < 3:
            raise ValueError("POSCAR 坐标行不完整: %r" % ln)
        coords.append([float(t[0]), float(t[1]), float(t[2])])
    if len(coords) != natoms:
        raise ValueError("POSCAR 坐标行数 (%d) != 原子数 (%d)" % (len(coords), natoms))
    if mode in ("C", "K"):
        inv = _inv3(lattice)            # frac = cart · lat^-1（行向量）
        cart = [[scale * x for x in c] for c in coords]
        coords = [[_dot(c, [inv[0][j], inv[1][j], inv[2][j]]) for j in range(3)]
                  for c in cart]

    if species is None:                 # VASP4：找 POTCAR 兜底
        if potcar is None:
            potcar = p.with_name("POTCAR")
        pot_symbols = _read_potcar_symbols(potcar) if potcar else None
        if pot_symbols and len(pot_symbols) == len(counts):
            species = pot_symbols
    if species is None:
        raise ValueError("POSCAR 无元素符号行，且没有可用的 POTCAR 兜底：%s" % p)
    if len(species) != len(counts):
        raise ValueError("元素种数 %d != 计数个数 %d" % (len(species), len(counts)))
    symbols = []
    for s, c in zip(species, counts):
        symbols += [s.capitalize()] * int(c)
    return {"lattice": lattice, "frac": coords, "symbols": symbols,
            "title": title, "counts": list(counts),
            "species": [s.capitalize() for s in species], "path": str(p)}


def write_poscar(path, lattice, frac, symbols, title="symmetrized"):
    """按 VASP5 写 POSCAR（Direct，按元素分组，组内保持输入顺序）。"""
    order = []
    for s in symbols:
        if s not in order:
            order.append(s)
    counts = [symbols.count(s) for s in order]
    out = [str(title)[:120], "    1.0000000000000000"]
    for row in lattice:
        out.append(" %21.16f %21.16f %21.16f" % (float(row[0]), float(row[1]), float(row[2])))
    out.append("  " + "  ".join(order))
    out.append("  " + "  ".join(str(c) for c in counts))
    out.append("Direct")
    fmap = {}
    for s, f in zip(symbols, frac):
        fmap.setdefault(s, []).append(f)
    for s in order:
        for f in fmap[s]:
            out.append(" %21.16f %21.16f %21.16f" % (float(f[0]) % 1.0,
                                                     float(f[1]) % 1.0,
                                                     float(f[2]) % 1.0))
    Path(path).write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


# ---------------------------------------------------------------------------
# 几何量：晶格度量 / 真空
# ---------------------------------------------------------------------------
def lattice_metrics(lattice):
    a = _norm(lattice[0]); b = _norm(lattice[1]); c = _norm(lattice[2])

    def ang(u, v):
        d = _norm(u) * _norm(v)
        return math.degrees(math.acos(max(-1.0, min(1.0, _dot(u, v) / d)))) if d else 0.0
    return {"a": a, "b": b, "c": c,
            "alpha": ang(lattice[1], lattice[2]),
            "beta": ang(lattice[0], lattice[2]),
            "gamma": ang(lattice[0], lattice[1])}


def _perp_height(lattice, axis):
    j, k = (axis + 1) % 3, (axis + 2) % 3
    vol = abs(_det3(lattice))
    area = _norm(_cross(lattice[j], lattice[k]))
    return vol / area if area > 1e-12 else 0.0


def vacuum_per_axis(lattice, frac):
    """沿 a/b/c 三方向的真空厚度 (Å)，与 dim_common.vacuum_per_axis 同义。"""
    out = []
    for i in range(3):
        height = _perp_height(lattice, i)
        f = sorted((c[i] % 1.0) for c in frac)
        if len(f) == 1:
            gap = 1.0
        else:
            gaps = [f[n + 1] - f[n] for n in range(len(f) - 1)]
            gaps.append(f[0] + 1.0 - f[-1])
            gap = max(gaps)
        out.append(gap * height)
    return out


def vacuum_info(lattice, frac, vacuum_min=VACUUM_MIN):
    """返回真空轴信息：axis 下标、真空厚度、垂直胞高、真空中心分数坐标。"""
    vacs = vacuum_per_axis(lattice, frac)
    axis = max(range(3), key=lambda i: vacs[i])
    height = _perp_height(lattice, axis)
    f = sorted((c[axis] % 1.0) for c in frac)
    gaps, centers = [], []
    for n in range(len(f) - 1):
        gaps.append(f[n + 1] - f[n])
        centers.append((f[n] + f[n + 1]) / 2.0)
    gaps.append(f[0] + 1.0 - f[-1])
    centers.append(((f[-1] + f[0] + 1.0) / 2.0) % 1.0)
    gi = max(range(len(gaps)), key=lambda i: gaps[i])
    return {"axis": int(axis), "axis_name": AXIS_NAMES[axis],
            "vacuum_A": gaps[gi] * height, "height_A": height,
            "center_frac": centers[gi], "vacuums_A": vacs,
            "is_2d": vacs[axis] >= vacuum_min}


# ---------------------------------------------------------------------------
# spglib 懒加载
# ---------------------------------------------------------------------------
def _import_numpy():
    try:
        import numpy as np
        return np
    except ImportError as exc:          # noqa: BLE001
        raise RuntimeError("需要 numpy：%s" % exc)


def _import_spglib():
    try:
        import spglib
        return spglib
    except ImportError as exc:          # noqa: BLE001
        raise RuntimeError("需要 spglib（pip install spglib）：%s" % exc)


def spglib_available():
    try:
        _import_spglib()
        return True
    except RuntimeError:
        return False


def _spglib_version():
    try:
        import spglib
        return getattr(spglib, "__version__", None)
    except Exception:                   # noqa: BLE001
        return None


def spacegroup_audit(lattice, frac, symbols, symprecs=AUDIT_SYMPRECS):
    """逐容差求空间群。返回 (entries, micro_distortion)。

    micro_distortion = 两个容差给出的 (空间群号, 符号) 不一致。"""
    spg = _import_spglib()
    np = _import_numpy()
    numbers = symbol_to_z(symbols)
    cell = (np.array(lattice, dtype=float), np.array(frac, dtype=float), numbers)
    entries = []
    for sp in symprecs:
        ds = spg.get_symmetry_dataset(cell, symprec=float(sp))
        if ds is None:
            entries.append({"symprec": float(sp), "number": None,
                            "international": None, "n_ops": 0, "ok": False})
            continue
        entries.append({
            "symprec": float(sp),
            "number": int(ds.number),
            "international": str(ds.international),
            "n_ops": int(len(ds.rotations)),
            "ok": True,
        })
    ids = {(e["number"], e["international"]) for e in entries if e["ok"]}
    return entries, len(ids) > 1


def symmetrize(lattice, frac, symbols, symprec=DEFAULT_SYMPREC):
    """用 spglib.refine_cell 产生精确对称的胞。

    返回 (new_lattice, new_frac, new_symbols, used_symprec)。refine_cell 会
    理想化晶格与原子位置；它可能改变原胞选取/原子顺序，调用方必须做三项确认。"""
    spg = _import_spglib()
    np = _import_numpy()
    numbers = symbol_to_z(symbols)
    cell = (np.array(lattice, dtype=float), np.array(frac, dtype=float), numbers)
    res = spg.refine_cell(cell, symprec=float(symprec))
    if res is None:
        raise RuntimeError("spglib.refine_cell 失败（symprec=%g）" % symprec)
    rl, rp, rn = res
    rl = np.asarray(rl, dtype=float)
    rp = np.mod(np.asarray(rp, dtype=float), 1.0)
    new_symbols = [z_to_symbol(z) for z in list(rn)]
    return rl, rp, new_symbols, float(symprec)


# ---------------------------------------------------------------------------
# 三项确认
# ---------------------------------------------------------------------------
def confirmations(before, after, energy_before=None, energy_after=None,
                  energy_tol_mev_per_atom=DEFAULT_ENERGY_TOL_MEV_PER_ATOM,
                  vac_center_tol=DEFAULT_VAC_CENTER_TOL):
    """三项确认。before/after = dict(lattice, frac, symbols)。energy_* 可为
    OUTCAR 路径或 eV 数值（None 表示"尚未跑，pending"）。

    返回 dict，其中：
      c1_stoichiometry : {ok, before, after}
      c2_energy        : {status: pass|fail|pending, ...}
      c3_vacuum_axis   : {ok, moved, ...}
      all_ok           : 几何两项通过 且 能量不是 fail（pending 不算 pass）
      fully_confirmed  : 几何两项通过 且 能量 pass
    """
    c1_before = dict(Counter(before["symbols"]))
    c1_after = dict(Counter(after["symbols"]))
    c1_ok = (c1_before == c1_after)
    c1 = {"ok": bool(c1_ok), "before": c1_before, "after": c1_after,
          "natoms_before": len(before["symbols"]),
          "natoms_after": len(after["symbols"])}

    vi_b = vacuum_info(before["lattice"], before["frac"])
    vi_a = vacuum_info(after["lattice"], after["frac"])
    axis_moved = (vi_b["axis"] != vi_a["axis"])
    # 真空中心分数位移（只在同轴时才有意义）。真空方向的【整体平移】在 PBC 下物理
    # 等价（真空隙不变、只是 slab 在胞内换了个位置），故中心位移只记录、不判失败；
    # 真正要拦的是真空轴被改（axis_moved）——refine_cell 可能重定向晶格矢量。
    dcenter = None
    if not axis_moved:
        d = abs(vi_a["center_frac"] - vi_b["center_frac"])
        dcenter = min(d, 1.0 - d)
    c3_ok = (not axis_moved)
    c3 = {"ok": bool(c3_ok), "axis_moved": bool(axis_moved),
          "before": vi_b, "after": vi_a,
          "center_shift_frac": dcenter, "center_tol": vac_center_tol,
          "center_shift_note": "真空方向平移不计入判据（PBC 下物理等价）"}

    en_b = _as_energy(energy_before)
    en_a = _as_energy(energy_after)
    natoms = len(after["symbols"])
    if en_b is None or en_a is None:
        c2 = {"status": "pending", "energy_before_eV": en_b, "energy_after_eV": en_a,
              "note": "需要两个单点（对称化前后各一个）才能完成本项；成本 = 2 个单点。"}
    else:
        de = en_a - en_b
        de_per_atom_meV = de / max(1, natoms) * 1000.0
        c2_ok = abs(de_per_atom_meV) <= energy_tol_mev_per_atom
        c2 = {"status": "pass" if c2_ok else "fail",
              "energy_before_eV": en_b, "energy_after_eV": en_a,
              "delta_eV": de, "delta_meV_per_atom": de_per_atom_meV,
              "tol_meV_per_atom": energy_tol_mev_per_atom,
              "note": "对称化前后单点能量差在数值噪声以内" if c2_ok
                      else "能量差超过噪声阈值 —— 不是单纯拉平微畸变，禁止静默通过"}

    all_ok = bool(c1["ok"] and c3["ok"] and c2["status"] != "fail")
    fully = bool(c1["ok"] and c3["ok"] and c2["status"] == "pass")
    return {"c1_stoichiometry": c1, "c2_energy": c2, "c3_vacuum_axis": c3,
            "all_ok": all_ok, "fully_confirmed": fully}


def _as_energy(value):
    """数值直接返回；字符串当 OUTCAR 路径解析；None 返回 None。"""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return parse_outcar_energy(value)


# ---------------------------------------------------------------------------
# OUTCAR 解析：能量 / 应力
# ---------------------------------------------------------------------------
_ENERGY_RE = re.compile(r"energy\(sigma->0\)\s*=\s*([-+0-9.eEdD]+)")
_ENERGY_FB_RE = re.compile(r"free\s+energy\s+TOTEN\s*=\s*([-+0-9.eEdD]+)")


def parse_outcar_energy(path):
    """OUTCAR 最后一个 energy(sigma->0)（eV）；退而求 TOTEN。读不到返回 None。"""
    try:
        txt = Path(path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    hits = _ENERGY_RE.findall(txt)
    if not hits:
        hits = _ENERGY_FB_RE.findall(txt)
    if not hits:
        return None
    try:
        return float(hits[-1].replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


_STRESS_ROW_RE = re.compile(r"^\s*in\s+kB[ \t]+([-+0-9.eEdD \t]+)$", re.M)
_PRESS_RE = re.compile(r"external pressure\s*=\s*([-+0-9.eEdD]+)\s*kB")


def parse_inplane_stress_kb(path):
    """OUTCAR 最后一行 "in kB  xx yy zz xy yz zx" 的面内分量 max(|xx|,|yy|,|xy|) (kB)。

    这与 relax_common.run_relax.sh 的 _last_inplane 完全同义（列序 XX YY ZZ XY YZ ZX）。"""
    try:
        txt = Path(path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    rows = _STRESS_ROW_RE.findall(txt)
    if not rows:
        return None
    toks = rows[-1].split()
    if len(toks) < 4:
        return None
    try:
        vals = [float(x.replace("D", "E").replace("d", "e")) for x in toks[:6]]
    except ValueError:
        return None
    # rows 正则抓了整行数字，toks 可能含 6 个数：xx yy zz xy yz zx
    if len(vals) < 4:
        return None
    return max(abs(vals[0]), abs(vals[1]), abs(vals[3]))


def parse_external_pressure_kb(path):
    try:
        txt = Path(path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    hits = _PRESS_RE.findall(txt)
    if not hits:
        return None
    try:
        return float(hits[-1].replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def energy_check(before, after, natoms=None, tol_mev_per_atom=DEFAULT_ENERGY_TOL_MEV_PER_ATOM):
    """两个单点能量对比。before/after 可为 OUTCAR 路径或 eV 数值。"""
    e_b, e_a = _as_energy(before), _as_energy(after)
    if e_b is None or e_a is None:
        return {"status": "pending", "energy_before_eV": e_b, "energy_after_eV": e_a,
                "note": "缺能量值/OUTCAR（需要 2 个单点）"}
    n = int(natoms) if natoms else 1
    de = e_a - e_b
    de_mev = de / max(1, n) * 1000.0
    ok = abs(de_mev) <= float(tol_mev_per_atom)
    return {"status": "pass" if ok else "fail", "energy_before_eV": e_b,
            "energy_after_eV": e_a, "delta_eV": de, "delta_meV_per_atom": de_mev,
            "tol_meV_per_atom": float(tol_mev_per_atom), "natoms": n}


def stress_recheck(before_outcar, after_outcar, tol_kb=0.2):
    """对称化后重新检查面内残余应力（拉平 a/b 引入的应力应很小）。

    tol_kb 用 S1 判据的胞口径阈值 INPLANE_STRESS_TOL_KB（默认 0.2 kB）。
    """
    ip_b = parse_inplane_stress_kb(before_outcar)
    ip_a = parse_inplane_stress_kb(after_outcar)
    p_b = parse_external_pressure_kb(before_outcar)
    p_a = parse_external_pressure_kb(after_outcar)
    if ip_a is None:
        return {"status": "pending", "inplane_before_kB": ip_b, "inplane_after_kB": None,
                "note": "对称化后的单点 OUTCAR 缺失/无应力行 —— 需跑 1 个单点复检面内应力"}
    ok = abs(ip_a) <= float(tol_kb)
    return {"status": "pass" if ok else "fail",
            "inplane_before_kB": ip_b, "inplane_after_kB": ip_a,
            "external_pressure_before_kB": p_b, "external_pressure_after_kB": p_a,
            "tol_kB": float(tol_kb),
            "note": "面内残余应力在阈值内" if ok else "面内残余应力超阈值，需再弛豫/复核"}


# ---------------------------------------------------------------------------
# 开关：CLI > env > step.conf [params] > 默认
# ---------------------------------------------------------------------------
_DEFAULT_SWITCHES = {
    "SYMMETRY_AUDIT": "warn",
    "SYMMETRY_SYMMETRIZE": "off",
    "SYMMETRY_AUDIT_SYMPREC": str(DEFAULT_SYMPREC),
    "SYMMETRY_ENERGY_TOL_MEV": str(DEFAULT_ENERGY_TOL_MEV_PER_ATOM),
}
_CONF_KV = re.compile(r"^\s*([A-Z][A-Z0-9_]*)\s*=\s*([^#\n]*?)\s*(?:#.*)?$")


def read_switches(step_conf=None, env=None):
    """读开关。step_conf 可为路径或 None（自动在 cwd 找 step.conf）。"""
    env = os.environ if env is None else env
    out = dict(_DEFAULT_SWITCHES)
    path = None
    if step_conf:
        path = Path(step_conf)
    else:
        cand = Path.cwd() / "step.conf"
        if cand.is_file():
            path = cand
    if path and path.is_file():
        for line in path.read_text(encoding="utf-8-sig", errors="ignore").splitlines():
            m = _CONF_KV.match(line)
            if m and m.group(1) in out:
                out[m.group(1)] = m.group(2).strip().strip("'\"")
    for key in list(out):
        envname = "AUTOZT_" + key
        if env.get(envname):
            out[key] = env[envname]
    return out


def symmetrize_enabled(switches=None):
    """供下一步 gen 消费：只有 SYMMETRY_SYMMETRIZE=on 才自动对称化（默认 off）。

    CLI 的 --mode symmetrize 是显式触发，不依赖本开关；本函数用于"进 S4 前的那个
    gen 想自动判断是否对称化"的场景，避免每个调用方各写一套字符串解析。"""
    sw = read_switches() if switches is None else switches
    return str(sw.get("SYMMETRY_SYMMETRIZE", "off")).strip().lower() in (
        "on", "true", "1", "yes")


# ---------------------------------------------------------------------------
# 顶层审计
# ---------------------------------------------------------------------------
def audit_structure(lattice, frac, symbols, symprecs=AUDIT_SYMPRECS,
                    symprec=DEFAULT_SYMPREC):
    """只读审计：晶格度量 + 真空 + 双容差空间群 + 微畸变判定。"""
    entries, micro = spacegroup_audit(lattice, frac, symbols, symprecs)
    m = lattice_metrics(lattice)
    vi = vacuum_info(lattice, frac)
    ref = next((e for e in entries if abs(e["symprec"] - symprec) < 1e-12), entries[-1])
    return {
        "micro_distortion": bool(micro),
        "audit_symprecs": list(symprecs),
        "spacegroup": entries,
        "micro_distortion_detail": {
            "a_minus_b_A": m["a"] - m["b"],
            "a_over_b_minus_1": (m["a"] / m["b"] - 1.0) if m["b"] else None,
            "gamma_minus_120_deg": m["gamma"] - 120.0,
            "alpha_minus_90_deg": m["alpha"] - 90.0,
            "beta_minus_90_deg": m["beta"] - 90.0,
        },
        "lattice_metrics": m,
        "vacuum": vi,
        "n_atoms": len(symbols),
        "composition": dict(Counter(symbols)),
        "symmetrize_symprec": float(ref["symprec"]) if ref else float(symprec),
        "spglib_version": _spglib_version(),
    }


def audit_structure_file(path, symprecs=AUDIT_SYMPRECS, symprec=DEFAULT_SYMPREC):
    """从文件审计；spglib 不可用返回 {available: False, ...}（不抛）。"""
    if not spglib_available():
        return {"available": False, "reason": "spglib 不可用（未安装）",
                "path": str(path)}
    st = read_poscar(path)
    res = audit_structure(st["lattice"], st["frac"], st["symbols"],
                          symprecs=symprecs, symprec=symprec)
    res["available"] = True
    res["path"] = str(path)
    return res


def format_audit_report(audit, tag="audit"):
    """把审计结果排成人类可读多行文本（验收要求：打印 a/b/gamma、空间群、原子数、真空轴）。"""
    m = audit["lattice_metrics"]
    L = []
    L.append("[symmetry_audit:%s] n_atoms=%d composition=%s"
             % (tag, audit["n_atoms"], audit["composition"]))
    L.append("  a=%.12f  b=%.12f  c=%.12f" % (m["a"], m["b"], m["c"]))
    L.append("  alpha=%.12f  beta=%.12f  gamma=%.12f"
             % (m["alpha"], m["beta"], m["gamma"]))
    L.append("  |a-b|=%.3e Å  gamma-120=%.3e deg"
             % (abs(m["a"] - m["b"]), m["gamma"] - 120.0))
    for e in audit["spacegroup"]:
        L.append("  symprec=%-8g -> %s (#%s)  %s ops"
                 % (e["symprec"], e["international"], e["number"], e["n_ops"]))
    v = audit["vacuum"]
    L.append("  真空轴=%s(下标%d)  真空=%.4f Å  垂直胞高=%.4f Å  中心(frac)=%.6f"
             % (v["axis_name"], v["axis"], v["vacuum_A"], v["height_A"],
                v["center_frac"]))
    L.append("  micro_distortion=%s" % audit["micro_distortion"])
    return "\n".join(L)


# ---------------------------------------------------------------------------
# 对称化主流程 + 留档
# ---------------------------------------------------------------------------
def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def symmetrize_structure_file(poscar, out=None, archive_dir=None, symprec=DEFAULT_SYMPREC,
                              energy_before=None, energy_after=None,
                              energy_tol_mev_per_atom=DEFAULT_ENERGY_TOL_MEV_PER_ATOM,
                              inplace=False, allow_energy_pending=False,
                              stress_before=None, stress_after=None, stress_tol_kb=0.2,
                              switch=None, source_note=None, write_prov=True):
    """审计 + 对称化 + 三项确认 + 留档。返回 result dict；确认不过会 raise SystemExit。

    留档（archive_dir，默认 <poscar 同目录>/symmetry_audit/）：
      POSCAR.before / POSCAR.symmetrized / symmetry_audit.json / symmetry_audit.log
    inplace=True 时把对称化结构写回 POSCAR，原结构另存 POSCAR.pre_symmetry。"""
    poscar = Path(poscar)
    st = read_poscar(poscar)
    audit_before = audit_structure(st["lattice"], st["frac"], st["symbols"], symprec=symprec)
    audit_before["path"] = str(poscar)

    used_symprec = audit_before["symmetrize_symprec"]
    if not audit_before["micro_distortion"]:
        # 空间群一致：不触发对称化（但仍可显式因 SYMMETRY_SYMMETRIZE=on 继续，见 CLI）
        result = {
            "triggered": False,
            "reason": "symprec=1e-5 与 1e-4 空间群一致，未检出数值微畸变",
            "audit_before": audit_before, "used_symprec": used_symprec,
            "confirmations": None,
        }
        if out:
            write_poscar(out, st["lattice"], st["frac"], st["symbols"],
                         title="%s (unchanged: no micro-distortion)" % st["title"][:80])
            result["out"] = str(out)
        return result

    rl, rp, rsyms, used_symprec = symmetrize(st["lattice"], st["frac"], st["symbols"],
                                             symprec=symprec)
    after = {"lattice": rl.tolist() if hasattr(rl, "tolist") else rl,
             "frac": rp.tolist() if hasattr(rp, "tolist") else rp,
             "symbols": rsyms}
    audit_after = audit_structure(after["lattice"], after["frac"], after["symbols"],
                                  symprec=symprec)

    conf = confirmations({"lattice": st["lattice"], "frac": st["frac"], "symbols": st["symbols"]},
                         after, energy_before=energy_before, energy_after=energy_after,
                         energy_tol_mev_per_atom=energy_tol_mev_per_atom)

    stress = None
    if stress_before or stress_after:
        stress = stress_recheck(stress_before, stress_after, tol_kb=stress_tol_kb)

    # ---- 留档 ----
    if archive_dir is None:
        archive_dir = poscar.parent / "symmetry_audit"
    archive_dir = Path(archive_dir)
    archive_dir.mkdir(parents=True, exist_ok=True)
    before_path = archive_dir / "POSCAR.before"
    after_path = archive_dir / "POSCAR.symmetrized"
    write_poscar(before_path, st["lattice"], st["frac"], st["symbols"],
                 title="%s (BEFORE symmetrization)" % st["title"][:80])
    write_poscar(after_path, after["lattice"], after["frac"], after["symbols"],
                 title="%s (AFTER spglib.refine_cell symprec=%g)" % (st["title"][:80], used_symprec))
    write_poscar(poscar.with_name(poscar.name + ".pre_symmetry"),
                 st["lattice"], st["frac"], st["symbols"],
                 title="%s (pre-symmetry backup)" % st["title"][:80]) if inplace else None

    if out is None and not inplace:
        out = archive_dir / "POSCAR.symmetrized"
    if inplace:
        write_poscar(poscar, after["lattice"], after["frac"], after["symbols"],
                     title="%s (symmetrized symprec=%g)" % (st["title"][:80], used_symprec))
        out = poscar
    elif out:
        write_poscar(out, after["lattice"], after["frac"], after["symbols"],
                     title="%s (symmetrized symprec=%g)" % (st["title"][:80], used_symprec))

    result = {
        "triggered": True,
        "reason": "symprec=1e-5 与 1e-4 空间群不一致 -> 数值微畸变",
        "when": _now(),
        "source": str(poscar) if not source_note else source_note,
        "used_symprec": used_symprec,
        "audit_before": audit_before,
        "audit_after": audit_after,
        "confirmations": conf,
        "stress_recheck": stress,
        "files": {"before": str(before_path), "after": str(after_path),
                  "json": str(archive_dir / "symmetry_audit.json"),
                  "log": str(archive_dir / "symmetry_audit.log")},
    }
    if out:
        result["files"]["out"] = str(out)
    result["fit_for_use"] = bool(conf["fully_confirmed"]
                                 and (stress is None or stress["status"] == "pass"))

    (archive_dir / "symmetry_audit.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    with open(archive_dir / "symmetry_audit.log", "a", encoding="utf-8") as fh:
        fh.write("%s  %s  used_symprec=%g  c1=%s c2=%s c3=%s  fit_for_use=%s\n"
                 % (_now(), result["reason"], used_symprec,
                    conf["c1_stoichiometry"]["ok"], conf["c2_energy"]["status"],
                    conf["c3_vacuum_axis"]["ok"], result["fit_for_use"]))

    # provenance.json 合并（best-effort）
    if write_prov:
        try:
            sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
            import provenance_common as P
            P.record(step=archive_dir.parent.name or None, skill=None,
                     inputs=[str(poscar)], outputs=[result["files"]["after"]],
                     extra={"symmetry_audit": {
                         "triggered": True, "used_symprec": used_symprec,
                         "reason": result["reason"],
                         "micro_distortion_detail": audit_before["micro_distortion_detail"],
                         "confirmations": conf, "stress_recheck": stress}},
                     note="结构体检第五条：数值微畸变对称化", quiet=True, outdir=str(archive_dir))
        except Exception as exc:        # noqa: BLE001
            print("[symmetry_audit] provenance 合并跳过：%s" % exc, file=sys.stderr)

    # ---- 不静默通过 ----
    if not conf["all_ok"]:
        raise SystemExit("[symmetry_audit] 三项确认未过：%s"
                         % json.dumps(conf, ensure_ascii=False))
    if conf["c2_energy"]["status"] == "pending" and not allow_energy_pending:
        raise SystemExit("[symmetry_audit] 能量确认 pending —— 必须提供对称化前后两个单点，"
                         "或显式 --allow-energy-pending（结果 fit_for_use=false）")
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_parser():
    p = argparse.ArgumentParser(description="结构体检第五条：数值微畸变审计 + spglib 对称化")
    p.add_argument("--poscar", help="输入结构（POSCAR/CONTCAR）")
    p.add_argument("--mode", choices=["audit", "symmetrize", "energy-check"],
                   default="audit")
    p.add_argument("--out", help="对称化结构输出路径（默认 archive-dir/POSCAR.symmetrized）")
    p.add_argument("--archive-dir", help="留档目录（默认 <poscar 同目录>/symmetry_audit）")
    p.add_argument("--inplace", action="store_true",
                   help="把对称化结构写回 POSCAR（原结构备份 POSCAR.pre_symmetry）")
    p.add_argument("--symprec", type=float, default=None,
                   help="对称化容差（默认取两容差里识别出高对称的那个，通常 1e-4）")
    p.add_argument("--energy-before", help="对称化前单点能量：OUTCAR 路径或 eV 数值")
    p.add_argument("--energy-after", help="对称化后单点能量：OUTCAR 路径或 eV 数值")
    p.add_argument("--energy-tol-mev", type=float, default=None,
                   help="能量确认阈值 meV/atom（默认 1.0）")
    p.add_argument("--natoms", type=int, default=None, help="energy-check 的原子数")
    p.add_argument("--stress-before", help="对称化前单点 OUTCAR（面内应力复检）")
    p.add_argument("--stress-after", help="对称化后单点 OUTCAR（面内应力复检）")
    p.add_argument("--stress-tol-kb", type=float, default=0.2)
    p.add_argument("--allow-energy-pending", action="store_true",
                   help="能量确认缺单点时仍写出结果（fit_for_use=false）")
    p.add_argument("--allow-no-distortion", action="store_true",
                   help="空间群一致时仍执行对称化（默认跳过）")
    p.add_argument("--switch", choices=["off", "warn", "error"], default=None,
                   help="体检模式覆盖（warn 只告警；error 检出畸变即退出非 0）")
    p.add_argument("--strict-exit", action="store_true",
                   help="audit 模式检出畸变时退出码 2")
    p.add_argument("--json", action="store_true", help="stdout 输出 JSON")
    p.add_argument("--json-out", help="把结果 JSON 写到该路径")
    p.add_argument("--step-conf", help="step.conf 路径（默认 cwd/step.conf）")
    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)
    switches = read_switches(args.step_conf)
    switch = (args.switch or switches.get("SYMMETRY_AUDIT", "warn")).strip().lower()
    symprec = args.symprec if args.symprec is not None else float(
        switches.get("SYMMETRY_AUDIT_SYMPREC", DEFAULT_SYMPREC))
    etol = args.energy_tol_mev if args.energy_tol_mev is not None else float(
        switches.get("SYMMETRY_ENERGY_TOL_MEV", DEFAULT_ENERGY_TOL_MEV_PER_ATOM))

    if args.mode == "energy-check":
        res = energy_check(args.energy_before, args.energy_after,
                           natoms=args.natoms, tol_mev_per_atom=etol)
        _dump(res, args)
        return 0 if res["status"] == "pass" else 4

    if not args.poscar:
        print("[ERROR] --poscar 必填", file=sys.stderr)
        return 5
    if not spglib_available():
        print("[symmetry_audit] spglib 不可用：跳过（不改变现有行为）", file=sys.stderr)
        return 5

    try:
        if args.mode == "audit":
            res = audit_structure_file(args.poscar, symprec=symprec)
            print(format_audit_report(res))
            if args.json_out:
                Path(args.json_out).write_text(json.dumps(res, ensure_ascii=False, indent=2),
                                               encoding="utf-8")
            _dump(res, args) if args.json else None
            if res["micro_distortion"]:
                print("[WARN] 结构体检第五条：检出数值微畸变 "
                      "(symprec=1e-5 -> %s vs 1e-4 -> %s)，建议进 S4 前对称化"
                      % (res["spacegroup"][0]["international"],
                         res["spacegroup"][-1]["international"]))
                if switch == "error" or args.strict_exit:
                    return 2
            return 0
        # symmetrize
        res = symmetrize_structure_file(
            args.poscar, out=args.out, archive_dir=args.archive_dir,
            symprec=symprec, energy_before=args.energy_before,
            energy_after=args.energy_after, energy_tol_mev_per_atom=etol,
            inplace=args.inplace, allow_energy_pending=args.allow_energy_pending,
            stress_before=args.stress_before, stress_after=args.stress_after,
            stress_tol_kb=args.stress_tol_kb)
        print("[symmetry_audit] triggered=%s used_symprec=%g fit_for_use=%s"
              % (res.get("triggered"), res.get("used_symprec"), res.get("fit_for_use")))
        if res.get("triggered"):
            print(format_audit_report(res["audit_before"], "before"))
            print(format_audit_report(res["audit_after"], "after"))
            conf = res["confirmations"]
            print("  ① 原子数/计量比不变: %s  %s -> %s"
                  % (conf["c1_stoichiometry"]["ok"],
                     conf["c1_stoichiometry"]["before"], conf["c1_stoichiometry"]["after"]))
            print("  ② 单点能量差: %s  (%s)"
                  % (conf["c2_energy"]["status"], conf["c2_energy"].get("note")))
            print("  ③ 真空轴未移动: %s (axis %d -> %d, center_shift=%s)"
                  % (conf["c3_vacuum_axis"]["ok"],
                     conf["c3_vacuum_axis"]["before"]["axis"],
                     conf["c3_vacuum_axis"]["after"]["axis"],
                     conf["c3_vacuum_axis"]["center_shift_frac"]))
            if res.get("stress_recheck"):
                print("  ④ 面内残余应力复检: %s" % res["stress_recheck"]["status"])
        if args.json_out:
            Path(args.json_out).write_text(json.dumps(res, ensure_ascii=False, indent=2),
                                           encoding="utf-8")
        _dump(res, args) if args.json else None
        return 0
    except SystemExit:
        raise
    except Exception as exc:            # noqa: BLE001
        print("[ERROR] symmetry_audit 失败：%s" % exc, file=sys.stderr)
        return 5


def _dump(res, args):
    if args.json:
        print(json.dumps(res, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    raise SystemExit(main())
