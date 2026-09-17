#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""skill/_common/thickness_2d.py —— 2D 层厚与 κ/σ 归一化因子（唯一真源）。

背景
----
phono3py / ShengBTE / FourPhonon / AMSET 都用【含真空的晶胞体积】做分母，
所以 2D 材料算出来的面内 κ（以及 σ）都被胞高稀释了。要还原成"层本身的量"
必须乘 h⊥/d：

    κ_inplane(层) = κ_raw(面内分量) × h⊥ / d

h⊥ = 周期性胞高，d = 层有效厚度。**kl 与 ke 两条链必须用同一个 d**，否则 zT
里 σ 与 κ_L 的厚度口径错开，d 约不掉。本模块就是那个唯一来源：kl-dft-cpu 的
S6_kappa 与 ke-dft-cpu 的 AMSET 步都读它（ke 侧接入前，kl 侧已按同一口径出数）。

三个量的定义
------------
h⊥ = V / |a_i × a_j|
    周期性胞高（体积 / 面内面积），**不是 |c|**。phono3py 与 ShengBTE 的分母
    其实是 V 与 A，所以归一化因子必须用 V/A；c 轴一旦倾斜（α/β≠90°），
    |c| ≠ V/A，用 |c| 就把因子算错了。

z  = 以"最大间隙（真空）"为切口把分数坐标展开后再乘 h⊥
    直接对 frac∈[0,1) 取 max−min 会出错：层若跨越周期边界（一半原子 z≈0.02、
    另一半 z≈0.98），算出的 span 接近整个胞高 → d > h⊥ → 因子小于 1。
    做法与 dim_common.vacuum_per_axis 的最大间隙法一致。

d  = atomic_span + r_vdW(顶层原子) + r_vdW(底层原子)
    这是文献里常见的一种取法（起伏高度 + 两侧 vdW 半径）。d 的取法本身在
    文献中不统一（层间距+vdW 修正 / 统一厚度 / 直接报面热导 W/K），关键不是
    取值而是两侧同源；要换口径就把 KAPPA_2D_THICKNESS 设成固定数值 Å。

只依赖 numpy + 同目录 vdw_radii.py（公共池保证两者一起推送；缺表时全部元素
退到 VDW_FALLBACK 并告警，不静默换口径）。
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from vdw_radii import VDW_RADII, VDW_FALLBACK
except ImportError:                                     # pragma: no cover
    print("[WARN] thickness_2d 找不到 vdw_radii.py：所有元素退到 r_vdW=2.0 Å，"
          "层厚口径不可比（请检查 gen_need 是否漏推 vdw_radii.py）")
    VDW_RADII = {}
    VDW_FALLBACK = 2.00

_CELL_MODES = ("cell", "lz", "none", "")


def slab_geometry(lat, frac, species=None, vac_axis=2, mode="vdw", fallback_r=None):
    """2D 层几何 + κ 归一化因子。

    参数
    ----
    lat      : 3x3 晶格（Å，行=a 向量）
    frac     : Nx3 分数坐标
    species  : 长度 N 的元素符号列表；None / 长度不符 → 退回不归一（并说明）
    vac_axis : 真空轴（0/1/2）
    mode     : "vdw"（默认）| "cell"/"lz"/"none"（d=h⊥，即不归一）| 数值字符串（固定 Å）
    fallback_r: 表里查不到的元素半径，默认 VDW_FALLBACK

    返回 dict（键名与 kappa_summary.json 里落地的字段一致）
    """
    lat = np.asarray(lat, float)
    frac = np.asarray(frac, float)
    ax = int(vac_axis if vac_axis is not None else 2)
    i, j = [k for k in range(3) if k != ax]

    area = float(np.linalg.norm(np.cross(lat[i], lat[j])))
    vol = abs(float(np.linalg.det(lat)))
    h = vol / area if area > 1e-12 else 0.0          # 周期性胞高 V/A
    lz = float(np.linalg.norm(lat[ax]))              # |c|（仅作对照，不参与归一化）

    n = int(len(frac))
    if n:
        f = np.mod(frac[:, ax], 1.0)
        fs = np.sort(f)
        # 最大间隙 = 真空；切口取在间隙之后，层就被连续展开（不跨越周期边界）
        gaps = np.diff(np.r_[fs, fs[0] + 1.0])
        k = int(np.argmax(gaps))
        start = float(fs[(k + 1) % n])
        z = np.mod(f - start, 1.0) * h
        itop, ibot = int(np.argmax(z)), int(np.argmin(z))
        span = float(z[itop] - z[ibot])
    else:
        itop = ibot = 0
        span = 0.0

    fb = VDW_FALLBACK if fallback_r is None else float(fallback_r)
    sp_ok = species is not None and len(species) == n
    m = str(mode).strip().lower()
    try:
        d = float(m)
        conv = "fixed %.3f A" % d
        top = bot = None
    except ValueError:
        if m in _CELL_MODES:
            d = h
            conv = "cell h_perp (no norm)"
            top = bot = None
        elif not sp_ok:
            d = h
            conv = "cell h_perp (no species -> fallback)"
            top = bot = None
        else:
            top, bot = str(species[itop]), str(species[ibot])
            r_top = float(VDW_RADII.get(top, fb))
            r_bot = float(VDW_RADII.get(bot, fb))
            d = span + r_top + r_bot
            conv = "span %.2f + vdW(%s %.2f + %s %.2f)" % (span, top, r_top, bot, r_bot)

    factor = (h / d) if d else 1.0
    gap = h - span if h else 0.0
    return {
        "h_perp_A": round(h, 4),
        "Lz_A": round(lz, 4),
        "atomic_span_A": round(span, 4),
        "vacuum_gap_A": round(gap, 4),
        "thickness_d_A": round(d, 4) if d else None,
        "kappa_2d_norm_factor": round(factor, 6),
        "thickness_convention": conv,
        "top_species": top,
        "bottom_species": bot,
    }


def write_material_level(matdir, meta, tol=0.05):
    """把层厚口径写到【材料根目录】的 thickness_2d.json —— kl / ke 共用契约。

    规则（文档 P0-2(b)）：谁先跑谁写；后跑的一方读它、并用自己弛豫好的结构重算
    对照，差值 > tol（默认 0.05 Å）就告警。
      - 目标不存在 / 已存在但来源就是自己（source 以 "kl-dft-cpu" 开头）→ 覆盖写。
      - 已存在且来自别的来源（如 ke 的 AMSET）→ 只校验不回写，冲突则告警，
        保持"先写入者生效"，避免两条链互相覆盖。
    返回一行可打印的状态说明；任何异常由调用方吞掉（这是 best-effort 契约，不该
    让作业失败）。
    """
    import json
    matdir = Path(matdir)
    p = matdir / "thickness_2d.json"
    src = str(meta.get("source", ""))
    old = None
    if p.is_file():
        try:
            old = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            old = None
    # 先写入者生效：已存在且来源【不是自己】时只校验不回写（避免两条链互相覆盖）。
    same_writer = old is not None and str(old.get("source", "")) == src
    if old is not None and not same_writer:
        d_old, d_new = old.get("thickness_d_A"), meta.get("thickness_d_A")
        if d_old and d_new and abs(float(d_old) - float(d_new)) > float(tol):
            return ("[WARN] 材料级 thickness_2d.json 已由 %s 写入 d=%.3f Å，本次用自己弛豫的结构算得 "
                    "d=%.3f Å（差 %.3f Å > %.2f）——两条链口径可能不一致，zT 里 d 约不掉，"
                    "请核对；本次不覆盖。"
                    % (old.get("source"), float(d_old), float(d_new),
                       abs(float(d_old) - float(d_new)), float(tol)))
        return ("[..] 材料级 thickness_2d.json 已存在（%s，d=%.3f Å），与本次计算一致，"
                "不覆盖。" % (old.get("source"), float(old.get("thickness_d_A") or 0.0)))
    p.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    _what = "覆盖更新（同一来源重跑）" if same_writer else "新建"
    return "[OK] 材料级 thickness_2d.json %s %s（source=%s, d=%.3f Å）" % (
        _what,
        p, meta.get("source"), float(meta.get("thickness_d_A") or 0.0))


def slab_geometry_from_poscar(poscar, vac_axis=2, mode="vdw", fallback_r=None,
                              reader=None):
    """读 POSCAR（或 CONTCAR）算层几何。reader 缺省用 dim_common.read_poscar_cell_frac
    （公共池随各技能一起推送）；没有 dim_common 时退到 ASE。"""
    if reader is None:
        try:
            from dim_common import read_poscar_cell_frac as reader
        except ImportError:
            from ase.io import read as _ase_read

            def reader(path):
                at = _ase_read(str(path), format="vasp")
                return at.cell[:], at.get_scaled_positions()

    lat, frac = reader(poscar)
    species = _species_from_poscar(poscar, len(frac))
    return slab_geometry(lat, frac, species, vac_axis=vac_axis, mode=mode,
                         fallback_r=fallback_r)


def _species_from_poscar(poscar, n):
    """从 POSCAR 第 6/7 行读元素符号展开成逐原子列表；VASP4（无元素行）返回 None。"""
    import re
    try:
        lines = Path(poscar).read_text(encoding="utf-8-sig").splitlines()
        line6 = lines[5].split()
        if not line6 or re.fullmatch(r"[+-]?\d+", line6[0]):
            return None
        counts = [int(x) for x in lines[6].split()]
        out = []
        for sym, c in zip(line6, counts):
            out += [sym] * c
        return out if len(out) == n else None
    except Exception:
        return None
