#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fermi_window_check.py —— 费米能级窗口自检（2026-09-28 用户批准；含 scissor 修正）。

为什么要有这一步：
  掺杂浓度单位一旦写错（2D 的 DOPING 是**面密度 cm^-2**，填成体浓度 cm^-3 就错 8 个量级），
  AMSET 会去解一个根本不存在的费米能级。2026-09-27 MoS2 S8.4 那次解出 p 型 -153.9 eV
  （当时 DOS 窗口只有 [-6.88, 1.51] eV），过程里只留下一句 zero-size array，
  "离谱"这件事要人翻日志才看得出来。

本检查**不依赖 WRITE_MESH、也不依赖后处理**，直接读作业自己产出的两样东西：
  1. 运行目录里最新的 transport*.json 的 fermi_levels（eV，与 vasprun 同一参考零点）；
  2. 同目录 vasprun.xml 的 VBM / CBM；
  3. 窗口 = [VBM - win, CBM + **Δ** + win]，Δ 是 AMSET 对导带的整体上移量；
  4. 任一费米能级落在窗口外 -> 打印明细并 **退出码 1**（作业链失败）；
  5. 金属 / 零带隙 / 取不到带边 -> 打印跳过，退出码 0（不误杀金属体系）。

★ Δ 为什么必须算（2026-09-28 用户指正）：S8 默认 REQUIRE_BANDGAP=True，一律写 bandgap:，
  AMSET 把导带整体上移 Δ = bandgap - DFT带隙（源码 amset/interpolation/bandstructure.py:716
  `scissor = bandgap - interp_bandgap`）。n 型费米能级因此跑到 CBM_DFT + Δ 附近，
  简并掺杂还会再进导带若干 —— 不算 Δ 会把**跑完好几个小时的正常作业误判为失败**。
  GaAs Δ≈0.86 eV、GaN Δ≈1.5 eV，正好是这个量级。

★ 带边怎么取（2026-09-28 实测修正）：**不能**用 pymatgen 的 bs.get_vbm()/get_cbm() ——
  它们要先 is_metal()，而 is_metal() 依赖 bs.efermi；实测 GaAs/GaN 的 step3 vasprun.xml
  **没有 <efermi> 标签**，efermi=None -> TypeError。改用本征值第 2 列（占据数）：
  VBM = 占据态最高能、CBM = 空态最低能。实测 GaAs 0.418 eV、GaN 1.901 eV，符合 PBE。

用法（作业链里 cwd 就是运行目录）：
    python fermi_window_check.py                 # win=2.0 eV，Δ 自动从 settings.yaml 算
    python fermi_window_check.py --win 3.0 --settings settings.yaml
"""
import argparse
import json
import sys
from pathlib import Path

WIN_DEFAULT = 2.0


def newest_transport(run_dir):
    """运行目录里最新的 transport*.json（排除 intrinsic_transport.json）。"""
    cands = [p for p in Path(run_dir).glob("transport*.json")
             if not p.name.startswith("intrinsic")]
    return max(cands, key=lambda p: p.stat().st_mtime) if cands else None


def load_fermi_levels(path):
    """读 transport json 的 fermi_levels，摊平成 [((i, j), eV), ...]。"""
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    out = []

    def walk(x, idx):
        if isinstance(x, (list, tuple)):
            for k, y in enumerate(x):
                walk(y, idx + [k])
        elif x is not None:
            out.append((tuple(idx), float(x)))

    walk(d.get("fermi_levels"), [])
    return out


def scissor_delta(settings_path, dft_gap):
    """AMSET 对导带的整体上移量 Δ（eV）。

    优先级与 AMSET 一致（amset/interpolation/bandstructure.py:711-716）：
      写了 bandgap  ->  Δ = bandgap - DFT带隙；
      只写了 scissor ->  Δ = scissor；
      都没有        ->  0。
    读不到 settings.yaml 时返回 0（退化成旧行为，不拦错）。
    """
    p = Path(settings_path)
    if not p.is_file():
        return 0.0
    try:
        import yaml
        s = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    except Exception:                                          # noqa: BLE001
        return 0.0
    bg = s.get("bandgap")
    if bg:
        return max(0.0, float(bg) - float(dft_gap))
    sc = s.get("scissor")
    if sc:
        return max(0.0, float(sc))
    return 0.0


def check_window(pairs, vbm, cbm, win=WIN_DEFAULT, delta=0.0):
    """纯函数：返回 (是否全在窗口内, 下限, 上限, 越界的 [(索引, 值), ...])。

    窗口 = [VBM - win, CBM + delta + win]（价带不上移，导带上移 Δ）。
    """
    lo, hi = vbm - win, cbm + delta + win
    bad = [(ix, f) for ix, f in pairs if not (lo <= f <= hi)]
    return (len(bad) == 0), lo, hi, bad


def band_edges(vasprun):
    """从 vasprun.xml 取 (VBM, CBM)（eV）；金属/取不到返回 (None, None)。

    用本征值的第 2 列（占据数）：VBM = 占据态最高能，CBM = 空态最低能。
    这样不依赖 <efermi>（实测这些 vasprun 里根本没有该标签）。
    """
    import numpy as np
    from pymatgen.io.vasp.outputs import Vasprun
    v = Vasprun(str(vasprun), parse_dos=False)
    vb, cb = -np.inf, np.inf
    any_v = any_c = False
    for _spin, ev in v.eigenvalues.items():
        a = np.asarray(ev, float)
        if a.ndim != 3 or a.shape[-1] < 2:
            raise ValueError("vasprun 本征值里没有占据数列（shape=%s）" % (a.shape,))
        e, occ = a[:, :, 0], a[:, :, 1]
        if np.any(occ > 0.5):
            vb = max(vb, float(e[occ > 0.5].max())); any_v = True
        if np.any(occ < 0.5):
            cb = min(cb, float(e[occ < 0.5].min())); any_c = True
    if not (any_v and any_c) or cb <= vb:
        return None, None
    return vb, cb


def _fmt_idx(ix):
    if not ix:
        return "?"
    parts = ["doping#%d" % ix[0]]
    if len(ix) > 1:
        parts.append("T#%d" % ix[1])
    return " ".join(parts)


def main(argv=None):
    ap = argparse.ArgumentParser(description="费米能级窗口自检")
    ap.add_argument("--win", type=float, default=WIN_DEFAULT,
                    help="允许超出带边的窗口 (eV)，默认 2.0")
    ap.add_argument("--file", default=None, help="指定 transport json（默认取最新）")
    ap.add_argument("--vasprun", default="vasprun.xml", help="vasprun.xml 路径")
    ap.add_argument("--settings", default="settings.yaml", help="settings.yaml 路径（取 bandgap/scissor）")
    ap.add_argument("--delta", type=float, default=None, help="手动指定 Δ（默认从 settings 算）")
    args = ap.parse_args(argv)

    run_dir = Path.cwd()
    tf = Path(args.file) if args.file else newest_transport(run_dir)
    if tf is None or not tf.is_file():
        print("[WARN] 运行目录里没有 transport*.json —— 跳过费米能级窗口检查")
        return 0
    try:
        pairs = load_fermi_levels(tf)
    except Exception as e:
        print("[WARN] 读 %s 失败（%s: %s）—— 跳过" % (tf, type(e).__name__, e))
        return 0
    if not pairs:
        print("[..] %s 里没有 fermi_levels（金属/本征？）—— 跳过" % tf.name)
        return 0
    vp = Path(args.vasprun)
    if not vp.is_file():
        print("[WARN] 找不到 %s —— 跳过费米能级窗口检查" % vp)
        return 0
    try:
        vbm, cbm = band_edges(vp)
    except Exception as e:
        print("[WARN] 读 %s 的带边失败（%s: %s）—— 跳过" % (vp, type(e).__name__, e))
        return 0
    if vbm is None or cbm is None or cbm <= vbm:
        print("[..] 金属/零带隙（VBM=%s CBM=%s）—— 跳过费米能级窗口检查" % (vbm, cbm))
        return 0
    delta = args.delta if args.delta is not None else scissor_delta(args.settings, cbm - vbm)
    ok, lo, hi, bad = check_window(pairs, vbm, cbm, args.win, delta)
    print("[..] 费米能级窗口：VBM=%.4f CBM=%.4f DFT带隙=%.4f Δ(scissor)=%.4f "
          "-> 允许 [%.4f, %.4f] eV；%d 个费米能级（%s）"
          % (vbm, cbm, cbm - vbm, delta, lo, hi, len(pairs), tf.name))
    if ok:
        print("[OK] 全部费米能级都在带边窗口内")
        return 0
    print("[ERROR] %d 个费米能级落在带边窗口之外：" % len(bad))
    for ix, f in bad[:8]:
        print("        %s: %.4f eV（超出 %.4f）" % (_fmt_idx(ix), f, lo if f < lo else hi))
    if len(bad) > 8:
        print("        ... 其余 %d 个省略" % (len(bad) - 8))
    print("        常见根因：① 掺杂浓度单位写错（2D 的 DOPING 是面密度 cm^-2）；"
          "② vasprun 与 h5 不同源；③ energy_cutoff/带隙设置把带边切掉了；"
          "④ settings.yaml 的 bandgap 写错（Δ 会跟着错）。"
          "        这种情况下的迁移率/Seebeck 不可用，请勿采信。")
    return 1


if __name__ == "__main__":
    sys.exit(main())
