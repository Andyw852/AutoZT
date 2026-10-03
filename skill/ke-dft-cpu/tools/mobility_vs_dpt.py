#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mobility_vs_dpt.py —— 按载流子列出 AMSET 分机制迁移率并与 DPT 对照（V143；只读，秒级，不需要 AMSET）。

**读 transport.json 的迁移率一律用本工具，不要手写脚本。**
起因（2026-10-02）：手写脚本把掺杂符号对调，MoSe2 的电子/空穴 ADP/DPT 整张表标反，来回核对了两轮。
AMSET 的约定：负掺杂 = n 型（电子），正掺杂 = p 型（空穴）—— 唯一真源是 ke_common.carrier_of_doping。
本工具另用 Seebeck 的符号自检（n 型 S < 0、p 型 S > 0）：对不上就标 ★，说明载流子判反了或数据有问题。

输出（指定温度，默认 300 K）：每个掺杂一行 —— 载流子、Seebeck、overall 与各机制迁移率（2D 取面内 xx/yy 平均，
3D 取迹/3），以及 DPT（step8.2_dpt/dpt_result.json 同一载流子）和 ADP/DPT。
ADP/DPT 超出 [1/3, 3]：标 ⚠，提示用 tools/dp_valley_probe.py 查是不是多谷（V141/V142）。
  例外（只注明、不提示多谷）：n/N_eff ≥ 0.5（接近/进入简并，DPT 的非简并统计不成立，V146）；
  DPT 的 m* 取自简并带边（单带公式只是近似，V145）；DPT 的 m* 网格分辨不出带边曲率（V147）。

用法（在材料目录）：
    python <skill>/ke-dft-cpu/tools/mobility_vs_dpt.py [--run step8.4_amset2d] [--T 300] [--json out.json]
    --run 默认：有 step8.4_amset2d/transport.json 用它，否则 step8_amset。
退出码：0 正常；1 有 ★（Seebeck 符号与载流子不符）；2 输入问题。
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in (str(HERE.parent), str(HERE.parent.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402

LABEL = {"electron": "电子(n)", "hole": "空穴(p)"}
RATIO_OK = (1.0 / 3.0, 3.0)
# V146：DPT 用的是非简并（Boltzmann）统计。载流子浓度接近有效态密度时（n/N_eff ≳ 0.5）它就不成立，
#   ADP/DPT 超界不能拿来判多谷。N_eff 用 DPT 自己的 m*（单谷、含自旋）估：
#   3D Nc = 2.509e19·m*^1.5·(T/300)^1.5 cm⁻³；2D N = 1.0799e13·m*·(T/300) cm⁻²（体浓度 × c 换回面密度）。
#   GaAs 电子 m* 0.052 -> Nc ≈ 3e17：S8 的 3D 掺杂 1e18–1e21 全在简并区。
DEGEN_RATIO = 0.5


def _avg(t, is_2d):
    t = np.asarray(t, float)
    return float((t[0, 0] + t[1, 1]) / 2.0) if is_2d else float(np.trace(t) / 3.0)


def dpt_by_carrier(path):
    """dpt_result.json -> {"electron": μ, "hole": μ}（没有 / 算不出为 None）。"""
    return _dpt(path)[0]


def dpt_degenerate(path):
    """dpt_result.json 里 m* 取自简并带边的载流子（V145：m_provenance 含"重简并"）。"""
    return _dpt(path)[1]


def _dpt(path):
    out, deg = {"electron": None, "hole": None}, set()
    p = Path(path)
    if not p.is_file():
        return out, deg
    for r in json.loads(p.read_text(encoding="utf-8")).get("results") or []:
        if r.get("carrier") in out:
            out[r["carrier"]] = r.get("mobility_cm2_Vs")
            if "重简并" in str((r.get("inputs") or {}).get("m_provenance", "")):
                deg.add(r["carrier"])
    return out, deg


def dpt_grid_unresolved(path):
    """V147：DPT 的 m* 网格分辨不出带边曲率的载流子（m_provenance 含"网格分辨不出"）。"""
    out = set()
    p = Path(path)
    if p.is_file():
        for r in json.loads(p.read_text(encoding="utf-8")).get("results") or []:
            if "网格分辨不出" in str((r.get("inputs") or {}).get("m_provenance", "")):
                out.add(r.get("carrier"))
    return out


def dpt_masses(path):
    """dpt_result.json -> {"electron": m*, "hole": m*}（m0 为单位；没有为 None）。"""
    out = {"electron": None, "hole": None}
    p = Path(path)
    if p.is_file():
        for r in json.loads(p.read_text(encoding="utf-8")).get("results") or []:
            if r.get("carrier") in out:
                out[r["carrier"]] = (r.get("inputs") or {}).get("m_eff_m0")
    return out


def degeneracy_ratio(doping, m_eff, T, is_2d, c_A=None):
    """n / N_eff（单谷、含自旋）。拿不到 m* 或（2D）c 时返回 None。"""
    if not m_eff or m_eff <= 0:
        return None
    n = abs(float(doping))
    if is_2d:
        if not c_A:
            return None
        return n * float(c_A) * 1e-8 / (1.0799e13 * float(m_eff) * (T / 300.0))
    return n / (2.509e19 * float(m_eff) ** 1.5 * (T / 300.0) ** 1.5)


def rows(transport, T=300.0, is_2d=True, dpt=None, masses=None, c_A=None):
    d = transport
    dop = np.asarray(d["doping"], float)
    temps = np.asarray(d["temperatures"], float)
    it = int(np.argmin(np.abs(temps - T)))
    mob = d.get("mobility") or {}
    mechs = [m for m in mob if m != "overall" and mob[m] is not None]
    S = np.asarray(d["seebeck"], float)
    dpt = dpt or {}
    out = []
    for i, x in enumerate(dop):
        car = kc.carrier_of_doping(x)
        s = _avg(S[i, it], is_2d)
        r = {"doping": float(x), "carrier": car, "T": float(temps[it]), "seebeck_uV_K": round(s, 2)}
        r["sign_ok"] = (car is None) or (car == "electron" and s < 0) or (car == "hole" and s > 0)
        if mob.get("overall") is not None:
            r["overall"] = round(_avg(np.asarray(mob["overall"], float)[i, it], is_2d), 2)
        for m in mechs:
            r[m] = round(_avg(np.asarray(mob[m], float)[i, it], is_2d), 2)
        mu_dpt = dpt.get(car) if car else None
        if car and masses:
            q_deg = degeneracy_ratio(x, masses.get(car), r["T"], is_2d, c_A)
            if q_deg is not None:
                r["n_over_Neff"] = round(q_deg, 3)
        r["DPT"] = mu_dpt
        if mu_dpt and r.get("ADP"):
            r["ADP_over_DPT"] = round(r["ADP"] / mu_dpt, 3)
        out.append(r)
    return out, mechs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("path", nargs="?", default=".", help="材料目录（默认当前目录）")
    ap.add_argument("--run", default=None, help="step8.4_amset2d / step8_amset（默认自动）")
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    mat = Path(a.path)
    run = a.run or ("step8.4_amset2d" if (mat / "step8.4_amset2d" / "transport.json").is_file() else "step8_amset")
    tj = mat / run / "transport.json"
    if not tj.is_file():
        print("[ERROR] 找不到 %s" % tj, file=sys.stderr)
        return 2
    is_2d = run.startswith("step8.4") or (mat / run / "2d_correction.json").is_file()
    dpt = dpt_by_carrier(mat / "step8.2_dpt" / "dpt_result.json")
    deg = dpt_degenerate(mat / "step8.2_dpt" / "dpt_result.json")
    masses = dpt_masses(mat / "step8.2_dpt" / "dpt_result.json")
    coarse = dpt_grid_unresolved(mat / "step8.2_dpt" / "dpt_result.json")
    c_A = None
    if is_2d and (mat / run / "2d_correction.json").is_file():
        try:
            c_A = json.loads((mat / run / "2d_correction.json").read_text(encoding="utf-8")).get("cell_c_A")
        except ValueError:
            c_A = None
    try:
        rs, mechs = rows(json.loads(tj.read_text(encoding="utf-8")), a.T, is_2d, dpt, masses, c_A)
    except (KeyError, ValueError, IndexError) as e:
        print("[ERROR] 读 %s 失败：%s: %s" % (tj, type(e).__name__, e), file=sys.stderr)
        return 2
    print("%s（%s，%g K，迁移率 cm²/Vs，%s；载流子按 AMSET 约定：负掺杂 = 电子）"
          % (tj, "2D 面内 xx/yy 平均" if is_2d else "3D 迹/3", rs[0]["T"] if rs else a.T,
             "DPT 来自 step8.2_dpt" if any(dpt.values()) else "没有 DPT 结果"))
    cols = ["overall"] + mechs
    print("  %-12s %-8s %9s " % ("doping/cm⁻³", "载流子", "S/μV·K⁻¹") + " ".join("%9s" % c for c in cols)
          + " %9s %9s" % ("DPT", "ADP/DPT"))
    bad = 0
    for r in rs:
        flag = ""
        if not r["sign_ok"]:
            flag = "  ★ Seebeck 符号与载流子不符"
            bad += 1
        q = r.get("ADP_over_DPT")
        if q is not None and not (RATIO_OK[0] <= q <= RATIO_OK[1]):
            nq = r.get("n_over_Neff")
            if nq is not None and nq >= DEGEN_RATIO:
                flag += ("  （n/N_eff = %.2g：接近/进入简并，DPT 的非简并统计不成立，ADP/DPT 不作判据）" % nq)
            elif r["carrier"] in deg:
                flag += "  （DPT 的带边简并，单带公式只是近似，ADP/DPT 不作判据）"
            elif r["carrier"] in coarse:
                flag += ("  （DPT 的 m* 网格分辨不出带边曲率、偏重，ADP/DPT 不作判据；"
                         "要比就用 amset eff-mass 的质量填 MANUAL 重跑 S8.2）")
            else:
                flag += "  ⚠ ADP/DPT 超出 [1/3, 3]：用 tools/dp_valley_probe.py 查是否多谷"
        print("  %-12.4g %-8s %9.1f " % (r["doping"], LABEL.get(r["carrier"], "本征"), r["seebeck_uV_K"])
              + " ".join("%9s" % (r.get(c, "")) for c in cols)
              + " %9s %9s" % (r["DPT"] if r["DPT"] is not None else "-", q if q is not None else "-") + flag)
    if a.json:
        Path(a.json).write_text(json.dumps(rs, ensure_ascii=False, indent=1), encoding="utf-8")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
