#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step12c_odp_read.py —— [0035] 读冻声子单点，算谷内光学形变势（step8.21_odp_read，S8.21_odp）。

run: gen 步骤（登录节点直接跑，秒级到分钟级）。输入：
  step7c_odp_frozen/{odp-ref, odp-mNNp, odp-mNNm}/vasprun.xml + LOCPOT + odp_manifest.json（S7.2）
  step8.2_dpt/dpt_result.json 的多谷谷表（每个谷的自旋、带号、各成员谷在 S7 网格上的 k_S7）
产出：step8.21_odp_read/odp_params.json（done_marker）+ odp_summary.txt。S8.4 的 ODP_FROZEN=auto 读它。

每个谷 × 每组简并光学模：D = d(E − V_vac)/du（eV/Å），各成员谷平均（odp_frozen.py 有推导与归一化）。
核对（不过就报错，绝不静默出数）：
  · S7.2 照抄的 S7 undeformed 输入（INCAR 物理键、KPOINTS、POTCAR）与 S7 现在的一致；
  · 参考构型 odp-ref 在各谷 k 点的本征值与 S7 undeformed 一致（< 5 meV）—— 带号/自旋/网格都对上了；
  · 每个构型的真空区平面平均势平坦（std < 1 meV）；
  · S5 现在的频率与生成位移时的一致（差 > 2% 告警）。
"""
import json
import os
import sys
import time as _time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, os.getcwd())
import ke_common as kc  # noqa: E402

of = None                     # odp_frozen：main() 里切好 python 环境后再 import（它要 numpy）

OUTDIR_NAME = "step8.21_odp_read"
SRC_DIR = "step7c_odp_frozen"
DEFORM_DIR = "step7_deform"
DIELECT_DIR = "step5_dielect"
DPT_RESULT = "step8.2_dpt/dpt_result.json"
MANIFEST = "odp_manifest.json"
REF_NAME = "odp-ref"
_SKILL_REV = "2026-10-11-0035-r1"
REF_TOL_EV = 0.005            # odp-ref 与 S7 undeformed 同一 k、同一带的本征值差上限
SNAP_TOL_INVA = 0.02          # 旧 S8.2（没有 k_S7）时把 k_grid 吸附到本步网格的容差（Å⁻¹）
SHIFT_WARN_EV = 0.1           # 带边移动超过它：出线性区，建议减小 ODP_DISP_A
FREQ_WARN = 0.02              # S5 现在的频率与生成时相差超过它 -> 告警
FREQ_FAIL = 0.05              # ……超过它 -> 报错（位移用的是旧本征矢，要 retry S7.2）
VAC_CACHE = "vacuum_levels.json"   # 各构型真空能级缓存（按 vasprun.xml 的 sha256 认）：LOCPOT 清掉后 S8.21 重读照样能跑
OUTCAR_DONE = "General timing and accounting informations"


def _die(msg):
    sys.exit("[ERROR] " + msg)


def _vasprun(d):
    for n in ("vasprun.xml", "vasprun.xml.gz"):
        if (d / n).is_file():
            return d / n
    return None


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Config(object):
    """一个单点构型：本征值（按自旋）、k 列表、真空能级。"""

    def __init__(self, d, vac_axis=2, cache=None):
        import numpy as np
        from pymatgen.io.vasp import Locpot, Vasprun
        self.name = d.name
        oc = d / "OUTCAR"
        if not oc.is_file() or OUTCAR_DONE not in oc.read_text(errors="ignore")[-400000:]:
            _die("%s/%s 的 OUTCAR 没跑完（S7.2 还有单点没算完？）" % (SRC_DIR, d.name))
        vr = _vasprun(d)
        if vr is None:
            _die("缺 %s/%s/vasprun.xml" % (SRC_DIR, d.name))
        v = Vasprun(str(vr), parse_dos=False, parse_potcar_file=False)
        if not v.converged_electronic:
            _die("%s/%s 的电子步没收敛（跑到 NELM 了）：本征值不可信，D 会是噪声。查 OSZICAR，调 NELM/ALGO 后 retry S7.2"
                 % (SRC_DIR, d.name))
        self.kpts = np.asarray(v.actual_kpoints, float)
        self.eig = {}
        for sp, arr in v.eigenvalues.items():
            self.eig[str(getattr(sp, "name", sp)).lower()] = np.asarray(arr)[:, :, 0]
        self.structure = v.final_structure
        lp = d / "LOCPOT"
        key = _sha256(vr)
        ent = (cache or {}).get(d.name) or {}
        if lp.is_file():
            planar = np.asarray(Locpot.from_file(str(lp)).get_average_along_axis(vac_axis))
            zf = np.asarray(self.structure.frac_coords)[:, vac_axis]
            try:
                self.vac, self.vac_std = of.vacuum_level(planar, zf)
            except ValueError as e:
                _die("%s/%s：%s" % (SRC_DIR, d.name, e))
            if cache is not None:
                cache[d.name] = {"vasprun_sha256": key, "vac_eV": self.vac, "std_eV": self.vac_std}
        elif ent.get("vasprun_sha256") == key:
            self.vac, self.vac_std = float(ent["vac_eV"]), float(ent["std_eV"])
            print("[..] %s：LOCPOT 已清，用上次读出的真空能级（%s，vasprun.xml 未变）" % (d.name, VAC_CACHE))
        else:
            _die("缺 %s/%s/LOCPOT，也没有这次运行的真空能级缓存（%s/%s）：真空对齐要它。INCAR 要 LVHAR=.TRUE.；"
                 "第一次跑 S8.21 之前别清 LOCPOT（跑过一次之后可以清）" % (SRC_DIR, d.name, OUTDIR_NAME, VAC_CACHE))

    def kidx(self, kf):
        import numpy as np
        for sg in (1.0, -1.0):
            dk = self.kpts - sg * np.asarray(kf, float)
            dk -= np.round(dk)
            dd = np.linalg.norm(dk, axis=1)
            i = int(np.argmin(dd))
            if dd[i] <= 1e-4:
                return i
        return None

    def bands(self, kf, spin):
        i = self.kidx(kf)
        if i is None:
            _die("k=%s 不在 %s 的 k 网格（含 −k）上 —— S7.2 与 S7 的网格不同？retry S7.2"
                 % ([round(float(x), 5) for x in kf], self.name))
        if spin not in self.eig:
            _die("%s 没有自旋道 %s（ISPIN 与 S3/S7 不同？）" % (self.name, spin))
        return self.eig[spin][i]


def check_inputs(cwd, man):
    """S7.2 生成时照抄的 S7 undeformed 输入，与 S7 现在的输入比（INCAR 只比物理键）。"""
    und = cwd / DEFORM_DIR / "undeformed"
    ref = cwd / SRC_DIR / REF_NAME
    bad = []
    if not (ref / "INCAR").is_file():
        _die("缺 %s/%s/INCAR（S7.2 没生成完整？retry S7.2_odp）" % (SRC_DIR, REF_NAME))
    for n in ("KPOINTS", "POTCAR"):
        a, b = und / n, ref / n
        if not a.is_file():
            bad.append("%s/undeformed/%s 不在了" % (DEFORM_DIR, n))
        elif not b.is_file() or kc._input_key(n, a.read_text(errors="ignore")) != kc._input_key(n, b.read_text(errors="ignore")):
            bad.append(n)
    if (und / "INCAR").is_file():
        ka = dict(kc._input_key("INCAR", (und / "INCAR").read_text(errors="ignore")))
        kb = dict(kc._input_key("INCAR", (ref / "INCAR").read_text(errors="ignore")))
        diff = sorted(k for k in set(ka) | set(kb) if ka.get(k) != kb.get(k))
        if diff:
            bad.append("INCAR（%s）" % ", ".join("%s: %s -> %s" % (k, kb.get(k), ka.get(k)) for k in diff[:6]))
    if bad:
        _die("S7.2 的单点用的输入与 S7 undeformed 现在的不同：%s。S7 重新生成过 -> 先 retry S7.2_odp（重算单点），"
             "再重跑本步。（step.conf 给 S7.2 单独改了 INCAR 也会这样：S7.2 必须与 S7 同一套参数）" % "；".join(bad))


def check_s5(cwd, man):
    """S5 现在的频率 vs 生成位移时记下的：差 > FREQ_WARN 告警，> FREQ_FAIL 或模数不同报错（位移用的是旧本征矢）。
    S5 读不了（被清理）只告警。"""
    import numpy as np
    try:
        from pymatgen.io.vasp import Vasprun
        vr = _vasprun(cwd / DIELECT_DIR)
        v = Vasprun(str(vr), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
        ev = -np.asarray(v.normalmode_eigenvals)[::-1]
        f = np.sqrt(np.abs(ev)) * np.sign(ev)
        if int(str(v.vasp_version).split(".")[0]) < 6:
            f = f * 15.633302
        old = np.asarray(man.get("s5_freq_THz") or [], float)
        f = np.sort(f)
        if len(old) != len(f):
            _die("S5 的模数与生成 S7.2 时不同（%d vs %d）：结构变了，retry S7.2_odp" % (len(f), len(old)))
        opt = np.abs(old) > 1.0
        rel = float(np.max(np.abs(f[opt] - old[opt]) / np.abs(old[opt]))) if opt.any() else 0.0
        if rel > FREQ_FAIL:
            _die("S5 现在的光学模频率与生成 S7.2 时差 %.1f%%（> %.0f%%）：位移是按旧本征矢做的、ħω 也是旧的，"
                 "retry S7.2_odp（重算单点）再重跑本步" % (100 * rel, 100 * FREQ_FAIL))
        if rel > FREQ_WARN:
            print("[WARN] S5 现在的光学模频率与生成 S7.2 时差 %.1f%%（> %.0f%%）：S5 重算过，建议 retry S7.2"
                  % (100 * rel, 100 * FREQ_WARN))
        return rel
    except Exception as e:                                      # noqa: BLE001
        print("[WARN] 核对 S5 频率失败（%s: %s），跳过" % (type(e).__name__, e))
        return None


def valley_kpoints(v, ref):
    """S8.2 谷表的一个星 -> 各成员谷在本步网格上的 k。优先 k_S7（S8.2 已吸附到 S7 网格），旧版没有就吸附 k_grid。"""
    import numpy as np
    ks = v.get("k_S7") or []
    if ks:
        return [np.asarray(k, float) for k in ks], 0.0
    recip = np.asarray(ref.structure.lattice.reciprocal_lattice.matrix)
    out, worst = [], 0.0
    cand = np.vstack([ref.kpts, -ref.kpts])
    for k in v.get("k_grid") or v.get("k_frac") or []:
        dk = cand - np.asarray(k, float)
        dk -= np.round(dk)
        dc = np.linalg.norm(dk @ recip, axis=1)
        i = int(np.argmin(dc))
        if dc[i] > SNAP_TOL_INVA:
            _die("谷 %s 的 k=%s 离本步网格最近点 %.3f Å⁻¹（> %.2f）" % (v.get("label"), k, dc[i], SNAP_TOL_INVA))
        out.append(np.asarray(k, float) + dk[i])
        worst = max(worst, float(dc[i]))
    return out, worst


def carrier_params(car, valleys, cfg, man, s7ref, warns):
    import numpy as np
    sel = {m["index"]: m for m in man["selected"]}
    groups = [g for g in man["groups"] if g.get("use")]
    name = {}
    for nm, c in man["configs"].items():
        if c["mode"] is not None:
            name[(c["mode"], c["sign"])] = nm
    ref = cfg[REF_NAME]
    out = []
    for v in valleys:
        spin = str(v.get("spin", "up")).lower()
        band = int(v["band_index"])
        ks, snap = valley_kpoints(v, ref)
        if snap:
            warns.append("%s %s：S8.2 没有 k_S7（旧版），k_grid 吸附到本步网格（最大偏 %.3f Å⁻¹）" % (car, v["label"], snap))
        # 参考构型 vs S7 undeformed：同一 k、同一带
        if s7ref is not None:
            for k in ks:
                e_a = float(ref.bands(k, spin)[band])
                e_b = float(s7ref.bands(k, spin)[band])
                if abs(e_a - e_b) > REF_TOL_EV:
                    _die("%s 谷 %s：odp-ref 的带 %d 能量 %.4f eV ≠ S7 undeformed 的 %.4f eV（差 %.1f meV > %.0f）——"
                         "两边输入/结构不一致，retry S7.2" % (car, v["label"], band, e_a, e_b, 1000 * abs(e_a - e_b),
                                                           1000 * REF_TOL_EV))
        grows = []
        for g in groups:
            d2, d2raw, info = {}, {}, []
            for i in g["modes"]:
                cp, cm = name[(i, 1)], name[(i, -1)]
                u = float(man["configs"][cp]["u_A"])
                dvp, dvm = cfg[cp].vac - ref.vac, cfg[cm].vac - ref.vac
                d2[i], d2raw[i] = [], []
                for k in ks:
                    e0, ep, em = ref.bands(k, spin), cfg[cp].bands(k, spin), cfg[cm].bands(k, spin)
                    where = "%s 谷 %s（%s 带 %d，k=%s）模 %d" % (car, v["label"], spin, band,
                                                            [round(float(x), 4) for x in k], i)
                    try:
                        r = of.edge_coupling(e0, ep, em, band, u, dvac_p=dvp, dvac_m=dvm)
                    except of.BandCrossing as e:
                        dm = float(man.get("disp_max_A") or 0.01)
                        _die("%s：%s。把 S7.2 的 ODP_DISP_A 调小（现在 %.3g Å，可试 %.3g）后 retry S7.2_odp"
                             % (where, e, dm, dm / 2))
                    except ValueError as e:
                        _die("%s：%s" % (where, e))
                    try:
                        rr = of.edge_coupling(e0, ep, em, band, u)
                    except ValueError:                 # 不对齐的值只作参考：解不出来就记 nan
                        rr = {"D2": float("nan")}
                    d2[i].append(r["D2"])
                    d2raw[i].append(rr["D2"])
                    info.append(r)
                    if r["shift_max_eV"] > SHIFT_WARN_EV:
                        warns.append("%s %s 模 %d：带边移动 %.3f eV（> %.2f），可能出线性区，可把 ODP_DISP_A 调小"
                                     % (car, v["label"], i, r["shift_max_eV"], SHIFT_WARN_EV))
                    if r["widened"]:
                        print("[..] %s：近邻带 %s 已核对，没有交叉" % (where, r["near"]))
            D, per, spread = of.combine_members(d2)
            Draw, _p, _s = of.combine_members(d2raw)
            if spread > of.SPREAD_WARN and D > of.DMIN_EV_A:
                warns.append("%s %s 组 %d（模 %s）：各成员谷 ΣD² 相差 %.0f%%（对称等价的谷应相同）"
                             % (car, v["label"], g["group"], g["modes"], 100 * spread))
            grows.append({"group": g["group"], "modes": g["modes"], "hw_meV": round(float(g["hw_meV"]), 4),
                          "freq_THz": round(float(g["freq_THz"]), 5), "polar": g.get("polar"),
                          "D_eV_A": round(D, 5), "D_raw_eV_A": round(Draw, 5) if Draw == Draw else None,
                          "D2_members": [round(x, 6) for x in per], "spread": round(spread, 4),
                          "method": sorted({r["method"] for r in info}),
                          "shift_max_eV": round(max(r["shift_max_eV"] for r in info), 5),
                          "u_A": [round(float(man["configs"][name[(i, 1)]]["u_A"]), 6) for i in g["modes"]]})
        out.append({"label": v["label"], "spin": spin, "band_index": band, "g": len(ks),
                    "dE_eV": v.get("dE_eV"), "k_used": [np.round(k, 6).tolist() for k in ks], "groups": grows})
    return out


def main():
    global of
    cwd = Path.cwd()
    kc.reexec_with_pymatgen(__file__)
    import odp_frozen as of                                    # noqa: F811
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    src = cwd / SRC_DIR
    try:
        man = json.loads((src / MANIFEST).read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        _die("读不了 %s/%s（S7.2 没生成？）：%s" % (SRC_DIR, MANIFEST, e))
    try:
        dpt = json.loads((cwd / DPT_RESULT).read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        _die("读不了 %s（S8.2 没跑完？）：%s" % (DPT_RESULT, e))
    check_inputs(cwd, man)
    s5_rel = check_s5(cwd, man)
    dim, vac_axis = kc.resolve_dim_for(src / "POSCAR", "auto")
    if dim != "2d":
        _die("本步只给 2D 用（%s 判为 %s）" % (SRC_DIR, dim))
    cfg = {}
    try:
        cache = json.loads((out / VAC_CACHE).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        cache = {}
    for nm in man["configs"]:
        d = src / nm
        if not d.is_dir():
            _die("缺 %s/%s（manifest 里有，目录没了）" % (SRC_DIR, nm))
        cfg[nm] = Config(d, vac_axis, cache)
        print("[..] %-12s 真空能级 %.5f eV（std %.1e）" % (nm, cfg[nm].vac, cfg[nm].vac_std))
    (out / VAC_CACHE).write_text(json.dumps(cache, indent=1), encoding="utf-8")
    s7ref = None
    und = cwd / DEFORM_DIR / "undeformed"
    if _vasprun(und) is not None:                       # 只要本征值（LOCPOT 在 S7.1 之后可能已清理）
        from pymatgen.io.vasp import Vasprun
        import numpy as np
        v = Vasprun(str(_vasprun(und)), parse_dos=False, parse_potcar_file=False)
        s7ref = object.__new__(Config)
        s7ref.name, s7ref.kpts = "%s/undeformed" % DEFORM_DIR, np.asarray(v.actual_kpoints, float)
        s7ref.eig = {str(getattr(sp, "name", sp)).lower(): np.asarray(a)[:, :, 0] for sp, a in v.eigenvalues.items()}
    else:
        print("[WARN] 没有 %s/undeformed/vasprun.xml：跳过 odp-ref 与 S7 的本征值核对" % DEFORM_DIR)

    warns, carriers = [], {}
    for r in dpt.get("results", []):
        car = r.get("carrier")
        valleys = ((r.get("multi_valley") or {}).get("valleys")) or []
        if not valleys:
            carriers[car] = {"status": "no_valleys",
                             "note": "S8.2 没有 %s 的多谷谷表（multi_valley.valleys）：本载流子不出冻声子参数" % car}
            continue
        grid = all(v.get("k_grid") for v in valleys)
        carriers[car] = {"status": "ok",
                         "valley_sig": {"use_grid": grid, "sig": json.loads(json.dumps(kc._valley_sig(valleys, grid)))},
                         "valleys": carrier_params(car, valleys, cfg, man, s7ref, warns)}
    params = {
        "schema": of.SCHEMA, "skill_rev": _SKILL_REV, "created": _time.strftime("%Y-%m-%d %H:%M:%S"),
        "valleys_source": DPT_RESULT, "manifest": "%s/%s" % (SRC_DIR, MANIFEST),
        "disp_max_A": man.get("disp_max_A"), "polar": man.get("polar"),
        "groups": man["groups"], "modes": man["modes"],
        "vacuum_levels_eV": {k: round(c.vac, 6) for k, c in cfg.items()},
        "s5_freq_rel_change": s5_rel, "warnings": warns, "carriers": carriers,
        "units": {"D_eV_A": "eV/Å（= 1e8 eV/cm），与 step.conf DP_OPTICAL_* 同单位", "hw_meV": "meV"},
        "definition": ("D = d(E−V_vac)/du，位移 Δ_κ = u·√(M/m_κ)·e_κ（M 原胞质量、e 质量加权归一本征矢）；"
                       "简并组各伙伴 ΣD²、星内各成员谷平均；简并带 D² = Σλ²/d²（见 odp_frozen.py）"),
    }
    (out / "odp_params.json").write_text(json.dumps(params, ensure_ascii=False, indent=1), encoding="utf-8")
    lines = ["# 冻声子 ODP 参数（%s）" % params["created"],
             "# D 单位 eV/Å（= 1e8 eV/cm），ħω 单位 meV；真空对齐；位移 %.3f Å；极性模 %s" % (
                 float(man.get("disp_max_A") or 0), man.get("polar"))]
    for car, c in carriers.items():
        if c["status"] != "ok":
            lines.append("%s: %s" % (car, c.get("note", c["status"])))
            continue
        for v in c["valleys"]:
            for g in v["groups"]:
                raw = "—" if g["D_raw_eV_A"] is None else "%6.3f" % g["D_raw_eV_A"]
                lines.append("%-8s %-10s ×%d  组%-2d 模 %-7s ħω=%6.2f meV  D=%6.3f eV/Å（不对齐 %s）  成员离散 %.1f%%  %s"
                             % (car, v["label"], v["g"], g["group"], ",".join(map(str, g["modes"])), g["hw_meV"],
                                g["D_eV_A"], raw, 100 * g["spread"], "/".join(g["method"])))
    for w in warns:
        lines.append("[WARN] " + w)
        print("[WARN] " + w)
    (out / "odp_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[2:]))
    kc.invalidate_odp_consumer(cwd, params)                 # S8.4 用过的参数变了 / 刚可用 -> 让它重排
    print("[DONE] %s/odp_params.json（S8.4 的 ODP_FROZEN=auto 会读它）" % OUTDIR_NAME)


if __name__ == "__main__":
    main()
