#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step12c_odp_read.py —— [0035] 读冻声子单点，算谷内光学形变势（step8.21_odp_read，S8.21_odp）。

run: gen 步骤（登录节点直接跑，秒级到分钟级）。输入：
  step7c_odp_frozen/{odp-ref, odp-mNNp, odp-mNNm}/vasprun.xml + LOCPOT + odp_manifest.json（S7.2）
  step8.2_dpt/dpt_result.json 的多谷谷表（每个谷的自旋、带号、各成员谷在 S7 网格上的 k_S7）
产出：step8.21_odp_read/odp_params.json（done_marker）+ odp_summary.txt。S8.4 的 ODP_FROZEN=auto 读它。

[0036] 三维（S8 的 ODP_FROZEN=auto 读它）：
  · 不对齐（没有真空，D = dE/du 直接用本征值差分；理由见 odp_frozen.py 文件头），不要 LOCPOT；
  · 谷表：材料 step.conf 手填了 DP_OPTICAL_VALLEYS_E/H 就用它（中心吸附到 S7 网格），否则在 odp-ref 的网格
    （= S7 undeformed，ISYM=0）上自动找带边极值（odp_frozen.find_valleys_3d，窗口 ODP_VALLEY_WINDOW_EV）。
    S8 直接用这里记下的谷中心（k_frac）作 ODP 的谷表，并按 valley_sig 核对手填谷表没改过。不需要 S8.2。

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
_SKILL_REV = "2026-10-12-0036"   # [0036] 三维（不对齐、自动/手填谷表）；[0035-r3]
VALLEY_SPEC = {"DP_OPTICAL_VALLEYS_E": ("", "str"), "DP_OPTICAL_VALLEYS_H": ("", "str"),
               "ODP_VALLEY_WINDOW_EV": (None, "float")}     # [0036] 3D 谷表（材料 step.conf，与 S8 同名同义）
REF_TOL_EV = 0.005            # odp-ref 与 S7 undeformed 同一 k、同一带的本征值差上限
REF_TOL_WARN_EV = 1e-4        # 更严告警（0.1 meV）：超过就提示。多半是两台机器的 VASP 编译版本不同
                              #   （如 AOCC/AOCL vs Intel/MKL）；D=(E₊−E₋)/2u、u≤0.01 Å，会放大 ~50–100×
SNAP_TOL_INVA = 0.02          # 旧 S8.2（没有 k_S7）时把 k_grid 吸附到本步网格的容差（Å⁻¹）
SHIFT_WARN_EV = 0.1           # 带边移动超过它：出线性区，建议减小 ODP_DISP_A
FREQ_WARN = 0.02              # S5 现在的频率与生成时相差超过它 -> 告警
FREQ_FAIL = 0.05              # ……超过它 -> 报错（位移用的是旧本征矢，要 retry S7.2）
DIPOLE_NOTE_EV = 0.005        # [0035-r3] 真空区斜坡（每胞）超过它就提示该构型有面外偶极
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

    def __init__(self, d, vac_axis=2, cache=None, use_vac=True):
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
        self.eig, self.occ = {}, {}
        for sp, arr in v.eigenvalues.items():
            self.eig[str(getattr(sp, "name", sp)).lower()] = np.asarray(arr)[:, :, 0]
            self.occ[str(getattr(sp, "name", sp)).lower()] = np.asarray(arr)[:, :, 1]
        self.structure = v.final_structure
        if not use_vac:                     # [0036] 3D：没有真空，不对齐
            self.vac, self.vac_std, self.vac_slope = 0.0, 0.0, 0.0
            return
        lp = d / "LOCPOT"
        key = _sha256(vr)
        ent = (cache or {}).get(d.name) or {}
        if lp.is_file():
            planar = np.asarray(Locpot.from_file(str(lp)).get_average_along_axis(vac_axis))
            zf = np.asarray(self.structure.frac_coords)[:, vac_axis]
            try:
                self.vac, self.vac_std, self.vac_slope = of.vacuum_level(planar, zf, return_slope=True)
            except ValueError as e:
                _die("%s/%s：%s" % (SRC_DIR, d.name, e))
            if cache is not None:
                cache[d.name] = {"vasprun_sha256": key, "vac_eV": self.vac, "std_eV": self.vac_std,
                                 "slope_eV": self.vac_slope}
        elif ent.get("vasprun_sha256") == key:
            self.vac, self.vac_std = float(ent["vac_eV"]), float(ent["std_eV"])
            self.vac_slope = float(ent.get("slope_eV", 0.0))
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
            continue
        if not b.is_file():
            bad.append(n)
            continue
        _ta, _tb = a.read_text(errors="ignore"), b.read_text(errors="ignore")
        if kc._input_key(n, _ta) == kc._input_key(n, _tb):
            continue
        # [kz=1] 2D：放行「面内网格相同、仅真空方向 kz 被降到 1」这一种差异（见 gen_step9c_odp.py）。
        #   真空方向能带平坦已由 LS 实测确认（带边 ±1 eV 内 |dE| ~0.07 meV）。
        # [kz=1-r2] 只在 S7.2 的 manifest 记了「真的降过」时放行，且按 manifest 记的真空轴比；
        #   旧 manifest（没有这个标记）KPOINTS 不同照旧报错。
        if n == "KPOINTS" and str((man or {}).get("dim")) == "2d" and (man or {}).get("kz_vacuum_reduced"):
            _ax = ((man.get("kz_vacuum") or {}).get("axis"))
            _ax = 2 if _ax is None else int(_ax)
            _ga, _gb = of.kpoints_mesh(_ta), of.kpoints_mesh(_tb)
            if (_ga and _gb and _gb[_ax] == 1
                    and all(_ga[i] == _gb[i] for i in range(3) if i != _ax)):
                print("[..] KPOINTS：面内一致（%s），真空方向 %d -> 1（kz=1 优化）放行"
                      % ([_gb[i] for i in range(3) if i != _ax], _ga[_ax]))
                continue
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


def _read_magnetism(d):
    """[kz=1-r2] 一个构型的磁矩：OSZICAR 最后一个 mag=（总磁矩，非共线是 3 个分量）与 OUTCAR 最后一个
    "magnetization (x)" 表的 tot 列（各原子）。非自旋极化（没有 mag=）返回 None。"""
    tot = None
    osz = d / "OSZICAR"
    if osz.is_file():
        for ln in osz.read_text(errors="ignore").splitlines():
            if "mag=" in ln:
                try:
                    tot = [float(x) for x in ln.split("mag=", 1)[1].split()]
                except ValueError:
                    pass
    if tot is None:
        return None
    sites = None
    oc = d / "OUTCAR"
    if oc.is_file():
        with open(oc, "rb") as fh:                             # 只读文件尾：磁矩表在单点 OUTCAR 的末尾
            fh.seek(0, 2)
            fh.seek(max(0, fh.tell() - (8 << 20)))
            tail = fh.read().decode(errors="ignore").splitlines()
        start = max((i for i, ln in enumerate(tail) if ln.strip().startswith("magnetization (x)")), default=None)
        if start is not None:
            vals = []
            for ln in tail[start + 1:]:
                f = ln.split()
                if not f or f[0].startswith(("#", "-")):
                    if vals and f and f[0].startswith("-"):
                        break
                    continue
                if f[0] == "tot":
                    break
                if f[0].isdigit():
                    try:
                        vals.append(float(f[-1]))
                    except ValueError:
                        break
            sites = vals or None
    return {"total": tot, "sites": sites}


def check_magnetism(src, names):
    """[kz=1-r2] 各构型与 odp-ref 是否停在同一个磁性解上（odp_frozen.magnetic_mismatch）。不一致就报错：
    换了磁性解，带边能量整体改变，D = ΔE/2u 会变成一个假的大耦合。非自旋极化直接跳过。返回核对记录。"""
    ref = _read_magnetism(src / REF_NAME)
    if ref is None:
        print("[..] 非自旋极化（odp-ref 的 OSZICAR 没有 mag=）：跳过磁性解一致性核对")
        return {"checked": False}
    bad = {}
    for nm in names:
        if nm == REF_NAME:
            continue
        m = _read_magnetism(src / nm)
        if m is None:
            bad[nm] = ["读不到磁矩（OSZICAR 没有 mag=），而 odp-ref 是自旋极化的"]
            continue
        probs = of.magnetic_mismatch(ref, m)
        if probs:
            bad[nm] = probs
    if bad:
        _die("下列冻声子构型与 odp-ref 不在同一个磁性解上（位移只有 ~0.01 Å，同一磁性解下磁矩只变 ~1e-3 μB）：\n  "
             + "\n  ".join("%s：%s" % (k, "；".join(v)) for k, v in sorted(bad.items()))
             + "\n这些点的 D 会是假的大耦合。处理：查它们的 OSZICAR/OUTCAR 磁矩；用 odp-ref 收敛后的各原子磁矩作 "
               "MAGMOM（写进 S7.2 的 step.conf [incar]）重算这些点，再重跑本步")
    print("[OK] 磁性解一致：%d 个构型与 odp-ref 的总磁矩/各原子磁矩差都在 %.2f/%.2f μB 内（odp-ref 总磁矩 %s）"
          % (len(names) - 1, of.MAG_TOTAL_TOL, of.MAG_SITE_TOL, [round(x, 3) for x in ref["total"]]))
    return {"checked": True, "ref_total": ref["total"], "tol_total": of.MAG_TOTAL_TOL, "tol_site": of.MAG_SITE_TOL}


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
                # [kz=1] kz 被降过时（odp-ref 用 N1 N2 1、S7 undeformed 用 N1 N2 3）本身就有约
                #   0.07 meV 的网格差（LS 实测），阈值放宽一倍，免得误报。
                _wr = REF_TOL_WARN_EV * (2.0 if (man or {}).get("kz_vacuum_reduced") else 1.0)
                if _wr < abs(e_a - e_b) <= REF_TOL_EV:
                    warns.append("%s 谷 %s：odp-ref vs S7 undeformed 的带 %d 差 %.3f meV（> %.1f meV 严告警，"
                                 "仍 < %.0f meV 判据）%s—— 若两者在不同机器/不同 VASP 编译版本上跑，"
                                 "D=(E₊−E₋)/2u 会放大该偏差几十倍，建议同机重算"
                                 % (car, v["label"], band, 1000 * abs(e_a - e_b),
                                    1000 * _wr, 1000 * REF_TOL_EV,
                                    "（kz=1 网格本身贡献约 0.07 meV，阈值已放宽）"
                                    if (man or {}).get("kz_vacuum_reduced") else ""))
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


def carriers_3d(cwd, cfg, man, s7ref, warns):
    """[0036] 三维：谷表（手填优先，否则自动）-> 每个载流子的冻声子参数。"""
    import numpy as np
    import stepconf
    p = {}
    if (cwd / stepconf.CONF_NAME).is_file():
        p = stepconf.load(VALLEY_SPEC, OUTDIR_NAME, str(cwd), strict=False)
    win = float(p["ODP_VALLEY_WINDOW_EV"]) if p and p["ODP_VALLEY_WINDOW_EV"] is not None else of.VALLEY_WINDOW_EV
    ref = cfg[REF_NAME]
    spins = list(ref.eig)
    nocc = {}
    for sp in spins:
        n_k = np.count_nonzero(ref.occ[sp] > 0.5, axis=1)
        if int(n_k.min()) != int(n_k.max()):
            _die("odp-ref 自旋 %s 各 k 点占据带数不同（%d..%d）：金属/半金属，ODP 常数形变势模型不适用"
                 % (sp, int(n_k.min()), int(n_k.max())))
        nocc[sp] = int(n_k[0])
    recip = np.asarray(ref.structure.lattice.reciprocal_lattice.matrix)
    out = {}
    for car, key, sign in (("electron", "DP_OPTICAL_VALLEYS_E", 1), ("hole", "DP_OPTICAL_VALLEYS_H", -1)):
        spec = str((p[key] if p else "") or "").strip()
        if spec:
            import dp_optical as dpo
            try:
                vs = dpo.parse_valleys(spec, full3d=True)
            except ValueError as e:
                _die("%s：%s" % (key, e))
            src = "step.conf %s" % key
            sp0 = spins[0]
            band = nocc[sp0] if sign > 0 else nocc[sp0] - 1
            for v in vs:
                v.update(spin=sp0, band_index=band, dE_eV=None)
        else:
            try:
                vs = of.find_valleys_3d(ref.kpts, ref.eig, nocc, sign, recip, window=win)
            except ValueError as e:
                _die("%s 自动找谷失败：%s。处理：在 step.conf 手填 %s" % (car, e, key))
            src = "auto（odp-ref 网格，窗口 %.2f eV）" % win
            print("[..] %s 谷表（自动）：%s" % (car, "; ".join(
                "%s×%d dE=%.3f eV %s 带 %d" % (v["label"], v["g"], v["dE_eV"], v["spin"], v["band_index"]) for v in vs)))
        if len(vs) > 1 or any(int(v.get("g", 1)) > 1 for v in vs):     # [0042] 单星多成员（Si X×6）也是多谷
            warns.append("%s：%d 组谷（%s）—— 只有谷内 ODP；谷间 IVS 不自动算，要的话按文献在 step.conf 写 DP_OPTICAL_%s 的 "
                         "X>X'（同星其它谷）或 X>Y 通道（谷标签用本文件的）" % (car, len(vs), ",".join("%s×%d" % (v["label"], int(v.get("g", 1))) for v in vs),
                                                         "E" if sign > 0 else "H"))
        vals = carrier_params(car, vs, cfg, man, s7ref, warns)
        for row, v in zip(vals, vs):
            row["k_frac"] = [np.round(np.asarray(k, float), 6).tolist() for k in v["k_frac"]]
        out[car] = {"status": "ok", "valleys_source": src,
                    "valley_sig": {"dim": "3d", "sig": of.valley_sig_3d(vs)}, "valleys": vals}
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
    dim, vac_axis = kc.resolve_dim_for(src / "POSCAR", "auto")
    if dim not in ("2d", "3d"):
        _die("本步只给 2D/3D 用（%s 判为 %s）" % (SRC_DIR, dim))
    is3d = dim == "3d"
    dpt = None
    if not is3d:
        try:
            dpt = json.loads((cwd / DPT_RESULT).read_text(encoding="utf-8"))
        except (OSError, ValueError) as e:
            _die("读不了 %s（S8.2 没跑完？）：%s" % (DPT_RESULT, e))
    check_inputs(cwd, man)
    s5_rel = check_s5(cwd, man)
    cfg = {}
    try:
        cache = json.loads((out / VAC_CACHE).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        cache = {}
    for nm in man["configs"]:
        d = src / nm
        if not d.is_dir():
            _die("缺 %s/%s（manifest 里有，目录没了）" % (SRC_DIR, nm))
        cfg[nm] = Config(d, vac_axis, cache, use_vac=not is3d)
        if is3d:
            continue
        slope = float(getattr(cfg[nm], "vac_slope", 0.0))
        print("[..] %-12s 真空能级 %.5f eV（去斜坡残差 %.1e，真空区斜坡 %+.1f meV/胞）"
              % (nm, cfg[nm].vac, cfg[nm].vac_std, 1000 * slope))
        if abs(slope) > DIPOLE_NOTE_EV:
            print("[..] %s 有面外偶极（真空区斜坡 %+.1f meV/胞）：取真空区中点为参考（两侧真空能级的平均）"
                  % (nm, 1000 * slope))
    (out / VAC_CACHE).write_text(json.dumps(cache, indent=1), encoding="utf-8")
    mag_rec = check_magnetism(src, list(man["configs"]))        # [kz=1-r2] 各构型同一个磁性解
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
    if is3d:
        carriers = carriers_3d(cwd, cfg, man, s7ref, warns)
    for r in ([] if is3d else dpt.get("results", [])):
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
        "dim": dim, "aligned": "none (3D)" if is3d else "vacuum",
        "valleys_source": "3D：见 carriers.<载流子>.valleys_source" if is3d else DPT_RESULT,
        "manifest": "%s/%s" % (SRC_DIR, MANIFEST),
        "disp_max_A": man.get("disp_max_A"), "polar": man.get("polar"),
        "groups": man["groups"], "modes": man["modes"],
        "vacuum_levels_eV": None if is3d else {k: round(c.vac, 6) for k, c in cfg.items()},
        "s5_freq_rel_change": s5_rel, "magnetism_check": mag_rec, "warnings": warns, "carriers": carriers,
        "units": {"D_eV_A": "eV/Å（= 1e8 eV/cm），与 step.conf DP_OPTICAL_* 同单位", "hw_meV": "meV"},
        "definition": (("D = dE/du（3D 不对齐），" if is3d else "D = d(E−V_vac)/du，") + "位移 Δ_κ = u·√(M/m_κ)·e_κ（M 原胞质量、e 质量加权归一本征矢）；"
                       "简并组各伙伴 ΣD²、星内各成员谷平均；简并带 D² = Σλ²/d²（见 odp_frozen.py）"),
    }
    (out / "odp_params.json").write_text(json.dumps(params, ensure_ascii=False, indent=1), encoding="utf-8")
    lines = ["# 冻声子 ODP 参数（%s）" % params["created"],
             "# D 单位 eV/Å（= 1e8 eV/cm），ħω 单位 meV；%s；位移 %.3f Å；极性模 %s" % (
                 "3D 不对齐" if is3d else "真空对齐", float(man.get("disp_max_A") or 0), man.get("polar"))]
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
    if is3d:                                                # [0036] 3D 的消费者是 S8（记录在 dp_optical_3d.json）
        kc.invalidate_odp_consumer(cwd, params, consumer="step8_amset", record="dp_optical_3d.json")
    else:
        kc.invalidate_odp_consumer(cwd, params)             # S8.4 用过的参数变了 / 刚可用 -> 让它重排
    print("[DONE] %s/odp_params.json（%s 的 ODP_FROZEN=auto 会读它）" % (OUTDIR_NAME, "S8" if is3d else "S8.4"))


if __name__ == "__main__":
    main()
