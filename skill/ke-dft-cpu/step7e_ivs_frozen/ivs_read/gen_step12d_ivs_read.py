#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step12d_ivs_read.py —— [0050] 读谷间冻声子单点，算 K↔K′ 零阶谷间形变势（step8.22_ivs_read，S8.22_ivs）。

run: gen 步骤（登录节点直接跑，秒级）。输入：
  step7e_ivs_frozen/{ivs-ref, ivs-prim, ivs-mNNp, ivs-mNNm}/vasprun.xml + ivs_manifest.json（S7.4）
  step7d_ivs_phonon/vasprun.xml（S7.3，只用来核对频率没变）
  step8.2_dpt/dpt_result.json 的多谷谷表（找 K 星：自旋、带号）
产出：step8.22_ivs_read/ivs_params.json（done_marker）+ ivs_summary.txt。S8.4 的 IVS_FROZEN=auto 读它。

做法（推导与归一化见 ivs_frozen.py）：
  1. ivs-prim（原胞，网格含全部折叠 k）给出 K、K′、Γ 等折叠 k 的能级；它们并起来排序 = 超胞 Γ̃ 的能级（差一个常数），
     据此给 ivs-ref 的 Γ̃ 能级标来源，找到 K/K′ 的那条带在超胞里的二重态（ivs_frozen.label_folded / edge_group）；
  2. 每个区边界模 ±u：二重态劈裂 ±|V|u，Σλ² = 组内去掉共同平移后的迹（与参考能级无关，不要真空对齐）；
  3. 简并伙伴模的 Σλ² 相加：D² = ΣΣλ²/(g(g−1))（g = 2）。
核对（不过就报错，绝不静默出数）：
  · S7.4 记下的 S7 undeformed 输入（INCAR/POTCAR/POSCAR 的 sha256）与现在一致；S7.3 的频率与生成位移时一致（差 >5% 报错）；
  · ivs-prim 在 K 点的带边与 S7 undeformed 一致（< 30 meV：带号/自旋/结构对上了）；
  · 超胞 Γ̃ 能级与原胞折叠能级逐条对上（< 50 meV）、K/K′ 二重态简并且附近没有别的来源的能级；
  · 非 SOC（[0050] A 方案）。
只算 K↔K′ 零阶：K↔Q、K↔Γ 与一阶谷间耦合不在内。离 K 星 IVS_NEAR_STAR_EV 以内还有别的星（如 MoS₂ 导带的 Q）-> 告警。
"""
import hashlib
import json
import os
import sys
import time as _time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, os.getcwd())
import ke_common as kc  # noqa: E402

iv = None                     # ivs_frozen：main() 里切好 python 环境后再 import（它要 numpy）

OUTDIR_NAME = "step8.22_ivs_read"
SRC_DIR = "step7e_ivs_frozen"
PH_DIR = "step7d_ivs_phonon"
DEFORM_DIR = "step7_deform"
DPT_RESULT = "step8.2_dpt/dpt_result.json"
MANIFEST = "ivs_manifest.json"
REF_NAME = "ivs-ref"
PRIM_NAME = "ivs-prim"
_SKILL_REV = "2026-10-14-0050"   # [0050] 谷间 IVS（K↔K′ 零阶）冻声子：读取
PRIM_TOL_EV = 0.03            # ivs-prim 与 S7 undeformed 在 K 点同一带的能量差上限（两边 k 网格不同，SCF 密度略不同）
SHIFT_WARN_EV = 0.1           # 二重态劈裂超过它：出线性区，建议减小 IVS_DISP_A
FREQ_WARN = 0.02              # S7.3 现在的频率与生成位移时相差超过它 -> 告警
FREQ_FAIL = 0.05              # ……超过它 -> 报错（位移用的是旧本征矢）
NEAR_STAR_EV = 0.2            # 离 K 星这么近的别的星（K↔Q、K↔Γ 不在本方法内）-> 告警
STAR_TOL = 1e-3               # S8.2 的 k_grid / k_S7 与 K/K′ 的分数坐标容差
STAR_TOL_FIT = 0.02           # 只有拟合中心 k_frac 时的容差
OUTCAR_DONE = "General timing and accounting informations"


def _die(msg):
    sys.exit("[ERROR] " + msg)


def _vasprun(d):
    for n in ("vasprun.xml", "vasprun.xml.gz"):
        if (d / n).is_file():
            return d / n
    return None


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Config(object):
    """一个单点构型：本征值（按自旋）与 k 列表。"""

    def __init__(self, d):
        import numpy as np
        from pymatgen.io.vasp import Vasprun
        self.name = d.name
        oc = d / "OUTCAR"
        if not oc.is_file() or OUTCAR_DONE not in oc.read_text(errors="ignore")[-400000:]:
            _die("%s/%s 的 OUTCAR 没跑完（S7.4 还有单点没算完？）" % (SRC_DIR, d.name))
        vr = _vasprun(d)
        if vr is None:
            _die("缺 %s/%s/vasprun.xml" % (SRC_DIR, d.name))
        v = Vasprun(str(vr), parse_dos=False, parse_potcar_file=False)
        if not v.converged_electronic:
            _die("%s/%s 的电子步没收敛（跑到 NELM 了）：本征值不可信。查 OSZICAR，调 NELM/ALGO 后 retry S7.4_ivs"
                 % (SRC_DIR, d.name))
        self.kpts = np.asarray(v.actual_kpoints, float)
        self.eig = {str(getattr(sp, "name", sp)).lower(): np.asarray(arr)[:, :, 0] for sp, arr in v.eigenvalues.items()}
        self.structure = v.final_structure

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
            _die("k=%s 不在 %s 的 k 网格（含 −k）上 —— S7.4 生成后网格被改过？retry S7.4_ivs"
                 % ([round(float(x), 5) for x in kf], self.name))
        if spin not in self.eig:
            _die("%s 没有自旋道 %s（ISPIN 与 S3/S7 不同？）" % (self.name, spin))
        return self.eig[spin][i]


def check_inputs(cwd, man):
    """S7.4 生成时记下的 S7 undeformed 输入（sha256）与现在比。"""
    und = cwd / DEFORM_DIR / "undeformed"
    old = (man.get("inputs_from") or {}).get("sha256") or {}
    bad = []
    for n, h in sorted(old.items()):
        f = und / n
        if not f.is_file():
            bad.append("%s 不在了" % n)
        elif _sha(f) != h:
            bad.append(n)
    if bad:
        _die("S7 undeformed 的 %s 与生成 S7.4 时不同：S7 重新生成过 -> 先 retry S7.3_ivsph（超胞 DFPT），"
             "再 retry S7.4_ivs（单点），然后重跑本步" % "、".join(bad))


def check_s73(cwd, man):
    """S7.3 现在的频率 vs 生成位移时记下的：vasprun.xml 没变直接过；差 > FREQ_WARN 告警，> FREQ_FAIL 或模数不同报错。
    S7.3 读不了（被清理）只告警。"""
    import numpy as np
    rec = man.get("s73") or {}
    vr = _vasprun(cwd / PH_DIR)
    if vr is None:
        print("[WARN] 没有 %s/vasprun.xml：跳过 S7.3 频率核对" % PH_DIR)
        return None
    if _sha(vr) == rec.get("vasprun_sha256"):
        return 0.0
    try:
        from pymatgen.io.vasp import Vasprun
        v = Vasprun(str(vr), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
        ev = -np.asarray(v.normalmode_eigenvals)[::-1]
        f = np.sqrt(np.abs(ev)) * np.sign(ev)
        if int(str(v.vasp_version).split(".")[0]) < 6:
            f = f * 15.633302
    except Exception as e:                                      # noqa: BLE001
        print("[WARN] 核对 S7.3 频率失败（%s: %s），跳过" % (type(e).__name__, e))
        return None
    old = np.sort(np.asarray(rec.get("freq_THz") or [], float))
    f = np.sort(f)
    if len(old) != len(f):
        _die("S7.3 的模数与生成 S7.4 时不同（%d vs %d）：结构变了，retry S7.4_ivs" % (len(f), len(old)))
    opt = np.abs(old) > 1.0
    rel = float(np.max(np.abs(f[opt] - old[opt]) / np.abs(old[opt]))) if opt.any() else 0.0
    if rel > FREQ_FAIL:
        _die("S7.3 现在的频率与生成 S7.4 时差 %.1f%%（> %.0f%%）：位移是按旧本征矢做的、ħω 也是旧的，"
             "retry S7.4_ivs（重算单点）再重跑本步" % (100 * rel, 100 * FREQ_FAIL))
    if rel > FREQ_WARN:
        print("[WARN] S7.3 现在的频率与生成 S7.4 时差 %.1f%%（> %.0f%%）：S7.3 重算过，建议 retry S7.4_ivs"
              % (100 * rel, 100 * FREQ_WARN))
    return rel


def _same_mod1(a, b, tol):
    import numpy as np
    d = np.asarray(a, float)[:2] - np.asarray(b, float)[:2]
    d -= np.round(d)
    return bool(np.all(np.abs(d) < tol))


def match_star(v, members):
    """S8.2 谷表的一个星是不是 K 星：各成员（模倒格矢）正好是 K、K′。优先 k_grid / k_S7（网格上的精确点），否则拟合中心。"""
    for key, tol in (("k_grid", STAR_TOL), ("k_S7", STAR_TOL), ("k_frac", STAR_TOL_FIT)):
        ks = v.get(key) or []
        if not ks:
            continue
        if len(ks) != len(members):
            return False
        return all(any(_same_mod1(k, m, tol) for m in members) for k in ks) and \
            all(any(_same_mod1(k, m, tol) for k in ks) for m in members)
    return False


def carrier_params(car, valleys, cfg, man, s7ref, warns):
    """一个载流子：在 S8.2 谷表里找 K 星（可能按自旋/带分成几个同中心的星），逐个算各组区边界模的 D。"""
    import numpy as np
    members = man["members_frac"]
    star = list(man["star"])
    folded = man["folded"]
    groups = [g for g in man["groups"] if g.get("use")]
    name = {}
    for nm, c in man["configs"].items():
        if c["mode"] is not None:
            name[(c["mode"], c["sign"])] = nm
    ref, prim_cfg = cfg[REF_NAME], cfg[PRIM_NAME]
    hits = [v for v in valleys if match_star(v, members)]
    if not hits:
        return None, "S8.2 的 %s 谷表里没有 K 星（%s）：本载流子不出 K↔K′ 谷间参数" % (
            car, ", ".join("%s×%d" % (v.get("label"), int(v.get("g", 1))) for v in valleys) or "空")
    for v in valleys:                         # [0050] 只算 K↔K′：离 K 星近的别的星（K↔Q、K↔Γ）提醒一句
        if v in hits or v.get("dE_eV") is None:
            continue
        dk = min(abs(float(v["dE_eV"]) - float(h.get("dE_eV") or 0.0)) for h in hits)
        if dk < NEAR_STAR_EV:
            warns.append("%s：%s 星离 K 星只有 %.0f meV（< %.0f）—— K↔%s 谷间散射不在本方法内（只算 K↔K′ 零阶），"
                         "该机制缺失、迁移率偏高；要的话按文献在 step.conf 写 DP_OPTICAL_%s 的 K>%s、%s>K 通道"
                         % (car, v.get("label"), 1000 * dk, 1000 * NEAR_STAR_EV, v.get("label"),
                            "E" if car == "electron" else "H", v.get("label"), v.get("label")))
    out = []
    for v in hits:
        spin = str(v.get("spin", "up")).lower()
        band = int(v["band_index"])
        prim = {}
        for kn in sorted(folded):
            prim[kn] = np.asarray(prim_cfg.bands(folded[kn], spin), float)
        if band >= len(prim[star[0]]) - 1:
            _die("%s 谷 %s 的带号 %d 超出 ivs-prim 的 NBANDS（%d）：加大 IVS_NBANDS_EXTRA 后 retry S7.4_ivs"
                 % (car, v.get("label"), band, len(prim[star[0]])))
        edge = float(prim[star[0]][band])
        if s7ref is not None:                 # ivs-prim vs S7 undeformed：同一 K、同一带（带号/自旋/结构对上了）
            i7 = s7ref.kidx(folded[star[0]])
            if i7 is not None and spin in s7ref.eig:
                e7 = float(s7ref.eig[spin][i7][band])
                if abs(edge - e7) > PRIM_TOL_EV:
                    _die("%s 谷 %s：ivs-prim 在 K 点带 %d 的能量 %.4f eV ≠ S7 undeformed 的 %.4f eV（差 %.0f meV > %.0f）"
                         "—— 两边输入/结构不一致，retry S7.4_ivs" % (car, v.get("label"), band, edge, e7,
                                                                1000 * abs(edge - e7), 1000 * PRIM_TOL_EV))
        e0 = np.asarray(ref.bands([0.0, 0.0, 0.0], spin), float)
        try:
            labels, shift = iv.label_folded(e0, prim, edge)
            grp = iv.edge_group(labels, e0, star, band)
        except ValueError as e:
            hint = "加大 step.conf 的 IVS_NBANDS_EXTRA 后 retry S7.4_ivs" if "NBANDS" in str(e) else \
                "ivs-ref 与 ivs-prim 的输入/结构是否一致（retry S7.4_ivs）"
            _die("%s 谷 %s（%s 带 %d）给超胞 Γ̃ 能级标来源失败：%s。处理：%s" % (car, v.get("label"), spin, band, e, hint))
        g = len(grp)
        grows = []
        for gr in groups:
            s2, info = {}, []
            for i in gr["modes"]:
                cp, cm = name[(i, 1)], name[(i, -1)]
                u = float(man["configs"][cp]["u_A"])
                ep = np.asarray(cfg[cp].bands([0.0, 0.0, 0.0], spin), float)
                em = np.asarray(cfg[cm].bands([0.0, 0.0, 0.0], spin), float)
                where = "%s 谷 %s（%s 带 %d）模 %d" % (car, v.get("label"), spin, band, i)
                try:
                    r = iv.split_sq(e0, ep, em, grp, u)
                except iv.BandCrossing as e:
                    dm = float(man.get("disp_max_A") or 0.01)
                    _die("%s：%s。把 S7.4 的 IVS_DISP_A 调小（现在 %.3g Å，可试 %.3g）后 retry S7.4_ivs" % (where, e, dm, dm / 2))
                except ValueError as e:
                    _die("%s：%s" % (where, e))
                s2[i] = r["S2"]
                info.append(r)
                if r["shift_max_eV"] > SHIFT_WARN_EV:
                    warns.append("%s %s 模 %d：二重态移动 %.3f eV（> %.2f），可能出线性区，可把 IVS_DISP_A 调小"
                                 % (car, v.get("label"), i, r["shift_max_eV"], SHIFT_WARN_EV))
                if r["widened"]:
                    print("[..] %s：近邻能级 %s 已核对，没有交叉" % (where, r["near"]))
            D = iv.group_D(s2, g)
            grows.append({"group": gr["group"], "modes": gr["modes"], "hw_meV": round(float(gr["hw_meV"]), 4),
                          "freq_THz": round(float(gr["freq_THz"]), 5), "D_eV_A": round(D, 5),
                          "S2_modes": {str(i): round(float(x), 6) for i, x in s2.items()},
                          "shift_max_eV": round(max(r["shift_max_eV"] for r in info), 5),
                          "near_checked": any(r["widened"] for r in info),
                          "u_A": [round(float(man["configs"][name[(i, 1)]]["u_A"]), 6) for i in gr["modes"]]})
        out.append({"label": v.get("label"), "spin": spin, "band_index": band, "g": g, "dE_eV": v.get("dE_eV"),
                    "star": star, "edge_prim_eV": round(edge, 5), "label_shift_eV": round(float(shift), 5),
                    "sc_levels": [int(i) for i in grp], "groups": grows})
    return out, None


def main():
    global iv
    cwd = Path.cwd()
    kc.reexec_with_pymatgen(__file__)
    import ivs_frozen as iv                                    # noqa: F811
    import numpy as np
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    src = cwd / SRC_DIR
    try:
        man = json.loads((src / MANIFEST).read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        _die("读不了 %s/%s（S7.4 没生成？）：%s" % (SRC_DIR, MANIFEST, e))
    dim = kc.read_method_dim(cwd / "step1_opt" / kc.METHOD_FILE) or \
        kc.resolve_dim_for(src / PRIM_NAME / "POSCAR", "auto")[0]
    if dim != "2d":
        _die("本步只给 2D（S8.4）用（判为 %s）" % dim)
    if kc.incar_is_ncl(src / REF_NAME / "INCAR") or kc.detect_ncl(cwd):
        _die("本材料是 SOC 能带（LSORBIT/LNONCOLLINEAR）。[0050] 谷间 IVS 冻声子只做 A 方案（非 SOC）：SOC 下 K/K′ 二重态"
             "被自旋劈裂、谷间散射只能落到同自旋子带（Xiao 等 PRL 108, 196802），劈裂法与 dp_optical 的星平均态密度都不适用。"
             "本材料 tf_*.yaml 写 ivs_frozen: false 关掉本组")
    try:
        dpt = json.loads((cwd / DPT_RESULT).read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        _die("读不了 %s（S8.2 没跑完？）：%s" % (DPT_RESULT, e))
    check_inputs(cwd, man)
    s73_rel = check_s73(cwd, man)
    cfg = {}
    for nm in man["configs"]:
        d = src / nm
        if not d.is_dir():
            _die("缺 %s/%s（manifest 里有，目录没了）" % (SRC_DIR, nm))
        cfg[nm] = Config(d)
    if [str(x.specie) for x in cfg[REF_NAME].structure] != list(man["species"]):
        _die("%s 的原子顺序与 %s 不同：retry S7.4_ivs" % (REF_NAME, MANIFEST))
    s7ref = None
    und = cwd / DEFORM_DIR / "undeformed"
    if _vasprun(und) is not None:                       # 只要本征值
        from pymatgen.io.vasp import Vasprun
        v = Vasprun(str(_vasprun(und)), parse_dos=False, parse_potcar_file=False)
        s7ref = object.__new__(Config)
        s7ref.name, s7ref.kpts = "%s/undeformed" % DEFORM_DIR, np.asarray(v.actual_kpoints, float)
        s7ref.eig = {str(getattr(sp, "name", sp)).lower(): np.asarray(a)[:, :, 0] for sp, a in v.eigenvalues.items()}
    else:
        print("[WARN] 没有 %s/undeformed/vasprun.xml：跳过 ivs-prim 与 S7 的带边核对" % DEFORM_DIR)

    warns, carriers = [], {}
    for r in dpt.get("results", []):
        car = r.get("carrier")
        valleys = ((r.get("multi_valley") or {}).get("valleys")) or []
        if not valleys:
            carriers[car] = {"status": "no_valleys",
                             "note": "S8.2 没有 %s 的多谷谷表（multi_valley.valleys）：本载流子不出谷间参数" % car}
            continue
        vals, note = carrier_params(car, valleys, cfg, man, s7ref, warns)
        if vals is None:
            carriers[car] = {"status": "no_k_star", "note": note}
            print("[..] " + note)
            continue
        grid = all(v.get("k_grid") for v in valleys)
        carriers[car] = {"status": "ok",
                         "valley_sig": {"use_grid": grid, "sig": json.loads(json.dumps(kc._valley_sig(valleys, grid)))},
                         "valleys": vals}
    params = {
        "schema": iv.SCHEMA, "skill_rev": _SKILL_REV, "created": _time.strftime("%Y-%m-%d %H:%M:%S"),
        "dim": dim, "scope": "K↔K′ 零阶（q = K 的区边界声子）；非 SOC；不含 K↔Q、K↔Γ 与一阶谷间耦合",
        "valleys_source": DPT_RESULT, "manifest": "%s/%s" % (SRC_DIR, MANIFEST),
        "M": man["M"], "det": man["det"], "star": man["star"], "members_frac": man["members_frac"],
        "disp_max_A": man.get("disp_max_A"), "groups": man["groups"], "modes": man["modes"],
        "s73_freq_rel_change": s73_rel, "warnings": warns, "carriers": carriers,
        "units": {"D_eV_A": "eV/Å（= 1e8 eV/cm），与 step.conf DP_OPTICAL_* 同单位", "hw_meV": "meV"},
        "definition": ("超胞（K、K′ 折叠到 Γ̃）冻结区边界模 ±u：K/K′ 二重态劈裂 ±|V|u；位移 Δ_j = u·√(M_s/m_j)·ẽ_j（超胞口径 = 原胞口径）；"
                       "D² = Σ_简并伙伴模 Σλ²/(g(g−1))，Σλ² 为组内去掉共同平移后的迹（与参考能级无关）；"
                       "dp_optical 的 X>X' 通道：1/τ = πD²/(Mω)·G_星·(g−1)/g（见 ivs_frozen.py）"),
        "reference": "Piscanec et al. PRL 93, 185503 (2004)（K 点冻声子劈裂）；Kaasbjerg et al. PRB 85, 115317 (2012)",
    }
    (out / "ivs_params.json").write_text(json.dumps(params, ensure_ascii=False, indent=1), encoding="utf-8")
    lines = ["# 冻声子谷间 IVS 参数（%s）" % params["created"],
             "# 只算 K↔K′ 零阶（非 SOC）；D 单位 eV/Å（= 1e8 eV/cm），ħω 单位 meV；超胞 M=%s；位移 %.3f Å"
             % (man["M"], float(man.get("disp_max_A") or 0))]
    for car, c in carriers.items():
        if c["status"] != "ok":
            lines.append("%s: %s" % (car, c.get("note", c["status"])))
            continue
        for v in c["valleys"]:
            for g in v["groups"]:
                lines.append("%-8s %-10s %s 带%-3d 组%-2d 模 %-7s ħω=%6.2f meV  D=%6.3f eV/Å  劈裂 %.1f meV"
                             % (car, v["label"], v["spin"], v["band_index"], g["group"],
                                ",".join(map(str, g["modes"])), g["hw_meV"], g["D_eV_A"], 1000 * g["shift_max_eV"]))
    for w in warns:
        lines.append("[WARN] " + w)
        print("[WARN] " + w)
    (out / "ivs_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[2:]))
    kc.invalidate_ivs_consumer(cwd, params)             # S8.4 用过的参数变了 / 刚可用 -> 让它重排
    print("[DONE] %s/ivs_params.json（S8.4 的 IVS_FROZEN=auto 会读它）" % OUTDIR_NAME)


if __name__ == "__main__":
    main()
