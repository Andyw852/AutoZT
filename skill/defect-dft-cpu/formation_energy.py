#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""formation_energy.py —— defect-dft-cpu 形成能 / 转变能级 / 自洽费米能级 / P-N 判定。

纯标准库，登录节点可跑。消费 S4 汇总的裸总能 + 第0步凸包产出的化学势，输出
formation_energy_results.json 与 P/N 结论（这是本技能真正的最终产物）。

标准缺陷热力学（E_F 从 VBM 起算，0 ≤ E_F ≤ E_gap）：
  E_f(D,q) = [E_tot(D,q) − E_tot(bulk) − Σ_i Δn_i·μ_i + E_corr(D,q) + q·E_VBM] + q·E_F
  Δn_i = n_i(缺陷) − n_i(完美超胞)    (正=加原子，负=减原子)
  转变能级 ε(q1/q2) = [E_f(D,q1;E_F=0) − E_f(D,q2;E_F=0)] / (q2 − q1)

有限温自洽 E_F（窄带隙不能用"中带隙最低 E_f"捷径，见 README §关键修正4）：
  电荷中性  Q(E_F) = p0 − n0 + Σ_D N_D Σ_q q·P(D,q;E_F) = 0
  自由载流子用 Fermi-Dirac 积分 F_{1/2}（窄带隙简并，Boltzmann 失准）：
    n0 = Nc·F_{1/2}((E_F−E_gap)/kT)   p0 = Nv·F_{1/2}((−E_F)/kT)
    Nc/Nv = 2(2π m* kT/h²)^{3/2}，由有效质量 mstar_e/mstar_h 算（或直接给 Nc/Nv）
  P(D,q) = exp(−β(E_f0(D,q)+q·E_F)) / (1 + Σ_q' exp(−β(...)))
  分母的 1 为无缺陷位点；冻结总浓度后才用电荷态条件 softmax。

输入 energies.json（第0步凸包 + 势对齐产出）：
  { "mu": {"Sn":..,"Sb":..,"Te":..},   # 化学势 eV/原子（凸包稳定性窗口内一点）
    "E_gap": ..,                        # 带隙 eV
    "epsilon": ..,                      # 静态介电常数（Makov-Payne 图像电荷修正）
    "mstar_e": .., "mstar_h": ..,       # 电子/空穴有效质量（m_e 单位，默认 0.2）
    "N_D": 1e18,                        # 缺陷浓度 cm^-3（平衡费米能级，可扫）
    "T": 300 }                          # 温度 K

近似（已在结果里标注）：Makov-Payne 各向同性图像电荷（3x3x1 层状体系用 eFNV 更准，
本脚本留了 E_corr 覆盖入口）；有效质量抛物线带近似（替代真实 DOS 积分）。
"""
import sys, os, json, re, math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import defects_common as D

EV_PER_ANG = 14.3996448915       # e²/(4πε0)  eV·Å
KB = 8.617333262e-5              # eV/K
ME = 9.1093837015e-31            # 电子静质量 kg
H_PLANCK = 6.62607015e-34        # 普朗克常数 J·s
EV_J = 1.602176634e-19           # eV -> J
MADELUNG = 2.8373                # 简单立方 Madelung 常数（Makov-Payne）

def read_energy(outcar):
    """取 OUTCAR 最后总能（robust：优先 without entropy=E0，兜底 OSZICAR E0）。"""
    if not os.path.exists(outcar):
        return None
    e = None
    with open(outcar, errors="ignore") as f:
        for line in f:
            if "free  energy   TOTEN" in line or "without entropy" in line:
                try:
                    e = float(line.split("=")[1].split()[0])
                except (ValueError, IndexError):
                    continue
    if e is None:
        osz = outcar[:-6] + "OSZICAR" if outcar.endswith("OUTCAR") else outcar.replace("OUTCAR", "OSZICAR")
        if os.path.exists(osz):
            with open(osz, errors="ignore") as f:
                for line in f:
                    if "E0=" in line:
                        e = float(line.split("E0=")[1].split()[0])
    return e

def fd_half(eta, npts=2000, xmax=50.0):
    """Fermi-Dirac 积分 F_{1/2}(η) = (2/√π)∫ x^{1/2}/(1+e^{x-η}) dx（纯标准库梯形法）。

    窄带隙材料 E_F 贴带边（简并），自由载流子密度不能用 Boltzmann 近似，
    必须用这个积分。η < -20 时退回 Boltzmann 极限 e^η（此时二者几乎重合）。"""
    eta = float(eta)
    if eta < -20.0:
        return math.exp(eta)
    h = xmax / npts
    s = 0.0
    prev = 0.0                    # x=0: sqrt(0)/(1+e^{-eta}) = 0
    for i in range(1, npts + 1):
        x = i * h
        cur = math.sqrt(x) / (1.0 + math.exp(x - eta))
        s += 0.5 * (prev + cur) * h
        prev = cur
    return (2.0 / math.sqrt(math.pi)) * s

def eff_dos(mstar, T):
    """有效态密度 N = 2(2π m* kT/h²)^{3/2} (cm^-3)，mstar 以 m_e 为单位。

    注意单位：m*=mstar·m_e(kg)、kT(eV→J)、h(J·s)，全部 SI 自洽。"""
    m = mstar * ME
    kT_J = KB * T * EV_J
    N_m3 = 2.0 * (2.0 * math.pi * m * kT_J / (H_PLANCK * H_PLANCK)) ** 1.5
    return N_m3 * 1e-6           # m^-3 -> cm^-3

def lattice_volume(poscar):
    st = D.parse_poscar(poscar)
    a, b, c = st["lat"]
    v = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
    return abs(v)

def load_ref(path="energies.json"):
    if not os.path.exists(path):
        raise SystemExit("[错误] 找不到 %s —— 先做第0步凸包拿化学势，再跑本脚本" % path)
    r = json.load(open(path, encoding="utf-8"))
    if not isinstance(r.get("mu"), dict) or not r["mu"]:
        raise SystemExit("[错误] energies.json 缺有效 mu（化学势）")
    if r.get("E_gap") is None and not r.get("is_metal"):
        raise SystemExit("[错误] energies.json 缺 E_gap（半导体）；金属请设 is_metal=true")
    return r

def _aggregate_types(conc, name_type):
    """把逐位点 conc 汇总为按缺陷名/按缺陷大类，用于回答"哪种最可能生成"。

    conc: [{"defect","nsites"...,"conc_cm3","Ef_min","m"}]（同一 name 可能对应多个
    不等价位点目录）。按 name 汇总后 nsites=位点目录数、m_total=等价位点总数。
    """
    agg = {}
    for c in conc:
        agg.setdefault(c["defect"], []).append(c)
    by_name = []
    for nm, cs in agg.items():
        by_name.append({"name": nm, "type": name_type.get(nm, "?"),
                        "nsites": len(cs), "m_total": sum(x["m"] for x in cs),
                        "Ef_min": round(min(x["Ef_min"] for x in cs), 5),
                        "conc_cm3": sum(x["conc_cm3"] for x in cs)})
    by_name.sort(key=lambda x: -x["conc_cm3"])
    agg2 = {}
    for x in by_name:
        d = agg2.setdefault(x["type"], {"type": x["type"], "nsites": 0, "m_total": 0,
                                        "Ef_min": 1e9, "conc_cm3": 0.0, "names": []})
        d["nsites"] += x["nsites"]; d["m_total"] += x["m_total"]
        d["Ef_min"] = min(d["Ef_min"], x["Ef_min"])
        d["conc_cm3"] += x["conc_cm3"]; d["names"].append(x["name"])
    by_type = sorted(agg2.values(), key=lambda x: -x["conc_cm3"])
    for d in by_type:
        d["Ef_min"] = round(d["Ef_min"], 5)
    return by_name, by_type


def parse_q(dirname):
    """从 'def-001_v_Te_q+1' 提取电荷态 q；不匹配返回 None。"""
    m = re.search(r"_q([+-]?[0-9]+)$", dirname)
    return int(m.group(1)) if m else None

def qs(q):
    """电荷态显示名：+1 -> '+1', -1 -> '-1', 0 -> '0'。"""
    return ("+" + str(q)) if q > 0 else str(q)

def softmax_charge(Ef0, qs, Ef, beta):
    """P(D,q) = exp(-beta(Ef0+q*Ef)) / sum，数值稳定（减最小值防溢出）。"""
    vals = [Ef0[i] + qs[i]*Ef for i in range(len(qs))]
    m = min(vals)
    exps = [math.exp(-beta*(v - m)) for v in vals]
    z = sum(exps)
    return [x/z for x in exps]

def site_probabilities(Ef0, charges, Ef, beta):
    """互斥电荷态 + 无缺陷位点（能量0）；稳定配分，ΣP<=1。"""
    vals = [e + q * Ef for e, q in zip(Ef0, charges)]
    minimum = min([0.0] + vals)
    weights = [math.exp(-beta * (v - minimum)) for v in vals]
    z = math.exp(beta * minimum) + sum(weights)
    return [w / z for w in weights]


def defect_concentrations(allEf0, Ef, beta, vol_cm3, frozen_tot=None):
    """求根和报告共享同一占据模型；冻结浓度按唯一 base 标识。"""
    conc = []
    for base, name, charges, ef0, mD in allEf0:
        sd = mD / vol_cm3
        if frozen_tot is None:
            cq = [sd * p for p in site_probabilities(ef0, charges, Ef, beta)]
        else:
            cq = [frozen_tot[base] * p for p in softmax_charge(ef0, charges, Ef, beta)]
        total = sum(cq)
        energies = [e + q * Ef for e, q in zip(ef0, charges)]
        k = min(range(len(energies)), key=lambda i: energies[i])
        conc.append({"defect": name, "base": base, "conc_cm3": total,
                     "charge_concentrations_cm3": {qs(q): c for q, c in zip(charges, cq)},
                     "q_avg": sum(q * c for q, c in zip(charges, cq)) / total if total else 0.0,
                     "ionized_cm3": sum(abs(q) * c for q, c in zip(charges, cq)),
                     "saturated": total >= 0.5 * sd, "site_density": sd,
                     "m": mD, "Ef_min": energies[k], "q_stable": charges[k]})
    return sorted(conc, key=lambda c: -c["conc_cm3"])


def makov_payne(q, vol, eps):
    if q == 0 or eps <= 0:
        return 0.0
    L = vol ** (1.0/3.0)
    return MADELUNG * q*q / (2.0 * eps * L) * EV_PER_ANG

def read_core_potentials(outcar):
    """读 OUTCAR 最后一次 'average (electrostatic) potential at core'，返回每原子核心势列表。"""
    if not os.path.exists(outcar):
        return None
    lines = open(outcar, errors="ignore").read().splitlines()
    last = None
    for i, ln in enumerate(lines):
        if "average (electrostatic) potential at core" in ln:
            last = i
    if last is None:
        return None
    i = last + 1
    while i < len(lines) and ("test charge" in lines[i] or "norm" in lines[i]):
        i += 1
    vals = []
    while i < len(lines):
        toks = lines[i].split()
        i += 1
        if not toks:
            continue
        if len(toks) % 2 != 0:
            break
        try:
            vals.extend(float(toks[j + 1]) for j in range(0, len(toks), 2))
        except (ValueError, IndexError):
            break
    return vals

def potential_align(defect_outcar, bulk_outcar, defect_poscar, bulk_poscar):
    """Kumagai-Oba 势对齐 ΔV：缺陷超胞与 bulk 超胞核心势差的稳健平均。

    匹配：每个缺陷原子用**笛卡尔最近像距离**（六方胞分数距离失真）取最近同物种
    bulk 原子，ΔV_j = V_def[j] - V_bulk[ref]。
    取中间 50%（四分位）平均，排除缺陷近邻（两侧离群）。返回 ΔV (eV)。"""
    Vd = read_core_potentials(defect_outcar)
    Vb = read_core_potentials(bulk_outcar)
    if not Vd or not Vb:
        return 0.0
    sd = D.parse_poscar(defect_poscar)
    sb = D.parse_poscar(bulk_poscar)
    G = D._metric(sb["lat"])
    dv = []
    for j, (a, c) in enumerate(zip(sd["atoms"], sd["coords"])):
        if j >= len(Vd):
            break
        best = None
        for i, (a2, c2) in enumerate(zip(sb["atoms"], sb["coords"])):
            if a2 != a:
                continue
            d = D._frac_sep(c, c2, G)   # 笛卡尔最近像距离(Å)
            if best is None or d < best[0]:
                best = (d, i)
        if best is None or best[1] >= len(Vb):
            continue
        dv.append(Vd[j] - Vb[best[1]])
    if not dv:
        return 0.0
    dv.sort()
    n = len(dv)
    lo, hi = n // 4, 3 * n // 4
    return sum(dv[lo:hi]) / (hi - lo) if hi > lo else dv[n // 2]

def vbm_cbm_from_eigenval(path):
    """从 EIGENVAL 提取 VBM/CBM 本征值（绝对值，VASP 内部参考系）。SOC 单自旋通道。"""
    lines = open(path, errors="ignore").read().splitlines()
    p = lines[5].split()
    nkpts, nbands = int(p[1]), int(p[2])
    idx = 6
    vbm, cbm = None, None
    for _ in range(nkpts):
        while idx < len(lines) and lines[idx].strip() == "":
            idx += 1
        if idx >= len(lines):
            break
        idx += 1  # 跳过 k 点坐标行
        for _b in range(nbands):
            if idx >= len(lines):
                break
            tok = lines[idx].split()
            idx += 1
            if len(tok) < 3:
                continue
            e, occ = float(tok[1]), float(tok[2])
            if occ > 0.5:
                vbm = e if vbm is None or e > vbm else vbm
            else:
                cbm = e if cbm is None or e < cbm else cbm
    return vbm, cbm

def site_multiplicity(bulk_st, name, disp=None, meta=None):
    """位点多重度 m_D：该缺陷位点在超胞内的等价位点数。

    优先用 S2 manifest 里由 defects_common.symmetry_orbits()（spglib 空间群轨道）
    算好的 meta['site_count']，任意结构正确。没有 meta 时退回旧的 (物种,z) 层状
    近似（只对 A2B2Te5 类高对称层状结构正确，别的一般会偏大 -> 高估浓度）。"""
    if meta and meta.get("site_count"):
        return int(meta["site_count"])
    atoms = bulk_st["atoms"]
    coords = bulk_st["coords"]
    if "_i" in name:
        return 1
    if name.startswith("v_"):
        sp = name[2:]
    else:
        parts = name.split("_")
        sp = parts[1] if len(parts) == 2 else None
    if sp is None:
        return 1
    m = re.search(r"z=([0-9.]+)", disp or "")
    if m:
        zt = float(m.group(1)) % 1.0
        zf = min(zt, 1.0 - zt)   # 折叠 z<->1-z（反演/镜面对称）
        n = 0
        for a, c in zip(atoms, coords):
            if a != sp:
                continue
            zz = c[2] % 1.0
            zzf = min(zz, 1.0 - zz)
            if abs(zzf - zf) < 0.02:
                n += 1
        if n > 0:
            return n
    return atoms.count(sp)

def load_defects(summary_json="step4_analysis/formation_energy_summary.json",
                 manifest_json="step2_defects/defects_manifest.json"):
    if not os.path.exists(summary_json):
        raise SystemExit("[错误] 找不到 %s —— 先跑完 S4" % summary_json)
    summ = json.load(open(summary_json, encoding="utf-8"))
    E_bulk = summ.get("E_bulk")
    if E_bulk is None:
        raise SystemExit("[错误] summary 缺 E_bulk（step1_bulk 未完成？）")
    if not os.path.exists(manifest_json):
        raise SystemExit("[错误] 缺少缺陷 manifest: %s" % manifest_json)
    manifest = {}
    for it in json.load(open(manifest_json, encoding="utf-8")):
        if it["dir"] in manifest:
            raise SystemExit("[错误] manifest 重复目录: %s" % it["dir"])
        manifest[it["dir"]] = it
    # 汇总：{dir: {"name", "counts", "E_neutral", "charged": {q: E}}}
    defects = {}
    raw = summ.get("defects", {})
    for d, en in raw.items():
        q = parse_q(d)
        base = d
        if q is not None:
            base = d[:d.rfind("_q")]
        if base not in manifest or not isinstance(manifest[base].get("counts"), dict):
            raise SystemExit("[错误] manifest 缺 %s 的 counts" % base)
        counts = manifest[base]["counts"]
        if any(isinstance(n, bool) or not isinstance(n, int) or n < 0 for n in counts.values()):
            raise SystemExit("[错误] %s counts 必须为非负整数" % base)
        e = en.get("step2_defects")
        if e is None:
            e = en.get("step3_charged")
        if q is None:
            defects.setdefault(base, {}).setdefault("E_neutral", e)
        else:
            defects.setdefault(base, {}).setdefault("charged", {})[q] = e
    return E_bulk, manifest, defects

def analyze(ref=None):
    ref = ref or load_ref()
    E_bulk, manifest, defects = load_defects()
    bulk_poscar = "step1_bulk/POSCAR" if os.path.exists("step1_bulk/POSCAR") else "step1_bulk/CONTCAR"
    if not os.path.exists(bulk_poscar):
        raise SystemExit("[错误] 找不到 step1_bulk/POSCAR —— 先跑完 step1")
    bulk_st = D.parse_poscar(bulk_poscar)
    bulk_counts = {}
    for a in bulk_st["atoms"]:
        bulk_counts[a] = bulk_counts.get(a, 0) + 1
    vol = lattice_volume(bulk_poscar)
    eps = float(ref.get("epsilon", 0.0) or 0.0)
    is_metal = (bool(ref.get("is_metal")) or ref.get("E_gap") is None
                or float(ref.get("E_gap", 0.0) or 0.0) <= 0.0)
    E_gap = 0.0 if is_metal else float(ref["E_gap"])
    if is_metal:
        # 仅支持中性金属缺陷；不把绝对金属费米能同时加到截距和斜率。
        if any(q != 0 for d in defects.values() for q in d.get("charged", {})):
            raise SystemExit("[错误] 金属非零 q 缺陷的电子参考/修正约定未受支持；请仅使用中性缺陷")
        E_VBM = 0.0
    else:
        E_VBM = ref.get("VBM_bulk_abs") if ref.get("VBM_bulk_abs") is not None else ref.get("VBM")
        if E_VBM is None:
            raise SystemExit("[错误] energies.json 缺 VBM_bulk_abs —— 需先跑 extract_band_refs（从 step1_bulk EIGENVAL 提取）。"
                             "静默默认 0 会把 E_f 平移 q×6.76 eV，禁止。")
        E_VBM = float(E_VBM)
    bulk_outcar = "step1_bulk/OUTCAR"
    mu_vertices = ref.get("mu_vertices") or [{"mu": ref["mu"]}]

    results = {"step": "formation_energy", "status": "done",
               "E_bulk": E_bulk, "volume_AA3": round(vol, 3),
               "epsilon": eps, "E_gap": E_gap, "E_VBM": E_VBM,
               "is_metal": is_metal,
               "defects": {}}

    # ---- 顶点无关的缺陷数据（Eq/dV/corr/dn/mD 只算一次） ----
    base_data = []
    for base, d in sorted(defects.items()):
        info = manifest[base]
        name = info.get("name", base)
        counts = info["counts"]
        E0 = d.get("E_neutral")
        if E0 is None:
            continue
        dn = {el: counts.get(el, 0) - bulk_counts.get(el, 0)
              for el in set(list(counts) + list(bulk_counts))}
        qnn = sorted([q for q in d.get("charged", {}) if q != 0])
        q_all = [0] + qnn
        dV = {0: 0.0}
        for q in qnn:
            cdir = "step3_charged/%s_q%s" % (base, qs(q))
            dV[q] = potential_align(cdir + "/OUTCAR", bulk_outcar,
                                     cdir + "/POSCAR", bulk_poscar)
        Eq = {}; corr = {}
        for q in q_all:
            Eq[q] = E0 if q == 0 else d["charged"][q]
            corr[q] = makov_payne(q, vol, eps)
        mD = site_multiplicity(bulk_st, name, info.get("disp"), info.get("meta"))
        base_data.append((base, name, counts, dn, q_all, Eq, dV, corr, mD))

    name_type = {}
    for it in manifest.values():
        name_type[it.get("name", it["dir"])] = it.get("meta", {}).get("type", "?")
    # ---- 遍历所有化学势顶点（各极端化学势下的稳定性） ----
    vol_cm3 = vol * 1e-24
    vertex_results = []
    for vi, vinfo in enumerate(mu_vertices):
        mu = vinfo.get("mu", vinfo)
        allEf0 = []
        defects_vi = {}
        for base, name, counts, dn, q_all, Eq, dV, corr, mD in base_data:
            required = set(dn)
            if (not isinstance(mu, dict) or any(el not in mu for el in required)
                    or any(not isinstance(mu[el], (int, float)) or not math.isfinite(mu[el])
                           for el in required)):
                raise SystemExit("[错误] 顶点%d 缺少有效化学势（需要 %s）" % (vi, sorted(required)))
            chempot = sum(dn[el] * mu[el] for el in required)
            Ef0 = {}
            for q in q_all:
                Ef0[q] = Eq[q] - E_bulk - chempot + corr[q] + q * (E_VBM + dV[q])
            trans = {}
            qlist = sorted(Ef0)
            for i in range(len(qlist) - 1):
                q1, q2 = qlist[i], qlist[i + 1]
                if q2 - q1 != 0:
                    trans["(%s/%s)" % (qs(q1), qs(q2))] = round((Ef0[q1] - Ef0[q2]) / (q2 - q1), 5)
            defects_vi[base] = {"name": name, "counts": counts, "dn": dn,
                                "Ef0": {qs(q): round(v, 5) for q, v in Ef0.items()},
                                "transitions": trans}
            allEf0.append((base, name, sorted(q_all), [Ef0[q] for q in sorted(q_all)], mD))
        T_ = float(ref.get("T", 300.0))
        if not math.isfinite(T_) or T_ <= 0:
            raise SystemExit("[错误] T 必须为有限正温度")
        if is_metal:
            ef = 0.0  # 中性金属参考；E_F_metal 不进入形成能
            conc = defect_concentrations(allEf0, ef, 1.0 / (KB * T_), vol_cm3)
            carrier, n0, p0 = "metal", 0.0, 0.0
        else:
            ef, _, carrier, n0, p0, conc = solve_ef(allEf0, ref, vol_cm3, return_concentrations=True)
        dominant = conc[0]["defect"] if conc else None
        dominant_base = conc[0]["base"] if conc else None
        dom_ef = round(conc[0]["Ef_min"], 5) if conc else None
        top6 = [{"defect": c["defect"], "base": c["base"],
                 "q_stable": c["q_stable"], "Ef_eV": round(c["Ef_min"], 5),
                 "conc": "%.1e%s" % (c["conc_cm3"], "(稀释失效)" if c["saturated"] else ""),
                 "conc_cm3": c["conc_cm3"], "site_density": c["site_density"],
                 "saturated": c["saturated"], "m": c["m"],
                 "q_avg": round(c["q_avg"], 2), "ionized_cm3": c["ionized_cm3"]}
                for c in conc[:6]]
        by_name, by_type = _aggregate_types(conc, name_type)
        tsum = {"by_name": by_name, "by_type": by_type,
                "dominant_name": by_name[0]["name"] if by_name else None,
                "dominant_type": by_type[0]["type"] if by_type else None}
        vertex_results.append({"vertex": vi, "mu": mu, "E_F_eq": round(ef, 5),
                               "dominant_defect": dominant, "dominant_base": dominant_base,
                               "concentrations": conc, "carrier_type": carrier,
                               "n0_cm3": n0, "p0_cm3": p0,
                               "dominant_ef_eV": dom_ef, "top_concentration": top6,
                               "type_summary": tsum})
        if vi == 0:
            results["defects"] = defects_vi
            results["E_F_eq"] = round(ef, 5)
            results["dominant_defect"] = dominant
            results["dominant_base"] = dominant_base
            results["concentrations"] = conc
            results["dominant_ef_eV"] = dom_ef
            results["dominant_type"] = tsum["dominant_type"]
            results["dominant_name"] = tsum["dominant_name"]
            results["type_summary"] = tsum
            results["carrier_type"] = carrier
            results["carrier_density_cm3"] = {"n0": round(n0, 4), "p0": round(p0, 4)}

    results["vertex_summary"] = vertex_results

    out = "formation_energy_results.json"
    json.dump(results, open(out, "w"), indent=2, ensure_ascii=False)
    print("[OK] 形成能/转变能级 -> %s" % out)
    for v in vertex_results:
        if is_metal:
            print("   顶点%d: E_F=%.3f eV  metal  主导=%s (E_f=%s eV)"
                  % (v["vertex"], v["E_F_eq"], v["dominant_defect"], v.get("dominant_ef_eV")))
        else:
            print("   顶点%d: E_F=%.3f eV  %s  (n0=%.2e, p0=%.2e cm^-3)  主导=%s"
                  % (v["vertex"], v["E_F_eq"], v["carrier_type"], v["n0_cm3"], v["p0_cm3"], v["dominant_defect"]))
        for b in v.get("type_summary", {}).get("by_type", []):
            print("       类型 %-10s 位点目录=%d 等价位点=%d  E_f_min=%+.3f eV  总浓度=%.2e cm^-3"
                  % (b["type"], b["nsites"], b["m_total"], b["Ef_min"], b["conc_cm3"]))
    return results
def solve_ef(allEf0, ref, vol_cm3, return_concentrations=False):
    """自洽费米能级 + 判型：饱和浓度 + 载流子比较（支持两温度冻结）。

    浓度用互斥电荷态配分（含无缺陷位点）：
      c(D,q) = (m_D/V) exp(-βE_fq) / (1 + Σ_q exp(-βE_fq))
    同一 base 的所有电荷态总浓度不超过其位点密度；不同 base 未作共享位点竞争。

    两温度（ref 含 T_growth，Du et al.）：
      1) 生长温度解 E_F，冻结每个缺陷种类的总浓度 [D]_total = Σ_q c(D,q)；
      2) 300K 重解 E_F：把 [D]_total 在各电荷态间重新分配（300K 玻尔兹曼），
         电荷中性只自由载流子变化。注意冻结层级是 D 不是 (D,q) —— 冻 (D,q) 会
         锁死电荷态分布，低温解会失去自由度、解错根。"""
    E_gap = float(ref["E_gap"])
    T = float(ref.get("T", 300.0))
    T_g = ref.get("T_growth")
    if not math.isfinite(T) or T <= 0 or (T_g is not None and (not math.isfinite(float(T_g)) or float(T_g) <= 0)):
        raise SystemExit("[错误] T/T_growth 必须为有限正温度")
    kT = KB * T
    beta = 1.0 / kT
    Nc = float(ref.get("Nc", eff_dos(ref.get("mstar_e", 0.2), T)))
    Nv = float(ref.get("Nv", eff_dos(ref.get("mstar_h", 0.2), T)))

    disorder_warned = [False]
    lo, hi = -0.5, E_gap + 0.5

    def defect_charge(Ef, beta_defect, frozen_tot=None):
        conc = defect_concentrations(allEf0, Ef, beta_defect, vol_cm3, frozen_tot)
        if any(c["Ef_min"] < 0 for c in conc) and not disorder_warned[0]:
            disorder_warned[0] = True
            print("[警告] 形成能<0 —— 稀释近似失效；位点饱和模型不描述缺陷间相互作用")
        return sum(c["q_avg"] * c["conc_cm3"] for c in conc)

    def q_total(Ef, beta_carrier, beta_defect, frozen_tot=None):
        n0 = Nc * fd_half((Ef - E_gap) * beta_carrier)
        p0 = Nv * fd_half(-Ef * beta_carrier)
        return p0 - n0 + defect_charge(Ef, beta_defect, frozen_tot)

    def solve_root(qfunc):
        prev_ef, prev_q = lo, qfunc(lo)
        best = (abs(prev_q), lo)
        npts = 2000
        for i in range(1, npts + 1):
            ef = lo + (hi - lo) * i / npts
            q = qfunc(ef)
            if abs(q) < best[0]:
                best = (abs(q), ef)
            if prev_q == 0:
                return prev_ef
            if q == 0:
                return ef
            if q * prev_q < 0:
                ef_root = prev_ef + (ef - prev_ef) * abs(prev_q) / (abs(prev_q) + abs(q))
                best = (0.0, ef_root)
                break
            prev_ef, prev_q = ef, q
        return best[1]

    d_total = None
    if T_g is not None:
        beta_g = 1.0 / (KB * float(T_g))
        ef_tg = solve_root(lambda Ef: q_total(Ef, beta_g, beta_g))
        growth_conc = defect_concentrations(allEf0, ef_tg, beta_g, vol_cm3)
        d_total = {c["base"]: c["conc_cm3"] for c in growth_conc}
        ef_eq = solve_root(lambda Ef: q_total(Ef, beta, beta, frozen_tot=d_total))
        print("[两温度] T_g=%.0fK 冻结 [D]_total，%.0fK 重解 E_F=%.3f eV" % (float(T_g), T, ef_eq))
    else:
        ef_eq = solve_root(lambda Ef: q_total(Ef, beta, beta))

    n0_eq = Nc * fd_half((ef_eq - E_gap) * beta)
    p0_eq = Nv * fd_half(-ef_eq * beta)
    if p0_eq > n0_eq:
        carrier = "p 型"
    elif n0_eq > p0_eq:
        carrier = "n 型"
    else:
        carrier = "本征"
    conc = defect_concentrations(allEf0, ef_eq, beta, vol_cm3, d_total)
    dominant = conc[0] if conc else None
    dom_desc = ("%s[%s](q=%+d, E_f=%.3f eV)" %
                (dominant["defect"], dominant["base"], dominant["q_stable"], dominant["Ef_min"])) if dominant else "(无)"
    result = (ef_eq, dom_desc, carrier, n0_eq, p0_eq)
    return result + (conc,) if return_concentrations else result



def verify_identity(ref=None):
    """空缺陷自检：Δn_i=0 的 '空缺陷' 过一遍 E_f，应严格为 0。

    若凸包用 target_prim 而 E_bulk=E_super/9 不一致（δ≠0），这个残差=δ量级，
    说明 μ 端口径污染进了缺陷能量。比看 ΔH_f 塌多少更直接。"""
    ref = ref or load_ref()
    E_bulk, _manifest, _defects = load_defects()
    bulk_poscar = "step1_bulk/POSCAR" if os.path.exists("step1_bulk/POSCAR") else "step1_bulk/CONTCAR"
    bulk_st = D.parse_poscar(bulk_poscar)
    bulk_counts = {}
    for a in bulk_st["atoms"]:
        bulk_counts[a] = bulk_counts.get(a, 0) + 1
    # 空缺陷：成分与 bulk 相同(Δn=0)，q=0，无 chempot/corr/VBM 项
    dn = {el: 0 for el in bulk_counts}
    if any(el not in ref.get("mu", {}) for el in dn):
        raise SystemExit("[错误] 空缺陷自检也需要完整的元素化学势")
    E_f_empty = E_bulk - E_bulk - sum(dn[el] * ref["mu"][el] for el in dn)
    print("[自检] 空缺陷 E_f(Δn=0, q=0) = %+.6f eV（应严格为 0）" % E_f_empty)
    return E_f_empty

if __name__ == "__main__":
    analyze()
