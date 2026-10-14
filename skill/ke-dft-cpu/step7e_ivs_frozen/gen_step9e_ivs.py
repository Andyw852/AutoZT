#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step9e_ivs.py —— [0050] 谷间 IVS（K↔K′）冻声子第二步：区边界模 ±u 单点，扇出（step7e_ivs_frozen，S7.4_ivs）。

做什么：
  1. 读 S7.3（step7d_ivs_phonon）超胞 Γ̃ 点 DFPT 的频率与质量加权本征矢；按"原胞平移不变"子空间的投影把模分成
     原胞 Γ 来的（权重 ≈ 1，不要）与区边界折叠来的（≈ 0，K 点声子）；区边界模按频率简并分组；
  2. 每个区边界模生成 ±u 两个位移构型：Δ_j = u·√(M_s/m_j)·ẽ_j（超胞口径，u 取到最大原子位移 = IVS_DISP_A；
     超胞口径的 D 就是原胞口径的 D，推导见 ivs_frozen.py）；
  3. 参考构型 ivs-ref（超胞不位移）+ ivs-prim（原胞，Γ 中心网格含全部折叠 k，给 S8.22 标超胞 Γ̃ 能级的来源）。
  INCAR 照抄 step7_deform/undeformed 的单点设置（同一泛函/ENCUT/PREC、ISYM=0），NBANDS 按同一规则显式给
  （超胞与原胞折叠后空带数相当，见 ivs_frozen.nbands_for）；不要 LOCPOT（劈裂与参考能级无关）。
产出：step7e_ivs_frozen/{ivs-ref, ivs-prim, ivs-mNNp, ivs-mNNm, ...}（autozt 扇出各自提交）+ ivs_manifest.json。
读取在 S8.22（gen_step12d_ivs_read.py）。
"""
import hashlib
import json
import math
import shutil
import sys
import time as _time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import ke_common as kc  # noqa: E402
import stepconf  # noqa: E402
from dim_common import require_dim, resolve_tpl  # noqa: E402

iv = None

# =========================== 可改参数区（step.conf [params] 可覆盖） ===========================
OUTDIR_NAME = "step7e_ivs_frozen"
PH_DIR = "step7d_ivs_phonon"
DEFORM_DIR = "step7_deform"
STEP_LABEL = "S7.4_ivs"
SUB_PREFIX = "ivs-"                # 与 skill.yaml 的 fanout "ivs-*" 一致
REF_NAME = "ivs-ref"
PRIM_NAME = "ivs-prim"
MANIFEST = "ivs_manifest.json"
PH_MANIFEST = "ivs_supercell.json"
_SKILL_REV = "2026-10-14-0050"   # [0050] 谷间 IVS（K↔K′ 零阶）冻声子：区边界模 ±u 单点
IVS_DISP_A = 0.01                  # 最大原子位移（Å），与 S7.2 的 ODP_DISP_A 同一线性区
IVS_FMIN_THZ = 0.5                 # 区边界模频率低于它不算（虚频/太软）
IVS_MAX_JOBS = 40                  # 单点数上限（2 + 2×模数）：MoS₂ 类 3 原子原胞 = 18 个区边界模 -> 38
IVS_NBANDS_EXTRA = 4               # 每个原胞多给的空带（NBANDS = ⌈0.6·NELECT⌉ + 它×原胞数）
SPEC = {"IVS_DISP_A": (None, "float"), "IVS_FMIN_THZ": (None, "float"), "IVS_MAX_JOBS": (None, "int"),
        "IVS_NBANDS_EXTRA": (None, "int")}
SP_TAGS = {"IBRION": "-1", "NSW": "0", "ISYM": "0", "LVHAR": None, "LWAVE": ".FALSE.", "LCHARG": ".FALSE.",
           "LEPSILON": None, "ISTART": "0", "ICHARG": "2"}
# =================================================================================================


def apply_conf(cwd):
    got = {}
    if (Path(cwd) / stepconf.CONF_NAME).is_file():
        p = stepconf.load(SPEC, OUTDIR_NAME, str(cwd), strict=False)
        for k in SPEC:
            if p[k] is not None:
                got[k] = p[k]
                globals()[k] = p[k]
    if not (0.0 < float(IVS_DISP_A) <= 0.05):
        sys.exit("[ERROR] IVS_DISP_A=%r 不合理（应在 (0, 0.05] Å）" % IVS_DISP_A)
    return got


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_dfpt_sc(d):
    """S7.3 的 vasprun.xml -> (频率 THz, 质量加权本征矢 (3N, N, 3), 结构)。换算同 gen_step9c_odp.read_dfpt。"""
    import numpy as np
    from pymatgen.io.vasp import Vasprun
    vr = next((d / n for n in ("vasprun.xml", "vasprun.xml.gz") if (d / n).is_file()), None)
    if vr is None:
        sys.exit("[ERROR] 缺 %s/vasprun.xml（S7.3 超胞 DFPT 没跑完？）" % PH_DIR)
    v = Vasprun(str(vr), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
    if getattr(v, "normalmode_eigenvals", None) is None or getattr(v, "normalmode_eigenvecs", None) is None:
        sys.exit("[ERROR] %s/vasprun.xml 里没有 Γ̃ 点声子（dynmat）—— S7.3 要 IBRION=8 的 DFPT" % PH_DIR)
    ev = -np.asarray(v.normalmode_eigenvals)[::-1]
    vec = np.asarray(v.normalmode_eigenvecs)[::-1]
    f = np.sqrt(np.abs(ev)) * np.sign(ev)
    if int(str(v.vasp_version).split(".")[0]) < 6:
        f = f * 15.633302
    return f, vec, v.final_structure, vr


def supercell_from_manifest(man):
    import numpy as np
    s = man["supercell"]
    M3 = np.eye(3)
    M3[:2, :2] = np.asarray(man["M"], float)
    return {"lattice": np.asarray(s["lattice"], float), "frac": np.asarray(s["frac"], float),
            "species": list(s["species"]), "translations": [np.asarray(t, int) for t in s["translations"]],
            "origin": [(int(i), np.asarray(t, int)) for i, t in s["origin"]], "M3": M3, "det": int(man["det"])}


def prim_kmesh(folded, lat_prim, lat_sc, mesh_sc):
    """原胞网格：每个面内方向是折叠 k 分母的公倍数（全部折叠 k 都在网格上），密度不低于超胞网格；kz 照抄。"""
    import numpy as np
    den = 1
    for k in folded:
        for x in list(k)[:2]:
            q = Fraction(float(x)).limit_denominator(12).denominator
            den = den * q // math.gcd(den, q)
    rp = 2 * np.pi * np.linalg.inv(np.asarray(lat_prim, float)).T
    rs = 2 * np.pi * np.linalg.inv(np.asarray(lat_sc, float)).T
    s = min(np.linalg.norm(rs[i]) / max(int(mesh_sc[i]), 1) for i in range(2))
    out = []
    for i in range(2):
        n = max(1, int(math.ceil(np.linalg.norm(rp[i]) / s - 1e-6)))
        out.append(int(math.ceil(n / float(den))) * den)
    return out + [int(mesh_sc[2])], den


def write_kpoints(path, mesh, comment):
    Path(path).write_text("%s\n0\nGamma\n%d %d %d\n0 0 0\n" % (comment, mesh[0], mesh[1], mesh[2]),
                          encoding="utf-8", newline="\n")


def main():
    import numpy as np
    global iv
    cwd = Path.cwd()
    kc.reexec_with_pymatgen(__file__)
    import ivs_frozen as iv                    # noqa: F811
    from pymatgen.core import Element, Structure
    from pymatgen.io.vasp import Poscar
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    # 本步重新生成 = 单点要重算：S8.22 读出的 ivs_params.json 作废（S8.4 是否重排由 S8.22 比指纹决定）
    kc.invalidate_downstream(cwd, OUTDIR_NAME, "%s 重新生成（谷间冻声子单点重算）" % OUTDIR_NAME)
    overrides = apply_conf(cwd)
    ph = cwd / PH_DIR
    if not (ph / PH_MANIFEST).is_file():
        sys.exit("[ERROR] 缺 %s/%s：S7.3_ivsph 没生成" % (PH_DIR, PH_MANIFEST))
    phm = json.loads((ph / PH_MANIFEST).read_text(encoding="utf-8"))
    dim = kc.read_method_dim(cwd / "step1_opt" / kc.METHOD_FILE) or "2d"
    require_dim(dim, ("2d",), OUTDIR_NAME, why="谷间 IVS（K↔K′）冻声子只给 S8.4（二维）用")
    und = cwd / DEFORM_DIR / "undeformed"
    for n in ("INCAR", "POTCAR", "POSCAR"):
        if not (und / n).is_file():
            sys.exit("[ERROR] 缺 %s/undeformed/%s" % (DEFORM_DIR, n))
    if kc.incar_is_ncl(und / "INCAR") or kc.detect_ncl(cwd):
        sys.exit("[ERROR] 本材料是 SOC 能带：[0050] 谷间 IVS 冻声子只支持非 SOC（A 方案），tf_*.yaml 写 ivs_frozen: false")
    now = {n: _sha(und / n) for n in ("INCAR", "POTCAR", "POSCAR")}
    old = (phm.get("inputs_from") or {}).get("sha256") or {}
    diff = [n for n in now if old.get(n) != now[n]]
    if diff:
        sys.exit("[ERROR] S7 undeformed 的 %s 与生成 S7.3 时不同：S7 重新生成过，先 retry S7.3_ivsph（超胞 DFPT 重算）"
                 % "/".join(diff))
    sc = supercell_from_manifest(phm)
    freqs, eig, st_ph, vr = read_dfpt_sc(ph)
    if [str(x.specie) for x in st_ph] != sc["species"]:
        sys.exit("[ERROR] S7.3 DFPT 的原子顺序与 %s 不同：retry S7.3_ivsph" % PH_MANIFEST)
    dfr = np.asarray(st_ph.frac_coords) - sc["frac"]
    dfr -= np.round(dfr)
    if float(np.abs(dfr @ sc["lattice"]).max()) > 1e-4:
        sys.exit("[ERROR] S7.3 DFPT 的结构与 %s 记的超胞不同：retry S7.3_ivsph" % PH_MANIFEST)
    perms = iv.translation_perms(sc)
    try:
        modes, groups, notes = iv.zone_boundary_groups(freqs, eig, perms, fmin_thz=float(IVS_FMIN_THZ))
    except ValueError as e:
        sys.exit("[ERROR] S7.3 的声子数据不能用：%s" % e)
    for n in notes:
        print("[WARN] %s" % n)
    print("[..] 区边界（K 点折叠来的）模：")
    for g in groups:
        print("     组 %d 模 %-10s %7.3f THz = %6.2f meV  %s" % (
            g["group"], ",".join(map(str, g["modes"])), g["freq_THz"], g["hw_meV"],
            "算" if g["use"] else "不算：%s" % g["why"]))
    use = [i for g in groups if g["use"] for i in g["modes"]]
    if not use:
        sys.exit("[ERROR] 没有可算的区边界模（全是虚频/软模？查 S7.3 的频率）")
    njobs = 2 + 2 * len(use)
    if njobs > int(IVS_MAX_JOBS):
        sys.exit("[ERROR] 要 %d 个单点（参考 + 原胞 + %d 个模 × 2），超过 IVS_MAX_JOBS=%d：确认后在 step.conf 调大"
                 % (njobs, len(use), int(IVS_MAX_JOBS)))
    masses = [float(Element(s).atomic_mass) for s in sc["species"]]
    by_index = {m["index"]: m for m in modes}
    prim = Structure.from_file(str(und / "POSCAR"))
    counts_p = [sum(1 for s in prim if str(s.specie) == el) for el in dict.fromkeys(str(s.specie) for s in prim)]
    potcar = (und / "POTCAR").read_text(encoding="utf-8", errors="ignore")
    ne_p = iv.nelect_from_potcar(potcar, counts_p)
    nb_p = iv.nbands_for(ne_p, 1, int(IVS_NBANDS_EXTRA))
    nb_s = iv.nbands_for(ne_p * sc["det"], sc["det"], int(IVS_NBANDS_EXTRA))
    folded = [phm["folded"][n] for n in sorted(phm["folded"])]
    mesh_sc = [int(x) for x in phm["kmesh_sc"]]
    mesh_p, den = prim_kmesh(folded, phm["prim"]["lattice"], sc["lattice"], mesh_sc)

    configs = [(REF_NAME, None, 0), (PRIM_NAME, None, 0)]
    for i in use:
        configs += [("%sm%02dp" % (SUB_PREFIX, i), i, +1), ("%sm%02dm" % (SUB_PREFIX, i), i, -1)]
    keep = {c[0] for c in configs}
    stamp = _time.strftime("%Y%m%d%H%M%S")
    for d in sorted(p for p in out.glob(SUB_PREFIX + "*") if p.is_dir()):
        if d.name not in keep:
            d.rename(out / ("_stale_%s_%s" % (d.name, stamp)))
            print("[..] %s 这一轮不算（选模变了），改名为 _stale_%s_%s" % (d.name, d.name, stamp))
    submit_body = resolve_tpl(Path(__file__).resolve().parent, "submit_std", dim).read_text(encoding="utf-8")
    base_incar = kc.parse_incar((und / "INCAR").read_text(encoding="utf-8", errors="ignore"))
    st0 = Structure(sc["lattice"], sc["species"], sc["frac"])
    rec_cfg = {}
    for name, idx, sign in configs:
        d = out / name
        d.mkdir(exist_ok=True)
        snap = kc.snapshot_inputs(d, ("INCAR", "POSCAR", "KPOINTS", "POTCAR"))
        u = 0.0
        if name == PRIM_NAME:
            st, nb, mesh, kind = prim, nb_p, mesh_p, "prim"
        elif idx is None:
            st, nb, mesh, kind = st0, nb_s, mesh_sc, "ref"
        else:
            delta, u = iv.displacement(eig[by_index[idx]["_src"]], masses, float(IVS_DISP_A))
            st = Structure(st0.lattice, st0.species, np.asarray(st0.cart_coords) + sign * delta,
                           coords_are_cartesian=True)
            nb, mesh, kind = nb_s, mesh_sc, "disp"
        Poscar(st).write_file(str(d / "POSCAR"), significant_figures=12)
        write_kpoints(d / "KPOINTS", mesh, "IVS %s Gamma-centred" % kind)
        shutil.copyfile(und / "POTCAR", d / "POTCAR")
        tags = dict(SP_TAGS)
        tags["NBANDS"] = str(nb)
        inc = kc.merge_incar(base_incar, tags)
        (d / "INCAR").write_text(kc.incar_text(inc, system="%s ivs %s" % (cwd.name, name)), encoding="utf-8",
                                 newline="\n")
        _ic_log = []
        stepconf.apply_incar_file(d / "INCAR", log=_ic_log)
        for _m in _ic_log:
            print("[..] %s: %s" % (name, _m))
        changed = kc.inputs_changed(d, snap)
        if changed:
            n = kc.archive_stale_outputs(d, "stale-input-" + stamp)
            if n:
                print("[..] %s：输入变了（%s）-> 归档旧产物 %d 个，会重算" % (name, "/".join(changed), n))
        else:
            kc.archive_if_run_mismatch(d, label=name)
        (d / "submit.sh").write_text(submit_body.replace(
            "{{JOBNAME}}", "%s-ke-dft-cpu-%s-%s" % (cwd.name, STEP_LABEL, name)), encoding="utf-8", newline="\n")
        stepconf.apply_submit(d / "submit.sh", stepconf.read_submit(stepconf.CONF_NAME, used_incar=True))
        rec_cfg[name] = {"kind": kind, "mode": idx, "sign": sign, "u_A": u, "nbands": nb, "kmesh": mesh}

    sel = []
    for i in use:
        m = by_index[i]
        _d, u = iv.displacement(eig[m["_src"]], masses, float(IVS_DISP_A))
        sel.append({"index": i, "group": m["group"], "freq_THz": m["freq_THz"], "hw_meV": m["hw_meV"], "u_A": u,
                    "gamma_weight": m["gamma_weight"], "eigvec": np.round(eig[m["_src"]], 10).tolist(),
                    "max_atom_disp_A": float(np.max(np.linalg.norm(_d, axis=1)))})
    manifest = {
        "schema": "ivs_manifest/1", "skill_rev": _SKILL_REV, "created": _time.strftime("%Y-%m-%d %H:%M:%S"),
        "disp_max_A": float(IVS_DISP_A), "fmin_THz": float(IVS_FMIN_THZ), "overrides": overrides,
        "M": phm["M"], "det": phm["det"], "folded": phm["folded"], "star": phm["star"],
        "members_frac": phm["members_frac"], "masses_amu": masses, "species": sc["species"],
        "modes": [{k: v for k, v in m.items() if not k.startswith("_")} for m in modes],
        "groups": groups, "selected": sel, "configs": rec_cfg,
        "nbands": {"supercell": nb_s, "prim": nb_p, "nelect_prim": ne_p},
        "kmesh": {"supercell": mesh_sc, "prim": mesh_p, "fold_denominator": den},
        "s73": {"vasprun_sha256": _sha(vr), "freq_THz": [float(x) for x in sorted(freqs)]},
        "inputs_from": {"dir": "%s/undeformed" % DEFORM_DIR, "sha256": now},
        "notes": notes,
    }
    (out / MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    print("[DONE] %s：%d 个单点（超胞参考 + 原胞 + %d 个区边界模 × ±%.3f Å），超胞网格 %s、原胞网格 %s；autozt 会各自提交"
          % (OUTDIR_NAME, njobs, len(use), float(IVS_DISP_A), mesh_sc, mesh_p))


if __name__ == "__main__":
    main()
