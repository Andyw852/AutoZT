#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step9c_odp.py —— [0035] 冻声子法 ODP 单点，扇出（step7c_odp_frozen，S7.2_odp）。

做什么：给 S8.4（2D）/ S8（3D，[0036]）的谷内光学形变势（ODP）自动出参数（D、ħω），不再依赖文献。
  1. 读 S5（step5_dielect，DFPT）的 Γ 点声子：频率、质量加权本征矢、Born 有效电荷；
  2. 选模：去掉 3 个声学模、虚频/太软的模；极性（红外活性）模默认不算（长程耦合在 POP 里，ODP_POLAR=include 才算）；
     简并的模整组一起算；
  3. 每个选中的模生成 ±u 两个位移构型：Δ_κ = u·√(M/m_κ)·e_κ，u 取到最大原子位移 = ODP_DISP_A（odp_frozen.py 有推导）；
     再加一个不位移的参考构型 odp-ref（简并带要用 u = 0 的本征值）；
  4. INCAR / KPOINTS / POTCAR 照抄 step7_deform/undeformed（与 S7/S8.2 同一套 k 网格、泛函、ISYM=0、LVHAR）——
     S8.2 谷表里的 k_S7 点在这个网格上，S8.21 直接在这些 k 点取本征值；LVHAR 给真空能级对齐。
产出：step7c_odp_frozen/{odp-ref, odp-mNNp, odp-mNNm, ...}（autozt 扇出各自 sbatch）+ odp_manifest.json（选模与位移记录）。
读取在 S8.21（gen_step12c_odp_read.py）。
"""
import hashlib
import json
import math
import os
import shutil
import sys
import time as _time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ke_common as kc  # noqa: E402
of = None                     # odp_frozen：main() 里切好 python 环境后再 import（它要 numpy）
import stepconf  # noqa: E402
from dim_common import require_dim, resolve_tpl  # noqa: E402

# =========================== 可改参数区（step.conf [params] 可覆盖） ===========================
OUTDIR_NAME = "step7c_odp_frozen"
PREV_CANDS = ["step1_opt", "step1_std_opt"]
DEFORM_DIR = "step7_deform"
DIELECT_DIR = "step5_dielect"
STEP_LABEL = "S7.2_odp"
SUB_PREFIX = "odp-"            # 与 skill.yaml 的 fanout "odp-*" 一致
REF_NAME = "odp-ref"
MANIFEST = "odp_manifest.json"
_SKILL_REV = "2026-10-13-0042"   # [0042] 面外极性判据只对有 σh 的 2D 材料（Janus/翘曲不剔同极类模）   # 2026-10-12-0036   # [0036] 三维也能跑（S8 的 ODP）；[0035-r3]
# 最大原子位移（Å）。0.01 Å 时强耦合的带边移动 ~10–30 meV（线性区），EDIFF=1E-8 下本征值噪声 ~1e-5 eV。
ODP_DISP_A = 0.01
# 极性（红外活性）模：exclude（出厂：长程 Fröhlich 耦合由 POP 管）/ include（把它们的短程形变势也算进 ODP）。
ODP_POLAR = "exclude"
ODP_POLAR_TOL = 0.1                  # 模有效电荷（e）大于它算极性（= odp_frozen.POLAR_TOL_E）
ODP_POLAR_Z_TOL = 0.02               # [0035-r3] 面外模有效电荷（e）大于它也算极性（面外极性，= odp_frozen.POLAR_Z_TOL_E）
ODP_FMIN_THZ = 0.5                   # 低于它的光学模不算（虚频/太软）（= odp_frozen.FMIN_THZ）
# auto = 全部合格的模；或者逗号分隔的模编号（按频率升序从 1 数，看 odp_manifest.json 的 modes），简并伙伴自动补全。
ODP_MODES = "auto"
# 单点数上限（1 + 2×模数）。大胞光学模多，超过就报错，确认后在 step.conf 调大。
ODP_MAX_JOBS = 40

SPEC = {
    "ODP_DISP_A": (None, "float"),
    "ODP_POLAR": (None, "str"),
    "ODP_POLAR_TOL": (None, "float"),
    "ODP_POLAR_Z_TOL": (None, "float"),
    "ODP_FMIN_THZ": (None, "float"),
    "ODP_MODES": (None, "str"),
    "ODP_MAX_JOBS": (None, "int"),
}
# =================================================================================================


def apply_conf(cwd):
    """step.conf 覆盖写进全局参数，返回 {键: 值}（写进 manifest）。"""
    global ODP_DISP_A, ODP_POLAR, ODP_POLAR_TOL, ODP_POLAR_Z_TOL, ODP_FMIN_THZ, ODP_MODES, ODP_MAX_JOBS
    got = {}
    if (Path(cwd) / stepconf.CONF_NAME).is_file():
        p = stepconf.load(SPEC, OUTDIR_NAME, str(cwd), strict=False)
        for k in SPEC:
            if p[k] is not None:
                got[k] = p[k]
                globals()[k] = p[k]
    ODP_POLAR = str(ODP_POLAR).strip().lower()
    if ODP_POLAR not in ("exclude", "include"):
        sys.exit("[ERROR] ODP_POLAR=%r 无效，只允许 exclude / include" % ODP_POLAR)
    if not (0.0 < float(ODP_DISP_A) <= 0.05):
        sys.exit("[ERROR] ODP_DISP_A=%r 不合理（应在 (0, 0.05] Å；太大出线性区，太小被本征值噪声淹没）" % ODP_DISP_A)
    return got


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _same_structure(a, b, label, tol=1e-4):
    """两份 pymatgen Structure：元素顺序、晶格、分数坐标（模 1）都一致，否则报错退出。"""
    import numpy as np
    if [str(x.specie) for x in a] != [str(x.specie) for x in b]:
        sys.exit("[ERROR] %s 的原子顺序/元素与 step1_opt CONTCAR 不同 —— 本步按原子序号对位移，必须一致" % label)
    dl = float(np.abs(np.asarray(a.lattice.matrix) - np.asarray(b.lattice.matrix)).max())
    df = np.asarray(a.frac_coords) - np.asarray(b.frac_coords)
    df -= np.round(df)
    dp = float(np.abs(df @ np.asarray(a.lattice.matrix)).max())
    if dl > tol or dp > tol:
        sys.exit("[ERROR] %s 的结构与 step1_opt CONTCAR 不同（晶格差 %.2e Å、位置差 %.2e Å）——"
                 "它是旧结构算的，先 rerun/retry 它再生成本步" % (label, dl, dp))


def read_dfpt(cwd):
    """S5 的 vasprun.xml/OUTCAR -> (频率 THz, 质量加权本征矢, Born 电荷 或 None, 结构)。
    换算与 step8.4_amset2d/amset2d_pop_freq.py 相同（它的 Δr 两条路线 5 位有效数字一致，验证过这套口径）。"""
    import numpy as np
    from pymatgen.io.vasp import Outcar, Vasprun
    d = Path(cwd) / DIELECT_DIR
    vr = next((d / n for n in ("vasprun.xml", "vasprun.xml.gz") if (d / n).is_file()), None)
    oc = next((d / n for n in ("OUTCAR", "OUTCAR.gz") if (d / n).is_file()), None)
    if vr is None or oc is None:
        sys.exit("[ERROR] 缺 %s/vasprun.xml 或 OUTCAR（S5 DFPT 没跑完？）" % DIELECT_DIR)
    v = Vasprun(str(vr), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
    if getattr(v, "normalmode_eigenvals", None) is None or getattr(v, "normalmode_eigenvecs", None) is None:
        sys.exit("[ERROR] %s/vasprun.xml 里没有 Γ 点声子（dynmat）—— S5 要 IBRION=7/8 的 DFPT" % DIELECT_DIR)
    ev = -np.asarray(v.normalmode_eigenvals)[::-1]
    vec = np.asarray(v.normalmode_eigenvecs)[::-1]
    f = np.sqrt(np.abs(ev)) * np.sign(ev)
    if int(str(v.vasp_version).split(".")[0]) < 6:
        f = f * 15.633302
    born = None
    try:
        o = Outcar(str(oc))
        o.read_lepsilon()
        if getattr(o, "born", None) is not None and len(o.born):
            born = np.asarray(o.born, float)
    except Exception as e:                                     # noqa: BLE001
        print("[WARN] 读 S5 的 Born 有效电荷失败（%s: %s）" % (type(e).__name__, e))
    return f, vec, born, v.final_structure


def select(groups, modes):
    """ODP_MODES 显式列表时：包含任一列出模的组整组算（简并伙伴补全），其余不算。返回实际要算的模 index 列表。"""
    spec = str(ODP_MODES or "auto").strip().lower()
    if spec == "auto":
        return [i for g in groups if g["use"] for i in g["modes"]]
    try:
        want = {int(x) for x in spec.replace(" ", "").split(",") if x}
    except ValueError:
        sys.exit("[ERROR] ODP_MODES=%r 应为 auto 或逗号分隔的模编号" % ODP_MODES)
    by = {m["index"]: m for m in modes}
    bad = sorted(i for i in want if i not in by or by[i]["acoustic"])
    if bad:
        sys.exit("[ERROR] ODP_MODES 里 %s 不是光学模编号（看 %s 的 modes）" % (bad, MANIFEST))
    out = []
    for g in groups:
        if want & set(g["modes"]):
            if set(g["modes"]) - want:
                print("[..] ODP_MODES：模 %s 与 %s 简并，一起算" % (sorted(want & set(g["modes"])), g["modes"]))
            if not g["use"]:
                print("[WARN] ODP_MODES 显式要了组 %s（%s），照算" % (g["modes"], g["why"]))
                g["use"], g["why"] = True, None
            out += g["modes"]
        else:
            g["use"], g["why"] = False, g["why"] or "ODP_MODES 没列"
    return out


def write_poscar(path, structure):
    from pymatgen.io.vasp import Poscar
    Poscar(structure).write_file(str(path), significant_figures=12)       # write_file 新旧 pymatgen 都有


def main():
    import numpy as np
    global of
    cwd = Path.cwd()
    kc.reexec_with_pymatgen(__file__)          # 先切环境（exec 重跑整个脚本），再做任何有副作用的事
    import odp_frozen as of                    # noqa: F811
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    # 本步重新生成 = 单点要重算：S8.21 读出的 odp_params.json 作废（S8.4 是否重排由 S8.21 比指纹决定）
    kc.invalidate_downstream(cwd, OUTDIR_NAME, "%s 重新生成（冻声子单点重算）" % OUTDIR_NAME)
    overrides = apply_conf(cwd)
    prev = kc.find_prev_dir(cwd, PREV_CANDS)
    if prev is None:
        sys.exit("[ERROR] 找不到含 CONTCAR 的上一步：%s" % PREV_CANDS)
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR", "step1_opt")
    dim = kc.read_method_dim(prev / kc.METHOD_FILE) or kc.resolve_dim_for(out / "POSCAR", "auto")[0]
    require_dim(dim, ("2d", "3d"), OUTDIR_NAME,
                why="冻声子 ODP 参数给 S8.4（2D）/ S8（3D，[0036] amset_dpo3d 插件）的 ODP 机制用；0D 没有能带")
    kc.write_method(out / kc.METHOD_FILE, dim, "冻声子 ODP 单点（扇出）")
    und = cwd / DEFORM_DIR / "undeformed"
    for n in ("INCAR", "KPOINTS", "POTCAR", "POSCAR"):
        if not (und / n).is_file():
            sys.exit("[ERROR] 缺 %s/undeformed/%s：先让 S7_deform 生成输入（本步照抄它的 INCAR/KPOINTS/POTCAR，"
                     "保证与 S7/S8.2 同一套 k 网格与参数）" % (DEFORM_DIR, n))
    from pymatgen.core import Structure

    base = Structure.from_file(str(out / "POSCAR"))
    _same_structure(base, Structure.from_file(str(und / "POSCAR")), "%s/undeformed/POSCAR" % DEFORM_DIR)
    freqs, eig, born, s5 = read_dfpt(cwd)
    _same_structure(base, s5, "%s（DFPT）" % DIELECT_DIR)
    masses = [float(s.specie.atomic_mass) for s in base]
    # [0042] 面外极性判据只对有水平镜面 σh 的 2D 材料成立：有 σh 时面外偶极模对谷内零阶耦合是对称性禁戒的；
    #   没有 σh（Janus、翘曲）时同极类 A1 模本身就带面外偶极，按面外极性剔掉会丢掉主导的谷内 ODP -> 不按面外判。
    z_tol, sigma_h = float("inf"), None
    if dim == "2d":
        z_tol = float(ODP_POLAR_Z_TOL)
        sigma_h = kc.za_coupling_check(cwd).get("sigma_h")
        if sigma_h is False:
            z_tol = float("inf")
            print("[..] 没有水平镜面 σh（Janus/翘曲）：面外偶极模不按面外极性剔除（谷内零阶耦合不被对称性禁戒），"
                  "只按总模有效电荷判极性")
        elif sigma_h is None:
            print("[WARN] 判不出有没有水平镜面 σh：仍按面外极性判据（面外模有效电荷 > %.2f e 剔除）；"
                  "Janus/翘曲材料请在 step.conf 写 ODP_POLAR_Z_TOL = 1e9 关掉它" % z_tol)
    try:
        modes, groups, notes = of.mode_table(freqs, eig, masses, born, polar_tol=float(ODP_POLAR_TOL),
                                             fmin_thz=float(ODP_FMIN_THZ), include_polar=(ODP_POLAR == "include"),
                                             # [0036] 面外极性只对 2D 有意义（层法向 z）；3D 只按总模有效电荷判
                                             # [0042] 2D 也只对有 σh 的材料判（见上）
                                             polar_z_tol=z_tol, two_d=(dim == "2d"))
    except ValueError as e:
        sys.exit("[ERROR] S5 的声子数据不能用：%s" % e)
    for n in notes:
        print("[WARN] %s" % n)
    use = select(groups, modes)
    by_index_ = {m["index"]: m for m in modes}
    print("[..] Γ 点光学模（频率升序编号；Z = 模有效电荷，总 / 面外%s，e）：" % ("" if dim == "2d" else "（3D 不用）"))
    for g in groups:
        zt = max((by_index_[i]["Z_mode_e"] or 0.0) for i in g["modes"]) if born is not None else None
        zz = max((by_index_[i]["Z_mode_z_e"] or 0.0) for i in g["modes"]) if born is not None else None
        print("     组 %d 模 %-8s %7.3f THz = %6.2f meV  Z=%s  %s%s" % (
            g["group"], ",".join(map(str, g["modes"])), g["freq_THz"], g["hw_meV"],
            "—" if zt is None else "%.3f/%.3f" % (zt, zz),
            ("%s极性 " % (by_index_[g["modes"][0]]["polar_kind"] or "")) if g["polar"] else "",
            "算" if g["use"] else "不算：%s" % g["why"]))
    if not use and born is None and ODP_POLAR != "include":
        sys.exit("[ERROR] 读不到 S5 的 Born 有效电荷（OUTCAR 里没有 BORN EFFECTIVE CHARGES？S5 要 LEPSILON=.TRUE.），"
                 "判不了哪些模是极性的：先查/重跑 S5；确认要把全部光学模都算上（含极性模的短程形变势）就在 S7.2 的 "
                 "step.conf 写 ODP_POLAR = include")
    if not use:
        sys.exit("[ERROR] 没有要算的光学模（全是声学/极性/软模）。本材料不需要冻声子 ODP：在 tf_*.yaml 关掉 odp_frozen，"
                 "或 step.conf 写 ODP_POLAR = include 把极性模的短程形变势也算上")
    njobs = 1 + 2 * len(use)
    if njobs > int(ODP_MAX_JOBS):
        sys.exit("[ERROR] 要 %d 个单点（%d 个模 × 2 + 参考），超过 ODP_MAX_JOBS=%d。确认后在 step.conf 调大 ODP_MAX_JOBS，"
                 "或用 ODP_MODES 只算一部分模" % (njobs, len(use), int(ODP_MAX_JOBS)))

    by_index = {m["index"]: m for m in modes}
    configs = [(REF_NAME, None, 0)]
    for i in use:
        configs += [("%sm%02dp" % (SUB_PREFIX, i), i, +1), ("%sm%02dm" % (SUB_PREFIX, i), i, -1)]
    # 上一轮生成过、这一轮不算的子目录：改名移出扇出范围（不删），免得 autozt 把它当成本步的一部分提交/判完成
    keep = {c[0] for c in configs}
    stamp = _time.strftime("%Y%m%d%H%M%S")
    for d in sorted(p for p in out.glob(SUB_PREFIX + "*") if p.is_dir()):
        if d.name not in keep:
            d.rename(out / ("_stale_%s_%s" % (d.name, stamp)))
            print("[..] %s 这一轮不算（选模变了），改名为 _stale_%s_%s" % (d.name, d.name, stamp))

    submit_body = resolve_tpl(Path(__file__).resolve().parent, "submit_std", dim).read_text(encoding="utf-8")
    incar_s7 = (und / "INCAR").read_text(encoding="utf-8", errors="ignore")
    rec_cfg = {}
    for name, idx, sign in configs:
        d = out / name
        d.mkdir(exist_ok=True)
        snap = kc.snapshot_inputs(d, ("INCAR", "POSCAR", "KPOINTS", "POTCAR"))
        if idx is None:
            st, u = base.copy(), 0.0
        else:
            delta, u = of.displacement(eig[by_index[idx]["_src"]], masses, float(ODP_DISP_A))
            st = Structure(base.lattice, base.species, np.asarray(base.cart_coords) + sign * delta,
                           coords_are_cartesian=True)
        write_poscar(d / "POSCAR", st)
        for n in ("KPOINTS", "POTCAR"):
            shutil.copyfile(und / n, d / n)
        lines = [("SYSTEM = %s odp %s" % (cwd.name, name)) if ln.strip().upper().startswith("SYSTEM") else ln
                 for ln in incar_s7.splitlines()]
        (d / "INCAR").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
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
        rec_cfg[name] = {"mode": idx, "sign": sign, "u_A": u}

    sel = []
    for i in use:
        m = by_index[i]
        _d, u = of.displacement(eig[m["_src"]], masses, float(ODP_DISP_A))
        sel.append({"index": i, "group": m["group"], "freq_THz": m["freq_THz"], "hw_meV": m["hw_meV"], "u_A": u,
                    "eigvec": np.round(eig[m["_src"]], 10).tolist(),
                    "max_atom_disp_A": float(np.max(np.linalg.norm(_d, axis=1)))})
    manifest = {
        "schema": "odp_manifest/1", "skill_rev": _SKILL_REV, "created": _time.strftime("%Y-%m-%d %H:%M:%S"),
        "disp_max_A": float(ODP_DISP_A), "polar": ODP_POLAR, "polar_tol_e": float(ODP_POLAR_TOL),
        "polar_z_tol_e": (z_tol if math.isfinite(z_tol) else None), "dim": dim, "sigma_h": sigma_h,
        "fmin_THz": float(ODP_FMIN_THZ), "overrides": overrides,
        "masses_amu": masses, "species": [str(s.specie) for s in base],
        "modes": [{k: v for k, v in m.items() if not k.startswith("_")} for m in modes],
        "groups": groups, "selected": sel, "configs": rec_cfg,
        "inputs_from": {"dir": "%s/undeformed" % DEFORM_DIR,
                        "sha256": {n: _sha(und / n) for n in ("INCAR", "KPOINTS", "POTCAR", "POSCAR")}},
        "s5_freq_THz": [float(x) for x in sorted(freqs)],
        "notes": notes,
    }
    (out / MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    print("[DONE] %s：%d 个单点（参考 + %d 个模 × ±%.3f Å），INCAR/KPOINTS/POTCAR 照抄 %s/undeformed，autozt 会各自提交"
          % (OUTDIR_NAME, njobs, len(use), float(ODP_DISP_A), DEFORM_DIR))


if __name__ == "__main__":
    main()
