#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_kappa.py -- lattice thermal conductivity (step2_kappa).

The last step of fit-fc-thermal: take the fc2/fc3 that step1_fit produced and run
the phono3py three-phonon BTE, writing kappa_summary.json (the step marker).

gen only prepares the submitted job: it validates the step1_fit artifacts, copies
fc2/fc3 + a phono3py YAML into step2_kappa/, works out the BTE mesh, and renders
submit.sh.  The compute node runs kappa_driver.py.

产物：step2_kappa/{kappa_summary.json, kappa-m<mesh>.hdf5, queue.out, ...}
"""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# 脱离 autozt 直接运行（python skill/fit-fc-thermal/gen_*.py）时，gen_need 里
# 来自 skill/_common 的 stepconf/dim_common/thickness_2d 等不在脚本旁边——
# 追加到 sys.path 末尾兜底（autozt 推送的同目录拷贝仍优先）。
for _d in ("_common/opt", "_common"):
    _p = Path(__file__).resolve().parent.parent / _d
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.append(str(_p))
import fc_common as fc
import stepconf

OUTDIR = "step2_kappa"
STEP = "step2_kappa"
FIT_DIR = "step1_fit"

SPEC = {
    # BTE q 网格：auto（按倒格矢长度自动匹配体系，phonopy 长度约定
    #   n_i = round(MESH_LENGTH*|b_i|)，对称性自洽）| 一个数（长度 L，Å）|
    #   "n n n"（显式网格，不做收敛）。2D 且 MESH_2D_VACUUM=on 时真空轴压 1。
    "MESH":             ("auto", "str"),
    "MESH_LENGTH":      (40.0, "float"),
    "MESH_2D_VACUUM":   ("on", "str"),
    # q 网格收敛（仅 auto/长度网格）：L 逐次 ×MESH_CONV_FACTOR 加密，直到
    #   MESH_CONV_T 下整个 κ 张量相邻变化 < MESH_CONV_TOL_PCT；上限
    #   MESH_CONV_MAX_LENGTH / MESH_CONV_MAX_POINTS（总 q 点数）。off = 只跑起始网格。
    "MESH_CONV":            ("auto", "str"),
    # per_component（默认）：xx/yy/zz 各用自己的 κ 归一，全部 < tol 才算收敛——
    #   各向异性体系慢收敛的面外 zz 不会被大的面内 κ 掩盖。分母有下限
    #   MESH_CONV_FLOOR×max κ_ii：小于最大分量 10% 的分量（2D 真空轴、极弱层间）
    #   按 0.1×max 归一，不会因为自身相对噪声卡住收敛。MESH_CONV_FLOOR=0 = 严格逐分量。
    #   max_ii：旧的整张量判据 max|Δκ_ij|/max|κ_ii|（zz 没收敛只告警，不影响判定）。
    "MESH_CONV_MODE":       ("per_component", "str"),
    "MESH_CONV_FLOOR":      (0.1, "float"),
    "MESH_CONV_TOL_PCT":    (3.0, "float"),
    "MESH_CONV_FACTOR":     (1.25, "float"),
    "MESH_CONV_MAX_LENGTH": (150.0, "float"),
    "MESH_CONV_MAX_POINTS": (64000, "int"),
    "MESH_CONV_T":          (300.0, "float"),
    # 温度扫描 (K)
    "T_MIN":            (100.0, "float"),
    "T_MAX":            (800.0, "float"),
    "T_STEP":           (100.0, "float"),
    "ISOTOPE":          (True, "bool"),
    # NAC：auto = 有 BORN 就开；on/off 强制
    "NAC":              ("auto", "str"),
    "BTE_METHOD":       ("rta", "str"),      # rta | lbte
    "P3PY_OMP_THREADS": (48, "int"),
    "SBATCH_QOS":       ("regular", "str"),
    "CONDA_SH":         ("", "str"),
    "CONDA_ENV":        ("", "str"),
    # 逐档截断 κ 扫描 + 三判据选截断：需要 step1_fit 用 CUT3_SCAN=auto 拟出
    #   cutoff_scan/cut3_<tag>/（每档 fc2/fc3）。auto = 有 ≥2 档就逐档跑 κ 再选；
    #   off = 只跑标称档。判据同 kl-dft-cpu S6：①残差 1-SE ②逐壳层稳定 ③κ 平台。
    "CUT3_KAPPA_SCAN":    ("auto", "str"),
    "CUT3_KAPPA_TOL_PCT": (5.0, "float"),
    "CUT3_CV_1SE_MULT":   (1.0, "float"),
    "CUT3_STABILITY_THR": (0.3, "float"),
    "CUT3_PICK":          ("smallest", "str"),
    # 2D 面内 κ 的层厚归一化（与 kl-dft-cpu S6 同一口径，真源 _common/thickness_2d.py）：
    #   phono3py 按整胞体积（含真空）归一，面内 κ 被 h⊥/d 稀释。kappa_summary.json 另写
    #   kappa_2d_normalized_* = 原始 κ × h⊥/d。vdw（默认：原子层跨度 + 上下 vdW 半径）|
    #   cell（不归一）| 数值（固定层厚 Å，如 MoS2 体相层间距 6.15）。
    "KAPPA_2D_THICKNESS": ("vdw", "str"),
}


def _make_disp_yaml(fit_dir, out):
    """Synthesise phono3py_disp.yaml from POSCAR + fc_dataset.json supercell."""
    from phonopy.interface.vasp import read_vasp
    from phono3py import Phono3py
    import numpy as np
    scm = None
    p = fit_dir / "fc_dataset.json"
    if p.is_file():
        try:
            scm = json.loads(p.read_text(encoding="utf-8")).get("supercell_matrix")
        except Exception:
            scm = None
    if not scm:
        sys.exit("[ERROR] %s 里没有 phono3py YAML，fc_dataset.json 也没有 "
                 "supercell_matrix —— 无法确定超胞" % fit_dir)
    uc = read_vasp(str(fit_dir / "POSCAR"))
    ph3 = Phono3py(uc, supercell_matrix=np.asarray(scm, float).tolist())
    ph3.save(str(out / "phono3py_disp.yaml"))
    print("[OK] 由 POSCAR + supercell_matrix 生成 phono3py_disp.yaml", flush=True)
    return "phono3py_disp.yaml"


def _find_disp_yaml(fit_dir, out):
    cands = [fit_dir / "phono3py" / "phono3py_disp.yaml",
             fit_dir / "phono3py" / "phono3py_params.yaml",
             fit_dir / "phono3py_disp.yaml",
             fit_dir / "phono3py_params.yaml"]
    for c in cands:
        if c.is_file():
            shutil.copyfile(str(c), str(out / c.name))
            print("[..] phono3py YAML <- %s" % c, flush=True)
            return c.name
    return _make_disp_yaml(fit_dir, out)


def _mesh_cfg(conf, fit_dir):
    """Validate MESH and read is_2d; the mesh itself is resolved on the compute
    node on phono3py's primitive cell (kappa_driver._apply_mesh), so an explicit
    "n n n" and the 2D vacuum axis are both taken in that basis."""
    m = str(conf["MESH"]).strip()
    parts = m.split()
    if m.lower() == "auto":
        pass
    elif len(parts) == 3:
        try:
            [int(x) for x in parts]
        except ValueError:
            sys.exit("[ERROR] MESH 需要 auto、一个长度或三个整数，得到 %r" % m)
    elif len(parts) == 1:
        try:
            if float(parts[0]) <= 0:
                raise ValueError
        except ValueError:
            sys.exit("[ERROR] MESH 需要 auto、一个长度或三个整数，得到 %r" % m)
    else:
        sys.exit("[ERROR] MESH 需要 auto、一个长度或三个整数，得到 %r" % m)
    if float(conf["MESH_LENGTH"]) <= 0:
        sys.exit("[ERROR] MESH_LENGTH 必须 > 0")
    if float(conf["MESH_CONV_FACTOR"]) <= 1.0:
        sys.exit("[ERROR] MESH_CONV_FACTOR 必须 > 1")
    is_2d = False
    ps = fit_dir / "phonon_summary.json"
    if ps.is_file():
        try:
            is_2d = bool(json.loads(ps.read_text(encoding="utf-8")).get("is_2d"))
        except Exception:
            is_2d = False
    vac = str(conf["MESH_2D_VACUUM"]).lower() in ("on", "true", "1", "yes", "auto")
    conv = str(conf["MESH_CONV"]).strip().lower() not in (
        "off", "false", "0", "no", "none", "")
    if len(parts) == 3 and conv:
        print("[..] MESH 是显式网格 %r —— 不做网格收敛（要收敛改 MESH=auto）" % m,
              flush=True)
    return {
        "mesh": m.lower() if m.lower() == "auto" else m,
        "mesh_length": float(conf["MESH_LENGTH"]),
        "is_2d": is_2d,
        "mesh_2d_vacuum": vac,
        "mesh_conv": "auto" if conv else "off",
        "mesh_conv_mode": str(conf["MESH_CONV_MODE"] or "per_component").strip().lower(),
        "mesh_conv_floor": (0.1 if conf["MESH_CONV_FLOOR"] is None
                            else float(conf["MESH_CONV_FLOOR"])),
        "mesh_conv_tol_pct": float(conf["MESH_CONV_TOL_PCT"]),
        "mesh_conv_factor": float(conf["MESH_CONV_FACTOR"]),
        "mesh_conv_max_length": float(conf["MESH_CONV_MAX_LENGTH"]),
        "mesh_conv_max_points": int(conf["MESH_CONV_MAX_POINTS"]),
        "mesh_conv_t": float(conf["MESH_CONV_T"]),
    }


def _fit_identity(fit):
    """{'method_label': 'ALASSO', 'nominal_cut3_A': 5.86} from the fit."""
    try:
        fc = json.loads((Path(fit) / "fit_config.json").read_text(encoding="utf-8"))
    except Exception:
        return {"method_label": None, "nominal_cut3_A": None}
    eng = str(fc.get("engine") or "")
    meth = {"phono3py": fc.get("fc_calc"), "pheasy": fc.get("pheasy_method"),
            "hiphive": fc.get("hiphive_fit_method")}.get(eng)
    raw = {"phono3py": fc.get("fc3_cutoff"), "pheasy": fc.get("pheasy_c3_cutoff"),
           "hiphive": fc.get("hiphive_cutoff3")}.get(eng)
    try:
        cut = float(raw) if raw not in (None, "") else None
    except (TypeError, ValueError):
        cut = None
    # 图上只标方法本身（不带引擎前缀），例如 ALASSO / symfc
    _label = str(meth).strip() if meth not in (None, "") else None
    return {"method_label": _label, "nominal_cut3_A": cut}


def two_d_norm(fit_dir, mode):
    """2D thickness normalisation for kappa_driver (None for a 3D cell or when
    the geometry cannot be read).  Same source as kl-dft-cpu: factor = h_perp/d."""
    pos = Path(fit_dir) / "POSCAR"
    if str(mode or "vdw").strip().lower() in ("cell", "lz", "none", "off"):
        return None
    try:
        import dim_common
        from thickness_2d import slab_geometry_from_poscar
        dim, axis, _vac = dim_common.detect_dimension(str(pos))
        if dim != "2d":
            return None
        g = slab_geometry_from_poscar(str(pos), vac_axis=axis, mode=str(mode or "vdw"))
    except Exception as e:                               # noqa: BLE001
        print("[WARN] 2D 层厚归一化算不出来（只写原始 κ，面内 κ 被真空稀释、不可直接"
              "使用）：%s" % e, flush=True)
        return None
    g = dict(g)
    g["vac_axis"] = int(axis)
    g["note"] = ("kappa_2d_normalized = kappa_raw * h_perp/d (h_perp = V/A); only the "
                 "in-plane components are physical for a 2D layer")
    print("[..] 2D κ 归一化：h⊥=%.3f Å 层厚 d=%.3f Å（%s）→ factor=%.4f"
          % (g["h_perp_A"], g["thickness_d_A"] or 0.0, g["thickness_convention"],
             g["kappa_2d_norm_factor"]), flush=True)
    return g


METHODS_TAIL = """
# ---- FIT_METHODS: kappa for the other fitting methods, then the comparison ----
set +e
for d in %s; do
    ( cd "methods/$d" && python kappa_driver.py kappa_config.json ) \\
        > "methods/$d/run.log" 2>&1 || echo "[WARN] kappa for $d failed -- see methods/$d/run.log"
done
python kappa_driver.py compare
set -e
"""


def main():
    """S2_kappa gen: the primary fit (step1_fit/) -> step2_kappa/; every other
    FIT_METHODS fit (step1_fit/methods/<tag>/, listed in step1_fit/methods.json)
    -> step2_kappa/methods/<tag>/, run in the same job, then
    `kappa_driver.py compare` writes methods_compare.{json,png,pdf}."""
    cwd = Path.cwd()
    out = cwd / OUTDIR
    conf = stepconf.load(SPEC, STEP, strict="warn")
    fit = cwd / FIT_DIR
    _gen_one(conf, fit, out, cwd, write_submit=True)
    mj = fit / "methods.json"
    entries = []
    if mj.is_file():
        try:
            entries = json.loads(mj.read_text(encoding="utf-8")).get("methods") or []
        except Exception:
            entries = []
    done = []
    for e in entries:
        if e.get("primary"):
            continue
        src = fit / e["dir"]
        if not ((src / "fc2.hdf5").is_file() and (src / "fc3.hdf5").is_file()):
            print("[WARN] 方法 %s 没有 fc2/fc3（S1 里失败了？见 %s/run.log）—— 跳过"
                  % (e["tag"], src), flush=True)
            e["skipped"] = "no fc2/fc3 from S1"
            continue
        print("\n[..] ==== 方法 %s ====" % e["tag"], flush=True)
        _gen_one(conf, src, out / e["dir"], cwd, write_submit=False)
        done.append(e["dir"].split("/", 1)[1])
    if entries:
        (out / "methods.json").write_text(
            json.dumps({"methods": entries}, indent=2, ensure_ascii=False),
            encoding="utf-8", newline="\n")
    if done:
        with open(out / "submit.sh", "a", encoding="utf-8", newline="\n") as fh:
            fh.write(METHODS_TAIL % " ".join(done))
        print("[..] 另 %d 个方法的 κ 在同一作业内串行算，最后出 methods_compare"
              % len(done), flush=True)


def _gen_one(conf, fit, out, cwd, write_submit=True):
    """One S2 recipe: kappa from the fc2/fc3 in `fit` into `out`."""
    out.mkdir(parents=True, exist_ok=True)
    if not (fit / "fc2.hdf5").is_file() or not (fit / "fc3.hdf5").is_file():
        sys.exit("[ERROR] %s 缺 fc2.hdf5/fc3.hdf5 —— step1_fit 还没拟合完" % fit)
    enable_fc = 3
    for f in ("fc2.hdf5", "fc3.hdf5"):
        shutil.copyfile(str(fit / f), str(out / f))

    nac_cfg = str(conf["NAC"]).strip().lower()
    born = fit / "BORN"
    if nac_cfg in ("auto", ""):
        nac = born.is_file()
    else:
        nac = nac_cfg in ("on", "true", "1", "yes")
    if nac and born.is_file():
        shutil.copyfile(str(born), str(out / "BORN"))
    elif nac and not born.is_file():
        print("[WARN] NAC=on 但 %s 没有 BORN，改为 nonac" % fit, flush=True)
        nac = False
    print("[..] NAC=%s（BORN %s）" % (nac, "有" if born.is_file() else "无"), flush=True)

    disp_name = _find_disp_yaml(fit, out)
    mcfg = _mesh_cfg(conf, fit)

    method = str(conf["BTE_METHOD"] or "rta").lower()
    if method not in ("rta", "lbte"):
        sys.exit("[ERROR] BTE_METHOD must be rta or lbte")
    t0, t1, dt = float(conf["T_MIN"]), float(conf["T_MAX"]), float(conf["T_STEP"])
    if dt <= 0 or t1 < t0:
        sys.exit("[ERROR] T_MIN/T_MAX/T_STEP 非法：%s %s %s" % (t0, t1, dt))
    temps = []
    t = t0
    while t <= t1 + 1e-9:
        temps.append(round(t, 6))
        t += dt

    # ---- 逐档截断：把 step1_fit/cutoff_scan/cut3_<tag>/ 的力常数搬进本步 ----
    #   只有 S1_fit 用 CUT3_SCAN=auto 时才有这些目录；够 2 档才做 κ 扫描。
    _want = str(conf["CUT3_KAPPA_SCAN"]).strip().lower() not in (
        "off", "false", "0", "no", "none", "")
    cut_scans = []
    scan_root = fit / "cutoff_scan"
    if _want and scan_root.is_dir():
        for cd in sorted(scan_root.glob("cut3_*")):
            if (cd / "fc2.hdf5").is_file() and (cd / "fc3.hdf5").is_file():
                tag = cd.name[len("cut3_"):]
                dst = out / cd.name
                dst.mkdir(exist_ok=True)
                shutil.copyfile(str(cd / "fc2.hdf5"), str(dst / "fc2.hdf5"))
                shutil.copyfile(str(cd / "fc3.hdf5"), str(dst / "fc3.hdf5"))
                cut_scans.append(tag)
    scan = len(cut_scans) >= 2
    if scan:
        shutil.copyfile(str(fit / "cutoff_scan.json"), str(out / "cutoff_scan.json"))
        print("[..] 逐档 κ 扫描：%d 档 %s" % (len(cut_scans), cut_scans), flush=True)
    else:
        cut_scans = []
        if _want:
            print("[WARN] step1_fit 没有逐档截断（CUT3_SCAN=off 或候选不足 2 档）—— 只跑标称截断，"
                  "画不出 κ 随截断的变化图。要扫描：conf --set params.CUT3_SCAN=auto 后 retry S1_fit",
                  flush=True)
    cfg = {
        "disp_yaml": disp_name,
        "scan": scan,
        "cut_scans": cut_scans,
        "cut3_kappa_tol_pct": float(conf["CUT3_KAPPA_TOL_PCT"]),
        "cut3_cv_1se_mult": float(conf["CUT3_CV_1SE_MULT"]),
        "cut3_stability_thr": float(conf["CUT3_STABILITY_THR"]),
        "cut3_pick": str(conf["CUT3_PICK"]),
        **mcfg,
        "temperatures": temps,
        "isotope": bool(conf["ISOTOPE"]),
        "nac": nac,
        "bte_method": method,
        "enable_fc": enable_fc,
        "kappa_2d_norm": (two_d_norm(fit, conf["KAPPA_2D_THICKNESS"])
                          if mcfg["is_2d"] else None),
        "source_fc": str(fit),
        # the fit's own cutoff and method: the per-method kappa figure is drawn
        # even without a cutoff scan (one point at the nominal cutoff)
        **_fit_identity(fit),
    }
    (out / "kappa_config.json").write_text(
        json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8",
        newline="\n")

    here = Path(__file__).resolve().parent
    for _f in ("kappa_driver.py", "cut3_select.py"):
        if not (here / _f).is_file():
            sys.exit("[ERROR] %s missing -- is it listed in gen_need?" % _f)
        shutil.copyfile(str(here / _f), str(out / _f))
    # the fit's own summary (method, gate, force RMSE) travels with the kappa so
    # the method comparison and the RMSE panel need no path back into step1_fit
    if (fit / "fc_fit_summary.json").is_file():
        shutil.copyfile(str(fit / "fc_fit_summary.json"), str(out / "fc_fit_summary.json"))
    if not write_submit:
        print("[DONE] %s: kappa_config.json ready" % out.relative_to(cwd), flush=True)
        return

    tpl = fc.resolve_submit(here, "submit_kappa")
    subs = {
        "JOBNAME": fc.new_jobname(cwd, "S2kappa"),
        "CONDA_SH": str(conf["CONDA_SH"] or ""),
        "CONDA_ENV": str(conf["CONDA_ENV"] or ""),
        "CPUS_PER_TASK": str(int(conf["P3PY_OMP_THREADS"] or 48)),
        "QOS": str(conf["SBATCH_QOS"] or "regular"),
    }
    fc.write_submit(tpl, out / "submit.sh", subs)
    stepconf.apply_submit(out / "submit.sh", dict(conf.submit or {}))

    print("[..] BTE=%s mesh=%s%s conv=%s T=%g..%g/%g isotope=%s 2D=%s"
          % (method, mcfg["mesh"],
             " (L=%g A)" % mcfg["mesh_length"] if mcfg["mesh"] == "auto" else "",
             mcfg["mesh_conv"], t0, t1, dt, bool(conf["ISOTOPE"]), mcfg["is_2d"]))
    print("[DONE] %s: kappa_config.json + submit.sh + kappa_driver.py ready; "
          "the compute node writes kappa_summary.json (BTE done marker)."
          % OUTDIR)


if __name__ == "__main__":
    main()
