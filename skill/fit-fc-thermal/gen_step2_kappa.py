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
import fc_common as fc
import stepconf

OUTDIR = "step2_kappa"
STEP = "step2_kappa"
FIT_DIR = "step1_fit"

SPEC = {
    # BTE q 网格（"24 24 24"）。2D 且 MESH_2D_VACUUM=on 时把真空轴压 1。
    "MESH":             ("24 24 24", "str"),
    "MESH_2D_VACUUM":   ("on", "str"),
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


def _mesh_str(conf, fit_dir):
    import numpy as np
    mesh = [int(x) for x in str(conf["MESH"]).split()]
    if len(mesh) != 3:
        sys.exit("[ERROR] MESH 需要三个整数，得到 %r" % conf["MESH"])
    is_2d = False
    ps = fit_dir / "phonon_summary.json"
    if ps.is_file():
        try:
            is_2d = bool(json.loads(ps.read_text(encoding="utf-8")).get("is_2d"))
        except Exception:
            is_2d = False
    on = str(conf["MESH_2D_VACUUM"]).lower() in ("on", "true", "1", "yes", "auto")
    if is_2d and on:
        try:
            lat, _ = fc.read_poscar_cell_frac(fit_dir / "POSCAR")
            vac = int(np.argmax([np.linalg.norm(v) for v in np.asarray(lat, float)]))
            mesh[vac] = 1
            print("[..] 2D：真空轴 %d 的网格压 1 → %s" % (vac, mesh), flush=True)
        except Exception as e:  # noqa: BLE001
            print("[WARN] 2D 真空轴判定失败，MESH 原样：%s" % e, flush=True)
    return " ".join(str(x) for x in mesh)


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    fit = cwd / FIT_DIR

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
    mesh = _mesh_str(conf, fit)

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
            print("[..] step1_fit 只有标称一档（未开 CUT3_SCAN=auto）—— 只跑标称截断",
                  flush=True)
    cfg = {
        "disp_yaml": disp_name,
        "scan": scan,
        "cut_scans": cut_scans,
        "cut3_kappa_tol_pct": float(conf["CUT3_KAPPA_TOL_PCT"]),
        "cut3_cv_1se_mult": float(conf["CUT3_CV_1SE_MULT"]),
        "cut3_stability_thr": float(conf["CUT3_STABILITY_THR"]),
        "cut3_pick": str(conf["CUT3_PICK"]),
        "mesh": mesh,
        "temperatures": temps,
        "isotope": bool(conf["ISOTOPE"]),
        "nac": nac,
        "bte_method": method,
        "enable_fc": enable_fc,
        "source_fc": str(fit),
    }
    (out / "kappa_config.json").write_text(
        json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8",
        newline="\n")

    here = Path(__file__).resolve().parent
    for _f in ("kappa_driver.py", "cut3_select.py"):
        if not (here / _f).is_file():
            sys.exit("[ERROR] %s missing -- is it listed in gen_need?" % _f)
        shutil.copyfile(str(here / _f), str(out / _f))

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

    print("[..] BTE=%s mesh=%s T=%g..%g/%g isotope=%s"
          % (method, mesh, t0, t1, dt, bool(conf["ISOTOPE"])))
    print("[DONE] %s: kappa_config.json + submit.sh + kappa_driver.py ready; "
          "the compute node writes kappa_summary.json (BTE done marker)."
          % OUTDIR)


if __name__ == "__main__":
    main()
