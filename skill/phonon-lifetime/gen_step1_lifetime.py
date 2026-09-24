#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_lifetime.py -- 声子寿命 / 线宽分析（登录节点，run: gen）。

数据来源（优先级）：
  1) SOURCE=kappa_hdf5 或 KAPPA_HDF5 指定：直接读 phono3py 的 kappa-m*.hdf5/gamma-m*.hdf5；
  2) SOURCE=auto：在材料目录树里找上游 kl-dft-cpu（或指定技能）的 kappa hdf5；
  3) 找不到 hdf5 但找到 fc2.hdf5+fc3.hdf5+phono3py_disp.yaml：
       RUN_PHONO3PY=true 时现跑 phono3py --br 生成 gamma；否则给出明确指引。

核心物理：tau[ps] = 1/(2*2*pi*gamma_total[THz])，
          gamma_total = gamma_anh + gamma_isotope（phono3py 官方约定，见 lifetime_common）。
产出：lifetime_summary.json（marker LIFETIME_DONE）、lifetime_data.csv、lifetime_T.csv、
      lifetime_vs_frequency.png / lifetime_vs_period.png / ioffe_regel.png /
      lifetime_vs_T.png（可选 lifetime_vs_q.png）。
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

OUTDIR = "step1_lifetime"
STEP = "step1_lifetime"


def _ensure_python():
    """当前 python 缺 numpy/h5py 时，re-exec 到带这些包的 python（一次性）。

    优先用 step.conf 的 CONDA_SH/CONDA_ENV；否则扫常见的 conda envs。
    标记 LIFETIME_REEXEC 防止递归。
    """
    try:
        import numpy  # noqa: F401
        import h5py   # noqa: F401
        return
    except Exception:
        pass
    if os.environ.get("LIFETIME_REEXEC"):
        sys.exit("[ERROR] 已 re-exec 过仍缺 numpy/h5py；请检查 CONDA_SH/CONDA_ENV。")
    cands = []
    env_py = os.environ.get("LIFETIME_PYTHON")
    if env_py:
        cands.append(env_py)
    # 从 step.conf 读 CONDA_ENV（若可解析）
    try:
        import stepconf as _sc  # noqa
        _cf = _sc.load({"CONDA_SH": ("", "str"), "CONDA_ENV": ("", "str")},
                       None, strict=False)
        sh = str(_cf["CONDA_SH"] or "")
        env = str(_cf["CONDA_ENV"] or "")
        if env:
            p = Path(os.path.expanduser(env))
            if (p / "bin" / "python").is_file():
                cands.append(str(p / "bin" / "python"))
            if sh and os.path.isfile(os.path.expanduser(sh)):
                # source conda.sh 后取 env 的 python
                try:
                    r = subprocess.run(
                        ["bash", "-lc",
                         "source %s && conda activate %s && command -v python"
                         % (os.path.expanduser(sh), env)],
                        capture_output=True, encoding="utf-8", timeout=120)
                    if r.returncode == 0 and r.stdout.strip():
                        cands.append(r.stdout.strip().splitlines()[-1])
                except Exception:
                    pass
    except Exception:
        pass
    import glob
    home = os.path.expanduser("~")
    for pat in ("miniconda3/envs/*/bin/python", "anaconda3/envs/*/bin/python",
                "miniforge3/envs/*/bin/python", "mambaforge/envs/*/bin/python"):
        cands.extend(sorted(glob.glob(os.path.join(home, pat))))
    cands.append(sys.executable)
    seen = set()
    for p in cands:
        if not p or p in seen or not os.path.exists(p):
            continue
        seen.add(p)
        if os.path.realpath(p) == os.path.realpath(sys.executable):
            continue
        try:
            ok = subprocess.run([p, "-c", "import numpy, h5py"],
                                capture_output=True, timeout=180).returncode == 0
        except Exception:
            continue
        if ok:
            print("[..] 当前解释器缺 numpy/h5py，改用 %s 重跑" % p, flush=True)
            os.execve(p, [p, os.path.abspath(__file__)] + sys.argv[1:],
                      dict(os.environ, LIFETIME_REEXEC="1"))


_ensure_python()

import numpy as np  # noqa: E402

import stepconf  # noqa: E402
import lifetime_common as LC  # noqa: E402
import lifetime_plots as LP  # noqa: E402

SPEC = {
    # ---- 数据来源 ----
    "SOURCE":          ("auto", "str"),        # auto | kappa_hdf5 | fc
    "KAPPA_HDF5":      ("", "str"),            # 显式指定 hdf5（优先）
    "UPSTREAM_SKILL":  ("kl-dft-cpu", "str"),
    "UPSTREAM_STEP":   ("step6_kappa", "str"),
    # ---- 分析 ----
    "TEMPERATURES":    ("", "str"),            # 逗号/空格分隔；空=全部
    "FOCUS_T":         ("300", "float"),
    "IMAG_THR":        ("0.10", "float"),      # THz，低于此频率的模排除
    "ACOUSTIC_FMAX":   ("", "str"),            # THz；空=自动
    "OMEGA_TAU_THR":   ("1.0", "float"),       # 过阻尼判据阈值
    "FMAX_PLOT":       ("20.0", "float"),
    "MAX_PLOT_T":      ("8", "int"),
    "Q_DIRECTION":     ("", "str"),            # 如 "0 0 1"；空=不出沿 q 的图
    "Q_DIRECTION_LABEL": ("", "str"),
    "MAKE_PLOT":       ("true", "bool"),
    "INCLUDE_ISOTOPE": ("true", "bool"),       # 寿命含同位素散射（与 phono3py kappa 一致）
    # ---- 可选：现跑 phono3py ----
    "RUN_PHONO3PY":    ("false", "bool"),
    "PHONO3PY_MESH":   ("15 15 15", "str"),
    "PHONO3PY_TS":     ("300", "str"),
    "PHONO3PY_CMD":    ("phono3py-load", "str"),
    "PHONO3PY_TIMEOUT": ("7200", "int"),
    # ---- 环境（由 setting/<hpc>.yaml 注入或项目 step.conf 覆盖）----
    "CONDA_SH":        ("", "str"),
    "CONDA_ENV":       ("", "str"),
}


def _parse_nums(s):
    if s is None:
        return []
    return [float(x) for x in re.split(r"[,\s]+", str(s).strip()) if x]


def _parse_direction(s):
    v = re.split(r"[,\s]+", str(s or "").strip())
    v = [x for x in v if x]
    return [float(x) for x in v] if v else None


def _has_gamma(p):
    return LC.hdf5_has_gamma(p)


def _resolve_dataset(conf, cwd, outdir):
    """-> (kind, payload)  kind in {kappa_hdf5, fc}。"""
    explicit = conf["KAPPA_HDF5"]
    if explicit:
        p = Path(os.path.expanduser(str(explicit)))
        if not p.is_absolute():
            p = (cwd / p).resolve()
        if not p.is_file():
            sys.exit("[ERROR] KAPPA_HDF5=%s 不存在" % explicit)
        if not _has_gamma(p):
            sys.exit("[ERROR] %s 里没有 gamma/frequency 数据集，不是 phono3py kappa/gamma hdf5" % p)
        return "kappa_hdf5", p

    # 本技能目录 / 步骤目录里已有的 hdf5
    for d in (outdir, cwd):
        for pat in ("kappa-m*.hdf5", "gamma-m*.hdf5", "kappa*.hdf5", "*.hdf5"):
            for h in sorted(Path(d).glob(pat)):
                if _has_gamma(h):
                    return "kappa_hdf5", h

    src = str(conf["SOURCE"] or "auto").lower()
    if src in ("auto", "upstream", "kappa_hdf5"):
        up = LC.find_kappa_hdf5(cwd, conf["UPSTREAM_SKILL"], conf["UPSTREAM_STEP"])
        if up:
            return "kappa_hdf5", up
    fc = LC.find_fc_dataset(cwd, conf["UPSTREAM_SKILL"], conf["UPSTREAM_STEP"])
    if fc:
        return "fc", fc
    return None, None


def _run_phono3py(conf, cwd, outdir, fc):
    """从 fc 数据集现跑 phono3py --br，返回生成的 hdf5 路径。"""
    base, fc2, fc3, disp = fc
    for f in (fc2, fc3, disp):
        shutil.copyfile(str(f), str(outdir / f.name))
    mesh = str(conf["PHONO3PY_MESH"] or "15 15 15")
    ts = str(conf["PHONO3PY_TS"] or "300")
    exe = str(conf["PHONO3PY_CMD"] or "phono3py-load")
    cmd = ("%s phono3py_disp.yaml --mesh %s --br --ts=\"%s\" "
           "2>&1 | tee phono3py_lifetime.log" % (exe, mesh, ts))
    py = sys.executable
    bindir = os.path.dirname(py)
    env_prefix = ("export PATH=%s:$PATH; " % bindir) if bindir else ""
    full = "cd %s && %s( %s )" % (str(outdir), env_prefix, cmd)
    print("[..] 现跑 phono3py：%s" % cmd, flush=True)
    try:
        rc = subprocess.run(["bash", "-lc", full],
                            timeout=int(conf["PHONO3PY_TIMEOUT"])).returncode
    except subprocess.TimeoutExpired:
        sys.exit("[ERROR] phono3py 超过 PHONO3PY_TIMEOUT=%ss 未完成" % conf["PHONO3PY_TIMEOUT"])
    if rc != 0:
        sys.exit("[ERROR] phono3py 退出码 %d，见 %s/phono3py_lifetime.log" % (rc, outdir))
    for pat in ("kappa-m*.hdf5", "gamma-m*.hdf5", "*.hdf5"):
        for h in sorted(Path(outdir).glob(pat)):
            if _has_gamma(h):
                return h
    sys.exit("[ERROR] phono3py 跑完但没找到含 gamma 的 hdf5")


def _iter_grid_kappa(data, focus_T):
    """若 hdf5 里有 kappa，返回最接近 focus_T 的对角张量。"""
    kappa = data.get("kappa")
    temps = data.get("temperatures")
    if kappa is None or temps is None or len(temps) == 0:
        return None
    i = int(np.argmin(np.abs(np.asarray(temps) - float(focus_T))))
    row = np.asarray(kappa)[i]
    return {"T_K": float(temps[i]), "kappa_xx_yy_zz": [float(x) for x in row[:6]]}


def main():
    cwd = Path.cwd()
    outdir = cwd / OUTDIR
    outdir.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)

    kind, payload = _resolve_dataset(conf, cwd, outdir)
    if kind is None:
        sys.exit(
            "[ERROR] 找不到声子寿命数据源。请任选其一：\n"
            "  a) 先跑上游 kl-dft-cpu（或 kl-mlff-*）到出 kappa-m*.hdf5；\n"
            "  b) 在 step.conf 设 params.KAPPA_HDF5=<kappa/gamma hdf5 路径>；\n"
            "  c) 若已有 fc2.hdf5/fc3.hdf5/phono3py_disp.yaml，设 params.RUN_PHONO3PY=true。")

    if kind == "fc":
        if not conf["RUN_PHONO3PY"]:
            base = payload[0]
            sys.exit("[ERROR] 只找到力常数数据集 %s，但没有可直接读的 gamma hdf5。"
                     "设 params.RUN_PHONO3PY=true 让本步现跑 phono3py --br。" % base)
        h5 = _run_phono3py(conf, cwd, outdir, payload)
        src_kind = "phono3py_run"
    else:
        h5 = Path(payload)
        src_kind = "upstream_hdf5" if str(h5).startswith(str(cwd.parent)) else "kappa_hdf5"

    data = LC.read_kappa_hdf5(h5)
    temperatures = data["temperatures"]
    frequencies = data["frequency"]
    gamma = data["gamma"]
    gamma_iso = data.get("gamma_isotope")
    if not conf["INCLUDE_ISOTOPE"]:
        print("[..] INCLUDE_ISOTOPE=false：只用三声子 gamma，寿命会偏长", flush=True)
        gamma_iso = None
    # 总和线宽：phono3py 的 kappa/mode_kappa 用的就是 gamma + gamma_isotope，
    # 只取 gamma 会系统性高估寿命（Si 300K 中位约 +19%）。
    gsum = LC.effective_gamma(gamma, gamma_iso)
    qpoints = data.get("qpoints")
    weights = data.get("weights")
    if qpoints is None or weights is None:
        sys.exit("[ERROR] %s 缺 qpoint/weight 数据集" % h5)

    tmap = LC.pick_temperature_indices(temperatures, _parse_nums(conf["TEMPERATURES"]))
    t_indices = sorted(tmap.values(), key=lambda i: temperatures[i])
    focus = float(conf["FOCUS_T"])
    focus_idx = int(np.argmin(np.abs(np.asarray(temperatures) - focus)))

    ac_fmax = conf["ACOUSTIC_FMAX"]
    ac_fmax = float(ac_fmax) if str(ac_fmax).strip() else None
    imag_thr = float(conf["IMAG_THR"])
    omega_thr = float(conf["OMEGA_TAU_THR"])

    per_T = [LC.summarize_temperature(frequencies, gamma, temperatures, ti,
                                      acoustic_fmax=ac_fmax, imag_thr=imag_thr,
                                      omega_tau_threshold=omega_thr,
                                      gamma_isotope=gamma_iso)
             for ti in t_indices]

    # ---- 逐模 CSV ----
    csv_mode = outdir / "lifetime_data.csv"
    LC.write_mode_csv(str(csv_mode), temperatures, frequencies, gamma, qpoints, weights,
                      group_velocity=data.get("group_velocity"), acoustic_fmax=ac_fmax,
                      imag_thr=imag_thr, t_indices=t_indices, gamma_isotope=gamma_iso)

    # ---- 每温度汇总 CSV ----
    csv_T = outdir / "lifetime_T.csv"
    with open(csv_T, "w", encoding="utf-8", newline="") as fh:
        import csv as _csv
        w = _csv.writer(fh)
        w.writerow(["T_K", "n_all", "tau_median_all_ps", "tau_median_all_anharm_ps",
                    "tau_median_acoustic_ps", "tau_median_optical_ps",
                    "omega_tau_median_acoustic", "overdamped_frac_acoustic",
                    "acoustic_fmax_THz"])
        for d in per_T:
            a = d["acoustic"]
            w.writerow(["%.4f" % d["T_K"], d["all"].get("n_modes"),
                        "%.6g" % d["all"].get("tau_median_ps", float("nan")),
                        "%.6g" % d["all"].get("tau_median_anharmonic_ps", float("nan")),
                        "%.6g" % a.get("tau_median_ps", float("nan")),
                        "%.6g" % d["optical"].get("tau_median_ps", float("nan")),
                        "%.6g" % a.get("omega_tau_median", float("nan")),
                        "%.6g" % a.get("overdamped_frac", float("nan")),
                        "%.4f" % d["acoustic_fmax_THz"]])

    # ---- 图 ----
    plots = []
    if conf["MAKE_PLOT"] and LP.HAVE_MPL:
        plot_T = t_indices[:int(conf["MAX_PLOT_T"])]
        fmax_plot = float(conf["FMAX_PLOT"])
        for fn in (
            LP.plot_lifetime_vs_frequency(str(outdir), temperatures, frequencies, gsum,
                                          plot_T, fmax_plot=fmax_plot, imag_thr=imag_thr),
            LP.plot_lifetime_vs_period(str(outdir), temperatures, frequencies, gsum,
                                       plot_T, fmax_plot=fmax_plot, imag_thr=imag_thr),
            LP.plot_ioffe_regel(str(outdir), temperatures, frequencies, gsum,
                                plot_T, fmax_plot=fmax_plot, imag_thr=imag_thr),
            LP.plot_lifetime_T(str(outdir), per_T),
        ):
            if fn:
                plots.append(os.path.basename(fn))
        qdir = _parse_direction(conf["Q_DIRECTION"])
        if qdir and qpoints is not None:
            fn = LP.plot_lifetime_vs_q(str(outdir), temperatures, frequencies, gsum,
                                       qpoints, qdir, [focus_idx],
                                       direction_label=str(conf["Q_DIRECTION_LABEL"] or ""),
                                       imag_thr=imag_thr)
            if fn:
                plots.append(os.path.basename(fn))
    elif conf["MAKE_PLOT"] and not LP.HAVE_MPL:
        print("[warn] 无 matplotlib，跳过出图（JSON/CSV 照常产出）")

    summary = {
        "LIFETIME_DONE": True,
        "step": STEP,
        "source": {
            "kind": src_kind,
            "file": str(h5),
            "mesh": [int(x) for x in np.asarray(data.get("mesh", [])).ravel().tolist()],
            "n_grid_points": int(frequencies.shape[0]),
            "n_bands": int(frequencies.shape[1]),
        },
        "tau_convention": "tau_ps = 1/(2*2*pi*gamma_total_THz), gamma_total = "
                          "gamma_anh + gamma_isotope  [phono3py kaccum/kdeplot + "
                          "conductivity/base.py:_get_main_diagonal]",
        "tau_factor": LC.TAU_FACTOR,
        "isotope_included": bool(gamma_iso is not None and np.any(np.asarray(gamma_iso) > 0)),
        "boundary_mfp_ang": (float(np.asarray(data["boundary_mfp"]).ravel()[0])
                             if data.get("boundary_mfp") is not None else None),
        "temperatures_all_K": [float(x) for x in temperatures],
        "temperatures_analyzed_K": [float(temperatures[i]) for i in t_indices],
        "focus_T_K": float(temperatures[focus_idx]),
        "imag_thr_THz": imag_thr,
        "omega_tau_threshold": omega_thr,
        "per_temperature": per_T,
        "kappa_from_hdf5": _iter_grid_kappa(data, focus),
        "plots": plots,
        "files": {
            "per_mode_csv": os.path.basename(csv_mode),
            "per_temperature_csv": os.path.basename(csv_T),
            "summary": "lifetime_summary.json",
        },
    }
    LC.dump_json(str(outdir / "lifetime_summary.json"), summary)

    foc = (min(per_T, key=lambda d: abs(d["T_K"] - float(temperatures[focus_idx])))
           if per_T else None)
    f = foc["acoustic"] if foc else {}
    print("[DONE] %s：%d 个温度，%d 个模；%.0fK 声学支 tau_median=%s ps，omega*tau_median=%s"
          % (OUTDIR, len(per_T), frequencies.size, foc["T_K"] if foc else float("nan"),
             ("%.3g" % f.get("tau_median_ps")) if f.get("tau_median_ps") else "n/a",
             ("%.3g" % f.get("omega_tau_median")) if f.get("omega_tau_median") else "n/a"))


if __name__ == "__main__":
    main()
