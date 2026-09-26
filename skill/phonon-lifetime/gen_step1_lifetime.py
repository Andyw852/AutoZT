#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_lifetime.py -- 声子寿命 / 线宽分析（登录节点，run: gen）。

数据来源（优先级）：
  1) SOURCE=kappa_hdf5 或 KAPPA_HDF5 指定：直接读 phono3py 的 kappa-m*.hdf5/gamma-m*.hdf5；
  2) SOURCE=auto：在材料目录树里找上游 kl-dft-cpu（或指定技能）的 kappa hdf5；
  3) 找不到 hdf5 但找到 fc2.hdf5+fc3.hdf5+phono3py_disp.yaml：
       RUN_PHONO3PY=true 时现跑 phono3py --br --isotope 生成 gamma；否则给出明确指引。

核心物理：tau[ps] = 1/(2*2*pi*gamma_total[THz])，
          gamma_total = gamma_anh + gamma_isotope + gamma_bd（见 lifetime_common）。
统计量按不可约网格 weight 加权；声学支默认取最低 3 条支、不加频率窗。
产出：lifetime_summary.json（marker LIFETIME_DONE）、lifetime_data.csv、lifetime_T.csv、
      lifetime_vs_frequency.png / lifetime_vs_period.png / ioffe_regel.png /
      lifetime_vs_T.png（可选 lifetime_vs_q.png）。
"""
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
    "ACOUSTIC_FMAX":   ("", "str"),            # THz；空=只用最低 3 条支，不加频率窗
    "OMEGA_TAU_THR":   ("1.0", "float"),       # 过阻尼判据阈值
    "FMAX_PLOT":       ("20.0", "float"),
    "MAX_PLOT_T":      ("8", "int"),
    "Q_DIRECTION":     ("", "str"),            # 如 "0 0 1"；空=不出沿 q 的图
    "Q_DIRECTION_LABEL": ("", "str"),
    "MAKE_PLOT":       ("true", "bool"),
    "INCLUDE_ISOTOPE": ("true", "bool"),       # 寿命含同位素散射（与 phono3py kappa 一致）
    "CSV_ALL_TEMPERATURES": ("false", "bool"), # 逐模 CSV 是否写全部温度（默认只写 FOCUS_T）
    # ---- 2D κ 归一化（只影响报告的 κ，不影响寿命）----
    "DIM":             ("auto", "str"),        # auto | 2d | 3d
    "KAPPA_2D_THICKNESS": ("vdw", "str"),      # 无上游因子时现算用：vdw | 固定数值 Å | cell
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


def _warn(msg):
    print("[warn] " + msg, flush=True)


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


def _fmt(v, spec="%.6g"):
    """CSV 单元格式化：None/NaN -> 空串（旧版直接 % 到 None 会崩）。"""
    if v is None:
        return ""
    try:
        fv = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not np.isfinite(fv):
        return ""
    return spec % fv


def _local_hdf5_candidates(outdir, cwd):
    files = []
    for d in (outdir, cwd):
        for pat in ("kappa-m*.hdf5", "gamma-m*.hdf5", "kappa*.hdf5", "*.hdf5"):
            files.extend(sorted(Path(d).glob(pat)))
    return files


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

    # 本技能目录 / 步骤目录里已有的 hdf5（多候选时选主文件，绝不按字符串序取第一个）
    picked, _ignored = LC.select_kappa_hdf5(_local_hdf5_candidates(outdir, cwd),
                                            log=lambda m: print(m, flush=True),
                                            note="[local] ")
    if picked is not None and LC.is_main_hdf5(picked):
        return "kappa_hdf5", picked

    # 本地只有 --gp 单点文件时不急着用：先看看上游有没有多 q 点主文件，
    # 都没有再退回本地那个单点文件（旧版会直接拿本地单点文件去读，然后报缺 qpoint）。
    src = str(conf["SOURCE"] or "auto").lower()
    if src in ("auto", "upstream", "kappa_hdf5"):
        up = LC.find_kappa_hdf5(cwd, conf["UPSTREAM_SKILL"], conf["UPSTREAM_STEP"],
                                log=lambda m: print(m, flush=True))
        if up is not None:
            return "kappa_hdf5", up
    if picked is not None:
        return "kappa_hdf5", picked
    fc = LC.find_fc_dataset(cwd, conf["UPSTREAM_SKILL"], conf["UPSTREAM_STEP"])
    if fc:
        return "fc", fc
    return None, None


def _phono3py_cmd(conf, mesh, ts):
    """生成 phono3py --br 命令字符串（单独成函数便于回归测试）。

    含同位素时加 --isotope；NAC 不在这里加——phono3py-load 没有 --nac 选项
    （只有 --nonac），它会在目录里发现 BORN 时自动启用 NAC，所以由调用方负责
    把 BORN 拷进来（见 _run_phono3py）。
    """
    exe = str(conf["PHONO3PY_CMD"] or "phono3py-load")
    parts = [exe, "phono3py_disp.yaml", "--mesh", str(mesh), "--br",
             "--ts=\"%s\"" % str(ts)]
    if conf["INCLUDE_ISOTOPE"]:
        parts.append("--isotope")
    return " ".join(parts)


def _run_phono3py(conf, cwd, outdir, fc):
    """从 fc 数据集现跑 phono3py --br，返回生成的 hdf5 路径。"""
    outdir = Path(outdir)
    base = Path(fc[0])
    fc2, fc3, disp = (Path(fc[1]), Path(fc[2]), Path(fc[3]))
    for f in (fc2, fc3, disp):
        shutil.copyfile(str(f), str(outdir / f.name))
    born = Path(base) / "BORN"
    if born.is_file():
        shutil.copyfile(str(born), str(outdir / "BORN"))
        print("[..] 检测到 BORN，已拷入；phono3py-load 会自动启用 NAC（极性材料 LO 支依赖它）",
              flush=True)
    else:
        _warn("fc 目录没有 BORN：极性材料（如 225 相）的 LO 支会算错；非极性材料可忽略。")
    mesh = str(conf["PHONO3PY_MESH"] or "15 15 15")
    ts = str(conf["PHONO3PY_TS"] or "300")
    if not conf["INCLUDE_ISOTOPE"]:
        _warn("INCLUDE_ISOTOPE=false：现跑 phono3py 不加 --isotope，寿命不含同位素散射")
    cmd = _phono3py_cmd(conf, mesh, ts)
    print("[..] 现跑 phono3py：%s" % cmd, flush=True)
    cmd = cmd + " 2>&1 | tee phono3py_lifetime.log"
    py = sys.executable
    bindir = os.path.dirname(py)
    env_prefix = ("export PATH=%s:$PATH; " % bindir) if bindir else ""
    full = "cd %s && %s( %s )" % (str(outdir), env_prefix, cmd)
    try:
        rc = subprocess.run(["bash", "-lc", full],
                            timeout=int(conf["PHONO3PY_TIMEOUT"])).returncode
    except subprocess.TimeoutExpired:
        sys.exit("[ERROR] phono3py 超过 PHONO3PY_TIMEOUT=%ss 未完成" % conf["PHONO3PY_TIMEOUT"])
    if rc != 0:
        sys.exit("[ERROR] phono3py 退出码 %d，见 %s/phono3py_lifetime.log" % (rc, outdir))
    picked, _ = LC.select_kappa_hdf5(sorted(Path(outdir).glob("*.hdf5")),
                                     log=lambda m: print(m, flush=True), note="[phono3py run] ")
    if picked is not None:
        return picked
    sys.exit("[ERROR] phono3py 跑完但没找到含 gamma 的 hdf5")


def _src_kind(h5, cwd):
    """区分数据源：本技能目录内的本地文件 / 上游技能目录 / 外部显式路径。

    0.1 只用 startswith(cwd.parent)，把本技能 outdir 里的本地 hdf5 也标成 upstream。
    """
    hp = Path(h5).resolve()
    cwd_r = Path(cwd).resolve()
    if hp.parent == cwd_r or cwd_r in hp.parents:
        return "local_hdf5"
    if str(hp).startswith(str(cwd_r.parent)):
        return "upstream_hdf5"
    return "kappa_hdf5"


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
        src_kind = _src_kind(h5, cwd)

    data = LC.read_kappa_hdf5(h5)
    temperatures = data["temperatures"]
    frequencies = data["frequency"]
    gamma = data["gamma"]
    gamma_iso = data.get("gamma_isotope")
    if conf["INCLUDE_ISOTOPE"]:
        if gamma_iso is None or not np.any(np.asarray(gamma_iso) > 0):
            _warn("INCLUDE_ISOTOPE=true，但数据里没有 gamma_isotope（或全为 0）："
                  "寿命只含三声子散射；若上游跑 phono3py 时没加 --isotope，需要重跑上游。")
    else:
        print("[..] INCLUDE_ISOTOPE=false：只用三声子 gamma，寿命会偏长", flush=True)
        gamma_iso = None

    qpoints = data.get("qpoints")
    weights = data.get("weights")
    if qpoints is None or weights is None:
        sys.exit("[ERROR] %s 缺 qpoint/weight 数据集" % h5)

    tmap = LC.pick_temperature_indices(temperatures, _parse_nums(conf["TEMPERATURES"]),
                                       log=_warn)
    t_indices = sorted(tmap.values(), key=lambda i: temperatures[i])

    # 总和线宽 = gamma_anh + gamma_isotope (+ gamma_bd)。边界散射是否计入，由 hdf5
    # 里 kappa 的重建误差决定；phono3py 默认 boundary_mfp=1e6 µm 就是“关”（单位微米）。
    g_base = LC.effective_gamma(gamma, gamma_iso)
    bd_mfp = data.get("boundary_mfp")
    bd_mfp = float(np.asarray(bd_mfp).ravel()[0]) if bd_mfp is not None else None
    g_bd = LC.boundary_gamma(data.get("group_velocity"), bd_mfp)
    bd_use, bd_info = LC.decide_boundary_gamma(data, g_base, g_bd, t_indices=t_indices)
    g_bd_use = g_bd if bd_use else None
    gsum = (g_base + g_bd_use) if g_bd_use is not None else g_base
    if bd_info.get("boundary_mfp_um") is not None:
        _g = bd_info.get("gamma_bd_median_THz")
        if bd_use:
            bd_note = "计入边界散射 gamma_bd=|v|/(4*pi*L)"
            if _g is not None:
                bd_note += "，gamma_bd 中位=%.3g THz" % _g
        elif bd_info.get("method") == "mfp_sentinel":
            bd_note = "不计入（phono3py 默认 1e6 µm，视为关）"
        else:
            bd_note = "不计入（计入后 kappa 误差更大）"
        print("[..] boundary_mfp=%.6g µm：%s；kappa 重建误差 %s -> %s"
              % (bd_info["boundary_mfp_um"], bd_note,
                 ("%.2e" % bd_info["err_no_boundary"]) if bd_info["err_no_boundary"] is not None else "n/a",
                 ("%.2e" % bd_info["err_with_boundary"]) if bd_info["err_with_boundary"] is not None else "n/a"),
              flush=True)
    focus = float(conf["FOCUS_T"])
    focus_idx = int(np.argmin(np.abs(np.asarray(temperatures) - focus)))
    if abs(float(temperatures[focus_idx]) - focus) > 1.0:
        _warn("FOCUS_T=%.6g K 不在可用温度里，改用最近的 %.6g K（可用：%s）"
              % (focus, temperatures[focus_idx],
                 ", ".join("%g" % t for t in temperatures)))

    ac_fmax = conf["ACOUSTIC_FMAX"]
    ac_fmax = float(ac_fmax) if str(ac_fmax).strip() else None
    imag_thr = float(conf["IMAG_THR"])
    omega_thr = float(conf["OMEGA_TAU_THR"])

    per_T = [LC.summarize_temperature(frequencies, gamma, temperatures, ti,
                                      acoustic_fmax=ac_fmax, imag_thr=imag_thr,
                                      omega_tau_threshold=omega_thr,
                                      gamma_isotope=gamma_iso, gamma_extra=g_bd_use,
                                      weights=weights)
             for ti in t_indices]

    # ---- 逐模 CSV（默认只写 FOCUS_T；大胞密网格 × 全部温度会写出上千万行）----
    csv_mode = outdir / "lifetime_data.csv"
    csv_t = t_indices if conf["CSV_ALL_TEMPERATURES"] else [focus_idx]
    LC.write_mode_csv(str(csv_mode), temperatures, frequencies, gamma, qpoints, weights,
                      group_velocity=data.get("group_velocity"), acoustic_fmax=ac_fmax,
                      imag_thr=imag_thr, t_indices=csv_t, gamma_isotope=gamma_iso,
                      gamma_extra=g_bd_use)

    # ---- 每温度汇总 CSV ----
    csv_T = outdir / "lifetime_T.csv"
    with open(csv_T, "w", encoding="utf-8", newline="") as fh:
        import csv as _csv
        w = _csv.writer(fh)
        w.writerow(["T_K", "n_all", "n_all_weighted", "tau_median_all_ps",
                    "tau_median_all_anharm_ps", "tau_median_acoustic_ps",
                    "tau_median_optical_ps", "omega_tau_median_acoustic",
                    "overdamped_frac_acoustic", "acoustic_fmax_THz"])
        for d in per_T:
            a = d["acoustic"]
            w.writerow(["%.4f" % d["T_K"], d["all"].get("n_modes"),
                        _fmt(d["all"].get("n_modes_weighted"), "%.0f"),
                        _fmt(d["all"].get("tau_median_ps")),
                        _fmt(d["all"].get("tau_median_anharmonic_ps")),
                        _fmt(a.get("tau_median_ps")),
                        _fmt(d["optical"].get("tau_median_ps")),
                        _fmt(a.get("omega_tau_median")),
                        _fmt(a.get("overdamped_frac")),
                        _fmt(d["acoustic_fmax_THz"], "%.4f")])

    # ---- κ 重建自检（gamma 口径对不对，用上游存的 kappa 兜底验证）----
    kappa_check = LC.kappa_reconstruction_error(data, gsum, t_indices=t_indices)
    if kappa_check:
        print("[..] κ 重建自检：最大相对误差 %.2e（%d 个温度）"
              % (kappa_check["max_rel_err"], kappa_check["n_temperature"]), flush=True)
    if kappa_check and kappa_check["max_rel_err"] > 1e-3:
        _warn("κ 重建最大相对误差达 %.2e：上游 kappa 与 gamma 的口径可能不一致"
              "（例如它开了边界散射而我们没计入）" % kappa_check["max_rel_err"])

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
            rotations = None
            yp = LC.find_phonopy_yaml(h5, cwd, conf["UPSTREAM_SKILL"], conf["UPSTREAM_STEP"])
            if yp is not None:
                rotations = LC.symmetry_rotations_from_yaml(yp)
            if rotations is None:
                _warn("没找到可用的晶体对称操作（缺 phono3py.yaml / pyyaml / spglib）："
                      "沿 q 图只取几何上落在该方向的不可约点，等价方向会漏。")
            fn = LP.plot_lifetime_vs_q(str(outdir), temperatures, frequencies, gsum,
                                       qpoints, qdir, [focus_idx],
                                       direction_label=str(conf["Q_DIRECTION_LABEL"] or ""),
                                       imag_thr=imag_thr, rotations=rotations)
            if fn:
                plots.append(os.path.basename(fn))
    elif conf["MAKE_PLOT"] and not LP.HAVE_MPL:
        print("[warn] 无 matplotlib，跳过出图（JSON/CSV 照常产出）")

    # ---- κ 报告：raw + 2D 归一化（寿命与厚度无关，这里只管 κ 口径）----
    norm2d = LC.two_d_normalization(
        h5, data.get("mesh"), dim=str(conf["DIM"] or "auto"),
        thickness_mode=str(conf["KAPPA_2D_THICKNESS"] or "vdw"),
        yaml_path=LC.find_phonopy_yaml(h5, cwd, conf["UPSTREAM_SKILL"], conf["UPSTREAM_STEP"]),
        log=_warn)
    kblock = LC.kappa_block(data, focus, norm2d)
    if kblock and norm2d["dim"] == "2d":
        if kblock.get("kappa_2d_inplane_avg") is not None:
            print("[..] 2D：κ_raw(面内平均)=%.4g -> κ_2d=%.4g W/mK @%.0fK"
                  "（×h⊥/d=%.4f，d=%s Å，来源 %s）"
                  % (kblock["kappa_raw_inplane_avg"], kblock["kappa_2d_inplane_avg"],
                     kblock["T_K"], norm2d["factor"],
                     ("%.4g" % float(norm2d["thickness_d_A"]))
                     if norm2d.get("thickness_d_A") is not None else "?",
                     norm2d["source"]), flush=True)
        else:
            _warn("2D 材料只拿到 raw κ=%.4g W/mK（含真空层），不能直接对文献"
                  % kblock["kappa_raw_inplane_avg"])

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
                          "gamma_anh + gamma_isotope + gamma_bd  [phono3py kaccum/kdeplot + "
                          "conductivity/base.py:_get_main_diagonal / _get_boundary_scattering]",
        "tau_factor": LC.TAU_FACTOR,
        "isotope_included": bool(gamma_iso is not None and np.any(np.asarray(gamma_iso) > 0)),
        "boundary_mfp_um": bd_mfp,
        "boundary_gamma_included": bd_use,
        "boundary_decision": bd_info,
        "weights_applied": bool(per_T and per_T[0].get("weights_applied")),
        "acoustic_definition": ("lowest 3 branches by band index"
                               + (", f<=%.4g THz" % ac_fmax if ac_fmax is not None else "")),
        "temperatures_all_K": [float(x) for x in temperatures],
        "temperatures_analyzed_K": [float(temperatures[i]) for i in t_indices],
        "csv_temperatures_K": [float(temperatures[i]) for i in csv_t],
        "focus_T_K": float(temperatures[focus_idx]),
        "imag_thr_THz": imag_thr,
        "omega_tau_threshold": omega_thr,
        "per_temperature": per_T,
        "kappa_from_hdf5": kblock,
        "dim": norm2d["dim"],
        "kappa_2d_normalization": norm2d,
        "kappa_reconstruction": kappa_check,
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
             ("%.3g" % f["tau_median_ps"]) if f.get("tau_median_ps") else "n/a",
             ("%.3g" % f["omega_tau_median"]) if f.get("omega_tau_median") else "n/a"))


if __name__ == "__main__":
    main()
