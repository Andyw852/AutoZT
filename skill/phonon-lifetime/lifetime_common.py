#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lifetime_common.py -- 声子寿命/线宽分析公共库（供 phonon-lifetime 技能各步 import）。

物理约定（与 phono3py 官方后处理工具一致，勿改）
------------------------------------------------
phono3py 的 kappa-m*.hdf5 / gamma-m*.hdf5 里的 gamma 是声子自能的虚部
（linewidth），单位 THz（普通频率）。其内部把它当作 2*pi*gamma 的角频线宽使用，
因此寿命为::

    tau[ps] = 1 / (2 * 2 * pi * gamma_total[THz])

其中 gamma_total = gamma_anh + gamma_isotope（phono3py 的 kappa 用的就是总和，
见 conductivity/base.py:_get_main_diagonal；只取 gamma 会把寿命算长）。
即 tau = 1/(4*pi*gamma_total)。三处官方实现均为该式：
  * phono3py/scripts/phono3py_kdeplot.py:  tau = 1/g/(2*2*pi)
  * phono3py/other/kaccum.py:              mfp = v/(2*2*pi*g)（mfp = v*tau）
  * phono3py/conductivity/base.py:         get_unit_to_WmK() 里的 1/(2*pi)
    （注释 "2pi comes from definition of lifetime"），配合 RTA 的 /(g*2)。

注意 1 THz = 1e12/s = 1/ps，所以上式的 tau 直接以 ps 计。

寿命与“振动周期”的比较（论文 Kim et al. Nature 2021 图 2e 的判据）
---------------------------------------------------------------
声子周期 T_period[ps] = 1/f[THz]；无量纲 omega*tau = 2*pi*f*tau。
omega*tau ~ 1（或 tau ~ T_period）表示强散射/近过阻尼，声子在完成一个振荡
前就失去了相干性（Ioffe-Regel 判据）。
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

# phono3py gamma[THz] -> tau[ps] 的转换因子:  tau = 1/(TAU_FACTOR*gamma)
TAU_FACTOR = 2.0 * 2.0 * np.pi

# 判据默认阈值
DEFAULT_IMAG_THR = 0.10        # THz，低于此频率视为数值虚频/声学近零模，排除
DEFAULT_ACOUSTIC_FMAX = 2.0    # THz，声学支分析的频率上限兜底


# ---------------------------------------------------------------------------
# 数据集定位
# ---------------------------------------------------------------------------
def _iter_candidate_dirs(cwd, skill, step):
    """给出可能存放上游技能产物的目录（由近及远）。

    兼容两种真实布局：
      * 远端：<matdir>/<skill>/<step>/...
      * 本地拉回副本：<matdir>/<skill>/result/<step>/...
    以及技能步目录 cwd 位于 <matdir>/<skill>/<step> 或 <matdir>/<step> 的情况。
    """
    cwd = Path(cwd).resolve()
    seen = []
    for base in (cwd, cwd.parent, cwd.parent.parent, cwd.parent.parent.parent):
        cands = [
            base / skill / step,
            base / skill / "result" / step,
            base / skill,
            base / "result" / step,
            base,
        ]
        for c in cands:
            c = c.resolve()
            if c not in seen:
                seen.append(c)
    return seen


def find_kappa_hdf5(cwd, skill="kl-dft-cpu", step="step6_kappa"):
    """在材料目录树里找 phono3py 的 kappa/gamma hdf5（含 gamma 数据集）。

    优先级：kappa-m*.hdf5 > gamma-m*.hdf5 > kappa*.hdf5 > *.hdf5（含 gamma）。
    返回 Path 或 None。
    """
    for d in _iter_candidate_dirs(Path(cwd), skill, step):
        if not d.is_dir():
            continue
        for pat in ("kappa-m*.hdf5", "gamma-m*.hdf5", "kappa*.hdf5", "*.hdf5"):
            hits = sorted(d.glob(pat))
            if hits and hdf5_has_gamma(hits[-1]):
                return hits[-1]
    return None


def find_fc_dataset(cwd, skill="kl-dft-cpu", step="step6_kappa"):
    """找 phono3py 力常数数据集（fc2.hdf5 + fc3.hdf5 + phono3py_disp.yaml）。

    返回 (dir, fc2, fc3, disp_yaml) 或 None。
    """
    for d in _iter_candidate_dirs(Path(cwd), skill, step):
        if not d.is_dir():
            continue
        for sub in (".", "phono3py"):
            base = d / sub
            fc2 = base / "fc2.hdf5"
            fc3 = base / "fc3.hdf5"
            disp = base / "phono3py_disp.yaml"
            if fc2.is_file() and fc3.is_file() and disp.is_file():
                return base, fc2, fc3, disp
    return None


def hdf5_has_gamma(path):
    try:
        import h5py

        with h5py.File(path, "r") as f:
            return "gamma" in f and "frequency" in f
    except Exception:
        return False


# ---------------------------------------------------------------------------
# 读取与转换
# ---------------------------------------------------------------------------
def read_kappa_hdf5(path):
    """读 phono3py kappa/gamma hdf5，返回 dict（numpy 数组）。

    关键字段：
      temperatures (nT,), frequency (n_gp, n_band) [THz], gamma (nT, n_gp, n_band) [THz],
      qpoints (n_gp, 3) [fractional], weights (n_gp,), mesh (3,),
      group_velocity (n_gp, n_band, 3) [THz·Å], heat_capacity (nT, n_gp, n_band) [eV/K],
      kappa (nT, 6) [W/mK], kappa_unit_conversion (scalar)
    """
    import h5py

    out = {}
    with h5py.File(path, "r") as f:
        for key in ("temperature", "frequency", "gamma", "gamma_isotope", "qpoint",
                    "weight", "mesh", "group_velocity", "gv_by_gv", "heat_capacity",
                    "kappa", "mode_kappa", "grid_point", "kappa_unit_conversion",
                    "boundary_mfp", "version"):
            if key in f:
                out[key] = f[key][()]
    out["_file"] = str(path)
    if "temperature" in out:
        out["temperatures"] = np.asarray(out.pop("temperature"), dtype=float)
    if "qpoint" in out:
        out["qpoints"] = np.asarray(out.pop("qpoint"), dtype=float)
    if "weight" in out:
        out["weights"] = np.asarray(out.pop("weight"), dtype=float)
    out["frequency"] = np.asarray(out["frequency"], dtype=float)
    out["gamma"] = np.asarray(out["gamma"], dtype=float)
    return out


def effective_gamma(gamma, gamma_isotope=None):
    """总和线宽 g_total[THz] = gamma_anh + gamma_isotope。

    这不是可选装饰：phono3py 计算 kappa 时，碰撞矩阵对角元用的是
    gamma + gamma_isotope（conductivity/base.py:_get_main_diagonal），
    所以 hdf5 里的 kappa 含同位素散射。若只取 gamma，寿命会偏长。
    gamma_isotope 为空或无正元素时原样返回 gamma。
    """
    g = np.asarray(gamma, dtype=float)
    if gamma_isotope is None:
        return g
    gi = np.asarray(gamma_isotope, dtype=float)
    if gi.size == 0 or not np.any(gi > 0):
        return g
    return g + gi


def gamma_to_tau_ps(gamma_thz):
    """gamma[THz] -> tau[ps]（phono3py 约定，见模块 docstring）。gamma<=0 -> +inf。"""
    g = np.asarray(gamma_thz, dtype=float)
    safe = np.where(g > 0.0, g, 1.0)
    with np.errstate(divide="ignore"):
        tau = np.where(g > 0.0, 1.0 / (TAU_FACTOR * safe), np.inf)
    return tau


def tau_to_gamma_thz(tau_ps):
    tau_ps = np.asarray(tau_ps, dtype=float)
    return 1.0 / (TAU_FACTOR * tau_ps)


def mean_free_path_ang(group_velocity, gamma_thz):
    """mfp[Å] = |v|[THz·Å] * tau[ps]（= |v|/(2*2*pi*gamma)），与 phono3py kaccum 一致。"""
    gv = np.asarray(group_velocity, dtype=float)
    g = np.asarray(gamma_thz, dtype=float)
    vnorm = np.sqrt((gv ** 2).sum(axis=-1))
    tau = gamma_to_tau_ps(g)
    return np.where(np.isfinite(tau), vnorm * tau, 0.0)


# ---------------------------------------------------------------------------
# 温度 / 声学支选择
# ---------------------------------------------------------------------------
def pick_temperature_indices(temperatures, wanted=None):
    """给定温度数组，返回 {T: index}。wanted 为空取全部；否则取最接近的索引（去重）。"""
    temperatures = np.asarray(temperatures, dtype=float)
    if not wanted:
        return {float(t): i for i, t in enumerate(temperatures)}
    chosen = {}
    for t in wanted:
        i = int(np.argmin(np.abs(temperatures - float(t))))
        chosen[float(temperatures[i])] = i
    return chosen


def acoustic_mask(frequencies, gamma, acoustic_fmax=None, imag_thr=DEFAULT_IMAG_THR):
    """声学支掩码。

    简单且稳健：最低 3 条支（index 0/1/2）视为声学支，再叠加频率窗口。
    2D/含真空材料真空方向有近零支，仍落在最低 3 条里。若 acoustic_fmax 为空，
    用“最低 3 条支中正频模的中位频率 * 2”兜底。
    返回 (mask, fmax_used)。
    """
    f = np.asarray(frequencies, dtype=float)   # (n_gp, n_band)
    g = np.asarray(gamma, dtype=float)
    n_band = f.shape[1]
    n_ac = min(3, n_band)
    branch = np.zeros_like(f, dtype=bool)
    branch[:, :n_ac] = True
    good = (f > imag_thr) & (g > 0)
    if acoustic_fmax is None:
        vals = f[branch & good]
        fmax = float(np.median(vals) * 2.0) if vals.size else DEFAULT_ACOUSTIC_FMAX
        fmax = max(fmax, DEFAULT_ACOUSTIC_FMAX)
    else:
        fmax = float(acoustic_fmax)
    return branch & good & (f <= fmax), fmax


# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
def _stats(tau, omega_tau):
    tau = np.asarray(tau, dtype=float)
    tau = tau[np.isfinite(tau)]
    if tau.size == 0:
        return {"n_modes": 0}
    ot = np.asarray(omega_tau, dtype=float)
    ot = ot[np.isfinite(ot)]
    return {
        "n_modes": int(tau.size),
        "tau_median_ps": float(np.median(tau)),
        "tau_mean_ps": float(np.mean(tau)),
        "tau_p10_ps": float(np.percentile(tau, 10)),
        "tau_p90_ps": float(np.percentile(tau, 90)),
        "omega_tau_median": float(np.median(ot)) if ot.size else None,
    }


def summarize_temperature(frequencies, gamma, temperatures, t_index,
                          acoustic_fmax=None, imag_thr=DEFAULT_IMAG_THR,
                          omega_tau_threshold=1.0, gamma_isotope=None):
    """对单个温度汇总：全体 / 声学 / 光学 分组的 tau、omega*tau 与过阻尼比例。

    gamma_isotope 非空时，tau 取总和线宽 (gamma + gamma_isotope)，与 phono3py 的
    kappa 一致；同时给出只看三声子（anharmonic）的 tau_median_anharmonic_ps。
    """
    f = np.asarray(frequencies, dtype=float)
    g_anh = np.asarray(gamma[t_index], dtype=float)
    g = effective_gamma(g_anh, gamma_isotope)
    tau = gamma_to_tau_ps(g)
    tau_anh = gamma_to_tau_ps(g_anh)
    omega_tau = 2.0 * np.pi * f * tau
    isotope = bool(gamma_isotope is not None and np.any(np.asarray(gamma_isotope) > 0))

    good = (f > imag_thr) & (g > 0)
    ac_mask, fmax = acoustic_mask(f, g, acoustic_fmax, imag_thr)
    opt_mask = good & ~ac_mask

    def pack(mask):
        d = _stats(tau[mask], omega_tau[mask])
        if mask.sum():
            finite = np.isfinite(omega_tau[mask])
            d["overdamped_frac"] = float(np.mean(omega_tau[mask][finite] < omega_tau_threshold)) \
                if finite.any() else None
            d["f_min_THz"] = float(f[mask].min())
            d["f_max_THz"] = float(f[mask].max())
            ta = tau_anh[mask]
            ta = ta[np.isfinite(ta)]
            d["tau_median_anharmonic_ps"] = float(np.median(ta)) if ta.size else None
        else:
            d["overdamped_frac"] = None
            d["f_min_THz"] = None
            d["f_max_THz"] = None
            d["tau_median_anharmonic_ps"] = None
        return d

    return {
        "T_K": float(temperatures[t_index]),
        "isotope_included": isotope,
        "all": pack(good),
        "acoustic": pack(ac_mask),
        "optical": pack(opt_mask),
        "acoustic_fmax_THz": float(fmax),
        "omega_tau_threshold": float(omega_tau_threshold),
    }


# ---------------------------------------------------------------------------
# 沿指定 q 方向取样（论文 Γ-A 那类“寿命 vs q”）
# ---------------------------------------------------------------------------
def select_direction_qpoints(qpoints, direction, n_points=21, tol=0.06):
    """从网格 q 点里挑出最接近 t*direction (t∈[0,1]) 的点。

    返回 list[(t, gp_index)]，按 t 升序，已对 t 去重。
    """
    q = np.asarray(qpoints, dtype=float)
    d = np.asarray(direction, dtype=float)
    nrm = np.linalg.norm(d)
    d = d / nrm if nrm else d
    qq = q - np.round(q)          # 映射到 [-0.5, 0.5)
    dist = np.linalg.norm(qq - np.outer(qq @ d, d), axis=1)
    proj = qq @ d
    picked, used = [], set()
    for k in range(n_points + 1):
        t = k / float(n_points)
        mask = (dist < tol) & (np.abs(proj - t) < 0.5 / n_points)
        idx = np.where(mask)[0]
        if idx.size == 0:
            continue
        j = int(idx[int(np.argmin(np.abs(proj[idx] - t)))])
        key = round(float(t), 6)
        if key in used:
            continue
        used.add(key)
        picked.append((float(t), j))
    picked.sort()
    return picked


# ---------------------------------------------------------------------------
# 落盘
# ---------------------------------------------------------------------------
def write_mode_csv(path, temperatures, frequencies, gamma, qpoints, weights,
                   group_velocity=None, acoustic_fmax=None,
                   imag_thr=DEFAULT_IMAG_THR, t_indices=None, gamma_isotope=None):
    """逐模写 CSV。

    tau_ps 用总和线宽 (gamma + gamma_isotope)，与 phono3py 的 kappa 一致；
    tau_anh_ps 只用三声子 gamma。mfp 用逐模群速度（不是 q 点平均）。
    """
    f = np.asarray(frequencies, dtype=float)
    g_all = np.asarray(gamma, dtype=float)
    gi_all = None if gamma_isotope is None else np.asarray(gamma_isotope, dtype=float)
    q = np.asarray(qpoints, dtype=float)
    w = np.asarray(weights, dtype=float)
    gv = None if group_velocity is None else np.asarray(group_velocity, dtype=float)
    n_gp, n_band = f.shape
    if t_indices is None:
        t_indices = range(len(temperatures))

    with open(path, "w", newline="", encoding="utf-8") as fh:
        csvw = csv.writer(fh)
        csvw.writerow(["T_K", "gp", "qx", "qy", "qz", "weight", "band",
                       "f_THz", "gamma_THz", "gamma_iso_THz", "gamma_total_THz",
                       "tau_ps", "tau_anh_ps", "omega_tau", "period_ps",
                       "v_THzA", "mfp_A", "acoustic"])
        for ti in t_indices:
            g_anh = g_all[ti]
            if gi_all is None:
                gi = None
            elif gi_all.ndim == 2:
                gi = gi_all
            else:
                gi = gi_all[ti]
            g = effective_gamma(g_anh, gi)
            tau = gamma_to_tau_ps(g)
            tau_anh = gamma_to_tau_ps(g_anh)
            omega_tau = 2.0 * np.pi * f * tau
            ac_mask, _ = acoustic_mask(f, g, acoustic_fmax, imag_thr)
            for gp in range(n_gp):
                vnorm_b = None
                if gv is not None and gv.ndim == 3:
                    # 逐条支的群速度模，不能对支取平均（旧版 bug：mfp 用 q 点平均 v）
                    vnorm_b = np.sqrt((gv[gp, :, :] ** 2).sum(axis=-1))
                for b in range(n_band):
                    period = 1.0 / f[gp, b] if f[gp, b] > 0 else 0.0
                    vb = float(vnorm_b[b]) if vnorm_b is not None else None
                    mfp = (vb * tau[gp, b]) if (vb is not None and g[gp, b] > 0) else 0.0
                    csvw.writerow([
                        "%.4f" % temperatures[ti], gp,
                        "%.6f" % q[gp, 0], "%.6f" % q[gp, 1], "%.6f" % q[gp, 2],
                        "%.6f" % w[gp], b,
                        "%.6f" % f[gp, b], "%.6e" % g_anh[gp, b],
                        "%.6e" % (float(gi[gp, b]) if gi is not None else 0.0),
                        "%.6e" % g[gp, b],
                        ("%.6g" % tau[gp, b]) if np.isfinite(tau[gp, b]) else "inf",
                        ("%.6g" % tau_anh[gp, b]) if np.isfinite(tau_anh[gp, b]) else "inf",
                        ("%.6g" % omega_tau[gp, b]) if np.isfinite(omega_tau[gp, b]) else "inf",
                        "%.6g" % period,
                        ("%.6f" % vb) if vb is not None else "",
                        "%.6g" % mfp,
                        int(bool(ac_mask[gp, b])),
                    ])


def dump_json(path, obj):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, default=_json_default)


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


__all__ = [
    "TAU_FACTOR", "DEFAULT_IMAG_THR", "DEFAULT_ACOUSTIC_FMAX",
    "find_kappa_hdf5", "find_fc_dataset", "hdf5_has_gamma",
    "read_kappa_hdf5", "effective_gamma", "gamma_to_tau_ps", "tau_to_gamma_thz",
    "mean_free_path_ang",
    "pick_temperature_indices", "acoustic_mask", "summarize_temperature",
    "select_direction_qpoints", "write_mode_csv", "dump_json",
]
