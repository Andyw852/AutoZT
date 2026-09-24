#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lifetime_plots.py -- phonon-lifetime 技能的出图（matplotlib Agg，无则跳过）。

所有图都以 phono3py 约定的 tau[ps] = 1/(2*2*pi*gamma_total[THz]) 为准；
gamma_total = gamma_anh + gamma_isotope（调用方传进来的 gamma 数组即总和线宽）。
"""
from __future__ import annotations

import os

import numpy as np

import lifetime_common as LC

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:            # pragma: no cover - 无 matplotlib 环境
    HAVE_MPL = False


def _ensure():
    if not HAVE_MPL:
        return False
    return True


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def plot_lifetime_vs_frequency(outdir, temperatures, frequencies, gamma,
                               t_indices, acoustic_fmax=None, fmax_plot=20.0,
                               imag_thr=LC.DEFAULT_IMAG_THR):
    if not _ensure():
        return None
    f = np.asarray(frequencies, dtype=float)
    g = np.asarray(gamma, dtype=float)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    cmap = plt.get_cmap("viridis")
    n = max(1, len(t_indices) - 1)
    for k, ti in enumerate(t_indices):
        tau = LC.gamma_to_tau_ps(g[ti])
        mask = (f > imag_thr) & (g[ti] > 0) & (f <= fmax_plot) & np.isfinite(tau)
        ax.scatter(f[mask], tau[mask], s=6, alpha=0.55,
                   color=cmap(k / n), label="%.0f K" % temperatures[ti])
    ax.set_yscale("log")
    ax.set_xlabel("frequency (THz)")
    ax.set_ylabel(r"lifetime $\tau$ (ps)")
    ax.set_title("Phonon lifetime vs frequency")
    ax.grid(True, which="both", alpha=0.25)
    if len(t_indices) <= 12:
        ax.legend(fontsize=8, ncol=2)
    return _save(fig, os.path.join(outdir, "lifetime_vs_frequency.png"))


def plot_lifetime_vs_period(outdir, temperatures, frequencies, gamma,
                            t_indices, fmax_plot=20.0, imag_thr=LC.DEFAULT_IMAG_THR):
    """tau vs 振动周期 1/f，含 tau = T_period（omega*tau=2*pi）参考线。"""
    if not _ensure():
        return None
    f = np.asarray(frequencies, dtype=float)
    g = np.asarray(gamma, dtype=float)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    cmap = plt.get_cmap("plasma")
    n = max(1, len(t_indices) - 1)
    for k, ti in enumerate(t_indices):
        tau = LC.gamma_to_tau_ps(g[ti])
        mask = (f > imag_thr) & (g[ti] > 0) & (f <= fmax_plot) & np.isfinite(tau)
        period = 1.0 / f[mask]
        ax.scatter(period, tau[mask], s=6, alpha=0.55,
                   color=cmap(k / n), label="%.0f K" % temperatures[ti])
    xs = np.array([1.0 / fmax_plot, 1.0 / max(imag_thr, 1e-3)])
    ax.plot(xs, xs, "k--", lw=1.0, label=r"$\tau = T_{\rm period}$ ($\omega\tau=2\pi$)")
    ax.plot(xs, xs / (2 * np.pi), "k:", lw=1.0, label=r"$\omega\tau=1$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"vibrational period $T=1/f$ (ps)")
    ax.set_ylabel(r"lifetime $\tau$ (ps)")
    ax.set_title("Lifetime vs vibrational period (overdamping chart)")
    ax.grid(True, which="both", alpha=0.25)
    if len(t_indices) <= 12:
        ax.legend(fontsize=8)
    return _save(fig, os.path.join(outdir, "lifetime_vs_period.png"))


def plot_ioffe_regel(outdir, temperatures, frequencies, gamma, t_indices,
                     fmax_plot=20.0, imag_thr=LC.DEFAULT_IMAG_THR):
    if not _ensure():
        return None
    f = np.asarray(frequencies, dtype=float)
    g = np.asarray(gamma, dtype=float)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    cmap = plt.get_cmap("cividis")
    n = max(1, len(t_indices) - 1)
    for k, ti in enumerate(t_indices):
        tau = LC.gamma_to_tau_ps(g[ti])
        ot = 2.0 * np.pi * f * tau
        mask = (f > imag_thr) & (g[ti] > 0) & (f <= fmax_plot) & np.isfinite(ot)
        ax.scatter(f[mask], ot[mask], s=6, alpha=0.55,
                   color=cmap(k / n), label="%.0f K" % temperatures[ti])
    ax.axhline(1.0, color="r", ls="--", lw=1.0, label=r"$\omega\tau=1$ (Ioffe-Regel)")
    ax.set_yscale("log")
    ax.set_xlabel("frequency (THz)")
    ax.set_ylabel(r"$\omega\tau$")
    ax.set_title("Ioffe-Regel parameter")
    ax.grid(True, which="both", alpha=0.25)
    if len(t_indices) <= 12:
        ax.legend(fontsize=8)
    return _save(fig, os.path.join(outdir, "ioffe_regel.png"))


def plot_lifetime_T(outdir, per_temperature):
    if not _ensure():
        return None
    Ts = [d["T_K"] for d in per_temperature]
    series = [("all", "tau_median_ps", "tab:blue"),
              ("acoustic", "tau_median_ps", "tab:green"),
              ("optical", "tau_median_ps", "tab:orange")]
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    for key, field, color in series:
        ys = [d.get(key, {}).get(field) for d in per_temperature]
        ys = [np.nan if y is None else y for y in ys]   # 空分组不能让 matplotlib 收到 None
        ax.plot(Ts, ys, "o-", color=color, label=key)
    ax.set_yscale("log")
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel(r"median lifetime $\tau$ (ps)")
    ax.set_title("Phonon lifetime vs temperature")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=9)
    return _save(fig, os.path.join(outdir, "lifetime_vs_T.png"))


def plot_lifetime_vs_q(outdir, temperatures, frequencies, gamma, qpoints,
                       direction, t_indices, direction_label="",
                       imag_thr=LC.DEFAULT_IMAG_THR, rotations=None):
    """沿指定 q 方向的逐支寿命（论文 Γ-A 图 2e 的 DFT 版）。

    rotations（分数坐标下的对称操作）给了才能把不可约点的等价方向映射回目标线上。
    """
    if not _ensure():
        return None
    picked = LC.select_direction_qpoints(qpoints, direction, n_points=30,
                                         rotations=rotations)
    if not picked:
        return None
    f = np.asarray(frequencies, dtype=float)
    g = np.asarray(gamma, dtype=float)
    n_band = f.shape[1]
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    cmap = plt.get_cmap("tab10")
    ti = t_indices[0]
    for b in range(n_band):
        xs, ys = [], []
        for t, gp in picked:
            if f[gp, b] > imag_thr and g[ti][gp, b] > 0:
                xs.append(t)
                ys.append(LC.gamma_to_tau_ps(g[ti][gp, b]))
        if xs:
            ax.plot(xs, ys, "o-", ms=3, lw=1.0, color=cmap(b % 10), label="band %d" % (b + 1))
    ax.set_yscale("log")
    ax.set_xlabel("reduced q along %s" % (direction_label or "direction"))
    ax.set_ylabel(r"lifetime $\tau$ (ps)")
    ax.set_title("Lifetime vs q  (T=%.0f K)" % temperatures[ti])
    ax.grid(True, which="both", alpha=0.25)
    if n_band <= 12:
        ax.legend(fontsize=7, ncol=2)
    return _save(fig, os.path.join(outdir, "lifetime_vs_q.png"))


__all__ = ["HAVE_MPL", "plot_lifetime_vs_frequency", "plot_lifetime_vs_period",
           "plot_ioffe_regel", "plot_lifetime_T", "plot_lifetime_vs_q"]
