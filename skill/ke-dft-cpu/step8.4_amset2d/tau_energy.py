#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tau_energy.py —— [0035] 从 AMSET 的 mesh.h5 输出电子弛豫时间 τ（S8.4 作业链里跑，写 tau_energy.json）。

两样东西：
1. τ(E)：逐机制、按能量分箱的弛豫时间。载流子能量从带边量起（电子 x = E − CBM，空穴 x = VBM − E），
   每个能量箱里所有 k 点、所有能带的散射率取平均，τ_m(E) = 1/⟨Γ_m⟩_E（即箱内 τ 的调和平均）；
   总的 τ(E) = 1/⟨Σ_m Γ_m⟩_E。只用 AMSET 真正算过散射的态（FD 截断窗口内），窗口外的率是 AMSET 补的，不进统计。
2. 迁移率等效 τ（τ_eff）：让"常数 τ 的玻尔兹曼输运"复现 AMSET 迁移率的那个 τ。
       τ_eff,m = τ0 · μ_m / μ_CRTA(τ0)
   μ_m 是 AMSET 对机制 m 的迁移率（与 transport.json 里的单机制值同一个函数 _calculate_mobility），
   μ_CRTA(τ0) 是同一套插值能带、同一费米能级、所有态 τ = τ0 时的迁移率（同一个函数）。
   这正是 S8.1（BoltzTraP2 常数 τ）里用 S8.2 DPT 的 τ 所扮演的角色 —— τ_eff,ADP 与 S8.2 的 τ_DPT 是同口径对比。
   分方向给 xx、yy 与面内（(μxx+μyy)/(μc,xx+μc,yy)）。
核对（写进 check 段）：
   · 用 mesh 散射率重算的 overall 迁移率与 transport.json 一致（同 postprocess_intrinsic 的复现检查）；
   · 由逐 k 点 τ 直接算的电导加权平均 τ_σ = Σ τ v²(−∂f/∂E) / Σ v²(−∂f/∂E) 与 τ_eff,overall 相近
     （BTE 走四面体 + 能量网格积分，差几个百分点正常；差得大说明载流子/能带选错了）。

用法（S8.4 运行目录里，作业链自动跑；人工补跑同样）：
    python tau_energy.py [运行目录] [--mesh mesh_xxx.h5] [--out tau_energy.json]
失败时也写 tau_energy.json（status = error + 原因），退出码非 0；作业链对它是"尽力而为"，不影响 transport.json。
"""
import argparse
import json
import os
import sys
import traceback
from pathlib import Path

import numpy as np

__version__ = "2026-10-11-tau-v1"
E_MAX_EV = 0.5           # τ(E) 从带边往上统计多宽
DE_EV = 0.01             # 能量箱宽
TAU0_S = 1.0e-14         # 常数 τ 参考值（10 fs）；μ 与 τ 成正比，取值不影响结果
KB_HA = 3.166811563e-6   # Boltzmann 常数，Hartree/K
HA_EV = 27.211386245988
SIGMA_REL_WARN = 0.10    # τ_σ 与 τ_eff 相对差超过它告警


# =====================================================================
# 纯函数（不依赖 amset，可单测）
# =====================================================================
def carrier_states(energies, vb_idx, n_type):
    """energies：{自旋: (nb, nk) Hartree}；vb_idx：{自旋: 最高价带下标}。
    返回 ({自旋: 该载流子的能带下标数组}, 带边能量 Ha)。n 型取导带、带边 = CBM；p 型取价带、带边 = VBM（与 AMSET 的
    _get_band_idx 同一选法）。"""
    sel, edge = {}, None
    for s, e in energies.items():
        nb = len(e)
        cb = int(vb_idx[s]) + 1
        idx = np.arange(cb, nb) if n_type else np.arange(0, cb)
        sel[s] = idx
        if len(idx):
            v = float(np.min(e[idx])) if n_type else float(np.max(e[idx]))
            edge = v if edge is None else (min(edge, v) if n_type else max(edge, v))
    return sel, edge


def bin_tau(x_ev, rates, edges):
    """x_ev (N,)：载流子能量（从带边量起，eV）；rates (nm, N)：各机制散射率（1/s）；edges：箱边界（eV）。
    返回 (τ_m (nm, nbin) 秒、总 τ (nbin,) 秒、各箱态数 (nbin,))；空箱为 nan。"""
    x = np.asarray(x_ev, float)
    r = np.atleast_2d(np.asarray(rates, float))
    nb = len(edges) - 1
    ib = np.searchsorted(edges, x, side="right") - 1
    ok = (ib >= 0) & (ib < nb)
    cnt = np.bincount(ib[ok], minlength=nb).astype(float)
    tau = np.full((r.shape[0], nb), np.nan)
    tot = np.full(nb, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        for m in range(r.shape[0]):
            ssum = np.bincount(ib[ok], weights=r[m, ok], minlength=nb)
            mean = ssum / cnt
            tau[m] = np.where(cnt > 0, 1.0 / mean, np.nan)
        ssum = np.bincount(ib[ok], weights=r[:, ok].sum(axis=0), minlength=nb)
        tot = np.where(cnt > 0, 1.0 / (ssum / cnt), np.nan)
    return tau, tot, cnt


def tau_sigma(tau_s, v2, w):
    """电导加权平均 τ：Σ τ v² w / Σ v² w（w = f(1−f) ∝ −∂f/∂E）。分母为 0 -> nan。"""
    num = float(np.sum(np.asarray(tau_s) * np.asarray(v2) * np.asarray(w)))
    den = float(np.sum(np.asarray(v2) * np.asarray(w)))
    return num / den if den > 0 else float("nan")


def tau_eff(mu, mu_c, tau0=TAU0_S):
    """mu、mu_c：(..., 3, 3) 迁移率张量。返回 {"xx","yy","inplane"}：τ0·μ/μ_c（秒），面内 = τ0·(μxx+μyy)/(μc,xx+μc,yy)。"""
    mu = np.asarray(mu, float)
    mc = np.asarray(mu_c, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return {"xx": tau0 * mu[..., 0, 0] / mc[..., 0, 0],
                "yy": tau0 * mu[..., 1, 1] / mc[..., 1, 1],
                "inplane": tau0 * (mu[..., 0, 0] + mu[..., 1, 1]) / (mc[..., 0, 0] + mc[..., 1, 1])}


def fermi_weight(e_ha, mu_ha, T):
    """f(1−f)，数值稳定写法。"""
    z = (np.asarray(e_ha, float) - mu_ha) / (KB_HA * max(float(T), 1e-9))
    z = np.clip(z, -200.0, 200.0)
    return 0.25 / np.cosh(z / 2.0) ** 2


def _js(a, nd=6):
    """nan/inf -> None，其余按有效数字取整（json 友好）。"""
    a = np.asarray(a, float)
    out = np.where(np.isfinite(a), a, np.nan)
    lst = out.tolist()

    def fix(v):
        if isinstance(v, list):
            return [fix(x) for x in v]
        return None if v != v else float("%.*g" % (nd, v))
    return fix(lst)


# =====================================================================
# 主流程（要 amset）
# =====================================================================
def analyze(ad, rates_full, labels, mobility_ref=None):
    """核心计算（要 amset 的 AmsetData，已设好费米能级/温度；不读文件）。
    rates_full：{自旋: (机制, 掺杂, 温度, 带, 全网格 k)} 散射率（1/s），labels：机制名；
    mobility_ref：原运行的 overall 迁移率 (掺杂, 温度, 3, 3)（给了就做复现核对，不一致抛错）。返回结果字典（不含文件元信息）。"""
    from amset.constants import bohr_to_cm
    from amset.core.transport import _calculate_mobility

    if ad.is_metal:
        raise RuntimeError("AMSET 判为金属（is_metal），不给迁移率，τ_eff 无从定义")
    labels = list(labels)
    nm = len(labels)
    ndop, nT = ad.fermi_levels.shape
    doping = np.asarray(ad.doping, float) * (1 / bohr_to_cm) ** 3
    temps = np.asarray(ad.temperatures, float)
    ad.set_scattering_rates({s: rates_full[s] for s in rates_full}, labels)

    # 1) 迁移率：overall（复现核对）+ 各机制
    mu_all = np.asarray(_calculate_mobility(ad, np.arange(nm), pbar_label=None))
    check = {}
    if mobility_ref is not None:
        ref = np.asarray(mobility_ref, float)
        if ref.shape != mu_all.shape:
            raise RuntimeError("原运行迁移率形状 %s 与重算 %s 不同" % (ref.shape, mu_all.shape))
        rel = float(np.nanmax(np.abs(mu_all - ref) / np.maximum(np.abs(ref), 1e-12)))
        check["mobility_reproduce_max_rel"] = rel
        if rel > 1e-6:
            raise RuntimeError("用 mesh 散射率重算的 overall 迁移率与原运行不一致（最大相对差 %.2e）——"
                               "mesh/settings/vasprun 不同源，τ 不可信" % rel)
    zero = {lab for i, lab in enumerate(labels) if all(np.all(np.asarray(rates_full[s])[i] == 0) for s in rates_full)}
    mu_m = {lab: np.asarray(_calculate_mobility(ad, i, pbar_label=None))
            for i, lab in enumerate(labels) if lab not in zero}
    # 2) 常数 τ 的迁移率（同一函数、同一费米能级）
    const = {s: np.full((1,) + np.asarray(rates_full[s]).shape[1:], 1.0 / TAU0_S) for s in rates_full}
    ad.set_scattering_rates(const, ["CRTA"])
    mu_c = np.asarray(_calculate_mobility(ad, [0], pbar_label=None))
    ad.set_scattering_rates({s: rates_full[s] for s in rates_full}, labels)
    teff = {"overall": tau_eff(mu_all, mu_c)}
    for lab in labels:
        if lab in mu_m:
            teff[lab] = tau_eff(mu_m[lab], mu_c)

    # 3) τ(E) 与 τ_σ
    edges = np.arange(0.0, E_MAX_EV + DE_EV / 2, DE_EV)
    lo, hi = ad.fd_cutoffs if ad.fd_cutoffs is not None else (-np.inf, np.inf)
    tauE = {lab: np.full((ndop, nT, len(edges) - 1), np.nan) for lab in labels}
    tauE_tot = np.full((ndop, nT, len(edges) - 1), np.nan)
    counts, carrier, tsig = {}, [], np.full((ndop, nT), np.nan)
    for n in range(ndop):
        n_type = bool(doping[n] < 0)
        carrier.append("electron" if n_type else "hole")
        sel, edge = carrier_states(ad.energies, ad.vb_idx, n_type)
        xs, rs, es, v2s = [], [], [], []
        for s in ad.spins:
            b = sel[s]
            if not len(b):
                continue
            e = np.asarray(ad.energies[s])[b]                               # (nbs, nk) Ha
            inwin = (e >= lo) & (e <= hi)
            x = ((e - edge) if n_type else (edge - e)) * HA_EV
            r = np.asarray(rates_full[s])[:, n][:, :, b]                    # (nm, nT, nbs, nk)
            v = np.asarray(ad.velocities[s])[b]                              # (nbs, nk, 3)
            xs.append(x[inwin])
            rs.append(r[:, :, inwin])
            es.append(e[inwin])
            v2s.append((v[..., 0] ** 2 + v[..., 1] ** 2)[inwin])
        x = np.concatenate(xs)
        r = np.concatenate(rs, axis=2)                                        # (nm, nT, N)
        e = np.concatenate(es)
        v2 = np.concatenate(v2s)
        for t in range(nT):
            tm, tt, cnt = bin_tau(x, r[:, t], edges)
            for i, lab in enumerate(labels):
                tauE[lab][n, t] = tm[i]
            tauE_tot[n, t] = tt
            counts[carrier[-1]] = cnt
            w = fermi_weight(e, float(ad.fermi_levels[n, t]), temps[t])
            with np.errstate(divide="ignore"):
                tsig[n, t] = tau_sigma(1.0 / r[:, t].sum(axis=0), v2, w)
    rel_sig = np.abs(tsig - teff["overall"]["inplane"]) / np.abs(teff["overall"]["inplane"])
    check["tau_sigma_vs_eff_max_rel"] = float(np.nanmax(rel_sig))
    i300 = int(np.argmin(np.abs(temps - 300.0)))
    check["tau_sigma_vs_eff_max_rel_%dK" % round(temps[i300])] = float(np.nanmax(rel_sig[:, i300]))
    notes = ["τ(E) = 1/⟨Γ⟩_E：箱内各态散射率的平均再取倒数（= τ 的调和平均）；只统计 FD 截断窗口内（AMSET 真算过散射的态）",
             "τ_eff = τ0·μ/μ_CRTA(τ0)：常数 τ 的玻尔兹曼输运要复现 AMSET 迁移率所需的 τ；τ_eff,ADP 与 S8.2 DPT 的 τ 同口径",
             "τ_σ = Σ τ v²(−∂f/∂E)/Σ v²(−∂f/∂E)（逐 k 点、面内 v²）只作核对：BTE 走四面体积分，与 τ_eff 差几个百分点正常"]
    if check["tau_sigma_vs_eff_max_rel"] > SIGMA_REL_WARN:
        notes.append("★ τ_σ 与 τ_eff,overall 最大相对差 %.0f%%（> %.0f%%）：多半在低温/高掺杂端（BTE 能量网格比 kT 粗），"
                     "核对 300 K 的那几行" % (100 * check["tau_sigma_vs_eff_max_rel"], 100 * SIGMA_REL_WARN))
    to_fs = 1e15
    return {
        "status": "ok", "version": __version__,
        "temperatures_K": temps.tolist(), "doping_cm3": doping.tolist(), "carrier": carrier,
        "mechanisms": labels, "zero_mechanisms": sorted(zero),
        "energy_bins_eV": {"edges": _js(edges), "centers": _js((edges[:-1] + edges[1:]) / 2)},
        "n_states": {k: v.astype(int).tolist() for k, v in counts.items()},
        "tau_E_fs": {lab: _js(tauE[lab] * to_fs, 5) for lab in labels},
        "tau_E_total_fs": _js(tauE_tot * to_fs, 5),
        "tau_eff_fs": {k: {d: _js(np.asarray(v[d]) * to_fs, 6) for d in ("xx", "yy", "inplane")}
                       for k, v in teff.items()},
        "tau_sigma_fs": _js(tsig * to_fs, 6),
        "mobility_crta_cm2_Vs": {"tau0_fs": TAU0_S * to_fs,
                                 "inplane": _js((mu_c[..., 0, 0] + mu_c[..., 1, 1]) / 2.0, 6)},
        "check": check, "notes": notes,
    }


def run(run_dir, mesh_arg=None):
    """S8.4 运行目录 -> 结果字典：用 postprocess_intrinsic 的同一套重建（与作业同一插值、不可约映射、费米能级），
    注入 mesh 的逐机制散射率，再 analyze。"""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import postprocess_intrinsic as pp

    mesh_file = pp.find_mesh(run_dir, mesh_arg)
    mesh = pp.read_mesh_h5(mesh_file)
    ir_fix = pp.mesh_ir_fix_state(mesh, pp._settings_symprec(run_dir))
    ad, settings, amset_version = pp.build_amset_data(run_dir, progress_bar=False, ir_fix=ir_fix)
    labels = list(mesh["scattering_labels"])
    perm = pp.apply_mesh_metadata(ad, mesh)
    rates_full, labels = pp.inject_mesh_rates(ad, mesh, labels, perm=perm)
    ref, ref_name = None, None
    tf = sorted(Path(run_dir).glob("transport*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
    if tf:
        ref = (json.loads(tf[0].read_text()).get("mobility") or {}).get("overall")
        ref_name = tf[0].name
    res = analyze(ad, rates_full, labels, mobility_ref=ref)
    res["check"]["mobility_reproduce_file"] = ref_name
    res.update(amset_version=amset_version, mesh_file=Path(mesh_file).name,
               transport_fingerprint=pp.transport_fingerprint(run_dir))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description="AMSET mesh -> τ(E) 与迁移率等效 τ（tau_energy.json）")
    ap.add_argument("run_dir", nargs="?", default=".")
    ap.add_argument("--mesh", default=None)
    ap.add_argument("--out", default="tau_energy.json")
    a = ap.parse_args(argv)
    run_dir = Path(a.run_dir).resolve()
    out = Path(a.out)
    if not out.is_absolute():
        out = run_dir / out
    try:
        res = run(run_dir, a.mesh)
        rc = 0
    except Exception as e:                                       # noqa: BLE001
        res = {"status": "error", "version": __version__, "error": "%s: %s" % (type(e).__name__, e),
               "traceback": traceback.format_exc()[-4000:]}
        rc = 1
        print("[WARN] τ(E) 后处理失败：%s（transport.json 不受影响）" % res["error"], file=sys.stderr)
    out.write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    if rc == 0:
        i300 = int(np.argmin(np.abs(np.asarray(res["temperatures_K"]) - 300.0)))
        print("[DONE] τ 写出：%s（τ_σ/τ_eff 最大相对差 %.1f%%）" % (out, 100 * res["check"]["tau_sigma_vs_eff_max_rel"]))
        for n, car in enumerate(res["carrier"]):
            te = res["tau_eff_fs"]["overall"]["inplane"][n][i300]
            print("[..] %.0f K  %-8s doping=%.3e cm^-3  τ_eff,overall(面内)=%s fs"
                  % (res["temperatures_K"][i300], car, res["doping_cm3"][n], "—" if te is None else "%.2f" % te))
    return rc


if __name__ == "__main__":
    sys.exit(main())
