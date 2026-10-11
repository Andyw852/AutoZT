#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_tau_energy_e2e.py —— [0035] tau_energy.analyze 的真 AMSET 端到端自检（需 amset + pymatgen 环境）。

用法（在 skill/ke-dft-cpu/step8.4_amset2d 目录里，约 3–10 分钟）：
    python test_tau_energy_e2e.py [--nworkers 2]
退出码 0 = 全部 PASS。

模型与 test_dp_optical_e2e.py 相同（六方单原子 2D 晶格，导带/价带极值在 K、K′）：ADP + ODP + IVS，电子开 ADP 谷内掩码。
跑一次 AMSET（内存里），再用同一个 AmsetData 调 tau_energy.analyze，检查：
  1. 复现：analyze 用 AMSET 的 _calculate_mobility 重算的 overall 迁移率 = AMSET 写的 transport json（相对差 < 1e-9）；
  2. 常数率精确性：两机制各自常数 τ_A、τ_B -> τ_eff,A = τ_A、τ_eff,B = τ_B、overall = 1/(1/τ_A+1/τ_B)、
     τ(E) 每个箱都等于 τ_A、τ_σ = overall（相对差只剩 json 的 6/5 位有效数字取整）；
  3. τ_σ（逐 k 点电导加权）与 τ_eff,overall 在 300 K 相差 < 10%（两条独立路径）；
  4. 物理：电子（掩码开 = 只谷内）的 τ_eff,ADP 与 Bardeen–Shockley τ = ħ³C₂D/(k_BT m_d E1²) 相差 < 25%
     （m_d 由紧束缚带在 K 点的曲率算）；空穴（掩码关，K↔K′ 谷间也算进 ADP）τ_eff,ADP/τ_BS ≈ 0.5；
  5. ODP 单机制 τ(E)：电子在 ħω 以下只能吸收声子、以上还能发射 -> τ_ODP(E) 在 ħω 之上明显变短。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FAILS = []
HBAR = 1.054571817e-34
KB = 1.380649e-23
EV = 1.602176634e-19


def ok(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _tb(kf, scale):
    S = (np.cos(2 * np.pi * kf[..., 0]) + np.cos(2 * np.pi * kf[..., 1])
         + np.cos(2 * np.pi * (kf[..., 0] + kf[..., 1])))
    return scale * (S + 1.5)


def mass_at_K(lat, scale):
    """紧束缚带 scale·(S+1.5) 在 K 的有效质量（kg，六方 K 点各向同性），有限差分。"""
    rec = np.asarray(lat.reciprocal_lattice.matrix)          # 行 = b_i（含 2π），Å⁻¹
    K = np.array([1 / 3, 1 / 3, 0.0])
    h = 1e-3
    def E(q):
        return float(_tb(K + np.linalg.solve(rec.T, q), scale))
    d2 = (E(np.array([h, 0, 0])) + E(np.array([-h, 0, 0])) - 2 * E(np.zeros(3))) / h ** 2   # eV·Å²
    return HBAR ** 2 / (abs(d2) * EV * 1e-20)


def child(run, nworkers):
    os.makedirs(run, exist_ok=True)
    for f in ("amset2d_plugin.py", "dp_optical.py", "amset_ir_fix.py", "amset_fermi_fix.py"):
        shutil.copy(HERE / f, run)
    os.chdir(run)
    sys.path.insert(0, run)
    sys.path.insert(1, str(HERE))
    import test_dp_optical_e2e as E2E
    Path("2d_correction.json").write_text(json.dumps(E2E.record(True)))
    import amset2d_plugin  # noqa: F401
    import tau_energy as TE
    from amset.core.run import Runner
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.bandstructure import BandStructure
    from pymatgen.electronic_structure.core import Spin

    lat = Lattice.hexagonal(3.18, 20.0)
    st = Structure(lat, ["Mo"], [[0, 0, 0]])
    n = 30
    k1, k2, k3 = np.meshgrid(np.arange(n) / n, np.arange(n) / n, [0.0, 0.5], indexing="ij")
    kf = np.stack([k1.ravel(), k2.ravel(), k3.ravel()], 1)
    cb = 1.8 + _tb(kf, 1.0)
    vb = -_tb(kf, 0.8)
    bs = BandStructure(kf, {Spin.up: np.array([-6.0 + 0.0 * cb, vb, cb, cb + 4.0])}, lat.reciprocal_lattice,
                       efermi=0.9, structure=st)
    Cv = np.zeros((6, 6))
    Cv[0, 0] = Cv[1, 1] = 60
    Cv[0, 1] = Cv[1, 0] = 15
    Cv[5, 5] = 22.5
    Cv[2, 2] = 1.0
    Cv[3, 3] = Cv[4, 4] = 0.5
    settings = {"scattering_type": ["ADP", "ODP", "IVS"], "doping": [-1e17, 1e17], "temperatures": [200, 300],
                "deformation_potential": (2.0, 4.0), "elastic_constant": Cv.tolist(),
                "unity_overlap": True, "write_mesh": False, "nworkers": nworkers, "interpolation_factor": 4,
                "free_carrier_screening": False, "file_format": "json", "print_log": False,
                "energy_cutoff": 1.5}
    ad = Runner(bs, 4, settings).run(directory=".")
    ad = ad[0] if isinstance(ad, tuple) else ad
    tj = json.loads(Path(sorted(Path(".").glob("transport_*.json"))[0]).read_text())
    rates = {s: np.array(ad.scattering_rates[s]) for s in ad.scattering_rates}
    labels = list(ad.scattering_labels)
    res = TE.analyze(ad, rates, labels, mobility_ref=tj["mobility"]["overall"])
    tA, tB = 2.0e-14, 5.0e-14
    rc = {s: np.stack([np.full(rates[s].shape[1:], 1 / tA), np.full(rates[s].shape[1:], 1 / tB)]) for s in rates}
    res_c = TE.analyze(ad, rc, ["A", "B"])
    out = {"res": res, "res_c": res_c, "tA_fs": tA * 1e15, "tB_fs": tB * 1e15,
           "m_e": mass_at_K(lat, 1.0), "m_h": mass_at_K(lat, 0.8),
           "C2D": 60e9 * 20.0e-10, "E1_e": 4.0, "E1_h": 2.0}
    Path("tau_check.json").write_text(json.dumps(out))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nworkers", type=int, default=2)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--_child", default=None, help=argparse.SUPPRESS)
    a = ap.parse_args()
    if a._child:
        return child(a._child, a.nworkers)
    tmp = tempfile.mkdtemp(prefix="tau_e2e_")
    run = os.path.join(tmp, "run")
    log = os.path.join(tmp, "run.log")
    print("[..] AMSET -> %s" % run, flush=True)
    with open(log, "w") as fh:
        rc = subprocess.call([sys.executable, "-W", "ignore", __file__, "--nworkers", str(a.nworkers), "--_child", run],
                             stdout=fh, stderr=subprocess.STDOUT)
    ok(rc == 0, "AMSET + analyze 退出码 %d（日志 %s）" % (rc, log))
    if rc != 0:
        print(Path(log).read_text(errors="ignore")[-3000:])
        print("\nFAIL（运行目录保留在 %s）" % tmp)
        return 1
    d = json.loads(Path(run, "tau_check.json").read_text())
    res, rc_ = d["res"], d["res_c"]
    # 1) 复现
    ok(res["check"]["mobility_reproduce_max_rel"] < 1e-9,
       "复现：重算 overall μ 与 AMSET transport json 最大相对差 %.1e" % res["check"]["mobility_reproduce_max_rel"])
    # 2) 常数率
    tA, tB = d["tA_fs"], d["tB_fs"]
    tO = 1 / (1 / tA + 1 / tB)

    def arr(x):
        return np.array([[np.nan if v is None else v for v in row] for row in x], float)
    worst = 0.0
    for key, want in (("A", tA), ("B", tB), ("overall", tO)):
        for dd in ("xx", "yy", "inplane"):
            worst = max(worst, float(np.nanmax(np.abs(arr(rc_["tau_eff_fs"][key][dd]) / want - 1))))
    ok(worst < 2e-6, "常数率：τ_eff,A/B/overall 精确等于 %g/%g/%.4g fs（最大相对差 %.1e；json 存 6 位有效数字）"
       % (tA, tB, tO, worst))
    tEA = np.array([[[np.nan if v is None else v for v in b] for b in t] for t in rc_["tau_E_fs"]["A"]], float)
    ok(np.nanmax(np.abs(tEA / tA - 1)) < 2e-5 and np.isfinite(tEA).sum() > 10,
       "常数率：τ_A(E) 每个非空箱都等于 τ_A（%d 个箱）" % np.isfinite(tEA).sum())
    ok(np.nanmax(np.abs(arr(rc_["tau_sigma_fs"]) / tO - 1)) < 2e-6, "常数率：τ_σ = overall")
    # 3) τ_σ vs τ_eff（真实率）
    temps = res["temperatures_K"]
    i300 = temps.index(300.0) if 300.0 in temps else int(np.argmin(np.abs(np.array(temps) - 300)))
    rel = res["check"]["tau_sigma_vs_eff_max_rel_300K"]
    ok(rel < 0.10, "真实率：300 K τ_σ 与 τ_eff,overall 最大相对差 %.1f%%（< 10%%）" % (100 * rel))
    # 4) Bardeen–Shockley
    for car, m, E1, want, tol in (("electron", d["m_e"], d["E1_e"], 1.0, 0.25), ("hole", d["m_h"], d["E1_h"], 0.5, 0.25)):
        n = res["carrier"].index(car)
        T = temps[i300]
        t_bs = HBAR ** 3 * d["C2D"] / (KB * T * m * (E1 * EV) ** 2) * 1e15
        t_adp = res["tau_eff_fs"]["ADP"]["inplane"][n][i300]
        r = t_adp / t_bs
        ok(abs(r - want) < tol * want, "%s：τ_eff,ADP = %.2f fs，Bardeen–Shockley τ = %.2f fs（m_d = %.3f m0），比值 %.3f（应 ≈%.1f%s）"
           % (car, t_adp, t_bs, m / 9.1093837e-31, r, want, "，掩码开=只谷内" if car == "electron" else "，K↔K′ 也在 ADP 里"))
    # 5) ODP 发射阈值
    n = res["carrier"].index("electron")
    x = np.array(res["energy_bins_eV"]["centers"])
    t = np.array([np.nan if v is None else v for v in res["tau_E_fs"]["ODP"][n][i300]], float)
    lo = np.nanmedian(t[(x > 0.005) & (x < 0.04)])
    hi = np.nanmedian(t[(x > 0.07) & (x < 0.12)])
    ok(hi < 0.7 * lo, "电子 τ_ODP(E)：ħω=50 meV 以下 %.1f fs、以上 %.1f fs（发射通道打开 -> 变短）" % (lo, hi))
    if not a.keep and not FAILS:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\n%s" % ("全部 PASS" if not FAILS else "FAIL %d 项（运行目录保留在 %s）" % (len(FAILS), tmp)))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
