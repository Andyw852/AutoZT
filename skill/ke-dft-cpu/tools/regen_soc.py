import json, sys, csv
import numpy as np
from pathlib import Path

def load(p): return json.loads(Path(p).read_text())

def at(d, T=300):
    dop = np.array(d["doping"], dtype=float)
    temps = np.array(d["temperatures"], dtype=float)
    it = int(np.argmin(np.abs(temps - T)))
    S = np.array(d["seebeck"])[:, it]
    C = np.array(d["conductivity"])[:, it]
    K = np.array(d["electronic_thermal_conductivity"])[:, it]
    M = np.array(d["mobility"]["overall"])[:, it]
    s = np.trace(S, axis1=1, axis2=2) / 3.0
    c = np.trace(C, axis1=1, axis2=2) / 3.0
    k = np.trace(K, axis1=1, axis2=2) / 3.0
    mu = np.trace(M, axis1=1, axis2=2) / 3.0
    pf = np.array([sum(S[i][j][j]**2 * C[i][j][j] * 1e-9 for j in range(3)) / 3.0 for i in range(len(dop))])
    return dop, s, c, k, mu, pf

def pair(a, b):
    """相对百分数按 |a| 归一（保留差值符号）—— 与归档 CSV 口径一致。"""
    if abs(a) < 1e-30: return float("nan"), float("nan")
    return b - a, (b - a) / abs(a) * 100.0

MATS = ["Pb2Bi2Te5", "Pb2Sb2Te5", "Sn2Bi2Te5", "Sn2Sb2Te5"]
BASE = "/public/home/<user>/ke_work/%s/ke-dft-cpu/step8_amset/transport.json"
SOC  = "/public/home/<user>/ke_soc_20260910/%s/ke-dft-cpu/step8_amset/transport.json"
HEAD = ["material","temperature_K","doping_cm3","direction",
        "mobility_noSOC","mobility_SOC","mobility_difference","mobility_relative_pct",
        "conductivity_noSOC","conductivity_SOC","conductivity_difference","conductivity_relative_pct",
        "seebeck_noSOC","seebeck_SOC","seebeck_difference","seebeck_relative_pct",
        "kappa_e_noSOC","kappa_e_SOC","kappa_e_difference","kappa_e_relative_pct",
        "PF_noSOC","PF_SOC","PF_difference","PF_relative_pct"]
rows = []
for m in MATS:
    dn, sn, cn, kn, mun, pn = at(load(BASE % m), 300)
    ds, ss, cs, ks, mus, ps = at(load(SOC % m), 300)
    for target in (-1e19, 1e19):
        i = int(np.argmin(np.abs(dn - target)))
        j = int(np.argmin(np.abs(ds - target)))
        r = {"material": m, "temperature_K": 300, "doping_cm3": repr(float(dn[i])), "direction": "mean"}
        for key, a, b in (("mobility", mun[i], mus[j]), ("conductivity", cn[i], cs[j]),
                          ("seebeck", sn[i], ss[j]), ("kappa_e", kn[i], ks[j]), ("PF", pn[i], ps[j])):
            d_, r_ = pair(float(a), float(b))
            r[key + "_noSOC"] = repr(float(a)); r[key + "_SOC"] = repr(float(b))
            r[key + "_difference"] = repr(float(d_)); r[key + "_relative_pct"] = repr(float(r_))
        rows.append(r)
with open(sys.argv[1], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=HEAD); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote", sys.argv[1], len(rows), "rows")