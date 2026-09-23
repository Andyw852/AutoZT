import json, sys, csv
import numpy as np
from pathlib import Path

def load(p): return json.loads(Path(p).read_text())
def pair(a, b):
    if abs(a) < 1e-30: return float("nan"), float("nan")
    return b - a, (b - a) / abs(a) * 100.0

MATS = ["Pb2Bi2Te5", "Pb2Sb2Te5", "Sn2Bi2Te5", "Sn2Sb2Te5"]
BASE = "/public/home/.../ke_work/%s/ke-dft-cpu/step8_amset/transport.json"
SOCP = "/public/home/.../ke_soc_20260910/%s/ke-dft-cpu/step8_amset/transport.json"
KEYS = ["mobility", "conductivity", "seebeck", "kappa_e", "PF"]
SRC = {"mobility": ("mobility", "overall"), "conductivity": ("conductivity",),
       "seebeck": ("seebeck",), "kappa_e": ("electronic_thermal_conductivity",)}
HEAD = ["material", "temperature_K", "doping_cm3", "direction"]
for k in KEYS: HEAD += [k+"_noSOC", k+"_SOC", k+"_difference", k+"_relative_pct"]

def get(d, key, i, it):
    if key == "PF": return None
    x = d
    for p in SRC[key]: x = x[p]
    return np.array(x)[i, it]

def pf_of(S, C, axis=None):
    if axis is None:
        return float(sum(S[j][j]**2 * C[j][j] * 1e-9 for j in range(3)) / 3.0)
    return float(S[axis][axis]**2 * C[axis][axis] * 1e-9)

rows = []
for m in MATS:
    bn, so = load(BASE % m), load(SOCP % m)
    dopn = np.array(bn["doping"], float); T = np.array(bn["temperatures"], float)
    Sn_all = np.array(bn["seebeck"]); Cn_all = np.array(bn["conductivity"])
    Ss_all = np.array(so["seebeck"]); Cs_all = np.array(so["conductivity"])
    for it in range(len(T)):
        for i in range(len(dopn)):
            for dname, di in [("x", 0), ("y", 1), ("z", 2), ("mean", None)]:
                r = {"material": m, "temperature_K": int(T[it]),
                     "doping_cm3": repr(float(dopn[i])), "direction": dname}
                for key in KEYS:
                    if key == "PF":
                        a = pf_of(Sn_all[i, it], Cn_all[i, it], di)
                        b = pf_of(Ss_all[i, it], Cs_all[i, it], di)
                    else:
                        An = get(bn, key, i, it); Bs = get(so, key, i, it)
                        if di is None:
                            a = float(np.trace(An) / 3.0); b = float(np.trace(Bs) / 3.0)
                        else:
                            a = float(An[di][di]); b = float(Bs[di][di])
                    d_, r_ = pair(a, b)
                    r[key+"_noSOC"] = repr(float(a)); r[key+"_SOC"] = repr(float(b))
                    r[key+"_difference"] = repr(float(d_)); r[key+"_relative_pct"] = repr(float(r_))
                rows.append(r)
with open(sys.argv[1], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=HEAD); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote", len(rows), "rows")