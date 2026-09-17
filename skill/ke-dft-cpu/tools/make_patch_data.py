import json, sys, re
import numpy as np
from pathlib import Path

MINUS = "\u2212"          # 报告用 U+2212 减号
def fmtS(v):
    """报告风格：正数带 +，负数用 U+2212，两位小数"""
    s = "%.2f" % abs(v)
    return ("+" + s) if v >= 0 else (MINUS + s)

def load(p): return json.loads(Path(p).read_text())

def parse(d):
    return (np.array(d["doping"], float), np.array(d["temperatures"], float),
            np.array(d["seebeck"]), np.array(d["conductivity"]))

def rows300(d):
    dop, T, S, C = parse(d)
    it = int(np.argmin(np.abs(T - 300)))
    out = []
    for i in range(len(dop)):
        sd = [S[i, it, k, k] for k in range(3)]; cd = [C[i, it, k, k] for k in range(3)]
        out.append((dop[i], sum(sd)/3.0, sum(cd)/3.0,
                    sum(sd[k]**2 * cd[k] * 1e-9 for k in range(3))/3.0))
    return out

def at_T(d, conc, Tt, axis=None):
    dop, T, S, C = parse(d)
    i = int(np.argmin(np.abs(dop - conc)))
    it = int(np.argmin(np.abs(T - Tt)))
    if axis is None:
        return float(np.trace(S[i, it]) / 3.0)
    return float(S[i, it, axis, axis])

OLD, NEW, CONCS = sys.argv[1], sys.argv[2], [float(x) for x in sys.argv[3].split(",")]
TAGS = ["v0", "v1", "v2", "v3"]
o, n = load(OLD), load(NEW)
ro, rn = rows300(o), rows300(n)

rows = []
for tag, c in zip(TAGS, CONCS):
    i = int(np.argmin(np.abs(np.array(o["doping"], float) - c)))
    j = int(np.argmin(np.abs(np.array(n["doping"], float) - c)))
    rows.append([tag, rn[j][2], rn[j][1], rn[j][3]])

s300 = []
for c in CONCS:
    old_v = fmtS(at_T(o, c, 300)); new_v = fmtS(at_T(n, c, 300))
    s300.append([old_v, new_v])

v3 = CONCS[3]
temp_v3 = []
for Tt in (750.0, 800.0, 900.0):
    temp_v3.append([fmtS(at_T(o, v3, Tt)), fmtS(at_T(n, v3, Tt))])

dir800 = []
for ax in (0, 1, 2):
    ov = fmtS(at_T(o, v3, 800.0, ax)); nv = fmtS(at_T(n, v3, 800.0, ax))
    if ax == 0:
        dir800.append([ov, nv])          # Sx=Sy 共用一处，只替换一次
    elif ax == 2:
        dir800.append([ov, nv])

print(json.dumps({"rows": rows, "s300": s300, "temp_v3": temp_v3, "dir800": dir800},
                 ensure_ascii=False, indent=1))