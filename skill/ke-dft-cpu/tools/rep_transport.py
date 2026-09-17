import json, sys
import numpy as np
from pathlib import Path

def load(p): return json.loads(Path(p).read_text())
def parse(d):
    dop = np.array(d["doping"], dtype=float)
    T = np.array(d["temperatures"], dtype=float)
    S = np.array(d["seebeck"]); C = np.array(d["conductivity"])
    return dop, T, S, C

def rows300(d):
    dop, T, S, C = parse(d)
    it = int(np.argmin(np.abs(T - 300)))
    out = []
    for i in range(len(dop)):
        sd = [S[i, it, k, k] for k in range(3)]
        cd = [C[i, it, k, k] for k in range(3)]
        s = sum(sd) / 3.0; c = sum(cd) / 3.0
        pf = sum(sd[k]**2 * cd[k] * 1e-9 for k in range(3)) / 3.0
        out.append((dop[i], s, c, pf))
    return out

mode = sys.argv[1]
if mode == "table":
    print("%-13s %10s %12s %10s" % ("doping", "S(uV/K)", "sigma(S/m)", "PF"))
    for d, s, c, pf in rows300(load(sys.argv[2])):
        print("%-13.5g %10.2f %12.4g %10.3f" % (d, s, c, pf))
elif mode == "cmp":
    ro = rows300(load(sys.argv[2])); rn = rows300(load(sys.argv[3]))
    print("%-13s %9s %9s %6s | %11s %11s %6s | %8s" % (
        "doping", "S_old", "S_new", "rS", "sig_old", "sig_new", "rC", "rPF"))
    for (d, s0, c0, p0), (_, s1, c1, p1) in zip(ro, rn):
        rS = s1 / s0 if abs(s0) > 1e-12 else float("nan")
        rC = c1 / c0 if abs(c0) > 1e-12 else float("nan")
        rP = p1 / p0 if abs(p0) > 1e-12 else float("nan")
        print("%-13.5g %9.2f %9.2f %6.3f | %11.4g %11.4g %6.3f | %8.3f" % (d, s0, s1, rS, c0, c1, rC, rP))
elif mode == "temp":
    d = load(sys.argv[2]); targets = [float(x) for x in sys.argv[3].split(",")]
    dop, T, S, C = parse(d)
    print("%-6s %s" % ("T(K)", " ".join("%9.4g" % t for t in targets)))
    for it in range(len(T)):
        vals = []
        for tv in targets:
            i = int(np.argmin(np.abs(dop - tv)))
            vals.append(float(np.trace(S[i, it]) / 3.0))
        print("%-6.0f %s" % (T[it], " ".join("%9.2f" % v for v in vals)))
elif mode == "pos":
    d = load(sys.argv[2]); targets = [float(x) for x in sys.argv[3].split(",")]
    dop, T, S, C = parse(d)
    print("%-13s %10s %10s %10s  %s" % ("doping", "S_min", "S_max", "T@S_min", "全程为正"))
    allpos = True
    for tv in targets:
        i = int(np.argmin(np.abs(dop - tv)))
        sv = np.array([float(np.trace(S[i, it]) / 3.0) for it in range(len(T))])
        j = int(np.argmin(sv))
        pos = bool((sv > 0).all())
        allpos = allpos and pos
        print("%-13.5g %10.2f %10.2f %10.0f  %s" % (dop[i], sv.min(), sv.max(), T[j], "是" if pos else "**否**"))
    print("结论: 全部浓度" + ("全程为正" if allpos else "**存在非正采样点**"))