
import re, sys, json
from pathlib import Path

TAG_INF = "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"
TAG_ION = "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION"

def grab(txt, tag):
    i = txt.rfind(tag)
    if i < 0:
        return None
    rows = []
    for ln in txt[i:].splitlines()[1:]:
        n = re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|NaN|nan", ln)
        if len(n) >= 3:
            rows.append([float(v) for v in n[:3]])
        if len(rows) == 3:
            break
    return rows if len(rows) == 3 else None

d = Path(sys.argv[1])
txt = (d / "OUTCAR").read_text(errors="ignore")
ei, io = grab(txt, TAG_INF), grab(txt, TAG_ION)
if ei is None or io is None:
    print(json.dumps({"ok": False, "reason": "eps_inf or ionic block missing"}))
    sys.exit(1)
e0 = [[ei[r][c] + io[r][c] for c in range(3)] for r in range(3)]
m = re.search(r"NKPTS\s*=\s*(\d+)", txt)
print(json.dumps({
    "ok": True, "dir": str(d), "nkpts": int(m.group(1)) if m else None,
    "eps_inf": ei, "eps_ion": io, "eps_static": e0,
    "eps_inf_diag": [ei[i][i] for i in range(3)],
    "eps_static_diag": [e0[i][i] for i in range(3)],
}, indent=1))

