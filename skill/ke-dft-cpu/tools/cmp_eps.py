import re, sys, subprocess, json
from pathlib import Path
TAG = "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"
def grab(txt):
    i = txt.rfind(TAG)
    if i < 0: return None
    rows = []
    for ln in txt[i:].splitlines()[1:]:
        n = re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?", ln)
        if len(n) >= 3: rows.append([float(v) for v in n[:3]])
        if len(rows) == 3: break
    return rows if len(rows) == 3 else None
print("%-12s %-22s %-22s %s" % ("材料", "生产 9x9x1 eps_inf(xx,zz)", "新 9x9x2 eps_inf(xx,zz)", "变化"))
print("-" * 88)
for m in ["Pb2Sb2Te5", "Sn2Sb2Te5", "Pb2Bi2Te5", "Sn2Bi2Te5"]:
    old = Path("/public/home/<user>/ke_work/%s/ke-dft-cpu/step5_dielect_soc_dfpt/OUTCAR" % m)
    new = Path("/public/home/<user>/ke_work/%s/ke-dft-cpu/soc_dfpt_v2/OUTCAR" % m)
    if m == "Pb2Sb2Te5":
        new = Path("/public/home/<user>/ke_work/%s/ke-dft-cpu/kconv_soc_dfpt/9x9x2/OUTCAR" % m)
    a = grab(old.read_text(errors="ignore"))
    b = grab(new.read_text(errors="ignore")) if new.is_file() else None
    if not a or not b:
        print("%-12s %s" % (m, "缺数据")); continue
    o = (a[0][0], a[2][2]); n2 = (b[0][0], b[2][2])
    print("%-12s %-22s %-22s xx %+.1f%%  zz %+.1f%%" % (m, "%.3f / %.3f" % o, "%.3f / %.3f" % n2,
          (n2[0]-o[0])/o[0]*100, (n2[1]-o[1])/o[1]*100))