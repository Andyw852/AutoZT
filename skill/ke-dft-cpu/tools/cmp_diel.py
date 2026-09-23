import json, subprocess, sys
def get(d):
    out = subprocess.run(["python3", "/tmp/extract_diel.py", d], capture_output=True, text=True)
    return json.loads(out.stdout)
A = get("/public/home/.../ke_work/Pb2Sb2Te5/ke-dft-cpu/step5_dielect_soc_dfpt")
B = get("/public/home/.../ke_work/Pb2Sb2Te5/ke-dft-cpu/kconv_soc_dfpt/12x12x1")
print("%-14s %-26s %-26s %8s" % ("量", "生产 9x9x1 (NKPTS=%d)" % A["nkpts"], "新 12x12x1 (NKPTS=%d)" % B["nkpts"], "变化"))
print("-" * 84)
for label, ka, kb in [("eps_inf", "eps_inf_diag", "eps_inf_diag"), ("eps_ion", None, None), ("eps_static", "eps_static_diag", "eps_static_diag")]:
    if ka:
        a, b = A[ka], B[kb]
    else:
        a = [A["eps_static"][i][i] - A["eps_inf"][i][i] for i in range(3)]
        b = [B["eps_static"][i][i] - B["eps_inf"][i][i] for i in range(3)]
    for j, nm in enumerate(["xx", "yy", "zz"]):
        ch = (b[j] - a[j]) / a[j] * 100 if a[j] else float("nan")
        print("%-14s %-26.4f %-26.4f %+7.2f%%" % (label + "_" + nm, a[j], b[j], ch))
froh = lambda e_inf, e_0: 1.0 / e_inf - 1.0 / e_0
fa = froh(A["eps_inf_diag"][0], A["eps_static_diag"][0])
fb = froh(B["eps_inf_diag"][0], B["eps_static_diag"][0])
print()
print("Frohlich 耦合 (xx): 旧 %.6f -> 新 %.6f  (%+.1f%%)" % (fa, fb, (fb - fa) / fa * 100))