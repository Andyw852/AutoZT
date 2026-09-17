import sys, yaml
for f in sys.argv[1:]:
    d = yaml.safe_load(open(f))
    hi = d["high_frequency_dielectric"]
    st = d["static_dielectric"]
    ok_shape = (len(hi) == 3 and all(len(r) == 3 for r in hi))
    print("%-26s eps_inf 对角 %s  形状OK=%s  pop=%s" % (
        f.split("/")[-1], [hi[i][i] for i in range(3)], ok_shape, d.get("pop_frequency")))
    print("    eps_static 对角", [st[i][i] for i in range(3)])