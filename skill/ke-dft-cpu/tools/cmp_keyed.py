import csv, sys
def read(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
A = read(sys.argv[1]); B = read(sys.argv[2])
def key(r):
    return (r["material"], int(float(r["temperature_K"])), round(float(r["doping_cm3"]), 6), r["direction"])
ka = {}; kb = {}
for r in A: ka[key(r)] = r
for r in B: kb[key(r)] = r
sa = set(ka.keys()); sb = set(kb.keys())
common = sa.intersection(sb)
print("keys: archive %d, regen %d, common %d" % (len(sa), len(sb), len(common)))
print("archive-only:", len(sa - sb), " regen-only:", len(sb - sa))
for k in list(sa - sb)[:3]: print("   archive-only key:", k)
for k in list(sb - sa)[:3]: print("   regen-only key:", k)
cols = [c for c in A[0].keys() if c not in ("material", "temperature_K", "doping_cm3", "direction")]
bad = {}; chk = 0
for k in sorted(common):
    for c in cols:
        try:
            fa = float(ka[k].get(c)); fb = float(kb[k].get(c))
        except (TypeError, ValueError):
            continue
        chk += 1
        if abs(fa - fb) > 1e-9 * max(abs(fa), 1e-12):
            bad[c] = bad.get(c, 0) + 1
print("compared %d values; per-column diffs:" % chk)
if bad:
    for c, n in sorted(bad.items(), key=lambda x: -x[1]): print("   %-28s %d" % (c, n))
else:
    print("   ALL MATCH")
if common:
    k0 = sorted(common)[0]
    print(); print("sample key:", k0)
    for c in cols[:6]: print("   %-24s archive=%s  regen=%s" % (c, ka[k0].get(c), kb[k0].get(c)))