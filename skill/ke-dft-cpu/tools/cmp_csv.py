import csv, sys
from pathlib import Path

def read(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

A = read(sys.argv[1]); B = read(sys.argv[2])
print("行数: 归档 %d, 重生成 %d" % (len(A), len(B)))
keys = [k for k in A[0] if not k.startswith("\ufeff")]
bad = 0; checked = 0
for ra, rb in zip(A, B):
    for k in keys:
        va, vb = ra.get(k), rb.get(k)
        if va is None or vb is None: continue
        try:
            fa, fb = float(va), float(vb)
            rel = abs(fa - fb) / max(abs(fa), 1e-12)
            checked += 1
            if rel > 1e-9:
                bad += 1
                if bad <= 6: print("  差异 %s %s: 归档=%s 重生=%s (rel=%.2e)" % (ra["material"], k, va, vb, rel))
        except ValueError:
            if str(va) != str(vb): print("  文本差异 %s: %s vs %s" % (k, va, vb))
print("比对数值 %d 个，超出容差 %d 个" % (checked, bad))
print("结论:", "完全一致" if bad == 0 else "存在差异")