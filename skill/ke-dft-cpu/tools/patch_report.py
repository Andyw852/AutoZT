import json, re, sys
from pathlib import Path

SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

def read_raw(p):
    with open(p, "r", encoding="utf-8", newline="") as f: return f.read()
def write_raw(p, s):
    with open(p, "w", encoding="utf-8", newline="") as f: f.write(s)
def fmt_sigma(v):
    if v == 0: return "0"
    e = 0; x = abs(v)
    while x >= 10: x /= 10.0; e += 1
    while x < 1: x *= 10.0; e -= 1
    return "%.2f×10%s" % (x, str(e).translate(SUP))

ANCHOR = "300 K平均σ/S·m⁻¹"      # 输运表表头锚点：改动限定在该表内

def patch(report, data, out=None):
    s = read_raw(report); n0 = len(s); log = []
    # ---- 1) 表 10 数据行（限定在输运表内）----
    mt = re.search(r"<table>(?:(?!</table>).)*?" + re.escape(ANCHOR) + r"(?:(?!</table>).)*?</table>", s, re.S)
    if not mt: raise SystemExit("[ERROR] 找不到输运表（锚点 %r）" % ANCHOR)
    tbl = mt.group(0); nrow = 0
    for tag, sigma, S, PF in data["rows"]:
        pat = re.compile(r"(<tr><td>" + tag + r"</td>(?:<td>[^<]*</td>){2}<td>)([^<]*)(</td><td>)([^<]*)(</td><td>)([^<]*)(</td></tr>)")
        m = pat.search(tbl)
        if not m:
            log.append((tag, "ROW-NOT-FOUND", "")); continue
        old = (m.group(2), m.group(4), m.group(6))
        new = (fmt_sigma(sigma), "%.2f" % S, "%.2f" % PF)
        if old == new:
            log.append((tag, old, "未变")); continue
        tbl = tbl[:m.start()] + m.group(1) + new[0] + m.group(3) + new[1] + m.group(5) + new[2] + m.group(7) + tbl[m.end():]
        nrow += 1; log.append((tag, old, new))
    s = s[:mt.start()] + tbl + s[mt.end():]
    log.append(("表内改动行数", "", nrow))
    # ---- 2) 其余文本组：任意 key（除 rows）都是 (old,new) 列表 ----
    for grp in [k for k in data if k != "rows"]:
        for pair in data[grp]:
            old_v, new_v = pair[0], pair[1]
            if old_v == new_v: continue
            c = s.count(old_v)
            if c != 1:
                log.append((grp, old_v, "[WARN] 匹配 %d 处，跳过" % c)); continue
            s = s.replace(old_v, new_v)
            log.append((grp, old_v, new_v))
    if out: write_raw(out, s)
    return log, len(s) - n0

if __name__ == "__main__":
    rep = sys.argv[1]; data = json.loads(read_raw(sys.argv[2]))
    out = sys.argv[3] if len(sys.argv) > 3 else None
    log, dl = patch(rep, data, out)
    for tag, old, new in log:
        print("  %-14s %-26s -> %s" % (tag, str(old), str(new)))
    print("  字符数变化:", dl)