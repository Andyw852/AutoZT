import re, sys, json
from pathlib import Path
NL = chr(10)

def compact(M):
    rows = ", ".join("[%s]" % ", ".join(repr(float(M[r][c])) for c in range(3)) for r in range(3))
    return "[%s]" % rows

def block(key, M):
    out = [key + ":"]
    for r in range(3):
        for c in range(3):
            out.append(("- - " if c == 0 else "  - ") + repr(float(M[r][c])))
    return out

def replace_key(text, key, M):
    lines = text.splitlines()
    out, i, hit = [], 0, False
    pat = re.compile(r"^" + re.escape(key) + r":(.*)$")
    while i < len(lines):
        m = pat.match(lines[i])
        if not m:
            out.append(lines[i]); i += 1; continue
        hit = True
        rest = m.group(1).strip()
        if rest:
            out.append(key + ": " + compact(M))
            i += 1
        else:
            out.extend(block(key, M))
            i += 1
            while i < len(lines) and (lines[i].lstrip().startswith("-") or lines[i].startswith(" ")):
                i += 1
    if not hit:
        raise SystemExit("key not found: " + key)
    return NL.join(out) + NL

settings = Path(sys.argv[1])
data = json.loads(Path(sys.argv[2]).read_text())
t = settings.read_text(encoding="utf-8")
t = replace_key(t, "high_frequency_dielectric", data["eps_inf"])
t = replace_key(t, "static_dielectric", data["eps_static"])
if len(sys.argv) > 3 and sys.argv[3] != "-":
    t = re.sub(r"^pop_frequency:.*$", "pop_frequency: " + str(sys.argv[3]), t, count=1, flags=re.M)
settings.write_text(t, encoding="utf-8", newline=NL)
print("updated:", settings)