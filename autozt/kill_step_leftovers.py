import os, signal, sys, time
d = os.path.realpath(os.path.expanduser(sys.argv[1]))
dry = len(sys.argv) > 2 and sys.argv[2] == "--dry-run"
if not os.path.isdir(d):                             # a wrong/local path must not read
    print("NO_DIR %s" % d)                           # as "nothing left, alive=0"
    sys.exit(2)
me, my_pgid, uid = {os.getpid(), os.getppid()}, os.getpgrp(), os.getuid()

def info(p):
    try:
        st = open("/proc/%d/stat" % p).read()
        f = st[st.rindex(")") + 2:].split()          # comm may contain spaces
        argv = [a.decode("utf-8", "replace")
                for a in open("/proc/%d/cmdline" % p, "rb").read().split(b"\0") if a]
        return dict(ppid=int(f[1]), pgid=int(f[2]), argv=argv,
                    wd=os.readlink("/proc/%d/cwd" % p), uid=os.stat("/proc/%d" % p).st_uid)
    except (OSError, ValueError):
        return None                                  # gone, or not ours to read

def in_dir(wd):
    return wd == d or wd.startswith(d + "/")

def is_pheasy(argv):                                 # the executable, not any mention
    return any(os.path.basename(a) in ("pheasy", "pheasy-gpu", "run_pheasy.py")
               for a in argv[:3])

procs = {int(e): r for e in os.listdir("/proc") if e.isdigit() for r in [info(int(e))] if r}
mine = {p: r for p, r in procs.items() if r["uid"] == uid and p not in me}
targets = [p for p, r in mine.items() if is_pheasy(r["argv"]) and in_dir(r["wd"])]
groups, pids = set(), set()
for p in targets:
    g = procs[p]["pgid"]
    members = [q for q, r in procs.items() if r["pgid"] == g]
    if g > 1 and g != my_pgid and all(q in mine and in_dir(procs[q]["wd"]) for q in members):
        groups.add(g)                                # whole group lives in this step dir
    else:
        pids.add(p)                                  # shared group: only the process and
        pp = procs[p]["ppid"]                        # its driver script in this dir
        if pp in mine and in_dir(procs[pp]["wd"]) and not is_pheasy(procs[pp]["argv"]):
            pids.add(pp)
kids = {}
for q, r in mine.items():
    kids.setdefault(r["ppid"], []).append(q)
stack = list(pids)                                   # a killed process's own children
while stack:                                         # (workers, helpers) go with it
    for c in kids.get(stack.pop(), []):
        if c not in pids:
            pids.add(c); stack.append(c)
print("TO_KILL dir=%s pgids=%s pids=%s" % (d, sorted(groups), sorted(pids)))
for p in sorted(pids | {q for q, r in procs.items() if r["pgid"] in groups}):
    print("  %d %s" % (p, " ".join(procs[p]["argv"])[:160]))
if dry:
    sys.exit(0)

def running(p):                                     # zombies are dead, just unreaped
    try:
        st = open("/proc/%d/stat" % p).read()
        return st[st.rindex(")") + 2] != "Z"
    except OSError:
        return False

def send(sig):
    for g in groups:
        try: os.killpg(g, sig)
        except ProcessLookupError: pass
    for p in pids:
        try: os.kill(p, sig)
        except ProcessLookupError: pass

watch = set(pids) | {q for q, r in procs.items() if r["pgid"] in groups}
send(signal.SIGTERM)
for _ in range(20):                                  # up to 10 s to exit cleanly
    if not any(running(p) for p in watch): break
    time.sleep(0.5)
send(signal.SIGKILL)
time.sleep(0.5)
alive = [p for p in watch if running(p)]
print("LEFTOVER killed=%d alive=%d" % (len(watch), len(alive)))
