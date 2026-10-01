#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fakeslurm -- just enough SLURM for AutoZT on a machine without SLURM.

AutoZT submits with ``sbatch``, watches with ``squeue``/``sacct`` and stops
with ``scancel``.  On a workstation (setting/local.yaml, ssh_host empty) these
four commands are provided by this file through the thin wrappers next to it;
AutoZT puts this directory on PATH automatically when no real ``sbatch`` is
found (autozt/collect.py::_effective_remote_path_prefix).

    sbatch [opts] script.sh   -> "Submitted batch job N"; runs in the background
    squeue -u USER -h -o FMT  -> pending/running jobs; %i %j %T %M %D %R %Z %u %C
    sacct -j N -n -P -X -o JobID,State
    scancel N [N ...]

Scheduling: first come first served on a CPU budget -- FAKESLURM_MAX_CPUS
(default: every core of the machine).  A job asks for #SBATCH
--cpus-per-task (or -c) x --ntasks; requests larger than the budget are
clipped to it so they never wait forever.  The script runs in its submit
directory with SLURM_JOB_ID / SLURM_SUBMIT_DIR / SLURM_CPUS_PER_TASK /
SLURM_JOB_NAME set, stdout/stderr to #SBATCH --output/--error (default
slurm-%j.out).  State lives in $FAKESLURM_HOME (default ~/.fakeslurm/).
"""
import fcntl
import getpass
import json
import os
import re
import shlex
import signal
import subprocess
import sys
import time

HOME = os.path.expanduser(os.environ.get("FAKESLURM_HOME") or "~/.fakeslurm")
JOBS = os.path.join(HOME, "jobs")
LOCK = os.path.join(HOME, "lock")
POLL_S = float(os.environ.get("FAKESLURM_POLL") or 3.0)
ACTIVE = ("PENDING", "RUNNING")
SHORT = {"PENDING": "PD", "RUNNING": "R", "COMPLETED": "CD", "FAILED": "F",
         "CANCELLED": "CA"}


def _max_cpus():
    v = os.environ.get("FAKESLURM_MAX_CPUS")
    try:
        n = int(v) if v else 0
    except ValueError:
        n = 0
    return n if n > 0 else (os.cpu_count() or 1)


class _Lock(object):
    def __enter__(self):
        os.makedirs(JOBS, exist_ok=True)
        self.fh = open(LOCK, "a+")
        fcntl.flock(self.fh, fcntl.LOCK_EX)
        return self

    def __exit__(self, *a):
        fcntl.flock(self.fh, fcntl.LOCK_UN)
        self.fh.close()


def _path(jid):
    return os.path.join(JOBS, "%s.json" % jid)


def _load(jid):
    try:
        with open(_path(jid), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _save(job):
    tmp = _path(job["id"]) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(job, fh)
    os.replace(tmp, _path(job["id"]))


def _all():
    out = []
    try:
        names = os.listdir(JOBS)
    except OSError:
        return out
    for n in names:
        if n.endswith(".json"):
            j = _load(n[:-5])
            if j:
                out.append(j)
    out.sort(key=lambda j: int(j["id"]))
    return out


def _alive(pid):
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def _reap(jobs):
    """A RUNNING job whose runner died without writing a state is FAILED."""
    for j in jobs:
        if j["state"] not in ACTIVE:
            continue
        if not j.get("runner_pid") and time.time() - j["submit"] < 30:
            continue                            # sbatch still recording the runner pid
        if not _alive(j.get("runner_pid")):
            j["state"] = "FAILED"
            j["end"] = time.time()
            _save(j)
    return jobs


# ------------------------------------------------------------------ sbatch
def _directives(script):
    """#SBATCH options of a script -> dict (long names, = form normalised)."""
    d = {}
    alias = {"-c": "cpus-per-task", "-n": "ntasks", "-J": "job-name",
             "-o": "output", "-e": "error", "-N": "nodes", "-t": "time",
             "-p": "partition", "-q": "qos"}
    with open(script, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if not s.startswith("#SBATCH"):
                if s and not s.startswith("#"):
                    break                      # directives end at the first command
                continue
            for tok in shlex.split(s[len("#SBATCH"):], comments=True):
                if tok.startswith("--"):
                    k, _, v = tok[2:].partition("=")
                    d[k] = v
                elif tok in alias:
                    d["_pending"] = alias[tok]
                elif "_pending" in d:
                    d[d.pop("_pending")] = tok
    d.pop("_pending", None)
    return d


def _int(v, default):
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return default


def cmd_sbatch(argv):
    script, cli = None, {}
    it = iter(argv)
    for a in it:
        if a.startswith("--") and "=" in a:
            k, _, v = a[2:].partition("=")
            cli[k] = v
        elif a in ("-J", "--job-name"):
            cli["job-name"] = next(it, "")
        elif a in ("-c", "--cpus-per-task"):
            cli["cpus-per-task"] = next(it, "")
        elif a.startswith("-"):
            continue                               # ignored (partition, qos, ...)
        else:
            script = a
            break
    if not script or not os.path.isfile(script):
        sys.stderr.write("sbatch: error: Unable to open file %s\n" % script)
        return 1
    d = _directives(script)
    d.update(cli)
    cpus = max(1, _int(d.get("cpus-per-task"), 1)) * max(1, _int(d.get("ntasks"), 1))
    cpus = min(cpus, _max_cpus())
    wd = os.path.abspath(os.path.dirname(script) or ".")
    if os.path.abspath(os.getcwd()) != wd:
        wd = os.path.abspath(os.getcwd())     # SLURM runs in the submit directory
    with _Lock():
        ids = [int(j["id"]) for j in _all()]
        jid = str(max(ids) + 1 if ids else 1000)
        job = {"id": jid, "name": d.get("job-name") or os.path.basename(script),
               "script": os.path.abspath(script), "workdir": wd, "cpus": cpus,
               "output": d.get("output") or "slurm-%j.out",
               "error": d.get("error") or "", "user": getpass.getuser(),
               "state": "PENDING", "submit": time.time(), "start": None,
               "end": None, "runner_pid": None, "pgid": None, "rc": None}
        _save(job)
    runner = subprocess.Popen(
        [sys.executable, os.path.abspath(__file__), "_run", jid],
        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, start_new_session=True, close_fds=True)
    with _Lock():
        job = _load(jid)
        job["runner_pid"] = runner.pid
        _save(job)
    print("Submitted batch job %s" % jid)
    return 0


def _expand(pat, job):
    return pat.replace("%j", job["id"]).replace("%x", job["name"])


def cmd_run(jid):
    """Background runner: wait for CPUs (FIFO), run the script, record the state."""
    while True:
        with _Lock():
            jobs = _reap(_all())
            me = next((j for j in jobs if j["id"] == jid), None)
            if me is None or me["state"] != "PENDING":
                return 0                        # cancelled while waiting
            used = sum(j["cpus"] for j in jobs if j["state"] == "RUNNING")
            older = [j for j in jobs if j["state"] == "PENDING"
                     and int(j["id"]) < int(jid)]
            if not older and used + me["cpus"] <= _max_cpus():
                me["state"] = "RUNNING"
                me["start"] = time.time()
                _save(me)
                break
        time.sleep(POLL_S)
    out = os.path.join(me["workdir"], _expand(me["output"], me))
    err = os.path.join(me["workdir"], _expand(me["error"], me)) if me["error"] else out
    env = dict(os.environ)
    env.update({"SLURM_JOB_ID": jid, "SLURM_JOBID": jid,
                "SLURM_SUBMIT_DIR": me["workdir"], "SLURM_JOB_NAME": me["name"],
                "SLURM_CPUS_PER_TASK": str(me["cpus"]), "SLURM_NTASKS": "1",
                "SLURM_JOB_NUM_NODES": "1", "SLURM_NODELIST": "localhost"})
    with open(out, "a") as fo:
        fe = open(err, "a") if err != out else fo
        try:
            p = subprocess.Popen(["bash", me["script"]], cwd=me["workdir"], env=env,
                                 stdout=fo, stderr=fe, stdin=subprocess.DEVNULL,
                                 start_new_session=True)
            with _Lock():
                j = _load(jid)
                j["pgid"] = p.pid
                _save(j)
            rc = p.wait()
        finally:
            if fe is not fo:
                fe.close()
    with _Lock():
        j = _load(jid)
        if j["state"] == "RUNNING":
            j["state"] = "COMPLETED" if rc == 0 else "FAILED"
        j["rc"] = rc
        j["end"] = time.time()
        _save(j)
    return 0


# ------------------------------------------------------------------ squeue
def _elapsed(j):
    if not j.get("start"):
        return "0:00"
    s = int(time.time() - j["start"])
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    return "%d:%02d:%02d" % (h, m, s) if h else "%d:%02d" % (m, s)


def _field(j, code):
    return {"i": j["id"], "A": j["id"], "j": j["name"], "T": j["state"],
            "t": SHORT.get(j["state"], j["state"]), "M": _elapsed(j), "D": "1",
            "R": "localhost" if j["state"] == "RUNNING" else "(Resources)",
            "Z": j["workdir"], "u": j["user"], "C": str(j["cpus"]),
            "P": "local"}.get(code, "")


def cmd_squeue(argv):
    fmt, user, header = "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R", None, True
    it = iter(argv)
    for a in it:
        if a in ("-o", "--format"):
            fmt = next(it, fmt)
        elif a.startswith("--format="):
            fmt = a.split("=", 1)[1]
        elif a in ("-u", "--user"):
            user = next(it, None)
        elif a in ("-h", "--noheader"):
            header = False
    with _Lock():
        jobs = [j for j in _reap(_all()) if j["state"] in ACTIVE
                and (not user or j["user"] == user)]
    pat = re.compile(r"%\.?\d*([A-Za-z])")
    if header:
        print(pat.sub(lambda m: {"i": "JOBID", "j": "NAME", "T": "STATE", "t": "ST",
                                 "M": "TIME", "D": "NODES", "R": "NODELIST(REASON)",
                                 "Z": "WORK_DIR", "u": "USER", "C": "CPUS",
                                 "P": "PARTITION"}.get(m.group(1), ""), fmt))
    for j in jobs:
        print(pat.sub(lambda m: _field(j, m.group(1)), fmt))
    return 0


# ------------------------------------------------------------------ sacct
def cmd_sacct(argv):
    ids, it = [], iter(argv)
    for a in it:
        if a in ("-j", "--jobs"):
            ids += [x for x in next(it, "").split(",") if x]
        elif a.startswith("--jobs="):
            ids += [x for x in a.split("=", 1)[1].split(",") if x]
    with _Lock():
        jobs = {j["id"]: j for j in _reap(_all())}
    for jid in ids:
        j = jobs.get(jid)
        if j:
            print("%s|%s" % (jid, j["state"]))
    return 0


# ------------------------------------------------------------------ scancel
def cmd_scancel(argv):
    rc = 0
    for jid in [a for a in argv if not a.startswith("-")]:
        with _Lock():
            j = _load(jid)
            if not j:
                sys.stderr.write("scancel: error: Invalid job id %s\n" % jid)
                rc = 1
                continue
            if j["state"] not in ACTIVE:
                continue
            pg = j.get("pgid")
            j["state"] = "CANCELLED"
            j["end"] = time.time()
            _save(j)
        if pg:
            try:
                os.killpg(int(pg), signal.SIGTERM)
            except OSError:
                pass
    return rc


def main(argv):
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        return 2
    cmd, rest = argv[1], argv[2:]
    if cmd == "_run":
        return cmd_run(rest[0])
    return {"sbatch": cmd_sbatch, "squeue": cmd_squeue, "sacct": cmd_sacct,
            "scancel": cmd_scancel}.get(cmd, lambda a: 2)(rest)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
