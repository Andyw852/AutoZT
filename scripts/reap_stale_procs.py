#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""reap_stale_procs.py —— 回收 agent 遗留的临时后台进程（孤儿 + /tmp 脚本 + 超龄）。

背景
----
agent 常用 `ssh host 'nohup python /tmp/xxx.py > /tmp/xxx.log 2>&1 &'` 启动临时诊断
进程；ssh 断开后进程变成 PPID=1 的孤儿。若忘记回收，它会长期空转占核
（2026-10-20 实测：/tmp/diag_imag.py 空转 17 小时、累计约 244 核·小时）。
本脚本只处理【本人、孤儿、命令行指向 /tmp 下 .py、且存活超过阈值】的进程；
其它一律不动。

安全规则（缺一不可）
--------------------
1. 进程属主 = 当前登录用户；
2. PPID == 1（已被收养的孤儿；有父进程的一律不动）；
3. 命令行同时含 "/tmp/"、".py"、"python"；
4. 存活时长 >= --min-age-hours（默认 12 小时）；
5. 命令行不含白名单关键词（pytest/unittest/pip/conda/pheasy/vasp 等）；
6. 杀前在远端用 ps 复核 PID 的命令行仍匹配，防止 PID 复用误杀。

用法
----
  python3 scripts/reap_stale_procs.py                      # 只报告（dry-run，默认）
  python3 scripts/reap_stale_procs.py --kill               # 真正回收
  python3 scripts/reap_stale_procs.py --hosts jzzn,wangchao_3090,A800
  python3 scripts/reap_stale_procs.py --min-age-hours 6 --json

由 monitor.sh 每 10 分钟自动调用（AUTOZT_REAP_STALE=0 可关闭）。
本机加 --local 一并检查（默认只查远端）。
"""
import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG = os.path.join(HERE, "tmp", "reap_stale_procs.log")

DEFAULT_HOSTS = ["jzzn", "wangchao_3090", "A800"]

WHITELIST = (
    "pytest", "unittest", "pip ", "pip3", "conda", "pheasy", "vasp",
    "mpirun", "mpi ", "rsync", "scp ", "squeue", "sbatch", "jupyter",
)

# 远端采集：先报告远端登录名（各集群用户名不同：jzzn=wangchao、3090=wangchaoyue852、
# A800=wangcch），再输出全量进程；分类时按【远端】用户名过滤。
REMOTE_PS = ('echo "ME=$(id -u)"; '
             'ps -eo pid=,ppid=,etimes=,uid=,args= 2>/dev/null || ps -eo pid,ppid,etime,uid,args')


def _ssh(host, command, timeout=60):
    return subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host, command],
        capture_output=True, text=True, timeout=timeout)


def parse_ps(text, uid):
    """按【数字 UID】过滤：ps 的 user 列默认截断成 8 字符（wangcha+），不能直接比用户名。"""
    procs = []
    for line in text.splitlines():
        parts = line.split(None, 4)
        if len(parts) < 5:
            continue
        pid, ppid, etimes, puid, args = parts
        try:
            pid_i, ppid_i, age_i, uid_i = int(pid), int(ppid), int(etimes), int(puid)
        except ValueError:
            continue
        if uid_i != int(uid):
            continue
        procs.append({"pid": pid_i, "ppid": ppid_i, "age": age_i, "args": args.strip()})
    return procs


def classify(procs, min_age_hours):
    out = []
    min_age = float(min_age_hours) * 3600.0
    for p in procs:
        a = p["args"].lower()
        if p["ppid"] != 1:
            continue
        if "/tmp/" not in p["args"] or ".py" not in a or "python" not in a:
            continue
        if p["age"] < min_age:
            continue
        # 白名单只匹配 /tmp/ 下的脚本名，不匹配解释器路径（miniconda3/... 会误命中 "conda"）
        script_part = p["args"].split("/tmp/", 1)[-1].lower()
        if any(w in script_part for w in WHITELIST):
            continue
        out.append(p)
    return out


def kill_one(host, cand):
    """远端复核命令行后 TERM→KILL。返回 'killed' / 'gone' / 'mismatch'。"""
    script = cand["args"].split("/tmp/")[-1]
    needle = "/tmp/" + script.split()[0].strip("'\"")
    shell = (
        'if ps -p %d >/dev/null 2>&1; then '
        'cur=$(ps -p %d -o args= 2>/dev/null); '
        'case "$cur" in *%s*) kill -TERM %d 2>/dev/null; sleep 2; '
        'kill -KILL %d 2>/dev/null; echo killed;; '
        '*) echo mismatch;; esac; else echo gone; fi'
        % (cand["pid"], cand["pid"], needle, cand["pid"], cand["pid"]))
    try:
        r = _ssh(host, shell, timeout=50)
        return (r.stdout or "").strip().splitlines()[-1] if r.stdout.strip() else "error"
    except Exception as e:  # noqa: BLE001
        return "error:%s" % e


def discover_hosts():
    hosts = []
    setting = os.path.join(HERE, "setting")
    try:
        for fn in sorted(os.listdir(setting)):
            if not fn.endswith(".yaml") or fn.startswith("template"):
                continue
            try:
                txt = open(os.path.join(setting, fn), encoding="utf-8").read()
            except OSError:
                continue
            m = re.search(r"^ssh_host:\s*(\S+)", txt, re.M)
            if m and m.group(1) and m.group(1) not in ('""', "''"):
                hosts.append(m.group(1).strip('"\''))
    except OSError:
        pass
    seen, uniq = set(), []
    for h in hosts:
        if h not in seen:
            seen.add(h)
            uniq.append(h)
    return uniq or list(DEFAULT_HOSTS)


def main():
    ap = argparse.ArgumentParser(description="回收遗留的 /tmp 临时 python 进程（孤儿+超龄）")
    ap.add_argument("--hosts", default=None, help="逗号分隔的 ssh 别名；默认读 setting/*.yaml")
    ap.add_argument("--min-age-hours", type=float, default=12.0, help="存活龄阈值（小时），默认 12")
    ap.add_argument("--kill", action="store_true", help="真正回收（默认只报告）")
    ap.add_argument("--local", action="store_true", help="同时检查本机")
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    ap.add_argument("--quiet", action="store_true", help="无候选时不输出")
    args = ap.parse_args()

    try:
        me = int(subprocess.run(["id", "-u"], capture_output=True, text=True).stdout.strip())
    except (OSError, ValueError):
        me = os.getuid()

    hosts = [h.strip() for h in (args.hosts.split(",") if args.hosts else discover_hosts()) if h.strip()]
    report = {"time": datetime.now().isoformat(timespec="seconds"), "min_age_hours": args.min_age_hours,
              "kill": args.kill, "hosts": {}, "total": 0}
    loglines = []

    for host in hosts:
        try:
            r = _ssh(host, REMOTE_PS)
            if r.returncode != 0:
                report["hosts"][host] = {"error": (r.stderr or "ssh/ps failed").strip()[-200:]}
                loglines.append("%s: ssh/ps failed: %s" % (host, (r.stderr or "").strip()[-160:]))
                continue
            remote_uid, body = me, r.stdout
            for ln in r.stdout.splitlines():
                if ln.startswith("ME="):
                    try:
                        remote_uid = int(ln[3:].strip())
                    except ValueError:
                        pass
                    body = r.stdout.split(ln, 1)[1]
                    break
            cands = classify(parse_ps(body, remote_uid), args.min_age_hours)
        except Exception as e:  # noqa: BLE001
            report["hosts"][host] = {"error": str(e)[:200]}
            loglines.append("%s: error %s" % (host, e))
            continue
        killed = []
        if args.kill:
            for c in cands:
                res = kill_one(host, c)
                killed.append({"pid": c["pid"], "age_h": round(c["age"] / 3600.0, 1),
                               "result": res, "args": c["args"][:120]})
        report["hosts"][host] = {"candidates": cands, "killed": killed}
        report["total"] += len(cands)

    if args.local:
        try:
            r = subprocess.run(["ps", "-eo", "pid=,ppid=,etimes=,uid=,args="],
                               capture_output=True, text=True)
            cands = classify(parse_ps(r.stdout, me), args.min_age_hours)
            report["hosts"]["@local"] = {"candidates": cands, "killed": []}
            report["total"] += len(cands)
        except Exception as e:  # noqa: BLE001
            report["hosts"]["@local"] = {"error": str(e)}

    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=1))
    elif not (args.quiet and report["total"] == 0):
        print("[reap] %s  min_age=%.1fh  kill=%s  候选=%d" %
              (report["time"], args.min_age_hours, args.kill, report["total"]))
        for host, d in report["hosts"].items():
            if "error" in d:
                print("  %-16s ERROR %s" % (host, d["error"]))
                continue
            if not d["candidates"] and not d["killed"]:
                print("  %-16s clean" % host)
                continue
            for c in d["candidates"]:
                print("  %-16s pid=%-8d age=%5.1fh %s" % (host, c["pid"], c["age"] / 3600.0, c["args"][:110]))
            for k in d.get("killed", []):
                print("  %-16s KILL pid=%-8d %s" % (host, k["pid"], k["result"]))

    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(report, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
