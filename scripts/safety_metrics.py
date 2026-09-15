#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""safety_metrics.py —— PhonoAgent 的三项安全/复现指标（论文里可引用的数字）。

指标定义
  1. interception_rate  破坏性/变更动作在 agent 会话下的**拦截率**（应 100%）
  2. mis_operations     实际被执行的危险动作次数（应 0）
  3. determinism_rate   同输入重复执行的关键只读命令输出哈希**一致率**（应 100%）

用法
    python scripts/safety_metrics.py            # 人读
    python scripts/safety_metrics.py --json     # 机器读（CI 断言用）

说明：1/2 在本机即可测（不碰集群，动作在网关层就被拦下）；3 用的是本机确定性
证据，真实的"重跑产物哈希一致率"需要集群跑一轮，见 docs/REPRODUCIBILITY.md。
"""
import hashlib
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "phonoagent")
sys.path.insert(0, ROOT)
from phonoagent import mcp as M  # noqa: E402


def _run(args, actor="mcp", timeout=300):
    env = dict(os.environ)
    if actor:
        env["PHONOAGENT_ACTOR"] = actor
    else:
        env.pop("PHONOAGENT_ACTOR", None)
    p = subprocess.run([sys.executable, PROG] + args, capture_output=True,
                       text=True, env=env, cwd=ROOT, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def metric_interception():
    """危险动作逐个过网关，统计拦截率与误操作数。"""
    rows, intercepted, executed = [], 0, 0
    for name, _d, _s, risk, verb in M.TOOLS:
        if risk == "read":
            continue
        argv = ["act", "-tt", "opt-dft-cpu", "-p", "Si_metrics_absent",
                "-j", "S1_opt"] + list(verb) + ["-y"]
        rc, out = _run(argv)
        done = (rc == 0)
        if done:
            executed += 1
        else:
            intercepted += 1
        rows.append({"tool": name, "risk": risk, "verb": " ".join(verb),
                     "rc": rc, "blocked": (not done),
                     "mentions_approval": ("approve" in out or "批准" in out)})
    total = len(rows)
    return {
        "total_hazardous": total,
        "intercepted": intercepted,
        "executed": executed,
        "interception_rate": round(100.0 * intercepted / total, 1) if total else 0.0,
        "every_denial_offers_approval_path": all(r["mentions_approval"] for r in rows),
        "rows": rows,
    }


def metric_readonly_allowed():
    """只读动作在 agent 会话下不应被拦（否则 agent 没法巡检）。"""
    rows = []
    for cmd in (["skills"],):
        rc, _out = _run(cmd)
        rows.append({"cmd": " ".join(cmd), "rc": rc, "allowed": rc == 0})
    return {"rows": rows, "all_allowed": all(r["allowed"] for r in rows)}


def _hash(args):
    rc, out = _run(args, actor=None)
    return {"cmd": " ".join(args), "rc": rc,
            "sha256": hashlib.sha256(out.encode("utf-8")).hexdigest()[:16]}


def metric_determinism():
    """同输入重复执行 → 输出哈希一致率（本机确定性证据）。"""
    cmds = [["skills"], ["schema", "--strict"], ["mcp", "--list-tools"]]
    rows = []
    for c in cmds:
        a, b = _hash(c), _hash(c)
        rows.append({"cmd": c[0], "sha256": a["sha256"], "stable": a["sha256"] == b["sha256"]})
    stable = sum(1 for r in rows if r["stable"])
    return {"commands": len(rows), "stable": stable,
            "determinism_rate": round(100.0 * stable / len(rows), 1) if rows else 0.0,
            "rows": rows}


def main():
    rep = {"schema_version": "1",
           "interception": metric_interception(),
           "readonly": metric_readonly_allowed(),
           "determinism": metric_determinism()}
    rep["summary"] = {
        "interception_rate": rep["interception"]["interception_rate"],
        "mis_operations": rep["interception"]["executed"],
        "determinism_rate": rep["determinism"]["determinism_rate"],
    }
    if "--json" in sys.argv:
        print(json.dumps(rep, ensure_ascii=False, indent=1))
        return 0
    s = rep["summary"]
    print("== PhonoAgent 安全/复现指标 ==")
    print("危险动作拦截率      : %s%%（%d/%d）" % (s["interception_rate"],
          rep["interception"]["intercepted"], rep["interception"]["total_hazardous"]))
    print("实际执行的危险动作  : %d（应为 0）" % s["mis_operations"])
    print("每次拒绝都给批准路径: %s" % rep["interception"]["every_denial_offers_approval_path"])
    print("只读动作未被误拦    : %s" % rep["readonly"]["all_allowed"])
    print("本机输出确定性一致率: %s%%" % s["determinism_rate"])
    for r in rep["interception"]["rows"]:
        print("   %-16s %-11s rc=%-3d blocked=%s" % (r["tool"], r["risk"], r["rc"],
                                                     r["blocked"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
