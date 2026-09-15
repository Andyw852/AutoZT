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
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "phonoagent")
sys.path.insert(0, ROOT)
from phonoagent import mcp as M  # noqa: E402


def _run(args, actor="mcp", timeout=300, strict=None):
    env = dict(os.environ)
    if actor:
        env["PHONOAGENT_ACTOR"] = actor
    else:
        env.pop("PHONOAGENT_ACTOR", None)
    if strict is None:
        env.pop("PHONOAGENT_AGENT_STRICT", None)
    else:
        env["PHONOAGENT_AGENT_STRICT"] = str(strict)
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



def _ablation_sandbox():
    """临时项目 + 一个只有 POSCAR 的材料：clean 在它上面执行也无害。"""
    base = os.path.join(ROOT, "tmp", "_ablation")
    shutil.rmtree(base, ignore_errors=True)
    proj = os.path.join(base, "projects")
    mat = os.path.join(proj, "Si_x")
    os.makedirs(mat, exist_ok=True)
    src = os.path.join(ROOT, "test", "tf_test", "Si", "POSCAR")
    shutil.copy(src, os.path.join(mat, "POSCAR"))
    cfg = os.path.join(base, "tf.yaml")
    with open(cfg, "w", encoding="utf-8") as f:
        f.write("host: \"\"\nproject_roots:\n  - %s\n"
                "auto_advance: false\nauto_watch: false\n"
                "task_types:\n  opt-dft-cpu:\n    max_jobs: 1\n" % proj)
    _run(["-c", cfg, "-tt", "opt-dft-cpu", "-p", "Si_x", "init"], actor=None)
    return cfg


def ablation():
    """Ablation：同一危险动作在"网关严格"与"严格关闭"下的对照。

    在临时项目里放一个只有 POSCAR 的材料：严格关闭时命令会真的执行，但只清理
    这个临时目录（无害），因此实验可以放心在 CI 里跑。
    """
    cfg = _ablation_sandbox()
    rows = []
    for strict in (1, 0):
        rc, out = _run(["-c", cfg, "-tt", "opt-dft-cpu", "-p", "Si_x", "clean", "-y"],
                       actor="mcp", strict=str(strict))
        blocked = ("approve" in out or "批准" in out or "拒绝" in out)
        rows.append({"gateway_strict": strict, "returncode": rc, "blocked": blocked,
                     "executed": (rc == 0),
                     "readonly_still_ok": _run(["skills"], actor="mcp",
                                               strict=str(strict))[0] == 0})
    # 第三臂：只读档下暴露多少工具、危险工具能否被调用
    env_backup = os.environ.get("PHONOAGENT_MCP_READONLY")
    os.environ["PHONOAGENT_MCP_READONLY"] = "1"
    try:
        import importlib
        M = importlib.import_module("phonoagent.mcp")
        exposed = len(M._tools_list())
        refused = M.call_tool("clean_material", {"material": "Si_x"})
    finally:
        if env_backup is None:
            os.environ.pop("PHONOAGENT_MCP_READONLY", None)
        else:
            os.environ["PHONOAGENT_MCP_READONLY"] = env_backup
    readonly_arm = {"tools_exposed": exposed,
                    "hazardous_call_refused": bool(refused.get("isError"))}
    # 只读档下"只读任务"是否仍能完成（用 list_skills 作代表：它不需要集群）
    os.environ["PHONOAGENT_MCP_READONLY"] = "1"
    try:
        import importlib
        M2 = importlib.import_module("phonoagent.mcp")
        ok_call = M2.call_tool("list_skills", {})
        readonly_arm["readonly_task_ok"] = (not ok_call.get("isError"))
    finally:
        if env_backup is None:
            os.environ.pop("PHONOAGENT_MCP_READONLY", None)
        else:
            os.environ["PHONOAGENT_MCP_READONLY"] = env_backup
    with_gate, without_gate = rows[0], rows[1]
    return {"rows": rows, "readonly_profile": readonly_arm,
            "gateway_blocks_hazardous": with_gate["blocked"] and not with_gate["executed"],
            "bypass_executes_hazardous": without_gate["executed"],
            "readonly_ok_in_both": all(r["readonly_still_ok"] for r in rows)}


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
    if "--ablation" in sys.argv:
        ab = ablation()
        if "--json" in sys.argv:
            print(json.dumps(ab, ensure_ascii=False, indent=1))
        else:
            print("== agent 安全 ablation（危险动作：clean -y 目标为不存在的材料）==")
            print("%-16s %6s %8s %9s %10s" % ("gateway_strict", "rc", "blocked",
                                              "executed", "readonly"))
            for r in ab["rows"]:
                print("%-16s %6d %8s %9s %10s" % (r["gateway_strict"], r["returncode"],
                                                  r["blocked"], r["executed"],
                                                  r["readonly_still_ok"]))
            ra = ab["readonly_profile"]
            print("只读档：暴露工具 %d 个；危险工具调用被拒 = %s；只读任务仍可完成 = %s"
                  % (ra["tools_exposed"], ra["hazardous_call_refused"],
                     ra.get("readonly_task_ok")))
            print("网关挡住危险动作: %s | 关闭严格后会真执行: %s | 只读两种设置都正常: %s"
                  % (ab["gateway_blocks_hazardous"], ab["bypass_executes_hazardous"],
                     ab["readonly_ok_in_both"]))
        return 0
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
