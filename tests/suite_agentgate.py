# -*- coding: utf-8 -*-
"""v1.0 自测：agent 动作网关（autozt act / autozt approve）+ 审计流水（P0-1）。

全程本地沙盒：tmp/_agentgate_test/ 里的假项目 + 假配置（host 为空 = 只在本机跑），
不连任何超算、不碰真实材料。用 pty 模拟「人在交互终端里批准」。
"""
import json
import os
import pty
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "tmp", "_agentgate_test")
CFG = os.path.join(BASE, "tf.yaml")
PROJ = os.path.join(BASE, "projects")
LOG = os.path.join(BASE, ".tf_agent_log.jsonl")
APPR = os.path.join(BASE, ".tf_approvals.json")
PROG = os.path.join(ROOT, "bin", "autozt")
FAILS = []


def ok(cond, msg):
    if cond:
        print("  ✓ %s" % msg)
    else:
        FAILS.append(msg)
        print("  ✗ %s" % msg)


def setup():
    shutil.rmtree(BASE, ignore_errors=True)
    mat = os.path.join(PROJ, "Si_ag")
    os.makedirs(mat, exist_ok=True)
    shutil.copy(os.path.join(ROOT, "test", "tf_test", "Si", "POSCAR"),
                os.path.join(mat, "POSCAR"))
    with open(CFG, "w", encoding="utf-8") as f:
        f.write("host: \"\"\n"
                "project_roots:\n  - %s\n"
                "auto_advance: false\nauto_watch: false\n"
                "task_types:\n  opt-dft-cpu:\n    max_jobs: 2\n" % PROJ)
    env = dict(os.environ)
    env.pop("AUTOZT_ACTOR", None)
    subprocess.run([sys.executable, PROG, "-c", CFG, "-tt", "opt-dft-cpu",
                    "-p", "Si_ag", "init"], capture_output=True, text=True, env=env)


def run(args, actor=None, strict=None, ttl=None, tty=False, answer="y\n"):
    env = dict(os.environ)
    env.pop("AUTOZT_AGENT_STRICT", None)
    env.pop("AUTOZT_APPROVE_TTL", None)
    if actor:
        env["AUTOZT_ACTOR"] = actor
    else:
        env.pop("AUTOZT_ACTOR", None)
    if strict is not None:
        env["AUTOZT_AGENT_STRICT"] = str(strict)
    if ttl is not None:
        env["AUTOZT_APPROVE_TTL"] = str(ttl)
    argv = [sys.executable, PROG, "-c", CFG] + [str(a) for a in args]
    if not tty:
        p = subprocess.run(argv, capture_output=True, text=True, env=env,
                           timeout=180)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    m, s = pty.openpty()                      # 真 pty：让 isatty() 为真
    p = subprocess.Popen(argv, stdin=s, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, env=env)
    os.close(s)
    time.sleep(0.6)
    try:
        os.write(m, answer.encode())
    except OSError:
        pass
    out = p.communicate(timeout=180)[0].decode("utf-8", "replace")
    os.close(m)
    return p.returncode, out


def log_recs():
    if not os.path.isfile(LOG):
        return []
    out = []
    with open(LOG, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    pass
    return out


def approvals():
    if not os.path.isfile(APPR):
        return []
    with open(APPR, encoding="utf-8") as f:
        return (json.load(f) or {}).get("approvals") or []


print("[0] 沙盒准备")
setup()
ok(os.path.isfile(CFG), "假配置就绪：tmp/_agentgate_test/tf.yaml")
ok(not os.path.exists(LOG), "还没审计日志（干净起点）")

print("[1] 零影响：不设 AUTOZT_ACTOR 时，一切照旧（不拦、不记）")
rc, out = run(["-tt", "opt-dft-cpu", "-p", "Si_ag", "clean", "-y"])
ok(rc == 0 and "已清理" in out, "直连 clean -y 正常执行（rc=%d）" % rc)
ok(not os.path.exists(LOG), "非 agent 会话不写审计日志")

print("[2] agent 直连 + 严格（缺省）：破坏性动作被拒，且**没**执行")
n0 = len(log_recs())
rc, out = run(["-tt", "opt-dft-cpu", "-p", "Si_ag", "clean", "-y"], actor="claude")
ok(rc == 3, "被拒绝（rc=3）")
ok("已清理" not in out, "命令没有被执行（输出里没有实际动作）")
ok("autozt approve" in out, "提示里给了人工批准的命令")
recs = log_recs()
ok(len(recs) == n0 + 1 and recs[-1]["decision"] == "deny-need-approval"
   and recs[-1]["risk"] == "destructive", "审计记了一条 deny-need-approval")

print("[3] AUTOZT_AGENT_STRICT=0：agent 直连放行（给用户自己的脚本留的后门）")
rc, out = run(["-tt", "opt-dft-cpu", "-p", "Si_ag", "clean", "-y"],
              actor="claude", strict=0)
ok(rc == 0 and "已清理" in out, "关闭严格后放行且真执行（rc=%d）" % rc)
ok(log_recs()[-1]["decision"] == "direct", "审计记 decision=direct")

print("[4] 网关只读命令：放行 + 记账（含退出码）")
rc, out = run(["act", "-p", "Si_ag", "list"], actor="claude")
ok(rc == 0 and "Si_ag" in out, "autozt act list 正常输出（rc=%d）" % rc)
r = log_recs()[-1]
ok(r["risk"] == "read" and r["decision"] == "allow" and r["exit_code"] == 0,
   "审计：read / allow / exit_code=0")

print("[5] 网关破坏性命令：无批准 → 拒绝（exit 3）")
rc, out = run(["act", "-tt", "opt-dft-cpu", "-p", "Si_ag", "clean", "-y"],
              actor="claude")
ok(rc == 3 and "拒绝执行" in out, "拒绝并说明理由")
ok("命令签名" in out and "autozt approve" in out, "给出签名与批准命令")
ok(log_recs()[-1]["decision"] == "deny-need-approval", "审计记拒绝")

print("[6] autozt approve 必须交互终端（非 TTY 一律拒绝）")
rc, out = run(["approve", "-tt", "opt-dft-cpu", "-p", "Si_ag", "clean", "-y"])
ok(rc == 3 and "交互终端" in out, "非 TTY 拒绝批准（rc=3）")
ok(approvals() == [], "批准文件里没有条目")

print("[7] autozt approve 在 pty 里确认 y → 写入一次性令牌")
rc, out = run(["approve", "-tt", "opt-dft-cpu", "-p", "Si_ag", "clean", "-y"],
              tty=True, answer="y\n")
ok(rc == 0 and "已批准" in out, "批准成功（rc=%d）" % rc)
ap = approvals()
ok(len(ap) == 1 and ap[0]["used"] is False and len(ap[0]["sig"]) == 12,
   "令牌已入库（签名 12 位、未使用）")

print("[8] agent 拿令牌跑同一条命令（参数顺序颠倒也算同一条）→ 放行")
rc, out = run(["act", "clean", "-y", "-tt", "opt-dft-cpu", "-p", "Si_ag"],
              actor="claude")
ok(rc == 0 and "已清理" in out, "参数顺序不同、签名一致 → 放行并真执行")
r = log_recs()[-1]
ok(r["decision"] == "allow" and r["approved_by"] and r["exit_code"] == 0,
   "审计：allow / approved_by=%s / exit_code=0" % r["approved_by"])
ok(all(x["used"] for x in approvals()), "令牌已销（一次性）")

print("[9] 令牌用完即废：同一命令再来一次 → 又拒绝")
rc, out = run(["act", "clean", "-y", "-tt", "opt-dft-cpu", "-p", "Si_ag"],
              actor="claude")
ok(rc == 3 and "拒绝执行" in out, "无令牌可用 → 拒绝（不可重放）")

print("[10] 令牌会过期（AUTOZT_APPROVE_TTL=1s）")
rc, _ = run(["approve", "-tt", "opt-dft-cpu", "-p", "Si_ag", "clean", "-y"],
            tty=True, answer="y\n", ttl=1)
ok(rc == 0, "先批准（TTL=1s）")
time.sleep(2.2)
rc, out = run(["act", "clean", "-y", "-tt", "opt-dft-cpu", "-p", "Si_ag"],
              actor="claude", ttl=1)
ok(rc == 3, "过期后拒绝（rc=3）")

print("[11] --dry-run 降级为只读：破坏性命令也能安全排练")
step_dir = os.path.join(PROJ, "Si_ag", "opt-dft-cpu", "step1_std_opt")
os.makedirs(step_dir, exist_ok=True)
keep = os.path.join(step_dir, "KEEP.txt")
with open(keep, "w") as f:
    f.write("x\n")
rc, out = run(["act", "-tt", "opt-dft-cpu", "-p", "Si_ag", "rerun", "--dry-run"],
              actor="claude")
ok(rc == 0 and "dry-run" in out, "dry-run 放行（rc=%d）" % rc)
ok(log_recs()[-1]["risk"] == "read", "审计把 dry-run 记为 read")
ok(os.path.isfile(keep), "排练没有删掉任何东西")

print("[12] -y / -f 一律算破坏性（需批准）")
rc, _ = run(["act", "-tt", "opt-dft-cpu", "-p", "Si_ag", "retry", "-y"],
            actor="claude")
ok(rc == 3, "retry -y 被拦（rc=3）")
rc, _ = run(["act", "-tt", "opt-dft-cpu", "-p", "Si_ag", "init", "-f"],
            actor="claude")
ok(rc == 3, "init -f 被拦（rc=3）")

print("[13] 参数位置：-p 写在 act 之前 = 拒绝（防止批准与执行不是同一条）")
rc, out = run(["-p", "Si_ag", "act", "clean", "-y"], actor="claude")
ok(rc == 2 and "act 之后" in out, "拒绝并指路（rc=2）")

print("[14] act 子命令：policy / log")
rc, out = run(["act", "policy"])
ok(rc == 0 and "破坏性" in out and "只读" in out and "推进" in out,
   "act policy 打出三档")
rc, out = run(["act", "log"], actor="claude")
ok(rc == 0 and "agent 审计" in out and "退出码" in out, "act log 打出审计流水")
ok("approve" not in out.split("时间")[0], "log 与 policy 不串")

print("[15] 不设 AUTOZT_ACTOR 也能用网关（actor 退化成 USER 环境变量）")
rc, out = run(["act", "-p", "Si_ag", "summary"])
ok(rc == 0, "act summary 正常（rc=%d）" % rc)
ok(log_recs()[-1]["actor"] == (os.environ.get("USER") or ""),
   "actor 退化为 $USER（%s）" % os.environ.get("USER"))

print("[16] act 不能嵌套，签名对 -y 不敏感")
rc, out = run(["act", "act", "list"], actor="claude")
ok(rc == 2 and "嵌套" in out, "act act 被拒（rc=2）")
sys.path.insert(0, ROOT)
from autozt import agent_signature as _sig          # noqa: E402
ok(_sig("clean", ["-p", "A", "clean", "-y"]) == _sig("clean", ["clean", "-p", "A"]),
   "签名忽略参数顺序与 -y")

print("")
if FAILS:
    print("FAILED %d 项：" % len(FAILS))
    for m in FAILS:
        print("  - %s" % m)
    sys.exit(1)
print("ALL PASS（v1.0 agentgate 网关与审计自测）")
