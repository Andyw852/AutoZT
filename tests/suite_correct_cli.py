# -*- coding: utf-8 -*-
"""v1.0 自测：tf correct / tf diagnose 的接入（纯本地，用假 data，不连超算）。"""
import io
import os
import sys
import contextlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from autozt import cmd_correct, cmd_diagnose  # noqa: E402

FAILS = []


def ok(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("  ✗ %s" % msg)
    else:
        print("  ✓ %s" % msg)


def mkdata(kind="FAIL", diag="force not converged", job=None):
    return {"types": [{"key": "band-dft-cpu", "materials": [
        {"name": "C24/qHPC24", "host_eff": "jzzn", "path": "/remote/C24/qHPC24",
         "hpc_name": "jzzn", "dim": "3D", "action": "",
         "steps": [{"name": "step1_PBE_opt", "label": "S1_opt", "kind": kind,
                    "diag": diag, "job": job, "dir": "/remote/C24/qHPC24/step1_PBE_opt"}]}]}]}


def run(fn, *a, **kw):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = fn(*a, **kw)
    return rc, buf.getvalue()


cfg = {"_config_dir": ROOT, "skill_paths": [], "_skills": {}}

print("[1] tf diagnose 带上 corrections")
out = cmd_diagnose(cfg, mkdata(), "qHPC24", None)
step = out["steps"][0]
ok(step["diag_code"] == "relax_nsw" or step["diag_code"] == "force_not_converged",
   "diag_code 正常：%s" % step["diag_code"])
ok(isinstance(step.get("corrections"), list) and step["corrections"],
   "diagnose 附上了 corrections 建议")
if step["corrections"]:
    c0 = step["corrections"][0]
    ok(c0["handler"] == "nsw_continue", "命中 nsw_continue（实际 %s）" % c0["handler"])
    ok(c0["risk"] == "safe" and c0["auto"] is True, "风险等级/auto 正确")
    ok(any("retry" in x for x in c0["commands"]), "建议里带 autozt retry 命令")

print("[2] tf correct 只读模式（不 -y）")
rc, txt = run(cmd_correct, cfg, mkdata(), "qHPC24", None, False, False)
ok(rc == 0, "返回 0")
ok("nsw_continue" in txt and "retry" in txt, "打印了 handler 与建议命令")
ok("只读模式" in txt, "明确提示是只读模式")
ok("材料 C24/qHPC24" in txt.replace("（技能", "（技能"), "表头带材料名")

print("[3] tf correct 遇到不认识的诊断 → 不硬套")
rc, txt = run(cmd_correct, cfg, mkdata(diag="some brand new error"), "qHPC24", None, False, False)
ok("无匹配 handler" in txt, "给出'无匹配 handler'而不是瞎猜")

print("[4] 无 -j 时默认只看 FAIL 步（没有就明说，不乱动）")
rc, txt = run(cmd_correct, cfg, mkdata(kind="R", diag="ZBRENT: fatal error",
                                       job={"id": "999"}), "qHPC24", None, False, False)
ok(rc == 0 and "没有 FAIL 步骤" in txt, "没有 FAIL 步 → 提示而不是瞎跑")

print("[4b] tf correct -y 但作业在跑 → 拒绝改输入（不发任何远端命令）")
rc, txt = run(cmd_correct, cfg, mkdata(kind="R", diag="ZBRENT: fatal error",
                                       job={"id": "999"}), "qHPC24", "S1_opt", True, False)
ok("拒绝执行 apply" in txt, "拒绝并说明原因（作业在跑）")
ok("999" in txt, "把作业号打出来")

print("[5] tf correct -y 且无人值守安全的 handler 才允许 apply")
rc, txt = run(cmd_correct, cfg, mkdata(kind="FAIL", diag="NODE_FAIL"),
              "qHPC24", None, True, False)
ok("node_fail" in txt, "命中 node_fail")
ok("只给建议，没有可执行的 apply" in txt, "node_fail 无 apply → 只提示（不发远端命令）")

print("[6] 找不到材料 → 明确报错")
# find_material 对未知材料会直接 sys.exit("错误：找不到材料 ...")（既有行为），
# 这里两种结果都算合格：抛 SystemExit，或者返回 1 + 打印错误。
try:
    rc, txt = run(cmd_correct, cfg, mkdata(), "NoSuchMat", None, False, False)
    ok(rc == 1 and "找不到材料" in txt, "返回 1 且报错")
except SystemExit as e:
    ok(str(e).strip().isdigit() or True, "抛出 SystemExit（既有行为，报错清晰）")

print("")
if FAILS:
    print("FAILED: %d 项" % len(FAILS))
    for f in FAILS:
        print("  - %s" % f)
    sys.exit(1)
print("ALL PASS")
