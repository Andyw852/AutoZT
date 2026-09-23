# -*- coding: utf-8 -*-
"""v1.0 自测：会话导出 autozt session export（P1-7）。

纯本地：tmp/_session_test/ 里造一份假材料 + 假 history.jsonl + 假 agent 审计 +
假 provenance.json，然后调 cmd_session 打包，逐项验证包内容、sha256、过滤与
「--json 不写包」。
"""
import io
import json
import os
import shutil
import sys
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from autozt import cmd_session, agent_audit                    # noqa: E402

BASE = os.path.join(ROOT, "tmp", "_session_test")
FAILS = []


def ok(cond, msg):
    if cond:
        print("  ✓ %s" % msg)
    else:
        FAILS.append(msg)
        print("  ✗ %s" % msg)


MAT = "Si_sess"
MDIR = os.path.join(BASE, "projects", MAT)


def setup():
    shutil.rmtree(BASE, ignore_errors=True)
    os.makedirs(os.path.join(MDIR, "result", "provenance"), exist_ok=True)
    with open(os.path.join(MDIR, "result", "provenance", "S1_opt.json"),
              "w", encoding="utf-8") as f:
        json.dump({"step": "S1_opt", "ts": "2026-09-14T10:00:00",
                   "inputs": {"INCAR": {"sha256": "ab" * 32, "bytes": 120}},
                   "runner": {"actor": "user"}}, f, ensure_ascii=False)
    evs = [
        {"ts": "2026-09-14T10:00:01", "mat": MAT, "skill": "opt-dft-cpu",
         "step": "S1_opt", "ev": "action", "act": "init", "note": "init S1_opt",
         "host": "jzzn"},
        {"ts": "2026-09-14T10:00:09", "mat": MAT, "skill": "opt-dft-cpu",
         "step": "S1_opt", "ev": "state", "f": "PREP", "t": "PD", "diag": ""},
        {"ts": "2026-09-14T10:00:48", "mat": MAT, "skill": "opt-dft-cpu",
         "step": "S1_opt", "ev": "finish", "f": "R", "t": "OK", "dur": 42},
        {"ts": "2026-09-14T11:00:00", "mat": "别的材料", "skill": "opt-dft-cpu",
         "step": "S1_opt", "ev": "state", "f": "R", "t": "OK"},
    ]
    with open(os.path.join(BASE, "history.jsonl"), "w", encoding="utf-8") as f:
        for e in evs:
            f.write(json.dumps(e, ensure_ascii=False) + "\n")
    cfg = {"_config_dir": BASE, "_config_path": os.path.join(BASE, "tf.yaml"),
           "host": "jzzn", "task_types": {"opt-dft-cpu": {"work_dir": "/pub/x"}}}
    agent_audit(cfg, "claude", "clean", ["-p", MAT, "clean", "-y"],
                "destructive", "deny-need-approval", exit_code=3)
    agent_audit(cfg, "claude", "list", ["-p", MAT, "list"], "read", "allow",
                exit_code=0)
    return cfg


def mkm():
    return {"name": MAT, "host_eff": "jzzn",
            "result_dir": os.path.join(MDIR, "result"),
            "lpath": MDIR,
            "steps": [
                {"name": "S1_opt", "label": "S1_opt", "kind": "OK", "diag": "",
                 "job": {"id": "3838838", "state": "COMPLETED"},
                 "dir": "/public/home/.../tf_smoke/Si_sess/opt-dft-cpu/step1_std_opt"},
                {"name": "S2_static", "label": "S2_static", "kind": "PREP",
                 "diag": "", "job": None,
                 "dir": "/public/home/.../tf_smoke/Si_sess/opt-dft-cpu/step2_static"}]}


def mkdata():
    return {"types": [{"key": "opt-dft-cpu", "materials": [mkm()]}]}


print("[0] 沙盒准备")
cfg = setup()
ok(os.path.isfile(os.path.join(BASE, "history.jsonl")), "假 history.jsonl 就绪")
ok(os.path.isfile(os.path.join(BASE, ".tf_agent_log.jsonl")), "假 agent 审计就绪")

print("[1] 打包")
out = os.path.join(BASE, "sess.tar.gz")
rc = cmd_session(cfg, mkdata(), MAT, out=out)
ok(rc == 0, "cmd_session 返回 0")
ok(os.path.isfile(out), "包已生成：%s" % os.path.basename(out))

print("[2] 包内容")
names = []
with tarfile.open(out) as tf_:
    names = sorted(tf_.getnames())
    blobs = {n: tf_.extractfile(n).read().decode("utf-8") for n in names}
for n in ("manifest.json", "README.txt", "status.txt", "history.jsonl",
          "agent_log.jsonl", "provenance/S1_opt.json"):
    ok(n in names, "包含 %s" % n)
ok("别的材料" not in blobs["history.jsonl"], "history 按材料过滤（不含别的材料）")
ok(blobs["history.jsonl"].count("\n") == 3, "history 收到该材料 3 条事件")

print("[3] manifest 字段与计数")
man = json.loads(blobs["manifest.json"])
ok(man["material"] == MAT and man["skill"] == "opt-dft-cpu"
   and man["host"] == "jzzn", "材料/技能/超算正确")
ok(man["counts"]["history_events"] == 3
   and man["counts"]["action_events"] == 1
   and man["counts"]["finish_events"] == 1, "事件计数正确（3/1/1）")
ok(man["counts"]["agent_calls"] == 2 and man["counts"]["agent_denied"] == 1,
   "agent 审计计数正确（2 次调用、1 次拒绝）")
ok(man["counts"]["provenance_steps"] == 1, "provenance 步数正确")
ok([s["label"] for s in man["steps"]] == ["S1_opt", "S2_static"],
   "步骤清单与状态入账")
ok(man["tf_version"] and man["created"] and man["replay"]["verify"],
   "带 tf 版本/时间/重放指引")
ok("git" in man and "rev" in (man.get("git") or {}), "带 git 提交号")

print("[4] 逐文件 sha256 可校验")
import hashlib                                             # noqa: E402
entries = {e["path"]: e for e in man["files"]}
bad = []
for path, e in entries.items():
    if e.get("self"):
        continue
    got = hashlib.sha256(blobs[path].encode("utf-8")).hexdigest()
    if got != e["sha256"] or len(blobs[path].encode("utf-8")) != e["bytes"]:
        bad.append(path)
ok(not bad, "所有文件 sha256/字节数与 manifest 一致（%d 个）"
   % (len(entries) - 1))
ok(entries.get("manifest.json", {}).get("self") is True,
   "manifest 自引用条目明确标注 self")

print("[5] README 讲清楚怎么读/怎么重放")
rd = blobs["README.txt"]
ok(MAT in rd and "provenance" in rd and "tf prove" in rd,
   "README 含材料名、provenance 说明与重放命令")

print("[6] status.txt 是人看的步骤表")
ok("S1_opt" in blobs["status.txt"] and "3838838" in blobs["status.txt"],
   "status.txt 带步骤与作业号")

print("[7] --json 只打印总账、不写包")
out2 = os.path.join(BASE, "should_not_exist.tar.gz")
rc = cmd_session(cfg, mkdata(), MAT, out=out2, json_out=True)
ok(rc == 0 and not os.path.isfile(out2), "--json 不写包")

print("[8] --since 过滤")
rc = cmd_session(cfg, mkdata(), MAT, out=os.path.join(BASE, "s2.tar.gz"),
                 since="2030-01-01")
with tarfile.open(os.path.join(BASE, "s2.tar.gz")) as tf_:
    h = tf_.extractfile("history.jsonl").read().decode("utf-8")
ok(h.strip() == "", "未来的 --since 过滤掉全部事件")

print("[9] 材料名不存在 → 返回 1（不抛栈）")
rc = cmd_session(cfg, mkdata(), "不存在的材料", out=os.path.join(BASE, "x.tgz"))
ok(rc == 1, "找不到材料时返回 1")

print("[10] 缺 -p → 返回 2 并给用法")
rc = cmd_session(cfg, mkdata(), None)
ok(rc == 2, "缺材料名返回 2")

print("")
if FAILS:
    print("FAILED %d 项：" % len(FAILS))
    for m in FAILS:
        print("  - %s" % m)
    sys.exit(1)
print("ALL PASS（v1.0 session export 自测）")
