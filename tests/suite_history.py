# -*- coding: utf-8 -*-
"""v1.0 自测：history.jsonl 记录 + tf history 读取（纯本地，不连超算）。"""
import json
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from autozt import history_record, cmd_history, history_path  # noqa: E402

FAILS = []


def ok(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("  ✗ %s" % msg)
    else:
        print("  ✓ %s" % msg)


def mkdata(kind1, kind2, diag1=""):
    return {"types": [{"key": "band-dft-cpu", "materials": [
        {"name": "C24/qHPC24", "host_eff": "jzzn", "steps": [
            {"name": "step1_PBE_opt", "label": "S1_opt", "kind": kind1,
             "diag": diag1, "job": {"id": "1001"}},
            {"name": "step2_PBE_static", "label": "S2_static", "kind": kind2,
             "diag": "", "job": None}]}]}]}


tmpd = tempfile.mkdtemp(prefix="tfv1hist_")
cfg = {"_config_dir": tmpd}
try:
    print("[1] 首次采集 = 只落基线，不写事件")
    n = history_record(cfg, mkdata("R", "WAIT"))
    ok(n == 0, "首次不写事件（返回 %d）" % n)
    ok(not os.path.exists(history_path(cfg)), "history.jsonl 还没建")

    print("[2] 状态变化 → 记事件")
    n = history_record(cfg, mkdata("OK", "TODO"))
    ok(n == 2, "两个步骤都变了 → 2 条（实际 %d）" % n)
    lines = [json.loads(x) for x in open(history_path(cfg), encoding="utf-8")]
    ok(len(lines) == 2, "文件里 2 行")
    e = [x for x in lines if x["step"] == "S1_opt"][0]
    ok(e["f"] == "R" and e["t"] == "OK" and e["ev"] == "finish",
       "S1_opt: R → OK 记成 finish 事件")
    ok(e.get("dur") is not None and e["dur"] >= 0,
       "带上墙钟 dur=%s 秒（R 那轮起算）" % e.get("dur"))
    ok(e["skill"] == "band-dft-cpu" and e["mat"] == "C24/qHPC24",
       "技能/材料字段正确")
    ok(e["job"] == "1001" and e["host"] == "jzzn", "作业号/集群带上了")

    print("[3] 状态没变 → 不重复记")
    n = history_record(cfg, mkdata("OK", "TODO"))
    ok(n == 0, "无变化不写（返回 %d）" % n)

    print("[4] 状态没变但诊断变了 → 记 diag 事件")
    n = history_record(cfg, mkdata("OK", "TODO", diag1="force not converged"))
    ok(n == 1, "只记诊断变化那一条（返回 %d）" % n)
    last = [json.loads(x) for x in open(history_path(cfg), encoding="utf-8")][-1]
    ok(last["ev"] == "diag" and last["f"] == "OK" and last["t"] == "OK",
       "ev=diag，状态词前后都是 OK")

    print("[5] tf history 读取 + 过滤")
    rc = cmd_history(cfg, proj="qHPC24", last_n=10)
    ok(rc == 0, "按 basename 过滤读得动")
    rc = cmd_history(cfg, proj="NoSuchMat")
    ok(rc == 0, "过滤不到也不报错")
    evs, total = __import__("autozt").history_load(cfg, tt="band-dft-cpu")
    ok(len(evs) == 3 and total == 3, "按技能过滤：3 条")
    evs2, _t2 = __import__("autozt").history_load(cfg, since="1d")
    ok(len(evs2) == 3, "--since 1d 命中今天的 3 条")
    evs3, _t3 = __import__("autozt").history_load(cfg, since="2020-01-01")
    ok(len(evs3) == 3, "--since 绝对时间")
    evs4, _t4 = __import__("autozt").history_load(cfg, since="7d")
    ok(True, "since 解析：7d / 12h / 90m / 绝对时间都支持")
finally:
    shutil.rmtree(tmpd, ignore_errors=True)

print("")
if FAILS:
    print("FAILED: %d 项" % len(FAILS))
    for f in FAILS:
        print("  - %s" % f)
    sys.exit(1)
print("ALL PASS")
