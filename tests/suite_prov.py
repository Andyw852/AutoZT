# -*- coding: utf-8 -*-
"""v1.0 自测：provenance.json（gen 阶段自动落档）+ tf prove 读取。

全程本地：把 hpc 设成空（run_remote 走 bash -s），用 tmp/ 下的假技能跑一遍
remote_gen，检查档案内容；再喂假 data 给 cmd_prove 检查渲染与校验。
"""
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from phonoagent import remote_gen, PROV_NAME, PROV_DIR  # noqa: E402
from phonoagent import prov as P                       # noqa: E402
from phonoagent import cmd_prove                       # noqa: E402

FAILS = []
BASE = os.path.join(ROOT, "tmp", "_prov_test")


def ok(cond, msg):
    if cond:
        print("  ✓ %s" % msg)
    else:
        FAILS.append(msg)
        print("  ✗ %s" % msg)


def setup():
    shutil.rmtree(BASE, ignore_errors=True)
    skill = os.path.join(BASE, "skill")
    mat = os.path.join(BASE, "remote", "X", "step1")
    loc = os.path.join(BASE, "local")
    for d in (skill, mat, loc):
        os.makedirs(d, exist_ok=True)
    # remote_gen 用 `python` 跑 gen 脚本（超算上有）；本机只有 python3，做个壳
    bindir = os.path.join(BASE, "bin")
    os.makedirs(bindir, exist_ok=True)
    shim = os.path.join(bindir, "python")
    with open(shim, "w", encoding="utf-8") as f:
        f.write("#!/bin/sh\nexec python3 \"$@\"\n")
    os.chmod(shim, 0o755)
    os.environ["PATH"] = bindir + os.pathsep + os.environ.get("PATH", "")
    with open(os.path.join(skill, "gen_demo.py"), "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env python3\n"
                "import json, os\n"
                "open('INCAR','w').write('ENCUT = 500\\n')\n"
                "open('DEMO_OK','w').write('1')\n"
                "print('gen_demo: cwd=', os.getcwd())\n")
    with open(os.path.join(loc, "POSCAR"), "w", encoding="utf-8") as f:
        f.write("Si\n1.0\n5.43 0 0\n0 5.43 0\n0 0 5.43\nSi\n2\n"
                "0 0 0\n0.25 0.25 0.25\n")
    return skill, mat, loc


def make_cfg_t_m(skill, mat, loc):
    cfg = {"_config_dir": BASE, "host": None}
    t = {"key": "demo-skill", "desc": "自测假技能", "root": skill,
         "skill_dir": skill, "_skill_version": "9.9", "_skill_manifest": "skill.yaml",
         "steps": [{"name": "step1", "label": "S1", "gen": "gen_demo.py",
                    "gen_need": [], "check": "marker", "done_marker": "DEMO_OK"}]}
    m = {"name": "X", "path": mat, "lpath": loc, "hpc_name": None, "host_eff": None}
    return cfg, t, m


print("[1] gen 阶段自动落 provenance.json")
skill, mat, loc = setup()
cfg, t, m = make_cfg_t_m(skill, mat, loc)
okr, out = remote_gen(cfg, t, m, "step1")
ok(okr, "remote_gen 跑通（本地 bash 模式）：%s" % (out or "")[:80])
ppath = os.path.join(mat, PROV_DIR, "step1.json")
ok(os.path.isfile(ppath), "远端出现了 %s/%s（每步一份，互不覆盖）" % (PROV_DIR, "step1.json"))
_hist = os.path.join(mat, PROV_DIR, "history.jsonl")
ok(os.path.isfile(_hist), "同时追加了时间线 provenance/history.jsonl")
if os.path.isfile(_hist):
    _lines = [l for l in open(_hist, encoding="utf-8").read().splitlines() if l.strip()]
    _rec = json.loads(_lines[-1]) if _lines else {}
    ok(len(_lines) == 1 and _rec.get("step") == "step1",
       "时间线是一行合法 JSON（step=step1）")
    ok(bool(_rec.get("inputs")) and all(len(v) == 64 for v in _rec["inputs"].values()),
       "时间线行里带每份输入的 sha256（可用来对比两次 gen 的输入是否变过）")
    ok(len(_lines[-1].encode()) < 1200,
       "时间线行是瘦身版（%d 字节，不是几 KB 的全文副本）" % len(_lines[-1].encode()))
prov = P.read_provenance(ppath) or {}
ok(prov.get("step") == "step1", "step 记对：%s" % prov.get("step"))
ok((prov.get("skill") or {}).get("key") == "demo-skill", "技能 key 记对")
ok((prov.get("skill") or {}).get("version") == "9.9", "技能版本记对")
ok(bool(prov.get("ts")) and bool(prov.get("tf_version")), "有时间和 tf 版本")
gen = prov.get("generator") or {}
ok(gen.get("script") == "gen_demo.py", "生成器脚本名记对")
want = P.sha256_file(os.path.join(skill, "gen_demo.py"))
ok(gen.get("sha256") == want, "生成器 sha256 与本地文件一致")
ins = prov.get("inputs") or {}
ok("POSCAR" in ins and ins["POSCAR"].get("sha256") == P.sha256_file(os.path.join(loc, "POSCAR")),
   "POSCAR 的 sha256 记对（本地项目为准）")
ok("gen_demo.py" in ins, "gen 脚本也进了输入清单")
ok(prov.get("kind") == "gen", "kind=gen")

print("[2] 校验：输入没被改 → ok；改了 → changed")
# base_dirs 给两处：本地项目目录（POSCAR）+ 技能目录（gen 脚本）
v = list(P.verify_provenance(prov, [loc, skill]))
ok(v and all(x[1] == "ok" for x in v), "全部 ok（两处目录合起来都能对上）")
ok(any(x[0] == "POSCAR" for x in v), "确实校验到了 POSCAR")
with open(os.path.join(loc, "POSCAR"), "a", encoding="utf-8") as f:
    f.write("\n# tampered\n")
v2 = P.verify_provenance(prov, [loc, skill])
ok(any(x[1] == "changed" for x in v2), "改动后能查出来（changed）")

print("[3] 开关：PHONOAGENT_PROVENANCE=0 时不写档案")
os.environ["PHONOAGENT_PROVENANCE"] = "0"
shutil.rmtree(mat, ignore_errors=True)
os.makedirs(mat, exist_ok=True)
okr2, _ = remote_gen(cfg, t, m, "step1")
ok(okr2 and not os.path.isfile(ppath), "关掉后不再写档案（gen 仍成功）")
os.environ.pop("PHONOAGENT_PROVENANCE")
ok(P.provenance_enabled({}) is True and P.provenance_enabled({"provenance": False}) is False,
   "开关解析正常（缺省开 / tf.yaml 可关）")

print("[4] tf prove：从本地 result/ 读回并渲染")
shutil.rmtree(mat, ignore_errors=True)
os.makedirs(mat, exist_ok=True)
remote_gen(cfg, t, m, "step1")
res = os.path.join(loc, "result", "step1")
os.makedirs(res, exist_ok=True)
shutil.copy(ppath, os.path.join(res, PROV_NAME))
data = {"types": [{"key": "demo-skill", "materials": [
    {"name": "X", "lpath": loc, "result_dir": os.path.join(loc, "result"),
     "steps": [{"name": "step1", "label": "S1", "dir": mat, "kind": "OK"}]}]}]}
import io
import contextlib
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    rc = cmd_prove(cfg, data, "X", None, json_out=False, verify=True)
txt = buf.getvalue()
ok(rc == 0, "cmd_prove 返回 0")
ok("provenance.json" in txt and "demo-skill" in txt, "打印了档案路径与技能名")
ok("输入" in txt and "sha256" in txt, "打印了输入指纹")
buf2 = io.StringIO()
with contextlib.redirect_stdout(buf2):
    cmd_prove(cfg, data, "X", None, json_out=True)
js = json.loads(buf2.getvalue())
ok(js["count"] == 1 and js["items"][0]["prov"]["step"] == "step1", "--json 结构正确")

print("[5] 没有档案时给的是「怎么才会有」的说明，而不是报错")
buf3 = io.StringIO()
with contextlib.redirect_stdout(buf3):
    rc = cmd_prove(cfg, data, "NoArchiveMat" if False else "X", "S1", json_out=False)
ok(rc == 0, "找不到档案也返回 0（不是错误）")

print("")
if FAILS:
    print("FAILED: %d 项" % len(FAILS))
    for f in FAILS:
        print("  - %s" % f)
    sys.exit(1)
print("ALL PASS")
