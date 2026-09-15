# -*- coding: utf-8 -*-
"""v1.0 自测：新技能 drop-in（P1-9）。

把 skill/_template/ 原样复制到两个"投放点"——① 配置目录旁的 skill/，
② tf.yaml 的 skill_paths 指向的另一个目录——然后断言：

  · 被自动发现（autozt skills 认得，版本/步骤数正确）；
  · tf skill show 能渲染出卡片（论文图 2 那套版式，零代码改动）；
  · tf schema --strict 不报错（自描述三段合法）；
  · 它的 step1 真能跑（gen → done_marker 落地 + provenance 自动落档）；
  · **核心代码零改动**：跑完前后 autozt/*.py 与 bin/autozt 的 sha256 一个没变。

全程本地（type hpc 指向不存在的 local 档 + ssh_host 置空 → gen 走 bash），
不连超算、不碰任何真实材料。
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "tmp", "_dropin_test")
CFG = os.path.join(BASE, "tf.yaml")
PROJ = os.path.join(BASE, "projects")
SHIM = os.path.join(BASE, "bin")
PROG = os.path.join(ROOT, "bin", "autozt")
FAILS = []


def ok(cond, msg):
    if cond:
        print("  ✓ %s" % msg)
    else:
        FAILS.append(msg)
        print("  ✗ %s" % msg)


def core_hashes():
    """核心代码指纹：autozt/*.py + bin/autozt（技能是数据，核心是代码——要分清）。"""
    out = {}
    d = os.path.join(ROOT, "autozt")
    for fn in sorted(os.listdir(d)):
        if fn.endswith(".py"):
            with open(os.path.join(d, fn), "rb") as f:
                out["autozt/" + fn] = hashlib.sha256(f.read()).hexdigest()
    with open(os.path.join(ROOT, "bin", "autozt"), "rb") as f:
        out["bin/autozt"] = hashlib.sha256(f.read()).hexdigest()
    return out


def install_skill(dest_root, name):
    """把模板复制成新技能：删掉 enabled: false、改 name（其余原样，模拟用户搬运）。"""
    dst = os.path.join(dest_root, name)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(os.path.join(ROOT, "skill", "_template"), dst)
    p = os.path.join(dst, "skill.yaml")
    with open(p, encoding="utf-8") as f:
        txt = f.read()
    txt = re.sub(r"^enabled:\s*false.*$", "enabled: true", txt, flags=re.M)
    txt = re.sub(r"^name:\s*template.*$", "name: %s" % name, txt, flags=re.M)
    txt = re.sub(r"^desc:\s*.*$", "desc: drop-in 自测技能 %s" % name, txt, flags=re.M)
    with open(p, "w", encoding="utf-8") as f:
        f.write(txt)
    return dst


def run(args, timeout=300):
    env = dict(os.environ)
    env["PATH"] = SHIM + os.pathsep + env.get("PATH", "")
    p = subprocess.run([sys.executable, PROG, "-c", CFG] + [str(x) for x in args],
                       capture_output=True, text=True, timeout=timeout, env=env)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def find_named(name, root=BASE):
    hits = []
    for dp, _dn, fns in os.walk(root):
        if name in fns:
            hits.append(os.path.join(dp, name))
        if os.path.basename(dp) == name:
            hits.append(dp)
    return hits


def setup():
    shutil.rmtree(BASE, ignore_errors=True)
    for d in (os.path.join(BASE, "skill"), os.path.join(BASE, "skill2"), SHIM):
        os.makedirs(d, exist_ok=True)
    mat = os.path.join(PROJ, "Si_new")
    os.makedirs(mat, exist_ok=True)
    shutil.copy(os.path.join(ROOT, "test", "tf_test", "Si", "POSCAR"),
                os.path.join(mat, "POSCAR"))
    install_skill(os.path.join(BASE, "skill"), "newsill")
    install_skill(os.path.join(BASE, "skill2"), "newsill2")
    # 本机只有 python3；tf 的 gen 用 python 跑脚本——做个壳（和 tmp/test_v1_prov.py 一样）
    shim = os.path.join(SHIM, "python")
    with open(shim, "w", encoding="utf-8") as f:
        f.write("#!/bin/sh\nexec python3 \"$@\"\n")
    os.chmod(shim, 0o755)
    with open(CFG, "w", encoding="utf-8") as f:
        f.write("host: \"\"\n"
                "project_roots:\n  - %s\n"
                "skill_paths:\n  - %s\n"
                "auto_advance: false\nauto_watch: false\n"
                "task_types:\n"
                "  newsill:\n    hpc: local\n    work_dir: %s\n"
                "  newsill2:\n    hpc: local\n    work_dir: %s\n"
                % (PROJ, os.path.join(BASE, "skill2"),
                   os.path.join(BASE, "remote"), os.path.join(BASE, "remote2")))


print("[0] 准备：核心代码指纹 + 两个投放点")
H0 = core_hashes()
setup()
ok(len(H0) >= 15, "核心指纹已取（%d 个文件）" % len(H0))
ok(os.path.isfile(os.path.join(BASE, "skill", "newsill", "skill.yaml")),
   "投放点①：配置目录旁的 skill/newsill/")
ok(os.path.isfile(os.path.join(BASE, "skill2", "newsill2", "skill.yaml")),
   "投放点②：tf.yaml 的 skill_paths 指向的 skill2/newsill2/")

print("[1] 自动发现：autozt skills 认得新技能（无需改任何配置）")
rc, out = run(["skills"])
ok(rc == 0, "autozt skills 正常退出")
ok("newsill" in out and "newsill2" in out, "两个新技能都在列表里")
line = next((x for x in out.splitlines() if x.strip().startswith("newsill ")), "")
ok("0.1" in line and re.search(r"\b1\b", line) is not None,
   "版本/步骤数正确：%s" % line.strip()[:70])
ok("模板本体" not in out, "模板本体（enabled: false）没被当成技能")

print("[2] tf skill show：新技能渲染出卡片（零代码改动）")
rc, out = run(["skill", "show", "newsill"])
ok(rc == 0, "渲染正常退出（rc=%d）" % rc)
ok("drop-in 自测技能 newsill" in out, "卡片抬头含技能描述")
ok("S1_summary" in out and "gen_step1_summary.py" in out,
   "卡片里有步骤与生成器（P0-3 声明式 I/O 对新技能同样生效）")

print("[3] tf schema --strict：自描述三段合法")
rc, out = run(["schema", "newsill", "--strict"])
ok(rc == 0, "newsill schema 无 [错误]（rc=%d）" % rc)
rc, out = run(["schema", "--strict"])
ok(rc == 0, "全库 schema --strict 仍然 0 错误")

print("[4] 新材料 + 项目配置（本地模式：ssh_host 置空）")
rc, out = run(["-tt", "newsill", "-p", "Si_new", "init"])
ok(rc == 0, "tf -tt newsill -p Si_new init 成功")
ps = os.path.join(PROJ, "Si_new", "newsill", "project_setting")
hpc = os.path.join(ps, "hpc.yaml")
ok(os.path.isfile(hpc), "生成了 project_setting/hpc.yaml")
if os.path.isfile(hpc):                     # 本地模式：ssh_host 置空 + work_dir 落本地
    with open(hpc, encoding="utf-8") as f:
        t = f.read()
    t = re.sub(r"^ssh_host:.*$", "ssh_host: \"\"", t, count=1, flags=re.M)
    with open(hpc, "w", encoding="utf-8") as f:
        f.write(t)
    with open(os.path.join(ps, "setting.yaml"), "a", encoding="utf-8") as f:
        f.write("\nwork_dir: \"%s\"\n" % os.path.join(BASE, "remote"))

print("[5] 新技能的 step1 真跑（gen → done_marker）")
rc, out = run(["-tt", "newsill", "-p", "Si_new", "-j", "S1_summary", "init"])
ok(rc == 0, "gen 成功（rc=%d）" % rc)
found = find_named("template_summary.json")
ok(bool(found), "done_marker 落地：%s"
   % (os.path.relpath(found[0], BASE) if found else "-"))
if found:
    with open(found[0], encoding="utf-8") as f:
        d = json.load(f)
    ok(d.get("step") == "step1_summary", "产物内容正确（step=step1_summary）")

print("[6] 新技能自动获得可追溯（provenance）")
prov = find_named("step1_summary.json", os.path.join(BASE, "remote"))
prov = [p for p in find_named("S1_summary.json") + find_named("step1_summary.json")
        if "provenance" in p or os.path.basename(p) == "step1_summary.json"]
ok(bool([p for p in prov if "provenance" in p]) or bool(prov),
   "provenance 由 tf 自动写入：%s"
   % (os.path.relpath(prov[0], BASE) if prov else "-"))
rc, out = run(["-tt", "newsill", "-p", "Si_new", "list"])
ok(rc == 0 and "Si_new" in out, "autozt list 能列出新材料")

print("[7] 核心代码零改动（技能是数据，不是代码）")
H1 = core_hashes()
diff = [k for k in H0 if H0[k] != H1.get(k)]
ok(not diff, "autozt/*.py 与 bin/autozt 的 sha256 一个没变")
ok(set(H0) == set(H1), "没有新增/删除核心文件")

print("")
if FAILS:
    print("FAILED %d 项：" % len(FAILS))
    for m in FAILS:
        print("  - %s" % m)
    sys.exit(1)
print("ALL PASS（v1.0 新技能 drop-in 自测）")
