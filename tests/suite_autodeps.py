#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""公共模块依赖补全（common_autodeps）自测。

背景：2026-09-08 往 skill/_common/opt/ 加了 method_select.py 并让 relax_common.py
import 它，但 5 个技能的 gen_need 都没声明 -> 永远推不到超算 -> 任何新材料首次
gen 直接 ModuleNotFoundError（autozt retry 也救不回来）。

本脚本验证新的"公共模块自动带同目录依赖"能在【清单故意漏写】时把它补上，
并且：已声明时不重复、技能目录优先、开关可关、gen 能真跑通（本地 bash 模式）。
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
FAILED = []


def check(name, ok, detail=""):
    print("  %s %s%s" % ("✓" if ok else "✗", name, ("  ← " + detail) if detail else ""))
    if not ok:
        FAILED.append(name)


def section(t):
    print("\n[%s]" % t)


from autozt import load_config, find_asset, _common_dep_closure  # noqa: E402
from autozt.workflow import (common_autodeps_enabled, _py_sibling_imports,  # noqa: E402
                            remote_gen)

CFG = load_config("tmp/tf_jzzn_smoke.yaml")[0]

T = {"key": "kl-dft-cpu", "skill_dir": "skill/kl-dft-cpu",
     "steps": [{"name": "step1_std_opt"}], "template_layout": "per_step"}


def _m():
    return {"name": "Si_auto", "template_map": {}, "ps": {}, "_seg": {}}


section("1. 解析：能找到公共池模块，并看出它的同目录依赖")
relax = find_asset(CFG, T, _m(), "relax_common.py", "step1_std_opt")
check("relax_common.py 解析到公共池", bool(relax) and "_common" in relax, relax or "None")
sib = _py_sibling_imports(relax) if relax else {}
check("看出它 import 了 method_select", "method_select" in sib)
check("也看得到 dim_common/mol_common/stepconf",
      {"dim_common", "mol_common", "stepconf"} <= set(sib), ",".join(sorted(sib)))
check("也看得到 symmetry_audit（结构体检第五条，2026-09-20）",
      "symmetry_audit" in sib)
check("也看得到 structure_health（结构体检四条，2026-09-21）",
      "structure_health" in sib)

section("2. 补全：清单漏写时自动补上")
need = ["relax_common.py", "mol_common.py", "dim_common.py", "stepconf.py",
        "check_common.py"]
extra = _common_dep_closure(CFG, T, _m(), need, "step1_std_opt")
# 2026-09-20/21：relax_common 新增 import symmetry_audit（第五条）+ structure_health（四条），
# 自动补推的闭包随之多两项 —— 这正是本套件要覆盖的"公共池加文件别忘清单"。
check("漏写的 method_select.py/symmetry_audit.py/structure_health.py 被补上",
      set(extra) == {"method_select.py", "symmetry_audit.py", "structure_health.py"},
      str(extra))

section("3. 幂等与优先级")
check("清单已声明时不重复",
      _common_dep_closure(CFG, T, _m(),
                          need + ["method_select.py", "symmetry_audit.py",
                                  "structure_health.py"],
                          "step1_std_opt") == [])
only = _common_dep_closure(CFG, T, _m(), ["relax_common.py"], "step1_std_opt")
check("只写 relax_common 时补齐全部同目录依赖",
      set(only) == {"dim_common.py", "method_select.py", "mol_common.py",
                    "stepconf.py", "symmetry_audit.py", "structure_health.py"},
      str(sorted(only)))
check("技能目录里已有的同名文件不补（技能优先）",
      "check_common.py" not in only)

section("4. 开关")
check("默认开", common_autodeps_enabled(CFG) is True)
os.environ["AUTOZT_COMMON_AUTODEPS"] = "0"
check("AUTOZT_COMMON_AUTODEPS=0 可关", common_autodeps_enabled(CFG) is False)
del os.environ["AUTOZT_COMMON_AUTODEPS"]
check("tf.yaml common_autodeps: false 可关",
      common_autodeps_enabled({"common_autodeps": False}) is False)

section("5. 真跑一次 gen（本地 bash 模式，自造小技能树）：清单【故意漏写】依赖也能过")
# 造一棵自足的小树：gen 脚本 import helper.py，helper.py 又 import 同目录的
# method_select.py，但 gen_need 只写 helper.py —— 修复前这就是"新材料必崩"的翻版。
tmp = tempfile.mkdtemp(prefix="tf_autodeps_")
work = os.path.join(tmp, "remote", "X", "step1")
skill = os.path.join(tmp, "skill")
pool = os.path.join(skill, "_common", "opt")
demo = os.path.join(skill, "demo")
for d in (work, pool, demo):
    os.makedirs(d, exist_ok=True)


def _w(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


_w(os.path.join(pool, "method_select.py"),
   "def default_method(*a, **k):\n    return 'pbe', 'default'\n")
_w(os.path.join(pool, "helper.py"),
   "import method_select\n\ndef run():\n    return method_select.default_method()\n")
_w(os.path.join(demo, "gen_demo.py"),
   "#!/usr/bin/env python3\nimport helper\n"
   "open('DEMO_OK','w').write(str(helper.run()[0]))\n")
# remote_gen 用 `python` 跑 gen 脚本（超算上有）；本机只有 python3，做个壳
bindir = os.path.join(tmp, "bin")
os.makedirs(bindir, exist_ok=True)
_shim = os.path.join(bindir, "python")
_w(_shim, "#!/bin/sh\nexec python3 \"$@\"\n")
os.chmod(_shim, 0o755)
os.environ["PATH"] = bindir + os.pathsep + os.environ.get("PATH", "")
cfg2 = {"_config_dir": tmp, "host": None, "_skill_root": skill}
# skill_dir 用绝对路径指向"技能目录"（demo/），与 tf 装配后的值一致
t2 = {"key": "demo", "root": demo, "skill_dir": demo, "template_layout": "shared",
      "steps": [{"name": "step1", "label": "S1", "gen": "gen_demo.py",
                 # ★ 故意只写 helper.py —— method_select.py 漏写
                 "gen_need": ["helper.py"]}]}
m2 = {"name": "X", "path": work, "hpc_name": None, "host_eff": None}
ok, msg = remote_gen(cfg2, t2, m2, "step1", host=None)
check("清单漏写 method_select.py 时 gen 仍成功", bool(ok),
      "" if ok else str(msg)[:600].replace("\n", " | "))
check("method_select.py 被自动补推到生成目录",
      os.path.isfile(os.path.join(work, "method_select.py")))
marker = os.path.join(work, "DEMO_OK")
check("gen 脚本真的 import 到了它（产物 DEMO_OK=pbe）",
      os.path.isfile(marker) and open(marker, encoding="utf-8").read().strip() == "pbe",
      open(marker, encoding="utf-8").read().strip() if os.path.isfile(marker) else "无产物")
# 关掉开关应当【不】补推（证明开关真的生效）
os.environ["AUTOZT_COMMON_AUTODEPS"] = "0"
work2 = os.path.join(tmp, "remote2", "X", "step1")
os.makedirs(work2, exist_ok=True)
m3 = dict(m2, path=work2, lpath=None)
remote_gen(cfg2, t2, m3, "step1", host=None)
del os.environ["AUTOZT_COMMON_AUTODEPS"]
check("AUTOZT_COMMON_AUTODEPS=0 时不补推（开关有效）",
      not os.path.isfile(os.path.join(work2, "method_select.py")))
shutil.rmtree(tmp, ignore_errors=True)

section("6. ZA 2D 共享模块（三副本合并）")
# 低层 ZA 判定收敛到 _common/za_2d.py；kl 的 S5_fc、fc-fit 的 S2_plot、以及全部 6 个
# mlff 技能（kl/opt/phonon × cpu/gpu，2026-09-22 第三副本迁移）都要能解析到它，
# 否则首次上机就是 ModuleNotFoundError（公共池文件没推到生成目录）。
for _tt, _sd, _step in (("kl-dft-cpu", "skill/kl-dft-cpu", "step5_fc"),
                        ("fc-fit", "skill/fc-fit", "step2_phonon_plot"),
                        ("kl-mlff-cpu", "skill/kl-mlff-cpu", "step3_fc"),
                        ("kl-mlff-gpu", "skill/kl-mlff-gpu", "step3_fc"),
                        ("opt-mlff-cpu", "skill/opt-mlff-cpu", "step1_mlff_relax"),
                        ("opt-mlff-gpu", "skill/opt-mlff-gpu", "step1_mlff_relax"),
                        ("phonon-mlff-cpu", "skill/phonon-mlff-cpu", "step3_phonon"),
                        ("phonon-mlff-gpu", "skill/phonon-mlff-gpu", "step3_phonon")):
    _Tz = {"key": _tt, "skill_dir": _sd, "steps": [{"name": _step}],
           "template_layout": "per_step"}
    _az = find_asset(CFG, _Tz, _m(), "za_2d.py", _step)
    check("%s/%s 解析到 _common/za_2d.py" % (_tt, _step),
          bool(_az) and "_common" in _az, _az or "None")

print()
if FAILED:
    print("FAILED: %d 项 -> %s" % (len(FAILED), ", ".join(FAILED)))
    sys.exit(1)
print("ALL PASS")
