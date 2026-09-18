# -*- coding: utf-8 -*-
"""自测：kl-dft-cpu S6_kappa 的 ShengBTE 完成判据（P49）不能用 small-grain limit 冒充。

为什么要有这条闸（2026-09-18 定）：
    BTE.KappaTensorVsT_sg 是 ShengBTE 的 **small-grain limit**（源码 ShengBTE.f90 里由
    kappasg() 写出，程序自印 "kappa in the small-grain limit"）—— 纯谐性、完全不含
    三声子散射的最小 κ，而且是在 calculate_Vp 之前就写完的。曾经的实现把它也算作
    "完成"，于是"三声子段崩了、只剩它"的作业被判 KAPPA_DONE=true，等于把最小 κ 当
    热导率交付（Mg4C60 7x7x7 实测 100 K：_sg=0.0782 vs _RTA=0.4740 W/m/K，差 6.06 倍）。
    本套件直接抽 submit_shengbte.tpl 里的汇总脚本，用假 BTE 产物跑三种情形，锁住：
      · 只有 _sg            → KAPPA_DONE=false，且写出 small_grain_limit_present 诊断
      · 有 _RTA（含 _sg）    → KAPPA_DONE=true，source=_RTA（不被 _sg 顶替）
      · _CONV 发散被过滤 <N → 回退 _RTA，仍判完成

跑法（仓库根）：python3 tests/suite_s6_marker.py；全绿 exit 0。
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TPL = Path(ROOT) / "skill/kl-dft-cpu/templates/step6_kappa/submit_shengbte.tpl"
FAILS = []


def ok(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("  x %s" % msg)
    else:
        print("  . %s" % msg)


def kt_rows(n=8, base=1.0):
    """ShengBTE KappaTensorVsT 行：col0=T，col1..9 张量（xx/yy/zz 有效即被接受）。"""
    out = []
    for i in range(n):
        out.append(" ".join("%.6f" % v for v in [100.0 + 100.0 * i] + [base] * 9))
    return "\n".join(out) + "\n"


def extract_summary_script():
    tpl = TPL.read_text(encoding="utf-8")
    m = re.search(r"python - <<'PY'\n(.*?)\nPY\n", tpl, re.S)
    if not m:
        sys.exit("模板里找不到汇总 heredoc（submit_shengbte.tpl 结构变了？）")
    src = m.group(1)
    src = src.replace("{{KAPPA_2D_FACTOR}}", "1.0")
    src = src.replace("{{KAPPA_2D_META}}", "{}")
    src = src.replace("{{KAPPA_2D_THICK2D}}", "{}")
    p = Path(ROOT) / "tmp" / "_s6_summary_extracted.py"
    p.parent.mkdir(exist_ok=True)
    p.write_text(src, encoding="utf-8")
    return p


def run_case(script, files):
    d = Path(tempfile.mkdtemp(prefix="s6marker_"))
    for fn, txt in files.items():
        (d / fn).write_text(txt, encoding="utf-8")
    p = subprocess.run([sys.executable, str(script)], cwd=str(d), capture_output=True, text=True)
    summ = {}
    if (d / "kappa_summary.json").is_file():
        summ = json.loads((d / "kappa_summary.json").read_text(encoding="utf-8"))
    shutil.rmtree(d, ignore_errors=True)
    return p, summ


print("[1] 抽取 submit_shengbte.tpl 的汇总脚本")
script = extract_summary_script()
ok(script.is_file(), "汇总脚本抽取成功：%s" % os.path.relpath(script, ROOT))

print("[2] P49：cand 不认 small-grain limit")
src = Path(TPL).read_text(encoding="utf-8")
ok("BTE.KappaTensorVsT_sg" not in re.search(r"(?m)^cand = \[.*\]$", src).group(0),
   "模板 cand 行不含 _sg：%s" % re.search(r"(?m)^cand = \[.*\]$", src).group(0).strip())

print("[3] 三种情形")
p, s = run_case(script, {"BTE.KappaTensorVsT_sg": kt_rows(8, 0.0782)})
ok(p.returncode == 0, "只有 _sg：脚本正常退出（rc=%d）" % p.returncode)
ok(s.get("KAPPA_DONE") is False, "只有 _sg：KAPPA_DONE=false（不再冒充完成）")
ok(s.get("small_grain_limit_present") is True, "只有 _sg：写出 small_grain_limit_present 诊断")
ok(s.get("small_grain_limit_rows") == 8, "只有 _sg：诊断记录行数 = %s" % s.get("small_grain_limit_rows"))
ok(bool(s.get("warning")), "只有 _sg：warning 说明不可当 κ 引用")
ok("KAPPA_DONE" not in p.stdout.split(), "只有 _sg：stdout 不打 KAPPA_DONE（marker 不出现）")

p, s = run_case(script, {"BTE.KappaTensorVsT_sg": kt_rows(8, 0.0782),
                         "BTE.KappaTensorVsT_RTA": kt_rows(8, 0.4740)})
ok(s.get("KAPPA_DONE") is True, "有 _RTA：KAPPA_DONE=true")
ok(str(s.get("source", "")).endswith("_RTA"), "有 _RTA：source=%s（不被 _sg 顶替）" % s.get("source"))
ok(s.get("small_grain_limit_present") is None, "有 _RTA：不写 small-grain 诊断")

p, s = run_case(script, {"BTE.KappaTensorVsT_CONV": kt_rows(2, 1e147),
                         "BTE.KappaTensorVsT_RTA": kt_rows(8, 0.4740),
                         "BTE.KappaTensorVsT_sg": kt_rows(8, 0.0782)})
ok(s.get("KAPPA_DONE") is True and str(s.get("source", "")).endswith("_RTA"),
   "_CONV 发散被过滤：回退 _RTA 判完成（source=%s）" % s.get("source"))

print("[4] gen 侧重写的那一行")
gen = (Path(ROOT) / "skill/kl-dft-cpu/gen_step6_kappa.py").read_text(encoding="utf-8")
ok("BTE.KappaTensorVsT_sg" not in gen.split("_cand = ")[-1].split("\n")[0],
   "gen_step6_kappa.py 的 _cand 行不含 _sg（渲染出来的 submit.sh 也安全）")

if FAILS:
    print("\nFAIL %d 项：" % len(FAILS))
    for m in FAILS:
        print("  - %s" % m)
    sys.exit(1)
print("\n全部通过。")
