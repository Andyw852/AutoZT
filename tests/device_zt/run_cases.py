#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""device-zt 回归用例运行器 —— 自包含，不联网、不连超算、不提交 SLURM。

用法:
    python3 tests/device_zt/run_cases.py [-v] [用例名 ...]

对 cases/ 下每个用例:
  1. 把 upstream/ 摊成"材料目录" <tmp>/<case>/<技能>/<步骤>/<文件> —— 这正是
     find_upstream 在超算上看到的第一候选布局 <matdir>/<技能>/<步骤>/<文件>;
  2. 用 autozt 自己的 report.find_asset 解析 device-zt 的 gen_need 并平铺进
     <tmp>/<case>/device-zt/（顺带验证 skill/_common/<组>/ 公共池解析）;
  3. 按 case.json 的 params 覆盖技能模板 step.conf;
  4. 清空 PYTHONPATH，在 <tmp>/<case>/device-zt/ 原地跑三个 gen 脚本
     （与 autozt 在超算上的运行方式一致）;
  5. 读产物，与 case.json 的 expected 逐项比对（数值用相对容差）。
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
import yaml                             # noqa: E402
from autozt import report               # noqa: E402

SKILL = "device-zt"
TOL = 1e-6
FAILS = []


def build(tmpmat, case):
    """按 autozt 的推送规则，把技能跑起来所需的文件平铺进 <tmpmat>/device-zt/。"""
    t = yaml.safe_load(open(os.path.join(ROOT, "skill", SKILL, "skill.yaml")))
    t["skill_dir"] = "skill/" + SKILL
    cfg = {"_config_dir": ROOT}
    m = {"hpc_name": None, "_skill_dir_local": None, "_seg": {}}
    dest = os.path.join(tmpmat, SKILL)
    os.makedirs(dest, exist_ok=True)
    for s in t["steps"]:
        need = list(s.get("gen_need") or [])
        gs = s["gen"].split()[0]
        if gs not in need:
            need.append(gs)
        for f in need:
            src = report.find_asset(cfg, t, m, f, s["name"])
            if not src:
                raise RuntimeError("find_asset 解析不到 gen_need: %s (步骤 %s)" % (f, s["name"]))
            shutil.copy2(src, os.path.join(dest, f))
    # step.conf = 技能模板 + 用例参数覆盖
    conf = os.path.join(dest, "step.conf")
    lines = open(conf, encoding="utf-8").read().split("\n")
    for k, v in (case.get("params") or {}).items():
        hit, out = False, []
        for ln in lines:
            st = ln.strip()
            if st.startswith(k + " ") or st.startswith("# " + k + " "):
                indent = ln[:len(ln) - len(ln.lstrip())]
                out.append(indent + "%s = %s" % (k, v))
                hit = True
            else:
                out.append(ln)
        if not hit:
            out.insert(out.index("[params]") + 1, "%s = %s" % (k, v))
        lines = out
    open(conf, "w", encoding="utf-8").write("\n".join(lines))
    return t


def run(tmpmat):
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    logs = []
    for gs in ("gen_step1_elec.py", "gen_step2_selftrap.py", "gen_step3_report.py"):
        p = subprocess.run([sys.executable, gs], cwd=os.path.join(tmpmat, SKILL),
                           env=env, capture_output=True, text=True)
        logs.append("--- %s exit=%d ---\n%s%s" % (gs, p.returncode, p.stdout, p.stderr))
        if p.returncode != 0:
            return False, "\n".join(logs)
    return True, "\n".join(logs)


def cmp(exp, got, path, verbose):
    if isinstance(exp, dict):
        for k, v in exp.items():
            if k not in got:
                FAILS.append("%s.%s 缺失" % (path, k))
            else:
                cmp(v, got[k], path + "." + k, verbose)
    elif isinstance(exp, list):
        if not isinstance(got, list) or len(got) != len(exp):
            FAILS.append("%s 长度不符 %s vs %s" % (path, len(exp), len(got)))
        else:
            for i, v in enumerate(exp):
                cmp(v, got[i], "%s[%d]" % (path, i), verbose)
    elif isinstance(exp, bool) or exp is None:
        if exp != got:
            FAILS.append("%s 期望 %r 实得 %r" % (path, exp, got))
    elif isinstance(exp, (int, float)):
        rel = abs(got - exp) / max(abs(exp), 1e-30)
        if rel > TOL:
            FAILS.append("%s 期望 %.10g 实得 %.10g 相对差 %.2e" % (path, exp, got, rel))
        elif verbose:
            print("      ok %s = %.10g (差 %.1e)" % (path, got, rel))
    else:
        if exp != got:
            FAILS.append("%s 期望 %r 实得 %r" % (path, exp, got))


def main(argv):
    verbose = "-v" in argv
    want = [a for a in argv if not a.startswith("-")]
    cases_dir = os.path.join(HERE, "cases")
    names = sorted(d for d in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir, d)))
    if want:
        names = [n for n in names if n in want]
    rc = 0
    for name in names:
        case = json.load(open(os.path.join(cases_dir, name, "case.json"), encoding="utf-8"))
        print("[%s] %s (%s, %s)" % (name, case.get("material"), case.get("dim"), "实例"))
        del FAILS[:]
        tmp = tempfile.mkdtemp(prefix="device_zt_case_")
        try:
            tmpmat = os.path.join(tmp, name)
            shutil.copytree(os.path.join(cases_dir, name, "upstream"), tmpmat)
            build(tmpmat, case)
            ok, log = run(tmpmat)
            if not ok:
                print("  gen 失败:\n" + log)
                rc = 1
                continue
            got = {}
            for key, rel in (("S1", "step1_elec/elec_props.json"),
                             ("S2", "step2_selftrap/device_zt_solution.json"),
                             ("S3", "step3_report/device_zt_report.json")):
                got[key] = json.load(open(os.path.join(tmpmat, SKILL, rel), encoding="utf-8"))
            exp = case["expected"]
            cmp(exp["S1"], got["S1"], "S1", verbose)
            cmp(exp["S2"], got["S2"], "S2", verbose)
            best = max(got["S3"]["zt"]["doping_scan"]["rows"], key=lambda x: x["zt_estimate"])
            cmp(exp["S3"]["zt_material"], got["S3"]["zt"]["zt_material"], "S3.zt_material", verbose)
            cmp(exp["S3"]["electronic_share"], got["S3"]["thermal"]["electronic_share"],
                "S3.electronic_share", verbose)
            cmp(exp["S3"]["best_doping_cm3"], best["doping_cm3"], "S3.best_doping_cm3", verbose)
            cmp(exp["S3"]["best_zt"], best["zt_estimate"], "S3.best_zt", verbose)
            if FAILS:
                print("  FAIL (%d 项):" % len(FAILS))
                for f in FAILS[:15]:
                    print("    -", f)
                rc = 1
            else:
                print("  PASS  dT=%.4f K  ZT=%.6g  V_max_safe=%.2f V"
                      % (got["S2"]["dT_peak_K"], got["S2"]["zt_material"], got["S2"]["V_max_safe_V"]))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    print("\n%s" % ("全部通过" if rc == 0 else "有用例失败"))
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
