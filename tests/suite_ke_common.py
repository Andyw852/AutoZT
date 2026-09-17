# -*- coding: utf-8 -*-
"""自测（A 批）：ke_common.py 单份化 + gen_need 声明完整性。

**这个套件只依赖技能清单与技能目录本身**，不依赖 autozt 的未提交改动，
所以它所属的那次提交（ke_common 单份化）能独立跑通它。
（跨步骤 src 兜底、step.conf strict 那两条依赖 autozt/_common 的改动，放在
 suite_asset_lookup.py 里，与那批一起提交。）

背景（2026-09-16 MoS2 实测）：ke 里曾有 3 份同名 ke_common.py，内容不一致，
远端只留最后 gen 的那份 → step7b 报 AttributeError: resolve_strain_pairs；
另外 S4/S8/S8.4 都 import ke_common 却没在 gen_need 里声明，靠上一步推过的
残留文件才跑得通。
"""
import glob
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from autozt import find_asset          # noqa: E402

FAILS = []


def ok(cond, msg):
    print(("  ✓ " if cond else "  ✗ ") + msg)
    if not cond:
        FAILS.append(msg)


def fake(t_base, sld, sname, fname, steps=None):
    t = {"key": "ke-dft-cpu", "skill_dir": sld, "template_dir": ".",
         "template_layout": "per_step", "_base_dir": t_base,
         "steps_cfg": steps or []}
    m = {"ps": None, "template_map": {}, "_skill_dir_local": None,
         "hpc_name": None, "_seg": {}}
    return find_asset({"_config_dir": None}, t, m, fname, sname)


def test_single_copy():
    print("[1] ke_common.py 只能有一份，且在技能根目录")
    base = os.path.join(ROOT, "skill")
    sk = os.path.join(base, "ke-dft-cpu")
    hits = glob.glob(os.path.join(sk, "**", "ke_common.py"), recursive=True)
    ok(len(hits) == 1 and os.path.dirname(hits[0]) == sk,
       "只有技能根目录一份（实际 %d 份）" % len(hits))
    if hits:
        p = fake(base, "ke-dft-cpu", "step7b_deform_read", "ke_common.py")
        ok(p is not None and "resolve_strain_pairs" in open(p).read(),
           "能按文件名解析到，且含 resolve_strain_pairs")


def test_gen_need_declared():
    print("[2] 凡是 import ke_common 的 gen 脚本，步骤必须声明它")
    import yaml
    sk = os.path.join(ROOT, "skill", "ke-dft-cpu")
    man = yaml.safe_load(open(os.path.join(sk, "skill.yaml")))
    steps = list(man.get("steps") or [])
    for grp in (man.get("optional_steps") or {}).values():
        steps += (grp or {}).get("steps") or []
    need = {}
    for s in steps:
        need[str(s.get("gen") or "").strip()] = [str(x) for x in (s.get("gen_need") or [])]
    users = []
    for root, _d, files in os.walk(sk):
        for fn in files:
            if fn.endswith(".py"):
                if "import ke_common" in open(os.path.join(root, fn), errors="ignore").read():
                    users.append(os.path.relpath(os.path.join(root, fn), sk))
    ok(bool(users), "找到 %d 个 import ke_common 的脚本" % len(users))
    missing = []
    for u in sorted(users):
        b = os.path.basename(u)
        declared = need.get(b) or need.get(u)
        if declared is None:
            missing.append("%s（skill.yaml 里没有这条 gen）" % u)
        elif "ke_common.py" not in declared:
            missing.append("%s 的 gen_need 没写 ke_common.py" % b)
    ok(not missing, "全部已声明" + ("" if not missing else "；缺：" + "; ".join(missing)))


def test_runtime_resolution():
    """运行期冒烟：对**每个** import ke_common 的步骤真的调用一次 find_asset，
    确认解析得到技能根目录那份 —— 只查"声明了没有"查不出"声明了但解析不到"。"""
    print("[3] 运行期解析：每个步骤 find_asset('ke_common.py') 都要命中根目录那份")
    import yaml
    sk = os.path.join(ROOT, "skill", "ke-dft-cpu")
    man = yaml.safe_load(open(os.path.join(sk, "skill.yaml")))
    steps = list(man.get("steps") or [])
    for grp in (man.get("optional_steps") or {}).values():
        steps += (grp or {}).get("steps") or []
    steps_cfg = [{"name": s.get("name"), "src": s.get("src"), "gen": s.get("gen")}
                 for s in steps]
    t = {"key": "ke-dft-cpu", "skill_dir": "ke-dft-cpu", "template_dir": ".",
         "template_layout": "per_step", "_base_dir": os.path.join(ROOT, "skill"),
         "steps_cfg": steps_cfg}
    m = {"ps": None, "template_map": {}, "_skill_dir_local": None,
         "hpc_name": None, "_seg": {}}
    exp = os.path.join(sk, "ke_common.py")
    users = []
    for s in steps_cfg:
        g = str(s.get("gen") or "").strip()
        if not g:
            continue
        cands = [os.path.join(sk, str(s.get("src") or ""), g), os.path.join(sk, g)]
        p = next((c for c in cands if os.path.isfile(c)), None)
        if p and "import ke_common" in open(p, errors="ignore").read():
            users.append((s["name"], g))
    ok(bool(users), "找到 %d 个需要解析 ke_common.py 的步骤" % len(users))
    bad = []
    for name, g in users:
        got = find_asset({"_config_dir": None}, t, m, "ke_common.py", name)
        if got is None or os.path.realpath(got) != os.path.realpath(exp):
            bad.append("%s(%s)->%s" % (name, g, got))
    ok(not bad, "全部命中 skill/ke-dft-cpu/ke_common.py"
       + ("" if not bad else "；未命中：" + "; ".join(bad)))


def main():
    print("== suite_ke_common：ke_common 单份化 / gen_need 声明 ==")
    test_single_copy()
    test_gen_need_declared()
    test_runtime_resolution()
    print()
    if FAILS:
        print("FAIL %d 项：" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
