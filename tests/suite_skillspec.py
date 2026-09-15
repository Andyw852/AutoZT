# -*- coding: utf-8 -*-
"""v1.0 自测：skillspec 校验 + _corrections/ handler 匹配（纯本地，不连超算）。

跑法（在 taskflow-v1.0 仓库根）：
    python3 tmp/test_v1_skillspec.py
全绿 = exit 0；任一断言失败 = exit 1（最后打印 FAIL 清单）。
"""
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import autozt  # noqa: E402
from autozt import (load_correction_handlers, correction_handler_names,  # noqa: E402
                   match_corrections, make_ctx, suggest_for_diag,
                   validate_skill_spec, spec_stats, issues_fatal)

FAILS = []


def ok(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("  ✗ %s" % msg)
    else:
        print("  ✓ %s" % msg)


CFG = {"_config_dir": ROOT, "skill_paths": [], "_skills": {}}

print("[1] 纠错库加载")
names = correction_handler_names(CFG)
ok(set(names) >= {"zbrent", "nsw_continue", "node_fail", "dir_missing"},
   "内置 4 个 handler 都注册了：%s" % names)


def hit(diag, skill="band-dft-cpu", label="S1_opt", step="step1_PBE_opt"):
    ctx = make_ctx(CFG, material="C24/qHPC24", skill=skill, label=label,
                   step=step, diag=diag)
    return [h.name for h, _s in match_corrections(CFG, ctx, skill=skill)]


print("[2] 诊断 → handler 匹配")
ok(hit("force not converged") == ["nsw_continue"], "NSW 打满 → nsw_continue")
ok(hit("ZBRENT: fatal error in subroutine") == ["zbrent"], "ZBRENT → zbrent")
ok(hit("NODE_FAIL") == ["node_fail"], "NODE_FAIL → node_fail")
ok(hit("dir missing") == ["dir_missing"], "目录缺失 → dir_missing")
ok(hit("unmapped weird thing") == [], "陌生诊断 → 不硬套 handler")

print("[3] 步骤适用性（steps 前缀过滤）")
ok(hit("force not converged", label="S9_zzz", step="step9_zzz") == [],
   "nsw_continue 不套到无关步骤（S9）")
ok(hit("ZBRENT", label="S9_zzz", step="step9_zzz") == ["zbrent"],
   "zbrent 是全步骤 handler，任何步骤都命中")

print("[4] 建议内容与风险等级")
sug = suggest_for_diag(CFG, material="C24/qHPC24", skill="band-dft-cpu",
                       step="step1_PBE_opt", label="S1_opt", diag="dir missing")
ok(sug and sug[0]["risk"] == "destructive", "dir_missing 标成 destructive")
ok(sug and sug[0]["apply_available"] is False,
   "destructive handler 不给 apply（只许人工 rerun）")
ok(sug and any("rerun" in c for c in sug[0]["commands"]), "建议里给了 rerun 命令")
sug2 = suggest_for_diag(CFG, material="C24/qHPC24", skill="band-dft-cpu",
                        step="step1_PBE_opt", label="S1_opt",
                        diag="ZBRENT: fatal error")
ok(sug2 and sug2[0]["apply_available"] is True, "zbrent 提供可执行的 apply")
ok(sug2 and any("retry" in c for c in sug2[0]["commands"]), "zbrent 建议里给 retry 命令")

print("[5] 技能私有纠错目录（_corrections/ 隔离）")
tmpd = tempfile.mkdtemp(prefix="tfv1test_")
try:
    fake_skill = os.path.join(tmpd, "fake")
    os.makedirs(os.path.join(fake_skill, "_corrections"))
    shutil.copy(os.path.join(ROOT, "skill", "_common", "_corrections", "base.py"),
                os.path.join(fake_skill, "_corrections", "base.py"))
    with open(os.path.join(fake_skill, "_corrections", "only_fake.py"), "w",
              encoding="utf-8") as f:
        f.write("from .base import CorrectionHandler, Suggestion\n"
                "class OnlyFake(CorrectionHandler):\n"
                "    name = 'only_fake'\n"
                "    title = 'Fake 技能专属纠错'\n"
                "    matches = ('zbrent',)\n"
                "    steps = ('*',)\n"
                "    risk = 'safe'\n"
                "    def suggest(self, ctx):\n"
                "        return Suggestion(handler=self.name, title=self.title,\n"
                "                          risk=self.risk, reason='fake')\n"
                "HANDLER = OnlyFake()\n")
    cfg2 = {"_config_dir": ROOT, "skill_paths": [],
            "_skills": {"fake": {"skill_dir": fake_skill}}}
    nm2 = correction_handler_names(cfg2)
    ok("only_fake" in nm2, "技能私有 _corrections/ 被加载：%s" % nm2)
    ctx_fake = make_ctx(cfg2, material="X", skill="fake", label="S1_opt", diag="ZBRENT")
    ok("only_fake" in [h.name for h, _s in match_corrections(cfg2, ctx_fake, skill="fake")],
       "fake 技能命中自己的 handler")
    ctx_band = make_ctx(cfg2, material="X", skill="band-dft-cpu", label="S1_opt", diag="ZBRENT")
    ok("only_fake" not in [h.name for h, _s in match_corrections(cfg2, ctx_band,
                                                                skill="band-dft-cpu")],
       "fake 的私有 handler 不误伤 band-dft-cpu")
finally:
    shutil.rmtree(tmpd, ignore_errors=True)

print("[6] skillspec 校验：合法清单")
good = {
    "schema": 2,
    "io_schema": {
        "inputs": [{"name": "POSCAR", "from": "user", "required": True}],
        "outputs": [{"name": "CONTCAR", "step": "S1_opt"}],
        "params": [{"name": "MODE", "values": ["a", "b"], "default": "a"}],
    },
    "flow": {"summary": "x", "stages": [{"name": "s", "steps": ["S1_opt"]}],
             "next_skills": [{"skill": "band-dft-cpu"}]},
    "corrections": ["zbrent", {"name": "node_fail", "desc": "节点故障"}],
}
skel_good = {"steps": [{"name": "step1_opt", "label": "S1_opt", "seq": 1}]}
iss = validate_skill_spec("demo", good, skel=skel_good,
                          known_skills={"demo", "band-dft-cpu"},
                          handler_names=set(names))
ok(iss == [], "合法清单 0 问题（实际：%s）" % iss)
ok(spec_stats(good)["n_in"] == 1 and spec_stats(good)["n_corr"] == 2, "统计正确")

print("[7] skillspec 校验：故意写错")
bad = {
    "schema": 1,                       # 用了扩展段却还是 schema 1 → 警告
    "io_schema": {
        "inputs": [{"name": "", "oops": 1}],          # 缺 name(错误) + 未知键(警告)
        "outputs": [{"name": "X", "step": "S99_nope"}],  # 步骤不存在 → 警告
        "params": [{"name": "P", "values": "a|b"},              # values 不是列表 → 警告
                   {"name": "P2", "values": ["a", "b"], "default": "z"}],  # 默认值不在集合 → 警告
    },
    "flow": {"bad_key": 1, "stages": [{"steps": ["S1_opt"]}]},   # 未知键 + 缺 name
    "corrections": [{"desc": "缺 name"}, "not_registered"],       # 缺 name + 未注册
}
iss2 = validate_skill_spec("demo", bad, skel=skel_good,
                           known_skills={"demo"}, handler_names=set(names))
txt = "\n".join(iss2)
ok(issues_fatal(iss2), "--strict 会失败（存在 [错误]）")
for want in ("缺少 name", "不认识的键", "S99_nope 在本技能步骤里找不到",
             "values 建议写成列表", "default=z 不在 values 里",
             "schema=1", "在 _corrections/ 库里没找到"):
    ok(want in txt, "抓到：%s" % want)
ok(spec_stats({})["has_io"] is False, "空 spec 统计不炸")


print("[9] io_schema.steps（按步 I/O）+ tf skill show 卡片渲染")
from autozt.skillspec import (_validate_io_steps, step_io_map, render_skill_card,
                             card_dict, _disp_width, _fit)  # noqa: E402

skel_card = {
    "desc": "演示技能", "_skill_version": "1.0", "steps": [
        {"seq": 1, "name": "step1_a", "label": "S1_a", "gen": "gen_step1_a.py",
         "check": "outcar"},
        {"seq": 2, "name": "step2_b", "label": "S2_b", "gen": "gen_step2_b.py",
         "check": "marker", "marker": 'x.json:"DONE": true'},
    ],
    "optional_steps": {"extra": {"default": False, "steps": [
        {"seq": 3, "name": "step3_c", "label": "S3_c", "gen": "plot.py",
         "check": "plot", "run": "gen", "done_marker": "plot.json"}]}},
}
spec_card = {
    "schema": 2,
    "io_schema": {"inputs": [{"name": "POSCAR"}], "params": [{"name": "P1"}],
                  "outputs": [{"name": "out.json"}],
                  "steps": [
                      {"step": "S1_a", "title": "Relax structure",
                       "tool": "VASP (vasp_std)", "outputs": ["CONTCAR"],
                       "validator": "力收敛"},
                      {"step": "step2_b", "title": "Static SCF",
                       "tool": "VASP", "outputs": ["OUTCAR"]},
                  ]},
    "flow": {"summary": "演示流程", "next_skills": [{"skill": "band-dft-cpu", "via": "CONTCAR"}]},
    "corrections": ["zbrent"],
}
iss = validate_skill_spec("demo", spec_card, skel=skel_card,
                          known_skills={"demo", "band-dft-cpu"},
                          handler_names=set(names))
ok(not iss, "合法的按步 I/O 不报问题：%s" % iss)
m_io = step_io_map(spec_card)
ok(set(m_io) == {"S1_a", "step2_b"}, "step_io_map 用 name/label 都能查到")
bad_io = _validate_io_steps([{"title": "缺 step"}, {"step": "S99_x", "oops": 1},
                             {"step": "S1_a", "outputs": "不是列表"}], {"S1_a", "S2_b"})
bt = "\n".join(bad_io)
ok("缺少 step" in bt, "抓到：缺 step")
ok("S99_x] 不是本技能的步骤" in bt, "抓到：步骤引用不存在")
ok("不认识的键 'oops'" in bt, "抓到：未知键")
ok("outputs 建议写成列表" in bt, "抓到：类型不是列表")
ok(not _validate_io_steps(None, set()), "没写 io_schema.steps 不报问题")

# 图 2 版式：外框 + Skill:/tagline 抬头 + 一行 4 个盒子（Step N / 标题 /
# Generator / Tool / Validator / Output）+ 两条页脚
card = render_skill_card("demo", skel_card, spec_card)
lines = card.splitlines()
for want in ("Skill: 演示技能", "Fixed validated workflow; configurable templates",
             "Skill definition: skill.yaml + templates + generators + validators",
             "Execution backend: VASP", "Step 1", "Step 2", "Step 3",
             "Relax structure", "Static SCF", "Generator", "gen_step1_a.py",
             "Tool", "VASP (vasp_std)", "Validator", "Output", "CONTCAR"):
    ok(want in card, "卡片里有：%s" % want)
ok(lines[0].startswith("┌") and lines[-1].startswith("└"), "外框完整（首尾是 ┌…┐/└…┘）")
ok(all(l.startswith("│") and l.endswith("│") for l in lines[1:-1]),
   "框内每行都有左右边界")
ok("→" in card, "步骤之间用箭头连接")
ok("(可选)" in card, "可选步骤标了 (可选)")
ok(all(_disp_width(l) == _disp_width(lines[0]) for l in lines), "所有行等宽（排版对齐）")
ok(_disp_width(lines[0]) <= 110, "不超过设定宽度（%d 列）" % _disp_width(lines[0]))
_row = [l for l in lines if "Step 1" in l and "Step 2" in l]
ok(bool(_row), "前两步在同一行（一行多列，照图 2 版式）")

card_full = render_skill_card("demo", skel_card, spec_card, full=True)
for want in ("Inputs: POSCAR", "Parameters: P1", "Next skills: band-dft-cpu",
             "Corrections: zbrent"):
    ok(want in card_full, "--full 里有：%s" % want)
ok("Next skills" not in card, "默认（图 2 版）不带 --full 的额外行")

cd = card_dict("demo", skel_card, spec_card)
ok(cd["steps"][0]["generator"] == "gen_step1_a.py"
   and cd["steps"][0]["tool"] == "VASP (vasp_std)"
   and cd["steps"][2]["optional"] is True, "card_dict 结构稳定（含 optional 标记）")
ok(cd["title"] and cd["tagline"], "card_dict 带 title/tagline（画图脚本要用）")
ok(_fit("abcdef", 4).endswith("…") and _disp_width(_fit("中文中文", 5)) <= 5,
   "截断按显示宽度算（中文不被切坏）")
ok(bool(render_skill_card("demo", {"steps": []}, {})), "没步骤的技能也不炸")
_spec_t = dict(spec_card, flow={"title": "Lattice thermal conductivity",
                                "backend": "VASP / MACE / Pheasy"})
ok("Skill: Lattice thermal conductivity" in render_skill_card("demo", skel_card, _spec_t)
   and "Execution backend: VASP / MACE / Pheasy" in render_skill_card("demo", skel_card, _spec_t),
   "flow.title / flow.backend 生效（论文图要英文抬头与工具链）")

print("")
if FAILS:
    print("FAILED: %d 项" % len(FAILS))
    for f in FAILS:
        print("  - %s" % f)
    sys.exit(1)
print("ALL PASS")
