#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""StepConf dict 协议 + S4 数据集四重比对（_dataset_matches）回归。

为什么单独立套（2026-09-22）：
  ① StepConf 原来只实现 __getitem__，而下游"通用比对/留档"代码按 dict 用
     （gen_step4_disp.py 写 disp_plan 时 conf.get("ALM_CUT2")）—— S4 生成位移
     数据集直接 AttributeError，崩在写完 POSCAR-*/yaml 之后，留下没有留档的
     半成品数据集，幂等检查随后也认不出它（实测 MoS₂ 踩到）。
  ② _dataset_matches 用 _old is None 当"老 plan 没这个键"，把"记着不截断"和
     "根本没这个键"混为一谈 ⇒ ALM_CUT2 从空改成 4.0 这类变化不会被拦，
     "前后步骤参数对不上就不要运行"的闸形同虚设。
本套件把这两条钉成硬闸。
"""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skill" / "_common" / "opt"))
sys.path.insert(0, str(ROOT / "skill" / "kl-dft-cpu"))

import stepconf  # noqa: E402

FAIL = []


def ck(cond, msg):
    if cond:
        print("  PASS  %s" % msg)
    else:
        print("  FAIL  %s" % msg)
        FAIL.append(msg)


def mkconf(values):
    """按 StepConf 的真实入参形状造 conf（merged = parse() 的输出）。"""
    spec = {"OVERSAMPLE": (3, "int"), "ALM_CUT2": (None, "float"),
            "ALM_CUT3": (6.0, "float"), "FD_DISTANCE": (0.03, "float"),
            "FC3_CUTOFF_PAIR": (None, "float")}
    rows = [(k, ("" if values.get(k) is None else str(values[k])), i + 1)
            for i, k in enumerate(spec) if k in values]
    return stepconf.StepConf({"params": rows}, spec)


print("== A. StepConf 必须表现得像只读 dict ==")
c = mkconf({"OVERSAMPLE": 3, "ALM_CUT2": None, "ALM_CUT3": 6.0,
            "FD_DISTANCE": 0.03, "FC3_CUTOFF_PAIR": None})
ck(c["OVERSAMPLE"] == 3, "__getitem__ 可用")
ck(c.get("ALM_CUT2") is None, "get() 存在（原崩溃点 conf.get('ALM_CUT2')）")
ck(c.get("NOPE") is None and c.get("NOPE", 7) == 7, "get(未声明键) 落到 default")
ck("ALM_CUT3" in c and "NOPE" not in c, "__contains__ 正确")
ck(sorted(c.keys()) == sorted(c.params.keys()), "keys() 与 params 一致")
ck(dict(c.items()) == dict(c.params), "items() 与 params 一致")
ck(list(c.values()) == [c.params[k] for k in c.keys()], "values() 顺序与 keys() 对齐")
ck(len(c) == len(c.params), "__len__ 正确")
ck(set(iter(c)) == set(c.params), "__iter__ 可迭代")
ck(c.copy() == dict(c.params) and isinstance(c.copy(), dict), "copy() 返回真 dict")
ck(c.get("ALM_CUT3") == 6.0, "get() 返回已转好类型（float）")
try:
    c["NEW"] = 1
    ck(False, "不该能直接写入 StepConf")
except TypeError:
    ck(True, "写入仍被拒绝（只读）")


print("== B. 两个 StepConf 副本行为一致 ==")
import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location(
    "stepconf_mlff", str(ROOT / "skill" / "mlff" / "stepconf.py"))
_m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m)
c2 = _m.StepConf({"params": [("OVERSAMPLE", "3", 1)]}, {"OVERSAMPLE": (3, "int")})
ck(c2.get("OVERSAMPLE") == 3 and c2.get("NOPE", "d") == "d",
   "skill/mlff/stepconf.py 的 StepConf 同样有 get()")

print("== C. _dataset_matches 的四重比对（含 None 语义）==")
try:
    import json
    import gen_step4_disp as g
except Exception as e:                                      # noqa: BLE001
    print("  [SKIP] 无法导入 gen_step4_disp（%s）" % e)
else:
    (ROOT / "tmp").mkdir(exist_ok=True)
    tmp = Path(tempfile.mkdtemp(dir=str(ROOT / "tmp"), prefix="stepconf-"))
    lat = [[3.14, 0.0, 0.0], [-1.57, 2.72, 0.0], [0.0, 0.0, 25.0]]
    pos = "Mo S\n1.0\n" + "".join("  %.10f %.10f %.10f\n" % tuple(r) for r in lat) \
        + "Mo S\n1 2\nDirect\n0 0 0\n0.33 0.66 0.1\n0.33 0.66 0.6\n"
    (tmp / "POSCAR").write_text(pos, encoding="utf-8")

    def write_yaml(sm):
        (tmp / "phono3py_disp.yaml").write_text(
            json.dumps({"unit_cell": {"lattice": lat}, "supercell_matrix": sm}),
            encoding="utf-8")
    write_yaml([[5, 0, 0], [0, 5, 0], [0, 0, 1]])
    reps = [5, 5, 1]

    def plan(**kw):
        d = {"method": "alm", "oversample": 3, "fd_distance": 0.03,
             "alm_cut2": None, "alm_cut3": 6.0, "fc3_cutoff_pair": None}
        d.update(kw)
        (tmp / "disp_plan.json").write_text(json.dumps(d), encoding="utf-8")

    plan()
    base = mkconf({"OVERSAMPLE": 3, "ALM_CUT2": None, "ALM_CUT3": 6.0,
                   "FD_DISTANCE": 0.03, "FC3_CUTOFF_PAIR": None})
    ck(g._dataset_matches(tmp, reps, "alm", base) is True, "完全一致 => 匹配（幂等跳过）")
    ck(g._dataset_matches(tmp, reps, "findiff", base) is False, "method 变了 => 重建")
    ck(g._dataset_matches(tmp, [4, 4, 1], "alm", base) is False, "超胞系数变了 => 重建")

    # 核心回归：plan 记"不截断"，conf 改成 4.0 —— 必须重建
    ck(g._dataset_matches(tmp, reps, "alm", mkconf({"ALM_CUT2": 4.0})) is False,
       "ALM_CUT2: plan=None, conf=4.0 => 重建（旧代码漏判）")
    plan(alm_cut2=4.0)
    ck(g._dataset_matches(tmp, reps, "alm", base) is False,
       "ALM_CUT2: plan=4.0, conf=None => 重建（旧代码漏判）")
    ck(g._dataset_matches(tmp, reps, "alm", mkconf({"ALM_CUT2": 4.0})) is True,
       "ALM_CUT2: plan=4.0, conf=4.0 => 匹配")
    slim = {"method": "alm", "oversample": 3}
    (tmp / "disp_plan.json").write_text(json.dumps(slim), encoding="utf-8")
    ck(g._dataset_matches(tmp, reps, "alm", mkconf({"ALM_CUT2": 4.0})) is True,
       "老 plan 缺 ALM_CUT2 键 => 不据此重建（向后兼容）")
    plan(alm_cut2=None)
    ck(g._dataset_matches(tmp, reps, "alm", base) is True, "base 复测仍匹配")

print("")
if FAIL:
    print("suite_stepconf: FAILED %d 项" % len(FAIL))
    sys.exit(1)
print("suite_stepconf: ALL PASS")
