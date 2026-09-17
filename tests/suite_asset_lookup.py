# -*- coding: utf-8 -*-
"""自测（C 批）：gen_need 跨步骤 src 兜底 + step.conf 的 strict 开关（2026-09-16）。

**只放依赖 autozt / _common 未提交改动的检查**：这样它与那批改动同一次提交能跑通。
ke_common 单份化与 gen_need 声明完整性在 tests/suite_ke_common.py（A 批）。

两件事都是 MoS2 实测暴露的（ke-dft-cpu 在 S3_uniform / S5_dielect 上 gen 失败）：

1. find_asset 原来只按【当前步骤的 src 子目录】找 gen_need 里的文件。
   ke 的 discriminant_common.py 住在 step2_bandgap/step2.15_discriminant/，
   而 step5_dielect / step8_amset / step8.1 / step8.4 的 gen 都 import 它 ——
   永远找不到，gen 直接失败。现在加了"技能目录内按 basename 唯一命中"的兜底。
2. stepconf.load 严格模式下会因**别的步骤**的键（FUNC）报错。材料级 step.conf
   是全技能共用的一份，所以只认一两个键的脚本需要 strict=False。
"""
import glob
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "skill", "_common", "opt"))
from autozt import find_asset          # noqa: E402
import stepconf                        # noqa: E402

FAILS = []


def ok(cond, msg):
    print(("  ✓ " if cond else "  ✗ ") + msg)
    if not cond:
        FAILS.append(msg)


def fake(t_base, sld, sname, fname, steps=None):
    t = {"key": "ke-dft-cpu", "skill_dir": sld, "template_dir": ".",
         "template_layout": "per_step", "_base_dir": t_base,
         "steps_cfg": steps or [{"name": "step5_dielect", "src": "step5_dielect"}]}
    m = {"ps": None, "template_map": {}, "_skill_dir_local": None,
         "hpc_name": None, "_seg": {}}
    return find_asset({"_config_dir": None}, t, m, fname, sname)


def test_find_asset():
    print("[1] find_asset：别的步骤 src 里的模块要能找到")
    base = os.path.join(ROOT, "skill")
    p = fake(base, "ke-dft-cpu", "step5_dielect", "discriminant_common.py")
    ok(p is not None and p.endswith("step2_bandgap/step2.15_discriminant/discriminant_common.py"),
       "ke：step5_dielect 找到了 step2.15 目录里的 discriminant_common.py → %s"
       % (p.replace(ROOT + "/", "") if p else None))

    p2 = fake(base, "ke-dft-cpu", "step5_dielect", "gen_step8_dielect.py")
    ok(p2 is not None and p2.endswith("step5_dielect/gen_step8_dielect.py"),
       "回归：正常按 src 命中的文件不受兜底影响（仍取本步骤那份）")

    d = tempfile.mkdtemp(prefix="tf_asset_")
    try:
        os.makedirs(os.path.join(d, "stepA"))
        os.makedirs(os.path.join(d, "stepB"))
        open(os.path.join(d, "stepA", "dup.py"), "w").write("A\n")
        open(os.path.join(d, "stepB", "dup.py"), "w").write("B\n")
        p3 = fake(os.path.dirname(d), os.path.basename(d), "stepC", "dup.py")
        ok(p3 is None, "多份同名且内容不同时兜底不猜（返回 None，保留原来的报错）")
        open(os.path.join(d, "stepA", "unique.py"), "w").write("x")
        p4 = fake(os.path.dirname(d), os.path.basename(d), "stepC", "unique.py")
        ok(p4 is not None and p4.endswith("stepA/unique.py"),
           "唯一命中时兜底生效（stepC 拿到 stepA/unique.py）")

        # 多份同名但内容完全一致 → 任取一份（ke_common.py 这种按步骤各放一份的公共模块）
        for sub in ("stepA", "stepB"):
            open(os.path.join(d, sub, "same.py"), "w").write("identical\n")
        p5 = fake(os.path.dirname(d), os.path.basename(d), "stepC", "same.py")
        ok(p5 is not None and p5.endswith("same.py"),
           "多份同名但内容相同 → 任取一份（内容无歧义）")
        open(os.path.join(d, "stepB", "same.py"), "w").write("different\n")
        p6 = fake(os.path.dirname(d), os.path.basename(d), "stepC", "same.py")
        ok(p6 is None, "多份同名且内容不同 → 仍不猜")

        # 真实技能：discriminant_common.py 住在别的步骤的 src 目录里，
        # 而 S5/S8/S8.4 的 gen 都 import 它 —— 兜底必须能找到
        p7 = fake(base, "ke-dft-cpu", "step5_dielect", "discriminant_common.py")
        ok(p7 is not None and p7.endswith("discriminant_common.py"),
           "ke：S5 的 discriminant_common.py 由兜底解析到")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_stepconf_strict():
    print("[2] stepconf：共用 step.conf 里别的步骤的键不该打死只用一两个键的脚本")
    d = tempfile.mkdtemp(prefix="tf_sc_")
    try:
        with open(os.path.join(d, "step.conf"), "w") as fh:
            fh.write("[params]\nBANDGAP = hse\nFUNC = auto\nVACUUM_KZ_MIN = 3\n")
        spec = {"VACUUM_KZ_MIN": (1, "int")}
        ok(stepconf.load(spec, "step3_uniform", d, strict=False)["VACUUM_KZ_MIN"] == 3,
           "strict=False：忽略 FUNC/BANDGAP，取到 VACUUM_KZ_MIN=3")
        raised = False
        try:
            stepconf.load(spec, "step3_uniform", d)
        except SystemExit:
            raised = True
        ok(raised, "strict=True（默认）：仍是老行为，报「不认识的键」")
        ok(stepconf.load(spec, "step3_uniform", d, strict=False).params.get("VACUUM_KZ_MIN") == 3,
           "strict 参数是加载项（默认 True），不影响既有调用方")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    print("== suite_asset_lookup：gen_need 查找兜底 / step.conf strict ==")
    test_find_asset()
    test_stepconf_strict()
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