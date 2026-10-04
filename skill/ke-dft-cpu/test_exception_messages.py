#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_exception_messages.py —— V153：兜底 except 不许只写异常类型名（纯源码检查，不需要 VASP/AMSET）。

用法：python test_exception_messages.py      退出码 0 = 全部 PASS。

起因：WS2 的 S8.2 只留下"二次型拟合异常：IndexError"，定位又绕了一轮（V150）。全技能扫一遍还有 12 处同样
只写 type(e).__name__ 的告警（S8.1 / S8.3 / S8.4 / overlap_preflight），另有几处会悄悄改变结果的 except-pass：
S7/S7.1 读不到 step.conf 的 amset 环境就按主机猜一个，S8/S8.4 写 interpolation_info.json 失败不吭声。
1) 本技能的 .py（AMSET 插件与 tools/ 之外）：凡出现 type(X).__name__，同一条语句里必须也带上 X 本身
   （"%s: %s" % (type(e).__name__, e)）或用 _exc_brief；
2) S7/S7.1 的 amset 环境回退、S8/S8.4 写 interpolation_info.json 失败都要告警，不许 except-pass。
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SKIP = {"amset2d_plugin.py"}                 # 跑在 AMSET 进程里，不动


def _sources():
    for p in sorted(ROOT.rglob("*.py")):
        rel = p.relative_to(ROOT)
        if rel.parts[0] == "tools" or p.name.startswith("test_") or p.name in SKIP or "__pycache__" in rel.parts:
            continue
        yield rel, p.read_text(encoding="utf-8")


class TypeOnlyTests(unittest.TestCase):
    def test_no_type_name_only_messages(self):
        bad = []
        for rel, src in _sources():
            lines = src.splitlines()
            for i, ln in enumerate(lines):
                for m in re.finditer(r"type\((\w+)\)\.__name__", ln):
                    var = m.group(1)
                    stmt = " ".join(lines[max(0, i - 2):i + 2])
                    if re.search(r"type\(%s\)\.__name__\s*,\s*(str\()?%s\b" % (var, var), stmt) \
                            or "_exc_brief(" in stmt or re.search(r"str\(%s\)" % var, stmt):
                        continue
                    bad.append("%s:%d  %s" % (rel, i + 1, ln.strip()))
        self.assertEqual(bad, [], "只写了异常类型名：\n" + "\n".join(bad))

    def test_silent_fallbacks_now_warn(self):
        for rel in ("step7_deform/gen_step9_deform.py", "step7_deform/step7b_read/gen_step9b_deform_read.py"):
            src = (ROOT / rel).read_text(encoding="utf-8")
            body = src[src.index("def _amset_env_src():"):src.index("AMSET_ENV_SRC")]
            self.assertIn("主机探测回退", body, rel)
            self.assertNotRegex(body, r"except Exception:\s*\n\s*pass", rel)
        for rel in ("step8_amset/gen_step10_amset.py", "step8.4_amset2d/gen_step14_amset2d.py"):
            src = (ROOT / rel).read_text(encoding="utf-8")
            i = src.index('(out / "interpolation_info.json").write_text(')
            self.assertIn("写 interpolation_info.json 失败", src[i:i + 400], rel)


if __name__ == "__main__":
    unittest.main(verbosity=2)
