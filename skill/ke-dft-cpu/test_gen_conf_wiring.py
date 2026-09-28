#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_gen_conf_wiring.py —— step.conf 覆盖通道的**静态接线**自检（不跑 gen、不需要任何数据）。

用法：python test_gen_conf_wiring.py      退出码 0 = 全部 PASS。

2026-09-28 查出的两处"静默失效"都是接线错误，运行时不报错、只是覆盖不生效：
  · gen_step10_amset.py：main() 给 EPS_INF_OVERRIDE 赋值却没声明 global —— 赋到局部变量，
    _apply_eps_inf_override() 读到的仍是模块级 None，step.conf 的 ε∞ 覆盖**永远不生效**；
  · gen_step14_amset2d.py：SPEC 里丢了 "INTERPOLATION_FACTOR" —— main() 里 _p[...] 抛 KeyError，
    被 except 吞掉，**它之后的全部覆盖**（WAVEFUNCTION_FULL / UNITY_OVERLAP / MESH_* /
    SCATTERING / WRITE_MESH / DOPING / TEMPERATURES）一起静默失效。
本测试对技能里每个用 stepconf.load 的 gen 脚本检查两条：
  1) 用 stepconf.load 的返回值取的键，必须都在 SPEC 里声明；
  2) 任一函数里给**模块级大写常量**赋值，必须先 global 声明（否则是局部变量，覆盖丢失）。
"""
import ast
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _gen_files():
    for p in sorted(ROOT.rglob("*.py")):
        if p.name.startswith("test_"):
            continue
        s = p.read_text(encoding="utf-8", errors="ignore")
        if "stepconf.load(" in s and "SPEC" in s:
            yield p, s


def _spec_keys(tree):
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "SPEC"
                                                for t in node.targets):
            if isinstance(node.value, ast.Dict):
                return {k.value for k in node.value.keys
                        if isinstance(k, ast.Constant) and isinstance(k.value, str)}
    return None


def _conf_vars(func):
    """函数里 `x = stepconf.load(...)` 的变量名。"""
    out = set()
    for n in ast.walk(func):
        if isinstance(n, ast.Assign) and isinstance(n.value, ast.Call):
            f = n.value.func
            if isinstance(f, ast.Attribute) and f.attr == "load" and \
                    isinstance(f.value, ast.Name) and f.value.id == "stepconf":
                out |= {t.id for t in n.targets if isinstance(t, ast.Name)}
    return out


def _module_upper_names(tree):
    names = set()
    for node in tree.body:
        targets = []
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        for t in targets:
            if isinstance(t, ast.Name) and t.id.isupper():
                names.add(t.id)
    return names


class ConfWiringTests(unittest.TestCase):
    def test_keys_declared_in_spec(self):
        bad = []
        n_files = 0
        for p, s in _gen_files():
            tree = ast.parse(s)
            spec = _spec_keys(tree)
            if spec is None:
                continue
            n_files += 1
            for func in [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]:
                cv = _conf_vars(func)
                if not cv:
                    continue
                for n in ast.walk(func):
                    if isinstance(n, ast.Subscript) and isinstance(n.value, ast.Name) \
                            and n.value.id in cv:
                        k = n.slice
                        k = k.value if isinstance(k, ast.Index) else k       # py<3.9
                        if isinstance(k, ast.Constant) and isinstance(k.value, str) \
                                and k.value not in spec and k.value != "STEP":
                            bad.append("%s:%d %s[%r] 不在 SPEC 里"
                                       % (p.relative_to(ROOT), n.lineno, n.value.id, k.value))
        self.assertGreater(n_files, 3)
        self.assertEqual(bad, [], "\n" + "\n".join(bad))

    def test_module_constants_declared_global(self):
        bad = []
        for p, s in _gen_files():
            tree = ast.parse(s)
            upper = _module_upper_names(tree)
            for func in [n for n in tree.body if isinstance(n, ast.FunctionDef)]:
                glob = set()
                for n in ast.walk(func):
                    if isinstance(n, ast.Global):
                        glob |= set(n.names)
                for n in ast.walk(func):
                    tgts = []
                    if isinstance(n, ast.Assign):
                        tgts = n.targets
                    elif isinstance(n, (ast.AugAssign, ast.AnnAssign)):
                        tgts = [n.target]
                    for t in tgts:
                        if isinstance(t, ast.Name) and t.id in upper and t.id not in glob:
                            bad.append("%s:%d %s() 给模块常量 %s 赋值但没声明 global"
                                       % (p.relative_to(ROOT), n.lineno, func.name, t.id))
        self.assertEqual(bad, [], "\n" + "\n".join(bad))


if __name__ == "__main__":
    unittest.main(verbosity=2)
