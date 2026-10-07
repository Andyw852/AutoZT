#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_preflight_elastic.py —— overlap_preflight 第 5 项（弹性张量）在 S8.4 插件路径上不再误报（V180）。

用法：python test_preflight_elastic.py     退出码 0 = 全部 PASS。

起因（CrS₂ S8.4，2026-10-07）：preflight 报"通过（有告警）"，唯一的告警是完整 3x3 Christoffel 最小特征值
-0.171（方向 (1,0,0)，即面外剪切 C55 为负）。S8.4 的 amset2d_plugin 里 ADP/PIE 只用面内 2x2，这个负值进不了计算；
告警还叫人"retry 重生成，或走 step8.4_amset2d 插件"——本步就是插件。每个面外剪切为负的二维材料都会被这样误报。
张量在测试里构造（量级照 CrS₂ 的原始 slab 值），不读任何真实文件。
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT), str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import overlap_preflight as pf

from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.inputs import Poscar


def voigt(c11=70.0, c12=20.0, c33=1.0, c44=-0.171, c66=25.0):
    """六方二维 slab 的标准 Voigt 张量（XX YY ZZ YZ XZ XY），面外剪切 C44=C55 可为负。"""
    m = [[0.0] * 6 for _ in range(6)]
    m[0][0] = m[1][1] = c11
    m[0][1] = m[1][0] = c12
    m[2][2] = c33
    m[3][3] = m[4][4] = c44
    m[5][5] = c66
    return m


def run_dir(elastic, plugin_file=True, plugin_key="amset2d_plugin.py"):
    base = Path(tempfile.mkdtemp())
    (base / "step3_uniform").mkdir()
    st = Structure(Lattice.hexagonal(3.0, 20.0), ["Cr", "S", "S"],
                   [[0, 0, .5], [1 / 3, 2 / 3, .57], [1 / 3, 2 / 3, .43]])
    Poscar(st).write_file(str(base / "step3_uniform" / "POSCAR"))
    out = base / "step8.4_amset2d"
    out.mkdir()
    (out / "settings.yaml").write_text("unity_overlap: false\nelastic_constant: %s\n" % json.dumps(elastic))
    rec = {"plugin": plugin_key} if plugin_key else {}
    (out / "2d_correction.json").write_text(json.dumps(rec))
    if plugin_file:
        (out / "amset2d_plugin.py").write_text("# stub\n")
    return out


def item5(lines):
    i = next(n for n, ln in enumerate(lines) if ln.startswith("  5) "))
    return lines[i:]


class ElasticLinesTests(unittest.TestCase):
    def test_plugin_path_negative_outofplane_is_not_a_warning(self):
        st, lines = pf.elastic_check_lines(voigt(), True, True)
        txt = "\n".join(lines)
        self.assertEqual(st, "ok", txt)
        self.assertIn("-0.1710", txt)                                   # 数照样记下
        self.assertIn("(1, 0, 0)", txt)
        self.assertNotIn("[WARN]", txt)
        self.assertIn("不影响结果", txt)

    def test_standard_2d_path_still_warns_with_current_wording(self):
        st, lines = pf.elastic_check_lines(voigt(), True, False)
        txt = "\n".join(lines)
        self.assertEqual(st, "warn", txt)
        self.assertIn("取绝对值", txt)
        self.assertNotIn("置零", txt)                                   # 旧说法（V76/V77 已改成取绝对值）

    def test_inplane_misorder_blocks_even_on_plugin_path(self):
        # VASP 顺序没重排：C66 槽里是接近 0 的面外量，C44 槽里是面内剪切
        st, lines = pf.elastic_check_lines(voigt(c44=25.0, c66=0.2), True, True)
        self.assertEqual(st, "error", "\n".join(lines))
        self.assertIn("★ 拦截", "\n".join(lines))

    def test_3d_negative_warns_without_2d_wording(self):
        st, lines = pf.elastic_check_lines(voigt(c33=60.0), False, False)
        txt = "\n".join(lines)
        self.assertEqual(st, "warn", txt)
        self.assertNotIn("面外", txt)
        self.assertIn("ELASTIC_GUARD", txt)

    def test_positive_definite_is_one_line(self):
        st, lines = pf.elastic_check_lines(voigt(c44=0.5), True, False)
        self.assertEqual((st, len(lines)), ("ok", 1), "\n".join(lines))

    def test_missing_elastic_skips(self):
        st, lines = pf.elastic_check_lines(None, True, True)
        self.assertEqual(st, "ok")
        self.assertIn("跳过", lines[0])


class PluginDetectTests(unittest.TestCase):
    def test_needs_file_and_record(self):
        self.assertTrue(pf.plugin_2d_active(run_dir(voigt())))
        self.assertFalse(pf.plugin_2d_active(run_dir(voigt(), plugin_file=False)))
        self.assertFalse(pf.plugin_2d_active(run_dir(voigt(), plugin_key=None)))     # S8 标准路径的记录没有 plugin 键
        self.assertFalse(pf.plugin_2d_active(run_dir(voigt(), plugin_key="other.py")))


class RunWiringTests(unittest.TestCase):
    def test_run_uses_plugin_aware_check(self):
        _, lines = pf.run(run_dir(voigt()), None, False, desym_fix=False)
        self.assertNotIn("[WARN]", "\n".join(item5(lines)))
        _, lines = pf.run(run_dir(voigt(), plugin_file=False, plugin_key=None), None, False, desym_fix=False)
        self.assertIn("[WARN]", "\n".join(item5(lines)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
