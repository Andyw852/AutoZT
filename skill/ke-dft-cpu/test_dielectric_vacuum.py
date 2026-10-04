#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_dielectric_vacuum.py —— V164：2D 板的真空方向不进 DFPT/独立粒子交叉比值（不需要 VASP）。

用法：python test_dielectric_vacuum.py      退出码 0 = 全部 PASS。

起因：WS2 的 S5.1 告警"DFPT 与独立粒子 eps 差 103.3%（>25%）—— 检查 NBANDS / k 网格"。实为 z（真空）方向：
DFPT 含局域场（层与真空串联 1/(f/εs+1-f) ≈ 1.21），独立粒子是体积平均 1+f(εs-1) ≈ 2.4 —— 物理如此，
每个 2D 材料都会误报。面内照旧按 25% 判。
"""
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "step5_dielect"))
import validate_dielectric as V  # noqa: E402

WS2 = """WS2 slab
1.0
3.18 0.0 0.0
-1.59 2.754 0.0
0.0 0.0 30.0
W S
1 2
Direct
0.0 0.0 0.5
0.3333 0.6667 0.5524
0.3333 0.6667 0.4476
"""
SI = """Si
5.43
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
Si
2
Direct
0.0 0.0 0.0
0.25 0.25 0.25
"""
WS2_CART = """WS2 slab cart
1.0
3.18 0.0 0.0
-1.59 2.754 0.0
0.0 0.0 30.0
1 2
Selective dynamics
Cartesian
0.0 0.0 15.0 T T T
0.0 1.836 16.57 T T T
0.0 1.836 13.43 T T T
"""


def _txt(ip):
    return ("HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE, excluding Hartree and local field effects)\n"
            " ------------------------------------------------------\n"
            + "".join("  %s\n" % "  ".join("%.6f" % v for v in row) for row in ip))


class VacuumTests(unittest.TestCase):
    def _dir(self, text):
        d = Path(tempfile.mkdtemp())
        (d / "POSCAR").write_text(text)
        return d

    def test_axes(self):
        self.assertEqual(V.vacuum_axes(self._dir(WS2)), (2,))
        self.assertEqual(V.vacuum_axes(self._dir(WS2_CART)), (2,))       # VASP4 + Selective + Cartesian
        self.assertEqual(V.vacuum_axes(self._dir(SI)), ())

    def test_ws2_z_not_counted(self):
        ei = [[3.56, 0, 0], [0, 3.56, 0], [0, 0, 1.204]]
        txt = _txt([[3.60, 0, 0], [0, 3.60, 0], [0, 0, 2.447]])
        phys, warns = V.born_asr_and_crosscheck(txt, ei, (2,))
        self.assertLess(phys["cross_ratio"], 0.25)
        self.assertGreater(phys["cross_ratio_vacuum_axis"], 1.0)
        self.assertFalse([w for w in warns if "独立粒子" in w])
        _p, warns_old = V.born_asr_and_crosscheck(txt, ei)                   # 不知道真空方向：照旧告警
        self.assertTrue([w for w in warns_old if "独立粒子" in w])

    def test_in_plane_still_checked(self):
        ei = [[3.56, 0, 0], [0, 3.56, 0], [0, 0, 1.204]]
        txt = _txt([[5.0, 0, 0], [0, 5.0, 0], [0, 0, 2.447]])
        _p, warns = V.born_asr_and_crosscheck(txt, ei, (2,))
        self.assertTrue([w for w in warns if "独立粒子" in w])


if __name__ == "__main__":
    unittest.main(verbosity=2)
