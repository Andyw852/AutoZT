#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_stale_and_mem.py —— ke_common 的"输入变了/上游重算 -> 失效"与内存粗估自检（纯本地、无数据）。

用法：python test_stale_and_mem.py      退出码 0 = 全部 PASS。
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
import ke_common as kc  # noqa: E402

POSCAR = """Si
1.0
0.0 2.715 2.715
2.715 0.0 2.715
2.715 2.715 0.0
Si
2
Direct
0.0 0.0 0.0
0.25 0.25 0.25
"""
INCAR = "ENCUT = 400\nISYM = 2\nNCORE = 4\n# comment\n"
KPOINTS = "auto\n0\nGamma\n 34 34 34\n"


def _s3(d):
    s3 = d / "step3_uniform"
    s3.mkdir(parents=True)
    (s3 / "POSCAR").write_text(POSCAR)
    (s3 / "INCAR").write_text(INCAR)
    (s3 / "KPOINTS").write_text(KPOINTS)
    for f in ("WAVECAR", "vasprun.xml", "OUTCAR"):
        (s3 / f).write_text("x")
    return s3


def _s4_s8(d):
    s4 = d / "step4_wave"
    s4.mkdir()
    for f in ("WAVECAR", "vasprun.xml"):
        (s4 / f).symlink_to(os.path.relpath(d / "step3_uniform" / f, s4))
    (s4 / "wavefunction.h5").write_text("h5")
    s8 = d / "step8_amset"
    s8.mkdir()
    (s8 / "wavefunction.h5").symlink_to(os.path.relpath(s4 / "wavefunction.h5", s8))
    (s8 / "transport.json").write_text("{}")
    s84 = d / "step8.4_amset2d"                      # 走全网格分支：不链 step4_wave
    s84.mkdir()
    (s84 / "transport.json").write_text("{}")
    return s4, s8, s84


class StaleTests(unittest.TestCase):
    def test_cosmetic_changes_do_not_count(self):
        d = Path(tempfile.mkdtemp())
        s3 = _s3(d)
        snap = kc.snapshot_inputs(s3)
        # 注释、空白、并行键（换宿主机）变化 -> 不算输入变化
        (s3 / "INCAR").write_text("ISYM   = 2\nENCUT=400\nNCORE = 1\nKPAR = 2\n")
        (s3 / "POSCAR").write_text(POSCAR.replace("Si\n1.0", "Si aligned\n1.0"))
        self.assertEqual(kc.inputs_changed(s3, snap), [])

    def test_real_change_archives_and_invalidates_linked_downstream(self):
        d = Path(tempfile.mkdtemp())
        s3 = _s3(d)
        s4, s8, s84 = _s4_s8(d)
        snap = kc.snapshot_inputs(s3)
        (s3 / "POSCAR").write_text(POSCAR.replace("0.25 0.25 0.25", "0.26 0.25 0.25"))
        changed, n, down = kc.finalize_stale_inputs(d, s3, "step3_uniform", snap)
        self.assertEqual(changed, ["POSCAR"])
        self.assertEqual(n, 3)                             # WAVECAR / vasprun.xml / OUTCAR
        self.assertFalse((s3 / "WAVECAR").exists())
        self.assertFalse((s4 / "wavefunction.h5").exists())   # S4 链着 S3 -> 失效
        self.assertFalse((s8 / "transport.json").exists())    # S8 链着 S4 -> 传递失效
        self.assertTrue((s84 / "transport.json").exists())    # S8.4 没链 S3/S4 -> 不动
        self.assertTrue(any(x.name.startswith("transport.json.stale-upstream-")
                            for x in s8.iterdir()))
        self.assertIn(("step4_wave", "wavefunction.h5"), down)

    def test_first_generation_is_noop(self):
        d = Path(tempfile.mkdtemp())
        (d / "step3_uniform").mkdir()
        snap = kc.snapshot_inputs(d / "step3_uniform")
        s3 = d / "step3_uniform"
        (s3 / "POSCAR").write_text(POSCAR)
        self.assertEqual(kc.finalize_stale_inputs(d, s3, "step3_uniform", snap), ([], 0, []))

    def test_grid_archive_alone_invalidates(self):
        d = Path(tempfile.mkdtemp())
        s3 = _s3(d)
        s4, s8, _ = _s4_s8(d)
        snap = kc.snapshot_inputs(s3)
        kc.finalize_stale_inputs(d, s3, "step3_uniform", snap, n_grid_archived=5)
        self.assertFalse((s8 / "transport.json").exists())


class StaleDeformationTests(unittest.TestCase):
    """patch_stale_dp（V122）：S7.1 重新生成 -> 软链它的 S8/S8.4 与 S8.2（DPT）失效，S8.3 递归失效；
    没链它的（例如走 MANUAL 形变势的 S8）不动。"""

    def test_s71_regen_invalidates_consumers(self):
        d = Path(tempfile.mkdtemp())
        s71 = d / "step7b_deform_read"
        s71.mkdir()
        (s71 / "deformation_vac.h5").write_text("h5")
        (s71 / "band_edges.json").write_text("{}")
        s84 = d / "step8.4_amset2d"
        s84.mkdir()
        (s84 / "deformation.h5").symlink_to(os.path.relpath(s71 / "deformation_vac.h5", s84))
        (s84 / "transport.json").write_text("{}")
        (s84 / "intrinsic_transport.json").write_text("{}")
        s8 = d / "step8_amset"                              # 不链 S7.1 -> 不动
        s8.mkdir()
        (s8 / "transport.json").write_text("{}")
        s82 = d / "step8.2_dpt"
        s82.mkdir()
        (s82 / "dpt_result.json").write_text("{}")
        s83 = d / "step8.3_output"
        s83.mkdir()
        (s83 / "comparison_300K.png").write_text("png")
        down = kc.invalidate_downstream(d, "step7b_deform_read", "test")
        self.assertFalse((s84 / "transport.json").exists())
        self.assertFalse((s84 / "intrinsic_transport.json").exists())
        self.assertTrue((s8 / "transport.json").exists())
        self.assertFalse((s82 / "dpt_result.json").exists())
        self.assertFalse((s83 / "comparison_300K.png").exists())     # 经 S8.2 递归
        self.assertIn(("step8.4_amset2d", "transport.json"), down)

    def test_s71_gen_calls_invalidation_before_regenerating(self):
        src = (_ROOT / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py").read_text(
            encoding="utf-8")
        i = src.index('kc.invalidate_downstream(cwd, OUTDIR_NAME')
        j = src.index("amset deform read undeformed")
        self.assertLess(i, j)


class MemTests(unittest.TestCase):
    def test_calibration_points_inside_band(self):
        lo, hi = kc.estimate_amset_mem_gib([263, 263, 21])       # MoS2 实测 91.1 GiB
        self.assertLessEqual(lo, 91.1 * 1.01)
        self.assertGreaterEqual(hi, 91.1)
        lo, hi = kc.estimate_amset_mem_gib([109, 109, 109])      # GaAs 实测 150.9 GiB
        self.assertLessEqual(lo, 150.9)
        self.assertGreaterEqual(hi, 150.9 * 0.99)

    def test_levels_and_exclusive(self):
        self.assertEqual(kc.mem_guard([100, 100, 100], warn_gib=300)[2], "ok")
        self.assertEqual(kc.mem_guard([329, 329, 27], warn_gib=300)[2], "warn")
        self.assertEqual(kc.mem_guard([400, 400, 400], warn_gib=300, limit_gib=500)[2], "error")
        d = Path(tempfile.mkdtemp())
        sub = d / "submit.sh"
        sub.write_text("#!/bin/bash\n#SBATCH -J x\n#SBATCH --ntasks=24\necho hi\n")
        self.assertTrue(kc.add_exclusive(sub))
        self.assertFalse(kc.add_exclusive(sub))                  # 幂等
        lines = sub.read_text().splitlines()
        self.assertEqual(lines[3], "#SBATCH --exclusive")


if __name__ == "__main__":
    unittest.main(verbosity=2)
