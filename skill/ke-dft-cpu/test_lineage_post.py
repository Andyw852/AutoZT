#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_lineage_post.py —— V163：lineage_check 的事后复核覆盖 S8/S8.4/S8.2/S8.1 本身（不需要 VASP/AMSET）。

用法：python test_lineage_post.py      退出码 0 = 全部 PASS。

起因：DOWNSTREAM 里有 24 条非 dir 失效边，lineage_check 一条都不复核 —— P1_Mo-MoS2 "4 项全 OK" 其实没验证
transport.json 是用最新的 S4/S5/S6/S7.1 算的；V162 的 S8.1 假完成正是失效没生效、事后又没查出来。
1) 每条非 dir 的 DOWNSTREAM 边：要么 LINEAGE_DERIVED/POST_DERIVED 有新鲜度规则，要么是 link 边且下游在
   LINKED_PRODUCTS（按软链目标比），要么是 gap 边（带隙一致性另查）；
2) S8 的软链输入比 transport.json 新 -> 报；断链不报；
3) S5/S6 比 transport.json 新 -> 报。
"""
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


def _w(m, rel, age_s, text="x"):
    f = m / rel
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(text)
    t = time.time() - age_s
    os.utime(f, (t, t))
    return f


class CoverageTests(unittest.TestCase):
    def test_every_non_dir_edge_has_a_freshness_rule(self):
        rules = {(s, src) for s, _m, src, _p in list(kc.LINEAGE_DERIVED) + list(kc.POST_DERIVED)}
        linked = {s for s, _m in kc.LINKED_PRODUCTS}
        miss = []
        for up, lst in kc.DOWNSTREAM.items():
            for down, mode in lst:
                if mode == "dir" or (down, up) in rules:
                    continue
                if mode == "link" and down in linked:
                    continue
                if mode == "gap" and down in ("step8_amset", "step8.4_amset2d"):   # post_lineage 查带隙一致性
                    continue
                if mode == "link" and down in ("step4_wave", "step4b_wave_full"):   # LINEAGE_DERIVED 按 vasprun/WAVECAR
                    continue
                miss.append("%s -> %s (%s)" % (up, down, mode))
        self.assertEqual(miss, [], "这些失效边 lineage_check 不复核：\n  " + "\n  ".join(miss))


class LinkedTests(unittest.TestCase):
    def setUp(self):
        self.m = Path(tempfile.mkdtemp())

    def test_newer_linked_input_reported(self):
        m = self.m
        _w(m, "step8_amset/transport.json", 3600)
        _w(m, "step4b_wave_full/wavefunction.h5", 60)                      # S4b 后来重算
        _w(m, "step4_wave/wavefunction.h5", 7200)
        (m / "step8_amset" / "wavefunction.h5").symlink_to("../step4b_wave_full/wavefunction.h5")
        probs = kc.post_lineage(m)
        self.assertEqual([s for s, _w2 in probs], ["step8_amset"])
        self.assertIn("step4b_wave_full/wavefunction.h5", probs[0][1])

    def test_old_linked_input_and_dangling_link_ok(self):
        m = self.m
        _w(m, "step8_amset/transport.json", 60)
        _w(m, "step4_wave/wavefunction.h5", 3600)
        (m / "step8_amset" / "wavefunction.h5").symlink_to("../step4_wave/wavefunction.h5")
        (m / "step8_amset" / "deformation.h5").symlink_to("../step7b_deform_read/deformation.h5")   # 本地没拉：断链
        self.assertEqual(kc.post_lineage(m), [])

    def test_s5_s6_newer_than_transport(self):
        m = self.m
        _w(m, "step8.4_amset2d/transport.json", 3600)
        _w(m, "step6_elastic/OUTCAR", 60)
        probs = kc.post_lineage(m)
        self.assertEqual([s for s, _w2 in probs], ["step8.4_amset2d"])
        self.assertIn("step6_elastic/OUTCAR", probs[0][1])


if __name__ == "__main__":
    unittest.main(verbosity=2)
