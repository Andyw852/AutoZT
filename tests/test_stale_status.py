# -*- coding: utf-8 -*-
"""V153：上游重算后被归档的步骤显示 STALE（kind=PREP），auto-advance 会重新生成，而不是显示 FAIL。

起因：ke-dft-cpu 的 gen 把下游完成标记改名 *.stale-upstream-<上游>-<时间>（V140/V148），整目录归档
<目录>.stale-upstream-<上游>-<时间>（V153）。采集端以前把"目录在、标记不在"判成 FAIL（画图步 plot_error、
作业步留着 slurm-*.out），而 auto-advance 不重试 FAIL —— "会重新排队"其实要人手 retry
（P1_Mo-MoS2 的任务图里 6 个 FAIL 全是这种）。
离线：只用临时目录，不连超算。
"""
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from autozt import workflow as wf  # noqa: E402


def _load_collector():
    """_collector_remote.py 是推到远端直接跑的脚本（导入即执行 main）：只取 main 之前的定义。"""
    import types
    src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
    mod = types.ModuleType("_collector_defs")
    exec(compile(src[:src.index("\ndef main():")], "_collector_remote.py", "exec"), mod.__dict__)
    return mod


C = _load_collector()


def _tag(ts):
    return time.strftime("%Y%m%d%H%M%S", time.localtime(ts))


def _touch(p, ts):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("x")
    os.utime(p, (ts, ts))


class StaleDetectTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.root = Path(self.td.name)
        self.now = time.time()

    def test_archived_marker_is_stale(self):
        d = self.root / "step8.1_boltztrap"
        _touch(d / "slurm-1.out", self.now - 7200)
        _touch(d / ("boltztrap_crta.json.stale-upstream-step8.2_dpt-" + _tag(self.now - 60)), self.now - 7200)
        up, ts = C.stale_upstream(str(d))
        self.assertEqual(up, "step8.2_dpt")
        self.assertAlmostEqual(ts, self.now - 60, delta=2)

    def test_regenerated_or_failed_regen_is_not_stale(self):
        d = self.root / "step8_amset"
        _touch(d / ("transport.json.stale-upstream-step7b_deform_read-" + _tag(self.now - 600)), self.now - 7200)
        _touch(d / "INCAR", self.now - 10)                                      # 之后重新生成过
        self.assertIsNone(C.stale_upstream(str(d)))
        d2 = self.root / "step8.3_output"
        _touch(d2 / ("comparison_300K.png.stale-upstream-step8_amset-" + _tag(self.now - 600)), self.now - 7200)
        _touch(d2 / ".autozt_gen_failed", self.now - 5)                         # 重新生成失败
        self.assertIsNone(C.stale_upstream(str(d2)))

    def test_whole_dir_archived(self):
        par = self.root / "step2_bandgap"
        (par / ("step2.3_hse.stale-upstream-step2_bandgap_step2.2_pbe-" + _tag(self.now - 30))).mkdir(parents=True)
        up, _ts = C.stale_upstream(str(par / "step2.3_hse"))
        self.assertEqual(up, "step2_bandgap_step2.2_pbe")
        self.assertIsNone(C.stale_upstream(str(par / "step2.2_pbe")))
        self.assertIsNone(C.stale_upstream(str(self.root / "nope")))

    def test_other_archive_kinds_ignored(self):
        d = self.root / "step3_uniform"
        _touch(d / "vasprun.xml.stale-grid-47x47x3", self.now - 7200)          # 本步自己换网格，不是等上游
        _touch(d / ("vasprun.xml.stale-input-" + _tag(self.now - 60)), self.now - 7200)
        self.assertIsNone(C.stale_upstream(str(d)))


class StepStateTests(unittest.TestCase):
    def _s(self, **kw):
        s = {"done": False, "has_outcar": False, "has_slurm_out": True, "has_incar": True, "job": None}
        s.update(kw)
        return s

    def test_stale_beats_fail(self):
        self.assertEqual(wf.step_state(self._s(stale=True, plot_error=True), False), ("STALE", "PREP"))
        self.assertEqual(wf.step_state(self._s(stale=True), True), ("----", "WAIT"))
        self.assertEqual(wf.step_state(self._s(), False), ("FAIL", "FAIL"))         # 没归档：照旧 FAIL
        self.assertEqual(wf.step_state(self._s(stale=True, done=True), False), ("OK", "OK"))
        self.assertEqual(wf.step_state(self._s(stale=True, scancel=True), False), ("scancel", "SCANCEL"))

    def test_failed_regen_marks(self):
        calls = []
        with patch("autozt.run_remote", lambda cfg, line, host=None: calls.append(line) or (0, "")):
            wf._mark_gen_failed({}, {"stale": True, "dir": "/w/m/step8.3_output", "name": "step8.3_output"})
            wf._mark_gen_failed({}, {"stale": False, "dir": "/w/m/x", "name": "x"})
        self.assertEqual(len(calls), 1)
        self.assertIn(".autozt_gen_failed", calls[0])


class CollectorWiringTests(unittest.TestCase):
    def test_probe_sets_stale(self):
        src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
        i = src.index("_st = stale_upstream(d)")
        self.assertLess(i, src.index("steps.append(f)", i))
        self.assertIn('f["plot_error"] = False', src[i:i + 400])


if __name__ == "__main__":
    unittest.main(verbosity=2)
