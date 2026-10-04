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


class RegenPendingTests(unittest.TestCase):
    """V158：归档后在原目录重新生成了输入、还没提交 -> 输入就绪（TODO），不是 FAIL。
    起因：P1_Mo-MoS2 的 S4_wave —— 作业槽满时 auto-advance 只预生成输入，submit.sh 比归档新（不再算 STALE），
    旧 slurm-*.out 还在 -> FAIL，auto-advance 不重试，一直卡着。"""

    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.d = Path(self.td.name) / "step4_wave"
        self.now = time.time()
        _touch(self.d / "slurm-3918000.out", self.now - 86400)
        _touch(self.d / "amset.log", self.now - 86400)
        _touch(self.d / ("wavefunction.h5.stale-upstream-step3_uniform-" + _tag(self.now - 3600)),
               self.now - 86400)

    def test_pregenerated_after_archive(self):
        self.assertIsNotNone(C.stale_upstream(str(self.d)))                    # 还没重新生成：STALE
        _touch(self.d / "submit.sh", self.now - 60)                             # 预生成
        self.assertIsNone(C.stale_upstream(str(self.d)))
        self.assertTrue(C.regen_pending(str(self.d)))

    def test_ran_or_gen_failed_after_regen(self):
        _touch(self.d / "submit.sh", self.now - 600)
        _touch(self.d / "slurm-3919500.out", self.now - 30)                     # 交了、跑过（失败）
        self.assertFalse(C.regen_pending(str(self.d)))
        (self.d / "slurm-3919500.out").unlink()
        _touch(self.d / ".autozt_gen_failed", self.now - 30)                    # 再生成失败
        self.assertFalse(C.regen_pending(str(self.d)))

    def test_needs_archive_and_submit_script(self):
        self.assertFalse(C.regen_pending(str(self.d)))                          # 没有提交脚本
        _touch(self.d / "submit.sh", self.now - 60)
        self.assertFalse(C.regen_pending(str(self.d), "run.sh"))
        plain = Path(self.td.name) / "step5_dielect"                           # 普通失败步骤手动重新生成：照旧 FAIL
        _touch(plain / "slurm-1.out", self.now - 7200)
        _touch(plain / "OUTCAR", self.now - 7200)
        _touch(plain / "submit.sh", self.now - 60)
        self.assertFalse(C.regen_pending(str(plain)))

    def test_archive_after_regen_is_stale_not_pending(self):
        _touch(self.d / "submit.sh", self.now - 7200)                           # 先生成、后被归档
        self.assertIsNotNone(C.stale_upstream(str(self.d)))
        self.assertFalse(C.regen_pending(str(self.d)))


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

    def test_regen_ready_is_todo(self):
        self.assertEqual(wf.step_state(self._s(regen_ready=True, has_incar=False), False), ("TODO", "TODO"))
        self.assertEqual(wf.step_state(self._s(regen_ready=True, has_outcar=True), False), ("TODO", "TODO"))
        self.assertEqual(wf.step_state(self._s(regen_ready=True), True), ("----", "WAIT"))
        self.assertEqual(wf.step_state(self._s(regen_ready=True, done=True), False), ("OK", "OK"))

    def test_failed_regen_marks(self):
        calls = []
        with patch("autozt.run_remote", lambda cfg, line, host=None: calls.append(line) or (0, "")):
            wf._mark_gen_failed({}, {"stale": True, "dir": "/w/m/step8.3_output", "name": "step8.3_output"})
            wf._mark_gen_failed({}, {"stale": False, "dir": "/w/m/x", "name": "x"})
            wf._mark_gen_failed({}, {"regen_ready": True, "dir": "/w/m/step4_wave", "name": "step4_wave"})
        self.assertEqual(len(calls), 2)
        self.assertIn(".autozt_gen_failed", calls[0])
        self.assertIn("step4_wave/.autozt_gen_failed", calls[1])


class CollectorWiringTests(unittest.TestCase):
    def test_probe_sets_stale(self):
        src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
        i = src.index("_st = stale_upstream(d)")
        self.assertLess(i, src.index("steps.append(f)", i))
        self.assertIn('f["plot_error"] = False', src[i:i + 400])
        k = src.index("regen_pending(d, f[\"submit\"])", i)
        self.assertLess(k, src.index("steps.append(f)", i))
        self.assertIn('f["regen_ready"] = True', src[k:k + 200])


if __name__ == "__main__":
    unittest.main(verbosity=2)
