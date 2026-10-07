# -*- coding: utf-8 -*-
"""V164：① retry 对本地即时步（画图/读取/校验）就地重跑；② strict_done_marker 之前的旧格式标记判 STALE（重新校验），
不再永久 FAIL。

起因（WS2 的 S5.1_dievalid）：dielectric_check.json 是 09-23 的旧格式（只有 ok: true），V136 的严格判据判
"invalid completion marker" -> FAIL；auto-advance 不重试 FAIL，retry 只打印"待 start 触发"，start 又因 FAIL 要 -f。
离线。
"""
import json
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
    import types
    src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
    mod = types.ModuleType("_collector_defs")
    exec(compile(src[:src.index("\ndef main():")], "_collector_remote.py", "exec"), mod.__dict__)
    return mod


C = _load_collector()
SC = {"check": "plot", "done_marker": "dielectric_check.json", "strict_done_marker": True}


class RetryLocalStepTests(unittest.TestCase):
    def _retry(self, run):
        calls = []
        s = {"name": "step5_dielect_validate", "label": "S5.1_dievalid", "job": None}
        with patch("autozt.step_cfg", lambda t, n, m: {"run": run} if run else {}), \
                patch.object(wf, "_retry_archive_scheduler_logs", lambda cfg, s, tag: True), \
                patch.object(wf, "_regen_marker_write", lambda *a, **k: calls.append("marker")), \
                patch.object(wf, "do_submit", lambda cfg, t, m, s2, force, gen_first, contcar_cp, tag, submit=True:
                             calls.append(("submit", submit)) or True):
            ok = wf.retry_submit({}, {}, {"name": "WS2"}, s, False, "retry WS2|S5.1")
        return ok, calls

    def test_gen_step_runs_now(self):
        ok, calls = self._retry("gen")
        self.assertTrue(ok)
        self.assertEqual(calls, [("submit", True)])          # 就地跑，不再"只生成、待 start"

    def test_job_step_unchanged(self):
        ok, calls = self._retry(None)
        self.assertEqual(calls, [("submit", False), "marker"])


class LegacyMarkerTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.d = Path(self.td.name) / "step5_dielect_validate"
        self.d.mkdir()

    def test_legacy_marker_is_not_fail(self):
        (self.d / "dielectric_check.json").write_text(json.dumps({"ok": True, "eps_inf": [[3.56]]}))
        done, diag = C.ck_plot(str(self.d), SC)
        self.assertFalse(done)
        self.assertEqual(diag, C.LEGACY_MARKER_DIAG)
        self.assertFalse(C._gen_failed_after(str(self.d), "dielectric_check.json"))

    def test_failed_new_format_still_invalid(self):
        (self.d / "dielectric_check.json").write_text(json.dumps({"ok": False, "status": "failed", "input_files": {}}))
        self.assertEqual(C.ck_plot(str(self.d), SC), (False, "invalid completion marker"))

    def test_regen_failure_falls_back_to_fail(self):
        m = self.d / "dielectric_check.json"
        m.write_text(json.dumps({"ok": True}))
        os.utime(m, (time.time() - 600, time.time() - 600))
        (self.d / ".autozt_gen_failed").write_text("")
        self.assertTrue(C._gen_failed_after(str(self.d), "dielectric_check.json"))

    def test_probe_wiring(self):
        src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
        i = src.index('elif f.get("diag") == LEGACY_MARKER_DIAG')
        self.assertIn('f["stale"] = True', src[i:i + 300])
        self.assertLess(i, src.index("steps.append(f)", i))
        self.assertEqual(wf.step_state({"done": False, "stale": True, "plot_error": False, "job": None,
                                        "has_outcar": False, "has_slurm_out": False, "has_incar": False}, False),
                         ("STALE", "PREP"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
