# -*- coding: utf-8 -*-
"""V170：显式 -p X -j S start 不再重交已完成（OK）的步骤。

起因（CrS₂ S8.4，2026-10-06）：S8.4 跑了 4h41m 正常完成，agent 想拉结果又执行了一次 start ——
这条路不看 OK，按现有输入原样重交，作业开头 rm -f transport.json 把结果删了，步骤变 FAIL，
只能手工从 transport_*.json 拷回。批量 start / step init / rerun 早就跳过 OK。
例外：retry 刚重新生成过输入、还没被 start 交过（指纹一致）—— retry→start 照常放行一次。
离线。
"""
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import autozt                      # noqa: E402
from autozt import workflow as wf  # noqa: E402


def _mat(tmpdir, kind="OK"):
    m = {"name": "CrS2_hex", "tt": "ke-dft-cpu", "lpath": tmpdir, "ps": {"dir": tmpdir},
         "host_eff": "mock", "steps": []}
    s = {"name": "step8.4_amset2d", "label": "S8.4_amset2d", "kind": kind,
         "dir": "/remote/step8.4_amset2d", "job": None, "has_incar": True, "has_outcar": False,
         "done": kind == "OK", "exists": True, "submit": "submit.sh", "_host": "mock"}
    m["steps"] = [s]
    return m, s


class StartOkStepTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.m, self.s = _mat(self.td.name)

    def _marker(self, fp="fpA", **extra):
        ent = {"version": 1, "time": "x", "ts": time.time(), "fingerprint": fp, "files": ["submit.sh"],
               "material": self.m["name"], "tt": self.m["tt"], "step": self.s["name"],
               "label": self.s["label"], "grants": []}
        ent.update(extra)
        wf._regen_marker_save(self.m, {wf._regen_key(self.m, self.s): ent})

    def _start(self, force=False, fp="fpA"):
        with patch.object(autozt, "find_material", return_value=({"key": "ke-dft-cpu"}, self.m)), \
                patch.object(wf, "guard_predecessors", return_value=True), \
                patch.object(wf, "_regen_input_fingerprint", return_value=(fp, ["submit.sh"])), \
                patch.object(autozt, "step_cfg", return_value={}), \
                patch.object(wf, "_gen_step_input", side_effect=AssertionError("OK 步骤不能重新生成输入")), \
                patch.object(wf, "do_submit", return_value=True) as ds, \
                patch("builtins.print") as pr:
            rc = wf.cmd_start({}, {"types": []}, "CrS2_hex", "S8.4_amset2d", force)
        msg = " ".join(str(c.args[0]) for c in pr.call_args_list if c.args)
        return rc, ds.called, msg

    def test_ok_step_is_not_resubmitted(self):
        rc, called, msg = self._start()
        self.assertEqual(rc, 0)
        self.assertFalse(called)
        self.assertIn("已完成（OK）", msg)
        self.assertIn("fetch", msg)

    def test_force_still_resubmits(self):
        rc, called, _ = self._start(force=True)
        self.assertEqual(rc, 0)
        self.assertTrue(called)

    def test_retry_then_start_allowed_once(self):
        self._marker("fpA")
        rc, called, msg = self._start()
        self.assertTrue(called)
        self.assertIn("retry", msg)
        ent = wf._regen_marker_load(self.m)[wf._regen_key(self.m, self.s)]
        self.assertTrue(ent.get("consumed_ts"))
        self.assertEqual(ent.get("grants"), [])                 # 不占 FAIL 的放行额度
        rc, called, _ = self._start()                           # 第二次：标记已用过
        self.assertFalse(called)

    def test_stale_marker_fingerprint_denies(self):
        self._marker("fpA")
        rc, called, _ = self._start(fp="fpB")                   # retry 之后输入又被改过
        self.assertFalse(called)

    def test_fail_path_consumes_marker_too(self):
        self.m, self.s = _mat(self.td.name, kind="FAIL")
        self._marker("fpA")
        rc, called, _ = self._start()
        self.assertTrue(called)
        self.s["kind"], self.s["done"] = "OK", True               # 跑完之后再误 start 一次
        rc, called, _ = self._start()
        self.assertFalse(called)


if __name__ == "__main__":
    unittest.main(verbosity=2)
