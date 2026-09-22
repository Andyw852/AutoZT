# -*- coding: utf-8 -*-
"""方案A（瘦身版）离线回归：retry 重生成标记 + cmd_start FAIL 免 -f 豁免。

全部离线：不打超算、不动主树、不读真实项目配置；指纹用 mock 或直接 patch
_regen_input_fingerprint。覆盖任务要求：标记缺失/指纹不符/24h 第 3 次不放行、
满足条件放行且写审计、retry scancel+regenerate 语义未变、do_submit 未被包装。
"""
import ast
import base64
import json
import os
import re
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


def _mat(tmpdir, tt="band-dft-cpu", name="C60/qHPC60"):
    m = {"name": name, "tt": tt, "lpath": tmpdir, "ps": {"dir": tmpdir},
         "host_eff": "mock", "steps": []}
    s = {"name": "step1_opt", "label": "S1_opt", "kind": "FAIL",
         "dir": "/remote/step1_opt", "job": None, "has_incar": True,
         "has_outcar": True, "done": False, "exists": True,
         "submit": "submit.sh", "_host": "mock"}
    m["steps"] = [s]
    return m, s


class RegenMarkerTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.m, self.s = _mat(self.td.name)

    # ---------- 放行条件 ----------
    def test_marker_missing_denies(self):
        with patch.object(wf, "_regen_input_fingerprint",
                          return_value=("fpA", ["INCAR"])):
            self.assertFalse(wf.regen_allow_without_force({}, {}, self.m, self.s))
        self.assertFalse(os.path.exists(wf._regen_path(self.m, wf.REGEN_AUDIT)))

    def test_fingerprint_mismatch_denies(self):
        self._write_marker("fpA")
        with patch.object(wf, "_regen_input_fingerprint",
                          return_value=("fpB", ["INCAR"])):
            self.assertFalse(wf.regen_allow_without_force({}, {}, self.m, self.s))
        ent = wf._regen_marker_load(self.m)[wf._regen_key(self.m, self.s)]
        self.assertEqual(ent.get("grants"), [])          # 不消耗额度

    def test_fingerprint_unavailable_denies(self):
        self._write_marker("fpA")
        with patch.object(wf, "_regen_input_fingerprint",
                          return_value=(None, [])):
            self.assertFalse(wf.regen_allow_without_force({}, {}, self.m, self.s))

    def test_first_two_allow_then_third_within_24h_denied(self):
        self._write_marker("fpA")
        with patch.object(wf, "_regen_input_fingerprint",
                          return_value=("fpA", ["INCAR"])):
            self.assertTrue(wf.regen_allow_without_force({}, {}, self.m, self.s))
            self.assertTrue(wf.regen_allow_without_force({}, {}, self.m, self.s))
            self.assertFalse(wf.regen_allow_without_force({}, {}, self.m, self.s))
        audit = self._audit_lines()
        self.assertEqual([a["n"] for a in audit], [1, 2])
        for a in audit:
            for k in ("time", "material", "tt", "step", "fingerprint", "n"):
                self.assertIn(k, a)
            self.assertEqual(a["fingerprint"], "fpA")
        ent = wf._regen_marker_load(self.m)[wf._regen_key(self.m, self.s)]
        self.assertEqual(len(ent["grants"]), 2)          # 第 3 次未记账

    def test_grants_older_than_24h_do_not_count(self):
        self._write_marker("fpA")
        key = wf._regen_key(self.m, self.s)
        marks = wf._regen_marker_load(self.m)
        marks[key]["grants"] = [time.time() - 25 * 3600, time.time() - 26 * 3600]
        wf._regen_marker_save(self.m, marks)
        with patch.object(wf, "_regen_input_fingerprint",
                          return_value=("fpA", ["INCAR"])):
            self.assertTrue(wf.regen_allow_without_force({}, {}, self.m, self.s))
            self.assertTrue(wf.regen_allow_without_force({}, {}, self.m, self.s))
            self.assertFalse(wf.regen_allow_without_force({}, {}, self.m, self.s))

    def test_grants_survive_new_retry_marker(self):
        self._write_marker("fpA")
        key = wf._regen_key(self.m, self.s)
        marks = wf._regen_marker_load(self.m)
        marks[key]["grants"] = [time.time()]
        wf._regen_marker_save(self.m, marks)
        with patch.object(wf, "_regen_input_fingerprint",
                          return_value=("fpA", ["INCAR"])), \
             patch.object(autozt, "step_cfg", return_value={}):
            self.assertTrue(wf._regen_marker_write(
                {}, {"key": "band-dft-cpu"}, self.m, self.s, "retry x"))
        self.assertEqual(
            len(wf._regen_marker_load(self.m)[key]["grants"]), 1)

    def test_marker_write_records_identity_and_files(self):
        with patch.object(wf, "_regen_input_fingerprint",
                          return_value=("FP1", ["INCAR", "POSCAR"])), \
             patch.object(autozt, "step_cfg", return_value={}):
            self.assertTrue(wf._regen_marker_write(
                {}, {}, self.m, self.s, "retry t"))
        ent = wf._regen_marker_load(self.m)[wf._regen_key(self.m, self.s)]
        self.assertEqual(ent["fingerprint"], "FP1")
        self.assertEqual(ent["material"], "C60/qHPC60")
        self.assertEqual(ent["tt"], "band-dft-cpu")
        self.assertEqual(ent["step"], "step1_opt")
        self.assertEqual(ent["label"], "S1_opt")
        self.assertIn("time", ent)

    def test_marker_write_skips_run_gen(self):
        with patch.object(autozt, "step_cfg", return_value={"run": "gen"}):
            self.assertFalse(wf._regen_marker_write(
                {}, {}, self.m, self.s, "retry x"))
        self.assertFalse(os.path.exists(wf._regen_path(self.m, wf.REGEN_MARK)))

    # ---------- 指纹 ----------
    def test_fingerprint_stable_and_content_sensitive(self):
        hashes = {"INCAR": "a" * 64, "POSCAR": "b" * 64}

        def fake_run(cfg, cmd, **kw):
            b64 = re.search(r"echo (\S+) \| base64", cmd).group(1)
            script = base64.b64decode(b64).decode()
            out = ""
            for name, h in hashes.items():
                if name in script:
                    out += "%s  %s\n" % (h, name)
            return 0, out

        with patch.object(autozt, "run_remote", side_effect=fake_run):
            fp1, files = wf._regen_input_fingerprint({}, self.m, self.s)
            fp2, _ = wf._regen_input_fingerprint({}, self.m, self.s)
        self.assertEqual(fp1, fp2)
        self.assertEqual(sorted(files), ["INCAR", "POSCAR"])
        hashes["POSCAR"] = "c" * 64
        with patch.object(autozt, "run_remote", side_effect=fake_run):
            fp3, _ = wf._regen_input_fingerprint({}, self.m, self.s)
        self.assertNotEqual(fp1, fp3)

    # ---------- retry_submit 语义 ----------
    def test_retry_scancel_regenerate_then_marker(self):
        order = []
        s_job = dict(self.s)
        s_job["job"] = {"id": "999"}
        with patch.object(wf, "kill_if_queued",
                          side_effect=lambda *a: (order.append("scancel") or True)), \
             patch.object(wf, "_retry_archive_scheduler_logs",
                          side_effect=lambda *a: (order.append("archive") or True)), \
             patch.object(wf, "do_submit",
                          side_effect=lambda *a, **k: (order.append(("gen", k)) or True)), \
             patch.object(autozt, "step_cfg",
                          return_value={"contcar_to_poscar": True}), \
             patch.object(wf, "_regen_marker_write",
                          side_effect=lambda *a: order.append("marker")):
            ok = wf.retry_submit({}, {}, self.m, s_job, False, "retry x")
        self.assertTrue(ok)
        self.assertEqual(order[0], "scancel")            # 先 scancel
        self.assertEqual(order[1], "archive")            # 归档旧调度日志
        self.assertEqual(order[2][0], "gen")             # 再重新生成
        self.assertTrue(order[2][1]["gen_first"])
        self.assertFalse(order[2][1]["submit"])
        self.assertEqual(order[3], "marker")             # 成功后写标记

    def test_retry_gen_failure_writes_no_marker(self):
        s_job = dict(self.s)
        s_job["job"] = {"id": "999"}
        with patch.object(wf, "kill_if_queued", return_value=True), \
             patch.object(wf, "_retry_archive_scheduler_logs", return_value=True), \
             patch.object(wf, "do_submit", return_value=False), \
             patch.object(autozt, "step_cfg", return_value={}), \
             patch.object(wf, "_regen_marker_write") as mw:
            self.assertFalse(wf.retry_submit({}, {}, self.m, s_job, False, "retry x"))
        mw.assert_not_called()

    def test_retry_archive_failure_stops_before_gen(self):
        with patch.object(wf, "_retry_archive_scheduler_logs", return_value=False), \
             patch.object(wf, "do_submit") as ds:
            self.assertFalse(wf.retry_submit({}, {}, self.m, self.s, False, "retry x"))
        ds.assert_not_called()

    # ---------- cmd_start FAIL 分支集成 ----------
    def test_cmd_start_denies_fail_without_marker(self):
        with patch.object(autozt, "find_material",
                          return_value=({"key": "band-dft-cpu"}, self.m)), \
             patch.object(wf, "guard_predecessors", return_value=True), \
             patch.object(wf, "do_submit") as ds:
            rc = wf.cmd_start({}, {"types": []}, "C60/qHPC60", "S1_opt", False)
        self.assertEqual(rc, 1)
        ds.assert_not_called()

    def test_cmd_start_allows_fail_with_valid_marker(self):
        self._write_marker("fpA")
        seen = {}

        def fake_submit(cfg, t, m, s, force, **k):
            seen["force"] = force
            seen["gen_first"] = k.get("gen_first")
            return True

        with patch.object(autozt, "find_material",
                          return_value=({"key": "band-dft-cpu"}, self.m)), \
             patch.object(wf, "guard_predecessors", return_value=True), \
             patch.object(wf, "_regen_input_fingerprint",
                          return_value=("fpA", ["INCAR"])), \
             patch.object(autozt, "step_cfg", return_value={}), \
             patch.object(wf, "do_submit", side_effect=fake_submit):
            rc = wf.cmd_start({}, {"types": []}, "C60/qHPC60", "S1_opt", False)
        self.assertEqual(rc, 0)
        self.assertFalse(seen["force"])          # -f 语义不变：仍以 force=False 提交
        self.assertFalse(seen["gen_first"])
        self.assertEqual([a["n"] for a in self._audit_lines()], [1])

    def test_cmd_start_force_still_bypasses(self):
        with patch.object(autozt, "find_material", return_value=({}, self.m)), \
             patch.object(wf, "guard_predecessors", return_value=True), \
             patch.object(autozt, "step_cfg", return_value={}), \
             patch.object(wf, "do_submit", return_value=True) as ds:
            rc = wf.cmd_start({}, {"types": []}, "C60/qHPC60", "S1_opt", True)
        self.assertEqual(rc, 0)
        self.assertTrue(ds.called)

    # ---------- do_submit 未被包装 ----------
    def test_do_submit_not_wrapped_or_reassigned(self):
        src = (ROOT / "autozt" / "workflow.py").read_text(encoding="utf-8")
        tree = ast.parse(src)
        defs = [n for n in ast.walk(tree)
                if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
                and n.name == "do_submit"]
        self.assertEqual(len(defs), 1)
        self.assertEqual(defs[0].decorator_list, [])
        stores = [n for n in ast.walk(tree)
                  if isinstance(n, ast.Name) and n.id == "do_submit"
                  and isinstance(n.ctx, ast.Store)]
        self.assertEqual(stores, [])
        self.assertEqual(wf.do_submit.__code__.co_name, "do_submit")
        self.assertIs(wf.do_submit, autozt.do_submit)

    # ---------- helpers ----------
    def _write_marker(self, fp):
        key = wf._regen_key(self.m, self.s)
        wf._regen_marker_save(self.m, {key: {
            "version": 1, "time": "x", "ts": time.time(), "fingerprint": fp,
            "files": ["INCAR"], "material": self.m["name"], "tt": self.m["tt"],
            "step": self.s["name"], "label": self.s["label"], "grants": []}})

    def _audit_lines(self):
        p = wf._regen_path(self.m, wf.REGEN_AUDIT)
        if not p or not os.path.isfile(p):
            return []
        with open(p, encoding="utf-8") as f:
            return [json.loads(x) for x in f if x.strip()]


if __name__ == "__main__":
    unittest.main(verbosity=2)
