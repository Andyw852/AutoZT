# -*- coding: utf-8 -*-
"""V171：同一步骤的 gen 不许并发（远端 flock），被挡下的那次不算 gen 失败。

起因（CrSe₂ S7.1，2026-10-06）：agent 第一次 retry 本地 2 分钟超时、远端 gen 照跑；第二次 retry 又起一个。
前一个最后 os.replace 把后一个刚写出的 deformation.h5 挪进 step7b，后一个读 h5 报 FileNotFoundError。
本地跑：remote_gen 拼出的远端命令原样交给本机 bash（不碰集群）。
"""
import os
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import autozt  # noqa: E402
from autozt import workflow as wf  # noqa: E402


def _local_run_remote(cfg, shell_line, host="__default__", use_stdin=False):
    r = subprocess.run(["bash", "-s"], input=shell_line, capture_output=True, text=True, timeout=60)
    return r.returncode, (r.stdout + r.stderr).strip()


class GenLockTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.m = {"name": "CrSe2_hex", "path": self.td.name, "tt": "ke-dft-cpu", "lpath": None}
        self.t = {"key": "ke-dft-cpu", "root": self.td.name}

    def _gen(self, cmd, out):
        sc = {"gen": cmd, "gen_need": []}
        with patch.object(autozt, "step_cfg", return_value=sc), \
                patch.object(autozt, "run_remote", side_effect=_local_run_remote), \
                patch.object(autozt, "provenance_enabled", return_value=False), \
                patch.object(wf, "common_autodeps_enabled", return_value=False), \
                patch.object(wf, "resolve_cores", return_value=None):
            out.append(wf.remote_gen({}, self.t, self.m, "step7_deform/step7b_read"))

    def test_second_concurrent_gen_is_refused(self):
        if subprocess.run(["bash", "-c", "command -v flock"], capture_output=True).returncode != 0:
            self.skipTest("本机没有 flock")
        a, b = [], []
        th = threading.Thread(target=self._gen, args=("sleep 3; echo A-done", a))
        th.start()
        time.sleep(1.0)
        self._gen("echo B-ran", b)                      # 第一个还在跑
        th.join()
        self.assertTrue(a[0][0], a)
        self.assertIn("A-done", a[0][1])
        self.assertFalse(b[0][0])
        self.assertTrue(wf._gen_busy(b[0][1]))
        self.assertNotIn("B-ran", b[0][1])               # 被挡下的那次一行都没执行
        c = []
        self._gen("echo C-ran", c)                       # 前一个跑完，锁随进程退出释放
        self.assertTrue(c[0][0], c)
        self.assertIn("C-ran", c[0][1])

    def test_other_steps_not_blocked(self):
        if subprocess.run(["bash", "-c", "command -v flock"], capture_output=True).returncode != 0:
            self.skipTest("本机没有 flock")
        a, b = [], []
        th = threading.Thread(target=self._gen, args=("sleep 3", a))
        th.start()
        time.sleep(1.0)
        sc = {"gen": "echo other", "gen_need": []}
        with patch.object(autozt, "step_cfg", return_value=sc), \
                patch.object(autozt, "run_remote", side_effect=_local_run_remote), \
                patch.object(autozt, "provenance_enabled", return_value=False), \
                patch.object(wf, "common_autodeps_enabled", return_value=False), \
                patch.object(wf, "resolve_cores", return_value=None):
            b.append(wf.remote_gen({}, self.t, self.m, "step8.4_amset2d"))
        th.join()
        self.assertTrue(b[0][0], b)

    def test_busy_is_not_marked_as_gen_failure(self):
        s = {"name": "step7_deform/step7b_read", "label": "S7.1_read", "dir": self.td.name, "stale": True,
             "job": None}
        with patch.object(wf, "remote_gen", return_value=(False, "TF_GEN_BUSY 同一步骤已有一个 gen 在跑")), \
                patch.object(wf, "_mark_gen_failed") as mk, patch("builtins.print"):
            self.assertFalse(wf.do_run_gen_step({}, self.t, self.m, s, "retry x"))
        mk.assert_not_called()
        with patch.object(wf, "remote_gen", return_value=(False, "Traceback ...")), \
                patch.object(wf, "_mark_gen_failed") as mk, patch("builtins.print"):
            self.assertFalse(wf.do_run_gen_step({}, self.t, self.m, s, "retry x"))
        mk.assert_called_once()


if __name__ == "__main__":
    unittest.main(verbosity=2)
