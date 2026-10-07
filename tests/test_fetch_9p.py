# -*- coding: utf-8 -*-
"""V179：本地结果目录在 9p / drvfs 挂载上，tar 只是设不了时间戳 / 权限（内容已落盘）时，fetch 不再整步失败。

起因（CrS₂ S8.4，2026-10-07）：autozt fetch 在 S1_opt 第一步就报 tar "Cannot utime: Operation not permitted" 并 abort，
后面的步骤一个都拉不回；transport.json 等只好手工 scp。
离线：真 tar、本地目录冒充远端；"9p" 用包一层 subprocess.run 模拟（真解包，再补上 GNU tar 在 9p 上打的那几行、退出码 2）。
"""
import io
import subprocess
import sys
import tempfile
import time
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from autozt import workflow as wf  # noqa: E402
from autozt import prov  # noqa: E402

GNU_9P = ("tar: ./transport.json: Cannot utime: Operation not permitted\n"
          "tar: .: Cannot utime: Operation not permitted\n"
          "tar: Exiting with failure status due to previous errors\n")
_REAL_RUN = subprocess.run


def _fake_9p(extra):
    def run(args, *a, **kw):
        r = _REAL_RUN(args, *a, **kw)
        if extra and list(args[:2]) in (["tar", "xvf"], ["tar", "xf"]) and r.returncode == 0:
            return subprocess.CompletedProcess(args, 2, r.stdout, (r.stderr or "") + extra)
        return r
    return run


class HelperTests(unittest.TestCase):
    def test_meta_only(self):
        self.assertTrue(wf.tar_meta_only_errors(GNU_9P))
        self.assertTrue(wf.tar_meta_only_errors("x transport.json\ntar: transport.json: Can't restore time\n"))  # bsdtar
        self.assertTrue(wf.tar_meta_only_errors("tar: ./a: Cannot change mode to rw-r--r--: Operation not permitted\n"))
        self.assertFalse(wf.tar_meta_only_errors(""))
        self.assertFalse(wf.tar_meta_only_errors("tar: Exiting with failure status due to previous errors\n"))
        self.assertFalse(wf.tar_meta_only_errors(
            "tar: ./transport.json: Cannot write: No space left on device\n" + GNU_9P))      # 真错误混在里面
        self.assertFalse(wf.tar_meta_only_errors("tar: This does not look like a tar archive\n"))


def _material(base):
    steps, names = [], ("step1_opt", "step8.4_amset2d")
    for n in names:
        d = base / "remote" / n
        d.mkdir(parents=True)
        (d / "transport.json").write_text('{"step": "%s"}\n' % n)
        steps.append({"name": n, "label": n, "dir": str(d), "exists": True, "done": True})
    return {"name": "CrS2_hex", "result_dir": str(base / "result"), "steps": steps,
            "fetch_files": ["transport.json"], "_seg": {"steps_cfg": []}}


class FetchTests(unittest.TestCase):
    def _fetch(self, extra):
        td = tempfile.TemporaryDirectory()
        self.addCleanup(td.cleanup)
        base = Path(td.name)
        m = _material(base)
        buf = io.StringIO()
        with patch.object(wf.subprocess, "run", side_effect=_fake_9p(extra)), \
                patch("autozt.log_action"), patch("autozt.fetch_provenance_dir", create=True), \
                redirect_stdout(buf):
            ok = wf.fetch_material({}, m, quiet=False)
        return ok, buf.getvalue(), base

    def test_9p_utime_does_not_abort(self):
        ok, out, base = self._fetch(GNU_9P)
        self.assertTrue(ok, out)
        for n in ("step1_opt", "step8.4_amset2d"):                    # 第一步之后的步骤也拉回来了
            self.assertEqual((base / "result" / n / "transport.json").read_text(), '{"step": "%s"}\n' % n)
        self.assertEqual(out.count("不允许设置时间戳"), 1)              # 只告警一次
        self.assertIn("已拉回 2 个步骤", out)

    def test_real_error_still_fails(self):
        ok, out, _ = self._fetch("tar: ./transport.json: Cannot write: No space left on device\n" + GNU_9P)
        self.assertFalse(ok)
        self.assertIn("fetch step1_opt 失败", out)

    def test_clean_fetch_has_no_warning(self):
        ok, out, _ = self._fetch("")
        self.assertTrue(ok)
        self.assertNotIn("不允许设置时间戳", out)


class OtherSitesTests(unittest.TestCase):
    def test_provenance_fetch_tolerates_9p(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            (base / "remote" / prov.PROV_DIR).mkdir(parents=True)
            (base / "remote" / prov.PROV_DIR / "step1_opt.json").write_text("{}")
            (base / "remote" / "step1_opt").mkdir()
            m = {"name": "CrS2_hex", "result_dir": str(base / "result"),
                 "steps": [{"name": "step1_opt", "dir": str(base / "remote" / "step1_opt")}]}
            with patch("subprocess.run", side_effect=_fake_9p(GNU_9P)):        # prov 在函数里 import subprocess
                self.assertTrue(prov.fetch_provenance_dir({}, m, quiet=True))
            with patch("subprocess.run",
                       side_effect=_fake_9p("tar: ./x: Cannot open: Permission denied\n" + GNU_9P)):
                self.assertFalse(prov.fetch_provenance_dir({}, m, quiet=True))

    def test_fanout_input_site_wired(self):
        src = (ROOT / "autozt" / "workflow.py").read_text(encoding="utf-8")
        i = src.index('return False, "扇出输入回拉失败')
        self.assertIn("tar_meta_only_errors(p2.stderr)", src[i - 200:i])


if __name__ == "__main__":
    unittest.main(verbosity=2)
