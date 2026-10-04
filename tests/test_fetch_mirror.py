# -*- coding: utf-8 -*-
"""V162：① 画图步骤的完成标记被归档后，旧 png 不再算"完成"；② 本地回拉镜像不再留着远端已没有的上一代文件，
扇出步骤的子目录一起拉回。

起因（Si_diamond，lineage_check 退出码 1）：
  · S8.1_boltztrap：S3 重算时 boltztrap_crta.json 被归档，可目录里旧 png 让 ck_plot 照样判完成 ->
    表上 19/19，S8.1 一直是 08-30 的结果，S8.2 10-04 已重算；
  · S2.3_hse（扇出 p1of4…p4of4）：fetch 只拉顶层，顶层本来就没有 POSCAR；本地那份是 08-30 旧一代的残留，
    lineage_check 在本地报"用的是旧结构"，集群上其实是新的。
离线：真 tar、本地目录冒充远端，不连超算。
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
    import types
    src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
    mod = types.ModuleType("_collector_defs")
    exec(compile(src[:src.index("\ndef main():")], "_collector_remote.py", "exec"), mod.__dict__)
    return mod


C = _load_collector()


def _touch(p, ts, text="x"):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    os.utime(p, (ts, ts))


class PlotMarkerTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.d = Path(self.td.name) / "step8.1_boltztrap"
        self.now = time.time()
        self.sc = {"check": "plot", "done_marker": "boltztrap_crta.json"}
        _touch(self.d / "mobility_300K.png", self.now - 86400 * 35)

    def test_archived_marker_with_old_png_is_not_done(self):
        tag = time.strftime("%Y%m%d%H%M%S", time.localtime(self.now - 86400 * 2))
        _touch(self.d / ("boltztrap_crta.json.stale-upstream-step3_uniform-" + tag), self.now - 86400 * 35)
        done, diag = C.ck_plot(str(self.d), self.sc)
        self.assertFalse(done)
        self.assertIn("已归档", diag)
        self.assertIsNotNone(C.stale_upstream(str(self.d)))          # -> STALE，auto-advance 会重新生成

    def test_marker_present_or_legacy_png(self):
        self.assertTrue(C.ck_plot(str(self.d), self.sc)[0])          # 没有标记也没归档：照旧按 png 判
        _touch(self.d / ("boltztrap_crta.json.stale-upstream-step8.2_dpt-20261001000000"), self.now - 9e5)
        _touch(self.d / "boltztrap_crta.json", self.now)              # 重新生成过
        self.assertTrue(C.ck_plot(str(self.d), self.sc)[0])


class TarNamesTests(unittest.TestCase):
    def test_gnu_and_bsd_listing(self):
        got = wf._tar_names("./\n./POSCAR\np1of4/POSCAR\n")
        self.assertEqual(got, {"POSCAR", "p1of4/POSCAR"})
        self.assertEqual(wf._tar_names("x INCAR\nx ./KPOINTS\n"), {"INCAR", "KPOINTS"})


class PruneTests(unittest.TestCase):
    def test_only_requested_names_not_received(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            for n in ("POSCAR", "INCAR", "notes.txt"):
                (d / n).write_text("old")
            moved = wf._prune_unreceived(str(d), "./INCAR\n", ["POSCAR", "INCAR", "*.png", "p1of4/POSCAR", None])
            self.assertEqual(moved, ["POSCAR"])
            self.assertFalse((d / "POSCAR").exists())
            self.assertEqual(len(list(d.glob("POSCAR.stale-remote-*"))), 1)
            self.assertTrue((d / "INCAR").is_file())
            self.assertTrue((d / "notes.txt").is_file())               # 不在清单里的本地文件不碰


class FetchFanoutTests(unittest.TestCase):
    """真跑 fetch_material（远端 = 本地目录，bash -c tar）：Si_diamond 的 S2.3_hse 情形。"""

    def test_fanout_subdirs_fetched_and_stale_top_pruned(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            remote = base / "remote" / "step2_bandgap" / "step2.3_hse"
            for k in range(1, 5):
                _touch(remote / ("p%dof4" % k) / "POSCAR", time.time(), "new 2.71784\n")
                _touch(remote / ("p%dof4" % k) / "OUTCAR", time.time(), "done\n")
            dest = base / "result" / "step2_bandgap" / "step2.3_hse"
            _touch(dest / "POSCAR", time.time() - 86400 * 35, "old 2.73313\n")   # 上一代残留
            step = {"name": "step2_bandgap/step2.3_hse", "label": "S2.3_hse", "dir": str(remote),
                    "exists": True, "done": True, "fanout": "*p*of*"}
            m = {"name": "Si_diamond", "result_dir": str(base / "result"), "steps": [step],
                 "fetch_files": ["POSCAR", "OUTCAR", "INCAR"],
                 "_seg": {"steps_cfg": [{"name": step["name"], "fanout": "*p*of*"}]}}
            with patch.object(wf, "fetch_provenance_dir", create=True), \
                    patch("autozt.log_action"), patch("autozt.fetch_provenance_dir", create=True):
                self.assertTrue(wf.fetch_material({}, m, quiet=True))
            self.assertFalse((dest / "POSCAR").exists())
            self.assertEqual(len(list(dest.glob("POSCAR.stale-remote-*"))), 1)
            for k in range(1, 5):
                self.assertEqual((dest / ("p%dof4" % k) / "POSCAR").read_text(), "new 2.71784\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
