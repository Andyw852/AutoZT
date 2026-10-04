# -*- coding: utf-8 -*-
"""V163：上游失效归档的文件，必须真能让采集端判"未完成"（否则失效不生效，表上照样 OK）。

起因：V162 的 S8.1 —— invalidate_downstream 归档 boltztrap_crta.json，可 ck_plot 还认旧 png，失效形同虚设。
这里对 ke-dft-cpu / zt-dft-cpu 里每个非 dir 的失效目标逐个造现场（标记 + 旧 png + slurm 输出），
先确认判"完成"，再按 invalidate_downstream 的命名归档 DONE_MARKERS，确认判"未完成"且识别为 STALE。
以后有人改 skill.yaml 的 check/done_marker 或 DONE_MARKERS，对不上就在这里挂。
离线。
"""
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _p in (str(ROOT), str(KE), str(ROOT / "skill" / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


def _load_collector():
    import types
    src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
    mod = types.ModuleType("_collector_defs")
    exec(compile(src[:src.index("\ndef main():")], "_collector_remote.py", "exec"), mod.__dict__)
    return mod


C = _load_collector()


def _steps():
    out = {}

    def walk(o):
        if isinstance(o, dict):
            if "name" in o and ("check" in o or "gen" in o):
                out.setdefault(o["name"], o)
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    for f in (KE / "skill.yaml", ROOT / "skill" / "zt-dft-cpu" / "skill.yaml"):
        if f.is_file():
            walk(yaml.safe_load(f.read_text(encoding="utf-8")))
    return out


class InvalidationTakesEffectTests(unittest.TestCase):
    def test_archived_markers_make_step_not_done(self):
        steps = _steps()
        targets = sorted({d for lst in kc.DOWNSTREAM.values() for d, mode in lst if mode != "dir"})
        self.assertTrue(targets)
        now = time.time()
        tag = "stale-upstream-up-" + time.strftime("%Y%m%d%H%M%S", time.localtime(now - 60))
        bad = []
        for name in targets:
            sc = steps.get(name)
            if not sc:
                continue
            ck = C.CHECKERS.get(sc.get("check", "outcar"))
            markers = kc.DONE_MARKERS.get(name, ())
            if ck is None or not markers:
                bad.append("%s：check=%s / DONE_MARKERS=%s" % (name, sc.get("check"), markers))
                continue
            with tempfile.TemporaryDirectory() as td:
                d = Path(td) / name
                d.mkdir(parents=True)
                body = (sc.get("marker") or ":").split(":", 1)[1] or "x"
                for f in set(markers) | {(sc.get("marker") or "").split(":", 1)[0], sc.get("done_marker")} - {"", None}:
                    (d / f).write_text(body)
                if sc.get("strict_done_marker"):                                # V136：status=done + 输入 sha256
                    import hashlib
                    import json
                    src = Path(td) / "step5_dielect" / "OUTCAR"
                    src.parent.mkdir(parents=True, exist_ok=True)
                    src.write_text("outcar")
                    (d / sc["done_marker"]).write_text(json.dumps({"status": "done", "input_files": {
                        "step5_dielect/OUTCAR": hashlib.sha256(b"outcar").hexdigest()}}))
                for leftover in ("old_300K.png", "slurm-1.out"):              # 旧图、旧作业输出
                    (d / leftover).write_text("old")
                for f in d.iterdir():
                    os.utime(f, (now - 86400, now - 86400))
                if not ck(str(d), sc)[0]:
                    bad.append("%s：造好的完成现场判不成完成（测试现场不对）" % name)
                    continue
                for f in markers:
                    if (d / f).exists():
                        (d / f).rename(d / (f + "." + tag))
                done, diag = ck(str(d), sc)
                if done:
                    bad.append("%s：归档 %s 后仍判完成（%s）" % (name, "/".join(markers), diag))
                elif not C.stale_upstream(str(d)):
                    bad.append("%s：归档后未识别为 STALE" % name)
        self.assertEqual(bad, [], "\n  ".join([""] + bad))


if __name__ == "__main__":
    unittest.main(verbosity=2)
