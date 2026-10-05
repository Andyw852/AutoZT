# -*- coding: utf-8 -*-
"""V168：扇出步骤有子作业在跑时，start 只补交没交上的子目录；-f 杀掉的子作业要重交；部分提交要记账。

起因（CrS₂ S7_deform，2026-10-06）：retry 后 start 只交上 5/10 个子目录，再 start 报"已有作业，先 stop 或加 -f"。
  · 没有"补交"这条路：扇出步骤任一子作业在跑，kill_if_queued 就整步拒绝；
  · -f 会 scancel 全部在跑的子作业，可提交用的 fan_todo 是杀之前采的（不含它们）-> 杀了不交，
    fan_todo 为空时还被"待补清单为空"拒交，一个也没交上；
  · 按作业名去重（squeue 不支持 %Z）时，同一步任一子作业在跑就整批拒交；
  · sbatch 交到一半被拒（QOS 上限）时，已交上的 jobid 不记账、不说剩下的怎么补。
离线：run_remote / sbatch 全部打桩或用假调度器。
"""
import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import autozt  # noqa: E402
from autozt import workflow as wf  # noqa: E402


def _step(**over):
    s = {"name": "step7_deform", "label": "S7_deform", "dir": "/remote/step7_deform",
         "fanout": "*deform*", "has_incar": True, "_host": "__default__", "submit": "submit.sh",
         "job": {"id": "3919547", "state": "R", "info": "0/10 5R 0PD"},
         "fan_jobids": ["3919547", "3919548", "3919549", "3919550", "3919551"],
         "fan_running": ["deform-01", "deform-02", "deform-03", "deform-04", "deform-09"],
         "fan_todo": ["deform-05", "deform-06", "deform-07", "deform-08", "undeformed"],
         "fan_done": 0}
    s.update(over)
    return s


M = {"name": "CrS2/ke-dft-cpu", "tt": "ke-dft-cpu"}


def _common_patches(guarded):
    return [patch.object(autozt, "step_cfg", return_value={"run": "sbatch"}),
            patch.object(autozt, "log_action"),
            patch.object(autozt, "agent_event"),
            patch.object(wf, "_local_env_ready", return_value=True),
            patch.object(wf, "_fetch_stamp_clear"),
            patch.object(wf, "_scancel_clear"),
            patch.object(wf, "remote_gen", side_effect=AssertionError("补交不能重新生成输入")),
            patch.object(wf, "_sbatch_guarded", side_effect=guarded)]


class _Ctx:
    def __init__(self, ps):
        self.ps = ps

    def __enter__(self):
        return [p.__enter__() for p in self.ps]

    def __exit__(self, *a):
        for p in reversed(self.ps):
            p.__exit__(*a)


class TopUpTests(unittest.TestCase):
    def test_start_tops_up_without_force_and_without_killing(self):
        calls = []

        def guarded(cfg, sdir, host, jobname=None, fanout=None, only=None, force=False,
                    submit=None, empty_fanout_manifest=None):
            calls.append(list(only or []))
            return True, "", "3919560,3919561,3919562,3919563,3919564"

        s = _step()
        with _Ctx(_common_patches(guarded)), \
                patch.object(wf, "remote_scancel", side_effect=AssertionError("补交不能 scancel 在跑的")):
            ok = wf.do_submit({}, {}, M, s, False, False, False, "start CrS2|S7")
        self.assertTrue(ok)
        self.assertEqual(calls, [["deform-05", "deform-06", "deform-07", "deform-08", "undeformed"]])

    def test_nothing_to_top_up_explains_instead_of_suggesting_force(self):
        s = _step(fan_todo=[], fan_done=5)
        with _Ctx(_common_patches(lambda *a, **k: self.fail("不该提交"))), \
                patch("builtins.print") as pr:
            ok = wf.do_submit({}, {}, M, s, False, False, False, "start CrS2|S7")
        self.assertFalse(ok)
        msg = " ".join(str(c.args[0]) for c in pr.call_args_list if c.args)
        self.assertIn("已判完成", msg)
        self.assertIn("retry", msg)

    def test_force_resubmits_the_killed_subdirs_too(self):
        calls = []

        def guarded(cfg, sdir, host, jobname=None, fanout=None, only=None, force=False,
                    submit=None, empty_fanout_manifest=None):
            calls.append(list(only or []))
            return True, "", "1,2,3,4,5"

        s = _step(fan_todo=[])                   # CrS₂：其余 5 个被旧 OUTCAR 判完成，fan_todo 为空
        with _Ctx(_common_patches(guarded)), \
                patch.object(wf, "remote_scancel", return_value=(True, "")) as sc, \
                patch.object(wf, "_remote_submit_preflight", return_value=(True, "")):
            ok = wf.do_submit({}, {}, M, s, True, False, False, "start -f CrS2|S7")
        self.assertTrue(ok)
        sc.assert_called_once()
        self.assertEqual(calls, [["deform-01", "deform-02", "deform-03", "deform-04", "deform-09"]])

    def test_gen_first_never_tops_up(self):
        s = _step()
        self.assertFalse(wf._fanout_topup(s, False, True, "auto"))   # auto-advance / retry 会重新生成
        self.assertFalse(wf._fanout_topup(s, True, False, "-f"))
        self.assertFalse(wf._fanout_topup(dict(s, fanout=None), False, False, "plain"))

    def test_partial_sbatch_is_accounted(self):
        def guarded(cfg, sdir, host, jobname=None, fanout=None, only=None, force=False,
                    submit=None, empty_fanout_manifest=None):
            return False, "sbatch 失败：QOSMaxSubmitJobPerUserLimit", "3919547,3919548"

        s = _step(job=None, fan_jobids=[], fan_running=[])
        with _Ctx(_common_patches(guarded)), \
                patch.object(autozt, "ledger_append") as la, \
                patch("builtins.print") as pr:
            ok = wf.do_submit({}, {}, M, s, False, False, False, "start CrS2|S7")
        self.assertFalse(ok)
        la.assert_called_once()
        self.assertEqual(la.call_args.args[3], "3919547,3919548")
        msg = " ".join(str(c.args[0]) for c in pr.call_args_list if c.args)
        self.assertIn("已交上 2 个", msg)
        self.assertIn("不要加 -f", msg)


def _load_fake_sched():
    spec = importlib.util.spec_from_file_location("_sbatch_dedup", ROOT / "test" / "test_sbatch_dedup.py")
    mod = importlib.util.module_from_spec(spec)
    cwd = os.getcwd()
    try:
        spec.loader.exec_module(mod)
    finally:
        os.chdir(cwd)
    return mod


class GuardNameModeTests(unittest.TestCase):
    """squeue 不支持 %Z 时按作业名去重：要逐子目录比，不能同一步任一子作业在跑就整批拒交。"""

    def test_name_mode_per_subdir(self):
        D = _load_fake_sched()
        sched = D._FAKE_SCHED.replace(
            "    use_wd = '%Z' in fmt\n",
            "    use_wd = '%Z' in fmt\n"
            "    if use_wd and os.environ.get('AUTOZT_FAKE_NO_WD') == '1': sys.exit(1)\n")
        self.assertNotEqual(sched, D._FAKE_SCHED)
        os.makedirs(ROOT / "tmp", exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(ROOT / "tmp")) as tmp:
            bin_dir, state = D._make_scheduler(tmp)
            with open(os.path.join(bin_dir, "fake_sched.py"), "w") as f:
                f.write(sched)
            s, _ = D._mk_step(tmp, fanout=True, subdirs=("d1", "d2", "d3"))
            mock = D._make_mock(bin_dir, state, {"AUTOZT_FAKE_NO_WD": 1})
            with patch.object(autozt, "run_remote", side_effect=mock):
                ok1, out1, _ = wf._sbatch_guarded({}, s["dir"], "fakehost", jobname="J", fanout="d*", only=["d1"])
                ok2, out2, jid2 = wf._sbatch_guarded({}, s["dir"], "fakehost", jobname="J", fanout="d*",
                                                     only=["d2", "d3"])
                ok3, out3, _ = wf._sbatch_guarded({}, s["dir"], "fakehost", jobname="J", fanout="d*", only=["d1"])
            self.assertTrue(ok1, out1)
            self.assertTrue(ok2, out2)                       # d1 在跑不挡 d2/d3
            self.assertEqual(len(jid2.split(",")), 2)
            self.assertFalse(ok3)                            # d1 自己仍去重
            self.assertIn("已有作业", out3)
            self.assertEqual(D._read_state(state)["sbatch_count"], 3)


class CollectorFanRunningTests(unittest.TestCase):
    def test_fan_running_lists_subdirs_with_jobs(self):
        import types
        src = (ROOT / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
        C = types.ModuleType("_collector_defs")
        exec(compile(src[:src.index("\ndef main():")], "_collector_remote.py", "exec"), C.__dict__)
        with tempfile.TemporaryDirectory() as td:
            st = Path(td) / "CrS2" / "step7_deform"
            for n in ("deform-01", "deform-02", "undeformed"):
                (st / n).mkdir(parents=True)
                (st / n / "INCAR").write_text("ISTART = 0\n")
            jobs = {os.path.normpath(str(st / "deform-02")): {"id": "77", "state": "R"}}
            t = {"key": "ke-dft-cpu", "root": td, "materials": ["CrS2"],
                 "steps": [{"name": "step7_deform", "label": "S7_deform", "check": "deform",
                            "fanout": "*deform*"}]}
            res = C.collect_type(t, jobs)
        s = res["materials"][0]["steps"][0]
        self.assertEqual(s["fan_running"], ["deform-02"])
        self.assertEqual(s["fan_todo"], ["deform-01", "undeformed"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
