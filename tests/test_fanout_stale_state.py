"""Offline regression test: stale fanout state must not block a fresh gen.

历史故障（Si_ke34 / S7_deform，2026-09-28）—— 路径为**推断**，原始终端输出未留存:
  `autozt -p Si_ke34 start` 这一次调用里，采集发生在 gen **之前** —— 那时远端
  S7_deform 下还留着上一轮的 4 个子目录 + 一个旧的 undeformed，且都被判成"完成"，
  于是采集写出 fan_todo=[]。同一次 start 随后 gen 出新输入、准备提交，
  `remote_sbatch_fanout` 读到**过期的** fan_todo=[] 就拒绝了：
      "扇出步骤的待补清单为空（fan_todo=[]）：远端 N 个子目录都已判定完成，已拒绝提交。"
  结果：刚生成的新子目录一个都没交出去。

修复（autozt/workflow.py::do_submit，gen 成功前后）:
  gen 一旦真的运行，上一轮采到的扇出状态整体作废 —— 删掉 fan_todo / subs /
  fan_done / fan_jobids，提交落回"首次 gen 之后提交全部"这条路。

本测试不需要集群：run_remote / remote_gen / _sbatch_guarded 全部打桩。
"""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import autozt
from autozt import workflow as wf


STALE = {
    "fan_todo": [],
    "subs": ["deform-01", "deform-02", "deform-03", "deform-04", "undeformed"],
    "fan_done": ["deform-01", "deform-02", "deform-03", "deform-04", "undeformed"],
    "fan_jobids": "3910001,3910002",
}


def _step(**over):
    s = {"name": "S7_deform", "label": "S7_deform", "dir": "/remote/step7_deform",
         "fanout": "deform-*", "has_incar": False, "job": None, "_host": "mock",
         "submit": None, "tt": "ke-dft-cpu"}
    s.update(STALE)
    s.update(over)
    return s


class StaleFanoutStateTests(unittest.TestCase):
    # ---- 1. 复现：过期的 fan_todo=[] 会被闸门拒绝 -------------------------
    def test_stale_empty_todo_is_refused(self):
        """旧状态的 fan_todo=[] + 远端确有子目录 -> 拒绝提交（这就是那次事故）。"""
        s = _step()
        with patch.object(autozt, "run_remote", return_value=(0, "5\n")) as rr, \
             patch.object(wf, "_sbatch_guarded") as sbatch:
            ok, out, jid = wf.remote_sbatch({}, s, jobname="J")
        self.assertFalse(ok)
        self.assertIsNone(jid)
        self.assertIn("待补清单为空", out)
        self.assertIn("5", out)                     # 远端子目录数写进提示
        sbatch.assert_not_called()                  # 一个 sbatch 都没发生
        self.assertEqual(rr.call_count, 1)          # 只数了一次子目录

    # ---- 2. 修复：同一次 do_submit 里 gen 之后，旧名单作废 -> 能提交 ----
    def test_gen_then_submit_not_blocked_by_stale_state(self):
        """gen 成功后 fan_todo/subs/fan_done/fan_jobids 必须消失，提交全部。"""
        s = _step()
        calls = []

        def fake_guarded(cfg, sdir, host, jobname=None, fanout=None, only=None,
                         force=False, submit=None, empty_fanout_manifest=None):
            calls.append({"dir": sdir, "fanout": fanout, "only": only})
            return (True, "", "3914906,3914907")

        m = {"name": "Si_ke34/ke-dft-cpu", "host_eff": "mock", "tt": "ke-dft-cpu"}
        with patch.object(autozt, "step_cfg", return_value={"run": "sbatch"}), \
             patch.object(autozt, "log_action"), \
             patch.object(autozt, "run_remote", return_value=(0, "3\n")), \
             patch.object(wf, "remote_gen", return_value=(True, "gen ok")), \
             patch.object(wf, "_fetch_stamp_clear"), \
             patch.object(wf, "_relay_prev_across_host"), \
             patch.object(wf, "_scancel_clear"), \
             patch.object(wf, "_remote_submit_preflight", return_value=(True, "")), \
             patch.object(wf, "_sbatch_guarded", side_effect=fake_guarded):
            # 反证（修复前的真实结局）：同样这份过期状态若原样交给提交闸门 -> 被拒。
            ok_before, out_before, _ = wf.remote_sbatch({}, _step(), jobname="J")
            self.assertFalse(ok_before)
            self.assertIn("待补清单为空", out_before)
            calls.clear()                            # 上面那次拒绝不该产生 sbatch
            ok = wf.do_submit({}, "ke-dft-cpu", m, s, False, False, False, "Step")
        self.assertTrue(ok)
        self.assertEqual(len(calls), 1)
        self.assertIsNone(calls[0]["only"])          # 修好后 = 提交全部（不是旧名单）
        for k in ("fan_todo", "subs", "fan_done", "fan_jobids"):
            self.assertNotIn(k, s, "%s 在 gen 之后必须被清掉" % k)

    # ---- 3. 半新半旧：旧名单里的 deform-04 已不存在，新子目录不能被漏掉 ----
    def test_partial_stale_list_does_not_skip_new_subdirs(self):
        s = _step(fan_todo=["deform-04"], subs=["deform-04"],
                  fan_done=["deform-04"], fan_jobids="3910001")
        captured = {}

        def fake_guarded(cfg, sdir, host, jobname=None, fanout=None, only=None,
                         force=False, submit=None, empty_fanout_manifest=None):
            captured["only"] = only
            return (True, "", "3914906,3914907,3914908,3914909")

        m = {"name": "Si_ke34/ke-dft-cpu", "host_eff": "mock", "tt": "ke-dft-cpu"}
        with patch.object(autozt, "step_cfg", return_value={"run": "sbatch"}), \
             patch.object(autozt, "log_action"), \
             patch.object(autozt, "run_remote", return_value=(0, "4\n")), \
             patch.object(wf, "remote_gen", return_value=(True, "gen ok")), \
             patch.object(wf, "_fetch_stamp_clear"), \
             patch.object(wf, "_relay_prev_across_host"), \
             patch.object(wf, "_scancel_clear"), \
             patch.object(wf, "_remote_submit_preflight", return_value=(True, "")), \
             patch.object(wf, "_sbatch_guarded", side_effect=fake_guarded):
            ok = wf.do_submit({}, "ke-dft-cpu", m, s, False, False, False, "Step")
        self.assertTrue(ok)
        self.assertIsNone(captured["only"])          # 不能只剩 [deform-04]

    # ---- 4. 反向保护：没有 fanout 的普通步骤不受影响 ----------------------
    def test_non_fanout_step_untouched(self):
        s = _step(fanout=None, has_incar=True, fan_todo=[])
        captured = {}

        def fake_guarded(cfg, sdir, host, jobname=None, fanout=None, only=None,
                         force=False, submit=None, empty_fanout_manifest=None):
            captured["fanout"] = fanout
            return (True, "", "1")

        m = {"name": "Si_ke34/ke-dft-cpu", "host_eff": "mock", "tt": "ke-dft-cpu"}
        with patch.object(autozt, "step_cfg", return_value={"run": "sbatch"}), \
             patch.object(autozt, "log_action"), \
             patch.object(autozt, "run_remote", return_value=(0, "0\n")), \
             patch.object(wf, "remote_gen", return_value=(True, "gen ok")), \
             patch.object(wf, "_fetch_stamp_clear"), \
             patch.object(wf, "_relay_prev_across_host"), \
             patch.object(wf, "_scancel_clear"), \
             patch.object(wf, "_remote_submit_preflight", return_value=(True, "")), \
             patch.object(wf, "_sbatch_guarded", side_effect=fake_guarded):
            ok = wf.do_submit({}, "ke-dft-cpu", m, s, False, False, False, "Step")
        self.assertTrue(ok)
        self.assertIsNone(captured["fanout"])
        self.assertIn("fan_todo", s)                 # 非扇出步骤的键不该被乱删

    # ---- 5. 提交被拒也要落 tf.log（2026-09-28 用户要求） --------------------
    def test_refused_submit_is_logged(self):
        s = _step(has_incar=True)          # 不 gen，直接走提交闸门
        logged = []
        m = {"name": "Si_ke34/ke-dft-cpu", "host_eff": "mock", "tt": "ke-dft-cpu"}
        with patch.object(autozt, "step_cfg", return_value={"run": "sbatch"}), \
             patch.object(autozt, "log_action",
                          side_effect=lambda _m, line: logged.append(line)), \
             patch.object(autozt, "run_remote", return_value=(0, "5\n")), \
             patch.object(wf, "_remote_submit_preflight", return_value=(True, "")), \
             patch.object(wf, "_sbatch_guarded") as sbatch:
            ok = wf.do_submit({}, "ke-dft-cpu", m, s, False, False, False, "Step")
        self.assertFalse(ok)
        sbatch.assert_not_called()
        self.assertTrue(any("提交失败" in x and "待补清单为空" in x for x in logged), logged)


if __name__ == "__main__":
    unittest.main(verbosity=2)
