# -*- coding: utf-8 -*-
"""V160：有可开始步骤、但 auto-advance 不会推进的材料，summary 要点名并说明原因。

起因：Si_diamond 的 project_setting/ 里没有 setting.yaml（没写 auto_advance: true），6 个就绪步骤一直 PREP、
不提交；表里只显示"可开始"，summary 只算进 wait，大家以为它在"等自动重排"。auto-advance 的跳过提示也一律说
"auto_advance: false"，其实是没写。显式 true 才推进是有意的（防误提交），不改；改的是让它看得见。
离线：不连超算。
"""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from autozt import report as rp  # noqa: E402
from autozt import workflow as wf  # noqa: E402


def _mat(name, setting, kinds, ps_dir="/p/ps"):
    steps = [{"name": "s%d" % i, "label": "S%d" % i, "kind": k} for i, k in enumerate(kinds)]
    return {"name": name, "tt": "ke-dft-cpu", "steps": steps,
            "ps": {"dir": ps_dir, "setting": setting}}


class ReasonTests(unittest.TestCase):
    def test_reasons(self):
        self.assertIsNone(rp.auto_off_reason(_mat("a", {"auto_advance": True}, [])))
        self.assertIn("没写", rp.auto_off_reason(_mat("a", {}, [])))
        self.assertIn("没写", rp.auto_off_reason(_mat("a", None, [])))
        self.assertEqual(rp.auto_off_reason(_mat("a", {"auto_advance": False}, [])), "auto_advance: false")
        self.assertIn("未 init", rp.auto_off_reason(_mat("a", {}, [], ps_dir=None)))


class SummaryTests(unittest.TestCase):
    def test_manual_line(self):
        data = {"types": [{"key": "ke-dft-cpu", "materials": [
            _mat("Si_diamond", {}, ["OK", "PREP", "PREP", "WAIT"]),
            _mat("P1_Mo-MoS2", {"auto_advance": True}, ["OK", "PREP"]),
            _mat("Mo2S3", {"auto_advance": False}, ["OK", "WAIT"]),          # 没有可开始的步骤：不报
            _mat("GaAs", {}, ["OK", "OK"]),
        ]}]}
        lines = rp._summary_lines(data)
        self.assertTrue(lines[0].startswith("ke-dft-cpu: 4 材料 done=1"))
        man = [l for l in lines if l.lstrip().startswith("手动")]
        self.assertEqual(len(man), 1)
        self.assertIn("Si_diamond", man[0])
        self.assertIn("没写 auto_advance", man[0])
        self.assertIn("S1,S2", man[0])
        self.assertIn("-p Si_diamond auto on", man[0])


class AutoAdvanceSkipTests(unittest.TestCase):
    def test_skip_message_names_reason(self):
        data = {"types": [{"key": "ke-dft-cpu", "materials": [
            _mat("Si_diamond", {}, ["PREP"]),
            _mat("Mo2S3", {"auto_advance": False}, ["PREP"]),
        ]}]}
        out = []
        with patch.object(wf, "_SkillGate", lambda cfg, data: None), \
                patch("builtins.print", lambda *a, **k: out.append(" ".join(str(x) for x in a))):
            wf.auto_advance({"auto_advance": True}, data)
        msg = "\n".join(out)
        self.assertIn("Si_diamond[ke-dft-cpu]（setting.yaml 没写 auto_advance）", msg)
        self.assertIn("Mo2S3[ke-dft-cpu]（auto_advance: false）", msg)
        self.assertNotIn("项目级 auto_advance: false）：", msg)


if __name__ == "__main__":
    unittest.main(verbosity=2)
