# -*- coding: utf-8 -*-
"""单步换超算对可选组步骤生效（Mg2C60 S2.3_hse 切 jzzn 不生效的回归）。

  · 项目 steps 里给可选组步骤写的 hpc：expand_optional_steps 重新插入组定义时要带上（原来静默丢失）；
  · 项目 step_hpc: {label|name|seq: 集群}：不用整段抄 steps，可选组步骤、被关掉的组定义都生效；
    写错的键告警；步骤 dict 在同技能各段之间共享，不能原地改。
"""
import copy
import io
import contextlib
import tempfile
import unittest
from pathlib import Path

import yaml

import autozt

REPO = Path(__file__).resolve().parents[1]
SK = yaml.safe_load((REPO / "skill" / "ke-dft-cpu" / "skill.yaml").read_text(encoding="utf-8"))
S23 = "step2_bandgap/step2.3_hse"


def _t(steps=None, **kw):
    t = {"steps": copy.deepcopy(steps if steps is not None else SK["steps"]),
         "optional_steps": copy.deepcopy(SK["optional_steps"])}
    t.update(kw)
    return t


def _step(t, label):
    return next(s for s in t["steps"] if s.get("label") == label)


class ExpandCarriesHpcTests(unittest.TestCase):
    def test_project_steps_hpc_on_optional_step_survives(self):
        grp = next(d for d in SK["optional_steps"]["bandgap_hse"]["steps"] if d["name"] == S23)
        t = _t(SK["steps"] + [dict(grp, hpc="jzzn")])
        autozt.expand_optional_steps(t)
        self.assertEqual(_step(t, "S2.3_hse").get("hpc"), "jzzn")

    def test_group_default_hpc_kept_and_project_override_wins(self):
        grp = {"g": {"default": True, "steps": [{"name": "stepX", "label": "SX", "after": "step1_opt",
                                                   "hpc": "clusterA"}]}}
        base = [{"name": "step1_opt", "label": "S1_opt"}]
        t = {"steps": list(base), "optional_steps": copy.deepcopy(grp)}
        autozt.expand_optional_steps(t)
        self.assertEqual(_step(t, "SX")["hpc"], "clusterA")          # 技能组定义自带的保留
        t = {"steps": base + [{"name": "stepX", "hpc": "clusterB"}], "optional_steps": copy.deepcopy(grp)}
        autozt.expand_optional_steps(t)
        self.assertEqual(_step(t, "SX")["hpc"], "clusterB")          # 项目写的优先


class StepHpcTests(unittest.TestCase):
    def _apply(self, t):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            autozt.expand_optional_steps(t)
            autozt.apply_step_hpc(t)
        return err.getvalue()

    def test_label_name_and_typo(self):
        shared = SK["steps"]                                  # 同技能各段共享的步骤 dict
        before = copy.deepcopy(shared)
        t = _t(steps=shared, step_hpc={"S2.3_hse": "jzzn", "step3_uniform": "jzzn",
                                       "step2.3_hse": "jzzn"}, _from="tf_Mg2C60_ke-dft-cpu.yaml")
        t["steps"] = shared                                    # 不拷贝：验证不原地改
        msg = self._apply(t)
        self.assertEqual(_step(t, "S2.3_hse")["hpc"], "jzzn")
        self.assertEqual(_step(t, "S3_uniform")["hpc"], "jzzn")
        self.assertNotIn("hpc", _step(t, "S6_elastic"))
        self.assertIn("'step2.3_hse' 对不上任何步骤", msg)
        self.assertEqual(shared, before)                      # 共享的步骤 dict 没被改

    def test_disabled_group_definition_gets_hpc(self):
        t = _t(bandgap_hse=False, step_hpc={"S2.3_hse": "jzzn"})
        self._apply(t)
        self.assertFalse(any(s.get("label") == "S2.3_hse" for s in t["steps"]))
        d = next(d for ds in t["_optional_off"].values() for d in ds if d.get("label") == "S2.3_hse")
        self.assertEqual(d["hpc"], "jzzn")                    # -j S2.3_hse start 按需启用时也走 jzzn

    def test_get_types_segment_only(self):
        """两个项目挂同一技能：只有写了 step_hpc 的那个段生效。"""
        a, b = tempfile.mkdtemp(), tempfile.mkdtemp()
        base = dict(copy.deepcopy(SK), local_root=a)
        base["_segments"] = [{"local_root": b, "step_hpc": {"S2.3_hse": "jzzn"}, "_from": "tf_b.yaml"}]
        with contextlib.redirect_stderr(io.StringIO()):
            types = autozt.get_types({"task_types": {"ke-dft-cpu": base}}, quiet=True)
        by_root = {t["local_root"]: t for t in types}
        self.assertEqual(_step(by_root[b], "S2.3_hse").get("hpc"), "jzzn")
        self.assertIsNone(_step(by_root[a], "S2.3_hse").get("hpc"))


if __name__ == "__main__":
    unittest.main()
