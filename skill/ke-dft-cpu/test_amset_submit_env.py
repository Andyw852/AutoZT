#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_amset_submit_env.py —— V159：AMSET 提交脚本激活的环境必须是 step.conf 的 AMSET_ENV（不需要 VASP/AMSET）。

用法：python test_amset_submit_env.py      退出码 0 = 全部 PASS。

起因：P1_Mo-MoS2 的 S4_wave 作业报 EnvironmentNameNotFound: amset_clean -> amset: command not found。
8 月 init 时拷进 project_setting/templates/ 的 submit_amset.tpl 是 09-22 之前的版本，写死
conda activate amset_clean（AMSET 0.4.19，已删），没有 {{AMSET_ENV}} 占位符；项目级模板优先，
gen 只查"占位符有没有残留"，照常生成、提交。环境还在的集群上会悄悄用 0.4.19 算（形变势减半，V75）。
"""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
COMMON = ROOT.parent / "_common" / "opt"
for _p in (str(ROOT), str(COMMON), str(ROOT / "tools")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402
import template_drift as TD  # noqa: E402

OLD_TPL = """#!/bin/bash
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
source /home/x/miniconda3/etc/profile.d/conda.sh
conda activate amset_clean
{{AMSET_CMD}}
"""
NEW_TPL = OLD_TPL.replace("conda activate amset_clean", "conda activate {{AMSET_ENV}}")
LOCAL_TPL = """#!/bin/bash
AMSET_ENV="{{AMSET_ENV}}"
if command -v conda >/dev/null 2>&1; then
    conda activate "${AMSET_ENV}" || echo "[WARN] conda env ${AMSET_ENV} not found" >&2
fi
{{AMSET_CMD}}
"""


def _render(raw, env="amset051"):
    return raw.replace("{{JOBNAME}}", "j").replace("{{AMSET_CMD}}", "amset wave").replace("{{AMSET_ENV}}", env)


class CheckSubmitTests(unittest.TestCase):
    def test_old_template_refused(self):
        with self.assertRaises(SystemExit) as cm:
            kc.check_amset_submit(OLD_TPL, _render(OLD_TPL), "amset051")
        msg = str(cm.exception)
        self.assertIn("amset_clean", msg)
        self.assertIn("0.4.19", msg)
        self.assertIn("stale", msg)

    def test_current_templates_pass(self):
        self.assertEqual(kc.check_amset_submit(NEW_TPL, _render(NEW_TPL), "amset051"), [])
        self.assertEqual(kc.check_amset_submit(LOCAL_TPL, _render(LOCAL_TPL), "amset051"), [])

    def test_hardcoded_other_env_refused(self):
        raw = OLD_TPL.replace("amset_clean", "amset")
        with self.assertRaises(SystemExit):
            kc.check_amset_submit(raw, _render(raw), "amset051")
        raw = OLD_TPL.replace("conda activate amset_clean", "conda run -n amset_clean amset wave")
        with self.assertRaises(SystemExit):
            kc.check_amset_submit(raw, _render(raw), "amset051")

    def test_hardcoded_same_env_warns(self):
        raw = OLD_TPL.replace("amset_clean", "amset051")
        w = kc.check_amset_submit(raw, _render(raw), "amset051")
        self.assertEqual(len(w), 1)
        self.assertIn("{{AMSET_ENV}}", w[0])

    def test_comments_ignored(self):
        raw = NEW_TPL + "# 旧：conda activate amset_clean\n"
        self.assertEqual(kc.check_amset_submit(raw, _render(raw), "amset051"), [])
        self.assertEqual(kc.amset_submit_envs("source activate old_env && x\n  conda activate  'a1' ; y"),
                         ["old_env", "a1"])


class EnvNameTests(unittest.TestCase):
    def _conf(self, val):
        d = Path(tempfile.mkdtemp())
        (d / "step.conf").write_text("[params]\nAMSET_ENV = %s\n" % val, encoding="utf-8")
        return d

    def test_old_env_in_step_conf_refused(self):
        with self.assertRaises(SystemExit) as cm:
            kc.amset_env_name(self._conf("amset_clean"))
        self.assertIn("0.4.19", str(cm.exception))
        self.assertEqual(kc.amset_env_name(self._conf("amset051")), "amset051")


GENS = ("step4_wave/gen_step6_wave.py", "step4b_wave_full/gen_step6b_wave_full.py",
        "step8_amset/gen_step10_amset.py", "step8.4_amset2d/gen_step14_amset2d.py")


class WiringTests(unittest.TestCase):
    def test_all_amset_gens_check_after_render(self):
        for g in GENS:
            src = (ROOT / g).read_text(encoding="utf-8")
            i = src.index('.replace("{{AMSET_ENV}}"')
            k = src.index("kc.check_amset_submit(_tpl_raw, text", i)
            self.assertLess(k, src.index("submit.write_text(text", i), g)
            self.assertLess(src.index("_tpl_raw = text"), i, g)


class S4GenTests(unittest.TestCase):
    """真跑 S4 的 gen：旧模板 -> 拒绝、不写 submit.sh；新模板 -> conda activate amset051。"""

    def _mat(self, tpl):
        d = Path(tempfile.mkdtemp())
        (d / "step.conf").write_text("[params]\nAMSET_ENV = amset051\n", encoding="utf-8")
        (d / "step3_uniform").mkdir()
        for f in ("vasprun.xml", "WAVECAR"):
            (d / "step3_uniform" / f).write_text("x")
        (d / "submit_amset.tpl").write_text(tpl, encoding="utf-8")
        return d

    def _run(self, d):
        env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT), str(COMMON)]))
        return subprocess.run([sys.executable, str(ROOT / GENS[0])], cwd=str(d), env=env,
                              capture_output=True, text=True, timeout=120)

    def test_old_template(self):
        d = self._mat(OLD_TPL)
        r = self._run(d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("amset_clean", r.stderr)
        self.assertFalse((d / "step4_wave" / "submit.sh").exists())

    def test_new_template(self):
        d = self._mat(NEW_TPL)
        r = self._run(d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        sub = (d / "step4_wave" / "submit.sh").read_text(encoding="utf-8")
        self.assertIn("conda activate amset051", sub)
        self.assertNotIn("amset_clean", sub)


class DriftEnvTests(unittest.TestCase):
    def test_flags_old_submit_template(self):
        root = Path(tempfile.mkdtemp())
        for name, text in (("P1_Mo-MoS2", OLD_TPL), ("GaAs", NEW_TPL)):
            t = root / name / "project_setting" / "templates"
            t.mkdir(parents=True)
            (t / "submit_amset.tpl").write_text(text, encoding="utf-8")
        rows = TD.audit(root, skill_dir=ROOT, repo=None)
        env = [r for r in rows if r["kind"] == "ENV"]
        self.assertEqual(len(env), 1)
        self.assertIn("P1_Mo-MoS2", env[0]["file"])
        self.assertIn("amset_clean", env[0]["note"])
        self.assertEqual(TD.main([str(root)]), 1)
        self.assertIsNone(TD.submit_env_issue(LOCAL_TPL))
        self.assertIn("占位符", TD.submit_env_issue(OLD_TPL.replace("amset_clean", "amset051")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
