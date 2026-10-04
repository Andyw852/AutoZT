# -*- coding: utf-8 -*-
"""V161：init 种下、没人改过的模板副本不再盖住技能后来的修复。

起因：init 把技能整套 *.tpl/*.conf 复制进 project_setting/templates/，查找链里项目副本排在技能前面，
没人改过的副本就把当时的技能版本永远钉住 —— P1 的 LPEAD=.TRUE.（DFPT NaN）、Mo2S3 的松 EDIFFG、
P1 S4 的 conda activate amset_clean 都是这么来的；Si_diamond 刚 init 又种了 19 个。
规则：副本 = 技能同路径文件的某个 git 历史旧版 -> 当它不存在；手改过的照旧优先；与现行版相同的照常用。
离线：临时 git 仓库，不连超算。
"""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from autozt import report as rp  # noqa: E402

V1 = "LEPSILON = .TRUE.\nLPEAD    = .TRUE.\n"
V2 = "LEPSILON = .TRUE.\nLPEAD    = .FALSE.\n"
OLD_SUB = "#!/bin/bash\nconda activate amset_clean\n{{AMSET_CMD}}\n"
CONF1 = "[params]\nMESH_MIN = 1\n"
CONF2 = "[params]\nMESH_MIN = 2\n"


def _git(repo, *a):
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t"] + list(a),
                   check=True, capture_output=True)


class SeedTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        base = Path(self.td.name)
        self.repo = base / "repo"
        self.skill = self.repo / "skill" / "ke-dft-cpu"
        (self.skill / "step5_dielect").mkdir(parents=True)
        (self.skill / "step4_wave").mkdir(parents=True)
        _git(self.repo, "init", "-q")
        (self.skill / "step5_dielect" / "incar_dfpt_2d.tpl").write_text(V1)
        (self.skill / "step4_wave" / "submit_amset.tpl").write_text(OLD_SUB)
        (self.skill / "step.conf").write_text(CONF1)
        _git(self.repo, "add", ".")
        _git(self.repo, "commit", "-q", "-m", "v1")
        (self.skill / "step5_dielect" / "incar_dfpt_2d.tpl").write_text(V2)
        (self.skill / "step.conf").write_text(CONF2)
        _git(self.repo, "rm", "-q", "skill/ke-dft-cpu/step4_wave/submit_amset.tpl")   # 技能里删掉（挪去集群模板）
        _git(self.repo, "add", ".")
        _git(self.repo, "commit", "-q", "-m", "v2")
        self.ps = base / "Si_diamond" / "project_setting"
        (self.ps / "templates" / "step5_dielect").mkdir(parents=True)
        (self.ps / "templates" / "step4_wave").mkdir(parents=True)
        self.home = base / "home"
        self.t = {"key": "ke-dft-cpu", "skill_dir": str(self.skill), "template_dir": ".", "steps": []}
        rp._SEED_HIST.clear()
        rp._SEED_NOTED.clear()

    def _m(self, setting=None, hpc=None):
        return {"name": "Si_diamond", "ps": {"dir": str(self.ps), "setting": setting or {}},
                "template_map": {}, "hpc_name": hpc}

    def _find(self, fname, sname, m=None):
        with patch.dict(os.environ, {"HOME": str(self.home)}):
            return rp.find_asset({}, self.t, m or self._m(), fname, sname)

    def _copy(self, rel, text):
        p = self.ps / "templates" / rel
        p.write_text(text)
        return str(p)

    def test_old_seed_skipped(self):
        self._copy("step5_dielect/incar_dfpt_2d.tpl", V1)
        got = self._find("incar_dfpt_2d.tpl", "step5_dielect")
        self.assertEqual(Path(got).resolve(), (self.skill / "step5_dielect" / "incar_dfpt_2d.tpl").resolve())

    def test_edited_copy_wins(self):
        c = self._copy("step5_dielect/incar_dfpt_2d.tpl", V1 + "NBANDS = 64\n")
        self.assertEqual(self._find("incar_dfpt_2d.tpl", "step5_dielect"), c)

    def test_current_copy_used(self):
        c = self._copy("step5_dielect/incar_dfpt_2d.tpl", V2)
        self.assertEqual(self._find("incar_dfpt_2d.tpl", "step5_dielect"), c)

    def test_opt_out(self):
        c = self._copy("step5_dielect/incar_dfpt_2d.tpl", V1)
        m = self._m({"templates_follow_skill": False})
        self.assertEqual(self._find("incar_dfpt_2d.tpl", "step5_dielect", m), c)

    def test_deleted_skill_file_falls_to_cluster(self):
        self._copy("step4_wave/submit_amset.tpl", OLD_SUB)
        cl = self.home / ".config" / "autozt" / "setting" / "tcl" / "templates"
        cl.mkdir(parents=True)
        (cl / "submit_amset.tpl").write_text("conda activate {{AMSET_ENV}}\n")
        got = self._find("submit_amset.tpl", "step4_wave", self._m(hpc="tcl"))
        self.assertEqual(Path(got).resolve(), (cl / "submit_amset.tpl").resolve())

    def test_no_git_history_keeps_copy(self):
        c = self._copy("step5_dielect/incar_dfpt_2d.tpl", "LPEAD = .TRUE.\n# 比仓库历史更早的快照\n")
        self.assertEqual(self._find("incar_dfpt_2d.tpl", "step5_dielect"), c)

    def test_step_conf_old_seed_dropped(self):
        (self.ps / "templates" / "step.conf").write_text(CONF1)
        with patch.dict(os.environ, {"HOME": str(self.home)}):
            src = [p for p, _n in rp.step_conf_sources({}, self.t, self._m(), "step5_dielect")]
        self.assertNotIn(str(self.ps / "templates" / "step.conf"), src)
        self.assertIn(os.path.normpath(str(self.skill / "step.conf")), src)
        (self.ps / "templates" / "step.conf").write_text(CONF1 + "AMSET_ENV = amset051\n")   # 手改过
        with patch.dict(os.environ, {"HOME": str(self.home)}):
            src = [p for p, _n in rp.step_conf_sources({}, self.t, self._m(), "step5_dielect")]
        self.assertIn(str(self.ps / "templates" / "step.conf"), src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
