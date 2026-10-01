#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_gap_warn_drift.py —— V137：小带隙告警的 TypeError、项目级模板副本漂移审计、S5 的 LPEAD 告警。

用法：python test_gap_warn_drift.py      退出码 0 = 全部 PASS。不需要 VASP / AMSET。

1) S8 / S8.4 的 _warn_small_pbe_gap：PBE < 0.5×HSE（P1_Al-AlN 0.11 vs 0.94、GaAs 0.4175 vs ~1.4）时
   原来 print 的格式串残留 %.4f/%.2f、只传一个字符串 -> TypeError，gen 直接崩。现在正常告警；
   没有 HSE 时按 0.3 eV 绝对阈值；带隙够大不告警。
2) tools/template_drift.py：项目 templates/ 里的副本
   · 等于技能模板的某个 git 历史版本 -> [STALE]，列出 ★ 关键键（LPEAD）；
   · 不是任何历史版本 -> [CUSTOM]；与现版完全相同 -> 不报；submit*.tpl 默认跳过；没有 git -> [DIFF]。
3) S5 的 gen 在 LPEAD=.TRUE. + IBRION=7/8 时告警（位置在 step.conf 的 [incar] 覆盖之后）。
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "tools"),
           str(ROOT / "step8_amset"), str(ROOT / "step8.4_amset2d")):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _quiet_import(name):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return __import__(name)


def _bands(pbe, hse=None):
    d = Path(tempfile.mkdtemp())
    for sub, g in (("step2.2_pbe_plot", pbe), ("step2.3_hse_plot", hse)):
        if g is None:
            continue
        p = d / "step2_bandgap" / sub / "band_summary.json"
        p.parent.mkdir(parents=True)
        p.write_text(json.dumps({"gap_eV": g}))
    return d


class SmallGapWarnTests(unittest.TestCase):
    def test_both_gens(self):
        for mod in ("gen_step10_amset", "gen_step14_amset2d"):
            G = _quiet_import(mod)
            for pbe, hse, want in ((0.11, 0.94, "0.1100 eV = HSE 带隙 0.9400 eV 的 12%"),
                                   (0.4175, 1.42, "0.4175 eV = HSE 带隙 1.4200 eV 的 29%"),
                                   (0.2, None, "0.2000 eV（< 0.3 eV 绝对阈值"),
                                   (1.0, 1.5, None)):
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    g = G._warn_small_pbe_gap(_bands(pbe, hse))       # 原来在这里 TypeError
                self.assertAlmostEqual(g, pbe)
                out = buf.getvalue()
                if want is None:
                    self.assertNotIn("[WARN]", out, mod)
                else:
                    self.assertIn("带隙只有 " + want, out, mod)


V1 = "SYSTEM = {{SYSTEM}}\nIBRION   = 8\nLEPSILON = .TRUE.\nLPEAD    = .TRUE.\nNCORE = 1\n"
V2 = "SYSTEM = {{SYSTEM}}\nIBRION   = 8\nLEPSILON = .TRUE.\n# 默认改成 .FALSE.\nLPEAD    = .FALSE.\nNCORE = 1\n"


def _git(repo, *a):
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t"] + list(a),
                   check=True, capture_output=True)


class TemplateDriftTests(unittest.TestCase):
    def setUp(self):
        import template_drift as TD
        self.TD = TD
        self.repo = Path(tempfile.mkdtemp())
        _git(self.repo, "init", "-q")
        self.skill = self.repo / "skill" / "ke-dft-cpu"
        tpl = self.skill / "step5_dielect" / "incar_dfpt_2d.tpl"
        tpl.parent.mkdir(parents=True)
        tpl.write_text(V1)
        _git(self.repo, "add", ".")
        _git(self.repo, "commit", "-q", "-m", "v1")
        tpl.write_text(V2)
        (self.skill / "step5_dielect" / "submit_std_2d.tpl").write_text("#!/bin/bash\n")
        _git(self.repo, "add", ".")
        _git(self.repo, "commit", "-q", "-m", "v2")

    def _project(self, files):
        root = Path(tempfile.mkdtemp())
        for name, text in files.items():
            p = root / name / "project_setting" / "templates" / "step5_dielect"
            p.mkdir(parents=True)
            (p / "incar_dfpt_2d.tpl").write_text(text)
            (p / "submit_std_2d.tpl").write_text("#!/bin/bash\n#SBATCH -p site\n")
        return root

    def test_classify(self):
        root = self._project({"P1_Mo-MoS2": V1, "custom": V2.replace("NCORE = 1", "NCORE = 4"), "same": V2})
        rows = {Path(r["file"]).parts[-5]: r for r in self.TD.audit(root, self.skill, self.repo)}
        self.assertEqual(set(rows), {"P1_Mo-MoS2", "custom"})            # same 不报；submit 默认跳过
        self.assertEqual(rows["P1_Mo-MoS2"]["kind"], "STALE")
        self.assertTrue(rows["P1_Mo-MoS2"]["critical"])
        self.assertIn(("LPEAD", ".TRUE.", ".FALSE."), rows["P1_Mo-MoS2"]["diff"])
        self.assertEqual(rows["custom"]["kind"], "CUSTOM")
        self.assertFalse(rows["custom"]["critical"])                     # 只改了 NCORE
        self.assertEqual(rows["custom"]["diff"], [("NCORE", "4", "1")])
        rows_all = self.TD.audit(root, self.skill, self.repo, include_submit=True)
        self.assertTrue(any(r["file"].endswith("submit_std_2d.tpl") for r in rows_all))

    def test_no_git(self):
        skill = Path(tempfile.mkdtemp())
        (skill / "step5_dielect").mkdir()
        (skill / "step5_dielect" / "incar_dfpt_2d.tpl").write_text(V2)
        root = self._project({"P1": V1})
        rows = self.TD.audit(root, skill, repo=skill)
        self.assertEqual([r["kind"] for r in rows], ["DIFF"])

    def test_cli_exit_code(self):
        root = self._project({"P1_Mo-MoS2": V1})
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            rc = self.TD.main([str(root), "--skill-dir", str(self.skill), "--repo", str(self.repo)])
        self.assertEqual(rc, 1)
        self.assertIn("[STALE]", buf.getvalue())
        self.assertIn("★LPEAD", buf.getvalue())


class S5LpeadWarnTests(unittest.TestCase):
    def test_position(self):
        src = (ROOT / "step5_dielect" / "gen_step8_dielect.py").read_text(encoding="utf-8")
        i_ov = src.index("stepconf.apply_incar_file(out / \"INCAR\", log=_ic_log)")
        i_w = src.index("[patch_lpead_warn V137]")
        i_sub = src.index("submit_tpl = resolve_tpl(")
        self.assertLess(i_ov, i_w)
        self.assertLess(i_w, i_sub)
        self.assertIn('_inc.get("IBRION", "").strip() in ("7", "8")', src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
