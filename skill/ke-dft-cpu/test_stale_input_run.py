#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_stale_input_run.py —— V168：已有产物不是现在这份 INCAR 跑的 -> 归档、判未完成（纯本地、无数据）。

起因（CrS₂ S7_deform，2026-10-06）：retry 换上带 LVHAR 的新 INCAR，undeformed / deform-05..08 的旧 OUTCAR
照样被 ck_deform 判完成、不进 fan_todo，start 只交了 5/10 个；S7.1 读不到这 5 个的 LOCPOT，
deformation_vac.h5 出不来。patch_stale_grid 只比 KPOINTS，snapshot_inputs 只在 gen 覆盖输入前拍快照——
INCAR 已被上一次 gen 换掉时，只有 vasprun.xml 里 VASP 回显的 <incar>/<parameters> 能说明问题。

用法：python test_stale_input_run.py      退出码 0 = 全部 PASS。
"""
import ast
import re
import sys
import tempfile
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
import ke_common as kc  # noqa: E402

# 格式照 VASP 6 的真实 vasprun.xml：<incar> 只列 INCAR 里写了的键；<parameters> 含默认值（PREC 截成 "accura"）
VASPRUN = """<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <generator>
  <i name="program" type="string">vasp </i>
 </generator>
 <incar>
  <i type="string" name="SYSTEM">CrS2 deform deform-05</i>
  <i type="string" name="PREC">accurate</i>
  <i name="ENCUT">    520.00000000</i>
  <i name="EDIFF">      0.00000001</i>
  <i type="int" name="NELM">   200</i>
  <i type="int" name="ISYM">     0</i>
  <i type="logical" name="LREAL"> F  </i>
  <i type="string" name="GGA">PE</i>
  <i type="string" name="ALGO">Normal</i>
  <i type="int" name="NCORE">     4</i>
  <v name="MAGMOM">      3.00000000      3.00000000      0.00000000</v>
 </incar>
 <kpoints>
  <varray name="kpointlist" >
   <v>       0.00000000       0.00000000       0.00000000 </v>
  </varray>
 </kpoints>
 <parameters>
  <separator name="general" >
   <i type="string" name="PREC">accura</i>
   <i type="int" name="ISTART">     0</i>
  </separator>
  <separator name="writing" >
   <i type="logical" name="LVHAR"> F  </i>
   <i type="logical" name="ADDGRID"> F  </i>
  </separator>
 </parameters>
 <calculation/>
</modeling>
"""
INCAR_OLD = """SYSTEM = CrS2 deform deform-05
PREC   = Accurate
ENCUT  = 520
EDIFF  = 1E-8        # 带边
NELM   = 200
ISYM   = 0
LREAL  = .FALSE.
GGA    = PE
ALGO   = N
NCORE  = 4
MAGMOM = 2*3.0 0
"""
INCAR_NEW = INCAR_OLD.replace("NCORE  = 4", "NCORE  = 6\nKPAR   = 2") + "LVHAR  = .TRUE.\n"


def _deform_dir(td, incar):
    d = Path(td) / "step7_deform" / "deform-05"
    d.mkdir(parents=True)
    (d / "INCAR").write_text(incar)
    (d / "POSCAR").write_text("x\n")
    (d / "KPOINTS").write_text("auto\n0\nGamma\n 46 46 3\n")
    (d / "vasprun.xml").write_text(VASPRUN)
    (d / "OUTCAR").write_text("... General timing and accounting informations for this job ...\n")
    return d


class RunIncarMismatchTests(unittest.TestCase):
    def test_same_physics_is_not_flagged(self):
        # NCORE/KPAR（换宿主机会变）、注释、n*x、ALGO 缩写、1E-8 vs 0.00000001 都不算变化
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD.replace("NCORE  = 4", "NCORE  = 6\nKPAR   = 2"),
                                               _write(VASPRUN)), [])

    def test_lvhar_added_is_flagged(self):
        self.assertEqual(kc.run_incar_mismatch(INCAR_NEW, _write(VASPRUN)), ["LVHAR: F -> .TRUE."])

    def test_value_change_and_removal(self):
        diff = kc.run_incar_mismatch(INCAR_OLD.replace("EDIFF  = 1E-8", "EDIFF  = 1E-6")
                                     .replace("ISYM   = 0\n", ""), _write(VASPRUN))
        self.assertEqual(diff, ["EDIFF: 0.00000001 -> 1E-6", "ISYM: 0 -> (删除)"])
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD.replace("2*3.0 0", "2*3.0 1"), _write(VASPRUN)),
                         ["MAGMOM: 3.00000000      3.00000000      0.00000000 -> 2*3.0 1"])
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD + "ADDGRID = .TRUE.\n", _write(VASPRUN)),
                         ["ADDGRID: F -> .TRUE."])          # 旧 INCAR 没写，按 <parameters> 的生效值比

    def test_output_only_tags(self):
        v = VASPRUN.replace(' </incar>', '  <i type="logical" name="LWAVE"> T  </i>\n </incar>')
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD + "LWAVE = .FALSE.\n", _write(v)), [])  # 关掉输出：不重算
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD, _write(v)), [])                         # 删掉输出键：同上
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD + "LVHAR = .FALSE.\n", _write(VASPRUN)), [])
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD + "LORBIT = 11\n", _write(VASPRUN)), [])  # 两处都没回显

    def test_tiny_numbers_and_unknown_tags(self):
        v = VASPRUN.replace('<i name="EDIFF">      0.00000001</i>', '<i name="EDIFF">      0.00000000</i>')
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD.replace("1E-8", "1E-9"), _write(v)), [])  # 8 位小数舍入
        self.assertEqual(kc.run_incar_mismatch(INCAR_OLD + "FOOBAR = 3\n", _write(VASPRUN)), [])   # VASP 不认识的键

    def test_undeterminable_returns_none(self):
        self.assertIsNone(kc.run_incar_mismatch(INCAR_OLD, "/nonexistent/vasprun.xml"))
        self.assertIsNone(kc.run_incar_mismatch(INCAR_OLD, _write("<modeling>\n <generator/>\n")))  # 跑到一半


def _write(text, _keep=[]):
    td = tempfile.TemporaryDirectory()
    _keep.append(td)
    p = Path(td.name) / "vasprun.xml"
    p.write_text(text)
    return p


class ArchiveTests(unittest.TestCase):
    def test_archive_makes_deform_not_done(self):
        import types
        src = (_ROOT.parents[1] / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
        C = types.ModuleType("_collector_defs")
        exec(compile(src[:src.index("\ndef main():")], "_collector_remote.py", "exec"), C.__dict__)
        with tempfile.TemporaryDirectory() as td:
            d = _deform_dir(td, INCAR_NEW)
            self.assertTrue(C.ck_deform(str(d), {})[0])                 # 修复前：旧 OUTCAR 判完成
            self.assertEqual(kc.archive_if_run_mismatch(d), ["LVHAR: F -> .TRUE."])
            self.assertFalse(C.ck_deform(str(d), {})[0])                # 归档后：进 fan_todo
            names = sorted(p.name for p in d.iterdir())
            self.assertIn("INCAR", names)                               # 输入不动
            self.assertTrue(any(n.startswith("OUTCAR.stale-input-") for n in names))
            self.assertTrue(any(n.startswith("vasprun.xml.stale-input-") for n in names))
            self.assertIsNone(kc.archive_if_run_mismatch(d))            # 再跑一次：没有回显可比，不重复归档

    def test_matching_run_is_kept(self):
        with tempfile.TemporaryDirectory() as td:
            d = _deform_dir(td, INCAR_OLD)
            self.assertEqual(kc.archive_if_run_mismatch(d), [])
            self.assertTrue((d / "OUTCAR").is_file())


class GenWiringTests(unittest.TestCase):
    SRC = (_ROOT / "step7_deform" / "gen_step9_deform.py").read_text(encoding="utf-8")

    def test_main_dirs_checked_after_incar_final(self):
        s = self.SRC
        i = s.index("kc.archive_if_run_mismatch(d)")
        self.assertLess(s.index("stepconf.apply_incar_file(d / \"INCAR\""), i)   # INCAR 定稿之后才比
        self.assertLess(i, s.index('sub = d / "submit.sh"'))
        self.assertLess(s.index('_snap_in = kc.snapshot_inputs(d, ("INCAR",))'), s.index("kc.render_tpl(tpl"))

    def test_ionrelax_reuse_requires_same_static_incar(self):
        s = self.SRC
        i = s.index("_why = kc.run_incar_mismatch(")
        self.assertIn('_ionrelax_incars((out / _name / "INCAR")', s[i:i + 200])
        self.assertLess(i, s.index("_skipped.append(_name)"))
        # 静态段继承主单点的 LVHAR（真空对齐要 LOCPOT），弛豫段关掉
        tree = ast.parse(s)
        fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_ionrelax_incars")
        ns = {"re": re, "IONRELAX_NSW": "60", "IONRELAX_EDIFFG": "-0.02"}
        exec(compile(ast.Module(body=[fn], type_ignores=[]), "gen", "exec"), ns)
        relax, stat = ns["_ionrelax_incars"](INCAR_NEW + "IBRION = -1\nNSW = 0\nLWAVE = .FALSE.\n")
        self.assertIn("LVHAR  = .TRUE.", stat)
        self.assertIn("LVHAR  = .FALSE.", relax)
        self.assertIn("ISTART = 1", stat)


if __name__ == "__main__":
    unittest.main(verbosity=2)
