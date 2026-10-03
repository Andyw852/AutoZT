#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_walltime.py —— V144：大 VASP 作业的墙时预估（S6 gen）与在跑作业的外推（tools/vasp_eta.py）。

用法：python test_walltime.py      退出码 0 = 全部 PASS。不需要 VASP。

起因：Mo2S3 S6（IBRION=6，20 原子，NFREE=4，ISYM 关）= 265 次 SCF，每步约 6.5 min ≈ 29 h，
在 regular 的 24 h 被杀（3918831 TIMEOUT 24:00:21）。
1) SLURM 时间格式、submit.sh 的 --time / --qos 推断、LOOP+ 解析、IBRION=6 步数、不等价原子数；
2) vasp_eta：用 Mo2S3 的实测（157 步 × 390 s）复现——regular 判超时、premium 来得及；非 IBRION=6 只报每步；
3) S6 gen 的 _walltime_estimate：写 walltime_estimate.json，超时给 ★ 告警，并且是 gen 输出的最后一段。
"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "tools"), str(ROOT / "step6_elastic")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402
import vasp_eta as E  # noqa: E402


def _poscar(n_mo=8, n_s=12, seed=0):
    rng = np.random.default_rng(seed)
    lat = [[9.0, 0, 0], [0, 6.0, 0], [0, 0, 20.0]]
    xyz = rng.random((n_mo + n_s, 3)) * [1, 1, 0.2] + [0, 0, 0.4]
    lines = ["Mo2S3", "1.0"] + ["  %.6f %.6f %.6f" % tuple(r) for r in lat] + ["Mo S", "%d %d" % (n_mo, n_s), "Direct"]
    lines += ["  %.6f %.6f %.6f" % tuple(r) for r in xyz]
    return "\n".join(lines) + "\n"


SI = """Si
1.0
0.0 2.715 2.715
2.715 0.0 2.715
2.715 2.715 0.0
Si
2
Direct
0 0 0
0.25 0.25 0.25
"""


def _outcar(n, t):
    return "".join("      LOOP+:  cpu time  %.2f: real time  %.2f\n" % (t - 5, t) for _ in range(n))


def _step(qos="regular", n_done=157, t=390.0, incar="IBRION = 6\nISIF = 3\nNFREE = 4\nISYM = 0\n", extra=""):
    d = Path(tempfile.mkdtemp())
    (d / "INCAR").write_text(incar)
    (d / "POSCAR").write_text(_poscar())
    (d / "OUTCAR").write_text(_outcar(n_done, t))
    (d / "submit.sh").write_text("#!/bin/bash\n#SBATCH --job-name=x\n#SBATCH --qos=%s\n%s" % (qos, extra))
    return d


class HelperTests(unittest.TestCase):
    def test_hms(self):
        for v, h in (("48:00:00", 48.0), ("2-00:00:00", 48.0), ("1-12", 36.0), ("90", 1.5), ("30:00", 0.5),
                     ("1-06:30", 30.5)):
            self.assertAlmostEqual(kc._hms_to_h(v), h, places=6, msg=v)
        self.assertIsNone(kc._hms_to_h("abc"))

    def test_limit(self):
        d = Path(tempfile.mkdtemp())
        f = d / "submit.sh"
        f.write_text("#SBATCH --qos=premium\n#SBATCH --time=36:00:00\n")
        self.assertEqual(kc.walltime_limit_h(f)[0], 36.0)                 # --time 优先
        f.write_text("#SBATCH --qos=premium\n")
        self.assertEqual(kc.walltime_limit_h(f)[0], 48.0)
        f.write_text("#SBATCH --qos=weird\n")
        h, why = kc.walltime_limit_h(f)
        self.assertEqual(h, 24.0)
        self.assertIn("weird", why)
        self.assertEqual(kc.walltime_limit_h(d / "nope")[0], 24.0)

    def test_loop_and_steps(self):
        d = Path(tempfile.mkdtemp())
        (d / "OUTCAR").write_text("junk\n" + _outcar(3, 100.0) + "LOOP:  cpu time 1: real time 2\n")
        self.assertEqual(kc.outcar_loop_times(d / "OUTCAR"), [100.0, 100.0, 100.0])
        self.assertEqual(kc.ibrion6_steps(20, 4), 265)                    # Mo2S3
        self.assertEqual(kc.ibrion6_steps(20, 2), 1 + 120 + 12)
        self.assertEqual(kc.ibrion6_steps(2, 2, isif=0), 13)

    def test_inequivalent(self):
        d = Path(tempfile.mkdtemp())
        (d / "POSCAR").write_text(SI)
        self.assertEqual(kc.n_inequivalent_atoms(d / "POSCAR"), 1)
        (d / "bad").write_text("x\n")
        self.assertIsNone(kc.n_inequivalent_atoms(d / "bad"))


class EtaTests(unittest.TestCase):
    def test_mo2s3_regular_times_out(self):
        r = E.eta(_step("regular"))
        self.assertEqual((r["total"], r["done"]), (265, 157))
        self.assertAlmostEqual(r["estimate_h"], 265 * 390 / 3600.0, places=1)   # ≈ 28.7 h
        self.assertTrue(r["over"])
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertEqual(E.main([str(_step("regular"))]), 1)
        self.assertIn("★", buf.getvalue())

    def test_premium_ok(self):
        r = E.eta(_step("premium"))
        self.assertFalse(r["over"])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(E.main([str(_step("premium"))]), 0)

    def test_dfpt_unknown_total(self):
        r = E.eta(_step(incar="IBRION = 8\nLEPSILON = .TRUE.\n"))
        self.assertIsNone(r["total"])
        self.assertNotIn("estimate_h", r)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertEqual(E.main([str(_step(incar="IBRION = 8\n"))]), 0)
        self.assertIn("总步数事先不知道", buf.getvalue())
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(E.main([tempfile.mkdtemp()]), 2)


class GenTests(unittest.TestCase):
    def _dirs(self, qos, loop_s):
        m = Path(tempfile.mkdtemp())
        s1, s2 = m / "step1_opt", m / "step6_elastic"
        s1.mkdir()
        s2.mkdir()
        (s1 / "OUTCAR").write_text(_outcar(20, loop_s))
        (s2 / "INCAR").write_text("IBRION = 6\nISIF = 3\nNFREE = 4\nISYM = 0\n")
        (s2 / "POSCAR").write_text(_poscar())
        (s2 / "submit.sh").write_text("#!/bin/bash\n#SBATCH --qos=%s\n" % qos)
        return s1, s2

    def test_estimate(self):
        import gen_step2_elastic as G
        s1, s2 = self._dirs("regular", 195.0)                          # ×2 = 390 s/步
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            rec = G._walltime_estimate(s1, s2)
        self.assertEqual(rec["scf_steps_upper_bound"], 265)
        self.assertTrue(rec["over"])
        self.assertIn("★", buf.getvalue())
        self.assertIn("submit.qos=premium", buf.getvalue())
        self.assertEqual(json.loads((s2 / "walltime_estimate.json").read_text())["scf_steps_upper_bound"], 265)
        s1, s2 = self._dirs("premium", 195.0)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertFalse(G._walltime_estimate(s1, s2)["over"])
        self.assertIn("[OK]", buf.getvalue())

    def test_no_s1_timing_and_wiring(self):
        import gen_step2_elastic as G
        s1, s2 = self._dirs("regular", 1.0)
        (s1 / "OUTCAR").unlink()
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            rec = G._walltime_estimate(s1, s2)
        self.assertNotIn("estimate_h", rec)
        self.assertIn("vasp_eta", buf.getvalue())
        src = (ROOT / "step6_elastic" / "gen_step2_elastic.py").read_text(encoding="utf-8")
        i_done = src.index('print("\\n[DONE] step2_elastic 已生成')
        self.assertLess(i_done, src.index("_walltime_estimate(step1, step2)"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
