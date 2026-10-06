#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_s3_grid_gate.py —— V173：S3 网格不合现行规则（旧规则下生成）时，S8/S8.4 拒绝生成（纯本地）。

起因（2026-10-07）：CrS₂_hex 的 S3 一直是 08-26 的 15×15×1（面内间距 0.16 Å⁻¹ = 规则 0.05 的 3.2 倍、kz 一层），
CrSe₂ 是 48×48×3。S7 照抄 S3 网格、S8.4 照用，谁都不查；两者 S8.4 的 ADP/DPT 一个 1.7、一个 0.5。

用法：python test_s3_grid_gate.py      退出码 0 = 全部 PASS。
"""
import math
import sys
import tempfile
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
import ke_common as kc  # noqa: E402


def _mat(a, c, mesh, dim="2D", hexagonal=True):
    d = Path(tempfile.mkdtemp())
    s3 = d / "step3_uniform"
    s3.mkdir()
    if hexagonal:
        lat = "%f 0 0\n%f %f 0\n0 0 %f\n" % (a, -a / 2, a * math.sqrt(3) / 2, c)
    else:
        lat = "%f 0 0\n0 %f 0\n0 0 %f\n" % (a, a, c)
    (s3 / "POSCAR").write_text("x\n1.0\n" + lat + "Cr S\n1 2\nDirect\n0 0 0.5\n0.3333 0.6667 0.57\n0.3333 0.6667 0.43\n")
    (s3 / "KPOINTS").write_text("auto\n0\nGamma\n %d %d %d\n" % tuple(mesh))
    (s3 / kc.METHOD_FILE).write_text("DIM=%s\n" % dim)
    return d


class S3GridTests(unittest.TestCase):
    def test_crs2_old_grid_is_rejected(self):
        iss = kc.s3_grid_issues(_mat(3.04, 20, (15, 15, 1)))
        self.assertEqual([lv for lv, _ in iss], ["error", "error", "error"])
        self.assertIn("3.2 倍", iss[0][1])
        self.assertIn("kz = 1", iss[2][1])

    def test_current_grid_passes(self):
        self.assertEqual(kc.s3_grid_issues(_mat(3.21, 20, (48, 48, 3))), [])

    def test_kz1_alone_is_rejected(self):
        iss = kc.s3_grid_issues(_mat(3.04, 20, (48, 48, 1)))
        self.assertEqual([lv for lv, _ in iss], ["error"])

    def test_slightly_coarse_only_warns(self):
        iss = kc.s3_grid_issues(_mat(3.04, 20, (42, 42, 3)))          # 0.057 Å⁻¹
        self.assertEqual([lv for lv, _ in iss], ["warn", "warn"])

    def test_3d_rule(self):
        d = _mat(5.43, 5.43, (8, 8, 8), dim="3D", hexagonal=False)    # 0.145 Å⁻¹ > 2×0.06
        self.assertEqual([lv for lv, _ in kc.s3_grid_issues(d)], ["error"] * 3)
        d = _mat(5.43, 5.43, (20, 20, 20), dim="3D", hexagonal=False)
        self.assertEqual(kc.s3_grid_issues(d), [])

    def test_missing_s3_is_silent(self):
        self.assertEqual(kc.s3_grid_issues(Path(tempfile.mkdtemp())), [])

    def test_gate_exits_unless_allowed(self):
        d = _mat(3.04, 20, (15, 15, 1))
        with self.assertRaises(SystemExit):
            kc.s3_grid_gate(d, "S8.4")
        self.assertEqual(len(kc.s3_grid_gate(d, "S8.4", allow=True)), 3)
        self.assertEqual(kc.s3_grid_gate(_mat(3.21, 20, (48, 48, 3)), "S8.4"), [])


class WiringTests(unittest.TestCase):
    def test_gens(self):
        for f, label in (("step8.4_amset2d/gen_step14_amset2d.py", "S8.4"), ("step8_amset/gen_step10_amset.py", "S8")):
            s = (_ROOT / f).read_text(encoding="utf-8")
            self.assertIn('"ALLOW_COARSE_S3": (False, "bool")', s, f)
            self.assertIn('ALLOW_COARSE_S3 = bool(_p["ALLOW_COARSE_S3"])', s, f)
            i = s.index('kc.s3_grid_gate(cwd, "%s"' % label)
            self.assertLess(s.index('ALLOW_COARSE_S3 = bool(_p["ALLOW_COARSE_S3"])'), i, f)
            self.assertLess(i, s.index('kc.check_lineage(cwd, enabled=STRUCTURE_GUARD, label="%s")' % label), f)

    def test_s3_defaults_and_s7_warning(self):
        for f in ("step3_uniform/gen_step5_uniform.py", "step3b_uniform_full/gen_step5b_uniform_full.py"):
            self.assertRegex((_ROOT / f).read_text(encoding="utf-8"), r"(?m)^VACUUM_KZ_MIN = 3\b", f)
        s7 = (_ROOT / "step7_deform" / "gen_step9_deform.py").read_text(encoding="utf-8")
        self.assertLess(s7.index("_need = list(_s3)"), s7.index("kc.s3_grid_issues(out.parent, dim)"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
