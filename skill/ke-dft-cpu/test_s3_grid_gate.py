#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_s3_grid_gate.py —— V173：S3 网格不合现行规则（旧规则下生成）时，S8/S8.4 拒绝生成（纯本地）。

起因（2026-10-07）：CrS₂_hex 的 S3 一直是 08-26 的 15×15×1（面内间距 0.16 Å⁻¹ = 规则 0.05 的 3.2 倍、kz 一层），
CrSe₂ 是 48×48×3。S7 照抄 S3 网格、S8.4 照用，谁都不查；两者 S8.4 的 ADP/DPT 一个 1.7、一个 0.5。
V174：S8.2（m* 直接拟合 S3 网格点）同样拦、S8.1（BoltzTraP2）告警；两步只查面内。CrS₂ 换 48×48×3 后
m* 0.967/1.009 -> 0.866/0.883、DPT +25%/+31%。
V175：S7 和 S3 不是同一张网格且 S7 不合规（S3 重算后 S7 没跟着重跑）-> S8/S8.4 拦、S7.1 告警；
tools/lineage_check.py 查已经算完的旧结果。

用法：python test_s3_grid_gate.py      退出码 0 = 全部 PASS。
"""
import contextlib
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
sys.path.insert(0, str(_ROOT / "step8.2_dpt"))
sys.path.insert(0, str(_ROOT / "step8.1_boltztrap"))
sys.path.insert(0, str(_ROOT / "tools"))
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

    def test_inplane_only(self):                                       # [V174] S8.1 / S8.2
        self.assertEqual(kc.s3_grid_issues(_mat(3.04, 20, (48, 48, 1)), check_kz=False), [])
        iss = kc.s3_grid_issues(_mat(3.04, 20, (15, 15, 1)), check_kz=False)
        self.assertEqual([lv for lv, _ in iss], ["error", "error"])

    def test_missing_s3_is_silent(self):
        self.assertEqual(kc.s3_grid_issues(Path(tempfile.mkdtemp())), [])

    def test_gate_exits_unless_allowed(self):
        d = _mat(3.04, 20, (15, 15, 1))
        with self.assertRaises(SystemExit):
            kc.s3_grid_gate(d, "S8.4")
        self.assertEqual(len(kc.s3_grid_gate(d, "S8.4", allow=True)), 3)
        self.assertEqual(kc.s3_grid_gate(_mat(3.21, 20, (48, 48, 3)), "S8.4"), [])


def _s7(d, mesh, dirs=("undeformed", "deform-01", "deform-02"), ionrelax_mesh=None):
    for n in dirs:
        (d / "step7_deform" / n).mkdir(parents=True, exist_ok=True)
        (d / "step7_deform" / n / "KPOINTS").write_text("auto\n0\nGamma\n %d %d %d\n" % tuple(mesh))
    if ionrelax_mesh:
        ir = d / "step7_deform" / dirs[1] / "ionrelax"
        ir.mkdir(parents=True, exist_ok=True)
        (ir / "KPOINTS").write_text("auto\n0\nGamma\n %d %d %d\n" % tuple(ionrelax_mesh))
    return d


class S7GridTests(unittest.TestCase):                                   # [V175]
    def test_s3_redone_but_s7_stale(self):
        d = _s7(_mat(3.04, 20, (48, 48, 3)), (15, 15, 1))
        iss = kc.s7_grid_issues(d)
        self.assertEqual([lv for lv, _ in iss], ["error"])
        self.assertIn("15×15×1", iss[0][1])
        self.assertIn("S7 没跟着重跑", iss[0][1])
        self.assertEqual(kc.s3_grid_issues(d), [])                     # S3 闸门本身是通过的
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as cm:
            kc.s3_grid_gate(d, "S8.4", with_s7=True)
        self.assertIn("retry step7_deform", str(cm.exception.code))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(len(kc.s3_grid_gate(d, "S8.4", allow=True, with_s7=True)), 1)
            self.assertEqual(kc.s3_grid_gate(d, "S8.2", check_kz=False), [])          # S8.2 不看 S7

    def test_same_or_compliant_mesh(self):
        self.assertEqual(kc.s7_grid_issues(_s7(_mat(3.04, 20, (48, 48, 3)), (48, 48, 3))), [])
        iss = kc.s7_grid_issues(_s7(_mat(3.04, 20, (48, 48, 3)), (54, 54, 3)))
        self.assertEqual([lv for lv, _ in iss], ["warn"])
        with contextlib.redirect_stdout(io.StringIO()):
            kc.s3_grid_gate(_s7(_mat(3.04, 20, (48, 48, 3)), (54, 54, 3)), "S8.4", with_s7=True)   # 不退出

    def test_ionrelax_and_missing(self):
        d = _s7(_mat(3.04, 20, (48, 48, 3)), (48, 48, 3), ionrelax_mesh=(15, 15, 1))
        iss = kc.s7_grid_issues(d)
        self.assertEqual(len(iss), 1)
        self.assertIn("deform-01/ionrelax", iss[0][1])
        self.assertEqual(kc.s7_grid_issues(_mat(3.04, 20, (48, 48, 3))), [])           # 没有 S7

    def test_coarse_s3_message_wins(self):
        d = _s7(_mat(3.04, 20, (15, 15, 1)), (15, 15, 1))                 # S7 照抄了同一张粗网格
        self.assertEqual(kc.s7_grid_issues(d), [])
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as cm:
            kc.s3_grid_gate(d, "S8.4", with_s7=True)
        self.assertIn("step3_uniform 的网格是旧规则下生成的", str(cm.exception.code))


class GridLineageTests(unittest.TestCase):                              # [V175]
    def test_grid_lineage(self):
        probs = kc.grid_lineage(_mat(3.04, 20, (15, 15, 1)))
        self.assertEqual([st for st, _ in probs], ["step3_uniform"] * 3)
        self.assertEqual(sum("S8/S8.4/S8.1/S8.2" in w for _, w in probs), 2)     # 面内两条
        self.assertEqual(sum("AMSET 结果受影响" in w for _, w in probs), 1)       # kz 一条
        self.assertEqual(kc.grid_lineage(_mat(3.04, 20, (42, 42, 3))), [])        # 略粗只告警，不算问题
        probs = kc.grid_lineage(_s7(_mat(3.04, 20, (48, 48, 3)), (15, 15, 1)))
        self.assertEqual([st for st, _ in probs], ["step7_deform"])

    def test_cli_flags_old_results(self):
        import lineage_check as L
        root = Path(tempfile.mkdtemp())
        for name, s3, s7 in (("CrS2_hex", (48, 48, 3), (15, 15, 1)), ("CrSe2_hex", (48, 48, 3), (48, 48, 3))):
            m = _s7(_mat(3.04, 20, s3), s7)
            (m / "step1_opt").mkdir()
            m.rename(root / name)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertEqual(L.main([str(root)]), 1)
        out = buf.getvalue()
        self.assertIn("[★]   %s" % (root / "CrS2_hex"), out)
        self.assertIn("[OK]  %s" % (root / "CrSe2_hex"), out)
        self.assertIn("retry S7 -> S7.1", out)


class DptGateTests(unittest.TestCase):                                  # [V174]
    def setUp(self):
        import importlib
        import gen_step12_dpt as D
        self.D = importlib.reload(D)

    def test_coarse_s3_blocks_dpt(self):
        d = _mat(3.04, 20, (15, 15, 1))
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as cm:
            self.D._s3_grid_gate(d, "2d")
        self.assertIn("m* 直接拟合", str(cm.exception.code))
        self.assertIn("ALLOW_COARSE_S3", str(cm.exception.code))
        self.assertEqual(self.D._s3_grid_gate(_mat(3.04, 20, (48, 48, 1)), "2d"),    # kz 一层不拦 DPT
                         {"mesh": [48, 48, 1], "issues": [], "allow_coarse": False})

    def test_allow_via_step_conf(self):
        d = _mat(3.04, 20, (15, 15, 1))
        (d / "step.conf").write_text("[params]\nALLOW_COARSE_S3 = true\n")
        with contextlib.redirect_stdout(io.StringIO()):
            rec = self.D.apply_conf(d)
            g = self.D._s3_grid_gate(d, "2d")
        self.assertTrue(rec["ALLOW_COARSE_S3"]["value"])
        self.assertEqual((g["mesh"], len(g["issues"]), g["allow_coarse"]), ([15, 15, 1], 2, True))

    def test_boltztrap_warns_only(self):
        import gen_step11_boltztrap as B
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            g = B._s3_grid_note(_mat(3.04, 20, (15, 15, 1)), "2d")
        self.assertEqual((g["mesh"], len(g["issues"])), ([15, 15, 1], 2))
        self.assertIn("BoltzTraP2 插值的就是这张网格", buf.getvalue())


class WiringTests(unittest.TestCase):
    def test_gens(self):
        for f, label in (("step8.4_amset2d/gen_step14_amset2d.py", "S8.4"), ("step8_amset/gen_step10_amset.py", "S8")):
            s = (_ROOT / f).read_text(encoding="utf-8")
            self.assertIn('"ALLOW_COARSE_S3": (False, "bool")', s, f)
            self.assertIn('ALLOW_COARSE_S3 = bool(_p["ALLOW_COARSE_S3"])', s, f)
            i = s.index('kc.s3_grid_gate(cwd, "%s"' % label)
            self.assertLess(s.index('ALLOW_COARSE_S3 = bool(_p["ALLOW_COARSE_S3"])'), i, f)
            self.assertLess(i, s.index('kc.check_lineage(cwd, enabled=STRUCTURE_GUARD, label="%s")' % label), f)

    def test_dpt_and_boltztrap_wiring(self):                           # [V174]
        s = (_ROOT / "step8.2_dpt" / "gen_step12_dpt.py").read_text(encoding="utf-8")
        body = s[s.index("\ndef main("):]
        self.assertIn('"ALLOW_COARSE_S3": (False, "bool")', s)
        self.assertLess(body.index("overrides = apply_conf(cwd)"), body.index("s3_grid = _s3_grid_gate(cwd, dim)"))
        self.assertLess(body.index("s3_grid = _s3_grid_gate(cwd, dim)"), body.index('"results": [_one_carrier('))
        self.assertIn('"s3_grid": s3_grid', body)
        b = (_ROOT / "step8.1_boltztrap" / "gen_step11_boltztrap.py").read_text(encoding="utf-8")
        body = b[b.index("\ndef main("):]
        self.assertLess(body.index("s3_grid = _s3_grid_note(cwd, dim)"), body.index("run_boltztrap_crta("))
        self.assertIn('res["s3_grid"] = s3_grid', body)

    def test_s7_wiring(self):                                          # [V175]
        for f, label in (("step8.4_amset2d/gen_step14_amset2d.py", "S8.4"), ("step8_amset/gen_step10_amset.py", "S8")):
            s = (_ROOT / f).read_text(encoding="utf-8")
            i = s.index('kc.s3_grid_gate(cwd, "%s"' % label)
            self.assertIn("with_s7=True", s[i:i + 120], f)
        s = (_ROOT / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py").read_text(encoding="utf-8")
        self.assertLess(s.index("kc.s7_grid_issues(cwd"), s.index("amset deform read undeformed"))
        s = (_ROOT / "tools" / "lineage_check.py").read_text(encoding="utf-8")
        self.assertIn("kc.grid_lineage(mat)", s)

    def test_s3_defaults_and_s7_warning(self):
        for f in ("step3_uniform/gen_step5_uniform.py", "step3b_uniform_full/gen_step5b_uniform_full.py"):
            self.assertRegex((_ROOT / f).read_text(encoding="utf-8"), r"(?m)^VACUUM_KZ_MIN = 3\b", f)
        s7 = (_ROOT / "step7_deform" / "gen_step9_deform.py").read_text(encoding="utf-8")
        self.assertLess(s7.index("_need = list(_s3)"), s7.index("kc.s3_grid_issues(out.parent, dim)"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
