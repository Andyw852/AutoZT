#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_elastic_decompose.py —— V154：tools/elastic_ionic_decompose.py（纯 numpy，合成 OUTCAR，不需要 VASP）。

用法：python test_elastic_decompose.py      退出码 0 = 全部 PASS。

模型：4 原子，力常数矩阵 K = U·diag(λ)·Uᵀ，U 的前 3 列是均匀平移（λ_T = 1e-4 eV/Å²，数值噪声级），
内应变 Λ = U·P。VASP 的离子弛豫项 = Σ −(1/Ω)PᵀP/λ（全部模）。OUTCAR 里 SECOND DERIVATIVES 按 −K 打印。
  A（Mo2S3 式）：平移模带一点 XY 耦合 -> XY 离子项被拉成大负数，TOTAL 不正定；扣掉平移模后正定 -> artifact；
  B（真实软模）：平移模不耦合，一个光学模 λ=0.05 强耦合 XY -> 扣掉平移模仍不正定 -> unstable；
  C：Λ 的剪切列按 ½ 打印（约定不同）-> 自动定出剪切因子 2，照样复现；
  D：离子项块被篡改 -> 复现不了 -> unreliable（退出码 2）；
  E：VASP 已去掉平移模（离子项 = 非平移模之和）-> 认出来，按真实值判；
  F：VASP 6 的实际格式（Mo2S3）：OUTCAR 无 SECOND DERIVATIVES / 无离子项块、内应变两段，力常数在 vasprun.xml。
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
sys.path.insert(0, str(ROOT / "tools"))
import elastic_ionic_decompose as E  # noqa: E402

N = 4
VOL = 300.0
CLAMPED = np.diag([600.0, 550.0, 5.0, 268.0, 1.0, 1.0])
CLAMPED[0, 1] = CLAMPED[1, 0] = 190.0


def _model(case):
    rng = np.random.default_rng(1)
    t = np.zeros((3 * N, 3))
    for ax in range(3):
        t[ax::3, ax] = 1.0 / np.sqrt(N)
    rest = rng.standard_normal((3 * N, 3 * N - 3))
    rest -= t @ (t.T @ rest)
    q, _ = np.linalg.qr(rest)
    u = np.hstack([t, q])
    lam = np.concatenate([[1e-4] * 3, np.linspace(2.0, 20.0, 3 * N - 3)])
    p = rng.standard_normal((3 * N, 6)) * 0.3
    p[:3] = 0.0
    if case == "A":
        p[0, 3] = 0.114                                    # 平移模 x 的 XY 残差 -> 约 −700 kBar
    if case == "B":
        lam[3] = 0.05
        p[3] = 0.0
        p[3, 3] = 2.6                                      # 软光学模强耦合 XY -> 约 −720 kBar
    k = u @ np.diag(lam) @ u.T
    lmat = u @ p
    ionic = sum(-np.outer(p[m], p[m]) / lam[m] for m in range(3 * N)) * E.EV_A3_TO_KBAR / VOL
    proj = sum(-np.outer(p[m], p[m]) / lam[m] for m in range(3, 3 * N)) * E.EV_A3_TO_KBAR / VOL
    return k, lmat, ionic, proj


def _block(title, m):
    out = [" %s (kBar)" % title, " Direction    XX          YY          ZZ          XY          YZ          ZX",
           " " + "-" * 80]
    out += [" %-5s" % lab + "".join("%12.4f" % v for v in row) for lab, row in zip(E.VASP_ORDER, m)]
    return out


def _outcar(k, lmat, ionic, shear_scale=1.0, tamper=False):
    lines = [" NIONS =      %d" % N, "  volume of cell :      %.4f" % VOL, "",
             " SECOND DERIVATIVES (NOT SYMMETRIZED)", " ------------------------------------",
             "        " + "".join("%10s" % ("%d%s" % (a + 1, x)) for a in range(N) for x in "XYZ")]
    labs = ["%d%s" % (a + 1, x) for a in range(N) for x in "XYZ"]
    for i, lab in enumerate(labs):
        lines.append("  %-5s" % lab + "".join("%12.6f" % v for v in -k[i]))
    lines.append("")
    lines += _block("ELASTIC MODULI", CLAMPED)
    lines += [""] + _block("SYMMETRIZED ELASTIC MODULI", CLAMPED) + [""]
    lm = lmat.copy()
    lm[:, 3:] *= shear_scale
    for a in range(N):
        lines += [" INTERNAL STRAIN TENSOR FOR ION    %d for displacements in x,y,z  (eV/Angst):" % (a + 1),
                  " Direction    XX          YY          ZZ          XY          YZ          ZX",
                  " " + "-" * 80]
        lines += [" %s   " % ax + "".join("%12.5f" % v for v in lm[3 * a + i]) for i, ax in enumerate("xyz")]
        lines.append("")
    ion = ionic.copy()
    if tamper:
        ion[3, 3] *= 0.3
    lines += _block("ELASTIC MODULI CONTR FROM IONIC RELAXATION", ion) + [""]
    lines += _block("TOTAL ELASTIC MODULI", CLAMPED + ion) + [""]
    return "\n".join(lines) + "\n"


def _run(case, **kw):
    k, lmat, ionic, proj = _model(case)
    d = Path(tempfile.mkdtemp())
    s6 = d / "step6_elastic"
    s6.mkdir()
    (s6 / "OUTCAR").write_text(_outcar(k, lmat, ionic, **kw))
    with contextlib.redirect_stdout(io.StringIO()) as buf, contextlib.redirect_stderr(io.StringIO()):
        rc = E.main([str(d), "--dim", "2d"])
    js = s6 / "elastic_ionic_decompose.json"
    r = json.loads(js.read_text(encoding="utf-8")) if js.is_file() else None
    return rc, r, buf.getvalue(), ionic, proj


MASS = [95.95, 95.95, 32.06, 32.06]                        # 2 Mo + 2 S


def _vasp6(k, lmat, ionic):
    """VASP 6 / IBRION=6 + ISIF=3 的实际样子（Mo2S3）：OUTCAR 没有 SECOND DERIVATIVES、没有离子项块；
    内应变分 FROM STRAINED CELLS / FROM DISPLACED ATOMS 两段（表头 X Y Z XY YZ ZX）；力常数只在 vasprun.xml 的
    <dynmat> hessian 里（质量加权、取负）。这里让 DISPLACED ATOMS 是真值，STRAINED CELLS 偏 10%。"""
    lines = [" NIONS =      %d" % N, "  volume of cell :      %.4f" % VOL, ""]
    lines += _block("ELASTIC MODULI", CLAMPED) + [""] + _block("SYMMETRIZED ELASTIC MODULI", CLAMPED) + [""]
    for title, lm in (("INTERNAL STRAIN TENSORS FROM STRAINED CELLS", lmat * 1.1),
                      ("INTERNAL STRAIN TENSORS FROM DISPLACED ATOMS", lmat)):
        lines += [" " + title, " " + "-" * 60]
        for a in range(N):
            lines += [" INTERNAL STRAIN TENSOR FOR ION    %d for displacements in x,y,z  (eV/Angst):" % (a + 1),
                      "          X           Y           Z          XY          YZ          ZX"]
            lines += [" %s " % ax + "".join("%12.5f" % v for v in lm[3 * a + i]) for i, ax in enumerate("xyz")]
            lines.append("")
    lines += _block("TOTAL ELASTIC MODULI", CLAMPED + ionic) + [""]
    m = np.repeat(MASS, 3)
    hess = -k / np.sqrt(np.outer(m, m))
    xml = ['<?xml version="1.0" encoding="ISO-8859-1"?>', "<modeling>", " <atominfo>",
           '  <array name="atoms"><set>']
    xml += ["   <rc><c>%s</c><c>%d</c></rc>" % (el, t) for el, t in (("Mo", 1), ("Mo", 1), ("S", 2), ("S", 2))]
    xml += ["  </set></array>", '  <array name="atomtypes"><set>',
            "   <rc><c>2</c><c>Mo</c><c>%.4f</c><c>14.0</c><c>PAW_PBE Mo_pv</c></rc>" % MASS[0],
            "   <rc><c>2</c><c>S</c><c>%.4f</c><c>6.0</c><c>PAW_PBE S</c></rc>" % MASS[2],
            "  </set></array>", " </atominfo>", " <calculation>", "  <dynmat>", '   <varray name="hessian">']
    xml += ["    <v>" + " ".join("%.10e" % v for v in row) + "</v>" for row in hess]
    xml += ["   </varray>", "  </dynmat>", " </calculation>", "</modeling>"]
    return "\n".join(lines) + "\n", "\n".join(xml) + "\n"


class DecomposeTests(unittest.TestCase):
    def test_translation_artifact(self):
        rc, r, out, ionic, proj = _run("A")
        self.assertEqual(rc, 0, out)
        self.assertTrue(r["reproduced"], r)
        self.assertEqual(r["n_translation_modes"], 3)
        self.assertEqual(r["worst_component"], "XY")
        self.assertFalse(r["total_pd"])
        self.assertTrue(r["projected_pd"])
        self.assertEqual(r["verdict"], "artifact")
        self.assertLess(ionic[3, 3] + CLAMPED[3, 3], 0)                     # 模型里 TOTAL 确实翻负了
        got = np.array(r["projected_total_vasp_order_kbar"])
        np.testing.assert_allclose(got, CLAMPED + proj, atol=0.05)
        self.assertIn("数值假象", out)

    def test_real_soft_mode(self):
        rc, r, out, _i, _p = _run("B")
        self.assertEqual(rc, 1, out)
        self.assertEqual(r["verdict"], "unstable")
        self.assertAlmostEqual(r["translation_part"], 0.0, delta=0.5)
        self.assertEqual(r["top_other_modes"][0]["eigenvalue_eV_A2"], 0.05)
        self.assertIn("真实的内应变耦合", out)

    def test_shear_convention_calibrated(self):
        rc, r, _o, _i, _p = _run("A", shear_scale=0.5)
        self.assertEqual(rc, 0)
        self.assertEqual(r["shear_factor"], 2.0)
        self.assertTrue(r["reproduced"])

    def test_vasp6_format_mo2s3(self):
        """Mo2S3 的实际格式：力常数从 vasprun.xml 拿，内应变选 DISPLACED ATOMS，离子项 = TOTAL − SYMMETRIZED。"""
        k, lmat, ionic, proj = _model("A")
        d = Path(tempfile.mkdtemp())
        s6 = d / "step6_elastic"
        s6.mkdir()
        oc, vr = _vasp6(k, lmat, ionic)
        (s6 / "OUTCAR").write_text(oc)
        (s6 / "vasprun.xml").write_text(vr)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            rc = E.main([str(d), "--dim", "2d"])
        r = json.loads((s6 / "elastic_ionic_decompose.json").read_text(encoding="utf-8"))
        self.assertEqual(rc, 0, buf.getvalue())
        self.assertEqual(r["force_constants_from"], "vasprun hessian × √(m_i m_j)")
        self.assertEqual(r["internal_strain_from"], "displaced atoms")
        self.assertEqual(r["ionic_reference"], "TOTAL − SYMMETRIZED")
        self.assertEqual(r["verdict"], "artifact")
        np.testing.assert_allclose(np.array(r["projected_total_vasp_order_kbar"]), CLAMPED + proj, atol=0.05)
        (s6 / "vasprun.xml").unlink()                                       # 没有力常数 -> 明确报错
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(E.main([str(d), "--dim", "2d"]), 2)
        self.assertIn("拿不到力常数矩阵", err.getvalue())

    def test_vasp_already_projected(self):
        """VASP 若已去掉平移模：它的离子项 = 非平移模之和；本工具认出来，结论按真实值（这里仍正定 -> ok）。"""
        k, lmat, _ionic, proj = _model("A")
        d = Path(tempfile.mkdtemp())
        (d / "step6_elastic").mkdir()
        (d / "step6_elastic" / "OUTCAR").write_text(_outcar(k, lmat, proj))
        with contextlib.redirect_stdout(io.StringIO()):
            rc = E.main([str(d), "--dim", "2d"])
        r = json.loads((d / "step6_elastic" / "elastic_ionic_decompose.json").read_text(encoding="utf-8"))
        self.assertEqual(rc, 0)
        self.assertFalse(r["vasp_includes_translations"])
        self.assertTrue(r["reproduced"])
        self.assertEqual(r["verdict"], "ok")

    def test_cannot_reproduce(self):
        rc, r, out, _i, _p = _run("A", tamper=True)
        self.assertEqual(rc, 2)
        self.assertEqual(r["verdict"], "unreliable")
        self.assertIn("复现不了", out)

    def test_missing_inputs(self):
        d = Path(tempfile.mkdtemp())
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(E.main([str(d)]), 2)
        (d / "OUTCAR").write_text(" NIONS = 2\n volume of cell : 10.0\n")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(E.main([str(d)]), 2)
        self.assertIn("SECOND DERIVATIVES", err.getvalue())


if __name__ == "__main__":
    unittest.main(verbosity=2)
