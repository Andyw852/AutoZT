#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_edge_mass_check.py —— V176：tools/edge_mass_check.py（带边细 k 点非自洽量 m*，核对 S8.2 的 m* 收敛没有）。纯本地。

模型能带：六方 2D（a = 3.04 Å），带边在 K；E = Eg + ħ²q²/2m + w q³ cos3θ + c₄q⁴（导带），价带取负。
三角翘曲 w q³cos3θ 是奇次项，±q 平均后应完全消掉；q⁴ 由拟合吸收。

用法：python test_edge_mass_check.py      退出码 0 = 全部 PASS。
"""
import contextlib
import io
import json
import math
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT.parent / "_common" / "opt"))
sys.path.insert(0, str(_ROOT / "tools"))
import ke_common as kc  # noqa: E402
import edge_mass_check as E  # noqa: E402

A, C = 3.04, 20.0
LAT = np.array([[A, 0, 0], [-A / 2, A * math.sqrt(3) / 2, 0], [0, 0, C]])
REC = 2 * np.pi * np.linalg.inv(LAT).T
K = np.array([1 / 3, 1 / 3, 0.0])
ME, MH, EG = 0.8663, 0.8829, 0.94


def _model(kf, me=ME, mh=MH, w=0.6, c4=-3.0):
    d = (np.asarray(kf) - K) @ REC
    q, th = math.hypot(d[0], d[1]), math.atan2(d[1], d[0])
    warp = w * q ** 3 * math.cos(3 * th)
    return (-(E.HB2M * q * q / mh) - warp - c4 * q ** 4,          # 价带
            EG + E.HB2M * q * q / me + warp + c4 * q ** 4)       # 导带


def _vasprun(kpts, eig, occ, decoy=True):
    ks = "\n".join("    <v> %.8f %.8f %.8f </v>" % tuple(k) for k in kpts)
    def block(e, o):
        return "\n".join("     <set comment=\"kpoint %d\">\n%s\n     </set>" % (
            i + 1, "\n".join("      <r> %10.4f %8.4f </r>" % (e[i][b], o[i][b]) for b in range(len(e[i]))))
            for i in range(len(e)))
    eig_xml = ("  <eigenvalues>\n   <array>\n    <set>\n     <set comment=\"spin 1\">\n%s\n     </set>\n    </set>\n"
               "   </array>\n  </eigenvalues>\n" % block(eig, occ))
    proj = ("  <projected>\n" + eig_xml.replace(" %10.4f" % eig[0][0], " %10.4f" % 99.0) + "  </projected>\n") if decoy else ""
    return ('<?xml version="1.0" encoding="ISO-8859-1"?>\n<modeling>\n <kpoints>\n  <varray name="kpointlist" >\n'
            + ks + "\n  </varray>\n </kpoints>\n <calculation>\n" + proj + eig_xml + " </calculation>\n</modeling>\n")


def _material(dim="2D"):
    m = Path(tempfile.mkdtemp()) / "CrS2_hex"
    s3 = m / "step3_uniform"
    s3.mkdir(parents=True)
    (s3 / "POSCAR").write_text("CrS2\n1.0\n" + "\n".join("%.10f %.10f %.10f" % tuple(r) for r in LAT)
                               + "\nCr S\n1 2\nDirect\n0 0 0.5\n0.3333333 0.6666667 0.57\n0.3333333 0.6666667 0.43\n")
    (s3 / "INCAR").write_text("SYSTEM = CrS2 uniform\nENCUT = 520\nEDIFF = 1E-6\nALGO = Fast\nISYM = 2\n"
                              "LVHAR = .TRUE.\nKPAR = 4\nNCORE = 6\n")
    (s3 / "KPOINTS").write_text("auto\n0\nGamma\n 48 48 3\n")
    (s3 / "POTCAR").write_text("PAW_PBE Cr\nPAW_PBE S\n")
    (s3 / "CHGCAR").write_text("chg\n")
    (s3 / kc.METHOD_FILE).write_text("DIM=%s\n" % dim)
    kpts = [[0, 0, 0], list(K), [0.5, 0, 0]]
    eig = [[-1.0, 2.0], list(_model(K)), [-0.8, 1.6]]                  # 模型只在 K 附近成立：Γ、M 直接给远离带边的值
    (s3 / "vasprun.xml").write_text(_vasprun(kpts, eig, [[1.0, 0.0]] * 3))
    (m / "step8.2_dpt").mkdir()
    return m


def _dpt(m, me, mh):
    (m / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": "electron", "inputs": {"m_eff_m0": me}}, {"carrier": "hole", "inputs": {"m_eff_m0": mh}}]}))


def _eigenval(out, vasp6=True, **model):
    ln = (out / "KPOINTS").read_text().splitlines()
    kf = [[float(x) for x in l.split()[:3]] for l in ln[3:3 + int(ln[1])]]
    rows = ["    3    3    1    1", "  0.1E+02", "  1.0E-04", "  CAR", " CrS2", "     18  %d  2" % len(kf)]
    for k in kf:
        ev, ec = _model(k, **model)
        rows += ["", "  %.7E  %.7E  %.7E  %.7E" % (k[0], k[1], k[2], 1.0)]
        rows += ["    1  %.6f%s" % (ev, "  1.000000" if vasp6 else ""), "    2  %.6f%s" % (ec, "  0.000000" if vasp6 else "")]
    (out / "EIGENVAL").write_text("\n".join(rows) + "\n")


def _run(argv):
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        rc = E.main(argv)
    return rc, buf.getvalue()


class GeometryAndFitTests(unittest.TestCase):
    def test_star_geometry(self):
        qs = [0.005 * (i + 1) for i in range(10)]
        pts = E.star_points(REC, K, qs, 6)
        self.assertEqual(len(pts), 1 + 6 * 2 * 10)
        for kf, j, s, q in pts[1:]:
            d = (kf - K) @ REC
            self.assertAlmostEqual(np.linalg.norm(d), q, places=10)      # 笛卡尔距离 = q（含 2π）
            self.assertAlmostEqual(kf[2], 0.0, places=12)                # 只在面内
        th = sorted({round(math.degrees(math.atan2(*(((kf - K) @ REC)[[1, 0]]))), 6) % 180
                     for kf, j, s, q in pts[1:]})
        self.assertEqual(th, [0.0, 30.0, 60.0, 90.0, 120.0, 150.0])

    def test_warping_cancels_and_mass_recovered(self):
        qs = [0.005 * (i + 1) for i in range(10)]
        pts = E.star_points(REC, K, qs, 6)
        ev = np.array([_model(p[0])[0] for p in pts])
        ec = np.array([_model(p[0])[1] for p in pts])
        for e, carrier, m in ((ec, "electron", ME), (ev, "hole", MH)):
            dirs = E.fit_directions(e, pts, 6, carrier, 0.05)
            for d in dirs:
                self.assertAlmostEqual(d["m"] / m, 1.0, delta=1e-6, msg=d)
            self.assertAlmostEqual(E.summarize(dirs)["anisotropy"], 1.0, places=6)
        bad = E.fit_directions(-ec, pts, 6, "electron", 0.05)            # 极大值当导带底 -> 报错
        self.assertTrue(all("error" in d for d in bad))

    def test_anisotropic(self):
        qs = [0.005 * (i + 1) for i in range(10)]
        pts = E.star_points(REC, K, qs, 6)
        e = np.array([E.HB2M * ((((p[0] - K) @ REC)[0]) ** 2 / 0.8 + (((p[0] - K) @ REC)[1]) ** 2 / 1.0)
                      for p in pts])
        dirs = E.fit_directions(e, pts, 6, "electron", 0.05)
        self.assertAlmostEqual(dirs[0]["m"], 0.8, places=6)               # θ = 0 沿笛卡尔 x
        self.assertAlmostEqual(dirs[3]["m"], 1.0, places=6)               # θ = 90°
        self.assertAlmostEqual(E.summarize(dirs)["anisotropy"], 1.25, places=6)


class ParserTests(unittest.TestCase):
    def test_vasprun_ignores_projected(self):
        m = _material()
        kf, eig, occ = E.read_vasprun_eigen(m / "step3_uniform" / "vasprun.xml")
        self.assertEqual(eig.shape, (1, 3, 2))
        self.assertAlmostEqual(eig[0, 0, 0], -1.0, places=6)              # 不是 <projected> 里 decoy 的 99
        ed = E.band_edges(eig, occ)
        self.assertEqual((ed["vbm"][1], ed["cbm"][1]), (1, 1))            # 都在 K
        self.assertIsNone(E.band_edges(eig, np.ones_like(occ)))           # 全占据：没有带隙

    def test_eigenval_vasp5_vasp6(self):
        m = _material()
        out = m.parent / "emc"
        _run(["make", str(m), str(out)])
        for v6 in (True, False):
            _eigenval(out, vasp6=v6)
            kf, eig = E.read_eigenval(out / "EIGENVAL")
            self.assertEqual(eig.shape, (1, 121, 2))
            self.assertAlmostEqual(eig[0, 0, 1], round(EG, 6), places=6)


class MakeFitTests(unittest.TestCase):
    def test_make_inputs(self):
        m = _material()
        out = m.parent / "emc"
        rc, txt = _run(["make", str(m), str(out)])
        self.assertEqual(rc, 0, txt)
        inc = kc.parse_incar((out / "INCAR").read_text())
        self.assertEqual((inc["ICHARG"], inc["ISYM"], inc["NBANDS"], inc["EDIFF"], inc["ALGO"]),
                         ("11", "0", "2", "1E-8", "Normal"))
        self.assertNotIn("LVHAR", inc)
        self.assertEqual(inc["KPAR"], "4")                                # 并行设置照抄
        self.assertEqual(int((out / "KPOINTS").read_text().splitlines()[1]), 121)   # CBM、VBM 同在 K：一个星
        for f in ("CHGCAR", "POTCAR", "POSCAR"):
            self.assertFalse((out / f).is_symlink())
            self.assertNotEqual(os.stat(out / f).st_ino, os.stat(m / "step3_uniform" / f).st_ino)  # 不是硬链接
        meta = json.loads((out / E.META).read_text())
        self.assertEqual(meta["s3_mesh"], [48, 48, 3])
        self.assertEqual(meta["stars"][0]["edges"], ["cbm", "vbm"])

    def test_make_refuses(self):
        m = _material()
        with self.assertRaises(SystemExit):
            _run(["make", str(m), str(m / "emc")])                       # 不许放进材料目录
        out = m.parent / "busy"
        out.mkdir()
        (out / "x").write_text("x")
        with self.assertRaises(SystemExit):
            _run(["make", str(m), str(out)])                             # 不许写进非空目录
        with self.assertRaises(SystemExit):
            _run(["make", str(_material(dim="3D")), str(m.parent / "e3")])

    def test_fit_pass_and_fail(self):
        m = _material()
        _dpt(m, ME, MH)
        out = m.parent / "emc"
        _run(["make", str(m), str(out)])
        _eigenval(out)
        rc, txt = _run(["fit", str(out)])
        self.assertEqual(rc, 0, txt)
        res = json.loads((out / E.RESULT).read_text())
        for c, mm in (("electron", ME), ("hole", MH)):
            self.assertAlmostEqual(res["carriers"][c]["summary"]["m_light2"] / mm, 1.0, delta=2e-3)
            self.assertLess(abs(res["carriers"][c]["half_window_change"]), 0.02)
        _dpt(m, ME * 1.10, MH)                                           # S8.2 拟重 10%（粗网格那种）
        out2 = m.parent / "emc2"
        _run(["make", str(m), str(out2)])
        _eigenval(out2)
        rc, txt = _run(["fit", str(out2)])
        self.assertEqual(rc, 1)
        self.assertIn("S3 网格对 m* 还不够", txt)

    def test_fit_wrong_run(self):
        m = _material()
        out = m.parent / "emc"
        _run(["make", str(m), str(out)])
        (out / "EIGENVAL").write_text("    3    3    1    1\n x\n x\n x\n x\n 18 1 2\n\n 0 0 0 1\n 1 -1.0\n 2 1.0\n")
        rc, txt = _run(["fit", str(out)])
        self.assertEqual(rc, 2)
        self.assertIn("不是这次 make 的输入", txt)


if __name__ == "__main__":
    unittest.main(verbosity=2)
