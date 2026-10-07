#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_edge_mass_check.py —— V176：tools/edge_mass_check.py（带边细 k 点非自洽量 m*，核对 S8.2 的 m* 收敛没有）。纯本地。

模型能带：六方 2D（a = 3.04 Å），带边在 K；E = Eg + ħ²q²/2m + w q³ cos3θ + c₄q⁴（导带），价带取负。
三角翘曲 w q³cos3θ 是奇次项，±q 平均后应完全消掉；q⁴ 由拟合吸收。
V177：S3 的 CHGCAR 是 0 字节（出厂 LCHARG = .FALSE.）时走两段式（stage1 从 WAVECAR 自洽 -> stage2 非自洽）。
V178：每个载流子单独给"通过/不通过"和原因；窗口减半不稳时区分"能量噪声"（各方向离散大）和"高阶非抛物"（各方向一致）。
V184：scan —— 在一串窗口上逐方向做抛物线拟合；六方胞标出锯齿/扶手椅。

用法：python test_edge_mass_check.py      退出码 0 = 全部 PASS。
"""
import contextlib
import io
import json
import math
import os
import subprocess
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


def _model(kf, me=ME, mh=MH, w=0.6, c4=-3.0, c6h=0.0):
    d = (np.asarray(kf) - K) @ REC
    q, th = math.hypot(d[0], d[1]), math.atan2(d[1], d[0])
    warp = w * q ** 3 * math.cos(3 * th)
    return (-(E.HB2M * q * q / mh) - warp - c4 * q ** 4 + c6h * q ** 6,   # 价带（c6h：高阶非抛物）
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


S3_TOTEN = -20.12345678


def _material(dim="2D", chgcar="chg\n", wavecar="wave"):
    m = Path(tempfile.mkdtemp()) / "CrS2_hex"
    s3 = m / "step3_uniform"
    s3.mkdir(parents=True)
    (s3 / "POSCAR").write_text("CrS2\n1.0\n" + "\n".join("%.10f %.10f %.10f" % tuple(r) for r in LAT)
                               + "\nCr S\n1 2\nDirect\n0 0 0.5\n0.3333333 0.6666667 0.57\n0.3333333 0.6666667 0.43\n")
    (s3 / "INCAR").write_text("SYSTEM = CrS2 uniform\nENCUT = 520\nEDIFF = 1E-6\nALGO = Fast\nISYM = 2\n"
                              "LVHAR = .TRUE.\nKPAR = 4\nNCORE = 6\n")
    (s3 / "KPOINTS").write_text("auto\n0\nGamma\n 48 48 3\n")
    (s3 / "POTCAR").write_text("PAW_PBE Cr\nPAW_PBE S\n")
    (s3 / "CHGCAR").write_text(chgcar)                          # "" = 出厂 LCHARG = .FALSE. 的 0 字节占位
    if wavecar is not None:
        (s3 / "WAVECAR").write_text(wavecar)
    (s3 / "OUTCAR").write_text("  free  energy   TOTEN  =       -19.00000000 eV\n"
                               "  free  energy   TOTEN  =       %.8f eV\n General timing\n" % S3_TOTEN)
    (s3 / kc.METHOD_FILE).write_text("DIM=%s\n" % dim)
    kpts = [[0, 0, 0], list(K), [0.5, 0, 0]]
    eig = [[-1.0, 2.0], list(_model(K)), [-0.8, 1.6]]                  # 模型只在 K 附近成立：Γ、M 直接给远离带边的值
    (s3 / "vasprun.xml").write_text(_vasprun(kpts, eig, [[1.0, 0.0]] * 3))
    (m / "step8.2_dpt").mkdir()
    return m


def _dpt(m, me, mh):
    (m / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": "electron", "inputs": {"m_eff_m0": me}}, {"carrier": "hole", "inputs": {"m_eff_m0": mh}}]}))


def _eigenval(out, vasp6=True, noise_h=0.0, seed=0, **model):
    rng = np.random.default_rng(seed)
    out = out / "stage2_nscf" if (out / "stage2_nscf").is_dir() else out
    ln = (out / "KPOINTS").read_text().splitlines()
    kf = [[float(x) for x in l.split()[:3]] for l in ln[3:3 + int(ln[1])]]
    rows = ["    3    3    1    1", "  0.1E+02", "  1.0E-04", "  CAR", " CrS2", "     18  %d  2" % len(kf)]
    for k in kf:
        ev, ec = _model(k, **model)
        ev += noise_h * rng.standard_normal()                     # 价带加噪声（模拟非自洽本征值没收敛透）
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
        self.assertIn("-> electron 判定：不通过（和 S8.2 差 ≥ 3%）", txt)

    def test_fit_wrong_run(self):
        m = _material()
        out = m.parent / "emc"
        _run(["make", str(m), str(out)])
        (out / "EIGENVAL").write_text("    3    3    1    1\n x\n x\n x\n x\n 18 1 2\n\n 0 0 0 1\n 1 -1.0\n 2 1.0\n")
        rc, txt = _run(["fit", str(out)])
        self.assertEqual(rc, 2)
        self.assertIn("不是这次 make 的输入", txt)


class VerdictTests(unittest.TestCase):                                   # [V178]
    def _fit(self, **kw):
        m = _material()
        _dpt(m, ME, MH)
        out = m.parent / "emc"
        _run(["make", str(m), str(out)])
        _eigenval(out, **kw)
        rc, txt = _run(["fit", str(out)])
        return rc, txt, json.loads((out / E.RESULT).read_text())["carriers"]

    def test_clean_pass_is_explicit(self):
        rc, txt, c = self._fit()
        self.assertEqual(rc, 0, txt)
        self.assertIn("-> electron 判定：通过", txt)
        self.assertIn("-> hole 判定：通过", txt)
        self.assertIn("总判定：通过（退出码 0）", txt)
        self.assertLess(max(d["rms_meV"] for d in c["hole"]["directions"]), 1e-3)

    def test_unstable_half_window_is_a_fail_not_a_warn(self):
        # CrS₂ 空穴那种：窗口减半 > 2%。V177 只打 WARN，agent 当成了通过
        rc, txt, c = self._fit(noise_h=2e-5)
        self.assertEqual(rc, 1)
        self.assertIn("-> hole 判定：不通过（拟合不稳", txt)
        self.assertIn("-> electron 判定：通过", txt)
        self.assertIn("总判定：不通过 —— hole：拟合不稳", txt)
        self.assertEqual(c["hole"]["half_window_diagnosis"], "noise")
        self.assertIn("多半是小 q 点的能量噪声", txt)
        lo, hi = c["hole"]["m_range"]
        self.assertLess(lo, hi)
        self.assertIn("带边 m* 在 %.4f–%.4f 之间" % (lo, hi), txt)

    def test_higher_order_diagnosed(self):
        rc, txt, c = self._fit(c6h=6e4)
        self.assertEqual(rc, 1)
        self.assertEqual(c["hole"]["half_window_diagnosis"], "higher_order")
        self.assertIn("多半是带的高阶非抛物项，半窗口值更接近带边", txt)
        self.assertLess(c["hole"]["summary_half_window"]["anisotropy"], E.NOISE_ANISO)


class TwoStageTests(unittest.TestCase):                                  # [V177]
    def test_empty_chgcar_gives_two_stage(self):
        m = _material(chgcar="")
        out = m.parent / "emc"
        rc, txt = _run(["make", str(m), str(out)])
        self.assertEqual(rc, 0, txt)
        self.assertIn("两段式", txt)
        s1, s2 = out / "stage1_scf", out / "stage2_nscf"
        i1 = kc.parse_incar((s1 / "INCAR").read_text())
        self.assertEqual((i1["ISTART"], i1["ICHARG"], i1["LCHARG"], i1["LWAVE"], i1["NBANDS"], i1["ISYM"]),
                         ("1", "0", ".TRUE.", ".FALSE.", "2", "2"))   # ISYM 照 S3：WAVECAR 是 S3 的不可约 k 点
        self.assertNotIn("LVHAR", i1)
        self.assertEqual((s1 / "KPOINTS").read_text(), (m / "step3_uniform" / "KPOINTS").read_text())
        self.assertNotEqual(os.stat(s1 / "WAVECAR").st_ino, os.stat(m / "step3_uniform" / "WAVECAR").st_ino)
        self.assertEqual(kc.parse_incar((s2 / "INCAR").read_text())["ICHARG"], "11")
        self.assertFalse((s2 / "CHGCAR").exists())                    # 等 stage1 写出来
        self.assertEqual(int((s2 / "KPOINTS").read_text().splitlines()[1]), 121)
        meta = json.loads((out / E.META).read_text())
        self.assertEqual((meta["layout"], meta["s3_toten"]), ("two_stage", S3_TOTEN))

    def test_no_wavecar_starts_from_scratch(self):
        m = _material(chgcar="", wavecar=None)
        out = m.parent / "emc"
        rc, txt = _run(["make", str(m), str(out)])
        self.assertIn("成本约等于重跑一次 S3", txt)
        i1 = kc.parse_incar((out / "stage1_scf" / "INCAR").read_text())
        self.assertEqual((i1["ISTART"], i1["ICHARG"]), ("0", "2"))
        self.assertFalse((out / "stage1_scf" / "WAVECAR").exists())

    def test_run_two_stage_script(self):
        m = _material(chgcar="")
        out = m.parent / "emc"
        _run(["make", str(m), str(out)])
        fake = m.parent / "fakevasp.sh"                               # 假 VASP：写 OUTCAR；LCHARG = .TRUE. 时写 CHGCAR
        fake.write_text("#!/bin/bash\necho \"$PWD\" >> %s/calls\n"
                        "echo ' General timing' > OUTCAR\n"
                        "grep -q 'LCHARG *= *.TRUE.' INCAR && echo chg > CHGCAR\nexit 0\n" % m.parent)
        fake.chmod(0o755)
        r = subprocess.run(["bash", str(out / "run_two_stage.sh"), str(fake)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        calls = (m.parent / "calls").read_text().split()
        self.assertEqual([Path(c).name for c in calls], ["stage1_scf", "stage2_nscf"])
        self.assertEqual((out / "stage2_nscf" / "CHGCAR").read_text(), "chg\n")
        # stage1 没写出 CHGCAR -> 停下，不跑 stage2
        m2 = _material(chgcar="")
        out2 = m2.parent / "emc"
        _run(["make", str(m2), str(out2)])
        bad = m2.parent / "badvasp.sh"
        bad.write_text("#!/bin/bash\necho ' General timing' > OUTCAR\n")
        bad.chmod(0o755)
        r = subprocess.run(["bash", str(out2 / "run_two_stage.sh"), str(bad)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 1)
        self.assertIn("没写出 CHGCAR", r.stderr)
        self.assertFalse((out2 / "stage2_nscf" / "OUTCAR").exists())

    def test_fit_checks_stage1_energy(self):
        for toten, want in ((S3_TOTEN + 2e-5, 0), (S3_TOTEN + 2e-3, 1)):
            m = _material(chgcar="")
            _dpt(m, ME, MH)
            out = m.parent / "emc"
            _run(["make", str(m), str(out)])
            (out / "stage1_scf" / "OUTCAR").write_text("  free  energy   TOTEN  =       %.8f eV\n" % toten)
            _eigenval(out)
            rc, txt = _run(["fit", str(out)])
            self.assertEqual(rc, want, txt)
            self.assertIn("复现了 S3" if want == 0 else "电荷密度不是 S3 那一份", txt)


class ScanTests(unittest.TestCase):                                      # [V184]
    QS = [0.01 * (i + 1) for i in range(30)]

    def _pts(self):
        return E.star_points(REC, K, self.QS, 6)

    def test_isotropic_parabola_same_everywhere(self):
        pts = self._pts()
        e = np.array([EG + E.HB2M * np.linalg.norm(((p[0] - K) @ REC)[:2]) ** 2 / ME for p in pts])
        for w in (0.05, 0.3):
            for j in range(6):
                for model in ("fixed", "free"):
                    m, why = E.parabola_mass(e, pts, j, w, "electron", model)
                    self.assertAlmostEqual(m / ME, 1.0, delta=1e-9, msg=(w, j, model, why))

    def test_sixfold_quartic_gives_window_dependent_anisotropy(self):
        """q⁴cos6θ 是六方 K 谷里锯齿 / 扶手椅差别的最低阶来源：带边处两者相同，窗口越大差得越多。"""
        d = 2.0
        pts = self._pts()

        def en(p):
            v = ((p[0] - K) @ REC)[:2]
            q, th = np.linalg.norm(v), math.atan2(v[1], v[0])
            return EG + E.HB2M * q * q / ME + d * q ** 4 * math.cos(6 * th)
        e = np.array([en(p) for p in pts])
        ratios = []
        for w in (0.05, 0.15, 0.3):
            qs = np.array([q for q in self.QS if q <= w + 1e-12])
            shift = d * np.sum(qs ** 6) / np.sum(qs ** 4)                 # 固定顶点的最小二乘：c = c0 ± d·Σq⁶/Σq⁴
            m0, _ = E.parabola_mass(e, pts, 0, w, "electron")              # θ = 0：cos6θ = +1
            m90, _ = E.parabola_mass(e, pts, 3, w, "electron")             # θ = 90°：cos6θ = −1
            self.assertAlmostEqual(m0, E.HB2M / (E.HB2M / ME + shift), places=9)
            self.assertAlmostEqual(m90, E.HB2M / (E.HB2M / ME - shift), places=9)
            ratios.append(m90 / m0)
        self.assertTrue(1.0 < ratios[0] < ratios[1] < ratios[2], ratios)

    def test_one_sided_picks_up_warping(self):
        pts = self._pts()
        e = np.array([_model(p[0])[1] for p in pts])                       # 导带：含三角翘曲 w q³cos3θ
        mb, _ = E.parabola_mass(e, pts, 0, 0.1, "electron", side="both")
        mp, _ = E.parabola_mass(e, pts, 0, 0.1, "electron", side="plus")
        mm, _ = E.parabola_mass(e, pts, 0, 0.1, "electron", side="minus")
        self.assertLess(mp, mb)
        self.assertGreater(mm, mb)                                          # 两侧的三次项符号相反
        self.assertIsNone(E.parabola_mass(-e, pts, 0, 0.1, "electron")[0])  # 曲率符号不对

    def test_hex_labels(self):
        th = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0]
        self.assertEqual(E.hex_direction_labels(LAT, th), ["锯齿", "扶手椅", "锯齿", "扶手椅", "锯齿", "扶手椅"])
        rot = np.array([[0, A, 0], [-A * math.sqrt(3) / 2, -A / 2, 0], [0, 0, C]])   # a1 沿 y
        self.assertEqual(E.hex_direction_labels(rot, th)[:4], ["扶手椅", "锯齿", "扶手椅", "锯齿"])
        rect = np.array([[5.26, 0, 0], [0, 3.04, 0], [0, 0, C]])
        self.assertIsNone(E.hex_direction_labels(rect, th))

    def test_cli_scan(self):
        m = _material()
        out = m.parent / "emc"
        rc, txt = _run(["make", str(m), str(out), "--qmax", "0.3", "--nq", "30"])
        self.assertEqual(rc, 0, txt)
        _eigenval(out)
        rc, txt = _run(["scan", str(out), "--windows", "0.05,0.3,0.5"])
        self.assertEqual(rc, 0, txt)
        self.assertIn("超过 make 时的 qmax", txt)
        self.assertIn("0°（锯齿）", txt)
        self.assertIn("90°（扶手椅）", txt)
        self.assertIn("扶手椅/锯齿", txt)
        res = json.loads((out / E.SCAN).read_text())
        self.assertEqual(res["windows"], [0.05, 0.3, 0.5])
        r0 = res["carriers"]["electron"]["rows"][0]
        self.assertAlmostEqual(r0["m"][0] / ME, 1.0, delta=0.02)            # 小窗口 ≈ 带边质量
        self.assertIsNotNone(r0["armchair_over_zigzag"])
        rc, _ = _run(["fit", str(out)])                                     # fit 照旧能读同一次结果
        self.assertIn(rc, (0, 1))


if __name__ == "__main__":
    unittest.main(verbosity=2)
