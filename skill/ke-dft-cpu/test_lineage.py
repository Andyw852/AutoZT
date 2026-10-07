#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_lineage.py —— V140：rerun 之后的上游同源核对与下游失效。不需要 VASP / AMSET。

用法：python test_lineage.py      退出码 0 = 全部 PASS。

起因（Mo2S3，S1 重新弛豫后整条 rerun）：autozt rerun 先 rm -rf 步骤目录再 gen，V122/V136 的
"gen 见到旧产物 -> 下游失效"一个都不触发；S2.3 HSE 从还没重跑的 S2.2 拷了旧结构。
1) read_poscar_cell / cell_deviation：VASP5 Direct、Cartesian + scale、负 scale、Selective dynamics；
   CONTCAR 与 POSCAR 的文本精度差不算不同，0.3% 的晶格差算不同。
2) structure_lineage：结构不同的步骤（含 S2.3 的扇出子目录）、派生产物旧于来源（S4 h5、S7.1）都报；
   check_lineage 开着就退出、关着只告警。
3) 失效：S3 目录是新建的（rerun）-> S4 的 h5 与 S8 归档；S7 -> S7.1 三个产物 + S8；各 gen 无条件调用。
4) S8/S8.4 的 gen 接线；tools/lineage_check.py 的报告与 --invalidate-from。
5) V149：比对扣除整体平移（S3/S3b 的 align_origin 是刚性平移，WS2 被误报 4.45 Å）；平移之外的位移照报。
"""
import contextlib
import io
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "tools")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402

LAT = [[3.19, 0.0, 0.0], [-1.595, 2.762621, 0.0], [0.0, 0.0, 20.0]]
FRAC = [[0.0, 0.0, 0.5], [1 / 3, 2 / 3, 0.578], [1 / 3, 2 / 3, 0.422]]


def poscar(lat=LAT, frac=FRAC, scale=1.0, cart=False, sd=False, digits=8):
    import numpy as np
    lat = np.asarray(lat, float)
    xyz = np.dot(np.asarray(frac, float), lat) / scale if cart else np.asarray(frac, float)
    fmt = "%%.%df" % digits
    lines = ["MoS2", "%.10f" % scale]
    lines += ["  " + " ".join(fmt % x for x in row / scale) for row in lat]
    lines += ["Mo S", "1 2"]
    if sd:
        lines.append("Selective dynamics")
    lines.append("Cartesian" if cart else "Direct")
    lines += ["  " + " ".join(fmt % x for x in r) + (" T T T" if sd else "") for r in xyz]
    return "\n".join(lines) + "\n"


def _w(p, text, mtime=None):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    if mtime is not None:
        os.utime(p, (mtime, mtime))
    return p


class CellTests(unittest.TestCase):
    def test_formats_equal(self):
        d = Path(tempfile.mkdtemp())
        ref = kc.read_poscar_cell(_w(d / "a", poscar(digits=16)))
        for text in (poscar(digits=8), poscar(cart=True), poscar(scale=1.5), poscar(sd=True),
                     poscar(cart=True, scale=2.0)):
            got = kc.read_poscar_cell(_w(d / "b", text))
            self.assertLess(kc.cell_deviation(ref, got), 1e-6, text)
        # 负 scale = 体积
        import numpy as np
        vol = abs(np.linalg.det(np.asarray(LAT)))
        t = poscar().splitlines()
        t[1] = "-%.10f" % vol
        self.assertLess(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "c", "\n".join(t) + "\n"))), 1e-6)

    def test_differences(self):
        d = Path(tempfile.mkdtemp())
        ref = kc.read_poscar_cell(_w(d / "a", poscar()))
        lat2 = [[x * 1.003 for x in r] if i < 2 else r for i, r in enumerate(LAT)]
        self.assertGreater(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "b", poscar(lat=lat2)))), 5e-3)
        frac2 = [FRAC[0], FRAC[1], [1 / 3, 2 / 3, 0.4225]]
        self.assertGreater(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "c", poscar(frac=frac2)))), 5e-3)
        # 跨周期边界的同一位置不算不同
        frac3 = [[1.0, 0.0, 0.5], FRAC[1], FRAC[2]]
        self.assertLess(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "e", poscar(frac=frac3)))), 1e-6)
        bad = poscar().replace("Mo S\n1 2", "Mo S\n2 1")
        self.assertEqual(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "f", bad))), float("inf"))

    def test_rigid_translation_v149(self):
        """S3/S3b 的 align_origin 是刚性平移：WS2 三个原子同移 [0, 0, 4.4494] Å，不算不同；平移之外的位移照报。"""
        import numpy as np
        d = Path(tempfile.mkdtemp())
        ref = kc.read_poscar_cell(_w(d / "a", poscar()))
        for t in ([0.0, 0.0, 4.4494 / 20.0], [0.37, -0.81, 0.6], [0.5, 0.5, 0.5]):
            fr = (np.asarray(FRAC) + t) % 1.0
            self.assertLess(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "b", poscar(frac=fr)))), 1e-6, t)
        fr = (np.asarray(FRAC) + [0.0, 0.0, 0.2225]) % 1.0
        fr[2, 2] += 0.0005                                                   # 平移 + 一个 S 动 0.01 Å
        self.assertGreater(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "c", poscar(frac=fr)))), 5e-3)
        fr = np.asarray(FRAC) + [0.0, 0.0, 0.2225]
        lat2 = [[x * 1.003 for x in r] if i < 2 else r for i, r in enumerate(LAT)]
        self.assertGreater(kc.cell_deviation(ref, kc.read_poscar_cell(_w(d / "e", poscar(lat=lat2, frac=fr)))),
                           5e-3)

    def test_align_origin_output_is_same_structure_v149(self):
        """端到端：WS2 式 slab（原点不在高对称位置）经 align_origin 改写后，lineage 不报。"""
        try:
            from pymatgen.core import Lattice, Structure
        except ImportError:
            self.skipTest("需要 pymatgen + spglib")
        import shutil
        m = Path(tempfile.mkdtemp())
        st = Structure(Lattice.hexagonal(3.18, 20.0), ["W", "S", "S"],
                       [[0.11, 0.07, 0.13], [0.11 + 1 / 3, 0.07 + 2 / 3, 0.208], [0.11 + 1 / 3, 0.07 + 2 / 3, 0.052]])
        (m / "step1_opt").mkdir()
        (m / "step3_uniform").mkdir()
        st.to(filename=str(m / "step1_opt" / "CONTCAR"), fmt="poscar")
        shutil.copy(m / "step1_opt" / "CONTCAR", m / "step3_uniform" / "POSCAR")
        with contextlib.redirect_stdout(io.StringIO()):
            ao = kc.align_origin(m / "step3_uniform" / "POSCAR")
        self.assertIsNotNone(ao, "align_origin 应当平移（原点不在高对称位置）")
        a = kc.read_poscar_cell(m / "step1_opt" / "CONTCAR")
        b = kc.read_poscar_cell(m / "step3_uniform" / "POSCAR")
        import numpy as np
        self.assertGreater(float(np.abs(b[2] - a[2]).max()), 1e-3)          # 坐标确实变了
        self.assertLess(kc.cell_deviation(a, b), kc.LINEAGE_TOL_A)
        self.assertEqual(kc.structure_lineage(m)[1], [])


def _material(new_s1=True):
    """S1 重新弛豫过（CONTCAR = 新结构），下游有的跟上了、有的还是旧的。"""
    m = Path(tempfile.mkdtemp())
    now = time.time()
    old = poscar()
    lat2 = [[x * 1.004 for x in r] if i < 2 else r for i, r in enumerate(LAT)]
    new = poscar(lat=lat2, digits=16) if new_s1 else old
    _w(m / "step1_opt" / "CONTCAR", new)
    _w(m / "step3_uniform" / "POSCAR", new)
    _w(m / "step5_dielect" / "POSCAR", new)
    _w(m / "step6_elastic" / "POSCAR", new)
    _w(m / "step7_deform" / "undeformed" / "POSCAR", new)
    _w(m / "step2_bandgap" / "step2.1_static" / "POSCAR", new)
    _w(m / "step2_bandgap" / "step2.2_pbe" / "POSCAR", old)                 # 还没重跑
    _w(m / "step2_bandgap" / "step2.3_hse" / "p1of2" / "POSCAR", old)       # 从旧 S2.2 拷的
    _w(m / "step2_bandgap" / "step2.3_hse" / "p2of2" / "POSCAR", old)
    _w(m / "step3_uniform" / "vasprun.xml", "x", mtime=now)
    _w(m / "step4_wave" / "wavefunction.h5", "x", mtime=now - 7200)          # 用旧 S3 做的
    _w(m / "step7_deform" / "deform-001" / "vasprun.xml", "x", mtime=now - 7200)
    _w(m / "step7b_deform_read" / "deformation.h5", "x", mtime=now - 3600)   # 比 S7 新：没问题
    return m


class LineageTests(unittest.TestCase):
    def test_report(self):
        m = _material()
        ref, probs = kc.structure_lineage(m)
        self.assertEqual(ref, m / "step1_opt" / "CONTCAR")
        steps = [s for s, _ in probs]
        self.assertEqual(sorted(steps), sorted(["step2_bandgap/step2.2_pbe", "step2_bandgap/step2.3_hse",
                                                "step4_wave"]), probs)
        self.assertEqual(steps.count("step2_bandgap/step2.3_hse"), 1)           # 扇出目录只报一次
        why = dict(probs)
        self.assertIn("p1of2/POSCAR", why["step2_bandgap/step2.3_hse"])
        self.assertIn("step3_uniform/vasprun.xml", why["step4_wave"])

    def test_clean_and_no_s1(self):
        m = _material(new_s1=False)
        for p in ("step4_wave/wavefunction.h5",):
            os.utime(m / p, None)                                               # h5 比 S3 新
        self.assertEqual(kc.structure_lineage(m)[1], [])
        (m / "step1_opt" / "CONTCAR").unlink()
        ref, probs = kc.structure_lineage(m)
        self.assertIsNone(ref)
        self.assertEqual(probs, [])

    def test_check_lineage(self):
        m = _material()
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            with self.assertRaises(SystemExit) as cm:
                kc.check_lineage(m, enabled=True, label="S8.4")
        self.assertIn("STRUCTURE_GUARD = false", str(cm.exception))
        self.assertIn("★ 上游同源核对", buf.getvalue())
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            probs = kc.check_lineage(m, enabled=False)
        self.assertEqual(len(probs), 3)
        self.assertIn("[WARN]", buf.getvalue())


class InvalidateTests(unittest.TestCase):
    def test_s3_rerun_fresh_dir(self):
        m = Path(tempfile.mkdtemp())
        s3 = m / "step3_uniform"
        s3.mkdir()
        _w(m / "step4_wave" / "wavefunction.h5", "x")
        (m / "step4_wave" / "vasprun.xml").symlink_to("../step3_uniform/vasprun.xml")   # rerun 删掉后是断链
        _w(m / "step8.4_amset2d" / "transport.json", "{}")
        (m / "step8.4_amset2d" / "wavefunction.h5").symlink_to("../step4_wave/wavefunction.h5")
        snap = kc.snapshot_inputs(s3)                                       # rerun 后：全空
        _w(s3 / "POSCAR", poscar())
        _w(s3 / "INCAR", "ISTART = 0\n")
        with contextlib.redirect_stdout(io.StringIO()):
            changed, n, down = kc.finalize_stale_inputs(m, s3, "step3_uniform", snap)
        self.assertEqual((changed, n), ([], 0))
        self.assertIn(("step4_wave", "wavefunction.h5"), down)
        self.assertIn(("step8.4_amset2d", "transport.json"), down)
        # 原地重新 gen、输入没变：不失效（原行为）
        _w(m / "step4_wave" / "wavefunction.h5", "x")
        snap = kc.snapshot_inputs(s3)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(kc.finalize_stale_inputs(m, s3, "step3_uniform", snap), ([], 0, []))

    def test_s7_to_s71(self):
        m = Path(tempfile.mkdtemp())
        (m / "step7_deform").mkdir()
        for f in ("deformation.h5", "deformation_vac.h5", "band_edges.json"):
            _w(m / "step7b_deform_read" / f, "x")
        _w(m / "step8.4_amset2d" / "transport.json", "{}")
        (m / "step8.4_amset2d" / "deformation.h5").symlink_to("../step7b_deform_read/deformation_vac.h5")
        _w(m / "step8.2_dpt" / "dpt_result.json", "{}")
        with contextlib.redirect_stdout(io.StringIO()):
            down = kc.invalidate_downstream(m, "step7_deform", "test")
        self.assertEqual(sorted(f for s, f in down if s == "step7b_deform_read"),
                         ["band_edges.json", "deformation.h5", "deformation_vac.h5"])
        self.assertIn(("step8.4_amset2d", "transport.json"), down)
        self.assertIn(("step8.2_dpt", "dpt_result.json"), down)
        self.assertFalse(any((m / "step7b_deform_read").glob("deformation*.h5")))

    def test_gens_unconditional(self):
        src = {f: (ROOT / f).read_text(encoding="utf-8") for f in (
            "step4_wave/gen_step6_wave.py", "step4b_wave_full/gen_step6b_wave_full.py",
            "step5_dielect/gen_step8_dielect.py", "step6_elastic/gen_step2_elastic.py",
            "step7_deform/step7b_read/gen_step9b_deform_read.py", "step7_deform/gen_step9_deform.py")}
        for f, s in src.items():
            self.assertIn("invalidate_downstream(", s, f)
            self.assertNotIn('and (out / "wavefunction.h5").is_file()', s, f)
            self.assertNotIn('if (out / "OUTCAR").is_file():\n        kc.invalidate_downstream', s, f)
            self.assertNotIn('if (step2 / "OUTCAR").is_file():', s, f)
        s7 = src["step7_deform/gen_step9_deform.py"]
        self.assertLess(s7.index('kc.invalidate_downstream(cwd, OUTDIR_NAME'), s7.index("kc.relay_poscar("))


class WiringTests(unittest.TestCase):
    def test_gens(self):
        for f, label in (("step8_amset/gen_step10_amset.py", "S8"),
                         ("step8.4_amset2d/gen_step14_amset2d.py", "S8.4")):
            s = (ROOT / f).read_text(encoding="utf-8")
            self.assertIn('"STRUCTURE_GUARD": (True, "bool")', s)
            self.assertIn('STRUCTURE_GUARD = bool(_p["STRUCTURE_GUARD"])', s)
            i = s.index('kc.check_lineage(cwd, enabled=STRUCTURE_GUARD, label="%s")' % label)
            self.assertLess(s.index("except (KeyError, ValueError, TypeError) as _e:"), i)
            self.assertLess(i, s.index("_warn_small_pbe_gap(cwd)\n    elastic"))

    def test_cli(self):
        import lineage_check as L
        root = Path(tempfile.mkdtemp())
        m = _material()
        target = root / "Mo2S3" / "ke-dft-cpu"
        target.parent.mkdir()
        m.rename(target)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertEqual(L.main([str(root)]), 1)
        self.assertIn("[★]", buf.getvalue())
        self.assertIn("step2_bandgap/step2.3_hse", buf.getvalue())
        (target / "step4_wave" / "vasprun.xml").symlink_to("../step3_uniform/vasprun.xml")   # S4 软链 S3
        _w(target / "step8.4_amset2d" / "transport.json", "{}")
        (target / "step8.4_amset2d" / "wavefunction.h5").symlink_to("../step4_wave/wavefunction.h5")
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            L.main([str(target), "--invalidate-from", "step3_uniform"])
        self.assertIn("step4_wave/wavefunction.h5", buf.getvalue())
        self.assertFalse((target / "step8.4_amset2d" / "transport.json").exists())
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stdout(io.StringIO()):
                L.main([str(target), "--invalidate-from", "step99"])


class FindMaterialsV183Tests(unittest.TestCase):
    """V183：project_setting/templates/ 里的步骤目录模板不当成材料。"""

    def test_skips_templates_and_archives(self):
        import lineage_check as L
        root = Path(tempfile.mkdtemp())
        for rel in ("CrA/ke-dft-cpu/step1_opt", "CrB/ke-dft-cpu/step1_opt",
                    "project_setting/templates/ke-dft-cpu/step1_opt", "project_setting/templates/step1_opt",
                    "CrC.stale-20261001/ke-dft-cpu/step1_opt", ".trash/CrD/step1_opt"):
            (root / rel).mkdir(parents=True)
        got = [p.relative_to(root).as_posix() for p in L.find_materials(root)]
        self.assertEqual(got, ["CrA/ke-dft-cpu", "CrB/ke-dft-cpu"])
        # 直接指向一个材料目录时照旧
        self.assertEqual(L.find_materials(root / "CrA" / "ke-dft-cpu"), [root / "CrA" / "ke-dft-cpu"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
