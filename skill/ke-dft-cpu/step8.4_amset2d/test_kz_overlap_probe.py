#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_kz_overlap_probe.py —— V126：k_z 插值重叠假象的探针 + postprocess 的 kz_cap 适配（不需要 VASP）。

用法：python test_kz_overlap_probe.py      退出码 0 = 全部 PASS。

1) 紧束缚单层（MoS2 类，12×12×3 全网格）写成真 h5：
   · 每个 k 乘随机整体相位（VASP 式任意规范）-> 21 层集合的平均重叠明显下降、最小值掉到 ~0.5，
     预测 μ新/μ旧 < 0.97（与 MoS₂ 实测 n 型 −5.7% 同量级）；
   · 平滑规范 -> 预测差 < 1.5%；
   · 3 层集合全在 h5 网格层上 -> 结果与规范无关（两种规范相同到 h5 单精度）。
2) main() 端到端（能带窗口读 amset.log，能量走 vasprun 接口）。
3) postprocess_intrinsic._apply_run_plugins：2d_correction.json 开了 kz_cap_rmax 才加载插件；缺插件拒绝。
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT), str(ROOT / "tools"), str(ROOT.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import test_desym_fix as TDF
    import test_symmetry_gate as tsg
    from amset.wavefunction.io import write_coefficients
    from pymatgen.electronic_structure.core import Spin
    _HAS = TDF._HAS_AMSET
except Exception:                                              # noqa: BLE001
    _HAS = False

import kz_overlap_probe as P  # noqa: E402


@unittest.skipUnless(_HAS, "需要 amset + pymatgen")
class ProbeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(np, "string_"):
            np.string_ = np.bytes_
        cls.st = tsg.mos2_aligned()
        cls.kp = TDF.full_mesh([12, 12, 3])
        cls.gp = TDF.gsphere(cls.st)
        cls.cf, cls.en = TDF.tb_coeffs(cls.st, cls.kp, cls.gp)
        rng = np.random.default_rng(3)
        cls.dirs = {}
        for tag, ph in (("smooth", np.ones(len(cls.kp))), ("random", np.exp(2j * np.pi * rng.random(len(cls.kp))))):
            d = Path(tempfile.mkdtemp())
            write_coefficients({Spin.up: cls.cf * ph[None, :, None]}, cls.gp, cls.kp, cls.st,
                               filename=str(d / "wavefunction.h5"))
            nb = cls.cf.shape[0]
            (d / "amset.log").write_text("  Interpolating spin-up bands 1-%d\n" % nb)
            cls.dirs[tag] = d

    def _ratio(self, tag, band):
        calc = P.load_calculator(self.dirs[tag] / "wavefunction.h5", None, False)
        k0 = np.array([1 / 12, 2 / 12, 0.0])
        Mo = P.pair_matrix(calc, Spin.up, band, k0, P.kz_set(21))
        Mn = P.pair_matrix(calc, Spin.up, band, k0, P.kz_set(3))
        return Mo, Mn, P.summarize(Mn)["inv_tau_weight"] / P.summarize(Mo)["inv_tau_weight"]

    def test_random_gauge_artifact_and_grid_layers_gauge_free(self):
        worst_random, worst_smooth = 1.0, 1.0
        for b in range(self.cf.shape[0]):
            Mo_r, Mn_r, r_r = self._ratio("random", b)
            Mo_s, Mn_s, r_s = self._ratio("smooth", b)
            self.assertLess(np.abs(Mn_r - Mn_s).max(), 1e-5)           # 网格层上与规范无关（h5 单精度）
            worst_random, worst_smooth = min(worst_random, r_r), min(worst_smooth, r_s)
            self.assertLess(Mo_r.min(), Mo_s.min() + 1e-9)
        self.assertLess(worst_random, 0.97)
        self.assertGreater(worst_smooth, 0.985)

    def test_main_end_to_end(self):
        d = self.dirs["random"]
        en = np.sort(self.en, axis=0)
        gaps = [(en[b + 1].min() - en[b].max(), b) for b in range(len(en) - 1)]
        g, vb = max(gaps)
        self.assertGreater(g, 0)
        ef = float(en[vb].max())
        old = P._energies_from_vasprun
        P._energies_from_vasprun = lambda *_a, **_k: ({Spin.up: self.en}, self.kp, ef)
        try:
            rc = P.main(["--dir", str(d), "--no-desym-fix", "--emax", "0.5", "--nsample", "4",
                         "--json", str(d / "probe.json")])
        finally:
            P._energies_from_vasprun = old
        self.assertEqual(rc, 0)
        rep = json.loads((d / "probe.json").read_text(encoding="utf-8"))
        self.assertEqual(rep["band_window"], [1, self.cf.shape[0]])
        self.assertIn("electron", rep["probe"])
        self.assertIn("hole", rep["probe"])
        for car in ("electron", "hole"):
            r = rep["probe"][car]
            self.assertLessEqual(r["mean_overlap_old"], r["mean_overlap_new"] + 1e-9)
            for row in r["rows"]:
                self.assertAlmostEqual(row["k"][2], 0.0)

    def test_band_window_parse(self):
        self.assertEqual(P.band_window(self.dirs["smooth"]), (1, self.cf.shape[0]))
        self.assertEqual(P.band_window(self.dirs["smooth"], "11:17"), (11, 17))
        with self.assertRaises(ValueError):
            P.band_window(tempfile.mkdtemp())


class PostprocessKzCapTests(unittest.TestCase):
    CODE = textwrap.dedent('''
        import json, sys, numpy as np
        if not hasattr(np, "Inf"):
            np.Inf = np.inf
        if not hasattr(np, "trapz"):
            np.trapz = np.trapezoid
        import postprocess_intrinsic as PP
        try:
            loaded = PP._apply_run_plugins(sys.argv[1])
        except SystemExit as e:
            print(json.dumps({"exit": str(e)[:40]})); sys.exit(0)
        import amset.scattering.calculate as C
        import amset.interpolation.bandstructure as B
        print(json.dumps({"loaded": loaded, "kzcap": bool(getattr(C, "_amset2d_kzcap", False)),
                          "init": B.Interpolator.__init__.__name__}))
    ''')

    def _case(self, rmax, plugin=True):
        d = Path(tempfile.mkdtemp())
        rec = json.loads((HERE / "example_2d_correction.json").read_text())
        if rmax:
            rec["kz_cap_rmax"] = rmax
        (d / "2d_correction.json").write_text(json.dumps(rec))
        if plugin:
            shutil.copyfile(HERE / "amset2d_plugin.py", d / "amset2d_plugin.py")
        env = dict(os.environ)
        stub = None
        import importlib.util
        if importlib.util.find_spec("sumo") is None:
            stub = Path(tempfile.mkdtemp())
            (stub / "sumo").mkdir()
            (stub / "sumo" / "__init__.py").write_text("")
            (stub / "sumo" / "symmetry.py").write_text("class Kpath: pass\nclass PymatgenKpath(Kpath): pass\n")
        env["PYTHONPATH"] = os.pathsep.join(p for p in (str(HERE), str(stub) if stub else "",
                                                        env.get("PYTHONPATH", "")) if p)
        r = subprocess.run([sys.executable, "-c", self.CODE, str(d)], capture_output=True, text=True,
                           env=env, cwd=tempfile.mkdtemp(), timeout=600)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        return json.loads(r.stdout.strip().splitlines()[-1])

    def test_loads_plugin_when_capped(self):
        out = self._case(1)
        self.assertTrue(out["loaded"])
        self.assertTrue(out["kzcap"])
        self.assertEqual(out["init"], "_init_kz")

    def test_noop_when_not_capped(self):
        out = self._case(None)
        self.assertFalse(out["loaded"])
        self.assertFalse(out["kzcap"])

    def test_refuses_without_plugin(self):
        out = self._case(1, plugin=False)
        self.assertIn("exit", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
