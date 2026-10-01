#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_fermi_probe.py —— V131：tools/fermi_probe.py 的核心函数（合成 DOS，不需要 VASP）。

用法：python test_fermi_probe.py      退出码 0 = 全部 PASS。

1) 干净的二维阶梯 DOS（带隙 1.9 eV，WS2 报错里的浓度 −2.759e-7 bohr⁻³，100–1000 K）：
   AMSET 的 get_fermi 能解，二分法与它一致到 1e-4 eV，实现浓度误差 < 1e-6。
2) 本征费米能级为 NaN：AMSET 的 get_fermi 不报错、静默返回 NaN（复现 WS2 IR_FIX 开时 pinv 崩的路径），
   probe 把它标成 nan；二分法不依赖本征费米能级，照样解出。
3) 带隙中间有插值假态：dos_diagnostics 数得出带隙内的态数。
4) V133：3D 运行目录（band_structure_data.json）上整条 fermi_probe 能跑；AMSET 原式在 tol=1 静默接受
   偏差 ~100% 的解（含 S8 出厂网格的 −1e17 cm⁻³ @ 100 K），probe 标"宽松!"；挂上 amset_fermi_fix
   后这些点都走二分、carrier_guard 通过。
"""
import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(HERE.parent / "tools")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
    if not hasattr(np, _a):
        setattr(np, _a, getattr(np, _b))

try:
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Spin
    from amset.electronic_structure.dos import FermiDos
    from amset.constants import ev_to_hartree as EH
    import fermi_probe as FP
    _HAS = True
except Exception:                                  # noqa: BLE001
    _HAS = False


def _dos(gap=1.9, in_gap=0.0, efermi=None):
    st = Structure(Lattice.hexagonal(3.18, 20.0), ["W", "S", "S"],
                   [[0, 0, .5], [1 / 3, 2 / 3, .578], [1 / 3, 2 / 3, .422]])
    E = np.linspace(-4 * EH, 6 * EH, 1000)
    d = np.where(E <= 0, 30.0, 0.0) + np.where(E >= gap * EH, 25.0, 0.0)
    if in_gap:
        d[np.abs(E - 0.5 * gap * EH) < 0.02 * EH] += in_gap
    fd = FermiDos(0.5 * gap * EH if efermi is None else efermi, E, {Spin.up: d}, st, atomic_units=True,
                  dos_weight=2)
    return fd, st


@unittest.skipUnless(_HAS, "需要 amset + pymatgen")
class FermiProbeTests(unittest.TestCase):
    def test_clean_dos_agrees(self):
        fd, st = _dos()
        rep = FP.probe_dos(fd, [-2.759e-7 * (1 / 0.529177210903e-8) ** 3], [100, 300, 1000], 0.0, 1.9 * EH,
                           fd.structure.volume)
        for r in rep["rows"]:
            self.assertTrue(r["amset_ok"])
            self.assertLess(abs(r["amset_ef_eV"] - r["bisect_ef_eV"]), 1e-4)
            self.assertLess(r["bisect_rel_err"], 1e-6)
        self.assertAlmostEqual(rep["dos"]["gap_bands_eV"], 1.9, places=6)
        self.assertLess(rep["dos"]["states_in_gap_per_cell"], 1e-12)

    def test_nan_intrinsic_fermi_is_silent_in_amset(self):
        fd, st = _dos()
        fd.efermi = float("nan")
        conc = -2.759e-7
        with np.errstate(all="ignore"):
            ef = FP._original_get_fermi(fd)(fd, conc, 300.0, tol=1e-5, precision=10)   # AMSET 原式：不报错
            self.assertTrue(np.isnan(ef))
            a = FP.amset_ladder(fd, conc, 300.0)
            self.assertTrue(a["nan"])
            self.assertFalse(a["ok"])
            b = FP.bisect_fermi(fd, conc, 300.0)
        self.assertIsNotNone(b)
        self.assertLess(abs(fd.get_doping(b, 300.0) / conc - 1), 1e-6)

    def test_in_gap_states_counted(self):
        fd, st = _dos(in_gap=0.5)
        diag = FP.dos_diagnostics(fd, 0.0, 1.9 * EH, fd.structure.volume)
        self.assertGreater(diag["states_in_gap_per_cell"], 1e-4)


def _sub_env():
    """子进程环境：本目录 + tools 在 PYTHONPATH 上；没装 sumo 时给一个最小桩（amset 顶层会 import 它）。"""
    import os
    import tempfile
    import importlib.util
    e = dict(os.environ)
    paths = [str(HERE), str(HERE.parent / "tools")]
    if importlib.util.find_spec("sumo") is None:
        d = Path(tempfile.mkdtemp()) / "sumo"
        d.mkdir()
        (d / "__init__.py").write_text("")
        (d / "symmetry.py").write_text("class Kpath: pass\nclass PymatgenKpath(Kpath): pass\n")
        paths.append(str(d.parent))
    e["PYTHONPATH"] = os.pathsep.join(paths + [e.get("PYTHONPATH", "")])
    e.pop("AZ_FERMI_FIX", None)
    e["AZ_IR_FIX"] = "0"
    return e


class PluginNanGuardTests(unittest.TestCase):
    """加载 amset2d_plugin 后，get_fermi 得到 NaN 会抛 ValueError（AMSET 原式静默返回 NaN）。"""

    def test_nan_raises_after_plugin(self):
        import json
        import subprocess
        import tempfile
        code = '''
import json, numpy as np
for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
    if not hasattr(np, _a):
        setattr(np, _a, getattr(np, _b))
import amset2d_plugin
import test_fermi_probe as T
fd, st = T._dos()
ok = float(fd.get_fermi(-2.759e-7, 300.0, tol=1e-5, precision=10))
fd.efermi = float("nan")
try:
    with np.errstate(all="ignore"):
        fd.get_fermi(-2.759e-7, 300.0, tol=1e-5, precision=10)
    raised = False
except ValueError as e:
    raised = "NaN" in str(e)
print(json.dumps({"ok_finite": bool(np.isfinite(ok)), "raised": raised}))
'''
        e = _sub_env()
        e["AMSET2D_RECORD"] = str(HERE / "example_2d_correction.json")
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=e,
                           cwd=tempfile.mkdtemp(), timeout=600)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        self.assertEqual(json.loads(r.stdout.strip().splitlines()[-1]), {"ok_finite": True, "raised": True})


# V133：3D 运行目录（band_structure_data.json + settings.yaml，没有 2d_correction.json）。纤锌矿型四带模型、
# 带隙 2.2 eV，掺杂/温度含 S8 出厂网格的最低点 −1e17 cm⁻³ @ 100 K。
PROBE_3D = '''
import json, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
    if not hasattr(np, _a):
        setattr(np, _a, getattr(np, _b))
import spglib, yaml
from monty.serialization import dumpfn
from pymatgen.core import Lattice, Structure
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.electronic_structure.core import Spin

st = Structure(Lattice.hexagonal(3.19, 5.19), ["Ga", "Ga", "N", "N"],
               [[1/3, 2/3, 0], [2/3, 1/3, .5], [1/3, 2/3, .377], [2/3, 1/3, .877]])
m = [12, 12, 8]
mp, grid = spglib.get_ir_reciprocal_mesh(m, (st.lattice.matrix, st.frac_coords, [s.Z for s in st.species]))
k = grid[np.unique(mp)] / np.array(m, float)
kc = k @ st.lattice.reciprocal_lattice.matrix
L = st.lattice.matrix
s = sum(np.cos(kc @ (np.array(r, float) @ L)) for r in ((1, 0, 0), (0, 1, 0), (1, 1, 0))) + 0.7 * np.cos(kc @ L[2])
E = np.vstack([-0.3 * (3.7 - s) - 2, -0.3 * (3.7 - s), 2.2 + 0.35 * (3.7 - s), 4.7 + 0.35 * (3.7 - s)])
bs = BandStructure(k, {Spin.up: E}, st.lattice.reciprocal_lattice, 1.1, structure=st)
d = bs.as_dict()
d["structure"] = st.as_dict()                  # as_dict 不带 structure（无投影时），AMSET 插值要它
dumpfn({"nelect": 4, "band_structure": d}, "band_structure_data.json")
yaml.safe_dump({"doping": [-1e15, -1e16, -1e17, -1e18, 1e16, 1e18], "temperatures": [50, 100, 300],
                "interpolation_factor": 5, "energy_cutoff": 2.0, "scattering_type": ["ADP"],
                "deformation_potential": [3.0, 3.0], "elastic_constant": (np.eye(6) * 200).tolist(),
                "nworkers": 1, "dos_estep": 0.01, "symprec": 0.01}, open("settings.yaml", "w"))

import fermi_probe
fermi_probe.main(["--dir", ".", "--json", "probe.json"])
rows = json.load(open("probe.json"))["rows"]
import postprocess_intrinsic as P, amset_fermi_fix as FF
guard = "ok"
try:
    P.build_amset_data(".", ir_fix=False)       # set_doping_and_temperatures（挂着插件）+ carrier_guard（>1% 退出）
except SystemExit as e:
    guard = str(e)
print(json.dumps({"loose": [[r["doping_cm3"], r["T"], r["amset_rel_err"], r["bisect_rel_err"]] for r in rows
                            if r["amset_ok"] and r["amset_rel_err"] > 0.01],
                  "failed": [[r["doping_cm3"], r["T"]] for r in rows if not r["amset_ok"]],
                  "fallbacks": len(FF.STATE["fallbacks"]), "guard": guard}))
'''


class ThreeDRunDirTests(unittest.TestCase):
    """V133：fermi_probe 在 3D 运行目录可用；AMSET 原式在 tol=1 静默接受偏差 100% 的解（含出厂网格的
    −1e17 @ 100 K）；挂上 amset_fermi_fix 后这些点都走二分，carrier_guard 通过。"""

    @unittest.skipUnless(_HAS, "需要 amset + pymatgen")
    def test_probe_flags_silent_loose_and_fix_solves(self):
        import json
        import subprocess
        import tempfile
        r = subprocess.run([sys.executable, "-c", PROBE_3D], capture_output=True, text=True, env=_sub_env(),
                           cwd=tempfile.mkdtemp(), timeout=1200)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        self.assertIn("宽松!", r.stdout)
        out = json.loads(r.stdout.strip().splitlines()[-1])
        loose = {(d, T) for d, T, _, _ in out["loose"]}
        self.assertIn((-1e17, 100.0), loose)
        for d, T, err, berr in out["loose"]:
            self.assertGreater(err, 0.5)                     # 原式接受的解差了一倍上下（载流子几乎为 0）
            self.assertLess(berr, 1e-6)                      # 二分解照样准确
        self.assertEqual(out["failed"], [])                  # 这些点原式不报错 —— 是静默的
        self.assertEqual(out["fallbacks"], len(out["loose"]))
        self.assertEqual(out["guard"], "ok")


if __name__ == "__main__":
    unittest.main(verbosity=2)
