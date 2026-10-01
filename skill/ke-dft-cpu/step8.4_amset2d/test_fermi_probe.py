#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_fermi_probe.py —— V131：tools/fermi_probe.py 的核心函数（合成 DOS，不需要 VASP）。

用法：python test_fermi_probe.py      退出码 0 = 全部 PASS。

1) 干净的二维阶梯 DOS（带隙 1.9 eV，WS2 报错里的浓度 −2.759e-7 bohr⁻³，100–1000 K）：
   AMSET 的 get_fermi 能解，二分法与它一致到 1e-4 eV，实现浓度误差 < 1e-6。
2) 本征费米能级为 NaN：AMSET 的 get_fermi 不报错、静默返回 NaN（复现 WS2 IR_FIX 开时 pinv 崩的路径），
   probe 把它标成 nan；二分法不依赖本征费米能级，照样解出。
3) 带隙中间有插值假态：dos_diagnostics 数得出带隙内的态数。
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
            ef = fd.get_fermi(conc, 300.0, tol=1e-5, precision=10)       # AMSET：不报错
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


class PluginNanGuardTests(unittest.TestCase):
    """加载 amset2d_plugin 后，get_fermi 得到 NaN 会抛 ValueError（AMSET 原式静默返回 NaN）。"""

    def test_nan_raises_after_plugin(self):
        import json
        import os
        import subprocess
        import tempfile
        import importlib.util
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
        e = dict(os.environ)
        paths = [str(HERE), str(HERE.parent / "tools")]
        if importlib.util.find_spec("sumo") is None:
            d = Path(tempfile.mkdtemp()) / "sumo"
            d.mkdir()
            (d / "__init__.py").write_text("")
            (d / "symmetry.py").write_text("class Kpath: pass\nclass PymatgenKpath(Kpath): pass\n")
            paths.append(str(d.parent))
        e["PYTHONPATH"] = os.pathsep.join(paths + [e.get("PYTHONPATH", "")])
        e["AMSET2D_RECORD"] = str(HERE / "example_2d_correction.json")
        e["AZ_IR_FIX"] = "0"
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=e,
                           cwd=tempfile.mkdtemp(), timeout=600)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        self.assertEqual(json.loads(r.stdout.strip().splitlines()[-1]), {"ok_finite": True, "raised": True})


if __name__ == "__main__":
    unittest.main(verbosity=2)
