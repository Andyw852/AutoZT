#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_fermi_fix.py —— V132：amset_fermi_fix（求费米能级的稳健化）。合成 DOS，不需要 VASP。

用法：python test_fermi_fix.py      退出码 0 = 全部 PASS。

1) WS2 同款带边（本征 E_F −2.1749 eV、VBM −3.243、CBM −1.107）、二维单层量级的导带 DOS（0.2 态/eV/胞）：
   AMSET 原搜索在 n 型 100 K 的整条容差阶梯上都失败（复现作业报错）；挂上插件后
   AmsetData.set_doping_and_temperatures 在第一级 1e-5 就拿到二分解，实现浓度误差 < 1e-8。
2) 原搜索能解的点：插件返回值与原式逐位相同（含电子/空穴浓度）。
3) NaN：原式静默返回 NaN，插件抛 ValueError。
4) 浓度超出 DOS 能容纳的范围：原式在 tol=1 接受偏差几十 % 的解，插件拒绝（> 1%）。
   旧运行目录里 V131 的 amset2d_plugin 先套了 NaN 壳时，插件照样挂上。
5) AZ_FERMI_FIX=0：回到原行为（只留 NaN 守卫）。
6) 幂等：重复 apply 不叠壳。
7) 接线：S8 / S8.4 的 gen、amset2d_plugin、postprocess、skill.yaml（ke 两处、zt 一处）、zt 的软链。
"""
import os
import sys
import types
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
KE = HERE.parent
SKILL = KE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
    if not hasattr(np, _a):
        setattr(np, _a, getattr(np, _b))

try:
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Spin
    from amset.electronic_structure.dos import FermiDos
    from amset.core.data import AmsetData
    from amset.constants import ev_to_hartree as EH, cm_to_bohr
    import amset_fermi_fix as FF
    _HAS = True
except Exception:                                  # noqa: BLE001
    _HAS = False

VBM, CBM, EF = -3.243, -1.107, -2.1749


def _ws2_dos(g_cb=0.2, estep=0.005, cb_width=3.0):
    """WS2 带边；价带 2 eV 宽装 4 个电子，导带 g_cb 态/eV/胞（含自旋）。c = 30 Å（≈ 3.5e12 cm⁻² ↔ 1.16e19 cm⁻³）。"""
    st = Structure(Lattice.hexagonal(3.18, 30.0), ["W", "S", "S"],
                   [[0, 0, .5], [1 / 3, 2 / 3, .553], [1 / 3, 2 / 3, .447]])
    E = np.arange(VBM - 2.0, CBM + cb_width, estep) * EH
    d = (np.where((E <= VBM * EH) & (E > (VBM - 2.0) * EH), 2.0 / EH, 0.0)
         + np.where(E >= CBM * EH, g_cb / EH, 0.0))
    return FermiDos(EF * EH, E, {Spin.up: d / 2}, st, atomic_units=True, dos_weight=2)


def _conc(cm3):
    return cm3 * (1 / cm_to_bohr) ** 3


def _ladder(fd, dop_cm3, temps):
    ns = types.SimpleNamespace(dos=fd, num_electrons=float(fd.nelect))
    with np.errstate(all="ignore"):
        AmsetData.set_doping_and_temperatures(ns, np.asarray(dop_cm3, float), np.asarray(temps, float))
    return ns


class _Env:
    def __init__(self, value):
        self.value = value

    def __enter__(self):
        self.old = os.environ.get(FF.ENV_NAME)
        if self.value is None:
            os.environ.pop(FF.ENV_NAME, None)
        else:
            os.environ[FF.ENV_NAME] = self.value

    def __exit__(self, *exc):
        if self.old is None:
            os.environ.pop(FF.ENV_NAME, None)
        else:
            os.environ[FF.ENV_NAME] = self.old


@unittest.skipUnless(_HAS, "需要 amset + pymatgen")
class FermiFixTests(unittest.TestCase):
    def setUp(self):
        FF.apply()
        self.orig = FF.original_get_fermi()

    def test_ws2_like_n_100K_original_fails_fix_solves(self):
        fd = _ws2_dos()
        for cm3 in (-1.86e18, -3.31e17):
            c = _conc(cm3)
            for tol in np.logspace(-5, 0, 6):                      # 前提：原搜索整条阶梯都失败
                with np.errstate(all="ignore"), self.assertRaises(ValueError):
                    self.orig(fd, c, 100.0, tol=tol, precision=10)
        n0 = len(FF.STATE["fallbacks"])
        with _Env(None):
            ns = _ladder(fd, [-1.86e18, -3.31e17], [100.0, 300.0])
        self.assertEqual(len(FF.STATE["fallbacks"]) - n0, 2)        # 只有两个 100 K 点走了二分
        for i, cm3 in enumerate((-1.86e18, -3.31e17)):
            for j, T in enumerate((100.0, 300.0)):
                ef = ns.fermi_levels[i, j]
                self.assertTrue(np.isfinite(ef))
                self.assertLess(abs(fd.get_doping(ef, T) / _conc(cm3) - 1), 1e-8)
                self.assertLess(abs((ns.hole_conc[i, j] - ns.electron_conc[i, j]) / _conc(cm3) - 1), 1e-6)
            self.assertLess(ns.fermi_levels[i, 0] / EH, CBM)            # 非简并：E_F 在导带边以下
            self.assertGreater(ns.fermi_levels[i, 0] / EH, CBM - 0.2)
        for rec in FF.STATE["fallbacks"][n0:]:
            self.assertAlmostEqual(rec["tol"], 1e-5, delta=1e-12)        # 第一级就满足，不往下退

    def test_unchanged_where_original_solves(self):
        fd = _ws2_dos()
        for cm3 in (-1e20, -1.86e18, 3.31e17, 1e16):
            for T in (300.0, 600.0):
                c = _conc(cm3)
                with np.errstate(all="ignore"):
                    a = self.orig(fd, c, T, tol=1e-5, precision=10, return_electron_hole_conc=True)
                    b = fd.get_fermi(c, T, tol=1e-5, precision=10, return_electron_hole_conc=True)
                    self.assertEqual(self.orig(fd, c, T, tol=1e-5, precision=10),
                                     fd.get_fermi(c, T, tol=1e-5, precision=10))
                self.assertEqual(tuple(map(float, a)), tuple(map(float, b)))

    def test_nan_raises(self):
        fd = _ws2_dos()
        fd.efermi = float("nan")
        c = _conc(-1.86e18)
        with np.errstate(all="ignore"):
            self.assertTrue(np.isnan(self.orig(fd, c, 300.0, tol=1e-5, precision=10)))
            with self.assertRaises(ValueError) as cm:
                fd.get_fermi(c, 300.0, tol=1e-5, precision=10)
        self.assertIn("NaN", str(cm.exception))
        with _Env("0"), np.errstate(all="ignore"), self.assertRaises(ValueError):
            fd.get_fermi(c, 300.0, tol=1e-5, precision=10)            # NaN 守卫不受开关影响

    def test_loose_solution_refused(self):
        fd = _ws2_dos(cb_width=0.5)                                    # 导带只装得下 0.1 个电子/胞
        c = -0.3 / fd.structure.volume                                 # 要 0.3 个 -> 无解
        with np.errstate(all="ignore"):
            ef = float(self.orig(fd, c, 300.0, tol=1.0, precision=10))  # 原式在 tol=1 接受
            err = abs(fd.get_doping(ef, 300.0) / c - 1)
        self.assertGreater(err, 0.5)
        self.assertIsNone(FF.solve_fermi(fd, c, 300.0))
        with _Env(None), np.errstate(all="ignore"):
            for tol in np.logspace(-5, 0, 6):
                with self.assertRaises(ValueError):
                    fd.get_fermi(c, 300.0, tol=tol, precision=10)
        with _Env("off"), np.errstate(all="ignore"):
            self.assertEqual(float(fd.get_fermi(c, 300.0, tol=1.0, precision=10)), ef)

    def test_switch_off_restores_original(self):
        fd = _ws2_dos()
        c = _conc(-1.86e18)
        with _Env("0"), np.errstate(all="ignore"), self.assertRaises(ValueError):
            fd.get_fermi(c, 100.0, tol=1e-5, precision=10)
        with _Env(None), np.errstate(all="ignore"):
            self.assertTrue(np.isfinite(fd.get_fermi(c, 100.0, tol=1e-5, precision=10)))

    def test_on_top_of_v131_plugin_wrapper(self):
        """旧运行目录里的 amset2d_plugin（V131）先套了 NaN 壳：插件照样挂上，unwrap 找得到 AMSET 原函数。"""
        import json
        import subprocess
        import tempfile
        code = '''
import json, numpy as np
for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
    if not hasattr(np, _a):
        setattr(np, _a, getattr(np, _b))
from amset.electronic_structure.dos import FermiDos as _FD
_orig_get_fermi = _FD.get_fermi
def _get_fermi_no_nan(self, concentration, temperature, *a, **k):
    out = _orig_get_fermi(self, concentration, temperature, *a, **k)
    ef = out[0] if isinstance(out, tuple) else out
    if not np.isfinite(ef):
        raise ValueError("[amset2d] get_fermi NaN")
    return out
_get_fermi_no_nan._amset2d_nan = True
_FD.get_fermi = _get_fermi_no_nan
import amset_fermi_fix as FF
import test_fermi_fix as T
fd = T._ws2_dos()
c = T._conc(-1.86e18)
with np.errstate(all="ignore"):
    ef, ne, nh = fd.get_fermi(c, 100.0, tol=1e-5, precision=10, return_electron_hole_conc=True)
print(json.dumps({"applied": FF.STATE["applied"], "orig": FF.original_get_fermi() is _orig_get_fermi,
                  "err": abs(fd.get_doping(ef, 100.0) / c - 1)}))
'''
        e = dict(os.environ)
        e["PYTHONPATH"] = os.pathsep.join([str(HERE), e.get("PYTHONPATH", "")])
        e.pop(FF.ENV_NAME, None)
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=e,
                           cwd=tempfile.mkdtemp(), timeout=600)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        out = json.loads(r.stdout.strip().splitlines()[-1])
        self.assertTrue(out["applied"])
        self.assertTrue(out["orig"])
        self.assertLess(out["err"], 1e-8)

    def test_apply_idempotent(self):
        self.assertTrue(FF.apply())
        self.assertTrue(FF.apply())
        fn = FermiDos.get_fermi
        self.assertTrue(getattr(fn, FF.MARK, False))
        self.assertFalse(getattr(fn.__wrapped__, FF.MARK, False))
        self.assertIs(FF.original_get_fermi(), fn.__wrapped__)


class WiringTests(unittest.TestCase):
    def test_gens_import_plugin(self):
        for rel in ("step8_amset/gen_step10_amset.py", "step8.4_amset2d/gen_step14_amset2d.py"):
            src = (KE / rel).read_text(encoding="utf-8")
            self.assertIn('("import amset_fermi_fix; " if _has_fermi else "")', src, rel)
            self.assertIn("_has_fermi = _install_fermi_fix(out)", src, rel)
            self.assertIn("kc.FERMI_FIX_PLUGIN", src, rel)
        kc = (KE / "ke_common.py").read_text(encoding="utf-8")
        self.assertIn('FERMI_FIX_PLUGIN = "amset_fermi_fix.py"', kc)

    def test_plugin_and_postprocess_load_it(self):
        self.assertIn("import amset_fermi_fix", (HERE / "amset2d_plugin.py").read_text(encoding="utf-8"))
        self.assertNotIn("_amset2d_nan", (HERE / "amset2d_plugin.py").read_text(encoding="utf-8"))
        pp = (HERE / "postprocess_intrinsic.py").read_text(encoding="utf-8")
        self.assertIn("_apply_fermi_fix(run_dir)", pp)

    def test_skill_yaml_and_zt_link(self):
        ke = (KE / "skill.yaml").read_text(encoding="utf-8")
        self.assertEqual(ke.count("amset_ir_fix.py, amset_fermi_fix.py]"), 2)
        zt = (SKILL / "zt-dft-cpu" / "skill.yaml").read_text(encoding="utf-8")
        self.assertEqual(zt.count("amset_ir_fix.py, amset_fermi_fix.py]"), 1)
        link = SKILL / "zt-dft-cpu" / "amset_fermi_fix.py"
        self.assertTrue(link.is_file())
        self.assertEqual(link.resolve(), (HERE / "amset_fermi_fix.py").resolve())


if __name__ == "__main__":
    unittest.main(verbosity=2)
