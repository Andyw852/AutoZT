#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_star_audit.py —— V130：tools/star_audit.py 的拆分公式 + amset2d_plugin 的输运前守卫（不需要 VASP）。

用法：python test_star_audit.py      退出码 0 = 全部 PASS。

1) 合成六方 mesh 对（同一张 12×12×1 网格；旧 = AMSET 原式不可约映射，按 H 轨道给随机率，全网格顺序打乱；
   新 = IR_FIX 映射）：
   · 新 = 旧成员率的星平均 -> 仿真项 = 0，Jensen 项 < 0 且等于直接算的值；
   · 新 = 旧在代表点上的率 -> 新/旧 = 代表点/旧；
   · 网格不同 -> 报错。
2) 守卫：transport_report 找得到 NaN 率、零率、NaN 费米能级；包装后的输运入口遇到 LinAlgError 时
   带着坏值说明改抛 RuntimeError。
"""
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent / "tools"
for _p in (str(HERE), str(TOOLS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import h5py
    import spglib
    from pymatgen.core import Lattice, Structure
    import star_audit as SA
    _HAS = True
except ImportError:                                # noqa: BLE001
    _HAS = False


def _mesh_pair(tmp, rng, dims=(12, 12, 1), nb=2):
    st = Structure(Lattice.hexagonal(3.19, 20.0), ["Mo", "S", "S"],
                   [[0, 0, .5], [1 / 3, 2 / 3, .578], [1 / 3, 2 / 3, .422]])
    cell = (st.lattice.matrix, st.frac_coords, [s.Z for s in st.species])
    cell_t = (np.ascontiguousarray(st.lattice.matrix.T), st.frac_coords, [s.Z for s in st.species])
    m_new, grid = spglib.get_ir_reciprocal_mesh(list(dims), cell, symprec=0.01)
    m_old, grid_o = spglib.get_ir_reciprocal_mesh(list(dims), cell_t, symprec=0.01)
    assert np.array_equal(grid, grid_o)
    k = grid / np.array(dims, float)
    u_new, inv_new = np.unique(m_new, return_inverse=True)
    u_old, inv_old = np.unique(m_old, return_inverse=True)
    assert len(u_new) < len(u_old)
    g_old_ir = 1e13 * (1 + 0.4 * rng.random((nb, len(u_old))))           # 旧：每个 H 轨道一个率
    g_old_full = g_old_ir[:, inv_old]
    cnt = np.bincount(inv_new).astype(float)
    star_mean = np.stack([np.bincount(inv_new, g_old_full[b]) / cnt for b in range(nb)])
    e = 0.05 * rng.random((nb, len(u_new))) + np.array([[-0.02], [0.03]])
    v = rng.normal(size=(nb, len(u_new), 3)) * 1e-3

    def write(path, kpts, ir_idx, inv, rates_ir):
        with h5py.File(str(path), "w") as f:
            f["kpoints"] = kpts
            f["ir_kpoints"] = kpts[ir_idx]
            f["ir_to_full_kpoint_mapping"] = inv
            f["scattering_rates_up"] = rates_ir[None, None, None]
            f["energies_up"] = e if rates_ir.shape[1] == len(u_new) else np.zeros(rates_ir.shape)
            f["velocities_up"] = v if rates_ir.shape[1] == len(u_new) else np.zeros(rates_ir.shape + (3,))
            f["fermi_levels"] = np.array([[0.0]])
            f["temperatures"] = np.array([300.0])
            f["doping"] = np.array([1e13])
            f["scattering_labels"] = np.array([b"ADP"])

    # 旧文件：全网格顺序打乱
    sh = rng.permutation(len(k))
    back = np.argsort(sh)
    ir_old_full_idx = np.array([np.where(inv_old == i)[0][0] for i in range(len(u_old))])
    write(tmp / "old.h5", k[sh], back[ir_old_full_idx], inv_old[sh], g_old_ir)
    ir_new_full_idx = np.array([np.where(inv_new == i)[0][0] for i in range(len(u_new))])
    write(tmp / "new_star.h5", k, ir_new_full_idx, inv_new, star_mean)
    write(tmp / "new_rep.h5", k, ir_new_full_idx, inv_new, g_old_full[:, ir_new_full_idx])
    return dict(inv_new=inv_new, g_old_full=g_old_full, star_mean=star_mean, e=e, v=v, cnt=cnt)


@unittest.skipUnless(_HAS, "需要 h5py + spglib + pymatgen")
class AuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        cls.d = _mesh_pair(cls.tmp, np.random.default_rng(5))

    def _direct_jensen(self):
        d = self.d
        kt = SA.KB_EV * 300.0
        f = 1 / (np.exp(d["e"] / kt) + 1)
        w = f * (1 - f) / kt * (d["v"][..., 0] ** 2 + d["v"][..., 1] ** 2) / 2 * d["cnt"][None]
        minv = np.stack([np.bincount(d["inv_new"], 1 / d["g_old_full"][b]) / d["cnt"] for b in range(2)])
        return 100 * (np.sum(w / d["star_mean"]) / np.sum(w * minv) - 1)

    def test_perfect_star_mean(self):
        r = SA.audit(self.tmp / "old.h5", self.tmp / "new_star.h5")["rows"][0]
        self.assertAlmostEqual(r["emulation_pct"], 0.0, places=9)
        self.assertLess(r["jensen_pct"], 0.0)
        self.assertAlmostEqual(r["jensen_pct"], self._direct_jensen(), places=9)
        self.assertAlmostEqual(r["new_vs_old_pct"], r["jensen_pct"], places=9)
        self.assertAlmostEqual(r["gamma_new_over_star_mean"]["p50"], 1.0, places=9)
        self.assertGreater(r["star_spread_rel"], 0.01)

    def test_rep_only(self):
        r = SA.audit(self.tmp / "old.h5", self.tmp / "new_rep.h5")["rows"][0]
        self.assertAlmostEqual(r["new_vs_old_pct"], r["rep_vs_old_pct"], places=9)

    def test_main_and_mismatch(self):
        out = self.tmp / "a.json"
        self.assertEqual(SA.main([str(self.tmp / "old.h5"), str(self.tmp / "new_star.h5"), "--json", str(out)]), 0)
        self.assertEqual(json.loads(out.read_text())["mech"], "ADP")
        with h5py.File(str(self.tmp / "other.h5"), "w") as f:
            for k in ("ir_kpoints", "ir_to_full_kpoint_mapping", "fermi_levels", "temperatures", "doping",
                      "scattering_labels", "scattering_rates_up"):
                h5py.File(str(self.tmp / "old.h5"), "r").copy(k, f)
            f["kpoints"] = np.zeros((10, 3))
        self.assertEqual(SA.main([str(self.tmp / "other.h5"), str(self.tmp / "new_star.h5")]), 2)


def _env():
    import importlib.util
    e = dict(os.environ)
    paths = [str(HERE)]
    if importlib.util.find_spec("sumo") is None:
        d = Path(tempfile.mkdtemp()) / "sumo"
        d.mkdir()
        (d / "__init__.py").write_text("")
        (d / "symmetry.py").write_text("class Kpath: pass\nclass PymatgenKpath(Kpath): pass\n")
        paths.append(str(d.parent))
    e["PYTHONPATH"] = os.pathsep.join(paths + [e.get("PYTHONPATH", "")])
    e["AMSET2D_RECORD"] = str(HERE / "example_2d_correction.json")
    e["AZ_IR_FIX"] = "0"
    return e


class GuardTests(unittest.TestCase):
    CODE = textwrap.dedent('''
        import json, sys, types, numpy as np
        for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
            if not hasattr(np, _a):
                setattr(np, _a, getattr(np, _b))
        import amset2d_plugin as P
        import amset.core.transport as T
        from pymatgen.electronic_structure.core import Spin
        r = np.full((2, 1, 1, 2, 5), 1e13)
        r[0, 0, 0, 1, 3] = np.nan
        r[1, 0, 0, 0, 2] = 0.0
        ad = types.SimpleNamespace(fermi_levels=np.array([[np.nan]]), scattering_labels=["ADP", "IMP"],
                                   kpoints=np.zeros((5, 3)), scattering_rates={Spin.up: r},
                                   energies={Spin.up: np.zeros((2, 5))},
                                   velocities_product={Spin.up: np.zeros((2, 3, 3, 5))})
        rep = P.transport_report(ad)
        def boom(amset_data):
            raise np.linalg.LinAlgError("SVD did not converge")
        g = P._guard_transport(boom)
        try:
            g(ad)
            err = None
        except RuntimeError as e:
            err = str(e)
        print(json.dumps({"rep": rep, "err": err,
                          "wrapped": [getattr(getattr(T, n), "_amset2d_guard", False)
                                      for n in ("_calculate_transport_properties", "_calculate_mobility")]}))
    ''')

    def test_report_and_wrap(self):
        r = subprocess.run([sys.executable, "-c", self.CODE], capture_output=True, text=True, env=_env(),
                           cwd=tempfile.mkdtemp(), timeout=600)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        out = json.loads(r.stdout.strip().splitlines()[-1])
        txt = "\n".join(out["rep"])
        self.assertIn("fermi_levels", txt)
        self.assertIn("ADP[up]", txt)
        self.assertIn("IMP[up]", txt)
        self.assertIn("SVD did not converge", out["err"])
        self.assertIn("ADP[up]", out["err"])
        self.assertEqual(out["wrapped"], [True, True])


if __name__ == "__main__":
    unittest.main(verbosity=2)
