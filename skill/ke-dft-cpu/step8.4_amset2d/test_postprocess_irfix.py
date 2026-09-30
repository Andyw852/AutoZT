#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_postprocess_irfix.py —— V128：postprocess_intrinsic 按 mesh.h5 判定 IR_FIX，重建用同一套不可约映射。

用法：python test_postprocess_irfix.py      退出码 0 = 全部 PASS。

MoS₂ IR_FIX-only 实测：AMSET 本体跑完，postprocess 的复现检查 FAIL（"mesh/settings/vasprun 不同源"）。
原因：重积分的 transport DOS 先把积分量加到不可约点上再做四面体积分，重建时的不可约映射
（不挂 amset_ir_fix = AMSET 原式）与作业（IR_FIX）不同，结果差 ~1e-1，远超 1e-6 的复现容差。

1) mesh_ir_fix_state：六方单层按不可约点数判 True/False，立方 None，都不对 -> ValueError。
2) _apply_ir_fix：True -> 设 AZ_IR_FIX=1 并挂上运行目录的 amset_ir_fix；False -> AZ_IR_FIX=0 不挂；
   缺文件 -> 停。
3) 端到端（合成纤锌矿，AMSET 真跑 IR_FIX + write_mesh）：重建不挂修正 -> apply_mesh_metadata 拒绝
   （不可约点数不同）；挂上 -> 重积分与 transport.json 差 ≤ 1e-6。
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
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import postprocess_intrinsic as P  # noqa: E402

try:
    import spglib  # noqa: F401
    from pymatgen.core import Lattice, Structure
    _HAS = True
except ImportError:                                # noqa: BLE001
    _HAS = False


def mos2(c=20.0):
    z = 1.56 / c
    return Structure(Lattice.hexagonal(3.19, c), ["Mo", "S", "S"],
                     [[0, 0, .5], [1 / 3, 2 / 3, .5 + z], [1 / 3, 2 / 3, .5 - z]])


def si():
    return Structure(Lattice([[0, 2.73, 2.73], [2.73, 0, 2.73], [2.73, 2.73, 0]]), ["Si", "Si"],
                     [[0, 0, 0], [.25, .25, .25]])


def fake_mesh(st, dims, n_ir):
    g = np.array(np.meshgrid(*[np.arange(n) for n in dims], indexing="ij")).reshape(3, -1).T / np.array(dims, float)
    g -= np.rint(g)
    return {"structure": st, "kpoints": g, "ir_kpoints": np.zeros((n_ir, 3)), "soc": False}


def n_ir(st, dims, transposed):
    L = np.ascontiguousarray(st.lattice.matrix.T if transposed else st.lattice.matrix)
    m, _ = spglib.get_ir_reciprocal_mesh(dims, (L, st.frac_coords, [s.Z for s in st.species]), symprec=0.01)
    return len(np.unique(m))


@unittest.skipUnless(_HAS, "需要 spglib + pymatgen")
class StateTests(unittest.TestCase):
    def test_hexagonal(self):
        st, dims = mos2(), [24, 24, 3]
        nf, no = n_ir(st, dims, False), n_ir(st, dims, True)
        self.assertLess(nf, no)
        self.assertIs(P.mesh_ir_fix_state(fake_mesh(st, dims, nf), 0.01), True)
        self.assertIs(P.mesh_ir_fix_state(fake_mesh(st, dims, no), 0.01), False)
        with self.assertRaises(ValueError):
            P.mesh_ir_fix_state(fake_mesh(st, dims, nf + 1), 0.01)

    def test_cubic_none(self):
        st, dims = si(), [8, 8, 8]
        self.assertIsNone(P.mesh_ir_fix_state(fake_mesh(st, dims, n_ir(st, dims, False)), 0.01))


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
    e.pop("AZ_IR_FIX", None)
    return e


PRE = textwrap.dedent('''
    import json, os, sys, numpy as np
    for _a, _b in (("Inf", "inf"), ("trapz", "trapezoid"), ("string_", "bytes_")):
        if not hasattr(np, _a):
            setattr(np, _a, getattr(np, _b))
''')


def _py(code, cwd, timeout=900):
    r = subprocess.run([sys.executable, "-c", PRE + textwrap.dedent(code)], capture_output=True, text=True,
                       env=_env(), cwd=cwd, timeout=timeout)
    return r


class ApplyTests(unittest.TestCase):
    CODE = '''
        import postprocess_intrinsic as P
        try:
            P._apply_ir_fix(sys.argv[1] if len(sys.argv) > 1 else ".", %s)
        except SystemExit as e:
            print(json.dumps({"exit": str(e)[:30]})); sys.exit(0)
        import amset.electronic_structure.kpoints as K
        m = sys.modules.get("amset_ir_fix")
        print(json.dumps({"env": os.environ.get("AZ_IR_FIX"), "loaded": m is not None,
                          "proxied": type(K.spglib).__name__ == "_SpglibProxy"}))
    '''

    def _case(self, state, with_file=True):
        d = Path(tempfile.mkdtemp())
        if with_file:
            shutil.copyfile(HERE / "amset_ir_fix.py", d / "amset_ir_fix.py")
        r = _py(self.CODE % state, cwd=str(d))
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        return json.loads(r.stdout.strip().splitlines()[-1])

    def test_true_loads(self):
        self.assertEqual(self._case("True"), {"env": "1", "loaded": True, "proxied": True})

    def test_false_does_not(self):
        self.assertEqual(self._case("False"), {"env": "0", "loaded": False, "proxied": False})

    def test_true_without_file_stops(self):
        self.assertIn("exit", self._case("True", with_file=False))


E2E = PRE + '''
import glob, warnings
warnings.filterwarnings("ignore")
from pymatgen.core import Lattice, Structure
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.electronic_structure.core import Spin
import spglib


def setup():
    st = Structure(Lattice.hexagonal(3.19, 5.19), ["Ga", "Ga", "N", "N"],
                   [[1/3, 2/3, 0], [2/3, 1/3, .5], [1/3, 2/3, .377], [2/3, 1/3, .877]])
    m = [12, 12, 8]
    mp, grid = spglib.get_ir_reciprocal_mesh(m, (st.lattice.matrix, st.frac_coords, [s.Z for s in st.species]))
    k = grid[np.unique(mp)] / np.array(m, float)
    rl, L = st.lattice.reciprocal_lattice.matrix, st.lattice.matrix
    kc = k @ rl
    s = sum(np.cos(kc @ (np.array(r, float) @ L)) for r in ((1, 0, 0), (0, 1, 0), (1, 1, 0)))
    s = s + 0.7 * np.cos(kc @ L[2])
    s0 = 3.7
    E = np.vstack([-0.3 * (s0 - s) - 2, -0.3 * (s0 - s), 1.5 + 0.35 * (s0 - s), 4.0 + 0.35 * (s0 - s)])
    C = np.zeros((6, 6)); C[0, 0] = C[1, 1] = 390.; C[0, 1] = C[1, 0] = 145.
    C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = 106.; C[2, 2] = 398.; C[3, 3] = C[4, 4] = 105.; C[5, 5] = 122.5
    return BandStructure(k, {Spin.up: E}, st.lattice.reciprocal_lattice, 0.0, structure=st), C


SET = {"scattering_type": ["ADP"], "deformation_potential": (3.0, 3.0), "energy_cutoff": 10.0,
       "doping": [-1e18, 1e18], "temperatures": [300], "unity_overlap": True, "interpolation_factor": 3,
       "nworkers": 2, "write_input": False, "write_mesh": True, "file_format": "json", "print_log": False,
       "high_frequency_dielectric": np.eye(3) * 5, "static_dielectric": np.eye(3) * 9}
'''

E2E_RUN = E2E + '''
if __name__ == "__main__":
    os.environ["AZ_IR_FIX"] = "1"
    import amset_ir_fix
    import amset.core.run as R
    R._log_structure_information = lambda *a, **k: None
    bs, C = setup()
    R.Runner(bs, 8, dict(SET, elastic_constant=C.tolist())).run(directory=".")
'''

E2E_CHECK = E2E + '''
if __name__ == "__main__":
    import postprocess_intrinsic as P
    apply = sys.argv[1] == "1"
    mesh = P.read_mesh_h5(glob.glob("mesh*.h5")[0])
    state = P.mesh_ir_fix_state(mesh, 0.01)
    if apply:
        P._apply_ir_fix(".", state)
    else:
        os.environ["AZ_IR_FIX"] = "0"
    from amset.interpolation.bandstructure import Interpolator
    from amset.core.transport import solve_boltzman_transport_equation
    bs, C = setup()
    ad = Interpolator(bs, num_electrons=8, interpolation_factor=3).get_amset_data(
        energy_cutoff=10.0, symprec=0.01, nworkers=2)
    ad.calculate_dos(estep=0.01)
    ad.set_doping_and_temperatures(np.array(SET["doping"], float), np.array(SET["temperatures"], float))
    out = {"state": state}
    try:
        perm = P.apply_mesh_metadata(ad, mesh)
    except ValueError as e:
        out["rejected"] = "不可约点数" in str(e)
        print(json.dumps(out)); sys.exit(0)
    P.inject_mesh_rates(ad, mesh, perm=perm)
    sig, seeb, kap, mob = solve_boltzman_transport_equation(ad, calculate_mobility=True, separate_mobility=True,
                                                            progress_bar=False)
    j = json.load(open(glob.glob("transport*.json")[0]))
    out["max_diff"] = max(float(np.nanmax(np.abs(np.asarray(mob["overall"]) - np.asarray(j["mobility"]["overall"])))),
                          float(np.nanmax(np.abs(sig - np.asarray(j["conductivity"])))),
                          float(np.nanmax(np.abs(seeb - np.asarray(j["seebeck"])))))
    print(json.dumps(out))
'''


@unittest.skipUnless(_HAS, "需要 spglib + pymatgen")
class EndToEndTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = Path(tempfile.mkdtemp())
        shutil.copyfile(HERE / "amset_ir_fix.py", cls.d / "amset_ir_fix.py")
        (cls.d / "run_e2e.py").write_text(E2E_RUN)
        (cls.d / "check_e2e.py").write_text(E2E_CHECK)
        e = _env()
        e["PYTHONPATH"] = os.pathsep.join([str(cls.d), e["PYTHONPATH"]])
        cls.env = e
        r = subprocess.run([sys.executable, "run_e2e.py"], capture_output=True, text=True, env=e,
                           cwd=str(cls.d), timeout=1500)
        assert r.returncode == 0, r.stderr[-3000:]

    def _check(self, apply):
        r = subprocess.run([sys.executable, "check_e2e.py", apply], capture_output=True, text=True, env=self.env,
                           cwd=str(self.d), timeout=1500)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        return json.loads(r.stdout.strip().splitlines()[-1])

    def test_without_fix_rejected(self):
        out = self._check("0")
        self.assertIs(out["state"], True)
        self.assertTrue(out.get("rejected"))

    def test_with_fix_reproduces(self):
        out = self._check("1")
        self.assertIs(out["state"], True)
        self.assertLessEqual(out["max_diff"], 1e-6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
