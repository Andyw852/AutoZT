#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_ir_fix_kzcap.py —— V125 两个提速补丁的检验（不需要 VASP；AMSET 部分在子进程里跑，互不污染）。

用法：python test_ir_fix_kzcap.py      退出码 0 = 全部 PASS。

1) amset_ir_fix：六方单层 / 纤锌矿的不可约 k 点数 = spglib 用正确晶格的结果（原 AMSET 多 5–6 倍）；
   fcc 原胞不受影响；AZ_IR_FIX=0 时保持原行为；源码对不上时拒绝打补丁。
2) amset2d_plugin 的 kz_cap：equivalence 截断、真空轴判定、"沿 k_z 平"核对；Interpolator 端到端：
   能带平 -> k_z 3 层（面内不变），能带沿 k_z 有色散 -> 不截断。
3) gen：KZ_CAP_2D 取值解析、与插件同一真空轴判据、_interp_mesh 预测 k_z=3、2d_correction.json 写/删键；
   S8/S8.4 命令里 import 插件并 export AZ_IR_FIX；ke_common.ir_fix_setting 默认开；gen_need 与 zt 软链齐全。
"""
import importlib.util
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
ROOT = HERE.parent
for _p in (str(HERE), str(ROOT), str(ROOT / "step8_amset"), str(ROOT.parent / "_common" / "opt"),
           str(ROOT / "step2_bandgap" / "step2.15_discriminant")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ke_common as kc  # noqa: E402

PRELUDE = textwrap.dedent('''
    import sys, warnings
    warnings.filterwarnings("ignore")
    import numpy as np
    if not hasattr(np, "Inf"):
        np.Inf = np.inf
    if not hasattr(np, "trapz"):
        np.trapz = np.trapezoid
    from pymatgen.core import Structure, Lattice
    def wse2(c=23.3):
        z = 1.65 / c
        return Structure(Lattice.hexagonal(3.32, c), ["W", "Se", "Se"],
                         [[0, 0, .5], [1/3, 2/3, .5 + z], [1/3, 2/3, .5 - z]])
    def gan():
        return Structure(Lattice.hexagonal(3.19, 5.19), ["Ga", "Ga", "N", "N"],
                         [[1/3, 2/3, 0], [2/3, 1/3, .5], [1/3, 2/3, .377], [2/3, 1/3, .877]])
    def si():
        return Structure(Lattice([[0, 2.73, 2.73], [2.73, 0, 2.73], [2.73, 2.73, 0]]), ["Si", "Si"],
                         [[0, 0, 0], [.25, .25, .25]])
    def spg_ir(st, mesh):
        import spglib
        cell = (st.lattice.matrix, st.frac_coords, [s.Z for s in st.species])
        m, _ = spglib.get_ir_reciprocal_mesh(mesh, cell, symprec=0.01, is_time_reversal=True)
        return len(np.unique(m))
''')


def _stub_dir():
    """本地环境没有 sumo 时（amset.interpolation 需要它）放一个空壳；集群环境有真的就不用。"""
    if importlib.util.find_spec("sumo") is not None:
        return None
    d = Path(tempfile.mkdtemp()) / "sumo"
    d.mkdir()
    (d / "__init__.py").write_text("")
    (d / "symmetry.py").write_text("class Kpath: pass\nclass PymatgenKpath(Kpath): pass\n")
    return str(d.parent)


def _run(code, env=None, cwd=None):
    e = dict(os.environ)
    paths = [str(HERE)] + ([_stub_dir()] if _stub_dir() else []) + [e.get("PYTHONPATH", "")]
    e["PYTHONPATH"] = os.pathsep.join(p for p in paths if p)
    e.update(env or {})
    r = subprocess.run([sys.executable, "-c", PRELUDE + textwrap.dedent(code)], capture_output=True,
                       text=True, env=e, cwd=cwd, timeout=900)
    if r.returncode != 0:
        raise AssertionError("子进程失败：\n%s\n%s" % (r.stdout[-3000:], r.stderr[-3000:]))
    return json.loads(r.stdout.strip().splitlines()[-1])


class IrFixTests(unittest.TestCase):
    CODE = '''
        import json
        import amset_ir_fix
        from amset.electronic_structure.kpoints import get_kpoints_tetrahedral
        out = {"applied": amset_ir_fix.STATE["applied"]}
        for name, st, mesh in (("wse2", wse2(), [24, 24, 3]), ("gan", gan(), [12, 12, 8]), ("si", si(), [12, 12, 12])):
            r = get_kpoints_tetrahedral(mesh, st, symprec=0.01)
            out[name] = [int(len(r[3])), int(spg_ir(st, mesh))]
        print(json.dumps(out))
    '''

    def test_fix_gives_true_ibz(self):
        out = _run(self.CODE, env={"AZ_IR_FIX": "1"})
        self.assertTrue(out["applied"])
        for name in ("wse2", "gan", "si"):
            self.assertEqual(out[name][0], out[name][1], (name, out[name]))

    def test_default_on_and_off_keeps_original(self):
        on = _run(self.CODE, env={"AZ_IR_FIX": ""})
        self.assertTrue(on["applied"])                                      # 未设 = 开
        off = _run(self.CODE, env={"AZ_IR_FIX": "0"})
        self.assertFalse(off["applied"])
        self.assertGreater(off["wse2"][0], 4 * off["wse2"][1])              # 原 AMSET：六方多 5–6 倍
        self.assertGreater(off["gan"][0], 3 * off["gan"][1])
        self.assertEqual(off["si"][0], off["si"][1])                        # fcc 原胞不受影响

    def test_refuses_when_source_differs(self):
        out = _run('''
            import json, os
            os.environ["AZ_IR_FIX"] = "0"
            import amset_ir_fix
            import amset.electronic_structure.kpoints as K
            def get_kpoints_tetrahedral(mesh, structure, symprec=0.01):   # 假装 AMSET 已改实现
                return None
            K.get_kpoints_tetrahedral = get_kpoints_tetrahedral
            print(json.dumps({"applied": amset_ir_fix.apply(force=True)}))
        ''')
        self.assertFalse(out["applied"])


class KzCapTests(unittest.TestCase):
    def _record_dir(self, rmax=1):
        d = Path(tempfile.mkdtemp())
        rec = json.loads((HERE / "example_2d_correction.json").read_text())
        if rmax:
            rec["kz_cap_rmax"] = rmax
        (d / "2d_correction.json").write_text(json.dumps(rec))
        return d

    def test_helpers(self):
        out = _run('''
            import json, spglib
            import amset2d_plugin as P
            st = wse2()
            eq = [np.array([[0, 0, 0]]), np.array([[1, 0, 0], [0, 1, 0]]), np.array([[0, 0, 1], [0, 0, -1]]),
                  np.array([[1, 0, 2], [1, 0, -2]])]
            kept = P.kz_cap_filter(eq, 2, 1)
            cell = (st.lattice.matrix, st.frac_coords, [s.Z for s in st.species])
            mp, grid = spglib.get_ir_reciprocal_mesh([45, 45, 3], cell, is_shift=[0, 0, 0])
            k = grid[np.unique(mp)] / 45.0
            k[:, 2] = grid[np.unique(mp)][:, 2] / 3.0
            flat = np.vstack([np.cos(2 * np.pi * k[:, 0]), 1 + np.cos(2 * np.pi * k[:, 1])])
            disp = flat + 0.02 * np.cos(2 * np.pi * k[:, 2])[None]
            tilted = Lattice([[3.3, 0, 0.3], [-1.65, 2.86, 0], [0, 0, 23.0]]).matrix   # a 出面 -> 不行
            slanted = Lattice([[3.3, 0, 0], [-1.65, 2.86, 0], [0.5, 0, 23.0]]).matrix  # 只有 c 倾斜 -> 可以
            print(json.dumps({"kept": len(kept), "axis": P.vacuum_axis(st.lattice.matrix),
                              "tilted": P.vacuum_axis(tilted), "slanted": P.vacuum_axis(slanted),
                              "flat": P.kz_spread_ev(k, flat, 2, -5, 5),
                              "disp": P.kz_spread_ev(k, disp, 2, -5, 5), "active": P._KZ["rmax"]}))
        ''', env={"AMSET2D_RECORD": str(self._record_dir() / "2d_correction.json")})
        self.assertEqual(out["kept"], 3)
        self.assertEqual(out["axis"], 2)
        self.assertIsNone(out["tilted"])
        self.assertEqual(out["slanted"], 2)
        self.assertLess(out["flat"][0], 1e-9)
        self.assertGreater(out["flat"][1], 50)            # IBZ 列表里确有同一面内坐标、不同 k_z 的点
        self.assertGreater(out["disp"][0], 0.02)
        self.assertEqual(out["active"], 1)

    INTERP = '''
        import json, spglib
        import amset2d_plugin as P
        from amset.interpolation.bandstructure import Interpolator
        from pymatgen.electronic_structure.bandstructure import BandStructure
        from pymatgen.electronic_structure.core import Spin
        st = wse2()
        cell = (st.lattice.matrix, st.frac_coords, [s.Z for s in st.species])
        mp, grid = spglib.get_ir_reciprocal_mesh([18, 18, 3], cell, is_shift=[0, 0, 0])
        k = grid[np.unique(mp)] / np.array([18, 18, 3], float)
        sh = [(1, 0), (0, 1), (1, 1)]
        s = sum(np.cos(2 * np.pi * (k[:, 0] * a + k[:, 1] * b)) for a, b in sh)
        ev = -0.6 - 0.35 * (3 - s)
        ec = 0.6 + 0.30 * (3 - s) + DISP * np.cos(2 * np.pi * k[:, 2])
        bs = BandStructure(k, {Spin.up: np.vstack([ev - 2, ev, ec, ec + 2.5])},
                           st.lattice.reciprocal_lattice, 0.0, structure=st)
        it = Interpolator(bs, num_electrons=4, interpolation_factor=20)
        print(json.dumps({"mesh": [int(x) for x in it.interpolation_mesh], "active": P._KZ["active"]}))
    '''

    def test_interpolator_flat_is_capped(self):
        d = self._record_dir()
        capped = _run(self.INTERP.replace("DISP", "0.0"), env={"AMSET2D_RECORD": str(d / "2d_correction.json")})
        d0 = self._record_dir(rmax=None)
        plain = _run(self.INTERP.replace("DISP", "0.0"), env={"AMSET2D_RECORD": str(d0 / "2d_correction.json")})
        self.assertTrue(capped["active"])
        self.assertEqual(capped["mesh"][2], 3)
        self.assertEqual(capped["mesh"][:2], plain["mesh"][:2])           # 面内网格不变
        self.assertGreater(plain["mesh"][2], 3)

    def test_interpolator_dispersive_not_capped(self):
        d = self._record_dir()
        out = _run(self.INTERP.replace("DISP", "0.05"), env={"AMSET2D_RECORD": str(d / "2d_correction.json")})
        self.assertFalse(out["active"])
        self.assertGreater(out["mesh"][2], 3)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class GenWiringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.G14 = _load("g14_v125", HERE / "gen_step14_amset2d.py")
        cls.G10 = _load("g10_v125", ROOT / "step8_amset" / "gen_step10_amset.py")

    def test_kz_cap_value_parsing(self):
        f = self.G14.kz_cap_rmax
        self.assertIsNone(f("off"))
        self.assertIsNone(f(""))
        self.assertIsNone(f("0"))
        self.assertEqual(f("on"), 1)
        self.assertEqual(f("2"), 2)
        self.assertIsNone(f("maybe"))
        self.assertEqual(self.G14.KZ_CAP_2D, "off")                        # 出厂关

    def test_vacuum_axis_same_as_plugin_and_interp_mesh(self):
        from pymatgen.core import Lattice, Structure
        st = Structure(Lattice.hexagonal(3.32, 23.3), ["W", "Se", "Se"],
                       [[0, 0, .5], [1 / 3, 2 / 3, .57], [1 / 3, 2 / 3, .43]])
        self.assertEqual(self.G14.vacuum_axis(st.lattice.matrix), 2)
        self.assertIsNone(self.G14.vacuum_axis([[3.3, 0, 0.3], [-1.65, 2.86, 0], [0, 0, 23.0]]))
        self.assertEqual(self.G14.vacuum_axis([[3.3, 0, 0], [-1.65, 2.86, 0], [0.5, 0, 23.0]]), 2)
        old = self.G14._KZ_CAP_RMAX
        try:
            self.G14._KZ_CAP_RMAX = None
            plain = self.G14._interp_mesh(st, 384, 20)
            self.G14._KZ_CAP_RMAX = 1
            capped = self.G14._interp_mesh(st, 384, 20)
        finally:
            self.G14._KZ_CAP_RMAX = old
        self.assertEqual(list(capped[:2]), list(plain[:2]))
        self.assertEqual(int(capped[2]), 3)
        self.assertGreater(int(plain[2]), 3)

    def test_record_kz_cap_writes_and_removes(self):
        d = Path(tempfile.mkdtemp())
        (d / "2d_correction.json").write_text(json.dumps({"cell_c_A": 23.3}))
        old = self.G14._KZ_CAP_RMAX
        try:
            self.G14._KZ_CAP_RMAX = 1
            self.G14.record_kz_cap(d)
            rec = json.loads((d / "2d_correction.json").read_text())
            self.assertEqual(rec["kz_cap_rmax"], 1)
            self.assertIn("kz_flat_tol_ev", rec)
            self.G14._KZ_CAP_RMAX = None
            self.G14.record_kz_cap(d)
            rec = json.loads((d / "2d_correction.json").read_text())
            self.assertNotIn("kz_cap_rmax", rec)
            self.assertEqual(rec["cell_c_A"], 23.3)
        finally:
            self.G14._KZ_CAP_RMAX = old

    def test_ir_fix_setting_and_commands(self):
        old = os.environ.pop(kc.IR_FIX_ENV, None)
        try:
            self.assertTrue(kc.ir_fix_setting(None)[0])
            self.assertTrue(kc.ir_fix_setting("auto")[0])
            self.assertFalse(kc.ir_fix_setting("off")[0])
            os.environ[kc.IR_FIX_ENV] = "0"
            self.assertFalse(kc.ir_fix_setting("auto")[0])
            self.assertTrue(kc.ir_fix_setting("on")[0])
        finally:
            os.environ.pop(kc.IR_FIX_ENV, None)
            if old is not None:
                os.environ[kc.IR_FIX_ENV] = old
        self.assertEqual(kc.ir_fix_cmd_prefix(False), "export AZ_IR_FIX=0; ")
        for mod in (self.G10, self.G14):
            src = Path(mod.__file__).read_text(encoding="utf-8")
            self.assertIn('"import amset_ir_fix; " if _has_ir', src)
            self.assertIn("kc.ir_fix_cmd_prefix(_IR_FIX_ON)", src)
            self.assertIn("IR_FIX", mod.SPEC)
        self.assertIn("KZ_CAP_2D", self.G14.SPEC)

    def test_gen_need_and_zt_symlink(self):
        import yaml
        ke = (ROOT / "skill.yaml").read_text(encoding="utf-8")
        zt = (ROOT.parent / "zt-dft-cpu" / "skill.yaml").read_text(encoding="utf-8")
        yaml.safe_load(ke)
        yaml.safe_load(zt)
        line10 = next(ln for ln in ke.splitlines() if "gen_step10_amset.py" in ln and "gen_need" in ln)
        self.assertIn("amset_ir_fix.py", line10)
        blk = ke[ke.index("gen: gen_step14_amset2d.py"):]
        self.assertIn("amset_ir_fix.py", blk[:blk.index("fetch_all")])
        zblk = zt[zt.index("gen_step10_amset.py"):]
        self.assertIn("amset_ir_fix.py", zblk[:zblk.index("fetch_all")])
        link = ROOT.parent / "zt-dft-cpu" / "amset_ir_fix.py"
        self.assertTrue(link.is_symlink() and link.resolve() == (HERE / "amset_ir_fix.py").resolve())


if __name__ == "__main__":
    unittest.main(verbosity=2)
