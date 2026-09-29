#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_dp_symmetrize.py —— V121 形变势场 D(k) 点群对称化的检验（不需要 VASP）。

用法：python test_dp_symmetrize.py      退出码 0 = 全部 PASS。

1) 约定（最关键）：用紧束缚模型对晶格施加 6 个应变分量做有限差分，得到**物理上**满足
   D(r k) = r D(k) rᵀ 的形变势场；dp_symmetrize 作用在它上面必须几乎不改值（< 1e-8），
   而故意写错的约定（张量用 rᵀ 旋转）必须明显改值 —— 2D MoS2、3D GaN 各一例。
2) 注入 4% 的 xx/yy 噪声（模拟 ε_xx、ε_yy 独立有限差分）：K 点（有 C3）投影回 xx=yy、xy=0；幂等。
3) h5：原地修改、保留已有 attrs（如 nspin_norm_fix 的标记）、打标记后重复调用跳过。
4) band_edges 的 (xx, yy) 投影：六方取平均，正交晶系不变。
5) S7.1 的两段嵌入脚本能编译，且在 nspin 修正之后、读 h5 / 闸门之前调用对称化。
"""
import re
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(HERE / "step8.4_amset2d"), str(HERE.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import dp_symmetrize as P                                  # noqa: E402
import test_symmetry_gate as tsg                           # noqa: E402

try:
    import test_desym_fix as TDF                           # tb_coeffs / full_mesh
    _HAS_TB = True
except Exception:                                          # noqa: BLE001
    _HAS_TB = False


def tb_D(st, kp, delta=1e-5):
    """紧束缚模型有限差分形变势（带符号张量应变，刚性离子），shape (nb, nk, 3, 3)。
    delta 取小：±δ 下的能移要远小于带间距，否则按能量排序的带在应变下交叉。"""
    from pymatgen.core import Lattice, Structure
    g0 = np.zeros((1, 3), int)

    def energies(eps):
        s2 = Structure(Lattice(st.lattice.matrix @ (np.eye(3) + eps).T), st.species, st.frac_coords)
        return TDF.tb_coeffs(s2, kp, g0)[1]
    nb = energies(np.zeros((3, 3))).shape[0]
    D = np.zeros((nb, len(kp), 3, 3))
    for i in range(3):
        for j in range(i, 3):
            e = np.zeros((3, 3))
            e[i, j] = e[j, i] = delta
            dE = (energies(e) - energies(-e)) / (2 * delta)
            D[:, :, i, j] = D[:, :, j, i] = dE / (1 if i == j else 2)
    return D, energies(np.zeros((3, 3)))


@unittest.skipUnless(_HAS_TB, "需要 test_desym_fix 的紧束缚模型（pymatgen/amset 环境）")
class ConventionTests(unittest.TestCase):
    CASES = (("MoS2 2D", tsg.mos2_aligned(), [6, 6, 1]), ("GaN 3D", tsg.gan_std(), [6, 6, 4]))

    def test_physical_field_unchanged_wrong_convention_changed(self):
        for name, st, mesh in self.CASES:
            with self.subTest(name=name):
                kp = TDF.full_mesh(mesh)
                D, en = tb_D(st, kp)
                rots = P.laue_rotations(st)
                maps = P.kpoint_maps(kp, st, rots)
                self.assertEqual(int((maps < 0).sum()), 0)
                gap = np.min(np.abs(np.diff(np.sort(en, axis=0), axis=0)), axis=0)
                ok = gap > 1e-2
                for o in range(len(rots)):
                    ok &= gap[maps[o]] > 1e-2
                self.assertGreater(ok.sum(), len(kp) // 3)
                Ds, _ = P.symmetrize_field(D, maps, rots)
                rel = np.abs(Ds - D)[:, ok].max() / np.abs(D)[:, ok].max()
                self.assertLess(rel, 1e-8, "%s 物理场被改了 %.2e" % (name, rel))
                Db, _ = P.symmetrize_field(D, maps, [r.T for r in rots])
                relb = np.abs(Db - D)[:, ok].max() / np.abs(D)[:, ok].max()
                self.assertGreater(relb, 1e-2, "%s 错误约定没被区分出来" % name)

    def test_noise_projected_and_idempotent(self):
        st, mesh = tsg.mos2_aligned(), [6, 6, 1]
        kp = TDF.full_mesh(mesh)
        D, _ = tb_D(st, kp)
        D[..., 0, 0] *= 1.04                                  # ε_xx 独立差分的噪声
        rots = P.laue_rotations(st)
        maps = P.kpoint_maps(kp, st, rots)
        Ds, _ = P.symmetrize_field(D, maps, rots)
        K = [i for i, k in enumerate(kp) if np.allclose(np.mod(k[:2], 1), [1 / 3, 1 / 3])]
        self.assertTrue(K)
        self.assertLess(np.abs(Ds[:, K, 0, 0] - Ds[:, K, 1, 1]).max(), 1e-10)
        self.assertLess(np.abs(Ds[:, K, 0, 1]).max(), 1e-10)
        Ds2, _ = P.symmetrize_field(Ds, maps, rots)
        self.assertLess(np.abs(Ds2 - Ds).max(), 1e-12)


class H5AndPairTests(unittest.TestCase):
    def test_h5_inplace_keeps_attrs_and_is_idempotent(self):
        import h5py
        st = tsg.mos2_aligned()
        g = np.array(list(np.ndindex(6, 6, 1)), float) / np.array([6, 6, 1])
        kp = g - np.rint(g)
        rng = np.random.default_rng(0)
        D = rng.normal(size=(3, len(kp), 3, 3))
        D = 0.5 * (D + np.swapaxes(D, -1, -2))
        p = Path(tempfile.mkdtemp()) / "deformation.h5"
        with h5py.File(str(p), "w") as f:                     # 与 amset write_deformation_potentials 同布局
            f.create_dataset("deformation_potentials_up", data=D)
            f["structure"] = np.bytes_(st.to_json())
            f["kpoints"] = kp
            f.attrs["nspin_norm_fixed"] = 1                    # 模拟 nspin_norm_fix 的标记
        self.assertFalse(P.is_symmetrized(p))
        rep = P.symmetrize_h5(p)
        self.assertEqual(rep["k_missing_images"], 0)
        self.assertGreater(rep["datasets"]["deformation_potentials_up"], 0.1)
        self.assertTrue(P.is_symmetrized(p))
        with h5py.File(str(p), "r") as f:
            self.assertEqual(int(f.attrs["nspin_norm_fixed"]), 1)
            first = np.array(f["deformation_potentials_up"])
        self.assertIn("skipped", P.symmetrize_h5(p))
        with h5py.File(str(p), "r") as f:
            self.assertTrue(np.array_equal(first, np.array(f["deformation_potentials_up"])))

    def test_inplane_pair(self):
        from pymatgen.core import Lattice, Structure
        hexa = tsg.mos2_aligned()
        xs, ys = P.symmetrize_inplane_pair(8.69, 8.35, hexa)
        self.assertAlmostEqual(xs, 8.52, places=6)
        self.assertAlmostEqual(ys, 8.52, places=6)
        ortho = Structure(Lattice.orthorhombic(3.0, 4.0, 15.0), ["Mo", "S"], [[0, 0, .5], [.5, .5, .6]])
        xo, yo = P.symmetrize_inplane_pair(8.69, 8.35, ortho)
        self.assertAlmostEqual(xo, 8.69, places=6)
        self.assertAlmostEqual(yo, 8.35, places=6)


class S71WiringTests(unittest.TestCase):
    SRC = (HERE / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py").read_text(encoding="utf-8")

    def _block(self, var):
        m = re.search(r"%s = r'''(.*?)\n'''" % var, self.SRC, re.S)
        self.assertIsNotNone(m, var)
        return m.group(1)

    def test_embedded_scripts_compile_and_order(self):
        be = self._block("_be_code")
        vac = self._block("_vac_code")
        compile(be, "_be_code", "exec")
        compile(vac, "_vac_code", "exec")
        # core：nspin 修正 < 对称化 < 读 h5 算带边
        i_fix = be.index('nspin_norm_fix.fix_h5("deformation.h5"')
        i_sym = be.index('dp_symmetrize.symmetrize_h5("deformation.h5")')
        i_load = be.index('load_deformation_potentials("deformation.h5")')
        self.assertLess(i_fix, i_sym)
        self.assertLess(i_sym, i_load)
        self.assertIn("symmetrize_inplane_pair(vac[\"xx\"], vac[\"yy\"], structure)", be)
        # vac：nspin 修正 < 对称化 < 硬闸门
        j_fix = vac.index("nspin_norm_fix.fix_h5(h5_out")
        j_sym = vac.index("dp_symmetrize.symmetrize_h5(h5_out)")
        j_gate = vac.index("_gate_bad = []")
        self.assertLess(j_fix, j_sym)
        self.assertLess(j_sym, j_gate)

    def test_gen_need_declares_module(self):
        import yaml
        for sk in ("ke-dft-cpu", "zt-dft-cpu"):
            txt = (HERE.parent / sk / "skill.yaml").read_text(encoding="utf-8")
            line = next(ln for ln in txt.splitlines() if "gen_step9b_deform_read.py" in ln)
            self.assertIn("dp_symmetrize.py", line, sk)
            self.assertIn("nspin_norm_fix.py", line, sk)
            yaml.safe_load(txt)


if __name__ == "__main__":
    unittest.main(verbosity=2)
