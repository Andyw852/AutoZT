#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_dp_symmetrize.py —— V121 形变势场 D(k) 点群对称化的检验（不需要 VASP）。

用法：python test_dp_symmetrize.py      退出码 0 = 全部 PASS。

1) 约定（最关键）：用紧束缚模型对晶格施加 6 个应变分量做有限差分，得到**物理上**满足
   D(r k) = r D(k) rᵀ 的形变势场；dp_symmetrize 作用在它上面必须几乎不改值（< 1e-8），
   而故意写错的约定（张量用 rᵀ 旋转）必须明显改值 —— 2D MoS2、3D GaN 各一例。
2) 注入 4% 的 xx/yy 噪声（模拟 ε_xx、ε_yy 独立有限差分）：K 点（有 C3）投影回 xx=yy、xy=0；幂等。
3) h5：原地修改、保留已有 attrs（如 nspin_norm_fix 的标记）、打标记后重复调用跳过。
   [V188] h5 是 AMSET 的**幅值**（逐分量 |·|，非对角 = 2·D_xy）；带符号中心差分在 *_signed.h5 里；
   V121 旧法标记（attrs = 1）的 h5 不再动、提示重跑。
6) [V188] 真实 AMSET 函数端到端（需要 amset 环境）：紧束缚模型 + amset 的 get_strain_mapping /
   calculate_deformation_potentials；旧法把非对角幅值抹掉，新法在非简并点上与线性响应真值一致。
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
    def _write(self, d, A, kp, st, attrs=None):
        import h5py
        p = Path(d) / "deformation.h5"
        with h5py.File(str(p), "w") as f:                     # 与 amset write_deformation_potentials 同布局
            f.create_dataset("deformation_potentials_up", data=A)
            f["structure"] = np.bytes_(st.to_json())
            f["kpoints"] = kp
            for k, v in (attrs or {}).items():
                f.attrs[k] = v
        return p

    def _data(self):
        st = tsg.mos2_aligned()
        g = np.array(list(np.ndindex(6, 6, 1)), float) / np.array([6, 6, 1])
        kp = g - np.rint(g)
        rng = np.random.default_rng(0)
        S = rng.normal(size=(3, len(kp), 3, 3))
        S = 0.5 * (S + np.swapaxes(S, -1, -2))
        return st, kp, S

    def test_h5_inplace_keeps_attrs_and_is_idempotent(self):
        import h5py
        from pymatgen.electronic_structure.core import Spin
        st, kp, S = self._data()
        d = tempfile.mkdtemp()
        p = self._write(d, np.abs(S), kp, st, {"nspin_norm_fixed": 1})   # 模拟 nspin_norm_fix 的标记
        P.write_signed(P.signed_name(p), {Spin.up: S}, kp,
                       gaps={Spin.up: np.ones(S.shape[:2])}, delta=0.005)
        self.assertFalse(P.is_symmetrized(p))
        rep = P.symmetrize_h5(p)
        self.assertEqual(rep["k_missing_images"], 0)
        self.assertEqual(rep["signed"], "used")
        info = rep["datasets"]["deformation_potentials_up"]
        self.assertEqual(info["mode"], "signed")
        self.assertEqual(info["n_fallback_abs"], 0)
        self.assertGreater(info["max_rel_change"], 0.1)
        self.assertTrue(P.is_symmetrized(p))
        with h5py.File(str(p), "r") as f:
            self.assertEqual(int(f.attrs["nspin_norm_fixed"]), 1)
            first = np.array(f["deformation_potentials_up"])
        self.assertGreaterEqual(first.min(), 0.0)              # 写回的仍是幅值
        self.assertIn("skipped", P.symmetrize_h5(p))
        with h5py.File(str(p), "r") as f:
            self.assertTrue(np.array_equal(first, np.array(f["deformation_potentials_up"])))

    def test_without_signed_file_uses_permutation_subgroup_only(self):
        st, kp, S = self._data()
        p = self._write(tempfile.mkdtemp(), np.abs(S), kp, st)
        rep = P.symmetrize_h5(p)
        self.assertIn("没有带符号文件", rep["signed"])
        self.assertEqual(rep["datasets"]["deformation_potentials_up"]["mode"], "abs-perm")
        self.assertEqual(rep["n_perm_ops"], 8)                # 6/mmm 里只有 mmm 子群能精确作用在幅值上

    def test_legacy_v121_marked_h5_is_left_alone(self):
        import h5py
        st, kp, S = self._data()
        p = self._write(tempfile.mkdtemp(), np.abs(S), kp, st, {P.ATTR: 1})
        rep = P.symmetrize_h5(p)
        self.assertIn("V121", rep["skipped"])
        self.assertFalse(P.is_symmetrized(p))
        with h5py.File(str(p), "r") as f:
            self.assertTrue(np.array_equal(np.abs(S), np.array(f["deformation_potentials_up"])))

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


def _amset_ok():
    try:
        import spglib
        if not hasattr(spglib, "spglib"):     # amset 0.5.1 写的是 from spglib import spglib（旧 spglib 的子模块名）
            spglib.spglib = spglib
        from amset.deformation.potentials import calculate_deformation_potentials  # noqa: F401
        return _HAS_TB
    except Exception:                                       # noqa: BLE001
        return False


@unittest.skipUnless(_amset_ok(), "需要 amset + pymatgen（集群 AMSET 环境里跑）")
class AmsetEndToEndTests(unittest.TestCase):
    """真实 AMSET 的形变势链路 + 紧束缚模型（跳跃加平滑截断，0.5% 应变下哈密顿量连续）。
    真值 = 紧束缚线性响应（±1e-5 中心差分）换成 AMSET 约定的幅值；只在非简并点（间隙 > 应变能移）上比较。"""

    @staticmethod
    def _smooth_tb():
        import inspect
        src = inspect.getsource(TDF.tb_coeffs).replace(
            "terms.append((i, j, R, -amp * np.exp(-d / lam)))",
            "terms.append((i, j, R, -amp * np.exp(-d / lam) * (1 - (d / rcut) ** 2) ** 3))")
        assert "(d / rcut) ** 2" in src
        ns = dict(TDF.__dict__)
        exec(src, ns)
        return ns["tb_coeffs"]

    def _run(self, st, mesh):
        from pymatgen.analysis.elasticity.strain import Deformation
        from pymatgen.electronic_structure.bandstructure import BandStructure
        from pymatgen.electronic_structure.core import Spin
        from amset.deformation.generation import get_strains
        from amset.deformation.potentials import get_strain_mapping, calculate_deformation_potentials
        tb = self._smooth_tb()
        kp = TDF.full_mesh(mesh)
        g0 = np.zeros((1, 3), int)

        def bs(s):
            e = tb(s, kp, g0)[1]
            return BandStructure(kp, {Spin.up: e}, s.lattice.reciprocal_lattice, efermi=0.0, structure=s)
        bulk = {"reference": 0.0, "bandstructure": bs(st)}
        calcs = [{"reference": 0.0,
                  "bandstructure": bs(Deformation(e.get_deformation_matrix()).apply_to_structure(st))}
                 for e in get_strains(0.005)]
        sm = get_strain_mapping(st, calcs)
        A = calculate_deformation_potentials(bulk, sm)[Spin.up]
        S = P.signed_deformation_potentials(bulk, sm)[Spin.up]
        gap = P.band_gaps(bulk["bandstructure"])[Spin.up]
        # 真值：同一个平滑紧束缚模型的线性响应
        orig = TDF.tb_coeffs
        TDF.tb_coeffs = tb
        try:
            D, _ = tb_D(st, kp)
        finally:
            TDF.tb_coeffs = orig
        X = np.abs(P._from_tensor(D))
        rots = P.laue_rotations(st)
        maps = P.kpoint_maps(kp, st, rots)
        new, info = P.combine(A, S, maps, rots, gap, P.strain_amount(sm))
        old, _ = P.symmetrize_field(A, maps, rots)
        ok = gap > 2 * 0.005 * X.max(axis=(-1, -2)) + 0.05
        off = ~np.eye(3, dtype=bool)
        return X, A, new, old, ok, off

    def test_hexagonal_and_cubic(self):
        for name, st, mesh in (("GaN 6/mmm", tsg.gan_std(), [4, 4, 2]),
                               ("zincblende m-3m", tsg.zincblende_std(), [4, 4, 4])):
            with self.subTest(name=name):
                X, A, new, old, ok, off = self._run(st, mesh)
                self.assertGreater(ok.sum(), 0)
                scale = X[ok].max()
                self.assertLess(np.abs(new - X)[ok].max() / scale, 1e-3, "新法应与线性响应真值一致")
                # 旧法：非对角幅值被平均掉
                self.assertGreater(np.abs(old - X)[ok][:, off].max() / scale, 0.3)
                # AMSET 自己的 |·| 再平均在非简并的小分量上有曲率底，新法去掉了它
                self.assertGreaterEqual(np.abs(A - X)[ok].max(), np.abs(new - X)[ok].max())

    def test_v189_contraction_matches_true_band_shift(self):
        """[V189] 沿随机 q̂、横向 û 施加声子应变 sym(q̂⊗û)·η，直接算带边移动斜率；AMSET 的 tensordot(q̂⊗û, D)
        用原生 h5（非对角 2·D_xy）差到几十个百分点，换成张量分量后一致。"""
        from pymatgen.analysis.elasticity.strain import Deformation
        from pymatgen.core import Lattice, Structure
        from pymatgen.electronic_structure.bandstructure import BandStructure
        from pymatgen.electronic_structure.core import Spin
        from amset.deformation.generation import get_strains
        from amset.deformation.potentials import get_strain_mapping
        tb = self._smooth_tb()
        st, kp, g0 = tsg.gan_std(), TDF.full_mesh([4, 4, 2]), np.zeros((1, 3), int)

        def E(s):
            return tb(s, kp, g0)[1]

        def bs(s):
            return BandStructure(kp, {Spin.up: E(s)}, s.lattice.reciprocal_lattice, efermi=0.0, structure=s)
        bulk = {"reference": 0.0, "bandstructure": bs(st)}
        calcs = [{"reference": 0.0, "bandstructure": bs(Deformation(e.get_deformation_matrix()).apply_to_structure(st))}
                 for e in get_strains(0.005)]
        X = P.signed_deformation_potentials(bulk, get_strain_mapping(st, calcs))[Spin.up]   # 原生约定，-dE/dε
        rng = np.random.default_rng(3)
        q = rng.normal(size=3); q /= np.linalg.norm(q)
        u = rng.normal(size=3); u -= u.dot(q) * q; u /= np.linalg.norm(u)
        eta = 1e-4

        def strained(sign):
            eps = sign * eta * 0.5 * (np.outer(q, u) + np.outer(u, q))
            return Structure(Lattice(st.lattice.matrix @ (np.eye(3) + eps).T), st.species, st.frac_coords)
        slope = (E(strained(+1)) - E(strained(-1))) / (2 * eta)
        gap = P.band_gaps(bulk["bandstructure"])[Spin.up]
        ok = (gap > 0.05) & (np.abs(slope) > 0.1 * np.abs(slope).max())
        self.assertGreater(ok.sum(), 5)
        native = -np.einsum("i,j,bnij->bn", q, u, X)
        tensor = -np.einsum("i,j,bnij->bn", q, u, P._to_tensor(X))
        rel = lambda y: np.max(np.abs(y - slope)[ok] / np.abs(slope[ok]))      # noqa: E731
        self.assertLess(rel(tensor), 0.01)
        self.assertGreater(rel(native), 0.2)


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
        # [V189] 对称化之后把非对角换成张量分量（core 在读带边之前、vac 在硬闸门之前；两处都只读对角）
        i_ten = be.index('dp_symmetrize.to_tensor_h5("deformation.h5")')
        self.assertLess(i_sym, i_ten)
        self.assertLess(i_ten, i_load)
        j_ten = vac.index("dp_symmetrize.to_tensor_h5(h5_out)")
        self.assertLess(j_sym, j_ten)
        self.assertLess(j_ten, j_gate)
        # [V188] 带符号文件：vac 在写幅值 h5 的同一处写出，且在对称化之前
        j_sig = vac.index("dp_symmetrize.write_signed(dp_symmetrize.signed_name(h5_out)")
        self.assertLess(j_sig, j_sym)
        self.assertIn("dp_symmetrize.signed_deformation_potentials(bulk, sm)", vac)
        # core：amset deform read 的包装里算带符号版、写 deformation_signed.h5
        nb = self._block("_NBANDS_FIX_CODE")
        compile(nb, "_NBANDS_FIX_CODE", "exec")
        self.assertIn("_DPS.signed_deformation_potentials(bulk, sm)", nb)
        self.assertIn("_DPS.write_signed(_DPS.signed_name(fn)", nb)
        self.assertLess(nb.index("_PT.calculate_deformation_potentials, _PT.extract_bands"),
                        nb.index('cli.main(args=["deform", "read"]'))
        # 旧的带符号文件在本次读取之前删掉（没算出来时不能让对称化读到上一次的）
        self.assertLess(self.SRC.index('dfm / "deformation_signed.h5", out / "deformation_signed.h5"'),
                        self.SRC.index('print("[..] amset deform read ...")'))

    def test_gen_need_declares_module(self):
        import yaml
        checked = 0
        for sk in ("ke-dft-cpu", "zt-dft-cpu"):
            txt = (HERE.parent / sk / "skill.yaml").read_text(encoding="utf-8")
            line = next((ln for ln in txt.splitlines() if "gen_step9b_deform_read.py" in ln), None)
            if line is None:                # zt-dft-cpu 现在不再自带 S7.1（只组合 ke 的结果），不在它那儿查
                continue
            checked += 1
            self.assertIn("dp_symmetrize.py", line, sk)
            self.assertIn("nspin_norm_fix.py", line, sk)
            yaml.safe_load(txt)
        self.assertGreaterEqual(checked, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
