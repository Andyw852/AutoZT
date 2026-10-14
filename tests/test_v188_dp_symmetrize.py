# -*- coding: utf-8 -*-
"""[V188] 形变势对称化：AMSET 的 h5 是逐分量幅值 |D|（非对角 = 2·D_xy），不是带符号张量（纯 numpy，不要 pymatgen）。

  · 旧法（V121：把幅值当张量做 rᵀ|D|r）在含 m_x/m_y/m_z 的 Laue 群里把非对角幅值平均成 0 —— 这里复现；
  · 置换类操作 + |r| 对幅值精确：精确对称的数据不被改；
  · 带符号 + 回退：干净数据不改；xx 噪声被投影掉；"近简并且 ± 符号相反"（应变劈裂简并能级）回退幅值；
    非简并点上 AMSET 的曲率/残差底（A = max(|D|, ½|K|δ)）被中心差分去掉；Kramers 型简并（符号一致）照用带符号。
"""
import itertools
import sys
import unittest
from pathlib import Path

import numpy as np

SKILL = Path(__file__).resolve().parents[1] / "skill" / "ke-dft-cpu"
if str(SKILL) not in sys.path:
    sys.path.insert(0, str(SKILL))
import dp_symmetrize as P  # noqa: E402

B = 2 * np.pi * np.eye(3)          # 简单立方倒格（行 = b_i）：笛卡尔操作直接作用在分数坐标上


def oh_ops():
    """m-3m 的 48 个笛卡尔操作（全是带符号置换）。"""
    ops = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            r = np.zeros((3, 3))
            for i, (j, sg) in enumerate(zip(perm, signs)):
                r[i, j] = sg
            ops.append(r)
    return ops


def d6h_ops():
    """6/mmm（x 沿 a）：绕 z 的 6 次旋转 × 垂直镜面 × 反演，共 24 个；其中只有 8 个是带符号置换。"""
    ops = []
    mz = np.diag([1.0, 1.0, -1.0])
    for n in range(6):
        c, s = np.cos(n * np.pi / 3), np.sin(n * np.pi / 3)
        rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
        mir = rot @ np.diag([1.0, -1.0, 1.0])
        for r in (rot, mir):
            for g in (r, -r):
                if not any(np.allclose(g, x) for x in ops):
                    ops.append(g)
    ops2 = list(ops)
    for g in ops:
        h = mz @ g
        if not any(np.allclose(h, x) for x in ops2):
            ops2.append(h)
    return ops2


def star_field(seed_k, D0, ops, B=B):
    """seed k 的对称星 + 物理上对称的带符号张量场 D(r k) = r D0 rᵀ（单条带）。返回 (kpoints, D[1, nk, 3, 3])。"""
    Binv = np.linalg.inv(B)
    kc = np.asarray(seed_k, float) @ B
    ks, Ds = [], []
    for r in ops:
        kf = (r @ kc) @ Binv
        kf = kf - np.rint(kf)
        if any(np.allclose(kf, k, atol=1e-8) for k in ks):
            continue
        ks.append(kf)
        Ds.append(r @ D0 @ r.T)
    return np.array(ks), np.array(Ds)[None]


def amset_conv(D):
    """张量 -> AMSET 约定（非对角 ×2）。"""
    return P._from_tensor(D)


# L 型谷（PbTe/Ge）：D = Ξd I + Ξu n nᵀ，n ∥ [111]；放在简单立方的 (¼,¼,¼)（8 个星点，小群 3m 保持 D_L）
XI_D, XI_U = 2.0, 15.0
N111 = np.ones(3) / np.sqrt(3)
D_L = XI_D * np.eye(3) + XI_U * np.outer(N111, N111)


class PermutationOpsTests(unittest.TestCase):
    def test_oh_all_perm_d6h_only_eight(self):
        self.assertEqual(sum(P.is_signed_permutation(r) for r in oh_ops()), 48)
        d6h = d6h_ops()
        self.assertEqual(len(d6h), 24)
        self.assertEqual(sum(P.is_signed_permutation(r) for r in d6h), 8)


class OldBugTests(unittest.TestCase):
    def test_v121_zeroes_off_diagonal_magnitudes(self):
        ops = oh_ops()
        kp, D = star_field([0.25, 0.25, 0.25], D_L, ops)
        A = np.abs(amset_conv(D))                          # AMSET h5：幅值
        maps = P.kpoint_maps_B(kp, B, ops)
        old, _ = P.symmetrize_field(A, maps, ops)
        self.assertLess(np.abs(old[..., 0, 1]).max(), 1e-12)          # Ξu/3 的剪切耦合被抹成 0
        self.assertAlmostEqual(A[0, 0, 0, 1], 2 * XI_U / 3, places=10)


class AbsPermTests(unittest.TestCase):
    def test_exact_abs_data_unchanged_cubic(self):
        ops = oh_ops()
        for seed in ([0.25, 0.25, 0.25], [0.13, 0.31, 0.07]):
            D0 = D_L if seed[0] == 0.25 else np.array([[3.0, -1.2, 0.4], [-1.2, -2.0, 0.9], [0.4, 0.9, 1.5]])
            kp, D = star_field(seed, D0, ops)
            A = np.abs(amset_conv(D))
            maps = P.kpoint_maps_B(kp, B, ops)
            out, cnt = P.symmetrize_field_abs(A, maps, ops)
            self.assertLess(np.abs(out - A).max(), 1e-12)
            self.assertTrue(np.all(cnt == 48))

    def test_abs_perm_cannot_equalize_hexagonal_xx_yy(self):
        """六方：xx/yy 只由 3/6 次轴联系，幅值做不了 —— 这正是需要带符号文件的原因。"""
        ops = d6h_ops()
        kp, D = star_field([1 / 3, 1 / 3, 0.0], np.diag([5.0, 5.0, 1.0]), ops)
        A = np.abs(amset_conv(D))
        A[..., 0, 0] *= 1.04
        out, _ = P.symmetrize_field_abs(A, P.kpoint_maps_B(kp, B, ops), ops)
        self.assertGreater(np.abs(out[..., 0, 0] - out[..., 1, 1]).max(), 0.1)


class CombineTests(unittest.TestCase):
    def _case(self, ops, seed, D0, noise=0.0):
        kp, D = star_field(seed, D0, ops)
        S = -amset_conv(D)                                 # AMSET 的 diff 是 E0 − E_def：符号整体取反无所谓
        A = np.abs(S)
        if noise:
            A = A.copy(); S = S.copy()
            A[..., 0, 0] *= 1 + noise
            S[..., 0, 0] *= 1 + noise
        gap = np.full(A.shape[:2], 1.0)                    # 非简并
        return kp, A, S, gap

    def test_clean_signed_data_unchanged_hex_and_cubic(self):
        for ops, seed, D0 in ((d6h_ops(), [0.21, 0.07, 0.0], np.array([[4.0, 0.7, 0.0], [0.7, 2.5, 0.0], [0.0, 0.0, 1.0]])),
                              (oh_ops(), [0.25, 0.25, 0.25], D_L)):
            kp, A, S, gap = self._case(ops, seed, D0)
            out, info = P.combine(A, S, P.kpoint_maps_B(kp, B, ops), ops, gap, 0.005)
            self.assertLess(np.abs(out - A).max(), 1e-10)
            self.assertEqual(info["n_fallback_abs"], 0)

    def test_hexagonal_xx_noise_projected_at_K(self):
        ops = d6h_ops()
        kp, A, S, gap = self._case(ops, [1 / 3, 1 / 3, 0.0], np.diag([5.0, 5.0, 1.0]), noise=0.04)
        out, _ = P.combine(A, S, P.kpoint_maps_B(kp, B, ops), ops, gap, 0.005)
        self.assertLess(np.abs(out[..., 0, 0] - out[..., 1, 1]).max(), 1e-10)
        self.assertAlmostEqual(out[0, 0, 0, 0], 5.0 * 1.02, places=10)    # (5.2 + 5.0)/2
        self.assertLess(np.abs(out[..., 0, 1]).max(), 1e-10)

    def test_degenerate_split_falls_back_to_abs(self):
        """简并能级被应变劈裂：S ≈ 0、A = |a|，且与相邻带间隙为 0 -> 用 AMSET 的幅值。"""
        ops = oh_ops()
        kp, A, S, gap = self._case(ops, [0.25, 0.25, 0.25], D_L)
        S = np.zeros_like(S)                               # ±ε 给出相反符号，平均为 0
        gap = np.zeros_like(gap)
        out, info = P.combine(A, S, P.kpoint_maps_B(kp, B, ops), ops, gap, 0.005)
        self.assertLess(np.abs(out - A).max(), 1e-10)
        self.assertEqual(info["n_fallback_abs"], A.shape[1])

    def test_curvature_floor_removed_when_nondegenerate(self):
        """非简并、D_xy 真值 0：AMSET 的 |·| 再平均给出 ½|K|δ 的底，中心差分是 0 -> 用 0。"""
        ops = oh_ops()
        kp, A, S, gap = self._case(ops, [0.25, 0.25, 0.25], np.diag([4.0, 4.0, 4.0]))
        A = A.copy()
        A[..., 0, 1] = A[..., 1, 0] = 1.5                  # 虚假的底（真值 0）
        out, info = P.combine(A, S, P.kpoint_maps_B(kp, B, ops), ops, gap, 0.005)
        self.assertLess(np.abs(out[..., 0, 1]).max(), 1e-10)
        self.assertEqual(info["n_fallback_abs"], 0)

    def test_kramers_like_degeneracy_keeps_signed(self):
        """简并但应变不劈裂（SOC + 反演的 Kramers 对）：符号一致 -> 照用带符号。"""
        ops = oh_ops()
        kp, A, S, gap = self._case(ops, [0.25, 0.25, 0.25], D_L, noise=0.04)
        gap = np.zeros_like(gap)
        out, info = P.combine(A, S, P.kpoint_maps_B(kp, B, ops), ops, gap, 0.005)
        self.assertEqual(info["n_fallback_abs"], 0)
        self.assertLess(np.abs(out[..., 0, 1] - 2 * XI_U / 3).max(), 1e-10)   # 剪切耦合保住

    def test_no_gap_data_is_conservative(self):
        ops = oh_ops()
        kp, A, S, gap = self._case(ops, [0.25, 0.25, 0.25], D_L)
        S = np.zeros_like(S)
        out, info = P.combine(A, S, P.kpoint_maps_B(kp, B, ops), ops)
        self.assertEqual(info["fallback_rule"], "sign-opposite (no gap data)")
        self.assertLess(np.abs(out - A).max(), 1e-10)


class ConventionTests(unittest.TestCase):
    def test_tensor_roundtrip(self):
        X = np.arange(9.0).reshape(3, 3)
        X = X + X.T
        T = P._to_tensor(X)
        self.assertEqual(T[0, 1], X[0, 1] / 2)
        self.assertEqual(T[1, 1], X[1, 1])
        self.assertTrue(np.allclose(P._from_tensor(T), X))

    def test_signed_name(self):
        self.assertEqual(P.signed_name("a/deformation.h5"), "a/deformation_signed.h5")
        self.assertEqual(P.signed_name("deformation_vac.h5"), "deformation_vac_signed.h5")
        self.assertNotIn("deformation_potentials", P.SIGNED_KEY)



# ============================ [V189] 非对角换成张量分量 ============================
try:
    import h5py  # noqa: F401
    HAVE_H5 = True
except ImportError:                                        # pragma: no cover
    HAVE_H5 = False


class TensorContractionTests(unittest.TestCase):
    """AMSET 运行期用 tensordot(q̂⊗û, D)。物理：ΔE = Σ D_ij ε_ij，ε = sym(q̂⊗û)，D 对称 => Σ D_ij q_i u_j。
    所以 D 必须是张量分量；AMSET 原生 h5 的非对角是 2·D_xy，剪切部分被算大一倍。"""

    def test_contraction_needs_tensor_components(self):
        rng = np.random.default_rng(7)
        A = rng.normal(size=(3, 3))
        D = 0.5 * (A + A.T)                                  # 张量
        native = P._from_tensor(D)                         # AMSET h5 里的值（非对角 ×2）
        for _ in range(20):
            q = rng.normal(size=3); q /= np.linalg.norm(q)
            u = rng.normal(size=3); u -= u.dot(q) * q; u /= np.linalg.norm(u)
            eps = 0.5 * (np.outer(q, u) + np.outer(u, q))
            dE = float(np.sum(D * eps))                        # 物理能移
            self.assertAlmostEqual(float(np.tensordot(np.outer(q, u), P._to_tensor(native))), dE, places=12)
            # 原生值直接缩并：差的正好是非对角那部分
            off = float(np.tensordot(np.outer(q, u), native)) - dE
            self.assertAlmostEqual(off, float(np.sum((D * ~np.eye(3, dtype=bool)) * eps)), places=12)

    def test_l_valley_la_mode(self):
        """L 谷 LA[111]：物理 Ξd + Ξu；原生 h5 直接缩并得 Ξd + 5Ξu/3。"""
        q = N111
        native = P._from_tensor(D_L)
        self.assertAlmostEqual(float(np.tensordot(np.outer(q, q), D_L)), XI_D + XI_U, places=10)
        self.assertAlmostEqual(float(np.tensordot(np.outer(q, q), native)), XI_D + 5 * XI_U / 3, places=10)


@unittest.skipUnless(HAVE_H5, "需要 h5py")
class ToTensorH5Tests(unittest.TestCase):
    def _h5(self, A, attrs=None):
        import h5py
        import tempfile
        p = Path(tempfile.mkdtemp()) / "deformation.h5"
        with h5py.File(str(p), "w") as f:
            f.create_dataset("deformation_potentials_up", data=A)
            f["kpoints"] = np.zeros((A.shape[1], 3))
            for k, v in (attrs or {}).items():
                f.attrs[k] = v
        return p

    def test_inplace_halves_offdiag_once_and_keeps_attrs(self):
        import h5py
        rng = np.random.default_rng(0)
        A = np.abs(rng.normal(size=(2, 5, 3, 3)))
        A = 0.5 * (A + np.swapaxes(A, -1, -2))
        p = self._h5(A, {P.ATTR: 2, "nspin_norm_fixed": 1})
        self.assertFalse(P.is_tensor_convention(p))
        P.to_tensor_h5(p)
        self.assertTrue(P.is_tensor_convention(p))
        self.assertIn("skipped", P.to_tensor_h5(p))          # 幂等
        with h5py.File(str(p), "r") as f:
            B = np.array(f["deformation_potentials_up"])
            self.assertEqual(int(f.attrs["nspin_norm_fixed"]), 1)
            self.assertEqual(int(f.attrs[P.ATTR]), 2)
        off = ~np.eye(3, dtype=bool)
        self.assertTrue(np.allclose(B[..., off], 0.5 * A[..., off]))
        self.assertTrue(np.allclose(B[..., ~off], A[..., ~off]))

    def test_symmetrize_refuses_tensor_file(self):
        p = self._h5(np.ones((1, 2, 3, 3)), {P.OFFDIAG_ATTR: "tensor"})
        rep = P.symmetrize_h5(p)
        self.assertIn("to_tensor_h5 之前", rep["skipped"])


if __name__ == "__main__":
    unittest.main()
