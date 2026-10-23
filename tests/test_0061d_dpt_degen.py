""" [0061-d] DPT 简并带边 E1 修正的单元测试。

覆盖：
  · d=1（带边不简并）：_e1_degen_3d 返回 None —— 调用方逐位沿用原 E1（Si 电子 8.72 eV 路径不变）；
  · d=3（Γ 三重简并价带）：用合成 deformation_matrix.h5 + 弹性 OUTCAR，验证 E1_eff = E1_naive·√r、
    r>1（带间耦合让 E1 变大、DPT μ 变小），且溯源字段齐全（provenance=adp_degen、d、ratio、skill_rev）；
  · 手填优先：MANUAL["E1_hole_eV"] 直接返回手填值（不做修正）。
纯本地，需要 numpy + h5py。
"""
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(KE / "step8.2_dpt")):
    if _d not in sys.path:
        sys.path.insert(0, _d)

try:
    import h5py
    import numpy as np
    HAVE = True
except ImportError:                                   # pragma: no cover
    HAVE = False


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _D():
    import gen_step12_dpt as D
    return D


OUTCAR_ELASTIC = """ TOTAL ELASTIC MODULI (kBar)
  Direction      XX        YY        ZZ      XY      YZ      ZX
  XX          1657.0    639.0    639.0     0.0     0.0     0.0
  YY           639.0   1657.0    639.0     0.0     0.0     0.0
  ZZ           639.0    639.0   1657.0     0.0     0.0     0.0
  XY             0.0      0.0      0.0   796.0     0.0     0.0
  YZ             0.0      0.0      0.0     0.0   796.0     0.0
  ZX             0.0      0.0      0.0     0.0     0.0   796.0
"""

PAIRS = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]


def _write_matrix(path, kpoints, eig, xblocks):
    """合成 deformation_matrix.h5：xblocks[ip] 是 (nk, nb, nb) 的应变算符（Hermitian）。"""
    with h5py.File(str(path), "w") as f:
        f["X_up"] = np.asarray(xblocks, complex)
        f["eig_up"] = np.asarray(eig, float)
        f["kpoints"] = np.asarray(kpoints, float)
        f["pairs"] = np.asarray(PAIRS, np.int32)
        f["bands_up"] = np.arange(eig.shape[1], dtype=np.int64)
        f["amounts"] = np.asarray([[-0.005, 0.005]] * 6, float)
        f.attrs["ncl"] = 0


def _make_cwd(root, kpoints, eig, xblocks):
    (root / "step8_amset").mkdir(parents=True, exist_ok=True)
    (root / "step6_elastic").mkdir(parents=True, exist_ok=True)
    _write_matrix(root / "step8_amset" / "deformation_matrix.h5", kpoints, eig, xblocks)
    (root / "step6_elastic" / "OUTCAR").write_text(OUTCAR_ELASTIC)


def _simple_x(nk, nb, diag=(2.0, 2.0, 2.0), offdiag=0.0):
    """构造 (6, nk, nb, nb) 的应变算符：xx/yy/zz 对角 + xy/yz/zx 的带间耦合。"""
    out = np.zeros((6, nk, nb, nb), complex)
    for j in range(nk):
        for a in range(min(3, nb)):
            out[0][j, a, a] = diag[a]          # xx
            out[3][j, a, a] = diag[a]          # yy
            out[5][j, a, a] = diag[a]          # zz
        for a in range(nb):
            for b in range(nb):
                if a < b:
                    out[1][j, a, b] = out[1][j, b, a] = offdiag   # xy
                    out[2][j, a, b] = out[2][j, b, a] = offdiag   # xz
                    out[4][j, a, b] = out[4][j, b, a] = offdiag   # yz
    return out


@unittest.skipUnless(HAVE, "需要 numpy + h5py")
class TestDptDegenerateE1(unittest.TestCase):
    def test_d1_returns_none(self):
        D = _D()
        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            kp = np.array([[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]])
            eig = np.array([[0.0], [-0.5]])                # 单带，带边 d=1
            _make_cwd(cwd, kp, eig, _simple_x(2, 1))
            self.assertIsNone(D._e1_degen_3d(cwd, "hole", 3.47))

    def test_d3_hole_correction(self):
        D = _D()
        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            kp = np.array([[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]])
            # Γ 点 3 重简并（VBM），其它 k 更低；没有 band_edges.json -> 走兜底极值检测
            eig = np.array([[0.0, 0.0, 0.0], [-0.1, -0.2, -0.3], [-0.1, -0.2, -0.3]])
            _make_cwd(cwd, kp, eig, _simple_x(3, 3, diag=(2.0, 2.0, 2.0), offdiag=1.0))
            r = D._e1_degen_3d(cwd, "hole", 3.47)
            self.assertIsNotNone(r)
            info = r["info"]
            self.assertEqual(info["provenance"], "adp_degen")
            self.assertEqual(info["d"], 3)
            self.assertGreater(info["ratio"], 1.0)          # 带间耦合 => 额外散射
            self.assertAlmostEqual(info["E1_eff_eV"], round(3.47 * (info["ratio"]) ** 0.5, 4), places=3)
            self.assertGreater(info["E1_eff_eV"], 3.47)
            self.assertIn("adp_degen", r["prov"])
            self.assertIn(D._SKILL_REV, r["prov"])
            for k in ("subspace", "thermal_nk", "window_eV", "matrix_source", "fold"):
                self.assertIn(k, info)

    def test_manual_override_wins(self):
        D = _D()
        old = dict(D.MANUAL)
        try:
            D.MANUAL["E1_hole_eV"] = 5.42
            e1, prov = D.get_E1(Path("/nonexistent"), "hole")
            self.assertEqual(e1, 5.42)
            self.assertIn("manual", prov.lower())          # 手填标记（_manual_prov 输出）
        finally:
            D.MANUAL.clear()
            D.MANUAL.update(old)

    def test_missing_matrix_returns_none(self):
        D = _D()
        with tempfile.TemporaryDirectory() as td:
            self.assertIsNone(D._e1_degen_3d(Path(td), "hole", 3.47))


if __name__ == "__main__":
    unittest.main()
