"""[0028] 非磁体系用 ISPIN=2 算时，多谷 DPT 不能把同一个谷按 up/down 各算一次。

CrS2/CrSe2 hex 实测：谷表里 K 与 K(down)、M 与 M(down) 的 ΔE 完全相同，单谷占比报"K 谷只占 50%"，
另一半其实是同一个 K 谷的另一自旋。两道能带逐点相同就只留一道：结果必须和 ISPIN=1 完全一样；
两道能带不同（自旋劈裂）时照旧分开算。复用 test_0022 的 WS2 型模型。
"""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests import test_0022_multi_valley_dpt as T22

HAVE = T22.HAVE


class SpinModel(T22.Model):
    """同一个模型，eigenvalues 给两个自旋道；split≠0 时 down 道整体上移 split（eV）。
    top_split 不为 None 时再加一条离导带底 5 eV 的最高空带，down 道这条带另外上移 top_split（[0029]：
    VASP 最高几条空带收敛最松，CrS2/CrSe2 hex 实测两道差 0.16–0.38 eV，其余带 < 1 meV）。"""

    def __init__(self, split=0.0, top_split=None, **kw):
        super().__init__(**kw)
        self.split, self.top_split = split, top_split

    def vasprun(self, kf, strain=None, nocc=2):
        import numpy as np
        from pymatgen.electronic_structure.core import Spin
        v = super().vasprun(kf, strain, nocc)
        up = v.eigenvalues[Spin.up]
        if self.top_split is not None:
            top = up[:, 2:3, :].copy()
            top[..., 0] += 5.0
            up = np.concatenate([up, top], axis=1)
        dn = up.copy()
        dn[..., 0] += self.split
        if self.top_split is not None:
            dn[:, -1, 0] += self.top_split
        v.parameters = {"ISPIN": 2}
        v.eigenvalues = {Spin.up: up, Spin.down: dn}
        return v


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class SpinDupTests(unittest.TestCase):
    def setUp(self):
        self.D = T22._D()
        self.D._MV_CACHE.clear()
        self._p = [mock.patch.dict(self.D.MANUAL_ANISO, {"C_2D_xy_N_per_m": (T22.C2D, T22.C2D)}),
                   mock.patch.object(self.D, "VALLEY_WINDOW_EV", 0.30),
                   mock.patch.object(self.D, "DPT_VALLEYS", "auto"),
                   mock.patch.object(self.D, "E1_SOURCE", "vac")]
        for p in self._p:
            p.start()

    def tearDown(self):
        for p in self._p:
            p.stop()

    def _mv(self, model, carrier="hole"):
        self.D._MV_CACHE.clear()
        m = Path(tempfile.mkdtemp())
        vr = T22.build(model, m)
        (m / self.D.UNIFORM_DIR / "OUTCAR").write_text(
            " number of electron      26.0000000 magnetization       0.0000000\n")
        with mock.patch("pymatgen.io.vasp.Vasprun", vr), contextlib.redirect_stdout(io.StringIO()):
            bd = self.D._aniso_block(m, True, carrier, 300.0)
            mv = self.D._multi_valley_block(m, True, carrier, 300.0, [300.0, 600.0], bd)
        return mv

    def test_identical_spins_same_as_ispin1(self):
        one = self._mv(T22.Model())
        two = self._mv(SpinModel())
        self.assertEqual(one["status"], "ok", one["status"])
        self.assertEqual(two["status"], "ok", two["status"])
        self.assertEqual([(v["label"], v["g"]) for v in two["valleys"]],
                         [(v["label"], v["g"]) for v in one["valleys"]])
        self.assertFalse(any("down" in v["label"] for v in two["valleys"]))
        for r1, r2 in zip(one["mobility_vs_T"], two["mobility_vs_T"]):
            self.assertAlmostEqual(r2["mobility_x_cm2_Vs"], r1["mobility_x_cm2_Vs"], places=6)
            for k, w in r1["weights"].items():
                self.assertAlmostEqual(r2["weights"][k], w, places=9)
        s1, s2 = self.D._single_valley_scope(one, 300.0), self.D._single_valley_scope(two, 300.0)
        self.assertAlmostEqual(s2["share"], s1["share"], places=9)
        self.assertGreater(s2["share"], 0.5)                 # 以前会被对半分成 < 0.5

    def test_top_band_mismatch_still_merged(self):
        """[0029] 只有远离带边的最高空带两道不一致：仍然合并，结果与 ISPIN=1 相同。"""
        for car in ("hole", "electron"):
            one = self._mv(T22.Model(), car)
            two = self._mv(SpinModel(top_split=0.3), car)
            self.assertEqual(two["status"], one["status"], car)
            self.assertEqual([(v["label"], v["g"]) for v in two.get("valleys") or []],
                             [(v["label"], v["g"]) for v in one.get("valleys") or []], car)

    def test_split_spins_kept_apart(self):
        two = self._mv(SpinModel(split=0.05))
        labels = [v["label"] for v in two.get("valleys") or []]
        self.assertTrue(any("down" in x for x in labels), (two["status"], labels))


if __name__ == "__main__":
    unittest.main()
