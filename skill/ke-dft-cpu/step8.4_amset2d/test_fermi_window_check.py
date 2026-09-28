#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_fermi_window_check.py —— fermi_window_check.py 的自检（离线，不需要 pymatgen）。

用法：python test_fermi_window_check.py     退出码 0 = 全部 PASS。
覆盖：fermi_levels 嵌套摊平 / scissor 的 Δ 计算 / 窗口边界 / 真实 GaAs S8 数值不误报 /
      -153.9 eV 必须失败 / 金属与缺产物必须跳过。
"""
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fermi_window_check as fwc

# 实测值（2026-09-28，GaAs S8 运行目录）：
#   vasprun 带边 VBM=3.7544 CBM=4.1720（DFT 带隙 0.4176）
#   settings.yaml bandgap=1.2735 -> Δ = 0.8559
#   transport.json 90 个费米能级，范围 3.043067 ~ 5.531336
GAAS_VBM, GAAS_CBM, GAAS_BG = 3.7544, 4.1720, 1.2735
GAAS_FERMI_MIN, GAAS_FERMI_MAX = 3.0430667867821493, 5.531336027855483


class FermiWindowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _write(self, name, obj):
        p = self.dir / name
        p.write_text(json.dumps(obj), encoding="utf-8")
        return p

    def _run(self, payload, vbm=-3.85, cbm=-2.03, argv=(), settings=None, delta=None):
        tf = self._write("transport_263x263x21.json", payload) if payload is not None else None
        vp = self.dir / "vasprun.xml"
        vp.write_text("<xml/>", encoding="utf-8")
        sp = self.dir / "settings.yaml"
        if settings is not None:
            sp.write_text(settings, encoding="utf-8")
        fwc.band_edges = lambda _p: (vbm, cbm)
        a = ["--file", str(tf)] if tf is not None else ["--file", str(self.dir / "nope.json")]
        a += ["--vasprun", str(vp), "--settings", str(sp)]
        if delta is not None:
            a += ["--delta", str(delta)]
        a += list(argv)
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = fwc.main(a)
        return rc, buf.getvalue()

    def test_flatten(self):
        p = self._write("transport_a.json", {"fermi_levels": [[-1.0], [-2.0]]})
        self.assertEqual(fwc.load_fermi_levels(p), [((0, 0), -1.0), ((1, 0), -2.0)])
        p = self._write("transport_b.json", {"fermi_levels": [-1.5, -2.5]})
        self.assertEqual(fwc.load_fermi_levels(p), [((0,), -1.5), ((1,), -2.5)])

    def test_scissor_delta(self):
        p = self.dir / "settings.yaml"
        p.write_text("bandgap: 1.2735\ninterpolation_factor: 10\n", encoding="utf-8")
        self.assertAlmostEqual(fwc.scissor_delta(p, 0.4176), 0.8559, places=4)
        p.write_text("scissor: 1.1\n", encoding="utf-8")
        self.assertAlmostEqual(fwc.scissor_delta(p, 0.4176), 1.1, places=6)
        p.write_text("scattering_type: [ADP]\n", encoding="utf-8")
        self.assertEqual(fwc.scissor_delta(p, 0.4176), 0.0)
        self.assertEqual(fwc.scissor_delta(self.dir / "nope.yaml", 0.5), 0.0)
        # bandgap 小于 DFT 带隙时不产生负位移
        p.write_text("bandgap: 0.1\n", encoding="utf-8")
        self.assertEqual(fwc.scissor_delta(p, 0.4176), 0.0)

    def test_window_boundaries(self):
        ok, lo, hi, bad = fwc.check_window([((0, 0), -1.85), ((1, 0), -3.61)], -3.85, -2.03, 2.0, 0.0)
        self.assertTrue(ok)
        self.assertAlmostEqual(lo, -5.85)
        self.assertAlmostEqual(hi, -0.03)
        ok, lo, hi, bad = fwc.check_window([((0, 0), 0.0)], -3.85, -2.03, 2.0, 0.0)
        self.assertFalse(ok)
        self.assertEqual(len(bad), 1)
        # Δ 抬高上限
        ok, lo, hi, _ = fwc.check_window([((0, 0), 0.5)], -3.85, -2.03, 2.0, 1.0)
        self.assertTrue(ok)
        self.assertAlmostEqual(hi, 0.97)

    def test_real_gaas_s8_passes(self):
        """真实 GaAs S8（±1e21 + scissor）不许误报 —— 用户 2026-09-28 明确要求实测。"""
        payload = {"fermi_levels": [[GAAS_FERMI_MAX], [GAAS_FERMI_MIN]]}
        rc, out = self._run(payload, vbm=GAAS_VBM, cbm=GAAS_CBM,
                            settings="bandgap: %s\n" % GAAS_BG)
        self.assertEqual(rc, 0, out)
        self.assertIn("[OK]", out)
        self.assertIn("0.8559", out)

    def test_scissor_prevents_false_positive(self):
        """Δ 不算是会误杀的：费米在 CBM+2.3，Δ=1.0 时应通过。"""
        payload = {"fermi_levels": [[GAAS_CBM + 2.3]]}
        rc, out = self._run(payload, vbm=GAAS_VBM, cbm=GAAS_CBM, delta=0.0)
        self.assertEqual(rc, 1, out)
        rc, out = self._run(payload, vbm=GAAS_VBM, cbm=GAAS_CBM, delta=1.0)
        self.assertEqual(rc, 0, out)

    def test_absurd_fermi_fails(self):
        rc, out = self._run({"fermi_levels": [[-1.8529], [-153.9]]}, settings="bandgap: 2.0\n")
        self.assertEqual(rc, 1, out)
        self.assertIn("[ERROR]", out)
        self.assertIn("-153.9000", out)

    def test_metal_skips(self):
        rc, out = self._run({"fermi_levels": [[-1.0]]}, vbm=0.0, cbm=0.0)
        self.assertEqual(rc, 0, out)
        self.assertIn("跳过", out)

    def test_missing_transport_skips(self):
        rc, out = self._run(None)
        self.assertEqual(rc, 0, out)
        self.assertIn("跳过", out)

    def test_no_fermi_levels_skips(self):
        rc, out = self._run({"doping": [-1e19]})
        self.assertEqual(rc, 0, out)
        self.assertIn("跳过", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
