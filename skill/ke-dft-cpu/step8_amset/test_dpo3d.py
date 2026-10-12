#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_dpo3d.py —— [0036] 三维 ODP/IVS 的单元自检（纯 numpy，不需要 amset）。

    cd skill/ke-dft-cpu/step8_amset && python test_dpo3d.py

1. dp_optical 的 2D 口径逐字节不变：full3d=False 的谷归属与改动前（git HEAD 的 dp_optical.py）完全相同；
2. 三维谷归属用 kz：中心 (0,0,0) 与 (0,0,½) 两个谷能分开（2D 口径分不开）；
3. 三维抛物带：carrier_rates(full3d=True) 的谷内 ODP 率 = Conwell 解析式（analytic_rate_3d），发射阈值精确；
4. parse_valleys(full3d=True) 拒绝 2 分量的谷中心；merge_coincident 三维判同看 kz；
5. gen_step10 的 resolve_dp_optical_3d：2D 忽略、3D 缺谷表报错、正常写记录并把 ODP/IVS 补进 SCATTERING、
   显式 SCATTERING 照办、ADP 谷内掩码 auto 规则；_install_dp_optical_3d 写/删记录。
"""
import importlib.util
import io
import json
import math
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (str(ROOT / "step8.4_amset2d"), str(ROOT), str(ROOT.parent / "_common" / "opt"),
           str(ROOT / "step2_bandgap" / "step2.15_discriminant")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import dp_optical as dpo  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _head_dp_optical():
    """改动前的 dp_optical.py（git HEAD）；不在 git 仓库里 -> None（跳过逐字节对照）。"""
    try:
        src = subprocess.run(["git", "show", "HEAD:skill/ke-dft-cpu/step8.4_amset2d/dp_optical.py"],
                             cwd=str(ROOT), capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    if "full3d" in src:          # HEAD 已含本补丁：没有"改动前"可比
        return None
    d = tempfile.mkdtemp()
    p = Path(d) / "dp_optical_head.py"
    p.write_text(src, encoding="utf-8")
    return _load("dp_optical_head", p)


class TestDpOptical3D(unittest.TestCase):
    def test_2d_unchanged(self):
        old = _head_dp_optical()
        if old is None:
            self.skipTest("没有改动前的 dp_optical.py 可比")
        rng = np.random.default_rng(1)
        rlat = np.array([[1.0, 0.577, 0], [0, 1.155, 0], [0, 0, 0.16]])
        k = rng.random((4000, 3)) - 0.5
        centers = np.array([[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0]])
        self.assertTrue(np.array_equal(dpo.assign_members(k, centers, rlat), old.assign_members(k, centers, rlat)))
        E = rng.random((2, 4000)) * 0.02
        vs = [{"label": "K", "g": 2, "k_frac": centers[:2].tolist()}, {"label": "Γ", "g": 1, "k_frac": [[0, 0, 0]]}]
        chs = dpo.parse_channels("K>K:3@50; K>K':2@40; K>Γ:1@30")
        a = dpo.carrier_rates(E, k, rlat, vs, chs, [300.0], 1e-25, {"intra", "inter", "inter_star"})
        b = old.carrier_rates(E, k, rlat, vs, chs, [300.0], 1e-25, {"intra", "inter", "inter_star"})
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(dpo.valley_report(E, k, rlat, vs), old.valley_report(E, k, rlat, vs))

    def test_3d_assignment_uses_kz(self):
        rlat = np.eye(3)
        centers = np.array([[0, 0, 0], [0, 0, 0.5]])
        k = np.array([[0.05, 0, 0.45], [0.05, 0, 0.05], [0, 0, -0.45]])
        self.assertEqual(dpo.assign_members(k, centers, rlat, full3d=True).tolist(), [1, 0, 1])
        self.assertEqual(dpo.assign_members(k, centers, rlat).tolist(), [0, 0, 0])   # 2D 口径：kz 被丢掉
        self.assertTrue(dpo.same_member([0, 0, 0.5], np.array([[0, 0, -0.48]]), centers, rlat, full3d=True)[0])

    def test_3d_parabolic_vs_conwell(self):
        a = 5.43e-10                                     # 立方原胞边长（m）
        n = 64
        g1 = (np.arange(n) - n // 2) / n
        kx, ky, kz = np.meshgrid(g1, g1, g1, indexing="ij")
        kf = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], 1)
        rlat = np.eye(3) * (2 * math.pi / a)
        kc = kf @ rlat
        m = 0.5
        HBAR, ME, HA = 1.054571817e-34, 9.1093837015e-31, dpo.HA_J
        E = (HBAR ** 2 * np.sum(kc ** 2, 1) / (2 * m * ME)) / HA        # Hartree
        vs = [{"label": "Γ", "g": 1, "k_frac": [[0, 0, 0]]}]
        chs = dpo.parse_channels("Γ>Γ:5@60")
        M = 28.0855 * 2 * dpo.AMU
        r = dpo.carrier_rates(E[None], kf, rlat, vs, chs, [300.0], M, {"intra"}, +1, full3d=True)[0, 0]
        V = a ** 3
        for x_eV in (0.02, 0.1, 0.2):
            sel = np.abs(E * dpo.HA_EV - x_eV) < 0.004
            num = float(np.median(r[sel]))
            ref = dpo.analytic_rate_3d(x_eV * dpo.E_C, 0.0, m, 1, V, 5.0, 60.0, 300.0, M)
            self.assertLess(abs(num / ref - 1), 0.06, "x=%.2f eV：数值 %.4g vs 解析 %.4g" % (x_eV, num, ref))
        # 发射阈值：x < ħω 时只有吸收项
        x = 0.03
        ref = dpo.analytic_rate_3d(x * dpo.E_C, 0.0, m, 1, V, 5.0, 60.0, 300.0, M)
        N = 1 / math.expm1(0.06 * dpo.E_C / (dpo.KB * 300))
        absorb_only = math.pi * (5e10 * dpo.E_C) ** 2 / (M * 0.06 * dpo.E_C / HBAR) * N \
            * (V * (2 * m * ME) ** 1.5 / (4 * math.pi ** 2 * HBAR ** 3)) * math.sqrt(0.09 * dpo.E_C)
        self.assertAlmostEqual(ref / absorb_only, 1.0, places=10)

    def test_parse_and_merge_3d(self):
        with self.assertRaises(ValueError):
            dpo.parse_valleys("L:0.5,0.5", full3d=True)
        v = dpo.parse_valleys("L:0.5,0.5,0.5|0,0,0.5; X:0.5,0,0.5", full3d=True)
        self.assertEqual([x["g"] for x in v], [2, 1])
        vs = [{"label": "A", "g": 1, "k_frac": [[0, 0, 0]]}, {"label": "B", "g": 1, "k_frac": [[0, 0, 0.5]]}]
        out, alias = dpo.merge_coincident(vs, full3d=True)
        self.assertEqual((len(out), alias), (2, {}))
        out2, alias2 = dpo.merge_coincident(vs)                 # 2D 口径：同一面内位置 -> 合并
        self.assertEqual((len(out2), alias2), (1, {"B": "A"}))


class TestGenResolve(unittest.TestCase):
    def setUp(self):
        self.g = _load("g10_dpo3d", HERE / "gen_step10_amset.py")
        self.g.SCATTERING = ["ADP", "IMP", "POP"]
        self.g._SCATTERING_EXPLICIT = False
        for k in ("DP_OPTICAL_E", "DP_OPTICAL_H", "DP_OPTICAL_VALLEYS_E", "DP_OPTICAL_VALLEYS_H"):
            setattr(self.g, k, "")
        self.g.ADP_INTRAVALLEY = "auto"

    def _run(self, is_2d=False):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rec = self.g.resolve_dp_optical_3d(is_2d)
        return rec, buf.getvalue()

    def test_off_by_default(self):
        rec, _ = self._run()
        self.assertIsNone(rec)
        self.assertEqual(self.g.SCATTERING, ["ADP", "IMP", "POP"])

    def test_spec_keys(self):
        for k in ("DP_OPTICAL_E", "DP_OPTICAL_H", "DP_OPTICAL_VALLEYS_E", "DP_OPTICAL_VALLEYS_H", "ADP_INTRAVALLEY"):
            self.assertIn(k, self.g.SPEC)
        self.assertTrue({"ODP", "IVS"} <= self.g._ALLOWED_SCATTERING)

    def test_2d_ignored(self):
        self.g.DP_OPTICAL_H = "K>K:3@50"
        self.g.SCATTERING = ["ADP", "ODP"]
        rec, log = self._run(is_2d=True)
        self.assertIsNone(rec)
        self.assertEqual(self.g.SCATTERING, ["ADP"])
        self.assertIn("S8.4", log)

    def test_missing_valleys(self):
        self.g.DP_OPTICAL_H = "G>G:6@63"
        with self.assertRaises(SystemExit) as cm:
            self._run()
        self.assertIn("DP_OPTICAL_VALLEYS_H", str(cm.exception))

    def test_two_component_valley_rejected(self):
        self.g.DP_OPTICAL_H = "G>G:6@63"
        self.g.DP_OPTICAL_VALLEYS_H = "G:0,0"
        with self.assertRaises(SystemExit):
            self._run()

    def test_record_and_scattering(self):
        self.g.DP_OPTICAL_H = "G>G:6@63"
        self.g.DP_OPTICAL_VALLEYS_H = "G:0,0,0"
        self.g.DP_OPTICAL_E = "X>X:1@60; X>X':0.8@19; X>X':1.1@47"
        self.g.DP_OPTICAL_VALLEYS_E = "X:0.425,0,0.425|-0.425,0,-0.425|0,0.425,0.425|0,-0.425,-0.425|0.425,0.425,0|-0.425,-0.425,0"
        rec, log = self._run()
        self.assertEqual(rec["mechanisms"], ["ODP", "IVS"])
        self.assertEqual(self.g.SCATTERING, ["ADP", "IMP", "POP", "ODP", "IVS"])
        d = rec["dp_optical"]
        self.assertEqual(d["adp_intravalley"], {"electron": True, "hole": False})   # auto：只有开 IVS 的载流子
        self.assertEqual(d["electron"]["valleys"][0]["g"], 6)
        self.assertEqual(d["geometry"], "3d")
        json.dumps(rec, ensure_ascii=False)
        self.assertIn("ODP/IVS 开（3D）", log)

    def test_explicit_scattering(self):
        self.g.DP_OPTICAL_H = "G>G:6@63"
        self.g.DP_OPTICAL_VALLEYS_H = "G:0,0,0"
        self.g.SCATTERING, self.g._SCATTERING_EXPLICIT = ["ADP"], True
        rec, log = self._run()
        self.assertIsNone(rec)
        self.assertEqual(self.g.SCATTERING, ["ADP"])
        self.assertIn("不生效", log)

    def test_scattering_without_channels(self):
        self.g.SCATTERING = ["ADP", "IVS"]
        self.g.DP_OPTICAL_H = "G>G:6@63"
        self.g.DP_OPTICAL_VALLEYS_H = "G:0,0,0"
        with self.assertRaises(SystemExit):
            self._run()

    def test_install_writes_and_cleans(self):
        self.g.DP_OPTICAL_H = "G>G:6@63"
        self.g.DP_OPTICAL_VALLEYS_H = "G:0,0,0"
        rec, _ = self._run()
        with tempfile.TemporaryDirectory() as d:
            with redirect_stdout(io.StringIO()):
                self.assertTrue(self.g._install_dp_optical_3d(d, rec))
            for f in ("dp_optical_3d.json", "dp_optical.py", "amset_dpo3d.py"):
                self.assertTrue((Path(d) / f).is_file(), f)
            with redirect_stdout(io.StringIO()):
                self.assertFalse(self.g._install_dp_optical_3d(d, None))
            self.assertFalse((Path(d) / "dp_optical_3d.json").exists())

    def test_job_command_imports_plugin_last(self):
        src = (HERE / "gen_step10_amset.py").read_text(encoding="utf-8")
        i = src.index('"import amset_desym_fix; " if _has_fix')
        j = src.index('"import amset_dpo3d; " if _has_dpo3d')
        self.assertLess(i, j)
        y = (ROOT / "skill.yaml").read_text(encoding="utf-8")
        line = next(ln for ln in y.splitlines() if "gen_step10_amset.py" in ln and "gen_need" in ln)
        self.assertIn("dp_optical.py", line)
        self.assertIn("amset_dpo3d.py", line)


if __name__ == "__main__":
    unittest.main(verbosity=2)
