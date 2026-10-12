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
        self.g.ADP_INTRAVALLEY, self.g.ODP_FROZEN = "auto", "off"     # 只测手填路径（冻声子路径见 TestGenFrozen3D）

    def _run(self, is_2d=False):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rec = self.g.resolve_dp_optical_3d(is_2d)
        return rec, buf.getvalue()

    def test_off_by_default(self):
        rec, _ = self._run()
        self.assertIsNone(rec["dp_optical"])
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
        self.g.ODP_FROZEN = "off"
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
        self.assertIsNone(rec["dp_optical"])
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


def _isym0_kpoints(N):
    """Γ 心 N³ 网格按 ±k 约化（VASP ISYM=0 的 k 列表）。"""
    seen, out = set(), []
    for i in range(N):
        for j in range(N):
            for l in range(N):
                key = (i, j, l)
                neg = tuple((-x) % N for x in key)
                if neg in seen:
                    continue
                seen.add(key)
                f = np.array(key, float) / N
                out.append(f - np.round(f))
    return np.array(out)


class TestValleys3D(unittest.TestCase):
    """[0036] S8.21 的三维自动谷表：类 Si 模型（导带 6 个 Δ 谷在 ±(3/8)·轴、价带 Γ 三重简并）。"""

    def setUp(self):
        import odp_frozen as of
        self.of = of
        N = 8
        self.k = _isym0_kpoints(N)
        kc = self.k * 2 * math.pi / 5.43
        cen = [s_ * 0.375 * np.eye(3)[a] for a in range(3) for s_ in (1, -1)]
        dmin = np.min([np.sum(((self.k - c) - np.round(self.k - c)) ** 2, 1) for c in cen], 0)
        cb = 1.0 + 8.0 * dmin
        cb2 = cb + 0.8
        vb = -np.sum(kc ** 2, 1) * 1.5
        e = np.stack([vb - 1.0, vb - 0.3, vb, vb, vb, cb, cb2], 1)          # 5 条占据带（后三条在 Γ 三重简并）
        self.eig = {"up": e}
        self.nocc = {"up": 5}

    def test_grid(self):
        N, off, table = self.of.grid_from_kpoints(self.k)
        self.assertEqual(N.tolist(), [8, 8, 8])
        self.assertEqual(len(table), 512)
        with self.assertRaises(ValueError):
            self.of.grid_from_kpoints(self.k[: len(self.k) // 2])

    def test_electron_six_valleys(self):
        vs = self.of.find_valleys_3d(self.k, self.eig, self.nocc, +1, np.eye(3))
        self.assertEqual(len(vs), 1)
        v = vs[0]
        self.assertEqual((v["label"], v["g"], v["band_index"]), ("V1", 6, 5))
        got = sorted(tuple(round(abs(x), 4) for x in k) for k in v["k_frac"])
        self.assertEqual(got, sorted([(0.375, 0, 0)] * 2 + [(0, 0.375, 0)] * 2 + [(0, 0, 0.375)] * 2))

    def test_hole_gamma_top_band(self):
        vs = self.of.find_valleys_3d(self.k, self.eig, self.nocc, -1, np.eye(3))
        self.assertEqual(len(vs), 1)
        self.assertEqual((vs[0]["label"], vs[0]["g"], vs[0]["band_index"]), ("Γ", 1, 4))

    def test_sig(self):
        a = [{"k_frac": [[0.5, 0, 0], [-0.5, 0, 0]]}]
        b = [{"k_frac": [[-0.5, 0, 0.0], [0.5, 0, 0]]}]
        self.assertEqual(self.of.valley_sig_3d(a), self.of.valley_sig_3d(b))
        self.assertNotEqual(self.of.valley_sig_3d(a), self.of.valley_sig_3d([{"k_frac": [[0.5, 0, 0.5]]}]))

    def test_bulk_polar_kind(self):
        import odp_frozen as of
        e = np.zeros((6, 2, 3))
        e[0, :, 0] = e[1, :, 1] = e[2, :, 2] = 1 / math.sqrt(2)
        e[3, 0, 0], e[3, 1, 0] = 1 / math.sqrt(2), -1 / math.sqrt(2)
        e[4, 0, 1], e[4, 1, 1] = 1 / math.sqrt(2), -1 / math.sqrt(2)
        e[5, 0, 2], e[5, 1, 2] = 1 / math.sqrt(2), -1 / math.sqrt(2)
        born = np.array([np.eye(3) * 2.0, -np.eye(3) * 2.0])
        modes, groups, _n = of.mode_table([0, 0, 0, 8, 8, 8], e, [30.0, 30.0], born, polar_z_tol=float("inf"))
        self.assertEqual({m["polar_kind"] for m in modes if m["polar"]}, {"体相"})
        self.assertFalse(groups[0]["use"])


def _frozen_params(of, valleys_hole, D=5.0, dim="3d", src="auto"):
    return {"schema": of.SCHEMA, "dim": dim,
            "carriers": {"hole": {"status": "ok", "valleys_source": src,
                                  "valley_sig": {"dim": "3d", "sig": of.valley_sig_3d(valleys_hole)},
                                  "valleys": [dict(v, groups=[{"group": 1, "modes": [4, 5, 6], "hw_meV": 63.0,
                                                               "D_eV_A": D}]) for v in valleys_hole]},
                         "electron": {"status": "ok", "valleys_source": src,
                                      "valley_sig": {"dim": "3d", "sig": [[[0.375, 0.0, 0.0]]]},
                                      "valleys": [{"label": "V1", "k_frac": [[0.375, 0, 0]],
                                                   "groups": [{"group": 1, "modes": [4], "hw_meV": 63.0,
                                                               "D_eV_A": 0.001}]}]}}}


class TestGenFrozen3D(unittest.TestCase):
    """[0036] S8 用 S8.21 的冻声子参数：出厂开 ODP、不开 IVS；手填谷内只作对照；谷表须同源。"""

    def setUp(self):
        import odp_frozen as of
        self.of = of
        self.g = _load("g10_frozen3d", HERE / "gen_step10_amset.py")
        self.g.SCATTERING = ["ADP", "IMP", "POP"]
        self.g._SCATTERING_EXPLICIT = False
        for k in ("DP_OPTICAL_E", "DP_OPTICAL_H", "DP_OPTICAL_VALLEYS_E", "DP_OPTICAL_VALLEYS_H"):
            setattr(self.g, k, "")
        self.g.ADP_INTRAVALLEY, self.g.ODP_FROZEN, self.g.ODP_FROZEN_DMIN = "auto", "auto", 0.05
        self.tmp = tempfile.TemporaryDirectory()
        self.cwd = Path(self.tmp.name)
        (self.cwd / "step8.21_odp_read").mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def _write(self, params):
        (self.cwd / "step8.21_odp_read" / "odp_params.json").write_text(json.dumps(params), encoding="utf-8")

    def _run(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rec = self.g.resolve_dp_optical_3d(False, self.cwd)
        return rec, buf.getvalue()

    def test_no_params_no_spec(self):
        rec, _ = self._run()
        self.assertIsNone(rec["dp_optical"])
        self.assertEqual(rec["dp_optical_frozen"]["frozen_use"], {"electron": "unavailable", "hole": "unavailable"})

    def test_default_frozen_odp_only(self):
        self._write(_frozen_params(self.of, [{"label": "Γ", "k_frac": [[0, 0, 0]]}]))
        rec, log = self._run()
        self.assertEqual(rec["mechanisms"], ["ODP"])
        self.assertEqual(self.g.SCATTERING, ["ADP", "IMP", "POP", "ODP"])
        h = rec["dp_optical"]["hole"]
        self.assertEqual([(c["from"], c["kind"], c["D_eV_A"], c["hw_meV"]) for c in h["channels"]],
                         [("Γ", "intra", 5.0, 63.0)])
        self.assertNotIn("electron", rec["dp_optical"])        # 电子 D < DMIN（对称性禁戒）-> 没有通道
        self.assertEqual(rec["dp_optical"]["adp_intravalley"], {"electron": False, "hole": False})
        self.assertEqual(rec["dp_optical_frozen"]["frozen_use"]["hole"], "used")
        self.assertTrue(rec["dp_optical_frozen"]["digest"]["hole"])
        self.assertIn("冻声子", log)

    def test_manual_intra_is_reference_only(self):
        self._write(_frozen_params(self.of, [{"label": "Γ", "k_frac": [[0, 0, 0]]}]))
        self.g.DP_OPTICAL_H = "G>G:9.9@60"
        rec, log = self._run()
        h = rec["dp_optical"]["hole"]
        self.assertEqual([c["D_eV_A"] for c in h["channels"]], [5.0])
        self.assertEqual(h["frozen"]["manual_reference"][0]["D_eV_A"], 9.9)
        self.assertIn("只作对照", log)

    def test_manual_ivs_with_frozen_valleys(self):
        vs = [{"label": "Γ", "k_frac": [[0, 0, 0]]}, {"label": "V1", "k_frac": [[0.5, 0.5, 0.5], [0.5, 0, 0]]}]
        self._write(_frozen_params(self.of, vs))
        self.g.DP_OPTICAL_H = "G>V1:1.5@40"
        rec, _ = self._run()
        self.assertEqual(rec["mechanisms"], ["ODP", "IVS"])
        self.assertEqual(rec["dp_optical"]["adp_intravalley"]["hole"], True)

    def test_off_uses_manual(self):
        self._write(_frozen_params(self.of, [{"label": "Γ", "k_frac": [[0, 0, 0]]}]))
        self.g.ODP_FROZEN = "off"
        self.g.DP_OPTICAL_H, self.g.DP_OPTICAL_VALLEYS_H = "G>G:9.9@60", "G:0,0,0"
        rec, _ = self._run()
        self.assertEqual([c["D_eV_A"] for c in rec["dp_optical"]["hole"]["channels"]], [9.9])
        self.assertEqual(rec["dp_optical_frozen"]["frozen_use"]["hole"], "off")

    def test_manual_valleys_must_match(self):
        self._write(_frozen_params(self.of, [{"label": "Γ", "k_frac": [[0, 0, 0]]}], src="step.conf"))
        self.g.DP_OPTICAL_VALLEYS_H = "G:0,0,0.5"
        with self.assertRaises(SystemExit) as cm:
            self._run()
        self.assertIn("S8.21", str(cm.exception))
        self.g.DP_OPTICAL_VALLEYS_H = "G:0,0,0"
        rec, _ = self._run()
        self.assertTrue(rec["dp_optical"]["hole"]["valleys_source"].startswith("manual"))

    def test_2d_params_ignored(self):
        self._write(_frozen_params(self.of, [{"label": "Γ", "k_frac": [[0, 0, 0]]}], dim=None))
        rec, log = self._run()
        self.assertIsNone(rec["dp_optical"])
        self.assertIn("不是三维口径", log)
        self.g.ODP_FROZEN = "on"
        with self.assertRaises(SystemExit):
            self._run()

    def test_explicit_scattering_excludes(self):
        self._write(_frozen_params(self.of, [{"label": "Γ", "k_frac": [[0, 0, 0]]}]))
        self.g.SCATTERING, self.g._SCATTERING_EXPLICIT = ["ADP", "IMP"], True
        rec, _ = self._run()
        self.assertIsNone(rec["dp_optical"])
        self.assertEqual(rec["dp_optical_frozen"]["frozen_use"]["hole"], "not_in_scattering")
        self.assertFalse(rec["dp_optical_frozen"]["odp_allowed"])

    def test_record_written_even_without_odp(self):
        rec, _ = self._run()
        d = self.cwd / "run"
        d.mkdir()
        with redirect_stdout(io.StringIO()):
            self.assertFalse(self.g._install_dp_optical_3d(d, rec))
        self.assertTrue((d / "dp_optical_3d.json").is_file())
        self.assertFalse((d / "amset_dpo3d.py").exists())

    def test_requeue_s8_when_params_appear(self):
        import ke_common as kc
        rec, _ = self._run()                                 # 参数出来之前生成的 S8
        d = self.cwd / "step8_amset"
        d.mkdir()
        with redirect_stdout(io.StringIO()):
            self.g._install_dp_optical_3d(d, rec)
        (d / "transport.json").write_text("{}", encoding="utf-8")
        params = _frozen_params(self.of, [{"label": "Γ", "k_frac": [[0, 0, 0]]}])
        with redirect_stdout(io.StringIO()):
            done = kc.invalidate_odp_consumer(self.cwd, params, consumer="step8_amset", record="dp_optical_3d.json")
        self.assertIn(("step8_amset", "transport.json"), done)
        self.assertFalse((d / "transport.json").exists())


class TestS821Read3D(unittest.TestCase):
    """[0036] S8.21 的三维读取（carriers_3d）：合成 ±u 构型，D 已知。
    导带 6 个 Δ 谷：带边随模 4 线性移动 c·u -> D = c；价带 Γ 三重简并劈裂 (+a, 0, −a)·u -> D² = 2a²/9。"""

    def test_known_D(self):
        from types import SimpleNamespace
        r = _load("g12c_3d", ROOT / "step7c_odp_frozen" / "odp_read" / "gen_step12c_odp_read.py")
        import odp_frozen as of
        r.of = of
        base = TestValleys3D("test_grid")
        base.setUp()
        k, e0 = base.k, base.eig["up"]
        u, c, a = 0.004, 3.0, 4.5
        gam = np.all(np.abs(k) < 1e-9, axis=1)

        def conf(name, sg):
            e = e0.copy()
            e[:, 5] += sg * c * u
            e[gam, 2] += sg * a * u
            e[gam, 4] -= sg * a * u
            e = np.sort(e, axis=1)
            cf = object.__new__(r.Config)
            cf.name, cf.kpts, cf.eig = name, k, {"up": e}
            cf.occ = {"up": (np.arange(e.shape[1])[None, :] < 5).astype(float) * np.ones((len(k), 1))}
            cf.structure = SimpleNamespace(lattice=SimpleNamespace(
                reciprocal_lattice=SimpleNamespace(matrix=np.eye(3) * 2 * math.pi / 5.43)))
            cf.vac = cf.vac_std = cf.vac_slope = 0.0
            return cf
        cfg = {"odp-ref": conf("odp-ref", 0), "odp-m04p": conf("p", 1), "odp-m04m": conf("m", -1)}
        man = {"selected": [{"index": 4}], "groups": [{"group": 1, "modes": [4], "hw_meV": 63.0, "freq_THz": 15.2,
                                                       "use": True, "polar": False}],
               "configs": {"odp-ref": {"mode": None, "sign": 0, "u_A": 0.0},
                           "odp-m04p": {"mode": 4, "sign": 1, "u_A": u}, "odp-m04m": {"mode": 4, "sign": -1, "u_A": u}}}
        warns = []
        with tempfile.TemporaryDirectory() as d, redirect_stdout(io.StringIO()):
            out = r.carriers_3d(Path(d), cfg, man, None, warns)
        e = out["electron"]
        self.assertEqual((e["status"], len(e["valleys"]), e["valleys"][0]["g"]), ("ok", 1, 6))
        self.assertAlmostEqual(e["valleys"][0]["groups"][0]["D_eV_A"], c, places=4)
        self.assertEqual(len(e["valleys"][0]["k_frac"]), 6)
        h = out["hole"]["valleys"][0]
        self.assertEqual((h["label"], h["band_index"]), ("Γ", 4))
        self.assertAlmostEqual(h["groups"][0]["D_eV_A"], math.sqrt(2 * a * a / 9), places=4)
        self.assertEqual(out["hole"]["valley_sig"]["sig"], of.valley_sig_3d([{"k_frac": [[0, 0, 0]]}]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
