#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_mobility_vs_dpt.py —— V143：载流子符号的唯一真源 + 分机制迁移率对照工具 + DPT 失败原因。

用法：python test_mobility_vs_dpt.py      退出码 0 = 全部 PASS。不需要 AMSET / 数据。

1) ke_common.carrier_of_doping：负 = 电子、正 = 空穴（AMSET FermiDos 的约定）。
2) tools/mobility_vs_dpt.py：用 MoSe2 S8.4 的实测数（300 K，±4.28e17）回归 —— 电子行必须是负掺杂、
   ADP 4639.8、ADP/DPT = 28.5（标 ⚠ 多谷），空穴 0.74；这正是 2026-10-02 手写脚本标反的那张表。
   Seebeck 符号与载流子不符 -> ★、退出码 1。3D 取迹/3。
3) S8.2 全部失败时的报错带上 m*/C/E1 的具体原因（WS2 的"K 离网"曾被误读成 vasprun 缺失）。
"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "tools")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402
import mobility_vs_dpt as M  # noqa: E402


def _t(v, zz=1.0):
    return [[v, 0, 0], [0, v, 0], [0, 0, zz]]


def _transport(swap_seebeck=False):
    """MoSe2 S8.4 300 K 的实测值（V141 表）；100 K 一档放占位数。"""
    vals = {  # doping -> (S, overall, ADP, POP, IMP)
        -4.28e17: (-250.0, 35.1, 4639.8, 334.7, 86.4),
        4.28e17: (260.0, 44.1, 1050.6, 342.7, 84.6),
    }
    dop = list(vals)
    d = {"doping": dop, "temperatures": [100.0, 300.0], "seebeck": [], "conductivity": [],
         "mobility": {k: [] for k in ("overall", "ADP", "POP", "IMP")}}
    for x in dop:
        s, ov, adp, pop, imp = vals[x]
        if swap_seebeck:
            s = -s
        d["seebeck"].append([_t(1.0), _t(s)])
        d["conductivity"].append([_t(1.0), _t(1e4)])
        for k, v in zip(("overall", "ADP", "POP", "IMP"), (ov, adp, pop, imp)):
            d["mobility"][k].append([_t(1.0), _t(v, zz=0.0)])
    return d


def _material(swap=False, dpt=True):
    m = Path(tempfile.mkdtemp())
    (m / "step8.4_amset2d").mkdir()
    (m / "step8.4_amset2d" / "transport.json").write_text(json.dumps(_transport(swap)))
    if dpt:
        (m / "step8.2_dpt").mkdir()
        (m / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
            {"carrier": "electron", "mobility_cm2_Vs": 162.9},
            {"carrier": "hole", "mobility_cm2_Vs": 1416.7}]}))
    return m


class CarrierTests(unittest.TestCase):
    def test_sign(self):
        self.assertEqual(kc.carrier_of_doping(-1e17), "electron")
        self.assertEqual(kc.carrier_of_doping(1e17), "hole")
        self.assertIsNone(kc.carrier_of_doping(0))

    def test_mose2_regression(self):
        m = _material()
        rs, mechs = M.rows(json.loads((m / "step8.4_amset2d" / "transport.json").read_text()), 300, True,
                           M.dpt_by_carrier(m / "step8.2_dpt" / "dpt_result.json"))
        self.assertEqual(sorted(mechs), ["ADP", "IMP", "POP"])
        e = next(r for r in rs if r["carrier"] == "electron")
        h = next(r for r in rs if r["carrier"] == "hole")
        self.assertLess(e["doping"], 0)                                  # 电子 = 负掺杂
        self.assertEqual((e["ADP"], e["DPT"]), (4639.8, 162.9))
        self.assertAlmostEqual(e["ADP_over_DPT"], 28.48, places=2)
        self.assertAlmostEqual(h["ADP_over_DPT"], 0.742, places=3)
        self.assertTrue(e["sign_ok"] and h["sign_ok"])
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertEqual(M.main([str(m)]), 0)
        out = buf.getvalue()
        line_e = next(ln for ln in out.splitlines() if "电子(n)" in ln)
        self.assertIn("4639.8", line_e)
        self.assertIn("⚠", line_e)                                        # 多谷提示
        self.assertNotIn("⚠", next(ln for ln in out.splitlines() if "空穴(p)" in ln))

    def test_seebeck_sign_mismatch(self):
        m = _material(swap=True)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertEqual(M.main([str(m)]), 1)
        self.assertEqual(buf.getvalue().count("★"), 2)

    def test_3d_and_no_dpt(self):
        m = Path(tempfile.mkdtemp())
        (m / "step8_amset").mkdir()
        d = _transport()
        for k in d["mobility"]:
            d["mobility"][k] = [[_t(1.0), [[3.0, 0, 0], [0, 6.0, 0], [0, 0, 9.0]]] for _ in d["doping"]]
        (m / "step8_amset" / "transport.json").write_text(json.dumps(d))
        rs, _ = M.rows(d, 300, False, M.dpt_by_carrier(m / "nope.json"))
        self.assertEqual(rs[0]["ADP"], 6.0)                               # 迹/3
        self.assertIsNone(rs[0]["DPT"])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(M.main([str(m)]), 0)
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(M.main([tempfile.mkdtemp()]), 2)


class OtherToolsTests(unittest.TestCase):
    def test_rep_transport_labels(self):
        import subprocess
        m = _material()
        r = subprocess.run([sys.executable, str(ROOT / "tools" / "rep_transport.py"), "table",
                            str(m / "step8.4_amset2d" / "transport.json")], capture_output=True, text=True,
                           cwd=tempfile.mkdtemp())
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = r.stdout.splitlines()
        self.assertIn("载流子", lines[0])
        self.assertIn("电子(n)", next(ln for ln in lines if ln.startswith("-4.28")))
        self.assertIn("空穴(p)", next(ln for ln in lines if ln.startswith("4.28")))


class DegenerateEdgeTests(unittest.TestCase):
    """V145：GaAs（无 SOC）Γ 点价带顶三重简并 —— 原来取到轻空穴，DPT 空穴 276 万 cm²/Vs。"""

    def _fake(self):
        import types
        import numpy as np
        from pymatgen.electronic_structure.core import Spin
        k = np.linspace(0.0, 0.2, 5)                             # k[0] = Γ
        vb = np.stack([-40.0 * k ** 2, -15.0 * k ** 2, -3.0 * k ** 2], axis=1)   # 按能量升序：轻 < 中 < 重
        cb = (1.0 + 20.0 * k ** 2)[:, None]
        ene = np.concatenate([vb, cb], axis=1)
        occ = np.concatenate([np.ones_like(vb), np.zeros_like(cb)], axis=1)
        return types.SimpleNamespace(eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})

    def test_heavy_band_picked(self):
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        v = self._fake()
        isp, k0, b0, ene = D._find_band_edge(v, "hole")
        self.assertEqual((k0, b0), (0, 2))                        # 重空穴（编号最大），不是轻空穴（0）
        self.assertEqual(D._EDGE_DEGEN["hole"], 3)
        isp, k0, b0, ene = D._find_band_edge(v, "electron")
        self.assertEqual((k0, b0), (0, 3))
        self.assertEqual(D._EDGE_DEGEN["electron"], 1)
        src = (ROOT / "step8.2_dpt" / "gen_step12_dpt.py").read_text(encoding="utf-8")
        self.assertIn("重简并，取重带", src)

    def test_mobility_tool_respects_degenerate_dpt(self):
        m = _material()
        f = m / "step8.2_dpt" / "dpt_result.json"
        js = json.loads(f.read_text())
        js["results"][1]["mobility_cm2_Vs"] = 2.76e6                 # 空穴 DPT 离谱
        js["results"][1]["inputs"] = {"m_provenance": "...；带边 3 重简并，取重带（单带 DPT 对简并带边只是近似）"}
        f.write_text(json.dumps(js))
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            M.main([str(m)])
        line_h = next(ln for ln in buf.getvalue().splitlines() if "空穴(p)" in ln)
        self.assertIn("带边简并", line_h)
        self.assertNotIn("dp_valley_probe", line_h)


class DegeneracyTests(unittest.TestCase):
    """V146：GaAs 电子 m* 0.052 -> Nc ≈ 3e17；S8 的 3D 掺杂 1e18–1e21 都在简并区，ADP/DPT 超界不是多谷。"""

    def test_ratio(self):
        self.assertAlmostEqual(M.degeneracy_ratio(-1e17, 0.0523, 300, False), 1e17 / (2.509e19 * 0.0523 ** 1.5), 6)
        self.assertLess(M.degeneracy_ratio(-1e17, 0.0523, 300, False), 0.5)
        self.assertGreater(M.degeneracy_ratio(-1e18, 0.0523, 300, False), 3)
        # 2D：MoSe2 ±4.28e17 cm⁻³、c = 20 Å -> 8.6e10 cm⁻²，远低于 N2D（m* 0.5 -> 5.4e12）
        self.assertLess(M.degeneracy_ratio(4.28e17, 0.5, 300, True, 20.0), 0.05)
        self.assertIsNone(M.degeneracy_ratio(4.28e17, 0.5, 300, True, None))
        self.assertIsNone(M.degeneracy_ratio(1e17, None, 300, False))

    def test_gaas_like_flags(self):
        m = Path(tempfile.mkdtemp())
        (m / "step8_amset").mkdir()
        dop = [-1e17, -1e18, -1e21]
        d = {"doping": dop, "temperatures": [300.0], "conductivity": [[_t(1.0)] for _ in dop],
             "seebeck": [[_t(-300.0)] for _ in dop],
             "mobility": {"overall": [[_t(9000.0)] for _ in dop], "ADP": [[_t(1.0e6)] for _ in dop]}}
        (m / "step8_amset" / "transport.json").write_text(json.dumps(d))
        (m / "step8.2_dpt").mkdir()
        (m / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
            {"carrier": "electron", "mobility_cm2_Vs": 145658.0, "inputs": {"m_eff_m0": 0.0523}},
            {"carrier": "hole", "mobility_cm2_Vs": 4527.0, "inputs": {"m_eff_m0": 0.623}}]}))
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            M.main([str(m)])
        lines = {ln.split()[0]: ln for ln in buf.getvalue().splitlines() if "电子(n)" in ln}
        self.assertIn("dp_valley_probe", lines["-1e+17"])               # 非简并：多谷提示保留
        for k in ("-1e+18", "-1e+21"):
            self.assertIn("接近/进入简并", lines[k])
            self.assertNotIn("dp_valley_probe", lines[k])


class DptReasonTests(unittest.TestCase):
    def test_failure_message_has_provenance(self):
        s = (ROOT / "step8.2_dpt" / "gen_step12_dpt.py").read_text(encoding="utf-8")
        self.assertIn('miss.add("m*←S3_uniform（%s）" % i.get("m_provenance"))', s)
        self.assertIn('miss.add("E1←S7.1_read（%s）" % i.get("E1_provenance"))', s)
        self.assertNotIn('miss.add("m*←S3_uniform(vasprun/网格)")', s)


if __name__ == "__main__":
    unittest.main(verbosity=2)
