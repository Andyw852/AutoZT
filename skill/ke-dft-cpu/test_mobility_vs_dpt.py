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


class GridResolutionTests(unittest.TestCase):
    """V147：GaAs（PBEsol 带隙 0.418 eV）的 Γ 谷很轻且非抛物，S3 最近网格点高出 CBM ~0.25 eV（~10 kT），
    DPT 拟合 0.0523，而 amset eff-mass 给带边 0.030。"""

    def _fake_vasprun(self, n=19, a=5.65, m_c=0.03, gap=0.418, m_v=0.6):
        import types
        import numpy as np
        from pymatgen.electronic_structure.core import Spin
        g = np.arange(-4, 5) / float(n)
        kf = np.array([[x, y, z] for x in g for y in g for z in g])
        recip = np.eye(3) * (2 * np.pi / a)
        q2 = np.sum((kf @ recip) ** 2, axis=1)
        alpha = 1.0 / gap
        rhs = 3.80998 * q2 / m_c
        ec = gap + (-1.0 + np.sqrt(1.0 + 4.0 * alpha * rhs)) / (2.0 * alpha)    # Kane 非抛物导带
        ev = -3.80998 * q2 / m_v
        ene = np.stack([ev, ec], axis=1)
        occ = np.stack([np.ones_like(ev), np.zeros_like(ec)], axis=1)
        lat = types.SimpleNamespace(reciprocal_lattice=types.SimpleNamespace(matrix=recip))
        return types.SimpleNamespace(final_structure=types.SimpleNamespace(lattice=lat),
                                     actual_kpoints=kf.tolist(),
                                     eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})

    def test_flagged_for_light_nonparabolic_band(self):
        from unittest import mock
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        m = Path(tempfile.mkdtemp())
        (m / D.UNIFORM_DIR).mkdir(parents=True)
        (m / D.UNIFORM_DIR / "vasprun.xml").write_text("<x/>")
        fake = self._fake_vasprun()
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake):
            for car in ("electron", "hole"):
                D._EDGE_NN_DE.pop(car, None)
                me, prov = D.get_effective_mass(m, car, False)
                self.assertIsNotNone(me, prov)
                with contextlib.redirect_stdout(io.StringIO()) as buf:
                    note = D._grid_resolution_note(car, me, prov, 300.0)
                if car == "electron":
                    self.assertGreater(D._EDGE_NN_DE[car], 0.2)          # ~0.27 eV，约 10 kT
                    self.assertIn("网格分辨不出带边曲率", note)
                    self.assertIn("amset eff-mass", buf.getvalue())
                else:
                    self.assertLess(D._EDGE_NN_DE[car], 0.04)            # 重空穴 ~0.02 eV < 4 kT
                    self.assertEqual(note, prov)
        self.assertEqual(D._grid_resolution_note("electron", 0.03, "manual", 300.0), "manual")

    def test_mobility_tool_respects_coarse_grid(self):
        m = _material()
        f = m / "step8.2_dpt" / "dpt_result.json"
        js = json.loads(f.read_text())
        js["results"][0]["inputs"] = {"m_provenance": "...；网格分辨不出带边曲率（最近网格点高出带边 266 meV = 10.3 kT）"}
        f.write_text(json.dumps(js))
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            M.main([str(m)])
        line_e = next(ln for ln in buf.getvalue().splitlines() if "电子(n)" in ln)
        self.assertIn("网格分辨不出带边曲率", line_e)
        self.assertNotIn("dp_valley_probe", line_e)


class DptReasonTests(unittest.TestCase):
    def test_failure_message_has_provenance(self):
        s = (ROOT / "step8.2_dpt" / "gen_step12_dpt.py").read_text(encoding="utf-8")
        self.assertIn('miss.add("m*←S3_uniform（%s）" % i.get("m_provenance"))', s)
        self.assertIn('miss.add("E1←S7.1_read（%s）" % i.get("E1_provenance"))', s)
        self.assertNotIn('miss.add("m*←S3_uniform(vasprun/网格)")', s)



class ExcBriefTests(unittest.TestCase):
    """V150：兜底 except 写出类型 + 信息 + 本脚本里出错的行（WS2 只看到"二次型拟合异常：IndexError"）。"""

    def test_brief_has_line(self):
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        try:
            D.mobility_dpt(True, None, 1.0, 1.0, 300.0)
        except Exception as e:                                           # noqa: BLE001
            msg = D._exc_brief(e)
        self.assertIn("TypeError", msg)
        self.assertRegex(msg, r"gen_step12_dpt\.py:\d+ `")

    def test_aniso_failure_reports_where(self):
        from unittest import mock
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        m = Path(tempfile.mkdtemp())
        (m / D.UNIFORM_DIR).mkdir(parents=True)
        (m / D.UNIFORM_DIR / "vasprun.xml").write_text("<x/>")
        fake = GridResolutionTests()._fake_vasprun()
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                mock.patch.object(D, "_band_edges_json", lambda cwd: {"electron": {"hits": [{"k_frac": [0, 0]}]}}):
            got, prov = D.get_effective_mass_aniso(m, "electron", True)
        self.assertIsNone(got)
        self.assertIn("二次型拟合异常：", prov)
        self.assertRegex(prov, r"gen_step12_dpt\.py:\d+")


class AnisoKzTests(unittest.TestCase):
    """V151：2D 的 S3 用 kz≥3（WS2 48×48×3）时，分方向拟合的面内掩码被拿去索引全 BZ 的 kfrac -> IndexError
    （6912 vs 2304），分方向 m* 整块失败；kz=1 时两者等长测不出来。合成六方单层，K 谷抛物 m_c=0.30、m_v=0.42。"""

    def _fake(self, n, nz):
        import types
        import numpy as np
        try:
            import spglib
            from pymatgen.core import Lattice, Structure
            from pymatgen.electronic_structure.core import Spin
        except ImportError:
            self.skipTest("需要 pymatgen + spglib")
        st = Structure(Lattice.hexagonal(3.18, 20.0), ["W", "S", "S"],
                       [[0, 0, 0.5], [1 / 3, 2 / 3, 0.578], [1 / 3, 2 / 3, 0.422]])
        mp, grid = spglib.get_ir_reciprocal_mesh([n, n, nz], (st.lattice.matrix, st.frac_coords, [74, 16, 16]),
                                                 is_shift=[0, 0, 0])
        kf = grid[np.unique(mp)] / np.array([n, n, nz], float)
        kf -= np.rint(kf)
        recip = st.lattice.reciprocal_lattice.matrix
        ks = [np.array(x + [0.0]) for x in ([1 / 3, 1 / 3], [-1 / 3, -1 / 3], [2 / 3, -1 / 3], [-2 / 3, 1 / 3],
                                             [1 / 3, -2 / 3], [-1 / 3, 2 / 3])]

        def band(e0, m, sgn):
            out = []
            for c in ks:
                d = kf - c
                d -= np.rint(d)
                out.append(e0 + sgn * 3.80998 * np.sum((d @ recip) ** 2, axis=1) / m)
            return np.min(out, axis=0) if sgn > 0 else np.max(out, axis=0)
        ene = np.stack([band(0.0, 0.42, -1) - 1.0, band(0.0, 0.42, -1), band(1.6, 0.30, +1)], axis=1)
        occ = np.stack([np.ones(len(kf))] * 2 + [np.zeros(len(kf))], axis=1)
        return types.SimpleNamespace(final_structure=st, actual_kpoints=kf.tolist(), parameters={"ISPIN": 1},
                                     eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})

    def test_kz3_grid(self):
        from unittest import mock
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        m = Path(tempfile.mkdtemp())
        (m / D.UNIFORM_DIR).mkdir(parents=True)
        (m / D.UNIFORM_DIR / "vasprun.xml").write_text("<x/>")
        for nz in (1, 3):
            fake = self._fake(12, nz)
            with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                    contextlib.redirect_stdout(io.StringIO()):
                for car, want in (("electron", 0.30), ("hole", 0.42)):
                    got, prov = D.get_effective_mass_aniso(m, car, True)
                    self.assertIsNotNone(got, "nz=%d %s：%s" % (nz, car, prov))
                    self.assertAlmostEqual(got[0], want, places=3)
                    self.assertAlmostEqual(got[1], want, places=3)

    def test_hex_offgrid_refused_v152(self):
        """V152：K 离网（N 不是 3 的倍数，WSe2 是 45）时分方向路径也要拒绝 —— 以前它给出 x≠y 且 status ok。"""
        from unittest import mock
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        m = Path(tempfile.mkdtemp())
        (m / D.UNIFORM_DIR).mkdir(parents=True)
        (m / D.UNIFORM_DIR / "vasprun.xml").write_text("<x/>")
        fake = self._fake(10, 3)
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            for car in ("electron", "hole"):
                got, prov = D.get_effective_mass_aniso(m, car, True)
                self.assertIsNone(got, prov)
                self.assertIn("K 离网", prov)
                self.assertIsNone(D.get_effective_mass(m, car, True)[0])

    def test_q_valley_not_blamed_on_grid_v152(self):
        """V152：WSe2 —— K 在网格上（45 是 3 的倍数），导带底在 Q/Λ 谷（0.178）。两条路径都拒绝，
        但说明是"非高对称谷、重做 S3 没用"，不是"K 离网"。"""
        from unittest import mock
        import numpy as np
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        from pymatgen.electronic_structure.core import Spin
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        m = Path(tempfile.mkdtemp())
        (m / D.UNIFORM_DIR).mkdir(parents=True)
        (m / D.UNIFORM_DIR / "vasprun.xml").write_text("<x/>")
        fake = self._fake(18, 1)
        kf = np.asarray(fake.actual_kpoints)
        recip = fake.final_structure.lattice.reciprocal_lattice.matrix
        ops = SpacegroupAnalyzer(fake.final_structure).get_point_group_operations(cartesian=False)
        star = {tuple(np.round(op.operate([4 / 18, 4 / 18, 0]) % 1.0, 6)) for op in ops}
        out = []
        for q in star:
            d = kf - np.array(q)
            d -= np.rint(d)
            out.append(1.5 + 3.80998 * np.sum((d @ recip) ** 2, axis=1) / 0.5)
        arr = fake.eigenvalues[Spin.up].copy()
        arr[:, 2, 0] = np.min(out, axis=0)                                  # 导带底在 Q 谷
        fake.eigenvalues = {Spin.up: arr}
        with mock.patch("pymatgen.io.vasp.Vasprun", lambda *a, **k: fake), \
                contextlib.redirect_stdout(io.StringIO()):
            for fn in (D.get_effective_mass, D.get_effective_mass_aniso):
                got, prov = fn(m, "electron", True)
                self.assertIsNone(got, prov)
                self.assertIn("非高对称谷", prov)
                self.assertNotIn("K 离网", prov)



class BandEdgesRevTests(unittest.TestCase):
    """V153：band_edges.json 的版本由 S7.1 写、S8.2 核对，两边共用 ke_common.BAND_EDGES_REV。
    以前 S8.2 拿自己的 _SKILL_REV 比，永远不一致，每次都告警。"""

    def _dpt(self):
        sys.path.insert(0, str(ROOT / "step8.2_dpt"))
        import gen_step12_dpt as D
        return D

    def test_current_file_no_warning_old_file_warns(self):
        import ke_common as kc
        D = self._dpt()
        m = Path(tempfile.mkdtemp())
        (m / D.DEFORM_READ_DIR).mkdir()
        f = m / D.DEFORM_READ_DIR / "band_edges.json"
        f.write_text(json.dumps({"skill_rev": kc.BAND_EDGES_REV, "electron": {}}))
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertIsNotNone(D._band_edges_json(m))
        self.assertNotIn("WARN", buf.getvalue())
        f.write_text(json.dumps({"skill_rev": "2026-08-01-old", "electron": {}}))
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            D._band_edges_json(m)
        self.assertIn("重跑 step7b_deform_read", buf.getvalue())

    def test_s71_writes_the_shared_constant(self):
        src = (ROOT / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py").read_text(encoding="utf-8")
        self.assertIn("_SKILL_REV = kc.BAND_EDGES_REV", src)

if __name__ == "__main__":
    unittest.main(verbosity=2)
