#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[0022] S8.2 多谷 DPT（自动找谷、逐谷质量张量与 E1、按布居加权的 μ(T)）+ S8.1/S8.3/mobility_vs_dpt 的接线。

S8.2 部分需要 pymatgen + spglib（没有就跳过），不需要 AMSET / 数据；接线部分纯本地。

合成六方单层（WS2 晶格），能带 = 各谷抛物面的 min（导带）/ max（价带），形变构型里每个谷按自己的形变势张量平移：
  导带：Q 星（Γ–K 线 0.6|K| 处，离 48 网格；纵/横质量 0.80/0.50，纵/横形变势 6.0/3.0 eV）为带边，
        K 星（m=0.40，D=4.0 eV）高 30 meV —— 像 WSe2：单谷 DPT 拒绝出数，多谷照算；
  价带：K 星（m=0.45，D=2.5 eV）为带边，Γ（m=1.20，D=5.0 eV）低 80 meV —— 像 MoS2。
期望值按 μ_αα = (e/m0)·Σ n_v τ_α,v W_αα,v / Σ n_v 独立算（不调被测代码），六方取 x/y 平均。
"""
import contextlib
import io
import json
import math
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1] / "skill" / "ke-dft-cpu"
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "step8.2_dpt"), str(ROOT / "tools")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import numpy as np
    import spglib  # noqa: F401
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Spin
    HAVE = True
except ImportError:                                   # pragma: no cover
    HAVE = False

E_C, HBAR, KB, M0 = 1.602176634e-19, 1.054571817e-34, 1.380649e-23, 9.1093837015e-31
HB = 3.80998
C2D = 120.0
DVAC = {"xx": 0.70, "yy": -0.40}                       # 真空势应变导数：本征值里带着它，E1 要减掉
G = 0.005


def _rot(t):
    c, s = math.cos(t), math.sin(t)
    return np.array([[c, -s], [s, c]])


class Model:
    def __init__(self, eK_cb=0.03, eG_vb=-0.08, with_q=True, with_gamma=True):
        self.st = Structure(Lattice.hexagonal(3.18, 20.0), ["W", "S", "S"],
                            [[0, 0, 0.5], [1 / 3, 2 / 3, 0.578], [1 / 3, 2 / 3, 0.422]])
        self.recip = self.st.lattice.reciprocal_lattice.matrix
        kc = np.array([1 / 3, 1 / 3, 0.0]) @ self.recip
        th0 = math.atan2(kc[1], kc[0])
        cb, vb = [], []
        for kf in ([1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]):                     # K、K'
            cb.append(self._v(kf, 1.6 + eK_cb, np.eye(2) / 0.40, 4.0 * np.eye(2), +1, "K"))
            vb.append(self._v(kf, 0.0, np.eye(2) / 0.45, 2.5 * np.eye(2), -1, "K"))
        if with_q:
            for n in range(6):                                                  # Q 星
                t = th0 + n * math.pi / 3
                qc = 0.6 * np.linalg.norm(kc[:2]) * np.array([math.cos(t), math.sin(t), 0.0])
                qf = np.linalg.solve(self.recip.T, qc)
                R = _rot(t)
                cb.append(self._v(qf, 1.6, R @ np.diag([1 / 0.80, 1 / 0.50]) @ R.T,
                                  R @ np.diag([6.0, 3.0]) @ R.T, +1, "Q"))
        if with_gamma:
            vb.append(self._v([0, 0, 0], eG_vb, np.eye(2) / 1.20, 5.0 * np.eye(2), -1, "Γ"))
        self.cb, self.vb = cb, vb

    @staticmethod
    def _v(kf, e0, W, D, sgn, lab):
        return {"k": np.array(kf, float), "e0": e0, "W": np.asarray(W, float), "D": np.asarray(D, float),
                "sgn": sgn, "lab": lab}

    def band(self, kf, valleys, strain=None):
        kf = np.asarray(kf, float)
        out = []
        for v in valleys:
            d = kf - v["k"]
            d -= np.rint(d)
            q = (d @ self.recip)[:, :2]
            e0 = v["e0"]
            if strain:
                ax, eps = strain
                i = {"xx": 0, "yy": 1}[ax]
                e0 = e0 + (v["D"][i, i] + DVAC[ax]) * eps
            out.append(e0 + v["sgn"] * HB * np.einsum("ni,ij,nj->n", q, v["W"], q))
        return np.min(out, axis=0) if valleys[0]["sgn"] > 0 else np.max(out, axis=0)

    def bands(self, kf, strain=None):
        vb = self.band(kf, self.vb, strain)
        return np.stack([vb - 1.0, vb, self.band(kf, self.cb, strain)], axis=1)

    def vasprun(self, kf, strain=None, nocc=2):
        import types
        ene = self.bands(kf, strain)
        occ = np.zeros_like(ene)
        occ[:, :nocc] = 1.0
        return types.SimpleNamespace(final_structure=self.st, actual_kpoints=np.asarray(kf).tolist(),
                                     parameters={"ISPIN": 1},
                                     eigenvalues={Spin.up: np.stack([ene, occ], axis=-1)})

    def s3_kpoints(self, n, nz):
        import spglib
        mp, grid = spglib.get_ir_reciprocal_mesh(
            [n, n, nz], (self.st.lattice.matrix, self.st.frac_coords, [74, 16, 16]), is_shift=[0, 0, 0])
        kf = grid[np.unique(mp)] / np.array([n, n, nz], float)
        return kf - np.rint(kf)

    @staticmethod
    def s7_kpoints(n):
        """VASP ISYM=0：全网格，但按时间反演只留 k / −k 之一。"""
        seen, out = set(), []
        for i in range(n):
            for j in range(n):
                k = np.array([i / n, j / n, 0.0])
                k -= np.rint(k)
                key = tuple(np.round(k % 1.0, 6) % 1.0)
                mkey = tuple(np.round((-k) % 1.0, 6) % 1.0)
                if mkey in seen:
                    continue
                seen.add(key)
                out.append(k)
        return np.array(out)

    # ---- 期望值（独立实现，不调被测代码）----
    def expect(self, carrier, T, window=0.30):
        vs = self.cb if carrier == "electron" else self.vb
        sg = 1.0 if carrier == "electron" else -1.0
        ref = min(sg * v["e0"] for v in vs)
        kT = KB * T / E_C
        sn, sntw, snw, per = 0.0, [0.0, 0.0], [0.0, 0.0], {}
        for v in vs:
            dE = sg * v["e0"] - ref
            if dE > window:
                continue
            md = 1.0 / math.sqrt(np.linalg.det(v["W"]))
            n = md * math.exp(-dE / kT)
            for a in (0, 1):
                tau = HBAR ** 3 * C2D / (KB * T * md * M0 * (abs(v["D"][a, a]) * E_C) ** 2)
                sntw[a] += n * tau * v["W"][a, a]
                snw[a] += n * v["W"][a, a]
            sn += n
            per[v["lab"]] = per.get(v["lab"], 0.0) + n
        mu = E_C / M0 * (sntw[0] + sntw[1]) / 2.0 / sn * 1e4
        return mu, (sntw[0] + sntw[1]) / (snw[0] + snw[1]), {k: x / sn for k, x in per.items()}


def build(model, m, n3=48, nz=1, n7=12, s7_nocc=2, e1_bias=1.0):
    """在 m 下建一个材料目录：S3/S7 的 vasprun（占位文件，内容由 mock 给）、band_edges.json。返回 Vasprun 分发器。"""
    D = _D()
    (m / D.UNIFORM_DIR).mkdir(parents=True, exist_ok=True)
    (m / D.UNIFORM_DIR / "vasprun.xml").write_text("<x/>")
    k3 = model.s3_kpoints(n3, nz)
    k7 = model.s7_kpoints(n7)
    pairs = {"xx": ["deform-01", "deform-02"], "yy": ["deform-03", "deform-04"]}
    strains = {"undeformed": None, "deform-01": ("xx", G), "deform-02": ("xx", -G),
               "deform-03": ("yy", G), "deform-04": ("yy", -G)}
    root = m / D.DEFORM_DIR
    read_paths = {}
    for fold in strains:
        d = root / fold / ("ionrelax" if fold != "undeformed" else "")
        d.mkdir(parents=True, exist_ok=True)
        (d / "vasprun.xml").write_text("<x/>")
        read_paths[fold] = {"vasprun.xml": "rigid" if fold == "undeformed" else "ionrelax"}
    fakes = {fold: model.vasprun(k7, st, nocc=(s7_nocc if fold == "undeformed" else 2))
             for fold, st in strains.items()}
    fakes["s3"] = model.vasprun(k3)

    def e1(ax, k, b):
        ep = model.bands([k], (ax, G))[0, b]
        em = model.bands([k], (ax, -G))[0, b]
        return abs((ep - em) / (2 * G) - DVAC[ax])
    be = {"skill_rev": _kc().BAND_EDGES_REV, "read_paths": read_paths,
          "vac_align": {"status": "ok", "pairs": pairs, "strain_mag": {"xx": G, "yy": G},
                        "dvac_eV_per_unit_strain": dict(DVAC)}}
    ene7 = model.bands(k7)
    for car, b, pick in (("electron", 2, np.argmin), ("hole", 1, np.argmax)):
        i = int(pick(ene7[:, b]))
        kf = k7[i]
        xx, yy = e1("xx", kf, b) * e1_bias, e1("yy", kf, b)
        s = (xx + yy) / 2
        be[car] = {"edge_energy_eV": float(ene7[i, b]), "n_hits": 1,
                   "hits": [{"spin": "up", "band_global": b, "band_h5": b, "k_h5": i,
                             "k_frac": np.round(kf, 6).tolist()}],
                   "E1_xx_eV": s, "E1_yy_eV": s, "E1_iso_eV": s,
                   "E1_vac_xx_raw_eV": round(xx, 4), "E1_vac_yy_raw_eV": round(yy, 4),
                   "E1_vac_xx_eV": round(s, 4), "E1_vac_yy_eV": round(s, 4), "E1_vac_iso_eV": round(s, 4)}
    (m / D.DEFORM_READ_DIR).mkdir(parents=True, exist_ok=True)
    (m / D.DEFORM_READ_DIR / "band_edges.json").write_text(json.dumps(be))

    def vr(path, *a, **k):
        p = str(path)
        if D.UNIFORM_DIR in p:
            return fakes["s3"]
        for fold in strains:
            if "/%s/" % fold in p:
                return fakes[fold]
        raise AssertionError("意外的 vasprun 路径 %s" % p)
    return vr


def _D():
    import gen_step12_dpt as D
    return D


def _kc():
    import ke_common as kc
    return kc


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class MultiValleyTests(unittest.TestCase):
    def setUp(self):
        self.D = _D()
        self.D._MV_CACHE.clear()
        self.m = Path(tempfile.mkdtemp())
        self._p = [mock.patch.dict(self.D.MANUAL_ANISO, {"C_2D_xy_N_per_m": (C2D, C2D)}),
                   mock.patch.object(self.D, "VALLEY_WINDOW_EV", 0.30),
                   mock.patch.object(self.D, "DPT_VALLEYS", "auto"),
                   mock.patch.object(self.D, "E1_SOURCE", "vac")]
        for p in self._p:
            p.start()

    def tearDown(self):
        for p in self._p:
            p.stop()

    def _run(self, model, carrier, temps=(300.0,), **kw):
        vr = build(model, self.m, **kw)
        with mock.patch("pymatgen.io.vasp.Vasprun", vr), contextlib.redirect_stdout(io.StringIO()):
            bd = self.D._aniso_block(self.m, True, carrier, temps[0])
            mv = self.D._multi_valley_block(self.m, True, carrier, temps[0], list(temps), bd)
        return bd, mv

    def test_q_edge_electron_analytic(self):
        """WSe2 型：导带底在 Q（离网），单谷拒绝；多谷找到 Q×6 + K×2，μ、τ、布居与解析值一致。"""
        model = Model()
        bd, mv = self._run(model, "electron", temps=(300.0, 600.0))
        self.assertNotEqual(bd.get("status"), "ok")
        self.assertEqual(mv["status"], "ok", mv["status"])
        labs = {v["label"]: v for v in mv["valleys"]}
        self.assertEqual(set(labs), {"Q", "K"})
        self.assertEqual(labs["Q"]["g"], 6)
        self.assertEqual(labs["K"]["g"], 2)
        self.assertAlmostEqual(labs["Q"]["dE_eV"], 0.0, places=4)
        self.assertAlmostEqual(labs["K"]["dE_eV"], 0.03, places=4)       # 拟合出的极值，不是网格点能量
        self.assertAlmostEqual(labs["Q"]["m_d_m0"], math.sqrt(0.8 * 0.5), places=3)
        self.assertAlmostEqual(labs["K"]["m_d_m0"], 0.40, places=3)
        for T, row in zip((300.0, 600.0), mv["mobility_vs_T"]):
            mu, tau, w = model.expect("electron", T)
            self.assertAlmostEqual(row["mobility_x_cm2_Vs"] / mu, 1.0, places=3)
            self.assertAlmostEqual(row["mobility_y_cm2_Vs"] / mu, 1.0, places=3)
            self.assertAlmostEqual(row["tau_x_s"] / tau, 1.0, places=3)
            self.assertAlmostEqual(row["weights"]["Q"], w["Q"], places=3)
        mu300, _t, _w = model.expect("electron", 300.0)
        self.assertAlmostEqual(mv["x"]["mobility_cm2_Vs"] / mu300, 1.0, places=3)
        # E1 逐谷：Q 成员各向异性（转过去的 6.0/3.0），K 各向同性 4.0（C₃ 位点对称化）
        self.assertTrue(any(abs(e[0] - e[1]) > 0.5 for e in labs["Q"]["E1_xy_eV"]))
        self.assertTrue(all(abs(e[0] - 4.0) < 1e-3 and abs(e[1] - 4.0) < 1e-3 for e in labs["K"]["E1_xy_eV"]))
        self.assertEqual(mv["E1_selfcheck"]["xx"]["recomputed_eV"], 4.0)

    def test_k_edge_hole_two_stars_and_single_valley_limit(self):
        """MoS2 型价带：K 带边 + Γ 低 80 meV。多谷 = 解析值；只留 K 星与 by_direction 一致（单谷极限）。"""
        model = Model()
        bd, mv = self._run(model, "hole", temps=(300.0,))
        self.assertEqual(bd.get("status"), "ok", bd.get("status"))
        self.assertEqual(mv["status"], "ok", mv["status"])
        self.assertEqual({v["label"]: v["g"] for v in mv["valleys"]}, {"K": 2, "Γ": 1})
        mu, _tau, _w = model.expect("hole", 300.0)
        self.assertAlmostEqual(mv["x"]["mobility_cm2_Vs"] / mu, 1.0, places=3)
        cons = mv["edge_star_vs_by_direction"]
        self.assertEqual(cons["star"], "K")
        self.assertAlmostEqual(cons["x"], 1.0, places=3)
        self.assertAlmostEqual(cons["y"], 1.0, places=3)
        self.assertEqual(mv["warnings"], [])

    def test_window_excludes_far_valleys(self):
        """窗口 0.05 eV：价带只剩 K 星 -> 多谷 μ 与 by_direction 完全相同，且 μ∝1/T。"""
        model = Model()
        with mock.patch.object(self.D, "VALLEY_WINDOW_EV", 0.05):
            bd, mv = self._run(model, "hole", temps=(300.0, 600.0))
        self.assertEqual(mv["status"], "ok", mv["status"])
        self.assertEqual([v["label"] for v in mv["valleys"]], ["K"])
        self.assertAlmostEqual(mv["x"]["mobility_cm2_Vs"] / bd["x"]["mobility_cm2_Vs"], 1.0, places=3)
        r300, r600 = mv["mobility_vs_T"]
        self.assertAlmostEqual(r300["mobility_x_cm2_Vs"] / r600["mobility_x_cm2_Vs"], 2.0, places=3)

    def test_temperature_dependence_not_simple_1_over_T(self):
        """两个星时布居随 T 变：μ·T 不是常数，高温时高能谷布居变大。"""
        model = Model()
        _bd, mv = self._run(model, "hole", temps=(100.0, 300.0, 900.0))
        rows = mv["mobility_vs_T"]
        muT = [r["mobility_x_cm2_Vs"] * r["T_K"] for r in rows]
        self.assertGreater(abs(muT[0] / muT[2] - 1.0), 0.05)
        self.assertLess(rows[0]["weights"]["Γ"], rows[2]["weights"]["Γ"])
        self.assertIn(900.0, mv["window_truncation_T_K"])

    def test_kz3_grid_same_result(self):
        """2D 的 S3 用 kz=3（WS2 48×48×3）：只取 kz=0 层，结果与 kz=1 相同。"""
        model = Model()
        _bd, mv1 = self._run(model, "electron")
        self.D._MV_CACHE.clear()
        self.m = Path(tempfile.mkdtemp())
        _bd, mv3 = self._run(model, "electron", nz=3)
        self.assertEqual(mv3["status"], "ok", mv3["status"])
        self.assertAlmostEqual(mv3["x"]["mobility_cm2_Vs"] / mv1["x"]["mobility_cm2_Vs"], 1.0, places=6)

    def test_band_count_mismatch_refused(self):
        """S7 undeformed 的占据带数与 S3 不同（NELECT/SOC 不一致）-> 拒绝出数，说清原因。"""
        _bd, mv = self._run(Model(), "electron", s7_nocc=1)
        self.assertNotEqual(mv["status"], "ok")
        self.assertIn("占据带数不同", mv["status"])

    def test_edge_e1_selfcheck_refused(self):
        """band_edges.json 的 E1_vac 与按同一定义重算的对不上 -> 拒绝（构型/带号/自旋读错了）。"""
        _bd, mv = self._run(Model(), "electron", e1_bias=1.2)
        self.assertNotEqual(mv["status"], "ok")
        self.assertIn("≠ S7.1", mv["status"])

    def test_missing_deform_vasprun_refused(self):
        model = Model()
        vr = build(model, self.m)
        (self.m / "step7_deform" / "deform-03" / "ionrelax" / "vasprun.xml").unlink()
        with mock.patch("pymatgen.io.vasp.Vasprun", vr), contextlib.redirect_stdout(io.StringIO()):
            mv = self.D._multi_valley_block(self.m, True, "electron", 300.0, [300.0], {})
        self.assertIn("缺 step7_deform/deform-03/ionrelax/vasprun.xml", mv["status"])

    def test_switches(self):
        with mock.patch.object(self.D, "DPT_VALLEYS", "off"):
            self.assertTrue(self.D._multi_valley_block(self.m, True, "electron", 300.0, [300.0])["status"]
                            .startswith("off"))
        self.assertIn("只实现 2D", self.D._multi_valley_block(self.m, False, "electron", 300.0, [300.0])["status"])
        build(Model(), self.m)
        with mock.patch.object(self.D, "E1_SOURCE", "amset"):
            st = self.D._multi_valley_block(self.m, True, "electron", 300.0, [300.0])["status"]
        self.assertIn("E1_SOURCE=vac", st)

    def test_q_message_points_to_multi_valley(self):
        """单谷对 Q 带边的拒绝说明：仍说"非高对称谷"，不再叫人手填 M_EFF_*，而是指向 multi_valley。"""
        vr = build(Model(), self.m)
        with mock.patch("pymatgen.io.vasp.Vasprun", vr), contextlib.redirect_stdout(io.StringIO()):
            got, prov = self.D.get_effective_mass_aniso(self.m, "electron", True)
        self.assertIsNone(got)
        self.assertIn("非高对称谷", prov)
        self.assertIn("multi_valley", prov)
        self.assertNotIn("K 离网", prov)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class CollectTests(unittest.TestCase):
    """按能量排序的两条导带相交：上面那条在交叉线上是 V 形"极小"，不是谷 -> 跳过并记下；真谷照收。"""

    def test_band_crossing_kink_skipped(self):
        D = _D()
        n = 24
        i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        low = 1.0 + 0.2 * (2 - np.cos(2 * np.pi * i / n) - np.cos(2 * np.pi * j / n))   # 极小在 (0,0)
        up = low + 0.3 * np.abs(np.sin(np.pi * (i - n / 2) / n))                          # i=n/2 线上与 low 相交
        idx = np.arange(n * n).reshape(n, n)
        ene = {"up": np.stack([np.zeros(n * n), low.ravel(), up.ravel()], axis=1)}
        valleys, seeds, skipped = D._mv_collect(ene, idx, (n, n), {"up": 1}, +1.0, 1.0, 0.5)
        got = {(b, pts[0]) for sp, b, pts, e in valleys}
        self.assertIn((1, (0, 0)), got)                       # 下面那条带的真谷
        self.assertIn((2, (0, 0)), got)                       # 上面那条带在 (0,0) 是真谷（离 low 0.3 eV）
        self.assertNotIn((2, (n // 2, 0)), got)               # 交叉点：跳过
        self.assertEqual(len(skipped), 1)
        self.assertIn("交叉", skipped[0])
        self.assertEqual(seeds, {k for k, (sp, b, pts, e) in enumerate(valleys) if e - 1.0 <= 0.5})


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class ExpandTests(unittest.TestCase):
    """pymatgen 点群展开：倒空间用 Rᵀ；区界 ±0.5 是同一点（以前 12×12 只得 97 点 / 修了转置后又多出 147 点）。"""

    def test_full_grid_even_odd(self):
        import types
        D = _D()
        model = Model()
        for n, nz in ((12, 1), (12, 3), (9, 1), (10, 2)):
            kf = model.s3_kpoints(n, nz)
            v = types.SimpleNamespace(final_structure=model.st, actual_kpoints=kf.tolist(), parameters={"ISPIN": 1})
            with mock.patch.dict(sys.modules, {"amset": None, "amset.electronic_structure": None,
                                               "amset.electronic_structure.symmetry": None}), \
                    mock.patch.dict(os.environ, {}), contextlib.redirect_stdout(io.StringIO()):
                full, kmap, how, err = D._expand_full_bz(v, Path(tempfile.mkdtemp()), "electron")
            self.assertIsNone(err, err)
            self.assertEqual(how, "pymatgen")
            self.assertEqual(len(full), n * n * nz, (n, nz))
            keys = {tuple(np.round(np.mod(k * [n, n, nz], [n, n, nz])).astype(int) % [n, n, nz]) for k in full}
            self.assertEqual(len(keys), n * n * nz)


@unittest.skipUnless(HAVE, "需要 pymatgen + spglib")
class MainEndToEndTests(unittest.TestCase):
    """main()：电子单谷拒绝（Q 带边）但多谷出数 -> 不报"全部未算出"；json/summary 带 multi_valley。"""

    def test_main_writes_multi_valley(self):
        D = _D()
        D._MV_CACHE.clear()
        m = Path(tempfile.mkdtemp())
        model = Model()
        vr = build(model, m)
        (m / "step1_opt").mkdir()
        (m / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\n")
        model.st.to(filename=str(m / D.UNIFORM_DIR / "POSCAR"), fmt="poscar")
        (m / D.UNIFORM_DIR / "KPOINTS").write_text("auto\n0\nGamma\n48 48 1\n0 0 0\n")
        h = abs(np.dot(model.st.lattice.matrix[2], np.cross(model.st.lattice.matrix[0], model.st.lattice.matrix[1]))) \
            / np.linalg.norm(np.cross(model.st.lattice.matrix[0], model.st.lattice.matrix[1]))
        c_kbar = C2D / (h * 1e-10) / 1e9 * 10.0
        rows = "\n".join(" XX %s" % " ".join("%.4f" % (c_kbar if i == j and i < 2 else 1.0) for j in range(6))
                         for i in range(6))
        (m / D.ELASTIC_DIR).mkdir()
        (m / D.ELASTIC_DIR / "OUTCAR").write_text(" TOTAL ELASTIC MODULI (kBar)\n Direction XX YY\n" + rows + "\n")
        cwd = os.getcwd()
        os.chdir(m)
        try:
            with mock.patch("pymatgen.io.vasp.Vasprun", vr), \
                    mock.patch.object(D, "TEMPERATURES", "300,600"), \
                    mock.patch.object(D, "MANUAL_ANISO", dict(D.MANUAL_ANISO)), \
                    contextlib.redirect_stdout(io.StringIO()) as buf:
                D.main()
        finally:
            os.chdir(cwd)
        j = json.loads((m / D.OUTDIR_NAME / "dpt_result.json").read_text())
        r = {x["carrier"]: x for x in j["results"]}
        self.assertIsNone(r["electron"]["mobility_cm2_Vs"])
        self.assertEqual(r["electron"]["multi_valley"]["status"], "ok")
        self.assertEqual(r["hole"]["multi_valley"]["status"], "ok")
        mu, _t, _w = model.expect("electron", 300.0)
        self.assertAlmostEqual(r["electron"]["multi_valley"]["x"]["mobility_cm2_Vs"] / mu, 1.0, places=2)
        self.assertEqual([t["T_K"] for t in r["electron"]["multi_valley"]["mobility_vs_T"]], [300.0, 600.0])
        txt = (m / D.OUTDIR_NAME / "dpt_summary.txt").read_text()
        self.assertIn("[多谷]", txt)
        self.assertIn("多谷 μ(T)", txt)
        self.assertIn("多谷 DPT 已算出", buf.getvalue())


# ============================ 下游接线（纯本地） ============================
def _load(name, path):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _mv_block(scale=1.0):
    rows = []
    for T, w in ((300.0, 0.94), (600.0, 0.88)):
        mx, my = 400.0 * scale * 300.0 / T * (1.0 if T == 300.0 else 0.95), 400.0 * scale * 300.0 / T
        rows.append({"T_K": T, "mobility_x_cm2_Vs": mx, "mobility_y_cm2_Vs": my,
                     "mobility_inplane_cm2_Vs": (mx + my) / 2.0,
                     "tau_x_s": 1.5e-13 * 300.0 / T, "tau_y_s": 1.6e-13 * 300.0 / T,
                     "weights": {"Q": w, "K": 1 - w}})
    return {"status": "ok", "T_K": 300.0, "window_eV": 0.3,
            "valleys": [{"label": "Q", "g": 6, "dE_eV": 0.0}, {"label": "K", "g": 2, "dE_eV": 0.03}],
            "x": {"mobility_cm2_Vs": rows[0]["mobility_x_cm2_Vs"], "tau_s": rows[0]["tau_x_s"], "m_c_m0": 0.60},
            "y": {"mobility_cm2_Vs": rows[0]["mobility_y_cm2_Vs"], "tau_s": rows[0]["tau_y_s"], "m_c_m0": 0.62},
            "mobility_vs_T": rows}


def _dpt_json(m):
    """WSe2 型：电子单谷算不出、只有多谷；空穴单谷 + 多谷都有。"""
    d = m / "step8.2_dpt"
    d.mkdir(parents=True, exist_ok=True)
    bd_h = {"status": "ok", "m_d_m0": 0.45,
            "x": {"mobility_cm2_Vs": 2000.0, "tau_s": 5.0e-13}, "y": {"mobility_cm2_Vs": 2000.0, "tau_s": 5.0e-13}}
    res = {"dim": "2d", "is_2d": True, "temperature_K": 300.0, "temperatures_K": [300.0, 600.0],
           "results": [{"carrier": "electron", "mobility_cm2_Vs": None, "inputs": {"m_eff_m0": None},
                        "by_direction": {"status": "缺输入：m*"}, "multi_valley": _mv_block()},
                       {"carrier": "hole", "mobility_cm2_Vs": 1900.0, "inputs": {"m_eff_m0": 0.45},
                        "by_direction": bd_h, "multi_valley": _mv_block(4.5)}]}
    (d / "dpt_result.json").write_text(json.dumps(res))
    return res


class S81TauTests(unittest.TestCase):
    def test_tau_falls_back_to_multi_valley(self):
        sys.path.insert(0, str(ROOT.parent / "_common" / "opt"))
        B = _load("bt2_0022", ROOT / "step8.1_boltztrap" / "gen_step11_boltztrap.py")
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            t600 = B._tau_from_json(m, 600.0, True)
            t450 = B._tau_from_json(m, 450.0, True)
            t1000 = B._tau_from_json(m, 1000.0, True)
        # 电子：多谷 τ 按温度直接取（不是 1/T 换算），两档之间 log–log 插值，网格外不外推
        self.assertAlmostEqual(t600["electron"]["x"], 1.5e-13 * 0.5, delta=1e-20)
        self.assertAlmostEqual(t450["electron"]["y"] / (1.6e-13 * 300 / 450), 1.0, places=9)
        self.assertNotIn("electron", t1000)
        # 空穴：单谷分方向优先（行为不变），×T0/T
        self.assertAlmostEqual(t600["hole"]["x"], 2.5e-13, delta=1e-20)
        self.assertIn("多谷", buf.getvalue())


class S83Tests(unittest.TestCase):
    def _mod(self):
        return _load("out_0022", ROOT / "step8.3_output" / "gen_step13_output.py")

    def test_load_dpt_uses_multi_valley_when_single_missing(self):
        O = self._mod()
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        with contextlib.redirect_stdout(io.StringIO()):
            d = O.load_dpt(m)
        self.assertAlmostEqual(d["electron"], 400.0)                  # 300 K 面内平均
        self.assertEqual(d["src_electron"], "multi_valley")
        self.assertEqual(d["dir_electron"]["x"]["tau_s"], 1.5e-13)
        self.assertAlmostEqual(O._dpt_tau_dir(d, "electron", "y"), 1.6e-13)
        self.assertAlmostEqual(d["hole"], 2000.0)                     # 空穴主栏仍是单谷
        self.assertAlmostEqual(d["mv_hole"]["mean"], 1800.0)
        rows = O.build_table(None, None, d)
        e = next(r for r in rows if r["carrier_conc_cm-3"] < 0)
        self.assertAlmostEqual(e["dpt_mv_mu_cm2/Vs"], 400.0)

    def test_mu_vs_T_has_mv_column(self):
        O = self._mod()
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        res = O.build_mu_vs_T(m, True)
        self.assertIn("mu_DPT_mv", O.MU_T_COLS)
        e600 = next(r for r in res["rows"] if r["carrier"] == "electron" and r["T_K"] == 600.0
                    and r["direction"] == "mean")
        self.assertIsNone(e600["mu_DPT"])
        self.assertAlmostEqual(e600["mu_DPT_mv"], (190.0 + 200.0) / 2.0)
        h300 = next(r for r in res["rows"] if r["carrier"] == "hole" and r["T_K"] == 300.0
                    and r["direction"] == "x")
        self.assertAlmostEqual(h300["mu_DPT"], 2000.0)
        self.assertAlmostEqual(h300["mu_DPT_mv"], 1800.0)
        self.assertIn("Q×6", res["meta"]["dpt_mv"]["electron"])


class ToolTests(unittest.TestCase):
    def test_rows_and_print_have_dpt_mv(self):
        import mobility_vs_dpt as T
        m = Path(tempfile.mkdtemp())
        _dpt_json(m)
        mvs = T.dpt_mv(m / "step8.2_dpt" / "dpt_result.json")
        diag = lambda v: [[v, 0, 0], [0, v, 0], [0, 0, 1]]
        tj = {"doping": [-1e17, 1e17], "temperatures": [300.0, 600.0],
              "seebeck": [[diag(-300.0)] * 2, [diag(300.0)] * 2],
              "mobility": {"overall": [[diag(30.0)] * 2, [diag(90.0)] * 2],
                           "ADP": [[diag(800.0), diag(380.0)], [diag(3600.0), diag(1700.0)]]}}
        rs, mechs = T.rows(tj, 600.0, True, {"electron": None, "hole": 2000.0}, None, None, dpt_T0=300.0, mv=mvs)
        e, h = rs
        self.assertIsNone(e["DPT"])
        self.assertAlmostEqual(e["DPT_mv"], 195.0)
        self.assertAlmostEqual(e["ADP_over_DPT_mv"], round(380.0 / 195.0, 3))
        self.assertAlmostEqual(h["DPT"], 1000.0)
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            T._print_block(Path("t.json"), True, rs, mechs, {"electron": None, "hole": 2000.0}, {}, 300.0,
                           set(), set(), 600.0)
        self.assertIn("DPT_mv", buf.getvalue())
        rs0, _ = T.rows(tj, 600.0, True, {"electron": None, "hole": 2000.0}, None, None, dpt_T0=300.0)
        self.assertNotIn("DPT_mv", rs0[0])                            # 不给 mv：输出与 0020 相同


if __name__ == "__main__":
    unittest.main()
