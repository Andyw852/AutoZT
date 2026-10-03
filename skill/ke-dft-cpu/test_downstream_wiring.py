#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_downstream_wiring.py —— V148：上游重算 -> 下游失效的接线（纯 Python，不需要 VASP/AMSET）。

用法：python test_downstream_wiring.py      退出码 0 = 全部 PASS。

起因：DOWNSTREAM 表里早有 S8.2 -> S8.1/S8.3，可没有一个 gen 调 invalidate_downstream —— GaAs 的 S8.2 重跑
4 次（V145 空穴 276 万 -> 4527），S8.1 一直用旧 τ；S2 画图重新生成后 S8/S8.4 的 settings.yaml 仍是旧带隙。
另：S8.2 的 MANUAL 只能改共享脚本，改着的时候所有材料的 S8.2 都会用上（gen 脚本总是从 skill 覆盖推送）。
1) skill.yaml 的每条 needs 边要么在 DOWNSTREAM、要么在 DOWNSTREAM_EXEMPT 写明理由；
2) DOWNSTREAM 的每个上游：它的 gen 调 invalidate_downstream / finalize_stale_inputs，gen_need 带 ke_common；
3) mode="gap"：只有下游用的是这份 band_summary、且带隙变了才失效；没跑完的下游给 ★；
4) S8.2：step.conf 覆盖只作用本材料并记进 overrides；脚本 MANUAL 被改过时 ★ 告警；重跑归档 S8.1/S8.3；
5) tools/lineage_check：S8.1/S8.3 旧于来源、S8/S8.4 的带隙与 S2 画图不一致（S8 的 gen 闸门不查）。
"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "step8.2_dpt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


def _steps(path):
    out = []

    def walk(o):
        if isinstance(o, dict):
            if "name" in o and ("gen" in o or "needs" in o):
                out.append(o)
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(yaml.safe_load(Path(path).read_text(encoding="utf-8")))
    return out


STEPS = _steps(ROOT / "skill.yaml")
BY_NAME = {s["name"]: s for s in STEPS}


def _gen_source(step):
    s = BY_NAME[step]
    gen = str(s["gen"]).split()[0]
    return (ROOT / s["src"] / gen).read_text(encoding="utf-8")


class WiringTests(unittest.TestCase):
    def test_every_needs_edge_is_covered(self):
        miss = []
        for s in STEPS:
            for up in s.get("needs") or []:
                if s["name"] in [d for d, _m in kc.DOWNSTREAM.get(up, ())]:
                    continue
                if (up, s["name"]) in kc.DOWNSTREAM_EXEMPT or (up, "*") in kc.DOWNSTREAM_EXEMPT:
                    continue
                miss.append("%s -> %s" % (up, s["name"]))
        self.assertEqual(miss, [], "needs 边既不在 DOWNSTREAM 也不在 DOWNSTREAM_EXEMPT")

    def test_every_upstream_gen_invalidates(self):
        bad = []
        for up in kc.DOWNSTREAM:
            self.assertIn(up, BY_NAME, up)
            src = _gen_source(up)
            if "invalidate_downstream(" not in src and "finalize_stale_inputs(" not in src:
                bad.append(up + "：gen 不调 invalidate_downstream")
            need = BY_NAME[up].get("gen_need")
            if need is not None and "ke_common.py" not in need:
                bad.append(up + "：gen_need 没有 ke_common.py")
            for d, _m in kc.DOWNSTREAM[up]:
                self.assertIn(d, kc.DONE_MARKERS, "%s 没有 DONE_MARKERS，失效不了" % d)
        self.assertEqual(bad, [])

    def test_zt_gen_need_matches(self):
        zt = {s["name"]: s for s in _steps(ROOT.parent / "zt-dft-cpu" / "skill.yaml")}
        for up in kc.DOWNSTREAM:
            if up in zt and zt[up].get("gen_need") is not None:
                self.assertIn("ke_common.py", zt[up]["gen_need"], "zt " + up)


def _mat():
    m = Path(tempfile.mkdtemp())
    for d in ("step8_amset", "step8.4_amset2d", "step8.1_boltztrap", "step8.3_output", "step8.2_dpt",
              "step2_bandgap/step2.3_hse_plot", "step2_bandgap/step2.2_pbe_plot"):
        (m / d).mkdir(parents=True)
    return m


def _settings(d, gap, src):
    (d / "settings.yaml").write_text("scissor: null\nbandgap: %s\n# bandgap_source: %s\nnworkers: 4\n" % (gap, src))


class GapTests(unittest.TestCase):
    def _run(self, m, step, gap):
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            done = kc.invalidate_downstream(m, step, "test", gap=gap)
        return done, buf.getvalue()

    def test_gap_changed_invalidates_consumer_and_s83(self):
        m = _mat()
        _settings(m / "step8_amset", 1.4200, "step2.3_hse_plot/band_summary.json（gap_eV）")
        (m / "step8_amset" / "transport.json").write_text("{}")
        (m / "step8.3_output" / "comparison_300K.png").write_text("x")
        (m / "step8.1_boltztrap" / "boltztrap_crta.json").write_text("{}")
        done, out = self._run(m, "step2_bandgap/step2.3_hse_plot", 1.3000)
        self.assertIn(("step8_amset", "transport.json"), done)
        self.assertIn(("step8.3_output", "comparison_300K.png"), done)     # 递归
        self.assertIn(("step8.1_boltztrap", "boltztrap_crta.json"), done)
        self.assertIn("1.4200 -> 1.3000", out)

    def test_same_gap_or_other_source_untouched(self):
        m = _mat()
        _settings(m / "step8_amset", 1.4200, "step2.3_hse_plot/band_summary.json（gap_eV）")
        (m / "step8_amset" / "transport.json").write_text("{}")
        _settings(m / "step8.4_amset2d", 1.4200, "step.conf BANDGAP_OVERRIDE")
        (m / "step8.4_amset2d" / "transport.json").write_text("{}")
        done, _ = self._run(m, "step2_bandgap/step2.3_hse_plot", 1.4205)       # 0.5 meV：不动
        self.assertNotIn("step8_amset", [d for d, _f in done])
        done, _ = self._run(m, "step2_bandgap/step2.2_pbe_plot", 0.80)         # S8 用的是 HSE 那份
        self.assertNotIn("step8_amset", [d for d, _f in done])
        done, _ = self._run(m, "step2_bandgap/step2.3_hse_plot", 1.30)         # S8.4 用的是覆盖值
        self.assertIn("step8_amset", [d for d, _f in done])
        self.assertNotIn("step8.4_amset2d", [d for d, _f in done])
        done, _ = self._run(m, "step2_bandgap/step2.3_hse_plot", None)         # 不给 gap：gap 下游不动
        self.assertNotIn("step8.4_amset2d", [d for d, _f in done])

    def test_running_consumer_warned(self):
        m = _mat()
        _settings(m / "step8.4_amset2d", 2.10, "step2.3_hse_plot/band_summary.json（gap_eV）")   # 无 transport.json
        done, out = self._run(m, "step2_bandgap/step2.3_hse_plot", 1.95)
        self.assertNotIn("step8.4_amset2d", [d for d, _f in done])
        self.assertIn("★ step8.4_amset2d", out)
        self.assertIn("旧带隙 2.1000", out)

    def test_plot_gen_wiring(self):
        for step in ("step2_bandgap/step2.2_pbe_plot", "step2_bandgap/step2.3_hse_plot"):
            src = _gen_source(step)
            i = src.index('result["invalidated_downstream"] = _invalidate_gap_consumers(')
            self.assertLess(i, src.index("emit(result, 0)\n\n\nif __name__"))
            self.assertIn("redirect_stdout(sys.stderr)", src)                 # stdout 只留给 JSON

    def test_plot_helper_runs(self):
        import importlib.util
        import os
        p = ROOT / "step2_bandgap" / "step2.3_hse_plot" / "gen_step4.1_plot_band.py"
        spec = importlib.util.spec_from_file_location("_plot41", str(p))
        G = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(G)
        m = _mat()
        _settings(m / "step8_amset", 1.42, "step2.3_hse_plot/band_summary.json（gap_eV）")
        (m / "step8_amset" / "transport.json").write_text("{}")
        cwd0 = os.getcwd()
        os.chdir(m)
        try:
            with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(io.StringIO()) as err:
                got = G._invalidate_gap_consumers(m / "step2_bandgap" / "step2.3_hse_plot", 1.30)
        finally:
            os.chdir(cwd0)
        self.assertIn("step8_amset/transport.json", got)
        self.assertEqual(out.getvalue(), "")
        self.assertIn("1.4200 -> 1.3000", err.getvalue())


class PostLineageTests(unittest.TestCase):
    def test_stale_s81_and_gap_mismatch_reported(self):
        import os
        import time
        m = _mat()
        (m / "step1_opt").mkdir()
        f81 = m / "step8.1_boltztrap" / "boltztrap_crta.json"
        f81.write_text("{}")
        old = time.time() - 3600
        os.utime(f81, (old, old))
        (m / "step8.2_dpt" / "dpt_result.json").write_text("{}")           # S8.2 后来重跑过（GaAs）
        _settings(m / "step8_amset", 1.4200, "step2.3_hse_plot/band_summary.json（gap_eV）")
        (m / "step2_bandgap" / "step2.3_hse_plot" / "band_summary.json").write_text(json.dumps({"gap_eV": 1.30}))
        probs = kc.post_lineage(m)
        self.assertEqual(sorted(s for s, _w in probs), ["step8.1_boltztrap", "step8_amset"])
        self.assertIn("1.4200", dict(probs)["step8_amset"])
        self.assertEqual(kc.structure_lineage(m)[1], [])                  # S8 的闸门不管这些
        sys.path.insert(0, str(ROOT / "tools"))
        import lineage_check as L
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            self.assertEqual(L.main([str(m)]), 1)
        self.assertIn("step8.1_boltztrap", buf.getvalue())
        (m / "step8_amset" / "transport.json").write_text("{}")
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            L.main([str(m), "--invalidate-from", "step2_bandgap/step2.3_hse_plot", "step8.2_dpt"])
        self.assertIn("step8_amset/transport.json", buf.getvalue())
        self.assertIn("step8.1_boltztrap/boltztrap_crta.json", buf.getvalue())


class DptConfTests(unittest.TestCase):
    def setUp(self):
        import importlib
        import gen_step12_dpt as D
        self.D = importlib.reload(D)

    def test_step_conf_override_is_per_material(self):
        D = self.D
        a, b = _mat(), _mat()
        (a / "step.conf").write_text("[params]\nFUNC = pbesol\nM_EFF_ELECTRON = 0.030\nE1_SOURCE = amset\n")
        with contextlib.redirect_stdout(io.StringIO()):
            rec = D.apply_conf(a)
        self.assertEqual(D.MANUAL["m_eff_electron"], 0.030)
        self.assertEqual(D.E1_SOURCE, "amset")
        self.assertEqual(rec["m_eff_electron"]["source"], "step.conf M_EFF_ELECTRON")
        self.assertEqual(D.get_effective_mass(a, "electron", False), (0.030, "manual(step.conf)"))
        self.assertEqual(D._grid_resolution_note("electron", 0.03, "manual(step.conf)", 300.0),
                         "manual(step.conf)")
        self.setUp()                                                        # 另一个材料：干净的脚本
        D = self.D
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(D.apply_conf(b), {})
        self.assertIsNone(D.MANUAL["m_eff_electron"])
        self.assertEqual(D.E1_SOURCE, "vac")

    def test_edited_script_warned(self):
        D = self.D
        D.MANUAL["m_eff_electron"] = 0.030
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            rec = D.apply_conf(_mat())
        self.assertIn("★ 共享脚本", buf.getvalue())
        self.assertIn("作用到所有材料", rec["m_eff_electron"]["source"])
        self.assertIn("脚本 MANUAL", D._manual_prov("m_eff_electron"))

    def test_bad_values_rejected(self):
        D = self.D
        m = _mat()
        (m / "step.conf").write_text("[params]\nE1_SOURCE = core\n")
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
            D.apply_conf(m)
        (m / "step.conf").write_text("[params]\nFORCE_NSTEP = 7\n")
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
            D.apply_conf(m)

    def test_rerun_invalidates_s81_s83(self):
        D = self.D
        m = _mat()
        (m / "step8.1_boltztrap" / "boltztrap_crta.json").write_text("{}")
        (m / "step8.3_output" / "comparison_300K.png").write_text("x")
        with contextlib.redirect_stdout(io.StringIO()):
            done = D.invalidate_consumers(m)
        self.assertEqual(sorted(done), [("step8.1_boltztrap", "boltztrap_crta.json"),
                                        ("step8.3_output", "comparison_300K.png")])
        src = _gen_source("step8.2_dpt")
        i = src.index("    invalidate_consumers(cwd)\n    overrides = apply_conf(cwd)")
        self.assertLess(i, src.index('"results": [_one_carrier('))
        self.assertIn('"overrides": overrides', src)

    def test_main_records_overrides(self):
        D = self.D
        m = _mat()
        (m / "step.conf").write_text("[params]\nM_EFF_HOLE = 0.5\nE1_HOLE_EV = 5\nC_3D_GPA = 100\n")
        import os
        cwd0 = os.getcwd()
        os.chdir(m)
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                try:
                    D.main()
                except SystemExit:
                    pass
        finally:
            os.chdir(cwd0)
        j = json.loads((m / "step8.2_dpt" / "dpt_result.json").read_text(encoding="utf-8"))
        self.assertEqual(j["overrides"]["m_eff_hole"]["value"], 0.5)
        hole = next(r for r in j["results"] if r["carrier"] == "hole")
        self.assertEqual(hole["inputs"]["m_provenance"], "manual(step.conf)")
        self.assertEqual(hole["inputs"]["E1_provenance"], "manual(step.conf)")
        self.assertEqual(hole["inputs"]["C_provenance"], "manual(step.conf)")
        self.assertIsNotNone(hole["mobility_cm2_Vs"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
