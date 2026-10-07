#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_functional_label.py —— 结果文件里写明泛函（V181）；lineage_check 按泛函分组列材料（V182）。

用法：python test_functional_label.py     退出码 0 = 全部 PASS。

起因（2026-10-07）：出厂泛函是 PBEsol（step.conf FUNC = pbesol -> GGA=PS），但 dpt_result.json、
2d_correction.json、8.3 汇总都不记泛函。agent 把 INCAR 的 GGA=PS 报成了"PBE"，和 PBE 文献逐项比时，
PBEsol 晶格更小带来的 C_2D 偏硬差点被归到别的原因上。材料目录在测试里构造，不读真实文件。
"""
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(HERE.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


def _load(rel, name):
    spec = importlib.util.spec_from_file_location(name, str(HERE / rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _mat(method=None, incar=None, step="step1_opt"):
    d = Path(tempfile.mkdtemp())
    (d / step).mkdir()
    if method is not None:
        (d / step / "workflow_method.txt").write_text(method)
    if incar is not None:
        (d / step / "INCAR").write_text(incar)
    return d


class MaterialFunctionalTests(unittest.TestCase):
    def test_method_file_pbesol(self):
        f = kc.material_functional(_mat("DIM=2D\nFUNC=pbesol\n# note\n"))
        self.assertEqual(f["func"], "pbesol")
        self.assertTrue(f["label"].startswith("PBEsol"))
        self.assertIn("GGA=PS", f["label"])
        self.assertEqual(f["source"], "step1_opt/workflow_method.txt")
        self.assertIn("不是 PBE", f["note"])

    def test_incar_fallback(self):
        f = kc.material_functional(_mat("DIM=2D\n", "SYSTEM = x\nGGA    = PS\nISIF = 3\n"))
        self.assertEqual((f["func"], f["source"]), ("pbesol", "step1_opt/INCAR（从 GGA/IVDW 反推）"))
        f = kc.material_functional(_mat(None, "GGA = PE\nIVDW = 12\n", step="step1"))
        self.assertEqual((f["func"], f["source"]), ("pbe-d3", "step1/INCAR（从 GGA/IVDW 反推）"))
        self.assertTrue(f["label"].startswith("PBE+D3"))

    def test_method_file_wins_over_incar(self):
        f = kc.material_functional(_mat("FUNC=pbe\n", "GGA = PS\n"))
        self.assertEqual(f["func"], "pbe")

    def test_unknown(self):
        f = kc.material_functional(_mat())
        self.assertIsNone(f["func"])
        self.assertIn("未知", f["label"])

    def test_labels_cover_supported(self):
        self.assertEqual(set(kc.FUNC_LABEL), set(kc.SUPPORTED_FUNCS))
        self.assertEqual(set(kc.FUNC_COMPARE_NOTE), set(kc.SUPPORTED_FUNCS))


class S83Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = _load("step8.3_output/gen_step13_output.py", "gen_step13_output_v181")

    def test_labels_match_ke_common(self):
        self.assertEqual(self.m._FUNC_LABEL, kc.FUNC_LABEL)

    def test_reads_s82_result_first(self):
        d = _mat("FUNC=pbe\n")
        (d / "step8.2_dpt").mkdir()
        (d / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(
            {"functional": kc.material_functional(_mat("FUNC=pbesol\n"))}))
        lab, src = self.m._functional(d)
        self.assertTrue(lab.startswith("PBEsol"))
        self.assertEqual(src, "step8.2_dpt/dpt_result.json")

    def test_falls_back_to_s84_then_method(self):
        d = _mat("FUNC=pbesol\n")
        lab, src = self.m._functional(d)
        self.assertEqual((lab, src), (kc.FUNC_LABEL["pbesol"], "step1_opt/workflow_method.txt"))
        (d / "step8.4_amset2d").mkdir()
        (d / "step8.4_amset2d" / "2d_correction.json").write_text(json.dumps(
            {"functional": {"func": "pbe", "label": kc.FUNC_LABEL["pbe"]}}))
        self.assertEqual(self.m._functional(d)[1], "step8.4_amset2d/2d_correction.json")
        self.assertEqual(self.m._functional(_mat()), (None, None))

    def test_summary_and_csv_carry_functional(self):
        d = _mat("FUNC=pbesol\n")
        out = d / "step8.3_output"
        out.mkdir()
        old = os.getcwd()
        os.chdir(d)
        try:
            self.m.write_table(out, [{"carrier_conc_cm-3": -1e19, "type": "n"}], None, None, None)
        finally:
            os.chdir(old)
        summ = (out / "comparison_summary.txt").read_text(encoding="utf-8")
        self.assertIn("[口径] 泛函 = PBEsol（GGA=PS，无色散修正）", summ)
        self.assertIn("不是 PBE", summ)
        head, row = (out / "comparison_300K.csv").read_text(encoding="utf-8").splitlines()[:2]
        self.assertIn("functional", head.split(","))
        self.assertIn("PBEsol", row)

    def test_summary_says_unknown(self):
        d = _mat()
        out = d / "step8.3_output"
        out.mkdir()
        old = os.getcwd()
        os.chdir(d)
        try:
            self.m.write_table(out, [], None, None, None)
        finally:
            os.chdir(old)
        self.assertIn("[口径] 泛函 = 未知", (out / "comparison_summary.txt").read_text(encoding="utf-8"))


class WiringTests(unittest.TestCase):
    def test_s82_writes_functional(self):
        src = (HERE / "step8.2_dpt" / "gen_step12_dpt.py").read_text(encoding="utf-8")
        self.assertIn('"functional": functional,', src)
        self.assertIn("functional = _functional(cwd)", src)
        self.assertIn('lines.append("# 泛函：%s —— %s"', src)
        m = _load("step8.2_dpt/gen_step12_dpt.py", "gen_step12_dpt_v181")
        f = m._functional(_mat("FUNC=pbesol\n"))
        self.assertEqual(f["func"], "pbesol")

    def test_s84_record_has_functional(self):
        src = (HERE / "step8.4_amset2d" / "gen_step14_amset2d.py").read_text(encoding="utf-8")
        i = src.index("    rec = {\n        \"cell_c_A\": round(c_len, 4),")
        self.assertIn('"functional": _func,', src[i:i + 200])
        self.assertIn("kc.material_functional(cwd)", src[i - 300:i])


class LineageReportTests(unittest.TestCase):
    """V182：lineage_check 末尾按泛函分组列材料；不止一种时提示，不改退出码。"""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(HERE / "tools"))
        import lineage_check
        cls.L = lineage_check

    def _project(self, funcs):
        root = Path(tempfile.mkdtemp())
        for name, func in funcs.items():
            d = root / name / "ke-dft-cpu" / "step1_opt"
            d.mkdir(parents=True)
            (d / "workflow_method.txt").write_text("DIM=2D\nFUNC=%s\n" % func)
        return root

    def _run(self, root):
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            rc = self.L.main([str(root)])
        return rc, buf.getvalue()

    def test_mixed_functionals_listed_and_flagged(self):
        rc, out = self._run(self._project({"CrA": "pbe", "CrB": "pbe", "MoA": "pbesol"}))
        self.assertEqual(rc, 0, out)                                      # 只是提示
        self.assertIn("PBE（GGA=PE，无色散修正） × 2：CrA、CrB", out)
        self.assertIn("PBEsol（GGA=PS，无色散修正） × 1：MoA", out)
        self.assertIn("[提示] 这些材料不是同一个泛函", out)

    def test_single_functional_no_flag(self):
        rc, out = self._run(self._project({"MoA": "pbesol", "WA": "pbesol"}))
        self.assertIn("PBEsol（GGA=PS，无色散修正） × 2：MoA、WA", out)
        self.assertNotIn("[提示]", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
