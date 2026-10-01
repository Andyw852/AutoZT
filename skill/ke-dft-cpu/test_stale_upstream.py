#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_stale_upstream.py —— V136：上游重算后下游结果的失效与核对。不需要 VASP / AMSET。

用法：python test_stale_upstream.py      退出码 0 = 全部 PASS。

1) ke_common.invalidate_downstream：S5（介电）重算 -> S5.1 校验、S8、S8.4 及其下游 S8.3 的完成标记归档；
   S6（弹性）重算 -> S8、S8.4、S8.2 归档。S5/S6 的 gen 在已有 OUTCAR 时调用它。
2) validate_dielectric.py 写 status / input_files（OUTCAR 的 sha256，路径相对材料目录）；
   autozt 的 strict_done_marker（ck_plot）：通过 -> 完成；ok=false -> 不算完成；S5 的 OUTCAR 变了 -> 作废。
3) tools/regate_projects.py：transport.json 比 step6_elastic/OUTCAR、形变势 h5 等旧 -> ⚠ STALE；
   S5.1 校验 ok=false -> ★ WRONG。
"""
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "tools"), str(ROOT / "step8.4_amset2d")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402

GOOD = """ MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)
 ------------------------------------------------------
            10.500     0.000     0.000
             0.000    10.500     0.000
             0.000     0.000    10.500
 ------------------------------------------------------
 MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION
 ------------------------------------------------------
             2.000     0.000     0.000
             0.000     2.000     0.000
             0.000     0.000     2.000
 ------------------------------------------------------
"""
BAD = GOOD.replace("10.500", "NaN").replace("2.000", "NaN")
POSCAR = """GaAs
1.0
0.0 2.83 2.83
2.83 0.0 2.83
2.83 2.83 0.0
Ga As
1 1
Direct
0.0 0.0 0.0
0.25 0.25 0.25
"""


def _touch(p, text="x", mtime=None):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    if mtime is not None:
        os.utime(p, (mtime, mtime))


def _collector():
    """只取 autozt/_collector_remote.py 里的 ck_plot（该文件是远端脚本，import 时会直接跑 main）。"""
    import ast
    import glob
    import types
    f = REPO / "autozt" / "_collector_remote.py"
    if not f.is_file():
        return None
    tree = ast.parse(f.read_text(encoding="utf-8"))
    fn = next((n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "ck_plot"), None)
    if fn is None:
        return None
    ns = {"json": json, "os": os, "glob": glob}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(f), "exec"), ns)
    return types.SimpleNamespace(ck_plot=ns["ck_plot"])


class InvalidateTests(unittest.TestCase):
    def _mat(self):
        d = Path(tempfile.mkdtemp())
        for f in ("step5_dielect_validate/dielectric_check.json", "step8_amset/transport.json",
                  "step8.4_amset2d/transport.json", "step8.4_amset2d/intrinsic_transport.json",
                  "step8.2_dpt/dpt_result.json", "step8.3_output/comparison_300K.png"):
            _touch(d / f)
        return d

    def test_s5_regen(self):
        d = self._mat()
        done = dict(kc.invalidate_downstream(d, "step5_dielect", "test"))
        for step in ("step5_dielect_validate", "step8_amset", "step8.4_amset2d", "step8.3_output"):
            self.assertIn(step, done, done)
        self.assertTrue((d / "step8.2_dpt" / "dpt_result.json").is_file())      # DPT 不用介电
        self.assertFalse((d / "step8_amset" / "transport.json").exists())

    def test_s6_regen(self):
        d = self._mat()
        done = dict(kc.invalidate_downstream(d, "step6_elastic", "test"))
        for step in ("step8_amset", "step8.4_amset2d", "step8.2_dpt"):
            self.assertIn(step, done, done)
        self.assertTrue((d / "step5_dielect_validate" / "dielectric_check.json").is_file())

    def test_gens_call_it(self):
        s5 = (ROOT / "step5_dielect" / "gen_step8_dielect.py").read_text(encoding="utf-8")
        i = s5.index('if (out / "OUTCAR").is_file():\n        kc.invalidate_downstream(cwd, OUTDIR_NAME')
        self.assertLess(i, s5.index("kc.relay_poscar(prev / \"CONTCAR\""))
        s6 = (ROOT / "step6_elastic" / "gen_step2_elastic.py").read_text(encoding="utf-8")
        i = s6.index('_kc.invalidate_downstream(Path.cwd(), STEP2_DIR')
        self.assertLess(i, s6.index('Path(step2 / "POSCAR").write_text('))
        for y in (ROOT / "skill.yaml", ROOT.parent / "zt-dft-cpu" / "skill.yaml"):
            t = y.read_text(encoding="utf-8")
            j = t.index("gen: gen_step2_elastic.py,")
            self.assertIn("ke_common.py", t[j:j + 200], y)
            k = t.index("name: step5_dielect_validate")
            self.assertIn("strict_done_marker: true", t[k:k + 300], y)


class StrictDievalidTests(unittest.TestCase):
    def _run(self, outcar_text):
        d = Path(tempfile.mkdtemp())
        _touch(d / "step5_dielect" / "OUTCAR", outcar_text)
        _touch(d / "step5_dielect" / "POSCAR", POSCAR)
        r = subprocess.run([sys.executable, str(ROOT / "step5_dielect" / "validate_dielectric.py"),
                            "--json", "step5_dielect_validate/dielectric_check.json"],
                           cwd=str(d), capture_output=True, text=True, timeout=120)
        js = json.loads((d / "step5_dielect_validate" / "dielectric_check.json").read_text())
        return d, r, js

    def test_marker_fields(self):
        d, r, js = self._run(GOOD)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(js["ok"])
        self.assertEqual(js["status"], "done")
        self.assertEqual(list(js["input_files"]), ["step5_dielect/OUTCAR"])
        d, r, js = self._run(BAD)
        self.assertEqual(r.returncode, 40)
        self.assertFalse(js["ok"])
        self.assertEqual(js["status"], "failed")

    def test_collector_strict(self):
        C = _collector()
        if C is None:
            self.skipTest("仓库里没有 autozt/_collector_remote.py")
        sc = {"done_marker": "dielectric_check.json", "strict_done_marker": True}
        d, _, _ = self._run(GOOD)
        ok, why = C.ck_plot(str(d / "step5_dielect_validate"), sc)
        self.assertTrue(ok, why)
        (d / "step5_dielect" / "OUTCAR").write_text(GOOD.replace("10.500", "10.600"))   # S5 重算了
        ok, why = C.ck_plot(str(d / "step5_dielect_validate"), sc)
        self.assertFalse(ok)
        self.assertIn("stale completion marker", why)
        d, _, _ = self._run(BAD)
        ok, why = C.ck_plot(str(d / "step5_dielect_validate"), sc)
        self.assertFalse(ok, why)
        ok, _ = C.ck_plot(str(d / "step5_dielect_validate"), {"done_marker": "dielectric_check.json"})
        self.assertTrue(ok)                                  # 旧判据（只看文件在不在）会放行 —— 就是这个 bug


class RegateStaleTests(unittest.TestCase):
    def test_audit(self):
        import yaml
        import regate_projects as R
        root = Path(tempfile.mkdtemp())
        now = time.time()
        cases = {"fresh": dict(t_tr=now, t_el=now - 3600, dfail=False, verdict="OK"),
                 "stale_el": dict(t_tr=now - 7200, t_el=now, dfail=False, verdict="⚠ STALE"),
                 "dfail": dict(t_tr=now, t_el=now - 3600, dfail=True, verdict="★ WRONG")}
        for name, c in cases.items():
            mat = root / name / "ke-dft-cpu"
            _touch(mat / "step3_uniform" / "POSCAR", POSCAR)
            _touch(mat / "step6_elastic" / "OUTCAR", "x", mtime=c["t_el"])
            run = mat / "step8_amset"
            _touch(run / "settings.yaml", yaml.safe_dump({"unity_overlap": True}))
            _touch(run / "transport.json", "{}", mtime=c["t_tr"])
            js = {"ok": not c["dfail"], "reasons": ["eps_inf 含 NaN/Inf"] if c["dfail"] else []}
            _touch(mat / "step5_dielect_validate" / "dielectric_check.json", json.dumps(js))
        rows = {r["material"].split("/")[0]: r for r in R.audit(str(root), "*/ke-dft-cpu")}
        for name, c in cases.items():
            self.assertEqual(rows[name]["verdict"], c["verdict"], rows[name]["why"])
        self.assertIn("step6_elastic/OUTCAR", rows["stale_el"]["stale_upstream"])
        self.assertIn("NaN", rows["dfail"]["why"])

    def test_symlinked_h5_newer(self):
        import regate_projects as R
        mat = Path(tempfile.mkdtemp())
        _touch(mat / "step7b_deform_read" / "deformation.h5", "x", mtime=time.time())
        run = mat / "step8_amset"
        _touch(run / "transport.json", "{}", mtime=time.time() - 7200)
        (run / "deformation.h5").symlink_to("../step7b_deform_read/deformation.h5")
        self.assertEqual(R.stale_upstreams(mat, run), ["形变势 h5"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
