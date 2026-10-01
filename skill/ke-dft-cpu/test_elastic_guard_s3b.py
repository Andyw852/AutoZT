#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_elastic_guard_s3b.py —— V134：弹性张量力学稳定性检查 + S3b 从 S3 的电荷密度起步。

用法：python test_elastic_guard_s3b.py      退出码 0 = 全部 PASS。不需要 VASP / AMSET。

1) ke_common.check_elastic_stability：
   · Mo2S3 / P1_Mo-MoS2 用户侧 settings.yaml 的张量（2D）：面内块正定，但 C66 只有 C11 的 0.7% -> 顺序告警；
   · 同一张量按"第 4 位其实是 XY"重排后：Mo2S3 的面内剪切 −40.16 -> 不正定 -> 退出；P1_Mo-MoS2 正定；
   · ELASTIC_GUARD 关 / 散射不含 ADP -> 只告警；3D 正定/不正定；标量。
2) S8 / S8.4 的 gen：ELASTIC_GUARD 进 SPEC、main 里 global 并读 step.conf、检查在 2D 处理之后、write_settings 之前。
3) gen_step5b_uniform_full.s3_density_plan：S3 完成 + 结构/INCAR 物理键一致 -> 可以从 S3 起步；
   缺 CHGCAR、结构不同、ENCUT 不同、ICHARG 被显式改过 -> 不行；只差并行/对称性/输出键 -> 可以。
   insert_before_launch 把复制命令插在第一条 mpirun/srun 前；_set_incar 改键保留其余行。
4) tools/regate_projects.py 扫已有运行目录：弹性张量不正定 -> ★ WRONG；2D 面内剪切近零 -> ⚠ ELASTIC。
"""
import contextlib
import io
import re
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for _p in (str(ROOT), str(ROOT.parent / "_common" / "opt"), str(ROOT / "step3b_uniform_full")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402

MO2S3 = [[155.653, 51.681, -0.202, -0.217, -0.003, 0.041], [51.681, 155.216, 0.050, -0.200, -0.009, -0.026],
         [-0.202, 0.050, 0.681, -0.408, -0.044, 0.068], [-0.217, -0.200, -0.408, -40.158, 0.032, -0.038],
         [-0.003, -0.009, -0.044, 0.032, 4.169, 0.026], [0.041, -0.026, 0.068, -0.038, 0.026, 1.021]]
P1_MOS2 = [[191.41, 45.154, 0.887, 0, -0.003, 0], [45.154, 171.939, 0.905, 0, 0.009, 0],
           [0.887, 0.905, 3.819, 0, 0.009, 0], [0, 0, 0, 46.035, 0, 0.015],
           [-0.003, 0.009, 0.009, 0, 3.939, 0], [0, 0, 0, 0.015, 0, 1.312]]
SI = [[153.0, 57.0, 57.0, 0, 0, 0], [57.0, 153.0, 57.0, 0, 0, 0], [57.0, 57.0, 153.0, 0, 0, 0],
      [0, 0, 0, 75.0, 0, 0], [0, 0, 0, 0, 75.0, 0], [0, 0, 0, 0, 0, 75.0]]


def _vasp_slot4_is_xy(m):
    p = (0, 1, 2, 4, 5, 3)
    return [[m[i][j] for j in p] for i in p]


def _quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = fn(*a, **k)
    return r, buf.getvalue()


class ElasticGuardTests(unittest.TestCase):
    def test_user_tensors_as_written_warn_order(self):
        for m, hint in ((MO2S3, True), (P1_MOS2, True)):
            r, out = _quiet(kc.check_elastic_stability, m, is_2d=True, label="S8")
            self.assertTrue(r["ok"])
            self.assertTrue(r["soft_c66"])
            self.assertEqual(r["misorder_hint"], hint)
            self.assertIn("已重排成标准 Voigt", out)

    def test_reordered_mo2s3_exits(self):
        with self.assertRaises(SystemExit) as cm:
            _quiet(kc.check_elastic_stability, _vasp_slot4_is_xy(MO2S3), is_2d=True, label="S8")
        self.assertIn("不正定", str(cm.exception))
        self.assertIn("-40.158", str(cm.exception))
        self.assertIn("ELASTIC_GUARD = false", str(cm.exception))
        r, _ = _quiet(kc.check_elastic_stability, _vasp_slot4_is_xy(P1_MOS2), is_2d=True)
        self.assertTrue(r["ok"])
        self.assertFalse(r["soft_c66"])

    def test_disabled_or_unused_only_warns(self):
        bad = _vasp_slot4_is_xy(MO2S3)
        r, out = _quiet(kc.check_elastic_stability, bad, is_2d=True, enabled=False)
        self.assertFalse(r["ok"])
        self.assertIn("[WARN]", out)
        r, out = _quiet(kc.check_elastic_stability, bad, is_2d=True, used=False)
        self.assertFalse(r["ok"])
        self.assertIn("不含 ADP", out)

    def test_3d(self):
        r, _ = _quiet(kc.check_elastic_stability, SI, is_2d=False)
        self.assertTrue(r["ok"])
        with self.assertRaises(SystemExit):
            _quiet(kc.check_elastic_stability, MO2S3, is_2d=False)       # 当 3D 看：C44 < 0
        near = [row[:] for row in SI]
        near[3][3] = 0.01
        r, out = _quiet(kc.check_elastic_stability, near, is_2d=False)
        self.assertTrue(r["near_singular"])
        self.assertIn("接近 0", out)

    def test_scalar_and_none(self):
        self.assertIsNone(kc.check_elastic_stability(None))
        r, _ = _quiet(kc.check_elastic_stability, 80.0)
        self.assertTrue(r["ok"])
        with self.assertRaises(SystemExit):
            _quiet(kc.check_elastic_stability, -5.0)


class GenWiringTests(unittest.TestCase):
    def test_s8_and_s84(self):
        for rel in ("step8_amset/gen_step10_amset.py", "step8.4_amset2d/gen_step14_amset2d.py"):
            src = (ROOT / rel).read_text(encoding="utf-8")
            self.assertIn('"ELASTIC_GUARD": (True, "bool")', src, rel)
            self.assertRegex(src, r"global [^\n]*ELASTIC_GUARD", rel)
            self.assertIn('ELASTIC_GUARD = bool(_p["ELASTIC_GUARD"])', src, rel)
            i_2d = src.index("is_2d, elastic, c_len = apply_2d_corrections(")
            i_chk = src.index("kc.check_elastic_stability(elastic, is_2d=bool(is_2d), enabled=ELASTIC_GUARD,")
            i_ws = src.index("    write_settings(out, eps_inf, eps_static, gap, elastic,")
            self.assertLess(i_2d, i_chk, rel)
            self.assertLess(i_chk, i_ws, rel)


POSCAR = """MoSe2
1.0
3.32 0.0 0.0
-1.66 2.875204 0.0
0.0 0.0 20.0
Mo Se
1 2
Direct
0.0 0.0 0.5
0.333333 0.666667 0.58
0.333333 0.666667 0.42
"""
INCAR_TPL = (ROOT / "step3b_uniform_full" / "incar_uniform_full_2d.tpl").read_text(encoding="utf-8")


def _material(chg=True, finished=True, poscar_s3b=POSCAR, s3b_edits=None):
    import gen_step5b_uniform_full as G
    d = Path(tempfile.mkdtemp())
    s3, s3b = d / "step3_uniform", d / "step3b_uniform_full"
    s3.mkdir()
    s3b.mkdir()
    sub = {"{{SYSTEM}}": "MoSe2 uniform", "{{ENCUT}}": "520", "{{GGA}}": "PE", "{{VDW_LINE}}": "IVDW = 12"}
    inc = INCAR_TPL
    for k, v in sub.items():
        inc = inc.replace(k, v)
    (s3 / "INCAR").write_text(inc.replace("ISYM   = -1", "ISYM   = 2"))
    b = inc
    for old, new in (s3b_edits or {}).items():
        b = b.replace(old, new)
    (s3b / "INCAR").write_text(b)
    (s3 / "POSCAR").write_text(POSCAR)
    (s3b / "POSCAR").write_text(poscar_s3b)
    if chg:
        (s3 / "CHGCAR").write_text("x" * 4096)
    (s3 / "OUTCAR").write_text("   NBANDS=     24\n" + ("General timing and accounting\n" if finished else ""))
    return G, d, s3b


class S3bFromS3Tests(unittest.TestCase):
    def test_ok(self):
        G, d, out = _material(s3b_edits={"KPAR   = 2": "KPAR   = 4"})       # 并行键不同不算
        ok, why, nb = G.s3_density_plan(d, out)
        self.assertTrue(ok, why)
        self.assertEqual(nb, 24)

    def test_refusals(self):
        cases = [(dict(chg=False), "CHGCAR"),
                 (dict(finished=False), "正常结束"),
                 (dict(poscar_s3b=POSCAR.replace("0.58", "0.59")), "POSCAR"),
                 (dict(s3b_edits={"ENCUT  = 520": "ENCUT  = 600"}), "ENCUT"),
                 (dict(s3b_edits={"ICHARG = 2": "ICHARG = 11"}), "ICHARG")]
        for kw, word in cases:
            G, d, out = _material(**kw)
            ok, why, _ = G.s3_density_plan(d, out)
            self.assertFalse(ok, kw)
            self.assertIn(word, why, kw)

    def test_submit_and_incar_patch(self):
        import gen_step5b_uniform_full as G
        tpl = "#!/bin/bash\n#SBATCH -N 1\nmodule load vasp\nmpirun -np $SLURM_NTASKS vasp_std\n"
        out = G.insert_before_launch(tpl, G.S3_CHG_COPY)
        self.assertLess(out.index("cp -f ../step3_uniform/CHGCAR CHGCAR"), out.index("mpirun"))
        self.assertIn("exit 1", out)
        self.assertIsNone(G.insert_before_launch("#!/bin/bash\necho hi\n", G.S3_CHG_COPY))
        self.assertIsNotNone(G.insert_before_launch("srun vasp_std\n", G.S3_CHG_COPY))
        _, d, s3b = _material()
        G._set_incar(s3b / "INCAR", {"ICHARG": "1", "NBANDS": "24"}, "V134")
        inc = kc.parse_incar((s3b / "INCAR").read_text())
        self.assertEqual(inc["ICHARG"], "1")
        self.assertEqual(inc["NBANDS"], "24")
        self.assertEqual(inc["ISYM"], "-1")
        self.assertEqual(len(re.findall(r"^ICHARG", (s3b / "INCAR").read_text(), re.M)), 1)

    def test_spec_and_skill_yaml(self):
        src = (ROOT / "step3b_uniform_full" / "gen_step5b_uniform_full.py").read_text(encoding="utf-8")
        self.assertIn('"START_FROM_S3": (START_FROM_S3, "bool")', src)
        self.assertIn('START_FROM_S3 = bool(_conf["START_FROM_S3"])', src)
        for y in (ROOT / "skill.yaml", ROOT.parent / "zt-dft-cpu" / "skill.yaml"):
            t = y.read_text(encoding="utf-8")
            i = t.index('name: "step3b_uniform_full"')
            self.assertIn('needs: ["step3_uniform"]', t[i:i + 800], y)


class RegateElasticTests(unittest.TestCase):
    """V134：regate_projects 从 settings.yaml 读弹性张量：不正定 -> ★ WRONG；2D 面内剪切近零 -> ⚠ ELASTIC。"""

    def test_audit(self):
        import yaml
        sys.path.insert(0, str(ROOT / "tools"))
        sys.path.insert(0, str(ROOT / "step8.4_amset2d"))
        import regate_projects as R
        root = Path(tempfile.mkdtemp())
        cases = {"as_written": (MO2S3, "⚠ ELASTIC"), "reordered": (_vasp_slot4_is_xy(MO2S3), "★ WRONG"),
                 "p1_fixed": (_vasp_slot4_is_xy(P1_MOS2), "OK"), "no_adp": (_vasp_slot4_is_xy(MO2S3), "OK")}
        for name, (el, _) in cases.items():
            mat = root / name / "ke-dft-cpu"
            (mat / "step3_uniform").mkdir(parents=True)
            (mat / "step3_uniform" / "POSCAR").write_text(POSCAR)
            run = mat / "step8.4_amset2d"
            run.mkdir()
            (run / "2d_correction.json").write_text("{}")
            (run / "settings.yaml").write_text(yaml.safe_dump({
                "unity_overlap": True, "elastic_constant": el,
                "scattering_type": ["IMP"] if name == "no_adp" else ["ADP", "IMP"]}))
        rows = {r["material"].split("/")[0]: r for r in R.audit(str(root), "*/ke-dft-cpu")}
        for name, (_, verdict) in cases.items():
            self.assertEqual(rows[name]["verdict"], verdict, rows[name]["why"])
        self.assertIn("疑似 VASP 顺序没重排", rows["as_written"]["why"])
        self.assertIn("-40.16", rows["reordered"]["why"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
