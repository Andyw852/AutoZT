# -*- coding: utf-8 -*-
"""[patch_ionrelax_locpot] ionrelax 缺 LOCPOT（纯本地、无数据）。

  · S7.1：本征值取自 ionrelax/vasprun.xml 时 LOCPOT 只能取 ionrelax/LOCPOT；缺了抛 RuntimeError
    （外层只把 OSError 当"合法跳过"），刚性构型的 LOCPOT 不会被用上；真空口径 h5 生不出来时报错、不退回 core；
    DEFORM_REF 可由 step.conf 显式写 core。
  · S7：只补静态段的触发矩阵；改名归档后 ck_deform 判未完成；submit.sh 整体重写成只跑静态段；
    INCAR.static 改 ISTART=0 / ICHARG=2，其余不动；补完后再 retry 不再触发、也不整个重建。
"""
import ast
import os
import re
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SKILL = REPO / "skill" / "ke-dft-cpu"
S7 = SKILL / "step7_deform"
S71 = S7 / "step7b_read" / "gen_step9b_deform_read.py"
for _p in (str(REPO), str(SKILL), str(S7)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

DONE = " General timing and accounting informations for this job:\n"


def _gen9():
    # 模块级要解析 AMSET 环境（只拼字符串，不执行）；测试里给个假的，免得 SystemExit
    os.environ.setdefault("AUTOZT_CONDA_SH", "/dev/null")
    os.environ.setdefault("AUTOZT_AMSET_ENV", "amset-test")
    import gen_step9_deform
    return gen_step9_deform


def _collector():
    src = (REPO / "autozt" / "_collector_remote.py").read_text(encoding="utf-8")
    C = types.ModuleType("_collector_defs")
    exec(compile(src[:src.index("\ndef main():")], "_collector_remote.py", "exec"), C.__dict__)
    return C


def _embedded(var):
    """S7.1 里嵌入脚本 _be_code / _vac_code 的源码文本。"""
    s = S71.read_text(encoding="utf-8")
    i = s.index("%s = r'''" % var) + len("%s = r'''" % var)
    return s[i:s.index("'''", i)]


def _func_from(code, name):
    """嵌入脚本里（含 try 块内）名为 name 的函数定义。"""
    for n in ast.walk(ast.parse(code)):
        if isinstance(n, ast.FunctionDef) and n.name == name:
            return n
    raise AssertionError("没找到 %s" % name)


def _exec_func(fn, ns):
    exec(compile(ast.Module(body=[fn], type_ignores=[]), "embedded", "exec"), ns)
    return ns[fn.name]


SUBMIT_TPL = """#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --job-name={{JOBNAME}}
#SBATCH --ntasks-per-node=24
cd $SLURM_SUBMIT_DIR
# 注释里也有 mpirun 字样：OpenMPI 的 openib 会让 mpirun 段错误
module purge
export OMP_NUM_THREADS=1
mpirun -np $SLURM_NTASKS /opt/vasp/bin/vasp_std
"""

INCAR_MAIN = """SYSTEM = MoS2 deform deform-01
ISTART = 0
ICHARG = 2
PREC   = Accurate
ENCUT  = 520
EDIFF  = 1E-8          # 带边能量差
IBRION = -1
NSW    = 0
ISYM   = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.
LVHAR  = .TRUE.
"""

VASPRUN = """<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <incar>
  <i type="string" name="SYSTEM">MoS2 deform deform-01</i>
  <i type="int" name="ISTART">     %(istart)s</i>
  <i type="int" name="ICHARG">     %(icharg)s</i>
  <i type="string" name="PREC">accurate</i>
  <i name="ENCUT">    520.00000000</i>
  <i name="EDIFF">      0.00000001</i>
  <i type="int" name="IBRION">    -1</i>
  <i type="int" name="NSW">     0</i>
  <i type="int" name="ISYM">     0</i>
  <i type="logical" name="LWAVE"> F  </i>
  <i type="logical" name="LCHARG"> F  </i>
  <i type="logical" name="LVHAR"> T  </i>
 </incar>
 <parameters/>
 <calculation/>
</modeling>
"""


def _finished_deform(root, name, g, locpot=False):
    """一个 ionrelax 已算完的面内形变目录（用真的 _build_ionrelax 建输入，再摆上产物）。"""
    d = Path(root) / "step7_deform" / name
    d.mkdir(parents=True)
    (d / "POSCAR").write_text("x\n1.0\n 3 0 0\n 0 3 0\n 0 0 20\nMo\n1\nDirect\n 0 0 0.5\n")
    (d / "KPOINTS").write_text("auto\n0\nGamma\n 48 48 1\n")
    (d / "POTCAR").write_text("potcar\n")
    (d / "INCAR").write_text(INCAR_MAIN)
    (d / "submit.sh").write_text(SUBMIT_TPL.replace("{{JOBNAME}}", "MoS2-S7-" + name))
    g._build_ionrelax(d, 520, {}, SUBMIT_TPL)
    (d / "OUTCAR").write_text("rigid run\n" + DONE)              # 刚性单点
    (d / "vasprun.xml").write_text("rigid\n")
    (d / "LOCPOT").write_text("rigid locpot\n")
    ir = d / "ionrelax"
    for nm in ("OUTCAR", "OSZICAR", "CONTCAR", "vasprun.xml"):
        (ir / (nm + ".relax")).write_text("relax %s\n" % nm)
    (ir / "OUTCAR").write_text("static\n" + DONE)
    (ir / "OSZICAR").write_text("static\n")
    (ir / "CONTCAR").write_text("relaxed\n")
    (ir / "vasprun.xml").write_text(VASPRUN % {"istart": 1, "icharg": 0})
    if locpot:
        (ir / "LOCPOT").write_text("ionrelax locpot\n")
    return d


def _launchers(text):
    return [ln.strip() for ln in text.splitlines()
            if ln.strip() and not ln.strip().startswith("#")
            and ln.strip().startswith(("mpirun", "srun", "mpiexec"))]


# ============================== S7.1 读取端 ==============================
class S71ResolvedTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.root = Path(self.td.name)
        self.paths = {}
        self.resolved = _exec_func(_func_from(_embedded("_be_code"), "_resolved"),
                                   {"_os": os, "_resolved_paths": self.paths})

    def _mk(self, name, files):
        d = self.root / name
        for f in files:
            (d / f).parent.mkdir(parents=True, exist_ok=True)
            (d / f).write_text("x\n")
        return str(d)

    def test_ionrelax_without_locpot_raises_runtime_error(self):
        d = self._mk("deform-01", ["vasprun.xml", "LOCPOT", "POSCAR",
                                   "ionrelax/vasprun.xml", "ionrelax/POSCAR"])
        self.assertTrue(self.resolved(d, "vasprun.xml").endswith("ionrelax/vasprun.xml"))
        with self.assertRaises(RuntimeError) as cm:
            self.resolved(d, "LOCPOT")
        self.assertNotIsInstance(cm.exception, OSError)          # 外层 except OSError 接不住
        self.assertIn("retry S7_deform 补静态段", str(cm.exception))
        self.assertNotIn("不平坦", str(cm.exception))            # 窗口扫描只吞"不平坦"
        self.assertNotIn("LOCPOT", self.paths[d])                # 刚性构型的 LOCPOT 没被登记、没被返回

    def test_ionrelax_with_locpot_uses_it(self):
        d = self._mk("deform-01", ["LOCPOT", "ionrelax/vasprun.xml", "ionrelax/LOCPOT"])
        self.assertEqual(self.resolved(d, "LOCPOT"), os.path.join(d, "ionrelax", "LOCPOT"))
        self.assertEqual(self.paths[d]["LOCPOT"], "ionrelax")

    def test_no_ionrelax_reads_rigid(self):
        d = self._mk("undeformed", ["vasprun.xml", "LOCPOT"])
        self.assertEqual(self.resolved(d, "LOCPOT"), os.path.join(d, "LOCPOT"))
        self.assertEqual(self.paths[d]["LOCPOT"], "rigid")
        # ionrelax/ 建了但还没有 vasprun.xml：本征值也走本级，同一构型
        d2 = self._mk("deform-02", ["vasprun.xml", "LOCPOT", "ionrelax/POSCAR", "ionrelax/LOCPOT"])
        self.assertEqual(self.resolved(d2, "LOCPOT"), os.path.join(d2, "LOCPOT"))

    def test_runtime_error_not_swallowed_by_outer_handlers(self):
        tree = ast.parse(_embedded("_be_code"))
        caught = set()
        for t in (n for n in tree.body if isinstance(n, ast.Try)):     # 真空对齐那个 try（带 _VacSkip）
            names = {getattr(e, "id", None) for h in t.handlers
                     for e in (h.type.elts if isinstance(h.type, ast.Tuple) else [h.type])}
            if "_VacSkip" in names:
                caught = names
        self.assertIn("OSError", caught)
        self.assertFalse(caught & {"RuntimeError", "Exception", "BaseException", None})

    def test_vac_code_ionrelax_without_locpot(self):
        fn = _exec_func(_func_from(_embedded("_vac_code"), "_vacuum_level"),
                        {"Path": Path, "Locpot": None, "np": None})
        d = self._mk("deform-01", ["LOCPOT", "ionrelax/vasprun.xml"])
        with self.assertRaises(RuntimeError) as cm:
            fn(os.path.join(d, "ionrelax"), None)
        self.assertIn("retry S7_deform 补静态段", str(cm.exception))
        # 嵌入脚本选目录的规则与 _resolved 一致：ionrelax/vasprun.xml 在 -> 只读 ionrelax/
        self.assertIn('dd = d / "ionrelax" if (d / "ionrelax" / "vasprun.xml").is_file() else d',
                      _embedded("_vac_code"))


class S71DeformRefTests(unittest.TestCase):
    def _mod(self):
        os.environ.setdefault("AUTOZT_CONDA_SH", "/dev/null")
        os.environ.setdefault("AUTOZT_AMSET_ENV", "amset-test")
        import importlib.util
        spec = importlib.util.spec_from_file_location("gen_step9b_deform_read_t", str(S71))
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_step_conf_explicit_core(self):
        m = self._mod()
        with tempfile.TemporaryDirectory() as t:
            (Path(t) / "step1_opt").mkdir()
            (Path(t) / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\n")
            self.assertEqual(m._deform_ref_setting(t), "auto")          # 没有 step.conf
            self.assertEqual(m._resolve_deform_ref(t), "vacuum")
            (Path(t) / "step.conf").write_text("[params]\nFUNC = pbe\nDEFORM_REF = core\n")
            self.assertEqual(m._resolve_deform_ref(t), "core")
            (Path(t) / "step.conf").write_text("[params]\nDEFORM_REF = bogus\n")
            with self.assertRaises(SystemExit):
                m._resolve_deform_ref(t)

    def test_vacuum_failure_is_fatal(self):
        s = S71.read_text(encoding="utf-8")
        i = s.index("if rc != 0 or not vac_h5.is_file():")
        blk = s[i:s.index("else:", i)]
        self.assertIn("sys.exit(", blk)
        self.assertNotIn("仍按 core 口径继续", s)
        self.assertIn("DEFORM_REF = core", blk)

    def test_gen_need_ships_step_conf(self):
        y = (SKILL / "skill.yaml").read_text(encoding="utf-8")
        line = next(ln for ln in y.splitlines() if "gen: gen_step9b_deform_read.py" in ln)
        self.assertIn("step.conf", line)
        self.assertIn("stepconf.py", line)

    def test_deform_ref_is_reserved(self):
        for p in (REPO / "stepconf.py", REPO / "skill" / "_common" / "opt" / "stepconf.py"):
            self.assertIn('"DEFORM_REF"', p.read_text(encoding="utf-8"), p)


# ============================== S7 生成端 ==============================
class S7TriggerMatrixTests(unittest.TestCase):
    def test_matrix(self):
        g = _gen9()
        with tempfile.TemporaryDirectory() as t:
            ir = Path(t) / "ionrelax"
            ir.mkdir()
            self.assertTrue(g.static_only_wanted(ir, "2d", "vacuum", False))     # 五条齐 -> 触发
            self.assertFalse(g.static_only_wanted(ir, "2d", "vacuum", True))     # 已有真空口径 h5
            self.assertFalse(g.static_only_wanted(ir, "3d", "vacuum", False))    # 3D
            self.assertFalse(g.static_only_wanted(ir, "2d", "core", False))      # DEFORM_REF=core
            (ir / "LOCPOT").write_text("x\n")
            self.assertFalse(g.static_only_wanted(ir, "2d", "vacuum", False))    # LOCPOT 在

    def test_deform_ref_resolution(self):
        g = _gen9()
        with tempfile.TemporaryDirectory() as t:
            self.assertEqual(g._deform_ref(t, "2d"), "vacuum")
            self.assertEqual(g._deform_ref(t, "3d"), "core")
            (Path(t) / "step.conf").write_text("[params]\nIONRELAX_3D = true\nDEFORM_REF = core\n")
            self.assertEqual(g._deform_ref(t, "2d"), "core")

    def test_vac_h5_checked_before_invalidate(self):
        # invalidate_downstream 会把 step7b_deform_read/deformation_vac.h5 改名归档；判"已有 h5"必须在它之前
        s = (S7 / "gen_step9_deform.py").read_text(encoding="utf-8")
        m = s.index("def main():")
        self.assertLess(s.index("_vac_h5_existed = (cwd / VAC_H5_REL).is_file()", m),
                        s.index("kc.invalidate_downstream(", m))
        self.assertIn("static_only_wanted(_irp, dim, _ref, _vac_h5_existed)", s)
        self.assertEqual(str(g_rel()), "step7b_deform_read/deformation_vac.h5")


def g_rel():
    return _gen9().VAC_H5_REL


class S7StaticOnlyTests(unittest.TestCase):
    def setUp(self):
        self.g = _gen9()
        self.C = _collector()
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)

    def test_prepare_static_only(self):
        g, C = self.g, self.C
        d = _finished_deform(self.td.name, "deform-01", g)
        other = _finished_deform(self.td.name, "deform-02", g, locpot=True)
        ir = d / "ionrelax"
        full = (d / "submit.sh").read_text()
        self.assertEqual(len(_launchers(full)), 3)                 # 刚性 + 弛豫 + 静态
        self.assertTrue(C.ck_deform(str(d), {})[0])                # 改名前：判完成，start 不会交
        stat_before = (ir / "INCAR.static").read_text()

        self.assertTrue(g.prepare_static_only(d))

        # ck_deform 判未完成 -> 进 fan_todo；没动的目录照旧判完成（start 只交这一个）
        ok, why = C.ck_deform(str(d), {})
        self.assertFalse(ok)
        self.assertIn("ionrelax/OUTCAR", why)
        self.assertTrue(C.ck_deform(str(other), {})[0])
        # 静态段产物改名，不删；弛豫段产物、刚性单点结果保留
        names = sorted(p.name for p in ir.iterdir())
        for nm in ("OUTCAR", "vasprun.xml", "OSZICAR"):
            self.assertNotIn(nm, names)
            self.assertTrue(any(re.match(r"%s\.nolocpot-\d{14}$" % re.escape(nm), n) for n in names), names)
        for nm in ("OUTCAR.relax", "CONTCAR.relax", "vasprun.xml.relax", "OSZICAR.relax"):
            self.assertIn(nm, names)
        self.assertEqual((d / "OUTCAR").read_text(), "rigid run\n" + DONE)
        # submit.sh：整体重写，只剩静态段的一条启动命令
        sub = (d / "submit.sh").read_text()
        self.assertEqual(len(_launchers(sub)), 1, sub)
        self.assertTrue(_launchers(sub)[0].startswith("mpirun -np $SLURM_NTASKS /opt/vasp/bin/vasp_std"))
        self.assertLess(sub.index("cd ionrelax"), sub.index("mpirun -np"))
        self.assertIn("cp CONTCAR.relax POSCAR", sub)
        self.assertIn("cp INCAR.static INCAR", sub)
        self.assertNotIn("INCAR.relax", sub)
        self.assertIn("#SBATCH --job-name=MoS2-S7-deform-01", sub)
        self.assertIn("export OMP_NUM_THREADS=1", sub)
        # INCAR.static：ISTART=0 / ICHARG=2，LVHAR / EDIFF / ISYM 不变
        stat = (ir / "INCAR.static").read_text()
        self.assertRegex(stat, r"(?m)^ISTART = 0$")
        self.assertRegex(stat, r"(?m)^ICHARG = 2$")
        self.assertNotRegex(stat, r"(?m)^ISTART\s*=\s*1")
        for key in ("LVHAR", "EDIFF ", "ISYM"):
            self.assertEqual([ln for ln in stat.splitlines() if ln.startswith(key)],
                             [ln for ln in stat_before.splitlines() if ln.startswith(key)], key)
        self.assertIn("LVHAR  = .TRUE.", stat)
        # 静态段还没跑：pending（再 retry 继续只补静态段，不整个重建）
        self.assertTrue(g.static_only_pending(ir))
        self.assertFalse(g.static_only_pending(other / "ionrelax"))

    def test_retry_after_static_rerun_does_not_trigger(self):
        g, C = self.g, self.C
        d = _finished_deform(self.td.name, "deform-01", g)
        ir = d / "ionrelax"
        stat_now = g._ionrelax_incars(INCAR_MAIN)[1]
        self.assertTrue(g.prepare_static_only(d))
        # 模拟静态段跑完：OUTCAR 收尾、vasprun 回显 ISTART=0 / ICHARG=2、LOCPOT 出来
        (ir / "OUTCAR").write_text("static rerun\n" + DONE)
        (ir / "vasprun.xml").write_text(VASPRUN % {"istart": 0, "icharg": 2})
        (ir / "LOCPOT").write_text("new locpot\n")
        self.assertTrue(C.ck_deform(str(d), {})[0])
        self.assertFalse(g.static_only_pending(ir))
        self.assertFalse(g.static_only_wanted(ir, "2d", "vacuum", False))
        # 主循环的 INCAR 一致性判据：按原静态段比不一致，按只补静态段的版本比一致 -> 不整个重建
        kc = g.kc
        self.assertTrue(kc.run_incar_mismatch(stat_now, ir / "vasprun.xml"))
        self.assertEqual(kc.run_incar_mismatch(g.static_only_incar(stat_now), ir / "vasprun.xml"), [])

    def test_static_only_incar_forces_lvhar(self):
        out = self.g.static_only_incar("ISTART = 1  # 接 WAVECAR\nICHARG = 0\nLVHAR = .FALSE.\nEDIFF = 1E-8\n")
        self.assertRegex(out, r"(?m)^LVHAR = \.TRUE\.$")
        self.assertRegex(out, r"(?m)^ISTART = 0$")
        self.assertIn("EDIFF = 1E-8", out)
        self.assertRegex(self.g.static_only_incar("EDIFF = 1E-8\n"), r"(?m)^ICHARG = 2$")   # 没写就补

    def test_old_layout_without_contcar_relax(self):
        g = self.g
        d = _finished_deform(self.td.name, "deform-01", g)
        ir = d / "ionrelax"
        (ir / "CONTCAR.relax").unlink()                             # 旧版 submit.sh 不留 *.relax
        self.assertTrue(g.prepare_static_only(d))
        self.assertEqual((ir / "CONTCAR.relax").read_text(), "relaxed\n")
        (ir / "CONTCAR.relax").unlink()
        (ir / "CONTCAR").unlink()
        self.assertFalse(g.prepare_static_only(d))                  # 弛豫结构都没有：调用方整个重建

    def test_main_loop_wiring(self):
        s = (S7 / "gen_step9_deform.py").read_text(encoding="utf-8")
        m = s[s.index("def main():"):]
        # pending 先于完成判据（否则 OUTCAR 被改名后会走 _build_ionrelax 整个重建）
        self.assertLess(m.index("static_only_pending(_irp)"), m.index("_done = False"))
        # 触发在全部完成判据（网格 / INCAR / 力判据）之后、"跳过重建"之前
        i = m.index("static_only_wanted(_irp, dim, _ref, _vac_h5_existed)")
        self.assertLess(m.index("IONRELAX_REBUILD_LOOSE=true -> 重建"), i)
        self.assertLess(i, m.index("_skipped.append(_name)"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
