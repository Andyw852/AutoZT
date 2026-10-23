"""[V191] ke-dft-cpu S7 ionrelax 静态段的起步磁态 + 磁态 / SCF 核对。

复现 2026-10-23 LS（8 个 Cr）：主单点 INCAR 继承上游的高自旋初猜 MAGMOM = 8*4 8*0 8*0；弛豫段收敛到非磁
（OSZICAR.relax mag=-0.0000）；静态段按初猜从头起算（ISTART=0 / ICHARG=2），ISYM=0 + EDIFF=1E-8 下在磁态
之间来回跳、跑不完被杀。修复后：
  · 新建的 ionrelax：弛豫段写 CHGCAR，静态段 ICHARG=1 从它接（电荷 + 磁化密度），收尾后删 CHGCAR；
  · 只补静态段：弛豫段 CHGCAR 在就 ICHARG=1，否则 ICHARG=2 + MAGMOM 换成弛豫段的收敛磁态；
  · 这些起步方式的变化不触发"旧 INCAR -> 整个重建"（不会删掉已算完的弛豫段）；
  · 已收尾的静态段磁态与弛豫段不符 / 跑满 NELM：只补一次静态段；
  · S7.1 跨构型核对总磁矩（含 S3）与 SCF 收敛，不一致直接报错。
"""
import os
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SKILL = REPO / "skill" / "ke-dft-cpu"
S7 = SKILL / "step7_deform"
S71 = S7 / "step7b_read" / "gen_step9b_deform_read.py"
for _p in (str(SKILL), str(S7)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

DONE = " General timing and accounting informations for this job:\n"


@pytest.fixture(scope="module")
def g9():
    os.environ.setdefault("AUTOZT_CONDA_SH", "/dev/null")
    os.environ.setdefault("AUTOZT_AMSET_ENV", "amset-test")
    import gen_step9_deform
    return gen_step9_deform


@pytest.fixture(scope="module")
def kc():
    import ke_common
    return ke_common


SUBMIT_TPL = """#!/bin/bash
#SBATCH --job-name={{JOBNAME}}
cd $SLURM_SUBMIT_DIR
mpirun -np $SLURM_NTASKS vasp_std
"""

# LS 形态：24 个原子（8 Cr + 8 + 8），主单点 INCAR 里是继承来的高自旋初猜
INCAR_LS = """SYSTEM = LS deform deform-001
ISTART = 0
ICHARG = 2
PREC   = Accurate
ENCUT  = 520
ALGO   = Normal
EDIFF  = 1E-8
NELM   = 200
IBRION = -1
NSW    = 0
ISIF   = 2
ISYM   = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.
LVHAR  = .TRUE.

# ---- 以下由 gen 脚本从 step2_bandgap/step2.1_static 继承（自旋/U 必须与弛豫时一致）----
ISPIN      = 2
MAGMOM     = 8*4  8*0  8*0
"""

POSCAR_LS = ("LS\n1.0\n 10 0 0\n 0 10 0\n 0 0 30\nCr S Se\n8 8 8\nDirect\n"
             + "".join(" 0.%02d 0.%02d 0.5\n" % (i, i) for i in range(24)))

OSZ_RELAX_NM = """       N       E                     dE             d eps       ncg     rms          rms(c)
DAV:   1     0.425437171107E+03    0.42544E+03   -0.37578E+04 10240   0.149E+03
DAV:   2    -0.150912345678E+03   -0.57635E+03   -0.12345E+03 10240   0.456E+02
RMM:   3    -0.150926400000E+03   -0.14054E-01   -0.14000E-02 10240   0.120E-01    0.333E-01
   1 F= -.15092641E+03 E0= -.15092641E+03  d E =-.150926E+03  mag=    -0.0000
DAV:   1    -0.150926410000E+03   -0.10000E-05   -0.10000E-05 10240   0.100E-03
   2 F= -.15092641E+03 E0= -.15092641E+03  d E =-.100000E-05  mag=     0.0000
"""

# 静态段从 32 μB 起步、在磁态之间跳、被墙钟杀掉（没有 F= 行）
OSZ_STATIC_KILLED = """       N       E                     dE             d eps       ncg     rms          rms(c)
DAV:   1     0.425437171107E+03    0.42544E+03   -0.37578E+04 10240   0.149E+03
DAV:   3    -0.151110000000E+03   -0.11000E+02   -0.10000E+01 10240   0.149E+01
DAV:   6    -0.183280000000E+03   -0.26000E+02   -0.10000E+01 10240   0.149E+01
DAV:   7    -0.145630000000E+03    0.38000E+02   -0.10000E+01 10240   0.149E+01
DAV:  11    -0.115220000000E+03    0.30000E+02   -0.10000E+01 10240   0.149E+01
"""


def _osz_done(mag, nelec):
    rows = "".join("DAV: %3d    -0.150926410000E+03   -0.10000E-08   -0.10000E-08 10240   0.100E-05\n" % i
                   for i in range(1, nelec + 1))
    return ("       N       E                     dE             d eps       ncg     rms          rms(c)\n" + rows
            + "   1 F= -.15092641E+03 E0= -.15092641E+03  d E =-.150926E+03  mag=    %.4f\n" % mag)


VASPRUN_ECHO = """<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <incar>
  <i type="string" name="SYSTEM">LS deform deform-001</i>
  <i type="int" name="ISTART">     %(istart)s</i>
  <i type="int" name="ICHARG">     %(icharg)s</i>
  <i type="string" name="PREC">accurate</i>
  <i name="ENCUT">    520.00000000</i>
  <i type="string" name="ALGO">Normal</i>
  <i name="EDIFF">      %(ediff)s</i>
  %(ediffg)s
  <i type="int" name="NELM">    200</i>
  <i type="int" name="IBRION">    %(ibrion)s</i>
  <i type="int" name="NSW">     %(nsw)s</i>
  <i type="int" name="ISIF">     2</i>
  <i type="int" name="ISYM">     %(isym)s</i>
  <i type="logical" name="LWAVE"> %(lwave)s  </i>
  <i type="logical" name="LCHARG"> F  </i>
  <i type="logical" name="LVHAR"> %(lvhar)s  </i>
  <i type="int" name="ISPIN">     2</i>
  <v name="MAGMOM">      4.00000000      4.00000000      4.00000000      4.00000000      4.00000000      4.00000000      4.00000000      4.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000</v>
 </incar>
 <parameters/>
</modeling>
"""
RELAX_ECHO = VASPRUN_ECHO % {"istart": 0, "icharg": 2, "ediff": "0.00000100",
                             "ediffg": '<i name="EDIFFG">     -0.00500000</i>', "ibrion": 2, "nsw": 60,
                             "isym": 2, "lwave": "T", "lvhar": "F"}
STATIC_ECHO_OLD = VASPRUN_ECHO % {"istart": 1, "icharg": 0, "ediff": "0.00000001", "ediffg": "",
                                  "ibrion": -1, "nsw": 0, "isym": 0, "lwave": "F", "lvhar": "T"}


def _ls_dir(root, g9):
    """LS deform-001：用（修复后的）_build_ionrelax 建输入，再摆上旧版跑出来的产物：
    弛豫段已收敛到非磁（旧 INCAR.relax，LCHARG=F，没有 CHGCAR），静态段被杀。"""
    d = Path(root) / "step7_deform" / "deform-001"
    d.mkdir(parents=True)
    (d / "POSCAR").write_text(POSCAR_LS)
    (d / "KPOINTS").write_text("auto\n0\nGamma\n 6 42 3\n")
    (d / "POTCAR").write_text("potcar\n")
    (d / "INCAR").write_text(INCAR_LS)
    (d / "submit.sh").write_text(SUBMIT_TPL.replace("{{JOBNAME}}", "LS-S7-deform-001"))
    g9._build_ionrelax(d, 520, {}, SUBMIT_TPL)
    (d / "OUTCAR").write_text("rigid\n" + DONE)
    ir = d / "ionrelax"
    (ir / "OUTCAR.relax").write_text("relax\n" + DONE)
    (ir / "OSZICAR.relax").write_text(OSZ_RELAX_NM)
    (ir / "CONTCAR.relax").write_text(POSCAR_LS)
    (ir / "vasprun.xml.relax").write_text(RELAX_ECHO)
    (ir / "OUTCAR").write_text("static running ...\n")          # 被墙钟杀掉：没有收尾横幅
    (ir / "OSZICAR").write_text(OSZ_STATIC_KILLED)
    return d


# ---------------------------------------------------------------- ke_common 解析
def test_oszicar_state(kc, tmp_path):
    p = tmp_path / "OSZICAR"
    p.write_text(OSZ_RELAX_NM)
    st = kc.oszicar_state(p)
    assert st["n_ionic"] == 2 and st["nelec"] == 1 and abs(st["mag"]) < 1e-9
    assert abs(st["F"] + 150.92641) < 1e-6
    p.write_text(OSZ_STATIC_KILLED)
    assert kc.oszicar_state(p) is None                              # 第一步都没跑完
    p.write_text(_osz_done(22.98, 200))
    (tmp_path / "INCAR").write_text("NELM = 200\n")
    st = kc.scf_state(tmp_path)
    assert st["nelm"] == 200 and not st["scf_ok"] and abs(st["mag"] - 22.98) < 1e-9
    p.write_text(_osz_done(0.0, 31))
    assert kc.scf_state(tmp_path)["scf_ok"]
    (tmp_path / "INCAR").write_text("EDIFF = 1E-8\n")               # 没写 NELM：VASP 默认 60
    p.write_text(_osz_done(0.0, 60))
    assert not kc.scf_state(tmp_path)["scf_ok"]
    # ISPIN=1：没有 mag=
    p.write_text("DAV:   1    -0.1E+03   -0.1E-08   -0.1E-08 100   0.1E-05\n   1 F= -.1E+03 E0= -.1E+03  d E =-.1E+03\n")
    assert kc.oszicar_state(p)["mag"] is None
    # 非共线：mag= 三个分量，取模
    p.write_text("DAV:   1    -0.1E+03   -0.1E-08   -0.1E-08 100   0.1E-05\n"
                 "   1 F= -.1E+03 E0= -.1E+03  d E =-.1E+03  mag=     3.0000     0.0000     4.0000\n")
    assert abs(kc.oszicar_state(p)["mag"] - 5.0) < 1e-9


def test_same_mag_state(kc):
    assert not kc.same_mag_state(0.0, 22.98)                        # 非磁 vs 铁磁
    assert not kc.same_mag_state(0.0, 2.0)
    assert kc.same_mag_state(24.0, 23.2)                            # 1% 应变下铁磁磁矩连续变一点
    assert kc.same_mag_state(0.0, 0.3)
    assert kc.same_mag_state(None, 24.0)                            # ISPIN=1 / 读不到：判不了，放行


def test_relaxed_magmom(kc, tmp_path):
    (tmp_path / "OSZICAR.relax").write_text(OSZ_RELAX_NM)
    new, why = kc.relaxed_magmom(tmp_path, "8*4  8*0  8*0", 24)
    assert new == "24*0" and "非磁" in why
    (tmp_path / "OSZICAR.relax").write_text(_osz_done(24.0, 20))
    new, _ = kc.relaxed_magmom(tmp_path, "8*4  8*0  8*0", 24)
    assert new == "8*3 16*0"                                        # 等比缩放到弛豫段的总磁矩
    new, why = kc.relaxed_magmom(tmp_path, "4*4 4*-4 16*0", 24)
    assert new is None and "反铁磁" in why                          # 总磁矩分不出 AFM 与非磁：不动
    assert kc.relaxed_magmom(tmp_path, "8*4  8*0  8*0", 23)[0] is None   # 个数对不上
    assert kc.relaxed_magmom(tmp_path, " ".join(["0 0 4"] * 24), 24)[0] is None   # 非共线
    # OUTCAR.relax 有逐离子磁矩（LORBIT>=10）就直接用
    table = ("\n magnetization (x)\n\n# of ion       s       p       d       tot\n"
             "------------------------------------------\n"
             + "".join("    %d        0.010   0.020   %.3f   %.3f\n" % (i + 1, m, m)
                       for i, m in enumerate([2.5, -2.5] * 4 + [0.0] * 16))
             + "--------------------------------------------------\ntot    0.0 0.0 0.0 0.0\n")
    (tmp_path / "OUTCAR.relax").write_text(table)
    new, why = kc.relaxed_magmom(tmp_path, "8*4  8*0  8*0", 24)
    assert new.startswith("2.5 -2.5 2.5") and new.endswith("16*0") and "逐离子" in why


def test_magmom_roundtrip(kc):
    assert kc.expand_magmom("8*4  8*0  8*0") == [4.0] * 8 + [0.0] * 16
    assert kc.compress_magmom([-0.0, 0.0, 1.25, 1.25, -3.0]) == "2*0 2*1.25 -3"
    assert kc.expand_magmom("a b") is None


def test_scf_state_issues(kc):
    nm = {"mag": 0.0, "nelec": 30, "nelm": 200, "scf_ok": True}
    fm = {"mag": 22.98, "nelec": 40, "nelm": 200, "scf_ok": True}
    nc = {"mag": 0.0, "nelec": 200, "nelm": 200, "scf_ok": False}
    assert kc.scf_state_issues({"undeformed": nm, "deform-001/ionrelax": nm}, "undeformed") == []
    iss = kc.scf_state_issues({"undeformed": nm, "deform-001/ionrelax": fm, "deform-002/ionrelax": nc},
                              "undeformed")
    assert len(iss) == 2
    assert any("deform-001/ionrelax" in x and "磁态" in x for x in iss)
    assert any("deform-002/ionrelax" in x and "NELM=200" in x for x in iss)
    iss = kc.scf_state_issues({"undeformed": nm}, "undeformed", {"step3_uniform": fm})
    assert len(iss) == 1 and "step3_uniform" in iss[0]
    assert kc.scf_state_issues({"undeformed": None, "deform-001": fm}, "undeformed") == []   # 参考判不了


# ---------------------------------------------------------------- S7 gen：新建的 ionrelax
def test_full_path_static_restarts_from_relax_chgcar(g9, tmp_path):
    d = _ls_dir(tmp_path, g9)
    ir = d / "ionrelax"
    relax = (ir / "INCAR.relax").read_text()
    stat = (ir / "INCAR.static").read_text()
    assert re.search(r"(?m)^LCHARG = \.TRUE\.", relax)
    assert re.search(r"(?m)^ISTART = 0$", stat) and re.search(r"(?m)^ICHARG = 1\b", stat)
    assert re.search(r"(?m)^LCHARG\s*=\s*\.FALSE\.", stat)          # 静态段不写 CHGCAR，不会盖掉弛豫段那份
    assert "ISTART = 1" not in stat
    sub = (d / "submit.sh").read_text()
    i_static = sub.index("cp INCAR.static INCAR")
    assert sub.index("[ -s CHGCAR ]") < i_static                     # 静态段之前核对 CHGCAR 在
    assert sub.index("rm -f CHGCAR") > i_static                      # 收尾后才删
    assert 'grep -q "General timing and accounting informations" OUTCAR 2>/dev/null && rm -f CHGCAR' in sub


# ---------------------------------------------------------------- S7 gen：LS 现场（弛豫段旧、静态段被杀）
def test_ls_static_killed_keeps_relax_and_starts_nonmagnetic(g9, tmp_path, capsys):
    d = _ls_dir(tmp_path, g9)
    ir = d / "ionrelax"
    assert g9.relax_done_static_dead(ir)
    # 主循环的弛豫段输入核对：旧弛豫段 LCHARG=F，现在 LCHARG=.TRUE. —— 只是输出开关，不能因此整个重建
    relax_now, stat_now = g9._ionrelax_incars(INCAR_LS)
    import ke_common as kc
    assert any(w.startswith("LCHARG:") for w in kc.run_incar_mismatch(relax_now, ir / "vasprun.xml.relax"))
    assert g9._mismatch_wo_lwave(relax_now, ir / "vasprun.xml.relax") == []
    # 物理输入真变了照样认
    assert g9._mismatch_wo_lwave(relax_now.replace("ENCUT  = 520", "ENCUT  = 600"), ir / "vasprun.xml.relax")
    (ir / "INCAR.static").write_text(stat_now)
    assert g9.prepare_static_only(d, False, "静态段没跑完（被杀/OOM），弛豫段已完成", tag=g9.STATIC_DEAD_TAG)
    out = capsys.readouterr().out
    stat = (ir / "INCAR.static").read_text()
    assert re.search(r"(?m)^ICHARG = 2$", stat)                      # 没有弛豫段 CHGCAR：从头起算
    assert re.search(r"(?m)^MAGMOM = 24\*0\b", stat), stat          # ……但按弛豫段的非磁态起步，不再是 8*4
    assert "非磁" in out
    names = {p.name for p in ir.iterdir()}
    assert {"OUTCAR.relax", "OSZICAR.relax", "CONTCAR.relax", "vasprun.xml.relax"} <= names
    sub = (d / "submit.sh").read_text()
    assert "[ -s CHGCAR ]" not in sub and sub.count("vasp_std") == 1
    assert g9.static_only_pending(ir)
    # 再 retry 一次（静态段还没跑）：MAGMOM 不会被反复改写
    assert g9.prepare_static_only(d, False, "静态段补算未完成")
    assert re.search(r"(?m)^MAGMOM = 24\*0\b", (ir / "INCAR.static").read_text())


def test_static_only_uses_relax_chgcar_when_present(g9, tmp_path):
    d = _ls_dir(tmp_path, g9)
    ir = d / "ionrelax"
    (ir / "CHGCAR").write_text("relax chgcar\n")                     # 新版弛豫段留下的
    assert g9.prepare_static_only(d, False, "静态段没跑完", tag=g9.STATIC_DEAD_TAG)
    stat = (ir / "INCAR.static").read_text()
    assert re.search(r"(?m)^ICHARG = 1$", stat)
    assert "24*0" not in stat                                        # ICHARG=1 时 MAGMOM 不用改（VASP 忽略它）
    sub = (d / "submit.sh").read_text()
    assert sub.index("[ -s CHGCAR ]") < sub.index("vasp_std")
    assert sub.rstrip().splitlines()[-1].startswith("rm -f CHGCAR")
    # use_chg=False（坏掉的静态段可能自己写过 CHGCAR）：一律按弛豫段磁态从头起算
    assert g9.prepare_static_only(d, False, "静态段 SCF 没收敛", tag=g9.BADSTATIC_TAG, use_chg=False)
    stat = (ir / "INCAR.static").read_text()
    assert re.search(r"(?m)^ICHARG = 2$", stat) and re.search(r"(?m)^MAGMOM = 24\*0\b", stat)
    assert g9.static_only_pending(ir)                                # badstatic 归档也认作"补算中"


def test_static_reuse_ignores_restart_tags(g9, tmp_path):
    """已算完的静态段（旧版 ISTART=1/ICHARG=0、MAGMOM 初猜）：起步方式变了不算输入变了，不整个重建。"""
    p = tmp_path / "vasprun.xml"
    p.write_text(STATIC_ECHO_OLD)
    stat_now = g9._ionrelax_incars(INCAR_LS)[1]
    assert g9._mismatch_wo_lwave(stat_now, p, static=True) == []
    stat_so = g9.static_only_incar(stat_now)                         # 只补过静态段的回显（ICHARG=2）
    assert g9._mismatch_wo_lwave(stat_so, p, static=True) == []
    stat_mag = g9._set_incar_tag(stat_now, "MAGMOM", "24*0")
    assert g9._mismatch_wo_lwave(stat_mag, p, static=True) == []
    assert any(w.startswith("MAGMOM:") for w in g9._mismatch_wo_lwave(stat_mag, p))   # 弛豫段照比 MAGMOM
    assert g9._mismatch_wo_lwave(g9._set_incar_tag(stat_now, "ICHARG", "11"), p, static=True)   # 非自洽是物理
    assert g9._mismatch_wo_lwave(stat_now.replace("EDIFF  = 1E-8", "EDIFF  = 1E-6"), p, static=True)


def test_static_state_issue(g9, tmp_path):
    ir = tmp_path / "ionrelax"
    ir.mkdir()
    (ir / "INCAR").write_text("NELM = 200\nISPIN = 2\n")
    (ir / "OSZICAR.relax").write_text(OSZ_RELAX_NM)
    (ir / "OSZICAR").write_text(_osz_done(0.0, 35))
    assert g9.static_state_issue(ir) is None
    (ir / "OSZICAR").write_text(_osz_done(22.98, 60))
    assert "磁态与弛豫段不符" in g9.static_state_issue(ir)
    (ir / "OSZICAR").write_text(_osz_done(0.0, 200))
    assert "NELM=200" in g9.static_state_issue(ir)
    (ir / "OSZICAR").write_text(OSZ_STATIC_KILLED)
    assert g9.static_state_issue(ir) is None                         # 没收尾的不归它管（[0067] 那条路）


def test_main_loop_wiring(g9):
    s = (S7 / "gen_step9_deform.py").read_text(encoding="utf-8")
    m = s[s.index("def main():"):]
    i = m.index("static_state_issue(_irp)")
    assert m.index("IONRELAX_REBUILD_LOOSE=true -> 重建") < i           # 要整个重建的先走重建
    assert i < m.index("static_only_wanted(_irp, dim, _ref, _vac_h5_existed)")
    assert i < m.index("_skipped.append(_name)")
    assert "tag=BADSTATIC_TAG, use_chg=False" in m
    assert '_mismatch_wo_lwave(_stat_now, _irp / "vasprun.xml", static=True)' in m


def test_s71_guard_wiring():
    s = S71.read_text(encoding="utf-8")
    m = s[s.index("def main():"):]
    i = m.index('kc.scf_state_issues(_scf, "undeformed", _scf_up)')
    assert m.index("kc.outcar_core_ref_kind") < i < m.index("amset deform read 在形变目录里跑")
    assert '"step3_uniform", "step3b_uniform_full"' in m
    assert '_be["scf_state_check"] = _scf_check' in s
