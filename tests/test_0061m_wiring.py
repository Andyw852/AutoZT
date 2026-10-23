"""[0061-m] DFT 形变势矩阵的接线（纯本地、不需要 amset / VASP）。

  · S7（gen_step9_deform）：LWAVE 只开在"本征值取自它"的那次计算（ionrelax 静态段 / 没有 ionrelax 的顶层单点），
    不改模板；LWAVE 不算"输入变了"；已完成但缺 deform_wave.h5 的 ionrelax 只补静态段（归档后缀 nowave）；
    submit.sh 末尾的抽取段（VASP 正常收尾才抽、失败不让作业判失败、相对根随目录层级）；step.conf 开关。
  · S7.1（gen_step9b_deform_read）：deform_wave build 真跑一遍（EPM 构型）；缺文件 / 没开 DEFORM_WAVE 只告警。
  · S8 / S8.4 gen：有 deformation_matrix.h5 才记 matrix_h5、软链、带上 deform_wave.py；ADP_DEGEN_MATRIX = off。
  · skill.yaml：S7 / S7.1 / S8 / S8.4 的 gen_need 都有 deform_wave.py；dp_symmetrize.write_signed 存 ibands。
"""
import io
import os
import re
import sys
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
S7 = KE / "step7_deform"
S71 = S7 / "step7b_read"
S84 = KE / "step8.4_amset2d"
for _p in (str(ROOT), str(KE), str(S7), str(S71), str(Path(__file__).resolve().parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

pytest.importorskip("h5py")


def _quiet(fn, *a, **kw):
    buf = io.StringIO()
    with redirect_stdout(buf):
        r = fn(*a, **kw)
    return r, buf.getvalue()


@pytest.fixture(scope="module")
def g9():
    os.environ.setdefault("AUTOZT_CONDA_SH", "/dev/null")
    os.environ.setdefault("AUTOZT_AMSET_ENV", "amset-test")
    import gen_step9_deform
    return gen_step9_deform


@pytest.fixture(scope="module")
def g9b():
    os.environ.setdefault("AUTOZT_CONDA_SH", "/dev/null")
    os.environ.setdefault("AUTOZT_AMSET_ENV", "amset-test")
    import gen_step9b_deform_read
    return gen_step9b_deform_read


INCAR = """SYSTEM = Si deform deform-001
ISTART = 0
ICHARG = 2
ENCUT  = 368
IBRION = -1
NSW    = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.
LVHAR  = .TRUE.
ISYM   = 0
"""


# ---------------------------------------------------------------------------
# S7
# ---------------------------------------------------------------------------
def test_s7_templates_untouched():
    for dim in ("2d", "3d"):
        assert re.search(r"(?m)^LWAVE\s*=\s*\.FALSE\.", (S7 / ("incar_deform_%s.tpl" % dim)).read_text(encoding="utf-8"))


def test_s7_wave_dirs(g9):
    names = ["undeformed", "deform-001", "deform-002", "deform-003"]
    assert g9.wave_top_dirs(names, {"deform-001", "deform-002"}, "3d") == {"undeformed", "deform-003"}
    assert g9.wave_top_dirs(names, set(names[1:]), "3d") == {"undeformed"}
    assert g9.wave_top_dirs(names, {"deform-001"}, "2d") == {"undeformed"}      # 2D 面外的离子固定构型不要


def test_s7_static_segment_lwave(g9):
    relax, stat = g9._ionrelax_incars(INCAR, wave=True)
    assert re.search(r"(?m)^LWAVE = \.TRUE\.$", stat)
    _r0, stat0 = g9._ionrelax_incars(INCAR)
    assert re.search(r"(?m)^LWAVE\s*=\s*\.FALSE\.", stat0)                       # 关着时与原来逐字相同
    assert stat0 == g9._ionrelax_incars(INCAR, wave=False)[1]
    so = g9.static_only_incar(stat0, wave=True)
    assert re.search(r"(?m)^LWAVE = \.TRUE\.$", so) and re.search(r"(?m)^ISTART = 0$", so)
    assert g9.static_only_incar(stat0) == g9.static_only_incar(stat0, wave=False)


VASPRUN = """<?xml version="1.0"?>
<modeling>
 <incar>
  <i type="string" name="SYSTEM">Si</i>
  <i name="ENCUT">    368.0</i>
  <i type="int" name="IBRION">   -1</i>
  <i type="int" name="NSW">    0</i>
  <i type="int" name="ISTART">    1</i>
  <i type="int" name="ICHARG">    0</i>
  <i type="logical" name="LWAVE"> F  </i>
  <i type="logical" name="LCHARG"> F  </i>
  <i type="logical" name="LVHAR"> T  </i>
  <i type="int" name="ISYM">    0</i>
 </incar>
 <parameters></parameters>
</modeling>
"""


def test_s7_lwave_not_an_input_change(g9, tmp_path):
    (tmp_path / "vasprun.xml").write_text(VASPRUN)
    stat = g9._ionrelax_incars(INCAR, wave=True)[1]
    import ke_common as kc
    assert any(w.startswith("LWAVE:") for w in kc.run_incar_mismatch(stat, tmp_path / "vasprun.xml"))
    assert g9._mismatch_wo_lwave(stat, tmp_path / "vasprun.xml") == []
    assert g9._mismatch_wo_lwave(stat.replace("368", "400"), tmp_path / "vasprun.xml")      # 物理输入变了照样认


SUBMIT = """#!/bin/bash
#SBATCH --job-name=x
cd $SLURM_SUBMIT_DIR
mpirun -np $SLURM_NTASKS vasp_std
"""


def test_s7_prepare_static_only_for_wave(g9, tmp_path):
    d = tmp_path / "deform-001"
    ir = d / "ionrelax"
    ir.mkdir(parents=True)
    (d / "INCAR").write_text(INCAR)
    (d / "submit.sh").write_text(SUBMIT)
    (ir / "INCAR.static").write_text(g9._ionrelax_incars(INCAR)[1])
    (ir / "CONTCAR.relax").write_text("relaxed\n")
    (ir / "OUTCAR").write_text(" General timing and accounting informations for this job:\n")
    (ir / "vasprun.xml").write_text(VASPRUN)
    (ir / "LOCPOT").write_text("x")
    ok, log = _quiet(g9.prepare_static_only, d, True, "缺 deform_wave.h5")
    assert ok and "缺 deform_wave.h5" in log
    assert list(ir.glob("OUTCAR.nowave-*")) and not list(ir.glob("OUTCAR.nolocpot-*"))
    st = (ir / "INCAR.static").read_text()
    assert re.search(r"(?m)^LWAVE = \.TRUE\.$", st) and re.search(r"(?m)^ICHARG = 2$", st)
    assert g9.static_only_pending(ir)                                   # nowave 归档也认作"补算中"
    sub = (d / "submit.sh").read_text()
    assert "cd ionrelax" in sub and sub.count("vasp_std") == 1


def test_s7_extract_block(g9, tmp_path):
    blk = g9.wave_extract_block("../../../", 6, False)
    assert "../../../deform_wave.py extract --margin 6 --delete" in blk
    assert 'grep -q "General timing and accounting informations" OUTCAR' in blk
    assert "|| echo" in blk                                             # 抽取失败不让作业判失败
    assert "--delete" not in g9.wave_extract_block("../../", 4, True)
    # bash 语法与"VASP 没收尾就不抽"：用假的 python 跑一遍
    d = tmp_path / "run"
    d.mkdir()
    fake = tmp_path / "deform_wave.py"
    fake.write_text("import sys; open('called', 'w').write(' '.join(sys.argv[1:]))\n")
    blk = g9.wave_extract_block("../", 6, False).replace(g9.AMSET_ENV_SRC, "CONDA_PREFIX=%s" % Path(sys.executable).parent.parent)
    blk = blk.replace('"${CONDA_PREFIX:?conda activate 没生效}/bin/python"', sys.executable)
    (d / "WAVECAR").write_text("w")
    import subprocess
    assert subprocess.run(["bash", "-c", blk], cwd=d).returncode == 0
    assert not (d / "called").exists()                                  # OUTCAR 没收尾 -> 不抽
    (d / "OUTCAR").write_text(" General timing and accounting informations for this job:\n")
    assert subprocess.run(["bash", "-c", blk], cwd=d).returncode == 0
    assert (d / "called").read_text() == "extract --margin 6 --delete"


def test_s7_step_conf_params(g9, tmp_path):
    assert g9._deform_wave_params(tmp_path) == (True, False, 6)
    (tmp_path / "step.conf").write_text("STEP = step7_deform\nDEFORM_WAVE = false\nDEFORM_WAVE_KEEP = true\n"
                                        "DEFORM_WAVE_MARGIN = 8\n")
    import stepconf
    os.chdir(tmp_path)
    try:
        assert g9._deform_wave_params(tmp_path) == (False, True, 8)
    finally:
        os.chdir(ROOT)
    _ = stepconf


def test_s7_main_wiring(g9):
    s = (S7 / "gen_step9_deform.py").read_text(encoding="utf-8")
    m = s[s.index("def main():"):]
    # matrix h5 的存在性在 invalidate_downstream 之前取（它会把下游产物归档）
    assert m.index("_mat_h5_existed = ") < m.index("kc.invalidate_downstream(")
    # ionrelax 目录在主循环之前选好（顶层要不要写 WAVECAR 取决于它）
    assert m.index("_ir_dirs = select_ionrelax_dirs(") < m.index("for d in subs:")
    # 顶层 LWAVE 在 step.conf 覆盖之后、V168 核对之前
    i = m.index("if d.name in _wave_top:")
    assert m.index("stepconf.apply_incar_file(d / \"INCAR\"") < i < m.index("kc.archive_if_run_mismatch(d)")
    # 缺波函数的只补静态段在全部完成判据之后、"跳过重建"之前；抽取段最后统一追加
    j = m.index("_need_wv = ")
    assert m.index("IONRELAX_REBUILD_LOOSE=true -> 重建") < j < m.index("_skipped.append(_name)")
    assert m.index("_skipped.append(_name)") < m.index("wave_extract_block(_rel")
    assert "_mismatch_wo_lwave(_stat_now" in m


# ---------------------------------------------------------------------------
# S7.1
# ---------------------------------------------------------------------------
def test_s71_build_skips_without_waves(g9b, tmp_path, monkeypatch):
    dfm = tmp_path / "step7_deform"
    (dfm / "undeformed").mkdir(parents=True)
    out = tmp_path / "step7b_deform_read"
    out.mkdir()
    monkeypatch.chdir(tmp_path)
    st, log = _quiet(g9b.build_deform_matrix, dfm, out, ["deform-001"])
    assert st.startswith("未生成（缺 deform_wave.py") and "gen_need" in log
    (tmp_path / "deform_wave.py").write_text("")
    st, log = _quiet(g9b.build_deform_matrix, dfm, out, ["deform-001"])
    assert st == "未生成（S7 没有波函数）" and not (out / g9b.MATRIX_H5).exists()


def test_s71_build_runs_deform_wave(g9b, tmp_path, monkeypatch):
    import test_0061m_deform_wave as TW
    import shutil
    root, dirs, _d0 = TW._s7(tmp_path / "w", n=2)
    for d in [root / "undeformed"] + [root / x for x in dirs]:
        (d / "WAVECAR").unlink()
    (root / dirs[0] / "deform_wave.h5").unlink()                       # 一个构型没抽到：跳过、照样 build
    proj = tmp_path / "w"
    shutil.copy(KE / "deform_wave.py", proj / "deform_wave.py")
    out = proj / "step7b_deform_read"
    out.mkdir()
    (out / g9b.MATRIX_H5).write_text("stale")
    monkeypatch.chdir(proj)
    monkeypatch.setattr(g9b, "AMSET_ENV_SRC", "true")
    monkeypatch.setattr(g9b, "_ENV_ECHO", "true")
    monkeypatch.setattr(g9b, "AMSET_PY", sys.executable)
    st, log = _quiet(g9b.build_deform_matrix, root, out, dirs)
    assert st == "%s（%d 个构型）" % (g9b.MATRIX_H5, len(dirs) - 1), log
    import deform_wave as DW
    m = DW.load_matrix(out / g9b.MATRIX_H5, "up")
    assert len(m["pairs"]) == 6 and m["X"].shape[0] == 6
    # build 失败（例如 dirs 里全是坏路径）：只告警、删掉半成品、不抛
    st, log = _quiet(g9b.build_deform_matrix, root, out, ["nonexistent"])
    assert st.startswith("未生成") or st.startswith("build 失败")


def test_s71_wiring(g9b):
    s = (S71 / "gen_step9b_deform_read.py").read_text(encoding="utf-8")
    m = s[s.index("def main():"):]
    assert "out / MATRIX_H5):" in m                                    # 旧矩阵先删
    assert m.index("build_deform_matrix(dfm, out, _rel)") < m.index("os.replace(str(h5), str(dst))")


def test_write_signed_stores_ibands(tmp_path):
    import h5py
    import dp_symmetrize as DPS

    class _S:
        def __init__(self, n):
            self.name = n
    up = _S("up")
    DPS.write_signed(tmp_path / "s.h5", {up: np.zeros((5, 2, 3, 3))}, np.zeros((2, 3)), ibands={up: np.arange(1, 4)})
    with h5py.File(tmp_path / "s.h5", "r") as f:
        assert list(f["ibands_up"][()]) == [1, 2, 3] and f["signed_up"].shape == (3, 2, 3, 3)


# ---------------------------------------------------------------------------
# S8 / S8.4 gen
# ---------------------------------------------------------------------------
from test_0061_adp_degen import _load, _resolve  # noqa: E402


def test_gen_s8_matrix(tmp_path):
    g = _load("g10_matrix", KE / "step8_amset" / "gen_step10_amset.py")
    g.SCATTERING = ["ADP", "IMP"]
    g.ADP_DEGEN = "auto"
    cwd = tmp_path / "ke"
    rd = cwd / g.READ_DIR
    rd.mkdir(parents=True)
    (rd / "deformation_signed.h5").write_bytes(b"x")
    rec, log = _resolve(g.resolve_adp_degen_3d, False, cwd)
    assert rec["status"] == "pending_runtime" and rec["matrix"] == "auto" and "matrix_h5" not in rec
    assert "从本征值拟合" in log
    (rd / "deformation_matrix.h5").write_bytes(b"x")
    rec, log = _resolve(g.resolve_adp_degen_3d, False, cwd)
    assert rec["matrix_h5"] == "deformation_matrix.h5" and "DFT 形变势矩阵" in log
    g.ADP_DEGEN_MATRIX = "off"
    rec, _log = _resolve(g.resolve_adp_degen_3d, False, cwd)
    assert rec["matrix"] == "off" and "matrix_h5" not in rec
    g.ADP_DEGEN_MATRIX = "bogus"
    with pytest.raises(SystemExit, match="ADP_DEGEN_MATRIX"):
        _resolve(g.resolve_adp_degen_3d, False, cwd)
    g.ADP_DEGEN_MATRIX = "auto"
    out = tmp_path / "out"
    out.mkdir()
    r = {"mechanisms": [], "dp_optical": None,
         "adp_degen": {"status": "pending_runtime", "matrix_h5": "deformation_matrix.h5"}}
    with redirect_stdout(io.StringIO()):
        assert g._install_dp_optical_3d(out, r) is True
    assert (out / "deform_wave.py").is_file() and (out / "adp_degen.py").is_file()
    out2 = tmp_path / "out2"
    out2.mkdir()
    with redirect_stdout(io.StringIO()):
        g._install_dp_optical_3d(out2, dict(r, adp_degen={"status": "pending_runtime"}))
    assert not (out2 / "deform_wave.py").exists()
    assert "ADP_DEGEN_MATRIX" in g.SPEC
    src = (KE / "step8_amset" / "gen_step10_amset.py").read_text(encoding="utf-8")
    assert "link(out, cwd / READ_DIR / DEFORM_MATRIX_H5, DEFORM_MATRIX_H5)" in src
    assert g._SKILL_REV[:10] >= "2026-10-21"          # [0072] 之后的版本号只会更新


def test_gen_s84_matrix(tmp_path, monkeypatch):
    for d in (ROOT, KE, S84):
        monkeypatch.syspath_prepend(str(d))
    g = _load("g14_matrix", S84 / "gen_step14_amset2d.py")
    monkeypatch.setattr(g, "SCATTERING", ["ADP", "IMP", "POP"])
    cwd = tmp_path / "ke2d"
    rd = cwd / g.READ_DIR
    rd.mkdir(parents=True)
    for f in ("deformation.h5", "deformation_vac.h5", "deformation_vac_signed.h5", "deformation_matrix.h5"):
        (rd / f).write_bytes(b"x")
    rec, _log = _resolve(g.resolve_adp_degen_2d, cwd)
    assert rec["status"] == "pending_runtime" and rec["matrix_h5"] == "deformation_matrix.h5"
    out = tmp_path / "out84"
    out.mkdir()
    monkeypatch.setattr(g, "_ADP_DEGEN_REC", rec)
    assert g._install_adp_degen(out) is True
    assert (out / "adp_degen.py").is_file() and (out / "deform_wave.py").is_file()
    src = (S84 / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert "link(out, cwd / READ_DIR / DEFORM_MATRIX_H5, DEFORM_MATRIX_H5)" in src
    assert "ADP_DEGEN_MATRIX" in g.SPEC and g._SKILL_REV.startswith("2026-10-21-0061m")


def test_skill_yaml_gen_need_deform_wave():
    import yaml
    d = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    found = {}

    def walk(x):
        if isinstance(x, dict):
            if x.get("gen") in ("gen_step9_deform.py", "gen_step9b_deform_read.py", "gen_step10_amset.py",
                                "gen_step14_amset2d.py"):
                found[x["gen"]] = "deform_wave.py" in (x.get("gen_need") or [])
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(d)
    assert found == {"gen_step9_deform.py": True, "gen_step9b_deform_read.py": True, "gen_step10_amset.py": True,
                     "gen_step14_amset2d.py": True}
