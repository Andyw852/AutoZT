# -*- coding: utf-8 -*-
"""[hse_two_stage] ke-dft-cpu S2.3 HSE 两段式（纯本地、不需要 VASP）。

  · 第 1 段 S2.28_hsescf（gen_step4_HSE.py --stage scf）：只有一个目录、没有 KPOINTS_OPT、带 NKRED（按网格能否整除）、
    WAVECAR 从 S2.2 PBE 拷；
  · 第 2 段 S2.3_hse（默认 --stage bands）：第 1 段就绪 -> 各片 WAVECAR / KPOINTS / NBANDS 取自第 1 段、不带 NKRED、
    写 hse_scf_source.json；第 1 段目录在但没跑完 -> 报错；第 1 段目录不存在（老材料）-> 退回旧口径并告警；
  · 判据 hse_scf：第 1 段收尾 + WAVECAR 才算完成；老材料 S2.3 已按单段算完 -> 判完成（不再交作业）；
  · skill.yaml：S2.3 依赖 S2.28；ke_common.DOWNSTREAM：S2.28 重算 -> S2.3 整目录归档。
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KE = REPO / "skill" / "ke-dft-cpu"
GEN = KE / "step2_bandgap" / "step2.3_hse" / "gen_step4_HSE.py"
DONE = " General timing and accounting informations for this job:\n"
DISC = KE / "step2_bandgap" / "step2.15_discriminant"

INCAR_PBE = """SYSTEM = MgC PBE
ISTART = 0
ICHARG = 2
ENCUT = 400
ISMEAR = 0
SIGMA = 0.05
GGA = PE
IVDW = 12
ISYM = 2
NCORE = 4
LWAVE = .TRUE.
"""
POSCAR = """MgC
1.0
 5.0 0.0 0.0
 0.0 5.0 0.0
 0.0 0.0 6.0
Mg C
1 2
Direct
 0.0 0.0 0.0
 0.5 0.5 0.3
 0.5 0.5 0.7
"""
SUBMIT = """#!/bin/bash
#SBATCH --job-name={{JOBNAME}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
mpirun -np $SLURM_NTASKS vasp_std
"""


def _material(td, kmesh="4 4 3"):
    m = Path(td)
    pbe = m / "step2_bandgap" / "step2.2_pbe"
    pbe.mkdir(parents=True)
    (pbe / "INCAR").write_text(INCAR_PBE)
    (pbe / "POSCAR").write_text(POSCAR)
    (pbe / "POTCAR").write_text("potcar\n")
    (pbe / "KPOINTS").write_text("auto\n0\nGamma\n %s\n" % kmesh)
    (pbe / "WAVECAR").write_bytes(b"P" * (2 * 1024 * 1024))
    (pbe / "OUTCAR").write_text("   NBANDS=     24\n" + DONE)
    (pbe / "workflow_method.txt").write_text("DIM=3d\nFUNC=pbe\n")
    labels = ["G"] + [""] * 8 + ["X"] + [""] * 8 + ["M"]
    pts = [[0.5 * i / 18, 0.5 * min(i, 9) / 9 * (i > 9), 0.0] for i in range(19)]
    (pbe / "kpath.json").write_text(json.dumps({"kpoints": pts, "point_labels": labels, "breaks": []}))
    (pbe / "KPOINTS_OPT").write_text("path\n%d\nReciprocal\n" % len(pts)
                                     + "".join("%.6f %.6f %.6f 1\n" % tuple(p) for p in pts))
    for t in ("submit_std_3d.tpl", "submit_ncl_3d.tpl"):
        (m / t).write_text(SUBMIT)
    return m


def _gen(mat, *args):
    # discriminant_common.py 由 autozt 按 gen_need 推到材料目录；这里用 PYTHONPATH 代替
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(KE), str(REPO), str(DISC)]))
    return subprocess.run([sys.executable, str(GEN), *args], cwd=str(mat),
                          capture_output=True, text=True, env=env, timeout=120)


def _incar(path):
    out = {}
    for ln in Path(path).read_text().splitlines():
        if "=" in ln:
            k, v = ln.split("=", 1)
            out[k.strip().upper()] = v.split("#")[0].strip()
    return out


# ------------------------------------------------------------------ 生成端
def test_stage1_scf_dir(tmp_path):
    m = _material(tmp_path)
    r = _gen(m, "--stage", "scf")
    assert r.returncode == 0, r.stdout + r.stderr
    d = m / "step2_bandgap" / "step2.28_hse_scf"
    assert d.is_dir() and not list(d.glob("*p*of*"))           # 单目录
    assert not (d / "KPOINTS_OPT").exists()                     # NKRED 与 KPOINTS_OPT 不兼容
    inc = _incar(d / "INCAR")
    assert inc["LHFCALC"] == ".TRUE." and inc["ISTART"] == "1" and inc["LWAVE"] == ".TRUE."
    assert (inc["NKREDX"], inc["NKREDY"], inc["NKREDZ"]) == ("2", "2", "1")   # 4 4 3：只有偶数轴降
    assert "NCORE" not in inc and "IVDW" not in inc
    assert (d / "WAVECAR").stat().st_size == 2 * 1024 * 1024    # 从 PBE 拷
    assert "_scf" in (d / "submit.sh").read_text()
    assert not (m / "step2_bandgap" / "step2.3_hse").exists()   # 第 1 段不碰第 2 段


def test_stage1_nkred_off_by_step_conf(tmp_path):
    m = _material(tmp_path)
    (m / "step.conf").write_text("[params]\nHSE_SCF_NKRED = off\n")
    assert _gen(m, "--stage", "scf").returncode == 0
    inc = _incar(m / "step2_bandgap" / "step2.28_hse_scf" / "INCAR")
    assert not any(k.startswith("NKRED") for k in inc)


def _finish_stage1(m, nbands=32):
    d = m / "step2_bandgap" / "step2.28_hse_scf"
    (d / "OUTCAR").write_text("   NBANDS=     %d\n" % nbands + DONE)
    (d / "WAVECAR").write_bytes(b"H" * (3 * 1024 * 1024))       # HSE 波函数（与 PBE 那份区分开）
    return d


def test_stage2_restarts_from_hse_wavecar(tmp_path):
    m = _material(tmp_path)
    assert _gen(m, "--stage", "scf").returncode == 0
    _finish_stage1(m)
    r = _gen(m, "--kpath-slice", "all/4")
    assert r.returncode == 0, r.stdout + r.stderr
    parts = sorted((m / "step2_bandgap" / "step2.3_hse").glob("p*of4"))
    assert len(parts) == 4
    for p in parts:
        assert (p / "WAVECAR").read_bytes()[:1] == b"H"             # 第 1 段的 HSE WAVECAR
        assert (p / "KPOINTS_OPT").is_file()
        inc = _incar(p / "INCAR")
        assert not any(k.startswith("NKRED") for k in inc)          # KPOINTS_OPT + NKRED 必崩
        assert inc["ISTART"] == "1" and int(inc["NBANDS"]) % 32 == 0 and int(inc["NBANDS"]) >= 32
        src = json.loads((p / "hse_scf_source.json").read_text())
        assert src["wavecar_from"] == "step2_bandgap/step2.28_hse_scf"
        assert (p / "KPOINTS").read_text() == (m / "step2_bandgap" / "step2.28_hse_scf" / "KPOINTS").read_text()


def test_stage2_waits_for_unfinished_stage1(tmp_path):
    m = _material(tmp_path)
    assert _gen(m, "--stage", "scf").returncode == 0                # 目录在、没跑完
    r = _gen(m, "--kpath-slice", "all/4")
    assert r.returncode != 0 and "S2.28_hsescf" in (r.stdout + r.stderr)
    assert not (m / "step2_bandgap" / "step2.3_hse").exists()


def test_stage2_legacy_fallback_without_stage1_dir(tmp_path):
    m = _material(tmp_path)
    r = _gen(m, "--kpath-slice", "all/4")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "退回旧口径" in r.stderr
    p = m / "step2_bandgap" / "step2.3_hse" / "p1of4"
    assert (p / "WAVECAR").read_bytes()[:1] == b"P"                 # PBE WAVECAR（旧口径）
    assert not (p / "hse_scf_source.json").exists()


# ------------------------------------------------------------------ [0066] ALGO 来源
def _eig_outcar(gap, nelect=20, nbands=24):
    """两个 k 点的本征值块：占据带顶 -5.0 eV（k1），导带底 -5.0+gap（k2）。"""
    nocc = nelect // 2
    out = ["   NELECT =      %d.0000    total number of electrons" % nelect,
           "   NBANDS=     %d" % nbands]
    for k, (kv, shift) in enumerate((("0.0000    0.0000    0.0000", 0.0),
                                     ("0.5000    0.0000    0.0000", 0.3)), 1):
        out.append(" k-point     %d :       %s" % (k, kv))
        for b in range(1, nbands + 1):
            if b <= nocc:
                e, f = -10.0 + 0.5 * b - shift, "2.00000"
            else:
                e, f = -5.0 + gap + 0.5 * (b - nocc - 1) + (0.3 - shift), "0.00000"
            out.append("     %d     %.4f      %s" % (b, e, f))
    return "\n".join(out) + "\n" + DONE


def _algo(m):
    return _incar(m / "step2_bandgap" / "step2.28_hse_scf" / "INCAR")["ALGO"]


def test_algo_falls_back_to_pbe_gap(tmp_path):
    m = _material(tmp_path)
    (m / "step2_bandgap" / "step2.2_pbe" / "OUTCAR").write_text(_eig_outcar(0.0645))   # Mg2C60 的 PBE 隙
    r = _gen(m, "--stage", "scf")
    assert r.returncode == 0, r.stdout + r.stderr
    assert _algo(m) == "All" and "S2.2 PBE" in r.stderr


def test_algo_pbe_overlap_stays_damped(tmp_path):
    m = _material(tmp_path)
    (m / "step2_bandgap" / "step2.2_pbe" / "OUTCAR").write_text(_eig_outcar(-0.2))
    assert _gen(m, "--stage", "scf").returncode == 0
    assert _algo(m) == "Damped"


def test_algo_discriminant_wins_over_pbe(tmp_path):
    m = _material(tmp_path)
    (m / "step2_bandgap" / "step2.2_pbe" / "OUTCAR").write_text(_eig_outcar(0.5))
    d = m / "step2_bandgap" / "step2.15_discriminant"
    d.mkdir(parents=True)
    (d / "OUTCAR").write_text(_eig_outcar(-0.1))                    # 密网格：带重叠
    r = _gen(m, "--stage", "scf")
    assert r.returncode == 0 and _algo(m) == "Damped" and "S2.15" in r.stderr


def test_algo_unfinished_discriminant_falls_back(tmp_path):
    m = _material(tmp_path)
    (m / "step2_bandgap" / "step2.2_pbe" / "OUTCAR").write_text(_eig_outcar(0.5))
    d = m / "step2_bandgap" / "step2.15_discriminant"
    d.mkdir(parents=True)
    (d / "OUTCAR").write_text("   NBANDS=     24\n running\n")        # 排队/OOM 留下的半截 OUTCAR
    r = _gen(m, "--stage", "scf")
    assert r.returncode == 0 and _algo(m) == "All"


# ------------------------------------------------------------------ 判据
@pytest.fixture()
def ck():
    src = (REPO / "skill" / "_common" / "opt" / "checks_relax.py").read_text(encoding="utf-8")
    ns = {}
    exec(compile(src, "checks_relax", "exec"), ns)
    return ns["CHECKERS"]["hse_scf"]


def test_check_stage1(tmp_path, ck):
    d = tmp_path / "step2_bandgap" / "step2.28_hse_scf"
    d.mkdir(parents=True)
    assert ck(str(d), {}) == (False, "OUTCAR missing")
    (d / "INCAR").write_text("x\n")
    (d / "OUTCAR").write_text("running\n")
    assert ck(str(d), {})[0] is False
    (d / "OUTCAR").write_text(DONE)
    assert ck(str(d), {})[0] is False                               # 收尾了但没有 WAVECAR
    (d / "WAVECAR").write_bytes(b"x" * (2 * 1024 * 1024))
    assert ck(str(d), {})[0] is True


def test_check_legacy_skip(tmp_path, ck):
    base = tmp_path / "step2_bandgap"
    d = base / "step2.28_hse_scf"
    d.mkdir(parents=True)
    for i in (1, 2):
        p = base / "step2.3_hse" / ("p%dof2" % i)
        p.mkdir(parents=True)
        (p / "OUTCAR").write_text(DONE)
    ok, why = ck(str(d), {})
    assert ok and "旧口径" in why                                    # 老材料：不再交一个昂贵的 SCF
    (base / "step2.3_hse" / "p2of2" / "OUTCAR").write_text("running\n")
    assert ck(str(d), {})[0] is False                               # 没算完的旧 S2.3 不算
    (base / "step2.3_hse" / "p2of2" / "OUTCAR").write_text(DONE)
    (base / "step2.3_hse" / "p1of2" / "hse_scf_source.json").write_text("{}")
    assert ck(str(d), {})[0] is False                               # 两段式产物不当旧口径
    (base / "step2.3_hse" / "p1of2" / "hse_scf_source.json").unlink()
    (d / "INCAR").write_text("x\n")
    assert ck(str(d), {})[0] is False                               # 本步已生成过：必须真跑完


# ------------------------------------------------------------------ 接线
def test_skill_yaml_and_cascade():
    import yaml
    man = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    steps = {s["label"]: s for s in man["optional_steps"]["bandgap_hse"]["steps"]}
    s1, s2 = steps["S2.28_hsescf"], steps["S2.3_hse"]
    assert s1["check"] == "hse_scf" and "--stage scf" in s1["gen"] and "fanout" not in s1
    assert "discriminant_common.py" in s1["gen_need"] and "discriminant_common.py" in s2["gen_need"]
    assert "step2_bandgap/step2.28_hse_scf" in s2["needs"]
    assert s1["seq"] < s2["seq"]
    sys.path.insert(0, str(KE))
    import ke_common as kc
    assert ("step2_bandgap/step2.3_hse", "dir") in kc.DOWNSTREAM["step2_bandgap/step2.28_hse_scf"]
    assert ("step2_bandgap/step2.28_hse_scf", "dir") in kc.DOWNSTREAM["step2_bandgap/step2.2_pbe"]
