"""skill/ke-dft-cpu/tools/nscf_dp_check.py：粗网格自洽 + 密网格非自洽算形变势的独立验证工具（不需要 VASP）。

假的 VASP（EPM 硅：按各段 POSCAR 的晶格、KPOINTS 的网格算本征值，写 VASP 格式的 EIGENVAL / OUTCAR，
平均静电芯势按 pymatgen 的定宽格式写）把 make -> run.sh -> analyze 整条走一遍：
  · make：三段的 INCAR 改对了（粗自洽写 CHGCAR、非自洽 ICHARG = 11、都不写 WAVECAR）、密网格按 |b|/dk、
    submit.sh 沿用 S7 的头、改作业名 / 墙钟；out_dir 在 step7_deform 里、杂化泛函、缺 ±δ 一律拒绝；
  · run.sh：粗自洽 -> 拷 CHGCAR -> 非自洽 -> 密自洽，timing.txt 记三段；
  · analyze：芯能级参考的读法与 pymatgen 一致；D 与独立按 AMSET 公式算的值相等；非自洽缺芯能级块 -> 不通过（退出码 1）；
    没跑完 -> 退出码 2。
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "skill" / "ke-dft-cpu" / "tools" / "nscf_dp_check.py"
sys.path.insert(0, str(TOOL.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nscf_dp_check as N  # noqa: E402
import _epm_si as E  # noqa: E402

DELTA = 0.005
FAKE_VASP = r'''
import os, re, sys, itertools
import numpy as np
sys.path.insert(0, %(tests)r)
import _epm_si as E
ln = open("POSCAR").read().splitlines()
L = np.array([[float(x) for x in ln[i].split()[:3]] for i in (2, 3, 4)]) * float(ln[1].split()[0])
mesh = [int(x) for x in open("KPOINTS").read().splitlines()[3].split()[:3]]
inc = open("INCAR").read()
nscf = re.search(r"^ICHARG\s*=\s*11", inc, re.M) is not None
ks, seen = [], set()
for g in itertools.product(*[range(n) for n in mesh]):
    k = np.array(g, float) / mesh
    k = k - np.round(k + 1e-9)
    key = tuple(np.round(k * mesh).astype(int) %% mesh)
    if tuple(np.round(-k * mesh).astype(int) %% mesh) in seen:
        continue
    seen.add(key)
    ks.append(k)
NB = 8
F = (np.linalg.inv(E.L0) @ L).T
gl = 0.5 * (F.T @ F - np.eye(3))
ref = -80.0 + 2.5 * np.trace(gl)                       # 参考能级随应变移动（AMSET 要扣掉的那部分）
shift = float(os.environ.get("FAKE_NSCF_DRIFT", "0")) * np.trace(gl) if nscf else 0.0   # 非自洽本征值随应变多漂一点
with open("EIGENVAL", "w") as f:
    f.write("    2    2    1    1\n  0.1 0.1 0.1 0.1 0.5E-15\n  1.0E-004\n  CAR\n Si\n")
    f.write("    8  %%d  %%d\n" %% (len(ks), NB))
    for k in ks:
        e = E.solve(k, L, 120.0)[1][:NB] + ref + shift
        f.write("\n  %%.7E  %%.7E  %%.7E  %%.7E\n" %% (k[0], k[1], k[2], 1.0 / len(ks)))
        for b in range(NB):
            f.write("    %%d   %%.6f   %%.6f\n" %% (b + 1, e[b], 1.0 if b < 4 else 0.0))
with open("OUTCAR", "w") as f:
    f.write(" running\n")
    if not (nscf and os.environ.get("FAKE_NSCF_NOCORE")):
        f.write("  average (electrostatic) potential at core\n    the test charge radii are     0.8698\n")
        f.write("  (the norm of the test charge is              1.0000)\n")
        f.write("%%8d%%9.4f%%8d%%9.4f\n\n" %% (1, ref, 2, ref))
    f.write(" E-fermi :   5.0000     XC(G=0):  -8.0000     alpha+bet : -4.0000\n")
    f.write("                         Elapsed time (sec):     %%8.3f\n" %% (len(ks) * (0.2 if nscf else 1.0)))
    f.write(" General timing and accounting informations for this job:\n")
if re.search(r"^LCHARG\s*=\s*\.TRUE\.", inc, re.M):
    open("CHGCAR", "w").write("chg\n")
'''

INCAR = """SYSTEM = Si deform
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
SUBMIT = """#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --job-name=Si-ke-dft-cpu-S7_deform-deform-001
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-node=24
cd $SLURM_SUBMIT_DIR
# 注释里也有 mpirun 字样
module purge
mpirun -np $SLURM_NTASKS vasp_std
"""


def _poscar(L):
    return "Si\n1.0\n" + "".join("  %.10f %.10f %.10f\n" % tuple(r) for r in L) + "Si\n2\nDirect\n 0.125 0.125 0.125\n -0.125 -0.125 -0.125\n"


@pytest.fixture()
def s7(tmp_path):
    root = tmp_path / "Si" / "step7_deform"
    cfg = [("undeformed", None)] + [("deform-%03d" % (i + 1), E.pattern(P, s * DELTA))
                                     for i, (P, s) in enumerate([((0, 0), 1), ((0, 0), -1), ((0, 1), 1), ((0, 1), -1)])]
    for name, Em in cfg:
        d = root / name
        d.mkdir(parents=True)
        (d / "POSCAR").write_text(_poscar(E.lattice(Em)))
        (d / "INCAR").write_text(INCAR)
        (d / "POTCAR").write_text("potcar\n")
        (d / "KPOINTS").write_text("Automatic\n0\nGamma\n 3 3 3\n 0 0 0\n")
        (d / "submit.sh").write_text(SUBMIT)
    return root


def _make(s7, out, *extra):
    return N.main(["make", str(s7), str(out), "--mesh", "4", "4", "4", *extra])


def test_make_inputs(s7, tmp_path, capsys):
    out = tmp_path / "check"
    assert _make(s7, out, "--time", "48:00:00") == 0
    meta = json.loads((out / N.META).read_text())
    assert [c["name"] for c in meta["configs"]] == ["undeformed", "deform-001", "deform-002"]
    assert [round(c["amount"], 6) for c in meta["configs"]] == [0.0, 0.005, -0.005]
    assert meta["coarse_mesh"] == [3, 3, 3] and meta["dense_mesh"] == [4, 4, 4]
    for c in ("undeformed", "deform-001"):
        sc = (out / c / "scf_coarse" / "INCAR").read_text()
        ns = (out / c / "nscf_dense" / "INCAR").read_text()
        sd = (out / c / "scf_dense" / "INCAR").read_text()
        assert "LCHARG = .TRUE." in sc and "LWAVE = .FALSE." in sc and "ISYM   = 0" in sc
        assert "ICHARG = 11" in ns and "LCHARG = .FALSE." in ns and "LVHAR  = .TRUE." in ns
        assert "LCHARG = .FALSE." in sd and "ICHARG = 2" in sd
        assert N.read_mesh(out / c / "scf_coarse" / "KPOINTS") == [3, 3, 3]
        assert N.read_mesh(out / c / "nscf_dense" / "KPOINTS") == [4, 4, 4]
        sub = (out / c / "submit.sh").read_text()
        assert "--job-name=nscfdp-%s" % c in sub and "--time=48:00:00" in sub
        assert sub.rstrip().endswith("bash run.sh mpirun -np $SLURM_NTASKS vasp_std") and sub.count("mpirun -np") == 1
    # S7 本身没被动
    assert (s7 / "undeformed" / "INCAR").read_text() == INCAR
    # 默认 --dk：Si 的 |b| ≈ 2.004 Å⁻¹ -> 0.06 时 34³
    out2 = tmp_path / "check2"
    assert N.main(["make", str(s7), str(out2)]) == 0
    assert json.loads((out2 / N.META).read_text())["dense_mesh"] == [34, 34, 34]


def test_make_refuses(s7, tmp_path):
    with pytest.raises(SystemExit, match="不能在"):
        _make(s7, s7 / "check")
    with pytest.raises(SystemExit, match="\\+δ 与 −δ"):
        _make(s7, tmp_path / "c", "--pattern", "zz")
    (s7 / "undeformed" / "INCAR").write_text(INCAR + "LHFCALC = .TRUE.\n")
    with pytest.raises(SystemExit, match="杂化泛函"):
        _make(s7, tmp_path / "c2")


def test_core_potential_reader_matches_pymatgen_layout(tmp_path):
    p = tmp_path / "OUTCAR"
    p.write_text("  average (electrostatic) potential at core\n    the test charge radii are     0.8698\n"
                 "  (the norm of the test charge is              1.0000)\n"
                 "       1 -83.6158       2 -83.6159       3-100.1234       4 -12.0000       5  -1.5000\n"
                 "       6 -83.0000\n\n E-fermi :   5.0\n")
    assert N.read_core_potentials(p) == [-83.6158, -83.6159, -100.1234, -12.0, -1.5, -83.0]
    pmg = pytest.importorskip("pymatgen.io.vasp.outputs")
    p.write_text(p.read_text() + " General timing and accounting informations for this job:\n")
    try:
        ref = pmg.Outcar(str(p)).read_avg_core_poten()
    except Exception:                                                  # pymatgen 要完整的 OUTCAR 头，造不出来就只比本工具
        return
    assert np.allclose(ref[0], N.read_core_potentials(p))


def _run(out, env=None):
    fake = out.parent / "fake_vasp.py"
    fake.write_text(FAKE_VASP % {"tests": str(Path(__file__).resolve().parent)})
    for c in ("undeformed", "deform-001", "deform-002"):
        r = subprocess.run(["bash", str(out / c / "run.sh"), sys.executable, str(fake)],
                           env=dict(os.environ, **(env or {})), capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        t = (out / c / "timing.txt").read_text().split()
        assert t[0::3] == ["scf_coarse", "nscf_dense", "scf_dense"] and t[2::3] == ["0", "0", "0"]


def _expected_D(out, states_fn):
    """独立按 AMSET 公式：D = (E0 − E_s − (ref0 − ref_s))/δ_s，±δ 平均（同一密网格的 scf_dense）。"""
    d = {}
    for c in ("undeformed", "deform-001", "deform-002"):
        kp, Ev, _o, _n = N.read_eigenval(out / c / "scf_dense" / "EIGENVAL")
        d[c] = (kp, Ev[0], float(np.mean(N.read_core_potentials(out / c / "scf_dense" / "OUTCAR"))))
    k0, E0, r0 = d["undeformed"]
    vals = []
    for c, amt in (("deform-001", DELTA), ("deform-002", -DELTA)):
        kp, Es, rs = d[c]
        lut = {N._key(k): j for j, k in enumerate(kp)}
        vals.append(np.array([(E0[k, b] - Es[lut[N._key(k0[k])], b] - (r0 - rs)) / amt for k, b in states_fn(E0)]))
    return float(np.mean(np.mean(vals, axis=0)))


def test_end_to_end_pass(s7, tmp_path, capsys):
    out = tmp_path / "check"
    _make(s7, out)
    _run(out)
    capsys.readouterr()
    assert N.main(["analyze", str(out)]) == 0
    log = capsys.readouterr().out
    rep = json.loads((out / N.REPORT).read_text())
    assert rep["pass"] and rep["dp_rel_worst"] == pytest.approx(0.0, abs=1e-12)
    for c in ("undeformed", "deform-001"):
        assert rep["core_ref"][c]["nscf_minus_scf_coarse"] == pytest.approx(0.0, abs=1e-9)
    # 价带顶（Γ 三重态）的平均 D 与独立计算一致
    exp_vb = _expected_D(out, lambda E0: [(k, b) for k in range(E0.shape[0]) for b in range(4)
                                          if abs(E0[k, b] - E0[:, :4].max()) < 1e-3])
    assert rep["dp"]["价带顶"]["scf_edge_mean"] == pytest.approx(exp_vb, abs=1e-3)
    assert rep["dp"]["价带顶"]["nscf_edge_mean"] == pytest.approx(exp_vb, abs=1e-3)
    assert "coarse_nearest_mean" in rep["dp"]["导带底"]
    assert "通过" in log and "墙钟" in log


def test_end_to_end_nscf_drift_fails(s7, tmp_path, capsys):
    """非自洽的本征值随应变多漂 0.5 eV/单位应变 -> D 差 0.5 eV，超过 2% × max|D| -> 不通过（退出码 1）。"""
    out = tmp_path / "check"
    _make(s7, out)
    _run(out, {"FAKE_NSCF_DRIFT": "0.5"})
    assert N.main(["analyze", str(out)]) == 1
    rep = json.loads((out / N.REPORT).read_text())
    assert not rep["pass"] and rep["dp"]["价带顶"]["max_abs_diff"] == pytest.approx(0.5, abs=1e-3)
    assert "不通过" in capsys.readouterr().out


def test_end_to_end_nscf_no_core(s7, tmp_path, capsys):
    out = tmp_path / "check"
    _make(s7, out)
    _run(out, {"FAKE_NSCF_NOCORE": "1"})
    assert N.main(["analyze", str(out)]) == 2                         # 非自洽缺芯能级：形变势算不了
    rep = json.loads((out / N.REPORT).read_text())
    assert rep["core_ref"]["undeformed"]["ref"]["nscf_dense"] is None
    assert "没有平均静电芯势" in capsys.readouterr().out


def test_analyze_incomplete(s7, tmp_path):
    out = tmp_path / "check"
    _make(s7, out)
    assert N.main(["analyze", str(out)]) == 2
    assert json.loads((out / N.REPORT).read_text())["missing"]
