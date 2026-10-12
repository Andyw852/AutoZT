"""[0040] 本征载流子（未掺杂，doping = 0）：
1. amset_fermi_fix 的本征 E_F 解 n = p（log-sum-exp），对抛物带与解析式一致；AMSET 原搜索在宽带隙低温下偏；
2. S8.3 transport_plots/intrinsic_carrier/<M>_intrinsic_carrier.{png,csv}：2D（σ、κe × c/t）与 3D；
   没有 0 掺杂点时跳过并清掉旧图；
3. S8 / S8.4 的 gen 有 UNDOPED_POINT（出厂开）。"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _load(name, path):
    for d in (ROOT, KE, KE / "step8.4_amset2d", ROOT / "skill" / "_common" / "opt", ROOT / "tests"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T32 = _load("t0032_for_0040", ROOT / "tests" / "test_0032_intrinsic_suite.py")


def _s83():
    return _load("s83_0040", KE / "step8.3_output" / "gen_step13_output.py")


def _set_undoped_sigma(cwd, value):
    p = cwd / "step8.4_amset2d" / "transport.json"
    j = json.loads(p.read_text())
    i0 = j["doping"].index(0.0)
    for k in range(len(j["temperatures"])):
        j["conductivity"][i0][k] = T32._t(value)
        j["seebeck"][i0][k] = T32._t(-150.0)
    txt = json.dumps(j)
    p.write_text(txt)
    import hashlib                      # 子集文件记的是 transport.json 指纹：改了内容要同步，否则被当成上一轮残留
    sha = hashlib.sha256(txt.encode()).hexdigest()
    for f in (cwd / "step8.4_amset2d").glob("intrinsic_*.json"):
        q = json.loads(f.read_text())
        q["transport_fingerprint"]["sha256"] = sha
        f.write_text(json.dumps(q))


def test_undoped_rows_2d_scaled_by_c_over_t(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path, doping0=True)
    _set_undoped_sigma(cwd, 2.0)
    rows, why = S.undoped_rows(cwd, True)
    assert why is None
    geom = S._geom_2d(cwd)
    ov = [r for r in rows if r["dataset"] == "overall"]
    assert len(ov) == len(T32.TEMPS)
    assert ov[0]["sigma_S/m"] == pytest.approx(2.0 * geom["c_over_t"])
    assert ov[0]["S_uV/K"] == pytest.approx(-150.0)
    assert ov[0]["PF_W/mK2"] == pytest.approx((150e-6) ** 2 * 2.0 * geom["c_over_t"])
    assert {r["dataset"] for r in rows} >= {"overall", "intrinsic", "ADP_only"}


def test_undoped_plot_files(tmp_path):
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = T32._material(tmp_path, doping0=True)
    _set_undoped_sigma(cwd, 2.0)
    plots = cwd.parent / "transport_plots"
    assert S.write_undoped_plots(cwd, plots, "MoS2", True)
    d = plots / S.UNDOPED_DIR
    assert (d / "MoS2_intrinsic_carrier.png").is_file()
    lines = [l for l in (d / "MoS2_intrinsic_carrier.csv").read_text().splitlines() if not l.startswith("#")]
    rows = list(csv.DictReader(lines))
    assert {r["dataset"] for r in rows} >= {"overall", "intrinsic", "ADP_only"}


def test_no_undoped_point_skips_and_cleans(tmp_path):
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = T32._material(tmp_path, doping0=False)
    d = cwd.parent / "transport_plots" / S.UNDOPED_DIR
    d.mkdir(parents=True)
    (d / "MoS2_intrinsic_carrier.png").write_bytes(b"old")
    assert not S.write_undoped_plots(cwd, cwd.parent / "transport_plots", "MoS2", True)
    assert not (d / "MoS2_intrinsic_carrier.png").exists()


def test_undoped_rows_3d(tmp_path):
    S = _s83()
    cwd = tmp_path / "Si" / "ke-dft-cpu"
    (cwd / "step8_amset").mkdir(parents=True)
    t3 = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 6.0]]
    (cwd / "step8_amset" / "transport.json").write_text(json.dumps({
        "doping": [-1e18, 0.0, 1e18], "temperatures": [300.0, 600.0],
        "seebeck": [[t3, t3]] * 3, "conductivity": [[t3, t3]] * 3,
        "electronic_thermal_conductivity": [[t3, t3]] * 3}))
    rows, why = S.undoped_rows(cwd, False)
    assert why is None and len(rows) == 2
    assert rows[0]["sigma_S/m"] == pytest.approx(4.0)          # (3+3+6)/3，3D 不乘 c/t


def test_gen_has_undoped_point():
    for p in (KE / "step8.4_amset2d" / "gen_step14_amset2d.py", KE / "step8_amset" / "gen_step10_amset.py"):
        src = p.read_text(encoding="utf-8")
        assert "UNDOPED_POINT = True" in src
        assert '"UNDOPED_POINT": (None, "bool")' in src
        assert "if UNDOPED_POINT and 0.0 not in _dop:" in src


def test_intrinsic_fermi_matches_parabolic_bands():
    np = pytest.importorskip("numpy")
    pytest.importorskip("amset")
    try:
        import spglib as _sg                     # 新 spglib 去掉了 spglib.spglib，老 amset 还在 import 它
        if not hasattr(_sg, "spglib"):
            _sg.spglib = _sg
            sys.modules["spglib.spglib"] = _sg
    except ImportError:
        pass
    F = _load("fermi_fix_0040", KE / "step8.4_amset2d" / "amset_fermi_fix.py")
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Spin
    from amset.constants import boltzmann_au, hartree_to_ev
    from amset.electronic_structure.dos import FermiDos
    st = Structure(Lattice.cubic(5.0), ["Si"], [[0, 0, 0]])
    ha = 1 / hartree_to_ev
    E = np.linspace(-3, 3, 60001) * ha
    Ev, Ec, Av, Ac = -0.9 * ha, 0.9 * ha, 3.0, 1.0
    dens = (np.where(E < Ev, Av * np.sqrt(np.clip(Ev - E, 0, None)), 0)
            + np.where(E > Ec, Ac * np.sqrt(np.clip(E - Ec, 0, None)), 0))
    dos = FermiDos(0.0, E, {Spin.up: dens / 2}, st, atomic_units=True)
    for T in (300.0, 900.0):
        ref = 0.5 * (Ev + Ec) + 0.5 * T * boltzmann_au * np.log(Av / Ac)
        assert F.intrinsic_fermi(dos, T) == pytest.approx(ref, abs=1e-6)
        assert dos.get_fermi_from_num_electrons(dos.nelect, T, tol=1e-8) == pytest.approx(ref, abs=1e-6)
    orig = F.unwrap(FermiDos.get_fermi_from_num_electrons)
    ef0 = orig(dos, dos.nelect, 300.0, tol=1e-11)
    assert abs(ef0 - (0.5 * (Ev + Ec) + 0.5 * 300 * boltzmann_au * np.log(3.0))) * hartree_to_ev > 0.05
    # 金属（无带隙）：照旧走原搜索
    metal = FermiDos(0.0, E, {Spin.up: np.full_like(E, 1.0)}, st, atomic_units=True)
    assert F.band_gap_edges(metal) is None
