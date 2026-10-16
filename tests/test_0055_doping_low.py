"""[0055] 低浓度段 DOPING_LOW（3D 1e10–1e16 cm^-3，2D 1e3–1e9 cm^-2）：
1. amset_fermi_fix.small_doping_fermi：n、p 分开 log-sum-exp 解 p − n = d，每胞 4e-13 个载流子也准；
2. get_fermi 外壳：每胞 < SMALL_DOPING_PER_CELL 走小掺杂解，否则原样交给 AMSET（出厂网格逐位不变）；
3. S8 / S8.4 gen 的 DOPING_LOW 出厂值与 off；
4. S8.3 浓度轴：没有低浓度点时与旧版相同，有时原点下探。"""
import importlib.util
import inspect
import json
import math
import sys
import types
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


def _ff():
    pytest.importorskip("scipy")
    return _load("fermi_fix_0055", KE / "step8.4_amset2d" / "amset_fermi_fix.py")


def _dos(F, Eg_ev, Av=3.0, Ac=1.0, vol=300.0):
    np = pytest.importorskip("numpy")
    ha = 1 / F.HA_EV
    E = np.linspace(-3, 3, 120001) * ha
    Ev, Ec = -Eg_ev / 2 * ha, Eg_ev / 2 * ha
    dens = (np.where(E < Ev, Av * np.sqrt(np.clip(Ev - E, 0, None)), 0)
            + np.where(E > Ec, Ac * np.sqrt(np.clip(E - Ec, 0, None)), 0))
    return types.SimpleNamespace(energies=E, tdos=dens, efermi=0.0, atomic_units=True,
                                 structure=types.SimpleNamespace(volume=vol)), Ev, Ec


@pytest.mark.parametrize("Eg,T", [(1.8, 300.0), (1.8, 900.0), (0.3, 300.0), (0.3, 900.0)])
@pytest.mark.parametrize("d_cell", [-4e-13, 4e-13, -4e-9, 4e-9])
def test_small_doping_fermi_exact_net_charge(Eg, T, d_cell):
    F = _ff()
    dos, Ev, Ec = _dos(F, Eg)
    ef, n, p = F.small_doping_fermi(dos, d_cell / 300.0, T)
    assert Ev < ef < Ec
    assert (p - n) * 300.0 / d_cell == pytest.approx(1.0, rel=1e-6)   # 窄带隙 900 K（n_i/|d| ~ 1e8）也成立


def test_small_doping_fermi_matches_boltzmann_analytic():
    F = _ff()
    dos, Ev, Ec = _dos(F, 1.8)
    T = 300.0
    kt = T * F.KB_AU
    Nc = 1.0 * math.sqrt(math.pi) / 2 * kt ** 1.5
    ef, n, p = F.small_doping_fermi(dos, -4e-13 / 300.0, T)
    assert (ef - (Ec + kt * math.log(4e-13 / Nc))) * F.HA_EV == pytest.approx(0.0, abs=1e-5)


def test_get_fermi_wrapper_routes_small_doping_only():
    """外壳：每胞 < 1e-7 走 log-sum-exp；大掺杂原样调用 AMSET 原搜索（不改出厂网格结果）。"""
    F = _ff()
    dos, Ev, Ec = _dos(F, 1.8)
    calls = []

    class FakeFermiDos:
        energies, tdos, efermi, atomic_units, structure = (dos.energies, dos.tdos, 0.0, True, dos.structure)

        def get_fermi(self, concentration, temperature, tol=0.01, nstep=50, step=0.1, precision=10,
                      return_electron_hole_conc=False):
            calls.append(concentration)
            return (0.123, 1.0, 2.0) if return_electron_hole_conc else 0.123

        def get_doping(self, ef, T, return_electron_hole_conc=False):
            return 0.0

    orig = FakeFermiDos.get_fermi
    FakeFermiDos.get_fermi = F._wrap(orig, inspect.signature(orig))
    fd = FakeFermiDos()
    small = -4e-13 / 300.0                               # 每胞 4e-13
    ef, n, p = fd.get_fermi(small, 300.0, tol=1e-5, return_electron_hole_conc=True)
    assert calls == [] and Ev < ef < Ec and (n - p) * 300.0 == pytest.approx(4e-13, rel=1e-6)
    big = -1e-5 / 300.0                                  # 每胞 1e-5（出厂 3D 1e17 量级）
    assert fd.get_fermi(big, 300.0, tol=1e-5) == 0.123 and calls == [big]


def test_gen_doping_low_defaults_and_off():
    for path, want in ((KE / "step8.4_amset2d" / "gen_step14_amset2d.py", [1e3, 1e5, 1e7, 1e9]),
                       (KE / "step8_amset" / "gen_step10_amset.py", [1e10, 1e12, 1e14, 1e16])):
        g = _load("gen_0055_%s" % path.stem, path)
        vals = g._doping_low_list(g.DOPING_LOW)
        assert sorted(v for v in vals if v > 0) == pytest.approx(want)
        assert sorted(-v for v in vals if v < 0) == pytest.approx(want)
        assert g._doping_low_list("off") == [] and g._doping_low_list("") == []
        src = path.read_text(encoding="utf-8")
        assert '"DOPING_LOW": (None, "str")' in src
        assert "_doping_low_list(DOPING_LOW)" in src


def test_fermi_fix_threshold_keeps_factory_grid_untouched():
    """出厂网格最小点每胞载流子 ≥ 1e-7（不进小掺杂分支）：3D 1e17 cm^-3 × 10 Å^3、2D 1e11 cm^-2 × 5 Å^2。"""
    F = _ff()
    assert 1e17 * 10e-24 >= F.SMALL_DOPING_PER_CELL
    assert 1e11 * 5e-16 >= F.SMALL_DOPING_PER_CELL


def _write_transport(cwd, dop, run="step8.4_amset2d"):
    (cwd / run).mkdir(parents=True, exist_ok=True)
    (cwd / run / "transport.json").write_text(json.dumps({"doping": dop, "temperatures": [300.0]}))


def test_doping_axis_unchanged_without_low_points(tmp_path):
    S = _load("s83_0055", KE / "step8.3_output" / "gen_step13_output.py")
    _write_transport(tmp_path, [-2.87e20, -2.87e17, 0.0, 2.87e17, 2.87e20])
    NLO, NFLOOR, XT, XTL = S._doping_axis(tmp_path, True)
    assert (NLO, NFLOOR) == (1e17, 5e17)
    assert XT == [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    assert XTL == ["1e21", "1e20", "1e19", "1e18", "1e17", "1e18", "1e19", "1e20", "1e21"]


def test_doping_axis_extends_with_low_points(tmp_path):
    S = _load("s83_0055b", KE / "step8.3_output" / "gen_step13_output.py")
    _write_transport(tmp_path, [-2.87e20, -2.87e9, 0.0, 2.87e9, 2.87e20])   # 2D 1e3 cm^-2 / c
    NLO, NFLOOR, XT, XTL = S._doping_axis(tmp_path, True)
    assert NLO == 1e9 and NFLOOR < 2.87e9
    assert XT[0] == -12 and XT[-1] == 12 and 0 in XT and XT == sorted(XT)
    assert XTL[0] == "1e21" and XTL[XT.index(0)] == "1e9"
