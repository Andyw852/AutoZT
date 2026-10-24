"""[0080] 用户 2026-10-24：横轴是掺杂浓度的图，近本征段（n ≤ 1e15）与掺杂段（n ≥ 1e16，三维等效）分开画；
化学势横轴的图画在一起；标题只写是什么 + 范围，图例不用 together。纯本地。"""
import importlib.util
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
sys.path.insert(0, str(ROOT / "tests"))
import test_0032_intrinsic_suite as T32  # noqa: E402

LOW_2D = [1e3, 1e5]                       # 面密度 cm^-2 -> 三维等效 ~1.5e10 / 1.5e12（近本征段）


@pytest.fixture
def mat(tmp_path, monkeypatch):
    monkeypatch.setattr(T32, "N2D", LOW_2D + T32.N2D)
    return T32._material(tmp_path)


def test_intrinsic_suite_split_by_range(mat):
    S = T32._s83()
    (mat / "step8.3_output").mkdir(exist_ok=True)
    files = S.write_intrinsic_suite(mat / "step8.3_output", mat, True)
    pngs = sorted(f for f in files if f.endswith(".png"))
    assert pngs == ["intrinsic_transport_300K.png", "intrinsic_transport_300K_low_doping.png",
                    "intrinsic_transport_600K.png", "intrinsic_transport_600K_low_doping.png"]


def _plotted_x(monkeypatch):
    import matplotlib.axes as ma
    xs, orig = [], ma.Axes.plot

    def rec(self, *a, **k):
        if a and not isinstance(a[0], str):
            try:
                xs.extend(float(v) for v in a[0])
            except TypeError:
                pass
        return orig(self, *a, **k)
    monkeypatch.setattr(ma.Axes, "plot", rec)
    return xs


def test_comparison_300k_split_by_range(mat, monkeypatch):
    S = T32._s83()
    out = mat / "step8.3_output"
    out.mkdir(exist_ok=True)
    geom = S._geom_2d(mat)
    am2 = S.load_amset(mat, True, S.AMSET2D_DIR)
    xs = _plotted_x(monkeypatch)
    S.make_figure(out, None, None, None, am2=am2, geom=geom, nrange=(0.0, S.PL_LOW_HI * 1.001),
                  fname="comparison_300K_low_doping.png", M="MoS2")
    t_cm = geom["t_A"] * 1e-8
    assert xs and max(xs) / t_cm <= S.PL_LOW_HI * 1.01                # x = n_2D；只有近本征段
    xs.clear()
    S.make_figure(out, None, None, None, am2=am2, geom=geom, nrange=(S.PL_TUNE_LO * 0.999, float("inf")), M="MoS2")
    assert xs and min(xs) / t_cm >= S.PL_TUNE_LO * 0.99
    assert (out / "comparison_300K.png").is_file() and (out / "comparison_300K_low_doping.png").is_file()


def test_titles_short_and_no_together():
    src = (KE / "step8.3_output" / "gen_step13_output.py").read_text(encoding="utf-8")
    body = src[src.index("def mech_tag"):]
    lines = [ln for ln in body.splitlines() if not ln.lstrip().startswith("#")]
    assert not any("together" in ln for ln in lines if "label" in ln or "suptitle" in ln or "_place" in ln)
    for old in ("4 models:", "%d models: ", "chemical-potential scan, %.0f K: BT2", "mechanism-resolved mobility%s",
                "Electronic transport comparison @ ~300 K", "AMSET transport by scattering mechanisms @"):
        assert old not in src, old
