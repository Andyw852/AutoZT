"""[0077] 掺杂网格两段（2026-10-24 用户定）：低浓度段（近本征）±1e11–1e15、掺杂段 ±1e16–1e21 cm^-3，每十年一个点。
3D（S8）直接用体浓度；2D（S8.4）同两段按三维等效浓度定义、按层厚 t 换算：n_2D = n_3D × t，AMSET doping = n_2D / c。
S8.3 的 2D 浓度轴同一口径（n_cell × c/t），文件夹 carrier_tuning_1e16-1e21 / low_doping_1e11-1e15。纯本地。"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
LOW = [1e11, 1e12, 1e13, 1e14, 1e15]
HIGH = [1e16, 1e17, 1e18, 1e19, 1e20, 1e21]
C_A, T_A = 20.0, 6.5389                     # MoS2 单层：含真空胞高 / vdW 层厚


def _load(name, path):
    for d in (ROOT, KE, KE / "step8.4_amset2d", ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _pm(vals):
    return sorted([-v for v in vals] + list(vals))


def test_s8_factory_grid_one_point_per_decade():
    g = _load("s8_0077", KE / "step8_amset" / "gen_step10_amset.py")
    assert sorted(g.expand_spec(g.DOPING, log=True)) == pytest.approx(_pm(HIGH))
    assert sorted(g._doping_low_list(g.DOPING_LOW)) == pytest.approx(_pm(LOW))


def _s84():
    g = _load("s84_0077", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")
    assert g.DOPING is None and g.DOPING_LOW is None       # 出厂：面密度不写死，由三维等效 × t 换算
    return g


def test_s84_factory_grid_is_3d_equivalent_times_t():
    g = _s84()
    t_cm = T_A * 1e-8
    main, low, src = g.doping_areal_grid(t_cm)
    assert sorted(v / t_cm for v in main) == pytest.approx(_pm(HIGH))
    assert sorted(v / t_cm for v in low) == pytest.approx(_pm(LOW))
    assert "DOPING_3D" in src and "DOPING_LOW_3D" in src
    # AMSET doping（元胞体积口径）× c/t 回到三维等效的十年刻度 —— S8.3 浓度轴就是这个口径
    dop = sorted(v / (C_A * 1e-8) * (C_A / T_A) for v in main + low)
    assert dop == pytest.approx(_pm(LOW + HIGH))


def test_s84_explicit_areal_override_and_off():
    g = _s84()
    g.DOPING, g.DOPING_LOW = "-1e14:-1e11:4, 1e11:1e14:4", "off"     # 旧项目 step.conf 写的面密度照旧生效
    main, low, src = g.doping_areal_grid(T_A * 1e-8)
    assert sorted(main) == pytest.approx(_pm([1e11, 1e12, 1e13, 1e14]))
    assert low == [] and "DOPING（面密度" in src and "DOPING_LOW（面密度" in src
    g.DOPING, g.DOPING_LOW, g.DOPING_LOW_3D = None, None, "off"
    assert g.doping_areal_grid(T_A * 1e-8)[1] == []
    g.DOPING_LOW_3D = "-1e16:-1e15:2, 1e15:1e16:2"                   # 与掺杂段重合的点不重复
    main, low, _ = g.doping_areal_grid(T_A * 1e-8)
    assert sorted(v / (T_A * 1e-8) for v in low) == pytest.approx([-1e15, 1e15])


def test_s84_without_thickness_keeps_3d_values():
    """不是 2D（拿不到 t）时三维等效值原样当体浓度。"""
    g = _s84()
    main, low, _ = g.doping_areal_grid(None)
    assert sorted(main) == pytest.approx(_pm(HIGH)) and sorted(low) == pytest.approx(_pm(LOW))


def test_s84_carrier_guard_accepts_factory_grid():
    """出厂最大点 1e21 × t × 面内面积（MoS2 ~8.8e-16 cm²、3 原子）远低于 0.5/原子，不拦。"""
    g = _s84()
    main, low, _ = g.doping_areal_grid(T_A * 1e-8)
    g.check_carriers_per_atom(list(main) + list(low), 8.8e-16, 3)


def _material(tmp_path):
    """2D 材料：AMSET 网格 = 出厂两段换算到元胞口径（n_3D × t / c），本征表同一网格。"""
    pytest.importorskip("matplotlib")
    cwd = tmp_path / "MoS2" / "ke-dft-cpu"
    out = cwd / "step8.3_output"
    (out / "intrinsic").mkdir(parents=True)
    run = cwd / "step8.4_amset2d"
    run.mkdir()
    (run / "2d_correction.json").write_text(json.dumps(
        {"cell_c_A": C_A, "elastic_rescale_factor_c_over_t": C_A / T_A, "num_atoms": 3,
         "inplane_cell_area_cm2": 8.8e-16, "layer_thickness": {"method": "vdw"}}))
    grid = sorted([v * T_A / C_A for v in _pm(LOW + HIGH)] + [0.0])
    temps = [300.0, 600.0]

    def _t(f):
        return [[[[f(n), 0, 0], [0, f(n), 0], [0, 0, 0]] for _ in temps] for n in grid]
    mob = {"ADP": _t(lambda n: 500.0), "POP": _t(lambda n: 300.0), "IMP": _t(lambda n: 1e6),
           "overall": _t(lambda n: 180.0)}
    (run / "transport.json").write_text(json.dumps({"doping": grid, "temperatures": temps, "mobility": mob}))
    cols = ["T_K", "dataset", "type", "n_cell_cm-3", "n_2D_cm-2", "n_lit_cm-3", "carriers_per_cell",
            "carriers_per_atom", "in_interpret_range", "S_uV/K", "sigma_S/m", "PF_W/mK2", "kappa_e_W/mK",
            "mu_cm2/Vs", "flag"]
    with open(out / "intrinsic" / "intrinsic_transport_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for T in temps:
            for ds in ("overall", "intrinsic", "ADP_only", "BT2_DPTtau"):
                for n in grid:
                    w.writerow([T, ds, "n" if n < 0 else "p", "%.6g" % n, abs(n) * C_A * 1e-8,
                                abs(n) * C_A / T_A, "", "", "yes", -200 if n < 0 else 200,
                                1e4 * (abs(n) + 1) / 1e18, 1e-3, 0.1, 100, ""])
    return cwd, out


def _n3d(p):
    rows = list(csv.DictReader(p.read_text().splitlines()))
    return sorted({abs(float(r["n_cell_cm-3"])) * C_A / T_A for r in rows})


def test_s83_2d_axis_is_3d_equivalent(tmp_path):
    S = _load("s83_0077", KE / "step8.3_output" / "gen_step13_output.py")
    cwd, out = _material(tmp_path)
    fac = S._n_scale(cwd, True)[0]
    assert fac == pytest.approx(C_A / T_A)
    assert S._n_scale(cwd, False)[0] == 1.0                     # 3D：元胞 = 体浓度
    assert S._doping_top(cwd, True, fac) == 1e21
    NLO, _, XT, XTL = S._doping_axis(cwd, True, 1e21, fac)
    assert NLO == 1e11 and XTL[XT.index(0)] == "1e11" and XTL[0] == "1e21"
    S.write_transport_plots(cwd, out, True)
    tp = cwd / "transport_plots"
    assert sorted(p.name for p in tp.glob("transport_models_*")) == [
        "transport_models_carrier_tuning_1e16-1e21", "transport_models_low_doping_1e11-1e15"]
    hi = tp / "transport_models_carrier_tuning_1e16-1e21" / "MoS2_transport_models_carrier_tuning_300K.csv"
    lo = tp / "transport_models_low_doping_1e11-1e15" / "MoS2_transport_models_low_doping_300K.csv"
    assert _n3d(hi) == pytest.approx(HIGH, rel=1e-5)            # 表里仍是 n_cell，按 c/t 换回正好是十年刻度
    assert _n3d(lo) == pytest.approx(LOW, rel=1e-5)
    assert "_nax" not in hi.read_text().splitlines()[0]         # 轴口径只用于画图，不进表


def test_s83_old_2d_grid_axis_top_follows_grid(tmp_path):
    """旧 2D 网格（面密度 1e11–1e14 cm^-2，三维等效最大 ~1.5e21）：轴上端放到 1e22，点不被截掉。"""
    S = _load("s83_0077b", KE / "step8.3_output" / "gen_step13_output.py")
    cwd = tmp_path / "m"
    (cwd / "step8.4_amset2d").mkdir(parents=True)
    (cwd / "step8.4_amset2d" / "2d_correction.json").write_text(json.dumps(
        {"cell_c_A": C_A, "elastic_rescale_factor_c_over_t": C_A / T_A}))
    (cwd / "step8.4_amset2d" / "transport.json").write_text(json.dumps(
        {"doping": [v / (C_A * 1e-8) for v in _pm([1e11, 1e14])] + [0.0], "temperatures": [300.0]}))
    fac = S._n_scale(cwd, True)[0]
    assert S._doping_top(cwd, True, fac) == 1e22
    assert S._phonon_limited_ranges(S._doping_axis(cwd, True, 1e22, fac)[0], 1e22)[0][0] == "carrier_tuning_1e16-1e22"
