"""[0041] 缺陷平衡载流子浓度 -> AMSET 掺杂点 -> S8.3 单独出图；odp_frozen 默认开。
1. ke_common.defect_doping_points：3D cm^-3（净 = p0 − n0）、2D 换面密度（× 缺陷超胞 V/A）、
   垂直能量 auto 不用 / on 强制、off、没有结果；
2. S8 / S8.4 gen：_add_defect_doping 追加掺杂、去重；
3. S8.3 defect_rows / write_defect_plots：按 defect_doping.json 的 amset_doping_cm3 取点出图，旧图清理；
4. skill.yaml：odp_frozen 默认 true。"""
import importlib.util
import json
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _load(name, path):
    for d in (ROOT, KE, KE / "step8.4_amset2d", ROOT / "skill" / "_common" / "opt", ROOT / "tests",
              KE / "step2_bandgap" / "step2.15_discriminant"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


KC = _load("ke_common_0041", KE / "ke_common.py")
T32 = _load("t0032_for_0041", ROOT / "tests" / "test_0032_intrinsic_suite.py")


def _defect(mat, thermo=True, verts=((-0.0, 1e17, 1e12, "n 型"), (1.0, 1e10, 3e16, "p 型")), vol=2000.0):
    d = mat / "defect-dft-cpu"
    (d / "step4_analysis").mkdir(parents=True)
    (d / "step1_bulk").mkdir()
    (d / "step1_bulk" / "POSCAR").write_text("X\n1.0\n10 0 0\n0 10 0\n0 0 20\nMo\n1\nDirect\n0 0 0\n")
    (d / "energies.json").write_text(json.dumps({"T": 500.0}))
    res = {"status": "done", "is_metal": False, "volume_AA3": vol, "thermodynamic_charge_states": thermo,
           "energy_interpretation": "thermodynamic" if thermo else "vertical",
           "vertex_summary": [{"vertex": i, "E_F_eq": ef, "n0_cm3": n0, "p0_cm3": p0, "carrier_type": ct}
                              for i, (ef, n0, p0, ct) in enumerate(verts)]}
    (d / "step4_analysis" / "formation_energy_results.json").write_text(json.dumps(res))
    cwd = mat / "ke-dft-cpu"
    cwd.mkdir(parents=True, exist_ok=True)
    return cwd


def test_points_3d(tmp_path):
    cwd = _defect(tmp_path / "GaN")
    pts, notes = KC.defect_doping_points(cwd, False)
    assert [p["unit"] for p in pts] == ["cm^-3", "cm^-3"]
    assert pts[0]["net"] == pytest.approx(1e12 - 1e17)            # n 型：负
    assert pts[1]["net"] == pytest.approx(3e16 - 1e10)            # p 型：正
    assert pts[0]["T_K"] == 500.0
    assert any("冻结" in n for n in notes)


def test_points_2d_areal(tmp_path):
    cwd = _defect(tmp_path / "MoS2")
    pts, _ = KC.defect_doping_points(cwd, True)
    # V/A = 2000 Å^3 / 100 Å^2 = 20 Å = 2e-7 cm
    assert pts[0]["unit"] == "cm^-2"
    assert pts[0]["net"] == pytest.approx((1e12 - 1e17) * 2e-7)


def test_vertical_auto_skips_on_forces_off_and_missing(tmp_path):
    cwd = _defect(tmp_path / "A", thermo=False)
    assert KC.defect_doping_points(cwd, False)[0] == []
    assert len(KC.defect_doping_points(cwd, False, "on")[0]) == 2
    assert KC.defect_doping_points(cwd, False, "off")[0] == []
    lone = tmp_path / "B" / "ke-dft-cpu"
    lone.mkdir(parents=True)
    pts, notes = KC.defect_doping_points(lone, False)
    assert pts == [] and "defect-dft-cpu" in notes[0]


@pytest.mark.parametrize("gen,is2d,scale", [
    (KE / "step8_amset" / "gen_step10_amset.py", False, 1.0),
    (KE / "step8.4_amset2d" / "gen_step14_amset2d.py", True, 1.0 / (20e-8)),
])
def test_gen_appends_doping(tmp_path, gen, is2d, scale):
    g = _load("gen_0041_%d" % is2d, gen)
    cwd = _defect(tmp_path / "M")
    dop, rec = g._add_defect_doping([-1e18, 1e18], cwd, is2d, scale)
    assert len(dop) == 4 and len(rec["points"]) == 2
    assert rec["points"][0]["amset_doping_cm3"] == pytest.approx(rec["points"][0]["net"] * scale)
    dop2, _ = g._add_defect_doping(dop, cwd, is2d, scale)             # 已有的不重复加
    assert len(dop2) == 4
    src = gen.read_text(encoding="utf-8")
    assert 'DEFECT_DOPING = "auto"' in src and '"DEFECT_DOPING": ("auto", "str")' in src


def test_s83_defect_plot(tmp_path):
    pytest.importorskip("matplotlib")
    S = _load("s83_0041", KE / "step8.3_output" / "gen_step13_output.py")
    cwd = T32._material(tmp_path, doping0=True)
    tj = json.loads((cwd / "step8.4_amset2d" / "transport.json").read_text())
    v = tj["doping"][1]                                               # 拿网格里已有的一档当缺陷点
    (cwd / "step8.4_amset2d" / "defect_doping.json").write_text(json.dumps({"points": [
        {"vertex": 0, "carrier_type": "n 型", "net": v * T32.C_A * 1e-8, "unit": "cm^-2",
         "E_F_eq_eV": 1.2, "T_K": 500.0, "amset_doping_cm3": v}]}))
    rows, why = S.defect_rows(cwd, True)
    assert why is None and {r["dataset"] for r in rows} >= {"overall", "intrinsic"}
    plots = cwd.parent / "transport_plots"
    assert S.write_defect_plots(cwd, plots, "MoS2", True)
    d = plots / S.DEFECT_DIR_OUT
    assert (d / "MoS2_defect_equilibrium.png").is_file()
    assert "vertex" in (d / "MoS2_defect_equilibrium.csv").read_text().splitlines()[1]
    (cwd / "step8.4_amset2d" / "defect_doping.json").unlink()            # 没记录 -> 跳过并清旧图
    assert not S.write_defect_plots(cwd, plots, "MoS2", True)
    assert not (d / "MoS2_defect_equilibrium.png").exists()


def test_odp_frozen_default_on():
    sk = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    assert sk["optional_steps"]["odp_frozen"]["default"] is True
