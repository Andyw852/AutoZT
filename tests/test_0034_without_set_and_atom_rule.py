"""[0034] S8.4 的 OVERALL_WITHOUT（overall 去掉某些机制，出厂去 ODP -> overall_noODP.json）与 S8.3 的接入；
comparison_300K 的可解读判据改每原子（与 S8.4、本征表同一判据）。纯本地，复用 0032 的合成 MoS₂ 材料。"""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
S84 = KE / "step8.4_amset2d"


def _load(name, path):
    for d in (ROOT, KE, S84, ROOT / "skill" / "_common" / "opt", ROOT / "tests"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T32 = _load("t0032_for_0034", ROOT / "tests" / "test_0032_intrinsic_suite.py")


def _s83():
    return _load("s83_0034", KE / "step8.3_output" / "gen_step13_output.py")


def _g14():
    return _load("s84_0034", S84 / "gen_step14_amset2d.py")


def _add_without(cwd, mu300=150.0, stale=False, dropped=("ODP",), name="overall_noODP.json"):
    """在 0032 材料上补一份 S8.4 的 overall_no*.json（与 transport.json 同一轮：同网格、同指纹）。"""
    run = cwd / "step8.4_amset2d"
    tj = (run / "transport.json").read_text()
    j = json.loads(tj)
    S, sg, ke, mu = T32._amset_grid(j["doping"], j["temperatures"], mu300)
    q = {"doping_cm3": j["doping"], "temperatures_K": j["temperatures"],
         "intrinsic_labels": ["ADP", "POP", "IMP", "IVS"], "intrinsic_dropped": list(dropped),
         "intrinsic": {"seebeck_uV_K": S, "conductivity_S_m": sg, "electronic_thermal_conductivity_W_mK": ke,
                       "mobility_cm2_Vs_s": {"overall": mu}},
         "transport_fingerprint": {"file": "transport_47x47x3.json",
                                   "sha256": hashlib.sha256((tj + ("x" if stale else "")).encode()).hexdigest()}}
    (run / name).write_text(json.dumps(q))


# ---------------------------------------------------------------- S8.4 gen
def test_overall_without_set_logic():
    g = _g14()
    f = g.overall_without_set
    full = ["ADP", "IMP", "POP", "ODP", "IVS"]
    assert f(None, full) == (["ODP"], "overall_noODP.json")                  # 出厂：去 ODP
    assert f("ODP,IVS", full) == (["ODP", "IVS"], "overall_noODP_IVS.json")
    assert f("odp ivs odp", full) == (["ODP", "IVS"], "overall_noODP_IVS.json")   # 大小写、空格、去重
    assert f("off", full)[0] is None and f("", full)[0] is None
    none, why = f(None, ["ADP", "IMP", "POP"])                                # 没开 ODP：就是 overall，不出
    assert none is None and "就是 overall" in why
    assert f("ADP,IMP,POP", ["ADP", "IMP", "POP"])[0] is None                 # 全剔除：不出
    with pytest.raises(SystemExit):
        f("ODQ", full)


def test_gen_chain_rm_spec_and_markers():
    src = (S84 / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert '"OVERALL_WITHOUT": (None, "str")' in src and '_p["OVERALL_WITHOUT"]' in src
    assert "\nOVERALL_WITHOUT = \"ODP\"\n" in src
    i = src.index('--drop IMP,POP,PIE,ODP,IVS"')
    j = src.index('_ow, _ow_name = overall_without_set(OVERALL_WITHOUT, SCATTERING)')
    k = src.index('" && python postprocess_intrinsic.py --drop %s --out %s" % (",".join(_ow), _ow_name)')
    assert i < j < k < src.index("_acmd = _acmd + AMSET_TAIL")              # 在落成 transport.json 之前
    assert "rm -f transport.json intrinsic_transport.json intrinsic_ADP.json overall_no*.json mesh_*.h5;" in src
    kc = _load("kc_0034", KE / "ke_common.py")
    assert {"overall_noODP.json", "overall_noODP_IVS.json"} <= set(kc.DONE_MARKERS["step8.4_amset2d"])


# ---------------------------------------------------------------- S8.3 本征套件 / μ(T)
def test_s83_suite_picks_up_overall_without(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path)
    _add_without(cwd)
    data, geom, rows = T32._rows(S, cwd)
    assert [d[0] for d in data["set_defs"]] == ["ADP_only", "intrinsic", "overall", "overall_noODP"]
    r = T32._pick(rows, 300.0, "overall_noODP", "n", T32.N2D[0])
    assert r["mu_cm2/Vs"] == pytest.approx(150.0, rel=1e-6)
    assert S._intr_legend("overall_noODP", data, "x") == "overall w/o ODP (ADP+POP+IMP+IVS)"
    L = S._intr_summary(tmp_path, rows, geom, data, ["a.png"])
    assert any("mu_noODP" in ln for ln in L) and any(ln.startswith("# overall_noODP overall_noODP.json") for ln in L)
    assert "w/o ODP: same mesh with ODP removed" in S._intr_caption(geom, data)


def test_s83_without_set_absent_or_stale(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path)
    data, *_ = T32._rows(S, cwd)
    assert [d[0] for d in data["set_defs"]] == ["ADP_only", "intrinsic", "overall"]      # 没有就不列
    _add_without(cwd, stale=True)
    data, *_ = T32._rows(S, cwd)
    assert "overall_noODP" not in data["sets"] and any("上一轮" in n for n in data["notes"])


def test_s83_figures_and_table_with_without_set(tmp_path):
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = T32._material(tmp_path)
    _add_without(cwd)
    out = cwd / "step8.3_output"
    out.mkdir()
    S.write_intrinsic_suite(out, cwd, True)
    lines = (out / "intrinsic" / "intrinsic_transport_table.csv").read_text(encoding="utf-8").splitlines()
    body = list(csv.DictReader([ln for ln in lines if not ln.startswith("#")]))
    assert {r["dataset"] for r in body} == {"ADP_only", "intrinsic", "overall", "overall_noODP", "BT2_DPTtau", "DPT"}


def test_mu_vs_T_has_without_column(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path)
    _add_without(cwd)
    am = S._amset_mu_T_sources(cwd)
    assert "mu_AMSET_all_without" in am["cols"] and am["without"] == ["ODP"]
    assert "overall_noODP.json（overall 去 ODP" in am["src"]["mu_AMSET_all_without"]
    assert S.MU_T_COLS[-1] == "mu_AMSET_all_without"


# ---------------------------------------------------------------- comparison_300K：每原子判据
def test_atom_rule_sources(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path)
    assert S._atom_rule(cwd)[:2] == (3, 0.02)                                   # 2D：2d_correction.json
    d3 = tmp_path / "Si" / "ke-dft-cpu" / "step3_uniform"
    d3.mkdir(parents=True)
    (d3 / "POSCAR").write_text("Si\n1.0\n0 2.7 2.7\n2.7 0 2.7\n2.7 2.7 0\nSi\n2\nDirect\n0 0 0\n.25 .25 .25\n")
    assert S._atom_rule(d3.parent)[:2] == (2, S.CARRIER_PER_ATOM_MAX)          # 3D：数结构文件（VASP5）
    (d3 / "POSCAR").write_text("Si\n1.0\n0 2.7 2.7\n2.7 0 2.7\n2.7 2.7 0\n2\nDirect\n0 0 0\n.25 .25 .25\n")
    assert S._atom_rule(d3.parent)[0] == 2                                      # VASP4
    assert S._atom_rule(tmp_path / "nothing")[0] is None


def test_geom_2d_counts_atoms_when_record_lacks_them(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path)
    rp = cwd / "step8.4_amset2d" / "2d_correction.json"
    rec = json.loads(rp.read_text())
    rec.pop("num_atoms")
    rec.pop("carrier_per_atom_threshold")
    rp.write_text(json.dumps(rec))
    g = S._geom_2d(cwd)
    assert g["num_atoms"] == 3 and g["thr_atom"] == 0.02                        # 不再退回每原胞 0.1


def test_comparison_300K_csv_is_per_atom(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = T32._material(tmp_path)
    rp = cwd / "step8.4_amset2d" / "2d_correction.json"
    rec = json.loads(rp.read_text())
    rec["carrier_per_atom_threshold"] = 0.005           # 阈值必须取自 S8.4 记录（同一判据）
    rp.write_text(json.dumps(rec))
    monkeypatch.chdir(cwd)
    S.main()
    out = cwd / "step8.3_output"
    rows = [r for r in csv.DictReader((out / "comparison_300K.csv").read_text(encoding="utf-8").splitlines())
            if r.get("carriers_per_cell")]
    assert rows and "carriers_per_atom" in rows[0]
    for r in rows:
        pc, pa = float(r["carriers_per_cell"]), float(r["carriers_per_atom"])
        assert pa == pytest.approx(pc / 3, rel=1e-2)
        assert r["in_interpret_range"] == ("yes" if pa <= 0.005 else "NO(>0.005/atom：超出刚性能带近似/实验范围)")
    # 至少有一档：旧的每原胞判据（≤ 0.1）会判 yes，每原子判 NO
    assert any(float(r["carriers_per_cell"]) <= 0.1 and r["in_interpret_range"].startswith("NO") for r in rows)
    assert any(r["in_interpret_range"] == "yes" for r in rows)
    txt = (out / "comparison_summary.txt").read_text(encoding="utf-8")
    assert "可解读上限 = 每原子 0.005 个载流子（3 原子/胞" in txt
