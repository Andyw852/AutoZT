"""[0032] 本征五参数图 + 输运参数表（S8.3）、S8.4 默认开 WRITE_MESH / INTRINSIC_ADP_ONLY、
autozt 整目录回拉的 fetch_exclude。纯本地，合成 MoS2 型 2D 材料。"""
import csv
import importlib.util
import json
import math
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"

E_C, M0 = 1.602176634e-19, 9.1093837015e-31
C_A, T_A = 34.88, 6.72                     # 胞高（含真空）与 vdW 层厚
N2D = [1e11, 5.62e11, 3.16e12, 1.78e13, 1e14]
TEMPS = [300.0, 600.0]
MU_E, MU_H, MD = 200.0, 1600.0, 0.5         # DPT 单谷 μ(300 K) 与 m_d（BT2 用 m*=MD，μ_BT2 应 = μ_DPT）
AREA_CM2 = math.sqrt(3) / 2 * 3.18 ** 2 * 1e-16


def _load(name, path):
    for d in (ROOT, KE, ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _s83():
    return _load("s83_0032", KE / "step8.3_output" / "gen_step13_output.py")


def _t(v, z=1e-9):
    return [[v, 0, 0], [0, v, 0], [0, 0, z]]


def _amset_grid(dop, temps, mu300):
    """σ = n e μ（元胞口径），μ ∝ 1/T；S 按非简并单带粗估；κe = L σ T（L=1.6e-8）。"""
    S, sig, ke, mu = [], [], [], []
    for n3 in dop:
        rs, rg, rk, rm = [], [], [], []
        for T in temps:
            m = mu300 * 300.0 / T
            s = abs(n3) * 1e6 * E_C * m * 1e-4
            seeb = ((-1 if n3 < 0 else 1) * 86.17 * (2.0 + math.log(1e19 / abs(n3)))) if n3 else 0.0
            rs.append(_t(seeb)), rg.append(_t(s)), rk.append(_t(1.6e-8 * s * T)), rm.append(_t(m))
        S.append(rs), sig.append(rg), ke.append(rk), mu.append(rm)
    return S, sig, ke, mu


def _material(tmp_path, adp_file=True, intr_file=True, bad_grid=False, stale=False, doping0=False):
    mat = tmp_path / "MoS2"
    cwd = mat / "ke-dft-cpu"
    for d in ("step1_opt", "step3_uniform", "step8.4_amset2d", "step8.1_boltztrap", "step8.2_dpt"):
        (cwd / d).mkdir(parents=True)
    (cwd / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\nFUNC=pbe\n")
    a = 3.18
    (cwd / "step3_uniform" / "POSCAR").write_text(
        "MoS2\n1.0\n%.6f 0 0\n%.6f %.6f 0\n0 0 %.4f\nMo S\n1 2\nDirect\n0 0 0.5\n"
        "0.3333 0.6667 0.545\n0.3333 0.6667 0.455\n" % (a, -a / 2, a * math.sqrt(3) / 2, C_A))
    c_cm = C_A * 1e-8
    dop = [-n / c_cm for n in N2D] + ([0.0] if doping0 else []) + [n / c_cm for n in N2D]
    run = cwd / "step8.4_amset2d"
    (run / "2d_correction.json").write_text(json.dumps({
        "cell_c_A": C_A, "elastic_rescale_factor_c_over_t": round(C_A / T_A, 4),
        "layer_thickness": {"method": "vdw"}, "areal_density_factor_cm": c_cm,
        "inplane_cell_area_cm2": AREA_CM2, "num_atoms": 3, "carrier_per_atom_threshold": 0.02}))
    S, sg, ke, mu = _amset_grid(dop, TEMPS, 100.0)
    _, _, _, mu_adp = _amset_grid(dop, TEMPS, 400.0)
    tj = json.dumps({"doping": dop, "temperatures": TEMPS, "seebeck": S, "conductivity": sg,
                     "electronic_thermal_conductivity": ke,
                     "mobility": {"overall": mu, "ADP": mu_adp, "IMP": mu, "POP": mu}})
    (run / "transport.json").write_text(tj)
    import hashlib
    sha = hashlib.sha256((tj + ("x" if stale else "")).encode()).hexdigest()   # stale：记录的是上一轮内容
    for name, mu300, kept, ok in (("intrinsic_ADP.json", 400.0, ["ADP"], adp_file),
                                  ("intrinsic_transport.json", 300.0, ["ADP", "POP"], intr_file)):
        if not ok:
            continue
        S_, sg_, ke_, mu_ = _amset_grid(dop, TEMPS, mu300)
        q = {"doping_cm3": dop if not bad_grid else dop[::-1], "temperatures_K": TEMPS,
             "intrinsic_labels": kept, "intrinsic_dropped": ["IMP"],
             "intrinsic": {"seebeck_uV_K": S_, "conductivity_S_m": sg_,
                           "electronic_thermal_conductivity_W_mK": ke_,
                           "mobility_cm2_Vs_s": {"overall": mu_}}}
        if kept == ["ADP", "POP"]:
            q["reproduction_check"] = {"passed": True, "max_abs_diff": {"overall": 1e-10}}
            q["doping_marks"] = {str(len(dop) - 1): ["超出刚带近似"]}
        q["transport_fingerprint"] = {"file": "transport_47x47x3.json", "sha256": sha}
        (run / name).write_text(json.dumps(q))
    # BT2：净浓度 ±1e15..1e22（元胞口径），σ/τ = |n| e² / (MD m0)，S、L 粗估
    ns = [-(10 ** (22 - 7 * i / 199)) for i in range(200)] + [10 ** (15 + 7 * i / 199) for i in range(200)]
    def row(T):
        sot = [abs(n) * 1e6 * E_C ** 2 / (MD * M0) for n in ns]
        return {"n": ns, "S": [(-1 if n < 0 else 1) * 86.17e-6 * (2 + math.log(1e19 / abs(n))) for n in ns],
                "L": [1.5e-8] * len(ns), "sot": sot, "kot": [1.5e-8 * v * T for v in sot]}
    rows = [row(T) for T in TEMPS]
    (cwd / "step8.1_boltztrap" / "boltztrap_crta.json").write_text(json.dumps({
        "temperatures_K": TEMPS, "carrier_conc_cm-3": [r["n"] for r in rows],
        "seebeck_V_per_K": [r["S"] for r in rows], "lorenz_WOhm_per_K2": [r["L"] for r in rows],
        "sigma_over_tau_S_per_m_s": [r["sot"] for r in rows],
        "kappa_e_over_tau_W_per_m_K_s": [r["kot"] for r in rows]}))
    res = []
    for car, m in (("electron", MU_E), ("hole", MU_H)):
        res.append({"carrier": car, "mobility_cm2_Vs": m, "inputs": {"m_eff_m0": MD},
                    "by_direction": {"status": "ok", "m_d_m0": MD,
                                     "x": {"mobility_cm2_Vs": m, "tau_s": m * 1e-4 * MD * M0 / E_C},
                                     "y": {"mobility_cm2_Vs": m, "tau_s": m * 1e-4 * MD * M0 / E_C}}})
    (cwd / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(
        {"temperature_K": 300.0, "is_2d": True, "results": res}))
    return cwd


def _rows(S, cwd):
    data = S.load_intrinsic_sets(cwd)
    geom = S._geom_2d(cwd)
    rows, _ = S.build_intrinsic_rows(cwd, data, geom)
    return data, geom, rows


def _pick(rows, T, ds, typ, n2d=None):
    for r in rows:
        if (r["T_K"] == T and r["dataset"] == ds and r["type"] == typ
                and (n2d is None or abs(r["n_2D_cm-2"] - n2d) <= 1e-6 * n2d)):
            return r
    raise AssertionError("no row %s %s %s %s" % (T, ds, typ, n2d))


# ---------------------------------------------------------------- S8.3
def test_geometry_and_concentration_conversions(tmp_path):
    S = _s83()
    cwd = _material(tmp_path)
    data, geom, rows = _rows(S, cwd)
    assert geom["c_A"] == C_A and geom["t_A"] == pytest.approx(T_A, rel=1e-4) and geom["t_method"] == "vdw"
    r = _pick(rows, 300.0, "intrinsic", "n", 1e11)
    assert r["n_cell_cm-3"] == pytest.approx(-1e11 / (C_A * 1e-8))          # AMSET 口径 2.87e17
    assert r["n_lit_cm-3"] == pytest.approx(1e11 / (T_A * 1e-8), rel=1e-4)  # 文献口径 1.49e18
    assert r["n_lit_cm-3"] == pytest.approx(abs(r["n_cell_cm-3"]) * C_A / T_A, rel=1e-4)
    assert r["in_interpret_range"] == "yes"
    # 判据与 S8.4 同口径（每原子 0.02）：MoS2 1e14 每胞 0.088、每原子 0.029 -> NO；1.78e13 每原子 0.005 -> yes
    r14 = _pick(rows, 300.0, "intrinsic", "p", 1e14)
    assert r14["carriers_per_cell"] == pytest.approx(0.0876, rel=1e-2)
    assert r14["carriers_per_atom"] == pytest.approx(0.0292, rel=1e-2)
    assert r14["in_interpret_range"] == "NO" and r14["S84_marks"] == "超出刚带近似"
    assert _pick(rows, 300.0, "intrinsic", "n", 1.78e13)["in_interpret_range"] == "yes"


def test_large_cell_not_flagged_per_atom():
    """C60 网络型大原胞（146 Å²、60 原子）在 1.78e13：每胞 0.26、每原子 0.0043 -> 可解读（09-28 修掉的误标）。"""
    S = _s83()
    g = {"c_A": 20.0, "t_A": 10.0, "c_over_t": 2.0, "area_cm2": 146e-16, "num_atoms": 60, "thr_atom": 0.02}
    d = S._dose_cols(-1.78e13 / (20e-8), g, None)
    assert d["carriers_per_cell"] == pytest.approx(0.26, rel=1e-2)
    assert d["in_interpret_range"] == "yes"
    g2 = dict(g, num_atoms=None, thr_atom=None)                      # 缺键：退回每原胞 0.1
    assert S._dose_cols(-1.78e13 / (20e-8), g2, None)["in_interpret_range"] == "NO"


def test_values_basis_and_lorenz(tmp_path):
    S = _s83()
    cwd = _material(tmp_path)
    data, geom, rows = _rows(S, cwd)
    r = _pick(rows, 600.0, "intrinsic", "p", 3.16e12)
    n3 = 3.16e12 / (C_A * 1e-8)
    sig_cell = n3 * 1e6 * E_C * 150.0 * 1e-4                               # μ(600)=300/2
    assert r["mu_cm2/Vs"] == pytest.approx(150.0)
    assert r["sigma_S/m"] == pytest.approx(sig_cell * C_A / T_A, rel=1e-4)  # 按层厚 t 归一
    assert r["sheet_sigma_S"] == pytest.approx(sig_cell * C_A * 1e-10)      # 与口径无关
    assert r["Lorenz_WOhm/K2"] == pytest.approx(1.6e-8)                     # 口径抵消
    assert r["PF_W/mK2"] == pytest.approx((r["S_uV/K"] * 1e-6) ** 2 * r["sigma_S/m"])
    assert r["mu_cm2/Vs"] == pytest.approx(r["sigma_S/m"] / (E_C * r["n_lit_cm-3"] * 1e6) * 1e4, rel=1e-4)
    assert _pick(rows, 600.0, "ADP_only", "p", 3.16e12)["mu_cm2/Vs"] == pytest.approx(200.0)
    assert _pick(rows, 600.0, "overall", "p", 3.16e12)["mu_cm2/Vs"] == pytest.approx(50.0)


def test_bt2_line_matches_dpt_and_dpt_rows(tmp_path):
    """BT2 用 m* = m_d 时 μ_BT2 = eτ/m* = μ_DPT（同一个 τ）；DPT 每个温度一行，2D ∝ 1/T。"""
    S = _s83()
    cwd = _material(tmp_path)
    _, _, rows = _rows(S, cwd)
    for T in TEMPS:
        for typ, mu in (("n", MU_E), ("p", MU_H)):
            b = _pick(rows, T, "BT2_DPTtau", typ, 1e11)
            assert b["mu_cm2/Vs"] == pytest.approx(mu * 300.0 / T, rel=2e-3)
            assert b["Lorenz_WOhm/K2"] == pytest.approx(1.5e-8)               # BT2 自己的 L
            d = _pick(rows, T, "DPT", typ)
            assert d["mu_cm2/Vs"] == pytest.approx(mu * 300.0 / T)


def test_table_summary_and_figures(tmp_path):
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = _material(tmp_path)
    out = cwd / "step8.3_output"
    out.mkdir()
    files = S.write_intrinsic_suite(out, cwd, True)
    d = out / "intrinsic"
    assert (d / "intrinsic_transport_table.csv").is_file() and (d / "intrinsic_summary.txt").is_file()
    assert {"intrinsic_transport_300K.png", "intrinsic_transport_600K.png"} <= set(files)
    for f in ("intrinsic_transport_300K.png", "intrinsic_transport_600K.png"):
        assert (d / f).stat().st_size > 20000
    lines = (d / "intrinsic_transport_table.csv").read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("#") and "c/t" in lines[1]
    body = list(csv.DictReader([l for l in lines if not l.startswith("#")]))
    sets = {r["dataset"] for r in body}
    assert sets == {"ADP_only", "intrinsic", "overall", "BT2_DPTtau", "DPT"}
    assert len([r for r in body if r["dataset"] == "intrinsic"]) == len(N2D) * 2 * len(TEMPS)
    txt = (d / "intrinsic_summary.txt").read_text(encoding="utf-8")
    assert "复现检查 PASS" in txt and "c/t = " in txt


def test_adp_fallback_and_grid_mismatch(tmp_path):
    S = _s83()
    cwd = _material(tmp_path, adp_file=False, bad_grid=True)
    data = S.load_intrinsic_sets(cwd)
    assert data["sets"]["ADP_only"]["partial"] is True                      # 只有 transport.json 的 μ_ADP
    assert "intrinsic" not in data["sets"]                                  # 网格不同 -> 不用
    assert any("只有 transport.json 的单机制" in n for n in data["notes"])
    assert any("网格与 transport.json 不同" in n for n in data["notes"])
    geom = S._geom_2d(cwd)
    rows, _ = S.build_intrinsic_rows(cwd, data, geom)
    r = _pick(rows, 300.0, "ADP_only", "n", 1e11)
    assert r["mu_cm2/Vs"] == pytest.approx(400.0) and r.get("sigma_S/m") is None


def test_stale_intrinsic_from_previous_run_rejected(tmp_path):
    """本征文件记录的 transport 指纹与当前 transport.json 不同 -> 上一轮残留，不用（不靠 mtime）。"""
    S = _s83()
    cwd = _material(tmp_path, stale=True)
    data = S.load_intrinsic_sets(cwd)
    assert "intrinsic" not in data["sets"]
    assert data["sets"]["ADP_only"]["partial"] is True              # 退回 transport.json 的 μ_ADP
    assert sum("上一轮 S8.4 留下的" in n for n in data["notes"]) == 2
    mu = S._amset_mu_T_sources(cwd)
    assert "mu_AMSET_intrinsic" not in mu["cols"]


def test_zero_doping_skipped_and_legends(tmp_path):
    S = _s83()
    cwd = _material(tmp_path, doping0=True)
    data, geom, rows = _rows(S, cwd)
    assert all(r.get("n_2D_cm-2") != 0 for r in rows if r["dataset"] != "DPT")
    assert S._intr_legend("intrinsic", data, "x") == "amset_ADP+POP"
    assert S._intr_legend("overall", data, "x") == "amset_ADP+IMP+POP"
    data["sets"]["intrinsic"]["labels"] = ["ADP"]                      # 非极性材料：没有 POP
    assert S._intr_legend("intrinsic", data, "x") == "amset_ADP"


def test_bt2_mobility_masked_in_bipolar_region(tmp_path):
    S = _s83()
    cwd = _material(tmp_path)
    p = cwd / "step8.1_boltztrap" / "boltztrap_crta.json"
    j = json.loads(p.read_text())
    j["seebeck_V_per_K"] = [[abs(v) for v in row] for row in j["seebeck_V_per_K"]]   # n 型也给正 S：双极
    p.write_text(json.dumps(j))
    _, _, rows = _rows(S, cwd)
    b = _pick(rows, 300.0, "BT2_DPTtau", "n", 1e11)
    assert b.get("mu_cm2/Vs") is None and b["flag"].startswith("bipolar")
    assert _pick(rows, 300.0, "BT2_DPTtau", "p", 1e11)["mu_cm2/Vs"] is not None


def test_dpt_multi_valley_labelled(tmp_path, monkeypatch):
    S = _s83()
    cwd = _material(tmp_path)
    real = S.load_dpt

    def fake(c, T=None, quiet=False):
        d = real(c, T=T, quiet=quiet)
        d["src_electron"] = "multi_valley"
        return d
    monkeypatch.setattr(S, "load_dpt", fake)
    _, _, rows = _rows(S, cwd)
    assert _pick(rows, 300.0, "DPT_mv", "n")["mu_cm2/Vs"] == pytest.approx(MU_E)
    assert _pick(rows, 300.0, "DPT", "p")["mu_cm2/Vs"] == pytest.approx(MU_H)


def test_old_outputs_cleared(tmp_path):
    S = _s83()
    cwd = _material(tmp_path)
    d = cwd / "step8.3_output" / "intrinsic"
    d.mkdir(parents=True)
    old = d / "intrinsic_transport_100K.png"
    old.write_bytes(b"old")
    (d / "intrinsic_transport_table.csv").write_text("old")
    (cwd / "step8.4_amset2d" / "2d_correction.json").unlink()        # 本轮跳过：旧文件也不能留着冒充
    assert S.write_intrinsic_suite(cwd / "step8.3_output", cwd, True) == []
    assert not old.exists() and not (d / "intrinsic_transport_table.csv").exists()


def test_suite_skips_3d_and_missing_geometry(tmp_path):
    S = _s83()
    cwd = _material(tmp_path)
    out = cwd / "step8.3_output"
    out.mkdir()
    assert S.write_intrinsic_suite(out, cwd, False) == []
    (cwd / "step8.4_amset2d" / "2d_correction.json").unlink()
    assert S.write_intrinsic_suite(out, cwd, True) == []


def test_main_end_to_end(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = _material(tmp_path)
    monkeypatch.chdir(cwd)
    S.main()
    assert (cwd / "step8.3_output" / "comparison_300K.png").is_file()
    assert (cwd / "step8.3_output" / "intrinsic" / "intrinsic_transport_table.csv").is_file()


def test_loaders_default_temperature_unchanged(tmp_path):
    """refactor：load_bt2 / load_dpt 不给 T 时仍在 TARGET_T（300 K）。"""
    S = _s83()
    cwd = _material(tmp_path)
    assert S.load_bt2(cwd)["T_used"] == 300.0 and S.load_bt2(cwd, T=600.0)["T_used"] == 600.0
    assert S.load_dpt(cwd)["electron"] == pytest.approx(MU_E)
    assert S.load_dpt(cwd, T=600.0, quiet=True)["electron"] == pytest.approx(MU_E / 2)


def test_postprocess_transport_fingerprint(tmp_path):
    pp = _load("pp_0032", KE / "step8.4_amset2d" / "postprocess_intrinsic.py")
    assert pp.transport_fingerprint(tmp_path) is None
    (tmp_path / "transport_47x47x3.json").write_text('{"a": 1}')
    fp = pp.transport_fingerprint(tmp_path)
    import hashlib
    assert fp == {"file": "transport_47x47x3.json", "sha256": hashlib.sha256(b'{"a": 1}').hexdigest()}


# ---------------------------------------------------------------- S8.4 开关
def test_s84_intrinsic_switches_default_on_and_explicit_off():
    src = (KE / "step8.4_amset2d" / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert "\nWRITE_MESH = True\n" in src and "\nINTRINSIC_ADP_ONLY = True\n" in src
    assert '"WRITE_MESH": (None, "bool")' in src and '"INTRINSIC_ADP_ONLY": (None, "bool")' in src
    g = _load("s84_0032", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")
    f = g.resolve_intrinsic_switches
    assert f(None, None) == (True, True)                 # 没写 = 出厂（开）
    assert f(False, None) == (False, False)              # 显式关 mesh -> 只 ADP 跟着关
    assert f(None, False) == (True, False)
    assert f(False, True) == (True, True)                # 显式要只 ADP -> 隐含 mesh
    assert f(True, None, False, False) == (True, False)  # 旧出厂值 + 只开 mesh
    # [0038] SUBSETS 驱动子集生成（subset_<机制>.json）；INTRINSIC_ADP_ONLY=false 摘掉"只 ADP"子集
    assert '"SUBSETS": (None, "str")' in src
    assert "def subset_specs(" in src
    i = src.index("_specs = subset_specs(SUBSETS, SCATTERING)")
    assert "for _keep, _drop, _name in _specs:" in src
    assert "--drop" in src and "subset_" in src
    assert "if not INTRINSIC_ADP_ONLY:" in src
    # 链首清掉上一轮的子集/本征产物与 mesh
    assert "rm -f transport.json subset_*.json intrinsic_*.json undoped_ni.json overall_no*.json tau_energy.json mesh_*.h5;" in src


# ---------------------------------------------------------------- autozt fetch_exclude
def test_fetch_exclude_helper_and_receipt():
    sys.path.insert(0, str(ROOT))
    from autozt import workflow as w
    assert w.fetch_excludes({"fetch_exclude": ["mesh_*.h5"]}) == ["mesh_*.h5"]
    assert w.fetch_excludes({"fetch_exclude": "mesh_*.h5"}) == ["mesh_*.h5"]
    assert w.fetch_excludes({}) == []
    m = {"_seg": {"steps_cfg": [{"name": "a", "fetch_all": True},
                                {"name": "b", "fetch_all": True, "fetch_exclude": ["mesh_*.h5"]}]},
         "fetch_files": []}
    ra = w._fetch_receipt({}, m, {"name": "a", "dir": "/x"})
    rb = w._fetch_receipt({}, m, {"name": "b", "dir": "/y"})
    assert "exclude" not in ra                           # 没配的步骤回执不变（不会被判失效重拉）
    assert rb["exclude"] == ["mesh_*.h5"]
    # 提交前输入检查：排除项（上一轮 mesh）不算输入，否则 retry 后 start 被"本地仍缺少 mesh"拦下
    assert w._drop_excluded(("INCAR", "mesh_35x35x5.h5", "d/mesh_1.h5", "submit.sh"),
                            ["mesh_*.h5"]) == ("INCAR", "submit.sh")
    src = (ROOT / "autozt" / "workflow.py").read_text(encoding="utf-8")
    assert '_exs = "" if all_files else' in src            # fetch --all 不套排除项


def test_skill_yaml_excludes_mesh_from_s84_fetch():
    yaml = pytest.importorskip("yaml")
    d = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    st = d["optional_steps"]["amset2d"]["steps"][0]
    assert st["name"] == "step8.4_amset2d" and st["fetch_exclude"] == ["mesh_*.h5"]
