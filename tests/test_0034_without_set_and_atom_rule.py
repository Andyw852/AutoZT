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
    assert '"SUBSETS": (None, "str")' in src and '_p["SUBSETS"]' in src            # [0038]
    assert "\nOVERALL_WITHOUT = \"ODP\"\n" in src
    i = src.index('_specs = subset_specs(SUBSETS, SCATTERING)')     # [0038] 子集循环
    j = src.index('_ow, _ow_name = overall_without_set(OVERALL_WITHOUT, SCATTERING)')
    k = src.index('" && python postprocess_intrinsic.py --drop %s --out %s" % (",".join(_ow), _ow_name)')
    assert i < j < k < src.index("_acmd = _acmd + AMSET_TAIL")              # 在落成 transport.json 之前
    assert "rm -f transport.json subset_*.json intrinsic_*.json overall_no*.json tau_energy.json mesh_*.h5;" in src
    assert "'rm -f subset_*.json intrinsic_*.json overall_no*.json tau_energy.json')" in src
    assert "POSTPROCESS_ONLY requires WRITE_MESH = true" in src
    assert "global INTRINSIC_ADP_ONLY, OVERALL_WITHOUT, TAU_ENERGY, POSTPROCESS_ONLY" in src
    kc = _load("kc_0034", KE / "ke_common.py")
    assert {"overall_noODP.json", "overall_noODP_IVS.json"} <= set(kc.DONE_MARKERS["step8.4_amset2d"])


# ---------------------------------------------------------------- S8.3 本征套件 / μ(T)
def test_s83_suite_picks_up_overall_without(tmp_path):
    S = _s83()
    assert S._intr_legend("ADP_POP", {"sets": {"ADP_POP": {"labels": ["ADP"]}}}, "ADP+POP") == "amset_ADP"
    cwd = T32._material(tmp_path)
    _add_without(cwd)
    data, geom, rows = T32._rows(S, cwd)
    assert [d[0] for d in data["set_defs"]] == ["ADP_only", "ADP_POP", "ADP_POP_ODP", "intrinsic", "overall", "overall_noODP"]
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
    assert [d[0] for d in data["set_defs"]] == ["ADP_only", "ADP_POP", "ADP_POP_ODP", "intrinsic", "overall"]      # extra(overall_no*) 没有就不列
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


# ---------------------------------------------------------------- [0034-r2] 复审修正
def _dpo():
    return _load("dpo_0034", S84 / "dp_optical.py")


def test_merge_coincident_stars_and_shadowed_channel_rejected(tmp_path, monkeypatch):
    D = _dpo()
    K = [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]
    vs = [{"label": "K", "g": 2, "k_frac": K}, {"label": "Q", "g": 1, "k_frac": [[0.17, 0.17, 0]]},
          {"label": "K(up带27)", "g": 2, "k_frac": [[2 / 3, 2 / 3, 0], [1 / 3 + 1, 1 / 3, 0.5]]}]   # 同中心（周期像、kz 无关）
    out, alias = D.merge_coincident(vs)
    assert [v["label"] for v in out] == ["K", "Q"] and alias == {"K(up带27)": "K"}
    with pytest.raises(ValueError):                                         # 只部分重合：谷表自相矛盾
        D.merge_coincident([{"label": "K", "g": 2, "k_frac": K}, {"label": "X", "g": 1, "k_frac": [K[0]]}])
    d = tmp_path / "step8.2_dpt"
    d.mkdir()
    (d / "dpt_result.json").write_text(json.dumps({"results": [{"carrier": "electron", "multi_valley": {"valleys": vs}}]}))
    g = _g14()
    monkeypatch.setattr(g, "SCATTERING", ["ADP", "IMP", "POP"])
    monkeypatch.setattr(g, "DP_OPTICAL_E", "K(up带27)>K(up带27):4@50")
    with pytest.raises(SystemExit):                                         # 用被并掉的星 -> 报错，不再静默 0 率
        g.resolve_dp_optical(tmp_path)
    monkeypatch.setattr(g, "DP_OPTICAL_E", "K>K:4@50; K>Q:2@45")
    rec = g.resolve_dp_optical(tmp_path)
    assert [v["label"] for v in rec["electron"]["valleys"]] == ["K", "Q"]


def test_explicit_scattering_is_honoured(tmp_path, monkeypatch):
    d = tmp_path / "step8.2_dpt"
    d.mkdir()
    K = [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]
    (d / "dpt_result.json").write_text(json.dumps({"results": [{"carrier": "electron", "multi_valley": {"valleys": [
        {"label": "K", "g": 2, "k_frac": K}]}}]}))
    g = _g14()
    monkeypatch.setattr(g, "_SCATTERING_EXPLICIT", True)
    monkeypatch.setattr(g, "SCATTERING", ["ADP"])                           # 严格对照：只开 ADP
    monkeypatch.setattr(g, "DP_OPTICAL_E", "K>K:4@50; K>K':2.6@46")
    assert g.resolve_dp_optical(tmp_path) is None and g.SCATTERING == ["ADP"]   # 不补 ODP/IVS、不掩码
    monkeypatch.setattr(g, "SCATTERING", ["ADP", "IVS"])
    rec = g.resolve_dp_optical(tmp_path)
    assert g.SCATTERING == ["ADP", "IVS"] and [c["kind"] for c in rec["electron"]["channels"]] == ["inter_star"]
    assert rec["adp_intravalley"]["electron"] is True
    monkeypatch.setattr(g, "_SCATTERING_EXPLICIT", False)                   # 没显式写：照旧自动补
    monkeypatch.setattr(g, "SCATTERING", ["ADP", "IMP", "POP"])
    g.resolve_dp_optical(tmp_path)
    assert g.SCATTERING == ["ADP", "IMP", "POP", "ODP", "IVS"]


def test_pie_without_piezo_dropped_from_settings_too(tmp_path, monkeypatch):
    g = _g14()
    monkeypatch.setattr(g, "SCATTERING", ["ADP", "IMP", "PIE"])
    monkeypatch.chdir(tmp_path)                                             # 没有 step6_elastic/OUTCAR
    assert g._amset_scatterers() == ["ADP", "IMP"]
    assert g._mechanisms_2d(tmp_path) == ["ADP", "IMP"]                      # 两边一致


def test_postprocess_note_follows_drop_list():
    src = (S84 / "postprocess_intrinsic.py").read_text(encoding="utf-8")
    assert '"%s 行已置零后重积分；overall=sum(%s)' in src and "IMP 行已置零后重积分；overall=sum(ADP,POP)" not in src


def test_plugin_symmetrizes_dpo_rates_and_warns_empty_star():
    src = (S84 / "amset2d_plugin.py").read_text(encoding="utf-8")
    i = src.index("def _dpo_rates(")
    body = src[i:src.index("\nclass _DPOptical2D", i)]
    assert 'r = r[..., np.asarray(ir)][..., np.asarray(i2f)]' in body           # 按 AMSET 不可约映射对称化
    assert body.index("r += dpo.RATE_FLOOR") < body.index("np.asarray(ir)") < body.index("np.broadcast_to")
    assert 'r["dk_rel"] is None and r["label"] in used' in src                 # 有通道却分不到 k 点 -> WARN


def test_mu_vs_T_skips_stale_first_without_file(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path)
    _add_without(cwd, stale=True)                                           # overall_noODP.json：上一轮残留
    _add_without(cwd, dropped=("ODP", "IVS"), name="overall_noODP_IVS.json")
    am = S._amset_mu_T_sources(cwd)
    assert am["without"] == ["ODP", "IVS"] and "overall_noODP_IVS.json" in am["src"]["mu_AMSET_all_without"]


def _s84_with_dpo(tmp_path, valleys, src="step8.2_dpt/dpt_result.json"):
    run = tmp_path / "step8.4_amset2d"
    run.mkdir(parents=True)
    (run / "2d_correction.json").write_text(json.dumps({"dp_optical": {"electron": {
        "valleys": valleys, "valleys_source": src, "channels": []}}}))
    (run / "transport.json").write_text("{}")
    (run / "intrinsic_transport.json").write_text("{}")
    s83 = tmp_path / "step8.3_output"
    s83.mkdir()
    (s83 / "comparison_300K.png").write_text("x")
    return run


def test_s82_valley_change_requeues_s84_only_when_used_and_changed(tmp_path, capsys):
    kc = _load("kc_0034b", KE / "ke_common.py")
    K = [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]
    same = {"results": [{"carrier": "electron", "multi_valley": {"valleys": [
        {"label": "K", "k_frac": [[0.333333, 0.333333, 0.0], [-0.333333, -0.333333, 0.0]]},
        {"label": "K(up带27)", "k_frac": [[2 / 3, 2 / 3, 0], [1 / 3, 1 / 3, 0]]}]}}]}     # 同中心重复星：不算变
    moved = {"results": [{"carrier": "electron", "multi_valley": {"valleys": [
        {"label": "K", "k_frac": K}, {"label": "Q", "k_frac": [[0.17, 0.17, 0]]}]}}]}
    run = _s84_with_dpo(tmp_path / "a", [{"label": "K", "g": 2, "k_frac": K}])
    assert kc.invalidate_dp_optical_valleys(tmp_path / "a", same) == []
    assert (run / "transport.json").is_file()
    done = kc.invalidate_dp_optical_valleys(tmp_path / "a", moved)
    assert ("step8.4_amset2d", "transport.json") in done and not (run / "transport.json").exists()
    assert any(f.name.startswith("transport.json.stale-upstream-step8.2_dpt-") for f in run.iterdir())
    assert ("step8.3_output", "comparison_300K.png") in done                  # 递归让 S8.3 失效
    run_b = _s84_with_dpo(tmp_path / "b", [{"label": "K", "g": 2, "k_frac": K}], src="step.conf DP_OPTICAL_VALLEYS_E")
    assert kc.invalidate_dp_optical_valleys(tmp_path / "b", moved) == []     # 手填谷表：与 S8.2 无关
    assert (run_b / "transport.json").is_file()
    assert kc.invalidate_dp_optical_valleys(tmp_path / "none", moved) == []  # 没有 S8.4


def test_s82_gen_calls_valley_invalidation_after_writing_result():
    src = (KE / "step8.2_dpt" / "gen_step12_dpt.py").read_text(encoding="utf-8")
    i = src.index('(out / "dpt_result.json").write_text(')
    assert src.index("invalidate_dp_optical_consumer(cwd, res)", i) - i < 200


# ---------------------------------------------------------------- [0034-r3] 第二轮复审修正
Q6 = [[0.17, 0.17, 0], [-0.17, -0.17, 0], [0.17, -0.34, 0], [-0.17, 0.34, 0], [0.34, -0.17, 0], [-0.34, 0.17, 0]]
Q6G = [[1 / 6, 1 / 6, 0], [-1 / 6, -1 / 6, 0], [1 / 6, -1 / 3, 0], [-1 / 6, 1 / 3, 0], [1 / 3, -1 / 6, 0], [-1 / 3, 1 / 6, 0]]


def _shift(ks, d):
    return [[k[0] + d, k[1] - d / 2, k[2]] for k in ks]


def test_merge_uses_k_grid_for_off_symmetry_duplicates():
    """Q 谷不钉在高对称点：两条带各自拟合的中心差零点几个网格步长（这里 2e-3），但网格极值点 k_grid 相同 -> 合并。"""
    D = _dpo()
    dpt = {"results": [{"carrier": "electron", "multi_valley": {"valleys": [
        {"label": "Q", "k_frac": Q6, "k_grid": Q6G},
        {"label": "Q(up带27)", "k_frac": _shift(Q6, 2e-3), "k_grid": Q6G}]}}]}
    vs = D.valleys_from_dpt(dpt, "electron")
    assert all("k_grid" in v for v in vs)
    out, alias = D.merge_coincident(vs)
    assert [v["label"] for v in out] == ["Q"] and alias == {"Q(up带27)": "Q"}
    # 没有 k_grid（手填 / 老结果）：照旧按拟合中心 1e-4 判，2e-3 的差不合并
    nog = [{k: v[k] for k in v if k != "k_grid"} for v in vs]
    assert len(D.merge_coincident(nog)[0]) == 2
    with pytest.raises(ValueError):                                     # 网格点只部分重合：谷表自相矛盾
        D.merge_coincident([vs[0], dict(vs[1], k_grid=Q6G[:3] + [[0.4, 0.1, 0]] * 3)])


def test_valley_signature_uses_k_grid(tmp_path, capsys):
    kc = _load("kc_0034c", KE / "ke_common.py")
    D = _dpo()
    raw = {"results": [{"carrier": "electron", "multi_valley": {"valleys": [
        {"label": "Q", "k_frac": Q6, "k_grid": Q6G},
        {"label": "Q(up带27)", "k_frac": _shift(Q6, 3e-5), "k_grid": Q6G}]}}]}
    merged, _ = D.merge_coincident(D.valleys_from_dpt(raw, "electron"))   # gen 记下的是合并后的
    run = _s84_with_dpo(tmp_path / "a", merged)
    assert kc.invalidate_dp_optical_valleys(tmp_path / "a", raw) == []      # 同一张原始谷表：不重排
    jit = json.loads(json.dumps(raw))
    jit["results"][0]["multi_valley"]["valleys"][0]["k_frac"] = _shift(Q6, 6e-5)   # 拟合中心抖动、网格点不变
    assert kc.invalidate_dp_optical_valleys(tmp_path / "a", jit) == []
    assert (run / "transport.json").is_file()
    moved = json.loads(json.dumps(raw))
    for v in moved["results"][0]["multi_valley"]["valleys"]:
        v["k_grid"] = [[k[0] + 1 / 12, k[1], 0] for k in Q6G]             # 网格极值点变了
    assert ("step8.4_amset2d", "transport.json") in kc.invalidate_dp_optical_valleys(tmp_path / "a", moved)


def test_valley_change_with_job_in_flight_warns_not_archives(tmp_path, capsys):
    kc = _load("kc_0034d", KE / "ke_common.py")
    K = [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]
    run = _s84_with_dpo(tmp_path, [{"label": "K", "g": 2, "k_frac": K}])
    (run / "slurm-123.out").write_text("running")                            # 刚更新过：作业在跑
    moved = {"results": [{"carrier": "electron", "multi_valley": {"valleys": [{"label": "Q", "k_frac": [[0.17, 0.17, 0]]}]}}]}
    assert kc.invalidate_dp_optical_valleys(tmp_path, moved) == []
    out = capsys.readouterr().out
    assert "★" in out and "scancel" in out and (run / "transport.json").is_file()
    run2 = _s84_with_dpo(tmp_path / "b", [{"label": "K", "g": 2, "k_frac": K}])
    for m in ("transport.json", "intrinsic_transport.json"):
        (run2 / m).unlink()
    (run2 / "submit.sh").write_text("#!/bin/bash\n")                         # 生成过、还没结果
    assert kc.invalidate_dp_optical_valleys(tmp_path / "b", moved) == []
    assert "还没有结果" in capsys.readouterr().out


def test_explicit_scattering_with_mask_on_only_masks(tmp_path, monkeypatch, capsys):
    d = tmp_path / "step8.2_dpt"
    d.mkdir()
    K = [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]
    (d / "dpt_result.json").write_text(json.dumps({"results": [{"carrier": "electron", "multi_valley": {"valleys": [
        {"label": "K", "g": 2, "k_frac": K}]}}]}))
    g = _g14()
    monkeypatch.setattr(g, "_SCATTERING_EXPLICIT", True)
    monkeypatch.setattr(g, "SCATTERING", ["ADP"])
    monkeypatch.setattr(g, "DP_OPTICAL_E", "K>K:4@50; K>K':2.6@46")
    monkeypatch.setattr(g, "ADP_INTRAVALLEY", "on")
    rec = g.resolve_dp_optical(tmp_path)
    assert g.SCATTERING == ["ADP"] and rec["electron"]["channels"] == [] and rec["adp_intravalley"]["electron"]
    out = capsys.readouterr().out
    assert "不加 ODP/IVS" in out and "ODP/IVS 开：无" not in out
