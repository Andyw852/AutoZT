"""[0053][0054] S8.3 transport_plots/intrinsic_carrier/<M>_mechanism_resolved_mobility.{png,csv}（原 dilute_mobility）：
稀薄极限（f ≈ 0）迁移率按散射机制拆开 —— 空穴取 n = p 本征点（AMSET 对 doping ≥ 0 只积价带），
电子取最低 n 型点；总值去 IMP；Matthiessen 倒数和；前因子为 0 的机制（μ ~ 1e47）不列。"""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _load(name, path):
    for d in (ROOT, KE, ROOT / "skill" / "_common" / "opt", ROOT / "tests"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T32 = _load("t0032_for_0053", ROOT / "tests" / "test_0032_intrinsic_suite.py")


def _s83():
    return _load("s83_0053", KE / "step8.3_output" / "gen_step13_output.py")


def _set_mobility(cwd, mech300, run="step8.4_amset2d"):
    """把 transport.json 的单机制 μ 换成 {机制: μ(300 K)}（∝ 1/T），并更新子集文件的同轮指纹。"""
    tp = cwd / run / "transport.json"
    j = json.loads(tp.read_text())
    temps = j["temperatures"]
    for k, v in mech300.items():
        j["mobility"][k] = [[T32._t(v * 300.0 / T) for T in temps] for _ in j["doping"]]
    for k in [k for k in j["mobility"] if k != "overall" and k not in mech300]:
        del j["mobility"][k]
    tj = json.dumps(j)
    tp.write_text(tj)
    sha = hashlib.sha256(tj.encode()).hexdigest()
    for f in (cwd / run).glob("intrinsic_*.json"):
        q = json.loads(f.read_text())
        q["transport_fingerprint"] = {"file": "transport_x.json", "sha256": sha}
        f.write_text(json.dumps(q))


def _pick(rows, car, role, mech, T=300.0):
    return next((r for r in rows if r["carrier"] == car and r["role"] == role
                 and r["mechanism"] == mech and r["T_K"] == T), None)


def test_dilute_rows_2d_hole_from_undoped_point(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path, doping0=True)
    _set_mobility(cwd, {"ADP": 400.0, "POP": 200.0, "IMP": 1e47})
    rows, why = S.mechanism_resolved_rows(cwd, True)
    assert why is None
    h = _pick(rows, "hole", "main", "ADP")
    assert h["source"].startswith("n = p") and h["doping_cm3"] == 0.0
    e = _pick(rows, "electron", "main", "ADP")
    assert e["doping_cm3"] < 0 and e["n_2D_cm-2"] == pytest.approx(1e11, rel=1e-6)   # 最低 n 型
    for car in ("electron", "hole"):
        assert _pick(rows, car, "main", "IMP") is None                               # IMP 不进稀薄极限
        mt = _pick(rows, car, "main", S.MATTHIESSEN)["mu_cm2/Vs"]
        assert mt == pytest.approx(1.0 / (1 / 400.0 + 1 / 200.0))
        sh = [r["share_1/mu"] for r in rows if r["carrier"] == car and r["role"] == "main"
              and r["T_K"] == 300.0 and r["share_1/mu"] is not None]
        assert sum(sh) == pytest.approx(1.0)
        assert _pick(rows, car, "main", "POP")["share_1/mu"] == pytest.approx(2 / 3)   # POP 限制
        assert _pick(rows, car, "main", S.TOTAL_NOIMP)["mu_cm2/Vs"] == pytest.approx(300.0)  # 去 IMP 子集
    assert _pick(rows, "hole", "check", S.TOTAL_NOIMP)["doping_cm3"] > 0              # 最低 p 型核对点
    assert _pick(rows, "electron", "check", S.TOTAL_NOIMP)["n_2D_cm-2"] == pytest.approx(5.62e11, rel=1e-6)


def test_zero_prefactor_mechanism_dropped(tmp_path):
    """非极性材料的 POP（前因子 0，μ ~ 1e47）不画、不进 Matthiessen。"""
    S = _s83()
    cwd = T32._material(tmp_path, doping0=True)
    _set_mobility(cwd, {"ADP": 400.0, "POP": 1e47})
    rows, _ = S.mechanism_resolved_rows(cwd, True)
    assert _pick(rows, "hole", "main", "POP") is None
    assert _pick(rows, "hole", "main", S.MATTHIESSEN)["mu_cm2/Vs"] == pytest.approx(400.0)
    assert _pick(rows, "hole", "main", "ADP")["share_1/mu"] == pytest.approx(1.0)


def test_no_undoped_point_falls_back_to_lowest_p(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path, doping0=False)
    _set_mobility(cwd, {"ADP": 400.0, "POP": 200.0})
    rows, _ = S.mechanism_resolved_rows(cwd, True)
    h = _pick(rows, "hole", "main", "ADP")
    assert h["doping_cm3"] > 0 and "lowest p-type" in h["source"]
    assert _pick(rows, "hole", "check", S.TOTAL_NOIMP) is None


def test_plot_files_and_cleanup(tmp_path):
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = T32._material(tmp_path, doping0=True)
    _set_mobility(cwd, {"ADP": 400.0, "POP": 200.0, "IMP": 1e47})
    plots = cwd.parent / "transport_plots"
    assert S.write_mechanism_resolved_plot(cwd, plots, "MoS2", True)
    d = plots / S.UNDOPED_DIR
    assert (d / "MoS2_mechanism_resolved_mobility.png").is_file()
    lines = [l for l in (d / "MoS2_mechanism_resolved_mobility.csv").read_text().splitlines() if not l.startswith("#")]
    rows = list(csv.DictReader(lines))
    assert {r["carrier"] for r in rows} == {"electron", "hole"}
    assert {r["mechanism"] for r in rows} == {"ADP", "POP", S.MATTHIESSEN, S.TOTAL_NOIMP}
    (cwd / "step8.4_amset2d" / "transport.json").unlink()                 # 没结果 -> 跳过并清旧图
    assert not S.write_mechanism_resolved_plot(cwd, plots, "MoS2", True)
    assert not (d / "MoS2_mechanism_resolved_mobility.png").exists()


def test_3d_without_subsets(tmp_path):
    """3D：只有 S8 transport.json（含 IMP、无去 IMP 子集）——电子不出总值（含 IMP 不算稀薄极限），
    空穴本征点 N_ii = 0，用 overall。"""
    S = _s83()
    cwd = tmp_path / "Si" / "ke-dft-cpu"
    (cwd / "step8_amset").mkdir(parents=True)
    dop, temps = [-1e17, 0.0, 1e17], [300.0, 600.0]
    mu = lambda v: [[T32._t(v * 300.0 / T, z=v * 300.0 / T) for T in temps] for _ in dop]
    (cwd / "step8_amset" / "transport.json").write_text(json.dumps({
        "doping": dop, "temperatures": temps,
        "mobility": {"overall": mu(150.0), "ADP": mu(500.0), "POP": mu(250.0), "IMP": mu(3000.0)}}))
    rows, why = S.mechanism_resolved_rows(cwd, False)
    assert why is None
    assert _pick(rows, "electron", "main", S.TOTAL_NOIMP) is None
    assert _pick(rows, "electron", "main", S.MATTHIESSEN)["mu_cm2/Vs"] == pytest.approx(1 / (1 / 500 + 1 / 250))
    assert _pick(rows, "hole", "main", S.TOTAL_NOIMP)["mu_cm2/Vs"] == pytest.approx(150.0)
    assert "n_2D_cm-2" not in rows[0] or rows[0]["n_2D_cm-2"] is None


# ---------------------------------------------------------------- [0054] 统一命名 / 本征点去重 / η 判据
def _write_ni(cwd, Ev, Ec, EF, run="step8.4_amset2d"):
    T = list(T32.TEMPS)
    (cwd / run / "undoped_ni.json").write_text(json.dumps(
        {"temperatures_K": T, "E_F_eV": [EF] * len(T), "n_i_cm-3": [1e5] * len(T),
         "E_v_eV": [Ev] * len(T), "E_c_eV": [Ec] * len(T)}))


def test_phys_label_naming():
    S = _s83()
    full = ["ADP", "IMP", "POP", "ODP"]
    assert S._phys_label(full, full) == "phonon + impurity (ADP+IMP+POP+ODP)"
    assert S._phys_label(["ADP", "POP", "ODP"], full) == "phonon-limited (ADP+POP+ODP)"
    assert S._phys_label(["ADP", "POP"], full) == "partial phonon (ADP+POP)"
    assert S._phys_label(["ADP"], full) == "ADP only (acoustic limit)"
    assert S._phys_label([], full) is None


def test_undoped_plot_dedupes_imp_and_reports_eta(tmp_path, capsys):
    """本征点 N_ii = 0：全机制（ADP+IMP+POP）与去 IMP（ADP+POP）重合 -> 只留去 IMP 那套；η 列按带边算。"""
    pytest.importorskip("matplotlib")
    S = _s83()
    cwd = T32._material(tmp_path, doping0=True)
    _write_ni(cwd, Ev=4.0, Ec=4.3, EF=4.15)            # 窄带隙 0.3 eV：600 K 时 η = 0.15/0.0517 ≈ 2.9 < 5
    plots = cwd.parent / "transport_plots"
    assert S.write_undoped_plots(cwd, plots, "MoS2", True)
    out = capsys.readouterr().out
    assert "η = min(E_c−E_F, E_F−E_v)/k_BT < 5" in out and "600K" in out
    lines = [l for l in (plots / S.UNDOPED_DIR / "MoS2_intrinsic_carrier.csv").read_text().splitlines()
             if not l.startswith("#")]
    rows = list(csv.DictReader(lines))
    assert "overall" not in {r["dataset"] for r in rows}             # 与 intrinsic 重合，去重
    assert "phonon-limited (ADP+POP)" in {r["label"] for r in rows}
    r600 = next(r for r in rows if float(r["T_K"]) == 600.0)
    assert float(r600["eta_c"]) == pytest.approx(0.15 / (S.KB_EV * 600.0), rel=1e-4)
    assert float(r600["f_CBM"]) == pytest.approx(1 / (1 + __import__("math").exp(0.15 / (S.KB_EV * 600.0))), rel=1e-3)


def test_mechanism_resolved_eta_on_undoped_rows(tmp_path):
    S = _s83()
    cwd = T32._material(tmp_path, doping0=True)
    _set_mobility(cwd, {"ADP": 400.0, "POP": 200.0, "IMP": 1e47})
    _write_ni(cwd, Ev=4.0, Ec=5.8, EF=4.9)             # 宽带隙 1.8 eV：η ≫ 5
    rows, _ = S.mechanism_resolved_rows(cwd, True)
    h = _pick(rows, "hole", "main", "ADP")
    assert h["eta"] == pytest.approx(0.9 / (S.KB_EV * 300.0), rel=1e-4)
    assert _pick(rows, "electron", "main", "ADP")["eta"] is None      # 掺杂点不记带边


def test_fermi_fix_writes_band_edges(tmp_path):
    F = _load("fermi_fix_0054", KE / "step8.4_amset2d" / "amset_fermi_fix.py")
    (tmp_path / "settings.yaml").write_text("x: 1\n")
    F.write_undoped_ni([(300.0, 4.9, 1e5, 4.0, 5.8)], run_dir=str(tmp_path))
    j = json.loads((tmp_path / "undoped_ni.json").read_text())
    assert j["E_v_eV"] == [4.0] and j["E_c_eV"] == [5.8]
    F.write_undoped_ni([(300.0, 4.9, 1e5)], run_dir=str(tmp_path))   # 旧格式记录（无带边）不崩
    assert json.loads((tmp_path / "undoped_ni.json").read_text())["E_c_eV"] == [None]


def test_old_names_gone():
    src = (KE / "step8.3_output" / "gen_step13_output.py").read_text(encoding="utf-8")
    assert '"%s_intrinsic.png" % M' not in src and '"%s_phonon_limited.png" % M' in src
    assert "amset (all but IMP)" not in src and "AMSET overall" not in src
