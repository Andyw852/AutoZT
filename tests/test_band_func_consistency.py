# -*- coding: utf-8 -*-
"""2026-10-14：ke/kl 能带与泛函设置修复的回归测试。

  [A] 2D 六方声子路径：γ=60° 原胞的 K 点（三份同源实现：kl-dft-cpu / phonon-dft-cpu / mlff 公共池）
  [B] S1 的 FUNC=auto：本技能已有 > 同材料配对技能（ke<->kl）> 键网络判据；决策写进 FUNC_REASON
  [C] ke S2.3 的 HYBRID（pbe0 不屏蔽、custom 校验、预设）
  [D] ke S8/S8.4 的能带一致性检查：SOC 影响闸门、杂化/半局域 m* 核对（只告警）
  [E] zt S20 的 ke/kl 泛函一致性核对
  [F] stepconf 保留键：HYBRID* / FUNC_ALIGN 写在材料共用 step.conf 不会让别的步骤报错

不依赖 VASP / phonopy / pymatgen（pymatgen 用假模块顶替）。
"""
import importlib.util
import json
import math
import sys
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
OPT = ROOT / "skill" / "_common" / "opt"
KE = ROOT / "skill" / "ke-dft-cpu"


def _load(name, path, extra_paths=()):
    for p in [str(OPT), str(ROOT / "skill" / "_common")] + [str(x) for x in extra_paths]:
        if p not in sys.path:
            sys.path.insert(0, p)
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- [A]
@pytest.mark.parametrize("path", [
    ROOT / "skill" / "kl-dft-cpu" / "kl_common.py",
    ROOT / "skill" / "phonon-dft-cpu" / "kl_common.py",
    ROOT / "skill" / "_common" / "mlff" / "klmlff_common.py",
])
@pytest.mark.parametrize("gamma", [120.0, 60.0])
@pytest.mark.parametrize("vac_first", [False, True])
def test_hex_2d_k_point(path, gamma, vac_first):
    kc = _load("klc_%s_%d_%d" % (path.parent.name.replace("-", "_"), gamma, vac_first),
               path, [path.parent])
    a = 3.2
    t = math.radians(gamma)
    cell = np.array([[a, 0, 0], [a * math.cos(t), a * math.sin(t), 0], [0, 0, 25.0]])
    vax = 2
    if vac_first:                       # 生产 MoS₂：真空在第 0 基矢
        cell, vax = cell[[2, 0, 1]], 0
    paths, labels, lat = kc.band_path_2d(cell, 51, vax)
    assert lat == "hexagonal"
    B = 2 * np.pi * np.linalg.inv(cell).T
    M = np.array(paths[0][1]) @ B
    K = np.array(paths[1][1]) @ B
    assert abs(np.linalg.norm(M) - 2 * np.pi / (math.sqrt(3) * a)) < 1e-9
    assert abs(np.linalg.norm(K) - 4 * np.pi / (3 * a)) < 1e-9       # 真 K（BZ 角）
    assert abs(np.linalg.norm(K - M) - 2 * np.pi / (3 * a)) < 1e-9   # 与 M 相邻


# --------------------------------------------------------------------------- [B]
@pytest.fixture
def rc(monkeypatch):
    mod = _load("relax_common_under_test", OPT / "relax_common.py")
    fake_pmg = types.ModuleType("pymatgen")
    fake_core = types.ModuleType("pymatgen.core")

    class _S:
        @staticmethod
        def from_file(_p):
            return object()
    fake_core.Structure = _S
    monkeypatch.setitem(sys.modules, "pymatgen", fake_pmg)
    monkeypatch.setitem(sys.modules, "pymatgen.core", fake_core)
    monkeypatch.setattr(mod, "detect_dimension", lambda *a, **k: ("3d", None, ""))
    import method_select
    monkeypatch.setattr(method_select, "default_method",
                        lambda *a, **k: ("pbe-d3", "键网络 2D < 结构 3D（假）"))
    return mod


def _mat(tmp_path, skill):
    d = tmp_path / "Bi2Te3" / skill
    d.mkdir(parents=True)
    (d / "POSCAR").write_text("x\n")
    return d


def _resolve(rc, cwd, monkeypatch, params=None):
    monkeypatch.chdir(cwd)
    monkeypatch.setattr(rc, "STEP_PARAMS", dict({"FUNC": "auto"}, **(params or {})))
    rc.FUNC_REASON = ""
    return rc.resolve_func(Path("incar_3d.tpl"), "step1_opt")


def test_auto_uses_bond_network_without_sibling(rc, tmp_path, monkeypatch):
    d = _mat(tmp_path, "ke-dft-cpu")
    func, why = _resolve(rc, d, monkeypatch)
    assert func == "pbe-d3" and "键网络" in rc.FUNC_REASON


def test_auto_follows_sibling_kl(rc, tmp_path, monkeypatch, capsys):
    d = _mat(tmp_path, "ke-dft-cpu")
    s = tmp_path / "Bi2Te3" / "kl-dft-cpu" / "step1_std_opt"
    s.mkdir(parents=True)
    (s / "workflow_method.txt").write_text("FUNC=pbesol\nGGA=PS\n")
    func, _ = _resolve(rc, d, monkeypatch)
    assert func == "pbesol"                         # 跟随配对技能
    assert "配对技能" in rc.FUNC_REASON
    assert "★" in capsys.readouterr().out           # 与键网络判据（pbe-d3）不同 -> 告警


def test_auto_sibling_in_local_result_dir(rc, tmp_path, monkeypatch):
    d = _mat(tmp_path, "kl-dft-cpu")
    s = tmp_path / "Bi2Te3" / "ke-dft-cpu" / "result" / "step1_opt"
    s.mkdir(parents=True)
    (s / "workflow_method.txt").write_text("FUNC=pbe-d3\n")
    assert _resolve(rc, d, monkeypatch)[0] == "pbe-d3"


def test_auto_align_can_be_disabled(rc, tmp_path, monkeypatch):
    d = _mat(tmp_path, "ke-dft-cpu")
    s = tmp_path / "Bi2Te3" / "kl-dft-cpu" / "step1_std_opt"
    s.mkdir(parents=True)
    (s / "workflow_method.txt").write_text("FUNC=pbesol\n")
    assert _resolve(rc, d, monkeypatch, {"FUNC_ALIGN": False})[0] == "pbe-d3"


def test_auto_keeps_own_previous_func(rc, tmp_path, monkeypatch):
    d = _mat(tmp_path, "ke-dft-cpu")
    (d / "step1_opt").mkdir()
    (d / "step1_opt" / "workflow_method.txt").write_text("FUNC=pbesol\n")
    s = tmp_path / "Bi2Te3" / "kl-dft-cpu" / "step1_std_opt"
    s.mkdir(parents=True)
    (s / "workflow_method.txt").write_text("FUNC=pbe-d3\n")
    assert _resolve(rc, d, monkeypatch)[0] == "pbesol"     # retry 不换泛函


def test_explicit_hybrid_func_rejected(rc, tmp_path, monkeypatch):
    d = _mat(tmp_path, "ke-dft-cpu")
    with pytest.raises(SystemExit) as e:
        _resolve(rc, d, monkeypatch, {"FUNC": "pbe0"})
    assert "HYBRID" in str(e.value)


def test_ke_and_kl_global_func_both_auto():
    def func_of(p):
        for ln in Path(p).read_text(encoding="utf-8").splitlines():
            s = ln.split("#", 1)[0].strip()
            if s.startswith("FUNC") and "=" in s and s.split("=")[0].strip() == "FUNC":
                return s.split("=", 1)[1].strip()
    assert func_of(KE / "step.conf") == "auto"
    assert func_of(ROOT / "skill" / "kl-dft-cpu" / "templates" / "step.conf") == "auto"


def test_root_and_pool_relax_common_in_sync():
    assert (ROOT / "relax_common.py").read_bytes() == (OPT / "relax_common.py").read_bytes()


# --------------------------------------------------------------------------- [C]
@pytest.fixture
def hse():
    return _load("gen_step4_hse_under_test",
                 KE / "step2_bandgap" / "step2.3_hse" / "gen_step4_HSE.py", [KE])


def _conf(tmp_path, body):
    (tmp_path / "step.conf").write_text("[params]\nFUNC = auto\n" + body, encoding="utf-8")
    return tmp_path


def test_hybrid_default_is_auto(hse, tmp_path):
    h = hse.resolve_hybrid(_conf(tmp_path, ""), None)
    assert h["name"] == "auto" and h["gga"] is None and h["hfscreen"] == 0.2


def test_hybrid_pbe0_unscreened(hse, tmp_path):
    h = hse.resolve_hybrid(_conf(tmp_path, "HYBRID = pbe0\n"), None)
    assert h["gga"] == "PE" and h["aexx"] == 0.25 and h["hfscreen"] is None
    assert hse.hybrid_label(h, "PE") == "PBE0"
    eff = hse.effective_incar([("LHFCALC", ".TRUE."), ("AEXX", "0.25")], ["HFSCREEN"], {})
    assert hse.resolve_hfrcut(eff)[0] == "-1"       # 非屏蔽 + KPOINTS_OPT -> HFRCUT=-1


def test_hybrid_custom(hse, tmp_path):
    h = hse.resolve_hybrid(_conf(tmp_path, "HYBRID = custom\nHYBRID_GGA = ps\nHYBRID_AEXX = 0.25\n"), None)
    assert (h["gga"], h["aexx"], h["hfscreen"]) == ("PS", 0.25, None)
    assert hse.hybrid_label(h, "PS") == "PBEsol0"
    h = hse.resolve_hybrid(_conf(tmp_path, "HYBRID = custom\nHYBRID_AEXX = 0.3\nHYBRID_HFSCREEN = 0.2\n"), None)
    assert h["gga"] is None and h["hfscreen"] == 0.2


@pytest.mark.parametrize("body", ["HYBRID = custom\n", "HYBRID = custom\nHYBRID_AEXX = 1.5\n",
                                  "HYBRID = b3lyp\n",
                                  "HYBRID = custom\nHYBRID_AEXX = 0.25\nHYBRID_HFSCREEN = -1\n"])
def test_hybrid_bad_input(hse, tmp_path, body):
    with pytest.raises(SystemExit):
        hse.resolve_hybrid(_conf(tmp_path, body), None)


def test_hybrid_cli_wins(hse, tmp_path):
    h = hse.resolve_hybrid(_conf(tmp_path, "HYBRID = pbe0\n"), "hse06")
    assert h["name"] == "hse06" and h["hfscreen"] == 0.2


# --------------------------------------------------------------------------- [D]
@pytest.fixture
def kc():
    return _load("ke_common_under_test", KE / "ke_common.py", [KE])


def _eigenval(path, nelect, kpts, bands, ispin=1):
    lines = ["  1  1  1  %d" % ispin, " 0", " 0", "  CAR", " x",
             "  %d  %d  %d" % (nelect, len(kpts), len(bands[0]))]
    for k, eb in zip(kpts, bands):
        lines += ["", "  %.7f  %.7f  %.7f  1.0" % tuple(k)]
        lines += ["  %d  %.6f  1.0" % (i + 1, e) for i, e in enumerate(eb)]
    Path(path).write_text("\n".join(lines) + "\n")


def _soc_case(tmp_path, split):
    """非 SOC：Γ 价带顶二重（各装 2 电子）；SOC：被劈开 split eV。nelect=4。"""
    m = tmp_path / "mat"
    b = m / "step2_bandgap" / "step2.2_pbe"
    u = m / "step3_uniform"
    b.mkdir(parents=True)
    u.mkdir(parents=True)
    (b / "INCAR").write_text("LSORBIT = .TRUE.\n")
    (u / "INCAR").write_text("ISPIN = 1\n")
    kp = [(0, 0, 0), (0.5, 0, 0)]
    _eigenval(u / "EIGENVAL", 4, kp + [(0.25, 0, 0)],
              [[-1.0, 0.0, 1.0], [-1.5, -0.4, 1.3], [-1.2, -0.2, 1.1]])
    _eigenval(b / "EIGENVAL", 4, kp,
              [[-1.0, -1.0, 0.0, -split, 1.0, 1.0], [-1.5, -1.5, -0.4, -0.4, 1.3, 1.3]])
    return m


def test_soc_gate_flags_large_splitting(kc, tmp_path):
    r = kc.soc_impact_check(_soc_case(tmp_path, 0.2))
    assert r["status"] == "warn" and r["n_common_k"] == 2
    assert abs(r["vb_level_shift_max_eV"] - 0.2) < 1e-6


def test_soc_gate_ok_small_splitting(kc, tmp_path):
    r = kc.soc_impact_check(_soc_case(tmp_path, 0.01))
    assert r["status"] == "ok", r


def test_soc_gate_consistent_when_both_nosoc(kc, tmp_path):
    m = _soc_case(tmp_path, 0.2)
    (m / "step2_bandgap" / "step2.2_pbe" / "INCAR").write_text("ISPIN = 1\n")
    assert kc.soc_impact_check(m)["status"] == "consistent"


def _plot(d, mass_v, mass_c, soc=False, nelect=2, kv=(0, 0, 0), xs=None, cv=0.25, alpha=0.0):
    """抛物带边：价带 E=-ħ²k²/2m_v，导带 1+ħ²k²/2m_c，两侧对称，k 距离 0..0.5 Å⁻¹，价带顶在 cv、导带底在 0.25。
    xs：路径点（默认步长 0.01）；alpha > 0：导带用 Kane 非抛物色散 E(1+αE) = ħ²k²/2m。"""
    d.mkdir(parents=True)
    xs = np.linspace(0.0, 0.5, 51) if xs is None else np.asarray(xs, float)
    vb = -kc_h2() * (xs - cv) ** 2 / mass_v
    g = kc_h2() * (xs - 0.25) ** 2 / mass_c
    cb = 1.0 + (g if alpha <= 0 else (-1 + np.sqrt(1 + 4 * alpha * g)) / (2 * alpha))
    lines = ["# x", "# y"] + ["%12.6f %10.6f %10.6f" % (x, v, c) for x, v, c in zip(xs, vb, cb)]
    (d / "band.dat").write_text("\n".join(lines) + "\n")
    (d / "band_summary.json").write_text(json.dumps({
        "status": "ok", "gap_eV": 1.0, "nelect": nelect, "soc": soc, "nspin": 1,
        "vbm": {"k_frac": list(kv)}, "cbm": {"k_frac": [0, 0, 0]},
        "files": {"dat": "/remote/x/band.dat"}}))


def kc_h2():
    return 3.80998


def test_mass_check_flags_heavier_hybrid_mass(kc, tmp_path):
    m = tmp_path / "mat"
    _plot(m / "step2_bandgap" / "step2.2_pbe_plot", 0.2, 0.1)
    _plot(m / "step2_bandgap" / "step2.3_hse_plot", 0.2, 0.15)     # 导带 m* 杂化重 50%
    r = kc.hybrid_mass_check(m)
    assert r["status"] == "warn"
    assert abs(r["fits"]["CB"]["semilocal"]["fwd"] - 0.1) < 1e-3
    assert abs(r["mass_ratio_hybrid_over_semilocal"]["CB_fwd"] - 1.5) < 1e-2
    assert abs(r["mass_ratio_hybrid_over_semilocal"]["VB_bwd"] - 1.0) < 1e-2


def test_mass_check_ok_and_valley_change(kc, tmp_path):
    m = tmp_path / "mat"
    _plot(m / "step2_bandgap" / "step2.2_pbe_plot", 0.2, 0.1)
    _plot(m / "step2_bandgap" / "step2.3_hse_plot", 0.21, 0.105)
    assert kc.hybrid_mass_check(m)["status"] == "ok"
    m2 = tmp_path / "mat2"
    _plot(m2 / "step2_bandgap" / "step2.2_pbe_plot", 0.2, 0.1)
    _plot(m2 / "step2_bandgap" / "step2.3_hse_plot", 0.2, 0.1, kv=(0.5, 0, 0), cv=0.45)   # 杂化价带顶真的挪了
    r = kc.hybrid_mass_check(m2)
    assert r["status"] == "warn" and any("谷次序" in x for x in r["reasons"])


def test_mass_check_downsampled_path_no_false_alarm(kc, tmp_path):
    """杂化路径是降采样的（每段几个点），带边在两个采样点之间：同一条能带不能报"谷次序变了"，
    非抛物（Kane）带也不能因为两边窗口不同报质量比。"""
    full = np.linspace(0.0, 0.5, 51)
    sub = full[::7]                                                    # 0, 0.07, ..., 0.49：不含 0.25
    for i, alpha in enumerate((0.0, 2.0)):
        m = tmp_path / ("mat%d" % i)
        _plot(m / "step2_bandgap" / "step2.2_pbe_plot", 0.2, 0.1, xs=full, cv=0.37, alpha=alpha,
              kv=(0.37, 0, 0))
        _plot(m / "step2_bandgap" / "step2.3_hse_plot", 0.2, 0.1, xs=sub, cv=0.37, alpha=alpha,
              kv=(0.35, 0, 0))
        r = kc.hybrid_mass_check(m)
        assert r["status"] in ("ok", "not_applicable"), r
        assert "同一组 k 点" in r["basis"]
        for k, v in r["mass_ratio_hybrid_over_semilocal"].items():
            assert abs(v - 1.0) < 1e-6, (alpha, k, v)


def test_soc_gate_inconclusive_when_edge_not_shared(kc, tmp_path):
    """S3 的真实带边（k=0.25）不在 S2.2 的自洽网格里：不能说"SOC 影响小"。"""
    m = _soc_case(tmp_path, 0.01)
    u = m / "step3_uniform"
    _eigenval(u / "EIGENVAL", 4, [(0, 0, 0), (0.5, 0, 0), (0.25, 0, 0)],
              [[-1.0, 0.0, 1.0], [-1.5, -0.4, 1.3], [-1.2, 0.3, 0.6]])         # 带边在 0.25
    r = kc.soc_impact_check(m)
    assert r["status"] == "inconclusive", r
    assert any("VBM" in x for x in r["missed_edges"]) and any("CBM" in x for x in r["missed_edges"])
    res = kc.band_consistency_checks(m, m / "step8_amset", "S8")
    assert res["soc_impact"]["status"] == "inconclusive"


def test_band_consistency_writes_json(kc, tmp_path):
    m = _soc_case(tmp_path, 0.2)
    out = m / "step8_amset"
    res = kc.band_consistency_checks(m, out, "S8")
    assert res["soc_impact"]["status"] == "warn"
    assert res["hybrid_mass"]["status"] == "not_applicable"
    assert json.loads((out / "band_consistency.json").read_text(encoding="utf-8"))


def test_s8_and_s84_call_band_consistency():
    for f in (KE / "step8_amset" / "gen_step10_amset.py",
              KE / "step8.4_amset2d" / "gen_step14_amset2d.py"):
        assert "kc.band_consistency_checks(cwd, out," in f.read_text(encoding="utf-8")


# --------------------------------------------------------------------------- [E]
@pytest.fixture
def zt():
    return _load("gen_step20_zt_under_test",
                 ROOT / "skill" / "zt-dft-cpu" / "step20_zt" / "gen_step20_zt.py",
                 [ROOT / "skill" / "zt-dft-cpu" / "step20_zt"])


def _zt_inputs(tmp_path, ke_func, kl_func):
    (tmp_path / "inputs" / "structure").mkdir(parents=True)
    (tmp_path / "inputs" / "kappa_L").mkdir(parents=True)
    (tmp_path / "inputs" / "structure" / "workflow_method.txt").write_text("FUNC=%s\n" % ke_func)
    if kl_func:
        (tmp_path / "inputs" / "kappa_L" / "kl_method.json").write_text(json.dumps({"func": kl_func}))
    return tmp_path


def test_zt_functional_check(zt, tmp_path):
    assert zt.functional_check(_zt_inputs(tmp_path / "a", "pbesol", "pbesol"))["ok"] is True
    r = zt.functional_check(_zt_inputs(tmp_path / "b", "pbesol", "pbe-d3"))
    assert r["ok"] is False and "★" in r["note"]
    assert zt.functional_check(_zt_inputs(tmp_path / "c", "pbesol", None))["ok"] is None


def test_kl_s6_writes_kl_method(tmp_path, monkeypatch):
    sys.modules.setdefault("gen_step2_kappa", types.ModuleType("gen_step2_kappa"))
    g6 = _load("gen_step6_kappa_under_test", ROOT / "skill" / "kl-dft-cpu" / "gen_step6_kappa.py",
               [ROOT / "skill" / "kl-dft-cpu"])
    (tmp_path / "step1_std_opt").mkdir()
    (tmp_path / "step1_std_opt" / "workflow_method.txt").write_text(
        "FUNC=pbe-d3\nGGA=PE\nFUNC_REASON=auto：层状\nDIM=3D\n")
    rec = g6.write_kl_method(tmp_path)
    assert rec["func"] == "pbe-d3"
    assert json.loads((tmp_path / "step6_kappa" / "kl_method.json").read_text())["func"] == "pbe-d3"


# --------------------------------------------------------------------------- [F]
@pytest.mark.parametrize("path", [ROOT / "stepconf.py", OPT / "stepconf.py"])
def test_reserved_keys(path):
    sc = _load("stepconf_%d" % (abs(hash(str(path))) % 10000), path)
    for k in ("HYBRID", "HYBRID_GGA", "HYBRID_AEXX", "HYBRID_HFSCREEN", "FUNC_ALIGN"):
        assert k in sc.RESERVED_PARAMS
