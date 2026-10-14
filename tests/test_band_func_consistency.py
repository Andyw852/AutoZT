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
    assert len(paths) == 1 and len(paths[0]) == len(labels) == 4   # 一条连通路径 Γ-M-K-Γ
    M = np.array(paths[0][1]) @ B
    K = np.array(paths[0][2]) @ B
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


def _plot(d, mass_v, mass_c, soc=False, nelect=2, kv=(0, 0, 0)):
    """抛物带边：价带 E=-ħ²k²/2m_v，导带 1+ħ²k²/2m_c，两侧对称，k 距离 0..0.5 Å⁻¹ 中心在 0.25。"""
    d.mkdir(parents=True)
    xs = np.linspace(0.0, 0.5, 51)   # 步长 0.01，带边 0.25 落在采样点上
    dk = xs - 0.25
    vb = -kc_h2() * dk ** 2 / mass_v
    cb = 1.0 + kc_h2() * dk ** 2 / mass_c
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
    assert abs(r["semilocal"]["CB"]["fwd"] - 0.1) < 1e-3
    assert abs(r["mass_ratio_hybrid_over_semilocal"]["CB_fwd"] - 1.5) < 1e-2
    assert abs(r["mass_ratio_hybrid_over_semilocal"]["VB_bwd"] - 1.0) < 1e-2


def test_mass_check_ok_and_valley_change(kc, tmp_path):
    m = tmp_path / "mat"
    _plot(m / "step2_bandgap" / "step2.2_pbe_plot", 0.2, 0.1)
    _plot(m / "step2_bandgap" / "step2.3_hse_plot", 0.21, 0.105)
    assert kc.hybrid_mass_check(m)["status"] == "ok"
    m2 = tmp_path / "mat2"
    _plot(m2 / "step2_bandgap" / "step2.2_pbe_plot", 0.2, 0.1)
    _plot(m2 / "step2_bandgap" / "step2.3_hse_plot", 0.2, 0.1, kv=(0.5, 0, 0))
    r = kc.hybrid_mass_check(m2)
    assert r["status"] == "warn" and any("谷次序" in x for x in r["reasons"])


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


# --------------------------------------------------------------------------- [G]
# 共享拟合引擎（fc_fit_driver 判据出图）与 fit-fc-thermal S2_plot：2D 走 za_2d 的 kz=0 路径
def _za():
    return _load("za_2d_under_test", ROOT / "skill" / "_common" / "za_2d.py")


@pytest.mark.parametrize("cell", [
    [[3.2, 0, 0], [-1.6, 2.771281, 0], [0, 0, 25.0]],          # 六方 120°
    [[3.2, 0, 0], [1.6, 2.771281, 0], [0, 0, 25.0]],           # 六方 60°
    [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 25.0]],                  # 正方
    [[3.0, 0, 0], [0, 4.5, 0], [0, 0, 25.0]],                  # 矩形
    [[3.0, 0, 0], [1.0, 2.8284, 0], [0, 0, 25.0]],             # 中心矩形（菱形原胞）
    [[3.0, 0, 0], [0.7, 4.1, 0], [0, 0, 25.0]],                # 斜方
])
@pytest.mark.parametrize("vac_first", [False, True])
def test_pool_band_path_matches_skill_copies(cell, vac_first):
    za = _za()
    cell = np.asarray(cell, float)
    vax = 2
    if vac_first:
        cell, vax = cell[[2, 0, 1]], 0
    ref = za.band_path_2d(cell, vax)
    for path in (ROOT / "skill" / "kl-dft-cpu" / "kl_common.py",
                 ROOT / "skill" / "phonon-dft-cpu" / "kl_common.py",
                 ROOT / "skill" / "_common" / "mlff" / "klmlff_common.py"):
        kc = _load("klc_cmp_%s" % path.parent.name.replace("-", "_"), path, [path.parent])
        got = kc.band_path_2d(cell, 51, vax)
        assert got[2] == ref[2] and got[1] == ref[1]
        assert np.allclose(np.asarray(got[0], float), np.asarray(ref[0], float))
    assert all(abs(p[vax]) < 1e-12 for seg in ref[0] for p in seg)    # kz=0


@pytest.mark.parametrize("gamma", [120.0, 60.0])
@pytest.mark.parametrize("vac_first", [False, True])
def test_write_band_yaml_2d_real_phonopy(tmp_path, gamma, vac_first):
    phonopy = pytest.importorskip("phonopy")
    from phonopy.structure.atoms import PhonopyAtoms
    za = _za()
    a = 2.5
    t = math.radians(gamma)
    cell = np.array([[a, 0, 0], [a * math.cos(t), a * math.sin(t), 0], [0, 0, 20.0]])
    pos = [[0, 0, 0.5], [1 / 3, 2 / 3, 0.5]] if gamma == 120.0 else [[0, 0, 0.5], [1 / 3, 1 / 3, 0.5]]
    scm = np.diag([3, 3, 1])
    cart = 2
    if vac_first:                                          # 真空放在第 0 个基矢（Cartesian 仍是 z）
        cell = cell[[2, 0, 1]]
        pos = [[p[2], p[0], p[1]] for p in pos]
        scm = np.diag([1, 3, 3])
    uc = PhonopyAtoms(symbols=["B", "N"], cell=cell, scaled_positions=pos)
    ph = phonopy.Phonopy(uc, supercell_matrix=scm)
    n = len(ph.supercell)
    ph.force_constants = np.zeros((n, n, 3, 3))
    lat = za.write_band_yaml_2d(ph, tmp_path / "band.yaml", npoints=11, cart_axis=cart)
    assert lat == "hexagonal" and (tmp_path / "band.yaml").is_file()
    bd = ph.get_band_structure_dict()
    B = 2 * np.pi * np.linalg.inv(np.asarray(ph.primitive.cell, float)).T
    qcart = np.concatenate([np.asarray(q, float) for q in bd["qpoints"]]) @ B
    assert np.max(np.abs(qcart[:, 2])) < 1e-9                          # 没有沿真空方向的线段
    kmax = max(np.linalg.norm(q) for q in qcart)
    assert abs(kmax - 4 * np.pi / (3 * a)) < 1e-6                     # 走到真 K（BZ 角）


def test_engine_and_fit_plot_use_2d_path():
    eng = (ROOT / "skill" / "_common" / "fcfit" / "fc_fit_driver.py").read_text(encoding="utf-8")
    i = eng.index("def _stability_gate(")
    body = eng[i:eng.index("\ndef ", i + 10)]
    assert "za_2d.write_band_yaml_2d(ph_nonac" in body
    assert body.index("if is2d:\n            # 2026-10-14") < body.index("auto_band_structure")
    plot = (ROOT / "skill" / "fit-fc-thermal" / "fc_plot_phonon.py").read_text(encoding="utf-8")
    assert "za_2d.write_band_yaml_2d(ph, dst" in plot


# --------------------------------------------------------------------------- [H]
# 2026-10-14 第三批：lattice_kappa / cutoff_conv_driver 的 2D 路径 + 标签不再被 phonopy 丢掉
def _hex_phonopy(gamma=60.0, vac_first=True, scm_inplane=3):
    phonopy = pytest.importorskip("phonopy")
    from phonopy.structure.atoms import PhonopyAtoms
    a = 2.5
    t = math.radians(gamma)
    cell = np.array([[a, 0, 0], [a * math.cos(t), a * math.sin(t), 0], [0, 0, 20.0]])
    pos = [[0, 0, 0.5], [1 / 3, 2 / 3, 0.5]] if gamma == 120.0 else [[0, 0, 0.5], [1 / 3, 1 / 3, 0.5]]
    scm = np.diag([scm_inplane, scm_inplane, 1])
    if vac_first:
        cell = cell[[2, 0, 1]]
        pos = [[p[2], p[0], p[1]] for p in pos]
        scm = np.diag([1, scm_inplane, scm_inplane])
    uc = PhonopyAtoms(symbols=["B", "N"], cell=cell, scaled_positions=pos)
    ph = phonopy.Phonopy(uc, supercell_matrix=scm, primitive_matrix="auto")
    n = len(ph.supercell)
    ph.force_constants = np.zeros((n, n, 3, 3))
    return ph, a


def _qcart(ph):
    bd = ph.get_band_structure_dict()
    B = 2 * np.pi * np.linalg.inv(np.asarray(ph.primitive.cell, float)).T
    return np.concatenate([np.asarray(q, float) for q in bd["qpoints"]]) @ B


def test_band_labels_survive_in_phonopy(tmp_path):
    ph, a = _hex_phonopy()
    za = _za()
    za.write_band_yaml_2d(ph, tmp_path / "band.yaml", npoints=11, cart_axis=2)
    assert ph.band_structure.labels == ["$\\Gamma$", "M", "K", "$\\Gamma$"]
    assert "labels:" in (tmp_path / "band.yaml").read_text(encoding="utf-8")


@pytest.fixture
def lk():
    pytest.importorskip("phonopy")
    return _load("lattice_kappa_under_test", ROOT / "skill" / "kl-dft-cpu" / "lattice_kappa.py",
                 [ROOT / "skill" / "kl-dft-cpu"])


@pytest.mark.parametrize("gamma", [120.0, 60.0])
@pytest.mark.parametrize("vac_first", [False, True])
def test_lattice_kappa_2d_path(lk, tmp_path, gamma, vac_first):
    ph, a = _hex_phonopy(gamma, vac_first)
    assert lk._vacuum_axis_2d(ph, {}) is not None
    lk._resolve_band_path(ph, tmp_path / "band.yaml", npoints=11, C={})
    q = _qcart(ph)
    assert np.max(np.abs(q[:, 2])) < 1e-9                       # Cartesian 真空 z：无 kz 线段
    assert abs(max(np.linalg.norm(q, axis=1)) - 4 * np.pi / (3 * a)) < 1e-6
    pos, lab = lk._band_ticks(ph.band_structure)                 # 原来 labels=None 会在这里出错
    assert lab == ["Γ", "M", "K", "Γ"]


def test_lattice_kappa_vacuum_axis_rules(lk):
    ph, _ = _hex_phonopy()
    assert lk._vacuum_axis_2d(ph, {"dim": "3d"}) is None           # 显式 3D 不改
    phonopy = pytest.importorskip("phonopy")
    from phonopy.structure.atoms import PhonopyAtoms
    uc = PhonopyAtoms(symbols=["Si", "Si"], cell=np.eye(3) * 5.43,
                      scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]])
    bulk = phonopy.Phonopy(uc, supercell_matrix=np.eye(3, dtype=int), primitive_matrix="P")
    assert lk._vacuum_axis_2d(bulk, {}) is None                    # 体相没有 8 Å 真空


def test_lattice_kappa_pymatgen_fallback_refuses_nonstandard(lk, tmp_path, monkeypatch):
    pytest.importorskip("pymatgen.core")
    phonopy = pytest.importorskip("phonopy")
    from phonopy.structure.atoms import PhonopyAtoms
    a, c, t = 3.0, 5.0, math.radians(60)                           # 体相六方，γ=60° 非标准原胞
    uc = PhonopyAtoms(symbols=["Si", "Si"],
                      cell=[[a, 0, 0], [a * math.cos(t), a * math.sin(t), 0], [0, 0, c]],
                      scaled_positions=[[0, 0, 0], [1 / 3, 1 / 3, 0.5]])
    ph = phonopy.Phonopy(uc, supercell_matrix=np.eye(3, dtype=int), primitive_matrix="P")
    n = len(ph.supercell)
    ph.force_constants = np.zeros((n, n, 3, 3))

    def _no_seekpath(*a, **k):
        raise ImportError("no seekpath")
    monkeypatch.setattr(ph, "auto_band_structure", _no_seekpath)
    with pytest.raises(RuntimeError, match="seekpath"):
        lk._resolve_band_path(ph, tmp_path / "band.yaml", npoints=11, C={})


def test_lattice_kappa_3d_uses_seekpath(lk, tmp_path):
    pytest.importorskip("seekpath")
    phonopy = pytest.importorskip("phonopy")
    from phonopy.structure.atoms import PhonopyAtoms
    uc = PhonopyAtoms(symbols=["Si"] * 2, cell=[[0, 2.715, 2.715], [2.715, 0, 2.715], [2.715, 2.715, 0]],
                      scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]])
    ph = phonopy.Phonopy(uc, supercell_matrix=np.eye(3, dtype=int) * 2, primitive_matrix="P")
    n = len(ph.supercell)
    ph.force_constants = np.zeros((n, n, 3, 3))
    lk._resolve_band_path(ph, tmp_path / "band.yaml", npoints=11, C={})
    assert any("X" in str(x) for x in ph.band_structure.labels)


def _run_driver(script, workdir, extra_paths):
    import os
    import subprocess
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(str(p) for p in extra_paths),
               PYTHONWARNINGS="ignore")
    return subprocess.run([sys.executable, str(script)], cwd=str(workdir), env=env,
                          capture_output=True, text=True, timeout=600)


def _bn_dataset(d, gamma=60.0):
    """γ=60° 六方 BN 单层（真空在第 0 基矢 = Cartesian x），1×4×4 超胞，LJ 势 6 帧随机位移。"""
    pytest.importorskip("ase")
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from ase.io import write
    a, t = 2.5, math.radians(gamma)
    cell = np.array([[20.0, 0, 0], [0, a, 0], [0, a * math.cos(t), a * math.sin(t)]])
    uc = Atoms("BN", cell=cell, scaled_positions=[[0.5, 0, 0], [0.5, 1 / 3, 1 / 3]], pbc=True)
    write(str(d / "POSCAR"), uc, format="vasp", direct=True)
    sc = uc.repeat((1, 4, 4))
    rng = np.random.default_rng(0)
    calc = LennardJones(sigma=1.3, epsilon=1.0, rc=4.0, smooth=True)
    D, F = [], []
    for _ in range(6):
        s = sc.copy()
        u = rng.normal(0, 0.01, (len(sc), 3))
        s.positions += u
        s.calc = calc
        D.append(u)
        F.append(s.get_forces())
    np.save(str(d / "disps.npy"), np.array(D))
    np.save(str(d / "forces.npy"), np.array(F))
    return a, len(sc)


def test_cutoff_conv_driver_2d_end_to_end(tmp_path):
    pytest.importorskip("symfc")
    pytest.importorskip("phonopy")
    import yaml
    skill = tmp_path / "mat" / "phonon-mlff-cpu"
    step = skill / "step4_cutoff_conv"
    step.mkdir(parents=True)
    for f in (ROOT / "skill" / "_common" / "za_2d.py", OPT / "dim_common.py"):  # gen_need 落在技能目录
        (skill / f.name).write_bytes(f.read_bytes())
    (step / "cutoff_conv_driver.py").write_bytes(
        (ROOT / "skill" / "phonon-mlff-cpu" / "cutoff_conv_driver.py").read_bytes())
    a, nsc = _bn_dataset(step)
    (step / "klmlff_params.txt").write_text("SUPERCELL = 1 4 4\nCUTOFF_LIST = 4.0\n")
    r = _run_driver(step / "cutoff_conv_driver.py", step, [])        # 不靠 PYTHONPATH，验证上一级入路径
    assert r.returncode == 0, r.stdout + r.stderr
    s = json.loads((step / "cutoff_conv_summary.json").read_text(encoding="utf-8"))
    assert s["dim"] == "2d" and s["supercell_atoms"] == nsc and "error" not in s["results"][0]
    d = yaml.safe_load((step / "band-cut4.yaml").read_text(encoding="utf-8"))
    q = np.array([p["q-position"] for p in d["phonon"]]) @ np.array(d["reciprocal_lattice"]) * 2 * np.pi
    assert np.max(np.abs(q[:, 0])) < 1e-9                             # 真空 = Cartesian x
    assert abs(np.linalg.norm(q, axis=1).max() - 4 * np.pi / (3 * a)) < 1e-6
    assert d.get("labels")


def test_phonon_mlff_s3_driver_runs(tmp_path):
    """phonon-mlff-cpu S3：原来引用未定义的 dim（NameError），n_sc 也恒为 0。"""
    pytest.importorskip("symfc")
    pytest.importorskip("phonopy")
    step = tmp_path / "step3_phonon"
    step.mkdir()
    (step / "phonon_fit_driver.py").write_bytes(
        (ROOT / "skill" / "phonon-mlff-cpu" / "phonon_fit_driver.py").read_bytes())
    a, nsc = _bn_dataset(step)
    (step / "klmlff_params.txt").write_text("SUPERCELL = 1 4 4\nFC_CUTOFF = 4.0\n")
    r = _run_driver(step / "phonon_fit_driver.py", step,
                    [ROOT / "skill" / "_common" / "mlff", ROOT / "skill" / "_common", OPT])
    out = r.stdout + r.stderr
    assert "NameError" not in out, out
    assert "%d 原子" % nsc in out, out
    assert "kz=0" in out                                              # 2D 路径分支
    # 真空在 a（下标 0）：原来 `vac_axis or 2` 把它当成 c，网格成了 60 60 1、路径判成矩形
    assert "hexagonal 格子" in out, out
    assert "mesh=1 60 60" in out, out


# --------------------------------------------------------------------------- [I]
# relax_common 的 2D 层厚：原来 `int(vac_axis or 2)`，真空在 a（下标 0）时按 c 算
def _mos2_like(path, vac_first):
    """六方单层（Mo + 上下 S），真空 20 Å。vac_first=True 时真空放在第 0 个基矢。"""
    a, c = 3.18, 20.0
    lat = [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, c]]
    fr = [[0, 0, 0.5], [1 / 3, 2 / 3, 0.5 + 1.56 / c], [1 / 3, 2 / 3, 0.5 - 1.56 / c]]
    if vac_first:
        lat = [lat[2], lat[0], lat[1]]
        fr = [[f[2], f[0], f[1]] for f in fr]
    txt = ["MoS2", "1.0"] + ["  %.8f %.8f %.8f" % tuple(v) for v in lat] + \
          ["Mo S", "1 2", "Direct"] + ["  %.8f %.8f %.8f" % tuple(f) for f in fr]
    path.write_text("\n".join(txt) + "\n")
    return 0 if vac_first else 2


def test_relax_common_layer_factor_vacuum_on_a(tmp_path):
    pytest.importorskip("ase")
    rc = _load("relax_common_layer_test", OPT / "relax_common.py")
    pc, pa = tmp_path / "c.vasp", tmp_path / "a.vasp"
    ax_c, ax_a = _mos2_like(pc, False), _mos2_like(pa, True)
    fc, fa = rc._layer_factor(pc, ax_c), rc._layer_factor(pa, ax_a)
    assert fc is not None and fa is not None
    assert abs(fc - fa) < 1e-6, (fc, fa)          # 同一个单层，真空放 a 或 c 结果必须相同


def test_relax_common_material_thickness_vacuum_on_a(tmp_path):
    pytest.importorskip("ase")
    rc = _load("relax_common_thick_test", OPT / "relax_common.py")
    res = {}
    for vac_first in (False, True):
        mat = tmp_path / ("mat_a" if vac_first else "mat_c")
        out = mat / "kl-dft-cpu" / "step1_std_opt"
        out.mkdir(parents=True)
        (mat / "POSCAR").write_text("x\n")
        ax = _mos2_like(out / "POSCAR", vac_first)
        rc.write_material_thickness(out, "2d", ax)
        js = json.loads((mat / "thickness_2d.json").read_text(encoding="utf-8"))
        assert js["vac_axis"] == ax                       # 记录的轴 = 实际算几何用的轴
        res[vac_first] = js
    for k in ("thickness_d_A", "h_perp_A"):
        assert abs(float(res[True][k]) - float(res[False][k])) < 1e-6, (k, res)
