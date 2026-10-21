"""[0062] S8.3 / S8.1 自动取 κ_L：2D 两链胞高不同（kl Lz ≠ ke c）时按片热导守恒 ×Lz/c 换到 ke 元胞口径。

以前只打印"元胞不一致"告警、原样拿 kl 胞的原始 κ_L 去配 ke 胞的 σ/κ_e，κ_L 放大 c/Lz 倍、
ZT 偏低（MoS2：Lz 25.0 vs c 34.879）；S8.1 则直接 SystemExit。纯本地。
"""
import contextlib
import importlib.util
import io
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (str(KE), str(ROOT / "skill" / "_common" / "opt"), str(ROOT)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

KL_RAW = 25.21          # kl 胞（Lz）上的原始 κ_L（W/mK）
C_KE = 34.879           # ke 胞高（Å）


def _g13():
    spec = importlib.util.spec_from_file_location("g13_0062", str(KE / "step8.3_output" / "gen_step13_output.py"))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


def _material(tmp, dim="2d", lz=25.0, c=C_KE):
    """<tmp>/ke（材料目录，带 S8.4 的 2d_correction.json）+ 同级 kl-dft-cpu 的 kappa_summary.json。"""
    m = Path(tmp) / "ke"
    (m / "step1_opt").mkdir(parents=True)
    (m / "step1_opt" / "workflow_method.txt").write_text("DIM=%s\n" % dim)
    (m / "step8.4_amset2d").mkdir()
    (m / "step8.4_amset2d" / "2d_correction.json").write_text(json.dumps({"cell_c_A": c}))
    kl = Path(tmp) / "kl-dft-cpu" / "step6_kappa"
    kl.mkdir(parents=True)
    summ = {"temperatures": [200.0, 400.0], "kappa_xx_yy_zz": [[KL_RAW, KL_RAW, 0.0]] * 2}
    if lz is not None:
        summ["Lz_ang"] = lz
    (kl / "kappa_summary.json").write_text(json.dumps(summ))
    return m


def test_2d_cell_mismatch_rescaled_by_lz_over_c(tmp_path):
    g = _g13()
    kxx, kyy, lines = g.resolve_kappa_L(_material(tmp_path))
    want = KL_RAW * 25.0 / C_KE
    assert abs(kxx - want) < 1e-9 and abs(kyy - want) < 1e-9
    assert any("[换算]" in ln for ln in lines)


def test_2d_rescale_keeps_sheet_conductance(tmp_path):
    """换算后 κ_L×c（ke 胞）= κ_raw×Lz（kl 胞）：片热导与真空无关。"""
    g = _g13()
    kxx, _, _ = g.resolve_kappa_L(_material(tmp_path))
    assert abs(kxx * C_KE - KL_RAW * 25.0) < 1e-9


def test_2d_zt_uses_rescaled_kappa(tmp_path):
    g = _g13()
    m = _material(tmp_path)
    row = {"amset_S_uV/K": -168.81, "amset_sigma_S/m": 9.733e4 / 5.190, "amset_kappa_e_W/mK": 0.6358 / 5.190}
    for d in ("xx", "yy"):
        row["amset_S_%s_uV/K" % d] = row["amset_S_uV/K"]
        row["amset_sigma_%s_S/m" % d] = row["amset_sigma_S/m"]
        row["amset_kappa_e_%s_W/mK" % d] = row["amset_kappa_e_W/mK"]
    old = os.getcwd()
    os.chdir(m)
    try:
        g._add_zt_columns(row)
    finally:
        os.chdir(old)
    kl = KL_RAW * 25.0 / C_KE
    want = (168.81e-6) ** 2 * row["amset_sigma_S/m"] * 300.0 / (row["amset_kappa_e_W/mK"] + kl)
    assert abs(row["amset_ZT"] - want) < 1e-12


def test_2d_same_cell_unchanged(tmp_path):
    g = _g13()
    kxx, _, lines = g.resolve_kappa_L(_material(tmp_path, lz=C_KE * 1.005))   # 差 0.5% < 1%
    assert kxx == KL_RAW
    assert not any("[换算]" in ln for ln in lines)


def test_3d_mismatch_only_warns(tmp_path):
    g = _g13()
    kxx, _, lines = g.resolve_kappa_L(_material(tmp_path, dim="3d"))
    assert kxx == KL_RAW
    assert not any("[换算]" in ln for ln in lines)
    assert any("CELL_POLICY" in ln for ln in lines)


def test_no_lz_metadata_unchanged(tmp_path):
    g = _g13()
    kxx, _, lines = g.resolve_kappa_L(_material(tmp_path, lz=None))
    assert kxx == KL_RAW
    assert not any("[换算]" in ln for ln in lines)


# ---------------------------------------------------------------- S8.1（BoltzTraP）
T_KE = 6.7201           # ke 的层厚（Å）


def _g11():
    spec = importlib.util.spec_from_file_location("g11_0062", str(KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py"))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


def _ninfo(norm="thickness", c=C_KE):
    f = c / T_KE if (c and norm == "thickness") else 1.0
    return {"norm": norm, "cell_c_A": c, "thickness_A": T_KE if c else None, "factor_applied": f}


def test_s81_2d_cell_mismatch_rescaled_then_thickness_factor(tmp_path):
    """以前 SystemExit；现在 κ_raw × Lz/c × (c/t) = κ_raw × Lz/t（按 ke 的层厚归一）。"""
    b = _g11()
    kxx, kyy, lines = b._resolve_kappa_L(_material(tmp_path), 300.0, _ninfo())
    want = KL_RAW * 25.0 / T_KE
    assert abs(kxx - want) < 1e-9 and abs(kyy - want) < 1e-9
    assert any("[换算]" in ln for ln in lines)


def test_s81_2d_same_cell_unchanged(tmp_path):
    b = _g11()
    kxx, _, lines = b._resolve_kappa_L(_material(tmp_path, lz=C_KE * 1.005), 300.0, _ninfo())
    assert abs(kxx - KL_RAW * C_KE / T_KE) < 1e-9
    assert any("元胞一致" in ln for ln in lines)
    assert not any("[换算]" in ln for ln in lines)


def test_s81_3d_without_cell_height_skips_gate(tmp_path):
    """3D：_norm_2d 不给 cell_c_A -> 闸门跳过，κ_L 原样（因子 1）。"""
    b = _g11()
    kxx, _, lines = b._resolve_kappa_L(_material(tmp_path, dim="3d"), 300.0, _ninfo(norm="cell(3D)", c=None))
    assert kxx == KL_RAW
    assert any("跳过元胞核对" in ln for ln in lines)


# ---------------------------------------------------------------- [0062-r2] step.conf 的 KL_ROOT
def _material_far_kl(tmp, kl_conf=None, step="step8.3_output", lz=25.0):
    """ke 在 <tmp>/MoS2/ke，kl 在另一个材料 <tmp>/MoS2_kl_v2/kl-dft-cpu（不是 ke 的同级目录）。"""
    m = Path(tmp) / "MoS2" / "ke"
    (m / "step1_opt").mkdir(parents=True)
    (m / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\n")
    (m / "step8.4_amset2d").mkdir()
    (m / "step8.4_amset2d" / "2d_correction.json").write_text(json.dumps({"cell_c_A": C_KE}))
    kl = Path(tmp) / "MoS2_kl_v2" / "kl-dft-cpu" / "step6_kappa"
    kl.mkdir(parents=True)
    (kl / "kappa_summary.json").write_text(json.dumps(
        {"temperatures": [200.0, 400.0], "kappa_xx_yy_zz": [[KL_RAW, KL_RAW, 0.0]] * 2, "Lz_ang": lz}))
    if kl_conf is not None:
        (m / "step.conf").write_text("[params]\nSTEP = %s\nKL_ROOT = %s\n" % (step, kl_conf))
    return m


def test_s83_kl_root_from_step_conf_relative(tmp_path):
    g = _g13()
    m = _material_far_kl(tmp_path, "../../MoS2_kl_v2/kl-dft-cpu")
    assert g.resolve_kappa_L(m)[0] is None                    # 不读 step.conf：同级找不到
    g = _g13()
    g._apply_conf(m)
    kxx, _, lines = g.resolve_kappa_L(m)
    assert abs(kxx - KL_RAW * 25.0 / C_KE) < 1e-9             # 找到了，并且照常 ×Lz/c
    assert any("MoS2_kl_v2" in ln or "step6_kappa" in ln for ln in lines)


def test_s83_kl_root_absolute(tmp_path):
    g = _g13()
    m = _material_far_kl(tmp_path, str(tmp_path / "MoS2_kl_v2" / "kl-dft-cpu"))
    g._apply_conf(m)
    assert g.resolve_kappa_L(m)[0] is not None


def test_s83_no_step_conf_keeps_default(tmp_path):
    g = _g13()
    m = _material_far_kl(tmp_path)
    g._apply_conf(m)                                          # 没有 step.conf：不动
    assert g.KL_ROOT is None


def test_s81_kl_root_from_step_conf(tmp_path):
    b = _g11()
    m = _material_far_kl(tmp_path, "../../MoS2_kl_v2/kl-dft-cpu", step="step8.1_boltztrap")
    b._apply_conf(m)
    kxx, _, lines = b._resolve_kappa_L(m, 300.0, _ninfo())
    assert kxx is not None and abs(kxx - KL_RAW * 25.0 / T_KE) < 1e-9
    assert any("[换算]" in ln for ln in lines)


def test_kl_root_in_common_step_conf_does_not_break_strict_steps(tmp_path):
    """[0063-r1] KL_ROOT / KL_KAPPA_DIRS 写在材料共用 step.conf 时会合并进每一步；
    S1 走 relax_common 的严格模式，必须放行这两个键，否则 gen 直接 SystemExit。"""
    sys.path.insert(0, str(ROOT / "skill" / "_common" / "opt"))
    import relax_common as R
    (tmp_path / "step.conf").write_text(
        "[params]\nSTEP = step1_opt\nKL_ROOT = /x/MoS2_kl_v2/kl-dft-cpu\n"
        "KL_KAPPA_DIRS = step6_kappa\n", encoding="utf-8")
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        assert R.load_step_params("step1_opt") is True
    finally:
        os.chdir(old)
