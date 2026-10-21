"""[0062] S8.3 自动取 κ_L：2D 两链胞高不同（kl Lz ≠ ke c）时按片热导守恒 ×Lz/c 换到 ke 元胞口径。

以前只打印"元胞不一致"告警、原样拿 kl 胞的原始 κ_L 去配 ke 胞的 σ/κ_e，κ_L 放大 c/Lz 倍、
ZT 偏低（MoS2：Lz 25.0 vs c 34.879）。纯本地。
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
