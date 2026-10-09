"""[0021] ① S8.1 的 DPT τ 换到扫描温度（以前生产路径总读 dpt_result.json 的 300 K τ，600 K 扫描 σ/PF/κe 偏大 2 倍）；
② S8.3 的 DPT μ/τ 换到 TARGET_T；③ S3/S7 网格闸门查 2D 面内高对称点对齐（K 的 1/3、M/X 的 1/2）；
④ S8.3 不再对 2D 真实重叠（出厂默认）判红。纯本地。"""
import contextlib
import importlib.util
import io
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (ROOT, KE, ROOT / "skill" / "_common" / "opt"):
    if str(_d) not in sys.path:
        sys.path.insert(0, str(_d))
import ke_common as kc  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- ③ 网格对齐
def _mat(tmp, lat, mesh, style="Gamma"):
    s3 = tmp / "step3_uniform"
    s3.mkdir(parents=True)
    (s3 / "POSCAR").write_text("x\n1.0\n" + "".join("%f %f %f\n" % tuple(v) for v in lat)
                               + "Cr S\n1 2\nDirect\n0 0 0.5\n0.3333 0.6667 0.57\n0.3333 0.6667 0.43\n")
    (s3 / "KPOINTS").write_text("auto\n0\n%s\n %d %d %d\n" % ((style,) + tuple(mesh)))
    (s3 / kc.METHOD_FILE).write_text("DIM=2D\n")
    return tmp


def _hex(a=3.2, c=20.0):
    return [(a, 0, 0), (-a / 2, a * math.sqrt(3) / 2, 0), (0, 0, c)]


RECT_SS = [(10.82, 0, 0), (0, 3.12, 0), (0, 0, 20.0)]           # SS 的矩形超胞：K 折叠到 y = 1/3、2/3


def _levels(iss, key):
    return [lv for lv, m in iss if key in m]


def test_hex_k_off_grid_is_error(tmp_path):
    iss = kc.s3_grid_issues(_mat(tmp_path, _hex(), (47, 47, 3)))
    assert _levels(iss, "K = (1/3, 1/3)") == ["error"]


def test_hex_odd_multiple_of_3_only_warns(tmp_path):            # WSe2 / MoSe2 的 45×45×3：K 在、M 不在
    iss = kc.s3_grid_issues(_mat(tmp_path, _hex(), (45, 45, 3)))
    assert [lv for lv, _ in iss] == ["warn"] and "M = (1/2, 0)" in iss[0][1]


def test_hex_multiple_of_6_passes(tmp_path):
    assert kc.s3_grid_issues(_mat(tmp_path, _hex(), (48, 48, 3))) == []


def test_rect_supercell_y41_is_error(tmp_path):                 # SS/LS 旧 S3
    iss = kc.s3_grid_issues(_mat(tmp_path, RECT_SS, (12, 41, 3)))
    assert any(lv == "error" and "第 2 轴 41 分" in m and "1/3" in m for lv, m in iss)


def test_rect_supercell_y42_passes(tmp_path):
    assert kc.s3_grid_issues(_mat(tmp_path, RECT_SS, (12, 42, 3))) == []


def test_rect_one_missing_point_warns(tmp_path):
    iss = kc.s3_grid_issues(_mat(tmp_path, RECT_SS, (12, 45, 3)))
    assert [lv for lv, _ in iss] == ["warn"]


def test_monkhorst_pack_not_checked(tmp_path):
    iss = kc.s3_grid_issues(_mat(tmp_path, _hex(), (47, 47, 3), style="Monkhorst-Pack"))
    assert iss == []


def test_gate_blocks_and_mentions_rule(tmp_path):
    d = _mat(tmp_path, RECT_SS, (12, 41, 3))
    with contextlib.redirect_stdout(io.StringIO()), pytest.raises(SystemExit) as e:
        kc.s3_grid_gate(d, "S8.2", check_kz=False)
    assert "6 的倍数" in str(e.value.code)


def test_s7_47_vs_s3_48_is_error(tmp_path):                     # WS2 的 S7
    d = _mat(tmp_path, _hex(), (48, 48, 3))
    sub = d / "step7_deform" / "deform-01"
    sub.mkdir(parents=True)
    (sub / "KPOINTS").write_text("auto\n0\nGamma\n 47 47 3\n")
    iss = kc.s7_grid_issues(d)
    assert [lv for lv, _ in iss] == ["error"] and "K = (1/3, 1/3)" in iss[0][1]


# ---------------------------------------------------------------- ① S8.1 τ 换温
def _dpt_json(mat, T0=300.0, with_bd=True, is_2d=True):
    res = []
    for car in ("electron", "hole"):
        r = {"carrier": car, "mobility_cm2_Vs": 100.0, "inputs": {"m_eff_m0": 0.5}}
        if with_bd:
            r["by_direction"] = {"status": "ok", "x": {"tau_s": 1.0e-14, "mobility_cm2_Vs": 100.0},
                                 "y": {"tau_s": 2.0e-14, "mobility_cm2_Vs": 120.0}, "m_d_m0": 0.5}
        res.append(r)
    (mat / "step8.2_dpt").mkdir(parents=True, exist_ok=True)
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(
        {"temperature_K": T0, "is_2d": is_2d, "dim": "2d" if is_2d else "3d", "results": res}))
    return mat


def _bt():
    return _load("bt_0021", KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py")


def test_s81_tau_scaled_to_scan_T(tmp_path):
    B = _bt()
    _dpt_json(tmp_path)
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        t = B._tau_aniso(tmp_path, True, 600.0)
    assert t["electron"]["x"] == pytest.approx(0.5e-14) and t["hole"]["y"] == pytest.approx(1.0e-14)
    assert "换算到 600 K" in buf.getvalue()


def test_s81_tau_isotropic_fallback_and_3d(tmp_path):
    B = _bt()
    _dpt_json(tmp_path, with_bd=False, is_2d=False)
    with contextlib.redirect_stdout(io.StringIO()):
        t = B._tau_from_json(tmp_path, 1200.0, False)
    t300 = (100.0 * 1e-4) * (0.5 * 9.1093837015e-31) / 1.602176634e-19
    assert t["electron"]["x"] == pytest.approx(t300 / 8.0)        # 3D ∝ T^-1.5


def test_s81_json_preferred_over_module(tmp_path, monkeypatch):
    B = _bt()
    _dpt_json(tmp_path)
    monkeypatch.setattr(B, "_dpt_module", lambda: (_ for _ in ()).throw(AssertionError("不该调模块")))
    with contextlib.redirect_stdout(io.StringIO()):
        assert B._tau_aniso(tmp_path, True, 300.0)["electron"]["x"] == pytest.approx(1e-14)


# ---------------------------------------------------------------- ② S8.3 换到 TARGET_T
def test_s83_dpt_scaled_to_target_T(tmp_path):
    S = _load("s83_0021", KE / "step8.3_output" / "gen_step13_output.py")
    _dpt_json(tmp_path, T0=600.0)
    with contextlib.redirect_stdout(io.StringIO()):
        d = S.load_dpt(tmp_path)
    assert d["T_used"] == 300.0 and d["T0"] == 600.0
    assert d["electron"] == pytest.approx(2 * 110.0)              # (100+120)/2 × 600/300
    assert S._dpt_tau_dir(d, "electron", "x") == pytest.approx(2e-14)
    _dpt_json(tmp_path, T0=300.0)
    d = S.load_dpt(tmp_path)
    assert d["electron"] == pytest.approx(110.0) and S._dpt_tau_dir(d, "hole", "y") == pytest.approx(2e-14)


# ---------------------------------------------------------------- ④ 重叠提示
def test_s83_no_longer_demands_unity_rerun():
    src = (KE / "step8.3_output" / "gen_step13_output.py").read_text(encoding="utf-8")
    assert "必须用 unity_overlap 重跑本步" not in src
    assert "2D + 真实重叠（出厂默认" in src
