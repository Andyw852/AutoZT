# -*- coding: utf-8 -*-
"""[metal_mode] ke-dft-cpu 金属模式（纯本地、无数据）。

  · S7.1：undeformed 能带 is_metal -> 写 metal.json、退出码 3，外层判金属：不写 band_edges.json / deformation_vac.h5，
    deformation.h5 照常放进本步目录（S7.1 判完成，S8.1 不被卡住）；
  · S8.1：金属模式不过体系判别闸门、不加剪刀、不读 DPT τ；在电中性 E_F(T) 处给 S、σ/τ、κ_e/τ、L；
    τ 只能来自 METAL_LAMBDA_TR（高温极限 ħ/τ = 2π λ_tr k_B T）或 METAL_TAU_FS（假设值），不给就不出绝对值；
  · S8.2 / S8.3：遇到 metal.json 直接报"不适用"并给出关组办法。
"""
import ast
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KE = REPO / "skill" / "ke-dft-cpu"
for _p in (str(REPO), str(KE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
S71 = KE / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py"
S81 = KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py"
S82 = KE / "step8.2_dpt" / "gen_step12_dpt.py"
S83 = KE / "step8.3_output" / "gen_step13_output.py"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _metal_dir(tmp_path):
    d = tmp_path / "step7b_deform_read"
    d.mkdir()
    (d / "metal.json").write_text(json.dumps({"is_metal": True, "reason": "test"}))
    return tmp_path


# ------------------------------------------------------------------ ke_common
def test_metal_status_and_not_applicable(tmp_path):
    import ke_common as kc
    assert kc.metal_status(tmp_path) is None
    kc.metal_not_applicable(tmp_path, "S8.2_dpt", "x")            # 不是金属：什么都不做
    _metal_dir(tmp_path)
    assert kc.metal_status(tmp_path)["reason"] == "test"
    with pytest.raises(SystemExit) as e:
        kc.metal_not_applicable(tmp_path, "S8.2_dpt", "x")
    assert "dpt: false" in str(e.value) and "output_compare: false" in str(e.value)
    assert "metal.json" in kc.DONE_MARKERS["step7b_deform_read"]  # S7 重算时旧判定一并归档


# ------------------------------------------------------------------ S7.1
def _s71_src():
    return S71.read_text(encoding="utf-8")


def test_s71_metal_block_writes_marker_and_exits_3(tmp_path):
    s = _s71_src()
    i = s.index("_be_code = r'''") + len("_be_code = r'''")
    code = s[i:s.index("'''", i)]
    node = next(n for n in ast.walk(ast.parse(code))
                if isinstance(n, ast.If) and isinstance(n.test, ast.Call)
                and getattr(n.test.func, "id", "") == "_is_metal")

    class BS:
        efermi = 1.25

        def is_metal(self):
            return True
    ns = {"json": json, "sys": type("S", (), {"argv": ["x", str(tmp_path), "rev-1"], "exit": sys.exit}),
          "outdir": str(tmp_path), "bs": BS(), "_is_metal": lambda b: b.is_metal()}
    with pytest.raises(SystemExit) as e:
        exec(compile(ast.Module(body=[node], type_ignores=[]), "be_code", "exec"), ns)
    assert e.value.code == 3
    rec = json.loads((tmp_path / "metal.json").read_text())
    assert rec["is_metal"] is True and rec["efermi_eV"] == 1.25 and rec["skill_rev"] == "rev-1"


def test_s71_outer_handles_metal():
    s = _s71_src()
    assert "_metal = (_rc == 3 and _metal_json.is_file())" in s
    m = s[s.index("_metal = (_rc == 3"):]
    assert m.index("if _metal:") < m.index("elif _rc != 0 or not be_out.is_file():")
    # 金属不走真空口径（否则 2D 金属会被"真空口径 h5 生不出来就报错"拦下），旧 vac h5 清掉
    v = m[m.index("ref = _resolve_deform_ref(cwd)"):]
    assert v.index("if _metal:") < v.index('elif ref == "vacuum":')
    assert '(out / "deformation_vac.h5").unlink()' in v[:v.index('elif ref == "vacuum":')]
    # deformation.h5 照常挪进本步目录（完成标记）
    assert "os.replace(str(h5), str(dst))" in v
    # 每次重跑先删旧 metal.json
    assert "_metal_json.unlink()" in s


# ------------------------------------------------------------------ S8.1
@pytest.fixture()
def s81(monkeypatch):
    monkeypatch.setenv("AUTOZT_CONDA_SH", "/dev/null")
    return _load(S81, "gen11_metal")


def test_metal_tau_formula(s81):
    tau, why = s81.metal_tau_s(300.0, lambda_tr=1.0, lambda_source="EPW 测试")
    assert tau == pytest.approx(1.054571817e-34 / (2 * math.pi * 1.380649e-23 * 300.0))
    assert tau * 1e15 == pytest.approx(4.05, abs=0.01) and "EPW 测试" in why
    t2, _ = s81.metal_tau_s(600.0, lambda_tr=1.0, lambda_source="x")
    assert t2 == pytest.approx(tau / 2)                           # τ ∝ 1/T
    assert s81.metal_tau_s(300, tau_fs=10)[0] == pytest.approx(1e-14)
    assert s81.metal_tau_s(300)[0] is None
    with pytest.raises(SystemExit):
        s81.metal_tau_s(300, tau_fs=10, lambda_tr=1, lambda_source="x")   # 二选一
    with pytest.raises(SystemExit):
        s81.metal_tau_s(300, lambda_tr=0.5)                                # λ 必须写出处


def _fake_res():
    mu = [-0.2, -0.1, 0.0, 0.1, 0.2]
    # 电中性点在 μ=+0.05（carrier 从 + 到 -）
    carr = [[3e20, 2e20, 1e20, -1e20, -2e20], [3e20, 2e20, 1e20, -1e20, -2e20]]
    sig = [[1e19] * 5, [2e19] * 5]
    kap = [[1e19 * 2.44e-8 * 300] * 5, [2e19 * 2.44e-8 * 600] * 5]
    return {"temperatures_K": [300.0, 600.0], "mu_rel_efermi_eV": mu, "carrier_conc_cm-3": carr,
            "sigma_over_tau_S_per_m_s": sig, "kappa_e_over_tau_W_per_m_K_s": kap,
            "seebeck_V_per_K": [[-1e-6, -2e-6, -3e-6, -4e-6, -5e-6]] * 2}


def test_metal_transport_tau_free_and_absolute(s81):
    rows, why = s81.metal_transport(_fake_res())
    assert "只输出" in why
    r = rows[0]
    assert r["mu_rel_efermi_eV"] == pytest.approx(0.05)
    assert r["seebeck_uV_per_K"] == pytest.approx(-3.5)
    assert r["lorenz_WOhm_per_K2"] == pytest.approx(2.44e-8)
    assert r["tau_fs"] is None and "sigma_S_per_m" not in r            # 不给 τ 就不编绝对值
    rows, _ = s81.metal_transport(_fake_res(), tau_fs=10)
    assert rows[0]["sigma_S_per_m"] == pytest.approx(1e5)
    assert rows[0]["resistivity_uOhm_cm"] == pytest.approx(1e8 / 1e5)
    assert rows[1]["kappa_e_W_per_m_K"] == pytest.approx(2e19 * 2.44e-8 * 600 * 1e-14)


def test_s81_main_metal_branch_wiring():
    s = S81.read_text(encoding="utf-8")
    m = s[s.index("def main():"):]
    assert m.index("metal = _metal_status(cwd)") < m.index("_disc_gate()")
    assert "scissor_gap_eV=(None if metal else _read_target_gap(cwd))" in m
    assert m.index("_write_metal(out, cwd, res, dim)") < m.index("if PAPER_SCAN:")
    for k in ("METAL_TAU_FS", "METAL_LAMBDA_TR", "METAL_LAMBDA_SOURCE"):
        assert '"%s"' % k in s[s.index("SPEC = {"):s.index("}", s.index("SPEC = {"))]


# ------------------------------------------------------------------ S8.2 / S8.3
def test_s82_guard_before_work():
    s = S82.read_text(encoding="utf-8")
    m = s[s.index("def main():"):]
    assert m.index('_kc.metal_not_applicable(cwd, "S8.2_dpt"') < m.index("invalidate_consumers(cwd)")


def test_s83_refuses_metal(tmp_path):
    _metal_dir(tmp_path)
    r = subprocess.run([sys.executable, "-c",
                        "import runpy,sys; sys.argv=['g']; runpy.run_path(%r, run_name='__main__')" % str(S83)],
                       cwd=tmp_path, capture_output=True, text=True, timeout=120)
    if "No module named" in r.stderr and "S8.3_cmp 不适用于金属" not in r.stderr:
        pytest.skip("S8.3 依赖的包不在本环境：%s" % r.stderr.strip().splitlines()[-1])
    assert r.returncode != 0 and "S8.3_cmp 不适用于金属" in r.stderr
    assert "output_compare: false" in r.stderr
