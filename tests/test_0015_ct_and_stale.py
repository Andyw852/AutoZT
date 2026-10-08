"""[0015] S8.1 的 c/t 不再静默退回元胞口径；S8 的 c/t 变了让 S8.1 重排；
组合步骤（S20_zt）的跨技能输入变了 -> 判过期重新生成。纯本地。"""
import contextlib
import importlib.util
import io
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from autozt import workflow as W  # noqa: E402

KE = ROOT / "skill" / "ke-dft-cpu"
MOS2 = ("MoS2\n1.0\n3.14 0 0\n-1.57 2.719330 0\n0 0 25\nMo S\n1 2\nDirect\n"
        "0.333333 0.666667 0.5\n0.666667 0.333333 0.5626\n0.666667 0.333333 0.4374\n")


def _load(path, name, extra=()):
    for d in (str(Path(path).parent),) + tuple(str(x) for x in extra):
        if d not in sys.path:
            sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


def _bt2():
    return _load(KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py", "g11_0015",
                 (KE, ROOT / "skill" / "_common", ROOT / "skill" / "_common" / "opt"))


# ---------------------------------------------------------------- S8.1 的 c/t
def test_s81_ct_prefers_s8_record(tmp_path):
    g = _bt2()
    (tmp_path / "step8_amset").mkdir()
    (tmp_path / "step8_amset" / "2d_correction.json").write_text(json.dumps(
        {"cell_c_A": 20.0, "elastic_rescale_factor_c_over_t": 2.9718,
         "layer_thickness": {"thickness_used_A": 6.73}}))
    i = g._ct_info(tmp_path)
    assert i["c_over_t"] == 2.9718 and i["thickness_A"] == 6.73
    assert i["source"] == "step8_amset/2d_correction.json"


def test_s81_ct_computed_when_s8_not_generated(tmp_path):
    g = _bt2()
    (tmp_path / "step3_uniform").mkdir()
    (tmp_path / "step3_uniform" / "POSCAR").write_text(MOS2)
    with contextlib.redirect_stdout(io.StringIO()):
        i = g._ct_info(tmp_path)
        fac, info = g._norm_2d(tmp_path, True)
    assert i["source"].startswith("vdw") and abs(i["thickness_A"] - 6.73) < 0.02
    assert abs(fac - 25.0 / i["thickness_A"]) < 1e-6 and info["norm"] == "thickness"


def test_s81_no_silent_cell_fallback(tmp_path):
    g = _bt2()
    with pytest.raises(SystemExit):
        g._norm_2d(tmp_path, True)                 # 2D 拿不到 c/t：报错，不再 cell(回退)
    fac, info = g._norm_2d(tmp_path, False)        # 3D 不需要 c/t
    assert fac == 1.0 and info["norm"] == "cell(3D)"


# ---------------------------------------------------------------- S8 写出 c/t 后让 S8.1 重排
def _kc():
    return _load(KE / "ke_common.py", "kc_0015", (ROOT / "skill" / "_common" / "opt",))


def _s81(cwd, **rec):
    d = cwd / "step8.1_boltztrap"
    d.mkdir(exist_ok=True)
    f = d / "boltztrap_crta.json"
    f.write_text(json.dumps(dict({"dim": "2d"}, **rec)))
    return f


def test_ct_consumers_invalidated_only_when_ct_changes(tmp_path):
    kc = _kc()
    f = _s81(tmp_path, twoD_c_over_t=2.9718, twoD_c_over_t_source="step8_amset/2d_correction.json")
    with contextlib.redirect_stdout(io.StringIO()):
        assert kc.invalidate_ct_consumers(tmp_path, 2.9719) == []        # 一样：不动
        assert f.is_file()
        done = kc.invalidate_ct_consumers(tmp_path, 2.8796)              # 变了：归档
    assert ("step8.1_boltztrap", "boltztrap_crta.json") in done and not f.is_file()


def test_ct_consumers_untouched_for_manual_thickness_or_3d(tmp_path):
    kc = _kc()
    f = _s81(tmp_path, twoD_c_over_t=3.5, twoD_c_over_t_source="THICKNESS_A")
    assert kc.invalidate_ct_consumers(tmp_path, 2.97) == [] and f.is_file()
    f = _s81(tmp_path, dim="3d")
    assert kc.invalidate_ct_consumers(tmp_path, 2.97) == [] and f.is_file()


def test_ct_consumers_legacy_fallback_record_is_invalidated(tmp_path):
    kc = _kc()
    f = _s81(tmp_path, twoD_c_over_t=None)          # 旧版本：S8.1 先跑、退回了元胞口径
    with contextlib.redirect_stdout(io.StringIO()):
        assert kc.invalidate_ct_consumers(tmp_path, 3.06)
    assert not f.is_file()


# ---------------------------------------------------------------- S20：跨技能输入变了 -> 过期
SPEC = {"transport": [{"path": "ke-dft-cpu/result/step8.4_amset2d/transport.json", "dim": "2d"},
                      {"path": "ke-dft-cpu/result/step8_amset/transport.json", "not_dim": "2d"}],
        "kappa_L": ["kl-dft-cpu/result/step6_kappa/kappa_summary.json"],
        "structure": ["ke-dft-cpu/result/step1_opt/workflow_method.txt"]}


def _fetched(lp, skill, step, fname, text="{}"):
    m = {"name": lp.name, "qualified_name": "p/" + lp.name, "tt": skill,
         "lpath": str(lp), "result_dir": str(lp / skill / "result")}
    d = lp / skill / "result" / step
    d.mkdir(parents=True, exist_ok=True)
    (d / fname).write_text(text)
    W.write_result_origin(m, {"name": step, "label": step})


def _zt(tmp_path):
    lp = tmp_path / "p" / "M"
    lp.mkdir(parents=True)
    (lp / "POSCAR").write_text(MOS2)
    _fetched(lp, "ke-dft-cpu", "step1_opt", "workflow_method.txt", "DIM=3d\n")
    _fetched(lp, "ke-dft-cpu", "step8_amset", "transport.json", '{"v": 1}')
    _fetched(lp, "kl-dft-cpu", "step6_kappa", "kappa_summary.json", '{"k": 1}')
    t = {"key": "zt-dft-cpu", "steps": [{"name": "step20_zt", "label": "S20_zt",
                                         "needs": [], "needs_results": SPEC}]}
    m = {"name": "M", "lpath": str(lp), "tt": "zt-dft-cpu",
         "result_dir": str(lp / "zt-dft-cpu" / "result")}
    return lp, t, m


def _done_step():
    return [{"name": "step20_zt", "label": "S20_zt", "job": None, "done": True,
             "has_outcar": False, "has_slurm_out": False, "has_incar": False}]


def test_s20_stale_when_upstream_result_changes(tmp_path):
    lp, t, m = _zt(tmp_path)
    chosen, miss = W.resolve_needs_results(m, SPEC)
    assert miss == []
    W.write_inputs_record(m, "step20_zt", W.needs_results_fingerprint(chosen))   # 模拟 gen 成功

    m["steps"] = _done_step()
    W._dag_recompute(t, m)
    assert m["steps"][0]["kind"] == "OK"

    _fetched(lp, "ke-dft-cpu", "step8_amset", "transport.json", '{"v": 1}')      # 重新 fetch 同一份
    m["steps"] = _done_step()
    W._dag_recompute(t, m)
    assert m["steps"][0]["kind"] == "OK"                                         # 只换来源戳：不重算

    _fetched(lp, "kl-dft-cpu", "step6_kappa", "kappa_summary.json", '{"k": 2}')  # κL 重算了
    m["steps"] = _done_step()
    ready = W._dag_recompute(t, m)
    s = m["steps"][0]
    assert s["label_txt"] == "STALE" and s in ready and "kappa_L" in s["diag"]


def test_s20_stale_when_route_switches(tmp_path):
    lp, t, m = _zt(tmp_path)
    chosen, _ = W.resolve_needs_results(m, SPEC)
    W.write_inputs_record(m, "step20_zt", W.needs_results_fingerprint(chosen))
    # 结构其实是 2D（S1 重新判维）：电子段改选 S8.4，旧 S20 用的是 S8
    _fetched(lp, "ke-dft-cpu", "step1_opt", "workflow_method.txt", "DIM=2d\n")
    _fetched(lp, "ke-dft-cpu", "step8.4_amset2d", "transport.json", '{"v": 2}')
    m["steps"] = _done_step()
    W._dag_recompute(t, m)
    s = m["steps"][0]
    assert s["label_txt"] == "STALE" and "step8.4_amset2d" in s["diag"]


def test_s20_without_record_is_not_touched(tmp_path):
    lp, t, m = _zt(tmp_path)                      # 0015 之前生成的 S20：没有记录，不判过期
    m["steps"] = _done_step()
    W._dag_recompute(t, m)
    assert m["steps"][0]["kind"] == "OK"
