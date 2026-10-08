"""[0016] 2D 新材料默认不跑 S8（amset3d 可选组 + auto_dim_off）；层厚走共用 step.conf；
S8.3 软依赖；σh / ZA 自动判断。纯本地，不需要 VASP / AMSET。"""
import contextlib
import copy
import importlib.util
import io
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from autozt import workflow as W  # noqa: E402

KE = ROOT / "skill" / "ke-dft-cpu"
OPT = ROOT / "skill" / "_common" / "opt"
H = "3.18 0 0\n-1.59 2.754 0\n0 0 25\n"
MOS2_2H = ("MoS2\n1.0\n" + H + "Mo S\n1 2\nDirect\n0.333333 0.666667 0.5\n"
           "0.666667 0.333333 0.5626\n0.666667 0.333333 0.4374\n")
TIS2_1T = ("TiS2\n1.0\n" + H + "Ti S\n1 2\nDirect\n0 0 0.5\n"
           "0.333333 0.666667 0.5626\n0.666667 0.333333 0.4374\n")
JANUS = ("MoSSe\n1.0\n" + H + "Mo S Se\n1 1 1\nDirect\n0.333333 0.666667 0.5\n"
         "0.666667 0.333333 0.5626\n0.666667 0.333333 0.4300\n")
SI = ("Si\n5.43\n0.0 0.5 0.5\n0.5 0.0 0.5\n0.5 0.5 0.0\nSi\n2\nDirect\n"
      "0.0 0.0 0.0\n0.25 0.25 0.25\n")


def _load(path, name, extra=()):
    for d in (str(Path(path).parent),) + tuple(str(x) for x in extra):
        if d not in sys.path:
            sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


def _ke():
    from autozt import load_config, apply_skills
    cfg, _ = load_config(None)
    return apply_skills(cfg)["task_types"]["ke-dft-cpu"]


def _expand(t0, flags):
    from autozt.bootstrap import expand_optional_steps
    t = copy.deepcopy(t0)
    t.update(flags)
    expand_optional_steps(t)
    return [s["name"] for s in t["steps"]]


# ---------------------------------------------------------------- S8 进可选组
def test_ke_default_steps_unchanged_and_2d_drops_only_s8():
    t0 = _ke()
    default = _expand(t0, {})
    assert "step8_amset" in default
    s8x = [n for n in default if n.startswith("step8")]
    assert s8x == ["step8_amset", "step8.1_boltztrap", "step8.2_dpt", "step8.3_output"]
    two_d = _expand(t0, {"amset3d": False, "amset2d": True})
    assert [n for n in two_d if n.startswith("step8")] == [
        "step8.1_boltztrap", "step8.2_dpt", "step8.3_output", "step8.4_amset2d"]
    assert set(_expand(t0, {"amset2d": True})) - set(two_d) == {"step8_amset"}


def test_init_writes_amset3d_false_for_2d(tmp_path, monkeypatch):
    from autozt import ops, report

    class _Dim:
        @staticmethod
        def detect_dimension(pos):
            return ("2d" if "MoS2" in Path(pos).read_text() else "3d", 2, 15.0)

    monkeypatch.setattr(report, "_dim_mod", lambda cfg, t: _Dim)
    sk = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    assert sk["optional_steps"]["amset3d"]["auto_dim_off"] == "2d"
    t = {"key": "ke-dft-cpu", "optional_steps": sk["optional_steps"]}
    content = "task_types:\n  ke-dft-cpu:\n    local_root: \"..\"\n"
    (tmp_path / "POSCAR").write_text(MOS2_2H)
    with contextlib.redirect_stdout(io.StringIO()):
        seg = yaml.safe_load(ops._auto_dim_flags({}, t, content, "ke-dft-cpu",
                                                 str(tmp_path)))["task_types"]["ke-dft-cpu"]
    assert seg["amset2d"] is True and seg["amset3d"] is False
    (tmp_path / "POSCAR").write_text(SI)
    assert ops._auto_dim_flags({}, t, content, "ke-dft-cpu", str(tmp_path)) == content


# ---------------------------------------------------------------- 软依赖
def _step(name, done=False):
    return {"name": name, "label": name, "job": None, "done": done,
            "has_outcar": False, "has_slurm_out": False, "has_incar": False}


def test_soft_dependency_waits_when_present_and_is_silent_when_absent(tmp_path):
    t = {"key": "x", "steps": [{"name": "a", "needs": []}, {"name": "b", "needs": []},
                               {"name": "c", "needs": ["a?", "b"]}]}
    m = {"name": "M", "lpath": str(tmp_path), "steps": [_step("b", True), _step("c")]}
    W._dag_recompute(t, m)
    c = m["steps"][-1]
    assert c["_missing_deps"] == [] and c["kind"] == "PREP"        # a 不在：不提示缺、不卡
    m["steps"] = [_step("a"), _step("b", True), _step("c")]
    W._dag_recompute(t, m)
    assert m["steps"][-1]["kind"] == "WAIT"                       # a 在、没完成：等它


# ---------------------------------------------------------------- S8.1 层厚
def _bt2():
    return _load(KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py", "g11_0016",
                 (KE, ROOT / "skill" / "_common", OPT))


def test_s81_reads_shared_layer_thickness_from_step_conf(tmp_path):
    g = _bt2()
    (tmp_path / "step3_uniform").mkdir()
    (tmp_path / "step3_uniform" / "POSCAR").write_text(MOS2_2H)
    (tmp_path / "step.conf").write_text("[params]\nLAYER_THICKNESS = 6.73\nFUNC = pbesol\n")
    with contextlib.redirect_stdout(io.StringIO()):
        g._apply_conf(tmp_path)
        i = g._ct_info(tmp_path)
    assert i["source"] == "LAYER_THICKNESS(step.conf)" and abs(i["thickness_A"] - 6.73) < 1e-9
    assert abs(i["c_over_t"] - 25.0 / 6.73) < 1e-6


def test_s81_falls_back_to_s84_record_without_s8(tmp_path):
    g = _bt2()
    (tmp_path / "step8.4_amset2d").mkdir()
    (tmp_path / "step8.4_amset2d" / "2d_correction.json").write_text(json.dumps(
        {"cell_c_A": 20.0, "elastic_rescale_factor_c_over_t": 3.0,
         "layer_thickness": {"thickness_used_A": 6.6667}}))
    i = g._ct_info(tmp_path)
    assert i["source"] == "step8.4_amset2d/2d_correction.json" and i["c_over_t"] == 3.0


def test_layer_thickness_is_a_reserved_shared_key():
    stepconf = _load(OPT / "stepconf.py", "stepconf_0016")      # 按路径加载：别的测试可能先缓存了另一份
    assert "LAYER_THICKNESS" in stepconf.RESERVED_PARAMS
    merged = stepconf.parse("[params]\nLAYER_THICKNESS = 6.73\n")
    stepconf.StepConf(merged, {"FUNC": ("pbesol", "str")}, strict=True)   # 严格模式不报错


# ---------------------------------------------------------------- c/t 失效的来源优先级
def _kc():
    return _load(KE / "ke_common.py", "kc_0016", (OPT,))


def test_s84_record_does_not_invalidate_s81_that_used_s8(tmp_path):
    kc = _kc()
    d = tmp_path / "step8.1_boltztrap"
    d.mkdir()
    f = d / "boltztrap_crta.json"
    f.write_text(json.dumps({"dim": "2d", "twoD_c_over_t": 2.9718,
                             "twoD_c_over_t_source": "step8_amset/2d_correction.json"}))
    with contextlib.redirect_stdout(io.StringIO()):
        assert kc.invalidate_ct_consumers(tmp_path, 2.87, step="step8.4_amset2d") == []
        assert f.is_file()                                         # S8 的记录优先，S8.4 不该动它
        assert kc.invalidate_ct_consumers(tmp_path, 2.87, step="step8_amset")
    assert not f.is_file()


# ---------------------------------------------------------------- σh / ZA
def test_za_coupling_check_by_horizontal_mirror(tmp_path):
    kc = _kc()
    for name, txt, want in (("2H", MOS2_2H, True), ("1T", TIS2_1T, False), ("Janus", JANUS, False)):
        d = tmp_path / name / "step3_uniform"
        d.mkdir(parents=True)
        (d / "POSCAR").write_text(txt)
        r = kc.za_coupling_check(tmp_path / name)
        assert r["sigma_h"] is want, (name, r)
    assert kc.za_coupling_check(tmp_path / "missing")["sigma_h"] is None


def test_za_recorded_by_s82_and_s84():
    s82 = (KE / "step8.2_dpt" / "gen_step12_dpt.py").read_text(encoding="utf-8")
    s84 = (KE / "step8.4_amset2d" / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert 'res["za_coupling"] = _za_check(cwd)' in s82
    assert '"za_coupling": _za_check(cwd)' in s84
