"""[0014] 2D 口径与代码一致：S8 弹性不乘 c/t、zt 对 2D 只认 S8.4、2D 自动打开 amset2d、
共享 κL 引擎 NAC=auto 对 2D 关。纯本地，不需要 VASP / AMSET / phono3py。"""
import contextlib
import importlib.util
import io
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from autozt import workflow as W  # noqa: E402

MOS2 = ("MoS2\n1.0\n3.14 0 0\n-1.57 2.719330 0\n0 0 25\nMo S\n1 2\nDirect\n"
        "0.333333 0.666667 0.5\n0.666667 0.333333 0.5626\n0.666667 0.333333 0.4374\n")
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


# ---------------------------------------------------------------- S8：2D 弹性不乘 c/t
def test_s8_2d_elastic_is_raw_slab(tmp_path):
    ke = ROOT / "skill" / "ke-dft-cpu"
    g = _load(ke / "step8_amset" / "gen_step10_amset.py", "g10_0014",
              (ke, ROOT / "skill" / "_common" / "opt", ROOT / "skill" / "_common",
               ke / "step8.4_amset2d"))
    (tmp_path / "step1_opt").mkdir()
    (tmp_path / "step1_opt" / "workflow_method.txt").write_text("DIM=2d\n")
    (tmp_path / "step1_opt" / "CONTCAR").write_text(MOS2)
    (tmp_path / "step8_amset").mkdir()
    el = [[130.0, 30.0, 1.0, 0, 0, 0], [30.0, 130.0, 1.0, 0, 0, 0], [1.0, 1.0, 2.0, 0, 0, 0],
          [0, 0, 0, 0.5, 0, 0], [0, 0, 0, 0, 0.5, 0], [0, 0, 0, 0, 0, 50.0]]
    with contextlib.redirect_stdout(io.StringIO()):
        is_2d, used, c_len = g.apply_2d_corrections(tmp_path, el)
    assert is_2d and abs(c_len - 25.0) < 1e-6
    assert used == el                                  # 原始 slab 值：不乘 c/t
    rec = json.loads((tmp_path / "step8_amset" / "2d_correction.json").read_text())
    assert rec["elastic_rescale_factor_c_over_t"] > 3.0   # c/t 照算（介电扣真空/厚度归一化用）
    assert rec["elastic_constant_used_GPa"] == el and rec["elastic_constant_rescaled_GPa"] is None


# ---------------------------------------------------------------- zt：2D 只认 S8.4
def _zt_spec():
    sk = yaml.safe_load((ROOT / "skill" / "zt-dft-cpu" / "skill.yaml").read_text(encoding="utf-8"))
    return sk["steps"][0]["needs_results"]


def _fetched(lp, skill, step, fname, text="{}"):
    m = {"name": lp.name, "qualified_name": "p/" + lp.name, "tt": skill,
         "lpath": str(lp), "result_dir": str(lp / skill / "result")}
    d = lp / skill / "result" / step
    d.mkdir(parents=True, exist_ok=True)
    (d / fname).write_text(text)
    W.write_result_origin(m, {"name": step, "label": step})


def _mat(tmp_path, dim):
    lp = tmp_path / "proj" / "M"
    lp.mkdir(parents=True)
    (lp / "POSCAR").write_text(MOS2 if dim == "2d" else SI)
    _fetched(lp, "ke-dft-cpu", "step1_opt", "workflow_method.txt", "DIM=%s\n" % dim)
    _fetched(lp, "kl-dft-cpu", "step6_kappa", "kappa_summary.json")
    return lp, {"name": "M", "lpath": str(lp)}


def test_zt_2d_waits_for_s84_not_s8(tmp_path):
    lp, m = _mat(tmp_path, "2d")
    _fetched(lp, "ke-dft-cpu", "step8_amset", "transport.json")   # S8 先跑完
    chosen, miss = W.resolve_needs_results(m, _zt_spec())
    assert "transport" not in chosen
    assert any("step8.4_amset2d" in x and "amset2d: true" in x and "DIM=2d" in x for x in miss), miss
    _fetched(lp, "ke-dft-cpu", "step8.4_amset2d", "transport.json")
    chosen, miss = W.resolve_needs_results(m, _zt_spec())
    assert miss == [] and chosen["transport"][1].endswith("step8.4_amset2d/transport.json")


def test_zt_3d_uses_s8(tmp_path):
    lp, m = _mat(tmp_path, "3d")
    _fetched(lp, "ke-dft-cpu", "step8_amset", "transport.json")
    chosen, miss = W.resolve_needs_results(m, _zt_spec())
    assert miss == [] and chosen["transport"][1].endswith("step8_amset/transport.json")


def test_plain_string_candidates_unchanged(tmp_path):
    lp, m = _mat(tmp_path, "2d")
    _fetched(lp, "ke-dft-cpu", "step8_amset", "transport.json")
    spec = {"transport": ["ke-dft-cpu/result/step8_amset/transport.json"]}
    chosen, miss = W.resolve_needs_results(m, spec)      # 不带 dim 的老写法不受维度影响
    assert miss == [] and "transport" in chosen


def test_s20_records_transport_route():
    g = _load(ROOT / "skill" / "zt-dft-cpu" / "step20_zt" / "gen_step20_zt.py", "g20_0014",
              (ROOT / "skill" / "_common" / "opt",))
    assert g.transport_route({"from": "ke-dft-cpu/result/step8.4_amset2d/transport.json"}, True) \
        == ("S8.4_amset2d", None)
    assert g.transport_route({"from": "ke-dft-cpu/result/step8_amset/transport.json"}, False) \
        == ("S8_kappa", None)
    step, warn = g.transport_route({"from": "ke-dft-cpu/result/step8_amset/transport.json"}, True)
    assert step == "S8_kappa" and "三维散射核" in warn


# ---------------------------------------------------------------- init：2D 自动 amset2d
def test_init_auto_dim_enables_amset2d(tmp_path, monkeypatch):
    from autozt import ops, report

    class _Dim:
        @staticmethod
        def detect_dimension(pos):
            return ("2d" if "MoS2" in Path(pos).read_text() else "3d", 2, 15.0)

    monkeypatch.setattr(report, "_dim_mod", lambda cfg, t: _Dim)
    sk = yaml.safe_load((ROOT / "skill" / "ke-dft-cpu" / "skill.yaml").read_text(encoding="utf-8"))
    assert sk["optional_steps"]["amset2d"]["auto_dim"] == "2d"
    t = {"key": "ke-dft-cpu", "optional_steps": sk["optional_steps"]}
    content = "task_types:\n  ke-dft-cpu:\n    local_root: \"..\"\n"
    (tmp_path / "POSCAR").write_text(MOS2)
    with contextlib.redirect_stdout(io.StringIO()):
        out = ops._auto_dim_flags({}, t, content, "ke-dft-cpu", str(tmp_path))
    seg = yaml.safe_load(out)["task_types"]["ke-dft-cpu"]
    assert seg["amset2d"] is True and seg["local_root"] == ".."
    (tmp_path / "POSCAR").write_text(SI)
    assert ops._auto_dim_flags({}, t, content, "ke-dft-cpu", str(tmp_path)) == content


# ---------------------------------------------------------------- κL：NAC=auto 对 2D 关
def test_fcfit_nac_auto_off_for_2d(tmp_path):
    g2 = _load(ROOT / "skill" / "_common" / "fcfit" / "gen_step2_kappa.py", "g2_0014",
               (ROOT / "skill" / "_common", ROOT / "skill" / "_common" / "opt"))
    fit = tmp_path / "fit"
    fit.mkdir()
    (fit / "POSCAR").write_text(SI)
    assert g2.resolve_nac("auto", fit) is False          # 没有 BORN
    (fit / "BORN").write_text("x\n")
    assert g2.resolve_nac("auto", fit) is True           # 3D + BORN
    (fit / "POSCAR").write_text(MOS2)
    with contextlib.redirect_stdout(io.StringIO()):
        assert g2.resolve_nac("auto", fit) is False      # 2D + BORN：不加 3D-NAC
    assert g2.resolve_nac("on", fit) is True             # 显式 on 仍可强制
    (fit / "phonon_summary.json").write_text(json.dumps({"is_2d": False}))
    assert g2.resolve_nac("auto", fit) is True           # phonon_summary 的 is_2d 优先
