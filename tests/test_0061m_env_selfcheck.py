"""作业内环境自检 + gen 端网格依赖降级的回归测试。

[2026-10-22] 缘起：yue 上 AMSET 作业因
  ① AMSET 环境缺 ase（pymatgen.io.ase.AseAtomsAdaptor 在算到中途才炸，白占节点 1:35）；
  ② gen 环境没有 BoltzTraP2（旧代码在 apply_mesh_min 里直接 ModuleNotFoundError 退出）。
两类都是"环境问题"，但都应在开头几秒测出来 / 优雅降级，而不是浪费排队。
"""
import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skill" / "ke-dft-cpu" / "step8.4_amset2d"))
sys.path.insert(0, str(ROOT / "skill" / "ke-dft-cpu" / "step8_amset"))


def _load_preflight():
    import overlap_preflight
    return overlap_preflight


def test_runtime_dep_problems_empty_when_all_present(monkeypatch):
    pf = _load_preflight()
    monkeypatch.setattr(pf, "_AMSET_RUNTIME_MODS", (("json", "x"),))
    assert pf._amset_runtime_dep_problems() == []


def test_runtime_dep_problems_reports_missing_with_hint(monkeypatch):
    pf = _load_preflight()
    monkeypatch.setattr(pf, "_AMSET_RUNTIME_MODS",
                        (("json", "x"), ("__definitely_missing_mod_0061__", "pip install it")))
    bad = pf._amset_runtime_dep_problems()
    assert len(bad) == 1
    mod, exc, hint = bad[0]
    assert mod == "__definitely_missing_mod_0061__"
    assert isinstance(exc, Exception)
    assert hint == "pip install it"


def test_job_env_check_flags_missing_ase_like_module(monkeypatch):
    """amset 版本对但依赖缺 -> job_env_check 必须 False 且带修法（不能只报版本 ok）。"""
    pf = _load_preflight()
    monkeypatch.setattr(pf, "REQUIRED_AMSET", "any")
    monkeypatch.setattr(pf, "_AMSET_RUNTIME_MODS",
                        (("__missing_ase_like_0061__", "pip install ase"),))

    class _FakeAmset:
        __version__ = "0.5.1"

    monkeypatch.setitem(sys.modules, "amset", _FakeAmset())
    ok, msg = pf.job_env_check()
    assert ok is False
    assert "运行期依赖不全" in msg and "pip install ase" in msg
