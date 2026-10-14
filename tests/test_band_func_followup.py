# -*- coding: utf-8 -*-
"""2026-10-14：能带/泛函设置修复（75c77e8）复审后的补丁回归测试。

  [A] FUNC=auto：rerun / 换集群后靠 autozt 推来的本地结果副本（.autozt_hints/）沿用 / 对齐；FUNC_REASON 留痕；
      空的 FUNC_ALIGN 当 true；step1*.bak 之类的备份目录不算
  [B] autozt push_paths 的 {from, to} 条目：改名落地、md5 相同不重写、本地没有就删远端那份
  [C] zt S20：kl_method.json 缺失时读 kl S1 的记录（只认 κ_L 来自 kl-dft-cpu）；kl_method.json 不是对象不崩
  [D] kl S5.1：缓存指纹含路径坐标与点数；seekpath 是可选依赖
  [E] ke：剪刀带隙来源标签读 band_summary 的 hybrid.label；hybrid_params 按 GGA 区分 HSE06 / HSEsol
（SOC 闸门 inconclusive、m* 同一组 k 点拟合的测试在 test_band_func_consistency.py）
"""
import json
import subprocess
import sys
import types
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_band_func_consistency import KE, OPT, ROOT, _load  # noqa: E402


# --------------------------------------------------------------------------- [A]
@pytest.fixture
def rc(monkeypatch):
    mod = _load("relax_common_followup", OPT / "relax_common.py")
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


def _skill_dir(tmp_path, skill="ke-dft-cpu"):
    d = tmp_path / "Bi2Te3" / skill
    d.mkdir(parents=True)
    (d / "POSCAR").write_text("x\n")
    return d


def _resolve(rc, cwd, monkeypatch, params=None):
    monkeypatch.chdir(cwd)
    monkeypatch.setattr(rc, "STEP_PARAMS", dict({"FUNC": "auto"}, **(params or {})))
    rc.FUNC_REASON = ""
    return rc.resolve_func(Path("incar_3d.tpl"), "step1_opt")


def _hint(d, name, text):
    h = d / ".autozt_hints"
    h.mkdir(exist_ok=True)
    (h / name).write_text(text)


def test_rerun_keeps_functional_from_local_result(rc, tmp_path, monkeypatch, capsys):
    """rerun 删了远端 step1、也没有 kl 链：老材料（pbesol）不能被键网络判据悄悄换成 pbe-d3。"""
    d = _skill_dir(tmp_path)
    _hint(d, "own_step1.txt", "FUNC=pbesol\nFUNC_REASON=step.conf 显式指定\n")
    func, _ = _resolve(rc, d, monkeypatch)
    assert func == "pbesol"
    assert "本地 result 回传" in rc.FUNC_REASON and rc.FUNC_REASON.endswith("当初：step.conf 显式指定")
    assert "远端没有本技能的 step1" in capsys.readouterr().out


def test_remote_step1_wins_over_hint_and_reason_not_nested(rc, tmp_path, monkeypatch):
    d = _skill_dir(tmp_path)
    (d / "step1_opt").mkdir()
    (d / "step1_opt" / "workflow_method.txt").write_text(
        "FUNC=pbe-d3\nFUNC_REASON=auto：沿用本技能 step1 上次的 FUNC（x）；当初：auto：键网络\n")
    _hint(d, "own_step1.txt", "FUNC=pbesol\n")
    assert _resolve(rc, d, monkeypatch)[0] == "pbe-d3"
    assert rc.FUNC_REASON.count("当初：") == 1 and rc.FUNC_REASON.endswith("当初：auto：键网络")


def test_sibling_hint_aligns_across_clusters(rc, tmp_path, monkeypatch):
    """ke、kl 在不同集群：远端没有 kl 目录，靠推来的 kl S1 记录对齐。"""
    d = _skill_dir(tmp_path)
    _hint(d, "sibling_kl-dft-cpu.txt", "FUNC=pbesol\n")
    assert _resolve(rc, d, monkeypatch)[0] == "pbesol"
    assert "配对技能" in rc.FUNC_REASON


def test_blank_func_align_means_true(rc, tmp_path, monkeypatch):
    d = _skill_dir(tmp_path)
    _hint(d, "sibling_kl-dft-cpu.txt", "FUNC=pbesol\n")
    assert _resolve(rc, d, monkeypatch, {"FUNC_ALIGN": None})[0] == "pbesol"
    assert _resolve(rc, d, monkeypatch, {"FUNC_ALIGN": False})[0] == "pbe-d3"


def test_backup_step1_dirs_ignored(rc, tmp_path, monkeypatch):
    d = _skill_dir(tmp_path)
    for name in ("step1_opt.bak", "step1_opt.stale-upstream-step1-20261014"):
        (d / name).mkdir()
        (d / name / "workflow_method.txt").write_text("FUNC=pbe\n")
    assert _resolve(rc, d, monkeypatch)[0] == "pbe-d3"          # 备份目录不算"本技能已有"


def test_reason_for_script_default_and_blank_func(rc, tmp_path, monkeypatch):
    d = _skill_dir(tmp_path)
    monkeypatch.chdir(d)
    monkeypatch.setattr(rc, "STEP_PARAMS", {})
    monkeypatch.setattr(rc, "load_step_params", lambda *a, **k: False)
    rc.resolve_func(Path("incar_3d.tpl"), "step1_opt")
    assert rc.FUNC_REASON.startswith("脚本默认值")
    monkeypatch.setattr(rc, "STEP_PARAMS", {"FUNC": ""})
    rc.resolve_func(Path("incar_3d.tpl"), "step1_opt")
    assert "为空" in rc.FUNC_REASON
    monkeypatch.setattr(rc, "STEP_PARAMS", {"FUNC": "pbe"})
    rc.resolve_func(Path("incar_3d.tpl"), "step1_opt")
    assert rc.FUNC_REASON == "step.conf 显式指定"


def test_root_and_pool_relax_common_in_sync():
    assert (ROOT / "relax_common.py").read_bytes() == (OPT / "relax_common.py").read_bytes()


@pytest.mark.parametrize("skill,own,sib", [
    ("ke-dft-cpu", "ke-dft-cpu/result/step1_opt/workflow_method.txt",
     "kl-dft-cpu/result/step1_std_opt/workflow_method.txt"),
    ("kl-dft-cpu", "kl-dft-cpu/result/step1_std_opt/workflow_method.txt",
     "ke-dft-cpu/result/step1_opt/workflow_method.txt"),
])
def test_s1_pushes_func_hints(skill, own, sib):
    sk = yaml.safe_load((ROOT / "skill" / skill / "skill.yaml").read_text(encoding="utf-8"))
    s1 = next(s for s in sk["steps"] if s["label"] == "S1_opt")
    other = "kl-dft-cpu" if skill == "ke-dft-cpu" else "ke-dft-cpu"
    assert {"from": own, "to": ".autozt_hints/own_step1.txt"} in s1["push_paths"]
    assert {"from": sib, "to": ".autozt_hints/sibling_%s.txt" % other} in s1["push_paths"]
    assert sib.split("/")[2] in ("step1_opt", "step1_std_opt")


# --------------------------------------------------------------------------- [B]
def _run(cmd, cwd):
    subprocess.run(["bash", "-c", cmd], cwd=str(cwd), check=True)


def test_push_file_entry_cmd(tmp_path):
    sys.path.insert(0, str(ROOT))
    from autozt.workflow import push_file_entry_cmd
    lp, step = tmp_path / "mat", tmp_path / "remote" / "ke-dft-cpu"
    (lp / "kl-dft-cpu" / "result" / "step1_std_opt").mkdir(parents=True)
    step.mkdir(parents=True)
    src = lp / "kl-dft-cpu" / "result" / "step1_std_opt" / "workflow_method.txt"
    src.write_text("FUNC=pbesol\n")
    entry = {"from": "kl-dft-cpu/result/step1_std_opt/workflow_method.txt",
             "to": ".autozt_hints/sibling_kl-dft-cpu.txt"}
    cmd, key, prov = push_file_entry_cmd(str(lp), str(step), entry)
    _run(cmd, step)
    dst = step / ".autozt_hints" / "sibling_kl-dft-cpu.txt"
    assert dst.read_text() == "FUNC=pbesol\n" and key == entry["to"] and prov["origin"] == "push_paths"
    src.unlink()                                                    # 本地没有了 -> 远端那份删掉
    cmd, key, _ = push_file_entry_cmd(str(lp), str(step), entry)
    _run(cmd, step)
    assert not dst.exists() and key is None
    for bad in ({"from": "x", "to": "/abs"}, {"from": "x", "to": "../up"}, {"to": "y"}, {}):
        assert push_file_entry_cmd(str(lp), str(step), bad) == ("", None, None)


# --------------------------------------------------------------------------- [C]
@pytest.fixture
def zt():
    return _load("gen_step20_zt_followup",
                 ROOT / "skill" / "zt-dft-cpu" / "step20_zt" / "gen_step20_zt.py",
                 [ROOT / "skill" / "zt-dft-cpu" / "step20_zt"])


def _zt(tmp_path, ke_func, kl_s1=None, kappa_from="kl-dft-cpu/result/step6_kappa/kappa_summary.json",
        kl_method=None):
    for k in ("structure", "kappa_L", "kl_step1"):
        (tmp_path / "inputs" / k).mkdir(parents=True, exist_ok=True)
    (tmp_path / "inputs" / "structure" / "workflow_method.txt").write_text("FUNC=%s\n" % ke_func)
    (tmp_path / "inputs" / "kappa_L" / "source.json").write_text(json.dumps({"from": kappa_from}))
    if kl_s1:
        (tmp_path / "inputs" / "kl_step1" / "workflow_method.txt").write_text(
            "FUNC=%s\nFUNC_REASON=auto：键网络\n" % kl_s1)
    if kl_method is not None:
        (tmp_path / "inputs" / "kappa_L" / "kl_method.json").write_text(kl_method)
    return tmp_path


def test_zt_falls_back_to_kl_s1(zt, tmp_path):
    r = zt.functional_check(_zt(tmp_path / "a", "pbesol", kl_s1="pbe-d3"))
    assert r["ok"] is False and r["kl_func"] == "pbe-d3" and "kl S1" in r["kl_reason"]
    r = zt.functional_check(_zt(tmp_path / "b", "pbesol", kl_s1="pbe-d3",
                                kappa_from="kl-mlff-cpu/result/step4_kappa/kappa_summary.json"))
    assert r["ok"] is None                                         # MLFF 来源不套用 kl S1 的泛函
    r = zt.functional_check(_zt(tmp_path / "c", "pbesol", kl_s1="pbesol", kl_method="[1, 2]"))
    assert r["ok"] is True                                         # 不是对象：不崩，改读 kl S1


def test_zt_pushes_kl_step1():
    sk = yaml.safe_load((ROOT / "skill" / "zt-dft-cpu" / "skill.yaml").read_text(encoding="utf-8"))
    s20 = sk["steps"][0]
    assert {"from": "kl-dft-cpu/result/step1_std_opt/workflow_method.txt",
            "to": "inputs/kl_step1/workflow_method.txt"} in s20["push_paths"]
    assert "kl_step1" not in s20["needs_results"]                  # 可选：缺了不卡 S20


# --------------------------------------------------------------------------- [D]
def test_s51_fingerprint_and_seekpath_optional():
    src = (ROOT / "skill" / "kl-dft-cpu" / "gen_step5.1_plot_phonon.py").read_text(encoding="utf-8")
    i = src.index("fp = json.dumps(")
    assert '"npoints": int(NPOINTS)' in src[i:i + 200] and '"segments": ends' in src[i:i + 200]
    assert "except ImportError as e:" in src and "seekpath 已安装" in src
    sk = yaml.safe_load((ROOT / "skill" / "kl-dft-cpu" / "skill.yaml").read_text(encoding="utf-8"))
    assert "seekpath" not in sk["requires"]["python"] and "seekpath" in sk["requires"]["optional_python"]


# --------------------------------------------------------------------------- [E]
def test_hybrid_label_by_gga(tmp_path):
    kc = _load("ke_common_followup", KE / "ke_common.py", [KE])
    for gga, label in (("PS", "HSEsol"), ("PE", "HSE06"), (None, "HSE06/HSEsol")):
        p = tmp_path / ("INCAR_%s" % gga)
        p.write_text("LHFCALC = .TRUE.\nAEXX = 0.25\nHFSCREEN = 0.2\n" + ("GGA = %s\n" % gga if gga else ""))
        assert kc.hybrid_params(p)["label"] == label


def test_scissor_label_reads_hybrid_label():
    for f in (KE / "step8_amset" / "gen_step10_amset.py", KE / "step8.4_amset2d" / "gen_step14_amset2d.py"):
        src = f.read_text(encoding="utf-8")
        assert '(d.get("hybrid") or {}).get("label")' in src
        assert '_fn = "HSE06" if "hse" in bs.parent.name' not in src
