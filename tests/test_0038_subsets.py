"""[0038] SUBSETS 散射机制子集：subset_<机制>.json 生成（gen 侧 subset_specs）与 S8.3 分类/读取。纯本地。"""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _load(name, path):
    for d in (ROOT, KE, ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _g14():
    return _load("g14_0038", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")


def _g10():
    return _load("g10_0038", KE / "step8_amset" / "gen_step10_amset.py")


def _s83():
    return _load("s83_0038", KE / "step8.3_output" / "gen_step13_output.py")


def test_subset_specs_degenerate_naming():
    g = _g14()
    # 出厂三套 + 默认 2D SCATTERING（ADP/IMP/POP，无 ODP/PIE）→ 退化为用户示例 ADP+POP+IMP / ADP+POP / ADP
    specs = g.subset_specs(None, ["ADP", "IMP", "POP"])
    assert specs == [
        (["ADP", "POP", "IMP"], [], "subset_ADP+POP+IMP.json"),
        (["ADP", "POP"], ["IMP"], "subset_ADP+POP.json"),
        (["ADP"], ["IMP", "POP"], "subset_ADP.json"),
    ]
    # 全机制（含 ODP/PIE）→ 用户要的四套里的前三套（第四套 BT2×DPT-τ 来自 S8.1）
    specs = g.subset_specs(None, ["ADP", "IMP", "POP", "ODP", "PIE"])
    assert [n for _, _, n in specs] == [
        "subset_ADP+POP+IMP+ODP+PIE.json",
        "subset_ADP+POP+ODP+PIE.json",
        "subset_ADP.json",
    ]
    assert [d for _, d, _ in specs] == [[], ["IMP"], ["IMP", "POP", "ODP", "PIE"]]


def test_subset_specs_off_bad_and_dedupe():
    g = _g14()
    assert g.subset_specs("off", ["ADP", "IMP", "POP"]) == []
    assert g.subset_specs("", ["ADP"]) == []
    assert g.subset_specs("none", ["ADP"]) == []
    # 重复子集去重；不认识的机制报错
    assert len(g.subset_specs("ADP; ADP", ["ADP"])) == 1
    with pytest.raises(SystemExit):
        g.subset_specs("ADP+ODQ", ["ADP"])


def test_subset_specs_3d_matches_2d():
    g14, g10 = _g14(), _g10()
    sc = ["ADP", "IMP", "POP", "ODP", "PIE"]
    assert g14.subset_specs(None, sc) == g10.subset_specs(None, sc)


def test_s83_subset_classification():
    S = _s83()
    assert S._subset_key({"intrinsic_labels": ["ADP", "POP", "IMP"], "intrinsic_dropped": []}) == "overall"
    assert S._subset_key({"intrinsic_labels": ["ADP", "POP", "ODP", "PIE"], "intrinsic_dropped": ["IMP"]}) == "intrinsic"
    assert S._subset_key({"intrinsic_labels": ["ADP"], "intrinsic_dropped": ["IMP", "POP"]}) == "ADP_only"
    assert S._subset_key({"intrinsic_labels": ["ADP", "POP"], "intrinsic_dropped": ["IMP", "PIE"]}) == "ADP_POP"


def test_s83_subset_file_defs_skips_overall(tmp_path):
    S = _s83()
    run = tmp_path / "step8.4_amset2d"
    run.mkdir()
    for name, kept, dropped in (
        ("subset_ADP+POP+IMP.json", ["ADP", "POP", "IMP"], []),       # 全机制 = overall，跳过
        ("subset_ADP+POP.json", ["ADP", "POP"], ["IMP"]),             # intrinsic
        ("subset_ADP.json", ["ADP"], ["IMP", "POP"]),                 # ADP_only
    ):
        (run / name).write_text(json.dumps({"intrinsic_labels": kept, "intrinsic_dropped": dropped}))
    defs = S._subset_file_defs(run)
    assert {d[0] for d in defs} == {"intrinsic", "ADP_only"}
    assert {d[1] for d in defs} == {"subset_ADP+POP.json", "subset_ADP.json"}
