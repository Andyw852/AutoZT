"""[0044] postprocess_intrinsic --check-reproduce：
1. 判据 = max |a − b| / max(|b|, 1)：单机制 μ 在零前因子处 ~1e47（本征点 IMP），浮点舍入（相对 ~1e-15）不再判 FAIL；
2. 小量仍按绝对差卡死；真不同源（相对 1e-3）照样 FAIL；
3. S8.4 的 POSTPROCESS_ONLY 只要 mesh_*.h5 + transport_*.json（后处理失败那一轮没有 transport.json）。"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
S84 = ROOT / "skill" / "ke-dft-cpu" / "step8.4_amset2d"


@pytest.fixture(scope="module")
def P():
    pytest.importorskip("pymatgen")
    for d in (ROOT, S84, S84.parent):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("pp_0044", str(S84 / "postprocess_intrinsic.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _case(scale_imp=3.875e47, rel_noise=2.46e-15, sig_err=0.0):
    t = np.ones((3, 2, 3, 3))
    ref = {"mobility": {"overall": (393.4 * t).tolist(), "IMP": (scale_imp * t).tolist()},
           "conductivity": (1e5 * t).tolist(), "seebeck": (-150.0 * t).tolist()}
    mob0 = {"overall": 393.4 * t * (1 + rel_noise), "IMP": scale_imp * t * (1 + rel_noise)}
    return mob0, 1e5 * t * (1 + sig_err), -150.0 * t, ref


def test_huge_single_mech_mobility_passes(P):
    mob0, sig0, seeb0, j = _case()
    diffs, rels = P.reproduce_diffs(mob0, sig0, seeb0, j, ["IMP"])
    assert diffs["IMP"] > 1e32                       # 旧的绝对判据：必 FAIL
    assert max(rels.values()) < 1e-12               # 新判据：PASS


def test_real_mismatch_still_fails(P):
    mob0, sig0, seeb0, j = _case(rel_noise=0.0, sig_err=1e-3)
    _, rels = P.reproduce_diffs(mob0, sig0, seeb0, j, ["IMP"])
    assert rels["conductivity_S_m"] == pytest.approx(1e-3)
    small = {"overall": np.full((1, 1, 3, 3), 0.5)}
    _, rels = P.reproduce_diffs({"overall": small["overall"] + 2e-4}, np.zeros((1, 1, 3, 3)), np.zeros((1, 1, 3, 3)),
                                {"mobility": {"overall": small["overall"].tolist()},
                                 "conductivity": np.zeros((1, 1, 3, 3)).tolist(),
                                 "seebeck": np.zeros((1, 1, 3, 3)).tolist()}, [])
    assert rels["overall"] == pytest.approx(2e-4)    # |ref| ≤ 1：仍是绝对差


def test_main_uses_relative_criterion():
    src = (S84 / "postprocess_intrinsic.py").read_text(encoding="utf-8")
    assert "diffs, rels = reproduce_diffs(mob0, sig0, seeb0, j, labels)" in src
    assert 'reproduce["passed"] = bool(worst <= args.reproduce_tol)' in src
    assert 'reproduce["max_rel_diff"] = rels' in src


def test_postprocess_only_does_not_need_transport_json():
    src = (S84 / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert "'ls mesh_*.h5 transport_*.json >/dev/null 2>&1 ||" in src
    assert "transport_*.json transport.json >/dev/null" not in src
