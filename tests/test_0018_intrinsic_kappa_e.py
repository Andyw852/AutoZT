"""[0018] postprocess_intrinsic.py：本征（剔除 IMP，或 --drop IMP,POP 只留 ADP）结果写出 κe。纯本地。"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
S84 = ROOT / "skill" / "ke-dft-cpu" / "step8.4_amset2d"


def _pp():
    for d in (S84, ROOT / "skill" / "ke-dft-cpu", ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location("ppi_0018", str(S84 / "postprocess_intrinsic.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _t(v):
    return np.array([[[[v, 0, 0], [0, v, 0], [0, 0, 1.0]]]])     # [doping=1][T=1][3][3]


def test_intrinsic_json_has_kappa_e(tmp_path):
    pp = _pp()
    tr = {"conductivity_S_m": _t(1e4), "seebeck_uV_K": _t(-200.0),
          "electronic_thermal_conductivity": _t(0.05),
          "mobility_cm2_Vs_s": {"overall": _t(80.0), "ADP": _t(80.0)}}
    res = pp.build_result(tmp_path, tmp_path / "mesh_1.h5", "0.5.1", {}, ["ADP", "IMP", "POP"],
                          ["IMP", "POP"], ["ADP"], tr, [-1e17], [300.0])
    k = res["intrinsic"]
    assert k["electronic_thermal_conductivity_inplane_W_mK"] == [[0.05]]
    assert np.asarray(k["electronic_thermal_conductivity_W_mK"]).shape == (1, 1, 3, 3)
    assert res["intrinsic_labels"] == ["ADP"] and pp.__version__.startswith("2026-10-08")


def test_s84_chain_runs_adp_only_postprocess_when_asked():
    src = (S84 / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    assert '"INTRINSIC_ADP_ONLY": (INTRINSIC_ADP_ONLY, "bool")' in src
    i = src.index("if INTRINSIC_ADP_ONLY:")
    assert "--drop IMP,POP,PIE" in src[i:i + 300] and "intrinsic_ADP.json" in src[i:i + 300]
    j = src.index('_acmd = _acmd + " && python postprocess_intrinsic.py --check-reproduce"')
    assert j < i < src.index("AMSET_TAIL", i)          # 在复现硬门之后、落成 transport.json 之前


def test_adp_only_json_is_archived_with_s84():
    sys.path.insert(0, str(ROOT / "skill" / "_common" / "opt"))
    spec = importlib.util.spec_from_file_location("kc_0018", str(ROOT / "skill" / "ke-dft-cpu" / "ke_common.py"))
    kc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kc)
    assert "intrinsic_ADP.json" in kc.DONE_MARKERS["step8.4_amset2d"]
