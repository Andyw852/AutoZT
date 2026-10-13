"""[0042] 2D 极性判据：面内极性只看面内模有效电荷；面外极性只对有 σh 的材料（gen 对没有 σh 的传 polar_z_tol=inf）。
用 MoS₂ 型三原子胞的 Γ 点声子（test_0035 的弹簧模型），Born 电荷改成 Janus 式（两个 S 层不等价：1.0 / 1.3 e）。"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


def _load(name, path):
    for d in (ROOT, KE, ROOT / "tests"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


OF = _load("odp_frozen_0042", KE / "odp_frozen.py")
T = _load("t0035_for_0042", ROOT / "tests" / "test_0035_odp_frozen_tau.py")


def _table(born, **kw):
    f, e = T._phonons()
    masses = list(T.MASS) if hasattr(T, "MASS") else [95.95, 32.06, 32.06]
    modes, groups, notes = OF.mode_table(f, e, masses, born, **kw)
    return modes, groups, e


JANUS = np.array([-2.3 * np.eye(3), 1.0 * np.eye(3), 1.3 * np.eye(3)])


def _a1_group(modes, groups, e):
    f, ev = T._phonons()
    idx = None
    for md in modes:
        v = ev[np.argsort(f, kind="stable")[md["index"] - 1]]
        if np.linalg.norm(v[0]) < 1e-6 and abs(v[1, 2]) > 0.5 and abs(v[1, 2] + v[2, 2]) < 1e-6:
            idx = md["index"]
    assert idx is not None
    md = next(m for m in modes if m["index"] == idx)
    return md, next(g for g in groups if idx in g["modes"])


def test_janus_a1_kept_without_sigma_h():
    modes, groups, e = _table(JANUS, polar_z_tol=float("inf"), two_d=True)     # 没有 σh
    md, g = _a1_group(modes, groups, e)
    assert md["Z_mode_inplane_e"] < 1e-6 and md["Z_mode_z_e"] > 0.1           # 只有面外偶极，且 > 0.1 e
    assert g["use"] and not g["polar"]                                        # 以前：总模有效电荷 > 0.1 -> 被当面内极性剔掉


def test_janus_a1_excluded_with_sigma_h_threshold():
    modes, groups, e = _table(JANUS, polar_z_tol=0.02, two_d=True)            # 有 σh（或判不出）：照旧按面外剔
    md, g = _a1_group(modes, groups, e)
    assert g["polar"] and not g["use"] and md["polar_kind"] == "面外"


def test_mos2_unchanged_and_3d_uses_total():
    born = np.array([-2 * np.eye(3), np.eye(3), np.eye(3)])
    m2, g2, e = _table(born, polar_z_tol=0.02, two_d=True)
    sel = [i for g in g2 if g["use"] for i in g["modes"]]
    assert sel == [4, 5, 6]                                                   # 与 test_0035 的 MoS₂ 选择相同
    m3, g3, _ = _table(JANUS, polar_z_tol=float("inf"), two_d=False)         # 3D：总模有效电荷
    md, g = _a1_group(m3, g3, e)
    assert g["polar"] and md["polar_kind"] == "体相"
