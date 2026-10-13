"""[V186] S7.1 amset deform read 的 NBANDS 对齐（_NBANDS_FIX_CODE）：
1. bulk 带多于某形变构型（VASP 默认 NBANDS 随并行布局取整）-> bulk 与已读构型统一裁到最小带数，不再跳过构型；
2. 形变构型带更多 -> 照 amset 原样裁形变构型；
3. 要裁的带落进能量窗口 -> 报错退出（exit 3），不静默裁；
4. gen 的 amset deform read 命令经这段代码调用。"""
import ast
import sys
from pathlib import Path

import numpy as np
import pytest

GEN = Path(__file__).resolve().parent / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py"


def _code():
    tree = ast.parse(GEN.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "_NBANDS_FIX_CODE":
            return node.value.value
    raise AssertionError("gen 里没有 _NBANDS_FIX_CODE")


@pytest.fixture()
def fix():
    pytest.importorskip("amset")
    try:
        import spglib as _sg                     # 新 spglib 去掉了 spglib.spglib，老 amset 还在 import 它
        if not hasattr(_sg, "spglib"):
            _sg.spglib = _sg
            sys.modules["spglib.spglib"] = _sg
    except ImportError:
        pass
    import amset.tools.deformation as T
    orig = T.check_calculation
    ns = {"__name__": "nbands_fix_test"}
    exec(_code(), ns)                                # noqa: S102  只装补丁，不调 CLI
    yield ns
    T.check_calculation = orig


def _calc(levels):
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.bandstructure import BandStructure
    from pymatgen.electronic_structure.core import Spin
    st = Structure(Lattice.hexagonal(3.2, 20.0), ["Mo"], [[0, 0, 0]])
    kf = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0.5, 0.5, 0]]
    bands = np.array([[e, e + 0.1, e + 0.2, e + 0.3] for e in levels])
    bs = BandStructure(kf, {Spin.up: bands}, st.lattice.reciprocal_lattice, 0.0, structure=st)
    return {"reference": 0.0, "bandstructure": bs}


def test_bulk_more_bands_trims_instead_of_skipping(fix):
    bulk = _calc([-1.0, 1.0, 4.0, 10.0])
    d1 = _calc([-1.0, 1.0, 4.0, 10.0])
    d2 = _calc([-1.0, 1.0, 4.0])                     # 少 1 条空带
    assert fix["check_calculation"](bulk, d1) is d1
    out = fix["check_calculation"](bulk, d2)
    assert out is d2                                 # 以前：False（跳过）
    for c in (bulk, d1, d2):
        assert c["bandstructure"].nb_bands == 3
        assert all(len(b) == 3 for b in c["bandstructure"].bands.values())


def test_deformation_more_bands_still_trimmed_by_amset(fix):
    bulk = _calc([-1.0, 1.0, 4.0])
    d = _calc([-1.0, 1.0, 4.0, 10.0])
    out = fix["check_calculation"](bulk, d)
    assert all(len(b) == 3 for b in out["bandstructure"].bands.values())
    assert bulk["bandstructure"].nb_bands == 3


def test_refuses_to_drop_band_inside_window(fix):
    bulk = _calc([-1.0, 1.0, 2.0])                   # 第 3 条在 CBM + 1.5 eV 窗口内
    d = _calc([-1.0, 1.0])
    with pytest.raises(SystemExit) as e:
        fix["check_calculation"](bulk, d)
    assert e.value.code == 3


def test_gen_routes_deform_read_through_fix():
    src = GEN.read_text(encoding="utf-8")
    assert "python3 -c %s undeformed %s" in src and "shlex.quote(_NBANDS_FIX_CODE)" in src
    assert "&& amset deform read undeformed" not in src


def test_cli_path_runs(fix):
    """python -c <code> --help 走到 amset 的 deform read 子命令（参数照传，退出码 0）。"""
    import os
    import subprocess
    r = subprocess.run([sys.executable, "-c", _code(), "--help"], capture_output=True, text=True,
                       env=dict(os.environ), timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "BULK_FOLDER" in r.stdout and "--energy-cutoff" in r.stdout
