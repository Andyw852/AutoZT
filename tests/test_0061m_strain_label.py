"""[0061-m 2026-10-22] AMSET get_symmetrized_strain_mapping 的应变标签改正（dp_symmetrize.symmetrized_strain_mapping）。

AMSET 把形变构型 E 的能带按 frac_op 搬到 k' = M·k（M = spglib Wᵀ），却把它标成 R·E·Rᵀ；这组能带属于 Rᵀ·E·R。
EPM 硅（金刚石 Fd-3m，只算 xx±、yz± 四个构型 = amset deform create 的约化集），桩 AMSET（只换 amset 的几个
k 点 / 对称工具，补应变的循环就是 dp_symmetrize 里那份）：改正版补出来的 8 个应变，能带与该应变下直接解的本征值一致；
把标签换回 AMSET 原版（R·E·Rᵀ），xy、xz 有对不上的（Si 实跑 S7.1 的 xy/xz 对齐离散度 2.7 / 3.9 eV 就是它）。
"""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "skill" / "ke-dft-cpu"), str(Path(__file__).resolve().parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
pytest.importorskip("pymatgen")
import _epm_si as E  # noqa: E402
import dp_symmetrize as DPS  # noqa: E402

EC, NB, D = 80.0, 6, 0.005
_G = np.arange(3) / 3.0
K = np.array([[a, b, c] for a in _G for b in _G for c in _G])
K = K - np.round(K)

_FAKE = {
    "__init__.py": "",
    "deformation/__init__.py": "",
    "deformation/potentials.py": "_mapping_tol = 0.002\n",
    "electronic_structure/__init__.py": "",
    "electronic_structure/kpoints.py": '''import numpy as np


def kpoints_to_first_bz(k, tol=1e-5, negative_zone_boundary=True):
    k = np.asarray(k, float)
    return k - np.round(k)


def get_kpoints_from_bandstructure(bs, cartesian=False, sort=False):
    k = np.asarray(bs.k, float)
    return k[np.lexsort(np.round(k, 6).T[::-1])] if sort else k
''',
    "electronic_structure/symmetry.py": '''from types import SimpleNamespace

import numpy as np
from pymatgen.core.operations import SymmOp

from amset.electronic_structure.kpoints import kpoints_to_first_bz


def get_symmops(structure, symprec=0.01):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    return [SymmOp.from_rotation_and_translation(rotation_matrix=o.rotation_matrix.T,
                                                 translation_vec=o.translation_vector)
            for o in SpacegroupAnalyzer(structure, symprec=symprec).get_symmetry_operations(cartesian=False)]


def expand_bandstructure(bs, symprec=0.01, time_reversal=True):
    return bs


def rotate_bandstructure(bs, frac_symop):
    k = np.dot(frac_symop.rotation_matrix, np.asarray(bs.k, float).T).T
    return SimpleNamespace(k=kpoints_to_first_bz(k, negative_zone_boundary=False), bands=bs.bands)
''',
}


@pytest.fixture()
def fake_amset(tmp_path, monkeypatch):
    for rel, txt in _FAKE.items():
        f = tmp_path / "amset" / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(txt)
    monkeypatch.syspath_prepend(str(tmp_path))
    for n in [x for x in sys.modules if x == "amset" or x.startswith("amset.")]:
        monkeypatch.delitem(sys.modules, n)
    yield
    for n in [x for x in sys.modules if x == "amset" or x.startswith("amset.")]:   # 桩模块别漏给后面的测试
        del sys.modules[n]


def _eig(strain, k):
    return E.solve(k, E.lattice(np.asarray(strain, float)), EC)[1][:NB]


@pytest.fixture(scope="module")
def si():
    from pymatgen.analysis.elasticity.strain import Strain
    from pymatgen.core import Lattice, Structure
    bulk = Structure(Lattice(E.L0), ["Si", "Si"], [E.TAU, 1 - E.TAU])
    comp = []                                   # [(Strain, (nb, nk) 本征值)]；Strain 在新 numpy 上不能当 dict 键
    for P in [(0, 0), (1, 2)]:
        for s in (1, -1):
            st = Strain(E.pattern(P, s * D))
            comp.append((st, np.array([_eig(st, k) for k in K]).T))
    return bulk, comp


def _mapping(comp):
    from pymatgen.core.tensors import TensorMapping
    sm = TensorMapping()
    for st, ev in comp:
        sm[st] = {"reference": 0.0, "bandstructure": SimpleNamespace(k=K.copy(), bands={"up": ev})}
    return sm


def _errors(sm, comp):
    """补出来的每个应变：它的能带与该应变下直接解的本征值之差（eV）。"""
    out = {}
    for st, calc in sm.items():
        if any(np.allclose(np.asarray(st), np.asarray(c)) for c, _ in comp):
            continue
        T = np.asarray(st, float)
        i, j = np.unravel_index(np.abs(T).argmax(), (3, 3))
        name = "%s%s%+d" % ("xyz"[min(i, j)], "xyz"[max(i, j)], int(np.sign(T[i, j])))
        bs = calc["bandstructure"]
        out[name] = max(float(np.abs(_eig(T, k) - bs.bands["up"][:, n]).max()) for n, k in enumerate(bs.k))
    return out


def test_fixed_mapping_labels_match_bands(si, fake_amset):
    bulk, comp = si
    rep = {}
    err = _errors(DPS.symmetrized_strain_mapping(bulk, _mapping(comp), report=rep), comp)
    assert sorted(err) == ["xy+1", "xy-1", "xz+1", "xz-1", "yy+1", "yy-1", "zz+1", "zz-1"]
    assert max(err.values()) < 1e-5, err
    assert rep["added"] == 8 and rep["relabeled"] > 0          # 有补出来的应变原版会标成别的


def test_amset_original_label_is_wrong(si, fake_amset, monkeypatch):
    """同一循环换回 AMSET 原版标签 strain.rotate(R)：剪切补出来的有对不上的（C3[111] 把 xy、xz 互换）。"""
    bulk, comp = si
    fixed = DPS.strain_label_rotation
    monkeypatch.setattr(DPS, "strain_label_rotation", lambda L, M: fixed(L, M).T)
    err = _errors(DPS.symmetrized_strain_mapping(bulk, _mapping(comp)), comp)
    bad = sorted(n for n, e in err.items() if e > 0.05)
    assert bad and all(n[:2] in ("xy", "xz") for n in bad), err
    assert max(err[n] for n in err if n[:2] in ("yy", "zz")) < 1e-5


def test_label_rotation_is_inverse_for_noninvolutory():
    """对合操作（R² = 1）两种标签相同；C3[111] 不同：原版 R·E_yz·Rᵀ 与改正 Rᵀ·E_yz·R 落在不同剪切分量上。"""
    L = E.L0
    W_c3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)          # fcc 原胞基矢的轮换 = 笛卡尔 C3[111]
    rt = DPS.strain_label_rotation(L, W_c3.T)
    R = rt.T
    assert np.allclose(R @ R.T, np.eye(3)) and np.allclose(np.linalg.matrix_power(R, 3), np.eye(3))
    Eyz = E.pattern((1, 2), D)
    new, old = rt @ Eyz @ rt.T, R @ Eyz @ R.T
    nz = lambda M: tuple(sorted(np.argwhere(np.abs(np.triu(M)) > 1e-12)[0]))   # noqa: E731
    assert nz(new) != nz(old) and {nz(new), nz(old)} == {(0, 1), (0, 2)}
    W_m = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], float)           # 交换 a1、a2：对合
    rt = DPS.strain_label_rotation(L, W_m.T)
    assert np.allclose(rt @ Eyz @ rt.T, rt.T @ Eyz @ rt)


GEN = ROOT / "skill" / "ke-dft-cpu" / "step7_deform" / "step7b_read" / "gen_step9b_deform_read.py"


def _embedded(var):
    import ast
    tree = ast.parse(GEN.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == var:
            return node.value.value
    raise AssertionError("gen 里没有 %s" % var)


def test_amset_read_wrapper_installs_fixed_mapping(si, fake_amset, tmp_path, monkeypatch, capsys):
    """S7.1 的 amset deform read 包装（_NBANDS_FIX_CODE）在桩 AMSET 里执行：read() 调用时取到的
    get_symmetrized_strain_mapping 是改正版，补出来的应变能带都对，日志写明原版会标错几个。"""
    extra = {
        "constants.py": "defaults = {'energy_cutoff': 1.5}\n",
        "tools/__init__.py": "",
        "tools/deformation.py": "def check_calculation(bulk, deformation):\n    return deformation\n",
        "electronic_structure/common.py": "def get_energy_cutoffs(ecut, bs):\n    return -ecut, ecut\n",
        "deformation/potentials.py": (_FAKE["deformation/potentials.py"]
                                      + "def calculate_deformation_potentials(b, sm):\n    return {}\n"
                                      + "def extract_bands(dp, ib):\n    return dp\n"
                                      + "def get_symmetrized_strain_mapping(*a, **k):\n    raise AssertionError('原版')\n"),
        "deformation/io.py": "def write_deformation_potentials(dp, k, s, filename='deformation.h5'):\n    return filename\n",
    }
    root = [p for p in sys.path if (Path(p) / "amset" / "deformation").is_dir()][0]
    for rel, txt in extra.items():
        f = Path(root) / "amset" / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(txt)
    monkeypatch.setattr(sys, "argv", ["deform_read.py"])
    exec(_embedded("_NBANDS_FIX_CODE"), {"__name__": "s71_wrapper_test"})       # noqa: S102  只装补丁，不调 CLI
    from amset.deformation.potentials import get_symmetrized_strain_mapping      # read() 在函数体里这样取
    bulk, comp = si
    sm = get_symmetrized_strain_mapping(bulk, _mapping(comp), symprec=0.01, symprec_deformation=1e-4)
    err = _errors(sm, comp)
    assert len(err) == 8 and max(err.values()) < 1e-5, err
    assert "[strain_label_fix] 对称补了 8 个应变" in capsys.readouterr().out


def test_vac_path_uses_fixed_mapping():
    """2D 真空参考那条路径（_vac_code）直接调 AMSET 的函数，也要换成改正版。"""
    src = _embedded("_vac_code")
    assert "get_symmetrized_strain_mapping(" not in src
    assert "dp_symmetrize.symmetrized_strain_mapping(bulk_struct, sm" in src
    assert src.index("import dp_symmetrize") < src.index("dp_symmetrize.symmetrized_strain_mapping(")
