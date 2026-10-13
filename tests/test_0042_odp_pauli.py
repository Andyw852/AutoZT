"""[0042] ODP/IVS 的泡利阻塞（简并统计）：
1. 费米能级远在带下（非简并）-> 与不传 fermi_levels 的结果逐点相同（Kaasbjerg Eq.31 极限）；
2. 简并电子：发射阈下 r/r0 = (N + f(x+ħω))/N；阈上 = (2N + 1 + f(x+ħω) − f(x−ħω))/(2N + 1)；
3. 空穴用 f_h = 1 − f_e（x = −E、x_F = −E_F）；
4. 形状 (ndop, nT, nb, nk)；温度数与 fermi_levels 不符时报错；
5. 真 AMSET（有就跑）：3D 插件 _dpo_rates 逐掺杂给出不同的率，非简并掺杂与旧结果相同。"""
import importlib.util
import json
import math
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
S84 = KE / "step8.4_amset2d"


def _load(name, path):
    for d in (ROOT, KE, S84):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def dpo():
    return _load("dp_optical_0042", S84 / "dp_optical.py")


A_LAT, M_D, CELL_KG, HW = 3.2, 0.5, 2.6e-25, 0.040
VAL = [{"label": "A", "g": 2, "k_frac": [[0.25, 0.25, 0.0], [0.75, 0.75, 0.0]]}]
KT300 = 8.617333262e-5 * 300.0


@pytest.fixture(scope="module")
def band():
    """方格子、两个等价抛物谷（m = 0.5），谷底 0 eV；返回 eV 能量、k、倒格子。"""
    n = 300
    rlat = 2 * np.pi / A_LAT * np.eye(3)
    rlat[2, 2] = 2 * np.pi / 20.0
    g = np.arange(n) / n
    k1, k2 = np.meshgrid(g, g, indexing="ij")
    kf = np.stack([k1.ravel(), k2.ravel(), np.zeros(n * n)], 1)
    hb2m = 3.80998212 / M_D
    E = np.full(len(kf), np.inf)
    for c in VAL[0]["k_frac"]:
        d = kf - np.array(c)
        d = d - np.rint(d)
        d[:, 2] = 0
        E = np.minimum(E, hb2m * np.sum((d @ rlat) ** 2, axis=1))
    return E, kf, rlat


def _rates(dpo, band, sign, fermi_eV=None):
    E, kf, rlat = band
    ch = dpo.parse_channels("A>A:3.0@%g" % (HW * 1e3))
    fl = None if fermi_eV is None else np.atleast_2d(np.asarray(fermi_eV, float)) / dpo.HA_EV
    r = dpo.carrier_rates((sign * E / dpo.HA_EV)[None], kf, rlat, VAL, ch, [300.0], CELL_KG, {"intra"}, sign,
                          fermi_levels=fl)
    return r


def test_nondegenerate_limit_equals_old(dpo, band):
    r0 = _rates(dpo, band, +1)[0, 0]
    rp = _rates(dpo, band, +1, [[-1.0]])[0, 0, 0]                 # E_F 在带底下 1 eV
    np.testing.assert_allclose(rp, r0, rtol=1e-9)


@pytest.mark.parametrize("sign", [+1, -1])
def test_degenerate_factors(dpo, band, sign):
    """电子：E_F = +0.05 eV（进带）；空穴：价带 E = −E_band，E_F = −0.05 eV（进带），载流子能量 x = E_band。"""
    E = band[0]
    ef = 0.05 * sign
    r0 = _rates(dpo, band, sign)[0, 0]
    rp = _rates(dpo, band, sign, [[ef]])[0, 0, 0]
    N = 1.0 / math.expm1(HW / KT300)
    f = lambda x: 1.0 / (math.exp((x - 0.05) / KT300) + 1.0)       # noqa: E731  载流子占据（x_F = 0.05）
    for x, emit in ((0.015, False), (0.160, True)):
        sel = np.abs(E - x) < 0.002
        assert sel.sum() > 20
        got = float(np.mean(rp[sel] / r0[sel]))
        want = ((N + f(x + HW)) / N) if not emit else ((2 * N + 1 + f(x + HW) - f(x - HW)) / (2 * N + 1))
        assert abs(got / want - 1) < 0.03, (sign, x, got, want)
    # 简并时发射被阻塞、吸收被增强：阈下（只有吸收）率升高
    sel = np.abs(E - 0.015) < 0.002
    assert np.mean(rp[sel]) > 1.05 * np.mean(r0[sel])


def test_shapes_and_errors(dpo, band):
    r = _rates(dpo, band, +1, [[-1.0], [0.0], [0.05]])
    assert r.shape == (3, 1, 1, len(band[0]))
    assert not np.allclose(r[0], r[2])
    with pytest.raises(ValueError):
        _rates(dpo, band, +1, [[0.0, 0.0]])                        # 2 个温度的费米能级 vs 1 个温度


def test_plugin_rates_depend_on_doping(tmp_path):
    """真 AMSET：3D 插件逐掺杂算率（以前广播同一份）。"""
    pytest.importorskip("amset")
    try:
        import spglib as _sg
        if not hasattr(_sg, "spglib"):
            _sg.spglib = _sg
            sys.modules["spglib.spglib"] = _sg
    except ImportError:
        pass
    for src in (KE / "step8_amset" / "amset_dpo3d.py", S84 / "dp_optical.py"):
        shutil.copy(src, tmp_path)
    rec = {"mechanisms": ["ODP"],
           "dp_optical": {"electron": {"valleys": [{"label": "Γ", "g": 1, "k_frac": [[0, 0, 0]]}],
                                       "channels": [{"from": "Γ", "to": "Γ", "D_eV_A": 3.0, "hw_meV": 40.0,
                                                     "kind": "intra"}]},
                          "adp_intravalley": {"electron": False, "hole": False}}}
    (tmp_path / "dp_optical_3d.json").write_text(json.dumps(rec))
    old = os.getcwd()
    os.chdir(tmp_path)
    sys.path.insert(0, str(tmp_path))
    try:
        for m in ("amset_dpo3d", "dp_optical"):
            sys.modules.pop(m, None)
        import amset_dpo3d as P
        from amset.constants import hartree_to_ev
        from pymatgen.core import Lattice, Structure
        from pymatgen.electronic_structure.core import Spin
        st = Structure(Lattice.cubic(5.43), ["Si", "Si"], [[0, 0, 0], [0.25, 0.25, 0.25]])
        n = 16
        g = (np.arange(n) - n // 2) / n
        kf = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
        kc = kf @ st.lattice.reciprocal_lattice.matrix
        Ev = -np.sum(kc ** 2, 1) * 0.5
        Ec = 1.0 + np.sum(kc ** 2, 1) * 0.5
        fl = np.array([[-1.0, -1.0], [0.2, 0.2], [1.05, 1.05]]) / hartree_to_ev   # 价带内 / 离导带 0.8 eV / 导带内（简并 n）
        ad = SimpleNamespace(structure=st, kpoints=kf, temperatures=np.array([300.0, 600.0]),
                             fermi_levels=fl, spins=[Spin.up],
                             energies={Spin.up: np.stack([Ev, Ec]) / hartree_to_ev},
                             is_metal=False, vb_idx={Spin.up: 0}, doping=np.array([1e17, 0.0, -1e21]))
        r = P._dpo_rates(ad, P._DPO_KINDS["ODP"])[Spin.up]
        assert r.shape == (3, 2, 2, len(kf))
        assert np.allclose(r[0, :, 1], r[1, :, 1], rtol=1e-5)       # 导带电子在非简并两档：相同（f < 1e-6）
        assert not np.allclose(r[2, :, 1], r[1, :, 1])               # 简并 n 型：泡利阻塞改变了率
    finally:
        os.chdir(old)
        sys.path.remove(str(tmp_path))
