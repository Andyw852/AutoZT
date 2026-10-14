"""[0050] 冻声子法谷间零阶形变势（IVS，K↔K′）核心模块 ivs_frozen.py：
1. fold_matrix：六角 K 星 -> √3×√3（det 3、最紧凑）；Q 这类不可公度的谷 -> 报错；
2. build_supercell / translation_perms / gamma_weight：原胞 Γ 来的模权重 1、区边界模权重 0；
3. 三角晶格紧束缚玩具模型（单轨道、导带底在 K/K′）：给定在位耦合 ε_j = α·Δz_j，整条流程
   （选模 -> displacement 归一化 -> ±u 能级 -> 标来源 -> 劈裂迹 -> 伙伴模求和）算出的 D 必须等于 α；
   只耦合 z 的模型里面内模的 D = 0；
4. label_folded / edge_group：常数平移下正确标出 K/K′ 二重态；输入对不上、来源附近有别的能级 -> 报错；
5. inter_channels：D² 按 ħω 合并，小于 dmin 丢掉，kind = inter_star。"""
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"


@pytest.fixture(scope="module")
def F():
    if str(KE) not in sys.path:
        sys.path.insert(0, str(KE))
    spec = importlib.util.spec_from_file_location("ivs_frozen_0050", str(KE / "ivs_frozen.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


A = 3.2
LAT = np.array([[A, 0, 0], [-A / 2, A * math.sqrt(3) / 2, 0], [0, 0, 20.0]])     # γ = 120°
K_STAR = [[1 / 3, 1 / 3, 0], [2 / 3, 2 / 3, 0]]                                    # K、K′
T_HOP = 0.4                     # eV；E(k) = 2t Σcos，t > 0 时带底在 K/K′（−3t），Γ 在 +6t
NEIGH = [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]                       # γ=120°：±a1、±a2、±(a1+a2)


def _band(k):
    return 2 * T_HOP * (math.cos(2 * math.pi * k[0]) + math.cos(2 * math.pi * k[1])
                        + math.cos(2 * math.pi * (k[0] + k[1])))


def test_fold_matrix_hexagonal_k(F):
    M = F.fold_matrix(K_STAR, LAT[:2, :2])
    assert int(round(np.linalg.det(M))) == 3
    v = M @ LAT[:2, :2]
    assert np.allclose(np.linalg.norm(v, axis=1), A * math.sqrt(3))        # √3×√3
    for k in K_STAR:
        assert np.allclose((M @ np.array(k[:2])) % 1.0, 0.0, atol=1e-9) or \
            np.allclose(np.abs((M @ np.array(k[:2])) - np.rint(M @ np.array(k[:2]))), 0.0, atol=1e-9)
    fk = F.folded_kpoints(M)
    assert len(fk) == 3 and np.allclose(fk[0], 0.0)


def test_fold_matrix_rejects_incommensurate(F):
    with pytest.raises(ValueError):
        F.fold_matrix([[0.2317, 0.2317, 0], [-0.2317, -0.2317, 0]], LAT[:2, :2])


def _supercell(F, nat_prim=1):
    M = F.fold_matrix(K_STAR, LAT[:2, :2])
    frac = [[0, 0, 0.5]] if nat_prim == 1 else [[0, 0, 0.5], [1 / 3, 2 / 3, 0.42], [1 / 3, 2 / 3, 0.58]]
    spec = ["X"] if nat_prim == 1 else ["Mo", "S", "S"]
    sc = F.build_supercell(LAT, frac, spec, M)
    return M, sc


def test_supercell_and_gamma_weight(F):
    M, sc = _supercell(F, nat_prim=3)
    assert len(sc["frac"]) == 9 and sc["species"] == ["Mo"] * 3 + ["S"] * 6
    perms = F.translation_perms(sc)
    assert len(perms) == 3 and list(perms[0]) == list(range(9))
    rng = np.random.default_rng(0)
    pat = rng.normal(size=(3, 3))                                   # 原胞里的任意位移图样（3 原子 × 3 方向）
    gam = np.zeros((9, 3))
    zb = np.zeros((9, 3))
    for j, (kap, t) in enumerate(sc["origin"]):
        phi = 2 * math.pi * (t[0] + t[1]) / 3.0                       # K·t
        gam[j] = pat[kap]
        zb[j] = pat[kap] * math.cos(phi)
    gam /= np.linalg.norm(gam)
    zb /= np.linalg.norm(zb)
    w = F.gamma_weight(np.array([gam, zb]), perms)
    assert w[0] == pytest.approx(1.0, abs=1e-12) and w[1] == pytest.approx(0.0, abs=1e-12)


def _toy(F, alpha, disp=0.01, shift=0.0):
    """返回 (D_z, D_xy, groups)：玩具模型走一遍完整流程。shift 加在超胞能级上（测对齐）。"""
    M, sc = _supercell(F)
    inv3 = np.linalg.inv(sc["M3"])
    frac = sc["frac"]

    def idx_of(t):
        fs = np.array([t[0], t[1], 0.0]) @ inv3
        d = frac - (np.array([fs[0], fs[1], 0.5]))
        d[:, :2] -= np.rint(d[:, :2])
        return int(np.where(np.all(np.abs(d) < 1e-6, axis=1))[0][0])

    H0 = np.zeros((3, 3))
    for i, (_k, t) in enumerate(sc["origin"]):
        for dl in NEIGH:
            H0[i, idx_of((t[0] + dl[0], t[1] + dl[1]))] += T_HOP
    # 9 个正交归一的超胞模：每个方向 平移（Γ）+ cos/sin 驻波（K）；频率：平移 0，面内驻波 5 THz，面外驻波 10 THz
    vecs, freqs = [], []
    for a in range(3):
        for kind in ("u", "c", "s"):
            e = np.zeros((3, 3))
            for j, (_k, t) in enumerate(sc["origin"]):
                phi = 2 * math.pi * (t[0] + t[1]) / 3.0
                e[j, a] = {"u": 1 / math.sqrt(3), "c": math.sqrt(2 / 3) * math.cos(phi),
                           "s": math.sqrt(2 / 3) * math.sin(phi)}[kind]
            vecs.append(e)
            freqs.append(0.0 if kind == "u" else (10.0 if a == 2 else 5.0))
    vecs, freqs = np.array(vecs), np.array(freqs)
    perms = F.translation_perms(sc)
    modes, groups, notes = F.zone_boundary_groups(freqs, vecs, perms)
    assert not notes
    masses = [50.0] * 3

    def eig(Hd):
        return np.sort(np.linalg.eigvalsh(H0 + np.diag(Hd))) + shift

    e0 = eig(np.zeros(3))
    prim = {"k%d" % i: [_band(k), _band(k) + 10.0]                       # 第二条 = 截断处的高能带（不进窗口）
            for i, k in enumerate(F.folded_kpoints(M))}
    star = [n for n, k in zip(prim, F.folded_kpoints(M))
            if any(np.allclose(np.asarray(k[:2]) % 1, np.asarray(s[:2]) % 1) for s in K_STAR)]
    labels, sh = F.label_folded(e0, prim, edge_eV=-3 * T_HOP)
    assert sh == pytest.approx(shift, abs=1e-12)
    grp = F.edge_group(labels, e0, star, 0)
    by = {m["index"]: m for m in modes}
    out = {}
    for g in groups:
        s2 = {}
        for i in g["modes"]:
            delta, u = F.displacement(vecs[by[i]["_src"]], masses, disp)
            ep, em = eig(alpha * delta[:, 2]), eig(-alpha * delta[:, 2])
            s2[i] = F.split_sq(e0, ep, em, grp, u)["S2"]
        out[round(g["freq_THz"], 3)] = F.group_D(s2, len(grp))
    return out[10.0], out[5.0], groups


@pytest.mark.parametrize("alpha", [0.8, 2.5])
def test_toy_model_recovers_coupling(F, alpha):
    Dz, Dxy, groups = _toy(F, alpha)
    assert Dz == pytest.approx(alpha, rel=1e-3)
    assert Dxy == pytest.approx(0.0, abs=1e-6)
    assert sorted(len(g["modes"]) for g in groups) == [2, 4]


def test_toy_model_independent_of_reference_shift(F):
    Dz, _, _ = _toy(F, 1.3, shift=0.37)
    assert Dz == pytest.approx(1.3, rel=1e-3)


def test_label_errors(F):
    M = F.fold_matrix(K_STAR, LAT[:2, :2])
    prim = {"k%d" % i: [_band(k), _band(k) + 3.0] for i, k in enumerate(F.folded_kpoints(M))}
    e_sc = sorted(x for v in prim.values() for x in v)
    e_bad = np.array(e_sc) + np.array([0, 0.2, 0, 0, 0, 0])               # 窗口里一条对不上
    with pytest.raises(ValueError, match="对不上"):
        F.label_folded(e_bad, prim, edge_eV=-3 * T_HOP, window=2.0)
    with pytest.raises(ValueError, match="NBANDS"):                         # 窗口顶端高过某折叠 k 的最高带
        F.label_folded(e_sc, prim, edge_eV=-3 * T_HOP, window=3.5)
    labels, _ = F.label_folded(e_sc, prim, edge_eV=-3 * T_HOP, window=2.0)
    names = [n for n, k in zip(prim, F.folded_kpoints(M)) if not np.allclose(k, 0)]
    assert len(F.edge_group(labels, e_sc, names, 0)) == 2
    prim2 = dict(prim)
    prim2["k0"] = [-3 * T_HOP + 0.004, 9.0]                                # Γ 的一条落在 K 二重态旁 4 meV
    e2 = sorted(x for v in prim2.values() for x in v)
    lab2, _ = F.label_folded(e2, prim2, edge_eV=-3 * T_HOP, window=2.0)
    with pytest.raises(ValueError):
        F.edge_group(lab2, e2, names, 0)


def test_inter_channels(F):
    params = {"schema": F.SCHEMA, "carriers": {"electron": {"status": "ok", "valleys": [
        {"label": "K", "groups": [{"D_eV_A": 2.0, "hw_meV": 46.0}, {"D_eV_A": 0.01, "hw_meV": 30.0},
                                  {"D_eV_A": 1.0, "hw_meV": 46.0004}]}]}}}
    chs, dropped, _ = F.inter_channels(params, "electron")
    assert len(chs) == 1 and chs[0]["kind"] == "inter_star" and chs[0]["from"] == chs[0]["to"] == "K"
    assert chs[0]["D_eV_A"] == pytest.approx(math.sqrt(5.0))
    assert len(dropped) == 1
    assert F.inter_channels(params, "hole") == ([], [], [])
    assert F.params_digest(params, "electron") and F.params_digest(params, "hole") is None
