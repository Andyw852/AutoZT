"""[0033] 光学形变势（ODP，谷内）与谷间散射（IVS）：常数形变势模型 dp_optical.py、S8.4 gen 的 step.conf
接线（DP_OPTICAL_E/H、谷表、ADP 谷内掩码）、插件/IR_FIX 的源码级接线。纯 numpy，不需要 amset。
真 AMSET 端到端（合成双谷模型）在 skill/ke-dft-cpu/step8.4_amset2d/test_dp_optical_e2e.py（需 amset 环境）。"""
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
S84 = KE / "step8.4_amset2d"


def _load(name, path):
    for d in (ROOT, KE, S84, ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def dpo():
    return _load("dp_optical_0033", S84 / "dp_optical.py")


# ---------------------------------------------------------------- 解析
def test_parse_channels_kinds_and_aliases(dpo):
    ch = dpo.parse_channels("K>K:4.1@50; K>K':2.6@46; K>Q:2.0@48\nG>K:1.5@30; Γ>Γ:3@40")
    assert [(c["from"], c["to"], c["kind"]) for c in ch] == [
        ("K", "K", "intra"), ("K", "K", "inter_star"), ("K", "Q", "inter"), ("Γ", "K", "inter"), ("Γ", "Γ", "intra")]
    assert ch[0]["D_eV_A"] == 4.1 and ch[0]["hw_meV"] == 50.0
    assert dpo.parse_channels("") == [] and dpo.parse_channels(None) == []
    assert dpo.parse_channels("k>k′:1@2")[0]["kind"] == "inter_star"      # 小写、′ 也认
    assert len(dpo.parse_channels("K>Q:1.2@24; K>Q:2.0@48")) == 2            # 同一对谷按声子支分列


@pytest.mark.parametrize("bad", ["K>Q:2.0", "K>Q:-1@40", "K>Q:1@0", "K>Q':1@40", "K>K:1@40; K>K:2@40", "KQ:1@4"])
def test_parse_channels_rejects(dpo, bad):
    with pytest.raises(ValueError):
        dpo.parse_channels(bad)


def test_parse_valleys_and_dpt_table(dpo):
    v = dpo.parse_valleys("K:0.333333,0.333333|-0.333333,-0.333333,0; Q:0.1667,0.1667,0")
    assert [(x["label"], x["g"]) for x in v] == [("K", 2), ("Q", 1)] and v[0]["k_frac"][0] == [0.333333, 0.333333, 0.0]
    with pytest.raises(ValueError):
        dpo.parse_valleys("K:0.3")
    dpt = {"results": [{"carrier": "electron", "multi_valley": {"status": "refused", "valleys": [
        {"label": "K", "g": 2, "dE_eV": 0.0, "k_frac": [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]},
        {"label": "Q", "g": 6, "dE_eV": 0.05, "k_frac": [[0.17, 0.17, 0]] * 6}]}},
        {"carrier": "hole", "multi_valley": {}}]}
    ve = dpo.valleys_from_dpt(dpt, "electron")
    assert [(x["label"], x["g"]) for x in ve] == [("K", 2), ("Q", 6)]          # 拒绝出数时谷表照样能用
    assert dpo.valleys_from_dpt(dpt, "hole") == []
    dpo.check_channels(dpo.parse_channels("K>Q:1@40; Q>Q':1@40"), ve)
    with pytest.raises(ValueError):
        dpo.check_channels(dpo.parse_channels("K>M:1@40"), ve)
    with pytest.raises(ValueError):
        dpo.check_channels(dpo.parse_channels("G>G':1@40"), [{"label": "Γ", "g": 1, "k_frac": [[0, 0, 0]]}])


# ---------------------------------------------------------------- 谷归属
def _hex_rlat(a=3.18, c=20.0):
    L = np.array([[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, c]])
    return 2 * np.pi * np.linalg.inv(L).T


def test_assign_members_hexagonal_K_Kprime(dpo):
    r = _hex_rlat()
    cen = np.array([[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]])
    k = np.array([[1 / 3, 1 / 3, 0], [2 / 3, 2 / 3, 0], [0.30, 0.36, 0.2], [1 / 3 + 1, 1 / 3 - 1, 0], [0.62, 0.70, 0]])
    m = dpo.assign_members(k, cen, r)
    assert m.tolist() == [0, 1, 0, 0, 1]         # 2/3,2/3 ≡ −1/3,−1/3（K′）；周期像、kz 无关
    keep = dpo.same_member([1 / 3, 1 / 3, 0], k, cen, r)
    assert keep.tolist() == [True, False, True, True, False]


# ---------------------------------------------------------------- 解析对照（双星、抛物带）
A_LAT = 3.18
M_D = 0.5
CELL_KG = 160.0 * 1.66053906660e-27
CH = "A>A:4.0@50; A>A':2.5@45; A>B:2.0@40; B>B:3.0@50; B>A:2.0@40"
VAL = [{"label": "A", "g": 2, "k_frac": [[0.25, 0.25, 0.0], [-0.25, -0.25, 0.0]]},
       {"label": "B", "g": 1, "k_frac": [[0.5, 0.0, 0.0]]}]
DE_B = 0.080


@pytest.fixture(scope="module")
def two_star(dpo):
    """方格子；星 A（2 个谷）与星 B（1 个谷，高 80 meV）都是各向同性抛物带 m=0.5；能带 = 各谷抛物面的下包络。"""
    n = 360
    rlat = 2 * np.pi / A_LAT * np.eye(3)
    rlat[2, 2] = 2 * np.pi / 20.0
    g = np.arange(n) / n
    k1, k2 = np.meshgrid(g, g, indexing="ij")
    kf = np.stack([k1.ravel(), k2.ravel(), np.zeros(n * n)], 1)
    hb2m = 3.80998212 / M_D                                       # ħ²/2m，eV·Å²
    E = np.full(len(kf), np.inf)
    for v, off in ((VAL[0], 0.0), (VAL[1], DE_B)):
        for c in v["k_frac"]:
            d = kf - np.array(c)
            d = d - np.rint(d)
            d[:, 2] = 0
            E = np.minimum(E, off + hb2m * np.sum((d @ rlat) ** 2, axis=1))
    E_H = (E / dpo.HA_EV)[None]
    rates = dpo.carrier_rates(E_H, kf, rlat, VAL, dpo.parse_channels(CH), [300.0], CELL_KG,
                              {"intra", "inter_star", "inter"}, +1)[0, 0]
    member = dpo.assign_members(kf, dpo.member_table(VAL)[0], rlat)
    return {"E": E, "rates": rates, "star": dpo.member_table(VAL)[1][member], "area": A_LAT ** 2 * 1e-20}


def _analytic(dpo, x, star, area):
    """该星里能量 x（eV，相对星 A 谷底）的态，各通道解析率之和。"""
    tot = 0.0
    lab = "AB"[star]
    for c in dpo.parse_channels(CH):
        if c["from"] != lab:
            continue
        tgt = c["to"]
        gt = {"intra": 1, "inter_star": 1, "inter": (2 if tgt == "A" else 1)}[c["kind"]]
        vmin = 0.0 if tgt == "A" else DE_B
        tot += dpo.analytic_rate(x * dpo.E_C, vmin * dpo.E_C, M_D, gt, area, c["D_eV_A"], c["hw_meV"], 300.0, CELL_KG)
    return tot


@pytest.mark.parametrize("star,x", [(0, 0.015), (0, 0.085), (0, 0.20), (1, 0.105), (1, 0.20)])
def test_rates_match_analytic_away_from_thresholds(dpo, two_star, star, x):
    sel = (two_star["star"] == star) & (np.abs(two_star["E"] - x) < 0.002)
    assert sel.sum() > 20
    num = float(np.mean(two_star["rates"][sel]))
    ref = _analytic(dpo, x, star, two_star["area"])
    assert ref > 1e11 and abs(num / ref - 1) < 0.06, (num, ref)


def test_thresholds_single_channel(dpo, two_star):
    """只开 A>B（ħω = 40 meV，B 谷底 80 meV）：x < 40 meV 吸收也到不了 -> 0；
    40–120 meV 只有吸收 N·G；> 120 meV 再加发射 (N+1)·G -> 两段之比 = (2N+1)/N。"""
    E, st = two_star["E"], two_star["star"]
    n = int(round(math.sqrt(len(E))))
    rlat = 2 * np.pi / A_LAT * np.eye(3)
    rlat[2, 2] = 2 * np.pi / 20.0
    g = np.arange(n) / n
    k1, k2 = np.meshgrid(g, g, indexing="ij")
    kf = np.stack([k1.ravel(), k2.ravel(), np.zeros(n * n)], 1)
    r = dpo.carrier_rates((E / dpo.HA_EV)[None], kf, rlat, VAL, dpo.parse_channels("A>B:2.0@40"), [300.0],
                          CELL_KG, {"inter"}, +1)[0, 0]
    at = lambda x: float(np.mean(r[(st == 0) & (np.abs(E - x) < 0.002)]))   # noqa: E731
    N = 1.0 / math.expm1(0.040 / (8.617333262e-5 * 300.0))
    assert at(0.015) < 0.01 * at(0.085)
    assert abs(at(0.160) / at(0.085) / ((2 * N + 1) / N) - 1) < 0.03
    assert abs(at(0.085) / _analytic_one(dpo, 0.085) - 1) < 0.05


def _analytic_one(dpo, x):
    return dpo.analytic_rate(x * dpo.E_C, DE_B * dpo.E_C, M_D, 1, A_LAT ** 2 * 1e-20, 2.0, 40, 300.0, CELL_KG)


def test_rate_floor_beats_amset_zero_rate_threshold(dpo):
    # AMSET 的 _interpolate_zero_rates 把 ≤ 1e6 s⁻¹ 当没算到、用最近邻顶替：解析率的真 0 不能落进去
    assert dpo.RATE_FLOOR > 1e6 and dpo.RATE_FLOOR < 1e8


def test_vacuum_cancels(dpo):
    """D²·G/M：原胞面积/质量不变，c 变化不影响（G 每原胞、M 每原胞，都与真空无关）。"""
    r1 = dpo.analytic_rate(0.1 * dpo.E_C, 0, 0.5, 1, 1e-19, 3.0, 50, 300, CELL_KG)
    assert r1 > 0 and dpo.analytic_rate(0.01 * dpo.E_C, 0, 0.5, 1, 1e-19, 3.0, 50, 300, CELL_KG) < r1


def test_valley_report(dpo, two_star):
    rlat = 2 * np.pi / A_LAT * np.eye(3)
    n = int(round(math.sqrt(len(two_star["E"]))))
    g = np.arange(n) / n
    k1, k2 = np.meshgrid(g, g, indexing="ij")
    kf = np.stack([k1.ravel(), k2.ravel(), np.zeros(n * n)], 1)
    E_H = (two_star["E"] / dpo.HA_EV)[None]
    ok = {r["label"]: r for r in dpo.valley_report(E_H, kf, rlat, VAL, +1)}
    assert ok["A"]["dk_rel"] < 0.01 and abs(ok["A"]["dE_eV"]) < 1e-6
    assert ok["B"]["dk_rel"] < 0.01 and abs(ok["B"]["dE_eV"] - DE_B) < 0.002
    wrong = [dict(VAL[0], k_frac=[[0.32, 0.18, 0], [-0.32, -0.18, 0]]), VAL[1]]   # 中心错位（换了晶胞设置）
    bad = {r["label"]: r for r in dpo.valley_report(E_H, kf, rlat, wrong, +1)}
    assert bad["A"]["dk_rel"] > 0.05


# ---------------------------------------------------------------- gen：step.conf -> 2d_correction.json
def _gen(tmp_path, monkeypatch, **conf):
    g = _load("s84_0033", S84 / "gen_step14_amset2d.py")
    for k, v in conf.items():
        monkeypatch.setattr(g, k, v)
    monkeypatch.setattr(g, "SCATTERING", list(conf.get("SCATTERING", ["ADP", "IMP", "POP"])))
    return g


def _dpt(tmp_path):
    d = tmp_path / "step8.2_dpt"
    d.mkdir(parents=True)
    d.joinpath("dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": "electron", "multi_valley": {"status": "ok", "valleys": [
            {"label": "K", "g": 2, "dE_eV": 0.0, "k_frac": [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]},
            {"label": "Q", "g": 6, "dE_eV": 0.06, "k_frac": [[0.17, 0.17, 0], [-0.17, -0.17, 0], [0.17, -0.34, 0],
                                                             [-0.17, 0.34, 0], [0.34, -0.17, 0], [-0.34, 0.17, 0]]}]}},
        {"carrier": "hole", "multi_valley": {"status": "ok", "valleys": [
            {"label": "K", "g": 2, "dE_eV": 0.0, "k_frac": [[1 / 3, 1 / 3, 0], [-1 / 3, -1 / 3, 0]]}]}}]}))


def test_gen_resolve_from_dpt(tmp_path, monkeypatch, capsys):
    _dpt(tmp_path)
    g = _gen(tmp_path, monkeypatch, DP_OPTICAL_E="K>K:4.1@50; K>K':2.6@46; K>Q:2.0@48; Q>Q:3@50; Q>K:2@48",
             DP_OPTICAL_H="K>K:4.0@48")
    rec = g.resolve_dp_optical(tmp_path)
    assert g.SCATTERING == ["ADP", "IMP", "POP", "ODP", "IVS"]
    assert rec["electron"]["valleys_source"] == "step8.2_dpt/dpt_result.json"
    assert [v["g"] for v in rec["electron"]["valleys"]] == [2, 6]
    assert rec["adp_intravalley"] == {"electron": True, "hole": False}       # auto：只有开了 IVS 的载流子
    assert "hole" in rec and rec["hole"]["channels"][0]["kind"] == "intra"
    json.dumps(rec)                                                          # 能落进 2d_correction.json


def test_gen_resolve_off_and_modes(tmp_path, monkeypatch):
    g = _gen(tmp_path, monkeypatch)
    assert g.resolve_dp_optical(tmp_path) is None and g.SCATTERING == ["ADP", "IMP", "POP"]   # 出厂：完全不变
    _dpt(tmp_path)
    g = _gen(tmp_path, monkeypatch, DP_OPTICAL_E="K>K:4@50", ADP_INTRAVALLEY="on")
    rec = g.resolve_dp_optical(tmp_path)
    assert g.SCATTERING[-1] == "ODP" and "IVS" not in g.SCATTERING
    assert rec["adp_intravalley"] == {"electron": True, "hole": False}
    g = _gen(tmp_path, monkeypatch, DP_OPTICAL_E="K>K':2@46", ADP_INTRAVALLEY="off")
    assert g.resolve_dp_optical(tmp_path)["adp_intravalley"]["electron"] is False


def test_gen_resolve_manual_valleys_without_s82(tmp_path, monkeypatch):
    g = _gen(tmp_path, monkeypatch, DP_OPTICAL_H="K>K:4@48; K>G:2@40",
             DP_OPTICAL_VALLEYS_H="K:0.3333,0.3333,0|-0.3333,-0.3333,0; G:0,0,0")
    rec = g.resolve_dp_optical(tmp_path)
    assert rec["hole"]["valleys_source"].startswith("step.conf") and [v["label"] for v in rec["hole"]["valleys"]] == ["K", "Γ"]


@pytest.mark.parametrize("conf", [
    {"DP_OPTICAL_E": "K>K:4@50"},                                       # 没跑 S8.2、也没手填谷表
    {"DP_OPTICAL_E": "K>M:4@50", "DP_OPTICAL_VALLEYS_E": "K:0.33,0.33,0"},   # 通道指向谷表里没有的谷
    {"DP_OPTICAL_E": "K>K:4@"},                                         # 写法错
    {"SCATTERING": ["ADP", "POP", "IVS"]},                              # 手写 IVS 却没有通道
    {"DP_OPTICAL_E": "K>K:4@50", "DP_OPTICAL_VALLEYS_E": "K:0.33,0.33,0", "ADP_INTRAVALLEY": "maybe"},
])
def test_gen_resolve_errors_are_loud(tmp_path, monkeypatch, conf):
    g = _gen(tmp_path, monkeypatch, **conf)
    with pytest.raises(SystemExit):
        g.resolve_dp_optical(tmp_path)


def test_gen_wiring_source():
    src = (S84 / "gen_step14_amset2d.py").read_text(encoding="utf-8")
    for k in ("DP_OPTICAL_E", "DP_OPTICAL_H", "DP_OPTICAL_VALLEYS_E", "DP_OPTICAL_VALLEYS_H", "ADP_INTRAVALLEY"):
        assert '"%s": (' % k in src                                    # SPEC 声明（否则 step.conf 写了不认）
        assert '_p["%s"]' % k in src                                   # main() 真的读了
    i = src.index("_DP_OPTICAL_REC = resolve_dp_optical(cwd)")
    assert i < src.index("is_2d, elastic, c_len = apply_2d_corrections(cwd, elastic, eps_inf, eps_static)")
    assert '"dp_optical": _DP_OPTICAL_REC' in src and "_install_dp_optical(out)" in src
    assert "--drop IMP,POP,PIE,ODP,IVS" in src                         # 只 ADP 也剔除 ODP/IVS
    assert '{"ADP", "IMP", "POP", "PIE", "ODP", "IVS"}' in src
    yml = (KE / "skill.yaml").read_text(encoding="utf-8")
    blk = yml[yml.index("gen: gen_step14_amset2d.py"):]
    assert "dp_optical.py" in blk[:blk.index("fetch_all")]


# ---------------------------------------------------------------- 插件 / IR_FIX 接线（源码级 + 纯 numpy）
def test_plugin_wiring_source():
    src = (S84 / "amset2d_plugin.py").read_text(encoding="utf-8")
    assert "class OpticalDeformation2D(_DPOptical2D)" in src and "class Intervalley2D(_DPOptical2D)" in src
    assert '"ODP": OpticalDeformation2D, "IVS": Intervalley2D' in src
    assert "_ODP_NORM" not in src                                       # 旧的弹性近似 ODP 已删
    i = src.index("class AcousticDeformation2D(_ADP)")
    adp = src[i:src.index("# ---------------- ODP / IVS", i)]
    for s in ("def to_reference(self)", "def from_reference(cls, *ref)", "def _iv_keep(", "def _factor_core("):
        assert s in adp                                                 # 掩码数据带进 spawn 子进程
    assert "dpo.RATE_FLOOR" in src and "np.broadcast_to" in src
    # FD 窗口外：AMSET 用全体态的对数平均填充；只给一种载流子写通道时另一种的下限率会把它拉到 ~1e9
    # （电子 IVS 单机制迁移率放大 ~40 倍）-> 填完后改回解析率
    assert "def _install_fill_keep_dpo()" in src and "_DPO_FULL[cls.name] = rates" in src
    i = src.index('if {"ODP", "IVS"} & set(_ACTIVE):')
    assert "_install_fill_keep_dpo()" in src[i:i + 80]
    assert "def _valley_check(amset_data)" in src and "_valley_check(amset_data)" in src[src.index("def _dpo_rates"):]


def test_ir_fix_star_average_masks_after_averaging():
    irf = _load("ir_fix_0033", S84 / "amset_ir_fix.py")

    class Interp:
        def interpolate(self, spin, b, k):
            return np.eye(3)[None]

    class S:
        name = "ADP"
        deformation_potential = Interp()

        def _factor_core(self, uq, nq, spin, b, kpoint, v):
            return np.ones((1, 1, len(nq))) * (1.0 + uq[:, 0] ** 2)[None, None]

        def factor(self, *a):
            raise AssertionError("星平均应调 _factor_core，不调带掩码的 factor")

        def _iv_keep(self, uq, nq, spin, b, kpoint):
            self.seen = np.array(uq)
            return np.array([True, False])

    s = S()
    R90 = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert irf.star_average_factor(s, [np.eye(3), R90], None)
    uq = np.array([[1.0, 0, 0], [0, 1.0, 0]])
    out = s.factor(uq, np.array([0.01, 0.01]), None, 0, np.zeros(3), np.zeros(3))
    assert np.allclose(out[0, 0], [1.5, 0.0])                          # 平均 (2+1)/2，第二个 q 被掩掉
    assert np.allclose(s.seen, uq)                                     # 掩码用的是未旋转的 q
