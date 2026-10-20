"""[0059] QE DFPT 谷间零阶形变势（eph-qe-cpu 可选组 ivs_dfpt）+ ke-dft-cpu S8/S8.4 的 IVS_DFPT。

不需要 QE：QE 的真实产物（Si，QE 7.4，SG15 模守恒、30 Ry、6³ k）放在 tests/fixtures/ivs_dfpt_si/，
  q-01 = g 过程（q = 0.3X），q-02 = f 过程。期望值即 M0 核实时的数（与 QE 7.0 逐位相同）：
  g：只有 LO（63.06 meV）耦合，D = 4.79 eV/Å；f：48.48 meV 2.73、58.40 meV 4.42，TA/LA 为零（零阶选择定则）。
核对：对称性（Si 48 / MoS₂ 12 个操作）、星展开、谷对等价类与 q 分组、AHC 二进制布局与大小核对、规范不变、
星通道折算（D²/g_Y、D²/(g_X−1)）、消费端 API；四个生成器在临时材料目录里真跑；作业内小工具；
S8（3D）/ S8.4（2D）读取与优先级（DFPT > 冻声子 K↔K′ > 手填）、S8.3 来源标注、技能开关组与集群注入。
"""
import copy
import importlib.util
import io
import json
import math
import shutil
import subprocess
import sys
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
IVSDIR = ROOT / "skill" / "_common" / "ivs"
EPH = ROOT / "skill" / "eph-qe-cpu"
KE = ROOT / "skill" / "ke-dft-cpu"
FIX = ROOT / "tests" / "fixtures" / "ivs_dfpt_si"
for _p in (str(IVSDIR), str(ROOT / "skill" / "_common" / "opt"), str(EPH), str(KE), str(KE / "step8.4_amset2d")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ivs_dfpt as iv  # noqa: E402

A_SI = 5.431
SI_LAT = np.array([[0, A_SI / 2, A_SI / 2], [A_SI / 2, 0, A_SI / 2], [A_SI / 2, A_SI / 2, 0]])
SI_FRAC = np.array([[0, 0, 0], [0.25, 0.25, 0.25]])
A_M, C_M = 3.18, 18.0
MOS2_LAT = np.array([[A_M, 0, 0], [-A_M / 2, A_M * math.sqrt(3) / 2, 0], [0, 0, C_M]])
MOS2_FRAC = np.array([[1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5 + 1.56 / C_M], [2 / 3, 1 / 3, 0.5 - 1.56 / C_M]])
MOS2_POSCAR = ("MoS2\n1.0\n%.6f 0 0\n%.6f %.6f 0\n0 0 %.1f\nMo S\n1 2\nDirect\n" % (A_M, -A_M / 2, A_M * math.sqrt(3) / 2, C_M)
               + "".join("%.10f %.10f %.10f\n" % tuple(f) for f in MOS2_FRAC))
K_STAR = [[1 / 3, 1 / 3, 0.0], [-1 / 3, -1 / 3, 0.0]]
QV = 0.515 / 3


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rops(lat, frac, sym):
    return iv.recip_ops(iv.point_group(lat, frac, sym))


def _q_star(rops):
    return iv.expand_star([[QV, QV, 0.0]], rops)[0]


# ---------------------------------------------------------------------------
# 核心：对称性、星、谷对
# ---------------------------------------------------------------------------
def test_point_group_and_stars():
    r_si = _rops(SI_LAT, SI_FRAC, ["Si", "Si"])
    assert len(r_si) == 48
    ks, norb, added = iv.expand_star([[0, 0.425, 0.425]], r_si, full3d=True)
    assert (len(ks), norb, added) == (6, 1, True)
    assert iv.detect_dim(SI_LAT, SI_FRAC) == "3d"
    r_m = _rops(MOS2_LAT, MOS2_FRAC, ["Mo", "S", "S"])
    assert len(r_m) == 12 and iv.detect_dim(MOS2_LAT, MOS2_FRAC) == "2d"
    assert len(iv.expand_star([K_STAR[0]], r_m)[0]) == 2
    assert len(_q_star(r_m)) == 6
    _ks, norb, _a = iv.expand_star([K_STAR[0], [QV, QV, 0]], r_m)
    assert norb == 2                                        # K 与 Q 不是同一个星
    full, norb, added = iv.expand_star(K_STAR, r_m)           # 已完整：不补
    assert (len(full), norb, added) == (2, 1, False)
    assert iv.auto_kgrid(SI_LAT, 0.2, "3d") == [11, 11, 11] and iv.auto_kgrid(MOS2_LAT, 0.2, "2d") == [12, 12, 1]


def test_plan_pairs_si_and_mos2():
    r_si = _rops(SI_LAT, SI_FRAC, ["Si", "Si"])
    X = iv.expand_star([[0, 0.425, 0.425]], r_si, True)[0]
    cl, runs = iv.plan_pairs({"electron": [{"label": "X", "k": X}]}, r_si)
    assert sorted(c["n_pairs"] for c in cl) == [1, 4] and len(runs) == 2       # g（对面 1 个）、f（垂直 4 个）
    assert all(c["kind"] == "inter_star" for c in cl)
    r_m = _rops(MOS2_LAT, MOS2_FRAC, ["Mo", "S", "S"])
    K = iv.expand_star([K_STAR[0]], r_m)[0]
    stars = [{"label": "K", "k": K}, {"label": "Q", "k": _q_star(r_m)}]
    cl, runs = iv.plan_pairs({"electron": stars}, r_m)
    assert len(cl) == 6 and len(runs) == 6                   # K↔K′、K↔Q 近/远、Q↔Q′ 三种
    kq = [c for c in cl if c["kind"] == "inter"]
    assert len(kq) == 2
    for c in kq:                                             # K→Q 与 Q→K（反向过程）归同一类
        assert {p["from"] for p in c["pairs"]} == {"K", "Q"}
        assert sum(p["from"] == "K" for p in c["pairs"]) == 3 and sum(p["from"] == "Q" for p in c["pairs"]) == 1
    qq = sorted(c["n_pairs"] for c in cl if c["from"] == "Q" and c["to"] == "Q")
    assert qq == [1, 2, 2]
    for c in cl:                                             # k_fin = k_init + q（模倒格矢）
        assert iv.same_q(np.asarray(c["k_init"]) + np.asarray(c["q"]), c["k_fin"])
        r = runs[c["run"]]
        assert iv.same_q(r["k"][c["ik"]], c["k_init"]) and iv.same_q(r["q"], c["q"])


# ---------------------------------------------------------------------------
# 核心：AHC 二进制、动力学矩阵、D
# ---------------------------------------------------------------------------
def _write_ahc(d, gkk, etk, etq):
    d.mkdir(parents=True, exist_ok=True)
    np.asarray(gkk, "<c16").tofile(str(d / "ahc_gkk_iq1.bin"))
    (np.asarray(etk) / iv.RY_EV).astype("<f8").tofile(str(d / "ahc_etk_iq1.bin"))
    (np.asarray(etq) / iv.RY_EV).astype("<f8").tofile(str(d / "ahc_etq_iq1.bin"))


def test_read_ahc_layout_and_size(tmp_path):
    nk, nat, nahc, nbnd = 2, 1, 2, 3
    g = (np.arange(nk * 3 * nahc * nbnd) + 1j).reshape(nk, 3, nahc, nbnd)
    e = np.arange(nk * nbnd, dtype=float).reshape(nk, nbnd)
    _write_ahc(tmp_path / "a", g, e, e + 1)
    gg, ek, eq = iv.read_ahc(str(tmp_path / "a"), nbnd, nahc, nat, nk)
    # Fortran ahc_gkk(ib, jb, imode) 每个 k 一条记录：ib 变得最快 -> numpy (k, mode, jb, ib)
    assert np.array_equal(gg, g) and np.allclose(ek, e) and np.allclose(eq, e + 1)
    with pytest.raises(ValueError, match="应为"):
        iv.read_ahc(str(tmp_path / "a"), nbnd + 1, nahc, nat, nk)


def test_mode_D2_normalization_and_gauge_invariance():
    M = 50.0 * iv.AMU_RY
    e = np.eye(3)
    gk = np.zeros((3, 2, 4), complex)
    gk[0, 0, 2] = 0.01                                        # ⟨m=2|∂V/∂u_x|n=0⟩ = 0.01 Ry/bohr
    d2 = iv.mode_D2(gk, e, [M], [0], [2])
    assert d2[0] == pytest.approx((0.01 * iv.GKK_TO_EV_A) ** 2) and d2[1] == d2[2] == 0.0   # 单原子：D = |g|
    rng = np.random.default_rng(3)
    g2 = rng.normal(size=(6, 2, 5)) + 1j * rng.normal(size=(6, 2, 5))
    m2 = [28.0 * iv.AMU_RY, 16.0 * iv.AMU_RY]
    ee = np.linalg.qr(rng.normal(size=(6, 6)))[0]
    base = iv.mode_D2(g2, ee, m2, [0, 1], [3, 4])
    U = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0]
    V = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0]
    g3 = g2.copy()
    g3[:, :, 3:5] = np.einsum("jnm,mk->jnk", g2[:, :, 3:5], V)                     # 终态简并子空间换基
    g3 = np.einsum("jnm,nk->jkm", g3, U)                                           # 初态简并子空间换基
    assert np.allclose(iv.mode_D2(g3, ee, m2, [0, 1], [3, 4]), base)


def test_class_couplings_and_window():
    hw = np.array([10.0, 10.02, 30.0])
    e = np.eye(3)
    gk = np.zeros((3, 2, 5), complex)
    gk[:, 1, 3] = [0.01, 0.0, 0.02]
    etk = np.array([-1.0, 0.0, 1.0, 2.0, 3.0])
    res = iv.class_couplings(gk, etk, etk, hw, e, [iv.AMU_RY * 30], nocc=3, ahc_skip=2, ahc_nbnd=2, carrier="electron")
    assert res["init_bands"] == [4] and res["final_bands"] == [4]
    assert [g["modes"] for g in res["groups"]] == [[1, 2], [3]]                   # 10.00/10.02 meV 简并
    assert res["groups"][0]["D_eV_A"] == pytest.approx(0.01 * iv.GKK_TO_EV_A)
    with pytest.raises(ValueError, match="AHC 窗口"):
        iv.class_couplings(gk, np.array([0, 0, 0, 0, 3.0]), etk, hw, e, [iv.AMU_RY], 3, 3, 1, "hole")


def test_fold_channels():
    cls = [{"id": 0, "groups": [{"hw_meV": 40.0, "D_eV_A": 2.0}, {"hw_meV": 0.1, "D_eV_A": 9.0}],
            "pairs": [{"from": "K", "to": "Q", "kind": "inter", "g_to": 6, "g_from": 2}] * 3},
           {"id": 1, "groups": [{"hw_meV": 40.05, "D_eV_A": 1.0}],
            "pairs": [{"from": "K", "to": "Q", "kind": "inter", "g_to": 6, "g_from": 2}] * 3},
           {"id": 2, "groups": [{"hw_meV": 46.0, "D_eV_A": 3.0}],
            "pairs": [{"from": "K", "to": "K", "kind": "inter_star", "g_to": 2, "g_from": 2}]}]
    chs, dropped = iv.fold_channels(cls)
    kq = [c for c in chs if c["kind"] == "inter"]
    assert len(kq) == 1 and kq[0]["D_eV_A"] == pytest.approx(math.sqrt(3 * 4 / 6 + 3 * 1 / 6))    # 40/40.05 合并
    kk = [c for c in chs if c["kind"] == "inter_star"][0]
    assert kk["D_eV_A"] == pytest.approx(3.0) and kk["hw_meV"] == 46.0                       # g−1 = 1
    assert len(dropped) == 1 and dropped[0]["hw_meV"] == 0.1


def test_read_dyn_real_qe():
    qr = json.loads((FIX / "step12_ivs_ph" / "q-01" / "qrun.json").read_text())
    dyn = iv.read_dyn(str(FIX / "step12_ivs_ph" / "q-01" / "dyn"))
    assert dyn["nat"] == 2 and dyn["masses_ry"][0] == pytest.approx(28.085 * iv.AMU_RY, rel=1e-6)
    assert np.allclose(dyn["q_list"][0], qr["xq_2pi_alat"], atol=1e-6)
    hw, _e = iv.phonon_modes(dyn["C"][0], dyn["masses_ry"])
    assert np.allclose(hw, [11.35, 11.35, 18.79, 61.31, 61.31, 63.06], atol=0.02)


# ---------------------------------------------------------------------------
# 消费端 API
# ---------------------------------------------------------------------------
def _params(channels, valleys, dim="2d"):
    return {"schema": iv.SCHEMA, "dim": dim,
            "carriers": {"electron": {"status": "ok", "valleys": valleys, "channels": channels}}}


def test_consumer_api(tmp_path):
    vs = [{"label": "K", "g": 2, "k_frac": K_STAR}]
    p = _params([{"from": "K", "to": "K", "kind": "inter_star", "D_eV_A": 1.66, "hw_meV": 40.6},
                 {"from": "K", "to": "K", "kind": "inter_star", "D_eV_A": 0.01, "hw_meV": 20.0}], vs)
    chs, drop, notes = iv.inter_channels(p, "electron", keep={"KK"}, alias={"K": "KK"})
    assert [(c["from"], c["D_eV_A"], c["source"]) for c in chs] == [("KK", 1.66, "dfpt_ivs")] and len(drop) == 1
    assert iv.inter_channels(p, "electron", keep={"Z"})[2]                       # 不在谷表里 -> 提示
    assert iv.inter_channels(p, "hole") == ([], [], [])
    d1 = iv.params_digest(p, "electron")
    p2 = copy.deepcopy(p)
    p2["carriers"]["electron"]["channels"][0]["D_eV_A"] = 1.7
    assert d1 and iv.params_digest(p2, "electron") != d1 and iv.params_digest(p, "hole") is None
    assert "K>K'" in iv.describe(p, "electron")
    cur = [{"label": "Kx", "g": 2, "k_frac": [[0.334, 0.333, 0], [-0.333, -0.333, 0]]}]
    assert iv.match_stars(p, "electron", cur) == ({"K": "Kx"}, None)
    assert iv.match_stars(p, "electron", cur, tol=1e-5)[1]                       # 容差外 -> 对不上
    extra = cur + [{"label": "G", "g": 1, "k_frac": [[0, 0, 0]]}]
    assert "没有 DFPT 参数" in iv.match_stars(p, "electron", extra)[1]
    f = tmp_path / "p.json"
    assert iv.load_consumer(str(f), "auto", "2d") == (None, None)
    with pytest.raises(ValueError, match="不存在"):
        iv.load_consumer(str(f), "on", "2d")
    f.write_text(json.dumps(p))
    assert iv.load_consumer(str(f), "auto", "2d")[0]["schema"] == iv.SCHEMA
    got, msg = iv.load_consumer(str(f), "auto", "3d")
    assert got is None and "口径" in msg
    with pytest.raises(ValueError):
        iv.load_consumer(str(f), "on", "3d")
    with pytest.raises(ValueError, match="无效"):
        iv.load_consumer(str(f), "maybe", "2d")


def test_valley_sig_same_as_ke_common():
    import ke_common as kc
    raw = [{"label": "K", "k_frac": K_STAR, "k_grid": K_STAR},
           {"label": "Q", "k_frac": [[0.17, 0.17, 0], [-0.17, -0.17, 0]], "k_grid": [[1 / 6, 1 / 6, 0], [-1 / 6, -1 / 6, 0]]}]
    for g in (False, True):
        assert iv.valley_sig(raw, use_grid=g) == json.loads(json.dumps(kc._valley_sig(raw, g)))


# ---------------------------------------------------------------------------
# eph-qe-cpu 生成器（在临时"材料目录"里真跑，不需要 QE）
# ---------------------------------------------------------------------------
@pytest.fixture()
def eph(tmp_path, monkeypatch):
    for f in ("ivs_qe_common.py", "ivs_qe_job.py", "eph_common.py", "submit_ivs_qe.tpl"):
        shutil.copy(EPH / f, tmp_path)
    shutil.copy(IVSDIR / "ivs_dfpt.py", tmp_path)
    shutil.copy(ROOT / "skill" / "_common" / "opt" / "stepconf.py", tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.syspath_prepend(str(tmp_path))
    for m in ("ivs_qe_common", "stepconf", "eph_common"):
        sys.modules.pop(m, None)
    import ivs_qe_common as c
    return tmp_path, c


def _conf(d, extra=""):
    base = (EPH / "templates" / "step.conf").read_text(encoding="utf-8")
    (d / "step.conf").write_text(base + "\n" + extra, encoding="utf-8")


def test_gen_read_si_fixture(eph, capsys):
    d, c = eph
    for sub in ("step10_ivs_qlist", "step11_ivs_scf", "step12_ivs_ph"):
        shutil.copytree(FIX / sub, d / sub)
    assert c.main("read") == 0
    p = json.loads((d / "step13_ivs_read" / "ivs_dfpt_params.json").read_text())
    assert p["schema"] == iv.SCHEMA and p["dim"] == "3d" and p["provenance"]["qe_version"] == "7.4"
    raw = {r["n_pairs"]: r for r in p["raw"]}
    g = {round(x["hw_meV"], 1): x["D_eV_A"] for x in raw[1]["groups"]}          # g 过程
    assert g[63.1] == pytest.approx(4.79, abs=0.01) and g[11.3] < 1e-3 and g[18.8] < 1e-3 and g[61.3] < 1e-3
    f = {round(x["hw_meV"], 1): x["D_eV_A"] for x in raw[4]["groups"]}          # f 过程
    assert f[48.5] == pytest.approx(2.73, abs=0.01) and f[58.4] == pytest.approx(4.42, abs=0.01)
    assert f[17.2] < 1e-3 and f[22.3] < 1e-3
    chs = {round(x["hw_meV"], 1): x["D_eV_A"] for x in p["carriers"]["electron"]["channels"] if x["D_eV_A"] > 0.05}
    assert chs == pytest.approx({48.5: 2.443, 58.4: 3.953, 63.1: 2.141}, abs=2e-3)   # f：√(4/5)·D，g：√(1/5)·D
    assert all(x["kind"] == "inter_star" for x in p["carriers"]["electron"]["channels"])
    assert "X>X':2.141@63.06" in (d / "step13_ivs_read" / "ivs_dfpt_table.txt").read_text()
    shutil.rmtree(d / "step12_ivs_ph" / "q-02" / "ahc_dir")
    (d / "step12_ivs_ph" / "q-02" / "ahc_dir").mkdir()
    with pytest.raises(SystemExit, match="retry S12_ivsph"):
        c.main("read")


def test_gen_chain_si(eph):
    d, c = eph
    shutil.copy(FIX / "POSCAR", d / "POSCAR")
    _conf(d, "IVS_VALLEYS_E = X:0,0.425,0.425\nIVS_ECUTWFC = 30\nIVS_SCF_KGRID = 6 6 6\nIVS_INPUT_DFT = none\n"
             "IVS_PSEUDO_FILES = Si:Si_ONCV_PBE-1.2.upf\nIVS_NP = 8\nIVS_SBATCH = --partition=cpu192; --qos=regular\n"
             "QE_BIN = /opt/qe/bin\nQE_PSEUDO_DIR = /opt/pseudo\nQE_SETUP = module load qe/7.3\n"
             "IVS_ODP = off\n")
    assert c.main("qlist") == 0
    q1 = json.loads((d / "step10_ivs_qlist" / "ivs_qlist.json").read_text())
    ref = json.loads((FIX / "step10_ivs_qlist" / "ivs_qlist.json").read_text())
    assert q1["runs"] == ref["runs"] and q1["n_ops"] == 48 and q1["carriers"]["electron"]["valleys"][0]["g"] == 6
    assert c.main("scf") == 0
    s = d / "step11_ivs_scf"
    sub = (s / "submit.sh").read_text()
    for frag in ("#SBATCH --ntasks=8", "#SBATCH --partition=cpu192", "#SBATCH --qos=regular", "module load qe/7.3",
                 "export PATH=/opt/qe/bin:", "IVS_PSEUDO_DIR=/opt/pseudo", 'mpirun -np 8 "$@"',
                 "run_mpi pw.x -in scf.in", "ivs_qe_job.py check", "printf 'OK\\n' > ivs_scf.ok"):
        assert frag in sub, frag
    assert "{{" not in sub and (s / "ivs_qe_job.py").is_file()
    assert "K_POINTS automatic\n  6 6 6 0 0 0" in (s / "scf.in").read_text()
    assert "nbnd = __NBND__" in (s / "valley.in.tpl").read_text()
    assert len(json.loads((s / "valley_ref.json").read_text())["points"]) == 6
    with pytest.raises(SystemExit, match="S11_ivsscf"):
        c.main("ph")
    for f in ("ivs_scf_info.json", "ivs_scf.ok"):
        shutil.copy(FIX / "step11_ivs_scf" / f, s / f)
    assert c.main("ph") == 0
    for q in ("q-01", "q-02"):
        qr = json.loads((d / "step12_ivs_ph" / q / "qrun.json").read_text())
        fx = json.loads((FIX / "step12_ivs_ph" / q / "qrun.json").read_text())
        assert np.allclose(qr["xq_2pi_alat"], fx["xq_2pi_alat"]) and (qr["nbnd"], qr["ahc_nbndskip"], qr["ahc_nbnd"]) == (12, 2, 4)
        ahc = (d / "step12_ivs_ph" / q / "ph_ahc.in").read_text()
        assert "electron_phonon = 'ahc'" in ahc and "trans = .false." in ahc and "skip_upperfan = .true." in ahc
        assert "cp -r ../../step11_ivs_scf/tmp/si2.save tmp/" in (d / "step12_ivs_ph" / q / "submit.sh").read_text()
    assert json.loads((d / "step12_ivs_ph" / "q-01" / "qrun.json").read_text())["tr2_ph"] == 1e-16   # [0059-r1] 不再舍成 0
    assert iv.jsonable({"a": 1e-16, "b": 0.1 + 0.2, "c": -0.0}) == {"a": 1e-16, "b": 0.3, "c": 0.0}
    (d / "step12_ivs_ph" / "q-01" / "ivs_ph.ok").write_text("OK\n")
    qf = d / "step12_ivs_ph" / "q-01" / "qrun.json"                   # 旧版写的 tr2_ph = 0.0：不算输入变了
    legacy = json.loads(qf.read_text())
    legacy["tr2_ph"] = 0.0
    qf.write_text(json.dumps(legacy))
    marker = d / "step12_ivs_ph" / "q-01" / "keep.me"
    marker.write_text("x")
    buf = io.StringIO()
    with redirect_stdout(buf):
        c.main("ph")
    assert marker.is_file() and "不动" in buf.getvalue()                       # 已完成且输入没变：不碰
    _conf(d, "IVS_VALLEYS_E = X:0,0.425,0.425\nIVS_ECUTWFC = 30\nIVS_SCF_KGRID = 6 6 6\nIVS_AHC_NB = 3\n")
    with redirect_stdout(io.StringIO()):
        c.main("ph")
    assert not marker.exists() and list((d / "step12_ivs_ph").glob(".old_q-01_*"))   # 输入变了：旧目录移开


def _dpt_json(q_members, extra=None):
    kv = {"label": "K", "dE_eV": 0.0, "k_frac": K_STAR, "k_grid": K_STAR}
    qv = {"label": "Q", "dE_eV": 0.13, "k_frac": q_members, "k_grid": q_members}
    res = [{"carrier": "electron", "multi_valley": {"status": "rejected", "valleys": [kv, qv]}},
           {"carrier": "hole", "multi_valley": {"status": "ok", "valleys": [{"label": "Γ", "dE_eV": 0.0, "k_frac": [[0, 0, 0]],
                                                                             "k_grid": [[0, 0, 0]]}]}}]
    return {"results": res + (extra or [])}


def test_gen_qlist_2d_from_s82(eph):
    d, c = eph
    (d / "POSCAR").write_text(MOS2_POSCAR)
    r_m = _rops(MOS2_LAT, MOS2_FRAC, ["Mo", "S", "S"])
    Q = [list(map(float, k)) for k in _q_star(r_m)]
    (d / "ivs_inputs").mkdir()
    (d / "ivs_inputs" / "dpt_result.json").write_text(json.dumps(_dpt_json(Q)))
    (d / "ivs_inputs" / "workflow_method.txt").write_text("FUNC=pbesol\nDIM=2D\n")
    _conf(d, "IVS_ODP = off\n")
    assert c.main("qlist") == 0
    q1 = json.loads((d / "step10_ivs_qlist" / "ivs_qlist.json").read_text())
    assert q1["dim"] == "2d" and q1["qe"]["kgrid"] == [12, 12, 1] and q1["qe"]["input_dft"] == "pbesol"
    assert list(q1["carriers"]) == ["electron"]                                 # 空穴单谷：跳过
    e = q1["carriers"]["electron"]
    assert e["use_grid"] is True and len(q1["runs"]) == 6
    raw = _dpt_json(Q)["results"][0]["multi_valley"]["valleys"]
    assert e["valley_sig"] == iv.valley_sig(raw, use_grid=True)
    (d / "ivs_inputs" / "dpt_result.json").write_text(json.dumps(_dpt_json(Q[:3])))
    with pytest.raises(SystemExit, match="没有多谷的载流子"):              # auto：不完整的载流子跳过 -> 一个都不剩
        c.main("qlist")
    _conf(d, "IVS_CARRIERS = electron\n")
    with pytest.raises(SystemExit, match="谷表不完整"):                   # 指名要电子：报错
        c.main("qlist")
    (d / "ivs_inputs" / "workflow_method.txt").write_text("FUNC=pbe-d3\n")
    P = c.load_params()
    assert c.input_dft(P)[0] == "pbe" and "D3" in c.input_dft(P)[1]


def test_gen_qlist_odp_gamma(eph):
    """[0059-r3] IVS_ODP=auto：q 列表末尾固定加 q=Γ 的谷内 ODP 类（kind=intra），每星一个初态 k。"""
    d, c = eph
    shutil.copy(FIX / "POSCAR", d / "POSCAR")
    _conf(d, "IVS_VALLEYS_E = X:0,0.425,0.425\nIVS_ECUTWFC = 30\nIVS_SCF_KGRID = 6 6 6\nIVS_INPUT_DFT = none\n"
             "IVS_PSEUDO_FILES = Si:Si_ONCV_PBE-1.2.upf\nIVS_NP = 8\n"
             "QE_BIN = /opt/qe/bin\nQE_PSEUDO_DIR = /opt/pseudo\n")
    assert c.main("qlist") == 0
    q1 = json.loads((d / "step10_ivs_qlist" / "ivs_qlist.json").read_text())
    runs = q1["runs"]
    assert runs[-1]["q"] == [0.0, 0.0, 0.0] and runs[-1]["name"] == "q-odp"
    assert len(runs[-1]["k"]) == 1                       # Si 电子只有一个星 -> 一个谷内 k
    intra = [c for c in q1["classes"] if c["kind"] == "intra"]
    assert len(intra) == 1 and intra[0]["from"] == intra[0]["to"] == "X" and intra[0]["run"] == runs[-1]["id"]
    # off：不加 Γ 点
    _conf(d, "IVS_ODP = off\nIVS_VALLEYS_E = X:0,0.425,0.425\nIVS_ECUTWFC = 30\nIVS_SCF_KGRID = 6 6 6\nIVS_INPUT_DFT = none\n"
             "IVS_PSEUDO_FILES = Si:Si_ONCV_PBE-1.2.upf\n")
    assert c.main("qlist") == 0
    q2 = json.loads((d / "step10_ivs_qlist" / "ivs_qlist.json").read_text())
    assert len(q2["runs"]) == 2 and not any(c["kind"] == "intra" for c in q2["classes"])


def test_jzzn_env_preset(eph):
    d, c = eph
    P = c.load_params()
    P.update(QE_SETUP="jzzn")
    blk = c.env_block(P)
    assert "module load gcc/10.2.0" in blk and "run_mpi()" in blk and 'IVS_PSEUDO_DIR="${PSEUDO_DIR:-}"' in blk


# ---------------------------------------------------------------------------
# 作业内小工具（只用标准库）
# ---------------------------------------------------------------------------
def _job(args, cwd, env=None):
    import os
    e = dict(os.environ)
    e.update(env or {})
    return subprocess.run([sys.executable, str(EPH / "ivs_qe_job.py")] + args, cwd=str(cwd), env=e,
                          capture_output=True, text=True)


def test_job_env(tmp_path):
    b = tmp_path / "bin"
    b.mkdir()
    (b / "pw.x").write_text("#!/bin/sh\n")
    (b / "ph.x").write_bytes(b"\x7fELF" + b"\0" * 200000 + b"ahc_nbnd" + b"\0" * 10)
    for x in ("pw.x", "ph.x"):
        (b / x).chmod(0o755)
    ps = tmp_path / "ps"
    ps.mkdir()
    (ps / "Si.upf").write_text('<UPF version="2.0.1">\n<PP_HEADER pseudo_type="NC" />\n')
    (ps / "C.upf").write_text('<UPF version="2.0.1">\n<PP_HEADER pseudo_type="PAW" is_paw="T" />\n')
    import os
    env = {"PATH": str(b) + os.pathsep + os.environ["PATH"]}
    r = _job(["env", "--pseudo-dir", str(ps), "--pseudos", "Si:Si.upf"], tmp_path, env)
    rep = json.loads((tmp_path / "ivs_env.json").read_text())
    assert r.returncode == 0 and rep["ok"] and rep["ahc"] == "yes" and rep["pseudos"]["Si"]["kind"] == "NC"
    r = _job(["env", "--pseudo-dir", str(ps), "--pseudos", "Si:Si.upf C:C.upf Ge:Ge.upf"], tmp_path, env)
    rep = json.loads((tmp_path / "ivs_env.json").read_text())
    assert r.returncode == 2 and any("PAW" in p for p in rep["problems"]) and any("不存在" in p for p in rep["problems"])
    (b / "ph.x").write_bytes(b"\x7fELF" + b"\0" * 200000)
    r = _job(["env", "--pseudo-dir", str(ps), "--pseudos", "Si:Si.upf"], tmp_path, env)
    assert r.returncode == 2 and "electron_phonon='ahc'" in r.stderr


def _xml(path, eigs_ha, nelec=8.0):
    ks = "".join('<ks_energies><k_point weight="1">0 0 0</k_point><eigenvalues size="%d">%s</eigenvalues></ks_energies>'
                 % (len(e), " ".join(repr(x) for x in e)) for e in eigs_ha)
    path.write_text(
        '<?xml version="1.0"?><qes:espresso xmlns:qes="x"><output><atomic_structure nat="2" alat="7.25">'
        "<cell><a1>0 5 5</a1><a2>5 0 5</a2><a3>5 5 0</a3></cell></atomic_structure><band_structure>"
        "<lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit><nbnd>%d</nbnd><nelec>%s</nelec>%s"
        "</band_structure></output></qes:espresso>" % (len(eigs_ha[0]), nelec, ks), encoding="utf-8")


def test_job_valley_and_check(tmp_path):
    ha = 1 / 27.211386245988
    _xml(tmp_path / "scf.xml", [[0, 0, 0, 0.0]])
    (tmp_path / "valley.in.tpl").write_text("nbnd = __NBND__\n")
    r = _job(["valley", "--xml", "scf.xml", "--tpl", "valley.in.tpl", "--out", "valley.in", "--extra", "4"], tmp_path)
    assert r.returncode == 0 and (tmp_path / "valley.in").read_text() == "nbnd = 8\n"
    assert json.loads((tmp_path / "ivs_scf_info.json").read_text())["nocc"] == 4
    # 两个星：K（k 0、1）与 Q（k 2），导带（第 5 条，下标 4）K 1.0 eV、Q 1.13 eV
    _xml(tmp_path / "v.xml", [[0, 0, 0, 0, 1.0 * ha], [0, 0, 0, 0, 1.0 * ha], [0, 0, 0, 0, 1.13 * ha]])
    ref = {"points": [[0, 0, 0]] * 3,
           "carriers": {"electron": [{"label": "K", "points": [0, 1], "dE_ref_eV": 0.0},
                                     {"label": "Q", "points": [2], "dE_ref_eV": 0.12}]}}
    (tmp_path / "ref.json").write_text(json.dumps(ref))
    r = _job(["check", "--xml", "v.xml", "--ref", "ref.json", "--tol-mev", "30"], tmp_path)
    assert r.returncode == 0 and json.loads((tmp_path / "valley_check.json").read_text())["status"] == "ok"
    r = _job(["check", "--xml", "v.xml", "--ref", "ref.json", "--tol-mev", "5"], tmp_path)
    assert r.returncode == 3 and "IVS_VALLEY_CHECK" in r.stderr
    r = _job(["check", "--xml", "v.xml", "--ref", "ref.json", "--tol-mev", "5", "--mode", "warn"], tmp_path)
    assert r.returncode == 0 and json.loads((tmp_path / "valley_check.json").read_text())["status"] == "mismatch"


def test_job_ahc_size(tmp_path):
    src = FIX / "step12_ivs_ph" / "q-01" / "ahc_dir"
    ok = ["ahc", "--dir", str(src), "--nbnd", "12", "--ahc-nbnd", "4", "--nat", "2", "--nk", "1"]
    assert _job(ok, tmp_path).returncode == 0
    assert _job(ok[:-1] + ["2"], tmp_path).returncode == 3


# ---------------------------------------------------------------------------
# ke-dft-cpu S8（3D）
# ---------------------------------------------------------------------------
@pytest.fixture()
def si_params(eph):
    d, c = eph
    for sub in ("step10_ivs_qlist", "step11_ivs_scf", "step12_ivs_ph"):
        shutil.copytree(FIX / sub, d / sub)
    with redirect_stdout(io.StringIO()):
        c.main("read")
    return json.loads((d / "step13_ivs_read" / "ivs_dfpt_params.json").read_text())


@pytest.fixture()
def g10(tmp_path):
    g = _load("g10_ivsdfpt", KE / "step8_amset" / "gen_step10_amset.py")
    g.SCATTERING = ["ADP", "IMP", "POP"]
    g._SCATTERING_EXPLICIT = False
    for k in ("DP_OPTICAL_E", "DP_OPTICAL_H", "DP_OPTICAL_VALLEYS_E", "DP_OPTICAL_VALLEYS_H"):
        setattr(g, k, "")
    g.ADP_INTRAVALLEY, g.ODP_FROZEN, g.IVS_DFPT = "auto", "off", "auto"
    cwd = tmp_path / "ke"
    (cwd / "ivs_dfpt").mkdir(parents=True)
    return g, cwd


def _run3(g, cwd):
    buf = io.StringIO()
    with redirect_stdout(buf):
        rec = g.resolve_dp_optical_3d(False, cwd)
    return rec, buf.getvalue()


def test_s8_uses_dfpt(si_params, g10):
    g, cwd = g10
    (cwd / "ivs_dfpt" / "ivs_dfpt_params.json").write_text(json.dumps(si_params))
    assert "IVS_DFPT" in g.SPEC and "ivs_dfpt.py" in (KE / "skill.yaml").read_text()
    rec, log = _run3(g, cwd)
    e = rec["dp_optical"]["electron"]
    assert rec["mechanisms"] == ["IVS"] and g.SCATTERING == ["ADP", "IMP", "POP", "IVS"]
    assert sorted(round(c["D_eV_A"], 2) for c in e["channels"]) == [2.14, 2.44, 3.95]
    assert {c["source"] for c in e["channels"]} == {"dfpt_ivs"} and e["valleys"][0]["g"] == 6
    assert rec["dp_optical"]["adp_intravalley"]["electron"] is True
    dr = rec["dp_optical_ivs_dfpt"]
    assert dr["dfpt_use"] == {"electron": "used", "hole": "unavailable"} and len(dr["sha256"]) == 64
    assert "DFPT 谷间 IVS" in log
    json.dumps(rec, ensure_ascii=False)
    # 手填谷间通道只作对照；手填谷表（不同标签、中心在 0.0083 内）按 k 中心配对
    g.DP_OPTICAL_E = "Δ>Δ':11@62"
    g.DP_OPTICAL_VALLEYS_E = "|".join("%.4f,%.4f,%.4f" % tuple(np.asarray(k) * 0.4167 / 0.425)
                                      for k in si_params["carriers"]["electron"]["valleys"][0]["k_frac"])
    g.DP_OPTICAL_VALLEYS_E = "Δ:" + g.DP_OPTICAL_VALLEYS_E
    rec, log = _run3(g, cwd)
    e = rec["dp_optical"]["electron"]
    assert {c["from"] for c in e["channels"]} == {"Δ"} and len(e["channels"]) == 3
    assert e["ivs_dfpt"]["manual_reference"][0]["D_eV_A"] == 11 and e["ivs_dfpt"]["label_map"] == {"X": "Δ"}
    assert "只作对照" in log
    g.DP_OPTICAL_VALLEYS_E = "L:0.5,0.5,0.5|0.5,0,0|0,0.5,0|0,0,0.5"
    with pytest.raises(SystemExit, match="对不上"):
        _run3(g, cwd)


def test_s8_modes_and_guards(si_params, g10):
    g, cwd = g10
    f = cwd / "ivs_dfpt" / "ivs_dfpt_params.json"
    rec, _ = _run3(g, cwd)                                      # 没有文件：与以前一样
    assert rec["dp_optical"] is None and rec["dp_optical_ivs_dfpt"]["dfpt_use"]["electron"] == "unavailable"
    g.IVS_DFPT = "on"
    with pytest.raises(SystemExit, match="不存在"):
        _run3(g, cwd)
    f.write_text(json.dumps(si_params))
    g.IVS_DFPT = "off"
    rec, _ = _run3(g, cwd)
    assert rec["dp_optical"] is None and rec["dp_optical_ivs_dfpt"]["dfpt_use"]["electron"] == "off"
    g.IVS_DFPT = "auto"
    f.write_text(json.dumps(dict(si_params, dim="2d")))
    rec, log = _run3(g, cwd)
    assert rec["dp_optical"] is None and "口径" in log
    f.write_text(json.dumps(si_params))
    g.SCATTERING, g._SCATTERING_EXPLICIT = ["ADP", "IMP"], True
    rec, _ = _run3(g, cwd)
    assert rec["dp_optical"] is None and rec["dp_optical_ivs_dfpt"]["dfpt_use"]["electron"] == "not_in_scattering"
    g._SCATTERING_EXPLICIT = False
    g.SCATTERING = ["ADP", "IMP", "POP"]
    (cwd / "step3_uniform").mkdir()
    (cwd / "step3_uniform" / "INCAR").write_text("LSORBIT = .TRUE.\n")
    with pytest.raises(SystemExit, match="SOC"):
        _run3(g, cwd)


# ---------------------------------------------------------------------------
# ke-dft-cpu S8.4（2D）：DFPT > 冻声子 K↔K′ > 手填
# ---------------------------------------------------------------------------
def _mos2_params(r_m):
    Q = [list(map(float, k)) for k in _q_star(r_m)]
    chans = [{"from": "K", "to": "K", "kind": "inter_star", "D_eV_A": 1.66, "hw_meV": 40.6},
             {"from": "K", "to": "Q", "kind": "inter", "D_eV_A": 3.5, "hw_meV": 49.1},
             {"from": "Q", "to": "Q", "kind": "inter_star", "D_eV_A": 2.45, "hw_meV": 23.0},
             {"from": "Q", "to": "K", "kind": "inter", "D_eV_A": 3.5, "hw_meV": 49.1}]
    vs = [{"label": "K", "g": 2, "k_frac": K_STAR}, {"label": "Q", "g": 6, "k_frac": Q}]
    return _params(chans, vs, "2d"), Q


def test_s84_dfpt_priority(tmp_path, monkeypatch):
    import ivs_frozen as ivf
    import ke_common as kc
    r_m = _rops(MOS2_LAT, MOS2_FRAC, ["Mo", "S", "S"])
    params, Q = _mos2_params(r_m)
    mat = tmp_path / "ke"
    (mat / "step8.2_dpt").mkdir(parents=True)
    dpt = _dpt_json(Q)
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(dpt))
    raw = dpt["results"][0]["multi_valley"]["valleys"]
    frozen = {"schema": ivf.SCHEMA, "scope": "K↔K′",
              "carriers": {"electron": {"status": "ok", "valley_sig": {"use_grid": True, "sig": kc._valley_sig(raw, True)},
                                        "valleys": [{"label": "K", "groups": [{"D_eV_A": 1.2, "hw_meV": 41.0}]}]}}}
    (mat / "step8.22_ivs_read").mkdir()
    (mat / "step8.22_ivs_read" / "ivs_params.json").write_text(json.dumps(frozen))
    g14 = _load("g14_ivsdfpt", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")
    monkeypatch.setattr(g14, "SCATTERING", ["ADP", "IMP", "POP"])
    monkeypatch.setattr(g14, "DP_OPTICAL_E", "K>Q:9.3@30")
    with redirect_stdout(io.StringIO()):
        rec = g14.resolve_dp_optical(mat)                    # 还没有 DFPT：冻声子 K↔K′ + 手填 K>Q
    chs = rec["electron"]["channels"]
    assert sorted((c["from"], c["to"], c.get("source")) for c in chs) == [("K", "K", "frozen_ivs"), ("K", "Q", None)]
    assert g14._DP_OPTICAL_IVS_DFPT_REC["dfpt_use"]["electron"] == "unavailable"
    (mat / "ivs_dfpt").mkdir()
    (mat / "ivs_dfpt" / "ivs_dfpt_params.json").write_text(json.dumps(params))
    monkeypatch.setattr(g14, "SCATTERING", ["ADP", "IMP", "POP"])
    buf = io.StringIO()
    with redirect_stdout(buf):
        rec = g14.resolve_dp_optical(mat)
    chs = rec["electron"]["channels"]
    assert {c["source"] for c in chs} == {"dfpt_ivs"} and len(chs) == 4
    assert "IVS" in g14.SCATTERING and rec["adp_intravalley"]["electron"] is True
    info = rec["electron"]["ivs_dfpt"]
    assert info["frozen_reference"][0]["D_eV_A"] == 1.2 and info["manual_reference"][0]["D_eV_A"] == 9.3
    dr = g14._DP_OPTICAL_IVS_DFPT_REC
    assert dr["dfpt_use"]["electron"] == "used" and dr["digest"]["electron"] == iv.params_digest(params, "electron")
    ir = g14._DP_OPTICAL_IVS_REC
    assert ir["frozen_use"]["electron"] == "superseded_dfpt" and "electron" not in ir["digest"]
    assert kc.ivs_consumer_stale(ir, frozen, full={"dp_optical": rec}) == []        # 冻声子参数再变也不重排
    assert "只作对照" in buf.getvalue()
    # S8.2 谷表变了（Q 换位置）-> 报错要求重算
    Q2 = [list(map(float, k)) for k in iv.expand_star([[0.2, 0.2, 0]], r_m)[0]]
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(_dpt_json(Q2)))
    monkeypatch.setattr(g14, "IVS_FROZEN", "off")
    with pytest.raises(SystemExit, match="对不上"):
        with redirect_stdout(io.StringIO()):
            g14.resolve_dp_optical(mat)
    # S8.3 汇总的来源标注
    g13 = _load("g13_ivsdfpt", KE / "step8.3_output" / "gen_step13_output.py")
    (mat / "step8.4_amset2d").mkdir()
    (mat / "step8.4_amset2d" / "2d_correction.json").write_text(json.dumps(
        {"dp_optical": rec, "dp_optical_ivs_dfpt": dr}))
    lines = g13._ivs_source_lines(mat)
    assert len(lines) == 1 and "DFPT" in lines[0] and "electron" in lines[0]


# ---------------------------------------------------------------------------
# 技能开关组、加载器、集群注入
# ---------------------------------------------------------------------------
def _eph_skel():
    from autozt import apply_skills, load_config
    cfg, _ = load_config(None)
    return apply_skills(cfg)["task_types"]["eph-qe-cpu"]


def _expand(t0, flags):
    from autozt.bootstrap import expand_optional_steps
    t = copy.deepcopy(t0)
    t.update(flags)
    expand_optional_steps(t)
    return [s["label"] for s in t["steps"]]


def test_eph_switch_groups():
    t0 = _eph_skel()
    assert _expand(t0, {}) == ["S0_env", "S1_scf", "S2_wannier", "S3_ph", "S4_qe2pert", "S5_ephmat", "S6_summary"]
    assert _expand(t0, {"ivs_dfpt": True, "perturbo": False}) == ["S10_ivsq", "S11_ivsscf", "S12_ivsph", "S13_ivs"]
    assert len(_expand(t0, {"ivs_dfpt": True})) == 11


def test_loader_empty_main_steps(tmp_path):
    from autozt.bootstrap import _load_manifest
    p = tmp_path / "sk" / "skill.yaml"
    p.parent.mkdir()
    p.write_text("schema: 2\nname: sk\nsteps: []\noptional_steps:\n  a:\n    default: true\n    steps:\n"
                 "      - {seq: 1, name: step1_x, label: S1_x, check: marker}\n", encoding="utf-8")
    key, skel = _load_manifest(str(p))
    assert key == "sk" and skel["steps"] == []
    p.write_text("schema: 2\nname: sk\nsteps: []\n", encoding="utf-8")
    assert _load_manifest(str(p)) == (None, "没有 steps")


def test_cluster_injects_qe_keys(tmp_path, monkeypatch):
    import autozt
    from autozt import report
    y = tmp_path / "qe3090.yaml"
    y.write_text("name: qe3090\nssh_host: user_3090\nqe_bin: ~/software/AutoZT/qe-7.4/bin\n"
                 "qe_pseudo_dir: ~/software/AutoZT/pseudo-nc\nqe_setup: source ~/software/AutoZT/qe-7.4/env.sh\n")
    monkeypatch.setattr(autozt, "pkg_setting_path", lambda n: str(y) if n == "qe3090.yaml" else None)
    txt = report._cluster_conda_step_conf({"hpc_name": "qe3090"})
    assert "QE_BIN = ~/software/AutoZT/qe-7.4/bin" in txt and "QE_PSEUDO_DIR = ~/software/AutoZT/pseudo-nc" in txt
    assert "QE_SETUP = source ~/software/AutoZT/qe-7.4/env.sh" in txt
    import stepconf
    assert {"QE_BIN", "QE_PSEUDO_DIR", "QE_SETUP"} <= stepconf.RESERVED_PARAMS


# ---------------------------------------------------------------------------
# [0059-r3] DFPT 顺带算谷内 ODP（q=Γ）：单谷载流子、按 q 放宽的 AHC 窗口、Γ 声学/极性模、消费端不重复计数
# ---------------------------------------------------------------------------
SI_E = "IVS_VALLEYS_E = X:0,0.425,0.425\nIVS_ECUTWFC = 30\nIVS_SCF_KGRID = 6 6 6\nIVS_INPUT_DFT = none\n" \
       "IVS_PSEUDO_FILES = Si:Si_ONCV_PBE-1.2.upf\nQE_BIN = /opt/qe/bin\nQE_PSEUDO_DIR = /opt/pseudo\n"


def _qlist(d, c, extra):
    _conf(d, SI_E + extra)
    with redirect_stdout(io.StringIO()) as buf:
        assert c.main("qlist") == 0
    return json.loads((d / "step10_ivs_qlist" / "ivs_qlist.json").read_text()), buf.getvalue()


def test_r3_qlist_single_valley_hole(eph):
    d, c = eph
    shutil.copy(FIX / "POSCAR", d / "POSCAR")
    q1, log = _qlist(d, c, "IVS_VALLEYS_H = G:0,0,0\n")
    assert [r["name"] for r in q1["runs"]] == ["q-01", "q-02", "q-odp"]          # 谷间 q 的编号不变
    odp = q1["runs"][-1]
    assert odp["q"] == [0.0, 0.0, 0.0] and len(odp["k"]) == 2                     # 电子 X、空穴 Γ 各一个初态
    intra = sorted((x["carrier"], x["from"], x["ik"]) for x in q1["classes"] if x["kind"] == "intra")
    assert intra == [("electron", "X", 0), ("hole", "Γ", 1)]
    assert q1["carriers"]["hole"]["valleys"][0]["g"] == 1 and "只算谷内 ODP" in log
    # 关掉 ODP：单谷的空穴照旧跳过
    q1, log = _qlist(d, c, "IVS_VALLEYS_H = G:0,0,0\nIVS_ODP = off\n")
    assert "hole" not in q1["carriers"] and len(q1["runs"]) == 2 and "跳过" in log
    # 只要空穴：ODP 开着 = 只有 q-odp；关掉 = 没东西可算 -> 报错
    q1, _ = _qlist(d, c, "IVS_VALLEYS_H = G:0,0,0\nIVS_CARRIERS = hole\n")
    assert [r["name"] for r in q1["runs"]] == ["q-odp"] and len(q1["classes"]) == 1
    _conf(d, SI_E + "IVS_VALLEYS_H = G:0,0,0\nIVS_CARRIERS = hole\nIVS_ODP = off\n")
    with pytest.raises(SystemExit, match="多谷"):
        with redirect_stdout(io.StringIO()):
            c.main("qlist")
    # 直接带隙（电子、空穴都在 Γ）：q-odp 里同一个 k 只算一次
    _conf(d, SI_E.replace("X:0,0.425,0.425", "G:0,0,0") + "IVS_VALLEYS_H = G:0,0,0\n")
    with redirect_stdout(io.StringIO()):
        assert c.main("qlist") == 0
    q1 = json.loads((d / "step10_ivs_qlist" / "ivs_qlist.json").read_text())
    assert len(q1["runs"]) == 1 and q1["runs"][0]["k"] == [[0.0, 0.0, 0.0]]
    assert [x["ik"] for x in q1["classes"]] == [0, 0]


def _xml_k(path, pts, eigs_ha, nelec):
    """谷能带 XML，k_point 按 QE 的口径写成笛卡尔 2π/alat（晶胞与 _xml 相同）。"""
    A = np.array([[0, 5.0, 5.0], [5.0, 0, 5.0], [5.0, 5.0, 0]])
    B = np.linalg.inv(A).T * 7.25
    ks = "".join('<ks_energies><k_point weight="1">%s</k_point><eigenvalues size="%d">%s</eigenvalues></ks_energies>'
                 % (" ".join(repr(float(x)) for x in np.asarray(k, float) @ B), len(e), " ".join(repr(x) for x in e))
                 for k, e in zip(pts, eigs_ha))
    path.write_text(
        '<?xml version="1.0"?><qes:espresso xmlns:qes="x"><output><atomic_structure nat="2" alat="7.25">'
        "<cell><a1>0 5 5</a1><a2>5 0 5</a2><a3>5 5 0</a3></cell></atomic_structure><band_structure>"
        "<lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit><nbnd>%d</nbnd><nelec>%s</nelec>%s"
        "</band_structure></output></qes:espresso>" % (len(eigs_ha[0]), nelec, ks), encoding="utf-8")


def _s11(d, c, ha_by_point, nocc=4):
    """S11 的产物：gen_scf 写 valley_ref.json，再按点给谷能带（Ha）写 valley_bands.xml。"""
    with redirect_stdout(io.StringIO()):
        assert c.main("scf") == 0
    s11 = d / "step11_ivs_scf"
    shutil.copy(FIX / "step11_ivs_scf" / "ivs_scf_info.json", s11)
    (s11 / "ivs_scf.ok").write_text("OK\n")
    pts = json.loads((s11 / "valley_ref.json").read_text())["points"]
    _xml_k(s11 / "valley_bands.xml", pts, [ha_by_point(k) for k in pts], nelec=2.0 * nocc)
    return pts


def test_r3_ahc_window_per_q(eph):
    d, c = eph
    shutil.copy(FIX / "POSCAR", d / "POSCAR")
    _qlist(d, c, "IVS_VALLEYS_H = G:0,0,0\n")
    ha = 1 / 27.211386245988
    gam = [x * ha for x in (-12.0, 0.0, 0.0, 0.0, 2.5, 2.5, 2.5, 3.3)]           # Γ：价带顶三重简并
    xpt = [x * ha for x in (-8.0, -3.0, -2.0, -2.0, 0.6, 0.6, 3.0, 4.0)]         # X：导带底二重简并
    _s11(d, c, lambda k: gam if max(abs(x) for x in k) < 1e-9 else xpt)
    with redirect_stdout(io.StringIO()) as buf:
        assert c.main("ph") == 0
    qr = {n: json.loads((d / "step12_ivs_ph" / n / "qrun.json").read_text()) for n in ("q-01", "q-02", "q-odp")}
    assert [(qr[n]["ahc_nbndskip"], qr[n]["ahc_nbnd"]) for n in ("q-01", "q-02")] == [(2, 4), (2, 4)]   # 不变
    assert (qr["q-odp"]["ahc_nbndskip"], qr["q-odp"]["ahc_nbnd"]) == (1, 6)        # 空穴三重简并 -> 带边上下各 3 条
    assert "epsil" not in (d / "step12_ivs_ph" / "q-odp" / "ph_dfpt.in").read_text()  # 单元素：不算 Born 电荷
    assert "AHC 带 2..7" in buf.getvalue()
    # 旧版 qrun.json（没有 epsil 键）+ 同样的输入：已完成的 q 不动
    for n in ("q-01", "q-02"):
        f = d / "step12_ivs_ph" / n / "qrun.json"
        j = json.loads(f.read_text())
        j.pop("epsil")
        f.write_text(json.dumps(j))
        (d / "step12_ivs_ph" / n / "ivs_ph.ok").write_text("OK\n")
    with redirect_stdout(io.StringIO()) as buf:
        assert c.main("ph") == 0
    assert buf.getvalue().count("已完成且输入没变") == 2 and not list((d / "step12_ivs_ph").glob(".old_*"))
    # S11 是按旧谷表跑的（谷能带里没有空穴的 Γ）-> 报错要求先重跑 S11
    s11 = d / "step11_ivs_scf"
    ref = json.loads((s11 / "valley_ref.json").read_text())
    keep = [i for i, k in enumerate(ref["points"]) if max(abs(x) for x in k) > 1e-9]
    ref["points"] = [ref["points"][i] for i in keep]
    (s11 / "valley_ref.json").write_text(json.dumps(ref))
    _xml_k(s11 / "valley_bands.xml", ref["points"], [xpt] * len(keep), nelec=8.0)
    with pytest.raises(SystemExit, match="不在 S11 谷能带"):
        with redirect_stdout(io.StringIO()):
            c.main("ph")
    # S11 重新生成了（valley_ref.json 换了点）但谷能带还是旧的 -> 报错，不按位置错配
    pts2 = [[0.1, 0.0, 0.0]] + ref["points"][1:]
    _xml_k(s11 / "valley_bands.xml", pts2, [xpt] * len(keep), nelec=8.0)
    with pytest.raises(SystemExit, match="对不上"):
        with redirect_stdout(io.StringIO()):
            c.main("ph")
    # 没有谷能带（旧版 S11 的结果）：按 IVS_AHC_NB，提示 S13 可能报错
    (s11 / "valley_bands.xml").unlink()
    with redirect_stdout(io.StringIO()) as buf:
        assert c.main("ph") == 0
    assert "IVS_AHC_NB" in buf.getvalue()


def test_r3_epsil_for_compounds(eph):
    d, c = eph
    (d / "POSCAR").write_text(MOS2_POSCAR)
    r_m = _rops(MOS2_LAT, MOS2_FRAC, ["Mo", "S", "S"])
    (d / "ivs_inputs").mkdir()
    Q = [list(map(float, k)) for k in _q_star(r_m)]
    (d / "ivs_inputs" / "dpt_result.json").write_text(json.dumps(_dpt_json(Q)))
    (d / "ivs_inputs" / "workflow_method.txt").write_text("FUNC=pbesol\nDIM=2D\n")
    _conf(d)
    with redirect_stdout(io.StringIO()):
        assert c.main("qlist") == 0
    s11 = d / "step11_ivs_scf"
    s11.mkdir()
    (s11 / "scf.in").write_text("x\n")
    (s11 / "ivs_scf.ok").write_text("OK\n")
    (s11 / "ivs_scf_info.json").write_text(json.dumps({"alat_bohr": 6.0, "nocc": 13,
                                                        "cell_bohr": [[6.0, 0, 0], [-3.0, 5.196, 0], [0, 0, 34.0]]}))
    with redirect_stdout(io.StringIO()):
        assert c.main("ph") == 0
    ph = d / "step12_ivs_ph"
    assert "epsil = .true." in (ph / "q-odp" / "ph_dfpt.in").read_text()
    assert "epsil" not in (ph / "q-01" / "ph_dfpt.in").read_text()
    assert json.loads((ph / "q-odp" / "qrun.json").read_text())["epsil"] is True
    _conf(d, "IVS_ODP_POLAR = include\n")                                         # include：不算 Born 电荷
    with redirect_stdout(io.StringIO()) as buf:
        assert c.main("ph") == 0
    assert "epsil" not in (ph / "q-odp" / "ph_dfpt.in").read_text() and "q-odp 的输入变了" in buf.getvalue()


def _diatomic(k_xy=0.3, k_z=0.5, m=(28.0, 16.0), noise=0.0):
    """双原子（各向异性弹簧）Γ 点力常数（Ry/bohr²）+ 质量（Ry）：3 个声学模（ħω≈0）、面内 2 个 + 面外 1 个光学模。"""
    K = np.diag([k_xy, k_xy, k_z])
    C = np.block([[K, -K], [-K, K]]).astype(complex)
    C += noise * np.eye(6)                                                       # 数值噪声：声学模不严格为 0
    masses = np.array(m) * iv.AMU_RY
    return C, masses


def test_r3_gamma_acoustic_and_polar_modes():
    C, masses = _diatomic(noise=1e-7)
    hw, e = iv.phonon_modes(C, masses)
    ac, tw = iv.acoustic_modes(e, masses)
    assert ac == [0, 1, 2] and min(tw[ac]) > 0.99 and max(np.delete(tw, ac)) < 1e-6
    assert max(abs(hw[ac])) < 2.0 and min(hw[3:]) > 20.0                          # 声学模 ħω 不严格为 0，但仍是声学
    Z = np.array([np.eye(3) * 2.0, -np.eye(3) * 2.0])                            # 电中性：ΣZ* = 0
    z3 = iv.mode_charges(e, masses, Z)
    assert max(z3[ac]) < 1e-6 and min(np.delete(z3, ac)) > 1.0                   # 3D：光学模全是极性
    z2 = iv.mode_charges(e, masses, Z, two_d=True)
    zmode = int(np.argmax(hw))                                                    # k_z 最硬：最高的是面外光学模
    assert z2[zmode] < 1e-6 and min(z2[[3, 4]]) > 1.0                            # 2D 只看面内：面外模不算极性
    # dyn 文件里的介电张量与 Born 电荷（QE write_epsilon_and_zeu 的格式：第 α 行 = 电场方向）
    lines = (FIX / "step12_ivs_ph" / "q-01" / "dyn").read_text().splitlines()
    zz = np.array([[[1.9, 0.1, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.1]], [[-1.9, -0.1, 0.0], [0.0, -2.0, 0.0],
                                                                          [0.0, 0.0, -2.1]]])
    blk = ["", "     Dielectric Tensor:", ""] + ["%24.12f%24.12f%24.12f" % (12.0, 0.0, 0.0)] * 3
    blk += ["", "     Effective Charges E-U: Z_{alpha}{s,beta}", ""]
    for na in range(2):
        blk.append("     atom #    %d" % (na + 1))
        blk += ["%24.12f%24.12f%24.12f" % tuple(r) for r in zz[na]]
    cut = next(i for i, ln in enumerate(lines) if "Diagonalizing" in ln)           # 在频率表之前（dynmatrix 的顺序）
    return_dyn = lines[:cut] + blk + lines[cut:]
    import tempfile
    with tempfile.TemporaryDirectory() as t:
        f = Path(t) / "dyn"
        f.write_text("\n".join(return_dyn) + "\n")
        dyn = iv.read_dyn(str(f))
    assert np.allclose(dyn["zstar"], zz) and np.allclose(dyn["epsilon"], [[12.0, 0.0, 0.0]] * 3)
    assert len(dyn["C"]) == len(iv.read_dyn(str(FIX / "step12_ivs_ph" / "q-01" / "dyn"))["C"])
    assert iv.read_dyn(str(FIX / "step12_ivs_ph" / "q-01" / "dyn"))["zstar"] is None


def test_r3_exclude_and_fold_intra():
    hw = np.array([0.6, 0.6, 0.7, 40.0, 40.0, 62.0])
    gk = np.zeros((6, 2, 5), complex)
    gk[:, 1, 3] = [0.05, 0.05, 0.05, 0.0, 0.0, 0.02]                            # 声学模上有数值噪声级的耦合
    etk = np.array([-1.0, 0.0, 1.0, 2.0, 3.0])
    ex = {0: "声学模", 1: "声学模", 2: "声学模", 5: "极性模"}
    res = iv.class_couplings(gk, etk, etk, hw, np.eye(6), [iv.AMU_RY * 30] * 2, nocc=3, ahc_skip=2, ahc_nbnd=2,
                             carrier="electron", exclude=ex)
    assert [g.get("exclude") for g in res["groups"]] == ["声学模", "声学模", None, "极性模"]
    assert res["modes"][5]["exclude"] == "极性模"
    res.update(id=7, pairs=[{"from": "X", "to": "X", "kind": "intra", "g_to": 6, "g_from": 6}])
    chs, dropped = iv.fold_channels([res])
    assert [(c["kind"], c["hw_meV"]) for c in chs] == [("intra", 40.0)]          # 谷内：D 不除星成员数
    assert chs[0]["D_eV_A"] == pytest.approx(res["groups"][2]["D_eV_A"])
    assert sorted(d["why"] for d in dropped) == ["声学模", "声学模", "极性模"]


def test_r3_odp_exclusions(eph):
    d, c = eph
    C, masses = _diatomic(noise=1e-7)
    hw, e = iv.phonon_modes(C, masses)
    R = {"e": e, "masses": masses, "hw": hw, "zstar": None, "qr": {"name": "q-odp"}}
    P = c.load_params()
    q1 = {"dim": "3d", "structure": {"types": ["Si"]}}
    ex, rec = c.odp_exclusions(R, q1, P)
    assert sorted(ex) == [0, 1, 2] and "单元素" in rec["born_source"]
    q1 = {"dim": "3d", "structure": {"types": ["Ga", "As"]}}
    with pytest.raises(SystemExit, match="Born"):
        c.odp_exclusions(R, q1, P)
    R["zstar"] = np.array([np.eye(3) * 2.0, -np.eye(3) * 2.0])
    ex, rec = c.odp_exclusions(R, q1, P)
    assert sorted(ex) == [0, 1, 2, 3, 4, 5] and rec["polar_modes"] == [4, 5, 6]
    ex, rec = c.odp_exclusions(R, dict(q1, dim="2d"), P)
    assert sorted(ex) == [0, 1, 2, 3, 4] and rec["polar_modes"] == [4, 5]       # 面外光学模不算极性
    P["IVS_ODP_POLAR"] = "include"
    ex, rec = c.odp_exclusions(R, q1, P)
    assert sorted(ex) == [0, 1, 2] and rec["polar_modes"] == []


def _with_odp(si_params, hole_D=6.0, e_intra=True):
    """si_params（真实 QE 电子谷间）+ 假造的 DFPT 谷内结果：空穴 Γ（单谷）D=hole_D，电子 X 谷内 ≈ 0（对称性禁戒）。"""
    p = copy.deepcopy(si_params)
    p["carriers"]["hole"] = {"status": "ok", "valleys": [{"label": "Γ", "g": 1, "k_frac": [[0.0, 0.0, 0.0]]}],
                             "valleys_source": "step.conf IVS_VALLEYS_H", "use_grid": False,
                             "channels": [{"from": "Γ", "to": "Γ", "kind": "intra", "D_eV_A": hole_D, "hw_meV": 62.5}]}
    p["raw"].append({"carrier": "hole", "kind": "intra", "from": "Γ", "to": "Γ"})
    if e_intra:
        p["carriers"]["electron"]["channels"].append({"from": "X", "to": "X", "kind": "intra", "D_eV_A": 0.002,
                                                      "hw_meV": 63.0})
        p["raw"].append({"carrier": "electron", "kind": "intra", "from": "X", "to": "X"})
    return p


def _frozen_hole(of, D=7.77):
    vs = [{"label": "Γ", "k_frac": [[0, 0, 0]]}]
    return {"schema": of.SCHEMA, "dim": "3d",
            "carriers": {"hole": {"status": "ok", "valleys_source": "auto",
                                  "valley_sig": {"dim": "3d", "sig": of.valley_sig_3d(vs)},
                                  "valleys": [dict(v, groups=[{"group": 1, "modes": [4, 5, 6], "hw_meV": 63.0,
                                                               "D_eV_A": D}]) for v in vs]}}}


def test_r3_s8_dfpt_odp_supersedes_frozen(si_params, g10):
    import ke_common as kc
    import odp_frozen as of
    g, cwd = g10
    g.ODP_FROZEN, g.ODP_FROZEN_DMIN = "auto", 0.05
    (cwd / "step8.21_odp_read").mkdir()
    frozen = _frozen_hole(of)
    (cwd / "step8.21_odp_read" / "odp_params.json").write_text(json.dumps(frozen))
    f = cwd / "ivs_dfpt" / "ivs_dfpt_params.json"
    f.write_text(json.dumps(si_params))                     # DFPT 没算空穴（IVS_ODP=off 的旧参数）：空穴用冻声子
    rec, _ = _run3(g, cwd)
    h = rec["dp_optical"]["hole"]
    assert [(c["kind"], c["D_eV_A"], c["source"]) for c in h["channels"]] == [("intra", 7.77, "frozen")]
    assert rec["dp_optical_frozen"]["frozen_use"]["hole"] == "used"
    # DFPT 算了谷内：空穴只留 DFPT 的一条（不与冻声子叠加），冻声子只作对照、不再触发重排
    g.SCATTERING = ["ADP", "IMP", "POP"]
    f.write_text(json.dumps(_with_odp(si_params)))
    rec, log = _run3(g, cwd)
    h = rec["dp_optical"]["hole"]
    assert [(c["kind"], c["D_eV_A"], c["source"]) for c in h["channels"]] == [("intra", 6.0, "dfpt_odp")]
    assert h["ivs_dfpt"]["odp_reference"] == [{"valley": "Γ", "D_eV_A": 7.77, "hw_meV": 63.0, "source": "frozen"}]
    fr, dr = rec["dp_optical_frozen"], rec["dp_optical_ivs_dfpt"]
    assert fr["frozen_use"]["hole"] == "superseded_dfpt" and "hole" not in fr["digest"]
    assert kc.odp_consumer_stale(fr, _frozen_hole(of, D=9.0), full=rec) == []
    assert dr["odp_use"] == {"electron": "used", "hole": "used"} and dr["dfpt_use"]["hole"] == "used"
    e = rec["dp_optical"]["electron"]                        # 电子谷内 ≈ 0（< DMIN）：只有 3 条谷间
    assert all(c["kind"] != "intra" for c in e["channels"]) and len(e["channels"]) == 3
    assert rec["mechanisms"] == ["ODP", "IVS"] and g.SCATTERING == ["ADP", "IMP", "POP", "ODP", "IVS"]
    assert "只作对照" in log and "谷内 ODP 参数来自 DFPT" in log
    # 手填的谷内通道（ODP_FROZEN=off）同样只作对照
    g.SCATTERING, g.ODP_FROZEN = ["ADP", "IMP", "POP"], "off"
    g.DP_OPTICAL_H, g.DP_OPTICAL_VALLEYS_H = "G>G:9.9@60", "G:0,0,0"
    rec, _ = _run3(g, cwd)
    h = rec["dp_optical"]["hole"]
    assert [c["D_eV_A"] for c in h["channels"]] == [6.0] and h["ivs_dfpt"]["odp_reference"][0]["source"] == "manual"
    # 显式 SCATTERING 不含 ODP：谷内不进，记录写明
    g.SCATTERING, g._SCATTERING_EXPLICIT = ["ADP", "IMP", "IVS"], True
    rec, _ = _run3(g, cwd)
    assert "hole" not in rec["dp_optical"] and rec["dp_optical_ivs_dfpt"]["odp_use"]["hole"] == "not_in_scattering"
    # S8.3 汇总的来源标注
    g._SCATTERING_EXPLICIT, g.SCATTERING = False, ["ADP", "IMP", "POP"]
    g.DP_OPTICAL_H = g.DP_OPTICAL_VALLEYS_H = ""
    rec, _ = _run3(g, cwd)
    g13 = _load("g13_r3", KE / "step8.3_output" / "gen_step13_output.py")
    (cwd / "step8_amset").mkdir()
    (cwd / "step8_amset" / "dp_optical_3d.json").write_text(json.dumps(rec, ensure_ascii=False))
    lines = g13._ivs_source_lines(cwd)
    assert any(x.startswith("# [ODP 参数来源] S8 hole：DFPT q=Γ") for x in lines)
    assert any(x.startswith("# [IVS 参数来源] S8 electron：DFPT") for x in lines)


def test_r3_s84_dfpt_odp_supersedes_frozen(tmp_path, monkeypatch):
    import ke_common as kc
    import odp_frozen as of
    r_m = _rops(MOS2_LAT, MOS2_FRAC, ["Mo", "S", "S"])
    params, Q = _mos2_params(r_m)
    params["carriers"]["electron"]["channels"].append({"from": "K", "to": "K", "kind": "intra", "D_eV_A": 3.1,
                                                       "hw_meV": 50.0})
    params["raw"] = [{"carrier": "electron", "kind": "intra", "from": "K", "to": "K"}]
    mat = tmp_path / "ke"
    (mat / "step8.2_dpt").mkdir(parents=True)
    dpt = _dpt_json(Q)
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps(dpt))
    raw = dpt["results"][0]["multi_valley"]["valleys"]
    frozen = {"schema": of.SCHEMA, "dim": "2d",
              "carriers": {"electron": {"status": "ok", "valley_sig": {"use_grid": True, "sig": kc._valley_sig(raw, True)},
                                        "valleys": [{"label": "K", "groups": [{"group": 1, "modes": [9], "hw_meV": 52.0,
                                                                               "D_eV_A": 4.4}]}]}}}
    (mat / "step8.21_odp_read").mkdir()
    (mat / "step8.21_odp_read" / "odp_params.json").write_text(json.dumps(frozen))
    g14 = _load("g14_r3", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")
    monkeypatch.setattr(g14, "SCATTERING", ["ADP", "IMP", "POP"])
    monkeypatch.setattr(g14, "IVS_FROZEN", "off")
    with redirect_stdout(io.StringIO()):
        rec = g14.resolve_dp_optical(mat)                    # 还没有 DFPT：冻声子 ODP
    assert [(c["kind"], c["D_eV_A"]) for c in rec["electron"]["channels"]] == [("intra", 4.4)]
    (mat / "ivs_dfpt").mkdir()
    (mat / "ivs_dfpt" / "ivs_dfpt_params.json").write_text(json.dumps(params))
    monkeypatch.setattr(g14, "SCATTERING", ["ADP", "IMP", "POP"])
    buf = io.StringIO()
    with redirect_stdout(buf):
        rec = g14.resolve_dp_optical(mat)
    intra = [c for c in rec["electron"]["channels"] if c["kind"] == "intra"]
    assert [(c["from"], c["D_eV_A"], c["source"]) for c in intra] == [("K", 3.1, "dfpt_odp")]   # 不与冻声子 4.4 叠加
    assert len(rec["electron"]["channels"]) == 5
    fr = g14._DP_OPTICAL_FROZEN_REC
    assert fr["frozen_use"]["electron"] == "superseded_dfpt" and "electron" not in fr["digest"]
    assert kc.odp_consumer_stale(fr, frozen, full={"dp_optical": rec}) == []
    dr = g14._DP_OPTICAL_IVS_DFPT_REC
    assert dr["odp_use"]["electron"] == "used" and dr["info"]["electron"]["odp_reference"][0]["D_eV_A"] == 4.4
    assert "谷间 IVS + 谷内 ODP" in buf.getvalue() and "只作对照" in buf.getvalue()
