"""[0050] 谷间 IVS（K↔K′ 零阶）冻声子的生成/读取脚本真跑一遍（需要 pymatgen；没有就跳过）。

把 S7.3、S7.4、S8.22 的 gen 连同 gen_need 复制进一个临时"材料目录"（与 autozt 推到超算的布局相同），用玩具模型代替 VASP：
  · 六角单原子层（γ = 120°，一个原子、两个轨道：s 轨道最近邻跳跃 t、在位能 α·Δz；另一条平带在 +8 eV 充当 NBANDS 截断处）；
  · S7.3 的超胞 Γ̃ 声子由最近邻弹簧（中心力 + 面外弯曲）给出（monkeypatch read_dfpt_sc）；
  · S7.4 各构型的"本征值"由同一个紧束缚模型按该构型 POSCAR 的原子位置、在它 KPOINTS 的网格上算出（monkeypatch Config）。
核对：S7.3 超胞 det 3、网格换算；S7.4 选出 6 个区边界模（面外 2 + 面内 4）、单点目录与位移、NBANDS；
S8.22 的面外组 D = α、面内组 D ≈ 0、K–M 近邻告警；S8.4 把它变成 K>K' 通道并打开 ADP 谷内掩码；
参数变了才让 S8.4 重排；SOC / 非六角 / S7 输入变了都拒绝。
"""
import importlib.util
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pymatgen")

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
COMMON = ROOT / "skill" / "_common" / "opt"

A = 3.2
C = 20.0
POSCAR = """toy hex
1.0
%.10f 0.0 0.0
%.10f %.10f 0.0
0.0 0.0 %.1f
Mo
1
Direct
0.0 0.0 0.5
""" % (A, -A / 2, A * math.sqrt(3) / 2, C)
INCAR_S7 = "SYSTEM = toy deform undeformed\nISTART = 0\nENCUT = 500\nISYM = 0\nLVHAR = .TRUE.\nNCORE = 6\n"
NK = 12
T_HOP = 0.4                       # eV；带底在 K/K′（−3t）
ALPHA = 1.7                       # eV/Å：在位能 = α·(z − z0)，只有面外位移耦合
Z0 = 0.5 * C
FLAT = 8.0                        # 第二条（平）带
KC, KZ = 4.0, 1.0                 # eV/Å²：中心力弹簧、面外弯曲弹簧
EV_A2_AMU = 1.602176634e-19 / (1e-20 * 1.66053906660e-27)       # (eV/Å²/amu) -> s⁻²


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def mat(tmp_path, monkeypatch):
    for f in (KE / "step7d_ivs_phonon" / "gen_step9d_ivs_ph.py",
              KE / "step7e_ivs_frozen" / "gen_step9e_ivs.py",
              KE / "step7e_ivs_frozen" / "ivs_read" / "gen_step12d_ivs_read.py",
              KE / "ke_common.py", KE / "ivs_frozen.py", KE / "odp_frozen.py",
              COMMON / "stepconf.py", COMMON / "dim_common.py",
              ROOT / "setting" / "example" / "templates" / "submit_std_2d.tpl"):
        shutil.copy(f, tmp_path)
    (tmp_path / "step.conf").write_text("[params]\n", encoding="utf-8")
    (tmp_path / "step1_opt").mkdir()
    (tmp_path / "step1_opt" / "CONTCAR").write_text(POSCAR)
    (tmp_path / "step1_opt" / "workflow_method.txt").write_text("DIM=2D\n")
    und = tmp_path / "step7_deform" / "undeformed"
    und.mkdir(parents=True)
    (und / "POSCAR").write_text(POSCAR)
    (und / "INCAR").write_text(INCAR_S7)
    (und / "KPOINTS").write_text("auto\n0\nGamma\n %d %d 1\n" % (NK, NK))
    (und / "POTCAR").write_text("  PAW_PBE Mo_pv\n   POMASS =   95.950; ZVAL   =   14.000    mass and valenz\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.syspath_prepend(str(tmp_path))
    for m in ("ke_common", "ivs_frozen", "odp_frozen", "stepconf", "dim_common"):
        sys.modules.pop(m, None)
    return tmp_path


# ---------------------------------------------------------------------------
# 玩具模型
# ---------------------------------------------------------------------------
def _pairs(lat, cart, cutoff=1.15 * A):
    """(i, j, R, 键矢量) —— 最近邻（含面内周期像）。"""
    out = []
    for i in range(len(cart)):
        for j in range(len(cart)):
            for a in range(-2, 3):
                for b in range(-2, 3):
                    d = cart[j] + a * lat[0] + b * lat[1] - cart[i]
                    if 1e-6 < np.linalg.norm(d) < cutoff:
                        out.append((i, j, (a, b), d))
    return out


def _eigs(st, kpts):
    """紧束缚：每个原子 s 轨道（在位 α·Δz、最近邻 t）+ 平带 FLAT。返回 (nk, 2·nat) 升序本征值。"""
    lat = np.asarray(st.lattice.matrix, float)
    cart = np.asarray(st.cart_coords, float)
    n = len(cart)
    pr = _pairs(lat, cart)
    out = []
    for k in kpts:
        H = np.zeros((n, n), complex)
        for i, j, R, _d in pr:
            H[i, j] += T_HOP * np.exp(2j * math.pi * (k[0] * R[0] + k[1] * R[1]))
        H += np.diag(ALPHA * (cart[:, 2] - Z0))
        e = np.concatenate([np.linalg.eigvalsh(H), np.full(n, FLAT)])
        out.append(np.sort(e))
    return np.array(out)


def _phonons(st):
    """超胞 Γ̃ 点：最近邻中心力 + 面外弯曲弹簧 -> (频率 THz, 质量加权归一本征矢 (3N, N, 3))。"""
    lat = np.asarray(st.lattice.matrix, float)
    cart = np.asarray(st.cart_coords, float)
    n = len(cart)
    m = np.array([float(s.specie.atomic_mass) for s in st])
    phi = np.zeros((n, 3, n, 3))
    for i, j, _R, d in _pairs(lat, cart):
        e = d / np.linalg.norm(d)
        blk = -KC * np.outer(e, e)
        blk[2, 2] += -KZ
        phi[i, :, j, :] += blk
        phi[i, :, i, :] -= blk
    D = phi.reshape(3 * n, 3 * n) / np.sqrt(np.repeat(m, 3)[:, None] * np.repeat(m, 3)[None, :])
    w, v = np.linalg.eigh(0.5 * (D + D.T))
    f = np.sign(w) * np.sqrt(np.abs(w) * EV_A2_AMU) / (2 * math.pi * 1e12)
    return f, v.T.reshape(3 * n, n, 3)


def _mesh(kp):
    ln = Path(kp).read_text().splitlines()
    n = [int(x) for x in ln[3].split()[:3]]
    return np.array([[i / n[0], j / n[1], 0.0] for i in range(n[0]) for j in range(n[1])])


def _gen73(mat):
    g = _load("gen73_0050", mat / "gen_step9d_ivs_ph.py")
    g.main()
    return g


def _gen74(mat, monkeypatch):
    from pymatgen.core import Structure
    g = _load("gen74_0050", mat / "gen_step9e_ivs.py")
    ph = mat / "step7d_ivs_phonon"
    st = Structure.from_file(str(ph / "POSCAR"))
    f, E = _phonons(st)
    (ph / "vasprun.xml").write_text("<fake dfpt/>")
    monkeypatch.setattr(g, "read_dfpt_sc", lambda d: (f, E, st, ph / "vasprun.xml"))
    g.main()
    return g, f, E


def _fake_reader(mat, monkeypatch, name="gen822_0050"):
    from pymatgen.core import Structure
    r = _load(name, mat / "gen_step12d_ivs_read.py")

    class FakeConfig(r.Config):
        def __init__(self, d):
            st = Structure.from_file(str(d / "POSCAR"))
            self.name, self.structure, self.kpts = d.name, st, _mesh(d / "KPOINTS")
            self.eig = {"up": _eigs(st, self.kpts)}
    monkeypatch.setattr(r, "Config", FakeConfig)
    return r


K_STAR = [[1 / 3, 1 / 3, 0.0], [2 / 3, 2 / 3, 0.0]]
M_STAR = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.5, 0.5, 0.0]]


def _dpt(mat, extra=True):
    kv = {"label": "K", "g": 2, "spin": "up", "band_index": 0, "dE_eV": 0.0,
          "k_frac": K_STAR, "k_grid": K_STAR, "k_S7": K_STAR}
    mv = {"label": "M", "g": 3, "spin": "up", "band_index": 0, "dE_eV": 0.15,
          "k_frac": M_STAR, "k_grid": M_STAR, "k_S7": M_STAR}
    (mat / "step8.2_dpt").mkdir(exist_ok=True)
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": "electron", "multi_valley": {"status": "ok", "valleys": [kv, mv] if extra else [kv]}},
        {"carrier": "hole", "multi_valley": {"status": "rejected"}}]}))


# ---------------------------------------------------------------------------
def test_full_pipeline(mat, monkeypatch, capsys):
    _gen73(mat)
    ph = mat / "step7d_ivs_phonon"
    man3 = json.loads((ph / "ivs_supercell.json").read_text())
    assert man3["det"] == 3 and man3["kmesh_sc"] == [7, 7, 1] and man3["kmesh_s7"] == [NK, NK, 1]
    inc = (ph / "INCAR").read_text()
    assert re.search(r"IBRION\s*=\s*8", inc) and "LVHAR" not in inc and re.search(r"NCORE\s*=\s*1\b", inc)
    assert len(man3["star"]) == 2 and (ph / "POSCAR.prim").read_text().split()[-3:] == POSCAR.split()[-3:]

    g, f, E = _gen74(mat, monkeypatch)
    out = mat / "step7e_ivs_frozen"
    man = json.loads((out / "ivs_manifest.json").read_text())
    grp = sorted(len(gr["modes"]) for gr in man["groups"] if gr["use"])
    assert grp == [2, 4]                                         # 面外 K/K′ 一对、面内两支各一对
    names = sorted(p.name for p in out.glob("ivs-*"))
    assert len(names) == 2 + 2 * 6 and "ivs-ref" in names and "ivs-prim" in names
    assert man["kmesh"]["prim"][0] % 3 == 0 and man["nbands"] == {"supercell": 38, "prim": 13, "nelect_prim": 14.0}
    for nm in names:
        assert re.search(r"NBANDS\s*=\s*%d\b" % (13 if nm == "ivs-prim" else 38), (out / nm / "INCAR").read_text())
        assert "{{JOBNAME}}" not in (out / nm / "submit.sh").read_text()
    from pymatgen.core import Structure
    base = np.asarray(Structure.from_file(str(out / "ivs-ref" / "POSCAR")).cart_coords)
    iv = _load("iv_pipe_0050", KE / "ivs_frozen.py")
    for s in man["selected"]:
        dlt, u = iv.displacement(np.array(s["eigvec"]), man["masses_amu"], man["disp_max_A"])
        assert u == pytest.approx(s["u_A"])
        for sign, tag in ((1, "p"), (-1, "m")):
            p = np.asarray(Structure.from_file(str(out / ("ivs-m%02d%s" % (s["index"], tag)) / "POSCAR")).cart_coords)
            assert np.allclose(p - base, sign * dlt, atol=1e-8)

    # ---- S8.22 ----
    _dpt(mat)
    r = _fake_reader(mat, monkeypatch)
    r.main()
    p = json.loads((mat / "step8.22_ivs_read" / "ivs_params.json").read_text())
    assert p["carriers"]["hole"]["status"] == "no_valleys"
    v = p["carriers"]["electron"]["valleys"]
    assert len(v) == 1 and v[0]["label"] == "K" and v[0]["g"] == 2
    got = {len(gr["modes"]): gr for gr in v[0]["groups"]}
    assert got[2]["D_eV_A"] == pytest.approx(ALPHA, rel=2e-3)      # 面外：D = α
    assert got[4]["D_eV_A"] == pytest.approx(0.0, abs=1e-3)         # 面内：不耦合
    assert any("M 星离 K 星只有 150 meV" in w for w in p["warnings"])
    log = capsys.readouterr().out
    assert "[DONE]" in log

    # ---- S8.4：K>K' 通道 ----
    sys.path.insert(0, str(KE / "step8.4_amset2d"))
    g14 = _load("gen14_pipe_0050", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")
    monkeypatch.setattr(g14, "SCATTERING", ["ADP", "IMP", "POP"])
    monkeypatch.setattr(g14, "DP_OPTICAL_E", "K>K':9.9@30")              # 手填的同星通道只作对照
    rec = g14.resolve_dp_optical(mat)
    chs = rec["electron"]["channels"]
    assert len(chs) == 1 and chs[0]["kind"] == "inter_star" and chs[0]["from"] == "K"
    assert chs[0]["D_eV_A"] == pytest.approx(ALPHA, rel=2e-3)
    assert chs[0]["hw_meV"] == pytest.approx(got[2]["hw_meV"], abs=1e-3)
    assert "IVS" in g14.SCATTERING and rec["adp_intravalley"]["electron"] is True
    ivr = g14._DP_OPTICAL_IVS_REC
    assert ivr["frozen_use"]["electron"] == "used" and ivr["info"]["electron"]["manual_reference"][0]["D_eV_A"] == 9.9
    assert "hole" not in rec

    # ---- 参数变了才让 S8.4 重排 ----
    d84 = mat / "step8.4_amset2d"
    d84.mkdir()
    (d84 / "2d_correction.json").write_text(json.dumps({"dp_optical": rec, "dp_optical_ivs": ivr}))
    (d84 / "transport.json").write_text("{}")
    kc = sys.modules["ke_common"]
    assert kc.invalidate_ivs_consumer(mat, p) == []
    p2 = json.loads(json.dumps(p))
    p2["carriers"]["electron"]["valleys"][0]["groups"][0]["D_eV_A"] += 0.3
    done = kc.invalidate_ivs_consumer(mat, p2)
    assert ("step8.4_amset2d", "transport.json") in done and not (d84 / "transport.json").exists()


def test_s822_refuses_changed_s7_inputs(mat, monkeypatch):
    _gen73(mat)
    _gen74(mat, monkeypatch)
    (mat / "step7_deform" / "undeformed" / "INCAR").write_text(INCAR_S7.replace("500", "520"))
    _dpt(mat)
    r = _fake_reader(mat, monkeypatch, "gen822b_0050")
    with pytest.raises(SystemExit) as e:
        r.main()
    assert "retry S7.3_ivsph" in str(e.value)


def test_s74_refuses_changed_s7_since_s73(mat, monkeypatch):
    _gen73(mat)
    (mat / "step7_deform" / "undeformed" / "POTCAR").write_text("  PAW_PBE Mo\n ZVAL = 6.000\n")
    with pytest.raises(SystemExit) as e:
        _gen74(mat, monkeypatch)
    assert "retry S7.3_ivsph" in str(e.value)


def test_soc_refused_everywhere(mat, monkeypatch):
    (mat / "step3_uniform").mkdir()
    (mat / "step3_uniform" / "INCAR").write_text("LSORBIT = .TRUE.\n")
    with pytest.raises(SystemExit) as e:
        _gen73(mat)
    assert "SOC" in str(e.value) and "ivs_frozen: false" in str(e.value)
    # S8.4：有参数而能带是 SOC -> 拒绝
    (mat / "step8.22_ivs_read").mkdir()
    (mat / "step8.22_ivs_read" / "ivs_params.json").write_text(json.dumps({"schema": "ivs_params/1", "carriers": {}}))
    sys.path.insert(0, str(KE / "step8.4_amset2d"))
    g14 = _load("gen14_soc_0050", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")
    with pytest.raises(SystemExit) as e:
        g14._load_ivs(mat)
    assert "SOC" in str(e.value) and "IVS_FROZEN = off" in str(e.value)
    monkeypatch.setattr(g14, "IVS_FROZEN", "off")
    assert g14._load_ivs(mat) == (None, "off")


def test_non_hexagonal_refused(mat):
    sq = POSCAR.replace("%.10f %.10f 0.0" % (-A / 2, A * math.sqrt(3) / 2), "0.0 %.10f 0.0" % A)
    (mat / "step1_opt" / "CONTCAR").write_text(sq)
    (mat / "step7_deform" / "undeformed" / "POSCAR").write_text(sq)
    with pytest.raises(SystemExit) as e:
        _gen73(mat)
    assert "不是六角" in str(e.value)


def test_no_k_star_and_skill_wiring(mat, monkeypatch):
    _gen73(mat)
    _gen74(mat, monkeypatch)
    (mat / "step8.2_dpt").mkdir()
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": "electron", "multi_valley": {"status": "ok", "valleys": [
            {"label": "M", "g": 3, "spin": "up", "band_index": 0, "dE_eV": 0.0,
             "k_frac": M_STAR, "k_grid": M_STAR, "k_S7": M_STAR}]}}]}))
    r = _fake_reader(mat, monkeypatch, "gen822c_0050")
    r.main()
    p = json.loads((mat / "step8.22_ivs_read" / "ivs_params.json").read_text())
    assert p["carriers"]["electron"]["status"] == "no_k_star"
    assert iv_digest(p) is None


def iv_digest(p):
    iv = _load("iv_dig_0050", KE / "ivs_frozen.py")
    return iv.params_digest(p, "electron")


def test_skill_yaml_group_default_off():
    import yaml
    sk = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    grp = sk["optional_steps"]["ivs_frozen"]
    assert grp["default"] is False
    by = {s["name"]: s for s in grp["steps"]}
    s73, s74, s822 = by["step7d_ivs_phonon"], by["step7e_ivs_frozen"], by["step8.22_ivs_read"]
    assert s74["fanout"] == "ivs-*" and s74["check"] == "outcar"
    for s in (s73, s74, s822):
        assert {"ivs_frozen.py", "odp_frozen.py", "ke_common.py"} <= set(s["gen_need"])
    assert s822["done_marker"] == "ivs_params.json" and s822["needs"] == ["step7e_ivs_frozen", "step8.2_dpt"]
    s84 = next(s for s in sk["optional_steps"]["amset2d"]["steps"] if s["name"] == "step8.4_amset2d")
    assert "step8.22_ivs_read?" in s84["needs"] and "ivs_frozen.py" in s84["gen_need"]
