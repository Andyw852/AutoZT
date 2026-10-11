"""[0035] 冻声子 ODP 的生成/读取脚本真跑一遍（需要 pymatgen；没有就跳过）。

把 S7.2 的 gen、S8.21 的 gen 连同 gen_need 复制进一个临时"材料目录"（与 autozt 推到超算的布局相同），用紧束缚玩具
代替 VASP：S5 的声子由 test_0035_odp_frozen_tau 的弹簧模型给出（monkeypatch read_dfpt），S7.2 各构型的"本征值"
由同一个紧束缚模型按该构型 POSCAR 的原子位置算出（monkeypatch 读取端的 Config），真空能级按 S 层间距人为平移。
核对：S7.2 的目录/位移/输入照抄/选模变化时的改名；S8.21 的 D 与 Hellmann–Feynman 直接算的 √M·⟨∂H/∂Q⟩ − dV_vac/du
一致；S8.4 把它变成谷内通道。
"""
import importlib.util
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pymatgen")

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
COMMON = ROOT / "skill" / "_common" / "opt"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_0035_odp_frozen_tau as T  # noqa: E402

POSCAR = """toy
1.0
3.2 0.0 0.0
0.0 3.2 0.0
0.0 0.0 20.0
Mo S
1 2
Direct
0.0 0.0 0.5
0.5 0.5 0.575
0.5 0.5 0.425
"""
INCAR_S7 = "SYSTEM = toy deform undeformed\nISTART = 0\nENCUT = 500\nISYM = 0\nLVHAR = .TRUE.\nNCORE = 6\n"
NK = 12
VAC_COEF = 0.3                    # 人为的真空能级：V = 0.3·(z_S1 − z_S2 − 3) eV


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def mat(tmp_path, monkeypatch):
    for f in (KE / "step7c_odp_frozen" / "gen_step9c_odp.py",
              KE / "step7c_odp_frozen" / "odp_read" / "gen_step12c_odp_read.py",
              KE / "ke_common.py", KE / "odp_frozen.py", COMMON / "stepconf.py", COMMON / "dim_common.py",
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
    (und / "POTCAR").write_text("dummy POTCAR Mo S\n")
    (tmp_path / "step5_dielect").mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.syspath_prepend(str(tmp_path))
    for m in ("ke_common", "odp_frozen", "stepconf", "dim_common"):
        sys.modules.pop(m, None)
    return tmp_path


def _gen72(mat, monkeypatch):
    from pymatgen.core import Structure
    g = _load("gen72_0035", mat / "gen_step9c_odp.py")
    f, E = T._phonons()
    born = np.array([-2 * np.eye(3), np.eye(3), np.eye(3)])
    monkeypatch.setattr(g, "read_dfpt", lambda cwd: (f, E, born, Structure.from_str(POSCAR, fmt="poscar")))
    return g, f, E


def _cart(st):
    """pymatgen 结构 -> 紧束缚模型的笛卡尔位置（模型只看相对位置）。"""
    return np.asarray(st.cart_coords, float)


def _vac(pos):
    return VAC_COEF * (pos[1, 2] - pos[2, 2] - 3.0)


def test_s72_gen_and_s821_read(mat, monkeypatch, capsys):
    g, f, E = _gen72(mat, monkeypatch)
    g.main()
    out = mat / "step7c_odp_frozen"
    man = json.loads((out / "odp_manifest.json").read_text())
    sel = [s["index"] for s in man["selected"]]
    assert sel == [4, 5, 6]                                     # 两个面内同极模（简并）+ A1′ 型；极性模不算
    names = sorted(p.name for p in out.glob("odp-*"))
    assert names == ["odp-m04m", "odp-m04p", "odp-m05m", "odp-m05p", "odp-m06m", "odp-m06p", "odp-ref"]
    und = mat / "step7_deform" / "undeformed"
    for nm in names:
        d = out / nm
        assert (d / "KPOINTS").read_text() == (und / "KPOINTS").read_text()
        assert (d / "POTCAR").read_text() == (und / "POTCAR").read_text()
        inc = (d / "INCAR").read_text().splitlines()
        assert inc[0] == "SYSTEM = %s odp %s" % (mat.name, nm) and inc[1:] == INCAR_S7.splitlines()[1:]
        assert "{{JOBNAME}}" not in (d / "submit.sh").read_text()
    from pymatgen.core import Structure
    base = _cart(Structure.from_file(str(out / "odp-ref" / "POSCAR")))
    masses = np.array(man["masses_amu"])
    for s in man["selected"]:
        e = np.array(s["eigvec"])
        dlt, u = T._load("of_pipe", KE / "odp_frozen.py").displacement(e, masses, man["disp_max_A"])
        assert u == pytest.approx(s["u_A"])
        for sign, tag in ((1, "p"), (-1, "m")):
            p = _cart(Structure.from_file(str(out / ("odp-m%02d%s" % (s["index"], tag)) / "POSCAR")))
            assert np.allclose(p - base, sign * dlt, atol=1e-8)

    # ---- S8.21：用紧束缚模型充当各构型的 VASP 结果 ----
    r = _load("gen821_0035", mat / "gen_step12c_odp_read.py")
    kgrid = np.array([[i / NK, j / NK, 0.0] for i in range(NK) for j in range(NK)])

    class FakeConfig(r.Config):
        def __init__(self, d, vac_axis=2):
            st = Structure.from_file(str(d / "POSCAR"))
            pos = _cart(st)
            self.name, self.structure, self.kpts = d.name, st, kgrid
            self.eig = {"up": np.array([np.linalg.eigvalsh(T._H(k, pos)) for k in kgrid])}
            self.vac, self.vac_std = _vac(pos), 0.0
    monkeypatch.setattr(r, "Config", FakeConfig)
    ka, kb = [3 / NK, 1 / NK, 0.0], [1 / NK, 3 / NK, 0.0]          # x↔y 镜面相关的两个成员谷
    valley = {"label": "X", "g": 2, "spin": "up", "band_index": 2, "dE_eV": 0.0,
              "k_frac": [ka, kb], "k_grid": [ka, kb], "k_S7": [ka, kb]}
    (mat / "step8.2_dpt").mkdir()
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": [
        {"carrier": "electron", "multi_valley": {"status": "ok", "valleys": [valley]}},
        {"carrier": "hole", "multi_valley": {"status": "rejected"}}]}))
    r.main()
    p = json.loads((mat / "step8.21_odp_read" / "odp_params.json").read_text())
    assert p["carriers"]["hole"]["status"] == "no_valleys"
    v = p["carriers"]["electron"]["valleys"][0]
    got = {tuple(gr["modes"]): gr for gr in v["groups"]}
    assert set(got) == {(4, 5), (6,)}
    # 参考：每个模 √M·⟨∂H/∂Q⟩ 减真空项，成员谷平均、简并伙伴 D² 相加
    msum = masses.sum()
    for modes, gr in got.items():
        per = []
        for k in (ka, kb):
            w0, V0 = np.linalg.eigh(T._H(k, base))
            tot = 0.0
            for i in modes:
                s = next(x for x in man["selected"] if x["index"] == i)
                pat = np.array(s["eigvec"]) / np.sqrt(masses)[:, None]
                dq = 1e-6
                dH = (T._H(k, base + dq * pat) - T._H(k, base - dq * pat)) / (2 * dq)
                lam = math.sqrt(msum) * float(np.real(V0[:, 2].conj() @ dH @ V0[:, 2]))
                dvdu = math.sqrt(msum) * (_vac(base + dq * pat) - _vac(base - dq * pat)) / (2 * dq)
                tot += (lam - dvdu) ** 2
            per.append(tot)
        assert gr["D_eV_A"] == pytest.approx(math.sqrt(np.mean(per)), rel=2e-3, abs=1e-4)
        assert gr["spread"] < 1e-3                                   # 对称等价的两个成员谷 ΣD² 相同
    a1 = got[(6,)]
    assert abs(a1["D_eV_A"] - a1["D_raw_eV_A"]) > 0.05                # 真空对齐确实起作用（A1′ 改变层厚）

    # ---- S8.4：冻声子参数 -> 谷内通道 ----
    sys.path.insert(0, str(KE / "step8.4_amset2d"))
    g14 = _load("gen14_pipe_0035", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")
    monkeypatch.setattr(g14, "SCATTERING", ["ADP", "IMP", "POP"])
    rec = g14.resolve_dp_optical(mat)
    chs = rec["electron"]["channels"]
    assert {c["kind"] for c in chs} == {"intra"} and all(c["from"] == "X" for c in chs)
    assert "ODP" in g14.SCATTERING and g14._DP_OPTICAL_FROZEN_REC["frozen_use"]["electron"] == "used"


def test_s72_mode_change_moves_stale_dirs(mat, monkeypatch):
    g, _f, _E = _gen72(mat, monkeypatch)
    g.main()
    out = mat / "step7c_odp_frozen"
    (out / "odp-m06p" / "OUTCAR").write_text("General timing and accounting informations\n")
    (mat / "step.conf").write_text("[params]\nODP_MODES = 6\n", encoding="utf-8")
    g, _f, _E = _gen72(mat, monkeypatch)
    g.main()
    assert sorted(p.name for p in out.glob("odp-*")) == ["odp-m06m", "odp-m06p", "odp-ref"]
    assert len(list(out.glob("_stale_odp-m04p_*"))) == 1
    assert (out / "odp-m06p" / "OUTCAR").is_file()                    # 输入没变：已算完的不归档
    man = json.loads((out / "odp_manifest.json").read_text())
    assert [s["index"] for s in man["selected"]] == [6] and man["overrides"] == {"ODP_MODES": "6"}


def test_s72_refuses_wrong_structure_and_missing_s7(mat, monkeypatch):
    (mat / "step7_deform" / "undeformed" / "POSCAR").write_text(POSCAR.replace("0.575", "0.580"))
    g, _f, _E = _gen72(mat, monkeypatch)
    with pytest.raises(SystemExit):
        g.main()
    (mat / "step7_deform" / "undeformed" / "KPOINTS").unlink()
    with pytest.raises(SystemExit):
        g.main()


def test_s821_refuses_changed_s7_inputs(mat, monkeypatch):
    g, _f, _E = _gen72(mat, monkeypatch)
    g.main()
    (mat / "step7_deform" / "undeformed" / "KPOINTS").write_text("auto\n0\nGamma\n 15 15 1\n")
    (mat / "step8.2_dpt").mkdir()
    (mat / "step8.2_dpt" / "dpt_result.json").write_text(json.dumps({"results": []}))
    r = _load("gen821b_0035", mat / "gen_step12c_odp_read.py")
    with pytest.raises(SystemExit) as e:
        r.main()
    assert "retry S7.2" in str(e.value)
