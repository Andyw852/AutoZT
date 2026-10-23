# -*- coding: utf-8 -*-
"""[0077] kl-dft-cpu NAC 链：S3_nac 的 vasprun.xml -> BORN -> S5/S6。

以前没有任何一步生成 BORN（S5 只复制现成的 step3_nac/BORN），且 S5 在 3090 上拿不到
别的集群的 step3_nac/，NAC=auto 永远是关的。这里覆盖：
  · gen 端 _stage_nac：BORN 优先、否则放 NAC_vasprun.xml、清掉上一轮遗留；
  · 引擎 prep 的 _born_from_vasprun：写出 phonopy 能读的 BORN，不物理 / 胞不一致则不写；
  · skill.yaml：S5_fc 软依赖 S3_nac（跨集群回传要靠它）。
"""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for _d in ("skill/kl-dft-cpu", "skill/_common/fcfit", "skill/_common/opt", "skill/_common"):
    if str(ROOT / _d) not in sys.path:
        sys.path.insert(0, str(ROOT / _d))


def _gen5():
    spec = importlib.util.spec_from_file_location(
        "kl_gen_step5_fc", str(ROOT / "skill/kl-dft-cpu/gen_step5_fc.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ gen 端
def test_stage_nac_prefers_born_then_vasprun_and_clears_stale(tmp_path):
    g = _gen5()
    disp = tmp_path / "step4_disp"
    nac = tmp_path / "step3_nac"
    disp.mkdir()
    # 无 step3_nac：清掉上一轮遗留的 BORN / NAC_vasprun.xml
    (disp / "BORN").write_text("old")
    (disp / g.NAC_VASPRUN).write_text("old")
    assert g._stage_nac(tmp_path, disp) is None
    assert not (disp / "BORN").exists() and not (disp / g.NAC_VASPRUN).exists()
    # 只有 vasprun.xml -> 放成 NAC_vasprun.xml，交给拟合作业生成 BORN
    nac.mkdir()
    (nac / "vasprun.xml").write_text("<vasprun/>")
    assert g._stage_nac(tmp_path, disp) == "vasprun"
    assert (disp / g.NAC_VASPRUN).read_text() == "<vasprun/>"
    assert not (disp / "BORN").exists()
    # 有现成 BORN -> 用 BORN，不再放 vasprun
    (nac / "BORN").write_text("born")
    assert g._stage_nac(tmp_path, disp) == "BORN"
    assert (disp / "BORN").read_text() == "born"
    assert not (disp / g.NAC_VASPRUN).exists()


def test_gen_and_engine_agree_on_the_staged_name():
    g = _gen5()
    src = (ROOT / "skill/_common/fcfit/fc_fit_driver.py").read_text(encoding="utf-8")
    assert 'NAC_VASPRUN = "%s"' % g.NAC_VASPRUN in src


def test_s5_soft_depends_on_s3_nac_for_cross_cluster_relay():
    import yaml
    y = yaml.safe_load((ROOT / "skill/kl-dft-cpu/skill.yaml").read_text(encoding="utf-8"))
    s5 = next(s for s in y["steps"] if s["name"] == "step5_fc")
    assert "step3_nac?" in s5["needs"] and "step4_disp" in s5["needs"]


# ------------------------------------------------------------------ 引擎 prep
def _rocksalt(a=5.64):
    from phonopy.structure.atoms import PhonopyAtoms
    lat = [[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]]
    return PhonopyAtoms(symbols=["Na", "Cl"], cell=lat,
                        scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]])


def _setup(tmp_path, monkeypatch, eps, z, cell=None):
    pytest.importorskip("phonopy")
    pytest.importorskip("ase")
    import numpy as np
    import phonopy.interface.vasp as pv
    from phonopy.interface.vasp import write_vasp
    import fc_fit_driver as d
    src, out = tmp_path / "src", tmp_path / "out"
    src.mkdir()
    out.mkdir()
    write_vasp(str(out / "POSCAR"), _rocksalt())
    (src / d.NAC_VASPRUN).write_text("<vasprun/>")

    class _Fake:                       # 只替换 vasprun 解析，张量约化走真 phonopy
        def __init__(self, fh):
            pass

        def parse(self):
            self.cell = cell if cell is not None else _rocksalt()
            self.born = np.asarray(z, float)
            self.epsilon = np.asarray(eps, float)
    monkeypatch.setattr(pv, "VasprunxmlExpat", _Fake)
    info = {"primitive_matrix": np.eye(3).tolist(), "supercell_matrix": (4 * np.eye(3)).tolist()}
    return d, src, out, info


def test_born_written_and_readable_by_phonopy(tmp_path, monkeypatch):
    import numpy as np
    eps = np.diag([2.4, 2.4, 2.4])
    z = [np.eye(3) * 1.1, -np.eye(3) * 1.1]
    d, src, out, info = _setup(tmp_path, monkeypatch, eps, z)
    st = d._born_from_vasprun(src, out, info)
    assert st["written"] and (out / "BORN").is_file()
    from phonopy import Phonopy
    from phonopy.file_IO import parse_BORN
    ph = Phonopy(_rocksalt(), supercell_matrix=np.eye(3) * 2)
    nac = parse_BORN(ph.primitive, filename=str(out / "BORN"))
    assert np.allclose(nac["dielectric"], eps)
    assert np.allclose(nac["born"][0], np.eye(3) * 1.1, atol=1e-6)
    assert np.allclose(nac["born"][1], -np.eye(3) * 1.1, atol=1e-6)


def test_dataset_born_wins_and_stale_out_born_is_rebuilt(tmp_path, monkeypatch):
    import numpy as np
    d, src, out, info = _setup(tmp_path, monkeypatch, np.diag([2.4] * 3),
                               [np.eye(3), -np.eye(3)])
    (out / "BORN").write_text("stale")
    assert d._born_from_vasprun(src, out, info)["written"]
    assert "stale" not in (out / "BORN").read_text()
    (src / "BORN").write_text("dataset")          # 数据集自带 BORN -> 不碰
    assert d._born_from_vasprun(src, out, info) is None


def test_metal_like_epsilon_gives_no_born(tmp_path, monkeypatch):
    import numpy as np
    d, src, out, info = _setup(tmp_path, monkeypatch, np.diag([-5.0, 2.0, 2.0]),
                               [np.eye(3), -np.eye(3)])
    st = d._born_from_vasprun(src, out, info)
    assert not st["written"] and not (out / "BORN").exists()
    assert "not physical" in st["reason"]


def test_different_cell_gives_no_born(tmp_path, monkeypatch):
    import numpy as np
    d, src, out, info = _setup(tmp_path, monkeypatch, np.diag([2.4] * 3),
                               [np.eye(3), -np.eye(3)], cell=_rocksalt(a=6.2))
    st = d._born_from_vasprun(src, out, info)
    assert not st["written"] and not (out / "BORN").exists()
    assert "lattice differs" in st["reason"]
