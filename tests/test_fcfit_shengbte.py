"""fit-fc-thermal 的 ShengBTE 导出：FC_2ND 必须是**超胞**（ShengBTE read2fc 要求
ntot = natoms*prod(scell)），3RD/POSCAR/CONTROL 是原胞；超胞顺序按 ShengBTE 的
split_index（x 最快、原胞原子最慢）。纯本地，不需要 ShengBTE 本体。"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skill" / "_common"))
ase = pytest.importorskip("ase")
from ase import Atoms  # noqa: E402
from ase.io import write  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "fcfd", ROOT / "skill" / "_common" / "fcfit" / "fc_fit_driver.py")
D = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(D)


def _split_index(index, nx, ny, nz):
    """ShengBTE Src/input.f90 split_index（1 起）。"""
    t, ix = divmod(index - 1, nx)
    t2, iy = divmod(t, ny)
    ia, iz = divmod(t2, nz)
    return ix, iy, iz, ia


def _setup(tmp_path, dims=(3, 3, 1), scramble=False, fc3_index=2, dim="2d"):
    uc = Atoms("MoS", cell=[[3.16, 0, 0], [-1.58, 2.7366, 0], [0, 0, 20.0]],
               scaled_positions=[[1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.56]], pbc=True)
    nx, ny, nz = dims
    rows = []                                   # phonopy/ShengBTE 顺序
    for ia in range(len(uc)):
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    rows.append((ia, (uc.get_scaled_positions()[ia] + [ix, iy, iz]) / dims))
    order = list(range(len(rows)))
    if scramble:                                # 像 ASE repeat 一样交错
        order = sorted(order, key=lambda k: (k % (nx * ny * nz), k))
    sc = Atoms([uc[rows[k][0]].symbol for k in order],
               cell=np.diag(dims) @ uc.cell, scaled_positions=[rows[k][1] for k in order],
               pbc=True)
    n = len(sc)
    # 可辨认的力常数：fc[a,b] = (原子 a 的 ShengBTE 索引, b 的索引) 编码
    true = np.zeros((n, n, 3, 3))
    for a in range(n):
        for b in range(n):
            true[a, b] = (a + 1) * 1000 + (b + 1)
    fc_in_sc_order = true[np.ix_(order, order)]
    sb = tmp_path / "shengbte"
    sb.mkdir()
    write(str(tmp_path / "POSCAR"), uc, format="vasp", direct=True)
    write(str(sb / "POSCAR"), uc, format="vasp", direct=True)
    write(str(tmp_path / "SPOSCAR"), sc, format="vasp", direct=True)
    D._write_fc2_text(fc_in_sc_order, sb / "FORCE_CONSTANTS_2ND")
    blk = ["1", "", "1", "0 0 0", "0 0 0", "1 %d 1" % fc3_index] + ["1 1 1 0.0"] * 27
    (sb / "FORCE_CONSTANTS_3RD").write_text("\n".join(blk) + "\n")
    cfg = {"export_shengbte": True, "dim": dim, "shengbte_t": "100 800 100"}
    return cfg, true


def test_supercell_fc2_is_the_correct_shengbte_layout(tmp_path):
    cfg, true = _setup(tmp_path)
    assert D._shengbte_finalize(cfg, tmp_path)
    man = json.loads((tmp_path / "shengbte" / "shengbte_manifest.json").read_text())
    assert man["usable"] and man["natoms"] == 2 and man["scell"] == [3, 3, 1]
    assert man["checks"]["fc2_header"] == 18 == 2 * 9      # 超胞，不是 "2 2"
    ctl = (tmp_path / "shengbte" / "CONTROL").read_text()
    assert "natoms=2," in ctl and "scell(:)=3 3 1" in ctl
    assert "T_min=100, T_max=800, T_step=100" in ctl
    assert man["ngrid"][2] == 1                              # 2D：真空方向 1 个 q 点
    assert "lfactor=0.1," in ctl and 'elements="Mo" "S"' in ctl


def test_scrambled_supercell_is_reordered_for_shengbte(tmp_path):
    cfg, true = _setup(tmp_path, scramble=True)
    assert D._shengbte_finalize(cfg, tmp_path)
    man = json.loads((tmp_path / "shengbte" / "shengbte_manifest.json").read_text())
    assert man["checks"].get("fc2_reordered")
    fc = D._read_fc2_text(tmp_path / "shengbte" / "FORCE_CONSTANTS_2ND")
    # 按 ShengBTE 的 split_index 读回来，每个块都对上真实的 (a,b)
    for a in (1, 5, 17):
        for b in (2, 9, 18):
            assert fc[a - 1, b - 1][0, 0] == a * 1000 + b
            assert _split_index(a, 3, 3, 1)[3] == (a - 1) // 9


def test_inconsistent_set_is_flagged_not_run(tmp_path):
    cfg, _ = _setup(tmp_path, fc3_index=5)                  # 3RD 里出现原胞外的原子号
    assert not D._shengbte_finalize(cfg, tmp_path)
    man = json.loads((tmp_path / "shengbte" / "shengbte_manifest.json").read_text())
    assert not man["usable"] and "natoms" in man["error"]
    assert not (tmp_path / "shengbte" / "CONTROL").exists()
    # 2ND 被折叠成原胞（"2 2"）——ShengBTE 会报 wrong number of force constants
    sb = tmp_path / "shengbte"
    D._write_fc2_text(np.zeros((2, 2, 3, 3)), sb / "FORCE_CONSTANTS_2ND")
    assert not D._shengbte_finalize(cfg, tmp_path)
    assert "expected natoms*prod(scell)" in json.loads(
        (sb / "shengbte_manifest.json").read_text())["error"]


def test_pheasy_bin_only_checked_for_pheasy_engine():
    src = (ROOT / "skill" / "_common" / "fcfit" / "gen_step1_fit.py").read_text(encoding="utf-8")
    i = src.index('p_bin = str(conf["PHEASY_BIN"]')
    block = src[i:i + 700]
    assert 'if engine == "pheasy":' in block and "ignored" in block
