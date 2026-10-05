"""pheasy 超胞序 / 输入换序 / 导出自检的回归防护。

锁住 Mn2In2Se5 (R-3m, 9-atom primitive, 3x3x1) 那个 Se 写成 x=1.0 的经典情形：
pheasy 的 create_supercell (ase get_scaled_positions 默认 wrap=True) 把 1.0 折回 0，
于是 pheasy 超胞第 54-62 号原子的 x 是 [0, 1/3, 2/3]，而 phonopy 的 dataset 序是
[1/3, 2/3, 0] —— 正好差 3 组 3-循环。job 5830 因为 .pheasy_order.json 幂等跳过而
在错位数据上拟合（6.89% 而不是 0.578%）。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common" / "fcfit"))
import fc_fit_driver as d  # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures" / "mn2in2se5"


def _expected_perm():
    """The 3 three-cycles on atoms 54..62 (z=0 Se images)."""
    perm = np.arange(81)
    perm[54], perm[55], perm[56] = 56, 54, 55
    perm[57], perm[58], perm[59] = 59, 57, 58
    perm[60], perm[61], perm[62] = 62, 60, 61
    return perm


def test_perm_from_pheasy_sposcar_locks_3_cycles(tmp_path):
    pytest.importorskip("ase")
    from ase.io import read, write
    ds = read(FIX / "SPOSCAR.dataset", format="vasp")
    perm = _expected_perm()
    # pheasy's SPOSCAR = dataset order permuted (pheasy row i = dataset row perm[i])
    ph = ds[perm]
    write(tmp_path / "SPOSCAR", ph, format="vasp", direct=True)
    (tmp_path / "SPOSCAR.dataset").write_text((FIX / "SPOSCAR.dataset").read_text())
    got = d._perm_from_pheasy_sposcar(tmp_path)
    assert got is not None
    np.testing.assert_array_equal(got, perm)


def test_perm_identity_when_same_order(tmp_path):
    pytest.importorskip("ase")
    from ase.io import write, read
    ds = read(FIX / "SPOSCAR.dataset", format="vasp")
    write(tmp_path / "SPOSCAR", ds, format="vasp", direct=True)
    (tmp_path / "SPOSCAR.dataset").write_text((FIX / "SPOSCAR.dataset").read_text())
    got = d._perm_from_pheasy_sposcar(tmp_path)
    np.testing.assert_array_equal(got, np.arange(len(ds)))


def test_apply_pheasy_order_content_idempotence(tmp_path):
    """dataset 序 -> 换序；已经 pheasy 序 -> 跳过；都不是 -> 报错。"""
    import pickle
    n, natom = 4, 81
    perm = _expected_perm()
    rng = np.random.default_rng(0)
    ds_disps = rng.standard_normal((n, natom, 3))
    ds_forces = rng.standard_normal((n, natom, 3))
    np.save(tmp_path / "dataset_disps.npy", np.concatenate([ds_disps, ds_disps[:1]]))
    np.save(tmp_path / "dataset_forces.npy", np.concatenate([ds_forces, ds_forces[:1]]))

    # pheasy 序的 SPOSCAR（让 _perm_from_pheasy_sposcar 命中非恒等 perm）
    pytest.importorskip("ase")
    from ase.io import write, read
    ds = read(FIX / "SPOSCAR.dataset", format="vasp")
    write(tmp_path / "SPOSCAR", ds[perm], format="vasp", direct=True)
    (tmp_path / "SPOSCAR.dataset").write_text((FIX / "SPOSCAR.dataset").read_text())

    cfg = {"n_frames": n}

    def _dump(a, b):
        pickle.dump(a, open(tmp_path / "disp_matrix.pkl", "wb"))
        pickle.dump(b, open(tmp_path / "force_matrix.pkl", "wb"))

    # 1) dataset 序 -> 换序
    _dump(ds_disps, ds_forces)
    assert d.apply_pheasy_order(cfg, tmp_path) is True
    got_d = np.asarray(pickle.load(open(tmp_path / "disp_matrix.pkl", "rb")), float)
    np.testing.assert_array_equal(got_d, ds_disps[:, perm, :])

    # 2) 已经 pheasy 序 -> 跳过（内容判断，不因 .pheasy_order.json 存在而跳过）
    assert d.apply_pheasy_order(cfg, tmp_path) is False
    got_d = np.asarray(pickle.load(open(tmp_path / "disp_matrix.pkl", "rb")), float)
    np.testing.assert_array_equal(got_d, ds_disps[:, perm, :])

    # 3) 既不是 dataset 序也不是 pheasy 序 -> 报错
    _dump(ds_disps + 1.0, ds_forces + 1.0)
    with pytest.raises(SystemExit):
        d.apply_pheasy_order(cfg, tmp_path)


def test_fc_reorder_roundtrip():
    """apply_pheasy_fc_order 的 inv 换序与 perm 互逆（换过去再换回来 = 原 fc）。"""
    perm = _expected_perm()
    inv = np.argsort(perm)
    fc = np.random.default_rng(1).standard_normal((81, 81, 3, 3))
    # pheasy order -> dataset order (apply_pheasy_fc_order 的方向)
    fc_ds = np.take(np.take(fc, inv, axis=0), inv, axis=1)
    # dataset order -> pheasy order (逆)
    fc_ph = np.take(np.take(fc_ds, perm, axis=0), perm, axis=1)
    np.testing.assert_allclose(fc_ph, fc)
    # 换序确实动了那 9 个原子
    assert not np.allclose(fc_ds, fc)


def test_symfc_fc2_reorder_roundtrip_and_asr():
    """用真实 symfc fc2（dataset 序）验证：换序往返不破坏 fc，且满足 ASR。

    这是"已知力常数往返"护栏的核心：一个正确序的 fc 经过 perm -> inv 往返
    必须原样回来，否则下游 phono3py 拿到的就是错位力常数。
    """
    pytest.importorskip("h5py")
    import h5py
    with h5py.File(FIX / "symfc_fc2.hdf5", "r") as h:
        fc = np.asarray(h["force_constants"][()], float)
    assert fc.shape == (81, 81, 3, 3)
    # ASR（平移不变性）：Σ_j Φ_ij ≈ 0
    assert np.abs(fc.sum(axis=1)).max() < 1e-6
    # 换序往返：dataset -> pheasy (perm) -> dataset (inv) 应原样回来
    perm = _expected_perm()
    inv = np.argsort(perm)
    fc_ph = np.take(np.take(fc, perm, axis=0), perm, axis=1)
    fc_ds = np.take(np.take(fc_ph, inv, axis=0), inv, axis=1)
    np.testing.assert_allclose(fc_ds, fc, atol=1e-12)
    # 换序确实动了那 9 个原子（perm 非恒等）
    assert not np.allclose(fc_ph, fc)
