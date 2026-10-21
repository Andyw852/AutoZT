"""skill/ke-dft-cpu/tools/adp_k0_probe.py：几秒钟重演 ADP_DEGEN 矩阵路径在 k0 的对齐与核对（EPM 硅，不需要 AMSET）。

正确的三份文件（带符号 S 含参考能级常数、幅值 A 张量约定、矩阵 X 的剪切只有 +δ 靠对称操作补）：对齐读回常数、
Γ 三重态由 X 重演的 S、A 与文件一致（剪切的 S = 0、A = |n|,0,|n|），退出码 0；
ibands 错一条带 -> 对齐离散度超门槛，退出码 1；幅值文件少了张量标记 -> A 差一倍（打印里看得出），S 照样对上。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "skill" / "ke-dft-cpu" / "tools"), str(ROOT / "skill" / "ke-dft-cpu"),
           str(Path(__file__).resolve().parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
h5py = pytest.importorskip("h5py")
import _epm_si as E  # noqa: E402
import adp_k0_probe as P  # noqa: E402
import deform_wave as DW  # noqa: E402
import test_0061m_deform_wave as TW  # noqa: E402


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("probe")
    root, dirs, (kpts, _dat) = TW._s7(tmp, n=4, plus_only=[(1, 2), (0, 2), (0, 1)])
    rd = tmp / "run"
    rd.mkdir()
    DW.build(root, dirs, rd / "deformation_matrix.h5", log=lambda m: None)
    NB = TW.NB
    S = np.zeros((NB, len(kpts), 3, 3))
    A = np.zeros((NB, len(kpts), 3, 3))
    e0 = np.array([E.solve(k, E.L0, TW.ENCUT)[1][:NB] for k in kpts])
    for ip, Pt in enumerate(E.PATS):
        v = []
        for s in (1, -1):
            es = np.array([E.solve(k, E.lattice(E.pattern(Pt, s * TW.DELTA)), TW.ENCUT)[1][:NB] for k in kpts])
            v.append((e0 - es) / (s * TW.DELTA) + TW.REF_SHIFT[ip])
        a, b = Pt
        S[:, :, a, b] = S[:, :, b, a] = (0.5 * (v[0] + v[1])).T
        A[:, :, a, b] = A[:, :, b, a] = (0.5 * (np.abs(v[0]) + np.abs(v[1]))).T / (1.0 if a == b else 2.0)
    with h5py.File(rd / "deformation_signed.h5", "w") as f:
        f["signed_up"] = S
        f["kpoints"] = kpts
        f["ibands_up"] = np.arange(NB)
        f.attrs["strain_amount"] = TW.DELTA
    with h5py.File(rd / "deformation.h5", "w") as f:
        f["deformation_potentials_up"] = A
        f["kpoints"] = kpts
        f.attrs["dp_offdiag"] = "tensor"
    return rd


def test_probe_consistent(run, capsys):
    assert P.main(["--dir", str(run)]) == 0
    log = capsys.readouterr().out
    assert "子空间 3 条" in log and "★" not in log
    xy = next(ln for ln in log.splitlines() if ln.strip().startswith("xy：eig(X)"))
    s0 = [float(x) for x in xy.split("S0 [")[1].split("]")[0].split(",")]
    a0 = [float(x) for x in xy.split("A0 [")[1].split("]")[0].split(",")]
    assert np.ptp(s0) < 0.05 and sorted(np.round(a0, 0))[0] == 0.0 and min(a0[0], a0[2]) > 0.3   # 剪切：三条带只差参考常数
    # 简并点上逐带的 S 带二阶污染（xx：首末两条带不对称），子空间内求和不带
    xx = next(ln for ln in log.splitlines() if ln.strip().startswith("xx：eig(X)"))
    s0 = [float(x) for x in xx.split("S0 [")[1].split("]")[0].split(",")]
    sp = [float(x) for x in xx.split("由 X 重演 [")[1].split("]")[0].split(",")]
    assert abs(s0[0] - s0[2]) > 0.5 and abs(sp[0] - sp[2]) < 0.05 and abs(sum(s0) - sum(sp)) < 0.05
    rel = log.split("k0 相对残差（子空间内求和")[1]
    assert float(rel.split("：")[1].split("%")[0]) < 1.0


def test_probe_flags_wrong_band_offset(run, tmp_path, capsys):
    import shutil
    rd = tmp_path / "bad"
    shutil.copytree(run, rd)
    with h5py.File(rd / "deformation_signed.h5", "r+") as f:
        del f["ibands_up"]
        f["ibands_up"] = np.arange(1, TW.NB + 1)
    assert P.main(["--dir", str(rd)]) == 1
    assert "超过门槛" in capsys.readouterr().out
