# -*- coding: utf-8 -*-
"""kl-dft-cpu S6：κ∞ 外推回归（build_extract 生成的真实脚本 + 合成 kappa-m*.hdf5）。

背景：外推曾用 mesh.split()[0] 当 N。phono3py 原胞重排后 2D 网格是 '1 N N'，
第一个数恒为 1 -> 外推退化（真值 100 给出 45.9）。两种轴序必须给出同一个 κ∞。
"""
import json
import os
import subprocess
import sys

import pytest

np = pytest.importorskip("numpy")
h5py = pytest.importorskip("h5py")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in ("skill/_common", "skill/_common/opt", "skill/kl-dft-cpu"):
    sys.path.insert(0, os.path.join(ROOT, _p))

g6 = pytest.importorskip("gen_step6_kappa")

KINF, A = 100.0, 900.0          # 合成 κ(N) = KINF - A/N（300 K 面内）
T = np.arange(100.0, 900.0, 100.0)


def _run(tmp_path, meshes):
    plan = []
    for m in meshes:
        n = max(int(x) for x in m.split())
        k300 = KINF - A / n
        K = np.zeros((len(T), 6))
        K[:, 0] = K[:, 1] = k300 * 300.0 / T
        K[:, 2] = 0.01
        with h5py.File(tmp_path / ("kappa-m%s.hdf5" % "".join(m.split())), "w") as h:
            h["temperature"] = T
            h["kappa"] = K
        plan.append((m, "rta"))
    script = g6.build_extract(1.0, {}, None, g6.extract_plan(plan, False),
                              meshes[-1] + "|rta")
    (tmp_path / "x.sh").write_text(script, encoding="utf-8")
    r = subprocess.run(["bash", "x.sh"], cwd=tmp_path, capture_output=True, text=True,
                       env=dict(os.environ, PATH=os.path.dirname(sys.executable)
                                + os.pathsep + os.environ.get("PATH", "")))
    assert "KAPPA_DONE" in r.stdout, r.stdout + r.stderr
    d = json.loads((tmp_path / "kappa_summary.json").read_text(encoding="utf-8"))
    return d["mesh_convergence"]["rta"]


@pytest.mark.parametrize("meshes", [
    ["90 90 1", "114 114 1", "138 138 1"],     # POSCAR 轴序
    ["1 90 90", "1 114 114", "1 138 138"],     # phono3py 原胞轴序（真空在 a1）
])
def test_kappa_inf_independent_of_axis_order(tmp_path, meshes):
    mc = _run(tmp_path, meshes)
    assert mc["kappa_inf_extrap_300K_inplane"] == pytest.approx(KINF, rel=1e-9)
    # 相邻档 <5% 会判"收敛"，但最细档离 κ∞ 仍 6.5% —— 新判据必须抓住
    assert mc["converged_5pct"] is True
    assert mc["finest_dev_from_kappa_inf_pct"] == pytest.approx(
        100.0 * (A / 138) / KINF, rel=1e-6)   # |κ(138)−κ∞|/κ∞
    assert mc["converged_extrap_3pct"] is False


def test_duplicate_n_is_not_extrapolated(tmp_path):
    # 两档面内 N 相同（仅轴序不同）：1/N 重复，不应给出外推值
    mc = _run(tmp_path, ["1 90 90", "90 90 1", "138 138 1"])
    assert "kappa_inf_extrap_300K_inplane" not in mc
    assert "extrap_error" in mc
