# -*- coding: utf-8 -*-
"""[0060] 三处"只影响打印/找谷"的修正（不碰物理计算）：

1. `postprocess_intrinsic.py` 打印"本征迁移率"取的是 mobility 数组索引 0，
   而索引 0 是 |掺杂| 最大的一档（−1e21），不是本征点 → 改取"最低非零掺杂档"。
2. S8.3 汇总 `[AMSET 机制]` 把 AMSET 用来表示"该机制散射率≈0"的极大 mobility
   （1e6–1e9）当实数印出来 → 超过阈值印 "—"，且不把它当"限制因素"。
3. `odp_frozen.find_valleys_3d` 的平台合并会跨 k 与 −k，把时间反演等价的两个谷
   （Si 的 ±Δ 谷在 11³ 网格上相邻、能量相等）并成一个 → 6 个谷只剩 3 个。
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
for _d in (str(ROOT), str(KE), str(KE / "step8.4_amset2d"), str(ROOT / "skill" / "_common" / "opt")):
    if _d not in sys.path:
        sys.path.insert(0, _d)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _pp():
    return _load("pp_0060", KE / "step8.4_amset2d" / "postprocess_intrinsic.py")


def _g13():
    return _load("g13_0060", KE / "step8.3_output" / "gen_step13_output.py")


# ---------------------------------------------------------------- 1) 掺杂档索引
def test_display_doping_index_picks_lowest_nonzero():
    P = _pp()
    dop = [-1e21, -1e20, 1e18, 0.0, -1e10, 1e10]
    assert P.display_doping_index(dop) == 4              # 1e10 绝对值最小；±1e10 并列时取先出现的 -1e10
    assert P.display_doping_index([-1e21, -1e10, 0.0]) == 1
    assert P.display_doping_index([0.0, 0.0]) == 0       # 全是 0（异常输入）退回 0
    assert P.display_doping_index([]) == 0
    assert P.display_doping_index([1e18, -1e17]) == 1    # 负档也按绝对值比


# ---------------------------------------------------------------- 2) 可忽略机制显示
def _mat(v):
    return [[v, 0.0, 0.0], [0.0, v, 0.0], [0.0, 0.0, v]]


def _series(vals):
    """[掺杂档] -> mobility 数组形状 (ndop, ntemp=1, 3, 3)。"""
    return [[_mat(float(v))] for v in vals]


def test_amset_mech_lines_hides_inactive(tmp_path):
    g = _g13()
    run = tmp_path / "step8_amset"
    run.mkdir()
    # 掺杂档顺序与 AMSET 一致：重掺杂在前，最低非零档在中间，p 侧在后
    tj = {
        "doping": [-1e21, -1e10, 1e10],
        "temperatures": [300.0],
        "mobility": {
            "overall": _series([64.0, 1459.0, 541.0]),
            "IVS": _series([483.0, 2747.0, 2.016e9]),      # p 侧没有 IVS 通道 -> 极大值
            "ADP": _series([2341.0, 3883.0, 7610.0]),
            "IMP": _series([80.0, 8.331e6, 6.356e6]),      # 低掺杂 -> 杂质散射可忽略
        },
    }
    (run / "transport.json").write_text(json.dumps(tj), encoding="utf-8")
    lines = g._amset_mech_lines(tmp_path, is_2d=False)
    assert len(lines) == 2
    el, ho = lines
    # 电子：n = 1e10 档，IVS 2747 最小 -> 限制因素 IVS；IMP 可忽略印 "—"
    assert "2747" in el and "3883" in el and "限制因素 IVS" in el and "IMP —" in el
    assert "可忽略/未参与" in el
    # 空穴：p 侧没有 IVS -> "IVS —"，限制因素 ADP
    assert "IVS —" in ho and "限制因素 ADP" in ho
    assert "2.016e+09" not in ho and "8.331e+06" not in el   # 极大值不再当实数印出


# ---------------------------------------------------------------- 3) 平台不跨 ±k
def test_find_valleys_3d_keeps_time_reversed_valleys(tmp_path):
    of = _load("of_0060", KE / "odp_frozen.py")
    N = 11
    grid = np.arange(N) / float(N)
    kpts = [[x, y, z] for x in grid for y in grid for z in grid]
    k = np.asarray(kpts, float)
    centers = np.asarray([[5 / 11, 0.0, 0.0], [6 / 11, 0.0, 0.0],
                          [0.0, 5 / 11, 0.0], [0.0, 6 / 11, 0.0],
                          [0.0, 0.0, 5 / 11], [0.0, 0.0, 6 / 11]], float)

    def _periodic_dist(kk, c):
        d = kk - c
        d -= np.round(d)
        return float(np.linalg.norm(d))

    e = np.asarray([min(_periodic_dist(kk, c) for c in centers) for kk in k])
    eig = {"up": np.column_stack([np.full(len(k), -1.0), e])}
    vs = of.find_valleys_3d(kpts, eig, {"up": 1}, +1, np.eye(3), window=0.30)
    # 6 个等价极小在同一颗星里（g=6）；旧版平台合并会只剩 3 个
    assert len(vs) == 1 and vs[0]["g"] == 6, vs
    tr = [tuple(np.round(x, 4)) for x in vs[0]["k_frac"]]
    # k_frac 采用 (−½, ½] 约定：5/11 → 0.4545，6/11 → −0.4545，两个时间反演伙伴都在
    assert (0.4545, 0.0, 0.0) in tr and (-0.4545, 0.0, 0.0) in tr
