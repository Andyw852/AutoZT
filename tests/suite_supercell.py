#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""超胞规格自测（v1.0）：对角 / 一般矩阵 / 校验 / 2D 约束 / phono3py 语义一致性。

纯标准库（本机没有 numpy/phonopy，数值部分在集群上单独核）。
"""
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OK, FAIL = [], []


def ck(cond, msg):
    (OK if cond else FAIL).append(msg)
    print("  %s %s" % ("PASS" if cond else "FAIL", msg))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def raises(fn, *a, **kw):
    try:
        fn(*a, **kw)
    except SystemExit:
        return True
    return False


# klmace_common / kl_common 会 import 同目录的公共模块（gen 时会一起推到超算）
for _d in ("skill/_common/mace", "skill/_common/opt", "skill/phonon-dft-cpu"):
    sys.path.insert(0, os.path.join(ROOT, _d))

kc = load("kc_mace", os.path.join(ROOT, "skill", "_common", "mace", "klmace_common.py"))
kp = load("kc_phonon", os.path.join(ROOT, "skill", "phonon-dft-cpu", "kl_common.py"))

for tag, m in (("mace.klmace_common", kc), ("phonon.kl_common", kp)):
    print("== %s ==" % tag)
    ck(m.parse_reps("3 3 3", "3d") == [3, 3, 3], "对角：'3 3 3' → [3,3,3]")
    ck(m.dim_str([3, 3, 3]) == "3 3 3", "对角回写：'3 3 3'")
    M = m.parse_reps("2 1 0 -1 2 0 0 0 1", "3d")
    ck(m.is_matrix(M) and M == [[2, 1, 0], [-1, 2, 0], [0, 0, 1]],
       "一般矩阵：9 个数 → 3×3（行主序）")
    ck(m.dim_str(M) == "2 1 0 -1 2 0 0 0 1",
       "矩阵回写与 phonopy/phono3py 的 --dim 写法逐字一致")
    ck(m.sc_matrix([4, 2, 1]) == [[4, 0, 0], [0, 2, 0], [0, 0, 1]],
       "sc_matrix：对角 → 对角阵")
    ck(m.sc_matrix(M) == M, "sc_matrix：矩阵原样透传（不转置）")
    ck(m._det3([[2, 1, 0], [-1, 2, 0], [0, 0, 1]]) == 5, "行列式 = 5（扩胞体积比）")
    ck(raises(m.parse_reps, "1 0 0 0 1 0 0 0 0", "3d"),
       "行列式为 0 → 报错（不是合法扩胞）")
    ck(raises(m.parse_reps, "2 2 2 2 2", "3d"), "个数不对（5 个）→ 报错")
    ck(raises(m.parse_reps, "2 x 2", "3d"), "非整数 → 报错")
    ck(m.parse_reps("3 3 3", "2d", 2) == [3, 3, 1], "2D 对角：真空方向强制 1")
    ck(raises(m.parse_reps, "2 1 0 -1 2 0 0 0 2", "2d", 2),
       "2D：矩阵动了真空方向 → 报错（不静默算错）")
    ck(m.parse_reps("2 1 0 -1 2 0 0 0 1", "2d", 2) == [[2, 1, 0], [-1, 2, 0], [0, 0, 1]],
       "2D：面内混、真空保持 [0,0,1] → 放行")

print("\n结果：PASS %d，FAIL %d" % (len(OK), len(FAIL)))
for x in FAIL:
    print("  ! %s" % x)
sys.exit(1 if FAIL else 0)
