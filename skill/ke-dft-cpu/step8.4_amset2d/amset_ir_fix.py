"""
amset_ir_fix.py —— AMSET 密网格不可约 k 点改用正确的晶格（运行时插件，不改 AMSET 源码）

用法（S8 / S8.4 的 gen 已把本文件复制进运行目录，并在 python -c 里最先 import）：
    AZ_IR_FIX=1 python -c "import amset_ir_fix; ...; Runner.from_directory('.').run()"

开关：环境变量 AZ_IR_FIX（未设 / 1/true/on/yes = 修；0/false/off/no = AMSET 原行为）。
    gen 按 step.conf 的 IR_FIX（auto/on/off，auto = on）在命令前写 export AZ_IR_FIX=…。

------------------------------------------------------------------------------
问题（V125）：amset/electronic_structure/kpoints.py::get_kpoints_tetrahedral 给 spglib 的是
    cell = (np.array(atoms.get_cell().T, ...),   # lattice
            atoms.get_scaled_positions(), ...)
ASE 的 get_cell() 按行存晶格矢量，spglib 的 cell[0] 也要求按行 —— 多出来的 .T 等于把晶格换成了
另一个。spglib 于是只找得到两种晶格都容许、又把原子映回原子的操作：
  · 晶格矩阵对称（立方、fcc/bcc 原胞、正交/四方常规胞）-> 没影响；
  · 六方/三方（2D TMD、纤锌矿 GaN/AlN…）-> 只剩 {E, σh}×时间反演。
找到的都是真对称操作，所以**结果是对的**，只是不可约 k 点多了 5.5–5.9 倍（实测：WSe2 单层
263×263×33 587945 vs 100232；GaN 61×61×37 35359 vs 6479）。散射按不可约点逐个算，于是白算 5–6 倍。
合成六方单层 AMSET 端到端（ADP+IMP，85×85×11）：不可约点 21678 -> 3870，散射 121 s -> 25 s，
迁移率一致到 0.05%。

修法：只替换 kpoints 模块里 get_ir_reciprocal_mesh 的调用 —— 把收到的 lattice 再转置一次（= 行矢量）。
该模块只有 get_kpoints_tetrahedral 这一处调用 spglib（护栏里核对）。
护栏：源码里必须能找到 "atoms.get_cell().T" 且 spglib 只被调用这一处；对不上（AMSET 已修或改了
实现）-> 不打补丁、照原样运行（原行为结果是对的，只是慢），并在 stderr 说明。
------------------------------------------------------------------------------
"""

import inspect
import os
import sys

import numpy as np

ENV_NAME = "AZ_IR_FIX"
_OFF = ("0", "false", "off", "no")
BUG_SNIPPET = "atoms.get_cell().T"
STATE = {"applied": False, "calls": 0, "last": None}


def enabled():
    return os.environ.get(ENV_NAME, "").strip().lower() not in _OFF


def _squash(s):
    return "".join(str(s).split())


class _SpglibProxy:
    """只改 get_ir_reciprocal_mesh：lattice 转回行矢量；其余属性原样转发。"""

    def __init__(self, real):
        self._real = real

    def __getattr__(self, name):
        return getattr(self._real, name)

    def get_ir_reciprocal_mesh(self, mesh, cell, *args, **kwargs):
        lat = np.ascontiguousarray(np.asarray(cell[0], dtype="double").T)
        fixed = (lat,) + tuple(cell[1:])
        mapping, grid = self._real.get_ir_reciprocal_mesh(mesh, fixed, *args, **kwargs)
        STATE["calls"] += 1
        try:
            m0, _ = self._real.get_ir_reciprocal_mesh(mesh, cell, *args, **kwargs)
            n_old = int(len(np.unique(m0)))
        except Exception:                                  # noqa: BLE001  只用于日志
            n_old = None
        n_new = int(len(np.unique(mapping)))
        STATE["last"] = {"mesh": [int(x) for x in mesh], "n_full": int(np.prod(mesh)),
                         "ir_fixed": n_new, "ir_amset_original": n_old}
        print("[ir_fix] mesh %s：不可约 k 点 %d（AMSET 原式 %s，少算 x%s）"
              % ("x".join(str(int(x)) for x in mesh), n_new, n_old,
                 ("%.1f" % (n_old / n_new)) if n_old and n_new else "?"),
              file=sys.stderr, flush=True)
        return mapping, grid


def apply(force=False):
    """给 amset.electronic_structure.kpoints 换上 spglib 代理。返回是否打了补丁。"""
    if not (force or enabled()):
        print("[ir_fix] %s=off -> 不打补丁（AMSET 原不可约 k 点）" % ENV_NAME, file=sys.stderr, flush=True)
        return False
    import amset
    import amset.electronic_structure.kpoints as _kp

    if isinstance(_kp.spglib, _SpglibProxy):
        return True
    try:
        fsrc = inspect.getsource(_kp.get_kpoints_tetrahedral)
        msrc = inspect.getsource(_kp)
    except (OSError, TypeError) as e:
        print("[ir_fix] 读不到 AMSET 源码（%s）-> 不打补丁（原行为：结果对，只是慢）" % e,
              file=sys.stderr, flush=True)
        return False
    n_calls = _squash(msrc).count("spglib.")
    if _squash(BUG_SNIPPET) not in _squash(fsrc) or n_calls != 1:
        print("[ir_fix] amset %s 的 get_kpoints_tetrahedral 与已核对的实现不同（%s 出现=%s，spglib 调用 %d 处）"
              " -> 不打补丁（原行为：结果对，只是可能慢）"
              % (getattr(amset, "__version__", "?"), BUG_SNIPPET, _squash(BUG_SNIPPET) in _squash(fsrc), n_calls),
              file=sys.stderr, flush=True)
        return False
    _kp.spglib = _SpglibProxy(_kp.spglib)
    STATE["applied"] = True
    print("[ir_fix] amset %s：密网格不可约 k 点改用正确晶格（行矢量）" % getattr(amset, "__version__", "?"),
          file=sys.stderr, flush=True)
    return True


apply()
