# -*- coding: utf-8 -*-
"""nspin_norm_fix.py —— AMSET <0.5.1 在 ISPIN=2 时形变势被自旋道数除的补丁。

背景（2026-09-21 确证）：
    AMSET 0.4.19 的 amset/deformation/potentials.py 里
    calculate_deformation_potentials() 把 norm += strain_loc 写在
    for spin, spin_deform in deform.items(): 循环**体内**。ISPIN=2 体系有
    两份自旋道，同一个应变条目被累加两次，随后
    deformation_potentials[spin] /= norm 把每个 D 除以「2 × 自旋道数」，
    也就是**所有形变势 D 整体减半**。上游 AMSET 0.5.1 已把 norm 累加移出
    自旋循环（修复）。

    判据是 VASP INCAR 的 ISPIN=2：只要算的是自旋极化（两份 h5 自旋道），
    哪怕磁矩 ~ 0 也必须修，因为两份道就是两份数据。实测 CrSe2_hex：
    deformation_vac.h5 的 E1_vac=2.4465 eV，而独立路径（原始本征值 +
    LOCPOT 真空对齐）得到的 band_edges.json E1_vac=4.8808 eV，差约 2 倍。

设计：把「判定」（纯函数，可离线用 stdlib 测试）与「改 h5」（h5py I/O）分开。
    * decide(amset_version, nspin, marked) 决定是否乘回自旋道数；
    * fix_h5(path, amset_version) 原地修正一份 h5，并在 attrs 写幂等标记。

幂等：第一次 fix_h5 后 attrs 有 nspin_norm_fixed，再次调用直接返回
(False, 'already-marked')，绝不重复乘。无论修不修，都会写
nspin_norm_fixed / nspin_norm_fix_version / nspin_norm_fix_reason，
方便回看每份 h5 是怎么处理的。
"""
import re

__all__ = [
    "parse_version", "version_lt", "decide", "fix_h5",
    "relative_diff", "gate_ok", "DATASET_PREFIX",
]

DATASET_PREFIX = "deformation_potentials_"


def parse_version(s):
    """把版本串解析成数字段元组（只比较数字段）。

    '0.4.19' -> (0, 4, 19)。每段只取开头连续数字，遇到非数字段就停：
    '0.5.1.dev0+g123' -> (0, 5, 1)、'0.5.1rc1' -> (0, 5, 1)。
    空串/None 返回 (0,)。
    """
    parts = []
    for seg in str(s if s is not None else "").strip().split("."):
        m = re.match(r"\s*(\d+)", seg)
        if not m:
            break
        parts.append(int(m.group(1)))
    return tuple(parts or [0])


def version_lt(a, b):
    """严格小于比较（只比较数字段，短的后补 0）：0.5 < 0.5.1 -> True。"""
    ta, tb = parse_version(a), parse_version(b)
    n = max(len(ta), len(tb))
    ta = ta + (0,) * (n - len(ta))
    tb = tb + (0,) * (n - len(tb))
    return ta < tb


def decide(amset_version, nspin, marked):
    """返回 (should_fix, reason)。

    * marked 为真                 -> (False, 'already-marked')
    * version < 0.5.1 且 nspin==2 -> (True,  'amset<0.5.1 & nspin=2')
    * version >= 0.5.1            -> (False, 'amset>=0.5.1')
    * 其余（nspin != 2）          -> (False, 'nspin=%d' % nspin)
    """
    if marked:
        return False, "already-marked"
    if version_lt(amset_version, "0.5.1") and int(nspin) == 2:
        return True, "amset<0.5.1 & nspin=2"
    if not version_lt(amset_version, "0.5.1"):
        return False, "amset>=0.5.1"
    return False, "nspin=%d" % int(nspin)


def fix_h5(path, amset_version, dataset_prefix=DATASET_PREFIX):
    """原地修正一份形变势 h5，返回 (fixed, reason)。

    以 'a' 打开。若 attrs 里已有 nspin_norm_fixed，直接返回
    (False, 'already-marked')（不重复乘）。否则统计
    <dataset_prefix>* datasets 数作为自旋道数 nspin，按 decide 判定；
    should_fix 时对每个 dataset 原地乘 nspin。

    总是写属性：nspin_norm_fixed（True 或字符串 'skipped'）、
    nspin_norm_fix_version=amset_version、nspin_norm_fix_reason=reason；
    修复时另写 nspin_norm_fix_factor。
    """
    import h5py

    with h5py.File(str(path), "a") as fh:
        if "nspin_norm_fixed" in fh.attrs:
            return False, "already-marked"
        names = [k for k in fh.keys() if k.startswith(dataset_prefix)]
        nspin = len(names)
        should_fix, reason = decide(amset_version, nspin, marked=False)
        if should_fix:
            for name in names:
                d = fh[name]
                d[...] = d[...] * nspin
            fh.attrs["nspin_norm_fixed"] = True
            fh.attrs["nspin_norm_fix_factor"] = nspin
        else:
            fh.attrs["nspin_norm_fixed"] = "skipped"
        fh.attrs["nspin_norm_fix_version"] = str(amset_version)
        fh.attrs["nspin_norm_fix_reason"] = reason
        return should_fix, reason


def relative_diff(a, b):
    """相对差 abs(a-b)/abs(b)；b==0 时 a==0 返回 0.0，否则 inf。"""
    a, b = float(a), float(b)
    denom = abs(b)
    if denom == 0.0:
        return 0.0 if a == 0.0 else float("inf")
    return abs(a - b) / denom


def gate_ok(h5_val, be_val, tol=0.05):
    """硬闸门判据：h5 值对 band_edges 参考值的相对差 <= tol（默认 5%）通过。"""
    return relative_diff(h5_val, be_val) <= tol
