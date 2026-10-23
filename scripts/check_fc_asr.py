#!/usr/bin/env python3
"""check_fc_asr.py —— 只读检查 phono3py fc2/fc3 的声学求和规则（ASR）满足程度。

用法：
    python3 scripts/check_fc_asr.py <目录> [<目录> ...]
其中每个 <目录> 下应有 phono3py 格式的 fc2.hdf5 / fc3.hdf5
（kl-dft-cpu 里就是 step5_fc/cutoff_scan/cut3_<tag>/ 或 step6_kappa/cut3_<tag>/）。

检查内容（fc 为 full 形式时）：
    fc2 ASR:  sum_j   Phi2[i,j,a,b]       == 0     （平移原子 j）
    fc3 ASR:  sum_k   Phi3[i,j,k,a,b,c]   == 0     （平移原子 k）
              sum_jk  Phi3[i,j,k,a,b,c]   == 0     （两个平移都消掉）

输出 rel = max|sum| / max|fc|（无量纲）：
    ~1e-10        正常（机器精度）
    1e-6 .. 1e-4  偏松，低频声学支可能被虚加散射
    >=1e-4        基本可判定 ASR 没满足好

依赖：numpy、h5py（不需要 phono3py）。
本脚本只读文件、不写任何东西。
"""
import os
import sys

import numpy as np
import h5py


def load_h5(path):
    """读 hdf5 里的全部数据集（小文件，2-3 MB）。"""
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            try:
                out[k] = f[k][()]
            except Exception:
                pass
    return out


def check(d2=None, d3=None):
    res = {}
    if d2 is not None:
        fc2 = np.asarray(d2)
        if fc2.ndim == 4:
            s = fc2.sum(axis=1)                     # (N,3,3)
            mx = max(float(np.abs(fc2).max()), 1e-30)
            res["fc2"] = dict(shape=list(fc2.shape),
                              absmax=float(np.abs(fc2).max()),
                              sum_absmax=float(np.abs(s).max()),
                              rel=float(np.abs(s).max() / mx))
        else:
            res["fc2"] = dict(shape=list(fc2.shape), note="unexpected ndim")
    if d3 is not None:
        fc3 = np.asarray(d3)
        if fc3.ndim == 6:
            mx = max(float(np.abs(fc3).max()), 1e-30)
            s_k = fc3.sum(axis=2)
            s_jk = fc3.sum(axis=(1, 2))
            res["fc3"] = dict(shape=list(fc3.shape), absmax=float(np.abs(fc3).max()),
                              rel_k=float(np.abs(s_k).max() / mx),
                              rel_jk=float(np.abs(s_jk).max() / mx))
        else:
            res["fc3"] = dict(shape=list(fc3.shape), note="unexpected ndim")
    return res


def main(argv):
    if not argv:
        sys.exit(__doc__)
    bad = 0
    for root in argv:
        print("=" * 84)
        print("DIR:", root)
        f2 = os.path.join(root, "fc2.hdf5")
        f3 = os.path.join(root, "fc3.hdf5")
        d2 = load_h5(f2).get("fc2") if os.path.isfile(f2) else None
        d3 = load_h5(f3).get("fc3") if os.path.isfile(f3) else None
        if d2 is None and d3 is None:
            print("  [WARN] 没找到 fc2.hdf5/fc3.hdf5")
            bad += 1
            continue
        r = check(d2, d3)
        if "fc2" in r:
            x = r["fc2"]
            if "rel" in x:
                print("  fc2 shape=%-22s max|fc2|=%.4e  max|sum_j|=%.4e  rel=%.3e"
                      % (x["shape"], x["absmax"], x["sum_absmax"], x["rel"]))
            else:
                print("  fc2 shape=%s  %s" % (x["shape"], x.get("note")))
        if "fc3" in r:
            x = r["fc3"]
            if "rel_k" in x:
                print("  fc3 shape=%-22s max|fc3|=%.4e  rel(sum_k)=%.3e  rel(sum_jk)=%.3e"
                      % (x["shape"], x["absmax"], x["rel_k"], x["rel_jk"]))
            else:
                print("  fc3 shape=%s  %s" % (x["shape"], x.get("note")))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
