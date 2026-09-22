#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_postprocess_intrinsic.py —— postprocess_intrinsic.py 的自检。

用法：
    python test_postprocess_intrinsic.py
退出码 0 = 全部 PASS。

覆盖（全部离线、不需要 amset/pymatgen）：
  · split_mesh_key：AMSET write_mesh 的自旋 key bug（down 存成 <name>_up_down）解析；
  · expand_ir_to_full：不可约 -> 全网格展开；
  · zero_mechanisms / drop_mechanisms / resolve_mechanism_indices：机制选择；
  · inplane_average / to_jsonable；
  · 关键等价性：置零 IMP 后的机制求和 == 删除 IMP 后的机制求和（overall 等价）。

可选集成测试（需要 amset + h5py，且设环境变量 INTRINSIC_TEST_RUN_DIR=<S8.4运行目录>）：
  · 读真实 mesh.h5，用 AMSET 重积分复现 transport.json（max|diff| 应 <= 1e-6）。
"""
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import postprocess_intrinsic as P  # noqa: E402

FAILS = []


def ok(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ---------------- split_mesh_key ----------------
ok(P.split_mesh_key("energies_up") == ("energies", "up"), "split: energies_up -> (energies, up)")
ok(P.split_mesh_key("energies_down") == ("energies", "down"), "split: energies_down")
ok(P.split_mesh_key("energies_up_down") == ("energies", "down"),
   "split: energies_up_down -> (energies, down)  [AMSET write_mesh bug]")
ok(P.split_mesh_key("scattering_rates_up_down") == ("scattering_rates", "down"),
   "split: scattering_rates_up_down  [bug]")
ok(P.split_mesh_key("doping") == ("doping", None), "split: 无自旋后缀")
ok(P.split_mesh_key("kpoints") == ("kpoints", None), "split: kpoints")
# 不能把普通名字里的 down 误判（如以 _down 结尾才判定）
ok(P.split_mesh_key("vb_idx_up") == ("vb_idx", "up"), "split: vb_idx_up")

# ---------------- expand_ir_to_full ----------------
r_ir = np.arange(2 * 4).reshape(2, 4).astype(float)  # (2 labels, 4 ir)
mapping = np.array([0, 1, 0, 1, 2, 3, 3])
full = P.expand_ir_to_full(r_ir, mapping)
ok(full.shape == (2, 7), "expand: 形状 (2 labels, 7 full)")
ok(np.allclose(full[0], [0, 1, 0, 1, 2, 3, 3]), "expand: 值按映射取列")

# ---------------- 机制选择 ----------------
labels = ["ADP", "IMP", "POP"]
# 直接构造 (n_labels, n_dop, n_temp, n_bands, n_k) = (3,2,1,1,2)
base = np.zeros((3, 2, 1, 1, 2))
base[0] = 1.0   # ADP
base[1] = 5.0   # IMP
base[2] = 3.0   # POP
rates = {1: base.copy()}
ok(P.resolve_mechanism_indices(labels, ["IMP"]) == [1], "resolve: IMP -> [1]")
ok(P.resolve_mechanism_indices(labels, ["imp"]) == [1], "resolve: 大小写不敏感")
ok(P.resolve_mechanism_indices(labels, ["PIE"]) == [], "resolve: 不存在的机制 -> []")

zr, dropped, kept = P.zero_mechanisms(rates, labels, ["IMP"])
ok(np.all(zr[1][1] == 0.0), "zero: IMP 行全为 0")
ok(np.all(zr[1][0] == 1.0) and np.all(zr[1][2] == 3.0), "zero: ADP/POP 行不变")
ok(dropped == ["IMP"] and kept == ["ADP", "POP"], "zero: dropped/kept 名单")
ok(list(rates[1][1].ravel()[:2]) == [5.0, 5.0], "zero: 不改原数组（拷贝）")

dr, kept2, dropped2 = P.drop_mechanisms(rates, labels, ["IMP"])
ok(dr[1].shape == (2, 2, 1, 1, 2), "drop: 机制维变 2")
ok(np.all(dr[1][0] == 1.0) and np.all(dr[1][1] == 3.0), "drop: 剩下 ADP/POP")
ok(kept2 == ["ADP", "POP"] and dropped2 == ["IMP"], "drop: 名单")

# 关键等价性：sum(置零) == sum(删除)
s_zero = sum(zr[1][i] for i in range(len(labels)))
s_drop = sum(dr[1][i] for i in range(len(kept2)))
ok(np.allclose(s_zero, s_drop), "等价性: sum(置零 IMP) == sum(删除 IMP)  [overall]")

# ---------------- inplane_average ----------------
t = np.array([[[1.0, 0, 0], [0, 3.0, 0], [0, 0, 9.0]],
              [[2.0, 0, 0], [0, 4.0, 0], [0, 0, 8.0]]])
ia = P.inplane_average(t)
ok(np.allclose(ia, [2.0, 3.0]), "inplane_average: (xx+yy)/2")
try:
    P.inplane_average(np.zeros((2, 2)))
    ok(False, "inplane_average: 非 3x3 应报错")
except ValueError:
    ok(True, "inplane_average: 非 3x3 报错")

ok(P.to_jsonable(np.array([[1.5, 2.5]])) == [[1.5, 2.5]], "to_jsonable")

# ---------------- 可选：真实 mesh.h5 的 key 解析 ----------------
try:
    import h5py  # noqa: F401
    import tempfile
    from pymatgen.electronic_structure.core import Spin
    tmp = tempfile.mkdtemp(prefix="intrinsic_test_")
    fp = os.path.join(tmp, "mesh_test.h5")
    with h5py.File(fp, "w") as f:
        f.create_dataset("energies_up", data=np.ones((2, 3)))
        f.create_dataset("energies_up_down", data=np.full((2, 3), 2.0))
        f.create_dataset("scattering_rates_up_down", data=np.full((1, 3), 3.0))
        f.create_dataset("scattering_labels", data=np.array([b"IMP"]))
        f.create_dataset("doping", data=np.array([1.0]))
    md = P.read_mesh_h5(fp)
    ok(md["energies"][Spin.up].shape == (2, 3), "read_mesh_h5: up 通道就位")
    ok(np.all(md["energies"][Spin.down] == 2.0), "read_mesh_h5: up_down(bug) 正确归到 down")
    ok(np.all(md["scattering_rates"][Spin.down] == 3.0), "read_mesh_h5: rates down 归位")
    ok(md["scattering_labels"] == ["IMP"], "read_mesh_h5: labels 解码")
    ok(float(np.asarray(md["doping"])[0]) == 1.0, "read_mesh_h5: 标量键")
except ImportError:
    print("[SKIP] h5py/pymatgen 不在 —— 跳过 read_mesh_h5 集成用例")
except Exception as e:  # noqa: BLE001
    ok(False, "read_mesh_h5 用例异常: %r" % e)

# ---------------- 可选：真实运行目录的复现集成测试 ----------------
run_dir = os.environ.get("INTRINSIC_TEST_RUN_DIR")
if run_dir:
    print("[..] 集成测试：%s" % run_dir)
    rc = P.main([run_dir, "--check-reproduce", "--out", "_intrinsic_test.json"])
    ok(rc == 0, "集成: 重积分复现 transport.json (rc=%s)" % rc)

print("")
if FAILS:
    print("FAILED %d:" % len(FAILS))
    for m in FAILS:
        print("  - " + m)
    sys.exit(1)
print("ALL PASS")
