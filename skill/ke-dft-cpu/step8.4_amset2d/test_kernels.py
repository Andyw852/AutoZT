#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_kernels.py —— amset2d 二维散射核的自检（需 amset + pymatgen 环境）。

用法（在 skill/ke-dft-cpu/step8.4_amset2d 目录里）：
    python test_kernels.py            # 用内置示例 2d_correction.json
退出码 0 = 全部 PASS。

覆盖：
  · ADP：F = D²/C11（各向同性）、与 C44/C55 符号无关、各向异性张量绕 z 旋转不变、
         未重排 VASP 顺序时被护栏拦下；
  · POP：真空不变性（r 不变时 F ∝ c）、q→0 有限且 = 2πc·Δr/ε_env²、屏蔽生效、旋转不变；
  · IMP：与手算公式一致、q→0 不发散（有屏蔽）；
  · PIE：q→0 有限、小 q 近似常数、旋转不变；
  · gen_step14 的 Voigt 重排：行列同时重排，C66(XY) 落到标准 Voigt 第 6 位。
"""
import importlib.util
import json
import os
import sys
import tempfile
import types
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FAILS = []


def ok(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ---------- 在临时目录里放示例 2d_correction.json，再 import 插件 ----------
tmp = tempfile.mkdtemp(prefix="amset2d_test_")
rec_src = HERE / "example_2d_correction.json"
if not rec_src.is_file():
    sys.exit("[ERROR] 缺 %s" % rec_src)
rec0 = json.loads(rec_src.read_text())
(tmp / "2d_correction.json") if False else None
Path(tmp, "2d_correction.json").write_text(json.dumps(rec0))
os.chdir(tmp)
sys.path.insert(0, str(HERE))
import amset2d_plugin as P  # noqa: E402

rng = np.random.default_rng(1)


def rand_q(n=500, zfrac=0.3):
    q = rng.normal(size=(n, 3))
    q[:, 2] *= zfrac
    nq = rng.uniform(1e-4, 0.3, n)          # |q| in bohr^-1
    return q / np.linalg.norm(q, axis=1)[:, None], nq ** 2


def Rz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rot4(Cf, R):
    return np.einsum("ia,jb,kc,ld,abcd->ijkl", R, R, R, R, Cf)


class Interp:
    def __init__(self, D):
        self.D = D

    def interpolate(self, spin, b, k):
        return self.D[None]


try:
    from pymatgen.core.tensors import Tensor
except Exception as e:                                            # noqa: BLE001
    sys.exit("[ERROR] 需要 pymatgen：%s" % e)

Cv = np.zeros((6, 6))
Cv[0, 0] = Cv[1, 1] = 130
Cv[0, 1] = Cv[1, 0] = 30
Cv[5, 5] = 50                       # 标准 Voigt：XY 在第 6 位
Cv[2, 2] = 0.5
Cv[3, 3] = Cv[4, 4] = -0.02         # 面外近零甚至为负（slab 的典型病态）
Cf = np.array(Tensor.from_voigt(Cv))

a = P.AcousticDeformation2D.__new__(P.AcousticDeformation2D)
a.elastic_constant = Cf
a.fermi_levels = np.zeros((2, 3))
a.is_metal = False
Dd = 5.0
a.deformation_potential = Interp(np.diag([Dd, Dd, 0.0]))
u, nsq = rand_q()
F = a.factor(u, nsq, None, 0, None, np.zeros(3))
ok(np.allclose(F, Dd ** 2 / 130),
   "ADP 各向同性：F = D²/C11 = %.5f，且与 C44/C55 的符号无关（最大偏差 %.1e）"
   % (Dd ** 2 / 130, np.abs(F - Dd ** 2 / 130).max()))

Ca = Cv.copy()
Ca[1, 1] = 80
Ca[0, 5] = Ca[5, 0] = 7
Ca[1, 5] = Ca[5, 1] = -4
Cfa = np.array(Tensor.from_voigt(Ca))
Da = np.array([[4.0, 1.2, 0], [1.2, 2.5, 0], [0, 0, 0]])
a.elastic_constant = Cfa
a.deformation_potential = Interp(Da)
F1 = a.factor(u, nsq, None, 0, None, np.zeros(3))
R = Rz(0.7)
a.elastic_constant = rot4(Cfa, R)
a.deformation_potential = Interp(R @ Da @ R.T)
F2 = a.factor(u @ R.T, nsq, None, 0, None, np.zeros(3))
ok(np.allclose(F1, F2), "ADP 各向异性张量：整体绕 z 旋转后结果不变（斜方/单斜可用）")

# 未重排的 VASP 顺序：C66(XY) 被换成近零的 C_zxzx -> 护栏必须拦下
Vv = Cv.copy()
Vv[3, 3] = Vv[4, 4] = 0.05
perm = [0, 1, 2, 5, 3, 4]          # 标准 -> VASP 顺序
Vv = Vv[np.ix_(perm, perm)]
a.elastic_constant = np.array(Tensor.from_voigt(Vv))
a.deformation_potential = Interp(Da)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        F3 = a.factor(u, nsq, None, 0, None, np.zeros(3))
        ok(False, "VASP 顺序未重排时护栏未触发（F 中位数 %.3g）" % np.median(F3))
    except ValueError as e:
        ok("标准 Voigt" in str(e) or "重排" in str(e),
           "VASP 顺序未重排被护栏拦下：%s" % str(e)[:60])

# ---------- POP ----------
def pop_obj(cA, eps_inf_ip, eps_0_ip, qtf=0.0, lo=True):
    P._C = cA * P.A2B
    P._LO_DISP = lo
    p = P.PolarOptical2D.__new__(P.PolarOptical2D)
    P._REC["eps_inf_slab"] = np.diag([eps_inf_ip] * 2 + [1]).tolist()
    P._REC["eps_static_slab"] = np.diag([eps_0_ip] * 2 + [1]).tolist()
    p.r_inf = P._r_tensor("eps_inf_slab")
    p.r_0 = P._r_tensor("eps_static_slab")
    p.n_po = np.ones((1, 1))
    p.q_tf = np.full((1, 1), qtf)
    return p


f0 = np.zeros((1, 1))
u, nsq = rand_q()
p20 = pop_obj(20, 1 + 0.9, 1 + 1.6)
F20 = p20.factor(u, nsq, False, f0)
p30 = pop_obj(30, 1 + 0.9 * 20 / 30, 1 + 1.6 * 20 / 30)
F30 = p30.factor(u, nsq, False, f0)
ok(np.allclose(F30 / F20, 1.5),
   "POP 真空不变性：r 不变时 F(c=30)/F(c=20) = 1.5（正好抵消 AMSET 积分里的 1/c）")

p = pop_obj(20, 1.9, 2.6)
tiny = np.array([[1, 0, 0.0]])
Fq0 = p.factor(tiny, np.array([1e-14]), False, f0)[0, 0, 0]
dr = 20 * P.A2B * (2.6 - 1.9) / 2
ok(np.isclose(Fq0, 2 * np.pi * 20 * P.A2B * dr),
   "POP q→0 有限且 = 2πc·Δr/ε_env²（%.4g）" % Fq0)

ps = pop_obj(20, 1.9, 2.6, qtf=0.05)
r_ = ps.factor(tiny, np.array([1e-14]), False, f0)[0, 0, 0] / Fq0
ok(r_ < 1e-6, "POP 自由载流子屏蔽：q→0 的散射被压到无屏蔽值的 %.1e" % r_)

P._REC["eps_inf_slab"] = [[2.2, 0.3, 0], [0.3, 1.5, 0], [0, 0, 1]]
P._REC["eps_static_slab"] = [[3.1, 0.5, 0], [0.5, 2.0, 0], [0, 0, 1]]
pa = pop_obj(20, 1, 1)
pa.r_inf = P._r_tensor("eps_inf_slab")
pa.r_0 = P._r_tensor("eps_static_slab")
G1 = pa.factor(u, nsq, False, f0)
R2 = Rz(1.1)[:2, :2]
pa.r_inf = R2 @ pa.r_inf @ R2.T
pa.r_0 = R2 @ pa.r_0 @ R2.T
G2 = pa.factor(u @ Rz(1.1).T, nsq, False, f0)
ok(np.allclose(G1, G2), "POP 各向异性 r 张量：旋转不变")

# ---------- IMP ----------
P._C = 20 * P.A2B
P._D_IMP = 0.0
im = P.IonizedImpurity2D.__new__(P.IonizedImpurity2D)
im.r_0 = np.eye(2) * 3.0
im.q_tf = np.full((1, 1), 0.02)
u1 = np.array([[1.0, 0, 0]])
qq = 0.05
Fi = im.factor(u1, np.array([qq ** 2]), None, 0, None, None)[0, 0, 0]
ref = P._C ** 2 * (2 * np.pi) ** 2 / ((qq * (1 + 3.0 * qq)) + 0.02) ** 2
ok(np.isclose(Fi, ref), "IMP 数值与手算公式 c²(2π)²/(q·ε0(q)+q_TF)² 一致")
Fq0_imp = im.factor(u1, np.array([1e-16]), None, 0, None, None)[0, 0, 0]
ok(np.isfinite(Fq0_imp), "IMP 在 q→0 不发散（有屏蔽，值 %.3g）" % Fq0_imp)

# ---------- PIE ----------
e = np.zeros((3, 3, 3))
e11 = 0.3
e[0, 0, 0] = e11
e[0, 1, 1] = -e11
e[1, 0, 1] = e[1, 1, 0] = -e11
pz = P.Piezoelectric2D.__new__(P.Piezoelectric2D)
pz.e_raw = e
pz.elastic_constant = Cf
pz.r_inf = np.eye(2) * 2.0
pz.q_tf = np.zeros((1, 1))
Fp0 = pz.factor(np.array([[1.0, 0, 0]]), np.array([1e-14]), None, 0, None, None)[0, 0, 0]
Fp1 = pz.factor(np.array([[1.0, 0, 0]]), np.array([1e-4]), None, 0, None, None)[0, 0, 0]
ok(np.isfinite(Fp0) and Fp0 > 0 and abs(Fp1 / Fp0 - 1) < 0.1,
   "PIE 无屏蔽时 q→0 有限（%.3g），小 q 近似常数" % Fp0)
Fa = pz.factor(u, nsq, None, 0, None, None)
pz.elastic_constant = rot4(Cf, Rz(np.pi / 3))
pz.e_raw = np.einsum("ia,jb,kc,abc->ijk", Rz(np.pi / 3), Rz(np.pi / 3), Rz(np.pi / 3), e)
Fb = pz.factor(u @ Rz(np.pi / 3).T, nsq, None, 0, None, None)
ok(np.allclose(Fa, Fb), "PIE 张量：旋转不变")

# ---------- gen 脚本的 Voigt 重排 ----------
stub = types.ModuleType("stepconf")
stub.load = lambda *a, **k: {}
stub.apply_submit = lambda *a, **k: None
stub.read_submit = lambda *a, **k: {}
stub.CONF_NAME = "step.conf"
sys.modules.setdefault("stepconf", stub)
gen_path = HERE / "gen_step14_amset2d.py"
if gen_path.is_file():
    spec = importlib.util.spec_from_file_location("gen_step14", gen_path)
    g = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(g)
    except Exception as e:                                        # noqa: BLE001
        g = None
        print("[WARN] 跳过 gen 重排自检（导入失败：%s）" % e)
    if g is not None:
        m = np.arange(36, dtype=float).reshape(6, 6)
        out = np.array(g.vasp_to_voigt_6x6(m))
        ok(np.allclose(out[np.ix_([0, 1, 2], [0, 1, 2])], m[np.ix_([0, 1, 2], [0, 1, 2])])
           and out[5, 5] == m[3, 3] and out[3, 3] == m[4, 4] and out[4, 4] == m[5, 5],
           "Voigt 重排：标准位 5(XY)=VASP 位 3(XY)、位 3(YZ)=VASP 位 4(YZ)、位 4(XZ)=VASP 位 5(ZX)")
        ok(g.vasp_to_voigt_piezo([[1, 2, 3, 4, 5, 6]])[0] == [1, 2, 3, 5, 6, 4],
           "压电张量列序重排：VASP 列 (XY YZ ZX) -> 标准 Voigt 列 (YZ XZ XY)")
        ent = np.arange(6, dtype=float)
        ok(np.allclose(g._reorder_voigt(np.diag(ent).tolist(),
                                        ["XX", "YY", "ZZ", "XY", "YZ", "ZX"]),
                       np.diag([0, 1, 2, 4, 5, 3])),
           "_reorder_voigt：按来源标签把对角元放到标准 Voigt 位置")

print()
if FAILS:
    print("%d 项 FAIL：" % len(FAILS))
    for f in FAILS:
        print("  -", f)
    sys.exit(1)
print("全部 PASS（示例记录：%s）" % rec_src.name)
