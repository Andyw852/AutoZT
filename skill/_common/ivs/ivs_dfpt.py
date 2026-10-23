# -*- coding: utf-8 -*-
"""ivs_dfpt.py —— QE DFPT 直接算谷间零阶形变势（IVS）的纯 numpy 部分（[0059]）。

生产端：eph-qe-cpu 可选组 ivs_dfpt（四步，方案 docs/PLAN-ivs-dfpt-20261018.md 里的 Q1–Q4）
  S10_ivsq   （Q1，gen_step10_ivs_qlist.py）：谷表 -> 谷对 -> 按晶体对称性分等价类 -> 按 q 分组（ivs_qlist.json）；
  S11_ivsscf （Q2，gen_step11_ivs_scf.py）  ：pw.x SCF（VASP 结构原样，不重新弛豫）+ 各谷能差自检；
  S12_ivsph  （Q3，gen_step12_ivs_ph.py）   ：每个 q 一个子目录：ph.x 单 q DFPT -> pw.x nscf（初态 k）-> ph.x AHC 矩阵元；
  S13_ivs    （Q4，gen_step13_ivs_read.py） ：ahc_gkk + 动力学矩阵 -> D0、ħω -> dp_optical 的星通道（ivs_dfpt_params.json）。
消费端：ke-dft-cpu 的 S8（3D）/ S8.4（2D）的 IVS_DFPT 开关：inter_channels / params_digest / describe / valley_sig。
本模块住在 skill/_common/ivs/，两个技能 gen_need 都写 ivs_dfpt.py（公共池按文件名找到它），只有一份、不会漂移。

QE 接口（M0 核实，QE 7.3 与 7.4 的 PHonon/PH/ahc.f90 输出部分逐行相同；记录见 METHODOLOGY [0059]）：
  ph.x 的 electron_phonon='ahc'（要 trans=.false.，读前一次 trans=.true. 存下的 dvscf 与位移模式）在 pw.x nscf 给的
  每个 k 上写 ahc_dir/ahc_gkk_iq1.bin：直接存取、无记录头，每个 k 一条记录，Fortran 次序 ahc_gkk(ib, jb, imode)：
      ahc_gkk(ib, jb, imode) = ⟨ψ_ib(k+q)| ∂V_SCF/∂u_{q,imode} |ψ_jb(k)⟩       单位 Ry/bohr
  ib = 1..nbnd（k+q 的全部带）、jb = ahc_nbndskip+1..ahc_nbndskip+ahc_nbnd（k 上的带）、imode = 3(原子−1)+方向（笛卡尔，
  一个原子沿一个方向按 e^{iq·R} 位移）。ahc_etk/etq_iq1.bin 是 k、k+q 上全部带的本征值（Ry）。
  setup_nscf 在 AHC 模式下跳过 irreducible_BZ —— nscf 给的 k 原样保留（M0：Si 的 k、k+q 与输入一致）。
  只支持模守恒赝势（ahc.f90：USPP/PAW/DFPT+U/磁性/2D 库仑截断未实现或未测试）。
  核对：q = Γ 时 ⟨n|∂V/∂u_{2x}|n⟩ 与冻声子有限差分 dε/du 四位有效数字一致（Si，0.00963/0.05539/0.24525/−0.20993 vs
  0.00962/0.05537/0.24520/−0.20987 eV/Å）；动力学矩阵按 fildyn 文件的 (na, nb) 块、行 = na 的方向解析，本征矢与 QE
  dynmat.x 重叠 = 1（取共轭则否，无反演中心的结构上验证）—— 与 postahc.x 的投影（g_ν = Σ_j g_j u_jν，u = e/√M）同一口径。

形变势（与 dp_optical / ivs_frozen 同一口径，|g|² = ħ D²/(2 N M_cell ω)）：
  g_ν(m,n) = Σ_j ahc_gkk(m, n, j) · e_jν / √M_j        （e 为质量加权、归一的本征矢，M 用 QE 的 amass，Ry 质量单位）
  D_ν² = M_cell · Σ_{n∈初态简并组} Σ_{m∈终态简并组} |g_ν(m,n)|² / (d_初 · d_终)      （规范不变；Ry/bohr -> eV/Å）
  声子简并组（ħω 相差 < PH_DEG_MEV）按 D² 相加当一支。
星通道（dp_optical 的态密度按星平均）：初态取星 X 的第一个成员 i0（等价成员的率相同），终态取所有星的所有成员 j：
  X>Y（别的星）：每个（成员 j、声子组）一条，D_ch² = D_{i0 j}² / g_Y；
  X>X'（同星其它成员）：D_ch² = D_{i0 j}² / (g_X − 1)；同一对星、ħω 相差 < HW_MERGE_MEV 的按 D² 相加（dp_optical 拒绝重复 ħω）。
对称性约化：谷对 (k_a, k_b) 与 (S k_a, S k_b)、(−k_a, −k_b)（时间反演）、(k_b, k_a)（反向过程，厄米性）给出同一 |g|，
  同一等价类只算一次；共用同一个 q 的类合进同一个 ph.x 作业（nscf 里放几个初态 k）。

不做：SOC、一阶谷间耦合（∝|q−K|）、EPW 插值；谷内（q ≈ 0）的 ADP/ODP 不在这里（小 q 要 2D 库仑截断，AHC 不支持）。
[0074] SOC 只备注不拦：电子、空穴都照算（非 SOC），soc_notes 给每个载流子写口径说明（无反演中心时时间反演对谷
  k↔−k 的通道在 SOC 下要自旋翻转 -> 上限），进 ivs_dfpt_params.json / 表 / S8、S8.4 日志与 S8.3 报告。
"""
import hashlib
import itertools
import json
import math
import re

import numpy as np

SCHEMA = "ivs_dfpt/1"
QLIST_SCHEMA = "ivs_qlist/1"

RY_EV = 13.605693122994            # eV / Ry
BOHR_A = 0.529177210903            # Å / bohr
GKK_TO_EV_A = RY_EV / BOHR_A       # Ry/bohr -> eV/Å
AMU_RY = 911.44424310865645        # 1 amu 在 Ry 质量单位（m_e/2 ... QE AMU_RY）
RY_MEV = RY_EV * 1000.0

DEG_TOL_EV = 0.001                 # 电子能级简并判据（QE ahc 的 e_degen_thr = 1e-4 Ry ≈ 1.4 meV）
PH_DEG_MEV = 0.05                  # 声子简并判据
HW_MERGE_MEV = 0.1                 # 同一对星的两条通道 ħω 相差小于它就合并
HW_MIN_MEV = 0.5                   # ħω 低于它（或虚频）的模不进通道（失稳/数值零），记进 dropped
DMIN_EV_A = 0.05                   # S8/S8.4 用：D 小于它的通道不进散射率（对称性禁戒）
POLAR_TOL_E = 0.1                  # [0059-r3] 模有效电荷（e）大于它 = 极性模（= ke-dft-cpu odp_frozen.POLAR_TOL_E）
N_ACOUSTIC = 3                     # [0059-r3] q=Γ 剔除平移权重最大的 3 个模（声学模，ħω≈0、耦合由 ADP 管）
SYMPREC_A = 0.01                   # 找对称操作的原子位置容差（Å，与 AMSET 的 symprec 同量级）
KEY_ND = 5                         # k 点比较键的小数位（模 1）
VALLEY_TOL_MEV = 30.0              # Q2：QE 的谷能差与 S8.2 的差多少以内算一致
CARRIERS = ("electron", "hole")


# ---------------------------------------------------------------------------
# 谷表（与 dp_optical 同一套写法与归一；这里复制最小必要部分，免得 eph-qe-cpu 依赖 ke-dft-cpu 的文件）
# ---------------------------------------------------------------------------
_ALIAS = {"G": "Γ", "GAMMA": "Γ", "Γ": "Γ", "GM": "Γ–M", "Γ-M": "Γ–M", "Γ–M": "Γ–M",
          "G-M": "Γ–M", "G–M": "Γ–M", "GAMMA-M": "Γ–M"}


def norm_label(s):
    """谷标签归一（与 dp_optical.norm_label 逐字相同）：Γ 的各种写法统一成 Γ；其余去空白、ASCII 字母转大写。"""
    t = str(s).strip()
    u = t.upper()
    if u in _ALIAS:
        return _ALIAS[u]
    return u if re.fullmatch(r"[A-Za-z0-9_()'.,+-]+", t) else t


def parse_valleys(spec, full3d=False):
    """'X:0,0.425,0.425|...; L:...' -> [{"label","g","k_frac"}]（与 dp_optical.parse_valleys 同一写法）。"""
    out = []
    for item in [p.strip() for p in str(spec or "").replace("\n", ";").split(";")]:
        if not item:
            continue
        if ":" not in item:
            raise ValueError("谷 %r 写法不对，应为 标签:kx,ky,kz|kx,ky,kz（分数坐标）" % item)
        lab, rest = item.split(":", 1)
        ks = []
        for kk in rest.split("|"):
            v = [float(x) for x in kk.replace(" ", "").split(",") if x != ""]
            if len(v) == 2 and full3d:
                raise ValueError("谷 %r 的 k 点 %r 只有 2 个分量：三维谷中心必须写 kx,ky,kz" % (item, kk))
            if len(v) == 2:
                v.append(0.0)
            if len(v) != 3:
                raise ValueError("谷 %r 的 k 点 %r 不是 2 或 3 个分量" % (item, kk))
            ks.append(v)
        out.append({"label": norm_label(lab), "g": len(ks), "k_frac": ks})
    return out


def valleys_from_dpt(dpt, carrier):
    """S8.2 dpt_result.json 某载流子的谷表（multi_valley.valleys），与 dp_optical.valleys_from_dpt 同口径，外加 raw 原表。"""
    for r in (dpt or {}).get("results", []):
        if r.get("carrier") != carrier:
            continue
        raw = ((r.get("multi_valley") or {}).get("valleys")) or []
        out = []
        for v in raw:
            ks = v.get("k_frac") or []
            if not ks:
                continue
            row = {"label": norm_label(v.get("label")), "g": len(ks),
                   "k_frac": [[float(x) for x in k] for k in ks], "dE_eV": v.get("dE_eV")}
            kg = v.get("k_grid") or []
            if len(kg) == len(ks):
                row["k_grid"] = [[float(x) for x in k] for k in kg]
            out.append(row)
        return out, raw, (r.get("multi_valley") or {}).get("status")
    return [], [], None


def grid_key(k, full3d=False):
    """分数坐标比较键：KEY_ND 位、模 1（与 dp_optical.grid_key 同口径；full3d 时三个分量都参与）。"""
    out = []
    for x in list(k)[:3 if full3d else 2]:
        y = round(round(float(x), KEY_ND) % 1.0, KEY_ND)
        out.append(0.0 if y >= 1.0 else y + 0.0)
    return tuple(out)


def _member_match(v, u, full3d=False):
    if v.get("k_grid") and u.get("k_grid"):
        uk = {grid_key(w, full3d) for w in u["k_grid"]}
        return [grid_key(k, full3d) in uk for k in v["k_grid"]]

    def same(a, b):
        d = np.asarray(a, float) - np.asarray(b, float)
        d = d - np.rint(d)
        if not full3d:
            d[2:] = 0.0
        return bool(np.all(np.abs(d) < 1e-4))
    return [any(same(k, w) for w in u["k_frac"]) for k in (v.get("k_frac") or [])]


def merge_coincident(valleys, full3d=False):
    """同一组 k 中心的重复星（S8.2 按能带/自旋分出来的）并成一个，与 dp_optical.merge_coincident 同一规则。
    返回 (新谷表, {被并掉的标签: 保留的标签})；部分重合 -> ValueError。"""
    out, alias = [], {}
    for v in valleys:
        ks = v.get("k_frac") or []
        hit = None
        for u in out:
            same = _member_match(v, u, full3d)
            if same and all(same) and len(ks) == len(u["k_frac"]):
                hit = u
                break
            if any(same):
                raise ValueError("谷表里 %s 与 %s 只有部分中心重合（%d/%d）——谷表自相矛盾"
                                 % (v.get("label"), u["label"], sum(same), len(ks)))
        if hit is None:
            out.append(dict(v))
        else:
            alias[v["label"]] = hit["label"]
    return out, alias


def valley_sig(valleys, use_grid=False, full3d=False):
    """谷表签名：每个星的中心集合（模 1），与标签、星顺序、重复星无关。2D 与 ke_common._valley_sig 相同
    （面内两个分量；use_grid 用 k_grid、5 位，否则 k_frac、4 位）；full3d 取三个分量。"""
    key, nd = ("k_grid", 5) if use_grid else ("k_frac", 4)
    stars = set()
    for v in valleys or []:
        ks = []
        for k in v.get(key) or []:
            f = []
            for x in list(k)[:3 if full3d else 2]:
                y = round(round(float(x), nd) % 1.0, nd)
                f.append(0.0 if y >= 1.0 else y + 0.0)
            ks.append(tuple(f))
        stars.add(tuple(sorted(ks)))
    return [list(map(list, s)) for s in sorted(stars)]


# ---------------------------------------------------------------------------
# 结构
# ---------------------------------------------------------------------------
def read_poscar_text(text):
    """VASP 5 POSCAR 文本 -> {"lattice" (3,3) Å 行矢量, "symbols" 每原子元素, "types" 元素顺序, "counts", "frac" (n,3)}。
    支持缩放因子（含负的体积写法）、Selective dynamics、Direct/Cartesian。"""
    ls = [x.strip() for x in str(text).splitlines() if x.strip()]
    if len(ls) < 8:
        raise ValueError("POSCAR 不完整")
    scale = float(ls[1].split()[0])
    raw = np.array([[float(x) for x in ls[i].split()[:3]] for i in range(2, 5)])
    types = ls[5].split()
    try:
        counts = [int(x) for x in ls[6].split()]
    except ValueError:
        raise ValueError("只支持 VASP 5 格式 POSCAR（第 6 行元素名、第 7 行原子数）")
    if len(types) != len(counts):
        raise ValueError("POSCAR 元素数与原子数行不匹配")
    det = float(np.linalg.det(raw))
    if abs(det) < 1e-12:
        raise ValueError("POSCAR 晶胞奇异")
    s = scale if scale > 0 else (abs(scale) / abs(det)) ** (1.0 / 3.0)
    lat = raw * s
    pos = 7
    if ls[pos].lower().startswith("s"):
        pos += 1
    direct = ls[pos].lower().startswith("d")
    if not ls[pos].lower().startswith(("d", "c", "k")):
        raise ValueError("POSCAR 缺 Direct/Cartesian 标记")
    pos += 1
    n = sum(counts)
    xyz = np.array([[float(x) for x in ls[i].split()[:3]] for i in range(pos, pos + n)])
    frac = xyz if direct else (xyz * s) @ np.linalg.inv(lat)
    symbols = [t for t, c in zip(types, counts) for _ in range(c)]
    return {"lattice": lat, "symbols": symbols, "types": types, "counts": counts, "frac": frac}


def detect_dim(lattice, frac, vacuum_A=6.0):
    """第三个晶格矢方向（c）上原子之间的最大空隙 > vacuum_A（Å）且 c 与面内基矢垂直 -> "2d"，否则 "3d"。"""
    lat = np.asarray(lattice, float)
    z = np.sort(np.mod(np.asarray(frac, float)[:, 2], 1.0))
    gaps = np.diff(np.concatenate([z, [z[0] + 1.0]]))
    c = lat[2]
    n = np.cross(lat[0], lat[1])
    n /= np.linalg.norm(n)
    height = abs(float(c @ n))
    perp = abs(float(c @ n)) / float(np.linalg.norm(c)) > 0.999
    return "2d" if perp and float(gaps.max()) * height > vacuum_A else "3d"


def recip_lattice(lattice):
    """倒格矢（行 = b1,b2,b3，含 2π，单位 1/Å）。"""
    return 2.0 * math.pi * np.linalg.inv(np.asarray(lattice, float)).T


def auto_kgrid(lattice, kspacing, dim):
    """VASP KSPACING 口径：N_i = max(1, ceil(|b_i| / kspacing))（|b| 含 2π，Å⁻¹）；2D 的 c 方向固定 1。"""
    b = np.linalg.norm(recip_lattice(lattice), axis=1)
    n = [max(1, int(math.ceil(float(x) / float(kspacing) - 1e-9))) for x in b]
    if dim == "2d":
        n[2] = 1
    return n


# ---------------------------------------------------------------------------
# 对称性（不依赖 spglib：晶格度规不变的 {−1,0,1} 整数矩阵 + 原子映射；约化过的胞足够，漏操作只会少约化、不会出错）
# ---------------------------------------------------------------------------
def lattice_rotations(lattice, tol=1e-3):
    """x' = W·x（直接晶格分数坐标）保持度规 G = A·Aᵀ 的整数矩阵 W（元素 −1/0/1，det = ±1）。"""
    A = np.asarray(lattice, float)
    G = A @ A.T
    cand = np.array(list(itertools.product((-1, 0, 1), repeat=9)), dtype=int).reshape(-1, 3, 3)
    det = np.rint(np.linalg.det(cand.astype(float))).astype(int)
    cand = cand[np.abs(det) == 1]
    gw = np.einsum("nji,jk,nkl->nil", cand, G, cand)
    ok = np.all(np.abs(gw - G[None]) <= tol * float(np.abs(G).max()), axis=(1, 2))
    return cand[ok]


def point_group(lattice, frac, symbols, symprec=SYMPREC_A):
    """结构的点群（直接空间分数坐标的整数旋转 W 列表；带分数平移的操作也算，k 空间只看旋转部分）。"""
    A = np.asarray(lattice, float)
    f = np.mod(np.asarray(frac, float), 1.0)
    sp = np.array([str(s) for s in symbols])
    same = sp[:, None] == sp[None, :]
    ops = []
    for W in lattice_rotations(A):
        fp = f @ W.T
        for j in np.where(sp == sp[0])[0]:
            t = f[j] - fp[0]
            d = (fp + t)[:, None, :] - f[None, :, :]
            d = d - np.rint(d)
            dist = np.linalg.norm(d @ A, axis=2)
            if np.all(np.min(np.where(same, dist, np.inf), axis=1) < symprec):
                ops.append(W.copy())
                break
    return ops


def recip_ops(rotations):
    """直接空间旋转 W -> 倒空间分数坐标上的作用 κ' = W⁻ᵀ κ（整数矩阵）。"""
    out = []
    for W in rotations:
        R = np.rint(np.linalg.inv(np.asarray(W, float)).T).astype(int)
        out.append(R)
    return out


def _k_images(k, rops):
    """k 在点群 + 时间反演下的全部像（分数坐标，未取模）。"""
    k = np.asarray(k, float)
    out = []
    for R in rops:
        kk = R @ k
        out.append(kk)
        out.append(-kk)
    return out


def expand_star(ks, rops, full3d=False):
    """星的成员 -> 对称轨道（点群 + 时间反演）。返回 (成员列表, 轨道数, 是否补了成员)。
    成员按给出的顺序在前，补出来的接在后面；轨道数 > 1 说明给的成员不是同一个星（调用方报错）。"""
    # 给出的成员原样在前（保留用户/S8.2 的坐标写法），补出来的取模到 (−0.5, 0.5]
    out, seen, roots = [], set(), []
    for k in [np.asarray(k, float) for k in ks]:
        if grid_key(k, full3d) not in seen:
            seen.add(grid_key(k, full3d))
            out.append(k)
    n_given = len(out)
    for k in list(out[:n_given]):
        imgs = [k] + _k_images(k, rops)
        ikeys = {grid_key(im, full3d) for im in imgs}
        if not any(r in ikeys for r in roots):          # 与已有的轨道代表都不等价 -> 新轨道
            roots.append(grid_key(k, full3d))
        for im in imgs:
            if grid_key(im, full3d) not in seen:
                seen.add(grid_key(im, full3d))
                out.append(wrap_half(im))
    return out, len(roots), len(out) > n_given


def wrap_half(v):
    """分数坐标取到 (−0.5, 0.5]。"""
    v = np.asarray(v, float)
    w = v - np.rint(v)
    w[np.isclose(w, -0.5, atol=1e-9)] = 0.5
    return w + 0.0


def pair_key(a, b, rops, full3d=True):
    """谷对 (a→b) 的等价类键：点群 S、时间反演（−）、反向（b→a）下所有像的最小键。"""
    best = None
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    for R in rops:
        for s in (1.0, -1.0):
            ra, rb = s * (R @ a), s * (R @ b)
            for x, y in ((ra, rb), (rb, ra)):
                key = (grid_key(x, full3d), grid_key(y, full3d))
                if best is None or key < best:
                    best = key
    return best


def pair_images(a, b, rops):
    """谷对的全部等价实现 [(k_初, k_终)]（含时间反演与反向），去重。"""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    out, seen = [], set()
    for R in rops:
        for s in (1.0, -1.0):
            ra, rb = s * (R @ a), s * (R @ b)
            for x, y in ((ra, rb), (rb, ra)):
                key = (grid_key(x, True), grid_key(y, True))
                if key not in seen:
                    seen.add(key)
                    out.append((x, y))
    return out


def same_q(q1, q2, tol=1e-4):
    d = np.asarray(q1, float) - np.asarray(q2, float)
    d = d - np.rint(d)
    return bool(np.all(np.abs(d) < tol))


# ---------------------------------------------------------------------------
# Q1：谷对 -> 等价类 -> q 作业
# ---------------------------------------------------------------------------
def plan_pairs(carriers, rops, full3d=True):
    """carriers = {载流子: [{"label","k": [k...]}, ...]}（星，成员已是完整轨道）。
    返回 (classes, runs)：
      classes：[{"id","carrier","from","to","kind","k_init","k_fin","q","run","ik","n_pairs"}]，每个等价类一条；
      runs：[{"id","q","k": [k...]}]，每个 q 一个 ph.x 作业，k = 该作业 nscf 里的初态 k（顺序即 ahc_gkk 的记录号）；
      另外每条 pair 记录 [{"carrier","from","to","kind","j","cls"}] 挂在 classes 上（"pairs"）。
    初态取每个星的第一个成员 i0；终态是全部星的全部成员（同星时跳过 i0 自己）。"""
    classes, index = [], {}
    for car in CARRIERS:
        stars = carriers.get(car) or []
        for X in stars:
            a = np.asarray(X["k"][0], float)
            for Y in stars:
                for j, b in enumerate(Y["k"]):
                    if Y is X and j == 0:
                        continue
                    b = np.asarray(b, float)
                    key = (car, pair_key(a, b, rops, full3d=True))
                    kind = "inter_star" if Y is X else "inter"
                    if key not in index:
                        index[key] = len(classes)
                        classes.append({"id": len(classes), "carrier": car, "from": X["label"], "to": Y["label"],
                                        "kind": kind, "a": a, "b": b, "pairs": []})
                    classes[index[key]]["pairs"].append({"from": X["label"], "to": Y["label"], "kind": kind,
                                                         "j": j, "g_to": len(Y["k"]), "g_from": len(X["k"])})
    runs = []
    for c in classes:
        placed = False
        for x, y in pair_images(c["a"], c["b"], rops):
            for r in runs:
                if same_q(y - x, r["q"]):
                    ik = next((i for i, kk in enumerate(r["k"]) if same_q(kk, x, 1e-6)), None)
                    if ik is None:
                        r["k"].append(wrap_half(x))
                        ik = len(r["k"]) - 1
                    c.update(run=r["id"], ik=ik, k_init=wrap_half(x), k_fin=wrap_half(y))
                    placed = True
                    break
            if placed:
                break
        if not placed:
            q = wrap_half(c["b"] - c["a"])
            runs.append({"id": len(runs), "q": q, "k": [wrap_half(c["a"])]})
            c.update(run=runs[-1]["id"], ik=0, k_init=wrap_half(c["a"]), k_fin=wrap_half(c["b"]))
    for c in classes:
        c["q"] = wrap_half(np.asarray(c["k_fin"]) - np.asarray(c["k_init"]))
        c["n_pairs"] = len(c["pairs"])
    return classes, runs


def jsonable(x):
    """numpy -> 纯 Python（json.dump 用）。"""
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, float):
        # [0059-r1] 保留 12 位有效数字（原来 round(x, 10) 按小数位舍入，把 tr2_ph = 1e-16 写成了 0.0）
        return float("%.12g" % x) + 0.0
    return x


# ---------------------------------------------------------------------------
# QE 输入
# ---------------------------------------------------------------------------
def xq_alat(q_frac, cell_bohr, alat_bohr):
    """分数坐标 q -> ph.x 的 xq（笛卡尔，2π/alat 单位）。cell_bohr：行矢量 a1..a3（bohr，取自 QE 的 XML）。"""
    A = np.asarray(cell_bohr, float)
    B = np.linalg.inv(A).T * float(alat_bohr)
    return np.asarray(q_frac, float) @ B


def qe_structure(lattice, frac, symbols, types, masses, pseudos):
    """CELL_PARAMETERS（Å）+ ATOMIC_SPECIES + ATOMIC_POSITIONS（crystal），与 VASP 结构逐位相同（不重新弛豫）。"""
    lines = ["CELL_PARAMETERS angstrom"]
    lines += ["  %.10f %.10f %.10f" % tuple(v) for v in np.asarray(lattice, float)]
    lines += ["", "ATOMIC_SPECIES"]
    lines += ["%s %.8g %s" % (t, float(masses[t]), pseudos[t]) for t in types]
    lines += ["ATOMIC_POSITIONS crystal"]
    lines += ["%s %.10f %.10f %.10f" % ((s,) + tuple(float(x) for x in f)) for s, f in zip(symbols, frac)]
    return "\n".join(lines) + "\n"


def kpoints_list(ks):
    """显式 k 列表（crystal，权重 1）。"""
    ks = [np.asarray(k, float) for k in ks]
    return "K_POINTS crystal\n%d\n" % len(ks) + "".join("  %.10f %.10f %.10f 1.0\n" % tuple(k) for k in ks)


# ---------------------------------------------------------------------------
# Q4：读 QE 产物
# ---------------------------------------------------------------------------
def read_ahc(dirpath, nbnd, ahc_nbnd, nat, nk):
    """ahc_dir/ahc_gkk_iq1.bin + ahc_etk/etq_iq1.bin -> (gkk (nk, 3nat, ahc_nbnd, nbnd) 复数 Ry/bohr,
    etk (nk, nbnd) eV, etq (nk, nbnd) eV)。文件大小与维度对不上 -> ValueError（不猜）。"""
    import os
    nm = 3 * int(nat)
    fg = os.path.join(dirpath, "ahc_gkk_iq1.bin")
    want = int(nk) * nm * int(ahc_nbnd) * int(nbnd) * 16
    have = os.path.getsize(fg) if os.path.isfile(fg) else -1
    if have != want:
        raise ValueError("%s 大小 %d 字节，按 nk=%d、3nat=%d、ahc_nbnd=%d、nbnd=%d 应为 %d" %
                         (fg, have, nk, nm, ahc_nbnd, nbnd, want))
    g = np.fromfile(fg, dtype="<c16").reshape(int(nk), nm, int(ahc_nbnd), int(nbnd))
    out = [g]
    for name in ("ahc_etk_iq1.bin", "ahc_etq_iq1.bin"):
        fe = os.path.join(dirpath, name)
        wante = int(nk) * int(nbnd) * 8
        if not os.path.isfile(fe) or os.path.getsize(fe) != wante:
            raise ValueError("%s 缺失或大小不是 %d 字节" % (fe, wante))
        out.append(np.fromfile(fe, dtype="<f8").reshape(int(nk), int(nbnd)) * RY_EV)
    return tuple(out)


_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"


def read_dyn(path):
    """ph.x 的 fildyn 文本（非 XML）-> {"masses_ry": (nat,), "q_list": [(3,)…], "C": [(3nat,3nat)…], "nat", "alat"}。
    每个 "Dynamical  Matrix in cartesian axes" 块一个 q（ph.x 写 q 的整个星，第一个是输入的 q）；
    块内 "na nb" 后三行 = 行方向 i（na 的笛卡尔分量）、每行 3 个复数 = 列方向 j（nb 的分量）。"""
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    head = lines[2].split()
    ntyp, nat = int(head[0]), int(head[1])
    alat = float(head[3])
    i = 3
    if lines[i].strip().lower().startswith("basis"):
        i += 4
    amass = {}
    for t in range(ntyp):
        parts = lines[i + t].split("'")
        amass[int(parts[0])] = float(parts[2])
    i += ntyp
    ityp = [int(lines[i + a].split()[1]) for a in range(nat)]
    masses = np.array([amass[t] for t in ityp], float)
    qs, Cs = [], []
    n = i + nat
    while n < len(lines):
        if "Dynamical  Matrix in cartesian axes" in lines[n]:
            q = [float(x) for x in re.findall(_NUM, lines[n + 2])][:3]
            C = np.zeros((3 * nat, 3 * nat), complex)
            p = n + 4
            for _ in range(nat * nat):
                a, b = [int(x) for x in lines[p].split()[:2]]
                p += 1
                for r in range(3):
                    v = [float(x) for x in lines[p].split()]
                    p += 1
                    C[3 * (a - 1) + r, 3 * (b - 1):3 * (b - 1) + 3] = np.array(v[0:6:2]) + 1j * np.array(v[1:6:2])
            qs.append(np.array(q))
            Cs.append(C)
            n = p
        else:
            n += 1
    if not Cs:
        raise ValueError("%s 里没有 Dynamical Matrix 块" % path)
    eps, zstar = None, None                    # [0059-r3] q=Γ 且 epsil=.true.：介电张量与 Born 有效电荷
    for n, ln in enumerate(lines):
        if "Dielectric Tensor:" in ln:
            rows = [r for r in lines[n + 1:n + 6] if r.strip()][:3]
            eps = np.array([[float(x) for x in re.findall(_NUM, r)[:3]] for r in rows])
        elif "Effective Charges E-U" in ln:
            # QE write_epsilon_and_zeu：每个原子 "atom #  na" + 3 行，第 α 行 = 电场方向 α，3 列 = 位移方向 β
            z, p = [], n + 1
            while len(z) < nat and p < len(lines):
                if lines[p].strip().startswith("atom #"):
                    z.append([[float(x) for x in re.findall(_NUM, lines[p + 1 + r])[:3]] for r in range(3)])
                    p += 4
                else:
                    p += 1
            if len(z) == nat:
                zstar = np.array(z, float)
    return {"masses_ry": masses, "q_list": qs, "C": Cs, "nat": nat, "alat": alat, "ityp": ityp,
            "epsilon": eps, "zstar": zstar}


def phonon_modes(C, masses_ry):
    """C(q)（Ry/bohr²）-> (ħω meV（虚频记负）, e 质量加权归一本征矢 (3nat, 3nat)，列 = 模，按 ħω 升序)。"""
    mw = np.repeat(1.0 / np.sqrt(np.asarray(masses_ry, float)), 3)
    D = np.asarray(C) * np.outer(mw, mw)
    D = 0.5 * (D + D.conj().T)
    w2, e = np.linalg.eigh(D)
    hw = np.sign(w2) * np.sqrt(np.abs(w2)) * RY_MEV
    return hw, e


def degenerate_group(energies_ev, idx, tol=DEG_TOL_EV):
    """与 idx 号能级简并（|ΔE| < tol，连成链）的全部能级下标（升序）。"""
    e = np.asarray(energies_ev, float)
    grp = {int(idx)}
    changed = True
    while changed:
        changed = False
        for j in range(len(e)):
            if j not in grp and any(abs(e[j] - e[i]) < tol for i in grp):
                grp.add(j)
                changed = True
    return sorted(grp)


def phonon_groups(hw_mev, tol=PH_DEG_MEV):
    """ħω 升序的模 -> 简并组 [[模号…]…]。"""
    out = []
    for nu, w in enumerate(hw_mev):
        if out and abs(w - hw_mev[out[-1][-1]]) < tol:
            out[-1].append(nu)
        else:
            out.append([nu])
    return out


def acoustic_modes(e, masses_ry, n=N_ACOUSTIC):
    """[0059-r3] q=Γ：平移子空间投影权重最大的 n 个模（0 基下标）+ 全部模的权重。
    平移 T_α = √m_κ ê_α / √ΣM（与 ke-dft-cpu odp_frozen.mode_table 同口径）；e 是质量加权本征矢（列 = 模）。"""
    m = np.repeat(np.asarray(masses_ry, float), 3)
    nat = len(m) // 3
    T = np.zeros((3 * nat, 3))
    for a in range(3):
        T[a::3, a] = np.sqrt(m[a::3] / m[a::3].sum())
    tw = np.sum(np.abs(T.T @ np.asarray(e)) ** 2, axis=0)
    return sorted(int(i) for i in np.argsort(-tw)[:n]), tw


def mode_charges(e, masses_ry, zstar, two_d=False):
    """[0059-r3] Born 有效电荷 zstar (nat, 3 电场, 3 位移) -> 每个模的模有效电荷 |Σ Z*·U|/|U|（e），U = e/√M。
    2D 只看面内偶极（二维 POP 只含面内极性模；与 odp_frozen 的面内判据同口径），3D 看总偶极。"""
    m = np.repeat(np.asarray(masses_ry, float), 3)
    U = np.asarray(e) / np.sqrt(m)[:, None]                               # (3nat, 模)
    Z = np.asarray(zstar, float)
    nat = Z.shape[0]
    Zm = np.zeros((3, 3 * nat))
    for k in range(nat):
        Zm[:, 3 * k:3 * k + 3] = Z[k]
    dip = Zm @ U                                                           # (3, 模)
    if two_d:
        dip = dip[:2]
    return np.sqrt(np.sum(np.abs(dip) ** 2, axis=0)) / np.sqrt(np.sum(np.abs(U) ** 2, axis=0))


def mode_D2(gkk_k, e, masses_ry, init_ahc, final_bands):
    """一个 k 的 gkk (3nat, ahc_nbnd, nbnd) -> 每个模的 D²（(eV/Å)²），已对初/终态简并组求和并除以 d_初·d_终。
    init_ahc：初态在 ahc 窗口里的下标（0 基，相对 ahc_nbndskip）；final_bands：终态在 k+q 全部带里的下标（0 基）。"""
    m = np.asarray(masses_ry, float)
    u = np.asarray(e) / np.repeat(np.sqrt(m), 3)[:, None]                  # 位移 u = e/√M
    sub = np.asarray(gkk_k)[:, init_ahc, :][:, :, final_bands]           # (3nat, d_i, d_f)
    gm = np.einsum("jab,jn->nab", sub, u)                                 # (模, d_i, d_f)
    s = np.sum(np.abs(gm) ** 2, axis=(1, 2)) / (len(init_ahc) * len(final_bands))
    return s * float(m.sum()) * GKK_TO_EV_A ** 2


def edge_band(nocc, carrier):
    """带边带号（0 基）：电子 = nocc（最低导带），空穴 = nocc − 1（最高价带）。"""
    return int(nocc) if carrier == "electron" else int(nocc) - 1


def class_couplings(gkk_k, etk_k, etq_k, hw, e, masses_ry, nocc, ahc_skip, ahc_nbnd, carrier,
                    deg_tol=DEG_TOL_EV, ph_tol=PH_DEG_MEV, exclude=None):
    """一个等价类（一个 k、一个 q）-> {"init_bands","final_bands","E_init_eV","E_final_eV","modes":[…],"groups":[…]}。
    groups：声子简并组，{"modes","hw_meV","D_eV_A"}（D² 为组内相加）。初态简并组越出 ahc 窗口 -> ValueError。
    [0059-r3] exclude = {模下标（0 基）: 原因}（q=Γ 的声学模、极性模）：模照样记下，所在的组标 "exclude": 原因，
    fold_channels 不让它进通道。"""
    nb = edge_band(nocc, carrier)
    ib = degenerate_group(etk_k, nb, deg_tol)
    fb = degenerate_group(etq_k, nb, deg_tol)
    lo, hi = int(ahc_skip), int(ahc_skip) + int(ahc_nbnd)
    if ib[0] < lo or ib[-1] >= hi:
        raise ValueError("初态简并组（带 %s）越出 AHC 窗口 %d..%d：加大 IVS_AHC_NB" % ([x + 1 for x in ib], lo + 1, hi))
    d2 = mode_D2(gkk_k, e, masses_ry, [x - lo for x in ib], fb)
    hw = np.asarray(hw, float)
    exclude = exclude or {}
    modes = [{"mode": int(nu) + 1, "hw_meV": float(hw[nu]), "D_eV_A": float(math.sqrt(max(d2[nu], 0.0)))}
             for nu in range(len(hw))]
    for nu, why in exclude.items():
        modes[int(nu)]["exclude"] = why
    groups = []
    for grp in phonon_groups(list(hw), ph_tol):
        g = {"modes": [x + 1 for x in grp], "hw_meV": float(np.mean(hw[grp])),
             "D_eV_A": float(math.sqrt(max(float(np.sum(d2[grp])), 0.0)))}
        why = sorted({exclude[x] for x in grp if x in exclude})
        if why:
            g["exclude"] = "；".join(why)
        groups.append(g)
    return {"init_bands": [x + 1 for x in ib], "final_bands": [x + 1 for x in fb],
            "E_init_eV": float(etk_k[nb]), "E_final_eV": float(etq_k[nb]), "modes": modes, "groups": groups}


# ---------------------------------------------------------------------------
# [0072] 谷内声学形变势（ADP）：小 q 声学支矩阵元的斜率（Li, Graziosi, Neophytou, PRB 104, 195201 (2021)）
# ---------------------------------------------------------------------------
# 带边能量随应变的移动（S7 形变势）要芯势对齐，静水分量 Ξd 有绝对参考的不确定；长波声学声子的 DFPT 矩阵元
# M_ν(q) = ⟨ψ(k+q)|∂V_SCF/∂u_ν|ψ(k)⟩ 自洽地含宏观势响应，不需要对齐（与 EPW 同口径）。对非简并带边：
#   G_ν = √M_cell · Σ_j gkk_j u_jν （u = e/√m）= 每 Å 质心位移的耦合（eV/Å），p_ν = Σ_κ √m_κ e_κν / √M_cell = 模 ν 的质心位移；
#   长波极限 G_ν = Σ_i p_ν,i w_i，w = i|q|·D·q̂ × 公共相位（电子波函数规范）——内应变（光学支混入，∝ q）随质心位移线性，已含在 w 里。
#   每个方向 n̂ 一个 q = dq·n̂：解 3×3 得 w，去掉公共相位得实矢量 v_n = s_n |q| D n̂（s_n = ±1 未知），
#   6 个方向（x、y、z、xy、xz、yz 对角）定对称张量 D（6 个分量，s_n 穷举），只差一个整体符号（ADP 率只用 D²）。
ADP_DQ_A = 0.03                    # Å⁻¹（含 2π）：小 q 的模长；Si 约为 ΓX 的 2.6%，ħω_LA ≈ 2 meV
ADP_DIRS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1))
ADP_FIT_TOL = 0.05                 # 拟合相对残差、虚部比例超过它就告警（status 仍 ok，记进 checks）
ADP_RUN_PREFIX = "adp"


def adp_q_fracs(lattice, dq=ADP_DQ_A, dirs=ADP_DIRS):
    """笛卡尔方向（结构的笛卡尔系，QE 同系）-> [(单位方向 (3,), q 分数坐标 (3,))]：q_cart = dq·n̂，κ_i = a_i·q/(2π)。"""
    A = np.asarray(lattice, float)
    out = []
    for d in dirs:
        n = np.asarray(d, float)
        n = n / np.linalg.norm(n)
        out.append((n, A @ (float(dq) * n) / (2.0 * math.pi)))
    return out


def adp_polarizations(e, masses_ry, acoustic):
    """声学模（下标列表）的质心位移 p_ν = Σ_κ √m_κ e_κν / √M_cell（复数，(模, 3)）。"""
    m = np.asarray(masses_ry, float)
    e = np.asarray(e)
    sw = np.repeat(np.sqrt(m), 3)[:, None] * e / math.sqrt(float(m.sum()))
    nat = len(m)
    return np.array([[np.sum(sw[a::3, nu]) for a in range(3)] for nu in acoustic])


def adp_vector(gkk_k, e, masses_ry, init_ahc, final_band, acoustic):
    """一个 k、一个小 q：w（复数 (3,)，eV/Å，含公共相位）= 质心位移 1 Å 沿 x/y/z 时的耦合。
    init_ahc：初态在 ahc 窗口里的下标（0 基）；final_band：k+q 上终态的带号（0 基）。"""
    m = np.asarray(masses_ry, float)
    u = np.asarray(e) / np.repeat(np.sqrt(m), 3)[:, None]
    g = np.asarray(gkk_k)[:, int(init_ahc), int(final_band)]                # (3nat,)
    G = math.sqrt(float(m.sum())) * GKK_TO_EV_A * (g @ u[:, list(acoustic)])  # (3,) 每个声学模
    P = adp_polarizations(e, masses_ry, acoustic)                           # (3 模, 3 方向)
    if abs(np.linalg.det(P)) < 1e-6:
        raise ValueError("三个声学模的质心位移线性相关（小 q 的声学模没分出来：dq 太大或结构不稳定）")
    return np.linalg.solve(P, G)


def adp_real(w):
    """去公共相位：v = Re(w·e^{−iθ})，θ 使 Σ w_i² 为实；返回 (v, |Im|/|w|)。"""
    w = np.asarray(w, complex)
    th = 0.5 * np.angle(np.sum(w * w))
    z = w * np.exp(-1j * th)
    nrm = float(np.linalg.norm(w))
    return z.real, (float(np.linalg.norm(z.imag)) / nrm if nrm > 0 else 0.0)


_SYM6 = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))


def adp_fit(dirs, vs, qnorm):
    """v_n = s_n |q| D n̂ -> 对称 D（eV）。s_1 = +1，其余穷举。返回 (D (3,3), 相对残差, signs)。"""
    dirs = np.asarray(dirs, float)
    vs = np.asarray(vs, float) / float(qnorm)
    rows = []
    for n in dirs:
        for i in range(3):
            r = np.zeros(6)
            for c, (a, b) in enumerate(_SYM6):
                if a == i:
                    r[c] += n[b]
                if b == i and a != b:
                    r[c] += n[a]
            rows.append(r)
    X = np.array(rows)
    best = None
    for tail in itertools.product((1.0, -1.0), repeat=len(dirs) - 1):
        s = np.array((1.0,) + tail)
        y = (vs * s[:, None]).ravel()
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = float(np.linalg.norm(X @ c - y))
        if best is None or res < best[0] - 1e-12:
            best = (res, c, s)
    res, c, s = best
    D = np.zeros((3, 3))
    for k, (a, b) in enumerate(_SYM6):
        D[a, b] = D[b, a] = c[k]
    scale = float(np.linalg.norm(vs))
    return D, (res / scale if scale > 0 else 0.0), [int(x) for x in s]


def cart_rotation(lattice, W):
    """直接空间分数坐标旋转 W -> 笛卡尔旋转 R = Aᵀ W A⁻ᵀ（A 行 = 晶格矢）。"""
    A = np.asarray(lattice, float)
    return A.T @ np.asarray(W, float) @ np.linalg.inv(A.T)


def adp_member_tensors(D, k0, members, lattice, rots, full3d=True):
    """参考谷 k0 的张量 D -> 星里每个成员 k_m 的张量 R D Rᵀ（R：把 k0 送到 ±k_m 的点群操作；时间反演不改 D）。
    找不到操作 -> ValueError（谷表不是一个对称轨道）。"""
    out = []
    k0 = np.asarray(k0, float)
    for km in members:
        km = np.asarray(km, float)
        hit = None
        for W in rots:
            kk = np.rint(np.linalg.inv(np.asarray(W, float)).T) @ k0
            if any(grid_key(s * kk, full3d) == grid_key(km, full3d) for s in (1.0, -1.0)):
                hit = W
                break
        if hit is None:
            raise ValueError("成员 %s 不在 %s 的对称轨道上" % (list(km), list(k0)))
        R = cart_rotation(lattice, hit)
        out.append(R @ np.asarray(D, float) @ R.T)
    return out


def adp_summary(D, k0, lattice):
    """张量的可读量：主值、沿谷轴（Γ→k0 的笛卡尔方向）的 â·D·â、静水 tr/3。k0 = Γ 时谷轴为 None。"""
    D = np.asarray(D, float)
    kc = np.asarray(k0, float) @ recip_lattice(lattice)
    ax = None
    if np.linalg.norm(kc) > 1e-8:
        a = kc / np.linalg.norm(kc)
        ax = float(a @ D @ a)
    return {"principal_eV": sorted(float(x) for x in np.linalg.eigvalsh(D)), "axis_eV": ax,
            "hydrostatic_eV": float(np.trace(D) / 3.0)}


def fold_channels(classes_out, hw_min=HW_MIN_MEV, merge_mev=HW_MERGE_MEV):
    """等价类结果 -> dp_optical 星通道（一个载流子）。
    classes_out：[{"pairs": [{"from","to","kind","g_to","g_from"}…], "groups": [{"hw_meV","D_eV_A"}…]}…]。
    返回 (channels, dropped)：channels = [{"from","to","kind","D_eV_A","hw_meV","n_terms"}]；
    dropped = ħω < hw_min（含虚频）的组、[0059-r3] 标了 exclude 的组（q=Γ 声学/极性模）。
    星平均：X>Y 除 g_Y；X>X' 除 g_X − 1；X>X（谷内 ODP，[0059-r3]）不除（每个谷各自的 D）。"""
    acc, dropped = {}, []
    for c in classes_out:
        for gr in c["groups"]:
            hw, D = float(gr["hw_meV"]), float(gr["D_eV_A"])
            if gr.get("exclude"):                  # [0059-r3] q=Γ 的声学模 / 极性模
                dropped.append({"class": c.get("id"), "hw_meV": hw, "D_eV_A": D, "why": gr["exclude"]})
                continue
            if not (hw >= hw_min) or not math.isfinite(D):
                dropped.append({"class": c.get("id"), "hw_meV": hw, "D_eV_A": D, "why": "ħω < %.2g meV 或虚频" % hw_min})
                continue
            for p in c["pairs"]:
                if p["kind"] == "inter_star":
                    gdiv = p["g_from"] - 1
                elif p["kind"] == "intra":          # [0059-r3] ODP：谷内光学耦合，每个谷的 D 不除星成员数
                    gdiv = 1
                else:
                    gdiv = p["g_to"]  # inter：X→Y，除目标星成员 g_Y
                acc.setdefault((p["from"], p["to"], p["kind"]), []).append((hw, D * D / float(gdiv)))
    out = []
    for (a, b, kind), terms in sorted(acc.items()):
        terms.sort()
        clusters = []
        for hw, d2 in terms:
            if clusters and hw - clusters[-1]["hw_last"] < merge_mev:
                cl = clusters[-1]
                cl["d2"] += d2
                cl["w"] += hw * d2
                cl["n"] += 1
                cl["hw_last"] = hw
                cl["hws"].append(hw)
            else:
                clusters.append({"d2": d2, "w": hw * d2, "n": 1, "hw_last": hw, "hws": [hw]})
        for cl in clusters:
            hw = cl["w"] / cl["d2"] if cl["d2"] > 0 else float(np.mean(cl["hws"]))
            out.append({"from": a, "to": b, "kind": kind, "D_eV_A": math.sqrt(cl["d2"]),
                        "hw_meV": round(hw, 2), "n_terms": cl["n"]})
    return out, dropped


# ---------------------------------------------------------------------------
# [0074] 非 SOC 备注：只说明口径，不改参数（指纹 params_digest 只看通道，加备注不让 S8/S8.4 重排）
# ---------------------------------------------------------------------------
SOC_TR_NOTE = ("非 SOC 口径：本结构无反演中心，%s 连的是时间反演对谷（k′ ≡ −k）。有 SOC 时这两个谷的带边自旋相反"
               "（自旋-谷锁定），谷间散射要自旋翻转：该带边 SOC 劈裂 ≫ k_BT 时（如 MoS₂ 等 TMD 单层价带 K，约 0.15 eV）"
               "该通道被强烈抑制，这里的值是上限；劈裂 ≪ k_BT 时（如 MoS₂ 导带 K，约 3 meV）影响小。定量要 SOC DFPT（本技能不做）")
SOC_HOLE_NOTE = "非 SOC 口径：空穴按非 SOC 价带算，价带 SOC 劈裂（如 Γ 简并价带的劈裂带）及其自旋选择未计入，数值是非 SOC 近似"


def has_inversion(rops):
    """点群（倒空间整数矩阵）里有没有反演 −I。"""
    m = -np.eye(3, dtype=int)
    return any(np.array_equal(np.asarray(R, int), m) for R in rops)


def tr_pair(a, b, tol=1e-4):
    """b ≡ −a（模倒格矢）且 a 不是时间反演不变点（2a ≢ 0）：时间反演对谷。"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return same_q(a + b, np.zeros(3), tol) and not same_q(2.0 * a, np.zeros(3), tol)


def soc_notes(classes_out, carriers, inversion):
    """[0074] 每个载流子的非 SOC 备注 -> {载流子: [备注…]}（没有就不出现）。
    无反演中心 + 有连时间反演对谷的谷间通道（MoS₂ 的 K↔K′、Q↔Q′）-> SOC_TR_NOTE（列出通道）；
    空穴另外（没有上一条时）给一句 SOC_HOLE_NOTE。classes_out：gen_read 的逐类结果（carrier/from/to/kind/k_init/k_fin）。"""
    out = {}
    for car in carriers:
        notes = []
        if not inversion:
            trp = sorted({"%s>%s%s" % (c["from"], c["to"], "'" if c["kind"] == "inter_star" else "")
                          for c in classes_out if c["carrier"] == car and c["kind"] != "intra"
                          and tr_pair(c["k_init"], c["k_fin"])})
            if trp:
                notes.append(SOC_TR_NOTE % "、".join(trp))
        if car == "hole" and not notes:
            notes.append(SOC_HOLE_NOTE)
        if notes:
            out[car] = notes
    return out


def carrier_notes(params, carrier):
    """消费端：ivs_dfpt_params.json 里某载流子的非 SOC 备注（旧文件没有 -> []）。"""
    return list((_carrier(params, carrier) or {}).get("soc_notes") or [])


# ---------------------------------------------------------------------------
# 消费端（ke-dft-cpu S8 / S8.4）
# ---------------------------------------------------------------------------
def _carrier(params, carrier):
    c = ((params or {}).get("carriers") or {}).get(carrier) or {}
    return c if c.get("status") == "ok" else None


def has_intra(params, carrier):
    """[0059-r3] 这份参数有没有算该载流子的谷内 ODP（q=Γ，IVS_ODP 开着）。算了但全部 < dmin（对称性禁戒，如 Si 电子
    的 X 谷）也算"有"：DFPT 的结论就是谷内 ODP ≈ 0，冻声子/手填的谷内通道照样只作对照。"""
    if not _carrier(params, carrier):
        return False
    return any(r.get("carrier") == carrier and r.get("kind") == "intra" for r in (params or {}).get("raw") or [])


def inter_channels(params, carrier, keep=None, norm=None, dmin=DMIN_EV_A, alias=None):
    """ivs_dfpt_params.json 某载流子 -> ([dp_optical 通道], [D < dmin 丢掉的], [提示])。
    通道：{"from","to","kind": inter|inter_star|intra,"D_eV_A","hw_meV","source": "dfpt_ivs"（谷间）/"dfpt_odp"（谷内，
    [0059-r3]）}——谷内的也在里面，调用方按 kind 分开（DFPT 谷内优先于冻声子/手填，不能两份都进散射率）；alias（S8.4 并掉的重复星）先换成
    保留的标签；keep 给了就只要两端都在 keep 里的（两端不在的记进提示，调用方据此报错）。同一对星换标签后 ħω 撞车的按 D² 相加。"""
    c = _carrier(params, carrier)
    if not c:
        return [], [], []
    norm = norm or (lambda s: s)
    alias = alias or {}
    acc, notes, dropped = {}, [], []
    for ch in c.get("channels") or []:
        a = alias.get(norm(ch["from"]), norm(ch["from"]))
        b = alias.get(norm(ch["to"]), norm(ch["to"]))
        if keep is not None and (a not in keep or b not in keep):
            notes.append("通道 %s>%s 的谷不在现在的谷表里（%s）" % (a, b, ",".join(sorted(keep))))
            continue
        kind = ch["kind"]
        key = (a, b, kind, round(float(ch["hw_meV"]), 2))
        acc[key] = acc.get(key, 0.0) + float(ch["D_eV_A"]) ** 2
    chs = []
    for (a, b, kind, hw), d2 in sorted(acc.items()):
        D = math.sqrt(d2)
        if D < dmin:
            dropped.append({"from": a, "to": b, "kind": kind, "hw_meV": hw, "D_eV_A": D})
            continue
        chs.append({"from": a, "to": b, "kind": kind, "D_eV_A": D, "hw_meV": hw,
                    "source": "dfpt_odp" if kind == "intra" else "dfpt_ivs"})
    return chs, dropped, notes


def adp_params(params, carrier):
    """[0072] ivs_dfpt_params.json 的 "adp" 块里某载流子可用的谷内 ADP 张量（status ok）-> 记录；否则 None。
    记录：{"status": "ok", "valleys": [{"label", "members": [{"k_frac", "D_eV"}]}], "summary": [...], "checks": {...}}。"""
    c = (((params or {}).get("adp") or {}).get("carriers") or {}).get(carrier) or {}
    return c if c.get("status") == "ok" and c.get("valleys") else None


def adp_digest(params, carrier):
    """谷内 ADP 张量的指纹（S8 记下来，张量变了才重排）；不可用 -> None。"""
    c = adp_params(params, carrier)
    if not c:
        return None
    rows = sorted([str(v["label"])] + [round(float(x), 4) for m in v["members"] for r in m["D_eV"] for x in r]
                  for v in c["valleys"])
    return hashlib.sha256(json.dumps(rows, ensure_ascii=False).encode("utf-8")).hexdigest()[:16]


def params_digest(params, carrier):
    """某载流子通道表的指纹（S8/S8.4 记下来，参数变了才重排）；该载流子不可用 -> None。"""
    c = _carrier(params, carrier)
    if not c:
        return None
    rows = sorted([str(ch["from"]), str(ch["to"]), str(ch["kind"]), round(float(ch["D_eV_A"]), 4),
                   round(float(ch["hw_meV"]), 2)] for ch in c.get("channels") or [])
    return hashlib.sha256(json.dumps(rows, ensure_ascii=False).encode("utf-8")).hexdigest()[:16]


def describe(params, carrier, norm=None, keep=None, dmin=DMIN_EV_A, alias=None):
    chs, _d, _n = inter_channels(params, carrier, keep=keep, norm=norm, dmin=dmin, alias=alias)
    if not chs:
        return "无（全部 < %.2f eV/Å）" % dmin
    by = {}
    for ch in chs:
        by.setdefault((ch["from"], ch["to"], ch["kind"]), []).append(ch)
    parts = []
    for (a, b, kind), lst in sorted(by.items()):
        tot = math.sqrt(sum(x["D_eV_A"] ** 2 for x in lst))
        parts.append("%s>%s%s %d 支(√ΣD² %.2f eV/Å)%s" % (a, b, "'" if kind == "inter_star" else "", len(lst), tot,
                                                         " 谷内 ODP" if kind == "intra" else ""))
    return "; ".join(parts)


def valleys_for_consumer(params, carrier):
    """3D 用：JSON 里的谷表（已展开成完整轨道）-> dp_optical 的谷表 [{"label","g","k_frac"}]。"""
    c = _carrier(params, carrier)
    if not c:
        return []
    return [{"label": norm_label(v["label"]), "g": len(v["k_frac"]), "k_frac": [list(map(float, k)) for k in v["k_frac"]]}
            for v in c.get("valleys") or []]


def match_stars(params, carrier, valleys, full3d=False, tol=0.08):
    """DFPT 的谷表 vs 现在 dp_optical 用的谷表（已 merge_coincident）：按 k 中心一一配对星。
    两个星对上 = 成员数相同、每个 DFPT 成员在 tol（分数坐标，模 1；2D 只看面内）内有一个成员对应（现在的谷表有 k_grid 就
    也拿它比，2D 的 S8.2 谷表与 DFPT 用的是同一组网格点，精确相等）。返回 ({DFPT 标签: 现在的标签}, 错误说明或 None)：
    DFPT 的星找不到对应、或现在的谷表多出 DFPT 没有的星，都算错误（谷表变了，参数要按新谷表重算）。"""
    c = _carrier(params, carrier)
    if not c:
        return {}, "%s 没有可用的 DFPT 参数" % carrier

    def _arr(ks):
        return [np.asarray(k, float) for k in ks]

    def _close(a, b):
        d = a - b
        d = d - np.rint(d)
        if not full3d:
            d[2:] = 0.0
        return bool(np.all(np.abs(d) < tol))

    def _same(A, B):
        if len(A) != len(B):
            return False
        used = set()
        for a in A:
            j = next((j for j, b in enumerate(B) if j not in used and _close(a, b)), None)
            if j is None:
                return False
            used.add(j)
        return True
    cur = [(v["label"], [_arr(v["k_frac"])] + ([_arr(v["k_grid"])] if v.get("k_grid") else [])) for v in valleys]
    out, left = {}, list(range(len(cur)))
    for dv in c.get("valleys") or []:
        dk = _arr(dv["k_frac"])
        hit = next((i for i in left if any(_same(dk, ks) for ks in cur[i][1])), None)
        if hit is None:
            return {}, ("DFPT 的星 %s（%d 个成员）在现在的谷表（%s）里找不到 k 中心对得上的星"
                        % (dv["label"], len(dk), ", ".join("%s×%d" % (lab, len(ks[0])) for lab, ks in cur)))
        out[norm_label(dv["label"])] = cur[hit][0]
        left.remove(hit)
    if left:
        return out, ("现在的谷表里 %s 没有 DFPT 参数（DFPT 谷表：%s）" %
                     (", ".join(cur[i][0] for i in left), ", ".join("%s×%d" % (v["label"], len(v["k_frac"]))
                                                                  for v in c.get("valleys") or [])))
    return out, None


def load_consumer(path, mode, want_dim):
    """S8 / S8.4 读 IVS_DFPT -> (params 或 None, 提示或 None)。
    mode：auto（文件在就用）/ on（必须有）/ off。文件缺、读不了、schema 或维度不对：on -> ValueError；auto -> (None, 提示)。"""
    import os
    mode = str(mode or "auto").strip().lower()
    if mode not in ("auto", "on", "off"):
        raise ValueError("IVS_DFPT=%r 无效，只允许 auto / on / off" % mode)
    if mode == "off":
        return None, None
    if not os.path.isfile(path):
        if mode == "on":
            raise ValueError("IVS_DFPT=on，但 %s 不存在：先跑完 eph-qe-cpu 的 ivs_dfpt 组（S10–S13），结果拉回本地后"
                             "（eph-qe-cpu/result/step13_ivs_read/ivs_dfpt_params.json）retry 本步，autozt 会推过来" % path)
        return None, None
    try:
        params = json.loads(open(path, encoding="utf-8").read())
    except ValueError as e:
        if mode == "on":
            raise ValueError("%s 读不了：%s" % (path, e))
        return None, "%s 读不了（%s）：本次不用" % (path, e)
    msg = None
    if params.get("schema") != SCHEMA:
        msg = "%s 的格式 %r 不是本脚本认的 %r：重跑 eph-qe-cpu S13_ivs" % (path, params.get("schema"), SCHEMA)
    elif params.get("dim") != want_dim:
        msg = "%s 是 %s 口径，本步是 %s：材料对不上（或 eph-qe-cpu 的 IVS_DIM 写错）" % (path, params.get("dim"), want_dim)
    if msg:
        if mode == "on":
            raise ValueError(msg)
        return None, msg + "；本次不用"
    return params, None


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()
