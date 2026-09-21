#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step7b_deform_read.py —— amset deform read → deformation.h5（step7b_deform_read）。

run:gen 步骤：在登录节点直接跑，不提交 SLURM。
把 step7_deform 的全部形变单点结果读进 deformation.h5。秒级完成。
产出：本步目录下的 deformation.h5（done_marker）。
"""
import glob
import os
import shlex
import subprocess
import sys
from pathlib import Path

# patch_inplane_kc（2026-09-20）：_inplane_dirs() 用 kc.read_lattice_matrix 判"面内/面外"，
# 但本生成器顶部从未 import ke_common —— 旧代码只在生成出来的 payload 里 import 了它。
# 于是判面内/面外**永远**抛 NameError -> except 把它当成"全部面内" -> 二维项目里
# "面内 ionrelax + 面外 clamped" 被误判成 mixed 口径而被**硬拦**（实测 CrSe2_hex：
# deform-01..04 有 ionrelax、05..09 是面外 clamped，本可放行）。
# 与 gen_step9_deform.py 的写法一致，补上 ke_common 导入。
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, os.getcwd())
import ke_common as kc  # noqa: E402

# =========================== 可改参数区 ===========================
# [SKILL_REV] 版本戳：写进 band_edges.json，与 gen_step12_dpt.py 交叉校验一致。
_SKILL_REV = "2026-09-16-deform-ref-vacuum"
OUTDIR_NAME  = "step7b_deform_read"
DEFORM_DIR   = "step7_deform"
# patch_deform_ref（2026-09-16）：形变势的**参考能级口径**。
#   "vacuum" = 以真空静电势为零点（二维文献通行做法；C2DB 给的也是真空口径）。
#   "core"   = AMSET 自带的"芯态平均势对齐"（官方文档写明，是为**体材料**设计的）。
#   "auto"   = 2D 用 vacuum、其余用 core（默认）。
# 实现：Δ_ij = -dE_vac/dε_ij 这一项**与能带无关**（真空势的应变导数只依赖构型），
# 所以把它加到 deformation.h5 里所有能带、所有 k 点的 D 上即可得到真空口径张量：
#     D_vac = D_core + Δ
# 另写 deformation_vac.h5（core 版原样保留，便于对比与回退）。
# 依据：二维文献普遍以真空为零点定义带能；AMSET 的芯势对齐在 slab 上会系统性偏小
# （实测 CrS2 单层：E1_core=3.53 eV vs E1_vac=5.86 eV，μ 差 2.76 倍）。
DEFORM_REF   = "auto"
# 只有这些应变分量有 LOCPOT 真空势可算；其余分量 Δ 留 0 并在 json 里标注。
_DELTA_COMP  = {"xx": (0, 0), "yy": (1, 1), "zz": (2, 2),
                "yz": (1, 2), "zy": (1, 2), "xz": (0, 2), "zx": (0, 2),
                "xy": (0, 1), "yx": (0, 1)}
# patch_deform_fix：只匹配形变目录。绝不能用 "*deform*"——它会把
# "undeformed" 自己也匹配进去，undeformed 必须单独作为 bulk 传给 read。
DEFORM_GLOB  = "deform-*"
# amset 环境：优先读 step.conf 里 tf 注入的集群 CONDA_SH/AMSET_ENV（setting/<集群>.yaml
# 里配 conda_sh/amset_env，每人按自己的机器配置，脚本不硬编码路径）；
# step.conf 没有（旧项目/直跑脚本）才回退主机探测。
import os as _os
def _amset_env_src():
    try:
        import stepconf as _sc
        _txt = open(_sc.CONF_NAME, encoding="utf-8-sig").read()
        _p = {k.upper(): v for k, v, _ in _sc.parse(_txt, _sc.CONF_NAME).get("params", [])}
        _sh, _env = _p.get("CONDA_SH"), _p.get("AMSET_ENV")
        if _sh and _env:
            return "source %s && conda activate %s" % (_sh, _env)
    except Exception:
        pass
    if _os.path.isdir("/home/wangchaoyue852/miniconda3"):
        return "source /home/wangchaoyue852/miniconda3/etc/profile.d/conda.sh && conda activate amset"
    return "source /public/home/wangchao/miniconda3/etc/profile.d/conda.sh && conda activate amset_clean"
AMSET_ENV_SRC = _amset_env_src()
# =================================================================


# patch_dim_guard：本步不跑 VASP、也不解析结构，所以没有 dim 变量可用。
# 直接从 step1 的 workflow_method.txt 读 DIM=，0D 就带原因退出，
# 免得 -f 强推时抛一句看不懂的"缺 xxx.h5"。
_STEP1_CANDS = ("step1_opt", "step1_std_opt",
                "step1c_PBE_opt", "step1b_PBE_opt", "step1a_PBE_opt")


def _guard_not_0d(cwd, step_name, why):
    from pathlib import Path as _P
    for name in _STEP1_CANDS:
        mf = _P(cwd) / name / "workflow_method.txt"
        if not mf.is_file():
            continue
        for ln in mf.read_text(errors="ignore").splitlines():
            if ln.strip().upper().startswith("DIM="):
                dim = ln.split("=", 1)[1].strip().lower()
                if dim == "0d":
                    sys.exit("[ERROR] %s 不支持 0D 体系。\n"
                             "        原因：%s\n"
                             "        支持的维度：2D, 3D\n"
                             "        若判定有误，检查 %s 的 DIM=。"
                             % (step_name, why, mf))
                return
        return


def _read_dim(cwd):
    """从 step1 的 workflow_method.txt 读 DIM=（2d/3d/0d）。"""
    for name in ("step1_opt", "step1_std_opt", "step1c_PBE_opt",
                 "step1b_PBE_opt", "step1a_PBE_opt"):
        mf = Path(cwd) / name / "workflow_method.txt"
        if not mf.is_file():
            continue
        for ln in mf.read_text(errors="ignore").splitlines():
            if ln.strip().upper().startswith("DIM="):
                return ln.split("=", 1)[1].strip().lower()
    return None


def _resolve_deform_ref(cwd):
    """DEFORM_REF: auto -> 2D 用 vacuum、其余 core；显式值原样返回。"""
    ref = str(DEFORM_REF).lower()
    if ref == "auto":
        return "vacuum" if _read_dim(cwd) == "2d" else "core"
    if ref not in ("vacuum", "core"):
        sys.exit("[ERROR] DEFORM_REF=%r 无效，只允许 auto / vacuum / core" % DEFORM_REF)
    return ref


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    _guard_not_0d(cwd, "step7b_deform_read",
                  "形变势是能带对应变的响应，孤立分子没有能带色散")
    dfm = cwd / DEFORM_DIR
    if not dfm.is_dir():
        sys.exit("[ERROR] 找不到 %s（形变单点没生成？）" % dfm)

    # 校验每个形变子目录都有 vasprun.xml，否则 read 会失败或给错结果
    subs = sorted(p for p in glob.glob(str(dfm / DEFORM_GLOB)) if os.path.isdir(p))
    und = str(dfm / "undeformed")
    if not os.path.isdir(und):
        sys.exit("[ERROR] 缺 %s —— 它是形变势的参考态，没有它 read 无法对齐"
                 % und)
    subs.append(und)
    if not [p for p in subs if p != und]:
        sys.exit("[ERROR] %s 下没有 deform-* 子目录（step7 没跑或跑挂了）" % dfm)
    missing = [os.path.basename(p) for p in subs
               if not (os.path.isfile(os.path.join(p, "vasprun.xml"))
                       or os.path.isfile(os.path.join(p, "vasprun.xml.gz")))]
    if missing:
        sys.exit("[ERROR] 以下形变单点还没算完（缺 vasprun.xml）：%s\n"
                 "        等 step7_deform 全部 done 再跑本步。"
                 % ", ".join(missing[:8]))

    # [PATCH-DP-IONRELAX] 形变势读**离子弛豫**构型（deform-*/ionrelax/，有则用）：
    #   ① 物理：q→0 声学声子 = 均匀应变 + 内坐标弛豫（长波声学模频率→0，原子坐标
    #      绝热跟随；石墨烯 q→0 声学声子的第一性原理处理同此）；
    #   ② 一致性（更硬）：ADP 核分母是 Christoffel 刚度＝含离子弛豫的弹性常数
    #      （VASP 的 TOTAL ELASTIC MODULI），band_edges/[8.2 DPT] 也一直读 ionrelax，
    #      分子分母必须同口径。
    #   实测 MoS2 单层 xx 0.5%：电子 E1_vac clamped 6.58 eV vs ionrelax 8.51 eV（差 29%），
    #   这是真实物理差别，论文方法部分要写明用哪一种。
    # ★ 回退不能静默（2026-09-16 用户要求）：没有 ionrelax/ 的构型是**离子固定**口径，
    #   与其它项目的离子弛豫口径不同源。跨项目的 ADP 数值本来就不可直接比，但"没人
    #   察觉自己用的是哪一套"是更糟的失效模式 —— 所以这里 WARN + 落盘实际口径。
    _geom_map = {}

    def _dp_folder(d):
        dd = Path(d)
        if (dd / "ionrelax" / "vasprun.xml").is_file():
            _geom_map[str(dd.relative_to(dfm))] = "ionrelax"
            return dd / "ionrelax"
        _geom_map[str(dd.relative_to(dfm))] = "clamped"
        return dd

    dp_subs = [_dp_folder(p) for p in subs if os.path.abspath(p) != os.path.abspath(und)]
    _rel = [str(p.relative_to(dfm)) if str(p).startswith(str(dfm)) else str(p)
            for p in dp_subs]
    print("[..] 形变构型（%d 个，ionrelax 优先）：%s"
          % (len(_rel), ", ".join(_rel[:6]) + (" ..." if len(_rel) > 6 else "")))
    # ---- 口径分级（2026-09-16 二次修订，按维度区分处理）------------------
    # 现实情况：step7 只对**面内**分量做了离子弛豫，所以一个项目里常常是
    # "面内 ionrelax + 面外 clamped"。这不是可以含糊过去的细节：
    #   · 二维核只用面内分量 → 面内全 ionrelax 即可放行，但标签必须写清楚；
    #   · 三维核要全部分量同源 → 出现混合口径直接报错（AMSET 会把不同来源的分量
    #     平均进同一个张量，结果没有物理意义）；
    #   · 面内只要有**一个**分量是离子固定口径，二维也不能放行。
    _clamped = sorted(k for k, v in _geom_map.items() if v == "clamped")
    _relaxed = sorted(k for k, v in _geom_map.items() if v == "ionrelax")

    def _inplane_dirs():
        """判每个形变目录的主导应变分量是不是面内（0/1 方向）。取不到返回 None。"""
        try:
            import numpy as _np
            und = _np.array(kc.read_lattice_matrix(dfm / "undeformed" / "POSCAR"))
            out = {}
            for k in _geom_map:
                L = _np.array(kc.read_lattice_matrix(dfm / k / "POSCAR"))
                F = _np.transpose(_np.dot(_np.linalg.inv(und), L))
                E = (_np.dot(F.T, F) - _np.eye(3)) / 2.0
                i, j = _np.unravel_index(_np.argmax(_np.abs(E)), E.shape)
                out[k] = (i, j) in ((0, 0), (1, 1), (0, 1), (1, 0))
            return out
        except Exception as _e:                                # noqa: BLE001
            print("[WARN] 应变分量判面内/面外失败（%s）；保守处理：全部当作面内" % _e)
            return None

    _ip = _inplane_dirs()
    _clamped_in = sorted(k for k in _clamped if (_ip is None or _ip.get(k, True)))
    _clamped_out = sorted(k for k in _clamped if not (_ip is None or _ip.get(k, True)))
    _dim_now = _read_dim(cwd)
    if not _clamped:
        _def_geom = "ionrelax"
    elif not _relaxed:
        _def_geom = "clamped"
    elif _clamped_in:
        _def_geom = "mixed"
    else:
        _def_geom = "ionrelax(in-plane)+clamped(out-of-plane)"

    if _def_geom == "mixed":
        sys.exit(
            "[ERROR] 形变势构型口径 mixed：面内分量里 %s 是【离子固定】、其余是离子弛豫。\n"
            "        二维核只用面内分量，面内口径必须一致；请给这些构型补 ionrelax/ "
            "（IBRION=2 弛豫内坐标）后重跑本步。" % ", ".join(_clamped_in))
    if _dim_now != "2d" and _clamped and _relaxed:
        sys.exit(
            "[ERROR] 三维体系的形变势构型口径不一致：%d 个 ionrelax、%d 个 clamped。\n"
            "        三维核需要全部分量同源（AMSET 会把它们平均进同一个张量），"
            "请补齐 ionrelax/ 后重跑本步；确实不需要弛豫就在 step7 里关掉。"
            % (len(_relaxed), len(_clamped)))
    if _def_geom == "clamped":
        print("[WARN] 全部 %d 个形变构型都是【离子固定】口径（没有 ionrelax/）：与弹性常数"
              "（TOTAL ELASTIC MODULI，含离子弛豫）不同源，ADP 绝对值跨项目不可直接比。"
              % len(_geom_map))
    elif _def_geom.startswith("ionrelax(in-plane)"):
        print("[..] 面内分量取自离子弛豫构型、面外分量取自离子固定构型：二维核只用面内分量，"
              "可用；三维核若复用同一份 h5 需注意口径（本体系 DIM=%s）。" % (_dim_now or "?"))
        print("     面外 clamped 构型：%s" % ", ".join(_clamped_out[:8]))
    print("[..] 形变势构型口径 deform_geometry=%s（落盘到 band_edges.json）" % _def_geom)
    # amset deform read 在形变目录里跑，产出 deformation.h5，再挪到本步目录
    # patch_deform_fix：amset 的签名是 read(bulk_folder, deformation_folders...)，
    # 【未形变的必须排第一个】。原来写成 `read *deform* undeformed`，既顺序
    # 颠倒，又因为 "*deform*" 会匹配到 "undeformed" 自己，把第一个形变目录
    # 当成了 bulk。这个错不报异常，只会安静地给出错误的形变势。
    cmd = ("%s && cd %s && amset deform read undeformed %s "
           ">> deform_read.log 2>&1"
           % (AMSET_ENV_SRC, str(dfm),
              " ".join(shlex.quote(r) for r in _rel)))
    print("[..] amset deform read ...")
    rc = subprocess.run(["bash", "-lc", cmd]).returncode
    h5 = dfm / "deformation.h5"
    if rc != 0 or not h5.is_file():
        sys.exit("[ERROR] amset deform read 失败，看 %s/deform_read.log" % dfm)
    dst = out / "deformation.h5"

    # ---- patch_dpt_edgefix：用 amset 权威带边定位，写 band_edges.json ----
    # 旧实现（gen_step12_dpt._band_edge_index）拿 step3_uniform（ISYM=2 IBZ 网格）
    # 的 (band,k) 索引直接去索引 deformation.h5（ISYM=0 全 BZ 网格 + amset 按
    # energy_cutoff 截断后的 band 轴），两套网格根本不是一套，索引必然错位；
    # 错位后静默走"全带中位数"兜底，electron/hole 取到同一个背景值。
    # 这里在 h5 生成的同时，用 amset 自己的 get_vbm/get_cbm 定位带边，band 轴
    # 经 get_ibands 映射到 h5 的截断轴，k 点按 frac-coords 匹配 h5 网格，
    # 把每个载流子的带边 E1（D_xx/D_yy 与面内均值）写进 band_edges.json，
    # step8.2_dpt 的 gen_step12_dpt.py 优先读它，不再自己猜索引。
    # 注意：必须在 h5 移入 out 之前跑（amset read 刚在 dfm 里产出 h5）；
    # json 用绝对路径写到 out，读的相对路径都相对于 dfm。
    _be_code = r'''
import json, numpy as np, glob
from amset.deformation.io import load_deformation_potentials, parse_calculation
from amset.electronic_structure.common import get_ibands
from amset.constants import defaults

import sys
import os as _os
from pathlib import Path as _P
# [PATCH-IONRELAX] 共用 ke_common.resolve_strain_pairs（与 gen_step9_deform.py 单一真源）
sys.path.insert(0, str(_P.cwd().resolve().parent))   # 项目根（本步 cwd 是 step7_deform）
import ke_common
outdir = sys.argv[1]

def _as_dict(edge):
    edge = edge if isinstance(edge, dict) else {
        "energy": edge.energy, "band_index": edge.band_index,
        "kpoint_index": edge.kpoint_index}
    edge = dict(edge)
    # [PATCH-BANDIDX] pymatgen 的 BandStructure.get_cbm/get_vbm 对**非自旋极化**
    #   体系把 band_index 返回成普通 list（只有 AMSET 自己的 Bandstructure 才是
    #   {Spin: [band_idx]}）。两种都归一成 {Spin: [band_idx]}，自旋键取自
    #   deformation.h5 的实际自旋道（与 ib / dp 的键同一套），否则下面的
    #   ib[spin.name] / dp[spin] 会 KeyError 或被静默错取。
    #   实测：2026-09-21 Mg4C60_relax3 S7.1_read 报
    #   AttributeError: 'list' object has no attribute 'items'。
    _bi = edge.get("band_index")
    if not isinstance(_bi, dict):
        _spins = list(dp.keys())
        if not _spins:
            raise RuntimeError("deformation.h5 里没有任何自旋道，无法归一 band_index")
        edge["band_index"] = {_spins[0]: [int(b) for b in np.atleast_1d(_bi)]}
    # [PATCH-BANDIDX] kpoint_index 同理摊平成一维 list（下面按 list 用：
    #   ki[0] 与 for k_idx in e["kpoint_index"]）。
    _ki = edge.get("kpoint_index")
    if isinstance(_ki, dict):
        _flat = []
        for _v in _ki.values():
            _flat.extend(int(x) for x in np.atleast_1d(_v))
        edge["kpoint_index"] = _flat
    elif _ki is not None:
        edge["kpoint_index"] = [int(x) for x in np.atleast_1d(_ki)]
    return edge

dp, kpoints, structure = load_deformation_potentials("deformation.h5")
bulk = parse_calculation("undeformed")
bs = bulk["bandstructure"]
ib = get_ibands(defaults["energy_cutoff"], bs)
ib = {s.name: np.asarray(v) for s, v in ib.items()}
# [PATCH-DPT-R4] band 轴断言：h5 的 band 轴由 amset deform read 当时的
# energy_cutoff 决定，这里用同一 defaults 重建 ibands；两者不一致说明
# 调用参数漂移了（如 read 传了 -e），必须显式失败而不是静默错位。
for _sp, _d in dp.items():
    if _d.shape[0] != len(ib[_sp.name]):
        raise RuntimeError(
            "band 轴与 h5 不一致（h5=%d, ibands[%s]=%d）：检查 deform read 的 "
            "-e/energy_cutoff" % (_d.shape[0], _sp.name, len(ib[_sp.name])))

# [PATCH-METAL] pymatgen 对金属/半金属的 get_cbm/get_vbm 返回 energy=None、
#   band_index=[]（bs.is_metal() 为真）；而带边形变势 E1/ADP 的前提就是"有带边"。
#   这里显式判定并给出可读结论，而不是让下游 float(None) 抛一个看不懂的
#   TypeError，也不是静默写一份 n_hits=0 的 band_edges.json（后者会被 step8.2
#   当成正常产物继续用）。实测：2026-09-21 Mg4C60_relax3（PBE gap=-0.114 eV、
#   HSE gap=-0.124 eV，两级别都判金属）。
def _is_metal(obj):
    f = getattr(obj, "is_metal", None)
    try:
        return bool(f()) if callable(f) else False
    except Exception:
        return False

if _is_metal(bs):
    raise RuntimeError(
        "体系为金属（pymatgen is_metal=True，undeformed bandstructure 无 VBM/CBM）："
        "带边形变势 E1 / ADP 对金属不适用 —— 本步不产出 band_edges.json，"
        "也不应继续 S8.2_dpt。若材料确实有隙，先查 step2 的带隙判别是否被"
        "零权重路径点或费米面穿过污染。")

out = {}
for carrier, edge in (("electron", bs.get_cbm()), ("hole", bs.get_vbm())):
    e = _as_dict(edge)
    ene = e["energy"]
    if ene is None:
        raise RuntimeError(
            "%s 带边能量为 None（无带边可定义）—— 带边形变势 E1/ADP 不适用" % carrier)
    hits = []
    for spin, bidxs in e["band_index"].items():
        for b_global in bidxs:
            loc = np.where(ib[spin.name] == b_global)[0]
            if len(loc) == 0:
                continue
            b_local = int(loc[0])
            for k_idx in e["kpoint_index"]:
                kf = np.array(bs.kpoints[k_idx].frac_coords)
                # [PATCH-DPT-R3] k 点周期回卷 + 容差：k=0.999 与 -0.001 是
                # 同一点；不回卷会 argmin 到无关 k（与旧 bug 同类失败）。
                dk = kpoints - kf
                dk -= np.round(dk)
                dd = np.linalg.norm(dk, axis=1)
                kh5 = int(np.argmin(dd))
                if dd[kh5] > 1e-4:
                    raise RuntimeError(
                        "带边 k 点 %s 在 h5 网格找不到（最近距离 %.3g）"
                        % (np.round(kf, 4), dd[kh5]))
                t = dp[spin][b_local, kh5]
                hits.append({"spin": spin.name, "band_global": int(b_global),
                             "band_h5": b_local, "k_h5": int(kh5),
                             "k_frac": np.round(kf, 6).tolist(),
                             "E1_xx_eV": round(float(t[0, 0]), 4),
                             "E1_yy_eV": round(float(t[1, 1]), 4),
                             "E1_iso_eV": round(float(abs((t[0, 0] + t[1, 1]) / 2)), 4)})
    if hits:
        # [PATCH-DPT-R2] 简并 spin×band×k 的命中做 |D| 平均（h5 值已非负），
        # 并报告离散度；不再取 hits[0]（依赖 dict 迭代顺序、不可复现）。
        xs = [h["E1_xx_eV"] for h in hits]
        ys = [h["E1_yy_eV"] for h in hits]
        iso = [h["E1_iso_eV"] for h in hits]
        avg = lambda v: round(float(sum(v)) / len(v), 4)
        spread = lambda v: (max(v) - min(v)) if v else 0.0
        out[carrier] = {
            "edge_energy_eV": round(float(ene), 4),
            "n_hits": len(hits),
            "E1_xx_eV": avg(xs), "E1_yy_eV": avg(ys), "E1_iso_eV": avg(iso),
            "spread_xx": round(spread(xs), 4), "spread_yy": round(spread(ys), 4),
            "hits": hits,
        }
        if spread(xs) > 0.5 or spread(ys) > 0.5:
            print("[WARN] %s 带边简并分量离散度大：xx spread=%.3f yy spread=%.3f"
                  % (carrier, spread(xs), spread(ys)))
    else:
        out[carrier] = {"edge_energy_eV": round(float(ene), 4),
                        "n_hits": 0, "hits": []}

# ---- [PATCH-DPT-R6] 真空能级对齐（Qiao 黑磷路线）----
# amset h5 的 DP 是 ⟨|D|⟩（先取模再平均），无符号、且参考是平均芯势——与文献
# Eq.(9) E1=dE_edge/dγ 口径不同。这里从原始数据重算：
#   E1_vac = |(ΔE_edge - ΔE_vac)/γ|      ← 生成端已取 |·|，写出的字段恒非负
# ΔE_vac 取 undeformed 与 deform-N 的 LOCPOT 真空区平面平均势之差（刚性平移）。
# [A2] 应变分量从结构反解（deform-* 的 POSCAR 晶格 vs undeformed），不硬编码编号。
# [A3/A4] k 点回卷+容差、自旋道按 band_index 选，与上方 amset 带边块共用同一套。
# [B1] 真空区取离所有原子最远的分数间隙中心 ±VAC_WINDOW_HALF，检查平坦度。
# [B2] 形变后校验带边未翻转，结果写进 json 的 edge_flip（不只 print）。
# [B5] 仅 2D 执行；3D 是合法跳过。
# ---- 本轮修复 ----
# [C1] 只有两类是合法跳过：LOCPOT 缺失/未生成（OSError）与 3D 无真空（_VacSkip）。
#      其余守卫（应变反解不全、k 点错位、真空不平坦、dim 读不到）一律抛出——
#      静默降级会产出「看起来正常但其实是 amset 芯势口径」的 dpt_result.json，
#      而 5 个材料自动跑时没人会逐个核 provenance 字符串。
# [C3] dim 从项目根的 step1*/workflow_method.txt 读，与 gen_step12_dpt._read_dim
#      同一套约定（本步 cwd 是 step7_deform，不能只 open 相对路径）。
# [C4] LOCPOT 平面平均势 / vasprun 本征值全部缓存，循环不变量提到内层循环外。
# [C7] 真空窗口按周期拼接（slab 居中时最大间隙跨 z=1 边界），并做半宽敏感性扫描。
# [C8] 每个构型用各自目录的 POSCAR 定位真空窗口。

VAC_WINDOW_HALF = 0.25          # 生产窗口：真空间隙中心 ± 间隙宽度×此值
VAC_WINDOW_SCAN = (0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)   # [C7] 敏感性扫描
VAC_FLAT_TOL_EV = 1e-3          # 真空区平面平均势平坦度上限（std）
_STEP1_CANDS = ("step1_opt", "step1_std_opt",
                "step1c_PBE_opt", "step1b_PBE_opt", "step1a_PBE_opt")


class _VacSkip(Exception):
    """[C1] 合法跳过真空对齐（3D 无真空），与「守卫失败」区分开。"""


_resolved_paths = {}   # 模块级兜底：3D 走 _VacSkip 提前跳过 try 内初始化，结尾 read_paths 仍可引用

try:
    from pymatgen.io.vasp.outputs import Locpot, Vasprun
    from pymatgen.core import Structure
    from pymatgen.analysis.elasticity.strain import Deformation
    import os as _os
    from pathlib import Path as _P

    # [C3] 维度守卫（[B5]）：从项目根的 step1*/workflow_method.txt 读 DIM=，
    # 与 gen_step12_dpt._read_dim 同一套约定（本步 cwd 是 step7_deform）
    _cwd = _P.cwd().resolve().parent              # 项目根
    _dim = None
    for _name in _STEP1_CANDS:
        _mf = _cwd / _name / "workflow_method.txt"
        if not _mf.is_file():
            continue
        for _ln in _mf.read_text(errors="ignore").splitlines():
            if _ln.strip().upper().startswith("DIM="):
                _dim = _ln.split("=", 1)[1].strip().lower()
                break
        if _dim:
            break
    if _dim is None:
        raise RuntimeError("读不到维度：检查项目根下 %s 的 workflow_method.txt"
                           % " / ".join(_STEP1_CANDS))
    if _dim != "2d":
        raise _VacSkip("dim=%s：真空对齐仅 2D 有意义" % _dim)

    def _k_index(kpts, kf, tol=1e-4):
        """[A3] k 点周期回卷 + 容差匹配（与 amset 带边块同款）。"""
        dk = kpts - kf
        dk -= np.round(dk)
        dd = np.linalg.norm(dk, axis=1)
        k = int(np.argmin(dd))
        if dd[k] > tol:
            raise RuntimeError("k 点 %.4s 不在网格（最近 %.3g）" % (kf, dd[k]))
        return k

    # [C4] vasprun 本征值缓存（key=folder），每构型只解析一次
    _vr_cache = {}

    _resolved_paths = {}   # folder -> {filename: "ionrelax"|"rigid"}，写进 json 可查

    def _resolved(folder, filename):
        """[PATCH-IONRELAX] 优先读 folder/ionrelax/filename（离子弛豫构型），
        否则读本级。E1 用弛豫后构型（IBRION=2 弛豫内坐标），amset h5 仍用
        本级刚性单点——两套口径并列，互不影响。实际读取口径记进 _resolved_paths。"""
        p = _os.path.join(folder, "ionrelax", filename)
        if _os.path.isfile(p):
            _resolved_paths.setdefault(folder, {})[filename] = "ionrelax"
            return p
        _resolved_paths.setdefault(folder, {})[filename] = "rigid"
        return _os.path.join(folder, filename)

    def _edge_raw(folder, kf, band, spin):
        """[A3/A4+C4] 读带边 (energy, occupancy)：回卷+容差匹配 k，自旋道显式选。
        
        [C4] 缓存 vasprun 解析结果（spin×band 双重循环里会反复读同一 folder）。"""
        if folder not in _vr_cache:
            v = Vasprun(_resolved(folder, "vasprun.xml"),
                        parse_dos=False, parse_potcar_file=False)
            _vr_cache[folder] = (np.array(v.actual_kpoints), v.eigenvalues)
        kpts, eigmap = _vr_cache[folder]
        k = _k_index(kpts, kf)
        ev = np.asarray(eigmap[spin])       # 自旋道本征值 (nk, nb, 2)
        return float(ev[k, band, 0]), float(ev[k, band, 1])

    def _vac_window(structure, axis=2, half_ratio=VAC_WINDOW_HALF):
        """[B1+C7] 沿真空轴找离所有原子最远的分数间隙，返回中心 ±half_ratio 窗口。
        
        返回契约（统一 4 元组）：
          不跨边界 → (a, b, None, None)，窗口 = 单段 [a, b]
          跨边界   → (a, 1.0, 0.0, b)，窗口 = 两段 [a, 1] + [0, b]
        判据在 _vac_level 里用 win[2] is None 区分（不是 len(win)）。
        回卷间隙（f[0]+1-f[-1]）也参与选最大，所以 slab 居中时能正确取跨 z=1 的真空。"""
        fracs = np.array([s.frac_coords for s in structure])
        f = np.sort(fracs[:, axis] % 1.0)
        gaps = np.diff(f).tolist() + [f[0] + 1.0 - f[-1]]   # n-1 个相邻 + 回卷
        imax = int(np.argmax(gaps))
        lo = f[imax]
        # 回卷间隙被选中时 hi = f[0]+1.0（> f[-1] = lo，故 lo<hi 恒成立）
        hi = f[imax + 1] if imax < len(f) - 1 else f[0] + 1.0
        c = (lo + hi) / 2.0
        w = (hi - lo) * half_ratio
        a, b = c - w, c + w
        if b > 1.0:            # 窗口右端越界，回卷到 [0, b-1]
            b -= 1.0
        if a < 0.0:
            a += 1.0
        if a > b:              # 回卷后左>右 → 两段
            return a, 1.0, 0.0, b
        return a, b, None, None

    # [C4] LOCPOT 平面平均势缓存（key=folder），每构型只解析一次
    _vz_cache = {}
    # [C8] POSCAR 结构缓存（由调用者管理，在 ax 循环里按需读取）
    _st_cache = {}

    def _vac_level(folder, structure, axis=2, half_ratio=VAC_WINDOW_HALF):
        """[B1+C7+C8] LOCPOT 真空区平面平均势，平坦度 std<VAC_FLAT_TOL_EV 否则报错。
        
        [C8] structure 参数必须是 folder 对应的结构（调用者负责读取和缓存）。
        [C7] 窗口可能跨周期边界（lo>1 或 hi<lo），周期拼接后取平均。"""
        if folder not in _vz_cache:
            lp = Locpot.from_file(_resolved(folder, "LOCPOT"))
            _vz_cache[folder] = np.asarray(lp.get_average_along_axis(axis))
        vz = _vz_cache[folder]
        n = len(vz)
        win = _vac_window(structure, axis, half_ratio)
        if win[2] is None:       # 不跨边界 [a, b]（单段）
            a = int(np.clip(np.floor(win[0] * n), 0, n - 1))
            b = int(np.clip(np.ceil(win[1] * n), a + 1, n))
            seg = vz[a:b]
        else:                    # 跨边界 [a, 1] + [0, b]，两段拼接
            a1 = int(np.clip(np.floor(win[0] * n), 0, n - 1))
            b1 = n
            a2 = 0
            b2 = int(np.clip(np.ceil(win[3] * n), 1, n))
            seg = np.concatenate([vz[a1:b1], vz[a2:b2]])
        if float(np.std(seg)) > VAC_FLAT_TOL_EV:
            raise RuntimeError("真空区平面平均势不平坦 std=%.3g eV（可能取到 slab 尾巴）"
                               % np.std(seg))
        return float(np.mean(seg))

    # [A2] 应变配对反解：与 gen_step9_deform.py 共用 ke_common.resolve_strain_pairs
    # （单一真源，避免 gen 建 ionrelax/ 与 step7b 找 ionrelax/ 各写一份反解而错位）
    _pairs, _strain_mag = ke_common.resolve_strain_pairs(".")
    if _pairs is None:
        raise RuntimeError("应变反解不完整：需 ±0.5% 面内形变各一对")

    for carrier, edge in (("electron", bs.get_cbm()), ("hole", bs.get_vbm())):
        e = _as_dict(edge)
        bi = e["band_index"]; ki = e["kpoint_index"]
        kf = np.array(bs.kpoints[ki[0]].frac_coords)
        # [C4] 循环不变量：undeformed 的能级与真空势（与 spin×band 无关）
        e0_cache = {}       # (spin, band) -> (energy, occ)
        for spin, bidxs in bi.items():
            for band in bidxs:
                e0_cache[(spin, band)] = _edge_raw("undeformed", kf, band, spin)
        v0 = _vac_level("undeformed", structure)
        # [C6] 占据态翻转累计（写进 json，不只 print）
        edge_flip = []
        vac = {}
        dvac_map = {}       # [patch_deform_ref] 带符号的 dE_vac/dε（与能带无关）
        vac_scan = {}       # [C7] 窗口敏感性扫描：{half_ratio: {ax: E1}}
        for ax, (dp_, dm_) in _pairs.items():
            g = _strain_mag[ax] if _strain_mag[ax] > 1e-6 else 0.005
            # [C8] 读 dp_/dm_ 各自的 POSCAR 并缓存（形变后原子位置不同；
            # ionrelax 存在时用弛豫后结构，真空窗口随离子位置移动）
            if dp_ not in _st_cache:
                _st_cache[dp_] = Structure.from_file(_resolved(dp_, "POSCAR"))
            if dm_ not in _st_cache:
                _st_cache[dm_] = Structure.from_file(_resolved(dm_, "POSCAR"))
            # [C4] 形变构型的真空势提到 spin×band 循环外（只依赖构型，不依赖带）
            v_p = _vac_level(dp_, _st_cache[dp_])
            v_m = _vac_level(dm_, _st_cache[dm_])
            dvac = ((v_p - v0) / g + (v_m - v0) / (-g)) / 2
            dvac_map[ax] = float(dvac)      # [patch_deform_ref] 带符号保存
            # [A4/B2] 遍历所有自旋道×带（简并时平均）；形变后校验该带边占据态未翻转
            vals = []
            for spin, bidxs in bi.items():
                for band in bidxs:
                    e0, o0 = e0_cache[(spin, band)]
                    e_p, o_p = _edge_raw(dp_, kf, band, spin)
                    e_m, o_m = _edge_raw(dm_, kf, band, spin)
                    # [B2+C6] 占据态翻转校验：undeformed 的带边占据态应与形变后一致
                    # （electron 应恒空 o≈0、hole 应恒满 o≈1）；翻转说明带序交换
                    if (o_p - 0.5) * (o0 - 0.5) < 0 or (o_m - 0.5) * (o0 - 0.5) < 0:
                        _msg = "%s %s band=%d spin=%s 形变后占据态翻转" % (
                            ax, carrier, band, spin.name)
                        print("[WARN] %s（带序交换）" % _msg)
                        edge_flip.append(_msg)
                    raw = ((e_p - e0) / g + (e_m - e0) / (-g)) / 2
                    vals.append(abs(raw - dvac))
            vac[ax] = round(float(sum(vals)) / len(vals), 4) if vals else None
            # [C7] 窗口半宽敏感性扫描（只在生产窗口之外试，生产值就是 vac[ax]）
            for _hr in VAC_WINDOW_SCAN:
                if abs(_hr - VAC_WINDOW_HALF) < 1e-9:
                    continue        # 已算过
                try:
                    # [Bug B] v0 也用同一 _hr 重算（否则扫描的参考窗口宽与 _vp/_vm 不一致，
                    # 会污染「窗口是否影响 dE_vac/dγ」的结论）
                    _v0 = _vac_level("undeformed", structure, half_ratio=_hr)
                    _vp = _vac_level(dp_, _st_cache[dp_], half_ratio=_hr)
                    _vm = _vac_level(dm_, _st_cache[dm_], half_ratio=_hr)
                    _dv = ((_vp - _v0) / g + (_vm - _v0) / (-g)) / 2
                    _scan_vals = []
                    for spin, bidxs in bi.items():
                        for band in bidxs:
                            e0, o0 = e0_cache[(spin, band)]
                            e_p, o_p = _edge_raw(dp_, kf, band, spin)
                            e_m, o_m = _edge_raw(dm_, kf, band, spin)
                            _r = ((e_p - e0) / g + (e_m - e0) / (-g)) / 2
                            _scan_vals.append(abs(_r - _dv))
                    if _hr not in vac_scan:
                        vac_scan[_hr] = {}
                    vac_scan[_hr][ax] = (round(float(sum(_scan_vals)) / len(_scan_vals), 4)
                                         if _scan_vals else None)
                except RuntimeError as _se:
                    # [C1] 只吞「该窗口平坦度不过」——这是扫描的正常结果之一；
                    # k 点错位等守卫仍必须炸出来，否则 C1 的漏洞在扫描路径复活
                    if "不平坦" not in str(_se):
                        raise
        if carrier in out and vac.get("xx") is not None and vac.get("yy") is not None:
            out[carrier]["E1_vac_xx_eV"] = vac["xx"]
            out[carrier]["E1_vac_yy_eV"] = vac["yy"]
            out[carrier]["E1_vac_iso_eV"] = round((vac["xx"] + vac["yy"]) / 2, 4)
            if edge_flip:       # [C6]
                out[carrier]["edge_flip"] = edge_flip
            if vac_scan:        # [C7]
                out[carrier]["vac_window_scan"] = {
                    str(_hr): {"xx": _v.get("xx"), "yy": _v.get("yy"),
                               # 用 is not None（0.0 是合法值，不能被 falsy 判掉）
                               "iso": (round((_v["xx"] + _v["yy"]) / 2, 4)
                                       if _v.get("xx") is not None
                                       and _v.get("yy") is not None else None)}
                    for _hr, _v in sorted(vac_scan.items())}
                # 立即计算该 carrier 的窗口扫描评判
                vals = [v.get("iso") for v in out[carrier]["vac_window_scan"].values()
                        if v.get("iso") is not None]
                if len(vals) >= 2:
                    avg = sum(vals) / len(vals)
                    std = (sum((x - avg) ** 2 for x in vals) / len(vals)) ** 0.5
                    cv = std / avg if avg > 1e-9 else 0.0
                    diffs = [vals[i+1] - vals[i] for i in range(len(vals)-1)]
                    sign_flips = sum(1 for i in range(len(diffs)-1)
                                     if diffs[i] * diffs[i+1] < 0)
                    _verdict = {
                        "n_windows": len(vals), "mean_iso_eV": round(avg, 4),
                        "std_eV": round(std, 4), "CV": round(cv, 4),
                        "sign_flips": sign_flips,
                        "recommend": ("窗口依赖弱（CV < 5%），生产值可靠" if cv < 0.05 else
                                      "窗口依赖中等（CV 5-15%），建议扩大扫描或取中位数" if cv < 0.15 else
                                      "窗口依赖强（CV > 15%），真空区可能不平坦/形变应变反解有误")}
                    if sign_flips >= 2:
                        _verdict["recommend"] += "；振荡（%d 次翻转，窗口可能伸进 slab）" % sign_flips
                    out[carrier]["scan_verdict"] = _verdict
                    _rec = _verdict["recommend"]
                else:
                    _verdict = {}
                    _rec = ""
            else:
                _verdict = {}
                _rec = ""
            print("[OK] %s 真空对齐 E1: xx=%.3f yy=%.3f (window_half=%.2f, scan=%d%s)"
                  % (carrier, vac["xx"], vac["yy"], VAC_WINDOW_HALF, len(vac_scan),
                     (" CV=%.1f%% %s" % (_verdict.get("CV", 0)*100, _rec.split("；")[0])
                      if _verdict else "")))
    # [C7] scan_verdict 已在 carrier 循环内计算并写进各 carrier 的字段
    out["vac_align"] = {"status": "ok", "window_half": VAC_WINDOW_HALF,
                        "flat_tol_eV": VAC_FLAT_TOL_EV,
                        "strain_mag": {k: round(v, 6) for k, v in _strain_mag.items()},
                        "pairs": {k: list(v) for k, v in _pairs.items()},
                        # [patch_deform_ref] 带符号 dE_vac/dε，供 deformation_vac.h5 用
                        "dvac_eV_per_unit_strain": {k: round(v, 6)
                                                    for k, v in dvac_map.items()}}
    # [patch_deform_ref] Δ_ij = -dE_vac/dε_ij（与能带无关，可平移到所有带）
    # ★ 这段跑在**嵌入的独立脚本**里（_be_code，见文件头注释），只能看到自己
    #   定义的名字；外层模块的 _DELTA_COMP 在这里是未定义名（NameError:
    #   name '_DELTA_COMP' is not defined，2026-09-16 MoS2 实测）。所以就地
    #   定义一份，与文件顶部的同名常量必须保持一致。
    _delta_comp = {"xx": (0, 0), "yy": (1, 1), "zz": (2, 2),
                   "yz": (1, 2), "zy": (1, 2), "xz": (0, 2), "zx": (0, 2),
                   "xy": (0, 1), "yx": (0, 1)}
    _delta = [[0.0] * 3 for _ in range(3)]
    _applied = []
    for _ax, _dv in dvac_map.items():
        _ij = _delta_comp.get(str(_ax).lower())
        if _ij is None:
            continue
        _delta[_ij[0]][_ij[1]] = _delta[_ij[1]][_ij[0]] = -float(_dv)
        _applied.append(str(_ax))
    out["deform_ref"] = {
        "delta_eV_per_unit_strain": _delta,
        "dvac_eV_per_unit_strain": {k: round(v, 6) for k, v in dvac_map.items()},
        "components_applied": _applied,
        "formula": "D_vac = D_core + delta,  delta_ij = -dE_vac/deps_ij",
        "note": ("Δ 与能带无关（真空势的应变导数只依赖构型），加到 deformation.h5 的"
                 "所有带×所有 k 上即得真空口径；未施加的应变分量 Δ 留 0。"),
    }
# [C1] 只有这两类是合法跳过；其余守卫失败必须炸出来（静默降级 = R5 要根除的失败模式）
except _VacSkip as _e:
    print("[SKIP] %s" % _e)
    out["vac_align"] = {"status": "skipped", "reason": str(_e)}
except (FileNotFoundError, OSError) as _e:
    # LOCPOT/POSCAR 没生成（LVHAR 未开、单点没跑完）——合法跳过，consumer 端会硬报错
    print("[WARN] 真空对齐跳过（%s: %s）：LOCPOT 缺失或 LVHAR 未开" % (
        type(_e).__name__, _e))
    out["vac_align"] = {"status": "skipped",
                        "reason": "%s: %s" % (type(_e).__name__, _e)}
except ImportError as _e:
    print("[WARN] 真空对齐跳过（缺依赖 %s）" % _e)
    out["vac_align"] = {"status": "skipped", "reason": "ImportError: %s" % _e}
# 注意：这里【故意】不接 RuntimeError / 其它 Exception ——
#   应变反解不完整、k 点不在网格、真空区不平坦、维度读不到，
#   全部是「结果不可信」而不是「本体系不适用」，必须让本步 FAIL、
#   让 tf 判 error，而不是写出一份少了 E1_vac_* 字段的 band_edges.json。

# [SKILL_REV] 版本戳 + 实际读取路径（ionrelax 还是 rigid），consumer 交叉校验
out["skill_rev"] = sys.argv[2] if len(sys.argv) > 2 else "UNKNOWN"
out["read_paths"] = {k: dict(v) for k, v in _resolved_paths.items()}

# [体系判别] 带重叠时"带边"只是形式上的极值点 —— 标注出来，避免被当物理带边用。
# 判据来自 step2_bandgap/step2.15_discriminant/discriminant.json（带指标法）。
try:
    from pathlib import Path as _DP
    import json as _DJ
    _dp = _DP(outdir).parents[1] / "step2_bandgap" / "step2.15_discriminant" / "discriminant.json"
    if _dp.is_file():
        _dd = _DJ.loads(_dp.read_text())
        _lab = _dd.get("label", "")
        _eff = _dd.get("effective") or {}
        _gap = _eff.get("gap_eV")
        out["discriminant"] = {
            "label": _lab,
            "functional": _eff.get("functional"),
            "gap_eV": _gap,
            "decisive": (_dd.get("block") or {}).get("decisive"),
            "source": str(_dp),
        }
        if _lab.startswith("SEMIMETAL") or _lab.startswith("METAL"):
            out["validity"] = ("INVALID: band overlap, edge is formal (%s, gap=%s eV). "
                               "带重叠时 VBM/CBM 不是真正的带边，E1/D 只是形式上的极值点，"
                               "不应作为物理带边引用。" % (_lab, _gap))
            print("[WARN] 体系判别 %s -> band_edges.json validity 标为 INVALID" % _lab)
        else:
            out["validity"] = "OK"
    else:
        out["validity"] = "UNKNOWN: no discriminant.json"
except Exception as _de:
    print("[WARN] 读体系判别失败，validity 未写：%s" % _de)

with open(outdir + "/band_edges.json", "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)
print("[OK] band_edges.json 已生成（amset 权威带边 E1 + 真空对齐 E1_vac）")
'''
    be_out = out / "band_edges.json"
    be_cmd = ("%s && cd %s && python3 -c %s %s %s > band_edges.log 2>&1"
              % (AMSET_ENV_SRC, shlex.quote(str(dfm)),
                 shlex.quote(_be_code), shlex.quote(str(out)),
                 shlex.quote(_SKILL_REV)))
    _rc = subprocess.run(["bash", "-lc", be_cmd]).returncode
    if _rc != 0 or not be_out.is_file():
        # [PATCH-DPT-R5] 硬失败：band_edges.json 是 step8.2_dpt 正确 E1 的唯一
        # 来源，缺失时回退旧逻辑 = 回退「全带中位数」bug。宁可显式报错。
        _tail = ""
        _log = dfm / "band_edges.log"
        if _log.is_file():
            _tail = _log.read_text(errors="ignore")[-400:]
        sys.exit("[ERROR] band_edges.json 生成失败（step8.2_dpt 将无法取到带边 E1，"
                 "请勿回退旧逻辑）。看 %s/band_edges.log：\n%s"
                 % (dfm, _tail))
    print("[DONE] %s：band_edges.json 已生成" % OUTDIR_NAME)

    # [2026-09-16] 把**实际使用的**形变势构型口径落进 band_edges.json：settings 注释、
    # 8.3 汇总表、provenance 都从这里取。回退到 clamped 时这里就是证据（配合上面的 WARN）。
    try:
        import json as _json
        _be = _json.loads(Path(be_out).read_text())
        _be["deform_geometry"] = _def_geom
        _be["deform_geometry_folders"] = _geom_map
        Path(be_out).write_text(_json.dumps(_be, indent=2, ensure_ascii=False) + "\n")
        print("[OK] band_edges.json 记录 deform_geometry=%s（%d 个构型）"
              % (_def_geom, len(_geom_map)))
    except Exception as _e:                                    # noqa: BLE001
        print("[WARN] 写 deform_geometry 失败（不影响计算）：%s" % _e)

    # ---- patch_deform_ref：真空参考口径的 deformation_vac.h5 ----
    ref = _resolve_deform_ref(cwd)
    _vac_state = "core（未生成 deformation_vac.h5）"
    print("[..] DEFORM_REF=%s（本体系 DIM=%s）-> %s 参考口径"
          % (DEFORM_REF, _read_dim(cwd) or "未知", ref))
    if ref == "vacuum":
        # 做法（用户 2026-09-16 校对 AMSET 0.4.19 源码后指定）：**直接替换参考能级**，
        # 而不是事后往 h5 上加 Δ。理由：potentials.py 里
        #     diff = E_bulk - E_deform - (ref_bulk - ref_deform)
        #     energy_diff = np.abs(diff / strain)          <- 取了绝对值
        # 所以 (a) h5 里的 D 不带符号，往它上面加带符号 Δ 只在 D_core>0 时才对；
        # (b) D_vac - D_core = d(E_core - E_vac)/dε 需要**两个**应变导数，
        # 只用 -dE_vac/dε 缺一项。改成把 parse_calculation 的 reference 换成各构型
        # LOCPOT 的真空能级，再调用 AMSET 自己的 calculate_deformation_potentials，
        # 这样符号、分量平均、噪声过滤与 AMSET 原版完全一致。
        _vac_code = r'''
import json
import sys
from pathlib import Path

import numpy as np
from pymatgen.io.vasp.outputs import Locpot

import amset.deformation.io as _io
from amset.constants import defaults
from amset.deformation.io import write_deformation_potentials
from amset.deformation.potentials import (
    calculate_deformation_potentials, extract_bands,
    get_strain_mapping, get_symmetrized_strain_mapping, strain_coverage_ok)
from amset.electronic_structure.common import get_ibands
from amset.electronic_structure.symmetry import expand_bandstructure
# bz_coverage_ok 在 tools/deformation 里，不在 deformation/potentials
#（0.4.19 与 0.5.1 都是；2026-09-16 MoS2 实测过 ImportError）。
from amset.tools.deformation import bz_coverage_ok, check_calculation

dfm, be_json, h5_out = sys.argv[1], sys.argv[2], sys.argv[3]
VAC_HALF, VAC_FLAT = 0.25, 0.001        # 与 band_edges 的真空窗口/平坦判据一致


def _vac_window(structure, axis=2, half_ratio=VAC_HALF):
    """沿真空轴选最大分数间隙，返回中心 ± half_ratio 窗口（跨边界用 4 元组）。"""
    fracs = np.array([s.frac_coords for s in structure])
    f = np.sort(fracs[:, axis] % 1.0)
    gaps = np.diff(f).tolist() + [f[0] + 1.0 - f[-1]]
    imax = int(np.argmax(gaps))
    lo = f[imax]
    hi = f[imax + 1] if imax < len(f) - 1 else f[0] + 1.0
    c = (lo + hi) / 2.0
    w = (hi - lo) * half_ratio
    a, b = c - w, c + w
    if b > 1.0:
        b -= 1.0
    if a < 0.0:
        a += 1.0
    if a > b:
        return a, 1.0, 0.0, b
    return a, b, None, None


def _vacuum_level(folder, structure):
    """该构型 LOCPOT 的真空区平面平均静电势（eV）。平坦度不达标直接报错。"""
    lp = Locpot.from_file(str(Path(folder) / "LOCPOT"))
    vz = np.asarray(lp.get_average_along_axis(2))
    n = len(vz)
    win = _vac_window(structure)
    if win[2] is None:
        a = int(np.clip(np.floor(win[0] * n), 0, n - 1))
        b = int(np.clip(np.ceil(win[1] * n), a + 1, n))
        seg = vz[a:b]
    else:
        a1 = int(np.clip(np.floor(win[0] * n), 0, n - 1))
        b2 = int(np.clip(np.ceil(win[3] * n), 1, n))
        seg = np.concatenate([vz[a1:n], vz[0:b2]])
    if seg.std() > VAC_FLAT:
        raise RuntimeError("真空区不平坦（std=%.4f eV > %.4f）：%s"
                           % (seg.std(), VAC_FLAT, folder))
    return float(seg.mean())


_orig_parse = _io.parse_calculation


def parse_calc_vac(folder, zero_weighted_kpoints=None):
    """AMSET 原版 parse_calculation + 把 reference 换成该构型的真空能级。"""
    if zero_weighted_kpoints is None:
        calc = _orig_parse(folder)
    else:
        calc = _orig_parse(folder, zero_weighted_kpoints=zero_weighted_kpoints)
    st = calc["bandstructure"].structure
    calc["reference"] = _vacuum_level(folder, st)
    calc["reference_kind"] = "vacuum"
    return calc


zwk = defaults["zero_weighted_kpoints"]
symprec = defaults["symprec"]
# symprec_deformation：constants.defaults 里**没有**这个键（0.4.19 与 0.5.1 都
# KeyError，2026-09-16 实测）。AMSET 自己的默认是 symprec/100
# （amset/deformation/potentials.py: get_symmetrized_strain_mapping 的形参默认值）；
# 运行目录若有 settings.yaml 且写了该键，则用项目值。
symprec_def = None
try:
    import yaml as _yaml
    _sp = Path("settings.yaml")
    if _sp.is_file():
        symprec_def = (_yaml.safe_load(_sp.read_text()) or {}).get("symprec_deformation")
except Exception as _e:                                   # noqa: BLE001
    print("[WARN] 读 settings.yaml 的 symprec_deformation 失败（用默认）：%s" % _e)
if symprec_def is None:
    symprec_def = float(symprec) / 100.0
print("[..] symprec=%s symprec_deformation=%s" % (symprec, symprec_def))

print("[..] 读未形变构型（参考 = 真空能级）")
bulk = parse_calc_vac("undeformed", zwk)
bulk_struct = bulk["bandstructure"].structure
print("     vacuum ref = %.4f eV" % bulk["reference"])

calcs = []
for d in sorted(Path(dfm).glob("deform-*")):
    # 与 core 口径 h5 保持一致：形变势读**离子弛豫**构型（有 ionrelax/ 就用）
    dd = d / "ionrelax" if (d / "ionrelax" / "vasprun.xml").is_file() else d
    c = parse_calc_vac(str(dd), zwk)
    print("     %-22s vacuum ref = %.4f eV" % (str(dd.relative_to(dfm)), c["reference"]))
    c = check_calculation(bulk, c)
    if c is not False:
        calcs.append(c)
if not calcs:
    raise RuntimeError("没有可用的形变构型")

sm = get_strain_mapping(bulk_struct, calcs)
bulk["bandstructure"] = expand_bandstructure(bulk["bandstructure"], symprec=symprec)
sm = get_symmetrized_strain_mapping(bulk_struct, sm, symprec=symprec,
                                    symprec_deformation=symprec_def)
if not strain_coverage_ok(list(sm.keys())):
    raise RuntimeError("应变没有覆盖完整张量")
if not bz_coverage_ok([bulk["bandstructure"]] + [c["bandstructure"] for c in sm.values()]):
    raise RuntimeError("k 网格没有覆盖整个 BZ")

print("[..] 用 AMSET 原版 calculate_deformation_potentials 计算（参考=真空）")
dp = calculate_deformation_potentials(bulk, sm)
ibands = get_ibands(defaults["energy_cutoff"], bulk["bandstructure"])
dp = extract_bands(dp, ibands)
kpts = np.array([k.frac_coords for k in bulk["bandstructure"].kpoints])
write_deformation_potentials(dp, kpts, bulk_struct, filename=h5_out)
print("[OK] 写出 %s（真空参考，AMSET 原版算法与符号约定）" % h5_out)

# 自检：新 h5（真空口径）与 core 版 h5 在带边的 D 必须**不同**——相同就说明参考
# 替换没生效（静默退回 core 口径，是本补丁最危险的失败模式）。两者都用 band_edges
# 里那套 (band,k) 索引，读的是同一批形变目录，可比。
# 另一条更强的内部自检：D_vac - D_core = d(E_core - E_vac)/dε，**与能带无关**，
# 所以电子与空穴的差值必须相同。旧版（clamped 且空穴芯势 D 为负）这条不成立：
# 空穴 1.435-0.263=1.172 vs 电子 6.595-4.897=1.698 —— AMSET 对 D 取了绝对值，
# 把空穴芯势口径的**负号**吃掉了。这正说明"事后往 h5 上加 Δ"的做法对空穴会错，
# 而"替换参考能级"（本实现）天然避开该坑。ionrelax 口径下两者都是 2.087（一致）。
import h5py
be = json.load(open(be_json))
core_h5 = str(Path(dfm) / "deformation.h5")
with h5py.File(h5_out, "r") as fh:
    up = np.array(fh["deformation_potentials_up"])
with h5py.File(core_h5, "r") as fh:
    up_core = np.array(fh["deformation_potentials_up"])
print("[check] 带边 D_xx：core h5 vs 真空 h5（同一批形变目录，ionrelax 优先）")
_any_diff = False
_shifts = {}
for car in ("electron", "hole"):
    hits = (be.get(car) or {}).get("hits") or []
    if not hits:
        continue
    b, k = int(hits[0]["band_h5"]), int(hits[0]["k_h5"])
    dv, dc = float(up[b, k, 0, 0]), float(up_core[b, k, 0, 0])
    _shifts[car] = dv - dc
    if abs(dv - dc) > 1e-6:
        _any_diff = True
    print("   %-8s D_xx: core=%8.4f  真空=%8.4f  比=%.2f"
          % (car, dc, dv, (dv / dc) if dc else float("nan")))
    # 两份 h5 现在都读 ionrelax，与 band_edges.json 同一套数据 → 应该对得上
    # （带边 k/带索引一致时差几个百分点以内；差很多说明哪一边读错了构型）。
    _be = be[car].get("E1_vac_xx_eV")
    if _be:
        print("            [对照] band_edges.json（同口径 ionrelax）E1_vac=%s → 偏差 %+.1f%%"
              % (_be, 100.0 * (dv - float(_be)) / float(_be)))
if not _any_diff:
    raise RuntimeError("真空口径 h5 与 core h5 完全相同 —— 参考能级替换没生效"
                       "（检查 parse_calculation 的 reference 是否被改掉）")
# D_vac - D_core 与能带无关（= d(E_core-E_vac)/dε）：两种载流子必须给同一个值。
# 不成立时通常是某一侧的 D 实际为负、被 AMSET 的 abs() 吃掉符号（口径敏感），
# 所以只告警不中止；但要写进记录，便于回看。
if len(_shifts) == 2:
    _a, _b = _shifts["electron"], _shifts["hole"]
    _rel = abs(_a - _b) / max(abs(_a), abs(_b), 1e-9)
    print("   [check] D_vac-D_core（应与能带无关）：electron %+.4f  hole %+.4f"
          "  相对差 %.2f%%" % (_a, _b, 100.0 * _rel))
    if _rel > 0.05:
        print("   [WARN] 两种载流子的参考能级偏移不一致（>5%%）—— 多半是某一侧"
              "形变势为负被 abs() 吃了符号，或两个 h5 用了不同构型；"
              "核对 deform 目录口径后再用。")
'''
        vac_h5 = out / "deformation_vac.h5"
        # 子进程输出必须落到 deform_vac.log：autozt 只回显 gen 的**最后一行**，
        # 不落盘就等于失败无迹可查（2026-09-16 MoS2 实测：这里 rc!=0 只打印了一行
        # WARN，被 autozt 吞掉，只能看到最后那句 "deformation.h5 已生成（core 口径）"，
        # 根本不知道真空口径没生成）。
        _vlog = dfm / "deform_vac.log"
        cmd = ("%s && cd %s && python3 -c %s %s %s %s >> %s 2>&1"
               % (AMSET_ENV_SRC, shlex.quote(str(dfm)), shlex.quote(_vac_code),
                  shlex.quote(str(dfm)), shlex.quote(str(be_out)),
                  shlex.quote(str(vac_h5)), shlex.quote(str(_vlog))))
        rc = subprocess.run(["bash", "-lc", cmd]).returncode
        if rc != 0 or not vac_h5.is_file():
            _tail = (_vlog.read_text(errors="ignore")[-600:]
                     if _vlog.is_file() else "(子进程没有留下日志)")
            _vac_state = "core（deformation_vac.h5 生成失败 rc=%d）" % rc
            print("[WARN] deformation_vac.h5 生成失败（rc=%d）-> 本步仍按 core 口径继续；"
                  "把 LOCPOT 留在 step7_deform 各子目录里再重跑本步即可。"
                  "看 %s：\n%s" % (rc, _vlog, _tail))
        else:
            _vac_state = "vacuum（deformation_vac.h5 已生成）"
            print("[DONE] %s：deformation_vac.h5（真空参考口径）已生成；"
                  "deformation.h5 保留 core 口径作对比" % OUTDIR_NAME)
    else:
        print("[..] core 口径：不生成 deformation_vac.h5（三维默认；"
              "AMSET 官方口径就是芯态对齐）")
    if dst.exists():
        dst.unlink()
    os.replace(str(h5), str(dst))
    # ★ 最后一行是 autozt 唯一回显的一行，必须把"真空口径到底有没有生成"说清楚。
    print("[DONE] %s：deformation.h5 已生成（core 口径）；DIM=%s DEFORM_REF=%s -> %s"
          % (OUTDIR_NAME, _read_dim(cwd) or "未知", DEFORM_REF, _vac_state))


if __name__ == "__main__":
    main()
