#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lifetime_common.py -- 声子寿命/线宽分析公共库（供 phonon-lifetime 技能各步 import）。

物理约定（与 phono3py 官方后处理工具一致，勿改）
------------------------------------------------
phono3py 的 kappa-m*.hdf5 / gamma-m*.hdf5 里的 gamma 是声子自能的虚部
（linewidth），单位 THz（普通频率）。其内部把它当作 2*pi*gamma 的角频线宽使用，
因此寿命为::

    tau[ps] = 1 / (2 * 2 * pi * gamma_total[THz])

其中 gamma_total = gamma_anh + gamma_isotope (+ 边界散射 gamma_bd)，
phono3py 的 kappa 用的就是 gamma + gamma_isotope（conductivity/base.py:
_get_main_diagonal）；只取 gamma 会把寿命算长。即 tau = 1/(4*pi*gamma_total)。
1 THz = 1/ps，故 tau 直接以 ps 计。

统计量必须按不可约网格权重加权
------------------------------
kappa-m*.hdf5 里的 qpoint 只是不可约点，weight 是它在全网格里的星数。
若按点等权统计，高对称点被过度代表（Si 15x15x15：120 个不可约点，权重 1~48），
中位数/P10/P90/过阻尼比例都会静默偏掉。本模块所有汇总统计都支持 weights。

寿命与“振动周期”的比较（论文 Kim et al. Nature 2021 图 2e 的判据）
------------------------------------------------------------------
声子周期 T_period[ps] = 1/f[THz]；无量纲 omega*tau = 2*pi*f*tau。
omega*tau ~ 1（或 tau ~ T_period）表示强散射/近过阻尼，声子在完成一个振荡
前就失去了相干性（Ioffe-Regel 判据）。
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np

# phono3py gamma[THz] -> tau[ps] 的转换因子:  tau = 1/(TAU_FACTOR*gamma)
TAU_FACTOR = 2.0 * 2.0 * np.pi

# 判据默认阈值
DEFAULT_IMAG_THR = 0.10        # THz，低于此频率视为数值虚频/声学近零模，排除
# 保留常量仅供外部引用；0.2 起 acoustic_mask 默认不再用频率窗（见其 docstring）
DEFAULT_ACOUSTIC_FMAX = 2.0    # THz

# phono3py 的 --boundary-mfp/--bmfp 与 hdf5 里的 boundary_mfp 单位都是**微米**
# （help 原文 "Boundary mean free path in micrometer"；默认 1e6 µm = 10 mm，即“关”），
# 而群速度是 THz·Å（1 THz·Å = 100 m/s），所以算 gamma_bd 时 L 要乘 1e4 换成 Å。
MICRON_TO_ANGSTROM = 1.0e4


# ---------------------------------------------------------------------------
# 数据集定位
# ---------------------------------------------------------------------------
# 上游技能 -> κ / fc 产物所在步骤目录（各技能步骤编号不同：
# kl-dft-cpu 是 step6_kappa / step5_fc，kl-mlff-* 是 step4_kappa / step3_fc）。
# 0.2 之前写死 step6_kappa，接 kl-mlff-* 时直接报“找不到数据源”。
_UPSTREAM_STEPS = {
    "kl-dft-cpu": {"kappa": ("step6_kappa",), "fc": ("step5_fc",)},
    "kl-mlff-cpu": {"kappa": ("step4_kappa",), "fc": ("step3_fc",)},
    "kl-mlff-gpu": {"kappa": ("step4_kappa",), "fc": ("step3_fc",)},
    "kl-mlff-3090": {"kappa": ("step4_kappa",), "fc": ("step3_fc",)},
}
_KAPPA_STEP_NAMES = ("step6_kappa", "step4_kappa", "step5_kappa")
_FC_STEP_NAMES = ("step5_fc", "step3_fc")


def upstream_steps(skill, step=None, kind="kappa"):
    """按上游技能给出候选步骤目录名：显式 step 优先，再技能映射，最后通用兜底。"""
    known = _UPSTREAM_STEPS.get(str(skill or ""), {})
    common = _KAPPA_STEP_NAMES if kind == "kappa" else _FC_STEP_NAMES
    out = []
    for s in (step,) + tuple(known.get(kind, ())) + tuple(common):
        if s and str(s) not in out:
            out.append(str(s))
    return out


def _iter_candidate_dirs(cwd, skill, step):
    """给出可能存放上游技能产物的目录（由近及远）。

    兼容两种真实布局：
      * 远端：<matdir>/<skill>/<step>/...
      * 本地拉回副本：<matdir>/<skill>/result/<step>/...
    以及技能步目录 cwd 位于 <matdir>/<skill>/<step> 或 <matdir>/<step> 的情况。
    """
    cwd = Path(cwd).resolve()
    seen = []
    for base in (cwd, cwd.parent, cwd.parent.parent, cwd.parent.parent.parent):
        cands = [
            base / skill / step,
            base / skill / "result" / step,
            base / skill,
            base / "result" / step,
            base,
        ]
        for c in cands:
            c = c.resolve()
            if c not in seen:
                seen.append(c)
    return seen


def _hdf5_info(path):
    """读一个 hdf5 的元信息（打不开/不是 phono3py 文件返回 has_gamma=False）。"""
    info = {"has_gamma": False, "n_gp": 0, "n_band": 0, "mesh": [], "has_qw": False,
            "gamma_ndim": 0, "has_mode_kappa": False, "has_boundary": False,
            "boundary_mfp": None, "datasets": []}
    try:
        import h5py

        with h5py.File(path, "r") as f:
            info["datasets"] = sorted(f.keys())
            info["has_gamma"] = ("gamma" in f) and ("frequency" in f)
            if "frequency" in f:
                shp = f["frequency"].shape
                if len(shp) == 2:
                    info["n_gp"], info["n_band"] = int(shp[0]), int(shp[1])
            info["gamma_ndim"] = int(f["gamma"].ndim) if "gamma" in f else 0
            info["has_qw"] = ((("qpoint" in f) or ("qpoints" in f))
                              and (("weight" in f) or ("weights" in f)))
            info["has_mode_kappa"] = "mode_kappa" in f
            if "mesh" in f:
                info["mesh"] = [int(x) for x in np.asarray(f["mesh"][()]).ravel().tolist()]
            if "boundary_mfp" in f:
                info["has_boundary"] = True
                try:
                    info["boundary_mfp"] = float(np.asarray(f["boundary_mfp"][()]).ravel()[0])
                except Exception:
                    info["boundary_mfp"] = None
    except Exception:
        pass
    return info


def hdf5_has_gamma(path):
    return bool(_hdf5_info(path)["has_gamma"])


# 展宽(smearing)法写的 kappa-m*-s<σ>.hdf5；四面体法（默认 --br）没有 -s 后缀
_SMEARING_RX = re.compile(r"-s[0-9]")

_NAME_PENALTY = (
    (re.compile(r"-g[0-9]+\.hdf5$"), 4),     # --write-gamma --gp 的单点文件
    (re.compile(r"mfp"), 2),                  # 累积 mfp / kappa-mfp
    (re.compile(r"(-|_)gp(-|_|\.|$)"), 1),   # 单网格点输出
    (_SMEARING_RX, 1),                        # 展宽法输出 -> 同 mesh 优先四面体法
)


def _is_smearing_file(path):
    """该文件是否来自 phono3py 展宽(smearing)法（四面体法无 -s 后缀）。"""
    return bool(_SMEARING_RX.search(Path(path).name.lower()))


def _file_penalty(path):
    name = Path(path).name.lower()
    return sum(pen for rx, pen in _NAME_PENALTY if rx.search(name))


def select_kappa_hdf5(paths, log=None, note=""):
    """从候选 hdf5 里选主文件，返回 (Path|None, [被忽略的候选 Path])。

    0.1 用 sorted()[-1]（字符串序）选文件，kappa-m999 会盖过 kappa-m151515。
    这里改成按物理量排序：
      1) 只保留含 gamma+frequency 的；
      2) 优先多 q 点主文件（gamma 三维、n_qpoint>1），单点/多套网格分别降级；
      3) 再比 mesh 乘积、q 点数、是否带 qpoint/weight/mode_kappa；
      4) 文件名里带 -g0 / mfp / gp 的额外降级。
    有多个候选时通过 log 输出“选了谁、忽略了谁”，绝不静默换网格。
    """
    cands = []
    for p in dict.fromkeys(Path(x) for x in paths):
        info = _hdf5_info(p)
        if info["has_gamma"]:
            cands.append((Path(p), info))
    if not cands:
        return None, []

    def key(item):
        p, info = item
        multi = 1 if (info["gamma_ndim"] >= 3 and info["n_gp"] > 1) else 0
        mesh_prod = int(np.prod(info["mesh"])) if info["mesh"] else 0
        return (multi, -_file_penalty(p), mesh_prod, info["n_gp"],
                1 if info["has_qw"] else 0, 1 if info["has_mode_kappa"] else 0)

    cands.sort(key=key, reverse=True)
    best, ignored = cands[0], cands[1:]
    if log and _is_smearing_file(best[0]):
        log("[warn] %s选中的 %s 是展宽(smearing)法输出；同 mesh 若还有四面体法主文件，"
            "应优先用后者（默认 --br 是四面体法）" % (note, best[0].name))
    if ignored and log:
        log("[warn] %s发现 %d 个含 gamma 的 hdf5，选 %s (mesh=%s, n_qpoint=%d)；忽略 %s"
            % (note, len(cands), best[0].name,
               "x".join(str(x) for x in best[1]["mesh"]) or "?",
               best[1]["n_gp"],
               ", ".join("%s(n_qpoint=%d)" % (q.name, i["n_gp"]) for q, i in ignored[:8])))
    return best[0], [q for q, _ in ignored]


def is_main_hdf5(path):
    """该 hdf5 是否是“多 q 点主文件”（而非 --write-gamma --gp 的单点/单条输出）。"""
    info = _hdf5_info(path)
    return bool(info["has_gamma"] and info["gamma_ndim"] >= 3 and info["n_gp"] > 1)


def find_kappa_hdf5(cwd, skill="kl-dft-cpu", step="step6_kappa", log=None):
    """在材料目录树里找 phono3py 的 kappa/gamma hdf5（含 gamma 数据集）。

    逐个候选目录找；同一目录里用 select_kappa_hdf5 选主文件（不按字符串序）。
    返回 Path 或 None。
    """
    fallback = None
    seen = set()
    for st in upstream_steps(skill, step, "kappa"):
        for d in _iter_candidate_dirs(Path(cwd), skill, st):
            if not d.is_dir() or str(d) in seen:
                continue
            seen.add(str(d))
            files = []
            for pat in ("kappa-m*.hdf5", "gamma-m*.hdf5", "kappa*.hdf5", "*.hdf5"):
                files.extend(sorted(d.glob(pat)))
            best, _ = select_kappa_hdf5(files, log=log, note="[%s] " % d)
            if best is None:
                continue
            if is_main_hdf5(best):
                return best
            # 只有 --gp 单点文件时先记住，继续看有没有别处的多 q 点主文件
            if fallback is None:
                fallback = best
    if fallback is not None:
        return fallback
    # 兜底：在上游技能目录下递归找主文件（步骤目录名与预期不符时）
    for base in (Path(cwd), Path(cwd).parent, Path(cwd).parent.parent):
        sk = base / skill
        if not sk.is_dir():
            continue
        hits = sorted(sk.rglob("kappa-m*.hdf5"))
        if hits:
            best, _ = select_kappa_hdf5(hits, log=log, note="[rglob %s] " % sk)
            if best is not None:
                return best
    return None


def find_fc_dataset(cwd, skill="kl-dft-cpu", step="step6_kappa"):
    """找 phono3py 力常数数据集（fc2.hdf5 + fc3.hdf5 + phono3py_disp.yaml）。

    步骤目录名同样按上游技能展开（kl-mlff-* 的 fc 在 step3_fc）。
    返回 (dir, fc2, fc3, disp_yaml) 或 None。
    """
    for st in upstream_steps(skill, step, "fc"):
        for d in _iter_candidate_dirs(Path(cwd), skill, st):
            if not d.is_dir():
                continue
            for sub in (".", "phono3py"):
                base = d / sub
                fc2 = base / "fc2.hdf5"
                fc3 = base / "fc3.hdf5"
                disp = base / "phono3py_disp.yaml"
                if fc2.is_file() and fc3.is_file() and disp.is_file():
                    return base, fc2, fc3, disp
    return None


# ---------------------------------------------------------------------------
# 对称操作（沿 q 方向取样要用；缺 yaml/pyyaml/spglib 时优雅降级）
# ---------------------------------------------------------------------------
_YAML_CANDIDATES = ("phono3py.yaml", "phono3py_disp.yaml", "phono3py_params.yaml")


def find_phonopy_yaml(h5path=None, cwd=None, skill="kl-dft-cpu", step="step6_kappa"):
    """找 hdf5 旁边的 phono3py/phonopy yaml（含 primitive_cell，可求对称操作）。"""
    dirs = []
    if h5path is not None:
        dirs.append(Path(h5path).resolve().parent)
    if cwd is not None:
        cwd = Path(cwd).resolve()
        for st in upstream_steps(skill, step, "kappa"):
            dirs.extend([d for d in _iter_candidate_dirs(cwd, skill, st) if d.is_dir()])
        dirs.extend([cwd / "phono3py", cwd.parent / "phono3py"])
    seen, uniq = set(), []
    for d in dirs:
        if str(d) not in seen:
            seen.add(str(d))
            uniq.append(d)
    for d in uniq:
        for name in _YAML_CANDIDATES:
            p = d / name
            if p.is_file():
                return p
    return None


def symmetry_rotations_from_yaml(path):
    """从 yaml 的 primitive_cell 用 spglib 求分数坐标下的对称操作 (n,3,3)。

    缺 pyyaml 或 spglib、或解析失败时返回 None（调用方据此降级）。
    """
    try:
        import yaml
        import spglib
    except Exception:
        return None
    try:
        doc = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
        pc = doc.get("primitive_cell") or {}
        lat = np.asarray(pc.get("lattice"), dtype=float)
        points = pc.get("points") or []
        pos = np.asarray([p["coordinates"] for p in points], dtype=float)
        if lat.shape != (3, 3) or pos.ndim != 2 or pos.shape[0] == 0:
            return None
        syms = [str(p.get("symbol")) for p in points]
        uniq = {s: i + 1 for i, s in enumerate(sorted(set(syms)))}
        nums = [uniq[s] for s in syms]
        sym = spglib.get_symmetry((lat, pos, nums))
        if not sym:
            return None
        rot = np.asarray(sym["rotations"], dtype=float)
        return rot if rot.ndim == 3 and rot.shape[1:] == (3, 3) else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# 读取与转换
# ---------------------------------------------------------------------------
def read_kappa_hdf5(path):
    """读 phono3py kappa/gamma hdf5，返回 dict（numpy 数组）。

    关键字段：
      temperatures (nT,), frequency (n_gp, n_band) [THz], gamma (nT, n_gp, n_band) [THz],
      qpoints (n_gp, 3) [fractional], weights (n_gp,), mesh (3,),
      group_velocity (n_gp, n_band, 3) [THz·Å], heat_capacity (nT, n_gp, n_band) [eV/K],
      kappa (nT, 6) [W/mK], kappa_unit_conversion (scalar)
    """
    import h5py

    out = {}
    with h5py.File(path, "r") as f:
        for key in ("temperature", "frequency", "gamma", "gamma_isotope", "qpoint", "qpoints",
                    "weight", "weights", "mesh", "group_velocity", "gv_by_gv", "heat_capacity",
                    "kappa", "mode_kappa", "grid_point", "kappa_unit_conversion",
                    "boundary_mfp", "version"):
            if key in f:
                out[key] = f[key][()]
        # 不同 phono3py 版本用单/复数两种拼写
        for plural, singular in (("qpoints", "qpoint"), ("weights", "weight")):
            if plural in out and singular not in out:
                out[singular] = out.pop(plural)
    out["_file"] = str(path)
    if "temperature" in out:
        out["temperatures"] = np.asarray(out.pop("temperature"), dtype=float)
    if "qpoint" in out:
        out["qpoints"] = np.asarray(out.pop("qpoint"), dtype=float)
    if "weight" in out:
        out["weights"] = np.asarray(out.pop("weight"), dtype=float)
    out["frequency"] = np.asarray(out["frequency"], dtype=float)
    out["gamma"] = np.asarray(out["gamma"], dtype=float)
    return out


def effective_gamma(gamma, gamma_isotope=None):
    """总和线宽 g_total[THz] = gamma_anh + gamma_isotope。

    这不是可选装饰：phono3py 计算 kappa 时，碰撞矩阵对角元用的是
    gamma + gamma_isotope（conductivity/base.py:_get_main_diagonal），
    所以 hdf5 里的 kappa 含同位素散射。若只取 gamma，寿命会偏长。
    gamma_isotope 为空或无正元素时原样返回 gamma。
    """
    g = np.asarray(gamma, dtype=float)
    if gamma_isotope is None:
        return g
    gi = np.asarray(gamma_isotope, dtype=float)
    if gi.size == 0 or not np.any(gi > 0):
        return g
    return g + gi


def boundary_gamma(group_velocity, boundary_mfp_um):
    """边界散射线宽 gamma_bd[THz] = |v|/(4*pi*L)，L 由**微米**换成 Å。

    与 phono3py conductivity/base.py:_get_boundary_scattering 同源：那里
    g_bd = |v|/(2*2*pi*L)，对应 tau_bd = L/|v|。上游若用 EXTRA_ARGS
    --boundary-mfp/--bmfp 打开边界散射，hdf5 里的 kappa 就含它，寿命必须一并计入
    才能和 kappa 对上。

    ⚠️ 单位：phono3py 的 boundary_mfp 是**微米**（help: "Boundary mean free path in
    micrometer"），群速度是 THz·Å，所以必须先 L[µm]*1e4 -> L[Å]。0.2 版本漏了这一步，
    gamma_bd 偏大 1e4 倍（0.1 µm 时算出 ~5 THz 的荒唐线宽），于是“计入边界散射”反而
    让 kappa 重建误差变大，自动判定就会静默地不计入——B8 等于没修。
    L<=0/非有限时返回 None。
    """
    if group_velocity is None or boundary_mfp_um is None:
        return None
    L_um = float(boundary_mfp_um)
    if not np.isfinite(L_um) or L_um <= 0:
        return None
    v = np.asarray(group_velocity, dtype=float)
    if v.ndim != 3:
        return None
    vnorm = np.sqrt((v ** 2).sum(axis=-1))
    return vnorm / (4.0 * np.pi * L_um * MICRON_TO_ANGSTROM)


def gamma_to_tau_ps(gamma_thz):
    """gamma[THz] -> tau[ps]（phono3py 约定，见模块 docstring）。gamma<=0 -> +inf。"""
    g = np.asarray(gamma_thz, dtype=float)
    safe = np.where(g > 0.0, g, 1.0)
    with np.errstate(divide="ignore"):
        tau = np.where(g > 0.0, 1.0 / (TAU_FACTOR * safe), np.inf)
    return tau


def tau_to_gamma_thz(tau_ps):
    tau_ps = np.asarray(tau_ps, dtype=float)
    return 1.0 / (TAU_FACTOR * tau_ps)


def mean_free_path_ang(group_velocity, gamma_thz):
    """mfp[Å] = |v|[THz·Å] * tau[ps]（= |v|/(2*2*pi*gamma)），与 phono3py kaccum 一致。"""
    gv = np.asarray(group_velocity, dtype=float)
    g = np.asarray(gamma_thz, dtype=float)
    vnorm = np.sqrt((gv ** 2).sum(axis=-1))
    tau = gamma_to_tau_ps(g)
    return np.where(np.isfinite(tau), vnorm * tau, 0.0)


def reconstruct_kappa(gv_by_gv, heat_capacity, gamma_total, kappa_unit_conversion, mesh=None):
    """按 RTA 逐模重建 kappa（自检用；不参与本技能产出）。

    mode_kappa = gv_by_gv * cv / (2*gamma_total) * kappa_unit_conversion
    kappa = mode_kappa.sum(axis=(0,1)) / prod(mesh)
    其中 gv_by_gv(n_gp,n_band,6) 已是 Voigt 六分量、且已带不可约网格权重。
    """
    gv2 = np.asarray(gv_by_gv, dtype=float)
    cv = np.asarray(heat_capacity, dtype=float)
    gt = np.asarray(gamma_total, dtype=float)
    if gv2.ndim != 3 or cv.ndim != 2 or gt.shape != cv.shape or gv2.shape[:2] != cv.shape:
        return None
    uc = float(kappa_unit_conversion)
    safe = np.where(gt > 0, gt, 1.0)
    mk = gv2 * cv[..., None] / (2.0 * safe[..., None]) * uc
    mk = np.where(gt[..., None] > 0, mk, 0.0)
    if mesh is None:
        n = gt.shape[0]
    elif np.isscalar(mesh):
        n = int(mesh)
    else:
        arr = np.asarray(mesh, dtype=float).ravel()
        n = int(np.prod(arr)) if arr.size else gt.shape[0]
    if n <= 0:
        n = gt.shape[0]
    return mk.sum(axis=(0, 1)) / float(n)


def kappa_reconstruction_error(data, gamma_total, t_indices=None):
    """把重建 kappa 与 hdf5 里存的 kappa 对比，返回 dict 或 None。

    自检项：能对上说明 gamma/权重/边界散射口径与上游一致；对不上就该怀疑
    统计口径（例如上游开了边界散射而我们没计入）。
    """
    need = ("gv_by_gv", "heat_capacity", "kappa", "kappa_unit_conversion")
    if any(k not in data for k in need):
        return None
    gt = np.asarray(gamma_total, dtype=float)
    if gt.ndim == 2:
        gt = gt[None, ...]
    kappa = np.asarray(data["kappa"], dtype=float)
    temps = np.asarray(data.get("temperatures", []), dtype=float)
    mesh = data.get("mesh")
    sel = list(range(gt.shape[0])) if t_indices is None else list(t_indices)
    per_T, errs = {}, []
    for ti in sel:
        if ti >= gt.shape[0] or ti >= kappa.shape[0] or ti >= len(temps):
            continue
        k = reconstruct_kappa(data["gv_by_gv"], np.asarray(data["heat_capacity"], dtype=float)[ti],
                              gt[ti], data["kappa_unit_conversion"], mesh=mesh)
        if k is None:
            return None
        ref = kappa[ti]
        # 以张量量级为分母（避免被 ~1e-14 的数值零分量放大相对误差）
        scale = max(float(np.max(np.abs(ref))), 0.0) if ref.size else 0.0
        denom = max(scale, 1e-12)
        e = float(np.max(np.abs(k - ref)) / denom)
        errs.append(e)
        per_T["%.4f" % float(temps[ti])] = e
    if not errs:
        return None
    return {"max_rel_err": max(errs), "n_temperature": len(errs), "rel_err_per_T": per_T}


def decide_boundary_gamma(data, gamma_total, g_bd, t_indices=None, sentinel_um=1e6):
    """决定要不要把边界散射线宽 gamma_bd 计入寿命（返回 (use, info)）。

    phono3py 的 hdf5 总是写 boundary_mfp，默认 1e6 µm 就是“关”（它自己的默认值，
    单位微米）。上游若真用 --boundary-mfp 打开了边界散射，hdf5 里的 kappa 就含它，
    此时计入 gamma_bd 会明显降低 kappa 重建误差；若只是默认占位，计入反而让误差变大。

    判定顺序：L >= sentinel_um 直接当“关”；否则能比对 kappa 就用重建误差决定，
    比不了再回退到阈值 L < sentinel_um。
    """
    info = {"method": "kappa_reconstruction", "err_no_boundary": None,
            "err_with_boundary": None, "boundary_mfp_um": None,
            "gamma_bd_median_THz": None}
    mfp = data.get("boundary_mfp")
    if mfp is not None:
        try:
            info["boundary_mfp_um"] = float(np.asarray(mfp).ravel()[0])
        except Exception:
            info["boundary_mfp_um"] = None
    if g_bd is not None:
        gb = np.asarray(g_bd, dtype=float)
        if gb.size:
            info["gamma_bd_median_THz"] = float(np.median(gb))
    if g_bd is None:
        info["method"] = "no_group_velocity"
        return False, info
    L = info["boundary_mfp_um"]
    if L is not None and np.isfinite(L) and L >= sentinel_um:
        info["method"] = "mfp_sentinel"
        return False, info
    err_no = kappa_reconstruction_error(data, gamma_total, t_indices)
    err_bd = kappa_reconstruction_error(data, gamma_total + g_bd, t_indices)
    if err_no and err_bd:
        info["err_no_boundary"] = err_no["max_rel_err"]
        info["err_with_boundary"] = err_bd["max_rel_err"]
        return bool(err_bd["max_rel_err"] < err_no["max_rel_err"]), info
    info["method"] = "mfp_threshold"
    return bool(L is not None and np.isfinite(L) and 0 < L < sentinel_um), info


# ---------------------------------------------------------------------------
# 温度 / 声学支选择
# ---------------------------------------------------------------------------
def pick_temperature_indices(temperatures, wanted=None, tol=1.0, log=None):
    """给定温度数组，返回 {T: index}。

    wanted 为空取全部；否则取最接近的索引（去重）。请求温度与实际可用温度
    相差超过 tol 时通过 log 输出警告（0.1 曾静默取最近温度，1000K 也能悄悄
    退化成 800K）。
    """
    temperatures = np.asarray(temperatures, dtype=float)
    if not wanted:
        return {float(t): i for i, t in enumerate(temperatures)}
    chosen = {}
    for t in wanted:
        i = int(np.argmin(np.abs(temperatures - float(t))))
        chosen[float(temperatures[i])] = i
        if log and abs(float(temperatures[i]) - float(t)) > tol:
            log("请求的温度 %.4g K 不在可用列表里，改用最近的 %.4g K" % (t, temperatures[i]))
    return chosen


def acoustic_mask(frequencies, gamma, acoustic_fmax=None, imag_thr=DEFAULT_IMAG_THR):
    """声学支掩码：默认取最低 3 条支（index 0/1/2），不加频率窗口。

    0.1 用“最低 3 支中位频率 x2”自动定频率上限，会把 X 点附近约 12 THz 的 LA
    归进光学组（Si 15x15x15 自动算出约 8.5 THz），声学/光学中位数因此失真。
    现在只有显式给出 acoustic_fmax 时才叠加频率窗（2D/含真空材料可用它剔除
    混在最低 3 支里的真空零模）。单原子胞（n_band<=3）没有光学支，光学组为空。
    返回 (mask, fmax_used)；fmax_used=None 表示未用频率窗。
    """
    f = np.asarray(frequencies, dtype=float)
    g = np.asarray(gamma, dtype=float)
    n_band = f.shape[1]
    n_ac = min(3, n_band)
    branch = np.zeros_like(f, dtype=bool)
    branch[:, :n_ac] = True
    good = (f > imag_thr) & (g > 0)
    if acoustic_fmax is None:
        return branch & good, None
    fmax = float(acoustic_fmax)
    return branch & good & (f <= fmax), fmax


# ---------------------------------------------------------------------------
# 加权统计
# ---------------------------------------------------------------------------
def _weights_for(weights, shape):
    """把 (n_gp,) 的不可约权重广播成 shape；不可用（None/形状不符/非正）返回 None。"""
    if weights is None:
        return None
    w = np.asarray(weights, dtype=float)
    if w.ndim == 1:
        if w.size != shape[0]:
            return None
        w = np.broadcast_to(w[:, None], shape).astype(float)
    elif w.shape != shape:
        return None
    if w.size == 0 or not np.all(np.isfinite(w)) or w.sum() <= 0:
        return None
    return w


_WEIGHT_EXPAND_CAP = 4_000_000


def _expand_by_weight(values, weights):
    """整数权重时把样本展开成等权样本（这样中位数/P10/P90 与“展开到全网格”
    的定义严格一致）；权重非整数或展开过大时返回 None。"""
    if weights is None:
        return None
    w = np.asarray(weights, dtype=float)
    wr = np.round(w)
    if not np.allclose(w, wr, atol=1e-9):
        return None
    total = int(wr.sum())
    if total <= 0 or total > _WEIGHT_EXPAND_CAP:
        return None
    return np.repeat(np.asarray(values, dtype=float), wr.astype(np.int64))


def _quantile(values, q, weights=None):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        rep = _expand_by_weight(v, w[np.isfinite(np.asarray(values, dtype=float))])
        if rep is not None and rep.size:
            return float(np.percentile(rep, q))
        order = np.argsort(v)
        vs, ws = v[order], np.asarray(weights, dtype=float)[np.isfinite(np.asarray(values, dtype=float))][order]
        cw = np.cumsum(ws)
        if cw[-1] > 0:
            idx = int(np.searchsorted(cw, q / 100.0 * cw[-1], side="left"))
            return float(vs[min(idx, vs.size - 1)])
    return float(np.percentile(v, q))


def _mean(values, weights=None):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        w = w[np.isfinite(np.asarray(values, dtype=float))]
        if w.size == v.size and w.sum() > 0:
            return float(np.average(v, weights=w))
    return float(np.mean(v))


def _stats(tau, omega_tau, weights=None):
    tau = np.asarray(tau, dtype=float)
    tau = tau[np.isfinite(tau)]
    if tau.size == 0:
        return {"n_modes": 0, "n_modes_weighted": 0.0}
    ot = np.asarray(omega_tau, dtype=float)
    ot = ot[np.isfinite(ot)]
    return {
        "n_modes": int(tau.size),
        "n_modes_weighted": (float(np.sum(weights)) if weights is not None else float(tau.size)),
        "tau_median_ps": _quantile(tau, 50, weights),
        "tau_mean_ps": _mean(tau, weights),
        "tau_p10_ps": _quantile(tau, 10, weights),
        "tau_p90_ps": _quantile(tau, 90, weights),
        "omega_tau_median": (_quantile(ot, 50, weights) if ot.size else None),
    }


# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
def summarize_temperature(frequencies, gamma, temperatures, t_index,
                          acoustic_fmax=None, imag_thr=DEFAULT_IMAG_THR,
                          omega_tau_threshold=1.0, gamma_isotope=None,
                          gamma_extra=None, weights=None):
    """对单个温度汇总：全体 / 声学 / 光学 分组的 tau、omega*tau 与过阻尼比例。

    参数
      gamma_isotope : (n_gp,n_band) 或 (nT,n_gp,n_band)；叠加到三声子线宽上。
      gamma_extra   : 额外的线宽项（如边界散射 gamma_bd），同样叠加。
      weights       : (n_gp,) 不可约网格权重；给出时所有统计量按权重加权。
    """
    f = np.asarray(frequencies, dtype=float)
    g_anh = np.asarray(gamma[t_index], dtype=float)
    g = effective_gamma(g_anh, gamma_isotope)
    extra = None
    if gamma_extra is not None:
        extra = np.asarray(gamma_extra, dtype=float)
        if extra.shape[-2:] == g.shape[-2:] and extra.ndim <= g.ndim:
            g = g + extra
        else:
            extra = None
    tau = gamma_to_tau_ps(g)
    tau_anh = gamma_to_tau_ps(g_anh)
    omega_tau = 2.0 * np.pi * f * tau
    isotope = bool(gamma_isotope is not None and np.any(np.asarray(gamma_isotope) > 0))
    w_all = _weights_for(weights, g.shape)

    good = (f > imag_thr) & (g > 0)
    ac_mask, fmax = acoustic_mask(f, g, acoustic_fmax, imag_thr)
    opt_mask = good & ~ac_mask

    def pack(mask):
        sel = mask & np.isfinite(tau) & np.isfinite(omega_tau)
        if not sel.any():
            return {"n_modes": 0, "n_modes_weighted": 0.0, "overdamped_frac": None,
                    "f_min_THz": None, "f_max_THz": None,
                    "tau_median_anharmonic_ps": None}
        w = w_all[sel] if w_all is not None else None
        d = _stats(tau[sel], omega_tau[sel], w)
        d["overdamped_frac"] = _mean((omega_tau[sel] < omega_tau_threshold).astype(float), w)
        d["f_min_THz"] = float(f[sel].min())
        d["f_max_THz"] = float(f[sel].max())
        d["tau_median_anharmonic_ps"] = _quantile(tau_anh[sel], 50, w)
        return d

    return {
        "T_K": float(temperatures[t_index]),
        "isotope_included": isotope,
        "boundary_included": bool(extra is not None and np.any(extra > 0)),
        "weights_applied": w_all is not None,
        "all": pack(good),
        "acoustic": pack(ac_mask),
        "optical": pack(opt_mask),
        "acoustic_fmax_THz": (float(fmax) if fmax is not None else None),
        "omega_tau_threshold": float(omega_tau_threshold),
    }


# ---------------------------------------------------------------------------
# 沿指定 q 方向取样（论文 Γ-A 那类“寿命 vs q”）
# ---------------------------------------------------------------------------
def select_direction_qpoints(qpoints, direction, n_points=21, tol=1e-5, rotations=None):
    """从（不可约）网格 q 点里挑出落在 Γ -> direction 线上的点，返回 [(t, gp), ...]。

    hdf5 存的是不可约点，等价方向的代表点未必几何上在这条线上（Si 15^3 的
    不可约楔子存 (x,0,0)，直接找 [001] 只会命中 Γ）。给了 rotations（分数坐标
    下的对称操作，见 symmetry_rotations_from_yaml）时，把每个不可约点用
    inv(R).T 映射一圈，命中目标线的像也算——τ 在对称操作下不变；另外把时间反演
    τ(q)=τ(-q) 的像（-imgs）也放进候选，没有反演中心的材料不可约楔里可能只存 -q。
    t 为该线上的约化坐标（0=Γ，0.5=该方向第一布里渊区边界），升序去重；同一个 t
    只保留离目标线最近的那个点（不是先遇到的那个）。

    ⚠️ 容差 tol 默认 1e-5（分数坐标）：网格点是精确有理数，31^3 的格点间距才
    1/31≈0.032，旧默认 0.06 会把紧贴目标线的离线点也当成线上的点（实测 31^3
    会整列取错）。方向与网格点不共线时会只剩 Γ，这是刻意的。
    n_points 只作上限（结果多于 n_points+1 个时按 t 均匀抽稀）。
    """
    q = np.asarray(qpoints, dtype=float)
    d = np.asarray(direction, dtype=float)
    nrm = np.linalg.norm(d)
    if nrm == 0 or q.size == 0:
        return []
    d = d / nrm
    qq = q - np.round(q)
    rinv_t = None
    if rotations is not None:
        rot = np.asarray(rotations, dtype=float)
        if rot.ndim == 3 and rot.shape[1:] == (3, 3) and rot.shape[0]:
            try:
                rinv_t = np.array([np.linalg.inv(r).T for r in rot])
            except np.linalg.LinAlgError:
                rinv_t = None
    tol = max(float(tol), 0.0)
    found = {}
    for i in range(qq.shape[0]):
        qi = qq[i]
        imgs = qi[None, :] if rinv_t is None else (rinv_t @ qi)
        imgs = imgs - np.round(imgs)
        imgs = np.concatenate([imgs, -imgs], axis=0)       # 时间反演 -q
        proj = imgs @ d
        dist = np.linalg.norm(imgs - np.outer(proj, d), axis=1)
        for k in np.where((dist <= tol) & (proj >= -tol) & (proj <= 0.5 + tol))[0]:
            tt = float(proj[k])
            if tt < 0.0:
                continue
            # 按实际 t 去重（不是按 n_points 量化后的格点），否则不在格点上的
            # t（如 1/15 的奇数倍）会在与邻居竞争时被丢掉；同一 t 取最近的
            key = round(tt, 6)
            prev = found.get(key)
            if prev is None or dist[k] < prev[2]:
                found[key] = (tt, int(i), float(dist[k]))
    picked = [(t, i) for t, i, _ in sorted(found.values())]
    if n_points and len(picked) > int(n_points) + 1:
        idx = np.linspace(0, len(picked) - 1, int(n_points) + 1).round().astype(int)
        picked = [picked[j] for j in dict.fromkeys(idx.tolist())]
    return picked


# ---------------------------------------------------------------------------
# 落盘
# ---------------------------------------------------------------------------
def write_mode_csv(path, temperatures, frequencies, gamma, qpoints, weights,
                   group_velocity=None, acoustic_fmax=None,
                   imag_thr=DEFAULT_IMAG_THR, t_indices=None, gamma_isotope=None,
                   gamma_extra=None):
    """逐模写 CSV。

    tau_ps 用总和线宽（gamma + gamma_isotope + gamma_extra），与 phono3py 的 kappa
    一致；tau_anh_ps 只用三声子 gamma。mfp 用逐模群速度（不是 q 点平均）。
    """
    f = np.asarray(frequencies, dtype=float)
    g_all = np.asarray(gamma, dtype=float)
    gi_all = None if gamma_isotope is None else np.asarray(gamma_isotope, dtype=float)
    ge_all = None if gamma_extra is None else np.asarray(gamma_extra, dtype=float)
    q = np.asarray(qpoints, dtype=float)
    w = np.asarray(weights, dtype=float)
    gv = None if group_velocity is None else np.asarray(group_velocity, dtype=float)
    n_gp, n_band = f.shape
    if t_indices is None:
        t_indices = range(len(temperatures))

    def _at(arr, ti):
        if arr is None:
            return None
        if arr.ndim == 2:
            return arr
        return arr[ti]

    with open(path, "w", newline="", encoding="utf-8") as fh:
        csvw = csv.writer(fh)
        csvw.writerow(["T_K", "gp", "qx", "qy", "qz", "weight", "band",
                       "f_THz", "gamma_THz", "gamma_iso_THz", "gamma_bd_THz",
                       "gamma_total_THz", "tau_ps", "tau_anh_ps", "omega_tau",
                       "period_ps", "v_THzA", "mfp_A", "acoustic"])
        for ti in t_indices:
            g_anh = g_all[ti]
            gi = _at(gi_all, ti)
            ge = _at(ge_all, ti)
            g = effective_gamma(g_anh, gi)
            if ge is not None and ge.shape == g.shape:
                g = g + ge
            else:
                ge = None
            tau = gamma_to_tau_ps(g)
            tau_anh = gamma_to_tau_ps(g_anh)
            omega_tau = 2.0 * np.pi * f * tau
            ac_mask, _ = acoustic_mask(f, g, acoustic_fmax, imag_thr)
            for gp in range(n_gp):
                vnorm_b = None
                if gv is not None and gv.ndim == 3:
                    vnorm_b = np.sqrt((gv[gp, :, :] ** 2).sum(axis=-1))
                for b in range(n_band):
                    period = 1.0 / f[gp, b] if f[gp, b] > 0 else 0.0
                    vb = float(vnorm_b[b]) if vnorm_b is not None else None
                    mfp = (vb * tau[gp, b]) if (vb is not None and g[gp, b] > 0) else 0.0
                    csvw.writerow([
                        "%.4f" % temperatures[ti], gp,
                        "%.6f" % q[gp, 0], "%.6f" % q[gp, 1], "%.6f" % q[gp, 2],
                        "%.6f" % w[gp], b,
                        "%.6f" % f[gp, b], "%.6e" % g_anh[gp, b],
                        "%.6e" % (float(gi[gp, b]) if gi is not None else 0.0),
                        "%.6e" % (float(ge[gp, b]) if ge is not None else 0.0),
                        "%.6e" % g[gp, b],
                        ("%.6g" % tau[gp, b]) if np.isfinite(tau[gp, b]) else "inf",
                        ("%.6g" % tau_anh[gp, b]) if np.isfinite(tau_anh[gp, b]) else "inf",
                        ("%.6g" % omega_tau[gp, b]) if np.isfinite(omega_tau[gp, b]) else "inf",
                        "%.6g" % period,
                        ("%.6f" % vb) if vb is not None else "",
                        "%.6g" % mfp,
                        int(bool(ac_mask[gp, b])),
                    ])


def dump_json(path, obj):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, default=_json_default)


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


__all__ = [
    "TAU_FACTOR", "DEFAULT_IMAG_THR", "DEFAULT_ACOUSTIC_FMAX",
    "find_kappa_hdf5", "find_fc_dataset", "find_phonopy_yaml", "hdf5_has_gamma",
    "select_kappa_hdf5", "is_main_hdf5", "upstream_steps", "symmetry_rotations_from_yaml",
    "read_kappa_hdf5", "effective_gamma", "boundary_gamma",
    "gamma_to_tau_ps", "tau_to_gamma_thz", "mean_free_path_ang",
    "reconstruct_kappa", "kappa_reconstruction_error", "decide_boundary_gamma",
    "pick_temperature_indices", "acoustic_mask", "summarize_temperature",
    "select_direction_qpoints", "write_mode_csv", "dump_json",
]
