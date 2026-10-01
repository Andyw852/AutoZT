#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""postprocess_intrinsic.py —— 从 AMSET 的 mesh.h5 严格重算"本征迁移率"。

背景（2026-09-22 用户方案）
-------------------------
S8.4_amset2d 一次运行里同时算 ADP / IMP / POP 三种机制，AMSET 把总散射率
（各机制求和）送进输运积分，得到的是"含电离杂质"的迁移率。要得到**严格本征**
（只含 ADP+POP、不含 IMP）的 μ、S、σ，正确做法不是把各机制的迁移率做
Matthiessen 倒数相加（那对张量/各向异性体系只是近似），而是：
    读取每机制每 k 点散射率 → 把 IMP 项置零 → 用 AMSET 自己的输运积分重算。

本脚本就做这件事：
    1. 读运行目录里的 settings.yaml / vasprun.xml / mesh_*.h5；
    2. 用 AMSET 自己的 Interpolator 重建 AmsetData（能带、速度、DOS、费米能级），
       这一步**不重算散射**（散射率直接来自 mesh.h5，是最贵的部分，跳过）；
    3. 把 mesh.h5 的（仅不可约 k 点的）每机制散射率展开回全网格；
    4. 把要剔除的机制（默认 IMP）整行置零；
    5. 调 amset.core.transport.solve_boltzman_transport_equation 重新积分，
       输出严格本征的 μ（overall + 保留机制单项）、S、σ 到 JSON。

等价性说明（重要）
----------------
· 各机制的散射率是彼此独立计算的（每个 scatterer 只用自己的 prefactor/factor），
  所以"含 IMP 的 3 机制运行里把 IMP 行置零" == "单独跑一次 [ADP, POP]"。
· 本脚本已在 0.4.19 上验证：用它重建+注入 mesh 散射率，**逐位复现**原运行的
  transport.json（max|diff| = 0.0）——见 --check-reproduce 与 tests.log。
· 因为 IMP 行被置零后，"单项 IMP 迁移率"会是 1/0=inf 并让 BoltzTraP2 的 SVD 不收敛，
  所以本脚本的 overall 用 separate_mobility=False，单项只对**保留的**机制逐一调用
  AMSET 的 _calculate_mobility；这样既满足"置零"，又能给出每机制数值。

AMSET 版本与已知坑
------------------
· 0.4.19 与 0.5.1 的 amset/core/data.py、amset/core/transport.py、amset/io.py
  **逐字节相同**；write_mesh / load_mesh / solve_boltzman_transport_equation
  的 API 在两版之间没有差异。差异只在 scattering/calculate.py 的重叠（umklapp
  g_diff）等处，与输运重积分无关。
· 两版的 amset/io.py:write_mesh 有一个**自旋极化 bug**：循环里把 key 变量改写，
  导致自旋向下的数据集被写成 "<name>_up_down"（如 energies_up_down），
  而 load_mesh 又把它误读成另一个顶层键 "<name>_up"。
  因此**不要用 amset.io.load_mesh** 读自旋极化体系的 mesh.h5；本脚本自带
  read_mesh_h5() 修正这个 key 解析，up/down 都能正确归位。

用法
----
    # 在含 settings.yaml + vasprun.xml + mesh_*.h5 的 S8.4 运行目录里：
    python postprocess_intrinsic.py                       # 默认剔除 IMP，写 intrinsic_transport.json
    python postprocess_intrinsic.py --drop IMP,PIE        # 剔除多个机制
    python postprocess_intrinsic.py --check-reproduce     # 额外校验能复现原 transport.json
    python postprocess_intrinsic.py --mesh mesh_47x47x3.h5 --out intrinsic.json

纯逻辑函数（split_mesh_key/expand_ir_to_full/zero_mechanisms/...）不依赖 amset/h5py，
可离线单测：python test_postprocess_intrinsic.py
"""
import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

__author__ = "AutoZT / DSH subagent"
__version__ = "2026-09-22-intrinsic-v1"

# AMSET 的 write_mesh 对自旋向下数据集的 bug 编码："<name>_up_down"。
# 顺序重要：先匹配最长的 "_up_down"，再 "_down"，最后 "_up"。
_SPIN_SUFFIXES = (("_up_down", "down"), ("_down", "down"), ("_up", "up"))

# 默认要从"本征"结果里剔除的机制（用户方案：IMP 置零）。
DEFAULT_DROP = ("IMP",)


# =====================================================================
# 纯逻辑（无 amset / h5py 依赖，可离线单测）
# =====================================================================
def split_mesh_key(name):
    """把一个 mesh.h5 原始数据集名拆成 (base_name, spin)。

    Args:
        name: h5 里的数据集名，如 "energies_up"、"scattering_rates_up_down"。

    Returns:
        (base, spin)：base 是去掉自旋后缀的名字；spin 是 "up"/"down"，无后缀时为 None。

    Notes:
        AMSET 0.4.19/0.5.1 的 write_mesh 有 key 改写 bug，自旋向下写成
        "<name>_up_down"，这里一并正确处理。
    """
    for suffix, spin in _SPIN_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)], spin
    return name, None


def norm_mech(name):
    """机制名归一化（大写、去空格）。"""
    return str(name).strip().upper()


def resolve_mechanism_indices(labels, drop):
    """返回 labels 中命中 drop 的下标列表（按 labels 顺序）。"""
    drop_set = {norm_mech(x) for x in drop}
    return [i for i, lab in enumerate(labels) if norm_mech(lab) in drop_set]


def zero_mechanisms(rates, labels, drop):
    """把指定机制的散射率整行置零。

    Args:
        rates: {spin: ndarray}，每个 ndarray 形如
               (n_labels, n_doping, n_temps, n_bands, n_kpoints)。
        labels: 机制名列表，长度 = rates[spin].shape[0]。
        drop: 要置零的机制名（可迭代）。

    Returns:
        (rates_zeroed, dropped_labels, kept_labels)：保持 labels 不变，
        仅把对应行拷贝后置 0；返回置零/保留的机制名列表（便于报告）。
    """
    idx = resolve_mechanism_indices(labels, drop)
    idx_set = set(idx)
    out = {}
    for spin, r in rates.items():
        r2 = np.array(r, copy=True)
        for i in idx:
            r2[i] = 0.0
        out[spin] = r2
    dropped = [labels[i] for i in idx]
    kept = [lab for i, lab in enumerate(labels) if i not in idx_set]
    return out, dropped, kept


def drop_mechanisms(rates, labels, drop):
    """把指定机制整行删除（等价于置零后的 overall，但避免 0 行导致的 1/0）。

    Returns:
        (rates_kept, kept_labels, dropped_labels)。
    """
    idx = set(resolve_mechanism_indices(labels, drop))
    keep_idx = [i for i in range(len(labels)) if i not in idx]
    out = {spin: np.asarray(r)[keep_idx] for spin, r in rates.items()}
    kept = [labels[i] for i in keep_idx]
    dropped = [labels[i] for i in sorted(idx)]
    return out, kept, dropped


def expand_ir_to_full(rates_ir, ir_to_full):
    """把不可约 k 点的散射率展开到全网格。

    mesh.h5 里散射率只存了 ir_kpoints_idx 对应的列（to_dict 里做 r[..., ir_idx]），
    ir_to_full 是"每个全网格 k 点 -> 其不可约代表下标"的映射（AMSET 的
    ir_to_full_kpoint_mapping），因此直接按最后一维取列即可。

    Args:
        rates_ir: (..., n_ir_kpoints)
        ir_to_full: (n_full_kpoints,) 整数映射

    Returns:
        (..., n_full_kpoints)
    """
    ir_to_full = np.asarray(ir_to_full)
    return np.asarray(rates_ir)[..., ir_to_full]


def inplane_average(tensor):
    """张量 x-y 面内平均 = (xx + yy) / 2（最后一维 3x3）。"""
    t = np.asarray(tensor, dtype=float)
    if t.shape[-2:] != (3, 3):
        raise ValueError("inplane_average 需要最后两维为 3x3，收到 %s" % (t.shape,))
    return 0.5 * (t[..., 0, 0] + t[..., 1, 1])


def to_jsonable(arr):
    """ndarray -> 可 json 序列化的嵌套 list（float 保精度）。"""
    return np.asarray(arr, dtype=float).tolist()


# =====================================================================
# mesh.h5 读取（修正 AMSET 的自旋 key bug）
# =====================================================================
def read_mesh_h5(path):
    """读 AMSET mesh.h5，返回与 amset.io.load_mesh 相同结构但修正自旋分组。

    返回 dict:
        标量/无自旋键 -> numpy 值（structure 为 pymatgen Structure）；
        自旋键 -> {Spin.up: ndarray, Spin.down: ndarray}。
    """
    import h5py
    from pymatgen.core.structure import Structure
    from pymatgen.electronic_structure.core import Spin

    spin_map = {"up": Spin.up, "down": Spin.down}

    def read_val(base, ds):
        if base == "structure":
            s = ds[()]
            s = s.decode() if isinstance(s, (bytes, np.bytes_)) else str(s)
            return Structure.from_str(s, fmt="json")
        if base == "scattering_labels":
            return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x)
                    for x in ds[()].astype("U")]
        if base == "vb_idx":
            d = ds[()]
            return None if d is False else d
        return ds[()]

    out = {}
    with h5py.File(path, "r") as f:
        for raw in f.keys():
            base, spin = split_mesh_key(raw)
            val = read_val(base, f[raw])
            if spin is None:
                out[base] = val
            else:
                out.setdefault(base, {})[spin_map[spin]] = val
    return out


# =====================================================================
# 用 AMSET 自己的代码重建 AmsetData（不重算散射）
# =====================================================================
def _settings_symprec(run_dir):
    """settings.yaml 的 symprec（AMSET 求不可约 k 点用的那个）；没写就用 AMSET 默认。"""
    try:
        import yaml
        s = yaml.safe_load((Path(run_dir) / "settings.yaml").read_text()) or {}
        if s.get("symprec"):
            return float(s["symprec"])
    except Exception:                                              # noqa: BLE001
        pass
    try:
        from amset.constants import defaults
        return float(defaults["symprec"] or 1e-8)
    except Exception:                                              # noqa: BLE001
        return 0.01


def mesh_ir_fix_state(mesh, symprec):
    """patch_ir_fix（V128）：mesh.h5 是用哪套不可约映射写的。
    True = IR_FIX（spglib 用正确晶格），False = AMSET 原式（转置晶格），None = 两者相同（晶格矩阵对称，如立方）。
    都对不上 -> ValueError。

    重积分用的是 AMSET 的对称约化四面体积分（transport DOS 先把积分量加到不可约点上），
    重建时的不可约映射必须与作业相同，否则复现检查在 1e-6 的容差下必然 FAIL（MoS₂ IR_FIX-only 实测）。"""
    import spglib
    st = mesh["structure"]
    k = np.asarray(mesh["kpoints"], float)
    dims = [len(np.unique(np.round(np.mod(k[:, i], 1.0), 6) % 1.0)) for i in range(3)]
    L = np.ascontiguousarray(st.lattice.matrix, dtype="double")
    pos = np.ascontiguousarray(st.frac_coords, dtype="double")
    num = np.asarray([s.Z for s in st.species], dtype="intc")
    tr = not bool(mesh.get("soc", False))

    def n_ir(lat):
        m, _ = spglib.get_ir_reciprocal_mesh(dims, (lat, pos, num), symprec=symprec, is_time_reversal=tr)
        return int(len(np.unique(m)))

    n_fix, n_orig = n_ir(L), n_ir(np.ascontiguousarray(L.T))
    n_mesh = int(len(mesh["ir_kpoints"]))
    if n_mesh == n_fix == n_orig:
        return None
    if n_mesh == n_fix:
        return True
    if n_mesh == n_orig:
        return False
    raise ValueError("mesh.h5 的不可约点数 %d 既不是 IR_FIX 的 %d，也不是 AMSET 原式的 %d（网格 %s，symprec %g）"
                     " —— 与当前运行目录不同源" % (n_mesh, n_fix, n_orig, "x".join(map(str, dims)), symprec))


def _apply_ir_fix(run_dir, ir_fix):
    """按 mesh 判定的 IR_FIX 状态设 AZ_IR_FIX，并在需要时加载运行目录的 amset_ir_fix（与作业同一套不可约映射）。
    必须在导入 amset2d_plugin 之前调用（V127 起插件会顺带导入 amset_ir_fix，开关看 AZ_IR_FIX）。"""
    if ir_fix is None:
        return False
    os.environ["AZ_IR_FIX"] = "1" if ir_fix else "0"
    if not ir_fix:
        mod = sys.modules.get("amset_ir_fix")
        if mod is not None and getattr(mod, "STATE", {}).get("applied"):
            raise SystemExit("[ERROR] mesh.h5 是 AMSET 原式不可约映射，但本进程已挂上 amset_ir_fix —— 请单独运行本脚本")
        print("[..] ir_fix：mesh.h5 用的是 AMSET 原式不可约映射（IR_FIX 关）")
        return False
    if not (Path(run_dir) / "amset_ir_fix.py").is_file():
        raise SystemExit("[ERROR] mesh.h5 用的是 IR_FIX 的不可约映射，但运行目录里没有 amset_ir_fix.py —— "
                         "重建的映射会与作业不同，重积分无法复现")
    rd = str(Path(run_dir).resolve())
    if rd not in sys.path:
        sys.path.insert(0, rd)
    import amset_ir_fix
    if not (amset_ir_fix.STATE.get("applied") or amset_ir_fix.apply(force=True)):
        raise SystemExit("[ERROR] mesh.h5 用的是 IR_FIX 的不可约映射，但 amset_ir_fix 在本机 AMSET 上拒绝打补丁")
    print("[..] ir_fix：mesh.h5 用的是 IR_FIX 的不可约映射，已加载 amset_ir_fix 重建")
    return True


def _apply_run_plugins(run_dir, ir_fix=None):
    """patch_kz_cap（V126）：运行目录的 2d_correction.json 开了 kz_cap_rmax 时，加载同目录的 amset2d_plugin，
    它按与作业相同的判据截断 k_z 方向的 equivalence。V125 之前这里用原版插值重建网格（263×263×21），
    与 KZ_CAP_2D 作业写出的 mesh.h5（263×263×3）对不上，apply_mesh_metadata 于是拒绝重积分。
    patch_ir_fix（V128）：ir_fix 为 True/False 时先按它设 AZ_IR_FIX（见 mesh_ir_fix_state、_apply_ir_fix）。
    返回是否加载了 amset2d_plugin。"""
    import json as _json
    import os as _os
    _apply_ir_fix(run_dir, ir_fix)
    rec = Path(run_dir) / "2d_correction.json"
    try:
        kz = _json.loads(rec.read_text()).get("kz_cap_rmax") if rec.is_file() else None
    except (OSError, ValueError):
        kz = None
    if not kz:
        return False
    if not (Path(run_dir) / "amset2d_plugin.py").is_file():
        raise SystemExit("[ERROR] 2d_correction.json 开了 kz_cap_rmax，但运行目录里没有 amset2d_plugin.py —— "
                         "重建的网格会与 mesh.h5 对不上")
    _os.environ["AMSET2D_RECORD"] = str(rec.resolve())
    sys.path.insert(0, str(Path(run_dir).resolve()))
    import amset2d_plugin  # noqa: F401
    print("[..] kz_cap：已按 2d_correction.json（kz_cap_rmax=%s）加载 amset2d_plugin 重建插值网格" % kz)
    return True


def _apply_fermi_fix(run_dir):
    """patch_fermi_fix（V132）：与作业同一套费米能级求解 —— 加载运行目录（或本脚本旁）的 amset_fermi_fix。
    返回是否加载上。"""
    for d in (Path(run_dir).resolve(), Path(__file__).resolve().parent):
        if (d / "amset_fermi_fix.py").is_file():
            if str(d) not in sys.path:
                sys.path.insert(0, str(d))
            break
    try:
        import amset_fermi_fix  # noqa: F401
    except ImportError:
        print("[WARN] 找不到 amset_fermi_fix.py —— 费米能级用 AMSET 原搜索（低温低掺杂可能失败）")
        return False
    return True


def build_amset_data(run_dir, nworkers=None, progress_bar=False, ir_fix=None, set_doping=True):
    """用 vasprun.xml + settings.yaml 重建 AmsetData（能带/速度/DOS/费米能级）。

    只做 interpolation + DOS + fermi levels，不做 overlap、不做 scattering；
    散射率随后从 mesh.h5 注入。ir_fix：mesh_ir_fix_state 的结果（None = 不管）。
    set_doping=False：到 DOS 为止，不求费米能级（tools/fermi_probe.py 用它诊断"求不到费米能级"）。
    """
    import amset
    from amset.core.run import Runner
    from amset.interpolation.bandstructure import Interpolator

    _apply_run_plugins(run_dir, ir_fix)  # patch_kz_cap / patch_ir_fix：与作业同一套插值与不可约映射
    _apply_fermi_fix(run_dir)            # patch_fermi_fix（V132）：与作业同一套费米能级求解
    runner = Runner.from_directory(run_dir)
    s = runner.settings
    nw = s["nworkers"] if nworkers is None else int(nworkers)

    interp = Interpolator(
        runner._band_structure,
        num_electrons=runner._num_electrons,
        interpolation_factor=s["interpolation_factor"],
        soc=s["soc"],
    )
    ad = interp.get_amset_data(
        energy_cutoff=s["energy_cutoff"],
        scissor=s["scissor"],
        bandgap=s["bandgap"],
        symprec=s["symprec"],
        nworkers=nw,
    )
    ad.calculate_dos(estep=s["dos_estep"], progress_bar=progress_bar)
    if not set_doping:
        return ad, s, amset.__version__
    ad.set_doping_and_temperatures(s["doping"], s["temperatures"])
    # patch_carrier_guard（2026-09-28 用户批准）：已实现 vs 请求载流子浓度（>1% 报错）。
    check_achieved_doping(ad)
    return ad, s, amset.__version__


def check_achieved_doping(ad, tol=0.01):
    """已实现载流子浓度 vs 请求掺杂（2026-09-28 用户批准）。

    AMSET 的 AmsetData.set_doping_and_temperatures 在解每个 (掺杂, 温度) 的费米能级时，
    会把 DOS 积分出来的**实际**电子/空穴浓度存进 ad.electron_conc / ad.hole_conc。
    正常网格下它应与请求值一致（AMSET 的容差阶梯最宽会接受 100% 偏差；V132 起 amset_fermi_fix 把关到 1%）；一旦偏差 > tol，说明 DOS/网格/
    费米求解没对上 —— 单位换算错、网格太粗、或带边被切掉，这种结果不可信，直接报错。
    （2026-09-27 MoS2：DOPING 单位事故 -> 费米能级 -153.9 eV，远在 DOS 窗口之外。）
    """
    import numpy as np
    from amset.constants import bohr_to_cm
    doping = getattr(ad, "doping", None)
    if doping is None:
        return
    doping = np.asarray(doping, float)
    ec = np.asarray(getattr(ad, "electron_conc", []), float)
    hc = np.asarray(getattr(ad, "hole_conc", []), float)
    if doping.size == 0 or ec.size == 0 or hc.size == 0:
        return
    to_cm3 = (1.0 / bohr_to_cm) ** 3
    temps = np.asarray(getattr(ad, "temperatures", []), float)
    # ★ 口径：AMSET 的 Fermi 求解目标就是**带符号的净浓度** p - n（get_doping 返回
    #   conc = vb_conc - cb_conc，负 = n 型），所以这里必须比 p-n 与 doping，
    #   不能比"多数载流子" —— 本征区 n_i 很大时 cb_conc 本来就不等于 |doping|。
    req = doping.copy()
    ach = hc - ec                                          # 带符号净浓度（负 = n 型）
    nz = req != 0
    if not np.any(nz):
        print("[..] 载流子浓度核对：doping 全为 0，跳过")
        return
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(ach / req[:, None] - 1.0)
    rel = np.where(nz[:, None], rel, np.nan)
    i, j = [int(x) for x in np.unravel_index(np.nanargmax(rel), rel.shape)]
    worst = float(rel[i, j])
    print("[..] 载流子浓度核对：max|achieved(p-n)/requested-1| = %.3e（索引 %d，请求 %.3e cm^-3 @ %g K，实际 %.3e cm^-3）" % (worst, i, req[i] * to_cm3, float(temps[j]), ach[i, j] * to_cm3))
    if worst > tol:
        raise SystemExit("[ERROR] 已实现载流子浓度与请求掺杂偏差 %.2f%% > %.0f%%（索引 %d @ %g K：请求 %.3e，实际 %.3e cm^-3）——DOS/费米求解没对上，本结果不可信。常见根因：① DOPING 单位写错（2D 应为面密度 cm^-2）；② 插值网格太粗（见 settings.yaml / interpolation_info.json）；③ vasprun 与 h5 不同源（S3 与 S3b 网格不一致）；④ energy_cutoff/带隙设置把带边切掉了。" % (worst * 100, tol * 100, i, float(temps[j]), req[i] * to_cm3, ach[i, j] * to_cm3))


def coordinate_lookup(query_kpts, all_kpts, decimals=6, tol=1e-5):
    """在 all_kpts 里查找 query_kpts 的下标（分数坐标，按周期性等价）。

    用于不依赖不可约 k 点选择顺序地核对能量、以及在不同 AMSET/spglib 版本
    下建立 k 点排列。找不到（或出现重复匹配）时抛错。

    ★ 2026-09-27 修（2D 后处理崩溃）：原实现只做"四舍五入后元组精确相等"，
    跨 AMSET/spglib 版本或 float32 存储时会因为末位差异（如 0.333333 vs
    0.33333333、-1e-9 vs 1e-9）**一个都匹配不上**而抛错。现在改为两级：
      1) 快路径：四舍五入元组精确匹配（与原行为一致）；
      2) 慢路径：**最小镜像最近邻 + 容差 tol**，唯一命中才接受。
    容差内命中多于一个仍然报错（不静默取第一个）——真正重复的网格不该被猜。
    """
    q = np.round(np.asarray(query_kpts, float), decimals)
    a = np.round(np.asarray(all_kpts, float), decimals)
    index = {}
    for i, row in enumerate(a):
        index.setdefault(tuple(row), []).append(i)

    a_raw = np.asarray(all_kpts, float)
    q_raw = np.asarray(query_kpts, float)
    out = []
    for qi in range(len(q)):
        row = q[qi]
        hits = index.get(tuple(row), [])
        if len(hits) == 1:
            out.append(hits[0])
            continue
        if len(hits) > 1:
            raise ValueError("k 点查找失败：%s 在重建网格里命中 %d 个" % (row, len(hits)))
        # 慢路径：最小镜像最近邻（分数坐标周期性，d -= rint(d)）
        d = a_raw - q_raw[qi]
        d -= np.rint(d)
        dist = np.linalg.norm(d, axis=1)
        m = int(np.argmin(dist))
        near = np.where(dist <= tol)[0]
        if dist[m] > tol:
            raise ValueError("k 点查找失败：%s 在重建网格里找不到"
                             "（最近距离 %.3g > tol %.3g）" % (row, dist[m], tol))
        if len(near) != 1:
            raise ValueError("k 点查找失败：%s 在容差 %g 内命中 %d 个" % (row, tol, len(near)))
        out.append(m)
    return np.asarray(out, dtype=int)


def kpoint_permutation(from_kpts, to_kpts, decimals=6):
    """求排列 perm，使 from_kpts[perm] == to_kpts（两者是同一套全网格）。"""
    perm = coordinate_lookup(to_kpts, from_kpts, decimals=decimals)
    return perm


def apply_mesh_metadata(ad, mesh, tol=1e-6):
    """把 mesh.h5 里的 fd_cutoffs / fermi_levels 覆盖回 ad，并校验一致性。

    校验项：全网格 kpoints（顺序可不同，用坐标匹配）、doping、temperatures、
    每个自旋在 mesh 不可约点上的能量。
    任何一项对不上都说明生成 mesh 时的 settings/interpolation_factor/vasprun
    与当前不一致，此时重积分没有意义 —— 直接报错，不静默出结果。

    Returns:
        perm: None（全网格顺序一致）或整数排列（mesh 顺序 -> ad 顺序）。
    """
    from amset.constants import hartree_to_ev, bohr_to_cm

    problems = []
    mesh_kpts = np.asarray(mesh["kpoints"])
    if mesh_kpts.shape != ad.kpoints.shape:
        problems.append("mesh.kpoints 形状 %s != 重建 %s"
                        % (mesh_kpts.shape, ad.kpoints.shape))
        perm = None
    elif np.allclose(mesh_kpts, ad.kpoints, atol=1e-8):
        perm = None
    else:
        # 全网格顺序不同（跨 spglib 版本常见），用坐标求排列
        perm = kpoint_permutation(mesh_kpts, ad.kpoints)

    if not np.allclose(np.asarray(mesh["doping"]),
                       ad.doping * (1 / bohr_to_cm) ** 3, rtol=1e-6):
        problems.append("mesh.doping 与 settings.doping 不一致")
    if not np.allclose(np.asarray(mesh["temperatures"]), ad.temperatures):
        problems.append("mesh.temperatures 与 settings.temperatures 不一致")

    # patch_ir_fix（V128）：重积分的对称约化依赖不可约映射，点数不同 = IR_FIX 状态与作业不同
    ad_ir = getattr(ad, "ir_kpoints_idx", None)
    if ad_ir is not None and len(ad_ir) != len(mesh["ir_kpoints"]):
        problems.append("不可约点数：mesh %d != 重建 %d（IR_FIX 状态与作业不同）"
                        % (len(mesh["ir_kpoints"]), len(ad_ir)))

    # mesh 的不可约点顺序 = mesh["ir_kpoints"] 的顺序（to_dict 里按 ir_kpoints_idx 取列）；
    # 用坐标在 ad 的全网格里找位置，逐点核对能量（与不可约选择无关）。
    ad_ir_idx = coordinate_lookup(mesh["ir_kpoints"], ad.kpoints)
    for spin in ad.spins:
        if spin not in mesh["energies"]:
            problems.append("mesh 缺少自旋 %s 的 energies" % spin)
            continue
        e_mesh = np.asarray(mesh["energies"][spin])
        e_fresh = ad.energies[spin][:, ad_ir_idx] * hartree_to_ev
        if e_mesh.shape != e_fresh.shape or not np.allclose(e_mesh, e_fresh, atol=1e-6):
            problems.append("mesh.energies[%s] 与重建能带不一致（maxdiff=%.3g）"
                            % (spin, np.max(np.abs(e_mesh - e_fresh))))
    if problems:
        raise ValueError("mesh.h5 与当前运行目录不同源，拒绝重积分：\n  - "
                         + "\n  - ".join(problems))

    # 用 mesh 的 fd_cutoffs / fermi_levels，保证与原始运行完全同口径
    if "fd_cutoffs" in mesh:
        ad.fd_cutoffs = tuple(np.asarray(mesh["fd_cutoffs"], float) / hartree_to_ev)
    fl = mesh.get("fermi_levels")
    if fl is not None:
        fl = np.asarray(fl, float) / hartree_to_ev
        if fl.shape == ad.fermi_levels.shape:
            ad.fermi_levels = fl
    return perm


def inject_mesh_rates(ad, mesh, labels=None, perm=None):
    """把 mesh.h5 的每机制散射率展开到全网格并 set 进 ad。

    展开用 **mesh 自带的** ir_to_full_kpoint_mapping（而不是重建的），这样即使
    跨 AMSET/spglib 版本导致不可约点选择不同，也能逐位还原原始全网格散射率；
    perm 负责把 mesh 的全网格顺序换算成 ad 的顺序。

    Returns:
        (rates_full, labels)
    """
    labels = list(labels if labels is not None else mesh["scattering_labels"])
    rates_full = {}
    for spin in ad.spins:
        if spin not in mesh["scattering_rates"]:
            raise ValueError("mesh.h5 缺少自旋 %s 的 scattering_rates" % spin)
        r_ir = np.asarray(mesh["scattering_rates"][spin])
        if r_ir.shape[0] != len(labels):
            raise ValueError("mesh 散射率机制数与 scattering_labels 不一致：%s vs %d"
                             % (r_ir.shape[0], len(labels)))
        r_full = expand_ir_to_full(r_ir, mesh["ir_to_full_kpoint_mapping"])
        if perm is not None:
            r_full = r_full[..., perm]
        expected = (len(labels),) + ad.fermi_levels.shape + ad.energies[spin].shape
        if r_full.shape != expected:
            raise ValueError("展开后的散射率形状 %s != 期望 %s"
                             % (r_full.shape, expected))
        rates_full[spin] = r_full
    ad.set_scattering_rates(rates_full, labels)
    return rates_full, labels


# =====================================================================
# 输运积分
# =====================================================================
def integrate(ad, labels, separate_labels=()):
    """调用 AMSET 自己的输运积分函数。

    Args:
        ad: 已注入散射率的 AmsetData（机制行按需置零）。
        labels: 全部机制名（与 ad.scattering_rates 轴 0 对应）。
        separate_labels: 需要单独给出迁移率的机制名（保证非零机制）。

    Returns:
        dict: sigma(S/m)、seebeck(uV/K)、kappa、mobility{...}（cm^2/Vs）。
    """
    from amset.core.transport import (
        solve_boltzman_transport_equation,
        _calculate_mobility,
    )

    # 使用 AMSET 的公开入口重算 σ/S/κ 与 overall μ。
    # separate_mobility=False：置零机制若单独求 mobility 会 1/0 -> SVD 不收敛。
    sigma, seebeck, kappa, mobility = solve_boltzman_transport_equation(
        ad,
        calculate_mobility=True,
        separate_mobility=False,
        progress_bar=False,
    )
    out = {
        "conductivity_S_m": np.asarray(sigma),
        "seebeck_uV_K": np.asarray(seebeck),
        "electronic_thermal_conductivity": np.asarray(kappa),
        "mobility_cm2_Vs_s": {},
    }
    if mobility is not None:
        out["mobility_cm2_Vs_s"]["overall"] = np.asarray(mobility["overall"])

    # 对保留机制逐一调用 AMSET 的 _calculate_mobility（公开入口的单项分支）
    for name in separate_labels:
        if name not in labels:
            continue
        i = labels.index(name)
        if all(np.all(np.asarray(ad.scattering_rates[s])[i] == 0.0) for s in ad.spins):
            continue  # 被置零的机制跳过
        out["mobility_cm2_Vs_s"][name] = np.asarray(_calculate_mobility(ad, i))
    return out


# =====================================================================
# 输出
# =====================================================================
def _doping_marks(run_dir, doping_cm3, temperatures_K, seebeck_inplane_uV_K, n_atoms=None):
    """返回 {doping_index: [标记...]}：每**原子**超阈值 -> "超出刚带近似"；S 符号与载流子类型不一致 -> "sign_anomaly"。

    ★ 2026-09-28 用户批准：阈值口径从"每原胞"改成"每**原子**"。每原胞数取决于晶胞大小，
      C60 网络那种大原胞（约 146 A^2、60 原子）在出厂掺杂上限下会被误标为"超出刚带近似"。
      优先读 2d_correction.json 的 carrier_per_atom_threshold（默认 0.02）；
      只有旧键 carrier_per_cell_threshold 时按 n_atoms 折算；n_atoms 未知则退回旧的每原胞判据。
    """
    marks = {}
    c_len_cm = area_cm2 = thr_atom = thr_cell = None
    corr = Path(run_dir) / "2d_correction.json"
    if corr.is_file():
        try:
            d = json.load(open(corr, encoding="utf-8"))
            c_len_cm = d.get("areal_density_factor_cm")
            area_cm2 = d.get("inplane_cell_area_cm2")
            n_atoms = n_atoms or d.get("num_atoms")
            thr_atom = d.get("carrier_per_atom_threshold")
            thr_cell = d.get("carrier_per_cell_threshold")
        except (OSError, ValueError, TypeError):
            pass
    if not thr_atom:
        thr_atom = (float(thr_cell) / float(n_atoms)) if (thr_cell and n_atoms) else 0.02
    if not thr_cell:
        thr_cell = 0.05
    # temperatures_K 可能是 numpy.ndarray（ad.temperatures 直传）或 list，
    # 统一转 list 再取索引，避免 ndarray 没有 .index() 报错。
    t_list = list(np.asarray(temperatures_K, dtype=float))
    if 300.0 in t_list:
        iT = t_list.index(300.0)
    elif t_list:
        iT = len(t_list) // 2
    else:
        iT = 0
    for i, n3 in enumerate(doping_cm3):
        m = []
        if c_len_cm and area_cm2:
            per_cell = abs(float(n3)) * float(c_len_cm) * float(area_cm2)
            if n_atoms:
                if per_cell / float(n_atoms) > float(thr_atom):
                    m.append("超出刚带近似")
            elif per_cell > float(thr_cell):
                m.append("超出刚带近似")
        try:
            S = float(seebeck_inplane_uV_K[i][iT])
        except (IndexError, TypeError, ValueError):
            S = 0.0
        if (float(n3) < 0 and S > 0) or (float(n3) > 0 and S < 0):
            m.append("sign_anomaly")
        if m:
            marks[str(i)] = m
    return marks


def build_result(run_dir, mesh_file, amset_version, settings, all_labels,
                 dropped, kept, transport, doping_cm3, temperatures,
                 reproduce=None, extra_notes=(), n_atoms=None):
    mob = transport.get("mobility_cm2_Vs_s", {})
    result = {
        "schema": "autozt.intrinsic_transport/1",
        "script": "postprocess_intrinsic.py",
        "script_version": __version__,
        "amset_version": amset_version,
        "run_dir": str(run_dir),
        "mesh_file": str(mesh_file),
        "settings_file": str(Path(run_dir) / "settings.yaml"),
        "doping_cm3": to_jsonable(doping_cm3),
        "temperatures_K": to_jsonable(temperatures),
        "scattering_labels_all": list(all_labels),
        "intrinsic_dropped": list(dropped),
        "intrinsic_labels": list(kept),
        "interpolation_factor": settings.get("interpolation_factor"),
        "notes": list(extra_notes),
        "intrinsic": {
            "mobility_cm2_Vs_s": {k: to_jsonable(v) for k, v in mob.items()},
            "mobility_inplane_cm2_Vs_s": {
                k: to_jsonable(inplane_average(v)) for k, v in mob.items()},
            "seebeck_uV_K": to_jsonable(transport["seebeck_uV_K"]),
            "seebeck_inplane_uV_K": to_jsonable(inplane_average(transport["seebeck_uV_K"])),
            "conductivity_S_m": to_jsonable(transport["conductivity_S_m"]),
            "conductivity_inplane_S_m": to_jsonable(
                inplane_average(transport["conductivity_S_m"])),
        },
        "doping_marks": _doping_marks(
            run_dir, doping_cm3, temperatures,
            inplane_average(transport["seebeck_uV_K"]), n_atoms=n_atoms),
    }
    if reproduce is not None:
        result["reproduction_check"] = reproduce
    return result


# =====================================================================
# CLI
# =====================================================================
def find_mesh(run_dir, explicit=None):
    if explicit:
        p = Path(explicit)
        if not p.is_absolute():
            p = Path(run_dir) / explicit
        if not p.is_file():
            raise FileNotFoundError("找不到 mesh 文件：%s" % p)
        return p
    cands = sorted(glob.glob(str(Path(run_dir) / "mesh*.h5")))
    if not cands:
        raise FileNotFoundError(
            "运行目录 %s 里没有 mesh*.h5 —— settings.yaml 要写 write_mesh: true "
            "并重跑 S8.4 才会产出。" % run_dir)
    # 多个时取最新
    return Path(max(cands, key=os.path.getmtime))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="从 AMSET mesh.h5 严格重算本征（剔除 IMP）μ/S/σ")
    ap.add_argument("run_dir", nargs="?", default=".",
                    help="S8.4 运行目录（含 settings.yaml/vasprun.xml/mesh*.h5），默认当前目录")
    ap.add_argument("--mesh", default=None, help="指定 mesh.h5（默认自动找 mesh*.h5）")
    ap.add_argument("--drop", default=",".join(DEFAULT_DROP),
                    help="要剔除的机制，逗号分隔（默认 IMP）")
    ap.add_argument("--mode", choices=["zero", "drop"], default="zero",
                    help="zero=置零（默认，用户方案）/ drop=整行删除（overall 等价）")
    ap.add_argument("--out", default="intrinsic_transport.json", help="输出 JSON 文件名")
    ap.add_argument("--nworkers", type=int, default=None, help="插值并行数（默认用 settings）")
    ap.add_argument("--check-reproduce", action="store_true",
                    help="先用全部机制重积分并与 transport*.json 对比，验证逐位复现")
    ap.add_argument("--reproduce-tol", type=float, default=1e-6,
                    help="复现检查的容差（绝对），默认 1e-6")
    args = ap.parse_args(argv)

    run_dir = Path(args.run_dir).resolve()
    mesh_file = find_mesh(run_dir, args.mesh)
    drop = [x for x in str(args.drop).replace(" ", ",").split(",") if x]

    # 1) 读 mesh，判定作业用的不可约映射（V128），再按同一套重建 AmsetData（不重算散射）
    mesh = read_mesh_h5(mesh_file)
    ir_fix = mesh_ir_fix_state(mesh, _settings_symprec(run_dir))
    ad, settings, amset_version = build_amset_data(
        run_dir, nworkers=args.nworkers, progress_bar=False, ir_fix=ir_fix)

    # 2) 口径校验
    labels = list(mesh["scattering_labels"])
    perm = apply_mesh_metadata(ad, mesh)
    if perm is not None:
        print("[..] mesh 全网格顺序与重建不同（跨 spglib 版本），已按坐标重排")
    rates_full, labels = inject_mesh_rates(ad, mesh, labels, perm=perm)
    print("[..] amset=%s  mechanisms=%s  spins=%s  n_full=%d"
          % (amset_version, labels, [str(s) for s in ad.spins], len(ad.kpoints)))

    # 3) 复现检查（可选）：用全部机制重积分，对比原 transport*.json
    reproduce = None
    if args.check_reproduce:
        from amset.core.transport import solve_boltzman_transport_equation
        ad.set_scattering_rates({s: rates_full[s] for s in rates_full}, labels)
        sig0, seeb0, kap0, mob0 = solve_boltzman_transport_equation(
            ad, calculate_mobility=not ad.is_metal, separate_mobility=True,
            progress_bar=False)
        # 多个 transport*.json 时取**最新**（与 find_mesh 同口径）：retry/重跑时目录里
        # 可能残留上一次的 transport_<旧网格>.json，按字典序可能取到旧文件 -> 复现假 FAIL。
        tf = sorted(glob.glob(str(run_dir / "transport*.json")),
                    key=os.path.getmtime, reverse=True)
        reproduce = {"transport_json": os.path.basename(tf[0]) if tf else None}
        if tf:
            j = json.load(open(tf[0]))
            diffs = {}
            for name in ["overall"] + labels:
                if name in j.get("mobility", {}) and name in mob0:
                    a = np.asarray(mob0[name]); b = np.asarray(j["mobility"][name])
                    n = min(a.shape[0], b.shape[0])
                    diffs[name] = float(np.nanmax(np.abs(a[:n] - b[:n])))
            diffs["conductivity_S_m"] = float(np.nanmax(np.abs(sig0 - np.asarray(j["conductivity"]))))
            diffs["seebeck_uV_K"] = float(np.nanmax(np.abs(seeb0 - np.asarray(j["seebeck"]))))
            reproduce["max_abs_diff"] = diffs
            worst = max(diffs.values()) if diffs else 0.0
            reproduce["passed"] = bool(worst <= args.reproduce_tol)
            print("[..] 复现检查：max|diff|=%s -> %s"
                  % ("%.3e" % worst, "PASS" if reproduce["passed"] else "FAIL"))
            if not reproduce["passed"]:
                raise SystemExit("[ERROR] 重积分无法复现原 transport.json —— "
                                 "mesh/settings/vasprun 不同源，勿信本征结果。")
        # 复原为原始（全机制）散射率，后续再置零
        ad.set_scattering_rates({s: rates_full[s] for s in rates_full}, labels)

    # 4) 置零/剔除目标机制，重积分得本征结果
    drop_idx = resolve_mechanism_indices(labels, drop)
    if not drop_idx:
        print("[WARN] 机制 %s 不在 mesh labels %s 里 —— 本征结果 = 原结果"
              % (drop, labels), file=sys.stderr)
        dropped, kept = [], list(labels)
        transport = integrate(ad, labels, separate_labels=kept)
    else:
        if args.mode == "drop":
            rates_mod, kept, dropped = drop_mechanisms(rates_full, labels, drop)
        else:
            rates_mod, dropped, kept = zero_mechanisms(rates_full, labels, drop)
        ad.set_scattering_rates(rates_mod, labels)
        transport = integrate(ad, labels, separate_labels=kept)
        print("[OK] 本征口径：剔除 %s，保留 %s" % (dropped, kept))

    # 5) 输出
    notes = []
    if args.mode == "zero":
        notes.append("IMP 行已置零后重积分；overall=sum(ADP,POP)，单项由 AMSET "
                     "_calculate_mobility 对保留机制逐一给出。")
    notes.append("各机制散射率独立，故本结果 == 单独运行 scattering_type=[%s]。"
                 % ",".join(kept))
    if ad.is_metal:
        notes.append("体系 is_metal=True，AMSET 拒绝给迁移率；仅 σ/S 有效。")
    from amset.constants import bohr_to_cm
    result = build_result(
        run_dir, mesh_file, amset_version, settings, labels,
        dropped, kept, transport,
        doping_cm3=ad.doping * (1 / bohr_to_cm) ** 3,
        temperatures=ad.temperatures, reproduce=reproduce, extra_notes=notes,
        n_atoms=int(ad.structure.num_sites),
    )
    out = Path(args.out)
    if not out.is_absolute():
        out = run_dir / out
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print("[DONE] 本征输运已写出：%s" % out)
    mi = result["intrinsic"]["mobility_inplane_cm2_Vs_s"]
    if "overall" in mi:
        t300 = list(np.asarray(result["temperatures_K"])).index(300) if 300 in list(np.asarray(result["temperatures_K"])) else 0
        print("[..] T=%g K 面内本征迁移率 overall = %.3f cm^2/Vs"
              % (result["temperatures_K"][t300], mi["overall"][0][t300]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
