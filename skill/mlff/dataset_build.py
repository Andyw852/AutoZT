#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dataset_build.py —— step6_dataset 引擎（在 venv 里跑，需要 numpy/ase/matplotlib）。

cwd = 材料目录。把 step5_label 全部已算完的 cfg-* OUTCAR 转成 extxyz 训练集：
    - 逐帧解析能量（energy(sigma->0) 优先）/ 力 / 应力（in kB 张量 → eV/Å³ Voigt）
    - 指纹校验（vs step1 的 DFT 设置；k 点密度用容差）
    - extend 模式并入 PRE_XYZ_FILES（指纹不一致的老数据整批丢弃 + WARN）
    - 离群过滤（ENERGY_LIMIT / FORCE_LIMIT；>10% WARN，>30% FAIL；filtered: true 保留）
    - FPS 全局排序 + 固定测试集划分（id 哈希 %10，跨代稳定）
    - e0s.json（孤立原子能量）、coverage_report.json、energy_forces.png 三联图
    - 写 step6_dataset/gen-<K>/{train.xyz, test.xyz, ...}

--coverage-only 模式（step4 extend 用）：只对 PRE_XYZ_FILES 做覆盖分析，把新采样
计划写进 coverage_plan.json（把新钱压到盲区：应变档/幅度档/元素）。
退出码 0 成功；非 0 = [ERROR]。
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dft_settings as ds  # noqa: E402
import mlff_common as mc  # noqa: E402

EV_A3_PER_KBAR = 6.24150907e-4   # 1 kB = 0.1 GPa；1 eV/Å³ = 160.217662 GPa


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gen", type=int, required=True)
    p.add_argument("--outdir", default=None)            # step6_dataset/gen-<K>
    p.add_argument("--step1-dir", default="step1_relax")
    p.add_argument("--step5-dir", default="step5_label")
    p.add_argument("--step4-dir", default="step4_genstruct")
    p.add_argument("--dim", required=True)
    p.add_argument("--energy-limit", type=float, default=0.005)
    p.add_argument("--force-limit", type=float, default=0.1)
    p.add_argument("--kspacing-tol", type=float, default=0.20)
    # ★ 2026-09-26 user 要求：k 点间距超容差 → FAIL（默认关 = 只 WARN，兼容旧行为）
    p.add_argument("--kspacing-tol-block", action="store_true")
    p.add_argument("--pre-xyz", default="")             # extend：逗号分隔
    p.add_argument("--fps-seed", type=int, default=42)
    p.add_argument("--vol-factors", default="0.97,1.00,1.03")
    p.add_argument("--n-per-cell", type=int, default=2)
    p.add_argument("--coverage-only", action="store_true")
    return p.parse_args()


# ============================================================ OUTCAR 解析
def outcar_done(path):
    p = Path(path)
    if not p.is_file():
        return False
    try:
        with open(p, "rb") as fh:
            try:
                fh.seek(-200000, 2)
            except OSError:
                fh.seek(0)
            tail = fh.read().decode("utf-8", "ignore")
    except OSError:
        return False
    return "General timing and accounting informations" in tail


def parse_outcar(path):
    """-> dict(E, F, stress) 或 None（没算完）。E 用 energy(sigma->0) 兜底 TOTEN。"""
    if not outcar_done(path):
        return None
    text = Path(path).read_text(errors="ignore")
    e = None
    for m in list(__import__("re").finditer(r"energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)", text)):
        e = float(m.group(1))
    if e is None:
        m = list(__import__("re").finditer(r"free\s+energy\s+TOTEN\s*=\s*([-+0-9.Ee]+)", text))
        if m:
            e = float(m[-1].group(1))
    if e is None:
        raise RuntimeError("OUTCAR 里找不到总能")

    # 力：最后一个 TOTAL-FORCE 块（外层整体捕获：内层 + 会让 group 只剩最后一行）
    fblocks = [b for b in __import__("re").finditer(
        r"TOTAL-FORCE \(eV/Angst\)\s*\n\s*-+\s*\n"
        r"((?:(?:\s*[-+0-9.eE]+){3,6}[^\n]*\n)+)", text)]
    if not fblocks:
        raise RuntimeError("OUTCAR 里找不到 TOTAL-FORCE 块")
    f = []
    for line in fblocks[-1].group(1).splitlines():
        t = line.split()
        # VASP 6.x：6 列 = 位置 x y z + 力 fx fy fz（力在后 3 列）；旧版只有 3 列力
        if len(t) >= 6:
            cols = t[3:6]
        elif len(t) == 3:
            cols = t
        else:
            continue
        try:
            f.append([float(cols[0]), float(cols[1]), float(cols[2])])
        except ValueError:
            continue
    if not f:
        raise RuntimeError("TOTAL-FORCE 块解析不出力")

    # 应力：最后一个 "in kB" 行（与 in kB 同行的 6 个数，列序 XX YY ZZ XY YZ ZX）
    # 符号：VASP 的 kB 张量以压缩为正（external pressure = -Tr/3 > 0 表示压缩）；
    # MACE/ASE 以拉伸为正 → 取负号。
    import re as _re
    srows = _re.findall(
        r"in\s+kB\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)", text)
    stress = None
    if srows:
        xx, yy, zz, xy, yz, zx = [float(x) for x in srows[-1]]
        stress = [-v * EV_A3_PER_KBAR for v in (xx, yy, zz, yz, zx, xy)]
    return {"E": e, "F": f, "stress": stress}


def outcar_vasp_version(path):
    """从 OUTCAR 首行提取 VASP 版本（如 "vasp.6.6.0" / "vasp.6.4.3"）。

    dft_fingerprint 只含 INCAR 键 + POTCAR + k 点密度，没有版本项——同一设置
    用不同 VASP 版本（如 CPU 6.4.3 vs GPU 6.6.0）算的帧指纹完全相同，混训时
    靠指纹发现不了。dataset_build 逐帧记版本，dataset_summary.json 汇总 +
    混用 WARN（README §15 补记）。"""
    try:
        import re as _re
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            first = fh.readline()
        m = _re.search(r"(vasp\.[0-9]+\.[0-9]+\.[0-9]+)", first)
        return m.group(1) if m else None
    except OSError:
        return None


# ============================================================ 特征 / 过滤 / FPS
def fps_order(features, seed=42):
    """贪心最远点采样，返回索引顺序（前缀天然嵌套）。"""
    n = len(features)
    if n == 0:
        return []
    rng = np.random.default_rng(seed)
    order = [int(rng.integers(0, n))]
    remain = list(set(range(n)) - {order[0]})
    dist = np.full(n, np.inf)
    while remain:
        d = np.sum((features[remain] - features[order[-1]]) ** 2, axis=1)
        dist[remain] = np.minimum(dist[remain], d)
        nxt = remain[int(np.argmax(dist[remain]))]
        order.append(nxt)
        remain.remove(nxt)
    return order


def frame_features(E, F, n_atom):
    f = np.array(F, dtype=float).ravel()
    f = f / (max(np.std(f), 1e-12) or 1.0)
    return np.concatenate([[E / n_atom], f])


def filter_frames(frames, energy_limit, force_limit):
    """两道过滤（§8.3）：
    - 力：max|F| ≥ FORCE_LIMIT（= autoplex data_distillation 的 force_max，源码
      默认 40.0 eV/Å，力分量上限——不是 0.1！0.1 会把大幅度 rattle 全滤掉）
    - 能量：|E/atom − 组中位数| ≥ ENERGY_LIMIT（组 = config_type + 应变 + 幅度档；
      组员 < 4 时跳过能量过滤——2 个种子定不出中位数，硬套会误杀）
    """
    groups = {}
    for fr in frames:
        key = (fr["config_type"], fr.get("strain_factor"), fr.get("strain_grun"),
               fr.get("rattle_std"))
        groups.setdefault(key, []).append(fr)
    n_force, n_energy = 0, 0
    for key, frs in groups.items():
        # rattle 帧跳过能量过滤：随机位移能量散布是物理的（幅度 0.27 Å 时同组可差
        # ~0.09 eV/atom >> ENERGY_LIMIT=0.005），跨代累积组员>=4 后固定阈值必误杀；
        # 坏结构由 FORCE_LIMIT 与 d_min 兜底。
        is_rattle = frs[0]["config_type"] == "rattle"
        use_energy = len(frs) >= 4 and not is_rattle
        med = float(np.median([f["E"] / len(f["F"]) for f in frs])) if use_energy else 0.0
        for fr in frs:
            ea = fr["E"] / len(fr["F"])
            fmax = float(np.max(np.abs(fr["F"])))
            bad_f = fmax >= force_limit
            bad_e = use_energy and abs(ea - med) >= energy_limit
            fr["filtered"] = bool(bad_f or bad_e)
            fr["e_vs_median_eV"] = float(ea - med) if use_energy else None
            fr["fmax_eV_A"] = fmax
            n_force += 1 if bad_f else 0
            n_energy += 1 if bad_e else 0
    return n_force, n_energy


# ============================================================ 主流程
def load_pre_xyz(paths, dim):
    """读 PRE_XYZ_FILES → 统一格式帧列表。读不出指纹的帧标记 fp=None。"""
    from ase.io import read as ase_read
    frames = []
    for p in paths:
        if not Path(p).is_file():
            sys.exit("[ERROR] PRE_XYZ_FILES 里的 %s 不存在" % p)
        try:
            atoms_list = ase_read(p, index=":")
        except Exception as e:
            sys.exit("[ERROR] 读 %s 失败：%s" % (p, e))
        for at in atoms_list:
            info = dict(at.info or {})
            forces = (at.arrays.get("REF_forces")
                      if "REF_forces" in at.arrays
                      else at.arrays.get("forces")
                      if "forces" in at.arrays
                      else info.get("REF_forces", info.get("forces", [])))
            frames.append({
                "E": float(info.get("REF_energy", info.get("energy", 0.0))),
                "F": np.array(forces),
                "stress": info.get("REF_stress", None),
                "atoms": at,
                "config_type": str(info.get("config_type", "pre")),
                "strain_factor": info.get("strain_factor"),
                "strain_grun": info.get("strain_grun"),
                "filtered": bool(info.get("filtered", False)),
                "source": Path(p).name,
                "fingerprint": info.get("mlff_fingerprint"),
            })
    return frames


def coverage_analysis(frames, dim):
    """覆盖分析：体积（2D 面积）/ RMS 位移 / 元素力分布。-> report dict。"""
    vols, rms_vals, elem_forces = [], [], {}
    for fr in frames:
        cell = np.array(fr["atoms"].get_cell(), dtype=float)
        if dim == "2d":
            v = float(np.linalg.norm(np.cross(cell[0], cell[1])))
        else:
            v = float(np.abs(np.linalg.det(cell)))
        vols.append(v)
        pos = fr["atoms"].get_positions()
        ref = pos.mean(axis=0)
        rms_vals.append(float(np.sqrt(np.mean(np.sum((pos - ref) ** 2, axis=1))))
                        if len(pos) > 0 else 0.0)
        for s, f in zip(fr["atoms"].get_chemical_symbols(), np.array(fr["F"])):
            elem_forces.setdefault(s, []).append(float(np.linalg.norm(f)))
    report = {
        "n_frames": len(frames),
        "volume_min": round(min(vols), 4) if vols else None,
        "volume_max": round(max(vols), 4) if vols else None,
        "volume_mean": round(float(np.mean(vols)), 4) if vols else None,
        "volume_std_rel": round(float(np.std(vols) / (np.mean(vols) or 1.0)), 4) if vols else None,
        "rms_min_A": round(min(rms_vals), 5) if rms_vals else None,
        "rms_max_A": round(max(rms_vals), 5) if rms_vals else None,
        "n_per_element": {k: len(v) for k, v in elem_forces.items()},
        "force_median_eV_A": {k: round(float(np.median(v)), 4) for k, v in elem_forces.items()},
    }
    return report


def make_coverage_plan(frames, dim, vol_factors, n_per_cell):
    """把新采样压到覆盖盲区。返回 {targets: [(strain_factor, rattle_std, n)], reason}。
    v1 策略：已有帧体积集中在单值附近（典型：固定胞 rattle）→ 1.00 档种子降 1，
    其余档平分；否则均匀。"""
    vols = []
    for fr in frames:
        cell = np.array(fr["atoms"].get_cell())
        if dim == "2d":
            vols.append(float(np.linalg.norm(np.cross(cell[0], cell[1]))))
        else:
            vols.append(float(np.abs(np.linalg.det(cell))))
    if not vols:
        return None
    rel = float(np.std(vols) / (np.mean(vols) + 1e-12))
    if rel < 1e-3:
        plan = {"targets": [], "reason": ("已有帧体积全部集中在单一晶胞（相对展宽 %.2e）"
                                          "——把新钱全花在应变构型上（1.00 档只留 1 个种子）"
                                          % rel)}
        n_edge = max(1, int(np.ceil((len(vol_factors) * n_per_cell - 1) / 2)))
        for f in vol_factors:
            if abs(f - 1.0) < 1e-9:
                plan["targets"].append({"strain": float(f), "n": 1})
            else:
                plan["targets"].append({"strain": float(f), "n": n_edge})
        return plan
    return None


def main():
    a = parse_args()
    t0 = time.time()
    cwd = Path.cwd()
    outdir = cwd / (a.outdir or "step6_dataset/gen-%d" % a.gen)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- 收集 step5 已算完的帧（跨代累积：遍历 step4 所有 gen-* 清单）----
    from collections import Counter as _Counter
    manifests = sorted((cwd / a.step4_dir).glob("gen-*/struct_manifest.json"))
    frames = []
    frame_fp = []                   # ★ 逐帧 DFT 设置指纹（id, config_type, fp, Δk, note）
    vasp_versions = _Counter()      # 逐帧 VASP 版本（fingerprint 不含版本项，见 outcar_vasp_version）
    for mp in manifests:
        man = json.loads(mp.read_text())
        for ent in man.get("frames", []):
            cfgdir = cwd / a.step5_dir / ent["id"]
            if not outcar_done(cfgdir / "OUTCAR"):
                continue
            try:
                d = parse_outcar(cfgdir / "OUTCAR")
            except Exception as e:
                print("[WARN] %s 解析失败：%s" % (ent["id"], e))
                continue
            vasp_versions[outcar_vasp_version(cfgdir / "OUTCAR") or "unknown"] += 1
            from ase.io import read as ase_read
            atoms = ase_read(str(cwd / a.step4_dir / ("gen-%d" % ent["gen"]) /
                                ent["file"]), format="vasp")
            if len(atoms) != len(d["F"]):
                print("[WARN] %s 原子数 %d ≠ 力行数 %d，跳过" %
                      (ent["id"], len(atoms), len(d["F"])))
                continue
            frames.append({
                "id": ent["id"], "config_type": ent["config_type"],
                "strain_factor": ent.get("strain_factor"),
                "strain_grun": ent.get("strain_grun"),
                "volume_factor": ent.get("volume_factor"),
                "rattle_std": ent.get("rattle_std"),
                "E": d["E"], "F": d["F"], "stress": d["stress"],
                "atoms": atoms, "filtered": False,
                "source": "step5_label",
            })
            # ★ 逐帧 DFT 设置指纹（2026-09-26，user 要求）：S5_label 自己产的帧此前
            #   【从不校验】（旧注释只是"假设本数据集所有 DFT 帧共用 step5 的生成设置"）。
            #   这里逐帧从 cfgdir 的 INCAR/POTCAR/KPOINTS + 该帧超胞晶格算指纹与
            #   胞不变 k 点间距 Δk，循环后统一做"帧间唯一 + 与 step1 参考一致"判定。
            #   Δk 必须用 kdelta_density()（胞不变量）：超胞帧与原胞参考只有它可比。
            try:
                _fi = ds.read_incar(cfgdir / "INCAR")
                _ffp, _ = ds.dft_fingerprint(
                    _fi, cfgdir / "POTCAR", atoms.cell.tolist(), cfgdir / "KPOINTS",
                    mc.read_kv(cwd / a.step1_dir / mc.METHOD_FILE).get("FUNC", "?"))
                _fdk, _fdkn = ds.kdelta_density(cfgdir / "KPOINTS", atoms.cell.tolist())
            except Exception as _e:                      # noqa: BLE001
                _ffp, _fdk, _fdkn = None, 0.0, "取不到（%s）" % _e
            frame_fp.append((ent["id"], ent.get("config_type"), _ffp, _fdk, _fdkn))
    if not frames:
        sys.exit("[ERROR] step5_label 里一帧算完的 OUTCAR 都没有 —— 先让 S5 跑完。")
    if len(vasp_versions) > 1:
        print("[WARN] VASP 版本混用：%s —— dft_fingerprint 不含版本项，靠指纹发现不了；"
              "跨版本帧混训可能带系统误差，请核对是否换过 VASP" % dict(vasp_versions))

    # ---- 指纹（本数据集所有 DFT 帧共用 step5 的生成设置）----
    incar = ds.read_incar(cwd / a.step1_dir / "INCAR")
    method = mc.read_kv(cwd / "step1_relax" / mc.METHOD_FILE)
    from dim_common import read_poscar_cell_frac
    prim_lat, _ = read_poscar_cell_frac(str(cwd / a.step1_dir / "CONTCAR"))
    fp, ksp = ds.dft_fingerprint(incar, cwd / a.step1_dir / "POTCAR",
                                 prim_lat, cwd / a.step1_dir / "KPOINTS",
                                 method.get("FUNC", "?"))

    # ---- ★ 帧间 DFT 设置一致性（2026-09-26，user 要求）----
    #   此前：S5_label 自己产的帧【从不校验】—— 指纹只从 step1_relax 算一次，而且只用于
    #   筛 PRE_XYZ_FILES。于是"某些帧被换了网格/INCAR"，甚至"整批 Γ-only"都能静默进
    #   训练集（2026-09 那个 243 原子 Γ-only 数据集就是这么来的：实测 Γ-only ≡ 原胞
    #   3×3×3，参考构型 max|F|=0.103 eV/Å，而 6×6×6 是 1e-4 —— 差 1000 倍）。
    #   iso 帧故意只有自身元素的 POTCAR（见 gen_step5_label.build_potcar_for_frame），
    #   其指纹天然与 step1 不同，故只对非 iso 帧做这两道判定。
    _bulk_fp = [t for t in frame_fp if t[1] != "iso"]
    if len(_bulk_fp) < len(frame_fp):
        print("[..] 指纹校验跳过 %d 个 iso 帧（单元素 POTCAR，天然与 step1 不同）"
              % (len(frame_fp) - len(_bulk_fp)))
    _groups = {}
    for _fid, _ct, _ffp, _fdk, _fdkn in _bulk_fp:
        _groups.setdefault(_ffp, []).append((_fid, _fdk, _fdkn))
    if len(_groups) > 1:
        print("[ERROR] S5_label 各帧的 DFT 设置指纹不唯一（共 %d 组）：" % len(_groups))
        for _g, _mem in sorted(_groups.items(), key=lambda kv: -len(kv[1])):
            print("          %d 帧，如 %s" % (len(_mem), ", ".join(m[0] for m in _mem[:4])))
            print("            指纹：%s" % ((_g or "(取不到)")[:200]))
        sys.exit(
            "[ERROR] 同一代里混了不止一种 DFT 设置（INCAR/POTCAR）—— 训练集内部不自洽，\n"
            "        混训会带系统误差。\n"
            "        处置：确认 step5_label 的 INCAR/POTCAR 是同一套后 retry 本步；\n"
            "        若中途改过生成设置，请 rerun S5_label 让全部帧用同一套设置重算。")
    _dks = [m[1] for _mem in _groups.values() for m in _mem if m[1] > 0]
    _kmax = max(_dks) if _dks else 0.0
    if _dks and min(_dks) > 0 and (_kmax / min(_dks) - 1.0) > a.kspacing_tol:
        _msg = ("S5_label 帧之间的 k 点间距 Δk 不一致：%.5f–%.5f 1/Å（相对差 %.0f%% > %.0f%%）"
                % (min(_dks), _kmax, (_kmax / min(_dks) - 1.0) * 100, a.kspacing_tol * 100))
        if a.kspacing_tol_block:
            sys.exit("[ERROR] %s\n        处置：统一 KPOINTS_GRID 后 rerun S5_label。"
                     "（KSPACING_TOL_BLOCK=false 可降级为 WARN，不推荐）" % _msg)
        print("[WARN] %s（KSPACING_TOL_BLOCK=false，放行）" % _msg)

    if _bulk_fp and _bulk_fp[0][2] is not None:
        _ref_dk, _ref_note = ds.kdelta_density(cwd / a.step1_dir / "KPOINTS", prim_lat)
        _blk = a.kspacing_tol_block
        if _ref_dk > 0 and _bulk_fp[0][3] > 0:
            print("[..] k 点间距：S5 帧 Δk=%.5f 1/Å vs step1 参考 Δk=%.5f 1/Å（比值 %.2f）"
                  % (_bulk_fp[0][3], _ref_dk, _bulk_fp[0][3] / _ref_dk))
        # ① k 点间距（胞不变量 Δk）：默认 WARN，KSPACING_TOL_BLOCK=on → FAIL
        if _bulk_fp[0][2] == fp:
            _ok, _msg = ds.check_fingerprint(
                _bulk_fp[0][2], _bulk_fp[0][3], fp, _ref_dk,
                kspacing_tol=a.kspacing_tol, kblock=_blk)
            if not _ok:
                sys.exit(
                    "[ERROR] S5_label 帧与 step1 参考的 k 点间距不一致：%s\n"
                    "        这条闸堵的就是「Γ-only 数据集静默进训练」：超胞 Γ-only ≡ 原胞"
                    " 3×3×3，\n"
                    "        实测参考构型 max|F| 会从 1e-4 涨到 0.103 eV/Å。\n"
                    "        处置：核对 step5_label 的 KPOINTS_GRID 后 rerun S5_label；\n"
                    "        确认要放行就把 KSPACING_TOL_BLOCK 设回 false。" % _msg)
            print("[OK] S5_label 帧 k 点间距 vs step1 参考：%s" % _msg)
        # ② INCAR/POTCAR 键级差异。★ ISMEAR/SIGMA 属【设计内合法差异】，不参与判定：
        #    step1 模板固定 ISMEAR=0，而 step5 走 decide_ismear(gap)（金属 → 1/0.2），
        #    拿它们当"不一致"会误杀所有金属体系的合法数据（2026-09-26 排查）。
        #    其余键不一致：默认 WARN，KSPACING_TOL_BLOCK=on → FAIL。
        _pa = dict(x.split("=", 1) for x in _bulk_fp[0][2].split("|") if "=" in x)
        _pb = dict(x.split("=", 1) for x in fp.split("|") if "=" in x)
        _diff = sorted(k for k in _pa if _pb.get(k) != _pa[k])
        _IGN = {"ISMEAR", "SIGMA"}
        _hard = [k for k in _diff if k not in _IGN]
        if _diff:
            print("[..] S5 帧 vs step1 参考的键差异：%s%s"
                  % (", ".join(_diff),
                     "（ISMEAR/SIGMA 为设计内合法差异，不判）" if set(_diff) & _IGN else ""))
        if _hard:
            _m2 = ("S5_label 帧与 step1 参考的 DFT 设置不一致，键 %s：帧 %s vs 参考 %s"
                   % (", ".join(_hard), {k: _pa[k] for k in _hard},
                      {k: _pb.get(k) for k in _hard}))
            if _blk:
                sys.exit("[ERROR] %s\n        处置：核对生成设置；确认要放行就把 "
                         "KSPACING_TOL_BLOCK 设回 false。" % _m2)
            print("[WARN] %s（KSPACING_TOL_BLOCK=false，放行）" % _m2)

    # ---- extend：并入 PRE_XYZ_FILES ----
    pre_paths = [x.strip() for x in (a.pre_xyz or "").split(",") if x.strip()]
    if pre_paths:
        pre_frames = load_pre_xyz(pre_paths, a.dim)
        kept, dropped = [], 0
        for fr in pre_frames:
            if fr["fingerprint"] is None:
                kept.append(fr)
                print("[WARN] %s（%s）没有 mlff_fingerprint，无法校验 DFT 设置，"
                      "默认接受（用户负责确认一致性）" % (fr["source"], fr["config_type"]))
            elif fr["fingerprint"] == fp:
                kept.append(fr)
            else:
                dropped += 1
        if dropped:
            print("[WARN] 指纹不一致的已有数据整批丢弃 %d 帧（不阻断流程）" % dropped)
        frames += kept

    # ---- 离群过滤 ----
    n_force_f, n_energy_f = filter_frames(
        [f for f in frames if f["source"] == "step5_label"],
        a.energy_limit, a.force_limit)
    n_filt = sum(1 for f in frames if f.get("filtered"))
    n_total = len(frames)
    frac = n_filt / max(n_total, 1)
    if n_force_f or n_energy_f:
        print("[..] 离群过滤分解：力上限滤 %d 帧，能量离群滤 %d 帧（总 %d/%d）"
              % (n_force_f, n_energy_f, n_filt, n_total))
    if frac > 0.3:
        sys.exit("[ERROR] 离群过滤比例 %.0f%% > 30%%：ENERGY_LIMIT/FORCE_LIMIT 过严，"
                 "有效数据被滤掉了。调大阈值后 retry 本步。" % (frac * 100))
    if frac > 0.1:
        print("[WARN] 离群过滤比例 %.0f%% > 10%%（阈值可能过严；被滤帧已标 filtered: true "
              "保留在 xyz 里，训练时排除）" % (frac * 100))

    # ---- 划分：固定测试集（id 哈希 %10），其余按 FPS 排序成训练集 ----
    import hashlib
    idx_nonfilt = [i for i, fr in enumerate(frames) if not fr["filtered"]]
    train_idx, test_idx, iso_idx = [], [], []
    for i in idx_nonfilt:
        fr = frames[i]
        if fr["config_type"] == "iso":
            iso_idx.append(i)                    # 孤立原子不参与 FPS（原子数不同），放训练集尾部
            continue
        if fr.get("id"):
            h = int(hashlib.md5(fr["id"].encode()).hexdigest(), 16)
            (test_idx if h % 10 == 0 else train_idx).append(i)
        else:
            train_idx.append(i)          # 无 id 的 PRE 帧全进训练集
    bulk = [i for i in idx_nonfilt if i not in iso_idx]
    feats = np.array([frame_features(frames[i]["E"], frames[i]["F"], len(frames[i]["F"]))
                      for i in bulk])
    order = fps_order(feats, a.fps_seed)          # bulk 内的位置
    train_set = set(train_idx)
    train_ordered = [bulk[i] for i in order if bulk[i] in train_set] + iso_idx

    # ---- 写 extxyz ----
    def write_xyz(path, idxs):
        from ase.io import write as ase_write
        outs = []
        for i in idxs:
            fr = frames[i]
            at = fr["atoms"].copy()
            info = {"config_type": fr["config_type"],
                    "REF_energy": float(fr["E"]),
                    "mlff_fingerprint": fp}
            if fr.get("strain_factor") is not None:
                info["strain_factor"] = float(fr["strain_factor"])
            if fr.get("volume_factor") is not None:
                info["volume_factor"] = float(fr["volume_factor"])
            if fr.get("strain_grun") is not None:
                info["strain_grun"] = float(fr["strain_grun"])
            if fr.get("rattle_std") is not None:
                info["rattle_std"] = float(fr["rattle_std"])
            if fr["filtered"]:
                info["filtered"] = True
            if fr["stress"] is not None and a.dim == "3d":
                info["REF_stress"] = np.array(fr["stress"], dtype=float)
            if a.dim == "2d":
                info["config_stress_weight"] = 0.0   # 2D 面外应力是垃圾，不训
            at.info = info
            # 力必须是 per-atom 数组（atoms.arrays）——mace 从 arrays_keys["forces"] 读，
            # 放进 info 里 mace 找不到，会静默退化成只训能量不训力（实测 RMSE_F=None、
            # 测试力 RMSE 322 meV/Å、声子 0.88 THz、体模量 0.1 GPa）。
            at.new_array("REF_forces", np.array(fr["F"], dtype=float))
            outs.append(at)
        ase_write(str(path), outs, format="extxyz")
        return len(outs)

    # all.xyz（含被过滤帧，供复查）与 train.xyz / test.xyz
    write_xyz(outdir / "all.xyz", list(range(n_total)))
    n_train = write_xyz(outdir / "train.xyz", train_ordered)
    n_test = write_xyz(outdir / "test.xyz", test_idx)

    # ---- e0s.json（孤立原子）----
    e0s = {}
    iso_man = [m for m in manifests]
    for mp in iso_man:
        man = json.loads(mp.read_text())
        for ent in man.get("iso_frames", []):
            cfgdir = cwd / a.step5_dir / ent["cfg_id"]
            if not outcar_done(cfgdir / "OUTCAR"):
                continue
            d = parse_outcar(cfgdir / "OUTCAR")
            from ase.data import chemical_symbols
            z = chemical_symbols.index(ent["element"])
            e0s[str(z)] = float(d["E"])
    if e0s:
        (outdir / "e0s.json").write_text(json.dumps(e0s, indent=2) + "\n",
                                         encoding="utf-8", newline="\n")
    else:
        print("[WARN] 没有孤立原子帧（extend 模式？）——e0s.json 留空，训练用基座 E0")

    # ---- 覆盖分析 + 三联图 ----
    report = coverage_analysis(frames, a.dim)
    plan = make_coverage_plan(frames, a.dim,
                              [float(x) for x in a.vol_factors.split(",") if x.strip()],
                              a.n_per_cell)
    report["sampling_plan"] = plan
    (outdir / "coverage_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8", newline="\n")

    try:
        _plot_parity(frames, train_ordered, test_idx, outdir / "energy_forces.png")
    except Exception as e:
        print("[WARN] energy_forces.png 画图失败：%s" % e)

    import hashlib as _hl
    _man_top = cwd / a.step4_dir / "struct_manifest.json"
    _dh = _hl.md5(_man_top.read_bytes()).hexdigest() if _man_top.is_file() else ""
    summary = {
        "generation": a.gen,
        "data_hash": _dh,
        "dim": a.dim,
        "n_frames_total": n_total,
        "n_train": n_train,
        "n_test": n_test,
        "n_filtered": n_filt,
        "filter_fraction": round(frac, 4),
        "energy_limit": a.energy_limit,
        "force_limit": a.force_limit,
        "fingerprint": fp,
        "vasp_versions": dict(sorted(vasp_versions.items())),   # {版本: 帧数}；fingerprint 不含版本，靠它发现换 VASP
        "kspacing_1_A": round(ksp, 4),
        "fps_seed": a.fps_seed,
        "e0s": e0s,
        "pre_xyz_files": pre_paths,
        "wall_time_s": round(time.time() - t0, 1),
    }
    (outdir / "dataset_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8", newline="\n")
    print("[DONE] dataset_summary.json：共 %d 帧（train %d / test %d / filtered %d），"
          "用时 %.1f s" % (n_total, n_train, n_test, n_filt, time.time() - t0))


CURRENT_FP = {}


def _plot_parity(frames, train_idx, test_idx, png_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    train_set = set(train_idx)
    test_set = set(test_idx)
    sets = {"训练集": [frames[i] for i in train_set],
            "测试集": [frames[i] for i in test_set],
            "过滤后": [f for f in frames if f["filtered"]]}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, (name, frs) in zip(axes, sets.items()):
        e = [f["E"] / len(f["F"]) for f in frs]
        f = [float(np.max(np.abs(f["F"]))) for f in frs]
        ax.scatter(e, f, s=8, alpha=0.6)
        ax.set_xlabel("E (eV/atom)")
        ax.set_ylabel("max|F| (eV/Å)")
        ax.set_title("%s (n=%d)" % (name, len(frs)))
    fig.tight_layout()
    fig.savefig(str(png_path), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
