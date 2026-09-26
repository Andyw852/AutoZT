#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step5_tbc.py -- (可选, GPU) 生成 GPUMD NEMD 界面热导输入（通用版）.

特点（相对旧版的关键修复）：
  * 热流由源/漏恒温器的累积传能直接给出（compute ... temperature 的最后两列），
    后处理不再依赖文献 kappa。
  * 通用：任意 extxyz 结构/元素、任意传输轴(x/y/z)、任意 GPUMD 势、pbc 可沿用结构文件。
  * deck 含 平衡 -> 烧入(heat_*) -> 测量 三段，测量段单独累积恒温器传能。

读 step.conf / 环境：
  TBC_STRUCTURE / TBC_NEP_MODEL / TBC_GPUMD_BIN / TBC_CUDA_VISIBLE_DEVICES
  TBC_AXIS (x|y|z) / TBC_PBC / TBC_INTERFACE_COORD / TBC_INTERFACE_Z(旧名)
  TBC_SOURCE_THICKNESS / TBC_SINK_THICKNESS / TBC_SOURCE_SIDE(high|low) / TBC_VACUUM_GAP_A
  TBC_T / TBC_DELTA_T / TBC_COUPLE / TBC_THERMOSTAT(heat_lan|heat_nhc|heat_bdp)
  TBC_TIME_STEP / TBC_SEED / TBC_BIN_WIDTH
  TBC_EQUIL_STEPS / TBC_BURN_STEPS / TBC_RUN_STEPS
  TBC_SAMPLE_INTERVAL / TBC_OUTPUT_INTERVAL / TBC_THERMO_INTERVAL
  TBC_KAPPA_A_W_mK / TBC_KAPPA_B_W_mK  (可选，仅后处理交叉核对)
产物: step5_tbc/{model.xyz, potential.txt, run.in, run_gpumd.sh, tbc_inputs.json}
"""
import json
import math
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stepconf  # noqa: E402

STEP = "step5_tbc"
OUTDIR = STEP

SPEC = {
    "TBC_STRUCTURE": ("", "str"),
    "TBC_NEP_MODEL": ("", "str"),
    "TBC_GPUMD_BIN": ("", "str"),
    "TBC_CUDA_VISIBLE_DEVICES": ("", "str"),
    "TBC_AXIS": ("z", "str"),
    "TBC_PBC": ("", "str"),
    "TBC_INTERFACE_COORD": (None, "float"),
    "TBC_INTERFACE_Z": (None, "float"),
    "TBC_SOURCE_THICKNESS": (5.0, "float"),
    "TBC_SINK_THICKNESS": (8.0, "float"),
    "TBC_SOURCE_SIDE": ("high", "str"),
    "TBC_VACUUM_GAP_A": (4.0, "float"),
    "TBC_T": (300.0, "float"),
    "TBC_DELTA_T": (40.0, "float"),
    "TBC_COUPLE": (100.0, "float"),
    "TBC_THERMOSTAT": ("heat_lan", "str"),
    "TBC_TIME_STEP": (1.0, "float"),
    "TBC_SEED": (0, "int"),
    "TBC_BIN_WIDTH": (2.0, "float"),
    "TBC_EQUIL_STEPS": (20000, "int"),
    "TBC_BURN_STEPS": (100000, "int"),
    "TBC_RUN_STEPS": (200000, "int"),
    "TBC_SAMPLE_INTERVAL": (10, "int"),
    "TBC_OUTPUT_INTERVAL": (100, "int"),
    "TBC_THERMO_INTERVAL": (1000, "int"),
    "TBC_KAPPA_A_W_mK": (None, "float"),
    "TBC_KAPPA_B_W_mK": (None, "float"),
}
AXES = {"x": 0, "y": 1, "z": 2}


def parse_extxyz(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    if len(lines) < 3:
        raise SystemExit("[ERROR] 结构文件内容不足: %s" % path)
    n = int(lines[0].split()[0])
    comment = lines[1]
    m = re.search(r'[Ll]attice\s*=\s*"([^"]+)"', comment)
    if not m:
        raise SystemExit('[ERROR] extxyz 第二行缺少 Lattice="...": %s' % path)
    lat = [float(x) for x in m.group(1).split()]
    if len(lat) != 9:
        raise SystemExit("[ERROR] Lattice 需要 9 个数（3x3 行主序）")
    mp = re.search(r'pbc\s*=\s*"([^"]+)"', comment)
    pbc = mp.group(1).strip() if mp else None
    mprop = re.search(r'[Pp]roperties\s*=\s*"?([A-Za-z0-9_:]+)"?', comment)
    props = mprop.group(1) if mprop else "species:S:1:pos:R:3"
    toks = re.findall(r'([A-Za-z_]+):([A-Za-z]):(\d+)', props)
    if not toks:
        toks = [("species", "S", "1"), ("pos", "R", "3")]
    layout, off = [], 0
    for name, typ, cnt in toks:
        layout.append((name, typ, int(cnt), off))
        off += int(cnt)
    col = {name: (typ, cnt, o) for name, typ, cnt, o in layout}
    if "species" not in col or "pos" not in col:
        raise SystemExit("[ERROR] extxyz Properties 需要含 species 与 pos: %s" % props)
    sym, pos = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        st, sc, so = col["species"]
        pt, pc, po = col["pos"]
        sym.append(p[so])
        pos.append(tuple(float(p[po + k]) for k in range(3)))
    if len(sym) != n:
        raise SystemExit("[ERROR] 原子数不符：声明 %d，实读 %d" % (n, len(sym)))
    return lat, pbc, sym, pos, props


def potential_elements(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().split()
    if head and head[0].lower().startswith("nep") and len(head) > 2:
        return set(head[2:])
    return None


def cross_area(lat, axis):
    a = lat[0:3]
    b = lat[3:6]
    c = lat[6:9]

    def cross(u, v):
        return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])

    def norm(u):
        return math.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2)

    if axis == 0:
        return norm(cross(b, c))
    if axis == 1:
        return norm(cross(c, a))
    return norm(cross(a, b))


def kappa_defaults(cwd, conf):
    ka, kb = conf["TBC_KAPPA_A_W_mK"], conf["TBC_KAPPA_B_W_mK"]
    p = os.path.join(cwd, "step2_props", "thermal_props.json")
    if (ka is None or kb is None) and os.path.isfile(p):
        with open(p, encoding="utf-8") as fh:
            d = json.load(fh)
        if ka is None:
            ka = d.get("channel", {}).get("ky_W_mK")
        if kb is None:
            kb = d.get("oxide", {}).get("kx_W_mK")
    return ka, kb


def auto_interface_by_composition(coords, species):
    """按元素组分沿轴的突变位置定界面；无法判定返回 None。

    用**原子窗**扫描，而不是固定 0.5 A bin + 正负 1 bin 平滑：当层间距 2~3.5 A 时，
    固定 bin 之间会出现整段空 bin，平滑后梯度被抹平，界面找不到（旧实现的失效区，
    恰好是二维材料/vdW 叠层的典型 3 A 层间距）。原子窗只看排序后连续的若干原子，
    与层间距无关。返回下层末原子与上层首原子的中点（界面键中心）。
    """
    uniq = sorted(set(species))
    if len(uniq) < 2:
        return None
    idx = {s: i for i, s in enumerate(uniq)}
    order = sorted(range(len(coords)), key=lambda i: coords[i])
    n = len(order)
    m = max(1, n // 4)
    c = [coords[i] for i in order]
    sidx = [idx[species[i]] for i in order]
    flo = sum(sidx[:m]) / m
    fhi = sum(sidx[-m:]) / m
    if abs(fhi - flo) < 0.5:
        return None
    w = max(3, n // 20)               # 原子窗宽
    if 2 * w >= n:
        w = max(1, n // 4)
    best = None
    for i in range(w, n - w + 1):
        left = sum(sidx[i - w:i]) / w
        right = sum(sidx[i:i + w]) / w
        j = abs(right - left)
        if best is None or j > best[0]:
            best = (j, i)
    if best is None or best[0] < 0.2 * abs(fhi - flo):
        return None
    i = best[1]
    if i <= 0 or i >= n:
        return None
    return 0.5 * (c[i - 1] + c[i])    # 界面键中心


def detect_interface(coords, vacuum_gap):
    zs = sorted(coords)
    gaps = [(zs[i + 1] - zs[i], 0.5 * (zs[i] + zs[i + 1])) for i in range(len(zs) - 1)]
    cand = [g for g in gaps if g[0] <= vacuum_gap]
    if cand:
        return max(cand)[1]
    return 0.5 * (zs[0] + zs[-1])


def align_coords(pos, lat, pbc):
    """把原子坐标对齐到 compute_chunk 的 bin 坐标系（原点为盒子原点）。

    周期轴：wrap 到 [0, L)。
    非周期轴：**不取模**。弛豫后的结构在非周期方向常有微小负坐标（如 -0.2 A），
    取模会把它挪到盒子另一端，既错位又可能被分进错误的源/漏组（真实案例：底层
    Ge 在 z=-0.2、pbc="T T F"，旧实现把它挪到 z=+69.8 并分进了源组）。这里改为
    整轴平移，使该轴最小值为 0，保持原子间相对构型不变。
    仅支持正交盒子（周期轴 wrap 按晶格对角元算）；只要存在周期轴且晶格有非对角分量
    就报错，避免静默按对角元 wrap 算错。
    返回 (new_pos, per_axis_is_periodic, warn_msgs)。
    """
    lens = [lat[0], lat[4], lat[8]]
    parts = pbc.split() if pbc else []
    is_T = [(len(parts) == 3 and parts[i].upper() == "T") for i in range(3)]
    if any(is_T):
        diag = max(abs(lat[0]), abs(lat[4]), abs(lat[8]), 1e-12)
        off = max(abs(lat[1]), abs(lat[2]), abs(lat[3]),
                  abs(lat[5]), abs(lat[6]), abs(lat[7]))
        if off > 1e-6 * diag:
            raise SystemExit(
                "[ERROR] 晶格非正交（非对角分量最大 %.4g A）：周期轴 wrap 只按对角元算，"
                "非正交盒子会静默算错。请先把结构转换/旋转到正交盒子，或自行 wrap 后把 "
                "pbc 设为全 F 交给本步只做整体平移。" % off)
    new = [[float(p[i]) for i in range(3)] for p in pos]
    for i in range(3):
        if is_T[i]:
            for q in new:
                q[i] = q[i] % lens[i]
    warn = []
    for i in range(3):
        if not is_T[i]:
            mn = min(q[i] for q in new)
            for q in new:
                q[i] -= mn
            mx = max(q[i] for q in new)
            if lens[i] > 0 and mx > lens[i] + 0.5:
                warn.append("[warn] %s 轴非周期，平移后跨度 %.2f A > 盒长 %.2f A；"
                            "compute_chunk 的 bin 可能不覆盖全部原子，建议把结构放进盒内"
                            % ("xyz"[i], mx, lens[i]))
    return [tuple(q) for q in new], is_T, warn


def main():
    conf = stepconf.load(SPEC, STEP, strict=False)
    cwd = os.getcwd()
    outdir = os.path.join(cwd, OUTDIR)
    os.makedirs(outdir, exist_ok=True)

    struct = (conf["TBC_STRUCTURE"] or "").strip()
    if not struct:
        sys.exit("[ERROR] 未提供 TBC_STRUCTURE（界面结构 extxyz）。见 GPU-TBC-RECIPE.md。")
    spath = struct if os.path.isabs(struct) else os.path.join(cwd, struct)
    if not os.path.isfile(spath):
        sys.exit("[ERROR] 找不到界面结构: %s" % spath)
    pot = os.path.expanduser((conf["TBC_NEP_MODEL"] or "").strip())
    if not pot or not os.path.isfile(pot):
        sys.exit("[ERROR] 找不到势文件 TBC_NEP_MODEL=%r" % pot)
    gbin = os.path.expanduser((conf["TBC_GPUMD_BIN"] or "").strip())
    if not gbin or not os.path.isfile(gbin):
        sys.exit("[ERROR] 找不到 gpumd 可执行文件 TBC_GPUMD_BIN=%r" % gbin)

    axis_name = (conf["TBC_AXIS"] or "z").strip().lower()
    if axis_name not in AXES:
        sys.exit("[ERROR] TBC_AXIS 只能是 x/y/z，收到 %r" % axis_name)
    ax = AXES[axis_name]

    lat, pbc_in, sym, pos, props = parse_extxyz(spath)
    els = potential_elements(pot)
    if els is not None:
        missing = sorted(set(sym) - els)
        if missing:
            sys.exit("[ERROR] 势不含元素 %s（势元素表: %s）" % (missing, sorted(els)[:15]))
    n = len(sym)
    # 坐标系对齐：compute_chunk 的 bin 坐标以盒子原点为 0 算；输入结构若没 wrap 到
    # [0, L)，两套坐标会错位。这里**只对周期轴 wrap**；非周期轴（通常是带真空的传输轴）
    # 不取模，改为整轴平移使最小值为 0——弛豫后出现 -0.2 A 这类负坐标很常见，取模会把
    # 它挪到盒子另一端并分错源/漏组。pbc 必须在 wrap 之前先定下来。
    _pbc_explicit = (conf["TBC_PBC"] or "").strip()
    if not _pbc_explicit and not pbc_in:
        sys.exit("[ERROR] 结构文件没有 pbc、TBC_PBC 也未给：不能默认 T T T（沿传输轴"
                 "周期会让源/漏隔着周期边界直接相邻，热流短路）。请显式设 "
                 "TBC_PBC='T T F'（传输轴取 F）或写进结构第二行。")
    pbc = _pbc_explicit or pbc_in
    pos, _pbc_T, _align_warn = align_coords(pos, lat, pbc)
    for _w in _align_warn:
        print(_w)
    _lens = [lat[0], lat[4], lat[8]]
    coord = [p[ax] for p in pos]
    cmin, cmax = min(coord), max(coord)
    span = cmax - cmin
    if span < 20.0:
        sys.exit("[ERROR] %s 向太薄（%.1f A），NEMD 需足够长度" % (axis_name, span))

    area = cross_area(lat, ax)
    vac_gap = float(conf["TBC_VACUUM_GAP_A"])
    cif = conf["TBC_INTERFACE_COORD"]
    if cif is None:
        cif = conf["TBC_INTERFACE_Z"]
    if cif is None:
        cif = auto_interface_by_composition(coord, sym)
        if cif is not None:
            auto = "composition"
        else:
            cif = detect_interface(coord, vac_gap)
            auto = "gap"
            print("[warn] 组分法判不出界面，回退最大间隙法（interface=%.2f A）：这只对"
                  "含真空层、且两侧同元素的堆叠可靠；共格异质界面请显式给 "
                  "TBC_INTERFACE_COORD。" % float(cif))
    else:
        auto = "explicit"
    cif = float(cif)

    src_th = float(conf["TBC_SOURCE_THICKNESS"])
    snk_th = float(conf["TBC_SINK_THICKNESS"])
    side = (conf["TBC_SOURCE_SIDE"] or "high").strip().lower()
    if side not in ("high", "low"):
        sys.exit("[ERROR] TBC_SOURCE_SIDE 只能是 high/low")
    groups = []
    for cc in coord:
        hi = cc >= (cmax - src_th)
        lo = cc <= (cmin + snk_th)
        if hi and cc > cif:
            groups.append(1 if side == "high" else 2)
        elif lo and cc < cif:
            groups.append(2 if side == "high" else 1)
        else:
            groups.append(0)
    ns, nk = groups.count(1), groups.count(2)
    if ns == 0 or nk == 0:
        sys.exit("[ERROR] 源/漏组为空（src=%d sink=%d）；请调 TBC_SOURCE/SINK_THICKNESS、"
                 "TBC_INTERFACE_COORD 或 TBC_SOURCE_SIDE" % (ns, nk))

    # 硬伤修复：沿传输轴若为周期性边界，源(高端)与漏(低端)隔着周期边界直接相邻，
    # 大量热流不经过界面直接从源流到漏（短路）→ q 被高估、G 虚高；"源/漏热流自洽"
    # 判据抓不到它（能量依然守恒）。只有界面本身是真空层(auto=="gap")时才可能打断
    # 这条短路；材料异质结(auto=="composition"/"explicit")在 pbc=T 下必然短路。
    _pbc_parts = pbc.split()
    _pbc_axis = _pbc_parts[ax].upper() if len(_pbc_parts) == 3 else "T"
    if _pbc_axis == "T" and auto != "gap":
        sys.exit("[ERROR] 沿 %s 轴 pbc=T 且界面非真空层(mode=%s)：源/漏隔着周期边界直接"
                 "相邻，热流不经界面短路，q/G 会被高估。请把 TBC_PBC 该轴设为 F（并在两端"
                 "加 fix 冻结原子），或改用含真空层的结构。" % (axis_name, auto))

    mp = os.path.join(outdir, "model.xyz")
    with open(mp, "w", encoding="utf-8") as fh:
        fh.write("%d\n" % n)
        fh.write('pbc="%s" Lattice="%s" Properties=species:S:1:pos:R:3:group:I:1\n'
                 % (pbc, " ".join("%.6f" % v for v in lat)))
        for s, (x, y, zz), g in zip(sym, pos, groups):
            fh.write("%-3s %12.6f %12.6f %12.6f %d\n" % (s, x, y, zz, g))
    shutil.copyfile(pot, os.path.join(outdir, "potential.txt"))

    T = float(conf["TBC_T"])
    dT = float(conf["TBC_DELTA_T"])
    coup = float(conf["TBC_COUPLE"])
    dt = float(conf["TBC_TIME_STEP"])
    binw = float(conf["TBC_BIN_WIDTH"])
    neq = int(conf["TBC_EQUIL_STEPS"])
    nburn = int(conf["TBC_BURN_STEPS"])
    nrun = int(conf["TBC_RUN_STEPS"])
    si = int(conf["TBC_SAMPLE_INTERVAL"])
    oi = int(conf["TBC_OUTPUT_INTERVAL"])
    ti = int(conf["TBC_THERMO_INTERVAL"])
    seed = int(conf["TBC_SEED"])
    thermo = (conf["TBC_THERMOSTAT"] or "heat_lan").strip().lower()
    if thermo not in ("heat_lan", "heat_nhc", "heat_bdp"):
        sys.exit("[ERROR] TBC_THERMOSTAT 只能是 heat_lan/heat_nhc/heat_bdp")

    L = []
    L.append("potential potential.txt")
    L.append("velocity %g%s" % (T, (" seed %d" % seed) if seed else ""))
    L.append("time_step %g" % dt)
    if neq > 0:
        L.append("ensemble nvt_nhc %g %g %g" % (T, T, coup))
        L.append("dump_thermo %d" % ti)
        L.append("run %d" % neq)
    # 烧入段：heat_* 恒温器建立温差（暂态，不测量；compute_* 留到测量段再声明）
    L.append("ensemble %s %g %g %g 1 2" % (thermo, T, coup, dT))
    L.append("dump_thermo %d" % ti)
    if nburn > 0:
        L.append("run %d" % nburn)
    # 测量段：GPUMD 每个 run 结束都会清掉 ensemble 与 compute 动作
    # （Integrate::finalize 的 ensemble_.reset()、Measure::post_run 的 actions_.clear()），
    # 所以必须重新声明 ensemble + 两个 compute 再 run；否则测量段会报
    # "An ensemble must be specified before each run." 或写不出 compute.out / compute_chunk.out。
    L.append("ensemble %s %g %g %g 1 2" % (thermo, T, coup, dT))
    L.append("compute 0 %d %d temperature" % (si, oi))
    L.append("compute_chunk %d %d bin/1d %s lower %g temperature density/number"
             % (si, oi, axis_name, binw))
    L.append("dump_thermo %d" % ti)
    L.append("run %d" % nrun)
    with open(os.path.join(outdir, "run.in"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    dev = (conf["TBC_CUDA_VISIBLE_DEVICES"] or "").strip()
    with open(os.path.join(outdir, "run_gpumd.sh"), "w", encoding="utf-8") as fh:
        fh.write("#!/bin/bash\n# 由 gen_step5_tbc.py 生成：运行 GPUMD NEMD\n")
        fh.write("set -e\ncd \"$(dirname \"$0\")\"\n")
        if dev:
            fh.write("export CUDA_VISIBLE_DEVICES=\"%s\"\n" % dev)
        fh.write("GPUMD=\"%s\"\n" % gbin)
        fh.write("[ -x \"$GPUMD\" ] || { echo \"gpumd 不可执行: $GPUMD\" >&2; exit 1; }\n")
        fh.write("echo \"GPU env: $CUDA_VISIBLE_DEVICES\"\n")
        fh.write("nvidia-smi --query-gpu=index,memory.used --format=csv,noheader 2>/dev/null | head -1 || true\n")
        fh.write("\"$GPUMD\" > gpumd.log 2>&1\n")
        fh.write("echo \"GPUMD done. 输出: compute.out compute_chunk.out thermo.out (见 gpumd.log)\"\n")
    os.chmod(os.path.join(outdir, "run_gpumd.sh"), 0o755)

    ka, kb = kappa_defaults(cwd, conf)
    meta = {"TBC_INPUTS_DONE": True, "n_atoms": n, "axis": axis_name,
            "interface_auto": auto, "interface_coord_A": cif,
            "coord_min_A": cmin, "coord_max_A": cmax,
            "n_source": ns, "n_sink": nk,
            "source_side": side, "source_thickness_A": src_th, "sink_thickness_A": snk_th,
            "vacuum_gap_A": vac_gap, "bin_width_A": binw, "area_A2": area, "pbc": pbc,
            "T_K": T, "delta_T_K": dT, "couple": coup, "thermostat": thermo,
            "time_step_fs": dt, "seed": seed,
            "sample_interval": si, "output_interval": oi, "thermo_interval": ti,
            "equil_steps": neq, "burn_steps": nburn, "run_steps": nrun,
            "kappa_a_W_mK": ka, "kappa_b_W_mK": kb,
            "structure": spath, "potential": pot, "gpumd_bin": gbin,
            "cuda_visible_devices": dev, "properties_in": props}
    with open(os.path.join(outdir, "tbc_inputs.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)

    print("[OK] %s  (%d atoms; axis=%s; src=%d sink=%d; interface %s=%.2f A, mode=%s)"
          % (os.path.join(outdir, "tbc_inputs.json"), n, axis_name, ns, nk,
             axis_name, cif, auto))
    print("    area=%.1f A^2 ; pbc=%s ; thermostat=%s" % (area, pbc, thermo))
    print("    run: bash %s" % os.path.join(OUTDIR, "run_gpumd.sh"))


if __name__ == "__main__":
    main()
