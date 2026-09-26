#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""overlap_preflight.py —— 重叠路径的运行前检查（VERIFICATION V24）。

把本次裁决中靠手工做的核对固化成自动检查，任何 2D 材料在 S8.4 生成输入时都会跑一遍。

【背景】裁决结论（V22/V23）：真实波函数重叠必须让 h5 走 AMSET 的 from_data 路径
（h5 的 k 点自身构成完整网格），否则 AMSET 会对系数做去对称化，实测把 ADP 迁移率
抬高 10~30 倍。unity_overlap 不受影响（它根本不读波函数）。

【在哪跑】（2026-09-18 补齐，V25.10）
  gen 时：cwd = 材料 result 目录。本地只有生成的输入、没有 vasprun.xml / wavefunction.h5，
          所以 ②③④ 会「跳过」——gen 时真正起作用的是 ⓪①。
  作业内（--in-job）：cwd = AMSET 运行目录（…/result/step8_amset 或 step8.4_amset2d），
          这里 vasprun.xml / wavefunction.h5 都是有效软链，②③④ 才真正生效。
          gen_step10_amset.py / gen_step14_amset2d.py 会把本文件复制进运行目录，
          并在 amset 命令之前调用它（失败即 exit 1，不出不可信的数）。

六项检查：
  ⓪ 2D + 非 unity -> **直接拦截**。二维必须 unity_overlap: true（V23 裁决）。
     这条不需要任何波函数数据，gen 时与作业内都生效 —— 也是能挡住原始 bug 的那条。
  ① 版本 + 重叠模式：<0.5.1 且真实重叠时告警（该版本缺 G 平移，见 V21）。
  ② h5 完整性：k 点数是否等于网格总点数（= 会走 from_data，不去对称化）。
     真实重叠 + 非完整网格：**2D 拦截**；3D 只告警（BLOCK_3D_REAL_OVERLAP=False，
     等 Si 全网格对照出结果再定，V25）。
  ③ S3 与 S3b 一致性：NBANDS 相同、共同 k 点上本征值偏差 < 1 meV（只看 AMSET 实际
     使用的带窗口，默认 11-17；最高的几条带没收敛会给出几十 meV，属正常）。
  ④ amset wave 的带窗口必须是显式的（自动选会得到不同带区间 -> 带索引错位，
     系数没法逐带比较）。
  ⑤ 弹性张量正定性：Christoffel 最小特征值 <= 0 时 **拦截** —— 与重叠无关的独立坑，
     AMSET 的形变势因子用负特征值 -> ADP 直接算废（V25.10 §九；Mo2S3/CrS2 实测）。
     只有运行目录里有 settings.yaml 时才生效（gen 写完 settings 之后 / 作业内都满足）。
"""
import glob
import os
import re
import sys
from pathlib import Path

# 兜底带窗口（1-based，与 vasprun/AMSET 日志同一编号）：读不到运行日志时用。
#   真实窗口优先从 AMSET 运行日志的 "Interpolating spin-up bands L-H" 里读
#   （2026-09-18 实测：MoS2 全网格 = 11-17、CrSe2_Hex = 9-15、Pb2Bi2Te5 = 36-51），
#   各材料不同 —— 钉死一个数会把 3D 材料误判（V25.10 修正）。
BAND_LO, BAND_HI = 11, 17
EIG_TOL_MEV = 1.0
# 3D 真实重叠 + 非完整网格：False = 只告警（当前），True = 与 2D 一样拦截。
# 为什么 3D 不拦：三维"是否真算错"要等 Si 的全网格对照（V25），而且 Si 有 6 个
# 等价谷，改 unity 最多低估约 6 倍 —— 结论出来之前不擅自改变三维的行为。
BLOCK_3D_REAL_OVERLAP = False
# 2D 判据：运行目录里有 2d_correction.json（gen 落盘的二维修正依据），
# 或 step1 的 workflow_method.txt 里 DIM=2d。
DIM_MARKERS = ("2d_correction.json",)


def resolve_layout(path):
    """返回 (base, out)：base = 材料 result 目录，out = AMSET 运行目录。

    · gen 时传 result 目录（里面没有 settings.yaml，名字也不是 step8*）-> base=它自己。
    · --in-job 时传运行目录（里面有 settings.yaml）-> base=上一级。
    """
    p = Path(path).resolve()
    if (p / "settings.yaml").is_file() or p.name.startswith("step8"):
        return p.parent, p
    out = None
    for name in ("step8.4_amset2d", "step8_amset"):
        if (p / name / "settings.yaml").is_file():
            out = p / name
            break
    return p, out


def is_2d_run(base, out):
    """是不是二维体系（用 gen 落盘的证据判断，两处任一命中即为 2D）。"""
    if out is not None:
        for m in DIM_MARKERS:
            if (out / m).is_file():
                return True
    for d in ("step1_opt", "step1"):
        f = base / d / "workflow_method.txt"
        if f.is_file():
            for ln in f.read_text(errors="ignore").splitlines():
                if ln.strip().upper().startswith("DIM="):
                    return ln.split("=", 1)[1].strip().lower() == "2d"
    return False


def h5_summary(h5_path):
    """返回 dict(nk, mesh, complete, ng, nb)；读不到返回 {"error": ...}。"""
    try:
        import h5py
        import numpy as np
        from amset.electronic_structure.kpoints import get_mesh_from_kpoint_numbers
        with h5py.File(h5_path, "r") as f:
            kp = np.array(f["kpoints"])
            ng = f["gpoints"].shape[0]
            nb = f["coefficients_up"].shape[0]
        mesh = get_mesh_from_kpoint_numbers(kp)
        prod = int(np.prod(mesh))
        return {"nk": len(kp), "mesh": tuple(int(x) for x in mesh), "mesh_prod": prod,
                "complete": prod == len(kp), "ng": ng, "nb": nb}
    except Exception as e:
        return {"error": "%s: %s" % (type(e).__name__, e)}


def elastic_min_christoffel(elastic):
    """弹性张量的最小 Christoffel 特征值（AMSET 同款），异常返回 None。

    AMSET 形变势散射因子 = deform²/c（各向异性为 Σ(strain·deform)²/c），c 是
    Christoffel 方程特征值（elastic.py:180-201）。弹性张量非正定 -> c 可为负
    -> 散射率/迁移率算废：Mo2S3 C44=-40.16 时 ADP 电子 6946，CrS2 C44=-0.52 时空穴 2e4 或 0。
    二维 slab 在真空方向没有层间剪切，负的 C44 基本都是 slab 归一化伪影（V25.10 §九）。
    """
    try:
        import numpy as np
        from pymatgen.core.tensors import Tensor
        from amset.scattering.elastic import (get_christoffel_tensors,
                                              solve_christoffel_equation)
    except Exception:
        return None
    try:
        C = np.array(Tensor.from_voigt(np.array(elastic, dtype=float)))
    except Exception:
        return None
    dirs = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1),
            (1, 1, 1), (1, -1, 0), (1, 0, -1), (0, 1, -1)]
    q = np.array(dirs, dtype=float)
    q = q / np.linalg.norm(q, axis=1)[:, None]
    try:
        ev = np.array(solve_christoffel_equation(get_christoffel_tensors(C, q))[0],
                      dtype=float)
    except Exception:
        return None
    n = int(np.argmin(ev))
    idx = np.unravel_index(n, ev.shape)
    return float(ev.ravel()[n]), dirs[idx[0]] if ev.ndim > 1 else dirs[0]


def inplane_eig_ratio(elastic):
    """二维**面内** Christoffel 特征值的 (最小/绝对最大) 比；异常返回 None。

    只看面内方向（q 的 z=0），因为二维 ADP/POP 的面内散射由面内 2x2 子块决定。
    ★ 典型错位：张量仍是 VASP 打印顺序 (XX YY ZZ XY YZ ZX)、没重排成标准 Voigt
    (XX YY ZZ YZ XZ XY) —— 面内剪切 C66 槽里放的是接近 0 的 C_zxzx，比值 << 1e-2。
    实测：SS 旧 settings 的 C66 槽 0.557 / C44 槽 56.337 -> 比值约 1e-2（拦下）；
    重排后 C66=56.337 -> 正常。
    """
    try:
        import numpy as np
        from pymatgen.core.tensors import Tensor
        from amset.scattering.elastic import (get_christoffel_tensors,
                                              solve_christoffel_equation)
    except Exception:
        return None
    try:
        C = np.array(Tensor.from_voigt(np.array(elastic, dtype=float)))
    except Exception:
        return None
    ang = np.linspace(0.0, np.pi, 9)[:-1]
    q = np.array([[np.cos(a), np.sin(a), 0.0] for a in ang])
    try:
        # ★ 只取**面内 2x2 子块**（与 2D 插件一致）：slab 的面外剪切 C44 本来就小，
        #   把它算进来会把所有二维材料都误判。
        g = get_christoffel_tensors(C, q)[:, :2, :2]
        ev = np.linalg.eigvalsh(g)
    except Exception:
        return None
    return float(ev.min() / (np.abs(ev).max() + 1e-30))


def _parse_vasprun(p):
    t = Path(p).read_text(errors="ignore")
    i = t.find('name="kpointlist"')
    seg = t[i:t.find("</varray>", i)]
    ks = [tuple(float(x) for x in m.groups())
          for m in re.finditer(r"<v>\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s*</v>", seg)]
    j = t.find("<eigenvalues>")
    j2 = t.find("</eigenvalues>", j)
    blocks = re.findall(r'<set\s+comment="kpoint[^"]*"\s*>(.*?)</set>', t[j:j2], re.S)
    eig = [[float(m.group(1)) for m in re.finditer(r"<r>\s*([-\d.Ee+]+)\s+[-\d.Ee+]+\s*</r>", b)]
           for b in blocks]
    return ks, eig


def check_s3_vs_s3b(s3_dir, s3b_dir, band_lo=BAND_LO, band_hi=BAND_HI):
    a, b = Path(s3_dir) / "vasprun.xml", Path(s3b_dir) / "vasprun.xml"
    if not (a.is_file() and b.is_file()):
        return None
    k1, e1 = _parse_vasprun(a)
    k2, e2 = _parse_vasprun(b)
    if not e1 or not e2:
        return {"error": "本征值解析失败"}
    n1, n2 = len(e1[0]), len(e2[0])
    idx = {tuple(round(v, 6) for v in kk): n for n, kk in enumerate(k2)}
    dmax, common = 0.0, 0
    for n in range(min(len(k1), len(e1))):
        m = idx.get(tuple(round(v, 6) for v in k1[n]))
        if m is None or m >= len(e2):
            continue
        common += 1
        for bi in range(max(band_lo - 1, 0), min(band_hi, n1, n2)):
            dmax = max(dmax, abs(e1[n][bi] - e2[m][bi]) * 1000.0)
    return {"nbands_s3": n1, "nbands_s3b": n2, "common": common, "max_dev_mev": dmax,
            "band_window": (band_lo, band_hi),
            "ok": (n1 == n2 and common > 0 and dmax < EIG_TOL_MEV)}


def run_band_window(out_dir, default=(BAND_LO, BAND_HI)):
    """从 AMSET 运行日志读它**实际插值**的带窗口（1-based，含两端）。

    AMSET 每次 run 都会打印 "Interpolating spin-up bands L-H" —— 这就是它从 h5 里
    取的带区间，也是 S3/S3b 该比对的那几条带。读不到就退回 default。
    """
    p = Path(out_dir) / "amset.log"
    if not p.is_file():
        return None          # 首次运行必然没有 amset.log —— 必须判"未知"，不能拿默认值冒充
    m = re.search(r"Interpolating spin-(?:up|down) bands\s+(\d+)\s*[-\u2014\u2013~]\s*(\d+)",
                  p.read_text(errors="ignore"))
    return (int(m.group(1)), int(m.group(2))) if m else None


def band_windows(log_paths):
    """读 amset wave 日志里的 "Including bands L-H" -> [(路径, lo, hi)]（去重保序）。"""
    out = []
    for p in log_paths:
        if not os.path.isfile(p):
            continue
        seen = set()
        t = Path(p).read_text(errors="ignore")
        for lo, hi in re.findall(r"Including bands\s+(\d+)[^\d]+(\d+)", t):
            if (lo, hi) in seen:
                continue
            seen.add((lo, hi))
            out.append((p, int(lo), int(hi)))
    return out


def run(cwd, out_dir=None, unity_overlap=False):
    two_d = False
    _passed = Path(cwd)
    base, out = resolve_layout(_passed)
    if out_dir is not None and not (_passed / "settings.yaml").is_file() \
            and not _passed.name.startswith("step8"):
        out = Path(out_dir)              # gen 时显式传进来的运行目录
    if out is None:
        return "warn", ["  0) 找不到 AMSET 运行目录（settings.yaml）—— 跳过（本步还没生成？）"]
    cwd, out_dir = base, out
    two_d = is_2d_run(cwd, out_dir)
    lines, err, warn = [], False, False
    # step.conf 显式覆盖 UNITY_OVERLAP 时，gen 会在 settings.yaml 里写
    # AZ_OVERLAP_CONTROLLED=1 —— ⓪/② 由"拦截"降为"告警"。
    # ★ 2026-09-20 澄清（用户）：这个标记只是**放行开关**，不改变结果的物理性质。
    #   若 h5 是**完整网格（走 from_data）**，这次 run 就是**生产数据**，不是"只能算比值的对照"；
    #   "受控对照"的定位只适用于**不完整 h5** 的真实重叠。⓪ 在 h5 检查之前执行、拿不到完整性，
    #   所以这里措辞偏保守——看到"不作生产结果"那句时，先看 ② 的 h5 是否完整：
    #   完整网格 + from_data 的真实重叠 = 可信生产数据（V23 的残差已解释：谷内通道也被压低，
    #   实测均值约 0.9，且 N_ch=3.13 > ADP 比值 2.52）。
    #   TODO（不急）：把 ⓪/② 的放行判据改成"完整网格即放行"，不再借用 AZ_OVERLAP_CONTROLLED。
    _st = Path(out_dir) / "settings.yaml"
    controlled = bool(_st.is_file() and "AZ_OVERLAP_CONTROLLED=1" in _st.read_text(errors="ignore"))

    # ---- ⓪ 二维默认 = 全网格 + 真实重叠（2026-09-26 用户决定，结论翻转）----
    #   旧逻辑（2026-09-17 V22/V23）要求"2D 必须 unity_overlap: true"，理由是真实重叠
    #   会被去对称化算坏。新证据把因果关系倒过来了：
    #     · unity 才是错的那一个 —— 它把谷间/大 q 散射当成完全耦合，实测 K→K' 真实
    #       |I|^2 ≈ 0.25，于是迁移率系统性**偏低**（只能看数量级，Seebeck 基本可用）；
    #     · 去对称化确实坏，但那是因为它本身就坏（无反演体系 38.9% 的 k 点错，
    #       TR 路径 5.2%；Si 对照 97.7%）—— 修法是改走全网格 from_data，而不是退回 unity。
    #   证据：tmp/amset2d/DESYM_FINDINGS.md、upstream_issue_tr_tau.md。
    lines.append("  0) 维度：%s；unity_overlap = %s"
                 % ("2D" if two_d else "3D/未知", unity_overlap))
    if two_d and unity_overlap is True:
        warn = True
        lines.append("     [WARN] 二维 unity_overlap=true：**仅作快速筛选**。"
                     "unity 把谷间/大 q 散射当成完全耦合（K->K' 真实 |I|^2≈0.25，"
                     "全网格实测），迁移率系统性**偏低**，只能看数量级；Seebeck 基本可用，"
                     "形变势(DPT)不受影响。正式结果请用 UNITY_OVERLAP=false + "
                     "WAVEFUNCTION_FULL=true（S3b ISYM=-1 -> S4b 全网格 -> 真实重叠）。")
    elif two_d:
        lines.append("     二维真实重叠：正确路径（要求全网格 h5，见下面第 2 项核对）。")

    h5s = [out_dir / "wavefunction.h5", cwd / "step4_wave" / "wavefunction.h5",
           cwd / "step4b_wave_full" / "wavefunction.h5"]
    info = None
    for h5 in h5s:
        if Path(h5).is_file():
            info = h5_summary(h5)
            info["path"] = str(h5)
            break
    nb = info.get("nb") if isinstance(info, dict) else None
    if info is None:
        lines.append("  2) h5：找不到 wavefunction.h5，跳过")
    elif "error" in info:
        lines.append("  2) h5：读取失败（%s）" % info["error"])
    else:
        complete = info["complete"]
        lines.append("  2) h5：nk=%d  网格=%s (总 %d)  -> %s   [G=%d, 带=%d]"
                     % (info["nk"], "x".join(map(str, info["mesh"])), info["mesh_prod"],
                        "完整网格，走 from_data（不去对称化）" if complete
                        else "非完整网格 -> AMSET 会去对称化",
                        info["ng"], info["nb"]))
        if not complete and not unity_overlap and two_d and controlled:
            warn = True
            lines.append("     [WARN] 2D + 真实重叠 + 非完整网格：受控对照放行（只用于算比值，"
                         "不作生产结果）。")
        elif not complete and not unity_overlap and two_d:
            err = True
            lines.append("     ★ 拦截：2D + 真实重叠 + **非完整网格**。真实重叠只有走完整网格的 "
                         "from_data 才对；非完整网格会触发 AMSET 的去对称化，"
                         "而它在本体系上会大面积算错。")
            lines.append("     正确做法：打开 wavefunction_full 分支（S3b 用 ISYM=-1 出全网格 "
                         "WAVECAR -> S4b 出全网格 wavefunction.h5）。")
            lines.append("     2026-09-26 量化（tmp/amset2d/DESYM_FINDINGS.md）：逐点真值裁决下 "
                         "MoS2 只有 38.9% 的 k 点正确，TR（R->-R）路径 5.2%，24 个操作里 15 个全错；"
                         "同一脚本 Si 对照 97.7%（Si 有反演，TR 被点群吸收、无此类操作）。"
                         "根因是 TR 分支 tau 被反号两次（symmetry.py:187-188 + common.py:123-124）。")
        elif not complete and not unity_overlap:
            warn = True
            lines.append("     [WARN] 3D + 真实重叠 + 非完整网格 -> AMSET 会对系数去对称化，"
                         "二维实测 ADP 被抬高 10~30 倍（V23）；三维 Si 全网格对照（V26）："
                         "ADP/overall 只差 3~9%、IMP 逐机制高 1.3~1.7 倍 -> 结论值可用，逐机制值需注明。"
                         "当前 BLOCK_3D_REAL_OVERLAP=%s -> 只告警不拦。"
                         % BLOCK_3D_REAL_OVERLAP)
            lines.append("     ★ 2026-09-19 补充（VERIFICATION V33）：MoS2 二维偏差的**主因**已定位到"
                         "去对称化的 **时间反演 x 旋转复合分支** —— 非完整网格里约 44% 的 TRxR 格点，"
                         "其系数子空间与直接 DFT 不符（100% 错），纯 TR 与非 TR 基本正确。"
                         "**凡空间群非中心对称（必须用时间反演补满网格）的材料，真实重叠 + 非完整网格"
                         "都不可信**；首选 S3b 用 ISYM=-1 出全网格波函数（走 from_data）。")

    ver = None
    try:
        import importlib.metadata as md
        ver = md.version("amset")
    except Exception:
        try:
            import amset
            ver = getattr(amset, "__version__", None)
        except Exception:
            ver = None
    if ver:
        lines.append("  1) AMSET 版本 = %s；unity_overlap = %s" % (ver, unity_overlap))
        if not unity_overlap and ver < "0.5.1":
            warn = True
            lines.append("     [WARN] <0.5.1 且用真实重叠：该版本缓存重叠路径缺 G 平移（V21）。"
                         "只有 h5 为完整网格时才可接受 —— 见上面第 2 条。")
    else:
        lines.append("  1) AMSET 版本：本机读不到（登录节点可能没装），跳过版本检查")

    # 窗口优先取 AMSET 运行日志里它实际插值的带区间。
    # ★ 读不到（首次运行/尚未跑 amset run）时必须判【未知】并跳过窗口一致性比较 ——
    #   拿兜底默认 11-17 去比会把任何窗口 != 11-17 的材料误拦
    #   （2026-09-18 Si 实测：step4_wave 真实窗口 2-6，被 11-17 误判"不一致"直接拦死 gen）。
    win = run_band_window(out_dir)
    win_known = win is not None
    if win_known:
        lines.append("  3/4) 本次 AMSET 实际插值的带窗口 = %d-%d" % win)
    else:
        lines.append("  3/4) 本次 AMSET 实际插值的带窗口：**运行日志缺失**（首次运行/尚未跑）"
                     " -> 跳过窗口一致性比较；等 amset.log 出来后再判定")
        win = (BAND_LO, BAND_HI)      # 仅作下面 S3/S3b 比对的默认口径
    r = check_s3_vs_s3b(cwd / "step3_uniform", cwd / "step3b_uniform_full",
                        band_lo=win[0], band_hi=win[1])
    if r is None:
        lines.append("  3) S3/S3b：没有 step3b_uniform_full，跳过（只在 wavefunction_full 分支下有意义）")
    elif "error" in r:
        lines.append("  3) S3/S3b：%s" % r["error"])
    else:
        lines.append("  3) S3/S3b：NBANDS %d vs %d；共同 k 点 %d；带 %d-%d 内最大偏差 %.4f meV  -> %s"
                     % (r["nbands_s3"], r["nbands_s3b"], r["common"],
                        r["band_window"][0], r["band_window"][1], r["max_dev_mev"],
                        "通过" if r["ok"] else "**不通过**"))
        if not r["ok"]:
            err = True
            lines.append("     ★ 拦截：两者能带不一致（NBANDS 不同或偏差 >= %.1f meV），"
                         "重叠与能带不能配套使用。" % EIG_TOL_MEV)

    logs = glob.glob(str(cwd / "step4_wave" / "amset.log")) + \
           glob.glob(str(cwd / "step4b_wave_full" / "amset.log"))
    wins = band_windows(logs)
    if wins and not win_known:
        # AMSET 运行日志缺失 -> 无从比较，只报出 wave 侧窗口（不当成"不一致"拦下）
        lines.append("  4) amset wave 带窗口：%s（AMSET 运行日志缺失，暂不判一致性）"
                     % "；".join("%s = %d-%d"
                                 % (os.path.basename(os.path.dirname(p)), lo, hi)
                                 for p, lo, hi in wins))
    elif wins:
        lines.append("  4) amset wave 带窗口：%s"
                     % "；".join("%s = %d-%d"
                                 % (os.path.basename(os.path.dirname(p)), lo, hi)
                                 for p, lo, hi in wins))
        bad = [(p, lo, hi) for p, lo, hi in wins if (lo, hi) != win]
        if bad:
            err = True
            for p, lo, hi in bad:
                lines.append("     ★ 拦截：%s 的带窗口 %d-%d 与 AMSET 实际插值的 %d-%d 不一致"
                             " -> 系数带索引错位，两者不能配套使用。"
                             "用同一个显式窗口 amset wave --bands %d:%d 重跑。"
                             % (os.path.basename(os.path.dirname(p)), lo, hi, win[0], win[1],
                                win[0], win[1]))
        elif nb is not None and wins[0][2] - wins[0][1] + 1 != nb:
            warn = True
            lines.append("     [WARN] h5 存了 %d 条带，而带窗口 %d-%d 是 %d 条 —— 对不上，"
                         "核对 amset wave --bands。" % (nb, wins[0][1], wins[0][2],
                                                        wins[0][2] - wins[0][1] + 1))
    elif logs:
        lines.append("  4) amset wave 带窗口：日志里没有 'Including bands' 行（老日志/未跑），跳过")

    # 5) 弹性张量正定性 —— 与重叠无关的独立坑（V25.10 §九）：非正定 -> ADP 直接算废
    st_path = out_dir / "settings.yaml"
    if st_path.is_file():
        try:
            import yaml
            _s = yaml.safe_load(st_path.read_text(errors="ignore")) or {}
        except Exception:
            _s = None
        if _s is None:
            lines.append("  5) 弹性张量：读 settings.yaml 失败（缺 yaml？），跳过")
        else:
            _el = _s.get("elastic_constant")
            me = elastic_min_christoffel(_el) if _el else None
            if not _el:
                lines.append("  5) 弹性张量：settings.yaml 里没有 elastic_constant，跳过")
            elif me is None:
                lines.append("  5) 弹性张量：装备缺 pymatgen/amset 或格式异常，跳过")
            else:
                val, d = me
                _r2 = inplane_eig_ratio(_el) if two_d else None
                lines.append("  5) 弹性张量：Christoffel 最小特征值 = %.4f（方向 %s）%s"
                             % (val, d, "" if _r2 is None else "；面内特征值比 = %.4f" % _r2))
                if two_d and _r2 is not None and (_r2 <= 0 or _r2 < 1e-2):
                    err = True
                    lines.append("     ★ 拦截（**面内**不可信）：二维面内 Christoffel 特征值比 = %.4f"
                                 "（<=0 或 <1e-2）—— 要么张量仍是 VASP 打印顺序没重排（C44/C66 互换；"
                                 "实测 SS 旧 settings C66 槽 0.557 / C44 槽 56.337），要么面内本身非正定。"
                                 "用当前 gen（gen_step10/gen_step14）重新生成 settings 后再跑。" % _r2)
                elif val <= 0:
                    warn = True
                    lines.append("     [WARN] 最小 Christoffel 特征值为负 —— 多为**面外**剪切（二维 slab "
                                 "的真空伪影），不是面内问题。gen 现已对二维面外剪切置零（原值记在 "
                                 "2d_correction.json 的 elastic_outofplane_shear_zeroed）；"
                                 "旧 settings 请 retry 重生成，或走 step8.4_amset2d 插件（面内 2x2）。")

    return ("error" if err else ("warn" if warn else "ok")), lines


def main():
    import argparse
    ap = argparse.ArgumentParser(description="AMSET 重叠路径运行前检查（VERIFICATION V24）")
    ap.add_argument("path", nargs="?", default=".",
                    help="gen 时传材料 result 目录；--in-job 时传 AMSET 运行目录")
    ap.add_argument("--in-job", action="store_true",
                    help="在 AMSET 作业里运行（cwd = 运行目录，h5/vasprun.xml 都在）")
    args = ap.parse_args()
    cwd = Path(args.path)
    base, out_dir = resolve_layout(cwd)
    unity = True
    if out_dir is None:
        out_dir = cwd / "step8.4_amset2d"
    st = out_dir / "settings.yaml"
    if st.is_file():
        m = re.search(r"^unity_overlap:\s*(\S+)", st.read_text(errors="ignore"), re.M)
        unity = (m.group(1).lower() == "true") if m else False
    v, lines = run(cwd, out_dir, unity)
    print("重叠路径运行前检查（VERIFICATION V24%s）：" % ("，作业内 --in-job" if args.in_job else ""))
    for l in lines:
        print(l)
    print("  结论：%s" % {"ok": "通过", "warn": "通过（有告警）", "error": "**拦截**"}[v])
    return 1 if v == "error" else 0


if __name__ == "__main__":
    sys.exit(main())
