#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step20_zt.py —— ZT 全流程汇总（zt-dft-cpu 的收尾步，run: gen）。

输入（cwd = <超算 work_dir>/<材料>/zt-dft-cpu/；autozt 按 skill.yaml 的 needs_results
从【本材料的本地目录】选出上游技能已拉回的结果，推到 inputs/<键>/，附 source.json）：
  · inputs/transport/transport.json      ke-dft-cpu：3D 取 S8_kappa，2D 取 S8.4_amset2d：S/σ/κ_e
  · inputs/kappa_L/kappa_summary.json    kl-dft-cpu S6_kappa（或 kl-mlff-* / fit-fc-thermal）：κ_L(T)
  · inputs/structure/workflow_method.txt ke-dft-cpu S1_opt：DIM=（2D/3D 口径）+ CONTCAR
  三份输入的 source.json 必须指向同一个材料目录、同一个 POSCAR 指纹，否则拒绝汇总——
  不再去远端按目录名找"兄弟技能"（同名材料会串）。

输出（done_marker = zt_summary.json）：
  · zt_summary.json   机读全表：温度网格 × 掺杂网格的 S/σ/κ_e/κ_L/κ_tot/PF/ZT、
                      峰值 ZT、口径说明、来源路径
  · zt_summary.txt    人读汇总（含峰值 ZT 与每档明细表）
  · zt_vs_T.png       各掺杂档 ZT(T)（n 型实线 / p 型虚线）
  · zt_vs_doping.png  各温度下 ZT(载流子浓度)
  · zt_components.png 最优档的 S/σ/κ/PF/ZT 五联图（诊断用）
  · per_temperature/<T>K/report.{txt,json}   每个温度一份全掺杂档明细
  · intrinsic_carrier/  [0052] 本征点（未掺杂，n = p = n_i）随温度的图 + 表：
                      zt_intrinsic.png / zt_intrinsic.csv / zt_intrinsic.txt（含 n_i、E_F、S/σ/κe/κL/PF/ZT）
  [0052] 每次 gen 先删上一轮的产物（温度网格/掺杂档变了不留旧文件）；本地拉回时 autozt 按
  skill.yaml 的 result_flat 直接放进 <材料>/zt-dft-cpu/result/（不再多一层 result/step20_zt/）。

判据：check: plot + done_marker zt_summary.json；画图失败只告警不失败。
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zt_common as zc          # noqa: E402

OUTDIR = "step20_zt"
INPUTS = "inputs"                 # autozt needs_results 推上来的上游结果
INTR_DIR = "intrinsic_carrier"    # [0052] 本征点（未掺杂）随温度的 ZT 图 + 表
AMSET_DIR = INPUTS + "/transport"
KL_DIR = INPUTS + "/kappa_L"
STEP1_CANDS = (INPUTS + "/structure",)
KL_STEP1 = INPUTS + "/kl_step1/workflow_method.txt"   # skill.yaml push_paths 推来的 kl S1 记录（可能没有）
T_TARGETS = (300.0, 500.0, 700.0, 900.0)   # zt_vs_doping 想画的目标温度

# ---- [0047] zT(n,T) 热图的【跨材料统一色板】 ---------------------------------------------
# 同一套颜色 + 同一个数值范围 -> 不同材料的 zt_heatmap.png 颜色含义一致，可直接并排比较。
# 单色相顺序色板「白 -> 浅粉 -> 酒红」（浅 = ZT 小，深 = ZT 大）。
# 两端按用户 2026-10-14 指定：最深 #B8244A、浅 #F6D3DE；更浅端用近白 #FFF7FA（比纯白略粉，
# 让 ZT≈0 的格子仍能看出边界，不和白格线融成一片）。中间为线性过渡。
ZT_CMAP_COLORS = ("#FFF7FA", "#F6D3DE", "#EAB0C0", "#DD8DA3", "#D16A85", "#C44768", "#B8244A")
ZT_CMAP = "autozt_zt"
ZT_INK = "#B8244A"                     # 标题用的同色系深色
ZT_AXIS = "#000000"                    # [0048] 坐标轴（带箭头的横/纵轴 + 刻度数字/轴名）用黑色
ZT_GRID = "#F2F2F2"                    # [0051] 档位格线：深色区读作白线，浅色区仍能看出边界
ZT_CMAP_VMAX_DEFAULT = "0.3"           # 统一色标上限（ZT）；写 auto = 退回每材料自适应的旧行为


def _params(cwd):
    """读 tf 合成并推下来的 step.conf 的 [params]（读不到就用默认值）。"""
    p = Path(cwd) / "step.conf"
    out = {}
    try:
        import stepconf as _sc
        txt = p.read_text(encoding="utf-8-sig")
        out = {k.upper(): v for k, v, _ in _sc.parse(txt, _sc.CONF_NAME).get("params", [])}
    except Exception:
        try:
            txt = p.read_text(encoding="utf-8-sig")
            for line in txt.splitlines():
                line = line.split("#", 1)[0].strip()
                if "=" in line and not line.startswith("["):
                    k, v = line.split("=", 1)
                    out[k.strip().upper()] = v.strip()
        except Exception:
            pass
    return out


def _read_dim(cwd):
    """从 step1 的 workflow_method.txt 读 DIM=；读不到返回 (None, 说明)。"""
    for d in STEP1_CANDS:
        f = Path(cwd) / d / "workflow_method.txt"
        if not f.is_file():
            continue
        try:
            for line in f.read_text(errors="ignore").splitlines():
                s = line.strip()
                if s.upper().startswith("DIM"):
                    v = s.split("=", 1)[-1].strip().strip('"').lower()
                    if v:
                        return v, "读自 %s/workflow_method.txt" % d
        except OSError:
            continue
    return None, "workflow_method.txt 里没有 DIM=（按 3D 处理）"


def _cell_c_axis(cwd):
    """读 step3_uniform / step1_opt 的结构，返回 c 轴长度（Å）；读不到返回 None。
    用途：与 kl 的 kappa_summary.json 的 Lz_ang 做元胞一致性闸门（同 ke 的 step8.3）。"""
    for rel in STEP1_CANDS:
        for name in ("CONTCAR", "POSCAR"):
            p = Path(cwd) / rel / name
            if not p.is_file():
                continue
            try:
                ln = p.read_text(errors="ignore").splitlines()
                s = float(ln[1].split()[0])
                v = [float(x) * s for x in ln[4].split()[:3]]
                return sum(x * x for x in v) ** 0.5
            except (OSError, IndexError, ValueError):
                continue
    return None


def _cell_caliber_check(cwd, kl):
    """κ_L 与电子段是否在【同一个元胞】上（只查 kl 提供了 Lz_ang 时）。
    两条链若用了不同的胞（例如手改过其中一段），κ_e+κ_L 直接相加就是错的。"""
    lz = kl.get("Lz_ang")
    c = _cell_c_axis(cwd)
    if lz is None or c is None:
        return None
    try:
        lz = float(lz)
    except (TypeError, ValueError):
        return None
    if lz <= 0:
        return None
    rel = abs(c - lz) / lz
    return {"cell_c_A": c, "kappa_L_Lz_A": lz, "rel_diff": rel,
            "ok": rel <= 0.02}


def _read_kv(path):
    kv = {}
    try:
        for ln in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
            if "=" in ln:
                k, v = ln.split("=", 1)
                kv[k.strip().upper()] = v.strip()
    except OSError:
        return None
    return kv


def functional_check(cwd):
    """2026-10-14：κ_e（ke）与 κ_L（kl）的几何是不是同一个泛函弛豫出来的。

    ke 的 FUNC 读 inputs/structure/workflow_method.txt；kl 的读 inputs/kappa_L/kl_method.json
    （kl-dft-cpu S6 起写；kl-mlff / fit-fc-thermal 的 κ_L 没有这个文件 -> 记"未知"，不拦）。
    不一致时只告警、写进 notes（不拦截）：两条链晶格常数不同，κ_e + κ_L 不在同一个晶格上。"""
    import json as _json
    ke = None
    for d in STEP1_CANDS:
        kv = _read_kv(Path(cwd) / d / "workflow_method.txt")
        if kv and kv.get("FUNC"):
            ke = {"func": kv["FUNC"].lower(), "reason": kv.get("FUNC_REASON")}
            break
    kl = None
    p = Path(cwd) / KL_DIR / "kl_method.json"
    if p.is_file():
        try:
            d = _json.loads(p.read_text(encoding="utf-8"))
            if isinstance(d, dict) and d.get("func"):
                kl = {"func": str(d["func"]).lower(), "reason": d.get("func_reason")}
        except (OSError, ValueError):
            kl = None
    if kl is None:
        # kl S6 早于 2026-10-14 没写 kl_method.json：κ_L 确实来自 kl-dft-cpu 时，改读 autozt 推来的 kl S1 记录
        #   （skill.yaml push_paths -> inputs/kl_step1/）。MLFF / fit-fc-thermal 来源的 κ_L 不套用 kl S1 的泛函。
        src = (_source(cwd, "kappa_L") or {}).get("from") or ""
        kv = _read_kv(Path(cwd) / KL_STEP1) if str(src).startswith("kl-dft-cpu/") else None
        if kv and kv.get("FUNC"):
            kl = {"func": kv["FUNC"].lower(),
                  "reason": "（kl S1 的 workflow_method.txt）%s" % (kv.get("FUNC_REASON") or "")}
    rec = {"ke_func": ke and ke["func"], "kl_func": kl and kl["func"],
           "ke_reason": ke and ke["reason"], "kl_reason": kl and kl["reason"]}
    if ke is None or kl is None:
        rec["ok"] = None
        rec["note"] = ("泛函一致性未核对：%s" % "；".join(
            x for x in (("ke 的 workflow_method.txt 没有 FUNC" if ke is None else ""),
                        ("κ_L 结果里没有 kl_method.json（MLFF 来源，或 kl S6 早于 2026-10-14）"
                         if kl is None else "")) if x))
    else:
        rec["ok"] = ke["func"] == kl["func"]
        rec["note"] = ("ke/kl 几何泛函一致：%s" % ke["func"] if rec["ok"] else
                       "★ ke 用 %s、kl 用 %s：两条链晶格常数不同，κ_e+κ_L 不在同一晶格上，ZT 不可直接采信。"
                       "处理：两边材料 step.conf 显式写同一个 FUNC，重跑泛函不对的那条链（S1 起）。"
                       % (ke["func"], kl["func"]))
    return rec


def _amset_settings(cwd):
    """读本技能 step8_amset/settings.yaml 的关键口径（散射机制 / 带隙），写进产物。

    Si 这类非极性体系 AMSET 会按物理口径剔除 POP（只剩 ADP/IMP），把这一条记进
    zt_summary.json，事后核对 ZT 数字口径时不用再去翻远端 settings.yaml。
    """
    p = Path(cwd) / AMSET_DIR / "settings.yaml"
    if not p.is_file():
        return None
    try:
        import yaml
        s = yaml.safe_load(p.read_text(encoding="utf-8", errors="ignore")) or {}
    except Exception:                                   # noqa: BLE001
        return None
    out = {}
    for k in ("scattering_type", "bandgap", "interpolation_factor",
              "free_carrier_screening", "doping", "temperatures"):
        if k in s:
            out[k] = s[k]
    return out or None


def _source(cwd, key):
    """inputs/<键>/source.json（autozt 推送时写）；没有返回 None。"""
    p = Path(cwd) / INPUTS / key / "source.json"
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def check_input_sources(cwd, keys=("transport", "kappa_L", "structure")):
    """三份输入必须来自同一个材料目录、同一个 POSCAR 指纹（防同名材料/旧结果串读）。

    返回 {键: source}；不一致直接退出。"""
    srcs = {k: _source(cwd, k) for k in keys}
    miss = [k for k, v in srcs.items() if not v]
    if miss:
        sys.exit("[ERROR] 缺 %s 的 source.json —— inputs/ 不是 autozt 按 needs_results 推上来的，"
                 "拒绝汇总（不再按目录名去找兄弟技能，避免同名材料串读）。" % "、".join(miss))
    dirs = {k: os.path.realpath(v.get("material_dir") or "") for k, v in srcs.items()}
    if len(set(dirs.values())) != 1:
        sys.exit("[ERROR] 输入来自不同的材料目录：%s" % json.dumps(dirs, ensure_ascii=False))
    shas = {k: (v.get("origin") or {}).get("poscar_sha256") for k, v in srcs.items()}
    if len({x for x in shas.values() if x}) > 1:
        sys.exit("[ERROR] 输入按不同的 POSCAR 算的（结构改过、有旧结果）：%s"
                 % json.dumps(shas, ensure_ascii=False))
    return srcs


def _find_transport(cwd):
    p = Path(cwd) / AMSET_DIR / "transport.json"
    return (p, "%s → %s" % ((_source(cwd, "transport") or {}).get("from", "ke-dft-cpu"), AMSET_DIR)) \
        if p.is_file() else (None, "未找到")


def transport_route(src, is_2d):
    """[0014] 电子段来自哪一步 + 2D 口径核对。返回 (步骤标签, 警告或 None)。

    2D 必须用 S8.4_amset2d（二维散射核）；S8_kappa 的 POP/IMP 是三维核，对 2D 只是近似。
    zt-dft-cpu 的 needs_results 对 2D 只认 S8.4，走到 S8 只可能是项目层改了 needs_results。"""
    frm = str((src or {}).get("from") or "")
    step = ("S8.4_amset2d" if "step8.4_amset2d" in frm else
            "S8_kappa" if "step8_amset" in frm else (frm or "?"))
    warn = None
    if is_2d and step != "S8.4_amset2d":
        warn = ("2D 材料的电子段来自 %s（不是 S8.4_amset2d）：POP/IMP 为三维散射核，"
                "ZT 只作对照，不作二维结论" % step)
    return step, warn


def _find_kappa(cwd):
    p = Path(cwd) / KL_DIR / "kappa_summary.json"
    return (p, "%s → %s" % ((_source(cwd, "kappa_L") or {}).get("from", "?"), KL_DIR)) \
        if p.is_file() else (None, "未找到")


def _have_matplotlib():
    try:
        import matplotlib            # noqa: F401
        return True
    except Exception:
        return False


def heatmap_data(grid):
    """zT(n, T) 热图数据：按载流子类型分开，|n| 升序为列、T 为行；ZT 缺值为 None。
    [0050] 用户 2026-10-14：**整行全空**的温度从热图里剔除（图上不留空白横带）。
    判据按实际数据：某个温度在 n/p 两栏都拿不到 ZT 才删 —— 例如 κ_L 只覆盖到 800 K 时
    900 K 整行是空的（ZT 记 null）就不画；**有 900 K 数据的材料照常保留那一行**，
    不写死温度、不按温度值硬筛。两栏共用同一套温度行，行距保持对齐。"""
    temps = list(grid["temperatures"])
    out = {}
    for kind in ("n", "p"):
        rows = [r for r in grid["rows"] if r.get("type") == kind]
        rows.sort(key=lambda r: abs(r["doping_cm-3"]))
        if not rows:
            continue
        out[kind] = {"doping_cm-3": [abs(r["doping_cm-3"]) for r in rows],
                     "temperatures": temps,
                     "zt": [[r["ZT"][j] for r in rows] for j in range(len(temps))]}
    if not out:
        return out
    keep = [j for j in range(len(temps))                      # 全栏都没值的温度 -> 剔除
            if any(any(v is not None for v in d["zt"][j]) for d in out.values())]
    for d in out.values():
        d["temperatures"] = [temps[j] for j in keep]
        d["zt"] = [d["zt"][j] for j in keep]
    return out


def _zt_cmap(plt):
    """[0047] 统一色板：浅粉 -> 酒红的单色相顺序 cmap（所有材料共用）。"""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(ZT_CMAP, list(ZT_CMAP_COLORS))


def _zt_vmax(par, data):
    """[0047] 统一色标上限：默认所有材料同一把尺子（ZT_CMAP_VMAX，出厂 1.0）；
    写成 auto 时退回「每材料取自身最大值」的旧行为。返回 (上限, 说明文字)。"""
    amax = max((v for d in data.values() for row in d["zt"] for v in row if v is not None),
               default=None)
    spec = str((par or {}).get("ZT_CMAP_VMAX", ZT_CMAP_VMAX_DEFAULT) or "").strip().lower()
    if spec and spec != "auto":
        try:
            v = float(spec)
            if v > 0:
                return v, "%.3g" % v
        except ValueError:
            pass
    return amax, "auto"


def _n_tick(x):
    """浓度档位标签（对数档位）：1e17 -> $10^{17}$；0 档（未掺杂点）标 0。"""
    try:
        if abs(float(x)) < 1e-12:
            return "0"
        return "$10^{%d}$" % int(("%.0e" % float(x)).split("e")[1])
    except (TypeError, ValueError, IndexError):
        return str(x)


def make_heatmap(out, grid, meta, plt, par=None):
    """[0047] zT(n, T) 热图：n 型 / p 型两栏，横轴 |n|（对数档位）、纵轴 T。
    色板与色标【跨材料统一】（ZT_CMAP / ZT_CMAP_VMAX）：浅粉->酒红、白格线分隔档位，
    峰值打星并标数值；超过统一上限的点按上限上色（图下注明）。"""
    import numpy as np
    data = heatmap_data(grid)
    if not data:
        return None
    kinds = [k for k in ("n", "p") if k in data]
    vmax, vsrc = _zt_vmax(par, data)
    if not vmax:
        return None
    clipped = any(v > vmax for d in data.values() for row in d["zt"] for v in row if v is not None)
    cmap = _zt_cmap(plt)
    fig, axes = plt.subplots(1, len(kinds), figsize=(4.7 * len(kinds), 4.6), squeeze=False)
    im = None
    for ax, kind in zip(axes[0], kinds):
        d = data[kind]
        nn, nt = len(d["doping_cm-3"]), len(d["temperatures"])
        z = np.array([[np.nan if v is None else min(float(v), vmax) for v in row]
                      for row in d["zt"]], float)
        im = ax.imshow(z, origin="lower", aspect="auto", cmap=cmap, vmin=0.0, vmax=vmax,
                       interpolation="nearest")
        ax.set_xticks(range(nn))
        ax.set_xticklabels([_n_tick(x) for x in d["doping_cm-3"]], fontsize=8.5)
        ax.set_yticks(range(nt))
        ax.set_yticklabels(["%g" % t for t in d["temperatures"]], fontsize=8.5)
        if nn > 6:                      # 档位多才转 45°，否则横排更清爽
            for lb in ax.get_xticklabels():
                lb.set_rotation(45)
                lb.set_ha("right")
                lb.set_rotation_mode("anchor")
        ax.set_xticks(np.arange(-.5, nn, 1), minor=True)        # 白格线：分隔相邻档位
        ax.set_yticks(np.arange(-.5, nt, 1), minor=True)
        ax.grid(which="minor", color=ZT_GRID, linewidth=1.8)      # [0051] 格线用极浅中性色：
        #   深色区读作白线（同参考图），浅色区（ZT≈0，浅端 #FFF7FA）仍能看出格子边界，不会融成一片
        ax.tick_params(which="minor", length=0)
        ax.tick_params(length=2.5, width=0.8, colors=ZT_AXIS)
        for lb in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            lb.set_color(ZT_AXIS)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xlabel(r"$n$  (cm$^{-3}$)", fontsize=10.5, color=ZT_AXIS, labelpad=3)
        ax.set_ylabel(r"$T$  (K)", fontsize=10.5, color=ZT_AXIS, labelpad=3)
        # [0048] 黑色坐标轴：左/下两条带末端箭头的轴（参考图风格）；刻度数字与轴名同黑
        _x0, _x1, _y0, _y1 = -0.5, nn - 0.5, -0.5, nt - 0.5
        ax.set_xlim(_x0 - 0.06, _x1 + 0.62)
        ax.set_ylim(_y0 - 0.06, _y1 + 0.62)
        for _hd, _tl in (((_x1 + 0.48, _y0), (_x0, _y0)), ((_x0, _y1 + 0.48), (_x0, _y0))):
            ax.annotate("", xy=_hd, xytext=_tl, annotation_clip=False,
                        arrowprops=dict(arrowstyle="-|>", color=ZT_AXIS, lw=1.3,
                                        shrinkA=0, shrinkB=0, mutation_scale=13),
                        zorder=6)
        ax.set_title("%s-type  zT($n$, $T$)" % kind, fontsize=13, pad=9,
                     color=ZT_INK, fontweight="bold")
        if np.isfinite(z).any():
            j, i = np.unravel_index(np.nanargmax(z), z.shape)
            ax.plot(i, j, marker="*", ms=13, mfc="white", mec=ZT_INK, mew=1.0, ls="none", zorder=3)
            ax.annotate("%.2f" % z[j, i], (i, j), xytext=(7, 6), textcoords="offset points",
                        color=("white" if z[j, i] > 0.55 * vmax else ZT_INK),
                        fontsize=9, fontweight="bold", zorder=4)
    cbar = fig.colorbar(im, ax=list(axes[0]), shrink=0.85, pad=0.035,
                        extend=("max" if clipped else "neither"))
    cbar.set_label("zT", fontsize=10, color=ZT_INK)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=2, width=0.7, colors=ZT_INK, labelsize=8.5)
    # [0049] 用户 2026-10-14：去掉最上面的总标题 "zT(n, T) — <材料>"，
    #   以及最下面那行 color scale 脚注（统一色标仍写进 zt_summary 的 notes，图面只留两栏面板）。
    fig.savefig(Path(out) / "zt_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {"vmax": vmax, "vsrc": vsrc, "clipped": clipped}


def make_per_temperature(out, grid, meta, peaks):
    """[0049] 每个温度点一个文件夹：report.txt + report.json。
    每个文件夹放该温度下【所有掺杂档】的 ZT 及电子/声子输运参数（S/σ/κe/κL/κtot/PF），
    并标出该温度的峰值 ZT（n/p 各一）。"""
    import json as _json
    temps = grid["temperatures"]
    kl_used = grid.get("kappa_L_used") or [None] * len(temps)
    rows = grid["rows"]
    root = Path(out) / "per_temperature"
    root.mkdir(exist_ok=True)
    amset_note = next((n for n in grid.get("notes", []) if "散射机制" in n), "")
    n_written = 0
    for jT, T in enumerate(temps):
        folder = root / ("%gK" % T)
        folder.mkdir(exist_ok=True)
        klT = kl_used[jT] if jT < len(kl_used) else None
        recs = []
        for r in rows:
            j = jT if jT < len(r["ZT"]) else None
            if j is None:
                continue
            recs.append({
                "doping_cm-3": r["doping_cm-3"],
                "type": r["type"],
                "seebeck_uV/K": r["seebeck_uV/K"][j],
                "sigma_S/m": r["sigma_S/m"][j],
                "kappa_e_W/mK": r["kappa_e_W/mK"][j],
                "kappa_L_W/mK": r["kappa_L_W/mK"][j],
                "kappa_tot_W/mK": r["kappa_tot_W/mK"][j],
                "PF_W/mK2": r["PF_W/mK2"][j],
                "ZT": r["ZT"][j],
            })
        # 该温度的峰值：n / p 各取 ZT 最大（只取 ZT 非 None）
        peak = {"n": None, "p": None}
        for r in recs:
            if r["ZT"] is None or r["type"] not in ("n", "p"):   # [0052] 本征点（i）不算 n/p 峰值
                continue
            k = r["type"]
            if peak[k] is None or r["ZT"] > peak[k]["ZT"]:
                peak[k] = r
        # ---- report.txt ----
        lines = []
        lines.append("# ZT 与输运参数 @ T = %.1f K  —  %s" % (T, meta["material"]))
        lines.append("# κ_L = %s W/m/K" % ("%.4g" % klT if klT is not None else "n/a"))
        if amset_note:
            lines.append("# %s" % amset_note)
        lines.append("#")
        lines.append("# %-14s %-4s %10s %11s %9s %9s %9s %10s %9s" % (
            "doping(cm-3)", "type", "S(uV/K)", "sigma(S/m)", "kappa_e", "kappa_L",
            "kappa_tot", "PF(W/mK2)", "ZT"))
        _hdr = "%-14.3e  %-4s %10s %11s %9s %9s %9s %10s %9s"
        for r in sorted(recs, key=lambda x: abs(x["doping_cm-3"])):
            lines.append(_hdr % (
                r["doping_cm-3"], r["type"],
                "n/a" if r["seebeck_uV/K"] is None else "%.4g" % r["seebeck_uV/K"],
                "n/a" if r["sigma_S/m"] is None else "%.4g" % r["sigma_S/m"],
                "n/a" if r["kappa_e_W/mK"] is None else "%.4g" % r["kappa_e_W/mK"],
                "n/a" if r["kappa_L_W/mK"] is None else "%.4g" % r["kappa_L_W/mK"],
                "n/a" if r["kappa_tot_W/mK"] is None else "%.4g" % r["kappa_tot_W/mK"],
                "n/a" if r["PF_W/mK2"] is None else "%.4g" % r["PF_W/mK2"],
                "n/a" if r["ZT"] is None else "%.4g" % r["ZT"]))
        for _k in ("n", "p"):
            if peak[_k] is not None:
                lines.append("# %s 型峰值 ZT = %.4g @ %s cm^-3" % (
                    _k, peak[_k]["ZT"], "%.3g" % peak[_k]["doping_cm-3"]))
        (folder / "report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        # ---- report.json ----
        (folder / "report.json").write_text(_json.dumps({
            "material": meta["material"], "T_K": T, "kappa_L_W_mK": klT,
            "rows": recs, "peak": peak,
        }, ensure_ascii=False, indent=2), encoding="utf-8")
        n_written += 1
    return n_written


OUT_PRODUCTS = ("zt_summary.json", "zt_summary.txt", "zt_heatmap.png", "zt_vs_T.png",
                "zt_vs_doping.png", "zt_components.png", "per_temperature", INTR_DIR)


def clear_old_products(out):
    """[0052] 先删上一轮的产物（只删本脚本自己写的那些名字；温度网格 / 掺杂档变了以后
    per_temperature/ 里的旧温度文件夹、旧图不再残留冒充当前结果）。返回删掉的项数。"""
    import shutil
    n = 0
    for name in OUT_PRODUCTS:
        f = Path(out) / name
        try:
            if f.is_dir() and not f.is_symlink():
                shutil.rmtree(f)
                n += 1
            elif f.exists() or f.is_symlink():
                f.unlink()
                n += 1
        except OSError as e:
            print("[WARN] 删旧产物 %s 失败：%s" % (f, e))
    return n


def make_intrinsic(out, grid, meta, ni, ni_src):
    """[0052] intrinsic_carrier/：本征点（AMSET doping = 0，未掺杂，E_F 由 n = p 定）随温度的
    n_i、E_F、S、σ、κe、κL、κtot、PF、ZT —— 表（csv / txt）+ 图（zt_intrinsic.png，n_i 逐点标数值）。
    transport.json 没有 0 掺杂点（ke-dft-cpu 的 UNDOPED_POINT 关了）-> 不建目录，返回 None。"""
    import csv
    rows = zc.intrinsic_rows(grid, ni)
    if not rows:
        return None
    d = Path(out) / INTR_DIR
    d.mkdir(exist_ok=True)
    cols = ["T_K", "n_i_cm-3", "E_F_eV", "seebeck_uV/K", "sigma_S/m", "kappa_e_W/mK", "kappa_L_W/mK",
            "kappa_tot_W/mK", "PF_W/mK2", "ZT"]
    head = ["# %s 本征载流子（未掺杂，AMSET doping = 0，E_F 由 n = p = n_i 定）随温度" % meta["material"],
            "# n_i 来源：%s" % (ni_src if ni else "无（%s；重跑 ke-dft-cpu S8/S8.4 后可得）" % ni_src),
            "# n_i 为每胞体积口径 cm^-3%s；ZT = S^2 sigma T / (kappa_e + kappa_L)"
            % ("（2D = 含真空元胞）" if meta.get("dim") == "2d" else "")]
    with open(d / "zt_intrinsic.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow(["" if r[c] is None else ("%.6g" % r[c]) for c in cols])
    L = head + ["#", "  %-7s %-11s %-9s %-10s %-11s %-10s %-10s %-10s %-11s %-8s" % (
        "T(K)", "n_i(cm-3)", "E_F(eV)", "S(uV/K)", "sigma(S/m)", "k_e", "k_L", "k_tot", "PF(W/mK2)", "ZT")]
    for r in rows:
        L.append("  %-7s %-11s %-9s %-10s %-11s %-10s %-10s %-10s %-11s %-8s" % tuple(
            zc._fmt(r[c], 0 if c == "T_K" else 4) for c in cols))
    (d / "zt_intrinsic.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    if not _have_matplotlib():
        return d
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [("seebeck_uV/K", "S (uV/K)", False), ("sigma_S/m", "sigma (S/m)", True),
              (("kappa_e_W/mK", "kappa_L_W/mK"), "kappa (W/m/K)", True),
              ("PF_W/mK2", "PF (W/m/K$^2$)", True), ("ZT", "ZT", False)]
    nip = [(r["T_K"], r["n_i_cm-3"]) for r in rows if r["n_i_cm-3"]]
    if nip:
        panels.insert(0, ("n_i_cm-3", "n_i = n = p (cm$^{-3}$)", True))
    fig, axs = plt.subplots(1, len(panels), figsize=(3.9 * len(panels), 3.6), squeeze=False)
    for a, (key, ttl, logy) in zip(axs[0], panels):
        for k, lab in ((key, None),) if isinstance(key, str) else zip(key, ("kappa_e", "kappa_L")):
            pts = [(r["T_K"], r[k]) for r in rows if r[k] is not None and (not logy or r[k] > 0)]
            if pts:
                a.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=3.5, lw=1.3,
                       color="k" if k == "n_i_cm-3" else None, label=lab)
            if k == "n_i_cm-3":
                for x, y in pts:
                    a.annotate("%.2g" % y, (x, y), textcoords="offset points", xytext=(0, 5),
                               ha="center", fontsize=6.5)
                a.margins(y=0.25)
        if logy:
            a.set_yscale("log")
        if not isinstance(key, str):
            a.legend(fontsize=7)
        a.set_title(ttl, fontsize=10)
        a.set_xlabel("T (K)", fontsize=8)
        a.grid(alpha=0.3)
    fig.suptitle("%s  intrinsic carriers (undoped, n = p = n_i)  —  ZT vs T%s"
                 % (meta["material"], "" if nip else "   [n_i not recorded: rerun ke-dft-cpu S8]"), fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(d / "zt_intrinsic.png", dpi=130)
    plt.close(fig)
    return d


def make_plots(out, grid, peaks, meta, par=None):
    """三张图；任何一张失败只告警，不影响 JSON/TXT 产物。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = grid["temperatures"]
    make_heatmap(out, grid, meta, plt, par=par)
    # ---- 1) ZT(T) ----
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for i, r in enumerate(grid["rows"]):
        ys = r["ZT"]
        if all(v is None for v in ys):
            continue
        xs = [t for t, v in zip(temps, ys) if v is not None]
        yy = [v for v in ys if v is not None]
        ls = {"n": "-", "p": "--"}.get(r["type"], ":")
        ax.plot(xs, yy, ls, marker="o", ms=3.5, lw=1.3,
                label=("intrinsic (undoped)" if r["type"] == "i" else
                       "%s %s cm$^{-3}$" % (r["type"], zc._fmt(abs(r["doping_cm-3"]), 1))))
    ax.set_xlabel("T (K)")
    ax.set_ylabel("ZT")
    ax.set_title("ZT(T) — %s   (S$^2\\sigma$T/($\\kappa_e$+$\\kappa_L$))" % meta["material"])
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(Path(out) / "zt_vs_T.png", dpi=140)
    plt.close(fig)

    # ---- 2) ZT(n) ----
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for T in T_TARGETS:
        xs, ys, marks = [], [], []
        for i, r in enumerate(grid["rows"]):
            if T not in temps:
                jT = min(range(len(temps)), key=lambda k: abs(temps[k] - T))
                if abs(temps[jT] - T) > 1e-6:
                    continue
            else:
                jT = temps.index(T)
            v = r["ZT"][jT]
            if v is None:
                continue
            xs.append(r["doping_cm-3"])
            ys.append(v)
        if not xs:
            continue
        order = sorted(range(len(xs)), key=lambda k: xs[k])
        ax.plot([xs[k] for k in order], [ys[k] for k in order], "o-", ms=3.5,
                lw=1.3, label="%g K" % T)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xscale("symlog", linthresh=1e16)
    ax.set_xlabel("carrier concentration (cm$^{-3}$, negative = n-type)")
    ax.set_ylabel("ZT")
    ax.set_title("ZT vs doping — %s" % meta["material"])
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(Path(out) / "zt_vs_doping.png", dpi=140)
    plt.close(fig)

    # ---- 3) 最优档五联图 ----
    best = peaks.get("p") or peaks.get("n")
    if best:
        i = grid["doping"].index(best["doping_cm-3"])
        r = grid["rows"][i]
        fig, axes = plt.subplots(1, 5, figsize=(19, 3.6))
        for ax, key, lab in zip(axes,
                                ("seebeck_uV/K", "sigma_S/m", "kappa_tot_W/mK",
                                 "PF_W/mK2", "ZT"),
                                ("S (uV/K)", "sigma (S/m)", "kappa_tot (W/m/K)",
                                 "PF (W/m/K^2)", "ZT")):
            ax.plot(temps, r[key], "o-", ms=3.5, lw=1.3)
            ax.set_xlabel("T (K)")
            ax.set_ylabel(lab)
            ax.grid(alpha=0.3)
        fig.suptitle("best %s-type doping %s cm^-3" %
                     (r["type"], zc._fmt(r["doping_cm-3"], 2)))
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fig.savefig(Path(out) / "zt_components.png", dpi=130)
        plt.close(fig)


def _material_name(cwd, argv):
    """材料名：优先 gen 命令行带进来的 --material（tf 的 {mat} 占位符展开），
    否则用技能目录的上一级目录名。"""
    for i, a in enumerate(argv):
        if a == "--material" and i + 1 < len(argv):
            return argv[i + 1]
        if a.startswith("--material="):
            return a.split("=", 1)[1]
    return cwd.parent.name or cwd.name


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    _nold = clear_old_products(out)               # [0052] 先删旧结果
    if _nold:
        print("[..] 已删除上一轮的 %d 项旧产物" % _nold)

    srcs = check_input_sources(cwd)
    tr_path, tr_src = _find_transport(cwd)
    if tr_path is None:
        sys.exit("[ERROR] 找不到 %s/transport.json —— 先把 ke-dft-cpu 的 S8_kappa（3D）"
                 "或 S8.4_amset2d（2D）跑完并拉回。" % AMSET_DIR)
    kl_path, kl_src = _find_kappa(cwd)
    if kl_path is None:
        sys.exit("[ERROR] 找不到 %s/kappa_summary.json —— 先把 kl-dft-cpu 的 S6_kappa 跑完并拉回。"
                 % KL_DIR)

    par = _params(cwd)
    ktemp_mode = str(par.get("KTEMP_MODE", "interp") or "interp").lower()
    const_T = float(par.get("KL_CONST_T", "300") or 300)

    dim, dim_note = _read_dim(cwd)
    is_2d = (dim == "2d")
    print("[..] %s" % dim_note)
    tr_step, tr_warn = transport_route(srcs.get("transport"), is_2d)
    if tr_warn:
        print("[WARN] " + tr_warn)

    try:
        tr = zc.load_transport(str(tr_path))
    except Exception as e:
        sys.exit("[ERROR] 解析 AMSET transport.json 失败：%s: %s"
                 % (type(e).__name__, e))
    try:
        kl = zc.load_kappa(str(kl_path))
    except Exception as e:
        sys.exit("[ERROR] 解析 %s 失败：%s: %s" % (kl_src, type(e).__name__, e))
    print("[..] 电子段 %s：%d 个掺杂档 × %d 个温度点"
          % (tr_src, len(tr["doping"]), len(tr["temperatures"])))
    if kl.get("single_T_fallback"):
        print("[WARN] κ_L 来源只有 300 K 单值（无惰温数组）——interp 模式下只有 300 K 能算 ZT；"
              "要全温曲线请把 step20_zt 的 KTEMP_MODE 设为 const300。")
    print("[..] 晶格段 %s：%d 个温度点，κ_L(300K 近似)=%.3f W/m/K"
          % (kl_src, len(kl["temperatures"]),
             zc.interp_T(kl["temperatures"], zc.kappa_scalar(kl, is_2d), 300.0)
             or float("nan")))

    cal = _cell_caliber_check(cwd, kl)
    if cal is not None:
        if cal["ok"]:
            print("[..] 元胞口径闸门：c=%.3f Å vs κ_L 的 Lz=%.3f Å（差 %.1f%%）——通过"
                  % (cal["cell_c_A"], cal["kappa_L_Lz_A"], cal["rel_diff"] * 100))
        else:
            print("[WARN] ★ 元胞口径闸门不通过：c=%.3f Å vs κ_L 的 Lz=%.3f Å（差 %.1f%%）"
                  % (cal["cell_c_A"], cal["kappa_L_Lz_A"], cal["rel_diff"] * 100))
            print("       κ_L 与 σ/κ_e 可能不在同一个胞上，κ_e+κ_L 的 ZT 不可直接采信。")

    fcheck = functional_check(cwd)
    print("[%s] %s" % ("..", fcheck["note"]) if fcheck["ok"] is not False
          else "[WARN] %s" % fcheck["note"])

    grid = zc.build_grid(tr, kl, is_2d=is_2d, ktemp_mode=ktemp_mode,
                         const_T=const_T)
    grid["notes"].append(fcheck["note"])
    _amset_s = _amset_settings(cwd)
    if _amset_s and _amset_s.get("scattering_type"):
        grid["notes"].append("电子段 AMSET 散射机制 = %s；带隙 = %s eV"
                             % (", ".join(map(str, _amset_s["scattering_type"])),
                                _amset_s.get("bandgap", "?")))
    grid["notes"].append("电子段来自 ke-dft-cpu %s" % tr_step)
    _hm = heatmap_data(grid)                      # [0047] 统一色板/色标也写进人读汇总
    if _hm:
        _v, _vsrc = _zt_vmax(par, _hm)
        grid["notes"].append("zT(n,T) 热图：跨材料统一色板（浅粉→酒红）；色标 0–%s%s"
                             % (_vsrc, "（各材料自身上限，未统一、不可横向比色）" if _vsrc == "auto"
                                else "（所有材料同一把尺子，颜色可直接横向比较）"))
    if tr_warn:
        grid["notes"].append("★ " + tr_warn)
    if cal is not None:      # 闸门结论也进人读汇总（notes 同时进 JSON 与 TXT）
        grid["notes"].append("元胞口径闸门：电子段胞 c=%.3f Å vs κ_L 的 Lz=%.3f Å"
                             "（差 %.1f%%）——%s"
                             % (cal["cell_c_A"], cal["kappa_L_Lz_A"],
                                cal["rel_diff"] * 100,
                                "通过" if cal["ok"] else "★ 不通过，ZT 不可直接采信"))
    peaks = zc.peak_zt(grid)

    # 健壮性闸门：一个可算的 ZT 都没有 = κ_L 温区与电子段温度网格完全不重叠
    # （或 κ_L 全为非有限值）。照常写出 zt_summary.json 会让步骤被判"完成"，
    # 下游拿到的却是空表 —— 直接报错，让人看见。
    if not any(v is not None for row in grid["grid_zt"] for v in row):
        _kt = grid["kappa_L_temperatures"] or []
        _et = grid["temperatures"] or []
        sys.exit("[ERROR] 没有任何温度点能算出 ZT：κ_L 温区 %s K 与电子段温度网格 %s K"
                 " 不重叠（或 κ_L 非有限）。\n"
                 "        可用 KTEMP_MODE=const300 固定取 κ_L(300K)，或检查 kl 的 "
                 "kappa_summary.json。"
                 % ("~%.0f" % _kt[0] + ("" if len(_kt) < 2 else "~%.0f" % _kt[-1]) if _kt else "?",
                    "~%.0f" % _et[0] + ("" if len(_et) < 2 else "~%.0f" % _et[-1]) if _et else "?"))

    mat = _material_name(cwd, sys.argv[1:])
    meta = {"material": mat, "skill": "zt-dft-cpu", "step": OUTDIR,
            "transport_src": tr_src,
            "kappa_src": (str(kl_path.relative_to(cwd.parent))
                          if cwd.parent in kl_path.parents else str(kl_path)),
            "dim": dim or "3d(默认)", "dim_note": dim_note,
            "ktemp_mode": ktemp_mode, "generated_by": "gen_step20_zt.py"}

    payload = {
        "ZT_DONE": True,
        "material": mat,
        "formula": "ZT = S^2 * sigma * T / (kappa_e + kappa_L)",
        "units": {"S": "uV/K", "sigma": "S/m", "kappa_e": "W/m/K",
                  "kappa_L": "W/m/K", "T": "K", "PF": "W/m/K^2"},
        "dim": dim or "3d",
        "temperatures": grid["temperatures"],
        "doping_cm-3": grid["doping"],
        "doping_type": [zc.carrier_type(d) for d in grid["doping"]],
        "kappa_L": {"source": kl_src, "temperatures": grid["kappa_L_temperatures"],
                    "scalar_W_mK": grid["kappa_L_scalar"],
                    "used_W_mK": grid["kappa_L_used"],
                    "xx": kl["k_xx"], "yy": kl["k_yy"], "zz": kl["k_zz"],
                    "caliber": "cell (kappa_xx_yy_zz)",
                    "Lz_ang": kl.get("Lz_ang"),
                    "thickness_d_ang": kl.get("thickness_d_ang"),
                    "thickness_convention": kl.get("thickness_convention"),
                    "kappa_300K_xx_yy_zz": kl.get("kappa_300K_xx_yy_zz")},
        "amset_settings": _amset_settings(cwd),
        "transport": {"source": tr_src, "step": tr_step, "warning": tr_warn},
        "kappa_e": {"source": meta["transport_src"],
                    "note": "AMSET electronic_thermal_conductivity（含双极项；未单独拆分）"},
        "grid_ZT": grid["grid_zt"],
        "rows": grid["rows"],
        "peak_ZT": peaks,
        "notes": grid["notes"] + [
            "κ_L 与 σ/κ_e 均为【元胞口径】，可相加；不用 kl 的 kappa_2d_normalized_*。",
            "ZT 只在 κ_L 温区内计算（不外推）；区域外的温度点 ZT 记 null。",
        ],
        "cell_caliber_check": cal,
        "functional_check": fcheck,
        "sources": {"transport_json": meta["transport_src"], "kappa_src": kl_src,
                    # 来源戳：哪个材料目录、哪个结构、哪个技能的哪一步（防同名串读）
                    "inputs": {k: {"from": v.get("from"),
                                   "material_dir": v.get("material_dir"),
                                   "origin": v.get("origin")} for k, v in srcs.items()}},
    }
    (out / "zt_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "zt_summary.txt").write_text(zc.text_report(grid, peaks, meta),
                                        encoding="utf-8")
    print("[OK] zt_summary.json / zt_summary.txt 已生成")
    for kind in ("n", "p"):
        p = peaks.get(kind)
        if p:
            print("[ZT] %s 型峰值 ZT=%.3f @ %.0f K, %.3g cm^-3"
                  % (kind, p["ZT"], p["T_K"], p["doping_cm-3"]))
        else:
            print("[ZT] %s 型无可用点（κ_L 温区没覆盖电子段温度？）" % kind)

    # [0049] 每个温度点一个文件夹（report.txt + report.json；不依赖 matplotlib）
    try:
        _npt = make_per_temperature(out, grid, meta, peaks)
        print("[OK] per_temperature/ 下已生成 %d 个温度文件夹（report.txt / report.json）" % _npt)
    except Exception as e:                            # noqa: BLE001
        print("[WARN] 每温度文件夹生成失败（%s: %s）——ZT 汇总仍正常" % (type(e).__name__, e))

    # [0052] intrinsic_carrier/：本征点随温度的图 + 表（n_i 读 S8 运行目录的 undoped_ni.json / 日志）
    try:
        _ni, _ni_src = zc.load_undoped_ni(str(tr_path.parent))
        _d = make_intrinsic(out, grid, meta, _ni, _ni_src)
        if _d is None:
            print("[..] transport.json 没有 0 掺杂点（ke-dft-cpu UNDOPED_POINT 关了）—— 不出 %s/" % INTR_DIR)
        else:
            print("[OK] %s/ 已生成（zt_intrinsic.png / .csv / .txt；n_i %s）"
                  % (INTR_DIR, "读自 " + _ni_src if _ni else "缺：" + _ni_src))
    except Exception as e:                            # noqa: BLE001
        print("[WARN] 本征点图表生成失败（%s: %s）——ZT 汇总仍正常" % (type(e).__name__, e))

    if _have_matplotlib():
        try:
            make_plots(out, grid, peaks, meta, par=par)
            print("[OK] zt_heatmap.png / zt_vs_T.png / zt_vs_doping.png / zt_components.png 已生成")
        except Exception as e:                          # noqa: BLE001
            print("[WARN] 画图失败（%s: %s）——JSON/TXT 已生成，图跳过"
                  % (type(e).__name__, e))
    else:
        print("[WARN] 没有 matplotlib —— 跳过画图（JSON/TXT 仍然是完整产物）")

    print("[DONE] %s：zt_summary.json / zt_summary.txt / zt_vs_T.png 已生成" % OUTDIR)


if __name__ == "__main__":
    main()
