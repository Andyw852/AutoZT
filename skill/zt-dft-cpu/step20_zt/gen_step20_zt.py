#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step20_zt.py —— ZT 全流程汇总（zt-dft-cpu 的收尾步，run: gen）。

输入（cwd = <超算 work_dir>/<材料>/zt-dft-cpu/；autozt 按 skill.yaml 的 needs_results
从【本材料的本地目录】选出上游技能已拉回的结果，推到 inputs/<键>/，附 source.json）：
  · inputs/transport/transport.json      ke-dft-cpu S8_kappa（或 S8.4_amset2d）：S/σ/κ_e
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
AMSET_DIR = INPUTS + "/transport"
KL_DIR = INPUTS + "/kappa_L"
STEP1_CANDS = (INPUTS + "/structure",)
T_TARGETS = (300.0, 500.0, 700.0, 900.0)   # zt_vs_doping 想画的目标温度


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
    return (p, "ke-dft-cpu → %s" % AMSET_DIR) if p.is_file() else (None, "未找到")


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
    """zT(n, T) 热图数据：按载流子类型分开，|n| 升序为列、T 为行；ZT 缺值为 None。"""
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
    return out


def make_heatmap(out, grid, meta, plt):
    """zT(n, T) 热图：n 型 / p 型两栏，横轴 |n|（对数刻度的档位）、纵轴 T。"""
    import numpy as np
    data = heatmap_data(grid)
    if not data:
        return
    kinds = [k for k in ("n", "p") if k in data]
    vmax = max((v for d in data.values() for row in d["zt"] for v in row
                if v is not None), default=None)
    if not vmax:
        return
    fig, axes = plt.subplots(1, len(kinds), figsize=(5.2 * len(kinds) + 1.0, 4.4),
                             squeeze=False)
    im = None
    for ax, kind in zip(axes[0], kinds):
        d = data[kind]
        z = np.array([[np.nan if v is None else v for v in row] for row in d["zt"]], float)
        im = ax.imshow(z, origin="lower", aspect="auto", cmap="magma", vmin=0, vmax=vmax,
                       interpolation="nearest")
        ax.set_xticks(range(len(d["doping_cm-3"])))
        ax.set_xticklabels(["%.0e" % x for x in d["doping_cm-3"]], rotation=45,
                           ha="right", fontsize=7)
        ax.set_yticks(range(len(d["temperatures"])))
        ax.set_yticklabels(["%g" % t for t in d["temperatures"]], fontsize=7)
        ax.set_xlabel("%s-type carrier concentration (cm$^{-3}$)" % kind)
        ax.set_ylabel("T (K)")
        ax.set_title("%s-type zT(n, T)" % kind)
        j, i = np.unravel_index(np.nanargmax(z), z.shape) if np.isfinite(z).any() else (None, None)
        if j is not None:
            ax.plot(i, j, marker="*", ms=12, mfc="white", mec="black", mew=0.8)
            ax.annotate("%.2f" % z[j, i], (i, j), xytext=(6, 6), textcoords="offset points",
                        color="black" if z[j, i] > 0.6 * vmax else "white",
                        fontsize=8, fontweight="bold")
    fig.colorbar(im, ax=list(axes[0]), label="zT", shrink=0.9)
    fig.suptitle("zT(n, T) — %s" % meta["material"])
    fig.savefig(Path(out) / "zt_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def make_plots(out, grid, peaks, meta):
    """三张图；任何一张失败只告警，不影响 JSON/TXT 产物。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = grid["temperatures"]
    make_heatmap(out, grid, meta, plt)
    # ---- 1) ZT(T) ----
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for i, r in enumerate(grid["rows"]):
        ys = r["ZT"]
        if all(v is None for v in ys):
            continue
        xs = [t for t, v in zip(temps, ys) if v is not None]
        yy = [v for v in ys if v is not None]
        ls = "-" if r["type"] == "n" else "--"
        ax.plot(xs, yy, ls, marker="o", ms=3.5, lw=1.3,
                label="%s %s cm$^{-3}$" % (r["type"], zc._fmt(abs(r["doping_cm-3"]), 1)))
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

    srcs = check_input_sources(cwd)
    tr_path, tr_src = _find_transport(cwd)
    if tr_path is None:
        sys.exit("[ERROR] 找不到 %s/transport.json —— 先把 ke-dft-cpu 的 S8_kappa 跑完并拉回。"
                 % AMSET_DIR)
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

    grid = zc.build_grid(tr, kl, is_2d=is_2d, ktemp_mode=ktemp_mode,
                         const_T=const_T)
    _amset_s = _amset_settings(cwd)
    if _amset_s and _amset_s.get("scattering_type"):
        grid["notes"].append("电子段 AMSET 散射机制 = %s；带隙 = %s eV"
                             % (", ".join(map(str, _amset_s["scattering_type"])),
                                _amset_s.get("bandgap", "?")))
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
        "transport": {"source": tr_src},
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

    if _have_matplotlib():
        try:
            make_plots(out, grid, peaks, meta)
            print("[OK] zt_heatmap.png / zt_vs_T.png / zt_vs_doping.png / zt_components.png 已生成")
        except Exception as e:                          # noqa: BLE001
            print("[WARN] 画图失败（%s: %s）——JSON/TXT 已生成，图跳过"
                  % (type(e).__name__, e))
    else:
        print("[WARN] 没有 matplotlib —— 跳过画图（JSON/TXT 仍然是完整产物）")

    print("[DONE] %s：zt_summary.json / zt_summary.txt / zt_vs_T.png 已生成" % OUTDIR)


if __name__ == "__main__":
    main()
