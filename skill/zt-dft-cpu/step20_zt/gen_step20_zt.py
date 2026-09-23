#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step20_zt.py —— ZT 全流程汇总（zt-dft-cpu 的收尾步，run: gen）。

输入（相对技能目录，cwd = <超算 work_dir>/<材料>/zt-dft-cpu/）：
  · step8_amset/transport.json          AMSET：S(μV/K) / σ(S/m) / κ_e(W/m/K)，[掺杂][温度][3][3]
  · step6_kappa/kappa_summary.json      phono3py BTE：κ_L(T)，W/m/K（元胞口径）
      找不到就退而求其次读【兄弟技能】的材料目录（../kl-dft-cpu/step6_kappa 等）——
      这样"只补跑 ZT 汇总"或"κ_L 用 MACE 链路算的"也能出 ZT；用它时产物里
      会记 sources.kappa_src，看得见来源。
  · step1_opt/workflow_method.txt       DIM=（决定张量约化口径 2D/3D）

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
AMSET_DIR = "step8_amset"
KL_DIR = "step6_kappa"
# 兄弟技能兜底（同一材料目录下、同一 work_dir）：(技能子目录, 步骤目录)
KL_SIBLINGS = (("kl-dft-cpu", "step6_kappa"),
               ("kl-mlff-cpu", "step4_kappa"),
               ("kl-mlff-gpu", "step4_kappa"))
# 电子段同理：若本技能的 S8_kappa_e 没跑（例如电子输运已由独立的 ke-dft-cpu
# 项目算过），允许读同级 ke-dft-cpu 目录的 transport.json —— 来源会写进产物。
AMSET_SIBLINGS = (("ke-dft-cpu", "step8_amset"),
                  ("ke-dft-cpu", "step8.4_amset2d"))
STEP1_CANDS = ("step1_opt", "step1_std_opt")
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
    for rel in ("step3_uniform", "step1_opt", "step1_std_opt"):
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


def _find_transport(cwd):
    """找 AMSET transport.json。返回 (path, src 说明)。本技能优先，其次兄弟技能。"""
    p = Path(cwd) / AMSET_DIR / "transport.json"
    if p.is_file():
        return p, "本技能 %s/transport.json" % AMSET_DIR
    base = Path(cwd).parent
    for sub, step in AMSET_SIBLINGS:
        q = base / sub / step / "transport.json"
        if q.is_file():
            return q, "兄弟技能 %s/%s/transport.json" % (sub, step)
    return None, "未找到"


def _find_kappa(cwd):
    """找 κ_L 汇总。返回 (path, src 说明)。同技能优先，其次兄弟技能目录。"""
    p = Path(cwd) / KL_DIR / "kappa_summary.json"
    if p.is_file():
        return p, "本技能 %s/kappa_summary.json" % KL_DIR
    base = Path(cwd).parent
    for sub, step in KL_SIBLINGS:
        q = base / sub / step / "kappa_summary.json"
        if q.is_file():
            return q, "兄弟技能 %s/%s/kappa_summary.json" % (sub, step)
    return None, "未找到"


def _have_matplotlib():
    try:
        import matplotlib            # noqa: F401
        return True
    except Exception:
        return False


def make_plots(out, grid, peaks, meta):
    """三张图；任何一张失败只告警，不影响 JSON/TXT 产物。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = grid["temperatures"]
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

    tr_path, tr_src = _find_transport(cwd)
    if tr_path is None:
        sys.exit("[ERROR] 找不到电子段产物：本技能 %s/transport.json 与兄弟技能 "
                 "%s 都没有。先把电子段（S8_kappa_e）跑完再汇总 ZT。"
                 % (AMSET_DIR, "、".join("%s/%s" % x for x in AMSET_SIBLINGS)))
    kl_path, kl_src = _find_kappa(cwd)
    if kl_path is None:
        sys.exit("[ERROR] 找不到 κ_L：本技能 %s/kappa_summary.json 与兄弟技能 "
                 "%s 都没有。先把晶格段（SK6_kappa）跑完。"
                 % (KL_DIR, "、".join("%s/%s" % x for x in KL_SIBLINGS)))

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
    if not str(tr_path).startswith(str(cwd)):
        print("[WARN] 电子段用的是兄弟技能目录：%s" % tr_src)
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
    _amset_s = _amset_settings(cwd) if tr_src.startswith("本技能") else None
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
        "transport": {"source": tr_src,
                      "sibling_fallback": not str(tr_path).startswith(str(cwd))},
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
        "sources": {"transport_json": meta["transport_src"], "kappa_src": kl_src},
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
            print("[OK] zt_vs_T.png / zt_vs_doping.png / zt_components.png 已生成")
        except Exception as e:                          # noqa: BLE001
            print("[WARN] 画图失败（%s: %s）——JSON/TXT 已生成，图跳过"
                  % (type(e).__name__, e))
    else:
        print("[WARN] 没有 matplotlib —— 跳过画图（JSON/TXT 仍然是完整产物）")

    print("[DONE] %s：zt_summary.json / zt_summary.txt / zt_vs_T.png 已生成" % OUTDIR)


if __name__ == "__main__":
    main()
