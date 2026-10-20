#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S6.1_plot: kl-dft-cpu 全流程出图归集（run: gen，登录节点，不交作业）。

★ 设计原则（2026-10-18 用户定）：
  1) **不改动**既有各步自动出的图（S5.1_plot 的声子谱 band_*.png 保持原样、仍在 step5_phonon_plot/）；
  2) 本步只做**归集**：把 S5.1 声子谱 + S6 κL 图/表，统一复制到
     `<材料>/kl-dft-cpu/phonon-transport_plot/`，作为交付给用户/论文的成品目录；
  3) 幂等：重复跑只覆盖同名文件，不删其它内容。

产物目录 `<材料>/kl-dft-cpu/phonon-transport_plot/`（与 stepN_* 同级）:
    phonon_band/            声子色散（S5.1 原图直拷，不改不重画）
    kappa/                  κL 图（ke 风格多面板）+ 表
    README.md               说明（自动生成）

完成标记: phonon-transport_plot/plot_manifest.json
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

OUT_NAME = "phonon-transport_plot"


def find_kappa_summary(root: Path) -> Path:
    cands = [root / "step6_kappa" / "kappa_summary.json"]
    cands += sorted(root.glob("*/step6_kappa/kappa_summary.json"))
    for c in cands:
        if c.is_file():
            return c
    sys.exit("[ERROR] 找不到 step6_kappa/kappa_summary.json（S6 未完成？）")


def find_phonon_plot_dir(root: Path):
    for c in [root / "step5_phonon_plot"] + sorted(root.glob("*/step5_phonon_plot")):
        if c.is_dir() and any(c.glob("band_*.png")):
            return c
    return None


def main():
    # autozt 的 run: gen 步骤在【技能目录】内执行：cwd = <材料>/kl-dft-cpu/
    # 产物必须写 cwd/phonon-transport_plot/；材料名取父目录。兼容手动从材料根/result 调起。
    cwd = Path.cwd()
    if (cwd / "step6_kappa" / "kappa_summary.json").is_file():
        skill_dir = cwd
    elif (cwd / "kl-dft-cpu" / "step6_kappa" / "kappa_summary.json").is_file():
        skill_dir = cwd / "kl-dft-cpu"
    else:
        skill_dir = cwd
    mat = skill_dir.parent.name if skill_dir.name in ("kl-dft-cpu", "result") else skill_dir.name

    out_root = skill_dir / OUT_NAME
    phonon_out = out_root / "phonon_band"
    kappa_out = out_root / "kappa"
    for d in (phonon_out, kappa_out):
        d.mkdir(parents=True, exist_ok=True)

    manifest = {"material": mat, "generator": "gen_step6.1_plot.py",
                "phonon_band": [], "kappa": [],
                "note": "S5.1 原图直拷、未重画；κL 图为 ke 风格多面板"}

    # 1) 声子谱：原样复制（不改、不重画）
    src = find_phonon_plot_dir(skill_dir)
    if src:
        n = 0
        for f in sorted(src.iterdir()):
            if f.is_file() and f.suffix.lower() in (".png", ".dat", ".yaml", ".path", ".json"):
                shutil.copy2(f, phonon_out / f.name)
                manifest["phonon_band"].append(f.name)
                n += 1
        print(f"[copy] 声子谱 {n} 个文件 <- {src}（原图直拷，未重画）")
    else:
        print("[warn] 未找到 step5_phonon_plot 声子谱产物")

    # 2) κL 图/表（ke 风格）
    sp = find_kappa_summary(skill_dir)
    d = json.loads(sp.read_text(encoding="utf-8"))
    if str(d.get("KAPPA_DONE", "")).lower() not in ("true", "1"):
        sys.exit(f"[ERROR] kappa_summary 未标记 KAPPA_DONE: {sp}")

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import kappa_plots as kp

    bte = d.get("bte_method") or "rta"
    mc = (d.get("mesh_convergence") or {}).get(bte, {})
    cs = d.get("cutoff_selection") or {}
    print("[converge] q-mesh     : %s (max rel %.3f%%)" % (
        "converged" if mc.get("converged_5pct") else "NOT converged",
        float(mc.get("max_rel_change_pct") or 0)))
    print("[converge] FC3 cutoff : %s" % (cs.get("status") or "unknown"))

    made = []
    p = kp.plot_overview(d, mat, kappa_out)
    if p:
        made.append(p); print("[fig]", p)
    p = kp.plot_anisotropy(d, mat, kappa_out)
    if p:
        made.append(p); print("[fig]", p)
    made += kp.write_tables(d, mat, kappa_out)
    for x in made:
        if str(x).endswith((".csv", ".json")):
            print("[table]", x)
    manifest["kappa"] = [Path(x).name for x in made]

    # 3) README
    q_ok = mc.get("converged_5pct")
    fc_ok = (cs.get("status") == "converged")
    readme = [
        f"# {mat} —— phonon + transport plots",
        "",
        "本目录由 `kl-dft-cpu` 的 **S6.1_plot** 步骤自动生成（`run: gen`，不交作业）。",
        "",
        "- `phonon_band/`：S5.1_plot 已算好的声子色散图**原样复制**（不改动、不重画）",
        "- `kappa/`：S6_kappa 结果的 κL 图与表（ke 风格多面板）",
        "",
        "## 收敛诊断（两类收敛分开看）",
        "",
        "| 项 | 状态 | 说明 |",
        "|---|---|---|",
        f"| q-mesh 收敛 | {'converged' if q_ok else '**NOT converged**'} | "
        f"最细 {d.get('mesh')}，最大相对变化 {float(mc.get('max_rel_change_pct') or 0):.3f}% |",
        f"| FC3 截断收敛 | {'converged' if fc_ok else '**not converged**'} | "
        f"采用 {d.get('chosen_cutoff_A')} Å；{cs.get('reason') or ''} |",
        "",
        "> 两者是**不同**的收敛判据：q-mesh 决定 BTE 积分精度；FC3 截断决定三阶力常数的空间范围。",
        "",
        "## 文件",
        "",
    ]
    for f in manifest["phonon_band"]:
        readme.append(f"- `phonon_band/{f}`")
    for f in manifest["kappa"]:
        readme.append(f"- `kappa/{f}`")
    readme += ["", "## 说明", "",
               "- 声子谱图来自 S5.1_plot，保持技能原有出图不变；",
               "- κL 原始值依赖 Lz；2D 归一化 κ·h/d 与 Lz 无关，跨算例可比；",
               "- 2D 的 κzz 无物理意义（真空方向），仅作参考。", ""]
    (out_root / "README.md").write_text("\n".join(readme), encoding="utf-8")

    manifest["done"] = True
    mjson = json.dumps(manifest, ensure_ascii=False, indent=2)
    (out_root / "plot_manifest.json").write_text(mjson, encoding="utf-8")
    # 同步写一份到步骤目录 step6.1_plot/，让 autozt 的 done_marker 判读
    #（do_run_gen_step 找 <步骤目录>/plot_manifest.json）能通过。
    #（2026-10-20：产物在 phonon-transport_plot/，但判读路径是 step6.1_plot/，需同步。）
    step_dir = skill_dir / "step6.1_plot"
    step_dir.mkdir(parents=True, exist_ok=True)
    (step_dir / "plot_manifest.json").write_text(mjson, encoding="utf-8")
    print(f"[DONE] 出图归集 -> {out_root}（marker 同步 -> {step_dir / 'plot_manifest.json'}）")


if __name__ == "__main__":
    main()