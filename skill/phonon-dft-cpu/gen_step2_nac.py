#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_nac.py —— NAC 介电/Born（step2_nac，可选步；对齐 kl-dft-cpu 的 S3_nac）。

在原胞上做 LEPSILON DFPT，得高频介电张量 ε∞ + Born 有效电荷，供声子谱的非解析项
修正（极性绝缘体 Γ 点 LO-TO 劈裂）。S3_phonon 的驱动从本步 vasprun.xml 用
phonopy-vasp-born 生成 BORN 并加载；没有本步（nac: false）或产物不物理 → 自动退回
无 NAC 出谱，不会污染结果。

金属：金属无 LO-TO 劈裂、DFPT 介电发散，请在项目配置写 nac: false 关掉本步
（autozt 就不再注入它）。
产出目录：step2_nac/
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kl_common as kc
from dim_common import require_dim  # noqa: E402
import stepconf

OUTDIR = "step2_nac"
STEP   = "step2_nac"
PREV   = ["step1_std_opt"]

SPEC = {
    "FUNC":         ("pbesol", "str"),
    "KSPACING":     ("0.03",   "str"),
    "KSCHEME":      ("2",      "str"),
    "ENCUT":        (None,     "int"),
    "ENCUT_FACTOR": (1.5,      "float"),
    "VASPKIT_EXE":  ("vaspkit", "str"),
}
# 结构体检第五条（B）：SYMMETRY_AUDIT/SYMMETRIZE 等键（见 kl_common.SYMMETRY_SPEC）
SPEC.update(kc.SYMMETRY_SPEC)


def warn_2d_nac(dim):
    """2D + NAC：phonopy 只有 3D 方案，对真 2D 是近似；LO-TO 在 q->0 应趋于零（V 形）。"""
    if str(dim).lower() != "2d":
        return
    print("[WARN] 2D 材料 + NAC：phonopy 只有 3D（Wang/Gonze）方案，对真 2D 是近似 ——")
    print("       2D 极性材料的 LO-TO 劈裂在 q->0 应趋于零（V 形），3D 方案却给出随真空")
    print("       厚度变化的伪劈裂。文献（Sohier et al., Nano Lett. 2017）的实用处方是：")
    print("       2D 极性材料在 phonopy 里【建议关掉 NAC】(nac: false)；严格 2D-NAC 需 QE 的")
    print("       2D 开边界 DFPT。要严格用 2D 方案请把本项目 nac 置 false 或改用 kl-dft 的")
    print("       处理口径。S3 仍会同时留一份无 NAC 的 band_noNAC.yaml 供对照。")


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)

    prev = kc.find_prev_dir(cwd, PREV)
    if prev is None:
        sys.exit("[ERROR] 找不到 step1_std_opt 的结构")
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR", "step1_std_opt")
    # 结构体检第五条：S1 弛豫后的数值微畸变审计 + 可选对称化
    kc.symmetry_gate(out / "POSCAR", conf)

    meth = kc.read_method(prev / kc.METHOD_FILE)
    dim = (meth.get("DIM", "").lower() or kc.resolve_dim(out / "POSCAR")[0])
    _, vac_axis = kc.resolve_dim(out / "POSCAR", dim)
    require_dim(dim, ("2d", "3d"), "step2_nac",
                why="NAC 修正的是 LO-TO 劈裂，孤立分子没有长程库仑的 q->0 行为")
    warn_2d_nac(dim)
    func = conf["FUNC"]
    if func in (None, "", "auto"):
        func = meth.get("FUNC", "pbesol").lower()
    # 泛函一致性硬检查：与 S1 workflow_method.txt 记录必须相同（同一势能面）
    _rec_func = str(meth.get("FUNC") or "").strip().lower()
    if _rec_func and func != _rec_func:
        sys.exit("[ERROR] S2_nac 的泛函与 step1 不一致：本步 FUNC=%s，workflow_method.txt "
                 "记 FUNC=%s。\n        两者必须相同（同一势能面）；拒绝生成。\n"
                 "        处置：把本步 FUNC 置 auto（推荐，自动继承），或先核对 S1 的泛函。"
                 % (func, _rec_func))

    kc.vaspkit_kpoints(out, conf["KSCHEME"], conf["KSPACING"],
                       conf["VASPKIT_EXE"], dim, vac_axis)
    kc.vaspkit_potcar(out, conf["VASPKIT_EXE"])
    encut = conf["ENCUT"] or kc.encut_from_potcar(out / "POTCAR", conf["ENCUT_FACTOR"])

    here = Path(__file__).resolve().parent
    incar_tpl = kc.resolve_submit(here, dim, "incar")
    subs = {"SYSTEM": "%s nac(LEPSILON)" % cwd.name, "ENCUT": encut,
            "GGA": kc.GGA_MAP.get(func, "PS"),
            "VDW_LINE": ("IVDW = %s" % kc.VDW_MAP[func]) if kc.VDW_MAP.get(func) else "# no vdW"}
    kc.render_tpl(incar_tpl, subs, out / "INCAR")

    # NAC 关键标签：LEPSILON DFPT 出 ε∞ + Born；LPEAD 提升数值稳定；NPAR/NCORE 与 LEPSILON 不兼容
    nac_lines = ["", "# ---- NAC：LEPSILON DFPT 介电 + Born 有效电荷（gen_step2_nac 注入）----",
                 "LEPSILON = .TRUE.", "LPEAD    = .TRUE.",
                 "IBRION   = -1", "NSW      = 0", "LWAVE    = .FALSE.", "LCHARG = .FALSE."]
    txt = (out / "INCAR").read_text(encoding="utf-8")
    txt = "\n".join(l for l in txt.splitlines()
                    if not l.strip().upper().startswith(("NPAR", "NCORE", "LWAVE", "LCHARG",
                                                         "LEPSILON", "LPEAD", "IBRION", "NSW")))
    (out / "INCAR").write_text(txt + "\n" + "\n".join(nac_lines) + "\n",
                               encoding="utf-8", newline="\n")
    print("[OK] INCAR 注入 LEPSILON/LPEAD，移除 NPAR/NCORE")

    submit_tpl = kc.resolve_submit(here, dim, "submit_std")
    kc.write_submit(submit_tpl, out / "submit.sh",
                    {"JOBNAME": kc.new_jobname(cwd, "S2nac")})
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    print("[DONE] %s：NAC(LEPSILON) 输入就绪，OUTCAR 出 MACROSCOPIC STATIC DIELECTRIC TENSOR"
          % OUTDIR)


if __name__ == "__main__":
    main()
