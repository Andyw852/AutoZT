#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_nac.py —— NAC 介电/Born（step3_nac，可选步）。

在原胞上做 LEPSILON DFPT，得高频介电张量 ε∞ + Born 有效电荷，供声子/热导率的
非解析项修正（极性绝缘体 Γ 点 LO-TO 劈裂）。BORN 文件在 step5 从本步 vasprun 生成。

金属：step2 的带隙≈0 时本步会打印提示——金属无 LO-TO 劈裂、DFPT 介电发散，
应在项目配置写 nac: false 把本步整个关掉（tf 就不注入它）。若仍跑，step5 会校验
ε∞/Born 的物理性，不物理则自动退回无 NAC 出谱，不会污染结果。
产出目录：step3_nac/
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kl_common as kc
from dim_common import require_dim  # noqa: E402
import stepconf

OUTDIR = "step3_nac"
STEP   = "step3_nac"
PREV   = ["step1_std_opt"]

SPEC = {
    "FUNC":         ("pbesol", "str"),
    "KSPACING":     ("0.03",   "str"),
    "KSCHEME":      ("2",      "str"),
    "ENCUT":        (None,     "int"),
    "ENCUT_FACTOR": (1.5,      "float"),
    "VASPKIT_EXE":  ("vaspkit", "str"),
    "SKIP_IF_METAL": (True,    "bool"),   # 仅提示；实际跳过请在项目配置置 nac: false
}
# 结构体检第五条（B）：SYMMETRY_AUDIT/SYMMETRIZE 等键（见 kl_common.SYMMETRY_SPEC）
SPEC.update(kc.SYMMETRY_SPEC)


def check_metal(cwd):
    """若能读到 step2 的 vasprun，用参考引擎判金属，仅打印提示。"""
    vr = cwd / "step2_static" / "vasprun.xml"
    if not vr.is_file():
        return
    try:
        import lattice_kappa as lk
        is_metal, gap = lk.detect_is_metal(str(vr))
        if is_metal:
            print("[WARN] step2 带隙≈%.3f eV，判为金属 —— 金属无 LO-TO 劈裂、DFPT 介电"
                  "发散。建议在项目配置写 nac: false 关掉本步（tf 就不再注入 step3_nac）。"
                  % (gap or 0.0))
        elif gap is not None:
            print("[..] step2 带隙 %.3f eV（绝缘体/半导体），NAC 有意义" % gap)
    except Exception as e:
        print("[..] 金属判定跳过（%s）" % e)


def warn_2d_nac(dim):
    """2D + NAC：phono3py 只有 3D 方案，对真 2D 是近似；文献建议 2D 关 NAC。"""
    if str(dim).lower() != "2d":
        return
    print("[WARN] 2D 材料 + NAC：phono3py 只有 Wang/Gonze 两种 3D 方案，对真 2D 是近似——")
    print("       2D 极性材料的 LO-TO 劈裂在 q->0 应趋于零（V 形），3D 方案却给出随真空")
    print("       厚度变化的伪劈裂。文献（Sohier et al., Nano Lett. 2017）的实用处方是：")
    print("       2D 极性材料在 phono3py/phonopy 里【建议关掉 NAC】(nac: false)；")
    print("       严格的 2D-NAC 需用 Quantum ESPRESSO 的 2D 开边界 DFPT，phono3py 不带。")
    print("       本流程 step6 已默认对 2D 不施加 NAC（KAPPA_NAC=auto）；要强制用设 KAPPA_NAC=on。")


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    check_metal(cwd)

    prev = kc.find_prev_dir(cwd, PREV)
    if prev is None:
        sys.exit("[ERROR] 找不到 step1_std_opt 的结构")
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR", "step1_std_opt")
    # 结构体检第五条（B，2026-09-20）：S1 弛豫后的数值微畸变审计 + 可选对称化
    # （默认 SYMMETRY_AUDIT=warn 只告警、SYMMETRY_SYMMETRIZE=off 不改结构）
    kc.symmetry_gate(out / "POSCAR", conf)

    meth = kc.read_method(prev / kc.METHOD_FILE)
    dim = (meth.get("DIM", "").lower() or kc.resolve_dim(out / "POSCAR")[0])
    _, vac_axis = kc.resolve_dim(out / "POSCAR", dim)
    require_dim(dim, ('2d', '3d'), "step3_nac",
                why="NAC 修正的是 LO-TO 劈裂，孤立分子没有长程库仑的 q->0 行为")
    warn_2d_nac(dim)
    func = conf["FUNC"]
    if func in (None, "", "auto"):
        func = meth.get("FUNC", "pbesol").lower()
    # ★ 泛函一致性硬检查（wangchao review 二）：本步泛函必须与 S1 记录在
    #   workflow_method.txt 里的 FUNC 相同 —— 泛函不同就不在同一步势能面上，
    #   S3_nac 的产物会带上不匹配的力/应力，且完全静默。
    _rec_func = str(meth.get("FUNC") or "").strip().lower()
    if _rec_func and func != _rec_func:
        sys.exit("[ERROR] S3_nac 的泛函与 step1 不一致：本步 FUNC=%s，workflow_method.txt "
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
    nac_lines = ["", "# ---- NAC：LEPSILON DFPT 介电 + Born 有效电荷（gen_step3 注入）----",
                 "LEPSILON = .TRUE.", "LPEAD    = .TRUE.",
                 "IBRION   = -1", "NSW      = 0", "LWAVE    = .FALSE.", "LCHARG = .FALSE."]
    txt = (out / "INCAR").read_text(encoding="utf-8")
    txt = "\n".join(l for l in txt.splitlines()
                    if not l.strip().upper().startswith(("NPAR", "NCORE", "LEPSILON",
                                                         "LPEAD", "IBRION", "NSW")))
    (out / "INCAR").write_text(txt + "\n" + "\n".join(nac_lines) + "\n",
                               encoding="utf-8", newline="\n")
    print("[OK] INCAR 注入 LEPSILON/LPEAD，移除 NPAR/NCORE")

    submit_tpl = kc.resolve_submit(here, dim, "submit_std")
    kc.write_submit(submit_tpl, out / "submit.sh",
                    {"JOBNAME": kc.new_jobname(cwd, "S3nac")})
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    print("[DONE] %s：NAC(LEPSILON) 输入就绪，OUTCAR 出 MACROSCOPIC STATIC DIELECTRIC TENSOR"
          % OUTDIR)


if __name__ == "__main__":
    main()
