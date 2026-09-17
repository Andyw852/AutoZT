#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step5_dielect.py —— DFPT 介电常数（step5_dielect）。

结构从优化结果接力，IBRION=8 + LEPSILON 一次微扰求 ε∞ 与 ε₀。
产出目录：step5_dielect/，判据看 OUTCAR 的 MACROSCOPIC STATIC DIELECTRIC TENSOR。
"""
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ke_common as kc
import stepconf  # noqa: E402
from dim_common import require_dim, resolve_tpl  # noqa: E402

# =========================== 可改参数区 ===========================
# ---------- 体系判别阻断（step2.15_discriminant）----------
# False = 默认拦截 SEMIMETAL/METAL@PBE（这套半导体框架不适用）。
# True  = 强制继续（金属体系也有人要算输运，或 PBE 误判而杂化还没来得及重判）。
# ★ 不是硬停：读不到 discriminant.json 时一律放行，行为与加这道闸门前完全一致。
FORCE_TRANSPORT = False


def _enforce_min_kz(kp, minimum=2):
    """3D 体相的 DFPT 介电不允许任何方向只取 1 个 k 点（2026-09-15 实测）。

    为什么必须卡：vaspkit 按 k-spacing 自动出网格时，面对 c 特别长的体相胞
    （例：A2B2Te5, a=b=4.395 A, c=17.44 A）会给出 9x9x1 —— 长轴方向只有 1 个 k 点。
    实测（Pb2Sb2Te5, PBE+SOC, ISYM=-1, 其余参数完全一致，只改网格）：
        9x9x1  (NKPTS= 81) -> eps_inf = 100.73 / 51.18   <- 生产用的就是这个
        9x9x2  (NKPTS=162) -> eps_inf =  91.80 / 55.41
       12x12x2 (NKPTS=288) -> eps_inf =  91.67 / 55.11   <- 与上一行差 0.1%/0.5%，已收敛
    即 kz=1 让面内 eps 高估 9.0%、面外低估 7.7%。**把 kz 从 1 提到 2 就收敛了**，
    面内网格不用动（9->9）。代价是 k 点数翻倍，远小于盲目加密面内网格。

    注意：真正的 2D 滑板（有真空层）kz=1 才是对的，那条路由 dim=2d 的
    force_kz1_2d 处理；本函数只作用于 3D 体相，不碰 2D。
    """
    try:
        lines = kp.read_text(encoding="utf-8").splitlines()
    except OSError:
        return
    for i, ln in enumerate(lines):
        nums = ln.split()
        if len(nums) == 3:
            try:
                mesh = [int(v) for v in nums]
            except ValueError:
                continue
            if min(mesh) == 1:
                new = [max(v, minimum) for v in mesh]
                lines[i] = "  %d   %d   %d" % tuple(new)
                kp.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
                print("[OK] KPOINTS %s -> %s：3D 体相 DFPT 介电强制 kz>=%d"
                      "（kz=1 实测使 eps 偏差 8~9%%）" % (mesh, new, minimum))
            return

def _fix_soc_isym(incar):
    """SOC 体系的 DFPT 必须 ISYM=-1（VASP 6.4.3/6.5.0 spinor 旋转缺失 bug）。

    实测（A2B2Te5 四体系, PBE+SOC, 窄隙 0.15~0.23 eV, 2026-09-15 复核）：
      ISYM=2  + LSORBIT=T + IBRION=8 + KPAR=8 -> DFPT 崩溃 / 张量非对角项异常
                                                 （本应对角，却达 0.1~0.4 量级）
      ISYM=-1 + LSORBIT=T + IBRION=7 + KPAR=1 -> 四材料全部产出正确 ε∞
    正确性证据：Born 有效电荷声学求和规则 Sigma Z* = 0 精确到 1e-5；
                与同一次运行的独立粒子 ε∞ 面内一致到 1.0~1.4%。

    本函数在 inherit_scf_tags 之后调用：继承来的 LSORBIT=.TRUE. 会把模板里的
    ISYM=2 顶掉（或模板值被继承值覆盖），所以必须在这一步之后做最终裁决。
    非 SOC 体系原样返回，不改变既有行为。
    """
    try:
        txt = incar.read_text(encoding="utf-8")
    except OSError:
        return
    if not re.search(r"^\s*LSORBIT\s*=\s*\.TRUE\.", txt, re.M | re.I):
        return
    new, n = re.subn(r"^\s*ISYM\s*=.*$", "ISYM   = -1", txt, count=1, flags=re.M)
    if n == 0:
        new = txt.rstrip("\n") + "\nISYM   = -1\n"
    incar.write_text(new, encoding="utf-8", newline="\n")
    print("[..] SOC 体系：ISYM 已强制为 -1（DFPT+SOC 的 spinor 旋转 bug 绕法）；"
          "代价是无对称约化、k 点数按全 BZ 计，内存与耗时相应增加。")

def _disc_gate():
    """体系判别闸门：默认拦截 SEMIMETAL/METAL@PBE，FORCE_TRANSPORT 可覆盖。"""
    import sys as _s
    from pathlib import Path as _P
    try:
        _s.path.insert(0, str(_P(__file__).resolve().parent))
        import discriminant_common as _dc
        if not _dc.gate(_P.cwd(), "step5_dielect", force=FORCE_TRANSPORT):
            _s.exit("[BLOCKED] 体系判别为 SEMIMETAL/METAL@PBE —— step5_dielect 已阻断。"
                    "完整提示见 step2_bandgap/step2.15_discriminant/discriminant.json "
                    "的 block.hint；强制继续请把本脚本顶部 FORCE_TRANSPORT 设为 True。")
    except SystemExit:
        raise
    except Exception as _e:
        print("[WARN] 体系判别闸门异常，放行：%s" % _e, file=_s.stderr)


OUTDIR_NAME  = "step5_dielect"
PREV_CANDS   = ["step1_opt", "step1_std_opt"]
DIMENSION    = "auto"
VASPKIT_EXE  = "vaspkit"
KSCHEME      = "2"
# KSPACING：DFPT 的 k 点网格间距（交给 vaspkit 的 KSCHEME=2 生成）。
# ★ 2026-09-14 由 0.04 收紧到 0.02，两条实测理由：
#   (1) LPEAD 已默认关掉（见 incar_dfpt_*.tpl）。LPEAD=.TRUE. 的本职就是**加速 k 收敛**，
#       关掉它等于用更密的网格换数值稳定性 —— 所以 KSPACING 必须比原来更严，不能一样严。
#       Si 实测：8x8x8 时 T 路偏 2.5%、F 路偏 8.8%；到 16x16x16 两路才收敛到同值
#       （13.390 vs 13.422，差 0.24%）—— 即 F 路需要明显更密的网格。
#   (2) ε∞ 的 k 收敛速度由带隙决定（响应分母是 ε_c−ε_v，隙越小费米面附近贡献越尖锐）。
#       Si 是 PBE 0.6 eV 的宽隙立方体系，属最容易收敛的一端；本链的真实目标是
#       Pb2Sb2Te5 0.227 eV / Sn2Sb2Te5 0.216 eV / Mg2C60 0.196 eV。**拿 Si 标定的值
#       搬到窄隙体系会不够**，故取更严的一档 0.02（网格约翻倍）。
#   若某体系仍偏，按带隙做自适应或继续收紧；不要因为"和 uniform 一样稀"而放松。
KSPACING     = "0.02"
FUNC         = "inherit"      # patch_ke_dag: inherit=继承 step1
# VASP 的 vdW 修正不进 DFPT 响应（见 wiki IVDW 词条），默认把 D3 剥掉；
# 想强行保留（只影响总能，不影响 ε）把下面改 True。
KEEP_D3_IN_DFPT = False
MANUAL_ENCUT = None
ENCUT_FACTOR = 1.5
STEP_LABEL   = "S5_dielect"
# =================================================================
GGA_MAP = {"pbe": "PE", "pbesol": "PS", "pbe-d3": "PE"}

def main():
    _disc_gate()
    cwd = Path.cwd(); out = cwd / OUTDIR_NAME; out.mkdir(exist_ok=True)
    prev = kc.find_prev_dir(cwd, PREV_CANDS)
    if prev is None:
        sys.exit("[ERROR] 找不到含 CONTCAR 的上一步：%s" % PREV_CANDS)
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR", "step1_opt")
    _func, _subs = kc.resolve_func(prev, FUNC, OUTDIR_NAME,
                                   drop_d3=not KEEP_D3_IN_DFPT)
    dim = kc.read_method_dim(prev / kc.METHOD_FILE) \
        or kc.resolve_dim_for(out / "POSCAR", DIMENSION)[0]
    _, vac_axis = kc.resolve_dim_for(out / "POSCAR", dim)
    require_dim(dim, ('2d', '3d'), "step5_dielect",
                why="DFPT 给的是介电张量；分子对应的是极化率，定义和量纲都不同")
    print("[..] 维度：%s" % dim.upper())
    kc.write_method(out / kc.METHOD_FILE, dim, "DFPT 介电常数",
                    func=_func)
    kc.vaspkit_kpoints(out, KSCHEME, KSPACING, VASPKIT_EXE, dim, vac_axis)
    if dim == "3d":
        _enforce_min_kz(out / "KPOINTS", minimum=2)
    kc.vaspkit_potcar(out, VASPKIT_EXE)
    encut = MANUAL_ENCUT or kc.encut_from_potcar(out / "POTCAR", ENCUT_FACTOR)
    tpl = Path(__file__).resolve().parent / ("incar_dfpt_%s.tpl" % dim)
    if not tpl.is_file():
        sys.exit("[ERROR] 找不到模板 %s" % tpl.name)
    _sub = {"SYSTEM": cwd.name + " DFPT", "ENCUT": encut}
    _sub.update(_subs)
    kc.render_tpl(tpl, _sub, out / "INCAR")
    # DFPT(IBRION=8) 在多数 VASP 版本不支持 LDA+U，故 with_u=False
    kc.inherit_scf_tags(out / "INCAR", cwd, with_u=False, label="dielect")
    _fix_soc_isym(out / "INCAR")
    # step.conf 的 [incar]/[incar.final]/[incar.delete] 覆盖。
    # 此前本技能所有 gen 脚本只调 read_submit()，这三节**写了没人读**，
    # 用户在 step.conf 里改 INCAR 会被静默忽略。2026-09-15 在 step2.3_hse 上
    # 暴露并修复；这里改用 stepconf 里的唯一实现，避免各脚本各写一套。
    _ic_log = []
    stepconf.apply_incar_file(out / "INCAR", log=_ic_log)
    for _m in _ic_log:
        print("[..] %s" % _m)

    submit_tpl = resolve_tpl(Path(__file__).resolve().parent, "submit_std", dim)
    submit = out / "submit.sh"
    submit.write_text(submit_tpl.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    kc.patch_submit_jobname(submit, kc.new_jobname(cwd, STEP_LABEL))
    stepconf.apply_submit(submit, stepconf.read_submit(stepconf.CONF_NAME, used_incar=True))
    print("[DONE] %s：DFPT 输入就绪（KPAR=NCORE=1），可提交" % OUTDIR_NAME)

if __name__ == "__main__":
    main()
