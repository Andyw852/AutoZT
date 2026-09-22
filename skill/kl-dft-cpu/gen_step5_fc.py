#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step5_fc.py —— 力常数拟合，提交计算节点（step5_fc）。

【结构变更】S5_fc 从"登录节点 gen 里裸跑 phono3py"改成"提交计算节点作业"：pheasy/symfc
拟合是重活（几十核 + 大内存 + 数小时），压不进登录节点。本 gen 只做准备：
  1. 校验 step4 产物（phono3py_disp.yaml + disp-*/vasprun.xml；alm 才有 SPOSCAR）
  2. 把 step.conf + kl_params 解析成 fit_config.json（作业里 kl_fc_backends.py 读它）
  3. 按 FIT_ENGINE 选提交模板渲染 submit.sh
tf 提交后，计算节点按 submit.sh 依次跑 kl_fc_backends 的 prep → 拟合 → post：
  拟合器（FIT_ENGINE）：phono3py(symfc/alm，默认) | pheasy(随机位移压缩感知，需 METHOD=alm)
  产出：step5_fc/phono3py/（fc2/fc3.hdf5 + phono3py_disp.yaml）
        step5_fc/shengbte/（FORCE_CONSTANTS_2ND/3RD，EXPORT_SHENGBTE=true 时）
        step5_fc/phonon_summary.json（虚频闸 marker：'"stable": true'）
产出目录：step5_fc/
"""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kl_common as kc
import stepconf

OUTDIR   = "step5_fc"
STEP     = "step5_fc"
DISP_DIR = "step4_disp"

SPEC = {
    "FUNC":        ("pbesol", "str"),   # 全局带入，本步不用
    # —— 拟合器选择 ——
    # auto：按维度选 —— 2D 用 pheasy，3D 用 phono3py（symfc 快且稳，3D 不需要旋转不变性）。
    #   ★ 限制在"必须走 pheasy"这一步，与用哪种回归方法无关：RASR 是在 pheasy 的 -c
    #   （零空间构造）步写进零空间的，之后 OLS / RFE / LASSO 都能用。所以引擎默认改成
    #   pheasy 才是 P0-1 的关键，PHEASY_FIT_METHOD 是另一个独立选择（见下）。
    "FIT_ENGINE":  ("auto", "str"),     # auto | phono3py | pheasy
    # phono3py 路（symfc/alm）
    "FC_CALC":     ("symfc",  "str"),   # symfc | alm
    "FC3_CUTOFF":  (None,     "str"),   # fc3 截断 Å（"5.0"）；空/None=不截断
    # pheasy 路
    # auto：2D → LASSO、3D → RFE。与 RASR 无关（RASR 只看引擎是不是 pheasy）；
    #   选 LASSO 是因为对照实验里 RFE/OLS 的 fc3 幅值本就贴近参考值，而 LASSO 必须靠
    #   --std + 去偏才能压到 0.13% 以内 —— 那两项已写死在 submit_fit_pheasy.tpl 里，
    #   用 LASSO 前务必确认模板没被项目级副本遮蔽。
    "PHEASY_FIT_METHOD": ("auto", "str"), # auto | LASSO | RFE | OLS（OLS 最吃内存）
    "PHEASY_C3_CUTOFF":  ("5.2", "str"), # pheasy fc3 截断 Å；None=不截断
    "PHEASY_ENABLE_FC":  (3,     "int"), # 2|3|4（热导率需 ≥3）
    # 旋转不变性/平衡条件（RASR）：auto = 2D 用 BHH、3D 不加（文献结论：对体材料可忽略，
    #   对 2D 是硬要求——不加则 ZA 近 Γ 线性化甚至出虚频）。BHH | BH | H | none 可强制。
    #   ★ 必须施加在 pheasy 的 -c（零空间构造）步；-f 步读的是 ns_*.npz，--rasr 在 -f 上无效。
    "PHEASY_RASR":       ("auto", "str"),
    "PHEASY_BIN":        ("pheasy", "str"), # pheasy 可执行名：pheasy | pheasy-gpu（GPU 版）
    "NULL_SPACE_EPS":    (0.001, "float"),
    # —— 导出 & 虚频闸 ——
    # —— 缺帧容错（仅随机位移 METHOD=alm 生效；findiff 必须帧帧齐全）——
    "MIN_SUCCESS_RATIO":  (0.9, "float"),  # 成功帧占比下限，低于它报错
    "MIN_SUCCESS_FRAMES": (0,   "int"),    # 成功帧绝对下限，0=不限
    # —— 导出 & 虚频闸 ——
    "EXPORT_SHENGBTE": (True, "bool"),   # 任一拟合器都产出 shengbte 力常数
    "BAND_POINTS":     (51,   "int"),
    "IMAG_THR":        (0.10, "float"),  # 虚频阈值(THz)
    # P2-2：2D 的 ZA 弯曲支二次性检查（ω ∝ q^p，要求 1.7<p<2.3）。
    #   auto = 2D 打开；只看最小频率不够 —— ZA 线性化时频率全是正的，照样过虚频闸门，
    #   但 κ 会整体错掉。检查结果写进 phonon_summary.json 的 za_exponent。
    "ZA_CHECK":        ("auto", "str"),   # auto | on | off
    "ZA_QMAX":         (0.05, "float"),   # 拟合用 q 上限（倒格子约化单位，0.05~0.1）
    # 作业资源（核数/qos/时长）走 step.conf 的 [submit] 段覆盖 #SBATCH，不在 [params] 里。
}


_RASR_VALUES = ("BHH", "BH", "H")


def _resolve_rasr(val, dim):
    """PHEASY_RASR 解析：auto → 2D=BHH / 3D=none；显式值原样（none 归一成小写）。

    pheasy 的 --rasr 只接受 BH / H / BHH（basic_io.py 的 choices），所以 "none" 不是
    一个能传的值 —— 关掉就是不传这个开关（模板里 RASR_FLAGS 为空）。
    """
    v = str(val if val is not None else "auto").strip().upper()
    if v == "AUTO":
        return "BHH" if dim == "2d" else "none"
    if v in _RASR_VALUES:
        return v
    if v in ("NONE", "OFF", "FALSE", "F", "0", ""):
        return "none"
    sys.exit("[ERROR] PHEASY_RASR=%r 非法：只允许 auto | BHH | BH | H | none" % val)


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    disp = cwd / DISP_DIR

    # ---- 校验 step4 产物 ----
    if not (disp / "phono3py_disp.yaml").is_file():
        sys.exit("[ERROR] %s 缺 phono3py_disp.yaml（step4 未生成位移）" % disp)
    if not list(disp.glob("disp-*/vasprun.xml")):
        sys.exit("[ERROR] %s 下无 disp-*/vasprun.xml，位移单点还没算完" % disp)

    # ---- 抽帧校验（2026-09-17，wangchao 要求）：力必须与位移对应 ----
    #   S4 是 fanout，retry 只补缺失帧；若 S1 换过结构而旧 disp-* 残留，会出现
    #   "旧结构的力 + 新位移"的静默错配。这里抽 3 帧把 vasprun.xml 的坐标与
    #   "phono3py_disp.yaml 的超胞 + 位移"逐原子比对（周期回绕后 < 1e-3 Å）。
    _ok, _note = kc.check_frames_match_displacements(disp)
    print("[%s] S4 帧一致性：%s" % ("OK" if _ok else "FAIL", _note))
    if not _ok:
        sys.exit("[ERROR] %s\n        请清空 step4_disp 的 disp-*/POSCAR-*/phono3py_disp.yaml/SPOSCAR "
                 "后重跑 S4（或 -j S4_disp rerun）。" % _note)

    # ---- SCF 收敛门禁（2026-09-19，wangchao 要求）：NELM 截断/未收敛的帧不能拟合 ----
    #   VASP 撞 NELM 时照样输出力、作业正常退出，只在 OUTCAR 留一段
    #   "number of steps (NELM) ... forces ... might not be reliable"；这种帧拿去拟合
    #   会让 fc2/fc3 与 κ 整体错掉，而下游所有检查都显示正常 —— 正是要拦的静默错误。
    #   反向：静态单点正常收敛时 OUTCAR 必有且只有 1 次 "aborting loop because EDIFF is reached"。
    _scf_ok, _scf = kc.check_outcar_scf_convergence(disp)
    print(kc.format_scf_report(_scf))
    (out / "scf_steps.json").write_text(
        json.dumps(_scf, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    if not _scf_ok:
        _bad = "、".join(b["frame"] for b in _scf["bad_frames"][:20])
        _more = "（共 %d 帧）" % _scf["n_bad"] if _scf["n_bad"] > 20 else ""
        sys.exit(
            "[ERROR] S4 有 %d/%d 帧 SCF 未正常收敛（NELM 截断 / 无 aborting loop / 帧未算完），"
            "力不可信，拒绝拟合：\n        作废帧：%s%s\n"
            "        这些帧的力是电子步没收敛的结果，拟合 fc2/fc3 会整体错掉而下游检查正常。\n"
            "        处理：对作废帧调 INCAR（提高 NELM / 换 ALGO / 调 AMIX,BMIX）后 "
            "autozt -tt kl-dft-cpu -p <材料> -j S4_disp retry 只补这些帧；\n"
            "        每帧电子步数见 %s。"
            % (_scf["n_bad"], _scf["n_frames"], _bad, _more, out / "scf_steps.json"))

    params = kc.read_kl_params(disp / kc.KL_PARAMS)
    method = (params.get("METHOD") or "alm").lower()
    supercell = params.get("SUPERCELL") or ""
    dim = (params.get("DIM") or "").lower()
    if dim not in ("2d", "3d") and (disp / "POSCAR").is_file():
        # kl_params 是 S4 写的（DIM 一定有）；真缺了就按结构现判，别把 2D 当 3D 静默放过
        try:
            dim = kc.resolve_dim(disp / "POSCAR", "auto")[0]
        except Exception:
            dim = ""

    engine = str(conf["FIT_ENGINE"] or "auto").strip().lower()
    if engine in ("auto", ""):
        engine = "pheasy" if dim == "2d" else "phono3py"
        print("[..] FIT_ENGINE=auto → %s（DIM=%s）" % (engine, dim or "?"))
    if engine not in ("phono3py", "pheasy"):
        sys.exit("[ERROR] FIT_ENGINE 只允许 auto / phono3py / pheasy")
    if engine == "pheasy" and method != "alm":
        sys.exit("[ERROR] FIT_ENGINE=pheasy 需要随机位移（step4 METHOD=alm）。\n"
                 "        findiff 请用 FIT_ENGINE=phono3py，或把 step4 改成 alm 重跑。")
    p_method = str(conf["PHEASY_FIT_METHOD"] or "auto").strip().upper()
    if p_method in ("AUTO", ""):
        p_method = "LASSO" if dim == "2d" else "RFE"
        print("[..] PHEASY_FIT_METHOD=auto → %s（DIM=%s）" % (p_method, dim or "?"))
    if p_method not in ("LASSO", "RFE", "OLS"):
        sys.exit("[ERROR] PHEASY_FIT_METHOD 只允许 auto / LASSO / RFE / OLS")
    p_bin = str(conf["PHEASY_BIN"] or "pheasy").lower()
    if p_bin not in ("pheasy", "pheasy-gpu"):
        sys.exit("[ERROR] PHEASY_BIN 只允许 pheasy / pheasy-gpu")
    if str(conf["FC_CALC"]).lower() not in ("symfc", "alm"):
        sys.exit("[ERROR] FC_CALC 只允许 symfc / alm")

    # ---- RASR（旋转不变性 + 零应力平衡条件）----
    # 2D 的 ZA 弯曲支 ω∝q² 由 Born-Huang 旋转不变性保证；不加时近 Γ 会线性化、
    #   常常还带小虚频，虚频闸可能过（频率全为正）但 κ 是错的，所以默认 2D 必加。
    rasr = _resolve_rasr(conf["PHEASY_RASR"], dim)
    if engine != "pheasy" and rasr != "none":
        print("[..] FIT_ENGINE=%s 不用 pheasy 的 RASR，忽略 PHEASY_RASR=%s" % (engine, rasr))
        rasr = "none"
    print("[..] RASR=%s（DIM=%s，PHEASY_RASR=%s）"
          % (rasr, dim or "?", conf["PHEASY_RASR"]))
    if dim == "2d" and rasr == "none":
        print("[WARN] 2D 体系关闭了 RASR（PHEASY_RASR=none）：ZA 近 Γ 可能线性化或出虚频，"
              "拟合出的力常数与 κ 不可信。除非在做对照实验，请设 PHEASY_RASR=auto/BHH。")
    elif dim != "2d" and rasr == "none":
        print("[..] 3D：按文献结论不施加 RASR（auto 的行为），要强制请设 PHEASY_RASR=BHH/BH/H")

    # kl_params 一并拷进 step5_fc（溯源/后续继承）
    for f in (kc.KL_PARAMS, kc.METHOD_FILE):
        if (disp / f).is_file():
            shutil.copyfile(disp / f, out / f)

    # ---- 写 fit_config.json（作业里读）----
    cfg = {
        "FIT_ENGINE": engine,
        "METHOD": method,
        "DIM": dim.upper(),
        "SUPERCELL": supercell,                       # 对角三整数，pheasy --dim / shengbte scell
        "FC_CALC": str(conf["FC_CALC"]).lower(),
        "FC3_CUTOFF": (None if conf["FC3_CUTOFF"] in (None, "", "None", "none")
                       else str(conf["FC3_CUTOFF"])),
        "PHEASY_RASR": rasr,
        "PHEASY_FIT_METHOD": p_method,
        "PHEASY_C3_CUTOFF": str(conf["PHEASY_C3_CUTOFF"]),
        "PHEASY_ENABLE_FC": int(conf["PHEASY_ENABLE_FC"]),
        "PHEASY_BIN": p_bin,
        "NULL_SPACE_EPS": float(conf["NULL_SPACE_EPS"]),
        "MIN_SUCCESS_RATIO": float(conf["MIN_SUCCESS_RATIO"]),
        "MIN_SUCCESS_FRAMES": int(conf["MIN_SUCCESS_FRAMES"]),
        "EXPORT_SHENGBTE": bool(conf["EXPORT_SHENGBTE"]),
        "BAND_POINTS": int(conf["BAND_POINTS"]),
        "IMAG_THR": float(conf["IMAG_THR"]),
        "ZA_CHECK": str(conf["ZA_CHECK"]),
        "ZA_QMAX": float(conf["ZA_QMAX"]),
    }
    (out / "fit_config.json").write_text(
        json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")

    # ---- 选模板渲染 submit.sh ----
    here = Path(__file__).resolve().parent
    # 关键：把作业运行时要执行的驱动脚本拷进产出目录。gen_need 只保证 gen 本地能用到它，
    #   不会自动进到发往计算节点的 step 目录；作业里要 `python kl_fc_backends.py`，
    #   必须显式拷过去（否则计算节点报 No such file）。
    shutil.copyfile(here / "kl_fc_backends.py", out / "kl_fc_backends.py")
    # kl_fc_backends 的 band_path_2d / 出图段是【懒】import kl_common（kl_common 又 import dim_common）。
    # 不把这两个也拷过去，2D 的 band-dft-cpu.yaml 会静默跳过（2026-09-21 实测：
    # "No module named 'kl_common'" —— 只影响出图、不影响 ZA/虚频判据，但那样 band_path_2d 的
    # 2D 路径修复就没被真正执行）。作业目录 = out，所以直接拷到 out。
    # za_2d.py：kl_fc_backends 现在【模块级】import za_2d（三副本合并后），
    # 作业目录少了它会在 import 阶段直接 ModuleNotFoundError。
    for _dep in ("kl_common.py", "dim_common.py", "za_2d.py"):
        if (here / _dep).is_file():
            shutil.copyfile(here / _dep, out / _dep)
    if engine == "pheasy":
        # GPU 拟合（PHEASY_BIN=pheasy-gpu）走独立模板 submit_fit_pheasy_gpu（--gres + 降核）。
        # 纯 CPU 集群（jzzn/hanhai25）没有该模板 → gen 期即报错，不排进队才失败（G2）。
        _kind = ("submit_fit_pheasy_gpu" if p_bin == "pheasy-gpu"
                 else "submit_fit_pheasy")
        try:
            tpl = kc.resolve_submit(here, dim or "3d", _kind)
        except SystemExit:
            if _kind == "submit_fit_pheasy_gpu":
                sys.exit("[ERROR] PHEASY_BIN=pheasy-gpu 需要 GPU 拟合模板 "
                         "submit_fit_pheasy_gpu.tpl（a800/3090 已配）。\n"
                         "        当前集群只有 CPU 模板——pheasy-gpu 在无 GPU 节点跑不了。\n"
                         "        改用 PHEASY_BIN=pheasy（CPU 拟合），或把材料 hpc 切到 "
                         "a800/3090。")
            raise
        subs = {"JOBNAME": kc.new_jobname(cwd, "S5fit"),
                "DIM": supercell or "1 1 1",
                "FIT_METHOD": cfg["PHEASY_FIT_METHOD"],
                "ENABLE_FC": str(cfg["PHEASY_ENABLE_FC"]),
                "PHEASY_BIN": p_bin,
                "RASR": rasr,
                "C3_CUTOFF": cfg["PHEASY_C3_CUTOFF"],
                "NULL_SPACE_EPS": str(cfg["NULL_SPACE_EPS"])}
    else:
        tpl = kc.resolve_submit(here, dim or "3d", "submit_fit_p3py")
        subs = {"JOBNAME": kc.new_jobname(cwd, "S5fit")}
    kc.write_submit(tpl, out / "submit.sh", subs)
    stepconf.apply_submit(out / "submit.sh", conf.submit)

    print("[..] 拟合器=%s 方法=%s 超胞=%s DIM=%s" % (engine, method, supercell, dim or "?"))
    print("[DONE] %s：submit.sh + fit_config.json 就绪。tf 提交后计算节点出 fc2/fc3，"
          "写 phonon_summary.json（'\"stable\": true' 为 marker）。" % OUTDIR)


if __name__ == "__main__":
    main()