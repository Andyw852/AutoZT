#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_phonon.py —— phonopy 收力 + fc2 拟合 + 声子谱（phonon-dft-cpu S3）。

接力 step2 的位移+力，写 submit.sh 把 phonopy 收力 + fc2 拟合 + 声子谱放到计算节点跑。
gen 期先做两道上游校验（对齐 kl-dft-cpu S5_fc）：
  · 抽帧一致性：vasprun.xml 坐标必须等于 phonopy_disp.yaml 的「超胞 + 位移」；
  · SCF 收敛门禁：NELM 截断 / 未收敛的帧力不可信，拒绝拟合。
并把 IMAG_*/ZA_* 阈值、CONDA_SH/CONDA_ENV 真正写进 kl_params / submit.sh（旧版
step.conf 的 IMAG_THR 因没写进 kl_params 而被驱动忽略）。
产出（作业跑完后）：phonon_summary.json（'"stable": true' 为 marker）+ band-dft-cpu.yaml。
"""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kl_common as kc
import stepconf

OUTDIR = "step3_phonon"
STEP   = "step3_phonon"
DISP   = "step2_disp"

SPEC = {
    "FUNC":        ("auto", "str"),   # 继承 step1 的方法卡
    "BAND_POINTS": (51,   "int"),
    # 虚频判据（唯一真源 skill/_common/imag_policy.py）：IMAG_THR 语义是「近 Γ 声学支上限」
    "IMAG_THR":          (0.10, "float"),
    "IMAG_THR_STRICT":   (0.10, "float"),
    "IMAG_QGAMMA":       (0.05, "float"),
    "IMAG_QGAMMA_GRACE": (1.2,  "float"),
    # 2D ZA 弯曲支二次性（ω∝q^p，要求 1.7<p<2.3）；auto = 2D 打开
    "ZA_CHECK":    ("auto", "str"),
    "ZA_QMAX":     (0.05, "float"),
    # 集群 conda（由 autozt 从 setting/<集群>.yaml 注入 step.conf；读不到报错）
    "CONDA_SH":    ("", "str"),
    "CONDA_ENV":   ("", "str"),
}


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    disp = cwd / DISP

    if not (disp / "phonopy_disp.yaml").is_file():
        sys.exit("[ERROR] %s 缺 phonopy_disp.yaml（step2 未生成位移）" % disp)
    if not list(disp.glob("disp-*/vasprun.xml")):
        sys.exit("[ERROR] %s 下无 disp-*/vasprun.xml，位移单点还没算完" % disp)

    # ---- 抽帧一致性：力必须与位移对应（对齐 kl-dft-cpu S5_fc）----
    _ok, _note = kc.check_frames_match_displacements(disp, disp_yaml="phonopy_disp.yaml")
    print("[%s] S2 帧一致性：%s" % ("OK" if _ok else "FAIL", _note))
    if not _ok:
        sys.exit("[ERROR] %s\n        请清空 step2_disp 的 disp-*/POSCAR-*/phonopy_disp.yaml/"
                 "SPOSCAR 后重跑 S2（或 -j S2_disp rerun）。" % _note)

    # ---- SCF 收敛门禁：NELM 截断 / 未收敛帧的力不可信，拒绝拟合 ----
    _scf_ok, _scf = kc.check_outcar_scf_convergence(disp)
    print(kc.format_scf_report(_scf))
    (out / "scf_steps.json").write_text(
        json.dumps(_scf, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    if not _scf_ok:
        _bad = "、".join(b["frame"] for b in _scf["bad_frames"][:20])
        _more = "（共 %d 帧）" % _scf["n_bad"] if _scf["n_bad"] > 20 else ""
        sys.exit(
            "[ERROR] S2 有 %d/%d 帧 SCF 未正常收敛（NELM 截断 / 无 aborting loop / 帧未算完），"
            "力不可信，拒绝拟合：\n        作废帧：%s%s\n"
            "        这些帧的力是电子步没收敛的结果，拟合 fc2 会整体错掉而下游检查正常。\n"
            "        处理：对作废帧调 INCAR（提高 NELM / 换 ALGO / 调 AMIX,BMIX）后 "
            "autozt -tt phonon-dft-cpu -p <材料> -j S2_disp retry 只补这些帧；\n"
            "        每帧电子步数见 %s。"
            % (_scf["n_bad"], _scf["n_frames"], _bad, _more, out / "scf_steps.json"))

    for f in ("POSCAR", "SPOSCAR", "phonopy_disp.yaml", kc.KL_PARAMS, kc.METHOD_FILE):
        src = disp / f
        if src.is_file():
            shutil.copyfile(str(src), str(out / f))

    # ---- kl_params 继承 + 把 S3 的阈值/网格真正写进去（旧版漏写 → step.conf 失效）----
    inherited = kc.read_kl_params(disp / kc.KL_PARAMS)
    inherited.update({
        "IMAG_THR": conf["IMAG_THR"],
        "IMAG_THR_STRICT": conf["IMAG_THR_STRICT"],
        "IMAG_QGAMMA": conf["IMAG_QGAMMA"],
        "IMAG_QGAMMA_GRACE": conf["IMAG_QGAMMA_GRACE"],
        "ZA_CHECK": conf["ZA_CHECK"],
        "ZA_QMAX": conf["ZA_QMAX"],
        "BAND_POINTS": int(conf["BAND_POINTS"]),
    })
    kc.write_kl_params(out / kc.KL_PARAMS, **inherited)

    here = Path(__file__).resolve().parent
    if not (here / "phonon_fit_driver.py").is_file():
        sys.exit("[ERROR] 缺 phonon_fit_driver.py —— gen_need 里漏了它？")
    shutil.copyfile(str(here / "phonon_fit_driver.py"), str(out / "phonon_fit_driver.py"))
    # 作业里的驱动 import kl_common / za_2d / imag_policy（kl_common 又 import dim_common）。
    # gen_need 只保证 gen 本地能用到，不会自动进到发往计算节点的 step 目录，必须显式拷。
    for _dep in ("kl_common.py", "dim_common.py", "za_2d.py", "imag_policy.py"):
        if (here / _dep).is_file():
            shutil.copyfile(str(here / _dep), str(out / _dep))
        else:
            print("[WARN] 找不到 %s（gen_need 里漏了？）—— 作业可能 import 失败" % _dep)

    # ---- 集群 conda：读不到就报错，不兜底个人路径/环境名 ----
    if not str(conf["CONDA_SH"] or "").strip():
        sys.exit("[ERROR] 未配置 CONDA_SH —— 请在 setting/<集群>.yaml 填 conda_sh"
                 "（conda.sh 绝对路径），autozt 会注入本步 step.conf")
    if not str(conf["CONDA_ENV"] or "").strip():
        sys.exit("[ERROR] 未配置 CONDA_ENV —— 请在 setting/<集群>.yaml 填 conda_env"
                 "（phonopy/atomate2_p_a 所在环境名）")

    dim = (inherited.get("DIM") or "").lower()
    tpl = kc.resolve_submit(here, dim or "3d", "submit_fit")
    kc.write_submit(tpl, out / "submit.sh",
                    {"JOBNAME": "%s-phonon-dft-cpu-S3" % cwd.name,
                     "CONDA_SH": str(conf["CONDA_SH"]),
                     "CONDA_ENV": str(conf["CONDA_ENV"])})
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    print("[DONE] %s：submit.sh 就绪（作业跑完写 phonon_summary.json + band-dft-cpu.yaml）"
          % OUTDIR)


if __name__ == "__main__":
    main()
