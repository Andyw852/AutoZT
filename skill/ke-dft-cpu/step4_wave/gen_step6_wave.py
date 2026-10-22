#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step4_wave.py —— amset wave → wavefunction.h5（step4_wave）。

amset wave 需要 uniform 步的 vasprun.xml 和 WAVECAR。WAVECAR 动辄几十 GB，
不复制，改用软链。提交脚本在计算节点 amset_clean 环境里跑 amset wave。
产出目录：step4_wave/，产物 wavefunction.h5
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stepconf  # noqa: E402
# 本步不跑 VASP，只软链 + 填提交脚本，不需要 ke_common 的 VASPKIT 部分
# 但复用它的 jobname/submit 工具
try:
    import ke_common as kc
    _HAS_KC = True
except Exception:
    _HAS_KC = False

# =========================== 可改参数区 ===========================
OUTDIR_NAME = "step4_wave"
UNIFORM_DIR = "step3_uniform"
STEP_LABEL  = "S4_wave"
# amset wave 命令（--planewave-cutoff 可按体系加大；见 amset 文档）
# patch_v63_port2（2026-09-23）：链首 rm -f wavefunction.h5 —— 判据是 marker(wavefunction.h5:)，
#   只看“文件在 + 含子串”；若 amset wave 失败而旧 h5 还在，整步会被判 OK、失败被隐藏（同 V63）。
# [patch_vasp_type-2026-10-22] 新版 pymatgen（实测 hfeshell amset 环境 2026.5.18）的 Wavecar
#   不再自动识别 vasp_type，而 amset 的 wavefunction/vasp.py 传的仍是 None，于是直接
#       ValueError: Incorrect vasp_type=None. Please open an issue if you are certain...
#   旧 pymatgen 环境不受影响。amset wave 本身有 --vasp-type std|ncl|gam 选项，按 uniform 步
#   INCAR 的 LSORBIT 判：.TRUE. -> ncl（vasp_ncl 写的非共线 WAVECAR），否则 std。
#   gam 不判：uniform 步一律 vasp_std / vasp_ncl，不会产出 vasp_gam 的 WAVECAR。
AMSET_CMD   = ("rm -f wavefunction.h5; amset wave --vasp-type %s "
               ">> amset.log 2>&1 && ls -l wavefunction.h5")
LINK_FILES  = ["vasprun.xml", "WAVECAR"]   # 从 uniform 步软链过来
# =================================================================


def _vasp_type(cwd):
    """[patch_vasp_type-2026-10-22] WAVECAR 的类型：看 uniform 步 INCAR 的 LSORBIT。

    新版 pymatgen 不再自动识别 vasp_type，必须显式给 amset wave --vasp-type。
    读不到 INCAR 时退回 "std"（保守：本技能 uniform 步默认就是 vasp_std）。
    """
    inc = Path(cwd) / UNIFORM_DIR / "INCAR"
    try:
        for ln in inc.read_text(errors="ignore").splitlines():
            s = ln.split("#", 1)[0]
            if "=" in s:
                k, v = s.split("=", 1)
                if k.strip().upper() == "LSORBIT" and v.strip().upper().startswith(".T"):
                    return "ncl"
    except Exception:
        pass
    return "std"


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    uni = cwd / UNIFORM_DIR
    # [patch_stale_input-2026-09-28] 本步重新生成 = 作业会 rm -f 旧 h5 再重算：
    #   下游 S8/S8.4（软链指向本步 h5 的那些）的 transport.json 一并失效、重新排队。
    #   V140：不再以"本步已有旧 h5"为条件 —— autozt rerun 先 rm -rf 本步目录，那样永远不触发。
    if _HAS_KC:
        kc.invalidate_downstream(cwd, OUTDIR_NAME, "%s 重新生成" % OUTDIR_NAME)

    for f in LINK_FILES:
        src = uni / f
        if not src.is_file():
            sys.exit("[ERROR] %s 缺失（uniform 步没跑完？）：%s" % (f, src))
        dst = out / f
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        # 相对软链，超算换挂载点也不断
        rel = os.path.relpath(src, out)
        dst.symlink_to(rel)
        print("[OK] 软链 %s -> %s" % (f, rel))

    # tf 把 submit_amset.tpl 与本脚本一起推到 gen 运行目录，但按原名推、不会改成
    # submit.sh；本步自己把它渲染成 out/submit.sh（维度步靠各自 gen 的 render，这里同理）。
    here = Path(__file__).resolve().parent
    tpl = next((p for p in (here / "submit_amset.tpl", cwd / "submit_amset.tpl")
                if p.is_file()), None)
    if tpl is None:
        sys.exit("[ERROR] 找不到 submit_amset.tpl（gen_need 里要有它，且应随 gen 脚本一起推送）")
    submit = out / "submit.sh"
    text = tpl.read_text(encoding="utf-8")
    _tpl_raw = text
    jobname = ("%s-ke-dft-cpu-%s" % (cwd.name, STEP_LABEL)) if not _HAS_KC \
        else kc.new_jobname(cwd, STEP_LABEL)
    _amset_env = kc.amset_env_name(cwd) if _HAS_KC else sys.exit("[ERROR] 无法 import ke_common 且 step.conf 缺 AMSET_ENV；请写 AMSET_ENV=<本集群 0.5.1 环境名>（jzzn/hanhai25=amset051、a800=amset_env、3090/hfeshell=amset）")
    _vt = _vasp_type(cwd)
    print("[OK] WAVECAR 类型：%s（--vasp-type %s）" % (_vt, _vt))
    text = (text.replace("{{JOBNAME}}", jobname)
                .replace("{{AMSET_CMD}}", AMSET_CMD % _vt)
                .replace("{{AMSET_ENV}}", _amset_env))
    if "{{AMSET_ENV}}" in text:
        sys.exit("[ERROR] submit_amset.tpl 的 {{AMSET_ENV}} 未填充（step.conf 缺 AMSET_ENV？）")
    if _HAS_KC:   # [V159] 激活的环境必须是 AMSET_ENV（旧项目级模板写死 amset_clean = 0.4.19）
        kc.check_amset_submit(_tpl_raw, text, _amset_env, tpl.name)
    submit.write_text(text, encoding="utf-8", newline="\n")
    stepconf.apply_submit(submit, stepconf.read_submit(stepconf.CONF_NAME))
    print("[OK] submit.sh 填好 amset 命令")
    print("[DONE] %s：软链就绪，提交后 amset wave 产出 wavefunction.h5" % OUTDIR_NAME)


if __name__ == "__main__":
    main()
