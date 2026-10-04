#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_fc.py —— kl-mlff S3_fc：拟合力常数（共享 fit-fc-thermal 引擎）。

把 step2_disp_force 的位移+力数据集交给共享引擎 _common/fcfit/gen_step1_fit，
写 fit_config.json + submit.sh；计算节点跑 fc_fit_driver.py prep|fit|post。
产出 step3_fc/{fc2,fc3}.hdf5 + phonon_summary.json + shengbte/（虚频闸 marker）。

导入模式（跨技能复用）：检测到推送来的 kl_bundle/（或 FC_IMPORT_DIR）→ 装配
step3_fc 布局、跳过拟合（见 reuse_fcfit.py）。
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# 共享引擎（_common/fcfit）+ stepconf/dim_common（_common/opt）+ _common
_base = Path(__file__).resolve().parent.parent
for _d in (_base / "fcfit", _base / "opt", _base):
    if _d.is_dir() and str(_d) not in sys.path:
        sys.path.append(str(_d))

import gen_step1_fit as g1  # noqa: E402
import stepconf             # noqa: E402

OUTDIR = "step3_fc"
STEP = "step3_fc"
SRC = "step2_disp_force"

# 共享引擎的参数面 + 本技能独有键。MACE_MODEL/DEVICE/DTYPE 本步不用，但全局
# step.conf 里有它们（step1/step2 用），声明出来避免 strict=warn 的未知键告警。
SPEC = dict(g1.SPEC)
SPEC.update({
    "MACE_MODEL": ("mace-mp:medium", "str"),
    "MACE_MODEL_DIR": ("", "str"),
    "DEVICE": ("auto", "str"),
    "DTYPE": ("float64", "str"),
    # MACE 给不出 Born 有效电荷：外部 DFPT 的 BORN 文件路径（空=不加 NAC）。
    "NAC_BORN": ("", "str"),
    # —— 由 fit-fc-thermal 产物导入（跨技能复用 fc2/fc3，对齐 kl-dft-cpu S5_fc）——
    "FC_IMPORT_DIR": ("", "str"),
    "FC_IMPORT_SUPERCELL": ("", "str"),
    "FC_IMPORT_DIM": ("", "str"),
    "FC_IMPORT_MESH": ("", "str"),
})


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP, strict="warn")

    # ---- fit-fc-thermal 产物导入模式（跨技能复用 fc2/fc3）----
    _fcimp = str(conf.get("FC_IMPORT_DIR") or "").strip()
    if _fcimp.lower() in ("off", "false", "0", "no", "none"):
        _fcimp = ""
    _bundle = cwd / "kl_bundle"
    if not _fcimp and _bundle.is_dir() and (_bundle / "fc2.hdf5").is_file():
        _fcimp = str(_bundle)
        print("[..] 检测到推送来的 kl_bundle/（fit-fc-thermal 力常数）→ 自动进入导入模式"
              "（想改为本技能自拟：FC_IMPORT_DIR = off）", flush=True)
    if _fcimp:
        _src = Path(_fcimp)
        if not _src.is_absolute():
            _src = cwd / _src
        if not (_src / "fc2.hdf5").is_file() or not (_src / "fc3.hdf5").is_file():
            sys.exit("[ERROR] FC_IMPORT_DIR=%s 里没有 fc2.hdf5/fc3.hdf5 —— 请指向 "
                     "fit-fc-thermal 的 step1_fit 产物目录（已拉回/已放到本集群）。" % _src)
        print("[..] S3_fc 导入模式：从 fit-fc-thermal 产物 %s 装配 step3_fc/（不重拟）" % _src)
        import reuse_fcfit
        _dim = str(conf.get("FC_IMPORT_DIM") or "").strip().lower() or None
        _mesh = str(conf.get("FC_IMPORT_MESH") or "").strip() or None
        _sc = str(conf.get("FC_IMPORT_SUPERCELL") or "").strip() or None
        try:
            reuse_fcfit.assemble(_src, out, poscar=None, supercell=_sc, dim=_dim,
                                 mesh=_mesh)
        except SystemExit:
            raise
        except Exception as _e:  # noqa: BLE001
            sys.exit("[ERROR] 从 fit-fc-thermal 产物装配 step3_fc 失败：%s\n"
                     "        （生成 phono3py_disp.yaml 需要 phono3py；本机 python=%s）"
                     % (_e, sys.executable))
        # 产物已就绪；写一个最小 submit.sh，autozt 收尾时提交的作业是 no-op。
        (out / "submit.sh").write_text(
            "#!/bin/bash\n"
            "# FC_IMPORT_DIR 导入模式：step3_fc/ 已在 gen 阶段装配完成，无需计算。\n"
            "#SBATCH --job-name=fcimport\n"
            "#SBATCH --output=queue.out\n"
            "#SBATCH --error=queue.err\n"
            "echo \"[S3_fc] imported from fit-fc-thermal; nothing to compute\"\n"
            "exit 0\n", encoding="utf-8", newline="\n")
        print("[DONE] step3_fc 已由 fit-fc-thermal 产物装配（导入模式），无需拟合作业")
        return

    # ---- 数据集默认指向本技能的 step2_disp_force（用户可 FIT_INPUT_DIR 覆盖）----
    if str(conf.get("FIT_INPUT_DIR") or "auto").strip().lower() == "auto":
        conf = g1.Overlay(conf, {"FIT_INPUT_DIR": SRC})

    # ---- NAC：MACE 给不出 Born，外部 DFPT 的 BORN 拷进数据集目录供 prep 拾取 ----
    _born = str(conf.get("NAC_BORN") or "").strip()
    if _born:
        _bp = Path(_born).expanduser()
        if not _bp.is_file():
            sys.exit("[ERROR] NAC_BORN 指向的文件不存在：%s" % _bp)
        _dst = cwd / SRC / "BORN"
        _dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(str(_bp), str(_dst))
        print("[OK] BORN <- %s（外部 DFPT 结果，MACE 自己给不出 Born 电荷）" % _bp)

    g1.main(conf=conf, outdir=OUTDIR, job_label="S3fit")


if __name__ == "__main__":
    main()
