#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step5_fc.py —— kl-dft-cpu S5_fc：拟合力常数（共享 fit-fc-thermal 引擎）。

保留 DFT 专属的【帧一致性 + SCF 收敛门禁】两个预检，然后把 step4_disp 的
位移+力数据集交给共享引擎 _common/fcfit/gen_step1_fit，写 fit_config.json +
submit.sh；计算节点跑 fc_fit_driver.py prep|fit|post。

导入模式（跨技能复用）：检测到推送来的 kl_bundle/（或 FC_IMPORT_DIR）→ 装配
step5_fc 布局、跳过拟合（见 reuse_fcfit.py）。

注意：换共享引擎后，本技能旧引擎的 PHEASY_ENABLE_FC=4（四声子）与 2D ZA 二次性
尚未并入共享 fc_fit_driver，后续单独补（3-3b / 四声子）。
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
import kl_common as kc      # noqa: E402  (帧一致性/SCF 门禁)

OUTDIR = "step5_fc"
STEP = "step5_fc"
DISP_DIR = "step4_disp"

SPEC = dict(g1.SPEC)
SPEC.update({
    "FUNC": ("pbesol", "str"),
    # —— 由 fit-fc-thermal 产物导入（跨技能复用 fc2/fc3）——
    "FC_IMPORT_DIR": ("", "str"),
    "FC_IMPORT_SUPERCELL": ("", "str"),
    "FC_IMPORT_DIM": ("", "str"),
    "FC_IMPORT_MESH": ("", "str"),
    "FC_IMPORT_FUNCTIONAL": ("", "str"),
    # —— imag_policy 近 Γ 判据阈值（唯一真源 _common/imag_policy.py）——
    "IMAG_THR_STRICT": (0.10, "float"),
    "IMAG_QGAMMA": (0.05, "float"),
    "IMAG_QGAMMA_GRACE": (1.2, "float"),
    # —— 2D ZA 二次性（3-3b 再并入共享 fc_fit_driver，这里先保留键）——
    "ZA_CHECK": ("auto", "str"),
    "ZA_QMAX": (0.05, "float"),
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
        print("[..] S5_fc 导入模式：从 fit-fc-thermal 产物 %s 装配 step5_fc/（不重拟）" % _src)
        import reuse_fcfit
        _dim = str(conf.get("FC_IMPORT_DIM") or "").strip().lower() or None
        _mesh = str(conf.get("FC_IMPORT_MESH") or "").strip() or None
        _sc = str(conf.get("FC_IMPORT_SUPERCELL") or "").strip() or None
        _func = str(conf.get("FC_IMPORT_FUNCTIONAL") or conf.get("FUNC")
                    or "pbesol").strip()
        try:
            reuse_fcfit.assemble(_src, out, poscar=None, supercell=_sc, dim=_dim,
                                 mesh=_mesh, nac=False, functional=_func)
        except SystemExit:
            raise
        except Exception as _e:  # noqa: BLE001
            sys.exit("[ERROR] 从 fit-fc-thermal 产物装配 step5_fc 失败：%s\n"
                     "        （生成 phono3py_disp.yaml 需要 phono3py；本机 python=%s）"
                     % (_e, sys.executable))
        (out / "submit.sh").write_text(
            "#!/bin/bash\n"
            "# FC_IMPORT_DIR 导入模式：step5_fc/ 已在 gen 阶段装配完成，无需计算。\n"
            "#SBATCH --job-name=fcimport\n"
            "#SBATCH --output=queue.out\n"
            "#SBATCH --error=queue.err\n"
            "echo \"[S5_fc] imported from fit-fc-thermal; nothing to compute\"\n"
            "exit 0\n", encoding="utf-8", newline="\n")
        print("[DONE] step5_fc 已由 fit-fc-thermal 产物装配（导入模式），无需拟合作业")
        return

    # ---- 校验 step4 产物 ----
    disp = cwd / DISP_DIR
    if not (disp / "phono3py_disp.yaml").is_file():
        sys.exit("[ERROR] %s 缺 phono3py_disp.yaml（step4 未生成位移）" % disp)
    if not list(disp.glob("disp-*/vasprun.xml")):
        sys.exit("[ERROR] %s 下无 disp-*/vasprun.xml，位移单点还没算完" % disp)

    # ---- 帧一致性（力必须与位移对应）----
    _ok, _note = kc.check_frames_match_displacements(disp)
    print("[%s] S4 帧一致性：%s" % ("OK" if _ok else "FAIL", _note))
    if not _ok:
        sys.exit("[ERROR] %s\n        请清空 step4_disp 的 disp-*/POSCAR-*/phono3py_disp.yaml/SPOSCAR "
                 "后重跑 S4（或 -j S4_disp rerun）。" % _note)

    # ---- SCF 收敛门禁（NELM 截断/未收敛的帧不能拟合）----
    import json as _json
    _scf_ok, _scf = kc.check_outcar_scf_convergence(disp)
    print(kc.format_scf_report(_scf))
    (out / "scf_steps.json").write_text(
        _json.dumps(_scf, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    if not _scf_ok:
        _bad = "、".join(b["frame"] for b in _scf["bad_frames"][:20])
        sys.exit("[ERROR] S4 有 %d/%d 帧 SCF 未正常收敛（NELM 截断 / 无 aborting loop），"
                 "力不可信，拒绝拟合：\n        作废帧：%s\n"
                 "        处理：对作废帧调 INCAR 后 autozt -tt kl-dft-cpu -p <材料> "
                 "-j S4_disp retry 只补这些帧。"
                 % (_scf["n_bad"], _scf["n_frames"], _bad))

    # ---- NAC：共享引擎的 prep 从【数据集目录】读 BORN；本技能的 BORN 在 step3_nac ----
    _born = cwd / "step3_nac" / "BORN"
    if _born.is_file():
        shutil.copyfile(str(_born), str(disp / "BORN"))
        print("[OK] BORN <- step3_nac/BORN（供共享引擎 prep 拾取）", flush=True)

    # ---- 委托共享引擎（数据集 = 本技能的 step4_disp）----
    if str(conf.get("FIT_INPUT_DIR") or "auto").strip().lower() == "auto":
        conf = g1.Overlay(conf, {"FIT_INPUT_DIR": DISP_DIR})
    g1.main(conf=conf, outdir=OUTDIR, job_label="S5fit")


if __name__ == "__main__":
    main()
