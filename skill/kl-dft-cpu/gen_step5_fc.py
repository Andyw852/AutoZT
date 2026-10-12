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
    # 导入模式来源校验：bundle 的 dataset_dir 指向非本技能（kl-mlff/phonon-mlff 等）时，
    # 默认拒绝导入；确要跨来源复用设 true（数据来源校验仍保留结构/超胞一致性）。
    "FC_IMPORT_ALLOW_FOREIGN": (False, "bool"),
    # —— imag_policy 近 Γ 判据阈值（唯一真源 _common/imag_policy.py）——
    "IMAG_THR_STRICT": (0.10, "float"),
    "IMAG_QGAMMA": (0.05, "float"),
    "IMAG_QGAMMA_GRACE": (1.2, "float"),
    # —— 2D ZA 二次性（3-3b 再并入共享 fc_fit_driver，这里先保留键）——
    "ZA_CHECK": ("auto", "str"),
    "ZA_QMAX": (0.05, "float"),
})


def _bundle_structure_signature(path):
    """(lattice 3x3, symbols, wrapped scaled_positions Nx3)，解析不了返回 None。"""
    try:
        from phonopy.interface.vasp import read_vasp
        import numpy as np
        a = read_vasp(str(path))
        frac = np.asarray(a.scaled_positions, float)
        frac = frac - np.floor(frac)
        return np.asarray(a.cell, float), list(a.symbols), frac
    except Exception:
        return None


def _structures_match(sig_a, sig_b, tol=1e-3):
    """比较两份结构签名 -> (bool|None, note)。None = 解析不了，跳过但不拦。"""
    if sig_a is None or sig_b is None:
        return None, "结构无法解析（缺 phonopy 或 POSCAR 异常）"
    import numpy as np
    la, sa, fa = sig_a
    lb, sb, fb = sig_b
    if sa != sb:
        return False, "元素种类/顺序不一致：%s vs %s" % (sa, sb)
    if fa.shape != fb.shape:
        return False, "原子数不一致：%d vs %d" % (len(fa), len(fb))
    _scale = max(float(np.abs(la).max()), 1.0)
    if float(np.abs(la - lb).max()) > tol * _scale:
        return False, "晶格矩阵不一致"
    diff = np.abs(fa - fb)
    diff = np.minimum(diff, 1.0 - diff)          # 周期等价
    if float(diff.max()) > tol:
        return False, "分数坐标不一致（max |Δ| = %.4f）" % diff.max()
    return True, "结构一致"


def _yaml_supercell_matrix(path):
    try:
        import numpy as np
        import yaml
        y = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        m = y.get("supercell_matrix") if isinstance(y, dict) else None
        return None if m is None else np.asarray(m, float)
    except Exception:
        return None


def _bundle_provenance(src, cwd, allow_foreign):
    """kl_bundle 来源校验（结构 / 超胞 / 数据来源），通过则返回 provenance dict。"""
    import json
    import numpy as np
    prov = {"source": str(src), "allow_foreign": bool(allow_foreign)}
    # 1) 数据来源：fc_dataset.json.dataset_dir(_abs) 必须指向本技能 step4_disp
    fc_ds = {}
    if (src / "fc_dataset.json").is_file():
        try:
            fc_ds = json.loads((src / "fc_dataset.json").read_text(encoding="utf-8"))
        except Exception:
            fc_ds = {}
    ds_dir = str(fc_ds.get("dataset_dir_abs") or fc_ds.get("dataset_dir") or "")
    prov["dataset_dir"] = ds_dir
    is_kl = "kl-dft-cpu" in ds_dir
    if not is_kl and not allow_foreign:
        sys.exit("[ERROR] kl_bundle 的数据来源不是本技能的 step4_disp（dataset_dir=%s）。\n"
                 "        它可能来自 kl-mlff / phonon-mlff 等其它技能，力常数与当前 DFT 数据集无关。\n"
                 "        确要跨来源复用：设 FC_IMPORT_ALLOW_FOREIGN=true；否则核对数据来源后重推 bundle。"
                 % (ds_dir or "（缺 fc_dataset.json）"))
    prov["foreign_source"] = (not is_kl)
    # 2) 结构一致：bundle POSCAR vs step4_disp POSCAR
    sig_b = _bundle_structure_signature(src / "POSCAR")
    sig_k = _bundle_structure_signature(cwd / "step4_disp" / "POSCAR")
    ok, note = _structures_match(sig_b, sig_k)
    prov["structure_check"] = note
    if ok is False:
        sys.exit("[ERROR] kl_bundle 的结构与 step4_disp 不一致：%s。\n"
                 "        通常是 S1 换了结构 / bundle 来自另一材料，导入会得到错误力常数。\n"
                 "        处置：按当前结构重跑 fit-fc-thermal 后再推 bundle。" % note)
    if ok is None:
        print("[WARN] %s（不拦）" % note, flush=True)
    # 3) 超胞一致：bundle fc_dataset.json.supercell_matrix vs step4_disp phono3py_disp.yaml
    sc_b = fc_ds.get("supercell_matrix")
    sc_k = _yaml_supercell_matrix(cwd / "step4_disp" / "phono3py_disp.yaml")
    prov["supercell_matrix"] = sc_b
    if sc_b is not None and sc_k is not None:
        _match = bool(np.allclose(np.asarray(sc_b, float), sc_k, atol=1e-6))
        prov["supercell_check"] = "一致" if _match else "不一致"
        if not _match:
            sys.exit("[ERROR] kl_bundle 的超胞矩阵与 step4_disp 不一致（%s vs %s）。\n"
                     "        超胞不同，力常数与当前数据集不对应，拒绝导入。"
                     % (sc_b, sc_k.tolist()))
    else:
        prov["supercell_check"] = "跳过（缺 supercell_matrix 数据）"
    return prov


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
        # ---- 来源校验（结构 / 超胞 / 数据来源），通过才装配（第 4 条）----
        _prov = _bundle_provenance(_src, cwd, bool(conf.get("FC_IMPORT_ALLOW_FOREIGN")))
        print("[..] 导入来源校验通过：结构=%s；超胞=%s；来源=%s"
              % (_prov["structure_check"], _prov["supercell_check"],
                 _prov["dataset_dir"] or "（无 fc_dataset.json）"), flush=True)
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
        # ---- NAC：导入不重拟、bundle 不带 BORN；_clean_fit_artifacts 已把旧 step5_fc/BORN
        #   删掉（旧结构残留），这里按当前 step3_nac 补回（结构一致性前面已校验过）。
        #   3D 极性材料 NAC=auto 需要它，否则 κ 少 LO-TO 修正。
        _born = cwd / "step3_nac" / "BORN"
        if _born.is_file():
            shutil.copyfile(str(_born), str(out / "BORN"))
            print("[OK] BORN <- step3_nac/BORN（导入模式补回，NAC=auto 才会开 NAC）", flush=True)
        import json as _json
        (out / "import_provenance.json").write_text(
            _json.dumps(_prov, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
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
    _born_dst = disp / "BORN"
    if _born.is_file():
        shutil.copyfile(str(_born), str(_born_dst))
        print("[OK] BORN <- step3_nac/BORN（供共享引擎 prep 拾取）", flush=True)
    elif _born_dst.is_file():
        # NAC 已关但 step4_disp 里还留着上一轮的旧 BORN —— 移除，避免共享引擎误加 NAC
        _born_dst.unlink()
        print("[WARN] step3_nac/BORN 不存在（NAC 已关），移除 step4_disp/ 里的旧 BORN", flush=True)

    # ---- 委托共享引擎（数据集 = 本技能的 step4_disp）----
    if str(conf.get("FIT_INPUT_DIR") or "auto").strip().lower() == "auto":
        conf = g1.Overlay(conf, {"FIT_INPUT_DIR": DISP_DIR})
    g1.main(conf=conf, outdir=OUTDIR, job_label="S5fit")


if __name__ == "__main__":
    main()
