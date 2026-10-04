#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""reuse_fcfit.py —— 把 fit-fc-thermal 的 step1_fit 产物装配成 kl-mlff 的 step3_fc/ 布局。

fit-fc-thermal（S1_fit）与 kl-mlff（S3_fc）拟合出的力常数内容相同，只是目录布局不同：
本脚本做纯搬运 + 生成 kl-mlff 侧需要的 phono3py_disp.yaml / klmlff_params.txt，
让 kl-mlff 的 S4_kappa 直接吃 fit-fc-thermal 的拟合结果，不必在 mlff 侧再重拟。

映射（fit-fc-thermal step1_fit/  ->  kl-mlff step3_fc/）：
    fc2.hdf5, fc3.hdf5            ->  fc2.hdf5, fc3.hdf5（顶层，mlff 布局）
    phono3py_disp.yaml (若有)      ->  phono3py_disp.yaml
    phono3py_params.yaml (若有)    ->  phono3py_params.yaml
    都没有则用 POSCAR + supercell_matrix 现生成 phono3py_disp.yaml（需要 phono3py）
    shengbte/FORCE_CONSTANTS_*     ->  shengbte/FORCE_CONSTANTS_*
    phonon_summary.json            ->  phonon_summary.json（虚频闸 marker）
    fc_fit_summary.json            ->  fc_fit_summary.json（溯源）
    POSCAR                         ->  POSCAR
并写 klmlff_params.txt：DIM / SUPERCELL / FC2_SUPERCELL / MESH / METHOD
和 workflow_method.txt（溯源）。

用法（在 mlff 集群上、有 phono3py 的环境里跑；gen_step3_fc 在 FC_IMPORT_DIR 非空时
会自动调用同一个 assemble()）：
    python reuse_fcfit.py --fcfit-dir <fit-fc-thermal step1_fit> --out <step3_fc> \
        [--poscar POSCAR] [--supercell "3 3 1"] [--dim 3d] [--mesh "24 24 24"]
"""
import argparse
import json
import shutil
import sys
from pathlib import Path


def _supercell_from_fcdataset(src):
    """fc_dataset.json 里的 supercell_matrix（3x3），拿不到返回 None。"""
    p = Path(src) / "fc_dataset.json"
    if not p.is_file():
        return None
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None
    scm = d.get("supercell_matrix")
    if not scm:
        return None
    try:
        import numpy as np
        m = np.asarray(scm, float)
        if m.shape == (3, 3):
            return m
    except Exception:
        return None
    return None


def _parse_supercell(text):
    import numpy as np
    vals = [float(x) for x in str(text).replace(",", " ").split()]
    if not vals:
        return None
    a = np.asarray(vals, float)
    if a.size == 3:                     # "3 3 1" -> 对角矩阵
        return np.diag(a)
    return a.reshape(3, 3)


def _supercell_str(m):
    # klmlff_params.txt 的 SUPERCELL 认对角三整数（gen_step4 直接 split）。
    import numpy as np
    m = np.asarray(m, float)
    if m.shape == (3, 3):
        d = np.diag(m)
        if np.abs(m - np.diag(d)).max() < 1e-8:
            return " ".join(str(int(round(float(x)))) for x in d)
    return " ".join(str(int(round(float(x)))) for x in m.reshape(-1))


def _default_mesh(dim):
    return "24 24 1" if str(dim).lower() == "2d" else "24 24 24"


def _make_phono3py_disp(poscar, scm, out_yaml):
    """用原胞 POSCAR + 超胞矩阵生成 phono3py_disp.yaml（需要 phono3py）。"""
    from phonopy.interface.vasp import read_vasp
    from phono3py import Phono3py
    uc = read_vasp(str(poscar))
    ph3 = Phono3py(uc, supercell_matrix=scm.tolist())
    ph3.save(str(out_yaml))
    return len(uc), len(ph3.supercell)


def _copy_tree_files(srcdir, dstdir):
    srcdir, dstdir = Path(srcdir), Path(dstdir)
    if not srcdir.is_dir():
        return 0
    dstdir.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in srcdir.iterdir():
        if f.is_file():
            shutil.copyfile(f, dstdir / f.name)
            n += 1
    return n


def assemble(fcfit_dir, out, poscar=None, supercell=None, dim=None,
             mesh=None, verbose=True):
    """把 fit-fc-thermal 的 step1_fit 目录装配成 kl-mlff 的 step3_fc 目录（幂等）。"""
    src = Path(fcfit_dir).expanduser()
    out = Path(out).expanduser()
    if not src.is_absolute():
        src = (Path.cwd() / src).resolve()
    else:
        src = src.resolve()
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    missing = [f for f in ("fc2.hdf5", "fc3.hdf5") if not (src / f).is_file()]
    if missing:
        sys.exit("[ERROR] %s 缺 %s —— 不是完整的 fit-fc-thermal step1_fit 产物"
                 % (src, "、".join(missing)))

    # ---- 力常数：mlff step3_fc 顶层布局（与 kl-dft-cpu 的 phono3py/ 子目录不同）----
    for f in ("fc2.hdf5", "fc3.hdf5"):
        shutil.copyfile(src / f, out / f)

    # ---- 力常数：ShengBTE 格式 ----
    sb = src / "shengbte"
    if sb.is_dir():
        _copy_tree_files(sb, out / "shengbte")

    # ---- POSCAR ----
    pos = Path(poscar) if poscar else (src / "POSCAR")
    if not pos.is_file():
        pos = src / "POSCAR"
    if pos.is_file():
        shutil.copyfile(pos, out / "POSCAR")

    # ---- phono3py YAML（最好从源带过来，否则现生成）----
    scm = _parse_supercell(supercell) if supercell else None
    if scm is None:
        scm = _supercell_from_fcdataset(src)
    for _y in ("phono3py_params.yaml", "phono3py_disp.yaml"):
        if (src / _y).is_file():
            shutil.copyfile(src / _y, out / _y)
            break
    else:
        if not (out / "phono3py_params.yaml").is_file()                 and not (out / "phono3py_disp.yaml").is_file():
            if scm is None:
                sys.exit("[ERROR] 源里没有 phono3py YAML，fc_dataset.json 也没有 "
                         "supercell_matrix —— 请用 --supercell 'n n n' 指定超胞")
            n_uc, n_sc = _make_phono3py_disp(out / "POSCAR", scm,
                                             out / "phono3py_disp.yaml")
            if verbose:
                print("[OK] 生成 phono3py_disp.yaml（%d 原胞原子 -> %d 超胞原子）"
                      % (n_uc, n_sc))

    # ---- marker / 溯源记录 ----
    for f in ("phonon_summary.json", "fc_fit_summary.json", "fc_dataset.json"):
        if (src / f).is_file():
            shutil.copyfile(src / f, out / f)
    if not (out / "phonon_summary.json").is_file():
        sys.exit("[ERROR] %s 缺 phonon_summary.json（S3_fc 虚频闸 marker）" % src)

    # ---- klmlff_params.txt ----
    if dim is None:
        dim = "3d"
        try:
            _ps = json.loads((out / "phonon_summary.json").read_text(encoding="utf-8"))
            if _ps.get("is_2d"):
                dim = "2d"
        except Exception:
            pass
    if supercell:
        sc_str = _supercell_str(_parse_supercell(supercell))
    elif scm is not None:
        sc_str = _supercell_str(scm)
    else:
        sc_str = ""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import klmlff_common as klc
    klc.write_kl_params(out / klc.KL_PARAMS, DIM=str(dim).lower(),
                        SUPERCELL=sc_str, FC2_SUPERCELL=sc_str,
                        MESH=(mesh or _default_mesh(dim)), METHOD="import")
    # workflow_method.txt：仅溯源（S4 只拷贝不读）。MODEL 用 import 标记，避免误当
    # MACE 弛豫产物；DIM 供其它步骤继承。
    klc.write_method(out / klc.METHOD_FILE, FUNC="pbesol", DIM=str(dim).lower(),
                     MODEL="fit-fc-thermal-import")
    if verbose:
        print("[OK] 装配完成：%s" % out)
        print("     fc2.hdf5 / fc3.hdf5 / phono3py YAML / shengbte/ / "
              "phonon_summary.json / klmlff_params.txt")
    return out


def main():
    ap = argparse.ArgumentParser(description="fit-fc-thermal -> kl-mlff step3_fc 布局适配")
    ap.add_argument("--fcfit-dir", required=True,
                    help="fit-fc-thermal 的 step1_fit 目录（含 fc2/fc3.hdf5）")
    ap.add_argument("--out", required=True, help="kl-mlff 的 step3_fc 目录")
    ap.add_argument("--poscar", default=None)
    ap.add_argument("--supercell", default=None, help='超胞矩阵，如 "3 3 1"')
    ap.add_argument("--dim", default=None, help="2d | 3d（默认按 phonon_summary.json 的 is_2d 判）")
    ap.add_argument("--mesh", default=None, help='BTE 网格，如 "24 24 24"')
    a = ap.parse_args()
    assemble(a.fcfit_dir, a.out, a.poscar, a.supercell, a.dim, a.mesh)


if __name__ == "__main__":
    main()
