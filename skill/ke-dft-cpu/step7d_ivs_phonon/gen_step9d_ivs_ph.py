#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step9d_ivs_ph.py —— [0050] 谷间 IVS（K↔K′）冻声子第一步：√3×√3 超胞的 Γ̃ 点 DFPT（step7d_ivs_phonon，S7.3_ivsph）。

做什么：原胞的 K 点声子在 √3×√3 超胞里折叠成 Γ̃ 点的模（ivs_frozen.py 有推导）。本步建超胞、算超胞 Γ̃ 点 DFPT
（IBRION=8，只要力常数，不要介电），给 S7.4 挑区边界模、做 ±u 位移。
  · 超胞矩阵：六角面内晶格的 K、K′ 折叠到 Γ̃ 的最小超胞（ivs_frozen.fold_matrix，det = 3）。不是六角晶格 -> 报错
    （没有 K 谷；本方法只做 K↔K′）。
  · INCAR：照抄 step7_deform/undeformed 的单点设置（泛函、ENCUT、PREC 与 S7/S7.2 同源），改成 DFPT：
    IBRION=8、NSW=1、ISYM=2、NCORE=1、KPAR=1（与 S5 DFPT 模板同一组保守并行设置），不要 LEPSILON/LVHAR/LWAVE/LCHARG。
  · KPOINTS：Γ 中心，面内按 S7 网格的 k 间距换算到超胞（IVS_KMESH 可覆盖），kz 照抄 S7。
  · [0050] A 方案只支持非 SOC：S7 undeformed 或 S3/S3b 是 SOC（LSORBIT/LNONCOLLINEAR）-> 报错。
产出：step7d_ivs_phonon/{POSCAR（超胞）, POSCAR.prim, INCAR, KPOINTS, POTCAR, submit.sh, ivs_supercell.json}。
"""
import hashlib
import json
import math
import shutil
import sys
import time as _time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import ke_common as kc  # noqa: E402
import stepconf  # noqa: E402
from dim_common import require_dim, resolve_tpl  # noqa: E402

iv = None                     # ivs_frozen：main() 里切好 python 环境后再 import（它要 numpy）

# =========================== 可改参数区（step.conf [params] 可覆盖） ===========================
OUTDIR_NAME = "step7d_ivs_phonon"
PREV_CANDS = ["step1_opt", "step1_std_opt"]
DEFORM_DIR = "step7_deform"
STEP_LABEL = "S7.3_ivsph"
MANIFEST = "ivs_supercell.json"
_SKILL_REV = "2026-10-14-0050"   # [0050] 谷间 IVS（K↔K′ 零阶）冻声子：超胞 DFPT
# auto = 按 S7 undeformed 的面内 k 间距换算到超胞；或 "n1 n2 n3"（Γ 中心）
IVS_KMESH = "auto"
SPEC = {"IVS_KMESH": (None, "str")}
DFPT_TAGS = {"IBRION": "8", "NSW": "1", "ISYM": "2", "NCORE": "1", "KPAR": "1", "LEPSILON": None, "LVHAR": None,
             "LWAVE": ".FALSE.", "LCHARG": ".FALSE.", "LORBIT": None, "NBANDS": None, "ISTART": "0", "ICHARG": "2"}
# =================================================================================================


def apply_conf(cwd):
    global IVS_KMESH
    got = {}
    if (Path(cwd) / stepconf.CONF_NAME).is_file():
        p = stepconf.load(SPEC, OUTDIR_NAME, str(cwd), strict=False)
        for k in SPEC:
            if p[k] is not None:
                got[k] = p[k]
                globals()[k] = p[k]
    return got


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def k_star_hexagonal(lat):
    """面内六角晶格（|a1| = |a2|、夹角 60°/120°）的 K、K′（原胞倒格分数坐标）；不是六角 -> None。
    取 i/3、j/3 的点里到 Γ 最近像距离 = |b|/√3（布里渊区角）的两类，各取一个代表。"""
    import numpy as np
    a = np.asarray(lat, float)
    la, lb = np.linalg.norm(a[0]), np.linalg.norm(a[1])
    ang = math.degrees(math.acos(float(a[0] @ a[1]) / (la * lb)))
    if abs(la - lb) > 1e-3 * la or min(abs(ang - 60.0), abs(ang - 120.0)) > 0.5:
        return None
    rec = 2 * np.pi * np.linalg.inv(a).T
    corner = np.linalg.norm(rec[0]) / math.sqrt(3.0)
    pts = []
    for i in range(3):
        for j in range(3):
            if (i, j) == (0, 0):
                continue
            k = np.array([i / 3.0, j / 3.0, 0.0])
            best = min(np.linalg.norm((k + np.array([p, q, 0])) @ rec) for p in (-1, 0, 1) for q in (-1, 0, 1))
            if abs(best - corner) < 1e-6 * corner:
                pts.append(k)
    if len(pts) < 2:
        return None
    K = pts[0]
    Kp = (-K) % 1.0
    if not any(np.allclose(Kp, p) for p in pts):
        return None
    return [K.tolist(), Kp.tolist()]


def supercell_kmesh(lat_prim, mesh_prim, lat_sc):
    """S7 网格的面内 k 间距（取两方向较密者）换算到超胞，kz 照抄。"""
    import numpy as np
    rp = 2 * np.pi * np.linalg.inv(np.asarray(lat_prim, float)).T
    rs = 2 * np.pi * np.linalg.inv(np.asarray(lat_sc, float)).T
    s = min(np.linalg.norm(rp[i]) / max(int(mesh_prim[i]), 1) for i in range(2))
    return [max(1, int(math.ceil(np.linalg.norm(rs[i]) / s - 1e-6))) for i in range(2)] + [int(mesh_prim[2])]


def write_kpoints(path, mesh, comment):
    Path(path).write_text("%s\n0\nGamma\n%d %d %d\n0 0 0\n" % (comment, mesh[0], mesh[1], mesh[2]),
                          encoding="utf-8", newline="\n")


def main():
    import numpy as np
    global iv
    cwd = Path.cwd()
    kc.reexec_with_pymatgen(__file__)
    import ivs_frozen as iv                    # noqa: F811
    from pymatgen.core import Structure
    from pymatgen.io.vasp import Poscar
    out = cwd / OUTDIR_NAME
    out.mkdir(exist_ok=True)
    overrides = apply_conf(cwd)
    prev = kc.find_prev_dir(cwd, PREV_CANDS)
    if prev is None:
        sys.exit("[ERROR] 找不到含 CONTCAR 的上一步：%s" % PREV_CANDS)
    kc.relay_poscar(prev / "CONTCAR", out / "POSCAR.prim", "step1_opt")
    dim = kc.read_method_dim(prev / kc.METHOD_FILE) or kc.resolve_dim_for(out / "POSCAR.prim", "auto")[0]
    require_dim(dim, ("2d",), OUTDIR_NAME, why="谷间 IVS（K↔K′）冻声子只给 S8.4（二维）用，K 谷要六角面内晶格")
    und = cwd / DEFORM_DIR / "undeformed"
    for n in ("INCAR", "KPOINTS", "POTCAR", "POSCAR"):
        if not (und / n).is_file():
            sys.exit("[ERROR] 缺 %s/undeformed/%s：先让 S7_deform 生成输入（本步照抄它的泛函/ENCUT/POTCAR 与 k 间距）"
                     % (DEFORM_DIR, n))
    if kc.incar_is_ncl(und / "INCAR") or kc.detect_ncl(cwd):
        sys.exit("[ERROR] 本材料是 SOC 能带（LSORBIT/LNONCOLLINEAR）。[0050] 谷间 IVS 冻声子只做 A 方案（非 SOC）：SOC 下 K/K′ "
                 "二重态被自旋劈裂、谷间散射只能落到同自旋子带（Xiao 等 PRL 108, 196802），劈裂法与 dp_optical 的星平均"
                 "态密度都不适用。本材料 tf_*.yaml 写 ivs_frozen: false 关掉本组")
    base = Structure.from_file(str(out / "POSCAR.prim"))
    s7 = Structure.from_file(str(und / "POSCAR"))
    if [str(x.specie) for x in base] != [str(x.specie) for x in s7] or \
            float(np.abs(np.asarray(base.lattice.matrix) - np.asarray(s7.lattice.matrix)).max()) > 1e-4:
        sys.exit("[ERROR] step1_opt CONTCAR 与 %s/undeformed/POSCAR 不同：S7 是旧结构算的，先 retry S7" % DEFORM_DIR)
    lat = np.asarray(base.lattice.matrix, float)
    members = k_star_hexagonal(lat)
    if members is None:
        sys.exit("[ERROR] 面内晶格不是六角（|a|≠|b| 或夹角不是 60°/120°）：没有 K/K′ 谷，本方法（K↔K′ 零阶 IVS）不适用。"
                 "本材料 tf_*.yaml 写 ivs_frozen: false 关掉本组")
    M = iv.fold_matrix(members, lat[:2, :2])
    species = [str(s.specie) for s in base]
    sc = iv.build_supercell(lat, base.frac_coords, species, M)
    st = Structure(sc["lattice"], sc["species"], sc["frac"])
    snap = kc.snapshot_inputs(out, ("INCAR", "POSCAR", "KPOINTS", "POTCAR"))
    Poscar(st).write_file(str(out / "POSCAR"), significant_figures=12)
    mesh_p = kc.read_kpoints_mesh(und / "KPOINTS")
    if mesh_p is None or kc.read_kpoints_gamma(und / "KPOINTS") is not True:
        sys.exit("[ERROR] %s/undeformed/KPOINTS 不是 Γ 中心的自动网格：本步按它的 k 间距换算超胞网格" % DEFORM_DIR)
    spec = str(IVS_KMESH or "auto").strip().lower()
    if spec == "auto":
        mesh = supercell_kmesh(lat, mesh_p, sc["lattice"])
    else:
        try:
            mesh = [int(x) for x in spec.split()]
            assert len(mesh) == 3 and min(mesh) >= 1
        except (ValueError, AssertionError):
            sys.exit("[ERROR] IVS_KMESH=%r 应为 auto 或三个正整数" % IVS_KMESH)
    write_kpoints(out / "KPOINTS", mesh, "IVS supercell Gamma-centred (S7 spacing)")
    shutil.copyfile(und / "POTCAR", out / "POTCAR")
    d = kc.merge_incar(kc.parse_incar((und / "INCAR").read_text(encoding="utf-8", errors="ignore")), DFPT_TAGS)
    (out / "INCAR").write_text(kc.incar_text(d, system="%s ivs supercell DFPT" % cwd.name), encoding="utf-8",
                               newline="\n")
    _ic_log = []
    stepconf.apply_incar_file(out / "INCAR", log=_ic_log)
    for _m in _ic_log:
        print("[..] %s" % _m)
    stamp = _time.strftime("%Y%m%d%H%M%S")
    changed = kc.inputs_changed(out, snap)
    if changed:
        n = kc.archive_stale_outputs(out, "stale-input-" + stamp)
        if n:
            print("[..] 输入变了（%s）-> 归档旧产物 %d 个，会重算" % ("/".join(changed), n))
    else:
        kc.archive_if_run_mismatch(out, label=OUTDIR_NAME)
    submit_body = resolve_tpl(Path(__file__).resolve().parent, "submit_std", dim).read_text(encoding="utf-8")
    (out / "submit.sh").write_text(submit_body.replace("{{JOBNAME}}", "%s-ke-dft-cpu-%s" % (cwd.name, STEP_LABEL)),
                                   encoding="utf-8", newline="\n")
    stepconf.apply_submit(out / "submit.sh", stepconf.read_submit(stepconf.CONF_NAME, used_incar=True))
    folded = iv.folded_kpoints(M)
    names = ["k%d" % i for i in range(len(folded))]
    star = [names[i] for i, k in enumerate(folded)
            if any(np.allclose(np.asarray(k[:2]) % 1.0, np.asarray(m[:2]) % 1.0, atol=1e-6) for m in members)]
    manifest = {
        "schema": "ivs_supercell/1", "skill_rev": _SKILL_REV, "created": _time.strftime("%Y-%m-%d %H:%M:%S"),
        "M": M.tolist(), "det": int(sc["det"]), "members_frac": members,
        "folded": {n: np.round(k, 8).tolist() for n, k in zip(names, folded)}, "star": star,
        "prim": {"lattice": lat.tolist(), "frac": np.asarray(base.frac_coords).tolist(), "species": species},
        "supercell": {"lattice": np.asarray(sc["lattice"]).tolist(), "species": sc["species"],
                      "frac": np.asarray(sc["frac"]).tolist(),
                      "origin": [[int(i), [int(t[0]), int(t[1])]] for i, t in sc["origin"]],
                      "translations": [[int(t[0]), int(t[1])] for t in sc["translations"]]},
        "kmesh_sc": mesh, "kmesh_s7": mesh_p, "overrides": overrides,
        "inputs_from": {"dir": "%s/undeformed" % DEFORM_DIR,
                        "sha256": {n: _sha(und / n) for n in ("INCAR", "KPOINTS", "POTCAR", "POSCAR")}},
    }
    (out / MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    print("[DONE] %s：超胞 M=%s（%d 个原胞、%d 个原子），K/K′（%s）折叠到 Γ̃；DFPT 网格 %s（S7 %s）"
          % (OUTDIR_NAME, M.tolist(), sc["det"], len(sc["species"]), "/".join(star), mesh, mesh_p))


if __name__ == "__main__":
    main()
