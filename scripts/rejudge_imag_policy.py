#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rejudge_imag_policy.py —— 只读重判：对已有 fc2 按新 imag_policy 重算虚频结论。

用途（2026-09-22）：虚频判据统一到 skill/_common/imag_policy.py 后，对**已有材料**的
S5 产物按新规则重判，列出结论发生变化的材料。**本脚本只读**：不改任何文件、不提交作业。

依赖 phono3py + phonopy + numpy（在集群的 atomate2_p_a 环境里跑）；imag_policy.py 与
za_2d.py 需在 PYTHONPATH 或同目录。

用法：
    python3 rejudge_imag_policy.py --root /public/home/wangchao/Fullerene_Network/work \
        [--pattern '**/kl-dft-cpu/step5_fc/phono3py/fc2.hdf5'] [--json out.json]
"""
import argparse
import glob
import json
import os
import sys


def _load_cfg(step_dir):
    """读该材料 S5 的 step.conf 里的 imag 阈值（未写就用 imag_policy 默认值）。"""
    conf = {}
    p = os.path.join(step_dir, "step.conf")
    if not os.path.isfile(p):
        return conf
    try:
        for ln in open(p, encoding="utf-8", errors="replace").read().splitlines():
            ln = ln.split("#", 1)[0].strip()
            if "=" not in ln:
                continue
            k, v = [x.strip() for x in ln.split("=", 1)]
            if k in ("IMAG_THR", "IMAG_THR_STRICT", "IMAG_QGAMMA"):
                try:
                    conf[k] = float(v.strip().strip('"').strip("'"))
                except ValueError:
                    pass
    except OSError:
        pass
    return conf


def _samples(ph, is2d, vaxp, mesh_n):
    """返回 (freqs, qpts)：mesh + band 两套都塞进去（与 S5 的口径一致）。"""
    freqs, qpts = [], []
    try:
        if is2d:
            mesh = [mesh_n, mesh_n, mesh_n]
            mesh[int(vaxp)] = 1
            ph.run_mesh(mesh=mesh, with_eigenvectors=False, is_mesh_symmetry=True,
                        is_gamma_center=True)
        else:
            ph.run_mesh(mesh=float(mesh_n), with_eigenvectors=False, is_mesh_symmetry=True)
        md = ph.get_mesh_dict()
        freqs += [[float(x) for x in row] for row in md["frequencies"]]
        qpts += [[float(x) for x in q] for q in md["qpoints"]]
    except Exception as e:                         # noqa: BLE001
        print("  [WARN] mesh 失败：%s" % e)
    try:
        import kl_common as kc
        from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
        # 与 S5 的口径一致用 101 点/段（31 点会漏掉很局部的软模，2026-09-22 实测 P1_Mo2S3）
        paths, labels, _lat = kc.band_path_2d(ph.primitive.cell, 101, int(vaxp)) if is2d \
            else (None, None, None)
        if is2d:
            bands, conn = get_band_qpoints_and_path_connections(paths, npoints=101)
            ph.run_band_structure(bands, path_connections=conn, labels=labels)
            bd = ph.get_band_structure_dict()
            freqs += [[float(x) for x in f]
                      for seg in (bd.get("frequencies") or []) for f in seg]
            qpts += [[float(x) for x in q]
                     for seg in (bd.get("qpoints") or []) for q in seg]
    except Exception as e:                         # noqa: BLE001
        print("  [WARN] band 失败（只按 mesh 判）：%s" % e)
    return freqs, qpts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/public/home/wangchao/Fullerene_Network/work")
    ap.add_argument("--pattern", default="**/kl-dft-cpu/step5_fc/phono3py/fc2.hdf5")
    ap.add_argument("--mesh-n", type=int, default=60)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    import numpy as np
    import phono3py
    from phono3py.file_IO import read_fc2_from_hdf5
    from phonopy import Phonopy
    import imag_policy
    import za_2d

    hits = sorted(glob.glob(os.path.join(args.root, args.pattern), recursive=True))
    print("扫描到 %d 个 fc2.hdf5（root=%s）" % (len(hits), args.root))
    rows = []
    for fc2_path in hits:
        p3 = os.path.dirname(fc2_path)
        step_dir = os.path.dirname(p3)                 # .../step5_fc
        mat_dir = os.path.dirname(os.path.dirname(step_dir))
        name = os.path.basename(mat_dir)
        old = {}
        osp = os.path.join(step_dir, "phonon_summary.json")
        if os.path.isfile(osp):
            try:
                old = json.load(open(osp, encoding="utf-8"))
            except Exception:                      # noqa: BLE001
                old = {}
        try:
            ph3 = phono3py.load(os.path.join(p3, "phono3py_disp.yaml"),
                                produce_fc=False, is_nac=False, log_level=0)
            ph = Phonopy(ph3.unitcell, supercell_matrix=ph3.supercell_matrix,
                         primitive_matrix=ph3.primitive_matrix)
            ph.force_constants = np.asarray(read_fc2_from_hdf5(filename=fc2_path))
            ph.nac_params = None
            norms = [float(np.linalg.norm(v)) for v in ph.primitive.cell]
            is2d = max(norms) > 10.0
            vaxp = za_2d.vacuum_axis_in_primitive(ph, None) if is2d else None
            freqs, qpts = _samples(ph, is2d, vaxp, args.mesh_n)
            conf = _load_cfg(step_dir)
            r = imag_policy.classify_imag(freqs, qpts, is2d, vaxp, conf or None)
        except Exception as e:                     # noqa: BLE001
            rows.append({"material": name, "path": fc2_path, "error": str(e)})
            print("!! %-46s 重判失败：%s" % (name, e))
            continue
        old_stable = old.get("stable")
        old_min = old.get("min_frequency_THz", old.get("min_freq_THz"))
        new_fail = (r["verdict"] == "fail")
        changed = (old_stable is not None and bool(old_stable) == new_fail)
        rows.append({"material": name, "fc2": fc2_path, "is_2d": is2d,
                     "old_stable": old_stable, "old_min_freq_THz": old_min,
                     "new_verdict": r["verdict"], "new_imag_class": r["imag_class"],
                     "new_min_freq_THz": r["min_freq_THz"],
                     "conf_used": conf or "defaults", "changed": bool(changed)})
        print("%-46s old(stable=%s,min=%s)  ->  new(%s/%s,min=%s)%s"
              % (name, old_stable, old_min, r["verdict"], r["imag_class"],
                 None if r["min_freq_THz"] is None else round(r["min_freq_THz"], 4),
                 "   ← 结论变化" if changed else ""))

    chg = [x for x in rows if x.get("changed")]
    print("\n=== 结论发生变化的材料（%d / %d）===" % (len(chg), len(rows)))
    for x in chg:
        print("  %-46s old_stable=%s  new=%s/%s  min=%s"
              % (x["material"], x["old_stable"], x["new_verdict"], x["new_imag_class"],
                 x["new_min_freq_THz"]))
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(rows, fh, ensure_ascii=False, indent=2)
        print("\n已写 %s" % args.json)


if __name__ == "__main__":
    main()
