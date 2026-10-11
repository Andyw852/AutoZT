#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_dp_optical_e2e.py —— [0033] ODP/IVS 与 ADP 谷内掩码的真 AMSET 端到端自检（需 amset + pymatgen 环境）。

用法（在 skill/ke-dft-cpu/step8.4_amset2d 目录里，约 5–15 分钟）：
    python test_dp_optical_e2e.py              # 默认 nworkers=2（顺带检验 spawn 子进程里的掩码）
    python test_dp_optical_e2e.py --nworkers 4
退出码 0 = 全部 PASS。

合成模型：六方单原子 2D 晶格（a=3.18 Å，c=20 Å），导带/价带都是三角格子紧束缚带，极值在 K、K′（两谷，g=2）。
电子开 ODP（K>K）+ IVS（K>K′），空穴只开 ODP；电子开 ADP 谷内掩码。同一模型跑两次（掩码开/关），检查：
  1. AMSET 落进 mesh 的 ODP/IVS 率 = dp_optical.carrier_rates 直接算的（窗口内外每个 k 点，相对差 < 1e-9）；
  2. 掩码开/关的 ADP 率之比：电子 = 0.5（两等价谷、unity overlap、q 无关的核 -> 谷间恰好占一半），空穴 = 1；
  3. 电子单机制迁移率 IVS/ODP 落在 [0.9, 1.6]（由两组参数直接估算 ≈ 1.2；修 FD 窗口外填充之前是 ~40）；
  4. 日志里有横幅、谷表核对（谷底离中心 < 0.01|b1|）、"FD 窗口外保留解析率"；
  5. 作业里的率 = 按 AMSET 不可约映射展开的率（[0034-r2]：不一致会让 postprocess 的复现检查失败、整步 FAIL）。
     另跑一个矩形胞两星模型（只开 ODP/IVS，几分钟）专测这一项：单星模型里谷归属不影响率，这一项在那里不可能失败；
     矩形两星模型在修复前有 1044 处不一致。
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
K = [[1 / 3, 1 / 3, 0.0], [-1 / 3, -1 / 3, 0.0]]
FAILS = []


def ok(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def record(mask):
    return {"cell_c_A": 20.0, "layer_normal": [0, 0, 1], "kz_cap_rmax": 1,
            "eps_inf_slab": [[3, 0, 0], [0, 3, 0], [0, 0, 1.2]],
            "eps_static_slab": [[3.2, 0, 0], [0, 3.2, 0], [0, 0, 1.2]],
            "mechanisms_2d": ["ADP", "ODP", "IVS"],
            "dp_optical": {
                "electron": {"valleys": [{"label": "K", "g": 2, "k_frac": K}],
                             "channels": [{"from": "K", "to": "K", "D_eV_A": 5.0, "hw_meV": 50.0, "kind": "intra"},
                                          {"from": "K", "to": "K", "D_eV_A": 3.0, "hw_meV": 40.0,
                                           "kind": "inter_star"}]},
                "hole": {"valleys": [{"label": "K", "g": 2, "k_frac": K}],
                         "channels": [{"from": "K", "to": "K", "D_eV_A": 4.0, "hw_meV": 48.0, "kind": "intra"}]},
                "adp_intravalley": {"electron": mask, "hole": False}}}


def record_rect():
    """[0034-r2] 矩形胞两星模型（电子 Γ + V(±0.4,0)，空穴 Γ）：谷归属在 Voronoi 边界上有等距平局，
    检验"作业率 = 按不可约映射展开的率"。单星模型里归属不影响率，这一项在那里不可能失败（复审指出）。"""
    return {"cell_c_A": 20.0, "layer_normal": [0, 0, 1], "kz_cap_rmax": 1,
            "eps_inf_slab": [[3, 0, 0], [0, 3, 0], [0, 0, 1.2]],
            "eps_static_slab": [[3.2, 0, 0], [0, 3.2, 0], [0, 0, 1.2]],
            "mechanisms_2d": ["ODP", "IVS"],
            "dp_optical": {
                "electron": {"valleys": [{"label": "Γ", "g": 1, "k_frac": [[0, 0, 0]]},
                                         {"label": "V", "g": 2, "k_frac": [[0.4, 0, 0], [-0.4, 0, 0]]}],
                             "channels": [{"from": "V", "to": "V", "D_eV_A": 4.0, "hw_meV": 48.0, "kind": "intra"},
                                          {"from": "V", "to": "Γ", "D_eV_A": 2.0, "hw_meV": 45.0, "kind": "inter"},
                                          {"from": "Γ", "to": "Γ", "D_eV_A": 6.0, "hw_meV": 50.0, "kind": "intra"},
                                          {"from": "Γ", "to": "V", "D_eV_A": 2.0, "hw_meV": 45.0, "kind": "inter"}]},
                "hole": {"valleys": [{"label": "Γ", "g": 1, "k_frac": [[0, 0, 0]]}],
                         "channels": [{"from": "Γ", "to": "Γ", "D_eV_A": 5.0, "hw_meV": 50.0, "kind": "intra"}]},
                "adp_intravalley": {"electron": False, "hole": False}}}


def run_one(run, mask, nworkers):
    """子进程入口：在 run 目录里放插件 + 记录，跑一次 AMSET（插件按 2d_correction.json 在导入时定型，所以一次一进程）。"""
    os.makedirs(run, exist_ok=True)
    for f in ("amset2d_plugin.py", "dp_optical.py", "amset_ir_fix.py", "amset_fermi_fix.py"):
        shutil.copy(HERE / f, run)
    os.chdir(run)
    sys.path.insert(0, run)
    rect = mask == "rect"
    Path("2d_correction.json").write_text(json.dumps(record_rect() if rect else record(mask == "on")))
    import amset2d_plugin  # noqa: F401
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.bandstructure import BandStructure
    from pymatgen.electronic_structure.core import Spin
    from amset.core.run import Runner

    lat = Lattice.orthorhombic(3.2, 5.5, 20.0) if rect else Lattice.hexagonal(3.18, 20.0)
    st = Structure(lat, ["Mo"], [[0, 0, 0]])
    n = 30
    k1, k2, k3 = np.meshgrid(np.arange(n) / n, np.arange(n) / n, [0.0, 0.5], indexing="ij")
    kf = np.stack([k1.ravel(), k2.ravel(), k3.ravel()], 1)
    if rect:
        S = np.cos(2 * np.pi * kf[:, 0]) + np.cos(2 * np.pi * kf[:, 1])
        cb = 1.8 + 0.3 * (2 - S)
        vb = -1.0 * (2 + S)
    else:
        S = (np.cos(2 * np.pi * kf[:, 0]) + np.cos(2 * np.pi * kf[:, 1])
             + np.cos(2 * np.pi * (kf[:, 0] + kf[:, 1])))
        cb = 1.8 + 1.0 * (S + 1.5)
        vb = -0.8 * (S + 1.5)
    bs = BandStructure(kf, {Spin.up: np.array([-6.0 + 0.0 * S, vb, cb, cb + 4.0])}, lat.reciprocal_lattice,
                       efermi=0.9, structure=st)
    Cv = np.zeros((6, 6))
    Cv[0, 0] = Cv[1, 1] = 60
    Cv[0, 1] = Cv[1, 0] = 15
    Cv[5, 5] = 22.5
    Cv[2, 2] = 1.0
    Cv[3, 3] = Cv[4, 4] = 0.5
    settings = {"scattering_type": ["ADP", "ODP", "IVS"], "doping": [-1e17, 1e17], "temperatures": [300],
                "deformation_potential": (2.0, 4.0), "elastic_constant": Cv.tolist(),
                "unity_overlap": True, "write_mesh": True, "nworkers": nworkers, "interpolation_factor": 4,
                "free_carrier_screening": False, "file_format": "json", "print_log": False,
                "energy_cutoff": 1.5}
    if rect:                     # 只有 basic 机制：不起 worker，几分钟
        settings.update(scattering_type=["ODP", "IVS"], doping=[-1e17, 1e17, 1e19], temperatures=[300, 900],
                        nworkers=1)
        for k in ("deformation_potential", "elastic_constant"):
            settings.pop(k)
    ad = Runner(bs, 4, settings).run(directory=".")
    # [0034-r2] 作业里的率必须等于按不可约映射展开的率（mesh 只存代表点，postprocess --check-reproduce 按它重积分）
    ad = ad[0] if isinstance(ad, tuple) else ad
    n_bad = 0
    for sp, r in ad.scattering_rates.items():
        r_ir = r[..., ad.ir_kpoints_idx][..., ad.ir_to_full_kpoint_mapping]
        n_bad += int((~np.isclose(r, r_ir, rtol=1e-9, atol=0)).sum())
    Path("sym_check.json").write_text(json.dumps({"n_differ": n_bad}))


def load(run):
    sys.path.insert(0, str(HERE))
    import postprocess_intrinsic as PP
    m = PP.read_mesh_h5(sorted(glob.glob(os.path.join(run, "mesh_*.h5")))[0])
    t = json.loads(Path(sorted(glob.glob(os.path.join(run, "transport_*.json")))[0]).read_text())
    return m, t


def mob(t, mech, i):
    a = np.asarray(t["mobility"][mech])
    return float(np.mean(np.diagonal(a[i, 0])[:2]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nworkers", type=int, default=2)
    ap.add_argument("--keep", action="store_true", help="保留临时运行目录")
    ap.add_argument("--_child", nargs=2, metavar=("RUN", "MASK"), help=argparse.SUPPRESS)
    a = ap.parse_args()
    if a._child:
        return run_one(a._child[0], a._child[1], a.nworkers)

    sys.path.insert(0, str(HERE))
    import dp_optical as D
    from pymatgen.electronic_structure.core import Spin

    tmp = tempfile.mkdtemp(prefix="dpo_e2e_")
    runs = {}
    for tag in ("on", "off", "rect"):
        run = os.path.join(tmp, "mask_" + tag)
        log = os.path.join(tmp, tag + ".log")
        print("[..] AMSET（ADP 谷内掩码 %s）-> %s" % (tag, run), flush=True)
        with open(log, "w") as fh:
            rc = subprocess.call([sys.executable, "-W", "ignore", __file__, "--nworkers", str(a.nworkers),
                                  "--_child", run, tag], stdout=fh, stderr=subprocess.STDOUT)
        txt = Path(log).read_text(errors="ignore")
        ok(rc == 0, "AMSET 运行（掩码 %s）退出码 %d（日志 %s）" % (tag, rc, log))
        if rc != 0:
            print(txt[-3000:])
            continue
        nd = json.loads(Path(run, "sym_check.json").read_text())["n_differ"]
        ok(nd == 0, "%s：作业里的散射率 = 按不可约映射展开的率（不一致 %d 处；不一致会让 --check-reproduce 失败）"
           % (tag, nd))
        if tag == "rect":        # 矩形两星模型只为对称化这一项（修复前：1044 处不一致）
            continue
        runs[tag] = load(run)
        if tag == "on":
            ok("ADP 谷内掩码：电子 开，空穴 关" in txt and "K>K':3@40" in txt, "横幅写明通道与掩码")
            ok("谷表核对（electron）：K 谷底离中心 0.000|b1|" in txt, "谷表核对：谷底落在谷中心")
            ok("FD 窗口外保留解析率" in txt, "FD 窗口外保留 ODP/IVS 解析率")
    if len(runs) == 2:
        (m, t), (m0, t0) = runs["on"], runs["off"]
        lab = m["scattering_labels"]
        rec = record(True)["dp_optical"]
        mp = m["ir_to_full_kpoint_mapping"]
        E_ir = m["energies"][Spin.up] / D.HA_EV                 # mesh 里是 eV
        E = E_ir[:, mp]
        st = m["structure"]
        rlat = st.lattice.reciprocal_lattice.matrix
        mass = float(st.composition.weight) * D.AMU
        vb = int(m["vb_idx"][Spin.up])
        R, R0 = m["scattering_rates"][Spin.up], m0["scattering_rates"][Spin.up]
        cnt = np.bincount(mp)
        for car, b, sign in (("electron", vb + 1, +1), ("hole", vb, -1)):
            c = rec[car]
            for mech, kinds in (("ODP", {"intra"}), ("IVS", {"inter_star", "inter"})):
                ref = D.carrier_rates(E[b:b + 1], m["kpoints"], rlat, c["valleys"], c["channels"], [300.0], mass,
                                      kinds, sign)[0, 0] + D.RATE_FLOOR
                ref_ir = np.bincount(mp, weights=ref) / cnt
                got = R[lab.index(mech), 0, 0, b]
                rel = np.abs(got / ref_ir - 1)
                ok(rel.max() < 1e-9, "%s %s：mesh 里的率 = dp_optical 直接算（全部 %d 个不可约点，最大相对差 %.1e）"
                   % (car, mech, len(rel), rel.max()))
            x = sign * E_ir[b]
            near = (x - x.min()) * D.HA_EV < 0.15
            r = R[lab.index("ADP"), 0, 0, b][near] / R0[lab.index("ADP"), 0, 0, b][near]
            want = 0.5 if car == "electron" else 1.0
            ok(abs(np.median(r) - want) < 0.02 and np.abs(r - want).max() < 0.05,
               "%s ADP 掩码开/关率之比 = %.3f（应 %.1f；范围 %.3f–%.3f）" % (car, np.median(r), want, r.min(), r.max()))
        mo, mi = mob(t, "ODP", 0), mob(t, "IVS", 0)
        ok(0.9 < mi / mo < 1.6, "电子单机制迁移率 IVS/ODP = %.0f/%.0f = %.2f（应 ≈1.2）" % (mi, mo, mi / mo))
        ra = mob(t, "ADP", 0) / mob(t0, "ADP", 0)
        ok(1.8 < ra < 2.2, "电子 μ_ADP 掩码开/关 = %.2f（应 ≈2）" % ra)
        ok(abs(mob(t, "ADP", 1) / mob(t0, "ADP", 1) - 1) < 0.01, "空穴 μ_ADP 不受掩码影响")
    if not a.keep and not FAILS:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\n%s" % ("全部 PASS" if not FAILS else "FAIL %d 项（运行目录保留在 %s）" % (len(FAILS), tmp)))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
