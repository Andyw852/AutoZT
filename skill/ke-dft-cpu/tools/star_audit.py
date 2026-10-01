#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""star_audit.py —— IR_FIX 改变迁移率的来源拆分（V130，只读两个 mesh*.h5，不需要 AMSET）。

同一张密网格上，一份是旧跑法（AMSET 原式不可约映射：星里的成员分别算），一份是 IR_FIX 跑法（代表点，
或 V127 的星平均）。按 IR_FIX 的轨道（完整点群下的星）逐轨道比较散射率，把迁移率差拆成：

    旧：        μ_old   ∝ Σ_k w_k · mean_star(1/Γ_old)        （每个成员各自 τ）
    Jensen：    μ_J     ∝ Σ_k w_k / mean_star(Γ_old)          （先平均率再取 τ —— 星平均能做到的最好）
    代表点：    μ_rep   ∝ Σ_k w_k / Γ_old(代表点)              （IR_FIX 不带星平均的预测）
    新：        μ_new   ∝ Σ_k w_k / Γ_new(代表点)

    Jensen 项 = μ_J/μ_old − 1（≤ 0）；仿真误差 = μ_new/μ_J − 1（星平均的代表点率 vs 旧成员率的平均）。
    w_k = f(1−f)/kT · |v|²（面内平均，--dims 3 时取 tr/3），同一颗星内相同。

用"代表点"一行与实测的 IR_FIX-only 结果（MoS₂：n −4.28%、p −0.04%）对照，可以检验这个近似靠不靠谱。

用法：
    python star_audit.py <旧 S8.4 目录或 mesh.h5> <新 S8.4 目录或 mesh.h5> [--mech ADP] [--temp 300]
                         [--dims 2] [--json out.json]
两份 mesh 必须是同一张密网格（KZ_CAP_2D 一致）。--mech 可以是机制名或 overall（全部机制求和）。
"""
import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "step8.4_amset2d"))

KB_EV = 8.617333262e-5


def find_mesh(p):
    p = Path(p)
    if p.is_file():
        return p
    c = sorted(glob.glob(str(p / "mesh*.h5")), key=os.path.getmtime)
    if not c:
        raise FileNotFoundError("%s 里没有 mesh*.h5（settings.yaml 要 write_mesh: true）" % p)
    return Path(c[-1])


class Mesh:
    """按需切片读 mesh.h5（散射率数组很大，只读用到的 (机制, 掺杂, 温度)）。"""

    def __init__(self, path):
        import h5py
        from postprocess_intrinsic import split_mesh_key
        self.path = Path(path)
        self.f = h5py.File(str(path), "r")
        self.keys = {}
        for raw in self.f.keys():
            base, spin = split_mesh_key(raw)
            self.keys[(base, spin)] = raw
        lab = np.asarray(self.f[self.keys[("scattering_labels", None)]][()])
        self.labels = [str(x) for x in lab.astype("U")]
        self.kpoints = np.asarray(self.f[self.keys[("kpoints", None)]][()], float)
        self.ir_kpoints = np.asarray(self.f[self.keys[("ir_kpoints", None)]][()], float)
        self.ir_to_full = np.asarray(self.f[self.keys[("ir_to_full_kpoint_mapping", None)]][()], int)
        self.temperatures = np.asarray(self.f[self.keys[("temperatures", None)]][()], float)
        self.doping = np.asarray(self.f[self.keys[("doping", None)]][()], float)
        self.fermi = np.asarray(self.f[self.keys[("fermi_levels", None)]][()], float)
        self.spins = sorted({s for (b, s) in self.keys if b == "scattering_rates" and s})

    def rates(self, spin, mech, n, t):
        ds = self.f[self.keys[("scattering_rates", spin)]]
        if mech.lower() == "overall":
            return np.asarray(ds[:, n, t], float).sum(axis=0)
        up = [x.strip().upper() for x in self.labels]
        if mech.strip().upper() not in up:
            raise KeyError("机制 %s 不在 %s 里" % (mech, self.labels))
        return np.asarray(ds[up.index(mech.strip().upper()), n, t], float)

    def ir(self, base, spin):
        return np.asarray(self.f[self.keys[(base, spin)]][()], float)


def mesh_dims(kpoints):
    return [len(np.unique(np.round(np.mod(kpoints[:, i], 1.0), 6) % 1.0)) for i in range(3)]


def grid_index(kpoints, dims):
    g = np.mod(np.rint(np.asarray(kpoints) * np.asarray(dims)), dims).astype(np.int64)
    return (g[:, 0] * dims[1] + g[:, 1]) * dims[2] + g[:, 2]


def align(old_k, new_k):
    """perm：old 全网格数组[..., perm] -> new 全网格顺序。"""
    dims = mesh_dims(new_k)
    if mesh_dims(old_k) != dims or len(old_k) != len(new_k):
        raise ValueError("两份 mesh 不是同一张密网格：%s vs %s" % (mesh_dims(old_k), dims))
    io, inew = grid_index(old_k, dims), grid_index(new_k, dims)
    pos_old = np.empty(int(np.prod(dims)), np.int64)
    pos_old[io] = np.arange(len(io))
    return pos_old[inew]


def star_sums(old_full, new_ir, orbit, rep_full, e_ir, v2_ir, efermi, T):
    """核心：全部是 numpy。old_full (nb, nk) 旧成员率（新顺序），new_ir (nb, n_orb) 新代表点率，
    orbit (nk,) 全网格点 -> 轨道号，rep_full (n_orb,) 每个轨道代表点的全网格下标，
    e_ir (nb, n_orb) eV，v2_ir (nb, n_orb) |v|²，efermi eV，T K。
    返回四种迁移率的（未归一）和、星内离散的加权和、以及 Γ新/星均 的逐点值与权重。"""
    n_orb = new_ir.shape[1]
    cnt = np.bincount(orbit, minlength=n_orb).astype(float)
    kt = KB_EV * T
    x = np.clip((e_ir - efermi) / kt, -200, 200)
    f = 1.0 / (np.exp(x) + 1.0)
    w = f * (1 - f) / kt * v2_ir * cnt[None]
    s = {"old": 0.0, "jensen": 0.0, "rep": 0.0, "new": 0.0, "spread": 0.0, "w": 0.0}
    ratios, rw = [], []
    for b in range(old_full.shape[0]):
        g = old_full[b]
        m1 = np.bincount(orbit, g, n_orb) / cnt
        minv = np.bincount(orbit, 1.0 / g, n_orb) / cnt
        m2 = np.bincount(orbit, g * g, n_orb) / cnt
        sd = np.sqrt(np.maximum(m2 - m1 * m1, 0.0))
        wb = w[b]
        s["old"] += float(np.sum(wb * minv))
        s["jensen"] += float(np.sum(wb / m1))
        s["rep"] += float(np.sum(wb / g[rep_full]))
        s["new"] += float(np.sum(wb / new_ir[b]))
        s["spread"] += float(np.sum(wb * sd / m1))
        s["w"] += float(np.sum(wb))
        ratios.append(new_ir[b] / m1)
        rw.append(wb)
    return s, np.concatenate(ratios), np.concatenate(rw)


def summarize(s, ratios, rw):
    keep = rw > 0
    order = np.argsort(ratios[keep])
    cw = np.cumsum(rw[keep][order]) / max(rw[keep].sum(), 1e-300)

    def pct(p):
        return float(ratios[keep][order][min(np.searchsorted(cw, p), len(cw) - 1)]) if keep.any() else float("nan")
    o = s["old"]
    return {"jensen_pct": 100 * (s["jensen"] / o - 1), "emulation_pct": 100 * (s["new"] / s["jensen"] - 1),
            "new_vs_old_pct": 100 * (s["new"] / o - 1), "rep_vs_old_pct": 100 * (s["rep"] / o - 1),
            "star_spread_rel": s["spread"] / max(s["w"], 1e-300),
            "gamma_new_over_star_mean": {"p5": pct(0.05), "p50": pct(0.5), "p95": pct(0.95)}}


def audit(old_path, new_path, mech="ADP", temp=300.0, dims=2):
    old, new = Mesh(find_mesh(old_path)), Mesh(find_mesh(new_path))
    perm = align(old.kpoints, new.kpoints)
    orbit = new.ir_to_full
    d = mesh_dims(new.kpoints)
    gi_full = grid_index(new.kpoints, d)
    lut = np.empty(int(np.prod(d)), np.int64)
    lut[gi_full] = np.arange(len(gi_full))
    rep_full = lut[grid_index(new.ir_kpoints, d)]              # 每个轨道代表点在全网格里的位置
    t = int(np.argmin(np.abs(new.temperatures - temp)))
    rows = []
    for n in range(len(new.doping)):
        tot, rs, ws = None, [], []
        for spin in new.spins:
            old_full = old.rates(spin, mech, n, t)[:, old.ir_to_full][:, perm]
            v = new.ir("velocities", spin)
            v2 = (v[..., 0] ** 2 + v[..., 1] ** 2) / 2 if dims == 2 else np.sum(v ** 2, axis=-1) / 3
            s, r, w = star_sums(old_full, new.rates(spin, mech, n, t), orbit, rep_full, new.ir("energies", spin),
                                v2, float(new.fermi[n, t]), float(new.temperatures[t]))
            tot = s if tot is None else {k: tot[k] + s[k] for k in tot}
            rs.append(r)
            ws.append(w)
        rows.append({"doping_cm3": float(new.doping[n]), "T": float(new.temperatures[t]),
                     **summarize(tot, np.concatenate(rs), np.concatenate(ws))})
    return {"old": str(old.path), "new": str(new.path), "mech": mech, "n_orbits": int(len(new.ir_kpoints)),
            "n_full": int(len(new.kpoints)), "rows": rows}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--mech", default="ADP")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--dims", type=int, default=2, choices=(2, 3))
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    try:
        rep = audit(a.old, a.new, a.mech, a.temp, a.dims)
    except (OSError, KeyError, ValueError) as e:
        print("[ERROR] %s" % e, file=sys.stderr)
        return 2
    print("旧：%s\n新：%s\n机制 %s，轨道 %d / 全网格 %d" % (rep["old"], rep["new"], rep["mech"], rep["n_orbits"],
                                                   rep["n_full"]))
    print("%12s %6s %10s %10s %12s %12s %10s %22s" % ("掺杂 cm^-3", "T", "Jensen%", "仿真%", "新/旧%", "代表点/旧%",
                                                      "星内离散", "Γ新/星均 p5/p50/p95"))
    for r in rep["rows"]:
        g = r["gamma_new_over_star_mean"]
        print("%12.3e %6.0f %10.3f %10.3f %12.3f %12.3f %10.3f %8.3f/%.3f/%.3f"
              % (r["doping_cm3"], r["T"], r["jensen_pct"], r["emulation_pct"], r["new_vs_old_pct"],
                 r["rep_vs_old_pct"], r["star_spread_rel"], g["p5"], g["p50"], g["p95"]))
    print("解读：新/旧 ≈ Jensen + 仿真。仿真 ≈ 0 -> 星平均忠实，剩下的只是对率平均的代价（改成逐成员 τ 才能消掉）；"
          "仿真明显 ≠ 0 -> 星平均没模仿到旧跑法（如重叠也随成员变）。代表点/旧 应接近 IR_FIX-only 实测。")
    if a.json:
        Path(a.json).write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
