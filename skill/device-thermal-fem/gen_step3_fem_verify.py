#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step3_fem_verify.py -- 交叉验证：FEM vs 解析 1D vs device-thermal 的 FVM。

产出 step3_fem_verify/fem_verification.json（marker: "FEM_VERIFY_DONE"）。
四类检查：
  1) 后端自检：sink=False（侧壁绝热）时，FEM 必须复现解析 1D 严格解；
  2) 源功率审计：实际注入功率必须等于名义 Q*L（防单元重复计数）；
  3) 网格收敛：粗/细网格的 dT 差异；
  4) 与 device-thermal 二维 FVM 逐点比对（含中轴剖面），并给出源项比值诊断。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fem_common as FC  # noqa: E402
import stepconf  # noqa: E402

OUTDIR = "step3_fem_verify"
STEP = "step3_fem_verify"
SPEC = {
    "FEM_BACKEND": ("skfem", "str"),
    "FEM_1D_TOL": (0.01, "float"),
    "FEM_POWER_TOL": (0.02, "float"),
    "FEM_FVM_TOL": (0.05, "float"),
    "FEM_MESH_TOL": (0.01, "float"),
    "FEM_MESH_CHECK": (True, "bool"),
}


def _interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i]) if xs[i + 1] > xs[i] else 0.0
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    sp = os.path.join(cwd, "step2_fem_solve", "fem_summary.json")
    if not os.path.isfile(sp):
        sys.exit("[ERROR] 找不到 %s，先跑 step2_fem_solve" % sp)
    S = FC.load_json(sp)
    model = S["model"]
    main_res = S["main"]
    checks, findings = [], []

    # 1) 后端自检（解析 1D 严格解）
    r1 = S.get("ref_1d")
    if r1:
        an = model["analytic"]["dT_distributed_K"]
        rel = abs(r1["dT_peak_K"] / an - 1.0) if an else None
        ok = rel is not None and rel <= float(conf["FEM_1D_TOL"])
        checks.append({"name": "backend_vs_analytic_1d", "pass": bool(ok),
                       "dT_fem_K": r1["dT_peak_K"], "dT_analytic_K": an,
                       "rel_err": rel, "tol": float(conf["FEM_1D_TOL"])})
        if not ok:
            findings.append("后端在侧壁绝热极限下偏离解析 1D（%.3f%%），先查后端实现"
                            % (100.0 * (rel or 0.0)))
    else:
        checks.append({"name": "backend_vs_analytic_1d", "pass": None,
                       "note": "FEM_DO_1D_CHECK=false，未做"})

    # 2) 源功率审计
    pn = main_res["P_nominal_W_per_m"]
    ps = main_res["P_source_W_per_m"]
    rel_p = abs(ps / pn - 1.0) if pn else None
    ok_p = rel_p is not None and rel_p <= float(conf["FEM_POWER_TOL"])
    checks.append({"name": "source_power_audit", "pass": bool(ok_p),
                   "P_source_W_per_m": ps, "P_nominal_W_per_m": pn, "rel_err": rel_p,
                   "tol": float(conf["FEM_POWER_TOL"])})
    if not ok_p:
        findings.append("FEM 实际注入源功率 %.6g W/m != 名义 Q*L=%.6g W/m（相对差 %s），"
                        "源项定义与器件功率口径不一致"
                        % (ps, pn, ("%.4f" % rel_p) if rel_p is not None else "n/a"))

    # 3) 网格收敛
    if bool(conf["FEM_MESH_CHECK"]):
        backend = FC.get_backend(main_res["backend"])
        seq = [{"nx": 30, "n_oxide": 50}, {"nx": 60, "n_oxide": 100},
               {"nx": 120, "n_oxide": 200}, {"nx": 240, "n_oxide": 400}]
        sweep = []
        for ms in seq:
            r = backend.solve(model, workdir=cwd, opts={"mesh": ms})
            sweep.append({"mesh": ms, "dT_peak_K": r["dT_peak_K"]})
        fine = sweep[-1]["dT_peak_K"]
        rel_m = abs(fine / main_res["dT_peak_K"] - 1.0) if main_res["dT_peak_K"] else None
        ok_m = rel_m is not None and rel_m <= float(conf["FEM_MESH_TOL"])
        checks.append({"name": "mesh_convergence", "pass": bool(ok_m), "sweep": sweep,
                       "dT_finest_K": fine, "rel_err_vs_main": rel_m,
                       "tol": float(conf["FEM_MESH_TOL"])})
        if not ok_m:
            findings.append("网格未收敛：最细网格与原网格 dT 差 %s"
                            % (("%.4f%%" % (100.0 * rel_m)) if rel_m is not None else "n/a"))
    else:
        checks.append({"name": "mesh_convergence", "pass": None, "note": "FEM_MESH_CHECK=false"})

    # 4) 与 device-thermal 二维 FVM 比对
    fvm_p = None
    matdir = os.path.dirname(os.path.abspath(cwd))
    for sub in ("step3_solve", os.path.join("result", "step3_solve")):
        c = os.path.join(matdir, "device-thermal", sub, "device_thermal_summary.json")
        if os.path.isfile(c):
            fvm_p = c
            break
    fvm = None
    if fvm_p:
        F = FC.load_json(fvm_p)
        dT_fvm, dT_fem = F.get("dT_peak_K"), main_res["dT_peak_K"]
        ratio = (dT_fvm / dT_fem) if dT_fem else None
        ok_f = ratio is not None and abs(ratio - 1.0) <= float(conf["FEM_FVM_TOL"])
        prof = None
        y_nm, T_fvm = F.get("y_nm"), F.get("T_profile_mid_K")
        if y_nm and T_fvm and main_res.get("y_nodes_m"):
            ys = [v * 1e9 for v in main_res["y_nodes_m"]]
            Ts = main_res["T_mid_x_K"]
            diffs = [_interp(ys, Ts, float(y)) - float(t) for y, t in zip(y_nm, T_fvm)]
            prof = {"n": len(diffs),
                    "max_abs_K": max(abs(d) for d in diffs),
                    "rms_K": (sum(d * d for d in diffs) / len(diffs)) ** 0.5}
        fvm = {"summary_json": fvm_p, "dT_fvm_K": dT_fvm, "dT_fem_K": dT_fem,
               "ratio_fvm_over_fem": ratio, "profile_diff": prof,
               "fvm_analytic_1d_K": F.get("dT_peak_1D_network_K")}
        checks.append({"name": "fem_vs_device_thermal_fvm", "pass": bool(ok_f),
                       **{k: fvm[k] for k in ("dT_fvm_K", "dT_fem_K", "ratio_fvm_over_fem")},
                       "tol": float(conf["FEM_FVM_TOL"]), "profile_diff": prof})
        if ratio is not None and ratio > 1.3:
            near2 = abs(ratio - 2.0) < 0.3
            findings.append(
                "device-thermal 的二维 FVM 比独立 FEM 高 %.3f 倍%s；若源项按「沟道每个单元各加 "
                "Q*dx」实现，注入功率会被沟道单元数重复计入（本例沟道 2 层 -> 约 2 倍）。"
                "请核对 device_common.solve_device 的源项是否按层厚权重分摊。"
                % (ratio, "（接近 2 倍）" if near2 else ""))
    else:
        checks.append({"name": "fem_vs_device_thermal_fvm", "pass": None,
                       "note": "未找到 device-thermal 的 step3_solve/device_thermal_summary.json"})

    hard = [c for c in checks if c.get("pass") is False]
    verdict = "pass" if not hard else "review"
    res = {"FEM_VERIFY_DONE": True, "verdict": verdict,
           "backend": main_res["backend"], "checks": checks, "findings": findings,
           "fem": {"dT_peak_K": main_res["dT_peak_K"],
                   "source_ratio": main_res.get("source_ratio"),
                   "mesh": main_res["mesh"]},
           "analytic": model["analytic"], "fvm": fvm}
    p = FC.dump_json(os.path.join(cwd, OUTDIR, "fem_verification.json"), res)

    print("[OK] %s" % p)
    print("    verdict=%s" % verdict)
    for c in checks:
        mark = {True: "PASS", False: "FAIL", None: "SKIP"}[c.get("pass")]
        print("    [%s] %s" % (mark, c["name"]))
    for f in findings:
        print("    [!] %s" % f)


if __name__ == "__main__":
    main()
