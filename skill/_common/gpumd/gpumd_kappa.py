#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gpumd_kappa.py —— GPUMD MD-κ：NVT 平衡 → HNEMD（x/y/z 三次驱动）(+ 可选 EMD) → kappa_summary.json

HNEMD 约定（GPUMD v5.6 文档）：compute_hnemd <输出间隔> <Fe_x> <Fe_y> <Fe_z>，Fe 单位 Å^-1，
结果写 kappa.out：驱动方向 μ 给 5 列 (κ_μx^in, κ_μx^out, κ_μy^in, κ_μy^out, κ_μz^tot)。
对角元因此需要三次运行：x 驱动取 κ_xx、y 驱动取 κ_yy、z 驱动取 κ_zz。
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpumd_common as gc  # noqa: E402


def build_cell(step_dir, p):
    from ase.io import read
    ref = Path(step_dir, "structure_reference.xyz")
    if not ref.is_file():
        sys.exit("[ERROR] 缺 structure_reference.xyz（先跑 S2 的数据转换）")
    atoms = read(str(ref))
    rep = [int(x) for x in p.get("CELL_REPLICATE", [1, 1, 1])]
    if any(x > 1 for x in rep):
        atoms = atoms.repeat(rep)
        print("[cell] 扩胞 %s -> %d 原子" % (rep, len(atoms)))
    return atoms


def write_xyz_in(path, atoms):
    lat = atoms.cell.array
    head = 'Lattice="%s" Properties=species:S:1:pos:R:3' % " ".join(
        "%.10f" % x for row in lat for x in row)
    lines = [str(len(atoms)), head]
    for s, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append("%-3s %18.10f %18.10f %18.10f" % (s, pos[0], pos[1], pos[2]))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def relax_lines(p):
    """可选：MD 之前先在该势的极小点上弛豫（含变胞）。
    GPUMD v5.6 语法：minimize fire <力容差> <步数> <box_change> <hydrostatic_strain>。
    box_change=1 变胞；hydrostatic_strain=0 为全胞（各向异性），1 为静水压。"""
    if not int(p.get("RELAX_ENABLED") or 0):
        return []
    return ["minimize fire %g %d %d %d" % (float(p.get("RELAX_TOL", 1e-3)),
                                           int(p.get("RELAX_STEPS", 5000)),
                                           int(p.get("RELAX_BOX_CHANGE", 1)),
                                           int(p.get("RELAX_HYDROSTATIC", 0)))]


def t_final(p, T):
    """T_FINAL 留空或缺省 = 与 T_LIST 相同（恒温）。
    SPEC 里 T_FINAL 默认是空串，StepConf 会给成 '' 或 None，直接 float('') 会崩。"""
    v = p.get("T_FINAL")
    return float(v) if v not in (None, "") else float(T)


def run_in_hnemd(p, T, direction):
    fe = float(p["DRIVING_FORCE"])
    vec = {"x": (fe, 0.0, 0.0), "y": (0.0, fe, 0.0), "z": (0.0, 0.0, fe)}[direction]
    return "\n".join(["potential nep.txt"] + relax_lines(p) + [
        "velocity %g" % T,
        "time_step %g" % float(p.get("TIME_STEP", 1.0)),
        "ensemble nvt_ber %g %g %g" % (T, t_final(p, T),
                                             float(p.get("T_COUP", 100))),
        "dump_thermo 1000",
        "run %d" % int(p.get("EQUIL_STEPS", 200000)),
        "compute_hnemd %d %g %g %g" % (int(p.get("OUTPUT_INTERVAL", 1000)), *vec),
        "dump_thermo 1000",
        "run %d" % int(p.get("PROD_STEPS", 2000000)),
    ]) + "\n"


def run_in_emd(p, T):
    return "\n".join(["potential nep.txt"] + relax_lines(p) + [
        "velocity %g" % T,
        "time_step %g" % float(p.get("TIME_STEP", 1.0)),
        "ensemble nvt_ber %g %g %g" % (T, t_final(p, T),
                                             float(p.get("T_COUP", 100))),
        "dump_thermo 1000",
        "run %d" % int(p.get("EQUIL_STEPS", 200000)),
        "compute_hac %d %d %d" % (int(p.get("EMD_SAMPLE", 10)),
                                  int(p.get("EMD_CORR_STEPS", 100000)),
                                  int(p.get("EMD_OUTPUT", 100))),
        "dump_thermo 1000",
        "run %d" % int(p.get("EMD_PROD_STEPS", 2000000)),
    ]) + "\n"


def run_in_shc(p, T, direction):
    """谱热流（SHC / Green-Kubo 谱分解）：先 NVT 平衡，再 NVE + compute_shc。
    compute_shc <采样间隔 1-50> <Nc 100-1000> <方向 0/1/2=x/y/z> <num_omega> <max_omega(THz)>。
    产物 shc.out：ω(THz) shc_i shc_o（谱 κ(ω)）。"""
    return "\n".join(["potential nep.txt"] + relax_lines(p) + [
        "velocity %g" % T,
        "time_step %g" % float(p.get("TIME_STEP", 1.0)),
        "ensemble nvt_ber %g %g %g" % (T, t_final(p, T),
                                             float(p.get("T_COUP", 100))),
        "dump_thermo 1000",
        "run %d" % int(p.get("EQUIL_STEPS", 200000)),
        "ensemble nve",
        "dump_thermo 1000",
        "compute_shc %d %d %d %d %g" % (int(p.get("SHC_SAMPLE", 10)),
                                        int(p.get("SHC_NC", 200)),
                                        {"x": 0, "y": 1, "z": 2}[direction],
                                        int(p.get("SHC_NUM_OMEGA", 200)),
                                        float(p.get("SHC_MAX_OMEGA", 50))),
        "run %d" % int(p.get("SHC_PROD_STEPS", 2000000)),
    ]) + "\n"


def run_in_nve(p, T, refinement=1):
    """Dump only the NVE segment; dt/2 run has the same physical duration."""
    dt = float(p.get("TIME_STEP", 1.0)) / refinement
    interval = int(p.get("NVE_OUTPUT_INTERVAL", 100)) * refinement
    return "\n".join(["potential nep.txt"] + relax_lines(p) + [
        "velocity %g" % T, "time_step %g" % dt,
        "ensemble nvt_ber %g %g %g" % (T, T, float(p.get("T_COUP", 100)) * refinement),
        "run %d" % (int(p.get("EQUIL_STEPS", 200000)) * refinement),
        "ensemble nve", "dump_thermo %d" % interval,
        "run %d" % (int(p["NVE_STEPS"]) * refinement),
    ]) + "\n"


def nve_metrics(path, natoms, dt_fs, interval, expected_rows):
    """GPUMD thermo columns 2/3 are total kinetic/potential energies (eV)."""
    import numpy as np
    a = np.loadtxt(path, ndmin=2)
    if expected_rows < 10 or a.shape[0] != expected_rows or a.shape[1] < 3:
        raise ValueError("incomplete NVE thermo output")
    if natoms <= 0 or dt_fs <= 0 or interval <= 0 or not np.isfinite(a).all():
        raise ValueError("invalid NVE values")
    kinetic = float(a[:, 1].mean()) / natoms
    if kinetic <= 0 or float(a[:, 0].min()) <= 0:
        raise ValueError("nonpositive kinetic energy or temperature")
    energy = (a[:, 1] + a[:, 2]) / natoms
    time = np.arange(len(a)) * dt_fs * interval / 1000.0
    centered = time - time.mean()
    slope = float(np.dot(centered, energy - energy.mean()) / np.dot(centered, centered))
    return dict(samples=len(a), observed_duration_ps=float(time[-1]),
                energy_drift_eV_atom_ps=slope,
                drift_fraction_of_kinetic=abs(slope) * float(time[-1]) / kinetic,
                energy_span_fraction_of_kinetic=float(np.ptp(energy)) / kinetic,
                mean_temperature_K=float(a[:, 0].mean()),
                min_temperature_K=float(a[:, 0].min()),
                max_temperature_K=float(a[:, 0].max()))


def run_case(step_dir, p, case_dir, run_in, label):
    d = Path(step_dir, case_dir)
    d.mkdir(parents=True, exist_ok=True)
    previous = [d / name for name in ("thermo.out", "kappa.out", "hac.out", "shc.out", "gpumd.log", "run.in")
                if (d / name).is_file()]
    if previous:
        archive = Path(tempfile.mkdtemp(prefix="previous_run_", dir=str(d)))
        for path in previous:
            shutil.move(str(path), str(archive / path.name))
    shutil.copyfile(str(Path(step_dir, "nep.txt")), str(d / "nep.txt"))
    shutil.copyfile(str(Path(step_dir, "model.xyz")), str(d / "model.xyz"))
    (d / "run.in").write_text(run_in, encoding="utf-8")
    log = d / "gpumd.log"
    with open(log, "w") as fh:
        rc = subprocess.call(["stdbuf", "-o0", "-e0", p["GPUMD_BIN"]], cwd=str(d), stdout=fh,
                             stderr=subprocess.STDOUT)
    print("[%s] rc=%d dir=%s" % (label, rc, d))
    return d, rc


def parse_emd(d):
    """EMD/Green-Kubo：hac.out col7-11 = κ_x^in/κ_x^out/κ_y^in/κ_y^out/κ_z^tot。
    返回 tau_ps / kxx / kyy / kzz 曲线（κ_xx=col7+8, κ_yy=col9+10, κ_zz=col11）。"""
    f = Path(d, "hac.out")
    if not f.is_file():
        return None
    return gc.parse_emd(str(f))


def incomplete_hnemd_results(out, temperatures):
    """Require finite results for every requested temperature and direction.

    This checks execution completeness, not physical convergence or accuracy.
    """
    missing = []
    for temperature in temperatures:
        key = "%g" % float(temperature)
        results = out.get("temperatures", {}).get(key, {})
        for direction in ("x", "y", "z"):
            result = results.get(direction, {})
            try:
                valid = (not result.get("error")
                         and math.isfinite(float(result["kappa"]))
                         and math.isfinite(float(result["stderr"]))
                         and float(result["stderr"]) >= 0
                         and int(result["n_samples"]) > 1)
            except (KeyError, TypeError, ValueError, OverflowError):
                valid = False
            if not valid:
                missing.append("T%s_%s" % (key, direction))
    return missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step-dir", default=".")
    ap.add_argument("--summary-only", action="store_true")
    args = ap.parse_args()
    p = gc.load_params(args.step_dir)
    if not args.summary_only:
        gc.write_summary(args.step_dir, "kappa_summary.json", dict(
            KAPPA_DONE=False, reason="current calculation not completed"))
        from ase.io import read  # noqa: F401  (提前暴露 ase 缺失)
        atoms = build_cell(args.step_dir, p)
        write_xyz_in(Path(args.step_dir, "model.xyz"), atoms)
        out = dict(temperatures={})
        if int(p.get("NVE_STEPS", 0)) > 0:
            checks = {}
            for T in [float(x) for x in p.get("T_LIST", [300])]:
                for refinement in (1, 2):
                    label = "T%g_nve_dt%s" % (T, refinement)
                    d, rc = run_case(args.step_dir, p, label,
                                     run_in_nve(p, T, refinement), label)
                    try:
                        if rc != 0:
                            raise ValueError("GPUMD rc=%d" % rc)
                        result = nve_metrics(d / "thermo.out", len(atoms),
                            float(p.get("TIME_STEP", 1.0)) / refinement,
                            int(p.get("NVE_OUTPUT_INTERVAL", 100)) * refinement,
                            int(p["NVE_STEPS"]) // int(p.get("NVE_OUTPUT_INTERVAL", 100)))
                        result["time_step_fs"] = float(p.get("TIME_STEP", 1.0)) / refinement
                        result["max_drift_fraction"] = float(p.get("NVE_MAX_DRIFT_FRACTION", 0.01))
                        result["max_span_fraction"] = float(p.get("NVE_MAX_SPAN_FRACTION", 0.05))
                        result["passed"] = (
                            result["drift_fraction_of_kinetic"] <= float(p.get("NVE_MAX_DRIFT_FRACTION", 0.01))
                            and result["energy_span_fraction_of_kinetic"] <= float(p.get("NVE_MAX_SPAN_FRACTION", 0.05)))
                    except (OSError, ValueError) as exc:
                        result = dict(passed=False, error=str(exc))
                    checks[label] = result
                    gc.write_summary(args.step_dir, "nve_summary.json", dict(
                        cases=checks, physical_accuracy_validated=False))
                    if not result["passed"]:
                        gc.write_summary(args.step_dir, "kappa_summary.json", dict(
                            KAPPA_DONE=False, reason="NVE preflight failed", nve_case=label))
                        sys.exit("[ERROR] NVE preflight failed: %s" % label)
        for T in [float(x) for x in p.get("T_LIST", [300])]:
            kk = {}
            for direction in ("x", "y", "z"):
                d, rc = run_case(args.step_dir, p, "T%g_hnemd_%s" % (T, direction),
                                 run_in_hnemd(p, T, direction), "HNEMD-%s" % direction)
                kf = d / "kappa.out"
                if rc != 0 or not kf.is_file():
                    kk[direction] = dict(error="rc=%d or no kappa.out" % rc)
                    continue
                k = gc.parse_kappa_out(str(kf))
                # 对角元：x 驱动 -> kmux(=κ_xx), y 驱动 -> kmuy(=κ_yy), z 驱动 -> kmuz(=κ_zz)
                diag = {"x": "kmux", "y": "kmuy", "z": "kmuz"}[direction]
                mean, err = gc.block_average(k[diag])
                kk[direction] = dict(kappa=mean, stderr=err,
                                     n_samples=int(len(k["kmux"])))
            if p.get("EMD_ENABLED"):
                d, rc = run_case(args.step_dir, p, "T%g_emd" % T, run_in_emd(p, T), "EMD")
                if rc == 0:
                    emd = parse_emd(d)
                    kk["emd"] = (dict(tau_ps_last=float(emd["tau_ps"][-1]),
                                      kxx_tail=float(emd["kxx"][-1]),
                                      kyy_tail=float(emd["kyy"][-1]),
                                      kzz_tail=float(emd["kzz"][-1]))
                                 if emd is not None else dict(error="no hac.out"))
                else:
                    kk["emd"] = dict(error="rc=%d" % rc)
            if p.get("SHC_ENABLED"):
                # 谱量子修正（严格做法）：跑 SHC 得到 κ(ω)，用 f(ω)=x²eˣ/(eˣ-1)² 逐频修正。
                shc = {}
                for direction in ("x", "y", "z"):
                    d, rc = run_case(args.step_dir, p, "T%g_shc_%s" % (T, direction),
                                     run_in_shc(p, T, direction), "SHC-%s" % direction)
                    sf = d / "shc.out"
                    if rc != 0 or not sf.is_file():
                        shc[direction] = dict(error="rc=%d or no shc.out" % rc)
                        continue
                    omega, kspec = gc.parse_shc(str(sf))
                    corr = gc.spectral_quantum_corrected_kappa(omega, kspec, T)
                    shc[direction] = dict(spectral_correction_ratio=corr["ratio"],
                                          kappa_classical_spectral=corr["kappa_classical_spectral"],
                                          kappa_quantum_spectral=corr["kappa_quantum_spectral"])
                    if "kappa" in kk.get(direction, {}):
                        kk[direction]["kappa_spectral_quantum_corrected"] = \
                            kk[direction]["kappa"] * corr["ratio"]
                kk["shc"] = shc
            qc = gc.cv_quantum_over_classical(T, p.get("DEBYE_T"))
            kk["quantum_correction_factor"] = qc
            for ax in ("x", "y", "z"):
                v = kk.get(ax, {})
                if "kappa" in v:
                    v["kappa_quantum_corrected"] = v["kappa"] * qc
            out["temperatures"]["%g" % T] = kk
        requested = p.get("T_LIST", [300])
        out["incomplete_hnemd_cases"] = incomplete_hnemd_results(out, requested)
        out["KAPPA_DONE"] = bool(requested) and not out["incomplete_hnemd_cases"]
        out["cell_atoms"] = int(len(atoms))
        gc.write_summary(args.step_dir, "kappa_summary.json", out)
        if not out["KAPPA_DONE"]:
            sys.exit("[ERROR] HNEMD 结果不完整或非有限值：%s" %
                     ", ".join(out["incomplete_hnemd_cases"]))
    else:
        print("[summary-only] kappa_summary.json 由主流程写出，这里只做存在性检查")
        if not Path(args.step_dir, "kappa_summary.json").is_file():
            sys.exit("[ERROR] 缺 kappa_summary.json")


if __name__ == "__main__":
    main()
