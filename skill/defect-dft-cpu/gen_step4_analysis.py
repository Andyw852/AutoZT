#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S4: summary is intermediate; analysis_complete.json is published last."""
import sys, os, json, math, uuid, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import formation_energy as FE
import convex_hull as CH
import extract_band_refs as EB
import defects_common as D

STEP = "step4_analysis"

def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")

def main():
    os.makedirs(STEP, exist_ok=True)
    marker = Path(STEP) / "analysis_complete.json"
    # Retain prior evidence, but never let it certify this attempt.
    if marker.exists():
        marker.rename(marker.with_name(marker.name + ".previous-" + uuid.uuid4().hex))
    conf = D.load_stepconf()
    manifest = json.loads(Path("step2_defects/defects_manifest.json").read_text())
    if not manifest:
        raise SystemExit("[错误] 缺陷 manifest 为空")
    def energy(path, relaxed=False):
        e = FE.read_energy(str(path))
        if e is None or not math.isfinite(e):
            raise SystemExit("[错误] 缺有效能量: " + str(path))
        if relaxed and "reached required accuracy" not in Path(path).read_text(errors="replace"):
            raise SystemExit("[错误] 中性结构未弛豫收敛: " + str(path))
        return e
    summary = {"step": STEP, "status": "intermediate",
               "E_bulk": energy("step1_bulk/OUTCAR", relaxed=True), "defects": {}}
    inputs = ["step.conf", "step2_defects/defects_manifest.json",
              "step0_references/references_energy.json", "step1_bulk/OUTCAR", "energies.json"]
    modes = set()
    for item in manifest:
        base = item["dir"]
        inputs.append("step2_defects/%s/OUTCAR" % base)
        summary["defects"][base] = {"step2_defects": energy(Path("step2_defects") / base / "OUTCAR", relaxed=True)}
        for q in D.charge_states(item["name"], conf, sorted(item["counts"])):
            if q == 0:
                continue
            name = "%s_q%+d" % (base, q)
            directory = Path("step3_charged") / name
            inputs.extend(str(directory / f) for f in ("INCAR", "OUTCAR"))
            # Actual INCAR, not current step.conf, defines energy interpretation.
            import re
            incar = (directory / "INCAR").read_text()
            nsw = re.search(r"(?im)^\s*NSW\s*=\s*(\d+)", incar)
            if not nsw:
                raise SystemExit("[错误] 缺 NSW: " + str(directory))
            relaxed = int(nsw.group(1)) > 0
            if relaxed and "reached required accuracy" not in (directory / "OUTCAR").read_text(errors="replace"):
                raise SystemExit("[错误] 带电态未弛豫收敛: " + name)
            modes.add("thermodynamic" if relaxed else "vertical")
            summary["defects"][name] = {"step3_charged": energy(directory / "OUTCAR")}
    interpretation = "neutral_only" if not modes else (next(iter(modes)) if len(modes) == 1 else "mixed")
    summary["energy_interpretation"] = interpretation
    dump(Path(STEP) / "formation_energy_summary.json", summary)
    if not Path("step0_references/references_energy.json").is_file():
        raise SystemExit("[错误] 缺 step0 参考相能量")
    # No swallowed SystemExit; a stale energies/results file cannot certify success.
    CH.run_from_references()
    ref = json.loads(Path("energies.json").read_text(encoding="utf-8"))
    if str(conf.get("IS_METAL", "")).strip():
        value = str(conf["IS_METAL"]).strip().lower()
        if value not in ("1", "true", "yes", "on", "0", "false", "no", "off"):
            raise SystemExit("[错误] IS_METAL 必须为布尔值")
        ref["is_metal"] = value in ("1", "true", "yes", "on")
        ref["is_metal_source"] = "step.conf:IS_METAL"
        dump("energies.json", ref)
    if ref.get("is_metal") is not True:
        EB.main()
    else:
        print("[信息] energies.json 显式 is_metal=true：不要求半导体带边参考")
    results = FE.analyze()
    if not isinstance(results, dict) or results.get("status") != "done" or not results.get("defects"):
        raise SystemExit("[错误] 形成能分析未返回完整结果")
    results["energy_interpretation"] = interpretation
    results["thermodynamic_charge_states"] = interpretation in ("neutral_only", "thermodynamic")
    if interpretation in ("vertical", "mixed"):
        results["interpretation_warning"] = "NSW=0 charged energies are vertical; transition levels and carrier predictions are not thermodynamic equilibrium results."
        print("[警告] 带电态包含垂直能量，不可作为热力学平衡电荷转变/P-N 结论")
    # FE's root output remains compatible; fetched output lives under the step.
    dump("formation_energy_results.json", results)
    dump(Path(STEP) / "formation_energy_results.json", results)
    temporary = marker.with_name(marker.name + ".pending-" + uuid.uuid4().hex)
    # Include geometries/potentials used by volume and charge alignment, when present.
    for directory in [Path("step1_bulk")] + list(Path("step2_defects").glob("def-*")) + list(Path("step3_charged").glob("def-*")):
        for name in ("POSCAR", "CONTCAR", "LOCPOT", "EIGENVAL"):
            if (directory / name).is_file():
                inputs.append(str(directory / name))
    inputs.append(STEP + "/formation_energy_results.json")
    fingerprints = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs}
    dump(temporary, {"schema": 1, "status": "done", "energy_interpretation": interpretation,
                     "input_files": fingerprints, "result": "formation_energy_results.json"})
    temporary.replace(marker)
    print("[OK] " + str(marker))

if __name__ == "__main__":
    main()
