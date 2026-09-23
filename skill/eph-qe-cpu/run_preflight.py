#!/usr/bin/env python3
"""Record the QE/Perturbo job environment and fail clearly if it is incomplete."""
from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

required = ["pw.x", "ph.x", "pw2wannier90.x", "wannier90.x", "qe2pert.x", "perturbo.x"]
executables = {name: shutil.which(name) for name in required}
mpi_launcher = shutil.which("mpirun")
pseudo_dir = os.environ.get("PSEUDO_DIR", "")
manifest_path = Path("eph_parameters.json")
manifest = {}
manifest_error = ""
if manifest_path.is_file():
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        manifest_error = f"eph_parameters.json: {exc}"
else:
    manifest_error = "缺少 eph_parameters.json"
if not isinstance(manifest, dict):
    manifest_error = manifest_error or "manifest 根对象不是 JSON 映射"
    manifest = {}
pseudo_files = manifest.get("pseudo_files", {}) if isinstance(manifest, dict) else {}
if not isinstance(pseudo_files, dict):
    manifest_error = manifest_error or "manifest 的 pseudo_files 不是映射"
    pseudo_files = {}
required_pseudo = {
    str(element): str(filename) for element, filename in pseudo_files.items()
}
pseudo_paths = {
    element: str(Path(pseudo_dir) / filename) if pseudo_dir else ""
    for element, filename in required_pseudo.items()
}
report = {
    "executables": executables,
    "mpirun": mpi_launcher,
    "PSEUDO_DIR": pseudo_dir,
    "prefix": manifest.get("prefix"),
    "species": manifest.get("species", []),
    "pseudo_files": required_pseudo,
    "pseudo_paths": pseudo_paths,
    "manifest_error": manifest_error,
    "ok": (
        all(executables.values()) and bool(mpi_launcher) and not manifest_error
        and bool(pseudo_paths) and all(Path(path).is_file() for path in pseudo_paths.values())
    ),
}
Path("qe_env.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
if not report["ok"]:
    sys.exit(2)
Path("preflight.ok").write_text("OK\n", encoding="utf-8")
