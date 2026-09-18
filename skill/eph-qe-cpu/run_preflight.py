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
si_upf = str(Path(pseudo_dir) / "Si.upf") if pseudo_dir else ""
report = {
    "executables": executables,
    "mpirun": mpi_launcher,
    "PSEUDO_DIR": pseudo_dir,
    "Si_upf": si_upf,
    "ok": all(executables.values()) and bool(mpi_launcher) and bool(si_upf) and Path(si_upf).is_file(),
}
Path("qe_env.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
if not report["ok"]:
    sys.exit(2)
Path("preflight.ok").write_text("OK\n", encoding="utf-8")
