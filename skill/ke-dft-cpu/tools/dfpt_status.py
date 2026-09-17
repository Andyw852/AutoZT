import sys, os
from pathlib import Path

def status(d):
    p = Path(d)
    oc = p / "OUTCAR"
    if not oc.is_file():
        return {"dir": str(p), "state": "no-outcar"}
    txt = oc.read_text(errors="ignore")
    r = {
        "dir": str(p),
        "eps_inf_block": txt.count("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"),
        "ionic_block": txt.count("MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION"),
        "born_summary": txt.count("BORN EFFECTIVE CHARGES (including local field effects)"),
        "finished": txt.count("General timing"),
        "dof_done": txt.count("BORN EFFECTIVE CHARGE FOR ION") // 9,
    }
    if r["finished"] and r["ionic_block"] and r["eps_inf_block"]:
        r["state"] = "COMPLETE"
    elif r["eps_inf_block"]:
        r["state"] = "eps-inf-only (linear response / ionic pending)"
    else:
        r["state"] = "running (no tensor yet)"
    return r

for d in sys.argv[1:]:
    r = status(d)
    print("%-58s %s" % (os.path.basename(os.path.dirname(d)) + "/" + os.path.basename(d), r["state"]))
    print("    eps_inf=%d ionic=%d born_summary=%d finished=%d 方向完成=%d/3" % (
        r.get("eps_inf_block", 0), r.get("ionic_block", 0), r.get("born_summary", 0),
        r.get("finished", 0), r.get("dof_done", 0)))