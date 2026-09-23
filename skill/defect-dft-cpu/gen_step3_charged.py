#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S3: charged relaxation by default; CHARGED_GEOMETRY=vertical opts into fixed-ion energies."""
import sys, os, json, shutil, hashlib, uuid
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import defects_common as D

STEP = "step3_charged"

# 电荷态清单改用 defects_common.charge_states（通用）：
# charge_states() 原先写死 Te/Sb/Bi/Pb/Sn 的价电子规则，已删除。
# 现在由 step.conf 的 CHARGE_MODE / QMIN / QMAX 决定（默认 window -2..+2）。

def nelect_neutral(order, counts, potcar_path):
    """由 POTCAR ZVAL 求中性态总电子数。"""
    zvals = D.read_nelect_from_potcar(potcar_path)
    if len(zvals) != len(order):
        raise SystemExit("[错误] POTCAR 元素数(%d) 与 POSCAR 元素数(%d) 不一致" % (len(zvals), len(order)))
    return int(round(sum(z * counts[el] for z, el in zip(zvals, order))))

def main():
    os.makedirs(STEP, exist_ok=True)
    old_manifest = Path(STEP) / "charged_manifest.json"
    if old_manifest.exists():
        old_manifest.rename(old_manifest.with_name(old_manifest.name + ".previous-" + uuid.uuid4().hex))
    conf = D.load_stepconf()
    if not os.path.isdir("step2_defects"):
        raise SystemExit("[错误] 找不到 step2_defects/ —— 先跑完 step2")
    # 中性态 NELECT 从 step2 任一 POTCAR 算
    manifest = json.loads(Path("step2_defects/defects_manifest.json").read_text(encoding="utf-8"))
    os.makedirs(STEP, exist_ok=True)
    mode = str(conf.get("CHARGED_GEOMETRY", "relax")).strip().lower()
    if mode not in ("vertical", "relax"):
        raise SystemExit("[错误] CHARGED_GEOMETRY 必须为 vertical 或 relax")
    nsw = int(conf.get("CHARGED_NSW", "60")) if mode == "relax" else 0
    if mode == "relax" and nsw <= 0:
        raise SystemExit("[错误] CHARGED_NSW 必须 > 0")
    if not manifest:
        raise SystemExit("[错误] 中性缺陷 manifest 为空")
    jobs = []
    n_jobs = 0
    for item in manifest:
        ddir = item["dir"]
        src = None
        for cand in ["step2_defects/%s/CONTCAR" % ddir]:
            if os.path.exists(cand):
                src = cand; break
        if src is None:
            raise SystemExit("[错误] %s 缺中性弛豫 CONTCAR" % ddir)
        struct = D.parse_poscar(src)
        order = [a for a in dict.fromkeys(struct["atoms"])]
        counts = {a: struct["atoms"].count(a) for a in order}
        species = sorted(set(struct["atoms"]))
        charges = [q for q in D.charge_states(item["name"], conf, species) if q != 0]
        if not charges:
            continue  # q=0 的弛豫能量来自 S2，不重复提交 VASP
        pot = "step2_defects/%s/POTCAR" % ddir
        n0 = nelect_neutral(order, counts, pot)
        for q in charges:
            qdir = "%s_q%+d" % (ddir, q)
            ne = n0 - q          # q=+1 -> 去一个电子 -> NELECT-1
            outdir = os.path.join(STEP, qdir)
            # 从中性 CHGCAR 续算（ICHARG=1）；旧流水线无 CHGCAR 时退回 ICHARG=2 从头算
            chg_src = "step2_defects/%s/CHGCAR" % ddir
            icharg = "ICHARG = 2"
            if os.path.exists(chg_src) and os.path.getsize(chg_src) > 0:
                icharg = "ICHARG = 1"
            # 奇数电子数带电态：seed 磁矩打破自旋对称（否则被算成非磁闭壳）
            natoms = len(struct["atoms"])
            soc = str(conf.get("SOC", "1")).strip().lower() not in ("0", "false", "no", "off")
            ncomp = 3 if soc else 1
            magmom = (D.spin_seed_magmom(natoms, ncomp=ncomp) if (ne % 2 != 0)
                      else "%d*0" % (ncomp * natoms))
            D.build_job(outdir, struct, conf, {
                "SYSTEM": "defect-dft-cpu step3 %s q=%+d" % (item["disp"], q),
                "IBRION": "2" if mode == "relax" else "-1",
                "ISIF": "2", "NSW": str(nsw),
                "EDIFFG_LINE": ("EDIFFG = %s" % conf.get("CHARGED_EDIFFG", "-0.02")) if mode == "relax" else "",
                "LVHAR_LINE": "LVHAR = .TRUE." if conf.get("LVHAR_CHARGED", "1") == "1" else "",
                "NELECT_LINE": "NELECT = %d" % ne,
                "ICHARG_LINE": icharg,
                "SIGMA_LINE": "SIGMA = 0.01",
                "MAGMOM": magmom,
            }, qdir)
            if icharg == "ICHARG = 1":
                shutil.copy(chg_src, os.path.join(outdir, "CHGCAR"))
            jobs.append(qdir)
            n_jobs += 1
    unexpected = {p.name for p in Path(STEP).glob("def-*") if p.is_dir()} - set(jobs)
    if unexpected:
        raise SystemExit("[错误] S3 存在旧/非本次电荷态目录，需人工检查（未删除）: " + ", ".join(sorted(unexpected)))
    if not jobs and (str(conf.get("CHARGE_MODE", "window")).strip().lower() != "window"
                     or int(conf.get("QMIN", "-2")) != 0 or int(conf.get("QMAX", "2")) != 0):
        raise SystemExit("[错误] 无带电任务时须显式设置 CHARGE_MODE=window QMIN=QMAX=0")
    inputs = ["step.conf", "step2_defects/defects_manifest.json"]
    for item in manifest:
        inputs.extend("step2_defects/%s/%s" % (item["dir"], name) for name in ("CONTCAR", "OUTCAR"))
    fingerprints = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs}
    with open(os.path.join(STEP, "charged_manifest.json"), "w", encoding="utf-8") as f:
        json.dump({"schema": 1, "status": "neutral_only" if not jobs else "generated",
                   "jobs": jobs, "geometry": mode, "neutral_source": "step2_defects",
                   "input_files": fingerprints}, f, indent=2)
    print("[OK] 生成 %s/ 下 %d 个带电态子目录" % (STEP, n_jobs))

if __name__ == "__main__":
    main()
