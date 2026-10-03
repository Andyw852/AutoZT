#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""vasp_eta.py —— 在跑的 VASP 作业还要多久、会不会撞墙时（V144；只读，秒级，不需要 VASP/AMSET）。

起因：Mo2S3 S6（IBRION=6，20 原子，NFREE=4）跑满 regular 的 24 h 才被 SLURM 按 TIMEOUT 杀掉，
IBRION=6 不能续算，白跑一天。按已完成的步数外推，几个小时内就能看出要超时。

做法：读步骤目录的 INCAR / POSCAR / OUTCAR / submit.sh：
  · 已完成步数 = OUTCAR 里 LOOP+ 的个数，已用时间 = 它们 real time 之和，每步耗时取最近 10 步的中位数；
  · 总步数：IBRION=6 = 1 + NFREE×3×不等价原子数 + NFREE×6（ISIF≥3），IBRION=5 = 1 + NFREE×3N（上界）；
    其它（弛豫、DFPT IBRION=7/8、LCALCEPS）总步数事先不知道，只报每步耗时与已用时间；
  · 墙时上限：submit.sh 的 --time，否则按 --qos 查 ke_common.WALLTIME_BY_QOS_H，否则 24 h。
用法：
    python <skill>/ke-dft-cpu/tools/vasp_eta.py [步骤目录 ...]（默认当前目录） [--limit-h 48]
退出码：0 来得及或无法判断；1 预计超过上限的 80%；2 输入问题。
"""
import argparse
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE.parent), str(HERE.parent.parent / "_common" / "opt")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ke_common as kc  # noqa: E402


def _incar(path):
    out = {}
    try:
        for ln in Path(path).read_text(errors="ignore").splitlines():
            ln = ln.split("#", 1)[0].split("!", 1)[0]
            for part in ln.split(";"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    out[k.strip().upper()] = v.strip()
    except OSError:
        pass
    return out


def expected_steps(inc, poscar):
    """返回 (总步数或 None, 说明)。"""
    ib = inc.get("IBRION", "").strip()
    if ib not in ("5", "6"):
        return None, "IBRION=%s：总步数事先不知道" % (ib or "未设")
    try:
        _l, counts, _f, _s = kc.read_poscar_cell(poscar)
    except Exception as e:                                           # noqa: BLE001
        return None, "读不了 POSCAR（%s）" % e
    n = int(sum(counts))
    nfree = inc.get("NFREE", "2")
    if ib == "5":
        return kc.ibrion6_steps(n, nfree, isif=0, n_ineq=n), "IBRION=5：1 + NFREE×3N（%d 原子）" % n
    isym = inc.get("ISYM", "2").strip()
    ni = n if isym in ("0", "-1") else (kc.n_inequivalent_atoms(poscar) or n)
    isif = inc.get("ISIF", "2")
    return (kc.ibrion6_steps(n, nfree, isif, ni),
            "IBRION=6：1 + NFREE×3×不等价原子 + %s（%d 原子，不等价 %d，NFREE=%s，ISYM=%s；上界）"
            % ("NFREE×6 个应变" if int(isif) >= 3 else "无应变", n, ni, nfree, isym))


def eta(step_dir, limit_h=None):
    d = Path(step_dir)
    oc = d / "OUTCAR"
    if not oc.is_file():
        return None
    times = kc.outcar_loop_times(oc)
    inc = _incar(d / "INCAR")
    total, why = expected_steps(inc, d / "POSCAR")
    if limit_h is None:
        limit_h, lsrc = kc.walltime_limit_h(d / "submit.sh")
    else:
        lsrc = "--limit-h"
    r = {"dir": str(d), "done": len(times), "elapsed_h": round(sum(times) / 3600.0, 2), "total": total,
         "total_note": why, "limit_h": limit_h, "limit_source": lsrc}
    if times:
        r["t_step_s"] = round(statistics.median(times[-10:]), 1)
    if total and times:
        est_h, over = kc.walltime_verdict(total, r["t_step_s"], limit_h, done_steps=len(times),
                                          elapsed_s=sum(times))
        r.update({"estimate_h": round(est_h, 1), "remaining_h": round(est_h - sum(times) / 3600.0, 1),
                  "over": bool(over)})
    return r


def _report(r):
    print("== %s" % r["dir"])
    print("   已完成 %d 步，已用 %.1f h%s；墙时上限 %.0f h（%s）"
          % (r["done"], r["elapsed_h"], "，最近每步 %.1f min" % (r["t_step_s"] / 60.0) if "t_step_s" in r else "",
             r["limit_h"], r["limit_source"]))
    if r.get("total"):
        print("   总步数 %d（%s）" % (r["total"], r["total_note"]))
    else:
        print("   %s" % r["total_note"])
    if "estimate_h" in r:
        print("   %s 预计总 %.1f h，还要 %.1f h -> %s"
              % ("[WARN] ★" if r["over"] else "[OK]", r["estimate_h"], r["remaining_h"],
                 "会超过上限的 %.0f%%：现在就改 QoS/墙时重交，别等它被杀（IBRION=5/6 不能续算）"
                 % (100 * kc.WALLTIME_WARN_FRAC) if r["over"] else "来得及"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dirs", nargs="*", default=["."])
    ap.add_argument("--limit-h", type=float, default=None, help="墙时上限（小时），覆盖 submit.sh 的推断")
    a = ap.parse_args(argv)
    over = False
    seen = False
    for d in a.dirs:
        r = eta(d, a.limit_h)
        if r is None:
            print("[..] %s：没有 OUTCAR（还没开始？）" % d)
            continue
        seen = True
        _report(r)
        over |= bool(r.get("over"))
    if not seen:
        return 2
    return 1 if over else 0


if __name__ == "__main__":
    sys.exit(main())
