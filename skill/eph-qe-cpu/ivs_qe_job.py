#!/usr/bin/env python3
"""ivs_qe_job.py —— eph-qe-cpu 可选组 ivs_dfpt（[0059]）在**作业里**跑的小工具，只用标准库（计算节点不一定有 numpy）。

  python3 ivs_qe_job.py env    --pseudo-dir D --pseudos "Mo:Mo.upf S:S.upf"     # 环境探针 -> ivs_env.json（缺东西 exit 2）
  python3 ivs_qe_job.py valley --xml tmp/P.save/data-file-schema.xml --scf-out scf.out \
                               --tpl valley.in.tpl --out valley.in                   # 按 SCF 的电子数写谷能带输入
  python3 ivs_qe_job.py check  --xml tmp_valley/P.save/data-file-schema.xml --ref valley_ref.json \
                               --tol-mev 30 --mode error                            # 各谷能差 vs S8.2 -> valley_check.json
  python3 ivs_qe_job.py ahc    --dir ahc_dir --nbnd N --ahc-nbnd M --nat A --nk K     # AHC 二进制大小核对（exit 3 = 不对）

只读、只写本步骤目录里的小 JSON；任何一步失败都以非 0 退出，作业随之失败（不写完成标记）。
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET

HA_EV = 27.211386245988


def _upf_kind(path):
    """UPF 的赝势类型：NC / SL（半局域模守恒）/ US / PAW / ?。只读开头 200 KB。"""
    with open(path, "rb") as fh:
        head = fh.read(200000).decode("utf-8", "replace")
    m = re.search(r'pseudo_type\s*=\s*"\s*(\w+)\s*"', head)
    if m:
        return m.group(1).upper()
    if re.search(r'is_paw\s*=\s*"\s*(T|\.?true\.?)\s*"', head, re.I):
        return "PAW"
    if re.search(r'is_ultrasoft\s*=\s*"\s*(T|\.?true\.?)\s*"', head, re.I):
        return "US"
    m = re.search(r"<PP_HEADER>\s*\n\s*\S+.*\n\s*\S+\s.*\n\s*(\S+)", head)       # UPF v1：第三行是 NC/US/PAW
    if m:
        tok = m.group(1).upper()
        return {"NC": "NC", "SL": "SL", "US": "US", "PAW": "PAW"}.get(tok, tok)
    return "?"


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def cmd_env(a):
    exes = {n: shutil.which(n) for n in ("pw.x", "ph.x")}
    rep = {"executables": exes, "mpirun": shutil.which("mpirun"), "pseudo_dir": a.pseudo_dir, "pseudos": {},
           "problems": []}
    for n, p in exes.items():
        if not p:
            rep["problems"].append("找不到 %s（setting/<集群>.yaml 的 qe_bin，或 step.conf 的 QE_BIN / QE_SETUP）" % n)
    ph = exes.get("ph.x")
    if ph:
        try:
            with open(os.path.realpath(ph), "rb") as fh:
                blob = fh.read()
            if b"ahc_nbnd" in blob:
                rep["ahc"] = "yes"
            elif blob[:2] == b"#!" or len(blob) < 100000:
                rep["ahc"] = "unknown（ph.x 是脚本包装，作业里由 ph.x 自己报错）"
            else:
                rep["ahc"] = "no"
                rep["problems"].append("ph.x 不支持 electron_phonon='ahc'（要 QE ≥ 6.7；推荐 7.3/7.4）")
        except OSError as e:
            rep["ahc"] = "unknown（%s）" % e
    for tok in str(a.pseudos or "").replace(",", " ").split():
        el, fn = tok.split(":", 1)
        p = os.path.join(a.pseudo_dir or "", fn)
        row = {"file": fn, "path": p, "exists": os.path.isfile(p)}
        if row["exists"]:
            row["kind"] = _upf_kind(p)
            row["sha256"] = _sha(p)
            if row["kind"] not in ("NC", "SL"):
                rep["problems"].append("%s 的赝势 %s 是 %s：AHC 只支持模守恒（NC，如 PseudoDojo/SG15）" % (el, fn, row["kind"]))
        else:
            rep["problems"].append("赝势不存在：%s（setting/<集群>.yaml 的 qe_pseudo_dir 或 step.conf 的 QE_PSEUDO_DIR）" % p)
        rep["pseudos"][el] = row
    rep["ok"] = not rep["problems"]
    with open("ivs_env.json", "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=2, ensure_ascii=False)
    for p in rep["problems"]:
        print("[ERROR] " + p, file=sys.stderr)
    return 0 if rep["ok"] else 2


def _xml_info(path):
    root = ET.parse(path).getroot()
    out = root.find("output")
    bs = out.find("band_structure")

    def _b(tag):
        t = bs.find(tag)
        return t is not None and t.text.strip().lower() in ("true", "t", ".true.")
    st = out.find("atomic_structure")
    cell = st.find("cell")
    a = [[float(x) for x in cell.find(t).text.split()] for t in ("a1", "a2", "a3")]
    ks = []
    for kse in bs.findall("ks_energies"):
        kp = [float(x) for x in kse.find("k_point").text.split()]
        ev = [float(x) * HA_EV for x in kse.find("eigenvalues").text.split()]
        ks.append({"k_cart_2pi_alat": kp, "eig_eV": ev})
    nelec = float(bs.find("nelec").text)
    nbnd = int(bs.find("nbnd").text)
    return {"alat_bohr": float(st.get("alat")), "cell_bohr": a, "nat": int(st.get("nat")), "nelec": nelec,
            "nbnd": nbnd, "lsda": _b("lsda"), "noncolin": _b("noncolin"), "spinorbit": _b("spinorbit"),
            "ks": ks}


def _qe_version(out_file):
    try:
        with open(out_file, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                m = re.search(r"Program PWSCF v\.?\s*([0-9][0-9.A-Za-z-]*)", line)
                if m:
                    return m.group(1)
    except OSError:
        pass
    return None


def cmd_valley(a):
    info = _xml_info(a.xml)
    if info["lsda"] or info["noncolin"]:
        print("[ERROR] SCF 是自旋极化/非共线（lsda=%s noncolin=%s）：ivs_dfpt 只做非磁、非 SOC" % (info["lsda"], info["noncolin"]),
              file=sys.stderr)
        return 2
    nocc = info["nelec"] / 2.0
    if abs(nocc - round(nocc)) > 1e-6:
        print("[ERROR] 电子数 %.4f 是奇数：不是闭壳层半导体" % info["nelec"], file=sys.stderr)
        return 2
    nocc = int(round(nocc))
    rec = {k: info[k] for k in ("alat_bohr", "cell_bohr", "nat", "nelec", "lsda", "noncolin")}
    rec.update(nocc=nocc, nbnd_scf=info["nbnd"], qe_version=_qe_version(a.scf_out))
    with open("ivs_scf_info.json", "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2)
    with open(a.tpl, encoding="utf-8") as fh:
        txt = fh.read()
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(txt.replace("__NBND__", str(nocc + int(a.extra))))
    print("[OK] nelec=%.1f nocc=%d -> %s nbnd=%d" % (info["nelec"], nocc, a.out, nocc + int(a.extra)))
    return 0


def cmd_check(a):
    info = _xml_info(a.xml)
    with open(a.ref, encoding="utf-8") as fh:
        ref = json.load(fh)
    with open("ivs_scf_info.json", encoding="utf-8") as fh:
        nocc = int(json.load(fh)["nocc"])
    ks = info["ks"]
    if len(ks) != len(ref["points"]):
        print("[ERROR] 谷能带的 k 点数 %d ≠ 输入 %d" % (len(ks), len(ref["points"])), file=sys.stderr)
        return 2
    rep = {"tol_meV": float(a.tol_mev), "mode": a.mode, "nocc": nocc, "carriers": {}, "problems": [], "warnings": []}
    for car, stars in ref["carriers"].items():
        band = nocc if car == "electron" else nocc - 1
        sign = 1.0 if car == "electron" else -1.0
        rows = []
        for st in stars:
            es = [ks[i]["eig_eV"][band] for i in st["points"]]
            rows.append({"label": st["label"], "E_eV": sum(es) / len(es), "spread_meV": (max(es) - min(es)) * 1000.0,
                         "dE_ref_eV": st.get("dE_ref_eV")})
        best = min(sign * r["E_eV"] for r in rows)
        for r in rows:
            r["dE_qe_eV"] = sign * r["E_eV"] - best
            if r["spread_meV"] > 5.0:
                rep["warnings"].append("%s %s 星各成员能量相差 %.1f meV（应简并）：谷表与这个结构的对称性对不上？"
                                       % (car, r["label"], r["spread_meV"]))
        have = [r for r in rows if r["dE_ref_eV"] is not None]
        if have:
            ref0 = min(float(r["dE_ref_eV"]) for r in have)
            qe0 = min(r["dE_qe_eV"] for r in have)
            for r in have:
                r["diff_meV"] = ((r["dE_qe_eV"] - qe0) - (float(r["dE_ref_eV"]) - ref0)) * 1000.0
                if abs(r["diff_meV"]) > float(a.tol_mev):
                    rep["problems"].append("%s %s 谷：QE 能差 %.0f meV、S8.2 %.0f meV，差 %.0f meV > %.0f"
                                           % (car, r["label"], (r["dE_qe_eV"] - qe0) * 1000,
                                              (float(r["dE_ref_eV"]) - ref0) * 1000, r["diff_meV"], float(a.tol_mev)))
        rep["carriers"][car] = rows
    rep["status"] = "mismatch" if rep["problems"] else "ok"
    with open("valley_check.json", "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=2, ensure_ascii=False)
    for w in rep["warnings"]:
        print("[WARN] " + w)
    for p in rep["problems"]:
        print(("[ERROR] " if a.mode == "error" else "[WARN] ") + p, file=sys.stderr)
    if rep["problems"] and a.mode == "error":
        print("[ERROR] 谷能差与 S8.2 不一致（赝势/泛函/结构差异）。确认可以接受就在 step.conf 写 IVS_VALLEY_CHECK = warn "
              "或放宽 IVS_VALLEY_TOL_MEV，再 retry 本步", file=sys.stderr)
        return 3
    print("[OK] 谷能差自检：%s" % rep["status"])
    return 0


def cmd_ahc(a):
    nm = 3 * int(a.nat)
    want = {"ahc_gkk_iq1.bin": int(a.nk) * nm * int(a.ahc_nbnd) * int(a.nbnd) * 16,
            "ahc_etk_iq1.bin": int(a.nk) * int(a.nbnd) * 8, "ahc_etq_iq1.bin": int(a.nk) * int(a.nbnd) * 8}
    bad = []
    for f, n in want.items():
        p = os.path.join(a.dir, f)
        have = os.path.getsize(p) if os.path.isfile(p) else -1
        if have != n:
            bad.append("%s：%d 字节，应为 %d" % (p, have, n))
    for b in bad:
        print("[ERROR] " + b, file=sys.stderr)
    return 3 if bad else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    p = sub.add_parser("env")
    p.add_argument("--pseudo-dir", default="")
    p.add_argument("--pseudos", default="")
    p = sub.add_parser("valley")
    p.add_argument("--xml", required=True)
    p.add_argument("--scf-out", default="scf.out")
    p.add_argument("--tpl", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--extra", default="4")
    p = sub.add_parser("check")
    p.add_argument("--xml", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--tol-mev", default="30")
    p.add_argument("--mode", default="error", choices=("error", "warn"))
    p = sub.add_parser("ahc")
    for k in ("--dir", "--nbnd", "--ahc-nbnd", "--nat", "--nk"):
        p.add_argument(k, required=True)
    a = ap.parse_args(argv)
    fn = {"env": cmd_env, "valley": cmd_valley, "check": cmd_check, "ahc": cmd_ahc}.get(a.cmd)
    if fn is None:
        ap.print_help()
        return 1
    return fn(a)


if __name__ == "__main__":
    sys.exit(main())
