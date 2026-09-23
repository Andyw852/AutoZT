#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_props.py -- 汇总材料热物性：上游 kl-dft-cpu 实测 kappa + 材料库兜底 + TBC。

产出 step2_props/thermal_props.json（判据 marker: "PROPS_DONE"）。
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import device_common as DC  # noqa: E402
import stepconf  # noqa: E402

OUTDIR = "step2_props"
STEP = "step2_props"
PREV_SPEC = "step1_device_spec"
SPEC = {
    "KAPPA_SOURCE": ("auto", "str"),
    "UPSTREAM_SKILL": ("kl-dft-cpu", "str"),
    "UPSTREAM_STEP": ("step6_kappa", "str"),
}


def _layer(db, name, thickness_override=None):
    m = (db.get("materials") or {}).get(name)
    if not m:
        raise SystemExit("[ERROR] 材料库没有 %s；可用：%s"
                         % (name, ", ".join(sorted(db.get("materials") or {}))))
    t = float(m["thickness_m"]) if m.get("thickness_m") else None
    if thickness_override:
        t = float(thickness_override)
    return {"material": name, "kind": m.get("kind"),
            "kx_W_mK": float(m["kappa_inplane_W_mK"]),
            "ky_W_mK": float(m["kappa_cross_W_mK"]),
            "thickness_m": t, "rho_cp_J_m3K": float(m["rho_cp_J_m3K"]),
            "note": m.get("note", ""), "source": "database"}


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    sp = os.path.join(cwd, PREV_SPEC, "device_spec.json")
    if not os.path.isfile(sp):
        sys.exit("[ERROR] 找不到 %s，先跑 %s" % (sp, PREV_SPEC))
    with open(sp, encoding="utf-8") as f:
        spec = json.load(f)
    db = DC.load_materials()
    T = float(spec["T_amb_K"])
    ch_name, ox_name = spec["channel_material"], spec["oxide_material"]
    src = str(conf["KAPPA_SOURCE"] or "auto").lower()
    up_skill = conf["UPSTREAM_SKILL"] or "kl-dft-cpu"
    up_step = conf["UPSTREAM_STEP"] or "step6_kappa"

    ch = _layer(db, ch_name)
    ox = _layer(db, ox_name)
    upstream_info = None
    if src in ("auto", "upstream", "dft"):
        up = DC.find_upstream(cwd, up_skill, up_step)
        if up:
            k = DC.read_kappa(up, T)
            if k:
                ch["kx_W_mK"] = k["inplane"]
                if k["cross"]:
                    ch["ky_W_mK"] = k["cross"]
                ch["source"] = "upstream:" + up
                ch["kappa_key"] = k.get("kappa_key")
                ch["kappa_note"] = k.get("note")
                ch["ky_source"] = ("upstream 跨面张量" if k["cross"]
                                   else "材料库值（%s）" % (k.get("cross_note") or "上游无有效跨面值"))
                th, thp = DC.read_thickness(cwd, up_skill)
                if th:
                    ch["thickness_m"] = th
                    ch["thickness_source"] = thp
                upstream_info = {"file": up, "T_K": k["T"], "tensor": k["tensor"],
                                 "kappa_key": k.get("kappa_key"), "note": k.get("note")}
        else:
            print("[warn] KAPPA_SOURCE=%s 但未找到上游 %s/%s 产物，改用材料库"
                  % (src, up_skill, up_step))
    elif src != "literature":
        raise SystemExit("[ERROR] KAPPA_SOURCE 只能是 auto/upstream/literature，收到 %r" % src)

    if spec.get("channel_thickness_m"):
        ch["thickness_m"] = float(spec["channel_thickness_m"])
        ch["thickness_source"] = "step.conf CHANNEL_THICKNESS"
    if ch.get("thickness_m") is None:
        sys.exit("[ERROR] 沟道 %s 厚度未知：请设 params.CHANNEL_THICKNESS，"
                 "或让上游 kl-dft-cpu 产出 thickness_2d.json" % ch_name)
    if ox.get("thickness_m") is None:
        ox["thickness_m"] = float(spec["oxide_thickness_m"])
        ox["thickness_source"] = "step.conf OXIDE_THICKNESS"

    tbc, tbc_ref = DC.resolve_tbc(db, ch_name, ox_name)
    if spec.get("tbc_override_W_m2K"):
        tbc = float(spec["tbc_override_W_m2K"])
        tbc_ref = "step.conf TBC_OVERRIDE (e.g. GPUMD/MD or experiment)"
    props = {"PROPS_DONE": True, "T_amb_K": T, "channel": ch, "oxide": ox,
             "tbc_W_m2K": tbc, "tbc_ref": tbc_ref,
             "upstream": upstream_info, "kappa_source": src,
             "substrate_material": spec["substrate_material"]}
    out = os.path.join(cwd, OUTDIR)
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, "thermal_props.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(props, f, ensure_ascii=False, indent=2)
    print("[OK] %s" % p)
    print("    channel %s: kx=%.4g ky=%.4g W/mK t=%.3f nm  [%s]"
          % (ch_name, ch["kx_W_mK"], ch["ky_W_mK"], ch["thickness_m"] * 1e9, ch["source"]))
    print("    oxide %s: k=%.4g W/mK t=%.1f nm ; TBC=%.4g MW/m2K  [%s]"
          % (ox_name, ox["kx_W_mK"], ox["thickness_m"] * 1e9, tbc / 1e6, tbc_ref))
    if ch.get("kappa_note"):
        print("    κ 口径: %s" % ch["kappa_note"])
    if ch.get("ky_source"):
        print("    跨面 κ: %s" % ch["ky_source"])


if __name__ == "__main__":
    main()
