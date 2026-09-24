#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step2_props.py -- 汇总材料热物性：上游 kl 实测 kappa + 材料库兜底 + TBC。

本步的参数（材料/厚度/TBC/来源）优先读本步自己的 step.conf，没设才回退 S1
固化在 device_spec.json 里的值——这样在 S2_props 上 conf --set 才真正生效，
同时沿用 S1 的默认。热物性组装统一走公共池 device_common.resolve_thermal_props，
不再本地抄一份（旧版本地副本漏了警告分支、且 upstream 模式静默回退材料库）。

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
# 只有下面这几个键允许在 S2 上覆盖；材料/温度/几何由 S1 拥有，直接读 device_spec.json。
SPEC = {
    "CHANNEL_THICKNESS": (None, "float"),
    "OXIDE_THICKNESS": (None, "float"),
    "TBC_OVERRIDE": (None, "float"),
    "KAPPA_SOURCE": (None, "str"),
    "UPSTREAM_SKILL": (None, "str"),
    "UPSTREAM_STEP": (None, "str"),
}


def _p(conf, ckey, spec, skey, default=None):
    """本步 step.conf 显式值 > S1 固化值 > 代码默认。"""
    v = conf.get(ckey)
    if v is not None:
        return v
    v = spec.get(skey) if skey else None
    return default if v is None else v


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    sp = os.path.join(cwd, PREV_SPEC, "device_spec.json")
    if not os.path.isfile(sp):
        sys.exit("[ERROR] 找不到 %s，先跑 %s" % (sp, PREV_SPEC))
    with open(sp, encoding="utf-8") as f:
        spec = json.load(f)

    ch_name = spec.get("channel_material") or "MoS2"
    ox_name = spec.get("oxide_material") or "SiO2"
    T = float(spec.get("T_amb_K", 300.0))
    ch_thick = _p(conf, "CHANNEL_THICKNESS", spec, "channel_thickness_m")
    ox_thick = _p(conf, "OXIDE_THICKNESS", spec, "oxide_thickness_m")
    tbc_override = _p(conf, "TBC_OVERRIDE", spec, "tbc_override_W_m2K")
    ksrc = _p(conf, "KAPPA_SOURCE", spec, "kappa_source", "auto")
    up_skill = _p(conf, "UPSTREAM_SKILL", spec, "upstream_skill", "kl-dft-cpu")
    up_step = _p(conf, "UPSTREAM_STEP", spec, "upstream_step", "step6_kappa")

    ch, ox, tbc, tbc_ref, upstream_info, src = DC.resolve_thermal_props(
        cwd, ch_name, ox_name, T,
        up_skill=up_skill, up_step=up_step, kappa_source=ksrc,
        thickness_override=ch_thick, tbc_override=tbc_override,
        oxide_thickness=ox_thick)

    props = {"PROPS_DONE": True, "T_amb_K": T, "channel": ch, "oxide": ox,
             "tbc_W_m2K": tbc, "tbc_ref": tbc_ref,
             "upstream": upstream_info, "kappa_source": src,
             "substrate_material": spec.get("substrate_material")}
    out = os.path.join(cwd, OUTDIR)
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, "thermal_props.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(props, f, ensure_ascii=False, indent=2)
    print("[OK] %s" % p)
    print("    channel %s: kx=%.4g ky=%.4g W/mK t=%.3f nm  [%s]"
          % (ch_name, ch["kx_W_mK"], ch["ky_W_mK"], ch["thickness_m"] * 1e9, ch["source"]))
    print("    oxide %s: kx=%.4g ky=%.4g W/mK t=%.1f nm ; TBC=%.4g MW/m2K  [%s]"
          % (ox_name, ox["kx_W_mK"], ox["ky_W_mK"], ox["thickness_m"] * 1e9,
             tbc / 1e6, tbc_ref))
    if ch.get("kappa_note"):
        print("    κ 口径: %s" % ch["kappa_note"])
    if ch.get("ky_source"):
        print("    跨面 κ: %s" % ch["ky_source"])
    if ch.get("kx_scale_note"):
        print("    厚度换算: %s" % ch["kx_scale_note"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
