#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step1_elec.py -- 读 ke 的 transport.json，固化成电学输入 elec_props.json

autozt 约定：run: gen 步骤的 cwd = 材料/device-zt/ ；产物 step1_elec/elec_props.json。

★ 口径：ke 的 transport.json 里 sigma/kappa_e 是【含真空超胞口径】
  （BoltzTraP 步自己的注释：「σ、κ_e 仍按超胞体积口径（与 amset 一致）」）。
  器件用面电导 G_s = sigma_cell * h_perp，与"物理层厚口径 sigma_2D*d"完全等价。
  h_perp 从 transport.json 同目录的 2d_correction.json 的 cell_c_A 取。
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stepconf          # noqa: E402
import zt_common as ZC   # noqa: E402

OUTDIR = "step1_elec"
STEP = "step1_elec"
SPEC = {
    "DEVICE_NAME": ("device", "str"),
    "KE_UPSTREAM_SKILL": ("ke-dft-cpu", "str"),
    "KE_UPSTREAM_STEP": ("", "str"),        # 空 = 自动依次找 step8_amset / step8.4_amset2d / step8.1_boltztrap
    "DOPING_CM3": (1e20, "float"),          # signed；负= p 型。按就近选 transport 里的掺杂曲线
    "T_REF": (300.0, "float"),              # 参考温度（仅用于报告；自洽用 T_eff）
    "INPLANE_MODE": ("avg", "str"),         # avg=(xx+yy)/2 | xx
    "SIGMA_T_CONV_M": (None, "float"),      # 可选：显式覆盖面电导口径厚度(m)。空=自动按胞口径
    "UPSTREAM_SKILL": ("kl-dft-cpu", "str"),   # 晶格热导上游（同 device-thermal）
    "UPSTREAM_STEP": ("step6_kappa", "str"),
}


def _cell_height_A(tp):
    """从 transport.json 同目录的 2d_correction.json 取超胞高 cell_c_A（Å）。"""
    p = os.path.join(os.path.dirname(os.path.abspath(tp)), "2d_correction.json")
    if not os.path.isfile(p):
        return None, None
    try:
        with open(p, encoding="utf-8") as f:
            rec = json.load(f)
        c = rec.get("cell_c_A")
        return (float(c) if c else None), p
    except (OSError, ValueError, TypeError):
        return None, None


def main():
    cwd = os.getcwd()
    conf = stepconf.load(SPEC, STEP, strict=False)
    out = os.path.join(cwd, OUTDIR)
    os.makedirs(out, exist_ok=True)

    ke_step = ZC.blank(conf["KE_UPSTREAM_STEP"])   # step.conf 的 "" 会解析成两个字面量引号
    tr_path, tr_step = ZC.find_transport(
        cwd, conf["KE_UPSTREAM_SKILL"], (ke_step,) if ke_step else None)
    if not tr_path:
        sys.exit("[ERROR] 找不到 ke 的 transport.json（技能=%s，步骤=%s）。"
                 "先跑 ke-dft-cpu 的 step8_amset 或 step8.4_amset2d。"
                 % (conf["KE_UPSTREAM_SKILL"], ke_step or "auto"))

    tr = ZC.read_transport(tr_path, float(conf["DOPING_CM3"]), conf["INPLANE_MODE"])

    # ---- 面电导口径厚度 t_conv ----
    t_conv, conv, conv_src = None, None, None
    if conf["SIGMA_T_CONV_M"]:
        t_conv, conv, conv_src = float(conf["SIGMA_T_CONV_M"]), "override(step.conf SIGMA_T_CONV_M)", "step.conf"
    else:
        c_A, cp = _cell_height_A(tr_path)
        if c_A:
            t_conv, conv, conv_src = c_A * 1e-10, "cell(含真空超胞高 h_perp)", cp
    if t_conv is None:
        sys.exit("[ERROR] 定不出面电导口径厚度：transport.json 同目录没有 2d_correction.json"
                 "（取 cell_c_A），也没设 SIGMA_T_CONV_M。3D 材料请显式设 SIGMA_T_CONV_M=1（"
                 "此时 sigma 是体电导，G_s=sigma*1 即单位面积的体电导）。")

    rec = {
        "ELEC_DONE": True,
        "device_name": conf["DEVICE_NAME"],
        "transport_source": tr_path,
        "transport_step": tr_step,
        "doping_cm3": tr["doping_cm3"],
        "doping_requested_cm3": tr["doping_requested_cm3"],
        "doping_all_cm3": tr["doping_all_cm3"],
        "T_ref_K": float(conf["T_REF"]),
        "inplane_mode": conf["INPLANE_MODE"],
        "T_K": tr["T_K"],
        "sigma_S_m": tr["sigma_S_m"],
        "seebeck_uV_K": tr["seebeck_uV_K"],
        "seebeck_V_K": tr["seebeck_V_K"],
        "kappa_e_W_mK": tr["kappa_e_W_mK"],
        "seebeck_unit": tr["seebeck_unit"],
        "t_conv_m": t_conv,
        "sigma_convention": conv,
        "sigma_convention_source": conv_src,
        "G_sheet_S": [ZC.sheet_conductance(s, t_conv) for s in tr["sigma_S_m"]],
        "R_sheet_ohm_sq": [ZC.sheet_resistance(s, t_conv) for s in tr["sigma_S_m"]],
        "upstream_skill_kl": conf["UPSTREAM_SKILL"],
        "upstream_step_kl": conf["UPSTREAM_STEP"],
        # 全掺杂扫描（10 档都在 transport.json 里）：器件 ZT 对掺杂极敏感，
        # 重掺下 S 塌到几 uV/K、ZT 归零，所以要把"哪档掺杂最好"一并报出来。
        "doping_scan": ZC.read_all_dopings(tr_path, float(conf["T_REF"]), conf["INPLANE_MODE"]),
    }
    p = os.path.join(out, "elec_props.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(rec, f, ensure_ascii=False, indent=2)

    i = ZC._nearest(tr["T_K"], float(conf["T_REF"]))
    print("[OK] %s" % p)
    print("    transport : %s (step=%s)" % (tr_path, tr_step))
    print("    掺杂      : %.4g cm^-3（请求 %.4g；可选 %d 档）"
          % (tr["doping_cm3"], tr["doping_requested_cm3"], len(tr["doping_all_cm3"])))
    print("    口径      : %s  t_conv=%.4e m  [%s]" % (conv, t_conv, conv_src))
    print("    T=%.0fK   : sigma=%.4g S/m  S=%.4g uV/K  kappa_e=%.4g W/mK"
          % (tr["T_K"][i], tr["sigma_S_m"][i], tr["seebeck_uV_K"][i], tr["kappa_e_W_mK"][i]))
    print("                G_sheet=%.4g S (=%.3f kOhm/sq)  面内热导占比 kappa_e/(kappa_L+kappa_e) 见 S2"
          % (rec["G_sheet_S"][i], rec["R_sheet_ohm_sq"][i] / 1000.0))


if __name__ == "__main__":
    main()
