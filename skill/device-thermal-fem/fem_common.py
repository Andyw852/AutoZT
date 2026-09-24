#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""device-thermal-fem 公共库：复用 device-thermal 的器件定义做独立 FEM 求解与交叉验证。

设计
----
* 不重复造器件参数：直接读 skill:device-thermal 的
    step1_device_spec/device_spec.json  （几何/热源/边界）
    step2_props/thermal_props.json      （各层 kappa/厚度/热容/界面 TBC）
* 求解器可插拔：FEM_BACKEND = skfem（默认，开源可 pip 装）| comsol（接口已留）。
* 后端契约：solve(model, workdir, opts) -> dict，字段见 backend_skfem.solve 文档。
"""
import json
import os

BACKENDS = ("skfem", "comsol")


def _first(paths):
    for p in paths:
        if p and os.path.isfile(p):
            return p
    return None


def find_artifacts(cwd, skill="device-thermal", override=None):
    """定位 device-thermal 的 device_spec.json / thermal_props.json。

    查找链：override 目录 > <材料>/<skill>/{stepdir, result/stepdir}/
    cwd 约定为 <材料>/<本技能目录>，故 <材料> = dirname(cwd)。
    """
    matdir = os.path.dirname(os.path.abspath(cwd))
    roots = []
    if override:
        roots.append(os.path.abspath(os.path.expanduser(str(override))))
    roots.append(os.path.join(matdir, skill))
    spec_c, prop_c = [], []
    for r in roots:
        for sub in ("step1_device_spec", os.path.join("result", "step1_device_spec")):
            spec_c.append(os.path.join(r, sub, "device_spec.json"))
        for sub in ("step2_props", os.path.join("result", "step2_props")):
            prop_c.append(os.path.join(r, sub, "thermal_props.json"))
    return _first(spec_c), _first(prop_c)


def _contact_is_sink(value):
    """与 device_common.contact_bc_is_dirichlet 同一套取值；拼错直接报错。"""
    v = str(value or "").strip().lower()
    if v in ("sink", "fixed", "contact"):
        return True
    if v in ("insulator", "adiabatic", "neumann"):
        return False
    raise SystemExit("[ERROR] CONTACT_BC 只能是 sink/insulator"
                     "（也接受 fixed/contact/adiabatic/neumann），收到 %r" % value)


def build_model(spec, props, backend="skfem", mesh=None):
    """把 device-thermal 的 spec+props 组装成与后端无关的 FEM 模型 dict。"""
    ch, ox = props["channel"], props["oxide"]
    g = {"L_m": float(spec["L_channel_m"]),
         "t_channel_m": float(ch["thickness_m"]),
         "t_oxide_m": float(ox["thickness_m"])}
    g["H_m"] = g["t_channel_m"] + g["t_oxide_m"]
    mesh = dict(mesh or {})
    # 界面用等效薄层表示：厚度 delta，导热系数 delta*G（层热阻 = 1/G）。
    # 薄层占据氧化物一侧的空间，故等效氧化层厚度 = t_oxide - delta。
    if "interface_layer_m" not in mesh:
        mesh["interface_layer_m"] = g["t_oxide_m"] / 1000.0
    model = {
        "source_skill": "device-thermal",
        "geometry": g,
        "layers": [
            {"name": ch["material"], "y_lo_m": 0.0, "y_hi_m": g["t_channel_m"],
             "kx_W_mK": float(ch["kx_W_mK"]), "ky_W_mK": float(ch["ky_W_mK"]),
             "rho_cp_J_m3K": ch.get("rho_cp_J_m3K")},
            {"name": ox["material"], "y_lo_m": g["t_channel_m"], "y_hi_m": g["H_m"],
             "kx_W_mK": float(ox["kx_W_mK"]), "ky_W_mK": float(ox["ky_W_mK"]),
             "rho_cp_J_m3K": ox.get("rho_cp_J_m3K")},
        ],
        "interface_G_W_m2K": float(props["tbc_W_m2K"]),
        "Q_W_m2": float(spec["Q_joule_W_m2"]),
        "T_amb_K": float(spec["T_amb_K"]),
        "contact_bc": str(spec["contact_bc"]),
        "sink": _contact_is_sink(spec["contact_bc"]),
        "W_device_m": float(spec["W_device_m"]),
        "backend": backend,
        "mesh": mesh,
    }
    model["analytic"] = analytic_1d(model)
    return model


def _r_series(model, t_ch_eff, t_ox_eff):
    ch, ox = model["layers"][0], model["layers"][1]
    # 氧化层是各向异性（hBN）时必须用跨面 ky；用面内 kx 会把 R_th 低估上百倍。
    return (t_ch_eff / ch["ky_W_mK"] + 1.0 / model["interface_G_W_m2K"]
            + t_ox_eff / ox["ky_W_mK"])


def analytic_1d(model):
    """一维串联解析（侧壁绝热）。

    series      : 技能现有 1D 网络公式，沟道热阻取整厚 t_ch/ky。
    distributed : 热源沿沟道厚度均匀分布时的严格解，沟道热阻为 t_ch/(2*ky)。
                  侧壁绝热时它就是精确解，可直接校验任何后端。
    氧化物厚度按等效薄层修正为 t_oxide - delta。
    """
    Q = model["Q_W_m2"]
    tch = model["geometry"]["t_channel_m"]
    tox_eff = model["geometry"]["t_oxide_m"] - model["mesh"]["interface_layer_m"]
    r_s = _r_series(model, tch, tox_eff)
    r_d = _r_series(model, tch / 2.0, tox_eff)
    return {"R_series_m2K_W": r_s, "R_distributed_m2K_W": r_d,
            "dT_series_K": Q * r_s, "dT_distributed_K": Q * r_d,
            "t_oxide_eff_m": tox_eff,
            "note": "distributed 为侧壁绝热的精确解；series 为技能现有 1D 网络公式"}


def analytic_1d_profile(model, y):
    """侧壁绝热 + 源在沟道内均匀分布时的解析温度（相对 T_amb）。"""
    Q = model["Q_W_m2"]
    g = model["geometry"]
    tch = g["t_channel_m"]
    d = model["mesh"]["interface_layer_m"]
    ch, ox = model["layers"][0], model["layers"][1]
    G = model["interface_G_W_m2K"]
    t0 = (Q * (g["t_oxide_m"] - d) / ox["ky_W_mK"] + Q / G
          + Q * tch / (2.0 * ch["ky_W_mK"]))
    if y <= tch:
        return t0 - Q * y * y / (2.0 * tch * ch["ky_W_mK"])
    if y <= tch + d:
        return t0 - Q * tch / (2.0 * ch["ky_W_mK"]) - Q * (y - tch) / (d * G)
    return t0 - Q * tch / (2.0 * ch["ky_W_mK"]) - Q / G - Q * (y - tch - d) / ox["ky_W_mK"]


def nominal_power_W_per_m(model):
    """名义源功率（每单位宽度）：Q 按器件面积计，总功率 = Q * L。"""
    return model["Q_W_m2"] * model["geometry"]["L_m"]


def get_backend(name):
    name = (name or "skfem").strip().lower()
    if name not in BACKENDS:
        raise SystemExit("[ERROR] FEM_BACKEND 只能是 %s，收到 %r"
                         % (" / ".join(BACKENDS), name))
    return __import__("backend_" + name)


def load_json(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def dump_json(p, obj):
    d = os.path.dirname(p)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2)
    return p
