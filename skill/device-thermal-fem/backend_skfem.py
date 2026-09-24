#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scikit-fem 后端：2D 截面各向异性稳态传热（P1 三角元）。

后端契约 solve(model, workdir=None, opts=None) -> dict
----------------------------------------------------
返回字段（S3 交叉验证只用其中一部分）：
  backend             : "skfem"
  dT_peak_K           : 热点温升（T_max - T_amb）
  T_top_max_K         : 顶面最高温（绝对温度）
  P_source_W_per_m    : 实际注入的源功率（每单位宽度），应等于 Q*L
  P_nominal_W_per_m   : 名义源功率 Q*L
  sink                : 本次求解是否用了恒温接触
  mesh                : 网格统计与参数
  y_nodes_m, T_mid_x_K: 中轴（x=L/2）剖面，供与 device-thermal 的 T_profile_mid_K 对比
  analytic_check      : 仅当 sink=False 时给出与解析 1D 的相对误差

物理：div(k grad T) + Q = 0，k 为各向异性 (kx, ky)；
界面热导 G 用等效薄层（厚度 delta，导热 delta*G，热阻 1/G）表示；
底面 T=T_amb（恒温）；左右接触 sink=恒温 T_amb / insulator=绝热。
"""
import numpy as np

DEFAULTS = {"nx": 60, "n_channel": 8, "n_interface": 4, "n_oxide": 200}


def available():
    try:
        import skfem  # noqa: F401
        return True, None
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)


def _require():
    try:
        import skfem
        return skfem
    except ImportError as exc:
        raise SystemExit(
            "[ERROR] FEM_BACKEND=skfem 需要 scikit-fem：pip install scikit-fem\n"
            "        本机报错：%s" % exc)


def solve(model, workdir=None, opts=None):
    opts = dict(opts or {})
    mesh = dict(DEFAULTS)
    mesh.update(opts.get("mesh") or {})
    skfem = _require()
    from skfem import (Basis, BilinearForm, LinearForm, ElementTriP1, MeshTri,
                       condense)
    from skfem import solve as sk_solve

    g = model["geometry"]
    L, tch, tox = g["L_m"], g["t_channel_m"], g["t_oxide_m"]
    H = g["H_m"]
    ch, ox = model["layers"][0], model["layers"][1]
    G, Q, Tamb = model["interface_G_W_m2K"], model["Q_W_m2"], model["T_amb_K"]
    sink = bool(opts.get("sink", model["sink"]))
    d_if = float(model["mesh"]["interface_layer_m"])
    k_if = d_if * G
    nx, n_ch, n_if, n_ox = (int(mesh["nx"]), int(mesh["n_channel"]),
                            int(mesh["n_interface"]), int(mesh["n_oxide"]))

    y_ch = np.linspace(0.0, tch, n_ch + 1)[:-1]
    y_if = np.linspace(tch, tch + d_if, n_if + 1)[:-1]
    y_ox = np.linspace(tch + d_if, H, n_ox + 1)
    y = np.concatenate([y_ch, y_if, y_ox])
    x = np.linspace(0.0, L, nx + 1)
    m = MeshTri.init_tensor(x, y)
    basis = Basis(m, ElementTriP1())

    @BilinearForm
    def diffusion(u, v, w):
        yy = w.x[1]
        kx = np.where(yy < tch, ch["kx_W_mK"],
                      np.where(yy < tch + d_if, k_if, ox["kx_W_mK"]))
        ky = np.where(yy < tch, ch["ky_W_mK"],
                      np.where(yy < tch + d_if, k_if, ox["ky_W_mK"]))
        return kx * u.grad[0] * v.grad[0] + ky * u.grad[1] * v.grad[1]

    @LinearForm
    def source(v, w):
        return np.where(w.x[1] < tch, Q / tch, 0.0) * v

    A = diffusion.assemble(basis)
    b = source.assemble(basis)

    nd = basis.nodal_dofs[0]
    # 注意：np.isclose 默认 atol=1e-8 m = 10 nm，在纳米器件上会误钉一片节点。
    at = max(1e-6 * H, 1e-18)
    D = nd[np.where(np.abs(m.p[1] - H) < at)[0]]
    if sink:
        # 金属接触只搭在沟道上（y < t_ch）——氧化层侧壁是绝热的。旧实现把 x=0/L
        # 的整列节点（含 300 nm 氧化层）都钉成 T_amb，等于给氧化层侧壁加了热沉。
        ch_nodes = m.p[1] < tch - 1e-18
        D = np.unique(np.concatenate([
            D,
            nd[np.where((np.abs(m.p[0] - 0.0) < at) & ch_nodes)[0]],
            nd[np.where((np.abs(m.p[0] - L) < at) & ch_nodes)[0]]]))
    T = sk_solve(*condense(A, b, x=np.full(basis.N, Tamb), D=D))

    pv = m.p[:, m.t]
    v0 = pv[:, 1, :] - pv[:, 0, :]
    v1 = pv[:, 2, :] - pv[:, 0, :]
    area = 0.5 * np.abs(v0[0] * v1[1] - v0[1] * v1[0])
    in_ch = pv.mean(axis=1)[1] < tch - 1e-18
    P_src = float(np.sum((Q / tch) * area[in_ch]))

    loc = basis.doflocs
    sel = np.where(np.abs(loc[0] - L / 2.0) < at)[0]
    ys = loc[1, sel]
    Ts = T[sel]
    o = np.argsort(ys)
    ys, Ts = ys[o], Ts[o]

    out = {
        "backend": "skfem",
        "dT_peak_K": float(T.max() - Tamb),
        "T_top_max_K": float(T.max()),
        "P_source_W_per_m": P_src,
        "P_nominal_W_per_m": Q * L,
        "source_ratio": (P_src / (Q * L) if Q * L else None),
        "sink": sink,
        "mesh": {"nx": nx, "n_channel": n_ch, "n_interface": n_if, "n_oxide": n_ox,
                 "n_nodes": int(basis.N), "n_elements": int(m.nelements),
                 "interface_layer_m": d_if},
        "y_nodes_m": [float(v) for v in ys],
        "T_mid_x_K": [float(v) for v in Ts],
    }
    if not sink:
        an = model["analytic"]["dT_distributed_K"]
        out["analytic_check"] = {
            "dT_analytic_K": an,
            "rel_err": (out["dT_peak_K"] / an - 1.0) if an else None,
        }
    return out


def solve_with_mesh_sweep(model, workdir=None, meshes=None):
    """按给定网格序列求解，返回收敛信息（用于 S2 自检）。"""
    meshes = meshes or [{"nx": 30, "n_oxide": 50}, {"nx": 60, "n_oxide": 100},
                        {"nx": 120, "n_oxide": 200}]
    out = []
    for ms in meshes:
        r = solve(model, workdir, opts={"mesh": ms})
        out.append({"mesh": ms, "dT_peak_K": r["dT_peak_K"]})
    return out
