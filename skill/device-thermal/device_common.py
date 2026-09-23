#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""device_common.py -- 器件级传热公共库（纯 numpy 红黑 SOR 有限体积法）

坐标系：x = 沟道方向，y = 深度（j=0 顶面/沟道；j=Ny-1 底/衬底热沉）。
按单位宽度（W=1 m）组装，最后按 W_DEVICE 折算器件总功率。

物理：div(k(x,y) grad T) + Q_joule = 0
  - 每层各向异性热导 kx/ky；层间界面热导 TBC [W/m^2K]
  - 源项 Q = 焦耳功率密度(W/m^2) 落在沟道层
  - 底边界(Si/衬底)恒温 T_AMB；左右(金属接触)可恒温 sink / 绝热 insulator
不依赖 scipy，登录节点 numpy 即可跑。
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------- 数据
def load_materials(path=None):
    p = path or os.path.join(HERE, "device_materials.json")
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def find_upstream(cwd, skill="kl-dft-cpu", step="step6_kappa",
                  fname="kappa_summary.json"):
    """材料目录下找上游技能产物：<材料>/<skill>/<step>/<fname>（远端）
    或 <材料>/<skill>/result/<step>/<fname>（本地拉回副本）。"""
    matdir = os.path.dirname(os.path.abspath(cwd))
    for c in (os.path.join(matdir, skill, step, fname),
              os.path.join(matdir, skill, "result", step, fname),
              os.path.join(matdir, skill, fname)):
        if os.path.isfile(c):
            return c
    return None


def read_kappa(path, T):
    """读 kl-dft-cpu 的 kappa_summary.json，取最接近 T 的 κ 张量。"""
    with open(path, encoding="utf-8") as f:
        js = json.load(f)
    Ts = list(js.get("temperatures") or [])
    Ks = js.get("kappa_xx_yy_zz") or []
    if not Ts or not Ks:
        return None
    i = int(np.argmin([abs(float(t) - T) for t in Ts]))
    row = [float(x) for x in Ks[i]]
    return {"inplane": 0.5 * (row[0] + row[1]),
            "cross": (row[2] if row[2] > 0 else None),
            "T": float(Ts[i]), "tensor": row, "source": path}


def read_thickness(cwd, skill="kl-dft-cpu"):
    """读上游 2D 有效厚度 thickness_d_A（Å）-> m。"""
    matdir = os.path.dirname(os.path.abspath(cwd))
    for step in ("step6_kappa", "step4_disp",
                 "result/step6_kappa", "result/step4_disp"):
        p = os.path.join(matdir, skill, step, "thickness_2d.json")
        if os.path.isfile(p):
            with open(p, encoding="utf-8") as f:
                js = json.load(f)
            for k in ("thickness_d_A", "thickness_A", "effective_thickness_A"):
                if js.get(k):
                    return float(js[k]) * 1e-10, p
    return None, None


def resolve_tbc(db, a, b):
    t = db.get("tbc_W_m2K", {})
    for k in ("%s|%s" % (a, b), "%s|%s" % (b, a)):
        if k in t:
            return float(t[k]["G"]), t[k].get("ref", k)
    d = t.get("default")
    if d:
        return float(d["G"]), d.get("ref", "default")
    return 10e6, "hardcoded default"


# ----------------------------------------------------------------- 组装
def build_gy(Ny, dy, ky, bnd, interfaces, dx):
    """y 方向面导（单位宽度 W=1）：界面用它自己的 G*dx，其余用调和平均。"""
    gy = np.zeros(Ny - 1)
    for j in range(Ny - 1):
        if j in bnd:
            gy[j] = interfaces[bnd.index(j)] * dx
        else:
            gy[j] = dx / (dy[j] / (2.0 * ky[j]) + dy[j + 1] / (2.0 * ky[j + 1]))
    return gy


def build_coeffs(Nx, Ny, dx, dy, kx, ky, gy, q, T_amb, dirichlet_x, T_contact):
    """返回 (gx, diag, b)。gx 为 x 方向面导（只依赖 j）。"""
    gx = kx * dy / dx
    diag = np.zeros((Nx, Ny))
    b = np.zeros((Nx, Ny))
    diag[:-1, :] += gx[None, :]
    diag[1:, :] += gx[None, :]
    diag[:, :-1] += gy[None, :]
    diag[:, 1:] += gy[None, :]
    gb = 2.0 * ky[-1] * dx / dy[-1]
    diag[:, -1] += gb
    b[:, -1] += gb * T_amb
    if dirichlet_x:
        gl = 2.0 * gx
        diag[0, :] += gl
        b[0, :] += gl * T_contact
        diag[-1, :] += gl
        b[-1, :] += gl * T_contact
    b += q
    return gx, diag, b


def sor(gx, gy, diag, b, Nx, Ny, T_init=None, T0=300.0,
        omega=1.7, tol=1e-10, maxit=200000):
    """红黑排序 SOR，向量化。返回 (T, iterations, last_delta)。"""
    T = np.full((Nx, Ny), float(T0)) if T_init is None else T_init.copy()
    par = (np.arange(Nx)[:, None] + np.arange(Ny)[None, :]) % 2
    gxb = gx[None, :]
    it = 0
    dmax = 0.0
    for it in range(1, maxit + 1):
        dmax = 0.0
        for color in (0, 1):
            nb = np.zeros_like(T)
            nb[1:, :] += gxb * T[:-1, :]
            nb[:-1, :] += gxb * T[1:, :]
            nb[:, 1:] += gy[None, :] * T[:, :-1]
            nb[:, :-1] += gy[None, :] * T[:, 1:]
            tstar = (b + nb) / diag
            m = (par == color)
            upd = omega * (tstar - T)
            dmax = max(dmax, float(np.abs(upd[m]).max()))
            T = T + np.where(m, upd, 0.0)
        if dmax < tol:
            break
    return T, it, dmax


# ----------------------------------------------------------------- 器件解
def solve_device(spec, props, n_grid_x=41, n_grid_sub=60, do_transient=True,
                 sor_omega=1.7, sor_tol=1e-10, sor_maxit=200000):
    """稳态 FVM + 可选瞬态。返回 (summary dict, T[Nx,Ny], dy, lay)。"""
    ch, ox = props["channel"], props["oxide"]
    layers = [
        dict(name=ch["material"], t=ch["thickness_m"], kx=ch["kx_W_mK"],
             ky=ch["ky_W_mK"], rho_cp=ch["rho_cp_J_m3K"]),
        dict(name=ox["material"], t=ox["thickness_m"], kx=ox["kx_W_mK"],
             ky=ox["ky_W_mK"], rho_cp=ox["rho_cp_J_m3K"]),
    ]
    interfaces = [props["tbc_W_m2K"]]
    n_cells = [2, int(n_grid_sub)]
    L = float(spec["L_channel_m"])
    Nx = int(n_grid_x)
    dx = L / (Nx - 1)
    dys, lay = [], []
    for Li, (lyr, n) in enumerate(zip(layers, n_cells)):
        for _ in range(int(n)):
            dys.append(lyr["t"] / n)
            lay.append(Li)
    dy = np.array(dys)
    lay = np.array(lay)
    Ny = len(dy)
    kx = np.array([layers[l]["kx"] for l in lay])
    ky = np.array([layers[l]["ky"] for l in lay])
    rho = np.array([layers[l]["rho_cp"] for l in lay])
    bnd, c = [], 0
    for n in n_cells[:-1]:
        c += int(n)
        bnd.append(c - 1)
    gy = build_gy(Ny, dy, ky, bnd, interfaces, dx)

    q = np.zeros((Nx, Ny))
    # 体热源按沟道各单元厚度分摊：单列总功率 = Q*dx（当量于整层厚度 t_ch），
    # 而不是给沟道每个单元各加一份 Q*dx —— 后者在沟道被离散成多层时会按层数
    # 重复计入源功率（n_cells[0]=2 时实际注入 2 倍 Q*L，dT 随之偏高约 2 倍）。
    q[:, lay == 0] = (float(spec["Q_joule_W_m2"]) * dx
                      * dy[lay == 0] / ch["thickness_m"])
    dirichlet_x = str(spec["contact_bc"]).lower() in ("sink", "fixed", "contact")
    gx, diag, b = build_coeffs(Nx, Ny, dx, dy, kx, ky, gy, q,
                               float(spec["T_amb_K"]), dirichlet_x,
                               float(spec["T_amb_K"]))
    T, it, dmax = sor(gx, gy, diag, b, Nx, Ny, T0=float(spec["T_amb_K"]),
                      omega=sor_omega, tol=sor_tol, maxit=sor_maxit)
    dT = float(T.max() - float(spec["T_amb_K"]))

    r_pp = (ch["thickness_m"] / ch["ky_W_mK"] + 1.0 / props["tbc_W_m2K"]
            + ox["thickness_m"] / ox["kx_W_mK"])
    dT_1d = float(spec["Q_joule_W_m2"]) * r_pp

    out = {
        "SOLVE_DONE": True,
        "method": "FVM red-black SOR (pure numpy)",
        "Nx": Nx, "Ny": Ny, "unknowns": Nx * Ny,
        "sor_iterations": int(it), "sor_last_delta": float(dmax),
        "T_amb_K": float(spec["T_amb_K"]),
        "dT_peak_K": dT, "dT_peak_1D_network_K": dT_1d,
        "dT_ratio_FVM_over_1D": (dT / dT_1d if dT_1d else None),
        "R_th_eff_m2K_per_W": r_pp, "G_eff_W_m2K": 1.0 / r_pp,
        "power_uW": float(spec["Q_joule_W_m2"]) * L * float(spec["W_device_m"]) * 1e6,
        "contact_bc": spec["contact_bc"], "dirichlet_x": dirichlet_x,
        "layer_names": [l["name"] for l in layers],
        "layer_thickness_nm": [l["t"] * 1e9 for l in layers],
    }
    jm = Nx // 2
    out["y_nm"] = [float(np.cumsum(dy)[j] - dy[j] / 2) * 1e9 for j in range(Ny)]
    out["T_profile_mid_K"] = [float(T[jm, j]) for j in range(Ny)]
    out["T_x_mid_K"] = [float(T[i, 0]) for i in range(Nx)]

    if do_transient:
        # 器件自热有【两个尺度】的模态：
        #   fast : 沟道自身热容 + 跨 TBC 的局部 RC   tau_fast ~ (t_ch/k_ch + 1/G) * rho_cp_ch*t_ch （亚 ns）
        #   slow : 整个氧化层/衬底热容被加热           tau_RC  ~ R_th_total * C_total        （百 ns 量级）
        # 用【几何时间步】(dt_{k+1}=dt_k*RATIO) 一张网格同时分辨两者；后向欧拉无条件稳定。
        cap = rho * dy * dx
        cap2 = np.tile(cap, (Nx, 1))
        R_th = r_pp / L
        C_tot = float(np.sum(rho * dy) * L)
        tau_rc = max(R_th * C_tot, 1e-15)
        R_fast = ch["thickness_m"] / ch["ky_W_mK"] + 1.0 / props["tbc_W_m2K"]
        C_ch = ch["rho_cp_J_m3K"] * ch["thickness_m"]
        tau_fast = max(R_fast * C_ch, 1e-15)
        ratio = 1.2
        dt = tau_fast / 100.0
        T0 = np.full((Nx, Ny), float(spec["T_amb_K"]))
        Tcur = T0.copy()
        ts, dts = [0.0], [0.0]
        t = 0.0
        for _ in range(2000):
            de = diag + cap2 / dt
            be = b + cap2 / dt * Tcur
            Tcur, _i, _d = sor(gx, gy, de, be, Nx, Ny, T_init=Tcur,
                               T0=float(spec["T_amb_K"]), omega=1.0,
                               tol=1e-7, maxit=20000)
            t += dt
            ts.append(t)
            dts.append(float(Tcur.max() - float(spec["T_amb_K"])))
            dt *= ratio
            if t > 3.0 * tau_rc:
                break
        dtt = np.array(dts)
        tt = np.array(ts)
        idx63 = int(np.argmax(dtt >= 0.632 * dtt[-1])) if dtt[-1] > 0 else 0
        out["tau_fast_s"] = tau_fast
        out["tau_RC_s"] = tau_rc
        out["transient_dt0_s"] = tau_fast / 100.0
        out["transient_dt_ratio"] = ratio
        out["transient_nsteps"] = int(len(ts) - 1)
        out["tau_63_s"] = float(tt[idx63])
        out["heating_t_s"] = tt.tolist()
        out["heating_dT_K"] = dtt.tolist()
    return out, T, dy, lay
