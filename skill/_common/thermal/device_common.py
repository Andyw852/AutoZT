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
    """读上游 kappa_summary.json，取最接近 T 的 κ 张量。

    ★ 2D 口径（关键）：phono3py/ShengBTE 算出的 κ 是【含真空的胞值】κ_cell（体积里
    带真空）。器件 FVM 用的是真实层厚 d（thickness_2d.json），面内热导 = κ×d，所以
    必须换成物理 2D 值 κ_2D = κ_cell × h⊥/d。若直接用 κ_cell，面内热导被低估
    h⊥/d 倍——单层 MoS2（h⊥≈20 Å、d≈6.7 Å）约 3×，实测 dT 会算高 ~50%。
    上游 kl-dft-cpu / kl-mlff 把这一档写在 `kappa_2d_normalized_xx_yy_zz` 里，本函数
    优先取它；缺失时退到 `kappa_2d_norm_factor` 现场归一化；再缺就原样用并警告。
    返回 dict 的 `kappa_key` / `note` 说明实际采用了哪一档。
    """
    with open(path, encoding="utf-8") as f:
        js = json.load(f)
    if not (js.get("temperatures") or []):
        # 兼容只有 runs[] 的 summary（旧版/精简产物）：优先 RTA，其次含 2D 归一化
        # 字段的那次运行；顶层只有 kappa_2d_norm_factor 时把它并进选中的 run。
        _top = js
        _runs = [r for r in (js.get("runs") or []) if r.get("temperatures")]
        if _runs:
            _prim = next((r for r in _runs
                          if str(r.get("method", "")).lower() == "rta"), None)
            if _prim is None:
                _prim = next((r for r in _runs
                              if r.get("kappa_2d_normalized_xx_yy_zz")), None)
            _prim = dict(_prim or _runs[-1])
            for _k in ("kappa_2d_norm_factor", "dim"):
                if _prim.get(_k) is None and _top.get(_k) is not None:
                    _prim[_k] = _top[_k]
            js = _prim
    Ts = list(js.get("temperatures") or [])
    if not Ts:
        return None
    dim = str(js.get("dim") or "").lower()
    i = int(np.argmin([abs(float(t) - T) for t in Ts]))
    kn = js.get("kappa_2d_normalized_xx_yy_zz") or []
    if kn:
        i = min(i, len(kn) - 1)
        row = [float(x) for x in kn[i]]
        key = "kappa_2d_normalized_xx_yy_zz"
        note = "2D 物理值 κ_2D = κ_cell × h⊥/d（与 thickness_2d.json 的 d 配套）"
    else:
        kr = js.get("kappa_xx_yy_zz") or []
        if not kr:
            return None
        i = min(i, len(kr) - 1)
        row = [float(x) for x in kr[i]]
        key = "kappa_xx_yy_zz"
        fac = js.get("kappa_2d_norm_factor")
        if dim == "2d" and fac and abs(float(fac) - 1.0) > 1e-9:
            row = [v * float(fac) for v in row]
            key = "kappa_xx_yy_zz × kappa_2d_norm_factor"
            note = ("2D：文件缺 kappa_2d_normalized_*，用 kappa_2d_norm_factor=%.4g 现场归一化"
                    % float(fac))
        elif dim == "2d":
            note = ("⚠️ 2D 材料，但文件既无 kappa_2d_normalized_* 也无 kappa_2d_norm_factor："
                    "当前是含真空的胞值，面内热导可能被低估 h⊥/d 倍（单层可达 ~3×）；"
                    "请重跑上游 kl 步补上归一化字段。")
        else:
            note = "3D 材料，κ 可直接用（无 2D 归一化）"
    inplane = 0.5 * (row[0] + row[1])
    # 跨面 κ 可用性：真空胞的 zz 是【数值零】（真实 2D 产物实测 4.6e-30），
    # 不能用 >0 当判据——否则 ky 被设成 1e-30，沟道跨面热阻炸到 1e19 m2K/W
    # （氧化层的 ~1e26 倍），器件近似绝热、dT 直接失控。
    cross, cross_note = float(row[2]), None
    if not (cross > 1e-6):
        cross, cross_note = None, "跨面 κ 数值为零/过小，保留材料库值"
    elif dim == "2d" and cross < 1e-3 * max(inplane, 1e-30):
        cross_note = ("2D 跨面 κ=%.3g 远小于面内 %.3g，判定为真空残留，保留材料库值"
                      % (float(row[2]), inplane))
        cross = None
    return {"inplane": inplane, "cross": cross,
            "T": float(Ts[i]), "tensor": row, "source": path,
            "kappa_key": key, "note": note, "cross_note": cross_note}


def read_thickness(cwd, skill="kl-dft-cpu"):
    """读上游 2D 有效厚度 thickness_d_A（Å）-> m。

    查找顺序：<材料>/<skill>/{step6_kappa, step4_disp, step4_kappa}/（各含 result/ 副本），
    最后兜底 <材料>/thickness_2d.json（kl/ke 共用契约的材料级文件）。
    step4_kappa 对应 kl-mlff-* 那条链（它的 κ 步就叫 step4_kappa）。
    """
    matdir = os.path.dirname(os.path.abspath(cwd))
    cands = [os.path.join(matdir, skill, step, "thickness_2d.json")
             for step in ("step6_kappa", "step4_disp", "step4_kappa",
                          "result/step6_kappa", "result/step4_disp",
                          "result/step4_kappa")]
    cands.append(os.path.join(matdir, "thickness_2d.json"))   # 材料级兜底
    for p in cands:
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
    """y 方向面导（单位宽度 W=1），按两层串联热阻取串联平均：

        R_j = dy[j]/(2*k_j) + [界面加 1/G] + dy[j+1]/(2*k_{j+1})
        gy[j] = dx / R_j

    界面只写 G*dx 会漏掉两侧各半个单元的热阻：绝热边界下的精确 1D 解是
    Q*(t_ox/ky_ox + 1/G + t_ch/(2*ky_ch))，旧实现低约 0.6%。"""
    gy = np.zeros(Ny - 1)
    for j in range(Ny - 1):
        r = dy[j] / (2.0 * ky[j]) + dy[j + 1] / (2.0 * ky[j + 1])
        if j in bnd:
            r += 1.0 / interfaces[bnd.index(j)]
        gy[j] = dx / r
    return gy


def build_coeffs(Nx, Ny, dx, dy, kx, ky, gy, q, T_amb, dirichlet_x, T_contact,
                 contact_rows=None):
    """返回 (gx, diag, b)。gx 为 x 方向面导（只依赖 j）。

    contact_rows：施加 x 向恒温（Dirichlet）的行掩码（长度 Ny 的布尔数组）。
    真实器件里金属接触只搭在沟道上，氧化层侧壁是绝热的；把整个侧壁都设成
    恒温会凭空多出一个氧化层热沉，改变结果。None = 所有行（旧行为，仅对照）。
    """
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
        rows = (np.ones(Ny, dtype=bool) if contact_rows is None
                else np.asarray(contact_rows, dtype=bool))
        diag[0, rows] += gl[rows]
        b[0, rows] += gl[rows] * T_contact
        diag[-1, rows] += gl[rows]
        b[-1, rows] += gl[rows] * T_contact
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


# ----------------------------------------------------------------- 输入校验
_CONTACT_SINK = ("sink", "fixed", "contact")
_CONTACT_INSULATOR = ("insulator", "adiabatic", "neumann")


def contact_bc_is_dirichlet(value):
    """CONTACT_BC 的合法取值 -> bool；拼错的值直接报错，而不是静默当绝缘。

    旧实现用 in ("sink","fixed","contact") 判断，于是 CONTACT_BC=sinks 这种
    拼写错误会被静默当成 insulator（结果完全不同）。
    """
    v = str(value or "").strip().lower()
    if v in _CONTACT_SINK:
        return True
    if v in _CONTACT_INSULATOR:
        return False
    raise SystemExit("[ERROR] CONTACT_BC 只能是 sink/insulator"
                     "（也接受 fixed/contact/adiabatic/neumann），收到 %r" % value)


def validate_device_inputs(spec, n_grid_x=None, n_grid_sub=None):
    """求解前挡下非法输入：Q<0、L<=0、Nx<2、拼错的 CONTACT_BC。"""
    try:
        L = float(spec["L_channel_m"])
    except (KeyError, TypeError, ValueError):
        raise SystemExit("[ERROR] L_CHANNEL 缺失或非数值：%r" % spec.get("L_channel_m"))
    if L <= 0:
        raise SystemExit("[ERROR] L_CHANNEL 必须 > 0，收到 %r" % spec.get("L_channel_m"))
    try:
        q = float(spec["Q_joule_W_m2"])
    except (KeyError, TypeError, ValueError):
        raise SystemExit("[ERROR] Q_JOULE 缺失或非数值：%r" % spec.get("Q_joule_W_m2"))
    if q < 0:
        raise SystemExit("[ERROR] Q_JOULE 必须 >= 0，收到 %r" % spec.get("Q_joule_W_m2"))
    contact_bc_is_dirichlet(spec.get("contact_bc"))
    if n_grid_x is not None and int(n_grid_x) < 2:
        raise SystemExit("[ERROR] N_GRID_X 必须 >= 2，收到 %r" % n_grid_x)
    if n_grid_sub is not None and int(n_grid_sub) < 1:
        raise SystemExit("[ERROR] N_GRID_SUB 必须 >= 1，收到 %r" % n_grid_sub)


# ----------------------------------------------------------------- 器件解
def solve_device(spec, props, n_grid_x=41, n_grid_sub=60, do_transient=True,
                 sor_omega=1.7, sor_tol=1e-10, sor_maxit=200000):
    """稳态 FVM + 可选瞬态。返回 (summary dict, T[Nx,Ny], dy, lay)。

    离散：x 为**单元中心**（Nx 个单元、间距 dx=L/Nx，壁面在 x=0/L 的半单元外），
    因此注入功率严格等于 Q*L（旧实现节点按顶点排 dx=L/(Nx-1)，注入 Q*L*Nx/(Nx-1)）。
    y 为节点（按层分格），底边界恒温、金属接触只搭在沟道行上。
    """
    validate_device_inputs(spec, n_grid_x, n_grid_sub)
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
    dx = L / Nx            # 单元中心间距；壁面在 0/L，距边界单元中心 dx/2
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
    dirichlet_x = contact_bc_is_dirichlet(spec["contact_bc"])
    # 金属接触只搭在沟道（lay==0）上；氧化层侧壁绝热。旧实现把整列（含 300 nm
    # 氧化层）都钉成 T_amb，相当于给氧化层侧壁凭空加了个热沉，改变了结果。
    contact_rows = (lay == 0)
    gx, diag, b = build_coeffs(Nx, Ny, dx, dy, kx, ky, gy, q,
                               float(spec["T_amb_K"]), dirichlet_x,
                               float(spec["T_amb_K"]), contact_rows=contact_rows)
    T, it, dmax = sor(gx, gy, diag, b, Nx, Ny, T0=float(spec["T_amb_K"]),
                      omega=sor_omega, tol=sor_tol, maxit=sor_maxit)
    dT = float(T.max() - float(spec["T_amb_K"]))

    # 能量守恒审计（稳态）：注入功率 = 底面流出 + 侧面流出。两项都按 FVM 的
    # 边界导纳算，收敛解应精确成立（相对误差 ~1e-12），可作为求解器自检。
    Tamb = float(spec["T_amb_K"])
    gb = 2.0 * ky[-1] * dx / dy[-1]
    p_inj = float(q.sum())
    p_bottom = float(np.sum(gb * (T[:, -1] - Tamb)))
    p_side = 0.0
    if dirichlet_x:
        glv = 2.0 * gx
        p_side = float(np.sum(glv[contact_rows] * (T[0, contact_rows] - Tamb))
                       + np.sum(glv[contact_rows] * (T[-1, contact_rows] - Tamb)))
    resid = p_inj - p_bottom - p_side

    # 1D 串联网络（沟道取整厚）用于对照；氧化层是各向异性（如 hBN）时必须用
    # 跨面 ky，不能用面内 kx —— 否则 R_th 被低估 ~200 倍。
    r_pp = (ch["thickness_m"] / ch["ky_W_mK"] + 1.0 / props["tbc_W_m2K"]
            + ox["thickness_m"] / ox["ky_W_mK"])
    dT_1d = float(spec["Q_joule_W_m2"]) * r_pp
    # 精确 1D 解析（侧壁绝热、沟道内均匀体热源）：沟道项是 t_ch/(2 ky)。
    r_1d_dist = (ch["thickness_m"] / (2.0 * ch["ky_W_mK"])
                 + 1.0 / props["tbc_W_m2K"]
                 + ox["thickness_m"] / ox["ky_W_mK"])
    dT_1d_dist = float(spec["Q_joule_W_m2"]) * r_1d_dist

    out = {
        "SOLVE_DONE": True,
        "method": "FVM red-black SOR (pure numpy)",
        "Nx": Nx, "Ny": Ny, "unknowns": Nx * Ny,
        "sor_iterations": int(it), "sor_last_delta": float(dmax),
        "T_amb_K": float(spec["T_amb_K"]),
        "dT_peak_K": dT, "dT_peak_1D_network_K": dT_1d,
        "dT_peak_1D_distributed_K": dT_1d_dist,
        "dT_ratio_FVM_over_1D": (dT / dT_1d if dT_1d else None),
        # R_th_eff 历史上指的是解析 1D 串联网络阻值（不是 2D 解的 dT/Q），
        # 这里显式给出 2D 有效值并保留旧键，避免下游读不到。
        "R_th_eff_m2K_per_W": r_pp, "G_eff_W_m2K": 1.0 / r_pp,
        "R_pp_1D_network_m2K_per_W": r_pp,
        "R_th_2D_m2K_per_W": (dT / float(spec["Q_joule_W_m2"])
                              if float(spec["Q_joule_W_m2"]) else None),
        "G_eff_2D_W_m2K": (float(spec["Q_joule_W_m2"]) / dT if dT else None),
        "power_uW": float(spec["Q_joule_W_m2"]) * L * float(spec["W_device_m"]) * 1e6,
        "contact_bc": spec["contact_bc"], "dirichlet_x": dirichlet_x,
        "P_injected_W_per_m": p_inj, "P_bottom_W_per_m": p_bottom,
        "P_side_W_per_m": p_side,
        "energy_balance_rel": (abs(resid) / p_inj if p_inj else 0.0),
        "layer_names": [l["name"] for l in layers],
        "layer_thickness_nm": [l["t"] * 1e9 for l in layers],
    }
    jm = Nx // 2
    out["x_nm"] = [float((i + 0.5) * dx) * 1e9 for i in range(Nx)]
    out["contact_rows"] = [int(v) for v in np.nonzero(contact_rows)[0]]
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
        unc = 0                    # 内层 SOR 撞到 maxit 仍未收敛的步数
        for _ in range(2000):
            de = diag + cap2 / dt
            be = b + cap2 / dt * Tcur
            Tcur, _i, _d = sor(gx, gy, de, be, Nx, Ny, T_init=Tcur,
                               T0=float(spec["T_amb_K"]), omega=1.0,
                               tol=1e-7, maxit=20000)
            if _i >= 20000 and _d >= 1e-7:
                unc += 1
            t += dt
            ts.append(t)
            dts.append(float(Tcur.max() - float(spec["T_amb_K"])))
            dt *= ratio
            if t > 3.0 * tau_rc:
                break
        dtt = np.array(dts)
        tt = np.array(ts)
        # dT 接近 0（如 Q=0）时 tau_63 无意义，给 None 而不是编一个数值。
        idx63 = int(np.argmax(dtt >= 0.632 * dtt[-1])) if dtt[-1] > 1e-6 else None
        out["tau_fast_s"] = tau_fast
        out["tau_RC_s"] = tau_rc
        out["transient_dt0_s"] = tau_fast / 100.0
        out["transient_dt_ratio"] = ratio
        out["transient_nsteps"] = int(len(ts) - 1)
        out["transient_n_unconverged"] = int(unc)
        out["transient_sor_tol"] = 1e-7
        out["transient_sor_maxit"] = 20000
        out["tau_63_s"] = float(tt[idx63]) if idx63 is not None else None
        out["heating_t_s"] = tt.tolist()
        out["heating_dT_K"] = dtt.tolist()
        out["transient_dT_final_K"] = float(dtt[-1])
        if dtt[-1] > 0:
            out["transient_final_vs_steady_rel"] = abs(dtt[-1] - dT) / dT
        if unc:
            print("[warn] 瞬态有 %d 个时间步在内层 SOR 20000 次迭代上限时仍未收敛；"
                  "dT 可能略偏低，可提高 SOR_MAXIT 或放宽内层 tol。" % unc)
    return out, T, dy, lay


# --------------------------------------------------- 热物性组装（多技能共用）
def layer_from_db(db, name, thickness_override=None):
    """按材料名从内置库取一层（channel / oxide / substrate ...）。"""
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


def resolve_thermal_props(cwd, ch_name, ox_name, T_amb,
                          up_skill="kl-dft-cpu", up_step="step6_kappa",
                          kappa_source="auto", thickness_override=None,
                          tbc_override=None, oxide_thickness=None):
    """上游 kl 实测 + 内置材料库兜底 -> (ch, ox, tbc, tbc_ref, upstream_info, src)。

    channel 面内 kx 走 read_kappa（已按 2D 口径取物理值 kappa_2d_normalized_*）；
    跨面 kx 由 read_kappa 判定有效性，无效则保留材料库值；厚度走 read_thickness。
    device-thermal 与 device-zt 共用本函数，避免两处各写一份导致口径漂移。
    """
    db = load_materials()
    ch = layer_from_db(db, ch_name)
    ox = layer_from_db(db, ox_name)
    src = str(kappa_source or "auto").lower()
    if src == "dft":
        src = "upstream"
    if src not in ("auto", "upstream", "literature"):
        raise SystemExit("[ERROR] KAPPA_SOURCE 只能是 auto/upstream/literature，收到 %r"
                         % kappa_source)
    upstream_info = None
    # 与当前 kx 配套的参考厚度：CHANNEL_THICKNESS 换算 κ·t 时要用它。
    kx_ref_t = ch.get("thickness_m")
    if src in ("auto", "upstream"):
        up = find_upstream(cwd, up_skill, up_step)
        k = read_kappa(up, float(T_amb)) if up else None
        if k:
            ch["kx_W_mK"] = k["inplane"]
            if k["cross"]:
                ch["ky_W_mK"] = k["cross"]
            ch["source"] = "upstream:" + up
            ch["kappa_key"] = k.get("kappa_key")
            ch["kappa_note"] = k.get("note")
            ch["ky_source"] = ("upstream 跨面张量" if k["cross"]
                               else "材料库值（%s）" % (k.get("cross_note") or "上游无有效跨面值"))
            th, thp = read_thickness(cwd, up_skill)
            if th:
                ch["thickness_m"] = th
                ch["thickness_source"] = thp
                kx_ref_t = th
            else:
                kx_ref_t = None      # 上游 kappa 没有配套厚度，无法做 κ·t 换算
            upstream_info = {"file": up, "T_K": k["T"], "tensor": k["tensor"],
                             "kappa_key": k.get("kappa_key"), "note": k.get("note")}
        else:
            reason = ("上游 %s 存在但读不出 kappa（缺顶层 temperatures/kappa_xx_yy_zz）" % up
                      if up else "未找到上游 %s/%s 产物" % (up_skill, up_step))
            if src == "upstream":
                # upstream/dft 是"必须用实测"的强口径，缺数据要炸掉而不是静默用库值，
                # 否则下游把库值当 DFT 结果引用（2026-09 用户 review）。
                raise SystemExit("[ERROR] KAPPA_SOURCE=upstream，但%s；"
                                 "请先跑上游，或改 KAPPA_SOURCE=auto/literature" % reason)
            print("[warn] KAPPA_SOURCE=auto，但%s，改用材料库" % reason)
    if thickness_override is not None:
        t_new = float(thickness_override)
        if t_new <= 0:
            raise SystemExit("[ERROR] CHANNEL_THICKNESS 必须 >0，收到 %r" % thickness_override)
        if (str(ch.get("kind") or "").lower() == "2d" and kx_ref_t
                and abs(t_new - float(kx_ref_t)) > 0.0):
            factor = float(kx_ref_t) / t_new
            ch["kx_W_mK"] = ch["kx_W_mK"] * factor
            ch["kx_scale_note"] = ("2D 口径：面内 kappa 从 t=%.4g nm 换算到 %.4g nm"
                                   "（kx *= %.5g，保持 κ·t 不变）"
                                   % (float(kx_ref_t) * 1e9, t_new * 1e9, factor))
        ch["thickness_m"] = t_new
        ch["thickness_source"] = "step.conf CHANNEL_THICKNESS"
    if ch.get("thickness_m") is None:
        raise SystemExit("[ERROR] 沟道 %s 厚度未知：请设 CHANNEL_THICKNESS，"
                         "或让上游 kl 产出 thickness_2d.json" % ch_name)
    if oxide_thickness is not None:
        # 显式 OXIDE_THICKNESS 一律覆盖材料库厚度（含 hBN 的 0.34 nm 层厚）；
        # 旧实现只在库厚度为空时才用，导致设了 OXIDE_THICKNESS 也不生效。
        t_ox = float(oxide_thickness)
        if t_ox <= 0:
            raise SystemExit("[ERROR] OXIDE_THICKNESS 必须 >0，收到 %r" % oxide_thickness)
        ox["thickness_m"] = t_ox
        ox["thickness_source"] = "step.conf OXIDE_THICKNESS"
    if ox.get("thickness_m") is None:
        raise SystemExit("[ERROR] 氧化层 %s 厚度未知：请设 OXIDE_THICKNESS" % ox_name)
    tbc, tbc_ref = resolve_tbc(db, ch_name, ox_name)
    if tbc_override:
        tbc = float(tbc_override)
        tbc_ref = "step.conf TBC_OVERRIDE (e.g. GPUMD/MD or experiment)"
    return ch, ox, tbc, tbc_ref, upstream_info, src
