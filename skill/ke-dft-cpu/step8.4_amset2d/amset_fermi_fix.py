"""
amset_fermi_fix.py —— AMSET 0.5.1 求费米能级的稳健化（运行时插件，不改 AMSET 源码；V132）

用法（S8 / S8.4 的 gen 已把本文件复制进运行目录，并在 python -c 里 import；amset2d_plugin 与
postprocess_intrinsic 也会导入它）：
    python -c "import amset_fermi_fix; ...; Runner.from_directory('.').run()"

开关：环境变量 AZ_FERMI_FIX。未设 = 开；0/false/off/no = 关（只用于复现旧结果做对照）。
    NaN 守卫（V131）与散射率插值死片保护（[0043]）不受开关影响，始终生效。

------------------------------------------------------------------------------
问题（V132，WS2 n 型 100 K 报 "Could not find fermi within 100.0% of concentration=-2.759e-07"）：
FermiDos.get_fermi 是贪心网格搜索 —— 从本征 E_F 出发，步长 0.1 Ha（2.72 eV）、0.272 eV、27 meV …
每层取 |doping/conc − 1| 最小的格点，再在它周围 ±50 步细分。可 doping = (nelect − Σ f·DOS·dE)/V
是两个 ≈ nelect 的数相减，双精度下每胞少于 ~nelect·1e-16 的载流子就是 0。低温时第二层（0.272 eV）
离目标最近、又没越过目标的格点在带边下 ~0.25 eV，n/目标 ~ exp(−0.2 eV/kT) ≈ 1e-10（100 K），
已被舍入吃掉：整个带隙的误差都"等于 1"，argmin 取第一个（最靠价带的），下一层的窗口（±1.36 eV）
够不到导带边，最后误差 100%。中不中招取决于"本征 E_F + j·0.272 eV"与带边怎么对齐，所以同配置的
MoS₂/WSe2 没事、WS2 n 型 100 K 失败；300 K 时 exp(−0.2/kT) ≈ 5e-4，看得见，不出事。3D 是同一段代码。
fermi_probe 实测（WS2）：DOS 正常（带隙内 0 态、本征 E_F 有限），失败的两行贪心误差 100%，二分法照样解。

更隐蔽的一面：AmsetData.set_doping_and_temperatures 按 tol = 1e-5, 1e-4, …, 1 依次试，前面的失败了
就接受后面的 —— 最后一级接受的是**掺杂偏差至多 100%** 的 E_F，不报错、不记日志。

修法：
  · 先跑 AMSET 原搜索；它在请求的 tol 下抛 ValueError 时，对 doping(E_F) − conc 做 brentq
    （doping 随 E_F 单调减；区间取 DOS 能量范围外扩 1 eV，浓度在 DOS 能容纳的范围内就有唯一根），
    达到 tol 就按原格式返回（含 return_electron_hole_conc）。阶梯第一级 1e-5 就满足，不再往下退。
    AMSET 原搜索能解的情况：结果逐位不变。
  · 宽松容差的上限：tol > 1% 时 AMSET 原搜索给出的解若实际偏差 > 1%（= postprocess 的
    carrier_guard 口径），不接受、报错 —— 走到这一步说明二分法也无解（浓度超出 DOS 能容纳的范围）。
  · NaN 守卫（V131 从 amset2d_plugin 移来，3D 也生效）：原搜索返回 NaN（中间量 NaN 时
    "NaN > tol" 为假，不报错）-> 抛 ValueError。
护栏：get_fermi 的签名里要有 tol / return_electron_hole_conc，对不上就不打补丁（stderr 说明）。

[0040] 本征点（doping = 0）：AmsetData 对 doping=0 走 get_fermi_from_num_electrons —— 按"总电子数相对误差"
贪心搜索。有带隙时整个带隙里的电子数都等于 nelect（本征载流子每胞 ~1e-17 个，远低于 1e-8 的容差），
argmin 取到网格上第一个带隙内的点（靠价带一侧），E_F 位置任意，S/σ 差几个数量级；electron_conc/hole_conc
也不填（= 0），迁移率除以 0。修法：DOS 在本征 E_F 处有带隙时，改解 n(E_F) = p(E_F)（log-sum-exp，
不受相消/下溢影响），并把 n、p 填进 electron_conc / hole_conc。无带隙（金属）照旧走 AMSET 原搜索。

[0043] 本征点 + IMP 的 "SVD did not converge"：不是本征浓度下溢（n_i ∝ exp(−E_g/2kT) 本来就小，是物理值）。
AMSET IMP 的 N_ii = |n − p|/Z，本征点 = 0 -> IMP 率整片为 0 -> _interpolate_zero_rates 的"任一掺杂/温度
≤ 1e6 就标记"把该机制全部 k 点标成待插值、又无处可插，其它掺杂里 MRTA 留下的负率不再被替掉 ->
fill_rates_outside_cutoffs 取 log 得 NaN -> pinv(L11) 的 SVD 不收敛。修法：本征点 n = p = √(np) 严格相等
（N_ii 精确为 0），插值时把这种死片单独处理（见 _wrap_interp）。
------------------------------------------------------------------------------
"""

import functools
import inspect
import os
import sys

import numpy as np

ENV_NAME = "AZ_FERMI_FIX"
_OFF = ("0", "false", "off", "no")
HA_EV = 27.211386245988
MAX_REL_ERR = 0.01          # 与 postprocess_intrinsic.check_achieved_doping 的 1% 口径一致
# [0055] 小掺杂（每胞净载流子 |d|·V < 1e-7）：AMSET 的 get_doping = (nelect − Σ f·DOS·dE)/V、hole_conc = (nelect − Σ_VB)/V
#   都是两个 ≈ nelect 的数相减，|d|·V 接近 nelect·1e-16 时 E_F 与空穴浓度都失真（3D 1e10 cm^-3 ≈ 每胞 4e-13）。
#   低于这个阈值改走 small_doping_fermi（n、p 分开 log-sum-exp）。出厂网格最小点（3D 1e17、2D 1e11 cm^-2）
#   每胞 ≥ ~1e-6，不受影响、结果逐位不变。
SMALL_DOPING_PER_CELL = 1e-7
KB_AU = 3.166811563455e-06  # Ha/K（= amset.constants.boltzmann_au；单测在没装 amset 时用）
MARK = "_az_fermi_fix"
STATE = {"applied": False, "fallbacks": [], "refused": []}


def enabled():
    """未设 = 开；只有显式 0/false/off/no 才关。"""
    return os.environ.get(ENV_NAME, "").strip().lower() not in _OFF


def _kb_au():
    try:
        from amset.constants import boltzmann_au
        return float(boltzmann_au)
    except Exception:                                              # noqa: BLE001
        return KB_AU


def solve_fermi(dos, concentration, temperature, pad_ev=1.0, xtol_ha=1e-13):
    """doping(E_F) = concentration 的 brentq 解（Ha）。doping 随 E_F 单调减（N(E_F) 单调增）。
    区间端点上不变号或不是有限数 -> None。"""
    from scipy.optimize import brentq
    E = np.asarray(dos.energies, dtype=float)
    lo = float(E.min()) - pad_ev / HA_EV
    hi = float(E.max()) + pad_ev / HA_EV

    def g(e):
        return dos.get_doping(e, temperature) - concentration

    with np.errstate(over="ignore", invalid="ignore"):
        glo, ghi = g(lo), g(hi)
        if not (np.isfinite(glo) and np.isfinite(ghi)) or glo * ghi > 0:
            return None
        return float(brentq(g, lo, hi, xtol=xtol_ha, maxiter=500))


def _dos_arrays(dos):
    """(E[Ha], g·dE[每胞态数], 本征 E_F[Ha])。"""
    E = np.asarray(dos.energies, dtype=float)
    ef0 = float(dos.efermi)
    if not getattr(dos, "atomic_units", True):
        E, ef0 = E / HA_EV, ef0 / HA_EV
    g = np.asarray(dos.tdos, dtype=float) * abs(float(E[1] - E[0]))
    return E, g, ef0


def band_gap_edges(dos, rel_thr=1e-10, min_gap_ev=0.01):
    """本征 E_F 两侧有态的最近能量 (VBM, CBM)（Ha）；E_F 不在带隙里（金属）或带隙 < min_gap_ev -> None。"""
    E, g, ef0 = _dos_arrays(dos)
    occ = g > rel_thr * max(float(g.max()), 1e-300)
    below, above = E[occ & (E <= ef0)], E[occ & (E > ef0)]
    if below.size == 0 or above.size == 0:
        return None
    vbm, cbm = float(below.max()), float(above.min())
    return (vbm, cbm) if (cbm - vbm) * HA_EV >= min_gap_ev else None


def log_carriers(dos, ef, temperature):
    """(ln n, ln p)：每胞导带电子数 / 价带空穴数的自然对数（log-sum-exp，不相消、不下溢）。"""
    E, g, ef0 = _dos_arrays(dos)
    kt = float(temperature) * _kb_au()
    pos = g > 0
    cb, vb = pos & (E > ef0), pos & (E <= ef0)

    def lse(v):
        m = float(v.max())
        return m + float(np.log(np.exp(v - m).sum()))

    ln_n = lse(np.log(g[cb]) - np.logaddexp(0.0, (E[cb] - ef) / kt))
    ln_p = lse(np.log(g[vb]) - np.logaddexp(0.0, (ef - E[vb]) / kt))
    return ln_n, ln_p


def intrinsic_fermi(dos, temperature, xtol_ha=1e-12):
    """本征 E_F（Ha）：解 ln n = ln p（单调增），区间 [VBM, CBM]。无带隙 -> None。"""
    from scipy.optimize import brentq
    edges = band_gap_edges(dos)
    if edges is None or not temperature > 0:
        return None

    def h(e):
        ln_n, ln_p = log_carriers(dos, e, temperature)
        return ln_n - ln_p

    lo, hi = edges
    if not (h(lo) < 0 < h(hi)):
        return None
    return float(brentq(h, lo, hi, xtol=xtol_ha, maxiter=500))


def small_doping_fermi(dos, concentration, temperature):
    """[0055] 小掺杂的 E_F：解 n − p = |d|（n 型，concentration < 0）或 p − n = d（p 型），d = |concentration|·V
    为每胞净载流子数。n、p 由 log_carriers 分开算（log-sum-exp），不经 nelect − Σ 的相减。
    写成 ln n − ln(p + d) = 0 / ln(n + d) − ln p = 0，两者都随 E_F 单调增，区间 [VBM, CBM]。
    返回 (E_F[Ha], n, p)，n、p 为 1/bohr³（与 AMSET electron_conc / hole_conc 同单位）；无带隙 / 不变号 -> None。"""
    import math
    from scipy.optimize import brentq
    edges = band_gap_edges(dos)
    if edges is None or not temperature > 0 or concentration == 0:
        return None
    vol = float(dos.structure.volume)
    ln_d = math.log(abs(float(concentration)) * vol)

    def h(e):
        ln_n, ln_p = log_carriers(dos, e, temperature)
        if concentration < 0:
            return ln_n - np.logaddexp(ln_p, ln_d)
        return np.logaddexp(ln_n, ln_d) - ln_p

    lo, hi = edges
    with np.errstate(over="ignore", invalid="ignore"):
        hlo, hhi = h(lo), h(hi)
        if not (np.isfinite(hlo) and np.isfinite(hhi)) or hlo * hhi > 0:
            return None
        # rtol 取到双精度极限：n_i ≫ |d|（窄带隙高温）时 n − p 要在 E_F 上分辨 ~kT·d/n_i
        ef = float(brentq(h, lo, hi, xtol=1e-16, rtol=4 * np.finfo(float).eps, maxiter=500))
    ln_n, ln_p = log_carriers(dos, ef, temperature)
    return ef, math.exp(ln_n) / vol, math.exp(ln_p) / vol


def rel_error(dos, ef, concentration, temperature):
    with np.errstate(over="ignore", invalid="ignore"):
        return float(abs(dos.get_doping(ef, temperature) / concentration - 1.0))


def _to_cm3(conc):
    try:
        from amset.constants import cm_to_bohr
        return conc * cm_to_bohr ** 3
    except Exception:                                              # noqa: BLE001
        return conc


def _result(dos, ef, temperature, want_eh):
    if want_eh:
        _, n_elec, n_hole = dos.get_doping(ef, temperature, return_electron_hole_conc=True)
        return ef, n_elec, n_hole
    return ef


def unwrap(fn):
    """剥掉本插件与 V131 amset2d_plugin（旧运行目录里的副本）套的壳，返回 AMSET 原 get_fermi。"""
    for _ in range(8):
        if getattr(fn, MARK, False):
            fn = fn.__wrapped__
        elif getattr(fn, "_amset2d_nan", False):          # V131：NaN 壳，*a/**k 透传，原函数在模块全局里
            fn = getattr(fn, "__globals__", {}).get("_orig_get_fermi", fn)
        else:
            break
    return fn


def _wrap(orig, sig):
    @functools.wraps(orig)
    def get_fermi(self, concentration, temperature, *args, **kwargs):
        ba = sig.bind(self, concentration, temperature, *args, **kwargs)
        ba.apply_defaults()
        tol = float(ba.arguments["tol"])
        want_eh = bool(ba.arguments["return_electron_hole_conc"])
        if (enabled() and concentration != 0 and np.isfinite(self.efermi)
                and abs(float(concentration)) * float(self.structure.volume) < SMALL_DOPING_PER_CELL):
            try:                                                   # [0055] 小掺杂：不走 nelect − Σ 的相减
                res = small_doping_fermi(self, concentration, temperature)
            except Exception:                                      # noqa: BLE001
                res = None
            if res is not None:
                ef, ne, nh = res
                key = float(concentration)
                if key not in {x["conc"] for x in STATE.setdefault("small", [])} and not _quiet():
                    print("[fermi_fix] 小掺杂 %.3e cm^-3（每胞 %.1e < %g）：n、p 分开 log-sum-exp 解 E_F"
                          "（%g K：E_F = %.4f eV，n = %.3e，p = %.3e cm^-3）"
                          % (_to_cm3(concentration), abs(key) * float(self.structure.volume),
                             SMALL_DOPING_PER_CELL, temperature, ef * HA_EV, _to_cm3(ne), _to_cm3(nh)), flush=True)
                STATE["small"].append({"conc": key, "T": float(temperature), "ef_ha": ef})
                return (ef, ne, nh) if want_eh else ef
        try:
            out = orig(self, concentration, temperature, *args, **kwargs)
        except ValueError:
            if not enabled() or concentration == 0 or not np.isfinite(self.efermi):
                raise                    # 本征 E_F 非有限时电子/空穴的划分没有意义，不兜底
            ef = solve_fermi(self, concentration, temperature)
            if ef is None:
                raise
            err = rel_error(self, ef, concentration, temperature)
            if not err <= tol:
                raise
            STATE["fallbacks"].append({"conc": float(concentration), "T": float(temperature), "tol": tol,
                                       "ef_ha": ef, "rel_err": err})
            print("[fermi_fix] 掺杂 %.3e cm^-3 @ %g K：AMSET 贪心搜索在 tol=%g 下失败 -> 二分法 E_F = %.4f eV"
                  "（相对误差 %.1e）" % (_to_cm3(concentration), temperature, tol, ef * HA_EV, err), flush=True)
            return _result(self, ef, temperature, want_eh)
        ef = out[0] if isinstance(out, tuple) else out
        if not np.isfinite(ef):
            raise ValueError("[fermi_fix] get_fermi 得到 NaN（本征费米能级 %r，浓度 %g）—— AMSET 原式会静默返回 NaN；"
                             "用 tools/fermi_probe.py 诊断" % (self.efermi, concentration))
        if enabled() and tol > MAX_REL_ERR and concentration != 0:
            err = rel_error(self, ef, concentration, temperature)
            if not err <= MAX_REL_ERR:
                STATE["refused"].append({"conc": float(concentration), "T": float(temperature), "tol": tol,
                                         "ef_ha": float(ef), "rel_err": err})
                raise ValueError("[fermi_fix] 掺杂 %.3e cm^-3 @ %g K：AMSET 只在 tol=%g 下找到 E_F = %.4f eV，"
                                 "实际掺杂偏差 %.1f%% > %.0f%%，不接受（二分法也无解：浓度超出 DOS 能容纳的范围？"
                                 "用 tools/fermi_probe.py 看）"
                                 % (_to_cm3(concentration), temperature, tol, ef * HA_EV, err * 100,
                                    MAX_REL_ERR * 100))
        return out

    setattr(get_fermi, MARK, True)
    return get_fermi


def _wrap_intrinsic(orig):
    @functools.wraps(orig)
    def get_fermi_from_num_electrons(self, num_electrons, temperature, *args, **kwargs):
        if enabled():
            try:
                ef = intrinsic_fermi(self, temperature)
            except Exception:                                      # noqa: BLE001
                ef = None
            if ef is not None:
                STATE.setdefault("intrinsic", []).append({"T": float(temperature), "ef_ha": ef})
                return ef
        return orig(self, num_electrons, temperature, *args, **kwargs)

    setattr(get_fermi_from_num_electrons, MARK, True)
    return get_fermi_from_num_electrons


UNDOPED_NI_JSON = "undoped_ni.json"   # [0052] 本征点 n_i(T) 记录（S8.3 的 intrinsic_carrier/ 图表标注 n_i 用）
#   [0054] 另记带边 E_v / E_c（与 E_F 同一能量参考、含 scissor），S8.3 据此算 η = (E_c − E_F)/k_BT，判本征点是否仍是 f ≈ 0


def write_undoped_ni(rec, run_dir="."):
    """[0052] 把本征点 (T, E_F, n_i) 写进运行目录的 undoped_ni.json —— transport.json 里没有载流子浓度
    字段，S8.3 / ZT 要标注本征浓度只能靠这份记录（旧结果退回解析 queue.out 里的 [fermi_fix] 本征点 行）。
    只在 AMSET 运行目录（有 settings.yaml）里写，单测/探针在别处调用时不落盘。n_i 是每胞体积口径 cm^-3
    （2D = 含真空的元胞体积，换算面密度/文献口径由 S8.3 按 c、t 做）。"""
    import json
    if not os.path.isfile(os.path.join(run_dir, "settings.yaml")):
        return None
    rec = sorted(dict((r[0], r) for r in rec).values())
    path = os.path.join(run_dir, UNDOPED_NI_JSON)
    try:
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump({"schema": "autozt-undoped-ni/1",
                       "note": "AMSET doping = 0 点：E_F 由 n = p 求得，n = p = n_i（每胞体积口径 cm^-3）",
                       "temperatures_K": [r[0] for r in rec], "E_F_eV": [r[1] for r in rec],
                       "n_i_cm-3": [r[2] for r in rec],
                       "E_v_eV": [r[3] if len(r) > 4 else None for r in rec],
                       "E_c_eV": [r[4] if len(r) > 4 else None for r in rec]}, fh, indent=1)
        os.replace(tmp, path)
    except OSError as e:
        print("[fermi_fix] 写 %s 失败（%s）—— 本征浓度只在日志里" % (UNDOPED_NI_JSON, e), file=sys.stderr)
        return None
    return path


def _wrap_set_doping(orig):
    @functools.wraps(orig)
    def set_doping_and_temperatures(self, doping, temperatures):
        out = orig(self, doping, temperatures)
        if enabled() and getattr(self, "dos", None) is not None and band_gap_edges(self.dos) is not None:
            vol = float(self.dos.structure.volume)
            vbm, cbm = band_gap_edges(self.dos)          # [0054] 带边（Ha），写进 undoped_ni.json
            rec = []
            for n, t in np.ndindex(self.fermi_levels.shape):
                if self.doping[n] == 0:
                    ln_n, ln_p = log_carriers(self.dos, self.fermi_levels[n, t], self.temperatures[t])
                    # [0043] n = p = n_i = √(np)（质量作用定律）严格相等：AMSET IMP 的 N_ii = |n − p|/Z 才是 0，
                    # 而不是 E_F 求解容差留下的残差（那点残差乘 q→0 处 1/(q_TF)⁴ 量级的因子 = 正负垃圾率）。
                    n_i = np.exp(0.5 * (ln_n + ln_p)) / vol
                    self.electron_conc[n, t] = n_i
                    self.hole_conc[n, t] = n_i
                    rec.append((float(self.temperatures[t]), float(self.fermi_levels[n, t] * HA_EV),
                                float(_to_cm3(self.electron_conc[n, t])), vbm * HA_EV, cbm * HA_EV))
                    if not _quiet():
                        print("[fermi_fix] 本征点 %g K：E_F = %.4f eV，n = p = %.3e cm^-3（n=p 求解，log-sum-exp）"
                              % (self.temperatures[t], self.fermi_levels[n, t] * HA_EV,
                                 _to_cm3(self.electron_conc[n, t])), flush=True)
            if rec and not _quiet():
                write_undoped_ni(rec)
        return out

    setattr(set_doping_and_temperatures, MARK, True)
    return set_doping_and_temperatures


RATE_FLAG = 1e6            # = amset.scattering.calculate._interpolate_zero_rates 的"非零"阈值


def _dead_slices(rates, masks):
    """[0043] 某机制在某 (掺杂, 温度, 能带) 窗口内没有一个 > 1e6 的率（前因子为 0：本征点的 IMP，N_ii = 0）。"""
    dead = []
    for spin, r in rates.items():
        for s, d, t, b in np.ndindex(r.shape[:-1]):
            win = np.ones(r.shape[-1], bool) if masks is None else ~masks[spin][s, d, t, b]
            if win.any() and not np.any(r[s, d, t, b, win] > RATE_FLAG):
                dead.append((spin, (s, d, t, b), win))
    return dead


def _wrap_interp(orig):
    """[0043] AMSET 的 _interpolate_zero_rates 按"任一掺杂/温度 ≤ 1e6"标记要插值的 k 点：只要某个
    (掺杂, 温度) 的率整片为 0（本征点 IMP），这个机制所有 k 点都被标记 → "全零、无法插值" 分支只加 1e-32，
    **其它掺杂**里 MRTA 因子留下的负率（正常情况下会被最近邻插值替掉）原样保留 →
    fill_rates_outside_cutoffs 对负数取 log = NaN → 寿命 NaN → BoltzTraP2 pinv(L11) "SVD did not converge"。
    修法：死片先垫占位值不参与标记，原函数照常插值其余各片；死片再还原，≤ 0 / 非有限的值换成 small_val
    （= AMSET 对前因子为 0 的机制本来的处理）。没有死片时原样调用，结果逐位不变。"""
    @functools.wraps(orig)
    def _interpolate_zero_rates(rates, kpoints, masks=None, *args, **kwargs):
        dead = _dead_slices(rates, masks)
        if not dead:
            return orig(rates, kpoints, masks, *args, **kwargs)
        try:
            from amset.constants import small_val
        except Exception:                                          # noqa: BLE001
            small_val = 1e-32
        saved = []
        for spin, idx, win in dead:
            row = rates[spin][idx]
            saved.append(row[win].copy())
            row[win] = 1e14
        out = orig(rates, kpoints, masks, *args, **kwargs)
        for (spin, idx, win), val in zip(dead, saved):
            out[spin][idx][win] = np.where(np.isfinite(val) & (val > 0), val, small_val)
        if not _quiet():
            labs = sorted({"%d/%d/%d" % idx[:3] for _, idx, _ in dead})
            print("[fermi_fix] 散射率插值：%d 个 (机制/掺杂/温度) 片窗口内全 ≤ 1e6（前因子 0，如本征点 IMP），"
                  "单独处理、不再拖垮其它掺杂的插值：%s" % (len(labs), ", ".join(labs[:12])), flush=True)
        return out

    setattr(_interpolate_zero_rates, MARK, True)
    return _interpolate_zero_rates


def apply_rate_guard():
    """[0043] 散射率插值的死片保护（幂等；与 NaN 守卫一样不受 AZ_FERMI_FIX 开关影响：
    没有死片时逐位不变，有死片时原版必然崩）。"""
    try:
        import amset.scattering.calculate as C
    except Exception as e:                                         # noqa: BLE001
        print("[fermi_fix] 导入 amset.scattering.calculate 失败（%s）—— 散射率插值不打补丁" % e, file=sys.stderr)
        return False
    cur = getattr(C, "_interpolate_zero_rates", None)
    if cur is None:
        print("[fermi_fix] AMSET 没有 _interpolate_zero_rates（版本不同）—— 散射率插值不打补丁", file=sys.stderr)
        return False
    if not getattr(cur, MARK, False):
        C._interpolate_zero_rates = _wrap_interp(cur)
    return True


def apply_intrinsic():
    """[0040] 本征点（doping = 0）的 E_F 与 n、p（幂等）。"""
    try:
        from amset.electronic_structure.dos import FermiDos
        from amset.core.data import AmsetData
    except Exception as e:                                         # noqa: BLE001
        print("[fermi_fix] 导入 amset 失败（%s）—— 本征点不打补丁" % e, file=sys.stderr)
        return False
    if not getattr(FermiDos.get_fermi_from_num_electrons, MARK, False):
        FermiDos.get_fermi_from_num_electrons = _wrap_intrinsic(FermiDos.get_fermi_from_num_electrons)
    if not getattr(AmsetData.set_doping_and_temperatures, MARK, False):
        AmsetData.set_doping_and_temperatures = _wrap_set_doping(AmsetData.set_doping_and_temperatures)
    return True


def apply():
    """给 FermiDos.get_fermi 套上外壳（幂等）。返回是否已套上。"""
    try:
        from amset.electronic_structure.dos import FermiDos
    except Exception as e:                                         # noqa: BLE001
        print("[fermi_fix] 导入 amset 失败（%s）—— 不打补丁" % e, file=sys.stderr)
        return False
    cur = FermiDos.get_fermi
    if getattr(cur, MARK, False):
        STATE["applied"] = True
        return True
    sig = inspect.signature(unwrap(cur))
    if not {"tol", "return_electron_hole_conc"} <= set(sig.parameters):
        print("[fermi_fix] FermiDos.get_fermi 的签名与 AMSET 0.5.1 不同（%s）—— 不打补丁"
              % list(sig.parameters), file=sys.stderr)
        return False
    FermiDos.get_fermi = _wrap(cur, sig)
    STATE["applied"] = True
    return True


def original_get_fermi():
    """未套壳的 FermiDos.get_fermi（fermi_probe 用它复现 AMSET 原搜索）。"""
    from amset.electronic_structure.dos import FermiDos
    return unwrap(FermiDos.get_fermi)


def _quiet():
    try:
        import multiprocessing
        # spawn 的散射子进程不刷屏（同 amset_ir_fix._quiet）。
        return (multiprocessing.parent_process() is not None
                or bool(getattr(multiprocessing.current_process(), "_inheriting", False)))
    except Exception:                                              # noqa: BLE001
        return False


apply_intrinsic()
apply_rate_guard()
if apply() and not _quiet():
    print("[fermi_fix] 已挂上：AMSET 求不到费米能级时改用二分法、宽松解偏差 > 1%% 报错、NaN 报错%s"
          % ("" if enabled() else "（%s=%s -> 只留 NaN 守卫）" % (ENV_NAME, os.environ.get(ENV_NAME))),
          file=sys.stderr, flush=True)
