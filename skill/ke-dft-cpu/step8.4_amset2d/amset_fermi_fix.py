"""
amset_fermi_fix.py —— AMSET 0.5.1 求费米能级的稳健化（运行时插件，不改 AMSET 源码；V132）

用法（S8 / S8.4 的 gen 已把本文件复制进运行目录，并在 python -c 里 import；amset2d_plugin 与
postprocess_intrinsic 也会导入它）：
    python -c "import amset_fermi_fix; ...; Runner.from_directory('.').run()"

开关：环境变量 AZ_FERMI_FIX。未设 = 开；0/false/off/no = 关（只用于复现旧结果做对照）。
    NaN 守卫（V131）不受开关影响，始终生效。

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
MARK = "_az_fermi_fix"
STATE = {"applied": False, "fallbacks": [], "refused": []}


def enabled():
    """未设 = 开；只有显式 0/false/off/no 才关。"""
    return os.environ.get(ENV_NAME, "").strip().lower() not in _OFF


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


if apply() and not _quiet():
    print("[fermi_fix] 已挂上：AMSET 求不到费米能级时改用二分法、宽松解偏差 > 1%% 报错、NaN 报错%s"
          % ("" if enabled() else "（%s=%s -> 只留 NaN 守卫）" % (ENV_NAME, os.environ.get(ENV_NAME))),
          file=sys.stderr, flush=True)
