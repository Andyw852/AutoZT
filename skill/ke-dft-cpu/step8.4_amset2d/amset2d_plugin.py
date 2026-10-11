"""
amset2d_plugin.py —— AMSET 的二维散射核（运行时插件，不改 AMSET 源码）

用法（在 AMSET 运行目录内，目录里必须有 2d_correction.json）：
    python -c "import amset2d_plugin; from amset.log import initialize_amset_logger as L; L(); \
               from amset.core.run import Runner; Runner.from_directory('.').run()"

本文件由 gen_step14_amset2d.py 复制到 step8.4_amset2d/ 运行目录，只在那一次运行里生效，
不改 AMSET 安装，也不影响三维项目（三维走原版 amset run）。

运行目录需要 2d_correction.json，至少含：
    cell_c_A            超胞 c 长度（Å）
    layer_normal        层法向（必须 ∥ 笛卡尔 z）
    eps_inf_slab        VASP DFPT 原始 slab 高频介电张量 3×3（不扣真空）
    eps_static_slab     VASP DFPT 原始 slab 静态介电张量 3×3（不扣真空）
可选：eps_env（默认 1.0）、imp_distance_A（默认 0.0）、pop_lo_dispersion（默认 true）、
      mechanisms_2d（默认 ["ADP","POP","IMP","PIE"]）
settings.yaml 里的 elastic_constant / piezoelectric_constant 必须是【原始 slab 值、标准 Voigt 顺序】，
即不能乘 c/t，且已由 gen 从 VASP 的 (XX YY ZZ XY YZ ZX) 重排成 (XX YY ZZ YZ XZ XY)。

核心映射（推导见 skill/ke-dft-cpu/METHODOLOGY.md 第 4 节）：
    AMSET 的积分是三维的，真空方向的 q_z 积分贡献因子 2π/c，而 slab 的能带沿 k_z 平；
    因此把二维耦合 G_2D 放进 AMSET 的 prefactor×factor 时要乘 c：G_s = c·G_2D。
    对 ADP 的等价说法：弹性常数直接用**原始 slab 值**（= C_2D/c），不能再乘 c/t。
"""

import json
import os
import sys

import numpy as np

import amset
import amset.scattering.calculate as C
from amset.constants import (angstrom_to_bohr as A2B, boltzmann_au,
                            coulomb_to_au, m_to_bohr)
from amset.electronic_structure.fd import fd
from amset.scattering.inelastic import PolarOpticalScattering
from amset.scattering.basic import AbstractBasicScattering as _ABS_BASIC
from amset.scattering.elastic import (
    AcousticDeformationPotentialScattering as _ADP,
    IonizedImpurityScattering as _IMP,
    PiezoelectricScattering as _PIE,
    get_christoffel_tensors,
)

# V127：spawn 的散射子进程只导入本插件（worker 入口在这里），IR_FIX 的 ADP 星平均
# 也必须在子进程里生效 —— 运行目录里有 amset_ir_fix.py 就一并导入（开关仍看 AZ_IR_FIX）。
try:
    import amset_ir_fix  # noqa: F401
except ImportError:
    pass

# ---------------- 版本与接口护栏 ----------------
# 集群 amset_clean 环境是 0.4.19；本地/新版是 0.5.1。两者的 scattering/elastic.py、
# inelastic.py、common.py 逐字节相同，calculate.py 只差重叠因子的 umklapp 处理，
# 所以同一套替换注册表在两边都成立。这里只对**用到的接口**做能力校验，
# 不做"必须等于某个版本号"的硬断言（旧版硬断言 0.5.1 会在集群上直接拒跑）。
_SUPPORTED_AMSET = ("0.4.19", "0.5.1")
if amset.__version__ not in _SUPPORTED_AMSET:
    print("[amset2d][WARN] 未在 %s 上核对过本插件（已核对：%s），"
          "将继续运行但请自行核对散射接口。"
          % (amset.__version__, "/".join(_SUPPORTED_AMSET)), file=sys.stderr)
for _n in ("_scattering_mechanisms", "scattering_worker",
           "AcousticDeformationPotentialScattering"):
    assert hasattr(C, _n), f"amset 内部结构已变化：amset.scattering.calculate 缺少 {_n}"
assert hasattr(C, "_get_cutoff_pad") or True     # _get_cutoff_pad 在 core.run 里，仅提示

# ---------------- 读取 2D 参数 ----------------
_REC_PATH = os.environ.get("AMSET2D_RECORD", "2d_correction.json")
if not os.path.isfile(_REC_PATH):
    raise ValueError("amset2d：运行目录里没有 %s —— 2D 插件不能在三维运行里使用"
                     "（请用 gen_step14_amset2d.py 生成的 step8.4_amset2d 目录）" % _REC_PATH)
with open(_REC_PATH) as _f:
    _REC = json.load(_f)

_N = np.asarray(_REC["layer_normal"], float)
_N = _N / np.linalg.norm(_N)
if abs(_N[2]) < 0.999:
    raise ValueError("amset2d：层法向必须平行于笛卡尔 z（请把 a、b 转到 xy 平面后再算），"
                     "当前 layer_normal=%s" % (_REC["layer_normal"],))
_C = float(_REC["cell_c_A"]) * A2B                    # bohr
_EPS_ENV = float(_REC.get("eps_env", 1.0) or 1.0)
_D_IMP = float(_REC.get("imp_distance_A", 0.0) or 0.0) * A2B
_LO_DISP = bool(_REC.get("pop_lo_dispersion", True))
_MECHS = set(_REC.get("mechanisms_2d", ["ADP", "POP", "IMP", "PIE"]))


def _r_tensor(key):
    """二维极化长度张量 r = c(ε_slab,∥ − 1)/2（2×2，bohr），真空收敛后与真空长度无关。

    输入的 ε_slab 是 VASP DFPT 原始 slab 张量（含真空稀释）的完整 3×3，
    这里只取面内 2×2（斜方/单斜需要非对角分量，所以不能只取对角元）。
    """
    if key not in _REC:
        raise KeyError("amset2d：2d_correction.json 缺少 %s（原始 slab 介电张量）" % key)
    e = np.asarray(_REC[key], float)[:2, :2]
    return _C * (e - np.eye(2)) / 2.0


def _inplane(unit_q, norm_q_sq):
    """把三维 q 投影到面内：返回 q∥（n×2, bohr⁻¹）、|q∥|、面内单位向量（n×3, z=0）。"""
    q = np.sqrt(norm_q_sq)[:, None] * unit_q
    qv = q[:, :2]
    qn = np.linalg.norm(qv, axis=1)
    u = np.zeros_like(unit_q)
    u[:, :2] = qv / np.maximum(qn, 1e-12)[:, None]
    return qv, qn, u


def _q_eps(qv, qn, r):
    """q·ε(q) = ε_env·|q| + q·r·q（张量 Keldysh 介电函数，q→0 无奇点）。"""
    return _EPS_ENV * qn + np.einsum("ni,ij,nj->n", qv, r, qv)


def _qtf(amset_data):
    """二维 Thomas–Fermi 波矢 q_TF = 2π·c·∂n/∂μ，形状 (ndop, ntemp)。

    ∂n/∂μ 用 AMSET 自己的总态密度积分（与 amset.scattering.common 同口径），
    再乘超胞 c 换成**面**密度；2π 来自二维屏蔽库仑（高斯单位，e²=1）。
    与 AMSET 三维 β² 的关系：q_TF = (c·ε_avg/2)·β²_3D。
    """
    trap = getattr(np, "trapezoid", None) or np.trapz
    tdos, en = amset_data.dos.tdos, amset_data.dos.energies
    vol = amset_data.structure.volume
    out = np.zeros(amset_data.fermi_levels.shape)
    for n, t in np.ndindex(out.shape):
        kt = amset_data.temperatures[t] * boltzmann_au
        f = fd(en, amset_data.fermi_levels[n, t], kt)
        out[n, t] = 2 * np.pi * _C * trap(tdos * f * (1 - f), x=en) / (kt * vol)
    return out


def _inplane_modes(elastic, u):
    """面内 2×2 Christoffel 子块：返回刚度 w (n,2) 与三维极化向量 v (n,2,3)。

    u 的 z 分量为 0，所以完整 3×3 Christoffel 的 xy 子块就是纯面内问题，
    两支本征值/本征矢都是面内极化的 LA/TA —— 从根上排除面外支（slab 的 C33/C44≈0
    会把原版 ADP 除爆）。
    """
    if np.isscalar(elastic) or np.asarray(elastic).ndim < 4:
        raise ValueError("amset2d：2D 的 ADP 需要完整 6×6 原始 slab 弹性张量"
                         "（settings.yaml 里的 elastic_constant 不能是各向同性标量）")
    # q 完全落在真空方向（面内投影恰好为 0）时面内方向没有定义：这些点在二维积分里是
    # 零测度，取任一固定面内方向即可（ADP/PIE 的分子本身与 q̂ 无关地有限/趋零）。
    u = np.array(u, copy=True)
    deg = np.linalg.norm(u, axis=1) < 0.5
    if np.any(deg):
        _note_qz_only(int(np.count_nonzero(deg)))
        u[deg] = (1.0, 0.0, 0.0)
    g = get_christoffel_tensors(elastic, u)[:, :2, :2]
    w, vec = np.linalg.eigh(g)
    ok = ~deg
    if np.any(w[ok] <= 0) or np.any(w[ok, 0] < 1e-2 * w[ok, 1]):
        raise ValueError(
            "amset2d：面内弹性张量非正定，或面内剪切刚度异常小"
            "（最小/最大 < 1%）——几乎总是 elastic_constant 没从 VASP 顺序 "
            "(XX YY ZZ XY YZ ZX) 重排成标准 Voigt (XX YY ZZ YZ XZ XY) 导致的："
            "不重排会把面内剪切 C66(XY) 换成接近零的 C_zxzx。"
            "请用 gen_step14_amset2d.py 重新生成 settings.yaml。")
    v = np.zeros((len(u), 2, 3))
    v[:, :, :2] = np.transpose(vec, (0, 2, 1))
    return w, v


def _zero_screening(amset_data):
    return np.zeros(amset_data.fermi_levels.shape)


# [AMSET2D_DEBUG] 诊断开关（AMSET2D_DEBUG=1 时打印散射核内部的量）
import os as _os_dbg
_DEBUG = _os_dbg.environ.get("AMSET2D_DEBUG", "") not in ("", "0", "false")
_DBG_N = [0]

_QZ_ONLY_SEEN = [False]


def _note_qz_only(n):
    """只报一次：有多少 q 点的面内投影恰好为 0（这些点用固定面内方向兜底）。"""
    if not _QZ_ONLY_SEEN[0]:
        _QZ_ONLY_SEEN[0] = True
        print("[amset2d] 提示：本次出现 %d 个 q∥=0（q 沿真空方向）的采样点，"
              "已用固定面内方向兜底（二维积分里是零测度）" % n, flush=True)


# ---------------- ADP ----------------
class AcousticDeformation2D(_ADP):
    """只保留面内传播、面内极化的 LA/TA 支；弹性常数用原始 slab 值（= C_2D/c）。

    [0033] 开了谷间散射（IVS）的载流子，ADP 只留**谷内**（k 与 k+q 落在同一个成员谷）：
    AMSET 的 ADP 对整个布里渊区的 q 积分，大 q 的 K→K′/K→Q 用长波声速和带边形变势算，
    既不物理、又与 IVS 的谷间通道重复计数。掩码数据在 _ivm（主进程 from_amset_data 建，
    经 to_reference/from_reference 带进散射子进程）；没有 _ivm 时与以前逐字节相同。
    """
    name = "ADP"

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        obj = super().from_amset_data(materials_properties, amset_data)
        obj._ivm = _adp_mask_setup(amset_data)
        return obj

    def to_reference(self):
        return super().to_reference() + (getattr(self, "_ivm", None),)

    @classmethod
    def from_reference(cls, *ref):
        obj = super().from_reference(*ref[:-1])
        obj._ivm = ref[-1]
        return obj

    def _iv_keep(self, unit_q, norm_q_sq, spin, band_idx, kpoint):
        """谷内掩码：k+q 与 k 在同一个成员谷 -> True。该载流子没开掩码 -> None（不掩）。"""
        ivm = getattr(self, "_ivm", None)
        if not ivm or kpoint is None:
            return None
        vb = (getattr(self, "vb_idx", None) or {}).get(spin)
        centers = ivm.get("electron" if (vb is not None and band_idx > vb) else "hole")
        if centers is None:
            return None
        q = np.sqrt(norm_q_sq)[:, None] * np.asarray(unit_q)
        kp = q @ ivm["rinv"] + np.asarray(kpoint, float)
        return _DPO.same_member(kpoint, kp, centers, ivm["rlat"])

    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        f = self._factor_core(unit_q, norm_q_sq, spin, band_idx, kpoint, velocity)
        keep = self._iv_keep(unit_q, norm_q_sq, spin, band_idx, kpoint)
        return f if keep is None else f * keep

    def _factor_core(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        # （amset_ir_fix 的 ADP 星平均对本方法做旋转平均，谷内掩码在平均之后按未旋转的 q 乘）
        _, _, u = _inplane(unit_q, norm_q_sq)
        w, v = _inplane_modes(self.elastic_constant, u)
        if hasattr(self.deformation_potential, "interpolate"):
            D = np.abs(self.deformation_potential.interpolate(spin, [band_idx], [kpoint])[0])
            D = D + np.outer(velocity, velocity)                   # 与 AMSET 原版一致
            f = np.zeros(len(unit_q))
            for b in range(2):                                     # 只求和面内两支
                strain = u[:, :, None] * v[:, b, None, :]
                f = f + np.tensordot(strain, D) ** 2 / w[:, b]
        elif self.is_metal:
            f = self.deformation_potential ** 2 / w[:, 1]           # 金属：只留 LA 支
        else:
            idx = 1 if band_idx > self.vb_idx[spin] else 0
            f = self.deformation_potential[idx] ** 2 / w[:, 1]
        # [AMSET2D_DEBUG] 诊断用：把核里的量打出来，便于与解析参考逐项对齐
        # （默认关；只打前 3 次调用，不影响生产运行）
        if _DEBUG and _DBG_N[0] < 3:
            _DBG_N[0] += 1
            _k = 0
            print("[amset2d][DEBUG] call=%d |q|=%.6e a^-1  u=(%.4f,%.4f,%.4f)  "
                  "w_TA=%.6e  w_LA=%.6e  D_proj=%s  f=%.6e"
                  % (_DBG_N[0], float(np.sqrt(norm_q_sq[_k])),
                     float(u[_k, 0]), float(u[_k, 1]), float(u[_k, 2]),
                     float(w[_k, 0]), float(w[_k, 1]),
                     (np.round(np.abs(D), 4).tolist() if np.size(D) <= 9 else "array"),
                     float(f[_k])), flush=True)
        return f[None, None] * np.ones(self.fermi_levels.shape + norm_q_sq.shape)


# ---------------- ODP / IVS（[0033] 常数形变势，默认不开） ----------------
# 谷内光学形变势（ODP）与谷间散射（IVS）：Jin 2014 / Kaasbjerg 2012 的常数形变势模型，率由
# dp_optical.py 解析给出（发射/吸收阈值精确、各向同性 -> 动量弛豫率 = 散射率），走 AMSET 的
# basic 接口（每 k 点直接给率，1/s），不经过四面体积分与 worker。
# 参数在 2d_correction.json 的 dp_optical 记录里（gen 从 step.conf 的 DP_OPTICAL_E/H 与 S8.2 谷表写入）：
#   {"electron": {"valleys": [...], "channels": [...]}, "hole": {...},
#    "adp_intravalley": {"electron": bool, "hole": bool}}
# mechanisms_2d 不含 ODP/IVS、adp_intravalley 全关时，本段对运行没有任何影响。
# （旧的"弹性近似 ODP"实验类（2N+1、端到端标定系数 0.88）已由本段取代。）
try:
    import dp_optical as _DPO
except ImportError:          # 老运行目录没有 dp_optical.py：不开 ODP/IVS/谷内掩码就用不到
    _DPO = None
_DPO_REC = _REC.get("dp_optical") or {}
_CARRIERS = ("electron", "hole")
_DPO_KINDS = {"ODP": frozenset({"intra"}), "IVS": frozenset({"inter_star", "inter"})}


def _need_dpo():
    if _DPO is None:
        raise ImportError("amset2d：2d_correction.json 开了 ODP/IVS 或 ADP 谷内掩码，但运行目录里没有"
                          " dp_optical.py（请用 gen_step14_amset2d.py 重新生成 step8.4_amset2d）")
    return _DPO


def _split_bands(amset_data, spin):
    """该自旋的能带 -> (电子带, 空穴带) 布尔数组。常数形变势模型按载流子分谷，金属没有意义。"""
    if amset_data.is_metal or not amset_data.vb_idx:
        raise ValueError("amset2d：ODP/IVS/ADP 谷内掩码需要半导体（AMSET 的 vb_idx），本体系被判为金属")
    cb = np.arange(len(amset_data.energies[spin])) > amset_data.vb_idx[spin]
    return cb, ~cb


def _adp_mask_setup(amset_data):
    """ADP 谷内掩码的数据（{"electron"/"hole": 成员谷中心 (n,3) 或 None, "rlat", "rinv"}）；都不开 -> None。
    rlat 与 calculate_rate 里的同一个（amset_data.structure 的倒格子，bohr⁻¹）：q 就是在这个笛卡尔系里给的。"""
    flags = _DPO_REC.get("adp_intravalley") or {}
    if not any(flags.get(c) for c in _CARRIERS):
        return None
    dpo = _need_dpo()
    rlat = np.asarray(amset_data.structure.lattice.reciprocal_lattice.matrix, float)
    out = {"rlat": rlat, "rinv": np.linalg.inv(rlat)}
    for car in _CARRIERS:
        vs = (_DPO_REC.get(car) or {}).get("valleys") or []
        out[car] = dpo.member_table(vs)[0] if (flags.get(car) and vs) else None
    if out["electron"] is not None or out["hole"] is not None:
        _split_bands(amset_data, amset_data.spins[0])          # 金属直接拒绝
        return out
    return None


_VALLEY_CHECKED = [False]
VALLEY_DK_WARN = 0.05          # |Δk|/|b1|：谷区域里的能量最低点离谷中心超过这么远就告警


def _valley_check(amset_data):
    """每次运行打一次：谷表与 AMSET 网格上的能带对不对得上（谷底是否落在给定的谷中心）。"""
    if _VALLEY_CHECKED[0]:
        return
    _VALLEY_CHECKED[0] = True
    if getattr(amset_data, "_soc", False):
        print("[amset2d][WARN] ODP/IVS：SOC 能带（旋量带）—— 目标谷态密度把两种自旋都算进去了；"
              "通道若应守恒自旋（如 MoS₂ 价带 K↔K′ 的自旋-谷锁定），率会偏大到 2 倍、谷间可能根本不该有。"
              "参数要与文献的 SOC 口径一致，结论慎用。", flush=True)
    rlat = amset_data.structure.lattice.reciprocal_lattice.matrix
    spin = amset_data.spins[0]
    for car, sel, sign in zip(_CARRIERS, _split_bands(amset_data, spin), (+1, -1)):
        vs = (_DPO_REC.get(car) or {}).get("valleys") or []
        if not vs or not sel.any():
            continue
        rows = _DPO.valley_report(np.asarray(amset_data.energies[spin])[sel], amset_data.kpoints, rlat, vs, sign)
        txt = "; ".join("%s 谷底离中心 %s、高于带边 %s" % (
            r["label"], "—" if r["dk_rel"] is None else "%.3f|b1|" % r["dk_rel"],
            "—" if r["dE_eV"] is None else "%.3f eV" % r["dE_eV"]) for r in rows)
        bad = [r["label"] for r in rows if r["dk_rel"] is not None and r["dk_rel"] > VALLEY_DK_WARN]
        used = {lab for h in (_DPO_REC.get(car) or {}).get("channels") or [] for lab in (h["from"], h["to"])}
        bad += [r["label"] for r in rows if r["dk_rel"] is None and r["label"] in used]   # 有通道却分不到 k 点
        print("[amset2d]%s 谷表核对（%s）：%s%s" % (
            "[WARN]" if bad else "", car, txt,
            ("  ——%s 的谷底不在给定中心或分不到 k 点：谷表与 AMSET 的晶胞/网格不同源？（S8.2 谷表来自 S3 晶胞）"
             % ",".join(bad))
            if bad else ""), flush=True)


def _dpo_rates(amset_data, kinds):
    """{spin: (ndop, ntemp, nbands, nk)} 的率（1/s），按载流子各用自己的谷表与通道；
    没有通道的态加 RATE_FLOOR（> AMSET 的 1e6 零率阈值，免得被最近邻插值顶替）。与掺杂无关 -> 广播视图。"""
    dpo = _need_dpo()
    _valley_check(amset_data)
    rlat = np.asarray(amset_data.structure.lattice.reciprocal_lattice.matrix, float)
    mass = float(amset_data.structure.composition.weight) * dpo.AMU           # 原胞质量 kg
    kf = np.asarray(amset_data.kpoints, float)
    temps = np.asarray(amset_data.temperatures, float)
    ndop = amset_data.fermi_levels.shape[0]
    out = {}
    for spin in amset_data.spins:
        E = np.asarray(amset_data.energies[spin], float)
        r = np.zeros((len(temps),) + E.shape)
        for car, sel, sign in zip(_CARRIERS, _split_bands(amset_data, spin), (+1, -1)):
            c = _DPO_REC.get(car) or {}
            if sel.any() and c.get("channels"):
                r[:, sel] = dpo.carrier_rates(E[sel], kf, rlat, c.get("valleys") or [], c["channels"],
                                              temps, mass, kinds, sign)
        r += dpo.RATE_FLOOR
        # 按 AMSET 的不可约映射对称化：每个 k 点取它代表点的率。窗口内 AMSET 本来就只用代表点的率，
        # 窗口外由 _fill_keep_dpo 写回这份；mesh 只存代表点、postprocess 再按映射展开 —— 三处必须同一份。
        # 谷归属按"最近中心"，等距（Voronoi 边界）时的取舍不随对称操作协变，不对称化会让窗口外的率与 mesh
        # 不一致，--check-reproduce 失败、整步 FAIL（复审实跑：矩形胞 1044 个 k 点，重积分差 0.14）。
        ir, i2f = getattr(amset_data, "ir_kpoints_idx", None), getattr(amset_data, "ir_to_full_kpoint_mapping", None)
        if ir is not None and i2f is not None:
            r = r[..., np.asarray(ir)][..., np.asarray(i2f)]
        out[spin] = np.broadcast_to(r[None], (ndop,) + r.shape)
    return out


class _DPOptical2D(_ABS_BASIC):
    """常数形变势光学声子散射的公共部分（AMSET basic 接口）。"""
    required_properties = ()          # 参数走 2d_correction.json，不进 settings 校验
    kinds = frozenset()

    def __init__(self, properties, doping, temperatures, nbands, rates):
        super().__init__(properties, doping, temperatures, nbands)
        self._rates = rates

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        rates = _dpo_rates(amset_data, cls.kinds)
        _DPO_FULL[cls.name] = rates          # FD 窗口外的态也要用解析率（见下面 _fill_keep_dpo）
        return cls(cls.get_properties(materials_properties), amset_data.doping,
                   amset_data.temperatures, cls.get_nbands(amset_data), rates)

    @property
    def rates(self):
        return self._rates


# FD 窗口外的态：AMSET 的 fill_rates_outside_cutoffs 把每个机制窗口外的率统一填成"窗口内全部态（电子带+空穴带）
# 的对数平均"。ODP/IVS 的率是解析的、全网格都有，没必要填；而且只给一种载流子写了通道时，另一种载流子的
# RATE_FLOOR（2e6）会把对数平均拉到 ~1e9，窗口外的态于是成了"几乎不散射"的通道，单机制迁移率被放大几十倍
# （合成双谷实测：电子 IVS 11778 vs 应有 ~350 cm²/Vs）。这里在 AMSET 填完之后把 ODP/IVS 窗口外的率改回解析值。
_DPO_FULL = {}


def _install_fill_keep_dpo():
    import amset.core.data as _AD
    if getattr(_AD.AmsetData.fill_rates_outside_cutoffs, "_amset2d_dpo", False):
        return
    _orig = _AD.AmsetData.fill_rates_outside_cutoffs

    def _fill_keep_dpo(self, fill_value=None):
        out = _orig(self, fill_value)
        labels = list(getattr(self, "scattering_labels", None) or [])
        hit = [lab for lab in labels if lab in _DPO_FULL]
        if hit and fill_value is None:
            lo, hi = self.fd_cutoffs
            for spin, e in self.energies.items():
                mask = (e < lo) | (e > hi)
                for lab in hit:
                    full = _DPO_FULL[lab].get(spin)
                    arr = self.scattering_rates[spin][labels.index(lab)]
                    if full is not None and full.shape == arr.shape:
                        arr[:, :, mask] = full[:, :, mask]
            print("[amset2d] %s：FD 窗口外保留解析率（不用 AMSET 的对数平均填充）" % "/".join(hit), flush=True)
        return out

    _fill_keep_dpo._amset2d_dpo = True
    _AD.AmsetData.fill_rates_outside_cutoffs = _fill_keep_dpo


class OpticalDeformation2D(_DPOptical2D):
    """谷内光学形变势（X>X 通道）。"""
    name = "ODP"
    kinds = _DPO_KINDS["ODP"]


class Intervalley2D(_DPOptical2D):
    """谷间散射（X>X′ 同星其它谷、X>Y 到别的星）。"""
    name = "IVS"
    kinds = _DPO_KINDS["IVS"]


# ---------------- POP ----------------
class PolarOptical2D(PolarOpticalScattering):
    """二维 Fröhlich：F = 2πc·[1/ε̃∞(q) − 1/ε̃0(q)]/q，ε(q) 为张量 Keldysh 形式。

    q→0 趋于有限值 2πc·Δr/ε_env²，大 q 按 1/q 衰减（与 Sohier 2016 的两个极限一致）。
    pop_lo_dispersion=True 时乘 √(ε̃0(q)/ε̃∞(q))（二维 LST：ω_LO(q)²/ω_TO² = ε̃0/ε̃∞），
    与玻恩电荷写法等价；此时 settings 的 pop_frequency 取 Γ 点**面内极性 TO** 模频率。
    """
    name = "POP"                                  # 保留原名：_get_cutoff_pad 只认 "POP"

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        obj = super().from_amset_data(materials_properties, amset_data)
        obj.r_inf = _r_tensor("eps_inf_slab")
        obj.r_0 = _r_tensor("eps_static_slab")
        if materials_properties["free_carrier_screening"]:
            obj.q_tf = _qtf(amset_data)
        else:
            obj.q_tf = _zero_screening(amset_data)
        return obj

    def factor(self, unit_q, norm_q_sq, emission, f):
        occ = self.n_po + 1 - f if emission else self.n_po + f
        qv, qn, _ = _inplane(unit_q, norm_q_sq)
        a_i, a_0 = _q_eps(qv, qn, self.r_inf), _q_eps(qv, qn, self.r_0)
        num = np.einsum("ni,ij,nj->n", qv, self.r_0 - self.r_inf, qv)
        A_i = a_i[None, None] + self.q_tf[..., None]
        A_0 = a_0[None, None] + self.q_tf[..., None]
        F = 2 * np.pi * _C * num[None, None] / (A_i * A_0)
        if _LO_DISP:      # ω_LO(q)/ω_TO = sqrt(ε̃0(q)/ε̃∞(q))
            F = F * np.sqrt((a_0 + 1e-30) / (a_i + 1e-30))[None, None]
        return occ[..., None] * F


# ---------------- IMP ----------------
class IonizedImpurity2D(_IMP):
    """二维屏蔽库仑杂质：F = c²(2π)² e^{-2qd} / (q·ε0(q) + q_TF)²。

    映射：G_s = c·G_2D，G_2D = N_2D|V_2D|²，N_2D = c·N_3D，V_2D = 2πZe^{-qd}/(q·ε0(q)+q_TF)
    ⇒ 在 AMSET 的三维 prefactor（= Z²N_3D）上乘 c²(2π)²/(…)。必须带屏蔽，否则 q→0 发散。
    """
    name = "IMP"

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        obj = super().from_amset_data(materials_properties, amset_data)
        obj.r_0 = _r_tensor("eps_static_slab")
        obj.q_tf = _qtf(amset_data)
        return obj

    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        qv, qn, _ = _inplane(unit_q, norm_q_sq)
        a_0 = _q_eps(qv, qn, self.r_0)
        num = _C ** 2 * (2 * np.pi) ** 2 * np.exp(-2 * qn * _D_IMP)
        return num[None, None] / (a_0[None, None] + self.q_tf[..., None]) ** 2


# ---------------- PIE ----------------
class Piezoelectric2D(_PIE):
    """二维面内压电：F = c²(2π)² Σ_b (q̂q̂:e·v_b)²/w_b · q²/(q·ε0(q) + q_TF)²。

    只处理面内压电（Janus 结构的面外偶极层在这里被忽略）；e 与 w 都用原始 slab 值。
    [0031] e 取 VASP 的 **including ionic contribution**（离子弛豫）总张量，而压电势是
    静态势，屏蔽必须用**静态** ε0（r0）而不是 ε∞；旧代码配 r∞ 属于张量口径错配
    （默认不开 PIE 时无影响；MoS2 这类弱压电材料也接近可忽略）。
    """
    name = "PIE"

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        e_au = np.array(materials_properties["piezoelectric_constant"], float)
        e_au = e_au * coulomb_to_au / m_to_bohr ** 2      # C/m² -> a.u.（原版会就地改 settings）
        obj = super().from_amset_data(materials_properties, amset_data)
        obj.e_raw = e_au
        obj.r_0 = _r_tensor("eps_static_slab")
        if materials_properties["free_carrier_screening"]:
            obj.q_tf = _qtf(amset_data)
        else:
            obj.q_tf = _zero_screening(amset_data)
        return obj

    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        qv, qn, u = _inplane(unit_q, norm_q_sq)
        w, v = _inplane_modes(self.elastic_constant, u)
        e2 = self.e_raw[:2, :2, :2]
        s = np.zeros(len(unit_q))
        for b in range(2):
            cpl = np.einsum("ijk,ni,nj,nk->n", e2, u[:, :2], u[:, :2], v[:, b, :2])
            s = s + cpl ** 2 / w[:, b]
        a_0 = _q_eps(qv, qn, self.r_0)
        return (_C ** 2 * (2 * np.pi) ** 2 * (s * qn ** 2)[None, None]
                / (a_0[None, None] + self.q_tf[..., None]) ** 2)


# ---------------- 注册 ----------------
_CLASSES = {"ADP": AcousticDeformation2D, "POP": PolarOptical2D,
            "IMP": IonizedImpurity2D, "PIE": Piezoelectric2D,
            "ODP": OpticalDeformation2D, "IVS": Intervalley2D}
_ACTIVE = [k for k in ("ADP", "POP", "IMP", "PIE", "ODP", "IVS") if k in _MECHS]
for _k in ("ODP", "IVS"):    # [0033] 开了就必须有对应种类的通道（gen 保证；手写 json 在这里拦）
    if _k in _MECHS and not any(h.get("kind") in _DPO_KINDS[_k]
                                for _c in _CARRIERS for h in (_DPO_REC.get(_c) or {}).get("channels") or []):
        raise ValueError("amset2d：mechanisms_2d 含 %s，但 2d_correction.json 的 dp_optical 里没有 %s 通道"
                         % (_k, "/".join(sorted(_DPO_KINDS[_k]))))
for _k in _ACTIVE:
    C._scattering_mechanisms[_k] = _CLASSES[_k]
if {"ODP", "IVS"} & set(_ACTIVE):
    _install_fill_keep_dpo()
if "ADP" in _MECHS:
    C.AcousticDeformationPotentialScattering = AcousticDeformation2D   # worker 用此名重建 ADP

_BANNER = ("[amset2d] amset %s | c=%.3f A | r_inf=%s | mechanisms=%s"
           % (amset.__version__, _C / A2B, np.round(_r_tensor("eps_inf_slab").tolist(), 3),
              ",".join(_ACTIVE)))
print(_BANNER, flush=True)
if _DPO_REC and ({"ODP", "IVS"} & set(_ACTIVE) or any((_DPO_REC.get("adp_intravalley") or {}).values())):
    _need_dpo()
    _ivf = _DPO_REC.get("adp_intravalley") or {}
    print("[amset2d] ODP/IVS（常数形变势，D eV/Å @ ħω meV）：%s | ADP 谷内掩码：电子 %s，空穴 %s"
          % (_DPO.describe(_DPO_REC), "开" if _ivf.get("electron") else "关",
             "开" if _ivf.get("hole") else "关"), flush=True)
# r∞/Δr 自检（gen_step14 落盘的值）。曾把体相 ε0/ε∞ 外推的 Δr=43.5 Å 当 Sohier 的
# 0.53 Å 用，POP 单机制迁移率差约两个量级；这里每次运行都把真正喂进核的两个数打出来。
_RCHK = _REC.get("r_inf_delta_r_check") or {}
if _RCHK.get("available"):
    print("[amset2d] POP: r_inf=%.3f A | delta_r=%.4f A | delta_r/r_inf=%.2f%%"
          " (G_2D = 2*pi*c*Delta_r/[(1+r_inf*q)(1+r_0*q)])"
          % (_RCHK["r_inf_A"], _RCHK["delta_r_A"],
             100.0 * _RCHK["delta_r_over_r_inf"]), flush=True)
else:
    print("[amset2d][WARN] 2d_correction.json 里没有 r_inf_delta_r_check ——"
          " 无法核对 POP 的 Δr 量级（老记录或手写 json）", flush=True)

# spawn 子进程：ADP 以 reference 传给子进程、并在子进程里按上面的类名重建，
# 所以子进程也必须导入本插件（否则静默退回原版 ADP）。这里包装 worker 打印实际生效的类。
if not getattr(C, "_amset2d_patched", False):
    _orig_worker = C.scattering_worker

    def worker_2d(*args, **kwargs):
        print("[amset2d] worker: " + " ".join(
            "%s=%s" % (k, C._scattering_mechanisms[k].__name__)
            for k in ("ADP", "POP", "IMP", "PIE")), flush=True)
        return _orig_worker(*args, **kwargs)

    C.scattering_worker = worker_2d
    C._amset2d_patched = True


# =====================================================================
# 输运前守卫（V130）—— 不改数值，只把 NaN/inf/零率说清楚
#
# WS2（KZ_CAP + 全网格 h5）两次在 BoltzTraP2 calc_Onsager_coefficients 的 pinv(L11) 报
# "SVD did not converge"：L11 里有 NaN/inf。可能的来源：某机制在 FD 窗口内的率 ≤ 0 或非有限
# （τ = 1/Γ -> ∞，单机制迁移率 1/0）、费米能级 NaN、速度 NaN。这里在 AMSET 的两个输运入口
# （_calculate_transport_properties / _calculate_mobility）前扫一遍，找到就打印前几处（机制/自旋/带/
# 掺杂/温度/k 点/离费米能多远），pinv 失败时再打一次并附上下文重新抛出。
# =====================================================================
def transport_report(amset_data, max_lines=8):
    """返回问题描述列表（空 = 没发现）。逐 (机制, 掺杂, 温度) 切片扫，避免整块布尔数组占内存。"""
    from amset.constants import hartree_to_ev
    out = []
    fl = np.asarray(amset_data.fermi_levels, float)
    if not np.all(np.isfinite(fl)):
        out.append("fermi_levels 有非有限值：%s" % np.argwhere(~np.isfinite(fl)).tolist()[:max_lines])
    labels = list(getattr(amset_data, "scattering_labels", []) or [])
    kp = np.asarray(amset_data.kpoints)
    for spin, r in amset_data.scattering_rates.items():
        r = np.asarray(r)
        e = np.asarray(amset_data.energies[spin])
        vv = np.asarray(amset_data.velocities_product[spin])
        if not np.all(np.isfinite(vv)):
            out.append("velocities_product[%s] 有非有限值 %d 个" % (spin.name, int(np.count_nonzero(~np.isfinite(vv)))))
        for m, n, t in np.ndindex(r.shape[:3]):
            sl = r[m, n, t]
            bad = ~np.isfinite(sl) | (sl <= 0)
            if not bad.any():
                continue
            nb = int(np.count_nonzero(bad))
            b, k = np.argwhere(bad)[0]
            de = (e[b, k] - fl[n, t]) * hartree_to_ev if np.isfinite(fl[n, t]) else float("nan")
            out.append("%s[%s] 掺杂#%d 温度#%d：%d 个率 ≤0 或非有限；首处 带 %d k=%s 值 %s，E−E_F = %+.3f eV"
                       % (labels[m] if m < len(labels) else "mech%d" % m, spin.name, n, t, nb, b,
                          np.round(kp[k], 4).tolist(), sl[b, k], de))
            if len(out) >= max_lines:
                out.append("……（只列前 %d 处）" % max_lines)
                return out
    return out


def _guard_transport(fn):
    def wrapped(amset_data, *a, **k):
        if not getattr(amset_data, "_amset2d_guarded", False):
            amset_data._amset2d_guarded = True
            rep = transport_report(amset_data)
            if rep:
                print("[amset2d][WARN] 输运前发现坏值（BoltzTraP2 的 pinv 可能因此 SVD 不收敛）：\n  - "
                      + "\n  - ".join(rep), flush=True)
        try:
            return fn(amset_data, *a, **k)
        except np.linalg.LinAlgError as e:
            rep = transport_report(amset_data)
            raise RuntimeError("[amset2d] 输运积分 %s 失败（%s）。坏值：%s"
                               % (fn.__name__, e, "; ".join(rep) if rep else "散射率/费米能级/速度都有限 —— "
                                  "请把 amset.log 与本段报错一并反馈")) from e
    wrapped.__name__ = fn.__name__
    wrapped._amset2d_guard = True
    return wrapped


import amset.core.transport as _T  # noqa: E402
for _fn in ("_calculate_transport_properties", "_calculate_mobility"):
    if hasattr(_T, _fn) and not getattr(getattr(_T, _fn), "_amset2d_guard", False):
        setattr(_T, _fn, _guard_transport(getattr(_T, _fn)))

# V131/V132：求费米能级的稳健化在 amset_fermi_fix（3D 的 S8 也挂它）：AMSET 贪心搜索失败时改用二分法、
# 宽松容差下偏差 > 1% 不接受、NaN 报错（V131 的 NaN 守卫从这里移过去）。gen 把它复制进运行目录；
# 这里再导入一次，postprocess_intrinsic / tools 只加载本插件时也生效。
try:
    import amset_fermi_fix  # noqa: F401,E402
except ImportError:
    print("[amset2d] [WARN] 找不到 amset_fermi_fix.py —— 费米能级用 AMSET 原搜索（低温低掺杂可能失败，"
          "宽松容差的解不拦）", file=sys.stderr)

# =====================================================================
# patch_kz_cap（V125）—— 插值网格在真空方向只留 2·r+1 层（2d_correction.json 的 kz_cap_rmax；
# gen 的 KZ_CAP_2D 出厂 on（V139），写 off 时没有这个键、不截断）
#
# 单层的能带沿 k_z 是平的（真空隔开，层间跳跃 ~0），本插件的 2D 核也只用 q∥ —— 沿 k_z 的积分
# 只贡献常数因子 2π/c。可 BoltzTraP2 的 equivalence 是实空间**球**，面内要 263 个格点时真空方向
# 也跟着要 ~33 层（WSe2/WS2 的 c≈23 Å），而 AMSET 的散射成本 ≈ (不可约点数 ∝ N_z) × (等能面上的
# 四面体数 ∝ N_z) ∝ N_z²。把 |R_z| > r 的 equivalence 去掉：面内半径不变（面内网格不变），
# k_z 只剩 2r+1 层。合成六方单层（2D ADP+IMP）实测：85×85×11 -> 85×85×3，迁移率差 < 0.15%。
# 真实重叠下多出的差别只有 slab 波函数沿 z 的形状因子在 k_z 上的平均（估计 ~0.1%）。
# 护栏：运行时用输入能带核对"沿 k_z 平"（同一面内坐标、不同 k_z 的能量差 ≤ kz_flat_tol_ev，
# 默认 5 meV，只看带边 ±1 eV 内的带）；核不了或不平 -> 不截断、照原样跑并告警。
# =====================================================================
KZ_FLAT_TOL_EV_DEFAULT = 0.005
_KZ = {"rmax": None, "axis": None, "spread_ev": None, "n_groups": 0, "active": False}


def kz_cap_filter(equivalences, axis, rmax):
    """只保留真空轴分量 |R_axis| <= rmax 的 equivalence（整个星都满足才留）。"""
    return [e for e in equivalences if np.abs(np.asarray(e)[:, axis]).max() <= rmax]


def vacuum_axis(lattice_matrix, normal=(0.0, 0.0, 1.0), tol=1e-3):
    """出面的晶格矢量下标。条件只是**另两个矢量 ⟂ 层法向**（此时 b_axis ∥ 法向，k_z 就是这一轴；
    c 本身可以倾斜）。不满足 -> None。"""
    L = np.asarray(lattice_matrix, float)
    n = np.asarray(normal, float)
    n = n / np.linalg.norm(n)
    cos = np.abs(L @ n) / np.linalg.norm(L, axis=1)
    ax = int(np.argmax(cos))
    if any(cos[i] > tol for i in range(3) if i != ax):
        return None
    return ax


def kz_spread_ev(kpoints, energies, axis, emin, emax):
    """同一面内坐标、不同 k_z 的 k 点之间，窗口 [emin, emax] 内各带的最大能量差（eV）。
    kpoints (nk, 3) 分数坐标；energies (nb, nk) eV。返回 (spread, 可比较的面内坐标组数)。"""
    k = np.mod(np.round(np.asarray(kpoints, float), 6), 1.0)
    k[k > 1 - 1e-6] = 0.0
    inp = [i for i in range(3) if i != axis]
    keys = np.round(k[:, inp], 5)
    kz = np.round(k[:, axis], 5)
    E = np.asarray(energies, float)
    bands = np.where(((E >= emin) & (E <= emax)).any(axis=1))[0]
    groups = {}
    for i, key in enumerate(map(tuple, keys)):
        groups.setdefault(key, []).append(i)
    spread, n = 0.0, 0
    for idx in groups.values():
        if len(set(kz[idx])) < 2:
            continue
        n += 1
        sub = E[np.ix_(bands, idx)]
        spread = max(spread, float((sub.max(axis=1) - sub.min(axis=1)).max()) if len(bands) else 0.0)
    return spread, n


def _kz_prepare(band_structure):
    """在 Interpolator 构造前核对：真空轴、沿 k_z 平。通过才让 equivalence 截断生效。"""
    from pymatgen.electronic_structure.core import Spin
    _KZ["active"] = False
    ax = vacuum_axis(band_structure.structure.lattice.matrix, _N)
    if ax is None:
        print("[amset2d][WARN] kz_cap：晶格的面内两个矢量不垂直于层法向 -> 不截断 k_z（照原样跑）", flush=True)
        return
    kp = np.array([k.frac_coords for k in band_structure.kpoints])
    ef = band_structure.efermi
    # 没有 efermi（vasprun 缺标签）又取不到带边 -> 全部带都核（只会更保守：多核的带不平就不截断）
    emin, emax = (ef - 1.0, ef + 1.0) if ef is not None else (-np.inf, np.inf)
    try:
        vbm = band_structure.get_vbm()["energy"]
        cbm = band_structure.get_cbm()["energy"]
        if vbm is not None and cbm is not None:
            emin, emax = vbm - 1.0, cbm + 1.0
    except Exception:                                          # noqa: BLE001
        pass
    spread, ng = 0.0, 0
    for spin in band_structure.bands:
        s, n = kz_spread_ev(kp, band_structure.bands[spin], ax, emin, emax)
        spread, ng = max(spread, s), ng + n
    tol = float(_REC.get("kz_flat_tol_ev") or KZ_FLAT_TOL_EV_DEFAULT)
    _KZ.update(axis=ax, spread_ev=spread, n_groups=ng)
    if ng == 0:
        print("[amset2d][WARN] kz_cap：输入 k 点里找不到同一面内坐标、不同 k_z 的点（S3 真空方向只有 1 层？）"
              " -> 核不了能带是否沿 k_z 平，不截断（照原样跑）", flush=True)
        return
    if spread > tol:
        print("[amset2d][WARN] kz_cap：能带沿 k_z 不平（最大差 %.1f meV > %.1f meV，%d 组）—— 真空不够或不是单层 ->"
              " 不截断 k_z（照原样跑）" % (1e3 * spread, 1e3 * tol, ng), flush=True)
        return
    _KZ["active"] = True
    print("[amset2d] kz_cap：沿 k_z 最大能量差 %.2f meV（%d 组，阈值 %.1f meV）-> equivalence 截到 |R_%s| <= %d"
          % (1e3 * spread, ng, 1e3 * tol, "abc"[ax], _KZ["rmax"]), flush=True)


_KZ_RMAX = _REC.get("kz_cap_rmax")
if _KZ_RMAX and not getattr(C, "_amset2d_kzcap", False):
    import inspect as _insp
    import types as _types
    import amset.interpolation.bandstructure as _bsm

    _KZ["rmax"] = int(_KZ_RMAX)
    if "sphere.get_equivalences(" not in "".join(_insp.getsource(_bsm.Interpolator.__init__).split()):
        print("[amset2d][WARN] kz_cap：amset %s 的 Interpolator 不再调 sphere.get_equivalences -> 不截断"
              % amset.__version__, flush=True)
    else:
        _orig_init_kz = _bsm.Interpolator.__init__
        _orig_eq_kz = _bsm.sphere.get_equivalences

        def _init_kz(self, band_structure, *a, **kw):
            _kz_prepare(band_structure)
            return _orig_init_kz(self, band_structure, *a, **kw)

        def _eq_kz(atoms, magmom, nkpt):
            eq = _orig_eq_kz(atoms, magmom, nkpt)
            if not _KZ["active"]:
                return eq
            out = kz_cap_filter(eq, _KZ["axis"], _KZ["rmax"])
            print("[amset2d] kz_cap：equivalence %d -> %d（k_z 层数 %d）"
                  % (len(eq), len(out), 2 * _KZ["rmax"] + 1), flush=True)
            return out

        _bsm.Interpolator.__init__ = _init_kz
        _ns = _types.SimpleNamespace(**{k: getattr(_bsm.sphere, k) for k in dir(_bsm.sphere)
                                        if not k.startswith("__")})
        _ns.get_equivalences = _eq_kz
        _bsm.sphere = _ns
        C._amset2d_kzcap = True

# =====================================================================
# patch_valley_overlap（2026-09-17）—— 诊断专用：**按谷区分**的重叠
#
# 目的（VERIFICATION V23.5 的下一步）：把"残差来自哪里"彻底分开。
#   同谷 |I|^2 = 1、不同谷 |I|^2 = 0，其余输入完全不变，跑一次。
#   若结果回到单谷 Takagi（电子 205 / 空穴 1370），说明散射核与参考值都没问题、
#   残差完全来自真实重叠的分布；若不吻合，残差来自参考模型本身。
#
# ★ 默认关闭，只由环境变量打开（生产绝不触发）：
#     AMSET2D_VALLEY_OVERLAP=1
#     AMSET2D_VALLEY_CENTERS="0.33333,0.33333,0;-0.33333,-0.33333,0"
#
# 实现：用缓存系数时重叠是在 C._get_overlap 里由系数内积算的，那里看不到 k 点；
#   所以同时包装 C.calculate_rate 把当前 tbs.kpoints 存到全局，再在 C._get_overlap
#   里按起点/终点的谷归属把跨谷通道清零。
# =====================================================================
if os.environ.get("AMSET2D_VALLEY_OVERLAP") == "1" and not getattr(C, "_amset2d_valley", False):
    import numpy as _npv

    _VC = []
    for _part in os.environ.get("AMSET2D_VALLEY_CENTERS", "").split(";"):
        _part = _part.strip()
        if _part:
            _VC.append(_npv.array([float(x) for x in _part.split(",")]))
    if not _VC:
        print("[amset2d] valley_overlap: 没给 AMSET2D_VALLEY_CENTERS，退回正常重叠", file=sys.stderr)
    else:
        _VAL = {"kpts": None}

        def _valley_of(kf):
            d = kf - _npv.rint(kf)
            best, bi = 1e9, 0
            for i, c in enumerate(_VC):
                dd = d - c
                dd -= _npv.rint(dd)
                n = float(_npv.dot(dd, dd))
                if n < best:
                    best, bi = n, i
            return bi

        _orig_cr_v = C.calculate_rate

        def _cr_v(*a, **kw):
            try:
                _VAL["kpts"] = a[0].kpoints
            except Exception:
                pass
            return _orig_cr_v(*a, **kw)

        C.calculate_rate = _cr_v
        _orig_go_v = C._get_overlap

        def _go_v(spin_coeffs, mapping, b_idx, k_idx, band_mask, kpoint_mask):
            out = _orig_go_v(spin_coeffs, mapping, b_idx, k_idx, band_mask, kpoint_mask)
            kp = _VAL.get("kpts")
            if kp is not None and len(out):
                va = _valley_of(kp[k_idx])
                for i in range(len(out)):
                    if _valley_of(kp[kpoint_mask[i]]) != va:
                        out[i] = 0.0
            return out

        C._get_overlap = _go_v
        print("[amset2d] **valley_overlap 已启用**：同谷 1、谷间 0；谷心 %s"
              % [_npv.round(c, 4).tolist() for c in _VC], flush=True)
    C._amset2d_valley = True

