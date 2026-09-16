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
from amset.constants import angstrom_to_bohr as A2B, boltzmann_au, coulomb_to_au, m_to_bohr
from amset.electronic_structure.fd import fd
from amset.scattering.inelastic import PolarOpticalScattering
from amset.scattering.elastic import (
    AcousticDeformationPotentialScattering as _ADP,
    IonizedImpurityScattering as _IMP,
    PiezoelectricScattering as _PIE,
    get_christoffel_tensors,
)

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


_QZ_ONLY_SEEN = [False]


def _note_qz_only(n):
    """只报一次：有多少 q 点的面内投影恰好为 0（这些点用固定面内方向兜底）。"""
    if not _QZ_ONLY_SEEN[0]:
        _QZ_ONLY_SEEN[0] = True
        print("[amset2d] 提示：本次出现 %d 个 q∥=0（q 沿真空方向）的采样点，"
              "已用固定面内方向兜底（二维积分里是零测度）" % n, flush=True)


# ---------------- ADP ----------------
class AcousticDeformation2D(_ADP):
    """只保留面内传播、面内极化的 LA/TA 支；弹性常数用原始 slab 值（= C_2D/c）。"""
    name = "ADP"

    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
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
        return f[None, None] * np.ones(self.fermi_levels.shape + norm_q_sq.shape)


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
    """二维面内压电：F = c²(2π)² Σ_b (q̂q̂:e·v_b)²/w_b · q²/(q·ε∞(q) + q_TF)²。

    只处理面内压电（Janus 结构的面外偶极层在这里被忽略）；e 与 w 都用原始 slab 值。
    """
    name = "PIE"

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        e_au = np.array(materials_properties["piezoelectric_constant"], float)
        e_au = e_au * coulomb_to_au / m_to_bohr ** 2      # C/m² -> a.u.（原版会就地改 settings）
        obj = super().from_amset_data(materials_properties, amset_data)
        obj.e_raw = e_au
        obj.r_inf = _r_tensor("eps_inf_slab")
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
        a_i = _q_eps(qv, qn, self.r_inf)
        return (_C ** 2 * (2 * np.pi) ** 2 * (s * qn ** 2)[None, None]
                / (a_i[None, None] + self.q_tf[..., None]) ** 2)


# ---------------- 注册 ----------------
_CLASSES = {"ADP": AcousticDeformation2D, "POP": PolarOptical2D,
            "IMP": IonizedImpurity2D, "PIE": Piezoelectric2D}
_ACTIVE = [k for k in ("ADP", "POP", "IMP", "PIE") if k in _MECHS]
for _k in _ACTIVE:
    C._scattering_mechanisms[_k] = _CLASSES[_k]
if "ADP" in _MECHS:
    C.AcousticDeformationPotentialScattering = AcousticDeformation2D   # worker 用此名重建 ADP

_BANNER = ("[amset2d] amset %s | c=%.3f A | r_inf=%s | mechanisms=%s"
           % (amset.__version__, _C / A2B, np.round(_r_tensor("eps_inf_slab").tolist(), 3),
              ",".join(_ACTIVE)))
print(_BANNER, flush=True)

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
