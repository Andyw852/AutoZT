"""
amset_pie3d.py —— [0037] 三维 AMSET（S8）的压电散射（PIE）核（运行期插件，不改 AMSET 安装）

原版 AMSET 0.4.19/0.5.1 的 PiezoelectricScattering 与教科书三维压电核（Ridley; Rode; Kaasbjerg 2013 的
三维对应式）有三处口径差：
  ① 屏蔽：原版 factor ∝ 1/(q² + β²)，相当于只屏蔽一次；压电势是长程库仑势，屏蔽后耦合
     φ(q) ∝ 1/ε(q) = q²/(q² + β²)，|M|² ∝ (1/q²)·q⁴/(q² + β²)² = q²/(q² + β²)²（q≫β 时两者一致，
     q→0 原版趋 1/β²、正确值趋 0 —— 原版高估小角散射，屏蔽越强偏得越多）；
  ② 介电：原版用 h = ε∞⁻¹·e（且屏蔽也用 ε∞）。e 取 VASP 的弛豫离子总张量（声学声子频率远低于光学支，
     离子跟得上），静态压电势应配**静态** ε0（与 2D 插件 [0031] 同一口径）；
  ③ 方向：纵场只沿 q̂，正确的投影是 q̂·e:(q̂v)/(q̂·ε0·q̂)；原版 (ε⁻¹q̂)·e:(q̂v) 只在 ε 各向同性时相同。
  原版 free_carrier_screening=false 时屏蔽填 1e-12（≈不屏蔽）；本插件的屏蔽单独由 pie_3d.json 控制
  （默认开），不影响 POP（全局 free_carrier_screening 不动）。

本插件的核（Hartree 原子单位，高斯制，e = ħ = 1；prefactor 与原版相同 = k_B T）：

    factor = (4π)² Σ_b [q̂·e:(q̂ v_b)]² / (c_b (q̂·ε0·q̂)²) · q² / (q² + κ/(q̂·ε0·q̂))²
    κ = 4π ∂n/∂μ（自由载流子 Thomas–Fermi，与 AMSET calculate_inverse_screening_length_sq 同一积分，ε=1）

b 为 Christoffel 方程的三支声学模（c_b 刚度、v_b 极化）。ε 各向同性、κ→0 时与原版逐点相同
（test_pie3d.py 用它锁住绝对归一化）。

用法（gen_step10_amset.py 自动接线）：运行目录里有 pie_3d.json（gen 写）时，python -c 入口里
`import amset_pie3d` 把注册表里的 "PIE" 换成 Piezoelectric3D；elastic 散射体按实例 pickle 进 spawn 子进程，
子进程按类引用自动导入本模块，不用包装 worker。
pie_3d.json: {"kernel": "az", "screening": true, ...}；kernel="amset" 或没有该文件 -> 什么都不做（原版 PIE）。
"""

import json
import os
import sys

import numpy as np

import amset
import amset.scattering.calculate as C
from amset.constants import boltzmann_au, coulomb_to_au, m_to_bohr
from amset.scattering.elastic import PiezoelectricScattering as _PIE
from amset.scattering.elastic import get_christoffel_tensors

from amset.electronic_structure.fd import fd

_TAG = "[amset_pie3d]"
_SUPPORTED_AMSET = ("0.4.19", "0.5.1")
_REC_PATH = os.environ.get("AMSET_PIE3D_RECORD", "pie_3d.json")
if os.path.isfile(_REC_PATH):
    with open(_REC_PATH) as _f:
        _REC = json.load(_f)
else:
    _REC = None


def screening_kappa(amset_data):
    """κ = 4π ∂n/∂μ（bohr⁻²·ε 单位，形状 (ndop, ntemp)）；∂n/∂μ 用 AMSET 自己的总态密度积分。"""
    trap = getattr(np, "trapezoid", None) or np.trapz
    tdos, en = amset_data.dos.tdos, amset_data.dos.energies
    vol = amset_data.structure.volume
    out = np.zeros(amset_data.fermi_levels.shape)
    for n, t in np.ndindex(out.shape):
        kt = amset_data.temperatures[t] * boltzmann_au
        f = fd(en, amset_data.fermi_levels[n, t], kt)
        out[n, t] = 4 * np.pi * trap(tdos * f * (1 - f), x=en) / (kt * vol)
    return out


def pie_factor(e_au, eps0, elastic_au, kappa, unit_q, norm_q_sq):
    """三维压电核（见文件头），返回 (ndop, ntemp, nq)。e_au: 3×3×3（a.u.），eps0: 3×3 相对介电，
    elastic_au: 3×3×3×3（a.u.），kappa: (ndop, ntemp)。"""
    g = get_christoffel_tensors(elastic_au, unit_q)
    w, vec = np.linalg.eigh(g)                       # w (n,3)，vec[:, :, b] 为第 b 支极化
    eq = np.einsum("ij,ni,nj->n", eps0, unit_q, unit_q)
    s = np.zeros(len(unit_q))
    for b in range(3):
        cpl = np.einsum("ijk,ni,nj,nk->n", e_au, unit_q, unit_q, vec[:, :, b])
        s = s + cpl ** 2 / w[:, b]
    amp = (4 * np.pi) ** 2 * s / eq ** 2
    qq = norm_q_sq[None, None, :]
    return amp[None, None] * qq / (qq + kappa[..., None] / eq[None, None]) ** 2


class Piezoelectric3D(_PIE):
    """三维压电散射：静态 ε0、纵场 q̂ 投影、屏蔽 q²/(q²+β²)²（见文件头）。"""
    name = "PIE"
    required_properties = ("piezoelectric_constant", "elastic_constant", "static_dielectric")

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        # 拷贝再换单位：原版会就地改 settings 里的数组
        e_au = np.array(materials_properties["piezoelectric_constant"], float) * coulomb_to_au / m_to_bohr ** 2
        eps0 = np.array(materials_properties["static_dielectric"], float)
        if eps0.shape != (3, 3) or np.any(np.linalg.eigvalsh(eps0) <= 0):
            raise ValueError("%s static_dielectric 必须是正定 3×3 张量（拿到 %s）" % (_TAG, eps0.tolist()))
        if (_REC or {}).get("screening", True):
            kappa = screening_kappa(amset_data)
        else:
            kappa = np.zeros(amset_data.fermi_levels.shape)
        obj = cls(cls.get_properties(materials_properties), amset_data.doping,
                  amset_data.temperatures, cls.get_nbands(amset_data), e_au, kappa)
        obj.eps0 = eps0
        return obj

    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        # piezoelectric_constant 槽位存 e（a.u.），inverse_screening_length_sq 槽位存 κ
        return pie_factor(self.piezoelectric_constant, self.eps0, self.elastic_constant,
                          self.inverse_screening_length_sq, unit_q, norm_q_sq)


def _install():
    if _REC is None or _REC.get("kernel", "az") != "az":
        print("%s 没有 %s 或 kernel=amset：PIE 用原版 AMSET 核" % (_TAG, _REC_PATH), flush=True)
        return
    if amset.__version__ not in _SUPPORTED_AMSET:
        print("%s[WARN] 未在 %s 上核对过本插件（已核对：%s），继续运行但请自行核对 elastic 散射接口。"
              % (_TAG, amset.__version__, "/".join(_SUPPORTED_AMSET)), file=sys.stderr)
    assert hasattr(C, "_scattering_mechanisms"), "amset 内部结构已变化：缺 _scattering_mechanisms"
    C._scattering_mechanisms["PIE"] = Piezoelectric3D
    print("%s amset %s | PIE = 三维压电核（静态 ε0、q̂ 纵场投影、屏蔽 %s）| |e|max = %s C/m²"
          % (_TAG, amset.__version__, "开" if _REC.get("screening", True) else "关",
             _REC.get("e_max_C_m2")), flush=True)


if not getattr(C, "_az_pie3d", False):
    C._az_pie3d = True
    _install()
