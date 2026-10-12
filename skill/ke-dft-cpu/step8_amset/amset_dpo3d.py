"""
amset_dpo3d.py —— [0036] 三维 AMSET（S8，原版散射核）的 ODP / IVS 运行期插件（不改 AMSET 安装）

AMSET 只有 ADP / POP / IMP / PIE，没有非极性光学声子（光学形变势 ODP）与谷间散射（IVS）。2D 由
step8.4_amset2d/amset2d_plugin.py 补上（[0033]）；本文件把同一套常数形变势模型接到三维 S8：

    1/τ_c(x) = π D_c² / (M ω_c) · [ N·G_c(x + ħω_c) + (N + 1)·G_c(x − ħω_c) ]

公式与维数无关（推导见 step8.4_amset2d/dp_optical.py 文件头；三维即 Conwell 的非极性光学声子散射率），
率由 dp_optical.carrier_rates(full3d=True) 解析给出，走 AMSET 的 basic 接口（每 k 点直接给率，1/s）。
三维与 2D 唯一的区别是谷的几何：k 点按完整三维距离（27 个周期像）归属最近的谷中心。

用法（gen_step10_amset.py 自动接线）：运行目录里有 dp_optical_3d.json（gen 写）+ dp_optical.py，
python -c 入口里 `import amset_dpo3d`（放在其它插件之后：本插件包装散射 worker，子进程按引用导入本模块）。
dp_optical_3d.json：
    {"mechanisms": ["ODP", "IVS"],                    # 本次要注册的机制（settings.yaml 的 scattering_type 里也有）
     "dp_optical": {"electron": {"valleys": [...], "channels": [...]}, "hole": {...},
                    "adp_intravalley": {"electron": bool, "hole": bool}, ...}}
可选：ADP 谷内掩码（开了 IVS 的载流子，ADP 只留 k 与 k+q 在同一个谷的散射，免得与 IVS 重复计数），
与 2D 插件同一做法：替换 AMSET 的 ADP 类（核不变，只乘掩码），worker 子进程里按引用重建。
没有 dp_optical_3d.json -> 什么都不做（打印一行）。
"""

import json
import os
import sys

import numpy as np

import amset
import amset.scattering.calculate as C
from amset.scattering.basic import AbstractBasicScattering as _ABS_BASIC
from amset.scattering.elastic import AcousticDeformationPotentialScattering as _ADP

# spawn 的散射子进程只导入 worker 所在模块：本插件包装 worker 时，其它插件（IR_FIX 的 ADP 星平均）也要在子进程里生效
try:
    import amset_ir_fix  # noqa: F401
except ImportError:
    pass

_TAG = "[amset_dpo3d]"
_SUPPORTED_AMSET = ("0.4.19", "0.5.1")
_REC_PATH = os.environ.get("AMSET_DPO3D_RECORD", "dp_optical_3d.json")
_CARRIERS = ("electron", "hole")
_DPO_KINDS = {"ODP": frozenset({"intra"}), "IVS": frozenset({"inter_star", "inter"})}

if os.path.isfile(_REC_PATH):
    with open(_REC_PATH) as _f:
        _REC = json.load(_f)
else:
    _REC = None
_DPO_REC = (_REC or {}).get("dp_optical") or {}
_MECHS = [m for m in ("ODP", "IVS") if m in ((_REC or {}).get("mechanisms") or [])]
_MASK_ON = any((_DPO_REC.get("adp_intravalley") or {}).get(c) for c in _CARRIERS)

if _REC is not None:
    try:
        import dp_optical as _DPO
    except ImportError:
        raise ImportError("%s %s 开了 ODP/IVS 或 ADP 谷内掩码，但运行目录里没有 dp_optical.py"
                          "（请用 gen_step10_amset.py 重新生成 step8_amset）" % (_TAG, _REC_PATH))
    if amset.__version__ not in _SUPPORTED_AMSET:
        print("%s[WARN] 未在 %s 上核对过本插件（已核对：%s），将继续运行但请自行核对 basic 散射接口。"
              % (_TAG, amset.__version__, "/".join(_SUPPORTED_AMSET)), file=sys.stderr)
    for _n in ("_scattering_mechanisms", "scattering_worker", "AcousticDeformationPotentialScattering"):
        assert hasattr(C, _n), "amset 内部结构已变化：amset.scattering.calculate 缺少 %s" % _n
else:
    _DPO = None


def _split_bands(amset_data, spin):
    """该自旋的能带 -> (电子带, 空穴带) 布尔数组。常数形变势模型按载流子分谷，金属没有意义。"""
    if amset_data.is_metal or not amset_data.vb_idx:
        raise ValueError("%s ODP/IVS/ADP 谷内掩码需要半导体（AMSET 的 vb_idx），本体系被判为金属" % _TAG)
    cb = np.arange(len(amset_data.energies[spin])) > amset_data.vb_idx[spin]
    return cb, ~cb


def _adp_mask_setup(amset_data):
    """ADP 谷内掩码的数据（{"electron"/"hole": 成员谷中心 (n,3) 或 None, "rlat", "rinv"}）；都不开 -> None。"""
    flags = _DPO_REC.get("adp_intravalley") or {}
    if not any(flags.get(c) for c in _CARRIERS):
        return None
    rlat = np.asarray(amset_data.structure.lattice.reciprocal_lattice.matrix, float)
    out = {"rlat": rlat, "rinv": np.linalg.inv(rlat)}
    for car in _CARRIERS:
        vs = (_DPO_REC.get(car) or {}).get("valleys") or []
        out[car] = _DPO.member_table(vs)[0] if (flags.get(car) and vs) else None
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
        print("%s[WARN] ODP/IVS：SOC 能带（旋量带）—— 目标谷态密度把两种自旋都算进去了；"
              "自旋守恒的通道率会偏大到 2 倍，参数要与文献的 SOC 口径一致，结论慎用。" % _TAG, flush=True)
    rlat = amset_data.structure.lattice.reciprocal_lattice.matrix
    spin = amset_data.spins[0]
    for car, sel, sign in zip(_CARRIERS, _split_bands(amset_data, spin), (+1, -1)):
        vs = (_DPO_REC.get(car) or {}).get("valleys") or []
        if not vs or not sel.any():
            continue
        rows = _DPO.valley_report(np.asarray(amset_data.energies[spin])[sel], amset_data.kpoints, rlat, vs, sign,
                                  full3d=True)
        txt = "; ".join("%s 谷底离中心 %s、高于带边 %s" % (
            r["label"], "—" if r["dk_rel"] is None else "%.3f|b1|" % r["dk_rel"],
            "—" if r["dE_eV"] is None else "%.3f eV" % r["dE_eV"]) for r in rows)
        bad = [r["label"] for r in rows if r["dk_rel"] is not None and r["dk_rel"] > VALLEY_DK_WARN]
        used = {lab for h in (_DPO_REC.get(car) or {}).get("channels") or [] for lab in (h["from"], h["to"])}
        bad += [r["label"] for r in rows if r["dk_rel"] is None and r["label"] in used]
        print("%s%s 谷表核对（%s）：%s%s" % (
            _TAG, "[WARN]" if bad else "", car, txt,
            ("  ——%s 的谷底不在给定中心或分不到 k 点：谷表（DP_OPTICAL_VALLEYS_*，S3 晶胞分数坐标）"
             "与 AMSET 的晶胞不同源？" % ",".join(bad)) if bad else ""), flush=True)


def _dpo_rates(amset_data, kinds):
    """{spin: (ndop, ntemp, nbands, nk)} 的率（1/s），按载流子各用自己的谷表与通道；没有通道的态加 RATE_FLOOR。"""
    _valley_check(amset_data)
    rlat = np.asarray(amset_data.structure.lattice.reciprocal_lattice.matrix, float)
    mass = float(amset_data.structure.composition.weight) * _DPO.AMU           # 原胞质量 kg
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
                r[:, sel] = _DPO.carrier_rates(E[sel], kf, rlat, c.get("valleys") or [], c["channels"],
                                               temps, mass, kinds, sign, full3d=True)
        r += _DPO.RATE_FLOOR
        # 按 AMSET 的不可约映射对称化（与 2D 插件同理：最近中心归属在等距边界上的取舍不随对称操作协变）
        ir, i2f = getattr(amset_data, "ir_kpoints_idx", None), getattr(amset_data, "ir_to_full_kpoint_mapping", None)
        if ir is not None and i2f is not None:
            r = r[..., np.asarray(ir)][..., np.asarray(i2f)]
        out[spin] = np.broadcast_to(r[None], (ndop,) + r.shape)
    return out


class _DPOptical3D(_ABS_BASIC):
    """常数形变势光学声子散射的公共部分（AMSET basic 接口）。"""
    required_properties = ()          # 参数走 dp_optical_3d.json，不进 settings 校验
    kinds = frozenset()

    def __init__(self, properties, doping, temperatures, nbands, rates):
        super().__init__(properties, doping, temperatures, nbands)
        self._rates = rates

    @classmethod
    def from_amset_data(cls, materials_properties, amset_data):
        rates = _dpo_rates(amset_data, cls.kinds)
        _DPO_FULL[cls.name] = rates          # FD 窗口外也用解析率（见 _install_fill_keep_dpo）
        return cls(cls.get_properties(materials_properties), amset_data.doping,
                   amset_data.temperatures, cls.get_nbands(amset_data), rates)

    @property
    def rates(self):
        return self._rates


class OpticalDeformation3D(_DPOptical3D):
    """谷内光学形变势（X>X 通道）。"""
    name = "ODP"
    kinds = _DPO_KINDS["ODP"]


class Intervalley3D(_DPOptical3D):
    """谷间散射（X>X′ 同星其它谷、X>Y 到别的星）。"""
    name = "IVS"
    kinds = _DPO_KINDS["IVS"]


# FD 窗口外：AMSET 把每个机制窗口外的率填成窗口内的对数平均；只给一种载流子写了通道时另一种载流子的
# RATE_FLOOR 会把平均拉低几个量级（2D 合成双谷实测单机制迁移率被放大几十倍，V196）。填完后改回解析率。
_DPO_FULL = {}


def _install_fill_keep_dpo():
    import amset.core.data as _AD
    if getattr(_AD.AmsetData.fill_rates_outside_cutoffs, "_az_dpo", False):
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
            print("%s %s：FD 窗口外保留解析率（不用 AMSET 的对数平均填充）" % (_TAG, "/".join(hit)), flush=True)
        return out

    _fill_keep_dpo._az_dpo = True
    # 2D 插件也装同一个补丁（标记 _amset2d_dpo）；两者不会同进程出现，但保持幂等
    _AD.AmsetData.fill_rates_outside_cutoffs = _fill_keep_dpo


class AcousticDeformation3DIntravalley(_ADP):
    """原版三维 ADP 核，开了掩码的载流子只留谷内（k 与 k+q 落在同一个成员谷）。
    核与 AMSET 原版逐字相同（_factor_core = 原版 factor）；amset_ir_fix 的星平均只平均 _factor_core，
    掩码在平均之后按未旋转的 k、q 乘（与 2D 插件同一约定）。没有 _ivm 时与原版逐字节相同。"""
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
        ivm = getattr(self, "_ivm", None)
        if not ivm or kpoint is None:
            return None
        vb = (getattr(self, "vb_idx", None) or {}).get(spin)
        centers = ivm.get("electron" if (vb is not None and band_idx > vb) else "hole")
        if centers is None:
            return None
        q = np.sqrt(norm_q_sq)[:, None] * np.asarray(unit_q)
        kp = q @ ivm["rinv"] + np.asarray(kpoint, float)
        return _DPO.same_member(kpoint, kp, centers, ivm["rlat"], full3d=True)

    def _factor_core(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        return _ADP.factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity)

    def factor(self, unit_q, norm_q_sq, spin, band_idx, kpoint, velocity):
        f = self._factor_core(unit_q, norm_q_sq, spin, band_idx, kpoint, velocity)
        keep = self._iv_keep(unit_q, norm_q_sq, spin, band_idx, kpoint)
        return f if keep is None else f * keep


def worker(*args, **kwargs):
    """散射子进程入口：spawn 按引用反序列化本函数 -> 子进程导入本模块（ADP 谷内掩码类在子进程里生效）。"""
    return _ORIG_WORKER[0](*args, **kwargs)


_ORIG_WORKER = [None]


def _install():
    if _REC is None:
        print("%s 没有 %s：不注册 ODP/IVS（与原版 AMSET 相同）" % (_TAG, _REC_PATH), flush=True)
        return
    for k in _MECHS:     # 开了就必须有对应种类的通道（gen 保证；手写 json 在这里拦）
        if not any(h.get("kind") in _DPO_KINDS[k]
                   for c in _CARRIERS for h in (_DPO_REC.get(c) or {}).get("channels") or []):
            raise ValueError("%s mechanisms 含 %s，但 %s 的 dp_optical 里没有 %s 通道"
                             % (_TAG, k, _REC_PATH, "/".join(sorted(_DPO_KINDS[k]))))
    classes = {"ODP": OpticalDeformation3D, "IVS": Intervalley3D}
    for k in _MECHS:
        C._scattering_mechanisms[k] = classes[k]
    if _MECHS:
        _install_fill_keep_dpo()
    if _MASK_ON:
        C._scattering_mechanisms["ADP"] = AcousticDeformation3DIntravalley
        C.AcousticDeformationPotentialScattering = AcousticDeformation3DIntravalley   # worker 用此名重建 ADP
        if C.scattering_worker is not worker and _ORIG_WORKER[0] is None:
            _ORIG_WORKER[0] = C.scattering_worker
            C.scattering_worker = worker
    ivf = _DPO_REC.get("adp_intravalley") or {}
    print("%s amset %s | 三维常数形变势 %s：%s | ADP 谷内掩码：电子 %s，空穴 %s"
          % (_TAG, amset.__version__, ",".join(_MECHS) or "（无，只掩码）", _DPO.describe(_DPO_REC),
             "开" if ivf.get("electron") else "关", "开" if ivf.get("hole") else "关"), flush=True)


if not getattr(C, "_az_dpo3d", False):
    C._az_dpo3d = True
    _install()
