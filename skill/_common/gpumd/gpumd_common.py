#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gpumd_common.py —— kl-gpumd-3090 的公共工具（gen 侧 + 作业侧共用）。

约定：gen 脚本把该步的全部参数写进 <step_dir>/gpumd_params.json，作业侧只读它，
不在计算节点重新解析 step.conf（避免三层合并在两个地方实现）。
"""
import json
import os
import re
import sys
from pathlib import Path

PARAMS_JSON = "gpumd_params.json"
METHOD_FILE = "workflow_method.txt"
DEFAULT_CONDA_SH = "/home/user_3090/miniconda3/etc/profile.d/conda.sh"
DEFAULT_CONDA_ENV = "wc"


def save_params(step_dir, params):
    Path(step_dir, PARAMS_JSON).write_text(
        json.dumps(params, indent=2, ensure_ascii=False), encoding="utf-8")
    return Path(step_dir, PARAMS_JSON)


def load_params(step_dir="."):
    return json.loads(Path(step_dir, PARAMS_JSON).read_text(encoding="utf-8"))


def write_method(path, **kv):
    lines = ["%s=%s" % (k, v) for k, v in kv.items()]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary(step_dir, name, payload):
    Path(step_dir, name).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print("[summary] %s" % json.dumps(payload, ensure_ascii=False)[:800])


def write_submit(tpl_path, out_path, subs):
    text = Path(tpl_path).read_text(encoding="utf-8")
    for k, v in subs.items():
        text = text.replace("{{%s}}" % k, str(v))
    left = re.findall(r"\{\{([A-Z_0-9]+)\}\}", text)
    if left:
        sys.exit("[ERROR] 模板 %s 还有未填占位符：%s"
                 % (Path(tpl_path).name, ", ".join(sorted(set(left)))))
    Path(out_path).write_text(text, encoding="utf-8", newline="\n")
    os.chmod(str(out_path), 0o755)


def resolve_submit(here, kind):
    here = Path(here)
    for cand in (here / ("%s.tpl" % kind), here.parent / ("%s.tpl" % kind)):
        if cand.is_file():
            return cand
    sys.exit("[ERROR] 找不到 %s.tpl（gen_need 里漏了？）" % kind)


def new_jobname(cwd, label):
    return "%s-kgp-%s" % (Path(cwd).name, label)


def extxyz_frame(lattice, symbols, positions, forces, energy, virial=None,
                 config_type=""):
    """一行 extxyz：GPUMD 的 nep 从 Properties 里的 force 与 energy= 读取标签。"""
    lat = " ".join("%.10f" % x for row in lattice for x in row)
    props = "species:S:1:pos:R:3:force:R:3"
    head = 'Lattice="%s" Properties=%s energy=%.10f' % (lat, props, float(energy))
    if virial is not None:
        head += ' virial="%s"' % " ".join("%.10f" % x for x in virial)
    if config_type:
        head += " config_type=%s" % config_type
    lines = [str(len(symbols)), head]
    for s, p, f in zip(symbols, positions, forces):
        lines.append("%-3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f"
                     % (s, p[0], p[1], p[2], f[0], f[1], f[2]))
    return "\n".join(lines) + "\n"


def parse_nep_log(text):
    """从 nep 的训练输出里抽 RMSE（不同版本措辞略有差异，尽量宽松）。"""
    out = {}
    pat = re.compile(
        r"RMSE\s*\(?\s*(energy|force|virial|total)?\s*\)?\s*[:=]\s*"
        r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")
    for m in pat.finditer(text):
        try:
            out[m.group(1) or "total"] = float(m.group(2))
        except ValueError:
            continue
    gen = re.findall(r"Generation\s+([0-9]+)", text)
    if gen:
        out["generation"] = int(gen[-1])
    return out


def parse_loss_out(path):
    """loss.out（NEP 非电荷）：10 列 = generation, loss_total, loss_L1, loss_L2,
    rmse_energy_train, rmse_force_train, rmse_virial_train, rmse_energy_test,
    rmse_force_test, rmse_virial_test。返回最后一行的指标。"""
    import numpy as np
    a = np.loadtxt(path)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    last = a[-1].astype(float)
    out = dict(generation=int(last[0]))
    names = ["rmse_energy_train", "rmse_force_train", "rmse_virial_train",
             "rmse_energy_test", "rmse_force_test", "rmse_virial_test"]
    for i, n in enumerate(names):
        idx = 4 + i
        out[n] = float(last[idx]) if last.size > idx else None
    return out


def force_file_rms(path):
    """force_*.out：6 列 = NEP 力(fx fy fz) + DFT 力(fx fy fz)。返回 DFT 力逐分量 RMS (eV/Å)。"""
    import numpy as np
    a = np.loadtxt(path)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.shape[1] < 6:
        return None
    dft = a[:, 3:6]
    return float(np.sqrt(np.mean(dft ** 2)))


def parse_kappa_out(path):
    """kappa.out：5 列 (κ_μx^in, κ_μx^out, κ_μy^in, κ_μy^out, κ_μz^tot)，相对**驱动方向 μ**。
    返回 kmux=κ_μx=col1+col2、kmuy=κ_μy=col3+col4、kmuz=κ_μz=col5。
    对角元：x 驱动取 kmux、y 驱动取 kmuy、z 驱动取 kmuz（各跑一次）。"""
    import numpy as np
    a = np.loadtxt(path)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.shape[1] < 5:
        raise SystemExit("[ERROR] kappa.out 少于 5 列：%s" % (a.shape,))
    return dict(kmux=a[:, 0] + a[:, 1], kmuy=a[:, 2] + a[:, 3], kmuz=a[:, 4], raw=a)


def block_average(arr, skip_frac=0.5):
    """丢掉前 skip_frac 的瞬态，返回 (mean, stderr)。"""
    import numpy as np
    tail = np.asarray(arr[int(len(arr) * skip_frac):], dtype=float)
    if tail.size == 0:
        return float("nan"), float("nan")
    blocks = [b for b in np.array_split(tail, 5) if b.size]
    means = np.array([b.mean() for b in blocks])
    return float(tail.mean()), float(means.std(ddof=1) / max(1.0, np.sqrt(len(means))))


def parse_emd(path):
    """hac.out（EMD/Green-Kubo）：col1=τ(ps)，col7-11=κ_x^in,κ_x^out,κ_y^in,κ_y^out,κ_z^tot。
    所以一次 EMD 运行就能拿到三个对角元：κ_xx=col7+8, κ_yy=col9+10, κ_zz=col11。"""
    import numpy as np
    a = np.loadtxt(path)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.shape[1] < 11:
        raise SystemExit("[ERROR] hac.out 少于 11 列：%s" % (a.shape,))
    return dict(tau_ps=a[:, 0], kxx=a[:, 6] + a[:, 7], kyy=a[:, 8] + a[:, 9],
                kzz=a[:, 10], raw=a)


def parse_shc(path):
    """shc.out：带 # 头的两段数据（先 time/ki/ko 的相关段，再 omega_THz/shc_i/shc_o 谱段）。
    返回 (omega_THz, kappa_spec) = (ω, shc_i+shc_o)。"""
    import numpy as np
    ncorr = None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if "num_correlation_rows" in line:
                    ncorr = int(line.split()[-1])
                continue
            parts = line.split()
            if len(parts) >= 3:
                rows.append([float(x) for x in parts[:3]])
    a = np.asarray(rows)
    if a.size == 0:
        raise SystemExit("[ERROR] shc.out 没有数据行")
    spec = a[ncorr:] if (ncorr is not None and ncorr < len(a)) else a
    return spec[:, 0], spec[:, 1] + spec[:, 2]


def parse_dos(path):
    """dos.out：宽松解析（取前两列 = ω(THz), DOS）。"""
    import numpy as np
    a = np.loadtxt(path, comments="#")
    if a.ndim == 1:
        a = a.reshape(1, -1)
    return a[:, 0], a[:, 1]


def quantum_factor_cv(omega_THz, T, hbar=1.054571817e-34, kB=1.380649e-23):
    """谱量子修正因子 f(x) = x^2 e^x /(e^x - 1)^2, x = hbar*omega/(kB*T)。
    x->0 时 f->1（经典极限），x>>1 时按 e^-x 衰减。"""
    import numpy as np
    w = 2.0 * np.pi * np.asarray(omega_THz, dtype=float) * 1e12   # THz -> rad/s
    x = hbar * w / (kB * float(T))
    x = np.clip(x, 1e-12, 700.0)
    ex = np.exp(x)
    return x ** 2 * ex / (ex - 1.0) ** 2


def spectral_quantum_corrected_kappa(omega_THz, kappa_spec, T):
    """用谱 κ(ω) 做严格量子修正：κ_q(T) = ∫κ_cl(ω) f(ω,T)dω；同时返回经典积分。"""
    import numpy as np
    om = np.asarray(omega_THz, dtype=float)
    ks = np.asarray(kappa_spec, dtype=float)
    domega = np.gradient(om)
    k_cl = float(np.sum(ks * domega))
    k_q = float(np.sum(ks * quantum_factor_cv(om, T) * domega))
    return dict(kappa_classical_spectral=k_cl, kappa_quantum_spectral=k_q,
                ratio=(k_q / k_cl if k_cl else float("nan")))


def cv_quantum_over_classical(T, theta_debye):
    """简版量子修正因子 c_v^q/c_v^cl（Debye 近似）；θ_D 未给时返回 1（不修正）。"""
    import numpy as np
    if not theta_debye or float(theta_debye) <= 0:
        return 1.0
    x = float(theta_debye) / float(T)
    if x > 200:
        return 0.0
    trapz = getattr(np, "trapezoid", None) or np.trapz
    n = 400
    t = np.linspace(1e-9, x, n)
    integ = trapz(t ** 3 / (np.exp(t) - 1.0), t)
    d3 = 3.0 / x ** 3 * integ
    cv_q = 3.0 * (4.0 * d3 - 3.0 * x / (np.exp(x) - 1.0))
    return float(max(0.0, cv_q) / 3.0)
