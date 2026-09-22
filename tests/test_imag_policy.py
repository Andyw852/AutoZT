#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tests/test_imag_policy.py —— 虚频判据唯一真源 imag_policy 的自测（2026-09-22）。

覆盖：真实 MoS₂ 数据快照、6 组合成用例、|q|_frac 边界、3D 输入、以及
"单一裁判"源码守卫（S5.1 / lattice_kappa / mlff driver 不得再自带阈值）。
"""
import os
import re
import sys
import math

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "skill", "_common"))
import imag_policy as ip          # noqa: E402

FAILED = []


def check(name, cond, detail=""):
    if cond:
        print("  ✓ %s%s" % (name, ("  ← " + str(detail)) if detail else ""))
    else:
        print("  ✗ %s   %s" % (name, detail))
        FAILED.append(name)


# ---------------------------------------------------------------------------
# 真实数据快照：MoS₂ 的 S5.1 生产输出 step5_phonon_plot/band_nonac.dat
#   （该文件由 step5_fc/phono3py/fc2.hdf5 算出；q_red = 原胞倒格基分数坐标，kz=0）
#   最小那一点的精确值取自 production band_nonac.yaml: -0.055377 THz
# ---------------------------------------------------------------------------
_MOS2_SOFT = [
    # q_red(原胞基矢), 9 支频率(THz, 升序)
    ([0.0, 0.0555556, 0.0555556], [0.1180, 0.8665, 1.3730, 8.5821, 8.6614,
                                   11.5078, 11.5655, 12.1855, 14.0700]),
    ([0.0, 0.0444444, 0.0444444], [0.0306, 0.6952, 1.1011, 8.5722, 8.6297,
                                   11.5301, 11.5694, 12.1757, 14.0986]),
    ([0.0, 0.0333333, 0.0333333], [-0.055377, 0.5225, 0.8274, 8.5645, 8.6045,
                                   11.5480, 11.5725, 12.1675, 14.1209]),
    ([0.0, 0.0222222, 0.0222222], [-0.050294, 0.3489, 0.5524, 8.5588, 8.5862,
                                   11.5612, 11.5746, 12.1613, 14.1370]),
    ([0.0, 0.0111111, 0.0111111], [-0.028460, 0.1746, 0.2764, 8.5554, 8.5752,
                                   11.5692, 11.5759, 12.1574, 14.1467]),
    ([0.0, 0.0, 0.0],             [0.0, 0.0, 0.0, 8.5543, 8.5714,
                                   11.5719, 11.5763, 12.1561, 14.1499]),
]
# 真空轴在【原胞基矢】里是第 0 个（生产 primitive_matrix 把 25 Å 真空放第 0 基矢）
_MOS2_VAC = 0


def _cls(freqs, q, is2d=True, vac=_MOS2_VAC):
    return ip.classify_imag([freqs], [q], is2d, vac)


def _acoustic(nu, extra=0.5):
    """最低 3 支里放一个负频，其余正当频率。"""
    return [nu, extra, extra + 0.4, 8.0, 8.1, 11.0, 11.1, 12.0, 14.0]


def _optical(nu):
    return [0.1, 0.2, 0.3, nu, 8.0, 8.1, 11.0, 11.1, 12.0, 14.0]


def main():
    print("[1. 真实数据：MoS₂ fc2 的软模区快照]")
    qs = [q for q, _ in _MOS2_SOFT]
    fs = [f for _, f in _MOS2_SOFT]
    r = ip.classify_imag(fs, qs, True, _MOS2_VAC)
    check("verdict = warn", r["verdict"] == "warn", r["verdict"])
    check("imag_class = near_gamma_acoustic", r["imag_class"] == "near_gamma_acoustic",
          r["imag_class"])
    check("min_freq_THz = -0.055377", abs(r["min_freq_THz"] + 0.055377) < 1e-9,
          r["min_freq_THz"])
    check("min_freq_cm1 ≈ -1.847", abs(r["min_freq_cm1"] + 1.847) < 0.01,
          round(r["min_freq_cm1"], 4))
    qn = ip.q_norm_frac([0.0, 1 / 30.0, 1 / 30.0], True, _MOS2_VAC)
    check("|q|_frac(1/30,1/30) ≈ 0.04714", abs(qn - 0.0471404) < 1e-6, round(qn, 6))
    check("q_norm_at_min ≈ 0.04714", abs(r["q_norm_at_min"] - qn) < 1e-5,
          round(r["q_norm_at_min"], 6))
    check("branch_at_min = 0（声学）", r["branch_at_min"] == 0, r["branch_at_min"])
    check("q_at_min_frac = (0, 1/30, 1/30)",
          r["q_at_min_frac"] is not None and abs(r["q_at_min_frac"][2] - 1 / 30.0) < 1e-7,
          r["q_at_min_frac"])
    check("n_neg_qpoints = 3", r["n_neg_qpoints"] == 3, r["n_neg_qpoints"])
    check("policy_version 有值", bool(r["policy_version"]), r["policy_version"])
    check("thresholds 三键齐全",
          all(k in r["thresholds"] for k in ("IMAG_THR", "IMAG_THR_STRICT", "IMAG_QGAMMA")),
          list(r["thresholds"].keys())[:3])

    print("[2. 合成：近 Γ 声学支]")
    check("近 Γ 声学 -0.03 → pass", _cls(_acoustic(-0.03), [0, 0.01, 0.01])["verdict"] == "pass")
    r2 = _cls(_acoustic(-0.03), [0, 0.01, 0.01])
    check("近 Γ 声学 -0.03 → imag_class noise", r2["imag_class"] == "noise", r2["imag_class"])
    check("近 Γ 声学 -0.20 → fail", _cls(_acoustic(-0.20), [0, 0.01, 0.01])["verdict"] == "fail")
    r2 = _cls(_acoustic(-0.20), [0, 0.01, 0.01])
    check("近 Γ 声学 -0.20 → near_gamma_large",
          r2["imag_class"] == "near_gamma_large", r2["imag_class"])

    print("[3. 合成：近 Γ 光学支]")
    r3 = _cls(_optical(-0.06), [0, 0.01, 0.01])
    check("近 Γ 光学 -0.06 → fail", r3["verdict"] == "fail")
    check("近 Γ 光学 -0.06 → optical", r3["imag_class"] == "optical", r3["imag_class"])

    print("[4. 合成：近 Γ 区外]")
    check("区外 -0.03 → pass", _cls(_acoustic(-0.03), [0, 0.2, 0.2])["verdict"] == "pass")
    r4 = _cls(_acoustic(-0.06), [0, 0.2, 0.2])
    check("区外 -0.06 → fail", r4["verdict"] == "fail")
    check("区外 -0.06 → off_gamma", r4["imag_class"] == "off_gamma", r4["imag_class"])

    print("[5. |q|_frac 边界 + grace 带（IMAG_QGAMMA_GRACE=1.2，用户 2026-09-22 批准）]")
    qb = [0.0, 0.05, 0.0]
    nb = ip.q_norm_frac(qb, True, _MOS2_VAC)
    check("|q|_frac(0.05) >= 0.05", nb >= 0.05, nb)
    rb = _cls(_acoustic(-0.20), qb)
    check("|q|=0.05 落在 grace 带内 → near_gamma_large（不再 off_gamma）",
          rb["imag_class"] == "near_gamma_large" and rb["verdict"] == "fail", rb["imag_class"])
    rb2 = _cls(_acoustic(-0.20), [0.0, 0.049, 0.0])
    check("|q|=0.049 → near_gamma_large", rb2["imag_class"] == "near_gamma_large",
          rb2["imag_class"])
    # 真实 MoS₂ 的边界点：ν=−0.0503 @ q=(0, 0.03667, 0.03667)，|q|=0.0519（只超 Petretto 窗口 3.8%）
    _qreal = [0.0, 0.0366667, 0.0366667]
    _nreal = ip.q_norm_frac(_qreal, True, _MOS2_VAC)
    check("真实边界点 |q|≈0.0519", abs(_nreal - 0.0518513) < 1e-4, round(_nreal, 5))
    rr = _cls(_acoustic(-0.0503), _qreal)
    check("真实边界点 → near_gamma_acoustic / warn（grace 生效）",
          rr["imag_class"] == "near_gamma_acoustic" and rr["verdict"] == "warn",
          "%s/%s" % (rr["verdict"], rr["imag_class"]))
    r0 = ip.classify_imag([_acoustic(-0.0503)], [_qreal], True, _MOS2_VAC,
                          {"IMAG_QGAMMA_GRACE": 1.0})
    check("grace=1.0 → 同一点退回 off_gamma / fail",
          r0["imag_class"] == "off_gamma" and r0["verdict"] == "fail",
          "%s/%s" % (r0["verdict"], r0["imag_class"]))
    rb3 = _cls(_acoustic(-0.06), [0.0, 0.065, 0.0])
    check("|q|=0.065（grace 带外）→ off_gamma / fail",
          rb3["imag_class"] == "off_gamma" and rb3["verdict"] == "fail", rb3["imag_class"])

    print("[6. 3D 输入：位置规则相同、不触发 ZA]")
    r6 = ip.classify_imag([_acoustic(-0.06)], [[0.02, 0.02, 0.02]], False, None)
    check("3D 近 Γ 声学 -0.06 → warn", r6["verdict"] == "warn", r6["verdict"])
    check("3D 用三分量（|q|=0.0346 < 0.05）",
          abs(ip.q_norm_frac([0.02, 0.02, 0.02], False, None) - 0.034641) < 1e-6)
    check("3D 无 ZA → za_verdict pass", ip.za_verdict(None) == "pass")
    check("combine(pass, warn) = warn", ip.combine_verdict("pass", "warn") == "warn")
    check("combine(warn, needs_review) = needs_review",
          ip.combine_verdict("warn", "needs_review") == "needs_review")
    check("combine(needs_review, fail) = fail",
          ip.combine_verdict("needs_review", "fail") == "fail")
    check("is_stable: warn 放行 / needs_review 放行 / fail 挡",
          ip.is_stable("warn") and ip.is_stable("needs_review") and not ip.is_stable("fail"))

    print("[7. 一致性 + 单一裁判源码守卫]")
    # 同一份真实快照 → 唯一期望值（四条路径都应拿到这个）
    expect = ("warn", "near_gamma_acoustic", round(r["min_freq_THz"], 6))
    got = []
    for _ in range(4):
        rr = ip.classify_imag(fs, qs, True, _MOS2_VAC)
        got.append((rr["verdict"], rr["imag_class"], round(rr["min_freq_THz"], 6)))
    check("重复调用给出同一 verdict/class/min", all(g == expect for g in got), expect)

    def _src(rel):
        with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
            return fh.read()

    s51 = _src("skill/kl-dft-cpu/gen_step5.1_plot_phonon.py")
    lk = _src("skill/kl-dft-cpu/lattice_kappa.py")
    m_cpu = _src("skill/phonon-mlff-cpu/phonon_fit_driver.py")
    m_gpu = _src("skill/phonon-mlff-gpu/phonon_fit_driver.py")
    kfb = _src("skill/kl-dft-cpu/kl_fc_backends.py")
    check("S5.1 不再有 IMAG_TOL", "IMAG_TOL" not in s51)
    check("S5.1 引用 imag_policy", "imag_policy" in s51)
    check("lattice_kappa 不再有 imag_thr = 0.5", not re.search(r"imag_thr\s*=\s*0\.5", lk))
    check("lattice_kappa 引用 imag_policy", "imag_policy" in lk)
    check("mlff-cpu driver 不再硬编码 -0.10", "-0.10" not in m_cpu)
    check("mlff-gpu driver 不再硬编码 -0.10", "-0.10" not in m_gpu)
    check("mlff 两个 driver 都引用 imag_policy",
          "imag_policy" in m_cpu and "imag_policy" in m_gpu)
    check("S5（kl_fc_backends）调用 classify_imag", "classify_imag" in kfb)

    print()
    if FAILED:
        print("FAILED: %d 项 -> %s" % (len(FAILED), ", ".join(FAILED)))
        sys.exit(1)
    print("ALL PASS")


def test_imag_policy():
    """pytest 入口（main() 内部失败会 sys.exit(1)）。"""
    main()


if __name__ == "__main__":
    main()
