# -*- coding: utf-8 -*-
"""phonon_stability.py —— 声子稳定性判据的【公共编排单一真源】（2026-09-27）。

## 为什么有它
2 阶声子的「虚频判据 + 2D ZA 二次性」此前在 kl-dft-cpu(S5_fc)、phonon-dft-cpu(S3)、
phonon-mlff-cpu/gpu、fc-fit 各处各写一份，改一处漏一处（kl_fc_backends.py 顶部 TODO
早记了这事）。本模块把【编排】收成一份：网格选择 / mesh+band 采样 / imag_policy 判定 /
ZA 二次性；阈值与分类规则仍在 _common/imag_policy.py 与 _common/za_2d.py（各自唯一真源）。

调用方需保证 ph 是已拟合 fc2 的 phonopy.Phonopy 对象；本模块不碰拟合。
IMAG_* / ZA_* 阈值从 cfg（dict / StepConf，支持 .get 或 []）读取，缺省走模块默认。

## 编排
  mesh_numbers(ph, is2d, cart_vac_axis, mesh_n)  → (mesh, vax_prim)
  imag_samples(ph)                               → (freqs, qpoints)  # mesh + band 合并
  classify(ph, is2d, cart_vac_axis, cfg)         → (imag_dict, vax_prim)
  za_check(cfg, ph, cart_vac_axis, is2d)         → dict | None（2D 的 ZA 二次性）
  finalize(imag, za)                             → (stability_verdict, stable)
"""
import imag_policy
from za_2d import (ZA_P_RANGE, ZA_ZFRAC_MIN, ZA_MARGIN_MIN, ZA_R2_MIN,
                   ZA_SYMPREC_DEFAULT, ZA_SYMPREC_RELAX,
                   za_power_law, za_power_law_eig,
                   inplane_qdirs, vacuum_axis_in_primitive,
                   cell_symmetry, relaxed_symprec, clone_phonopy)

IMAG_MESH_N = 60        # 虚频门禁面内网格（显式整数；float mesh 会漏近 Γ 软模）


def _cfg_get(cfg, key, default=None):
    """cfg 可能是 dict 或 StepConf：优先 []，没有 .get 就退化。"""
    if cfg is None:
        return default
    try:
        v = cfg.get(key, default)
    except AttributeError:
        try:
            v = cfg[key]
        except (KeyError, TypeError):
            return default
    return default if v is None else v


def mesh_numbers(ph, is2d, cart_vac_axis, mesh_n=IMAG_MESH_N):
    """虚频门禁的显式整数网格：2D 真空轴 1、面内 mesh_n（Γ 中心）。

    ★ 不要用 float mesh：phonopy 把 float 当【面间距长度】(N_i=nint(l/d_i)) 并强制 Γ 中心，
      实测 mesh=60.0 → [2,22,22]，面内最近 q 只有 0.105 Å⁻¹，会漏掉近 Γ 的软模。
      3D 返回 (None, None)，由调用方用 run_mesh(mesh=60.0)（3D 绝不能把某轴设 1）。
    """
    if not is2d:
        return None, None
    vax_prim = vacuum_axis_in_primitive(ph, cart_vac_axis)
    mesh = [int(mesh_n)] * 3
    mesh[int(vax_prim)] = 1
    return mesh, int(vax_prim)


def imag_samples(ph):
    """mesh + band 上全部 q 点的频率与分数坐标（一一对应）。缺失项按空跳过。"""
    freqs, qpts = [], []
    try:
        md = ph.get_mesh_dict()
        freqs += [[float(x) for x in row] for row in md["frequencies"]]
        qpts += [[float(x) for x in q] for q in md["qpoints"]]
    except Exception:                                   # noqa: BLE001
        pass
    try:
        bd = ph.get_band_structure_dict()
        freqs += [[float(x) for x in row]
                  for seg in (bd.get("frequencies") or []) for row in seg]
        qpts += [[float(x) for x in row]
                 for seg in (bd.get("qpoints") or []) for row in seg]
    except Exception:                                   # noqa: BLE001
        pass
    return freqs, qpts


def classify(ph, is2d, cart_vac_axis, cfg):
    """mesh+band 合并频率交 imag_policy.classify_imag。返回 (imag_dict, vax_prim)。"""
    vax_prim = None
    if is2d:
        try:
            vax_prim = vacuum_axis_in_primitive(ph, cart_vac_axis)
        except Exception:                               # noqa: BLE001
            vax_prim = None
    freqs, qpts = imag_samples(ph)
    imag = imag_policy.classify_imag(freqs, qpts, is2d, vax_prim, cfg)
    return imag, vax_prim


def finalize(imag, za):
    """stability_verdict = combine(imag verdict, ZA verdict)；只有 fail 挡下一步。"""
    za_v = imag_policy.za_verdict(za)
    verdict = imag_policy.combine_verdict(imag["verdict"], za_v)
    return verdict, imag_policy.is_stable(verdict)


def za_check(cfg, ph, cart_vac_axis, is2d):
    """P2-2：2D 的 ZA 二次性闸门。返回写进 phonon_summary.json 的 dict（3D → None）。

    ★ 与 skill/fc-fit/fc_plot_phonon.py 的 _za_summary 同一套判据（同一 bug 的两份
    副本，见文件顶部 TODO）：
      * vac 轴映射：Cartesian 真空轴先映射成 原胞基矢 下标（vacuum_axis_in_primitive），
        否则生产原胞（真空在第 0 基矢）会量到"真空方向 + Γ-M"；
      * 面内第二方向：非正交原胞用 _k_sum_sign 取 e_i + s*e_j，否则 60° 生产原胞给出
        两个 Γ-M；
      * ZA 支识别：本征矢量面外占比 >= ZA_ZFRAC_MIN 的"最低频支"（旧代码取全局最低支）；
      * 质量字段：za_r2_* / za_log_resid_rms_* / za_margin_* / za_needs_review / 支号 /
        面外占比 / 混合说明；旧最低支结果并列保留但不作为判据；
      * 对称性审计：symprec 1e-5 与 1e-4 空间群不一致时，用 1e-4 的克隆做 ZA 拟合，
        默认构造的 p 记录在 za_exponent_default_qX 供对照。
    """
    mode = str(cfg.get("ZA_CHECK", "auto")).strip().lower()
    if mode in ("off", "false", "0", "no") or (not is2d and mode not in ("on", "true", "1", "yes")):
        return None
    p_lo, p_hi = ZA_P_RANGE
    qmax = float(cfg.get("ZA_QMAX", 0.05))
    res = {
        "mode": mode, "qmax": qmax, "p_range": [p_lo, p_hi], "is_2d": bool(is2d),
        "vacuum_axis_cartesian": (None if cart_vac_axis is None else int(cart_vac_axis)),
        "vacuum_axis_in_primitive": None,
        "dirs": [], "p": [], "min_freq": [], "p_lowest": [],
        "za_r2_q1": None, "za_r2_q2": None,
        "za_log_resid_rms_q1": None, "za_log_resid_rms_q2": None,
        "za_margin_q1": None, "za_margin_q2": None,
        "za_margin_min": ZA_MARGIN_MIN, "za_r2_min": ZA_R2_MIN,
        "za_needs_review": None,
        "za_branch_index_q1": None, "za_branch_index_q2": None,
        "za_zfrac_q1": None, "za_zfrac_q2": None,
        "za_zfrac_min_q1": None, "za_zfrac_min_q2": None,
        "za_branch_note_q1": "", "za_branch_note_q2": "",
        "za_lowest_exponent_q1": None, "za_lowest_exponent_q2": None,
        "za_lowest_r2_q1": None, "za_lowest_r2_q2": None,
        "za_lowest_log_resid_rms_q1": None, "za_lowest_log_resid_rms_q2": None,
        "za_symprec": ZA_SYMPREC_DEFAULT, "za_symprec_relax": ZA_SYMPREC_RELAX,
        "za_spacegroup_default": None, "za_spacegroup_relax": None,
        "za_symmetry_relaxed": None,
        "za_exponent_default_q1": None, "za_exponent_default_q2": None,
        "ok": False, "note": "",
    }
    # ★ 方向 bug 修复：Cartesian 真空轴 → 原胞基矢下标
    try:
        ax = vacuum_axis_in_primitive(ph, cart_vac_axis)
    except Exception as e:
        res["error"] = "真空轴 → 原胞基矢映射失败：%s" % e
        return res
    res["vacuum_axis_in_primitive"] = int(ax)
    try:
        qdirs = inplane_qdirs(ph, ax)
    except Exception as e:
        res["error"] = "面内方向选择失败：%s" % e
        return res
    res["dirs"] = [list(d) for d in qdirs]

    # --- 对称性审计（数值微畸变使默认 symprec 低估空间群 → ZA 假线性）---
    notes = []
    res["za_spacegroup_default"] = cell_symmetry(ph, ZA_SYMPREC_DEFAULT)
    res["za_spacegroup_relax"] = cell_symmetry(ph, ZA_SYMPREC_RELAX)
    relaxed = relaxed_symprec(ph, ZA_SYMPREC_RELAX)
    ph_za = ph
    if relaxed is not None:
        try:
            ph_za = clone_phonopy(ph, relaxed)
            res["za_symprec"] = relaxed
        except Exception as e:
            relaxed = None
            notes.append("无法按 symprec=%.0e 重建 phonopy（%s）" % (ZA_SYMPREC_RELAX, e))
    res["za_symmetry_relaxed"] = bool(ph_za is not ph)
    if ph_za is not ph:
        notes.append(
            "phonopy 默认 symprec=%.0e 看到 %s，而 symprec=%.0e 看到 %s：原胞只是数值上"
            "微畸变，ZA 拟合改用放宽容差（默认构造的 p 见 za_exponent_default_q1/q2）"
            % (ZA_SYMPREC_DEFAULT, res["za_spacegroup_default"],
               ZA_SYMPREC_RELAX, res["za_spacegroup_relax"]))

    had_nac = getattr(ph, "nac_params", None)
    tags = ("q1", "q2")
    ps, wmins, fits = [], [], []
    try:
        ph.nac_params = None
        try:
            ph_za.nac_params = None
        except Exception:
            pass
        for i, d in enumerate(qdirs):
            eig = None
            try:
                eig = za_power_law_eig(ph_za, qdir=d, qmax=qmax, n=12,
                                       cart_vac_axis=cart_vac_axis)
            except Exception as e:
                notes.append("qdir %s 本征矢量 ZA 拟合失败：%s" % (list(d), e))
            try:
                lp, lw, lf = za_power_law(ph_za, qdir=d, qmax=qmax, n=12)
            except Exception as e:
                lp, lw, lf = None, None, None
                notes.append("qdir %s 最低支拟合失败：%s" % (list(d), e))
            if eig is not None:
                p, wmin, fit, extra = eig
                res["za_branch_index_%s" % tags[i]] = extra["branch_index"]
                res["za_zfrac_%s" % tags[i]] = extra["zfrac_first"]
                res["za_zfrac_min_%s" % tags[i]] = extra["zfrac_min"]
                res["za_branch_note_%s" % tags[i]] = extra["note"]
                if extra["note"]:
                    notes.append("%s %s" % (tags[i], extra["note"]))
            else:
                # 无本征矢量（假 phonopy）→ 退回最低支
                p, wmin, fit = lp, lw, lf
            ps.append(p)
            wmins.append(wmin)
            fits.append(fit)
            res["p_lowest"].append(lp)
            res["za_lowest_exponent_%s" % tags[i]] = lp
            res["za_lowest_r2_%s" % tags[i]] = lf["r2"] if lf else None
            res["za_lowest_log_resid_rms_%s" % tags[i]] = (
                lf["log_resid_rms"] if lf else None)
            if ph_za is not ph:
                dflt = None
                try:
                    dflt = za_power_law_eig(ph, qdir=d, qmax=qmax, n=12,
                                            cart_vac_axis=cart_vac_axis)
                except Exception:
                    dflt = None
                res["za_exponent_default_%s" % tags[i]] = (
                    dflt[0] if dflt is not None else None)
            if p is None and wmin is not None:
                notes.append("qdir %s 无正频（omega_min=%.4f THz）" % (list(d), wmin))
    finally:
        try:
            ph.nac_params = had_nac
        except Exception:
            pass
    if ph_za is ph:
        res["za_exponent_default_q1"] = ps[0]
        res["za_exponent_default_q2"] = ps[1]
    res["p"] = list(ps)
    res["min_freq"] = list(wmins)
    res["za_r2_q1"], res["za_r2_q2"] = [(f["r2"] if f else None) for f in fits]
    res["za_log_resid_rms_q1"], res["za_log_resid_rms_q2"] = [
        (f["log_resid_rms"] if f else None) for f in fits]
    margins = [None if p is None else float(min(p - p_lo, p_hi - p)) for p in ps]
    res["za_margin_q1"], res["za_margin_q2"] = margins
    res["ok"] = bool(all(p is not None and p_lo < p < p_hi for p in ps))
    pstr = ["%.3f" % p if p is not None else "None" for p in ps]
    lowstr = ["%.3f" % p if p is not None else "None" for p in res["p_lowest"]]

    # 人工复核提示：p 贴近判据边界、或 log-log 拟合太脏时，粗判 ok 可能骗人
    review = False
    if is2d:
        for tag, p, fit, m in zip(tags, ps, fits, margins):
            if p is None or m is None:
                continue
            if m < ZA_MARGIN_MIN:
                review = True
                notes.append("%s p=%.3f 距 [%.1f, %.1f] 判据边界仅 %.3f（< %.2f）"
                             % (tag, p, p_lo, p_hi, m, ZA_MARGIN_MIN))
            r2 = fit.get("r2") if fit else None
            if r2 is not None and r2 < ZA_R2_MIN:
                review = True
                notes.append("%s log-log 拟合噪声过大（R2=%.4f < %.2f，log 残差 RMS=%.4f），"
                             "p 不可信" % (tag, r2, ZA_R2_MIN, fit.get("log_resid_rms")))
    res["za_needs_review"] = bool(review)
    if res["ok"]:
        notes.append("两个面内方向都有 %.1f < p < %.1f" % (p_lo, p_hi))
    else:
        notes.append("弯曲支在每个面内方向都不是二次色散")
    if review:
        notes.append("建议人工复核，不要仅凭 ok 下结论")
    if not res["ok"]:
        notes.append("弯曲支被线性化/有虚频时 2D 的 κ 会整体失真，S6 不启动。先查 pheasy 的 "
                     "RASR 是否真在 -c 步施加了 BHH"
                     "（pheasy_c.log 里的 Imposing rotational invariance and equilibrium conditions）、"
                     "结构是否还有残余面内应力。")
    res["note"] = ("ZA 指数 p(q1,q2)=%s（新判据=本征矢量面外占比≥%.1f 的最低频支，"
                   "vac 轴原胞下标=%s/笛卡尔=%s）；旧最低支 p=%s；%s"
                   % (pstr, ZA_ZFRAC_MIN, res["vacuum_axis_in_primitive"],
                      res["vacuum_axis_cartesian"], lowstr, "；".join(notes)))
    return res
