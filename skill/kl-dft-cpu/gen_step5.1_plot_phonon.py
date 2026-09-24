#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step5.1_plot_phonon.py —— S5.1 声子谱绘图（run: gen，不提交作业）。

读 step5_fc/phono3py 里拟合好的 fc2.hdf5 画声子谱（参考用户 plot_phonon_band.py）。
结构/超胞/原胞矩阵从 phono3py_disp.yaml 直接取（与拟合、虚频闸完全一致，避免 SPOSCAR 反推）。

【kl-dft-cpu-s51 路径修复】tf 的 remote_gen 是 `cd <材料>/kl-dft-cpu && python gen_step5.1_plot_phonon.py`，
CWD = 材料的技能目录（step5_fc/ 的兄弟层），不是步骤目录。原来写死 `../step5_fc/phono3py`
会跳到上一级，必然找不到 fc2.hdf5；而产物又写在 CWD，tf 去 step5_phonon_plot/ 找
done_marker 也必然落空（两个 bug 叠在一起）。
现在改成：
  · 输入目录 FCDIR 在几个候选位置里探测 fc2.hdf5（兼容 CWD=材料目录 / 步骤目录 / step5_fc）；
  · 自建 step5_phonon_plot/ 并 chdir 进去，所有产物（png/yaml/summary）都落在步骤目录里。
与 band-dft-cpu 技能的 gen_step3.1_plot_band.py 同一套约定：裸名访问兄弟目录、产物写自建目录。

NAC 规则（按用户要求）：
  · 目录里有 BORN（材料考虑了 NAC）→ 画两张，命名区分：band_nac.* 和 band_nonac.*；
  · 没有 BORN → 只画 band_nonac.*。
注意：2D 材料的 3D-NAC 在近 Γ 会有 LO-TO 假象，两张对照正好能看出 NAC 的影响
（虚频闸判据对 2D 用的是无 NAC，见 kl_fc_backends._stability_gate）。

产物（全部写入 step5_phonon_plot/）：
  band_nac.png/yaml、band_nonac.png/yaml（视有无 BORN）、phonon_plot_summary.json（done_marker）。
"""
import json
import os
import sys
from pathlib import Path

import numpy as np   # 画低频放大图/写数据表用（原来本文件没导入 numpy）
import imag_policy  # 虚频判据唯一真源（阈值常量与 THz→cm-1；gen_need 会带过来）

OUT_NAME = "step5_phonon_plot"       # 步骤目录名，必须与 skill.yaml 的 name 一致
FCSUB    = "phono3py"                # S5 拟合产物子目录
NPOINTS  = int(os.environ.get("KAPPA_PLOT_NPOINTS", "31"))
# 每段路径采样点。128 原子胞在集群登录节点上算一条完整 band structure 要几分钟
# （jzzn 实测 101 点/段 > 400 s，超过 tf 给 gen 的预算，步骤会被判失败），
# 31 点/段足够看清色散与软模，也留出余量；要更细就设 KAPPA_PLOT_NPOINTS=101。
# ★ 虚频判据不再在本步自行判定（2026-09-22）：唯一裁判是 S5，结论写在
#   ../step5_fc/phonon_summary.json 的 stability_verdict / imag_class / min_freq_THz。
#   本步只用它决定【图注文案】；summary 缺失时才用同一函数 imag_policy.classify_imag
#   配合 S5 的 step.conf 现算（并打 WARN）。阈值语义见 skill/_common/imag_policy.py。
S5_SUMMARY = "phonon_summary.json"
LOWF_MAX = 10.0                      # 低频放大图的上限（THz）：声学支/软模区
NL = chr(10)                         # 行尾常量（避免字符串里的 \n 被工具链吃掉）

FCDIR = None                         # 由 _locate() 定为绝对路径

# ==========================================================================
# gen 模式：只写 submit.sh，真正的算谱+出图交给计算节点上的 --work（2026-09-17
#   从 taskflow-v2.0 同步）。
#   原因见 templates/step5_phonon_plot/submit_plot.tpl 顶部：登录节点上一条完整
#   band structure 要几分钟到十几分钟，超过 autozt 给 gen 的预算会被中途杀掉。
# ==========================================================================
PLOT_SPEC = {
    "SBATCH_QOS": ("regular", "str"),
    "PLOT_CORES": (16, "int"),          # 该作业要几个核（画图/建动力学矩阵是 CPU）
    "CONDA_SH": ("", "str"),            # 由集群配置注入
    "CONDA_ENV": ("atomate2_p_a", "str"),
}
STEP = "step5_phonon_plot"


def _write_submit(outdir):
    """gen 模式：只写 submit.sh（作业脚本），谱的计算交给计算节点上的 --work。"""
    conf = {k: v[0] for k, v in PLOT_SPEC.items()}
    if Path("step.conf").is_file():
        try:
            import stepconf
            conf = stepconf.load(PLOT_SPEC, STEP)
        except BaseException as e:   # stepconf 用 sys.exit("[ERROR] ...") 报错
            print("[..] step.conf 读取失败（%s），用默认值" % e)
    else:
        print("[..] cwd 没有 step.conf（gen_need 未带上），用默认值：qos=%s cores=%s"
              % (conf["SBATCH_QOS"], conf["PLOT_CORES"]))
    here = Path(__file__).resolve().parent
    tpl = None
    for cand in (Path.cwd() / "submit_plot.tpl",
                 Path.cwd() / "templates" / "step5_phonon_plot" / "submit_plot.tpl",
                 here / "submit_plot.tpl",
                 here / "templates" / "step5_phonon_plot" / "submit_plot.tpl"):
        if cand.is_file():
            tpl = cand
            break
    if tpl is None:
        _emit({"status": "error",
               "reason": "缺 submit_plot.tpl（技能模板没推到远端？查过：cwd、"
                         "cwd/templates/step5_phonon_plot、脚本同级）"}, 40)
    print("[..] 模板：%s" % tpl)
    text = tpl.read_text(encoding="utf-8")
    subs = {"JOBNAME": "S51plot",
            "QOS": str(conf["SBATCH_QOS"] or "premium"),
            "NTASKS": str(int(conf["PLOT_CORES"] or 16)),
            "CONDA_SH": str(conf["CONDA_SH"]
                            or "/public/home/.../miniconda3/etc/"
                               "profile.d/conda.sh"),
            "CONDA_ENV": str(conf["CONDA_ENV"] or "atomate2_p_a"),
            # 作业的 cwd 是【步骤目录】，而本脚本在它的上一级（材料技能目录），
            # 所以必须写绝对路径——否则作业里 "can't open file" 直接失败。
            "PLOT_CMD": "python -u %s --work" % Path(__file__).resolve()}
    for k, v in subs.items():
        text = text.replace("{{%s}}" % k, str(v))
    (outdir / "submit.sh").write_text(text, encoding="utf-8", newline=NL)
    print("[OK] submit.sh 已写出（%s 核, qos=%s）→ 由 autozt 提交到计算节点"
          % (subs["NTASKS"], subs["QOS"]))


def gen_mode():
    """登录节点：只写 submit.sh，不真算谱。"""
    global FCDIR
    FCDIR, outdir = _locate()
    outdir.mkdir(parents=True, exist_ok=True)
    os.chdir(str(outdir))
    _write_submit(outdir)
    _emit({"status": "ok", "stage": "gen", "submit": "submit.sh",
           "fc_dir": str(FCDIR)}, 0)

def _emit(result, code):
    print(json.dumps(result, ensure_ascii=False), flush=True)
    sys.exit(code)


def _locate():
    """定位 S5 产物目录与本步输出目录，返回 (fcdir, outdir)，都是绝对路径。

    候选按 CWD 的三种可能排列（第一个能看到 fc2.hdf5 的胜出）：
      A. CWD = 材料技能目录 <mat>/kl-dft-cpu      —— tf remote_gen 的实际行为
      B. CWD = <mat>/kl-dft-cpu/step5_phonon_plot —— 手动进步骤目录里跑
      C. CWD = <mat>/kl-dft-cpu/step5_fc          —— 手动进 S5 目录里跑
    """
    cwd = Path.cwd().resolve()
    cands = [(cwd / "step5_fc" / FCSUB,        cwd / OUT_NAME),
             (cwd.parent / "step5_fc" / FCSUB, cwd if cwd.name == OUT_NAME
                                               else cwd.parent / OUT_NAME),
             (cwd / FCSUB,                     cwd.parent / OUT_NAME)]
    for fc, out in cands:
        if (fc / "fc2.hdf5").is_file():
            return fc.resolve(), out
    _emit({"status": "error",
           "reason": "缺 fc2.hdf5，先把 S5_fc 跑完（已查找：%s）"
                     % " | ".join(str(fc / "fc2.hdf5") for fc, _ in cands)}, 40)


def _build_phonopy():
    """从 phono3py_disp.yaml 建 Phonopy + 读 fc2.hdf5（每次新建，避免 NAC 状态串台）。"""
    import numpy as np
    import phono3py
    from phono3py.file_IO import read_fc2_from_hdf5
    from phonopy import Phonopy
    ph3 = phono3py.load(str(FCDIR / "phono3py_disp.yaml"),
                        produce_fc=False, is_nac=False, log_level=0)
    ph = Phonopy(ph3.unitcell, supercell_matrix=ph3.supercell_matrix,
                 primitive_matrix=ph3.primitive_matrix)
    ph.force_constants = np.asarray(read_fc2_from_hdf5(filename=str(FCDIR / "fc2.hdf5")))
    return ph


def _set_nac(ph, on):
    if on:
        from phonopy.file_IO import parse_BORN
        nac = parse_BORN(ph.primitive, filename=str(FCDIR / "BORN"))
        if isinstance(nac, dict) and not nac.get("factor"):
            nac["factor"] = 14.399652     # VASP NAC 单位换算因子（phonopy-vasp-born 常缺）
        ph.nac_params = nac
    else:
        ph.nac_params = None              # 显式建"无 NAC"动力学矩阵


def _dim_axis():
    """读 kl_params 的 DIM；读不到就用 POSCAR 现判。返回 (dim, vac_axis)。"""
    import kl_common as kc
    dim, ax = "", 2
    for cand in (FCDIR.parent / kc.KL_PARAMS, FCDIR / kc.KL_PARAMS,
                 Path.cwd() / kc.KL_PARAMS):
        if cand.is_file():
            dim = (kc.read_kl_params(cand).get("DIM") or "").lower()
            break
    pos = next((p for p in (FCDIR / "POSCAR", FCDIR.parent / "POSCAR",
                            Path.cwd() / "POSCAR") if p.is_file()), None)
    if pos is not None:
        try:
            d2, ax2 = kc.resolve_dim(pos, dim or "auto")
            dim = dim or d2
            ax = int(ax2 if ax2 is not None else 2)
        except Exception:
            pass
    return dim, ax


def _read_s5_imag():
    """读 S5 唯一裁判的结论（<step5_fc>/phonon_summary.json）；缺失/未写入返回 None。"""
    f = FCDIR.parent / S5_SUMMARY
    if not f.is_file():
        return None
    try:
        d = json.loads(f.read_text(encoding="utf-8"))
    except Exception:                              # noqa: BLE001
        return None
    return d if d.get("stability_verdict") else None


def _imag_annotation(verdict, imag_class, min_thz):
    """图注文案（规则见 skill/_common/imag_policy.py）。

    near_gamma_acoustic → "近Γ声学支软化 (−x.xxx THz)"；fail → "含虚频"；noise/none 不标。
    """
    if imag_class == "near_gamma_acoustic":
        try:
            return "近Γ声学支软化 (%.4f THz)" % float(min_thz)
        except (TypeError, ValueError):
            return "近Γ声学支软化"
    if str(verdict) == "fail":
        return "含虚频"
    return ""


def _s5_conf():
    """从 S5 的 step.conf 取三个 imag 阈值；读不到返回 None（→ imag_policy 默认值）。"""
    conf = {}
    for cand in (FCDIR.parent / "step.conf", Path.cwd() / "step.conf"):
        if not cand.is_file():
            continue
        try:
            for ln in cand.read_text(encoding="utf-8", errors="replace").splitlines():
                ln = ln.split("#", 1)[0].strip()
                if "=" not in ln:
                    continue
                k, v = [x.strip() for x in ln.split("=", 1)]
                if k in ("IMAG_THR", "IMAG_THR_STRICT", "IMAG_QGAMMA"):
                    try:
                        conf[k] = float(v.strip().strip('"').strip("'"))
                    except ValueError:
                        pass
        except OSError:
            pass
    return conf or None


def _imag_fallback(ph):
    """phonon_summary.json 缺失时现算：与 S5 同一函数 + 同一套 conf（打 WARN）。"""
    import imag_policy
    print("[WARN] 找不到 S5 的 %s —— 用 imag_policy 现算（与 S5 同一判据）"
          % (FCDIR.parent / S5_SUMMARY))
    dim, vax = _dim_axis()
    is2d = str(dim).lower().startswith("2")
    vaxp, freqs, qpts = None, [], []
    try:
        if is2d:
            import kl_common as kc
            vaxp = kc.vacuum_axis_in_primitive(ph.primitive.cell, vax)
            mesh = [60, 60, 60]
            mesh[int(vaxp)] = 1
            ph.run_mesh(mesh=mesh, with_eigenvectors=False, is_mesh_symmetry=True,
                        is_gamma_center=True)
        else:
            ph.run_mesh(mesh=60.0, with_eigenvectors=False, is_mesh_symmetry=True)
        md = ph.get_mesh_dict()
        freqs += [[float(x) for x in row] for row in md["frequencies"]]
        qpts += [[float(x) for x in q] for q in md["qpoints"]]
    except Exception as e:                         # noqa: BLE001
        print("[WARN] 现算 mesh 失败：%s" % e)
    try:
        bsd = ph.get_band_structure_dict()
        freqs += [[float(x) for x in f]
                  for seg in (bsd.get("frequencies") or []) for f in seg]
        qpts += [[float(x) for x in q]
                 for seg in (bsd.get("qpoints") or []) for q in seg]
    except Exception as e:                         # noqa: BLE001
        print("[WARN] 现算 band 失败：%s" % e)
    r = imag_policy.classify_imag(freqs, qpts, is2d, vaxp, _s5_conf())
    print("[WARN] 现算结论：verdict=%s imag_class=%s min=%.4f THz (%.2f cm-1)"
          % (r["verdict"], r["imag_class"], r["min_freq_THz"] or 0.0,
             r["min_freq_cm1"] or 0.0))
    return {"stability_verdict": r["verdict"], "imag_class": r["imag_class"],
            "min_freq_THz": r["min_freq_THz"], "verdict": r["verdict"]}


_S5_IMAG_CACHE = None


def _s5_imag(ph):
    """S5 的虚频结论（读 summary；缺失则现算一次并缓存）。"""
    global _S5_IMAG_CACHE
    if _S5_IMAG_CACHE is None:
        _S5_IMAG_CACHE = _read_s5_imag() or _imag_fallback(ph)
        print("[..] 虚频结论（S5 唯一裁判）：stability_verdict=%s imag_class=%s min=%.4f THz"
              % (_S5_IMAG_CACHE.get("stability_verdict"), _S5_IMAG_CACHE.get("imag_class"),
                 _S5_IMAG_CACHE.get("min_freq_THz") or 0.0))
    return _S5_IMAG_CACHE


def _band_path(ph):
    """高对称路径：2D 走 kz=0 的二维路径表；3D 优先 seekpath，失败退回六方通用路径。

    P2-1：seekpath/VASPKIT-301 是 3D 工具，用在 2D 上会给出 Γ-A 这类 kz 线段
    （真空方向的平凡色散，白占图幅还误导判读）。
    """
    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
    dim, vax = _dim_axis()
    if dim == "2d":
        import kl_common as kc
        try:
            # ★ vax 是 POSCAR 的 Cartesian 真空轴；primitive_matrix 常把真空换到别的原胞基矢
            #   （生产 MoS₂：a1=(0,0,-25)，真空=基矢 0），必须映射成原胞基矢下标 ——
            #   否则 band_path_2d 会拿"25 Å 真空基矢 + 一个面内基矢"判成 rectangular、
            #   把六方路径画成 Γ-X-S-Y-Γ（2026-09-21 真实数据实测）。
            vax_prim = kc.vacuum_axis_in_primitive(ph.primitive.cell, vax)
            paths, labels, lat = kc.band_path_2d(ph.primitive.cell, NPOINTS, vax_prim)
            bands, connections = get_band_qpoints_and_path_connections(paths, npoints=NPOINTS)
            print("[..] 2D 高对称路径：%s 格子，kz=0：%s"
                  % (lat, "-".join(l.replace("$\\Gamma$", "Γ") for l in labels)))
            return bands, connections, labels, "2d-%s" % lat
        except Exception as e:
            print("[WARN] 2D 路径生成失败，退回六方通用路径：%s" % e)
    try:
        from phonopy.phonon.band_structure import get_band_qpoints_by_seekpath
        bands, labels, connections = get_band_qpoints_by_seekpath(ph.primitive, NPOINTS)
        return bands, connections, labels, "seekpath"
    except Exception as e:
        print("[..] seekpath 不可用（%s），用通用路径 Γ-M-K-Γ" % type(e).__name__)
        paths = [[[0, 0, 0], [0.5, 0, 0]], [[0.5, 0, 0], [1./3, 1./3, 0]],
                 [[1./3, 1./3, 0], [0, 0, 0]]]
        bands, connections = get_band_qpoints_and_path_connections(paths, npoints=NPOINTS)
        return bands, connections, ["$\\Gamma$", "M", "K", "$\\Gamma$"], "fallback"


def _cached_band(tag):
    """已算好的 band_<tag>.yaml（上一次运行被超时打断留下的）→ dict，否则 None。

    登录节点上算不动整条路径时，谱本身在被打断前就写进 yaml 了；下一步重跑
    直接读它，几秒钟就能出图和数据表，不必再算一遍（tf 的 gen 预算有限，
    "算谱" 与 "画图存数据" 必须能分开重试）。
    """
    p = Path("band_%s.yaml" % tag)
    if not p.is_file():
        return None
    try:
        import yaml
        d = yaml.safe_load(p.read_text(encoding="utf-8"))
        # phonopy 的 band yaml 把【所有】q 点摊平成一个 list（每项：q-position /
        # distance / band[分支]），而 get_band_structure_dict() 是按【路径段】分组的
        # —— 这里用顶层的 segment_nqpoint 还原成段布局，两种来源才对得上。
        segs = d["phonon"]
        seg_n = [int(x) for x in (d.get("segment_nqpoint") or [len(segs)])]
        dists, freqs, i = [], [], 0
        for nq in seg_n:
            chunk = segs[i:i + nq]
            i += nq
            if not chunk:
                continue
            dists.append(np.array([float(c.get("distance", 0.0)) for c in chunk],
                                  dtype=float))
            freqs.append(np.array([[float(b["frequency"]) for b in c["band"]]
                                   for c in chunk], dtype=float))
        if not freqs or not any(float(dd.max()) > 0.0 for dd in dists):
            raise ValueError("yaml 里没有可用的 distance/band 数据")
        return {"distances": dists, "frequencies": freqs,
                "labels": d.get("labels") or []}
    except Exception as e:
        print("[..] band_%s.yaml 缓存不可用（%s），重新计算" % (tag, e))
        return None


def _plot(tag):
    """tag='nac' | 'nonac'：画声子谱（全频段 + 0..LOWF_MAX 放大）+ 数据表。"""
    import matplotlib
    matplotlib.use("Agg")
    # 先按【当前代码】定路径，并写一个路径指纹 sidecar：vax→原胞基矢的修复会改变路径
    # （2026-09-21 实测：MoS₂ 未修前画成 Γ-X-S-Y-Γ，修复后才是 Γ-M-K-Γ），旧 band_*.yaml
    # 必须失效重算，否则断点续画会把错误路径的谱一直缓存下去。
    ph = _build_phonopy()
    _set_nac(ph, tag == "nac")
    bands, conn, labels, src = _band_path(ph)
    fp = json.dumps({"source": src, "labels": labels}, ensure_ascii=False, sort_keys=True)
    side = Path("band_%s.path" % tag)
    bs = None
    if side.is_file():
        try:
            if side.read_text(encoding="utf-8").strip() == fp:
                bs = _cached_band(tag)
        except OSError:
            bs = None
    if bs is not None:
        print("[..] 复用已算好的 band_%s.yaml（路径指纹一致：%s）" % (tag, src))
    else:
        if side.is_file():
            old_fp = "?"
            try:
                old_fp = json.loads(side.read_text(encoding="utf-8")).get("source", "?")
            except Exception:                              # noqa: BLE001
                pass
            print("[..] 路径指纹变化（旧 %s → 新 %s）：丢弃缓存重算" % (old_fp, src))
        ph.run_band_structure(bands, path_connections=conn, labels=labels,
                              with_eigenvectors=False)
        # 先落盘谱数据：即使后面绘图阶段被打断，数据也已保住
        ph.write_yaml_band_structure(filename="band_%s.yaml" % tag)
        try:
            side.write_text(fp, encoding="utf-8")
        except OSError:
            pass
        bs = ph.get_band_structure_dict()
        plt = ph.plot_band_structure()
        plt.savefig("band_%s.png" % tag, dpi=200, bbox_inches="tight")
        try:
            plt.close("all")
        except Exception:
            pass
    bs = {"distances": [np.asarray(d, float) for d in bs["distances"]],
          "frequencies": [np.asarray(fr, float) for fr in bs["frequencies"]],
          "labels": bs.get("labels") or []}
    fmin = float(min(fr.min() for fr in bs["frequencies"]))
    # 虚频结论来自 S5（唯一裁判）；本步只据此写图注，不再自己判阈值
    _im = _s5_imag(ph)
    _annot = _imag_annotation(_im.get("stability_verdict"), _im.get("imag_class"),
                              _im.get("min_freq_THz"))
    # 数据文件：距离-频率表，便于复算/画自己的图（不只留 phonopy 的 yaml）
    ndata = _dump_band_data(bs, tag)
    # 低频放大图（0..LOWF_MAX THz）：声学支与软模/虚频一眼可见
    low = _draw_band(bs, tag, fmin, LOWF_MAX, "_lowfreq", _annot)
    # 全频段图：走缓存时（上一次被超时打断）phonopy 那张没画出来，这里补上
    if not Path("band_%s.png" % tag).is_file():
        _draw_band(bs, tag, fmin, None, "", _annot)
    print("[OK] band_%s.png / band_%s.yaml / band_%s.dat / %s  路径=%s  "
          "最低频率=%.3f THz (%.2f cm-1)%s"
          % (tag, tag, tag, low, src, fmin, fmin * imag_policy.THZ_TO_CM1,
             ("  ⚠️" + _annot) if _annot else ""))
    return {"tag": tag, "nac": tag == "nac", "png": "band_%s.png" % tag,
            "yaml": "band_%s.yaml" % tag, "dat": "band_%s.dat" % tag,
            "lowfreq_png": low, "n_band_rows": ndata,
            "min_freq_THz": round(fmin, 4),
            "stability_verdict": _im.get("stability_verdict"),
            "imag_class": _im.get("imag_class"),
            "s5_min_freq_THz": _im.get("min_freq_THz"),
            "imaginary": str(_im.get("stability_verdict")) == "fail",
            "path_source": src}


def _dump_band_data(bs, tag):
    """band_<tag>.dat：段号/段内序号/路径距离/各支频率（THz）的纯文本表。"""
    n = 0
    with open("band_%s.dat" % tag, "w", encoding="utf-8") as fh:
        fh.write("# segment q_index distance_1overA  frequencies_THz(...)" + NL)
        for si, (d, fr) in enumerate(zip(bs["distances"], bs["frequencies"])):
            d = np.asarray(d, float)
            fr = np.asarray(fr, float)
            for qi in range(fr.shape[0]):
                fh.write("%d %d %.6f %s" % (si, qi, d[qi],
                         " ".join("%.6f" % v for v in fr[qi])) + NL)
                n += 1
    return n


def _draw_band(bs, tag, fmin, ymax, suffix, annot=""):
    """画一张声子谱：ymax=None → 全频段；ymax=LOWF_MAX → 0..10 THz 放大图。

    整条谱的上限 40+ THz 会把声学支压成一条线，软模/虚频看不见，所以低频那张
    单独出：声学支色散 + fmin<0 的虚频位置。返回文件名。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for d, fr in zip(bs["distances"], bs["frequencies"]):
        d = np.asarray(d, float)
        fr = np.asarray(fr, float)
        for b in range(fr.shape[1]):
            ax.plot(d, fr[:, b], "-", lw=0.9, color="#1f4e79")
    ax.axhline(0.0, color="0.55", lw=0.8, ls="--")
    if ymax is None:                       # 全频段
        fmax = float(max(fr.max() for fr in bs["frequencies"]))
        ax.set_ylim(min(-1.0, fmin - 3.0), fmax + 3.0)
        ax.set_title("Phonon dispersion (%s%s)"
                     % (tag, (", " + annot) if annot else ""), fontsize=10)
    else:                                  # 0..ymax THz 放大
        ax.set_ylim(min(-1.0, fmin - 0.5), float(ymax))
        ax.set_title("Phonon dispersion (low frequency <= %g THz, %s%s)"
                     % (ymax, tag, (", " + annot) if annot else ""), fontsize=10)
    ax.set_xlim(0.0, float(np.asarray(bs["distances"][-1], float)[-1]))
    labels = bs.get("labels") or []
    if labels:
        try:
            ax.set_xticks([float(np.asarray(bs["distances"][i], float)[0])
                           for i, lb in enumerate(labels) if lb])
            ax.set_xticklabels([str(lb).replace("$", "") for lb in labels if lb])
        except Exception:
            pass
    ax.set_xlabel("Wave vector")
    ax.set_ylabel("Frequency (THz)")
    name = "band_%s%s.png" % (tag, suffix)
    fig.savefig(name, dpi=200, bbox_inches="tight")
    try:
        plt.close(fig)
    except Exception:
        pass
    return name


def main():
    global FCDIR
    FCDIR, outdir = _locate()
    outdir.mkdir(parents=True, exist_ok=True)
    os.chdir(str(outdir))                 # 之后所有相对写入都落在步骤目录里
    out = Path.cwd()
    print("[..] fc 源：%s\n[..] 产物目录：%s" % (FCDIR, out))

    has_born = (FCDIR / "BORN").is_file()
    tags = ["nac", "nonac"] if has_born else ["nonac"]
    print("[..] NAC=%s → 输出 %s" % ("有(BORN)" if has_born else "无",
                                     " + ".join("band_%s" % t for t in tags)))

    plots = []
    for t in tags:
        try:
            plots.append(_plot(t))
        except Exception as e:
            print("[WARN] band_%s 绘制失败：%s" % (t, e))
    if not plots:
        _emit({"status": "error", "reason": "声子谱全部失败（见上 WARN）"}, 40)

    summary = {"status": "ok", "nac_considered": has_born,
               "fc_dir": str(FCDIR), "out_dir": str(out),
               "note": ("有 NAC：band_nac + band_nonac 两张对照"
                        if has_born else "无 NAC：仅 band_nonac"),
               "lowfreq_max_THz": LOWF_MAX,
               "data_files": [p["dat"] for p in plots],
               "lowfreq_figures": [p["lowfreq_png"] for p in plots],
               "plots": plots}
    (out / "phonon_plot_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    print("[DONE] S5.1 声子谱：%s" % ", ".join(p["png"] for p in plots))
    _emit(summary, 0)


if __name__ == "__main__":
    # 两种运行方式（2026-09-17 从 taskflow-v2.0 同步）：
    #   登录节点 gen：只写 submit.sh（gen_mode）
    #   计算节点作业：python gen_step5.1_plot_phonon.py --work → 真算谱+出图（main）
    if "--work" in sys.argv:
        main()
    else:
        gen_mode()
