# -*- coding: utf-8 -*-
"""kl_common.py —— 晶格热导率技能各步骤公共工具（step.conf 版）。

依赖：标准库 + dim_common（同目录）。重物理（ALM 定位移数、NAC/Born 提取、
拟合）由各 gen 脚本按需 import lattice_kappa（同目录，随 gen_need 推送）复用参考
引擎里已验证的函数，本文件只放轻量编排/字符串/子进程工具，保证登录节点系统
python 也能跑通编排部分。

放置：随每个用它的步骤 gen_need 推送到材料目录（skill 库里就一份在技能根，
per_step 布局下 find_asset 会从技能根回落取到）。
"""
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
# dim_common 在集群上随 gen_need 与 kl_common.py 同目录；仓库里它在 skill/_common/opt。
# 独立 CLI（python kl_common.py --scan-nelm …）可能在别处执行，逐级回退找它。
for _cand in (_HERE.parent / "_common" / "opt", _HERE.parent / "_common"):
    if (_cand / "dim_common.py").is_file():
        sys.path.insert(0, str(_cand))
        break
from dim_common import (detect_dimension, force_kz1, read_poscar_cell_frac,  # noqa: E402
                        validate_poscar, resolve_tpl, VACUUM_MIN,
                        _norm, _cross, _det3)

METHOD_FILE = "workflow_method.txt"     # step1 写的 FUNC/GGA/DIM/MAG（继承泛函/维度）
KL_PARAMS   = "kl_params.txt"           # 本技能跨步共享：DIM/SUPERCELL/MESH/METHOD/NAC

# phono3py / alm 所在 conda 环境（按集群改；与 submit_*.tpl 保持一致）。
# 非交互 shell 里 conda activate 前必须 source conda.sh。
PHONO3PY_ENV_SRC = ("source /opt/miniconda3/etc/profile.d/conda.sh "
                    "&& conda activate atomate2_p_a")


# ==========================================================================
# 维度 / 泛函继承（读 step1 的 workflow_method.txt）
# ==========================================================================
def read_method(method_file):
    """读 step1 的 workflow_method.txt → dict（FUNC/GGA/IVDW/DIM/MAG 等，键大写）。"""
    d = {}
    p = Path(method_file)
    if not p.is_file():
        return d
    for ln in p.read_text(errors="ignore").splitlines():
        s = ln.strip()
        if s and not s.startswith("#") and "=" in s:
            k, v = s.split("=", 1)
            d[k.strip().upper()] = v.strip()
    return d


def read_method_dim(method_file):
    v = read_method(method_file).get("DIM", "").lower()
    return v if v in ("2d", "3d") else None


def resolve_dim(poscar, dimension="auto", vacuum_min=VACUUM_MIN):
    """返回 (dim, vac_axis)。dim ∈ {'2d','3d'}；vac_axis 仅 2D 有意义。"""
    mode = str(dimension).lower()
    if mode in ("2d", "3d"):
        return mode, (2 if mode == "2d" else None)
    dim, axis, _ = detect_dimension(poscar, vacuum_min)
    if dim == "2d" and axis != 2:
        sys.exit("[ERROR] 检测到 2D 但真空不在 c 轴（在 %d 轴）。请把结构旋转成"
                 "真空沿第 3 个晶格矢量再重跑（step1 会自动轮换，从 step1 重跑即可）。" % axis)
    return dim, (axis if dim == "2d" else None)


GGA_MAP = {"pbe": "PE", "pbesol": "PS", "pbe-d3": "PE"}
VDW_MAP = {"pbe": None, "pbesol": None, "pbe-d3": "12"}


# ==========================================================================
# 超胞倍数 / phono3py --dim / --mesh 字符串（2D 真空方向恒 1）
# ==========================================================================
def _supercell_geometry(poscar, reps):
    """超胞几何（纯标准库）：返回 (lengths, perp, insphere)。
    lengths = 三条边模长；perp = 每方向垂直胞高 V/|a_j×a_k|（周期性最短镜像距离，
    恒 ≤ 对应边长）；insphere = min(perp) = 内切球直径。"""
    lat, _ = read_poscar_cell_frac(poscar)
    sc = [[reps[i] * lat[i][k] for k in range(3)] for i in range(3)]
    lengths = [_norm(sc[i]) for i in range(3)]
    vol = abs(_det3(sc))
    perp = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        area = _norm(_cross(sc[j], sc[k]))
        perp.append(vol / area if area > 1e-12 else 0.0)
    return lengths, perp, min(perp)


def warn_cutoff_vs_supercell(poscar, reps, cutoff, label="cutoff", margin=0.1,
                             dim="3d", vac_axis=2):
    """截断半径必须 ≤ 超胞安全截断（0.5×内切球直径 − margin），否则周期性镜像让
    力常数对重复计数（参考 lattice_kappa._max_safe_cutoff）。只 WARN 不拦截。
    返回安全截断 (Å)；cutoff 为 None 时返回 None。

    P1-5：2D 只按【面内】方向判。真空方向的"最短镜像距离"是层间距加真空（很大），
    把它算进 min(perp) 没有意义；而跨真空的相互作用另有专门判据（见 check_2d_vacuum
    的"截断 ≥ 真空隙"硬闸）。"""
    if cutoff is None:
        return None
    _, perp, insphere = _supercell_geometry(poscar, reps)
    if str(dim).lower() == "2d":
        ax = int(vac_axis if vac_axis is not None else 2)
        inplane = [perp[i] for i in range(3) if i != ax]
        if inplane:
            insphere = min(inplane)
    max_safe = 0.5 * insphere - margin
    if float(cutoff) > max_safe:
        print("[WARN] %s=%.2f Å > 超胞安全截断 %.2f Å（内切球直径 %.2f Å）——"
              "周期镜像会污染力常数！建议 %s ≤ %.2f Å，或扩超胞（MIN_SC_LEN/SUPERCELL）。"
              % (label, float(cutoff), max_safe, insphere, label, max_safe))
    return max_safe


def check_2d_vacuum(poscar, vac_axis, cut3=None, min_vac=15.0):
    """2D 真空层检查（P1-5）。返回 slab 几何 dict，出错直接 sys.exit。

    为什么是硬闸：pheasy 的 Born-Huang 约束用 Wigner-Seitz 最近镜像的 r_ij
    （symmetry_constraints.py 的 ws_offsets），2D 超胞在真空方向倍数为 1。
    真空隙一旦 ≤ 层厚，层内某些原子对的"最近镜像"就会跨过真空，r_z 是错的，
    力常数整体不可信；而三阶截断若 ≥ 真空隙，还会直接引入跨真空的三体簇。

    真空 ≥ min_vac（默认 15 Å）只 WARN：够不够跟体系有关，文献里石墨烯/WS2 用
    20 Å、MoS2 用 25 Å（那些工作还配了 2D 库仑截断），VASP 侧只能靠加厚真空。
    """
    from thickness_2d import slab_geometry_from_poscar
    ax = vac_axis if vac_axis is not None else 2
    g = slab_geometry_from_poscar(poscar, vac_axis=ax, mode="vdw")
    gap, span = float(g["vacuum_gap_A"]), float(g["atomic_span_A"])
    if gap <= span:
        sys.exit("[ERROR] 2D 真空隙 %.2f Å ≤ 层厚（原子起伏）%.2f Å：Wigner-Seitz 最近镜像"
                 "会跨真空，Born-Huang 约束里的 r_ij 是错的，力常数不可信。\n"
                 "        请把真空加到 ≥ 15 Å（最好 20 Å）后从 S1 重跑结构弛豫。"
                 "\n        当前值由 _common/thickness_2d.slab_geometry 算出（h⊥=%.2f Å）。"
                 % (gap, span, float(g["h_perp_A"])))
    if gap < float(min_vac):
        print("[WARN] 2D 真空隙 %.2f Å < %.1f Å —— 层间虚假相互作用不可忽略；"
              "文献里石墨烯/WS2 用 20 Å、MoS2 用 25 Å。当前 h⊥=%.2f Å，层厚 %.2f Å。"
              % (gap, float(min_vac), float(g["h_perp_A"]), span))
    if cut3 is not None and float(cut3) >= gap:
        sys.exit("[ERROR] 三阶截断 %.2f Å ≥ 真空隙 %.2f Å：会在真空方向引入跨层的三体簇，"
                 "力常数被污染。\n        请把截断收到 < %.2f Å（step.conf 的 ALM_CUT3 / "
                 "FC3_CUTOFF_PAIR），或加厚真空后重跑。"
                 % (float(cut3), gap, gap))
    print("[..] 2D 层几何：h⊥=%.2f Å  层起伏=%.2f Å  真空隙=%.2f Å  层厚 d=%.2f Å（%s）"
          % (float(g["h_perp_A"]), span, gap, float(g["thickness_d_A"]),
             g.get("thickness_convention", "?")))
    return g


_STRESS_RE = re.compile(
    r"^\s*in\s+kB\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", re.M)


def last_stress_kbar(outcar):
    """读 OUTCAR 里【最后一次】"in kB" 应力行，返回 (xx,yy,zz,xy,yz,zx) kB 或 None。

    VASP 的列序是 XX YY ZZ XY YZ ZX。ISIF=2 不写应力行（返回 None）；分段弛豫的
    ISIF=3 变胞段会写，取最后一次就是弛豫完的末态应力。
    """
    p = Path(outcar)
    if not p.is_file():
        return None
    txt = p.read_text(errors="ignore")
    ms = _STRESS_RE.findall(txt)
    if not ms:
        return None
    return tuple(float(x) for x in ms[-1])


def check_2d_stress(outcar, h_perp_A, d_A, thr_kbar=0.5, quiet=False):
    """P1-4：2D 残余应力门禁。返回诊断 dict（也写进 summary）；超阈值 sys.exit。

    Huang 条件本质是零应力条件：只要面内有外应力，ZA 色散里就会出现线性项、压制
    二次项（Lin 等对 MoS2 用的判据是 1e-3 kbar 量级）。VASP 报的应力是对整个胞
    （含真空）平均的，层内真实应力 ≈ σ_VASP × h⊥/d —— 真空越厚这个换算越重要，
    直接用 σ_VASP 会低估真实层内应力一个量级。

    σ_zz 只记录不判（c 轴固定）。thr<=0 表示显式关掉门禁。
    """
    import math
    sig = last_stress_kbar(outcar)
    fac = (float(h_perp_A) / float(d_A)) if float(d_A) else 1.0
    if sig is None:
        if not quiet:
            print("[WARN] 读不到 OUTCAR 的应力（%s 里没有 'in kB' 行——ISIF=2 不算应力？）。"
                  "2D 的残余面内应力无法核验；Huang 条件不满足时 ZA 会线性化、κ 失真。"
                  % outcar)
        return {"available": False, "source": str(outcar)}
    xx, yy, zz, xy, yz, zx = sig
    # Voigt 剪切分量是工程应变口径的一半，2D 面内剪切只看 σ_xy
    vals = {"xx": xx, "yy": yy, "xy": xy}
    layer = {k: v * fac for k, v in vals.items()}
    worst = max(layer, key=lambda k: abs(layer[k]))
    res = {"available": True, "source": str(outcar), "cell_kbar": vals,
           "layer_kbar": layer, "sigma_zz_cell_kbar": zz, "sigma_zz_layer_kbar": zz * fac,
           "h_perp_over_d": fac, "threshold_kbar": float(thr_kbar),
           "worst": worst, "worst_layer_kbar": layer[worst],
           "ok": float(thr_kbar) <= 0 or all(abs(v) < float(thr_kbar) for v in layer.values())}
    if not quiet:
        print("[..] 2D 残余应力（σ_VASP × h⊥/d，h⊥/d=%.4f）：xx=%.3f yy=%.3f xy=%.3f kbar"
              "（层内）；σ_zz(原胞)=%.3f kbar 仅记录"
              % (fac, layer["xx"], layer["yy"], layer["xy"], zz))
    if float(thr_kbar) > 0 and not res["ok"]:
        sys.exit("[ERROR] 2D 残余面内应力过大：|σ_%s|=%.3f kbar（层内口径）≥ %.2f kbar。\n"
                 "        Huang 条件是零应力条件 —— 有外应力时 ZA 弯曲支会出现线性项、"
                 "压制二次项，声子谱与 κ 都不可信。\n"
                 "        先回 S1 把结构/晶格收干净（VASP 对 2D 建议用 "
                 "LATTICE_CONSTRAINTS 或 OPTCELL 只放面内自由度），或显式放宽 "
                 "STRESS_2D_THR。"
                 % (worst, layer[worst], float(thr_kbar)))
    return res


def enforce_incar_tags(path, tags, label=""):
    """把 tags（{标签: 值}）强制写进 INCAR：删掉同名旧行再统一追加。

    为什么需要：项目/集群级的 incar_*.tpl 优先级高于技能模板（find_asset），老项目里
    那些副本没有新加的占位符（{{FORCE_PREC}} / {{DIPOLE_LINE}}…），新修的东西会被
    **静默吃掉**。渲染后强制落一遍，缺了就打 WARN，保证物理设置真的进了 INCAR。
    返回被修正的标签列表（空 = 模板已经写对）。
    """
    p = Path(path)
    if not p.is_file():
        return []
    text = p.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()
    fixed = []
    for k, v in tags.items():
        pat = re.compile(r"^\s*%s\s*=" % re.escape(k), re.I)
        hit = [i for i, ln in enumerate(lines)
               if pat.match(ln.split("#")[0].split("!")[0])]
        if not hit:
            lines.append("%s = %s" % (k, v))
            fixed.append(k)
    if fixed:
        p.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        print("[WARN] %sINCAR 缺 %s（多半是项目/集群级模板没有对应占位符，被 find_asset "
              "遮蔽了）→ 已强制补写。" % (label, ", ".join(fixed)))
    return fixed


def _vasprun_basis_positions(path):
    """读 vasprun.xml 第一个 <structure> 的晶格与（分数）坐标。返回 (basis, positions) 或 None。"""
    import re as _re
    try:
        txt = Path(path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    def _arr(name):
        m = _re.search(r'<varray name="%s"\s*>(.*?)</varray>' % name, txt, _re.S)
        if not m:
            return None
        return [[float(x) for x in row.split()]
                for row in _re.findall(r"<v>\s*([^<]+?)\s*</v>", m.group(1))]
    b, p = _arr("basis"), _arr("positions")
    return (b, p) if b and p else None


def check_frames_match_displacements(step4_dir, n_sample=3, tol=1e-3):
    """抽几帧：vasprun.xml 的坐标必须等于 phono3py_disp.yaml 的「超胞 + 位移」。

    为什么需要（2026-09-17，user 要求）：S4 是 fanout，retry 只补缺失帧；若 S1 换了
    结构而旧 disp-* 残留，就会出现"旧结构的力 + 新位移"的静默错配 —— 检查晶格能拦住
    大部分，但同胞不同位移（换过 rattle 种子/帧数）仍可能漏过，所以这里**逐原子比对坐标**。
    返回 (ok, note)。
    """
    import yaml as _yaml
    import numpy as np
    d4 = Path(step4_dir)
    y = d4 / "phono3py_disp.yaml"
    if not y.is_file():
        return True, "无 phono3py_disp.yaml，跳过抽帧校验"
    try:
        dd = _yaml.safe_load(y.read_text(encoding="utf-8"))
        sc = dd["supercell"]
        lat = np.array(sc["lattice"], dtype=float)
        # phonopy/phono3py 的 points 是 [{symbol, coordinates, mass, reduced_to}, ...]，
        # 分数坐标在 coordinates 里（2026-09-17 实测；直接 np.array(points) 会报
        # "float() argument must be ... not 'dict'"）。
        pts = np.array([p["coordinates"] if isinstance(p, dict) else p for p in sc["points"]],
                       dtype=float)
        # 随机位移（type-2 dataset）的位移记录在 dataset.displacements，
        # 顶层 displacements 可能是 None。
        _ds = dd.get("dataset") or {}
        disps = dd.get("displacements") or _ds.get("displacements") or []
    except Exception as e:                              # noqa: BLE001
        return True, "phono3py_disp.yaml 解析失败（%s），跳过抽帧校验" % e
    sc_cart = pts @ lat
    # disp-00000 是 S4 特意生成的【平衡帧】（未位移的完美超胞，拟合时用来扣残余力），
    # 它不是 phono3py 的位移帧：位移帧是 POSCAR-00001..N <-> disps[0..N-1]。
    # 把它算进来会让抽帧整体错位一帧，并在末尾 disp-(N) 上越界。
    frames = [f for f in sorted(d4.glob("disp-*")) if f.name != "disp-00000"]
    if not frames or not disps:
        return True, "无帧或无位移记录，跳过抽帧校验"
    idx = sorted({0, len(frames) // 2, len(frames) - 1})[:n_sample]
    bad = []
    for k in idx:
        if k >= len(disps):          # 目录数与位移数不一致时不要崩
            continue
        f = frames[k]
        vp = f / "vasprun.xml"
        if not vp.is_file():
            continue
        # 期望坐标 = 超胞 + 该帧位移（displacements[k] 是 {atom: 序号, displacement: 笛卡尔} 列表）
        exp = sc_cart.copy()
        # 位移有两种记录格式（2026-09-17 实测）：
        #   ① [[dx,dy,dz], ...] —— 按原子序，每原子一个矢量（kleq/rattle 的 type-2 dataset）
        #   ② [{"atom": i, "displacement": [...]}, ...] —— 只列被位移的原子
        rec = disps[k]
        if rec and isinstance(rec[0], dict):
            for it in rec:
                try:
                    exp[int(it["atom"])] += np.array(it["displacement"], dtype=float)
                except Exception:                       # noqa: BLE001
                    pass
        else:
            arr = np.array(rec, dtype=float)
            if arr.shape == exp.shape:
                exp = exp + arr
        got = _vasprun_basis_positions(vp)
        if got is None:
            bad.append("%s: vasprun.xml 解析失败" % f.name)
            continue
        b, p = np.array(got[0], dtype=float), np.array(got[1], dtype=float)
        cart = p @ b
        if cart.shape != exp.shape:
            bad.append("%s: 原子数 %s != %s" % (f.name, cart.shape, exp.shape))
            continue
        # ★ 必须按周期回绕后再比：VASP 把分数坐标写在 [0,1) 内，而"理想超胞+位移"
        #   可能落在区间外 —— 不回绕的话差值会凭空多出一整个晶格矢量（实测 12.34 Å）。
        fdiff = (cart - exp) @ np.linalg.inv(b)
        fdiff -= np.round(fdiff)
        dev = float(np.abs(fdiff @ b).max())
        if dev > tol:
            bad.append("%s: 坐标最大偏差 %.4f Å" % (f.name, dev))
    if bad:
        return False, ("抽帧校验失败（%d/%d 帧）：%s —— 力数据与位移不对应，"
                       "S5 会把两者混拟合。请清空 step4_disp 后重跑 S4。"
                       % (len(bad), len(idx), "；".join(bad[:3])))
    return True, "抽帧校验通过（%d 帧，逐原子比对 < %.0e Å）" % (len(idx), tol)


# ==========================================================================
# SCF 收敛门禁（2026-09-19）：NELM 截断 / 未收敛的帧不能拿去拟合
# ==========================================================================
# VASP 在 SCF 撞 NELM 时**照样输出力、作业正常退出**，只在 OUTCAR 留一段：
#     The electronic self-consistency was not achieved in the given
#     number of steps (NELM). The forces and other quantities evaluated ...
# 这种帧力不可信，拿去拟合会让 fc2/fc3 与 κ 整体错掉，而下游所有检查都显示正常。
# 反向的情形：静态单点正常收敛时 OUTCAR 必有且只有 1 次
#     "aborting loop because EDIFF is reached"；计数为 0 说明电子步没正常收住，
#     != 1（>1 通常是多离子步）也不符合 S4 单点帧的预期。
NELM_WARNING = "number of steps (NELM)"
EDIFF_ABORT_MARK = "aborting loop because EDIFF is reached"
_DAV_RE = re.compile(r"^\s*DAV:\s*(\d+)", re.M)
# SCF 行字段：<DAV|RMM|CG|EDDAV>: n E dE deps ncg rms [rms(c)]
_DE_RE = re.compile(r"^\s*(?:DAV|RMM|CG|EDDAV):\s*\d+\s+\S+\s+(\S+)", re.M)
_DEPS_RE = re.compile(r"^\s*(?:DAV|RMM|CG|EDDAV):\s*\d+\s+\S+\s+\S+\s+(\S+)", re.M)
# 本帧 INCAR 的 EDIFF —— 判据用它，而不是写死阈值
_EDIFF_RE = re.compile(r"^\s*EDIFF\s*=\s*([0-9.eEdD+-]+)", re.M)
# 撞 NELM 但**确实已收敛**的帧，力可用（wangchao review II.2；口径 2026-09-23 收紧）：
#   最后 SCF_DE_MIN_STEPS(5) 个电子步的 |dE| 全部 ≤ EDIFF，且末步 d_eps ≤ 10×EDIFF。
# ★ 只看末步有风险 —— 振荡的 SCF 可能恰好在某一步落到很小的 dE，所以要求"连续 5 步都达标"。
# 实测 WS₂ 6 帧：末 5 步 dE ≈ ±(0.6~1.0)e-9（EDIFF=1E-8 全达标），d_eps≈1.01e-8 ≤ 10×EDIFF；
# 它们没能触发 abort 只因 VASP 的第二判据 d_eps 比 EDIFF 高一丁点。
# 读不到 EDIFF 时回落到 SCF_DE_OK_EV（AUTOZT_SCF_DE_OK 可覆盖）并跳过 d_eps 判据。
SCF_DE_OK_EV = float(os.environ.get("AUTOZT_SCF_DE_OK", "1e-6"))
SCF_DEPS_FACTOR = 10.0
SCF_DE_MIN_STEPS = 5


def scan_frame_scf(frame_dir, de_ok_ev=SCF_DE_OK_EV):
    """只读扫一帧（disp-XXXXX/）的 OUTCAR + OSZICAR，返回诊断 dict。

    作废判据：
      - OUTCAR 含 NELM 警告（"number of steps (NELM)"）→ SCF 未在 NELM 内收敛；
      - OUTCAR 存在但 "aborting loop because EDIFF is reached" 出现次数 != 1；
      - OUTCAR 不存在 → 帧未算完，无法核验，同样作废。
    电子步数：OSZICAR 里 "DAV:" 行数与最后一个 DAV 序号（静态单点两者相等）；
    以最后一个 DAV 序号为准（兼容多离子步的极端情况），同时记录行数。
    """
    d = Path(frame_dir)
    row = {"frame": d.name, "dir": str(d), "has_outcar": False, "has_oszicar": False,
           "nelm_warn": False, "n_aborting": 0, "n_dav_lines": 0, "last_dav": 0,
           "scf_steps": 0, "dE_last": None, "dE_last5": [], "deps_last": None,
           "ediff": None, "nelm_converged": False,
           "ok": True, "reasons": []}
    outcar = d / "OUTCAR"
    if outcar.is_file():
        row["has_outcar"] = True
        try:
            txt = outcar.read_text(errors="ignore")
        except OSError as e:                                  # pragma: no cover
            txt = ""
            row["reasons"].append("OUTCAR 读取失败：%s" % e)
        row["nelm_warn"] = NELM_WARNING in txt
        row["n_aborting"] = txt.count(EDIFF_ABORT_MARK)
    else:
        row["reasons"].append("OUTCAR 不存在（帧未算完）")
    osz = d / "OSZICAR"
    if osz.is_file():
        row["has_oszicar"] = True
        try:
            _txt = osz.read_text(errors="ignore")
            idx = [int(m) for m in _DAV_RE.findall(_txt)]
            _de = [float(x) for x in _DE_RE.findall(_txt)]
            _deps = [float(x) for x in _DEPS_RE.findall(_txt)]
            row["dE_last"] = _de[-1] if _de else None
            row["dE_last5"] = _de[-int(SCF_DE_MIN_STEPS):]
            row["deps_last"] = _deps[-1] if _deps else None
            _inc = d / "INCAR"
            if _inc.is_file():
                _m = _EDIFF_RE.search(_inc.read_text(errors="ignore"))
                if _m:
                    row["ediff"] = float(_m.group(1).replace("D", "E").replace("d", "e"))
        except (OSError, ValueError):                         # pragma: no cover
            idx = []
        row["n_dav_lines"] = len(idx)
        row["last_dav"] = max(idx) if idx else 0
    row["scf_steps"] = row["last_dav"] or row["n_dav_lines"]
    # 分级：撞 NELM 但"连续 5 步 |dE| ≤ EDIFF 且 d_eps ≤ 10×EDIFF" ⇒ 力可用，放行只留痕。
    if row["ediff"] is not None:
        _tol = row["ediff"]
        _ok_de = (len(row["dE_last5"]) >= int(SCF_DE_MIN_STEPS)
                  and all(abs(x) <= _tol for x in row["dE_last5"]))
        _ok_deps = (row["deps_last"] is not None
                    and abs(row["deps_last"]) <= SCF_DEPS_FACTOR * _tol)
        _crit = "末 %d 步 |dE| ≤ EDIFF=%g 且 d_eps ≤ %g×EDIFF" % (
            SCF_DE_MIN_STEPS, _tol, SCF_DEPS_FACTOR)
    else:                                    # 读不到 EDIFF：回落旧阈值，跳过 d_eps 判据
        _tol = de_ok_ev
        _ok_de = (row["dE_last"] is not None and abs(row["dE_last"]) <= _tol)
        _ok_deps = True
        _crit = "末步 |dE| ≤ %g（未读到 EDIFF，回落阈值；未查 d_eps）" % _tol
    row["nelm_converged"] = bool(row["nelm_warn"] and _ok_de and _ok_deps)
    if row["nelm_warn"] and not row["nelm_converged"]:
        row["reasons"].append(
            "SCF 未在 NELM 内收敛（NELM 警告；不满足 %s；实测末 5 步 dE=%s、d_eps=%s）"
            % (_crit, row["dE_last5"], row["deps_last"]))
    if row["has_outcar"] and row["n_aborting"] != 1 and not row["nelm_converged"]:
        row["reasons"].append("'aborting loop because EDIFF is reached' 出现 %d 次（应为 1）"
                              % row["n_aborting"])
    row["ok"] = not row["reasons"]
    return row


def check_outcar_scf_convergence(step4_dir, n_outliers=8):
    """逐帧扫 OUTCAR 的 SCF 收敛门禁。返回 (ok, info)。

    info 含每帧明细（frames）、作废帧（bad_frames）、电子步数分布（steps_stats）
    与按步数排序的异常大值（step_outliers）。只读、不写盘；任何一帧作废则 ok=False。
    """
    d4 = Path(step4_dir)
    frames = sorted([f for f in d4.glob("disp-*") if f.is_dir()])
    rows = [scan_frame_scf(f) for f in frames]
    bad = [r for r in rows if not r["ok"]]
    steps = [r["scf_steps"] for r in rows if r["scf_steps"] > 0]
    stats = {}
    if steps:
        s = sorted(steps)
        n = len(s)
        stats = {"n": n, "min": s[0], "max": s[-1],
                 "mean": round(sum(s) / float(n), 1), "median": s[n // 2],
                 "p90": s[int(round(0.9 * (n - 1)))]}
    outliers = sorted([r for r in rows if r["scf_steps"] > 0],
                      key=lambda r: -r["scf_steps"])[:max(1, int(n_outliers))]
    info = {"step4_dir": str(d4), "n_frames": len(rows),
            "n_with_outcar": sum(1 for r in rows if r["has_outcar"]),
            "n_bad": len(bad),
            "n_nelm": sum(1 for r in rows if r["nelm_warn"]),
            "n_nelm_converged": sum(1 for r in rows if r["nelm_converged"]),
            "n_abort_ne_1": sum(1 for r in rows
                                if r["has_outcar"] and r["n_aborting"] != 1),
            "n_no_oszicar": sum(1 for r in rows if not r["has_oszicar"]),
            "steps_stats": stats,
            "step_outliers": [{"frame": r["frame"], "scf_steps": r["scf_steps"],
                               "nelm_warn": r["nelm_warn"],
                               "nelm_converged": r["nelm_converged"],
                               "dE_last": r["dE_last"],
                               "n_aborting": r["n_aborting"]} for r in outliers],
            "bad_frames": [{"frame": r["frame"], "reasons": r["reasons"],
                            "scf_steps": r["scf_steps"]} for r in bad],
            "frames": rows, "ok": not bad}
    return (not bad), info


def format_scf_report(info, top=6):
    """把 check_outcar_scf_convergence 的 info 格式化成一屏可读的日志/报告。"""
    lines = []
    st = info.get("steps_stats") or {}
    lines.append("[SCF] %s：帧 %d（有 OUTCAR %d）｜NELM 警告 %d 帧｜"
                 "'aborting loop'!=1 共 %d 帧｜无 OSZICAR %d 帧"
                 % (info.get("step4_dir", "?"), info.get("n_frames", 0),
                    info.get("n_with_outcar", 0), info.get("n_nelm", 0),
                    info.get("n_abort_ne_1", 0), info.get("n_no_oszicar", 0)))
    if st:
        lines.append("[SCF] 电子步数分布：min=%s 中位=%s mean=%s p90=%s max=%s"
                     % (st.get("min"), st.get("median"), st.get("mean"),
                        st.get("p90"), st.get("max")))
    if info.get("step_outliers"):
        lines.append("[SCF] 电子步数最多：" + "，".join(
            "%s=%s%s" % (o["frame"], o["scf_steps"],
                         "（NELM 警告!）" if o["nelm_warn"] else "")
            for o in info["step_outliers"][:max(1, int(top))]))
    if info.get("n_nelm_converged"):
        lines.append("[SCF] 其中 %d 帧撞 NELM 但末步 |dE| 已 ≤ %.0e eV —— 能量已收敛、力可用，"
                     "按放行处理（只留痕不作废；见 wangchao review II.2）。"
                     % (info["n_nelm_converged"], SCF_DE_OK_EV))
    if info.get("bad_frames"):
        lines.append("[SCF] 作废帧 %d 个：" % len(info["bad_frames"]))
        for b in info["bad_frames"][:max(1, int(top)) * 4]:
            lines.append("        %s：%s（电子步 %s）"
                         % (b["frame"], "；".join(b["reasons"]), b["scf_steps"]))
        if len(info["bad_frames"]) > top * 4:
            lines.append("        … 其余 %d 帧见 scf_steps.json"
                         % (len(info["bad_frames"]) - top * 4))
    else:
        if info.get("n_nelm_converged"):
            lines.append("[SCF] 无作废帧：%d 帧撞 NELM 但 |dE| 已收敛（分级放行），"
                         "其余帧各 1 次 aborting loop。" % info["n_nelm_converged"])
        else:
            lines.append("[SCF] 所有帧 SCF 正常收敛（各 1 次 aborting loop，无 NELM 警告）。")
    return "\n".join(lines)


def _scan_main(argv=None):
    """只读 CLI：python kl_common.py --scan-nelm <step4_disp 目录> [--json] [--out PATH]。

    退出码：0 = 全部帧 SCF 正常收敛；1 = 有帧作废（NELM / 未收敛 / OUTCAR 缺失）。
    """
    import argparse
    import json as _json
    ap = argparse.ArgumentParser(
        prog="kl_common.py",
        description="只读扫描 S4 各帧 OUTCAR 的 SCF 收敛情况（NELM 截断 / 未正常收住）。")
    ap.add_argument("--scan-nelm", metavar="STEP4_DIR", required=True,
                    help="含 disp-*/OUTCAR 的 step4_disp 目录")
    ap.add_argument("--json", action="store_true", help="打印完整 JSON（默认打印汇总表）")
    ap.add_argument("--out", metavar="PATH", help="把完整 JSON 另写到该文件（只写本地）")
    args = ap.parse_args(argv)
    ok, info = check_outcar_scf_convergence(args.scan_nelm)
    if args.out:
        Path(args.out).write_text(_json.dumps(info, ensure_ascii=False, indent=2),
                                  encoding="utf-8", newline="\n")
    if args.json:
        print(_json.dumps(info, ensure_ascii=False, indent=2))
    else:
        print(format_scf_report(info))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_scan_main())


DIPOLE_TAGS = ("LDIPOL", "IDIPOL", "DIPOL")


def dipole_line_from_incar(*paths):
    """从上游步骤的 INCAR 抄 LDIPOL/IDIPOL/DIPOL，返回 (INCAR 行列表, 来源路径)。

    P2-3：偶极修正要求 S1(弛豫)/S2(静态)/S4(取力) 处在【同一个静电边界条件】下。
    与其在三处各自判一遍（判据漂移就前功尽弃），不如让下游直接继承 S1 的写法。
    """
    for p in paths:
        p = Path(p)
        if not p.is_file():
            continue
        on, extra = False, []
        for ln in p.read_text(errors="ignore").splitlines():
            s = ln.split("#")[0].split("!")[0].strip()
            if not s:
                continue
            for k in DIPOLE_TAGS:
                if re.match(r"^%s\s*=" % k, s, re.I):
                    if k == "LDIPOL":
                        if re.search(r"\.TRUE\.|^\s*LDIPOL\s*=\s*T", s, re.I):
                            on = True
                    else:
                        extra.append(s)
        if on:
            return (["# 偶极修正：从 %s 继承（P2-3，三步必须一致）" % p.parent.name,
                     "LDIPOL = .TRUE."] + extra), str(p)
        return (["# 偶极修正：%s 未开（LDIPOL），本步保持一致" % p.parent.name]), str(p)
    return (["# 偶极修正：找不到上游 INCAR，未施加"]), None


def supercell_matrix(poscar, dim, min_len=15.0, max_multiple=6, vac_axis=2,
                     cutoff=None, margin=0.1):
    """按"每个非真空方向胞长 ≥ min_len"定对角超胞倍数 [na,nb,nc]。
    2D 真空方向恒 1。
    cutoff 给定时（三阶/二阶截断半径 Å）再加内切球判据：垂直胞高 ≥ 2×(cutoff+margin)
    （垂直胞高是周期性最短镜像距离，必须 ≥ 2×截断半径，否则镜像虚假相互作用污染力常数）。
    被 max_multiple 截断时只 WARN 不拦截（对齐 reference engine validate_user_supercell）。
    纯标准库（读 POSCAR 晶格），不依赖 ASE。"""
    lat, _ = read_poscar_cell_frac(poscar)
    lens = [_norm(lat[i]) for i in range(3)]
    vol = abs(_det3(lat))
    perp = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        area = _norm(_cross(lat[j], lat[k]))
        perp.append(vol / area if area > 1e-12 else 0.0)
    need = 2.0 * float(cutoff) + 2.0 * float(margin) if cutoff is not None else None
    reps = []
    for i in range(3):
        if dim == "2d" and i == (vac_axis if vac_axis is not None else 2):
            reps.append(1)
            continue
        n = max(1, int(math.ceil(min_len / max(lens[i], 1e-6))))
        if need is not None:
            n = max(n, int(math.ceil(need / max(perp[i], 1e-6))))
        if n > max_multiple:
            print("[WARN] 方向 %d 需 %d 倍（边长 %.3f Å，垂直胞高 %.3f Å%s），"
                  "被 MAX_MULTIPLE=%d 截断 → 该方向实际仅 %.2f Å"
                  % (i, n, lens[i], perp[i],
                     "，截断 %.2f Å 要求垂直胞高 ≥ %.2f Å" % (float(cutoff), need)
                     if need is not None else "",
                     max_multiple, max_multiple * perp[i]))
            n = max_multiple
        reps.append(n)
    return reps


def dim_str(reps):
    return " ".join(str(int(x)) for x in reps)


def mesh_str(mesh, dim, vac_axis=2):
    m = [mesh, mesh, mesh] if isinstance(mesh, int) else [int(x) for x in mesh]
    if dim == "2d":
        m[vac_axis if vac_axis is not None else 2] = 1
    elif dim == "3d" and 1 in m:
        # 3D 材料：某方向网格=1 是 2D 模板残留，phono3py 会报 "Grid symmetry is
        # broken"。自动对齐到三个方向的最大值（各向同性），避免手写 2D 网格静默出错。
        n = max(m)
        if n > 1:
            m = [n, n, n]
    return " ".join(str(int(x)) for x in m)


# ---- 2D 高对称路径（P2-1）--------------------------------------------------
# seekpath / VASPKIT-301 的路径生成只对 3D 有效：用在 2D 上会给出 Γ-A 这类 kz 方向的
# 线段（真空方向的平凡色散，白占图幅还会误导判读），而原来的回退路径又写死六方
# Γ-M-K-Γ（正方/矩形/斜方格子全错）。这里按 2D 的五种布拉维格子给固定路径，
# 全部 kz=0；坐标是"与 fc 同一个原胞"的倒格子约化坐标。
_HS_2D = {
    "hexagonal": [("Γ", (0, 0)), ("M", (0.5, 0)), ("K", (1.0 / 3, 1.0 / 3)), ("Γ", (0, 0))],
    "square": [("Γ", (0, 0)), ("X", (0.5, 0)), ("M", (0.5, 0.5)), ("Γ", (0, 0))],
    "rectangular": [("Γ", (0, 0)), ("X", (0.5, 0)), ("S", (0.5, 0.5)),
                    ("Y", (0, 0.5)), ("Γ", (0, 0))],
    "centered_rectangular": [("Γ", (0, 0)), ("X", (0.5, 0)), ("S", (0.5, 0.5)),
                             ("Y", (0, 0.5)), ("Γ", (0, 0))],
    "oblique": [("Γ", (0, 0)), ("X", (0.5, 0)), ("H", (0.5, 0.5)),
                ("C", (0, 0.5)), ("Γ", (0, 0))],
}


def vacuum_axis_in_primitive(cell, cart_axis):
    """把 Cartesian 真空轴映射成【原胞基矢】下标（面外分量占比最大的基矢）。

    primitive_matrix 可以任意排列原胞基矢（生产 MoS₂ 的 a1=(0,0,-25)，真空在第 0 基矢）。
    band_path_2d / ZA 方向判据索引的都是"原胞基矢"，把 POSCAR 的 Cartesian 下标直接传进去
    会取错方向 —— 这是 kl 链方向 bug 的根源。cart_axis=None 退化为"最长基矢即真空"。
    单一真源：kl_fc_backends._vacuum_axis_in_primitive 与
    gen_step5.1_plot_phonon._band_path 都调它（2026-09-21）。
    """
    import numpy as np
    C = np.asarray(cell, float)
    norms = np.linalg.norm(C, axis=1)
    if cart_axis is None:
        return int(np.argmax(norms))
    vhat = np.zeros(3, float)
    vhat[int(cart_axis)] = 1.0
    frac = np.abs(C @ vhat) / np.maximum(norms, 1e-300)
    return int(np.argmax(frac))


def classify_2d_lattice(cell, vac_axis=2, tol=0.02):
    """按面内两个原胞矢量判 2D 布拉维格子。cell 是 3×3（行=晶格矢量）。

    中心矩形格子的原胞一定是菱形（|a1|=|a2|，γ 不是 60/90/120），据此与斜方区分。
    """
    idx = [i for i in range(3) if i != int(vac_axis)]
    a1 = [float(x) for x in cell[idx[0]]]
    a2 = [float(x) for x in cell[idx[1]]]
    l1, l2 = _norm(a1), _norm(a2)
    cosg = sum(a1[i] * a2[i] for i in range(3)) / max(l1 * l2, 1e-12)
    cosg = max(-1.0, min(1.0, cosg))
    eq = abs(l1 - l2) / max(l1, l2) < tol
    if eq and (abs(cosg - 0.5) < 0.03 or abs(cosg + 0.5) < 0.03):
        return "hexagonal"
    if abs(cosg) < 0.03:
        return "square" if eq else "rectangular"
    return "centered_rectangular" if eq else "oblique"


def band_path_2d(cell, npoints=101, vac_axis=2):
    """2D 高对称路径（kz=0 恒成立），返回 (paths, labels, lattice)。

    paths 形如 [[[0,0,0],[0.5,0,0]], ...]，可直接喂
    phonopy.phonon.band_structure.get_band_qpoints_and_path_connections。
    """
    lat = classify_2d_lattice(cell, vac_axis)
    ax = int(vac_axis)
    hs = _HS_2D[lat]
    paths, labels = [], []
    for (n1, q1), (n2, q2) in zip(hs[:-1], hs[1:]):
        p1 = [0.0, 0.0, 0.0]
        p2 = [0.0, 0.0, 0.0]
        inplane = [i for i in range(3) if i != ax]
        p1[inplane[0]], p1[inplane[1]] = q1
        p2[inplane[0]], p2[inplane[1]] = q2
        paths.append([p1, p2])
    labels = [("$\\Gamma$" if n == "Γ" else n) for n, _ in hs]
    return paths, labels, lat


def _hex_or_ortho(lat):
    """晶格是否"正交 或 2D 六方"—— 这两种情形下实空间周期口径 N·|a_i| 已足够。

    正交：三轴两两垂直（此时倒空间式与实空间式恒等）。
    2D 六方：|a1|≈|a2|、夹角 60°/120°，且 a3 同时垂直 a1/a2。
    其余（单斜/三斜等斜方晶格）改用倒空间口径 —— 见 auto_mesh（wangchao review 一.4）。
    """
    import numpy as np
    L = np.asarray(lat, float)
    n = np.linalg.norm(L, axis=1)
    if min(n) <= 1e-9:
        return True
    c01 = float(np.dot(L[0], L[1]) / (n[0] * n[1]))
    c02 = float(np.dot(L[0], L[2]) / (n[0] * n[2]))
    c12 = float(np.dot(L[1], L[2]) / (n[1] * n[2]))
    ortho = abs(c01) < 0.02 and abs(c02) < 0.02 and abs(c12) < 0.02
    hex2d = (abs(n[0] - n[1]) / max(n[0], n[1]) < 0.01
             and abs(abs(c01) - 0.5) < 0.02
             and abs(c02) < 0.02 and abs(c12) < 0.02)
    return bool(ortho or hex2d)


def _recip_norm(lat, i):
    """倒格矢 |b_i| = 2π|a_j × a_k| / |a1·(a2×a3)|（无 2π 因子取 1 亦可，见调用处）。"""
    import numpy as np
    L = np.asarray(lat, float)
    j, k = [x for x in range(3) if x != i]
    vol = abs(float(np.dot(L[0], np.cross(L[1], L[2]))))
    if vol <= 1e-12:
        return 0.0
    return 2.0 * math.pi * float(np.linalg.norm(np.cross(L[j], L[k]))) / vol


def auto_mesh(spec, dim, vac_axis, poscar, q_len, nmin=20, dflt3d="15 15 15"):
    """q 网格：显式值照用；auto/空 时按"倒空间长度"估算（P1-1）。

    N_i = max(nmin, ceil(Q_LEN / |a_i|))：Q_LEN 是在倒空间里希望覆盖的长度（Å），
    物理含义是"每个方向至少采这么长的倒格矢"，比写死 15 15 15 靠谱得多 ——
    15 15 15 对 a≈3 Å 的 2D 材料等于只覆盖 ~47 Å，kappa 远没收敛。

    2D 只在面内估（真空轴恒 1）；3D 保持旧默认（dflt3d），避免老材料因为这一改
    突然多跑几十倍的网格 —— 3D 想用同一套判据就把 KAPPA_MESH 显式写成 auto3d。

    返回 (mesh_str, note)；note 为空表示用的是显式值。
    """
    s = str(spec if spec is not None else "").strip()
    if s and s.lower() not in ("auto", "auto3d"):
        return mesh_str(s.split(), dim, vac_axis), ""
    if dim != "2d" and s.lower() != "auto3d":
        return mesh_str(dflt3d.split(), dim, vac_axis), "（3D 用默认 %s）" % dflt3d
    lat, _ = read_poscar_cell_frac(poscar)
    ax = vac_axis if vac_axis is not None else 2
    m = []
    # ★ 口径（wangchao review 一.4）：Q_LEN 的物理含义是**实空间 Born–von Kármán 周期**
    #   下界（N·|a_i| ≥ Q_LEN），对正交/2D 六方已足够；但斜方（单斜/三斜）晶格下
    #   "实空间周期"与"倒空间采样密度"不成比例，改用倒空间口径：
    #     |b_i| / N_i ≤ 2π / Q_LEN   ⇒   N_i ≥ Q_LEN·|b_i| / 2π
    #   （正交时两式恒等；2D 六方若改用后者会让 WS₂ 从 88 变 ~101，故按 review 保持旧式。）
    _use_recip = not _hex_or_ortho(lat)
    for i in range(3):
        if dim == "2d" and i == ax:
            m.append(1)
            continue
        li = max(_norm(lat[i]), 1e-6)
        if _use_recip:
            bi = _recip_norm(lat, i)
            m.append(max(int(nmin), int(math.ceil(float(q_len) * bi / (2.0 * math.pi)))))
        else:
            m.append(max(int(nmin), int(math.ceil(float(q_len) / li))))
    _how = "斜方→倒空间口径 N_i≥Q_LEN·|b_i|/2π" if _use_recip else "正交/六方→实空间口径 N_i≥Q_LEN/|a_i|"
    return (" ".join(str(x) for x in m),
            "（auto：Q_LEN=%.0f Å，|a|=%.2f/%.2f/%.2f Å，下限 %d，%s）"
            % (float(q_len), _norm(lat[0]), _norm(lat[1]), _norm(lat[2]), int(nmin), _how))


def write_kl_params(path, **kv):
    """把 DIM/SUPERCELL/MESH/METHOD/NAC 等落盘，供后续步骤严格继承。"""
    lines = ["%s=%s" % (k.upper(), v) for k, v in kv.items() if v is not None]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def read_kl_params(path):
    d = {}
    p = Path(path)
    if p.is_file():
        for ln in p.read_text(errors="ignore").splitlines():
            s = ln.strip()
            if s and not s.startswith("#") and "=" in s:
                k, v = s.split("=", 1)
                d[k.strip().upper()] = v.strip()
    return d


# ==========================================================================
# phono3py / 通用 conda-env 命令
# ==========================================================================
def run_in_env(cmd_body, cwd, logname=None, env_src=None):
    """在 conda 环境里跑一段命令（可含多条 && 串联），可 tee 日志。返回 returncode。"""
    env = env_src if env_src is not None else PHONO3PY_ENV_SRC
    tee = " 2>&1 | tee -a %s" % logname if logname else ""
    full = "%s && cd %s && ( %s )%s" % (env, str(cwd), cmd_body, tee)
    return subprocess.run(["bash", "-lc", full]).returncode


def run_phono3py(args_str, cwd, logname="phono3py.log", env_src=None):
    print("[..] phono3py %s" % args_str)
    return run_in_env("phono3py %s" % args_str, cwd, logname, env_src)


# ==========================================================================
# VASPKIT：KPOINTS / POTCAR / ENCUT
# ==========================================================================
def vaspkit_kpoints(outdir, kscheme="2", kspacing="0.04", exe="vaspkit",
                    dim="3d", vac_axis=2, force_kz1_2d=True):
    print("[..] VASPKIT KPOINTS：1 -> 102 -> %s -> %s" % (kscheme, kspacing))
    subprocess.run([exe], input="1\n102\n%s\n%s\n" % (kscheme, kspacing),
                   text=True, cwd=outdir, check=True)
    if dim == "2d" and force_kz1_2d:
        changed, note = force_kz1(Path(outdir) / "KPOINTS",
                                  axis=vac_axis if vac_axis is not None else 2)
        print("[%s] 2D KPOINTS 真空方向细分：%s" % ("OK" if changed else "..", note))


def vaspkit_potcar(outdir, exe="vaspkit"):
    if (Path(outdir) / "POTCAR").exists():
        return
    print("[..] VASPKIT POTCAR：1 -> 103")
    subprocess.run([exe], input="1\n103\n", text=True, cwd=outdir, check=True)


def encut_from_potcar(potcar, factor=1.5, fallback=400):
    vals = []
    for ln in Path(potcar).read_text(errors="ignore").splitlines():
        m = re.search(r"ENMAX\s*=\s*([\d.]+)", ln)
        if m:
            vals.append(float(m.group(1)))
    if not vals:
        return int(fallback)
    return int(math.ceil(max(vals) * factor / 10.0) * 10)


# ==========================================================================
# 模板渲染 / 提交 / 结构接力
# ==========================================================================
def render_tpl(tpl_path, subs, out_path):
    text = Path(tpl_path).read_text(encoding="utf-8")
    for k, v in subs.items():
        text = text.replace("{{%s}}" % k, str(v))
    left = re.findall(r"\{\{([A-Z_]+)\}\}", text)
    if left:
        sys.exit("[ERROR] 模板 %s 还有未填占位符：%s"
                 % (Path(tpl_path).name, ", ".join(sorted(set(left)))))
    Path(out_path).write_text(text, encoding="utf-8", newline="\n")
    print("[OK] %s" % Path(out_path).name)


def _strip_doc_placeholders(text):
    """把说明注释里的 {{X}} 中和成 X。

    kls6：模板头部的"占位符 {{JOBNAME}} {{P3PY_CMD}}"这类说明行也会被
    全局 replace 命中，而 P3PY_CMD 是多行命令块 —— 一旦塞进注释行，它就排在
    所有 #SBATCH 之前，SLURM 遇到第一条可执行语句即停止解析指令，partition/
    ntasks/qos/output 全部失效（作业按默认队列跑，且后处理先于主程序执行）。
    真正的占位符写在 #SBATCH 行或正文里，不受影响。
    """
    out = []
    for ln in text.split("\n"):
        s = ln.lstrip()
        if s.startswith("#") and not s.startswith("#SBATCH") and "{{" in ln:
            ln = re.sub(r"\{\{(\w+)\}\}", r"\1", ln)
        out.append(ln)
    return "\n".join(out)


def _check_sbatch_order(path):
    """#SBATCH 必须全部位于第一条可执行语句之前，否则 SLURM 静默忽略。"""
    lines = Path(path).read_text(encoding="utf-8").split("\n")
    first_cmd = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        first_cmd = i
        break
    if first_cmd is None:
        return
    bad = [i + 1 for i, ln in enumerate(lines[first_cmd + 1:], start=first_cmd + 1)
           if ln.strip().startswith("#SBATCH")]
    if bad:
        sys.exit("[ERROR] %s 第 %s 行的 #SBATCH 排在可执行语句(第 %d 行)之后，"
                 "SLURM 会全部忽略它们（分区/核数/时限/输出都会失效）。"
                 "多半是多行占位符被填进了头部注释——检查模板。"
                 % (Path(path).name, ",".join(str(b) for b in bad), first_cmd + 1))


def write_submit(tpl_path, out_path, subs, require=(), label=""):
    """渲染提交模板。require 里的字面量在渲染结果里必须出现，否则硬失败。

    为什么需要：项目级 project_setting/templates/ 的优先级高于技能模板
    （见 autozt find_asset）。项目里那份副本往往是技能模板的旧拷贝 —— 旧副本
    没有新加的 {{占位符}}，渲染出来"语法完全正常、作业也能跑"，只是新功能被
    静默吃掉。enforce_incar_tags() 只能兜 INCAR，提交脚本这一侧要靠 require。
    """
    text = _strip_doc_placeholders(Path(tpl_path).read_text(encoding="utf-8"))
    for k, v in subs.items():
        text = text.replace("{{%s}}" % k, str(v))
    missing = [t for t in require if t not in text]
    if missing:
        sys.exit("[ERROR] %s提交模板 %s 渲染后缺少 %s —— 该项目/集群级的同名模板副本"
                 "遮蔽了技能模板（project_setting/templates/ 优先级最高），新功能会"
                 "被静默吃掉。请把技能模板 %s 同步覆盖过去（先备份旧副本）。"
                 % (label, Path(tpl_path).name, "、".join(missing),
                    Path(tpl_path).name))
    Path(out_path).write_text(text, encoding="utf-8", newline="\n")
    _check_sbatch_order(out_path)
def relay_poscar(prev_contcar, dst_poscar, label="上一步"):
    if not Path(prev_contcar).is_file():
        sys.exit("[ERROR] %s 的 CONTCAR 不存在：%s\n"
                 "        请确认上一步已完成再生成本步。" % (label, prev_contcar))
    bad = validate_poscar(prev_contcar)
    if bad:
        sys.exit("[ERROR] %s 的 CONTCAR 残缺（%s）：%s" % (label, bad, prev_contcar))
    shutil.copyfile(prev_contcar, dst_poscar)
    print("[OK] POSCAR ← %s" % prev_contcar)


# ---------------------------------------------------------------------------
# 结构体检第五条（B，2026-09-20）：S1 弛豫结构的数值微畸变审计 + spglib 对称化
# ---------------------------------------------------------------------------
#   S2/S3/S4 的 gen 都从 S1 的 CONTCAR 接力原胞；relay 之后调用 symmetry_gate。
#   * 默认只审计告警（SYMMETRY_AUDIT=warn）；spglib 不可用则静默跳过，不阻断 gen。
#   * SYMMETRY_SYMMETRIZE=off（默认）时【绝不】改写结构 —— 默认不改变现有行为。
#   * 对称化必须过三项确认，其中②要求对称化前后各一个单点能量
#     （SYMMETRY_ENERGY_BEFORE/AFTER 给 OUTCAR 路径）；缺省即拒绝静默通过。
#   开关写在【本步的】step.conf（templates/<step>/step.conf），不要写进全局
#   templates/step.conf（会漏进别的步骤，触发"本脚本不认识的键"）。
SYMMETRY_SPEC = {
    "SYMMETRY_AUDIT":          ("warn", "str"),   # off | warn | error
    "SYMMETRY_SYMMETRIZE":     ("off",  "str"),   # off | on（on 才改写结构）
    "SYMMETRY_AUDIT_SYMPREC":  (1e-4,   "float"),
    "SYMMETRY_ENERGY_TOL_MEV": (1.0,    "float"),
    "SYMMETRY_ENERGY_BEFORE":  (None,   "str"),   # 对称化前单点 OUTCAR 路径
    "SYMMETRY_ENERGY_AFTER":   (None,   "str"),   # 对称化后单点 OUTCAR 路径
    "SYMMETRY_STRESS_BEFORE":  (None,   "str"),
    "SYMMETRY_STRESS_AFTER":   (None,   "str"),
    "SYMMETRY_ALLOW_PENDING":  (False,  "bool"),  # true=无能量也先对称化(标 fit_for_use=false)
}


def _cget(conf, key, default=None):
    """从 StepConf / dict / None 里安全取一个键（None 或缺失返回 default）。"""
    if conf is None:
        return default
    try:
        v = conf[key]
    except (KeyError, TypeError):
        return default
    return default if v is None else v


def _switch_on(value):
    return str(value).strip().lower() in ("on", "true", "1", "yes")


def symmetry_gate(poscar, conf=None):
    """结构体检第五条：审计（可选对称化）relay 过来的原胞 POSCAR。

    返回审计结果 dict；关闭/不可用时 None。默认只在检出数值微畸变时告警，
    绝不改结构；SYMMETRY_SYMMETRIZE=on 才调用 spglib.refine_cell 就地对称化，
    且确认②（能量）未提供时拒绝静默通过（除非 SYMMETRY_ALLOW_PENDING=true）。
    """
    mode = str(_cget(conf, "SYMMETRY_AUDIT", "warn")).strip().lower()
    if mode in ("off", "false", "0", "none", ""):
        return None
    try:
        import symmetry_audit as SA
    except Exception as exc:                      # noqa: BLE001
        print("[..] 结构体检第五条跳过：symmetry_audit.py 不可用（%s）" % exc)
        return None
    if not SA.spglib_available():
        print("[..] 结构体检第五条跳过：spglib 不可用"
              "（gen 需在 phono3py/atomate2_p_a 环境里跑）")
        return None

    poscar = Path(poscar)
    if not poscar.is_file():
        print("[..] 结构体检第五条跳过：%s 不存在" % poscar)
        return None
    symprec = float(_cget(conf, "SYMMETRY_AUDIT_SYMPREC", 1e-4) or 1e-4)
    try:
        aud = SA.audit_structure_file(poscar, symprec=symprec)
    except Exception as exc:                      # noqa: BLE001
        print("[..] 结构体检第五条跳过（审计失败）：%s" % exc)
        return None
    if not aud.get("available"):
        print("[..] 结构体检第五条跳过：%s" % aud.get("reason"))
        return aud

    print(SA.format_audit_report(aud, tag="S4前"))
    sg = aud["spacegroup"]
    if not aud["micro_distortion"]:
        print("[..] 结构体检第五条：两容差空间群一致（%s #%s，%d ops），无需对称化"
              % (sg[-1]["international"], sg[-1]["number"], sg[-1]["n_ops"]))
        return aud

    print("[WARN] 结构体检第五条：检出数值微畸变 —— symprec=%g -> %s(#%s)/%d ops，"
          "symprec=%g -> %s(#%s)/%d ops（ZA 被线性化的几何根因之一）"
          % (sg[0]["symprec"], sg[0]["international"], sg[0]["number"], sg[0]["n_ops"],
             sg[-1]["symprec"], sg[-1]["international"], sg[-1]["number"], sg[-1]["n_ops"]))

    if not _switch_on(_cget(conf, "SYMMETRY_SYMMETRIZE", "off")):
        msg = ("进 S4 前建议对称化：python symmetry_audit.py --poscar %s "
               "--mode symmetrize --inplace（或在本步 step.conf 设 "
               "SYMMETRY_SYMMETRIZE=on）" % poscar)
        if mode == "error":
            sys.exit("[ERROR] 结构体检第五条：数值微畸变（SYMMETRY_AUDIT=error）。" + msg)
        print("[WARN] " + msg)
        return aud

    # ---- 对称化：确认②（能量）必须显式提供，否则拒绝静默通过 ----
    allow_pending = bool(_cget(conf, "SYMMETRY_ALLOW_PENDING", False))
    en_before = _cget(conf, "SYMMETRY_ENERGY_BEFORE")
    en_after = _cget(conf, "SYMMETRY_ENERGY_AFTER")
    if not allow_pending and (en_before is None or en_after is None):
        # 只把【候选】对称化结构留档（不改本步 POSCAR），停下来要单点能量 ——
        # 避免 gen 失败时留下一个"未确认"的结构被下游静默用掉。
        cand = poscar.parent / "symmetry_audit" / "POSCAR.symmetrized"
        try:
            SA.symmetrize_structure_file(
                poscar, inplace=False, symprec=symprec,
                allow_energy_pending=True, write_prov=False,
                source_note="kl-dft-cpu %s (candidate, energy pending)"
                            % poscar.parent.name)
        except Exception as exc:                  # noqa: BLE001
            sys.exit("[ERROR] 结构体检第五条：对称化失败（%s）" % exc)
        sys.exit("[ERROR] 结构体检第五条：检出数值微畸变，但确认②（单点能量）还没给 —— "
                 "已把【候选】对称化结构写到 %s（本步 POSCAR 未改）。\n"
                 "        请对【对称化前/后】各跑一个单点，再把两个 OUTCAR 路径写进本步 "
                 "step.conf 的 SYMMETRY_ENERGY_BEFORE / SYMMETRY_ENERGY_AFTER，然后 retry。\n"
                 "        确实要先出结构再补单点，可设 SYMMETRY_ALLOW_PENDING=true"
                 "（结果标 fit_for_use=false）。" % cand)
    res = SA.symmetrize_structure_file(
        poscar,
        inplace=True,
        symprec=symprec,
        energy_before=en_before,
        energy_after=en_after,
        energy_tol_mev_per_atom=float(_cget(conf, "SYMMETRY_ENERGY_TOL_MEV", 1.0) or 1.0),
        stress_before=_cget(conf, "SYMMETRY_STRESS_BEFORE"),
        stress_after=_cget(conf, "SYMMETRY_STRESS_AFTER"),
        allow_energy_pending=allow_pending,
        source_note="kl-dft-cpu %s" % poscar.parent.name,
    )
    if res.get("triggered"):
        c2 = (res.get("confirmations") or {}).get("c2_energy", {}).get("status")
        print("[OK] 结构体检第五条：已对称化（used_symprec=%g，能量确认=%s，"
              "fit_for_use=%s），原结构备份 %s.pre_symmetry"
              % (res.get("used_symprec", symprec), c2,
                 res.get("fit_for_use"), poscar.name))
    return res


def find_prev_dir(cwd, candidates):
    for name in candidates:
        d = Path(cwd) / name
        if (d / "CONTCAR").is_file() or (d / "POSCAR").is_file():
            return d
    return None


def new_jobname(cwd, step_label):
    return "%s-kl-dft-cpu-%s" % (Path(cwd).name, step_label)


def resolve_submit(base_dir, dim, kind="submit_std"):
    """按维度找提交模板：<kind>_<dim>.tpl，回退 <kind>.tpl / resolve_tpl。"""
    base = Path(base_dir)
    for name in ("%s_%s.tpl" % (kind, dim), "%s.tpl" % kind):
        if (base / name).is_file():
            return base / name
    try:
        return resolve_tpl(base, kind, dim)
    except SystemExit:
        sys.exit("[ERROR] 找不到 %s 的提交模板（%s_%s.tpl）" % (dim.upper(), kind, dim))
