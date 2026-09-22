#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_step6_kappa.py —— BTE 晶格热导率，提交计算节点（step6_kappa）。

从 step5_fc 拷 fc2/fc3/phono3py_disp.yaml/BORN，按 kl_params 的 MESH 组 BTE 命令，
渲染提交模板 → submit.sh，tf 提交到计算节点。成功后把 κ 张量写进 kappa_summary.json
并落 KAPPA_DONE（marker 判据）。
求解器（step.conf 的 SOLVER）：
  phono3py   : phono3py-load --br（完整支持 findiff/alm + NAC，默认）
  shengbte   : 写 ShengBTE CONTROL（复用参考引擎例程）。注意 fc3→ShengBTE 导出仅
               random/hiphive 路线可靠，findiff 的 compact fc3 无稳定导出口——solver=shengbte
               建议配 METHOD=alm，且需在集群装好 ShengBTE、把 exe 填进 step.conf。
  fourphonon : ShengBTE 同源引擎，跑 3ph RTA（four_phonon=F，不需 fc4）。输入格式与
               shengbte 相同（FORCE_CONSTANTS_2ND/3RD 拷自 S5_fc/shengbte/，需
               EXPORT_SHENGBTE=true）。适合 GPU 机多卡加速（FOURPHONON_NGPU 张卡、
               rank=卡、acc_set_device_num 自动分卡）；CPU 机请用 shengbte 别用
               fourphonon（CPU 单进程枚举慢、AOCC 版易卡死，见 README）。
产出目录：step6_kappa/
"""
import json
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kl_common as kc
import stepconf
# 2D 层厚/归一化唯一真源（与 ke-dft-cpu 的 AMSET 步共用同一份口径；
#   见 skill/_common/thickness_2d.py，随 gen_need 推送到本步目录）
from thickness_2d import slab_geometry

OUTDIR = "step6_kappa"
STEP   = "step6_kappa"
FC_DIR = "step5_fc"

SPEC = {
    "FUNC":        ("pbesol", "str"),  # 全局 step.conf 带入，本步不用
    "SOLVER":       ("phono3py", "str"),  # phono3py | shengbte（= 热导率计算软件）
    "BTE_METHOD":   ("rta",     "str"),   # phono3py 路：rta(--br) | lbte(--lbte)
    # phono3py 的单节点并行（2026-09-18 实测）：本环境的 phono3py 是 **OpenMP-only**
    #   —— C 扩展只链 libgomp，包内 0 处 mpi4py/MPI，所以不能 mpirun（那会跑 N 份全量
    #   进程互相覆盖输出）。唯一有效并行 = 1 进程 × N 个 OMP 线程，N 建议 = 节点物理核数
    #   （jzzn 192 逻辑核 = 96 物理核）。改这一个键即可调；0/空回落 96。
    "P3PY_OMP_THREADS": (96,    "int"),
    "MESH_OVERRIDE": (None,     "str"),   # 空=用 step4 写入 kl_params 的 MESH
    # q 网格收敛扫描（P1-1）：""=只跑一套 | auto=三档(N/1.25N/1.5625N) | "a b c; d e f"
    #   判据：相邻档 300K 面内 κ 变化 < 5%（写进 kappa_summary.json 的 mesh_convergence）
    "MESH_SCAN":    ("auto",    "str"),   # auto=2D 三档扫描/3D 单套；""=强制单套
    "MESH_MIN":     (20,        "int"),   # auto 时每方向下限
    # RTA vs 完整解对照（P1-2）：auto = 2D 跑一次 --lbte（正规过程主导，RTA 会低估 κ）。
    #   两种解法的输出文件名相同，需跑完 RTA 先改名再跑 LBTE；比值写进 rta_over_full。
    "COMPARE_LBTE": ("auto",    "str"),   # auto | on | off
    "T_MIN":        (100,       "int"),
    "T_MAX":        (800,       "int"),
    "T_STEP":       (100,       "int"),
    "ISOTOPE":      (True,      "bool"),
    "SCALEBROAD":   (1.0,       "float"), # shengbte/fourphonon 展宽（ShengBTE 默认 1.0）
    "SHENGBTE_EXE": ("ShengBTE", "str"),
    # ShengBTE 3ph 求解：True=ShengBTE 默认的迭代自洽解（写 BTE.KappaTensorVsT_CONV），
    # False=只算 RTA（写 ..._sg）。迭代段在 128 原子/7x7x7 这种大胞上极贵：实测 Mg8C120
    # 96 核跑 1000 次迭代 15.1 h 后 SIGSEGV（CONV 全 NaN、白烧 13 h），而 RTA 1.5 h 就出全
    # 温区结果。RTA-only 时 kappa_summary.json 会自动以 ..._sg 为 source（模板已回退）。
    "KAPPA_CONVERGENCE": (True, "bool"),
    # ---- ShengBTE 的 MPI/OMP 布局（见 submit_shengbte.tpl 文件头的两个坑）----
    # ShengBTE 靠 MPI 按 q 点并行，mpirun -n 1 = 串行（实测 11.8 h 零产物）。
    # 正确布局 = 一个 MPI rank 占一个 NUMA 域：
    #     NTASKS = TOTAL_CORES / CORES_PER_NUMA，CPUS_PER_TASK = CORES_PER_NUMA
    # jzzn 计算节点（cpu192）实测：192 核 = 2 socket x 96，**8 个 NUMA 节点 x 24 核**
    #   （lscpu / numactl -H；早先误记为 24 NUMA x 8 核）。mpirun 的 --map-by numa 是
    #   "1 rank 占 1 个 NUMA 节点"：rank 数 > NUMA 节点数时会循环复用同一节点、把多个
    #   rank 绑到同一组核上 —— 实测 24 rank 时 3 rank 挤 8 核，每线程只拿到 33% 的核，
    #   整作业只用 64/192 核（其余 128 核闲置）。故 CORES_PER_NUMA 必须是**真实每 NUMA
    #   核数 24**：192 核 -> 8 rank x 24 线程；96 核 -> 4 rank x 24 线程。
    #   换集群只改这两个数，模板不用动。SHENGBTE_NTASKS 留 "auto" 即按上面公式算。
    "SHENGBTE_TOTAL_CORES":    (96,     "int"),
    "SHENGBTE_CORES_PER_NUMA": (24,     "int"),
    "SHENGBTE_NTASKS":         ("auto", "str"),
    # SLURM QoS：默认 premium（每用户 5 个并发作业）。jzzn 上 premium 槽位常被本账号
    # 其它技能占满，S6_kappa 会以 QOSMaxJobsPerUserLimit 排队等很久；regular 允许
    # 每用户 50 个并发、墙钟 1 天，RTA 级别的 κ 计算足够。见 submit_shengbte.tpl。
    "SBATCH_QOS":              ("premium", "str"),
    "FOURPHONON_EXE": ("", "str"),          # fourphonon(multi-GPU) 可执行文件绝对路径
    "FOURPHONON_NGPU": (4, "int"),          # fourphonon 用几张 GPU(rank=卡)
    "FOURPHONON_CPUS_PER_GPU": (8, "int"),  # 每 GPU 配几个 CPU 核(cpus-per-task+OMP)
    # ★ 安全闸：FourPhonon v1.3 **GPU** 版对 phonopy/ShengBTE 格式 fc2 有官方 OpenACC
    #   数值 bug（κ 错 67~71×，只有 espresso 路径正确、官方未修）；本技能喂的正是
    #   phonopy 格式，故 SOLVER=fourphonon 默认**拒绝**。确知风险仍要跑（对照/复现 bug）
    #   时显式置 true。详见 skill/kl-dft-cpu/README.md 的「fourphonon」节。
    "ALLOW_FOURPHONON_GPU_PHONOPY": (False, "bool"),
    # 集群 conda.sh（tf 从 setting/<集群>.yaml 的 conda_sh 注入 step.conf，切集群自动跟着走）
    "CONDA_SH": ("", "str"),
    # 2D κ 厚度归一化：phono3py 用含真空的原胞体积做分母，2D 面内 κ 被胞高稀释，
    #   需乘 h⊥/d（h⊥=V/A 周期胞高，d=层有效厚度）。
    #   取法：vdw=原子层跨度+两侧vdW半径 | cell=用 h⊥(即不归一) | 数值=固定Å
    "KAPPA_2D_THICKNESS": ("vdw",  "str"),
    # 2D NAC 覆盖：auto=2D默认不用3D-NAC(LO-TO在2D应趋零)/3D随BORN；on=强制用；off=强制不用
    "KAPPA_NAC":          ("auto", "str"),
}

# 层厚/归一化口径全部下沉到 skill/_common/thickness_2d.py（vdw 半径表在
#   skill/_common/vdw_radii.py，ke-dft-cpu 读同一份）——本文件不再自带半径表副本。

def _poscar_species(poscar):
    """从 POSCAR 读每个原子的元素符号（VASP5）；VASP4 无元素行则返回 None。"""
    import re as _re
    L = Path(poscar).read_text(encoding="utf-8-sig").splitlines()
    line6 = L[5].split()
    if not line6 or _re.fullmatch(r"[+-]?\d+", line6[0]):
        return None
    counts = [int(x) for x in L[6].split()]
    out = []
    for sym, n in zip(line6, counts):
        out += [sym] * n
    return out


def two_d_norm_factor(poscar, vac_axis, mode):
    """2D κ 厚度归一化因子 factor=h⊥/d 与元数据。mode: vdw | cell | 数值(Å)。

    真源是 skill/_common/thickness_2d.slab_geometry（ke 侧读同一份），本函数只负责
    从 POSCAR 取晶格 + 元素表再转成 kappa_summary 的字段名。

    修正的三个点（旧实现）：
      ① 旧版用 Lz=|c| 当分母基准；phono3py/ShengBTE 的分母其实是 V 与面内面积 A，
         所以基准必须是 h⊥=V/A。c 轴倾斜（α/β≠90°）时 |c|≠V/A，因子会算错。
      ② 旧版直接取 frac 投影的 max−min；层跨越 z=0/1 周期边界时 span 被算成接近
         整个胞高 → d>h⊥ → 因子<1（把 κ 又缩小一遍）。现按"最大间隙=真空"切口展开。
      ③ 旧版没把因子交给 shengbte/fourphonon 模板（那两条路只写原始 κ）。
    """
    from dim_common import read_poscar_cell_frac
    lat, frac = read_poscar_cell_frac(poscar)
    ax = vac_axis if vac_axis is not None else 2
    sp = _poscar_species(poscar)
    if sp is not None and len(sp) != len(frac):
        sp = None
    g = slab_geometry(lat, frac, sp, vac_axis=ax, mode=mode)
    factor = float(g["kappa_2d_norm_factor"])
    meta = {k: g[k] for k in ("h_perp_A", "Lz_A", "atomic_span_A", "vacuum_gap_A",
                              "thickness_d_A", "thickness_convention",
                              "kappa_2d_norm_factor")}
    meta["note"] = ("kappa_2d_normalized = kappa_raw * h_perp/d"
                    "（h_perp=V/A；面内分量才有物理意义，2D 的 zz 分量无物理意义）")
    # ★ 兼容字段（不要删）：ke-dft-cpu/step8.1_boltztrap/gen_step11_boltztrap.py 直接读
    #   kappa_summary.json 的 Lz_ang / thickness_d_ang / thickness_convention 做两条链的
    #   元胞与厚度一致性闸门，gen_step13_output 也读 Lz_ang。这里的 Lz_ang 保持旧含义
    #   =|c|（闸门要跟 ke 的 cell_c_A 比），h⊥ 另有 h_perp_A 字段，两者不再混用。
    meta["Lz_ang"] = g["Lz_A"]
    meta["thickness_d_ang"] = g["thickness_d_A"]
    meta["atomic_zspan_ang"] = g["atomic_span_A"]
    return factor, meta


# 材料级 thickness_2d.json（文档 P0-2(b) 的 kl/ke 共用契约）：谁先跑谁写，对方读它做对照。
#   作业 cwd 是 <材料>/<步骤>，材料根 = ".."；先确认 ".." 里真有别的步骤目录（说明
#   它就是 gen 运行的那个步骤根）才写，避免路径猜错写到别处。best-effort，失败不拦作业。
_MAT_LEVEL = [
    "if THICK2D:\n",
    "    try:\n",
    "        import glob as _g, os as _os, re as _re, sys as _s\n",
    # ★ 材料根 = 技能目录的上一级（与 S1 的 write_material_thickness 完全对齐）：
    #   作业 cwd = <材料>/<技能>/<步骤>；'..' = 技能目录，'../..' = 材料根。
    #   老代码把 '..'（技能目录）当材料根，写到 <材料>/<技能>/thickness_2d.json，
    #   与 S1 写到 <材料>/thickness_2d.json 不一致（ke-dft-cpu 的 AMSET 步读材料根，
    #   会读不到 kl 链的层厚 → zT 里 d 约不掉）。2026-09-21 修。
    "        _skill = _os.path.abspath('..')\n",
    "        _mat = _os.path.dirname(_skill)\n",
    "        if (_re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)+', _os.path.basename(_skill))\n",
    "                and (_os.path.isfile(_os.path.join(_mat, 'POSCAR'))\n",
    "                     or _g.glob(_os.path.join(_mat, '*', '*', 'step*')))):\n",
    # thickness_2d.py 被 gen_need 推到 <材料>/<技能>/（= '..'），不在 cwd 里，
    # 必须把技能目录加进 sys.path，否则 ModuleNotFoundError 被 except 吞成 [WARN]。
    "            _s.path.insert(0, _skill)\n",
    "            from thickness_2d import write_material_level\n",
    "            print(write_material_level(_mat, THICK2D))\n",
    "        else:\n",
    "            print('[..] 材料根识别失败（skill=%s mat=%s），跳过材料级 thickness_2d.json'\n",
    "                  % (_skill, _mat))\n",
    "    except Exception as _e:\n",
    "        print('[WARN] 材料级 thickness_2d.json 写入跳过: %s' % _e)\n",
]

def build_extract(factor, meta, thick2d=None, plan=None, primary=None):
    """phono3py 跑完后就地抽 κ 到 kappa_summary.json（计算节点 conda 里执行）。

    factor!=1 时额外写 kappa_2d_normalized_*（原始 κ × h⊥/d）。2D 里 zz 分量没有
    物理意义，所以面内 xx/yy 单独也出一份 kappa_*_inplane_xx_yy 方便下游直接用
    （kappa_xx_yy_zz 保持原样，不动下游已有读取口径）。
    thick2d 非空时再落一份 thickness_2d.json（本步算出的层厚口径留档，
    供 ke-dft-cpu 的 AMSET 步对照/复用同一口径）。
    """
    import json as _json
    plan = plan or []
    body = [
        "import glob, json, h5py, numpy as np",
        "FACTOR=%r" % float(factor),
        "META=json.loads(%r)" % _json.dumps(meta, ensure_ascii=False),
        "THICK2D=json.loads(%r)" % _json.dumps(thick2d or {}, ensure_ascii=False),
        "PLAN=json.loads(%r)" % _json.dumps(plan, ensure_ascii=False),
        "PRIMARY=%r" % (primary or ""),
        "if THICK2D:",
        "    json.dump(THICK2D,open('thickness_2d.json','w'),ensure_ascii=False,indent=2)",
    ] + [ln.rstrip("\n") for ln in _MAT_LEVEL] + [
        "def _read(f):",
        "    with h5py.File(f, 'r') as h:",
        "        T = np.array(h['temperature']); K = np.array(h['kappa'])",
        "    # phono3py 的 kappa 数据集是 (n_T, 6) Voigt，这里只取对角 xx/yy/zz",
        "    raw = [[float(K[i,0]), float(K[i,1]), float(K[i,2])] for i in range(len(T))]",
        "    return [float(t) for t in T], raw",
        "d = {'KAPPA_DONE': False}",
        "d.update(META)",
        "if THICK2D: d['thickness_2d_source'] = 'kl-dft-cpu:step6_kappa'",
        "runs, missing, used = [], [], set()",
        "for p in PLAN:",
        "    _f0 = p['file']",
        "    _cand = [x for x in glob.glob(_f0) if x not in used]",
        "    if not _cand:",
        "        # 文件名兜底：主文件是 kappa-m<digits>.hdf5，别把 -sg/-gv/mfp 当它",
        "        _g = [x for x in glob.glob('kappa-m*.hdf5') if x not in used",
        "              and 'mfp' not in x and '-sg' not in x and '-gv' not in x]",
        "        _d = [x for x in _g if p['digits'] in x]",
        "        if p.get('tagged'):",
        "            _d = [x for x in _d if x.endswith('.' + p['method'] + '.hdf5')]",
        "        _cand = _d or _g",
        "    if not _cand:",
        "        missing.append(_f0)",
        "        d.setdefault('errors', []).append(",
        "            '缺 %s（%s / %s）：该次 phono3py 没写出 kappa，看 phono3py_kappa.log'",
        "            % (_f0, p['mesh'], p['method']))",
        "        continue",
        "    f = sorted(_cand)[0]",
        "    used.add(f)",
        "    T, raw = _read(f)",
        "    j = int(np.argmin(np.abs(np.array(T) - 300.0)))",
        "    rec = {'file': f, 'mesh': p['mesh'], 'method': p['method'], 'temperatures': T,",
        "           'kappa_xx_yy_zz': raw,",
        "           'kappa_inplane_xx_yy': [[v[0], v[1]] for v in raw],",
        "           'kappa_300K_xx_yy_zz': raw[j],",
        "           'kappa_inplane_300K_xx_yy': [raw[j][0], raw[j][1]],",
        "           'kappa_inplane_300K': 0.5 * (raw[j][0] + raw[j][1])}",
        "    if abs(FACTOR - 1.0) > 1e-9:",
        "        rec['kappa_2d_normalized_xx_yy_zz'] = [[v * FACTOR for v in r] for r in raw]",
        "        rec['kappa_2d_normalized_inplane_xx_yy'] = [[v[0]*FACTOR, v[1]*FACTOR] for v in raw]",
        "        rec['kappa_2d_normalized_300K_xx_yy_zz'] = [v * FACTOR for v in raw[j]]",
        "        rec['kappa_2d_normalized_inplane_300K'] = rec['kappa_inplane_300K'] * FACTOR",
        "    runs.append(rec)",
        "d['runs'] = runs",
        "d['mesh_run_count'] = len(PLAN)",
        "d['meshes'] = [p['mesh'] for p in PLAN]",
        "# 顶层字段仍取【主口径】那一档（step.conf 的 BTE_METHOD + 最细网格），",
        "#   这样 ke-dft-cpu / te-screen 读 kappa_xx_yy_zz 的口径完全不变。",
        "prim = next((r for r in runs if (r['mesh'] + '|' + r['method']) == PRIMARY), None)",
        "if prim is None and runs: prim = runs[-1]",
        "if prim is not None:",
        "    for k in ('file','mesh','temperatures','kappa_xx_yy_zz','kappa_inplane_xx_yy',",
        "              'kappa_300K_xx_yy_zz','kappa_inplane_300K_xx_yy','kappa_inplane_300K',",
        "              'kappa_2d_normalized_xx_yy_zz','kappa_2d_normalized_inplane_xx_yy',",
        "              'kappa_2d_normalized_300K_xx_yy_zz','kappa_2d_normalized_inplane_300K'):",
        "        if k in prim: d[k] = prim[k]",
        "    d['bte_method'] = prim['method']",
        "# ---- 网格收敛（P1-1）：同方法相邻档的 300K 面内 κ 相对变化 ----",
        "_by = {}",
        "for r in runs: _by.setdefault(r['method'], []).append(r)",
        "for _m, _rs in _by.items():",
        "    if len(_rs) < 2: continue",
        "    _rs = sorted(_rs, key=lambda r: [int(x) for x in r['mesh'].split()])",
        "    _kv = [r['kappa_inplane_300K'] for r in _rs]",
        "    _rel = [abs(_kv[i+1]-_kv[i])/abs(_kv[i])*100.0",
        "            for i in range(len(_kv)-1) if abs(_kv[i]) > 1e-12]",
        "    d.setdefault('mesh_convergence', {})[_m] = {",
        "        'meshes': [r['mesh'] for r in _rs], 'kappa300_inplane': _kv,",
        "        'rel_change_pct': _rel, 'max_rel_change_pct': (max(_rel) if _rel else None),",
        "        'converged_5pct': bool(_rel) and max(_rel) < 5.0}",
        "    if _rel and max(_rel) >= 5.0:",
        "        d.setdefault('warnings', []).append(",
        "            '%s 网格还没收敛：相邻档 300K 面内 κ 变化 %.1f%% ≥ 5%%，请继续加密网格'",
        "            % (_m, max(_rel)))",
        "# ---- RTA vs 完整解（P1-2）：正规过程主导的 2D 里 RTA 会低估 κ ----",
        "_rta = [r for r in runs if r['method'] == 'rta']",
        "_full = [r for r in runs if r['method'] == 'lbte']",
        "_common = sorted({r['mesh'] for r in _rta} & {r['mesh'] for r in _full},",
        "                 key=lambda s: [int(x) for x in s.split()])",
        "if _common:",
        "    _m0 = _common[0]",
        "    _a = next(r for r in _rta if r['mesh'] == _m0)",
        "    _b = next(r for r in _full if r['mesh'] == _m0)",
        "    _ra, _rb = _a['kappa_inplane_300K'], _b['kappa_inplane_300K']",
        "    d['rta_over_full'] = {'mesh': _m0, 'kappa300_inplane_rta': _ra,",
        "                          'kappa300_inplane_full': _rb,",
        "                          'ratio': (_ra / _rb if abs(_rb) > 1e-12 else None),",
        "                          'note': 'ratio<1 说明 RTA 低估 κ（2D 里正规过程主导时很常见）'}",
        "    if _rb and _ra / _rb < 0.9:",
        "        d.setdefault('warnings', []).append(",
        "            'RTA/LBTE = %.3f，RTA 明显低估，生产值建议用 lbte' % (_ra / _rb))",
        "# 成功判据：主口径（PRIMARY 的 method，通常是 RTA）的各档网格齐全即可；",
        "#   对照腿（LBTE）缺失只告警、不判失败 —— 否则 LBTE 太重/超时会把整步拖成 NO_KAPPA。",
        "_pmethod = PRIMARY.split('|')[1]",
        "_need_files = {p['file'] for p in PLAN if p['method'] == _pmethod}",
        "_missing_primary = [f for f in missing if f in _need_files]",
        "d['KAPPA_DONE'] = bool(runs) and prim is not None and not _missing_primary",
        "if missing and d['KAPPA_DONE']:",
        "    d.setdefault('warnings', []).append(",
        "        '对照运行缺失（%s）：主口径 %s 仍有效，但该项未验证' % (', '.join(missing), _pmethod))",
        "json.dump(d, open('kappa_summary.json','w'), ensure_ascii=False, indent=2)",
        "for _e in d.get('errors', []): print('  ' + _e)",
        "for _w in d.get('warnings', []): print('[WARN] ' + _w)",
        "print('KAPPA_DONE' if d['KAPPA_DONE'] else 'NO_KAPPA')",
    ]
    return "python - <<'PY'\n" + "\n".join(body) + "\nPY"


def _prim_lattice(out):
    """从 step6 目录的 phono3py_disp.yaml 读原胞晶格（3x3 行主序）；失败返回 None。"""
    p = Path(out) / "phono3py_disp.yaml"
    if not p.is_file():
        return None
    try:
        import yaml
        d = yaml.safe_load(p.read_text(encoding="utf-8"))
        return d["primitive_cell"]["lattice"]
    except Exception:                              # noqa: BLE001
        return None


def _mesh_prim_remap(mesh_str, poscar_lat, prim_lat, cart_vac_axis):
    """把 [POSCAR/Cartesian 轴序] 的网格串重排成 [原胞基矢轴序]（phono3py --mesh 的口径）。

    primitive_matrix 会重排原胞基矢（生产 MoS₂：真空从 Cartesian c 变成原胞基矢 0）。
    phono3py 的 --mesh 索引的是原胞基矢 —— 不重排就会出现"两个面内对称相关轴网格不等"
    （88 88 1 把 88 给真空基矢 a1 与面内 a2、1 给面内 a3），phono3py 直接
    `Grid symmetry is broken`（2026-09-21 MoS₂ 实测）。
    规则：原胞真空基矢恒 1；其余原胞基矢取"与它最平行的 POSCAR 非真空轴"的网格值。
    """
    import numpy as np
    vals = [int(x) for x in str(mesh_str).split()]
    if len(vals) != 3:
        return mesh_str
    P = np.asarray(prim_lat, float)
    L = np.asarray(poscar_lat, float)
    prim_vac = kc.vacuum_axis_in_primitive(P, cart_vac_axis)
    out = [1, 1, 1]
    used = set()
    if cart_vac_axis is not None:
        used.add(int(cart_vac_axis))
    for i in range(3):
        if i == prim_vac:
            out[i] = 1
            continue
        best_j, best_cos = None, -1.0
        for j in range(3):
            if j in used:
                continue
            den = float(np.linalg.norm(P[i]) * np.linalg.norm(L[j]))
            c = abs(float(np.dot(P[i], L[j])) / den) if den else 0.0
            if c > best_cos:
                best_cos, best_j = c, j
        if best_j is None:
            best_j = min(range(3),
                         key=lambda j: abs(float(np.linalg.norm(P[i]) - np.linalg.norm(L[j]))))
        out[i] = vals[best_j]
        used.add(best_j)
    return " ".join(str(x) for x in out)


def meshes(conf, params, dim, vac_axis):
    """本步要跑的 q 网格列表（P1-1）。

    MESH_SCAN 未设/auto → 2D 自动三档收敛扫描，3D 单套（旧行为）；
                          网格收敛性没法从单次结果看出来，而 2D 的 ZA 支在 Γ 附近
                          发散最慢、恰恰最需要这个判据，所以默认值按维度分。
    MESH_SCAN = ""     → 显式只跑 kl_params（S4 按 Q_LEN 估好的）那一套
    MESH_SCAN = on/auto→ 强制三档：N、ceil(1.25N)、ceil(1.5625N)（2D 真空轴恒 1）
    MESH_SCAN = off    → 强制单套（同 ""）
    MESH_SCAN = "a b c; d e f" → 显式多套
    代价：三档 = 3 次完整 κ 求解（3D 不受影响）。
    """
    import math
    ax = vac_axis if vac_axis is not None else 2
    base = kc.mesh_str((conf["MESH_OVERRIDE"] or params.get("MESH") or "20 20 20").split(),
                       dim, ax)
    scan = conf["MESH_SCAN"]
    scan = "auto" if scan is None else str(scan).strip()
    if scan == "":                       # 显式关掉（step.conf 写 MESH_SCAN = ）
        return [base]
    if scan.lower() in ("auto",):
        # 默认 auto：只有 2D 展开扫描；3D 保持单套，不动既有 3D 项目的成本
        if dim != "2d":
            return [base]
        scan = "on"
    if scan.lower() in ("on", "auto3d", "auto"):
        n = [int(x) for x in base.split()]
        out = []
        for f in (1.0, 1.25, 1.5625):
            # 第一档就是 step4 定的那套（原样，不再被 MESH_MIN 抬高 —— 扫描要以
            #   配置的网格为锚点）；只有加密档才吃 MESH_MIN 下限。
            m = [1 if (dim == "2d" and i == ax)
                 else (int(v) if f == 1.0
                       else max(int(conf["MESH_MIN"]), int(math.ceil(v * f))))
                 for i, v in enumerate(n)]
            out.append(" ".join(str(x) for x in m))
        # 去重 + 保持升序（网格小的时候 1.25 倍可能四舍五入撞档）
        seen, uniq = set(), []
        for m in out:
            if m not in seen:
                seen.add(m)
                uniq.append(m)
        return uniq
    return [kc.mesh_str(s.split(), dim, ax) for s in scan.split(";") if s.strip()]


def extract_plan(plan, tagged):
    """PLAN：extract 要读哪几个 kappa 文件、分别是哪档网格/哪种解法。
    tagged=True 时文件名是 kappa-m<digits>.<method>.hdf5（见 build_phono3py_cmd）。"""
    out = []
    for m, meth in plan:
        d = "".join(str(m).split())
        f = ("kappa-m%s.%s.hdf5" % (d, meth)) if tagged else ("kappa-m%s.hdf5" % d)
        out.append({"file": f, "mesh": m, "method": meth, "digits": d, "tagged": bool(tagged)})
    return out


def norm_subs(factor, meta, thick2d=None):
    """shengbte / fourphonon 模板的 2D 归一化占位符（那两条路原先只写原始 κ）。"""
    return {"KAPPA_2D_FACTOR": repr(float(factor)),
            "KAPPA_2D_META": json.dumps(meta, ensure_ascii=False),
            "KAPPA_2D_THICK2D": json.dumps(thick2d or {}, ensure_ascii=False)}


def build_phono3py_cmd(plan, ts, isotope, use_nac, extract, tagged=False, ts_override=None):
    """按 plan=[(mesh, method)] 依次跑 phono3py，最后一次收尾做 extract。

    tagged=True（同一网格既跑 RTA 又跑 LBTE 时）：每跑完一次就把主 kappa 文件改名成
    kappa-m<digits>.<method>.hdf5 —— phono3py 两种解法的输出文件名完全相同，不改名
    后者会直接覆盖前者，对照就白跑了。文件名仍以 kappa-m 开头，不影响 skill.yaml 里
    声明的 outputs 通配。
    """
    steps = []
    for _mesh, _method in plan:
        _digits = "".join(str(_mesh).split())
        _flag = "--lbte" if str(_method).lower() == "lbte" else "--br"
        _nac = "" if use_nac else " --nonac"
        # 逐方法温度覆盖：LBTE 对照只跑 300 K（见调用处），其余走 ts。
        _ts = (ts_override or {}).get(str(_method).lower(), ts)
        # fc2/fc3 已在 step5_fc 拟好并拷到本目录，phono3py-load 默认读 cwd 的
        # fc2.hdf5/fc3.hdf5（--no-read-fc2/--no-read-fc3 关闭），不会再从 disp.yaml 重拟。
        steps.append('phono3py-load phono3py_disp.yaml %s --mesh %s --ts="%s"%s%s '
                     '2>&1 | tee phono3py_kappa.log'
                     % (_flag, _mesh, _ts, " --isotope" if isotope else "", _nac))
        if tagged:
            steps.append('[ -f "kappa-m%s.hdf5" ] && mv "kappa-m%s.hdf5" '
                         '"kappa-m%s.%s.hdf5" || true'
                         % (_digits, _digits, _digits, _method))
    return "\n".join(steps + [extract])


def _unused_build_phono3py_cmd(mesh, ts, isotope, use_nac, extract, bte="rta"):
    method = "--lbte" if str(bte).lower() == "lbte" else "--br"
    # NAC：有 nac_params/BORN 就默认启用；无 --nac 开关（会被当 --nac-method），要关才 --nonac。
    #   phono3py 3.24/4.x 行为一致。use_nac=True 就默认带上（不加开关），否则显式 --nonac。
    nac_flag = "" if use_nac else " --nonac"
    # fc2/fc3 已在 step5_fc 由 symfc/alm 拟好并拷到本目录（fc2.hdf5/fc3.hdf5）。
    # 用 phono3py-load（phono3py>=3.x 的载入入口）：它默认读 cwd 的 fc2.hdf5/fc3.hdf5
    # （--no-read-fc2/--no-read-fc3 关闭），不会再从 disp.yaml 重新拟合（否则报
    # "Forces were not found"——forces 已在上一步拟进 hdf5，不在 disp.yaml 里）。
    p3 = ('phono3py-load phono3py_disp.yaml %s --mesh %s --ts="%s"%s%s'
          % (method, mesh, ts, " --isotope" if isotope else "", nac_flag))
    return "%s 2>&1 | tee phono3py_kappa.log\n%s" % (p3, extract)


def prepare_shengbte(cwd, out, sb_src, conf, mesh, use_nac):
    """S5 已把 FORCE_CONSTANTS_2ND/3RD 导到 step5_fc/shengbte/；这里拷进 step6，
    再按本步运行参数（T/mesh/scalebroad/NAC）写 CONTROL。CONTROL 依赖运行参数，
    故归 S6 生成，不在 S5 固化。"""
    try:
        import lattice_kappa as lk
        from ase.io import read as ase_read
    except Exception as e:
        sys.exit("[ERROR] shengbte 需要 lattice_kappa/ase：%s" % e)

    need = ("FORCE_CONSTANTS_2ND", "FORCE_CONSTANTS_3RD")
    missing = [f for f in need if not (sb_src / f).is_file()]
    if missing:
        sys.exit("[ERROR] SOLVER=shengbte 但 %s 缺 %s。\n"
                 "        请确认 S5_fc 的 step.conf 里 EXPORT_SHENGBTE=true（且 hiphive 可用），"
                 "重跑 S5_fc 后再来。" % (sb_src, ", ".join(missing)))
    for f in need:
        shutil.copyfile(sb_src / f, out / f)
    # POSCAR：优先用 shengbte 导出时的原胞（与力常数同源）
    src_pos = sb_src / "POSCAR" if (sb_src / "POSCAR").is_file() else (out / "POSCAR")
    atoms = ase_read(str(src_pos), format="vasp")
    shutil.copyfile(src_pos, out / "POSCAR")

    params = kc.read_kl_params(out / kc.KL_PARAMS)
    sc = [int(x) for x in (params.get("SUPERCELL") or "2 2 2").split()]
    C = {"kappa_mesh": [int(x) for x in mesh.split()],
         "kappa_t_min": conf["T_MIN"], "kappa_t_max": conf["T_MAX"],
         "kappa_t_step": conf["T_STEP"], "kappa_scalebroad": conf["SCALEBROAD"],
         "kappa_isotope": conf["ISOTOPE"],
         "kappa_convergence": bool(conf["KAPPA_CONVERGENCE"])}
    lk._write_shengbte_control(C, atoms, sc, out / "CONTROL", use_nac)
    print("[OK] ShengBTE 输入就绪：FORCE_CONSTANTS_2ND/3RD（拷自 S5）+ CONTROL")
    if use_nac:
        print("[WARN] CONTROL 已置 nonanalytic=T，但未自动写 born/epsilon；"
              "极性材料请手动在 CONTROL 补 Born 有效电荷与介电张量。")


def prepare_fourphonon(cwd, out, sb_src, conf, mesh, use_nac, ngpu):
    """fourphonon 输入准备 = ShengBTE 同源（FORCE_CONSTANTS_2ND/3RD 格式一致），
    复用 prepare_shengbte 后把 CONTROL 的 &flags 改成 fourphonon 3ph RTA：
      - convergence=F（RTA；fourphonon 迭代解未移植 GPU、且大胞易发散）
      - four_phonon=F（3ph only；4ph 需 fc4，S5 不产 fc4，勿开）
    运行时用 multi-GPU 提交模板（rank=GPU，acc_set_device_num 自动分卡）。"""
    prepare_shengbte(cwd, out, sb_src, conf, mesh, use_nac)
    ctl = out / "CONTROL"
    txt = ctl.read_text(encoding="utf-8")
    # 确保 flags 块含 four_phonon=F；convergence 强制 F（RTA）
    if "four_phonon" not in txt:
        txt = txt.replace("nanowires=F,", "nanowires=F,\n  four_phonon=F,", 1)
    txt = txt.replace("convergence=T,", "convergence=F,", 1)
    txt = txt.replace("convergence=F,", "convergence=F,", 1)
    ctl.write_text(txt, encoding="utf-8")
    # 若 kappa_convergence 默认 True 生成的是 CONV，这里显式对齐 RTA
    print("[OK] fourphonon(%d-GPU) 输入就绪：ShengBTE 格式 fc2/fc3 + CONTROL"
          "(four_phonon=F, convergence=F RTA)" % ngpu)


def main():
    cwd = Path.cwd()
    out = cwd / OUTDIR
    out.mkdir(exist_ok=True)
    conf = stepconf.load(SPEC, STEP)
    fcd = cwd / FC_DIR                 # step5_fc
    p3d = fcd / "phono3py"            # S5 产出的 phono3py 格式子目录
    sbd = fcd / "shengbte"           # S5 产出的 shengbte 力常数子目录
    if not fcd.is_dir():
        sys.exit("[ERROR] 找不到 step5_fc")
    if not p3d.is_dir():
        sys.exit("[ERROR] 找不到 step5_fc/phono3py（S5 拟合未完成或为旧版布局）")

    for f in ("fc2.hdf5", "fc3.hdf5", "phono3py_disp.yaml"):
        if not (p3d / f).is_file():
            sys.exit("[ERROR] %s 缺 %s（step5 力常数没建成）" % (p3d, f))
        shutil.copyfile(p3d / f, out / f)
    if (p3d / "BORN").is_file():
        shutil.copyfile(p3d / "BORN", out / "BORN")
    for f in ("POSCAR", kc.KL_PARAMS):
        src = p3d / f if (p3d / f).is_file() else fcd / f
        if src.is_file():
            shutil.copyfile(src, out / f)
    use_nac = (out / "BORN").is_file()

    # 维度 + 2D NAC 门槛：phono3py 只有 3D-Wang/Gonze 方案，对真 2D 是近似
    #   （2D 极性材料 LO-TO 在 q->0 应趋零，3D 方案给的是随真空变化的伪劈裂）。
    #   auto（默认）：2D 不用 NAC，3D 随 BORN；on/off 强制。正确的 2D-NAC 在 QE。
    params0 = kc.read_kl_params(out / kc.KL_PARAMS)
    dim = (params0.get("DIM") or "").lower()
    try:
        _d, vac_axis = kc.resolve_dim(out / "POSCAR", dim or "auto")
        dim = dim or _d
    except Exception:
        vac_axis = None
    nac_mode = str(conf["KAPPA_NAC"] or "auto").strip().lower()
    if nac_mode in ("on", "true", "1", "yes"):
        pass
    elif nac_mode in ("off", "false", "0", "no"):
        use_nac = False
    elif dim == "2d" and use_nac:
        use_nac = False
        print("[WARN] 2D + KAPPA_NAC=auto：默认不对 2D 施加 phono3py 的 3D-NAC")
        print("       （LO-TO 在 2D 应趋零，3D 方案是随真空变化的伪劈裂）。")
        print("       要强制用请设 KAPPA_NAC=on；正确的 2D-NAC 需用 QE 的 2D-DFPT。")

    # 稳定性闸：step5 判过虚频才该到这（phonon_summary.json 在 step5_fc 根）
    ps = fcd / "phonon_summary.json"
    if ps.is_file():
        import json
        try:
            if not json.loads(ps.read_text()).get("stable", True):
                sys.exit("[ERROR] step5 判定声子谱有虚频（不稳定），热导率无物理意义，已中止。")
        except Exception:
            pass

    params = kc.read_kl_params(out / kc.KL_PARAMS)
    solver = str(conf["SOLVER"]).lower()
    vac_ax = vac_axis if vac_axis is not None else 2
    # MESH_OVERRIDE / kl_params 的 MESH 也过 mesh_str：3D 被误写成 "N N 1" 时自动纠正。
    mesh_list = meshes(conf, params, dim, vac_ax)
    # ★ phono3py --mesh 是【原胞基矢轴序】，而 meshes() 产出的是 POSCAR/Cartesian 轴序。
    #   primitive_matrix 重排基矢后必须重排网格，否则 phono3py 报 Grid symmetry is broken
    #   （2026-09-21 MoS₂ S6 实测）。见 _mesh_prim_remap。
    if dim == "2d":
        try:
            _pl = kc.read_poscar_cell_frac(out / "POSCAR")[0]
            _prl = _prim_lattice(out)
            if _pl and _prl:
                _before = list(mesh_list)
                mesh_list = [_mesh_prim_remap(m, _pl, _prl, vac_ax) for m in mesh_list]
                if mesh_list != _before:
                    print("[..] q 网格按原胞基矢重排（POSCAR 序 -> 原胞序）：%s -> %s"
                          % (_before, mesh_list))
        except Exception as _e:                    # noqa: BLE001
            print("[WARN] q 网格原胞重排失败（%s），按原样使用" % _e)
    mesh = mesh_list[-1]
    if solver != "phono3py" and len(mesh_list) > 1:
        sys.exit("[ERROR] MESH_SCAN 多套网格只有 phono3py 支持（shengbte/fourphonon 单套）；"
                 "当前 SOLVER=%s，请把 MESH_SCAN 留空或只写一套。" % solver)
    if dim == "2d" and len(mesh_list) == 1:
        print("[..] 2D 提示：MESH_SCAN=auto 可一次跑三档网格做收敛判据（P1-1，5%% 判据）")
    bte_primary = str(conf["BTE_METHOD"] or "rta").lower()
    # RTA vs 完整解（P1-2）：auto = 2D 打开（正规过程主导，RTA 会低估 κ）
    _cm = str(conf["COMPARE_LBTE"] or "auto").strip().lower()
    if _cm in ("on", "true", "1", "yes"):
        compare_lbte = True
    elif _cm in ("off", "false", "0", "no"):
        compare_lbte = False
    else:
        # auto：不自动跑 LBTE 对照，只提示。
        #   实测（2026-09-16，Mo2S3 2D，32×32×1）：RTA 三档网格合计 ~5 分钟，
        #   而 LBTE 的碰撞矩阵是 (nq×nb×3)² 稠密矩阵、逐温度 dsyev 对角化，
        #   同一网格上跑了 >50 分钟仍未结束（量级 1~2 个数量级）。默认静默开这种
        #   开销不可接受，所以 auto = 关；文档的"RTA 会低估"结论改成显式提示。
        compare_lbte = False
        if dim == "2d":
            print("[WARN] 2D 的正规过程可能让 RTA 低估 κ（文档结论），但本步默认不做 LBTE "
                  "对照：实测 LBTE 比 RTA 贵 1~2 个数量级（稠密碰撞矩阵对角化）。\n"
                  "       需要对照就在 step.conf 里设 COMPARE_LBTE = on（只跑最粗那档网格 + 仅 300 K）。")
    plan = [(m, bte_primary) for m in mesh_list]
    ts_override = None
    if compare_lbte and solver == "phono3py":
        # 对照跑在【最粗那档】网格上、且【只跑 300 K】（wangchao 2026-09-21）：
        #   LBTE 的碰撞矩阵按 (nq·nb)² 稠密、逐温度对角化，比同网格 RTA 贵 1~2 个数量级；
        #   "RTA 是否低估"只需一个温度点 + 最省的那档网格就够。
        #   （同一档网格既跑 RTA 又跑 LBTE → 输出文件名会撞，tagged 改名，见 build_phono3py_cmd）
        _other = "lbte" if bte_primary == "rta" else "rta"
        plan.append((mesh_list[0], _other))
        ts_override = {_other: "300"}
    tagged = len({m for m, _ in plan}) != len(plan)
    primary = "%s|%s" % (mesh_list[-1], bte_primary)
    # P1-1/P1-2 的代价必须显式可见：一次 κ 求解很贵，"默认开了个对照"不该是惊喜。
    _extra = []
    if len(mesh_list) > 1:
        _extra.append("网格扫描 %d 档" % len(mesh_list))
    if len(plan) > len(mesh_list):
        _extra.append("%s 对照" % plan[-1][1].upper())
    if _extra:
        print("[..] 本步共 %d 次完整 κ 求解（%s）—— 每次都是一整套 phono3py 计算，"
              "不想跑就在 step.conf 里关：MESH_SCAN = / COMPARE_LBTE = off"
              % (len(plan), " + ".join(_extra)))
    if compare_lbte:
        print("[WARN] 已开启 LBTE 对照：该次比同网格 RTA 慢 1~2 个数量级"
              "（碰撞矩阵按 (nq×nb×3)² 稠密对角化），请留足墙钟时间。")
    ts = " ".join(str(t) for t in range(conf["T_MIN"], conf["T_MAX"] + 1, conf["T_STEP"]))
    print("[..] 求解器=%s mesh=%s 温度=%s K NAC=%s DIM=%s"
          % (solver, " -> ".join(mesh_list), ts, use_nac, dim or "?"))

    # 2D κ 厚度归一化因子（3D 时 factor=1、不归一）
    # REQ_2D：2D 时提交模板渲染结果里必须出现的字面量。项目级副本遮蔽技能模板时，
    #   老副本没有这些字段 → 归一化被静默吃掉（实测 2026-09-16：本项目 step6_kappa/
    #   submit_shengbte.tpl 就是这种 0 钩子的旧拷贝）。这里硬拦。
    REQ_2D = ("kappa_2d_normalized",) if dim == "2d" else ()
    factor, meta = 1.0, {"dim": dim or "?"}
    thick2d = None
    if dim == "2d":
        try:
            factor, m2 = two_d_norm_factor(out / "POSCAR", vac_axis, conf["KAPPA_2D_THICKNESS"])
            meta.update(m2)
            thick2d = dict(m2)
            thick2d.update({"dim": "2d",
                            "vac_axis": int(vac_axis if vac_axis is not None else 2),
                            "source": "kl-dft-cpu:step6_kappa",
                            "poscar": "POSCAR"})
            print("[..] 2D κ 归一化：h⊥=%.3f Å（|c|=%.3f Å）层厚 d=%.3f Å "
                  "（span=%.3f，真空隙=%.3f）→ factor=h⊥/d=%.4f  [%s]"
                  % (m2["h_perp_A"], m2["Lz_A"], m2["thickness_d_A"], m2["atomic_span_A"],
                     m2["vacuum_gap_A"], factor, m2["thickness_convention"]))
            if thick2d["vacuum_gap_A"] <= thick2d["atomic_span_A"]:
                print("[WARN] 真空隙 %.2f Å ≤ 层厚 %.2f Å：Wigner-Seitz 最近镜像可能跨真空，"
                      "Born-Huang 约束里的 r_ij 会取错；建议回 S1 加厚真空重跑。"
                      % (thick2d["vacuum_gap_A"], thick2d["atomic_span_A"]))
        except Exception as e:
            print("[WARN] 2D 归一化因子算失败，只出原始 κ —— 面内 κ 仍是含真空胞的体积"
                  "口径（被 h⊥ 稀释，数值不可直接使用）：%s" % e)

    here = Path(__file__).resolve().parent
    if solver == "phono3py":
        cmd = build_phono3py_cmd(
            plan, ts, conf["ISOTOPE"], use_nac,
            build_extract(factor, meta, thick2d, extract_plan(plan, tagged), primary),
            tagged=tagged, ts_override=ts_override)
        print("[..] BTE 方法=%s%s" % (bte_primary,
              ("；另跑 %d 次（%s 对照，P1-2）" % (len(plan) - len(mesh_list),
               "lbte" if bte_primary == "rta" else "rta"))
              if len(plan) > len(mesh_list) else ""))
        tpl = kc.resolve_submit(here, "3d", "submit_p3py")   # 单节点，无 2D/3D 之分
        # submit_p3py.tpl 的 {{CONDA_SH}}/{{CONDA_ENV}} 必须补传，否则残留字面占位符
        # （运行时 source {{CONDA_SH}} 报 No such file）。kl-dft 的 phono3py 环境是
        # atomate2_p_a（见 skill.yaml 的 conda 字段），与 MACE 技能的 mace_cpu 不同，
        # 故 CONDA_ENV 固定 atomate2_p_a、CONDA_SH 用集群 conda.sh（与 S5_fc 模板一致）。
        # 并行布局（2026-09-18 标定）：本环境 phono3py 是 OpenMP-only（C 扩展只链
        # libgomp，无 mpi4py/MPI），单节点只能是「1 进程 × N 个 OMP 线程」：
        #   --ntasks=1 --cpus-per-task=N，OMP_NUM_THREADS=N（见 submit_p3py.tpl 文件头）。
        _p3_threads = int(conf["P3PY_OMP_THREADS"] or 96)
        print("[..] phono3py 布局：1 进程 x %d OMP 线程（OpenMP-only，不能用 mpirun）"
              % _p3_threads)
        kc.write_submit(tpl, out / "submit.sh",
                        {"JOBNAME": kc.new_jobname(cwd, "S6kappa"),
                         "CONDA_SH": (conf["CONDA_SH"]
                                      or "/public/home/wangchao/miniconda3/etc/profile.d/conda.sh"),
                         "CONDA_ENV": "atomate2_p_a",
                         "CPUS_PER_TASK": str(_p3_threads),
                         "QOS": str(conf["SBATCH_QOS"] or "premium"),
                         "P3PY_CMD": cmd},
                        require=REQ_2D + ("--ntasks=1",),
                        label="phono3py 提交模板：")
    elif solver == "shengbte":
        prepare_shengbte(cwd, out, sbd, conf, mesh, use_nac)
        tpl = kc.resolve_submit(here, "3d", "submit_shengbte")
        # MPI/OMP 布局：一个 rank 一个 NUMA 域。rank=1 会退化成串行（ShengBTE 靠
        # MPI 按 q 点并行），所以这里必须显式算出来，不能沿用模板里的常量。
        _cpus_per_numa = int(conf["SHENGBTE_CORES_PER_NUMA"] or 24)
        _total = int(conf["SHENGBTE_TOTAL_CORES"] or 96)
        _nt_raw = str(conf["SHENGBTE_NTASKS"] or "auto").strip()
        if _nt_raw and _nt_raw.lower() != "auto":
            try:
                _ntasks = int(_nt_raw)
            except ValueError:
                sys.exit("[ERROR] SHENGBTE_NTASKS=%r 不是整数也不是 auto" % _nt_raw)
        else:
            _ntasks = max(1, _total // _cpus_per_numa)
        if _ntasks <= 1:
            sys.exit("[ERROR] SHENGBTE_NTASKS 算出来是 %d —— ShengBTE 靠 MPI 按 q 点"
                     "并行，单进程等于串行（实测 11.8 h 零产物）。请检查 step.conf 的"
                     " SHENGBTE_TOTAL_CORES=%d / SHENGBTE_CORES_PER_NUMA=%d。"
                     % (_ntasks, _total, _cpus_per_numa))
        if _ntasks * _cpus_per_numa != _total:
            print("[WARN] SHENGBTE 布局 %d rank x %d 线程 = %d 核，与 "
                  "SHENGBTE_TOTAL_CORES=%d 不一致（按 rank 数为准）"
                  % (_ntasks, _cpus_per_numa, _ntasks * _cpus_per_numa, _total))
        print("[..] ShengBTE 布局：%d MPI ranks x %d OMP threads = %d 核"
              % (_ntasks, _cpus_per_numa, _ntasks * _cpus_per_numa))
        # NTASKS/CPUS_PER_TASK/QOS：2026-09-16 从 taskflow-v2.0 同步。
        #   ShengBTE 靠 MPI 按 q 点并行，模板里写死 mpirun -n 1 = 串行（实测 11.8 h
        #   零产物），所以这三个必须由 gen 按集群 NUMA 拓扑算好填进去；模板里另有
        #   SLURM_NTASKS<=1 的退化保护兜底。
        kc.write_submit(tpl, out / "submit.sh",
                        {"JOBNAME": kc.new_jobname(cwd, "S6kappa"),
                         "SHENGBTE_EXE": conf["SHENGBTE_EXE"],
                         "NTASKS": str(_ntasks),
                         "CPUS_PER_TASK": str(_cpus_per_numa),
                         "QOS": str(conf["SBATCH_QOS"] or "premium"),
                         **norm_subs(factor, meta, thick2d)},
                        require=REQ_2D + ("--map-by numa",), label="2D 归一化：")
        # [FIX P48] SLURM QoS comes from step.conf, not from whichever copy of the
        # template already sits on the cluster.  autozt never overwrites an existing
        # gen_need asset ("材料目录已有的文件不覆盖"), and every project keeps its own
        # project_setting/templates/step6_kappa/submit_shengbte.tpl from the day it
        # was initialised -- so a project copy carrying --qos=premium pins S6_kappa
        # to premium's 5-jobs-per-user limit forever, no matter what the skill's
        # template says (measured 2026-09-16: Mg4C60 kappa jobs stuck at
        # QOSMaxJobsPerUserLimit for 40+ min while the local template said
        # {{QOS}}).  Rewrite the RENDERED submit.sh from conf["SBATCH_QOS"], which
        # is always fresh because the gen script itself is always re-pushed.
        #
        # 同一段里重写"完成标记"那一行（项目级旧模板会遮蔽技能模板）。
        #
        # [FIX P49 2026-09-18，移植自 taskflow-v2.0] 原先这里把
        # BTE.KappaTensorVsT_sg 也算作"完成"，理由是有条注释以为 _sg 是 ShengBTE 的
        # RTA 名字。**物理上反了**：ShengBTE 源码（Src/ShengBTE.f90）里 _sg 由
        # kappasg(energy, velocity) 写出，程序自己打印 "kappa in the small-grain
        # limit" —— 那是纯谐性、完全不含三声子散射的量（最小/边界限 κ），而且在
        # calculate_Vp 开始之前就写完了。真正的 RTA 产物是 BTE.KappaTensorVsT_RTA
        # （unit 303，温度循环里 calculate_Vp 之后写出）；BTE.KappaTensorVsT_CONV
        # （unit 403）是迭代自洽解。Mg4C60 7x7x7 实测（2026-09-17）100 K：
        # _sg = 0.0782 vs _RTA = 0.4740 W/m/K，**差 6.06 倍**；而崩掉的作业留下的
        # 恰好是完整的 _sg（8 个温度 15 分钟写完）与只有一行的 _RTA（13.7 h 后）。
        # 因此接受 _sg 等于把"三声子段根本没跑完"判成绿色的 KAPPA_DONE ——
        # 交付的是最小 κ 而不是热导率。这里只认 _CONV/_RTA；_sg 的存在改由模板写成
        # kappa_summary.json 里的诊断字段（small_grain_limit_present / _rows / warning）。
        _qos = str(conf["SBATCH_QOS"] or "premium").strip()
        _sh = out / "submit.sh"
        _txt = _sh.read_text(encoding="utf-8")
        _txt = re.sub(r"(?m)^#SBATCH --qos=.*$", "#SBATCH --qos=%s" % _qos, _txt)
        _cand = 'cand = ["BTE.KappaTensorVsT_CONV", "BTE.KappaTensorVsT_RTA"]'
        _txt, _ncand = re.subn(r"(?m)^cand = \[.*\]$", _cand, _txt)
        if _ncand == 0:
            print("[WARN] submit.sh 没有 cand = [...] 行（旧模板？）：kappa_summary 的"
                  "完成标记可能仍只认 _CONV/_RTA，RTA-only 运行会被误判 KAPPA_DONE=false")
        _sh.write_text(_txt, encoding="utf-8", newline="\n")
        print("[..] SLURM QoS = %s（step.conf 的 SBATCH_QOS；模板里的旧值会被覆盖）"
              % _qos)
    elif solver == "fourphonon":
        # StepConf 没有 .get()（只有 __getitem__）—— 原来写 conf.get(...) 会直接
        # AttributeError，SOLVER=fourphonon 根本进不来。用 [] 取（缺键回落 SPEC 默认）。
        fp_exe = str(conf["FOURPHONON_EXE"] or "").strip()
        if not fp_exe:
            sys.exit("[ERROR] SOLVER=fourphonon 但 step.conf 没填 FOURPHONON_EXE"
                     "（multi-GPU 版绝对路径）。见 README「fourphonon」节。")
        ngpu = int(conf["FOURPHONON_NGPU"] or 4)

        if not bool(conf["ALLOW_FOURPHONON_GPU_PHONOPY"]):
            sys.exit("[ERROR] SOLVER=fourphonon 已停用：FourPhonon v1.3 **GPU** 版对"
                     " phonopy/ShengBTE 格式 fc2 有官方 OpenACC 数值 bug（κ 错 67~71×，"
                     "只有 espresso 路径正确、官方未修），本技能喂的正是该格式，结果不可信。\n"
                     "         确知风险仍要跑（做对照 / 复现 bug）：把 step.conf 的"
                     " ALLOW_FOURPHONON_GPU_PHONOPY 设为 true。见 README「fourphonon」节。")
        print("[WARN] SOLVER=fourphonon：GPU 版对 phonopy fc2 有已知数值 bug（κ 错 67~71×），"
              "本次结果只能用于对照 / 复现，不可当真值。")
        prepare_fourphonon(cwd, out, sbd, conf, mesh, use_nac, ngpu)
        tpl = kc.resolve_submit(here, "3d", "submit_fourphonon")
        kc.write_submit(tpl, out / "submit.sh",
                        {"JOBNAME": kc.new_jobname(cwd, "S6kappa"),
                         "FOURPHONON_EXE": fp_exe,
                         "FOURPHONON_NGPU": str(ngpu),
                         "FOURPHONON_CPUS_PER_GPU": str(conf["FOURPHONON_CPUS_PER_GPU"] or 8),
                         **norm_subs(factor, meta, thick2d)},
                        require=REQ_2D + ("--map-by numa",), label="2D 归一化：")
    else:
        sys.exit("[ERROR] SOLVER 只允许 phono3py / shengbte / fourphonon")
    stepconf.apply_submit(out / "submit.sh", conf.submit)
    print("[DONE] %s：submit.sh 就绪，提交后计算节点出 κ，写 kappa_summary.json(KAPPA_DONE)"
          % OUTDIR)


if __name__ == "__main__":
    main()