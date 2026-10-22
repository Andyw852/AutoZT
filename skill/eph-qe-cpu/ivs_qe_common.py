#!/usr/bin/env python3
"""eph-qe-cpu 可选组 ivs_dfpt（[0059]）的生成器：QE DFPT 直接算谷间零阶形变势，给 ke-dft-cpu 的 S8/S8.4 当 IVS 通道。

  S10_ivsq    gen_step10_ivs_qlist.py  （run: gen）谷表 -> 谷对 -> 等价类 -> q 作业            -> step10_ivs_qlist/ivs_qlist.json
  S11_ivsscf  gen_step11_ivs_scf.py    （作业）    环境探针 + pw.x SCF + 各谷能差自检          -> step11_ivs_scf/ivs_scf.ok
  S12_ivsph   gen_step12_ivs_ph.py     （作业，扇出 q-*）ph.x 单 q DFPT -> nscf -> ph.x AHC     -> step12_ivs_ph/q-*/ivs_ph.ok
  S13_ivs     gen_step13_ivs_read.py   （run: gen）D0、ħω -> dp_optical 星通道                -> step13_ivs_read/ivs_dfpt_params.json

数学与 QE 接口在公共池 skill/_common/ivs/ivs_dfpt.py（ke-dft-cpu 读结果也用它）。gen 在 <材料>/eph-qe-cpu/ 下运行。
输入（autozt 按 skill.yaml 的 push_paths 从本地 result 推过来，本地没有就不推、远端旧的删掉）：
  ivs_inputs/CONTCAR            ke-dft-cpu S1 的弛豫结构（IVS_STRUCTURE=auto 优先用它；没有就用材料 POSCAR 并告警）
  ivs_inputs/dpt_result.json    ke-dft-cpu S8.2 的多谷谷表（2D 默认读它；3D 没有，用 IVS_VALLEYS_E/H 手填）
  ivs_inputs/workflow_method.txt ke-dft-cpu S1 的 FUNC（IVS_INPUT_DFT=auto 据此设 QE 的 input_dft）
集群环境：setting/<集群>.yaml 的 qe_bin / qe_pseudo_dir / qe_setup -> step.conf 的 QE_BIN / QE_PSEUDO_DIR / QE_SETUP（autozt 注入）。
"""
import hashlib
import json
import shlex
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _d in (HERE, Path.cwd()):
    if str(_d) not in sys.path:
        sys.path.insert(0, str(_d))

STEP = {"qlist": "step10_ivs_qlist", "scf": "step11_ivs_scf", "ph": "step12_ivs_ph", "read": "step13_ivs_read"}
INPUTS = Path("ivs_inputs")
FANOUT_PREFIX = "q-"

# ---- step.conf [params]（技能默认在 templates/step.conf；项目级覆盖）----
PARAMS = {
    "IVS_STRUCTURE": ("auto", "str"),        # auto | ke | poscar
    "IVS_DIM": ("auto", "str"),              # auto | 2d | 3d
    "IVS_CARRIERS": ("auto", "str"),         # auto | electron | hole | both
    "IVS_VALLEYS_E": ("", "str"),            # 手填谷表（同 DP_OPTICAL_VALLEYS_E 写法；每星写一个成员即可，按对称性补全）
    "IVS_VALLEYS_H": ("", "str"),
    "IVS_PSEUDO_FILES": ("auto", "str"),     # Element:file；auto = 沿用 PSEUDO_FILES（再缺省 Element.upf）
    "PSEUDO_FILES": ("auto", "str"),
    "ATOMIC_MASSES": ("auto", "str"),
    "IVS_INPUT_DFT": ("auto", "str"),        # auto（按 ke S1 的 FUNC）| none | QE 的 input_dft 名
    "IVS_ECUTWFC": (60.0, "float"),
    "IVS_ECUTRHO": (0.0, "float"),           # 0 = 4 × IVS_ECUTWFC（模守恒）
    "IVS_KSPACING": (0.2, "float"),          # Å⁻¹（VASP KSPACING 口径，含 2π）
    "IVS_SCF_KGRID": ("auto", "str"),
    "IVS_SCF_CONV_THR": (1.0e-12, "float"),
    "IVS_TR2_PH": (1.0e-14, "float"),        # [0068] 原 1e-16；D₀/ħω 到 1e-14 已远低于模型误差，与 Perturbo 路线一致
    "IVS_AHC_NB": (2, "int"),                # AHC 窗口：带边上下各几条带
    "IVS_NBND_EXTRA": (6, "int"),
    "IVS_VALLEY_TOL_MEV": (30.0, "float"),
    "IVS_VALLEY_CHECK": ("error", "str"),    # error | warn
    "IVS_SYMPREC": (0.01, "float"),
    "IVS_NP": (4, "int"),                    # 每个作业的 MPI 进程数
    "IVS_SBATCH": ("", "str"),               # 额外 #SBATCH 选项，分号分隔
    "IVS_MPIRUN": ("auto", "str"),           # auto | none | 自定义（{np} 处填进程数）
    "IVS_KEEP_TMP": (False, "bool"),
    # [0059-r3] ODP：q 列表加一个 q=Γ 作业（q-odp），AHC 顺带给每个谷的谷内光学形变势；单谷载流子也算。
    #   auto/on = 算，off = 不算（只算谷间 IVS）
    "IVS_ODP": ("auto", "str"),     # auto | on | off
    # 极性（红外活性）模：exclude（出厂，与冻声子 ODP_POLAR 同口径：长程 Fröhlich 由 POP 管）= q-odp 的 DFPT 加 epsil 算
    #   Born 有效电荷判极性、不进谷内通道；include = 极性模的短程形变势也算进 ODP（不算 epsil，金属/零带隙也能跑）
    "IVS_ODP_POLAR": ("exclude", "str"),
    "IVS_ODP_POLAR_TOL": (0.1, "float"),   # 模有效电荷（e）大于它算极性（= ke-dft-cpu 的 ODP_POLAR_TOL）
    "QE_BIN": ("", "str"),                   # autozt 从 setting/<集群>.yaml 的 qe_bin 注入
    "QE_PSEUDO_DIR": ("", "str"),            # qe_pseudo_dir
    "QE_SETUP": ("", "str"),                 # qe_setup：作业开头执行的 shell（module load / source …）；jzzn = 用 eph-qe-cpu 的 jzzn 环境
}


def load_params():
    if not Path("step.conf").is_file():
        return {k: v for k, (v, _t) in PARAMS.items()}
    import stepconf
    return dict(stepconf.load(PARAMS, cwd=".", strict=False).params)


def _iv():
    import ivs_dfpt
    return ivs_dfpt


def _sha_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _write(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def _read_json(path, what):
    p = Path(path)
    if not p.is_file():
        raise SystemExit("[ERROR] 缺 %s（%s）" % (p, what))
    return json.loads(p.read_text(encoding="utf-8"))


def _elem_map(text, name, cast=str):
    """'Mo:Mo.upf S:S.upf' -> {Mo: Mo.upf, S: S.upf}；auto/空 -> {}。"""
    t = str(text or "").strip()
    if not t or t.lower() == "auto":
        return {}
    out = {}
    for tok in t.replace(",", " ").replace(";", " ").split():
        if ":" not in tok:
            raise SystemExit("[ERROR] %s 要写成 元素:值，收到 %r" % (name, tok))
        k, v = tok.split(":", 1)
        k = k.strip()
        k = k[:1].upper() + k[1:].lower()
        out[k] = cast(v.strip())
    return out


# ---------------------------------------------------------------------------
# 作业环境与提交脚本
# ---------------------------------------------------------------------------
def env_block(P):
    """作业开头的环境：QE_SETUP（或 jzzn 预设）+ QE_BIN 进 PATH + run_mpi + IVS_PSEUDO_DIR。"""
    setup = str(P.get("QE_SETUP") or "").strip()
    lines = ["# ---- QE 环境（[0059] ivs_dfpt）----", "export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}"]
    np_ = int(P.get("IVS_NP") or 1)
    if setup.lower() == "jzzn":
        import eph_common
        lines.append(eph_common.JZZN_ENVIRONMENT.rstrip())
    else:
        if setup:
            lines.append(setup)
        mp = str(P.get("IVS_MPIRUN") or "auto").strip()
        if mp.lower() == "none" or np_ <= 1:
            lines.append('run_mpi() { "$@"; }')
        elif mp.lower() == "auto":
            lines.append('run_mpi() { if command -v mpirun >/dev/null 2>&1; then mpirun -np %d "$@"; '
                         'else echo "[WARN] 没有 mpirun，串行运行" >&2; "$@"; fi; }' % np_)
        else:
            lines.append('run_mpi() { %s "$@"; }' % mp.replace("{np}", str(np_)))
    qb = str(P.get("QE_BIN") or "").strip()
    if qb:
        lines.append('export PATH=%s:"$PATH"' % shlex.quote(qb))
    pd = str(P.get("QE_PSEUDO_DIR") or "").strip()
    lines.append('IVS_PSEUDO_DIR=%s' % (shlex.quote(pd) if pd else '"${PSEUDO_DIR:-}"'))
    lines.append('[ -n "$IVS_PSEUDO_DIR" ] || { echo "[ERROR] 没有赝势目录：setting/<集群>.yaml 写 qe_pseudo_dir，'
                 '或 step.conf 写 QE_PSEUDO_DIR" >&2; exit 2; }')
    return "\n".join(lines) + "\n"


def render_submit(dest, job, command, P):
    tpl = Path("submit_ivs_qe.tpl")
    if not tpl.is_file():
        tpl = HERE / "submit_ivs_qe.tpl"
    if not tpl.is_file():
        raise SystemExit("[ERROR] 缺 submit_ivs_qe.tpl（gen_need 里要有它）")
    extra = [x.strip() for x in str(P.get("IVS_SBATCH") or "").split(";") if x.strip()]
    sb = "\n".join("#SBATCH %s" % x for x in extra)
    txt = (tpl.read_text(encoding="utf-8").replace("{{JOBNAME}}", job).replace("{{NP}}", str(int(P["IVS_NP"])))
           .replace("{{SBATCH_EXTRA}}", sb).replace("{{COMMAND}}", env_block(P) + command))
    _write(Path(dest) / "submit.sh", txt)
    (Path(dest) / "submit.sh").chmod(0o755)


def _pseudo_links(types, pseudos):
    out = ["mkdir -p pseudo"]
    for t in types:
        f = shlex.quote(pseudos[t])
        out.append('ln -sfn "$IVS_PSEUDO_DIR/%s" "pseudo/%s"' % (f, f))
    return "\n".join(out)


def _copy_job_helper(dest):
    src = Path("ivs_qe_job.py") if Path("ivs_qe_job.py").is_file() else HERE / "ivs_qe_job.py"
    if not src.is_file():
        raise SystemExit("[ERROR] 缺 ivs_qe_job.py（gen_need 里要有它）")
    _write(Path(dest) / "ivs_qe_job.py", src.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# 结构、赝势、泛函
# ---------------------------------------------------------------------------
def pick_structure(P):
    iv = _iv()
    mode = str(P.get("IVS_STRUCTURE") or "auto").strip().lower()
    if mode not in ("auto", "ke", "poscar"):
        raise SystemExit("[ERROR] IVS_STRUCTURE=%r 无效，只允许 auto / ke / poscar" % mode)
    ke = INPUTS / "CONTCAR"
    if mode in ("auto", "ke") and ke.is_file():
        f, src = ke, "ke-dft-cpu S1 CONTCAR（ivs_inputs/CONTCAR）"
    elif mode == "ke":
        raise SystemExit("[ERROR] IVS_STRUCTURE=ke，但没有 ivs_inputs/CONTCAR：先跑完 ke-dft-cpu 的 S1_opt（本地 result 里要有"
                         " ke-dft-cpu/result/step1_opt/CONTCAR），再 retry 本步")
    else:
        f, src = Path("POSCAR"), "材料 POSCAR"
        if not f.is_file():
            raise SystemExit("[ERROR] 没有 ivs_inputs/CONTCAR，材料目录也没有 POSCAR")
        if mode == "auto":
            print("[WARN] 没有 ke-dft-cpu S1 的 CONTCAR（ivs_inputs/CONTCAR）：用材料 POSCAR —— 谷位置要与 AMSET 一致，"
                  "正式计算应先跑完 ke-dft-cpu S1（IVS_STRUCTURE=ke 强制）")
    try:
        st = iv.read_poscar_text(f.read_text(encoding="utf-8", errors="replace"))
    except ValueError as e:
        raise SystemExit("[ERROR] %s 解析失败：%s" % (f, e))
    st["source"] = src
    st["file"] = str(f)
    st["sha256"] = _sha_file(f)
    return st


def species_maps(P, types):
    import eph_common
    pse = _elem_map(P.get("IVS_PSEUDO_FILES"), "IVS_PSEUDO_FILES") or _elem_map(P.get("PSEUDO_FILES"), "PSEUDO_FILES")
    pse = {t: pse.get(t, "%s.upf" % t) for t in types}
    ms = _elem_map(P.get("ATOMIC_MASSES"), "ATOMIC_MASSES", float)
    masses = {}
    for t in types:
        if t in ms:
            masses[t] = ms[t]
        elif t in eph_common.PERIODIC_MASSES:
            masses[t] = eph_common.PERIODIC_MASSES[t]
        else:
            raise SystemExit("[ERROR] 不认识的元素 %r（ATOMIC_MASSES 里给它的质量）" % t)
    return pse, masses


FUNC_TO_QE = {"pbesol": "pbesol", "pbe": "pbe", "pbe-d3": "pbe"}


def input_dft(P):
    """IVS_INPUT_DFT -> (QE input_dft 或 None, 说明)。"""
    v = str(P.get("IVS_INPUT_DFT") or "auto").strip()
    if v.lower() in ("none", ""):
        return None, "赝势自带的泛函（IVS_INPUT_DFT=none）"
    if v.lower() != "auto":
        return v, "step.conf IVS_INPUT_DFT"
    f = INPUTS / "workflow_method.txt"
    if f.is_file():
        for ln in f.read_text(encoding="utf-8", errors="replace").splitlines():
            if ln.strip().upper().startswith("FUNC="):
                func = ln.split("=", 1)[1].strip().lower()
                if func in FUNC_TO_QE:
                    note = "ke-dft-cpu S1 FUNC=%s" % func
                    if func == "pbe-d3":
                        note += "（DFPT 不含 D3 色散：结构仍是 VASP PBE-D3 弛豫的，只差声子频率里的色散贡献）"
                    return FUNC_TO_QE[func], note
                print("[WARN] ke-dft-cpu S1 的 FUNC=%s 没有对应的 QE 泛函：用赝势自带的" % func)
    print("[WARN] IVS_INPUT_DFT=auto，但没有 ivs_inputs/workflow_method.txt：用赝势自带的泛函（应与 VASP 一致，见 README）")
    return None, "赝势自带的泛函（没有 ke S1 的 FUNC）"


def pw_input(calc, P, q1, prefix, outdir, kblock, nbnd=None, extra_sys=""):
    iv = _iv()
    st = q1["structure"]
    dft = q1["qe"]["input_dft"]
    ecut = float(P["IVS_ECUTWFC"])
    ecutrho = float(P["IVS_ECUTRHO"] or 0) or 4.0 * ecut
    sys_lines = ["  ibrav = 0", "  nat = %d" % len(st["symbols"]), "  ntyp = %d" % len(st["types"]),
                 "  ecutwfc = %.8g" % ecut, "  ecutrho = %.8g" % ecutrho, "  occupations = 'fixed'"]
    if dft:
        sys_lines.append("  input_dft = '%s'" % dft)
    if nbnd:
        sys_lines.append("  nbnd = %s" % (nbnd if isinstance(nbnd, str) else int(nbnd)))
    if extra_sys:
        sys_lines.append(extra_sys.rstrip())
    if calc == "nscf":
        # [0059-r4] AHC 的 nscf 必须原样保留给定 k 列表：pw.x 默认按对称性扩增 k，
        #   会让 ahc_gkk/etk/etq 的记录数变成输入 k 数的整数倍（MoS2 q=Gamma：2 -> 4 个 k），
        #   ivs_qe_job.py 的尺寸校验和 ivs_dfpt.py 的逐 k 映射都会错位（Si 单 k 不暴露）。
        sys_lines.append("  nosym = .true.")
        sys_lines.append("  noinv = .true.")
    el = ["  conv_thr = %.3e" % float(P["IVS_SCF_CONV_THR"]), "  mixing_beta = 0.4"]
    if calc != "scf":
        el.append("  diago_full_acc = .true.")
    return ("&CONTROL\n  calculation = '%s'\n  prefix = '%s'\n  outdir = '%s'\n  pseudo_dir = './pseudo'\n"
            "  verbosity = 'high'\n/\n&SYSTEM\n%s\n/\n&ELECTRONS\n%s\n/\n" % (calc, prefix, outdir, "\n".join(sys_lines),
                                                                         "\n".join(el))
            + iv.qe_structure(st["lattice"], st["frac"], st["symbols"], st["types"], q1["qe"]["masses"],
                              q1["qe"]["pseudos"])
            + kblock)


# ---------------------------------------------------------------------------
# S10：谷表 -> q 作业
# ---------------------------------------------------------------------------
ODP_RUN = FANOUT_PREFIX + "odp"     # [0059-r3] q=Γ 作业的子目录名（与 q-01、q-02… 同一个扇出）


def add_odp_gamma(carriers, classes, runs):
    """[0059-r3] 每个载流子的每个星加一个谷内类（kind="intra"，q=Γ，初态 = 星的第一个成员），全部放进一个 q=Γ 作业。
    同一个 k（直接带隙时电子、空穴都在 Γ）只算一次。返回加了几个类。"""
    iv = _iv()
    run = {"id": len(runs), "q": [0.0, 0.0, 0.0], "k": [], "name": ODP_RUN}
    n = 0
    for car in iv.CARRIERS:
        for X in carriers.get(car) or []:
            k0 = iv.wrap_half(X["k"][0])
            ik = next((i for i, kk in enumerate(run["k"]) if iv.same_q(kk, k0, 1e-6)), None)
            if ik is None:
                run["k"].append(k0)
                ik = len(run["k"]) - 1
            g = len(X["k"])
            classes.append({"id": len(classes), "carrier": car, "from": X["label"], "to": X["label"], "kind": "intra",
                            "k_init": k0, "k_fin": k0, "q": [0.0, 0.0, 0.0], "run": run["id"], "ik": ik, "n_pairs": 1,
                            "pairs": [{"from": X["label"], "to": X["label"], "kind": "intra", "j": 0,
                                       "g_to": g, "g_from": g}]})
            n += 1
    if n:
        runs.append(run)
    return n


def gen_qlist(P):
    iv = _iv()
    st = pick_structure(P)
    dim = str(P.get("IVS_DIM") or "auto").strip().lower()
    if dim == "auto":
        dim = iv.detect_dim(st["lattice"], st["frac"])
    if dim not in ("2d", "3d"):
        raise SystemExit("[ERROR] IVS_DIM=%r 无效，只允许 auto / 2d / 3d" % dim)
    full3d = dim == "3d"
    rots = iv.point_group(st["lattice"], st["frac"], st["symbols"], symprec=float(P["IVS_SYMPREC"]))
    rops = iv.recip_ops(rots)
    print("[..] 结构：%s（%s，%d 原子，点群 %d 个操作）" % (st["source"], dim.upper(), len(st["symbols"]), len(rots)))
    want = str(P.get("IVS_CARRIERS") or "auto").strip().lower()
    if want not in ("auto", "electron", "hole", "both"):
        raise SystemExit("[ERROR] IVS_CARRIERS=%r 无效，只允许 auto / electron / hole / both" % want)
    # [0059-r3] ODP：q=Γ 的谷内光学耦合。开着时单谷载流子（Si 空穴 Γ）也保留：它没有谷间散射，但有谷内 ODP
    odp_mode = str(P.get("IVS_ODP") or "auto").strip().lower()
    if odp_mode not in ("auto", "on", "off"):
        raise SystemExit("[ERROR] IVS_ODP=%r 无效，只允许 auto / on / off" % odp_mode)
    odp = odp_mode != "off"
    dpt, dpt_sha = None, None
    f = INPUTS / "dpt_result.json"
    if f.is_file():
        dpt = json.loads(f.read_text(encoding="utf-8"))
        dpt_sha = _sha_file(f)
    carriers, cinfo = {}, {}
    for car, key in (("electron", "IVS_VALLEYS_E"), ("hole", "IVS_VALLEYS_H")):
        if want in ("electron", "hole") and car != want:
            continue
        spec = str(P.get(key) or "").strip()
        notes = []
        if spec:
            try:
                vs = iv.parse_valleys(spec, full3d=full3d)
            except ValueError as e:
                raise SystemExit("[ERROR] %s：%s" % (key, e))
            raw, src, use_grid, mv = vs, "step.conf %s" % key, False, None
        elif dim == "2d" and dpt is not None:
            vs, raw, mv = iv.valleys_from_dpt(dpt, car)
            src = "S8.2 dpt_result.json（sha256 %s）" % (dpt_sha or "")[:12]
            use_grid = bool(vs) and all(v.get("k_grid") for v in vs)
        else:
            if want in (car, "both"):
                raise SystemExit("[ERROR] %s 没有谷表：%s" % (car, "3D 在 step.conf 手填 %s（每个星写一个成员即可）" % key
                                                       if dim == "3d" else "没有 ivs_inputs/dpt_result.json（先跑完 "
                                                       "ke-dft-cpu S8.2_dpt），或手填 %s" % key))
            continue
        if not vs:
            if want in (car, "both"):
                raise SystemExit("[ERROR] %s 的谷表是空的（%s）" % (car, src))
            continue
        try:
            vs, alias = iv.merge_coincident(vs, full3d=full3d)
        except ValueError as e:
            raise SystemExit("[ERROR] %s 谷表（%s）：%s" % (car, src, e))
        if alias:
            notes.append("同中心的星并成一个：%s" % ", ".join("%s→%s" % kv for kv in alias.items()))
        stars, bad = [], None
        for v in vs:
            ks = v["k_grid"] if use_grid else v["k_frac"]
            full, norb, added = iv.expand_star(ks, rops, full3d=full3d)
            if norb > 1:
                bad = ("%s 谷表（%s）的星 %s 的成员不在同一个对称轨道上（%d 个轨道）：一个星只能是一组等价谷"
                       % (car, src, v["label"], norb))
            elif added and not spec:
                bad = ("%s 谷表（%s）的星 %s 只有 %d 个成员，按晶体对称性应有 %d 个：谷表不完整（S8.4 的态密度按这份谷表的"
                       "成员数平均，补成员会对不上）。手填 IVS_VALLEYS_%s 覆盖"
                       % (car, src, v["label"], len(ks), len(full), "E" if car == "electron" else "H"))
            if bad:
                break
            if added:
                notes.append("%s：按对称性补全成员 %d -> %d" % (v["label"], len(ks), len(full)))
            stars.append({"label": v["label"], "k": [list(map(float, k)) for k in full],
                          "dE_ref_eV": v.get("dE_eV")})
        if bad:
            if want == "auto":                     # 自动挑载流子：这个载流子跳过（另一个照算），不让整步失败
                print("[WARN] %s：本载流子跳过" % bad)
                continue
            raise SystemExit("[ERROR] " + bad)
        nmem = sum(len(s["k"]) for s in stars)
        if nmem <= 1:
            if not odp:
                print("[..] %s 只有一个谷（%s）：没有谷间散射，跳过" % (car, ", ".join(s["label"] for s in stars)))
                continue
            notes.append("只有一个谷（%s）：没有谷间散射，只算谷内 ODP" % ", ".join(s["label"] for s in stars))
        carriers[car] = stars
        full_sig = iv.valley_sig([{"k_frac": s["k"]} for s in stars], use_grid=False, full3d=full3d)
        cinfo[car] = {"valleys_source": src, "use_grid": use_grid, "multi_valley_status": mv,
                      "valley_sig": iv.valley_sig(raw, use_grid=use_grid, full3d=full3d),
                      "valley_sig_full": full_sig, "notes": notes,
                      "valleys": [{"label": s["label"], "g": len(s["k"]), "k_frac": s["k"], "dE_ref_eV": s["dE_ref_eV"]}
                                  for s in stars]}
        for n_ in notes:
            print("[..] %s：%s" % (car, n_))
    if not carriers:
        raise SystemExit("[ERROR] 没有%s的载流子：不需要算（IVS_CARRIERS=%s，IVS_ODP=%s）"
                         % ("可算" if odp else "多谷", want, odp_mode))
    multi = {car: stars for car, stars in carriers.items() if sum(len(s["k"]) for s in stars) > 1}
    classes, runs = iv.plan_pairs(multi, rops)
    if odp:
        n_intra = add_odp_gamma(carriers, classes, runs)
        print("[OK] ODP（Γ 点谷内光学耦合）：%d 个星，q = (0,0,0)，作业 %s" % (n_intra, ODP_RUN))
    if not runs:
        raise SystemExit("[ERROR] 没有要算的 q：没有多谷载流子，IVS_ODP=off")
    pse, masses = species_maps(P, st["types"])
    dft, dft_note = input_dft(P)
    prefix = "".join("%s%d" % (t.lower(), c) for t, c in zip(st["types"], st["counts"]))
    kg = str(P.get("IVS_SCF_KGRID") or "auto").strip().lower()
    kgrid = iv.auto_kgrid(st["lattice"], float(P["IVS_KSPACING"]), dim) if kg == "auto" else \
        [int(x) for x in kg.replace("x", " ").split()]
    if len(kgrid) != 3 or min(kgrid) < 1:
        raise SystemExit("[ERROR] IVS_SCF_KGRID=%r 应为 auto 或三个正整数" % kg)
    out = {"schema": iv.QLIST_SCHEMA, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
           "structure": {"lattice": st["lattice"], "symbols": st["symbols"], "types": st["types"], "counts": st["counts"],
                         "frac": st["frac"], "source": st["source"], "sha256": st["sha256"]},
           "dim": dim, "n_ops": len(rots), "symprec_A": float(P["IVS_SYMPREC"]),
           "qe": {"prefix": prefix, "pseudos": pse, "masses": masses, "input_dft": dft, "input_dft_note": dft_note,
                  "ecutwfc_Ry": float(P["IVS_ECUTWFC"]),
                  "ecutrho_Ry": float(P["IVS_ECUTRHO"] or 0) or 4.0 * float(P["IVS_ECUTWFC"]),
                  "kgrid": kgrid, "tr2_ph": float(P["IVS_TR2_PH"])},
           "inputs": {"dpt_result_sha256": dpt_sha},
           "carriers": cinfo,
           "classes": [{k: c[k] for k in ("id", "carrier", "from", "to", "kind", "k_init", "k_fin", "q", "run", "ik",
                                          "n_pairs", "pairs")} for c in classes],
           "runs": [{"id": r["id"], "q": r["q"], "k": r["k"],
                      "name": r.get("name") or "%s%02d" % (FANOUT_PREFIX, r["id"] + 1)} for r in runs]}
    dest = Path(STEP["qlist"])
    _write(dest / "ivs_qlist.json", json.dumps(iv.jsonable(out), indent=1, ensure_ascii=False) + "\n")
    for car in carriers:
        nc = [c for c in classes if c["carrier"] == car]
        print("[OK] %s：%s；%d 个等价谷对" % (car, ", ".join("%s×%d" % (s["label"], len(s["k"])) for s in carriers[car]),
                                           len(nc)))
    for r in runs:
        print("[..]   %s  q = (%s)  初态 k %d 个" % (r.get("name") or "%s%02d" % (FANOUT_PREFIX, r["id"] + 1),
                                                    ", ".join("%.4f" % x for x in r["q"]), len(r["k"])))
    print("[DONE] %s/ivs_qlist.json：%d 个 q（每个一个 ph.x 作业，S12_ivsph 扇出）；QE 泛函：%s"
          % (dest, len(runs), dft or dft_note))
    return 0


# ---------------------------------------------------------------------------
# S11：SCF + 谷能差自检
# ---------------------------------------------------------------------------
def gen_scf(P):
    iv = _iv()
    q1 = _read_json(Path(STEP["qlist"]) / "ivs_qlist.json", "先跑 S10_ivsq")
    if q1.get("schema") != iv.QLIST_SCHEMA:
        raise SystemExit("[ERROR] ivs_qlist.json 格式 %r 不是 %r：重跑 S10_ivsq" % (q1.get("schema"), iv.QLIST_SCHEMA))
    pre = q1["qe"]["prefix"]
    dest = Path(STEP["scf"])
    kg = q1["qe"]["kgrid"]
    scf = pw_input("scf", P, q1, pre, "./tmp", "K_POINTS automatic\n  %d %d %d 0 0 0\n" % tuple(kg))
    pts, ref = [], {"points": [], "carriers": {}}
    for car, c in q1["carriers"].items():
        rows = []
        for v in c["valleys"]:
            idx = []
            for k in v["k_frac"]:
                for i, p in enumerate(pts):
                    if iv.same_q(p, k, 1e-6):
                        break
                else:
                    pts.append(k)
                    i = len(pts) - 1
                idx.append(i)
            rows.append({"label": v["label"], "points": idx, "dE_ref_eV": v.get("dE_ref_eV")})
        ref["carriers"][car] = rows
    ref["points"] = pts
    valley = pw_input("bands", P, q1, pre, "./tmp_valley", iv.kpoints_list(pts), nbnd="__NBND__")
    types = q1["structure"]["types"]
    pses = " ".join("%s:%s" % (t, q1["qe"]["pseudos"][t]) for t in types)
    cmd = "\n".join([
        "rm -f ivs_scf.ok",
        "python3 ivs_qe_job.py env --pseudo-dir \"$IVS_PSEUDO_DIR\" --pseudos %s" % shlex.quote(pses),
        _pseudo_links(types, q1["qe"]["pseudos"]),
        "run_mpi pw.x -in scf.in > scf.out",
        "grep -q 'JOB DONE' scf.out",
        "rm -rf tmp_valley && mkdir -p tmp_valley && cp -r tmp/%s.save tmp_valley/" % pre,
        "python3 ivs_qe_job.py valley --xml tmp/%s.save/data-file-schema.xml --scf-out scf.out "
        "--tpl valley.in.tpl --out valley.in --extra 4" % pre,
        "run_mpi pw.x -in valley.in > valley.out",
        "grep -q 'JOB DONE' valley.out",
        "python3 ivs_qe_job.py check --xml tmp_valley/%s.save/data-file-schema.xml --ref valley_ref.json "
        "--tol-mev %.6g --mode %s" % (pre, float(P["IVS_VALLEY_TOL_MEV"]), _valley_mode(P)),
        "cp -f tmp_valley/%s.save/data-file-schema.xml valley_bands.xml" % pre,
        "rm -rf tmp_valley",
        "printf 'OK\\n' > ivs_scf.ok"])
    _write(dest / "scf.in", scf)
    _write(dest / "valley.in.tpl", valley)
    _write(dest / "valley_ref.json", json.dumps(iv.jsonable(ref), indent=1, ensure_ascii=False) + "\n")
    _copy_job_helper(dest)
    render_submit(dest, "%s-ivs-scf" % pre, cmd, P)
    print("[DONE] %s：SCF %s k 网格、ecutwfc %.0f Ry、泛函 %s；谷能带 %d 个 k（自检容差 %.0f meV，%s）"
          % (dest, "×".join(map(str, kg)), float(P["IVS_ECUTWFC"]), q1["qe"]["input_dft"] or "赝势自带", len(pts),
             float(P["IVS_VALLEY_TOL_MEV"]), _valley_mode(P)))
    return 0


def _valley_mode(P):
    m = str(P.get("IVS_VALLEY_CHECK") or "error").strip().lower()
    if m not in ("error", "warn"):
        raise SystemExit("[ERROR] IVS_VALLEY_CHECK=%r 无效，只允许 error / warn" % m)
    return m


# ---------------------------------------------------------------------------
# S12：每个 q 一个子目录
# ---------------------------------------------------------------------------
def ahc_window(nocc, P, need=0):
    """AHC 带窗：带边上下各 nb 条（nb = max(IVS_AHC_NB, need)）。[0059-r3] need = 这个 q 的初态带边简并组要的条数
    （Si 空穴 Γ 三重简并 -> 3），只放宽需要的那个 q，别的 q 的窗口不变（已算完的不因此重算）。"""
    nb = int(P["IVS_AHC_NB"])
    if nb < 1:
        raise SystemExit("[ERROR] IVS_AHC_NB 必须 ≥ 1")
    nb = max(nb, int(need or 0))
    skip = max(0, int(nocc) - nb)
    nahc = int(nocc) - skip + nb
    nbnd = int(nocc) + nb + max(2, int(P["IVS_NBND_EXTRA"]))
    return skip, nahc, nbnd


def edge_needs(nocc):
    """[0059-r3] S11 的谷能带（valley_bands.xml + valley_ref.json 的 k 点）-> [(k, 电子要的 nb, 空穴要的 nb)]：
    带边简并组（|ΔE| < 1 meV 连成链）有几条带。没有这两个文件（旧版 S11 的结果没拉回来）-> None。"""
    iv = _iv()
    d = Path(STEP["scf"])
    if not ((d / "valley_bands.xml").is_file() and (d / "valley_ref.json").is_file()):
        return None
    import numpy as np
    import ivs_qe_job as job
    pts = json.loads((d / "valley_ref.json").read_text(encoding="utf-8"))["points"]
    info = job._xml_info(str(d / "valley_bands.xml"))
    ks = info["ks"]
    # 谷能带的 k（笛卡尔 2π/alat）必须就是 valley_ref.json 的点：S11 重新生成后还没跑完时两者会对不上
    B = np.linalg.inv(np.asarray(info["cell_bohr"], float)).T * float(info["alat_bohr"])
    if len(ks) != len(pts) or any(np.max(np.abs(np.asarray(k, float) @ B - np.asarray(row["k_cart_2pi_alat"]))) > 1e-4
                                  for k, row in zip(pts, ks)):
        raise SystemExit("[ERROR] %s/valley_bands.xml 的 k 点与 valley_ref.json 对不上：S11 重新生成后还没跑完（或结果没拉回来）。"
                         "先 start S11_ivsscf 跑完，再生成本步" % STEP["scf"])
    out = []
    for k, row in zip(pts, ks):
        e = row["eig_eV"]
        ge = iv.degenerate_group(e, nocc)
        gh = iv.degenerate_group(e, nocc - 1)
        out.append((k, max(ge) - nocc + 1, nocc - min(gh)))
    return out


def run_need(q1, run, needs):
    """[0059-r3] 一个 q 作业的初态带边要多宽的 AHC 窗口；needs 为 None -> 0（用 IVS_AHC_NB）。
    初态 k 不在 S11 的谷能带 k 点里 = S11 是按旧谷表跑的（S10 后来加了载流子/谷）-> 报错要求 retry S11。"""
    if needs is None:
        return 0
    iv = _iv()
    need = 0
    for c in q1["classes"]:
        if c["run"] != run["id"]:
            continue
        hit = next((n for n in needs if iv.same_q(n[0], c["k_init"], 1e-5)), None)
        if hit is None:
            raise SystemExit("[ERROR] %s 的初态 k (%s)（%s %s）不在 S11 谷能带的 k 点里：S11 是按旧的 S10 谷表跑的"
                             "（后来加了载流子或谷，如 IVS_ODP 的单谷空穴）。先 retry + start S11_ivsscf（SCF + 谷能差自检，"
                             "已算完的 q 不受影响），再生成本步"
                             % (run["name"], ", ".join("%.4f" % x for x in c["k_init"]), c["carrier"], c["from"]))
        need = max(need, hit[1] if c["carrier"] == "electron" else hit[2])
    return need


def odp_polar_mode(P):
    m = str(P.get("IVS_ODP_POLAR") or "exclude").strip().lower()
    if m not in ("exclude", "include"):
        raise SystemExit("[ERROR] IVS_ODP_POLAR=%r 无效，只允许 exclude / include" % m)
    return m


def gen_ph(P):
    iv = _iv()
    q1 = _read_json(Path(STEP["qlist"]) / "ivs_qlist.json", "先跑 S10_ivsq")
    info = _read_json(Path(STEP["scf"]) / "ivs_scf_info.json", "S11_ivsscf 的作业要先跑完")
    if not (Path(STEP["scf"]) / "ivs_scf.ok").is_file():
        raise SystemExit("[ERROR] S11_ivsscf 还没完成（没有 ivs_scf.ok）")
    pre = q1["qe"]["prefix"]
    nocc = int(info["nocc"])
    needs = edge_needs(nocc)
    if needs is None:
        print("[WARN] 没有 S11 的谷能带（%s/valley_bands.xml、valley_ref.json）：AHC 带窗按 IVS_AHC_NB，"
              "带边简并多于它时 S13 会报错" % STEP["scf"])
    polar = odp_polar_mode(P)
    nat = len(q1["structure"]["symbols"])
    types = q1["structure"]["types"]
    dest = Path(STEP["ph"])
    dest.mkdir(parents=True, exist_ok=True)
    names = []
    for r in q1["runs"]:
        name = r["name"]
        names.append(name)
        skip, nahc, nbnd = ahc_window(nocc, P, run_need(q1, r, needs))
        gamma = all(abs(x) < 1e-9 for x in r["q"])
        # [0059-r3] q=Γ 的谷内 ODP 要判极性模：DFPT 顺带算 Born 有效电荷（epsil；要有带隙）。单元素体系按非极性处理，不算
        epsil = bool(gamma and polar == "exclude" and len(types) > 1
                     and any(c["kind"] == "intra" for c in q1["classes"] if c["run"] == r["id"]))
        xq = iv.xq_alat(r["q"], info["cell_bohr"], info["alat_bohr"])
        qrun = {"id": r["id"], "name": name, "q_frac": r["q"], "xq_2pi_alat": xq, "k_frac": r["k"], "nk": len(r["k"]),
                "nbnd": nbnd, "ahc_nbndskip": skip, "ahc_nbnd": nahc, "nat": nat, "nocc": nocc, "prefix": pre,
                "classes": [c["id"] for c in q1["classes"] if c["run"] == r["id"]],
                "qlist_sha256": _sha_file(Path(STEP["qlist"]) / "ivs_qlist.json"),
                "scf_sha256": _sha_file(Path(STEP["scf"]) / "scf.in"),
                "ecutwfc": float(P["IVS_ECUTWFC"]), "tr2_ph": float(P["IVS_TR2_PH"]), "epsil": epsil}
        d = dest / name
        key = {k: qrun[k] for k in ("q_frac", "k_frac", "nbnd", "ahc_nbndskip", "ahc_nbnd", "scf_sha256", "ecutwfc", "tr2_ph",
                                    "epsil")}
        old = d / "qrun.json"
        if old.is_file():
            try:
                prev = json.loads(old.read_text(encoding="utf-8"))
                pk = {k: prev.get(k) for k in key}
                if pk.get("tr2_ph") == 0.0:       # [0059-r1] 旧版把 1e-16 舍成 0.0 写进 qrun.json：当作没记录，不因此判输入变了
                    pk["tr2_ph"] = key["tr2_ph"]
                # [0068] 旧作业收敛得比现在要求的更严（出厂 1e-16 -> 1e-14）：结果照样有效，不因此移走重算
                #   （包括还在跑的 q：重新生成会把它的目录挪走）；比现在松才算输入变了。
                try:
                    if pk.get("tr2_ph") is not None and float(pk["tr2_ph"]) <= key["tr2_ph"]:
                        pk["tr2_ph"] = key["tr2_ph"]
                except (TypeError, ValueError):
                    pass
                if "epsil" not in prev:           # [0059-r3] 旧版没有这个键 = 没开 epsil
                    pk["epsil"] = False
                same = iv.jsonable(pk) == iv.jsonable(key)
            except ValueError:
                same = False
            if same and (d / "ivs_ph.ok").is_file():
                print("[..] %s 已完成且输入没变：不动" % name)
                continue
            if not same:
                bak = dest / (".old_%s_%s" % (name, time.strftime("%Y%m%d%H%M%S")))
                d.rename(bak)
                print("[WARN] %s 的输入变了（q/k/带窗/SCF）：旧目录移到 %s，重新生成" % (name, bak.name))
        d.mkdir(parents=True, exist_ok=True)
        xqs = "%.12f %.12f %.12f" % tuple(xq)
        ph = ("&inputph\n  prefix = '%s'\n  outdir = './tmp'\n  fildyn = 'dyn'\n  fildvscf = 'dvscf'\n"
              "  tr2_ph = %.3e\n  ldisp = .false.\n%s/\n%s\n"
              % (pre, float(P["IVS_TR2_PH"]), "  epsil = .true.\n" if epsil else "", xqs))
        ahc = ("&inputph\n  prefix = '%s'\n  outdir = './tmp'\n  fildyn = 'dyn.ahc'\n  fildvscf = 'dvscf'\n"
               "  ldisp = .false.\n  trans = .false.\n  electron_phonon = 'ahc'\n  ahc_nbnd = %d\n  ahc_nbndskip = %d\n"
               "  skip_upperfan = .true.\n  ahc_dir = 'ahc_dir/'\n/\n%s\n" % (pre, nahc, skip, xqs))
        nscf = pw_input("nscf", P, q1, pre, "./tmp", iv.kpoints_list(r["k"]), nbnd=nbnd)
        chk = ("python3 ivs_qe_job.py ahc --dir ahc_dir --nbnd %d --ahc-nbnd %d --nat %d --nk %d"
               % (nbnd, nahc, nat, len(r["k"])))
        scf_save = "../../%s/tmp/%s.save" % (STEP["scf"], pre)
        cmd = "\n".join([
            "rm -f ivs_ph.ok",
            _pseudo_links(types, q1["qe"]["pseudos"]),
            "if ! %s >/dev/null 2>&1; then" % chk,
            "  if [ ! -s dyn ] || ! grep -q 'JOB DONE' ph_dfpt.out 2>/dev/null || [ ! -d tmp/_ph0/%s.phsave ]; then" % pre,
            "    rm -rf tmp && mkdir -p tmp && cp -r %s tmp/" % scf_save,
            "    run_mpi ph.x -in ph_dfpt.in > ph_dfpt.out",
            "    grep -q 'JOB DONE' ph_dfpt.out",
            "  fi",
            "  run_mpi pw.x -in nscf.in > nscf.out",
            "  grep -q 'JOB DONE' nscf.out",
            "  rm -rf ahc_dir",
            "  run_mpi ph.x -in ph_ahc.in > ph_ahc.out",
            "  grep -q 'JOB DONE' ph_ahc.out",
            "  %s" % chk,
            "fi",
            "" if P.get("IVS_KEEP_TMP") else "rm -rf tmp",
            "printf 'OK\\n' > ivs_ph.ok"])
        _write(d / "ph_dfpt.in", ph)
        _write(d / "ph_ahc.in", ahc)
        _write(d / "nscf.in", nscf)
        _write(d / "qrun.json", json.dumps(iv.jsonable(qrun), indent=1) + "\n")
        _copy_job_helper(d)
        render_submit(d, "%s-ivs-%s" % (pre, name), cmd, P)
        print("[DONE] %s：q = (%s)，xq = (%s) 2π/alat，初态 k %d 个，nbnd %d，AHC 带 %d..%d%s"
              % (name, ", ".join("%.4f" % x for x in r["q"]), xqs, len(r["k"]), nbnd, skip + 1, skip + nahc,
                 "，epsil（Born 有效电荷判极性模）" if epsil else ""))
    for d in sorted(dest.glob(FANOUT_PREFIX + "*")):
        if d.is_dir() and d.name not in names:
            bak = dest / (".old_%s_%s" % (d.name, time.strftime("%Y%m%d%H%M%S")))
            d.rename(bak)
            print("[WARN] %s 不在现在的 q 列表里：移到 %s" % (d.name, bak.name))
    return 0


# ---------------------------------------------------------------------------
# S13：读结果
# ---------------------------------------------------------------------------
def gen_read(P):
    import numpy as np
    iv = _iv()
    q1 = _read_json(Path(STEP["qlist"]) / "ivs_qlist.json", "先跑 S10_ivsq")
    info = _read_json(Path(STEP["scf"]) / "ivs_scf_info.json", "S11_ivsscf 没跑完")
    vchk = json.loads((Path(STEP["scf"]) / "valley_check.json").read_text(encoding="utf-8")) \
        if (Path(STEP["scf"]) / "valley_check.json").is_file() else None
    env = json.loads((Path(STEP["scf"]) / "ivs_env.json").read_text(encoding="utf-8")) \
        if (Path(STEP["scf"]) / "ivs_env.json").is_file() else {}
    by_run, problems = {}, []
    for r in q1["runs"]:
        d = Path(STEP["ph"]) / r["name"]
        if not (d / "ivs_ph.ok").is_file():
            problems.append("%s 没完成（没有 ivs_ph.ok）" % r["name"])
            continue
        qr = json.loads((d / "qrun.json").read_text(encoding="utf-8"))
        try:
            gkk, etk, etq = iv.read_ahc(str(d / "ahc_dir"), qr["nbnd"], qr["ahc_nbnd"], qr["nat"], qr["nk"])
            dyn = iv.read_dyn(str(d / "dyn"))
        except (OSError, ValueError) as e:
            problems.append("%s：%s" % (r["name"], e))
            continue
        dq = np.asarray(dyn["q_list"][0]) - np.asarray(qr["xq_2pi_alat"])
        if np.max(np.abs(dq)) > 1e-5:
            problems.append("%s：dyn 第一个 q %s 与输入 xq %s 不同" % (r["name"], dyn["q_list"][0].tolist(), qr["xq_2pi_alat"]))
            continue
        hw, e = iv.phonon_modes(dyn["C"][0], dyn["masses_ry"])
        by_run[r["id"]] = {"qr": qr, "gkk": gkk, "etk": etk, "etq": etq, "hw": hw, "e": e, "masses": dyn["masses_ry"],
                           "zstar": dyn.get("zstar"), "epsilon": dyn.get("epsilon")}
    if problems:
        for p in problems:
            print("[ERROR] " + p)
        raise SystemExit("[ERROR] S12_ivsph 的产物不全/不对（见上）：retry S12_ivsph 只补没完成的 q")
    odp_info = {}
    for rid, R in by_run.items():                  # [0059-r3] q=Γ（谷内 ODP）：剔声学模与极性模
        if not any(c["kind"] == "intra" for c in q1["classes"] if c["run"] == rid):
            continue
        R["exclude"], odp_info = odp_exclusions(R, q1, P)
    out_classes = []
    for c in q1["classes"]:
        R = by_run[c["run"]]
        qr = R["qr"]
        try:
            res = iv.class_couplings(R["gkk"][c["ik"]], R["etk"][c["ik"]], R["etq"][c["ik"]], R["hw"], R["e"], R["masses"],
                                     qr["nocc"], qr["ahc_nbndskip"], qr["ahc_nbnd"], c["carrier"],
                                     exclude=R.get("exclude") if c["kind"] == "intra" else None)
        except ValueError as e:
            raise SystemExit("[ERROR] 等价类 %d（%s %s→%s）：%s" % (c["id"], c["carrier"], c["from"], c["to"], e))
        res.update({k: c[k] for k in ("id", "carrier", "from", "to", "kind", "k_init", "k_fin", "q", "n_pairs", "pairs")})
        res["run"] = qr["name"]
        out_classes.append(res)
    carriers = {}
    for car, ci in q1["carriers"].items():
        cl = [c for c in out_classes if c["carrier"] == car]
        chs, dropped = iv.fold_channels(cl)
        carriers[car] = {"status": "ok", "valleys": ci["valleys"], "valleys_source": ci["valleys_source"],
                         "use_grid": ci["use_grid"], "valley_sig": ci["valley_sig"], "valley_sig_full": ci["valley_sig_full"],
                         "channels": chs, "dropped": dropped,
                         "valley_check": (vchk or {}).get("carriers", {}).get(car)}
    params = {"schema": iv.SCHEMA, "created": time.strftime("%Y-%m-%d %H:%M:%S"), "dim": q1["dim"],
              "structure_sha256": q1["structure"]["sha256"], "structure_source": q1["structure"]["source"],
              "valleys_source_sha256": (q1.get("inputs") or {}).get("dpt_result_sha256"),
              "scope": ("零阶谷间形变势（DFPT 直接矩阵元，非 SOC；不含一阶谷间耦合）"
                        + ("＋谷内零阶光学形变势 ODP（q=Γ）" if odp_info else "")),
              "provenance": {"qe_version": info.get("qe_version"), "input_dft": q1["qe"]["input_dft"],
                             "input_dft_note": q1["qe"]["input_dft_note"], "ecutwfc_Ry": q1["qe"]["ecutwfc_Ry"],
                             "ecutrho_Ry": q1["qe"]["ecutrho_Ry"], "kgrid": q1["qe"]["kgrid"], "tr2_ph": q1["qe"]["tr2_ph"],
                             "pseudos": {t: {"file": q1["qe"]["pseudos"][t],
                                             "sha256": ((env.get("pseudos") or {}).get(t) or {}).get("sha256"),
                                             "kind": ((env.get("pseudos") or {}).get(t) or {}).get("kind")}
                                         for t in q1["structure"]["types"]},
                             "method": "ph.x electron_phonon='ahc'（trans=.false.，单 q），|g| 按初/终态简并组平均"},
              "checks": {"valley_check": (vchk or {}).get("status"), "valley_check_mode": (vchk or {}).get("mode")},
              "odp": odp_info or None, "carriers": carriers, "raw": out_classes}
    dest = Path(STEP["read"])
    _write(dest / "ivs_dfpt_params.json", json.dumps(iv.jsonable(params), indent=1, ensure_ascii=False) + "\n")
    _write(dest / "ivs_dfpt_table.txt", table_text(params))
    print(table_text(params))
    print("[DONE] %s/ivs_dfpt_params.json（schema %s）。ke-dft-cpu 要用它：结果拉回本地后对 S8_kappa（3D）/ S8.4_amset2d（2D）"
          "执行 retry（gen 时自动推过去，IVS_DFPT=auto 读取）" % (dest, iv.SCHEMA))
    return 0


def odp_exclusions(R, q1, P):
    """[0059-r3] q=Γ 作业 -> ({模下标（0 基）: 不进谷内通道的原因}, 记录)。
    声学模：平移权重最大的 3 个（ħω≈0，谷内声学耦合由 ADP 管；数值上不严格为 0，r3 初版在 Si 上漏进一条 0.58 meV 的）。
    极性模：IVS_ODP_POLAR=exclude 时按 Born 有效电荷（q-odp 的 DFPT 开 epsil）的模有效电荷 > IVS_ODP_POLAR_TOL 判，
    长程 Fröhlich 耦合由 POP 管（与冻声子 ODP_POLAR 同口径）。单元素体系没有 Born 电荷时按非极性处理。"""
    iv = _iv()
    ac, tw = iv.acoustic_modes(R["e"], R["masses"])
    ex = {i: "声学模（q=Γ，平移权重 %.2f）" % tw[i] for i in ac}
    polar = odp_polar_mode(P)
    tol = float(P.get("IVS_ODP_POLAR_TOL") or iv.POLAR_TOL_E)
    rec = {"acoustic_modes": [i + 1 for i in ac], "polar": polar, "polar_tol_e": tol, "polar_modes": [],
           "zmode_e": None, "born_source": None}
    if polar == "exclude":
        if R.get("zstar") is not None:
            two_d = q1["dim"] == "2d"
            zm = iv.mode_charges(R["e"], R["masses"], R["zstar"], two_d=two_d)
            rec["zmode_e"] = [round(float(z), 4) for z in zm]
            rec["born_source"] = "ph.x epsil（%s）" % R["qr"]["name"]
            for i, z in enumerate(zm):
                if i not in ex and z > tol:
                    ex[i] = "%s极性模（%s模有效电荷 %.2f e > %.2f；长程耦合由 POP 管，IVS_ODP_POLAR=include 才算）" % (
                        "面内" if two_d else "", "面内" if two_d else "", z, tol)
                    rec["polar_modes"].append(i + 1)
        elif len(q1["structure"]["types"]) == 1:
            rec["born_source"] = "单元素体系，按非极性处理（没有 Born 有效电荷）"
        else:
            raise SystemExit("[ERROR] %s 的 dyn 里没有 Born 有效电荷（DFPT 没开 epsil：旧版生成的输入？），判不了哪些模是极性的。"
                             "retry + start S12_ivsph（只重算 %s），或 step.conf 写 IVS_ODP_POLAR = include（极性模的短程"
                             "形变势也算进 ODP）" % (R["qr"]["name"], R["qr"]["name"]))
    return ex, rec


def table_text(params):
    lines = ["# ivs_dfpt_params.json 摘要（D 单位 eV/Å，ħω 单位 meV）", "# %s" % params["scope"]]
    for car, c in params["carriers"].items():
        lines.append("## %s：%s（谷表 %s）" % (car, ", ".join("%s×%d" % (v["label"], v["g"]) for v in c["valleys"]),
                                          c["valleys_source"]))
        lines.append("#   逐谷对（初态取每个星的第一个成员；D 已按初/终态简并组平均；声子简并组 D² 相加）")
        for r in [x for x in params["raw"] if x["carrier"] == car]:
            gr = ", ".join("%.2f@%.1f" % (g["D_eV_A"], g["hw_meV"]) for g in r["groups"])
            lines.append("  %s→%s (%s) q=(%s) ×%d 对  带 %s→%s：%s" % (r["from"], r["to"], r["run"],
                                                                 ", ".join("%.3f" % x for x in r["q"]), r["n_pairs"],
                                                                 r["init_bands"], r["final_bands"], gr))
        lines.append("#   dp_optical 星通道（X>Y：D²/g_Y；X>X'：D²/(g_X−1)；X>X 谷内 ODP 不除；同 ħω 合并）")
        for ch in c["channels"]:
            lines.append("  %s>%s%s:%.3f@%.2f" % (ch["from"], ch["to"], "'" if ch["kind"] == "inter_star" else "",
                                                ch["D_eV_A"], ch["hw_meV"]))
        if c.get("dropped"):
            lines.append("#   不进通道：%s" % "; ".join("%.1f meV（%s）" % (d["hw_meV"], d["why"]) for d in c["dropped"]))
    return "\n".join(lines) + "\n"


def main(kind):
    P = load_params()
    fn = {"qlist": gen_qlist, "scf": gen_scf, "ph": gen_ph, "read": gen_read}.get(kind)
    if fn is None:
        raise SystemExit("unknown kind %r" % kind)
    return fn(P)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
