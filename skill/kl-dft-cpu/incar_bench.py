#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""incar_bench.py —— 取力 INCAR 基准对照（S4_disp 的隔离 side-car）。

目的
----
正式 S4 扇出（disp-* 位移取力）**之前**，用【同一批位移】的极少数几帧
（默认 2 帧），比较几种 INCAR 变体的力一致性，从而为目标体系定下生产取力
INCAR。

判据（2026-09-19 升级为相对量）
------------------------------
默认：任意变体相对**参考设置**的相对力差
      |DF|_rms / |F|_rms < 1%（--threshold-rel，默认 1.0）
算通过。绝对判据 max|DF| < 0.5 meV/Å 保留为可选（--threshold 0.5），
两者可同时给出，相对判据优先。

参考设置（baseline）
-------------------
默认 = 变体中**精度最高**的那个，按规则选取：
  1) LREAL=.FALSE. + ADDGRID=.TRUE. 且 ENCUT 最大；
  2) 没有 ADDGRID=.TRUE. 的变体时，退化为 LREAL=.FALSE. + ENCUT 最大；
  3) 再没有则退化为 ENCUT 最大者。
可用 --baseline NAME 显式指定（可配置）。

变体（默认，见 VARIANTS）
-------------------------
    lreal_false_noaddgrid   LREAL=.FALSE.，无 ADDGRID，生产 ENCUT
    lreal_false_addgrid     LREAL=.FALSE. + ADDGRID=.TRUE.
    lreal_auto_noaddgrid    LREAL=Auto   ，无 ADDGRID
    lreal_auto_addgrid      LREAL=Auto   + ADDGRID=.TRUE.
    encut_1.5x              LREAL=.FALSE. + ADDGRID，ENCUT=1.5×max(ENMAX)
    encut_2.0x              LREAL=.FALSE. + ADDGRID，ENCUT=2.0×max(ENMAX)
    encut_650               LREAL=.FALSE. + ADDGRID，ENCUT=650（第三点，与 2.0x 只差 ENCUT）
共 7 个目录 × N 帧的极小 VASP 作业。

per-frame 记录字段（analyze report）
------------------------------------
    max_abs_delta_meV_per_A   相对参考设置的逐分量最大绝对差 max|DF|
    rms_abs_delta_meV_per_A   相对参考设置的逐分量差 RMS |DF|_rms
    force_rms_meV_per_A       该帧自身的力模 RMS |F|_rms
    ref_force_rms_meV_per_A   参考帧的力模 RMS
    rel_error_pct             100 * |DF|_rms / |F|_rms(ref)
    cpu_time_s                OUTCAR "Total CPU time used (sec)"
    elapsed_s                 OUTCAR "Elapsed time (sec)"
    node                      作业节点（jzzn：sacct/squeue -o %N，scontrol NodeList 为 null）
    jobid                     提交时记录或从 sacct 发现

★ 隔离性（硬要求）
------------------
bench 帧一律放 step4_disp/bench/<variant>/disp-XXXXX/，**不占用**生产
disp-* 命名空间（生产帧在 step4_disp/disp-*，bench 在下一层 bench/ 里）。
生产扇出与 S5 帧校验用的都是**非递归** glob：

    - autozt fanout:      glob(os.path.join(step_dir, "disp-*"))      # workflow.py
    - kl_common.check_frames_match_displacements: d4.glob("disp-*")   # 非递归
    - gen_step4_disp.check_existing_matches_input: glob(out/"disp-*") # 非递归
    - kl_fc_backends.collect_vaspruns: glob(disp_dir/"disp-*")        # 非递归
    - gen_step5_fc: disp.glob("disp-*/vasprun.xml")                   # 非递归

仓库里**没有任何** rglob / ** 递归遍历 kl 的 disp-*（已 grep 确认），所以
bench/ 天然落在生产扇出与 S5 校验之外。本模块只写 bench/ 子树。

用法
----
    python incar_bench.py generate [--dir step4_disp] [--frames 2] [--force]
    python incar_bench.py submit   [--dir step4_disp]      # 提交 bench 作业
    python incar_bench.py nodes    [--dir step4_disp]      # (重新)采集 jobid/节点
    python incar_bench.py analyze  [--dir step4_disp] [--baseline NAME]
                                   [--threshold-rel 1.0] [--threshold 0.5]
    python incar_bench.py          # == analyze

gen_step4_disp.py 在 INCAR_BENCH=on 时会直接 import 本模块调 generate()；
submit / nodes / analyze 走本文件的独立命令。
"""
import argparse
import datetime
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

BENCH_SUBDIR        = "bench"
PLAN_NAME           = "incar_bench_plan.json"
REPORT_NAME         = "incar_bench_report.json"
JOBS_NAME           = "incar_bench_jobs.json"
DEFAULT_BASELINE    = None          # None = 自动取精度最高的变体
DEFAULT_THRESHOLD   = 0.5           # meV/Å（绝对判据，可选）
DEFAULT_THRESHOLD_REL_PCT = 1.0     # %（相对判据，默认启用）
PROD_ENCUT_FACTOR   = 1.5          # ENCUT_FACTOR 的口径（生产取力）
FACTOR_LOW          = 1.5          # encut_1.5x
FACTOR_HIGH         = 2.0          # encut_2.0x（CELL_STAGE_ENCUT_FACTOR 的口径）
FIXED_ENCUT_HIGH    = 650          # encut_650 的固定 ENCUT（第三点）

_TAG_SYSTEM  = "SYSTEM"
_TAG_ENCUT   = "ENCUT"
_TAG_LREAL   = "LREAL"
_TAG_ADDGRID = "ADDGRID"

# 7 个变体：前 4 个是 LREAL × ADDGRID 的 2×2，后 3 个是 ENCUT 对照。
#   encut_*  固定 LREAL=.FALSE. + ADDGRID=.TRUE.（生产口径），只动 ENCUT；
#   encut_1.5x 对 MoS2（max ENMAX=260）刚好 = 生产 ENCUT=390，与
#   lreal_false_addgrid 的 INCAR 物理设置相同，作为"同输入复现性对照"；
#   encut_650 是第三点，用于判断 520 是否已收敛。
VARIANTS = [
    {"name": "lreal_false_noaddgrid", "lreal": ".FALSE.", "addgrid": None,
     "encut": "prod", "desc": "LREAL=.FALSE.，无 ADDGRID，生产 ENCUT"},
    {"name": "lreal_false_addgrid", "lreal": ".FALSE.", "addgrid": ".TRUE.",
     "encut": "prod", "desc": "LREAL=.FALSE. + ADDGRID=.TRUE."},
    {"name": "lreal_auto_noaddgrid", "lreal": "Auto", "addgrid": None,
     "encut": "prod", "desc": "LREAL=Auto，无 ADDGRID"},
    {"name": "lreal_auto_addgrid", "lreal": "Auto", "addgrid": ".TRUE.",
     "encut": "prod", "desc": "LREAL=Auto + ADDGRID=.TRUE."},
    {"name": "encut_1.5x", "lreal": ".FALSE.", "addgrid": ".TRUE.",
     "encut": "factor", "factor": FACTOR_LOW,
     "desc": "生产 LREAL/ADDGRID，ENCUT=1.5×max(ENMAX)（ENCUT_FACTOR 口径）"},
    {"name": "encut_2.0x", "lreal": ".FALSE.", "addgrid": ".TRUE.",
     "encut": "factor", "factor": FACTOR_HIGH,
     "desc": "生产 LREAL/ADDGRID，ENCUT=2.0×max(ENMAX)（CELL_STAGE_ENCUT_FACTOR 口径）"},
    {"name": "encut_650", "lreal": ".FALSE.", "addgrid": ".TRUE.",
     "encut": FIXED_ENCUT_HIGH,
     "desc": "生产 LREAL/ADDGRID，ENCUT=650（第三点，与 encut_2.0x 只差 ENCUT）"},
]

_STRESS_RE = re.compile(
    r"^\s*in\s+kB\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", re.M)


# ==========================================================================
# 通用小工具
# ==========================================================================
def enabled(val):
    """INCAR_BENCH 开关解析：on/true/1/yes/.true. 都算开，其余（含 off/空）算关。"""
    return str(val if val is not None else "").strip().lower() in (
        "on", "true", ".true.", "yes", "1", "enable", "enabled")


def _conf_get(conf, key, default=None):
    """同时兼容 stepconf.StepConf 对象与普通 dict；缺键返回 default。"""
    if conf is None:
        return default
    try:
        v = conf[key]
    except Exception:                                       # noqa: BLE001
        try:
            v = conf.get(key)
        except Exception:                                   # noqa: BLE001
            v = default
    return default if v is None else v


def _incar_tag(text, key):
    """从 INCAR 文本取某个标签的值（忽略 # / ! 注释），找不到返回 None。"""
    for ln in text.splitlines():
        code = ln.split("#")[0].split("!")[0]
        if "=" not in code:
            continue
        k, v = code.split("=", 1)
        if k.strip().upper() == key.upper():
            return v.strip()
    return None


def _incar_int_tag(text, key, default=None):
    v = _incar_tag(text, key)
    if v is None:
        return default
    m = re.match(r"[-+]?\d+", v)
    return int(m.group(0)) if m else default


def _is_lreal_false(v):
    return str(v if v is not None else "").strip().strip(".").lower() in (
        "false", "f", "no", "0")


def _is_addgrid_true(v):
    return str(v if v is not None else "").strip().strip(".").lower() in (
        "true", "t", "yes", "1")


def select_baseline(metas):
    """按"精度最高"规则选参考变体名。

    规则：① LREAL=.FALSE. + ADDGRID=.TRUE. 中 ENCUT 最大者；
          ② 没有 ADDGRID=.TRUE. 的变体时，退化为 LREAL=.FALSE. + ENCUT 最大者；
          ③ 再没有则退化为 ENCUT 最大者；返回变体名或 None。
    metas: [{"name":..., "LREAL":..., "ADDGRID":..., "ENCUT":...}, ...]
    """
    items = [m for m in metas if m.get("name")]

    def _pick(pred):
        cand = [m for m in items if pred(m) and m.get("ENCUT") is not None]
        if not cand:
            return None
        return max(cand, key=lambda m: int(m["ENCUT"]))

    picked = _pick(lambda m: _is_lreal_false(m.get("LREAL"))
                   and _is_addgrid_true(m.get("ADDGRID")))
    if picked is None:
        picked = _pick(lambda m: _is_lreal_false(m.get("LREAL")))
    if picked is None:
        picked = _pick(lambda m: True)
    return picked.get("name") if picked else None


def _load_json(p):
    try:
        if Path(p).is_file():
            return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:                                       # noqa: BLE001
        pass
    return {}


def _max_enmax(potcar):
    vals = []
    for ln in Path(potcar).read_text(errors="ignore").splitlines():
        m = re.search(r"ENMAX\s*=\s*([\d.]+)", ln)
        if m:
            vals.append(float(m.group(1)))
    return max(vals) if vals else None


def _encut_from_potcar(potcar, factor):
    """ceil(max(ENMAX) * factor / 10) * 10 —— 与 kl_common.encut_from_potcar 同口径。"""
    mx = _max_enmax(potcar)
    if mx is None:
        return None
    return int(math.ceil(mx * float(factor) / 10.0) * 10)


def _render_incar(base_text, system, encut, lreal, addgrid):
    """在【生产帧 INCAR】基础上改 SYSTEM/ENCUT/LREAL/ADDGRID，其余一字不动。

    这样偶极修正、VDW、精度、并行参数等全部继承生产口径，变体只差目标标签。
    ADDGRID=None 表示显式【删除】该行（"不加 ADDGRID"）。
    """
    drop = {_TAG_SYSTEM, _TAG_ENCUT, _TAG_LREAL, _TAG_ADDGRID}
    keep = []
    for ln in base_text.splitlines():
        code = ln.split("#")[0].split("!")[0]
        if "=" in code:
            key = code.split("=", 1)[0].strip().upper()
            if key in drop:
                continue
        keep.append(ln.rstrip())
    while keep and not keep[-1].strip():
        keep.pop()
    lines = keep + ["",
                    "SYSTEM  = %s" % system,
                    "ENCUT   = %s" % encut,
                    "LREAL   = %s" % lreal]
    if addgrid:
        lines.append("ADDGRID = %s" % addgrid)
    return "\n".join(lines) + "\n"


def _set_jobname(text, name):
    """把提交脚本的 #SBATCH --job-name= 换成 bench 专用名；没有该行就插在首个 #SBATCH 后。"""
    new, n = re.subn(r"^(#SBATCH\s+--job-name=).*$", r"\g<1>" + name,
                     text, flags=re.M)
    if n:
        return new
    lines = text.splitlines()
    idx = next((i for i, ln in enumerate(lines) if ln.startswith("#SBATCH")), None)
    if idx is not None:
        lines.insert(idx, "#SBATCH --job-name=%s" % name)
        return "\n".join(lines) + "\n"
    return text


def _production_frames(out, n_frames):
    """生产位移帧（disp-00000 是平衡帧，排除）。返回 [(num, Path), ...]，按编号排序。

    只 glob out/disp-*（非递归），bench/ 不会进来。
    """
    fr = []
    for d in sorted(Path(out).glob("disp-*")):
        m = re.fullmatch(r"disp-(\d+)", d.name)
        if not m or int(m.group(1)) == 0:
            continue
        if (d / "POSCAR").is_file() and (d / "INCAR").is_file():
            fr.append((int(m.group(1)), d))
    fr.sort(key=lambda x: x[0])
    if n_frames and int(n_frames) > 0:
        fr = fr[:int(n_frames)]
    return fr


# ==========================================================================
# generate：从生产帧生成 bench 输入
# ==========================================================================
def _variant_encut(v, prod_encut, potcar):
    mode = v.get("encut")
    if mode == "prod":
        return int(prod_encut) if prod_encut else None
    if isinstance(mode, (int, float)) and not isinstance(mode, bool):
        return int(mode)
    fac = float(v["factor"])
    if potcar is not None and Path(potcar).is_file():
        e = _encut_from_potcar(potcar, fac)
        if e:
            return e
    if prod_encut:
        # 没有 POTCAR：按生产 ENCUT_FACTOR=1.5 的口径折算
        return int(round(float(prod_encut) * fac / PROD_ENCUT_FACTOR / 10.0) * 10)
    return None


def _apply_overrides(text, overrides):
    """在生成出的 INCAR 文本上替换/追加 KEY=VALUE（overrides: {KEY: VALUE}）。"""
    if not overrides:
        return text
    keys = {str(k).upper() for k in overrides}
    out = []
    for ln in text.splitlines():
        code = ln.split("#")[0].split("!")[0]
        if "=" in code and code.split("=", 1)[0].strip().upper() in keys:
            continue
        out.append(ln.rstrip())
    while out and not out[-1].strip():
        out.pop()
    for k, v in overrides.items():
        out.append("%-8s= %s" % (str(k).upper(), v))
    return "\n".join(out) + "\n"


def generate(out, conf=None, dim="3d", n_frames=2, force=False, verbose=True,
             overrides=None, only=None, suffix=""):
    """在 <out>/bench/<variant><suffix>/disp-XXXXX 下生成基准帧。

    overrides: 额外/覆盖的 INCAR 标签，如 {"NELM": 300, "EDIFF": "1E-7"}（诊断腿用）。
    only:      只生成这些变体名（逗号串或 list）；缺省全部 VARIANTS。
    suffix:    变体目录/名字后缀（如 "_nelm300"），与主对照隔离命名。
    多次调用会【合并】进已有 plan，不会丢掉此前生成的变体。
    位移完全复用生产帧（逐位一致）；幂等：已有 INCAR+POSCAR 且非 force 则跳过。
    """
    out = Path(out)
    overrides = dict(overrides or {})
    if isinstance(only, str):
        only = [x.strip() for x in only.split(",") if x.strip()]
    only_set = set(only) if only else None

    def log(msg):
        if verbose:
            print(msg)

    if not out.is_dir():
        log("[WARN] INCAR_BENCH：%s 不存在，跳过" % out)
        return None
    prod = _production_frames(out, n_frames)
    if not prod:
        log("[WARN] INCAR_BENCH：%s 下找不到生产位移帧（disp-NNNNN/POSCAR+INCAR），跳过" % out)
        return None

    vs = [v for v in VARIANTS if (only_set is None or v["name"] in only_set)]
    if not vs:
        log("[WARN] INCAR_BENCH：--only=%s 没匹配到任何变体，跳过" % (only_set,))
        return None

    base_text = (prod[0][1] / "INCAR").read_text(encoding="utf-8", errors="ignore")
    potcar = prod[0][1] / "POTCAR"
    potcar = potcar if potcar.is_file() else None
    factor = float(_conf_get(conf, "ENCUT_FACTOR", PROD_ENCUT_FACTOR))
    prod_encut = _conf_get(conf, "ENCUT")
    if prod_encut is None and potcar is not None:
        prod_encut = _encut_from_potcar(potcar, factor)
    if prod_encut is None:
        prod_encut = _incar_int_tag(base_text, _TAG_ENCUT)
    if prod_encut is None:
        log("[WARN] INCAR_BENCH：拿不到生产 ENCUT（无 POTCAR / 无 step.conf / INCAR 里没有），跳过")
        return None

    bench = out / BENCH_SUBDIR
    bench.mkdir(parents=True, exist_ok=True)
    plan = _load_json(bench / PLAN_NAME) or {}
    plan.setdefault("skill", "kl-dft-cpu")
    plan.setdefault("step", "step4_disp")
    plan.setdefault("variants", [])
    plan["generated_at"] = datetime.datetime.now().isoformat(timespec="seconds")
    plan["baseline_rule"] = ("LREAL=.FALSE.+ADDGRID=.TRUE. 中 ENCUT 最大"
                             "（无 ADDGRID=.TRUE. 则 LREAL=.FALSE.+ENCUT 最大）")
    plan["criterion_mode"] = "relative"
    plan.setdefault("threshold_rel_pct", DEFAULT_THRESHOLD_REL_PCT)
    plan.setdefault("threshold_meV_per_A", DEFAULT_THRESHOLD)
    plan["dim"] = str(dim).lower()
    plan["prod_encut_eV"] = int(prod_encut)
    plan["generated_frames"] = ["disp-%05d" % n for n, _ in prod]
    plan["note"] = "位移与生产集逐位一致（直接 copy 生产 disp-NNNNN/POSCAR）"
    existing = {v["name"]: v for v in plan["variants"] if v.get("name")}
    seen, dup_note = {}, []
    for v in vs:
        encut = _variant_encut(v, prod_encut, potcar)
        if encut is None:
            log("[WARN] INCAR_BENCH：变体 %s 的 ENCUT 无法确定，跳过" % v["name"])
            continue
        key = (v["lreal"], v.get("addgrid"), int(encut),
               json.dumps(overrides, sort_keys=True))
        dup = seen.get(key)
        seen.setdefault(key, v["name"])
        if dup:
            dup_note.append("%s==%s" % (v["name"], dup))
        dname = v["name"] + suffix
        vinfo = {"name": dname, "base_variant": v["name"], "desc": v.get("desc", ""),
                 "LREAL": v["lreal"], "ADDGRID": v.get("addgrid"),
                 "ENCUT": int(encut), "overrides": dict(overrides),
                 "frames": [], "duplicate_of": dup}
        for num, pd in prod:
            fd = bench / dname / ("disp-%05d" % num)
            fd.mkdir(parents=True, exist_ok=True)
            incar_p, poscar_p = fd / "INCAR", fd / "POSCAR"
            framename = "disp-%05d" % num
            if incar_p.is_file() and poscar_p.is_file() and not force:
                vinfo["frames"].append(framename)
                continue
            shutil.copyfile(pd / "POSCAR", poscar_p)       # ★ 逐位一致
            for f in ("KPOINTS", "POTCAR"):
                if (pd / f).is_file():
                    shutil.copyfile(pd / f, fd / f)
            incar_text = _apply_overrides(
                _render_incar(base_text,
                              "bench %s %s %s" % (out.parent.name, dname, framename),
                              int(encut), v["lreal"], v.get("addgrid")),
                overrides)
            incar_p.write_text(incar_text, encoding="utf-8", newline="\n")
            if (pd / "submit.sh").is_file():
                sub = (pd / "submit.sh").read_text(encoding="utf-8", errors="ignore")
                (fd / "submit.sh").write_text(
                    _set_jobname(sub, "bench-%s-%s" % (dname, framename)),
                    encoding="utf-8", newline="\n")
            vinfo["frames"].append(framename)
        existing[dname] = vinfo

    plan["variants"] = sorted(existing.values(), key=lambda x: x.get("name", ""))
    plan["baseline"] = select_baseline(plan["variants"]) or (
        plan["variants"][0]["name"] if plan["variants"] else None)
    (bench / PLAN_NAME).write_text(
        json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    log("[OK] INCAR_BENCH：本次生成 %d 变体 × %d 帧（plan 合计 %d 变体）→ %s/"
        % (len(vs), len(prod), len(plan["variants"]), bench))
    log("[..] INCAR_BENCH：参考设置（精度最高）= %s" % plan["baseline"])
    if dup_note:
        log("[..] INCAR_BENCH：INCAR 物理设置重复的变体（可作同输入复现性对照）：%s"
            % "; ".join(dup_note))
    return plan


# ==========================================================================
# OUTCAR 解析
# ==========================================================================
def parse_forces(outcar):
    """读 OUTCAR 里【最后一次】TOTAL-FORCE 块，返回 [[fx,fy,fz], ...]（eV/Å）或 None。"""
    p = Path(outcar)
    if not p.is_file():
        return None
    text = p.read_text(errors="ignore")
    i = text.rfind("TOTAL-FORCE")
    if i < 0:
        return None
    lines = text[i:].splitlines()
    data, started = [], False
    for ln in lines[1:]:
        s = ln.strip()
        if s and set(s) <= set("-="):
            if started:
                break
            started = True
            continue
        if not started:
            continue
        parts = s.split()
        if len(parts) < 6:
            break
        try:
            data.append([float(x) for x in parts[-3:]])
        except ValueError:
            break
    return data or None


def parse_elapsed(outcar):
    """OUTCAR 最后一次 Elapsed time (sec)（秒）；没有返回 None。"""
    p = Path(outcar)
    if not p.is_file():
        return None
    ms = re.findall(r"Elapsed time \(sec\):\s*([\d.]+)", p.read_text(errors="ignore"))
    return float(ms[-1]) if ms else None


def parse_cpu_time(outcar):
    """OUTCAR 最后一次 "Total CPU time used (sec)"（秒）；没有返回 None。"""
    p = Path(outcar)
    if not p.is_file():
        return None
    ms = re.findall(r"Total CPU time used \(sec\):\s*([\d.]+)",
                    p.read_text(errors="ignore"))
    return float(ms[-1]) if ms else None


def parse_electronic_steps(outdir):
    """电子步数（DAV 迭代计数）+ 来源文件名；拿不到返回 (None, None)。

    DAV 行只出现在 VASP 标准输出里（submit.sh 的 --output=queue.out），不在 OUTCAR。
    依次尝试 queue.out / slurm-*.out / stdout，最后兜底在 OUTCAR 里找。
    """
    p = Path(outdir)
    if p.is_file():
        p = p.parent
    cands = [p / "queue.out"] + sorted(p.glob("slurm-*.out")) + [p / "stdout"]
    for f in cands:
        if f.is_file():
            n = len(re.findall(r"^DAV:", f.read_text(errors="ignore"), re.M))
            if n:
                return n, f.name
    oc = p / "OUTCAR"
    if oc.is_file():
        n = len(re.findall(r"^DAV:", oc.read_text(errors="ignore"), re.M))
        if n:
            return n, "OUTCAR"
    return None, None


_NELM_WARN = "number of steps (NELM)"
_ABORT_EDIFF = "aborting loop because EDIFF is reached"


def scf_converged(outcar):
    """电子 SCF 是否真收敛。返回 (bool|None, note)。

    收敛 = OUTCAR 含 "aborting loop because EDIFF is reached" 且不含 VASP 的
    "number of steps (NELM)" 警告（后者说明打满 NELM、力可能不可靠）。
    """
    p = Path(outcar)
    if not p.is_file():
        return None, "OUTCAR 缺失"
    txt = p.read_text(errors="ignore")
    reached = _ABORT_EDIFF in txt
    hit_nelm = _NELM_WARN in txt
    if reached and not hit_nelm:
        return True, ""
    if hit_nelm and not reached:
        return False, "打满 NELM，电子步未收敛（VASP 警告力可能不可靠）"
    if not reached and not hit_nelm:
        return False, "无收敛标志（未跑完 / 被取消）"
    return False, "同时含收敛与 NELM 警告标志（异常）"


def last_stress_kbar(outcar):
    """OUTCAR 最后一次 "in kB" 应力行 (xx,yy,zz,xy,yz,zx)，单位 kbar；没有返回 None。

    与 kl_common.last_stress_kbar 同口径（列序 XX YY ZZ XY YZ ZX）。
    """
    p = Path(outcar)
    if not p.is_file():
        return None
    ms = _STRESS_RE.findall(p.read_text(errors="ignore"))
    return tuple(float(x) for x in ms[-1]) if ms else None


def max_abs_delta(fa, fb, ev_per_a=True):
    """两条力表的逐分量最大绝对差。返回 (max|ΔF|, 分量)。单位 eV/Å 或 meV/Å。"""
    if not fa or not fb or len(fa) != len(fb):
        return None, None
    worst, which = 0.0, None
    for ia, (ra, rb) in enumerate(zip(fa, fb)):
        for c in range(3):
            d = abs(ra[c] - rb[c])
            if d > worst:
                worst, which = d, (ia, "xyz"[c])
    if ev_per_a:
        worst *= 1000.0                                    # eV/Å -> meV/Å
    return worst, which


def rms_delta(fa, fb):
    """两条力表的逐分量差 RMS（eV/Å）；不一致返回 None。"""
    if not fa or not fb or len(fa) != len(fb):
        return None
    s = n = 0
    for ra, rb in zip(fa, fb):
        for c in range(3):
            d = ra[c] - rb[c]
            s += d * d
            n += 1
    return math.sqrt(s / n) if n else None


def force_rms(fa):
    """力模 RMS：sqrt(mean(|F_i|²))（eV/Å）；空返回 None。"""
    if not fa:
        return None
    s = n = 0
    for r in fa:
        for c in range(3):
            s += r[c] * r[c]
            n += 1
    return math.sqrt(s / n) if n else None


# ==========================================================================
# 作业 / 节点信息（jzzn：scontrol show job 的 NodeList 为 (null)，
# 所以用 squeue -o "%i %T %N" 与 sacct 的 NodeList）
# ==========================================================================
def _slurm_lookup(jobids):
    """按 jobid 查 {jid: {"node","state","elapsed"}}；sacct + squeue 双通道，只读。"""
    info = {}
    ids = [str(j) for j in jobids if j]
    if ids:
        try:
            p = subprocess.run(
                ["sacct", "-n", "-P", "--format=JobID,NodeList,State,Elapsed",
                 "-j", ",".join(ids)],
                capture_output=True, text=True, timeout=120)
            for ln in (p.stdout or "").splitlines():
                parts = ln.split("|")
                if len(parts) < 4 or parts[0] != parts[0].split(".")[0]:
                    continue
                jid = parts[0]
                d = info.setdefault(jid, {})
                node = parts[1].strip()
                if node and node not in ("(null)", "None", "N/A"):
                    d["node"] = node
                if parts[2].strip():
                    d["state"] = parts[2].strip()
                if parts[3].strip():
                    d["elapsed"] = parts[3].strip()
        except Exception:                                   # noqa: BLE001
            pass
    try:
        p = subprocess.run(["squeue", "-h", "-o", "%i %T %N"],
                           capture_output=True, text=True, timeout=60)
        for ln in (p.stdout or "").splitlines():
            parts = ln.split(None, 2)
            if not parts:
                continue
            jid = parts[0].split("_")[0]
            if jid not in info:
                continue
            if len(parts) > 1 and parts[1]:
                info[jid]["state"] = parts[1]
            if len(parts) > 2:
                node = parts[2].strip()
                if node not in ("", "(null)", "None", "N/A"):
                    info[jid]["node"] = node
    except Exception:                                       # noqa: BLE001
        pass
    return info


def collect_nodes(bench, verbose=True):
    """采集 bench 作业的 jobid/节点，写 incar_bench_jobs.json。

    先读已有 jobs 文件；再从 sacct（近 30 天，按 job-name bench-<variant>-disp-NNNNN
    发现历史作业）+ squeue（running/pending，%N 取节点）补齐 node/state；最后只保留
    bench 下真实存在的帧。返回 jobs dict。
    """
    bench = Path(bench)

    valid = set()
    if bench.is_dir():
        for vdir in bench.iterdir():
            if vdir.is_dir():
                for f in vdir.glob("disp-*"):
                    valid.add("%s/%s" % (vdir.name, f.name))

    jobs = _load_json(bench / JOBS_NAME) or {}
    # 历史作业发现（sacct）
    try:
        p = subprocess.run(
            ["sacct", "-n", "-P", "--starttime=now-30days",
             "--format=JobID,JobName,NodeList,State,Elapsed"],
            capture_output=True, text=True, timeout=120)
        for ln in (p.stdout or "").splitlines():
            parts = ln.split("|")
            if len(parts) < 5 or parts[0] != parts[0].split(".")[0]:
                continue
            m = re.match(r"^bench-(.+)-disp-(\d+)$", parts[1].strip())
            if not m:
                continue
            key = "%s/disp-%s" % (m.group(1), m.group(2))
            d = jobs.setdefault(key, {})
            d.setdefault("jobid", parts[0])
            d.setdefault("job_name", parts[1].strip())
            node = parts[2].strip()
            if node and node not in ("(null)", "None", "N/A"):
                d["node"] = node
            if parts[3].strip():
                d["state"] = parts[3].strip()
            if parts[4].strip():
                d["elapsed"] = parts[4].strip()
    except Exception:                                       # noqa: BLE001
        pass
    # 运行/排队作业补节点（squeue -o %N）
    try:
        p = subprocess.run(["squeue", "-h", "-o", "%i %T %N"],
                           capture_output=True, text=True, timeout=60)
        for ln in (p.stdout or "").splitlines():
            parts = ln.split(None, 2)
            if not parts:
                continue
            jid = parts[0].split("_")[0]
            for d in jobs.values():
                if str(d.get("jobid")) == jid:
                    if len(parts) > 1 and parts[1]:
                        d["state"] = parts[1]
                    if len(parts) > 2:
                        node = parts[2].strip()
                        if node not in ("", "(null)", "None", "N/A"):
                            d["node"] = node
    except Exception:                                       # noqa: BLE001
        pass
    if valid:
        jobs = {k: v for k, v in jobs.items() if k in valid}
    if bench.is_dir() and jobs:
        (bench / JOBS_NAME).write_text(
            json.dumps(jobs, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    if verbose:
        got = sum(1 for v in jobs.values() if v.get("node"))
        print("[..] INCAR_BENCH：记录 %d 帧 jobid，%d 帧有节点名 → %s/%s"
              % (len(jobs), got, bench, JOBS_NAME))
    return jobs


def _frame_jobid(jobs, name, fname):
    return (jobs.get("%s/%s" % (name, fname)) or {}).get("jobid")


def _frame_node(bench, name, fname, jobs):
    d = jobs.get("%s/%s" % (name, fname)) or {}
    if d.get("node"):
        return d["node"]
    p = Path(bench) / name / fname / "node.txt"
    if p.is_file():
        t = p.read_text(errors="ignore").strip().split()
        return t[0] if t else None
    return None


# ==========================================================================
# analyze：读 bench OUTCAR → 表 + JSON
# ==========================================================================
def _load_plan(bench):
    p = bench / PLAN_NAME
    if p.is_file():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:                                   # noqa: BLE001
            pass
    return {}


def _variant_dir_meta(vdir, plan_meta):
    meta = dict(plan_meta.get(vdir.name, {}))
    meta["name"] = vdir.name
    inc = vdir / next((f.name for f in sorted(vdir.glob("disp-*"))), "") / "INCAR"
    if inc.is_file():
        txt = inc.read_text(encoding="utf-8", errors="ignore")
        meta.setdefault("LREAL", _incar_tag(txt, _TAG_LREAL))
        ag = _incar_tag(txt, _TAG_ADDGRID)
        meta.setdefault("ADDGRID", ag if ag else None)
        meta.setdefault("ENCUT", _incar_int_tag(txt, _TAG_ENCUT))
        meta.setdefault("NELM", _incar_int_tag(txt, "NELM"))
        meta.setdefault("ALGO", _incar_tag(txt, "ALGO"))
        meta.setdefault("EDIFF", _incar_tag(txt, "EDIFF"))
    return meta


def _criterion_text(threshold_mev, threshold_rel_pct):
    if threshold_rel_pct is not None:
        txt = ("任意变体相对参考设置的 |ΔF|_rms / |F|_rms < %.2f%%（相对判据）"
               % float(threshold_rel_pct))
        if threshold_mev is not None:
            txt += "；绝对参考 max|ΔF| < %.2f meV/Å" % float(threshold_mev)
        return txt
    return "任意变体相对参考设置的最大力分量差 max|ΔF| < %.2f meV/Å" % float(threshold_mev)


def _judge(rec, threshold_mev, threshold_rel_pct):
    if threshold_rel_pct is not None:
        r = rec.get("rel_error_pct")
        if r is None:
            return "PENDING"
        return "PASS" if r < float(threshold_rel_pct) else "FAIL"
    if threshold_mev is not None:
        d = rec.get("max_abs_delta_meV_per_A")
        if d is None:
            return "PENDING"
        return "PASS" if d < float(threshold_mev) else "FAIL"
    return "PENDING"


def analyze(out, baseline=None, threshold_mev=None,
            threshold_rel_pct=DEFAULT_THRESHOLD_REL_PCT, verbose=True):
    """统计各变体相对【参考设置】的 max|ΔF|、RMS|ΔF|、相对误差、CPU time、节点。

    默认判据：|ΔF|_rms / |F|_rms < threshold_rel_pct%（默认 1%）。
    threshold_rel_pct=None 时退回绝对判据 max|ΔF| < threshold_mev meV/Å。
    baseline=None（或 "auto"）时自动选精度最高的变体。写 report JSON。
    """
    out = Path(out)
    bench = out / BENCH_SUBDIR
    plan = _load_plan(bench)

    if baseline in (None, "", "auto"):
        baseline = None
    if threshold_mev is None and threshold_rel_pct is None:
        threshold_rel_pct = DEFAULT_THRESHOLD_REL_PCT

    if not bench.is_dir():
        rep = {"baseline": baseline, "threshold_meV_per_A": threshold_mev,
               "threshold_rel_pct": threshold_rel_pct,
               "criterion_mode": ("relative" if threshold_rel_pct is not None else "absolute"),
               "variants": [], "verdict": "PENDING",
               "recommendation": "还没有 %s/ —— 先在 INCAR_BENCH=on 下跑一次 S4 gen。" % BENCH_SUBDIR}
        if verbose:
            print("[WARN] %s 不存在" % bench)
        return rep

    jobs = collect_nodes(bench, verbose=False)

    frames = {}                                             # variant -> {frame: record}
    for vdir in sorted([d for d in bench.iterdir() if d.is_dir()]):
        fdirs = sorted([f for f in vdir.glob("disp-*") if re.fullmatch(r"disp-\d+", f.name)])
        if not fdirs:
            continue
        frames[vdir.name] = {}
        for f in fdirs:
            oc = f / "OUTCAR"
            rec = {"frame": f.name, "n_atoms": None,
                   "max_abs_delta_meV_per_A": None, "rms_abs_delta_meV_per_A": None,
                   "force_rms_meV_per_A": None, "ref_force_rms_meV_per_A": None,
                   "rel_error_pct": None, "delta_atom_axis": None,
                   "elapsed_s": None, "cpu_time_s": None,
                   "electronic_steps": None, "electronic_steps_source": None,
                   "scf_converged": None, "scf_note": "",
                   "node": _frame_node(bench, vdir.name, f.name, jobs),
                   "jobid": _frame_jobid(jobs, vdir.name, f.name),
                   "stress_kbar": None,
                   "status": "PENDING", "note": "OUTCAR 缺失"}
            if oc.is_file():
                forces = parse_forces(oc)
                rec["n_atoms"] = len(forces) if forces else None
                rec["_forces"] = forces
                if forces:
                    frec = force_rms(forces)
                    rec["force_rms_meV_per_A"] = frec * 1000.0 if frec is not None else None
                rec["elapsed_s"] = parse_elapsed(oc)
                rec["cpu_time_s"] = parse_cpu_time(oc)
                rec["electronic_steps"], rec["electronic_steps_source"] = \
                    parse_electronic_steps(f)
                rec["scf_converged"], rec["scf_note"] = scf_converged(oc)
                rec["stress_kbar"] = last_stress_kbar(oc)
                rec["status"] = "OK" if forces else "PENDING"
                rec["note"] = "" if forces else "OUTCAR 里没有 TOTAL-FORCE"
            frames[vdir.name][f.name] = rec

    # 参考设置：显式 baseline 优先；否则按"精度最高"规则自动选
    metas = [_variant_dir_meta(bench / name, plan.get("variants") and
                               {v.get("name"): v for v in plan.get("variants", []) if v.get("name")} or {})
             for name in sorted(frames)]
    if baseline is None or baseline not in frames:
        with_forces = [m for m in metas
                       if any(rec.get("_forces") for rec in frames.get(m["name"], {}).values())]
        pool = with_forces or metas
        converged = [m for m in pool
                     if any(frames[m["name"]][fr].get("scf_converged") is True
                            for fr in frames.get(m["name"], {}))]
        baseline = (select_baseline(converged or pool) or baseline
                    or (sorted(frames)[0] if frames else None))

    base = frames.get(baseline) or {}
    have_base = any(r.get("_forces") for r in base.values())

    rows = []
    for name, fmap in frames.items():
        meta = next((m for m in metas if m.get("name") == name), {"name": name})
        for fname, rec in sorted(fmap.items()):
            br = base.get(fname)
            if name == baseline and rec.get("_forces"):
                rec["max_abs_delta_meV_per_A"] = 0.0
                rec["rms_abs_delta_meV_per_A"] = 0.0
                rec["rel_error_pct"] = 0.0
                rec["ref_force_rms_meV_per_A"] = rec.get("force_rms_meV_per_A")
                rec["status"] = "OK"
                rec["note"] = "参考设置"
            elif have_base and rec.get("_forces") and br and br.get("_forces"):
                d, wh = max_abs_delta(br["_forces"], rec["_forces"])
                rd = rms_delta(br["_forces"], rec["_forces"])
                ref = force_rms(br["_forces"])
                rec["max_abs_delta_meV_per_A"] = d
                rec["rms_abs_delta_meV_per_A"] = rd * 1000.0 if rd is not None else None
                rec["ref_force_rms_meV_per_A"] = ref * 1000.0 if ref is not None else None
                rec["rel_error_pct"] = (100.0 * rd / ref) if (rd is not None and ref) else None
                rec["delta_atom_axis"] = ("#%s-%s" % (wh[0] + 1, wh[1])) if wh else None
                if d is None:
                    rec["status"] = "PENDING"
                    rec["note"] = "原子数不一致 %s vs %s" % (
                        br.get("n_atoms"), rec.get("n_atoms"))
                else:
                    rec["status"] = _judge(rec, threshold_mev, threshold_rel_pct)
                    rec["note"] = ""
            elif rec.get("_forces"):
                # 该帧有力、但参考设置的同一帧缺力（如调研期间只跑完 frame-1）
                rec["status"] = "PENDING"
                rec["note"] = ("缺参考设置 %s 的力" % baseline) if not have_base                     else ("参考设置缺 %s 帧" % fname)
            rows.append((name, meta, rec))

    # ---- 汇总 ----
    variants = []
    for name in sorted(frames):
        meta = next((m for m in metas if m.get("name") == name), {"name": name})
        recs = [r for (n, m, r) in rows if n == name]
        deltas = [r["max_abs_delta_meV_per_A"] for r in recs
                  if r["max_abs_delta_meV_per_A"] is not None]
        rmses = [r["rms_abs_delta_meV_per_A"] for r in recs
                 if r["rms_abs_delta_meV_per_A"] is not None]
        rels = [r["rel_error_pct"] for r in recs if r["rel_error_pct"] is not None]
        frms = [r["force_rms_meV_per_A"] for r in recs
                if r["force_rms_meV_per_A"] is not None]
        elapsed = [r["elapsed_s"] for r in recs if r["elapsed_s"] is not None]
        cpus = [r["cpu_time_s"] for r in recs if r["cpu_time_s"] is not None]
        nodes = sorted({r["node"] for r in recs if r.get("node")})
        esteps = [r["electronic_steps"] for r in recs
                  if r.get("electronic_steps") is not None]
        convs = [r.get("scf_converged") for r in recs]
        v_conv = (all(c is True for c in convs)
                  if convs and all(c is not None for c in convs) else None)
        worst = max(deltas) if deltas else None
        statuses = {r["status"] for r in recs}
        any_notconv = any(c is False for c in convs)
        if "FAIL" in statuses:
            vstatus = "FAIL"
        elif "PENDING" in statuses or not recs:
            vstatus = "PENDING"
        elif any_notconv:
            vstatus = "NOTCONV"
        else:
            vstatus = "PASS"
        variants.append({
            "name": name, "LREAL": meta.get("LREAL"), "ADDGRID": meta.get("ADDGRID"),
            "ENCUT": meta.get("ENCUT"), "NELM": meta.get("NELM"),
            "ALGO": meta.get("ALGO"), "EDIFF": meta.get("EDIFF"),
            "desc": meta.get("desc", ""), "overrides": meta.get("overrides"),
            "duplicate_of": meta.get("duplicate_of"),
            "frames": [{k: v for k, v in r.items() if not k.startswith("_")} for r in recs],
            "worst_max_abs_delta_meV_per_A": worst,
            "worst_rms_abs_delta_meV_per_A": max(rmses) if rmses else None,
            "mean_force_rms_meV_per_A": (sum(frms) / len(frms)) if frms else None,
            "worst_rel_error_pct": max(rels) if rels else None,
            "mean_rel_error_pct": (sum(rels) / len(rels)) if rels else None,
            "mean_elapsed_s": (sum(elapsed) / len(elapsed)) if elapsed else None,
            "mean_cpu_time_s": (sum(cpus) / len(cpus)) if cpus else None,
            "electronic_steps": esteps,
            "scf_converged": v_conv,
            "nodes": nodes,
            "status": vstatus,
        })

    if not have_base:
        verdict = "PENDING"
    elif any(v["status"] == "FAIL" for v in variants):
        verdict = "FAIL"
    elif any(v["status"] == "NOTCONV" for v in variants):
        verdict = "NOTCONV"
    elif all(v["status"] == "PASS" for v in variants) and variants:
        verdict = "PASS"
    else:
        verdict = "PENDING"

    if threshold_rel_pct is not None:
        crit_desc = "相对 |ΔF|_rms/|F|_rms < %.2f%%" % float(threshold_rel_pct)
    else:
        crit_desc = "绝对 max|ΔF| < %.2f meV/Å" % float(threshold_mev or 0.0)

    # ---- 建议 ----
    if verdict == "PASS":
        passing = [v for v in variants if v["status"] == "PASS"
                   and v["mean_elapsed_s"] is not None]
        fastest = min(passing, key=lambda v: v["mean_elapsed_s"]) if passing else None
        rec = ("全部变体相对参考设置 %s —— 可自由选生产口径。" % crit_desc)
        if fastest is not None:
            rec += " 最快：%s（均 %.1f s/帧，worst rel=%.3f%%，max|ΔF|=%.3f meV/Å）；" % (
                fastest["name"], fastest["mean_elapsed_s"],
                fastest.get("worst_rel_error_pct") or 0.0,
                fastest["worst_max_abs_delta_meV_per_A"] or 0.0)
        rec += ("若不舍精度，建议维持 LREAL=.FALSE.+ADDGRID 的生产口径，"
                "除非速度收益明显。")
    elif verdict == "NOTCONV":
        bad = [v for v in variants if v["status"] == "NOTCONV"]
        rec = ("以下变体电子步未收敛（打满 NELM / 无 EDIFF 收敛标志），力可能不可靠，"
               "**不能**作为干净的 ENCUT/口径证据：%s。"
               "建议提高 NELM 或放宽 EDIFF（保持单一变量重跑）后再比较。"
               % "；".join("%s e-steps=%s" % (
                   v["name"], "/".join(str(x) for x in v.get("electronic_steps", [])) or "-")
                   for v in bad))
    elif verdict == "FAIL":
        bad = [v for v in variants if v["status"] == "FAIL"]
        detail = "；".join(
            "%s rel=%.3f%% RMS|ΔF|=%.3f max|ΔF|=%.3f meV/Å"
            % (v["name"], v.get("worst_rel_error_pct") or 0.0,
               v.get("worst_rms_abs_delta_meV_per_A") or 0.0,
               v["worst_max_abs_delta_meV_per_A"] or 0.0)
            for v in bad)
        rec = ("超阈值变体：%s。参考设置 %s 与它们不等价 —— **不要**把这些口径用于生产；"
               "生产 INCAR 维持能过阈值的组合（优先 LREAL=.FALSE.，"
               "ADDGRID 视 lreal_false_addgrid 行是否 PASS）。可增大超胞/收紧 EDIFF 后再测。"
               % (detail, baseline))
    else:
        missing = [v["name"] for v in variants if v["status"] == "PENDING"]
        rec = ("数据不完整：%s。先提交 bench 作业并等 OUTCAR 齐全后再 analyze。"
               % (", ".join(missing) if missing else "缺参考设置"))

    rep = {
        "skill": "kl-dft-cpu", "step": "step4_disp",
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "bench_dir": str(bench),
        "baseline": baseline,
        "baseline_rule": plan.get("baseline_rule",
                                  "LREAL=.FALSE.+ADDGRID=.TRUE. 中 ENCUT 最大"),
        "criterion_mode": ("relative" if threshold_rel_pct is not None else "absolute"),
        "threshold_rel_pct": (float(threshold_rel_pct)
                              if threshold_rel_pct is not None else None),
        "threshold_meV_per_A": (float(threshold_mev)
                                if threshold_mev is not None else None),
        "criterion": _criterion_text(threshold_mev, threshold_rel_pct),
        "prod_encut_eV": plan.get("prod_encut_eV"),
        "verdict": verdict, "recommendation": rec, "variants": variants,
    }
    (bench / REPORT_NAME).write_text(
        json.dumps(rep, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    if verbose:
        print_table(rep)
    return rep


def _disp_width(s):
    """显示宽度（中文/全角算 2 列），表头对齐用。"""
    return sum(2 if ord(c) > 0x2E80 else 1 for c in str(s))


def _pad(s, w, right=False):
    pad = " " * max(0, w - _disp_width(s))
    return (pad + str(s)) if right else (str(s) + pad)


def _fmt(x, pat="%.4f"):
    return "-" if x is None else (pat % x)


def print_table(rep):
    rel = rep.get("threshold_rel_pct")
    mev = rep.get("threshold_meV_per_A")
    crit = ("|ΔF|_rms/|F|_rms < %.2f%%" % float(rel)) if rel is not None \
        else ("max|ΔF| < %.2f meV/Å" % float(mev or 0.0))
    print("\n=== INCAR_BENCH 对照（参考=%s，判据 %s）==="
          % (rep["baseline"], crit))
    hdr = [("变体", 26), ("帧", 11), ("原子", 5), ("LREAL", 8), ("ADDGRID", 9),
           ("ENCUT", 6), ("max|ΔF|", 11), ("RMS|ΔF|", 11), ("|F|rms", 10),
           ("rel%", 8), ("CPU(s)", 10), ("Elapsed(s)", 10), ("node", 9),
           ("e-steps", 8), ("SCF", 8), ("判定", 8)]
    line = " ".join(_pad(h, w) for h, w in hdr)
    print(line)
    print("-" * _disp_width(line))
    for v in rep["variants"]:
        for r in v["frames"]:
            sig = r.get("stress_kbar")
            sig = max(abs(x) for x in sig) if sig else None
            d = r.get("max_abs_delta_meV_per_A")
            rd = r.get("rms_abs_delta_meV_per_A")
            fr = r.get("force_rms_meV_per_A")
            relr = r.get("rel_error_pct")
            print(" ".join([
                _pad(v["name"], 26),
                _pad(r["frame"], 11),
                _pad(r.get("n_atoms") if r.get("n_atoms") is not None else "-", 5),
                _pad(v.get("LREAL"), 8),
                _pad(v.get("ADDGRID") or "-", 9),
                _pad(v.get("ENCUT") if v.get("ENCUT") is not None else "-", 6),
                _pad(_fmt(d), 11),
                _pad(_fmt(rd), 11),
                _pad(_fmt(fr, "%.1f"), 10),
                _pad(_fmt(relr, "%.3f"), 8),
                _pad(_fmt(r.get("cpu_time_s"), "%.1f"), 10),
                _pad(_fmt(r.get("elapsed_s"), "%.1f"), 10),
                _pad(r.get("node") or "-", 9),
                _pad(r.get("electronic_steps")
                     if r.get("electronic_steps") is not None else "-", 8),
                _pad(("CONV" if r.get("scf_converged") is True
                      else ("NOTCONV" if r.get("scf_converged") is False else "-")), 8),
                _pad(r.get("status"), 8),
            ]))
    print("-" * _disp_width(line))
    for v in rep["variants"]:
        _esteps = "/".join(str(x) for x in v.get("electronic_steps", [])) or "-"
        _conv = ("CONV" if v.get("scf_converged") is True
                 else ("NOTCONV" if v.get("scf_converged") is False else "-"))
        print("  [%s] %-26s worst max|ΔF|=%-10s worst RMS|ΔF|=%-10s "
              "worst rel=%-9s mean CPU=%-8s e-steps=%-9s SCF=%-8s %s"
              % (v["status"], v["name"],
                 _fmt(v.get("worst_max_abs_delta_meV_per_A")),
                 _fmt(v.get("worst_rms_abs_delta_meV_per_A")),
                 (("%.3f%%" % v["worst_rel_error_pct"])
                  if v.get("worst_rel_error_pct") is not None else "-"),
                 _fmt(v.get("mean_cpu_time_s"), "%.1f"),
                 _esteps, _conv, v.get("desc", "")))
    print("\n结论：%s" % rep["verdict"])
    print("建议：%s" % rep["recommendation"])
    print("报告：%s/%s\n" % (rep["bench_dir"], REPORT_NAME))


# ==========================================================================
# submit：提交 bench 作业（独立命令，不进 autozt DAG）
# ==========================================================================
def submit(out, verbose=True, only=None, suffix=None):
    """对 bench 下每个还没 OUTCAR 的帧跑 sbatch submit.sh。返回 [(frame, jobid)]。

    only:   只提交这些变体目录名（逗号串或 list）；缺省全部。
    suffix: 只提交目录名以该后缀结尾的变体（如 "_ediiff1e7"）。
    提交后把 jobid 记进 incar_bench_jobs.json，并立即用 sacct/squeue 补节点名。
    ★ 这是【独立诊断 side-car】的提交口，故意不进 S4 的生产 fanout（否则会被
    autozt 当生产帧统计，违反隔离要求）。在 jzzn/A800 登录节点直接 sbatch 即可。
    """
    bench = Path(out) / BENCH_SUBDIR
    if not bench.is_dir():
        print("[ERROR] %s 不存在；先在 INCAR_BENCH=on 下跑 S4 gen" % bench)
        return []
    if isinstance(only, str):
        only = [x.strip() for x in only.split(",") if x.strip()]
    only_set = set(only) if only else None
    todo = []
    for vdir in sorted(d for d in bench.iterdir() if d.is_dir()):
        if only_set is not None and vdir.name not in only_set:
            continue
        if suffix and not vdir.name.endswith(suffix):
            continue
        for f in sorted(vdir.glob("disp-*")):
            if (f / "submit.sh").is_file() and not (f / "OUTCAR").is_file():
                todo.append(f)
    jobs = _load_json(bench / JOBS_NAME) or {}
    done = []
    for f in todo:
        p = subprocess.run(["sbatch", "submit.sh"], cwd=str(f),
                           capture_output=True, text=True)
        m = re.search(r"Submitted batch job\s+(\d+)", p.stdout or "")
        if p.returncode == 0 and m:
            rel = str(f.relative_to(bench))
            jid = m.group(1)
            d = jobs.setdefault(rel, {})
            d["jobid"] = jid
            d["job_name"] = "bench-%s" % rel.replace("/", "-")
            d["submitted_at"] = datetime.datetime.now().isoformat(timespec="seconds")
            done.append((rel, jid))
            if verbose:
                print("[OK] sbatch %s -> job %s" % (rel, jid))
        else:
            print("[FAIL] sbatch %s（rc=%d）：%s%s"
                  % (f, p.returncode, (p.stdout or "").strip(), (p.stderr or "").strip()))
    if done and bench.is_dir():
        (bench / JOBS_NAME).write_text(
            json.dumps(jobs, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
        collect_nodes(bench, verbose=verbose)
    if verbose:
        print("[..] 共提交 %d 个 bench 作业（总待提交 %d）" % (len(done), len(todo)))
    return done


# ==========================================================================
# CLI
# ==========================================================================
def _load_conf(material_dir):
    """只读 step.conf 的 ENCUT/ENCUT_FACTOR，不 import gen_step4_disp（避免拖 numpy）。"""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_common" / "opt"))
        import stepconf
        spec = {"ENCUT": (None, "int"), "ENCUT_FACTOR": (PROD_ENCUT_FACTOR, "float")}
        return stepconf.load(spec, None, cwd=str(material_dir), strict=False)
    except Exception as e:                                  # noqa: BLE001
        print("[WARN] 读不到 step.conf（%s），ENCUT 回落到 POTCAR/INCAR 推断" % e)
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="kl-dft-cpu S4 取力 INCAR 基准对照（隔离 side-car）")
    ap.add_argument("action", nargs="?", default="analyze",
                    choices=["generate", "submit", "nodes", "analyze"])
    ap.add_argument("--dir", default="step4_disp", help="step4_disp 目录（默认 cwd 下）")
    ap.add_argument("--frames", type=int, default=2, help="每个变体取生产前几帧（默认 2）")
    ap.add_argument("--force", action="store_true", help="generate 时强制重写已有帧")
    ap.add_argument("--baseline", default=None,
                    help="参考变体名；缺省=自动取精度最高者（LREAL=.FALSE.+ADDGRID+.最高 ENCUT）")
    ap.add_argument("--threshold-rel", type=float, default=None,
                    help="相对判据阈值（%%）；缺省 1.0（当未给 --threshold 时）")
    ap.add_argument("--threshold", type=float, default=None,
                    help="绝对判据阈值（meV/Å，可选；给了且未给 --threshold-rel 时用绝对判据）")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="generate 时对 INCAR 追加/覆盖标签（可重复，如 --set NELM=300）")
    ap.add_argument("--only", default=None,
                    help="generate 只生成这些变体（逗号分隔，如 encut_2.0x,encut_650）")
    ap.add_argument("--suffix", default="",
                    help="generate 的变体目录后缀（隔离命名，如 _nelm300）")
    a = ap.parse_args(argv)
    out = Path(a.dir)
    if a.action == "generate":
        conf = _load_conf(out.resolve().parent)
        overrides = {}
        for kv in a.set:
            if "=" in kv:
                k, v = kv.split("=", 1)
                v = v.strip()
                overrides[k.strip()] = int(v) if re.fullmatch(r"[-+]?\d+", v) else v
        plan = generate(out, conf=conf, n_frames=a.frames, force=a.force,
                        overrides=overrides, only=a.only, suffix=a.suffix)
        if plan is None:
            sys.exit(1)
    elif a.action == "submit":
        submit(out, only=a.only, suffix=(a.suffix or None))
    elif a.action == "nodes":
        collect_nodes(Path(out) / BENCH_SUBDIR)
    else:
        rel = a.threshold_rel
        mev = a.threshold
        if rel is None and mev is None:
            rel = DEFAULT_THRESHOLD_REL_PCT
        rep = analyze(out, baseline=a.baseline, threshold_mev=mev,
                      threshold_rel_pct=rel)
        if rep["verdict"] == "FAIL":
            sys.exit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
