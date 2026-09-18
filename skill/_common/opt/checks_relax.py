#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
checks_relax.py —— 公共技能池的私有判据：作业内分段弛豫（STAGE_MODE=in_job）

用法：技能的 skill.yaml 里写
    checks: ../_common/checks_relax.py
步骤里写
    {seq: 1, name: step1_PBE_opt, label: S1_opt, check: relax_injob, gen: gen_step1_PBE_opt.py}

tf 会把本文件的源码随 payload 下发到远端 exec，然后把 CHECKERS 里的判据注册进去
（见 tf 的"技能私有判据"一段）。所以本文件必须自包含：不要 import 技能里的其它模块。
可以直接用 payload 里已有的 tail_text/os/glob——下面做了兜底，缺了也能跑。

判据语义
--------
done  : OUTCAR 出现 "reached required accuracy"（收敛总闸）
        或者兼容老材料：同级 step1{a,b,c}_PBE_opt 里任意一段已收敛
fail  : 所有段都跑完了但没收敛 / 某段异常退出（.sN.started 有而 .sN.done 无）
running 之外的"还没跑" -> OUTCAR missing

注意：作业内分段只有一个目录，所以 tf 看不到段间状态；进度信息从
.sN.done 标记和 OUTCAR.sN 存档里读，判据的 note 会带上"已完成 N/M 段"。

.sN.done 只表示【这一段真的跑完了】。因力判据被提前跳过的段写 .sN.skipped
（见 relax_common.EARLY_EXIT：变胞段禁止被跳过），不计入已完成，note 里单独报。
"""

import os
import re
import sys
import glob as _glob


def _tail(path, n=200):
    fn = globals().get("tail_text")
    if fn:
        return fn(path)
    try:
        with open(path, "rb") as fh:
            try:
                fh.seek(-200000, os.SEEK_END)
            except OSError:
                fh.seek(0)
            return fh.read().decode("utf-8", "ignore")
    except OSError:
        return ""


def _conv(d):
    p = os.path.join(d, "OUTCAR")
    return os.path.isfile(p) and "reached required accuracy" in _tail(p)


def _stage_progress(d):
    """返回 (已完成段数, 计划段数, 有没有段跑了一半没完成, 被跳过的段数)。

    .sN.done = 该段真跑完；.sN.skipped = 该段因"上一段力已收敛"被跳过（只有
    不改晶胞的段会这样）。两者都算"占了一个位"，但只有 done 算完成。"""
    planned = len(_glob.glob(os.path.join(d, "INCAR.s*_*")))
    done = len(_glob.glob(os.path.join(d, ".s*.done")))
    started = len(_glob.glob(os.path.join(d, ".s*.started")))
    skipped = len(_glob.glob(os.path.join(d, ".s*.skipped")))
    return done, planned, started > done, skipped


# 末态压力门禁（kB）：脚本自身稳定判据是 PRESS_TOL_KB=1.0，这里放宽一倍留 Pulay 余量。
#   调大/调小请改这里，或在 step.conf 用同样的键覆盖（见下 STRESS_GATE_KB 读取）。
_STRESS_GATE_KB = 2.0
try:                                    # step.conf [params] 里可覆盖
    _sc_txt = open(os.path.join(os.getcwd(), "step.conf"), encoding="utf-8-sig").read()
    _m = re.search(r"(?m)^\s*STRESS_2D_THR\s*=\s*([\d.]+)", _sc_txt)
    if _m:
        _STRESS_GATE_KB = float(_m.group(1)) * 4.0   # 层内口径阈值 → 胞口径粗换算
except Exception:                       # noqa: BLE001 —— 读不到就用默认
    pass

_ISIF_RE = re.compile(r"^\s*ISIF\s*=\s*(\S+)", re.I | re.M)   # noqa: E305
_IOPT_RE = re.compile(r"^\s*IOPTCELL\s*=", re.I | re.M)


def _incar_moves_cell(text):
    """这一段会不会动晶胞（ISIF>=3，或写了 IOPTCELL）。与 relax_common.stage_changes_cell 同义。"""
    m = _ISIF_RE.search(text)
    if m and m.group(1)[:1].isdigit() and int(m.group(1)[:1]) >= 3:
        return True
    return bool(_IOPT_RE.search(text))


def _cell_stage_states(d):
    """-> [(段文件名, 是否动胞, 状态)]，状态 ∈ {done, skipped, started, ''}。"""
    out = []
    for p in sorted(_glob.glob(os.path.join(d, "INCAR.s*_*"))):
        tag = os.path.basename(p)[len("INCAR."):]          # 形如 s2_b
        key = tag.split("_")[0]                            # s2
        try:
            moves = _incar_moves_cell(open(p, encoding="utf-8", errors="ignore").read())
        except OSError:
            moves = False
        st = ""
        for s in ("done", "skipped", "started"):
            if os.path.isfile(os.path.join(d, ".%s.%s" % (key, s))):
                st = s
                break
        out.append((tag, moves, st))
    return out


_PRESS_RE = re.compile(r"external pressure\s*=\s*(-?[\d.]+)\s*kB")


def _last_external_pressure(d):
    """末态 external pressure (kB)；读不到返回 None。

    自包含实现（不 import 技能模块）：只认 OUTCAR 里最后一条
    "external pressure = X kB"。它是应力张量迹的 1/3，与"晶胞是否到位"直接对应。
    """
    p = os.path.join(d, "OUTCAR")
    if not os.path.isfile(p):
        return None
    hits = _PRESS_RE.findall(_tail(p, 400))
    try:
        return float(hits[-1]) if hits else None
    except (TypeError, ValueError):
        return None


_ZBRENT_RE = re.compile(r"ZBRENT:\s*fatal error in bracketing"
                        r"|I REFUSE TO CONTINUE WITH THIS SICK JOB")


def _stage_crash_kind(d):
    """这一段是【怎么】中断的 —— 'zbrent' / 'watchdog' / ''。

    为什么值得单独判：旧文案一律说"可能撞墙钟或被看门狗杀掉"。实测（Ti2S3 Z4-3-1，
    jzzn jobid 3847708，2026-09-17）是 VASP 自己 ZBRENT 崩的，与墙钟/看门狗完全无关；
    照旧文案去查墙钟和 queue.err 会查错方向，白耗一轮。只读各段存档 OUTCAR* 与 queue.err。
    """
    for p in sorted(_glob.glob(os.path.join(d, "OUTCAR*"))):
        try:
            if _ZBRENT_RE.search(_tail(p, 200)):
                return "zbrent"
        except OSError:
            pass
    p = os.path.join(d, "queue.err")
    if os.path.isfile(p):
        try:
            if "[watchdog]" in open(p, encoding="utf-8", errors="ignore").read():
                return "watchdog"
        except OSError:
            pass
    return ""


_STRESS_GATE_MARKERS = ("S1_STRESS_GATE.off", "S1_STRESS_GATE_OFF")


def _stress_gate_exempt(d):
    """应力门禁豁免：步骤目录里放 S1_STRESS_GATE.off（内容写原因）即跳过末态压力判据。

    为什么需要：末态压力判据针对的是【2D 生产材料】——面内张力会把 ZA 线性化、
    Huang 零应力条件不成立。但基准/验证类结构（Si 金刚石、LiCoO2 等 shengbte-gpu
    对照线）本来就是拿官方 example 的结构与力常数做对照，冻结晶格很可能是有意为之，
    重新弛豫反而会让它们与官方基准对不上。豁免必须写明原因，留痕可审计。
    """
    for name in _STRESS_GATE_MARKERS:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            try:
                txt = open(p, encoding="utf-8", errors="ignore").read().strip()
            except OSError:
                txt = ""
            return txt.splitlines()[0][:160] if txt else "（未写原因）"
    return None


def _job_reported_unconverged_cell(d):
    """作业自报变胞段未稳定（run_relax.sh 末尾的 RELAX_CELL_UNCONVERGED）。"""
    for f in ("queue.out", "queue.err"):
        p = os.path.join(d, f)
        if not os.path.isfile(p):
            continue
        try:
            if "RELAX_CELL_UNCONVERGED" in open(p, encoding="utf-8", errors="ignore").read():
                return True
        except OSError:
            pass
    return False


def ck_relax_injob(d, sc):
    """作业内分段弛豫的判据（check: relax_injob）。"""
    mat = os.path.dirname(d)

    # 老材料兼容：tf 分段时代留下的 step1a/b/c，任意一段收敛即认账
    for legacy in ("step1a_PBE_opt", "step1b_PBE_opt", "step1c_PBE_opt",
                   "step1a_std_opt", "step1b_std_opt", "step1c_std_opt"):
        if _conv(os.path.join(mat, legacy)):
            return True, "旧分段 %s 已收敛，跳过" % legacy

    if _conv(d):
        # ★ 2026-09-16 加固（这是"前面遇到过、后面不许再遇到"的那一类）：
        #   ① 变胞段被"上一段力已收敛"跳过 —— 力判据说明不了晶胞/应力。实测 jzz 批次：
        #      段2(b)/段3(c) 被跳过并写成 .sN.done，晶格 16 位不变、面内层内口径应力
        #      −2~−62 kbar（全拉伸），而判据报 "converged（3/3 段）"，kappa 照跑。
        #   ② 作业自己报了 RELAX_CELL_UNCONVERGED（力收敛但晶格变化/|P| 没到位）——
        #      以前这个标记写了没人读，等于白写。
        #   两条都在"看起来收敛"之前拦，绝不让带残余应力的结构流到 S2/S4。
        bad = [t for t, moves, st in _cell_stage_states(d) if moves and st == "skipped"]
        if bad:
            return False, ("变胞段 %s 被跳过 —— 力收敛不代表晶胞/应力收敛，晶胞没有弛豫，"
                           "残余应力会原样传到声子/热导步。清掉 .s?.done/.s?.skipped 后 "
                           "autozt retry 重投（新代码的 gen 会自动清）。" % "、".join(bad))
        if _job_reported_unconverged_cell(d):
            return False, ("作业自报 RELAX_CELL_UNCONVERGED —— 变胞段的晶格变化量或 |external "
                           "pressure| 没到位（力判据过了也没用）。看 queue.out 末态 P 值，"
                           "按提示 cp CONTCAR POSCAR 后重投。")
        # ③ 末态压力后置校验（物理量判据，不依赖标记语义）：有变胞段却没把应力收到位，
        #    说明晶胞没真弛豫。实测故障件：段 b/c 被跳过 → 晶格 16 位不变、末态
        #    external pressure = -7.93 kB，而旧判据报 "converged（3/3 段）"。
        #    阈值取脚本自身稳定判据 PRESS_TOL_KB=1.0 的两倍，给 Pulay 应力留余量。
        _exempt = _stress_gate_exempt(d)
        if _exempt is not None:
            # ★ 必须写 stderr：本判据的 stdout 会被上层采集器当作【一整份 JSON】解析，
            #   往 stdout 打任何东西都会让 json.loads 从第 1 个字符就失败
            #   （2026-09-18 实测：jzzn 的 kl-dft-cpu 组因此整组"无材料"，
            #    因为同目录下的豁免材料把提示打进了采集器 stdout）。
            print("[..] 变胞应力门禁已豁免（S1_STRESS_GATE.off：%s）" % _exempt,
                  file=sys.stderr)
        if _exempt is None and any(m for _, m, _ in _cell_stage_states(d)):
            _p = _last_external_pressure(d)
            if _p is not None and abs(_p) > _STRESS_GATE_KB:
                return False, ("有变胞段但末态 external pressure = %.2f kB（阈值 %.1f）—— "
                               "晶胞没有弛豫到位（力收敛说明不了应力收敛）。这正是"
                               "'变胞段被跳过'故障的指纹：晶格不变 + 残余应力留下。"
                               "autozt retry 重投（gen 会清掉旧阶段标记，变胞段不会再被跳过）。"
                               % (_p, _STRESS_GATE_KB))
        done, planned, _, skipped = _stage_progress(d)
        if planned:
            # 变胞多遍循环会让同一个 INCAR 被跑多遍（.s2r1/.s2r2/…），于是
            # done > planned。原来直接印 "converged（5/3 段）"，看着像 bug
            # （2026-09-17 实测就是这个），所以分开表述。
            if done == planned:
                note = "converged（%d/%d 段" % (done, planned)
            else:
                note = "converged（%d 个阶段 INCAR、共 %d 次段运行" % (planned, done)
            if skipped:
                note += "，跳过 %d 段" % skipped
            return True, note + "）"
        return True, "converged"

    if not os.path.isfile(os.path.join(d, "OUTCAR")):
        return False, "OUTCAR missing"

    done, planned, half, skipped = _stage_progress(d)
    if planned and done >= planned:
        return False, ("%d 段全部跑完但未收敛 —— 看 OUTCAR.s* / OSZICAR.s* 定位是哪一段"
                       "开始震荡，调完 step.conf 后 autozt retry" % planned)
    if half:
        # 是哪一段没跑完（.sN.started 且无 .sN.done），报名字比报序号有用
        _half_tags = [os.path.basename(t)[len("INCAR."):].split("_")[0]
                      for t, _, st in _cell_stage_states(d) if st == "started"]
        _who = "、".join(_half_tags) if _half_tags else "第 %d 段" % (done + 1)
        kind = _stage_crash_kind(d)
        if kind == "zbrent":
            return False, ("%s 段中断（有 .started 无 .done）：VASP 报 ZBRENT: fatal error in "
                           "bracketing —— 变胞段的 CG 线搜索夹不到极小值，最常见的原因是"
                           "【该段起点晶胞已接近平衡、没有梯度可走】（例如接手了上一次作业"
                           "遗留的旧 CONTCAR）。旧代码会把这种上报成\"可能撞墙钟或被看门狗"
                           "杀掉\"，方向是错的。重投即可：gen 会清掉旧阶段标记、从当前 POSCAR "
                           "重跑整条段序（2026-09-17 起 run_relax 还有 STAGES_RUN 闸门，"
                           "不会再用上次遗留的 OUTCAR 把第一段跳过）；若反复出现再考虑"
                           "收紧 EDIFF 或减小 POTIM。" % _who)
        if kind == "watchdog":
            return False, ("%s 段中断（有 .started 无 .done）：被 run_relax 看门狗判卡死杀掉"
                           "（queue.err 里有 [watchdog] 判定快照），重投会从 CONTCAR 续跑"
                           % _who)
        return False, ("%s 段中断（有 .started 无 .done）：可能撞墙钟或被看门狗杀掉；"
                       "重投会从 CONTCAR 续跑" % _who)
    note = "已完成 %d/%d 段" % (done, planned or 1)
    if skipped:
        note += "（另有 %d 段因力判据被跳过）" % skipped
    return False, note


CHECKERS = {"relax_injob": ck_relax_injob}
