# -*- coding: utf-8 -*-
"""progress —— 机器可读进度文件 .tf_progress.json + autozt progress 命令。

为什么要它：2026-09 的 295 材料实测里，所有麻烦都集中在「拿状态」——summary 要
ssh 采集、和 monitor 抢 ssh/9p、1500s 跑不完；只有 count 没有逐材料状态；想看
单个项目只能 -p 列 295 个名字。

做法：每次**真正采集**（monitor 每轮、以及 list/summary/status 等不走缓存的命令）
之后，把结果压成一份小 JSON 落在配置目录。AI / 脚本 / 人直接读文件：秒级、
不连超算、不受 9p 影响，也不会和 monitor 抢 ssh。

文件结构（schema = autozt-progress/1，只追加字段、不改含义）：
  ts / time           写入时刻（epoch / ISO）
  writer              monitor | cli
  monitor             {pid, interval, round}（monitor 写入时）
  round               本轮摘要：提交/新完成/新失败数、耗时（monitor 写入时）
  counts              全部步骤按状态计数（OK/R/PD/FAIL/TODO/PREP/WAIT/…）
  projects            每个项目（tf_<项目>.yaml）× 技能的材料级计数
  materials[]         每材料一行：id（<项目名>/<完整名>，稳定主键）、project、tt、
                      state（done 或活动步骤的状态）、active、steps{label: 状态}、
                      job（活动步骤作业）、fails[]（code/class/action/fail_count）、
                      since（当前状态从何时起）、eta_s（按历史同类步骤中位耗时估计）、
                      busy（在跑/排队作业数，扇出按子作业计；与 max_jobs 闸门同口径）
  coverage            {技能: 上次「整技能」采集的时刻}——只有这里够新，单材料命令才敢
                      只采目标材料、其余材料的并发占用从本文件 + 提交账本补（见 busy_baseline）

范围：采集只覆盖一部分（-p/-tt/--project）时按材料 id 合并进已有文件，不会把
monitor 写的全量视图截断成局部视图。
"""

import datetime
import json
import os
import sys
import time

PROGRESS_NAME = ".tf_progress.json"
SCHEMA = "autozt-progress/1"
_ETA_MIN_SAMPLES = 3
# 提交账本：每次成功 sbatch 追加一行。进度文件只在「采集」后写，提交之后到下一次采集
# 之间它看不到新作业；只采目标材料的命令靠账本把这段空档补上，max_jobs 不会被低估。
LEDGER_NAME = ".tf_submit_ledger.jsonl"
_LEDGER_KEEP_S = 3 * 86400
_LEDGER_PRUNE_BYTES = 1 << 20


def progress_path(cfg):
    return os.path.join((cfg or {}).get("_config_dir") or os.getcwd(), PROGRESS_NAME)


def _iso(ts):
    return datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%dT%H:%M:%S")


def _parse_iso(text):
    try:
        return time.mktime(datetime.datetime.strptime(
            str(text)[:19], "%Y-%m-%dT%H:%M:%S").timetuple())
    except (TypeError, ValueError):
        return None


def load_progress(cfg):
    try:
        with open(progress_path(cfg), encoding="utf-8") as f:
            d = json.load(f)
        return d if isinstance(d, dict) and d.get("schema") == SCHEMA else None
    except (OSError, ValueError):
        return None


def load_progress_view(cfg):
    """给人/AI 看的进度：原文件 + 采集之后的提交（账本）叠加成 PD。

    进度文件只在采集后写；monitor 更是先采集、再推进提交、最后才写文件——刚提交的
    作业在文件里还是 TODO，AI 下一轮 inspect 会把它们再提议一遍。这里读时把账本里
    比该材料观测时刻更新的提交标成 PD（submitted_after_collect），文件本身不改。"""
    doc = load_progress(cfg)
    if doc is None:
        return None
    mats = doc.get("materials") or []
    since = min([float(e.get("obs") or 0) for e in mats] or [0])
    led = ledger_recent(cfg, since) if mats else []
    if not led:
        return doc
    by_id = {}
    for r in led:
        by_id.setdefault((r.get("tt"), r.get("id")), []).append(r)
    out = []
    for e in mats:
        recs = [r for r in by_id.get((e.get("tt"), e.get("id")), ())
                if float(r.get("ts") or 0) > float(e.get("obs") or 0)]
        if not recs:
            out.append(e)
            continue
        e = dict(e)
        e["steps"] = dict(e.get("steps") or {})
        for r in recs:
            lab = str(r.get("step"))
            if lab in e["steps"]:
                e["steps"][lab] = "PD"
            if lab == str(e.get("active")):
                e["state"] = "PD"
                e["job"] = {"id": r.get("jid"), "state": "SUBMITTED"}
            e["busy"] = int(e.get("busy") or 0) + int(r.get("jobs") or 1)
            e["fails"] = [f for f in e.get("fails") or [] if str(f.get("step")) != lab]
        if not e["fails"]:
            e.pop("fails", None)
        e["submitted_after_collect"] = True
        out.append(e)
    doc = dict(doc)
    doc["materials"] = out
    return doc


def _history_stats(cfg):
    """从 history.jsonl 读出 (各 技能|步骤 的完成耗时列表, 各 材料|技能|步骤 的失败次数)。"""
    from autozt.history import history_path
    durs, fails = {}, {}
    path = history_path(cfg)
    if not os.path.isfile(path):
        return durs, fails
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if '"FAIL"' not in line and '"finish"' not in line:
                    continue
                try:
                    e = json.loads(line)
                except ValueError:
                    continue
                if e.get("ev") == "finish" and e.get("dur"):
                    durs.setdefault((e.get("skill"), e.get("step")), []).append(
                        int(e["dur"]))
                elif e.get("t") == "FAIL" and e.get("f") != "FAIL" \
                        and e.get("ev") in ("state", "fail"):
                    k = (e.get("mat"), e.get("skill"), e.get("step"))
                    fails[k] = fails.get(k, 0) + 1
    except OSError:
        pass
    return durs, fails


def _history_started(cfg):
    """history 状态快照里记下的「第一次看到 R/PD」时刻：材料\\t技能\\t步骤名 -> epoch。"""
    from autozt.history import history_state_path
    try:
        with open(history_state_path(cfg), encoding="utf-8") as f:
            st = json.load(f)
    except (OSError, ValueError):
        return {}
    out = {}
    for k, v in (st or {}).items():
        ts = _parse_iso((v or {}).get("s")) if isinstance(v, dict) else None
        if ts:
            out[k] = ts
    return out


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def _material_state(m):
    kinds = [s.get("kind") for s in m.get("steps") or []]
    if kinds and all(k == "OK" for k in kinds):
        return "done"
    a = m.get("active")
    if isinstance(a, dict):
        return a.get("kind") or "?"
    return "?"


def build_entries(cfg, data, now=None):
    """collect_data 的结果 → {材料 id: 进度条目}。"""
    from autozt import _diag_class, _diag_code, _material_busy_jobs, _suggested_action
    now = now or time.time()
    durs, fails = _history_stats(cfg)
    started = _history_started(cfg)
    try:
        max_retries = int((cfg or {}).get("max_auto_retries", 2))
    except (TypeError, ValueError):
        max_retries = 2
    out = {}
    for t in (data or {}).get("types") or []:
        key = t.get("key")
        for m in t.get("materials") or []:
            mid = m.get("qualified_name") or m.get("name")
            a = m.get("active") if isinstance(m.get("active"), dict) else None
            steps = m.get("steps") or []
            e = {"id": mid, "name": m.get("name"), "project": m.get("project") or "",
                 "tt": key, "hpc": m.get("hpc_name") or "",
                 "state": _material_state(m),
                 "active": (a or {}).get("label"),
                 "action": m.get("action") or "",
                 "done_steps": sum(1 for s in steps if s.get("kind") == "OK"),
                 "total_steps": len(steps),
                 "busy": _material_busy_jobs(m),
                 "steps": {str(s.get("label")): s.get("kind") for s in steps}}
            j = (a or {}).get("job") or {}
            if j:
                e["job"] = {k: j.get(k) for k in ("id", "state", "info", "time")
                            if j.get(k) not in (None, "")}
            fl = []
            for s in steps:
                if s.get("kind") != "FAIL":
                    continue
                code = _diag_code(s.get("diag") or "")
                act, reason = _suggested_action(code)
                n = fails.get((m.get("name"), key, s.get("label") or s.get("name")), 0)
                item = {"step": s.get("label"), "code": code, "class": _diag_class(code),
                        "action": act, "reason": reason, "fail_count": n,
                        "diag": (s.get("diag") or "")[:200]}
                if act == "retry" and n > max_retries:
                    # AGENTS.md 决策表：同一步骤 retry 过 2 次仍 FAIL 就停手找人。
                    item["action"] = "human_review"
                    item["reason"] = ("同一步骤已失败 %d 次（>max_auto_retries=%d），"
                                      "停止自动重投，需人工介入" % (n, max_retries))
                fl.append(item)
            if fl:
                e["fails"] = fl
            if a and a.get("kind") in ("R", "PD"):
                sk = "%s\t%s\t%s" % (m.get("name"), key, a.get("name"))
                t0 = started.get(sk)
                if t0:
                    e["since"] = _iso(t0)
                    samples = durs.get((key, a.get("label") or a.get("name"))) or []
                    if len(samples) >= _ETA_MIN_SAMPLES:
                        e["eta_s"] = max(0, int(_median(samples) - (now - t0)))
                        e["eta_basis"] = "median of %d finished %s runs" % (
                            len(samples), a.get("label"))
            out[mid] = e
    return out


def _rollup(entries):
    counts, projects, mstates = {}, {}, {}
    for e in entries.values():
        for k in (e.get("steps") or {}).values():
            counts[k] = counts.get(k, 0) + 1
        st = e.get("state") or "?"
        mstates[st] = mstates.get(st, 0) + 1
        pj = projects.setdefault(e.get("project") or "(global)", {})
        row = pj.setdefault(e.get("tt") or "?", {"materials": 0, "states": {}})
        row["materials"] += 1
        row["states"][st] = row["states"].get(st, 0) + 1
    return counts, mstates, projects


def write_progress(cfg, data, writer="cli", full_scope=False, round_info=None,
                   monitor=None, covered_keys=None, observed=None):
    """把本次采集写进 .tf_progress.json（原子替换）。失败只告警，绝不影响主流程。

    covered_keys：本次采集完整覆盖了哪些技能（没按 -p/--project/目录裁剪）；记进
    coverage，供 busy_baseline 判断本文件的并发占用数据是否可信。
    observed：采集**开始**时刻（monitor 在推进/提交之后才写文件，必须用采集时刻，
    否则本轮刚提交的作业会被当成「文件里已经看到」而漏算）。"""
    try:
        now = time.time()
        _col = sys.modules.get("autozt.collect")
        if (getattr(_col, "COLLECT_TIMING", None) or {}).get("failed_groups"):
            # 有采集组失败被跳过：这轮缺材料，按局部合并，也不标覆盖
            full_scope, covered_keys = False, None
        obs = float(observed or now)
        fresh = build_entries(cfg, data, now)
        prev = load_progress(cfg) or {}
        old = {e.get("id"): e for e in prev.get("materials") or [] if e.get("id")}
        merged = {} if full_scope else dict(old)
        for mid, e in fresh.items():
            o = old.get(mid)
            fp = (e.get("state"), e.get("active"), json.dumps(e.get("steps"), sort_keys=True))
            ofp = ((o.get("state"), o.get("active"), json.dumps(o.get("steps"), sort_keys=True))
                   if o else None)
            if "since" not in e:
                # 没有作业起始时刻（非 R/PD）：状态没变就沿用上次的 since
                e["since"] = (o.get("since") if o and ofp == fp and o.get("since")
                              else _iso(now))
            e["updated"] = _iso(now)
            e["obs"] = round(obs, 3)
            merged[mid] = e
        mats = sorted(merged.values(), key=lambda e: (e.get("project") or "",
                                                      e.get("tt") or "", e.get("id") or ""))
        counts, mstates, projects = _rollup({e["id"]: e for e in mats})
        out = {"schema": SCHEMA, "ts": int(now), "time": _iso(now), "writer": writer,
               "full_scope": bool(full_scope or prev.get("full_scope")),
               "counts": counts, "material_states": mstates, "projects": projects,
               "queue": (data or {}).get("queue") or prev.get("queue") or {},
               "coverage": dict(prev.get("coverage") or {}),
               "materials": mats}
        if full_scope:
            covered_keys = [t.get("key") for t in (data or {}).get("types") or []]
        for key in covered_keys or ():
            if key:
                out["coverage"][key] = round(obs, 3)
        if monitor:
            out["monitor"] = monitor
        elif prev.get("monitor"):
            out["monitor"] = prev["monitor"]
        if round_info:
            out["round"] = round_info
        elif prev.get("round"):
            out["round"] = prev["round"]
        path = progress_path(cfg)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        tmp = path + ".tmp.%d" % os.getpid()
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(out, f, ensure_ascii=False, separators=(",", ":"))
        os.replace(tmp, path)
        return path
    except Exception as exc:                      # noqa: BLE001
        print("警告：进度文件写入失败（忽略）：%s" % exc, file=sys.stderr)
        return None


# =============================================================================
# 提交账本 + 只采目标材料时的并发基线
# =============================================================================
def ledger_path(cfg):
    return os.path.join((cfg or {}).get("_config_dir") or os.getcwd(), LEDGER_NAME)


def ledger_append(cfg, m, s, jid):
    """成功 sbatch 后追加一行（O_APPEND 小行写，多进程安全）。失败静默：账本只用于
    补并发计数，写不进去时只采目标材料的命令会因为没有账本而更保守（见 busy_baseline）。
    没有配置目录（单元测试里的裸 cfg）就不写——绝不往当前目录乱落文件。"""
    if not (cfg or {}).get("_config_dir"):
        return
    try:
        n = len([x for x in str(jid or "").split(",") if x.strip()]) or 1
        rec = {"ts": round(time.time(), 3), "tt": m.get("tt"),
               "id": m.get("qualified_name") or m.get("name"),
               "project": m.get("project") or "",
               "from": ((m.get("_seg") or {}).get("_from")) or "",
               "step": s.get("label"), "jobs": n, "jid": str(jid or "")}
        path = ledger_path(cfg)
        with open(path, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False, sort_keys=True) + "\n")
        if os.path.getsize(path) > _LEDGER_PRUNE_BYTES:
            keep = ledger_recent(cfg, time.time() - _LEDGER_KEEP_S)
            tmp = path + ".tmp.%d" % os.getpid()
            with open(tmp, "w", encoding="utf-8") as f:
                for r in keep:
                    f.write(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n")
            os.replace(tmp, path)
    except Exception:                             # noqa: BLE001
        pass


def ledger_recent(cfg, since_ts):
    out = []
    try:
        with open(ledger_path(cfg), encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                if isinstance(r, dict) and float(r.get("ts") or 0) >= since_ts:
                    out.append(r)
    except OSError:
        pass
    return out


def agent_event(kind, **fields):
    """给 agent 执行器看的机器可读行（仅 AUTOZT_AGENT_EVENTS=1 时输出，人看的输出不变）。"""
    if not os.environ.get("AUTOZT_AGENT_EVENTS"):
        return
    fields["event"] = kind
    print("[agent-event] " + json.dumps(fields, ensure_ascii=False, sort_keys=True),
          flush=True)


def narrow_max_age(cfg):
    """进度文件里「整技能覆盖」多久以内算新（秒）：AUTOZT_NARROW_MAX_AGE > tf.yaml
    narrow_max_age > 默认 7200；0 = 关闭只采目标材料（一律全量采集）。"""
    from autozt import tuning_value
    return tuning_value("narrow_max_age", cfg or {})


def narrow_keys(cfg, keys, now=None):
    """哪些技能可以只采目标材料：进度文件在 narrow_max_age 内对该技能做过整技能采集。
    返回 (可裁剪的技能集合, 进度文档, {技能: 不能裁剪的原因})。"""
    now = now or time.time()
    limit = narrow_max_age(cfg)
    doc = load_progress(cfg) if limit > 0 else None
    cov = (doc or {}).get("coverage") or {}
    ok, why = set(), {}
    for key in keys:
        ts = cov.get(key)
        if limit <= 0:
            why[key] = "narrow_max_age=0"
        elif doc is None:
            why[key] = "no progress file"
        elif not ts:
            why[key] = "no full-skill collection recorded"
        elif now - float(ts) > limit:
            why[key] = "last full-skill collection %ds ago > %ds" % (now - float(ts), limit)
        else:
            ok.add(key)
    return ok, doc, why


def busy_baseline(cfg, doc, data, keys, from_by_project):
    """只采了目标材料时，没采到的同技能材料的并发占用（给 _SkillGate 加底数）。

    来源：进度文件每材料的 busy（采集时刻的在跑/排队作业数）+ 账本里该材料在那次
    采集之后的提交。已完成但文件里还是 R/PD 的会被多算——宁多勿少（多算只会让提交
    等一等，少算会超 max_jobs）。返回 {"skill": {key: n}, "project": {"key\tfrom": n}}。"""
    present = set()
    for t in (data or {}).get("types") or []:
        for m in t.get("materials") or []:
            present.add((t.get("key"), m.get("qualified_name") or m.get("name")))
    cov = (doc or {}).get("coverage") or {}
    oldest = min([float(cov.get(k) or 0) for k in keys] or [0])
    led = ledger_recent(cfg, oldest)
    skill, proj = {}, {}

    def _add(key, pfrom, n):
        if n <= 0:
            return
        skill[key] = skill.get(key, 0) + n
        if pfrom:
            pk = "%s\t%s" % (key, pfrom)
            proj[pk] = proj.get(pk, 0) + n

    seen = set()
    for e in (doc or {}).get("materials") or []:
        key, mid = e.get("tt"), e.get("id")
        if key not in keys or (key, mid) in present:
            continue
        seen.add((key, mid))
        n = e.get("busy")
        if n is None:   # 旧文件没有 busy：按步骤状态数（扇出会少算子作业，仅兜底）
            n = sum(1 for k in (e.get("steps") or {}).values() if k in ("R", "PD"))
        obs = float(e.get("obs") or _parse_iso(e.get("updated")) or 0)
        n += sum(int(r.get("jobs") or 1) for r in led
                 if r.get("tt") == key and r.get("id") == mid
                 and float(r.get("ts") or 0) > obs)
        _add(key, from_by_project.get(e.get("project") or ""), int(n))
    for r in led:   # 进度文件里还没有的材料：覆盖之后的提交全算
        key, mid = r.get("tt"), r.get("id")
        if key not in keys or (key, mid) in present or (key, mid) in seen:
            continue
        if float(r.get("ts") or 0) > float(cov.get(key) or 0):
            _add(key, r.get("from") or from_by_project.get(r.get("project") or ""),
                 int(r.get("jobs") or 1))
    return {"skill": skill, "project": proj}


# =============================================================================
# 读取（autozt progress）
# =============================================================================
def filter_progress(doc, project=None, tt=None, material=None, status=None):
    """按项目/技能/材料/状态过滤进度文档，返回新文档（rollup 重算）。"""
    from autozt.agent_protocol import material_matches, split_scope, status_kinds
    projs = set(split_scope(project))
    tts = set(split_scope(tt))
    mats = split_scope(material)
    kinds = status_kinds(status)
    keep = []
    for e in (doc or {}).get("materials") or []:
        if projs and (e.get("project") or "") not in projs:
            continue
        if tts and e.get("tt") not in tts:
            continue
        if mats and not material_matches(e.get("name"), mats, e.get("id")):
            continue
        if kinds and not any(k in kinds for k in (e.get("steps") or {}).values()):
            continue
        keep.append(e)
    out = dict(doc or {})
    counts, mstates, projects = _rollup({e["id"]: e for e in keep})
    out.update({"materials": keep, "counts": counts, "material_states": mstates,
                "projects": projects, "fail_summary": fail_summary(keep)})
    return out


def fail_summary(entries, max_ids=30):
    """FAIL 按建议动作分桶：{action: {count, codes{code: n}, ids[…]}}。
    给 AI 的分诊入口：retry 桶可提议重投，human_review 桶直接汇报给人。"""
    out = {}
    for e in entries or []:
        for f in e.get("fails") or []:
            b = out.setdefault(f.get("action") or "human_review",
                               {"count": 0, "codes": {}, "ids": []})
            b["count"] += 1
            code = f.get("code") or "unknown"
            b["codes"][code] = b["codes"].get(code, 0) + 1
            if len(b["ids"]) < max_ids:
                b["ids"].append("%s#%s" % (e.get("id"), f.get("step")))
    return out


def monitor_liveness(cfg, doc):
    """后台 monitor 是否还活着、进度文件有多旧（读时现算，不依赖文件内容）。"""
    from autozt import _watch_running_pid
    pid, _pf = _watch_running_pid(cfg)
    age = int(time.time() - float((doc or {}).get("ts") or 0)) if doc else None
    interval = ((doc or {}).get("monitor") or {}).get("interval")
    # 一轮本身可能比间隔长得多（9p 上 300+ 材料一轮 ~50 分钟）：按「间隔 + 上一轮耗时」算
    # 周期，超过 3 个周期没更新才算过期——只按间隔算会在健康的慢 monitor 上误报。
    try:
        dur = float(((doc or {}).get("round") or {}).get("dur_s") or 0)
    except (TypeError, ValueError):
        dur = 0.0
    period = int(interval) + dur if interval else 0
    stale = bool(doc) and age is not None and bool(period) and age > 3 * period
    return {"alive": bool(pid), "pid": pid, "age_s": age,
            "stale": bool(stale) if doc else True}


def cmd_progress(cfg, project=None, tt=None, material=None, status=None,
                 json_out=False, limit=None):
    """autozt progress [--project P] [-tt T] [-p M] [-status S] [--json]

    只读本地进度文件：不采集、不连超算、不扫盘。"""
    doc = load_progress_view(cfg)
    if doc is None:
        msg = ("还没有进度文件 %s。它由 monitor 每轮、或任何一次真正采集"
               "（autozt list --refresh / summary --refresh / status）自动写出。"
               % progress_path(cfg))
        if json_out:
            print(json.dumps({"ok": False, "error": msg, "path": progress_path(cfg)},
                             ensure_ascii=False))
        else:
            print(msg)
        return 1
    view = filter_progress(doc, project, tt, material, status)
    view["path"] = progress_path(cfg)
    view["liveness"] = monitor_liveness(cfg, doc)
    if json_out:
        if limit:
            view["materials_total"] = len(view["materials"])
            view["materials"] = view["materials"][:int(limit)]
        print(json.dumps(view, ensure_ascii=False, separators=(",", ":")))
        return 0
    live = view["liveness"]
    print("进度  %s（%s 写于 %s，%s 秒前）  monitor：%s"
          % (view["path"], view.get("writer"), view.get("time"), live.get("age_s"),
             ("在跑 PID %s" % live["pid"]) if live.get("alive") else "未运行"))
    if live.get("stale"):
        print("  ⚠ 进度文件已过期（超过 3 个监控周期未更新）——monitor 可能挂了：autozt monitor -d")
    rnd = view.get("round")
    if rnd:
        print("  上一轮：%s" % json.dumps(rnd, ensure_ascii=False, sort_keys=True))
    for proj, by_tt in sorted((view.get("projects") or {}).items()):
        for key, row in sorted(by_tt.items()):
            st = row.get("states") or {}
            print("  %-32s %-16s %4d 材料  done=%d R=%d PD=%d FAIL=%d 其它=%d"
                  % (proj[:32], key[:16], row["materials"], st.get("done", 0),
                     st.get("R", 0), st.get("PD", 0), st.get("FAIL", 0),
                     row["materials"] - sum(st.get(k, 0) for k in
                                            ("done", "R", "PD", "FAIL"))))
    for act, b in sorted((view.get("fail_summary") or {}).items()):
        top = sorted(b["codes"].items(), key=lambda kv: -kv[1])[:4]
        print("  FAIL 分诊 %-13s %3d  %s" % (act, b["count"],
                                            ", ".join("%s×%d" % kv for kv in top)))
    fails = [(e, f) for e in view["materials"] for f in e.get("fails") or []]
    shown = fails if limit is None else fails[:int(limit)]
    for e, f in shown[:40] if limit is None else shown:
        print("  FAIL %s %s  [%s/%s → %s]%s"
              % (e["id"], f.get("step"), f.get("class"), f.get("code"), f.get("action"),
                 ("  已失败 %d 次" % f["fail_count"]) if f.get("fail_count") else ""))
    if limit is None and len(fails) > 40:
        print("  … 还有 %d 条 FAIL（--json 或 -status error 看全部）" % (len(fails) - 40))
    return 0
