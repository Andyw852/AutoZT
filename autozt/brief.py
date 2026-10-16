# -*- coding: utf-8 -*-
"""autozt brief —— 材料基础信息卡（给 LLM 回复开头用）。

    autozt [-tt T] -p MAT [-j STEP] brief [--json]

一张卡回答"这个材料现在跑在哪、怎么配的、作业到哪了"：

  1. hpc       实际用的集群 + 来自哪一层配置（材料/<技能>/hpc.yaml、project_setting/hpc.yaml、技能默认）
  2. work_dir  远端工作目录 + 来源（切集群后最容易踩的坑：setting.yaml 里钉着旧集群路径）
  3. 参数      项目层相对技能出厂默认【改过】的 step.conf 键（如 FIT_ENGINE=phono3py）
  4. 作业      当前活动步骤、最近一次提交的作业号与时间、排队/运行/失败状态
  5. 自动推进  auto_advance 生效值及来源

只读、秒级：配置走 conf/dir 同一条"只读配置"快路径（不 ssh），作业状态读 monitor
写的进度文件 .tf_progress.json 与提交账本 .tf_submit_ledger.jsonl，不连超算。
卡片里的 Markdown（brief_md）就是 LLM 回复开头应原样输出的那几行。
"""

import os
import time

_SECTIONS = ("params", "submit", "incar", "incar.final")
_MAX_PARAMS = 12
_LEDGER_LOOKBACK_S = 30 * 86400

_STATE_TEXT = {
    "R": ("▶", "运行中"), "PD": ("⏳", "排队中"), "TODO": ("·", "输入已生成，待提交"),
    "PREP": ("·", "未生成输入"), "WAIT": ("·", "等上游"), "FAIL": ("❌", "失败"),
    "SCANCEL": ("⏹", "已取消"), "OK": ("✅", "已完成"), "DONE": ("✅", "全部完成"),
    "IMAG": ("⚠", "完成但有虚频"),
}


def _rel(m, path):
    """配置文件路径写成相对材料目录的形式（卡片里短、一眼看出是哪一层）。"""
    if not path:
        return path
    base = m.get("lpath") or ""
    try:
        if base and os.path.commonpath([os.path.abspath(path), os.path.abspath(base)]) \
                == os.path.abspath(base):
            return os.path.relpath(path, base)
    except ValueError:
        pass
    return path


def _ago(ts, now=None):
    try:
        sec = int((now or time.time()) - float(ts))
    except (TypeError, ValueError):
        return ""
    if sec < 90:
        return "%d 秒前" % max(sec, 0)
    if sec < 5400:
        return "%d 分钟前" % (sec // 60)
    if sec < 172800:
        return "%.1f 小时前" % (sec / 3600.0)
    return "%d 天前" % (sec // 86400)


def _fmt_ts(ts):
    try:
        return time.strftime("%m-%d %H:%M", time.localtime(float(ts)))
    except (TypeError, ValueError, OverflowError):
        return ""


# ---------------------------------------------------------------- hpc / work_dir
def hpc_info(cfg, t, m):
    from autozt import _load_yaml_file, step_cfg
    from autozt.collect import LOCAL_HOST
    ps = m.get("ps") or {}
    src = None
    sk = m.get("_skill_dir_local")
    if sk and (_load_yaml_file(os.path.join(sk, "hpc.yaml")) or {}).get("name"):
        src = _rel(m, os.path.join(sk, "hpc.yaml"))
    elif (ps.get("hpc") or {}).get("name") and ps.get("dir"):
        src = _rel(m, os.path.join(ps["dir"], "hpc.yaml"))
    else:
        src = "技能默认（task_types.%s.hpc）" % t.get("key")
    host = m.get("host_eff")
    if host == LOCAL_HOST:
        host_txt = "本机"
    elif host:
        host_txt = str(host)
    else:
        host_txt = str((cfg or {}).get("host") or "") or "本机"
    per_step = []
    for s in m.get("steps") or []:
        sh = step_cfg(t, s.get("name"), m).get("hpc")
        if sh and str(sh) != str(m.get("hpc_name")):
            per_step.append({"step": s.get("label"), "hpc": str(sh)})
    return {"name": m.get("hpc_name") or "", "source": src, "ssh_host": host_txt,
            "step_overrides": per_step}


def work_dir_info(m):
    return {"path": m.get("work_dir_eff") or "", "source": m.get("work_dir_src") or "",
            "material_dir": m.get("rpath") or m.get("path") or ""}


# ---------------------------------------------------------------- 改过的参数
def changed_params(cfg, t, m, steps):
    """项目层 step.conf 相对技能出厂默认改过的键（init 原样拷来的默认值不算）。

    返回 [{section, key, value, default, steps:[label...], source}]，同一 (节,键,值)
    跨步骤合并（写在共用 step.conf 里的键会出现在每一步）。"""
    from autozt import _stepconf_mod, step_conf_sources
    mod = _stepconf_mod(cfg, t, m)
    if not mod:
        return []
    merged_rows = {}
    order = []
    for s in steps:
        sname = s.get("name")
        try:
            srcs = step_conf_sources(cfg, t, m, sname)
        except Exception:                       # noqa: BLE001
            continue
        skill, full = [], []
        for path, note in srcs:
            try:
                with open(path, encoding="utf-8-sig") as fh:
                    txt = fh.read()
            except OSError:
                continue
            item = (_rel(m, path), txt)
            full.append(item)
            if str(note).startswith("skill"):
                skill.append(item)
        if len(full) == len(skill):
            continue                            # 没有项目层
        base, _bp = mod.merge(skill)
        final, prov = mod.merge(full)
        for sec in _SECTIONS:
            dflt = {k.upper(): v for k, v, _l in base.get(sec, []) if v is not None}
            for k, v, _l in final.get(sec, []):
                if v is None or dflt.get(k.upper()) == v:
                    continue
                rk = (sec, k.upper(), v)
                row = merged_rows.get(rk)
                if row is None:
                    row = {"section": sec, "key": k, "value": v,
                           "default": dflt.get(k.upper()), "steps": [],
                           "source": prov.get((sec, k))}
                    merged_rows[rk] = row
                    order.append(rk)
                row["steps"].append(s.get("label") or sname)
    return [merged_rows[k] for k in order]


# ---------------------------------------------------------------- 作业
def _progress_entry(doc, t, m):
    if not doc:
        return None
    ids = {m.get("qualified_name"), m.get("name")} - {None}
    for e in doc.get("materials") or []:
        if e.get("tt") == t.get("key") and (e.get("id") in ids or e.get("name") in ids):
            return e
    return None


def _last_submit(cfg, t, m, step_label=None):
    from autozt.progress import ledger_recent
    ids = {m.get("qualified_name"), m.get("name")} - {None}
    best = None
    for r in ledger_recent(cfg, time.time() - _LEDGER_LOOKBACK_S):
        if r.get("tt") != t.get("key") or r.get("id") not in ids:
            continue
        if step_label and str(r.get("step")) != str(step_label):
            continue
        if best is None or float(r.get("ts") or 0) >= float(best.get("ts") or 0):
            best = r
    return best


def job_info(cfg, t, m, doc, step_label=None, now=None):
    e = _progress_entry(doc, t, m)
    out = {"progress": bool(e)}
    if e:
        out.update({k: e.get(k) for k in ("state", "active", "done_steps", "total_steps",
                                           "since", "eta_s") if e.get(k) is not None})
        if e.get("job"):
            out["job"] = e["job"]
        if e.get("fails"):
            out["fails"] = [{k: f.get(k) for k in ("step", "code", "action", "diag")}
                            for f in e["fails"]][:3]
        if e.get("obs"):
            out["observed"] = _ago(e["obs"], now)
        if step_label:
            out["step"] = step_label
            out["step_state"] = (e.get("steps") or {}).get(str(step_label))
    last = _last_submit(cfg, t, m, step_label)
    if last:
        out["last_submit"] = {"step": last.get("step"), "jid": last.get("jid"),
                              "ts": last.get("ts"), "at": _fmt_ts(last.get("ts")),
                              "ago": _ago(last.get("ts"), now)}
    return out


# ---------------------------------------------------------------- 自动推进
def auto_info(cfg, m):
    ps = m.get("ps") or {}
    st = ps.get("setting") or {}
    if not ps.get("dir"):
        return {"on": False, "source": "未 init（不会自动推进）"}
    if "auto_advance" in st:
        return {"on": st.get("auto_advance") is True,
                "source": _rel(m, os.path.join(ps["dir"], "setting.yaml"))}
    return {"on": bool((cfg or {}).get("auto_advance")), "source": "全局 tf.yaml"}


# ---------------------------------------------------------------- 组装 + 渲染
def build_brief(cfg, data, mname, jname=None, doc=None, liveness=None, now=None):
    from autozt import find_material, find_step
    t, m = find_material(data, mname)
    steps = list(m.get("steps") or [])
    focus = None
    if jname:
        focus = find_step(m, jname)
        steps = [focus]
    out = {
        "material": m.get("name"), "id": m.get("qualified_name") or m.get("name"),
        "tt": t.get("key"),
        "hpc": hpc_info(cfg, t, m),
        "work_dir": work_dir_info(m),
        "params": changed_params(cfg, t, m, steps),
        "job": job_info(cfg, t, m, doc, focus.get("label") if focus else None, now),
        "auto_advance": auto_info(cfg, m),
        "_nsteps": len(steps),
    }
    if liveness is not None:
        out["monitor"] = {"alive": bool(liveness.get("alive")),
                          "progress_age_s": liveness.get("age_s"),
                          "stale": bool(liveness.get("stale"))}
    out["brief_md"] = render_md(out)
    out.pop("_nsteps", None)
    return out


def _job_line(j, mid):
    last = j.get("last_submit") or {}
    if not j.get("progress"):
        if last:
            return ("%s 已提交（%s，%s）；当前状态未采集——`autozt -p %s list --refresh` 刷新"
                    % (last.get("step"), last.get("jid"), last.get("at"), mid))
        return "暂无作业记录（未采集；`autozt -p %s list --refresh` 刷新）" % mid
    if j.get("state") == "DONE" and not j.get("step"):
        return "✅ 全部 %s/%s 步完成" % (j.get("done_steps"), j.get("total_steps"))
    step = j.get("step") or j.get("active") or last.get("step")
    state = j.get("step_state") if j.get("step") else j.get("state")
    job = j.get("job") if not j.get("step") or j.get("step") == j.get("active") else None
    job = job or {}
    jid = job.get("id")
    if not jid and str(last.get("step")) == str(step):
        jid = last.get("jid")
    icon, word = _STATE_TEXT.get(str(state), ("·", str(state or "?")))
    s = str(step or "?")
    if jid:
        at = last.get("at") if str(last.get("jid")) == str(jid) else ""
        s += " 已提交（%s%s）" % (jid, "，" + at if at else "")
    s += " %s %s" % (icon, word)
    if state == "PD":
        reason = (job.get("info") or "").strip("() ")
        s += "（%s）" % reason if reason else "，等排到"
    elif state == "R":
        extra = [x for x in ((job.get("info") and "@%s" % job["info"]), job.get("time")) if x]
        if extra:
            s += "（%s）" % " ".join(extra)
    elif state == "FAIL":
        f = next((x for x in j.get("fails") or [] if str(x.get("step")) == str(step)),
                 (j.get("fails") or [{}])[0])
        if f.get("code"):
            s += "：%s，建议 %s" % (f.get("code"), f.get("action") or "?")
    if j.get("done_steps") is not None and not j.get("step"):
        s += "；进度 %s/%s" % (j.get("done_steps"), j.get("total_steps"))
    if j.get("observed"):
        s += "（数据采于 %s）" % j["observed"]
    return s


def render_md(b):
    h, w, a = b["hpc"], b["work_dir"], b["auto_advance"]
    lines = ["**【AutoZT 基础信息】%s（%s）**" % (b["id"], b["tt"])]
    hl = "1. hpc → %s ✅（来源：%s；ssh：%s）" % (h["name"] or "?", h["source"], h["ssh_host"])
    if h.get("step_overrides"):
        hl += "；分步：" + "、".join("%s@%s" % (x["step"], x["hpc"])
                                   for x in h["step_overrides"][:6])
    lines.append(hl)
    lines.append("2. work_dir → `%s` %s（来源：%s）"
                 % (w["path"] or "(未解析)", "✅" if w["path"] else "⚠", w["source"] or "?"))
    ps = b.get("params") or []
    if ps:
        items = []
        for p in ps[:_MAX_PARAMS]:
            key = p["key"] if p["section"] == "params" else "%s.%s" % (p["section"], p["key"])
            n = len(p["steps"])
            where = ("全部步骤" if n > 1 and n >= b.get("_nsteps", 0)
                     else "/".join(p["steps"][:4]) + ("…" if n > 4 else ""))
            items.append("%s=%s（%s）" % (key, p["value"], where))
        more = "；另 %d 项" % (len(ps) - _MAX_PARAMS) if len(ps) > _MAX_PARAMS else ""
        lines.append("3. 参数 → %s ✅%s" % ("，".join(items), more))
    else:
        lines.append("3. 参数 → %s全部为技能默认值" % ("本步" if b["job"].get("step") else ""))
    lines.append("4. 作业 → " + _job_line(b["job"], b["id"]))
    al = "5. 自动推进 → %s（%s）" % ("开" if a["on"] else "关", a["source"])
    mon = b.get("monitor")
    if mon is not None:
        al += "；monitor %s" % ("在跑" if mon["alive"] else "未运行")
    lines.append(al)
    return "\n".join(lines)


def cmd_brief(cfg, data, projs, jname=None, json_out=False):
    import json
    from autozt import progress as _progress
    doc = _progress.load_progress_view(cfg)
    try:
        live = _progress.monitor_liveness(cfg, doc)
    except Exception:                           # noqa: BLE001
        live = None
    outs = [build_brief(cfg, data, pj, jname, doc=doc, liveness=live) for pj in projs]
    if json_out:
        print(json.dumps(outs[0] if len(outs) == 1 else outs, ensure_ascii=False, indent=2))
    else:
        print("\n\n".join(o["brief_md"] for o in outs))
    return 0
