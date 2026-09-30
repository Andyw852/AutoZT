# -*- coding: utf-8 -*-
"""agent 读状态 / 执行动作的共享层（``autozt agent`` 与 ``autozt mcp`` 共用）。

为什么要它（2026-09 生产实测）：inspect/snapshot/propose/cycle 每调一次都要
``autozt list --json`` 把整个技能现采一遍——328 个材料在 WSL 9p 上几十分钟；执行时
每个动作再各起一条 ``act``、各自再全量采集一遍。结果是 AI 看一眼状态要半小时，
cycle 一轮 20 个动作要十几个小时，而且大范围下两次现采之间总有作业变状态，
整批 cursor 校验几乎必然判「过期」。

这里的做法：

* **读**：``source=auto``（默认）优先读 monitor 每轮写的 ``.tf_progress.json``——
  秒级、不连超算；文件太旧/没有/不含所问范围时自动回退现采。``live`` 强制现采，
  ``progress`` 强制读文件。返回里带 ``source`` 说明数据来自哪、多旧。
* **执行**：同一 (动作, 技能, 步骤) 的动作合成**一条** ``act -p A,B,C`` —— 一次采集、
  一个共享的 max_jobs 闸门；每个动作带 ``--expect-state`` 前置条件，由 CLI 按
  **现采**状态逐个核对（计划来自旧文件也不会把已算完的步骤重交、把在跑的作业
  retry 掉）。CLI 的 ``[agent-event]`` 行回填每个动作的结果：submitted / deferred
  （达 max_jobs）/ skipped_stale（状态已变）/ ok / failed。
"""

from __future__ import annotations

import json
import os
import time
from collections import OrderedDict
from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple

from autozt import agent_protocol as protocol

SOURCES = ("auto", "progress", "live")
VERBS = {"start_step": "start", "prepare_step": "retry",
         "sync_results": "fetch", "run_ready_steps": "advance"}
# 计划里的状态 → 执行前现采必须仍是这些状态之一（与 proposals 的生成条件一致）
EXPECT = {"start_step": "TODO,PREP,SCANCEL", "prepare_step": "FAIL"}
EVENT_PREFIX = "[agent-event] "


def load_cfg(config_path: Optional[str] = None) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """与主 CLI 相同的配置解析（只读 tf.yaml，不扫项目、不采集）。"""
    try:
        from autozt import load_config
        path = (config_path or os.environ.get("AUTOZT_CONFIG")
                or os.environ.get("AUTOZT_CFG") or "").strip() or None
        cfg, found = load_config(path)
    except SystemExit as exc:
        return None, str(exc.code)
    except Exception as exc:                      # noqa: BLE001
        return None, "cannot load config: %s" % exc
    cfg["_config_dir"] = os.path.dirname(os.path.abspath(found)) if found else os.getcwd()
    return cfg, None


def compact_from_progress(doc: Mapping[str, Any], scope: Optional[Mapping[str, Any]] = None
                          ) -> Dict[str, Any]:
    """进度文件 → 与 ``agent_protocol.compact_state`` 同形的紧凑状态（proposals 直接可用）。"""
    from autozt import progress as _progress
    scope = scope or {}
    view = _progress.filter_progress(doc, scope.get("project"), scope.get("tt"),
                                     scope.get("material"), scope.get("status"))
    out: Dict[str, Any] = {"types": [], "queue": doc.get("queue") or {}}
    by_tt: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
    for e in view.get("materials") or []:
        ct = by_tt.setdefault(e.get("tt") or "", {"key": e.get("tt"), "materials": []})
        fails = {str(f.get("step")): f for f in e.get("fails") or []}
        cm: Dict[str, Any] = {"name": e.get("name"), "id": e.get("id") or e.get("name"),
                              "project": e.get("project") or "", "hpc": e.get("hpc") or "",
                              "dim": "", "active": e.get("active"),
                              "action": e.get("action") or "", "steps": []}
        for label, kind in (e.get("steps") or {}).items():
            cs: Dict[str, Any] = {"label": label, "kind": kind}
            f = fails.get(str(label))
            if f:
                for src, dst in (("diag", "diag"), ("code", "diag_code"),
                                 ("class", "diag_class"), ("action", "suggested_action"),
                                 ("reason", "action_reason")):
                    if f.get(src):
                        cs[dst] = f[src]
                if f.get("fail_count"):
                    cs["fail_count"] = f["fail_count"]
            if label == e.get("active") and e.get("job"):
                cs["job"] = {k: e["job"].get(k) for k in ("id", "state", "info")
                             if e["job"].get(k) not in (None, "")}
            cm["steps"].append(cs)
        ct["materials"].append(cm)
    out["types"] = [t for t in by_tt.values() if t["materials"]]
    return out


def _fresh_limit(cfg: Mapping[str, Any]) -> int:
    from autozt import narrow_max_age
    return narrow_max_age(cfg) or 7200


def choose_source(requested: Optional[str], cfg: Optional[Mapping[str, Any]],
                  scope: Optional[Mapping[str, Any]] = None
                  ) -> Tuple[str, Optional[Dict[str, Any]], Dict[str, Any]]:
    """决定这次读状态用进度文件还是现采。返回 (kind, 进度文档或 None, source 说明)。
    没指定时取 AUTOZT_AGENT_SOURCE，再默认 auto。"""
    requested = (requested or os.environ.get("AUTOZT_AGENT_SOURCE") or "auto").strip().lower()
    if requested not in SOURCES:
        requested = "auto"
    if requested == "live" or cfg is None:
        return "live", None, {"kind": "live", "requested": requested}
    from autozt import progress as _progress
    doc = _progress.load_progress_view(cfg)   # 含采集后的提交（账本叠加成 PD）
    if doc is None:
        info = {"kind": "live", "requested": requested,
                "fallback_reason": "no progress file at %s" % _progress.progress_path(cfg)}
        return ("progress" if requested == "progress" else "live"), None, info
    live = _progress.monitor_liveness(cfg, doc)
    age = live.get("age_s")
    info = {"kind": "progress", "requested": requested, "age_s": age,
            "time": doc.get("time"), "writer": doc.get("writer"),
            "monitor_alive": live.get("alive"), "stale": live.get("stale"),
            "path": _progress.progress_path(cfg),
            "note": "状态来自进度文件（不是现采）；执行动作时 CLI 会按现采状态逐个复核"}
    if requested == "progress":
        return "progress", doc, info
    limit = _fresh_limit(cfg)
    fresh = (live.get("alive") and not live.get("stale")) or (age is not None and age <= limit)
    if not fresh:
        return "live", None, {"kind": "live", "requested": requested,
                              "fallback_reason": "progress file is %ss old (limit %ss) and "
                                                 "no live monitor" % (age, limit)}
    probe = _progress.filter_progress(doc, (scope or {}).get("project"),
                                      (scope or {}).get("tt"), (scope or {}).get("material"))
    if not probe.get("materials"):
        return "live", None, {"kind": "live", "requested": requested,
                              "fallback_reason": "progress file has no materials in this scope"}
    return "progress", doc, info


def load_compact(scope: Optional[Mapping[str, Any]], source: Optional[str],
                 config_path: Optional[str],
                 live_loader: Callable[[Mapping[str, Any]], Tuple[Optional[Dict[str, Any]],
                                                                  Optional[str], int]]
                 ) -> Tuple[Optional[Dict[str, Any]], Optional[str], int, Dict[str, Any]]:
    """按 source 读紧凑状态。config_path 与主 CLI 的 -c 同义（None = AUTOZT_CONFIG/默认搜索）；
    live_loader(scope) 是各传输层自己的现采函数。配置读不到时只能现采。"""
    cfg = None
    if (source or os.environ.get("AUTOZT_AGENT_SOURCE") or "auto") != "live":
        cfg, _err = load_cfg(config_path)
    kind, doc, info = choose_source(source, cfg, scope)
    if kind == "progress":
        if doc is None:
            return None, info.get("fallback_reason") or "no progress file", 1, info
        return compact_from_progress(doc, scope), None, 0, info
    compact, error, rc = live_loader(scope or {})
    return compact, error, rc, info


# ---------------------------------------------------------------------------
# 执行：分组 + 前置条件 + 事件回填
# ---------------------------------------------------------------------------
def group_actions(actions: List[Mapping[str, Any]]) -> "OrderedDict[tuple, List[Mapping[str, Any]]]":
    groups: "OrderedDict[tuple, List[Mapping[str, Any]]]" = OrderedDict()
    for a in actions:
        mat = str(a.get("material") or "")
        # 名字里带逗号的材料不能拼进 -p A,B,C，单独一组
        solo = mat if "," in mat else ""
        key = (a.get("action"), a.get("tt") or "", a.get("step") or "",
               a.get("project") or "", solo)
        groups.setdefault(key, []).append(a)
    return groups


def group_argv(key: tuple, acts: List[Mapping[str, Any]]) -> List[str]:
    action, tt, step, project, _solo = key
    argv = ["act"]
    if project:
        argv += ["--project", str(project)]
    if tt:
        argv += ["-tt", str(tt)]
    mats: List[str] = []
    for a in acts:
        m = str(a.get("material") or "")
        if m and m not in mats:
            mats.append(m)
    if mats:
        argv += ["-p", ",".join(mats)]
    if step:
        argv += ["-j", str(step)]
    if action in EXPECT and mats:
        argv.append("--expect-state=" + EXPECT[action])
    argv.append(VERBS[action])
    return argv


def parse_events(text: str) -> List[Dict[str, Any]]:
    out = []
    for line in (text or "").splitlines():
        if line.startswith(EVENT_PREFIX):
            try:
                ev = json.loads(line[len(EVENT_PREFIX):])
            except ValueError:
                continue
            if isinstance(ev, dict):
                out.append(ev)
    return out


def _event_matches(ev: Mapping[str, Any], action: Mapping[str, Any]) -> bool:
    want = str(action.get("material") or "")
    if not want:
        return True
    for key in ("material", "id", "name"):
        val = ev.get(key)
        if val and protocol.material_matches(val, [want], ev.get("id")):
            return True
    return False


def _outcome(action: Mapping[str, Any], events: List[Mapping[str, Any]], rc: int) -> str:
    # 一组只有一个步骤（-j 相同），按材料对上即可；-j 可能写序号而事件里是 label，不比步骤。
    kinds = [ev.get("event") for ev in events if _event_matches(ev, action)]
    if "skipped" in kinds:
        return "skipped_stale"
    if "submitted" in kinds:
        return "submitted"
    if "deferred" in kinds:
        return "deferred"
    return "ok" if rc == 0 else "failed"


def execute_grouped(actions: List[Mapping[str, Any]],
                    run: Callable[[List[str]], Tuple[int, str, str]],
                    output_limit: int = 4000) -> Dict[str, Any]:
    """按组执行；遇到失败组即停（与以前「第一个失败即停」一致）。"""
    results: List[Dict[str, Any]] = []
    groups_out: List[Dict[str, Any]] = []
    for key, acts in group_actions(actions).items():
        if key[0] not in VERBS:
            return {"results": results, "groups": groups_out, "actions": list(actions),
                    "failed": True, "error": "unsupported or destructive action: %s" % key[0]}
        argv = group_argv(key, acts)
        t0 = time.time()
        rc, out, err = run(argv)
        text = ((out or "") + (("\n" + err) if (err or "").strip() else "")).strip()
        events = parse_events(text)
        human = "\n".join(l for l in text.splitlines() if not l.startswith(EVENT_PREFIX))
        groups_out.append({"argv": argv, "returncode": rc, "seconds": round(time.time() - t0, 1),
                           "actions": len(acts), "output": human[-output_limit:]})
        for a in acts:
            outcome = _outcome(a, events, rc)
            row = {"action": a.get("action"),
                   "arguments": {k: a[k] for k in ("tt", "material", "step", "project")
                                 if a.get(k) not in (None, "")},
                   "ok": rc == 0 and outcome != "failed", "outcome": outcome,
                   "returncode": rc, "group": len(groups_out) - 1,
                   "output": human[-min(output_limit, 1000):]}
            if outcome == "failed":
                row["error"] = human[-min(output_limit, 1000):] or "AutoZT command failed"
            results.append(row)
        if rc != 0:
            return {"results": results, "groups": groups_out, "actions": list(actions),
                    "failed": True, "error": "batch stopped after first failed action"}
    counts: Dict[str, int] = {}
    for r in results:
        counts[r["outcome"]] = counts.get(r["outcome"], 0) + 1
    return {"results": results, "groups": groups_out, "actions": list(actions),
            "outcomes": counts, "failed": False}


def event_env(env: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    out = dict(os.environ if env is None else env)
    out["AUTOZT_AGENT_EVENTS"] = "1"
    return out
