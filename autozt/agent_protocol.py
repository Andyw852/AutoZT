# -*- coding: utf-8 -*-
"""Shared model-facing protocol helpers.

This module contains no transport code and never submits a job.  The MCP server
and ``autozt agent`` CLI both use it so that scope filtering, compact state,
skill contracts, snapshots, and deterministic proposals cannot drift apart.
"""

from __future__ import annotations

import hashlib
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set


# Shared by the JSON CLI, MCP capabilities document, and the service layer.
# The MCP transport has its own result envelope version; this is the model-facing
# contract version that callers should negotiate and cache.
PROTOCOL_VERSION = "agent/6"
ALLOWED_ACTIONS = {
    "start_step", "prepare_step", "sync_results", "run_ready_steps",
}
MAX_ACTIONS = 20


def _value(args: Any, key: str, default: Any = None) -> Any:
    if isinstance(args, Mapping):
        return args.get(key, default)
    return getattr(args, key, default)


def compact_json(value: Any) -> str:
    import json

    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"))


def state_cursor(compact: Mapping[str, Any]) -> str:
    """Return the stable identity of one compact, scoped state snapshot.

    The cursor deliberately hashes only the compact state sent to the model.
    It is therefore cheap to compare, independent of transport, and safe to
    persist in a plan without storing remote files or command output.
    """
    raw = compact_json(compact).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:20]


def split_scope(value: Any) -> List[str]:
    return [x.strip() for x in str(value or "").split(",") if x.strip()]


def status_kinds(spec: Any) -> Set[str]:
    aliases = {
        "done": "OK", "ok": "OK", "running": "R", "run": "R",
        "pd": "PD", "queued": "PD", "error": "FAIL", "fail": "FAIL",
        "waiting": "WAIT", "wait": "WAIT", "scancel": "SCANCEL",
        "scancelled": "SCANCEL", "cancelled": "SCANCEL", "todo": "TODO",
        "prep": "PREP", "other": "OTHER",
    }
    kinds: Set[str] = set()
    for raw in split_scope(spec):
        key = raw.lower()
        if key == "ready":
            kinds.update(("TODO", "PREP"))
        else:
            kinds.add(aliases.get(key, raw.upper()))
    return kinds


def material_matches(name: Any, wanted: Iterable[str]) -> bool:
    wanted = list(wanted)
    if not wanted:
        return True
    text = str(name or "")
    base = text.rstrip("/").rsplit("/", 1)[-1]
    return any(text == item or base == item or text.endswith("/" + item)
               for item in wanted)


def scope_data(data: Any, args: Any = None) -> Any:
    """Apply model-requested type/material/status scope to a collected payload."""
    if not isinstance(data, dict) or not isinstance(data.get("types"), list):
        return data
    wanted_tt = set(split_scope(_value(args, "tt")))
    wanted_mat = split_scope(_value(args, "material"))
    wanted_status = status_kinds(_value(args, "status"))
    out = dict(data)
    types = []
    for task_type in data.get("types") or []:
        if wanted_tt and task_type.get("key") not in wanted_tt:
            continue
        ct = dict(task_type)
        mats = []
        for material in task_type.get("materials") or []:
            if not material_matches(material.get("name"), wanted_mat):
                continue
            steps = material.get("steps") or []
            if wanted_status and not any(s.get("kind") in wanted_status for s in steps):
                continue
            mats.append(material)
        ct["materials"] = mats
        types.append(ct)
    out["types"] = types
    return out


def active_label(material: Mapping[str, Any]) -> Any:
    active = material.get("active")
    if isinstance(active, Mapping):
        return active.get("label") or active.get("name")
    return active


def compact_state(data: Any, args: Any = None) -> Dict[str, Any]:
    """Keep only fields needed for planning, display, and state comparison."""
    data = scope_data(data, args)
    out: Dict[str, Any] = {"types": [], "queue": data.get("queue") or {}}
    for task_type in data.get("types") or []:
        compact_type = {"key": task_type.get("key"), "materials": []}
        for material in task_type.get("materials") or []:
            cm: Dict[str, Any] = {
                "name": material.get("name"),
                "hpc": material.get("hpc_name") or material.get("hpc") or "",
                "dim": material.get("dim") or "",
                "active": active_label(material),
                "action": material.get("action") or "",
                "steps": [],
            }
            for step in material.get("steps") or []:
                cs: Dict[str, Any] = {
                    "label": step.get("label"), "kind": step.get("kind")}
                for key in ("diag", "diag_code", "suggested_action", "action_reason"):
                    if step.get(key):
                        cs[key] = step[key]
                job = step.get("job") or {}
                if job:
                    cs["job"] = {k: job.get(k) for k in ("id", "state", "info")
                                 if job.get(k) not in (None, "")}
                cm["steps"].append(cs)
            compact_type["materials"].append(cm)
        if compact_type["materials"]:
            out["types"].append(compact_type)
    return out


def compact_contract(contract: Any) -> Any:
    """Drop generator paths and prose while preserving the model contract."""
    if not isinstance(contract, dict):
        return contract
    out = {k: contract.get(k) for k in
           ("skill", "version", "desc", "schema", "corrections", "issues")
           if contract.get(k) not in (None, [], "")}
    steps: List[Dict[str, Any]] = []
    for step in contract.get("steps") or []:
        if isinstance(step, dict):
            item = {k: step.get(k) for k in
                    ("name", "label", "seq", "check", "run", "needs", "produces")
                    if step.get(k) not in (None, [], "")}
            if item:
                steps.append(item)
    if steps:
        out["steps"] = steps
    io = contract.get("io_schema")
    if isinstance(io, dict):
        out["io_schema"] = {k: io.get(k) for k in
                            ("inputs", "outputs", "params", "notes")
                            if io.get(k) not in (None, [], "")}
    elif io not in (None, [], ""):
        out["io_schema"] = io
    flow = contract.get("flow")
    if isinstance(flow, dict):
        out["flow"] = {k: flow.get(k) for k in
                        ("summary", "stages", "requires", "backend", "title", "tagline")
                        if flow.get(k) not in (None, [], "")}
    elif flow not in (None, [], ""):
        out["flow"] = flow
    if contract.get("stats"):
        out["stats"] = contract["stats"]
    missing = [key for key in ("io_schema", "flow", "corrections")
               if not contract.get(key)]
    io_decl = contract.get("io_schema") if isinstance(contract.get("io_schema"), dict) else {}
    outputs = io_decl.get("outputs") or []
    normalized = any(isinstance(x, dict) and any(k in x for k in ("type", "unit", "validator", "provenance"))
                     for x in outputs)
    out["contract_status"] = {
        "complete": not missing, "missing": missing,
        "declarative_only": not normalized,
        "input_binding": any(isinstance(x, dict) and x.get("from") for x in (io_decl.get("inputs") or [])),
        "output_normalization": normalized,
        "readiness": "executable_contract" if not missing and normalized else ("declared_only" if not missing else "incomplete"),
    }
    if missing and steps:
        out["inferred"] = {"step_graph": [
            {k: step.get(k) for k in ("label", "seq", "needs", "check")
             if step.get(k) not in (None, [], "")}
            for step in steps
        ]}
    return out


def compact_skill_list(payload: Any) -> Dict[str, Any]:
    skills = []
    for item in (payload.get("skills") if isinstance(payload, dict) else []) or []:
        if not isinstance(item, dict):
            continue
        row = {k: item.get(k) for k in ("skill", "version", "desc", "schema", "issues")
               if item.get(k) not in (None, [], "")}
        row["step_count"] = len(item.get("steps") or [])
        if item.get("stats"):
            row["coverage"] = item["stats"]
        missing = [key for key in ("io_schema", "flow", "corrections")
                   if not item.get(key)]
        io_decl = item.get("io_schema") if isinstance(item.get("io_schema"), dict) else {}
        normalized = any(isinstance(x, dict) and any(k in x for k in ("type", "unit", "validator", "provenance"))
                         for x in (io_decl.get("outputs") or []))
        row["contract_status"] = {
            "complete": not missing, "missing": missing,
            "declarative_only": not normalized,
            "input_binding": any(isinstance(x, dict) and x.get("from") for x in (io_decl.get("inputs") or [])),
            "output_normalization": normalized,
            "readiness": "executable_contract" if not missing and normalized else ("declared_only" if not missing else "incomplete"),
        }
        skills.append(row)
    return {"schema": payload.get("schema") if isinstance(payload, dict) else None,
            "count": len(skills), "skills": skills}


def state_counts(compact: Mapping[str, Any]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for task_type in compact.get("types") or []:
        for material in task_type.get("materials") or []:
            for step in material.get("steps") or []:
                kind = step.get("kind") or "unknown"
                counts[kind] = counts.get(kind, 0) + 1
    return counts


def snapshot_flat(compact: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for task_type in compact.get("types") or []:
        tt = task_type.get("key") or ""
        for material in task_type.get("materials") or []:
            name = material.get("name") or ""
            key = "%s\x1f%s\x1f__material__" % (tt, name)
            out[key] = {"type": tt, "material": name, "step": "__material__",
                        "active": material.get("active") or "",
                        "hpc": material.get("hpc") or "", "dim": material.get("dim") or "",
                        "action": material.get("action") or ""}
            for step in material.get("steps") or []:
                key = "%s\x1f%s\x1f%s" % (tt, name, step.get("label") or "")
                job = step.get("job") or {}
                out[key] = {"type": tt, "material": name, "step": step.get("label"),
                            "kind": step.get("kind"),
                            "diag_code": step.get("diag_code") or "",
                            "queue_reason": job.get("info") or ""}
    return out


def snapshot_changes(before: Mapping[str, Any], after: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Describe changed material/step rows between two compact snapshots."""
    old_flat, new_flat = snapshot_flat(before), snapshot_flat(after)
    changes: List[Dict[str, Any]] = []
    for key in sorted(set(old_flat) | set(new_flat)):
        if old_flat.get(key) != new_flat.get(key):
            current = new_flat.get(key) or old_flat.get(key) or {}
            changes.append({"type": current.get("type"),
                            "material": current.get("material"),
                            "step": current.get("step"),
                            "old": old_flat.get(key),
                            "new": new_flat.get(key)})
    return changes


def filter_snapshot(compact: Mapping[str, Any], view: str = "attention") -> Dict[str, Any]:
    if view == "all":
        return dict(compact)
    out: Dict[str, Any] = {"types": [], "queue": compact.get("queue") or {}}
    for task_type in compact.get("types") or []:
        ct = {"key": task_type.get("key"), "materials": []}
        for material in task_type.get("materials") or []:
            active = material.get("active")
            if view == "active":
                steps = [s for s in material.get("steps") or []
                         if s.get("label") == active or s.get("kind") == "FAIL"]
            else:
                steps = [s for s in material.get("steps") or []
                         if s.get("kind") not in ("OK", "WAIT") or s.get("label") == active]
            if steps:
                row = dict(material)
                row["steps"] = steps
                ct["materials"].append(row)
        if ct["materials"]:
            out["types"].append(ct)
    return out


def proposals(compact: Mapping[str, Any], include_monitoring: bool = False,
              max_items: int = MAX_ACTIONS) -> List[Dict[str, Any]]:
    """Turn stable status fields into bounded, explainable candidate actions.

    Large projects can have thousands of ready or failed steps. Keep a model-facing
    response bounded even when the caller does not scope by skill or material.
    """
    try:
        limit = max(1, min(int(max_items), MAX_ACTIONS))
    except (TypeError, ValueError):
        limit = MAX_ACTIONS
    proposals: List[Dict[str, Any]] = []
    tool_for = {"retry": "prepare_step", "start": "start_step", "rerun": "rebuild_step"}
    risk_for = {"prepare_step": "mutate", "start_step": "mutate", "rebuild_step": "destructive"}
    for task_type in compact.get("types") or []:
        tt = task_type.get("key")
        for material in task_type.get("materials") or []:
            active = material.get("active")
            for step in material.get("steps") or []:
                kind = step.get("kind")
                suggested = step.get("suggested_action")
                if kind == "FAIL":
                    tool = tool_for.get(suggested)
                    risk = risk_for.get(tool, "read")
                    proposals.append({
                        "tool": tool,
                        "risk": risk,
                        "tt": tt,
                        "material": material.get("name"),
                        "step": step.get("label"),
                        "reason_code": step.get("diag_code") or "unknown",
                        "reason": step.get("action_reason") or step.get("diag") or "",
                        "requires_approval": risk == "destructive",
                        "requires_human_review": tool is None,
                    })
                elif step.get("label") == active and kind in ("TODO", "PREP", "SCANCEL"):
                    proposals.append({"tool": "start_step", "risk": "mutate", "tt": tt,
                                      "material": material.get("name"), "step": step.get("label"),
                                      "reason_code": "ready", "reason": "active step is ready",
                                      "requires_approval": False})
                elif include_monitoring and step.get("label") == active and kind in ("R", "PD"):
                    proposals.append({"tool": None, "risk": "read", "tt": tt,
                                      "material": material.get("name"), "step": step.get("label"),
                                      "reason_code": "monitor", "reason": "job is running or queued"})
                if len(proposals) >= limit:
                    return proposals
    return proposals


def executable_actions(items: Iterable[Mapping[str, Any]], max_actions: int = MAX_ACTIONS,
                       include_retry: bool = True) -> List[Dict[str, Any]]:
    """Select unique non-destructive actions from proposals, bounded for one cycle."""
    out: List[Dict[str, Any]] = []
    seen = set()
    for item in items:
        tool = item.get("tool") or item.get("action")
        action = tool
        if action not in ALLOWED_ACTIONS:
            continue
        if action == "prepare_step" and not include_retry:
            continue
        if item.get("requires_approval") or item.get("requires_human_review"):
            continue
        row = {"action": action}
        for key in ("tt", "material", "step"):
            if item.get(key) not in (None, ""):
                row[key] = str(item[key])
        if action != "run_ready_steps" and not row.get("material"):
            continue
        key = tuple(sorted(row.items()))
        if key not in seen:
            out.append(row)
            seen.add(key)
        if len(out) >= max_actions:
            break
    return out
