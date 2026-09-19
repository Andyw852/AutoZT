# -*- coding: utf-8 -*-
"""autozt agent —— 面向脚本和 LLM 的稳定 JSON CLI。

这层是 AutoZT 的命令行 façade，不实现新的工作流语义：

* 技能契约来自现有 ``autozt schema --json``；
* 状态来自现有只读接口 ``autozt list --json``；
* 变更动作通过现有 ``autozt act`` 网关执行；
* 快照只在本地配置目录保存一个小型状态文件，避免每次把全量 JSON 交给模型。

因此即使没有 MCP 客户端，人或普通脚本也能使用同一套“发现、观察、规划、执行”接口。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from autozt import agent_protocol as _protocol
from autozt import agent_service as _service
from autozt import science as _science


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "autozt")
# Keep the outer envelope on the same version as the shared model-facing
# contract. MCP has a separate transport result envelope version.
SCHEMA_VERSION = _protocol.PROTOCOL_VERSION
SNAPSHOT_VERSION = "1"
SNAPSHOT_LIMIT = 20
ALLOWED_ACTIONS = {
    "start_step": "start",
    "prepare_step": "retry",
    "sync_results": "fetch",
    "run_ready_steps": "advance",
}
REQUEST_OPS = {
    "capabilities", "schema", "skills", "contract", "snapshot", "inspect", "plan",
    "cycle", "run", "evidence", "propose", "apply", "research_plan", "preflight", "results",
}
REQUEST_FIELDS = {
    "op", "scope", "tt", "material", "status", "step", "view",
    "include_monitoring", "include_retry", "execute", "dry_run", "max_actions",
    "cursor", "full", "skill", "plan", "actions", "goal", "result_dir", "property",
    "temperature", "carrier", "direction", "dimension", "thickness",
}


def _json_compact(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"))


def _envelope(command: str, data: Any = None, *, ok: bool = True,
              error: Optional[str] = None, returncode: int = 0) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "ok": bool(ok),
        "command": command,
        "returncode": int(returncode),
    }
    if data is not None:
        out["data"] = data
        if isinstance(data, dict) and data.get("conversation_state"):
            # Keep the full scientific payload under data while making the next
            # dialogue turn discoverable without schema-specific parsing.
            out["conversation"] = {
                key: data.get(key) for key in (
                    "schema_version", "conversation_state", "message", "next_action",
                    "requires_user_confirmation", "confirmation_payload", "plan_id")
                if key in data
            }
    if error:
        out["error"] = str(error)
    return out


def _print(value: Dict[str, Any], pretty: bool = False) -> None:
    if pretty:
        print(json.dumps(value, ensure_ascii=False, indent=2))
    else:
        print(_json_compact(value))


def _config_path(explicit: Optional[str]) -> Optional[str]:
    """Resolve the same config path used by the main CLI when possible."""
    path = explicit or os.environ.get("AUTOZT_CONFIG") or os.environ.get("AUTOZT_CFG")
    if path:
        return os.path.abspath(os.path.expanduser(path))
    try:
        from autozt import load_config

        _cfg, found = load_config(None)
        return os.path.abspath(found) if found else None
    except Exception:
        return None


def _config_dir(explicit: Optional[str]) -> str:
    path = _config_path(explicit)
    return os.path.dirname(path) if path else os.getcwd()


def _cfg_env(explicit: Optional[str]) -> Dict[str, str]:
    env = dict(os.environ)
    if explicit:
        env["AUTOZT_CONFIG"] = os.path.abspath(os.path.expanduser(explicit))
    # Machine-readable output must be UTF-8 even when the parent is a Windows
    # console whose legacy code page cannot encode skill descriptions.
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    # Direct agent-cli calls are audited by the same gate as MCP calls. Mutating
    # actions still go through `autozt act`; this actor label distinguishes the source.
    env.setdefault("AUTOZT_ACTOR", "agent-cli")
    return env


def _run(argv: Iterable[str], config: Optional[str], timeout: int = 1800
         ) -> Tuple[int, str, str]:
    cmd = [sys.executable, PROG]
    if config:
        cmd += ["-c", os.path.abspath(os.path.expanduser(config))]
    cmd += [str(x) for x in argv]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, env=_cfg_env(config))
    except subprocess.TimeoutExpired:
        return 1, "", "timeout after %ss" % timeout
    except OSError as exc:
        return 127, "", "cannot execute autozt: %s" % exc
    return p.returncode, p.stdout or "", p.stderr or ""


def _parse_json_output(out: str, err: str = "") -> Tuple[Optional[Any], Optional[str]]:
    try:
        return json.loads((out or "").strip()), None
    except (TypeError, ValueError) as exc:
        detail = ((out or "") + (("\n" + err) if err else "")).strip()
        return None, "AutoZT did not return valid JSON: %s%s" % (
            exc, ("; output: " + detail[:1000]) if detail else "")


def _cli_error(rc: int, out: str, err: str) -> str:
    text = ((out or "") + (("\n" + err) if err.strip() else "")).strip()
    return text or "AutoZT command failed with exit code %d" % rc


def _context_argv(args: argparse.Namespace) -> List[str]:
    out: List[str] = []
    if getattr(args, "tt", None):
        out += ["-tt", str(args.tt)]
    if getattr(args, "material", None):
        out += ["-p", str(args.material)]
    if getattr(args, "status", None):
        out += ["-status", str(args.status)]
    if getattr(args, "step", None):
        out += ["-j", str(args.step)]
    return out


def _call_json(argv: List[str], config: Optional[str]) -> Tuple[Optional[Any], Optional[str], int]:
    rc, out, err = _run(argv, config)
    if rc:
        return None, _cli_error(rc, out, err), rc
    data, parse_error = _parse_json_output(out, err)
    return data, parse_error, 1 if parse_error else 0


def _compact_contract(contract: Any) -> Any:
    return _protocol.compact_contract(contract)


def _compact_skill_list(payload: Any) -> Dict[str, Any]:
    return _protocol.compact_skill_list(payload)


def _active_label(material: Dict[str, Any]) -> Any:
    active = material.get("active")
    return (active.get("label") or active.get("name")) if isinstance(active, dict) else active


def _compact_state(data: Any, args: Any = None) -> Dict[str, Any]:
    scope = vars(args) if hasattr(args, "__dict__") else (args or {})
    return _protocol.compact_state(data, scope)


def _state_counts(compact: Dict[str, Any]) -> Dict[str, int]:
    return _protocol.state_counts(compact)


def _snapshot_flat(compact: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    return _protocol.snapshot_flat(compact)


def _filter_snapshot(compact: Dict[str, Any], view: str) -> Dict[str, Any]:
    if view == "all":
        return compact
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


def _snapshot_path(args: argparse.Namespace) -> str:
    scope = "|".join([str(getattr(args, key, "") or "")
                       for key in ("tt", "material", "status", "view")])
    digest = hashlib.sha256(scope.encode("utf-8")).hexdigest()[:16]
    return os.path.join(_config_dir(getattr(args, "config", None)),
                        ".tf_agent_snapshot_%s.json" % digest)


def _load_snapshot(path: str) -> Optional[Dict[str, Any]]:
    try:
        with open(path, encoding="utf-8") as handle:
            value = json.load(handle)
        return value if isinstance(value, dict) and value.get("version") == SNAPSHOT_VERSION else None
    except (OSError, ValueError):
        return None


def _save_snapshot(path: str, value: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, sort_keys=True)
    os.replace(tmp, path)


def _status(args: argparse.Namespace) -> Tuple[Optional[Dict[str, Any]], Optional[str], int]:
    data, error, rc = _call_json(_context_argv(args) + ["list", "--json"],
                                 getattr(args, "config", None))
    if error:
        return None, error, rc
    return _compact_state(data, args), None, 0


def _proposals(compact: Dict[str, Any], include_monitoring: bool = False,
               max_items: int = _protocol.MAX_ACTIONS) -> List[Dict[str, Any]]:
    """Use the shared proposal policy; MCP and CLI must never drift."""
    return _protocol.proposals(compact, include_monitoring, max_items=max_items)


def _read_plan(path: str) -> Tuple[Optional[Any], Optional[str]]:
    try:
        text = sys.stdin.read() if path == "-" else open(path, encoding="utf-8").read()
    except OSError as exc:
        return None, "cannot read plan %s: %s" % (path, exc)
    try:
        return json.loads(text), None
    except ValueError as exc:
        return None, "plan is not valid JSON: %s" % exc


def _normalise_actions(plan: Any) -> Tuple[Optional[List[Dict[str, Any]]], Optional[str]]:
    if isinstance(plan, dict) and isinstance(plan.get("data"), dict):
        plan = plan["data"]
    if isinstance(plan, dict):
        actions = plan.get("actions")
        if actions is None:
            proposals = plan.get("proposals")
            if isinstance(proposals, list):
                actions = proposals
    elif isinstance(plan, list):
        actions = plan
    else:
        actions = None
    if not isinstance(actions, list) or not actions:
        return None, "plan must contain a non-empty actions or proposals array"
    if len(actions) > SNAPSHOT_LIMIT:
        return None, "plan contains too many actions (maximum %d)" % SNAPSHOT_LIMIT
    out: List[Dict[str, Any]] = []
    for index, raw in enumerate(actions):
        if not isinstance(raw, dict):
            return None, "action[%d] must be an object" % index
        action = raw.get("action") or raw.get("tool")
        if action not in ALLOWED_ACTIONS:
            return None, "action[%d] is destructive or unsupported: %s" % (index, action)
        if raw.get("requires_approval") or raw.get("requires_human_review"):
            return None, "action[%d] is marked for human review and cannot be applied automatically" % index
        item = {"action": action}
        for key in ("tt", "material", "step"):
            if raw.get(key) not in (None, ""):
                item[key] = str(raw[key])
        if action != "run_ready_steps" and not item.get("material"):
            return None, "action[%d] requires material" % index
        out.append(item)
    return out, None


def _cmd_skills(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    payload, error, rc = _call_json(["schema", "--json"], getattr(args, "config", None))
    if error:
        return _envelope("skills", ok=False, error=error, returncode=rc), rc
    return _envelope("skills", _compact_skill_list(payload)), 0


def _cmd_capabilities(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    return _envelope("capabilities", _service.capabilities()), 0


def _cmd_schema(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    """Return the complete protocol contract without collecting cluster state."""
    return _envelope("schema", _service.schema()), 0


def _cmd_evidence(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    if not getattr(args, "material", None):
        return _envelope("evidence", ok=False,
                         error="evidence requires -p/--material", returncode=2), 2
    payload, error, rc = _call_json(_context_argv(args) + ["diagnose"],
                                    getattr(args, "config", None))
    if error:
        return _envelope("evidence", ok=False, error=error, returncode=rc), rc
    return _envelope("evidence", payload), 0


def _cmd_research_plan(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    payload, error, rc = _call_json(["schema", "--json"], getattr(args, "config", None))
    if error:
        return _envelope("research_plan", ok=False, error=error, returncode=rc), rc
    skills = (payload or {}).get("skills", []) if isinstance(payload, dict) else []
    data = _science.research_plan(getattr(args, "goal", ""), skills,
                                  dimension=getattr(args, "dimension", None),
                                  temperature=getattr(args, "temperature", None),
                                  carrier=getattr(args, "carrier", None),
                                  material=getattr(args, "material", None))
    return _envelope("research_plan", data), 0


def _cmd_preflight(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    root = getattr(args, "result_dir", None)
    if not root:
        return _envelope("preflight", ok=False, error="preflight requires --result-dir", returncode=2), 2
    data = _science.preflight(root, dimension=getattr(args, "dimension", None),
                              thickness=getattr(args, "thickness", None),
                              temperature=getattr(args, "temperature", None),
                              carrier=getattr(args, "carrier", None))
    return _envelope("preflight", data), 0


def _cmd_results(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    root = getattr(args, "result_dir", None)
    prop = getattr(args, "property", None)
    if not root or not prop:
        return _envelope("results", ok=False, error="results requires --result-dir and --property", returncode=2), 2
    data = _science.query_results(root, prop, temperature=getattr(args, "temperature", None),
                                  carrier=getattr(args, "carrier", None), direction=getattr(args, "direction", None))
    return _envelope("results", data), 0


def _cmd_contract(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    skill = getattr(args, "skill", None) or getattr(args, "tt", None)
    argv = ["schema"]
    if skill:
        argv.append(str(skill))
    argv.append("--json")
    payload, error, rc = _call_json(argv, getattr(args, "config", None))
    if error:
        return _envelope("contract", ok=False, error=error, returncode=rc), rc
    if skill:
        skills = payload.get("skills") if isinstance(payload, dict) else []
        if not skills:
            return _envelope("contract", ok=False, error="skill not found: %s" % skill,
                             returncode=1), 1
        item = skills[0] if args.full else _compact_contract(skills[0])
        return _envelope("contract", item), 0
    items = payload.get("skills") or []
    return _envelope("contract", {"count": len(items),
                                   "skills": items if args.full else
                                   [_compact_contract(x) for x in items]}), 0


def _cmd_snapshot(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    compact, error, rc = _status(args)
    if error:
        return _envelope("snapshot", ok=False, error=error, returncode=rc), rc
    cursor = _protocol.state_cursor(compact)
    path = _snapshot_path(args)
    previous = _load_snapshot(path)
    record = {"version": SNAPSHOT_VERSION, "cursor": cursor, "snapshot": compact,
              "updated": int(time.time())}
    _save_snapshot(path, record)
    if previous is None:
        data = {"cursor": cursor, "first_snapshot": True,
                "counts": _state_counts(compact),
                "snapshot": _filter_snapshot(compact, args.view)}
    else:
        before = previous.get("snapshot") or {}
        changes = _protocol.snapshot_changes(before, compact)
        data = {"cursor": cursor, "unchanged": not changes,
                "counts": _state_counts(compact), "changes": changes}
    if args.cursor and args.cursor != (previous or {}).get("cursor"):
        data["cursor_reset"] = True
    data["snapshot_file"] = path
    return _envelope("snapshot", data), 0


def _cmd_propose(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    compact, error, rc = _status(args)
    if error:
        return _envelope("propose", ok=False, error=error, returncode=rc), rc
    data = {"cursor": _protocol.state_cursor(compact),
            "scope": _agent_scope(args),
            "counts": _state_counts(compact),
            "proposals": _proposals(compact, args.include_monitoring,
                                     max_items=args.max_actions),
            "proposal_limit": args.max_actions}
    return _envelope("propose", data), 0


def _plan_field(plan: Any, field: str) -> Any:
    """Read metadata from an envelope, service result, or plain plan object."""
    if isinstance(plan, dict) and isinstance(plan.get("data"), dict):
        value = plan["data"].get(field)
        if value is not None:
            return value
        plan = plan["data"]
    return plan.get(field) if isinstance(plan, dict) else None


def _check_plan_cursor(args: argparse.Namespace, plan: Any,
                       actions: List[Dict[str, Any]]) -> Tuple[Optional[Dict[str, Any]],
                                                               Optional[str], int]:
    expected = _plan_field(plan, "cursor")
    if not expected:
        return None, None, 0
    scope = _plan_field(plan, "scope") or {}
    if not isinstance(scope, dict):
        return None, "plan scope must be an object", 2
    current_args = argparse.Namespace(
        config=getattr(args, "config", None),
        tt=scope.get("tt"), material=scope.get("material"),
        status=scope.get("status"))
    compact, error, rc = _status(current_args)
    if error:
        return None, error, rc
    actual = _protocol.state_cursor(compact or {})
    if str(expected) != actual:
        return {"stale_plan": True, "expected_cursor": str(expected),
                "current_cursor": actual, "actions": actions}, \
            "plan is stale; collect a new snapshot before applying actions", 1
    return None, None, 0


def _apply_payload(args: argparse.Namespace, plan: Any
                   ) -> Tuple[Dict[str, Any], int]:
    actions, error = _normalise_actions(plan)
    if error:
        return _envelope("apply", ok=False, error=error, returncode=2), 2
    assert actions is not None
    stale_data, stale_error, stale_rc = _check_plan_cursor(args, plan, actions)
    if stale_error:
        return _envelope("apply", stale_data, ok=False, error=stale_error,
                         returncode=stale_rc), stale_rc
    data = _execute_actions(actions, getattr(args, "config", None), bool(args.dry_run))
    data["plan_actions"] = actions
    data["dry_run"] = bool(args.dry_run)
    if data.get("failed"):
        return _envelope("apply", data, ok=False,
                         error="batch stopped after first failed action",
                         returncode=1), 1
    return _envelope("apply", data), 0


def _cmd_apply(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    plan, error = _read_plan(args.plan)
    if error:
        return _envelope("apply", ok=False, error=error, returncode=2), 2
    return _apply_payload(args, plan)


def _agent_scope(args: argparse.Namespace) -> Dict[str, Any]:
    return {key: getattr(args, key, None) for key in ("tt", "material", "status")
            if getattr(args, key, None) not in (None, "")}


def _service_loader(args: argparse.Namespace):
    def load(scope: Mapping[str, Any]):
        scoped = argparse.Namespace(
            config=getattr(args, "config", None),
            tt=scope.get("tt"), material=scope.get("material"),
            status=scope.get("status"))
        return _status(scoped)
    return load


def _cmd_inspect(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    data, error, rc = _service.inspect(
        _service_loader(args), _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS),
        include_retry=True)
    if error:
        return _envelope("inspect", ok=False, error=error, returncode=rc), rc
    return _envelope("inspect", data), 0


def _cmd_plan(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    data, error, rc = _service.plan(
        _service_loader(args), _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS))
    if error:
        return _envelope("plan", ok=False, error=error, returncode=rc), rc
    return _envelope("plan", data), 0


def _execute_actions(actions: List[Dict[str, Any]], config: Optional[str],
                     dry_run: bool = False,
                     expected_cursor: Optional[str] = None,
                     scope: Optional[Mapping[str, Any]] = None) -> Dict[str, Any]:
    """Execute only the bounded non-destructive action vocabulary through ``act``."""
    if dry_run:
        return {"results": [], "actions": actions, "dry_run": True, "failed": False}
    if expected_cursor:
        # The service has already rechecked once; this final read is the CAS
        # boundary immediately before invoking the mutation gateway.
        scoped = argparse.Namespace(config=config, tt=(scope or {}).get("tt"),
                                    material=(scope or {}).get("material"),
                                    status=(scope or {}).get("status"))
        compact, read_error, rc = _status(scoped)
        if read_error:
            return {"results": [], "actions": actions, "failed": True,
                    "error": read_error, "returncode": rc}
        current = _protocol.state_cursor(compact or {})
        if current is None or str(expected_cursor) != current:
            return {"results": [], "actions": actions, "failed": True,
                    "stale_plan": True, "expected_cursor": str(expected_cursor),
                    "current_cursor": current,
                    "error": "plan is stale; call inspect again"}
    results = []
    for action in actions:
        name = action.get("action")
        if name not in _protocol.ALLOWED_ACTIONS:
            return {"results": results, "actions": actions, "failed": True,
                    "error": "unsupported or destructive action: %s" % name}
        argv = ["act"]
        if action.get("tt"):
            argv += ["-tt", str(action["tt"])]
        if action.get("material"):
            argv += ["-p", str(action["material"])]
        if action.get("step"):
            argv += ["-j", str(action["step"])]
        verb = {"start_step": "start", "prepare_step": "retry",
                "sync_results": "fetch"}.get(name)
        argv += (["advance"] if name == "run_ready_steps" else [verb])
        rc, out, err = _run(argv, config)
        text = ((out or "") + (("\n" + err) if err.strip() else "")).strip()
        results.append({"action": name, "arguments": action, "ok": rc == 0,
                        "returncode": rc, "output": text[:4000]})
        if rc != 0:
            return {"results": results, "actions": actions, "failed": True,
                    "error": "batch stopped after first failed action"}
    return {"results": results, "actions": actions, "failed": False}


def _cmd_cycle(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    execute = bool(getattr(args, "execute", False)) and not bool(
        getattr(args, "dry_run", False))
    data, error, rc = _service.cycle(
        _service_loader(args),
        lambda actions, dry_run, cursor=None: _execute_actions(
            actions, getattr(args, "config", None), dry_run, cursor,
            _agent_scope(args)),
        _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS),
        execute=execute,
        include_retry=bool(getattr(args, "include_retry", False)))
    if error:
        return _envelope("cycle", data, ok=False, error=error, returncode=rc), rc
    return _envelope("cycle", data), 0


def _cmd_run(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    """Run the deterministic AutoZT loop; dry-run is the default."""
    execute = bool(getattr(args, "execute", False)) and not bool(
        getattr(args, "dry_run", False))
    data, error, rc = _service.run(
        _service_loader(args),
        lambda actions, dry_run, cursor=None: _execute_actions(
            actions, getattr(args, "config", None), dry_run, cursor,
            _agent_scope(args)),
        _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS),
        execute=execute,
        include_retry=bool(getattr(args, "include_retry", False)))
    if error:
        return _envelope("run", data, ok=False, error=error, returncode=rc), rc
    return _envelope("run", data), 0


def _read_request(path: str) -> Tuple[Optional[Any], Optional[str]]:
    try:
        text = sys.stdin.read() if path == "-" else open(path, encoding="utf-8").read()
    except OSError as exc:
        return None, "cannot read request %s: %s" % (path, exc)
    try:
        value = json.loads(text)
    except ValueError as exc:
        return None, "request is not valid JSON: %s" % exc
    if not isinstance(value, dict):
        return None, "request must be a JSON object"
    return value, None


def _request_args(base: argparse.Namespace, request: Mapping[str, Any]
                  ) -> argparse.Namespace:
    scope = request.get("scope") if isinstance(request.get("scope"), dict) else {}
    merged = dict(scope)
    for key in ("tt", "material", "status", "step"):
        if request.get(key) not in (None, ""):
            merged[key] = request[key]
    return argparse.Namespace(
        config=getattr(base, "config", None),
        tt=merged.get("tt"), material=merged.get("material"),
        status=merged.get("status"), step=merged.get("step"),
        view=request.get("view", "attention"),
        include_monitoring=bool(request.get("include_monitoring", False)),
        max_actions=request.get("max_actions", _protocol.MAX_ACTIONS),
        include_retry=bool(request.get("include_retry", False)),
        execute=bool(request.get("execute", False)),
        dry_run=bool(request.get("dry_run", False)),
        cursor=request.get("cursor"),
        full=bool(request.get("full", False)),
        skill=request.get("skill"),
        goal=request.get("goal", ""), result_dir=request.get("result_dir"),
        property=request.get("property"), temperature=request.get("temperature"),
        carrier=request.get("carrier"), direction=request.get("direction"),
        dimension=request.get("dimension"), thickness=request.get("thickness"),
        agent_command=str(request.get("op") or ""),
    )


def _dispatch_request(request: Any, base: argparse.Namespace
                      ) -> Tuple[Dict[str, Any], int]:
    """Dispatch one JSON object without reading or writing a transport stream."""
    if not isinstance(request, dict):
        return _envelope("request", ok=False,
                         error="request must be a JSON object", returncode=2), 2
    unknown = sorted(set(request) - REQUEST_FIELDS)
    if unknown:
        return _envelope("request", ok=False,
                         error="unknown request field(s): %s" % ", ".join(unknown),
                         returncode=2), 2
    op = str(request.get("op") or "").strip().lower()
    if not op:
        return _envelope("request", ok=False,
                         error="request requires op", returncode=2), 2
    if op not in REQUEST_OPS:
        return _envelope("request", ok=False,
                         error="unsupported request op: %s" % op,
                         returncode=2), 2
    scope = request.get("scope")
    if "scope" in request:
        if not isinstance(scope, dict):
            return _envelope("request", ok=False,
                             error="scope must be an object", returncode=2), 2
        bad_scope = sorted(set(scope) - {"tt", "material", "status", "step"})
        if bad_scope:
            return _envelope("request", ok=False,
                             error="unknown scope field(s): %s" % ", ".join(bad_scope),
                             returncode=2), 2
        for key, value in scope.items():
            if not isinstance(value, str):
                return _envelope("request", ok=False,
                                 error="scope field %s must be string" % key,
                                 returncode=2), 2
    for key in ("op", "tt", "material", "status", "step", "view", "cursor", "skill"):
        if key in request and not isinstance(request[key], str):
            return _envelope("request", ok=False,
                             error="request field %s must be string" % key,
                             returncode=2), 2
    if "view" in request and request["view"] not in ("attention", "active", "all"):
        return _envelope("request", ok=False,
                         error="request field view must be attention, active, or all",
                         returncode=2), 2
    if "actions" in request and not isinstance(request["actions"], list):
        return _envelope("request", ok=False,
                         error="request field actions must be array", returncode=2), 2
    for key in ("include_monitoring", "include_retry", "execute", "dry_run", "full"):
        if key in request and not isinstance(request[key], bool):
            return _envelope("request", ok=False,
                             error="request field %s must be boolean" % key,
                             returncode=2), 2
    if "max_actions" in request:
        value = request["max_actions"]
        if isinstance(value, bool) or not isinstance(value, int):
            return _envelope("request", ok=False,
                             error="request field max_actions must be integer",
                             returncode=2), 2
        if not 1 <= value <= _protocol.MAX_ACTIONS:
            return _envelope("request", ok=False,
                             error="request field max_actions is outside 1..%d"
                                   % _protocol.MAX_ACTIONS,
                             returncode=2), 2
    child = _request_args(base, request)
    if op == "capabilities":
        return _cmd_capabilities(child)
    if op == "schema":
        return _cmd_schema(child)
    if op == "skills":
        return _cmd_skills(child)
    if op == "contract":
        return _cmd_contract(child)
    if op == "inspect":
        return _cmd_inspect(child)
    if op == "plan":
        return _cmd_plan(child)
    if op == "cycle":
        return _cmd_cycle(child)
    if op == "run":
        return _cmd_run(child)
    if op == "evidence":
        return _cmd_evidence(child)
    if op == "research_plan":
        return _cmd_research_plan(child)
    if op == "preflight":
        return _cmd_preflight(child)
    if op == "results":
        return _cmd_results(child)
    if op == "propose":
        return _cmd_propose(child)
    if op == "snapshot":
        return _cmd_snapshot(child)
    if op == "apply":
        plan = request.get("plan")
        if plan is None:
            plan = {"actions": request.get("actions")}
        return _apply_payload(child, plan)
    # The op set is checked above; reaching this line means the dispatcher table
    # itself is incomplete and should fail loudly during development.
    return _envelope("request", ok=False,
                     error="request dispatcher has no handler for %s" % op,
                     returncode=2), 2


def _cmd_request(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    request, error = _read_request(args.request)
    if error:
        return _envelope("request", ok=False, error=error, returncode=2), 2
    return _dispatch_request(request, args)


def _cmd_serve(args: argparse.Namespace) -> int:
    """Run a persistent newline-delimited JSON transport for agents and scripts.

    Each non-empty input line receives exactly one compact JSON envelope. A bad
    line is isolated to that request; the process keeps serving later lines.
    No shell command is constructed from input, and mutations still go through
    the same ``autozt act`` path used by the one-shot CLI.
    """
    overall_rc = 0
    for line_number, line in enumerate(sys.stdin, 1):
        text = line.strip()
        if not text:
            continue
        try:
            request = json.loads(text)
        except (TypeError, ValueError) as exc:
            result, rc = _envelope(
                "request", ok=False,
                error="line %d is not valid JSON: %s" % (line_number, exc),
                returncode=2), 2
        else:
            result, rc = _dispatch_request(request, args)
        _print(result, pretty=False)
        sys.stdout.flush()
        if rc:
            overall_rc = rc
    return overall_rc


def _add_context_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-c", "--config", default=argparse.SUPPRESS,
                        help="配置文件；默认按 AutoZT 的搜索顺序")
    parser.add_argument("-tt", dest="tt", default=argparse.SUPPRESS,
                        help="技能类型")
    parser.add_argument("-p", "--material", dest="material", default=argparse.SUPPRESS,
                        help="材料名")
    parser.add_argument("-status", "--status", dest="status", default=argparse.SUPPRESS,
                        help="状态过滤，例如 error,running,pd")
    parser.add_argument("--pretty", action="store_true", default=argparse.SUPPRESS,
                        help="缩进输出 JSON")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="autozt agent",
        description="AutoZT 的稳定 JSON agent 接口：技能、快照、规划、执行。")
    _add_context_options(parser)
    sub = parser.add_subparsers(dest="agent_command")

    p = sub.add_parser("skills", help="列出紧凑技能目录")
    _add_context_options(p)
    p = sub.add_parser("capabilities", help="输出稳定的 agent 协议和动作能力")
    _add_context_options(p)
    p = sub.add_parser("schema", help="输出完整的 JSON 请求/响应/动作 Schema")
    _add_context_options(p)
    p = sub.add_parser("contract", help="输出一个技能的机器可读契约")
    _add_context_options(p)
    p.add_argument("skill", nargs="?", help="技能名；省略则输出全部契约")
    p.add_argument("--full", action="store_true",
                   help="保留生成器路径等完整契约字段；默认输出紧凑契约")

    p = sub.add_parser("snapshot", help="持久化增量状态快照")
    _add_context_options(p)
    p.add_argument("--cursor", help="上一轮 cursor；用于检测跨进程失效/重置")
    p.add_argument("--view", choices=("attention", "active", "all"), default="attention")

    p = sub.add_parser("evidence", help="读取一个材料/步骤的有界结构化诊断")
    _add_context_options(p)
    p.add_argument("--step", dest="step", help="步骤 label 或序号")

    p = sub.add_parser("research_plan", help="把研究目标转换为可审查的技能/阶段方案（只读）")
    _add_context_options(p)
    p.add_argument("--goal", required=True, help="研究目标，例如：比较二维材料的 zT")
    p.add_argument("--dimension", choices=("0D", "1D", "2D", "3D"))
    p.add_argument("--temperature", nargs="*", help="目标温度网格")
    p.add_argument("--carrier", nargs="*", help="目标载流子网格")

    p = sub.add_parser("preflight", help="检查已回收结果的跨文件科学输入一致性（只读）")
    _add_context_options(p)
    p.add_argument("--result-dir", required=True)
    p.add_argument("--dimension", choices=("0D", "1D", "2D", "3D"))
    p.add_argument("--thickness", type=float)
    p.add_argument("--temperature", nargs="*")
    p.add_argument("--carrier", nargs="*")

    p = sub.add_parser("results", help="查询带来源路径的结构化结果（只读）")
    _add_context_options(p)
    p.add_argument("--result-dir", required=True)
    p.add_argument("--property", required=True)
    p.add_argument("--temperature", type=float)
    p.add_argument("--carrier", type=float)
    p.add_argument("--direction")

    p = sub.add_parser("request", help="从 stdin 或 JSON 文件读取一个结构化请求")
    _add_context_options(p)
    p.add_argument("request", nargs="?", default="-",
                   help="JSON 文件；省略或写 - 表示 stdin")

    p = sub.add_parser("serve", help="以 JSONL 常驻服务处理多条结构化请求")
    _add_context_options(p)

    for command, help_text in (
            ("inspect", "一次读取当前关注状态并给出确定性候选动作"),
            ("plan", "一次生成可审阅的非破坏性动作计划"),
            ("cycle", "观察、规划，并可选执行一个安全工作流周期"),
            ("run", "运行确定性的观察-规划-执行闭环")):
        p = sub.add_parser(command, help=help_text)
        _add_context_options(p)
        p.add_argument("--view", choices=("attention", "active", "all"), default="attention")
        p.add_argument("--include-monitoring", action="store_true",
                       help="把正在运行/排队的步骤也列为观察项")
        p.add_argument("--max-actions", type=int, default=_protocol.MAX_ACTIONS,
                       help="单次计划最多动作数（默认 %(default)s）")
        if command in ("cycle", "run"):
            p.add_argument("--dry-run", dest="dry_run", action="store_true",
                           help="只观察和展示计划（默认行为）")
            p.add_argument("--include-retry", action="store_true",
                           help="把 retry 纳入本轮执行；默认留给模型/人工选择")
            p.add_argument("--execute", action="store_true",
                           help="执行计划中的非破坏性动作；默认只做 dry-run")

    p = sub.add_parser("propose", help="只读生成候选动作计划")
    _add_context_options(p)
    p.add_argument("--include-monitoring", action="store_true",
                   help="把正在运行/排队的步骤也列为观察项")
    p.add_argument("--max-actions", type=int, default=_protocol.MAX_ACTIONS,
                   help="最多返回候选动作数（默认 %(default)s）")

    p = sub.add_parser("apply", help="执行 JSON 计划；默认只允许非破坏性动作")
    _add_context_options(p)
    p.add_argument("plan", help="计划 JSON 文件；用 - 从 stdin 读取")
    p.add_argument("--dry-run", action="store_true", help="只验证并展示动作，不执行")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    # The CLI is a JSON transport. Force UTF-8 at the boundary so descriptions
    # containing Greek letters or Chinese remain valid on a GBK Windows console.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except AttributeError:
            pass
    parser = build_parser()
    raw = list(sys.argv[1:] if argv is None else argv)
    # bin/autozt passes the complete argv so options before `agent` survive.
    # Direct imports may pass only the subcommand; both forms are accepted.
    if "agent" in raw:
        raw.pop(raw.index("agent"))
    args = parser.parse_args(raw)
    if not args.agent_command:
        parser.print_help()
        return 2
    if args.agent_command == "skills":
        result, rc = _cmd_skills(args)
    elif args.agent_command == "capabilities":
        result, rc = _cmd_capabilities(args)
    elif args.agent_command == "schema":
        result, rc = _cmd_schema(args)
    elif args.agent_command == "contract":
        result, rc = _cmd_contract(args)
    elif args.agent_command == "snapshot":
        result, rc = _cmd_snapshot(args)
    elif args.agent_command == "evidence":
        result, rc = _cmd_evidence(args)
    elif args.agent_command == "research_plan":
        result, rc = _cmd_research_plan(args)
    elif args.agent_command == "preflight":
        result, rc = _cmd_preflight(args)
    elif args.agent_command == "results":
        result, rc = _cmd_results(args)
    elif args.agent_command == "request":
        result, rc = _cmd_request(args)
    elif args.agent_command == "serve":
        return _cmd_serve(args)
    elif args.agent_command == "inspect":
        result, rc = _cmd_inspect(args)
    elif args.agent_command == "plan":
        result, rc = _cmd_plan(args)
    elif args.agent_command == "cycle":
        result, rc = _cmd_cycle(args)
    elif args.agent_command == "run":
        result, rc = _cmd_run(args)
    elif args.agent_command == "propose":
        result, rc = _cmd_propose(args)
    elif args.agent_command == "apply":
        result, rc = _cmd_apply(args)
    else:
        result, rc = _envelope("agent", ok=False, error="unknown agent command",
                               returncode=2), 2
    _print(result, bool(getattr(args, "pretty", False)))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())


# Keep the historical private names importable for callers/tests while making
# the shared protocol module the single implementation used at runtime.
_compact_contract = _protocol.compact_contract
_compact_skill_list = _protocol.compact_skill_list
_compact_state = _protocol.compact_state
_state_counts = _protocol.state_counts
_snapshot_flat = _protocol.snapshot_flat
_filter_snapshot = _protocol.filter_snapshot
_proposals = _protocol.proposals
