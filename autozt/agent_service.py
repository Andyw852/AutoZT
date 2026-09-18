# -*- coding: utf-8 -*-
"""Transport-independent agent service for AutoZT.

The workflow engine remains usable without this module.  MCP and ``autozt
agent`` call these functions with their own read/execute callbacks; neither
transport owns planning policy or state-shaping logic.
"""

from __future__ import annotations

import inspect as _inspect
from typing import Any, Callable, Dict, Mapping, Optional, Tuple

from autozt import agent_protocol as protocol


LoadState = Callable[[Mapping[str, Any]], Tuple[Optional[Dict[str, Any]], Optional[str], int]]
ExecuteActions = Callable[[Any, bool, Optional[str]], Dict[str, Any]]


def _execute(execute_actions: ExecuteActions, actions: Any, dry_run: bool,
             expected_cursor: Optional[str]) -> Dict[str, Any]:
    """Call the mutation callback with the observed cursor.

    Keep compatibility with third-party two-argument callbacks while all
    built-in transports use the three-argument CAS-aware form.
    """
    try:
        parameters = _inspect.signature(execute_actions).parameters
        accepts_cursor = any(
            p.kind == _inspect.Parameter.VAR_POSITIONAL
            for p in parameters.values()
        ) or len(parameters) >= 3
    except (TypeError, ValueError):
        accepts_cursor = True
    if accepts_cursor:
        return execute_actions(actions, dry_run, expected_cursor)
    return execute_actions(actions, dry_run)  # type: ignore[misc]


def _bounded_max(value: Any) -> int:
    try:
        return max(1, min(int(value), protocol.MAX_ACTIONS))
    except (TypeError, ValueError):
        return protocol.MAX_ACTIONS


def inspect(
    load_state: LoadState,
    scope: Optional[Mapping[str, Any]] = None,
    *,
    view: str = "attention",
    include_monitoring: bool = False,
    max_actions: int = protocol.MAX_ACTIONS,
    include_retry: bool = True,
) -> Tuple[Optional[Dict[str, Any]], Optional[str], int]:
    """Collect one compact observation and derive deterministic candidates."""
    scope = dict(scope or {})
    compact, error, rc = load_state(scope)
    if error:
        return None, error, rc
    compact = compact or {"types": [], "queue": {}}
    proposal_limit = _bounded_max(max_actions)
    candidates = protocol.proposals(compact, include_monitoring,
                                    max_items=proposal_limit)
    actions = protocol.executable_actions(
        candidates,
        max_actions=_bounded_max(max_actions),
        include_retry=include_retry,
    )
    data = {
        "cursor": protocol.state_cursor(compact),
        "scope": {key: value for key, value in scope.items()
                  if value not in (None, "")},
        "counts": protocol.state_counts(compact),
        "snapshot": protocol.filter_snapshot(compact, view or "attention"),
        "proposals": candidates,
        "proposal_limit": proposal_limit,
        "actions": actions,
        "policy": {
            "allowed_actions": sorted(protocol.ALLOWED_ACTIONS),
            "cycle_default_actions": sorted(protocol.ALLOWED_ACTIONS - {"retry_step"}),
            "retry": {"allowed": True, "default": "review", "opt_in": "include_retry"},
            "destructive_actions": "human_review",
            "default_execution": "dry_run",
        },
    }
    return data, None, 0


def plan(
    load_state: LoadState,
    scope: Optional[Mapping[str, Any]] = None,
    *,
    view: str = "attention",
    include_monitoring: bool = False,
    max_actions: int = protocol.MAX_ACTIONS,
) -> Tuple[Optional[Dict[str, Any]], Optional[str], int]:
    data, error, rc = inspect(
        load_state,
        scope,
        view=view,
        include_monitoring=include_monitoring,
        max_actions=max_actions,
        include_retry=True,
    )
    if error:
        return None, error, rc
    assert data is not None
    return {
        "cursor": data["cursor"],
        "scope": data.get("scope") or {},
        "counts": data["counts"],
        "proposals": data["proposals"],
        "actions": data["actions"],
        "snapshot": data["snapshot"],
        "execution": "dry_run_first",
    }, None, 0


def cycle(
    load_state: LoadState,
    execute_actions: ExecuteActions,
    scope: Optional[Mapping[str, Any]] = None,
    *,
    view: str = "attention",
    include_monitoring: bool = False,
    max_actions: int = protocol.MAX_ACTIONS,
    execute: bool = False,
    include_retry: bool = False,
) -> Tuple[Optional[Dict[str, Any]], Optional[str], int]:
    """Run one observe-plan cycle, dry-running unless explicitly requested."""
    data, error, rc = inspect(
        load_state,
        scope,
        view=view,
        include_monitoring=include_monitoring,
        max_actions=max_actions,
        include_retry=include_retry,
    )
    if error:
        return None, error, rc
    assert data is not None
    if execute:
        # A queue transition or another worker may change the state between
        # inspection and the first mutation. Re-read before crossing the
        # execution boundary so a one-shot cycle cannot apply stale actions.
        fresh, fresh_error, fresh_rc = load_state(dict(scope or {}))
        if fresh_error:
            return data, fresh_error, fresh_rc
        fresh_cursor = protocol.state_cursor(fresh or {"types": [], "queue": {}})
        if fresh_cursor != data.get("cursor"):
            data = dict(data)
            data["mode"] = "execute"
            data["execution"] = {
                "results": [], "actions": data.get("actions") or [],
                "failed": True, "stale_plan": True,
                "expected_cursor": data.get("cursor"),
                "current_cursor": fresh_cursor,
                "error": "state changed before execution",
            }
            return data, "state changed before execution; collect a new plan", 1
    execution = _execute(execute_actions, data["actions"], not execute,
                         data.get("cursor"))
    failed = bool(execution.get("failed"))
    data = dict(data)
    data["execution"] = execution
    data["mode"] = "execute" if execute else "dry_run"
    return data, execution.get("error") if failed else None, 1 if failed else 0


def run(
    load_state: LoadState,
    execute_actions: ExecuteActions,
    scope: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> Tuple[Optional[Dict[str, Any]], Optional[str], int]:
    """Stable name for the deterministic observe/plan/execute entrypoint.

    ``cycle`` remains available for MCP compatibility; ``run`` is the clearer
    entrypoint for JSON-only callers and has identical safety semantics.
    """
    return cycle(load_state, execute_actions, scope, **kwargs)


def capabilities() -> Dict[str, Any]:
    """Stable discovery document for JSON-only callers.

    Skills remain data discovered from ``skill.yaml``; this document describes
    the fixed agent protocol and its safety boundary.
    """
    read_commands = ["capabilities", "schema", "skills", "contract", "snapshot", "inspect",
                     "plan", "evidence", "propose", "research_plan", "preflight", "results"]
    mutate_commands = ["cycle", "run", "apply"]
    commands = read_commands + mutate_commands + ["request"]
    request_ops = ["capabilities", "schema", "skills", "contract", "snapshot", "inspect",
                   "plan", "cycle", "run", "evidence", "propose", "research_plan",
                   "preflight", "results", "apply"]
    request_schema = _request_schema(request_ops)
    return {
        "schema_version": protocol.PROTOCOL_VERSION,
        "entrypoints": {
            "one_shot": "autozt agent request -",
            "jsonl": "autozt agent serve",
            "mcp": "autozt mcp",
            "ordinary_cli": "autozt <command> (works without an LLM)",
        },
        "transport": ["argv-json", "stdin-json", "jsonl-stdio", "mcp-stdio"],
        "commands": commands,
        "read": {
            "commands": read_commands,
            "state_source": "autozt list --json",
            "steady_state": "autozt summary --diff",
        },
        "mutate": {"commands": mutate_commands, "default": "dry_run"},
        "conversation": {
            "schema_version": "autozt/conversation/1",
            "sequence": ["research_plan", "preflight", "inspect",
                          "cycle(dry_run)", "user_confirmation", "cycle(execute)", "results"],
            "confirmation_boundary": "execution actions may submit jobs; read-only planning and preflight do not",
        },
        "cursor": {"scope": "session", "restart": "cursor_reset",
                    "usage": "pass the cursor returned by inspect/get_snapshot to cycle or apply"},
        "actions": {
            "allowed": sorted(protocol.ALLOWED_ACTIONS),
            "default": "dry_run",
            "retry": "explicit_review",
            "destructive": ["stop", "rerun", "clean", "conf_set", "hpc"],
            "gateway": "autozt act",
        },
        "request": {
            "ops": request_ops,
            "schema": request_schema,
            "output": "one JSON envelope per request; serve uses one request/response per line",
        },
        "schema_command": "autozt agent schema",
    }


def _request_schema(request_ops: Any = None) -> Dict[str, Any]:
    """Return the one JSON Schema shared by capabilities and ``schema``."""
    ops = list(request_ops or ["capabilities", "schema", "skills", "contract", "snapshot",
                               "inspect", "plan", "cycle", "run", "evidence", "propose",
                               "research_plan", "preflight", "results", "apply"])
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "AutoZT agent request", "type": "object", "required": ["op"],
        "additionalProperties": False,
        "properties": {
            "op": {"type": "string", "enum": ops},
            "scope": {"type": "object", "additionalProperties": False,
                      "properties": {key: {"type": "string"}
                                     for key in ("tt", "material", "status", "step")}},
            "view": {"type": "string", "enum": ["attention", "active", "all"]},
            "include_monitoring": {"type": "boolean"},
            "include_retry": {"type": "boolean"},
            "execute": {"type": "boolean", "default": False},
            "dry_run": {"type": "boolean", "default": True},
            "max_actions": {"type": "integer", "minimum": 1,
                            "maximum": protocol.MAX_ACTIONS},
            "cursor": {"type": "string"}, "full": {"type": "boolean"},
            "skill": {"type": "string"}, "plan": {},
            "actions": {"type": "array", "maxItems": protocol.MAX_ACTIONS},
            "goal": {"type": "string"}, "result_dir": {"type": "string"},
            "property": {"type": "string"}, "direction": {"type": "string"},
            "dimension": {"type": "string"}, "thickness": {"type": "number"},
            "temperature": {}, "carrier": {},
        },
    }


def schema() -> Dict[str, Any]:
    """Complete, model-readable protocol contract with examples and risks."""
    request_ops = ["capabilities", "schema", "skills", "contract", "snapshot", "inspect",
                   "plan", "cycle", "run", "evidence", "propose", "research_plan",
                   "preflight", "results", "apply"]
    action_schema = {
        "type": "object", "additionalProperties": False, "required": ["action"],
        "properties": {
            "action": {"type": "string", "enum": sorted(protocol.ALLOWED_ACTIONS)},
            "tt": {"type": "string"}, "material": {"type": "string"},
            "step": {"type": "string"},
        },
    }
    operations = {
        "capabilities": {"class": "read", "state": False},
        "schema": {"class": "read", "state": False},
        "skills": {"class": "read", "state": False},
        "contract": {"class": "read", "state": False},
        "snapshot": {"class": "read", "state": True, "cursor": True},
        "inspect": {"class": "read", "state": True, "cursor": True},
        "plan": {"class": "read", "state": True, "cursor": True},
        "evidence": {"class": "read", "state": True},
        "propose": {"class": "read", "state": True, "cursor": True},
        "research_plan": {"class": "read", "state": False},
        "preflight": {"class": "read", "state": False},
        "results": {"class": "read", "state": False},
        "cycle": {"class": "mutate", "default": "dry_run", "execute": True},
        "run": {"class": "mutate", "default": "dry_run", "execute": True,
                "description": "确定性 observe -> plan -> optional execute 闭环"},
        "apply": {"class": "mutate", "default": "dry_run", "cursor": True},
    }
    return {
        "schema_version": protocol.PROTOCOL_VERSION,
        "kind": "autozt-agent-protocol",
        "purpose": "让 LLM 或脚本以稳定 JSON 调用 AutoZT；技能是数据契约，不是工具集合",
        "transport": {"one_shot": "autozt agent request -", "jsonl": "autozt agent serve",
                      "mcp": "autozt mcp"},
        "request_schema": _request_schema(request_ops),
        "response_envelope": {
            "type": "object", "required": ["schema_version", "ok", "command", "returncode"],
            "properties": {"schema_version": {"const": protocol.PROTOCOL_VERSION},
                           "ok": {"type": "boolean"}, "command": {"type": "string"},
                           "returncode": {"type": "integer"}, "data": {},
                           "error": {"type": "string"}},
        },
        "operations": operations,
        "conversation": {
            "schema_version": "autozt/conversation/1",
            "states": ["needs_user_input", "awaiting_confirmation", "preflight_review",
                        "ready_to_execute", "executing", "completed", "blocked"],
            "next_actions": ["provide_inputs", "confirm_plan", "review_preflight",
                              "confirm_execution", "review_result", "inspect_results"],
            "sequence": ["research_plan", "preflight", "inspect", "cycle(dry_run)",
                          "user_confirmation", "cycle(execute)", "results"],
        },
        "actions": {"schema": action_schema,
                    "automatic": sorted(protocol.ALLOWED_ACTIONS - {"retry_step"}),
                    "retry": {"name": "retry_step", "requires": "include_retry=true or explicit selection",
                              "keeps_products": True},
                    "human_only": ["stop", "rerun", "clean", "conf_set", "hpc"],
                    "gateway": "autozt act"},
        "state": {"source": "autozt list --json", "steady_state": "autozt summary --diff",
                  "views": ["attention", "active", "all"],
                  "cursor": "sha256 of the compact scoped state; stale plans are rejected"},
        "examples": {"discover": {"op": "schema"},
                     "inspect_failures": {"op": "inspect", "scope": {"status": "error"}},
                     "dry_run": {"op": "run", "scope": {"tt": "band-dft-cpu"}, "execute": False},
                     "execute_safe": {"op": "run", "scope": {"material": "C24/qHPC24"},
                                      "execute": True, "include_retry": False}},
        "token_policy": {"no_change": "Use summary --diff; do not call an LLM",
                          "polling": "Pass cursor to snapshot/get_snapshot and request only changes",
                          "large_scope": "Use tt/material/status filters and max_actions"},
    }
