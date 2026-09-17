# -*- coding: utf-8 -*-
"""autozt mcp —— 把 AutoZT 以 MCP（Model Context Protocol）stdio 服务暴露给上层 agent。

设计三条硬约束（决定了"新技能仍然方便"和"agent 不能乱来"）：
  1. **只调 CLI，不 import 内部函数**：门禁/审计/档案只有一份事实来源，不可能与命令行漂移。
  2. **工具是通用动词，不认识任何具体技能**：技能清单/步骤清单都从 skill.yaml 动态读，
     所以加新技能 = 加目录，MCP 侧零改动、工具数量恒定（上限 21，超了就说明有人按技能加工具）。
  3. **风险在协议边界收口**：只读工具直连 CLI；变更/破坏性工具一律走 autozt act
     （风险分级 + 人工批准 + 审计日志），approve 只从真 TTY 生效，agent 无法自我批准。

协议：JSON-RPC 2.0 over stdio，支持 initialize / tools/list / tools/call，
以及只读 resources/list、resources/read、prompts/list、prompts/get。
"""
import glob
import json
import os
import subprocess
import sys
from collections import OrderedDict

from autozt import agent_protocol as _protocol
from autozt import agent_service as _service

PROTOCOL_VERSION = "2025-06-18"
SUPPORTED_PROTOCOL_VERSIONS = {"2024-11-05", "2025-06-18"}
SCHEMA_VERSION = "2"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "autozt")

# 通用动词工具表：读/写/破坏三类，与技能无关
TOOLS = [
    # ---- 只读 ----
    ("list_skills", "列出所有技能（skill.yaml 汇总：步骤数/默认集群/版本）",
     {"type": "object", "properties": {}}, "read", ["skills"]),
    ("list_materials", "列出材料及其步骤状态（只读，不提交、不拉文件）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "status": {"type": "string"}}}, "read", ["summary"]),
    ("get_summary", "极简汇总（每类型一行计数 + FAIL 清单 + 全局队列）",
     {"type": "object", "properties": {"tt": {"type": "string"}}}, "read", ["summary"]),
    ("get_status", "单材料详情（含每步诊断、hpc、work_dir 来源）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"}},
      "required": ["material"]}, "read", ["status"]),
    ("capabilities", "读取固定 agent 协议、动作和风险边界（只读）",
     {"type": "object", "properties": {}, "additionalProperties": False},
     "read", ["capabilities"]),
    ("schema", "读取完整 JSON 请求/响应/动作 Schema（只读，首次接入调用一次）",
     {"type": "object", "properties": {}, "additionalProperties": False},
     "read", ["schema"]),
    ("conf_get", "查看某步 step.conf 合并后的最终值（只读）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"}},
      "required": ["material"]}, "read", ["conf"]),
    ("describe_skill", "读取技能能力目录：步骤、依赖、判据和专属说明（只读）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                           "detail": {"type": "string",
                                                      "enum": ["compact", "full"]}},
      "required": ["tt"], "additionalProperties": False}, "read", ["schema"]),
    ("probe_step", "自适应步骤诊断：缺陷技能作业探测，其余技能读取通用诊断（只读）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"}},
      "required": ["material"], "additionalProperties": False}, "read", ["probe"]),
    # ---- 变更（走 act 网关）----
    ("start_step", "推进材料/步骤（输入没生成先 gen 再提交）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"}},
      "required": ["material"]}, "mutate", ["start"]),
    ("retry_step", "保留产物重新生成输入并重交（不删产物）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"}},
      "required": ["material"]}, "mutate", ["retry"]),
    ("fetch_results", "强制把已完成步骤的结果拉回本地",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"}},
      "required": ["material"]}, "mutate", ["fetch"]),
    ("advance_ready", "按依赖图推进所有就绪步骤（等价于 autozt auto on）",
     {"type": "object", "properties": {"tt": {"type": "string"}}},
     "mutate", ["auto", "on"]),
    ("inspect", "一次返回关注状态、失败诊断和确定性候选动作（只读）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                           "material": {"type": "string"},
                                           "status": {"type": "string"},
                                           "view": {"type": "string",
                                                     "enum": ["attention", "active", "all"]},
                                           "include_monitoring": {"type": "boolean"},
                                           "max_actions": {"type": "integer",
                                                           "minimum": 1, "maximum": 20}},
      "additionalProperties": False}, "read", ["inspect"]),
    # ---- 破坏性（走 act 网关，默认拒绝，需人工 approve）----
    ("stop_step", "取消该步骤的作业（破坏性）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"}},
      "required": ["material"]}, "destructive", ["stop"]),
    ("rerun_step", "删目录重新生成（破坏性：会删除该步全部产物）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"}},
      "required": ["material"]}, "destructive", ["rerun"]),
    ("clean_material", "删除生成物回到 PREP（破坏性）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"}},
      "required": ["material"]}, "destructive", ["clean"]),
    # ---- 可复现性 ----
    ("cycle", "观察、规则规划，并可选执行一轮非破坏性动作",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                           "material": {"type": "string"},
                                           "status": {"type": "string"},
                                           "view": {"type": "string",
                                                     "enum": ["attention", "active", "all"]},
                                           "include_monitoring": {"type": "boolean"},
                                           "max_actions": {"type": "integer",
                                                           "minimum": 1, "maximum": 20},
                                           "include_retry": {"type": "boolean"},
                                           "execute": {"type": "boolean"}},
      "additionalProperties": False}, "mutate", ["cycle"]),
    # ---- LLM 高层接口：内部仍只调用上面的 CLI 能力 ----
    ("get_snapshot", "紧凑状态快照；传 cursor 时只返回变化，降低上下文开销",
     {"type": "object", "properties": {
         "tt": {"type": "string"}, "material": {"type": "string"},
         "status": {"type": "string"}, "cursor": {"type": "string"},
         "view": {"type": "string", "enum": ["attention", "active", "all"]}},
      "additionalProperties": False}, "read", ["status"]),
    ("propose_actions", "按状态与稳定诊断码生成候选动作；只建议，不执行",
     {"type": "object", "properties": {
         "tt": {"type": "string"}, "material": {"type": "string"},
         "status": {"type": "string"},
         "include_monitoring": {"type": "boolean"},
         "max_actions": {"type": "integer", "minimum": 1, "maximum": 20}},
      "additionalProperties": False}, "read", ["status"]),
    ("apply_actions", "批量执行已选定的非破坏性动作；破坏性动作始终拒绝",
     {"type": "object", "properties": {
         "dry_run": {"type": "boolean"},
         "cursor": {"type": "string"},
         "tt": {"type": "string"}, "material": {"type": "string"},
         "status": {"type": "string"},
         "actions": {"type": "array", "minItems": 1, "maxItems": 20,
                     "items": {"type": "object", "properties": {
                         "action": {"type": "string", "enum": [
                             "start_step", "retry_step", "fetch_results", "advance_ready"]},
                         "tt": {"type": "string"}, "material": {"type": "string"},
                         "step": {"type": "string"}},
                         "required": ["action"], "additionalProperties": False}}},
      "required": ["actions"], "additionalProperties": False},
     "mutate", ["batch"]),
]
BY_NAME = {t[0]: t for t in TOOLS}

# 默认 full 暴露完整的固定工具表；新 LLM 接入推荐 workflow，只把高层闭环放进模型上下文。
COMPACT_TOOLS = {
    # Compatibility profile retained for existing clients.
    "schema", "list_skills", "describe_skill", "get_snapshot", "probe_step",
    "propose_actions", "apply_actions",
}
MONITOR_TOOLS = {
    # Smallest profile for a long-running monitor. It can observe a cursor and
    # run the already-safe cycle, while discovery and triage stay opt-in.
    "get_snapshot", "cycle",
}
WORKFLOW_TOOLS = {
    # One-shot observation/planning/execution for new LLM integrations.
    "schema", "capabilities", "list_skills", "describe_skill", "get_snapshot", "inspect",
    "probe_step", "cycle", "apply_actions",
}

TOOL_RESULT_SCHEMA = {
    "type": "object",
    "properties": {
        "ok": {"type": "boolean"},
        "risk": {"type": "string", "enum": ["read", "mutate", "destructive"]},
        "returncode": {"type": "integer"},
        "schema_version": {"type": "string"},
        "data": {},
        "error": {"type": "string"},
    },
    "required": ["ok", "risk", "returncode", "schema_version"],
}

_SNAPSHOTS = OrderedDict()
_SNAPSHOT_LIMIT = 128


def _cfg_args():
    cfg = (os.environ.get("AUTOZT_CONFIG")
           or os.environ.get("AUTOZT_CFG") or "").strip()
    return ["-c", cfg] if cfg else []


def _run(argv, timeout=1800):
    """只通过 CLI 执行——门禁/审计/档案全部沿用命令行那一套。"""
    cmd = [sys.executable, PROG] + _cfg_args() + argv
    env = dict(os.environ)
    env.setdefault("AUTOZT_ACTOR", "mcp")   # 让 act 网关认出这是 agent 会话
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return 1, "", "timeout after %ss" % timeout
    return p.returncode, p.stdout or "", p.stderr or ""


def _mat_args(a):
    out = []
    if a.get("tt"):
        out += ["-tt", str(a["tt"])]
    if a.get("material"):
        out += ["-p", str(a["material"])]
    if a.get("step"):
        out += ["-j", str(a["step"])]
    return out


def _validate_args(name, args):
    """Validate the small JSON-Schema subset used by this server.

    MCP clients normally validate tool inputs, but the protocol boundary must also reject
    malformed direct calls. Keep this dependency-free so AutoZT's core remains zero-deps.
    """
    if not isinstance(args, dict):
        return "arguments must be an object"
    schema = BY_NAME[name][2]
    props = schema.get("properties") or {}
    for key in schema.get("required") or []:
        if key not in args or args[key] is None or args[key] == "":
            return "missing required argument: %s" % key
    if schema.get("additionalProperties") is False:
        unknown = sorted(set(args) - set(props))
        if unknown:
            return "unknown argument(s): %s" % ", ".join(unknown)
    pytypes = {"string": str, "boolean": bool, "integer": int,
               "array": list, "object": dict}
    for key, value in args.items():
        spec = props.get(key) or {}
        typ = spec.get("type")
        if typ == "integer" and isinstance(value, bool):
            return "argument %s must be integer" % key
        if typ in pytypes and not isinstance(value, pytypes[typ]):
            return "argument %s must be %s" % (key, typ)
        if typ == "integer" and "minimum" in spec and value < spec["minimum"]:
            return "argument %s is below minimum %s" % (key, spec["minimum"])
        if typ == "integer" and "maximum" in spec and value > spec["maximum"]:
            return "argument %s is above maximum %s" % (key, spec["maximum"])
        if "enum" in spec and value not in spec["enum"]:
            return "argument %s must be one of: %s" % (
                key, ", ".join(str(x) for x in spec["enum"]))
        if typ == "array":
            if len(value) < spec.get("minItems", 0):
                return "argument %s has too few items" % key
            if len(value) > spec.get("maxItems", len(value)):
                return "argument %s has too many items" % key
            item_spec = spec.get("items") or {}
            for i, item in enumerate(value):
                if item_spec.get("type") == "object" and not isinstance(item, dict):
                    return "argument %s[%d] must be object" % (key, i)
                for req in item_spec.get("required") or []:
                    if req not in item or item[req] in (None, ""):
                        return "argument %s[%d] missing %s" % (key, i, req)
                iprops = item_spec.get("properties") or {}
                if item_spec.get("additionalProperties") is False:
                    extra = sorted(set(item) - set(iprops))
                    if extra:
                        return "argument %s[%d] has unknown field(s): %s" % (
                            key, i, ", ".join(extra))
                for ikey, ivalue in item.items():
                    enum = (iprops.get(ikey) or {}).get("enum")
                    if enum and ivalue not in enum:
                        suffix = " (destructive or unsupported)" if (
                            name == "apply_actions" and ikey == "action") else ""
                        return "argument %s[%d].%s is not allowed%s" % (
                            key, i, ikey, suffix)
    return None


def _result(risk, rc=0, data=None, error=None, text=None):
    """Build an MCP result with identical structured and text representations.

    structuredContent avoids reparsing CLI prose. The compact JSON text is retained for
    older MCP clients, as recommended by the MCP tools specification.
    """
    payload = {"ok": rc == 0 and not error, "risk": risk, "returncode": int(rc),
               "schema_version": SCHEMA_VERSION}
    if data is not None:
        payload["data"] = data
    if error:
        payload["error"] = str(error)
    profile = (os.environ.get("AUTOZT_MCP_PROFILE") or "").strip().lower()
    if text is not None:
        rendered = text
    elif data is not None and profile in ("compact", "monitor", "workflow", "agent"):
        # Modern MCP clients consume structuredContent directly. Avoid sending the same
        # potentially large snapshot twice in the model-facing text channel.
        brief = {k: payload[k] for k in ("ok", "risk", "returncode", "schema_version")}
        brief["data_in_structured_content"] = True
        rendered = json.dumps(brief, ensure_ascii=False, separators=(",", ":"))
    else:
        rendered = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    return {"isError": not payload["ok"],
            "content": [{"type": "text", "text": rendered[:20000]}],
            "structuredContent": payload}


def _json_stdout(out):
    try:
        return json.loads((out or "").strip()), None
    except (TypeError, ValueError) as e:
        return None, "AutoZT did not return valid JSON: %s" % e


def _compact_contract(contract):
    """Keep the skill contract useful while dropping generator paths and prose bulk."""
    if not isinstance(contract, dict):
        return contract
    out = {k: contract.get(k) for k in
           ("skill", "version", "desc", "schema", "corrections", "issues")
           if contract.get(k) not in (None, [], "")}
    steps = []
    for step in contract.get("steps") or []:
        if isinstance(step, dict):
            steps.append({k: step.get(k) for k in
                          ("name", "label", "seq", "check", "run", "needs", "produces")
                          if step.get(k) not in (None, [], "")})
    if steps:
        out["steps"] = steps
    io = contract.get("io_schema") or {}
    if io:
        out["io_schema"] = {k: io.get(k) for k in ("inputs", "outputs", "params", "notes")
                            if io.get(k) not in (None, [], "")}
    flow = contract.get("flow") or {}
    if flow:
        out["flow"] = {k: flow.get(k) for k in ("summary", "stages", "requires")
                        if flow.get(k) not in (None, [], "")}
    stats = contract.get("stats")
    if stats:
        out["stats"] = stats
    missing = [key for key in ("io_schema", "flow", "corrections")
               if not contract.get(key)]
    out["contract_status"] = {"complete": not missing, "missing": missing}
    if missing and steps:
        # A new skill may start with only executable steps. Expose a clearly marked
        # structural fallback without inventing physical inputs/outputs or policies.
        out["inferred"] = {"step_graph": [
            {k: step.get(k) for k in ("label", "seq", "needs", "check")
             if step.get(k) not in (None, [], "")}
            for step in steps
        ]}
    return out


def _state_argv(args):
    argv = _mat_args(args)
    if args.get("status"):
        argv += ["-status", str(args["status"])]
    # ``list`` is AutoZT's explicitly read-only collector. Keep the model-facing
    # observation path free of auto-fetch/auto-advance semantics even if the CLI's
    # human-oriented ``status`` command changes in a future release.
    argv += ["list", "--json"]
    return argv


def _load_state(args):
    rc, out, err = _run(_state_argv(args))
    if rc != 0:
        return rc, None, ((out or "") + (("\n" + err) if err else "")).strip()
    data, parse_error = _json_stdout(out)
    return (1, None, parse_error) if parse_error else (0, data, None)


def _active_label(material):
    active = material.get("active")
    if isinstance(active, dict):
        return active.get("label") or active.get("name")
    return active


def _split_scope(value):
    return [x.strip() for x in str(value or "").split(",") if x.strip()]


def _status_kinds(spec):
    aliases = {
        "done": "OK", "ok": "OK", "running": "R", "run": "R",
        "pd": "PD", "queued": "PD", "error": "FAIL", "fail": "FAIL",
        "waiting": "WAIT", "wait": "WAIT", "scancel": "SCANCEL",
        "scancelled": "SCANCEL", "cancelled": "SCANCEL", "todo": "TODO",
        "prep": "PREP", "other": "OTHER",
    }
    kinds = set()
    for raw in _split_scope(spec):
        key = raw.lower()
        if key == "ready":
            kinds.update(("TODO", "PREP"))
        else:
            kinds.add(aliases.get(key, raw.upper()))
    return kinds


def _material_matches(name, wanted):
    if not wanted:
        return True
    text = str(name or "")
    base = os.path.basename(text.rstrip("/"))
    return any(text == item or base == item or text.endswith("/" + item)
               for item in wanted)


def _scope_data(data, args=None):
    """Apply model-requested scope after CLI collection.

    The human CLI uses ``-p`` mainly to narrow collection and table rendering; JSON
    branches intentionally preserve the full schema. MCP must enforce its own scope
    so a material-specific call cannot leak an entire project into the context.
    """
    args = args or {}
    if not isinstance(data, dict) or not isinstance(data.get("types"), list):
        return data
    wanted_tt = set(_split_scope(args.get("tt")))
    wanted_mat = _split_scope(args.get("material"))
    wanted_status = _status_kinds(args.get("status"))
    out = dict(data)
    types = []
    for task_type in data.get("types") or []:
        if wanted_tt and task_type.get("key") not in wanted_tt:
            continue
        ct = dict(task_type)
        mats = []
        for material in task_type.get("materials") or []:
            if not _material_matches(material.get("name"), wanted_mat):
                continue
            steps = material.get("steps") or []
            if wanted_status and not any(s.get("kind") in wanted_status for s in steps):
                continue
            mats.append(material)
        ct["materials"] = mats
        types.append(ct)
    out["types"] = types
    return out


def _compact_state(data, args=None):
    """Reduce full status JSON to fields needed for planning and state comparison."""
    data = _scope_data(data, args)
    out = {"types": [], "queue": data.get("queue") or {}}
    for task_type in data.get("types") or []:
        compact_type = {"key": task_type.get("key"), "materials": []}
        for material in task_type.get("materials") or []:
            cm = {"name": material.get("name"),
                  "hpc": material.get("hpc_name") or material.get("hpc") or "",
                  "dim": material.get("dim") or "",
                  "active": _active_label(material),
                  "action": material.get("action") or "", "steps": []}
            for step in material.get("steps") or []:
                job = step.get("job") or {}
                cs = {"label": step.get("label"), "kind": step.get("kind")}
                for key in ("diag", "diag_code", "suggested_action", "action_reason"):
                    if step.get(key):
                        cs[key] = step[key]
                if job:
                    cs["job"] = {k: job.get(k) for k in ("id", "state", "info")
                                 if job.get(k) not in (None, "")}
                cm["steps"].append(cs)
            compact_type["materials"].append(cm)
        out["types"].append(compact_type)
    return out


def _state_counts(compact):
    counts = {}
    for task_type in compact.get("types") or []:
        for material in task_type.get("materials") or []:
            for step in material.get("steps") or []:
                kind = step.get("kind") or "unknown"
                counts[kind] = counts.get(kind, 0) + 1
    return counts


def _snapshot_flat(compact):
    out = {}
    for task_type in compact.get("types") or []:
        tt = task_type.get("key") or ""
        for material in task_type.get("materials") or []:
            name = material.get("name") or ""
            meta_key = "%s\x1f%s\x1f__material__" % (tt, name)
            out[meta_key] = {"type": tt, "material": name, "step": "__material__",
                             "active": material.get("active") or "",
                             "hpc": material.get("hpc") or "",
                             "dim": material.get("dim") or "",
                             "action": material.get("action") or ""}
            for step in material.get("steps") or []:
                key = "%s\x1f%s\x1f%s" % (tt, name, step.get("label") or "")
                job = step.get("job") or {}
                out[key] = {"type": tt, "material": name, "step": step.get("label"),
                            "kind": step.get("kind"),
                            "diag_code": step.get("diag_code") or "",
                            "queue_reason": job.get("info") or ""}
    return out


def _remember_snapshot(compact):
    cursor = _protocol.state_cursor(compact)
    _SNAPSHOTS[cursor] = compact
    _SNAPSHOTS.move_to_end(cursor)
    while len(_SNAPSHOTS) > _SNAPSHOT_LIMIT:
        _SNAPSHOTS.popitem(last=False)
    return cursor


def _filter_snapshot(compact, view):
    if view == "all":
        return compact
    out = {"types": [], "queue": compact.get("queue") or {}}
    for task_type in compact.get("types") or []:
        ct = {"key": task_type.get("key"), "materials": []}
        for material in task_type.get("materials") or []:
            active = material.get("active")
            if view == "active":
                steps = [s for s in material.get("steps") or []
                         if s.get("label") == active or s.get("kind") == "FAIL"]
            else:
                steps = [s for s in material.get("steps") or []
                         if s.get("kind") not in ("OK", "WAIT")
                         or s.get("label") == active]
            if steps:
                cm = dict(material)
                cm["steps"] = steps
                ct["materials"].append(cm)
        if ct["materials"]:
            out["types"].append(ct)
    return out


def _get_snapshot(args):
    rc, data, error = _load_state(args)
    if rc:
        return _result("read", rc=rc, error=error)
    compact = _compact_state(data, args)
    new_cursor = _remember_snapshot(compact)
    old_cursor = args.get("cursor")
    if old_cursor:
        old = _SNAPSHOTS.get(old_cursor)
        if old is None:
            payload = {"cursor": new_cursor, "cursor_reset": True,
                       "counts": _state_counts(compact),
                       "snapshot": _filter_snapshot(compact, args.get("view") or "attention")}
        else:
            changes = _protocol.snapshot_changes(old, compact)
            payload = {"cursor": new_cursor, "unchanged": not changes,
                       "counts": _state_counts(compact), "changes": changes}
    else:
        payload = {"cursor": new_cursor, "first_snapshot": True,
                   "counts": _state_counts(compact),
                   "snapshot": _filter_snapshot(compact, args.get("view") or "attention")}
    return _result("read", data=payload)


def _propose_actions(args):
    rc, data, error = _load_state(args)
    if rc:
        return _result("read", rc=rc, error=error)
    compact = _compact_state(data, args)
    include_monitoring = bool(args.get("include_monitoring"))
    limit = args.get("max_actions", _protocol.MAX_ACTIONS)
    try:
        proposal_limit = max(1, min(int(limit), _protocol.MAX_ACTIONS))
    except (TypeError, ValueError):
        proposal_limit = _protocol.MAX_ACTIONS
    proposals = _protocol.proposals(compact, include_monitoring,
                                    max_items=proposal_limit)
    return _result("read", data={"cursor": _protocol.state_cursor(compact),
                                 "scope": {key: args.get(key) for key in
                                            ("tt", "material", "status")
                                            if args.get(key) not in (None, "")},
                                 "counts": _state_counts(compact),
                                 "proposals": proposals,
                                 "proposal_limit": proposal_limit})


def _inspect(args):
    """Single read call for an LLM: state, bounded attention view, and plan."""
    data, error, rc = _service.inspect(
        _mcp_load_compact, args,
        view=args.get("view") or "attention",
        include_monitoring=bool(args.get("include_monitoring")),
        max_actions=args.get("max_actions", _protocol.MAX_ACTIONS),
        include_retry=bool(args.get("include_retry", True)))
    if error:
        return _result("read", rc=rc, error=error)
    return _result("read", data=data)


def _cycle(args):
    """Observe and plan once; execute only when the caller explicitly sets execute."""
    execute = bool(args.get("execute"))

    def apply(actions, dry_run):
        result = _apply_actions({"dry_run": dry_run, "actions": actions})
        structured = result.get("structuredContent") or {}
        return {
            **(structured.get("data") or {}),
            "failed": bool(result.get("isError")),
            "error": structured.get("error"),
        }

    data, error, rc = _service.cycle(
        _mcp_load_compact, apply, args,
        view=args.get("view") or "attention",
        include_monitoring=bool(args.get("include_monitoring")),
        max_actions=args.get("max_actions", _protocol.MAX_ACTIONS),
        execute=execute,
        include_retry=bool(args.get("include_retry", False)))
    return _result("mutate", rc=rc, data=data, error=error)


def _mcp_load_compact(scope):
    rc, raw, error = _load_state(scope or {})
    if rc:
        return None, error, rc
    return _compact_state(raw, scope or {}), None, 0


def _apply_actions(args):
    allowed = {"start_step", "retry_step", "fetch_results", "advance_ready"}
    actions = args.get("actions") or []
    for i, action in enumerate(actions):
        if action.get("action") not in allowed:
            return _result("mutate", rc=1,
                           error="action[%d] is destructive or unsupported" % i)
        if action.get("action") != "advance_ready" and not action.get("material"):
            return _result("mutate", rc=1,
                           error="action[%d] requires material" % i)
    expected = args.get("cursor")
    if expected:
        compact, error, rc = _mcp_load_compact(args)
        if error:
            return _result("mutate", rc=rc, error=error)
        actual = _protocol.state_cursor(compact or {})
        if str(expected) != actual:
            return _result("mutate", rc=1,
                           data={"stale_plan": True,
                                 "expected_cursor": str(expected),
                                 "current_cursor": actual,
                                 "actions": actions},
                           error="plan is stale; call get_snapshot or inspect again")
    if args.get("dry_run"):
        return _result("mutate", data={"dry_run": True, "actions": actions})
    results = []
    failed = False
    for action in actions:
        name = action["action"]
        call_args = {k: v for k, v in action.items() if k != "action"}
        # The profile controls the public tool surface.  A high-level cycle is
        # allowed to dispatch its bounded primitive actions internally, while
        # those actions still pass through the same validation and act gateway.
        got = call_tool(name, call_args, _internal=True)
        structured = got.get("structuredContent") or {}
        results.append({"action": name, "arguments": call_args,
                        "ok": not got.get("isError"),
                        "returncode": structured.get("returncode"),
                        "data": structured.get("data"),
                        "error": structured.get("error")})
        if got.get("isError"):
            failed = True
            break
    return _result("mutate", rc=1 if failed else 0,
                   data={"results": results, "stopped_on_error": failed},
                   error="batch stopped after first failed action" if failed else None)


def call_tool(name, args, _internal=False):
    """执行一个工具：只读直连，变更/破坏性走 autozt act（会返回需批准的错误）。"""
    if name not in BY_NAME:
        return _result("read", rc=1, error="unknown tool: %s" % name)
    _n, _d, _schema, risk, verb = BY_NAME[name]
    if not _internal and name not in _profile_tool_names():
        return _result(risk, rc=1,
                       error="tool %s is not exposed in the selected MCP profile"
                       % name)
    if readonly_profile() and risk != "read":
        return _result(risk, rc=1,
                       error="tool %s is not exposed in the read-only profile "
                             "(AUTOZT_MCP_READONLY=1)" % name)
    args = args or {}
    validation_error = _validate_args(name, args)
    if validation_error:
        return _result(risk, rc=1, error=validation_error)
    if name == "get_snapshot":
        return _get_snapshot(args)
    if name == "capabilities":
        return _result("read", data=_service.capabilities())
    if name == "schema":
        return _result("read", data=_service.schema())
    if name == "inspect":
        return _inspect(args)
    if name == "cycle":
        return _cycle(args)
    if name == "propose_actions":
        return _propose_actions(args)
    if name == "apply_actions":
        return _apply_actions(args)

    argv = _mat_args(args)
    json_result = False
    if name == "list_skills":
        argv = ["schema", "--json"]
        json_result = True
    elif name == "describe_skill":
        # schema is the machine-readable skill contract, not a terminal table.
        argv = ["schema", str(args["tt"]), "--json"]
        json_result = True
    elif name == "conf_get":
        argv += ["conf"]
    elif name == "conf_set":
        argv += ["conf", "--set", "%s=%s" % (args.get("key"), args.get("value"))]
    elif name == "session_export":
        argv += ["session", "export"]
        if args.get("out"):
            argv += ["--out", str(args["out"])]
    elif name == "probe_step":
        # The legacy probe script is defect-dft-cpu-specific.  Preserve that
        # behaviour for defect jobs and route every other skill to the generic,
        # read-only diagnose command instead of applying a hard-coded defect path.
        argv += ["probe" if args.get("tt") in (None, "defect-dft-cpu")
                 else "diagnose"]
        json_result = True
    elif name == "list_materials":
        if args.get("status"):
            argv += ["-status", str(args["status"])]
        argv += ["summary", "--json"]
        json_result = True
    elif name == "get_summary":
        argv += ["summary", "--json"]
        json_result = True
    elif name == "get_status":
        # Keep detailed reads on the same explicitly read-only collector as snapshots.
        argv += ["list", "--json"]
        json_result = True
    else:
        argv += list(verb)
    if risk != "read":
        # 变更/破坏性：交给动作网关（agent 会话下破坏性动作会被拒并给出批准命令）
        argv = ["act"] + argv
    rc, out, err = _run(argv)
    combined = ((out or "") + (("\n" + err) if err.strip() else "")).strip()
    if rc != 0:
        return _result(risk, rc=rc, error=combined or "AutoZT command failed")
    if json_result:
        parsed, parse_error = _json_stdout(out)
        if parse_error:
            return _result(risk, rc=1, error=parse_error,
                           data={"raw": combined[:4000]})
        if name in ("get_status",):
            parsed = _scope_data(parsed, args)
        if name == "list_skills":
            parsed = _protocol.compact_skill_list(parsed)
        elif name == "describe_skill" and args.get("detail", "compact") != "full":
            # The full contract remains available for explicit requests; compact is the
            # default so a model does not pay for generator paths and long prose every turn.
            if isinstance(parsed, dict) and isinstance(parsed.get("skills"), list):
                parsed = dict(parsed)
                parsed["skills"] = [_protocol.compact_contract(x) for x in parsed["skills"]]
            else:
                parsed = _protocol.compact_contract(parsed)
        return _result(risk, data=parsed)
    return _result(risk, data={"stdout": combined})


def readonly_profile():
    """只读档：AUTOZT_MCP_READONLY 为真时只暴露 read 档工具。

    给"只让 agent 看、不让它动"的场景用：变更/破坏性动词在协议层就不出现，
    比"出现但被拒"更省事，也更难被绕过。
    """
    v = (os.environ.get("AUTOZT_MCP_READONLY") or "").strip().lower()
    return v in ("1", "true", "yes", "on")


def _profile_tool_names():
    # Keep the default model-facing surface small.  The complete compatibility
    # surface remains available explicitly with AUTOZT_MCP_PROFILE=full.
    profile = (os.environ.get("AUTOZT_MCP_PROFILE") or "workflow").strip().lower()
    if profile == "compact":
        return COMPACT_TOOLS
    if profile == "monitor":
        return MONITOR_TOOLS
    if profile in ("workflow", "agent"):
        return WORKFLOW_TOOLS
    return {item[0] for item in TOOLS}


def _tools_list():
    out = []
    only_read = readonly_profile()
    exposed = _profile_tool_names()
    for name, desc, schema, risk, _verb in TOOLS:
        if only_read and risk != "read":
            continue
        if name not in exposed:
            continue
        s = dict(schema)
        out.append({"name": name, "description": "[%s] %s" % (risk, desc),
                    "inputSchema": s,
                    "outputSchema": TOOL_RESULT_SCHEMA,
                    "annotations": {"readOnlyHint": risk == "read",
                                    "destructiveHint": risk == "destructive",
                                    "idempotentHint": risk == "read"},
                    "x-risk": risk})
    return out


def handle(req):
    method = req.get("method")
    rid = req.get("id")
    if method == "initialize":
        requested = ((req.get("params") or {}).get("protocolVersion") or "").strip()
        selected = requested if requested in SUPPORTED_PROTOCOL_VERSIONS else PROTOCOL_VERSION
        result = {"protocolVersion": selected,
                  "capabilities": {"tools": {}, "resources": {"subscribe": False,
                                                                  "listChanged": False},
                                    "prompts": {}},
                  "serverInfo": {"name": "autozt", "version": SCHEMA_VERSION}}
    elif method in ("tools/list", "list_tools"):
        result = {"tools": _tools_list()}
    elif method == "resources/list":
        result = {"resources": _resources_list()}
    elif method == "resources/read":
        params = req.get("params") or {}
        got = _resource_read(params.get("uri"))
        if got is None:
            return {"jsonrpc": "2.0", "id": rid,
                    "error": {"code": -32602, "message": "resource not found: %s"
                              % params.get("uri")}}
        result = got
    elif method == "prompts/list":
        result = {"prompts": _prompts_list()}
    elif method == "prompts/get":
        params = req.get("params") or {}
        got = _prompt_get(params.get("name"), params.get("arguments"))
        if got is None:
            return {"jsonrpc": "2.0", "id": rid,
                    "error": {"code": -32602, "message": "prompt not found: %s"
                              % params.get("name")}}
        result = got
    elif method == "tools/call":
        params = req.get("params") or {}
        result = call_tool(params.get("name"), params.get("arguments") or {})
    elif method in ("notifications/initialized", "initialized"):
        return None
    else:
        return {"jsonrpc": "2.0", "id": rid,
                "error": {"code": -32601, "message": "method not found: %s" % method}}
    return {"jsonrpc": "2.0", "id": rid, "result": result}


def main(argv=None):
    # MCP stdio is a protocol stream, so never inherit a Windows legacy code page.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except AttributeError:
            pass
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--list-tools" in argv:
        print(json.dumps(_tools_list(), ensure_ascii=False, indent=1))
        return 0
    if "--call" in argv:
        i = argv.index("--call")
        name = argv[i + 1]
        raw = argv[i + 2] if len(argv) > i + 2 else "{}"
        print(json.dumps(call_tool(name, json.loads(raw)), ensure_ascii=False))
        return 0
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
        except ValueError:
            continue
        resp = handle(req)
        if resp is None:
            continue
        sys.stdout.write(json.dumps(resp, ensure_ascii=False) + "\n")
        sys.stdout.flush()
    return 0


# ---- MCP resources：把技能与项目文档作为只读资源暴露（agent 直接拿到上下文）----
def _skill_dirs():
    """Return discovered skill directories, including configured extensions.

    The CLI already supports ``skill_paths`` as a drop-in extension point.  Keep
    MCP resources aligned with that discovery path so a new skill's contract and
    README are visible to the model without changing the MCP server code.
    """
    bases = []
    cfg_path = (os.environ.get("AUTOZT_CONFIG")
                or os.environ.get("AUTOZT_CFG") or "").strip()
    try:
        from autozt import load_config

        cfg, found = load_config(cfg_path or None)
        for path in (cfg.get("skill_paths") or []):
            bases.append(os.path.expanduser(str(path)))
        if found:
            cdir = os.path.dirname(os.path.abspath(found))
            bases.extend((os.path.join(cdir, "skill"),
                          os.path.normpath(os.path.join(cdir, "..", "skill"))))
    except (Exception, SystemExit):
        # Resource discovery must never make the MCP handshake fail because a
        # user config is incomplete; the built-in skill pool remains available.
        pass
    bases.append(os.path.join(ROOT, "skill"))

    out, seen = [], set()
    for base in bases:
        for d in sorted(glob.glob(os.path.join(base, "*"))):
            if not os.path.isdir(d):
                continue
            rd = os.path.realpath(d)
            if rd not in seen:
                seen.add(rd)
                out.append(rd)
    return out


def _resources_list():
    out, seen_names = [], set()
    for d in _skill_dirs():
        name = os.path.basename(d)
        if name in seen_names:
            continue
        seen_names.add(name)
        for f in ("skill.yaml", "README.md"):
            p = os.path.join(d, f)
            if os.path.isfile(p):
                out.append({"uri": "autozt://skill/%s/%s" % (name, f),
                            "name": "%s/%s" % (name, f),
                            "mimeType": "text/markdown" if f.endswith(".md")
                            else "application/yaml"})
    for f in ("README.md", "README.en.md", "CHANGELOG.md", "docs/ACCEPTANCE.md",
              "docs/CONFIGURING.md", "docs/mcp.md"):
        if os.path.isfile(os.path.join(ROOT, f)):
            out.append({"uri": "autozt://doc/%s" % f, "name": f,
                        "mimeType": "text/markdown"})
    return out


def _resource_path(uri):
    """Resolve a built-in/configured skill resource or a repository document.

    Skill paths are explicitly configured extension roots; traversal and arbitrary
    files are still rejected so an agent cannot use the resource API as a disk reader.
    """
    if not isinstance(uri, str) or not uri.startswith("autozt://"):
        return None
    rest = uri[len("autozt://"):]
    if rest.startswith("skill/"):
        _name, _file = rest[len("skill/"):].split("/", 1) \
            if "/" in rest[len("skill/"):] else ("", "")
        if not _name or _file not in ("skill.yaml", "README.md"):
            return None
        for skill_dir in _skill_dirs():
            if os.path.basename(skill_dir) != _name:
                continue
            p = os.path.realpath(os.path.join(skill_dir, _file))
            root = os.path.realpath(skill_dir)
            if p.startswith(root + os.sep) and os.path.isfile(p):
                return p
        return None
    elif rest.startswith("doc/"):
        rel = rest[len("doc/"):]
    else:
        return None
    if ".." in rel.split("/") or rel.startswith("/"):
        return None
    p = os.path.realpath(os.path.join(ROOT, rel))
    root = os.path.realpath(ROOT)
    if not p.startswith(root + os.sep) or not os.path.isfile(p):
        return None
    return p


def _resource_read(uri):
    p = _resource_path(uri)
    if not p:
        return None
    try:
        text = open(p, encoding="utf-8", errors="replace").read()
    except OSError:
        return None
    mime = "text/markdown" if p.endswith(".md") else "application/yaml"
    return {"contents": [{"uri": uri, "mimeType": mime, "text": text[:200000]}]}

# ---- MCP prompts：两个"照做就不会踩坑"的现成提词 ----
PROMPTS = (
    {"name": "triage-failures",
     "description": "Look at failing steps and decide what to do, using the tool table",
     "arguments": []},
    {"name": "review-before-submit",
     "description": "Check generated inputs before submitting a step",
     "arguments": [{"name": "material", "description": "material name", "required": True},
                   {"name": "step", "description": "step label", "required": False}]},
)


def _prompts_list():
    return [{"name": p["name"], "description": p["description"],
             "arguments": p["arguments"]} for p in PROMPTS]


def _prompt_get(name, args=None):
    args = args or {}
    if name == "triage-failures":
        text = ("Call inspect with status=error to obtain the bounded failure snapshot, "
                "diagnosis codes, and deterministic proposals. For each failing step "
                "use the returned evidence to decide whether retry is appropriate; "
                "retry keeps products and must be explicitly selected before execution. "
                "Never try to perform stop, rerun, or clean through an automatic plan: "
                "report the evidence and let a human approve destructive actions.")
    elif name == "review-before-submit":
        mat = args.get("material") or "<material>"
        step = args.get("step") or "<step>"
        text = ("Call inspect with material %s and review the active step %s. If the "
                "contract or input evidence is insufficient, call describe_skill or "
                "get_snapshot for the same scope. Check the effective step settings, "
                "remote directory, and work_dir source before selecting start_step; "
                "execute only after the plan is explicit and its cursor is current."
                % (mat, step))
    else:
        return None
    return {"description": name, "messages": [{"role": "user",
                                               "content": {"type": "text", "text": text}}]}


# Shared protocol implementation.  The aliases preserve the historical private
# names used by tests and integrations while eliminating MCP/CLI policy drift.
_compact_contract = _protocol.compact_contract
_compact_state = _protocol.compact_state
_state_counts = _protocol.state_counts
_snapshot_flat = _protocol.snapshot_flat
_filter_snapshot = _protocol.filter_snapshot
_scope_data = _protocol.scope_data


if __name__ == "__main__":
    sys.exit(main())
