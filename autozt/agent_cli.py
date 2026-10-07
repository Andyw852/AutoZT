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
import shutil
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from autozt import agent_protocol as _protocol
from autozt import agent_service as _service
from autozt import agent_state as _agent_state
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
REQUEST_OPS = set(_service.REQUEST_OPS)
REQUEST_FIELDS = {
    "op", "scope", "project", "tt", "material", "status", "step", "view", "limit",
    "include_monitoring", "include_retry", "execute", "dry_run", "max_actions",
    "cursor", "full", "skill", "plan", "actions", "goal", "result_dir", "poscar", "property",
    "temperature", "carrier", "direction", "dimension", "thickness", "source",
}
# 最近一次 _status 的数据来源（progress / live、多旧）；随结果一起返回给模型
_LAST_SOURCE: Dict[str, Any] = {}


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
    # start/retry 输出机器可读的 [agent-event] 行（提交/达上限/状态已变），执行器据此回填
    env["AUTOZT_AGENT_EVENTS"] = "1"
    return env


def autozt_argv() -> List[str]:
    """调用 AutoZT CLI 的 argv 前缀。源码树里有 bin/autozt 就用它；pip 非 editable
    安装时 bin/ 不在包里，退回 `python -m autozt`（以前直接报 cannot execute）。"""
    if os.path.isfile(PROG):
        return [sys.executable, PROG]
    return [sys.executable, "-m", "autozt"]


def _run(argv: Iterable[str], config: Optional[str], timeout: int = 1800
         ) -> Tuple[int, str, str]:
    cmd = autozt_argv()
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
    if getattr(args, "project", None):
        out += ["--project", str(args.project)]
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
                       for key in ("tt", "material", "status", "view")]
                      + (["project=%s" % args.project] if getattr(args, "project", None)
                         else []))
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


def _live_status(args: argparse.Namespace) -> Tuple[Optional[Dict[str, Any]], Optional[str], int]:
    data, error, rc = _call_json(_context_argv(args) + ["list", "--json"],
                                 getattr(args, "config", None))
    if error:
        return None, error, rc
    return _compact_state(data, args), None, 0


def _status(args: argparse.Namespace) -> Tuple[Optional[Dict[str, Any]], Optional[str], int]:
    """读紧凑状态：source=auto（默认）先看进度文件（秒级、无 ssh），太旧/没有/不含
    该范围才现采；live 强制现采；progress 强制读文件。来源记在 _LAST_SOURCE。"""
    scope = {key: getattr(args, key, None) for key in ("project", "tt", "material", "status")}
    compact, error, rc, info = _agent_state.load_compact(
        scope, getattr(args, "source", None), getattr(args, "config", None),
        lambda _scope: _live_status(args))
    _LAST_SOURCE.clear()
    _LAST_SOURCE.update(info)
    return compact, error, rc


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


def _cmd_progress(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    """读本地进度文件（monitor 每轮写）：不起子进程、不采集、不连超算。"""
    try:
        from autozt import load_config
        from autozt import progress as _progress
        path = _config_path(getattr(args, "config", None))
        cfg, found = load_config(path)
        cfg["_config_dir"] = os.path.dirname(os.path.abspath(found)) if found else os.getcwd()
    except SystemExit as exc:
        return _envelope("progress", ok=False, error=str(exc.code), returncode=1), 1
    doc = _progress.load_progress_view(cfg)
    if doc is None:
        return _envelope("progress", ok=False, returncode=1,
                         error="no progress file at %s yet; it is written by `autozt monitor` "
                               "each round or by any collecting command (e.g. `autozt list "
                               "--refresh`)" % _progress.progress_path(cfg)), 1
    view = _progress.filter_progress(doc, getattr(args, "project", None),
                                     getattr(args, "tt", None),
                                     getattr(args, "material", None),
                                     getattr(args, "status", None))
    view["path"] = _progress.progress_path(cfg)
    view["liveness"] = _progress.monitor_liveness(cfg, doc)
    limit = getattr(args, "limit", None)
    if limit:
        view["materials_total"] = len(view["materials"])
        view["materials"] = view["materials"][:int(limit)]
    return _envelope("progress", view), 0


AGENT_RULES_MD = """# AutoZT 调用规则（给 LLM；由 `autozt agent setup --rules` 生成）

1. 看状态先用 `autozt agent progress [--project P] [-status error]`（MCP: get_progress）：
   读 monitor 写的本地文件，秒级、不连超算。不要 grep 日志，不要解析 summary/list 的表格文本。
2. FAIL 分诊看 `fail_summary` 和 `materials[].fails[].action`：retry 桶可以提议 prepare_step；
   human_review 桶直接汇报给用户，不要自己反复重投。
3. 要候选动作用 `autozt agent inspect [--project P]`（MCP: inspect）。它默认读进度文件
   （返回里的 source 写明数据多旧）；确需实时状态才加 `--source live`（大体系要几十分钟）。
4. 执行只用 `autozt agent cycle --execute`（MCP: cycle execute=true）或 apply：只会 start/retry/
   fetch/advance；同技能同步骤合成一条命令，按现采状态逐个复核。结果看 `results[].outcome`：
   submitted=已提交，deferred=达 max_jobs 在排队（正常，不是故障），skipped_stale=状态已变
   （重新 inspect），failed=看 output。
5. 破坏性动作（stop/rerun/clean、-f/-y）不要自己决定：直接调命令（或 MCP
   request_destructive_action）会被拦下，返回**审批卡片**（命令、后果、任务图里标出目标步骤）
   和请求号。把卡片原样给用户看并询问；用户明确同意后执行
   `autozt approve --request <号> --reply "<用户原话>"`（MCP: approve_request）——就地批准并
   执行，不要让用户自己去终端敲命令；用户拒绝就 `autozt approve --deny <号>`。客户端支持 MCP
   确认框（elicitation）时会直接弹框，用户点了才执行。conf --set 改项目配置，同样先说明、得到同意。
6. 批量操作前可跑 `autozt agent doctor`（max_jobs 生效值、被屏蔽项目、同名材料）。
7. 材料用稳定 id `<项目名>/<完整名>`（progress/inspect 返回的 id 字段）指代，避免同名歧义。
8. 不改全局开关（auto_advance、auto_watch 等）和项目配置，除非用户明确要求。
9. 新材料/只有一份数据集：`autozt -tt <技能> -p <材料> register [--dataset D] [--poscar F]
   [--cluster local]`（MCP: register_material），改参数用 MCP conf_set（= conf --set）——两者都
   改项目配置，先向用户说明再调用。不要自己建目录、手跑 gen_*.py、手写 submit.sh 或 sbatch；
   没有超算就 `--cluster local`（本机执行，无 SLURM 时自动用仓库自带 tools/fakeslurm）。
10. 命令失败看退出码：非 0 即失败；代理环境下错误也会以 `[autozt error] …`（--json 时为
   `{"ok": false, "error": …}`）出现在 stdout。`-p` 写错材料名会直接报"找不到材料"。
11. `autozt` 不在 PATH 时用 setup 返回的 `cli.absolute`（绝对路径调用），或让用户
   `pip install -e <仓库>` / 把 `<仓库>/bin` 加进 PATH。
12. 汇报进度、或请用户同意任何操作之前，先给用户看任务图：`autozt -p <材料> graph`
   （整技能各步骤计数：`autozt -tt <技能> graph`；MCP: task_graph）。
"""


MCP_CONFIG_FILES = {"claude": (".mcp.json",), "cursor": (".cursor", "mcp.json")}


def _write_mcp_configs(server: Dict[str, Any], target: str) -> List[Dict[str, Any]]:
    """把 autozt 的 MCP 服务配置合并进软件目录下的项目级配置文件。

    只改 mcpServers.autozt 这一项，文件里其它服务/字段原样保留；文件坏了（不是 JSON）
    就不碰并报告。路径是本机绝对路径，所以这两个文件已加进 .gitignore。"""
    out = []
    for name in (["claude", "cursor"] if target == "all" else [target]):
        path = os.path.join(ROOT, *MCP_CONFIG_FILES[name])
        doc: Dict[str, Any] = {}
        if os.path.isfile(path):
            try:
                with open(path, encoding="utf-8") as fh:
                    doc = json.load(fh) or {}
            except (OSError, ValueError) as exc:
                out.append({"client": name, "path": path, "ok": False,
                            "error": "现有文件不是合法 JSON，未修改：%s" % exc})
                continue
        if not isinstance(doc, dict):
            out.append({"client": name, "path": path, "ok": False,
                        "error": "现有文件顶层不是对象，未修改"})
            continue
        servers = doc.setdefault("mcpServers", {})
        changed = servers.get("autozt") != server
        servers["autozt"] = server
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, ensure_ascii=False, indent=2)
            fh.write("\n")
        os.replace(tmp, path)
        out.append({"client": name, "path": path, "ok": True, "changed": changed})
    return out


def _cmd_setup(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    """接入说明：MCP 服务配置（绝对路径，直接粘进客户端配置）+ 给 LLM 的规则卡。只读。"""
    config = _config_path(getattr(args, "config", None))
    argv = autozt_argv()
    env = {"AUTOZT_MCP_PROFILE": "core"}
    if config:
        env["AUTOZT_CONFIG"] = config
    data: Dict[str, Any] = {
        "mcp_server": {"mcpServers": {"autozt": {
            "command": argv[0], "args": argv[1:] + ["mcp"], "env": env}}},
        "mcp_profiles": {"core": "默认：12 个工具，按 发现/规划/执行/审批/汇报 五段组织（tools/list 的 x-stage）",
                         "workflow": "21 个：core + inspect/apply_actions/conf_get/conf_set/doctor 等",
                         "monitor": "最小：get_progress/get_snapshot/cycle",
                         "full": "全部工具（含需人工批准的破坏性工具）"},
        "cli": {"prefix": " ".join(argv + ["agent"]),
                "absolute": " ".join(argv),
                "on_path": shutil.which("autozt"),
                "install_hint": (None if shutil.which("autozt") else
                                 "autozt 不在 PATH：pip install -e %s，或 export PATH=%s:$PATH；"
                                 "也可直接用 cli.absolute" % (ROOT, os.path.join(ROOT, "bin"))),
                "first_calls": ["progress --project <P>", "doctor",
                                "inspect --project <P>", "cycle --project <P> --execute"]},
        "rules_md": AGENT_RULES_MD,
        "config": config,
    }
    target = getattr(args, "write", None)
    if target:
        data["written"] = _write_mcp_configs(data["mcp_server"]["mcpServers"]["autozt"],
                                             target)
    try:
        from autozt import progress as _progress
        cfg, _err = _agent_state.load_cfg(getattr(args, "config", None))
        if cfg is not None:
            doc = _progress.load_progress(cfg)
            data["progress"] = {"path": _progress.progress_path(cfg), "exists": doc is not None}
            if doc is not None:
                data["progress"].update(_progress.monitor_liveness(cfg, doc))
                data["progress"]["coverage"] = doc.get("coverage") or {}
    except Exception:                             # noqa: BLE001
        pass
    _bad_write = any(not w.get("ok") for w in data.get("written") or [])
    return _envelope("setup", data), (1 if _bad_write else 0)


def _cmd_doctor(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    argv = ["doctor", "--json"]
    if getattr(args, "project", None):
        argv = ["--project", str(args.project)] + argv
    if getattr(args, "tt", None):
        argv = ["-tt", str(args.tt)] + argv
    rc, out, err = _run(argv, getattr(args, "config", None))
    payload, parse_error = _parse_json_output(out, err)
    if parse_error:
        return _envelope("doctor", ok=False, error=_cli_error(rc, out, err) if rc else parse_error,
                         returncode=rc or 1), rc or 1
    # doctor 有 error 级发现时退出码 1，但报告本身是完整有效的——照样返回
    return _envelope("doctor", payload, ok=rc == 0, returncode=rc,
                     error=None if rc == 0 else "doctor reported error-level findings"), rc


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
    """--result-dir：结果回收后的校验；否则 --poscar / -p 材料：执行前的输入检查。"""
    root = getattr(args, "result_dir", None)
    common = {"dimension": getattr(args, "dimension", None),
              "thickness": getattr(args, "thickness", None),
              "temperature": getattr(args, "temperature", None),
              "carrier": getattr(args, "carrier", None)}
    if root:
        return _envelope("preflight", _science.preflight(root, **common)), 0
    poscar = getattr(args, "poscar", None)
    material = getattr(args, "material", None)
    if not poscar and material:
        argv = (["-tt", args.tt] if getattr(args, "tt", None) else [])
        argv += ["-p", material, "list", "--json"]
        listing, error, rc = _call_json(argv, getattr(args, "config", None))
        if error:
            return _envelope("preflight", ok=False, error=error, returncode=rc or 1), rc or 1
        poscar = _science.poscar_from_listing(listing, material)
        if not poscar:
            return _envelope("preflight", ok=False, returncode=1,
                             error="preflight: 材料 %s 的本地目录下没有 POSCAR" % material), 1
    if not poscar:
        return _envelope("preflight", ok=False, returncode=2,
                         error="preflight 需要 --result-dir（结果校验）或 --poscar / -p 材料"
                               "（执行前输入检查）"), 2
    return _envelope("preflight", _science.preflight_inputs(poscar, material=material, **common)), 0


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
    kind = _LAST_SOURCE.get("kind") or "live"
    record = {"version": SNAPSHOT_VERSION, "cursor": cursor, "snapshot": compact,
              "updated": int(time.time()), "source": kind}
    _save_snapshot(path, record)
    if previous is not None and previous.get("source", "live") != kind:
        # 进度文件与现采的紧凑状态字段不完全相同，跨来源比 diff 会全是假变化
        previous = None
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
    data["source"] = dict(_LAST_SOURCE)
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
            "proposal_limit": args.max_actions,
            "source": dict(_LAST_SOURCE)}
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
        config=getattr(args, "config", None), project=scope.get("project"),
        tt=scope.get("tt"), material=scope.get("material"),
        status=scope.get("status"), source=_plan_source(plan, args))
    compact, error, rc = _status(current_args)
    if error:
        return None, error, rc
    actual = _protocol.state_cursor(compact or {})
    if str(expected) != actual:
        return {"stale_plan": True, "expected_cursor": str(expected),
                "current_cursor": actual, "actions": actions}, \
            "plan is stale; collect a new snapshot before applying actions", 1
    return None, None, 0


def _plan_source(plan: Any, args: Any = None) -> str:
    """计划是从哪种来源做出来的（复核 cursor 必须用同一来源，否则必然对不上）。"""
    src = _plan_field(plan, "source")
    if isinstance(src, dict):
        src = src.get("kind")
    if src in ("progress", "live"):
        return str(src)
    return getattr(args, "source", None) or "live"


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
    return {key: getattr(args, key, None) for key in ("project", "tt", "material", "status")
            if getattr(args, key, None) not in (None, "")}


def _service_loader(args: argparse.Namespace):
    pinned: Dict[str, str] = {}

    def load(scope: Mapping[str, Any]):
        # 同一次 inspect/cycle 里的复核读必须和第一次同一来源（auto 第一次选定后钉住）
        scoped = argparse.Namespace(
            config=getattr(args, "config", None), project=scope.get("project"),
            tt=scope.get("tt"), material=scope.get("material"),
            status=scope.get("status"),
            source=pinned.get("kind") or getattr(args, "source", None))
        got = _status(scoped)
        if not got[1] and _LAST_SOURCE.get("kind"):
            pinned.setdefault("kind", _LAST_SOURCE["kind"])
        return got
    load.pinned = pinned  # type: ignore[attr-defined]
    return load


def _cmd_inspect(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    loader = _service_loader(args)
    data, error, rc = _service.inspect(
        loader, _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS),
        include_retry=True)
    if error:
        return _envelope("inspect", ok=False, error=error, returncode=rc), rc
    data["source"] = dict(_LAST_SOURCE)
    return _envelope("inspect", data), 0


def _cmd_plan(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    loader = _service_loader(args)
    data, error, rc = _service.plan(
        loader, _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS))
    if error:
        return _envelope("plan", ok=False, error=error, returncode=rc), rc
    data["source"] = dict(_LAST_SOURCE)
    return _envelope("plan", data), 0


def _execute_actions(actions: List[Dict[str, Any]], config: Optional[str],
                     dry_run: bool = False,
                     expected_cursor: Optional[str] = None,
                     scope: Optional[Mapping[str, Any]] = None,
                     source: Optional[str] = None) -> Dict[str, Any]:
    """Execute only the bounded non-destructive action vocabulary through ``act``.

    同一 (动作, 技能, 步骤) 合成一条 ``act -p A,B,C``：一次采集、一个共享 max_jobs
    闸门；start/retry 带 ``--expect-state``，CLI 按现采状态逐个复核（见 agent_state）。"""
    if dry_run:
        return {"results": [], "actions": actions, "dry_run": True, "failed": False}
    for action in actions:
        if action.get("action") not in _protocol.ALLOWED_ACTIONS:
            return {"results": [], "actions": actions, "failed": True,
                    "error": "unsupported or destructive action: %s" % action.get("action")}
    if expected_cursor:
        # The service has already rechecked once; this final read is the CAS
        # boundary immediately before invoking the mutation gateway.
        scoped = argparse.Namespace(config=config, tt=(scope or {}).get("tt"),
                                    project=(scope or {}).get("project"),
                                    material=(scope or {}).get("material"),
                                    status=(scope or {}).get("status"),
                                    source=source or "live")
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
    return _agent_state.execute_grouped(actions, lambda argv: _run(argv, config))


def _cmd_cycle(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    execute = bool(getattr(args, "execute", False)) and not bool(
        getattr(args, "dry_run", False))
    loader = _service_loader(args)
    data, error, rc = _service.cycle(
        loader,
        lambda actions, dry_run, cursor=None: _execute_actions(
            actions, getattr(args, "config", None), dry_run, cursor,
            _agent_scope(args), loader.pinned.get("kind")),
        _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS),
        execute=execute,
        include_retry=bool(getattr(args, "include_retry", False)))
    if isinstance(data, dict):
        data["source"] = dict(_LAST_SOURCE)
    if error:
        return _envelope("cycle", data, ok=False, error=error, returncode=rc), rc
    return _envelope("cycle", data), 0


def _cmd_run(args: argparse.Namespace) -> Tuple[Dict[str, Any], int]:
    """Run the deterministic AutoZT loop; dry-run is the default."""
    execute = bool(getattr(args, "execute", False)) and not bool(
        getattr(args, "dry_run", False))
    loader = _service_loader(args)
    data, error, rc = _service.run(
        loader,
        lambda actions, dry_run, cursor=None: _execute_actions(
            actions, getattr(args, "config", None), dry_run, cursor,
            _agent_scope(args), loader.pinned.get("kind")),
        _agent_scope(args),
        view=getattr(args, "view", "attention"),
        include_monitoring=bool(getattr(args, "include_monitoring", False)),
        max_actions=getattr(args, "max_actions", _protocol.MAX_ACTIONS),
        execute=execute,
        include_retry=bool(getattr(args, "include_retry", False)))
    if isinstance(data, dict):
        data["source"] = dict(_LAST_SOURCE)
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
    for key in ("project", "tt", "material", "status", "step"):
        if request.get(key) not in (None, ""):
            merged[key] = request[key]
    return argparse.Namespace(
        config=getattr(base, "config", None), project=merged.get("project"),
        limit=request.get("limit"),
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
        poscar=request.get("poscar"),
        property=request.get("property"), temperature=request.get("temperature"),
        carrier=request.get("carrier"), direction=request.get("direction"),
        dimension=request.get("dimension"), thickness=request.get("thickness"),
        source=request.get("source"),
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
        bad_scope = sorted(set(scope) - {"project", "tt", "material", "status", "step"})
        if bad_scope:
            return _envelope("request", ok=False,
                             error="unknown scope field(s): %s" % ", ".join(bad_scope),
                             returncode=2), 2
        for key, value in scope.items():
            if not isinstance(value, str):
                return _envelope("request", ok=False,
                                 error="scope field %s must be string" % key,
                                 returncode=2), 2
    for key in ("op", "project", "tt", "material", "status", "step", "view", "cursor",
                "skill", "source"):
        if key in request and not isinstance(request[key], str):
            return _envelope("request", ok=False,
                             error="request field %s must be string" % key,
                             returncode=2), 2
    if "source" in request and request["source"] not in _agent_state.SOURCES:
        return _envelope("request", ok=False,
                         error="request field source must be auto, progress, or live",
                         returncode=2), 2
    if "view" in request and request["view"] not in ("attention", "active", "all"):
        return _envelope("request", ok=False,
                         error="request field view must be attention, active, or all",
                         returncode=2), 2
    if "actions" in request and not isinstance(request["actions"], list):
        return _envelope("request", ok=False,
                         error="request field actions must be array", returncode=2), 2
    if "limit" in request and (isinstance(request["limit"], bool)
                               or not isinstance(request["limit"], int)
                               or request["limit"] < 1):
        return _envelope("request", ok=False,
                         error="request field limit must be a positive integer",
                         returncode=2), 2
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
    if op == "progress":
        return _cmd_progress(child)
    if op == "doctor":
        return _cmd_doctor(child)
    if op == "setup":
        return _cmd_setup(child)
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
                        help="材料名或稳定 id（<项目名>/<完整名>）")
    parser.add_argument("-proj", "--project", dest="project", default=argparse.SUPPRESS,
                        help="只看指定项目（tf_<项目>.yaml 的 <项目>，逗号分隔）")
    parser.add_argument("-status", "--status", dest="status", default=argparse.SUPPRESS,
                        help="状态过滤，例如 error,running,pd")
    parser.add_argument("--pretty", action="store_true", default=argparse.SUPPRESS,
                        help="缩进输出 JSON")


def _add_source_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--source", choices=_agent_state.SOURCES, default=None,
                        help="状态来源：auto（默认，或 AUTOZT_AGENT_SOURCE）=进度文件够新就读它"
                             "（秒级、无 ssh），否则现采；live=强制现采（慢）；progress=只读进度文件")


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

    p = sub.add_parser("progress", help="读本地进度文件：逐材料状态/FAIL 码/ETA（无 ssh，秒级）")
    _add_context_options(p)
    p.add_argument("--limit", type=int, help="最多返回多少个材料")
    p = sub.add_parser("doctor", help="配置预检：max_jobs 生效值/屏蔽项目/同名材料（无 ssh）")
    _add_context_options(p)
    p = sub.add_parser("setup", help="接入说明：MCP 服务配置（绝对路径）+ 给 LLM 的规则卡")
    _add_context_options(p)
    p.add_argument("--rules", action="store_true",
                   help="只输出规则卡（Markdown 原文），可追加进 agent 的 AGENTS.md")
    p.add_argument("--write", nargs="?", const="all", choices=["all", "claude", "cursor"],
                   help="把 MCP 服务配置合并写进软件目录：claude=.mcp.json（Claude Code 等"
                        "读项目级 .mcp.json 的客户端）、cursor=.cursor/mcp.json、all=两个都写。"
                        "已有的其它 MCP 服务不动；在软件目录启动客户端即自动接上 autozt mcp")

    p = sub.add_parser("snapshot", help="持久化增量状态快照")
    _add_context_options(p)
    _add_source_option(p)
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

    p = sub.add_parser("preflight", help="执行前检查输入（--poscar 或 -p 材料），"
                                         "或结果回收后校验结果（--result-dir）；只读")
    _add_context_options(p)
    p.add_argument("--result-dir", help="给出 = 结果校验；不给 = 执行前输入检查")
    p.add_argument("--poscar", help="执行前输入检查用的 POSCAR（或用 -p 材料）")
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
        _add_source_option(p)
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
    _add_source_option(p)
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
    # 只删裸命令词 agent（-p agent 这种取值不能删）。
    from autozt.cli import _bare_index
    _i = _bare_index(raw, "agent")
    if _i is not None:
        raw.pop(_i)
    args = parser.parse_args(raw)
    if not args.agent_command:
        parser.print_help()
        return 2
    if args.agent_command == "skills":
        result, rc = _cmd_skills(args)
    elif args.agent_command == "progress":
        result, rc = _cmd_progress(args)
    elif args.agent_command == "doctor":
        result, rc = _cmd_doctor(args)
    elif args.agent_command == "setup":
        if getattr(args, "rules", False):
            sys.stdout.write(AGENT_RULES_MD)
            return 0
        result, rc = _cmd_setup(args)
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
