# -*- coding: utf-8 -*-
"""phonoagent mcp —— 把 PhonoAgent 以 MCP（Model Context Protocol）stdio 服务暴露给上层 agent。

设计三条硬约束（决定了"新技能仍然方便"和"agent 不能乱来"）：
  1. **只调 CLI，不 import 内部函数**：门禁/审计/档案只有一份事实来源，不可能与命令行漂移。
  2. **工具是通用动词，不认识任何具体技能**：技能清单/步骤清单都从 skill.yaml 动态读，
     所以加新技能 = 加目录，MCP 侧零改动、工具数量恒定（上限 20，超了就说明有人按技能加工具）。
  3. **风险在协议边界收口**：只读工具直连 CLI；变更/破坏性工具一律走 phonoagent act
     （风险分级 + 人工批准 + 审计日志），approve 只从真 TTY 生效，agent 无法自我批准。

协议：JSON-RPC 2.0 over stdio，支持 initialize / tools/list / tools/call（MCP 最小可用集）。
"""
import glob
import json
import os
import subprocess
import sys

PROTOCOL_VERSION = "2024-11-05"
SCHEMA_VERSION = "1"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "phonoagent")

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
    ("get_json", "结构化全量状态（types→materials→steps），带 schema_version",
     {"type": "object", "properties": {}, "schema_version": SCHEMA_VERSION},
     "read", ["json"]),
    ("conf_get", "查看某步 step.conf 合并后的最终值（只读）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"}},
      "required": ["material"]}, "read", ["conf"]),
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
    ("conf_set", "写项目配置（会改材料配置，需批准）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "step": {"type": "string"},
                                       "key": {"type": "string"},
                                       "value": {"type": "string"}},
      "required": ["material", "key", "value"]}, "mutate", ["conf"]),
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
    ("session_export", "导出可复现会话归档（manifest + 档案 + 重放提示）",
     {"type": "object", "properties": {"tt": {"type": "string"},
                                       "material": {"type": "string"},
                                       "out": {"type": "string"}},
      "required": ["material"]}, "mutate", ["session", "export"]),
]
BY_NAME = {t[0]: t for t in TOOLS}


def _cfg_args():
    cfg = (os.environ.get("PHONOAGENT_CONFIG")
           or os.environ.get("PHONOAGENT_CFG") or "").strip()
    return ["-c", cfg] if cfg else []


def _run(argv, timeout=1800):
    """只通过 CLI 执行——门禁/审计/档案全部沿用命令行那一套。"""
    cmd = [sys.executable, PROG] + _cfg_args() + argv
    env = dict(os.environ)
    env.setdefault("PHONOAGENT_ACTOR", "mcp")   # 让 act 网关认出这是 agent 会话
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           env=env)
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


def call_tool(name, args):
    """执行一个工具：只读直连，变更/破坏性走 phonoagent act（会返回需批准的错误）。"""
    if name not in BY_NAME:
        return {"isError": True, "content": [{"type": "text",
                "text": "unknown tool: %s" % name}]}
    _n, _d, _schema, risk, verb = BY_NAME[name]
    args = args or {}
    argv = _mat_args(args)
    if name == "conf_get":
        argv += ["conf"]
    elif name == "conf_set":
        argv += ["conf", "--set", "%s=%s" % (args.get("key"), args.get("value"))]
    elif name == "session_export":
        argv += ["session", "export"]
        if args.get("out"):
            argv += ["--out", str(args["out"])]
    elif name == "list_materials":
        if args.get("status"):
            argv += ["-status", str(args["status"])]
        argv += ["summary"]
    elif name == "get_summary":
        argv += ["summary"]
    else:
        argv += list(verb)
    if risk != "read":
        # 变更/破坏性：交给动作网关（agent 会话下破坏性动作会被拒并给出批准命令）
        argv = ["act"] + argv
    rc, out, err = _run(argv)
    text = (out or "") + (("\n" + err) if err.strip() else "")
    return {"isError": rc != 0, "content": [{"type": "text", "text": text.strip()[:20000]}],
            "structuredContent": {"risk": risk, "returncode": rc,
                                  "schema_version": SCHEMA_VERSION}}


def _tools_list():
    out = []
    for name, desc, schema, risk, _verb in TOOLS:
        s = dict(schema)
        s["description"] = desc
        s["x-risk"] = risk
        out.append({"name": name, "description": "[%s] %s" % (risk, desc),
                    "inputSchema": s})
    return out


def handle(req):
    method = req.get("method")
    rid = req.get("id")
    if method == "initialize":
        result = {"protocolVersion": PROTOCOL_VERSION,
                  "capabilities": {"tools": {}},
                  "serverInfo": {"name": "phonoagent", "version": SCHEMA_VERSION}}
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


if __name__ == "__main__":
    sys.exit(main())

# ---- MCP resources：把技能与项目文档作为只读资源暴露（agent 直接拿到上下文）----
def _skill_dirs():
    return sorted(d for d in glob.glob(os.path.join(ROOT, "skill", "*"))
                  if os.path.isdir(d))


def _resources_list():
    out = []
    for d in _skill_dirs():
        name = os.path.basename(d)
        for f in ("skill.yaml", "README.md"):
            p = os.path.join(d, f)
            if os.path.isfile(p):
                out.append({"uri": "phonoagent://skill/%s/%s" % (name, f),
                            "name": "%s/%s" % (name, f),
                            "mimeType": "text/markdown" if f.endswith(".md")
                            else "application/yaml"})
    for f in ("README.md", "README.en.md", "CHANGELOG.md", "docs/ACCEPTANCE.md",
              "docs/CONFIGURING.md", "docs/mcp.md"):
        if os.path.isfile(os.path.join(ROOT, f)):
            out.append({"uri": "phonoagent://doc/%s" % f, "name": f,
                        "mimeType": "text/markdown"})
    return out


def _resource_path(uri):
    """phonoagent://skill/<name>/<file> 或 phonoagent://doc/<path> → 仓库内真实路径。

    只允许仓库内的常规文件：任何 .. 或绝对路径都拒绝（agent 不能借资源接口越权读盘）。
    """
    if not isinstance(uri, str) or not uri.startswith("phonoagent://"):
        return None
    rest = uri[len("phonoagent://"):]
    if rest.startswith("skill/"):
        rel = os.path.join(rest.split("/", 1)[0], rest.split("/", 1)[1])
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
        text = ("Run get_summary, then get_status for every material it lists as FAIL. "
                "For each failing step read the diagnosis, decide between retry (keep "
                "products) and rerun (delete and regenerate), and remember that "
                "destructive tools are refused for agent sessions: propose the command "
                "and let a human approve it.")
    elif name == "review-before-submit":
        mat = args.get("material") or "<material>"
        step = args.get("step") or "<step>"
        text = ("Call conf_get for material %s (step %s) and get_status for the same "
                "material. Check the effective step.conf, the remote directory reported "
                "by status, and the work_dir source line; only then call start_step."
                % (mat, step))
    else:
        return None
    return {"description": name, "messages": [{"role": "user",
                                               "content": {"type": "text", "text": text}}]}
