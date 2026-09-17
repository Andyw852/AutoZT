# MCP interface

AutoZT can expose itself as a Model Context Protocol (MCP) server so that an
LLM-driven planner (or any MCP-capable client) can drive the workflow without bespoke
glue code, while every dangerous action still passes through the approval gate.

The stdio server negotiates MCP `2025-06-18` and keeps `2024-11-05` as a fallback for
older clients.

    autozt mcp                 # JSON-RPC 2.0 over stdio: initialize / tools/list / tools/call
    autozt mcp --list-tools    # print the tool table (for humans and tests)
    autozt mcp --call NAME JSON

## Design rules

1. The server only calls the command-line interface. Approval, audit and provenance have
   a single source of truth and cannot drift from the CLI.
2. Tools are generic verbs, never per-skill commands. `skill.yaml` is the skill
   contract: its `io_schema`, `flow`, `corrections`, and step definitions are exposed
   by `describe_skill`; the README is available as an MCP resource. Adding or changing
   a skill therefore does not change the MCP tool surface. The table is capped at 21
   entries on purpose.
3. Risk is enforced at the protocol boundary. Read-only tools use AutoZT's explicit
   `list`/`summary`/diagnostic collectors; mutating and
   destructive tools are routed through the action gateway, so an agent session receives
   "needs approval" instead of executing, and approvals only work from a real TTY.
4. The MCP layer may aggregate and reduce CLI output, but it never imports skill
   generators or changes workflow semantics. Skill fixes and new skill development stay
   in their existing directories.

## Tools

| Tool | Risk | Purpose |
|---|---|---|
| list_skills | read | skills and their steps, from skill.yaml |
| list_materials | read | materials and step states |
| get_summary | read | per-skill counters, FAIL list, global queue |
| get_status | read | one material in detail (diagnosis, cluster, work_dir source) |
| capabilities | read | fixed agent protocol, actions, and risk boundary |
| conf_get | read | effective step.conf of one step |
| describe_skill | read | machine-readable skill contract (I/O, DAG, checks, corrections) |
| probe_step | read | adaptive probe for defect-dft-cpu, generic diagnosis for other skills |
| start_step | mutate | generate inputs if needed, then submit |
| retry_step | mutate | regenerate inputs, keep products, do not submit |
| fetch_results | mutate | pull finished results back |
| advance_ready | mutate | advance all dependency-ready steps for one skill or all skills |
| inspect | read | one call: scoped attention state + deterministic actions |
| stop_step | destructive | cancel the job of a step |
| rerun_step | destructive | delete the step directory and regenerate |
| clean_material | destructive | delete generated files, back to PREP |
| cycle | mutate | one dry-run or explicitly executed safe cycle |
| get_snapshot | read | compact state; cursor returns only changes |
| propose_actions | read | deterministic candidate actions from status/diagnosis |
| apply_actions | mutate | execute a bounded batch of non-destructive actions |

`tools/list` also returns the MCP tool annotations `readOnlyHint`,
`destructiveHint`, and `idempotentHint`, plus AutoZT's explicit `x-risk` field.
Clients can present the risk before a call without parsing natural-language descriptions.
Every tool declares an `outputSchema`; successful and failed calls return the same compact
JSON envelope in `structuredContent`. Full mode also serializes that envelope in text for
older clients; the compact, workflow, and monitor profiles keep only a short text pointer
to avoid duplicating large payloads.

## Workflow profile for new LLM integrations

新接入优先使用 `AUTOZT_MCP_PROFILE=workflow`。它把一次观察、规则规划和可选执行收敛为少量稳定工具：

```text
inspect → (model reviews only if needed) → cycle(execute=false) → cycle(execute=true)
```

`inspect` 一次返回 attention 状态、稳定诊断码、候选动作和可执行动作；模型不需要先分别调用
状态、提案和动作工具。`cycle` 默认 dry-run，只有显式传 `execute=true` 才会把
`start_step`、`fetch_results`、`advance_ready` 交给 `autozt act`；`retry_step` 还需要显式
传 `include_retry=true`。`stop`、
`rerun`、`clean`、配置修改和其它高风险操作仍留在普通 CLI/人工审批路径。

```json
{
  "mcpServers": {
    "autozt": {
      "command": "autozt",
      "args": ["mcp"],
      "env": {"AUTOZT_MCP_PROFILE": "workflow"}
    }
  }
}
```

stdio 客户端启动进程后按 MCP 顺序发送 `initialize`、`tools/list`，之后直接调用：

```json
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"inspect","arguments":{"tt":"band-dft-cpu","status":"error"}}}
{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"cycle","arguments":{"tt":"band-dft-cpu","execute":false}}}
```

返回值从 `result.structuredContent` 读取；`result.content[0].text` 只是兼容旧客户端的短指针。
模型确认计划后才把 `execute` 改为 `true`，并把 `include_retry` 明确打开。
执行型 `cycle` 在提交动作前会再次读取同一作用域；cursor 变化就返回
`stale_plan`，整轮不执行。

如果 MCP 客户端本身运行在 WSL 内，建议使用 AutoZT 的绝对路径，避免客户端工作目录
改变后找不到程序：

```json
{
  "mcpServers": {
    "autozt": {
      "command": "/home/wangchao/software/AutoZT/bin/autozt",
      "args": ["mcp"],
      "env": {
        "AUTOZT_MCP_PROFILE": "workflow",
        "AUTOZT_CONFIG": "/home/wangchao/software/AutoZT/setting/tf.yaml"
      }
    }
  }
}
```

如果客户端是 Windows 原生进程，则让它启动 WSL 的 stdio 子进程；环境变量要放进
`bash -lc` 内部：

```json
{
  "mcpServers": {
    "autozt": {
      "command": "wsl.exe",
      "args": [
        "-d", "Ubuntu", "--", "bash", "-lc",
        "export AUTOZT_MCP_PROFILE=workflow AUTOZT_CONFIG=/home/wangchao/software/AutoZT/setting/tf.yaml; exec /home/wangchao/software/AutoZT/bin/autozt mcp"
      ]
    }
  }
}
```

如果 MCP 客户端只负责长期巡检，使用更小的 `monitor` profile。它只暴露
`get_snapshot` 和 `cycle` 两个工具：前者配合 cursor 做增量轮询，后者默认 dry-run，
只有显式传 `execute=true` 才推进 AutoZT 已允许的安全动作。技能发现、失败探测和人工复核
在出现异常时再切换到 `workflow` profile，避免把完整工具表和技能目录长期放进模型上下文。

```json
{
  "mcpServers": {
    "autozt": {
      "command": "autozt",
      "args": ["mcp"],
      "env": {"AUTOZT_MCP_PROFILE": "monitor"}
    }
  }
}
```

## 本地 MCP 与 API 接入的边界

`autozt mcp` 是本地 **stdio** 服务，适合 Claude、Kimi、ChatGPT/Codex 等能启动本地 MCP
进程的客户端。它的配置是客户端启动 `autozt mcp`，不是把一个本地命令直接填进
Responses API 的 `server_url`。OpenAI Responses API 的 MCP 工具要求可访问的远程
Streamable HTTP 或 HTTP/SSE 服务；如果要从 API 调用 AutoZT，需要在受控环境中做一个
远程适配层，把 `autozt agent serve` 或同一套协议转成远程 MCP，并继续保留 `autozt act`
作为唯一变更入口。

API 适配层只应先暴露高层工具，例如 `inspect`、`get_snapshot` 和 `cycle`，并用
`allowed_tools` 限制导入集合。对支持 OpenAI Tool Search 的模型，可以把 AutoZT 作为一个
有清晰描述的 MCP/namespace 延迟加载面（`defer_loading: true`），让模型需要时再载入具体
工具；这不会改变 AutoZT 的技能契约，也不要求每个技能增加一个工具。若客户端不支持
Tool Search，直接使用 `AUTOZT_MCP_PROFILE=workflow` 即可达到同样的固定小工具面。

官方 OpenAI 文档说明：MCP 工具定义和工具输出都会进入模型上下文；保留会话中的工具列表
可以避免每轮重新拉取，`allowed_tools` 可以减少初始工具面，Tool Search 可以把大工具集
延迟到需要时加载。AutoZT 的本地实现对应关系是：`workflow`/`monitor` 减少工具定义，
`get_snapshot(cursor)` 只返回变化，`summary --diff` 在无变化时完全不调用模型；没有必要
为了“使用 MCP”而让模型参与每十分钟的普通巡检。

参考：

- [OpenAI：MCP and Connectors](https://developers.openai.com/api/docs/guides/tools-connectors-mcp)
- [OpenAI：Tool Search](https://developers.openai.com/api/docs/guides/tools-tool-search)

## Compact LLM profile

The default profile is `workflow`, which exposes the small model-facing control loop. The full
profile exposes the complete fixed generic tool table and must be selected explicitly with
`AUTOZT_MCP_PROFILE=full`. For an even smaller compatibility surface, set
`AUTOZT_MCP_PROFILE=compact` before starting the server. `tools/list` then
exposes only seven tools (including the one-shot `schema` discovery tool):

For a stdio MCP client, the configuration is equivalent to:

```json
{
  "mcpServers": {
    "autozt": {
      "command": "autozt",
      "args": ["mcp"],
      "env": {"AUTOZT_MCP_PROFILE": "compact"}
    }
  }
}
```

```text
schema (first connection only)
list_skills → describe_skill → get_snapshot → probe_step
                                ↘ propose_actions → apply_actions
```

The compact profile is an interface-size optimization, not a permission bypass. Read-only
tools remain read-only, and `apply_actions` accepts only `start_step`, `retry_step`,
`fetch_results`, and `advance_ready`; `stop`, `rerun`, `clean`, and force flags are
rejected before they reach the CLI.

In compact mode, large results are carried once in MCP `structuredContent`; the text
content contains only a small envelope pointing to that structured payload. Full mode
keeps the serialized result in both fields for older clients.

`get_snapshot` is designed for polling. The first call returns a short cursor and the
requested view (`attention` by default). Passing that cursor on the next call returns
only changed steps. The cursor is an in-memory cache key, so a new server process resets
it safely and returns `cursor_reset: true`. `tt`, `material`, and `status` scopes are
applied again inside the MCP layer, because the CLI's machine-readable JSON branches
may retain the full collection for compatibility.

`propose_actions` is deliberately advisory. It translates the existing stable
`diag_code` and `suggested_action` fields into candidate tool calls; it does not submit
anything. The model can inspect the evidence, then pass selected candidates to
`apply_actions`. The batch stops at the first failed action and returns per-action
structured results.

`inspect`, `propose_actions`, and `get_snapshot` return the same scoped state `cursor`.
When the model passes that cursor to `apply_actions`, the server re-collects the scoped
state before executing. A changed cursor returns `stale_plan` and executes nothing, so a
queued job completing between planning and execution cannot make the model apply an old
decision silently. Omitting the cursor remains valid for an explicit low-level call, but
model-generated plans should always carry it.

Typical split calls are therefore:

```json
{"name":"propose_actions","arguments":{"tt":"band-dft-cpu","status":"error"}}
{"name":"apply_actions","arguments":{"tt":"band-dft-cpu","cursor":"<cursor-from-propose_actions>","actions":[{"action":"retry_step","material":"C24/qHPC24","step":"S1_opt"}]}}
```

## Why This Boundary

AutoZT should remain the durable workflow engine: the state machine, dependency graph,
physical checks, retries, queue monitoring, and audit trail stay deterministic. This is
the same separation used by durable workflow systems such as
[Temporal Workflows](https://docs.temporal.io/workflows), which replay an event history
to recover a workflow, and by [Prefect flows](https://docs.prefect.io/v3/concepts/flows),
which keep retries and timeouts in the flow runtime. An LLM is useful above that layer for
failure triage, evidence interpretation, and selecting among already-declared actions.

The checkpoint and approval shape also follows the pattern documented by
[LangGraph persistence](https://docs.langchain.com/oss/python/langgraph/persistence):
persist a small cursor/checkpoint, interrupt for human review when needed, and resume the
same workflow instead of rebuilding the whole conversation. In AutoZT, `summary --diff`
and the cursor snapshot are the checkpoint, while `autozt act` and the TTY approval gate
are the authoritative mutation boundary.

## Skill contract for models

Every skill should keep its existing executable implementation and add or maintain only
the declarative fields in `skill.yaml`:

- `io_schema`: required inputs, declared outputs, and user-adjustable parameters;
- `flow`: stages, dependencies, optional branches, and produced artifacts;
- `corrections`: recognized failure families and supported recovery policies;
- `steps` / `optional_steps`: labels, checks, generators, and completion markers.

`resources/list` and `resources/read` include the bundled skill directories and any
directories named by the active config's `skill_paths`. A drop-in skill can therefore
publish its own `skill.yaml` and `README.md` to the model without adding an MCP tool or
changing `autozt/mcp.py`; the configured search path keeps the same precedence as CLI
skill discovery.

When a new skill has only executable `steps`, `describe_skill` reports the missing
contract fields and supplies an explicitly marked structural `inferred.step_graph`.
That fallback helps an LLM navigate an early skill without pretending to know its
physical inputs or recovery policy; adding the declarative fields later upgrades the
same skill in place.

The MCP server does not hard-code a list of skills or steps. A model first calls
`list_skills`, then `describe_skill(tt)` for the selected capability, and only then
uses the material-scoped state and action tools. This keeps the interface stable as
skills are revised or new ones are added.

## CLI counterpart

MCP is an optional transport, not the workflow engine. The same model-facing loop is
available without a persistent MCP client:

```bash
autozt agent capabilities
autozt agent schema              # 一次性读取完整请求/响应/风险 Schema
autozt agent skills
autozt agent contract band-dft-cpu
autozt agent inspect
autozt agent plan --status error
autozt agent cycle --dry-run
autozt agent run --dry-run       # 无需 MCP 的确定性闭环
autozt agent run --execute       # 显式执行允许的安全动作
autozt agent snapshot
autozt agent propose > plan.json
autozt agent apply plan.json --dry-run
autozt agent apply plan.json
autozt agent serve                 # 常驻 JSONL：一行请求对应一行 JSON 响应
```

`autozt agent` emits one stable JSON envelope per invocation and stores its snapshot on
disk, so the cursor survives a process restart. It reuses the same skill schema, state
normalisation and `autozt act` gateway as MCP. Plans only accept `start_step`,
`retry_step`, `fetch_results` and `advance_ready`; destructive operations remain outside
the automatic plan interface and require the normal human approval path.

## Token behaviour

MCP itself does not reduce tokens. The workflow profile reduces tool calls; the compact
profile reduces tool-definition tokens; the `monitor` profile reduces the steady-state
tool table to two tools; and `get_snapshot(cursor=...)` reduces repeated state tokens. The
server still performs the full AutoZT collection needed for correctness; it only sends a
smaller result to the model. Clients should avoid full-state reads during polling and use
`get_snapshot(cursor=...)` instead; reserve `view=all` for explicit batch analysis. The
10-minute steady-state monitor should run `autozt summary --diff` without an LLM call;
involve the model only when a new FAIL, queue anomaly, or ambiguous diagnosis appears.

With the current schemas, the compact JSON tool definitions measure about 13,455
characters in the full profile, 4,808 in compact, 6,243 in workflow, and 1,552 in
monitor (compact and workflow are about 64% and 54% smaller; monitor about 88% smaller
before the client adds its own wrapper). These are transport/context savings,
not a reduction in VASP or SSH work. Re-measure after changing a tool schema because the
values are not a protocol guarantee.

## Measured behaviour

safety_metrics.py reports the numbers that can be cited:

    interception rate of hazardous actions  100% (9/9)
    hazardous actions actually executed        0
    every denial offers an approval path     yes
    read-only actions still allowed           yes
    local output determinism                 100%

Run it with

    python scripts/safety_metrics.py          # human readable
    python scripts/safety_metrics.py --json   # machine readable (CI asserts on this)

Protocol references: [MCP Tools](https://modelcontextprotocol.io/specification/2025-06-18/server/tools),
[Resources](https://modelcontextprotocol.io/specification/2025-06-18/server/resources), and
[Prompts](https://modelcontextprotocol.io/specification/2025-06-18/server/prompts). The
implementation follows the optional `outputSchema`/`structuredContent` contract and
treats server annotations as hints; the action gateway remains the authoritative safety
boundary.
