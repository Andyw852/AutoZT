# `autozt agent` CLI

`autozt agent` 是 AutoZT 的稳定 JSON 接口。它和 MCP 使用同一套技能契约、状态字段和动作网关，但不需要 MCP 客户端；人、cron、Shell、Python 脚本或只会调用命令的模型都可以使用它。

技能开发规范统一维护在软件包根目录 `TASKFLOW.md` 第 7 章，尤其 7.15–7.18；入口与服务生命周期见 9.1。`io_schema` 是描述契约，不是任意科学输入输出的自动解析器。新增符合现有生命周期的技能通常无需修改本接口；新增执行语义应先扩展核心，再同步共享 Agent 协议。普通集成优先使用一次性 `request -`，`serve` 仅供需要连续 JSONL 管道的高级 wrapper。

```bash
# 第一次接入：只读协议能力，不采集超算状态
autozt agent capabilities
# 完整机器可读请求/响应/风险 Schema（建议首次只调用一次）
autozt agent schema

# 技能目录：只返回模型需要的字段
autozt agent skills

# 一个技能的输入、输出、步骤和纠错契约
autozt agent contract band-dft-cpu
# 需要 generator 等完整字段时再显式放大输出
autozt agent contract band-dft-cpu --full

# 第一次输出当前关注状态；以后只输出变化，快照保存在配置目录
autozt agent snapshot
autozt agent snapshot --view active

# 根据状态和稳定诊断码生成候选动作，不提交、不改配置
autozt agent propose > plan.json

# 先验证计划，不执行；确认后才执行
autozt agent apply plan.json --dry-run
autozt agent apply plan.json

# 新接入的单次闭环：一次观察 + 规则计划；默认不提交
autozt agent inspect
autozt agent plan --status error
autozt agent cycle --dry-run
autozt agent cycle --execute       # 只执行 start/retry/fetch/advance
# 无需 MCP 的确定性闭环；默认 dry-run，--execute 才执行
autozt agent run --dry-run
autozt agent run --execute

# 失败证据：只读，输出 diagnose 的结构化结果
autozt agent evidence -tt band-dft-cpu -p C24/qHPC24 --step S1_opt

# JSON request：避免模型拼接 shell 参数
echo '{"op":"inspect","scope":{"tt":"band-dft-cpu","status":"error"}}' |
  autozt agent request -

# 常驻 JSONL：每个非空输入行对应一个 JSON 输出行
printf '%s\n%s\n' '{"op":"capabilities"}' '{"op":"snapshot"}' |
  autozt agent serve
```

全局上下文可以放在 `agent` 前面，也可以放在子命令后面：

```bash
autozt -c ~/software/AutoZT/setting/tf.yaml -tt band-dft-cpu agent snapshot
autozt agent propose -tt band-dft-cpu -p C24/qHPC24
autozt agent contract -tt band-dft-cpu
```

默认输出一行 JSON；加 `--pretty` 只改变排版，不改变字段。`snapshot` 的本地状态文件名形如 `.tf_agent_snapshot_<hash>.json`，位于 AutoZT 配置目录，不写材料目录，也不保存远端文件内容。

## 计划格式

`apply` 接受 `propose` 的完整输出、只含 `proposals` 的对象，或只含 `actions` 的对象。动作是通用动词，与技能无关：

```json
{
  "actions": [
    {"action": "start_step", "tt": "band-dft-cpu", "material": "C24/qHPC24", "step": "S1_opt"},
    {"action": "advance_ready", "tt": "band-dft-cpu"}
  ]
}
```

CLI 会拒绝 `stop`、`rerun`、`clean`、`-f`、`-y` 以及标记为需要人工复核的建议。允许的动作逐条交给 `autozt act`，因此仍会进入 `.tf_agent_log.jsonl` 审计；`--dry-run` 不执行任何远端动作。

`inspect`、`plan`、`cycle` 是高层入口。它们读取同一个只读状态快照，先使用 AutoZT
已有的诊断码和依赖规则生成动作，再把需要人工判断的 FAIL 单独留在
`proposals` 中。`cycle` 默认 dry-run；`--execute` 只会调用 `autozt act` 的四个
允许的非破坏性动作，且默认不自动 retry；需要时显式加 `--include-retry`，便于 cron、脚本和 MCP 使用同一套安全边界。
执行型 `cycle` 会在跨过动作边界前再次读取状态；若 cursor 已变化，则整轮返回
`stale_plan`，不执行其中任何动作。

每个 `snapshot`、`inspect`、`plan`、`propose` 结果都带一个基于当前紧凑状态的
`cursor`。把完整计划交给 `apply` 时，AutoZT 会在执行前重新读取同一作用域的状态；
cursor 不一致就返回 `stale_plan`，不会提交任何动作。这样队列状态或另一个 agent 改变
状态后，模型会重新观察，而不会把旧计划当成事实继续执行。

`request -` 接受一个 JSON 对象，`op` 可以是 `capabilities`、`schema`、`skills`、`contract`、
`snapshot`、`inspect`、`plan`、`cycle`、`run`、`evidence`、`propose` 或 `apply`。`scope` 统一承载 `tt`、`material`
和 `status`，所有请求仍只返回一个 JSON envelope；这比让模型拼接多个可选参数更稳定。
`serve` 使用同一套请求对象，但保持进程常驻，适合只会读写 stdin/stdout 的 LLM wrapper、
脚本或编排器；坏 JSON 只影响当前行，后续请求仍会继续处理。

```json
{
  "op": "cycle",
  "scope": {"tt": "band-dft-cpu", "material": "C24/qHPC24"},
  "execute": false,
  "include_retry": false,
  "max_actions": 10
}
```

## 选择 MCP 还是 CLI

| 场景 | 入口 |
|---|---|
| Claude、ChatGPT、Kimi 等支持 MCP 的客户端 | `autozt mcp`，异常分析建议 `AUTOZT_MCP_PROFILE=workflow`；长期巡检用 `monitor` |
| 不支持 MCP 的单次调用、Shell、cron、CI、人工排查 | `autozt agent ...` 或 `autozt agent request -` |
| 需要保持进程、逐行收发 JSON 的 wrapper | `autozt agent serve` |
| 10 分钟稳定巡检 | `autozt summary --diff`，不调用模型 |
| 失败后的诊断和选择 retry/rerun | `agent propose`，必要时把单点证据交给模型 |

MCP 协议本身不会减少模型 token。节省来自三个地方：profile 少发工具定义、快照 cursor 不重复发送未变化步骤、稳定状态由 AutoZT 的规则和物理判据先筛选。真正稳定的巡检仍由 `autozt summary --diff` 完成，无变化时不调用模型。两条入口都通过只读的 `autozt list --json` 观察状态，并在入口层按材料/状态再次裁剪；若模型不能持久化 MCP cursor，CLI 的磁盘快照通常更可靠。

对于大规模项目，`inspect`、`propose` 和 `cycle` 的候选动作也受 `max_actions` 限制（默认
20，最大 20）；这只是模型上下文的上限，不会改变 AutoZT 的材料状态。先用 `-tt`、`-p`
或 `-status error,running,pd` 缩小范围，再让模型读取单点证据，通常比一次发送全量状态更省。
