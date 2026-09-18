# AutoZT 的定位与 AI 智能体接入

建议主定位：**面向热电材料、支持 AI 智能体接入的自动化计算平台**。完整 zT 工作流回答“计算什么”，2D 原生处理回答“领域规则如何落实”，智能体接口回答“人和模型如何可靠地使用这些计算能力”。

MCP 插件是 AutoZT 的一种接入形态，不是软件的整体身份。AutoZT 独立提供技能、输入生成、工作流推进、物理完成判据和结果归档；外部模型可以选择技能、解释异常并通过受控接口请求执行。没有模型时，既有自动化流程仍可运行。

## 定位边界与发展方向

目前可主张的是“支持智能体接入”（agent-ready）：技能发现、结构化契约、状态摘要、动作计划、执行前状态核对及审计已有接口基础。不能仅凭这些接口就主张“自主科学家”或完整的“AI 原生研究工作台”。契约字段齐全也不代表每个技能的物理语义和全部输入输出都经过了验证。

第三项贡献应表述为**面向 AI 智能体的科研工作流接口**，并用具体任务展示模型如何发现能力、读取所需契约、检查状态、请求动作和取得证据。MCP/CLI 是实现这一能力的接入协议，不单独作为科学创新。

若进一步发展为“目标驱动的热电研究工作台”，优先补齐并验收以下链路，而不是继续增加接口名称：

1. 研究目标到可审查的计算方案：列出所选技能、前置产物、方法假设和缺失参数；区分已有规则生成的运行期动作计划与尚需验证的跨技能研究规划。
2. 执行前的科学输入检查：结构、维度、温度/载流子网格、单位与厚度约定、外部程序和上游产物必须匹配；描述性 Schema 不能代替这些检查。
3. 面向结果的查询：按材料、性质、温度、载流子浓度和方向返回数值、单位、方法、验证状态及来源文件；不能仅返回“作业完成”。已有结果文件与 provenance 可作为基础，统一查询覆盖仍需逐技能核验。
4. 长任务的持续执行与恢复：关闭模型会话后，由已启动的 AutoZT monitor 继续推进；新会话从持久状态恢复观察。启动 MCP 或 `agent serve` 本身不等于启动后台监控。

验收应比较同一批任务在人工/脚本、通用模型直接操作命令、模型经 AutoZT 接口三种方式下的完成率、错误动作、人工介入和 token 消耗。模型、输入和任务范围需一致；没有测量前，不宣称显著降低 token 或优于所有工作流系统。

## 三层接口

```text
MCP 客户端（Claude / Kimi / ChatGPT / 本地 MCP agent）
                 │ MCP tools/list + tools/call
                 ▼
          autozt mcp（本地 stdio 适配器）

终端模型、Python、CI ── JSON ──► autozt agent request / serve
                 │
                 └──────────────► 共享 Agent 协议与服务层
                                      │
                                      ▼
                         autozt 普通 CLI + autozt act
                                      │
                                      ▼
                         skill.yaml / generators / checks
                                      │
                                      ▼
                         VASP / MACE / phonopy / HPC
```

普通 CLI 是执行核心；MCP 和 Agent CLI 是两个可替换的外部入口。MCP 客户端是否实际启动 `autozt mcp`，由客户端配置决定；仅安装 AutoZT 不会自动连接任何模型。

## 模型能调用什么

默认 `workflow` profile 面向新接入，提供少量通用工具：

```text
list_skills → describe_skill(tt) → inspect / get_snapshot
                                      ↓
                               cycle(dry-run)
                                      ↓
                               cycle(execute=true)
```

工具按动作抽象，而不是每个技能一个工具。模型可以读取：

- 技能的输入、输出、参数、步骤和依赖；
- 材料当前步骤、队列状态、集群和诊断码；
- 已变化的步骤和候选动作；
- 失败证据及是否建议 retry。

模型请求执行时，动作仍进入 `autozt act`。`stop`、`rerun`、`clean`、配置修改和强制参数不能由普通模型调用直接放行。`cycle` 默认只生成 dry-run 计划，执行型请求还会用 cursor 检查状态是否已经变化。

## 三个加强环节

### 1. 从研究目标到可审查方案

模型先调用 `research_plan`，提交目标、维度、温度网格和载流子网格。返回值列出候选技能、结构/电子输运/晶格输运/汇合阶段、缺失输入和人工复核项。它是方案草案，不会自动提交作业；用户确认后才进入 `inspect → cycle`。

```bash
autozt agent research_plan --goal "比较二维材料的 zT" --dimension 2D \
  --temperature 300 600 --carrier 1e18 2e18
```

### 2. 执行前科学预检查

结果回收后，`preflight` 只读检查结果目录、电子输运和晶格热导产物、zT 汇总产物，以及二维厚度约定。它明确区分“文件存在”和“技能物理 validator 已判定通过”；不会把描述性 Schema 冒充为收敛证明。

```bash
autozt agent preflight --result-dir /path/to/material/zt-dft-cpu/result \
  --dimension 2D --thickness 6.5 --temperature 300 600
```

### 3. 按条件查询带来源的结果

`results` 从已回收的 JSON 产物中按性质、温度、载流子浓度和方向查询，并返回每个数值的源文件和字段路径。模型不需要猜列名或把 `zT`、`κL` 和单位混在一起；返回结果仍应结合 `provenance` 和技能说明解释。

```bash
autozt agent results --result-dir /path/to/material/zt-dft-cpu/result \
  --property zt --temperature 600 --carrier 1e18
```

MCP 的 `workflow` profile 同样暴露 `research_plan`、`preflight` 和 `results` 三个只读工具。它们不增加每个技能一套工具，也不触发超算提交。

### 论文式的人机确认边界

面向用户的 agent 应把科学接口当作一段可暂停的对话，而不是连续盲调用：

```text
用户：为 MoS2 计算 300–800 K 的 zT(T,n)
Agent：我先生成研究计划，列出结构、电子输运、晶格热导和汇合步骤；当前不会提交作业。
      research_plan -> conversation_state=awaiting_confirmation
      [展示 plan_summary、plan_steps、缺失输入和 confirmation_payload]
用户：继续
Agent：执行只读 preflight 和 inspect，再给出 dry-run 动作清单；是否允许提交？
用户：允许提交
Agent：调用 cycle(execute=true)，动作经 autozt act，返回 job、集群、审计和结果来源。
```

`research_plan`、`preflight` 和 `results` 都带有统一的 `conversation_state`、`message`、
`next_action`、`requires_user_confirmation` 和 `confirmation_payload` 字段。缺少温度网格、
载流子网格或二维材料名时状态为 `needs_user_input`；计划完整但未确认时为
`awaiting_confirmation`；预检查失败时为 `preflight_review`；通过后为
`ready_to_execute`。这些字段只表达对话边界，不替代技能自身的物理 validator。

`research_plan.review_card` 是面向界面的计划卡：包含目标、编号步骤、每步输入/输出/判据、
后续工具调用顺序和 `Proceed` 状态。`action_surface` 会把普通 CLI/人工动作与 MCP 可执行动作
分开，模型不能把展示用的 `init_project` 或 `start_workflow` 误当成 `apply_actions` 的动作。

MCP 的 `content[0].text` 也会返回一行短摘要，便于只读取文本的客户端；完整证据仍在
`structuredContent.data`。CLI 则在同一 JSON envelope 中增加 `conversation`，因此 MCP 和
`autozt agent` 可以互换，模型不需要为两种传输重新学习流程。推荐同时加载
`plan-2d-zt`、`run-validated-workflow` 和 `explain-result-provenance` 三个 MCP prompt。

## 技能如何被模型发现和使用

符合 AutoZT 现有步骤生命周期的新技能，通常只需要加入：

```text
skill/<name>/skill.yaml
skill/<name>/gen_*.py
skill/<name>/checks.py（只有内置 marker 不足时）
skill/<name>/templates/ 或公共池依赖
skill/<name>/README.md
```

`skill.yaml` 的 `schema: 2` 自描述段应包含：

- `io_schema.inputs`：文件/结构/参数来源、类型、单位和必需性；
- `io_schema.outputs`：真实路径、格式、单位、产生步骤和回收方式；
- `io_schema.params`：模型可选择的参数、默认值和配置位置；
- `io_schema.steps`：每步实际使用的工具、输入、输出和验证器；
- `flow`：阶段、依赖、产生物以及可连接的下游技能；
- `corrections`：已实现的诊断处理器和人工复核边界。

这些字段是模型契约，不是自动生成器。声明不会自动解析任意科学文件，也不会自动把一个技能的输出转成另一个技能的输入。`autozt schema --strict` 只能检查结构和部分引用，仍需对照生成器、`CONF_SPEC`、模板、检查器和真实样本。

如果新技能仍能用 `start`、`retry`、`fetch`、`advance` 管理，MCP 和 Agent CLI 通常无需修改；如果引入新的生命周期、动作、状态类型、结果查询方式或风险等级，应先扩展 AutoZT 核心，再同步共享 Agent 协议和 MCP Schema。不要为每个技能写一套 MCP 工具。

## 分发和接入

软件包可以包含 MCP/Agent 入口和技能契约，但使用者仍需提供自己的：

- `autozt.yaml` 和项目根目录；
- SSH/集群账号、工作目录和队列配置；
- VASP、MACE、POTCAR、模型权重和外部程序；
- MCP 客户端配置。

MCP 本地 stdio 配置示例：

```json
{
  "mcpServers": {
    "autozt": {
      "command": "/home/user/AutoZT/bin/autozt",
      "args": ["mcp"],
      "env": {
        "AUTOZT_MCP_PROFILE": "workflow",
        "AUTOZT_CONFIG": "/home/user/AutoZT/setting/tf.yaml"
      }
    }
  }
}
```

这意味着 AutoZT 是一个“可接入的大模型工具包”，不是一个包含模型、超算账号和全部科学环境的即插即用云服务。

## 三个可验证亮点

1. **完整 zT 工作流**：一个结构可以进入电子输运和声子输运分支，最终汇合为 `zT(n,T)`，每一步有完成闸门。
2. **2D 原生处理**：维度、真空轴、面内弛豫、二维弹性单位、NAC/ZA 稳定性和输运厚度换算属于工作流规则，不只是材料标签。
3. **面向 AI 智能体的科研工作流接口**：通过契约发现能力、依据状态和证据请求动作；MCP 和 Agent CLI 负责标准化接入，AutoZT 负责物理判据、状态机、恢复语义和动作网关。模型可以替换，既有自动化执行仍可独立运行。

这三个亮点应与实验证据一起表述。MCP 本身不是科学创新，也不代表所有新技能都自动获得结果解析；它证明的是接入面稳定、可替换、可审计。
