# AutoZT（autozt）—— 面向热电材料的自动化计算平台

[![CI](https://github.com/Andyw852/AutoZT/actions/workflows/ci.yml/badge.svg)](https://github.com/Andyw852/AutoZT/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> 一句话：从 POSCAR 到 `zT(n,T)`，在多个超算上自动运行 VASP / MACE 热电计算，并让脚本、人工或 AI 智能体通过同一套可验证接口协作。

> **v2.0 重构**：本目录是重构版——原 7463 行单体脚本已拆成 `autozt/` 包（`bin/autozt` 为入口），代码导航见 [`CONTEXT.md`](CONTEXT.md)。原始单体保留在 `versions/v1.0/tf` 作对比基准。

## 能干什么

AutoZT 的三个主亮点是：

- **完整 zT 工作流**：从一个 POSCAR 串起结构优化、电子输运、声子/晶格热输运，汇合得到 `zT(n,T)`。
- **2D 原生处理**：自动识别维度并沿整条链路处理真空轴、面内弛豫、二维单位、ZA 稳定性和厚度换算。
- **面向 AI 智能体的科研工作流接口**：模型按需发现技能、读取输入输出契约、依据状态和证据请求动作；通过 MCP 或 Agent CLI 接入，物理判据、状态机、恢复和审计由 AutoZT 执行核心负责。
- **三步科学协作闭环**：`research_plan` 生成可审查方案，`preflight` 检查跨步骤输入，`results` 按条件返回带来源的结果；三者均只读，不绕过物理 validator。

- **多材料**：一个项目根下几千个结构（如 `C20/qHPC20`、`Ag/qHPC20_Ag1C20_s0`）统一管理，逐材料多步流水线（弛豫 → 静态 → 后处理）自动推进。
- **多技能**：技能从 `skill/*/skill.yaml` 动态发现，覆盖 VASP 能带/弹性/电热导/晶格热导/结构优化，以及 MACE 同类、声子、MLFF 训练和替代模型；遵循现有工作流抽象的新技能通常不增加 MCP 工具。
- **多超算**：jzzn（CPU 真 SLURM）、a800（A800 GPU 真 SLURM）、3090（无 SLURM 的 fakeslurm 垫片服务器），换超算只改一个 `hpc` 名。
- **全自动**：`auto_advance` + `autozt monitor` 后台监控，作业算完自动拉结果、自动提交下一步；挂死作业自动 `scancel`+续跑（`hang_check`）。
- **省心巡检**：`autozt summary --diff` 无变化输出 0 字节，有变化才吐几行——适合 AI / cron 定时巡检。
- **三种接入面**：MCP 默认使用 `workflow` 小工具面（完整面需显式设置 `AUTOZT_MCP_PROFILE=full`）；一次性脚本用 `autozt agent request -`；常驻 wrapper 用 `autozt agent serve` 的 JSONL。三者共用技能契约、快照和动作审计。
- **低上下文巡检**：长期监控可用 `AUTOZT_MCP_PROFILE=monitor autozt mcp`，只暴露 `get_snapshot`/`cycle` 两个通用工具；异常时再切换到 `workflow`。
- **过期计划保护**：`inspect`/`propose`/`snapshot` 返回状态 `cursor`；`apply_actions` 或 `autozt agent apply` 执行前重新核对，状态变化就拒绝旧计划并要求重新观察。

## 快速开始

```bash
# 1. 安装（把 autozt 放进 PATH）
git clone <repo> ~/software/AutoZT
ln -s ~/software/AutoZT-v2.0/bin/autozt ~/.local/bin/autozt   # 确保 ~/.local/bin 在 PATH（.bashrc/.profile 里）
autozt --version

# 2. 配置（全局 autozt.yaml：ssh 别名、项目根、技能 work_dir、auto_advance）
vi ~/software/AutoZT/setting/tf.yaml

# 3. 初始化一个材料项目（生成 project_setting）
cd 你的材料根目录 && autozt -tt opt-mlff-cpu init

# 4. 开跑 + 后台监控
autozt auto on                      # 开全局自动推进
autozt -tt opt-mlff-cpu monitor -d    # 后台监控：自动拉结果 + 自动提交下一步

# 5. 巡检
autozt summary --diff               # 首选巡检：无变化 0 字节
autozt -tt opt-mlff-cpu summary     # 只看某技能

# 给模型或脚本的稳定 JSON 接口（不需要 MCP）
autozt agent capabilities
autozt agent schema
autozt agent skills
autozt agent snapshot
autozt agent propose > plan.json
autozt agent apply plan.json --dry-run

# 新接入推荐：一次读取/规划；cycle 默认 dry-run
autozt agent inspect
autozt agent cycle --dry-run
# 明确需要时才执行非破坏性动作
autozt agent cycle --execute
autozt agent run --dry-run          # 确定性规则闭环，不需要 MCP/LLM
autozt agent run --execute          # 显式执行允许的安全动作

# 也可以用一个 JSON 请求代替命令行参数
echo '{"op":"inspect","scope":{"status":"error"}}' | autozt agent request -

# 多次调用时使用常驻 JSONL，逐行输入、逐行输出
printf '%s\n%s\n' '{"op":"capabilities"}' '{"op":"snapshot"}' | autozt agent serve
```

## 常用命令与接入入口

`autozt --help`（或 `autozt help`）显示简明帮助；`autozt --help-all` 显示高级命令、全部参数和示例。本仓库入口为 `bin/autozt`，可用 `python3 bin/autozt --help` 检查新版；`autozt --version` 可核对 PATH 实际指向的版本。

| 用途 | 推荐命令 | 行为 |
|---|---|---|
| 只读巡检 | `autozt summary --diff` | 无变化静默；总表用 `autozt list` |
| 提交 | `autozt -tt TT -p MAT start` | 生成缺失输入并提交 |
| 保留产物重生成 | `autozt -tt TT -p MAT retry` | 不提交，检查后执行 `start` |
| 删除后重生成 | `autozt -tt TT -p MAT rerun` | 删除旧步骤产物，不提交 |
| 只删除 | `autozt -tt TT -p MAT clean` | 不重新生成 |
| 自动推进 | `autozt -tt TT -p MAT auto on` | 开启并推进就绪步骤，保留取消标记 |
| 恢复取消步骤 | `autozt -tt TT -p MAT auto resume` | 开启推进并清除指定技能的取消标记，仍等待依赖 |
| 持续监控 | `autozt monitor -d` | `--stop` 停止，`--restart` 重启 |

`watch` 兼容 `monitor`；`restart`、`watch restart`、`monitor restart` 兼容 `monitor --restart`。旧入口保留，但新脚本统一使用推荐写法。`auto resume` 必须同时指定 `-tt` 和 `-p`，只在明确需要恢复取消步骤时执行，不放进周期巡检。`status` 会拉结果且可能自动提交，不属于只读查询。

## 三个超算

| 集群 | ssh 别名 | 类型 | 适用技能 |
|---|---|---|---|
| jzzn | `jzzn` | CPU，真 SLURM（cpu192） | VASP 各技能 + MACE CPU 类 |
| a800 | `A800` | A800 GPU，真 SLURM | GPU 类技能 |
| 3090 | `user_3090` | 8×RTX3090，无 SLURM（fakeslurm 垫片） | MACE GPU 类（注意：垫片无排队，需靠 max_jobs 节流） |

## 技能清单

技能数量和版本以 `autozt skills` / `autozt agent skills` 为准；每个技能的模型契约位于
对应的 `skill/<技能>/skill.yaml`，包括输入、输出、步骤依赖、参数和纠错策略。MCP
只暴露通用动词，不为技能或步骤单独增加工具。

## 大规模实战（2926 材料 · opt-mlff-cpu · jzzn）

一次跑完 2926 个金属掺杂富勒烯结构（38 金属 × 各笼型 × 各位点）的 MACE 结构优化 + 形成能，经验记录如下，供新项目参考：

- **本地模式（v3）**：输入文件以本地项目目录为准，超算只是算力；每个项目自己的 `project_setting/`。
- **每技能并发上限 `max_jobs`**：大体系批量提交靠它节流（`task_types.<技能>.max_jobs`），到上限后剩余任务先本地生成输入待命，自动补交——是正常状态，不是故障。
- **共享 project_setting（体系级布局）**：几千材料时用一份 `project_setting`（`local_root` 指向项目根），避免每材料一份配置把扫描拖慢；对应地 `autozt` 已修复共享布局下 `-p <材料>` 定位。
- **大体系采集自动分块**：`autozt` 会按 `AUTOZT_COLLECT_CHUNK`（默认 500）把同组材料分块 ssh 采集，避免 "Argument list too long"。
- **2D 结构真空轴规范**：工作流要求 2D 结构真空沿 c 轴（第 3 个晶格矢量），生成结构时注意；已算过的若真空在别的轴，需旋转晶格（保持笛卡尔坐标）再重跑。
- **形成能参考化学势 μ**：`scripts/mlff_mu/` 提供本地算 38 金属单质每原子能量的固化脚本（`setup_local.sh` + `run_mu.sh`），产出 `MU = ...` 单行填进 step3 配置。

## 关键修复记录（本仓库已含）

1. 共享（体系级）布局下 `-p <材料>` 定位失效 —— 已修复（预过滤保留共享段）。
2. `max_jobs` 并发计数把"已完成作业残留状态（OTHER）"误计入 —— 已修复（`_BUSY_KINDS` 只数 R/PD）。
3. 大体系采集 ssh 命令超长（Argument list too long）—— 已修复（`AUTOZT_COLLECT_CHUNK` 默认 500 自动分块）。
4. 3090 fakeslurm 垫片无排队会瞬间拉起几千作业压垮服务器 —— 该平台务必配 `max_jobs` 节流，或改用真 SLURM 集群（jzzn/a800）。

## 提交代码到 GitHub

改完代码，一键「脱敏 → 提交 → 推送」到远程（本地真实版 + 远程脱敏版双轨，真实用户名/超算路径不会泄露）：

```bash
bash scripts/tf-git-push.sh "你的提交说明"
```

- 敏感词映射在 `setting/git-sanitize.conf`（本地私有、不进仓库），新增敏感词往里加一行「源词=目标词」
- github 走 `ssh.github.com:443`（绕过 22 端口封锁），带 8 次重试防卡

## 文档

- 完整手册：[TASKFLOW.md](TASKFLOW.md)
- 技能开发与模型契约：[总手册第 7 章](TASKFLOW.md#7-技能开发规范原-skill_devmd-内容)，其中 7.15–7.18 说明输入输出、报错、新技能适配和验收边界。
- 模型调用与服务生命周期：[总手册第 9 章](TASKFLOW.md#9-给大语言模型用agent-接入)；支持本地 MCP 的客户端使用 [MCP](docs/mcp.md)，终端模型/脚本使用 [Agent CLI](docs/agent-cli.md)。两条路径可独立使用，`agent serve` 仅供高级 wrapper 集成。
- 软件定位与模型接入：[AutoZT 的定位与 AI 智能体接入](docs/llm-plugin.md)
- AI/监控接入规范：[AGENTS.md](AGENTS.md)
- 各技能细节：`skill/<技能>/README.md`
