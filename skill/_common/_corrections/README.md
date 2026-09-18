# _corrections/ —— 纠错 handler 库

> AutoZT-v1.0 建议 1.1。目标：**把"这个报错该怎么办"从人的脑子里搬进仓库**，
> 让新成员/新技能不必读完 AGENTS.md 决策表 + 问作者。

## 1. 它解决什么

老流程：FAIL → 人/AI 看 diag → 翻 AGENTS.md 第五节决策表 → 猜 → retry/rerun。
问题是：决策表是散文，各技能的"专属坑"散在 16 份 README 里，**加新技能时没人知道
要写什么**，也没有机器可判读的"报错→动作"映射。

handler 库把这件事变成**可枚举的代码**：

| 诊断特征 | handler | 风险 | 动作 |
|---|---|---|---|
| ZBRENT / EDDAV / SCF 空转 | `zbrent` | review | retry；连续失败则打 INCAR 补丁（AMIX/BMIX→ALGO=All→NELM≥200） |
| force not converged / relax_nsw | `nsw_continue` | safe | retry 续算（自动 cp CONTCAR） |
| NODE_FAIL / 被抢占 | `node_fail` | safe | retry 重交 |
| dir missing | `dir_missing` | destructive | 只建议；rerun 须人工执行 |

## 2. 怎么用

```bash
autozt diagnose -p C24/qHPC24              # 诊断里直接带 corrections 建议（只读）
autozt correct  -p C24/qHPC24 -j S1_opt    # 只打印命中的 handler + 建议命令（只读）
autozt correct  -p C24/qHPC24 -j S1_opt -y # 执行 handler.apply（改远端 INCAR，自动备份）
autozt schema   band-dft-cpu               # 看某技能声明了哪些 corrections
```

★ `autozt correct` **永不提交作业、永不删目录**：apply 只允许改输入文件（备份 + 原子写），
提交仍走 `autozt start`、重生成仍走 `autozt retry/rerun`。

## 2.1 已知修复与未知错误的分流

`autozt diagnose -p MAT [-j STEP] --json`（以及 Agent CLI 的 `evidence`）在每个失败步
返回 `repair` 分类：

| 分类 | 含义 | 下一步 |
|---|---|---|
| `known/auto_safe` | handler 明确声明 `auto: true` 且风险为 `safe` | 可由受控自动策略选择；仍须由调用方决定是否提交 |
| `known/approval_required` | 已知规则，但要改 INCAR/配置或影响计算 | 给出补丁/命令，用户批准后执行；执行后再检查并 `retry`/`start` |
| `unknown/llm_review` | 没有可靠 handler，或已有建议不足以自动执行 | 把 `llm_handoff` 证据交给模型；模型只能提出方案，不能直接改文件/提交 |

`llm_handoff` 是有界证据包，不是自动执行请求。它包含技能、材料、步骤、诊断码、诊断文本、
主机和作业号，并要求模型给出最小修改、验证和重交方案。新增 handler 前先证明匹配特征
足够窄；宁可转人工/LLM 复核，不要把未知文本加入宽泛匹配。

每条 `Suggestion` 还可以返回一份跨技能通用的修复计划：

| 字段 | 作用 |
|---|---|
| `changes` | 结构化声明要改什么（目标、路径、操作、值）；只是计划，不等于已经写入 |
| `verify` | 修改或重交后的验收条件 |
| `rollback` | 备份、恢复或停止继续尝试的方式 |
| `requires_approval` | 是否必须用户批准；默认按风险等级保守处理 |
| `execute_via` | 正式入口，例如 `autozt retry -> autozt start` |

模型只能读取这些字段并提出方案，不能把 `actions` 或自然语言命令自行转成文件修改。
只有已注册 handler 的 `apply()` 能执行受控输入修正；提交、取消和重生成仍必须经
`autozt` 命令和审批网关。技能 YAML 中的 `match/repair/verify/rollback` 也只是契约说明，
不会直接执行。

当前 `auto_safe` 只表示“规则允许安全处理”，并不绕过 AutoZT 的 `act` 网关，也不意味着
计算结果已经收敛。自动策略若要上线，还必须设置重试次数、批量提交上限和审计记录；本模块
本身只负责分类、建议和受控输入修正。

## 3. 怎么写一个新 handler（照 zbrent.py 复制改 20 分钟）

```python
# skill/_common/_corrections/my_case.py
from .base import CorrectionHandler, Suggestion

class MyCaseCorrection(CorrectionHandler):
    name = "my_case"                  # ★ 唯一名，skill.yaml 的 corrections: [my_case] 引用它
    title = "一句话说明"
    matches = ("特征子串1", "特征子串2")   # 小写子串，命中 diag/diag_code/日志即算命中
    steps = ("S3", "step3")            # 适用步骤前缀；("*",) = 所有步骤
    risk = "review"                    # safe | review | destructive

    def suggest(self, ctx):
        return Suggestion(
            handler=self.name, title=self.title, risk=self.risk,
            reason="为什么这么判",
            actions=["人看的动作清单"],
            commands=[self.retry_cmd(ctx), self.start_cmd(ctx)],
            changes=[], verify=["写出该 handler 的验收条件"],
            rollback=["写出备份或停止继续尝试的方式"],
            requires_approval=True,
            execute_via="autozt retry -> autozt start")

    # 可选：真要改输入文件才实现（写前必须备份 + 原子替换）
    # def apply(self, ctx, cfg=None): return rc, ["改了啥"]

HANDLER = MyCaseCorrection()
```

写完 `autozt schema <技能>` 就能看到它（`autozt diagnose` 会自动命中，不注册也行）。

## 4. 放哪 & 优先级

| 位置 | 作用范围 |
|---|---|
| `skill/_common/_corrections/` | **全局**：所有技能都能命中（推荐放这里） |
| `skill/<技能>/_corrections/` | **技能私有**：只对该技能生效，同名可覆盖全局 |

技能 `skill.yaml` 里写 `corrections: [名字, ...]` 表示"本技能声明会用到这些纠错"——
这只是**自描述**（`autozt schema` 展示给人看），不写也照样能被全局兜底命中。

## 5. 纪律（写 handler 前必读）

1. `apply()` 只改输入文件，**必须**先备份再原子写（参考 `autozt._hung_incar_fix`）；
2. **禁止**在 handler 里 sbatch / scancel / rm / mv 目录——提交与取消永远只走 autozt；
3. 风险等级别乱标：`safe` 才允许无人值守；`destructive` 只许给建议（`dir_missing` 就是范例）；
4. 判据要窄：宁可漏判（人再判断），不可误判（把正常推进的作业当成故障处理）。
