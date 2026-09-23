# AutoZT-v1.0 —— 隔离开发副本（ISOLATED DEV FORK）

> **【2026-09-17 状态变更 · 先读这段】** 本文件描述的"副本 / 生产仓在别处"的格局**已作废**：
> 用户已明确 **`~/software/AutoZT` 是唯一权威仓库**，旧的 `~/software/taskflow-v2.0`、
> `~/software/taskflow` 不再维护（AutoZT 已把它们的近期改动同步过来）。
> 下面第 1 节的 `auto_advance: false` / `auto_watch: false` **仍然保留**，但理由变了：
> 旧仓的后台监控目前还在跑（实测 v2.0 的 bin/tf monitor 有 5 个进程、crontab 里还有
> taskflow/monitor.sh），此时开自动推进会与它们抢提交、重复投作业；等旧仓监控停掉后再按需打开。
> 第 0 节的来源表、第 2 节的 `AutoZT-v1.0` / `AutoZT-v2.0` 路径均属历史记录，不必再照做。

> 本目录 **不是** 生产仓库。生产仓库是 `~/software/AutoZT-v2.0`（用户明确要求：
> "不要修改目前这个版本"）。这里是从 v2.0 复制出来的隔离副本，**所有新特性开发都在
> 本副本内进行**，v2.0 一行不改。

## 0. 来源与基线

| 项 | 值 |
|---|---|
| 复制来源 | `/home/wangchao/software/AutoZT-v2.0` |
| 复制时间 | 2026-09-14 |
| 上游提交 | `56fc99f`（feat(fc-fit): add the fc-fit skill） |
| 复制方式 | `rsync -a`（保留当时工作区 91 个未提交改动，原样带过来） |
| 基线快照提交 | `baseline: 隔离快照 …`（分支 `dev-v1.0`） |
| 排除内容 | `tmp/`（672M 临时区）、`__pycache__`、`setting/.tf_state_cache.json`、`setting/.tf_watch.log`、`setting/.tf_hung.json`、`setting/.tf_summary_*.txt`、`setting/.tf_watch.pid` |

git 远端已重命名为 **`origin-v2.0`**（避免误 push 到 GitHub 上的 v2.0 仓库）。
要 push 到新仓库时先 `git remote rename origin-v2.0 origin && git remote set-url origin <新地址>`。

## 1. ★ 安全约定（重要）

副本的 `setting/tf.yaml` 与 v2.0 **共用同一批 project_roots**（同一批真实项目、
同一批超算账号）。也就是说：**在副本里敲错命令，会真的往集群提交/取消作业。**

因此副本已被改成：

- `auto_advance: false` —— 副本内任何 autozt 命令都**不会**自动提交作业；
- `auto_watch: false` —— 副本绝不拉起后台监控（避免与 v2.0 的 monitor 抢提交）。

开发本副本时只允许跑**只读命令**：
`autozt skills` / `autozt schema` / `autozt history` / `autozt list` / `autozt summary` / `autozt -V` / `autozt config`。
**禁止**在副本里跑 `autozt start/stop/retry/rerun/clean`（会动真集群）。

## 2. 怎么跑副本

```bash
cd ~/software/AutoZT-v1.0
python3 bin/autozt skills          # 只读：列出技能（走副本自己的 skill/ 与 setting/）
python3 bin/autozt -V
```

`~/.local/bin/autozt` **仍指向旧仓库 `~/software/AutoZT/versions/v1.0/tf`**，未被本副本改动。

## 3. 本副本做了什么（v1.0 加技能友好化）

| 计划 | 状态 | 内容 |
|---|---|---|
| 建议 1.3 | ✅ 已完成 | `io_schema` 段（技能自报吃什么/吐什么/有哪些旋钮）+ 校验 |
| 建议 1.2 | ✅ 已完成 | `flow` 段（整条流程 + 产物能喂给谁） |
| 建议 1.1 | ✅ 已完成 | `_corrections/` handler 库（base + 4 handler + `autozt correct`） |
| 建议 3.1 | ✅ 已完成 | 零配置自动发现已在，本次补齐自描述的自动校验/展示 |
| W1–4 | ✅ 已完成 | `autozt schema [技能]`（`--json` / `--strict`） |
| W5–8 | ✅ 已完成 | `history.jsonl`（采集时自动记录）+ `autozt history` |
| 模板 | ✅ 已完成 | `skill/_template/`（复制即上线，实测零配置被发现） |
| **第二批 P0-3** | ✅ 已完成 | 每步声明式 I/O（`io_schema.steps`）+ `autozt skill show` 技能卡片（论文图 2 数据源） |
| **第二批 P0-2** | ✅ 已完成 | per-step provenance（gen 时自动落档）+ `autozt prove`（含 `--verify`） |
| W9–12 / P0-1/4/5 / P1-6…9 | ⏳ 未做 | 缓存、`autozt act` 网关、`autozt bench`、筛选报告、session export、评测集——见 ROADMAP 第 5.3 节 |

详见 `V1.0-ROADMAP.md`。

## 4. 离线自测（不碰集群，随时可跑）

```bash
cd ~/software/AutoZT-v1.0
python3 tmp/test_v1_skillspec.py     # 自描述校验 + 每步 I/O + 卡片渲染 + 纠错库匹配/隔离
python3 tmp/test_v1_history.py       # history 记录/过滤/--since
python3 tmp/test_v1_correct_cli.py   # autozt correct / autozt diagnose 接入（假数据）
python3 tmp/test_v1_prov.py          # provenance：真跑一遍 gen（本地 bash 模式）+ autozt prove
# 只读命令冒烟（都不 ssh）：
python3 bin/autozt skills ; python3 bin/autozt skill ; python3 bin/autozt schema --strict
python3 bin/autozt skill show kl-mlff-gpu        # 技能卡片（论文图 2）
# 本地沙盒主路径（tmp/tf_smoke.yaml 只指向仓库自带 test/tf_test，不碰真集群）：
python3 bin/autozt -c tmp/tf_smoke.yaml list
python3 bin/autozt -c tmp/tf_smoke.yaml history
```

四个测试脚本都是纯本地断言（不连超算、不写项目、不提交作业），失败返回非零。
`tmp/test_v1_prov.py` 会把 hpc 置空走本地 bash，真的执行一次 gen 并检查落档。

**真机验证记录**（jzzn 4 核 Si 小体系，隔离配置 `tmp/tf_jzzn_si.yaml`）：
`tmp/v1_si4c_verification.md`——含 4 核改动位置、作业号、provenance 证据链、
实测暴露的 6 个缺陷与修法、复现命令。
