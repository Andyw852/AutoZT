# <集群名> 接入 AutoZT —— 部署与验收记录（模板）

> 用法：复制本文件为 `docs/<HPC>-DEPLOY-<YYYYMMDD>.md`，按实测填写；不要留 <占位符>。
> `setting/*` 属本地文件（.gitignore）；本记录涉及个人路径时脱敏为 `<work>/...`、`<user>`。
> 通用步骤与坑见 `docs/DEPLOY-NEW-CLUSTER.md`；本文件只记录**这台集群的事实**。

---

## 0. 概要

| 项 | 值 |
|---|---|
| 集群 / 单位 | <如 SCNet 西安 / 自建> |
| 集群 id / 调度器 | <id / SLURM 22.05 / PBS ...> |
| 账号 | <user>（密钥到期 <YYYY-MM-DD>） |
| ssh 别名 | <alias>（端口 <port>；<直连 / 经 <jump> 跳板>） |
| 登录节点 | <hostname>，<n> 核 / <mem> GB |
| 完成日期 / 部署人 | <YYYY-MM-DD> / <name> |

## 1. 资源实测（全部自己验，不要只抄管理员口径）

### 1.1 分区可用性

| 分区 | 可提交 | 节点数 | 核/节点 | 内存/节点 | GRES | MaxTime |
|---|---|---|---|---|---|---|
| <partition> | ✅/❌ | <n> | <c> | <m> GB | <gpu:N / 无> | <time> |

> 判据：逐分区 `sbatch --test-only`；❌ 记下报错原文（如 `Invalid account or account/partition combination`）。
> 批量探测时每分区 sleep 1–2 s，避免 `sbatch: request too frequently` 把可用分区误判成不可用。

### 1.2 账户限额（实测反推）

| 项 | 实测值 | 验证方法 |
|---|---|---|
| 并发上限（CPU 核） | <N> | 超配提交报错 `group max tres limit ...` |
| 单作业上限 | <N 核 / M 节点> | `-N/-n` 逐档 `--test-only` |
| 任务数上限 | <无 / N> | 实提 N 个 1 核作业观察 |
| 单作业最长 | <X 天> | 模板 `--time` 上限 |
| 存储配额 | <work X TB / home Y GB> | `df -h`、`quota` |

### 1.3 调度行为实测（决定 AutoZT 并发策略）

- 整节点作业：<是否要等完全空闲节点；观察到的排队特征>
- 同时提交 n 个作业：<实际起几个、排队原因>
- 结论：<默认 32 核/作业、32 核 × 10 并发、...>

## 2. 软件与路径

| 组件 | 版本 | 路径（远端） | 备注 |
|---|---|---|---|
| 编译器 | <gcc 9.3.0 / Intel 2021> | module: <name> | |
| MPI | <HPC-X 2.11 / IntelMPI 2017> | <path> | 与编译器 ABI 匹配 |
| BLAS/LAPACK | <MKL / 系统 3.4.2> | <path> | 是否补 dgegv/dggsvd |
| FFT | <FFTW3 double> | <path> | |
| VASP 原版 | <6.4.3> | <.../bin> | std/gam/ncl |
| VASP 变胞版 | <同版本+optcell> | <.../bin> | 是否启用 relax_2d |
| POTCAR | <potpaw_PBE_64: N 元素> | <...> | 布局 A/B/C |
| VASPKIT | <1.2.4> | <...> | ~/.vaspkit 用空格分隔 |
| Python 环境 | <venv py3.11 / conda> | <...> | 关键包版本 |
| 可选 QE | <6.8 module> | <...> | |
| 可选 amset | <0.5.1 独立 env> | <...> | |
| 可选 MACE | <torch x.y + mace-torch z> | <...> | 模型清单 |
| 可选其它 | <COGITO / UniHamGNN> | <...> | |

## 3. 部署步骤与踩坑（只记本集群特殊项）

| # | 现象（报错关键行） | 根因 | 处理 | 验证 |
|---|---|---|---|---|
| 1 | <error line> | <原因> | <命令/改动> | <怎么确认好了> |
| 2 | | | | |

> 通用坑见 DEPLOY-NEW-CLUSTER.md 第 13 节；这里只记本集群特有的（特定位分区名、
> 某型号节点的 MPI 绑核问题、特定 module 冲突等）。

## 4. AutoZT 配置摘要（脱敏）

```yaml
name: <HPC>
ssh_host: <alias>
work_dir: <work>/software/AutoZT/work
partition: <partition>
max_cpus: <每节点核数>
conda_env: <work>/software/AutoZT/venvs/autozt
potcar_dir: <work>/software/AutoZT/vasp_pseudopotentials
remote_path_prefix: <work>/software/AutoZT/pybin:<work>/software/AutoZT/venvs/autozt/bin
```

- 提交模板：`setting/<HPC>/templates/`（核数 <N>，time=<...>，mpirun 参数 <...>）
- 与其它集群的差异：<必须 --bind-to none / 必须 module load / 需要跳板 / ...>

## 5. 验收（必须全部有输出或数字）

| # | 检查 | 结果 |
|---|---|---|
| 1 | ssh -o BatchMode=yes <alias> true | ✅ |
| 2 | 非交互 python：ssh <alias> 'python3 -V' | ✅ <版本> |
| 3 | autozt doctor | ✅ |
| 4 | autozt check --cluster <HPC> | ✅ <k>/25 |
| 5 | 冒烟 VASP（Si 2 原子，ENCUT 250） | 收敛 ✅，E=<...> eV（差 <...>） |
| 6 | 并发实测 | <n> 作业 × <m> 核结果 |
| 7 | 作业清理 | squeue 0 个 ✅ |

## 6. 已知限制 / TODO

- <例：技能默认 192 核，用 step.conf [submit] ntasks_per_node 覆盖>
- <例：某可选组件未装 / 某分区不可用>
- <例：跳板 ControlMaster 过期需重建>

## 7. 变更记录

| 日期 | 变更 | 影响 |
|---|---|---|
| <YYYY-MM-DD> | 初始部署 | — |
