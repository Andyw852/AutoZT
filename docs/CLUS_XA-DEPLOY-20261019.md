# clus_xa（SCNet 西安）接入与 VASP 6.4.3 部署记录

日期：2026-10-19　执行环境：eshell111（SCNet E-Shell 节点，用户 wangchao217）

## 1. 集群与配额

| 项 | 值 |
|---|---|
| 集群 | `clus_xa`（NSCC 西安，Slurm 22.05.8，节点 Hygon C86 7285H） |
| 登录 | `ssh eshell111`（eshell111.hpccube.com:65082，直连可达） |
| 可用分区 | **仅 `xahctest`**：433 节点 × 64 核 / 255 GB，**纯 CPU**，最长 3 天 |
| 账户配额 | 并发 **320 核**（`-N 5` 通过，512 核报 group max tres limit） |
| 存储 | `/work` 43 PB、`/public` 554 TB；家目录 `/work/home/wangchao217` |

### 1.1 并发与限额（2026-10-19 实测）

| 项目 | 值 | 依据 |
|---|---|---|
| 同时最大核数 | **320 核**（账户级 GrpTRES cpu=320） | `-N 5 -n 320` 通过；`-N 8 -n 512` 报 `group max tres limit 320 for account wangchao217` |
| 单作业最大核数 | **320 核** | 同上 |
| 节点规格 | 64 核 / 255 GB（默认 3.8 GB/核） | `scontrol show partition xahctest`：`TotalCPUs=27712 TotalNodes=433` |
| 最大整节点数 | **5 个**（320 / 64） | 账户并发上限换算 |
| 同时任务数 | **无 per-user 限制** | QOS/账户 `MaxJobsPU`、`MaxSubmitPU` 为空；实测 37 个 1 核任务并存 |
| 单作业时长 | 最长 **3 天** | `MaxTime=3-00:00:00` |
| 分区 MaxNodes | UNLIMITED（受 320 核额度约束） | `scontrol show partition` |

按每任务核数换算的并发任务数：1 核 → 320 个；8 核 → 40；16 核 → 20；32 核 → 10；**64 核（整节点）→ 5 个**。

> 实测注意：整节点（64 核）作业必须等**完全空闲**的节点；2026-10-19 分区只有约 652 空闲核
> （`alloc 26690 / idle 652`）时，同时提 6 个整节点作业只有 1 个立刻 RUNNING，其余
> `PENDING (Resources)`。想更快排上建议用 32 核左右的任务，或先 `sinfo -p xahctest` 看空闲节点。

## 2. 远端目录（`~/software/AutoZT/`，约 1.5 G）

```
vasp/vasp.6.4.3-pristine/bin/{vasp_std,vasp_gam,vasp_ncl}   # 官方版
vasp/vasp.6.4.3-optcell /bin/{vasp_std,vasp_gam,vasp_ncl}   # IOPTCELL 补丁版
vasp_pseudopotentials/potpaw_PBE_64/                       # 346 元素 POTCAR
venvs/autozt/      # py3.11: numpy1.26.4 scipy matplotlib spglib pymatgen2024.8.9 seekpath
pybin/vaspkit      # VASPKIT 1.2.4 软链；~/.vaspkit 用空格分隔格式配 PBE_PATH
work/              # AutoZT 项目根
```

## 3. 编译要点（gcc 9.3.0 + HPC-X + 系统 BLAS/LAPACK/FFTW）

1. gcc 7.3 在 `ml_reader.F` 触发内部编译器段错误 → 换 `module load compiler/gcc/9.3.0`。
2. FFTW 的 `fftw3.f` 需要 `INCS=-I/public/software/mathlib/fftw/3.3.8/double/gnu/include`。
3. VASP 6.4.3 链接变量是 `LLIBS/BLASPACK`（不是 BLAS/LAPACK），按 `arch/makefile.include.gnu` 写。
4. LAPACK 3.8 已删除 `dgegv_`（6.4.3 仍调用）→ 从 netlib 取回 `dgegv.f`/`dggsvd.f` 编进 `libdmy`。
5. BLAS/LAPACK/FFTW 用系统 `/usr/lib64`（3.4.2 / 3.3.2），运行期只需 `module load compiler/gcc/9.3.0`。
6. VASPKIT 1.2.4 的 `~/.vaspkit` 必须用空格分隔（`KEY value`）；写成 `KEY = value` 会静默不生成 POTCAR。

## 4. AutoZT 配置

- `setting/clus_xa.yaml`（`setting/*` 被 .gitignore 忽略，与其它集群一致）
- `setting/clus_xa/templates/`：`submit_std_{0d,2d,3d}.tpl`、`submit_ncl_{2d,3d}.tpl`、
  `step1_opt/submit_std_2d.tpl`（1 节点 64 核、最长 3 天、单节点 MPI 走共享内存）
- 远端额外：`~/.ssh/environment` 追加 PATH（AutoZT 采集需要 `python3`），
  `~/bin/python3` 指向 venv python

## 5. 验收结果（2026-10-19，测试材料 `Si_clxa`）

| 步骤 | 结果 |
|---|---|
| `autozt -tt band-dft-cpu check --cluster clus_xa` | ✅ 可运行 |
| Si 冒烟 SCF（8 核） | ✅ 12 步收敛，E = −5.4167 eV |
| optcell 版冒烟 + `IOPTCELL = 1 1 0 1 1 0 0 0 0` | ✅ 输出 `Constraining stress: 1 1 0 | 1 1 0 | 1 1 0` |
| S1_opt 三段式弛豫 | ✅ job 71901987 COMPLETED（7:42，64 核），RELAX_OK |
| S2_static 自洽 | ✅ job 71902847 COMPLETED（0:29） |
| S3_WAVECAR 波函数 | ✅ job 71903364 COMPLETED（1:03），WAVECAR 61 MB |
| S3.1_plot 能带 | ✅ 已出图；PBEsol 带隙 **0.4589 eV（间接）**、最小直接隙 2.508 eV、209 路径 k 点 |
| S3.2_vac 真空对齐 | ✅ vacuum_align_summary.json |
| S4_HSE | ⏸ 按项目 `BANDGAP = pbe` 关闭（原因见下） |

**结论：PBE 能带全链路（S1→S2→S3→S3.1→S3.2）在 clus_xa 上全部通过。**

### 未跑 HSE 的原因（技能默认值不适用于 64 核小节点）

`skill/band-dft-cpu/gen_step4_HSE.py:113\` 的 `NTASKS_PER_NODE_LARGE = 192`：当自洽 k 点 +
能带路径点总数 >100 时会硬写 `--ntasks-per-node=192`，而 xahctest 节点只有 64 核，
sbatch 直接报 `CPU count per node can not be satisfied`。脚本支持用 step.conf `[submit]`
覆盖，但 AutoZT 推送的合并 step.conf 不含该节，故未生效。

建议（需改 skill 源码，待批准）：
1. `gen_step4_HSE.py` 的 192 改为“集群可提供核数”封顶（如读 step.conf 的 `[submit]` 或
   `max(64, ...)`）。
2. `gen_step3.1_plot_band.py` / `gen_step4.1_plot_band.py` 的 `--prefix` 默认值
   `band-dft-cpu` 与 `skill.yaml` 的 `done_marker: band_summary.json` 名字不一致
   （当前靠 check 的通配匹配才判 OK，建议统一）。

## 6. 复现命令

```bash
ssh eshell111
cd ~/software/AutoZT/vasp && bash build3.sh      # 或按 §3 手动 make
cd ~/software/AutoZT && ./bin/autozt -tt band-dft-cpu check --cluster clus_xa
```