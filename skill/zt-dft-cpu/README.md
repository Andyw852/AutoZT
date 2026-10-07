# zt-dft-cpu —— 热电优值 ZT（组合技能：ke-dft-cpu + kl-dft-cpu → ZT 汇总）

ZT = S²σT / (κ_e + κ_L)。从 v0.2（2026-10）起本技能**只有一个步骤 S20_zt**：

| 量 | 谁算 | 在哪 |
|---|---|---|
| S、σ、κ_e | `ke-dft-cpu`（S8_kappa，2D 可用 S8.4_amset2d） | `<材料>/ke-dft-cpu/` |
| κ_L | `kl-dft-cpu`（S6_kappa；也接受 kl-mlff-* S4_kappa、fit-fc-thermal S2_kappa） | `<材料>/kl-dft-cpu/` |
| ZT(n, T) | 本技能 S20_zt | `<材料>/zt-dft-cpu/` |

ke / kl 是**完整的独立技能**，在同一个材料目录下各自计算；本技能不再抄它们的步骤定义
（v0.1 抄了一份 + 符号链接，kl 段一直停在旧力常数引擎上）。**改 ke-dft-cpu 或
kl-dft-cpu（换引擎、修 bug、加步骤），ZT 自动跟着用。**

## 1. 用法

```bash
# 给一个结构：建材料目录 + 同时 init zt / ke / kl 三个技能（不提交）
autozt -tt zt-dft-cpu -p Bi2Te3 register --poscar POSCAR [--cluster jzzn|local]
# 推进：不带 -tt 时 start 会依次推进这个材料的每个技能
autozt -p Bi2Te3 start
# 看三个技能的任务图（S20_zt 在 ke/kl 结果拉回前是 WAIT，写着在等哪份结果）
autozt -p Bi2Te3 graph
```

已有材料（之前单独跑过 ke-dft-cpu / kl-dft-cpu）：`autozt -tt zt-dft-cpu -p <材料> init` 后
直接 `start` S20_zt。电子/晶格参数改上游技能本身（`autozt -tt ke-dft-cpu -p M -j S conf --set …`）；
本技能的参数只有 `step20_zt/step.conf` 的 `KTEMP_MODE` / `KL_CONST_T`。

## 2. 跨技能结果依赖（needs_results）

`skill.yaml` 的 S20_zt 声明要哪些结果（路径相对**本材料的本地目录**，按顺序取第一个合格的）：

```yaml
needs_results:
  transport: [ke-dft-cpu/result/step8_amset/transport.json, ke-dft-cpu/result/step8.4_amset2d/transport.json]
  kappa_L:   [kl-dft-cpu/result/step6_kappa/kappa_summary.json, kl-mlff-cpu/..., kl-mlff-gpu/..., fit-fc-thermal/...]
  structure: [ke-dft-cpu/result/step1_opt/workflow_method.txt]
```

- 依赖没满足时 S20_zt 是 WAIT，任务图 / 状态表写明在等哪一份。
- gen 时 autozt 把选中的结果目录推到远端 `zt-dft-cpu/inputs/<键>/`，附 `source.json`
  （来自哪里、来源戳、每个文件 sha256）；S20 **只读 inputs/**。
- 想换 κ_L 来源的优先级：项目层覆盖 `needs_results` 的顺序即可。

## 3. 同名材料 / 旧结果不会被误读

- autozt 每次把某技能某步的结果拉回本地时，在结果目录写 `.autozt_origin.json`：
  材料目录（绝对路径）、材料 POSCAR 的 sha256、技能、步骤、远端目录、时间。
- `needs_results` 只认：**同一个材料目录** + **当前 POSCAR 指纹一致** 的结果。
  另一个项目里同名的材料、改过结构之后留下的旧结果、没有来源戳的旧版本结果，都不算
  （旧版本结果重新 `autozt -tt <技能> -p M fetch` 即可补上来源戳）。
- S20 再核一次：三份输入的 `source.json` 必须是同一个材料目录、同一个 POSCAR 指纹，
  否则拒绝汇总。不再去远端按目录名找"兄弟技能"（同名材料会串）。

## 4. 产物：`zt_summary.json` 关键字段

| 字段 | 含义 |
|---|---|
| `temperatures` / `doping_cm-3` / `doping_type` | 电子段的温度与掺杂网格（负 = n 型） |
| `rows[i]` | 第 i 个掺杂档：`seebeck_uV/K`、`sigma_S/m`、`kappa_e_W/mK`、`kappa_L_W/mK`、`kappa_tot_W/mK`、`PF_W/mK2`、`ZT`（按温度对齐的列表） |
| `grid_ZT` | [掺杂][温度] 的 ZT 矩阵 |
| `peak_ZT` | n / p 型的峰值 ZT 及其温度与掺杂（外加当时的 κ_L） |
| `kappa_L` | κ_L 来源、温区、xx/yy/zz、元胞口径标记、2D 厚度信息 |
| `notes` | 口径说明（约化方式、插值方式、不外推区间…） |

配套产物：`zt_summary.txt`（人读汇总 + 每档明细表）、`zt_vs_T.png`、`zt_vs_doping.png`、`zt_components.png`。

## 5. 历史记录（v0.1：在本技能里重抄 ke/kl 两段的版本）

下面的验证与修复记录针对 v0.1 的装配方式（本技能目录里的符号链接 + 抄写的 24 个步骤），
保留作为方法学记录；物理口径（张量约化、κ_L 插值、元胞口径闸门）在 v0.2 不变。

### 5. 验证记录（2026-09-18）

| 项 | 结果 |
|---|---|
| `autozt skills` / `schema --strict` | 认到技能；退出码 0；四模式展开 24/17/8/1 步 |
| 资源解析自检（autozt 自身 `find_asset`/`step_conf_sources`） | 四种模式**均 0 缺件** |
| 单测 `tests/test_zt_dft_cpu.py` | 18 项全过（公式/单位/口径/插值/取峰/解析/装配） |
| **3D 实跑（Si，autozt 全链路）** | `S20_zt` 在 jzzn 登录节点生成并回拉：`zt_summary.json` 43 KB + 3 张图；n 型峰值 ZT=0.394 @800 K/1e20，p 型 0.301 |
| 电子段实机（Si，经本技能装配） | S1_opt/S2.1_scf/S2.15/S2.2_pbe(+plot)/**S2.3_hse(4/4)**/S2.3_hseplot/S5_dielect/S6_elastic/S7_deform(4/4)/S7.1_read 全部 OK；**S3_uniform OK（WAVECAR 24.4 MB，26³）**、S4_wave OK、S8_kappa_e 已提交（上游 4 处缺陷修复后） |
| 晶格段实机（Si，经本技能装配） | SK1_opt/SK2_static/SK3_nac/**SK4_disp(11/11)**/**SK5_fc OK（stable, min_freq=0.000 THz）**；SK5.1_plot 与 SK6_kappa 已提交 |
| **MACE 交叉对照（Si）** | 同材料 `kl-mlff-gpu` 的 κ_L(300K)=**91.27 W/mK** vs 本技能 DFT phono3py 的 **88.8 W/mK**（差 3%）——两法在 300 K 一致；用 `KTEMP_MODE=const300` 汇总得 n 型峰值 ZT=0.162 @900 K、p 型 0.129（高温柔性差异来自"常数 κ_L 近似"，不是两法分歧） |
| S8 预检 in-job 预跑（2026-09-18） | 在 `amset_clean` 里于 `step8_amset/` 直跑 `overlap_preflight.py --in-job` → **结论：通过（有告警），exit=0**；带窗口比较按补丁 05 跳过，弹性 Christoffel 最小特征值 47.07 正定。两条 WARN 为设计内信息（3D+非完整网格去对称化、AMSET 0.4.19 G 平移） |
| 起跑前预检（2026-09-18） | S8_kappa_e 的提交命令调 `overlap_preflight.py` 而材料目录里没有（上游 WIP 漏声明）→ 已在技能侧补软链 + `gen_need`；SK6_kappa（自包含内联汇总）与 SK5.1_plot（脚本/依赖/目录探测均齐）预检通过 |
| **2D 实测（真实材料 P1_Mo-MoS2_…_Mo2S3）** | DIM=2D 正确识别；张量按**面内 (xx+yy)/2** 约化（σ 2734/289.9→1512 S/m、S −234.9/−167.5→−201.2 µV/K、κ_e 面内平均 0.00637）；κ_L 用元胞口径 `kappa_xx_yy_zz`（1.721/0.846→1.284 W/m/K，zz≈0 剔除）；**元胞口径闸门 c=20.0 Å vs Lz=20.0 Å 通过**；n 型峰值 ZT=0.379 @800 K |

#### 5.1 全网格对照（`wavefunction_full`）A/B 实跑结论（2026-09-21）

**A 臂**（默认，`step4_wave` 去对称化 h5）——已完成并归档：

- 集群产物 `result/step20_zt/`：`zt_summary.json` 43 KB + `zt_summary.txt` + 3 张图；峰 **n 型 0.3163 @800 K/1e21**、**p 型 0.2386 @800 K/1e20**；
  `sources` = 本技能 `step8_amset/transport.json` + `step6_kappa/kappa_summary.json`；
  `amset_settings = {scattering_type:[ADP,IMP], bandgap:1.0913, interpolation_factor:10}`。
- 归档（sha256 前 12 位）：`tmp/zt_arm_archive/A_step4_wave_transport.json` = `edb0023a9432`、
  `A_zt_summary.json` = `3c06c55c4515`、`kappa_summary_selfcomputed.json` = `3d915e8aa7e8`。
- 复核：当前 `result/step20_zt/zt_summary.json` 与归档逐字段比对，**值差异 0 处**（仅多出
  `amset_settings.doping/temperatures` 两个 provenance 字段）。

**B 臂**（`WAVEFUNCTION_FULL = true` → `step4b_wave_full` 的 26³ 全网格 h5）——**端到端未完成**：

- 两次提交：`3850661` NODE_FAIL；`3852782` 跑满 24 h 后
  `CANCELLED AT 2026-09-20T16:51:19 DUE TO TIME LIMIT`，`transport.json` 未更新（仍是 A 臂文件）。
- 成本实测（AMSET 0.4.19，24 核，`qos=regular` = 1 天硬顶）：日志 `Interpolating spin-up bands 2-6`；
  band 1 = 25 min、band 2 = 2 h41、band 3 = **7 h41**、band 4 到 42% 时已 4 h35，**band 5/6 未开始**。
  单带耗时随带序显著上升 → 24 h 内不可能完成；按此斜率即使 `qos=premium`（2 天）也不够。
- 并行参数已专项核对，**无低垂果实**：`NCORE=1`、`PRECFOCK=Fast` 已是最优；
  `KPAR=8` 在 96 核 / `NBANDS=324` 约束下已是合法上限（KPAR 层并行饱和）。
- ★ 操作提醒：项目 `templates/step8_amset/step.conf` 里的 `WAVEFUNCTION_FULL = true` 目前**仍停在 B 臂设置**；
  今后若要让 S8 回到 A 臂（去对称化）路线，需先把它翻回 `false`。

**科学裁定沿用 V26，不由本次端到端产出**（`skill/ke-dft-cpu/step8.4_amset2d/VERIFICATION.md` §11，2026-09-18）：

> 三维裁定（Si 全网格对照）已完成：同一 vasprun / settings / nworkers，只切 h5 ——
> ADP **1.01 / 0.92**、overall **1.02 / 0.97**、IMP **0.58 / 0.77**；重复跑**逐位相同**。
> 结论：三维真实重叠可用（不必改 unity）。

**结论**：`wavefunction_full` 分支的**装配与单变量开关**已端到端验证到"能产出全网格 h5 并让 S8 走
`from_data`"（S3b/S3c/S4b 全 OK、日志 0 次 `Desymmetrizing`）；**B 臂的 ZT 数值未由本技能端到端产出**，
其裁定引用 V26。若要补齐，需先把 AMSET 这一段搬上更长墙钟的 qos 或做插值规模 benchmark，两者都要先请示。

### 5.x 实跑中暴露并已修复的上游缺陷（2026-09-18，补丁在 `tmp/zt_upstream_patches_20260918/`）

用户授权「按提案修」后，以下 4 处上游缺陷已按「先 `git apply --check` → 应用 → `py_compile` → 集群实跑验证」的流程修复：

| 补丁 | 技能/文件 | 问题 | 实跑验证 |
|---|---|---|---|
| 01 | `kl-dft-cpu/kl_common.py` | 抽帧校验把平衡帧 `disp-00000` 也算进位移帧 → 末尾越界 `IndexError`（最后一次提交引入，**全仓范围**） | `-j SK5_fc init` → `gen 完成` |
| 02 | `ke-dft-cpu/step3_uniform/gen_step5_uniform.py` | `DK_MAX/UNIFORM_NMAX` 不是 `step.conf` 键，护栏提示的覆盖方式无效 | `DK_MAX_3D=0.08` → 网格 34³→**26³=17576**，`gen 完成` |
| 03 | `ke-dft-cpu/step2.3_hse_plot/gen_step4.1_plot_band.py` | 缺 `import os`（第 410 行 `os.path.join`）→ `NameError`，连带挡死 S8_amset | `-j S2.3_hseplot start` → 出 `band_summary.json` |
| 04 | `ke-dft-cpu/step8_amset/gen_step10_amset.py` | 非极性体系 ε_static ≡ ε_inf 被当 DFPT 失效拦下（下游本会按物理口径剔除 POP） | `-j S8_kappa_e start` → `gen 完成` 并提交 AMSET 作业 |
| 05 | `ke-dft-cpu/step8.4_amset2d/overlap_preflight.py`（在建工作） | "AMSET 实际插值带窗口"读不到时用兜底默认 11–17 去比 → 把任何窗口≠11–17 的材料（Si=2–6）拦死 gen 与 in-job | 行为级：verdict 由 `error` 变 `ok`，打印"运行日志缺失 → 跳过窗口一致性比较" |
| 06 | `ke-dft-cpu/skill.yaml` + 本技能同一声明 | S2.155 的 `done_marker` 与脚本真实落点不一致 → 状态表长期 error | marker 改 `../step2.15_discriminant/discriminant.json`（该 json 在集群上确实存在） |
| 07 | `ke-dft-cpu/step5_dielect/validate_dielectric.py` | 单元素非极性体系 ε_s≡ε_∞ 被误判 FAIL（与补丁 04 同一物理） | 行为级：单元素 `ok=True`+notes；双元素对照仍 `ok=False`；集群实测 Si 产物 `ok=True` |
| 08/09/10 | `validate_dielectric.py` / `decide_discriminant.py` / 两个 `skill.yaml` | 两个 `run:gen` 步的 `done_marker` 永远找不到：**步骤目录没人创建**，且校验脚本默认把 json 写进 `step5_dielect/`（不是本步目录） | 集群实测：`-j S5.1_dievalid` 与 `-j S2.155_discr` 都「已生成」并回拉 |

补丁 03/04 是执行中新暴露的（不在最初提案内），已逐条向用户披露。以下两处**未修**（不影响流程推进）：

1. **`S2.155_discr` 长期显示 error**：`decide_discriminant.py` 把 `discriminant.json` 写进
   `step2_bandgap/step2.15_discriminant/`，而该步目录是 `step2.155_discriminant_decide` → `done_marker` 找不到（下游直接读 2.15 的 OUTCAR，不读它，无实际影响）。
2. **`S5.1_dievalid` 对非极性体系报 FAIL**：`ε_s ≡ ε_∞`（Si 无 IR 活性声子的正确结果）被校验器要求人工确认。
3. **`S8_kappa_e` 被在建的「重叠预检」误报挡死**（2026-09-18 实测，属**未提交的在建工作**
   `ke-dft-cpu/step8.4_amset2d/overlap_preflight.py`）：它把「AMSET 实际插值带窗口」的**兜底默认
   11–17** 当成真值，与本材料 `step4_wave` 的真实窗口 **2–6**（`amset.log` 原文 `Including bands 2—6`）
   比较 → 报 `★ 拦截：… 带窗口不一致`。gen 阶段和作业内都会命中（运行时 `amset.log` 尚不存在，
   同样拿兜底值）。**本技能不擅自改这个在建文件**：当前 Si 测试改用 `electronic: false`，
   S8 等上游修好（读不到运行日志时应跳过该比较，或直接从 h5/step4_wave 日志取窗口）后再打开。

以下为修复前的原始记录（保留以便对照）：

1. **`ke-dft-cpu` S3_uniform 成本护栏挡住小胞材料**：`gen_step5_uniform.py` 的
   `DK_MAX_3D=0.06` 对 Si 原胞（|b|≈2.0 Å⁻¹）要求 34³=39304 点，超过脚本里的
   `UNIFORM_NMAX=20000` → 直接 `sys.exit`。**而 `DK_MAX`/`UNIFORM_NMAX` 是模块常量、
   不在 `SPEC` 里**，所以脚本自己提示的“在项目里覆盖 DK_MAX”改不动（项目 step.conf 写了会被
   `strict=False` 静默忽略）。脚本注释给的参考量级是 “Si 0.08 → 25³ ≈ 1.6e4”。
   本技能不绕过它：Si 测试改用 `electronic: false`。
2. **`kl-dft-cpu` S5_fc 抽帧校验 off-by-one**：`kl_common.check_frames_match_displacements`
   对 `sorted(glob("disp-*"))` 取 {0, mid, last} 三个下标去索引 `phono3py_disp.yaml` 的
   `displacements`；但 rattle 流程生成的目录是 `disp-00000…disp-00010`（11 个），
   位移表只有 10 条 → `disps[10]` 抛 `IndexError`，gen 直接失败（该步没有开关可跳过）。
   Si 测试因此也把 `lattice: false`：κ_L 取同材料 kl-dft-cpu 的既有结果。

两处修好后，把项目的 `electronic` / `lattice` 改回 `true` 即可跑完整 24 步。

### 5.y 依赖

`requires: python [pymatgen, numpy, matplotlib, amset, BoltzTraP2, phonopy, phono3py, spglib]`、
`conda: amset_clean`、`exe: [vasp_std, vaspkit, phono3py]`（与上游两段一致）。
`S20_zt` 本身只用标准库；没有 matplotlib 时自动跳过画图，JSON/TXT 照常产出。
