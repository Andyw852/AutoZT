# device-thermal —— 器件级传热仿真（技能介绍，v0.1）

> 把「材料筛选」推进到「真实器件热评估」的最后一环：用上游实测 κ + 界面热导 + 器件几何
> + 焦耳热源，算出器件热点温升与温度场。

---

## 1. 它回答什么问题

给定一个 2D 材料 FET 结构（沟道 / 栅介质 / 衬底）和它的工作功率：

- 热点温升 **ΔT** 多少？
- 有效热阻 **R_th** 多少？
- 温度在沟道和氧化层里怎么分布（二维温度场）？
- 升温要多久（热时间常数 **τ**，快/慢两个尺度）？
- 界面热导 **G** 从哪来？（可回灌 GPUMD/MD 或实验值）

覆盖材料：MoS2 / WS2 / MoSe2 / WSe2 / graphene / hBN 作沟道；SiO2 / Al2O3 / HfO2 / hBN 作栅介质。

**不覆盖**：需要弹道/波动效应（沟道尺度 ≲ 声子平均自由程）、非傅里叶、强非线性材料的场景，
那需要 S6_bte（OpenBTE/JAX-BTE，尚未安装）。

---

## 2. 算力与运行环境（重要）

| 部分 | 环境 | 规模 | 实测 |
|---|---|---|---|
| S1–S4 器件求解 | CPU 登录节点，**纯 numpy，不依赖 scipy** | 未知数 41×62 ≈ 2542 | 四步端到端 **wall ≈ 18 s**（含瞬态 2000 步） |
| S1 / S2 / S4 单步 | 同上 | JSON 读写 | < 0.1 s |
| S3 稳态 | 同上 | 红黑 SOR 到 1e-10 | 4443 次迭代 |
| **S5/S6（可选 GPU）** | GPU 主机 + GPUMD | 原子数 × 步数 | 1280 原子 × 320000 步 ≈ **485 s**（3090） |

**结论**：器件部分几乎不耗算力（登录节点秒级、不排队）；只有想用 MD 第一性拿界面热导时才需要 GPU，
成本随「原子数 × 步数」线性涨。

---

## 3. 输入（详细）

> 参数三处可写，优先级从低到高：技能 `templates/step.conf` → 项目 `templates/` → 项目 `step.conf`。
> 改：`autozt -tt device-thermal -p <材料> -j <步骤> conf --set params.KEY=值`（改配置，先请示）。
> 看合并后的最终值：`autozt -tt device-thermal -p <材料> -j <步骤> conf`（只读）。

### 3.1 器件几何（S1 读，写进 device_spec.json）

**`DEVICE_NAME`** · 字符串 · 默认 `device`
- 只作报告与出图标签，不进任何物理计算。
- 例：`--set params.DEVICE_NAME=MoS2_FET_A`。

**`CHANNEL_MATERIAL`** · 字符串 · 默认 `MoS2`
- 沟道材料在 `device_materials.json` 里的键，决定 4 个物性：面内 κx、跨面 κy、厚度 t、体积热容 ρc_p。
- 可取：`MoS2` `WS2` `MoSe2` `WSe2` `graphene` `hBN` 等；写错 S2 直接报错并列出可用键。
- **关键交互**：此键只提供「材料库兜底值」。当 `KAPPA_SOURCE=auto` 且上游产物存在时，κx/κy 与厚度会被
  上游实测覆盖，但 **ρc_p 永远来自材料库**（上游不提供热容）。
- 影响：κy 与 t 决定沟道跨面热阻 t/κy（→ ΔT）；t 与 ρc_p 决定 τ_fast。

**`CHANNEL_THICKNESS`** · float · 单位 **m** · 默认空
- 沟道厚度。**优先级最高**，设了就覆盖上游厚度与材料库厚度。
- 缺省回退链：本参数 → 上游 `thickness_2d.json` → 材料库 `thickness_m`；**三者都没有则 S2 直接报错退出**。
- 影响：同时进跨面热阻 t/κy 与热容 ρc_p·t，所以 ΔT 和 τ 都随它变。
- 例：`--set params.CHANNEL_THICKNESS=6.7177e-10`（= 6.7177 Å = 0.672 nm）。

**`OXIDE_MATERIAL`** · 字符串 · 默认 `SiO2`
- 栅介质/氧化层材料键，决定其 kx/ky/厚度/ρc_p。
- 可取：`SiO2` `Al2O3` `HfO2` `hBN` `SiC` 等。
- 影响：氧化层热阻 t/kx 是 R_th 的一项；其热容主导慢模态 τ_RC。

**`OXIDE_THICKNESS`** · float · 单位 **m** · 默认 `300e-9`
- **仅当**材料库该材料的 `thickness_m` 为 null 时才使用（S2 里的兜底）。材料库已带厚度时此参数无效。
- 影响：直接进 R_th（t/kx）与 C_total（→ τ_RC）。
- 例：`--set params.OXIDE_THICKNESS=90e-9`（90 nm 栅氧）。

**`SUBSTRATE_MATERIAL`** · 字符串 · 默认 `Si`
- 底部热沉材料名。**当前物理解里底边界恒为 T_AMB，与衬底材料无关**，此键只写进报告作标签。
- 想按衬底材料建热阻，需改 `device_common.solve_device`（属于改技能源码，先请示）。

**`L_CHANNEL`** · float · 单位 **m** · 默认 `200e-9`
- 沟道长度（沿 x，源漏方向），决定计算域长度与网格 dx = L/(N_GRID_X−1)，默认 dx = 5 nm。
- 影响：与 Q 一起定总功率 power_uW = Q·L·W；也影响 x 方向的温度分布与接触散热相对权重。

**`W_DEVICE`** · float · 单位 **m** · 默认 `1e-6`
- 器件宽度。**不进传热方程**（求解按单位宽度组装），只用于把结果折算成总功率 power_uW = Q·L·W。
- 例：`--set params.W_DEVICE=2e-6`。

### 3.2 热源与边界（S1 读）

**`Q_JOULE`** · float · 单位 **W/m²** · 默认 `1e8`
- 焦耳功率密度（按器件面积计），只加在沟道层的 2 个网格上。
- 换算：`1e8 W/m² = 100 mW/µm²`；`1 mW/µm² = 1e9 W/m²`。
- 影响：**ΔT 的第一驱动量**。1D 网络下 ΔT ≈ Q·R_pp，近似线性。
- 典型：2D FET 自热约 10–500 mW/µm² → 1e7–5e8 W/m²。
- 例：`--set params.Q_JOULE=5e8`。

**`T_AMB`** · float · 单位 **K** · 默认 `300.0` · **双重作用**
- ① S3 里作底边界固定温度（`CONTACT_BC=sink` 时左右接触也固定为它）；
- ② **S2 里用它从上游 κ 的 `temperatures` 列表里挑最接近的一行**。
- 影响：改它同时改边界温度和取用的 κ。例如上游算了 200/300/400 K，设 350 会取 400 K 那行。
- 例：`--set params.T_AMB=350`。

**`CONTACT_BC`** · `sink` | `insulator` · 默认 `sink`
- `sink`：左右金属接触为恒温（Dirichlet = T_AMB）→ 散热好、ΔT 低，对应良好接触。
- `insulator`：左右接触绝热（Neumann）→ 热量只能走底部，ΔT 高，是**最坏上限**。
- 代码里 `sink`/`fixed`/`contact` 都判为恒温，其余按绝热处理。
- 建议：同一器件两个都算，给 ΔT 的上下限。例：`--set params.CONTACT_BC=insulator`。

### 3.3 物性来源与界面（S2 读）

**`KAPPA_SOURCE`** · `auto` | `upstream` | `literature` · 默认 `auto`
- `auto`/`upstream`/`dft`：优先用上游 `kl-dft-cpu` 实测 κ（面内 = 0.5(κxx+κyy)，跨面取 κzz）与上游厚度；
  找不到就打 warning 回退材料库。
- `literature`：完全用材料库，不去找上游（复现文献口径时用）。
- 其它值：S2 直接报错。
- 例：`-j S2_props conf --set params.KAPPA_SOURCE=literature`。

**`UPSTREAM_SKILL`** · 字符串 · 默认 `kl-dft-cpu`
- 上游技能目录名；MLFF 路线可改 `kl-mlff-gpu` / `kl-mlff-cpu`。

**`UPSTREAM_STEP`** · 字符串 · 默认 `step6_kappa`
- `kappa_summary.json` 所在的步骤目录名。

**`TBC_OVERRIDE`** · float · 单位 **W/m²K** · 默认空
- 覆盖界面热导 G。空 → 查库里的 `沟道材料|氧化层材料` 键（没有则用 `default`）。
- 这是**把第一性/实验 G 灌进器件模型的正规入口**（例如 GPU 步 S6 算出的 G）。
- 影响：1/G 是 R_pp 的一项，G 越小界面温降越大、ΔT 越大。
- 例：MoS2/SiO2 G = 14 MW/m²K → `--set params.TBC_OVERRIDE=1.4e7`。

### 3.4 数值求解（S3 读）与绘图（S4 读）

**`N_GRID_X`** · int · 默认 `41`
- 沿沟道 x 的节点数；dx = L_CHANNEL/(N_GRID_X−1)。默认 200 nm / 40 = 5 nm。
- 影响：温度场 x 方向分辨率与迭代耗时（未知数 = N_GRID_X × (2 + N_GRID_SUB)）。
- 网格收敛证据见第 10 节（n_sub 15/30/60/120 → 19.294/19.714/19.934/20.046 K）。

**`N_GRID_SUB`** · int · 默认 `60`
- 氧化层 y 方向格数；氧化层 dy = t_oxide/N_GRID_SUB（默认 300 nm / 60 = 5 nm）。
- **沟道层固定 2 格**（dy = t_ch/2 ≈ 0.34 nm），保证界面前后都有节点。
- 影响：界面附近梯度分辨率；dT 对它不敏感。例：`--set params.N_GRID_SUB=120` 做收敛检查。

**`DO_TRANSIENT`** · bool · 默认 `true`
- 是否算瞬态升温曲线与 τ。设 `false` → summary 里没有 tau_* / heating_* 字段，S3 更快。
- 瞬态实现：几何时间步 dt₀ = τ_fast/100、每步 ×1.2，最多 2000 步，直到 t > 3·τ_RC；
  后向欧拉（无条件稳定）+ 每步热启动 SOR。实际步数随材料/几何自适应，不固定为 2000。

**`SOR_OMEGA`** · float · 默认 `1.7`
- 红黑 SOR 松弛因子，要求 1 < ω < 2。这类各向异性泊松问题 1.7 附近最优。
- 太接近 2 → 振荡/发散；太小 → 迭代次数暴涨。

**`SOR_TOL`** · float · 单位 **K** · 默认 `1e-10`
- 收敛判据：最大网格更新量 < tol。对本问题很严（默认网格 4443 次迭代）。
- 快速试算可放宽到 `1e-8`（dT 变化远小于 1e-6 K）。

**`SOR_MAXIT`** · int · 默认 `200000`
- 迭代上限（保险丝）。达到上限会用当前解继续，实际次数记在 `sor_iterations`。

**`SOR_*` 与网格参数从哪读**：S3 用自己的 step.conf；若 S3 没显式设置，会回退到 S1 固化的
`device_spec.json` 里的值。**改数值参数请指到 `-j S3_solve`**，改几何/热源指到 `-j S1_spec`，
改物性/界面指到 `-j S2_props`。

**`MAKE_PLOT`** · bool · 默认 `true`
- S4 是否输出 PNG。无 matplotlib 时自动跳过并在报告里记 `plot_error`，仍正常写 JSON。

### 3.5 可选 GPU 步参数（`TBC_*`，需先开 `optional.gpu_tbc=true`）

**路径与设备**

| 参数 | 默认 | 说明 |
|---|---|---|
| `TBC_STRUCTURE` | `tbc_interface.xyz` | 界面结构 extxyz（第二行须含 `Lattice=...`）。放技能目录或该材料的 `project_setting/skill_dir/`，随 S5 的 gen_need 推到 GPU 主机；也可直接给 GPU 主机上的绝对路径。传输轴跨度须 ≥ 20 Å |
| `TBC_NEP_MODEL` | 无（必填） | GPUMD 势文件在 **GPU 主机**上的路径（支持 `~`）。NEP 会自动核对结构元素是否都在势的元素表里 |
| `TBC_GPUMD_BIN` | 无（必填） | `gpumd` 可执行文件在 GPU 主机上的路径 |
| `TBC_CUDA_VISIBLE_DEVICES` | 空 | 指定显卡，写进 `run_gpumd.sh`。**强烈建议用 UUID**（3090 的 GPU4 故障会让数字索引错位） |

**几何与界面判定**

| 参数 | 默认 | 说明 |
|---|---|---|
| `TBC_AXIS` | `z` | 传输/堆叠方向 `x`/`y`/`z`。分组、`compute_chunk` 的 bin 方向、界面判定都按它走 |
| `TBC_PBC` | 空 | 覆盖 pbc（如 `T T F`）；空 = 沿用结构文件。写法为 3 个 T/F 字母 |
| `TBC_INTERFACE_COORD` | 空 | 界面沿 `TBC_AXIS` 的坐标（Å）。空 = 按元素组分沿轴突变自动判定（取界面键中心）。共格界面且需严格比较时建议显式给值 |
| `TBC_SOURCE_THICKNESS` | `5.0` | 源恒温器区厚度（Å），从坐标一端算起 |
| `TBC_SINK_THICKNESS` | `8.0` | 漏恒温器区厚度（Å），从另一端算起。太薄控温不稳，太厚会吃掉本就不长的体区 |
| `TBC_SOURCE_SIDE` | `high` | 源组在坐标高侧还是低侧（`high`/`low`），漏组在对侧 |
| `TBC_VACUUM_GAP_A` | `4.0` | 判真空的间隙阈值（Å），仅在组分法失败、回退到最大间隙法时使用 |

**热力学与积分**

| 参数 | 默认 | 说明 |
|---|---|---|
| `TBC_T` | `300.0` | 平均温度（K）；也作 NVT 平衡的目标温度与初速温度 |
| `TBC_DELTA_T` | `40.0` | **半温差**（K）：源 = T+ΔT、漏 = T−ΔT，实际源漏温差 = 2·ΔT（默认 80 K）。太大放大非线性，太小信噪比差 |
| `TBC_COUPLE` | `100.0` | 恒温器耦合强度，越大越强行把区域拉到目标温度 |
| `TBC_THERMOSTAT` | `heat_lan` | `heat_lan` / `heat_nhc` / `heat_bdp`。本技能的热流累积方案用 `heat_lan` |
| `TBC_TIME_STEP` | `1.0` | MD 时间步（fs）。须保证能量守恒，轻元素/高温可能要更小 |
| `TBC_SEED` | `0` | 0 = 不固定随机初速（不可复现）；非 0 → `velocity T seed <seed>`，同种子可复现 |
| `TBC_BIN_WIDTH` | `2.0` | T(轴) 剖面 bin 宽（Å），直接决定两侧各有多少拟合点 |

**三段积分（平衡 → 烧入 → 测量）**

| 参数 | 默认 | 说明 |
|---|---|---|
| `TBC_EQUIL_STEPS` | `20000` | 第一段 NVT 平衡步数；0 则跳过 |
| `TBC_BURN_STEPS` | `100000` | 第二段烧入（施加 heat_* 但不测量），释放初始应变。**强烈建议 > 0**：早期版本设 0 时源漏热流差达 3 倍 |
| `TBC_RUN_STEPS` | `200000` | 第三段测量步数；恒温器累积传能在这一段内计数（写 `compute.out`） |
| `TBC_SAMPLE_INTERVAL` | `10` | compute 采样间隔 |
| `TBC_OUTPUT_INTERVAL` | `100` | compute 输出间隔。**每行时间 = 两者相乘 × TBC_TIME_STEP**（默认 10×100×1 = 1000 fs/行） |
| `TBC_THERMO_INTERVAL` | `1000` | `dump_thermo` 间隔 |

**S6 后处理与质量判据**

| 参数 | 默认 | 说明 |
|---|---|---|
| `TBC_FIT_MARGIN_A` | `2.0` | 拟合体区时离界面与恒温器各留的余量（Å）。S6 会额外在 1/2/4 Å 上各算一次做敏感性检查 |
| `TBC_MIN_BINS` | `3` | 每侧最少拟合 bin 数，不够直接判 poor |
| `TBC_STEADY_TOL` | `0.2` | 源/漏热流相对差容差，超了判 poor（未达稳态） |
| `TBC_MIN_WINDOW_FRAC` | `0.2` | 稳态窗口最短占比；S6 在后缀里找最长的自洽窗口 |
| `TBC_DT_TOL` | `0.3` | ΔT_i 外推相对误差容差，超了判 poor |
| `TBC_MARGIN_TOL` | `1.5` | G 对拟合余量敏感度 max/min 上限，超了判 poor（引线太短/非扩散区） |
| `TBC_KAPPA_A_W_mK` | 空 | 可选，给出材料 A 的 κ 做交叉核对（`kappa_crosscheck`） |
| `TBC_KAPPA_B_W_mK` | 空 | 同上，材料 B。两者**不参与 quality 判据** |

---

### 3.6 上游输入（`kl-dft-cpu`，**可选**）

**`kappa_summary.json`**
- `temperatures`: [K] 温度列表；
- `kappa_xx_yy_zz`: [[κxx, κyy, κzz], ...] 每条对应一个温度。
- S2 取最接近 `T_AMB` 的一行：面内 κ = 0.5(κxx+κyy) 覆盖 κx；κzz > 0 时覆盖 κy。

**`thickness_2d.json`**
- `thickness_d_A`（Å，也接受 `thickness_A` / `effective_thickness_A`）→ S2 转成 m 覆盖沟道厚度。

**查找链**（相对材料目录）：`<材料>/<skill>/<step>/` → `<材料>/<skill>/result/<step>/` → `<材料>/<skill>/`。
两个文件都**可选**：缺失时 auto 模式回退材料库并打 warning；`KAPPA_SOURCE=literature` 时根本不查。

---

### 3.7 内置材料库 `device_materials.json`

**材料字段**（`materials.<键>`）：

| 字段 | 含义 |
|---|---|
| `kappa_inplane_W_mK` | 面内热导率 κx (W/mK) → 进 x 方向导热 |
| `kappa_cross_W_mK` | 跨面热导率 κy (W/mK) → 进跨面热阻 t/κy |
| `thickness_m` | 材料厚度 (m)，可为 null（此时用 step.conf 的 *_THICKNESS） |
| `rho_cp_J_m3K` | 体积热容 (J/m³K) → 决定热容与 τ |
| `kind` / `note` | 分类与备注 |

- **14 种材料**：MoS2、WS2、MoSe2、WSe2、graphene、hBN、SiO2、Si、SiC、Al2O3、HfO2、Au、Cu、Al
- **8 组界面 TBC**：`tbc_W_m2K.<"A|B">` = `{G: W/m²K, ref: 出处}`；7 个材料对 + `default` 兜底
- 来源：MoS2 面内 κ = 20.8 W/mK 取本仓库 `kl-dft-cpu` 实测（RTA raw, 300 K）；
  MoS2/SiO2 G = 14 MW/m²K 取 Yalon et al., Nano Lett. 17, 3429 (2017)。**文献值只作量级兜底**。

---

### 3.8 GPU 步额外输入

- 界面结构 extxyz（`TBC_STRUCTURE`）
- GPUMD 势文件（`TBC_NEP_MODEL`，GPU 主机本地路径）
- `gpumd` 可执行文件（`TBC_GPUMD_BIN`）

---

### 3.9 参数优先级与交互（最容易踩的 6 条）

1. **沟道厚度**：`CHANNEL_THICKNESS` > 上游 `thickness_2d.json` > 材料库 `thickness_m`；都没有 → S2 报错。
2. **沟道 κ**：上游 `kappa_summary.json`（`KAPPA_SOURCE=auto/upstream`）> 材料库。
3. **ρc_p 永远来自材料库**（上游不提供热容）。
4. **界面 G**：`TBC_OVERRIDE` > 库里的材料对 > `default`。
5. **`T_AMB` 同时决定边界温度和上游 κ 的取值温度**，改之前想清楚。
6. **按步骤设参数**：几何/热源/来源 → `-j S1_spec`；物性/界面 → `-j S2_props`；网格/SOR → `-j S3_solve`；
   绘图 → `-j S4_report`；GPU → `-j S5_tbc` / `-j S6_tbc_post`。S3 的数值参数若自身没设，会回退 S1 固化的值。

---

## 4. 输出

| 步骤 | 产物 | 关键字段 |
|---|---|---|
| S1 | `step1_device_spec/device_spec.json` | 固化后的几何/热源/边界 |
| S2 | `step2_props/thermal_props.json` | channel/oxide 的 kx/ky/厚度/ρc_p + `tbc_W_m2K` + 各自 source |
| S3 | `step3_solve/device_thermal_summary.json` | **dT_peak_K**、**R_th_eff_m2K_per_W**、`dT_peak_1D_network_K`、`sor_iterations`、`T_profile_mid_K`、瞬态 `tau_fast_s`/`tau_RC_s`/`tau_63_s`/`heating_t_s`/`heating_dT_K` |
| S3 | `step3_solve/T_field.json` | 二维温度场 `T[Nx][Ny]` + `x_nm`/`y_nm` |
| S4 | `step4_report/device_report.json` | 汇总报告 |
| S4 | `step4_report/device_thermal_Tmap.png` | 温度场图（无 matplotlib 则跳过，仍出 JSON） |
| **S5**（可选） | `step5_tbc/{model.xyz, potential.txt, run.in, run_gpumd.sh, tbc_inputs.json}` | NEMD deck；跑出 `compute.out`/`compute_chunk.out`/`thermo.out` |
| **S6**（可选） | `step6_tbc_post/tbc_result.json` | `G_W_m2K`、`R_interface_m2K_W`、`q_W_m2`、`dT_interface_K`、`quality`、`quality_reasons` |

autozt 回拉清单（`fetch_files`）：上述 8 个 JSON/PNG。

---

## 5. 物理模型

二维各向异性有限体积法（红黑排序 SOR，纯 numpy）：

    x = 沟道方向，y = 深度（j=0 顶面/沟道，j=Ny-1 底/衬底热沉）
    稳态：div(k(x,y) grad T) + Q_joule = 0
    瞬态：rho*cp dT/dt = div(k grad T) + Q

- 每层各向异性 kx（面内）/ ky（跨面）；层间用有限界面热导 G 连接。
- 底边界恒温 T_AMB；左右金属接触恒温（sink）或绝热（insulator）。
- 焦耳功率密度 Q 只加在沟道层。

一维串联网络（用于对照）：

    R_pp = t_ch/ky_ch + 1/G + t_ox/kx_ox        [m²K/W]
    ΔT_1D = Q * R_pp
    tau_fast ≈ (t_ch/ky_ch + 1/G) * rho_cp_ch * t_ch      (亚 ns)
    tau_RC   ≈ R_th_total * C_total                       (百 ns)

程序同时输出二维解与 1D 网络解，以及 `dT_ratio_FVM_over_1D`（二维修正比例）。

---

## 6. 四步流水线

| 步骤 | 名称 | 干什么 | 产物 | 判据 marker |
|---|---|---|---|---|
| S1 | step1_device_spec | 固化器件几何/热源/边界 | device_spec.json | `DEVICE_SPEC` |
| S2 | step2_props | 上游 κ+厚度 + 材料库兜底 + TBC | thermal_props.json | `PROPS_DONE` |
| S3 | step3_solve | 二维 FVM（稳态 + 瞬态） | device_thermal_summary.json + T_field.json | `SOLVE_DONE` |
| S4 | step4_report | 温度场图 + 汇总 | device_report.json (+ PNG) | `REPORT_DONE` |

全部 `run: gen`（登录节点执行，不提交 SLURM）。依赖链 S1→S2→S3→S4。

---

## 7. 用法

```bash
# 0) 可选：确认上游 kl-dft-cpu 已算出 kappa（KAPPA_SOURCE=auto 时自动取）
autozt -tt kl-dft-cpu -p <材料> status

# 1) 生成 + 推进（S1..S4 自动按依赖跑完）
autozt -tt device-thermal -p <材料> -j S1_spec init
autozt -tt device-thermal -p <材料> start

# 或只生成某一步检查输入
autozt -tt device-thermal -p <材料> -j S3_solve init
```

改参数（不要改脚本）：

```bash
autozt -tt device-thermal -p MoS2 -j S1_spec conf --set params.Q_JOULE=5e8
autozt -tt device-thermal -p MoS2 -j S1_spec conf --set params.CONTACT_BC=insulator
autozt -tt device-thermal -p MoS2 -j S2_props conf --set params.TBC_OVERRIDE=1.2e8

# 查看合并后的最终值（只读）
autozt -tt device-thermal -p MoS2 -j S1_spec conf
```

参数三层合并，后者覆盖前者：**技能 `templates/step.conf` → 项目 `templates/` → 项目 `step.conf`**。

**实践上最值得调的三个**：`Q_JOULE`（热源强度）、`TBC_OVERRIDE`（界面热导）、`CONTACT_BC`（边界上下限）。

---

## 8. 材料库与数据来源

- 默认 MoS2 面内 κ = 20.8 W/mK，取自 AutoZT `kl-dft-cpu` 实测（RTA raw, 300 K）。
- 默认 MoS2/SiO2 TBC = 14 MW/m²K，取自 Yalon et al., Nano Lett. 17, 3429 (2017)。
- **文献值只作量级兜底**；正式投稿请用上游实测 κ，并在文中注明口径。

---

## 9. GPU 可选路线（`optional_steps.gpu_tbc`，默认关）

用来给界面 TBC 提供第一性来源，再把 G 通过 `TBC_OVERRIDE` 灌回 S2：

- **S5_tbc**（GPU）：通用生成 GPUMD NEMD 输入 —— 任意 extxyz 结构/元素、任意传输轴 x/y/z、
  任意 GPUMD 势；界面按元素组分突变自动判定；deck 为「平衡 → 烧入 → 测量」三段。
- **S6_tbc_post**（GPU）：热流 q **直接来自源/漏恒温器的累积传能**（不依赖文献 κ），
  自动挑最长稳态窗口 + 两侧体区外推得 ΔT_i，输出 G 与 `quality` 标记。
- **S6_bte**（未安装）：OpenBTE / JAX-BTE 弹道修正。

```bash
# 回灌 GPU 算出的 G
autozt -tt device-thermal -p <材料> -j S2_props conf --set params.TBC_OVERRIDE=<G>
autozt -tt device-thermal -p <材料> start
```

> 可执行配方见 **GPU-TBC-RECIPE.md**。

---

## 10. 验证证据（2026-09-23）

- **求解器交叉验证**：本技能红黑 SOR 与独立 scipy spsolve 原型同一输入下 dT_peak 均为 19.93406 K，
  相对差 **1.19e-10**；与 1D 解析热阻串联一致到机器精度。
  ⚠️ 该交叉验证**只覆盖求解器，不覆盖源项口径**——两套实现共用同一源项约定，所以它没能发现下面的 bug。
- **网格收敛**（n_sub=15/30/60/120）：19.294 / 19.714 / 19.934 / 20.046 K，单调收敛。
- **源项修正（2026-09-23）**：独立 FEM 复核（技能 `device-thermal-fem`）发现 S3 源项
  `q[:, lay == 0] = Q*dx` 给沟道**每个**单元各加一份，而沟道离散为 2 层，
  实际注入功率达名义 `Q*L` 的 **2.05 倍**，dT 偏高 **2.018 倍**。
  已修正为按单元厚度分摊 `q[:, lay == 0] = Q*dx*dy[lay == 0]/t_channel`。
  **修正后**：FVM dT_peak = **9.8177 K** vs 独立 FEM 9.7307 K（差 0.89%），FEM 侧四类检查全 PASS。
- **端到端**（真实 MoS2 上游 κ=21.43 W/mK、厚度 6.7177 Å）：FVM dT_peak = **9.8177 K**
  （修正前 19.6355 K，即偏高 2 倍）vs 1D 网络 28.5864 K（ratio 0.343）；SOR 4295 次迭代；四步 wall ≈ 18 s。
- **参数覆盖**：`CHANNEL_THICKNESS=1e-9`、`TBC_OVERRIDE=2.5e7` 均生效且 source 正确标注。
- **autozt 集成**：`schema --strict` ✓（io_schema 入/出/参 3/7/49、flow ✓、0 问题）；py_compile OK。
- **GPU 端到端**：新版 S5 生成 deck → 3090 GPU5 跑 320000 步（485 s）→ S6 出结果。
  注意该 Si/Ge 例子的 G **未收敛**（两种子 0.92 / 2.21 GW/m²K，均被判 `quality=poor`）。

完整记录：`tmp/device_thermal_e2e/VALIDATION.md`。

---

## 11. 局限与注意事项

1. **`SUBSTRATE_MATERIAL` 目前只作标签**：底边界恒为 T_AMB，未按衬底材料建热阻。
2. **`W_DEVICE` 不进传热方程**：求解按单位宽度，此参数只用于折算总功率。
3. **层数固定两层**（沟道 + 氧化层）；更多层需改 `device_common.solve_device`。
4. **连续介质假设**：`L_CHANNEL` 接近或小于声子平均自由程时结果偏乐观，需弹道修正。
5. **材料库是文献兜底**，含 2D 材料有效厚度口径差异；正式结果请用上游实测。
6. GPU 步的 G 依赖引线长度与稳态：短引线会被新版的 `TBC_DT_TOL`/`TBC_MARGIN_TOL` 判为 `poor`，
   此时**不要引用具体数字**，先加长引线并多种子平均。

---

## 12. 自检

```bash
python3 bin/autozt skills                          # device-thermal 应出现在表里
python3 bin/autozt schema device-thermal --strict  # io_schema / flow / 0 问题
```
