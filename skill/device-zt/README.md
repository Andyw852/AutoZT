# device-zt —— 器件级电-热自洽 ZT 仿真

把 **ke-dft-cpu 的电子输运**（σ / Seebeck / κe）接进 **device-thermal 的传热 FVM**，
让「电流 → 焦耳热 → 温升 → 物性随 T 变 → 焦耳热变」自己迭代到**自洽**，
输出器件工作点（I、P、ΔT、T_eff）、**I-V / P-V 曲线**，以及**全掺杂档的 ZT 扫描**。

**定位**：device-thermal 只解「给定一个 `Q`，温度场是多少」；本技能解「给定一个偏压，**电流自己产生**的 `Q` 对应什么温度，而温度又反过来改电流」。

---

## 1. 三步（都是 `run: gen`，登录节点秒级，纯 numpy，不提交 SLURM）

| 步 | 名字 | 输入 | 产物 |
|---|---|---|---|
| S1_elec | 电学输入 | ke 的 `transport.json` + `2d_correction.json` | `elec_props.json` |
| S2_selftrap | ★ 自洽迭代 | S1 + 上游 kl 的 κ + 材料库 | `device_zt_solution.json` + `T_field.json` |
| S3_report | 汇总 | S2 | `device_zt_report.json` |

```bash
autozt -tt device-zt -p <材料> -j S1_elec init      # 只生成输入，不提交
autozt -tt device-zt -p <材料> start                # 依次跑 S1 → S2 → S3
autozt -tt device-zt -p <材料> -j S2_selftrap conf  # 看合并后的 step.conf
```

---

## 2. ★★ 口径（本技能最容易错的地方）

上游给的量**不是同一套口径**，混用就错：

| 量 | 上游口径 | 器件里要用的形式 | 换算 |
|---|---|---|---|
| 晶格 κ_L | **物理 2D 值**（`kappa_2d_normalized_*`） | 面内热导 | 无（`read_kappa` 已归一化） |
| 跨面 κ_⊥ | 物理值，但 2D 产物常是**数值零** | 跨面热阻 | κ_⊥ 无效时回落材料库 |
| **σ** | **含真空超胞值** | **面电导** `G_s = σ_cell × h⊥` | 见 2.1 |
| **κe** | **含真空超胞值** | 面内热导 | `κe_2D = κe_cell × h⊥/d`（见 2.2） |
| Seebeck | 与口径无关 | 直接用 | 无 |

### 2.1 为什么 σ 用面电导就"与口径无关"

```
σ_cell × h⊥  ==  σ_2D × d          （两者都等于面电导 G_s，单位 S / square）
```

所以电阻、焦耳热、ZT 都不依赖上游用的是胞体积还是层厚。本技能因此**只算 G_s**。

`h⊥` 取自 `transport.json` 同目录的 `2d_correction.json` 的 `cell_c_A`。
AMSET 与 BoltzTraP 两条路径的 σ **都是超胞口径**——BoltzTraP 步自己的注释写着
「σ、κ_e 仍按超胞体积口径（与 amset 一致）」。
另一条路：`autozt ... conf --set params.SIGMA_T_CONV_M=<m>` 显式指定。

> ★ 顺带的好处：`G_s = σ_cell × h⊥` 对**胞高不敏感**。ke 与 kl 两次 DFT 用的超胞
> 高度常常不同，用「σ_cell × 各自的 c」比用「σ_cell × c/t × 别人的 d」更稳。

### 2.2 κe 必须换算（**曾经的混口径 bug**）

σ 和 κe 是同一个超胞口径（实测 MoS2 掺杂 1e21 cm⁻³：κe_cell = 4.1 W/mK，
Wiedemann-Franz 数 L = κe/(σT) = 2.47e-8 ≈ L₀ = 2.44e-8，物理自洽）。
但它跟**物理 2D 值**的晶格 κ_L（63.7 W/mK）不是一套口径。若直接 κ_L + κe_cell，
等于把 4.1 加到 63.7 上，电子份额被低估 5.19 倍。

正确：`κe_2D = κe_cell × h⊥/d = 4.1 × 5.19 = 21.3 W/mK`，
面上总 κ = 63.7 + 21.3 = 85.0 W/mK，**电子占 26%**（忽略 κe 会让器件温度高估约 1/4）。

### 2.3 step.conf 里 `""` 的坑

`stepconf` **不剥引号**：写 `KE_UPSTREAM_STEP = ""` 解析出来是**两个字面量引号**
的字符串（真值），会被当成"指定了一个叫空串的步骤"→ 找不到上游。
本技能用 `zt_common.blank()` 统一把 `""` / `''` / `none` / `null` / `-` 当空处理。

---

## 3. 自洽算法（S2）

```
T_eff <- T_amb
重复：
    σ, S, κe <- 按 T_eff 从 transport 表线性插值（超出网格不外推，夹到端点）
    G_s = σ × h⊥
    R   = L / (G_s × W)                # 沟道体电阻（不含接触电阻）
    I   = V / R
    Q   = V² × G_s / L²                # 焦耳面功率密度 W/m²（按俯视面积分摊）
    κ_channel = κ_L + κe_2D            # 电子热导并入（同为物理口径）
    跑 device_common.solve_device(Q, …)  ->  T(x,y)
    T_eff' <- mean(T[:, 沟道层])        # 集总模型的温度代表值
    T_eff  <- T_eff + relax × (T_eff' - T_eff)
直到 |ΔT_eff| < SC_TOL_K 或达到 SC_MAXIT
```

S2 还对 `V_SWEEP_MIN…V_SWEEP_MAX` 的每个偏压**各自**做一遍上述自洽，
得到真实的 I-V / P-V 曲线（因为 R 随 T 变，曲线不是直线）。

---

## 4. ZT 与器件指标

- `zt_material = S² σ T_eff / (κ_L + κe_2D)` —— 材料品质因子，**与电导口径无关**。
- S3 还给**全掺杂档 ZT 扫描**（`transport.json` 里本来就有 10 档）：ZT 对掺杂极敏感，
  重掺下 S 塌到几 µV/K、ZT 归零，必须看折中在哪一档。
  扫描里 κ_L 固定为当前器件值（近似：κ_L 与掺杂关系弱），ZT 是比值故与口径无关。
- 器件侧的真实差异体现在 `dT_peak_K`、`R_electrical_ohm`、`P_max_W` 上：
  自热升温把 T_eff 抬高 → σ 变 → 电流变 → 焦耳热变，这是 device-thermal 给不出的。

### 4.1 接触电阻（简单串联）

设 `R_CONTACT_OHM > 0` 后：R_total = R_channel + R_contact，I = V/R_total。

```
R_channel = (L/W) * R_sheet = L/(G_s*W)
I         = V / (R_channel + R_contact)
Q_channel = I^2 * R_channel / (L*W)      <- 只有这部分进温度场
P_contact = I^2 * R_contact              <- 留在接触处（CONTACT_BC=sink 时接触即热沉）
```

为什么不把接触耗散也塞进沟道：本技能的接触边界默认就是**热沉**，接触处产生的热
沿最近的路径进热沉，加热沟道不是主导路径。把接触热算进沟道会**高估** dT。
两项在输出里单列（`P_channel_W` / `P_contact_W` / `contact_share`），便于核对。

### 4.2 温度/可靠性上限判据

S2 的每个扫描点都带 `within_limit` = (dT <= DT_LIMIT_K) and (T_channel_max <= T_LIMIT_K)`

- 默认 `DT_LIMIT_K = 100 K`、`T_LIMIT_K = 450 K`。
- S2/S3 报出 `V_max_safe_V`（扫描范围内仍满足判据的最大偏压）。
- **超限不等于算错**：只说明这个工作点超出热/可靠性预算，模型本身照常有效。
- 实测意义很大：MoS2 基线下 `P_max` 出现在 0.5 V，但 `V_max_safe` 只有 0.2 V——
  **P_max 那个工作点会把器件烧掉**，真正可用的是 0.2 V。

---

## 5. 模型局限（**明确声明，不要过度解读**）

1. **只做焦耳自热**：不含 Peltier 界面热、Seebeck 反向电动势、Thomson 效应。
2. **均匀掺杂 + 集总温度**：用沟道层**平均温度**查 σ(T)/S(T)/κe(T)。真实器件里
   沟道温度不均匀、掺杂有空间分布，这里都没有。
3. **接触电阻是简单串联近似**（`R_CONTACT_OHM`）：不做真实接触几何/电流拥挤（TLM），
   故接触电阻对温度场的影响是**保守**的（接触热不进沟道）。真实器件里接触电阻常占主导。
4. **掺杂是 3D 密度**：`DOPING_CM3` 是 cm⁻³，面密度 n_2D = n_3D × h⊥。默认已改为 **1e20**
   （MoS2 的 ZT 最优档，S≈114 µV/K）；**1e21 是重掺极限**，n_2D ≈ 3.5e14 cm⁻² 远超典型 FET，
   且 S 塌到 3.8 µV/K、ZT 归零——不要用它做 ZT 演示。
5. **单层材料、两层结构**（沟道 + 氧化层），与 device-thermal 同。
6. **σ/S/κe 未含温度以外的影响**（无栅压、无量子限域修正）。

---

## 6. 复用（不是复制）

传热 FVM（红黑 SOR）、材料库、上游定位、`read_kappa` 的 2D 口径判定、
`resolve_thermal_props` 全部来自**公共池** `skill/_common/thermal/device_common.py`。
技能目录里没有第二份实现——`find_asset` 的解析链是
`<技能>/templates/<步> → <技能>/templates → <技能>/ → _common/templates → **_common/<组>/**`，
所以 device-thermal 与 device-zt 引用的是**同一个文件**。

---

## 7. MoS2 基线（真实 ke + kl 数据，详见 tmp/device_zt_e2e/VALIDATION.md）

电学口径：cell（h⊥ = 34.879 Å），t_conv = 3.4879e-9 m；器件 L=200 nm、W=1 µm、SiO2 300 nm、TBC=1.4e7 W/m²K。

| 量 | **1e20（默认）** | 1e21（旧默认/重掺） |
|---|---|---|
| σ / S / κe_cell（300K） | 4.296e5 S/m / 114.2 µV/K / 2.28 W/mK | 5.529e5 / 3.769 / 4.10 |
| G_sheet / R_sheet | 1.498e-3 S / 667 Ω/sq | 1.925e-3 S / 519 Ω/sq |
| κ_L / κe_2D / 总（T_eff） | 63.7 / 11.7 / 75.4（电子 15.5%） | 63.7 / 22.4 / 86.1（电子 26.0%） |
| V=0.1 V：I / R | 7.22e-4 A / 138.5 Ω | 9.625e-4 A / 103.9 Ω |
| Q / dT_peak / T_eff | 3.61e8 W/m² / 21.38 K / 314.6 K | 4.813e8 / 26.38 K / 318.0 K |
| **ZT** | **0.0243** | 2.77e-5 |
| 自洽 | 9 次收敛，SOR 2292 it | 10 次收敛，SOR 2190 it |
| P_max（扫到 0.5 V） | 1.023e-3 W @ 0.5 V | 2.187e-3 W @ 0.5 V |
| **V_max_safe**（dT≤100K） | **0.2 V** | 0.2 V |

**两个关键结论**：

1. **默认掺杂从 1e21 改到 1e20 后 ZT 从 2.8e-5 提到 0.0243（547 倍）**——重掺把 S 压到
   几 µV/K，ZT 必然归零；1e20 才是 S 与 σ 的折中。掺杂扫描（10 档全报）显示
   n 型峰在 1e20、p 型峰在 1e19，正是 MoS2 众所周知的电子/空穴不对称。
2. **P_max 与 V_max_safe 不是同一个点**：P_max 在 0.5 V，但那里 dT 会远超 100 K；
   真正可用的是 **0.2 V**。只看 P_max 会给出会把器件烧掉的设计。

**接触电阻实测**（`R_CONTACT_OHM=100 Ω`，同工况）：

| 量 | 值 | 校核 |
|---|---|---|
| R_channel / R_contact / R_total | 135.1 / 100 / 235.1 Ω | I=0.1/235.1=4.253e-4 A ✓ |
| P_channel / P_contact | 2.444e-5 / 1.809e-5 W | I²R，占比 **42.5%** ✓ |
| Q_channel | 1.222e8 W/m² | =P_ch/(L·W) ✓；与无接触的 3.61e8 之比 0.3385 |
| dT_peak | 7.23 K | 比无接触的 21.38 K 低——电流小了，且接触热进热沉 |
| V_max_safe | 0.3 V（↑） | 更凉 → 可承受更高偏压 ✓ |

> 一致性校核：I 比 = 0.589，Q 比 = 0.589² × (135.1/138.5) = 0.3385 ✓ 与实测完全吻合。
