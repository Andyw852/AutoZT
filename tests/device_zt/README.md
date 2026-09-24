
# device-zt 回归用例与技能说明

- 技能目录：`skill/device-zt/`（文档 `skill/device-zt/README.md`）
- 用例目录：`tests/device_zt/cases/`
- 运行器：`tests/device_zt/run_cases.py`（自包含，不联网、不连超算、不提交 SLURM）

```bash
python3 tests/device_zt/run_cases.py          # 跑全部用例
python3 tests/device_zt/run_cases.py -v Si_3D # 只跑一个，逐项打印
```

运行器做的事（**刻意复刻 autozt 在超算上的运行方式**）：

1. 把 `case/upstream/` 摊成材料目录 `<tmp>/<case>/<技能>/<步骤>/<文件>` —— 这正是
   `find_upstream` 在超算上看到的第一候选布局 `<matdir>/<技能>/<步骤>/<文件>`；
2. 用 autozt 自己的 `report.find_asset` 解析技能 `gen_need` 并平铺进 `<tmp>/<case>/device-zt/`
   （顺带验证 `skill/_common/<组>/` 公共池解析）；
3. 用 `case.json` 的 `params` 覆盖技能模板 `step.conf`；
4. **清空 `PYTHONPATH`**，在 `<tmp>/<case>/device-zt/` 原地顺序跑三个 gen 脚本；
5. 读产物，与 `case.json` 的 `expected` 逐项比对（数值相对容差 1e-6）。

---

## 1. 用例一览

| 用例 | 维度 | 材料 | 关键参数 | dT_peak | 器件 ZT | V_max_safe |
|---|---|---|---|---|---|---|
| `MoS2_2D` | 2D | MoS2 | 默认（1e20，R_contact=0） | 21.3771 K | 0.0242778 | 0.20 V |
| `Si_3D` | 3D | Si | + CHANNEL_THICKNESS/SIGMA_T_CONV_M=100 nm | 2.6734 K | 0.0091542 | 0.50 V |

两个用例都另外通过了**真实 autozt CLI** 端到端（`init` → `-j S1_elec init` → `start` → `start`，
三步全 `OK` 并 fetch 回本地），数值与这里的期望值**逐项完全一致**。
完整记录见 `tmp/device_zt_e2e/VALIDATION.md` §8。

### 用例目录结构

```
cases/<用例名>/
├── case.json                  # 维度/材料/参数覆盖/期望标量/说明
├── upstream/                  # 技能的真实前置产物（原样保存）
│   ├── ke-dft-cpu/<步骤>/...
│   └── kl-dft-cpu/<步骤>/...
└── expected/                  # 三步的真实产物（原样保存）
    ├── step1_elec/elec_props.json
    ├── step2_selftrap/{device_zt_solution.json, T_field.json}
    └── step3_report/device_zt_report.json
```

---

## 2. 技能输入

### 2.1 前置产物（必须已存在于**同一台超算的同一材料目录**）

| 文件 | 位置 | 提供什么 |
|---|---|---|
| `transport.json` | `<mat>/ke-dft-cpu/step8_amset/`（也认 `step8.4_amset2d`、`step8.1_boltztrap`） | σ(T)、S(T)、κe(T) 三张温度表 × 10 档 signed 掺杂 |
| `2d_correction.json` | 同上（**2D 才有**） | `cell_c_A` = 含真空超胞高 h⊥ → 定面电导口径 |
| `kappa_summary.json` | `<mat>/kl-dft-cpu/step6_kappa/` | κ(T)；2D 优先读 `kappa_2d_normalized_xx_yy_zz`（物理 2D 值） |
| `thickness_2d.json` | 同上（2D 才有） | 物理层厚 d（单层厚度） |
| `device_materials.json` | 技能自带（公共池） | 材料库：κ_inplane/κ_cross、层厚、TBC 及其文献出处 |

另外材料目录里要有 `POSCAR`：**autozt 用它识别材料目录**（`init` 没它会报"找不到材料"）。
device-zt 三个 gen 脚本**都不读结构**（已 grep 确认），所以对 POSCAR 内容无要求。

### 2.2 ★ 口径约定（本技能最容易错的地方，也是核心设计）

- **σ 与 κe 同为「含真空超胞」口径**（AMSET / BoltzTraP 的输出就是如此）。
  故 **面电导 `G_s = σ_cell × h⊥`**，等价于 `σ_2D × d`，**与口径无关** → 器件电阻、
  焦耳热、ZT 全部不受影响。
- **κ_L 是物理 2D 值**（kl 已用 h⊥/d 归一化过），**不能**直接与 κe_cell 相加。
  必须先换算 `κe_2D = κe_cell × t_conv/d`（输出里保留 `ke_cell_to_layer_factor` 便于核对）。
- **3D 特例**：体材料 σ/κe 本就是 per-volume，取 `SIGMA_T_CONV_M = 沟道厚` 即可，
  此时因子 = 1，公式自然退化。
- **`step.conf` 不剥引号**：`KE_UPSTREAM_STEP = ""` 解析出来是**两个引号字符**的真值字符串。
  `zt_common.blank()` 统一把 `""`/`''`/`none`/`null`/`-` 当空。

### 2.3 可调参数（`step.conf`，项目级可覆盖）

器件几何：`DEVICE_NAME`、`CHANNEL_MATERIAL`、`CHANNEL_THICKNESS`、`OXIDE_MATERIAL`、
`OXIDE_THICKNESS`、`SUBSTRATE_MATERIAL`、`L_CHANNEL`、`W_DEVICE`、`T_AMB`、`CONTACT_BC`
（`sink` 接触恒温 / `insulator` 接触绝热）。
热学：`TBC_OVERRIDE`、`KAPPA_SOURCE`、`UPSTREAM_SKILL`、`UPSTREAM_STEP`。
电学：`KE_UPSTREAM_SKILL`、`KE_UPSTREAM_STEP`、`DOPING_CM3`、`T_REF`、`INPLANE_MODE`、
`SIGMA_T_CONV_M`。
驱动与自洽：`V_BIAS`、`V_SWEEP_MIN/MAX/N`、`SC_MAXIT`、`SC_TOL_K`、`SC_RELAX`。
可靠性：`R_CONTACT_OHM`、`DT_LIMIT_K`、`T_LIMIT_K`。
数值：`N_GRID_X`、`N_GRID_SUB`、`SOR_OMEGA`、`SOR_TOL`、`SOR_MAXIT`。

**3D 必须显式给三项**（`CHANNEL_THICKNESS`、`SIGMA_T_CONV_M`，以及 `CHANNEL_MATERIAL` 换成 3D 材料）——
材料库里 `Si` 的 `thickness_m = None`，不给会直接报"沟道厚度未知"。

---

## 3. 技能输出

三步都是 `run: gen`：登录节点、纯 numpy、秒级、**不提交 SLURM**。

### 3.1 S1_elec → `step1_elec/elec_props.json`

读上游、定口径、算面电导，并输出 10 档掺杂的电学表。

| 字段 | 含义 |
|---|---|
| `t_conv_m` | 面电导口径厚度（2D = `cell_c_A`；3D = `SIGMA_T_CONV_M`） |
| `sigma_convention` / `sigma_convention_source` | 口径名与来源（`2d_correction.json` 还是 `step.conf`） |
| `doping_cm3` / `doping_all_cm3` | 选中档与实际可选档（signed，负=p 型） |
| `sigma_S_m[]` / `seebeck_uV_K[]` / `kappa_e_W_mK[]` | 10 档 × 300 K 的 σ/S/κe（κe 为**胞口径**） |
| `G_sheet_S[]` / `R_sheet_ohm_sq[]` | 逐档面电导（S/square）与面电阻 |
| `upstream_skill_kl` / `upstream_step_kl` | 实际用的 kl 上游 |

### 3.2 S2_selftrap → `step2_selftrap/device_zt_solution.json` + `T_field.json`

焦耳自热闭环：解 αT 温度场 → 层平均温度回代查 σ(T)/S(T)/κe(T) → 迭代到收敛。

| 分组 | 字段 |
|---|---|
| 温度 | `dT_peak_K`、`T_eff_final_K`、`T_channel_mean_K`、`T_channel_max_K` |
| 电 | `I_A`、`R_channel_ohm`、`R_contact_ohm`、`R_electrical_ohm`、`P_W`、`P_channel_W`、`P_contact_W`、`contact_share`、`sigma_S_m`、`seebeck_V_K` |
| 热 | `Q_joule_W_m2`、`kappa_lattice_W_mK`、`kappa_e_cell_W_mK`、`kappa_e_W_mK`、`ke_cell_to_layer_factor`、`kappa_total_W_mK`、`tbc_W_m2K`、`tbc_source`、`sor_iterations` |
| 自洽 | `converged`、`iterations`、`history[]`（每轮的 T_eff/σ/Q/dT/κe） |
| 器件 ZT | `zt_material`（= S²σT/κ_total @ T_eff） |
| 可靠性上限 | `V_max_safe_V`、`P_max_W`、`V_at_P_max_V`、`dt_limit_K`、`T_limit_K`、`within_limit` |
| I–V 扫描 | `iv_sweep[]`：每点 `V_V/I_A/P_W/Q_joule_W_m2/dT_peak_K/T_eff_K/T_channel_max_K/sigma_S_m/zt_material/converged/within_limit` |

`T_field.json` 是 x–y 平面温度场网格（沟道 2 格 + 氧化层 `N_GRID_SUB` 格），供画图/自查用。

### 3.3 S3_report → `step3_report/device_zt_report.json`

汇总成四段：`geometry` / `electrical` / `thermal` / `zt`，并给出
`zt.doping_scan.rows[]` —— 10 档掺杂在 T=300 K、当前器件 κ_L 下的 ZT 估算，
用于**找最优掺杂**（不随器件口径变化）。

---

## 4. 结果怎么分析

### 4.1 八个独立校核（都可手工复算，不依赖代码路径）

```
G_sheet   = σ(T_eff) × t_conv
R_channel = L / (G_sheet × W)
I         = V / (R_channel + R_contact)
Q         = P_channel / (L × W)           P_channel = I² R_channel
dT_peak   = Q × R_th_eff
ZT        = S(T_eff)² × σ(T_eff) × T_eff / κ_total
κe_2D     = κe_cell × t_conv / d
κ_total   = κ_L + κe_2D
```

**Wiedemann–Franz 必须同口径**：`L = κe_cell / (σ_cell × T)`。
若误用 κe_2D 除 σ_cell，比值会虚高 h⊥/d 倍（MoS2 是 5.19 倍），这是**口径混用**不是物理问题。

### 4.2 物理判读要点

- **`L/L0` 随掺杂升高趋向 1**：非简并半导体 L < L0（0.7 左右），简并极限 → L0 = 2.443e-8。
  MoS2 实测：1e17 0.746 → 1e20 0.724 → **1e21 1.012**（进入简并）。趋势对才是对的。
- **κ 电子占比**：MoS2 1e20 约 15%，Si 约 1.4%（Si 声子主导），量级必须合理。
- **电流随偏压次线性**：自热把 σ 压下去（负反馈）。若 I–V 是直的，说明自洽没生效。
- **ZT 随偏压上升 ≠ 该工作点可用**：MoS2 在 0.5 V 时 ZT 升到 0.048，但 T_max = 609 K。
  **`P_max` 与 `V_max_safe` 必须分开看** —— 只看 P_max/ZT 会选出烧器件的工作点。
- **`within_limit`** = `dT ≤ DT_LIMIT_K` 且 `T_max ≤ T_LIMIT_K`；`V_max_safe_V` 是扫描内最大的合格偏压。

### 4.3 两个用例的分析

> **符号约定**：`transport.json` 的 doping **负 = n 型电子、正 = p 型空穴**
> （见 `ke-dft-cpu/step8_amset/gen_step10_amset.py` 的注释与 `step8.1_boltztrap` 的
> `branch = "electron" if n < 0 else "hole"`）。判据：负掺杂的 S 为负、且 Si 的
> 负掺杂分支迁移率更高（207 > 115 cm²/Vs），恰好符合电子比空穴易动的物理。

**MoS2 2D（1e20，V=0.1 V）**：I=7.220e-4 A，R=138.50 Ω，Q=3.610e8 W/m²，
dT=21.377 K，T_eff=314.61 K；κ_L=63.700 / κe_2D=11.660 / κ_tot=75.360（电子占 15.5%）；
**器件 ZT=0.02428**，9 次迭代收敛。上限：V_max_safe=0.20 V，P_max 出现在 0.5 V（不可用）。
掺杂扫描（**负=n 型电子、正=p 型空穴**）最优是 **p 型 1e20（ZT 0.02226）**；n 型峰在 1e19
（ZT 0.00946）。重掺 p 型 1e21 掉到 2.77e-5 —— 教科书级掺杂优化行为。
> 注意：ke 与 kl 来自**两次不同 DFT 超胞**（ke h⊥=34.879 Å，kl h⊥=19.98 Å）。各自归一化成物理
> 2D 值，`G_s = σ·h⊥` 对胞高不敏感，所以结果自洽；但这不是真实材料的自然配对。

**Si 3D（1e20，V=0.1 V）**：I=9.190e-3 A，R=10.881 Ω，Q=4.595e9 W/m²，
dT=2.673 K，T_eff=301.77 K；κ_L=88.809（取上游 `kappa_xx_yy_zz`，3D 不做 2D 归一化）/
κe=1.301 / κ_tot=90.110（电子占 **1.44%**）；**器件 ZT=0.009154**。
上限：**V_max_safe=0.50 V**（dT 才 65 K，Si 高 κ 使器件更凉）。
掺杂扫描最优是 **n 型 1e20（ZT 0.0407）**，p 型峰在 1e21（ZT 0.0121）—— 方向与"Si 电子比空穴
易动"一致（1e20 迁移率 207 vs 115 cm²/Vs）。**但这个数偏高**：见 §4.4。
> 注意：Si/SiO2 的 TBC 走的是材料库**通用 vdW 量级默认值 1e7 W/m²K**，不是 Si/SiO2 专用实测值
> （真实界面通常 1e8~1e9），所以 `dT` 绝对值只能作量级参考；要严肃用请设 `TBC_OVERRIDE`。

---

### 4.4 ⚠ 上游数据的两处保留意见（本次核对发现，都不是 device-zt 的问题）

**① Si 的 κ_L 明显欠收敛**：上游 `kappa_summary.json` 来自 `kappa-m151515.hdf5`，
即 **15×15×15 q 网格**，给出 300 K 的 κ = **88.81 W/mK**。体硅公认值 **~150 W/mK**
（材料库 `device_materials.json` 里 Si 写的正是 `kx = 150.0`）——15³ 网格对 Si 远未收敛。
而 `KAPPA_SOURCE = auto` **无条件优先上游**，于是用 88.81 而不是 150：
κ 低 41%，**ZT ∝ 1/κ 被放大 ~1.7 倍**，dT 同向偏大。
要可靠值：项目 `step.conf` 里设 `KAPPA_SOURCE = literature`（用库值 150.0），
或把上游 phono3py 的 q 网格加密后重算。

**② 空穴分支的 μ_ADP 比文献高约一个量级 —— 已定位到「不是重叠模式、不是分维问题」**

`transport.json` 的 `mobility.ADP` 在 **p 型（正掺杂）** 分支比 n 型分支高 **3~4 倍**：

| 材料 | n 型 μ_ADP（1e20） | p 型 μ_ADP（1e20） | 文献对照 |
|---|---|---|---|
| Si | 1244.6 | 4163.0 | Si 300 K 声子限：电子 ≈ 1400、空穴 ≈ 470 |
| MoS2 | 1361.8 | 5289.2 | MoS2 通常电子比空穴易动 |

n 型分支与文献吻合、p 型分支高约一个量级。**本技能只是读 σ**，且已核实 AMSET 的 σ
与其自带 mobility 完全一致（比值 1.0000），所以源头在上游 `ke-dft-cpu`。

**追查结论（2026-09-24 本地受控实测，未动集群）**：这**不是** 2D/3D 分类问题。
上游的分维处理是**显式且已验证**的 —— 2D 强制 `unity_overlap: true`、3D 显式 `false`
（真实重叠），且 3D 的 Si 全网格对照已完成（去对称化对 ADP 只有 1.01/0.92）。
我补做了**上游从未做过的 `unity_overlap` 单变量臂** —— 同一 vasprun / settings /
形变势，只多一行 `unity_overlap: true`：

| 机制（300 K，1e17） | 真实重叠（生产） | unity | unity / 真实 |
|---|---|---|---|
| 电子 ADP | 1608.4 | 363.4 | **0.23** |
| 电子 overall | 805.0 | 247.1 | **0.31** |
| 空穴 ADP | 5434.8 | 1875.8 | **0.35** |
| 空穴 overall | 1474.3 | 681.1 | **0.46** |

（参考臂逐位复现了上游 VERIFICATION §六 ③ 的记录 1608.4 / 5434.8 / 805.0，台子可靠。）

**结论：改 unity 只会把电子打到文献值的 1/4（363 vs ≈1400），空穴仍高 4 倍 ——
unity 不是解药，上游「3D 不改 unity」的决定被实测确认正确。**

**空穴分支偏高的真正原因（2026-09-24 定论）：AMSET 没有非极性光学声子（光学形变势 ODP）
散射。** AMSET 实现的机制只有 **ADP / POP / IMP / PIE**（源码 `amset/scattering/elastic.py` =
ADP+IMP+PIE，`inelastic.py` 只有 `PolarOpticalScattering`；官方文档 Summary of scattering
rates 也只列这四种）。Si 是非极性晶体，POP 前因子为零（`scattering/calculate.py:842` 注释
明说），于是**整条光学声子通道缺失** —— 而 Si **空穴**的声子限迁移率里光学声子贡献很大。

电子/空穴的双重实验锚点正好对上：

| 分支 | AMSET（1e17） | 本征声子限实验值 | 倍数 |
|---|---|---|---|
| 电子 ADP | 1622.7 | ≈1400 | 1.16 |
| **电子 overall** | **805.0** | **≈800** | **1.01** |
| 空穴 ADP | 5434.8 | ≈470 | **11.6** |
| **空穴 overall** | **1474.3** | **≈330** | **4.5** |

**电子几乎完全吻合**（1e17 实验值 ≈800 ✔），**空穴高 4.5 倍** —— 与「缺一项主要散射机制」
的预期一致。**所以不是形变势算错**：`deformation.h5` 里价带顶 D(Γ)=1.74 eV，接近文献
价带静水压形变势 a_v ≈ 2 eV，是合理值。**「价带形变势绝对归一化」不作为立项方向**
（2026-09-24 已结案）：AMSET 官方 Si 算例（matsci.org 论坛）空穴 ≈5400、电子 ≈1300，
与本流程 ADP 分支（空穴 5435、电子 1608）几乎一致 —— 空穴偏高是 AMSET 模型本身的特性
（无非极性光学声子散射），流程设置无误。证据：
https://matsci.org/t/error-in-the-hole-mobility-of-silicon-in-the-example-of-asme/59809
已写进上游 `step8.4_amset2d/README.md` §已知局限、`METHODOLOGY.md` §5.1 与
`VERIFICATION.md` §六点一。

**后果（收敛后的说法）**：`device-zt` 用的 σ 是 **overall** 迁移率。
- 在 **1e19~1e20**（热电用的掺杂档）IMP 主导，ADP 分支偏差对器件 ZT 的影响
  **没有分机制数字看起来那么大**；
- 在 **1e17** 附近 Si 仍是晶格散射主导，这条缺陷会实打实抬高 p 型迁移率；
- 因此 **Si「p 型不如 n 型」的方向性结论是在不利偏差下仍然成立的** ——
  缺陷把 p 型迁移率**高估**了（真实 p 型更低、差距更明显），反而在削弱这个结论；
  **结论在这种偏向下仍成立，所以更可信**；
- **MoS2 是极性材料，有 POP，不受这条影响。**
- 分机制 ADP 值在引用时要附上表；上游分维/重叠的结论（3D 不改 unity）经本轮实测确认不变。

**③ 因此 Si 用例的 ZT = 0.0407 只应作为管线自洽性证据，不能引用为体硅的器件 ZT。**
换成库值 κ = 150 重算约 **0.0246**；而文献 300 K 体硅重掺 ZT 约 **0.005~0.01**
（本数据集 S 与 μ 都偏大，所以仍偏高约 2~4 倍）。

## 5. 局限（不要过度解读）

1. **只做焦耳自热**：无 Peltier 界面热、无 Seebeck 反向电动势、无 Thomson 效应。
2. **均匀掺杂 + 层平均温度**：用沟道层平均温度查表，不解析沟道内温度/掺杂分布。
3. **接触电阻是简单串联**（`R_CONTACT_OHM`）：不做真实接触几何/电流拥挤（TLM）。
   设计选择：**只有沟道体内耗散进温度场**，接触耗散留在接触处（`CONTACT_BC=sink` 时接触即热沉）；
   两项在输出里单列（`P_channel_W`/`P_contact_W`/`contact_share`）可自行核对。
4. **`DOPING_CM3` 是 3D 密度**：面密度 n_2D = n_3D × h⊥。
5. **单层材料、两层结构**（沟道 + 氧化层），无栅压、无量子限域修正。
6. **上游必须同一台超算同一材料目录**：三个 gen 步骤在超算上执行（autozt 无本机执行模式），
   而 `find_upstream` 以 `dirname(cwd)` 作材料目录。跨集群或只在本地 `result/` 有回拉结果时，
   S1 会报 `[ERROR] 找不到 ke 的 transport.json`。**不回退本地 `result/`**。
7. **不提交 SLURM**：全流程是登录节点上的 numpy 后处理。
