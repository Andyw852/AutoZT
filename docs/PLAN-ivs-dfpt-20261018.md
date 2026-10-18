# 方案：用 QE DFPT 自动计算谷间散射（IVS）参数，接入 AMSET（2026-10-18）

> 状态：**方案稿，未动代码**。用户 2026-10-18 要求出总体方案、交给其他执行者。
> 按 AGENTS.md 铁律 10，给技能新增步骤/脚本属受控操作：执行方**动代码前**先把本方案里"待用户确认"的几项问给用户，得到同意再做。
> 背景与文献依据见 `skill/ke-dft-cpu/METHODOLOGY.md` [0058]。

---

## 1. 目标与边界

**目标**：多谷载流子的零阶谷间形变势 D0（eV/Å）和声子能量 ħω 由第一性原理自动算出，S8.4（2D）/ S8（3D）的 AMSET
自动接入为 IVS 通道，不再依赖手填文献值。

**做什么**
- 在"连接两个谷"的声子波矢 q 上做 QE DFPT，直接计算电声矩阵元 g，换算为 D0。这是 Filipovic, arXiv:2607.18895 (2026)
  直接 DFPT 那一层的做法（他在 34 个指定 q 上算，其中谷间部分是本方案的范围）。
- 覆盖 K↔K′、K↔Q、Q↔Q′、Q↔K（2D TMD 电子），以及 3D 的同类谷间过程（如 Si 的 X 谷 f/g 过程）。
- 输出一个 JSON，S8.4/S8 读取后写成 dp_optical 的 `X>X'`、`X>Y` 通道。

**不做什么**
- 不做 EPW/Wannier 全布里渊区插值，不改 AMSET 的 ADP/POP/IMP。
- 不做 SOC（dp_optical 与现有 ivs_frozen 都只支持非 SOC）；TMD 空穴不在范围内。
- 不做一阶谷间耦合（∝|q−K|；Filipovic Table 5：缺它 μ 偏高约 8%）。
- 谷内 ODP 仍以冻声子（S7.2/S8.21）为主；本方案可顺带给出谷内光学耦合作对照，不替换。

**预期效果（需验证）**：MoS₂ 电子通用流程（不手填参数）的声子限制 μ(300 K) 从约 377 回落到 DFPT 基准附近
（Filipovic 187–212），高温 ADP-only 异常消失。剩余误差来源：K 谷 ADP 用真空 E1（8.54 对 DFPT 6.5）、缺一阶 K↔K′、
常数 D0、非 SOC。预计与基准差 20–30%。

---

## 2. 现有基础（复用，不重写）

| 已有 | 位置 | 用途 |
|---|---|---|
| 常数形变势散射率（ODP/IVS），2D 与 3D | `ke-dft-cpu/step8.4_amset2d/dp_optical.py`、`step8_amset/amset_dpo3d.py` | 率公式 = Jin 2014 Eq.(4)，不改 |
| 冻声子 K↔K′ 链 | ke-dft-cpu 可选组 `ivs_frozen`（S7.3/S7.4/S8.22），`ivs_frozen.py`，产物 `step8.22_ivs_read/ivs_params.json`（schema `ivs_params/1`） | 接线方式照抄；K↔K′ 结果作交叉校验 |
| S8.4 读取冻声子 IVS | `gen_step14_amset2d.py` 的 `IVS_FROZEN`（`_load_ivs`、`inter_channels`、`invalidate_ivs_consumer`） | 新开关照此实现 |
| 2D 谷表 | S8.2 `dpt_result.json` 的 `multi_valley.valleys`（标签、成员 k 中心、`k_grid`、`dE_eV`） | 决定要算哪些 q |
| 3D 谷表 | step.conf `DP_OPTICAL_VALLEYS_E/H`（三个分量） | 3D 手填 |
| QE 环境 | `eph-qe-cpu`（S0_env 探针；jzzn 上 QE 7.3 + Perturbo 3.0：`~/software/AutoZT/qe-7.3-perturbo/`） | ph.x 已装 |

---

## 3. 总体架构

```
ke-dft-cpu                                   eph-qe-cpu（新可选组 ivs_dfpt）
 S1 CONTCAR ──────────────┐
 S8.2 谷表(2D) / 手填(3D) ─┼──> Q1 ivs_qlist   (gen)  定 k 初态、q 列表、对称约化
                           │    Q2 ivs_scf     (QE)   pw.x SCF（用 VASP 结构，不重新弛豫）
                           │    Q3 ivs_ph      (QE，按 q 扇出) ph.x 单 q DFPT + 指定 k、k+q 的矩阵元
                           │    Q4 ivs_read    (gen)  g → D0、ħω，按星归并 → ivs_dfpt_params.json
 S8.4 / S8  <──────────────┘     （IVS_DFPT=auto 读取，优先级高于冻声子 K↔K′ 和手填）
```

**待用户确认 A：放在哪个技能。**
- 推荐：eph-qe-cpu 新增可选组 `ivs_dfpt`（QE 环境、赝势、模板都在这里），ke-dft-cpu S8.4/S8 读其产物。
  需要定跨技能取文件的方式：step.conf `IVS_DFPT_SRC` 写路径 + 记录 sha256，或 Q4 把结果复制进 ke-dft-cpu 材料目录。
- 备选：放进 ke-dft-cpu（与 ivs_frozen 并列），需要在 ke-dft-cpu 里引入 QE 提交模板与环境探针。

---

## 4. 各步骤

### Q1 `ivs_qlist`（run: gen，本地）
- 输入：ke-dft-cpu S1 结构；2D 读 S8.2 谷表，3D 读 `DP_OPTICAL_VALLEYS_*`。
- 对每个载流子：初态取每个星的一个代表成员（对称等价成员给出相同 D0），终态取所有星的所有成员；
  q = k_终 − k_初，按晶体对称性约化、合并等价 q。
- 输出 `ivs_qlist.json`：`[{carrier, from_star, from_k, to_star, to_k, q_frac, n_members_same_D}]` 与去重后的 q 集合。
- MoS₂ 电子预期 5–8 个独立 q（K、M 型、Q 型（约 K/2，Q 谷在 0.515 ΓK，不是严格 K/2）、Q–Q′ 的几种）。
- 检查：谷表里非高对称谷（Q）的 k 中心用 S8.2 的网格极值点 `k_grid`，与 S8.4 的谷归属一致。

### Q2 `ivs_scf`（QE）
- 结构：直接用 VASP 弛豫结构（晶格常数与原子坐标原样），**不在 QE 里重新弛豫**，保证谷位置与 AMSET 一致。
- 泛函：尽量与 VASP 一致（当前 MoS₂ 是 PBEsol，`GGA=PS`）；赝势用 PseudoDojo NC（Filipovic 同款，Mo 含半芯）。
- 2D：真空 ≥ 20 Å。谷间耦合是短程、区边界的耦合，Filipovic 的直接 DFPT 也没加 2D 库仑截断。
  是否加 `assume_isolated='2D'`，取决于 Q3 用的矩阵元输出功能是否支持（见 Q3）。
- 收敛档：ECUT 60–90 Ry，12×12×1 k（2D）。
- 自检：QE 的能带里，各谷相对带边的能差与 S8.2 的 `dE_eV` 相差 ≤ 30 meV（例如 MoS₂ 的 Q 谷 129 meV），超出就报错。

### Q3 `ivs_ph`（QE，按 q 扇出，`fanout`）
- 每个 q 一个子目录：ph.x 单 q DFPT（9 支声子；3 原子原胞）。
- 矩阵元：在初态 k（谷中心）和 k+q（终态谷中心）上取 g_{mnν}(k,q)，导带取带边那条（简并时整个简并子空间）。
- **执行前核实（M0）**：QE 7.3 的 ph.x 能否在指定 k 列表上输出 g。候选：
  1. `electron_phonon='ahc'`（Allen–Heine–Cardona 输出，Filipovic 用的就是这个；需要配套的 nscf k 列表）；
  2. Perturbo `calc_mode='ephmat'`（仓库现成，但要先做 Wannier 和粗 q 网格 DFPT，成本高，只作交叉校验）。
- 计算量估算：每个 q 约 1–2 节点小时，MoS₂ 合计约 10–20 节点小时（需试跑确认）。

### Q4 `ivs_read`（run: gen，本地）
- 规范不变量：对初态、终态简并子空间求 Σ|g|²，再按简并度平均（与 `ivs_frozen` 的劈裂法同一口径）。
- 换算：D_ν = |g_ν| · √(2 M_cell ω_ν) / ħ（Filipovic Note S2 与 Kaasbjerg 2012 的定义；与 dp_optical 的
  |g|² = ħD²/(2NMω) 一致，实现时用 `ivs_frozen.py` 头部的归一化说明核对一遍）。
- 归并成 dp_optical 的星通道：
  - `X>Y`（到别的星，dp_optical 用整个 Y 星的态密度）：每个（目标成员 j，声子支 ν）一条通道，D_ch² = D_{ij,ν}² / g_Y；
  - `X>X'`（同星其它成员，用 (g−1)/g 的态密度）：D_ch² = D_{ij,ν}² / (g−1)；
  - 同一对星、同一 ħω（±0.1 meV）的几条按平方和合并（dp_optical 拒绝重复 ħω）；
  - D < `IVS_DFPT_DMIN`（建议 0.05 eV/Å）的丢弃，记录在案。
- 输出 `ivs_dfpt_params.json`，schema `ivs_dfpt/1`：
  ```json
  {"schema": "ivs_dfpt/1", "structure_sha256": "...", "valleys_source": "S8.2 dpt_result.json sha256 ...",
   "provenance": {"qe_version": "7.3", "functional": "PBEsol", "ecut_Ry": 80, "kgrid": [12,12,1], "pseudos": {...}},
   "carriers": {"electron": {"channels": [{"from": "K", "to": "Q", "kind": "inter", "D_eV_A": 3.50, "hw_meV": 49.1,
                                            "branch": 8, "q_frac": [..], "raw": [...]}]}},
   "dropped": [...], "checks": {"valley_dE_vs_S8.2_meV": {...}}}
  ```
  `raw` 保留未归并的逐成员、逐声子支数值，便于和 Filipovic Table S4/S8 逐条对照。

### S8.4（2D）/ S8（3D）接入（改 ke-dft-cpu）
- 新开关 `IVS_DFPT = auto | on | off`、`IVS_DFPT_SRC`、`IVS_DFPT_DMIN`，写法照 `IVS_FROZEN`。
- 优先级：DFPT > 冻声子 K↔K′ > 手填（手填的谷间通道在有 DFPT 结果时只打印对照、不进散射率，与 ODP_FROZEN 的规则一致）。
- 有谷间通道就自动把 IVS 补进 SCATTERING。显式写了 `SCATTERING` 时按 [0034] 规则处理并打 WARN。
- `ADP_INTRAVALLEY=auto` 随之生效，ADP 只留谷内。
- 记录进 `2d_correction.json` 的 `dp_optical_ivs`：mode、source、sha256、通道摘要。参数指纹变了，S8.4 失效重排
  （照抄 `invalidate_ivs_consumer`）。
- S8.3 图注与汇总写明 IVS 参数来源（DFPT / 冻声子 / 文献）。
- 3D：`step8_amset/gen_step10_amset.py` 同样加开关；谷表来自 `DP_OPTICAL_VALLEYS_*`（full3d=True）。

---

## 5. 验证与验收

| 层级 | 内容 | 验收标准 |
|---|---|---|
| 单元测试 | 合成 g → D0 → 星通道：成员占比、简并求和、同 ħω 合并、DMIN 丢弃、schema 校验；照 `tests/test_0050_*` 的风格 | 全部通过；现有测试不回归 |
| 内部交叉 | MoS₂ K↔K′：DFPT 对比本流程冻声子 `ivs_frozen`（若已开） | D 差 ≤ 15% |
| 文献逐通道 | MoS₂ 电子对照 Filipovic Table S4（K→K′ E′ 40.6 meV 1.66；K→Q：LA 22.2/1.15、LO 42.7/1.44、A1′ 49.3/1.04、TA 45.5/0.57）和 Table S8（Q→Q′、Q→K） | D > 1 eV/Å 的通道差 ≤ 30%，ħω 差 ≤ 3 meV |
| 端到端 2D | MoS₂ 用干净的通用配置重跑 S8.4：`ODP_FROZEN=auto`、不手填 `DP_OPTICAL_*`、不显式写 SCATTERING、PIE 与当前版一致 | 声子限制 μ(300 K) 在 150–280；γ（μ∝T^−γ，300–600 K）在 1.4–1.8；高温 ADP-only 不再反常上升 |
| 推广 2D | MoSe₂、WS₂、WSe₂ 电子（S8.2 判多谷"拒绝"的 4 个） | 与 MoS₂ 同样的端到端检查；有文献时对照（Jin 2014 Table IV 只作量级参考，其 K→Q 近邻声学偏大约 7 倍） |
| 3D | Si 电子 X 谷 f/g 过程 | D0 与文献同量级（Jacoboni & Reggiani RMP 1983 的拟合值，需核对原文）；电子 μ(300 K) 从 ADP-only 约 2381 降到约 1400 左右的量级 |

---

## 6. 风险与注意事项

- **矩阵元输出接口**（M0 必须先确认）：ph.x 在 QE 7.3 的电声输出功能、输出格式、是否兼容 2D 截断。确认不了就先停，报告用户。
- **谷位置一致性**：QE（NC 赝势）与 VASP（PAW）对 Q 谷能量可能差几十 meV。只要求矩阵元来自同一结构；
  谷能量仍用 VASP/AMSET 的。Q2 的自检超过 30 meV 时报错，交给用户判断。
- **简并与规范**：只用规范不变量（简并子空间的 Σ|g|²）；声子简并支（E′、E″）同样按简并组求和。
- **声学支在 q≈0 附近**：本方案只取区边界和谷间的 q，不碰小 q，避开 Filipovic 提到的周期 slab 长程假象。
- **SOC**：不支持。遇到 SOC 能带直接拒绝，与 `ivs_frozen` 一致。
- **AutoZT 纪律**：所有计算走 `autozt`（gen + start），不手动 sbatch 或 sed 远程文件；`conf --set`、`rerun`、带 `-f` 的命令都要用户审批。
  新步骤加进 skill.yaml 前要用户同意（铁律 10）。改 step8.x 脚本必须 bump `_SKILL_REV`，并在 METHODOLOGY 加 `[00xx]` 条目。

---

## 7. 里程碑

| 里程碑 | 内容 | 产出 |
|---|---|---|
| M0 | 核实 jzzn 上 QE 7.3 ph.x 的矩阵元输出功能、PBEsol NC 赝势是否可用 | 一页核实记录；不可行则停 |
| M1 | Q1 + Q4 的纯 numpy 部分加单元测试（不依赖集群） | 通过的测试 |
| M2 | Q2/Q3 接入 eph-qe-cpu（或 ke-dft-cpu，按确认 A），MoS₂ 只算 q = K（K↔K′） | 和冻声子、Filipovic 1.66 对照 |
| M3 | MoS₂ 全部 q，Q4 输出 `ivs_dfpt_params.json`，逐通道对照文献 | 对照表 |
| M4 | S8.4 接入（`IVS_DFPT`），MoS₂ 端到端重跑（审批后） | μ(T) 与基准对照 |
| M5 | MoSe₂、WS₂、WSe₂ | 同上 |
| M6 | 3D 接入（S8）+ Si 验证 | 同上 |

---

## 8. 待用户确认

- A. 新步骤放 eph-qe-cpu（推荐）还是 ke-dft-cpu；跨技能取文件的方式。
- B. 是否同意给技能新增步骤（铁律 10）。M0、M1 不碰集群也不改 skill.yaml，可以先做。
- C. 验收标准（第 5 节）是否认可，特别是"与 DFPT 基准差 20–30% 即通过"。
- D. 集群与核数：jzzn 4 核测试模板只适合冒烟测试，正式 DFPT 用多少核、哪台机器。
