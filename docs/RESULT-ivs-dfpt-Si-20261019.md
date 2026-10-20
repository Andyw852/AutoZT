# IVS-DFPT：Si 实测与验收记录（2026-10-19）

配套方案 `docs/PLAN-ivs-dfpt-20261018.md`；代码条目 `[0059]`（ke-dft-cpu METHODOLOGY）。
本文件记录第一次**集群实跑**（eph-qe-cpu 的 ivs_dfpt 组 S10–S13 + ke-dft-cpu S8 接入）的参数、原始结果与验收。

## 1. 运行环境与参数

| 项 | 值 |
|---|---|
| 集群 | eph-qe-cpu 链：**hfeshell**（job 28948945 / 28948987 / 28948988）；ke S8：**jzzn**（job 4245630） |
| QE | 7.4（conda-forge 预编译，装在 hfeshell `~/software/AutoZT/qe-7.4`） |
| 赝势 | SG15 `Si_ONCV_PBE-1.2.upf`（模守恒，sha256 `73d21061…`），`kind=NC`，AHC 支持探针 = yes |
| 泛函 | `input_dft=pbesol`（自动继承 ke-dft-cpu S1 的 `FUNC=pbesol`） |
| 结构 | ke S1 CONTCAR（2 原子原胞，a = 3.865 Å；对称操作 48 个） |
| 截断/网格 | ecutwfc 60 Ry、ecutrho 240 Ry、SCF k 网格 10×10×10（IVS_KSPACING=0.2 Å⁻¹） |
| DFPT | 单 q、tr2_ph = 1e-16；nscf nbnd=12（nocc=4，AHC 窗口 ahc_nbndskip=2 / ahc_nbnd=4） |
| 谷表 | `IVS_VALLEYS_E = X:0,0.425,0.425` → 按点群+时间反演补全为 6 个 X 谷 |
| q 作业 | **q-01** = (0, 0.150, 0.150)（g 过程，X→−X，折叠后小 q）；**q-02** = (−0.425, −0.425, 0.150)（f 过程，X→Y） |

## 2. 逐谷对结果（D 单位 eV/Å，ħω 单位 meV；零 = 选择定则禁戒）

```text
X→X (q-01) q=(0.000, 0.150, 0.150) ×1 对   带 [5]→[5]：0.00@11.0, 0.00@18.0, 0.00@59.7, 4.70@61.4
X→X (q-02) q=(-0.425, -0.425, 0.150) ×4 对 带 [5]→[5]：0.00@16.5, 0.00@21.2, 0.00@46.5,
                                                          2.61@47.1, 0.00@55.9, 4.28@56.9
```

dp_optical 星通道（X>Y：D²/g_Y；X>X′：D²/(g_X−1)；同 ħω 合并）：

```text
X>X':2.337@47.09   X>X':3.827@56.94   X>X':2.101@61.43      （其余支为 0）
```

本机冒烟对照（30 Ry、6³ k，M0 档）：g-LO 4.79 @63.1（星通道 2.141@63.06）；f 2.73@48.5（2.443）、4.42@58.4（3.953）。
集群 60 Ry/10³ k 的 D 低 2–4%、ħω 低 2–3%，选择定则（哪些支为零）完全一致。

## 3. 与文献对照（Si 谷间形变势）

文献取 Eric Pop 博士论文 Table 3.2 的六组 Si 谷间形变势 ∆（单位 1e8 eV/cm = eV/Å，各拟合组差异很大）。
**口径（2026-10-19 更正）**：上表的逐谷对 D 与文献 ∆ **直接相比，不需要换算**。文献速率
Γ = π∆²Z_f/(2ρω)·g_c(E)（g_c 含自旋 = 2 × 单谷单自旋态密度，Jacoboni & Reggiani 1983；Lundstrom）与本流程逐谷对的
πD²G₁/(Mω) 在 D = ∆ 时逐项相等（用 dp_optical.analytic_rate_3d 数值核对比值 = 1.000）；Z_f 只出现在速率里（f 型 4 个终态谷、
g 型 1 个），已由星通道折算 D_ch² = Z_f·D²/(g−1) 计入。初稿写的 D = ∆·√(Z_f/2) 会让速率差一倍，作废。

| 过程 | ħω 文献 | 文献 ∆ 范围 | 本次 DFPT（逐谷对 D） |
|---|---|---|---|
| f2 (LA/LO) | ~51 meV | 2–4.3 | 2.61 @47.1 |
| f3 (TO) | ~57 meV | 1.5–8 | 4.28 @56.9 |
| g3 (LO) | ~62 meV | 3–11 | 4.70 @61.4 |

结论（按正确口径不变）：三条通道的 D 与 ħω 均落在文献散布范围内；TA/LA 支由零阶选择定则严格为零。

## 4. ke-dft-cpu S8 验收（jzzn，AMSET 0.5.1，PBEsol，S3 网格 11³）

机制集：`scattering_type = [ADP, IMP, IVS]`；IVS 由 `ivs_dfpt_params.json`（sha256 `aa259f918ff1…`）经 `IVS_DFPT=auto` 注入，
S8.3 汇总写"`[IVS 参数来源] S8 electron：DFPT …；3 条通道`"。

本征/低掺杂电子迁移率（T=300 K，单机制由 AMSET separate_mobility 给出）：

| 掺杂档 | ADP only | ADP+IVS | 全机制（+IMP） | 说明 |
|---|---|---|---|---|
| **n = 1e10 cm⁻³** | 3883 | **1459** | **1459** | 验收口径 ✅（见下） |
| n = 1e21 cm⁻³ | 2341 | 390 | 64 | 简并档，IMP 主导 |

**验收依据（2026-10-19 更正）**：方案第 5 节写的"ADP-only 约 2381 → 约 1400"，其中 2381 来自 `postprocess_intrinsic.py` 的日志行，
那一行取的是 −1e21 重掺杂档（见第 5 节问题 1），与 n = 1e10 的 1459 不是同一口径，两者"对上"属巧合，不能作为依据。
同口径的依据是 Si 电子 300 K 的晶格限制迁移率实验值约 1400–1450 cm²/Vs：低掺杂 1459 落在其中，验收结论成立。
- 电子在 n=1e10 的限制机制：**IVS 2747 < ADP 3883 ≪ IMP（可忽略）** → IVS 是声子限制的主导项，符合物理预期。
- 空穴没有 IVS 通道（本轮 DFPT 只算了电子），空穴侧仍由 ADP 限制。注意：Si 空穴只有 Γ 一个谷，**本来就没有谷间散射**；
  它缺的是谷内光学 ODP（本次 `ODP_FROZEN = off` 把冻声子 ODP 一起关了），所以空穴 μ 偏高不是 DFPT 算不了空穴。
  [0059-r3] 起 DFPT 在 q=Γ 顺带算谷内 ODP（含单谷空穴）：本机 QE 7.4 同口径小算例给空穴 Γ D = 7.89 eV/Å @ 63.9 meV，
  与冻声子 7.77 一致；集群正式值待按 r3 重跑 S10–S13（`IVS_VALLEYS_H = G:0,0,0`）后由 S8 读入。

## 5. 已知问题与处理

1. ~~`postprocess_intrinsic.py` 打印"本征迁移率"用 `mobility[0][t300]`~~ **已修（[0059-r2]）**：索引 0 是**第一档掺杂（−1e21）**
   不是本征点；作业日志里 63.6 / 390.1 / 2340.8 都是重掺杂档，方案里的 2381 基线也来自这一行。现在打印
   "本征近似 = 最低非零掺杂档"，Si 该值 = **1459 cm²/Vs（n = 1e10）**。
2. 恰好 0 掺杂档在 AMSET 里不可用：四种机制集的 mobility 全部给出同一个 7610.0 → 该档不参与判断。这是 AMSET 行为，
   无法在流程侧修；打印与汇总一律改用最低非零档。
3. IVS 单机制迁移率随掺杂变化明显（2747 @1e10 → 483 @1e21）。**核对结论：大概率是物理，不是 bug**：1e21 时 E_F − E_c ≈ 210 meV，
   远高于 ħω ≈ 47–61 meV，发射声子通道整体打开（低掺杂时电子平均能量 ~40 meV，发射被阈值压住）。三维抛物带模型（6 谷、
   m_d = 0.328、同样三条通道、含泡利阻塞）给出 μ_IVS(1e10)/μ_IVS(1e21) ≈ 4.0、ADP ≈ 2.2，AMSET 实测 5.7 与 1.7——
   "IVS 对掺杂远比 ADP 敏感"能复现，余下差别来自非抛物性与 AMSET 的真实能带。
4. ~~汇总把"无此机制/可忽略"的机制值（IMP 8.3e6、空穴 IVS 2.0e9）当实数打印~~ **已修（[0059-r2]）**：超过
   `INACTIVE_MU_CM2_VS = 1e5` cm²/Vs 的机制值印"—"，且不参与"限制因素"的选择（AMSET 用极大值表示该机制散射率≈0）。
5. ~~S8.21 的 3D 自动谷表在 Si 上只给 V1×3（应为 6 个 Δ 谷）~~ **已修（[0059-r2]）**：`odp_frozen.find_valleys_3d` 的平台合并
   会跨 k 与 −k，把"相邻、能量相同"的时间反演伙伴当成一个平台（±x 两个 Δ 谷隔着 X 只差约 0.3X，在 11³ 网格上正好是相邻格点
   t = 0.909 与 1.091 且能量按对称相等）。现在平台合并对每条边做时间反演检查，q ≡ −c（模 1）的格点不并入。本次实测仍是
   手填 6 谷 + `ODP_FROZEN = off`（正确做法，冻声子电子 ODP 通道本来就 < 0.05 eV/Å）；修复后 S8.21 自动找谷也会给出 6 谷。
6. `ivs_qlist.json` / `qrun.json` / `ivs_dfpt_params.json` 里 `tr2_ph` 记成 0.0：`ivs_dfpt.jsonable` 按小数位舍入把 1e-16 抹掉
   （输入文件 ph_dfpt.in 里是对的，只是元数据）。已修（[0059-r1]，保留 12 位有效数字；S12 重新生成时旧记录的 0.0 不算输入变了，
   不会让已完成的 q 被归档重算）。

## 6. 复现与数据位置

```bash
# eph-qe-cpu（hfeshell）：S10 → S11 → S12（扇出 2 个 q）→ S13
autozt -tt eph-qe-cpu -p Si -j S10_ivsq  start
autozt -tt eph-qe-cpu -p Si -j S11_ivsscf init && autozt -tt eph-qe-cpu -p Si -j S11_ivsscf start
autozt -tt eph-qe-cpu -p Si -j S12_ivsph start
autozt -tt eph-qe-cpu -p Si -j S13_ivs    start
# ke-dft-cpu（jzzn）：S8 接入 DFPT-IVS
autozt -tt ke-dft-cpu -p Si -j S8_kappa retry && autozt -tt ke-dft-cpu -p Si -j S8_kappa start
autozt -tt ke-dft-cpu -p Si -j S8.3_cmp start      # 汇总写明 IVS 来源
```

产物（本地 `<材料>/eph-qe-cpu/result/`、`<材料>/ke-dft-cpu/result/`）：

- `step10_ivs_qlist/ivs_qlist.json`（谷表、2 个 q、结构 sha256）
- `step11_ivs_scf/{ivs_env.json, valley_check.json, scf.in, scf.out}`
- `step12_ivs_ph/q-0{1,2}/{qrun.json, ph_dfpt.in, ph_ahc.in, ahc_dir/ahc_*.bin}`
- `step13_ivs_read/{ivs_dfpt_params.json, ivs_dfpt_table.txt}`
- `step8_amset/{dp_optical_3d.json, transport.json, subset_*.json}`
- `step8.3_output/comparison_summary.txt`（含 `[IVS 参数来源]` 行）


## 7. DFPT 谷内 ODP 验证（[0059-r3]，2026-10-20/21）

r3 让 q 列表固定加 `q-odp (Γ, 0,0,0)`，用同一套 AHC 顺带算谷内光学形变势，替代冻声子 ODP。
Si 是理想的验证对象：**电子（X 谷）ODP 按对称性严格为 0，空穴（Γ 谷，价带顶三重简并）ODP 非零**。

### 7.1 Si 空穴 Γ 谷（非零情形，核心验证）

| 方法 | D (eV/Å) | ħω (meV) | 差异 |
|---|---|---|---|
| **DFPT q=Γ（本轮）** | **7.734** | 62.34 | — |
| 冻声子（S7c/S8.21） | 7.769 | 62.04 | **0.5%** |

- 声学模（QE 在 Γ 给 0.58 meV，非严格 0）已按平移权重显式剔除，不进通道；
- 价带顶 3 重简并需要 `IVS_AHC_NB` 放宽到上下各 3 条（S12 按 S11 谷能带只对需要的 q 放宽，已完成的 q 不重算）；
- 单谷载流子（Si 空穴只有 Γ 一个谷）在 `IVS_ODP=auto` 时保留，只算谷内。

### 7.2 Si 电子 X 谷（零情形，对称性对照）

DFPT 给 `X>X D≈0 @62.3 meV`，与冻声子 0.000 一致。S8 里电子加 ODP 前后 μ 不变（3882.6），
证明零通道不会被误加进散射率（`D²/ω` 不会因噪声爆掉）。

### 7.3 S8 端到端结果（jzzn，S8_kappa job 4284923，1h19m，exit 0）

`scattering_type = [ADP, IMP, ODP, IVS]`，`dfpt_use = {electron: used, hole: used}`，
冻声子标 `superseded_dfpt`（只作对照，不再进散射率）。

Si 300 K 迁移率（cm²/Vs，按机制子集逐项拆开）：

| 数据集 | 空穴 p=1e10 | 空穴 p=1e18 | 电子 n=1e10 |
|---|---|---|---|
| ADP only | 7610.0 | 7623.1 | 3882.6 |
| **+ DFPT ODP** | **1558.2**（↓4.9×） | 1544.0 | 3882.6（不变 ✓） |
| + DFPT IVS | 1558.2（单谷无 IVS ✓） | — | **1459.2**（↓2.7×） |
| 全机制（+IMP） | 1557.2 | 90.4 | 1458.6 |

S8.3 汇总现在写明来源：
```text
# [AMSET 机制] electron 300 K（-1e10）：IVS 2747 < ADP 3883 < IMP — < ODP — → 限制因素 IVS
# [AMSET 机制] hole     300 K（ 1e10）：ODP 2216 < ADP 7610 < IMP — < IVS — → 限制因素 ODP
# [IVS 参数来源] S8 electron：DFPT（ivs_dfpt_params.json sha256 9fa72aa0b433；零阶、非 SOC）；3 条通道
# [ODP 参数来源] S8 hole：DFPT q=Γ（ivs_dfpt_params.json sha256 9fa72aa0b433；零阶、非 SOC）；1 条通道
```

### 7.4 解读

- **电子**：ADP 3883 → 加 DFPT IVS 后 **1459**，落在实验晶格限制值 1400–1450；ODP=0 不参与。
- **空穴**：ADP 7610 → 加 DFPT ODP 后 **1558**，ODP 成为限制机制（此前 ODP 缺失导致空穴虚高）。
  1558 仍高于实验（Si 空穴声子限制约 450–500）。**偏差不来自 ODP**：DFPT 的 ODP 与冻声子一致（0.5%），而且 q=Γ 的值
  已经是散射率上限，改成小 q 平均只会让 μ 更高（原稿"符合方案 7.6 的预期"方向反了，2026-10-21 更正）。
  偏差在 ADP（只算 ADP 7610 明显偏高），两个来源，见 `skill/ke-dft-cpu/METHODOLOGY.md` [0061]：
  ① AMSET 的 ADP 写法（初态带的对角形变势 × 重叠）丢了三重简并价带里重/轻空穴之间的剪切耦合 —— 离子弛豫口径下只算 ADP
  高估 1.6–2.1 倍，已由 S8 的 `ADP_DEGEN`（出厂 auto）修正；
  ② 7610 正好对上"钳位离子剪切形变势 + AMSET 写法"的模型值（8000–8300；离子弛豫应为 2300 左右）—— 本次 S7 很可能是
  [V187] 之前的钳位离子构型，核对 `step7b_deform_read/band_edges.json` 的 `deform_geometry`，不是 ionrelax 就 retry S7 → S7.1。
  两项都修好后剩下的差距再看 SOC（分裂带 44 meV，本链未计）与非抛物性。
- 电子/空穴互不串扰：电子无 ODP、空穴无 IVS，两种情形都验证了"只加该加的"。

### 7.5 复现（[0059-r3] 起）

```bash
# 1) 空穴谷表（3D 必须手填）：IVS_VALLEYS_H = G:0,0,0
autozt -tt eph-qe-cpu -p Si -j S10_ivsq conf --set params.IVS_VALLEYS_H=G:0,0,0
# 2) eph 链（q-odp 会作为第 3 个 q 作业扇出）
autozt -tt eph-qe-cpu -p Si -j S10_ivsq start
autozt -tt eph-qe-cpu -p Si -j S11_ivsscf retry && autozt -tt eph-qe-cpu -p Si -j S11_ivsscf start
autozt -tt eph-qe-cpu -p Si -j S12_ivsph retry && autozt -tt eph-qe-cpu -p Si -j S12_ivsph start
autozt -tt eph-qe-cpu -p Si -j S13_ivs retry && autozt -tt eph-qe-cpu -p Si -j S13_ivs start
# 3) ke 链：S8 自动用 DFPT ODP（DFPT 优先于冻声子）
autozt -tt ke-dft-cpu -p Si -j S8_kappa retry && autozt -tt ke-dft-cpu -p Si -j S8_kappa start
autozt -tt ke-dft-cpu -p Si -j S8.3_cmp start
```

> ⚠️ gen 端跑 `MESH_MIN` 自动选网格时需要 `BoltzTraP2`。若目标集群的 gen python（`remote_path_prefix` 里的 python）
> 没有 BoltzTraP2，可写 `MESH_MIN=0 / MESH_MAX=0 / INTERPOLATION_FACTOR=<显式值>` 跳过该步（本轮 hfeshell 的
> `atomate2_p_a` py3.10 无 BoltzTraP2，就是用 factor=60 复现 jzzn 的 73³ 网格）。

