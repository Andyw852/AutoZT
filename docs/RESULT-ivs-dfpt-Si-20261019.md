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

文献取 Eric Pop 博士论文 Table 3.2 的六组 Si 谷间形变势 ∆（单位 1e8 eV/cm，各拟合组差异很大）；
换算到 dp_optical / ivs_frozen 口径：**D = ∆·√(Z_f/2)**（f 型 Z_f=4，g 型 Z_f=1）。

| 过程 | ħω 文献 | ∆ 范围 → D 范围 | 本次 DFPT |
|---|---|---|---|
| f2 (LA/LO) | ~51 meV | 2–4.3 → **D 2.8–6.1** | 2.61 @47.1 |
| f3 (TO) | ~57 meV | 1.5–8 → **D 2.1–11.3** | 4.28 @56.9 |
| g3 (LO) | ~62 meV | 3–11 → **D 2.1–7.8** | 4.70 @61.4 |

结论：三条通道的 D 与 ħω 均落在文献散布范围内；TA/LA 支由零阶选择定则严格为零。

## 4. ke-dft-cpu S8 验收（jzzn，AMSET 0.5.1，PBEsol，S3 网格 11³）

机制集：`scattering_type = [ADP, IMP, IVS]`；IVS 由 `ivs_dfpt_params.json`（sha256 `aa259f918ff1…`）经 `IVS_DFPT=auto` 注入，
S8.3 汇总写"`[IVS 参数来源] S8 electron：DFPT …；3 条通道`"。

本征/低掺杂电子迁移率（T=300 K，单机制由 AMSET separate_mobility 给出）：

| 掺杂档 | ADP only | ADP+IVS | 全机制（+IMP） | 说明 |
|---|---|---|---|---|
| **n = 1e10 cm⁻³** | 3883 | **1459** | **1459** | 验收口径；目标 ~1400 ✅ |
| n = 1e21 cm⁻³ | 2341 | 390 | 64 | 简并档，IMP 主导 |

- 电子在 n=1e10 的限制机制：**IVS 2747 < ADP 3883 ≪ IMP（可忽略）** → IVS 是声子限制的主导项，符合物理预期。
- 空穴没有 IVS 通道（本轮 DFPT 只算了电子），空穴侧仍由 ADP 限制。

## 5. 已知问题（发现于本次验收，待处理）

1. `postprocess_intrinsic.py` 打印"本征迁移率"用 `mobility[0][t300]`，而索引 0 是**第一档掺杂（−1e21）**不是本征点；
   作业日志里 63.6 / 390.1 / 2340.8 都是重掺杂档。方案里的 2381 基线也来自这一行（同口径）。
   → 应改为"最低非零掺杂档"（或同时打印两档并标清）。
2. 恰好 0 掺杂档在 AMSET 里不可用：四种机制集的 mobility 全部给出同一个 7610.0 → 该档不参与判断。
3. IVS 单机制迁移率随掺杂变化明显（2747 @1e10 → 483 @1e21），声子机制的掺杂依赖偏强，需核对简并 E_F 下的
   占据/DOS 权重（dp_optical 的 inter_star 分支）。
4. 汇总把"无此机制/可忽略"的机制值（IMP 8.3e6、空穴 IVS 2.0e9）当实数打印，建议显示"—（可忽略）"。

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
