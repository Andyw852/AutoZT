# MoS2（单层 2H，2D TMD）—— ke / eph 关键结果

- 材料：MoS₂ 单层（2H 相，P-6m2），3 原子/胞，真空层 c = 34.879 Å
- 技能链：eph-qe-cpu（QE DFPT：谷间 IVS + 谷内 ODP）→ ke-dft-cpu（AMSET 电子输运 + BT2 + DPT）
- **导出日期：2026-10-24（快照；空穴 DFPT 与最终 AMSET 重跑尚未完成，见下方"进行中"）**
- 溯源：
  - 集群 jzzn：`/public/home/wangchao/Fullerene_Network/work/MoS2/{ke-dft-cpu,eph-qe-cpu}`
  - 集群 clus_xa：`/work/home/wangchao217/software/AutoZT/work/MoS2/eph-qe-cpu`（epf DFPT 与 AHC）
  - 本地完整数据：`/mnt/d/tf_data/MoS2/`
  - Windows 归档：`D:\joint_research_project\AutoZT\基准测试\MoS2\`

## 方法要点

| 项 | 设置 |
|---|---|
| 泛函 | PBEsol（GGA=PS，无色散修正） |
| 带隙（剪刀） | HSE06/HSEsol **2.3239 eV**（源：step2.3_hse_plot/band_summary.json） |
| 电子-声子 | QE 7.4 DFPT + AHC：4 个 q（K↔K′、K→Q 近/远、q=Γ），tr2_ph=1e-16 |
| 形变势 | 冻声子（S7，9 个应变构型，ionrelax 面内 + clamped 面外） |
| 散射机制 | ADP + POP + ODP(DFPT) + IVS(DFPT) + IMP |

## 关键数值（300 K）

### 声子限制迁移率（AMSET，低掺杂）

| 机制 | 电子 μ (cm²/V·s) | 空穴 μ (cm²/V·s) |
|---|---:|---:|
| ADP | 594.3 | 1976 |
| POP | 1759 | 1462 |
| ODP | 1451 | 3420 **（冻声子）** |
| IVS | 7385 | — **（未算）** |
| **全声子限制** | **165** | **608** |

电子的限制因素是 **ADP**；空穴是 **POP**。

### 输运与 ZT（κ_L = 17.515 W/m·K，已按胞高换算）

| n (cm⁻³) | 类型 | μ | S (μV/K) | κ_e (W/m·K) | **ZT** |
|---|---|---:|---:|---:|---:|
| −1e20 | 电子 | 18.07 | −104.6 | 0.159 | 0.00538 |
| −1e19 | 电子 | 19.07 | −305.3 | 0.0139 | 0.00487 |
| −1e18 | 电子 | 24.94 | −515.2 | 0.0016 | 0.00181 |
| +1e18 | 空穴 | 32.07 | +562.6 | 0.0032 | 0.00278 |
| +1e19 | 空穴 | 20.23 | +349.0 | 0.0219 | 0.00675 |
| +1e20 | 空穴 | 21.38 | +132.7 | 0.222 | **0.01020** |

### DFPT 形变势参数（电子）

| 通道 | D₀ (eV/Å) | ħω (meV) | 说明 |
|---|---:|---:|---|
| ODP K→K | 5.760 | 50.54 | 谷内，DFPT q=Γ |
| ODP Q→Q | 3.909 | 50.54 | 谷内，DFPT q=Γ |
| IVS K→K′ | 1.731 | 41.89 | 最强谷间通道（文献 1.66 @ 40.6） |
| IVS K→Q | 3.845 / 3.431 | 51.04 / 29.31 | 两条 K→Q 分支 |
| IVS Q→Q′ | 2.638 | 23.41 | — |

共 **27 条 IVS 通道 + 2 条谷内 ODP**（`eph/ivs_dfpt_params.json`）。

### 交叉验证

| 项 | 冻声子 | DFPT | 差 |
|---|---:|---:|---:|
| ODP K | 5.803 @ 50.33 meV | **5.760 @ 50.54** | −0.7% |
| ODP Q | 3.818 @ 50.33 meV | **3.909 @ 50.54** | +2.4% |

## 进行中（本快照之后会更新）

1. **空穴 DFPT**（S12 跑中）：空穴谷表 K×2 + Γ×1，q-01/q-odp 的 DFPT 正在 clus_xa 重算（64 核）
   → 完成后空穴的 ODP/IVS 也从 DFPT 取（现在空穴 ODP 还是冻声子、IVS 空缺）
2. **S8.4 AMSET 重跑**：将用上 ①最新 ADP 提取（step7b 已按 `2026-10-22-0061m-strain-label` 重跑）
   ②最新 DPT 代码（`2026-10-23-0061d-dptdegen`）③空穴 DFPT 通道
3. 因此上表的**电子 μ/ZT 会小幅变化**，空穴侧变化较大

## 文件

```
ke/  dpt_result.json       S8.2 形变势迁移率（DPT，含 K/Q 谷表与 m*）
     transport.json        S8.4 AMSET 全结果（overall + 各机制 mobility）
     2d_correction.json    S8.4 散射参数来源记录（ODP/IVS 来源、通道数、sha256）
     band_edges.json       S7.1 带边形变势 E1（版本 0061m-strain-label）
     boltztrap_summary.txt BoltzTraP2 CRTA 摘要
     comparison_300K.csv/.png、comparison_summary.txt、comparison_vs_T.csv    S8.3 汇总
eph/ ivs_dfpt_params.json  DFPT 谷间/谷内参数（schema ivs_dfpt/1；当前为电子版）
     ivs_dfpt_table.txt     同上的人类可读表
     ivs_qlist.json         IVS q 列表（4 个 q 作业）
```

中间文件（wavefunction.h5 ~1.3 GB、WAVECAR、OUTCAR、DFPT 的 dyn/phsave 等）未收录。
