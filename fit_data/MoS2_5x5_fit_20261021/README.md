# MoS2 5×5×1 力常数拟合数据包（ALASSO vs RFE-OLS）

生成时间：2026-10-21
来源：AutoZT 技能 `fit-fc-thermal`，材料 `MoS2_5x5_fit`，集群 3090
（`/data/wangchaoyue852/taskflow/work/MoS2_5x5_fit/fit-fc-thermal/`）

## 输入数据

- 结构：MoS2 单层 1H，原胞 3 原子，5×5×1 超胞 = 75 原子，25 Å 真空。
- 数据集：36 帧随机位移 + 1 帧平衡态（`MoS2_kl_5x5_dataset`），位移 RMS 0.0171 Å。
- 拟合：pheasy，二阶 + 三阶力常数；`FIT_METHODS = pheasy:ALASSO,RFE-OLS`（ALASSO 为主方法）。

## 两种方法设置

| 项 | 值 |
|---|---|
| 三阶截断候选 | 3.57 / 4.52 / 5.27 / 5.90 / 6.54 Å |
| CUT3_SCAN | auto（每档单独重拟合） |
| CUT3_BOOTSTRAP | **0**（本包未做逐壳层 bootstrap，截断选择只按 κ 平台判据） |
| 公式/参数 | 见 `step1_fit/fit_config.json`、`step1_fit/methods/m-pheasy-RFE-OLS/fit_config.json` |

## 关键结果（300 K，面内）

| 截断 (Å) | ALASSO 拟合相对误差 | RFE-OLS 拟合相对误差 | ALASSO κ (W/m·K) | RFE-OLS κ (W/m·K) |
|---|---|---|---|---|
| 3.57 | 0.589% | 0.608% | 35.62 | 33.44 |
| 4.52 | 0.369% | 0.387% | 32.63 | 35.65 |
| 5.27 | 0.276% | 0.305% | 25.17 | 25.00 |
| 5.90 | 0.247% | 0.288% | **22.00** | **25.21** |
| 6.54 | 0.199% | 0.255% | 13.33 | 15.81 |

- ALASSO 的 CV 把 α\* 选为 0，**实际交付的是 OLS 解**（`fit_method_configured=ALASSO`、
  `fit_method_delivered=OLS`、`alpha_star_source=cv_conclusion`，见对应 `fc_fit_summary.json`）。
- RFE-OLS 交付稀疏解：5.90 Å 档 894/2304 个非零 IFC，6.54 Å 档 1177/3755。
- 两方法在 5.90 Å 档都被选为最终截断；ALASSO 判 `not_converged`（无 κ 平台），
  RFE-OLS 判 `ok`（5.27–5.90 Å 平台，+0.8%）。
- 2D 归一化 κ（h⊥/d = 3.7092）：ALASSO 81.60，RFE-OLS 93.53 W/m·K。

## 目录

```
step1_fit/                     # ALASSO（主方法）的拟合输入/配置与结果
  fit_config.json              # 拟合参数（超胞、截断候选、引擎、方法等）
  methods.json                 # 多方法清单（谁在 methods/ 下）
  fc_dataset.json              # 数据集元信息（帧数、原子数、平衡态等）
  disp_matrix.pkl              # 位移（特征）矩阵
  force_matrix.pkl             # 力矩阵
  neighbor_list.pkl            # 近邻表（特征定义）
  cutoff_scan.json             # 5 档截断扫描结果
  fc_fit_summary.json          # 拟合汇总（含配置/交付方法、α*、残差、稳定性）
step1_fit/m-pheasy-RFE-OLS/    # RFE-OLS 的同一套文件
step2_kappa/                   # ALASSO 的 κ 结果（kappa_summary、κ–截断图/数据、选中截断）
step2_kappa/m-pheasy-RFE-OLS/  # RFE-OLS 的 κ 结果
step2_kappa/figures/           # 两方法对比图 kappa_vs_cutoff_<tag>.{png,pdf}
step2_kappa/methods_compare.json
```

## 未包含（体积大，如需另行获取）

- `fc2.hdf5` / `fc3.hdf5` / `FORCE_CONSTANTS*`（各档力常数，单个 ~0.06–3.4 MB）
- `sm_dense.npy`（75 MB）、`sm_prime.npz`（95 MB）、`cs.pkl`（13 MB）等中间矩阵
- 各档 `cutoff_scan/cut3_*/` 下的 fc2/fc3

## 说明

- 本包由 `tmp/mos2_kappa_local` 的本地重算与 3090 作业 6206 交叉验证：
  κ 数值逐位一致（ALASSO 81.60117266011453、RFE-OLS 93.52665397475592，2D 归一）。
