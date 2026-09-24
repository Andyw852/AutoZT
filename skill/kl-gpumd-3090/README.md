# kl-gpumd-3090 —— 晶格热导率（GPUMD MD 版：NEP + HNEMD / Green-Kubo）

与 `kl-dft-cpu`（VASP 力 → fc2/fc3 → ShengBTE）和 `kl-mace-gpu/cpu`（MACE 力 → fc →
phono3py）**并列的第三条路线**：不经过力常数截断，用 NEP 势直接做 MD 算 κ。

适用场景（也是它存在的理由）：**强非谐 / 玻璃型导热材料**——ωτ ≲ 1、RTA 已经失效。
> 注："玻璃型"指**玻璃型热输运**（声子寿命 ≈ 最小 MFP、κ 逼近 κ_sg），**不是非晶结构**：
> Mg4C60（= Mg8C120，C2/c）是**晶体**。
判据是实测的：Mg4C60 的 `kl-dft-cpu` 结果里 κ_RTA(c3=4.5 Å) = 0.325 已经**低于**
小晶粒极限 κ_sg ≈ 0.455 W/mK，说明 IFC+BTE(RTA) 的声子寿命短于物理最小 MFP，
截断越大力常数越不确定，κ 只会被噪声推着往下掉。MD 版不受这个限制。

## 用法

```bash
tf -tt kl-gpumd-3090 -p <材料> init
tf -tt kl-gpumd-3090 -p <材料> -j step1_struct conf --set \
   params.DATA_DIR=/home/user_3090/Data/Mg4C60_rattle
tf -tt kl-gpumd-3090 -p <材料> -j step2_nep conf --set \
   params.PRETRAINED_NEP=nep89_20250409.txt \
   params.PRETRAINED_RESTART=nep89_20250409.restart
tf -tt kl-gpumd-3090 -p <材料> start        # S1 gen → S2 作业 → S3 作业 → S4 画图
```

四步：

| 步骤 | 位置 | 产物 | 判据 |
|---|---|---|---|
| S1_struct | 登录节点（gen） | `gpumd_params.json`、`struct_summary.json` | `"STRUCT_DONE": true` |
| S2_nep | 3090 CPU 队列 | `train.xyz/test.xyz`、`nep.txt`、`nep_summary.json` | `"NEP_DONE": true` 且力 RMSE（相对/绝对）≤ 门槛 |
| S3_kappa | 3090 CPU 队列 | `kappa.out`、`kappa_summary.json` | `"KAPPA_DONE": true` |
| S4_plot | 登录节点（gen） | `kappa_T.png`、`compare_summary.json` | — |

## κ 协议（S3）

每个温度跑三次 HNEMD（GPUMD 文档：驱动方向 μ 只给 κ_μx/κ_μy/κ_μz，所以对角元要三次）：

```
potential nep.txt
velocity 300
time_step 1
ensemble nvt_ber 300 300 100     # nvt_ber <T_init> <T_final> <τ/Δt>（v5.6 无 nvt_berendsen）
run 200000                          # 平衡
compute_hnemd 1000 1e-5 0 0         # 输出间隔 与 Fe(Å^-1)
run 2000000                         # 生产
```

- 对角元：x 驱动取 `κ_xx`（col1+col2）、y 驱动取 `κ_yy`（col3+col4）、z 驱动取 `κ_zz`（col5）。
- 误差棒：丢前 50% 瞬态后做 5 段块平均（`block_average`）。
- 收敛矩阵（务必做）：`CELL_REPLICATE`（512→2048 原子）、`PROD_STEPS`、`DRIVING_FORCE` 各扫一遍；
  HNEMD 的 κ 不应随 Fe 变化（线性响应区）。
- EMD 对照：`EMD_ENABLED=True` 时额外跑 `compute_hac`。`hac.out` 列映射已核对：
  col7-11 = κ_x^in/κ_x^out/κ_y^in/κ_y^out/κ_z^tot，对角元 **in+out 相加**
  （κ_xx=col7+8, κ_yy=col9+10, κ_zz=col11；只取 in 会让 running κ 开头为负）。
- 量子修正：默认用 Debye 近似因子 `c_v^q/c_v^cl`（给 `DEBYE_T` 才生效）。
  严格做法是谱分解修正：`SHC_ENABLED=True` 时额外跑 `compute_shc`（NVE + 谱热流），
  得到 κ(ω) 后用 f(ω)=x²eˣ/(eˣ−1)² 逐频修正（每方向一次，写 `kappa_spectral_quantum_corrected`）。
  语法已按 GPUMD 源码核对：`compute_shc <采样1-50> <Nc 100-1000> <方向 0/1/2> <num_omega> <max_omega THz>`。

## 前置条件（跑之前必须满足）

1. **GPUMD v5.6**：3090 机 `~/gpumd/src/{gpumd,nep}` 已编译（step.conf 里是绝对路径）。
2. **微调基座**：GPUMD 自带 `potentials/nep/nep89_20250409.{txt,restart}`，但**当前 3090 机的
   `~/gpumd/potentials/` 里只有两个纯碳 NEP4，没有它** → 需要先补上；或把
   `PRETRAINED_NEP` 留空做从头训练（数据量不大时精度会差一些）。
   架构键（`version/zbl/cutoff/n_max/basis_size/l_max/neuron`）**必须与基座一致、不可改**。
3. **训练数据**：VASP OUTCAR 目录（每个构型一个 `POSCAR-*/OUTCAR`）。Mg4C60 的 312 个
   MC-rattle 构型在 jzzn 的 `/public/home/.../joint_research_project/hyl/Mg4C60/new2/new2/new_mcrattle_structures`，
   每个 OUTCAR 只有 ~0.65 MB（**不要拷 CHGCAR/WAVECAR**）→ 一次拷到 3090 机再改 `DATA_DIR`。
4. **环境**：`CONDA_ENV=wc` 里要有 `ase`（读 OUTCAR 与建 xyz）、`matplotlib`（S4 画图）。

## 必须知道的物理边界（别只信力 RMSE）

MD 版换的是**估计量**（无 IFC 截断、含全阶与温度重正化），但**换不掉数据**：训练标签
仍然是那批位移 RMS ≈ 0.017 Å 的数据，长程三阶非谐在力上只有 ~1.5e-4 eV/Å（比残差
9.5e-4 小一个量级），**势对长程非谐的精度无法由这份数据约束**。所以：

- 力 RMSE ≤ 门槛（**优先相对门槛** `FORCE_RMSE_REL_MAX` = 力 RMSE / DFT 力 RMS，
  默认 5%；设 0 才回落到绝对门槛 `FORCE_RMSE_MAX` eV/Å）只是**必要**条件；
- 必须再做两项（当前**还没实现**，见下节 TODO）：
  1. **非谐一致性**：势做有限差分得到的 Φ³ 与 pheasy 拟合在**已确定区（≤3 Å）**的 Φ³ 对比
     （参考：mean|Φ³| = 47 / 21 / 0.9 eV/Å³ @ 0–1/1–2/2–3 Å，跨截断稳定 ≤1–20%）；
  2. **谐性一致性**：势的声子谱 vs 现成 fc2 声子谱（低频支/声速偏差 < 5%）。
- 只有这两条过了，MD 的 κ 才有资格和 IFC 路线的 κ 放在一张图上讨论。

## TODO / 已知未完成

- [ ] README 里说的 S3_valid 验证步（声子 vs fc2、Φ³ vs pheasy）——目前只做了力 RMSE 闸；
- [x] `hac.out`（EMD/Green-Kubo）列映射与文档逐列核对（col7-11，in+out 相加）；
- [x] 力 RMSE 相对门槛（`FORCE_RMSE_REL_MAX`，力 RMSE / DFT 力 RMS）；
- [x] 量子修正的谱分解版本（`compute_shc` → κ(ω) → f(ω) 逐频修正，已接入 gpumd_kappa.py 的 SHC_ENABLED）；
- [ ] 端到端 smoke test：用 `~/gpumd/potentials/C_2022_NEP4.txt` + 纯碳体系把 S3 跑通
      （只验流程/数值，不追物理值）；
- [ ] 2D 材料的真空方向与 κ 张量归一（本技能先按 3D 写）。
## NEP 侧的关键坑（2026-09-25 实跑踩到，务必先看）

`nep` 的报错默认被缓冲吞掉 —— **一律用 `stdbuf -o0 -e0` 跑**，否则只能看到 rc=1 而看不到原因。
本技能的 `nep_train.py` / `gpumd_kappa.py` 已经内置 stdbuf。依次踩到的四个坑：

1. `type` 需要 **≥2 个参数**：`type <类型数> <元素符号...>`（只给符号报 number of types should be integer，
   只给数字报 type should have at least 2 parameters）。
2. **微调时 `type` 必须与基座模型一致**：nep89 有 **89 种元素**，`type 2 Mg C` 会被拒；
   正确做法是从基座 `.txt` 首行读（`nep4 89 H He Li ...`）再原样写进 `nep.in` —— `nep_train.py` 已自动这么做。
3. **没有 train/test 文件参数**：`nep` 在工作目录里固定读 `train.xyz`（`test.xyz` 可选）；
   写 train_xyz 之类的关键字会报 invalid keyword。
4. `save_potential <间隔> <start> <end>`、`generation <N>`、`population`、`batch` 等照基座 `nep.in` 给即可。

## Si smoke test 结果（2026-09-25，已修 3 个真 bug）

用 `potentials/nep/Si_2022_NEP4_4body.txt` + 金刚石 Si 跑通了 S3 的 MD-κ 流程，
暴露出三个只有在真跑时才会显形的错误，都已修好：

1. **ensemble 语法**：v5.6 只认 `nvt_ber / nvt_nhc / nvt_bdp / nvt_lan / nvt_mttk / nvt_qtb / nve …`，
   **没有 `nvt_berendsen`**；且 `nvt_ber` 要 **3 个参数**：`ensemble nvt_ber <T_init> <T_final> <tau/dt>`。
   引擎现在写的是这个（`T_FINAL` 留空 = 恒温）。报错信息默认是缓冲的，看不到时用 `stdbuf -o0` 跑。
2. **初始结构文件名必须是 `model.xyz`**（不是 `xyz.in`）——否则报 "Failed to open model.xyz"。
3. **κ 张量的对角元要按驱动方向取**：`compute_hnemd` 输出的是相对驱动方向 μ 的 κ_μx/κ_μy/κ_μz，
   所以 x 驱动取 col0+col1、y 驱动取 col2+col3、z 驱动取 col4（第一版一律取 col0+col1，已修）。

环境相关：

- 3090 机上 `~/gpumd/src/gpumd` 是 **CUDA 版**（启动日志会列 GPU）→ 走 gpu 分区单卡；
  `submit_gpumd.tpl` 已按此写。CPU 版请用 step.conf 的 `[submit]` 段把 partition 改回 cpu192。
- 卡 4（`GPU0000:84:00.0`）当前有驱动错误（`nvidia-smi` 也认不到它），用 0–3 / 5。

**小胞不收敛（重要的物理判据）**：64 原子（L = 10.9 Å）跑 220 ps，得到
κ_xx = −995 ± 491、κ_yy = 672 ± 297 W/mK，首尾值相差几千 —— 完全在噪声里。
Si 的声子 MFP 长，必须 **≥ 512 原子 + ≥ 1 ns**，并做 Fe / cell / run-length 收敛才有意义。
> 对照值口径（重要）：文献 κ(Si, 300 K) ≈ 150 W/mK 是**实验值**（含量子效应 + 同位素散射），
> **不能当经典 MD 的标尺**。经典 MD 的 Si κ 应约 180–220 W/mK，正确对照是
> **同一个势**算出的 BTE κ、或已发表的同类势 MD 值。Si 这种高 κ 材料应直接上 **HNEMD**
> （`compute_hnemd`，收敛快一个量级），不要用 512 原子/1 ns 的 EMD（必然不收敛）。

