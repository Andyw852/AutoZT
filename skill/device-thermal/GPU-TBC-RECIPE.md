# GPU 扩展：GPUMD NEMD 界面热导（通用版，已在 3090 跑通）

本技能可选步 **S5_tbc / S5b_tbc_run / S6_tbc_post**（optional_steps.gpu_tbc，默认关）在 GPU 上用
GPUMD NEMD 算界面热导 G，再用 TBC_OVERRIDE 回灌器件模型（S2_props）。

## 0. 硬件与软件现状（实测）

- **3090（本步默认 hpc）**：GPUMD v5.6 已编译，二进制
  `~/software/AutoZT/GPUMD/build/gpumd`（cmake+CUDA 12.4, sm_86）。
  通用势 NEP89：`~/software/AutoZT/GPUMD/potentials/nep/nep89_20250409/nep89_20250409.txt`。
- **a800**：同样装有 GPUMD 与 NEP89。
- ⚠️ **3090 的 GPU4 处于 Unknown Error**（会拖累整机 `nvidia-smi` 全量枚举），
  GPU0-3 常被他人占用，GPU5 空闲。**必须用 UUID 指定 GPU**（数字索引会因故障卡错位）：
  `CUDA_VISIBLE_DEVICES=GPU-377470de-1983-1db4-5a7b-c8d04f5157ac`（= nvidia-smi 的 GPU5）。
  数字 `CUDA_VISIBLE_DEVICES=5` 会报 "no CUDA-capable device"。
  单卡查状态用 `nvidia-smi -i 5 ...`（全量 `nvidia-smi` 会被故障卡打断）。

## 1. 通用性（2026-09-23 修复）

本次把 S5/S6 从"只认 Si/SiO2 + 文献 kappa 反推热流"改成通用实现：

| 维度 | 支持 |
|---|---|
| 结构 | 任意 extxyz：自动解析 Lattice / Properties 列布局 / pbc |
| 元素 | 任意；若是 NEP 势会自动校验结构元素都在势的元素表里 |
| 势 | 任意 GPUMD 势文件（NEP / Tersoff / ...） |
| 传输轴 | x / y / z（`TBC_AXIS`） |
| 界面判定 | 按**元素组分沿轴的突变**自动定界面（`TBC_INTERFACE_COORD` 可显式覆盖） |
| 热流 | 由源/漏恒温器**累积传能**直接给出 -> 不再需要文献 kappa |
| 源/漏位置 | `TBC_SOURCE_SIDE=high\|low` 可选，厚度可调 |
| 可复现 | `TBC_SEED` 非 0 时固定初速随机种子 |

## 2. 用法

1. 打开可选组（改项目配置，需先请示）：在项目配置 `project_setting/tf_<材料>_device-thermal.yaml`
   的 `task_types.device-thermal:` 下加一行 `gpu_tbc: true`（**类型顶层**，不是 `optional:` 子段）。
   ⚠️ **不要用** `autozt -tt device-thermal -p <材料> conf --set optional.gpu_tbc=true`：
   `conf --set` 只把键写进**本步 step.conf**（变成 `[optional] gpu_tbc = true`），而可选组开关读的是
   类型顶层的 `t["gpu_tbc"]`（`autozt/bootstrap.py` 的 `expand_optional_steps`），写 step.conf
   **不会**展开 S5/S5b/S6（2026-09-26 实测确认：加 `gpu_tbc: true` 后状态表从 4 列变 7 列，
   而 `conf --set optional.gpu_tbc=true` 后仍是 4 列）。
2. 放界面结构：把 extxyz 放到该材料的 `project_setting/skill_dir/`（随 project_setting
   同步到远端），再把 `TBC_STRUCTURE` 写成该文件名；或直接把 `TBC_STRUCTURE` 写成 GPU
   主机上的绝对路径。**结构文件不在 `gen_need` 里**（属用户数据，不进 per-gen 推送）；
   缺文件时 S5 会报“找不到界面结构”，而不是含糊的 gen_need 报错。
3. 设 `TBC_NEP_MODEL / TBC_GPUMD_BIN / TBC_CUDA_VISIBLE_DEVICES`，必要时设
   `TBC_AXIS / TBC_INTERFACE_COORD / TBC_SOURCE|SINK_THICKNESS / TBC_*_STEPS`。
4. 生成输入：S5_tbc -> `step5_tbc/{model.xyz, potential.txt, run.in, run_gpumd.sh, tbc_inputs.json}`。
5. 运行：S5b_tbc_run 在 GPU 主机执行 `step5_tbc/run_gpumd.sh`（等价于手动 `bash step5_tbc/run_gpumd.sh`）。完成判据 = `compute_chunk.out` 块数 >= `TBC_RUN_STEPS // (TBC_SAMPLE_INTERVAL * TBC_OUTPUT_INTERVAL)`；达标才写 `step5b_run/run_done.json` 并放行 S6。
6. 提取 G：S6_tbc_post -> `step6_tbc_post/tbc_result.json`（G / dT_i / q / quality）。
7. 回灌：`autozt -tt device-thermal -p <材料> -j S2_props conf --set params.TBC_OVERRIDE=<G>` -> start。

## 3. 生成的 run.in（三段：平衡 -> 烧入 -> 测量）

```
potential potential.txt
velocity 300 seed 20260923        # TBC_SEED 默认 20260923；设 0 则不写 seed（不可复现）
time_step 1
ensemble nvt_nhc 300 300 100      # 平衡
dump_thermo 1000
run 20000                         # TBC_EQUIL_STEPS
ensemble heat_lan 300 100 40 1 2  # 烧入：源=T+40 漏=T-40（TBC_THERMOSTAT/TBC_DELTA_T）
dump_thermo 1000
run 100000                        # TBC_BURN_STEPS（暂态，不测量）
ensemble heat_lan 300 100 40 1 2  # 测量段必须**重新**声明 ensemble
compute 0 10 100 temperature      # 输出各组 T + 源/漏恒温器累积传能
compute_chunk 10 100 bin/1d z lower 2 temperature density/number
dump_thermo 1000
run 200000                        # TBC_RUN_STEPS（测量段）
```

- ★ **每个 `run` 前都要有自己的 `ensemble` 和 `compute*`**：
  - `Measure::finalize()`（每个 run 结束调用）会 `properties.clear()`：`compute` / `compute_chunk`
    只对紧接着的那一个 run 生效。只声明一次时 `compute_chunk.out` 只有它那一段（如烧入段）
    的剖面，与测量段的 q 不同窗 —— 实测本机 GPUMD 源码即如此。
  - 用户核对的 GPUMD **v5.8** 里 `Integrate::finalize()` 还会 `ensemble_.reset()`：测量段 run
    前没有 `ensemble` 会直接报 `An ensemble must be specified before each run.`（`run.cu`），
    后面根本跑不到。（本机 3090 的旧构建不报这个错，但不能依赖。）
  技能脚本已把 `ensemble + compute + compute_chunk` 都写进测量段（`tests/suite_device_thermal.py`
  的 `[s5-deck]` 回归测试守着这两点）。
- `ensemble heat_* <T> <T_coup> <delta_T> <src> <sink>`：源/漏温度分别为 T+delta_T / T-delta_T。
- `compute_chunk ... bin/1d <axis> lower <d> temperature`：逐 bin 温度 -> compute_chunk.out。
- model.xyz 第二行：`pbc="..." Lattice="..." Properties=species:S:1:pos:R:3:group:I:1`。
  - pbc 优先用 `TBC_PBC`，否则沿用结构文件的 pbc；**两者都没有时直接报错，不再默认 `T T T`**。
  - ★ 沿传输轴 `pbc=T` 时源(高端)/漏(低端)隔着周期边界直接相邻，热流不经界面短路、q 被高估
    （"源/漏热流自洽"判据抓不到它，因为能量两边仍守恒）。技能在 `pbc[传输轴]=T` 且界面非真空层
    时直接报错；正确做法：传输轴设 `F`，或结构两端留真空/改用对称布置。

## 4. 热流为什么不用 kappa（关键修复）

GPUMD 的 `compute <group> <s> <o> temperature` 会输出各组温度，**并在最后两列追加
源/漏恒温器的累积传能（eV）**（见 GPUMD `compute_out.rst`）。该累积在每个 `run` 段内独立计数，
所以"声明 compute -> 单独 run 测量段"就能拿到干净的 E_src(t)、E_sink(t)：

    q = |dE_thermostat/dt| / A       (A = 垂直于 TBC_AXIS 的横截面积)

旧版用 `q = -kappa * dT/dz`（kappa 取文献值），在引线偏短（弹道区）或未达稳态时
会差好几个量级 —— 这正是旧版 S6 把 quality 判成 poor 的根因。新版直接用恒温器能量，
并**用源/漏两侧热流是否相等作为稳态判据**。

## 5. 后处理（gen_step6_tbc_post.py）

1. 先丢掉前 1/3 暂态，再在剩余后缀里挑**最长**的源/漏热流自洽窗口（`TBC_STEADY_TOL` 容差、
   `TBC_MIN_WINDOW_FRAC` 最短占比）；温度剖面只平均窗口内的 block，保证 T(z) 与 q 同窗。
   dT_i 的统计误差用**时间分块 + 自相关修正**：compute_chunk 每块约 1 ps，而温度涨落
   的相关时间有几十 ps，相邻块不独立，直接 std/sqrt(N) 会低估（相关时间 20 块时约 8 倍）。
   现同时用积分自相关时间（N_eff=N/2tau）与 8~10 个超级块的批均值，取较大者（保守）；
   不再把空间上高度相关的相邻 bin 当独立样本。可用时间块少于 3 个时误差记 None，
   调用方必须据此判 `poor`，而不是跳过 `TBC_DT_TOL` 判据。
2. 对界面两侧体区（各剔除源/漏恒温器区与 `TBC_FIT_MARGIN_A` 余量）线性拟合 T(z)。
3. dT_i = 两侧外推到界面坐标的温差；**G = q / dT_i**。
4. 输出 `quality`：`good` 要求源/漏热流自洽、两侧 bin 数够、dT_i>0；否则 `poor`。

## 6. 质量判据与局限（G 不要当定量结论）

S6 输出 `quality`（good/poor）+ `quality_reasons`，判据：

- `q_consistency`（源/漏热流相对差）≤ `TBC_STEADY_TOL`（默认 0.2）
- 两侧 bin 数 ≥ `TBC_MIN_BINS`，且 `dT_i > 0`
- **dT_i 相对误差** ≤ `TBC_DT_TOL`（默认 0.3）——由稳态窗口内各时间 block 的 dT_i
  做**自相关修正**（N_eff=N/2tau）+ 超级块批均值给出，取较大者（时间块不足时误差为
  None，直接判 poor）
- **G 对拟合余量的敏感度** max/min ≤ `TBC_MARGIN_TOL`（默认 1.5）——对 1/2/4 A 余量各算一次

局限：

1. **引线长度**：体区只有几 nm 时处于弹道区，"体斜率"不是傅里叶梯度，两线外推失效；
   G 对拟合余量仍敏感（Si/Ge 三种子 439/443/475 MW/m2K，余量敏感度 max/min 1.37~1.43）。
   拿到收敛 G 必须**加长引线**到出现真正线性中段，并多种子平均。
2. **界面位置约定**：自动判定取界面键中心；改变位置会让 dT_i 变化 ~(slope_A-slope_B)*dz，
   需要严格比较时显式给 `TBC_INTERFACE_COORD`。
3. **共格、无混合界面 = TBC 上限**；真实界面有互混/缺陷会更低。
4. **通用势未针对界面验证**（NEP89 是 89 元素通用势），定量前建议用专用势或做收敛测试。
5. `kappa_crosscheck` 仅在扩散极限有意义，短引线偏离 1 属预期，不参与 quality。

## 7. 3090 端到端验证（2026-09-26，新 deck 三种子）

- 体系：共格 Si/Ge(001)，1280 原子，22.18 x 22.18 x 55.4 A，源=顶 8 A(Si)，漏=底 8 A(Ge)；
  界面自动定在 z=27.03 A（键中心），源/漏各 192 原子，面积 491.9 A²。
- 流程：S5 生成 deck -> S5b 在 3090 GPU5（UUID `GPU-377470de-…`）跑 GPUMD
  320000 步（20000 平衡 + 100000 烧入 + 200000 测量）-> S6 出 G。三个种子
  `TBC_SEED`=20260923/24/25，远端目录 `/data/user_3090/AutoZT_tbc_v3/seed_<seed>`。

| seed | q (W/m²) | 自洽度 | dT_i (K) | G (MW/m²K) | tau_int | n_eff | 余量敏感 | quality |
|---|---|---|---|---|---|---|---|---|
| 20260923 | 1.516e10 | 1.1% | 31.9 ± 5.6 | 475 | 2.3 | 29.0 | 1.43 | good |
| 20260924 | 1.515e10 | 1.9% | 34.5 ± 2.8 | 439 | 0.9 | 72.3 | 1.37 | good |
| 20260925 | 1.526e10 | 1.5% | 34.5 ± 5.6 | 443 | 3.6 | 18.5 | 1.41 | poor |

- q 三种子一致到 0.4%，G = 452 ± 20 MW/m²K（±4.4%），`n_blocks_averaged`=134
  （测量段 200 行丢前 1/3），符合预期。
- seed 20260925 判 poor 的唯一原因是 `tau_int`=3.6 块、超级块 16 < 5·tau_int=18.2
  （批均值仍相关）——新的误差可信度判据在起作用，不是链路故障。
- 结论：**链路、S5b 块数门控与质量/误差判据都工作正常；该短引线体系的 G 未收敛**，
  绝对 G 不可引用（引线 27 A 仍在弹道区，需加长引线 + 多种子平均）。
- ⚠️ **历史结论作废**：2026-09-23 那次 0.92 / 2.21 GW/m²K 用的是**旧 deck**
  （测量段 run 前没有重新声明 ensemble/compute，`compute_chunk` 只声明在烧入段）。
  `Measure::finalize()` 每个 run 后 `properties.clear()`（新旧版本一致），旧 deck 因此
  q 取自测量段、T(z) 取自烧入段，两者不同窗 —— 那两个 G 作废。
  现行脚本已把 `ensemble + compute + compute_chunk` 都写进测量段（见第 3 节）；
  2026-09-17 之后的 GPUMD（`17d06ad3`/`329d05ea`）旧 deck 会在测量段直接报错跑不到底，
  新 deck 在每个 run 前重新声明 ensemble，新旧版本都能跑。

## 8. 弹道 BTE（未安装）

- OpenBTE（https://github.com/romanodev/OpenBTE，arXiv:2106.02764）：CUDA 稳态声子 BTE，读 phono3py fc2/fc3。
- JAX-BTE（arXiv:2503.23657）：GPU 可微 BTE。两者都需在 GPU 主机安装，依赖 kl-dft-cpu 的力常数。

## 9. 一句话

**连续介质传热（S1-S4）= CPU 秒级；界面 TBC（S5/S6，GPUMD）= GPU（3090/a800）；弹道 BTE 仍需另装。**
