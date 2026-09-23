# GPU 扩展：GPUMD NEMD 界面热导（通用版，已在 3090 跑通）

本技能可选步 **S5_tbc / S6_tbc_post**（optional_steps.gpu_tbc，默认关）在 GPU 上用
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

1. 打开可选组（改项目配置，需先请示）：
   `autozt -tt device-thermal -p <材料> conf --set optional.gpu_tbc=true`
2. 放界面结构：把 extxyz 命名为 `tbc_interface.xyz`，放在**技能目录**
   `skill/device-thermal/` 或该材料的 `project_setting/skill_dir/`（S5 的 gen_need 会推送到
   GPU 主机）；也可把 `TBC_STRUCTURE` 直接写成 GPU 主机上的绝对路径。
3. 设 `TBC_NEP_MODEL / TBC_GPUMD_BIN / TBC_CUDA_VISIBLE_DEVICES`，必要时设
   `TBC_AXIS / TBC_INTERFACE_COORD / TBC_SOURCE|SINK_THICKNESS / TBC_*_STEPS`。
4. 生成输入：S5_tbc -> `step5_tbc/{model.xyz, potential.txt, run.in, run_gpumd.sh, tbc_inputs.json}`。
5. 运行：`bash step5_tbc/run_gpumd.sh`（已 pin CUDA_VISIBLE_DEVICES）-> compute.out / compute_chunk.out / thermo.out。
6. 提取 G：S6_tbc_post -> `step6_tbc_post/tbc_result.json`（G / dT_i / q / quality）。
7. 回灌：`autozt -tt device-thermal -p <材料> -j S2_props conf --set params.TBC_OVERRIDE=<G>` -> start。

## 3. 生成的 run.in（三段：平衡 -> 烧入 -> 测量）

```
potential potential.txt
velocity 300 seed 20260923        # TBC_SEED=0 时不写 seed
time_step 1
ensemble nvt_nhc 300 300 100      # 平衡
dump_thermo 1000
run 20000                         # TBC_EQUIL_STEPS
ensemble heat_lan 300 100 40 1 2  # 源=T+40 漏=T-40（TBC_THERMOSTAT/TBC_DELTA_T）
compute_chunk 10 100 bin/1d z lower 2 temperature density/number
dump_thermo 1000
run 100000                        # TBC_BURN_STEPS（烧入，释放初始应变，强烈建议 >0）
compute 0 10 100 temperature      # 输出各组 T + 源/漏恒温器累积传能
dump_thermo 1000
run 200000                        # TBC_RUN_STEPS（测量段）
```

- `ensemble heat_* <T> <T_coup> <delta_T> <src> <sink>`：源/漏温度分别为 T+delta_T / T-delta_T。
- `compute_chunk ... bin/1d <axis> lower <d> temperature`：逐 bin 温度 -> compute_chunk.out。
- model.xyz 第二行：`pbc="..." Lattice="..." Properties=species:S:1:pos:R:3:group:I:1`（pbc 沿用结构文件，可用 TBC_PBC 覆盖）。

## 4. 热流为什么不用 kappa（关键修复）

GPUMD 的 `compute <group> <s> <o> temperature` 会输出各组温度，**并在最后两列追加
源/漏恒温器的累积传能（eV）**（见 GPUMD `compute_out.rst`）。该累积在每个 `run` 段内独立计数，
所以"声明 compute -> 单独 run 测量段"就能拿到干净的 E_src(t)、E_sink(t)：

    q = |dE_thermostat/dt| / A       (A = 垂直于 TBC_AXIS 的横截面积)

旧版用 `q = -kappa * dT/dz`（kappa 取文献值），在引线偏短（弹道区）或未达稳态时
会差好几个量级 —— 这正是旧版 S6 把 quality 判成 poor 的根因。新版直接用恒温器能量，
并**用源/漏两侧热流是否相等作为稳态判据**。

## 5. 后处理（gen_step6_tbc_post.py）

1. 自动挑选**最长的**源/漏热流自洽窗口（后缀搜索，`TBC_STEADY_TOL` 容差、`TBC_MIN_WINDOW_FRAC` 最短占比）。
2. 对界面两侧体区（各剔除源/漏恒温器区与 `TBC_FIT_MARGIN_A` 余量）线性拟合 T(z)。
3. dT_i = 两侧外推到界面坐标的温差；**G = q / dT_i**。
4. 输出 `quality`：`good` 要求源/漏热流自洽、两侧 bin 数够、dT_i>0；否则 `poor`。

## 6. 质量判据与局限（G 不要当定量结论）

S6 输出 `quality`（good/poor）+ `quality_reasons`，判据：

- `q_consistency`（源/漏热流相对差）≤ `TBC_STEADY_TOL`（默认 0.2）
- 两侧 bin 数 ≥ `TBC_MIN_BINS`，且 `dT_i > 0`
- **dT_i 外推相对误差** ≤ `TBC_DT_TOL`（默认 0.3）——由两侧拟合的预测标准误合成
- **G 对拟合余量的敏感度** max/min ≤ `TBC_MARGIN_TOL`（默认 1.5）——对 1/2/4 A 余量各算一次

局限：

1. **引线长度**：体区只有几 nm 时处于弹道区，"体斜率"不是傅里叶梯度，两线外推失效；
   G 会随拟合余量显著变化（Si/Ge 例子里两次种子得 0.92 与 2.21 GW/m2K，都被判 poor）。
   拿到收敛 G 必须**加长引线**到出现真正线性中段，并多种子平均。
2. **界面位置约定**：自动判定取界面键中心；改变位置会让 dT_i 变化 ~(slope_A-slope_B)*dz，
   需要严格比较时显式给 `TBC_INTERFACE_COORD`。
3. **共格、无混合界面 = TBC 上限**；真实界面有互混/缺陷会更低。
4. **通用势未针对界面验证**（NEP89 是 89 元素通用势），定量前建议用专用势或做收敛测试。
5. `kappa_crosscheck` 仅在扩散极限有意义，短引线偏离 1 属预期，不参与 quality。

## 7. 3090 端到端验证（2026-09-23）

- 体系：共格 Si/Ge(001)，1280 原子，22.18 x 22.18 x 55.4 A，源=顶 8 A(Si)，漏=底 8 A(Ge)。
- 新版 S5 在 3090 生成 deck（自动定界面 z=27.03 A，源/漏各 192 原子）-> GPU5 跑 GPUMD
  （20000 平衡 + 100000 烧入 + 200000 测量 = 320000 步，485 s，exit 0）-> S6 出结果。
- q 两侧自洽 0.3%；但 dT_i 在两次独立种子间为 15.4 K 与 6.4 K（Si 引线非线性），
  G 分别为 924 与 2209 MW/m2K，都被新版判据判为 poor。
- 结论：**链路与判据都工作正常；该短引线体系的 G 未收敛**。细节见 tmp/device_thermal_e2e/。

## 8. 弹道 BTE（未安装）

- OpenBTE（https://github.com/romanodev/OpenBTE，arXiv:2106.02764）：CUDA 稳态声子 BTE，读 phono3py fc2/fc3。
- JAX-BTE（arXiv:2503.23657）：GPU 可微 BTE。两者都需在 GPU 主机安装，依赖 kl-dft-cpu 的力常数。

## 9. 一句话

**连续介质传热（S1-S4）= CPU 秒级；界面 TBC（S5/S6，GPUMD）= GPU（3090/a800）；弹道 BTE 仍需另装。**
