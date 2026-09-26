# phonon-lifetime —— 声子寿命 / 线宽（三声子，phono3py）

> 把上游 kl-dft-cpu（或 kl-mlff-*）已经算出来的 phono3py 声子自能 **gamma**
> 显式变成**声子寿命**：tau[ps] = 1/(2*2*pi*gamma_total[THz])，
> gamma_total = gamma_anh + gamma_isotope。
> 附带线宽、平均自由程、omega*tau（Ioffe-Regel）、寿命 vs 振动周期图。

## 1. 为什么需要它（不重复造轮子）

kl-dft-cpu 的 S6 跑 phono3py --br 得到 kappa-m*.hdf5，里面**同时存了**每模的
gamma（声子自能虚部/线宽）与 gamma_isotope，但技能只输出 κ 张量，
**没有把 gamma 变成寿命**。本技能就是那个下游分析层：

* **不重跑 DFT**：直接读上游已算好的 kappa-m*.hdf5（首选），
  或读 fc2.hdf5 / fc3.hdf5 / phono3py_disp.yaml 现跑一次 phono3py --br。
* **不碰 SLURM**：整个技能只有一步 run: gen（登录节点秒级出结果）。
* 与 device-thermal / device-zt 同属“消费上游产物”的分析技能，用 find_upstream
  同一套目录约定。

> 论文 Kim et al. Nature 597, 660 (2021) 里 **bulk MoS2** 的寿命来自 PHONO3PY
> （最低阶微扰），本技能算的就是这一类；论文里 r-MoS2 的**有限温度 MD 声子电流
> 关联函数拟合**（dynasor）是另一条路线，**不在本技能范围**。

## 2. 物理约定（关键，勿混淆）

phono3py 的 gamma 单位是 **THz（普通频率）**，内部按角频 2*pi*gamma 使用，因此

    tau[ps] = 1 / (2 * 2*pi * gamma_total[THz])  =  1 / (4*pi*gamma_total)

其中 **gamma_total = gamma_anh + gamma_isotope + gamma_bd**：phono3py 计算 κ 时，
碰撞矩阵对角元用的就是这个总和（conductivity/base.py 的 _get_main_diagonal），
所以 hdf5 里的 κ 含同位素散射。只取 gamma_anh 会把寿命系统性算长
（Si 300 K 全模中位 +19%，光学支 +24%）。用 INCLUDE_ISOTOPE 控制（默认 true）。
边界散射 gamma_bd = |v|/(4*pi*L) 只在**上游确实开了**边界散射时叠加——由 hdf5 里
κ 的重建误差自动判定（见 §5）；phono3py 默认的 boundary_mfp=1e6 **µm** 只是“关”的占位。
注意 L 的单位是微米（help 原文 "Boundary mean free path in micrometer"），换算式见 §6。

统计量一律按**不可约网格权重 weight** 加权（中位数/P10/P90/均值/过阻尼比例）：
kappa-m*.hdf5 里的 qpoint 只是不可约点，按点等权会让高对称点过度代表
（Si 15³ 只有 120 个不可约点、权重 1~48：等权中位 3.53 ps，加权 3.39 ps）。

这与 phono3py **官方后处理工具逐位一致**：

| 出处 | 写法 |
|---|---|
| phono3py/scripts/phono3py_kdeplot.py | tau = 1/g/(2*2*pi) |
| phono3py/other/kaccum.py | mfp = v/(2*2*pi*g)（mfp = v·tau） |
| phono3py/conductivity/base.py | get_unit_to_WmK() 的 1/(2*pi)（注释 “2pi comes from definition of lifetime”）+ RTA 的 /(g*2) + 边界散射 g=|v|/(4*pi*L) |

注意：直接写 tau = 1/(2*gamma) 会**偏大 2π 倍**（约 6.28×）。
本技能把这个因子固化成 lifetime_common.TAU_FACTOR = 2*2*pi，并写进
lifetime_summary.json 的 tau_convention 字段。

论文图 2e 的对照量：

* 振动周期 T_period[ps] = 1/f[THz]；
* omega*tau = 2*pi*f*tau；omega*tau < 1 即 Ioffe-Regel 判据下的强散射/过阻尼。

## 3. 用法

    # 1) 确认上游 kl 已经出 kappa-m*.hdf5（材料目录下，通常自动被 find_upstream 找到）
    autozt -tt kl-dft-cpu -p <材料> status

    # 2) 只看输入怎么生成（不跑）
    autozt -tt phonon-lifetime -p <材料> init
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime init

    # 3) 正式生成/运行（run: gen，登录节点）
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime start

常用参数（改参数不要改脚本）：

    # 只看 300/600 K
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.TEMPERATURES="300 600"
    # 显式指定 hdf5
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.KAPPA_HDF5=/path/kappa-m151515.hdf5
    # 沿 Gamma-A 出 tau-q 图
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.Q_DIRECTION="0 0 1"
    # 只看三声子（不含同位素）寿命
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.INCLUDE_ISOTOPE=false
    # 没有 kappa hdf5、但有 fc2/fc3：现跑 phono3py（自动带 --isotope；有 BORN 自动启用 NAC）
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.RUN_PHONO3PY=true
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.PHONO3PY_MESH="7 7 7"
    # 逐模 CSV 默认只写 FOCUS_T；要全温度再打开（大胞密网格会写出上千万行）
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.CSV_ALL_TEMPERATURES=true

### 数据源解析顺序（SOURCE=auto）

1. KAPPA_HDF5 显式路径；
2. 技能/步骤目录里已有的 kappa-m*.hdf5 / gamma-m*.hdf5；
3. 材料目录树里上游技能（默认 kl-dft-cpu / step6_kappa；kl-mlff-* 自动改试 step4_kappa，
   对不上再 rglob 兜底）的 kappa hdf5；
4. 上游 fc2.hdf5 + fc3.hdf5 + phono3py_disp.yaml → RUN_PHONO3PY=true 时现跑
   phono3py-load ... --br；否则报错并给出三种解决指引。

## 4. 产出

| 文件 | 内容 |
|---|---|
| lifetime_summary.json | marker LIFETIME_DONE；各温度全/声学/光学分组的 tau 中位/均值/P10/P90（含同位素）、只看三声子的 tau_median_anharmonic_ps、omega*tau 中位、过阻尼比例；数据源、mesh、约定、isotope_included、boundary_mfp_um / boundary_gamma_included / boundary_decision |
| lifetime_data.csv | 逐模：T,gp,q,band,f,gamma_anh,gamma_iso,gamma_total,tau,tau_anh,omega_tau,period,逐模 v,mfp,acoustic |
| lifetime_T.csv | 逐温度汇总一行（含 tau_median_all_anharm_ps） |
| lifetime_vs_frequency.png | tau vs 频率（多温度） |
| lifetime_vs_period.png | tau vs 振动周期 1/f，含 tau=T_period、omega*tau=1 参考线 |
| ioffe_regel.png | omega*tau vs 频率 |
| lifetime_vs_T.png | 中位寿命 vs 温度 |
| lifetime_vs_q.png | （可选，设 Q_DIRECTION）沿该方向的逐支寿命 |

## 5. 验证（2026-09-25，修 B1–B9 与 C1–C5 后重测）

真实数据自测（Si，上游 kl-dft-cpu 的 kappa-m151515.hdf5，15³ 网格、4×4×4 超胞、ALM）：

* **τ 因子与 κ 自洽（自检）**：按 mode_kappa = gv_by_gv * cv / (2*gamma_total) *
  kappa_unit_conversion、kappa = Σmode_kappa / Πmesh 重建，与 hdf5 里存的 kappa
  **最大相对误差 7.3e-6**（8 个温度）。这条自检写进 summary 的 kappa_reconstruction；
  口径不对（例如上游开了边界散射而我们没计入）会立刻报出来。
* **边界散射自动判定（单位微米，已修）**：默认 boundary_mfp=1e6 µm 直接判定“关”
  （phono3py 自己的日志也打印 Boundary mean free path (millimeter): 1000.000）。
  上游真开边界散射时才用 κ 重建误差判定：实测 '--bmfp 0.1'（100 nm）时“不计入”误差
  1.91e-01、“计入”**2.02e-16**（γ_bd 中位 5.15e-4 THz），声学 τ_median 14.4 → 13.2 ps。
* 300 K（含同位素、按权重）：全模中位 tau=**3.39 ps**、声学 **12.2 ps**、光学 **2.31 ps**，
  P10–P90 = 1.64–20.0 ps；纯三声子中位 4.22 ps。
* τ 近似 1/T：100→800 K 全模加权中位 9.88 / 4.97 / 3.39 / 2.64 / 2.19 / 1.89 / 1.63 / 1.46 ps。
* 声学支 = 最低 3 条支（不加频率窗）：Si 的 LA 到 11.76 THz 仍算声学，不再被旧版
  自动算出的 8.5 THz 频率窗误判进光学组。
* 沿 [001]：用 spglib 对称操作把不可约楔子 (x,0,0) 映射回目标线，取到 Γ→X **8 个点**
  （不给对称操作只剩 Γ ——旧版出的就是这种空图）。
* 同一数据源 κ(300 K)=88.8 W/mK；RUN_PHONO3PY=true 从 fc2/fc3 现跑 3³ 也走通，
  命令带 --isotope，κ 重建误差 2.0e-8（fc 目录有 BORN 时自动拷入启用 NAC）。

回归测试：tests/test_phonon_lifetime.py（**38 条**：phono3py 约定比对、加权统计、
选主文件、对称方向、边界散射判定、fc 命令构造、gen 端到端）：

    python -m pytest tests/test_phonon_lifetime.py -q

## 6. 修正记录

### 2026-09-26（0.4：2D κ 自动归一化）

* summary 的 `kappa_from_hdf5` 旧版直接吐 hdf5 的 raw κ（分母是含真空的胞体积），
  2D 材料会被当成文献值误用（MoS2：raw 27.9 vs 层内 103.9 W/mK）。
* 现在自动判定 2D（上游 summary 的 `dim`，或 q 网格恰有一轴为 1；`DIM` 可强制），
  输出 `kappa_2d_inplane_avg`（及逐温度、xx/yy 分量）= κ_raw × h⊥/d，并注明 d 与来源。
* 因子来源：优先上游 `kappa_summary.json` 的 `kappa_2d_norm_factor`（与 kl 报告严格同源）；
  读不到才用 `skill/_common/thickness_2d.slab_geometry`（kl/ke 共用的唯一真源）从
  phono3py yaml 的 primitive_cell 现算（`KAPPA_2D_THICKNESS`：vdw | 固定 Å | cell）；
  两路都有时互相对照，差 >1% 告警。两路都没有则只报 raw 并告警“不能直接对文献”。
* 寿命 tau 与厚度无关，本改动不影响任何寿命数值。
* 测试里 Si 参考数据路径改为读环境变量 `AUTOZT_SI_REF`（不再写死本机路径）。

### 2026-09-25（C1–C5；第二轮外部审计逐条对照后修）

* **C1 边界散射单位错了 1e4 倍（严重，B8 实质没修好）**：phono3py 的
  --boundary-mfp/--bmfp 与 hdf5 里的 boundary_mfp 单位都是**微米**
  （help 原文 "Boundary mean free path in micrometer"；默认 1e6 µm），而群速度是
  THz·Å，所以正确的是 γ_bd = |v|/(4π·L[µm]·1e4)。旧代码按 Å 算，γ_bd 偏大 1e4 倍：
  上游真开 '--bmfp 0.1'（100 nm）时会算出 ~5 THz 的荒唐线宽，计入后 κ 重建误差反而
  变大，“自动判定”就选了“不计入”——边界散射被静默丢掉，还被自动判定掩盖。
  现在：L >= 1e6 µm 直接判定“关”（phono3py 默认占位），真开了才用 κ 重建误差判。
  真跑验证（Si 3³、'--bmfp 0.1'）：不计入误差 1.91e-01、计入 2.02e-16。
* **C2 密网格混入离线点**：select_direction_qpoints 容差原是固定 0.06（分数坐标），
  31³ 格点间距 1/31≈0.032 比它还小，紧贴目标线的离线点会被当成线上点（实测整列取错）。
  现容差默认 1e-5，且同一个 t 只保留离目标线**最近**的候选。
* **C3 缺时间反演**：等价点原来只靠空间群旋转，没有 q→−q。没有反演中心的材料，
  不可约楔里若只存 −q 就会漏点（实测只剩 Γ）；现把 −imgs 也放进候选（τ(q)=τ(−q)）。
* **C4 接 kl-mlff 找不到上游**：UPSTREAM_STEP 写死 step6_kappa，而 kl-mlff-* 的 κ 在
  step4_kappa、fc2/fc3 在 step3_fc。现按 UPSTREAM_SKILL 自动展开候选步骤名，
  对不上再在上游技能目录下 rglob 兜底。
* **C5a 展宽法文件**：kappa-m*-s<σ>.hdf5（smearing）与四面体法主文件同分时先后不定；
  现给 -s 文件降级，并在选中它时告警。
* 仍未做：除 Si 外没有真跑过别的材料；过阻尼比例需要强非谐体系
  （MnIn2Se4 / 225 相）验证一次，集群上还没跑过。

### 2026-09-25（B1–B9；首轮外部审计用合成 hdf5 复现后修）

* **B1 按字符串序选文件**：find_kappa_hdf5 原用 sorted()[-1]，κ-m999 会盖过 κ-m151515。
  现按“多 q 点主文件 → mesh 乘积 → q 点数 → 文件名惩罚(-g0/mfp/gp)”排序，
  多候选时打印选了谁、忽略了谁；跨目录也优先主文件。
* **B2 本地目录取排序第一个**：改用同一套选主文件逻辑，新增 is_main_hdf5；
  本地只有 -g0 单点文件时不再挡掉上游主文件（都没有才退回单点文件）。
* **B3 空分组崩溃**：lifetime_T.csv / 出图遇 None 写空串（_fmt），不再 TypeError。
* **B4 统计量没乘权重**：中位数/P10/P90/均值/过阻尼比例全部按不可约 weight 加权；
  整数权重时展开到全网格再取分位数（与“展开后中位数”严格一致）；新增
  n_modes_weighted / weights_applied。
* **B5 声学支划分**：默认只取最低 3 条支、**不加频率窗**；只有显式给 ACOUSTIC_FMAX
  才叠加频率窗。
* **B6 沿 q 方向取点**：用 spglib 从 phono3py.yaml 的 primitive_cell 求对称操作，
  把不可约点的等价像映射到目标线；缺 yaml/pyyaml/spglib 时降级为几何取点并告警。
  同时修了“按量化格点取点会漏掉非格点 t”的问题。
* **B7 fc 路线缺同位素/NAC**：INCLUDE_ISOTOPE=true 时现跑命令加 --isotope；
  fc 目录有 BORN 时自动拷入（phono3py-load 没有 --nac 选项，靠 BORN 自动开 NAC）。
* **B8 边界散射**：新增 gamma_bd = |v|/(4*pi*L)，是否计入由 κ 重建误差自动判定；
  summary 记录 boundary_mfp_um / boundary_gamma_included / boundary_decision
  （单位与判定在 C1 修正，见上）。
* **B9 温度静默就近**：请求温度或 FOCUS_T 不在列表时告警，不再悄悄取最近值。
* 小项：src_kind 不再把本地文件标成 upstream；上游缺 gamma_isotope 时告警；
  逐模 CSV 默认只写 FOCUS_T（CSV_ALL_TEMPERATURES=false）。

### 2026-09-24（首版测试后）

* **漏了同位素散射**：首版只用 gamma，寿命偏长。现默认 gamma+gamma_isotope。
* **逐模 mfp 用错群速度**：首版对 q 点内 6 条支的 |v| 取平均再乘各支 tau，
  现改为逐模自己的 |v|。

## 7. 局限

1. 只做**三声子（最低阶微扰 / RTA）**；四声子、电子-声子寿命不在内。
2. 边界散射只在“上游确实开了”时计入：默认 1e6 µm 直接判为关，其余由 κ 重建误差判定；
   RTA 也不是全 BTE。L 的单位是微米，换算见 §6 的 C1。
3. 声学/光学按**支编号**分组（默认最低 3 条支），因此两组频率范围在支交叉处会重叠
   （Si 声学支最高 11.76 THz 略高于光学支最低 10.42 THz，属正常，不是错分）。
4. 沿 q 图只取**精确落在**目标线上的网格点（容差 1e-5，分数坐标）；方向与网格点不共线
   时只剩 Γ。对称操作依赖 phono3py.yaml / pyyaml / spglib，缺了会降级成几何取点并告警。
5. 同一 mesh 若同时有四面体法与展宽法文件，优先四面体法（-s 降级）并告警；本技能不重算，
   只是选文件。
6. 不含 MD 声子电流关联函数（dynasor）那条有限温度路线。
7. RUN_PHONO3PY=true 在登录节点同步跑，大胞/密网格请改用上游 kl 的 S6 产物。
8. 除 Si 外还没有真跑过别的材料；过阻尼比例（omega*tau<1）只在合成数据上验证过，
   真实强非谐体系（MnIn2Se4 / 225 相）与集群环境都还没跑。

## 8. 文件

    skill/phonon-lifetime/
    ├── skill.yaml                 # 技能定义（1 步 run: gen）
    ├── templates/step.conf        # 全部参数默认值（shared 布局）
    ├── gen_step1_lifetime.py      # 驱动：定位数据源 -> 读 gamma(+iso,+bd) -> tau -> CSV/图/JSON
    ├── lifetime_common.py         # 物理约定 + 读取 + 选文件 + 对称方向 + 加权统计（无 matplotlib）
    ├── lifetime_plots.py          # 出图（无 matplotlib 时自动跳过）
    └── README.md
    tests/test_phonon_lifetime.py  # 38 条回归（约定 + 加权 + 选文件 + 方向 + κ 自检 + 端到端）
