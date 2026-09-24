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

其中 **gamma_total = gamma_anh + gamma_isotope**：phono3py 计算 κ 时，碰撞矩阵
对角元用的就是这个总和（conductivity/base.py 的 _get_main_diagonal），
所以 hdf5 里的 κ 含同位素散射。只取 gamma_anh 会把寿命系统性算长
（Si 300 K 全模中位 +19%，光学支 +24%）。用 INCLUDE_ISOTOPE 控制（默认 true）。

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
    # 没有 kappa hdf5、但有 fc2/fc3：现跑 phono3py
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.RUN_PHONO3PY=true
    autozt -tt phonon-lifetime -p <材料> -j S1_lifetime conf --set params.PHONO3PY_MESH="7 7 7"

### 数据源解析顺序（SOURCE=auto）

1. KAPPA_HDF5 显式路径；
2. 技能/步骤目录里已有的 kappa-m*.hdf5 / gamma-m*.hdf5；
3. 材料目录树里上游技能（默认 kl-dft-cpu / step6_kappa）的 kappa hdf5；
4. 上游 fc2.hdf5 + fc3.hdf5 + phono3py_disp.yaml → RUN_PHONO3PY=true 时现跑
   phono3py-load ... --br；否则报错并给出三种解决指引。

## 4. 产出

| 文件 | 内容 |
|---|---|
| lifetime_summary.json | marker LIFETIME_DONE；各温度全/声学/光学分组的 tau 中位/均值/P10/P90（含同位素）、只看三声子的 tau_median_anharmonic_ps、omega*tau 中位、过阻尼比例；数据源、mesh、约定、isotope_included、boundary_mfp |
| lifetime_data.csv | 逐模：T,gp,q,band,f,gamma_anh,gamma_iso,gamma_total,tau,tau_anh,omega_tau,period,逐模 v,mfp,acoustic |
| lifetime_T.csv | 逐温度汇总一行（含 tau_median_all_anharm_ps） |
| lifetime_vs_frequency.png | tau vs 频率（多温度） |
| lifetime_vs_period.png | tau vs 振动周期 1/f，含 tau=T_period、omega*tau=1 参考线 |
| ioffe_regel.png | omega*tau vs 频率 |
| lifetime_vs_T.png | 中位寿命 vs 温度 |
| lifetime_vs_q.png | （可选，设 Q_DIRECTION）沿该方向的逐支寿命 |

## 5. 验证（2026-09-24，含修正后复测）

真实数据自测（Si，上游 kl-dft-cpu 的 kappa-m151515.hdf5，15³ 网格、4×4×4 超胞、ALM）：

* **寿命因子用 phono3py 自己的 κ 反推验证**：用
  mode_kappa = gv_by_gv * cv / (2*(gamma+gamma_iso)) * (get_unit_to_WmK()/V)
  重建，与 hdf5 里的 kappa/mode_kappa **逐模吻合到 1e-6（99.3% 的模）**，
  且 kappa=Σmode_kappa/N 精确成立 → 证明 tau=1/(4*pi*gamma_total) 就是 phono3py 的约定；
* 群速度量级自洽：近 Gamma 处 LA 的 |v|≈86 THz·Å ≈ 8600 m/s，与 Si 纵声速一致；
* 300 K（含同位素）全模中位 tau=3.53 ps、声学 12.9 ps、光学 2.17 ps；
  不含同位素则 4.35/13.9/2.87 ps（同位素使光学支 −24%，与 ω² 标度一致）；
* 随温度单调下降（100 K 全模中位 9.88 ps）；
* 同一数据源 κ(300 K)=88.8 W/mK；
* RUN_PHONO3PY=true 从 fc2/fc3 现跑 3³ 网格也走通（自动生成 kappa-m333.hdf5）。

回归测试：tests/test_phonon_lifetime.py（12 条，含约定比对、同位素、gen 端到端）：

    python -m pytest tests/test_phonon_lifetime.py -q

## 6. 修正记录（2026-09-24，首版测试后）

* **漏了同位素散射**：首版只用 gamma，寿命偏长。现默认用 gamma+gamma_isotope
  （与 phono3py κ 一致），并在 summary/CSV 同时保留纯三声子值，加开关 INCLUDE_ISOTOPE。
* **逐模 mfp 用错群速度**：首版对 q 点内 6 条支的 |v| 取平均再乘各支 tau，
  导致 mfp/v 列错误。现改为逐模自己的 |v|。

## 7. 局限

1. 只做**三声子（最低阶微扰 / RTA）**；四声子、电子-声子寿命不在内。
2. 边界散射：hdf5 的 boundary_mfp 若为有限值，phono3py 的 κ 会含边界散射；
   本技能默认只叠加 gamma_isotope（Si 这组 boundary_mfp=1e6 Å，可忽略），
   未把边界散射折进 tau。
3. 声学支用“最低 3 条支 + 频率窗口”近似；交叉点处与“真声学支”会有少量错分。
4. ACOUSTIC_FMAX 留空时自动取“最低 3 支正频模中位频率 ×2”，2D 材料建议显式给值。
5. 不含 MD 声子电流关联函数（dynasor）那条有限温度路线。
6. RUN_PHONO3PY=true 在登录节点同步跑，大胞/密网格请改用上游 kl 的 S6 产物。

## 8. 文件

    skill/phonon-lifetime/
    ├── skill.yaml                 # 技能定义（1 步 run: gen）
    ├── templates/step.conf        # 全部参数默认值（shared 布局）
    ├── gen_step1_lifetime.py      # 驱动：定位数据源 -> 读 gamma(+iso) -> tau -> CSV/图/JSON
    ├── lifetime_common.py         # 物理约定 + 读取 + 汇总（无 matplotlib 依赖）
    ├── lifetime_plots.py          # 出图（无 matplotlib 时自动跳过）
    └── README.md
    tests/test_phonon_lifetime.py  # 12 条回归（约定 + 同位素 + 端到端）
