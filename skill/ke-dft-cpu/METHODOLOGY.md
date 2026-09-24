# ke-dft-cpu 方法学（2D 电子输运复现）

> 复现 J. Appl. Phys. 139, 024301 (2026) doi:10.1063/5.0308349 Table I 的**方法学口径**。
> 本文件只含方法学：口径定义、根因分析、验收判据、参数表。
> 六材料对比表 / 相对误差 / AMSET 旁证等**未发表结果**见数据目录（不进 git）：
> `/mnt/d/tf_data/jzz/jap/comparison_vs_literature.md`。

## 1. 有效厚度（厚度口径）

二维体系的 σ、κ_e、κ_p 都 ∝ 1/t，t 为有效层厚。本项目取
vdW 口径 d = zspan + r_vdW(top) + r_vdW(bot)（S 1.80、Se 1.90 Å），
单层自动得 CrS₂ 6.54、CrSe₂ 6.94 Å，与文献的 6.53/6.94 差 0.15%
——说明文献单层用的是同一套口径。

**ZT 与 t 无关。** `step8.1` 的 `_resolve_kappa_L` 只读 kl 链的原始
元胞口径 κ，再乘 ke 侧的 c/t，使 σ、κ_e、κ_p 共用同一个 t；
分子分母的 1/t 抵消。kl 链自己的 `KAPPA_2D_THICKNESS` 不进入 ZT，
只影响 kl 单独报出的二维归一 κ_p。

**但绝对值依赖 t。** 文献超晶格取两单层平均 6.73 Å，而 vdW 自动
值为 6.97 Å（由较厚的 CrSe₂ 侧决定），两者差 3.4%。要与文献的
κ_p 37/53（SS）逐点对比，需在两处分别锁 6.73：
- `ke-dft-cpu` 的 `LAYER_THICKNESS`（决定 σ、κ_e，以及 8.1 报出的 κ_p）
- `kl-*` 的 `KAPPA_2D_THICKNESS`（决定 kl 自己报的 κ_p）

两者互不影响，各自锁各自的；漏锁其一只影响该侧的绝对值，不影响 ZT。

## 2. E1 真空对齐（根因）

- ionrelax 两段式：`INCAR.relax IBRION=2` + `INCAR.static IBRION=-1 LVHAR=.TRUE.`，
  取 band_edges.json 的 `E1_vac_*`。
- 成因：amset 的 `get_reference_energy` 用原子核处的平均静电势，对 2D slab
  不跟随真空能级——这正是本工作用 LVHAR + 真空对齐绕开的那件事。amset 没有原生
  真空对齐选项，**重跑 step8 也无法解决**（不是 ICORELEVEL 那次能修的）。
- 验收：真空对齐线性拟合的 `off ≈ 0`（对齐斜率残差）。这步修掉 amset 芯势口径
  约 24% 的 E1 误差。

## 3. m* 口径

文献用的是 **full-BZ 二次型拟合的 m_d**（DOS 有效质量），不是 3 点抛物拟合的 m*。
用后者 μ 会低一半。

## 4. ε 真空稀释与二维散射核（amset2d）

step5_dielect 的 VASP 输出是 slab-in-a-box 介电（含真空稀释）。skill 的
`gen_step10_amset.py` 里 `_dielectric_2d_inplane` 已扣真空：
`ε_m = 1 + (c/t)(ε_slab - 1)`，c/t 复用弹性同款。

**但散射核本身还是三维的** —— 这一条不再靠「只作比值」糊过去：2D 项目改走
`step8.4_amset2d`（可选步骤组 `amset2d: true`），运行期插件把四种散射核换成二维形式。
原版 `step8_amset` 保留给三维项目，两者用同一批上游输入、在 DAG 上并行。

### 4.1 弹性张量顺序（三维也中招，已修）

VASP `TOTAL ELASTIC MODULI` 的行列标签序是 `XX YY ZZ XY YZ ZX`，AMSET 经 pymatgen
按标准 Voigt `XX YY ZZ YZ XZ XY` 解析。直接照抄会把第 4 位（XY）与第 6 位（ZX）互换：

- 2D slab：面内剪切 `C66 = C_xyxy` 被换成接近零的 `C_zxzx`。实测 CrS2 单层
  （step6_elastic/OUTCAR）：`C66 = 22.287 GPa` 被写到 `YZ` 槽位、`C_zxzx = -0.104 GPa`
  被当成面内剪切，面内 TA 支刚度近零。
- 3D 低对称（225 相 tetradymite 的 `C44 != C66`、`C14 != 0`）同样读错；立方晶系不受影响。

修复：`gen_step10_amset.py` / `gen_step14_amset2d.py` 的 `read_elastic` 读表头标签后
按标签重排（`_reorder_voigt`），行列同时重排；压电张量同理（`vasp_to_voigt_piezo`）。
**已经跑完的三维 AMSET 结果（如 225 相四个合金）需要评估是否重算。**

### 4.2 二维散射核的映射：G_s = c * G_2D

AMSET 实现的是三维费米黄金规则。真空 slab 的能带沿 `k_z` 平，能量守恒面是柱面，
`q_z` 的积分只贡献 `2*pi/c`：

```
AMSET 结果 = (2*pi/hbar) * Int d2q/(2*pi)^2 * (G_s/c) * delta(dE)
二维真值   = (2*pi/hbar) * Int d2q/(2*pi)^2 * G_2D     * delta(dE)
=> G_s = c * G_2D
```

插件把 `c*G_2D` 塞进 AMSET 的 `prefactor x factor`。四种机制的二维形式：

| 机制 | G_2D | 插件 factor（= c*G_2D / prefactor） |
|---|---|---|
| ADP | `k_B T*(q.v:D)^2/C_2D` | `sum_面内两支 (q.v_b:D)^2/w_b`，w 取**原始 slab** 弹性张量的面内 2x2 Christoffel 本征值 |
| POP | `2*pi*[1/eps_inf(q) - 1/eps_0(q)]/q` | `2*pi*c*(q.dr.q)/[(q.eps_inf(q)+q_TF)(q.eps_0(q)+q_TF)]` |
| IMP | `N_2D*(2*pi*Z*exp(-qd))^2/(q.eps_0(q)+q_TF)^2` | `c^2*(2*pi)^2*exp(-2qd)/(q.eps_0(q)+q_TF)^2` |
| PIE | `k_B T*(2*pi)^2*(q q:e_2D.v)^2*q^2/[C_2D*(q.eps_inf(q)+q_TF)^2]` | `c^2*(2*pi)^2*sum_b (q q:e.v_b)^2*q^2/(w_b*(q.eps_inf(q)+q_TF)^2)` |

其中张量 Keldysh 介电函数 `q.eps(q) = eps_env*|q| + q.r.q`，`r = c*(eps_slab,|| - 1)/2`
（2x2 张量，bohr）；`q_TF = 2*pi*c*dn_3D/dmu`（二维 Thomas-Fermi）。

**POP 的两个输入量 r∞ 与 Δr（2026-09-16 新增自检）**：由**同一次** DFPT 的原始 slab 张量给出
`r∞ = c*(eps_inf,|| - 1)/2`、`r0 = c*(eps_0,|| - 1)/2`、`Δr = r0 - r∞ = c*(eps_0,|| - eps_inf,||)/2`。
gen_step14 现在会算出这三个数、打印并落盘到 `2d_correction.json` 的 `r_inf_delta_r_check`，
插件每次运行也在横幅里复述一遍（`[amset2d] POP: r_inf=... delta_r=...`）。闸门：
`Δr < 0` 报警（离子介电缺失/读错）；`Δr/r∞ > 0.6` 报警（典型病因是拿**体相** ε0/ε∞ 外推
当 slab 值用——早期验证台子正是因此拿到 Δr = 43.5 Å，比 Sohier 的 0.53 Å 大 ~80 倍）。
文献锚点：MoS2 r∞ ≈ 46.5 Å / Δr ≈ 0.53 Å（Δr/r∞ ≈ 1.1%）；h-BN r∞ ≈ 7.6 Å / Δr ≈ 2.8 Å（≈37%）。

**ADP 的等价说法（最要紧的一条）**：`C_2D = c*C_slab`，所以 `G_s = k_B T*D^2/C_slab`，
即 settings.yaml 里的 `elastic_constant` 必须写**原始 slab 值**，**不能再乘 c/t**。
旧做法乘了 c/t -> ADP 散射率偏小 `t/c` 倍 -> ADP 迁移率偏大 `c/t` 倍
（CrS2：c/t = 3.06；实测旧流程 mu_ADP/正确值 = 3.1，见 VERIFICATION.md 的 T2）。

### 4.3 pop_frequency 取面内极性模

AMSET 自带的 `amset phonon-frequency` 对**三维球面**方向做 Born 电荷/本征矢加权，
slab 里会把面外极性模也算进去。二维 Fröhlich 只用面内极性模，所以本步改用
`amset2d_pop_freq.py`：同一套公式，只把 q 方向限制在 xy 面内。实测 CrS2：
三维口径 12.45 THz（含 2.4% 的 15.07 THz 面外模），面内口径 12.38 THz（纯 E' 简并模）。
插件再乘 `sqrt(eps_0(q)/eps_inf(q))`（二维 LST）得到 `w_LO(q)`。

### 4.4 形变势的构型口径（ionrelax）与参考能级（2D 默认真空）

**构型口径（2026-09-16 修订，影响三维）**：step7b 生成的**两份** h5（芯势口径的
`deformation.h5` 与真空口径的 `deformation_vac.h5`）都改为读 `deform-*/ionrelax/`
（离子弛豫构型，有则用）。物理依据：q→0 声学声子 = 均匀应变 + 内坐标弛豫；一致性依据：
ADP 核分母是 Christoffel 刚度（VASP `TOTAL ELASTIC MODULI`，含离子弛豫）。实测 MoS2 单层
xx 0.5%：电子 E1_vac 由 6.58（离子固定）变 8.51 eV（离子弛豫），差 **29%**。

> ⚠ **这是行为变化，不是纯修 bug**：HEAD 里三维 `step8/gen_step10` 读的正是
> `deformation.h5`，所以此改动之后**所有三维项目新算的 ADP 都会变**（含内坐标的体系
> 幅度类似）。三维数值与改动前**不可直接比较**；225 相合金重做时的差异要分成
> "缺陷修复"与"口径切换"两部分说明。
> ⚠ **回退不静默**：没有 `ionrelax/` 的构型退回离子固定口径时，gen 会打 WARN 并把
> `deform_geometry`（ionrelax / clamped / mixed）写进 `band_edges.json`，进而出现在
> settings 注释、`2d_correction.json` 与 8.3 汇总表里，跨项目比较前先核对这一项。

二维工作的通行做法是以**远离单层处的真空静电势**为零点定义带能，
C2DB 给出的形变势也是真空口径。两者在本项目实测差 1.66 倍（CrS2 单层：芯势 3.53 eV、
真空 5.86 eV），迁移率因此差 **2.76 倍**。

实现（2026-09-16 二次修订，按用户对 AMSET 0.4.19 源码的校对）：**替换参考能级**，不是事后加 Δ。
AMSET 的公式是 `diff = E_bulk − E_deform − (ref_bulk − ref_deform)` 再取 `np.abs(diff/strain)`，
所以 ① h5 里的 D 不带符号，"往它上面加带符号 Δ"只在 D_core>0 时才对；② `D_vac − D_core`
需要**两个**应变导数，只补 −dE_vac/dε 缺一项。现做法：patch 掉 parse_calculation 的
reference（换成各构型 LOCPOT 的真空能级），再调用 AMSET 原版的
`calculate_deformation_potentials`，符号/分量平均/噪声过滤与官方一致，另写
`deformation_vac.h5`（core 版保留作对比）；开关 `DEFORM_REF = auto|vacuum|core`，2D 默认 vacuum。

**两条运行期自检**（写进 gen_step9b，MoS2 实测都通过）：① 真空 h5 与 core h5 在带边必须**不同**
（相同 = 参考替换静默失效，直接抛错）；② `D_vac − D_core` **与能带无关**，两种载流子必须给
同一个值（MoS2：电子/空穴均 +2.0865 eV，相对差 0.00%；>5% 告警 —— 旧 clamped 数据在空穴上
不成立，因为空穴芯势口径的 D 实为负、被 AMSET 的 abs() 吃了符号）。

### 4.5 验证状态（详见 `step8.4_amset2d/VERIFICATION.md`）

| 项 | 结论 |
|---|---|
| 四种核的极限行为、旋转不变、护栏 | 合成数据单元测试全 PASS（amset 0.4.19 / 0.5.1） |
| Voigt 重排 | 真实 OUTCAR 实测；护栏能拦下未重排的输入 |
| "不该乘 c/t"（ADP 归一化的**比值**） | **通过**：旧/新 μ_ADP = 3.08，c/t = 3.059（比值测试，与 m*、D 无关） |
| **POP 绝对归一化** | **通过**：与独立 TB 参考（\`pop2d_ref_tb.py\`）逐点比，散射率中位比值 1.014、标准差 0.025（63 点全部 ±20% 内）；μ 2226 vs 2445（-9%） |
| **ADP 绝对归一化** | **通过**：同一 TB 台子上与独立解析参考比，μ 435.0 vs 444.3（**-2%**），逐能量散射率比值 1.03-1.2 |
| **IMP 绝对归一化** | **基本一致**：μ 188.8 vs 224.9（**-16%**），差的来源是屏蔽口径（参考用非简并 q_TF = 2πn/kT，插件用态密度积分出的费米能处 q_TF） |
| PIE 绝对归一化 | **未做定量验证**（默认关闭；只有单元测试覆盖极限行为），论文中如实写明 |
| r∞ / Δr 量级自检 | 已加入 gen（落盘 + 打印 + 闸门）与插件横幅；单元测试覆盖正常/零离子/缺张量三种输入 |
| 真空不变性（映射系数 c） | **可控台子上通过到 ~7%**：同色散、同步标定下 μ(c=40)/μ(c=20) = 1.071，而固定 c 只改网格的散布 ±5%。真实 slab 数据上仍是"分辨率不足"（插值外推噪声 ±30-45%） |
| 插值网格收敛 | **新查出的流程级问题**：interpolation_factor 6/10/15 给出 ~30% 差异；建议 2D 项目加收敛检查、DFT 侧 uniform 步 kz>=3 |
| 真实材料端到端 | 6 个 2D 项目跑通 5 个（CrS2/CrSe2 单层与正交胞、SS；LS 因 wavefunction 版本不符未通） |
| MoS2 本身 | 2026-09-16 起在 jzzn 跑 S1-S7（`/mnt/d/tf_data/MoS2`，uniform 步 kz>=3，S7 保留 LOCPOT 供真空参考能级用）；S8.4 待 S7 完成后在本机/集群单独跑 |

## 5. c 统一（比值法前提）

各材料元胞 c 轴统一 20 Å（有意为之），AMSET σ 的跨体系比值不含 c 差异伪影。
IMP/POP 扣真空、σ 归一化都依赖此前提。

**8.4（amset2d）之后**：二维核把真空长度显式积掉了，但**真空不变性尚未验证通过**
（V3：固定 c 只改插值网格就有 ×1.29-1.83 的噪声，大于 c 效应）——在解析参考测试与
kz>=3 的网格收敛测试完成前，仍建议用同 c（20 Å）的跨体系比值。

## 5.1 适用边界（按文献核对后修订，2026-09-16）

1. **"去掉 ZA 支"只对存在水平镜面（σh）的材料成立**：缺少 σh 的翘曲结构（硅烯、锗烯）
   与 Janus 结构可以有**一阶 ZA 耦合**，对它们是模型遗漏（ZA 是二次色散，形变势模型
   本身也处理不了）。
2. **缺非极性光学声子（光学形变势 ODP）**：AMSET 实现的机制只有 **ADP / POP / IMP / PIE**
   （`amset/scattering/elastic.py` = ADP+IMP+PIE，`inelastic.py` 只有 `PolarOpticalScattering`；
   官方文档 Summary of scattering rates 同）。**非极性晶体 POP 前因子为零**
   （`scattering/calculate.py:842` 注释明说），于是**整条光学声子通道缺失**：
   - **2D**：MoS2 的同极 A1 模与谷间光学模耦合很重要 -> **总**迁移率系统偏高。
     与文献对照（T4）必须**分机制**做：POP 对 Sohier 2016 的 LO 模散射时间；ADP 需要一个
     **有出处**的 LA 声子限制迁移率目标值（"340 cm2/Vs"这个数找不到出处，已撤回）。
   - **3D 非极性半导体（Si、Ge、金刚石这类）**：**空穴迁移率高估约 4~5 倍，电子几乎不受影响。**
     Si 300 K / 1e17 实测（AMSET 0.4.19 本地受控，见 `step8.4_amset2d/VERIFICATION.md` §六点一）：
     电子 overall **805**（实验 ≈800 ✔，比值 1.01）、ADP 1623（≈1400，1.16）；
     空穴 overall **1474**（实验 ≈330，**4.5 倍**）、ADP 5435（≈470，**11.6 倍**）。
     原因就是空穴声子限迁移率里光学声子（光学形变势 d₀）贡献大。
     **价带形变势本身合理**（VBM D(Γ)=1.74 eV ≈ 文献 a_v ≈ 2 eV），不要往它上面找原因；
     「价带形变势绝对归一化」已按 2026-09-24 裁定**不立项**，改用 AMSET 原文对照。
     **报数时**：用 AMSET 报非极性材料的**空穴**迁移率必须标注此缺失；极性材料不受影响。
3. **偶极近似、无四极矩**：Poncé 等（PRL 2023）指出准确预测二维迁移率需规范协变四极矩，
   **MoS2 具体数字**：Poncé PRL 2023 —— 忽略四极矩时室温 Hall 迁移率误差为
   电子 23%、空穴 76%。
4. 二维 Fröhlich **不弱**：Sohier 2016 用截断库仑的 DFPT 表明它比早期第一性原理研究的
   假设强得多；Kaasbjerg 2012 把 MoS2 的室温限制归因于光学形变势 + Fröhlich。

## 5.2 每条命题的来源状态（按文献核对，2026-09-16）

| 命题 | 来源状态 |
|---|---|
| 二维弹性常数 = 三维值 × 晶胞长度；AMSET 超胞里应直接用未换算的 slab 值 | **有文献支持**（二维标准换算）+ 本流程数值佐证（3.08 vs 3.06） |
| 二维极化用 r = c(ε_slab − 1)/2；Keldysh 屏蔽长度 r0 = 2πα_2D | **有文献支持**（Royo–Stengel 强调二维响应必须与超胞长度无关） |
| 各向异性 Keldysh：1/(k + A kx² + B ky²)，与 q·ε(q) = ε_env q + q·r·q 同构 | **有文献支持** |
| POP 必须用二维核（周期镜像的虚假相互作用） | **有文献支持**（平面波周期边界下不能直接套三维核） |
| ω_LO(q)² = ω_TO² · eps0(q)/eps∞(q) 的**完整表达式** | 推导；两个极限与文献一致，但表达式本身非文献原文 |
| ADP 用方向相关的弹性刚度（面内 Christoffel 方程） | **有文献支持**（Lang 等的各向异性处理是其特例） |
| **G_s = c·G_2D**（含 IMP/PIE 核里的 c² 因子） | **纯推导，无文献**；数值上已有：ADP 比值测试（c/t）与绝对归一化（-2%）、POP 绝对归一化（散射率比值 1.01±0.03）、IMP 绝对归一化（-16%，差异归因于屏蔽口径）。**PIE 仍未验证** |
| "介电差形式"与"玻恩电荷形式"等价（取 ω_LO(q) 时） | 推导（与 Sohier 的定性结论吻合，无文献直接给出） |
| 二维压电核的具体表达式 | 部分文献（Kaasbjerg 2013 的物理结论；公式原文未核到）；D'Souza 未含长程压电（只说低温重要） |
| ZA 支不耦合 | **有条件**：仅对存在 σh 的材料成立（见 5.1） |

## 6. skill_rev 防呆

step8.x 的 `_SKILL_REV` 写进 comparison_summary.txt。改脚本即提交 GitHub 并 bump，
防陈旧副本（step12 / step9b / step13 三次被同一形状的 bug 咬过）。

---

> 结果文档（未发表，不进 git）：`/mnt/d/tf_data/jzz/jap/comparison_vs_literature.md`
