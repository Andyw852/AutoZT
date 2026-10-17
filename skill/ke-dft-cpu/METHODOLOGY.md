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
- `ke-dft-cpu` 的 `LAYER_THICKNESS`（决定 σ、κ_e，以及 8.1 报出的 κ_p）。2026-10（补丁 0016）起写在
  **本材料共用 step.conf**（`project_setting/templates/step.conf` 的 `[params]`），S8 / S8.1 / S8.4
  同读这一个；旧做法写在 `templates/step8_amset/step.conf` 只对 S8 生效，2D 默认不跑 S8 之后要搬过来。
  S8.1 的取值顺序：`THICKNESS_A` → 共用 `LAYER_THICKNESS` → S8 / S8.4 的 `2d_correction.json` →
  按结构现算 vdW；实际来源写在 `paper_mu_scan_*.csv` 首行的 `t_source`。
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
`gen_step10_amset.py` / `gen_step14_amset2d.py` 里 `_dielectric_2d_inplane`
按方向扣真空：**面内**并联 `ε_layer − 1 = (c/t)(ε_slab − 1)`，**面外 zz 串联**
`1/ε_slab = (t/c)/ε_layer + (1 − t/c)`（层与真空去极化；c/t 复用弹性同款）。
★ 口径（[0031]）：S8.4 的二维散射核真正"消除真空"靠的是插件里的
`r = c(ε_slab,∥ − 1)/2`（取 `2d_correction.json` 的原始 slab 张量）；settings.yaml 里这个
layer ε 只在旧的三维核（S8）里起作用，不进 S8.4 的数值。

**但散射核本身还是三维的** —— 这一条不再靠「只作比值」糊过去：2D 项目改走
`step8.4_amset2d`（可选步骤组 `amset2d`；`auto_dim: 2d`，2D 材料 init/register 时自动写
`amset2d: true`），运行期插件把四种散射核换成二维形式。原版 `step8_amset` 保留给三维项目：
它在可选步骤组 `amset3d` 里（`auto_dim_off: 2d`，补丁 0016），2D 材料新建时自动写 `amset3d: false`、
不跑 S8（三维核的 POP/IMP 对 2D 不适用，只能作对照）；要三维核对照就把它改回 true，或
`-j S8_kappa start` 按需启用。已有材料的配置里没有这个键，行为不变。`zt-dft-cpu` 的 S20 对 2D
只认 S8.4 的 transport.json。

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
| PIE | `k_B T*(2*pi)^2*(q q:e_2D.v)^2*q^2/[C_2D*(q.eps_0(q)+q_TF)^2]` | `c^2*(2*pi)^2*sum_b (q q:e.v_b)^2*q^2/(w_b*(q.eps_0(q)+q_TF)^2)`（e 为离子弛豫总张量，故屏蔽用静态 ε0） |

其中张量 Keldysh 介电函数 `q.eps(q) = eps_env*|q| + q.r.q`，`r = c*(eps_slab,|| - 1)/2`
（2x2 张量，bohr）；`q_TF = 2*pi*c*dn_3D/dmu`（二维 Thomas-Fermi）。

**PIE 的推导核对（[0037]）**：面内压电极化 `P = e_2D:S`，面电荷 `ρ = -iq.P = q^2 (q q:e_2D.v) A`，
二维泊松 `φ = 2*pi*ρ/(q.eps(q) + q_TF)`，均分 `|A|^2 = k_B T/(Area*C_2D*q^2)`，得表中 `G_2D`
（与 Kaasbjerg 2013 的 `e11^2 [cos^2(3θ)/C_LA + sin^2(3θ)/C_TA]` 角分布一致）；代入 `e_2D = c*e`、`C_2D = c*C` 即插件 factor。
不屏蔽、`r = 0` 时 PIE 核与 ADP 核同形（常数 |M|²，等效 `D = 2*pi*c*(q q:e.v)/eps_env`），绝对归一化继承 ADP 已验证的 `G_s = c*G_2D`；
有 Keldysh 屏蔽时仍未端到端验证。

**3D PIE（[0037]，S8，插件 `step8_amset/amset_pie3d.py`）**：原版 AMSET 的 PIE 有三处口径差——
屏蔽写成 `1/(q^2+β^2)`（应为 `q^2/(q^2+β^2)^2`：压电势是库仑势，屏蔽后 `φ ∝ q^2/(q^2+β^2)`，再乘 `|A|^2 ∝ 1/q^2`）、
用 ε∞（e 是弛豫离子总张量，应配静态 ε0）、纵场用 `ε^-1.q`（应为 `q.e:(q v)/(q.ε0.q)`，各向同性时相同）。插件核：
`factor = (4*pi)^2 Σ_b [q.e:(q v_b)]^2/(c_b (q.ε0.q)^2) * q^2/(q^2 + κ/(q.ε0.q))^2`，`κ = 4*pi*dn/dμ`，prefactor 同原版 `k_B T`；
屏蔽由 `PIE_SCREENING`（出厂开）单独控制，不碰全局 `free_carrier_screening`（POP 不变）。
ε 各向同性、不屏蔽时与原版逐点相同（`test_pie3d.py`，锁住绝对归一化）。S8 从 `step6_elastic/OUTCAR` 读
`PIEZOELECTRIC TENSOR (including ionic contribution)` 写 settings；**S6 出厂没开 `LCALCEPS`/`LEPSILON`，OUTCAR 里没有压电张量**，
此时 PIE 自动去掉并告警（2D 的 S8.4 同样如此）。中心对称（`|e|max < 1e-3 C/m^2`）也去掉。`PIE_KERNEL = amset` 回到原版核作对照。

**POP 的两个输入量 r∞ 与 Δr（2026-09-16 新增自检）**：由**同一次** DFPT 的原始 slab 张量给出
`r∞ = c*(eps_inf,|| - 1)/2`、`r0 = c*(eps_0,|| - 1)/2`、`Δr = r0 - r∞ = c*(eps_0,|| - eps_inf,||)/2`。
gen_step14 现在会算出这三个数、打印并落盘到 `2d_correction.json` 的 `r_inf_delta_r_check`，
插件每次运行也在横幅里复述一遍（`[amset2d] POP: r_inf=... delta_r=...`）。闸门：
`Δr < 0` 报警（离子介电缺失/读错）；`Δr/r∞ > 0.6` 报警（偏大；先核对 ε0/ε∞ 是否同源、
ε_slab 是否真空收敛——早期台子曾因拿**体相**外推移当 slab 值得到 Δr = 43.5 Å，比 Sohier 的 0.53 Å 大 ~80 倍）。
文献锚点：MoS2 r∞ ≈ 46.5 Å / Δr ≈ 0.53 Å（Δr/r∞ ≈ 1.1%）；h-BN r∞ ≈ 7.6 Å / Δr ≈ 2.8 Å（≈37%）。

**ADP 的等价说法（最要紧的一条）**：`C_2D = c*C_slab`，所以 `G_s = k_B T*D^2/C_slab`，
即 settings.yaml 里的 `elastic_constant` 必须写**原始 slab 值**，**不能再乘 c/t**。
旧做法乘了 c/t -> ADP 散射率偏小 `t/c` 倍 -> ADP 迁移率偏大 `c/t` 倍
（CrS2：c/t = 3.06；实测旧流程 mu_ADP/正确值 = 3.1，见 VERIFICATION.md 的 T2）。
ADP 与 q 无关，所以这条对原版 `step8_amset` 同样成立：2026-10（补丁 0014）起 S8 的 2D 路径
也改为原始 slab 值，c/t 只用于介电扣真空与 σ/κ 的厚度归一化（`2d_correction.json` 的
`elastic_constant_used_GPa`）。此前 S8 的 2D μ_ADP 偏大 c/t 倍，需要重算。

### 4.3 pop_frequency 取面内极性模

AMSET 自带的 `amset phonon-frequency` 对**三维球面**方向做 Born 电荷/本征矢加权，
slab 里会把面外极性模也算进去。二维 Fröhlich 只用面内极性模，所以本步改用
`amset2d_pop_freq.py`：同一套公式，只把 q 方向限制在 xy 面内。实测 CrS2：
三维口径 12.45 THz（含 2.4% 的 15.07 THz 面外模），面内口径 12.38 THz（纯 E' 简并模）。
插件再乘 `sqrt(eps_0(q)/eps_inf(q))`（二维 LST）得到 `w_LO(q)`。

### 4.4 形变势的构型口径（ionrelax）与参考能级（2D 默认真空）

**构型口径（2026-09-16 修订；2026-10 按代码核对）**：step7b 生成的**两份** h5（芯势口径的
`deformation.h5` 与真空口径的 `deformation_vac.h5`）都读 `deform-*/ionrelax/`
（离子弛豫构型，有则用）。物理依据：q→0 声学声子 = 均匀应变 + 内坐标弛豫；一致性依据：
ADP 核分母是 Christoffel 刚度（VASP `TOTAL ELASTIC MODULI`，含离子弛豫）。实测 MoS2 单层
xx 0.5%：电子 E1_vac 由 6.58（离子固定）变 8.51 eV（离子弛豫），差 **29%**。

**哪些构型做离子弛豫（[V187] 起按维度）**：step7（`gen_step9_deform.py` 的 `select_ionrelax_dirs`）给形变目录建
`ionrelax/`（固定晶胞 `ISIF=2`、`IBRION=2` 弛豫内坐标，不加选择性动力学，再接 `ISYM=0` 静态段）：

| 体系 | 形变构型 | `deform_geometry` | z 方向 |
|---|---|---|---|
| 2D | 面内分量（xx/yy/xy）ionrelax，面外分量（zz/yz/xz）clamped | `ionrelax(in-plane)+clamped(out-of-plane)` | z 是真空：面外应变拉/剪的是真空，二维核不用，不弛豫；但面内应变下原子的 **z 坐标放开**（层厚随应变变化，是二维离子弛豫响应的一部分） |
| 3D | 全部分量 ionrelax（[V187]；`IONRELAX_3D = false` 回到全 clamped） | `ionrelax` | z 是真实晶向，六个分量一视同仁 |

> [V187] 之前 3D 一律离子固定（`clamped`），与含离子弛豫的弹性常数不同源。受影响的是有内应变的结构：金刚石/闪锌矿的
> 剪切分量（Kleinman 内应变 → TA 支 ADP）、R-3m（Bi₂Te₃ 类，有自由 z 参数）连正应变也受影响；原子都在反演中心的结构
> （岩盐等）受力为零，弛豫段 1 步即停，结果不变。**已算完的 3D 项目不会自动重算**：retry S7_deform 才会给全部构型补
> `ionrelax/`（ck_deform 判未完成 → 重交），再 S7.1 → S8 / S8.2。
> ⚠ **回退不静默**：没有 `ionrelax/` 的构型按离子固定口径处理时，gen 会打 WARN 并把
> `deform_geometry`（ionrelax / clamped / mixed）写进 `band_edges.json`，进而出现在
> settings 注释、`2d_correction.json` 与 8.3 汇总表里，跨项目比较前先核对这一项。
> 2D 的 `ionrelax(in-plane)+clamped(out-of-plane)` 是正常口径：S8.4（二维核）与 8.3 不再告警；只有 S8（三维核，2D 只作对照）
> 会用到面外分量，那里单独注明。

**弛豫力判据（[V187]）**：`IONRELAX_EDIFFG` 由 -0.02 收紧到 **-0.005 eV/Å**（与 S1 为弹性收紧的同值）。带边本征值对内坐标是
**一阶**的（斜率 = 光学形变势 D_opt，1–10 eV/Å），残余位移 δu ≈ F_tol/k 落到 E1 上的误差上限 ≈ D_opt·δu/γ（γ = 0.005）：
-0.02 → ~1 eV，-0.005 → ~0.3 eV。旧注释"残余力对带边影响 ~1E-4 eV"是把总能（二阶）的量级套到了本征值上。MoS₂ 粗估：
±0.5% 应变诱导的内坐标力只有 ~0.04 eV/Å，与 -0.02 同量级。S7.1 逐构型写 `band_edges.json` 的 `ionrelax_check`
（首/末步最大力、离子步数、是否跑满 NSW、判据是否偏松；末/首 > 0.25 记 `weak`）。判据偏松的旧结果默认复用、只告警；
S7 的 step.conf 写 `IONRELAX_REBUILD_LOOSE = true` 再 retry 只重建这些构型（`IONRELAX_EDIFFG` 也可在那里改）。

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

**LOCPOT 与本征值同构型（[patch_ionrelax_locpot]，2026-10-15）**：某构型的本征值取自 `ionrelax/vasprun.xml`，
它的 LOCPOT（以及 POSCAR）就只从 `ionrelax/` 读；`ionrelax/LOCPOT` 缺了 S7.1 报 `RuntimeError`
（"ionrelax 缺 LOCPOT，retry S7_deform 补静态段"），不再退回刚性构型的 LOCPOT。真空口径 h5 生不出来时 S7.1 直接
FAIL，不再静默退回 core 口径；确实要 core 就在**材料共用** step.conf 显式写 `DEFORM_REF = core`（S7 与 S7.1 同读）。
retry S7_deform 时，2D + 真空口径 + ionrelax 已完成 + 缺 `ionrelax/LOCPOT` + gen 前还没有
`step7b_deform_read/deformation_vac.h5` 的构型**只补静态段**：静态段的 OUTCAR/vasprun.xml/OSZICAR 改名
`*.nolocpot-<时间>`，submit.sh 重写成只跑静态段（`CONTCAR.relax` → POSCAR、`INCAR.static` 改 `ISTART=0`/`ICHARG=2`），
刚性单点与弛豫段不重算。已有真空口径 h5 时不触发 —— S7.1 跑完后清 LOCPOT 是允许的，此后 retry S7 不会整批重算静态段
（但那次 retry 会让 S7.1 重排、重读时报缺 LOCPOT，届时再 retry S7 一次即只补静态段）。

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
| PIE 绝对归一化 | 2D：不屏蔽、r=0 极限下核与 ADP 同形，归一化继承 ADP 的 -2%；带 Keldysh 屏蔽的端到端**未验证**。3D（[0037]）：各向同性不屏蔽极限与原版 AMSET 逐点相同；端到端迁移率未与文献对照。默认关闭，论文中如实写明 |
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
   本身也处理不了）。补丁 0016 起自动判断（`ke_common.za_coupling_check`：spglib 对称操作里
   有没有笛卡尔部分为 I − 2nnᵀ 的操作，n 为层法向，含滑移镜面），结论写进 S8.2 的
   `dpt_result.json` 与 S8.4 的 `2d_correction.json`（键 `za_coupling`），没有 σh 时 gen 打 WARN。
   只能标注，补不上：DPT 与 AMSET 都不含 ZA 一阶耦合，文献 DPT 同样不含。
2. **缺非极性光学声子（光学形变势 ODP）**：AMSET 实现的机制只有 **ADP / POP / IMP / PIE**
   （`amset/scattering/elastic.py` = ADP+IMP+PIE，`inelastic.py` 只有 `PolarOpticalScattering`；
   官方文档 Summary of scattering rates 同）。**非极性晶体 POP 前因子为零**
   （`scattering/calculate.py:842` 注释明说），于是**整条光学声子通道缺失**：
   - **2D**：MoS2 的同极 A1 模与谷间光学模耦合很重要 -> **总**迁移率系统偏高。
     **[0033] 起 S8.4 可选补上**：常数形变势模型（Jin 2014 / Kaasbjerg 2012）的谷内 ODP 与谷间 IVS，
     step.conf 写 `DP_OPTICAL_E/H` 通道表（参数抄文献表），谷表取自 S8.2 多谷；开了 IVS 的载流子 ADP 只留谷内，
     避免重复计数。谷内 ODP 自 [0041] 起出厂由冻声子参数打开（`odp_frozen` 组默认开）；谷间 IVS 仍要手填，出厂不算。
     细节见 `step8.4_amset2d/README.md`「ODP 与谷间散射」。
     **[0035] 谷内 ODP 的参数可以自动算**：开 `odp_frozen` 组后，S7.2 对 S5 的 Γ 点非极性光学模做 ±u 冻声子单点，
     S8.21 在 S8.2 各谷的 k 点取真空对齐的带边能量算 D = dE/du（位移 Δ_κ = u·√(M/m_κ)·e_κ，与常数形变势模型的
     |g|² = ħD²/(2NMω) 同一归一化），S8.4 自动用上。只有零阶（q 无关）谷内耦合；谷间 IVS 仍要文献参数。
     细节见 `step8.4_amset2d/README.md`「冻声子 ODP 参数」。
     **[0036] 3D（S8）也可选补上**：同一常数形变势模型（公式本来与维数无关：|M|² 的 1/V 与末态求和的 V 抵消，
     G 取每原胞态密度；三维即 Conwell 的非极性光学声子率），插件 `step8_amset/amset_dpo3d.py` 把 ODP/IVS 注册进
     原版三维 AMSET，谷归属用完整三维距离（27 个周期像）。材料 step.conf 写 `DP_OPTICAL_E/H` +
     **必须手填** `DP_OPTICAL_VALLEYS_E/H`（三维没有 S8.2 多谷谷表，k 点写 3 个分量），如 Si 空穴
     `DP_OPTICAL_H = G>G:6.0@63`、`DP_OPTICAL_VALLEYS_H = G:0,0,0`；`ADP_INTRAVALLEY` 规则同 2D。
     transport.json 多出 ODP/IVS 单机制迁移率，S8.3 自动带上。
     **谷内 ODP 出厂用冻声子参数**（与 2D 同一规则）：开 `odp_frozen` 组后 S7.2 在 3D 上照样跑（照抄 S7 undeformed），
     S8.21 不对齐（没有真空）、谷表在 S7 网格（ISYM=0）上自动找带边极值（或用手填 `DP_OPTICAL_VALLEYS_E/H`），
     S8 的 `ODP_FROZEN=auto` 读 `odp_params.json`：算出来的优先，手填谷内值只作对照（`ODP_FROZEN = off` 才用手填）；
     **谷间 IVS 不自动算**，只有手填 `X>X'`/`X>Y` 才开（谷标签用 S8.21 的 `odp_summary.txt`）。
     自检 `step8_amset/test_dpo3d.py`（纯 numpy：2D 口径逐字节不变、3D 抛物带对 Conwell 解析式、3D 自动谷表、冻声子接线、S8 重排）。
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

## 5.3 波函数去对称化的相位（V115，2026-09-28）

AMSET 用 IBZ 波函数时要把系数按对称操作展开到全网格。其相位因子只在满足
非 TR：(R+I)τ ≡ 0、TR：(R−I)τ ≡ 0（mod 1）的操作上精确（逐操作判据，`ke_common.op_phase_exact`）；
不满足时真实重叠被算坏（与有无反演中心无关）。三种合法做法：① 全网格波函数（ISYM=-1，AMSET 走
from_data）；② S3 的 `align_origin` 把原点移到让坏操作归零的位置（简单空间群可行）；③ 相位补丁
`amset_desym_fix`（`DESYM_FIX=on`，任何结构都精确；实数据验证通过前默认关）。判据与 AMSET 用同一套
操作（symprec 0.01 Å）。推导、模型检验与接线见 `step8.4_amset2d/VERIFICATION.md` V115。

## 6. skill_rev 防呆

step8.x 的 `_SKILL_REV` 写进 comparison_summary.txt。改脚本即提交 GitHub 并 bump，
防陈旧副本（step12 / step9b / step13 三次被同一形状的 bug 咬过）。

---

> 结果文档（未发表，不进 git）：`/mnt/d/tf_data/jzz/jap/comparison_vs_literature.md`

## [0056] phonon-limited 掺杂扫描按浓度拆两段（2026-10-17 用户定）

- **为什么**：[0055] 加了低浓度段以后，`<M>_phonon_limited.png` 的 ± 浓度轴原点下探到 1e9，载流子调控区（1e17–1e21）
  只占每格两端一小截，看不清。
- **现在**：不再出整轴的 `<M>_phonon_limited.png` / `_phonon_limited_data_<T>K.csv`，改成两个文件夹，每个里面一张 9×5 图
  （温度 × μ/S/σ/PF/κe，基础 4 条曲线）+ 每温度一张表：
  - `transport_plots/phonon_limited_carrier_tuning_1e17-1e21/`：`<M>_phonon_limited_carrier_tuning.png` + `_<T>K.csv`，
    |n| ∈ [1e17, ∞)，轴与 [0055] 之前逐字相同（原点 1e17，±4 十年），S 纵轴 ±1000 μV/K；
  - `transport_plots/phonon_limited_low_doping_<NLO>-1e16/`：`<M>_phonon_limited_low_doping.png` + `_<T>K.csv`，
    |n| ∈ [NLO, 1e16]，原点 NLO（最小点所在十年），S 纵轴 ±2500 μV/K；没有 DOPING_LOW 时不出这个文件夹。
  - 两段上下限都含端点；(1e16, 1e17) 出厂网格上没有点（3D 低段止于 1e16、主网格起于 1e17；2D 换算后 ~2.9e15 / ~2.9e17）。
    本征点 n = 0 不进任何一段（在 `intrinsic_carrier/`）。
- **r2（同日用户要求）**：`<M>_transport.png`（基础 4 条 + 本轮机制对照）同样拆成
  `transport_carrier_tuning_1e17-1e21/<M>_transport_carrier_tuning.png + _<T>K.csv` 与
  `transport_low_doping_<NLO>-1e16/<M>_transport_low_doping.png + _<T>K.csv`，表含全数据集（只截浓度段）；
  顶层 `<M>_transport.png` / `_transport_data_<T>K.csv` 不再出。
- `datasets/`、`amset_scan/` 仍用整轴，不变。`tools/make_transport_plots.py` 同步（它的 transport 图是 DPT 路线 +
  AMSET phonon + impurity + ADP only；低浓度段 σ/PF/κe 纵轴改自动范围）。

## [0055] 低浓度段 DOPING_LOW（2026-10-16 用户要求）

- **为什么**：出厂掺杂网格 3D 1e17–1e21 cm⁻³、2D 1e11–1e14 cm⁻²，和本征点（n = p）之间空着好几个十年；
  电声文献算 phonon-limited 迁移率用的正是低浓度、非简并、单载流子口径（固定低浓度，μ 只积一种能带——AMSET 的
  `_calculate_mobility` 按掺杂符号只积导带/价带，本来就是分开的；σ、S、κe 对全部能带积分，含双极）。
- **出厂值**：3D `DOPING_LOW = "-1e16:-1e10:4, 1e10:1e16:4"`（cm⁻³，n/p 各 1e10、1e12、1e14、1e16）；
  2D `"-1e9:-1e3:4, 1e3:1e9:4"`（cm⁻²，n_2D ≈ n_3D × t、t ≈ 1 nm 取整；文献口径 n_lit = n_2D/t，MoS₂ t≈0.67 nm
  即 1.5e10–1.5e16 cm⁻³）。掺杂点 11 → 19，AMSET 计算量约 ×1.7。step.conf `DOPING_LOW = off` 关，或写自己的范围串。
  3D 的 S8 以前 step.conf 不能改 DOPING，DOPING_LOW 两边都能改。
- **E_F 必须换解法**：AMSET `FermiDos.get_doping` = (nelect − Σ f·DOS·dE)/V，`hole_conc` = (nelect − Σ_VB)/V，
  都是两个 ≈ nelect 的数相减。3D 1e10 cm⁻³ 每胞只有 ~4e-13 个载流子，理想合成 DOS 上相减式误差已 ~0.5%（逼近 1% 护栏），
  真实 DOS 的 nelect 归一误差只会更大。`amset_fermi_fix.small_doping_fermi`：每胞 |d|·V < `SMALL_DOPING_PER_CELL`
  （1e-7）时改解 ln n − ln(p + d) = 0（n 型）/ ln(n + d) − ln p = 0（p 型），n、p 分开 log-sum-exp，并把 n、p 直接填进
  electron_conc / hole_conc（迁移率分母）。合成抛物带验证：每胞 4e-13、窄带隙 900 K（n_i/|d| ~ 1e8）下 (p−n)/d 误差 ≤ 4e-8，
  E_F 与 Boltzmann 解析式差 ~1e-6 eV。出厂网格最小点每胞 ≥ ~1e-6，不进这条分支，结果逐位不变。
- **出图**：S8.3 浓度轴（`_doping_axis`）跟着 AMSET 实际掺杂网格：最小非零 |n| ≥ 1e17（元胞口径）时与旧版逐字相同
  （原点 1e17、下限 5e17、±4 十年）；有低浓度段时原点下探到最小点所在十年，全部 AMSET 点都画，每 2 个十年一个刻度，
  S 纵轴 ±1000 → ±2500 μV/K。`tools/make_transport_plots.py` 同算法。mechanism-resolved 图的电子主值、μ(T) 表、τ 对比
  都取 |掺杂| 最低的一档，有低浓度段后自动用到 1e10 cm⁻³ / 1e3 cm⁻² 这一档（更接近文献 f ≈ 0）。
- **低浓度下的 IMP**：N_ii = |n − p| 很小，IMP 率可能低于 AMSET 插值的"非零"阈值 1e6 s⁻¹，按 [0043] 的死片规则处理（≈ 无 IMP），
  物理上也可忽略。

## [0054] 命名与电声文献统一（2026-10-16 用户定）

规则：**"intrinsic" 只指本征载流子**（未掺杂，n = p = n_i，教科书 intrinsic carrier concentration）；散射一律按
物理含义 + 实际机制组合命名（`_phys_label`）：

| 含义 | 图例 / 对外名 | 文献对应 | 内部键 / 文件（不改，兼容旧结果与 zt-dft-cpu） |
|---|---|---|---|
| 全部声子机制、不含电离杂质 | `phonon-limited (ADP+POP+…)` | phonon-limited mobility（Cheng & Liu 2018 也叫 intrinsic mobility） | 数据集键 `intrinsic`、`intrinsic_transport.json`、`subset_<去 IMP>.json` |
| 含电离杂质 | `phonon + impurity (ADP+IMP+POP+…)` | — | 键 `overall`、`transport.json` |
| 只 ADP | `ADP only (acoustic limit)` | DPT / Bardeen–Shockley 口径 | 键 `ADP_only`、`subset_ADP.json` |
| 其余声子子集 | `partial phonon (ADP+POP)` | — | `subset_*.json` |
| 本征载流子 | `intrinsic carriers (undoped, n = p = n_i)` | intrinsic carrier concentration | `intrinsic_carrier/`、`undoped_ni.json` |

对外文件改名：`transport_plots/<M>_intrinsic.png` / `_intrinsic_data_<T>K.csv` → `<M>_phonon_limited.png` /
`_phonon_limited_data_<T>K.csv`；`intrinsic_carrier/<M>_dilute_mobility.*` → `<M>_mechanism_resolved_mobility.*`。
S8.3 每次重建 `transport_plots/`，旧名文件随之清掉。`step8.3_output/intrinsic/` 下的内部表、`comparison_vs_T.csv` 的
列名 `mu_AMSET_intrinsic` / `mu_AMSET_all` 不改（图例已改）。datasets/ 的文件名仍是纯机制组合（如 `<M>_ADP+POP.png`）。

- **本征点去重**：N_ii = |n − p| = 0，IMP 率为 0，含 IMP 的数据集与去掉 IMP 的同机制集逐点重合。
  `<M>_intrinsic_carrier` 按"去掉 IMP 后的机制组合"去重，优先留本来不含 IMP 的那套；只剩含 IMP 的（旧回退文件）
  图例标 `[IMP = 0 at n = p]`。出厂（IVS 默认关、PIE 需 step.conf 显式开）本征点通常剩：phonon-limited
  （ADP+POP+ODP，有 S8.21 冻声子 ODP 时）、ADP+POP（来自 overall_noODP）、ADP only。
- **f ≈ 0 判据 η**：`amset_fermi_fix` 在 `undoped_ni.json` 多记带边 `E_v_eV` / `E_c_eV`（与 E_F 同参考、含 scissor）。
  S8.3 算 η_c = (E_c − E_F)/k_BT、η_v = (E_F − E_v)/k_BT 与 f(CBM)，写进 CSV；min(η) < 5（f > ~0.7%）的温度打 WARN、
  图上灰底——窄带隙/高温时本征点不再是文献的稀薄极限。[0054] 之前的结果没有带边，重跑 S8/S8.4 后才有。

## [0053] 稀薄极限（f ≈ 0）迁移率的机制分解（2026-10-16）

- **本征点 n = p 等价于稀薄极限 f ≈ 0**：E_F = E_mid + (k_BT/2)·ln(N_v/N_c) 在带隙中部，带边占据
  f ≈ exp(−E_g/2k_BT)（MoS₂ 300 K 约 1e-15）——无泡利阻塞；N_ii = |n − p| = 0，无 IMP；n_i 极小，自由载流子屏蔽可忽略。
  这与 Cheng & Liu, JACS 2018, 140, 17895 的"费米能级放带隙中间、只看带边态"是同一个占据口径。
  前提是 E_g/2 ≫ k_BT；窄带隙材料高温下 n_i 变大，这条不再成立（核对点会报 WARN）。
- **AMSET 在 0 点给的是空穴迁移率**：`amset/core/transport.py` 的 `_get_band_idx` 对 doping ≥ 0 只取价带，
  `_calculate_mobility` 除以 `hole_conc`（`amset_fermi_fix` 已把它填成 n_i），所以 0 点的 `mobility[...]` = 空穴、只含价带。
  电子在 0 点拿不到，取 |掺杂| 最低的 n 型点（出厂 1e11 cm⁻²，非简并），总值用去 IMP 的子集。
- **出图** `transport_plots/intrinsic_carrier/<M>_mechanism_resolved_mobility.png/.csv`（[0054] 前叫 `dilute_mobility`）：只取这一个浓度，μ(T) 分电子/空穴两栏，
  画各单机制 μ_i（`transport.json` 的 `mobility[机制]`，即只开该机制）、去 IMP 总值、Matthiessen 倒数和 (Σ1/μ_i)⁻¹，
  CSV 带 `share_1/mu`（限制机制占比，对应文献"由 LO / LA 限制"的说法）。μ > 1e10 的机制（前因子 0，本征点 IMP、非极性 POP）不列。
  核对点：空穴另取最低 p 型、电子取次低 n 型（都去 IMP），与主值差 > 5% 打 WARN（未到非简并平台或屏蔽已起作用）。
- **与文献的口径差**：AMSET 是 RTA 对 −∂f/∂E 的热平均，单机制 μ 的倒数和只给上界（总值 ≤ (Σ1/μ_i)⁻¹，总值更高会打 WARN）；
  文献取带边单态 μ = eτ/m*，Matthiessen 严格成立。两者差在热平均（高于 LO 发射阈值的态散射更强），不在占据数。
  AMSET 的 ADP ≈ LA+TA、POP ≈ LO，其余光学支与谷间散射要开 ODP/IVS 才有。

## [0052] 本征浓度 n_i 标注、`<M>_intrinsic.png`（[0054] 起改名 `<M>_phonon_limited.png`）的含义、先删旧结果（2026-10-15）

- **`<M>_intrinsic.png` 横轴是掺杂浓度**：这里的 intrinsic 指**本征散射**（AMSET 去掉 IMP 的机制组合），
  与 `<M>_transport.png` 一样是掺杂扫描；只画基础 4 条（全机制 / 去 IMP / ADP only / BT2×DPT-τ）。
  图标题与横轴名已写明这一点。本征**载流子浓度**（未掺杂、n = p = n_i，随温度变）的图在 `intrinsic_carrier/`。
- **n_i 的来源**：AMSET 的 transport.json 不带载流子浓度。`amset_fermi_fix` 解本征点时把
  (T, E_F, n_i) 写进运行目录的 `undoped_ni.json`（只在有 settings.yaml 的 AMSET 运行目录写）；
  旧结果退回解析 `queue.out` 里的 `[fermi_fix] 本征点 <T> K：E_F = … eV，n = p = … cm^-3` 行。
  `intrinsic_carrier/<M>_intrinsic_carrier.png` 第一栏画 n_i(T) 并逐点标数值，CSV 每行带 `n_i_cm-3`、`E_F_eV`
  （2D 另有 `n_i_2D_cm-2 = n_i × c`、`n_i_lit_cm-3 = n_i × c/t`）。两处都没有时图标题注明"重跑 S8 后可得"。
- **先删旧结果**：`transport_plots/` 是 S8.3 独占的输出目录，每次 S8.3 先清空再重建；本地回拉
  （fetch_extra / fetch_all 整目录）按镜像处理，远端已不存在的旧文件本地也删除（清单 `.tf_fetched_manifest.json`）。

## [0040] 本征载流子点（UNDOPED_POINT）与"intrinsic"命名

- **命名**：S8.3 的图例与 `datasets/` 文件名一律用**实际散射机制组合** `amset_<机制+机制>`（如 `amset_ADP+POP+ODP`、
  `datasets/<M>_ADP+POP+ODP.png`），不再用 "intrinsic / no IMP / overall" 这类代称；"intrinsic" 只指**本征载流子**
  （未掺杂、E_F 由 n = p 决定），单独出 `transport_plots/intrinsic_carrier/<M>_intrinsic_carrier.png/.csv`
  （σ、S、PF、κe 随温度；2D 的 σ、κe、PF 按层厚 t 归一）。S8.3 内部表 `step8.3_output/intrinsic/` 与 `intrinsic_*.json`
  文件名不变（兼容）。
- **掺杂网格**：S8 / S8.4 出厂在 DOPING 后加一个 0（`UNDOPED_POINT = true`）。1e11–1e14 cm⁻²（2D）/ 1e17–1e21 cm⁻³（3D）
  是人为扫描范围，不是本征浓度；本征浓度 n_i = √(N_c N_v)·exp(−E_g/2k_BT) 因材料、温度而异（MoS₂ 300 K 约 1e-2 cm⁻²）。
- **AMSET doping = 0 的坑**：AmsetData 对 0 掺杂走 `get_fermi_from_num_electrons`（按总电子数相对误差贪心搜索）。
  有带隙时本征载流子每胞远少于容差，整个带隙"误差都为 0"，E_F 落在带隙内任意格点（合成抛物带 1.8 eV、300 K 实测偏 0.12 eV，
  σ、S 差约 2 个数量级）；electron_conc / hole_conc 也不填（迁移率除以 0）。`amset_fermi_fix.py` 在 DOS 有带隙时改解
  ln n = ln p（log-sum-exp，不相消不下溢），与解析式 E_F = E_mid + (k_BT/2)·ln(N_v/N_c) 一致到 1e-6 Ha，并填好 n、p；
  无带隙（金属）照旧走原搜索。本征点 |n − p| ≈ 0，IMP 不起作用（全机制 ≈ 去掉 IMP 的那套）。
  `<M>_intrinsic_carrier` 不画 μ（σ、S 是双极量）；本征点的单载流子迁移率见 [0053] 的 `<M>_mechanism_resolved_mobility`。

## [0041] 缺陷平衡载流子浓度（DEFECT_DOPING）与 odp_frozen 默认开

- **两种"未掺杂"浓度**：[0040] 的本征点是无缺陷理想晶体（n = p，完整 DFT 态密度）；defect-dft-cpu 给的是有本征缺陷时的
  平衡浓度（p0 − n0 + Σq[D^q] = 0，单抛物带有效态密度，含化学势顶点与可选生长温度冻结）——后者才对应实际未掺杂样品。
- **接线**：S8 / S8.4 的 gen 读 `../defect-dft-cpu/step4_analysis/formation_energy_results.json` 各顶点的 n0、p0，
  把净浓度 p0 − n0 加进 AMSET 掺杂（3D cm⁻³；2D 换面密度 = 净浓度 × 缺陷超胞 V/A），记录写 `defect_doping.json`；
  S8.3 出 `transport_plots/defect_equilibrium/<M>_defect_equilibrium.png/.csv`（σ、S、PF、κe 随温度，每顶点一色，
  实线全部机制、虚线去掉 IMP）。带电态是垂直能量（NSW=0）时 auto 不用。
- **局限**：缺陷平衡只在 defect 技能的温度下求；AMSET 各温度用同一净浓度（等价于缺陷冻结）。IMP 的杂质浓度
  **N_ii = compensation_factor × |n−p| / defect_charge**，其中 compensation_factor / defect_charge 是 AMSET 的设置项，
  默认 **2 / 1**（`amset/defaults.yaml`）；[0053] 起这两个量**显式写进** S8/S8.4 的 `settings.yaml`（可在 step.conf 覆盖），
  也就是实际按 **2×|n−p|** 计 —— 但它仍只是"净载流子浓度的倍数"，真实的补偿缺陷（施主受主成对）没算进去
  （实际电离杂质更多 -> 这里的 IMP 偏低）。两边带隙要同源（defect 的 E_gap 与 ke 的 HSE scissor）。
- **odp_frozen 默认开**（用户 2026-10-13 决定）：每材料多 1 + 2×非极性光学模数个 VASP 单点；已算完的材料补跑 S7.2/S8.21，
  参数出来后 S8/S8.4 重排一次。材料 tf_*.yaml 写 `odp_frozen: false` 关。

## [0042] ODP 审查后的修正与局限

五路独立审查（速率公式、冻声子 D、AMSET 接线、参数接线、数值复现）确认：速率公式与 Kaasbjerg Eq.31 / Conwell 一致到 1e-6，
位移归一化 Δ = u·√(M/m)·e（M = 原胞质量）与 Kaasbjerg Eq.3/A1 同一约定（D = dE/du 即 Γ 点耦合），本征矢取质量加权正交归一的
那一组，简并带迹恒等式与态密度配套，K/K′ 平均 D² 不会相消，AMSET 里的率与解析值逐点相同，2D 没有误乘 c。修正：
- **泡利阻塞**：ODP/IVS 率改为逐掺杂的线性化出射率 (N + f(x+ħω))·G(x+ħω) + (N + 1 − f(x−ħω))·G(x−ħω)，与 AMSET 的 POP
  同一口径（以前非简并极限、对所有掺杂广播同一份；简并区 PF 偏 20–30%、1e14 cm⁻² 总迁移率 300 K 偏 ~28%）。
- **2D 极性判据**：面内极性看面内模有效电荷；面外极性只对有 σh 的材料剔除（Janus 的同极类模以前会被误删）。
未改、需在报数时注明的局限：常数 D 取 Γ 点裸值（速率上限，MoS₂ 约高 1.37 倍）；IVS 只认手填；极性模的短程零阶耦合缺失；
G 的 6 meV 展宽抹圆阈值（只看 ODP 的 ⟨τ⟩ 偏低约 10%）；轻重带质量不同的简并带有偏差。

## [V186] S7.1 读形变势：各构型 NBANDS 不一致

- **现象**（LS，2026-10-13）：`Bulk calculation has more bands than deformation. Skipping deformation ...` →
  `ERROR: Strains do not cover full tensor`。undeformed NBANDS = 132，deform-NN/ionrelax = 131。
- **根因**：S7 的 INCAR 不写 NBANDS，VASP 默认值要向上取整到"能带并行组数"（核数 /(KPAR·NCORE)）的倍数，换机器或换核数
  （LS 的 ionrelax 在 hanhai25 跑）就会差几条最顶上的空带。amset 的 `check_calculation` 对"形变构型带多"会裁掉多余的，
  对"bulk 带多"却把整个构型跳过，面内分量一丢，应变张量就不全。
- **修法**：S7.1 不直接调 amset CLI，经 `_NBANDS_FIX_CODE` 包一层：bulk 与已读构型统一裁到最小带数。被裁的带若落进 amset
  的能量窗口（VBM − cutoff, CBM + cutoff），就报错退出（exit 3），不静默裁；各构型带数本来一致时什么都不做。差的只是窗口外的空带，
  对形变势、band_edges.json、deformation.h5 的带轴都没有影响，所以**不用重算 ionrelax**，S7.1 重跑即可。
- **不在 S7 的 INCAR 里补 NBANDS**：INCAR 一变，[V168] 的输入核对会把所有已算完的形变单点判成过期、整批重算。

## [V187] 形变势构型：3D 也做离子弛豫；弛豫力判据收紧；z 方向按维度处理

- **3D 离子弛豫**：S7 对 3D 的全部形变构型建 `ionrelax/`（默认开，S7 的 step.conf `IONRELAX_3D = false` 关）。
  理由与 2D 相同：q→0 声学声子含内坐标弛豫，ADP 分母是含离子弛豫的 `TOTAL ELASTIC MODULI`。细节见 4.4。
- **力判据**：`IONRELAX_EDIFFG` -0.02 → -0.005；S7.1 的 `ionrelax_check` 核对回收是否充分；旧结果按
  `IONRELAX_REBUILD_LOOSE` 由用户决定是否重算（不自动重算）。
- **z 方向按维度（全技能核对后修的几处）**：
  - S7.1 `band_edges.json`：带边 hits 加 `E1_zz_eV`；3D 的 `E1_iso_eV` 改为 (D_xx+D_yy+D_zz)/3（以前照搬 2D 只取面内两项），
    `E1_iso_def` 写明口径；`BAND_EDGES_REV` 随之升级（S8.2 对旧文件告警"重跑 S7.1"）。S8.2 的 3D 主路径本来就直接读 h5 的
    三个对角（[0025]），只改了回退路径的说明。
  - S8.4 / S8.3：2D 的 `ionrelax(in-plane)+clamped(out-of-plane)` 不再误报"非离子弛豫口径"；8.3 的张量约化标签按 step1 的
    DIM 判（2D 默认不跑 S8 时以前会被标成 "[3D] (xx+yy+zz)/3"）。
  - S8（三维核）：读 `deform_geometry`，3D 非离子弛豫时告警；2D 注明面外分量离子固定、S8 只作对照。
- **未改**：2D 的面外应变仍是离子固定单点 —— 只为让 AMSET 读到完整张量，二维核不用；未形变参考构型仍是 S1 CONTCAR 的单点
  （中心差分里抵消；AMSET 逐应变取 |·| 后对 |D| 小于 ~s/γ 的分量会留下一个噪声底，s 为 S1/S7 设置差造成的残余力位移所致的带边漂移）。

## [V188] 形变势对称化改正：AMSET 的 h5 是幅值，不是带符号张量

- **错在哪**：V121 的 `dp_symmetrize` 以为 h5 存的是带符号张量 D，对它做 rᵀ D(rk) r 的群平均。实际上 AMSET（0.4.19 与 0.5.1
  相同）在 `get_strain_deformation_potential` 里 `np.abs(diff / strain_amount)`，**逐分量取绝对值**再对各应变平均 ——
  h5 里是幅值 |D_ij|。把幅值当张量旋转：Laue 群里只要有垂直 x、y、z 的镜面（mmm、4/mmm、m-3m、6/mmm，MoS₂ 也是），
  非对角幅值带着 ± 号被平均成 **0**（剪切形变势 → 横声学支 ADP 被抹掉；PbTe/SnTe/Ge 的 L 谷 Ξu/3 整个没了），
  3/6 次轴还把非对角混进对角。另外 AMSET 的非对角存的是 **2·D_xy**（剪切应变取张量分量 ε_xy = δ，ΔE = 2D_xy δ），
  不是张量分量，旧法的混合变换这一点也没处理。
- **AMSET 自己的 |·| 再平均也有偏差**：±应变给 d± = −D ∓ ½Kδ（二阶曲率）或 −D ± s/δ（未形变参考的残差），
  先取绝对值再平均得 max(|D|, ½|K|δ, |s|/δ) —— 小分量（多是非对角）被垫上一个虚假的底。带符号的中心差分把它精确抵消。
  只有**简并能级被应变劈裂**处（Si 的 Γ 价带顶：排序后 E = E0 ± |a|δ），导数没有定义，AMSET 的 |a| 才是合理的幅值。
- **改法**：
  1. S7.1 在同一次读取里（core：包一层 `amset deform read`；vac：自己的脚本）用与 AMSET 相同的 k 映射与分量归属
     再算一份**带符号中心差分**，连同各带与相邻带的能隙、应变幅度写进 `deformation_signed.h5` /
     `deformation_vac_signed.h5`（数据集名不含 `deformation_potentials`，AMSET 不会误读）；
  2. 对称化：非对角先 ÷2 换成真张量，对带符号场做全 Laue 群平均，再 ×2、取 |·| 写回 h5（AMSET 约定的幅值；
     [V189] 之后紧接着再把非对角换成张量分量，见下一节）；
  3. 只在"近简并（能隙 < 2δ·max|D| + 1 meV）**且** ±应变符号相反"的 (带, k) 及其全部对称像上退回 AMSET 的幅值，
     并只用"带符号置换"类操作（对幅值精确：|r D rᵀ| = |r| |D| |r|ᵀ）；Kramers 简并（应变不劈裂）照用带符号；
  4. 没有带符号文件（旧产物）时只做第 3 步的置换子群对称化（精确，但六方的 xx/yy 拉不齐）。
- **标记**：h5 的 `dp_symmetrized` = 2；= 1 是 V121 旧法，S8/S8.4 会告警、`tools/regate_projects.py` 记 STALE-DP，
  `dp_symmetrize` 不会再动它（无法复原）—— **重新 gen S7.1（秒级、不提交作业）再重跑 S8/S8.4/S8.2**。
- **影响面**：带符号可信的点上，h5 的值从 AMSET 的 ⟨|d±|⟩ 变成 |中心差分|：主分量几乎不变，小分量的虚假底被去掉；
  非对角（剪切）形变势恢复（V121 之后的所有结果里它们都是 0）。立方/六方材料的 ADP 迁移率会变，L 谷/偏轴谷材料变化最大。
- 验证见 `step8.4_amset2d/VERIFICATION.md` 的 V188；测试 `test_dp_symmetrize.py`（含真实 AMSET 函数的端到端）、
  `tests/test_v188_dp_symmetrize.py`（纯 numpy）。

## [V189] 形变势 h5 的非对角换成张量分量（AMSET 运行期把剪切耦合算大一倍）

- **问题**：AMSET 原生 h5 的非对角是 2·D_xy（[V188]）。运行期的 ADP（`amset/scattering/elastic.py`）与 S8.4 的
  `amset2d_plugin` 都用 `tensordot(q̂⊗û, D)`（q̂⊗û 不对称化）。物理上 ΔE = Σ D_ij ε_ij，ε = sym(q̂⊗û)·A，D 对称，
  所以 ΔE ∝ Σ D_ij q_i u_j —— D 必须是张量分量。用原生值缩并，剪切部分的耦合振幅大一倍（散射率里那部分 ×4）。
  例：L 谷 LA[111] 应为 Ξd + Ξu，原生 h5 给 Ξd + 5Ξu/3。V121 把非对角抹成 0 时这个问题被掩盖；V188 恢复非对角后必须一起修。
- **验证（第一性原理）**：紧束缚 GaN，沿随机 q̂、横向 û 施加声子应变 sym(q̂⊗û)·η，直接算带边移动斜率：
  原生 h5 缩并最大偏差 **53%**，张量分量 **0.1%**（`test_dp_symmetrize.py` 的 test_v189_contraction_matches_true_band_shift）。
- **改法（数据层，不碰运行期）**：S7.1 在对称化之后调 `dp_symmetrize.to_tensor_h5`，把 `deformation.h5` /
  `deformation_vac.h5` 的非对角原地 ÷2，打 `dp_offdiag = "tensor"`。原版 AMSET、2D 插件、`amset_ir_fix` 的星平均
  （对 D 做 R D Rᵀ，本来就要求张量）都自动变对；S8/S8.4 继续软链 S7.1 的文件，上游失效检测不受影响。
  对角不变：band_edges.json、S8.2 的 E1 只读对角，不受影响。不选运行期补丁：AMSET 的散射子进程是 spawn，
  补丁要靠包装 worker 链式导入（ir_fix / dpo3d / 2D 插件都在抢同一个 worker），脆弱。
- **检查**：`ke_common.deformation_h5_issues` —— S8/S8.4 生成时 h5 没对称化（< 2）或不是张量约定都会告警；
  `tools/regate_projects.py` 记 STALE-DP；`tools/compare_deformation_h5.py` 先把两份 h5 统一到张量约定再比，
  旧 h5 不是张量约定且新场非对角权重 > 2% 时判"需要重跑"（旧运行的剪切耦合口径不同）。
- **h5 不再是 AMSET 原生约定**：拿它喂原版 AMSET 正是想要的（物理正确）；但不要再把它当 `amset deform read`
  的直接输出去和别的 AMSET 结果比非对角。带符号文件 `*_signed.h5` 仍是 AMSET 约定（只给 dp_symmetrize 用）。


## [0050] 冻声子谷间 IVS：K↔K′ 零阶（可选组 `ivs_frozen`，默认关）

- **做什么**：自动算六角 2D 材料 K↔K′ 的零阶谷间形变势 D（eV/Å）与 ħω，S8.4 的 `IVS_FROZEN=auto` 把它变成 dp_optical 的
  `K>K'` 通道（手填的同星 `X>X'` 只作对照；手填的 `X>Y`，如 K↔Q，照用；ADP 谷内掩码随之打开）。用户 2026-10-14 决定默认关，
  材料 tf_*.yaml 写 `ivs_frozen: true` 开。
- **方法（劈裂法）**：√3×√3 超胞让 K、K′ 都折叠到 Γ̃（`ivs_frozen.fold_matrix`，M·k ∈ Z²）。S7.3 算超胞 Γ̃ 点 DFPT，按"原胞平移
  不变子空间"上的投影挑出区边界模（权重 ≈ 0，即原胞 q = K 的声子）；S7.4 对每个模做 ±u 单点。区边界模没有 q = 0 分量，K/K′ 二重态
  的对角元为 0，只被谷间矩阵元 V 劈成 ±|V|u —— 线性、与参考能级无关（同一构型同一 k 的两条能级相减），不要真空对齐。
  Σλ² 取组内去掉共同平移后的迹，简并伙伴模相加：D² = ΣΣλ²/(g(g−1))，g = 2。位移按超胞口径 Δ_j = u·√(M_s/m_j)·ẽ_j，
  超胞口径的 D 就是原胞口径的 D（N_s M_s = N M）；实驻波 C/S 两个模的 |V|² 之和 = 复声子的 |V_A|²。与 dp_optical 的
  `inter_star` 率 πD²/(Mω)·G_星·(g−1)/g 配套。同类做法：Piscanec 等 PRL 93, 185503 (2004) 用超胞冻声子求石墨烯 K 点电声耦合。
- **标来源**：原胞 `ivs-prim`（网格含全部折叠 k）各折叠 k 的能级并起来排序 = 超胞 Γ̃ 的能级（差一个常数），按序配对、中位数对齐
  （残差 > 50 meV 报错），再核对 K/K′ 二重态简并、±10 meV 内没有别的来源的能级；窗口碰到 NBANDS 截断就报错（`IVS_NBANDS_EXTRA`）。
- **验证**：三角晶格紧束缚玩具模型（在位能 α·Δz）——核心函数与整条 gen 流程（S7.3→S7.4→S8.22→S8.4，`tests/test_0050_*`）
  算出的面外组 D = α（相对误差 < 0.2%）、面内组 D = 0，与参考能级平移无关。
- **范围与局限（报数时注明）**：
  - 只做非 SOC（A 方案）：SOC 下 K/K′ 的二重态按自旋劈开（自旋-谷锁定，Xiao 等 PRL 108, 196802 (2012)），谷间散射只能落到
    同自旋子带，劈裂法与星平均态密度都不成立 —— S7.3/S7.4/S8.22 与 S8.4 遇到 SOC 一律拒绝并说明（B 方案另议）。
  - 只算零阶（q 正好在 K）：一阶谷间耦合（∝ |q − K|，Kaasbjerg 等 PRB 85, 115317 (2012) 里 K 点声学支以一阶为主）不在内。
  - 只算 K↔K′：K↔Q、K↔Γ 不在内。S8.22 对离 K 星 < 0.2 eV 的别的星（如 MoS₂ 导带的 Q、价带的 Γ）告警，迁移率仍偏高。
- **代价**：MoS₂ 类 3 原子原胞：9 原子超胞 DFPT 一次 + 18 个区边界模 × 2 + 2 = 38 个超胞/原胞单点。
- **接线**：S1 → S7.3/S7.4 整目录归档（dir）；S7.4 → S8.22、S8.2 → S8.22 一律失效；S7 → S7.3/S7.4、S7.3 → S7.4 不进 DOWNSTREAM，
  由 sha256 / 频率核对报错要求 retry；S8.22 → S8.4 有条件重排（`invalidate_ivs_consumer`：只在用过的参数指纹变了或刚可用时）。
