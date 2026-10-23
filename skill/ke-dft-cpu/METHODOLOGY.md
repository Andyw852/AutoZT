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
     避免重复计数。**谷内 ODP / 谷间 IVS 现在的正路是 DFPT 提取（[0059] eph-qe-cpu 的 `ivs_dfpt`）**；
     冻声子那条自 [0075]（2026-10-23 用户决定）起**已隐藏**：`odp_frozen` 组默认关，只有显式写 `odp_frozen: true`
     的历史材料（当前仅 MoS2）才会注出 S7.2/S8.21。
     细节见 `step8.4_amset2d/README.md`「ODP 与谷间散射」。
     **[0035][0075] 谷内 ODP 的参数（冻声子路，已隐藏）**：显式开 `odp_frozen` 后，S7.2 对 S5 的 Γ 点非极性光学模做 ±u 冻声子单点，
     S8.21 在 S8.2 各谷的 k 点取真空对齐的带边能量算 D = dE/du（位移 Δ_κ = u·√(M/m_κ)·e_κ，与常数形变势模型的
     |g|² = ħD²/(2NMω) 同一归一化），S8.4 自动用上。只有零阶（q 无关）谷内耦合；谷间 IVS 仍要文献参数。
     细节见 `step8.4_amset2d/README.md`「冻声子 ODP 参数」。
     **[0036] 3D（S8）也可选补上**：同一常数形变势模型（公式本来与维数无关：|M|² 的 1/V 与末态求和的 V 抵消，
     G 取每原胞态密度；三维即 Conwell 的非极性光学声子率），插件 `step8_amset/amset_dpo3d.py` 把 ODP/IVS 注册进
     原版三维 AMSET，谷归属用完整三维距离（27 个周期像）。材料 step.conf 写 `DP_OPTICAL_E/H` +
     **必须手填** `DP_OPTICAL_VALLEYS_E/H`（三维没有 S8.2 多谷谷表，k 点写 3 个分量），如 Si 空穴
     `DP_OPTICAL_H = G>G:6.0@63`、`DP_OPTICAL_VALLEYS_H = G:0,0,0`；`ADP_INTRAVALLEY` 规则同 2D。
     transport.json 多出 ODP/IVS 单机制迁移率，S8.3 自动带上。
     **谷内 ODP（冻声子路，[0075] 起隐藏，改走 DFPT）**：显式开 `odp_frozen` 组后 S7.2 在 3D 上照样跑（照抄 S7 undeformed），
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
     **[0061] 更正（2026-10-21）**：原因不止光学声子一个。"ADP 11.6 倍"拿的是只算 ADP 的值比实验总值，不能说明 ADP 本身对；
     AMSET 的 ADP 写法（初态带的对角形变势 × 重叠）还丢了简并价带之间的剪切耦合，Si 只算 ADP 时偏高约 1.6–2 倍（离子弛豫口径），
     已由 `ADP_DEGEN` 修正（见 [0061]）。"静水形变势 a 合理"这一句仍然成立；剪切形变势要看 S7 是不是离子弛豫口径。
     光学声子由 [0035]/[0036] 冻声子或 [0059-r3] DFPT 的谷内 ODP 补上。
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

## [metal_mode] 金属：没有实验数据时怎么算电子输运（2026-10-17 用户批准）

**为什么半导体那条链不能用**：DPT（Bardeen–Shockley，S8.2）和 AMSET 的 ADP 都建立在**非简并带边**上（带边形变势 E1、
有效质量、按掺杂浓度定费米能）。金属没有带边，费米面上是简并电子气，散射由整个费米面上的电声耦合决定。
所以 S7.1 判金属后写 `metal.json`、不写 `band_edges.json`；S8.2 / S8.3 遇到它直接报"不适用"。

**文献里不依赖实验数据的做法（由严到宽）**：

1. **完整第一性原理电声 + 玻尔兹曼方程**（现在的标准做法）：DFPT 算电声矩阵元，Wannier 插值到细网格，
   解 IBTE 或用 SERTA 得到逐态 τ_nk，再算 ρ(T)、σ、κ_e。工具 EPW、Perturbo；综述 Giustino, Rev. Mod. Phys. 89,
   015003 (2017)。代价最高，大胞（如 Mg4C60）很难做；本技能**不包含**这一步（仓库的 `eph-qe-cpu` 只产出电声矩阵元
   数据链，不产出 τ 和输运量）。
2. **输运电声耦合常数 λ_tr + 最低阶变分（Allen）公式**：由输运 Eliashberg 函数 α²_tr F(ω) 得 λ_tr
   （常用近似 λ_tr ≈ λ），高温极限 ħ/τ = 2π λ_tr k_B T，ρ(T) = 1/[(σ/τ)·τ]。Allen, PRB 17, 3725 (1978)；
   用法与金属元素的验证见 Savrasov & Savrasov, PRB 54, 16487 (1996)。只在高温区（T 高于 Debye 温度的几分之一）
   成立，低温（Bloch–Grüneisen 区）会高估散射。**本技能实现的就是这一档**：`METAL_LAMBDA_TR`。
   λ_tr 本身仍要一次电声计算（或引用他人的第一性原理值），`METAL_LAMBDA_SOURCE` 必须写明出处。
3. **常数弛豫时间（CRTA）**：BoltzTraP2 给 σ/τ、κ_e/τ，再取一个**假设**的 τ（文献常用 10 fs 量级）。
   Madsen, Carrete & Verstraete, Comput. Phys. Commun. 231, 140 (2018)。这不是"算出来的"τ，只能给量级；
   论文里必须写明是假设值。本技能的 `METAL_TAU_FS`。

**与 τ 无关、可以直接报的量**（任何一档都一样，S8.1 金属模式总是输出）：Seebeck 系数 S（CRTA 下 τ 抵消）、
σ/τ、κ_e/τ，以及 Lorenz 数 L = κ_e/(σT)（与 Sommerfeld 值 2.44×10⁻⁸ WΩ/K² 比，看 Wiedemann–Franz 是否成立）。
要 κ_e 的绝对值，最省事的组合是 κ_e = L·σ·T，但 σ 仍要 τ。

**本技能的金属模式（S8.1_bt2）**：

- 触发：S7.1 写了 `step7b_deform_read/metal.json`（undeformed 能带 `is_metal`）。不过体系判别闸门、不加剪刀、
  不读 S8.2 的 DPT τ、不做按掺杂的文献口径扫描。
- 取值点：电中性的 E_F(T)（BoltzTraP2 净载流子 = 0 处）。
- 输出：`metal_transport.json` + `boltztrap_summary.txt` 末尾一张表（T、S、σ/τ、κ_e/τ、L，有 τ 时再加 τ、ρ、κ_e）。
- step.conf（本步或材料共用）：
  - `METAL_LAMBDA_TR = 0.6` + `METAL_LAMBDA_SOURCE = "EPW, 2026-xx, 8x8x8 q 网格"` → τ(T) = ħ/(2πλ_tr k_B T)；
  - 或 `METAL_TAU_FS = 10`（假设值）；
  - 两个都不写 → 只出与 τ 无关的量，**不编 τ**；两个都写 → 报错。
- 金属材料在 `tf_<项目>.yaml` 里关掉半导体专用的组：`amset3d / dpt / odp_frozen / ivs_frozen / output_compare: false`
  （S8.2 / S8.3 遇到金属会报错并提示这一行）。S5 介电（DFPT）对金属也没有意义，只被 S8 用；关掉 `amset3d` 后它
  不再被下游需要，已在跑的可以按需停掉（`stop` 先请示）。

## [0071] 2D 单谷 DPT：面内各向异性谷的表头 m* 改用 m_d；面内拟合的 cond 按列归一化（S8.2；2026-10-22 用户批准）

- **表头 m* 偏轻**：2D 表头 m* 取方向分解里"最轻两支"的几何平均（[PATCH-2026-09-23]，为六方 K 谷设计）。面内各向异性谷
  里最轻的两支都在轻轴上（+x 与 −x），取到的就是轻轴质量。磷烯型谷 0.17/1.1 给出 0.170，态密度质量应为
  √(m_x·m_y) = 0.432。2D 表头公式 μ ∝ 1/(m*·m_d) 用的是 m* = m_d，所以 μ 虚高约 (0.432/0.170)² ≈ 6 倍。3D 的同一
  问题已在 [0025] 修过，2D 一直没改。`by_direction` 用的是面内二次型拟合的 m_d，不受影响。
- **改法**：带边 k 点的小群不强制面内质量张量各向同性时（`_site_isotropic(..., dim=2)`：矩形胞的 Γ、带边在 Y–Γ 线上的
  SS/LS；六方胞里带边在 Q 谷的情况本来就被 [V152] 的闸拒绝、交给多谷 DPT，不经过这里），表头 m* 改用 `by_direction` 的 m_d = √(m_x·m_y)，与 [0025] 的 3D 同一口径；`m_provenance` 写明
  旧值。六方 K/Γ、四方 Γ 这类各向同性谷两种做法一致，表头不变（MoS2、CrS2/CrSe2 不受影响）。手填 M_EFF_* 或
  MANUAL_ANISO 质量时也不变。**SS/LS 这类矩形胞材料的表头 μ 会变**：以前走方向分解（偏向轻轴质量）或径向单一曲率拟合
  （接近轻重两向的调和平均），都比几何平均轻；现在用 m_d，表头 μ 一般会下降。
- **cond 误报**：`_quad_fit_2d` 的设计矩阵各列量纲不同（Δk²、Δk、1），不归一化的 cond 随网格加密而变大：带边不在高对称点
  （带线性项）时，48、72 网格上 cond = 1.5e5 / 7.7e5，被当成"病态"直接报错，拿不到数。现在先按列归一化再与 1e5 比较，
  cond 只反映取点几何。只在两条坐标轴上取到点（xy 交叉项无约束）这类真正的病态仍然拒绝。多谷 DPT 共用这个函数，
  一并生效。3D 的拟合没有改：它在实测网格上 cond 都在 1e3 量级，改了会让 [0070] 的窗口选择变化。
- **测试**：`tests/test_0071_2d_headline_md.py` 覆盖以下情形：
  - 面内小群判断（六方 K/Γ、四方 Γ 各向同性；矩形胞 Γ、Q 谷不是）；
  - 磷烯型谷表头取 m_d，μ 与公式一致；
  - 六方 K 谷表头不变；
  - 手填分方向质量时表头不变；
  - 48/72 网格带线性项能出数，主质量误差 < 6%；
  - 十字形取点仍判病态。

## [0070] 3D 带边质量：各向异性谷改用能量窗口取点拟合（S8.2；2026-10-22）

- **错在哪**：[0025] 的三维二次型拟合按同一个球半径 |Δk| ≤ R 取点，再加一个各向同性的 q⁴ 项吸收非抛物性。长椭球谷
  （m_l ≫ m_t）在同一 |Δk| 上，横向能量是纵向的 m_l/m_t 倍（Si ≈ 4.8 倍）：球内横向的点已经离开抛物区很远，纵向的点
  还贴着带边；q⁴ 是各向同性的，分不开两个方向的非抛物性，横向的弯曲就被算进了纵向曲率。Si 实跑（S3 34³）拟出
  m_l = 0.386（实验 0.916）；EPM Si 34³ 同样的做法给出 0.405（EPM 真值 0.873）。m_l 偏轻让 m_c、m_d 都偏轻，
  [0069] 之后的 DPT 电子仍是 4477 cm²/Vs；换文献质量约 2600。
- **文献怎么做**：有效质量是带边极值点的曲率（Ashcroft & Mermin, p. 227），抛物近似只在带边附近成立。三种常规做法：
  1. 带边处有限差分：EMC（Fonari & Sutton）、QuantumATK（五点差分，步长 0.001 Å⁻¹）都要在带边周围另算一组小步长的非自洽 k 点；
  2. DFPT 直接给 Hessian：Laflamme Janssen et al., PRB 93, 205147 (2016)。文中指出有限差分质量需要做步长收敛研究；
  3. 对带边附近的能带做最小二乘拟合：sumo 用抛物与非抛物拟合；effmass（Whalley, JOSS 3, 797 (2018)）让用户选拟合的
     能量范围和多项式阶数，并用 Kane 准线性色散 E(1+αE) = ħ²k²/2m* 量化非抛物性。

  S3 只有均匀网格，能用的是第 3 种；所以按**能量**而不是按 |Δk| 划拟合范围。Kane 色散的偏离只取决于 E：同一个能量窗口
  里，各个方向偏离抛物的程度相同，拟合不会偏向某个方向。Si 导带 α ≈ 0.5 eV⁻¹（Jacoboni & Reggiani, RMP 55, 645
  (1983)），0.15 eV 处偏离抛物约 7.5%。
- **改法**（`_quad_fit_3d_energy`，`MASS_EW_*`）：
  - **谷的小群强制质量张量各向同性**（立方晶体的 Γ、R 等，`_site_isotropic`）时，仍用 [0025] 的球形取点 + q⁴，因为 q⁴
    正好描述各向同性的非抛物性，GaAs 型轻谷用旧做法更准。其余（Si Δ、Ge/PbTe L、Mg₂Si X、非立方晶体的 Γ、一般位置）
    都走能量窗口；
  - 只取离带边 |ΔE| ≤ Ew 的网格点，Ew 从 0.05 eV 逐级放宽到 0.20 eV（`MASS_EW_STEPS_EV`）。拟合二次型；带边没被对称性
    钉住时加线性项和常数项。Δk 取最短等价矢量，所以区边界上的谷（X、L）两侧的点都能取到；
  - **只取本谷的 Voronoi 区域**（`_valley_basin`）：离本谷谷心比离任何其它等价谷谷心都更近的点才算。Si 的 Δ 谷离 X 只有
    0.15·2π/a，X 另一侧属于相邻的等价谷，能量窗口会把它们包进来（EPM Si 不限区域时 m_l 拟合成 6–19）；
  - 每个窗口要满足：点数 ≥ 2 倍参数、三个倒格矢方向都有非零步、cond ≤ 1e5、相对残差 ≤ `MASS_EW_REL` = 0.10（真正的
    抛物区是 0.02–0.04）；
  - **窗口收敛**：与下一个可用窗口的主曲率相差 ≤ `MASS_EW_AGREE` = 15%，才取较小的那个窗口。这对应有限差分的步长收敛
    检查。只剩一个窗口可用时验不了收敛，按失败处理；
  - 失败（网格相对带边曲率太粗）时退回球形取点，并在 `dpt_result.json` 的 `caveats` 里写"m* 不可靠，加密 step3_uniform
    网格后重跑"。`m_provenance` 写明用了哪种取点、哪个窗口、与下一个窗口差多少。
- **通用性**（合成 Kane 能带，fcc 原胞，误差 = 三个主质量的最大相对误差）：

  | 谷 | 真值 m（α, eV⁻¹） | 24³ 旧 / 新 | 34³ 旧 / 新 | 48³ 旧 / 新 |
  |---|---|---|---|---|
  | Si Δ（偏格点） | 0.19/0.92（0.5） | 28% / 退回 | 24% / 6% | 18% / 3% |
  | Mg₂Si 型 X | 0.2/0.5（0.3） | 9% / 退回 | 6% / 2% | 3% / 1% |
  | Ge 型 L | 0.08/1.6（0.5） | 77% / 退回 | 70% / 7% | 59% / 4% |
  | 扁椭球，一般位置 | 1.0/1.0/0.2（0.5） | 23% / 7% | 13% / 5% | 9% / 3% |
  | 三轴，一般位置 | 0.2/0.5/1.5（0.5） | 27% / 7% | 21% / 3% | 15% / 3% |
  | 三轴，强非抛物 | 0.2/0.5/1.5（1.0） | 42% / 18% | 32% / 9% | 26% / 6% |
  | 空穴，各向异性 Γ | 0.3/0.6/1.2（0.5） | 11% / 3% | 6% / 2% | 3% / 2% |
  | 空穴，一般位置 | 0.4/0.4/2.0（0.5） | 28% / 6% | 21% / 3% | 13% / 3% |
  | 纯抛物三轴 | 0.2/0.5/1.5（0） | 0% / 0% | 0% / 0% | 0% / 0% |
  | 立方 Γ 重谷（各向同性，走旧做法） | 0.5（0.3） | 0% | 0% | 0% |
  | GaAs 型 Γ 轻谷（各向同性，走旧做法） | 0.067（0.64） | 11% | 4% | 1% |
  | PbTe 型 L 轻谷 | 0.024/0.24（2） | 99% / 退回 | 66% / 退回 | 67% / 退回 |

  EPM Si（真实能带，34³ 即 Si 实跑的 S3 网格）：能量窗口给 m_t/m_t/m_l = 0.183/0.183/0.854（0.10 eV 窗口），真值
  0.873；48³ 给 0.867。旧做法 0.405。
- **局限**：带边很轻、非抛物性很强的谷（PbTe 型）在常规网格上 0.05 eV 内只有几个点，能量窗口拿不到可靠的拟合，退回旧
  做法并标注 m* 不可靠。这时 DPT 的 m* 要靠加密 S3 网格，或手填文献值 `M_EFF_*`。
- **测试**：`tests/test_0070_mass_energy_window.py` 覆盖以下情形：
  - 合成 Kane 谷（Si Δ、Ge L、扁椭球、三轴、空穴）的精度；
  - 旧做法对长椭球谷的偏差；
  - 纯抛物精确；
  - PbTe 型轻谷拒绝出数；
  - 小群各向同性判断；
  - EPM Si 34³；
  - 端到端：网格够细走能量窗口、粗网格回退并标注、只剩一个窗口时的提示、立方 Γ 空穴仍走旧做法。

  `test_0069` 的 X 谷夹具改用 32³ 网格（12³ 上能量窗口点不够，会被标注）。

## [0069] 3D 多谷 DPT：E1 取谷轴形变势、C 取沿谷轴的纵向模量（S8.2；2026-10-22 用户批准）

- **错在哪**：[0025] 起 3D 分方向的 E1_α 取带边张量的对角；带边谷有对称等价谷时，先对等价谷把张量对角**平均**，再进
  μ ∝ 1/E1²。等价谷的张量互为旋转，平均后只剩静水部分：Si 导带 Δ 谷 D = Ξ_d·I + Ξ_u·ââ，平均后每个方向都是
  Ξ_d + Ξ_u/3。Si 实跑：`band_edges.json` 的 D_xx/D_yy/D_zz = 0.13/0.13/8.72 eV，旧 E1 = 2.99 eV，谷轴值 8.72 eV，
  μ 高估 (8.72/2.99)² ≈ 8.5 倍（DPT 电子 36135 cm²/Vs，实验约 1450）。一个谷里的电子按它自己的耦合散射，散射率
  ∝ ⟨|D:ε|²⟩；先平均张量再平方，就把各向异性的耦合抵消掉了。
- **改法**（`_valley_axis_3d`）：带边谷有对称等价谷、带边不简并时，谷轴 â 取 Γ 指向谷心（最短等价 k）的方向；
  - E1 = |â·D·â|：沿谷轴传播的 LA 声子的形变势。Δ、L 谷就是 Ξ_d + Ξ_u，与文献 Si 电子声学形变势约 9 eV 同一口径
    （Jacoboni & Reggiani, Rev. Mod. Phys. 55, 645 (1983)）；
  - C = C_L(â) = Σ C_ijkl â_i â_j â_k â_l：沿谷轴的纵向模量，取自 step6 全弹性张量。立方 [100] 即 C11，[111] 为
    (C11 + 2C12 + 4C44)/3；
  - 三个方向共用 E1 与 C（τ 各向同性近似），各向异性只来自已对等价谷平均的电导质量 m_α；表头 μ = (μx+μy+μz)/3，
    其中 μ_α ∝ 1/(m_α m_d^{3/2})，不再用 m_d^{5/2}；
  - D 取带边谷自己的张量：优先 S7.1 的带符号中心差分 `deformation_signed.h5`（换成 dE/dε，非对角 ÷2）；没有就用
    `deformation.h5` 的幅值并按 |â| 投影（Ξ_u > 0 且 Ξ_d + Ξ_u/3 ≥ 0 时与带符号结果相同）。
- **不变**：单谷（含立方 Γ 谷，仍是 [0025] 的分方向口径）、简并带边（立方 Γ 价带顶，单带 DPT 只是近似）、2D、手填的
  E1/C/m*。step.conf 写 `DPT_VALLEY_AXIS = off` 复现旧口径。
- **近似**：谷轴 LA 口径是经典多谷 DPT 的常用约定，不是完整的 Herring–Vogt 各向异性弛豫时间：没有显式处理 TA 支
  （Ξ_u 的剪切耦合），也没有处理谷内散射波矢的方向分布。强扁长谷（m_l ≫ m_t，如 Si、Ge）里散射波矢集中在谷轴附近，
  这个口径与 Herring–Vogt 同量级；弱各向异性谷或三轴谷（谷轴不是 D 的主轴）上，它只是一个约定。
- **可信度标注**：`dpt_result.json` 每个载流子加 `caveats` 字段，两种情况会标：带边简并；三维二次型拟合失败、表头 m*
  用了回退值。Si 空穴两条都中（价带顶三重简并，三维拟合残差 rel = 0.44），所以 DPT 空穴值只作量级参考。
- **测试**：`tests/test_0069_dpt_valley_axis.py` 覆盖以下情形：
  - 简立方 X 谷；
  - 面心立方 L 谷，带符号与幅值两路；
  - 纵向模量与 Voigt 重排；
  - 简并带边与单谷保持不变；
  - 表头取分方向平均；
  - 回退值的标注。

  `test_0025` 的等价谷用例改为 `DPT_VALLEY_AXIS = off`。

## [0061] 简并 / 近简并带边的 ADP 带间修正（S8 / S8.4 的 `ADP_DEGEN`；2026-10-21 用户批准，同日通用化）

**问题**：AMSET 的 ADP 矩阵元是 |ê:|D_n(k)||²·|⟨m,k′|n,k⟩|² —— 只用初态那条带的对角形变势（S7 形变构型里按带号排序的能移，
逐分量取 |·|），再乘波函数重叠。带边多重简并（或几条带只差几个 kT）时，应变的主要作用常常是把这几条带**互相耦合**：
应变算符在这个子空间里的非对角矩阵元（立方 p 型价带的 Bir–Pikus b、d；六方的 D3–D6；二维 Γ 点 p_x/p_y 双重态的剪切项；
SOC 下 Γ8 与 Γ7 之间……）。初态与末态是子空间里不同的组合，重叠一乘这部分就没了 → ADP 偏弱、迁移率偏高。只有静水部分
（算符 ∝ 单位阵）时两种写法完全相同，差别全部来自带间耦合。

**定量**（三带 k·p + Bir–Pikus + Christoffel 的模型迁移率，弹性均分、非简并统计，绝对值对 Bardeen–Shockley 解析式 1%；
[0061] 首版的 `model_mobility`，比值部分与通用版逐位一致）。Si 300 K、只算 ADP（cm²/Vs）：

| 剪切形变势 | 完整矩阵元 | AMSET 写法 | 比 |
|---|---|---|---|
| 离子弛豫（\|n\| = 9.46 eV，本地 QE 实测） | 1090–1500 | 2290–2360 | 1.6–2.1 |
| 钳位离子（\|n\| = 2.31 eV） | 2300–2440 | 8000–8300 | 3.3–3.6 |

（范围 = 静水 a 取 ±1.74 eV。）AMSET 实算的 Si 空穴 ADP-only 7610 落在"钳位离子 + AMSET 写法"上 —— 那次 S7 很可能是 [V187]
之前的钳位离子构型：**先核对 S7.1 `band_edges.json` 的 `deform_geometry`**，不是 `ionrelax` 就 retry S7_deform → S7.1。

**约定核对（本地 QE 7.4 PBE Si，SG15，30 Ry）**：k·p L/M/N = −5.57/−3.60/−8.46（Luttinger 换算 −5.64/−3.61/−8.68）；
Γ 点 ε_xx 劈裂给 b = (l−m)/3 = −2.31 eV；离子弛豫 ε_xy 劈裂给 |n| = 9.46，|n|/√3 = 5.46 eV（文献 |d| ≈ 4.85 实验 / 5.3 LDA），
即 H_xy = n ε_xy（张量分量）、n = √3 d；钳位离子 |n| 只有 2.31（Si 的剪切形变势几乎全来自 Kleinman 内应变）；n 与 N 同号。

**做法（`skill/ke-dft-cpu/adp_degen.py`，不预设对称性，全部从本材料的 DFT 数据里读）**：
1. **谷**：DFT 全网格（S3 的 vasprun，与 wavefunction.h5 同源）上带边带的局部极值（离带边 `ADP_DEGEN_EFIT` = 0.15 eV 以内，
   ±2 步网格内比较，免得很平的方向出网格假极值），按点群 + 时间反演归成谷星。谷中心 k0 处与带边带逐条相隔 < 50 meV
   （且大于 k0 处应变造成的能移 2δ·max|S|）的带构成子空间，d 条；d = 1 不需要修正。价带、导带都找。
2. **k·p 矩阵**：k0 附近（|δk| < 0.25 最短（面内）倒格矢，且不超过到别的谷一半的距离）的 DFT 网格点上，用 AMSET 的波函数
   系数（`WavefunctionOverlapCalculator`，网格点上的插值 = 节点值；隔着倒格矢时按 `g_maps` 平移，与 AMSET 的 overlap 同法）
   把各带投影到 k0 处的子空间基，Löwdin 正交化：H_DFT(δk) = Σ_n (E_n − E_edge)|p_n⟩⟨p_n|。与各点的波函数相位、简并子空间里的
   组合方式都无关。投影权重 < 0.8 的点（子空间外的带混进来）不用。
3. **对称性**：先求 {H_DFT(δk)} 的公共不变子空间，可分块就分块（例：二维 p_z 与 p_x/p_y 在面内 k 下不耦合；Kramers 对处处
   ∝ 单位阵 -> 一维块，不需要修正）。每块对 k0 小群的每个操作（点群操作 R 与时间反演 ∘R）解 H(Mδk) = U H(δk)[*] U† 的 U
   （零空间唯一到相位才用），闭合成群，求系数的对称允许子空间（群平均投影）。Si：96 个操作，k·p 3 个参数（L、M、N）、
   应变 3 个（l、m、n）——就是 DKK + Bir–Pikus。
4. **k·p 模型** H(δk) = D0 + Σ δk_a B_a + Σ δk_a δk_b Q_ab（只在对称允许的子空间里）：网格点 H_DFT 按载流子能量加权线性拟合，
   再用 AMSET 自己的（插值后更密的）能带做本征值最小二乘精修（网格点锚定，窗口不够点放宽 2、4 倍）；相对残差 ≤ 15%。
5. **应变算符** H_ε = Σ_ij ε_ij Y_ij（同样只取对称允许的部分：Bir–Pikus 原理，ε_ij 与 δk_i δk_j 同样变换）：
   (a) k0 点严格重演 S7.1：±δ 应变 → 按能量排序 → 能移/δ，对上带符号中心差分 S（`*_signed.h5`）与 AMSET 幅值 A（运行期 h5）；
   (b) k0 附近的 S7 网格点：不交换次序的带（能隙 > 1.5·δ·|ΔS|）单独用 ⟨c_n|X^(P)|c_n⟩ − Σ P_ab δk_b ∂_aE_n = −S_n,P
       （固定分数坐标的 k 随应变移动），会交换的取组内之和；离 k0 越远、越靠近交换阈值（三阶误差 (δ|X|/能隙)²）权重越小；
   多起点 Levenberg–Marquardt（纯 numpy）。k0 点的相对残差 ≤ 20%，雅可比满秩。数据分不开的几个解（S7 太粗，Γ 附近没有
   非简并点时立方 n 的符号）：各解的放大系数差 ≤ 10% 取平均，否则不修正并报告"加密 S7 网格"——不靠猜。
6. **比值表**：同一个模型里，每个载流子能量（对数网格 5 meV – efit）、方向（3D 240 / 2D 96）、带：
   r = Σ_λ∫|M_full|²/C_λ·w ÷ Σ_λ∫(ε_λ:D_AMSET)²|⟨m,k′|n,k⟩|²/C_λ·w（等能面积分；3D 三支声学模，2D 只取面内 2×2 Christoffel
   子块，与 amset2d_plugin 同口径）。两边**照 AMSET 弹性散射的实际算法**（[0061-g] 起；首版两边都只用 H_ε、总散射率）：
   - 末态权重 w 含动量弛豫因子 1 − v·v′/|v|²（AMSET `calculate_rate` 对弹性机制总是乘）；
   - D_AMSET = |Y_n − sym(δk⊗∇E_n)|·(非对角约定跟 `dp_offdiag`) + m_e v⊗v：S7 文件里是固定分数坐标的形变势，运行期取 |·| 再加 v⊗v；
   - M_full = ⟨m,k′|H_ε|n,k⟩ + ε:(m_e v⊗v − sym(δk⊗∇E))_n·⟨m,k′|n,k⟩（运动学项与 v⊗v 按初态带内项，与 AMSET 同口径）。
   这样比值只差在"带间应变耦合"（以及子空间里不再对形变势取 |·|）这一处。Si 上这两项让各带平均放大系数变几个百分点
   （首版 1.93/3.21/2.26 -> 0.01 eV 处 1.89/3.11/2.25、0.1 eV 处 1.81/2.78/2.04）；二维双重态影响更大（动量弛豫因子对带间
   跳到更快的带是负的，个别态放大系数 < 1）——这是 AMSET 动量弛豫近似本身的性质，比值照它算才能乘到 AMSET 的率上。
7. **运行期**：插件（3D `amset_dpo3d.py`，2D `amset2d_plugin.py`）的 ADP 类在主进程 `from_amset_data` 里调
   `adp_degen.runtime_tables` 拟合一次，表随 reference 带进散射子进程。
   - **逐对替换矩阵元（[0061-p]，出厂 `ADP_DEGEN_PAIR = auto`）**：在 AMSET 自己的四面体积分里换。AMSET 的 `calculate_rate`
     （0.4.19 / 0.5.1 同构）对每个初态依次算 `tbs.get_masks(spin, tet_mask)`（tet_mask[0] = 各四面体的末态带号）、
     `tet_overlap = get_cross_section_values(vert_overlap, …)`（截面上的重叠）、`get_fine_mesh_qpoints`（细 q 点 → 四面体），
     然后才调各机制的 `factor()`，再乘 `tet_overlap[mapping] * weights * mrta_factor`。插件包装这三个函数记下中间量（导入时装，
     散射子进程里也装；源码结构对不上就不装、退回比值并告警），ADP 的 `factor()` 对子空间内的 (初态, 末态) 返回
     F_full / 该点的 AMSET 重叠 —— AMSET 随后乘的正是这个重叠，精确抵消，留下
     F_full = Σ_λ |⟨m,k′|H_ε|n,k⟩ + ε:(m_e v⊗v − sym(δk⊗v))_n⟨m,k′|n,k⟩|² / C_λ（v 用 AMSET 自己的速度，2D 乘 `ADP_D_SCALE`）。
     初态、末态的子空间分量来自 k·p 模型（谷星成员按对称操作映回代表谷），四面体权重、动量弛豫因子、态密度、末态集合全是 AMSET 的。
     子空间外的末态、离谷 rcut 以外的点照 AMSET 原样；rcut 前 20% 线性过渡；AMSET 重叠 < 1e-14 的点无法用"除以重叠"表达，照原样。
   - **按初态比值**（`ADP_DEGEN_PAIR = off` 或钩子装不上）：在 `_iv_keep` 里乘第 6 步的比值（`amset_ir_fix` 星平均之后，同谷内掩码）。
   报告 `adp_degen_fit.json`（`method` 写明用的哪种；每个谷：子空间带数、分块、对称操作数、参数数、两段残差、各带平均放大系数——
   逐对模式下这是按初态比值口径的诊断值），S8.3 汇总每个修正过的谷一行 `[ADP 简并修正]`。报告 `adp_degen_fit.json`（每个谷：子空间带数、分块、对称操作数、参数数、两段残差、
   各带平均放大系数），S8.3 汇总每个修正过的谷一行 `[ADP 简并修正]`。

**gen**：`ADP_DEGEN = auto`（出厂：有简并谷就修正，没有/数据定不住照原版）/ `on`（不适用报错）/ `off`。gen 只判开关、
SCATTERING 有没有 ADP、有没有 S7.1 的带符号文件（与本步实际用的 h5 同口径：S8.4 真空口径配 `deformation_vac_signed.h5`），
S8.4 还要求真实重叠（`UNITY_OVERLAP = true` 时做不了）；`deform_geometry` 不是 ionrelax 时告警。

**合成数据验证**（`tests/test_0061_adp_degen.py`：已知 k·p + 应变算符生成能量、随机规范的波函数系数、AMSET 密网格、S7.1 文件，
比值表对上真模型）：立方 Γ 三重态（误差 < 1e-3，与 [0061] 首版立方实现逐位一致）、Si + SOC 的 Γ8 + Γ7（6 条带，< 2%）、二维六方
Γ 双重态与 K/K′ 谷星、正交晶系 20 meV 近简并、二维 p_z + p_x/p_y 自动分块、Kramers 对与不简并谷不修正、符号分不开时取平均；
4×4 真模型（远带让三带 k·p 残差约 10%）时各带平均放大系数仍在 4% 以内。

**近似（相对 AMSET 自己的 ADP 而言）**：
- 弹性、均分：AMSET 的 ADP 本来就是弹性 + 均分（率 ∝ k_BT），比值用同一口径，不额外引入；低温（≲ 30–50 K）要非弹性声学散射，
  不在 AMSET ADP 的范围里。
- 只算谷内：形变势（长波应变）理论只对小 q 成立，谷间大 q 本来就该由 IVS（[0059] DFPT / 冻声子）给，AMSET 原版对大 q 用长波
  公式才是不物理的那部分（`ADP_INTRAVALLEY` 掩码）。
- "按初态乘比值"的近似（比值在 k·p 模型的等能面 / 态密度上算，乘到 AMSET 在 DFT 插值能带上的积分里；合成数据里远带让 k·p
  残差 10% 时 ≤ 4%）在逐对模式下没有了：末态集合、权重、态密度、动量弛豫都是 AMSET 自己的，只有矩阵元换成完整的。
  剩下的近似只是初态 / 末态在子空间里的分量取自 k·p 模型（谷附近成立，rcut 外不修正）。逐对与比值两条路在同一组末态上逐位
  一致（测试 `test_pairwise_matches_ratio_table`，3D Si 与 2D 双重态，相对 1e-8；把运动学项符号改错就对不上）。

**局限**：
- 谷中心要落在 DFT（S3）网格上（K 点要 3 的倍数网格；落不上时最近的网格点上简并已被线性项劈开，常常定不住，报告写明）；
- k·p 与"应变算符不随 k 变"只在谷附近成立（截断外不修正）；与 k 有关的形变势、非抛物性进残差，只经由 AMSET 自己的能带/态密度进入；
- 主进程一次性开销：3D d = 3 约 20 s / 谷，SOC d = 6 约 1–2 min（240 个方向）；要读一遍 vasprun.xml；
- 弹性用 Christoffel、均分近似，与 AMSET ADP 核同口径；谷间散射不在这里（IVS / ADP 谷内掩码另管）。

### [0061-m] 应变算符直接取 DFT 波函数：形变势矩阵（方案 A1；2026-10-21 用户选定）

**为什么**：上面第 5 步从 S7 的**本征值**反推应变算符，简并点上本征值只给"劈裂"，非对角元的大小与符号要靠 k·p 模型 +
邻近非简并点去猜（S7 太粗时定不住、只能取平均）。A1 把形变构型的**波函数**投影回未形变的能带基，应变一阶哈密顿量
X^(P)_mn(k) = ⟨m,k|∂H/∂ε_P|n,k⟩（固定分数坐标，即约化坐标下的应变 DFPT 一阶哈密顿量，Hamann–Wu–Rabe–Vanderbilt
PRB 71, 035117 (2005)）在能带基上的**整个矩阵**直接出来，没有歧义解。
**说明**："把 DFT 形变势矩阵接进 AMSET 的四面体积分"这个组合**没有文献先例**，是本项目自己的实现（单独的部分各有出处：
Bir & Pikus 1974；Khan & Allen PRB 29, 3341 (1984)；Li, Graziosi, Neophytou PRB 104, 195201 (2021) 的 Si 价带 3×3 声学耦合矩阵；
Fang et al. PRB 98, 075106 (2018) 的应变 Wannier 哈密顿量）。

**三段（`skill/ke-dft-cpu/deform_wave.py`，只依赖 numpy + h5py）**：
1. **S7 作业里**（`gen_step9_deform.py`，step.conf `DEFORM_WAVE = true` 出厂）：只在"本征值取自它"的那次计算写 WAVECAR ——
   ionrelax 的静态段（`INCAR.static` 加 `LWAVE`）、没有 ionrelax 的顶层单点（undeformed；3D 且 `IONRELAX_3D = false` 的形变目录）；
   2D 面外的离子固定构型不要。**不改模板的 LWAVE**：`run_incar_mismatch` 把 LWAVE .F.→.T. 当"输入变了"，会整批重算；LWAVE
   不算进 ionrelax 的复用核对。已完成但缺 `deform_wave.h5` 的 ionrelax 只补静态段（`prepare_static_only`，归档后缀 `nowave`）；
   undeformed 的顶层 INCAR 变了，V168 会归档重算这一个单点。submit.sh 末尾 VASP 正常收尾才 `deform_wave.py extract`：
   带边 ± 0.8 eV 内有态的 k 点、能带窗口（带边 ± 1.2 eV 的带，两侧各多抽 `DEFORM_WAVE_MARGIN` = 6 条）、系数截到 300 eV，
   写 `deform_wave.h5` 后删 WAVECAR（`DEFORM_WAVE_KEEP = true` 不删）；抽取失败只告警。
2. **S7.1**（`gen_step9b_deform_read.py`）：`amset deform read` 之后 `deform_wave.py build` → `step7b_deform_read/deformation_matrix.h5`。
   窗口只由未形变构型定（在各形变构型两侧仍剩 2 条余量的范围内，把窗口边界扩到不切开近简并组；不再各构型各自按能量扩 ——
   导带到处交叉，应变把简并劈得各构型不一样，会连锁缩窗口）。每个 k 点两组系数各自 Löwdin 正交化，
   H^s = S^-½ O diag(E^s) O† S^-½（O = ⟨未形变|形变⟩），X^(P) = (H^+ − H^−)/(δ+ − δ−)。存未对齐的 X、未形变的正交化系数、
   覆盖度 q_min 与赝波函数非正交度。缺文件 / build 失败只告警（S8 退回第 5 步的拟合）。带符号 h5 同时存 `ibands_<spin>`。
   **AMSET 按对称性约化掉的图样与一侧**（立方 Si 只算 xx± 与 yz+）：对缺的每个图样、每一侧，从 undeformed/POSCAR 找出把该图样
   找把某个真算构型的应变 a·E_P 变成它的空间群操作 g（x → Wx + t，应变变成 a·RᵀE_P R：yy/zz 由 xx 转来，xy/xz 由 yz 转来，
   另一侧由镜面），新晶体的态 = 来源构型在 k W⁻¹（或其时间反演）处的态经 g 变换（系数乘 e^{2πiκ·t}，波矢 κ → κW），
   全部中心差分；某个 k 点两侧的像都没抽到记为缺格点；找不到操作（非共线、没有 POSCAR）才退回前向差分（O(δ)）并点名。
   **AMSET 自己的对称补应变标签错（2026-10-22 查实并改正，影响 AMSET 基线形变势，不只是本节）**：`get_symmetrized_strain_mapping`
   把构型 E 的能带搬到 k' = M·k（M = spglib 的 Wᵀ），这组能带属于 Rᵀ·E·R，它却标成 R·E·Rᵀ（`strain.rotate(r_cart)`）。对合操作
   （镜面、2 次轴、反演）两者相同；3/4/6 次轴不同 —— 立方 Si 的 C3[111] 把 yz 补出来的 xy、xz 互换（AMSET 0.4.19 / 0.5.1 源码在 EPM Si
   上重演：xy+、xz± 能带与标签差 0.22 eV）。S7.1 后面的 Laue 群对称化平均不掉它：EPM Si 走完整个 S7.1 流程（AMSET 幅值 + 带符号
   + `dp_symmetrize.combine`），剪切分量的 BZ 平均 |D| 价带顶 3.45 eV（12 个应变全算 7.23）、导带底 1.37（5.04），对角分量不受影响。
   Si 实跑的信号：带符号文件的 xy、xz 与本节的 X 对齐离散度 2.7 / 3.9 eV（`tools/adp_k0_probe.py` 的 1b 节能认出"S_xz 实为 X_xy"）。
   改正：`dp_symmetrize.symmetrized_strain_mapping`（照搬 AMSET 的循环，只把标签换成 strain.rotate(Rᵀ)），S7.1 的
   `amset deform read` 包装与 2D 真空参考路径都用它，日志 `[strain_label_fix] 对称补了 N 个应变，其中 M 个 AMSET 原版会标成别的应变`；
   改后与 12 个全算逐位一致。`BAND_EDGES_REV` 随之改版。受影响的：立方（剪切 xy/xz 互换；视 AMSET 取到的操作，yy/zz 也可能）、
   四方（C4 把 xz/yz 的正负号换掉）——要重跑 S7.1（只是重读，不重算 S7）；六方 / 三方 / 正交及更低对称只靠对合操作补出单分量应变，
   M = 0、结果不变（MoS₂ P-6m2、SnSe Pnma 逐个操作核过）。
3. **S8 / S8.4 运行期**（`adp_degen.py`，`ADP_DEGEN_MATRIX = auto` 出厂，有 `deformation_matrix.h5` 才用）：
   - **对齐**：X 里混着参考能级（芯势 / 真空）的应变导数。非简并态的 X_nn + S_n = dref/dδ（S = 带符号形变势），逐图样取中位数
     作常数 c_P；离散度 > max(0.1 eV, 5% × 对角形变势中位数) 就不用（带号错位一条时离散度立刻变大）。
   - **网格算符**：在 S7 网格上对**算符本身**三线性插值 D̂(k) = Σ_i α_i Σ_ab |a,k_i⟩X_ab⟨b,k_i|（与规范无关）；ISYM = 0 只给半个
     网格，缺的格点用时间反演补（SOC 的旋量 (c↑, c↓) → (c↓*, −c↑*)），跨布里渊区边界按 Miller 平移。作用在 AMSET 的 Miller 列表上。
   - **k0（谷中心）**：Z = Bᴴ D̂(k0) B（B = AMSET 在 k0 的子空间系数 × 分块酉阵 W）。先按 B 在网格基上的投影范数²（≥ 0.9）认出
     AMSET 系数用的是 k0 的哪个倒格矢表示；分块之间的元素置零；能投影到对称允许子空间（残差 ≤ 10%）就用投影后的；k0 在 S7 网格上
     时与 S 核对（与拟合同一口径，≤ 20%）。任一不过 → 退回拟合，报告 `matrix_fallback` 写原因。
   - **逐对替换**：rcut 以内 AMSET 网格点的 Zk = Bᴴ D̂(k_r) B 在主进程预算（散射子进程只拿小数组）；它已含固定分数坐标的
     运动学项（带内 = 拟合模式的 −ε:T，外加带间部分），所以不再另减 T（测试：k0 处 Zk = Z，离开 k0 对角元 = Z 的对角元 − T）。
   - **网格公式**（核心区外、有简并谷的那种载流子的全部初态、全部末态）：AMSET 缓存了波函数系数时，在它的四面体积分里对每个末态
     顶点算 A = ⟨f|[D̂^(P)(k)|i⟩…, |i⟩]（末态按 g_diff 平移，与 0.5.1 `_get_overlap` 逐字相同），顶点双线性量 Re(A_μA_ν*) 用 AMSET
     原来的 `get_cross_section_values` 插到截面上，F = Σ_λ e_λᵀ G e_λ / C_λ，e_λ = [ε_λ 的图样分量（剪切取 ε_ab，不乘 2）/Hartree,
     ε_λ:(m_e v⊗v)]；返回 F / AMSET 截面重叠（再乘谷内掩码）。|⟨f|i⟩|² 那一项与 AMSET 的截面重叠逐位核对，不过就整体关掉网格公式。
     核心区 min(rcut, 3 个 DFT 网格间距) 里子空间内的末态照用 k·p 态（AMSET 的系数线性插值在简并点附近不可信），外 20% 线性过渡。
     StrainGrid 的大数组经 `amset.util.create_shared_array` 进子进程（complex64 按 float32 存）。

**验证（`tests/test_0061m_*.py`，EPM 测试台 `tests/_epm_si.py`：金刚石结构经验赝势，∂H/∂ε 直接对平面波矩阵求导 = 精确值）**：
WAVECAR 写读逐位一致（真实 VASP 文件的读法另对照过 pymatgen：标准 / 自旋极化 / 非共线）；投影 X 与精确值 < 1%（含 Γ 三重态
非对角元与符号）；对齐读回混进去的参考常数，带号错位离散度变大；网格算符离格点误差随网格加密二阶下降（6³ 价带中位 3%，
Si 实际 S7 网格比它细约 1.8 倍）；Γ 三重态 Bᴴ D̂ B 与精确 ∂H/∂ε < 1%；合成 k·p 模型 + 精确算符时比值表对上真模型、
拟合分不开符号的粗 S7 直接定下来；网格公式与"逐顶点求和再插值"逐位相等，自检能抓住倒格矢表示错误。

**近似与风险（汇报时要讲）**：
- ADP 框架本身（长波应变、弹性、均分）与 AMSET 原版相同；只有 EPW 那条路（方案 B，未选）能去掉。
- PAW 赝波函数系数（不含补偿项，与 AMSET 的波函数重叠同口径）；`build` 报告赝波函数非正交度。
- 中心差分 O(δ²)（δ = 0.5%）；截断球随应变变（Pulay），S7 ENCUT = 1.5 × ENMAX 时较小。
- 窗口截断：D̂ 只含窗口里的带，远带贡献 O(q·v/ΔE)；EPM 上窗口从 3³/1-7 带到 6³/0-9 带误差降一个量级。
- 网格插值：S7 网格越粗误差越大（3D 默认 vaspkit 网格，Si 11³；设了 `DK_MAX_3D` 才与 S3 相同）。
- AMSET 自己的系数插值在简并点附近不可信 → 核心区仍用 k·p 态；核心区半径取 3 个 DFT 网格间距（经验值）。
- 成本：每个形变构型多一个静态 SCF 的 WAVECAR（作业里抽完就删）；`deform_wave.h5` 只含带边附近的 k 点，Si 每个几十 MB，
  2D 大真空胞会到几百 MB（S7.1 的 `deformation_matrix.h5` 也随 fetch_all 回拉本地）。

## [0059] DFPT 直接算谷间零阶形变势（eph-qe-cpu 可选组 `ivs_dfpt`；S8 / S8.4 的 `IVS_DFPT`）（2026-10-18 用户批准）

**决定**（方案 `docs/PLAN-ivs-dfpt-20261018.md` 的待确认 A–D）：A 放进 eph-qe-cpu（QE 环境在那里），ke-dft-cpu 的 S8/S8.4 读它的产物；
B 同意新增步骤；C 与 DFPT 基准差 20–30% 以内算通过；D 用 3090、先拿 Si 测。这条取代 [0058] 里"流程本身不再扩展谷间散射的自动计算"
那一句：[0058] 的文献参数入口（`DP_OPTICAL_E/H`）照旧可用，有 DFPT 结果时手填谷间通道只作对照。

- **四步**（eph-qe-cpu，`ivs_dfpt: true` 开；只算 IVS 的材料再写 `perturbo: false` 关掉 S0–S6）：
  S10_ivsq（谷表 → 谷对 → 对称等价类 → q 作业）→ S11_ivsscf（环境探针 + pw.x SCF，VASP 结构原样 + 谷能差自检）→
  S12_ivsph（每个 q 一个作业：ph.x 单 q DFPT → pw.x nscf（初态 k）→ ph.x `electron_phonon='ahc'`）→ S13_ivs（D0、ħω → 星通道，
  `ivs_dfpt_params.json`，schema `ivs_dfpt/1`）。数学在公共池 `skill/_common/ivs/ivs_dfpt.py`（生产、消费两端共用）。
- **谷表**：2D 读本材料 ke-dft-cpu S8.2 的 `multi_valley.valleys`（用 `k_grid` 网格极值点，与 S8.4 的谷归属一致；重复星按
  `merge_coincident` 合并）；3D 在 eph-qe-cpu 的 step.conf 手填 `IVS_VALLEYS_E/H`，每个星写一个成员即可，按晶体点群 + 时间反演补全
  （Si：`X:0,0.425,0.425` → 6 个）。S8.2 的星成员不全（按对称性应更多）就拒绝：S8.4 的态密度按谷表成员数平均，补了会对不上。
- **QE 接口（M0 核实，本机 conda QE 7.0 与 7.4 结果逐位相同；7.3/7.4 的 `PHonon/PH/ahc.f90` 输出部分只差 GPU 指令）**：
  `ahc_dir/ahc_gkk_iq1.bin` 直接存取、无记录头，每个 k 一条记录，Fortran 次序 `ahc_gkk(ib, jb, imode)` =
  ⟨ψ_ib(k+q)|∂V_SCF/∂u_{q,imode}|ψ_jb(k)⟩（Ry/bohr；ib 是 k+q 的全部带，jb 是 AHC 窗口里 k 的带，imode = 笛卡尔位移）；
  `trans=.false.` 时读前一次 DFPT 存下的位移模式（`u_from_file`），`setup_nscf` 在 AHC 模式下跳过 `irreducible_BZ`，
  nscf 给的 k 原样保留。只支持模守恒赝势（USPP/PAW/DFPT+U/磁性/2D 库仑截断未实现或未测试）。
  核对 ①：q = Γ 时 ⟨n|∂V/∂u_{2x}|n⟩ 与冻声子有限差分 dε/du 四位有效数字一致（Si 一般 k 点，带 3–6：0.00963/0.05539/0.24525/−0.20993
  vs 0.00962/0.05537/0.24520/−0.20987 eV/Å）—— 量纲与按原胞归一确认；② 动力学矩阵按 fildyn 的 (na, nb) 块、行 = na 的方向解析，
  无反演中心的结构上本征矢与 QE `dynmat.x` 重叠 = 1（取共轭只有 0.46–0.95）—— 与 postahc.x 的投影 g_ν = Σ_j g_j u_jν 同一口径。
- **形变势**：g_ν(m,n) = Σ_j ahc_gkk(m,n,j)·e_jν/√M_j，D_ν² = M_cell·Σ_{n∈初态简并组, m∈终态简并组}|g_ν|²/(d_初 d_终)，
  与 dp_optical / ivs_frozen 的 |g|² = ħD²/(2NMω) 同一口径；声子简并组按 D² 相加当一支。规范不变（测试里对简并子空间做随机幺正变换不变）。
- **星通道**：初态取星 X 的第一个成员（等价成员率相同），终态遍历所有星的所有成员：`X>Y` 每条 D² 除以 g_Y，`X>X'` 除以 g_X − 1，
  同一对星 ħω 相差 < 0.1 meV 的按 D² 合并。对称性约化：(k_a, k_b) 与 (S k_a, S k_b)、(−k_a, −k_b)、(k_b, k_a) 同类（后两条来自时间反演
  与厄米性），每类一个 ph.x 作业；MoS₂ 电子 6 类（K↔K′、K↔Q 近/远、Q↔Q′ 三种），Si 电子 2 类（g、f）。
- **Si 验证（M0，SG15 PBE、30 Ry、6³ k，冒烟精度）**：g 过程（q = 0.3X）只有 LO（63.1 meV）耦合，D = 4.79 eV/Å，TA/LA = 0；
  f 过程 48.5 meV 2.73、58.4 meV 4.42 eV/Å，TA = 0（零阶选择定则）。ħω 与 Jacoboni & Reggiani RMP 55, 645 (1983) 的 g/f 声子
  （12/18.5/62、19/47/59 meV）一致；D 与其拟合值（g-LO 11、f-LA 2、f-TO 2 × 10⁸ eV/cm）同量级 —— 拟合值吸收了一阶过程，不要求逐条相等。
  折成星通道：`X>X':2.141@63.06; X>X':2.443@48.48; X>X':3.953@58.40`（f 4 个终态 × D²/5，g 1 个 × D²/5）。
- **S8 / S8.4 接线**：autozt 按 push_paths 把本地 `eph-qe-cpu/result/step13_ivs_read/ivs_dfpt_params.json` 推到 ke 材料目录的
  `ivs_dfpt/`（本地没有就删掉远端旧的）；`IVS_DFPT = auto`（出厂）有就用。优先级 DFPT > 冻声子 K↔K′（`IVS_FROZEN`）> 手填：
  用上时冻声子 `K>K'` 与手填 `X>X'`、`X>Y` 只打印对照、记进 `dp_optical_ivs_dfpt.info`，不进散射率（冻声子记 `superseded_dfpt`，
  参数再变也不让 S8.4 重排）。谷表按 k 中心配对星（`IVS_DFPT_KTOL`：S8 0.08，粗网格上自动找的谷中心可差半个网格步长；S8.4 0.005，2D 两边都是 S8.2 的网格点；DFPT 与本步谷标签可以不同）；对不上（S8.2 谷表变了、
  手填谷表不同）报错。3D 没有手填谷表、也没有冻声子 ODP 谷表时直接用 DFPT 的谷表。`ADP_INTRAVALLEY=auto` 随之把 ADP 收成谷内。
  记录：S8.4 的 `2d_correction.json`、S8 的 `dp_optical_3d.json` 都写 `dp_optical_ivs_dfpt`（mode、sha256、逐载流子 used/指纹/摘要）；
  S8.3 的 comparison_summary.txt 写一行 `[IVS 参数来源]`。参数是跨技能来的，S13 不会自动让 S8/S8.4 重排：结果拉回本地后对
  S8_kappa / S8.4_amset2d 执行 retry（gen 时推过去）。SOC 能带一律拒绝。
- **范围与局限**：零阶、常数 D（不含一阶 K↔K′，Filipovic Table 5：缺它 μ 偏高约 8%）；非 SOC（W 化合物 Q 谷的 SOC 劈裂大，WSe₂ 电子
  可信度最低，见 [0058]）；K 谷 ADP 仍用真空 E1（8.54 vs DFPT 6.5）。验收（用户定 C）：与 DFPT 基准差 20–30% 以内算通过。
- **3090 部署**（铁律 11：装到 `~/software/AutoZT/`）：QE 7.3/7.4（AHC 自 6.7 起；conda-forge 有 7.4：
  `micromamba create -p ~/software/AutoZT/qe-7.4 -c conda-forge qe=7.4`）+ 模守恒赝势目录（PseudoDojo/SG15）；
  `setting/3090.yaml` 写 `qe_bin` / `qe_pseudo_dir`（可选 `qe_setup`），autozt 注入 step.conf 的 `QE_BIN` / `QE_PSEUDO_DIR` / `QE_SETUP`。
  通用提交模板 `submit_ivs_qe.tpl` 按 `IVS_NP` 申请 `--ntasks`（fakeslurm 同样按它限核）；S11 的作业先跑环境探针（pw.x/ph.x、
  ph.x 是否支持 AHC、赝势是不是模守恒），缺东西直接失败并写 `ivs_env.json`。
- **Si 集群实测（2026-10-19，hfeshell + jzzn）**：完整记录见 `docs/RESULT-ivs-dfpt-Si-20261019.md`。要点：QE 7.4 + SG15
  `Si_ONCV_PBE-1.2.upf` + pbesol + 60 Ry + 10³ k；g 过程（q-01）只有 61.4 meV LO 耦合、D=4.70 eV/Å，f 过程（q-02）47.1/56.9 meV、
  D=2.61/4.28 eV/Å，TA/LA 严格为 0；三条通道的 D/ħω 都落在文献散布范围内（逐谷对 D 与 Jacoboni/Pop 的 ∆ **直接比**：
  Γ = π∆²Z_f/(2ρω)·g_c 与 πD²G₁/(Mω) 在 D = ∆ 时相等；初稿的 D = ∆·√(Z_f/2) 作废）。ke S8（AMSET 0.5.1）低掺杂电子
  μ(300 K)=1459（ADP-only 3883），与 Si 电子晶格限制迁移率实验值约 1400–1450 一致，IVS 单机制 2747 为限制项。
  方案里"ADP-only 2381 → ~1400"的 2381 是 `postprocess_intrinsic.py` 日志取到的 −1e21 档，不能作基线。
  IVS 单机制的掺杂依赖经抛物带模型核对属物理（简并后发射通道打开）。
- **遗留项处理（[0059-r2]，2026-10-19）**：日志"本征迁移率"改取最低非零掺杂档（`display_doping_index`，原来取索引 0 = −1e21）；
  S8.3 汇总对 AMSET 的"机制关闭"极大值（>1e5 cm²/Vs）印 "—" 且不参与限制因素；`odp_frozen.find_valleys_3d` 的平台合并
  不再跨 k 与 −k（Si ±Δ 谷在 11³ 网格上相邻且简并，旧版 6 谷并成 V1×3）。恰好 0 掺杂档仍是 AMSET 自身行为，只记录不改。
  测试：`tests/test_0060_ivs_fixes.py`。
- **DFPT 顺带算谷内 ODP（[0059-r3]，2026-10-20）**：eph-qe-cpu `IVS_ODP = auto`（出厂开）在 q 列表末尾加一个 q=Γ 作业 `q-odp`，
  每个谷取星的第一个成员做初态，同一套 AHC 给出谷内零阶光学形变势（kind `intra`，D 不除星成员数）。
  - **单谷载流子也算**：开 ODP 时只有一个谷的载流子（Si 空穴 Γ）保留（只算谷内），3D 要手填 `IVS_VALLEYS_H = G:0,0,0`；
    `IVS_ODP = off` 时照旧跳过。
  - **AHC 带窗按 q 放宽**：S12 读 S11 的谷能带（`valley_bands.xml`），某个 q 的初态带边简并组多于 `IVS_AHC_NB` 时只放宽这个 q
    （Si 空穴 Γ 三重简并 → 上下各 3 条），已算完的谷间 q 不重算。S11 的谷能带里没有新 q 的初态 k（S10 后来加了载流子/谷）就报错，
    要求先重跑 S11（谷能差自检要覆盖新载流子）。
  - **Γ 的声学模显式剔除**：按平移子空间权重取 3 个声学模，不进谷内通道。QE 在 Γ 给的声学频率不严格为 0
    （本机 Si 0.82 meV，r3 初版集群上 0.58 meV），只靠 ħω < 0.5 meV 的阈值会漏进去。
  - **极性模**：`IVS_ODP_POLAR = exclude`（出厂，与冻声子 `ODP_POLAR` 同口径，长程 Fröhlich 由 POP 管）时 q-odp 的 DFPT
    加 `epsil`，读 dyn 里的 Born 有效电荷，模有效电荷 > `IVS_ODP_POLAR_TOL`（0.1 e；2D 只看面内）的模不进通道；单元素体系按
    非极性处理，不加 epsil。`include` 则极性模的短程形变势也算（不加 epsil，零带隙体系也能跑）。
  - **消费端优先级**：S8 / S8.4 的谷内 ODP 按 DFPT > 冻声子 > 手填。DFPT 算了某载流子的谷内（`raw` 里有 intra 类；
    全部 < `IVS_DFPT_DMIN` 也算，结论就是 ≈ 0），冻声子/手填的谷内通道只作对照（`ivs_dfpt.odp_reference`），
    冻声子记 `frozen_use = superseded_dfpt`、指纹去掉（S8.21 参数再变也不让 S8/S8.4 重排）。r3 初版两份都进了散射率
    （Si 空穴冻声子 7.77 + DFPT 6.0 同时出现），本次修正。记录 `dp_optical_ivs_dfpt.odp_use`；S8.3 汇总多一行 `[ODP 参数来源]`。
  - **本机 QE 7.4 实测**（Si，SG15，30 Ry，6³ k，与 `tests/fixtures/ivs_dfpt_si` 同一套）：空穴 Γ 谷内 D = 7.89 eV/Å @ 63.9 meV
    （冻声子 7.77，差 1.5%）；电子 Δ 谷谷内 D = 0（对称性禁戒，与冻声子一致）；谷间三条通道与之前逐位相同。
  - 局限：q=Γ 的 D 是 q→0 的值，没有在热载流子 q 范围内平均（散射率偏上限）；2D 没有库仑截断（谷内光学模与面外偶极相关的
    部分有 σh 时对称性禁戒，影响有限）。
  测试：`tests/test_0059_ivs_dfpt.py` 的 `test_r3_*`（8 个）。

## [0058] 形变势口径、多谷 TMD 电子的局限与文献参数入口（2026-10-18 用户定：保持现状，不上 EPW）

**决定**：2D 继续用真空参考的形变势（`DEFORM_REF=auto` → vacuum）；ODP/IVS 的常数形变势通道（`DP_OPTICAL_E/H`）
定位为**文献参数的接入入口**，流程本身不再扩展谷间散射的自动计算。

- **真空参考实现正确**：MoS₂ 的 S7.1 真空形变势与 Wiktor & Pasquarello PRB 94, 245411 (2016) 的绝对形变势一致
  （他们是单轴应力，按 ν ≈ 0.25 除以 1−ν 换成单轴应变）：K 导带 −8.1 / 本流程 −8.54，K 价带 −2.48 / −2.48，
  Γ 价带 +4.17 / +4.23，Q 导带都 < 1 eV。芯能级 h5 = 真空 h5 + 2.09 eV（带符号常数，S7.1 自检②）。
- **芯能级不当作修正**：芯能级下 MoS₂ 电子 K 6.44 / Q 2.19 eV 恰与 Filipovic 2026 的 DFPT（Ξ_LA 6.5；Q 谷 Ξ_L 3.0 / Ξ_T 2.4）接近，
  但同一口径下空穴 K 只有 0.38 eV（Jin 2014 DFPT 2.5）。没有一个统一参考能同时对上两种载流子，属巧合，不改默认。
- **本流程里的 DFPT 不能校准电声耦合**：VASP `IBRION=7/8` 只给力常数、ε、Born 电荷（S5、S7.3），不给电声矩阵元。
  冻声子劈裂法（S7.2 ODP、S7.4 IVS）只适用于 Γ 点光学模与**简并**谷间（K↔K′、Q↔Q′）；K↔Q 不简并（MoS₂ 差 129 meV），
  一阶劈裂为零，劈裂法提不出 D，需要直接算矩阵元（DFPT/EPW）。
- **多谷 TMD 电子的已知局限（S8.2 多谷判定"拒绝"即此情形）**：真空参考下 Q 谷形变势几乎抵消（MoS₂ 0.36、MoSe₂ 0.27、
  WS₂ 0.16、WSe₂ 0.11 eV，K 谷 7.4–9.4 eV），Q 谷电子的 ADP 接近 0；而 Q 谷的主要散射是谷间（按 Filipovic 的耦合估算，
  300 K 时 Q→Q′ + Q→K ≈ 1.2×10¹⁴ s⁻¹，约为 Q 谷内声学的 10 倍），出厂不算。所以这类载流子的声子限制迁移率**偏高**：
  MoS₂ 电子在冻声子 ODP、无 IVS 下约 377 cm²/Vs（300 K），DFPT 基准 187–212（Filipovic 2026）。
  单机制"只 ADP"在高温上升（MoS₂ 600 K 734）是同一原因，看 overall 不看它。空穴（S8.2 判 ok）不受此影响，但没开 SOC
  （MoS₂ 价带 K 劈裂 148 meV），也没有含 SOC 的 DFPT 基准可对照。
- **文献参数入口**（写法见 `step8.4_amset2d/README.md`「ODP 与谷间散射」）：
  - `DP_OPTICAL_E/H` 的 `X>X'`、`X>Y` 通道接谷间散射（IVS），`X>X` 接谷内 ODP；`ODP_FROZEN=off` 时用手填的谷内值。
  - step.conf 显式写了 `SCATTERING` 时不自动补 ODP/IVS（[0034]），要生效须把 `IVS` 写进去。
  - **参数选源**：MoS₂ 电子谷间推荐 Filipovic, arXiv:2607.18895 (2026) Table S4/S8（DFPT、2D 库仑截断，K–Q 21 条通道）。
    Jin 2014 Table IV 的 K→Q 近邻声学 9.3 eV/Å 约为 Filipovic（≈1.3）的 7 倍，会把电子 μ 压到 ~130，不要直接用。
    按 dp_optical 的星平均态密度换算（D × √(成员占比)：K→Q 近/远各 3/6，Q→Q′ 2/5、2/5、1/5，Q→K 1/2；同 ħω 两条按平方和合并），
    非 SOC 口径：
    ```
    DP_OPTICAL_E = K>K':1.66@40.6;
     K>Q:0.81@22.2; K>Q:1.02@42.7; K>Q:0.74@49.3; K>Q:0.40@45.5; K>Q:0.16@16.5;
     K>Q:3.50@49.1; K>Q:3.28@28.5; K>Q:2.02@43.9; K>Q:1.10@19.0;
     Q>Q':2.45@23.0; Q>Q':1.31@49.2; Q>Q':0.89@45.4; Q>Q':0.49@17.4; Q>Q':0.47@42.4;
     Q>Q':1.99@42.8; Q>Q':1.88@49.1; Q>Q':2.59@28.5; Q>Q':1.31@22.5; Q>Q':1.09@47.4; Q>Q':0.64@40.6;
     Q>K:0.81@22.2; Q>K:1.02@42.7; Q>K:0.73@49.3; Q>K:0.41@45.5; Q>K:0.16@16.5;
     Q>K:3.50@49.1; Q>K:3.28@28.5; Q>K:2.02@43.9; Q>K:1.10@19.0
    ```
    一阶 K↔K′ 通道（∝|q−K|）本模型表示不了，缺它 μ 约偏高 8%（Filipovic Table 5）。用了文献参数的结果须在报告里注明来源。
- **MoS₂ 文献基准（300 K，声子限制）**：电子 Filipovic 2026 EPW 187–188 / EMC 212；Cheng & Liu JACS 140, 17895 (2018)
  EPW 395（只取带边态、无 2D 截断）；Li 2013 / Kaasbjerg 2012 只算 K 谷 320 / 410（Kaasbjerg 的 LA 耦合约为 DFPT 的一半）；
  Jin 2014 含 K→Q 130。空穴 Cheng & Liu 219、Jin 270（均非 SOC）。各家 DFPT 之间本身差约 2 倍。
- **2026-10-18 MoS₂ 的 S8.4 是对照配置**：`ODP_FROZEN=off`、手填 K>K 4.1 / Q>Q 7.1（Kaasbjerg/Jin）、显式
  `SCATTERING=ADP,IMP,POP,ODP`（写在 `DP_OPTICAL_E` 里的 Jin 谷间通道因此未生效）。300 K 电子 overall 171 与基准接近，
  部分来自误差抵消（K 谷 ADP 偏强、Q 谷 ADP≈0、Q 谷手填 ODP 偏大），不代表出厂流程的精度。

## [0057] transport_plots 改版：目录改名、μ 按机制拆、温度蓝→红、μ − E_VBM 扫描（2026-10-18 用户定）

| 新目录 | 原目录 | 内容 |
|---|---|---|
| `transport_models_<段>/` | `phonon_limited_<段>/` | 不变：基础 4 条（phonon + impurity / phonon-limited / ADP only / BT2×DPT-τ）× μ、S、σ、PF、κe |
| `mobility_mechanisms_<段>/` | `transport_<段>/` | **只画 μ**，每温度一个面板（3×3）：AMSET 每个散射机制单独一条（transport.json 的单机制 μ，虚线）+ phonon + impurity 总值（黑实线）+ DPT 带边横线；不再画 phonon-limited 线、不再重复 S/σ/PF/κe 四列 |
| `temperature_series/` | `datasets/` | 每个数据集一张、每温度一条线；颜色低温蓝 → 高温红；按浓度两段 `<M>_<数据集>_{carrier_tuning,low_doping}.png/.csv`，同一个文件夹 |
| `mu_scan/` | `bt2_muscan/` + `amset_scan/` | 600 K，横轴 μ − E_VBM：BT2×DPT-τ 连续曲线 + AMSET 点 |

- **为什么**：`transport_<段>` 与 `phonon_limited_<段>` 的 S/σ/PF/κe 四列几乎一样（transport 只多了部分机制子集，
  MoS2 的 ADP+POP+ODP 子集就是 phonon-limited 本身）；μ 列又没把每个机制拆开。
- **mobility_mechanisms**：单机制 μ ≥ 1e10（该机制前因子为 0）不要；超过「总值 / DPT 最大值 × 100」的不画
  （低浓度 IMP 可到 1e7+，会把总值压成一条线），图注写明、值全在表里（`mu_<机制>_only_cm2/Vs`、`mu_total_cm2/Vs`、
  `mu_DPT_band_edge_cm2/Vs`）。黑线是 AMSET 全部机制一起算（每个态的散射率先相加再热平均），比单机制曲线的
  Matthiessen 倒数和低 —— 两者不一样是 RTA 下的正常现象，不是错。
- **mu_scan 的换算**（两边能量零点不同，必须各减各的 VBM）：
  - BT2：x = μ_rel − VBM_rel，μ_rel 是 S8.1 网格（零点 = S8.1 的 `data.fermi`）。VBM_rel 优先读 boltztrap_crta.json 新增的
    `band_edges_rel_efermi_eV`（S8.1 [0057] 起写，剪刀后的能带）；旧结果做过剪刀时 S8.1 把 E_F 放在带隙正中
    （`data.fermi = VBM + Eg/2`），所以 VBM_rel = −Eg/2；两者都没有（没做剪刀，零点 = VASP E_F）就**不画 BT2**，不猜。
  - AMSET：x = fermi_levels − E_v，fermi_levels 是 transport.json 里 AMSET 对每个掺杂档自己解的 E_F（eV）；
    E_v 读 undoped_ni.json（amset_fermi_fix 用 AMSET 自己的 DOS 求的带边，剪刀后，与 fermi_levels 同零点）；
    缺时退回运行目录 vasprun.xml 的 DFT 带边 − Δ/2（AMSET 剪刀对半分）。顺带自检：transport.json 本征点 E_F
    与 undoped_ni.json 记录的 E_F 差 > 1 meV 打 WARN。
  - BT2 的 σ、κe = BoltzTraP2 σ/τ、κe/τ × S8.3 同一个 DPT τ(T)（n < 0 电子 τ、n > 0 空穴 τ），面内平均、层厚归一，
    与本征表 BT2 行同一口径（不再用 paper_mu_scan 的 x/y 分方向、WF κe）；双极区不给 μ。
- **mechanism_resolved_mobility.png**：去掉核对圆圈（相邻掺杂档的同一总值，CSV 里 role = check 仍保留，
  相差 > 5% 照旧打 WARN）；标题拆两行，写明黑线 / 点线的含义。
- 只改出图：已有结果用 `tools/make_transport_plots.py <材料>/ke-dft-cpu` 重出即可，不用重跑 S8.3；mu_scan 的 BT2 曲线
  对没做剪刀的旧 S8.1 结果要重跑 S8.1 才有。
- **r1（2026-10-18 用户本地重出后报）**：MoS2 的 mu_scan 一个 AMSET 点都没有 —— 本征表按 `%.6g` 写浓度
  （−2.86705e20），transport.json 是全精度（−2.867034…e20），配对容差写成了 1e-9。改为 `_same_n`：同号且相对差
  ≤ 1e-4（网格间距至少 ×1.5，不会配错）；mobility_mechanisms 表的 n_2D / n_lit / in_interpret_range 列是同一个
  问题（配不上就空着），一并修。另：BT2 在带隙深处 σ/τ 下溢成 0，S 是 0/0 占位值（实测一段 S ≡ 0），这些点不画；
  分机制图对数轴只留主网格。
- **r2（2026-10-18 用户「温度图更好看、类似阴影」）**：temperature_series 每个温度一条细线，相邻两温度之间按色带半透明填充，
  整组读成一条蓝（低温）→ 红（高温）的渐变带；9 项图例换成横向温度色标。
- **r3（2026-10-18 用户看 MoS2 真实 mu_scan「不太对」）**：BT2 的 μ = σ_总/(e|n_净|) 在带隙中央两侧竖直掉下去
  （n_净 是 N + nelect 的相减残差、少子不可忽略），κe 在同一处掉到 1e-21（κ0 − L1²/L0 相减抵消）。现在 BT2 μ 只在离
  带隙中央 > 5 kT 处给，κe 在 Lorenz < 0.1·L0 时不画；图注加一行：掺杂低于 n_i(T) 的 AMSET 点 E_F 都钉在带隙中央，
  点会叠在一起（S 正负都有 = 双极），这是物理，不是配错。
- **r4（2026-10-18 用户问）**：① 分机制图低浓度 IMP「只有两个点连成线」—— 超过上限的点不画后线在最后一个点就断了；
  现在从该点按 log μ 线性延到上限处、末端画 ▲（= 继续往上，μ_IMP ~ 1/N_ii）。② mu_scan 的 AMSET 只在掺杂档上算
  （每档一个 E_F，没有连续 μ 网格），仍画点，但同一载流子类型的点按 μ 连细线；带隙中央 ±5 kT（E_F 钉住、叠在一起）的不连。

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

- **r3（2026-10-18 用户要求）图注**：四张分段图底部加图注 —— 每条曲线 = 哪种计算 + 用了哪些散射机制
  （phonon + impurity = 本轮 transport.json 全部机制；phonon-limited = 去 IMP 的子集；ADP only；BT2×DPT-τ =
  BoltzTraP2 常弛豫时间、τ 取 DPT 的 Bardeen–Shockley 值，只含声学形变势），并给出现过的机制缩写全称
  （ADP/POP/ODP/IVS/IMP/PIE）。刚性能带近似不成立的浓度（本征表 `in_interpret_range = NO`，每原子载流子数 >
  阈值，与 S8.4 同判据）灰底标出；以前 9×5 图不标，最外侧（最高浓度）那几档看起来像"数据有问题"，其实是这条判据外。
  载流子调控段上限由 ∞ 改为轴端 1e21（轴外的点不再画出边框）。独立脚本的旧图例 `DPT route` / `DPT band-edge e- h+`
  改为 `BT2×DPT-τ` / `DPT μ, electrons|holes (band edge)`，图注说明横线 = Bardeen–Shockley 带边迁移率
  μ = eτ/m*，每个载流子、每个温度一个数，与浓度无关。
- **化学势 μ 与浓度 n**：BoltzTraP2 扫化学势，每个 μ 由插值能带的 DOS 算净载流子 n(μ, T) =
  (N(μ, T) + nelect)/V_cell（N 为 −∫DOS·f dE，负 = 电子）；固定 T 下一一对应、单调。S8.3 的 BT2×DPT-τ 曲线把
  BT2 的 (n, 输运量) 在同号支上对 log|n| 插值到 AMSET 的掺杂档（不外推）。

- **r4/r5（2026-10-18 用户要求）统一 + DPT 横线**：`tools/make_transport_plots.py` 不再有自己的算法，直接调用
  `gen_step13_output.write_transport_plots(cwd, out, is_2d, plots=<输出目录>, M=<材料名>)`，产物与 S8.3 逐字相同
  （测试 test_0056 逐个 CSV 比对）。以前它自己算的 BT2×DPT-τ 与 S8.3 不是一回事：

  | | 以前的独立脚本（已删） | S8.3（现在两边都用这套） |
  |---|---|---|
  | BT2 横轴 | BT2 每个化学势点自己的 n(μ)，一直扫到 1e21 以上 | BT2 的 (n, 量) 按 log\|n\| 插值到 AMSET 掺杂档，不外推 |
  | σ/τ、S、AMSET | 只取 xx | 面内平均 (xx+yy)/2（3D 迹/3） |
  | τ | DPT x 方向 τ(300 K) × 300/T | 每个温度的 DPT τ = μ·m*/e |
  | BT2 κe | Wiedemann–Franz L0·σ·T（L0 = 2.44e-8） | BT2 自己的 κe/τ × τ |

  AMSET 方向相反：给定掺杂浓度求费米能级（`amset_fermi_fix`），没有 μ → n 换算。
  transport 图的 μ 列加 DPT 带边迁移率横线（本征表的 `DPT` / `DPT_mv` 行）：电子只画 n 型（左）半边、空穴只画
  p 型（右）半边；μ 轴顶部留 22% 空白，数值标注写在空白里（左上电子、右上空穴），不压曲线。
  目录名为 `transport_plots` 才清空旧结果（独立脚本指定别的输出目录时不删）。
- **S8.1 低温锯齿**：BT2 自己的 μ 扫描在 100–200 K、高浓度出锯齿，可能是 DOS 直方图格宽（`DOS_NPTS = 4000`
  覆盖全部能带）≳ kT。统一后 BT2 只在 AMSET 掺杂档取插值、不画 1e21 以外，影响已很小；要消掉需改脚本常量并
  重跑 S8.1 → S8.3（全部材料），暂不做（用户 2026-10-18：没必要）。

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
