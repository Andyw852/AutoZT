# step8.4_amset2d -- AMSET 二维散射核（插件式）

> 只在 2D 项目里打开；三维项目走原版 `step8_amset`，完全不受影响。
> 开关：项目 `project_setting` 里写 `amset2d: true`（技能里该组 `default: false`、`auto_dim: 2d`：
> 2D 材料 init/register 时自动写这一行；老材料手写）。zt-dft-cpu 的 S20 对 2D 只认本步的 transport.json。
> 补丁 0016 起 2D 新材料默认不跑原版 S8（`amset3d: false`）；层厚写在本材料共用 step.conf 的 `LAYER_THICKNESS`。
> `2d_correction.json` 的 `za_coupling` 记录结构有没有水平镜面 σh（没有时本步不含的 ZA 一阶耦合可能不为零）。

## 这个步骤解决什么

原版 AMSET 的散射核是三维的。把它直接用在真空 slab 上会有三个问题：

1. **弹性张量顺序**（三维也中招，已在 gen_step10 / gen_step14 一起修）
   VASP `TOTAL ELASTIC MODULI` 的行列序是 `XX YY ZZ XY YZ ZX`，而 AMSET 经 pymatgen
   按标准 Voigt `XX YY ZZ YZ XZ XY` 解析。不重排就把面内剪切 `C66(XY)` 塞进 `YZ` 槽位。
   实测 CrS2 单层：`C66 = 22.29 GPa` 被换成 `C_zxzx = -0.10 GPa`（面内 TA 支刚度近零）。
2. **ADP 不该乘 c/t**。AMSET 的散射积分在超胞里做，二维核的正确输入是**原始 slab
   弹性常数**（`= C_2D/c`）。乘了 c/t 会让 ADP 散射率偏小 `t/c` 倍（CrS2 约 1/3，
   对应迁移率偏大 3.06 倍）。
3. **散射核的 q 依赖**。二维 Fröhlich / 压电 / 带电杂质 / 声学形变势在 `q->0` 与
   大 `q` 的行为都与三维不同；真空方向的 `q_z` 在二维里必须积分掉（映射系数见下）。

## 映射关系（一句话）

AMSET 的三维积分里，真空方向 `q_z` 的积分贡献因子 `2*pi/c`，因此

```
G_s = c * G_2D
```

`G_s` 是插件塞进 AMSET `prefactor x factor` 的量。对 ADP 的等价说法：直接用原始
slab 弹性常数（`C_2D/c`）即可，不要再乘 c/t。推导与各机制的式子见
`../METHODOLOGY.md` 第 4 节；实测验证见 `VERIFICATION.md`。

## 目录里的文件

| 文件 | 作用 |
|---|---|
| `gen_step14_amset2d.py` | 登录节点 gen：写 settings.yaml + 2d_correction.json + submit.sh，复制插件 |
| `amset2d_plugin.py` | 运行期插件：把 ADP/POP/IMP/PIE 四类换成二维核（同名注册表替换，不改 AMSET 安装）；[0033] 可选 ODP/IVS 两个新机制与 ADP 谷内掩码 |
| `dp_optical.py` | [0033] 光学形变势（ODP）与谷间散射（IVS）的常数形变势模型：解析散射率、谷归属、step.conf 解析（纯 numpy） |
| `amset2d_pop_freq.py` | 取 Gamma 点**面内极性模**有效频率（AMSET 自带口径会把面外 ZO 模算进去） |
| `amset_desym_fix.py` | 运行期插件（S8 与 S8.4 共用）：修 AMSET 去对称化的相位因子（V115）；环境变量 `AZ_DESYM_FIX` 未设时不打补丁 |
| `amset_fermi_fix.py` | 运行期插件（S8 与 S8.4 共用）：AMSET 求不到费米能级时改用二分法、宽松解偏差 > 1% 报错、NaN 报错（V132）；`AZ_FERMI_FIX=0` 关 |
| `test_desym_fix.py` | 相位补丁的模型检验（紧束缚真值，12 构型；需要 amset） |
| `test_fermi_fix.py` | 费米能级插件的检验（WS2 同款 DOS 复现与修复、原式能解时逐位不变、接线） |
| `test_gen_desym_wiring.py` | S8/S8.4 gen 的逐操作判据、命令、覆盖通道接线自检 |
| `test_kernels.py` | 四类核的自检（旋转不变、极限行为、护栏），`python test_kernels.py` 应全 PASS |
| `test_dp_optical_e2e.py` | [0033] ODP/IVS 与 ADP 谷内掩码的真 AMSET 端到端自检（合成双谷模型，5–15 分钟；需要 amset） |
| `tau_energy.py` | [0035] 作业链尾部：从 mesh 出逐机制 τ(E) 与迁移率等效 τ，写 `tau_energy.json` |
| `test_tau_energy_e2e.py` | [0035] τ 的真 AMSET 端到端自检（常数率精确性、复现、τ_σ 交叉核对、与 Bardeen–Shockley 对照；需要 amset） |
| `example_2d_correction.json` | 自检用的最小 2d_correction.json 示例 |
| `VERIFICATION.md` | 在真实 2D 数据上做的端到端验证记录（T1/T2 + 三材料对比） |

## 运行前自检

```bash
cd skill/ke-dft-cpu/step8.4_amset2d
python test_kernels.py          # 需要 amset + pymatgen 环境（如 conda amset_clean）
```

集群上确认 `python -c "import amset; print(amset.__version__)"`：插件在
**0.4.19 与 0.5.1** 上都核对过（两者的 `scattering/elastic.py`、`inelastic.py`、
`common.py` 逐字节相同）。其它版本会打印 WARN 但继续跑，请自行核对散射接口。

跑完后在 `amset.log` 里应看到：

```
[amset2d] amset 0.4.19 | c=20.000 A | r_inf=[[...]] | mechanisms=ADP,POP,IMP
[amset2d] POP: r_inf=46.500 A | delta_r=0.5300 A | delta_r/r_inf=1.14% (G_2D = 2*pi*c*Delta_r/[(1+r_inf*q)(1+r_0*q)])
[amset2d] worker: ADP=AcousticDeformation2D POP=PolarOptical2D IMP=IonizedImpurity2D PIE=Piezoelectric2D
```

没有这两行说明插件没生效（子进程尤其要看第二行 -- ADP 是以 reference 传给 spawn
子进程的，插件没被导入就会**静默**退回原版）。

## 去除 IMP 的无杂质结果（`WRITE_MESH`）

S8.4 默认 `SCATTERING = ADP,IMP,POP`。因此默认 `transport.json` 含电离杂质散射，对应给定载流子浓度下包含 IMP 的输运结果。
`intrinsic_transport.json` 是把 IMP 散射率移除后、对其余已启用机制重新积分的结果：默认机制下为 ADP+POP；若启用 PIE、ODP 或 IVS，这些机制仍会保留。
它表示当前模型和给定温度、载流子浓度下的无电离杂质结果。若 ODP、IVS、PIE 或其他重要耦合未纳入，不能把它解释成真实晶体迁移率的绝对上限。
**[0032] 起 `WRITE_MESH` 与 `INTRINSIC_ADP_ONLY` 出厂即开**（S8.4 只跑 2D）：每个材料跑完自动有
`intrinsic_transport.json`（移除 IMP，保留其余已启用机制）和 `intrinsic_ADP.json`（只保留 ADP；若 IVS 已启用，ADP 按默认掩码限制为谷内）。启用 ADP、POP 时还会生成 `intrinsic_ADP_POP.json`；同时启用 ODP 时再生成 `intrinsic_ADP_POP_ODP.json`。S8.3 会读取可用且与本轮输运结果匹配的文件，绘制分级散射对照图表
（`step8.3_output/intrinsic/`，见 VERIFICATION V195）。mesh_*.h5 不回拉本地（skill.yaml 的 `fetch_exclude`）。
要关掉，在 `project_setting/templates/step8.4_amset2d/step.conf` 写：

```ini
[params]
WRITE_MESH = false
```

gen（`patch_write_mesh` + `patch_intrinsic_auto`）会：

1. 在 `settings.yaml` 写 `write_mesh: true` → AMSET 落 `mesh_*.h5`（每机制每 k 点散射率 + 能带/速度/DOS 元数据）；
2. 作业链尾部**自动**接 `python postprocess_intrinsic.py --check-reproduce` → 产出 `intrinsic_transport.json`。

`postprocess_intrinsic.py` 做三件事：

- 用 **mesh 自带**的 `ir_to_full_kpoint_mapping` 把不可约点散射率展开到全网格；`read_mesh_h5` 自己解析
  AMSET `write_mesh` 的**自旋 key bug**（down 被写成 `<name>_up_down`），不依赖 `amset.io.load_mesh`；
- `--check-reproduce`：**先用全部机制重积分**，与原始输运结果核对；默认最大绝对误差容差为 `1e-4`。检查迁移率、电导率和 Seebeck 系数，不是逐位比较，也不覆盖 κe。**不通过就整步失败、
  不写 `intrinsic_transport.json`**，防止 mesh/settings/vasprun 不同源时静默出错；
- 置零（默认）IMP 后重积分，写 `intrinsic_transport.json`（schema `autozt.intrinsic_transport/1`）；按实际启用机制理解该文件，默认机制以外的 PIE、ODP、IVS 不会被自动删除。

**注意**：`WRITE_MESH` 只负责产出 mesh；**正确性取决于这个后处理被套用**——直接换 AMSET 环境不会自动对。

`POSTPROCESS_ONLY = true` 可跳过 AMSET、复用已有 mesh 和输运结果来补做后处理。仅当 `mesh_*.h5`、`transport_*.json`、`transport.json`、`settings.yaml`、vasprun、散射机制、掺杂和温度网格都来自同一次计算且配置一致时使用；文件存在本身不能证明来源一致。该模式会清理旧的派生 JSON，但保留 mesh 和输运源文件。

**overall 去掉某些机制（[0034]）**：`OVERALL_WITHOUT`（出厂 `ODP`）列的机制本次开着时，作业链再跑
`postprocess_intrinsic.py --drop <列表> --out overall_no<列表>.json`（如 `overall_noODP.json`），用来看该机制对总结果的影响。

**只含 ADP 的一套（[0018]）**：`INTRINSIC_ADP_ONLY`（[0032] 起 2D 出厂开，隐含 `WRITE_MESH = true`）。作业链在上面
的后处理之后用**同一份 mesh** 再跑 `postprocess_intrinsic.py --drop IMP,POP,PIE,ODP,IVS --out intrinsic_ADP.json`
（[0033] 起把 ODP/IVS 也剔除；没开的机制自动忽略）。各机制散射率彼此独立，等价于单独跑一次只有 ADP 的 AMSET
（开了 IVS 时这个 ADP 是谷内 ADP，见下面「ODP 与谷间散射」），不用第二次 AMSET。
**递增机制对照（[0037]）**：当 `SCATTERING` 同时包含 ADP、POP 时，另用同一 mesh 生成
`intrinsic_ADP_POP.json`（移除 IMP、PIE、ODP、IVS）；若 ODP 也启用，则生成
`intrinsic_ADP_POP_ODP.json`（移除 IMP、PIE、IVS）。因此后两份分别保留 ADP+POP 和 ADP+POP+ODP，便于分辨加入 POP、ODP 后的变化；它们不含 IMP、PIE 或 IVS。启用条件不满足时不生成相应文件。
`intrinsic_transport.json` 则保留所有实际启用机制、只移除 IMP；`transport.json` 保留全部实际启用机制。默认机制下分别是 ADP+POP 和 ADP+IMP+POP。`OVERALL_WITHOUT`（[0034]）可额外生成 overall 去除指定机制的结果。
本征 JSON 含 μ、S、σ、κe（κe 自 [0018] 起写出）；PF 用 S²σ 自算。都在本步目录，不另建步骤目录。
0.4.19 与 0.5.1 在真实 2 自旋 / 不可约网格数据上的端到端复现均 PASS（见 `VERIFICATION.md` V89/V90/V91）。

## ODP 与谷间散射（[0033]；谷内 ODP 自 [0041] 起默认由冻声子参数开启，谷间 IVS 仍需手填）

AMSET 只有 ADP / POP / IMP / PIE，**没有非极性光学声子（ODP）和谷间散射**。MoS₂ 这类 TMD 的同极 A1′ 模与 K→K′、K→Q
谷间声子很重要（Kaasbjerg 2012、Jin 2014），缺了它们总迁移率系统偏高。[0033] 用文献的**常数形变势模型**补上：

    1/τ_c(x) = π D_c² / (M ω_c) · [ (N + f(x + ħω_c))·G_c(x + ħω_c) + (N + 1 − f(x − ħω_c))·G_c(x − ħω_c) ]

- `D_c`：通道 c 的形变势（eV/Å）；`ħω_c`：该通道声子能量；`N`：玻色占据；`M`：原胞质量。
- `G_c`：目标谷的态密度（每原胞、每自旋），**直接用 AMSET 插值网格上的真实能带**统计（不假设抛物带），按星平均。
- `f`：末态的载流子占据（[0042] 起逐掺杂计算泡利阻塞，与 AMSET 自己的 POP 同一口径；等价于 Kaasbjerg Eq.30 的
  Σ P(1−f′)/(1−f)）。以前没有 f，相当于非简并极限（Kaasbjerg Eq.31），简并掺杂（热电区间）下 ODP 有 5–30% 量级的偏差。
- 发射/吸收阈值：目标谷在 x ± ħω 处没有态就是 0；但 `G` 做了 6 meV 高斯展宽，阈值附近 ±2σ 被抹圆
  （2D 只看 ODP 的 ⟨τ⟩ 约偏低 10%，总迁移率影响更小）。矩阵元与 q 无关 → 各向同性，动量弛豫率 = 散射率。
- 体积（含真空）在 D²·G/M 里抵消，与真空长度无关。
- 走 AMSET 的 basic 接口（每 k 点直接给率），不经四面体积分，几乎不增加耗时。

**参数从哪来**：谷内 ODP 用冻声子法自动算（[0035]，见下一节「冻声子 ODP 参数」，不用文献；有算出来的值时手填的谷内值
只作对照）；谷间 IVS 没有自动参数，只能手填文献值。

**手填的开法**：材料的 step.conf 写通道表（`起>止:D@ħω`，D 单位 eV/Å = 1e8 eV/cm，ħω 单位 meV），**数值从文献表里抄**
（Jin et al. PRB 90, 045422 (2014) 的 Table IV/V，arXiv:1406.4569；或 Kaasbjerg PRB 85, 115317 (2012)）：

```ini
[params]
DP_OPTICAL_E = K>K:<D>@<ħω>; K>K':<D>@<ħω>; K>Q:<D>@<ħω>; Q>Q:<D>@<ħω>; Q>Q':<D>@<ħω>; Q>K:<D>@<ħω>
DP_OPTICAL_H = K>K:<D>@<ħω>; K>K':<D>@<ħω>; K>Γ:<D>@<ħω>; Γ>K:<D>@<ħω>
# 可选：
# DP_OPTICAL_VALLEYS_E = K:0.333333,0.333333,0|-0.333333,-0.333333,0; Q:...   # 手填谷中心（S3 晶胞分数坐标）
# ADP_INTRAVALLEY = auto        # auto / on / off
```

- `X>X`：谷内（同一个谷）→ **ODP**；`X>X'`：同一星的其它谷（K→K′）、`X>Y`：到别的星（K→Q）→ **IVS**。
  文献按声子支分列的（如 K→Q 的 TA、LA、LO 各有 D 与 ħω），**每支写一条**（`K>Q:1.2@24; K>Q:2.0@48`），各支的率相加。
- **谷表**默认读 S8.2 的 `step8.2_dpt/dpt_result.json`（`multi_valley.valleys`，拒绝出数时谷表照样在）；没跑 S8.2 或不对，
  用 `DP_OPTICAL_VALLEYS_E/H` 手填。通道两端的谷必须在谷表里，否则 gen 直接报错。
- 写了通道表 gen 就自动把 ODP / IVS 加进 SCATTERING（日志 `[OK] ODP/IVS 开：…`），并写进 2d_correction.json 的 `dp_optical`。
  **例外**（[0034]）：step.conf 显式写了 `SCATTERING`（如严格对照写 `SCATTERING = ADP`）时照它来，不补；
  要开就把 ODP、IVS 写进去（`SCATTERING = ADP,IMP,POP,ODP,IVS`），没写进去的那类通道打 WARN、不生效。
- S8.2 的谷表里，同一组谷可能因为"按能带/自旋分星"出现两次（第二个标成 `K(up带27)` 这样）。本模型按谷中心归属 k 点、
  态密度统计所有能带，它们就是同一个谷：gen 自动并成一个；通道里写被并掉的标签直接报错（请用 `K`）。
  判"同一组谷"用 S8.2 的网格极值点 `k_grid`（Q 这类非高对称谷的拟合中心按带各自浮动，不能拿来比）。
- S8.2 重算后谷变了（网格极值点集合变了）、而 S8.4 的 ODP/IVS 用的是 S8.2 的谷表：S8.2 自动把 S8.4 的完成标记归档，
  S8.4 重新生成、重新排队。谷没变、或用的是手填谷表就不动。S8.4 的作业正在排队/运行时不归档，打 ★ 提示先 scancel 再 retry。
- **ADP 谷内掩码**（`ADP_INTRAVALLEY = auto`，出厂）：开了 IVS 的载流子，ADP 只留谷内（k 与 k+q 落在同一个谷）。
  AMSET 的 ADP 对整个布里渊区的 q 积分，大 q 的 K→K′、K→Q 用长波声速和带边形变势算，既不物理、又和 IVS 重复计数。
  只开 ODP（没有谷间通道）时 auto 不掩。掩码同样作用于「只 ADP」那一套（同一个 ADP 机制）。
- **输运结果与分级对照**（以下派生 JSON 需 `WRITE_MESH = true`，同一次 AMSET、同一份 mesh 重积分；可用结果由本次 `SCATTERING` 机制决定，S8.3 读取有效文件）：

  | 文件 | 含哪些机制 | 用途 |
  |---|---|---|
  | `intrinsic_ADP.json` | 只 ADP；还需 `INTRINSIC_ADP_ONLY` | 声学散射参照，可与文献 ADP/DPT 对照 |
  | `intrinsic_ADP_POP.json` | ADP+POP；需同时启用 ADP、POP | 分离 POP 对声学参考结果的增量 |
  | `intrinsic_ADP_POP_ODP.json` | ADP+POP+ODP；需同时启用 ADP、POP、ODP | 分离 ODP 对 ADP+POP 结果的增量；不含 IVS |
  | `transport.json` | 本次启用的全部机制 | overall；默认 ADP+IMP+POP |
  | `intrinsic_transport.json` | 本次启用的全部机制减 IMP | 无电离杂质结果；默认 ADP+POP，不代表未纳入机制后的真实绝对上限 |
  | `overall_noODP.json` | 本次启用的全部机制减 ODP | 看 ODP 对 overall 的影响（[0034]，默认 `OVERALL_WITHOUT = ODP`；仅在 ODP 启用时生成） |

  "只 ADP"就是把其余机制的散射率全部置零后重积分；ODP/IVS 是新加的机制，所以也在置零之列，
  否则那一套就成了 ADP+ODP+IVS。
- `OVERALL_WITHOUT`（[0034]）：出厂 `ODP`。写 `ODP,IVS` 得到 `overall_noODP_IVS.json`（去掉 ODP 与 IVS）；
  写 `off` 不出。列出的机制本次都没开时这一套就是 overall，自动不出。
- **"只 ADP"与文献对得上吗**：
  - 没写 `DP_OPTICAL_*`：ADP 与以前逐字节相同（不掩码）。
  - 开了 IVS 的载流子，ADP 只留谷内，正是形变势理论（Bardeen–Shockley / Takagi、Kaasbjerg、Jin 2014 的声学项）
    和 S8.2 DPT 的定义 —— 谷间散射另由 IVS（区边界声子能量 ħω）给出。所以"只 ADP"反而更接近文献口径。
  - 量级参考：V25.5 在 MoS₂ 真实重叠上做过 K↔K′ 谷间置零的诊断（与本掩码同一思路，当时谷表只有 K/K′），
    ADP-only 电子 250 → 311（+24%）、空穴 1070 → 1162（+9%）。谷表里有 Q 时，掩码还会去掉 K↔Q 的声学散射。
  - 想对照开 ODP/IVS 之前的结果，直接用这个材料上一轮（没写 `DP_OPTICAL_*`）的 S8.4 结果。

**跑完看日志**（`amset.log` 前面几行）：

- `[amset2d] ODP/IVS（常数形变势…）：electron[K×2,Q×6] K>K:…` —— 实际生效的通道；
- `[amset2d] 谷表核对（electron）：K 谷底离中心 0.000|b1|、高于带边 0.000 eV; Q …` —— 谷底要落在谷中心（< 0.05|b1|），
  Q 比带边高多少应与 S8.2 的 `dE_eV` 一致；带 `[WARN]` 说明谷表与 AMSET 的晶胞/网格不同源，先查再用；
- `[amset2d] ODP/IVS：FD 窗口外保留解析率` —— ODP/IVS 在费米窗口外不用 AMSET 的对数平均填充（见 V196）；
- transport.json 里多出 ODP、IVS 两个单机制迁移率；某载流子没写通道时它的单机制值会极大（只有 2e6 s⁻¹ 的下限率），这是正常的。

**局限**：常数形变势（q 无关、各向同性）是文献的简化模型，参数质量决定结果；不含一阶（∝q）谷间耦合；
**模型局限（报数时注明）**：
- 常数 D 取的是 Γ 点（q = 0）的耦合。冻声子给的是 Γ 点裸值，而真实 |M_q| 离开 Γ 点会下降，所以速率是**上限**、
  迁移率偏低：Kaasbjerg MoS₂ 同极模 Γ 点 ≈4.8 eV/Å、拟合散射率的等效常数 4.1 eV/Å，速率差 ≈1.37 倍。
- 谷间散射（IVS）只认手填文献参数，出厂不算。多谷材料（MoS₂ 的 K/K′、Si 的 X×6）缺谷间光学散射，迁移率偏高
  （MoS₂ 仅 K→K′ 的 LO 通道就约占零阶 ODP 率的 30%）；gen 对多谷载流子会打 WARN。
- 极性模默认不算（长程部分在 POP）；它们对称性允许的短程零阶耦合（如闪锌矿空穴的 d₀）因此缺失。
- 简并带（d 条）用 D² = Σλ²/d² 配全部 d 条带的态密度，各带质量相同时精确；轻重带质量不同、或劈裂刚超过 5 meV
  （按单带算）时有偏差。

谷归属按「离谷中心最近」划分。SOC 能带（旋量带）下目标谷态密度含两种自旋，自旋守恒的通道（MoS₂ 价带 K↔K′）会偏大到 2 倍，日志有 WARN。通道里的谷标签要与 S8.2 谷表一致（Γ / K / M / Q / Γ–M；G、G-M 也认）。要发表的绝对值仍建议与 EPW/Perturbo 对照。

## 冻声子 ODP 参数（[0035]；[0041] 起出厂开）

谷内 ODP 的 D 与 ħω 不再需要文献：两个新步骤按材料自动算，S8.4 自动用上。

**开关**：[0041] 起 `odp_frozen` 组出厂开；某材料不要就在它的 `tf_*.yaml` 写一行 `odp_frozen: false`（或 S8/S8.4 的
step.conf 写 `ODP_FROZEN = off` 改用手填/文献值）。[0042] 极性判据：2D 面内极性只看面内模有效电荷；面外极性只对有
水平镜面 σh 的材料剔除（Janus/翘曲材料的同极类模保留）。多出两步：
（[0036] 3D 材料也能开：参数给 S8（原版三维 AMSET + `step8_amset/amset_dpo3d.py`）；不对齐、谷表由 S8.21 在 S7 网格上
自动找或用手填 `DP_OPTICAL_VALLEYS_E/H`，不需要 S8.2。规则同下：谷内 ODP 算出来的优先，谷间 IVS 只认手填。）

| 步骤 | 做什么 |
|---|---|
| S7.2_odp（`step7c_odp_frozen`，VASP 扇出 `odp-*`） | 读 S5 DFPT 的 Γ 点声子；每个**非极性光学模**做 ±u 两个位移单点，再加一个不位移的参考 `odp-ref`。INCAR/KPOINTS/POTCAR 照抄 `step7_deform/undeformed`（同一套 k 网格、泛函、ISYM=0、LVHAR） |
| S8.21_odp（`step8.21_odp_read`，登录节点） | 按 S8.2 的谷表，在每个谷的 k 点（`k_S7`）取带边能量，算 D，写 `odp_params.json` + `odp_summary.txt` |

S8.4 的 `ODP_FROZEN = auto`（出厂）读 `odp_params.json`，给每个谷加谷内通道 `X>X:D@ħω`（每组光学模一条），
把 ODP 加进 SCATTERING。满足机制条件时会生成 `intrinsic_ADP_POP_ODP.json`；默认 `OVERALL_WITHOUT = ODP` 也会生成 `overall_noODP.json`，可分别比较去杂质的 ODP 增量与包含杂质时 ODP 对 overall 的影响。

**算法**（推导在 `odp_frozen.py` 文件头）：
- 位移 Δ_κ = u·√(M/m_κ)·e_κ（e 为质量加权归一本征矢，M 原胞质量），u 取到最大原子位移 = `ODP_DISP_A`（出厂 0.01 Å）。
  这样 D = dE/du 正好对上 dp_optical 的 |g|² = ħD²/(2NMω)，不用再换算。
- 非简并带：中心差分 D = [E(+u) − E(−u)]/2u，二次项相消；对称性禁戒的模自动给 D ≈ 0。
- 简并带（u = 0 时相差 ≤ 5 meV 的 d 条带）：D² = Σλ²/d²，Σλ² 用迹恒等式算（与带的排序无关；与 dp_optical
  态密度统计全部能带的口径一致）。
- 组外邻带近到位移下可能交叉时做核对：排序配对与迹恒等式一致才用；不一致（交叉或强混合，这时单带的 D 随位移变）
  就报错，提示把 `ODP_DISP_A` 调小后 retry S7.2，不硬给一个数。
- 同一组简并声子的各伙伴 D² 相加；星内各成员谷取平均，离散度 > 10% 告警。
- 能量用真空能级对齐（LOCPOT 平面平均，与 S7.1/S8.2 的 E1 同一口径）；`odp_summary.txt` 同时列不对齐的 D。
  真空区若是一段直线（构型有面外偶极时周期性 slab 的补偿电场），取窗口中点（= 两侧真空能级的平均）；平坦度按扣掉直线后的残差判。
- 极性（红外活性）模默认不算，分两类（gen 的模表打印每组的模有效电荷"总/面外"）：
  - 面内极性（模有效电荷 > 0.1 e，如 MoS₂ 的 E′）：长程 Fröhlich 耦合由 POP 管；
  - 面外极性（面外模有效电荷 > 0.02 e，如 MoS₂ 的 A2″）：二维 POP 不含它，但有 σh 的材料里它对谷内一阶耦合是对称性
    禁戒的（D≈0），算了也是白花两个单点。[0035-r3] 之前只按总模有效电荷判，A2″ 这种面外偶极弱的模被当成非极性算了进去。
  要算：S7.2 的 step.conf 写 `ODP_POLAR = include`（阈值 `ODP_POLAR_TOL` / `ODP_POLAR_Z_TOL` 可改）。

**核对**（不过就报错，不静默出数）：
- S7.2 用的 INCAR（物理键）/KPOINTS/POTCAR 与 S7 undeformed 现在的一致（S7 重新生成过就要 retry S7.2）；
- 参考构型与 S7 undeformed 在各谷 k 点、同一条带的本征值相差 < 5 meV（带号、自旋、网格都对上了）；
- 每个构型真空区平面平均势平坦（std < 1 meV）；每个单点电子步收敛（没收敛的本征值是噪声）；
- S5 现在的频率与生成位移时的一致（差 > 2% 告警，> 5% 或模数变了报错要求 retry S7.2；S5 重算会触发 S8.21 重读）。

**和手填通道的关系**（[0035-r2] 起**算出来的优先**）：有冻声子参数时，`DP_OPTICAL_E/H` 里手填的谷内通道（`X>X`）
**只作对照**：gen 日志逐谷打印"手填 … | 冻声子 …"，记进 2d_correction.json 的 `dp_optical.<载流子>.frozen.manual_reference`，
不进散射率（谷名写错照样报错）。所以 step.conf 里以前抄的文献谷内值不用删。要用手填（文献）值：写 `ODP_FROZEN = off`。
手填的谷间通道（`X>X'`、`X>Y`）与冻声子谷内通道一起用。手填了谷表 `DP_OPTICAL_VALLEYS_*` 的载流子不用冻声子参数
（它是按 S8.2 谷表算的），照旧用手填通道。
S8.2 同一组 k 上按带/自旋分出的几个星（S8.4 并成一个谷）：各星的 D² 取平均，相差 > 10% 打 WARN。
`ODP_FROZEN = off` 不用，`on` 必须有（没有就报错）。`ODP_FROZEN_DMIN`（出厂 0.05 eV/Å）以下的模不进通道。

**自动重排**：S7.2 重新生成 → S8.21 重读；S8.2 重算 → S8.21 重读；S8.21 的参数指纹变了（或刚出来、S8.4 还没用上）
→ S8.4 的完成标记归档、重新排队。S8.4 正在排队/运行时只打 ★ 提示。

**计算量**：1 + 2 × 非极性光学模数个单点，单点成本与 S7 的一个形变构型相同。MoS₂：A1′ + E″（二重）= 3 个模，7 个单点
（E′ 面内极性、A2″ 面外极性不算）。
超过 `ODP_MAX_JOBS`（出厂 40）报错；可用 `ODP_MODES` 只算一部分模（编号见 `odp_manifest.json`）。

**局限**：只有零阶（q 无关）耦合；一阶（∝q）耦合与谷间参数都没有。SOC 下简并组按带能量 5 meV 判，劈裂更大的
自旋分裂带按单带算，目标谷态密度含两条带，自旋守恒通道可能偏大（同 [0033] 的 SOC 告警）。
**LOCPOT**：S7.2 各单点的 LOCPOT 在 **S8.21 第一次跑完之前**不能清。S8.21 把各构型的真空能级记进
`step8.21_odp_read/vacuum_levels.json`（按 vasprun.xml 的 sha256 认），之后 S8.2 重算触发的重读用缓存；
vasprun.xml 要留着（本征值从它读）。

## τ(E) 与迁移率等效 τ（[0035]，出厂开，要 WRITE_MESH）

作业链在 transport.json 落成**之后**跑 `tau_energy.py`（它超时或被杀不会丢掉 AMSET 结果），同一份 mesh 出 `tau_energy.json`：

- **τ(E)**：逐机制、按能量分箱（从带边量起，箱宽 10 meV，到 0.5 eV）。τ_m(E) = 1/⟨Γ_m⟩_E，即箱内各态 τ 的调和平均。
  只统计 FD 窗口内 AMSET 真正算过散射的态。
- **τ_eff**：常数 τ 的玻尔兹曼输运要复现 AMSET 迁移率所需的 τ，τ_eff = τ0·μ/μ_CRTA(τ0)，两个迁移率用 AMSET 同一个函数、
  同一费米能级算，分 xx / yy / 面内。**τ_eff,ADP 与 S8.2 DPT 的 τ（S8.1 BT2 用的那个）同口径**，可以直接比。
- **核对**：重算的 overall 迁移率必须复现 transport json（不复现就报错）；逐 k 点电导加权的 τ_σ 与 τ_eff 相差 > 10% 时写 ★ 注记。

失败只影响 `tau_energy.json`（写 status = error 和原因），不影响 transport.json。`TAU_ENERGY = false` 关掉。
S8.3 出 `tau/`：`tau_comparison.csv`（各温度、各掺杂档的 τ_eff 逐机制与 DPT τ、ADP/DPT 比值）、`tau_vs_E_300K.png`、
`tau_eff_vs_T.png`、`tau_summary.txt`。

## 已知局限（写论文要写进方法学）

- **（2026-09-28 汇总，第 9–13 项）未验证 / 模型局限**：
  - 插值不确定度：MoS₂ 在验证过的 263×263×21 网格上空穴仍有 **8–10%** 的插值不确定度（V112）；
    出生产结果前做 factor 收敛 f = 61 / 90 / 120（f=120 约 195 GB，gen 的内存粗估会告警，先确认节点内存）。
  - 2D POP 偏弱（μ_POP 约 1.9e4）**未验证**；kz=3 对照没做。
  - DPT 拟合有效质量的窗口随网格变（换网格要重看 m*）。
  - **没有 SOC 路径**（Bi₂Te₃ / PbTe / SnSe 这类都需要）。相位补丁的 SOC 分支用同一因子、未测；
    且 AMSET 的 ncl 去对称化对 TR 操作不取共轭、不乘 iσ_y —— SOC + 无反演一律要求全网格（V115 §5）。
  - AMSET 不含非极性光学声子与谷间散射（见下条），Si 空穴、MoS₂ 迁移率偏高；PBEsol 有效质量偏轻，
    GaAs 迁移率偏高。写进方法学即可，不是流程 bug。
- **"去掉 ZA 支"只对有水平镜面（σh）的材料成立**：翘曲结构（硅烯、锗烯）与 Janus
  结构存在一阶 ZA 耦合，对它们是模型遗漏（ZA 是二次色散，形变势模型也处理不了）。
- **缺非极性光学声子（光学形变势 ODP）**：AMSET 实现的机制只有
  **ADP / POP / IMP / PIE**（`amset/scattering/elastic.py` = ADP+IMP+PIE，`inelastic.py` 只有
  `PolarOpticalScattering`；官方文档 Summary of scattering rates 同）。**非极性晶体 POP 前因子为零**
  （`scattering/calculate.py:842` 注释明说 "for POP when studying non polar materials"），
  于是这类材料**整条光学声子通道缺失**：
  - **2D**：MoS2 的同极 A1 模与谷间光学模耦合很重要 —— **总**迁移率系统偏高，
    与文献对照（T4）必须**分机制**做，不比总迁移率。
  - **3D 非极性半导体（Si、Ge、金刚石这类）**：**空穴迁移率被显著高估，电子受影响小得多。**
    实测（Si，300 K，1e17，AMSET 0.4.19 本地受控，详见 `VERIFICATION.md` §六点一）：

    | 分支 | AMSET | 本征声子限实验值 | 倍数 |
    |---|---|---|---|
    | 电子 ADP | 1622.7 | ≈1400 | 1.16 |
    | **电子 overall** | **805.0** | **≈800** | **1.01** |
    | 空穴 ADP | 5434.8 | ≈470 | **11.6** |
    | **空穴 overall** | **1474.3** | **≈330** | **4.5** |

    电子整体几乎完全吻合（1e17 实验值 ≈800 ✔）；空穴高 4.5 倍 —— 与"缺一项主要散射机制"
    的预期一致。空穴的声子限迁移率里**光学声子（光学形变势 d₀）贡献很大**，正是缺的那一项。
    价带形变势本身没问题（`deformation.h5` 里 VBM D(Γ)=1.74 eV，接近文献 a_v ≈ 2 eV），
    **不要往形变势归一化上找原因**。
    > 上游对照（2026-09-24 首记，**2026-09-27 更正出处**）：AMSET 官方 `examples/Si` 的 `doping`
    > **全是负值（纯 n 型）**，它本身**不产出空穴分支**；官方算例复跑只得到电子 ADP=**1360.1**
    > （≈论坛里的 1300 ✔）。所以"官方算例 空穴 ≈5400"**不是官方算例的原样结果** ——
    > 那个 5400 是 [matsci.org 论坛帖](https://matsci.org/t/error-in-the-hole-mobility-of-silicon-in-the-example-of-asme/59809)
    > 里的讨论值。**本流程自己算出的 Si 空穴 ADP 是 5000–5435**
    > （全网格臂 5003.7、去对称化臂 5434.8），与论坛量级一致。
    > **结论不变**：空穴偏高是 AMSET 模型本身特性（缺非极性光学声子/光学形变势散射），
    > 不是本流程的设置问题。
    >
    > ★ 上表数字**必须标来源臂**：`1622.7 / 5434.8` 出自**去对称化 h5** 臂；
    > **全网格 h5** 臂是 `1618.9 / 5003.7`（同掺杂 1e17、300 K）。两者不可混用。
    > 另：曾有一版把 n/p 对调后得出"电子/空穴误差都恰好 +258%"——那是**写反了**，
    > 负掺杂才是电子（ADP 1611–1619）、正掺杂是空穴（5004–5435）。
  - **写论文/报数时**：用 AMSET 报**非极性**材料的**空穴**迁移率，必须标注缺少光学形变势散射；
    极性材料（MoS2 等，有 POP）不受此条影响。
- 偶极近似，**没有四极矩**：Poncé 等（PRL 2023）表明准确预测二维迁移率需规范协变四极矩，
  **MoS2 具体数字**：Poncé PRL 2023 —— 室温 Hall 迁移率误差 电子 23%、空穴 76%。
- 极性散射只用**单一有效声子模**。
- 忽略 Janus 结构的面外偶极耦合。
- 自由载流子屏蔽用局域场 + 二维 Thomas-Fermi（文献对比偏低 10-20%）。
- 要求层法向平行于笛卡尔 z（插件会检查）。
- **形变势参考口径**：2D 默认真空（`deformation_vac.h5`，与二维文献、C2DB 一致）；
  该文件不存在时退回 core 口径（AMSET 芯态对齐，为体材料设计），两者绝对值差约 2.8 倍。
- **绝对归一化的验证状态**（G_s = c·G_2D 本身仍是推导，靠独立解析参考支撑）：
  ADP **-2%**、POP **散射率中位比值 1.014±0.025（μ -9%）**、IMP **-16%**（差异归因于屏蔽口径：
  参考用非简并 q_TF = 2πn/kT，插件用态密度积分），**PIE 未做定量验证**（默认关闭；
  屏蔽长度用静态 ε0，配离子弛豫的压电总张量）。
  与文献比总迁移率时必须**分机制**，且记住模型缺 ODP / 四极矩，总迁移率会系统偏高。
- **★★ 并行度（nworkers）不影响结果 —— 早期"nw=24 偏高 15.5%"已撤回（2026-09-27）**
  同一套输入/掺杂/温度（Si 输入、f=10、117³、T=300、doping = -1.99e14 … -4.39e18），
  只变 nworkers，1→24 共 8 次（含 **nw=24 重复两次**）：

  | nworkers | 电子 ADP @ -1.99e14 | ADP @ -1e18 | overall |
  |---|---|---|---|
  | 1 / 2 / 8 / 12 / 16 / 23 / **24a** / **24b** | **1394.48** | **1379.31** | **1241.19** |

  **八次逐位相同**（驱动 tmp/v1_bench/run_nws2.sh，结果 tmp/v1_bench/NWS2/*）。
  若分块/并发真有 bug，不可能做到逐位相同。
  **撤回说明**：早期那个 "nw=24 -> 1610.6（+15.5%）" 用的是**另一套掺杂与温度**
  （10 个 doping 点 ±1e21…、9 个温度 100–900 K），而且那次跑**中途被中止**、没留下
  最终输运文件；与 nw<=23 那组（6 个 doping 点、T=300）**本来就不可比**。
  ⇒ **"用 nw=24 跑出的 3D 电子迁移率偏高 15%" 这条作废**；nworkers 维持"按提交分配核数"。
  教训：参数对照前必须先 diff 两份完整 settings.yaml（工具 tmp/amset2d/cmp_settings.py），
  只盯着要变的那个参数是不够的。
- 第二行横幅的 r∞/Δr 是输入参数量级自检（VERIFICATION V14）：锚点 MoS2 r∞≈46.5 Å / Δr≈0.53 Å（1.1%）、
  h-BN r∞≈7.6 Å / Δr≈2.8 Å（37%）；老记录缺该字段会打 [WARN]。

- **★★ 真实重叠的裁决已定论（2026-09-17，见 VERIFICATION V23/V24）：根因是去对称化**
  - **默认配置**：`settings.yaml` 用 **`unity_overlap: true`**（生成器现已**显式写入**，
    以前根本没写 → AMSET 默认 False=真实重叠，所有 2D 项目都默认踩坑，见 V24.1）。
    已知偏差方向：多能谷材料的 ADP 迁移率**最多被低估 N_v 倍**（MoS₂ K+K′ 为 2）。
  - **要真实重叠就必须开 `wavefunction_full` 分支**（S3b 用 `ISYM=-1` 出全网格 WAVECAR，
    S4b 出 6627 点的 h5），让 AMSET 走 `from_data`、跳过 `desymmetrize_coefficients`。
    实测（V23.1，输入逐字一致、唯一变量是 h5）：

    | 配置 | μ_e | μ_h | 对判据带 205/102、1370/685 |
    |---|---|---|---|
    | 旧 h5(416) + 真实重叠 —— 0.4.19 / 0.5.1 | 1056.8 / 3341.2 | 8595.8 / 15802.1 | **超 5~16 倍** |
    | **全网格 h5(6627) + 真实重叠 —— 0.4.19 / 0.5.1** | **250.4 / 280.3** | **1070.3 / 1148.2** | 电子 1.2~1.4 倍，**空穴带内** |
    | 两种 h5 + unity（对照） | 99.3 / 111.1 | 569.3 / 526.6 | 逐位相同 → 变量隔离干净 |

  - **插值已排除**（V23.3）：step3c 的 14 个离网格点（含 K 点本身、kz=1/6、Q 谷）直检，
    `|⟨c_AMSET|c_DFT⟩|²` **中位 0.9936** → 不需要加密网格。
  - **代价**：全网格分支的 k 点数 ×对称操作数（MoS₂ 16 倍），h5 85 MB → 1.3 GB。
  - **残差（未解释清）**：电子仍是单谷参考的 1.2~1.4 倍、ADP 比值 2.52 略超 N_v=2。
    在解释清楚前，真实重叠的**绝对值仍不宜直接写进论文**。
  - **检查与报警**（V24 / V25.10 / V26，2026-09-18 更新）：
    ① 运行前检查 `overlap_preflight.py`（`step8_amset` / `step8.4_amset2d` 的作业里都会先跑，不过就 `exit 1`、不出数）：
       ⓪ **2D 且 `unity_overlap != true` -> 拦截**；但若 `step.conf` **显式**写 `UNITY_OVERLAP = false`（受控对照），
          gen 会在 settings 里写 `AZ_OVERLAP_CONTROLLED=1`，此时降为**告警放行**（结果只作对照、不作生产）。
       ① AMSET 版本 + 重叠模式；
       ② h5 完整性：k 点数 == 网格乘积（= `from_data`、不去对称化）；真实重叠 + 非完整网格：
          **2D 拦截 / 3D 告警**（Si 全网格对照已完成：ADP/overall 只差 3~9%、IMP 逐机制 1.3~1.7x，见 VERIFICATION V26）；
       ③ S3 vs S3b：NBANDS + 本征值一致性（<1 meV，窗口取运行日志实际值）；
       ④ `amset wave` 带窗口必须等于 AMSET 实际插值的窗口；
       ⑤ 弹性张量：**面内** Christoffel 特征值比 `<=0` 或 `<1e-2` -> **拦截**（多为 VASP 顺序未重排、C44/C66 互换）；
          **仅面外**为负（二维 slab 真空伪影）：本步（插件，ADP/PIE 只用面内 2x2）**只记录、不告警**（V180）；
          `step8_amset` 标准路径用完整 3x3 -> **告警**（gen 会把**负的**面外剪切取绝对值，原值->新值记进
          `2d_correction.json` 的 `elastic_outofplane_shear_zeroed`；还看到负值多半是旧 settings）。
    ② `use_projections: true` 的运行**不适用**比值报警：投影重叠不满足 `|I|^2<=1`（实测比值可 0.71/0.41），guard 遇到它直接跳过。
    ③ 8.3 自动算比值（`overlap_ratio_guard.py`，2026-09-19 分级）：ADP **`<1` 红 / `[1, N_v]` 绿 / `>N_v` 黄（要求人工确认）**；
       **N_ch（按谷区分实测，同谷=1/谷间=0）只作解释打印，不再当放行阈值**（SS 实测 N_ch 电子 ~23、空穴 ~250，见 V29）；材料目录放 `valley_ratio.json`（`{"n_ch": 3.13}`）会被读到并打印。
       **POP 只记录、不参与红黄**（真实重叠也压低谷内通道，±20% 会误红）。实测：SS 修正弹性顺序后面内 1.06 -> 绿；MoS₂ 2.52 > N_v=2 -> 黄（人工确认）。
    ① 运行前检查 `overlap_preflight.py`：⓪ **2D 且 `unity_overlap != true` → 直接拦截**
    （这条不读波函数，gen 时即生效，正是能挡住原始 bug 的那条）；① AMSET 版本 + 重叠模式；
    ② h5 完整性（k 点数 == 网格乘积 = `from_data`；真实重叠 + 非完整网格：**2D 拦截、
    3D 只告警**——三维结论待 Si 全网格对照，`BLOCK_3D_REAL_OVERLAP=False`）；
    ③ S3 与 S3b 的 NBANDS + 本征值一致性（<1 meV）；④ `amset wave` 带窗口**必须等于
    AMSET 运行日志里实际插值的窗口**（从 `Interpolating spin-up bands L-H` 读，各材料不同：
    MoS₂ 11-17、CrSe₂_Hex 9-15、Pb₂Bi₂Te₅ 36-51）。不通过直接拦截。
    ② gen 是在**集群上**跑的（cwd 里 vasprun.xml / h5 / amset.log 都在），所以上面四项在
    gen 时就是实效的；`step8_amset` / `step8.4_amset2d` 的作业里再跑一次
    `--in-job` 才开算，检查不过就**不出数**。
    ③ 8.3 输出里自动算 `overlap_ratio_guard.py` 的比值，ADP 分级
    **[1, N_v] 绿 / >N_v 黄（人工确认）/ <1 红（拦截、不得进 zT 汇总）**，POP 只记录。
  - **现象**：2×2 实验（输入与集群 run 逐字一致）显示，真实重叠把 ADP 迁移率整体抬高
    **电子 ×10.6、空穴 ×15.1**，而插值因子 4→10 只值 ×1.06~1.11。
    同输入下 `unity_overlap: true` 给 99.3 cm²/Vs（与 Takagi/2 参考 102 一致），
    `false` 给 1056.8（集群 1067.6）。
  - **判据带（口径已修正，见 V22.1）**：二维 Takagi 单谷 **电子 205 / 空穴 1370**，
    两谷且重叠取 1 时再除以 2（102 / 685）。回归锚点：TB 台子算得 480.7，与
    `ref2d_tb.py` 的 ADP 自检 472.5 差 1.7%（`tmp/amset2d/takagi2d.py`）。
    → 真实重叠两版（1057 / 3341）都远超上限。
  - **根因：尚未定论（见 V22）**。曾怀疑 `desymmetrize_coefficients`，但那是**错误度量**造成的
    假阳性（h5 里同一 k 点各带系数并不正交，Gram 非对角达 0.18，连恒等路径都会给出 0.658）。
    改用 QR 正交化后的**主角度**重做：**原证据作废**。但**不能说"去对称化没问题"**——
    那 416 个点本身就是 DFT 算出的不可约点，AMSET 给它们几乎都选了**恒等操作**，
    该检查**没用到 TR 分支**。正确表述：**度量有误、证据作废；去对称化是否正确仍未确定，
    时间反演分支可疑**。依据：379 条非恒等路径失败 24%（全是 TR），而全网格里
    **46.9% 的格点由 TR 分支生成**（K 谷附近 40.3%）。
    非选中路径只有 75.7% 正确 → 该变换**不满足路径独立性**；
    6211 个非 DFT 格点是否可信**本地无法裁决**。
  - **最小复现**：`tmp/amset2d/desym_final.py`、`desym_altpaths_truth.py`、`desym_pairs.py`、`desym_validate.py`。
  - **下一步（已获用户批准动作；但步骤映射更正见 V22.8）**：
    **ISYM = -1 要设在 S3_uniform，不是 S4** —— S4 只跑 `amset wave`、没有 INCAR；
    S3 的 `KPOINTS` 已经是完整的 Γ 中心 47×47×3，只是 `ISYM = 2` 把它约化成了 416 点。
    `ISYM = -1`（**不能是 0**：0 仍用时间反演，只给半个网格）后 k 点 → 6627，
    h5 约 85 MB → 1.3 GB，内存需留意。裁决表见 V22.4。
    要在不覆盖原步骤的前提下落地，**需给技能加一对新步骤**（属受控改动，待确认）。
  - **升级到 0.5.1 不解决问题（V21，已实测）**：0.5.1 确实补上了 0.4.19 缺的
    `g_diff`/`g_maps` G 平移（源码对照成立），但**同输入 A/B** 显示真实重叠的迁移率
    被**再抬高 3.0~3.2 倍（电子）/ 1.8 倍（空穴）**——电子 1057→3341、空穴 8596→15802，
    离判据带更远；对照组（unity 重叠）两版一致到 ±12%。
    **因此不要为二维项目在集群上升级到 0.5.1。**
  - **「只是缺一个 G 平移」也不成立**：把平移搜索扩到 ±4，30 个失败样例只救回 10 个。
  - **后果**：真实重叠路径把散射率系统性压低约 10~30 倍 → 迁移率被抬高同样倍数。
  - **在 `desymmetrize_coefficients` 修好之前，S8.4 的绝对数值不能写进论文**
    （用户要求，继续有效）；二维项目固定 `unity_overlap: true`。

要作为定量结论发表，建议选 1-2 个材料用 EPW / Perturbo（开二维库仑截断）对照。
