# step8.4_amset2d -- AMSET 二维散射核（插件式）

> 只在 2D 项目里打开；三维项目走原版 `step8_amset`，完全不受影响。
> 开关：项目 `project_setting` 里写 `amset2d: true`（技能里该组 `default: false`）。

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
| `amset2d_plugin.py` | 运行期插件：把 ADP/POP/IMP/PIE 四类换成二维核（同名注册表替换，不改 AMSET 安装） |
| `amset2d_pop_freq.py` | 取 Gamma 点**面内极性模**有效频率（AMSET 自带口径会把面外 ZO 模算进去） |
| `test_kernels.py` | 四类核的自检（旋转不变、极限行为、护栏），`python test_kernels.py` 应全 PASS |
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

## 已知局限（写论文要写进方法学）

- **"去掉 ZA 支"只对有水平镜面（σh）的材料成立**：翘曲结构（硅烯、锗烯）与 Janus
  结构存在一阶 ZA 耦合，对它们是模型遗漏（ZA 是二次色散，形变势模型也处理不了）。
- **缺非极性光学形变势（ODP）**：MoS2 的同极 A1 模与谷间光学模耦合很重要，AMSET 没有；
  因此**总**迁移率会系统偏高 —— 与文献对照（T4）必须**分机制**做，不比总迁移率。
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
  参考用非简并 q_TF = 2πn/kT，插件用态密度积分），**PIE 未做定量验证**（默认关闭）。
  与文献比总迁移率时必须**分机制**，且记住模型缺 ODP / 四极矩，总迁移率会系统偏高。
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
  - **检查与报警**（V24 / V25.10）：
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
    **[1, N_v] 绿 / (N_v, 2N_v] 黄（人工确认）/ > 2N_v 红（拦截、不得进 zT 汇总）**，
    POP 另按 1 ± 20%。
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
