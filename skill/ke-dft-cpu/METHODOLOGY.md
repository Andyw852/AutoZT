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

### 4.4 验证状态（详见 `step8.4_amset2d/VERIFICATION.md`）

| 项 | 结论 |
|---|---|
| 四种核的极限行为、旋转不变、护栏 | 合成数据单元测试全 PASS（amset 0.4.19 / 0.5.1） |
| Voigt 重排 | 真实 OUTCAR 实测；护栏能拦下未重排的输入 |
| ADP 归一化（不乘 c/t） | 与闭式 2D Bardeen-Shockley 对照：电子 302 vs 241 cm2/Vs（+25%）；旧流程偏大 3.1 倍 |
| 真空不变性（映射系数 c） | c = 20/30/40 Å 三次实跑，ADP/总迁移率变化 < 10%（低掺杂电子档最大 26%，来自 BoltzTraP2 插值网格随 c 变化，无法完全隔离） |
| 真实材料端到端 | CrS2 单层 / CrS2 正交胞 / CrSe2 单层跑通（amset 0.4.19 = 集群 amset_clean 版本） |
| MoS2 本身 | 本机与集群都没有 MoS2 的 ke-dft-cpu DFT 数据（只有 MACE 优化/声子），未做；同族 TMD 单层已覆盖 |

## 5. c 统一（比值法前提）

各材料元胞 c 轴统一 20 Å（有意为之），AMSET σ 的跨体系比值不含 c 差异伪影。
IMP/POP 扣真空、σ 归一化都依赖此前提。

**8.4（amset2d）之后这条限制放宽**：二维核里的 `c` 因子已把真空长度显式积掉，
真空不变性实测 < 10%（同为 20 Å 的跨体系比值仍是更干净的比较口径，建议继续保留）。

## 6. skill_rev 防呆

step8.x 的 `_SKILL_REV` 写进 comparison_summary.txt。改脚本即提交 GitHub 并 bump，
防陈旧副本（step12 / step9b / step13 三次被同一形状的 bug 咬过）。

---

> 结果文档（未发表，不进 git）：`/mnt/d/tf_data/jzz/jap/comparison_vs_literature.md`
