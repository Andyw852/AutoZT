# AutoZT 论文组织方案（2026 秋学期）

> 写于 2026-09-22。用途：① 回复课题组的"本学期投稿梳理"；② 把图 1（系统总览）
> 与图 2（kl 技能 + 状态机）背后的实验补齐，锁定投稿窗口。
> 事实来源：`docs/EVALUATION.md`、`docs/RELATED-WORK.md`、`docs/RELATED-WORK-EXTENDED.md`、
> `docs/HANDOVER-2D-kl-20260917.md`、`README.md` 与两个 benchmark 数据集目录。

---

## 0. 一句话

以 **AutoZT** 为成果：把"一个 POSCAR → `zT(n,T)`"的 VASP/MACE 热电计算做成**跨异构超算的生产级执行底座**，
下一步由物理判据推进、动作经风险分档与人工闸门、产物按文件哈希归档、失败按状态分类恢复；
LLM 只在异常路径介入且可被替换。**2D/3D 验证材料集已建好但还没跑**——这是本学期剩下的主要工作量。

---

## 0.5 ★ 状态更正（2026-09-22 晚，逐项复验后改写）

初稿的"缺口"一节是照 `HANDOVER-2D-kl-20260917.md` §0/§28.1/§29.4 写的，而那几处**已经过期**。
本次实际复验结果：

| 原以为的"硬阻塞" | 实际状态 | 证据 |
|---|---|---|
| phonon-mlff 第三副本 ZA 判定是错的 | ✅ **已修**（2026-09-22 迁移到 `_common/za_2d.py`；三副本合并完成） | 本次独立复跑真实 MoS₂ fc2（md5 `acfd37b2…`）：`vac(2)→prim=0`、`dirs=[(0,1,0),(0,1,-1)]`、qmax=0.05 → `p=(1.9930749384721, 1.9801186415610)`，与 §28.2 验收数字**逐位一致**；`imag_mesh_numbers→[1,60,60]`；`band_path_2d_prim→hexagonal Γ-M-K-Γ` |
| S6 的"网格轴序"bug 未修（MoS₂ κ 出不来） | ✅ **已修**（`gen_step6_kappa._mesh_prim_remap`，commit `a0a8e57`） | 已产出**原胞轴序**网格的 κ：`kappa-m18888 / m1110110 / m1138138`（= "1 88 88"/"1 110 110"/"1 138 138"），RTA 三档收敛（§31.2） |
| BHH 390/520 真数据对照没做 | ✅ **已做** | `tmp/mos2_fc_compare/MoS2_ZA_refix_report.md` §3a：390 p=(1.99307,1.98012)、520 p=(1.99307,1.98014)，Δp≈1e-5；两套 pheasy 日志均含 `Imposing rotational invariance and equilibrium conditions` |

**因此 Tier A 的门槛不是"修 ZA"，而是下面三件（都已在交接文档里，需要人拍板或重跑）：**

1. **P0-2 厚度归一化端到端**：代码/逻辑已验证（§39.1），但线上 MoS₂ 的 `kappa_summary.json` 是**修复前的陈旧产物**（§39.2/§39.3）→ 重跑一次 S6 即闭环。
2. **MoS₂ 对称化链卡在 S5**：对称化后重拟合**比未对称化更差**（§46.9、§47）——这是科学口径问题，**需要你决定走哪条**（保留未对称化 / 接受对称化 / 换采样）。
3. **同一虚频量三套阈值**（S5 0.10 / S5.1 −0.05 / lattice_kappa 0.5）互相打架 → §45 第 4 项，**需要你给一个值**。

---

## 1. 群回复草稿（可直接粘贴）

> 本学期主投 **1 篇 leading first author**：AutoZT——面向热电材料的自动化计算与运维平台。
>
> **创新点**
> 1. 单一 POSCAR 到 `zT(n,T)` 的端到端流水线：电子（能带→散射→载流子输运）与声子（稳定性→力常数→声子输运）双链汇合，**3D 与 2D 走同一套工作流**（2D 原生处理真空轴、面内弛豫、二维单位、ZA 支与厚度归一化）。
> 2. **物理判据驱动的自治**：把"不报错但会给出合理错值"的静默错误逐条设成门禁（目前 11 类：BHH 旋转不变性、κ 厚度归一化、面内应力取代 |P|、Pulay 偏差、旧结构帧混拟合、SCF 撞 NELM 等），稳态推进 **0 次 LLM 调用**，成本可预测、每一步可追溯到物理量。
> 3. **LLM 只在失败路径介入**：经 MCP/CLI 操作，动作分只读/推进/破坏性三档，破坏性需人在终端发一次性令牌；planner 可换成脚本或人工而不改变安全边界（架构主张 + 消融验证）。
> 4. 跨异构超算：真 SLURM 的 CPU/GPU 集群 + 无 SLURM 的 GPU 主机，逐材料选集群、DFT 与 MLIP 步骤混跑；per-file SHA-256 provenance + `prove --verify` + 会话导出可重放。
> 5. MLIP 预筛 + DFT 分级 + GPU 力常数拟合（图 2），把"精度—机时"权衡做成可量化的筛选漏斗。
>
> **进度**
> - 平台：20 个技能、3 个集群实跑；单批 2926 个结构（金属掺杂富勒烯）MACE 优化 + 形成能跑通。
> - 硬指标（可复跑脚本）：危险动作拦截 9/9 = 100%、实际执行 0；跨集群 17/17 输入逐字节一致；18 个技能的输入生成 16/18（另 2 个按设计依赖上游）；真实提交的首个可算步骤 15/15 全部 OK。
> - 2D 链路：8 个材料的 S1 结构优化在生产上验证（面内应力全部 ≤ 0.114 kB）；旧结构力帧被门禁拦下（差 0.1263 Å）；力噪声底同节点≈0、跨节点 ≈1e-3 meV/Å。
> - 图 1（v6）、图 2（v13）已出；2D/3D 验证集（POSCAR + 带 DOI 的文献参考值 + 分档）已建好。
>
> **下一步（决定投稿时间的唯一路径）**：把验证集跑出来——3D 先 Si/BAs，2D 先 MoS₂/WS₂（WS₂ 有 pheasy 文献值 241.6 W/mK 可直接比），再扩到 Bi₂Te₃/PbTe/SnSe 的 `zT` 端到端。
> **预期投稿：2026 年 12 月中旬**（保守 2027 年 1 月上旬）；首选 *npj Computational Materials*，备选 *Nature Computational Science* / *Computer Physics Communications*。

---

## 2. 论文定位与题目候选

| 项 | 内容 |
|---|---|
| 定位 | 软件/方法学论文（production workflow + agent governance），不是单材料科学发现 |
| 差异化一句话 | DREAMS / AtomAgents 证明了 agent 栈能规划与校验；AutoZT 的保证**附着在执行底座**（判据、状态机、逐文件哈希、动作闸门、多集群），规划层被替换时仍然成立 |
| 题目候选 1 | AutoZT: physics-gated autonomous workflows from a single POSCAR to zT(n,T) across heterogeneous HPC |
| 题目候选 2 | Deterministic autonomy for thermoelectric discovery: criteria-driven execution with an LLM at the failure boundary |
| 目标期刊 | 首选 *npj Computational Materials*（`docs/AutoZT_software_overview_npj.docx` 已在按此格式准备）；备选 *Nature Comput. Sci.* / *Comput. Phys. Commun.* |

**论文中要避免的表述**（`RELATED-WORK-EXTENDED.md` §3）：不说 DREAMS 无 HPC/安全/provenance/确定性；
不说 AtomAgents 无审批；不把"文献未讨论"写成"系统不具备"；不把"可重放"写成"浮点 bit-exact"；
不把"支持 3 个后端"扩大成"适用于所有超算"。

---

## 3. 创新点 ↔ 图 ↔ 证据

| # | 创新点 | 图 | 现有证据 | 缺什么 |
|---|---|---|---|---|
| 1 | 单 POSCAR → `zT(n,T)` 双链汇合，3D/2D 统一 | 图 1a/1b | README 快速开始；技能契约；2D 厚度/真空/ZA 处理 | 2D 与 3D 各一组端到端 `zT` 实跑 |
| 2 | 物理判据门禁（静默错误拦截） | 图 2 状态机 + §14 表 | HANDOVER §14 的 11 行表，多数"已生产验证" | 表 1 行需真数据：BHH 的 390/520 对照、S6 厚度归一化端到端 |
| 3 | LLM 只在失败路径，动作分档 + 人工令牌 | 图 1c | `EVALUATION.md`：拦截 9/9、执行 0、三档 profile；`agentgate.py` + 40 项测试 | planner 替换消融（脚本 / 不同 LLM / 人工 三种 planner）|
| 4 | 跨异构超算 + 逐文件 provenance | 图 1a 底部 | 17/17 跨集群输入逐字节一致；`prove --verify`、`session export` 已有测试 | 规模—成本实验（材料数 × token × 墙钟 × 人工介入）|
| 5 | MLIP 预筛 + DFT 分级 + GPU 拟合 | 图 2a | fc-fit（Si，10 帧×250 原子）产出 fc2/fc3 + ShengBTE 导出；phonon 判据稳定 | pheasy-gpu 的实测加速比（图 2 上标注的 "×?" 还没填）；有限位移 57× vs 随机位移 21× 的对照 |

---

## 4. 论文结构（对应图号）

| 图 | 内容 | 数据来源 | 现状 |
|---|---|---|---|
| 图 1 | 系统总览：POSCAR→zT(n,T)、双链+验收门禁、3D/2D 统一、LLM 经 MCP/CLI + 动作三档 | 架构 | ✅ v6 已出 |
| 图 2 | kl 技能：MLIP 预筛→DFT 力→力常数拟合→输运 + 显式状态机/挂死恢复/风险分级 | `autozt skill show kl-dft-cpu`（卡片由 skill.yaml 自动渲染，不会过期） | ✅ v13 已出，**"×?" 加速比待填** |
| 图 3 | 安全与治理：动作分档、一次性令牌、审计流水；gateway 开/关消融 | `scripts/safety_metrics.py` | 数据已有（9/9、三档 profile），待作图 |
| 图 4 | 资源与精度基准：method × backend × fit 的机时/精度表 | 需新增 `autozt bench` 聚合（每步 `resource_summary.json`） | ❌ 需真实跑批数据 |
| 图 5 | 筛选漏斗与排名（te-screen 各阶段剩余数、淘汰原因） | `screening_report.json/csv` | ❌ 需 te-screen 实跑 |
| 表 1 | 技能 × 软件栈 × 硬件（CPU/GPU） | `docs/AutoZT_skill_and_software_tables*.docx` | ✅ 三版已出 |

---

## 5. 进度盘点（已做，可复跑）

| 类别 | 数字 | 出处 |
|---|---|---|
| 技能 | 20 个（VASP + MLFF + 辅助），3 个集群实跑 | `skill/`、README |
| 批量实战 | 2926 个金属掺杂富勒烯结构（opt-mlff-cpu） | README |
| 安全 | 危险动作拦截 100%（9/9），执行 0，只读仍可用；compact 7 / workflow 12 / monitor 2 / full 24 工具面 | `scripts/safety_metrics.py`、`EVALUATION.md` |
| 复现 | 跨集群 17/17 输入逐字节一致；本地输出确定性 100% | `pa -tt <技能> -p <材料> prove --verify` |
| 技能覆盖 | 18 技能输入生成 16/18；首个可算步骤真实提交 15/15 OK | `EVALUATION.md` §3 |
| 力常数 | fc-fit：Si 10 帧 × 250 原子 → fc2.hdf5 / fc3.hdf5 + ShengBTE 导出 | `EVALUATION.md` §3 |
| 2D kl | 8 材料 S1 生产验证（面内 ‖σ‖ ≤ 0.114 kB）；Mo₂S₃ Z4-3-1 旧结构被拦（0.1263 Å） | HANDOVER §10.1、§14 |
| 噪声底 | 同节点 ≈ 0；跨节点 ≈ 1e-3 meV/Å（比待分辨的口径效应小 2–3 个数量级） | HANDOVER §14 |
| ZA 判定 | 三副本合并到 `_common/za_2d.py`，真实 MoS₂ fc2 上 p=(1.99307, 1.98012)，390/520 一致 | 本次复验 + HANDOVER §38 |
| MoS₂ κ | RTA 三档网格（1×88×88 / 1×110×110 / 1×138×138）跑完且收敛 | HANDOVER §31.2 |
| BHH | 390/520 两套真实 fc 上均确认 RASR=BHH 已施加，Δp≈1e-5 | `MoS2_ZA_refix_report.md` §3a |
| 图 | 图 1 v6、图 2 v13 | drops |

---

## 6. ★ 2D/3D 验证方案（本学期核心工作量）

### 6.1 现成资产：验证集**已经建好，但一个都还没跑**

| 数据集 | 路径 | 内容 |
|---|---|---|
| 2D 一档（严格） | `/mnt/d/tf_data/work_AutoZT/2D_transport_benchmark/validation_tiers/01_strict_standard_benchmarks/` | graphene、hBN、phosphorene、MoS2_2H、MoSe2_2H、WS2_2H、WSe2_2H、SnS、SnSe、Bi2Te3_QL（每个含 POSCAR + 结构审计 + 文献 DOI） |
| 2D 二档（扩展） | `…/02_extended_literature_2d/` | GaSe、GeS、GeSe、InSe、ReS2、TiS3 |
| 2D 三档（回归） | `…/03_database_bulk_regression/JARVIS_dft2d/` | 1103 个数据库结构（只做解析器/数值稳定性回归） |
| 2D 参考值 | `2D_transport_benchmark/transport_properties.csv` | 24 条：kl 12、S 3、σ 3、ke 3、ZT 3（带 DOI + `value_note`） |
| 3D 声子标准 | `…/bulk_transport_benchmark/01_standard_phonon_benchmarks/` | Si、Ge、diamond、GaAs、LiF、MgO、NaCl、BAs |
| 3D 热电标准 | `…/02_standard_thermoelectric_benchmarks/` | Bi2Te3、CoSb3、Mg2Si、PbTe、Sb2Te3、SnSe |
| 3D 参考值 | `bulk_properties_long.csv` | 14 条：kl 8、ZT 6（带 DOI + evidence_level） |

**关键说明（写进论文的 validity 一节）**：2D 的 W/m/K 依赖各文献自定的有效厚度，
2D κ 值必须按来源口径归一化后再比；公开数据里**不存在**一套对全部一档材料用完全相同设置、
同时给出 S/σ/κe/κL/ZT 的统一结果表 → 因此精度验证采用"**统一程序重算 + 文献值交叉检查**"两轨制。

### 6.2 三档跑批计划（先窄后宽，每档都有明确验收）

**Tier A —— 打通与保底（必须先绿，约 2 周）**

| 材料 | 目的 | 验收 |
|---|---|---|
| 3D **Si** | 最简单、已有大量历史数据（opt/kl/phonon/ke/TE 全技能都跑过） | κ(300K) 对 148 W/mK；`zT` 链全程无人工干预 |
| 3D **BAs** | 超高 κ 极端值，检验拟合与网格是否被推到极限 | κ 与文献量级一致；记录失败模式 |
| 2D **MoS₂（2H 单层）** | 线 A 已启动；ZA 支必须打开 | 面内 κ + ZA 支判据；厚度归一化记录 |
| 2D **WS₂（2H 单层）** | pheasy 论文基准 **241.6 W/mK**（d=6.16 Å、90×90×1），可直接比 | 同设置重算 κ，误差报告 |

**Tier B —— 论文主精度表（3–4 周）**

| 组 | 材料 | 用途 |
|---|---|---|
| 3D 声子全量 | Ge、diamond、GaAs、MgO、LiF、NaCl（+ Si、BAs 共 8） | 小胞低成本，构成 κ 精度的主证据表 |
| 3D 热电端到端 | **Bi₂Te₃**（室温 ZT≈1）、**PbTe**（中温） | 电子（AMSET）+ 声子（phono3py/pheasy）双链汇合到 `zT(n,T)` |
| 2D 热电端到端 | **Bi₂Te₃_QL**、**SnSe**、**SnS** | 2D 原生链路（真空轴/面内弛豫/厚度换算）的 `zT` 验证 |

**Tier C —— 泛化与外推（可选，投稿前能跑多少跑多少）**

| 组 | 材料 | 用途 |
|---|---|---|
| 2D 二档 | GaSe、GeS、GeSe、InSe、ReS2、TiS3 | 覆盖更多对称性/化学环境；须先匹配文献的相与厚度 |
| 数据库回归 | JARVIS dft_2d 1103 结构 | 只跑便宜步骤（opt + 静态），验证解析器与数值稳定性、吞吐与失败率 |

### 6.3 成本与排序原则

- **先小胞、先声子、后电子**：声子基准（2–4 原子/胞）机时低、判决快；AMSET 电子侧最贵，放最后。
- **2D 的隐藏成本在真空与 q 网格**（层状 + 大胞 + ZA 支需要足够密的 q 网格），单材料比同原子数的 3D 贵；先各 2 个把口径定死，再放量。
- **每个材料都要落 `resource_summary.json`**（核数/墙钟/GPU 标志）→ 这是图 4 与"规模—成本"实验的唯一数据来源，顺手就采，不要事后补。

---

## 7. 时间表与预期投稿日期

| 周 | 日期 | 任务 | 里程碑 |
|---|---|---|---|
| W1–W2 | 09/22–10/05 | **MoS₂ 口径决策 + S4 回 alm（`retry`）+ S6 重跑**（P0-2 端到端闭环）；Tier A 四材料跑通；图 1/图 2 定稿 | 2D/3D 各一条链全绿 |
| W3–W5 | 10/06–10/26 | Tier B 精度批（3D 8 声子 + 2D 3 材料）；噪声底/复现性统计（ZA 与 BHH 已就绪，直接引用） | 精度表成形 |
| W6–W7 | 10/27–11/09 | `zT(n,T)` 端到端（3D Bi₂Te₃/PbTe + 2D）；te-screen 漏斗实跑；bench 聚合 | 图 4、图 5 有数据 |
| W8–W9 | 11/10–11/23 | 三组消融（planner 替换 / 安全 / 规模—成本）；补充材料（provenance、session export、Zenodo DOI） | 证据链闭合 |
| W10–W12 | 11/24–12/14 | 写作（正文 + 图 3–5 + 表 1 + 补充材料）；组内审阅 | 首稿 |
| W13–W15 | 12/15–01/04 | 投稿；arXiv 预印本 + `git tag -a v1.3.0-paper` + Zenodo 存档 | 投出 |

**预期投稿日期：2026 年 12 月中旬**（保守：2027 年 1 月上旬）。
若 Tier B 不完整，退路是"图 4/图 5 放补充材料、正文只留 Tier A + 3D 8 材料声子表"。

---

## 8. 缺口、风险与硬阻塞

| # | 事项 | 严重度 | 处置 |
|---|---|---|---|
| 1 | ~~phonom-mlff 的 ZA 判定仍是错的~~ **已解除**（三副本合并到 `_common/za_2d.py`，2026-09-22 复验逐位吻合验收数字） | ✅ | 无；论文里可把"ZA 方向/本征矢量/微畸变 symprec"作为门禁设计的一项证据 |
| 1b | **MoS₂ 对称化链卡在 S5**：对称化后重拟合比未对称化更差（§46.9/§47） | 🔴 挡住 Tier A | **需你决策**：保留未对称化口径 / 接受对称化 / 换采样方法 |
| 1c | **同一虚频量三套阈值**（S5 0.10 / S5.1 −0.05 / lattice_kappa 0.5） | 🟠 | **需你给一个统一值**（§45 第 4 项）；这一条正好是论文"门禁设计"一章要讲的 |
| 1d | S4 当前是 findiff 数据集（214 个 disp 目录、**0 帧计算**），findiff 对比已取消（§49） | 🟠 | 要回到 alm 口径需 `retry`（只重新生成输入，不占 VASP）；属破坏性前置操作，需你点头 |
| 2 | S6 厚度归一化端到端：代码/逻辑已验证，但线上 MoS₂ 用的是**旧代码**产出的 `kappa_summary.json`（§39.2/§39.3） | 🟠 | 重跑一次 S6 即闭环（在 1b 决策之后） |
| 3 | ~~BHH 390/520 对照未做~~ **已做**（Δp≈1e-5，RASR=BHH 日志已确认） | ✅ | 无 |
| 4 | 图 2 上 pheasy-gpu 的加速比仍是 "×?" | 🟠 | 同一体系 CPU/GPU 各跑一次，测墙钟与数值一致性（Δ~1e-10 量级） |
| 5 | 图 4 / 图 5 无真实跑批数据 | 🟡 | 依赖 Tier B/C 与 te-screen 实跑 |
| 6 | 结构体检四判据仍未实现（AlN / Zn₅O₃ / BeO 曾被手工拦下） | 🟡 | 可留作"未来工作"，但需在论文中如实写明 |
| 7 | 2D 二档材料的相/厚度尚未逐一匹配文献 | 🟡 | Tier C 前先做匹配表 |

---

## 9. 已定 / 待定

**已定（user 2026-09-22）**

1. **两篇**：Paper 1 = AutoZT 平台（本文，leading first author）；**Paper 2 = 富勒烯那批（2926 材料）单独成篇**，见 §10。
2. **目标期刊 = *npj Computational Materials***（定为首选；备选只在桌拒/转投时启用）。

**待定**

3. **投稿窗口**：本学期截止 / 毕业送审时间点？（现按"Paper 1 于 2026 年 12 月中投稿"排）
4. **MoS₂ 口径决策**（§8 第 1b 项）——Tier A 唯一的拦路石。
5. **虚频阈值统一值**（§8 第 1c 项）。
6. 是否要我把 Tier A 四材料**按 `autozt` 正式流程起流水线**（生成输入 → 提交 → 推进）？

---

## 10. Paper 2：金属掺杂富勒烯热电（2926 材料那批）

| 项 | 内容 |
|---|---|
| 定位 | **材料/发现类**论文（不是方法论文）：一批金属掺杂富勒烯 / 富勒烯网络的热电性能筛选与候选 |
| 现成弹药 | 2926 个结构（38 金属 × 笼型 × 位点）的 MACE 优化 + 形成能已跑通；`fc-fit` 的 Mg₂C₆₀ 超胞力常数（OLS/LASSO/RFE 多组对照）；Mg₄C₆₀ / Mg₄C₆₀ 单层的声子链；`te-screen` 替代筛选器 |
| 缺什么 | ① 形成能 → 稳定性漏斗的统计图（决定谁进 κ / ZT）；② 入选材料的 kl + ke 全链（DFT 级）；③ 外部锚点（实验或独立文献值） |
| 与 Paper 1 的边界 | Paper 1 讲"平台怎么保证不出错"，Paper 2 讲"用平台发现了什么"。Paper 2 里 AutoZT 只作方法段一句话 + 补充材料引用，不抢主线 |
| 排期 | 本学期先做"数据整理 + 漏斗图"（与 Paper 1 的 te-screen 实跑**共用同一批作业**）；Paper 1 投出后全面启动 |
| 风险 / 退路 | 若漏斗筛不出足量有竞争力的候选，退成"方法 + 数据库"型论文（*npj Comput. Mater.* 接受这类） |
---

## 11. 简洁表格版（群回复 / 汇报用，可直接粘贴）

**Paper 1｜AutoZT：从单个 POSCAR 到 `zT(n,T)` 的热电计算平台**　目标期刊 *npj Computational Materials*　预期投稿 **2026-12 中旬**（保守 2027-01 上旬）
**Paper 2｜金属掺杂富勒烯热电（2926 结构）**　同刊，Paper 1 投出后启动

**表 1 创新点**

| # | 创新点 | 图 | 状态 |
|---|---|---|---|
| 1 | 端到端双链（电子：能带→散射→载流子输运；声子：稳定性→力常数→输运）汇合 `zT(n,T)`；**3D/2D 同一套工作流**（2D 原生：真空轴、面内弛豫、二维单位、ZA 支、厚度归一化） | 图 1a–b | ✅ 链路通 |
| 2 | **物理判据驱动的自治**：11 类"不报错但给错值"的静默错误设成门禁；稳态推进 **0 次 LLM 调用** | 图 2b | ✅ 多数已生产验证 |
| 3 | **LLM 只在失败路径介入**（MCP/CLI），动作分只读/推进/破坏性三档，破坏性需终端一次性令牌；planner 可替换 | 图 1c | ✅ 拦截 9/9 |
| 4 | 跨异构超算（真 SLURM CPU/GPU + 无 SLURM GPU 主机）+ per-file SHA-256 provenance + 可重放会话导出 | 图 1a | ✅ 17/17 |
| 5 | MLIP 预筛 + DFT 分级 + GPU 力常数拟合，把"精度—机时"做成可量化漏斗 | 图 2a | 🟠 加速比待测 |

**表 2 进度（硬数据，可复跑脚本）**

| 项 | 数字 |
|---|---|
| 技能 / 集群 | 20 / 3 |
| 批量实战 | 2926 结构（MACE 优化 + 形成能） |
| 安全 | 危险动作拦截 9/9 = **100%**，实际执行 **0** |
| 复现 | 跨集群 **17/17** 输入逐字节一致 |
| 技能覆盖 | 输入生成 16/18；首步真实提交 **15/15 OK** |
| 2D 验证 | 8 材料 S1 生产验证，面内应力 ≤ 0.114 kB |
| ZA / BHH | ZA 三副本合并，p=(1.99307, 1.98012)；390 vs 520 一致（Δp≈1e-5） |

**表 3 图分工**

| 图 | 内容 | 状态 |
|---|---|---|
| 图 1 | 系统总览（双链+验收门禁 / 3D-2D 统一 / LLM 经 MCP-CLI + 动作三档） | ✅ v6 |
| 图 2 | kl 技能（MLIP 预筛→DFT 力→力常数拟合→输运）+ 显式状态机 | ✅ v13；pheasy-gpu 加速比 "×?" **待填** |
| 图 3 | 安全消融（动作分档 / 令牌 / 审计） | 数据已有，待作图 |
| 图 4 | 资源与精度基准（method × backend × fit） | ❌ 缺实跑数据 |
| 图 5 | te-screen 筛选漏斗与排名 | ❌ 缺实跑数据 |

**表 4 2D/3D 验证计划**（基准集已建好，位于 `/mnt/d/tf_data/work_AutoZT/`，全部还没跑）

| 档 | 材料 | 目的 / 验收 | 期限 |
|---|---|---|---|
| **A**（必须先绿） | 3D **Si、BAs**；2D **MoS₂、WS₂** | 打通 2D/3D 端到端；WS₂ 直接对 pheasy 文献 **241.6 W/mK** | 2 周 |
| **B**（论文主表） | 3D 声子 8 全量（Si/Ge/diamond/GaAs/LiF/MgO/NaCl/BAs）；3D **Bi₂Te₃、PbTe**；2D **Bi₂Te₃_QL、SnSe、SnS** | κ 精度主证据 + `zT` 端到端 | 3–4 周 |
| **C**（可选） | 2D 二档 6 种（GaSe/GeS/GeSe/InSe/ReS₂/TiS₃）；JARVIS 1103 结构 | 泛化与外推；只跑便宜步骤做回归 | 投稿前能跑多少跑多少 |

> 基准集规模：2D 一档 10 个（含 graphene/hBN/phosphorene…，带 DOI 参考值）；3D 声子 8 个、3D 热电 6 个；参考值表 2D 24 条 / 3D 14 条。
