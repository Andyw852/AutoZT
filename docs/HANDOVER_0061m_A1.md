# 交接：[0061] ADP_DEGEN —— 方案 A1（DFT 形变势矩阵）进度与后续

日期：2026-10-21  分支：`claude/exciting-ride-1wamyy`  基线：`origin/master = 9785e75`（只有立方版 ADP_DEGEN）

---

## 续（2026-10-22，第三个会话）—— Si 真机跑出来的问题与修复（先看这一节）

| 提交 | 内容 |
|---|---|
| `317a1f7` | 运行期把 AMSET 的剪刀挪回 DFT（`undo_scissor`，报告记 `scissor_eV`）。Si 价带顶差 0.314 eV = 剪刀 0.628 的一半，被误判成带号对不上。 |
| `dd53cf8` | k0 核对改用简并子空间内求和（逐带的中心差分在简并点带 O(δ) 二阶污染）；新增 `tools/adp_k0_probe.py`。 |
| `bae7ebb` | `build` 用空间群操作补齐 AMSET 约化掉的全部图样（Si 原来矩阵只有 xx、yz，运行期判"缺图样"）。 |
| （本提交） | **AMSET 自己的对称补应变标签错**的改正 + probe 1b 标签核对 + `BAND_EDGES_REV` 改版。 |

**bae7ebb 之后 Si 的 probe**：6 个图样齐、k0 残差 0.7% 过；但对齐 xy 离散度 2.69 eV、xz 3.87 eV 超门槛 → 矩阵仍被拒。
**根因（已用 AMSET 0.5.1 源码在 EPM Si 上端到端重演）**：`get_symmetrized_strain_mapping` 把构型 E 的能带搬到 k' = M·k
（M = spglib Wᵀ），这组能带属于 Rᵀ·E·R，它却标成 R·E·Rᵀ。立方 Si 由 yz 经 C3[111] 补出的 xy、xz 被互换（xy+、xz± 错 0.22 eV，
xy− 碰巧走了对合操作是对的 —— 与实跑 xz 离散度大于 xy 一致）。所以错的是 S7.1 的 `deformation_signed.h5` 与 `deformation.h5`
的剪切分量，DFT 矩阵是对的。S7.1 的 Laue 群对称化**平均不掉**：EPM Si 走完整个 S7.1 流程，剪切分量 BZ 平均 |D| 价带顶
3.45 eV（12 个应变全算 7.23）、导带底 1.37（5.04）；对角分量不受影响。**这也影响 AMSET 基线**（ADP_DEGEN 关时）的形变势：
剪切偏小 -> ADP 散射偏弱 -> 迁移率偏高。受影响的晶系：立方、四方（C4 换掉 xz/yz 的正负号）；六方 / 三方 / 正交及更低不受影响。

**改正**：`dp_symmetrize.symmetrized_strain_mapping`（照搬 AMSET 的循环，只把标签换成 strain.rotate(Rᵀ)），S7.1 的
`amset deform read` 包装（`_NBANDS_FIX_CODE` 里换 `amset.deformation.potentials.get_symmetrized_strain_mapping`）与 2D 真空参考
路径（`_vac_code`）都用它；与"12 个应变全算"逐位一致。测试 `tests/test_0061m_strain_label.py`（桩 AMSET + EPM：改正版 8 个补出的
应变全对、换回原版标签 xy/xz 错、S7.1 包装在桩 AMSET 里真装上补丁、真空路径不再调原版）；`test_adp_k0_probe.py` 加了模拟 AMSET
标签错的一条（1b 节认出 S_xz 实为 X_xy、S_xy 实为 ½(X_xy+X_xz)）。

**本地会话接下来（动集群 / 改状态都按 AGENTS.md 走审批卡片）**：
1. **先确认假设，不改任何东西**：把新版 `adp_k0_probe.py`（连同 `adp_degen.py`、`deform_wave.py`）拷到现有 S8 运行目录跑
   `python adp_k0_probe.py --dir .` —— 1b 节应出现 `S_xz：… 最好 ±X_xy … ★ 标签错` 和 `S_xy：… 最好 ½(…) … ★ 标签错`。
2. 用户同意后合并到 master（本分支 fast-forward）。
3. 重跑 S7.1（只是重读，几分钟；S7 不用重算）。日志里 grep `[strain_label_fix]`，Si 应是"对称补了 8 个应变，其中 M 个 AMSET 原版
   会标成别的应变"（M > 0）。
4. 在新 S7.1 产物上再跑 probe：1 节 6 个图样离散度都在门槛内，1b 全部"最好 = 自己"。
5. 然后才跑 S8。**提醒用户**：AMSET 基线（ADP_DEGEN 关）结果也会变；其它立方 / 四方材料的 S7.1 也要重读，哪些重跑 S8 由用户定
   （S8.2 会对旧 band_edges.json 告警，因为 `BAND_EDGES_REV` 改了）。

---

## 续（2026-10-21 第二个会话）

**A1 的代码已全部写完并通过本地测试**（第 1–5 节是上一次的交接，保留作背景；第 4 节的计划已逐条落实，偏差见下）。

| 提交 | 内容 |
|---|---|
| `9b7acd0` | 修窗口选取（3.3 的失败）：extract 只按能量选、MARGIN 2→6；build 以未形变构型为准、在各形变构型两侧仍剩 PAD=2 条的范围内扩简并组，不再连锁缩。6³ 离格点误差 0.16→0.024。7/7 过。 |
| `53c580c` | 4.2–4.4 接线 + 运行期核心区：S7 `DEFORM_WAVE*`、S7.1 build、`ibands` 进带符号 h5、S8/S8.4 `ADP_DEGEN_MATRIX` 与软链；`adp_degen` 的 `matrix_grid`（对齐）、k0 处 Z = Bᴴ D̂ B（覆盖度认倒格矢表示、对称化、与 S 核对，失败退回拟合）、rcut 内预算 Zk（不再减 T）。 |
| `570eba3` | 4.5 网格公式：钩子多记 calculate_rate 参数 / masks / 截面插值参数；`grid_pair_values`、`matrix_apply`（核心区 3 个 DFT 网格间距内混 k·p）；共享内存进子进程；两个插件改走 `matrix_apply`。 |
| （本提交） | 文档：METHODOLOGY [0061-m] 小节、CHANGELOG、step8 step.conf 的 `ADP_DEGEN_MATRIX`、本交接。 |

与第 4 节计划的偏差：
- `DEFORM_WAVE_MARGIN` 出厂 6（不是 2）；没有 `step7_deform/step.conf` 文件，三个键的说明写在 `gen_step9_deform.py` 里（与 IONRELAX_* 同）。
- 运行期没有 `ADP_DEGEN_SCOPE`：网格公式只处理"有简并谷的那种载流子"（计划里的默认），没做 `all`。
- Zk 不在子进程里现算，而是主进程对 rcut 内 AMSET 网格点预算好（子进程只拿小数组）；网格公式本身要 StrainGrid，走共享内存。
- 网格公式只取四面体顶点（不建 (带, 全部 k, 28) 的数组 —— 细网格上每次几百 MB）。

测试：新增 `tests/test_0061m_matrix_mode.py`（7）、`test_0061m_wiring.py`（15）、`test_0061m_grid_formula.py`（8），
`test_0061m_deform_wave.py` 7/7，原 `test_0061_adp_degen.py` 25/25（改了一条源码断言：2D 插件现在调 `matrix_apply`）；
`skill/ke-dft-cpu/test_stale_input_run.py` 改了一条源码断言（`_mismatch_wo_lwave`）。
`skill/ke-dft-cpu/test_downstream_wiring.py::test_every_upstream_gen_invalidates` 在 master 与 a9d90e4 上同样失败，与本次无关。
本地没有 AMSET：钩子对照的是下载的 0.4.19 / 0.5.1 源码（`calculate_rate`、`get_masks`、`get_cross_section_values`、`get_g_maps`、
`create_shared_array`），用 0.5.1 形状的假 `calculate_rate` 测过；**真 AMSET 跑一次之前，网格公式只算"源码对照 + 合成测试"过关**。

还没做 / 要用户定的：
1. **真数据验证**：Si 要 retry S7_deform（undeformed 重算一个单点；12 个 ionrelax 只补静态段）→ S7.1 → S8。动集群作业要按 AGENTS.md
   走审批卡片；另一个会话在跑 Si 的 S7 ionrelax（4298484–6），等它完再动。看 `adp_degen_fit.json` 的 `matrix`（对齐离散度、带号）、
   各谷 `strain.source` 是不是 `deformation_matrix.h5`、`coverage_min`、`k0_rel_res`，以及 `method` 里有没有"网格公式"。
2. 网格公式在散射子进程里的计数（`_PAIR["n_grid"]` / `n_grid_lowcov`）回不到主进程报告；现在每个进程第一次生效打一行
   `[adp_degen] 网格公式生效（进程 …）`，自检不过打 `[WARN] … 本进程关掉网格公式` —— 真跑时在作业输出里 grep 这两句。
3. pymatgen 对照 WAVECAR 读法的脚本仍没进仓库（真实 VASP 样例只有 1 个 k 点的）。
4. 合并 master：仍须用户明确同意。
5. **S7 网格加密的验证工具**（用户 2026-10-21 要求）：`skill/ke-dft-cpu/tools/nscf_dp_check.py`。用户的本机会话在能找到 Si 的 S7（gen 过即可）的集群上 `make` -> 在每个构型目录 `sbatch submit.sh`（独立诊断作业，用户已同意）-> `analyze`。看三项：非自洽 OUTCAR 的芯能级参考、带边 D（含 Si 电子谷的静水值 Ξ_d + Ξ_u/3，xx 图样在谷星上的平均）、耗时。通过了再谈改 S7（自洽段多写 CHGCAR 不触发整批重算要照 LWAVE 的办法处理）。

## 0. 一句话

简并 / 近简并带边的 ADP 修正，已经做完"通用 k·p 版 + 在 AMSET 四面体积分里逐对替换矩阵元"，并已推送到分支。
用户要求更精确，选了**方案 A1**：形变势矩阵直接从 DFT 波函数投影得到，不再靠本征值反推。
A1 只做了第一块（`deform_wave.py`：WAVECAR 读取、抽取、投影、对齐、网格算符），**S7/S7.1/S8 的接线和运行期都还没写**。
**Si 的 S8 要等 A1 做完再跑；分支没有合进 master，合并必须等用户明确同意。**

---

## 1. 用户的决定与约束

- **2026-10-21 用户的选择**
  - 方案 A（用 DFT 波函数算形变势矩阵）。
  - 范围只做 **A1**：不加"在 AMSET 密网格上补算非自洽"的新步骤（A2）。
  - Si 的 S8 **等 A1 做完再跑**；分支先不合 master。
- **A1 的含义**（已向用户说明并获同意）
  - 形变势矩阵全部来自 DFT 波函数，所有 k 点都用。S7 写 WAVECAR，作业里就地抽成小文件后删掉 WAVECAR，S7.1 做投影。
  - 离简并点稍远的态直接用 DFT 态（AMSET 自己的系数）。
  - 只在最内层（简并点附近一两个网格单元，那里 AMSET 的系数线性插值不可信）用 DFT 拟合的 k·p 态，报告里给逐点误差。
- **不能违反的规则**（AGENTS.md）
  - 铁律 2：破坏性操作要走审批卡片，且必须用户亲口同意。
  - 铁律 10：新建步骤 / 新 gen 脚本要先请示；`deform_wave.py` 属于已获批 A1 的辅助模块。
  - 临时文件放 `tmp/` 或 scratchpad。
  - 用中文回复。
  - 不碰集群作业（另一个会话在跑 Si 的 S7 ionrelax：4298484–6）。
- **Git**
  - 只在本分支开发。
  - 不开 PR。
  - commit 结尾加两行：`Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` 和 `Claude-Session: https://claude.ai/code/session_01KwgbCRMcAbxcek9mJ68naT`。
  - 产物里不写模型 ID。

---

## 2. 已完成并推送的（分支上，master 上没有）

| 提交 | 内容 |
|---|---|
| `55150d9` | [0061-g] 通用引擎：任意晶系、谷星、价带导带、SOC，2D（S8.4）接入。<br>谷识别、Löwdin 投影 k·p、对称性从数据解出、k·p 拟合（加密精修）、应变算符从 S7 本征值拟合，出比值表。 |
| `0a9c6b6` | 比值表照 AMSET 弹性散射的真实算法：动量弛豫因子 + 固定分数坐标形变势的 \|·\| + m_e v⊗v。 |
| `e4866ed` | [0061-p] 在 AMSET 四面体积分里逐对替换矩阵元，默认 `ADP_DEGEN_PAIR = auto`；源码结构对不上时退回比值。<br>同时修了 bohr 单位 bug（比值表里的 v⊗v 原来按 Å 算）。测试 25 条全过。 |

[0061-p] 的做法：
- 钩子包住 AMSET `calculate_rate` 的三个中间量：`tbs.get_masks` 给末态带号，`get_cross_section_values` 给截面重叠，`get_fine_mesh_qpoints` 给 q 点 → 四面体映射。
- ADP 的 `factor()` 对子空间内的 (初态, 末态) 返回 `F_full / AMSET 重叠`；AMSET 随后乘的正是这个重叠，两者精确抵消。
- 源码锚点已对照真实 AMSET 0.4.19 与 0.5.1，五处全部吻合。
- 相关代码：`skill/ke-dft-cpu/adp_degen.py` 的 `install_pair_hooks` / `pair_factor` / `pair_apply`，`step8_amset/amset_dpo3d.py`，`step8.4_amset2d/amset2d_plugin.py`。

文献（已核实）：
- Bir & Pikus 1974
- Khan & Allen, PRB 29, 3341 (1984)
- Hamann, Wu, Rabe, Vanderbilt, PRB 71, 035117 (2005)：约化坐标下的应变 DFPT
- Li, Graziosi, Neophytou, PRB 104, 195201 (2021)：Si 价带 3×3 声学耦合矩阵，HH↔LH 很大
- Kaasbjerg, Thygesen, Jauho, PRB 87, 235312 (2013)
- Poncé et al. PRB 97, 121201 (2018)；Brunin et al. PRL 2020；Park et al. PRB 2020（四极矩）
- Fang et al. PRB 98, 075106 (2018)：应变 Wannier 哈密顿量

注意："把 DFT 形变势矩阵接进 AMSET 四面体积分"这个组合**没有文献先例**，是本项目自己的实现，向用户说明时要讲清楚。

---

## 3. A1 已写但未提交的文件（本次一并提交为 WIP）

### 3.1 `skill/ke-dft-cpu/deform_wave.py`（新模块，只依赖 numpy + h5py，不依赖 AMSET）

- **`Wavecar`**：流式读 WAVECAR，按 (自旋, k, 带) 定位记录。
  - 支持标准版、ISPIN=2、非共线；gamma-only 会明确报错。
  - **已与 pymatgen v2022.0.8 的读法逐位对照**，用的是 pymatgen 仓库 `test_files/` 里的真实 VASP 文件：WAVECAR.N2、.N2.45210、.N2.spin、.H2_low_symm、.H2.ncl。
    - 下载地址：`https://raw.githubusercontent.com/materialsproject/pymatgen/v2022.0.8/test_files/<名字>`。
    - 对照脚本当时放在 scratchpad，没进仓库，建议补成测试夹具。
  - 局限：这些样例都只有 1 个 k 点；多 k 点时的记录偏移只由 EPM 写读往返测试覆盖。
- **`extract()` + 命令行 `extract`**：在 S7 作业里、VASP 跑完后执行，生成 `deform_wave.h5`。
  - 只抽带边 ± `EWIN_K`（0.8 eV）以内有态的 k 点。
  - 能带窗口 ± `EWIN_B`（1.2 eV），两侧再各加 `MARGIN`（2）条；窗口边界在 `DEG_TOL`（0.05 eV）以内不切开简并组。
  - 系数截到 `ENCUT_STORE`（300 eV），按 Miller 并集存。
  - `--delete` 表示抽完删 WAVECAR。
- **`build()` + 命令行 `build`**：在 S7.1 里执行，生成 `deformation_matrix.h5`。
  - 由晶格算 Green–Lagrange 应变，判断图样（与 AMSET 同口径）。
  - 每个 k 点：两组系数各自 Löwdin 正交化，然后 `H^s = S^-½ O diag(E^s) O† S^-½`（`O = ⟨未形变|形变⟩`），`X^(P) = (H^+ − H^−)/(δ+ − δ−)`。
  - 存的是原生约定（剪切 ε_ab = ε_ba = δ）、**未对齐**的 X，外加未形变的正交化系数和诊断量（覆盖度 `q_min`、赝波函数非正交度）。
  - 形变构型抽到的带若盖不住未形变窗口，会缩窗口，同样不切开简并组。
- **`alignment()`**：在运行期，用 S8/S8.4 实际读的那份带符号 h5，求每个图样的参考能级常数 `c_P = X_nn + S_n`。
  - 只用非简并态；同时返回离散度，带号错位时离散度会立刻变大。
- **`StrainGrid`**：在 S7 网格上对**算符本身**做三线性插值：`D̂(k) = Σ α_i Σ_ab |a,k_i⟩X_ab⟨b,k_i|`，与规范无关。
  - 缺的格点用时间反演补：`Θ|a⟩`，`X → X*`；SOC 时自旋分量做 (c↑, c↓) → (c↓*, −c↑*)。
  - 跨布里渊区边界时按 Miller 平移。
  - `apply(k, vecs)` 把 D̂ 作用到 AMSET Miller 列表上的系数，单位 eV。

### 3.2 `tests/_epm_si.py`（测试台）

- 金刚石结构的经验赝势模型：Cohen–Bergstresser 的三个壳层换成光滑的紧支撑鼓包，体积按 Ω0/Ω 缩放。
- `dH()` 在固定 Miller 基下直接对平面波矩阵求导，作为精确参照。
- `write_wavecar()` 按 VASP 格式写文件。
- 不用慢衰减形状因子的原因：它让平面波截断随应变变化（Pulay 效应），误差可达 25%。
- 改成紧支撑后，ENCUT 150 / 200 / 300 eV 时投影误差分别是 0.63% / 0.27% / 0.23%。

### 3.3 `tests/test_0061m_deform_wave.py`：7 条，**6 过 1 败**

| 测试 | 状态 | 说明 |
|---|---|---|
| WAVECAR 写读往返、带边 | 过 | |
| 应变图样识别、抽取 | 过 | |
| 投影 vs 精确 X | 过 | < 1%，含 Γ 三重态非对角元和符号 |
| 对齐读回参考常数 | 过 | 离散度小；带号错位时离散度变大 |
| 格点上精确（含时间反演点） | 过 | |
| 命令行 extract / --delete | 过 | |
| **离格点收敛（3³ → 6³）** | **败** | 原因还没查，见下 |

离格点收敛测试失败的线索：
- 最后一次改动是"形变构型抽到的带盖不住时缩窗口"。
- 日志显示窗口从 1-7 缩到 1-3、从 0-11 缩到 0-3（形变构型只抽到 0-6 / 0-7）。也就是**导带被缩掉了**。
- 根本问题：各构型**各自按能量独立选窗口**，`DEG_TOL` 的外扩在未形变（简并）和形变（简并被劈开）两边结果不一样。
- **建议修法**：窗口按"占据带数 ± 固定带数"选，或 S7 各构型统一用 undeformed 选出的带号范围，并把 `MARGIN` 放大。不要再靠缩窗口。
- 修完后，该测试的误差阈值也要按新窗口重新确认。

---

## 4. 剩下要做的（按顺序）

### 4.1 修 `deform_wave` 的窗口选取，让 7 条测试全过（见 3.3）

### 4.2 S7 接线（`step7_deform/gen_step9_deform.py`）—— 方案已定，没写

1. **新 step.conf 键**：`DEFORM_WAVE = on|off`（默认 on）、`DEFORM_WAVE_KEEP = false`、`DEFORM_WAVE_MARGIN = 2`。
2. **不要改 INCAR 模板的 LWAVE**。
   - `kc.run_incar_mismatch` 把 LWAVE 当输出类开关：旧结果 `.FALSE.`、新 INCAR `.TRUE.` 会被判"输入变了"，所有单点和离子弛豫会**整批重算**。
3. **按目录类型分别处理 LWAVE**：
   - **离子弛豫目录**：只在静态段 `INCAR.static` 里加 `LWAVE = .TRUE.`（`_ionrelax_incars` 加开关）。
     - 主循环"已完成"判定里，从 `_why` 中滤掉 `LWAVE:` 这一条。
     - 已完成但缺 `ionrelax/deform_wave.h5` 的目录，走现成的 `prepare_static_only`（只补静态段，不重做弛豫，ISTART=0/ICHARG=2）。给它加一个"原因"参数，打印别再写成"缺 LOCPOT"。
   - **没有离子弛豫的目录**：undeformed 总是要；3D 且 `IONRELAX_3D = false` 时各形变目录也要。在顶层 INCAR（stepconf 覆盖之后）加 `LWAVE = .TRUE.`。
     - V168 会自动归档旧结果并重算，undeformed 必然重算一次。
   - 2D 的面外 clamped 目录不需要波函数（2D 核只用 xx/yy/xy）。
   - 要在主循环之前就知道哪些目录做离子弛豫：把 `_ionrelax_params` 和 `select_ionrelax_dirs` 提前调用（后者无副作用，只打印一行）。
4. **submit.sh 末尾追加抽取段**（最后统一做一遍）：
   - 内容：`( <AMSET_ENV_SRC> && "${CONDA_PREFIX:?}/bin/python" <相对根>deform_wave.py extract --margin N [--delete] ) || echo "[deform_wave][WARN] …" >&2`。
   - 相对根：顶层目录用 `../../`；离子弛豫目录（脚本最后停在 `ionrelax/`）用 `../../../`。
   - 抽取失败不能让作业判失败。
5. `skill.yaml` 里 S7、S7.1、S8、S8.4 的 `gen_need` 都加上 `deform_wave.py`。

### 4.3 S7.1 接线（`step7b_read/gen_step9b_deform_read.py`）

- 在 `amset deform read` 成功之后、h5 挪进本步目录之前执行：`cd dfm && python <cwd>/deform_wave.py build --root . --out <out>/deformation_matrix.h5 <有 deform_wave.h5 的目录>`。
  - 目录用 `_rel` 清单（离子弛豫优先），在 AMSET 环境里跑。
- 先删旧的 `deformation_matrix.h5`。
- 缺文件或 build 失败**只告警、不让 S7.1 失败**（S8 会退回 k·p 拟合的应变算符）。
- 顺手：`dp_symmetrize.write_signed` 把 `ibands_<spin>` 也写进带符号 h5（现在没存），这样运行期对齐的带号映射有直接依据。目前靠运行期的带号偏移 `off` 推断。

### 4.4 S8 / S8.4 gen

- `gen_step10_amset.py` 的 `resolve_adp_degen_3d` 和 `gen_step14_amset2d.py` 的 `resolve_adp_degen_2d`：如果 `step7b_deform_read/deformation_matrix.h5` 存在，就软链进运行目录，记录 `rec["matrix_h5"]`。
- 可加 `ADP_DEGEN_MATRIX = auto|off` 开关。
- 改 `_SKILL_REV`。

### 4.5 运行期矩阵模式（`adp_degen.py` + 两个插件）—— 设计已定，没写

- **加载与准备**：
  - `runtime_tables` 里读 `matrix_h5`。
  - 用 `deform_wave.alignment(data, 带符号 S, kS, band_offset=off, delta)` 求对齐常数。`off` = DFT 带号 − AMSET 带号，`runtime_tables` 里已有，带符号 h5 的带轴与 AMSET 相同。
  - 离散度过大就拒绝矩阵模式、退回现状，并在报告里写明。
  - 建 `StrainGrid(data, overlap_calculator.gpoints, offsets)`。
  - 需要的图样：3D 六个，2D 是 xx/yy/xy。缺了就退回。
- **钩子扩展**：
  - 包住 `C.calculate_rate`，记下 `(tbs, overlap_calculator, coeffs, coeffs_mapping, spin, b_idx, k_idx)`。worker 是按名字调用的，替换模块属性就能生效；插件的 `worker` 包装保证子进程会导入插件。
  - `get_masks` 记下返回值：`property_mask`、`band_kpoint_mask`、`band_mask`、`kpoint_mask`。
  - `get_cross_section_values` 第一次调用时记下 `tet_contributions`（参数），后面用**原函数**按同一组权重插值。
- **网格公式**（对每个初态）：
  1. `c_n = coeffs[spin][coeffs_mapping[b,k]]`；`w^P = grid.apply(k, c_n)`，共 6（2D 为 3）个向量。
  2. 末态 `F = coeffs[spin][coeffs_mapping[band_mask, kpoint_mask]]`。0.5.1 要按 `g_diff = (k′−k) − pbc_diff(k′,k)` 用 `overlap_calculator.g_maps[−g+1]` 平移，与 AMSET 0.5.1 的 `_get_overlap` 逐字相同；0.4.19 不平移。
  3. 振幅：`A = conj(F) @ [w_1..w_P, c_n]ᵀ`，形状 (n_final, P+1)。最后一列的 \|·\|² 必须等于 AMSET 的顶点重叠，作为自检。
  4. 双线性量：`G = A Aᴴ`，转成实坐标（7×7 共 49 个）。
  5. 顶点值：`all_G[band_kpoint_mask] = G`，取 `vert = all_G[property_mask]`，用原 `get_cross_section_values(vert, *tet_contributions)` 插到截面上（它支持带尾部维度的实数数组）。
  6. 对每个细 q 点：`e_λ = [d_scale·Es_P/HARTREE (P), Es:(v⊗v)]`，`F = Σ_λ e_λᵀ G e_λ / C_λ`。
     - Es 用原生系数：正应变取 Es_aa，剪切取 Es_ab，**不乘 2**。
     - v 是 AMSET 的初态速度，原子单位。
     - 运动学项已包含在 X 里（固定分数坐标），**不要再减 T**。
  7. 返回 `keep × F / tet_overlap[mapping]`，`keep` 是插件的谷内掩码；重叠 < 1e-14 的点保持原样。
- **核心区（k·p 态 + DFT 算符）**：
  - 谷心附近 `r_core = min(rcut, CORE_CELLS × S3 网格间距)`，`CORE_CELLS` 可设为 3。
  - 初态、末态用现有 k·p 模型（代表谷坐标）。
  - 应变算符用 `Zk^P = Bᴴ D̂^P(k_r) B`：B 是模型基（AMSET 在 k0 处的系数 × 分块酉阵 W），`k_r` 是映回代表谷的 k。
  - 在 `pair_factor` 里用 `Σ_P Er_P Zk^P` 代替 `strain_H(Z, Er)`，**去掉 T 项**。
  - 在外 20% 与网格公式线性混合。
  - 核心区里不属于子空间的末态，用网格公式。
- **报告**：`ratio_mean` 改用 `Z0 = Bᴴ D̂(k0) B`，不再拟合，所以也没有歧义解；`_fit_strain` 只作为没有 matrix_h5 时的退路。
- **范围**：默认只处理"有简并谷的载流子类型"的全部带；`ADP_DEGEN_SCOPE = all` 时两种载流子都处理。
- **共享内存**：AMSET 用 spawn 起子进程，大数组要用 `amset.util.create_shared_array`。complex64 要先 view 成 float32 存。
- **单位提醒**：`amset_data.structure` 是 bohr；X 是 eV，要除以 `HARTREE_EV`；C 是原子单位。

### 4.6 测试

- 用 EPM 构造与 AMSET 同形的假 `calculate_rate`，比较网格公式的逐对振幅和精确的 `⟨m,k′|∂H(k)/∂ε|n,k⟩`（初态侧）。要量化窗口截断和插值误差。
- 测钩子、0.5.1 的 g_diff 平移、自检（最后一列 \|·\|² = AMSET 重叠）。
- 测 S7 gen：LWAVE 只出现在该出现的地方；已完成的离子弛豫只补静态段；抽取段写进 submit.sh。
- 测 S7.1 调 build，S8/S8.4 记录 `matrix_h5`。
- 现有 `tests/test_0061_adp_degen.py`（25 条）必须保持全过。

### 4.7 文档

- `skill/ke-dft-cpu/METHODOLOGY.md` 的 [0061] 一节：加 A1 的做法、近似和文献。注意说明"组合无先例"。
- `CHANGELOG.md` 加 [0061-m] 条目。
- `step7_deform/step.conf` 写 `DEFORM_WAVE*`，`step8_amset/step.conf` 写 `ADP_DEGEN_MATRIX` 等。

### 4.8 收尾

- 跑全量测试：`python -m pytest -q -p no:warnings tests skill/ke-dft-cpu/step8_amset`。
  - 环境要有 numpy、h5py、pytest。本会话用的是 scratchpad 里的 conda 环境，其他会话访问不到。
  - 已知只有 `tests/test_suites.py::test_local_suite[agentgate]` 失败，是环境问题，master 上也一样失败。
  - 全量大约要 15 分钟以上。
- 提交并推送到本分支。
- 向用户汇报，并再问一次合并 master 的事。

---

## 5. 已知的近似与风险（汇报时要讲）

- **ADP 框架本身**：长波应变、弹性散射、均分，与 AMSET 原版相同。只有 EPW 那条路（方案 B，用户没选）能去掉。
- **PAW 赝波函数**：系数不含补偿项，与 AMSET 的波函数重叠同口径；`build` 会报告赝波函数的非正交度。
- **中心差分误差**：O(δ²)，δ = 0.5%。真实 VASP 里截断球随应变变化（Pulay 效应）也有影响；S7 用 ENCUT = 1.5 × ENMAX，效应较小。
- **窗口截断**：插入的完备集只到窗口里的带。远带的贡献是 O(q·v/ΔE)，需要用 EPM 量化。
- **AMSET 自己的系数插值**：线性插值、相位随机。在简并点附近不可信，所以最内层保留 k·p。用户已认可这个取舍。
- **S7 的网格**：3D 默认用 vaspkit 生成的粗网格（Si 为 11³），只有设了 `DK_MAX_3D` 才与 S3 相同。A1 的 `StrainGrid` 不要求 S7 和 S3 共网格，但网格越粗，插值误差越大。
- **磁盘**：WAVECAR 只是过渡文件，作业里抽完就删；`deform_wave.h5` 只含带边附近的 k 点。Si 估计每个目录几十到一两百 MB。
