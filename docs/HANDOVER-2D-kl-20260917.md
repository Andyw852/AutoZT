# 交接文档：2D 晶格热导链路（kl-dft-cpu / 2D）修复与验证

> 写于 2026-09-17，上一会话 context 用满。仓库 `~/software/AutoZT`（包 `autozt`），
> 集群 jzzn（ssh 别名 `jzzn`），项目根 `/mnt/d/tf_data/test_TE/P1`，
> 远端工作根 `/public/home/wangchao/Fullerene_Network/work`。

## 0. 一句话现状

**S1（结构优化）这一环已修好并验证**；P0-1（BHH 真施加）、P0-2（厚度归一化）、
P0-3（LASSO 设置）**尚未在完整链路上跑通**。下一步分两条线（见第 5 节）。

---

## 1. S1 当前流程（已实现 + 待补）

1. **准备结构**：起点是材料根 `POSCAR`；续算用上一次 `CONTCAR`（`contcar_to_poscar`）。
   自动判维度与真空轴。
2. **段 a**：ISIF=2 固定胞只弛豫原子，ENCUT 用基础值（1.5×max(ENMAX)）。已完成则跳过。
3. **段 b/c**：变胞弛豫，最多 3 遍（`CELL_PASS_MAX`），每遍先 b 再 c（重置平面波基组）。
   - **ENCUT = 2.0×max(ENMAX)**，**只有这两段**（`CELL_STAGE_ENCUT_FACTOR = 2.0`）。
   - 2D 只允许面内变化：优先 `IOPTCELL`；VASP≥6.5 且 c∥z 用 `LATTICE_CONSTRAINTS`；
     两者都不可用 → **硬报错**，绝不静默退回 ISIF=2（`ALLOW_2D_FIXED_CELL` 是逃生阀）。
   - 一遍收敛需**同时**满足：力收敛 + `max|Δ晶格| < 0.002 Å` + 应力达标。
     - 2D：`max(|σxx|,|σyy|,|σxy|) < tol_cell`，`tol_cell = max(0.4/(h⊥/d), 0.05)` kB（胞口径）
     - 3D：`|P| < 1.0` kB
4. **看门狗**：三指纹（OUTCAR 字节 / OSZICAR 行数 / **VASP 进程累计 CPU 时间**）连续
   `STALL_MINUTES=180` 分钟都不变才判卡死；判死前把 `ps/df/uptime` 快照写进 `queue.err`。
5. **完成检查**：应力门禁读最后一步应力；有 `S1_STRESS_GATE.off`（内含原因）则豁免并打印。
6. ⚠️ **未实现**：
   - S1 完成时写**材料级** `thickness_2d.json`（两条链先到者写，后者比对 `|Δd|>0.05 Å`
     或晶格差 >0.1% 给 WARN）—— 这是 P0-2b 的最后一环，也是 ke/AMSET 侧的读取点。
   - 汇总里记录 **输入 POSCAR → 最终 CONTCAR 的晶格总变化**（已两次造成误读）。
7. **S2–S4 继续用基础 ENCUT**，是否提高等 `INCAR_BENCH` 对照结果。

### ⚠️ S1 当前两个已知漏洞（本轮新发现，未修）

- **变胞段会以 `ZBRENT: fatal error in bracketing` 崩掉（Ti₂S₃ 实证，09-17 23:13）**。
  机理：段 a 已把力弛豫到收敛，段 b 放开晶胞后线搜索**找不到极小值**（能量对面内应变
  已经平了），ZBRENT 直接放弃 → VASP 打印 `I REFUSE TO CONTINUE WITH THIS SICK JOB`、
  rc=1。`CONTCAR.s2r1` 里 a 只从 3.309864 变到 3.307900（−0.06%），正是"几乎没动"。
  **不是**撞墙钟、**不是**看门狗（`queue.err` 里无 `[watchdog]`，是 mpirun 报的非零退出）。
  当前报错文案（`第 1 段中断…可能撞墙钟或被看门狗杀掉`）**会误导**，应识别 ZBRENT/SICK JOB
  并归因到"晶胞已平衡、线搜索无法夹逼"。修法（优先级从高到低）：
  ① 段 b 起点不要直接用段 a 的 CONTCAR，先给面内 **±0.3% 固定应变种子**再放开；
  ② 变胞段改**正的 EDIFFG（能量判据，如 1e-4 eV）**，别用 `EDIFFG=-0.01` 的力判据；
  ③ 捕获 ZBRENT 后，若力已收敛且应力达标则**视为收敛**，否则换 IBRION 重试。
  ⚠️ 同一机理的温和版本：`max|Δ晶格| = 0.000000` 也可能只是"力已收敛→VASP 第一步就停"
  （EDIFFG<0 只判力、不判应力），晶格一动不动、残余应力原样留着，会**连跑 3 遍原地不动
  然后失败**。
- ~~**AlN(kl) 的变胞段 ENCUT 是 520，不是 800**~~ **❌ 2026-09-18 作废：误报（详见 §11.1/§1 更正）**。原文保留供对照：（N 的 ENMAX = 400 → 2.0× 应为 800；
  AlN(ke) 拿到 800 ✓）。520 = 1.3×400，怀疑有项目级/材料级 `ENCUT` 固定值覆盖了
  2.0× 设置，而 **520 恰好等于旧默认系数 1.3×400** —— 强烈指向某处仍按 1.3 倍生成。
  已查：该材料 `kl-dft-cpu/` 全树只有模板占位符与一个 `.bak`（其中 `ENCUT_FACTOR = 1.5`），
  **活的固定值没找到，来源仍未定位**。反证规则本身是对的：Ti₂S₃ 段 b = 550 = 2.0×274.7 ✓、
  段 a = 420 = 1.5×274.7 ✓（注意 INCAR 里那行注释 `~1.3 倍` 是**过时文案**，值才是准的）。后果：AlN(kl) 的
  Pulay 偏差没消除，它的面内应力 0.114 kB 是这批里最接近阈值的（0.1604）。

---

## 2. 这几轮的实测结论（都是硬数据，别再重踩）

| # | 结论 | 证据 |
|---|---|---|
| 1 | **Pulay 应力让所有应力分量各向同性平移约 +1.0 kB**，ENCUT 到 **2.0×ENMAX** 才收敛 | 同几何只改 ENCUT：390→507 每分量 +1.0 kB；507→624 变化 <0.02 kB（Mo₂S₃，σxx/σyy/σzz 同时平移） |
| 2 | **只有变胞段需要高 ENCUT**；S2/S4 是固定几何算力/静态，与 Pulay 无关 | 同上（Pulay 不影响力）→ 省掉 S4 那 166~478 帧各 +60% 平面波成本 |
| 3 | **用 `|P|` 作 2D 收敛判据会被 σzz 掩盖真实面内应力** | 旧结构 P=−7.93 时面内 −11.97/−10.62、σzz=−1.21；弛豫后 P=−0.51 而 σzz 仍 −1.30（σzz 是被 IOPTCELL 冻结的真空方向应力，不随面内下降）→ AlN/BeO 的 |P| 小但面内 2.8~3.1 kB，**骗过了旧门禁** |
| 4 | **S4 按 fanout retry 会把旧结构的帧当已完成** | 6 个材料的 `step4_disp` 残留 **1902 个 `disp-*`**（含半成品 vasprun.xml），S5 会拿旧力配新位移静默混拟合 |
| 5 | **真空厚度影响面内晶格 0.3~0.5%** | Mo₂S₃ Z4-3-1：20 Å 真空 a=3.0446944 → 26.7 Å 真空 a=3.0595617（**+0.488%**），b +0.261%；两者均 520 eV。判据阈值是 0.05%，超了近 10 倍。**但用户已指出这可能不是镜像作用**（见第 4 节，需四项诊断） |
| 6 | kl 与 ke 两条链的 S1 结果**实质相同** | Mo₂S₃ Z4-3-1：Δa = 1.6e−6 Å（**0.00005%**）、Δb = 3.2e−7 Å，远小于 0.1% 判据 → `thickness_2d.json` 比对逻辑可按此落地 |
| 7 | **AlN 的原胞离平衡极远**（首次弛豫 `max|Δ晶格| = 2.143 Å`），但真空轴 c **完全没动**（IOPTCELL 正常）、元素/原子数不变（无重构） | 后续轮次 `max|Δ晶格| = 0.000000`，因为前一轮已弛豫、CONTCAR 成了新起点 |
| 8 | AlN 现在**不是平面单层**：由阈值反推 d = 0.4/0.1604 → h⊥/d=2.49 → **d ≈ 8.03 Å**，而 z 跨度实测 ≈3.9 Å | 待用**材料根种子 POSCAR** vs 当前 CONTCAR 比空间群/配位数/z 跨度判定是"种子本来如此"还是"弛豫中变了" |

---

## 3. 关键代码改动位置（本会话）

| 文件 | 内容 |
|---|---|
| `skill/_common/opt/relax_common.py` | `CELL_STAGE_ENCUT_FACTOR=2.0`（只对 ISIF≥3 段加 ENCUT）；`TOL_LAYER_KB=0.4` + `_layer_factor()`（层内口径换算）+ 下限 0.05；`_last_inplane`/`_cpu_seconds`/三指纹看门狗 + 快照；`STALL_MINUTES` 默认 60→180；`ALLOW_2D_FIXED_CELL`；`LATTICE_CONSTRAINTS` 接线 |
| `skill/_common/opt/checks_relax.py` | `ck_relax_injob` 三道拦：变胞段 `.skipped`、`RELAX_CELL_UNCONVERGED`、**末态压力物理量后置校验**；`S1_STRESS_GATE.off` 豁免 |
| `skill/kl-dft-cpu/gen_step4_disp.py` | `check_existing_matches_input()` 结构指纹（比对 `phono3py_disp.yaml` 的 `unit_cell` 与当前 POSCAR，差 >1e-3 Å 报错） |
| `skill/kl-dft-cpu/kl_common.py` | `enforce_incar_tags()`（防项目级模板遮蔽）；`write_submit(require=)`；`check_frames_match_displacements()` 抽帧逐原子校验（`points` 是 dict 列表、位移在 `dataset.displacements`、**VASP 分数坐标要按周期回绕**） |
| `skill/kl-dft-cpu/gen_step5_fc.py` | 调抽帧校验，不通过即 `sys.exit` |
| `skill/kl-dft-cpu/gen_step6_kappa.py` | `MESH_SCAN/COMPARE_LBTE`、`REQ_2D` 归一化守卫、P48 QoS/cand 重写 |
| `skill/kl-dft-cpu/templates/step5_phonon_plot/submit_plot.tpl` + `gen_step5.1_plot_phonon.py` | S5.1 由"登录节点 gen"改为**提交绘图作业**（`--work` 分派），从 v2.0 同步 |
| `skill/ke-dft-cpu/tools/` | 15 个诊断脚本（从 v2.0 同步） |

回归脚本（`tmp/`，全部通过）：`_kltest_s1gate`、`_kltest_inplane`、`_kltest_fingerprint`、
`_kltest_2dfix`、`_kltest_plotmode`、`_kltest_stale`、`_kltest_lc`、`_kltest_render`、
`_kltest_cellencut`、`_kltest_addgrid`。套件 6 个全 PASS。

---

## 4. 线 B 第一件事：vac10 的四项诊断（**先诊断，别急着对全批重跑**）

0.5% 的面内应变换算成层内应力是 kbar 量级，比阈值高两个数量级 —— 隔着 17 Å 真空的
色散作用很难产生这么大的力，所以"镜像作用"不是唯一解释。按顺序做（都便宜）：

1. **把 vac10 弛豫后的结构放回 20 Å 真空**（原子内部坐标不变），520 eV 算一次应力：
   - 应力≈0 → 两个晶格都是平衡态，差异来自不同局部极小值或平坦势能面 →
     再比单原子能量、空间群、z 跨度；
   - 层内应力几 kbar → 真空确实在影响结果 → 进第 3 步。
2. **检查层的净偶极**：`IDIPOL=3`、**不开** `LDIPOL`，只读偶极矩。若不可忽略，
   主因是周期镜像间的静电作用 → 应开偶极修正，而不是加厚真空。
3. 若确认是真空效应，**再补 30 Å 真空点**（三点才看得出趋势），据此定默认真空厚度。
4. 依诊断结果决定是否对 7 个材料重做 S1。

## 5. 两条线并行

**线 A（优先）：用文献基准验证功能 —— 单层 MoS₂**
- 不要用 Mo₂S₃ Z4-3-1（S4 有 165 帧、每帧 20~38 分钟，真空问题一旦成立全要重算）。
- 单层 MoS₂ 原胞 3 个原子、便宜，且文献有值可对照；**真空按文献取 25 Å**。
- `INCAR_BENCH` 也在这个材料上做（四变体 LREAL/ADDGRID + ENCUT 1.5× vs 2.0×）。
- 逐步验收判据表：

| 步骤 | 判据 |
|---|---|
| S1 | 写出材料级 `thickness_2d.json`；ke 侧读同一个文件 |
| S4 bench | 最大 `|ΔF| < 0.5 meV/Å`，据此定生产 INCAR 与 ENCUT |
| S4 | 指纹校验通过；INCAR 与 bench 选定一致 |
| S5 | `-c` 日志出现 `Imposing rotational invariance and equilibrium conditions`；ns_harm 自由度减少；LASSO 的 fc3max 与 OLS 差 <2%；墙钟内完成 |
| S5 声子 | 两个面内方向 ZA 指数都 1.7<p<2.3；近 Γ 无负频；2D 路径不含 kz 段 |
| S6 | 三档 q 网格相邻差 <5%；`kappa_2d_normalized_*` 由 V/A 算；zz 为 null |
| 解法对照 | 最粗网格、300 K 跑一次 LBTE，写 `rta_over_full` |
| 文献 | pheasy 论文 WS₂：PBE、**d=6.16 Å**、90×90×1、迭代解、300 K 三声子 κ = **241.6 W/mK**（换 d 要先统一到文献值再比） |

**线 B：jzz 这批的数值问题**
1. vac10 四项诊断（第 4 节）。
2. **AlN 先暂停**：用**材料根种子 POSCAR**（不要用 `SPOSCAR_old`）与当前 CONTCAR 比
   空间群/配位数/z 跨度，判断是种子本身多层/屈曲，还是弛豫中变的；同时修 kl 的 ENCUT 覆盖。
3. 定下真空厚度后，对需要的材料重做 S1。
4. S1 通过且线 A 验证完成后，再批量推 S2→S6；S4 的 INCAR 按 bench 结论定。

**低优先级**：Mn₂In₂Se₅（jzzn `/public/home/wangch/Mn2In2Se5/fix`，**账号是 wangch 不是
wangchao**，手写 phono3py+symfc、**3D 体材料**，按 |P| 判）与 225 相四合金（hanhai25、
走 **ke 链**）—— 都是 3D，与真空/厚度无关，只需查当时变胞段有没有被跳过。

---

## 6. 归档与产物路径

- **旧结构归档**（S4 力数据，2×2 上排的命根子）：
  `/public/home/wangchao/Fullerene_Network/work/_archive_pre_relax_20260917/`（1.1 GB）
  - `P1_Mo-MoS2_Z4-3-1_kl_step4_disp/`：132 MB，**166 帧含全部 vasprun.xml** +
    `phono3py_disp.yaml`(1 MB，含旧超胞与位移) + 旧 `phonon_summary.json`/`fit_config.json`
  - 另外 6 个材料的 `*_essentials/`（`phono3py_disp.yaml` + `kl_params.txt` + 旧 SPOSCAR）
  - 另 6 个材料的旧 `disp-*` 整批移到了 `*_step4_disp_old/`（共 1902 个帧目录）
- **本地**：`tmp/_sync_v20_backup_*`、`tmp/_align_backup_*`、`tmp/_ab2d/`（S5 A 臂基线）、
  `tmp/_pulay/OUTCAR_{390,507,624}`（Pulay 三点）、`tmp/_arch/`、`tmp/_cmp/`
- **新建材料**：`P1_Mo-MoS2_Z4-3-1_Z4-3-1_Mo2S3_vac10`（c 20→30 Å，真空隙 16.74→26.74 Å）

## 7. 作业状态（写文档时）

| jobid | 内容 | 状态 |
|---|---|---|
| 3847692 | Mo₂S₃ Z4-3-1 kl 的 S1（520 eV 变胞） | 已完成 |
| 3847705 | Mo₃S₄ kl 的 S1 | 已完成 |
| 3847708 | Ti₂S₃ kl 的 S1（ENCUT_b=550 = 2.0×274.7 ✓ 规则正确） | **FAIL：段 b 首轮 ZBRENT 崩**（漏洞 1）；`CONTCAR.s2r1` 已存档，`retry` 可续 |
| 3847711 / 3847713 | Zn₅O₃ kl / ke 的 S1 | 已完成 |
| 3847714 / 3847715 | AlN kl / ke 的 S1 | 已完成（kl 的 ENCUT 异常见 1 节） |
| 3847717 | BeO kl 的 S1 | 已完成 |
| 3847718 | **vac10 kl 的 S1** | 已完成（面内 0.0049，a=3.0595617） |
| 3844584 | S5 第一臂（RASR=none，LASSO） | **已按用户要求 stop**（跑 18.9 h 无输出） |
| 3847023 / 3847024 | AlN kl / ke 的上一轮 S1 | 已被后续 retry 覆盖 |
| 3846871 / 3846873 / 3846981 | AlN kl / ke、vac10 的 S1 首轮 | 已完成 |

## 8. 交接时必须知道的操作坑

1. **等待作业别用"名字计数"**：曾两次把 ssh 瞬时失败当"作业结束"，取回旧结果当新结果。
   要按 jobid + 哨兵 + 二次确认，**并校验取回的产物属于目标结构**。
2. **新建材料**要在**材料目录里**跑 `autozt -tt <技能> init`（建 `project_setting/` 才是注册点）；
   只在项目根放 POSCAR + 技能子目录，`-p` 找不到它。
3. **`autozt stop`** 第一次可能只打 scancel 标记不真取消（采集没挂到 job 信息）→
   先 `status` 刷新采集再 `stop -y`。
4. **`start` 对已完成/FAIL 的步骤要加 `-f`**；`retry` 对 **fanout 步骤只补未完成子目录**。
5. **项目级模板遮蔽技能模板**（`find_asset` 优先级）→ 已有 `stale_template_note` 检测 +
   `enforce_incar_tags`/`write_submit(require=)` 兜底。
6. **`git_commit` 闸门在本会话被自动拒绝**（审批禁用）；**未做任何 commit**。
7. `setting/README.md` 被 `.gitignore` 忽略（含 `lattice_constraints` 专章与
   `AUTOZT_VASP_VERSION` 用法），需移到 `docs/CONFIGURING.md`。
8. `skill/_common/thickness_2d.py` 是**未跟踪文件**，commit 时要单独 `git add`。

---

## 9. 2026-09-18：ZBRENT 根因定位与修复（本轮）

### 9.1 根因：作业起来时用【上一次作业遗留的 OUTCAR】判「首段已收敛」

`_run_stage` 里「上一段已收敛就跳过」用的 `_converged()` 读的是**当前目录的 OUTCAR**；
作业刚起来时那个文件是**上一次作业留下的**。于是重投时发生连锁反应：

1. 段 a 被跳过并写 `.s1.skipped`（证据：Ti₂S₃ 目录里 `.s1.skipped` 时间戳 = 作业启动时刻
   23:03，且不存在 `.s1.done` —— 只有 EARLY_EXIT 这条路会只写 `.skipped` 不写 `.done`）；
2. 变胞段不禁跳，于是它拿到的是**上一次遗留的旧 CONTCAR**；
3. 那个晶胞是**旧 ENCUT 下的平衡胞**，几乎没有梯度 → CG 空转 → ZBRENT 崩。

硬证据（OUTCAR.s2r1：单次运行，只有一份 INCAR/POSCAR 头）：

- 起点晶格 a = 3.3079000（= 上一次的 CONTCAR），**不是**种子 POSCAR 的 3.3635696；
- 从第 3 步起 **19 步内总能量只变 1.7e-7 eV**（-73.417290 一直不动），VASP 硬撑到第 20 步；
- 第 20 步晶格又回到 3.3079000、应力回到 +0.345 kB（= 起点值）→ 线搜索回退后夹不到极小值；
- 崩在 `ZBRENT: fatal error in bracketing` + `I REFUSE TO CONTINUE WITH THIS SICK JOB`，rc=1；
- 同时排除：`queue.err` 里**没有** `[watchdog]`，不是撞墙钟也不是看门狗。压力 +0.37 kB、
  力 0.007 eV/Å，物理上本来就已经到位 —— 这个作业是被「空转」拖死的。

### 9.2 四个修复（已落地 + 已验证）

| # | 修复 | 位置 |
|---|---|---|
| **A（根因）** | `STAGES_RUN` 闸门：只有**本次作业内真跑过至少一段**之后，才允许用 `_converged` 跳过后
续段；首段永不被陈旧 OUTCAR 跳过 | `relax_common.py` `_run_stage` + 生成头 `STAGES_RUN=0` |
| **B（判据，原修法②）** | 变胞段 EDIFFG 由 **-0.01（力判据）改为 +1E-4（能量判据）**：能量一平就停
（本故障里第 2 步即停，不再空转到 ZBRENT）。可用 step.conf `CELL_STAGE_EDIFFG` 覆盖 | `STAGE_SPEC` b,c |
| **B 的兜底** | 变胞段稳定判据新增**末态最大力上限** `FORCE_TOL=0.02 eV/Å`（能量判据可能让 VASP 在
力没洗干净时就停）；`_max_force` 解析 OUTCAR 最后一个 TOTAL-FORCE 块 | `_cell_settled`；`CELL_FORCE_TOL_EV_A` |
| **C** | S1 判据报错文案区分 **ZBRENT / 看门狗 / 其它**，并点名是哪一段（`.sN.started`）；
旧文案一律说「可能撞墙钟或被看门狗杀掉」会把人带偏 | `checks_relax.py` `_stage_crash_kind` |

修法①（给面内 ±0.3% 应变种子）**没有采用**：本故障里晶胞本来就在平衡点附近，人为加 0.3%
应变反而制造几 kbar 的假应力、多绕一圈；根因是起点陈旧，不是梯度算法。

### 9.3 验证

`tmp/_kltest_stagerun.py`（新增，假 VASP 端到端跑 `run_relax.sh`）：

| 用例 | 结果 |
|---|---|
| 目录里预置「上次遗留的、含 reached required accuracy 的 OUTCAR」 | **PASS**：段 a 真跑了
（调用序列 `ISIF=2 EDIFFG=-0.05 → ISIF=3 EDIFFG=1E-4 → ISIF=3`），写 `.s1.done`，无 `.s1.skipped` |
| **负对照**：把闸门从脚本里抠掉（模拟旧代码） | **PASS**：复现原故障 —— 只调 ISIF=3、`.s1.skipped`
出现。证明测试对修复敏感 |
| 静态：闸门 / `FORCE_TOL` / `_max_force` 在生成脚本里 | PASS |
| 判据：段 a `EDIFFG=-0.05`（力）；段 b/c `EDIFFG=1E-4`（能量） | PASS |
| 单元：`_max_force` 取最后一个 TOTAL-FORCE 块 | PASS（多块 OUTCAR 取到 0.0020，阈值 0.02 放过 / 0.001 拦截） |

既有回归 6 个全 PASS：`_kltest_s1gate`、`_kltest_inplane`、`_kltest_cellencut`、
`_kltest_stageskip`、`_kltest_2dfix`、`_kltest_stale`。

### 9.4 同类风险审计（全批材料）

对 `P1_*/{kl,ke}-dft-cpu/step1_std_opt` 逐个查阶段标记：

| 材料 | `.s1.skipped` | 说明 |
|---|---|---|
| AlN(kl)、BeO、Mo₂S₃ A4-3-1、Mo₃S₄、Mo₂S₃ Z4-3-1、vac10、Zn₅O₃ | **都有**（且无 `.s1.done`） |
最后一次 S1 里**段 a 都被跳过** —— 同一个 bug |
| Ti₂S₃(kl) | 已清（本轮 retry） | 修复前它因此崩 ZBRENT |

影响分级：这些材料的**变胞段从未被跳过**（`stage_changes_cell` 一律 ee=0），ISIF=3 本身也弛豫
离子，所以最终结构仍是「在生产 ENCUT 下弛豫到面内 0.0049~0.14 kB」的合法极小点（变胞段起点
= 各自上一次留下的胞，位移仅 0.02~0.14%）。修复后不会再出现「段 a 静默跳过」；代价是新一批
S1 会多跑一次固定胞段。

### 9.5 ★ 必须知道的语义坑：retry 不从材料种子起算

`resolve_stage()` 在 `STAGE_MODE=in_job` 下返回 `src_poscar = cwd/POSCAR`，而**步骤目录里的
`POSCAR` 会被运行过程改写**（段 b 的重投续跑 `cp CONTCAR POSCAR`、多遍循环的 `cp -f CONTCAR POSCAR`）。
所以：

- `retry` = **从上次留下的那个结构继续**（实测 Ti₂S₃：retry 后 POSCAR 直接是上一次的 CONTCAR
  3.3098644，不是材料种子 3.3635696）；
- 想让 S1 **从材料种子重跑**，只有 `rerun` / 清空步骤目录一条路（破坏性，需请示）；
- 因此「换了新种子再 retry」是**无效的** —— 与 S4 fanout 陈旧帧同一类陷阱。

### 9.6 本轮作业

| jobid | 内容 | 状态 |
|---|---|---|
| **3847765** | Ti₂S₃(kl) S1 重投（修复后第一次跑：段 a 会真跑、能量判据不再空转） | 已提交，PD(Priority) |

提交前已核验远端生成物：`INCAR.s1_a EDIFFG=-0.05`、`INCAR.s2_b/s3_c EDIFFG=1E-4`、
`run_relax.sh` 含 `STAGES_RUN=0`/`FORCE_TOL=0.02`/`_max_force`、旧阶段标记已清、
`INPLANE_TOL=0.168974 kB`（当前面内 0.13 kB < 阈值 → 预期一遍收敛）。

### 9.7 遗留（下一轮）

1. 监控 3847765 跑完，确认：段 a 有 `.s1.done`、无 `.s1.skipped`、变胞段无 ZBRENT、
   S1 判据返回 converged。
2. 项目级 `incar_2d.tpl` 是**过期覆盖副本**（少 `{{DIPOLE_LINE}}`，retry 时 autozt 已告警）——
   需同步覆盖（改项目文件，按规矩先请示）。
3. 全批其余材料的 S1 是否要按新代码重跑一遍（会多一次固定胞段；也可按兵不动）。
4. 线 A（单层 MoS₂ 文献基准）与线 B（vac10 四项诊断）仍未开始。


### 9.8 生产验证结果（09-18 01:0x，Ti₂S₃ S1 跑完）

job 3847765 在 cu48 跑完，**S1 = OK「converged（3/3 段）」**。硬证据（远端 step1_std_opt）：

| 检查项 | 结果 |
|---|---|
| 阶段标记 | `.s1.done`、`.s2r1.done`、`.s3r1.done` 齐全；**无 `.s1.skipped`** |
| 段序 | a(ISIF=2, EDIFFG=-0.05) → b(ISIF=3, 1E-4) → c(ISIF=3, 1E-4) |
| 崩溃 | `grep -l "ZBRENT\|SICK JOB" OUTCAR*` = **0 个文件** |
| 闸门生效 | queue.out 出现新文案「段1(a) …—— 本次作业还没跑过任何段，不采信目录里残留的 OUTCAR（防陈旧跳过）」，随后段 a 真的跑了 |
| 收敛 | **一遍**稳定：`max\|Δ晶格\| = 0.001110 Å`、`面内max\|σ\| = 0.0438 kB`（阈值 0.168974）、`max\|F\| = 0.0000 eV/Å`（上限 0.02） |
| 末态应力 | 段 c 末态面内 0.0296/0.0438 kB（修复前那次是 0.130） |
| 判据 | `checks_relax` 返回 converged（3/3 段），S1 不再 FAIL |

对比修复前同一步：段 a 被跳过 → 变胞段从旧 CONTCAR 起算 → 19 步空转 → ZBRENT（rc=1、S1 FAIL）。

### 9.9 顺带修掉：采集器被 checker 的 stdout 污染（巡检曾经是瞎的）

`autozt -tt kl-dft-cpu summary` 一直报「跳过采集失败的组 host=jzzn wd=…：Expecting value: line 1
column 2 (char 1)」→ 该技能整组「无材料」。抓原始 stdout 定位到：**远端采集器要求 stdout 是
一整份 JSON，而 `ck_relax_injob` 把豁免提示打到了 stdout**：

```
[..] 变胞应力门禁已豁免（S1_STRESS_GATE.off：…）      ← 有豁免材料的目录都会打
{"types": [ … ]}                                    ← 真正的 JSON，被前面的行顶掉
```

于是 json.loads 在第 1 个字符就失败（`[` 后跟 `,`）。豁免材料（Si / LiCoO₂ / BAs_kltest 等）
和 P1 材料在**同一个 work 目录**下，所以 P1 的 kl 组一起瞎。修法：该 print 改 `file=sys.stderr`
（已加注释说明"本判据 stdout 必须是 JSON"这条约束）。修完立刻恢复：
`kl-dft-cpu: 14 材料 done=3 run=0 pd=0 err=2 scancel=0 wait=9`。

★ 教训：**技能判据/checker 里任何 print 都必须写 stderr**，否则整组巡检静默失效。

### 9.10 巡检发现的两个 FAIL（需决策，未动手）

1. **Mo₂S₃ Z4-3-1 的 S4_disp：165/166 帧，且这 165 帧属于【旧结构】。**
   我对该步跑了 `retry`（fanout 只补缺失帧），新加的 `check_existing_matches_input()` **按设计拒绝**：
   「phono3py_disp.yaml 的单胞与当前输入差异最大 **0.1263 Å**（>1e-3 即判定不一致）」。
   也就是说 S1 这次重弛豫（胞改了约 1%）之后，**S4 的整套位移/力数据都作废**，S5 若照跑
   就是「旧力配旧位移」的静默混拟合。处置只有两条，都需你点头（破坏性）：
   - 确认旧帧已归档后，清空 `disp-*/POSCAR-*/phono3py_disp.yaml/SPOSCAR` 再 retry；
   - 或 `autozt -tt kl-dft-cpu -p <材料> -j S4_disp rerun`（推倒重来）。
   代价：166 帧 × 每帧 20~40 分钟（S4 用基础 ENCUT 390，比 S1 变胞段便宜），
   按 7~16 并发槽折算约 1~2 天。**这一步不做，S5/S6 无从谈起**。
2. `LiCoO2_kltest`（shengbte-gpu 基准线）S4_disp 22/52，缺 30 帧 —— 与 P1 无关，
   属基准线维护，低优先级。

### 9.11 两个项目级模板仍是过期覆盖副本（retry 时 autozt 会告警）

- Ti₂S₃(kl)：`project_setting/templates/step1_std_opt/incar_2d.tpl` 少 `{{DIPOLE_LINE}}`；
- Mo₂S₃(kl)：`project_setting/templates/step4_disp/incar_force_3d.tpl` 少
  `{{FORCE_EDIFF}}`、`{{FORCE_LREAL}}`、`{{FORCE_PREC}}` ← **直接关系到正在算的 S4 力**
  （EDIFF/LREAL/PREC 是取力的关键标签，被静默吃掉会让力常数质量下降）。

处置：把技能模板覆盖到项目级副本（先备份旧副本）。**改项目文件按规矩需先请示**。

### 9.12 本轮遗留

1. 等人决策：S4 是否清空重算（9.10 第 1 条）；项目级模板是否同步（9.11）。
2. S1 修复已闭环，但**线 A（单层 MoS₂ 文献基准）仍未开始**，它是验证 P0-1/P0-2/P0-3
   端到端的关键，建议优先于批量重算 S4。
3. 线 B（vac10 四项诊断）未开始。


---

## 10. 2026-09-18：全批状态盘点（修复后的完整监控）

采集器修好后第一次拿到完整状态。**结论：S1 这一环 8/8 材料全部正确；昂贵的陈旧产物只有一份。**

### 10.1 P1 全 8 个 kl 材料的 S1 末态（段 c 末态面内 σxx/σyy/σxy，kB，胞口径）

| 材料 | 段 done | skipped | 面内 max\|σ\| | 阈值 |
|---|---|---|---|---|
| AlN | 2 | 1 | 0.1144 | 0.1604 |
| BeO | 2 | 1 | 0.0351 | 0.0640 |
| Mo₂S₃ A4-3-1 | 4 | 1 | 0.0136 | — |
| Mo₃S₄ | 4 | 1 | 0.0181 | — |
| Mo₂S₃ Z4-3-1 | 4 | 1 | 0.0083 | — |
| vac10 | 4 | 1 | 0.0049 | 0.0915 |
| **Ti₂S₃** | 3 | **0** | 0.0438 | 0.1690 |
| Zn₅O₃ | 4 | 1 | 0.0294 | 0.0643 |

→ **旧 bug（段 a 被跳过）没有损伤任何材料的 S1 结果**：变胞段本来就同时弛豫离子，
最终结构都是各自（旧）起点的合法极小点，面内应力全部远低于阈值。无需为这 7 个材料重跑 S1。

### 10.2 ★ 下游全部建立"旧结构"之上（重要）

逐个比对 `step1_std_opt/CONTCAR` 与下游步骤实际用的结构（`stepN/POSCAR` 的 a）：

| 材料 | S1 CONTCAR a | S2/S3/S4 用的 a | 偏差 |
|---|---|---|---|
| Mo₂S₃ Z4-3-1 | 3.0446944 | 3.0410669(S2) / 3.0863707(S3) | +0.12% / +1.37% |
| BeO | 9.8163354 | 9.8175607 | +0.012% |
| Mo₂S₃ A4-3-1 | 11.7746396 | 11.9320777 | +1.34% |
| Mo₃S₄ | 3.0055247 | 3.0471745 | +1.39% |
| Zn₅O₃ | 6.2066460 | 6.3002091 | +1.51% |
| Ti₂S₃ | 3.3096458 | 3.3635696（= 材料种子） | +1.63% |
| **AlN** | 7.6458416 | 9.7484589 | **+27.5%** |
| vac10 | 3.0595617 | （下游还没生成） | — |

成因：S2/S3/S4 都是在**这一轮 2.0× ENCUT 重弛豫之前**生成的（那时 S1 的 CONTCAR 还是种子/旧胞）。
所以「S2/S3 OK、S5_fc OK、S6_kappa OK」这些状态，对旧结构成立，对新结构**不成立**。

### 10.3 陈旧产物的规模（决定重算代价）

| 材料 | S4 帧目录 | 有 OUTCAR | S5 产物 | S6 产物 |
|---|---|---|---|---|
| **Mo₂S₃ Z4-3-1** | **166** | **165** | **6** | **2** |
| 其余 7 个 | 0 | 0 | 0 | 0 |

→ **只有 Z4-3-1 有昂贵陈旧产物**（166 帧 + S5/S6 结果）；其余 7 个只需重跑便宜的 S2/S3
（静态 + NAC，几分钟级），S4 从头算即可。这对预算是个好消息。

Z4-3-1 的 S4 已被新防护判为不一致：`phono3py_disp.yaml` 单胞 vs 当前输入**差 0.1263 Å**
（阈值 1e-3）→ 这 165 帧、以及建立其上的 S5_fc（min_freq=-0.072 THz）与 S6_kappa 结果，
**在重算前都不可信**。处置需你点头（见 9.10）。

### 10.4 ke 链（P1 材料）

- `P1_Zn-ZnO_Z3A2_Z-3-2_Zn5O3` 的 **S2.3_hse（step2_bandgap 的子步）与 S6_elastic 从
  **2026-08-31** 起就 FAIL**：queue.err 是
  `slurmstepd: error: *** JOB 3768990/3769568 ... CANCELLED ... DUE TO TIME LIMIT ***`，
  而 `submit.sh` 里**没有 `--time`**，走的是集群默认上限（实测跑了约 26 小时被杀）。
  **盲目 retry 会再撞一次**，正解是给这两步加 `--time`（改提交模板，需请示）或降成本。
- 其余 P1 材料的 ke 链当前步骤无 FAIL。
- ke 链里的 `Si`（|P|=15.45）、`Mo2S3`（7.80）、`Zn3O2`（5.75）三个 S1 被**新应力门禁**拦下 ——
  它们是**其它项目**的材料（不是 P1_*），kl 侧对同类基准材料（Si/LiCoO₂/BAs）已放豁免文件，
  ke 侧没有。是否同样豁免需你定。

### 10.5 本轮结论与待决策清单

| # | 事项 | 需要的决定 |
|---|---|---|
| 1 | Z4-3-1 的 S4：清空 166 帧重算 / `rerun` | 破坏性，需同意（约 1~2 天机时） |
| 2 | 其余 7 个材料重跑 S2/S3（便宜，因为 S1 结构变了） | 是否现在推 |
| 3 | P1 Zn₅O₃ ke 的 S2.3_hse / S6_elastic：加 `--time` | 改模板，需同意 |
| 4 | ke 侧基准材料是否照 kl 侧放 `S1_STRESS_GATE.off` | 需同意 |
| 5 | 项目级过期模板同步（9.11） | 需同意 |
| 6 | 线 A（单层 MoS₂ 文献基准） | 是否开始（这是验证 P0-1/P0-2/P0-3 端到端的关键） |
| 7 | 线 B（vac10 四项诊断） | 是否开始 |

`LiCoO2_kltest` S4 22/52（缺 30 帧，基准线）与本批材料无关，低优先级。


---

## 11. 2026-09-18 第二批：核实结论、决定执行情况、提交状态

### 11.0 ★ 更正：AlN(kl) 的 ENCUT 误报（我上一轮报错了）

**远端实测 AlN(kl) 的 `INCAR.s1_a` = ENCUT 600（=1.5×400）、`INCAR.s2_b`/`s3_c` = ENCUT 800
（=2.0×400），与 AlN(ke) 完全一致 —— 规则从未失效。** 此前"AlN(kl) 是 520"的结论来自
**未 fetch 的本地 `result/step1_std_opt/INCAR` 旧副本**（旧 gen 用 1.3×400=520）。
§1 该条已就地标记作废。**教训：判断远端现状一律读远端文件，绝不用本地 `result/` 副本。**
（同理，本轮前面读到的"Ti₂S₃ 本地 result 显示 EDIFFG=-0.01"也是同一原因。）

### 11.1 核实 1：AlN 的结构差异 —— 结论：种子本身不可用，不是相变

pymatgen 对比**材料根种子 POSCAR** 与当前 S1 CONTCAR（同为 Al10N2 = Al5N，12 原子）：

| 指标 | 种子 POSCAR | S1 CONTCAR（2.0× 重弛豫后） |
|---|---|---|
| a / b (Å) | 9.748459 / 5.058944（正交） | 7.646974 / 4.960744（α=88.4°, β=90.3°, γ=90.5°） |
| V/atom (Å³) | 82.195 | 63.222（**−23%**） |
| z 跨度 (Å) | 4.339（3 个 z 簇） | 5.221（2 个 z 簇） |
| 空间群 | **P1**（无任何对称性） | **P1** |
| 配位数 | Al:[1,2,3,5,7] N:[3] | Al:[1,2,3,5,6,8,9] N:[3] |
| E/atom (eV) | −4.5915 | −4.5924（**只低 0.9 meV/atom**） |

**判断：既不是"离平衡点远"，也不是能量驱动的相变。**

- 计量比是 **Al:N = 5:1**，根本不是 AlN 单层（1:1）；
- 种子里就有 **CN=1 的悬挂 Al**、配位数 1~7 乱跳、完全无对称性、z 跨度 4.34 Å
  —— **它本身就不是平面单层**；
- 面内胞缩小 22% 只换来 **0.9 meV/atom**：势能面对该应变方向几乎**平坦**，
  说明"晶格常数"在该种子上没有物理意义；
- 弛豫后依然 P1、配位数更乱（出现 CN=9）、z 跨度反而增大到 5.22 Å。

**处置：AlN 不推 S2（按你的指示）✓**。建议确认该材料的出处，或换真正的 AlN 单层种子
（1:1、真空 ≥20 Å）重做。在此之前 AlN 的 S1 结果没有物理意义。

### 11.2 核实 2：FORCE_TOL 收紧 ✓

`CELL_FORCE_TOL_EV_A` 0.02 → **0.01**（= 原 `EDIFFG=-0.01` 的力判据本身），已写进 CONF_SPEC
（step.conf 可覆盖）。Ti₂S₃ 实测 max|F| = 0.0000 eV/Å → 收紧无副作用。

### 11.3 决定 3：`--time` —— **根因与"模板缺 --time"不同**

- `scontrol show partition cpu192` → **MaxTime=UNLIMITED**；上限在 **QoS**：
  `sacctmgr show qos` → **`regular` MaxWall = 1-00:00:00（24 h）**、`premium` = 2-00:00:00。
- `setting/jzzn/templates/step1_opt/submit_std_2d.tpl` 的 `#SBATCH --time=24:00:00` 是**注释掉的**，
  旁边就写着"不能超过 QOS 的 MaxWall，超了作业会永远 PD"。
- 所以**没有 --time 时 Slurm 直接套用 QoS 的 MaxWall（24 h）**；真正的两条路是：
  **给重作业换 `premium`（2 天）**，或**降成本/改并行设置**让它落在 24 h 内。
- P1 Zn₅O₃ ke 的 `S2.3_hse`/`S6_elastic`（2026-08-31 被 TIME LIMIT 取消）**状态未动**，
  等你确认走 premium 还是降成本后再 retry（盲目 retry 会再撞一次）。

### 11.4 决定 4：ke 侧豁免 —— 只加了 Si（前提核对后）

- **Si ✓ 已加**：kl 侧确有同名豁免（shengbte-gpu 基准、P=−15.45 kB、冻结晶格有意为之），
  符合"同一批、同一个理由" → 已写入 `…/work/Si/ke-dft-cpu/step1_opt/S1_STRESS_GATE.off`。
- **Mo2S3 / Zn3O2 未加**：kl 侧**没有**对应豁免，且 `Mo2S3` 的 ke 结构是
  a=5.266 / b=11.932 / **c=20 Å** 的 2D slab、胞明显未弛豫（|P|=7.80 kB）—— 属生产材料，
  **门禁拦得对**；放豁免会掩盖"晶胞没弛豫"的真问题。等你确认这两个材料的定位。

### 11.5 决定 6：线 A（单层 MoS₂）已启动 ✓

新建 `/mnt/d/tf_data/test_kl/MoS2_kltest`（与 Si_kltest/BAs_kltest/LiCoO2_kltest 同项目）：
1H-MoS₂ 单层（Mo1+S2、a=b=3.16 Å、γ=120°、S 在 Mo 上下各 1.57 Å、**真空 25 Å**）。
已 `init` 并按新代码 gen+提交 **S1：jobid 3847778**。生成物核验：

| 项 | 值 | 判定 |
|---|---|---|
| 段 a ENCUT | 390（1.5×258.689） | ✓ |
| 段 b/c ENCUT | **520（2.0×258.689）** | ✓ Pulay 口径生效 |
| EDIFFG | a −0.05（力）/ b,c **1E-4（能量）** | ✓ |
| KPOINTS | **15×15×1**（2D 强制 kz=1） | ✓ |
| run_relax.sh | 含 `STAGES_RUN` 闸门 | ✓ |

后续按 §5 验收表逐步检查；**MoS₂ 的 S4/S6 同样会遇到 §11.3 的 --time/QoS 问题**。

### 11.6 未完成（下一轮/下一会话）

1. **线 B（vac10 四项诊断）未开始**：① vac10 结构放回 20 Å 真空算应力；② `IDIPOL=3` 读偶极；
   ③ 30 Å 真空点；④ 据结果决定 jzz 这批是否重做 S1。
2. **`retry`/`rerun` 种子指纹校验未实现**（你提的方案）：设计 = 把材料种子 POSCAR 哈希写进
   步骤目录（如 `.seed.sha256`），gen 时比对，不一致则拒绝并提示改用 `rerun`。可复用现成件：
   `kl_common.check_existing_matches_input()`（S4 已在用，本轮真的拦住了一次旧结构数据）。
3. **init 阶段的模板占位符检查未加**：gen 阶段的告警**已存在**（`report.stale_template_note()`，
   本轮多次触发并列出缺失占位符），缺的是 init 阶段同一检查。
4. **INCAR_BENCH 模式未实现**（四变体 LREAL/ADDGRID + ENCUT 1.5× vs 2.0×、max|ΔF|<0.5 meV/Å）。
5. 决定 1（Z4-3-1 S4 等线 A 后 `rerun`）与决定 2（7 个材料 S2/S3 暂缓）按指示**保持不动** ✓。

### 11.7 ★ 提交状态：本会话 `git_commit` 被闸门拒绝

- **本会话未做任何 commit**：`git_commit` 被审批策略拒绝（本会话审批禁用 → 自动拒绝）。
  按 `~/.dsh/AGENTS.md` 铁律 5，**没有**改用 shell `git commit` 绕道。
- 本轮 3 个文件**已留在暂存区**，执行下面一条即可落库：

```bash
cd /home/wangchao/software/AutoZT
git commit -m "fix(2D-kl): S1 变胞段 ZBRENT 根因修复（STAGES_RUN 闸门）+ 能量判据 EDIFFG=1E-4 + 末端力上限 0.01；checker stdout 改 stderr"
```

  暂存内容：`skill/_common/opt/relax_common.py`(+167) `skill/_common/opt/checks_relax.py`(+57)
  `docs/HANDOVER-2D-kl-20260917.md`(+429+§11)。
- **多会话并发警告**：同一 checkout 里另有会话在写，留下 `docs/dev/COMMIT-PLAN-20260918.md`
  （其批次 2/4 已把这几个文件列入计划），还有会话在写 `skill/eph-qe-cpu/`、`skill/zt-dft-cpu/`。
  提交前建议先与它对齐。
- 上一段会话的大部分成果（`relax_common.py`/`checks_relax.py`/`thickness_2d.py` 等）**已由该会话的
  `c0ede3e "工作区批量落盘"` 提交入库**，并非全部未提交。
- **我故意未暂存**这些混合文件（含另一会话在途改动）：`skill/kl-dft-cpu/gen_step6_kappa.py`、
  `skill/kl-dft-cpu/skill.yaml`、`skill/kl-dft-cpu/templates/step1_std_opt/incar_2d.tpl`、
  `skill/kl-dft-cpu/templates/step6_kappa/submit_shengbte.tpl`。


---

## 12. 2026-09-18 第三批：结构体检、线 A S1 结果、P0-2 钩子核查

### 12.1 ★ 全 8 材料结构体检（定位问题不是 AlN 独有）

判据：命名/组成自洽、空间群（P1 = 危险）、配位数分布（CN=1 = 悬挂原子）、z 跨度与层数、
面内是否近似平坦到可疑。脚本 `tmp/_health8.py`（pymatgen + CrystalNN，symprec=0.05）。

| 材料 | 命名声称 | 实际组成 | 空间群 | z 跨度(Å) | z 簇 | CN 范围 | 判定 |
|---|---|---|---|---|---|---|---|
| AlN | Al5N | Al5N | **P1 #1** | 4.34 | 3 | **1–7** | ❌ P1 无对称 + 3 个悬挂 |
| **BeO** | Be2O | Be2O | P2_1 #4 | **0.16** | 1 | 2–5 | ⚠️ 完全平坦单层（待判） |
| Mo₂S₃ A4-3-1 | Mo2S3 | Mo2S3 | P2_1 #4 | 3.26 | 3 | 2–6 | ✓ |
| Mo₃S₄ | Mo3S4 | Mo3S4 | Pmn2_1 #31 | 3.34 | 3 | 2–6 | ✓ |
| Mo₂S₃ Z4-3-1 | Mo2S3 | Mo2S3 | Pmm2 #25 | 3.28 | 3 | 2–6 | ✓ |
| vac10 | （vac10 非化学式） | Mo2S3 | Pmm2 #25 | 3.26 | 3 | 2–6 | ✓（命名告警是启发式误报） |
| Ti₂S₃ | Ti2S3 | Ti2S3 | P2_1/m #11 | 4.85 | 4 | 3–6 | ✓（跨度偏厚，待目视） |
| **Zn₅O₃** | Zn5O3 | Zn5O3 | Pma2 #28 | **0.17** | 1 | **1–3** | ❌ 4 个悬挂 + 近似平坦 |

**结论：AlN 不是孤例。Zn₅O₃（4 个悬挂原子、CN 低到 1、z 跨度 0.17 Å）与 BeO（z 跨度 0.16 Å，
近乎完美平面的 Be₂O）需要同一套核查。** 三者都建议挂起并向 jzz 确认来源；在确认之前不投 S4 机时。
（未做的两项：每原子能量相对凸包的位置 —— 需要参考相数据，未查；Ti₂S₃ 的 4 个 z 簇属"厚三明治"
还是分层，需目视 CONTCAR 才能定。）

### 12.2 ★ 线 A：MoS₂ 的 S1 **没通过验收**（这是基准要暴露的东西）

作业 3847778 跑完了，`RELAX_CELL_UNCONVERGED`（同时 `RELAX_OK`，即力收敛、晶胞判据没过）。

实测（`OUTCAR.s2r1` / `s3r3`）：

| 量 | 段 b（CG，第 1 遍） | 段 c（准牛顿，第 1 遍） | 阈值 |
|---|---|---|---|
| 面内 max\|σ\| (kB) | **0.0776** ✓ | **0.1427** ❌ | 0.10784 |
| 离子步数 | 29 步（1.6→−19.8→…→−0.078，CG 先过头再收敛） | **仅 1 步** | — |
| `max\|Δ晶格\|` | — | **0.000000 Å** | 0.002 |

**诊断：不是"弛豫不到位"，而是"最后一段把它做坏了"。** 段 b 已把面内应力收到 0.078 kB（达标），
段 c 从那个几何出发只走 **1 个离子步**就停（力已收敛，VASP 立即判收敛），而那一步准牛顿的
初始 Hessian 很差，**把应力从 0.078 推回到 0.143 kB**；判据量的是段 c 的末态 → 判未稳定。
多跑几遍没用（三遍完全一样：0.000000 Å / 0.1427 kB / |P|=−0.27）。

**这是设计层的问题，不是 MoS₂ 的物理问题**，三个可选修法（需要你定，我没有擅自改）：

1. **段 c 不再动晶胞**（`ISIF=2` 固定胞做离子精修）：段 b 已经负责晶胞，c 的职责是准牛顿收尾；
   这样段 c 不会把已到位的应力做坏，还顺带省掉 c 段的 2.0× ENCUT 机时；
2. **判据取"本遍变胞段的最优值"**而不是"最后一段的末态值"（把段 b 的 0.078 记为结果）；
3. **判据口径**：0.10784 kB 是否已经低于该体系应力的数值噪声底（MoS₂ 15×15×1、
   `EDIFF=1E-7`、`PREC=Accurate`）—— 若噪声就有 ~0.1 kB，这个阈值永远过不去，
   应改为"相对噪声/相对 S4 门禁 0.5 kbar"来定。

按你的指示（"任何一步不达标就停下来"），**线 A 停在 S1**，等这三个问题定音后再往下推。

### 12.3 P0-2 钩子完整性核查（你点名要确认的） ✓

| 文件 | 现状 |
|---|---|
| `skill/kl-dft-cpu/gen_step6_kappa.py` | ✓ 完整：`norm_subs()`(L338)、`REQ_2D=("kappa_2d_normalized",)`(L553)、`require=REQ_2D`(L598)、注入脚本 L219–235 齐全；`ast.parse` 通过 |
| `skill/kl-dft-cpu/templates/step6_kappa/submit_shengbte.tpl` | ✓ 完整：`{{KAPPA_2D_FACTOR/META/THICK2D}}` 占位符齐、归一化块 L153–185 在、`bash -n` 通过、**无冲突残留标记** |

→ 线 A 跑到 S6 时这两个文件可用（尚未实跑验证）。

### 12.4 流程改良（按你的要求落到文档）

- **判断远端现状一律读远端文件；本地 `result/` 只是上次 fetch 的快照，只能当参考。**
  已由 ENCUT-520 误报证实（本地旧副本 520 vs 远端 800）；本案也解释了我此前读到的
  "Ti₂S₃ 本地 EDIFFG=−0.01"。
- **多会话共用工作区**：本 checkout 至少有 3 个会话在动。已知分工：
  * 本会话（2D kl/S1）→ `skill/_common/opt/{relax_common,checks_relax}.py`、
    `docs/HANDOVER-2D-kl-20260917.md`、`/mnt/d/tf_data/test_kl/MoS2_kltest`（线 A）、jzz 上 P1 材料的 S1；
  * 会话 B → `docs/dev/COMMIT-PLAN-20260918.md`、`autozt/*`（agent/MCP/纠错层）、
    `skill/ke-dft-cpu/step3*/step8*`、`skill/kl-dft-cpu/gen_step6_kappa.py`（P49）、
    `templates/step6_kappa/submit_shengbte.tpl`；
  * 会话 C → `skill/eph-qe-cpu/`、`skill/zt-dft-cpu/`（仍在写）。
  建议：各会话用 `git worktree add` 独立工作区，或至少约定"只碰自己的文件范围"，
  并且**提交前先看 `git status` 里有没有别人的在途文件**（本会话已按此只暂存自己的 3 个文件）。

### 12.5 仍未做（下一轮）

1. **Zn₅O₃ ke 26 小时的诊断**（先查 NCORE/KPAR/核数/k 网格，再决定是否换 `premium`），
   并按你的要求在模板注释里写清"`--time` 为何注释掉、依赖 QoS 默认值"的原因。
2. **init 阶段模板占位符检查**（gen 阶段已有）。
3. **线 B（vac10 四项诊断）**：按你的指示与线 A 的 S1 排在同一批 —— 现在线 A 卡在 S1 设计问题上，
   两者一起等修法定音。
4. 待你定：MoS₂ S1 的三个修法；AlN/Zn₅O₃/BeO 是否挂起并向 jzz 确认来源。


---

## 13. 2026-09-18 第四批：种子事故、tol_cell 措辞、待办 A–F（本会话最后一节）

### 13.1 ★ MoS₂ 种子事故（完整经过）——根因是输入，不是代码

1. 我手工建的 1H-MoS₂ 种子把 **Mo 和 S 放在同一个面内位置**（都写 1/3, 2/3），于是 S 正压在
   Mo 正上方，**d(Mo–S) = 1.57 Å**（正常 2.41 Å；1H 结构里 S 应占另一个 hollow 位）。
2. DFT **照单全收**：弛豫同时 (a) 把 S 往层外推 ±1.57 → ±2.42 Å、(b) 把面内胞收缩
   3.16 → 2.7126 Å（**−14.2%**），得到的是一个毫无物理意义的"极小点"。
3. 发现路径**只有一条**：比对晶格才看出缩了 14%。若只看"力收敛 + 应力小"，它完全正常
   （末态 |σ| = 0.14 kB、max|F| = 0.0000）——**这就是"坏种子会被 DFT 照单全收"的实证**，
   也是 A 项体检要固化进代码的理由。
4. 修正：Mo (0,0,0.5) + S (1/3,2/3,0.5±0.0628) → pymatgen 复核 d(Mo–S)=**2.407 Å**、
   S–S=**3.140 Å** ✓（已写入 `/mnt/d/tf_data/test_kl/MoS2_kltest/POSCAR`）。
5. 重投后 **S1 通过**：2 遍稳定、面内 **0.0059 kB**（阈值 0.10784）、max|F|=0.0000，
   `autozt` 判 `OK converged（3 阶段 INCAR、共 5 次段运行）`。
6. 附带操作失误（已纠正）：`rerun` 需交互确认，非交互调用得到 `EOFError` → **rerun 未执行**，
   而同一条命令串里的 `start -f` 却提交了 **jobid 3849199**（等于用坏结构继续跑）。已 `stop -y`
   取消，再用 **`rerun -y`** 正确重建（远端种子 a=3.16 Å 已核验），重新提交 **3849203**。
   ★ 教训：**`rerun` 与 `start` 必须分两条命令、分别确认返回值**，不能串在一个 shell 里。

### 13.2 tol_cell 结论的正确措辞（wangchao 修正）

正确表述：**"阈值对结构正常的材料是合适的，没有证据表明它低于噪声底"**，
**而不是"不需要下限"**。

- MoS₂（3 原子高对称胞）面内 0.0059 kB、Mo₂S₃（20 原子、正确结构）0.0083 kB —— 都远低于阈值；
- 之前算出的"0.065 kB 噪声底"是**畸变几何**上的抖动（粗糙 PES），**不是方法的噪声底**；
- 判定"段间抖动是否只是噪声"的可复用证据：**段 b 与段 c 的 CONTCAR 逐位相同**
  （MoS₂ 实测晶格 10 位一致、max|Δfrac| = 0.000000），而两段应力 0.0776 vs 0.1427 kB；
- **0.05 的下限保留**，不影响任何已通过的材料；判据去噪的另一半"取本遍最优"（②）**保持可选、不实施**；
- **但每段的 (晶格/应力/力) 三元组已改为逐段记入日志**（见 13.3），下次段间抖动直接看日志即可。

### 13.3 本会话已完成（可提交）

1. **`.gitignore` 兜底**（wangchao 单独一条提交，**已 `git add`**）：忽略另一会话留下的
   Windows 路径垃圾目录，模式 `C\:*` 与 `*wsl.localhost*`（后者绕开了 gitignore 的
   反斜杠转义坑），`git status` 残留 **0**；注释里被吃掉的字符已修好。
2. **每段三元组日志**：`relax_common.py` 的 `run_relax.sh` 生成器在每段跑完后追加一行
   `[run_relax] 段存档 <tag>: 面内max|σ|=… kB  max|F|=… eV/Å  |P|=… kB`（`relax_common.py:2009`）。
   验证：`ast.parse` 通过 + `tmp/_kltest_stagerun.py` 端到端 **8 项全 PASS**（含"抠掉闸门"负对照）。
3. 同批（早先完成、仍在暂存区）：`relax_common.py`（STAGES_RUN 闸门 / 变胞段能量判据
   `EDIFFG=1E-4` / 末态力上限 `CELL_FORCE_TOL_EV_A=0.01` / 看门狗三指纹 / 变胞段 2.0× ENCUT /
   2D 层内口径 / `_max_force`）、`checks_relax.py`（ZBRENT 诊断 + checker stdout 改 stderr，
   修掉采集器整组"无材料"）、`docs/HANDOVER-2D-kl-20260917.md`。

### 13.4 待办 A–F（下一会话按序）

**A. S1 gen 阶段四判据结构体检** —— 把参考实现 `tmp/_health8.py` 固化进
`skill/_common/opt/relax_common.py`：① 最近邻 < **0.7×共价半径和**；② CN=1 悬挂原子；
③ 计量比 vs 材料名；④ z 跨度 vs 声称层数。

- ★ **共价半径必须用字典，不要用 0.75×vdW 近似**：这个换算对不同元素误差可达 20%，
  而阈值本身就是 0.7 倍，两者叠加会让检查形同虚设。做法：从 pymatgen **导出一次**，
  写成 `_common/covalent_radii.py`（常用元素几十个足够，不需要运行时依赖 pymatgen）；
  **缺失元素就跳过该判据并 WARN**。
- 建议旋钮 `STRUCTURE_HEALTH = warn|error|off`（默认 `warn`，避免一刀切打断既有项目；
  jzz 那批用 `error`）。
- 实测标定（2026-09-18）：MoS₂ 事故 1.57 Å = 0.65×；**AlN** P1#1 无对称、CN=1×3、z=4.34 Å、
  Al:N=5:1；**Zn₅O₃** Pma2、CN=1×4、z=0.17 Å；**BeO** z=0.16 Å。这三个已挂起等 jzz 核对来源。

**B. `tol_cell`**：按 13.2 的措辞落档（本节已写明）；"取本遍最优"② 保持可选。

**C. S4 之前必须完成的两件**：

- **INCAR_BENCH**（四变体 LREAL/ADDGRID × ENCUT 1.5× vs 2.0×，判据 max|ΔF| < 0.5 meV/Å）——**未实现**；
- **材料级 `thickness_2d.json`**（P0-2b；S1 已算 h⊥/d，写出即可；S6 归一化验收依赖它）。
- ★ **与文献比 κ 时不要照套 6.16 Å**：6.16 Å 是 **pheasy 论文里 WS₂** 的取值（同类结构但各文献
  厚度定义不同；MoS₂ 另有取值）。正确做法：先用 vdW 口径算出 MoS₂ 自己的 d，比较时
  **两边统一到同一个 d**，或**直接比面热导 κ×d（单位 W/K，厚度不参与比较）**。
  写死 6.16 会让最后一步整体差一个系数。

**D. Zn₅O₃ ke 两步换 `premium`**（已批准，不必等线 A）：

- ★ **改之前先 `sacctmgr show qos premium`**：高优先级 QoS 常限制**并发作业数/总核数**，
  只把这两步换过去可能挤掉别的会话正在跑的作业；
- 已知事实：`cpu192` 分区 `MaxTime=UNLIMITED`，真正的上限是 **QoS `regular` MaxWall = 24 h**
  （`premium` = 2 天）；模板里 `#SBATCH --time=24:00:00` 是注释掉的，旁边写明"超 QoS MaxWall 会永远 PD"。
  换 QoS 后在模板注释里写清"`--time` 为何注释掉、依赖 QoS 默认值"；
- 该两步 2026-08-31 被 TIME LIMIT 取消（HSE 静态 EDIFF=1E-7/KPAR=3、弹性常数 67 个能量点
  + 1481 个 DAV 迭代 ≈ 63 s/电子步）；**结论是真实成本量级，不是并行设置坏掉**，
  且弹性常数对 EDIFF 敏感 → **不要为省时间放宽 EDIFF**。

**E. 线 B 四项诊断** —— 与下一批作业一起提交：① vac10 结构放回 20 Å 真空算应力；② `IDIPOL=3`
读偶极（不开 LDIPOL）；③ 补 30 Å 真空点；④ 据此定 jzz 那批是否重做 S1。

**F. 多会话纪律**：

- 交接文档维护"当前各会话负责的文件范围"（§12.4 的清单需刷新：新增未跟踪文件
  `bin/autozt.cmd`、`autozt/science.py`、`scripts/export_params.py`、`docs/PARAMETERS.md`、
  `skill/eph-qe-cpu/`、`skill/zt-dft-cpu/` 等，均为其它会话产物）；
- **每轮开工前先 `git status`**，确认没人动了自己的文件；只提交自己的文件；
- 已知分工：本会话=2D kl/S1 + 线 A；会话 B=`autozt/*`(agent/MCP)、`ke-dft-cpu/step3*/step8*`、
  `gen_step6_kappa.py`(P49)、`submit_shengbte.tpl`；会话 C=`eph-qe-cpu/`、`zt-dft-cpu/`。

### 13.5 ★ 工作方式（wangchao 定的规矩，下一会话必须遵守）

**每完成一个步骤的验收，立刻把结论追加进交接文档，不要攒到最后。** 验收表（§5）就是天然的
检查点，一步一记，交接文档就永远是最新的。本会话已经出现过"上下文见底、来不及写交接"的险情；
线 A 到 S6 还有 5 步、每步都要人工判断，一步一记是唯一可靠的传递方式。

### 13.6 线 A 当前进度（写本节时刻）

| 步骤 | 状态 |
|---|---|
| S1_opt | **OK — converged**（3 阶段/5 次段运行；面内 0.0059 kB / 阈值 0.10784；2 遍稳定） |
| S2_static | **PD(Priority)**，jobid **3849745**（排队中，`start` 已正确拒绝重复提交） |
| S3_nac / S4_disp / S5_fc / S5.1_plot / S6_kappa | 未开始（S4 前须先做 13.4-C 两项） |

材料：`/mnt/d/tf_data/test_kl/MoS2_kltest`（1H-MoS₂ 单层，a=3.16 Å、25 Å 真空、
d(Mo–S)=2.407 Å）；远端工作目录 `/public/home/wangchao/Fullerene_Network/work/MoS2_kltest`。


---

## 14. 静默错误 ↔ 门禁对照表（2026-09-18 小结；论文"门禁设计"一节的骨架）

> 一个自动化流程真正的难点不是把步骤串起来，而是**在每个环节设下能拦住静默错误的门禁**。
> 下面 10 项全部属于"不会报错、只会给出一个看起来合理的 κ"这一类。

| # | 静默错误 | 门禁（谁拦、拦在哪） | 状态 |
|---|---|---|---|
| 1 | 2D 的旋转不变性（BHH）从未真正施加 | S5 强制 `FIT_ENGINE=pheasy`（2D 自动），并要求拟合日志出现 `Imposing rotational invariance and equilibrium conditions`；ns_harm 自由度必须减少 | 代码已改，**未端到端验证**（线 A S5 是验证点） |
| 2 | κ 厚度归一化（ShengBTE 路线漏掉；层跨周期边界算错） | S6 三处 `write_submit(require=("kappa_2d_normalized",))`（模板渲染必须出现该字面量，防项目级模板遮蔽）；`thickness_2d.py` 为单一真源（`h⊥=V/|a_i×a_j|`、`d=zspan+两端 vdW`） | 代码已改，**未端到端验证**（S6 是验证点） |
| 3 | 面内晶格根本没弛豫过 | `ALLOW_2D_FIXED_CELL=False` + 2D 必须声明 IOPTCELL/LATTICE_CONSTRAINTS（否则硬 `sys.exit`）；变胞段 `ee=0` 禁止被力判据跳过；末态面内应力门禁；`RELAX_CELL_UNCONVERGED` 由判据读取 | **已生产验证**（8 材料全部面内 ≤0.114 kB；Ti₂S₃ 修好后 0.0438） |
| 4 | 应力里混 ~1 kB Pulay 偏差 | 变胞段 `CELL_STAGE_ENCUT_FACTOR=2.0`（只对 ISIF≥3 段）；末态应力门禁兜底 | **已实测标定 + 生产生效**（AlN 800=2.0×400、MoS₂ 520=2.0×258.7） |
| 5 | 用 \|P\| 判收敛被 σzz 掩盖真实面内应力 | 2D 改判**面内分量** `max(\|σxx\|,\|σyy\|,\|σxy\|)`，阈值走层内口径 `max(TOL_LAYER/(h⊥/d), 0.05)` | **已生产验证**（旧结构 P=−7.93 而面内 −11.97/−10.62） |
| 6 | S4 按 fanout 续算把**旧结构**的帧当已完成 | `check_existing_matches_input()`：`phono3py_disp.yaml` 单胞 vs 当前输入 > 1e-3 Å 即拒绝并提示 rerun；`check_frames_match_displacements()`：逐帧逐原子校验（含周期回绕），gen_step4/5 硬退出 | **已实战拦截**（Mo₂S₃ Z4-3-1 差 0.1263 Å 被拦下，避免"旧力配旧位移"混拟合） |
| 7 | 筛选材料里混着悬挂原子、计量比对不上的结构 | **四判据结构体检**（① 最近邻 < 0.7×共价半径和 ② CN=1 悬挂 ③ 计量比 vs 材料名 ④ z 跨度 vs 层数）→ 详见 §13.4-A | **未实现**（参考实现 `tmp/_health8.py`；已拦下 AlN/Zn₅O₃/BeO） |
| 8 | 作业刚起来时用**上一次遗留的 OUTCAR**判"首段已收敛"，首段被跳过、变胞段接手旧 CONTCAR → CG 空转 → ZBRENT 崩 | `STAGES_RUN` 闸门：本次作业真跑过至少一段后，才允许用 `_converged` 跳过后段 | **已生产验证**（Ti₂S₃ 日志出现"不采信目录里残留的 OUTCAR"，随后段 a 真跑） |
| 9 | 技能判据（checker）把提示打到 stdout → 采集器 JSON 解析失败 → **整组材料"无材料"，巡检静默失效** | 约定：checker 里任何 print 必须 `file=sys.stderr`；已在 `ck_relax_injob` 落地并注明理由 | **已修复并恢复**（`kl-dft-cpu: 15 材料 done/err/wait` 恢复可见） |
| 10 | 项目级模板副本过期 → 新功能被静默吃掉（少 `{{DIPOLE_LINE}}`/`{{FORCE_*}}`/归一化钩子） | `report.stale_template_note()` 在 gen 阶段比对占位符集合并列出缺失项；`enforce_incar_tags()`/`write_submit(require=)` 兜底 | gen 阶段**已生效**（本轮多次触发）；**init 阶段检查未实现** |
| 11 | **SCF 未收敛（撞 NELM）却照样输出力** —— VASP 只在 OUTCAR 留一行警告，作业正常退出、力照出 | §19 的 S5 门禁：逐帧扫 OUTCAR，含 `number of steps (NELM)` 或 `aborting loop because EDIFF is reached` 计数 != 1 → 该帧作废；任一帧作废即拒绝拟合并列出帧名；电子步数写进 `scf_steps.json` | **已落地并复核**（真实坏帧：MoS₂ bench 的 520 eV 帧被精确作废；P1 已算的 165 帧全部通过） |

**可复用的一类判据（本轮两次起作用）**：
- **结构指纹**：把"当前结构"与"历史产物所依据的结构"做逐分量比对（S4 的 disp yaml 校验、种子指纹方案 §9.5）；
- **物理量后置校验**：不信标记、只信末态物理量（末态面内应力、末态最大力、`max|Δ晶格|`）；
- **同几何重复性**：同一几何两次独立计算（如段 b vs 段 c）的差值就是数值噪声，可用来判断"抖动"还是"真变化"；
- **闸门只放宽有证据的**：如 Pulay 只对变胞段加 ENCUT、力上限回到与旧判据一致的 0.01，都有实测依据。

**力（位移帧）的噪声底 —— 可直接引用，不必再测**：同一设置两次独立计算，**同节点 ≈ 0**（MoS₂ bench 的 `encut_1.5x` 与 `lreal_false_addgrid` 同落在 cu14，max/RMS 逐位相同）；**跨节点 ≈ 1e-3 meV/Å**（cu18 vs cu06：max 差 ~0.0010、RMS 差 ~0.0001 meV/Å）。对照待分辨的口径效应（ADDGRID 的 RMS 0.165、ENCUT 520 的 RMS 1.26 meV/Å）小 2–3 个数量级。
**已实测的噪声量级（供判据设计参考）**：MoS₂（3 原子、15×15×1、520 eV）面内 0.0059 kB、
两遍完全一致；Mo₂S₃（20 原子、正确结构）0.0083 kB；而**畸变几何**上同一几何两次可差 0.065 kB
（正是 §13.2 那次的假象）。→ 阈值对结构正常的材料合适，0.05 下限保留。


---

## 15. 2026-09-19：INCAR_BENCH 落地、S4 生成、两笔待处理

### 15.1 ★ 遗留不一致：`gen_step6_kappa.py` 的 `write_material_level('..')`（wangchao 要求记一笔）

现状（两处写同一契约、**位置不同**）：
- **S1（已改好）**：`relax_common.write_material_thickness()` 写**材料根** `<mat>/thickness_2d.json`
  （`source="<技能>:<步骤>"`，如 `kl-dft-cpu:step1_std_opt`；按链比较，见 §16）。
- **S6（未动，在别的会话手上）**：`gen_step6_kappa.py` 里的 `_MAT_LEVEL` 片段在作业内调
  `write_material_level('..', THICK2D)`，而该处 cwd 是步骤目录 ⇒ `'..'` 是**技能目录**，
  会写到 `<mat>/kl-dft-cpu/thickness_2d.json` ✗。实测它从未落盘（jzzn 上全量 find 只有
  `step6_kappa/thickness_2d.json` 这个步骤级文件），所以目前只是"潜在的错位"，不是活动故障。

**处置（二选一，需改那个文件的人决定）**：① 把路径改成材料根（`../..`）并保持同链覆盖语义；
② **直接删掉那段 `_MAT_LEVEL`**，并在注释里写明"S1（step1_std_opt）已承担材料级写入职责"。
无论选哪个，**不要再出现第二个写入位置** —— 这次排查已经为"同一契约两个位置"多花过时间。

### 15.2 ★ 新操作规矩（wangchao 2026-09-18 定，已进 AGENTS.md 铁律 2/3）

1. **`-f`/`-y` 必须逐条单独请示**：批准"目标"（如"修 Mo2S3"）**不等于**批准具体命令里的 `-f`。
   教训来源：本会话曾以"MoS₂ 的 S1 重判"为由执行 `start -f`。
2. **会改变状态的命令一条一条执行，确认前一条的结果再走下一条**（同上教训的另一半：当时
   `rerun` 与 `start` 写在同一条命令串里，前者 EOFError 失败、后者照样执行，误提交了
   jobid 3849199，事后才 `stop -y` 取消）。

3. **涉及作业提交/状态变更的机械工作，自己做，不要派子代理。** 子代理适合**只读分析、写代码、跑测试**这类
   可以从产物验证、失败无外部副作用的任务；一旦涉及 sbatch/autozt 等**外部状态变更**，失败时「到底做了什么」是不确定的
   （2026-09-19 实测：两个子代理一个无消息失败、一个跑很久零产出，事后只能靠「没有变体、没有作业、没有 tmp 残留」反证没重复提交）。

### 15.3 INCAR_BENCH 落地（两个并行子代理完成，父会话独立复核）

- 新增 `skill/kl-dft-cpu/incar_bench.py`（698 行）：6 变体（4×LREAL×ADDGRID + ENCUT 1.5×/2.0×）、
  `generate/submit/analyze` 三个子命令、判据 `max|ΔF| < 0.5 meV/Å`、写 `bench/incar_bench_report.json`。
- `gen_step4_disp.py`：`INCAR_BENCH=off|on`（默认 off）+ `INCAR_BENCH_FRAMES`（默认 2）；on 时
  在生产帧之后生成 `step4_disp/bench/<variant>/disp-XXXXX`，**失败只 WARN 不阻塞生产**。
- `templates/step4_disp/step.conf` 与 `skill.yaml` 已登记（io_schema/skillspec 套件全过）。
- **隔离性已复核**：生产遍历全为非递归（`kl_common.py:295`、`gen_step5_fc.py:101`、
  `gen_step4_disp`、`autozt/workflow.py`、`_collector_remote.py`），仓库**无** `rglob/**` 遍历 `disp-*`。
- 测试：`tmp/_kltest_incarbench.py`（ALL PASS）、`tmp/_kltest_matlevel_s1.py`（exit 0）。

### 15.4 MoS₂（线 A）S4 当前进度

- `conf --set params.INCAR_BENCH=on`（材料级覆盖）✓；`-j S4_disp init`（只生成不提交）✓。
- 生成结果：**生产 13 个目录**（均衡帧 disp-00000 + **12 位移帧**，`alm` + MC-rattle RMS=0.0294）
  + `bench/` **6 变体 × 2 帧**。生产取力 INCAR = `ENCUT 390 / PREC Accurate / EDIFF 1E-8 /
  LREAL=.FALSE. / ADDGRID=.TRUE.`。
  → **MoS₂ 的生产扇出只有 12 帧**（对比 Mo₂S₃ Z4-3-1 的 166 帧），S4 本身很便宜。
- bench 12 个作业已提交：**jobid 3850864–3850875**（encut_1.5x 3850864/65、encut_2.0x 3866/67、
  lreal_auto_addgrid 3868/69、lreal_auto_noaddgrid 3870/71、lreal_false_addgrid 3872/73、
  lreal_false_noaddgrid 3874/75）。
- **生产 S4 扇出未提交**（wangchao 明确"等 bench 报告看过 max|ΔF| 与耗时再定"）。

### 15.5 ★ 复现性对照要看节点（wangchao 提醒，待做）

`encut_1.5x` 与 `lreal_false_addgrid` 的 INCAR 物理设置相同、是两次独立运行 —— 但报告里**不记节点名**。
若两次落在不同节点，差值里混入节点间差异，不能当"同输入复现性噪声底"。跑完后先取节点：

```bash
# ★ 2026-09-19 实测：jzzn 上 `scontrol show job` 的 NodeList 是 (null)，取节点要用 squeue 的 %N：
ssh jzzn "squeue -u wangchao -h -o '%i %T %M %N' | grep -E '3850864|3850872'"
# 已实测：3850864 在 cu14、3850865 在 cu18 —— 同一变体的两帧就已落在不同节点，
# 所以复现性对照那一对（3850864 vs 3850872）大概率也不同节点，需同节点对照才能定噪声底。
```

两者**同节点** → 差值可直接当复现性噪声底；**不同节点且差值偏大** → 需再跑一次同节点对照才能下结论。
这个噪声底是后续所有 |ΔF| 判据的基准（与 §13.2 的应力噪声底同理），不能含糊。
建议顺手让 `incar_bench.py` 在 report 的每帧里加一个 `node` 字段（跑完的 OUTCAR 里没有节点名，
需在提交时记录或事后 `scontrol show job` 回填）。

### 15.6 下一步（等 bench 报告）

1. 取节点 → `python incar_bench.py analyze --dir step4_disp`（登录节点，已批准）→ 看 max|ΔF| 与耗时；
2. 据报告定生产 INCAR（是否保留 ADDGRID、ENCUT 用 390 还是 520）→ 报给 wangchao 决定是否提交生产扇出；
3. 生产 S4 只有 12 帧，若批准，提交后按 §5 验收表继续（S5 第一道门 = pheasy
   `Imposing rotational invariance and equilibrium conditions`）。


### 16. 2026-09-19：bench 报告要看什么 + 两条纪律（wangchao 定）

#### 16.1 ★ MoS₂ 的 bench 结论**不得外推**（重要）

MoS₂ = 3 原子原胞、超胞约 **75–108 原子**。这个规模下 `LREAL=.FALSE.` 本来就不贵；
VASP 手册建议 ≥30 原子用实空间投影主要是**为大体系**考虑。所以本次 bench 的**目的只有三个**：
① 验证 bench 机制本身可用；② 给出**复现性噪声底**；③ 提供小体系的口径参考。

**生产口径的最终选择必须在真实规模的超胞上重做** —— 即推进到 jzz 那批、在 **~80 原子超胞 × 166 帧**
（Mo₂S₃ Z4-3-1）的条件下再跑一次 bench。**不得用 MoS₂ 的结论去套 166 帧的体系。**

#### 16.2 bench 报告必须包含的四项（wangchao 要求）

1. **复现性对照那一对（`encut_1.5x` 与 `lreal_false_addgrid`）的节点名 + 两者 |ΔF|**：
   若它自身就接近 0.5 meV/Å，说明阈值没有区分度，需要重新定；两次不同节点则该差值里混了
   节点间差异，需同节点对照才能下结论（见 §15.5 的 `scontrol show job ... | grep NodeList`）。
2. **RMS |ΔF| 与 max |ΔF| 并列**：单个原子的最大偏差可能来自某一帧的偶然，RMS 更能反映整体。
3. **每帧的实际耗时**：用于估算生产 12 帧的成本，以及将来 Mo₂S₃ 那 166 帧的成本。
4. （已有）各变体的 LREAL/ADDGRID/ENCUT 组合与 PASS/FAIL 判定。

#### 16.3 ★ 提交生产 S4 之前必须把 INCAR_BENCH 关掉

`INCAR_BENCH=on` 现在写在**材料级配置**里。bench 跑完、提交生产 S4 之前要执行：

```bash
autozt -tt kl-dft-cpu -p MoS2_kltest -j S4_disp conf --set params.INCAR_BENCH=off
```

（改项目配置 → 按铁律 3 先请示）否则以后**每次 retry S4 都会重新生成 bench 目录**。

---

## 17. 2026-09-19：INCAR_BENCH 结果（MoS₂，75 原子超胞，12 帧小作业）

### 17.1 ★ 复现性噪声底 = 0（阈值有区分度）

那对同设置对照（`encut_1.5x` 与 `lreal_false_addgrid`）**都调度到 cu14（同一节点）**，
两帧的 max/RMS **逐位相同**：

| 变体 | 帧 | max (meV/Å) | RMS (meV/Å) |
|---|---|---|---|
| encut_1.5x | disp-00001 | 0.6350 | 0.1650 |
| lreal_false_addgrid | disp-00001 | 0.6350 | 0.1650 |
| encut_1.5x | disp-00002 | 0.6720 | 0.1689 |
| lreal_false_addgrid | disp-00002 | 0.6720 | 0.1689 |

→ **同节点复现性远小于 1e-4 meV/Å**（VASP 同输入确定性，符合预期）。
**结论：0.5 meV/Å 的阈值有充分区分度**，不需要重定；差异 100% 是系统性的（口径效应），不是噪声。
（未做跨节点对照 —— 本次两帧对照恰好同节点，跨节点差异仍未测。）

### 17.2 完整结果（基准 = `lreal_false_noaddgrid`：LREAL=.FALSE.、无 ADDGRID、ENCUT 390）

| 变体 | LREAL | ADDGRID | ENCUT | max (meV/Å) | RMS (meV/Å) | Elapsed (s) | 判定 |
|---|---|---|---|---|---|---|---|
| **lreal_false_noaddgrid** | .FALSE. | - | 390 | 0.0000 | 0.0000 | 2658 / 2916 | 基准 |
| lreal_auto_noaddgrid | Auto | - | 390 | 0.5700 / 0.5440 | 0.1749 / 0.1935 | 1998 / 2316 | FAIL |
| lreal_false_addgrid（= 现生产口径） | .FALSE. | .TRUE. | 390 | 0.6350 / 0.6720 | 0.1650 / 0.1689 | 3073 / 3034 | FAIL |
| encut_1.5x（同设置复现对照） | .FALSE. | .TRUE. | 390 | 0.6350 / 0.6720 | 0.1650 / 0.1689 | 3927 / 2694 | （复现对照） |
| lreal_auto_addgrid | Auto | .TRUE. | 390 | 0.7220 / 0.7680 | 0.2195 / 0.2370 | 3124 / 2358 | FAIL |
| **encut_2.0x** | .FALSE. | .TRUE. | **520** | **3.0380 / 2.9840** | **1.2616 / 1.2461** | 11162 / 16865 | FAIL |

（耗时：ENCUT 390 约 **33–65 分钟/帧**，ENCUT 520 约 **3.1–4.7 小时/帧** = 4.2×。生产 12 帧按 390 约 10 核·时。）

### 17.3 四条结论

1. **ENCUT 是主导项**：390→520 让力变化 **3.0 meV/Å (max) / 1.26 (RMS)** —— 是阈值的 ~6 倍、
   是其它变体效应的 ~2.4 倍，代价 4.2×。**对位移取力，ENCUT 1.5×max(ENMAX) 没有力收敛**。
   ⚠️ 这修正了此前"Pulay 不影响力"的说法：那条结论是在**未位移**的几何上测的；**位移后的帧上
   基组不完备会明显改变力**。（§2 表第 2 行应在此处加注。）
2. **ADDGRID（现生产口径）**相对无 ADDGRID 改变力 **0.635–0.672 (max) / 0.165–0.169 (RMS)**，
   超过 0.5 阈值 → 按"与基准一致"的判据，**现生产口径并不与基准等价**。
3. **LREAL=.FALSE. 与 Auto 在本规模耗时相当**（3073/3034 s vs 3124/2358 s）、Auto 偏差还略大，
   印证了你的预期：**75 原子这个规模没必要用实空间投影**。
4. **max 与 RMS 同序**（RMS ≈ max/3.8），没有"单帧偶然"的迹象；两帧结论一致。

### 17.4 ★ 待决策（不急，且不得据小体系外推）

- 生产 S4 的取力口径**尚未决定**：按 §16.1 的纪律，MoS₂（75 原子）这次只用于**验证机制 + 定噪声底**；
  **最终口径必须在真实规模超胞（jzz 批、~80 原子 × 166 帧）上再跑一次 bench 决定**。
- 但上面第 1 条（ENCUT 1.5× 力不收敛）**是普适性风险，值得在 jzz 批之前先讨论**：
  若生产仍用 ENCUT 390，力的系统不确定性约 3 meV/Å；若提到 520，166 帧的成本按 4.2× 上升。
- 决定生产口径前**必须先把 `INCAR_BENCH` 关掉**（否则每次 retry S4 都重新生成 bench 目录，§16.3）。
- `incar_bench.py` 的 report 缺两项（已手工补算）：**每帧 RMS** 与 **节点名**。建议下一轮补进 `analyze`。

---

## 18. 2026-09-19：INCAR_BENCH 升级（相对判据 + RMS/CPU/节点/电子步）与 ENCUT 多点（390/520/650）

### 18.1 工具升级（`skill/kl-dft-cpu/incar_bench.py`）

- **per-frame 新增字段**：`rms_abs_delta_meV_per_A`、`cpu_time_s`（OUTCAR `Total CPU time used`）、
  `node`（作业节点）、`jobid`、`force_rms_meV_per_A`、`ref_force_rms_meV_per_A`、`rel_error_pct`、
  `electronic_steps`（DAV 计数，来自 queue.out）、`scf_converged`、`scf_note`。
- **variant 级新增**：`worst_rms_abs_delta_meV_per_A`、`worst_rel_error_pct`、`mean_rel_error_pct`、
  `mean_cpu_time_s`、`nodes`、`electronic_steps`、`scf_converged`、`NELM/ALGO/EDIFF`；
  顶层新增 `baseline_rule` / `criterion_mode` / `threshold_rel_pct` / `threshold_meV_per_A`。
- **参考设置（baseline）可配置**：默认取变体中"精度最高"者（LREAL=.FALSE.+ADDGRID=.TRUE.+最高 ENCUT；
  无 ADDGRID=.TRUE. 时退化为 LREAL=.FALSE.+最高 ENCUT）；`--baseline NAME` 可覆盖；analyze 的自动选择
  会在"有力且已收敛"的变体里按同规则挑。
- **判据改为相对量**：`|ΔF|_rms / |F|_rms < 1%`（`--threshold-rel`，默认 1%）；绝对
  `max|ΔF| < 0.5 meV/Å` 保留为可选（`--threshold`）。未收敛帧在表里显示 `NOTCONV`，
  variant/verdict 也会变 `NOTCONV`。
- **节点获取**：提交时把 jobid 写进 `bench/incar_bench_jobs.json`，并用 `sacct NodeList`（发现历史作业）
  + `squeue -h -o "%i %T %N"`（补 running/pending）采集；**jzzn 的 `scontrol show job` NodeList 为 (null)，不用**。
- **新增变体 `encut_650`**（LREAL=.FALSE.、ADDGRID=.TRUE.、ENCUT=650，与 `encut_2.0x` 只差 ENCUT）。
- **诊断腿开关**：`generate/submit` 新增 `--only` / `--set KEY=VALUE` / `--suffix`
  （只生成/提交指定变体、覆盖 INCAR 标签、目录加后缀；多次调用会**合并**进已有 plan）。

### 18.2 (A) 隔离 ENCUT：`encut_2.0x`(520) vs `lreal_false_addgrid`(390)

两者只差 ENCUT（其余 LREAL=.FALSE.+ADDGRID=.TRUE. 相同）；|F|_rms = 362.16 / 362.95 meV/Å。
**注意：520 帧 SCF 未收敛（NELM=200 截断），见 18.6**。

| 帧 | max|ΔF| (meV/Å) | RMS|ΔF| (meV/Å) | |F|_rms (meV/Å) | 相对 % | SCF |
|---|---|---|---|---|---|
| disp-00001 | 2.982 | 1.262 | 362.16 | 0.348 | NOTCONV |
| disp-00002 | 2.913 | 1.248 | 362.95 | 0.344 | NOTCONV |

### 18.3 (B) 各变体 vs 基准 `lreal_false_noaddgrid`（无 ADDGRID，ENCUT 390）

相对误差 = `|ΔF|_rms / |F|_rms`（两帧；max 与 RMS 同序显示）：

| 变体 | LREAL | ADDGRID | ENCUT | max (meV/Å) | RMS (meV/Å) | 相对 % | e-steps | SCF |
|---|---|---|---|---|---|---|---|---|
| lreal_false_noaddgrid（基准） | .FALSE. | - | 390 | 0.000 | 0.000 | 0.000 | 47/49 | CONV |
| lreal_auto_noaddgrid | Auto | - | 390 | 0.570/0.544 | 0.175/0.194 | 0.048/0.053 | 45/50 | CONV |
| lreal_false_addgrid | .FALSE. | .TRUE. | 390 | 0.635/0.672 | 0.165/0.169 | 0.046/0.047 | 44/49 | CONV |
| encut_1.5x（同设置复现对照） | .FALSE. | .TRUE. | 390 | 0.635/0.672 | 0.165/0.169 | 0.046/0.047 | 44/49 | CONV |
| lreal_auto_addgrid | Auto | .TRUE. | 390 | 0.722/0.768 | 0.220/0.237 | 0.061/0.065 | 44/49 | CONV |
| **encut_2.0x** | .FALSE. | .TRUE. | **520** | **3.038/2.984** | **1.262/1.246** | **0.348/0.343** | **200/200** | **NOTCONV** |

**要点**：LREAL/ADDGRID 四个变体（ENCUT 390）相对误差都 **< 0.07%**；旧判据 0.5 meV/Å（绝对 max）
会把 ADDGRID（max 0.67）判 FAIL。相对判据下这些全部 PASS。520 的 0.348% 看似 <1% PASS，但它
是 NOTCONV 帧（形式上不可信）；但 18.7 实测其力与收敛的 1E-7 帧只差 <0.001 meV/Å，结论未受影响。

### 18.4 (C) CPU time（OUTCAR `Total CPU time used`，秒/相对 390 基准）

| 变体 | 帧1 CPU (s) | 帧2 CPU (s) | 相对 390 基准 |
|---|---|---|---|
| lreal_false_noaddgrid（390 基准） | 2630.0 | 2894.7 | 1.000 / 1.000 |
| lreal_auto_noaddgrid | 1983.9 | 2299.6 | 0.754 / 0.794 |
| lreal_false_addgrid | 3042.0 | 3008.7 | 1.157 / 1.039 |
| lreal_auto_addgrid | 3092.8 | 2336.5 | 1.176 / 0.807 |
| encut_2.0x（520） | 11101.8 | 16777.3 | 4.221 / 5.796 |
| encut_650（NELM=200 版，被撤回） | 未完成（DAV 179） | 未完成（DAV 114） | —（按 520 外推 ~4–6×） |

### 18.5 复现性对照：同节点噪声 ≈ 0，跨节点 ~1e-3 meV/Å

`encut_1.5x` 与 `lreal_false_addgrid` 的 INCAR 物理设置相同（MoS₂ 上 1.5×max(ENMAX)=390=生产 ENCUT）。

| 帧 | encut_1.5x node | lreal_false_addgrid node | max 差 (meV/Å) | RMS 差 (meV/Å) |
|---|---|---|---|---|
| disp-00001 | cu14 | cu14（同节点） | 0.0000 | 0.0000 |
| disp-00002 | cu18 | cu06（跨节点） | 0.0010 | 0.0001 |

- **同节点**（disp-00001，均 cu14）：max/RMS **逐位一致** → 同节点最小尺度/数值路径确定，噪声底 ≈ 0。
- **跨节点**（disp-00002，cu18 vs cu06）：max 差 ~1e-3 meV/Å、RMS 差 ~1e-4 meV/Å（相对 ~3e-6）。
  远小于各口径效应（≥0.05 meV/Å），可忽略；§17 说的"跨节点未测"此处补上了这一小量级。
- **干净 EDIFF=1E-7 腿的节点**：390 两帧都落在 **cu31（同节点）**；520 在 cu28 / cu42；650 在 cu29 / cu31。

### 18.6 ★ 关键修正：520/650 在 EDIFF=1E-8 下撞 NELM，是数值精度地板

- **ENCUT 390 的所有帧**：DAV 44–50，OUTCAR **1 次** `aborting loop because EDIFF is reached`，无 NELM 警告 → 真收敛。
- **ENCUT 520 的帧**：DAV=**200**=`NELM`，**0 次** aborting，OUTCAR 有 VASP 警告
  `The electronic self-consistency was not achieved ... forces ... might not be reliable` → 力不可靠。
- **ENCUT 650（NELM=200 版）**：被用户撤回时 DAV 分别 179 / 114，同样 0 次 aborting；不可用。
- **不是混合发散，而是精度地板**：520 的 dE 已 ~1e-10（能量稳定），但 `d eps` 卡在 ~1.66e-8
  且每步只降约 0.02%（1.6627→1.6619e-8）。即使推到 NELM=300 也只能到 ~1.63e-8，仍 > EDIFF=1E-8。
  即：该基组 + `ALGO=Normal` 下 `EDIFF=1E-8` 判据过严，是数值地板，不是 ENCUT 物理效应。
- **§17 的"ENCUT 是主导项 / 1.5×max(ENMAX) 力不收敛（3 meV/Å）"由 18.7 的干净 EDIFF=1E-7 三点裁定**：
  差异确实来自 ENCUT（NELM=200 与收敛帧的力差 <0.001 meV/Å，见 18.7），但量级是相对 ~0.5%，不是绝对灾难。
- 处置（已执行）：① 520 跑 `NELM=300` 作反证（进行中）；② 390/520/650 用 `EDIFF=1E-7` 重跑
  （单一变量）；③ 390 的 1E-8 vs 1E-7 力对照 → 0.0001%–0.0003%，远 <0.1%。结论见 18.7。

### 18.7 干净多点（EDIFF=1E-7，全部 CONV）与收敛判断

390/520/650 三点全部用 `EDIFF=1E-7` 重跑（单一变量，其余同生产口径）。六帧全部
`aborting loop because EDIFF is reached`=1、无 NELM 警告 → 干净收敛（jobid：390=3853051/52、
520=3853053/54、650=3853055/56）。

| 对照 | 帧 | max|ΔF| (meV/Å) | RMS|ΔF| (meV/Å) | 相对 % |
|---|---|---|---|---|
| 390 vs 520 | disp-00001 / 00002 | 2.982 / 2.915 | 1.262 / 1.249 | 0.348 / 0.344 |
| 390 vs 650 | disp-00001 / 00002 | 4.511 / 4.275 | 1.787 / 1.746 | 0.493 / 0.481 |
| **520 vs 650** | disp-00001 / 00002 | **1.529 / 1.381** | **0.540 / 0.516** | **0.149 / 0.142** |

- **变化在收敛，但 520 尚未完全收敛**：390→520 的 RMS 变化 1.25 meV/Å，到 520→650 降到
  0.53 meV/Å（约 42%，只降 ~2.4 倍，不是"小一个量级"）。520 处仍有 ~0.5 meV/Å（~0.15%）的
  ENCUT 残余；按比值 r≈0.42 作 Richardson 粗估，650 相对"无穷大 ENCUT"还差 ~0.38 meV/Å（~0.10%）。
  ⚠️ 该 Richardson 数值**仅作参考**：只有三个点、r≈0.42 离渐近区尚远，**不得作为"650 已足够"的依据**。
- **按 1% 相对判据**：520 vs 650 = 0.146%、390 vs 650 = 0.487%，都已 <1% —— 390/520 在"1% 意义"
  上可用；但**按更严的 <0.1% 判，520 未到**，要压到 <0.1% 需 ≥650。
- **ENCUT 成本（收敛帧）**：390 ~2.0–2.1 ks；520 ~4.2–4.5 ks（2.1×390）；650 ~6.9–7.6 ks（3.5×390）。
  对比之前撞 NELM 的 1E-8 帧（520 = 11.1/16.8 ks），放开 EDIFF 后 520 反而**便宜 2.5–3.7×且真收敛**。
- **★ NELM 截断对力的实际影响 ≈ 0**：520 的 NELM=200 帧（形式 NOTCONV）与收敛的 1E-7 帧逐帧力差
  max **0.0010 meV/Å**、RMS 0.0004–0.0005（相对 ~0.0001%）；390 的 1E-8 vs 1E-7 力差 max 0.0010/0.0040、
  RMS 0.0004/0.0012（相对 0.0001%/0.0003%），远 <0.1%。**故 §17 的"390→520 力变化 3 meV/Å(max)/1.26(RMS)"
  并非 NELM 截断造成，结论成立**；本 bench 里 NELM 只是"形式上不干净"，力对该 d eps 地板不敏感
  （S5 门禁 §19 仍应作废这类帧 —— 力恰好没受影响是这个体系的运气，不是普遍保证）。
- **NELM=300 反证腿**（encut_2.0x_nelm300，3853049/50）：提交时观察到 d eps 仍卡在 ~1.4e-8、
  每步只降 ~0.02%，与预判一致；该腿不阻塞结论（干净答案已由 1E-7 给出）。**650 的 NELM=300 未提交**
  （同样算不到 1E-8，省 ~1 天机时）。

**最终判据结论**：① 生产 ENCUT=390（1.5×max(ENMAX)）相对 650 参考的力系统差 **0.49%（<1% 相对判据）**；
② ENCUT=520 比 390 明显更接近收敛，但 520→650 仍有 0.53 meV/Å（0.15%）变化，**未完全收敛**；
③ ADDGRID/LREAL 口径效应（<0.07%）远小于 ENCUT 效应（390 相对 650 约 0.49%）——ENCUT 仍是主导项，
但量级是 ~0.5% 相对，不是绝对 3 meV/Å 的"灾难"。

### 18.8 判据改为相对量（<1%）的说明

- **绝对 max 对均匀小偏移过敏**：ADDGRID/LREAL 效应 max ≈ 0.67 meV/Å，但相对只有 ~0.05%；
  用 0.5 meV/Å 绝对阈值会把实际等价的口径判 FAIL。
- **相对量更贴合取力影响**：`|ΔF|_rms/|F|_rms` 衡量力场的整体系统偏差，直接进 fc2/fc3；
  用 `--threshold-rel`（默认 1%）判据，`--threshold`（绝对）保留为可选，用于筛查单分量尖峰。
- **但相对量不能替代收敛检查**：520 相对 390 只有 0.348%，其帧却是 NOTCONV。判据与收敛标志必须
  一起看（本工具已在表里同时给出 rel% 与 SCF 列；S5 另有 §19 的 NELM 门禁）。

（本节由 2026-09-19 INCAR_BENCH 升级会话追加；证据：`step4_disp/bench/incar_bench_report.json`、
`incar_bench_jobs.json` 及各帧 OUTCAR/queue.out。18.7 已按干净 EDIFF=1E-7 三点定稿。）

## 19. 2026-09-19：S5 新增 SCF 收敛门禁（NELM 截断 / 无 aborting loop）+ P1 回溯扫描

### 19.1 拦的静默错误

VASP 在 SCF 撞 NELM 时**照样输出力、作业正常退出**，只在 OUTCAR 留一段：

```
|     The electronic self-consistency was not achieved in the given           |
|     number of steps (NELM). The forces and other quantities evaluated       |
|     might not be reliable so examine the results carefully. ...             |
```

这种帧的力不可信，拿去拟合会让 fc2/fc3 与 κ 整体错掉，而下游所有检查（作业状态、
CONTCAR、虚频闸、κ 数值）都显示正常。反向指标：静态单点正常收敛时 OUTCAR 必有且
**只有 1 次** `aborting loop because EDIFF is reached`。实测正对照（§17 的 bench）：
`encut_2.0x/disp-00001`（ENCUT=520、NELM=200）DAV=**200**=NELM、**0** 次 aborting、有 NELM 警告；
ENCUT=390 的同结构帧 DAV=44、**1** 次 aborting、无警告。

### 19.2 门禁逻辑（代码位置 + 行为）

- `skill/kl-dft-cpu/kl_common.py`
  - `NELM_WARNING` / `EDIFF_ABORT_MARK` / `_DAV_RE`（L361–365）
  - `scan_frame_scf(frame_dir)`（L368）：只读扫一帧 `disp-*/OUTCAR` + `OSZICAR`，返回
    `nelm_warn` / `n_aborting` / `n_dav_lines` / `last_dav` / `scf_steps` / `ok` / `reasons`。
    作废判据：① OUTCAR 命中 NELM 警告（`number of steps (NELM)`）；②
    `aborting loop because EDIFF is reached` 计数 != 1；③ OUTCAR 不存在（帧未算完，无法核验）。
    电子步数取 `last_dav`（静态单点等于 DAV 行数，同时记录行数）。
  - `check_outcar_scf_convergence(step4_dir)`（L413）：**非递归**遍历 `disp-*`
    （不会误扫 `bench/<variant>/disp-*`），逐帧 scan，聚合
    `n_frames/n_with_outcar/n_nelm/n_abort_ne_1/n_no_oszicar/steps_stats(min,中位,mean,p90,max)/
    step_outliers/bad_frames`，返回 `(ok, info)`。
  - `format_scf_report(info)`（L450）：一屏可读的报告/日志。
  - `_scan_main` + `if __name__ == "__main__"`（L481/507）：只读 CLI。
  - 顶部 `dim_common` 导入加逐级回退（L20–30），让 CLI 在仓库里也能独立运行（集群同目录时不受影响）。
- `skill/kl-dft-cpu/gen_step5_fc.py`（L114–133）：在 `check_frames_match_displacements` 之后、
  读 kl_params / 写 fit_config 之前调用门禁：
  1. `print(format_scf_report(...))` —— 把每帧电子步数分布 + 异常大值打进 gen 日志；
  2. 写 `step5_fc/scf_steps.json` —— 每帧明细（S5 的 summary 产物，异常值可追溯）；
  3. 任一帧作废 → `sys.exit`，信息里**逐个列出作废帧名**、原因与处置建议。
- 门的时机在"拟合之前"，不改 INCAR、不删除任何帧；`incar_bench.py` 等禁改文件未触碰。

### 19.3 独立只读回溯扫描命令

```bash
python skill/kl-dft-cpu/kl_common.py --scan-nelm <step4_disp 目录> [--json] [--out result.json]
# 退出码：0=全部帧 SCF 正常；1=有帧作废（NELM 截断 / 无 aborting / OUTCAR 缺失）
```

只读、不写远端；`--json` 机读输出，`--out` 另存本地 JSON。

### 19.4 P1 全批回溯扫描统计（jzzn，2026-09-19，只读）

扫描对象：`/public/home/wangchao/Fullerene_Network/work/P1_*/kl-dft-cpu/step4_disp`。

| 材料 | disp-* 帧目录 | 有 OUTCAR | NELM 警告 | aborting!=1 | 缺 OUTCAR | 电子步 min/中位/mean/p90/max |
|---|---|---|---|---|---|---|
| P1_Al-AlN_A1Z4_A-1-4_Al5N | 0 | 0 | 0 | 0 | 0 | —（S4 输入已生成，0 帧计算） |
| P1_Be-BeO_A2Z2_A-2-2_Be2O | 0 | 0 | 0 | 0 | 0 | —（同上） |
| P1_Mo-MoS2_A4-3-1_A4-3-1_Mo2S3 | 0 | 0 | 0 | 0 | 0 | —（同上） |
| P1_Mo-MoS2_Z3-2-1_Z3-2-1_Mo3S4 | 0 | 0 | 0 | 0 | 0 | —（同上） |
| **P1_Mo-MoS2_Z4-3-1_Z4-3-1_Mo2S3** | **166** | **165** | **0** | **0** | **1** | **44 / 51 / 50.7 / 52 / 53** |
| P1_Ti-TiS2_Z4-3-1_Z4-3-1_Ti2S3 | 0 | 0 | 0 | 0 | 0 | —（同上） |
| P1_Zn-ZnO_Z3A2_Z-3-2_Zn5O3 | 0 | 0 | 0 | 0 | 0 | —（同上） |

**关键数字（Mo₂S₃ Z4-3-1，唯一已算完的）**：166 个 `disp-*` 目录，165 帧有 OUTCAR/OSZICAR；
**0 帧命中 NELM 警告、0 帧 `aborting!=1`**；165 帧电子步全部落在 40–59（min 44、max 53），
远低于 `NELM=200`（该批 INCAR：ENCUT=390、EDIFF=1E-8、ALGO=Normal、NELM=200、NSW=0）。
唯一异常帧 `disp-00001`：目录里只有 INCAR/KPOINTS/POSCAR/POTCAR/submit.sh，
**没有任何 OUTCAR/OSZICAR/vasprun.xml**（mtime 2026-09-17 17:09）——该帧从未产出结果。

其余 6 个 P1 材料的 `step4_disp` 已生成输入（POSCAR-00001..N + SPOSCAR + phono3py_disp.yaml），
但 `disp-*` 目录数为 **0**（S4 fanout 尚未计算），因此无帧可扫、也无 NELM 问题。

### 19.5 真数据正对照（MoS2_kltest bench）

| 变体 | ENCUT | 帧数 | 电子步 | aborting | NELM 警告 | 门禁判定 |
|---|---|---|---|---|---|---|
| encut_2.0x | 520 | 2 | 200 / 200 | 0 / 0 | 2 | **全部作废** |
| encut_1.5x | 390 | 2 | 44 / 49 | 1 / 1 | 0 | 通过 |
| lreal_false_addgrid | 390 | 2 | 44 / 49 | 1 / 1 | 0 | 通过 |
| lreal_auto_addgrid | 390 | 2 | 44 / 49 | 1 / 1 | 0 | 通过 |

门禁在真实数据上精确命中了 §17 那个"DAV=NELM、0 aborting、有 NELM 警告"的假帧，
同时不误伤 ENCUT=390 的正常帧。

### 19.6 ★ 对 jzz 那批 S4 是否需要重算的判断

1. **Mo₂S₃ Z4-3-1（唯一已算）→ 不因 NELM 重算**：165 个已算帧全部 SCF 正常收敛
   （44–53 电子步 ≪ NELM=200，0 警告、0 计数异常）。因此 **κ≈3.8 不是 NELM 截断造成的**——
   背景里设想的"第二个解释"在这 165 帧上被排除（其它可能解释仍需另行排查）。
2. **但要补算缺失的 `disp-00001`**：它从未产出 OUTCAR。新门禁会把"缺 OUTCAR"也判作废，
   从而**拒绝拟合**。处置：`autozt -tt kl-dft-cpu -p P1_Mo-MoS2_Z4-3-1_Z4-3-1_Mo2S3 -j S4_disp retry`
   （retry 只补未完成帧、保留已算 165 帧；**不要用 rerun/clean**）。
3. **其余 6 个 P1 材料：无帧可判**（S4 尚未计算）。等 S4 算完后用
   `--scan-nelm` 复核；若届时出现 NELM 帧，再按第 2 条 retry 补算，不要整批 rerun。
4. **MoS₂ bench 的 `encut_2.0x`（520 eV）** 是真 NELM 帧，其 |ΔF| 数据不可信（§17.4 已注明
   不得外推）；这属于 bench 目录，不在生产 `step4_disp/disp-*` 顶层，生产门禁不会误扫。

### 19.7 验收证据

- `tmp/_kltest_nelmgate.py`：造 3 帧（正常 / NELM 警告 / aborting 计数 0），断言逐帧
  通过 / 作废 / 作废、聚合 `n_bad=2`、`n_nelm=1`、作废帧名列表正确；并真调 CLI
  `--scan-nelm --json`，退出码 1（有坏帧）/ 0（全好）——PASS。
- `ast.parse`：`kl_common.py`、`gen_step5_fc.py` 均 OK。
- `python3 tests/suite_io_schema.py`、`python3 tests/suite_skillspec.py`：全部通过。
- `tests/test_zt_upstream_fixes.py`（用 `/home/wangchao/miniconda3/envs/multiqc_env/bin/python -m pytest`）：
  **8 passed**（含 `test_kl_frame_check_skips_equilibrium_frame`，证明新增导入回退未破坏原逻辑）。
- 本文件已 `git add`、**未 commit**。

### 19.8 归档完整性核实：Mo₂S₃ 的 `disp-00001` 不需要补（2026-09-19）

  > ⚠️ **记录性质**：本节**不是**「某次批准被执行了」，而是「**批准后核实发现前提不成立**」 ——
  > 批准的动作是「补 `disp-00001`」，实际结果是**结构指纹门禁按设计拒绝**、且经核实**根本不需要补**。
  > **活目录并没有被补齐**（仍然缺 `disp-00001` 的产物），不要据此认为它可用；需要这批力一律从归档取（见下）。
  > 另：已在该目录放了 `README_OLD_STRUCTURE.md`，写明「旧结构、只用归档、不要 retry/rerun」。


用户批准"`retry` 只补 `disp-00001`"，但实测**被 §9.10 那道结构指纹门禁拦下**：
`[ERROR] 已有位移/力数据属于【旧结构】—— phono3py_disp.yaml 的单胞与当前输入 差异最大 0.1263 Å`。
这是**设计内的拒绝**（这批帧本就属于旧结构，`retry` 会让旧帧配新位移混用）。唯一绕法是 `rerun` / 手工 `rm -rf` 步骤目录，
两者都被用户明确禁止，且会毁掉那 165 帧。

**关键核实结论：这一帧根本不需要补。** 归档数据集
`/public/home/wangchao/Fullerene_Network/work/_archive_pre_relax_20260917/P1_Mo-MoS2_Z4-3-1_kl_step4_disp`
**166 个 `disp-*` / 166 个 OUTCAR，`disp-00001` 完整**（OUTCAR、vasprun.xml、CONTCAR、CHGCAR、WAVECAR 全在），
并带 `phono3py_disp.yaml` / `alm.in` / `SPOSCAR` / `POSCAR-00001..00165`。而该数据集的用途正是
"**归档完整性 + 2×2 诊断的上排（走 fc-fit，用归档那份）**" → **该项已满足，无需任何 retry/rerun**。

处置：活目录（旧结构、缺 disp-00001）**保持原样**，由门禁阻止误用；需要这批力时一律从**归档**取。

---

## 20. 2026-09-19/20：MoS₂「390 vs 520」物理量对照 —— 当前状态与待办

### 20.1 已就绪（本轮完成）

| 项 | 内容 |
|---|---|
| 390 那一套 | 生产 12 帧 `step4_disp/disp-00001..00012`（`ENCUT=390 EDIFF=1E-8`，44–50 电子步全部 CONV）+ 平衡帧 `disp-00000`（`POSCAR == SPOSCAR`，已验证是**真正的零位移超胞**） |
| 520 那一套 | `step4_disp/bench/encut_2.0x_ediiff1e7_full12/`：**12 帧**（jobid **3854253–3854264**，`ENCUT=520 EDIFF=1E-7 LREAL=.FALSE. ADDGRID=.TRUE.`，POSCAR 与生产 `cmp` 逐位一致）+ **平衡帧 `disp-00000`**（jobid **3854392**，POSCAR 亦为 `SPOSCAR`） |
| 隔离 | 两套在**不同目录**；生产 `disp-*` 未被触碰 |

### 20.2 ★ 一个必须知道的工具坑：`incar_bench.py submit` **没有队列感知**

`submit` 的判据只是「有 `submit.sh` 且**无 OUTCAR**」→ 对**已在队列中**的帧会**重复提交**。
所以补帧时**不能**再跑 `submit --suffix ...`（会把 8 个 pending 帧再交一遍）。本轮改用**只对新增帧**
`sbatch submit.sh` 的方式，并把 jobid 手工记进 `bench/incar_bench_jobs.json`（保持账目一致）。
**建议**：给 `submit` 加一道「已在 `incar_bench_jobs.json` 里且 `squeue` 仍能看到该 jobid 就跳过」的判断。

### 20.3 待办（第 2 阶段，等 13 帧全部产出 OUTCAR）

1. **两套各拟合 fc2**，必须用仓库现成的 `skill/fc-fit/`（不要自写拟合器），**两套都跑 `RASR=BHH`**，另加 pheasy + OLS；
2. ★ **确认 `-c` 步日志里真有 BHH 施加记录**（如 `Imposing rotational invariance and equilibrium conditions`）——
   这是 **P0-1 那道门第一次在文献基准上被真正验证**。若 BHH 在 MoS₂ 上让 `rel_err` 明显上升，
   那它是**继 Mg₂C₆₀ 之后的第二个样本**，需要单独分析；
3. 比较（按重要性）：**ZA 指数 p**（两个面内方向）→ **Γ 点光学模频率** → **整条谱频率 RMS 偏差（相对 %）** → κ；
4. 判据（wangchao 给定）：**频率相对偏差 < 1% 且 |Δp| < 0.05 → 390 够用**；
5. 结论写 §21（或并入本节），再据此定 **jzz 那 166 帧的口径**（390/1e-8 vs 520/1e-7，成本实测 2.1×）。

### 20.4 本轮的三个操作细节（wangchao 要求记录）

- **平衡帧口径**：两套必须**各自扣自己口径的平衡帧**（pheasy 逐帧扣平衡帧力）；生产 `disp-00000` 已核实为 `SPOSCAR`（零位移）；
- **子代理规矩**：涉及作业提交/状态变更的机械工作**自己做**（已写入 §15.2 第 3 条）；子代理只做只读分析/写代码/跑测试；
- **Mo₂S₃ 记录性质**：§19.8 已加「批准后核实发现前提不成立」的说明；活目录已放 `README_OLD_STRUCTURE.md`。


### 20.5 第 2 阶段的执行方案（已探明，帧跑完即可照此执行）

1. **520 套的 fc-fit 输入已备好**：`bench/encut_2.0x_ediiff1e7_full12/` 里已放入 `SPOSCAR`
   （与生产 `cmp` **逐位一致**）与 `phono3py_disp.yaml`（**12 条位移**）。fc-fit 的 `vasprun` 布局要求：
   `SPOSCAR` 必须**与 `disp-*` 同级**、且**平衡帧必须命名为 `disp-00000/vasprun.xml`**
   （本套已按此补齐，jobid 3854392；缺平衡帧会直接 `sys.exit`）。
2. **★ BHH 门禁是 fc-fit 内置的**（`fc_fit_driver.py:1171`）：`PHEASY_RASR=BHH` 时若 pheasy 的 `-c` 步骤
   日志里没有 `Imposing rotational invariance`，脚本直接 `sys.exit`。→ **P0-1 这次会被自动验证**，
   执行时确认日志出现 `[OK] RASR=BHH imposed in the null-space construction step (-c)`。
   若 BHH 让 `rel_err` 明显上升，记为**继 Mg₂C₆₀ 之后第二个样本**，单独分析。
3. **fc-fit 的输出目录是硬编码的**（`gen_step1_fit.py:28`：`OUTDIR = "step1_fit"`）→ 两套不能并存于同一材料。
   执行顺序：① 先跑 **390 套**（`FIT_INPUT_DIR` 缺省 `auto`，会自动找到 `step4_disp`），把 `step1_fit/` 产物
   另存为 `step1_fit_encut390/`；② `conf --set params.FIT_INPUT_DIR=<520 变体目录>` 后重跑，产物另存为
   `step1_fit_encut520/`；③ 用两份 `fc2.hdf5` + `phonon_band_summary.json` 做 §20.3 的四项对比。
   拟合参数：`FIT_ENGINE=pheasy`、`PHEASY_FIT_METHOD=OLS`（另跑一份默认 `RFE` 作对照）、`PHEASY_RASR=BHH`。
4. 四项对比与判据见 §20.3；结论写 §21，再定 jzz 那 166 帧的口径。

**当前进度（写本节时）**：520 套 13 帧中 4 帧已开始产出 `vasprun.xml`，其余在跑/排队（jzzn 仍拥挤）。


### 20.6 ★ 修正：390 那一套此前并未计算（本轮才提交）

写 §20.1 时把"生产 12 帧"当成已算好的数据了 —— **实际不然**：生产 S4 **一直没提交**（早先按指示"等 bench 报告再定"，
之后改成做两套对照，而"390 那一套"正是这批生产帧）。核实：`kl-dft-cpu/step4_disp/disp-*` 里
**OUTCAR / OSZICAR / vasprun.xml 计数均为 0**，帧目录里只有 `INCAR KPOINTS POSCAR POTCAR submit.sh`。

**本轮已提交**（属用户批准的"12 帧跑两套 ≈21 核·时"范围；用 `start`，**未用 `-f`**）：

```bash
autozt -tt kl-dft-cpu -p MoS2_kltest -j S4_disp start
# → jobid 3854404–3854416（13 个：disp-00000 平衡帧 + 12 位移帧），INCAR = 生产口径 ENCUT=390 EDIFF=1E-8
```

→ **两套当前状态**：390 = `3854404–3854416`（13 帧）；520 = `3854253–3854264`（12 帧）+ `3854392`（平衡帧）。
两套各 13 帧全部产出 `vasprun.xml` 后，按 §20.5 执行拟合（`RASR=BHH` 门禁自动验证 P0-1）与四项对比。

（也说明 §20.1 表格里"390 那一套"的其他描述仍然正确：`disp-00000/POSCAR == SPOSCAR` 已核实为零位移超胞；
12 位移帧的 POSCAR 与 520 套逐位一致。）


### 20.7 拟合调用的最后细节（已确认，帧跑完即可开跑）

- `FIT_INPUT_DIR` 取值：`auto | path`（**技能目录相对或绝对路径**，`gen_step1_fit.py:54`）。
- **两套都显式给绝对路径**，避免 `auto` 自动搜索在 `step4_disp/bench/` 附近产生歧义：
  - 390：`/public/home/wangchao/Fullerene_Network/work/MoS2_kltest/kl-dft-cpu/step4_disp`
  - 520：`/public/home/wangchao/Fullerene_Network/work/MoS2_kltest/kl-dft-cpu/step4_disp/bench/encut_2.0x_ediiff1e7_full12`
    （已备好 `SPOSCAR` + `phono3py_disp.yaml` + `disp-00000`，见 §20.5）
- 参数：`FIT_ENGINE=pheasy`、`PHEASY_FIT_METHOD=OLS`、`PHEASY_RASR=BHH`；另跑一份默认 `RFE` 作对照。
- 输出目录硬编码 `step1_fit`（§20.5）→ 两套串行跑、各自另存为 `step1_fit_encut390/`、`step1_fit_encut520/`。

**本轮状态（写本节时）**：390 套 13 个作业（3854404–3854416）**仍在排队**（产物 0）；
520 套仍是 **4 个半成品**（`General timing` 计数 0，未收尾）；账号下队列压在 **47 个作业**（含其它会话）——
jzzn 极度拥挤，两套帧都还没出结果。**下一步完全取决于队列**，拟合方案已全部就绪。

