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

**可复用的一类判据（本轮两次起作用）**：
- **结构指纹**：把"当前结构"与"历史产物所依据的结构"做逐分量比对（S4 的 disp yaml 校验、种子指纹方案 §9.5）；
- **物理量后置校验**：不信标记、只信末态物理量（末态面内应力、末态最大力、`max|Δ晶格|`）；
- **同几何重复性**：同一几何两次独立计算（如段 b vs 段 c）的差值就是数值噪声，可用来判断"抖动"还是"真变化"；
- **闸门只放宽有证据的**：如 Pulay 只对变胞段加 ENCUT、力上限回到与旧判据一致的 0.01，都有实测依据。

**已实测的噪声量级（供判据设计参考）**：MoS₂（3 原子、15×15×1、520 eV）面内 0.0059 kB、
两遍完全一致；Mo₂S₃（20 原子、正确结构）0.0083 kB；而**畸变几何**上同一几何两次可差 0.065 kB
（正是 §13.2 那次的假象）。→ 阈值对结构正常的材料合适，0.05 下限保留。

