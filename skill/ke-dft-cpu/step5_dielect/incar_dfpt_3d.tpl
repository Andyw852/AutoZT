# =====================================================================
# incar_dfpt_3d.tpl —— DFPT 介电常数（3D，IBRION=8 + LEPSILON）
# 一次微扰求 ε∞ 和 ε₀（离子+电子）。占位符：{{SYSTEM}} {{ENCUT}} {{GGA}}
# 【注意】KPAR/NCORE 默认取 1（保守、已被四材料 SOC DFPT 验证）。
#   官方文档对 DFPT 的 k 点并行历来持保留态度；实测（2026-09-14/15）：
#     KPAR=8 + ISYM=2 + SOC 的 IBRION=8 能跑出 "8 groups" 并产出张量，
#     但那次同时改了 ISYM/网格/NBANDS，**结果正确性未确认**，故不作默认。
#   并且——提高 KPAR **并不降低内存总量**（见文末内存模型），别指望它救 OOM。
# =====================================================================
SYSTEM = {{SYSTEM}}

ISTART = 0
ICHARG = 2
GGA    = {{GGA}}
{{VDW_LINE}}

PREC   = Accurate
ENCUT  = {{ENCUT}}
LREAL  = .FALSE.
LASPH  = .TRUE.
ADDGRID= .TRUE.        # patch_addgrid：细化增广网格，减小力常数/声学求和规则数值误差（amset 推荐）

ALGO   = Normal
EDIFF  = 1E-8
NELM   = 200
NELMIN = 6
# ISMEAR/SIGMA：窄隙体系要留神。SIGMA 相对带隙偏大时占据数会糊开，
#   而 LPEAD 的 PEAD 重叠矩阵与 DFPT 的占据流形都依赖锐利的占据划分。
#   gen_step8_dielect.py 会读 step2.15 的 discriminant.json，带隙 < 0.30 eV 时告警。
ISMEAR = 0
SIGMA  = 0.05

# ---- DFPT 微扰 ----
IBRION   = 8           # DFPT + 对称约化
LEPSILON = .TRUE.      # 静态介电张量（含离子贡献）
# LPEAD：.TRUE. = 用有限差分(PEAD)求 |∇k u_nk>；.FALSE. = 解 Sternheimer 方程（纯 DFPT）。
# 两者都合法，官方 SiC 介电教程用的就是 LEPSILON=.TRUE.+LPEAD=.TRUE.+IBRION=8，
# 并指出 .TRUE. 对 k 点收敛往往更快 —— 所以这里**不是**"路线接错"，别据此禁用。
# 但官方同时声明：**LPEAD 不支持金属**。PEAD 的重叠矩阵 S 定义在**占据流形**上，
# 占据/空态划分一模糊，矩阵就退化。窄隙 / 强 SOC / 近金属体系正落在这个风险里。
# 实测（A2B2Te5 四体系, PBE+SOC, 带隙 0.15~0.23 eV, ISYM=2, KPAR=8）：
#   LPEAD=.TRUE.  -> IBRION=8 的 DFPT 直接 SIGSEGV（rc=139），死在 MACROSCOPIC 张量之前；
#   LPEAD=.FALSE. -> 正常产出 ε∞（含局域场）与 Born 有效电荷、离子贡献。
# 故默认取保守的 .FALSE.。体系明确是宽隙绝缘体时，可以自行开回 .TRUE. 换取更快的 k 收敛。
LPEAD    = .FALSE.
NSW      = 1
ISYM     = 2

LWAVE  = .FALSE.
LCHARG = .FALSE.

# ★DFPT 不能 k 点并行
NCORE  = 1
KPAR   = 1

# ★ KSPACING 加密不是免费的（2026-09-14 实测）
#   DFPT 的响应数组内存随 k 点数增长。实测（Pb2Sb2Te5, 9 原子, SOC, 96 核）：
#     KSPACING=0.04 (12x12x2, NKPTS=148) -> 正常跑完（B/C/D 三格，约 6-8 h）
#     KSPACING=0.02 (24x24x4, NKPTS=?)   -> 在 DFPT 响应段一开始就被 OOM 杀：
#         "Linear response G [H, r] |phi>, progress: Direction: 1"
#         -> "process rank 42 ... exited on signal 9 (Killed)"，rc=137
#   ★★ 收敛判据的陷阱：只看 eps_static 会被"假收敛"骗过（2026-09-15 实测）
#     A2B2Te5 (PBE+SOC, ISYM=-1, 只改 k 网格，其余参数逐字节一致)：
#        网格        NKPTS   eps_inf_xx   eps_ion_xx   eps_static_xx
#        9x9x1         81      100.727      153.640       254.367   <- 生产值
#       12x12x1       144       88.298      173.193       261.491
#        9x9x2        162       91.797        ...           ...     <- 已收敛
#       12x12x2       288       91.674        ...           ...     <- 与 9x9x2 差 0.13%
#     要点：**eps_ion 同样强烈依赖 k 网格，且变化方向与 eps_inf 相反**，
#       两者大幅抵消 —— eps_static 只变 2.8%，而两个分量各自变了约 12%。
#       所以"eps_static 稳定"**不能**证明收敛。
#     后果：AMSET 的 Frohlich 耦合是 1/eps_inf - 1/eps_static，两份量不抵消：
#       实测该体系耦合 0.005996 -> 0.007501（**+25%**），POP 散射显著增强。
#       即 eps_inf 未收敛会让迁移率被高估。
#     做法：看收敛性必须**分别**核对 eps_inf 与 eps_ion，不能只看 eps_static。
#     本体系实测：kz=1 是主因（c=17.44 A 只取 1 个 k 点），把 kz 提到 2 即收敛，
#       面内网格不必加密 —— 见 gen_step8_dielect.py 的 _enforce_min_kz()。
#   ★ 内存需求模型（2026-09-15 修正；此前"提高 KPAR"的建议是错的，别再照做）
#     DFPT 总内存需求 ≈ NBANDS x NKPTS x phi，实测标定 phi ≈ 5.5 MB/(k点·能带)。
#     关键性质：**该总量与 rank 数、与 KPAR 都无关**。
#       KPAR=1 + 24 ranks：每 rank 持 (NBANDS/24) x NKPTS
#       KPAR=8 + 96 ranks：每 rank 持 (NBANDS/12) x (NKPTS/8) —— 与上式同值
#     多开节点或提高 KPAR 只是把同一块内存切成更多份，**不降低总量**；
#     真正降总内存只有两条路：减少 NKPTS（粗网格）或减少 NBANDS。
#     要跑密网格，必须让「聚合节点内存 >= 总需求」——即**跨多节点**跑。
#     实测标定（96 ranks, NKPTS=1156, 每 rank 1 个能带）：6.41 GB/rank
#       => 总 615 GB，单节点 768 GB 撞顶 OOM；同规模放 2 节点(~1.5 TB)可容纳。
#     本类体系参考值（NBANDS=96）：81 k点 ~43 GB；288 k点 ~153 GB；
#       768 k点 ~408 GB（单节点尚可）；2304 k点 ~1.2 TB（必须多节点）。
#   ★ SOC 体系的额外内存代价：ISYM=-1 关掉对称约化，NKPTS 按全 BZ 计。
#     例：9x9x1 在 ISYM=2 下约 12 个 k 点，ISYM=-1 下是 81 个，差 6~7 倍，
#     内存与耗时同比例放大。这是 SOC + DFPT 绕 spinor bug 必须付的代价。
#   申请整节点内存：--mem=0（SLURM 表示"本节点全部内存"）或 --exclusive。
#     hanhai25 的 DefMemPerNode=UNLIMITED —— 调度器不预留内存，超订只在运行时炸。
#   ★ 提交时限：**不必显式写 --time**。jzzn 实测复核（2026-09-15）：fc-fit 等技能
#     的模板同样没写 --time，仍拿到 Timelimit=2:00:00，长作业正常跑完；
#     A2B2Te5 的生产 SOC DFPT 作业（无 --time）也跑满 6-8 h。默认值够用。
#   ★ 真正要防的是另一件事：**#SBATCH 指令块必须完整有效**。
#     Slurm 只解析文件最前部、且在第一个非注释命令之前的 #SBATCH。
#     一旦渲染后的命令/占位符插在 #SBATCH 之前，**整块指令被静默忽略**：
#     作业名回落成脚本名（submit.sh）、分区与时限全部回落到默认值。
#     实测症状：sacct 里 JobName=submit.sh 且 Timelimit=00:01:00，长作业秒死。
#     （这正是 2026-09-12 cohp-cogito 作业连挂两次的原因：占位符位置放错。）
#     做法：模板里命令块一律放在所有 #SBATCH 之后；提交后核对作业名与 TIME_LIMIT。
#   反面教训：网格加密与"DOF x NFREE"是两个独立的乘数，不要一起放大。
#   在 128 原子体相上 IBRION=6 + ISIF=3 有 102 个自由度，NFREE=2 即 ~204 次 SCF，
#   这一项与 k 点无关却往往是真正的成本大头。
