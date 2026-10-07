# ke-dft-cpu / tools —— 介电重算与输运复核工具集

2026-09-16 由四材料介电重算工作沉淀而来，全部经实测验证。

## 为什么有这套工具

ke 技能原本只覆盖到「生成 settings.yaml 并提交 amset」。当需要
**更换介电（或任何输入）后重跑 AMSET 并复核结果**时，缺少可复用的抓手，
于是产生了这些脚本。它们解决的是流程层面的重复劳动，不是新的物理方法。

## 工具清单

| 脚本 | 用途 |
|---|---|
| extract_diel.py | 从 DFPT 的 OUTCAR 提取 eps_inf / eps_static / NKPTS |
| update_amset_settings.py | 更新 AMSET settings.yaml 的介电块（自适应紧凑/多行两种格式）|
| val_yaml.py | 校验 settings.yaml 的 3x3 形状与数值 |
| rep_transport.py | 输运结果分析：table / cmp / temp / pos 四模式 |
| make_patch_data.py | 由新旧 transport 生成报告补丁 JSON |
| patch_report.py | 改写 HTML 报告的表与正文（限定表内、唯一匹配）|
| regen_soc.py | 重生成 300K SOC/无SOC 对比 CSV |
| regen_soc_full.py | 重生成全温度全浓度对比 CSV |
| dfpt_status.py | DFPT 完成判据（eps_inf / ionic / born / finished 四信号）|
| apply_new_dielectric.sh | 编排：数值校验 -> 提取 -> 更新 -> 提交（四道闸门）|
| cmp_diel.py / cmp_eps.py | 新旧介电张量对比 |
| cmp_csv.py / cmp_keyed.py | CSV 逐值比对 / 按键对齐比对 |
| desym_fix_validate.py | V115/V116：用同一份全网格 h5 检验去对称化相位公式（原式/修正式）。实数据用**对照组判据**（坏操作点 vs 好操作点 = AMSET 原生路径）；`--shift t1 t2 t3` 平移原点造坏操作 + 查平移不变性（标准原点 GaN 必须用它）。V117 起自动从 vasprun 取全部带能量识别被带窗口切开的简并组（纤锌矿 kz=π/c 粘连）、边界标签换到 AMSET 约定、查源/目标本征值一致性、打印平面波截断。判据见 VERIFICATION V115 §6、V116、V117 |
| template_drift.py | V137：项目级模板副本（init 时复制进 project_setting/templates/ 的 *.tpl，find_asset 优先用它）与技能模板的漂移审计：逐键列出不同的 INCAR 键（★ 关键键），用 git 历史判断副本是技能模板的原样旧版 [STALE]（可直接删，让技能模板生效）还是有手改 [CUSTOM]；只读 |
| vasp_eta.py | V144：在跑的 VASP 作业还要多久、会不会撞墙时。按 OUTCAR 的 LOOP+ 数与耗时外推；IBRION=6 的总步数 = 1 + NFREE×3×不等价原子 + NFREE×6（ISIF≥3），IBRION=5 = 1 + NFREE×3N；墙时上限取 submit.sh 的 --time，否则按 QoS 查表（regular 24 h / premium 48 h）。预计超过上限 80% 退出码 1。大作业跑起来几个小时就该看一次 |
| mobility_vs_dpt.py | V143：**读 transport.json 迁移率一律用它，不要手写脚本**。每个掺杂一行：载流子（按 AMSET 约定，负掺杂 = 电子，唯一真源 ke_common.carrier_of_doping）、Seebeck、overall 与各机制迁移率（2D 面内平均 / 3D 迹/3）、同载流子的 DPT 与 ADP/DPT。Seebeck 符号与载流子不符标 ★（退出码 1）；ADP/DPT 超出 [1/3, 3] 标 ⚠ 并提示 dp_valley_probe；n/N_eff ≥ 0.5（简并，V146）或 DPT 带边简并（V145）、DPT 的 m* 网格分辨不出带边曲率（V147）的行只注明、不提示多谷 |
| dp_valley_probe.py | V141：带边附近形变势 D(k) 的分布按谷拆开（读 S7.1 的 h5 + undeformed 能带，AMSET 环境，秒级）。报每个谷的 ΔE、载流子占比、对 ADP-only 迁移率的贡献占比、rms D，以及倍数 E1²·⟨1/D²⟩（≈ 只因 D 的分布，ADP-only 迁移率相对单谷 DPT 的倍数）。用于判断 ADP-only ≫ DPT 是多谷物理还是 D 场出错 |
| lineage_check.py | V140：上游同源核对 —— 各步 POSCAR 是否都等于 S1 当前的 CONTCAR（S1 重新弛豫后没全部重跑时，哪些步骤还是旧结构；S2.3 HSE 会从还没重跑的 S2.2 拷旧结构），派生产物（S4 h5、S7.1 形变势、S2 画图的 band_summary、S5.1）是否旧于来源的计算输出；S8/S8.4 的 gen 做同一核对（STRUCTURE_GUARD）。V148 起另查 S8.1/S8.3 是否旧于 S8.2/S8/S8.4、S8/S8.4 settings.yaml 的带隙是否与 S2 画图一致（这两项 S8 的 gen 不查）。V175 起另查 k 网格：S3 不合现行规则（面内间距、2D kz）、S7 和 S3 不是同一张网格且不合规——已经算完的旧结果只有这里能查出来。`--invalidate-from <步骤>`：已经 rerun 过的上游，把它下游的完成标记改名归档，交给 auto-advance 重做（S2 画图按当前带隙决定是否归档 S8/S8.4）；其余只读 |
| elastic_ionic_decompose.py | V154：S6 的 TOTAL ELASTIC MODULI 不正定时先跑它再决定要不要重算。把离子弛豫项 −(1/Ω)ΛᵀK⁻¹Λ 按力常数本征模拆开，先复现 VASP 自己的离子项（自检，复现不了就退出码 2），再单列 3 个平移模的贡献：翻负主要来自平移模 = 数值假象（扣掉后的张量写进 elastic_ionic_decompose.json，不用重算）；扣掉仍不正定 = 真实的内应变耦合（列出软模），不可用。Mo2S3 的 C66：刚性离子 +267.80、离子项 −430.7、TOTAL −162.94。V155：VASP 6 的 IBRION=6 不打印 SECOND DERIVATIVES / 离子项块，力常数改从 S6 目录的 vasprun.xml（<dynmat> hessian）读，离子项用 TOTAL − SYMMETRIZED，内应变两套口径都试。V156：自检改为"VASP − 本工具 的差值是否落在严格均匀平移子空间"（VASP 怎么处理近零模无从精确复现），另用 OUTCAR 的振动频率核对力常数单位；差值落在软模方向则判不可用 |
| regate_projects.py | V115：用逐操作判据回查已有项目的 S8/S8.4 是否"IBZ + 真实重叠 + 有坏操作 + 没打补丁"（静默算错）；V117：全网格 h5 的 k 点标签不在 (-0.5, 0.5] 也判 ★ WRONG；V134：settings.yaml 的弹性张量不正定（含 ADP）判 ★ WRONG，2D 面内剪切 C66 不到 C11 的 2% 判 ⚠ ELASTIC（疑似 Voigt 顺序没重排）；V135：settings.yaml 与 step6_elastic/OUTCAR 逐位对照，确定是 VASP 顺序（没重排）判 ★ WRONG（3D 非立方也查得出），重排后仍不正定的提示先重做 S6；V136：transport.json 比它用到的上游产物（形变势 h5 / vasprun / h5 / S5、S6 的 OUTCAR）旧判 ⚠ STALE，S5.1 介电校验 ok=false 判 ★ WRONG |
| compare_deformation_h5.py | V123：重新 gen S7.1 前后两份形变势 h5 的**整场**对比（⟨\|D\|²⟩ 新/旧 ≈ ADP 迁移率粗估因子的倒数、逐点相对变化），给 ⚠ STALE-DP 项目分诊：变化小的可不重跑 S8/S8.4。只看带边 E1 判断不了（MoS₂ 带边 iso 几乎不变，ADP 迁移率却 ×2） |
| compare_transport_json.py | V125：两份 AMSET transport.json 逐 (掺杂, 温度) 对比（迁移率 overall/各机制、电导率、Seebeck；2D 取面内平均），并报告新旧各自的 xx/yy 各向异性。用于 IR_FIX / KZ_CAP_2D 的真实数据验证：阈值默认 2%，退出码 0 = PASS、1 = FAIL、2 = 网格对不上 |
| kz_overlap_probe.py | V126：2D 运行目录里用 AMSET 自己的重叠计算，量化沿 k_z **插值**波函数系数的假象：比较旧 21 层 / 新 3 层集合上的平均重叠，给出预测的 μ新/μ旧（秒级、不提交作业）。用来判断 KZ_CAP_2D 与旧结果的差别是截断引入的误差，还是旧结果里的插值假象 |
| star_audit.py | V130：同一张密网格上，旧跑法（AMSET 原式映射，逐成员算）与 IR_FIX 跑法的 mesh.h5 逐星比较散射率，把迁移率差拆成 Jensen（对率平均 vs 对 τ 平均）、仿真误差（星平均没模仿到的部分）两项，并给出只用代表点的预测；只读 h5、不需要 AMSET。用来判断某材料能否开 IR_FIX |
| fermi_probe.py | V131/V132：AMSET 求不到费米能级（"Could not calculate Fermi level position"）或下游 pinv 报 SVD 不收敛时，在运行目录（S8.4 或 S8；V133 在合成 3D 目录上验证过）按作业同一套插值重建 DOS（不算散射），给出带隙内假态、本征 E_F、逐个掺杂/温度的 AMSET **原搜索**结果与其实际掺杂偏差（> 1% 标"宽松!" = 旧结果在该点不可信）、贪心搜索停点与二分解（与 amset_fermi_fix 同一函数）；几分钟、不提交作业 |

注意：DFPT 产物的**数值校验器**不在这里，它是技能正式步骤的一部分，
见 step5_dielect/validate_dielectric.py（含 NaN/单位矩阵/声学求和规则/
双算法交叉/对角性/SOC-ISYM 六项判据）。

## 典型用法：换介电后重跑 AMSET

    bash tools/apply_new_dielectric.sh <材料> <新DFPT目录> [--submit]

四道闸门依次把关，任一不过即中止：
1. validate_dielectric.py 校验新 DFPT 产物（不完整/NaN 直接拒绝）
2. 提取值的有限性与对角 >= 1 检查
3. 更新后 yaml.safe_load 校验
4. 归档旧 transport.json 后再提交

## 环境依赖

- 需能 import yaml（amset_clean 环境满足）
- apply_new_dielectric.sh 里的 amset phonon-frequency 需要 amset_clean

## 相关：本次发现的 AMSET nworkers 问题（已在技能内修复）

AMSET 源码 amset/interpolation/bandstructure.py:182 把 nworkers=-1 解释成
multiprocessing.cpu_count()，即用满节点全部核。技能此前不写 nworkers，
于是超订 submit 模板申请的核数（jzzn 节点 192 核 vs 申请 24 核）。

现已在 step8_amset/gen_step10_amset.py 显式写入 nworkers（默认 24，
可由 step.conf 的 NWORKERS 覆盖），并在 step.conf 里写明必须与
提交模板的 --ntasks-per-node 保持一致。
