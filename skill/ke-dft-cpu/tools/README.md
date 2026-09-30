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
| regate_projects.py | V115：用逐操作判据回查已有项目的 S8/S8.4 是否"IBZ + 真实重叠 + 有坏操作 + 没打补丁"（静默算错）；V117：全网格 h5 的 k 点标签不在 (-0.5, 0.5] 也判 ★ WRONG |
| compare_deformation_h5.py | V123：重新 gen S7.1 前后两份形变势 h5 的**整场**对比（⟨\|D\|²⟩ 新/旧 ≈ ADP 迁移率粗估因子的倒数、逐点相对变化），给 ⚠ STALE-DP 项目分诊：变化小的可不重跑 S8/S8.4。只看带边 E1 判断不了（MoS₂ 带边 iso 几乎不变，ADP 迁移率却 ×2） |
| compare_transport_json.py | V125：两份 AMSET transport.json 逐 (掺杂, 温度) 对比（迁移率 overall/各机制、电导率、Seebeck；2D 取面内平均），并报告新旧各自的 xx/yy 各向异性。用于 IR_FIX / KZ_CAP_2D 的真实数据验证：阈值默认 2%，退出码 0 = PASS、1 = FAIL、2 = 网格对不上 |

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
