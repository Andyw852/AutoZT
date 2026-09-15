# AutoZT 相关工作补充与优势分析

本文献清单用于定位 AutoZT，不把“文献没有讨论”写成“文献不具备”。DREAMS 与 AtomAgents 的修正结论见 `RELATED-WORK.md`。

## 1. 相关文献谱系

| 方向 | 代表文献 | 系统贡献 | AutoZT 可区分的切入点 |
|---|---|---|---|
| 工作流与计算基础设施 | Pizzi et al., **AiiDA**, *Comput. Mater. Sci.* 111, 218 (2016), DOI: 10.1016/j.commatsci.2015.09.013 | provenance-aware 的材料计算工作流、数据节点和可追溯执行 | AutoZT 面向 VASP/MACE 热输运等固定技能，提供物理判据驱动的推进和运维恢复 |
| 工作流与作业编排 | Jain et al., **FireWorks**, *Concurrency Computat.: Pract. Exper.* 27, 5037 (2015), DOI: 10.1002/cpe.3505 | 大规模计算任务的队列、依赖和故障恢复 | AutoZT 将步骤状态、集群切换、失败分类和 retry 语义绑定到材料工作流 |
| 自动化材料工作流 | Ong et al., **atomate**, *Comput. Mater. Sci.* 68, 314 (2013), DOI: 10.1016/j.commatsci.2012.10.023 | 基于 pymatgen/FireWorks 的可复用材料计算 recipes | AutoZT 的重点是跨 DFT/MLIP、跨 CPU/GPU 和生产吞吐治理 |
| 新一代工作流 | Montoya et al., **atomate2**, *Comput. Mater. Sci.* 223, 112169 (2023), DOI: 10.1016/j.commatsci.2023.112169 | 更模块化、可组合的材料工作流构建 | AutoZT 增加面向实际超算差异的 HPC profile、审批网关和逐文件校验 |
| 材料数据与计算平台 | Jain et al., **Materials Project**, *APL Mater.* 1, 011002 (2013), DOI: 10.1063/1.4812323 | 大规模材料性质数据库与计算基础设施 | AutoZT 是执行层与运行监督层，不是材料数据库 |
| 原子模拟接口 | Larsen et al., **ASE**, *J. Phys.: Condens. Matter* 29, 273002 (2017), DOI: 10.1088/1361-648X/aa680e | 统一原子模拟 calculator 与结构操作接口 | AutoZT 位于 calculator 之上，负责多步任务、队列、判据和恢复 |
| 声子计算 | Togo and Tanaka, **Phonopy**, *Scr. Mater.* 108, 1 (2015), DOI: 10.1016/j.scriptamat.2015.07.003 | 有限位移声子、动力学矩阵与热力学性质 | AutoZT 编排 phonon/force/fit/transport 链路并检查虚频、残差和收敛 |
| 晶格热输运 | Li et al., **ShengBTE**, *Comput. Phys. Commun.* 185, 1747 (2014), DOI: 10.1016/j.cpc.2014.02.015 | 三阶力常数与玻尔兹曼输运方程求解 | AutoZT 将结构优化、力生成、力常数拟合和输运求解变为可监控流水线 |
| 晶格热输运 | Carrete et al., **almaBTE**, *Comput. Phys. Commun.* 185, 277 (2014), DOI: 10.1016/j.cpc.2013.09.007 | 声子输运与复杂晶体热导率计算 | AutoZT 关注重复批量计算的输入一致性、失败恢复和结果归档 |
| 机器学习势 | Batatia et al., **MACE**, arXiv:2206.07697 (2022); 后续版本发表于 *NeurIPS* 2022 workshop | 等变消息传递原子模型与高效势能面 | AutoZT 将 MACE 与 VASP 放进同一技能/材料/HPC 状态模型，并支持 MLFF 迭代 |
| AI 科学智能体 | Bran et al., **ChemCrow**, *Nature Machine Intelligence* 8, 525 (2024), DOI: 10.1038/s42256-024-00832-8 | LLM 规划器调用化学工具完成开放式任务 | AutoZT 的稳态推进由物理判据和状态机完成，LLM 可作为失败路径的可替换上层 |
| AI 材料发现 | Zeni et al., **MatterGen**, *Nature* 626, 123 (2025), DOI: 10.1038/s41586-025-08628-5 | 生成式模型进行无机材料结构设计 | AutoZT 的目标是已知计算流程的可靠生产吞吐，而非生成新结构 |
| 多智能体材料模拟 | Wang et al., **DREAMS**, arXiv:2507.14267v2 (2025) | DFT agent、HPC agent、安全 guard 与 claim provenance | AutoZT 将保证附着于执行底座：判据、动作闸门、产物哈希、多集群和可重放会话 |
| 多模态材料发现 | Ghafarollahi and Buehler, **Automating alloy design and discovery with physics-aware multimodal multiagent AI**, PNAS 2025, DOI: 10.1073/pnas.2414074122 | planner/critic 与物理感知多模态合金设计 | AutoZT 不竞争开放式发现定位，而强调运行期确定性、批处理吞吐和运维可观测性 |

## 2. AutoZT 最稳妥的优势

### 2.1 执行底座独立于规划层

图 1 中 LLM 的角色是提出请求、选择技能、诊断异常；图 2 中实际推进由 taskflow 状态机、步骤依赖和物理验证器完成。因此可以把 LLM 换成脚本、不同模型或人工操作，而不改变“输入生成 - 提交 - 监控 - 验证 - 收集 - 推进”的安全边界。这是架构主张，需用替换 planner 的消融实验验证。

### 2.2 将“热输运工作流”作为一等公民

图 2 不是通用工具清单，而是把 KL 计算拆成结构弛豫、力生成、力常数拟合和输运求解，并为每步定义 backend、validator 和 output。相较仅提供 calculator 或 workflow primitives 的系统，这种粒度直接对应研究者判断“能否进入下一步”的物理条件。

### 2.3 面向异构超算，而非单一调度器假设

当前实现覆盖 jzzn 真 SLURM CPU、A800 真 SLURM GPU 和无 SLURM 的 3090 GPU 主机，材料可选择 HPC，VASP 与 MACE 步骤可混合。这个差异应表述为“支持异构后端和调度语义”，不要表述为“其他系统不能跨集群”。

### 2.4 失败是可分类、可恢复的状态

AutoZT 用退出码、OUTCAR/日志判据和步骤状态区分收敛困难、节点故障、排队、挂死和结构/参数错误；`retry` 保留已有产物，`rerun` 才从头重建。图 2 的 known failure 分支可对应到可审计动作，而不是把所有异常交给 LLM 猜测。

### 2.5 产物级 provenance 与会话级重放

每步可以记录输入文件 SHA-256、合并后的参数、生成器/工具信息、作业信息和时间线，并用 `prove --verify` 校验；会话导出把审计、历史、provenance 和 manifest 汇总为可复核档案。这里的卖点不是“首次拥有 provenance”，而是粒度、校验方式和与执行动作的绑定。

### 2.6 稳态成本可预测

自动推进、状态摘要、失败分类和挂死检测不需要每轮调用 LLM。可以把“稳态 0 LLM 调用”作为实现事实，把“比 agentic loop 更省 token”留给基准实验：相同材料数、相同失败率、相同队列条件下，比较 token、墙钟和人工干预次数。

## 3. 建议在论文中避免的表述

- 不要说 DREAMS 没有 HPC、安全、provenance 或确定性；它明确讨论了这些内容。
- 不要说 AtomAgents 没有任何审批；它有 critic agent 审批计划。
- 不要把“论文未讨论”写成“系统不具备”。
- 不要把“可重放”写成“所有 VASP 浮点结果 bit-exact”；应区分输入/产物哈希一致、环境固定和数值结果一致。
- 不要把“支持 3 个后端”扩大成“适用于所有超算”；应给出已验证的后端清单。

## 4. 建议的四组实验

| 实验 | 指标 | 证明什么 |
|---|---|---|
| planner 替换 | 脚本、不同 LLM、人工三种 planner 的完成率/错误率 | 执行底座与规划层解耦 |
| 动作安全 | destructive action 拦截率、误批准率、TTL 过期率、审计完整率 | 工具边界的安全闸门可度量 |
| 重跑一致性 | N 次相同输入的 manifest/hash agreement；环境变化下的差异 | provenance 和重放边界 |
| 规模与成本 | 材料数、失败率、LLM token、墙钟、人工介入次数 | 生产吞吐和稳态成本模型 |

## 5. 可直接使用的定位句

> DREAMS 和 AtomAgents 证明了 agent 栈可以规划、校验并组织材料研究；AutoZT 把已知的 VASP/MACE 计算流程变成跨异构超算的生产执行底座：下一步由物理判据推进，动作经过风险分级和人工闸门，产物按文件哈希归档，失败按状态分类恢复。其核心主张不是拥有某个单独功能，而是当规划层被替换时，执行确定性、可追溯性和运行治理仍然保持不变。

