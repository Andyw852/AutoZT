autozt --help 输出（自动生成）

    autozt [-tt <skill>] [-p <material>] [-j <step>] <command>

完整命令列表：

    用法：autozt [-tt 技能] [-p 材料] [-j 步骤] 命令
    
    查看：
      summary [--diff]    只读汇总；--diff 无变化静默
      list               只读总表；--refresh 强制刷新
      status             状态详情，同时可能拉结果、自动提交
      probe              单材料作业诊断（需 -p）
    
    计算：
      start              生成缺失输入并提交
      stop               取消作业并阻止自动重跑
      retry              保留产物重新生成输入，不提交；检查后 start
      rerun              删除旧步骤并重新生成，不提交；检查后 start
      clean              删除产物，不重新生成
    
    管理：
      auto on|off        开关自动推进；-tt/-p 限定范围
      auto resume        恢复指定项目取消步骤，必须配 -tt 和 -p
      monitor [-d]       持续监控；--stop 停止，--restart 重启
      conf [--set ...]   查看或修改步骤配置
      hpc 集群           切换指定项目集群（需 -p）
      init / skills      初始化项目 / 查看技能
    
    技能自描述与可追溯（v1.0，全部只读）：
      skill [show] [技能] 技能卡片：每步 生成器/工具/判据/产物（论文图 2 那种）
                          （--json 机器可读；不给技能名=列出全部技能一行摘要）
      schema [技能]      看技能吃什么/吐什么/有哪些旋钮/能接谁/怎么纠错
                         （--json 机器可读，--strict 有错误返回非零）
      correct -p 材料    把 FAIL 诊断喂给 _corrections/ 纠错 handler 库，
                         列出建议；-y 才执行 handler.apply（只改输入文件、不提交）
      history [-p 材料]  步骤状态的时间序列（history.jsonl，采集时自动记录）
      prove -p 材料      这一步"结果怎么来的"：输入 sha256 / step.conf 参数 /
                         工具版本 / 作业号（gen 时自动落档，--verify 校验输入没被改）
    
    AI 审计（v1.0 P0-1，agent 走网关：风险分档 + 每次调用留痕）：
      act <命令>         agent 的唯一入口：autozt act -p 材料 summary / autozt act start …
                         （只读与推进类放行并记账；stop/rerun/clean/-f/-y 需人工批准）
      act log           看审计流水（.tf_agent_log.jsonl）；act policy 看风险分档表
      approve <命令>     人工在**交互终端**批准一条破坏性动作（一次性令牌，默认 15 分钟）
    
    旧命令和别名继续兼容。高级命令、全部参数及示例：tf --help-all
    注意：status/auto/monitor 可提交作业；只看状态用 summary 或 list。
    
