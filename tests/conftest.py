# -*- coding: utf-8 -*-
"""测试默认让 agent/MCP 读状态走「现采」（被测试替换的假 CLI）。

source=auto 会优先读配置目录里的 .tf_progress.json；在装着真实配置的机器上跑测试时，
那会读到真实生产状态而不是测试构造的假数据。需要测进度文件来源的用例显式传 source。
"""
import os

os.environ["AUTOZT_AGENT_SOURCE"] = "live"

# 审批卡片默认会在子进程里采集目标材料的任务图：测试里一律关掉，免得在装着真实配置的
# 机器上为一条假命令去采集真实项目（ssh、最多 90 秒）。需要测任务图的用例显式打开。
os.environ.setdefault("AUTOZT_APPROVAL_GRAPH", "0")
