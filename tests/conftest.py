# -*- coding: utf-8 -*-
"""测试默认让 agent/MCP 读状态走「现采」（被测试替换的假 CLI）。

source=auto 会优先读配置目录里的 .tf_progress.json；在装着真实配置的机器上跑测试时，
那会读到真实生产状态而不是测试构造的假数据。需要测进度文件来源的用例显式传 source。
"""
import os

os.environ["AUTOZT_AGENT_SOURCE"] = "live"
