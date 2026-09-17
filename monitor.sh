#!/bin/bash
# =============================================================================
# monitor.sh —— 十分钟巡检：推进(静默) + 采集 diff(无变化 0 输出)
#
# 用法（crontab，每 10 分钟）：
#   */10 * * * * "$(cd "$(dirname "$0")" && pwd)"/monitor.sh >> "$(cd "$(dirname "$0")" && pwd)"/.tf_monitor.out 2>&1
#   ★ 2026-09-16 修正：原来这里写的是 'cd "$(dirname "$0")" <仓库绝对路径> pwd'，
#     cd 收到两个参数必然失败 → 顶层那行 '|| exit 1' 让整个脚本一行都不跑（HEAD 里
#     就是这个状态，属于某次改名替换替换坏了）。现在改成标准写法。
#
# 省 token 的关键：
#   1. autozt auto on 的推进日志静默到 .tf_monitor.log（不刷屏）；
#   2. autozt summary --diff 无变化时输出 0 字节，只有状态真变（有作业完成/新失败/
#      排队变化）才打印几行汇总——agent/人只看这最后几行。
# =============================================================================
cd "$(cd "$(dirname "$0")" && pwd)" || exit 1

# ★ 修复：cron 环境 PATH 极简（/usr/bin:/bin），不含 ~/.local/bin，
#   导致裸调 tf 报 "No such file or directory"。显式补上，让 tf 软链可被找到。
export PATH="/home/wangchao/.local/bin:$PATH"

# =============================================================================
# ★ 与 taskflow-v2.0 的隔离开关（2026-09-17 加）
#   下面两处会【真的提交作业】：内层块的 `auto on`（ke-dft-cpu 的 Mg4C60 两个材料）
#   和第 1) 步的全局 `autozt auto on`。它们面向的 /mnt/d/tf_data/work_taskflow 与
#   taskflow-v2.0 是【同一批真实项目、同一批超算账号】，而 v2.0 那边还有 5 个 monitor
#   在跑 —— 这里再推进一次就是重复投作业、互相改远端输入。
#   故默认【只跑只读巡检】（下面的 summary --diff）。确实要本仓推进时显式打开：
#       AUTOZT_MONITOR_SHARED=1 ./monitor.sh
# =============================================================================
if [ "${AUTOZT_MONITOR_SHARED:-0}" = "1" ]; then
(
    # ★ 2026-09-17 修正：原来这里硬编码 cd 到 /home/wangchao/software/taskflow-v2.0，
    #   而 v2.0 里只有 bin/tf、没有 bin/autozt → 每轮 python3 报 "can't open file"
    #   静默失败，ke-dft-cpu 的这两个材料实际从没被这一块推进过。
    #   改成指向【本脚本所在仓库】。
    cd "$(cd "$(dirname "$0")" && pwd)" || exit 1
    /home/wangchao/bin/hanhai25-connect >> tmp/ke_auto_monitor.log 2>&1 || exit 1
    AUTOZT_OP_WORKERS=8 timeout 600 python3 bin/autozt -tt ke-dft-cpu -p Mg4C60,Mg4C60_monolayer auto on >> tmp/ke_auto_monitor.log 2>&1
    timeout 600 python3 bin/autozt -tt ke-dft-cpu -p Mg4C60,Mg4C60_monolayer summary --diff
)

# 1) 推进流水线（DAG 自动推进：按依赖找就绪步骤，S0 FAIL 不再阻塞 S3；FAIL 只报告不动），日志静默
AUTOZT_OP_WORKERS=8 timeout 600 python3 bin/autozt auto on >> .tf_monitor.log 2>&1
else
    echo "[隔离模式] 已跳过会提交作业的推进（taskflow-v2.0 监控仍在跑同一批项目）。" \
         "需要本仓推进时：AUTOZT_MONITOR_SHARED=1 $0"
fi

# 2) 采集 + 变更检测：无变化 0 输出，有变化才打印汇总 + FAIL 清单
AUTOZT_OP_WORKERS=8 timeout 600 python3 bin/autozt summary --diff 2>&1
