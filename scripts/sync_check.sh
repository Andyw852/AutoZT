#!/usr/bin/env bash
# sync_check.sh —— 内容级对齐审计：本仓 skill/ + setting/ 与【参考仓】逐文件 diff。
#
# ★ 2026-09-17 更新：AutoZT 已是唯一权威仓库，旧仓 ~/software/taskflow（v1.0）不再维护，
#   默认参考仓改成仍在用的姊妹仓 ~/software/taskflow-v2.0 —— 两边共用同一批真实项目/
#   超算账号，skill/ 一旦漂移，同一个项目会被生成出不同输入（会互相踩）。
#   本脚本只读、只审计，不自动同步；差异由人工判断方向后再动手。
# 用法:
#   scripts/sync_check.sh                      # 对比默认参考仓 ~/software/taskflow-v2.0
#   scripts/sync_check.sh /path/to/reference   # 对比指定参考仓
# 退出码: 0=完全对齐；1=有差异（内容级，排除 __pycache__/运行时缓存）
set -uo pipefail

ORIG="${1:-/home/wangchao/software/taskflow-v2.0}"
V2="$(cd "$(dirname "$0")/.." && pwd)"

EXC='__pycache__|\.pyc$|\.tf_(state_cache|hung|summary|watch)'

dirty=0
for sub in skill setting; do
  diff_out=$(diff -rq "$V2/$sub" "$ORIG/$sub" 2>/dev/null | grep -vE "$EXC")
  if [ -n "$diff_out" ]; then
    echo "[$sub/] 与参考仓不一致："
    echo "$diff_out" | sed 's/^/  /'
    dirty=1
  fi
done

if [ "$dirty" -eq 0 ]; then
  echo "OK：skill/ 与 setting/ 已与参考仓内容级对齐。"
  exit 0
fi

echo ""
echo "发现差异。★ 先确认参考仓没有在跑 monitor（否则会同时提交作业），再选方向："
echo "  AutoZT 已是权威仓，正常方向是 本仓 → 参考仓（只同步 skill/，参考仓的 setting/ 是它自己的）："
echo "    rsync -ac --exclude='__pycache__/' --exclude='*.pyc' '$V2/skill/' '$ORIG/skill/'"
echo "  只有确认本仓落后时，才反向拉（会覆盖本仓，先自行备份）："
echo "    rsync -ac --exclude='__pycache__/' --exclude='*.pyc' '$ORIG/skill/' '$V2/skill/'"
exit 1
