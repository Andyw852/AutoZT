#!/usr/bin/env bash
# tf-push.sh —— 一步：本地真实版提交 + 脱敏快照推送（单历史，主工作树零污染）
#
# 用法: scripts/tf-git-push.sh [-p <路径> ...] [--dry-run] [--yes] "commit message"
#   -p <路径>   只提交/推送指定路径（可多个）；省略 = 全量 git add -A
#   --dry-run   只生成脱敏提交，不 push（安全预演）
#   --yes       跳过交互确认（等价 AUTOZT_PUSH_YES=1）
#
# 与旧版的区别：旧版切 push-sanitized 分支 + read-tree 强改主工作树 + stash + trap pop，
# 任一步失败就留 UU 冲突/stash 残局；新版在临时 git worktree 里干活，主工作树全程不动，
# 远端始终只收到基于 origin/<分支> 的一条脱敏线性历史（每次 fast-forward）。
set -euo pipefail

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "$SELF")" && pwd)"
REPO="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
cd "$REPO"

PATHS=(); MSG=""; DRY=0; ANS=""
while [ $# -gt 0 ]; do
  case "$1" in
    -p|--path) [ $# -gt 1 ] || { echo "缺 -p 参数"; exit 1; }; PATHS+=("$2"); shift 2;;
    --dry-run) DRY=1; shift;;
    --yes) AUTOZT_PUSH_YES=1; shift;;
    *) MSG="$1"; shift;;
  esac
done
[ -n "$MSG" ] || MSG="update"

# 词表优先级：仓库内 setting/git-sanitize.conf 优先
#   （与远端脱敏线一致；实测 ~/.config 旧词表会把已脱敏内容二次改写，造出大量伪差异）
#   > ~/.config/git-sanitize.conf 仅作后备
CONF=""
for c in "$REPO/setting/git-sanitize.conf" "$HOME/.config/git-sanitize.conf"; do
  [ -f "$c" ] && { CONF="$c"; break; }
done
[ -n "$CONF" ] || { echo "缺敏感词映射：setting/git-sanitize.conf 或 $HOME/.config/git-sanitize.conf"; exit 1; }

git remote get-url origin >/dev/null 2>&1 || { echo "无 origin 远程"; exit 1; }
CUR="$(git symbolic-ref --short HEAD 2>/dev/null || echo master)"

if ! grep -q 'ssh.github.com' "$HOME/.ssh/config" 2>/dev/null; then
  printf '\nHost github.com\n    HostName ssh.github.com\n    Port 443\n    User git\n    ConnectTimeout 10\n    ServerAliveInterval 30\n    ServerAliveCountMax 2\n' >> "$HOME/.ssh/config" 2>/dev/null || true
  echo "已配置 github SSH 走 443"
fi

echo "===== 待提交改动 ====="
if [ ${#PATHS[@]} -gt 0 ]; then git status --short -- "${PATHS[@]}"; else git status --short; fi
echo "======================"

if [ ${#PATHS[@]} -gt 0 ]; then
  mapfile -t _py < <(git status --porcelain -- "${PATHS[@]}" | sed 's/^...//; s/.* -> //' | grep '\.py$' || true)
else
  mapfile -t _py < <(git status --porcelain | sed 's/^...//; s/.* -> //' | grep '\.py$' || true)
fi
for f in "${_py[@]:-}"; do
  [ -n "$f" ] && [ -f "$f" ] || continue
  python3 -m py_compile "$f" || { echo "❌ py_compile 失败: $f"; exit 1; }
done
[ ${#_py[@]} -gt 0 ] && echo "py_compile 通过（${#_py[@]} 个 .py）"

if [ "${AUTOZT_PUSH_YES:-0}" != "1" ]; then
  read -r -p "确认 sanitize 并推到 GitHub（远端 $CUR）？[y/N] " ANS
  case "$ANS" in y|Y|yes) ;; *) echo "已取消"; exit 1;; esac
fi

# 1) 本地真实版提交（留在本地分支，永不推到远端）
if [ ${#PATHS[@]} -gt 0 ]; then git add -- "${PATHS[@]}"; else git add -A; fi
git commit -m "$MSG（本地真实版）" >/dev/null 2>&1 || echo "(本地无变更可提交，直接用当前树快照)"

# 2) 临时 worktree：基于远端最新，取真实版内容 + 脱敏 + 提交
git fetch origin --quiet
BASE=""
for b in "origin/$CUR" origin/main origin/master; do
  git rev-parse --verify "$b" >/dev/null 2>&1 && { BASE="$b"; break; }
done
[ -n "$BASE" ] || { echo "找不到远程基线分支"; exit 1; }

WT="$(mktemp -d /tmp/tf-push.XXXXXX)"
cleanup() { cd "$REPO"; git worktree remove --force "$WT" >/dev/null 2>&1 || true; }
trap cleanup EXIT

git worktree add --detach "$WT" "$BASE" >/dev/null
if [ ${#PATHS[@]} -gt 0 ]; then
  git -C "$WT" checkout "$CUR" -- "${PATHS[@]}"
else
  git -C "$WT" read-tree -u --reset "$CUR"
fi
cd "$WT"
git add -A
python3 "$SCRIPT_DIR/sanitize.py" "$CONF"
git add -A
if git diff --cached --quiet HEAD; then
  echo "远端已包含当前内容（脱敏后无差异），跳过 push。"
  exit 0
fi
git commit -m "$MSG（脱敏版）" >/dev/null
echo "已生成脱敏提交: $(git rev-parse --short HEAD)"

if [ "$DRY" = "1" ]; then
  echo "--dry-run：不 push，远端 $CUR 未改动。"
  exit 0
fi

ok=0
for i in 1 2 3 4 5; do
  echo "push 第 $i 次..."
  if timeout 60 git push origin HEAD:"$CUR"; then ok=1; break; fi
  sleep 4
done
[ "$ok" = "1" ] || { echo "❌ push 失败（5 次）"; exit 1; }
echo "✅ 完成：远端 $CUR 已是脱敏版，主工作树未改动。"
