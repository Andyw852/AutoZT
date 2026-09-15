# -*- coding: utf-8 -*-
"""i18n.py —— 运行期消息的语言开关。

设计（刻意保持零依赖、可增量迁移）：
  · t(zh, en) 按 PHONOAGENT_LANG / LANG 选择；en 时用英文，其余（含未翻译）用中文。
  · 已有 1286 条中文字面量，全部人工翻译是独立工程；因此这里只提供机制 +
    迁移规范，未迁移的字符串保持中文（不会坏）。用 scripts/i18n_report.py 盘点。

迁移规范：把面向用户打印的字符串包一层 t("中文", "English")。
"""
import os

_TRUE = ("1", "true", "yes", "on")


def lang():
    v = (os.environ.get("PHONOAGENT_LANG") or "").strip().lower()
    if not v:
        v = (os.environ.get("LANG") or "").strip().lower()
    return v


def is_en():
    return lang().startswith("en")


def t(zh, en):
    """英文环境返回 en，否则 zh。"""
    return en if is_en() else zh


# ---- 诊断短语的展示层翻译 ----
# 这些短语由远端采集器（_collector_remote.py）生成，随状态一起回传。为避免改推送脚本
# 带来的风险，英文环境下在**展示层**做映射（正则），未知短语原样返回。
_DIAG_RULES = (
    (r"^converged（(\d+)/(\d+) 段）$", r"converged (\1/\2 stages)"),
    (r"^converged（(.+)）$", r"converged (\1)"),
    (r"^dir missing$", "dir missing"),
    (r"^not started$", "not started"),
    (r"^OUTCAR missing$", "OUTCAR missing"),
    (r"^force not converged$", "force not converged"),
    (r"^pressure ([\d.]+) kB > ([\d.]+)$", r"pressure \1 kB > \2"),
    (r"^WAVECAR too small$", "WAVECAR too small"),
    (r"^job PD  (.+)$", r"job pending: \1"),
    (r"^scancel 标记", "marked as cancelled"),
)


def diag(text):
    """英文环境下把诊断短语映射成英文；未命中的原样返回（绝不丢信息）。"""
    if not text or not is_en():
        return text
    import re as _re
    for pat, rep in _DIAG_RULES:
        if _re.match(pat, text):
            return _re.sub(pat, rep, text)
    return text
