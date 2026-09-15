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
