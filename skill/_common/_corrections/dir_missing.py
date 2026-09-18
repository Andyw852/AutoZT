# -*- coding: utf-8 -*-
"""dir_missing —— 步骤目录缺失纠错 handler（示范 destructive 风险等级）。

现场：diag = "dir missing"。远端该步骤目录整个不见了（被 clean 过 / 手工删过 /
从没生成过）。恢复要重新生成输入，属于 **rerun 级** 操作。

★ 这类 handler **只给建议、不给 apply**：tf correct 永远只打印命令，
   真正的 rerun 必须由人执行（AGENTS.md 铁律：破坏性操作先请示）。
"""
from .base import CorrectionHandler, Suggestion


class DirMissingCorrection(CorrectionHandler):
    name = "dir_missing"
    title = "步骤目录缺失（需重新生成）"
    matches = ("dir missing", "目录缺失")
    steps = ("*",)
    risk = "destructive"

    def suggest(self, ctx):
        base = self.retry_cmd(ctx)
        rerun = base.rsplit(" ", 1)[0] + " rerun"
        return Suggestion(
            handler=self.name, title=self.title, risk=self.risk, auto=False,
            reason="远端目录不存在，retry 只能重生成输入文件、无法恢复已删产物。",
            actions=["先 retry 试一次（多数情况重生成输入足够）",
                     "确实要推倒重来 → rerun（★ 破坏性：删该步目录，先请示）",
                     "若该步是扇出步骤（如 step5_label）：只补未完成的用 retry，"
                     "绝不 rerun/clean"],
            commands=[base, rerun],
            changes=[{"target": "step_directory", "operation": "regenerate",
                      "path": ctx.workdir or ""}],
            verify=["确认目标步骤目录和输入文件已由 AutoZT 重新生成",
                    "确认没有误删除已有产物"],
            rollback=["rerun 会删除步骤目录，无法由 handler 自动回滚；执行前必须人工确认"],
            requires_approval=True,
            execute_via="autozt retry（优先）或经批准的 autozt rerun",
            notes="受保护材料（目录下有 AGENTS-PROTECTED.md）一律不许 rerun/clean。")


HANDLER = DirMissingCorrection()
