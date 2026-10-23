# -*- coding: utf-8 -*-
"""[V192] QOS 作业数/提交数上限排队的识别与提示。

背景（2026-10-23 实测）：jzzn 的 premium QoS MaxJobsPU=5。把 2D VASP 提交模板全局
切到 premium 后，LS 的 6 个 S7 帧全部 PD(QOSMaxJobsPerUserLimit)、squeue --start 直接
显示 N/A，看起来像"集群没资源"，实际是同一用户的其它项目（Mg2C60）占满了 5 个槽位；
而 regular 的 MaxJobsPU=50。以前程序只把原因原样显示成 PD(...)，不说为什么、怎么办。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from autozt import preflight as pf  # noqa: E402
from autozt import report  # noqa: E402


def _step(label, kind="PD", info="QOSMaxJobsPerUserLimit", qos="premium"):
    return {"label": label, "kind": kind,
            "job": {"state": kind, "info": info, "qos": qos} if kind == "PD" else None}


def _mats(name="LS", steps=None, tt="ke-dft-cpu"):
    return [({"key": tt}, {"tt": tt, "name": name, "steps": steps or []})]


# ---------------------------------------------------------------- 原因识别
def test_recognises_all_user_and_account_limits():
    for r in ("QOSMaxJobsPerUserLimit", "QOSMaxSubmitJobPerUserLimit",
              "QOSMaxJobsPerAccountLimit", "QOSMaxSubmitJobPerAccountLimit"):
        assert pf.is_qos_limit_reason(r), r


def test_tolerates_squeue_parentheses_and_spaces():
    assert pf.is_qos_limit_reason("(QOSMaxJobsPerUserLimit)")
    assert pf.is_qos_limit_reason("  QOSMaxJobsPerUserLimit  ")


def test_ignores_other_pending_reasons():
    for r in ("Priority", "PartitionConfig", "Dependency", "Resources",
              "ReqNodeNotAvail", "", None):
        assert not pf.is_qos_limit_reason(r), r


# ---------------------------------------------------------------- 文案
def test_message_explains_cause_and_remedies():
    msg = pf.qos_limit_message("QOSMaxJobsPerUserLimit", "premium")
    assert "QOSMaxJobsPerUserLimit" in msg
    assert "premium" in msg
    for kw in ("不是集群没资源", "等槽位", "换到限额更大的 QoS", "分批"):
        assert kw in msg, kw


def test_message_says_do_not_lower_cores():
    # 铁律 9：降核不解决槽位上限，还会改并行设置
    assert "别为这个降核数" in pf.qos_limit_message("QOSMaxJobsPerUserLimit")


def test_message_distinguishes_job_limit_from_submit_limit():
    assert "作业数" in pf.qos_limit_message("QOSMaxJobsPerUserLimit")
    assert "提交数" in pf.qos_limit_message("QOSMaxSubmitJobPerUserLimit")


def test_message_without_qos_name():
    msg = pf.qos_limit_message("QOSMaxJobsPerUserLimit")
    assert "QOSMaxJobsPerUserLimit" in msg
    assert "QoS=" not in msg


# ---------------------------------------------------------------- 分组
def test_groups_by_qos_and_skips_non_qos():
    jobs = [
        {"state": "PD", "info": "QOSMaxJobsPerUserLimit", "qos": "premium"},
        {"state": "PD", "info": "QOSMaxJobsPerUserLimit", "qos": "premium"},
        {"state": "PD", "info": "Priority", "qos": "regular"},      # 普通排队
        {"state": "R", "info": "QOSMaxJobsPerUserLimit", "qos": "premium"},  # 在跑的不算
    ]
    g = pf.qos_limit_groups(jobs)
    assert list(g) == ["premium"]
    assert len(g["premium"]) == 2


# ---------------------------------------------------------------- 提示行（纯函数）
def test_no_hint_when_only_normal_queueing():
    assert report.qos_limit_hint_lines(
        _mats(steps=[_step("S7", info="Priority", qos="regular")])) == []


def test_no_hint_when_steps_have_no_job():
    assert report.qos_limit_hint_lines(
        _mats(steps=[{"label": "S1", "kind": "OK", "job": None}])) == []


def test_hint_lists_qos_count_and_titles():
    steps = [_step("S7_deform-%d" % i) for i in range(3)]
    lines = report.qos_limit_hint_lines(_mats(steps=steps))
    text = "\n".join(lines)
    assert "QOS 槽位上限排队" in text
    assert "QoS=premium" in text
    assert "3 个步骤被挡" in text
    assert "S7_deform-0" in text


def test_hint_truncates_after_five_steps():
    steps = [_step("S7_deform-%d" % i) for i in range(8)]
    text = "\n".join(report.qos_limit_hint_lines(_mats(steps=steps)))
    assert "另有 3 个" in text


def test_hint_groups_multiple_qos_separately():
    steps = [_step("A", qos="premium"), _step("B", qos="gpu")]
    text = "\n".join(report.qos_limit_hint_lines(_mats(steps=steps)))
    assert "QoS=premium" in text and "QoS=gpu" in text


def test_hint_survives_missing_qos_field():
    steps = [{"label": "S7", "kind": "PD",
              "job": {"state": "PD", "info": "QOSMaxJobsPerUserLimit"}}]  # 老格式没有 qos
    text = "\n".join(report.qos_limit_hint_lines(_mats(steps=steps)))
    assert "未采集" in text


# ---------------------------------------------------------------- 扇出步骤
def test_fanout_uses_sub_job_reasons():
    """扇出步骤的 job.info 被压成 "4/10 0R 6PD"，原因要从 fan_pd_reasons 明细里取。"""
    step = {"label": "S7_deform", "kind": "PD",
            "job": {"state": "PD", "info": "4/10 0R 6PD", "qos": "premium"},
            "fan_pd_reasons": [{"reason": "QOSMaxJobsPerUserLimit",
                                "qos": "premium", "n": 6}]}
    text = "\n".join(report.qos_limit_hint_lines(_mats(steps=[step])))
    assert "QoS=premium" in text
    assert "QOSMaxJobsPerUserLimit" in text
    assert "S7_deform" in text


def test_fanout_with_normal_queueing_gets_no_hint():
    step = {"label": "S7_deform", "kind": "PD",
            "job": {"state": "PD", "info": "4/10 0R 6PD", "qos": "regular"},
            "fan_pd_reasons": [{"reason": "Priority", "qos": "regular", "n": 6}]}
    assert report.qos_limit_hint_lines(_mats(steps=[step])) == []


def test_fanout_mixed_reasons_only_reports_qos_ones():
    step = {"label": "S7_deform", "kind": "PD",
            "job": {"state": "PD", "info": "4/10 0R 6PD", "qos": "regular"},
            "fan_pd_reasons": [
                {"reason": "Priority", "qos": "regular", "n": 4},
                {"reason": "QOSMaxJobsPerUserLimit", "qos": "regular", "n": 2}]}
    text = "\n".join(report.qos_limit_hint_lines(_mats(steps=[step])))
    assert "QOSMaxJobsPerUserLimit" in text
    assert "Priority" not in text

