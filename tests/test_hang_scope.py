"""Hang recovery must not act on jobs outside the current AutoZT project set."""

import json
from unittest.mock import Mock

import autozt
import autozt.ops as ops


def _node_fail(jobid, workdir):
    return {
        "jobid": jobid,
        "wd": workdir,
        "running": False,
        "age": 0,
        "err_node": True,
        "err_disk": False,
        "bytes": 0,
        "lines": 0,
        "qbytes": 0,
        "last3": [],
    }


def _cfg(tmp_path):
    return {
        "_config_dir": str(tmp_path),
        "hang_check": True,
        "hang_dry_run": False,
        "hang_exclude_hosts": [],
        "hang_max_retries": 2,
        "task_types": {"demo": {}},
    }


def test_hang_recovery_ignores_unmanaged_cluster_job(monkeypatch, tmp_path):
    rec = _node_fail("9001", "/shared/other_project/step1")
    monkeypatch.setattr(ops, "_hung_scan", lambda _cfg, host: (0, json.dumps([rec])))
    cancel = Mock(return_value=(True, ""))
    resume = Mock(return_value=(0, "submitted"))
    monkeypatch.setattr(ops, "_hung_scancel_wait", cancel)
    monkeypatch.setattr(ops, "_hung_resume", resume)
    monkeypatch.setattr(ops, "_hung_state_load", lambda _cfg: {})
    autozt.auto_recover_hung(_cfg(tmp_path), {"types": []})
    cancel.assert_not_called()
    resume.assert_not_called()


def test_hang_recovery_still_handles_managed_job(monkeypatch, tmp_path):
    rec = _node_fail("9002", "/managed/Mn2In2Se5/step1")
    cancel_calls = []
    resume_calls = []
    monkeypatch.setattr(ops, "_hung_scan", lambda _cfg, host: (0, json.dumps([rec])))
    monkeypatch.setattr(
        ops, "_hung_scancel_wait",
        lambda *a, **k: cancel_calls.append(a) or (True, ""),
    )
    monkeypatch.setattr(
        ops, "_hung_resume",
        lambda *a, **k: resume_calls.append(a) or (0, "submitted"),
    )
    monkeypatch.setattr(ops, "_hung_state_load", lambda _cfg: {})
    monkeypatch.setattr(ops, "_hung_state_save", lambda *_a: None)
    monkeypatch.setattr(autozt, "log_action", lambda *_a, **_k: None)
    data = {"types": [{
        "key": "demo",
        "materials": [{
            "name": "Mn2In2Se5",
            "tt": "demo",
            "path": "/managed/Mn2In2Se5",
            "host_eff": "jzzn",
            "steps": [{"_host": "jzzn"}],
            "ps": {"setting": {}},
        }],
    }]}

    cfg = _cfg(tmp_path)
    autozt.auto_recover_hung(cfg, data)

    assert cancel_calls == [(cfg, "9002")]
    assert resume_calls == [(cfg, "/managed/Mn2In2Se5/step1")]
