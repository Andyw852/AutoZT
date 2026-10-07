"""`_pheasy_metrics` 读 fit_manifest.json 而非正则抓日志的回归防护。

锁住 job 5830 之后的另一个坑：pheasy 的 RFE 系方法在日志里先打印
``best CV RMSE: <x>`` 再打印 ``RMSE: <y>``，旧的 ``\\bRMSE:`` 正则第一次命中 CV 行，
把 CV 误差存进 ``pheasy_rmse_eV_per_A``（train 键）——方法对比表里 RFE 会显得比实际差。
修复后 train/CV 直接从 fit_manifest.json 的 metrics.rmse / metrics.re /
metrics.rmse_path_mean 读，日志只用来抓 manifest 没有的标记（corr/free_ifcs/alpha/GPU）。
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common" / "fcfit"))
import fc_fit_driver as d  # noqa: E402


# 真实 RFE 日志片段（5849 Mn2In2Se5 step1_fit/methods/m-pheasy-RFE-OLS/run.log）：
# "best CV RMSE:" 在 "RMSE:" 之前，旧正则第一次就命中 CV 值。
RFE_LOG = """
[RFE] argmin: n_active=745, CV=2.657390e-03
[RFE] selected: n_active=671, CV=2.664949e-03 (+-2.57e-05)
- best CV RMSE: 0.002664949041122533 eV/A
- RMSE: 0.0024863780637438097 eV/A
- Relative error: 0.017363005208632862
worst corr=0.9999
Free IFC terms: 3638
best alpha= 1.0e-11
"""


def _write_manifest(tmp_path):
    """pheasy 的 fit_manifest.json：metrics.rmse=训练、rmse_path_mean=CV、re=训练相对误差。"""
    (tmp_path / "fit_manifest.json").write_text(json.dumps({
        "metrics": {
            "rmse": 0.0009153184898178548,          # train RMSE（与日志 "RMSE:" 0.002486 不同！）
            "re": 0.006391899903723942,             # train relative error 0.639%
            "rmse_path_mean": 0.0011378242547216261,  # CV RMSE（5 折按构型分组）
            "mse_path_mean": 1.294644034632824e-06,
        },
    }), encoding="utf-8")


def test_pheasy_metrics_reads_manifest_not_log(tmp_path):
    """train/CV 必须来自 fit_manifest.json，而不是日志里那个歧义的 "RMSE:" 行。"""
    _write_manifest(tmp_path)
    m = d._pheasy_metrics(RFE_LOG, out=tmp_path)
    # train：来自 manifest（日志 "RMSE:" 是 0.002486，manifest metrics.rmse 是 0.0009153）
    assert m["pheasy_rmse_train"] == pytest.approx(0.0009153184898178548)
    assert m["pheasy_relative_error"] == pytest.approx(0.006391899903723942)
    # CV：来自 manifest metrics.rmse_path_mean
    assert m["pheasy_rmse_cv"] == pytest.approx(0.0011378242547216261)
    # 旧的歧义字段不应再出现
    assert "pheasy_rmse_eV_per_A" not in m


def test_pheasy_metrics_log_only_markers_still_regexed(tmp_path):
    """manifest 没有的标记（corr/free_ifcs/alpha）仍从日志抓。"""
    _write_manifest(tmp_path)
    m = d._pheasy_metrics(RFE_LOG, out=tmp_path)
    assert m["pheasy_worst_force_correlation"] == pytest.approx(0.9999)
    assert m["pheasy_free_ifcs"] == 3638
    assert m["pheasy_best_alpha"] == pytest.approx(1.0e-11)


def test_pheasy_metrics_no_manifest(tmp_path):
    """没有 manifest 时不报错，train/CV 字段留空（不伪造）。"""
    m = d._pheasy_metrics(RFE_LOG, out=tmp_path)
    assert "pheasy_rmse_train" not in m
    assert "pheasy_rmse_cv" not in m
    assert "pheasy_relative_error" not in m
    # 日志标记照常
    assert m["pheasy_free_ifcs"] == 3638


def test_rel_err_gate_pheasy_missing_manifest_fails(tmp_path):
    """pheasy 缺 manifest（rel=None）且门开着时，必须写 .fit_gate_fail。"""
    assert d._rel_err_gate("pheasy", None, 0.02, tmp_path) is True
    assert (tmp_path / ".fit_gate_fail").is_file()
    assert "unavailable" in (tmp_path / ".fit_gate_fail").read_text()


def test_rel_err_gate_non_pheasy_no_marker(tmp_path):
    """phono3py / hiphive 没有 pheasy_relative_error，绝不能因缺字段判失败。"""
    for eng in ("phono3py", "hiphive", ""):
        assert d._rel_err_gate(eng, None, 0.02, tmp_path) is False
        assert not (tmp_path / ".fit_gate_fail").exists()


def test_rel_err_gate_passes_below_threshold(tmp_path):
    """pheasy 且 rel 低于阈值时不写。"""
    assert d._rel_err_gate("pheasy", 0.005, 0.02, tmp_path) is False
    assert not (tmp_path / ".fit_gate_fail").exists()


