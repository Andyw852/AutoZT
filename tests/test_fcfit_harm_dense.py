"""HARM_DENSE 门闸的回归测试：稀疏方法没保留全部二阶列时必须判 FAIL。

背景：Mg4C60/MoS2 的 RFE-OLS 拟合在 PHEASY_HARM_DENSE 未开的情况下把二阶块
一起剪掉（Mg4C60 3.0 Å 只剩 36,044/74,017），却一路 PASS 到声子/κ。门闸要求：
  * manifest.environment.PHEASY_HARM_DENSE 为真；
  * manifest.results 里有 nnz_unpenalized / unpenalized_columns；
  * 且 nnz_unpenalized == unpenalized_columns。
缺字段 = FAIL（不跳过），字段不相等 = FAIL。
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common" / "fcfit"))

import fc_fit_driver as d  # noqa: E402


def _manifest(tmp_path, env=None, results=None):
    (tmp_path / "fit_manifest.json").write_text(
        json.dumps({"environment": env or {}, "results": results or {}}),
        encoding="utf-8")


def test_missing_fields_fails(tmp_path):
    """manifest 缺 nnz_unpenalized / unpenalized_columns：必须 FAIL，不能跳过。"""
    _manifest(tmp_path, {"PHEASY_HARM_DENSE": "1"}, {})
    metrics = d._harm_dense_gate("RFE-OLS", 3, tmp_path, {})
    assert metrics["harm_dense_env"] == "1"
    assert metrics["harm_kept"] is None and metrics["harm_free"] is None
    assert (tmp_path / ".fit_gate_fail").is_file()
    assert "lacks nnz_unpenalized" in (tmp_path / ".fit_gate_fail").read_text()


def test_partial_fails(tmp_path):
    """二阶列被剪（保留 < 自由）：必须 FAIL。"""
    _manifest(tmp_path, {"PHEASY_HARM_DENSE": "1"},
              {"nnz_unpenalized": 100, "unpenalized_columns": 200})
    metrics = d._harm_dense_gate("RFE-OLS", 3, tmp_path, {})
    assert (metrics["harm_kept"], metrics["harm_free"]) == (100, 200)
    assert (tmp_path / ".fit_gate_fail").is_file()
    assert "pruned" in (tmp_path / ".fit_gate_fail").read_text()


def test_full_kept_passes(tmp_path):
    """二阶全保留：不写 .fit_gate_fail，计数写进 metrics（post 会并入 summary）。"""
    _manifest(tmp_path, {"PHEASY_HARM_DENSE": "1"},
              {"nnz_unpenalized": 200, "unpenalized_columns": 200})
    metrics = d._harm_dense_gate("RFE-OLS", 3, tmp_path, {"pheasy_rmse_train": 1.0})
    assert (metrics["harm_kept"], metrics["harm_free"]) == (200, 200)
    assert metrics["pheasy_rmse_train"] == 1.0
    assert not (tmp_path / ".fit_gate_fail").exists()


def test_flag_off_fails_even_with_counts(tmp_path):
    """字段齐全但 HARM_DENSE 没开：FAIL（防止旧 pheasy 或环境漏传）。"""
    _manifest(tmp_path, {}, {"nnz_unpenalized": 200, "unpenalized_columns": 200})
    d._harm_dense_gate("LASSO", 3, tmp_path, {})
    assert (tmp_path / ".fit_gate_fail").is_file()


def test_fc2_only_and_dense_methods_skip(tmp_path):
    """纯 fc2（enable_fc<3）与非稀疏方法（OLS/RIDGE）不适用：不写 FAIL。"""
    _manifest(tmp_path, {"PHEASY_HARM_DENSE": "1"}, {})
    d._harm_dense_gate("RFE-OLS", 2, tmp_path, {})
    d._harm_dense_gate("OLS", 3, tmp_path, {})
    assert not (tmp_path / ".fit_gate_fail").exists()
