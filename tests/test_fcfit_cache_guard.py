"""pheasy cs.pkl/ns 缓存串档防护：每档截断前清缓存 + 校验 -s 日志里的 c3。"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common" / "fcfit"))
import fc_fit_driver as d  # noqa: E402


def test_clear_cache_removes_only_pheasy_cache(tmp_path):
    for n in ("cs.pkl", "cs.pkl.meta.json", "ns_harm.npz", "ns_anharm3.npz.meta.json",
              "fc3.hdf5", "disp_matrix.pkl"):
        (tmp_path / n).write_text("x")
    d._clear_pheasy_cache(tmp_path)
    left = sorted(p.name for p in tmp_path.iterdir())
    assert left == ["disp_matrix.pkl", "fc3.hdf5"]


def test_cutoff_check_accepts_match_and_rejects_stale():
    ok = "Cutoff distance (A) for 3-order IFCs: 5.12\nANHARM3 | cluster number: 182"
    d._check_cutoff_honoured(ok, 5.12, "t")
    stale = "Cutoff distance (A) for 3-order IFCs: 2.73\nANHARM3 | cluster number: 15"
    with pytest.raises(SystemExit):
        d._check_cutoff_honoured(stale, 5.12, "t")


def test_cutoff_check_skips_when_no_c3():
    d._check_cutoff_honoured("anything", None, "t")
