"""[0052] 本地结果回拉：
1. result_flat 步骤（zt-dft-cpu S20_zt）直接拉进 result/，旧版 result/<步骤>/ 子目录删除；
2. 整目录回拉（fetch_all）是镜像：上一轮拉回、远端这轮已没有的顶层项本地删除；
3. fetch_extra（transport_plots/）同样镜像，第一次（无清单）也按本地现有项清；远端什么都没有时不删。"""
import sys
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import autozt                      # noqa: E402
from autozt import workflow        # noqa: E402


def _mat(tmp_path, flat, extra=None):
    mat = tmp_path / "Si"
    remote = mat / "zt-dft-cpu" / "step20_zt"
    remote.mkdir(parents=True)
    sc = {"name": "step20_zt", "label": "S20_zt", "fetch_all": True, "done_marker": "zt_summary.json"}
    if flat:
        sc["result_flat"] = True
    if extra:
        sc["fetch_extra"] = extra
    s = {"name": "step20_zt", "label": "S20_zt", "dir": str(remote), "exists": True, "done": True}
    m = {"name": "Si", "tt": "zt-dft-cpu", "result_dir": str(mat / "zt-dft-cpu" / "result"),
         "fetch_files": ["zt_summary.json"], "steps": [s], "_seg": {"steps_cfg": [sc]}}
    return m, s, remote, Path(m["result_dir"])


def _fetch(m):
    with patch.object(autozt, "log_action"):
        assert workflow.fetch_material({}, m, quiet=True)


def test_flat_result_and_legacy_removed(tmp_path):
    m, s, remote, rd = _mat(tmp_path, flat=True)
    (remote / "zt_summary.json").write_text("{}")
    (remote / "per_temperature" / "300K").mkdir(parents=True)
    (remote / "per_temperature" / "300K" / "report.txt").write_text("x")
    (rd / "step20_zt").mkdir(parents=True)                     # 旧版布局
    (rd / "step20_zt" / "zt_summary.json").write_text("old")
    (rd / ".autozt_inputs").mkdir()
    (rd / ".autozt_inputs" / "step20_zt.json").write_text("{}")
    _fetch(m)
    assert (rd / "zt_summary.json").read_text() == "{}"
    assert (rd / "per_temperature" / "300K" / "report.txt").is_file()
    assert not (rd / "step20_zt").exists()
    assert (rd / ".autozt_inputs" / "step20_zt.json").is_file()  # 本地自己的记录不动
    assert workflow.step_result_dir(m, s) == str(rd)
    assert workflow._fetch_receipt_valid({}, m, s) in (True, False)  # 戳记路径可读


def test_mirror_prune_removes_remote_deleted(tmp_path):
    m, s, remote, rd = _mat(tmp_path, flat=False)
    (remote / "zt_summary.json").write_text("{}")
    (remote / "per_temperature" / "1000K").mkdir(parents=True)
    (remote / "per_temperature" / "1000K" / "r.txt").write_text("x")
    (remote / "zt_old.png").write_bytes(b"x")
    _fetch(m)
    dest = rd / "step20_zt"
    assert (dest / "zt_old.png").is_file()
    (dest / "local_note.txt").write_text("mine")                 # 不是拉回来的：不删
    (remote / "zt_old.png").unlink()
    import shutil
    shutil.rmtree(remote / "per_temperature")
    _fetch(m)
    assert not (dest / "zt_old.png").exists()
    assert not (dest / "per_temperature").exists()
    assert (dest / "zt_summary.json").is_file() and (dest / "local_note.txt").is_file()


def test_fetch_extra_pure_mirror(tmp_path):
    m, s, remote, rd = _mat(tmp_path, flat=False, extra=["../../transport_plots"])
    (remote / "zt_summary.json").write_text("{}")
    tp = remote.parent.parent / "transport_plots"
    tp.mkdir()
    (tp / "Si_transport.png").write_bytes(b"new")
    local = rd / "transport_plots"
    local.mkdir(parents=True)
    (local / "Si_transport_data_1000K.csv").write_text("old")    # 旧版残留（还没有清单）
    _fetch(m)
    assert (local / "Si_transport.png").read_bytes() == b"new"
    assert not (local / "Si_transport_data_1000K.csv").exists()
    # 远端目录没了：什么都没拉到 -> 本地不删
    import shutil
    shutil.rmtree(tp)
    _fetch(m)
    assert (local / "Si_transport.png").is_file()


def test_zt_skill_declares_result_flat():
    from autozt import load_config, apply_skills
    cfg, _ = load_config(None)
    t = apply_skills(cfg)["task_types"]["zt-dft-cpu"]
    assert t["steps"][0].get("result_flat") is True
