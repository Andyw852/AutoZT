# -*- coding: utf-8 -*-
"""[kz=1-r2] S7.2 冻声子单点：2D 真空方向 kz 降到 1 + S8.21 的放行/拒绝 + 磁性解一致性核对。

  · odp_frozen.kpoints_mesh / reduce_vacuum_kz：只认自动网格；显式列表、本来就是 1 时如实报告
  · S7.2 真跑（沿用 test_0035_gen_pipeline 的玩具材料，S7 undeformed 改成 kz=3）：
      reduce -> 各构型 N1 N2 1、manifest 记 reduced；keep -> 不动；显式列表 -> 不动且不谎报
  · S8.21 check_inputs：只在 manifest 记了 reduced 时放行「面内相同、真空=1」；面内不同 / 旧 manifest 照旧报错
  · 磁矩：OSZICAR mag= / OUTCAR magnetization 表解析；换了磁性解（翻转、跳变）报错，1e-3 μB 级变化放行
"""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"
sys.path.insert(0, str(Path(__file__).resolve().parent))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def of():
    pytest.importorskip("numpy")
    return _load("of_kz1", KE / "odp_frozen.py")


AUTO3 = "K-Spacing Value to Generate K-Mesh: 0.030\n0\nGamma\n  6  42  3\n0.0  0.0  0.0\n"


# --------------------------------------------------------------------------- 网格解析 / 降 kz
def test_reduce_vacuum_kz(of):
    assert of.kpoints_mesh(AUTO3) == (6, 42, 3)
    txt, st, old, new = of.reduce_vacuum_kz(AUTO3, 2)
    assert st == "reduced" and old == [6, 42, 3] and new == [6, 42, 1]
    assert of.kpoints_mesh(txt) == (6, 42, 1)
    assert txt.splitlines()[0] == AUTO3.splitlines()[0] and txt.splitlines()[4] == "0.0  0.0  0.0"
    # 与 bd37be3 的写法逐字相同：已用旧版生成过的材料重新生成时不会被判成"输入变了"
    assert txt.splitlines()[3] == "  6  42  1"
    assert of.reduce_vacuum_kz(txt, 2)[1] == "already_1"
    assert of.reduce_vacuum_kz(AUTO3, 0)[3] == [1, 42, 3]                      # 按真空轴降，不写死第 3 个
    mp = "x\n0\nMonkhorst-Pack\n 4 4 2\n 0 0 0\n"
    assert of.reduce_vacuum_kz(mp, 2)[1] == "reduced"


@pytest.mark.parametrize("text", [
    "explicit\n3\nReciprocal\n0 0 0 1\n0.5 0 0 1\n0 0.5 0 1\n",              # 显式列表
    "fully automatic\n0\nAuto\n 20\n",                                        # Auto
    "",
])
def test_reduce_vacuum_kz_not_mesh(of, text):
    out, st, old, new = of.reduce_vacuum_kz(text, 2)
    assert st == "not_mesh" and out == text and old is None and of.kpoints_mesh(text) is None


# --------------------------------------------------------------------------- S7.2 真跑 + S8.21 核对
T35 = pytest.importorskip("test_0035_gen_pipeline")
mat = T35.mat                     # 复用 [0035] 的玩具材料夹具（S7 undeformed、step1、S5 目录、gen 脚本）


def _set_s7_kpoints(mat, text):
    (mat / "step7_deform" / "undeformed" / "KPOINTS").write_text(text)


def _gen(mat, monkeypatch, conf="[params]\n"):
    (mat / "step.conf").write_text(conf, encoding="utf-8")
    g, _f, _E = T35._gen72(mat, monkeypatch)
    g.main()
    out = mat / "step7c_odp_frozen"
    return out, json.loads((out / "odp_manifest.json").read_text())


def _reader(mat):
    r = _load("gen821_kz1_%d" % id(mat), mat / "gen_step12c_odp_read.py")
    r.of = _load("of_kz1_reader", mat / "odp_frozen.py")
    return r


def test_s72_reduces_kz_and_s821_accepts(mat, monkeypatch, capsys):
    _set_s7_kpoints(mat, "auto\n0\nGamma\n 12 12 3\n")
    out, man = _gen(mat, monkeypatch)
    for d in out.glob("odp-*"):
        assert (d / "KPOINTS").read_text().splitlines()[3].split() == ["12", "12", "1"]
    assert man["kz_vacuum_reduced"] is True
    assert man["kz_vacuum"]["axis"] == 2 and man["kz_vacuum"]["status"] == ["reduced"]
    assert man["kz_vacuum"]["s7_grid"] == [12, 12, 3] and man["kz_vacuum"]["odp_grid"] == [12, 12, 1]
    r = _reader(mat)
    r.check_inputs(mat, man)                                       # 放行
    assert "kz=1 优化）放行" in capsys.readouterr().out
    # 面内也变了：照旧报错
    _set_s7_kpoints(mat, "auto\n0\nGamma\n 15 15 3\n")
    with pytest.raises(SystemExit):
        r.check_inputs(mat, man)
    # 旧 manifest（没有 reduced 标记）：网格不同一律报错
    _set_s7_kpoints(mat, "auto\n0\nGamma\n 12 12 3\n")
    old = {k: v for k, v in man.items() if not k.startswith("kz_vacuum")}
    with pytest.raises(SystemExit):
        r.check_inputs(mat, old)


def test_s72_regen_is_stable(mat, monkeypatch):
    """降过一次再重新生成：KPOINTS 逐字不变 -> 已算完的点不归档。"""
    _set_s7_kpoints(mat, "auto\n0\nGamma\n 12 12 3\n")
    out, _ = _gen(mat, monkeypatch)
    done = out / "odp-ref" / "OUTCAR"
    done.write_text("General timing and accounting informations\n")
    _gen(mat, monkeypatch)
    assert done.is_file() and not list(out.glob("odp-ref/*stale*"))


def test_s72_keep_and_explicit_list(mat, monkeypatch, capsys):
    _set_s7_kpoints(mat, "auto\n0\nGamma\n 12 12 3\n")
    out, man = _gen(mat, monkeypatch, "[params]\nODP_KZ_VACUUM = keep\n")
    assert (out / "odp-ref" / "KPOINTS").read_text().splitlines()[3].split() == ["12", "12", "3"]
    assert man["kz_vacuum_reduced"] is False and man["kz_vacuum"]["mode"] == "keep"
    _reader(mat).check_inputs(mat, man)                            # 与 S7 逐字相同，照常通过
    lst = "explicit\n2\nReciprocal\n0 0 0 1\n0.5 0 0 1\n"
    _set_s7_kpoints(mat, lst)
    out, man = _gen(mat, monkeypatch)
    assert (out / "odp-ref" / "KPOINTS").read_text() == lst
    assert man["kz_vacuum_reduced"] is False and man["kz_vacuum"]["status"] == ["not_mesh"]
    assert "不是自动网格" in capsys.readouterr().out


def test_s72_bad_kz_mode(mat, monkeypatch):
    with pytest.raises(SystemExit):
        _gen(mat, monkeypatch, "[params]\nODP_KZ_VACUUM = half\n")


# --------------------------------------------------------------------------- 磁性解一致性
OUTCAR_TAIL = """ some header
 magnetization (x)

# of ion       s       p       d       tot
------------------------------------------
    1        0.010   0.020   {m1:.3f}   {m1t:.3f}
    2        0.010   0.020   {m2:.3f}   {m2t:.3f}
    3       -0.001  -0.010   0.000  -0.011
--------------------------------------------------
tot          0.019   0.030   6.000   6.050

 General timing and accounting informations
"""


def _cfg(d, total, m1, m2):
    d.mkdir(parents=True, exist_ok=True)
    (d / "OSZICAR").write_text("DAV: 60 ...\n   1 F= -.1E+03 E0= -.1E+03  d E =-.1E+03  mag=   %.4f\n" % total)
    (d / "OUTCAR").write_text(OUTCAR_TAIL.format(m1=m1, m1t=m1 + 0.03, m2=m2, m2t=m2 + 0.03))


def test_read_magnetism_and_consistency(tmp_path, monkeypatch, of):
    monkeypatch.syspath_prepend(str(KE))
    r = _load("gen821_mag", KE / "step7c_odp_frozen" / "odp_read" / "gen_step12c_odp_read.py")
    r.of = of
    src = tmp_path / "step7c_odp_frozen"
    _cfg(src / "odp-ref", 6.0, 3.0, 3.0)
    _cfg(src / "odp-m04p", 6.001, 3.0005, 2.9995)                 # 同一磁性解，~1e-3 μB
    _cfg(src / "odp-m04m", 5.999, 2.9995, 3.0005)
    m = r._read_magnetism(src / "odp-ref")
    assert m["total"] == [6.0] and m["sites"] == pytest.approx([3.03, 3.03, -0.011])
    rec = r.check_magnetism(src, ["odp-ref", "odp-m04p", "odp-m04m"])
    assert rec["checked"] is True
    _cfg(src / "odp-m05p", 0.0, 3.0, -3.0)                         # 一个 Cr 翻转：AFM
    with pytest.raises(SystemExit) as e:
        r.check_magnetism(src, ["odp-ref", "odp-m04p", "odp-m05p"])
    assert "odp-m05p" in str(e.value) and "odp-m04p" not in str(e.value)


def test_magnetism_skipped_when_not_spin_polarized(tmp_path, monkeypatch, of, capsys):
    monkeypatch.syspath_prepend(str(KE))
    r = _load("gen821_mag2", KE / "step7c_odp_frozen" / "odp_read" / "gen_step12c_odp_read.py")
    r.of = of
    src = tmp_path / "step7c_odp_frozen"
    (src / "odp-ref").mkdir(parents=True)
    (src / "odp-ref" / "OSZICAR").write_text("   1 F= -.1E+03 E0= -.1E+03  d E =-.1E+03\n")
    assert r.check_magnetism(src, ["odp-ref", "odp-m04p"]) == {"checked": False}
    assert "跳过磁性解一致性核对" in capsys.readouterr().out


def test_magnetic_mismatch_pure(of):
    ref = {"total": [6.0], "sites": [3.0, 3.0]}
    assert of.magnetic_mismatch(ref, {"total": [6.05], "sites": [3.04, 2.97]}) == []
    assert of.magnetic_mismatch(ref, {"total": [2.0], "sites": [3.0, 3.0]})            # 总磁矩跳变
    assert of.magnetic_mismatch(ref, {"total": [6.0], "sites": [4.0, 2.0]})            # 总量不变、原子间转移
    assert of.magnetic_mismatch(ref, {"total": [6.0], "sites": None}) == []            # 缺 OUTCAR 表：只比总量
