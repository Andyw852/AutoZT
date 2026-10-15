"""[V124-fix] 插值 factor 支持小数：3D/2D 全网格 vasprun 按 IBZ 等价值换算。

背景（2026-10-15 复核）：AMSET 的 `interpolation_factor` 乘的是 **vasprun 的 k 点数**
（equivalence 数 = nk × factor，amset/interpolation/bandstructure.py:105）。
WAVEFUNCTION_FULL / S3b 的 vasprun 是整张 Γ/MP 全网格，点数比 IBZ 大十几到几十倍，
同一个出厂 factor 在全网格上会给出大 (n_full/n_ir)^(1/3) 倍的网格。

Si 34^3 全网格（39304 点）的 IBZ 是 1059 点（37 倍）：f0=10 换算后 f≈0.269，
旧代码 `round()` + `max(1, ·)` 把它钳到 **1** -> 等价值 10590 变 39304，
插值网格从 109^3 静默涨到 167^3（9 倍点、OOM 风险）。本文件锁住修复后的契约。

不依赖超算。纯算术用例任何环境都跑；真实网格用例需要 pymatgen + BoltzTraP2（缺则跳过）。
"""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
KE = ROOT / "skill" / "ke-dft-cpu"

# Si 34^3 全网格 / 它的 IBZ 点数（真实体系，见 tmp/amset2d/si_full/）
N_FULL, N_IR, F0 = 39304, 1059, 10.0
EQUIV = 10590          # = n_ir * f0 = 1059 * 10，两边换算后必须都等于它


def _load(name, path):
    for d in (ROOT, KE, ROOT / "skill" / "_common" / "opt"):
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _g10():
    return _load("g10_0124", KE / "step8_amset" / "gen_step10_amset.py")


def _g14():
    return _load("g14_0124", KE / "step8.4_amset2d" / "gen_step14_amset2d.py")


# ---------------------------------------------------------------- SPEC / 落盘
@pytest.mark.parametrize("which", ["3d", "2d"])
def test_spec_accepts_float_factor(which):
    """[V124-fix] SPEC 必须是 float —— 写 "int" 时 step.conf 写 0.269 会被 int() 截成 0。"""
    g = _g10() if which == "3d" else _g14()
    assert g.SPEC["INTERPOLATION_FACTOR"] == (None, "float")


@pytest.mark.parametrize("which", ["3d", "2d"])
def test_fnum_keeps_integer_as_int(which):
    """IBZ 路径（f0 = 10 / 61 这类整数）落盘仍是 int —— 现有输出逐字不变。"""
    g = _g10() if which == "3d" else _g14()
    for v in (10, 10.0, 61, 4.0):
        x = g._fnum(v)
        assert isinstance(x, int) and not isinstance(x, bool), (v, x)
        assert x == int(v)


@pytest.mark.parametrize("which", ["3d", "2d"])
def test_fnum_keeps_small_float(which):
    """小数 factor 不能被截断成 0（旧 int(f) 的行为）—— 0 会让 nk*factor 归零。"""
    g = _g10() if which == "3d" else _g14()
    f_start = F0 * (N_IR / N_FULL)
    assert 0.0 < f_start < 1.0
    assert g._fnum(f_start) == pytest.approx(f_start, abs=1e-6)
    assert g._fnum(f_start) > 0


# ---------------------------------------------------------------- 换算契约
@pytest.mark.parametrize("which", ["3d", "2d"])
def test_full_grid_conversion_preserves_equivalence(which):
    """全网格臂的换算必须保持「等价值 nk × factor」不变。

    apply_mesh_min 里的换算：f_start = f0 × n_ir / n_full，只保留 nk*factor >= 1 的保底。
    旧代码 round + max(1, ·) 会把 0.269 钳到 1 -> 等价值 39304（网格 167^3）。
    """
    g = _g10() if which == "3d" else _g14()
    f_start = F0 * (N_IR / N_FULL)
    assert f_start < 1.0                                    # 钳到 1 就是 bug 现场
    assert round(N_FULL * f_start) == round(N_IR * F0) == EQUIV
    floor = 1.0 / max(1, int(N_FULL))                       # 保底：nk*factor >= 1
    assert max(f_start, floor) == f_start                   # 保底不能顶掉正常值
    assert round(N_FULL * max(f_start, floor)) == EQUIV


def test_settings_line_uses_g_not_d(tmp_path, monkeypatch):
    """写 settings.yaml 时 factor 用 %g（旧 %d + int() 会把 0.269 写成 0）。"""
    g = _g10()
    monkeypatch.setattr(g, "check_carriers_per_atom", lambda *a, **k: None)
    monkeypatch.setattr(g, "_cell_volume_cm3", lambda *a, **k: 1e-22)
    monkeypatch.setattr(g, "_count_atoms", lambda *a, **k: 1)
    monkeypatch.setattr(g, "_add_defect_doping", lambda dop, *a, **k: (dop, None))
    monkeypatch.setattr(g, "read_pop_frequency", lambda *a, **k: None)
    monkeypatch.setattr(g, "NWORKERS", 6)                  # 出厂 None = 用满节点核数

    def _line(factor):
        monkeypatch.setattr(g, "INTERPOLATION_FACTOR", factor)
        g.write_settings(tmp_path, None, None, None, None)
        txt = (tmp_path / "settings.yaml").read_text(encoding="utf-8")
        return [ln for ln in txt.splitlines() if ln.startswith("interpolation_factor:")]

    assert _line(F0 * (N_IR / N_FULL)) == ["interpolation_factor: 0.269438"]
    assert _line(10) == ["interpolation_factor: 10"]        # 整数路径逐字不变
    assert _line(None) == []                                # 没给 factor 就不写这一行


# ---------------------------------------------------------------- 真实网格
def _si():
    pytest.importorskip("pymatgen")
    pytest.importorskip("BoltzTraP2")
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice([[0, 2.715, 2.715], [2.715, 0, 2.715], [2.715, 2.715, 0]]),
                     ["Si", "Si"], [[0, 0, 0], [0.25, 0.25, 0.25]])


@pytest.mark.parametrize("which", ["3d", "2d"])
def test_si_34cube_full_grid_and_ibz_land_on_the_same_mesh(which):
    """★ 核心契约：Si 34^3 全网格臂（f≈0.269）与 IBZ 臂（f=10）必须给出同一张 109^3。

    网格只由等价值 nk×factor 决定（_interp_mesh -> BoltzTraP2 sphere.get_equivalences）。
    """
    g = _g10() if which == "3d" else _g14()
    si = _si()
    m_full = g._interp_mesh(si, N_FULL, F0 * (N_IR / N_FULL))
    m_ir = g._interp_mesh(si, N_IR, F0)
    assert list(m_full) == list(m_ir) == [109, 109, 109]


def test_mesh_depends_only_on_equivalence_value():
    """网格是 nk×factor 的函数：等价值相同 -> 网格相同（换算成立的前提）。"""
    g = _g10()
    si = _si()
    a = g._interp_mesh(si, N_FULL, F0 * (N_IR / N_FULL))
    b = g._interp_mesh(si, N_IR, F0)
    c = g._interp_mesh(si, 2 * N_IR, F0 / 2)
    assert list(a) == list(b) == list(c)


def test_mesh_min_float_factor_end_to_end(tmp_path, monkeypatch):
    """apply_mesh_min 全链路：全网格 vasprun + f0=10 -> factor 落成小数、网格 109^3。

    同时覆盖上限循环的小数分支（MESH_MAX 设小 -> 必须按比例降，而不是 np.floor 到 0）。
    """
    g = _g10()
    si = _si()

    class _VR:
        actual_kpoints = [[0.0, 0.0, 0.0]] * N_FULL
        final_structure = si

    import pymatgen.io.vasp.outputs as _pvo          # noqa: PLC0415  (importorskip 已确保存在)
    monkeypatch.setattr(_pvo, "Vasprun", lambda *a, **k: _VR())
    monkeypatch.setattr(g.kc, "full_grid_ibz_count",
                        lambda st, kpts: (N_IR, N_FULL, (34, 34, 34)))
    monkeypatch.setattr(g, "INTERPOLATION_FACTOR", F0)
    monkeypatch.setattr(g, "INTERPOLATION_FACTOR_EXPLICIT", False)
    # 下限设 109 会让 need=[109,109,109] 正好被 f0 换算后的网格满足（不触发自动升 factor 循环）；
    # 全 0 则函数会在入口 `if not lo and not lo_z and not MESH_MAX: return None` 直接返回。
    monkeypatch.setattr(g, "MESH_MIN", 109)
    monkeypatch.setattr(g, "MESH_MIN_KZ", 0)
    monkeypatch.setattr(g, "MESH_MAX", None)

    f, mesh = g.apply_mesh_min("vasprun.xml", tmp_path)
    assert list(mesh) == [109, 109, 109]
    assert 0.0 < float(f) < 1.0                       # 换算后是小数，不是被钳到的 1
    assert round(N_FULL * float(f)) == EQUIV
    info = (tmp_path / "interpolation_info.json").read_text(encoding="utf-8")
    assert '"factor_used": 0.269438' in info          # 落盘存小数（不是 int）
    assert '"nk_ir": 1059' in info and '"full_grid": true' in info

    # 上限循环：MESH_MAX 远小于 109 -> 必须按比例降 factor（旧 int(np.floor(...)) 会变 0）
    monkeypatch.setattr(g, "INTERPOLATION_FACTOR", F0)
    monkeypatch.setattr(g, "MESH_MIN", 11)            # 小下限，避免降完又被升 factor 循环推回去
    monkeypatch.setattr(g, "MESH_MAX", 60)
    f2, mesh2 = g.apply_mesh_min("vasprun.xml", tmp_path)
    assert max(mesh2) <= 60
    assert 0.0 < float(f2) < float(f)                 # 降下来但仍是正小数
    assert (tmp_path / "interpolation_info.json").exists()
