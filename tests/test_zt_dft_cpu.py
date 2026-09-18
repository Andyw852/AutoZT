# -*- coding: utf-8 -*-
"""zt-dft-cpu 自测：ZT 汇总步的纯逻辑（zt_common.py）+ 组合技能的装配完整性。

分两部分：
  [A] zt_common 物理/数值：张量约化口径、ZT 公式与单位、κ_L 插值（不外推）、
      n/p 分型取峰、AMSET/kl 产物解析。
  [B] 装配完整性：24 个步骤（含 electronic 开关两种模式）用 autozt 自身的
      find_asset / step_conf_sources 全量核对，任何缺件都会在这里被抓住
      —— 本技能靠软链装到上游技能目录，上游改了目录/步骤名时这里应当先红。

不依赖超算、不依赖 matplotlib。
"""
import importlib.util
import json
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SKILL = os.path.join(ROOT, "skill", "zt-dft-cpu")


def _load_zt_common():
    p = os.path.join(SKILL, "step20_zt", "zt_common.py")
    spec = importlib.util.spec_from_file_location("ztc_under_test", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


zc = _load_zt_common()


# ---------------------------------------------------------------- 造数工具
def _t33(xx, yy, zz):
    return [[xx, 0.0, 0.0], [0.0, yy, 0.0], [0.0, 0.0, zz]]


def _write_json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh)
    return str(path)


def _fake_transport(tmp_path, doping=(-1e20, 1e20), temps=(300.0, 800.0),
                    sigma=(1.0e5, 2.0e5), seebeck=(-200.0, -100.0),
                    kappa_e=(1.0, 2.0)):
    """按 AMSET 的 [i_doping][i_T][3][3] 结构造一份 transport.json。"""
    nd, nt = len(doping), len(temps)

    def _q(base):
        # 每个掺杂档一个矩阵：xx=base, yy=base, zz=base（各向同性，便于手算）
        return [[_t33(base[i], base[i], base[i]) for _ in range(nt)]
                for i in range(nd)]

    obj = {"doping": list(doping), "temperatures": list(temps),
           "conductivity": _q(sigma), "seebeck": _q(seebeck),
           "electronic_thermal_conductivity": _q(kappa_e)}
    return obj


def _fake_kappa(tmp_path, temps=(300.0, 800.0), k=(100.0, 50.0)):
    return {"KAPPA_DONE": True, "temperatures": list(temps),
            "kappa_xx_yy_zz": [[v, v, v] for v in k]}


# ---------------------------------------------------------------- [A] 纯逻辑
def test_diag_avg_caliber():
    """3D 取对角平均；2D 取面内 (xx+yy)/2 —— zz 必须被剔除。"""
    t = _t33(2.0, 4.0, 100.0)          # zz 故意给个离谱值
    assert zc.diag_avg(t, is_2d=False) == pytest.approx((2.0 + 4.0 + 100.0) / 3.0)
    assert zc.diag_avg(t, is_2d=True) == pytest.approx(3.0)      # (2+4)/2


def test_zt_formula_and_units():
    """ZT = S²σT/(κ_e+κ_L)，S 从 μV/K 换算到 V/K；与 ke 的 step8.3 同式。"""
    S_uVK, sigma, ke, kl, T = 200.0, 1.0e5, 3.0, 2.0, 800.0
    want = (200.0e-6) ** 2 * sigma * T / (ke + kl)
    assert zc.zt_value(S_uVK, sigma, ke, kl, T) == pytest.approx(want)
    # 单位自洽：若 S 忘乘 1e-6，结果会差 1e12 倍（量级断言，防单位回退）
    naive = S_uVK ** 2 * sigma * T / (ke + kl)
    assert naive / zc.zt_value(S_uVK, sigma, ke, kl, T) == pytest.approx(1e12)
    assert zc.zt_value(S_uVK, sigma, ke, kl, T) < 1.0
    # PF = S²σ [W/m/K²]
    assert zc.power_factor(S_uVK, sigma) == pytest.approx((200e-6) ** 2 * sigma)


def test_zt_invalid_inputs():
    assert zc.zt_value(None, 1e5, 1.0, 1.0, 300.0) is None
    assert zc.zt_value(200.0, 1e5, 1.0, -1.0, 300.0) is None       # 分母 <= 0
    assert zc.zt_value(float("nan"), 1e5, 1.0, 1.0, 300.0) is None


def test_interp_T_no_extrapolation():
    xs, ys = [100.0, 200.0, 300.0], [10.0, 20.0, 30.0]
    assert zc.interp_T(xs, ys, 150.0) == pytest.approx(15.0)
    assert zc.interp_T(xs, ys, 100.0) == pytest.approx(10.0)
    assert zc.interp_T(xs, ys, 99.0) is None        # 不外推（下界）
    assert zc.interp_T(xs, ys, 301.0) is None       # 不外推（上界）
    assert zc.interp_T(list(reversed(xs)), list(reversed(ys)), 250.0) == pytest.approx(25.0)


def test_kappa_scalar_caliber():
    kl = {"temperatures": [300.0], "k_xx": [10.0], "k_yy": [20.0], "k_zz": [1000.0]}
    assert zc.kappa_scalar(kl, is_2d=False)[0] == pytest.approx((10 + 20 + 1000) / 3)
    assert zc.kappa_scalar(kl, is_2d=True)[0] == pytest.approx(15.0)   # 面内，剔除 zz


def test_build_grid_2d_inplane_and_zt(tmp_path):
    """2D：ZT 必须用面内平均；κ_L 线性插值到电子段温度网格。"""
    tr = _fake_transport(tmp_path)
    tr["conductivity"] = [[_t33(1e5, 3e5, 0.0), _t33(1e5, 3e5, 0.0)],
                          [_t33(1e5, 3e5, 0.0), _t33(1e5, 3e5, 0.0)]]
    kl = _fake_kappa(tmp_path)
    tp = _load_zt_common().load_transport(_write_json(tmp_path / "t.json", tr))
    kp = zc.load_kappa(_write_json(tmp_path / "k.json", kl))
    g = zc.build_grid(tp, kp, is_2d=True)
    r = g["rows"][1]                                    # p 型档
    assert r["sigma_S/m"][0] == pytest.approx(2.0e5)     # (1e5+3e5)/2，zz 被剔除
    assert r["kappa_L_W/mK"][0] == pytest.approx(100.0)
    assert r["ZT"][0] == pytest.approx(
        (r["seebeck_uV/K"][0] * 1e-6) ** 2 * r["sigma_S/m"][0] * 300.0
        / (r["kappa_e_W/mK"][0] + 100.0))
    # 3D 口径下同一个输入会得到不同结果（证明口径真的生效）
    g3 = zc.build_grid(tp, kp, is_2d=False)
    assert g3["rows"][1]["sigma_S/m"][0] == pytest.approx((1e5 + 3e5 + 0.0) / 3)


def test_build_grid_zt_null_outside_kappa_range(tmp_path):
    """κ_L 只有 300/800 K：900 K 那一点必须记 None，而不是外推。"""
    tr = _fake_transport(tmp_path, temps=(300.0, 800.0, 900.0))
    kl = _fake_kappa(tmp_path, temps=(300.0, 800.0), k=(100.0, 50.0))
    tp = zc.load_transport(_write_json(tmp_path / "t.json", tr))
    kp = zc.load_kappa(_write_json(tmp_path / "k.json", kl))
    g = zc.build_grid(tp, kp, is_2d=False)
    assert g["rows"][0]["ZT"][0] is not None
    assert g["rows"][0]["ZT"][2] is None
    assert any("覆盖不到" in n for n in g["notes"])


def test_const300_mode(tmp_path):
    tr = _fake_transport(tmp_path, temps=(300.0, 800.0))
    kl = _fake_kappa(tmp_path, temps=(300.0, 800.0), k=(100.0, 50.0))
    tp = zc.load_transport(_write_json(tmp_path / "t.json", tr))
    kp = zc.load_kappa(_write_json(tmp_path / "k.json", kl))
    g = zc.build_grid(tp, kp, is_2d=False, ktemp_mode="const300")
    assert g["rows"][0]["kappa_L_W/mK"] == [100.0, 100.0]      # 两个温度都用 300 K 的值


def test_peak_zt_split_by_type(tmp_path):
    tr = _fake_transport(tmp_path, doping=(-1e21, -1e20, 1e20, 1e21),
                         temps=(300.0, 800.0),
                         sigma=(1e5, 2e5, 1e5, 2e5),
                         seebeck=(-50.0, -300.0, 50.0, 300.0),
                         kappa_e=(1.0, 1.0, 1.0, 1.0))
    kl = _fake_kappa(tmp_path)
    tp = zc.load_transport(_write_json(tmp_path / "t.json", tr))
    kp = zc.load_kappa(_write_json(tmp_path / "k.json", kl))
    peaks = zc.peak_zt(zc.build_grid(tp, kp, is_2d=False))
    for kind in ("n", "p"):
        p = peaks[kind]
        assert p is not None and p["ZT"] > 0
        assert (p["doping_cm-3"] < 0) == (kind == "n")
        assert p["T_K"] in (300.0, 800.0)
    # 掺杂越大 σ 越大、|S| 越小 → 峰值出现在两组掺杂里 ZT 更大的一档
    assert peaks["n"]["doping_cm-3"] in (-1e21, -1e20)


def test_load_transport_handles_monty_and_np_dict(tmp_path):
    """AMSET 有时把数组包成 monty 的 {@class: array, data: [...]}；n/p 也会分键。"""
    tr = _fake_transport(tmp_path)
    wrapped = {"doping": tr["doping"], "temperatures": tr["temperatures"],
               "conductivity": {"@class": "array", "data": tr["conductivity"]},
               "seebeck": {"n": tr["seebeck"][:1], "p": tr["seebeck"][1:]},
               "electronic_thermal_conductivity": tr["electronic_thermal_conductivity"]}
    got = zc.load_transport(_write_json(tmp_path / "t.json", wrapped))
    assert len(got["doping"]) == 2 and len(got["seebeck"]) == 2


def test_load_kappa_requires_fields(tmp_path):
    """既没有惰温数组、也没有 300K 单值 → 必须报错（不能静默出空结果）。"""
    p = _write_json(tmp_path / "bad.json", {"KAPPA_DONE": True, "temperatures": [300.0]})
    with pytest.raises(ValueError):
        zc.load_kappa(p)


def test_load_kappa_accepts_single_temperature_fallback(tmp_path):
    """kl-mace-* 老格式只有 kappa_300K_xx_yy_zz：降级成单温度点并打标。"""
    p = _write_json(tmp_path / "macelike.json",
                    {"KAPPA_DONE": True, "kappa_300K_xx_yy_zz": [110.0, 90.0, 5.0]})
    k = zc.load_kappa(p)
    assert k["single_T_fallback"] is True
    assert k["temperatures"] == [300.0] and k["k_xx"] == [110.0]

    # 单点 κ_L + const300 模式 → 仍能给出全温 ZT 曲线
    tr = _fake_transport(tmp_path, temps=(300.0, 800.0))
    tp = zc.load_transport(_write_json(tmp_path / "t.json", tr))
    g = zc.build_grid(tp, k, is_2d=False, ktemp_mode="const300")
    _scalar = (110.0 + 90.0 + 5.0) / 3.0      # 3D 对角平均
    assert g["rows"][1]["kappa_L_W/mK"] == pytest.approx([_scalar, _scalar])
    assert any("只有 300 K 单值" in n for n in g["notes"])

    # interp 模式下：只有 300 K 落在 κ_L 温区内，800 K 记 null（不外推）
    g2 = zc.build_grid(tp, k, is_2d=False, ktemp_mode="interp")
    assert g2["rows"][1]["ZT"][0] is not None and g2["rows"][1]["ZT"][1] is None


def test_text_report_contains_key_sections(tmp_path):
    tr = _fake_transport(tmp_path)
    kl = _fake_kappa(tmp_path)
    tp = zc.load_transport(_write_json(tmp_path / "t.json", tr))
    kp = zc.load_kappa(_write_json(tmp_path / "k.json", kl))
    g = zc.build_grid(tp, kp, is_2d=False)
    txt = zc.text_report(g, zc.peak_zt(g), {"material": "X", "transport_src": "a",
                                            "kappa_src": "b", "dim": "3d",
                                            "dim_note": "n"})
    assert "峰值 ZT" in txt and "ZT = S^2" in txt and "3d" in txt


# -------------------------------------------- [C] gen 脚本级：温区不重叠必须报错
def _gen_fixture(root, tr, kl, dim="3d"):
    """搭一个最小的远端技能目录：step8_amset/transport.json + step6_kappa/"""
    (root / "step8_amset").mkdir(parents=True)
    (root / "step6_kappa").mkdir(parents=True)
    (root / "step1_opt").mkdir(parents=True)
    (root / "step8_amset" / "transport.json").write_text(json.dumps(tr), encoding="utf-8")
    (root / "step6_kappa" / "kappa_summary.json").write_text(json.dumps(kl), encoding="utf-8")
    (root / "step1_opt" / "workflow_method.txt").write_text("DIM = %s\n" % dim,
                                                             encoding="utf-8")
    src = os.path.join(SKILL, "step20_zt")
    for f in ("gen_step20_zt.py", "zt_common.py", "step.conf"):
        shutil.copyfile(os.path.join(src, f), root / f)


def test_gen_step_errors_when_kappa_range_disjoint(tmp_path):
    """κ_L 温区与电子段温度网格毫无重叠时：必须非 0 退出，且不能写出空结果。

    （否则步骤会被判完成、下游拿到一张全 null 的 ZT 表 —— 静默错比报错危险得多。）"""
    tr = {"doping": [-1e20, 1e20], "temperatures": [300.0, 400.0],
          "conductivity": [[_t33(1e5, 1e5, 1e5)] * 2] * 2,
          "seebeck": [[_t33(-200.0, -200.0, -200.0)] * 2] * 2,
          "electronic_thermal_conductivity": [[_t33(1.0, 1.0, 1.0)] * 2] * 2}
    kl = {"KAPPA_DONE": True, "temperatures": [1000.0, 1100.0],
          "kappa_xx_yy_zz": [[50.0, 50.0, 50.0], [40.0, 40.0, 40.0]]}
    root = tmp_path / "Si" / "zt-dft-cpu"
    _gen_fixture(root, tr, kl)
    r = subprocess.run([sys.executable, "gen_step20_zt.py", "--material", "Si"],
                       cwd=str(root), capture_output=True, text=True)
    assert r.returncode != 0
    assert "没有任何温度点能算出 ZT" in (r.stdout + r.stderr)
    # done_marker 是 zt_summary.json：它不能出现，否则步骤会被判"完成"
    assert not (root / "step20_zt" / "zt_summary.json").exists()


def test_gen_step_ok_on_overlapping_range(tmp_path):
    """温区重叠时：正常退出、写出 zt_summary.json 与 txt。"""
    tr = {"doping": [-1e20, 1e20], "temperatures": [300.0, 400.0],
          "conductivity": [[_t33(1e5, 1e5, 1e5)] * 2] * 2,
          "seebeck": [[_t33(-200.0, -200.0, -200.0)] * 2] * 2,
          "electronic_thermal_conductivity": [[_t33(1.0, 1.0, 1.0)] * 2] * 2}
    kl = {"KAPPA_DONE": True, "temperatures": [200.0, 500.0],
          "kappa_xx_yy_zz": [[50.0, 50.0, 50.0], [40.0, 40.0, 40.0]]}
    root = tmp_path / "Si" / "zt-dft-cpu"
    _gen_fixture(root, tr, kl)
    r = subprocess.run([sys.executable, "gen_step20_zt.py", "--material", "Si"],
                       cwd=str(root), capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    out = root / "step20_zt"
    assert (out / "zt_summary.json").is_file() and (out / "zt_summary.txt").is_file()
    d = json.loads((out / "zt_summary.json").read_text(encoding="utf-8"))
    assert d["ZT_DONE"] is True and d["material"] == "Si"
    assert d["peak_ZT"]["n"] is not None


# ------------------------------------------------------------- [B] 装配完整性
def _expand(electronic, lattice=True):
    from autozt import load_config, apply_skills, expand_optional_steps, _seq_sort_steps
    cfg, _ = load_config(None)
    cfg = apply_skills(cfg)
    t = dict(cfg["task_types"]["zt-dft-cpu"])
    t["key"] = "zt-dft-cpu"
    t["electronic"] = electronic
    t["lattice"] = lattice
    expand_optional_steps(t)
    _seq_sort_steps(t.get("steps") or [])
    seg = {"steps_cfg": t.get("steps"), "skill_dir": t.get("skill_dir"),
           "template_dir": t.get("template_dir"),
           "template_layout": t.get("template_layout")}
    m = {"name": "X", "hpc_name": "jzzn", "_seg": seg, "template_map": {}, "ps": {}}
    return cfg, t, m


@pytest.mark.parametrize("electronic,lattice,min_steps",
                         [(True, True, 24),      # 完整 ZT 全流程
                          (True, False, 17),     # 只电子段 + 汇总
                          (False, True, 8),      # 只晶格段 + 汇总
                          (False, False, 1)])    # 只出 ZT 汇总（两段都借兄弟技能结果）
def test_all_step_assets_resolve(electronic, lattice, min_steps):
    """四种模式下，gen 脚本 / 模板 / gen_need / step.conf 都必须找得到。"""
    from autozt import step_cfg, find_asset, step_conf_sources
    cfg, t, m = _expand(electronic, lattice)
    assert len(t["steps"]) == min_steps
    missing = []
    for s in t["steps"]:
        n = s["name"]
        sc = step_cfg(t, n, m)
        gen = (sc.get("gen") or "").split()[0]
        for f in [gen] + list(sc.get("gen_need") or []):
            if f and f != "step.conf" and not find_asset(cfg, t, m, f, n):
                missing.append((n, f))
        if not step_conf_sources(cfg, t, m, n):
            missing.append((n, "step.conf"))
    assert missing == []


def test_segment_groups_drop_only_their_own_steps():
    """两个段组各自只影响自己那一段；S20_zt 的缺失依赖必须被忽略。"""
    from autozt import _dag_needs
    _, t_on, _ = _expand(True, True)
    full = {s["name"] for s in t_on["steps"]}

    _, t_e, m_e = _expand(False, True)          # 关电子段
    e_off = {s["name"] for s in t_e["steps"]}
    assert e_off < full and "step8_amset" in full - e_off
    assert "step6_kappa" in e_off

    _, t_l, m_l = _expand(True, False)          # 关晶格段
    l_off = {s["name"] for s in t_l["steps"]}
    assert l_off < full and "step6_kappa" in full - l_off
    assert "step8_amset" in l_off

    # 关电子段：step8_amset 不在步骤表 → autozt 忽略它，只等 step6_kappa
    s20 = [s for s in t_e["steps"] if s["name"] == "step20_zt"][0]
    deps = _dag_needs(t_e, m_e, s20, "step6_kappa")
    assert deps == ["step8_amset", "step6_kappa"]
    assert [d for d in deps if d in e_off] == ["step6_kappa"]

    # 两段都关：只剩汇总步，且无有效依赖（立即就绪）
    _, t_z, m_z = _expand(False, False)
    assert [s["name"] for s in t_z["steps"]] == ["step20_zt"]
    s20z = t_z["steps"][0]
    assert [d for d in _dag_needs(t_z, m_z, s20z, None) if d in {"step20_zt"}] == []


def test_skill_dir_is_symlinks_not_copies():
    """装配纪律：除 skill.yaml / step.conf / step20_zt 外，技能目录里必须都是软链
    （防止有人图省事把上游脚本拷进来，产生副本漂移）。"""
    real = []
    for name in os.listdir(SKILL):
        if name in ("skill.yaml", "step.conf", "step20_zt", "README.md",
                    "__pycache__", "tests"):
            continue
        p = os.path.join(SKILL, name)
        if not os.path.islink(p):
            real.append(name)
        else:
            assert os.path.exists(p), "软链失效：%s" % name
    assert real == []
