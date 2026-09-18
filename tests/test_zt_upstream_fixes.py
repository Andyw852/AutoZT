# -*- coding: utf-8 -*-
"""zt-dft-cpu 的 Si 全流程实跑中暴露的 4 处**上游**缺陷 —— 回归钉子。

这些缺陷都不在 zt-dft-cpu 自己的代码里，而在 ke-dft-cpu / kl-dft-cpu 侧。
本文件把"修好之后必须成立"的行为/结构钉住，防止被悄悄改回去。
每条都对应 tmp/zt_upstream_patches_20260918/ 里的一份补丁，并都做过集群实跑验证：

  01 kl-dft-cpu 抽帧越界      -> 实跑 -j SK5_fc init 由 IndexError 变 gen 完成
  02 ke S3 护栏参数不可配置    -> 实跑 Si 网格 34^3 -> 26^3 并 gen 完成
  03 ke HSE 画图缺 import os  -> 实跑 -j S2.3_hseplot start 出 band_summary.json
  04 ke AMSET 非极性介电误拦    -> 实跑 -j S8_kappa_e start gen 完成并提交作业

标 [行为] 的是真的调用函数；标 [结构] 的只钉住源码里的关键行（03/04 要造
DFPT/AMSET 输入才能行为级复现，代价不划算，集群实跑已覆盖）。
"""
import importlib.util
import os
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(path, name):
    """按 autozt 推送 gen_need 的方式补齐搜索路径（脚本目录 + 技能根 + 公共池）。"""
    for d in (os.path.dirname(path), os.path.dirname(os.path.dirname(path)),
              os.path.join(ROOT, "skill", "_common", "opt"),
              os.path.join(ROOT, "skill", "_common")):
        if d not in sys.path:
            sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    mod.__dict__["__file__"] = path
    mod.__dict__["__name__"] = name
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- 01 [行为]
def test_kl_frame_check_skips_equilibrium_frame(tmp_path):
    """S4 会额外生成平衡帧 disp-00000；抽帧校验必须只数位移帧。

    修复前：目录 11 个（含平衡帧）而 displacements 只有 10 条 ->
    disps[len(frames)-1] 抛 IndexError，S5_fc 的 gen 直接失败。
    这里不造 vasprun.xml（每帧都会被 continue 跳过），所以只验证**不再越界**。
    """
    kl = _load(os.path.join(ROOT, "skill", "kl-dft-cpu", "kl_common.py"), "klc_pin")
    d4 = tmp_path / "step4_disp"
    d4.mkdir()
    n_disp = 10
    (d4 / "phono3py_disp.yaml").write_text(yaml.safe_dump({
        "supercell": {"lattice": [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]],
                      "points": [{"symbol": "Si", "coordinates": [0.0, 0.0, 0.0],
                                  "mass": 28.0855}]},
        "dataset": {"displacements": [[[0.0, 0.0, 0.0]] for _ in range(n_disp)]},
    }), encoding="utf-8")
    for i in range(n_disp + 1):          # disp-00000（平衡帧）+ disp-00001..10
        (d4 / ("disp-%05d" % i)).mkdir()
    ok, note = kl.check_frames_match_displacements(d4)
    assert ok is True, note

    # 反向钉子：真要越界时必须被防御住（目录数 > 位移数、且带 vasprun 也不该崩）
    for i in range(n_disp + 1):
        (d4 / ("disp-%05d" % i) / "vasprun.xml").write_text("<modeling/>", encoding="utf-8")
    ok2, note2 = kl.check_frames_match_displacements(d4)
    assert ok2 in (True, False), note2     # 允许判定失败，但**不允许抛异常**
    assert isinstance(note2, str) and note2


# ---------------------------------------------------------------- 02 [结构]
def test_ke_uniform_guard_params_are_step_conf_keys():
    """ke S3 护栏参数必须能由 step.conf 覆盖（Si 靠 DK_MAX_3D=0.08 才过护栏）。"""
    p = os.path.join(ROOT, "skill", "ke-dft-cpu", "step3_uniform", "gen_step5_uniform.py")
    src = open(p, encoding="utf-8").read()
    for key in ('"DK_MAX"', '"DK_MAX_2D"', '"DK_MAX_3D"', '"UNIFORM_NMAX"'):
        assert key in src, "SPEC 里缺 %s（补丁 02 被回退？）" % key
    assert "global DK_MAX, DK_MAX_2D, DK_MAX_3D, UNIFORM_NMAX" in src
    mod = _load(p, "ke_uniform_pin")   # 能被导入 = 语法/顶层没问题
    for key in ("DK_MAX", "DK_MAX_2D", "DK_MAX_3D", "UNIFORM_NMAX"):
        assert key in mod.SPEC, "SPEC 里缺 %s" % key


# ---------------------------------------------------------------- 03 [结构]
def test_ke_hse_plot_imports_os():
    """step2.3_hse_plot 的 gen 用了 os.path.join，必须 import os（否则 NameError）。"""
    p = os.path.join(ROOT, "skill", "ke-dft-cpu", "step2_bandgap",
                     "step2.3_hse_plot", "gen_step4.1_plot_band.py")
    src = open(p, encoding="utf-8").read()
    assert "\nimport os\n" in src, "缺 import os（补丁 03 被回退？）"
    assert "os.path.join" in src          # 确实用到了，才有意义


# ---------------------------------------------------------------- 04 [结构]
def test_ke_amset_nonpolar_dielectric_not_blocked():
    """非极性体系 eps_static 恒等于 eps_inf 是物理正确结果，不能当成 DFPT 失效拦下。"""
    p = os.path.join(ROOT, "skill", "ke-dft-cpu", "step8_amset", "gen_step10_amset.py")
    src = open(p, encoding="utf-8").read()
    anchor = "if max(abs(c) for c in _coup) < 1e-9:"
    i = src.find(anchor)
    assert i > 0, "找不到 Frohlich 耦合判据"
    block = src[i: i + 260]
    assert "_is_nonpolar()" in block, \
        "eps_static 恒等于 eps_inf 的分支没有先问 _is_nonpolar()（补丁 04 被回退？）"
