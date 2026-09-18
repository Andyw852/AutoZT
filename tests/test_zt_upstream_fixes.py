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
import json
import os
import subprocess
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


# ---------------------------------------------------------------- 05 [行为]
def _load_exec(path, name, extra_dirs=()):
    import importlib.util
    for d in (os.path.dirname(path),) + tuple(extra_dirs):
        if d not in sys.path:
            sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    mod.__dict__["__file__"] = path
    mod.__dict__["__name__"] = name
    spec.loader.exec_module(mod)
    return mod


def test_overlap_preflight_skips_bandwindow_when_runlog_absent(tmp_path):
    """修复前：拿"AMSET 插值带窗口"的兜底默认 11-17 与 step4_wave 的真实窗口比 ->
    任何窗口 != 11-17 的材料都会被拦死（Si 实测 2-6）。运行日志缺失时必须判"未知"并跳过。"""
    p = os.path.join(ROOT, "skill", "ke-dft-cpu", "step8.4_amset2d", "overlap_preflight.py")
    pf = _load_exec(p, "pf_pin")
    (tmp_path / "step8_amset").mkdir()
    (tmp_path / "step4_wave").mkdir()
    # AMSET 运行日志故意不存在；amset wave 日志写真实窗口 2—6（em dash）
    (tmp_path / "step4_wave" / "amset.log").write_text("Including bands 2\u20146\n",
                                                      encoding="utf-8")
    (tmp_path / "step8_amset" / "settings.yaml").write_text(
        "doping: [-1e20, 1e20]\ntemperatures: [100, 300]\nscattering_type: [ADP, IMP]\n",
        encoding="utf-8")
    assert pf.run_band_window(tmp_path / "step8_amset") is None      # 未知，不是 (11, 17)
    verdict, lines = pf.run(tmp_path, tmp_path / "step8_amset", False)
    assert verdict != "error", "运行日志缺失时仍然拦下了：%s" % lines
    assert any("运行日志缺失" in l for l in lines)


# ---------------------------------------------------------------- 07 [行为]
def _vd_fixture(root, symbols, eps_inf, eps_ion):
    root.mkdir(parents=True)
    (root / "POSCAR").write_text("\n".join([
        "c", "1.0", "5 0 0", "0 5 0", "0 0 5", symbols, "1", "Direct", "0 0 0"]) + "\n",
        encoding="utf-8")
    L = [" OUTCAR fixture", "",
         " MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"]
    for r in range(3):
        L.append("           " + " ".join("%.6f" % eps_inf[r][c] for c in range(3)))
    L += ["", " MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION"]
    for r in range(3):
        L.append("           " + " ".join("%.6f" % eps_ion[r][c] for c in range(3)))
    (root / "OUTCAR").write_text("\n".join(L) + "\n", encoding="utf-8")


def test_dielectric_validator_exempts_nonpolar_but_not_polar(tmp_path):
    """ε_static ≡ ε_inf：单元素（非极性）是物理正确 -> 不判 FAIL；极性体系仍判 FAIL。"""
    script = os.path.join(ROOT, "skill", "ke-dft-cpu", "step5_dielect",
                          "validate_dielectric.py")
    diag = [[13.0, 0, 0], [0, 13.0, 0], [0, 0, 13.0]]
    zero = [[0.0] * 3] * 3
    d1 = tmp_path / "step5_dielect_nonpolar"
    _vd_fixture(d1, "Si", diag, zero)
    d2 = tmp_path / "step5_dielect_polar"
    _vd_fixture(d2, "Na Cl", diag, zero)
    for d in (d1, d2):
        r = subprocess.run([sys.executable, script, "--step-dir", str(d)],
                           capture_output=True, text=True)
        assert (d / "dielectric_check.json").is_file(), r.stderr
    non = json.loads((d1 / "dielectric_check.json").read_text(encoding="utf-8"))
    pol = json.loads((d2 / "dielectric_check.json").read_text(encoding="utf-8"))
    assert non["ok"] is True and non.get("notes"), non
    assert pol["ok"] is False and pol["reasons"], pol


# ---------------------------------------------------------------- 06 [逻辑]
def test_discriminant_done_marker_resolves_to_script_output():
    """S2.155 的 done_marker 必须解析到 decide_discriminant.py 真正写出的位置。

    脚本写在 step2_bandgap/step2.15_discriminant/，而该步目录是
    step2_bandgap/step2.155_discriminant_decide —— 清单里必须用 .. 才找得到
    （ck_plot/do_run_gen_step 都是 os.path.join(step_dir, marker)）。"""
    sk = yaml.safe_load(open(os.path.join(ROOT, "skill", "ke-dft-cpu", "skill.yaml"),
                             encoding="utf-8"))
    defs = list(sk["steps"])
    for grp in (sk.get("optional_steps") or {}).values():
        defs += list((grp or {}).get("steps") or [])
    step = next(s for s in defs
                if s.get("name") == "step2_bandgap/step2.155_discriminant_decide")
    marker = step["done_marker"]
    resolved = os.path.normpath(
        os.path.join("step2_bandgap/step2.155_discriminant_decide", marker))
    expected = os.path.normpath("step2_bandgap/step2.15_discriminant/discriminant.json")
    assert resolved == expected, \
        "done_marker 没指到脚本真实写出的位置：%s" % resolved

    # 本技能（zt-dft-cpu）里那份声明也要一致
    zs = yaml.safe_load(open(os.path.join(ROOT, "skill", "zt-dft-cpu", "skill.yaml"),
                             encoding="utf-8"))
    zdefs = list(zs["steps"])
    for grp in (zs.get("optional_steps") or {}).values():
        zdefs += list((grp or {}).get("steps") or [])
    zstep = next(s for s in zdefs
                 if s.get("name") == "step2_bandgap/step2.155_discriminant_decide")
    assert zstep["done_marker"] == marker


# ---------------------------------------------------------------- 08/09/10 [结构]
def test_gen_steps_marker_lands_in_their_own_step_dir():
    """两个 run:gen 步的 done_marker 必须真的会出现在【本步目录】里：

    - S5.1：脚本默认写 step5_dielect/dielectric_check.json（DFPT 步的目录）→ 清单里必须
      用 --json 指到本步目录，且脚本要会建父目录；
    - S2.155：脚本只写 step2.15_discriminant/，清单里用 ../ 取它 → 脚本必须补建本步目录，
      否则路径里的 .. 解析不了（test -f / isfile 都要求中间目录存在）。
    """
    src = open(os.path.join(ROOT, "skill", "ke-dft-cpu", "step5_dielect",
                            "validate_dielectric.py"), encoding="utf-8").read()
    assert "out.parent.mkdir(parents=True, exist_ok=True)" in src

    dsrc = open(os.path.join(ROOT, "skill", "ke-dft-cpu", "step2_bandgap",
                             "step2.15_discriminant", "decide_discriminant.py"),
                encoding="utf-8").read()
    assert "step2.155_discriminant_decide" in dsrc and "mkdir" in dsrc

    for skill in ("ke-dft-cpu", "zt-dft-cpu"):
        sk = yaml.safe_load(open(os.path.join(ROOT, "skill", skill, "skill.yaml"),
                                 encoding="utf-8"))
        defs = list(sk["steps"])
        for grp in (sk.get("optional_steps") or {}).values():
            defs += list((grp or {}).get("steps") or [])
        st = next(s for s in defs if s.get("name") == "step5_dielect_validate")
        assert "--json step5_dielect_validate/dielectric_check.json" in st["gen"], \
            "%s 的 S5.1 gen 没把 json 指到本步目录：%s" % (skill, st["gen"])
