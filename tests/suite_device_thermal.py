#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""device-thermal / device-common 回归套件（纯本机、几秒~1 分钟，不碰集群）。

覆盖用户 review 列出的会改变结果/流程的问题：
  1. 侧壁边界：恒温只搭沟道，氧化层侧壁绝热；与独立 FEM 交叉验证；
  2. 参数流向：S2 上 TBC_OVERRIDE/CHANNEL_THICKNESS 生效；S3 回退 S1 的网格；
  3. 能量守恒、绝热极限对精确 1D、x/y 网格收敛 + Richardson；
  4. 各向异性氧化层（hBN）R 用 ky；CHANNEL_THICKNESS 换算保持 κ·t；
  5. KAPPA_SOURCE 三取值 × 上游四形态（normalized / factor / 都没有 / runs[]）；
  6. 非法输入（Q<0、L<=0、Nx=1、CONTACT_BC/KAPPA_SOURCE 拼错、Q=0）；
  7. S5 界面检测（层间距 0.7~3.5 A × x/y/z × 元素序）与无 pbc；
  8. S6 用合成数据回归，并确认稳态窗口丢掉初始暂态。

运行: python3 tests/suite_device_thermal.py   （rc=0 全过）
"""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
THERMAL = os.path.join(ROOT, "skill", "_common", "thermal")
OPT = os.path.join(ROOT, "skill", "_common", "opt")
DT = os.path.join(ROOT, "skill", "device-thermal")
FEM = os.path.join(ROOT, "skill", "device-thermal-fem")
for _p in (THERMAL, OPT, DT):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import device_common as DC  # noqa: E402

FAILS = []
NCHECK = [0]


def ok(cond, msg):
    NCHECK[0] += 1
    if cond:
        print("  . %s" % msg)
    else:
        print("  x %s" % msg)
        FAILS.append(msg)


def approx(got, want, rtol, msg, atol=0.0):
    NCHECK[0] += 1
    err = abs(got - want)
    good = err <= atol + rtol * abs(want)
    if good:
        print("  . %s (got=%.8g want=%.8g)" % (msg, got, want))
    else:
        print("  x %s (got=%.8g want=%.8g err=%.3g)" % (msg, got, want, err))
        FAILS.append(msg)


def _db_layer(name, th=None):
    db = DC.load_materials()
    m = db["materials"][name]
    return {"material": name, "kind": m["kind"],
            "kx_W_mK": float(m["kappa_inplane_W_mK"]),
            "ky_W_mK": float(m["kappa_cross_W_mK"]),
            "thickness_m": float(th) if th is not None else m["thickness_m"],
            "rho_cp_J_m3K": float(m["rho_cp_J_m3K"]), "source": "database"}


def make_props(ch="MoS2", ox="SiO2", ox_th=300e-9, tbc=1.4e7):
    p = {"channel": _db_layer(ch), "oxide": _db_layer(ox, ox_th), "tbc_W_m2K": tbc}
    if p["channel"]["thickness_m"] is None:
        p["channel"]["thickness_m"] = 6.7177e-10
    return p


def make_spec(**kw):
    s = {"DEVICE_SPEC": True, "device_name": "test", "channel_material": "MoS2",
         "channel_thickness_m": None, "tbc_override_W_m2K": None,
         "oxide_material": "SiO2", "oxide_thickness_m": 300e-9,
         "substrate_material": "Si", "contact_bc": "sink",
         "L_channel_m": 200e-9, "W_device_m": 1e-6, "Q_joule_W_m2": 1e8,
         "T_amb_K": 300.0, "kappa_source": "literature",
         "upstream_skill": "kl-dft-cpu", "upstream_step": "step6_kappa",
         "n_grid_x": 41, "n_grid_sub": 60, "do_transient": False,
         "sor_omega": 1.7, "sor_tol": 1e-10, "sor_maxit": 200000}
    s.update(kw)
    return s


def solve(spec, props=None, **kw):
    return DC.solve_device(spec, props or make_props(), do_transient=False, **kw)[0]


# ---------------------------------------------------------------- 求解器
def test_validation():
    print("[validation]")
    p = make_props()
    for kw, label in (({"Q_joule_W_m2": -1.0}, "Q_JOULE<0"),
                      ({"L_channel_m": 0.0}, "L_CHANNEL<=0")):
        try:
            DC.solve_device(make_spec(**kw), p, do_transient=False)
            ok(False, "validate %s 应报错" % label)
        except SystemExit:
            ok(True, "validate %s 报错" % label)
    for val in ("sinks", "sink_", "", "conducting"):
        try:
            DC.solve_device(make_spec(contact_bc=val), p, do_transient=False)
            ok(False, "CONTACT_BC=%r 拼错应报错" % val)
        except SystemExit:
            ok(True, "CONTACT_BC=%r 拼错报错" % val)
    try:
        DC.solve_device(make_spec(), p, n_grid_x=1, do_transient=False)
        ok(False, "N_GRID_X=1 应报错")
    except SystemExit:
        ok(True, "N_GRID_X=1 报错")
    s0 = solve(make_spec(Q_joule_W_m2=0.0))
    ok(s0["dT_peak_K"] == 0.0 and s0["dT_ratio_FVM_over_1D"] is None,
       "Q_JOULE=0 行为明确：dT=0、ratio=n/a（不抛 TypeError）")
    tmp = tempfile.mkdtemp(prefix="dt_val_")
    try:
        d = os.path.join(tmp, "mat", "device-thermal")
        os.makedirs(d)
        try:
            DC.resolve_thermal_props(d, "MoS2", "SiO2", 300.0, kappa_source="bogus")
            ok(False, "KAPPA_SOURCE=bogus 应报错")
        except SystemExit:
            ok(True, "KAPPA_SOURCE=bogus 报错")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_energy_conservation():
    print("[energy]")
    p = make_props()
    for bc in ("sink", "insulator"):
        s = solve(make_spec(contact_bc=bc), p)
        approx(s["P_injected_W_per_m"], 1e8 * 200e-9, 1e-12,
               "注入功率 = Q*L (%s)" % bc)
        approx(s["energy_balance_rel"], 0.0, 0.0, "能量守恒 rel (%s)" % bc, atol=1e-7)
        if bc == "insulator":
            approx(s["P_side_W_per_m"], 0.0, 0.0, "绝热侧壁 P_side=0", atol=1e-12)


def test_analytic_1d_insulator():
    print("[analytic-1d]")
    s = solve(make_spec(contact_bc="insulator"))
    rel = abs(s["dT_peak_K"] / s["dT_peak_1D_distributed_K"] - 1.0)
    ok(rel < 1e-3, "绝热 FVM 对齐精确 1D（含 t/(2ky)）rel=%.2e" % rel)


def test_x_convergence():
    print("[x-convergence]")
    dts = [solve(make_spec(), n_grid_x=n)["dT_peak_K"] for n in (21, 41, 81, 161)]
    ok(all(dts[i + 1] < dts[i] for i in range(3)),
       "x 网格单调收敛: %s" % ["%.4f" % v for v in dts])
    rel = abs(dts[1] - dts[3]) / dts[3]
    ok(rel < 0.01, "默认 Nx=41 vs Nx=161 差 %.3f%% < 1%%" % (100 * rel))
    rich = dts[2] - (dts[1] - dts[2]) / 3.0
    ok(abs(rich - dts[3]) / dts[3] < 1e-3,
       "Richardson 外推 %.5f 对 161 网格 rel=%.1e" % (rich, abs(rich - dts[3]) / dts[3]))


def test_y_convergence():
    print("[y-convergence]")
    dts = [solve(make_spec(), n_grid_sub=n)["dT_peak_K"] for n in (15, 30, 60, 120)]
    diffs = [abs(dts[i + 1] - dts[i]) for i in range(3)]
    ok(diffs[-1] < diffs[0], "y 网格加密变化量减小: %s" % ["%.4f" % v for v in dts])
    ok(abs(dts[3] - dts[2]) / dts[3] < 5e-3, "nsub=60 vs 120 差 < 0.5%")


def test_anisotropic_oxide():
    print("[hBN-oxide]")
    p = make_props(ox="hBN", ox_th=20e-9)
    s = solve(make_spec(oxide_material="hBN", oxide_thickness_m=20e-9), p)
    ch, ox = p["channel"], p["oxide"]
    want = (ch["thickness_m"] / ch["ky_W_mK"] + 1.0 / p["tbc_W_m2K"]
            + ox["thickness_m"] / ox["ky_W_mK"])
    approx(s["R_pp_1D_network_m2K_per_W"], want, 1e-12, "hBN 氧化层 R_pp 用跨面 ky")
    thick = make_props(ox="hBN", ox_th=2e-6)
    s2 = solve(make_spec(oxide_material="hBN", oxide_thickness_m=2e-6), thick)
    _ch2, _ox2 = thick["channel"], thick["oxide"]
    wrong = (_ch2["thickness_m"] / _ch2["ky_W_mK"] + 1.0 / thick["tbc_W_m2K"]
             + _ox2["thickness_m"] / _ox2["kx_W_mK"])
    ok(abs(s2["R_pp_1D_network_m2K_per_W"] / wrong - 1.0) > 10.0,
       "R_pp 确实不是用面内 kx（厚 hBN 差异 >10x）")
    tmp = tempfile.mkdtemp(prefix="dt_hbn_")
    try:
        d = os.path.join(tmp, "mat", "device-thermal")
        os.makedirs(d)
        _c, ox2, _t, _r, _u, _s = DC.resolve_thermal_props(
            d, "MoS2", "hBN", 300.0, kappa_source="literature", oxide_thickness=20e-9)
        approx(ox2["thickness_m"], 20e-9, 1e-12, "OXIDE_THICKNESS 覆盖 hBN 库厚度 0.34 nm")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_fem_crosscheck():
    print("[FEM-crosscheck]")
    try:
        import skfem  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        print("  . skip（skfem 不可用: %s）" % exc)
        return
    if FEM not in sys.path:
        sys.path.insert(0, FEM)
    import fem_common as FC
    import backend_skfem as BK
    for bc in ("insulator", "sink"):
        spec = make_spec(contact_bc=bc)
        p = make_props()
        fem = BK.solve(FC.build_model(spec, p), opts={"sink": (bc == "sink")})
        fvm = solve(spec, p)
        rel = abs(fvm["dT_peak_K"] / fem["dT_peak_K"] - 1.0)
        ok(rel < 0.01, "FVM vs 独立 FEM（%s）rel=%.3f%%" % (bc, 100 * rel))
    fem = BK.solve(FC.build_model(make_spec(), make_props()), opts={"sink": True})
    approx(fem["P_source_W_per_m"], 1e8 * 200e-9, 2e-2, "FEM 注入功率 = Q*L")


# ---------------------------------------------------------------- 参数流向
def _write_conf(d, params):
    with open(os.path.join(d, "step.conf"), "w", encoding="utf-8") as fh:
        fh.write("[params]\n")
        for k, v in params.items():
            fh.write("%s = %s\n" % (k, v))


def _copy_dt_assets(dest):
    for f in ("gen_step1_device_spec.py", "gen_step2_props.py", "gen_step3_solve.py"):
        shutil.copyfile(os.path.join(DT, f), os.path.join(dest, f))
    shutil.copyfile(os.path.join(OPT, "stepconf.py"), os.path.join(dest, "stepconf.py"))
    shutil.copyfile(os.path.join(THERMAL, "device_common.py"),
                    os.path.join(dest, "device_common.py"))
    shutil.copyfile(os.path.join(THERMAL, "device_materials.json"),
                    os.path.join(dest, "device_materials.json"))


def _run_gen(d, script):
    return subprocess.run([sys.executable, script], cwd=d, capture_output=True, text=True)


def test_param_flow():
    print("[param-flow]")
    tmp = tempfile.mkdtemp(prefix="dt_flow_")
    try:
        _copy_dt_assets(tmp)
        _write_conf(tmp, {"N_GRID_X": "21", "N_GRID_SUB": "12", "DO_TRANSIENT": "false"})
        r = _run_gen(tmp, "gen_step1_device_spec.py")
        ok(r.returncode == 0, "S1 生成 rc=0")
        spec = json.load(open(os.path.join(tmp, "step1_device_spec/device_spec.json")))
        ok(spec["n_grid_x"] == 21 and spec["n_grid_sub"] == 12, "S1 把 N_GRID_X=21 固化进 spec")

        _write_conf(tmp, {"TBC_OVERRIDE": "1.2e8", "CHANNEL_THICKNESS": "1e-9",
                          "KAPPA_SOURCE": "literature", "DO_TRANSIENT": "false"})
        r = _run_gen(tmp, "gen_step2_props.py")
        ok(r.returncode == 0, "S2 生成 rc=0")
        props = json.load(open(os.path.join(tmp, "step2_props/thermal_props.json")))
        approx(props["tbc_W_m2K"], 1.2e8, 1e-12, "S2 上设 TBC_OVERRIDE 真的生效")
        approx(props["channel"]["thickness_m"], 1e-9, 1e-12, "S2 上设 CHANNEL_THICKNESS 真的生效")
        approx(props["channel"]["kx_W_mK"] * props["channel"]["thickness_m"],
               63.7 * 6.7177e-10, 1e-9, "改厚度时 κ·t 保持（面内热导不变）")

        _write_conf(tmp, {"DO_TRANSIENT": "false"})
        r = _run_gen(tmp, "gen_step3_solve.py")
        ok(r.returncode == 0, "S3 生成 rc=0（网格回退 S1）")
        summ = json.load(open(os.path.join(tmp, "step3_solve/device_thermal_summary.json")))
        ok(summ["Nx"] == 21 and summ["Ny"] == 14, "S3 未设网格时回退 S1 的 Nx=21/Ny=14")

        _write_conf(tmp, {"N_GRID_X": "31", "N_GRID_SUB": "7", "DO_TRANSIENT": "false"})
        r = _run_gen(tmp, "gen_step3_solve.py")
        summ = json.load(open(os.path.join(tmp, "step3_solve/device_thermal_summary.json")))
        ok(summ["Nx"] == 31 and summ["Ny"] == 9, "S3 上显式 N_GRID_X=31 覆盖 S1")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- 物性来源
def _write_kappa(matdir, obj):
    up = os.path.join(matdir, "kl-dft-cpu", "step6_kappa")
    os.makedirs(up, exist_ok=True)
    with open(os.path.join(up, "kappa_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(obj, fh)
    with open(os.path.join(up, "thickness_2d.json"), "w", encoding="utf-8") as fh:
        json.dump({"thickness_d_A": 6.7177}, fh)


def test_kappa_source_modes():
    print("[kappa-source]")
    tmp = tempfile.mkdtemp(prefix="dt_ksrc_")
    try:
        mat = os.path.join(tmp, "mat")
        d = os.path.join(mat, "device-thermal")
        os.makedirs(d)
        _c, _o, _t, _r, _u, src = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, kappa_source="literature", oxide_thickness=300e-9)
        ok(src == "literature" and _c["source"] == "database", "literature：用材料库")
        try:
            DC.resolve_thermal_props(d, "MoS2", "SiO2", 300.0,
                                     kappa_source="upstream", oxide_thickness=300e-9)
            ok(False, "upstream 缺上游应报错（不能静默回退）")
        except SystemExit:
            ok(True, "upstream 缺上游直接报错")
        c_a, _o, _t, _r, _u, src = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, kappa_source="auto", oxide_thickness=300e-9)
        ok(c_a["source"] == "database", "auto 缺上游：回退材料库并继续")
        # 形态 1：normalized 张量
        _write_kappa(mat, {"dim": "2d", "temperatures": [300.0],
                           "kappa_xx_yy_zz": [[20.0, 22.0, 0.0]],
                           "kappa_2d_norm_factor": 2.973164,
                           "kappa_2d_normalized_xx_yy_zz": [[59.46328, 65.409608, 0.0]]})
        c1, _o, _t, _r, _u, _s = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, kappa_source="upstream", oxide_thickness=300e-9)
        approx(c1["kx_W_mK"], 0.5 * (59.46328 + 65.409608), 1e-12, "上游 normalized 张量")
        # 形态 2：只有归一化因子
        _write_kappa(mat, {"dim": "2d", "temperatures": [300.0],
                           "kappa_xx_yy_zz": [[20.0, 22.0, 0.0]],
                           "kappa_2d_norm_factor": 2.973164})
        c2, _o, _t, _r, _u, _s = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, kappa_source="upstream", oxide_thickness=300e-9)
        approx(c2["kx_W_mK"], 0.5 * (20.0 + 22.0) * 2.973164, 1e-12, "只有因子：现场归一化")
        # 形态 3：都没有（raw）—— 仍可用，但带警告
        _write_kappa(mat, {"dim": "2d", "temperatures": [300.0],
                           "kappa_xx_yy_zz": [[20.0, 22.0, 0.0]]})
        c3, _o, _t, _r, _u, _s = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, kappa_source="auto", oxide_thickness=300e-9)
        approx(c3["kx_W_mK"], 21.0, 1e-12, "都没有：原样用并给警告（auto）")
        # 形态 4：只有 runs[]
        _write_kappa(mat, {"runs": [{"method": "rta", "temperatures": [300.0],
                                     "kappa_xx_yy_zz": [[20.0, 22.0, 0.0]],
                                     "kappa_2d_normalized_xx_yy_zz": [[59.46328, 65.409608, 0.0]]}]})
        c4, _o, _t, _r, _u, _s = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, kappa_source="upstream", oxide_thickness=300e-9)
        approx(c4["kx_W_mK"], 0.5 * (59.46328 + 65.409608), 1e-12, "只有 runs[]：取 RTA 运行")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_thickness_scaling():
    print("[thickness-scaling]")
    tmp = tempfile.mkdtemp(prefix="dt_thick_")
    try:
        mat = os.path.join(tmp, "mat")
        d = os.path.join(mat, "device-thermal")
        os.makedirs(d)
        _write_kappa(mat, {"dim": "2d", "temperatures": [300.0],
                           "kappa_xx_yy_zz": [[20.0, 22.0, 0.0]],
                           "kappa_2d_norm_factor": 2.973164,
                           "kappa_2d_normalized_xx_yy_zz": [[59.46328, 65.409608, 0.0]]})
        kx_up = 0.5 * (59.46328 + 65.409608)
        t_up = 6.7177e-10
        c0, _o, _t, _r, _u, _s = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, up_skill="kl-dft-cpu", up_step="step6_kappa",
            kappa_source="upstream", oxide_thickness=300e-9)
        approx(c0["thickness_m"], t_up, 1e-12, "上游厚度采用")
        c1, _o, _t, _r, _u, _s = DC.resolve_thermal_props(
            d, "MoS2", "SiO2", 300.0, up_skill="kl-dft-cpu", up_step="step6_kappa",
            kappa_source="upstream", thickness_override=1e-9, oxide_thickness=300e-9)
        approx(c1["kx_W_mK"] * c1["thickness_m"], kx_up * t_up, 1e-12,
               "上游厚度 + CHANNEL_THICKNESS：κ·t 保持")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- S5 / S6
def _two_layer(spacing, axis=2, nA=10, nB=10):
    coords, species = [], []
    for i in range(nA):
        p = [0.0, 0.0, 0.0]
        p[axis] = i * spacing
        coords.append(p)
        species.append("Si")
    zbase = (nA - 1) * spacing + spacing
    for i in range(nB):
        p = [0.0, 0.0, 0.0]
        p[axis] = zbase + i * spacing
        coords.append(p)
        species.append("Ge")
    return [p[axis] for p in coords], species, zbase - spacing / 2.0


def _load_step5():
    spec = importlib.util.spec_from_file_location("dt_g5", os.path.join(DT, "gen_step5_tbc.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_step6():
    spec = importlib.util.spec_from_file_location("dt_g6", os.path.join(DT, "gen_step6_tbc_post.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_s5_align_coords():
    print("[s5-align-coords]")
    g5 = _load_step5()
    lat = [10.0, 0, 0, 0, 10.0, 0, 0, 0, 70.0]
    # 真实 bug：非周期 z 轴上一层原子在 -0.2 A，旧实现对三轴都取模，把它挪到 +69.8 并分进源组
    pos = [(0.3, 0.0, -0.2), (0.0, 0.0, 1.0), (0.0, 0.0, 69.8)]
    new, isT, _warn = g5.align_coords(pos, lat, "T T F")
    ok(isT == [True, True, False], "pbc='T T F' 解析为 [T,T,F]")
    zs = [round(p[2], 6) for p in new]
    ok(min(zs) == 0.0, "非周期 z 轴平移后最小值为 0（旧版取模把 -0.2 挪到 +69.8）")
    ok(zs == [0.0, 1.2, 70.0], "非周期 z 轴整体平移、保持相对构型 z=%s" % zs)
    ok(all(abs(n[0] - p[0]) < 1e-12 and abs(n[1] - p[1]) < 1e-12
           for n, p in zip(new, pos)), "非周期轴平移不改变周期轴坐标")
    new2, _i2, _w2 = g5.align_coords([(-0.5, 0.0, 1.0)], lat, "T T F")
    approx(new2[0][0], 9.5, 0, "周期 x 轴仍 wrap: -0.5 -> 9.5", atol=1e-12)
    new3, isT3, _w3 = g5.align_coords([(-0.5, 11.0, -1.0)], lat, "T T T")
    ok(all(isT3), "pbc='T T T' 解析为全周期")
    approx(new3[0][0], 9.5, 0, "x wrap", atol=1e-12)
    approx(new3[0][1], 1.0, 0, "y wrap", atol=1e-12)
    approx(new3[0][2], 69.0, 0, "z wrap", atol=1e-12)


def test_interface_detection():
    print("[interface-detection]")
    g5 = _load_step5()
    allok = True
    for spacing in (0.7, 1.0, 1.36, 2.0, 2.77, 3.2, 3.5):
        for axis in (0, 1, 2):
            coords, species, want = _two_layer(spacing, axis)
            got = g5.auto_interface_by_composition(coords, species)
            good = got is not None and abs(got - want) < 1e-9
            allok = allok and good
    ok(allok, "层间距 0.7~3.5 A × x/y/z 全部命中界面（旧版 2.0+ 失效）")
    coords, species, want = _two_layer(3.2, 2)
    species = ["Ge" if s == "Si" else "Si" for s in species]
    got = g5.auto_interface_by_composition(coords, species)
    ok(got is not None and abs(got - want) < 1e-9, "元素序颠倒（Ge 低 / Si 高）仍命中")
    # 无 pbc 的结构
    tmp = tempfile.mkdtemp(prefix="dt_pbc_")
    try:
        p = os.path.join(tmp, "s.xyz")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write('2\nLattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3\n'
                     "Si 0 0 0\nGe 0 0 1\n")
        _lat, pbc, _sym, _pos, _props = g5.parse_extxyz(p)
        ok(pbc is None, "无 pbc 的结构解析为 None（S5 会要求显式 TBC_PBC）")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_s6_regression():
    print("[S6-regression]")
    EV = 1.602176634e-4
    A2 = 1e-20
    dt_fs, area, slope, dti = 1000.0, 100.0, 0.5, 10.0
    q_exp = (slope / dt_fs) * EV / (area * A2)
    tmp = tempfile.mkdtemp(prefix="dt_s6_")
    try:
        shutil.copyfile(os.path.join(DT, "gen_step6_tbc_post.py"),
                        os.path.join(tmp, "gen_step6_tbc_post.py"))
        shutil.copyfile(os.path.join(OPT, "stepconf.py"), os.path.join(tmp, "stepconf.py"))
        with open(os.path.join(tmp, "step.conf"), "w", encoding="utf-8") as fh:
            fh.write("[params]\n")
        s5 = os.path.join(tmp, "step5_tbc")
        os.makedirs(s5)
        meta = {"axis": "z", "interface_coord_A": 50.0, "coord_min_A": 0.0,
                "coord_max_A": 100.0, "source_thickness_A": 10.0,
                "sink_thickness_A": 10.0, "source_side": "high", "area_A2": area,
                "time_step_fs": 1.0, "sample_interval": 10, "output_interval": 100}
        json.dump(meta, open(os.path.join(s5, "tbc_inputs.json"), "w"))
        # 前 40 行是不自洽的暂态，之后 80 行稳态
        with open(os.path.join(s5, "compute.out"), "w") as fh:
            for i in range(120):
                if i < 40:
                    fh.write("%.10f %.10f\n" % (2.0 * i, -0.2 * i))
                else:
                    fh.write("%.10f %.10f\n" % (2.0 * 40 + slope * (i - 40),
                                                -0.2 * 40 - slope * (i - 40)))
        # 前 40 个 block 用错误的剖面（若被平均进去 dT_i 会明显偏）
        with open(os.path.join(s5, "compute_chunk.out"), "w") as fh:
            for b in range(120):
                for cid in range(50):
                    z = 2.0 * cid + 1.0
                    if 12.0 < z < 48.0:
                        t = (9.0 * z + 4.0 * 273.75) if b < 40 else (0.625 * z + 273.75)
                    elif 52.0 < z < 88.0:
                        t = (9.0 * z + 4.0 * 308.75) if b < 40 else (0.125 * z + 308.75)
                    else:
                        t = 300.0
                    fh.write("%d %.6f 10.0 %.6f\n" % (cid, z, t))
        r = subprocess.run([sys.executable, "gen_step6_tbc_post.py"],
                           cwd=tmp, capture_output=True, text=True)
        ok(r.returncode == 0, "S6 合成数据运行 rc=0")
        if r.returncode != 0:
            print(r.stdout[-500:])
            print(r.stderr[-500:])
            return
        out = json.load(open(os.path.join(tmp, "step6_tbc_post/tbc_result.json")))
        approx(out["q_W_m2"], q_exp, 1e-6, "S6 提取 q")
        approx(out["dT_interface_K"], dti, 1e-6, "S6 dT_i=10 K（稳态窗口）")
        approx(out["G_W_m2K"], q_exp / dti, 1e-6, "S6 G = q/dT_i")
        ok(out["quality"] == "good", "S6 quality=good")
        approx(out["steady_window"]["rows_used"], 80, 1e-12, "稳态窗口丢弃前 40 行暂态")
        approx(out["n_blocks_averaged"], 80, 1e-12, "剖面只平均稳态 block")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_s6_temporal_stats():
    print("[s6-temporal-stats]")
    g6 = _load_step6()
    # 只有 2 个时间块：误差必须返回 None（调用方据此判 poor，而不是跳过 TBC_DT_TOL）
    blocks2 = [{0: (0.0, 300.0, 10.0), 1: (1.0, 301.0, 10.0)},
               {0: (0.0, 300.0, 10.0), 1: (1.0, 301.0, 10.0)}]
    se, _h, _c, n = g6.temporal_stats(blocks2, [0, 1], [0, 1], 0.5)
    ok(se is None and n == 2, "时间块 < min_blocks 时误差返回 None（旧版空间分块会整段失效）")
    # 6 个时间块，每块斜率有噪声：时间分块应给出 >0 的 dT_i 标准误
    blocks6 = []
    for b in range(6):
        d = float(0.3 * (b - 2.5))
        blocks6.append({0: (0.0, 300.0 + d, 10.0), 1: (1.0, 346.0 + d, 10.0),
                        2: (2.0, 388.0, 10.0), 3: (3.0, 434.0, 10.0)})
    se6, sdH, sdC, n6 = g6.temporal_stats(blocks6, [0, 1], [2, 3], 1.5)
    ok(n6 == 6 and se6 is not None and se6 > 0.0, "时间分块给出 dT_i 标准误 se=%.4g" % (se6 or -1))
    ok(sdH is not None and sdC is not None, "同时给出两侧斜率的时间离散度")


def test_s6_noise_is_poor():
    print("[s6-noise-poor]")
    EV = 1.602176634e-4
    A2 = 1e-20
    dt_fs, area, slope = 1000.0, 100.0, 0.5
    rng = np.random.default_rng(20260925)
    tmp = tempfile.mkdtemp(prefix="dt_s6n_")
    try:
        shutil.copyfile(os.path.join(DT, "gen_step6_tbc_post.py"),
                        os.path.join(tmp, "gen_step6_tbc_post.py"))
        shutil.copyfile(os.path.join(OPT, "stepconf.py"), os.path.join(tmp, "stepconf.py"))
        with open(os.path.join(tmp, "step.conf"), "w", encoding="utf-8") as fh:
            fh.write("[params]\n")
        s5 = os.path.join(tmp, "step5_tbc")
        os.makedirs(s5)
        # 源/漏各 30 A 厚 -> 拟合窗约 16 A，bin 宽 4 A -> 每侧只有 4 个 bin（用户场景）
        meta = {"axis": "z", "interface_coord_A": 50.0, "coord_min_A": 0.0,
                "coord_max_A": 100.0, "source_thickness_A": 30.0,
                "sink_thickness_A": 30.0, "source_side": "high", "area_A2": area,
                "time_step_fs": 1.0, "sample_interval": 10, "output_interval": 100}
        json.dump(meta, open(os.path.join(s5, "tbc_inputs.json"), "w"))
        with open(os.path.join(s5, "compute.out"), "w") as fh:
            for i in range(120):
                if i < 40:
                    fh.write("%.10f %.10f\n" % (2.0 * i, -0.2 * i))
                else:
                    fh.write("%.10f %.10f\n" % (2.0 * 40 + slope * (i - 40),
                                                -0.2 * 40 - slope * (i - 40)))
        with open(os.path.join(s5, "compute_chunk.out"), "w") as fh:
            for b in range(120):
                for cid in range(50):
                    z = 4.0 * cid + 1.0
                    if 12.0 < z < 48.0:
                        t = 0.625 * z + 273.75
                    elif 52.0 < z < 88.0:
                        t = 0.125 * z + 308.75
                    else:
                        t = 300.0
                    if b >= 40 and (12.0 < z < 48.0 or 52.0 < z < 88.0):
                        t += float(rng.normal(0.0, 8.0))
                    fh.write("%d %.6f 10.0 %.6f\n" % (cid, z, t))
        r = subprocess.run([sys.executable, "gen_step6_tbc_post.py"],
                           cwd=tmp, capture_output=True, text=True)
        ok(r.returncode == 0, "S6 少 bin + 大噪声数据运行 rc=0")
        if r.returncode != 0:
            print(r.stdout[-600:])
            print(r.stderr[-600:])
            return
        out = json.load(open(os.path.join(tmp, "step6_tbc_post/tbc_result.json")))
        unc = out.get("dT_interface_uncertainty_K")
        grel = out.get("G_relative_uncertainty")
        ok(out["quality"] == "poor",
           "少 bin + 大噪声判 poor（旧版空间分块误差为 None 会跳过判据误判 good）")
        ok(unc is None or grel is None or grel > 0.3,
           "误差判据生效 dT_unc=%r G_rel=%r" % (unc, grel))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    print("=" * 62)
    print("suite_device_thermal  (%s)" % ROOT)
    print("=" * 62)
    for fn in (test_validation, test_energy_conservation, test_analytic_1d_insulator,
               test_x_convergence, test_y_convergence, test_anisotropic_oxide,
               test_fem_crosscheck, test_param_flow, test_kappa_source_modes,
               test_thickness_scaling, test_s5_align_coords, test_interface_detection,
               test_s6_regression, test_s6_temporal_stats, test_s6_noise_is_poor):
        try:
            fn()
        except Exception as exc:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            FAILS.append("%s 抛异常: %s" % (fn.__name__, exc))
    print("-" * 62)
    if FAILS:
        print("suite_device_thermal: FAILED %d / %d 项" % (len(FAILS), NCHECK[0]))
        for f in FAILS[:20]:
            print("   -", f)
        return 1
    print("suite_device_thermal: ALL PASS (%d 项)" % NCHECK[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
