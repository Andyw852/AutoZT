"""tests/test_fcfit_imag_findiff.py —— 三处修复回归（纯本地，不碰集群）。

1. gen_step1_fit 把虚频阈值按【大写键】写进 fit_config.json：
   imag_policy._get 只查大写 IMAG_THR/IMAG_THR_STRICT/IMAG_QGAMMA/IMAG_QGAMMA_GRACE，
   之前写小写 imag_thr -> 静默回落默认值 0.10；且 kl S5 声明的 ZA_CHECK/ZA_QMAX 没写进去。
2. fc_fit_driver._dataset_from_yaml 识别 phono3py 有限位移（顶层 displacement_pairs）。
3. kl-dft-cpu reuse_fcfit 产出【顶层布局】（fc2/fc3 + cutoff_scan/cut3_*），
   共享 S6 gen_step2_kappa._gen_one 只认顶层文件。
"""
import contextlib
import importlib.util
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FCFIT = ROOT / "skill" / "_common" / "fcfit"
OPT = ROOT / "skill" / "_common" / "opt"
COMMON = ROOT / "skill" / "_common"

for _d in (str(COMMON), str(FCFIT), str(OPT)):
    if _d not in sys.path:
        sys.path.insert(0, _d)


def _load(path, name):
    for d in (str(Path(path).parent), str(OPT), str(COMMON)):
        if d not in sys.path:
            sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


def test_spec_declares_imag_and_za_keys():
    g1 = _load(FCFIT / "gen_step1_fit.py", "g1_imag_spec")
    for k in ("IMAG_THR", "IMAG_THR_STRICT", "IMAG_QGAMMA", "IMAG_QGAMMA_GRACE",
              "ZA_CHECK", "ZA_QMAX"):
        assert k in g1.SPEC, k


def test_fit_config_writes_uppercase_imag_keys(tmp_path, monkeypatch):
    import numpy as np
    g1 = _load(FCFIT / "gen_step1_fit.py", "g1_imag_e2e")
    import fc_common
    ds = tmp_path / "mat" / "step4_disp"
    sk = tmp_path / "mat" / "fit-fc-thermal"
    ds.mkdir(parents=True)
    sk.mkdir()
    ds.joinpath("POSCAR").write_text(
        "sc\n1.0\n3 0 0\n0 3 0\n0 0 3\nX\n1\nDirect\n0 0 0\n")
    pos = "".join("%g %g %g\n" % (i / 2, j / 2, k / 2)
                  for i in range(2) for j in range(2) for k in range(2))
    ds.joinpath("SPOSCAR").write_text(
        "sc\n1.0\n6 0 0\n0 6 0\n0 0 6\nX\n8\nDirect\n" + pos)
    dd = np.random.default_rng(1).normal(0, 0.03, (6, 8, 3))
    dd[-1] = 0
    np.save(ds / "dataset_disps.npy", dd)
    np.save(ds / "dataset_forces.npy", -dd)
    sk.joinpath("step.conf").write_text(
        "[params]\nSTEP = step1_fit\nFIT_ENGINE = pheasy\nCUT3_CANDIDATES = off\n"
        "IMAG_THR = 0.5\nIMAG_THR_STRICT = 0.2\nIMAG_QGAMMA = 0.1\n"
        "IMAG_QGAMMA_GRACE = 1.5\nZA_CHECK = on\nZA_QMAX = 0.08\n"
        "MIN_SUCCESS_RATIO = 0.7\nMIN_SUCCESS_FRAMES = 2\n")
    for t in ("submit_fcfit_p3py", "submit_fcfit_pheasy", "submit_fcfit_hiphive",
              "submit_fcfit_pheasy_gpu"):
        (sk / (t + ".tpl")).write_text(
            (FCFIT / "templates" / (t + ".tpl")).read_text(encoding="utf-8"))
    monkeypatch.chdir(sk)
    monkeypatch.setattr(fc_common, "resolve_submit",
                        lambda base, kind: sk / (kind + ".tpl"))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        g1.main()
    cfg = json.loads((sk / "step1_fit" / "fit_config.json").read_text())
    # 大写键写进 fit_config.json，旧的小写 imag_thr 不再出现
    assert cfg["IMAG_THR"] == 0.5
    assert "imag_thr" not in cfg
    assert cfg["IMAG_THR_STRICT"] == 0.2
    assert cfg["IMAG_QGAMMA"] == 0.1
    assert cfg["IMAG_QGAMMA_GRACE"] == 1.5
    assert cfg["ZA_CHECK"] == "on"
    assert cfg["ZA_QMAX"] == 0.08
    # 缺帧门禁（第 5 条）也写进 fit_config.json
    assert cfg["min_success_ratio"] == 0.7
    assert cfg["min_success_frames"] == 2
    # 下游 imag_policy 确实能从这份 cfg 读到（大写键生效，不再回落 0.10）
    import imag_policy
    thr = imag_policy.thresholds(cfg)
    assert thr["IMAG_THR"] == 0.5
    assert thr["IMAG_THR_STRICT"] == 0.2
    assert thr["IMAG_QGAMMA"] == 0.1
    assert thr["IMAG_QGAMMA_GRACE"] == 1.5


def test_dataset_from_yaml_recognizes_findiff():
    import fc_fit_driver as drv
    # phono3py 有限位移（type-1）：位移在顶层 displacement_pairs，没有 dataset 块
    findiff = {"displacement_pairs": [{"atom": 1, "displacement_id": 1}],
               "displacement_pair_info": {}}
    kind, disps, forces = drv._dataset_from_yaml(findiff)
    assert kind == "findiff" and disps is None and forces is None
    # 随机位移（type-2）：dataset.displacements / forces
    random_y = {"dataset": {"displacements": [[[0.01, 0, 0], [0, 0.01, 0]]],
                            "forces": [[[1.0, 0, 0], [0, 1.0, 0]]]}}
    kind, disps, forces = drv._dataset_from_yaml(random_y)
    assert kind == "random" and disps is not None and forces is not None
    # 空 / 未知结构保持 None（不抛）
    assert drv._dataset_from_yaml(None) == (None, None, None)
    assert drv._dataset_from_yaml({}) == (None, None, None)


def test_reuse_fcfit_top_level_layout(tmp_path):
    import fc_fit_driver as drv
    src = tmp_path / "step1_fit"
    src.mkdir()
    for f in ("fc2.hdf5", "fc3.hdf5", "phonon_summary.json", "fc_fit_summary.json",
              "cutoff_scan.json", "fc_dataset.json", "phono3py_disp.yaml"):
        (src / f).write_bytes(b"x")
    (src / "cutoff_scan").mkdir()
    (src / "cutoff_scan" / "cut3_5p00").mkdir()
    (src / "cutoff_scan" / "cut3_5p00" / "fc2.hdf5").write_bytes(b"z")
    (src / "cutoff_scan" / "cut3_5p00" / "fc3.hdf5").write_bytes(b"z")
    (src / "shengbte").mkdir()
    (src / "shengbte" / "FORCE_CONSTANTS_2ND").write_bytes(b"y")
    drv._write_kl_bundle(src, {})
    b = src / "kl_bundle"
    rf = _load(ROOT / "skill" / "kl-dft-cpu" / "reuse_fcfit.py", "rf_imag")
    out = tmp_path / "step5_fc"
    rf.assemble(b, out, supercell="3 3 1", dim="3d", mesh="24 24 24")
    # S6 _gen_one 只认顶层 fc2/fc3 与 cutoff_scan/cut3_*，旧的 phono3py/ 子目录必须消失
    assert (out / "fc2.hdf5").is_file() and (out / "fc3.hdf5").is_file()
    assert not (out / "phono3py").exists()
    assert (out / "cutoff_scan" / "cut3_5p00" / "fc2.hdf5").is_file()
    assert (out / "cutoff_scan" / "cut3_5p00" / "fc3.hdf5").is_file()
    assert (out / "fc_fit_summary.json").is_file()
    assert (out / "fc_dataset.json").is_file()
    assert (out / "shengbte" / "FORCE_CONSTANTS_2ND").is_file()


def test_structures_match_and_provenance(tmp_path):
    import numpy as np
    g5 = _load(ROOT / "skill" / "kl-dft-cpu" / "gen_step5_fc.py", "g5_prov")
    # 结构签名比对
    sig_a = (np.eye(3) * 3.0, ["Si"], np.array([[0.0, 0.0, 0.0]]))
    sig_b = (np.eye(3) * 3.0, ["Si"], np.array([[0.0, 0.0, 0.0]]))
    assert g5._structures_match(sig_a, sig_b)[0] is True
    sig_c = (np.eye(3) * 3.1, ["Si"], np.array([[0.0, 0.0, 0.0]]))   # 晶格不同
    assert g5._structures_match(sig_a, sig_c)[0] is False
    sig_d = (np.eye(3) * 3.0, ["Ge"], np.array([[0.0, 0.0, 0.0]]))   # 元素不同
    assert g5._structures_match(sig_a, sig_d)[0] is False
    assert g5._structures_match(None, None)[0] is None
    # 来源校验：dataset_dir 指向 kl-mlff → 默认拒绝，allow_foreign 放行
    src = tmp_path / "bundle"
    src.mkdir()
    (src / "fc_dataset.json").write_text(
        json.dumps({"dataset_dir": "/x/kl-mlff-cpu/mat/step2_disp_force"}))
    cwd = tmp_path / "mat"
    (cwd / "step4_disp").mkdir(parents=True)
    try:
        g5._bundle_provenance(src, cwd, False)
        assert False, "foreign source 应被拒绝"
    except SystemExit:
        pass
    prov = g5._bundle_provenance(src, cwd, True)
    assert prov["foreign_source"] is True

