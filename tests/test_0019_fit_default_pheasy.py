"""[0019] κL 力常数拟合默认 pheasy + ALASSO（kl-dft-cpu / kl-mlff-cpu / kl-mlff-gpu / fit-fc-thermal
统一）；旧模板留下的 FIT_ENGINE = auto、PHEASY_FIT_METHOD = auto 按它解析，不再报错退出。纯本地。"""
import contextlib
import importlib.util
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FCFIT = ROOT / "skill" / "_common" / "fcfit"
OPT = ROOT / "skill" / "_common" / "opt"


def _load(path, name):
    for d in (str(Path(path).parent), str(OPT), str(ROOT / "skill" / "_common")):
        if d not in sys.path:
            sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


def test_auto_resolves_to_pheasy_alasso():
    g1 = _load(FCFIT / "gen_step1_fit.py", "g1_0019")
    eng, meth, notes = g1.resolve_auto_recipe("auto", "auto")
    assert (eng, meth) == ("pheasy", "ALASSO") and len(notes) == 2
    assert g1.resolve_auto_recipe("", "")[:2] == ("pheasy", "ALASSO")
    assert g1.resolve_auto_recipe(None, None)[2] == []                     # 没写：不打印解析说明
    assert g1.resolve_auto_recipe("phono3py", "auto")[:2] == ("phono3py", "ALASSO")
    assert g1.resolve_auto_recipe("Pheasy", "LASSO") == ("pheasy", "LASSO", [])
    assert g1.normalize_pheasy_method(g1.resolve_auto_recipe("auto", "auto")[1]) in g1.PHEASY_METHODS


def test_templates_default_to_pheasy_alasso():
    stepconf = _load(OPT / "stepconf.py", "stepconf_0019")
    for rel in ("kl-dft-cpu/templates/step5_fc/step.conf",
                "kl-mlff-cpu/templates/step3_fc/step.conf",
                "kl-mlff-gpu/templates/step3_fc/step.conf",
                "fit-fc-thermal/templates/step1_fit/step.conf"):
        params = {k: v for k, v, _ in
                  stepconf.parse((ROOT / "skill" / rel).read_text(encoding="utf-8"))["params"]}
        assert params["FIT_ENGINE"].lower() == "pheasy", rel
        assert params["PHEASY_FIT_METHOD"].upper() == "ALASSO", rel


def test_legacy_auto_conf_generates_pheasy_alasso_job(tmp_path, monkeypatch):
    """旧 kl-dft-cpu 项目 step.conf 里的 auto/auto：gen 不再退出，按 pheasy + ALASSO 生成（GPU 模板）。"""
    import json
    import numpy as np
    g1 = _load(FCFIT / "gen_step1_fit.py", "g1_0019_e2e")
    import fc_common
    ds, sk = tmp_path / "mat" / "step4_disp", tmp_path / "mat" / "fit-fc-thermal"
    ds.mkdir(parents=True)
    sk.mkdir()
    ds.joinpath("POSCAR").write_text("sc\n1.0\n3 0 0\n0 3 0\n0 0 3\nX\n1\nDirect\n0 0 0\n")
    pos = "".join("%g %g %g\n" % (i / 2, j / 2, k / 2)
                  for i in range(2) for j in range(2) for k in range(2))
    ds.joinpath("SPOSCAR").write_text("sc\n1.0\n6 0 0\n0 6 0\n0 0 6\nX\n8\nDirect\n" + pos)
    dd = np.random.default_rng(1).normal(0, 0.03, (6, 8, 3))
    dd[-1] = 0
    np.save(ds / "dataset_disps.npy", dd)
    np.save(ds / "dataset_forces.npy", -dd)
    sk.joinpath("step.conf").write_text("[params]\nSTEP = step1_fit\nFIT_ENGINE = auto\n"
                                        "PHEASY_FIT_METHOD = auto\nCUT3_CANDIDATES = off\n")
    for t in ("submit_fcfit_p3py", "submit_fcfit_pheasy", "submit_fcfit_hiphive",
              "submit_fcfit_pheasy_gpu"):
        (sk / (t + ".tpl")).write_text((FCFIT / "templates" / (t + ".tpl")).read_text(encoding="utf-8"))
    monkeypatch.chdir(sk)
    monkeypatch.setattr(fc_common, "resolve_submit", lambda base, kind: sk / (kind + ".tpl"))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        g1.main()
    cfg = json.loads((sk / "step1_fit" / "fit_config.json").read_text())
    assert cfg["engine"] == "pheasy" and cfg["pheasy_method"] == "ALASSO"
    assert "FIT_ENGINE=auto -> pheasy" in buf.getvalue()
    assert "--gres=gpu" in (sk / "step1_fit" / "submit.sh").read_text()      # pheasy-gpu 走 GPU 模板


def test_pheasy_is_a_required_dependency():
    """[0019] 默认引擎是 pheasy → 三个技能的 pheasy 从可选依赖升为必需；symfc 仍可选（phono3py 路）。"""
    import yaml
    for rel in ("kl-dft-cpu", "kl-mlff-cpu", "kl-mlff-gpu"):
        req = yaml.safe_load((ROOT / "skill" / rel / "skill.yaml").read_text(encoding="utf-8"))["requires"]
        req_py = [str(x) for x in (req.get("python") or [])]
        opt_py = [str(x) for x in (req.get("optional_python") or [])]
        assert "pheasy" in req_py, (rel, req_py)
        assert "pheasy" not in opt_py, (rel, opt_py)
        assert "symfc" in opt_py, (rel, opt_py)

