"""[0020] fit-fc-thermal 默认 FIT_ENGINE=pheasy → pheasy 从可选升为必需依赖（symfc/hiphive 仍可选）。纯本地。"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_fit_fc_thermal_requires_pheasy():
    import yaml
    req = yaml.safe_load(
        (ROOT / "skill" / "fit-fc-thermal" / "skill.yaml").read_text(encoding="utf-8")
    )["requires"]
    req_py = [str(x) for x in (req.get("python") or [])]
    opt_py = [str(x) for x in (req.get("optional_python") or [])]
    assert "pheasy" in req_py, req_py
    assert "pheasy" not in opt_py, opt_py
    assert "symfc" in opt_py and "hiphive" in opt_py, opt_py
