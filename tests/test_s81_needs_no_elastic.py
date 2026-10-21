# -*- coding: utf-8 -*-
"""[0065] S8.1 不直接等 S6（金属关 dpt 后不被 IBRION=6 卡住）；S6 出厂 NFREE=2 / EDIFF=1E-6 / KSPACING=0.03。"""
import re
from pathlib import Path

import yaml

KE = Path(__file__).resolve().parents[1] / "skill" / "ke-dft-cpu"


def _steps():
    man = yaml.safe_load((KE / "skill.yaml").read_text(encoding="utf-8"))
    out = {s["label"]: s for s in man.get("steps", [])}
    for g in (man.get("optional_steps") or {}).values():
        out.update({s["label"]: s for s in g.get("steps", [])})
    return out


def test_s81_needs():
    st = _steps()
    assert "step6_elastic" not in st["S8.1_bt2"]["needs"]
    assert "step8.2_dpt" in st["S8.1_bt2"]["needs"]          # 半导体仍经 S8.2 间接等 S6
    assert "step6_elastic" in st["S8.2_dpt"]["needs"]


def test_s81_does_not_read_elastic():
    src = (KE / "step8.1_boltztrap" / "gen_step11_boltztrap.py").read_text(encoding="utf-8")
    assert "step6_elastic" not in src


def test_elastic_defaults():
    src = (KE / "step6_elastic" / "gen_step2_elastic.py").read_text(encoding="utf-8")
    val = lambda k: re.search(r'^%s\s*=\s*"([^"]+)"' % k, src, re.M).group(1)
    assert (val("NFREE"), val("EDIFF"), val("KSPACING")) == ("2", "1E-6", "0.03")
