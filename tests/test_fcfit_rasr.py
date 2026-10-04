# -*- coding: utf-8 -*-
"""fit-fc-thermal PHEASY_RASR: 'auto' must mean BHH on a 2D slab, none on bulk.

Regression guard for the bug that burned a full pheasy run on the Mn-In-Se
datasets: the skill default was PHEASY_RASR=BHH.  RASR is mandatory for 2D (the
ZA branch has to be quadratic near Gamma), but on a bulk crystal the same
constraints drove the MnIn2Se4 189-atom pheasy OLS fit from 0.62 % force error
(min frequency -0.03 THz, stable) to 8.9 % (min frequency -1.57 THz, spurious
imaginary modes).  The cutoff was irrelevant, so 'auto' keys off DIM.
"""
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SK = os.path.join(ROOT, "skill", "fit-fc-thermal")
ENG = os.path.join(ROOT, "skill", "_common", "fcfit")


def _load():
    sys.path.insert(0, SK)
    spec = importlib.util.spec_from_file_location(
        "fc_fit_driver_rasr_t", os.path.join(ENG, "fc_fit_driver.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_auto_is_bhh_on_2d():
    drv = _load()
    assert drv.resolve_rasr("auto", "2d") == "BHH"
    assert drv.resolve_rasr("AUTO", "2D") == "BHH"


def test_auto_is_none_on_bulk():
    drv = _load()
    for dim in ("3d", "3D", "", None, "auto"):
        assert drv.resolve_rasr("auto", dim) == "none"


def test_explicit_value_is_untouched():
    drv = _load()
    for v in ("BHH", "BH", "H", "none", "", "false", "off"):
        for dim in ("2d", "3d"):
            assert drv.resolve_rasr(v, dim) == v


def test_gen_default_is_auto():
    """The gen default must stay in step with the driver resolution."""
    sys.path.insert(0, SK)
    spec = importlib.util.spec_from_file_location(
        "gen_step1_fit_rasr_t", os.path.join(ENG, "gen_step1_fit.py"))
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    assert gen.SPEC["PHEASY_RASR"][0] == "auto"
