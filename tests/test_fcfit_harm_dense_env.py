"""HARM_DENSE 两套逻辑的交互：默认必须强制 1；显式写 0 也必须强制 1（不允许跑完才在检查步失败）。

本会话加了“按方法强制”（_pheasy_env），另一会话加了“按配置打开”（_apply_harm_dense）。
两者并存时，显式 pheasy_harm_dense=0 不能把强制值关掉。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skill" / "_common" / "fcfit"))
import fc_fit_driver as d  # noqa: E402

SPARSE = ("LASSO", "ALASSO", "RFE-OLS", "RFE-OLS-TSQR", "ARDR", "RVM")


def _env(method):
    return d._pheasy_env(method, "fit", 8, 10)


def test_default_forces_on_for_every_sparse_method():
    for m in SPARSE:
        env = _env(m)
        assert env.get("PHEASY_HARM_DENSE") == "1", (m, env.get("PHEASY_HARM_DENSE"))


def test_explicit_off_cannot_disable_it():
    for m in SPARSE:
        env = _env(m)
        d._apply_harm_dense(env, {"pheasy_harm_dense": "0"})
        assert env.get("PHEASY_HARM_DENSE") == "1", (m, env.get("PHEASY_HARM_DENSE"))


def test_explicit_on_still_on():
    for m in SPARSE:
        env = _env(m)
        d._apply_harm_dense(env, {"pheasy_harm_dense": "1"})
        assert env.get("PHEASY_HARM_DENSE") == "1"


def test_dense_methods_have_no_flag():
    for m in ("OLS", "RIDGE"):
        env = _env(m)
        assert "PHEASY_HARM_DENSE" not in env, (m, env.get("PHEASY_HARM_DENSE"))
