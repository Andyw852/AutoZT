"""AutoZT 冒烟测试（不需要集群）。pytest -q 即可；需要真集群的用例打 cluster 标记。"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "autozt")


def run(*args):
    return subprocess.run([sys.executable, PROG] + list(args),
                          capture_output=True, text=True, cwd=ROOT)


def test_version_banner_names_the_tool():
    p = run("--version")
    out = (p.stdout or "") + (p.stderr or "")
    assert "AutoZT" in out and "autozt" in out
    assert "taskflow" not in out and "(tf)" not in out


def test_help_uses_new_cli_name():
    out = (lambda p: (p.stdout or "") + (p.stderr or ""))(run("--help"))
    assert "autozt" in out and "用法：tf" not in out


def test_skills_and_schema():
    assert run("skills").returncode == 0
    assert run("schema", "--strict").returncode == 0


def test_package_has_version():
    sys.path.insert(0, ROOT)
    import autozt
    assert autozt.__version__.count(".") >= 1


def test_publication_metadata_present():
    for f in ("LICENSE", "pyproject.toml", "CITATION.cff", "CHANGELOG.md",
              "docs/CONFIGURING.md", "setting/template-cluster.yaml"):
        assert os.path.isfile(os.path.join(ROOT, f)), f
