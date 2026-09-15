"""i18n 单测（CI 可跑）：英文选择、回退、提示表。"""
import importlib
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run(code, lang=None):
    env = dict(os.environ)
    env.pop("PHONOAGENT_LANG", None)
    env.pop("LANG", None)
    if lang:
        env["PHONOAGENT_LANG"] = lang
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       env=env, cwd=ROOT)
    return p.stdout.strip()


def test_default_is_chinese_fallback():
    out = _run("import sys;sys.path.insert(0,'.');"
               "from phonoagent import i18n;print(i18n.t('中文','English'))")
    assert out == "中文"


def test_en_switches_to_english():
    out = _run("import sys;sys.path.insert(0,'.');"
               "from phonoagent import i18n;print(i18n.t('中文','English'))", lang="en_US.UTF-8")
    assert out == "English"


def test_hint_table_falls_back_and_translates():
    code = ("import sys;sys.path.insert(0,'.');"
            "from phonoagent import preflight as pf;"
            "print(pf.hint('resources','中文提示')[:6]);"
            "print(pf.hint('nope','中文提示'))")
    assert _run(code).splitlines() == ["中文提示", "中文提示"]          # zh：原样 + 未知种类回退
    en = _run(code, lang="en").splitlines()
    assert en[0].startswith("hint:") and en[1] == "中文提示"           # en：有译文；未知种类仍回退


def test_lang_env_alone_is_enough():
    out = _run("import sys;sys.path.insert(0,'.');"
               "from phonoagent import i18n;print(i18n.is_en())", lang="en")
    assert out == "True"
