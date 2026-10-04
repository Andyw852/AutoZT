#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""template_drift.py —— 项目级模板副本与技能模板的漂移审计（V137，只读，不改任何文件）。

背景：ke-dft-cpu 的 skill.yaml 写的是 template_dir: "."，autozt init 会把技能目录下**所有** *.tpl / *.conf
按相对路径复制进 project_setting/templates/（如 templates/step5_dielect/incar_dfpt_2d.tpl），之后 find_asset
永远优先用这份副本。于是技能模板后来的修复到不了老项目：
  · P1_Mo-MoS2 / P1_Al-AlN / P1_Zn-ZnO：副本是 8 月的 LPEAD=.TRUE.（技能 09-14 起改成 .FALSE.），DFPT 出 NaN；
  · Mo2S3：S1 的残余力 0.022 eV/Å，而现行 S1 模板是 EDIFFG = −0.005（副本多半是更松的旧版）。
autozt 自带的漂移检查只比 {{占位符}}，键值变了它看不出来。

本工具对 <root> 下每个 templates/ 目录里的 *.tpl：
  · 找技能里同一相对路径的模板（--skill-dir，默认本技能目录）；内容相同 -> 不报；
  · 不同：逐键比 INCAR 式的 KEY = value（忽略注释、空白、大小写），列出不同的键，关键键标 ★；
  · 再用 git 历史判断副本是不是技能模板的某个**历史版本原样**（--repo，默认本仓库）：
      [STALE]  是历史版本、没人改过 -> 删掉副本即可让技能模板生效（或 init -f 刷新，先备份）；
      [CUSTOM] 不是任何历史版本 -> 有人手改过；把要保留的改动移到 step.conf 的 [incar]，再删副本。
  · 提交脚本模板（submit*.tpl）是站点相关的，默认不比键值（--all 也查）；
  · [V159] 但 submit_amset*.tpl 一律查环境：写死旧环境（amset_clean = 0.4.19）或没有 {{AMSET_ENV}} 占位符
      -> [ENV]（有 [ENV] 退出码也是 1）。P1_Mo-MoS2 的 S4_wave 就是被 8 月的项目级副本钉在 amset_clean 上
      （作业报 EnvironmentNameNotFound）。

用法：
    python template_drift.py <root> [--skill-dir DIR] [--repo DIR] [--all] [--csv out.csv]
    例：python template_drift.py /mnt/d/tf_data/test_TE
退出码：有 [STALE]/[CUSTOM] 且含 ★ 键 -> 1，否则 0。
"""
import argparse
import csv
import os
import subprocess
import sys
import re
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SKILL_DIR = HERE.parent
# 与 ke_common.OLD_AMSET_ENVS / amset_submit_envs 同口径（本工具只读、不 import ke_common，免得拉起 dim_common 等依赖）
OLD_AMSET_ENVS = ("amset_clean",)
_ACTIVATE_RE = re.compile(r"\b(?:conda|mamba|micromamba)\s+(?:activate|run\s+(?:-n|--name))\s+[\"']?([^\s\"';&|)]+)"
                          r"|\bsource\s+activate\s+[\"']?([^\s\"';&|)]+)")
CRITICAL = {"LPEAD", "IBRION", "LEPSILON", "LCALCEPS", "EDIFF", "EDIFFG", "ISIF", "ISYM", "PREC", "LREAL",
            "ENCUT", "NSW", "POTIM", "ISMEAR", "SIGMA", "LORBIT", "NEDOS", "ICHARG", "LWAVE", "KSPACING",
            "IVDW", "LASPH", "ALGO", "NELM", "ADDGRID"}


def incar_keys(text):
    """KEY = value（去注释，值压空白、转大写）。{{占位符}} 原样保留。"""
    out = {}
    for ln in text.splitlines():
        ln = ln.split("#", 1)[0].split("!", 1)[0].strip()
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        k = k.strip().upper()
        if not k or " " in k:
            continue
        out[k] = " ".join(v.split()).upper()
    return out


def key_diff(proj_text, skill_text):
    a, b = incar_keys(proj_text), incar_keys(skill_text)
    rows = []
    for k in sorted(set(a) | set(b)):
        if a.get(k) != b.get(k):
            rows.append((k, a.get(k), b.get(k)))
    return rows


def _git(repo, *args):
    try:
        r = subprocess.run(["git", "-C", str(repo)] + list(args), capture_output=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout if r.returncode == 0 else None


def history_match(repo, skill_file, data):
    """副本是否等于技能模板 skill_file 在 git 历史里的某个版本。
    返回 (提交短哈希, 日期)；不是任何历史版本返回 None；没有 git / 文件不在仓库里返回 False。"""
    import hashlib
    top = _git(repo, "rev-parse", "--show-toplevel")
    if top is None:
        return False
    top = Path(top.decode().strip())
    try:
        rel = os.path.relpath(str(Path(skill_file).resolve()), str(top)).replace(os.sep, "/")
    except ValueError:
        return False
    if rel.startswith(".."):
        return False
    want = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()   # git 的 blob id
    log = _git(top, "log", "--all", "--format=%h %ad", "--date=short", "--", rel)
    if not log:
        return None
    for line in log.decode().splitlines():
        h, d = line.split(" ", 1)
        blob = _git(top, "rev-parse", "%s:%s" % (h, rel))
        if blob and blob.decode().strip() == want:
            return h, d
    return None


def submit_env_issue(text):
    """[V159] submit_amset*.tpl 的环境问题：返回说明文字；没问题返回 None。"""
    found = []
    for ln in text.splitlines():
        st = ln.strip()
        if st and not st.startswith("#"):
            for m in _ACTIVATE_RE.finditer(st):
                n = m.group(1) or m.group(2)
                if n and not n.startswith("$") and n not in found:
                    found.append(n)
    old = [n for n in found if n in OLD_AMSET_ENVS]
    if old:
        return "写死了 AMSET 0.4.19 的旧环境 %s（作业会 EnvironmentNameNotFound，环境还在时会悄悄用 0.4.19）" % "/".join(old)
    if "{{AMSET_ENV}}" not in text:
        return ("没有 {{AMSET_ENV}} 占位符（写死 %s）—— 不跟 step.conf 的 AMSET_ENV 走" % "/".join(found)
                if found else "没有 {{AMSET_ENV}} 占位符，也没有 conda activate 行")
    return None


def audit(root, skill_dir=SKILL_DIR, repo=None, include_submit=False):
    repo = Path(repo) if repo else Path(skill_dir)
    rows = []
    for tdir in sorted(p for p in Path(root).rglob("templates") if p.is_dir()):
        for f in sorted(tdir.rglob("*.tpl")):
            rel = f.relative_to(tdir)
            if rel.name.startswith("submit_amset"):
                why = submit_env_issue(f.read_text(encoding="utf-8", errors="ignore"))
                if why:
                    rows.append({"kind": "ENV", "file": str(f), "skill": "（集群 setting/<集群>/templates/submit_amset.tpl）",
                                 "note": why + "；改名为 submit_amset.tpl.stale-<日期> 让集群模板生效，或把环境名改成 "
                                 "{{AMSET_ENV}}；副本修改于 %s"
                                 % time.strftime("%Y-%m-%d", time.localtime(f.stat().st_mtime)),
                                 "critical": True, "diff": []})
            if not include_submit and rel.name.startswith("submit"):
                continue
            sk = Path(skill_dir) / rel
            if not sk.is_file():
                continue
            pdata, sdata = f.read_bytes(), sk.read_bytes()
            if pdata == sdata:
                continue
            diff = key_diff(pdata.decode("utf-8", "ignore"), sdata.decode("utf-8", "ignore"))
            hm = history_match(repo, sk, pdata)
            if hm:
                kind, note = "STALE", "= 技能模板 %s（%s）的原样旧版" % (hm[0], hm[1])
            elif hm is None:
                kind, note = "CUSTOM", ("不是技能模板在本仓库历史里的任何版本（有手改，或比仓库历史更早的旧快照 —— "
                                        "看下面的键差和副本日期判断）")
            else:
                kind, note = "DIFF", "读不到 git 历史，无法判断是旧版还是手改"
            note += "；副本修改于 %s" % time.strftime("%Y-%m-%d", time.localtime(f.stat().st_mtime))
            rows.append({"kind": kind, "file": str(f), "skill": str(sk), "note": note,
                         "critical": any(k in CRITICAL for k, _, _ in diff),
                         "diff": diff})
    return rows


def _fmt_diff(diff):
    out = []
    for k, a, b in diff:
        out.append("%s%s：副本 %s / 技能 %s" % ("★" if k in CRITICAL else "", k,
                                            "未设" if a is None else a, "未设" if b is None else b))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root")
    ap.add_argument("--skill-dir", default=str(SKILL_DIR))
    ap.add_argument("--repo", default=None, help="技能所在 git 仓库（默认 --skill-dir 所在仓库）")
    ap.add_argument("--all", action="store_true", help="连 submit*.tpl 也查")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args(argv)
    rows = audit(a.root, a.skill_dir, a.repo, a.all)
    if not rows:
        print("没有与技能模板不同的项目级模板副本（root=%s）" % a.root)
        return 0
    for r in rows:
        print("[%s] %s\n       %s\n       技能模板：%s" % (r["kind"], r["file"], r["note"], r["skill"]))
        if r["kind"] == "ENV":
            continue
        for ln in _fmt_diff(r["diff"]) or ["（键值相同，只有注释/格式不同）"]:
            print("         - " + ln)
    n_env = sum(r["kind"] == "ENV" for r in rows)
    rows_t = [r for r in rows if r["kind"] != "ENV"]
    n_stale = sum(r["kind"] == "STALE" for r in rows_t)
    n_crit = sum(r["critical"] for r in rows_t)
    if n_env:
        print("\n★ %d 份 submit_amset 模板的 AMSET 环境有问题（[ENV]）：写死的环境与 step.conf 的 AMSET_ENV 不同时，"
              "S4/S4b/S8/S8.4 的 gen 会拒绝生成；没有占位符的，换集群或换环境时不会跟着变" % n_env)
    print("\n共 %d 份副本与技能模板不同：[STALE] %d（可直接删掉，让技能模板生效）、其余 %d 需人工确认；含 ★ 关键键 %d 份"
          % (len(rows_t), n_stale, len(rows_t) - n_stale, n_crit))
    if a.csv:
        with open(a.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["kind", "file", "skill", "note", "critical", "diff"])
            for r in rows:
                w.writerow([r["kind"], r["file"], r["skill"], r["note"], r["critical"],
                            "；".join(_fmt_diff(r["diff"]))])
        print("写出 %s" % a.csv)
    return 1 if (n_crit or n_env) else 0


if __name__ == "__main__":
    sys.exit(main())
