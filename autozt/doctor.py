# -*- coding: utf-8 -*-
"""doctor —— 配置预检（autozt doctor）：提交前把「配置层」的坑一次查清。

2026-09 的 295 材料实测里踩的坑，没有一个是计算调度的问题，全在配置/发现层，
而且都能在动手前查出来：
  · 项目 tf_*.yaml 里写的 max_jobs 以前被静默忽略 → 这里逐项目列出生效上限；
  · AUTOZT_COLLECT_CHUNK 只设在交互 shell、cron 里漏了 → 列出每个调优旋钮的
    生效值与来源（env / 配置 / 默认），只在 env 里的给出「写进 tf.yaml」提示；
  · 陈旧项目被屏蔽、又与正常项目材料同名 → 列出屏蔽项目及原因、跨项目同名材料；
  · 旧技能名改名遗留 → 按「旧名→新名」汇总，给出要改的文件；
  · monitor 守护死了没人拉 → 看 pid、auto_watch、crontab 保活；
  · AI 会话未走审批网关 → 报告当前网关策略与是否识别为 agent 会话。

纯本地：不连超算、不提交、不写任何文件（只读配置与本地目录）。
"""

import json
import os
import subprocess
import sys
import time

_LEVELS = ("error", "warn", "info")


def _finding(level, code, msg, fix=None, **extra):
    d = {"level": level, "code": code, "message": msg}
    if fix:
        d["fix"] = fix
    d.update({k: v for k, v in extra.items() if v not in (None, [], {})})
    return d


def _proj_name(t):
    from autozt import _seg_proj_name
    return _seg_proj_name(t) or ""


def check_tuning(cfg):
    from autozt import TUNING_KNOBS, tuning_value
    rows, finds = [], []
    for key, (env, default, _lo) in sorted(TUNING_KNOBS.items()):
        val, src = tuning_value(key, cfg, with_source=True)
        rows.append({"key": key, "env": env, "value": val, "source": src,
                     "default": default})
        if src.startswith("env:"):
            in_cfg = (cfg.get(key) is not None
                      or (isinstance(cfg.get("tuning"), dict)
                          and cfg["tuning"].get(key) is not None))
            if not in_cfg:
                finds.append(_finding(
                    "warn", "knob_env_only",
                    "%s=%s 只来自环境变量 %s；cron/monitor 保活进程里很容易漏设"
                    % (key, val, env),
                    fix="写进 tf.yaml：%s: %s（环境变量仍可临时覆盖）" % (key, val)))
    return rows, finds


def check_max_jobs(cfg, types):
    from autozt import _skill_max_jobs
    rows, seen = [], set()
    for t in types or []:
        key = t.get("key")
        proj = _proj_name(t)
        if (key, proj) in seen:
            continue
        seen.add((key, proj))
        skill_cap = _skill_max_jobs(cfg, t)
        pcap = t.get("_project_max_jobs")
        try:
            pcap = int(pcap) if pcap is not None else None
        except (TypeError, ValueError):
            pcap = None
        eff = min([c for c in (skill_cap, pcap) if c]) if (skill_cap or pcap) else None
        rows.append({"skill": key, "project": proj or "(global)",
                     "skill_cap": skill_cap, "project_cap": pcap,
                     "effective_cap": eff,
                     "note": ("项目级 %s 与技能级 %s 同时生效，先到先卡" % (pcap, skill_cap))
                     if pcap and skill_cap else ""})
    return rows


def check_blocked():
    from autozt import bootstrap as _b
    finds = []
    for path, reason in sorted(_b._CONFIG_SKILL_ERRORS.items()):
        finds.append(_finding("warn", "project_blocked",
                              "项目配置 %s 已被屏蔽：%s" % (path, reason),
                              fix="修好/迁移该配置，或移入 _archive，或加进 project_root_excludes"))
    for name, paths in sorted(_b._CONFIG_CONFLICTS.items()):
        finds.append(_finding("warn", "config_name_conflict",
                              "tf_%s.yaml 有 %d 份同名配置，全部屏蔽：%s"
                              % (name, len(paths), ", ".join(sorted(paths))),
                              fix="重命名其中一份或移入 _archive"))
    return finds


def check_aliases():
    from autozt import bootstrap as _b
    by_pair = {}
    for src, old, new in sorted(_b._SKILL_MIGRATIONS):
        by_pair.setdefault("%s→%s" % (old, new), []).append(src)
    finds = []
    for pair, srcs in sorted(by_pair.items()):
        finds.append(_finding("info", "legacy_skill_name",
                              "%d 个配置仍用旧技能名 %s（已在内存映射）" % (len(srcs), pair),
                              fix="确认物理目录后把这些文件里的技能名改成新名",
                              files=srcs))
    return finds


def check_same_names(cfg, types):
    """跨项目同名材料 + 与被屏蔽项目撞名的材料（纯本地发现）。"""
    from autozt import discover_local
    from autozt import bootstrap as _b
    by_name, by_base = {}, {}
    for t in types or []:
        lr = t.get("local_root")
        if not lr:
            continue
        try:
            _r, mats = discover_local(lr, tt=t.get("key"))
        except Exception:                       # noqa: BLE001
            continue
        proj = _proj_name(t) or "(global)"
        for m in mats:
            rp = os.path.realpath(m["lpath"])
            by_name.setdefault((t.get("key"), m["name"]), {})[rp] = proj
            by_base.setdefault(os.path.basename(m["name"]), set()).add(rp)
    finds = []
    dups = {k: v for k, v in by_name.items() if len(v) > 1}
    if dups:
        sample = sorted(dups.items())[:10]
        finds.append(_finding(
            "warn", "cross_project_same_name",
            "%d 个材料名在多个项目里重复（-p 写 basename 会歧义）" % len(dups),
            fix="-p 写 <项目名>/<完整名>，或用 --project 限定；AI 用 JSON 里的 id 字段",
            examples=["%s %s ← %s" % (k[0], k[1], ", ".join(sorted(set(v.values()))))
                      for k, v in sample]))
    blocked_hits = []
    for owner, _keys, name, _aliases in _b._CONFIG_BLOCKS:
        for lp in _b._block_candidates(owner):
            base = os.path.basename(lp)
            live = [p for p in by_base.get(base, ()) if not (
                p == owner or p.startswith(owner + os.sep))]
            if live:
                blocked_hits.append("%s（被屏蔽项目 %s 与正常项目同名）" % (base, name))
    if blocked_hits:
        finds.append(_finding(
            "warn", "blocked_project_shadow",
            "%d 个材料与被屏蔽项目里的材料同名；按 basename 操作时会先撞上屏蔽项目"
            % len(blocked_hits),
            fix="把陈旧项目加进 project_root_excludes（或移入 _archive）；-p 用 <项目名>/<完整名>",
            examples=sorted(blocked_hits)[:10]))
    return finds


def check_clusters(cfg, types):
    from autozt import _PKG_ROOT, pkg_setting_path, _load_yaml_file
    import glob
    import re
    finds, names = [], set()
    for t in types or []:
        h = t.get("hpc")
        if isinstance(h, str) and h:
            names.add(h)
    for n in sorted(names):
        path = pkg_setting_path(n + ".yaml")
        if not path:
            finds.append(_finding("error", "cluster_config_missing",
                                  "技能/项目引用的集群 %s 没有 setting/%s.yaml" % (n, n),
                                  fix="照 setting/template-cluster.yaml 建一份"))
            continue
        c = _load_yaml_file(path) or {}
        declared = str(c.get("partition") or c.get("queue") or "").strip()
        if not declared:
            continue
        bad = set()
        for f in glob.glob(os.path.join(_PKG_ROOT, "setting", n, "templates", "*.tpl")):
            try:
                txt = open(f, encoding="utf-8").read()
            except OSError:
                continue
            for p in re.findall(r"^#SBATCH\s+--partition=(\S+)", txt, re.M):
                if p != declared:
                    bad.add("%s: %s" % (os.path.basename(f), p))
        if bad:
            finds.append(_finding("warn", "partition_mismatch",
                                  "集群 %s 声明分区 %s，但模板里写的是：%s"
                                  % (n, declared, "; ".join(sorted(bad))),
                                  fix="改模板分区或 setting/%s.yaml 的 partition" % n))
    return finds


def check_monitor(cfg):
    from autozt import _watch_files, _watch_running_pid
    pid, _pf = _watch_running_pid(cfg)
    pidfile, logfile = _watch_files(cfg)
    cron = None
    try:
        r = subprocess.run(["crontab", "-l"], capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=10)
        cron = (r.returncode == 0 and "tf-watch-keepalive" in (r.stdout or ""))
    except (OSError, subprocess.SubprocessError):
        cron = None
    info = {"pid": pid, "alive": bool(pid), "log": logfile,
            "auto_watch": bool(cfg.get("auto_watch")), "cron_keepalive": cron,
            "auto_advance": bool(cfg.get("auto_advance"))}
    finds = []
    if cfg.get("auto_advance") and not pid:
        finds.append(_finding("warn", "monitor_not_running",
                              "auto_advance 开着但后台 monitor 没在跑，流水线不会自己推进",
                              fix="autozt monitor -d；并 autozt monitor --install 装 cron 保活"))
    if pid and not cfg.get("auto_watch") and cron is False:
        finds.append(_finding("info", "monitor_no_keepalive",
                              "monitor 在跑，但既没开 auto_watch 也没装 cron 保活；进程死了不会被拉起",
                              fix="autozt monitor --install"))
    return info, finds


def check_fast_path(cfg, types):
    """单材料快路径（-p 只采目标材料、agent 读进度文件）对每个技能是否可用。"""
    from autozt import narrow_keys, narrow_max_age
    keys = sorted({t.get("key") for t in types if t.get("key")})
    ok, doc, why = narrow_keys(cfg, keys)
    cov = (doc or {}).get("coverage") or {}
    now = time.time()
    rows = [{"skill": k, "narrow": k in ok,
             "full_collect_age_s": int(now - float(cov[k])) if cov.get(k) else None,
             "reason": why.get(k, "")} for k in keys]
    finds = []
    off = [r["skill"] for r in rows if not r["narrow"]]
    if off and narrow_max_age(cfg) > 0:
        finds.append(_finding("info", "fast_path_off",
                              "%d 个技能暂不能只采目标材料（-p 命令会全量采集）：%s"
                              % (len(off), ", ".join(off[:8])),
                              fix="让 monitor 跑一轮（或 autozt -tt <技能> list --refresh）"
                                  "记下整技能采集；有效期 narrow_max_age（默认 7200 秒）"))
    return {"narrow_max_age": narrow_max_age(cfg), "skills": rows}, finds


def check_gate(cfg):
    from autozt import agent_detect, agent_gate_policy
    actor, src = agent_detect(cfg)
    return {"policy": agent_gate_policy(cfg), "actor": actor or "", "detected_by": src or "",
            "stdin_tty": _isatty()}


def _isatty():
    try:
        return bool(sys.stdin and sys.stdin.isatty())
    except (AttributeError, ValueError):
        return False


def run_doctor(cfg, types, deep=True):
    report = {"schema": "autozt-doctor/1"}
    finds = []
    knobs, f = check_tuning(cfg)
    report["tuning"] = knobs
    finds += f
    report["max_jobs"] = check_max_jobs(cfg, types)
    finds += check_blocked()
    finds += check_aliases()
    if deep:
        finds += check_same_names(cfg, types)
    finds += check_clusters(cfg, types)
    mon, f = check_monitor(cfg)
    report["monitor"] = mon
    finds += f
    report["agent_gate"] = check_gate(cfg)
    fp, f = check_fast_path(cfg, types)
    report["fast_path"] = fp
    finds += f
    report["collector_transport"] = "stdin（payload 不占 argv，无 E2BIG 上限）"
    finds.sort(key=lambda d: _LEVELS.index(d["level"]))
    report["findings"] = finds
    report["counts"] = {lv: sum(1 for d in finds if d["level"] == lv) for lv in _LEVELS}
    return report


def cmd_doctor(cfg, types, json_out=False, strict=False, deep=True):
    rep = run_doctor(cfg, types, deep=deep)
    if json_out:
        print(json.dumps(rep, ensure_ascii=False, indent=2))
    else:
        c = rep["counts"]
        print("autozt doctor：error %d / warn %d / info %d（纯本地检查，未连超算）"
              % (c["error"], c["warn"], c["info"]))
        print("\n[max_jobs 生效上限]（技能级=全局 task_types.<技能>.max_jobs；"
              "项目级=项目 tf_*.yaml 里写的；两者同时生效）")
        for r in rep["max_jobs"]:
            print("  %-16s %-32s 技能级=%-5s 项目级=%-5s 生效=%s"
                  % (r["skill"], r["project"][:32], r["skill_cap"], r["project_cap"],
                     r["effective_cap"]))
        print("\n[调优旋钮]（环境变量 > tf.yaml > 默认）")
        for r in rep["tuning"]:
            print("  %-16s = %-6s ← %s" % (r["key"], r["value"], r["source"]))
        m = rep["monitor"]
        print("\n[monitor] %s；auto_advance=%s auto_watch=%s cron 保活=%s"
              % (("在跑 PID %s" % m["pid"]) if m["alive"] else "未运行",
                 m["auto_advance"], m["auto_watch"],
                 {True: "已装", False: "未装", None: "未知"}[m["cron_keepalive"]]))
        fp = rep["fast_path"]
        print("[单材料快路径]（-p 只采目标材料；并发计数由进度文件+提交账本补；"
              "有效期 %ss）" % fp["narrow_max_age"])
        for r in fp["skills"]:
            print("  %-16s %s" % (r["skill"], ("可用（整技能采集 %s 秒前）" % r["full_collect_age_s"])
                                  if r["narrow"] else ("全量采集：%s" % r["reason"])))
        g = rep["agent_gate"]
        print("[agent 网关] 策略=%s；本会话%s"
              % (g["policy"], ("识别为 agent（%s，%s）" % (g["actor"], g["detected_by"]))
                 if g["actor"] else "未识别为 agent（破坏性命令在非交互终端仍需批准）"))
        if rep["findings"]:
            print("\n[发现]")
        for d in rep["findings"]:
            print("  %-5s %-26s %s" % (d["level"], d["code"], d["message"]))
            for ex in (d.get("examples") or [])[:5]:
                print("        · %s" % ex)
            if d.get("files") and len(d["files"]) <= 5:
                for ff in d["files"]:
                    print("        · %s" % ff)
            if d.get("fix"):
                print("        → %s" % d["fix"])
    if rep["counts"]["error"] or (strict and rep["counts"]["warn"]):
        return 1
    return 0
