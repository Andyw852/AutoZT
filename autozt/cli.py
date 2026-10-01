# -*- coding: utf-8 -*-
"""cli —— 命令入口（17_cli）。main + 分发 + dry-run。
对外接口：main。"""

import os

from autozt import i18n as _i18n  # noqa: E402
import sys
import re
import json
import time
import shlex
import hashlib
import base64
import collections
import functools
import itertools
import subprocess
import tempfile
import threading
import socket
import argparse
import glob
import math
import random
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

# -*- coding: utf-8 -*-
# 17_cli —— main() 入口与命令分发
#
# 本分片由 autozt/__init__.py 装配器在单一命名空间里按顺序执行；
# 函数之间的引用按名字解析（与原单文件一致），分片间无需 import。
# 内容清单（按原文件行号）：
#   L7078  main
#   （数据簇 _dbg_t/_state_cache_*/collect_data/apply_exclude/filter_projs/
#    filter_status/status_spec_has_scancel/_snapshot 已抽成真模块 autozt/data.py）

def _dry_run_steps_for(cmd, m, jb):
    from autozt import find_step_soft
    """--dry-run：某材料在命令 cmd、步骤筛选 jb 下实际会动的步骤。
    返回步骤列表；None 表示"整材料级"（rerun/clean 无 -j 时）。"""
    if jb:
        s = find_step_soft(m, jb)
        return [s] if s is not None else []
    if cmd == "retry":
        return [s for s in m["steps"] if s["kind"] == "FAIL"]
    if cmd in ("start", "advance"):
        ready = m.get("actives")
        if ready is None:
            ready = [m.get("active")] if m.get("active") is not None else []
        return [s for s in ready if s and s["kind"] in ("TODO", "PREP")]
    if cmd == "stop":
        return [s for s in m["steps"] if s.get("job")]
    if cmd == "fetch":
        return [s for s in m["steps"] if s["kind"] == "OK"]
    return None   # rerun / clean：无 -j 时整材料级


def _dry_run_report(cfg, data, cmd, projs, jobs):
    from autozt import find_material
    """--dry-run：打印真实目标对象（复用 _dry_run_steps_for 的语义）。"""
    mats = []
    if projs:
        for pj in projs:
            try:
                mats.append(find_material(data, pj)[1])
            except SystemExit:
                print("  材料 %s  （未能解析）" % pj)
    else:
        mats = [m for t in data["types"] for m in t["materials"]]
    shown = 0
    for m in mats:
        for jb in jobs:
            steps = _dry_run_steps_for(cmd, m, jb)
            if steps is None:
                print("  材料 %s  （整材料 %s）" % (m["name"], cmd))
                shown += 1
                continue
            if not steps:
                print("  材料 %s  %s（无需操作）"
                      % (m["name"], ("步骤 %s " % jb) if jb else ""))
                continue
            for s in steps:
                print("  材料 %s  步骤 %s [%s]  ->  %s"
                      % (m["name"], s.get("label", s.get("name", "?")),
                         s.get("kind", "?"), s.get("dir", "?")))
                shown += 1
    if shown == 0:
        print("  （共 %d 个材料，无任何目标）" % len(mats))


def normalize_monitor_command(command, positional, restart=False):
    positional = list(positional)
    if command == "restart":
        return "monitor", positional, True
    if command == "watch":
        command = "monitor"
    if command == "monitor" and "restart" in positional:
        positional.remove("restart")
        restart = True
    return command, positional, restart


# 取值型选项：找裸命令词时跳过它们后面的值（-p agent 里的 agent 是材料名，不是子命令）
_ROUTE_VALUE_FLAGS = {'-c', '--config', '-tt', '-p', '-j', '-job', '-status', '--status',
                      '-x', '--exclude', '--host', '-u', '--user', '-proj', '--project',
                      '--expect-state', '--dataset', '--poscar', '--root', '--cluster',
                      '--set', '--out'}

# 带 -p 时可以只采目标材料的命令（进度文件对该技能有足够新的整技能覆盖时）。
# 这些命令只动/只看 -p 指定的材料；并发闸门缺的那部分由 busy_baseline 补。
# stop/rerun/clean（破坏性）与 status（会顺带 auto-fetch 全部材料）保持全量采集。
_NARROW_CMDS = {"start", "retry", "fetch", "advance", "list", "json", "diagnose",
                "conf", "dir", "prove"}


def _bare_index(argv, name):
    skip = False
    for idx, tok in enumerate(argv):
        if skip:
            skip = False
            continue
        if tok in _ROUTE_VALUE_FLAGS:
            skip = True
            continue
        if tok == name:
            return idx
    return None


def route_subcommand(argv):
    """agent / mcp 子命令走各自的稳定 JSON 面，在主参数解析和任何采集之前分流。

    以前只有 bin/autozt 与 python -m autozt 做了这层拦截；pip 安装生成的 autozt
    命令直接进 main()：`autozt agent ...` 报"不是命令"，`autozt mcp` 要先把全部
    材料 ssh 采集一遍才启动。返回退出码；不是这两个子命令返回 None。"""
    argv = list(argv)
    i_mcp = _bare_index(argv, "mcp")
    i_agent = _bare_index(argv, "agent")
    if i_mcp is not None and (i_agent is None or i_mcp < i_agent):
        from autozt import mcp as _mcp
        return _mcp.main(argv) or 0
    if i_agent is not None:
        from autozt import agent_cli as _agent
        return _agent.main(argv)
    return None


def _errors_to_stdout(argv):
    """错误是否镜像到 stdout：--json、AUTOZT_JSON_ERRORS=1、或处在 AI 代理环境。

    很多代理/MCP 宿主只读 stdout（或把 stderr 截断/丢弃），以前 `sys.exit("错误：…")`
    只写 stderr，代理看到的是"空输出 + 非 0"，只能瞎猜。"""
    if "--json" in argv:
        return True
    flag = (os.environ.get("AUTOZT_JSON_ERRORS") or "").strip().lower()
    if flag in ("0", "false", "no", "off"):
        return False
    if flag in ("1", "true", "yes", "on"):
        return True
    try:
        from autozt.agentgate import AGENT_ENV_MARKERS
    except Exception:
        AGENT_ENV_MARKERS = ("CLAUDECODE", "GEMINI_CLI", "CODEX_SANDBOX")
    return any(os.environ.get(k) for k in AGENT_ENV_MARKERS)


def main():
    """入口：sys.exit("错误：…") 这类字符串退出照旧（解释器把它写 stderr、退出码 1）；
    代理/--json 场景再往 stdout 镜像一行，免得只读 stdout 的宿主看到"空输出"。"""
    argv = sys.argv[1:]
    try:
        return _main()
    except SystemExit as e:
        code = e.code
        if code is not None and not isinstance(code, int) and _errors_to_stdout(argv):
            try:
                if "--json" in argv:
                    sys.stdout.write(json.dumps({"ok": False, "error": str(code), "rc": 1},
                                                ensure_ascii=False) + "\n")
                else:
                    sys.stdout.write("[autozt error] %s\n" % code)
                sys.stdout.flush()
            except Exception:
                pass
        raise


def _main():
    from autozt.bootstrap import reset_config_conflicts, reject_config_conflict_targets, filter_config_conflicts
    reset_config_conflicts()
    # JSON/MCP/agent callers and Windows consoles must see the same UTF-8 text.
    # line_buffering：输出重定向到文件/管道时也逐行落盘（以前块缓冲，
    # `autozt ... > log` 只能靠文件 mtime 猜它还在不在动）。
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
        except (AttributeError, ValueError):
            pass
    _routed = route_subcommand(sys.argv[1:])
    if _routed is not None:
        sys.exit(_routed)
    from autozt import cmd_doctor, cmd_progress, write_progress, _watch_running_pid, _watch_files, _seg_proj_name
    from autozt import EXAMPLE_CONFIG, JSON_SCHEMA, AUTOZT_VERSION, USAGE, USAGE_EN, _PKG_ROOT, _add_diag_codes, _json_changes, _json_errors_only, _json_paginate, _dbg_t, _state_cache_load, _state_cache_save, _summary_json, _watch_cron, _watch_daemon, _watch_ensure, _watch_stop, apply_exclude, apply_hide_done, apply_skills, auto_advance, auto_fetch, auto_recover_hung, cmd_adopt, cmd_auto, cmd_auto_project, cmd_auto_skill, cmd_clean, cmd_conf, cmd_diagnose, cmd_fetch, cmd_hpc, cmd_init, cmd_level, cmd_migrate_subdir, cmd_rerun, cmd_retry, cmd_skills, cmd_start, cmd_status, cmd_step_init, cmd_stop, cmd_summary, cmd_watch, collect_data, fill_local_dim, filter_projs, filter_status, find_material, find_step, find_uninited, get_types, load_config, merge_project_configs, render_table, status_spec_has_scancel, cmd_schema, cmd_skill_show, cmd_correct, cmd_correct_usage, cmd_history, history_record, cmd_prove, set_active_cfg, cmd_act, cmd_approve, agent_direct_gate, agent_audit, cmd_session
    if "--help-all" in sys.argv[1:]:
        # 英文帮助：AUTOZT_LANG=en 或命令行 --lang en
        _lang = (os.environ.get('AUTOZT_LANG') or '').lower()
        if '--lang' in sys.argv:
            try:
                _lang = sys.argv[sys.argv.index('--lang') + 1].lower()
            except IndexError:
                pass
        print(USAGE_EN if _lang.startswith('en') else USAGE)
        return
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        from autozt import QUICK_USAGE
        print(QUICK_USAGE)
        return
    if any(a in ("-V", "--version") for a in sys.argv[1:]):
        print("AutoZT (autozt) version %s" % AUTOZT_VERSION)
        pkg = _PKG_ROOT
        _prog = os.path.realpath(globals().get("_PROG_PATH") or _PKG_ROOT)
        print("程序: %s" % _prog)
        print("包根: %s（setting/ = 默认模板，skill/ = 技能脚本）" % pkg)
        return
    if len(sys.argv) == 1:   # patch_auto：纯 tf = 只报版本，不采集/不提交
        print("AutoZT (autozt) version %s" % AUTOZT_VERSION)
        print("程序: %s" % os.path.realpath(globals().get("_PROG_PATH") or _PKG_ROOT))
        print("")
        print("  autozt list      只读总表（不拉取、不提交）")
        print("  autozt summary   只读极简汇总（巡检省 token，见 AGENTS.md）")
        print("  autozt status    刷新状态 + auto-fetch + auto-advance")
        print("  autozt monitor   后台监控（-i 秒，-d 后台，restart 重做；watch 为旧名）")
        print("  tf -h        常用命令；tf --help-all 查看完整帮助")
        return
    p = argparse.ArgumentParser(prog="autozt")
    p.add_argument("-tt", dest="tt")
    p.add_argument("-p", dest="proj")
    p.add_argument("-proj", "--project", dest="project", metavar="项目",
                   help="只看/只动指定项目（tf_<项目>.yaml 的 <项目>，逗号分隔）："
                        "采集前就按项目裁剪，不用把几百个材料名列进 -p")
    p.add_argument("-j", "-job", dest="job")
    p.add_argument("-c", "--config")
    p.add_argument("--host")
    p.add_argument("-u", "--user")
    p.add_argument("-x", "--exclude", dest="exclude",
                   help="跳过指定项目（逗号分隔，全名或 basename）")
    p.add_argument("-status", "--status", dest="status_f", metavar="ST",
                   help="只保留含指定状态步骤的材料：done/running/pd/error/"
                        "waiting/scancel（逗号分隔），对任意命令生效；"
                        "如 -status scancel start 重跑全部被 stop 取消的")
    p.add_argument("-i", "--interval", type=int, default=300,
                   help="monitor 刷新间隔秒数（默认 300）")
    p.add_argument("-d", "--daemon", action="store_true",
                   help="monitor 放后台运行（日志 .tf_watch.log）")
    p.add_argument("--stop", action="store_true",
                   help="停止后台运行的 autozt monitor")
    p.add_argument("--install", action="store_true",
                   help="写入 crontab 保活（autozt monitor 重启后自动恢复）")
    p.add_argument("--uninstall", action="store_true",
                   help="移除 crontab 保活")
    p.add_argument("--restart", action="store_true",
                   help="重启后台监控（先停旧的再起新的）")
    p.add_argument("-f", "--force", action="store_true")
    p.add_argument("--all", dest="all_files", action="store_true",
                   help="fetch 时拉回每个步骤的全部文件（不只 fetch_files 清单）")
    p.add_argument("--hide-done", dest="hide_done", action="store_true",
                   help="状态表隐藏全部步骤都完成的项目")
    p.add_argument("--show-done", dest="show_done", action="store_true",
                   help="取消 hide_done 配置的隐藏效果")
    p.add_argument("--diff", dest="diff", action="store_true",
                   help="summary：与上次快照对比，无变化不输出（省 token，巡检用）")
    p.add_argument("--changes", dest="changes", action="store_true",
                   help="json：只输出相对上次快照的步骤级变更（结构化，批分析用）")
    p.add_argument("--refresh", dest="refresh", action="store_true",
                   help="list/summary 强制跳过本地状态缓存，重新 ssh 采集")
    p.add_argument("-y", "--yes", action="store_true")
    p.add_argument("--dry-run", dest="dry", action="store_true",
                   help="只打印将影响的对象/计划，不执行（advance/start/stop/retry/rerun/clean/fetch/adopt/migrate-subdir）")
    p.add_argument("--json", dest="json_out", action="store_true",
                   help="list/summary/status/dir 以 JSON 输出（机器可读）")
    p.add_argument("--schema", dest="schema", action="store_true",
                   help="json：打印字段 schema 说明")
    p.add_argument("--verify", dest="verify", action="store_true",
                   help="prove：逐份校验输入的 sha256 是否与档案一致")
    p.add_argument("--out", dest="out", metavar="文件",
                   help="session export：包输出路径（默认 cwd/tmp/session_<材料>_<时间>.tar.gz）")
    p.add_argument("--full", dest="full", action="store_true",
                   help="skill show：卡片里再给输入/参数/可接技能/纠错")
    p.add_argument("--strict", dest="strict", action="store_true",
                   help="schema：自描述有 [错误] 级问题时返回非零（CI/检查用）")
    p.add_argument("--since", dest="since", metavar="时间",
                   help="history：只看该时间之后的记录（如 2026-09-14 或 7d）")
    p.add_argument("-n", dest="last_n", type=int, metavar="N",
                   help="history：只显示最近 N 条（默认 40）")
    p.add_argument("--write", dest="hist_write", action="store_true",
                   help="history：采集一次并把变化写进 history.jsonl（默认只读）")
    p.add_argument("-clean", "--clean", dest="clean", action="store_true")
    p.add_argument("--purge-config", dest="purge_config", action="store_true",
                   help="clean：连 project_setting 一起删（默认保留，重算需 autozt init）")
    p.add_argument("--from-skill", dest="from_skill", action="store_true",
                   help="rerun：忽略项目侧模板/step.conf，只用 skill 库出厂版生成")
    p.add_argument("--set", dest="sets", action="append", metavar="节.键=值",
                   help="conf：改本步 step.conf（如 --set incar.EDIFF=1E-6；"
                        "值留空=删除该键）")
    p.add_argument("--dataset", dest="reg_dataset", metavar="目录",
                   help="register：位移+力数据集目录（写进该技能的数据集参数，如 FIT_INPUT_DIR）")
    p.add_argument("--poscar", dest="reg_poscar", metavar="文件",
                   help="register：材料结构文件（缺省取 --dataset 目录里的 POSCAR）")
    p.add_argument("--root", dest="reg_root", metavar="目录",
                   help="register：材料建在哪个 project_roots 下（缺省第一个）")
    p.add_argument("--cluster", dest="reg_cluster", metavar="集群",
                   help="register：顺带切到该集群（setting/<集群>.yaml；本机用 local）")
    p.add_argument("--errors-only", dest="errors_only", action="store_true",
                   help="json：只保留含 FAIL 步骤的材料与 FAIL 步骤")
    p.add_argument("--limit", dest="limit", type=int, metavar="N",
                   help="json：只输出前 N 个材料（分页）")
    p.add_argument("--offset", dest="offset", type=int, metavar="M",
                   help="json：跳过前 M 个材料（配合 --limit 分页）")
    p.add_argument("--expect-state", dest="expect_state", metavar="状态",
                   help="start/retry 的前置条件（agent 执行计划用）：-p 目标步骤的当前状态"
                        "不在列表里就跳过该材料，如 --expect-state TODO,PREP,SCANCEL")
    p.add_argument("args", nargs="*")
    a, unknown = p.parse_known_args()  # v3.14：位置参数可穿插在选项中间
    a.args = list(a.args) + unknown    # （-p A B retry / -p A B -j 1 dir）
    if a.schema:   # v2.0：--schema 打印 json 字段说明（无需配置/采集）
        print(JSON_SCHEMA)
        return

    commands = {"status", "list", "summary", "advance", "start", "stop", "retry", "rerun",
                "json", "config", "dir", "fetch", "init", "clean", "watch", "mcp",
                "monitor", "restart", "help", "auto", "adopt", "migrate-subdir",
                "hpc", "skills", "conf", "level", "diagnose", "probe", "push",
                # v1.0（加技能友好化）：schema = 看技能自描述（io_schema/flow/corrections）
                #                      correct = 把 FAIL 诊断喂给 _corrections/ handler 库
                "schema", "skill", "correct", "history", "prove",
                # v1.0（P0-1）：act = agent 动作网关（风险分档+审计）；
                #                approve = 人工批准破坏性动作（仅交互终端）
                "act", "approve",
                # v1.0（P1-7）：session export = 把一个材料的操作历史/provenance/
                #                  审计打包成论文补充材料
                "session",
                # AI 接入：progress = 读本地进度文件（不采集、不连超算）；
                #          doctor   = 配置预检（不连超算）
                "progress", "doctor",
                # register = 新材料一次接入（建目录+POSCAR+init+可选集群/数据集）
                "register"}
    root, cmd, pos = None, "status", []
    for tok in a.args:  # v3.14：位置参数先收集，之后按"材料名/目录"消歧
        if tok == "help":
            from autozt import QUICK_USAGE
            print(QUICK_USAGE)
            return
        if tok in commands and cmd == "status":
            cmd = tok
        else:
            pos.append(tok)
    cmd, pos, a.restart = normalize_monitor_command(cmd, pos, a.restart)
    # 本地存在的目录 → 旧版 ROOT 语义；其余位置参数留给材料名消歧
    mat_toks = []
    for tok in pos:
        if os.path.isdir(os.path.expanduser(tok)):
            if root is not None:
                sys.exit(_i18n.t("错误：", "error: ") + "只能指定一个 ROOT（收到多个目录）。")
            root = tok
        else:
            mat_toks.append(tok)

    if cmd == "config":
        print(EXAMPLE_CONFIG)
        return
    from autozt.bootstrap import canonical_skill_name
    a.tt = canonical_skill_name(a.tt)
    if cmd == "probe":   # 只读探测：delegate 到 scripts/probe_jobs.py（不采集/不提交/不改文件）
        import subprocess as _sp
        _probe = os.path.join(_PKG_ROOT, "scripts", "probe_jobs.py")
        _mats = [m.strip() for m in (a.proj or "").split(",") if m.strip()] + mat_toks
        if not _mats:
            sys.exit(_i18n.t("错误：", "error: ") + "probe 需要 -p 材料名（如 tf -tt defect-dft-cpu -p Sn2Sb2Te5 probe）。")
        _probe_cfg, _ = load_config(a.config)
        merge_project_configs(apply_skills(_probe_cfg, verbose=True))
        reject_config_conflict_targets(",".join(_mats), a.tt)
        _argv = [sys.executable, _probe, "-p", ",".join(_mats)]
        if a.job:
            _argv += ["-j", a.job]
        if a.host:
            _argv += ["--host", a.host]
        sys.exit(_sp.call(_argv))
    if cmd == "push":   # git 提交：delegate 到 scripts/tf-git-push.sh（本地真实版 + 远端脱敏版）
        import subprocess as _sp
        _push = os.path.join(_PKG_ROOT, "scripts", "tf-git-push.sh")
        _msg = " ".join(mat_toks) if mat_toks else "update"
        sys.exit(_sp.call(["bash", _push, _msg]))
    if a.job and not a.proj and not mat_toks and cmd not in (
            "start", "stop", "retry", "rerun", "clean", "status"):
        sys.exit(_i18n.t("错误：", "error: ") + "-j 必须和 -p 一起用（start/stop/retry/rerun/clean "
                 "支持不带 -p，表示对全部材料只操作该步骤）。")

    cfg, cfg_path = load_config(a.config)
    cfg["_config_dir"] = (os.path.dirname(os.path.abspath(cfg_path))
                          if cfg_path else os.getcwd())
    cfg["_config_path"] = cfg_path
    set_active_cfg(cfg)    # v1.0：让 log_action 能把动作记进 history.jsonl
    if cmd in ("act", "approve"):   # v1.0（P0-1）：网关在**任何**采集之前短路
        sys.exit(cmd_act(cfg, sys.argv[1:]) if cmd == "act"
                 else cmd_approve(cfg, sys.argv[1:]))
    _gate = agent_direct_gate(cfg, cmd, sys.argv[1:])   # agent 会话直接调用 → 审计
    if _gate is not None:
        sys.exit(_gate)
    if cmd == "progress":   # 只读本地进度文件：不装技能、不扫项目、不连超算（秒级）
        sys.exit(cmd_progress(cfg, project=a.project, tt=a.tt,
                              material=",".join(filter(None, [a.proj] + mat_toks)) or None,
                              status=a.status_f, json_out=a.json_out, limit=a.limit))
    if cmd == "monitor" and (a.stop or a.install or a.uninstall
                             or (a.daemon and not a.restart)):
        # 控制类操作不需要项目配置。cron 保活每 10 分钟跑一次 monitor -d：以前先
        # 做完整的 apply_skills + merge_project_configs（9p 上 1–2 分钟、还会卡在
        # D 状态）才去看 pid——守护明明在跑也白扫一遍。现在先看 pid。
        if a.install:
            sys.exit(_watch_cron(True))
        if a.uninstall:
            sys.exit(_watch_cron(False))
        if a.stop:
            sys.exit(_watch_stop(cfg))
        _pid, _ = _watch_running_pid(cfg)
        if _pid:
            print("autozt monitor 已在后台运行 (PID %d)" % _pid)
            print("日志：%s（tail -f 查看）；停止：autozt monitor --stop"
                  % _watch_files(cfg)[1])
            return
    cfg = apply_skills(cfg, verbose=True)
    # 只读技能库的命令用不到项目配置，在扫描 project_roots 之前就分派：
    # merge_project_configs 要遍历全部项目目录，WSL 9p 上每次 1–2 分钟。以前
    # `autozt skills`/`schema`/`skill`——以及 agent skills/contract/research_plan、
    # MCP list_skills/describe_skill/research_plan（它们都调 `schema --json`）——
    # 每调一次都白扫一遍。这三个命令只读 cfg["_skills"]（apply_skills 已装好）。
    # 带了 -p 材料目标时不走捷径：显式目标必须先过下面的屏蔽项目检查（安全约定，
    # 见 tests/test_skill_aliases.py）。schema/skill 的位置参数是技能名，不是材料。
    _fast = not a.proj
    if cmd == "skills" and _fast and not mat_toks:
        return cmd_skills(cfg, tt=a.tt)
    if cmd == "schema" and _fast:   # v1.0：看技能自描述（纯本地、不采集、不提交）
        # 技能名既可用 -tt，也可直接当位置参数写：autozt schema band-dft-cpu
        _which = a.tt or (mat_toks[0] if mat_toks else None)
        sys.exit(cmd_schema(cfg, tt=_which, json_out=a.json_out, strict=a.strict))
    if cmd == "skill" and _fast:   # v1.0：技能卡片（论文图 2 的机器可读来源；纯本地）
        # tf skill show <技能> / tf skill show -tt <技能> / tf skill <技能> / tf skill
        _args = list(mat_toks)
        if _args and _args[0] in ("show", "list"):
            _args.pop(0)
        sys.exit(cmd_skill_show(cfg, which=a.tt or (_args[0] if _args else None),
                                json_out=a.json_out, full=a.full))
    if cmd == "history" and not a.hist_write and _fast and not mat_toks:
        # 不带目标的 history 只读 history.jsonl；带 -p 时仍走下面的屏蔽项目检查
        sys.exit(cmd_history(cfg, proj=None, tt=a.tt, since=a.since,
                             last_n=a.last_n or 40, json_out=a.json_out))
    cfg = merge_project_configs(cfg)

    def _types_for_block_check():   # 只在 basename 撞上屏蔽项目时才算（本地发现，有缓存）
        return get_types(cfg, tt=a.tt, quiet=True)
    reject_config_conflict_targets(a.proj, a.tt, types_fn=_types_for_block_check)
    if cmd == "monitor":
        reject_config_conflict_targets(",".join(mat_toks), a.tt,
                                       types_fn=_types_for_block_check)
    if cmd == "monitor":   # 控制类操作不采集状态，提前短路
        if a.install:
            sys.exit(_watch_cron(True))
        if a.uninstall:
            sys.exit(_watch_cron(False))
        if a.stop:
            sys.exit(_watch_stop(cfg))
        if a.restart:
            _watch_stop(cfg)
            _watch_daemon(a, mat_toks, root, cfg)
            return
        if a.daemon:
            _watch_daemon(a, mat_toks, root, cfg)
            return
    if cmd == "skills":   # 带 -p 目标时才会走到这里（已过屏蔽检查）
        return cmd_skills(cfg, tt=a.tt)
    if cmd == "schema":
        _which = a.tt or (mat_toks[0] if mat_toks else None)
        sys.exit(cmd_schema(cfg, tt=_which, json_out=a.json_out, strict=a.strict))
    if cmd == "skill":
        _args = list(mat_toks)
        if _args and _args[0] in ("show", "list"):
            _args.pop(0)
        sys.exit(cmd_skill_show(cfg, which=a.tt or (_args[0] if _args else None),
                                json_out=a.json_out, full=a.full))
    if cmd == "session" and not (a.proj or mat_toks):
        # 缺材料名不必采集（省一次 ssh）：直接给用法
        sys.exit(cmd_session(cfg, None, None))
    if cmd == "history" and not a.hist_write:
        # v1.0：直接读 history.jsonl（不采集、不连超算、不提交）。
        # 记录是自动的——任何一次真正采集都会追加；--write 时才先采集一轮再读。
        sys.exit(cmd_history(cfg,
                             proj=a.proj or (mat_toks[0] if mat_toks else None),
                             tt=a.tt, since=a.since, last_n=a.last_n or 40,
                             json_out=a.json_out))
    if a.host is not None:
        cfg["host"] = a.host or None
    if a.user:
        cfg["user"] = a.user
    types = get_types(cfg, tt=a.tt,
                      root_override=None if cmd in ("init", "register") else root,
                      quiet=(cmd in ("init", "register")))
    if a.project:   # --project：采集前按项目配置裁剪段（只采这个项目，快且不串项目）
        _wp = {x.strip() for x in a.project.split(",") if x.strip()}
        _avail = sorted({_seg_proj_name(t) for t in types if _seg_proj_name(t)})
        types = [t for t in types if _seg_proj_name(t) in _wp]
        if not types:
            sys.exit(_i18n.t("错误：", "error: ") + "没有项目 %s%s（可用项目：%s）。"
                     % (a.project, ("（在技能 %s 下）" % a.tt) if a.tt else "",
                        ", ".join(_avail[:40]) or "（无）"))
    if cmd == "doctor":   # 配置预检：纯本地，不采集、不连超算
        sys.exit(cmd_doctor(cfg, types, json_out=a.json_out, strict=a.strict))
    if cmd not in ("watch", "monitor"):
        _watch_ensure(cfg)   # v1.10：auto_watch 时顺带确保后台监控在跑
    if cmd in ("status", "json") and not types and not a.tt:
        sys.exit(_i18n.t("错误：", "error: ") + "没有任何任务类型"
                 "（在全局 tf.yaml 或项目 project_setting/tf_*.yaml 里定义）。")
    # v1.1：-tt 指定的类型有骨架但无项目段时 types 为空——前面已打印引导
    # 提示，这里放行，按空表/无目标处理（不算错误）。

    if cmd == "register":   # 新材料接入：纯本地（建目录/init/写 step.conf），不连超算
        _gate = agent_direct_gate(cfg, cmd, sys.argv[1:])
        if _gate is not None:
            sys.exit(_gate)
        from autozt import cmd_register
        sys.exit(cmd_register(cfg, types, a.proj, a.tt, dataset=a.reg_dataset,
                              poscar=a.reg_poscar, root=a.reg_root,
                              cluster=a.reg_cluster, sets=a.sets, step=a.job,
                              json_out=a.json_out))

    if cmd == "init" and not a.job:  # 项目配置初始化：纯本地，不连超算
        _gate = agent_direct_gate(cfg, cmd, sys.argv[1:])   # P0-1：init 也进审计
        if _gate is not None:
            sys.exit(_gate)
        if a.proj and mat_toks:
            print("提示：多个材料要用逗号分隔，如  -p %s,%s"
                  % (a.proj, ",".join(mat_toks)))
        sys.exit(cmd_init(cfg, types, a.proj, tt=a.tt,
                          name=(mat_toks[0] if mat_toks else root),
                          force=a.force, yes=a.yes))

    if cmd == "level":  # patch_level：设/查计算级别（纯本地改 step.conf）
        _arg = mat_toks[0] if mat_toks else None
        sys.exit(cmd_level(cfg, types, a.tt, a.proj, _arg))

    _force_advance = False   # autonow：仅 auto on 这一次允许推进
    if cmd == "auto":   # v1.5：一键开关 auto_advance（纯本地改 tf.yaml）
        # autonow2：从位置参数里挑 on/off 当开关，其余位置参数当材料名，
        # 并标记为已消费。原来固定取 mat_toks[0] 且不消费，导致
        #   `autozt auto on`           -> "on" 落到后面被当材料名解析而报错
        #   `tf -tt ke <材料> auto on` -> 材料名被当成了 on/off 参数
        _AUTO_WORDS = ("on", "off", "1", "0", "true", "false", "resume",
                       "\u5f00", "\u5173")
        _arg = next((x for x in mat_toks
                     if str(x).strip().lower() in _AUTO_WORDS), None)
        _rest = [x for x in mat_toks
                 if str(x).strip().lower() not in _AUTO_WORDS]
        _proj = a.proj or (",".join(_rest) if _rest else None)
        if _arg == "resume" and (not _proj or not a.tt):
            print(_i18n.t("错误：", "error: ") + "auto resume 必须用 -tt 指定技能和 -p 指定项目。")
            sys.exit(1)
        mat_toks, a.proj = [], _proj
        if _proj:       # v1.9.9：带 -p（或位置参数）就只改这些材料/技能
            _rc = cmd_auto_project(cfg, types, _proj, a.tt, _arg, dry=a.dry)
        elif a.tt:      # patch_auto：-tt 不带材料 = 该技能下全部材料
            _rc = cmd_auto_skill(cfg, types, a.tt, _arg, dry=a.dry)
        elif a.project:  # --project 不带 -tt：该项目挂的每个技能（types 已按项目裁剪）
            _rc = 0
            for _k in sorted({t["key"] for t in types}):
                _rc |= cmd_auto_skill(cfg, types, _k, _arg, dry=a.dry)
        else:
            if a.dry and _arg is not None:
                print("[dry-run] 将把全局 %s 的 auto_advance 改为 %s（未写入）。"
                      % (cfg.get("_config_path") or "tf.yaml", _arg))
                sys.exit(0)
            _rc = cmd_auto(cfg, _arg)
        if a.dry:       # --dry-run：只列出会改的文件，不写、不推进
            sys.exit(_rc)
        # patch_auto_now：原实现三条路都 sys.exit，走不到下面的 status 分支，
        # 而 auto_advance() 只在 status / watch 里调用 —— 于是 auto on 只翻
        # 开关不干活，还得再手敲一次 tf。这里改成：on 且全部成功 -> 不退出，
        # 落到 status 分支跑一轮采集+推进。off / 查询 / 有失败 -> 保持原行为。
        # autonow：on 且全部成功 → 不退出，置标志落到下面的 status 分支，
        # 当场跑一轮采集 + 推进。off / 无参查询 / 有失败 → 原样，只翻开关。
        if _rc != 0 or str(_arg or "").strip().lower() not in (
                "on", "1", "true", "开", "resume"):
            sys.exit(_rc)
        print("auto_advance 已开，下面立刻提交可开始的步骤"
              "（只想看不提交：autozt list）。")
        _force_advance = True
        cmd = "status"

    if cmd == "adopt":  # v1.5：接管手工整理的技能子目录（内部自做采集）
        sys.exit(cmd_adopt(cfg, types,
                           a.proj or (mat_toks[0] if mat_toks else None),
                           a.yes, a.dry, tt=a.tt))

    if cmd == "hpc":    # v1.7：指定项目分配到指定超算（纯本地改 hpc.yaml）
        sys.exit(cmd_hpc(cfg, types,
                         [x.strip() for x in (a.proj or "").split(",")
                          if x.strip()],
                         mat_toks[0] if mat_toks else None, a.tt, a.yes))

    # v3：含 local_root 的类型走本地发现 + 按项目超算分组采集；其余远端扫描
    # v3.1：同 key 的多个段（项目配置）合并进同一个类型条目
    import time as _time
    _t0 = _time.time()
    # -p 指定了材料时，收集前先把无关的段过滤掉——大体系（数千材料）下全量
    # 采集会逐材料读 yaml + ssh，极慢。skill_subdir 段的材料名 =
    # basename(local_root)，可零成本精确匹配；其余段原样保留（回退全量，行为不变）。
    _types_all = list(types)
    if a.proj:
        _want = {x.strip() for x in a.proj.split(",") if x.strip()}
        # 共享(体系级)布局：local_root 是项目目录（无 POSCAR），材料名 = 相对路径
        # 前缀，无法用 basename 精确匹配 → 保留该段回退全量收集。
        types = [t for t in types
                 if (not (t.get("skill_subdir") and t.get("local_root"))
                     or os.path.basename(os.path.realpath(t["local_root"])) in _want
                     or ("%s/%s" % (_seg_proj_name(t),
                                    os.path.basename(os.path.realpath(t["local_root"]))))
                     in _want
                     or not os.path.isfile(os.path.join(
                         os.path.realpath(t["local_root"]), "POSCAR")))]
    # 单材料快路径：-p 指定材料时，共享布局的段只采目标材料（以前每条
    # `start -p X` 都要把整个技能几百个材料在 9p 上解析 + ssh 探测一遍）。
    # 前提是进度文件对该技能有足够新的整技能采集（narrow_max_age）——并发闸门里
    # 没采到的材料的占用从进度文件 + 提交账本补（busy_baseline），否则回退全量。
    _nkeys, _ndoc = set(), None
    if a.proj and not mat_toks and not root:
        from autozt import narrow_keys
        _nkeys, _ndoc, _nwhy = narrow_keys(cfg, sorted({t["key"] for t in _types_all}))
        if cmd in _NARROW_CMDS and not a.clean:
            _toks = [x.strip() for x in a.proj.split(",") if x.strip()]
            for t in types:
                if t["key"] in _nkeys and t.get("local_root"):
                    t["_narrow"] = _toks
        if _nwhy and os.environ.get("AUTOZT_DEBUG_TIME"):
            print("[采集] 以下技能全量采集（不能只采目标材料）：%s"
                  % "; ".join("%s: %s" % kv for kv in sorted(_nwhy.items())),
                  file=sys.stderr)
    _narrowed = any(t.get("_narrow") for t in types)
    # patch_state_cache：list/summary 只读命令优先读本地缓存（跳过 ssh 采集），
    # --refresh 或 AUTOZT_CACHE_TTL=0 强制刷新；会改状态的命令一律现采。
    from autozt import tuning_value
    _ttl = tuning_value("cache_ttl", cfg)
    _cached = None
    if cmd in ("list", "summary") and not a.refresh and _ttl > 0:
        _cached = _state_cache_load(cfg, types, a.tt, root, _ttl)
    if _cached is not None:
        data = _cached
        if os.environ.get("AUTOZT_DEBUG_TIME"):
            print("[缓存] 命中本地状态缓存，跳过 ssh 采集（--refresh 强制刷新）",
                  file=sys.stderr)
    else:
        data = collect_data(cfg, types)
        fill_local_dim(cfg, data, types)
        if cmd in ("list", "summary") and not _narrowed:   # 只采了目标材料的不能当全量缓存
            _state_cache_save(cfg, data, types, a.tt, root)
        # v1.0（W5–8）：把本轮的状态转移追加进 history.jsonl。这样**任何技能**
        # 加进来就自动有历史，不必各技能自己写日志；走缓存那一支不重复记录。
        try:
            _n_hist = history_record(cfg, data)
            if _n_hist and os.environ.get("AUTOZT_DEBUG_TIME"):
                print("[history] 记录 %d 条状态转移" % _n_hist, file=sys.stderr)
        except Exception:
            pass
        # 机器可读进度文件（.tf_progress.json）：AI 之后用 autozt progress /
        # agent progress 直接读，不用再 ssh。局部采集按材料合并，不截断全量视图。
        filter_config_conflicts(data)
        _whole = not (a.proj or a.project or root or mat_toks)
        write_progress(cfg, data, writer="cli", full_scope=_whole and not a.tt,
                       covered_keys=[t["key"] for t in data["types"]] if _whole else None,
                       observed=_t0)
    if _nkeys and data is not _cached:
        # 没采到的同技能材料（被裁剪的段/材料）的在跑作业数，给 max_jobs 闸门加底数
        from autozt import busy_baseline
        data["_busy_baseline"] = busy_baseline(
            cfg, _ndoc, data, _nkeys,
            {_seg_proj_name(t): t.get("_from") for t in _types_all if _seg_proj_name(t)})
    _dbg_t("状态采集（ssh+远端扫描）", _t0)

    filter_config_conflicts(data)   # 包含旧缓存；必须先于状态过滤/动作分派
    if a.project:   # 缓存命中时 types 裁剪不生效，这里再按材料的 project 字段兜底
        _wp = {x.strip() for x in a.project.split(",") if x.strip()}
        for _t in data["types"]:
            _t["materials"] = [m for m in _t["materials"]
                               if (m.get("project") or "") in _wp]
    if a.proj:
        # -p 写错材料名以前静默返回"（没有任何任务类型）"/空 JSON 且 rc=0，
        # agent 会当成"该材料没事"。现在全部 -p 都找不到 → 报错退出（rc≠0）；
        # 部分找不到 → 警告后照常处理找得到的。
        from autozt import _name_matches
        _all_m = [m for t in data["types"] for m in t["materials"]]
        def _hit(x):
            # 材料名 / basename / <项目>/<名>；clean 等还接受体系目录（-p C20 → C20/*）
            return any(_name_matches(m, x) or (m.get("name") or "").startswith(x + "/")
                       or (m.get("qualified_name") or "").startswith(x + "/")
                       for m in _all_m)
        _miss = [x.strip() for x in a.proj.split(",") if x.strip() and not _hit(x.strip())]
        # 本地存在、只是本轮没采到（ssh 不通/采集组失败）的不算"找不到"——只告警，
        # 否则网络抖一下 agent 就会以为材料没了。
        _roots = [os.path.realpath(os.path.expanduser(str(r)))
                  for r in (cfg.get("project_roots") or [])]
        _lr = {os.path.basename(os.path.realpath(t["local_root"]))
               for t in _types_all if t.get("local_root")}

        def _local(x):
            return x in _lr or any(os.path.isfile(os.path.join(r, x, "POSCAR"))
                                   for r in _roots)
        _uncollected = [x for x in _miss if _local(x)]
        if _uncollected:
            print(_i18n.t("警告：", "warning: ")
                  + "材料 %s 在本地存在，但本轮没有采集到它的状态（见上方采集告警，"
                    "多半是 ssh/集群不可达）。" % ", ".join(_uncollected), file=sys.stderr)
        _miss = [x for x in _miss if x not in _uncollected]
        if _miss:
            import difflib
            # -p 预过滤可能已把 skill_subdir 段裁掉了：候选名再从未裁剪的段补
            _cands = sorted({m["name"] for m in _all_m}
                            | {os.path.basename(m["name"]) for m in _all_m}
                            | {os.path.basename(os.path.realpath(t["local_root"]))
                               for t in _types_all if t.get("local_root")})
            _hint = ["%s→%s" % (x, ",".join(difflib.get_close_matches(x, _cands, n=2,
                                                                   cutoff=0.6)))
                     for x in _miss if difflib.get_close_matches(x, _cands, n=1, cutoff=0.6)]
            _msg = (_i18n.t("找不到材料：", "material not found: ") + ", ".join(_miss)
                    + ((a.tt and "（技能 %s 内）" % a.tt) or "")
                    + (("；相近：" + "; ".join(_hint)) if _hint else "")
                    + "。新材料先 autozt -tt <技能> -p <材料> init（或 MCP register_material）。")
            if len(_miss) == len([x for x in a.proj.split(",") if x.strip()]):
                sys.exit(_i18n.t("错误：", "error: ") + _msg)
            print(_i18n.t("警告：", "warning: ") + _msg, file=sys.stderr)
    apply_exclude(data, a.exclude)   # v3.11：-x 跳过指定项目

    incl_sc = status_spec_has_scancel(a.status_f)   # v1.4
    if a.status_f:   # v1.4：-status 只保留含指定状态步骤的材料
        n0 = sum(len(t["materials"]) for t in data["types"])
        filter_status(data, a.status_f)
        n1 = sum(len(t["materials"]) for t in data["types"])
        # Machine-readable branches must remain a single JSON document.  The
        # human table keeps the helpful match-count line.
        if not a.json_out:
            print("（-status %s：%d/%d 个材料匹配）" % (a.status_f, n1, n0))

    # v3.14：多项目——-p 支持逗号分隔，位置参数里的材料名自动并入
    # v1.0：多步骤——-j 支持逗号分隔，位置参数里的步骤名/label 自动并入
    projs = [x.strip() for x in (a.proj or "").split(",") if x.strip()]
    jobs = [x.strip() for x in (a.job or "").split(",") if x.strip()]
    if mat_toks:
        names, bases, stepnames = set(), {}, set()
        for t in data["types"]:
            for m in t["materials"]:
                names.add(m["name"])
                bases.setdefault(os.path.basename(m["name"]), m["name"])
                for s in m["steps"]:
                    stepnames.add(s["name"])
                    stepnames.add(s["label"])
        for tok in mat_toks:
            if tok in names or tok in bases:
                projs.append(tok)
            elif tok in stepnames:
                jobs.append(tok)
            else:
                import difflib
                close = difflib.get_close_matches(
                    tok, sorted(names | set(bases) | stepnames | commands),
                    n=1, cutoff=0.6)
                sys.exit(_i18n.t("错误：", "error: ") + "'%s' 不是命令、材料或步骤%s"
                         % (tok, ("，你是不是想 '%s'？" % close[0]) if close else "。"))
    jobs = jobs or [None]

    if a.expect_state:
        # agent 执行计划的前置条件：计划可能来自几十分钟前的进度文件，这里按**现采**
        # 状态逐个核对，不符合的材料跳过（不会把已算完的步骤重交、把在跑的作业 retry 掉）。
        if cmd not in ("start", "retry") or not projs:
            sys.exit(_i18n.t("错误：", "error: ") + "--expect-state 只用于带 -p 的 start/retry。")
        from autozt import agent_event, find_step_soft
        from autozt.agent_protocol import status_kinds
        _want_k = status_kinds(a.expect_state)
        _keep = []
        for pj in projs:
            try:
                _t, _m = find_material(data, pj)
            except SystemExit as _e:
                print("跳过 %s：%s" % (pj, _e.code))
                agent_event("skipped", material=pj, reason=str(_e.code))
                continue
            _bad = []
            for jb in jobs:
                _s = find_step_soft(_m, jb) if jb else _m.get("active")
                _k = (_s or {}).get("kind")
                if _k not in _want_k:
                    _bad.append({"step": (_s or {}).get("label") or jb, "kind": _k})
            if _bad:
                print("跳过 %s：当前状态 %s 不满足 --expect-state %s（计划已过期）。"
                      % (_m["name"], ", ".join("%s=%s" % (b["step"], b["kind"]) for b in _bad),
                         a.expect_state))
                agent_event("skipped", material=pj,
                            id=_m.get("qualified_name") or _m.get("name"),
                            name=_m.get("name"), steps=_bad,
                            expected=sorted(_want_k), reason="state_changed")
                continue
            _keep.append(pj)
        if not _keep:
            print("--expect-state：没有材料满足前置条件，未执行任何操作。")
            return
        projs = _keep

    # v2.0：--dry-run 排练。破坏性/有副作用命令只打印将影响的对象，不执行。
    # 按命令语义打印【真实】目标：retry=FAIL 步，start=就绪步，stop=有作业步，
    # fetch=已完成可拉回步；rerun/clean 无 -j 时是整材料级（不再笼统"全部步骤"）。
    _eff_cmd = "clean" if (a.clean or cmd == "clean") else cmd
    if a.dry and _eff_cmd in ("advance", "start", "stop", "retry", "rerun", "clean",
                              "fetch"):
        print("【dry-run】命令 '%s' 将影响以下对象（未执行任何变更、未提交作业）："
              % _eff_cmd)
        _dry_run_report(cfg, data, _eff_cmd, projs, jobs)
        return

    if cmd == "monitor":   # 前台监控（控制标志已在前面短路）
        _ov = {}                        # v1.8：自动重载时命令行覆盖照旧生效
        if a.host is not None:
            _ov["host"] = a.host or None
        if a.user:
            _ov["user"] = a.user
        cmd_watch(cfg, types, projs, a.exclude, a.interval,
                  tt=a.tt, root=root, overrides=_ov, project=a.project)
        return

    if a.clean or cmd == "clean":  # 删除生成物回到 PREP（留 POSCAR）
        fails = 0
        for pj in (projs or [None]):
            for jb in jobs:
                fails += 0 if cmd_clean(cfg, data, pj, jb, a.yes,
                                        purge_config=a.purge_config) == 0 else 1
        sys.exit(0 if fails == 0 else 1)

    if cmd == "init":  # -p MAT -j STEP init：只生成该步骤输入，不提交
        if len(projs) > 1:
            sys.exit(_i18n.t("错误：", "error: ") + "步骤级 init 一次只支持一个材料。")
        fails = 0
        for jb in jobs:
            fails += 0 if cmd_step_init(cfg, data, projs[0] if projs else None,
                                        jb, a.force) == 0 else 1
        sys.exit(fails)

    if cmd == "history":   # v1.0：--write 已先采集记录一轮，这里再读出来
        sys.exit(cmd_history(cfg, proj=a.proj or None, tt=a.tt, since=a.since,
                             last_n=a.last_n or 40, json_out=a.json_out))
    if cmd == "list":   # v3.23：只读总览表格；不 auto_fetch/auto_advance（绝不提交）
        if a.hide_done or (cfg.get("hide_done") and not a.show_done):
            apply_hide_done(data)
        if a.json_out:
            print(json.dumps(_add_diag_codes(data), ensure_ascii=False, indent=2))
        else:
            render_table(data)
        return
    if cmd == "diagnose":   # v2.0：一键结构化诊断（只读；默认输出 FAIL 步）
        if not projs:
            sys.exit(_i18n.t("错误：", "error: ") + "diagnose 需要 -p 材料（如 tf -p C24/qHPC24 diagnose）。")
        _outs = [cmd_diagnose(cfg, data, pj, jobs[0]) for pj in projs]
        print(json.dumps(_outs[0] if len(_outs) == 1 else _outs,
                         ensure_ascii=False, indent=2))
        return
    if cmd == "correct":   # v1.0：把 FAIL 诊断喂给 _corrections/ handler 库
        if not projs:
            sys.exit(cmd_correct_usage())
        _fails = 0
        for pj in projs:
            _fails += cmd_correct(cfg, data, pj, jobs[0], yes=a.yes, dry=a.dry)
        sys.exit(1 if _fails else 0)
    if cmd == "prove":   # v1.0：看这一步"结果是怎么来的"（只读本地档案，不提交）
        if not projs:
            sys.exit(_i18n.t("错误：", "error: ") + "prove 需要 -p 材料（如 tf -p C24/qHPC24 prove）。")
        _rc = 0
        for pj in projs:
            _rc |= cmd_prove(cfg, data, pj, jobs[0], json_out=a.json_out,
                             verify=a.verify)
        sys.exit(_rc)
    if cmd == "session":   # v1.0（P1-7）：会话导出（只读本地，打论文补充材料包）
        _sargs = list(mat_toks)
        if _sargs and _sargs[0] in ("export", "bundle"):
            _sargs.pop(0)
        _sp = a.proj or (_sargs[0] if _sargs else None)
        sys.exit(cmd_session(cfg, data, _sp, job=jobs[0], out=a.out,
                             since=a.since, json_out=a.json_out))
    if cmd == "advance":
        # 一次性推进：采集并回传已完成结果，然后推进就绪步骤；不改写
        # 全局 auto_advance 配置。与 `auto on`（持久化开关）严格区分。
        # Respect an explicit material scope.  `auto on` persists scope through
        # project settings; one-shot advance must filter the collected snapshot
        # directly so it cannot advance unrelated materials.
        if projs:
            filter_projs(data, projs)
        _t1 = _time.time()
        auto_fetch(cfg, data)
        _dbg_t("advance：auto-fetch 拉回", _t1)
        _t1 = _time.time()
        auto_advance(cfg, data, force=True)
        _dbg_t("advance：一次性推进", _t1)
        if a.hide_done or (cfg.get("hide_done") and not a.show_done):
            apply_hide_done(data)
        for pj in (projs or [None]):
            for jb in jobs:
                cmd_status(cfg, data, pj, jb)
        return
    if cmd == "summary":   # 只读极简汇总；不 auto_fetch/auto_advance（绝不提交）
        if a.hide_done or (cfg.get("hide_done") and not a.show_done):
            apply_hide_done(data)
        _sp = None
        if a.diff:   # 快照按过滤范围分开存，不同 -tt/-status/-x/-p 互不干扰
            import hashlib as _hashlib
            _p = ",".join(sorted(x.strip() for x in (a.proj or "").split(",") if x.strip()))
            _scope = ",".join([a.tt or "", a.status_f or "", a.exclude or "", _p]
                              + (["project=" + a.project] if a.project else []))
            _h = _hashlib.md5(_scope.encode("utf-8")).hexdigest()[:8]
            _sp = os.path.join(cfg.get("_config_dir") or os.getcwd(),
                               ".tf_summary_%s.txt" % _h)
        if a.json_out:
            print(json.dumps(_summary_json(data), ensure_ascii=False, indent=2))
        else:
            cmd_summary(data, diff=a.diff, state_path=_sp)
        return
    if cmd == "status":
        if a.json_out:
            print(json.dumps(_add_diag_codes(data), ensure_ascii=False, indent=2))
            return
        _t1 = _time.time()
        auto_fetch(cfg, data)   # 算完的步骤自动保存到本地 result/
        _dbg_t("auto-fetch 拉回", _t1)
        # autonow：只有从 auto on 落过来时才推进；裸 tf / autozt status / tf list
        # 仍是 fixte⑤ 的只读语义（不提交任务）。
        if _force_advance:
            _t1 = _time.time()
            auto_recover_hung(cfg, data)   # v1.11: 挂死自动恢复（含 UCX 卡死），与 auto_advance 同轮
            _dbg_t("hang-check 挂死检测", _t1)
            _t1 = _time.time()
            auto_advance(cfg, data)
            _dbg_t("auto-advance 推进", _t1)
        if a.hide_done or (cfg.get("hide_done") and not a.show_done):
            apply_hide_done(data)   # v1.1：隐藏全部完成的项目
        for pj in (projs or [None]):
            for jb in jobs:
                cmd_status(cfg, data, pj, jb)
        if not projs:
            new = find_uninited(cfg)   # v3.11：新材料目录自动检测（v1.2 扫 project_roots）
            if new:
                print("发现 %d 个新材料目录未初始化：%s"
                      % (len(new), ", ".join(new)))
                print("→ autozt init 纳入管理；配 auto_advance: true 后下次 tf 自动开算")
    elif cmd == "conf":
        if not projs or not jobs or not jobs[0]:
            sys.exit(_i18n.t("错误：", "error: ") + "conf 需要 -p 材料 -j 步骤（如 tf -tt bd -p Mg2C60 -j 2 conf）。")
        fails = 0
        for pj in projs:
            for jb in jobs:
                fails += cmd_conf(cfg, data, pj, jb, a.sets)
        sys.exit(1 if fails else 0)
    elif cmd == "json":
        if a.schema:
            print(JSON_SCHEMA)
        else:
            if a.changes:                          # v2.0：只输出步骤级变更快照
                import hashlib as _hashlib
                _p = ",".join(sorted(x.strip() for x in
                                     (a.proj or "").split(",") if x.strip()))
                _scope = ",".join([a.tt or "", a.status_f or "", a.exclude or "",
                                   _p] + (["project=" + a.project] if a.project else []))
                _h = _hashlib.md5(_scope.encode("utf-8")).hexdigest()[:8]
                _sp = os.path.join(cfg.get("_config_dir") or os.getcwd(),
                                   ".tf_json_snapshot_%s.txt" % _h)
                _out = _json_changes(data, _sp)
                _out["schema_version"] = 2
                _out["tf_version"] = AUTOZT_VERSION
                print(json.dumps(_out, ensure_ascii=False, indent=2))
                return
            _out = {"schema_version": 2, "tf_version": AUTOZT_VERSION}
            _out.update(data)
            if a.errors_only:                       # v2.0：只留 FAIL 材料/步骤
                _out = _json_errors_only(_out)
            if a.limit is not None or a.offset is not None:   # v2.0：材料分页
                _out = _json_paginate(_out, a.offset or 0, a.limit)
            print(json.dumps(_out, ensure_ascii=False, indent=2))
    elif cmd == "dir":
        # 只输出路径本身，方便拼进 ssh/cd 命令：ssh jzzn "cd $(tf -p X -j 1 dir)"
        if a.json_out:
            _dirs = []
            if projs:
                for pj in projs:
                    t, m = find_material(data, pj)
                    for jb in jobs:
                        _dirs.append(find_step(m, jb)["dir"] if jb else m["path"])
            elif a.tt:
                _dirs.append(data["types"][0]["root"])
            else:
                sys.exit(_i18n.t("错误：", "error: ") + "dir 需要 -p（可配 -j）或 -tt 指定对象。")
            print(json.dumps(_dirs, ensure_ascii=False))
        elif projs:
            for pj in projs:
                t, m = find_material(data, pj)
                for jb in jobs:
                    print(find_step(m, jb)["dir"] if jb else m["path"])
        elif a.tt:
            print(data["types"][0]["root"])
        else:
            sys.exit(_i18n.t("错误：", "error: ") + "dir 需要 -p（可配 -j）或 -tt 指定对象。")
    elif cmd == "migrate-subdir":
        if not a.tt:
            sys.exit(_i18n.t("错误：", "error: ") + "migrate-subdir 需要 -tt 指定迁哪个技能"
                     "（如 tf -tt band migrate-subdir）。")
        fails = 0
        for pj in (projs or [None]):
            fails += cmd_migrate_subdir(cfg, data, pj, a.yes, a.dry)
        sys.exit(0 if fails == 0 else 1)
    else:
        fails = 0
        # -p A,B,C 逐个 start 时共用一个闸门：以前每个材料都从同一份提交前的计数
        # 起算，批量能冲过 max_jobs。
        from autozt import _SkillGate
        _start_gate = _SkillGate(cfg, data) if cmd == "start" else None
        for pj in (projs or [None]):
            for jb in jobs:
                if cmd == "start":
                    fails += cmd_start(cfg, data, pj, jb, a.force,
                                       incl_scancel=incl_sc, gate=_start_gate)
                elif cmd == "stop":
                    fails += cmd_stop(cfg, data, pj, jb, a.yes)
                elif cmd == "retry":
                    fails += cmd_retry(cfg, data, pj, jb, a.force,
                                       incl_scancel=incl_sc)
                elif cmd == "rerun":
                    fails += cmd_rerun(cfg, data, pj, jb, a.yes, a.force,
                                       from_skill=a.from_skill)
                elif cmd == "fetch":
                    fails += cmd_fetch(cfg, data, pj, all_files=a.all_files)
        if fails:
            sys.exit(1)
