# -*- coding: utf-8 -*-
"""agentgate —— LLM 动作网关 + 审计（v1.0 P0-1）。

**要解决的问题**：tf 分不清"人敲的"和"agent 敲的"。LLM 会话里一次手滑的
`autozt rerun` 就会 rm -rf 掉算完的步骤目录，事后连"谁、什么时候、下的什么手"
都查不到（材料目录里的 tf.log 只记成功动作，不记被拒的调用、没有 actor）。

**做法**（两层，都不改既有命令的行为）：

1. **网关 `autozt act <真实命令>`** —— agent 的唯一入口。命令按风险分三档：
   - `read`（list/summary/status/dir/skills/skill/schema/history/prove/probe/
     config/help/diagnose/session，以及任何带 `--dry-run` 的调用）：直接放行；
   - `mutate`（start/retry/fetch/init/adopt/level/hpc/auto/conf --set/correct/
     monitor…）：放行（本来就是 agent 该干的事）；
   - `destructive`（stop/rerun/clean/migrate-subdir/push/`correct -y`，以及任何带
     `-f`/`--force`/`-y`/`--yes`/`--purge-config` 的调用）：**必须用户同意**。
     被拦下时签发请求号 + 打印审批卡片（命令、后果、任务图里标出目标步骤）；
     chat 模式（默认）agent 把卡片给用户看、用户在对话里同意后，agent 执行
     `autozt approve --request <号> --reply "<用户原话>"` 就地批准并执行；
     tty 模式（approval_mode: tty）只认交互终端里的人工批准。批准按"命令签名"记账，
     默认 15 分钟、一次用完即销。

2. **审计 `.tf_agent_log.jsonl`**（落在配置目录）：每次 `autozt act` 调用都追加一条
   {ts, actor, cmd, argv, risk, decision, exit_code, dur, cwd, approved_by, …}
   ——放行、拒绝、失败都记。另外，agent 会话（环境变量 `AUTOZT_ACTOR` 已设）**直接**
   敲 tf 也会记一条（decision=direct），堵住"绕过网关就没人知道"。

**边界与开关**：
- 网关默认开启（`agent_gate: auto`，可在 tf.yaml 写或用 `AUTOZT_AGENT_GATE`）：
  以前只有设了 `AUTOZT_ACTOR` 才算 agent 会话，AI 忘了设就直接放行（2026-09 实测）。
  现在 agent 会话按下面顺序识别：`AUTOZT_ACTOR` → 常见 AI 代理留在子进程环境里的
  标记（CLAUDECODE / GEMINI_CLI / CODEX_SANDBOX …，可用 tf.yaml 的
  `agent_env_markers` 追加，如 DSH 自己的变量）→ 都没有时，**非交互终端里的破坏性
  命令**也按 agent 处理（人在终端里敲不受影响；approve 必须在真 TTY 里做）；
- `agent_gate: off`（或 `AUTOZT_AGENT_GATE=off`）→ 回到旧行为：只认 `AUTOZT_ACTOR`；
- 识别为 agent → 直接调用也审计；破坏性动作是否要令牌由 `AUTOZT_AGENT_STRICT`
  决定（缺省严格；`AUTOZT_AGENT_STRICT=0` 可关，给用户自己的定时脚本留后门）；
- 只写配置目录下的两个文件，不碰材料目录、不碰超算、不发任何网络请求。

铁律对齐：stop/rerun/clean/-f/-y 本来就要人工同意（~/.dsh/AGENTS.md 第 2 节、
仓库 AGENTS.md 铁律 2）；这里把"口头同意"变成可审计的一次性令牌。
"""

import os
import sys
import time
import json
import hashlib
import datetime
import subprocess
import threading

_AUDIT_WRITE_LOCK = threading.Lock()

# 审计与批准落盘位置（都在**配置目录**里，和 history.jsonl / .tf_hung.json 同级）
AGENT_LOG_NAME = ".tf_agent_log.jsonl"
AGENT_APPROVAL_NAME = ".tf_approvals.json"
AGENT_ACTOR_ENV = "AUTOZT_ACTOR"
AGENT_STRICT_ENV = "AUTOZT_AGENT_STRICT"
AGENT_GATEWAY_ENV = "AUTOZT_AGENT_GATEWAY"     # 网关子进程用它避免重复记账
AGENT_TTL_ENV = "AUTOZT_APPROVE_TTL"
AGENT_TTL_DEFAULT = 900                    # 批准有效期（秒）
AGENT_GATE_ENV = "AUTOZT_AGENT_GATE"       # auto（默认）/ off
# 批准方式：chat（默认）= agent 在对话里问、用户同意后 agent 凭请求号执行
#           （MCP 客户端支持 elicitation 时直接弹确认框）；tty = 只认交互终端里的人工批准
AGENT_APPROVAL_MODE_ENV = "AUTOZT_APPROVAL_MODE"
AGENT_CARD_GRAPH_TIMEOUT = 90               # 审批卡片取任务图的采集上限（秒）
# 常见 AI 编码代理在它执行的 shell 子进程里设置的环境变量；存在即视为 agent 会话。
# 其它代理（如 DSH）在 tf.yaml 里用 agent_env_markers: [变量名, ...] 追加。
AGENT_ENV_MARKERS = ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT", "GEMINI_CLI",
                     "CODEX_SANDBOX", "CODEX_SANDBOX_NETWORK_DISABLED")

# 三档风险：read / mutate / destructive
AGENT_READ_CMDS = {
    "list", "summary", "status", "json", "dir", "skills", "skill", "schema",
    "history", "prove", "probe", "config", "help", "diagnose", "session",
    "progress", "doctor", "check", "graph",
}
AGENT_MUTATE_CMDS = {
    "start", "retry", "fetch", "advance", "init", "adopt", "level", "hpc", "auto",
    "conf", "correct", "monitor", "watch", "restart", "register",
}
AGENT_DESTRUCTIVE_CMDS = {"stop", "rerun", "clean", "migrate-subdir", "push"}
AGENT_DESTRUCTIVE_FLAGS = {"-f", "--force", "-y", "--yes", "--purge-config"}

# 取值型选项：扫"命令词"时要跳过它们后面的值（autozt act -p Si stop → stop 才是命令）
AGENT_VALUE_FLAGS = {
    "-p", "-j", "-job", "-tt", "-c", "--config", "-x", "--exclude",
    "-status", "--status", "-i", "--interval", "-n", "--since", "--limit",
    "--offset", "--host", "-u", "--user", "--set", "--out",
    # 以前漏了 --project：`act --project P advance` 会把 P 当命令词 → 按未登记命令拒绝
    "-proj", "--project", "--expect-state",
    # register 的取值选项（值是路径/集群名，不能被当成命令词）
    "--dataset", "--poscar", "--root", "--cluster",
}


# ===== 基础工具 =====
def agent_actor():
    """当前动作的发起者：AUTOZT_ACTOR 优先；空串 = 不是 agent 会话。"""
    return (os.environ.get(AGENT_ACTOR_ENV) or "").strip()


def agent_strict():
    """破坏性动作要不要令牌。设了 AUTOZT_ACTOR 就默认严格；显式 0/false/off/关 可关。"""
    v = os.environ.get(AGENT_STRICT_ENV)
    if v is None or str(v).strip() == "":
        return True
    return str(v).strip().lower() not in ("0", "false", "no", "off", "关", "否")


def agent_gate_policy(cfg=None):
    """网关策略：auto（默认，AI 调用默认走网关）/ off（只认 AUTOZT_ACTOR）。"""
    v = (os.environ.get(AGENT_GATE_ENV) or (cfg or {}).get("agent_gate") or "auto")
    v = str(v).strip().lower()
    if v in ("0", "false", "no", "off", "关", "否"):
        return "off"
    return "auto"


def agent_detect(cfg=None):
    """识别 agent 会话：返回 (actor, 来源)；不是 agent 会话返回 ("", None)。"""
    actor = agent_actor()
    if actor:
        return actor, "env:" + AGENT_ACTOR_ENV
    if agent_gate_policy(cfg) == "off":
        return "", None
    extra = (cfg or {}).get("agent_env_markers") or []
    if isinstance(extra, str):
        extra = [x.strip() for x in extra.split(",")]
    for name in list(AGENT_ENV_MARKERS) + [str(x) for x in extra if str(x).strip()]:
        if (os.environ.get(name) or "").strip():
            return "agent:%s" % name.lower(), "env:" + name
    return "", None


def _stdin_isatty():
    try:
        return bool(sys.stdin and sys.stdin.isatty())
    except (AttributeError, ValueError):
        return False


def agent_ttl():
    try:
        return max(1, int(os.environ.get(AGENT_TTL_ENV) or AGENT_TTL_DEFAULT))
    except (TypeError, ValueError):
        return AGENT_TTL_DEFAULT


def agent_cfg_dir(cfg=None):
    cfg = cfg or {}
    return (cfg.get("_config_dir")
            or (os.path.dirname(os.path.abspath(cfg["_config_path"]))
                if cfg.get("_config_path") else os.getcwd()))


def agent_log_path(cfg=None):
    return os.path.join(agent_cfg_dir(cfg), AGENT_LOG_NAME)


def agent_approval_path(cfg=None):
    return os.path.join(agent_cfg_dir(cfg), AGENT_APPROVAL_NAME)


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _now_epoch():
    return time.time()


# ===== 命令切分与风险判定 =====
def agent_split(raw_argv):
    """把 sys.argv[1:] 切成 (verb, inner, outer)。

    verb  = "act" / "approve" 出现的位置（第一个裸词）；
    inner = verb 之后的**原样** token（真实命令，含它的 -p/-j/-y…）；
    outer = verb 之前允许的全局选项（-c/--config、--host、-u/--user）。

    verb 之前若混进了 -p/-j/-tt/-f/-y 之类，直接拒绝——那会造成"批准的命令"
    与"执行的命令"不是同一条（`tf -p Si act stop` 里 -p 会被丢掉）。
    """
    toks = [str(x) for x in (raw_argv or [])]
    idx = None
    for i, t in enumerate(toks):
        if t in ("act", "approve"):
            idx = i
            break
    if idx is None:
        return None, list(toks), []
    outer, i = [], 0
    keep_pair = {"-c", "--config", "--host", "-u", "--user"}
    while i < idx:
        t = toks[i]
        if t in keep_pair and i + 1 < idx:
            outer += [t, toks[i + 1]]
            i += 2
            continue
        if not t.startswith("-"):
            return ("act-error", list(toks[idx + 1:]), outer)
        i += 1
    return toks[idx], list(toks[idx + 1:]), outer


_OUTER_PAIRS = ("-c", "--config", "--host", "-u", "--user")


def _split_outer(tokens):
    """直接调用的 argv → (命令 token, 全局选项)。-c/--host/-u 与 act 网关同口径放进
    全局选项、不进签名：直接调用和 act 调用的同一条命令签名一致。"""
    inner, outer, i = [], [], 0
    while i < len(tokens):
        t = tokens[i]
        key = t.split("=", 1)[0]
        if key in _OUTER_PAIRS:
            if "=" in t:
                outer.append(t)
                i += 1
                continue
            if i + 1 < len(tokens):
                outer += [t, tokens[i + 1]]
                i += 2
                continue
        inner.append(t)
        i += 1
    return inner, outer


def agent_command(tokens):
    """从 token 流里挑出真正的命令词（跳过取值型选项的值）。"""
    skip = False
    for t in tokens or []:
        if skip:
            skip = False
            continue
        if t.startswith("-"):
            base = t.split("=", 1)[0]
            if "=" not in t and base in AGENT_VALUE_FLAGS:
                skip = True
            continue
        return t
    return None


def agent_classify(cmd, argv=()):
    """返回 (risk, why)：risk ∈ read / mutate / destructive，why 是中文理由清单。"""
    argv = [str(x) for x in (argv or [])]
    flags = set()
    for x in argv:
        if x.startswith("-"):
            flags.add(x.split("=", 1)[0])
    if "--dry-run" in flags:
        return "read", ["--dry-run 排练：只打印将影响的对象，无副作用"]
    if cmd == "migrate-subdir":
        return "destructive", ["会迁移/重排技能目录和配置，可能改变项目路径"]
    if cmd == "push":
        return "destructive", ["会向远端 Git 仓库写入提交，属于外部不可逆副作用"]
    if cmd in AGENT_DESTRUCTIVE_CMDS:
        return "destructive", ["%s 会取消作业 / 删除已算产物" % cmd]
    hit = sorted(flags & AGENT_DESTRUCTIVE_FLAGS)
    if hit:
        return "destructive", ["带强制/免确认开关 %s" % " ".join(hit)]
    if cmd == "conf" and "--set" not in flags:
        return "read", ["conf 不带 --set：只读展示"]
    if cmd in AGENT_MUTATE_CMDS:
        return "mutate", ["%s 会生成输入/推进状态" % cmd]
    if cmd in AGENT_READ_CMDS:
        return "read", ["只读命令"]
    # Fail closed: a newly added CLI command must be classified explicitly
    # before an agent can run it.  This keeps future destructive commands from
    # silently inheriting the mutate/allow policy.
    return "destructive", ["未登记命令，默认按破坏性动作处理，需人工批准"]


def agent_signature(cmd, argv):
    """命令签名：对"命令词 + 参数（忽略顺序与 -y/--yes）"取 sha256 前 12 位。

    人工 `autozt approve` 与 agent `autozt act` 各算一次，参数顺序不同也能对上；
    -y/--yes 只是"我知道自己在干什么"，不参与签名。
    """
    toks = [str(x) for x in (argv or []) if str(x) not in ("-y", "--yes")]
    payload = json.dumps({"cmd": cmd, "args": sorted(toks)},
                         ensure_ascii=False, sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def agent_targets(argv):
    """从参数里抠出材料/步骤/技能，便于审计与按材料过滤。"""
    out = {"mat": None, "step": None, "tt": None, "host": None}
    key = {"-p": "mat", "--proj": "mat", "-tt": "tt",
           "-j": "step", "-job": "step", "--host": "host"}
    toks = [str(x) for x in (argv or [])]
    i = 0
    while i < len(toks):
        t = toks[i]
        if t in key and i + 1 < len(toks):
            out[key[t]] = toks[i + 1]
            i += 2
            continue
        i += 1
    return out


# ===== 审计日志 =====
def agent_audit(cfg, actor, cmd, argv, risk, decision, why=None,
                exit_code=None, dur=None, approved_by=None, sig=None,
                gateway="act", ev="call", note=None, extra=None):
    """追加一条审计记录。整段包 try/except：审计失败绝不影响命令本身。"""
    rec = {
        "ts": _now(),
        "ev": ev,
        "actor": actor or "",
        "cmd": cmd or "",
        "argv": list(argv or []),
        "risk": risk or "",
        "gateway": gateway,
        "decision": decision,
        "exit_code": exit_code,
        "dur": round(dur, 2) if isinstance(dur, (int, float)) else None,
        "cwd": os.getcwd(),
        "approved_by": approved_by,
        "sig": sig,
        "why": list(why or []),
    }
    rec.update({k: v for k, v in agent_targets(argv).items() if v})
    if note:
        rec["note"] = str(note)[:200]
    for k, v in (extra or {}).items():
        if v not in (None, ""):
            rec[k] = v if not isinstance(v, str) else v[:500]
    try:
        path = agent_log_path(cfg)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        line = (json.dumps(rec, ensure_ascii=False) + "\n").encode("utf-8")
        # One append write prevents interleaving records between threads in
        # this process; the reader remains tolerant of pre-existing torn rows.
        with _AUDIT_WRITE_LOCK:
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
            try:
                os.write(fd, line)
            finally:
                os.close(fd)
    except Exception as _e:                     # noqa: BLE001
        print("警告：审计日志写不进去（忽略）：%s" % _e, file=sys.stderr)
    return rec


def agent_log_load(cfg, proj=None, since=None, risk=None, limit=None):
    """读审计日志并按材料/时间/风险过滤。返回 (记录列表, 文件总行数)。"""
    path = agent_log_path(cfg)
    out, total = [], 0
    if not os.path.isfile(path):
        return out, 0
    wants = {x.strip() for x in str(proj or "").split(",") if x.strip()}
    since = str(since or "").strip().replace(" ", "T")
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    e = json.loads(line)
                except ValueError:
                    continue
                total += 1
                if wants:
                    mat = str(e.get("mat") or "")
                    if not (mat in wants or os.path.basename(mat) in wants):
                        continue
                if since and str(e.get("ts") or "") < since:
                    continue
                if risk and e.get("risk") != risk:
                    continue
                out.append(e)
    except OSError:
        return [], 0
    if limit:
        out = out[-int(limit):]
    return out, total


# ===== 人工批准（一次性令牌）=====
def _approvals_read(path):
    try:
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}


def _approvals_write(path, d):
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(d, f, ensure_ascii=False, indent=2)
        os.replace(tmp, path)
        return True
    except OSError:
        return False


def agent_save_approval(cfg, sig, cmd, argv, by, ttl=None):
    """人工批准一条命令签名（一次性、带有效期）。"""
    path = agent_approval_path(cfg)
    d = _approvals_read(path)
    items = [x for x in (d.get("approvals") or [])
             if not x.get("used")
             and float(x.get("expires_epoch") or 0) > _now_epoch()]
    items.append({"sig": sig, "cmd": cmd, "argv": list(argv or []),
                  "by": by or "", "ts": _now(), "used": False,
                  "expires_epoch": _now_epoch() + float(ttl or agent_ttl())})
    d["approvals"] = items
    d["updated"] = _now()
    return _approvals_write(path, d)


def agent_take_approval(cfg, sig):
    """取用一条尚未用过的批准（用完置 used，防止重放）。返回 (ok, by)。"""
    path = agent_approval_path(cfg)
    d = _approvals_read(path)
    items = list(d.get("approvals") or [])
    now = _now_epoch()
    hit = None
    for x in items:
        if (not x.get("used") and x.get("sig") == sig
                and float(x.get("expires_epoch") or 0) > now):
            hit = x
            break
    if hit is None:
        return False, None
    hit["used"] = True
    hit["used_at"] = _now()
    d["approvals"] = items
    _approvals_write(path, d)
    return True, hit.get("by")


def agent_pending_approvals(cfg):
    d = _approvals_read(agent_approval_path(cfg))
    now = _now_epoch()
    return [x for x in (d.get("approvals") or [])
            if not x.get("used")
            and float(x.get("expires_epoch") or 0) > now]


# ===== 对话内审批：请求号 + 审批卡片（含任务图）=====
# 流程：agent 调破坏性命令 → 网关不执行，签发一个请求号并打印「审批卡片」（命令、后果、
# 任务图里标出目标步骤）→ agent 把卡片给用户看、问要不要做 → 用户同意 → agent 执行
#   autozt approve --request <请求号> --reply "<用户原话>"
# 就地批准并**立即执行**（不用用户自己去终端敲命令）。MCP 客户端支持 elicitation 时，
# 由 MCP 服务直接弹确认框（用户点「执行」才算）。请求号绑定命令签名：一次有效、
# 有有效期、改任何参数都要重新申请；同意的原话、通过哪条渠道批准都进审计日志。
# approval_mode: tty（tf.yaml 或 AUTOZT_APPROVAL_MODE）回到「只认交互终端」的严格模式。
_REQ_KEEP_S = 86400              # 请求记录保留一天（审计日志另有完整流水）
_NEGATIVE_REPLIES = {"不", "否", "no", "n", "不要", "不同意", "拒绝", "别", "算了",
                     "不行", "先别", "nope", "deny"}
_ACTION_WORD = {"stop": "取消", "rerun": "删除重算", "clean": "清理",
                "migrate-subdir": "迁移目录", "push": "推送", "start": "强制提交",
                "retry": "强制重生成"}
_CONSEQ = {
    "stop": "scancel 目标步骤排队/在跑的作业，并打 SCANCEL 标记（auto 不会自动重跑）；"
            "已有文件保留，之后 start 可续交。",
    "rerun": "先 scancel，再删除该步骤目录（已算产物全部删除、不可恢复），然后重新生成并提交。",
    "clean": "删除生成物回到 PREP（关联作业一并 scancel；材料级保留 POSCAR）。",
    "migrate-subdir": "迁移/重排技能目录和配置，可能改变项目路径。",
    "push": "向远端 Git 仓库写入提交（外部不可逆）。",
}
_CONFIRM_CMDS = ("stop", "rerun", "clean", "migrate-subdir")   # 自带 [y/N] 二次确认的命令


def agent_approval_mode(cfg=None):
    """chat（默认）：对话里征得同意后 agent 凭请求号执行；tty：只认交互终端。"""
    v = (os.environ.get(AGENT_APPROVAL_MODE_ENV) or (cfg or {}).get("approval_mode")
         or "chat")
    v = str(v).strip().lower()
    return "tty" if v in ("tty", "terminal", "strict", "终端", "严格") else "chat"


def agent_consequences(cmd, argv):
    """一句话后果（给人看的，不是给 agent 看的）。"""
    flags = {str(x).split("=", 1)[0] for x in (argv or []) if str(x).startswith("-")}
    out = []
    whole = not agent_targets(argv).get("step")
    if cmd == "rerun" and whole:
        out.append("不带 -j = 整材料重来：scancel 全部作业，清空整个工作目录（只留 POSCAR，"
                   "本地 result 也删），从第一步重新生成提交。不可恢复。")
    elif cmd == "clean" and whole:
        out.append("不带 -j = 整材料清理：删除全部生成物回到 PREP（关联作业一并 scancel，保留 POSCAR）。")
    elif cmd in _CONSEQ:
        out.append(_CONSEQ[cmd])
    if flags & {"-f", "--force"}:
        out.append("带 -f：越过安全闸门（max_jobs 并发上限 / 依赖未完成 / 本机环境检查等）强行执行。")
    if "--purge-config" in flags:
        out.append("--purge-config：连 project_setting 一起删（之后要重新 init）。")
    if not out:
        out.append("未登记的命令，按破坏性动作处理。")
    return " ".join(out)


def _split_list(v):
    return [x.strip() for x in str(v or "").split(",") if x.strip()]


def _node_hit(node, tok):
    """-j 写法（label / name / 末段名 / 序号）是否指这个节点。"""
    tok = str(tok).strip()
    lb, nm = str(node.get("label") or ""), str(node.get("name") or "")
    if tok in (lb, nm) or nm.rsplit("/", 1)[-1] == tok:
        return True
    if tok.replace(".", "", 1).isdigit():
        return lb.startswith("S%s_" % tok) or nm.rsplit("/", 1)[-1].startswith("step%s_" % tok)
    return False


def agent_card_marks(graph, cmd, argv):
    """在任务图上标出本次动作的目标步骤（◀）和受影响的下游（↳）。"""
    from autozt.taskgraph import descendants
    word = _ACTION_WORD.get(cmd, "执行")
    steps = _split_list(agent_targets(argv).get("step"))
    nodes = graph.get("nodes") or []
    if steps:
        targets = [n["label"] for n in nodes if any(_node_hit(n, j) for j in steps)]
    elif cmd == "stop":
        targets = [n["label"] for n in nodes if n.get("kind") in ("R", "PD")]
    elif cmd in ("clean", "rerun"):         # 材料级：整个工作目录清空，所有步骤从头来
        targets = [n["label"] for n in nodes]
    elif cmd in ("start", "retry"):
        targets = [n["label"] for n in nodes
                   if n.get("kind") in ("TODO", "PREP", "SCANCEL", "FAIL", "WAIT")][:1]
    else:
        targets = []
    marks = {lb: "◀ 本次：%s" % word for lb in targets}
    if cmd in ("stop", "rerun", "clean"):
        down = "↳ 上游被取消，不会推进" if cmd == "stop" else "↳ 上游产物将被删除，需重算"
        kinds = {n["label"]: n.get("kind") for n in nodes}
        for lb in descendants(graph, targets):
            # stop 只影响还没算完、也没被取消的下游；rerun/clean 删了上游产物，下游全要重算
            if cmd != "stop" or kinds.get(lb) not in ("OK", "SCANCEL"):
                marks.setdefault(lb, down)
    return marks


def _card_graph_text(cfg, cmd, argv, outer=None):
    """取目标材料的任务图（只读子进程，带超时）并渲染成带标注的文本。

    返回 (文本或 None, [带标注的图 dict], 失败原因或 None)。"""
    from autozt import _PKG_ROOT
    from autozt.taskgraph import fetch_graphs, render_text, with_renderings
    tg = agent_targets(argv)
    mats = _split_list(tg.get("mat"))
    if not mats and not tg.get("step"):
        return None, [], None              # 全局动作：没有具体对象可画
    if (os.environ.get("AUTOZT_APPROVAL_GRAPH") or "").strip().lower() in (
            "0", "false", "no", "off"):
        return None, [], None
    outer = list(outer or [])
    if not any(x in ("-c", "--config") or x.startswith("--config=") for x in outer):
        if not (cfg or {}).get("_config_path"):
            # 不知道是哪份配置就不去采（单元测试/嵌入调用），免得碰到别的项目
            return None, [], None
        outer = ["-c", cfg["_config_path"]] + outer
    prog = os.path.join(_PKG_ROOT, "bin", "autozt")
    base = ([sys.executable, prog] if os.path.isfile(prog)
            else [sys.executable, "-m", "autozt"]) + outer
    try:
        wait = int((cfg or {}).get("approval_graph_timeout") or AGENT_CARD_GRAPH_TIMEOUT)
    except (TypeError, ValueError):
        wait = AGENT_CARD_GRAPH_TIMEOUT
    graphs, err = fetch_graphs(base, tt=tg.get("tt"), materials=mats[:3] or None,
                               timeout=wait)
    if not graphs:
        return None, [], err
    texts, out = [], []
    for g in graphs[:3]:
        marks = agent_card_marks(g, cmd, argv) if g.get("view") == "material" else {}
        if g.get("view") == "skill":
            steps = _split_list(tg.get("step"))
            marks = {n["label"]: "◀ 本次：%s（全部材料）" % _ACTION_WORD.get(cmd, "执行")
                     for n in g.get("nodes") or [] if any(_node_hit(n, j) for j in steps)}
        texts.append(render_text(g, marks=marks, ascii_mode=False, legend=False,
                                 indent="    "))
        out.append(with_renderings(g, marks=marks))
    if len(mats) > 3:
        texts.append("    …… 另有 %d 个材料（同一命令一起执行）" % (len(mats) - 3))
    return "\n\n".join(texts), out, None


def _cmd_text(argv):
    return "autozt " + " ".join(_shell_quote(x) for x in argv)


def _shell_quote(x):
    x = str(x)
    if x and all(c.isalnum() or c in "-_./,=:@%+" for c in x):
        return x
    return "'" + x.replace("'", "'\\''") + "'"


def agent_request_create(cfg, sig, cmd, argv, why, actor, gateway, card="",
                         graphs=None, outer=None):
    """签发（或复用同签名未过期的）审批请求号。返回请求 dict；写不进文件返回 None。"""
    import secrets
    path = agent_approval_path(cfg)
    d = _approvals_read(path)
    now = _now_epoch()
    reqs = [r for r in (d.get("requests") or [])
            if float(r.get("created_epoch") or 0) > now - _REQ_KEEP_S]
    hit = next((r for r in reqs if r.get("status") == "pending" and r.get("sig") == sig
                and float(r.get("expires_epoch") or 0) > now), None)
    if hit is None:
        used = {r.get("id") for r in reqs}
        rid = "r" + secrets.token_hex(3)
        while rid in used:
            rid = "r" + secrets.token_hex(3)
        hit = {"id": rid, "sig": sig, "cmd": cmd, "argv": list(argv or []),
               "why": list(why or []), "actor": actor or "", "gateway": gateway,
               "outer": list(outer or []),
               "created": _now(), "created_epoch": now,
               "expires_epoch": now + float(agent_ttl()), "status": "pending"}
        reqs.append(hit)
    hit["card"] = card or hit.get("card") or ""
    if graphs is not None:
        hit["graphs"] = [{k: g.get(k) for k in ("view", "skill", "id", "text", "mermaid",
                                                 "marks") if g.get(k) is not None}
                         for g in graphs]
    d["requests"] = reqs
    d["updated"] = _now()
    return hit if _approvals_write(path, d) else None


def agent_request_get(cfg, rid):
    d = _approvals_read(agent_approval_path(cfg))
    return next((r for r in (d.get("requests") or []) if r.get("id") == rid), None)


def agent_request_update(cfg, rid, **fields):
    path = agent_approval_path(cfg)
    d = _approvals_read(path)
    for r in d.get("requests") or []:
        if r.get("id") == rid:
            r.update(fields)
            d["updated"] = _now()
            return _approvals_write(path, d)
    return False


def agent_requests_pending(cfg):
    now = _now_epoch()
    d = _approvals_read(agent_approval_path(cfg))
    return [r for r in (d.get("requests") or []) if r.get("status") == "pending"
            and float(r.get("expires_epoch") or 0) > now]


def agent_render_card(cmd, argv, why, req, graph_text, graph_err=None, mode="chat"):
    """审批卡片：给用户看的「要干什么、会怎样、在流程图的哪儿」。"""
    left = int(max(0, float(req.get("expires_epoch") or 0) - _now_epoch())) if req else 0
    L = ["━━ 需要你的同意（破坏性动作，AGENTS.md 铁律 2）━━",
         "命令   %s" % _cmd_text(argv),
         "后果   %s" % agent_consequences(cmd, argv)]
    if graph_text:
        L.append("任务图（◀ = 本次动作的目标，↳ = 受影响的下游）")
        L.append(graph_text)
    elif graph_err:
        L.append("任务图   暂时取不到（%s）；可先看：autozt -p <材料> graph" % graph_err)
    if req:
        L.append("请求号 %s（%d 秒内有效，只能用一次；改任何参数都要重新申请）"
                 % (req["id"], left))
        if mode == "chat":
            L.append("同意 → 回复「同意」即可，agent 会执行：")
            L.append("        autozt approve --request %s --reply \"<你的原话>\"" % req["id"])
            L.append("不同意 → 回复「不」（请求作废：autozt approve --deny %s）" % req["id"])
        else:
            L.append("本机是严格模式（approval_mode: tty）：请你在交互终端运行")
            L.append("        autozt approve --request %s" % req["id"])
    else:
        L.append("（请求号写不进配置目录，只能用旧方式：在交互终端 autozt approve %s）"
                 % " ".join(str(x) for x in argv))
    return "\n".join(L)


def agent_request_approval(cfg, cmd, inner, why, actor, gateway, outer=None):
    """破坏性命令被拦下时：签发请求号 + 打印审批卡片（含任务图）+ 记审计。返回 3。"""
    sig = agent_signature(cmd, inner)
    graph_text, graphs, gerr = _card_graph_text(cfg, cmd, inner, outer=outer)
    mode = agent_approval_mode(cfg)
    req = agent_request_create(cfg, sig, cmd, inner, why, actor, gateway, graphs=graphs,
                               outer=outer)
    card = agent_render_card(cmd, inner, why, req, graph_text, gerr, mode=mode)
    if req:
        agent_request_update(cfg, req["id"], card=card)
    agent_audit(cfg, actor, cmd, inner, "destructive", "deny-need-approval", why=why,
                exit_code=3, sig=sig, gateway=gateway,
                extra={"request": (req or {}).get("id"), "approval_mode": mode})
    print("✗ 拒绝执行（先征得用户同意）：tf %s 属破坏性动作。" % " ".join(inner))
    print(card)
    print("  命令签名：%s  审计日志：%s" % (sig, agent_log_path(cfg)))
    if os.environ.get("AUTOZT_AGENT_EVENTS") == "1":
        # MCP / 执行器用：结构化的审批请求（一行 JSON）
        ev = {"request_id": (req or {}).get("id"), "sig": sig, "cmd": cmd,
              "command": _cmd_text(inner), "approval_mode": mode,
              "consequences": agent_consequences(cmd, inner), "why": list(why or []),
              "expires_in_s": int(max(0, float((req or {}).get("expires_epoch") or 0)
                                      - _now_epoch())),
              "card": card, "graphs": graphs}
        print("[agent-approval] " + json.dumps(ev, ensure_ascii=False))
    sys.stdout.flush()
    return 3


# ===== 渲染 =====
AGENT_RISK_CN = {"read": "只读", "mutate": "推进", "destructive": "破坏性"}
AGENT_POLICY_ROWS = [
    ("read", "list summary status json dir skills skill schema history prove "
             "probe config help diagnose session progress doctor check graph / 任何 --dry-run", "放行"),
    ("mutate", "start retry fetch advance init adopt level hpc auto conf --set "
               "correct monitor restart", "放行（记账）"),
    ("destructive", "stop rerun clean migrate-subdir push / 未登记命令 / 任何 -f -y --yes --purge-config",
     "需用户同意（请求号 + 审批卡片）"),
]


def _wrap_words(text, width=78):
    """按空格把长串折行（纯展示用）。"""
    out, cur = [], ""
    for w in str(text).split():
        if cur and len(cur) + 1 + len(w) > width:
            out.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    if cur:
        out.append(cur)
    return out or [""]


def render_agent_policy():
    L = ["autozt act —— agent 动作网关（v1.0 P0-1）", ""]
    for risk, cmds, policy in AGENT_POLICY_ROWS:
        L.append("【%s】 %s" % (AGENT_RISK_CN.get(risk, risk), policy))
        L.append("    命令：%s" % _wrap_words(cmds, 74)[0])
        for extra in _wrap_words(cmds, 74)[1:]:
            L.append("          %s" % extra)
    L.append("")
    L.append("用法：")
    L.append("  autozt act <和 tf 一模一样的命令>     # agent 的唯一入口（自动记账）")
    L.append("  autozt approve --list                 # 看待你同意的请求（审批卡片 + 任务图）")
    L.append("  autozt approve --request <号> --reply \"同意\"   # 批准并立即执行")
    L.append("  autozt approve --deny <号>            # 拒绝/作废")
    L.append("  autozt act log [-n 40] [--json]       # 看审计流水")
    L.append("  autozt act policy                     # 看这张表")
    L.append("")
    L.append("说明：破坏性命令被拦下时签发一个请求号并打印审批卡片（命令、后果、任务图里标出目标步骤）。")
    L.append("      chat 模式（默认）：agent 把卡片给用户看，用户在对话里同意后 agent 执行")
    L.append("      approve --request（必须带 --reply 用户原话，进审计）；MCP 客户端支持 elicitation")
    L.append("      时直接弹确认框。tty 模式（tf.yaml approval_mode: tty）：只认交互终端里的人工批准。")
    L.append("      请求号绑定命令签名（参数顺序无关），默认 %d 秒内一次有效。" % AGENT_TTL_DEFAULT)
    L.append("      默认网关开启（agent_gate: auto）：认 AUTOZT_ACTOR、常见 AI 代理环境标记，")
    L.append("      以及非交互终端里的破坏性命令；tf.yaml agent_gate: off / AUTOZT_AGENT_GATE=off 关。")
    L.append("      环境变量：AUTOZT_ACTOR=名字（谁在操作）、AUTOZT_AGENT_STRICT=0（关令牌）、")
    L.append("                AUTOZT_APPROVE_TTL=秒（批准有效期）。")
    return "\n".join(L)


def agent_render_log(cfg, n=40, json_out=False, proj=None, since=None):
    path = agent_log_path(cfg)
    evs, total = agent_log_load(cfg, proj=proj, since=since)
    if json_out:
        print(json.dumps({"path": path, "total": total, "count": len(evs),
                          "events": evs[-int(n or 40):]},
                         ensure_ascii=False, indent=2))
        return 0
    if total == 0:
        print("还没有审计记录（%s 不存在或为空）。" % path)
        print("agent 侧：autozt act <命令> 自动记一条；设了 AUTOZT_ACTOR 的直接调用也记。")
        return 0
    show = evs[-int(n or 40):]
    print("agent 审计  %s" % path)
    print("共 %d 条（文件累计 %d 条），显示最近 %d 条：" % (len(evs), total, len(show)))
    print("")
    print("%-19s %-10s %-9s %-11s %-24s %-6s %s"
          % ("时间", "actor", "风险", "决定", "命令", "退出码", "批准人"))
    for e in show:
        argv = " ".join(str(x) for x in (e.get("argv") or []))
        line = ("%s %s" % (e.get("cmd") or "?", argv)).strip()
        print("%-19s %-10s %-9s %-11s %-24s %-6s %s"
              % (str(e.get("ts") or "").replace("T", " ")[:19],
                 str(e.get("actor") or "")[:10],
                 AGENT_RISK_CN.get(e.get("risk"), str(e.get("risk") or ""))[:9],
                 str(e.get("decision") or "")[:11],
                 line[:24],
                 "-" if e.get("exit_code") is None else e.get("exit_code"),
                 str(e.get("approved_by") or "-")))
    pend = agent_pending_approvals(cfg)
    if pend:
        print("")
        print("尚有 %d 条未使用的批准（签名 %s）。"
              % (len(pend), ", ".join(str(x.get("sig")) for x in pend)))
    if len(evs) > len(show):
        print("看更多：-n %d（或 --json 拿结构化数据）" % (len(evs) * 2))
    return 0


def agent_deny_message(cfg, cmd, inner, sig, why, prog=None):
    prog = prog or "bin/autozt"
    L = ["✗ 拒绝执行：tf %s %s 属**破坏性动作**，agent 不能自行决定（AGENTS.md 铁律 2）。"
         % (cmd, " ".join(str(x) for x in inner))]
    if why:
        L.append("  理由：%s" % "；".join(why))
    L.append("")
    L.append("  请**人工**在交互终端里批准（批准后一次有效，默认 %d 秒）："
             % agent_ttl())
    L.append("      python3 %s approve %s" % (prog, " ".join(str(x) for x in inner)))
    L.append("")
    L.append("  命令签名：%s（参数顺序无关；改了任何一个参数都要重新批准）" % sig)
    L.append("  审计日志：%s" % agent_log_path(cfg))
    return "\n".join(L)


# ===== 命令入口 =====
def cmd_act(cfg, raw_argv):
    """autozt act <真实命令> —— agent 的唯一入口：判定风险 + 记账 + 转发。

    返回子进程退出码（放行）或 3（拒绝）。
    """
    from autozt import _PKG_ROOT
    verb, inner, outer = agent_split(raw_argv)
    if verb == "act-error":
        print("错误：act 之前的选项里混进了 -p/-j/-tt/-f/-y 这类会影响目标的参数。\n"
              "      请把它们写在 act 之后（否则批准的命令和执行的命令不是同一条）：\n"
              "        autozt act -p <材料> [-j <步骤>] <命令>")
        return 2
    inner, _outer2 = _split_outer(inner)    # act 之后写的 -c/--host 同样不进签名
    outer = list(outer) + _outer2
    if not inner:
        print(render_agent_policy())
        return 0
    sub = agent_command(inner)
    if sub == "policy":
        print(render_agent_policy())
        return 0
    if sub == "log":
        return agent_render_log(cfg, n=40, json_out=False)
    if sub in ("act", "approve"):
        print("错误：act/approve 不能嵌套（网关只此一层）。")
        return 2
    if sub is None:
        print(render_agent_policy())
        return 0
    risk, why = agent_classify(sub, inner)
    actor = agent_detect(cfg)[0] or os.environ.get("USER") or "?"
    sig = agent_signature(sub, inner)
    approved_by = None
    prog = os.path.join(_PKG_ROOT, "bin", "autozt")
    if risk == "destructive":
        ok, approver = agent_take_approval(cfg, sig)
        if not ok:
            # 不执行：签发请求号 + 审批卡片（含任务图），等用户在对话里同意
            return agent_request_approval(cfg, sub, inner, why, actor, "act",
                                          outer=outer)
        approved_by = approver or "human"
    # pip 非 editable 安装时没有 bin/autozt：退回 python -m autozt（同一入口）
    child = ([sys.executable, prog] if os.path.isfile(prog)
             else [sys.executable, "-m", "autozt"]) + outer + inner
    env = dict(os.environ)
    env[AGENT_GATEWAY_ENV] = "act"        # 子进程不再重复记账/重复要令牌
    t0 = time.time()
    try:
        rc = subprocess.call(child, env=env)
    except OSError as _e:
        agent_audit(cfg, actor, sub, inner, risk, "error", why=why, exit_code=127,
                    sig=sig, approved_by=approved_by, note=str(_e), gateway="act")
        print("错误：无法执行 tf 子进程：%s" % _e)
        return 127
    dur = time.time() - t0
    agent_audit(cfg, actor, sub, inner, risk, "allow", why=why, exit_code=rc,
                dur=dur, sig=sig, approved_by=approved_by, gateway="act")
    return rc


_APPROVE_OPTS = {"--request": 1, "--reply": 1, "--deny": 1, "--via": 1,
                 "--list": 0, "--no-run": 0, "--json": 0}


def _approve_parse(inner):
    """把 approve 自己的选项（--request/--reply/--deny/--via/--list/--no-run/--json）
    从命令 token 里拆出来。返回 (opts, 其余 token)。"""
    opts, rest, i = {}, [], 0
    toks = [str(x) for x in (inner or [])]
    while i < len(toks):
        key, eq, val = toks[i].partition("=")
        if key in _APPROVE_OPTS:
            if _APPROVE_OPTS[key]:
                if eq:
                    opts[key], i = val, i + 1
                elif i + 1 < len(toks):
                    opts[key], i = toks[i + 1], i + 2
                else:
                    opts[key], i = "", i + 1
            else:
                opts[key], i = True, i + 1
            continue
        rest.append(toks[i])
        i += 1
    return opts, rest


def _approve_usage():
    return ("用法：\n"
            "  autozt approve --list                              # 看 agent 正在等你同意的请求（含任务图）\n"
            "  autozt approve --request <请求号> --reply \"同意\"    # 批准并立即执行（agent 在对话里征得同意后用）\n"
            "  autozt approve --deny <请求号>                     # 拒绝/作废\n"
            "  autozt approve -p <材料> [-j <步骤>] <命令>         # 旧写法：在交互终端批准一条命令（agent 再用 act 执行）")


def _approve_list(cfg, json_out=False):
    reqs = agent_requests_pending(cfg)
    if json_out:
        print(json.dumps({"pending": [{k: r.get(k) for k in (
            "id", "cmd", "argv", "created", "expires_epoch", "actor", "card", "graphs")}
            for r in reqs], "approval_mode": agent_approval_mode(cfg)},
            ensure_ascii=False))
        return 0
    if not reqs:
        print("没有待你同意的请求。\n")
        print(_approve_usage())
        return 0
    for r in reqs:
        print(r.get("card") or ("请求 %s：%s" % (r["id"], _cmd_text(r.get("argv")))))
        print()
    return 0


def _approve_deny(cfg, rid, reply=None):
    r = agent_request_get(cfg, rid)
    if not r:
        print("✗ 找不到请求号 %s（autozt approve --list 看待批准的）。" % rid)
        return 2
    if r.get("status") != "pending":
        print("请求 %s 已是 %s 状态，无需再拒绝。" % (rid, r.get("status")))
        return 0
    agent_request_update(cfg, rid, status="declined", declined=_now(), reply=reply or "")
    agent_audit(cfg, agent_detect(cfg)[0] or os.environ.get("USER") or "?", r.get("cmd"),
                r.get("argv"), "destructive", "declined-by-user", why=r.get("why"),
                sig=r.get("sig"), gateway="approve", ev="approve",
                extra={"request": rid, "reply": reply})
    print("✓ 已作废请求 %s（没有执行任何操作）。" % rid)
    return 0


def _approve_request(cfg, opts, outer):
    """凭请求号批准；默认批准后立即执行（--no-run 只发一次性令牌）。"""
    rid = str(opts.get("--request") or "").strip()
    r = agent_request_get(cfg, rid) if rid else None
    if not r:
        print("✗ 找不到请求号 %r（autozt approve --list 看待批准的）。" % rid)
        return 2
    if r.get("status") != "pending":
        print("✗ 请求 %s 已%s，不能再用；重新执行原命令会签发新的请求号。"
              % (rid, {"approved": "批准过", "used": "执行过", "declined": "被拒绝",
                       "expired": "过期"}.get(r.get("status"), r.get("status"))))
        return 3
    if float(r.get("expires_epoch") or 0) <= _now_epoch():
        agent_request_update(cfg, rid, status="expired")
        print("✗ 请求 %s 已过期；重新执行原命令会签发新的请求号（并重新展示任务图）。" % rid)
        return 3
    mode = agent_approval_mode(cfg)
    reply = str(opts.get("--reply") or "").strip()
    via = str(opts.get("--via") or "").strip()
    if _stdin_isatty() and not reply:
        # 人就在终端前：卡片里给 agent 的「同意 → …」提示不用再印
        print((r.get("card") or _cmd_text(r.get("argv"))).split("\n同意 →")[0])
        try:
            ans = input("确认批准并执行？[y/N] ")
        except EOFError:
            print("（没有输入，未批准。）")
            return 1
        if str(ans).strip().lower() not in ("y", "yes", "是", "好", "同意"):
            print("（已取消，未批准。）")
            return 1
        reply, via = str(ans).strip(), via or "tty"
    else:
        if mode != "chat":
            print("✗ 拒绝：本机是严格模式（approval_mode: tty），批准必须在**交互终端**里做：\n"
                  "      autozt approve --request %s" % rid)
            return 3
        if not reply:
            print("✗ 缺 --reply：把用户同意的原话原样传进来（写进审计日志），"
                  "如 --reply \"同意\"。没问过用户就不能批准。")
            return 2
        if reply.lower().rstrip("。.!！~ ") in _NEGATIVE_REPLIES:
            _approve_deny(cfg, rid, reply)
            print("✗ 用户的回复是拒绝（%r），没有执行。" % reply)
            return 3
        via = via or "chat"
    by = "user(%s)" % via
    if not agent_save_approval(cfg, r["sig"], r.get("cmd"), r.get("argv"), by):
        print("✗ 批准写不进 %s（检查配置目录权限）。" % agent_approval_path(cfg))
        return 1
    agent_request_update(cfg, rid, status="approved", approved=_now(), approved_by=by,
                         reply=reply, via=via)
    agent_audit(cfg, agent_detect(cfg)[0] or os.environ.get("USER") or "?", r.get("cmd"),
                r.get("argv"), "destructive",
                "approved-by-human" if via == "tty" else "approved-in-chat",
                why=r.get("why"), sig=r["sig"], approved_by=by, gateway="approve",
                ev="approve", extra={"request": rid, "reply": reply, "via": via})
    if opts.get("--no-run"):
        print("✓ 已批准（请求 %s，%d 秒内一次有效）。用同一条命令执行：autozt act %s"
              % (rid, agent_ttl(), " ".join(r.get("argv") or [])))
        return 0
    run_argv = list(r.get("argv") or [])
    if r.get("cmd") in _CONFIRM_CMDS and not ({"-y", "--yes"} & set(run_argv)):
        run_argv.append("-y")              # 用户已经在卡片上确认过，不再二次 [y/N]
    print("✓ 已批准（%s，原话：%s）。执行：%s" % (via, reply, _cmd_text(run_argv)))
    sys.stdout.flush()
    outer = list(r.get("outer") or outer or [])
    if "-c" not in outer and "--config" not in outer and (cfg or {}).get("_config_path"):
        outer = ["-c", cfg["_config_path"]] + outer
    rc = cmd_act(cfg, outer + ["act"] + run_argv)
    agent_request_update(cfg, rid, status="used", finished=_now(), exit_code=rc)
    return rc


def cmd_approve(cfg, raw_argv):
    """autozt approve —— 批准破坏性动作。

    --request <请求号>：批准 agent 被拦下的那条命令并立即执行（chat 模式下 agent 在对话里
    征得同意后可代为执行，必须带 --reply 用户原话；交互终端里会再问一次 y/N）；
    旧写法 approve <命令>：只在交互终端里发一次性令牌，agent 再用 act 执行。"""
    verb, inner, outer = agent_split(raw_argv)
    if verb == "act-error":
        print(_approve_usage())
        return 2
    opts, inner = _approve_parse(inner)
    if opts.get("--list") or (not inner and not opts):
        return _approve_list(cfg, json_out=bool(opts.get("--json")))
    if "--deny" in opts:
        return _approve_deny(cfg, str(opts["--deny"]).strip(), opts.get("--reply"))
    if "--request" in opts:
        return _approve_request(cfg, opts, outer)
    if not inner:
        print(_approve_usage())
        return 2
    sub = agent_command(inner)
    if sub is None:
        print(_approve_usage())
        return 2
    risk, why = agent_classify(sub, inner)
    if risk != "destructive":
        print("提示：tf %s 属「%s」档，本来就不需要批准（approve 只管破坏性动作）。"
              % (sub, AGENT_RISK_CN.get(risk, risk)))
        return 0
    if not sys.stdin.isatty():
        print("✗ 拒绝：按命令批准必须在**交互终端**里做（防止 agent 自己批准自己）。\n"
              "  agent 的正确做法：先执行原命令拿到请求号和审批卡片，给用户看、征得同意后\n"
              "      autozt approve --request <请求号> --reply \"<用户原话>\"")
        return 3
    sig = agent_signature(sub, inner)
    print("即将批准一次破坏性动作：")
    print("    命令    tf %s" % " ".join(str(x) for x in inner))
    print("    后果    %s" % agent_consequences(sub, inner))
    print("    签名    %s" % sig)
    print("    有效期  %d 秒（一次有效，用完即销）" % agent_ttl())
    try:
        ans = input("确认批准？[y/N] ")
    except EOFError:
        print("（没有输入，未批准。）")
        return 1
    if str(ans).strip().lower() not in ("y", "yes", "是", "好"):
        print("（已取消，未批准。）")
        return 1
    by = os.environ.get("USER") or "human"
    if not agent_save_approval(cfg, sig, sub, inner, by):
        print("✗ 批准写不进 %s（检查配置目录权限）。" % agent_approval_path(cfg))
        return 1
    agent_audit(cfg, by, sub, inner, risk, "approved-by-human", why=why,
                sig=sig, approved_by=by, gateway="approve", ev="approve")
    print("✓ 已批准（签名 %s，%d 秒内一次有效）。\n"
          "  agent 现在可以把**同一条命令**用 autozt act 重跑。" % (sig, agent_ttl()))
    return 0


def agent_direct_gate(cfg, cmd, raw_argv):
    """agent 会话（AUTOZT_ACTOR 已设）**直接**敲 tf 时的旁路钩子。

    返回 None（放行）或非零退出码（拒绝）。网关子进程 / 非 agent 会话直接放行，
    所以对现有用法是**零影响**。
    """
    if os.environ.get(AGENT_GATEWAY_ENV) == "act":
        return None
    if cmd in ("act", "approve"):
        return None
    inner, outer = _split_outer([str(x) for x in (raw_argv or [])])
    risk, why = agent_classify(cmd, inner)
    actor, _src = agent_detect(cfg)
    if not actor:
        # 没认出是 agent：人在交互终端里敲的一律放行（行为不变）；非交互终端里的
        # 破坏性命令默认也要批准令牌——AI 调用多半没有 TTY，不能靠它自觉设 AUTOZT_ACTOR。
        if (risk == "destructive" and agent_gate_policy(cfg) != "off"
                and not _stdin_isatty()):
            actor = "non-tty"
        else:
            return None
    sig = agent_signature(cmd, inner)
    if risk == "destructive" and agent_strict():
        ok, approver = agent_take_approval(cfg, sig)
        if not ok:
            rc = agent_request_approval(cfg, cmd, inner, why, actor, "direct",
                                        outer=outer)
            print("  （自己的定时脚本需要非交互执行时：AUTOZT_AGENT_STRICT=0，"
                  "或 tf.yaml 写 agent_gate: off）")
            return rc
        agent_audit(cfg, actor, cmd, inner, risk, "allow-direct-approved",
                    why=why, sig=sig, approved_by=approver or "human",
                    gateway="direct")
        return None
    agent_audit(cfg, actor, cmd, inner, risk, "direct", why=why, sig=sig,
                gateway="direct")
    return None
