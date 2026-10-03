"""对话内审批：请求号 + 审批卡片 → 用户同意 → agent 代为执行；MCP 确认框（elicitation）。

纯本地：cmd_act / MCP 的 CLI 调用都被替换掉，不采集、不连超算。
"""
import contextlib
import io
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from autozt import agentgate as G  # noqa: E402
from autozt import cli as C  # noqa: E402
from autozt import mcp as M  # noqa: E402


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    for k in ("AUTOZT_APPROVAL_MODE", "AUTOZT_APPROVE_TTL", "AUTOZT_AGENT_EVENTS",
              "AUTOZT_AGENT_GATEWAY"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("AUTOZT_ACTOR", "test-agent")
    monkeypatch.setattr(G, "_stdin_isatty", lambda: False)
    return {"_config_dir": str(tmp_path)}     # 没有 _config_path：卡片不去采任务图


def _quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = fn(*a, **k)
    return rc, buf.getvalue()


def _deny(cfg, argv):
    """agent 直接执行破坏性命令 → 被拦下，返回 (rc, 输出, 请求号)。"""
    rc, out = _quiet(G.agent_direct_gate, cfg, argv[-1] if argv[-1] != "-y" else argv[-2],
                     argv)
    reqs = G.agent_requests_pending(cfg)
    return rc, out, (reqs[-1]["id"] if reqs else None)


def test_blocked_command_issues_request_and_card(cfg):
    rc, out, rid = _deny(cfg, ["-tt", "fit-fc-thermal", "-p", "MnIn2Se4", "-j",
                               "step2_kappa", "stop", "-y"])
    assert rc == 3 and rid and rid.startswith("r")
    assert "需要你的同意" in out and "后果" in out and rid in out
    assert 'autozt approve --request %s --reply' % rid in out
    # 同一条命令再被拦一次：复用同一个请求号，不刷屏
    rc2, _out, rid2 = _deny(cfg, ["-tt", "fit-fc-thermal", "-p", "MnIn2Se4", "-j",
                                  "step2_kappa", "stop"])
    assert rid2 == rid and len(G.agent_requests_pending(cfg)) == 1
    # -c 等全局选项不进签名（直接调用与 act 的同一条命令是同一个签名）
    inner, outer = G._split_outer(["-c", "/x/tf.yaml", "-p", "A", "stop"])
    assert inner == ["-p", "A", "stop"] and outer == ["-c", "/x/tf.yaml"]


def test_chat_approval_runs_command_with_y_and_audits_reply(cfg, monkeypatch):
    _rc, _out, rid = _deny(cfg, ["-tt", "fit-fc-thermal", "-p", "Mn", "-j", "S2_kappa",
                                 "stop"])
    ran = []

    def fake_act(c, raw):
        ran.append(raw)
        ok, by = G.agent_take_approval(c, G.agent_signature("stop", raw[raw.index("act") + 1:]))
        assert ok and by == "user(chat)"      # 批准令牌和执行的是同一条命令
        return 0

    monkeypatch.setattr(G, "cmd_act", fake_act)
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "--request", rid])
    assert rc == 2 and "--reply" in out and not ran          # 没有用户原话：不批准
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "--request", rid, "--reply", "同意，停掉"])
    assert rc == 0 and ran and ran[0][-1] == "-y" and "stop" in ran[0]
    assert G.agent_request_get(cfg, rid)["status"] == "used"
    rec = [json.loads(x) for x in open(G.agent_log_path(cfg), encoding="utf-8")]
    hit = [r for r in rec if r["decision"] == "approved-in-chat"]
    assert hit and hit[-1]["reply"] == "同意，停掉" and hit[-1]["request"] == rid
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "--request", rid, "--reply", "同意"])
    assert rc == 3 and "执行过" in out and len(ran) == 1       # 一次有效，不能重放


def test_negative_reply_deny_expiry_and_tty_mode(cfg, monkeypatch):
    monkeypatch.setattr(G, "cmd_act", lambda c, raw: pytest.fail("不应执行"))
    _rc, _o, rid = _deny(cfg, ["-p", "A", "rerun"])
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "--request", rid, "--reply", "不"])
    assert rc == 3 and G.agent_request_get(cfg, rid)["status"] == "declined"

    _rc, _o, rid = _deny(cfg, ["-p", "B", "clean"])
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "--deny", rid])
    assert rc == 0 and G.agent_request_get(cfg, rid)["status"] == "declined"

    _rc, _o, rid = _deny(cfg, ["-p", "C", "clean"])
    G.agent_request_update(cfg, rid, expires_epoch=1)
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "--request", rid, "--reply", "同意"])
    assert rc == 3 and "过期" in out

    monkeypatch.setenv("AUTOZT_APPROVAL_MODE", "tty")
    _rc, out, rid = _deny(cfg, ["-p", "D", "clean"])
    assert "严格模式" in out
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "--request", rid, "--reply", "同意"])
    assert rc == 3 and "交互终端" in out
    # 旧写法：非交互终端照旧拒绝
    rc, out = _quiet(G.cmd_approve, cfg, ["approve", "-p", "D", "clean"])
    assert rc == 3 and "交互终端" in out


def test_approval_event_line_for_executors(cfg, monkeypatch):
    monkeypatch.setenv("AUTOZT_AGENT_EVENTS", "1")
    rc, out, rid = _deny(cfg, ["-p", "A", "-j", "S1", "stop"])
    ev = M._approval_event(out)
    assert ev["request_id"] == rid and ev["approval_mode"] == "chat"
    assert ev["command"] == "autozt -p A -j S1 stop" and "scancel" in ev["consequences"]


def test_route_does_not_mistake_approve_values_for_subcommands(monkeypatch):
    monkeypatch.setattr(M, "main", lambda argv=None: pytest.fail("误路由到 mcp"))
    assert C.route_subcommand(["approve", "--request", "r1", "--reply", "mcp",
                               "--via", "mcp"]) is None
    assert C.route_subcommand(["act", "-p", "agent", "stop"]) is None


# ---------------------------------------------------------------- MCP
_EV = {"request_id": "rabc123", "sig": "x", "cmd": "stop", "approval_mode": "chat",
       "command": "autozt -p A -j S1 stop -y", "consequences": "scancel …",
       "expires_in_s": 900, "card": "━━ 需要你的同意 ━━\n任务图 …",
       "graphs": [{"text": "A   0/1 完成\n▶ S1  运行   ◀ 本次：取消", "mermaid": "flowchart TD"}]}


@pytest.fixture
def mcp_calls(monkeypatch):
    calls = []

    def fake_run(argv, timeout=1800):
        calls.append(list(argv))
        if argv[:1] == ["act"]:
            return 3, "✗ 拒绝执行\n[agent-approval] " + json.dumps(_EV, ensure_ascii=False), ""
        if argv[:1] == ["approve"]:
            return 0, "✓ 已批准（mcp）。执行：autozt -p A -j S1 stop -y\nscancel 1 成功", ""
        if "graph" in argv:
            return 0, json.dumps({"view": "material", "text": "TEXT", "mermaid": "MER"}), ""
        return 1, "", "unexpected"

    monkeypatch.setattr(M, "_run", fake_run)
    monkeypatch.setattr(M, "_CLIENT", {"elicitation": False, "protocol": "2025-06-18"})
    monkeypatch.setattr(M, "_IO", {"q": None, "backlog": [], "seq": 0, "current": None})
    return calls


def test_mcp_request_returns_card_then_approve_executes(mcp_calls):
    got = M.call_tool("request_destructive_action",
                      {"action": "cancel", "material": "A", "step": "S1"})
    data = got["structuredContent"]["data"]
    assert got["isError"] and data["needs_approval"] and data["request_id"] == "rabc123"
    assert "需要你的同意" in got["content"][0]["text"] and "approve_request" in data["next"]
    assert mcp_calls[0] == ["act", "-p", "A", "-j", "S1", "stop", "-y"]
    ok = M.call_tool("approve_request", {"request_id": "rabc123", "reply": "同意"})
    assert not ok["isError"] and ok["structuredContent"]["data"]["approved"]
    assert mcp_calls[-1] == ["approve", "--request", "rabc123", "--reply", "同意",
                             "--via", "mcp"]
    M.call_tool("approve_request", {"request_id": "rabc123", "reply": "不要",
                                    "decision": "deny"})
    assert mcp_calls[-1][:3] == ["approve", "--deny", "rabc123"]
    g = M.call_tool("task_graph", {"material": "A", "format": "mermaid"})
    assert g["content"][0]["text"] == "```mermaid\nMER\n```"
    assert mcp_calls[-1] == ["-p", "A", "graph", "--json"]


def _elicit(monkeypatch, answer):
    import queue
    q = queue.Queue()
    sent = []
    monkeypatch.setattr(M, "_write", lambda msg: sent.append(msg))
    monkeypatch.setattr(M, "_CLIENT", {"elicitation": True, "protocol": "2025-06-18"})
    monkeypatch.setattr(M, "_IO", {"q": q, "backlog": [], "seq": 0, "current": 7})
    monkeypatch.setattr(M, "_approval_wait", lambda: 0.3)
    q.put(json.dumps({"jsonrpc": "2.0", "id": 99, "method": "tools/list"}) + "\n")
    q.put(json.dumps({"jsonrpc": "2.0", "id": "p", "method": "ping"}) + "\n")
    if answer is not None:
        q.put(json.dumps({"jsonrpc": "2.0", "id": "autozt-1", "result": answer}) + "\n")
    got = M.call_tool("request_destructive_action",
                      {"action": "cancel", "material": "A", "step": "S1"})
    return got, sent


def test_mcp_elicitation_accept_decline_and_timeout(mcp_calls, monkeypatch):
    got, sent = _elicit(monkeypatch, {"action": "accept", "content": {"confirm": True}})
    assert sent[0]["method"] == "elicitation/create"
    assert "◀ 本次：取消" in sent[0]["params"]["message"]
    assert {"jsonrpc": "2.0", "id": "p", "result": {}} in sent          # ping 当场回
    assert M._IO["backlog"] and "tools/list" in M._IO["backlog"][0]     # 其它请求排队
    assert not got["isError"]
    assert mcp_calls[-1][-2:] == ["--via", "mcp-elicitation"]

    got, _ = _elicit(monkeypatch, {"action": "decline"})
    assert mcp_calls[-1][:2] == ["approve", "--deny"]
    assert got["structuredContent"]["data"]["declined"]

    got, _ = _elicit(monkeypatch, {"action": "accept", "content": {"confirm": False}})
    assert mcp_calls[-1][:2] == ["approve", "--deny"]                   # 没勾「执行」

    n = len(mcp_calls)
    got, _ = _elicit(monkeypatch, None)                                   # 用户没理：超时
    data = got["structuredContent"]["data"]
    assert data["needs_approval"] and "timeout" in data["elicitation"]
    assert len(mcp_calls) == n + 1                                        # 只有被拦的那次


def test_initialize_records_elicitation_capability(monkeypatch):
    monkeypatch.setattr(M, "_CLIENT", {"elicitation": False, "protocol": None})
    M.handle({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2025-06-18",
                         "capabilities": {"elicitation": {}}}})
    assert M._CLIENT["elicitation"] is True
    M.handle({"jsonrpc": "2.0", "id": 2, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05",
                         "capabilities": {"elicitation": {}}}})
    assert M._CLIENT["elicitation"] is False                            # 旧协议没有确认框
