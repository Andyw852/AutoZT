"""MCP 接口测试（不需要集群）：工具表、协议握手、风险分级、网关拒绝。

MCP 是 PhonoAgent 对外给 agent 用的唯一入口，所以这层必须锁死：
  1. 工具是通用动词，数量封顶 20（不许按技能加工具）；
  2. 危险动作必须在协议边界被拒绝并给出批准命令（agent 不能自我批准）；
  3. 只读调用不触发任何提交。
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROG = os.path.join(ROOT, "bin", "phonoagent")
sys.path.insert(0, ROOT)
from phonoagent import mcp as M  # noqa: E402


def test_tool_table_is_generic_and_capped():
    assert 10 <= len(M.TOOLS) <= 20, "工具数应恒定且封顶 20"
    names = [t[0] for t in M.TOOLS]
    assert len(names) == len(set(names))
    # 不许出现"某技能专用"的工具名
    for n in names:
        assert "dft" not in n and "mace" not in n and "phonon" not in n, n


def test_risk_tiers_cover_three_levels():
    risks = {t[3] for t in M.TOOLS}
    assert risks == {"read", "mutate", "destructive"}


def test_every_tool_has_schema():
    for name, desc, schema, risk, verb in M.TOOLS:
        assert schema.get("type") == "object", name
        assert desc and verb, name


def test_initialize_and_tools_list():
    p = subprocess.run([sys.executable, PROG, "mcp"], input=
                       json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize"}) + "\n" +
                       json.dumps({"jsonrpc": "2.0", "id": 2, "method": "tools/list"}) + "\n",
                       capture_output=True, text=True, cwd=ROOT, timeout=120)
    lines = [json.loads(x) for x in p.stdout.splitlines() if x.strip()]
    assert lines[0]["result"]["serverInfo"]["name"] == "phonoagent"
    assert len(lines[1]["result"]["tools"]) == len(M.TOOLS)


def test_unknown_method_is_an_error():
    resp = M.handle({"jsonrpc": "2.0", "id": 9, "method": "nope"})
    assert resp["error"]["code"] == -32601


def test_destructive_tool_goes_through_gateway():
    # 不连集群：只断言这条路走的是 act 网关（agent 会话下会被拒）
    rc, out, err = M._run(["act", "-tt", "opt-dft-cpu", "-p", "Si_none", "clean", "-y"])
    text = out + err
    assert rc != 0, "破坏性动作在 agent 会话下必须被拒"
    assert "approve" in text or "批准" in text, text[:200]


def test_resources_list_exposes_skill_docs():
    items = M._resources_list()
    uris = [r["uri"] for r in items]
    assert items and any(u.startswith("phonoagent://skill/") for u in uris)
    assert any(u.endswith("/skill.yaml") for u in uris)


def test_resources_read_returns_file_content():
    uri = "phonoagent://doc/README.en.md"
    got = M._resource_read(uri)
    assert got and got["contents"][0]["text"].startswith("# PhonoAgent")


def test_resources_read_refuses_traversal_and_unknown():
    assert M._resource_read("phonoagent://doc/../../etc/passwd") is None
    assert M._resource_read("phonoagent://doc/does-not-exist.md") is None
    assert M._resource_read("http://example.com/x") is None


def test_handle_serves_resource_methods():
    resp = M.handle({"jsonrpc": "2.0", "id": 7, "method": "resources/list"})
    assert resp["result"]["resources"]
    bad = M.handle({"jsonrpc": "2.0", "id": 8, "method": "resources/read",
                    "params": {"uri": "phonoagent://doc/nope.md"}})
    assert bad["error"]["code"] == -32602


def test_prompts_list_and_get():
    items = M._prompts_list()
    assert {p["name"] for p in items} == {"triage-failures", "review-before-submit"}
    got = M._prompt_get("review-before-submit", {"material": "Si_demo", "step": "S1_opt"})
    assert "Si_demo" in got["messages"][0]["content"]["text"]
    assert M._prompt_get("nope") is None
    resp = M.handle({"jsonrpc": "2.0", "id": 9, "method": "prompts/get",
                     "params": {"name": "nope"}})
    assert resp["error"]["code"] == -32602


def test_readonly_profile_hides_hazardous_tools(monkeypatch=None):
    import os
    os.environ["PHONOAGENT_MCP_READONLY"] = "1"
    try:
        names = [t["name"] for t in M._tools_list()]
        assert names, "只读档不该为空"
        assert all(n in {x[0] for x in M.TOOLS if x[3] == "read"} for n in names)
        res = M.call_tool("clean_material", {"material": "Si_x"})
        assert res["isError"], "只读档下危险工具必须被拒"
    finally:
        os.environ.pop("PHONOAGENT_MCP_READONLY", None)
    full = [t["name"] for t in M._tools_list()]
    assert len(full) > len(names), "默认档应比只读档暴露更多工具"
