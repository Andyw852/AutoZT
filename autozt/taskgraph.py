# -*- coding: utf-8 -*-
"""taskgraph —— 技能任务图：步骤依赖 + 每步当前状态，一屏看清「整条流程、算到哪、卡在哪」。

通用实现：只读采集结果里每个材料的 steps[]（kind / label_txt / diag / job / _deps），
依赖来自 skill.yaml 的 needs（没写 needs 的步骤 = 上一步，与 _dag_recompute 同一口径），
不认识任何具体技能——加新技能 = 加目录，任务图零改动。

三种视图：
  material  单材料：每步一个节点，按依赖画成树；节点带状态、作业号/排队原因/诊断；
  skill     整技能（不带 -p）：同一棵步骤树，节点上是各状态的材料计数；
  card      审批卡片：在 material 视图上标出本次动作的目标步骤（◀ 本次：取消），
            受影响的下游步骤标 ↳（agentgate 用它把"要干什么"摊在用户眼前）。

输出：text（终端；不支持 UTF-8 时自动退回 ASCII 图标）/ json（nodes+edges+text+mermaid）/
mermaid（给能渲染 Markdown 的 agent 客户端）。
"""

import json
import os
import subprocess
import sys

# 状态 → (图标, ASCII 图标, 中文短词)。图标都选单宽字符，终端里对得齐。
KINDS = (
    ("OK", "✔", "+", "完成"),
    ("R", "▶", ">", "运行"),
    ("PD", "◷", "~", "排队"),
    ("FAIL", "✖", "x", "失败"),
    ("TODO", "◇", "o", "可提交"),
    ("PREP", "○", ".", "可开始"),
    ("WAIT", "·", "-", "等上游"),
    ("SCANCEL", "⊘", "/", "已取消"),
    ("IMAG", "△", "!", "虚频"),
    ("OTHER", "?", "?", "其它"),
)
_ICON = {k: u for k, u, _a, _w in KINDS}
_ICON_ASCII = {k: a for k, _u, a, _w in KINDS}
_WORD = {k: w for k, _u, _a, w in KINDS}
_ORDER = [k for k, _u, _a, _w in KINDS]
_MERMAID_CLASS = {"OK": "ok", "R": "run", "PD": "pd", "FAIL": "fail", "TODO": "ready",
                  "PREP": "ready", "WAIT": "wait", "SCANCEL": "cancel", "IMAG": "warn",
                  "OTHER": "wait"}


def _unicode_ok(stream=None):
    stream = stream or sys.stdout
    try:
        "✔├─◷⊘△◀↳".encode(getattr(stream, "encoding", None) or "ascii")
        return True
    except (UnicodeEncodeError, LookupError):
        return False


def _disp_width(s):
    """终端显示宽度（中文/全角算 2）。"""
    w = 0
    for ch in str(s):
        o = ord(ch)
        if (0x1100 <= o <= 0x115F or 0x2E80 <= o <= 0xA4CF or 0xAC00 <= o <= 0xD7A3
                or 0xF900 <= o <= 0xFAFF or 0xFE30 <= o <= 0xFE6F
                or 0xFF00 <= o <= 0xFF60 or 0xFFE0 <= o <= 0xFFE6
                or 0x1F300 <= o <= 0x1FAFF):
            w += 2
        else:
            w += 1
    return w


def _pad(s, w):
    return s + " " * max(0, w - _disp_width(s))


def _kind(s):
    k = str((s or {}).get("kind") or "OTHER")
    return k if k in _ICON else "OTHER"


def _short(text, n=70):
    t = " ".join(str(text or "").split())
    return t if len(t) <= n else t[:n - 1] + "…"


# ===== 图结构 =====
def _layout(nodes):
    """给 nodes（有序，含 label / deps）算 depth、主父节点和子节点表。

    主父节点 = 深度最大的依赖（同深度取靠后的）；其余依赖记在 extra_deps 里，
    文本视图写成"另需 X、Y"。返回 (roots, children)。"""
    by = {n["label"]: n for n in nodes}
    order = {n["label"]: i for i, n in enumerate(nodes)}
    depth = {}

    def _depth(lb, seen=()):
        if lb in depth:
            return depth[lb]
        if lb in seen:          # 环（配置错误）：截断，别死循环
            return 0
        ps = [d for d in by[lb]["deps"] if d in by]
        v = 0 if not ps else 1 + max(_depth(p, seen + (lb,)) for p in ps)
        depth[lb] = v
        return v

    for n in nodes:
        n["depth"] = _depth(n["label"])
    children = {n["label"]: [] for n in nodes}
    roots = []
    for n in nodes:
        ps = [d for d in n["deps"] if d in by and d != n["label"]]
        if not ps:
            roots.append(n["label"])
            n["extra_deps"] = []
            continue
        main = max(ps, key=lambda p: (depth.get(p, 0), order.get(p, 0)))
        n["parent"] = main
        n["extra_deps"] = [p for p in ps if p != main]
        children[main].append(n["label"])
    return roots, children


def _status_text(n):
    """节点的人话状态：作业号/节点/排队原因/诊断/缺哪个上游。"""
    k = n.get("kind") or "OTHER"
    j = n.get("job") or {}
    lt = str(n.get("label_txt") or "")
    fan = ""
    if "/" in lt and lt.split("/", 1)[0].strip().isdigit():   # 扇出步骤 "3/5 2R 0PD"
        fan = " [%s]" % lt
    if k == "R":
        bits = ["#%s" % j["id"]] if j.get("id") else []
        if j.get("info"):
            bits.append("@%s" % j["info"])
        if j.get("time"):
            bits.append("已跑 %s" % j["time"])
        return ("运行 " + " ".join(bits)).strip() + fan
    if k == "PD":
        why = str(j.get("info") or "").strip().strip("()")
        return ("排队 %s%s" % (("#%s" % j["id"]) if j.get("id") else "",
                               (" (%s)" % why) if why else "")).strip() + fan
    if k == "FAIL":
        return "失败：%s" % _short(n.get("diag") or lt or "?", 60) + fan
    if k == "WAIT":
        need = [d for d in (n.get("pending_deps") or [])]
        return ("等上游 %s" % "、".join(need)) if need else "等上游"
    if k == "TODO":
        return "可提交（输入已生成）" + fan
    if k == "PREP":
        return "可开始" + fan
    if k == "SCANCEL":
        return "已取消（不会自动重跑；start 续交）"
    if k == "IMAG":
        return "算完但有虚频（不算失败）"
    if k == "OK":
        return "完成" + fan
    return (lt or "?") + fan


def material_graph(t, m):
    """采集结果中的一个材料 → 任务图（dict，可直接 json.dumps）。"""
    steps = [s for s in (m.get("steps") or []) if isinstance(s, dict)]
    name2label = {s.get("name"): s.get("label") or s.get("name") for s in steps}
    nodes, prev = [], None
    for s in steps:
        lb = s.get("label") or s.get("name")
        raw = s.get("_deps")
        if raw is None:                         # 没有 DAG 信息（旧缓存）：回退上一步
            raw = [prev] if prev else []
        deps = [name2label.get(d, os.path.basename(str(d))) for d in raw if d]
        j = s.get("job") or None
        node = {"name": s.get("name"), "label": lb, "kind": _kind(s),
                "label_txt": s.get("label_txt"), "diag": _short(s.get("diag"), 200),
                "deps": deps, "missing_deps": list(s.get("_missing_deps") or [])}
        if j:
            node["job"] = {k: j.get(k) for k in ("id", "state", "info", "time")
                           if j.get(k) not in (None, "")}
        nodes.append(node)
        prev = s.get("name")
    kinds = {n["label"]: n["kind"] for n in nodes}
    for n in nodes:
        n["pending_deps"] = [d for d in n["deps"] if kinds.get(d) not in (None, "OK")]
        n["status"] = _status_text(n)
    counts = {}
    for n in nodes:
        counts[n["kind"]] = counts.get(n["kind"], 0) + 1
    _layout(nodes)
    return {"view": "material", "skill": (t or {}).get("key") or m.get("tt"),
            "material": m.get("name"), "id": m.get("qualified_name") or m.get("name"),
            "project": m.get("project") or "", "hpc": m.get("hpc_name") or "",
            "done": counts.get("OK", 0), "total": len(nodes), "counts": counts,
            "nodes": nodes,
            "edges": [[d, n["label"]] for n in nodes for d in n["deps"]
                      if d in kinds]}


def skill_graph(t):
    """整技能视图：步骤树 + 每步各状态的材料计数 + 前几条失败。"""
    mats = [m for m in (t.get("materials") or [])]
    nodes, seen = [], {}
    fails = []
    for m in mats:
        g = material_graph(t, m)
        for n in g["nodes"]:
            hit = seen.get(n["label"])
            if hit is None:
                hit = {"name": n["name"], "label": n["label"], "deps": list(n["deps"]),
                       "counts": {}}
                seen[n["label"]] = hit
                nodes.append(hit)
            hit["counts"][n["kind"]] = hit["counts"].get(n["kind"], 0) + 1
            if n["kind"] == "FAIL" and len(fails) < 8:
                fails.append({"material": g["id"], "step": n["label"],
                              "diag": _short(n.get("diag"), 80)})
    _layout(nodes)
    done = sum(1 for m in mats if m.get("steps")
               and all(_kind(s) == "OK" for s in m["steps"]))
    labels = {n["label"] for n in nodes}
    return {"view": "skill", "skill": t.get("key"), "desc": t.get("desc") or "",
            "materials": len(mats), "materials_done": done, "nodes": nodes,
            "edges": [[d, n["label"]] for n in nodes for d in n["deps"] if d in labels],
            "fails": fails}


def descendants(graph, labels):
    """labels 的全部下游步骤（不含自身），按图中顺序。"""
    want = set(labels or [])
    out = []
    changed = True
    while changed:
        changed = False
        for n in graph.get("nodes") or []:
            lb = n["label"]
            if lb in want or lb in out:
                continue
            if any(d in want or d in out for d in n.get("deps") or []):
                out.append(lb)
                changed = True
    order = [n["label"] for n in graph.get("nodes") or []]
    return [x for x in order if x in out]


# ===== 渲染 =====
def _tree_rows(nodes, ascii_mode=False):
    """[(前缀, 节点)]：按主父节点画树（├─ └─ │）。"""
    if not nodes:
        return []
    for n in nodes:
        n.setdefault("deps", [])
    roots, children = _layout(nodes)
    by = {n["label"]: n for n in nodes}
    tee, ell, bar = ("|-", "`-", "| ") if ascii_mode else ("├─", "└─", "│ ")
    rows = []

    def walk(lb, first, rest):
        rows.append((first, by[lb]))
        kids = children.get(lb) or []
        for i, k in enumerate(kids):
            last = i == len(kids) - 1
            walk(k, rest + (ell if last else tee), rest + ("  " if last else bar))

    for r in roots:
        walk(r, "", "")
    return rows


def render_text(graph, marks=None, ascii_mode=None, legend=True, indent=""):
    """任务图 → 终端文本。marks = {步骤 label: 标注文字}（审批卡片用）。"""
    if ascii_mode is None:
        ascii_mode = not _unicode_ok()
    icon = _ICON_ASCII if ascii_mode else _ICON
    marks = marks or {}
    lines = []
    view = graph.get("view")
    nodes = [dict(n) for n in graph.get("nodes") or []]
    if view == "skill":
        lines.append("%s%s：%d 个材料，全部完成 %d 个"
                     % (indent, graph.get("skill"), graph.get("materials", 0),
                        graph.get("materials_done", 0)))
    else:
        c = graph.get("counts") or {}
        extra = "，".join("%s %d" % (_WORD[k], c[k]) for k in ("R", "PD", "FAIL", "SCANCEL")
                         if c.get(k))
        where = " · ".join(x for x in (graph.get("skill"), graph.get("id") or graph.get("material"),
                                       graph.get("hpc")) if x)
        lines.append("%s%s   %d/%d 完成%s" % (indent, where, graph.get("done", 0),
                                            graph.get("total", 0),
                                            ("，" + extra) if extra else ""))
    rows = _tree_rows(nodes, ascii_mode=ascii_mode)
    heads = []
    for pre, n in rows:
        if view == "skill":
            heads.append("%s%s" % (pre, n["label"]))
        else:
            heads.append("%s%s %s" % (pre, icon.get(n.get("kind"), "?"), n["label"]))
    w = max([_disp_width(h) for h in heads] + [0]) + 2
    for (pre, n), h in zip(rows, heads):
        if view == "skill":
            cnt = n.get("counts") or {}
            body = " ".join("%s%d" % (icon[k], cnt[k]) for k in _ORDER if cnt.get(k))
        else:
            body = n.get("status") or _status_text(n)
        extra = []
        if n.get("extra_deps") and n.get("kind") != "WAIT":   # WAIT 的状态里已列出未完成的上游
            extra.append("（另依赖 %s）" % "、".join(n["extra_deps"]))
        if n.get("missing_deps"):
            extra.append("（缺 %s：可选步骤已关）" % "、".join(n["missing_deps"]))
        mk = marks.get(n["label"])
        line = "%s%s%s%s" % (indent, _pad(h, w), body, "".join(extra))
        if mk:
            line += "   " + mk
        lines.append(line.rstrip())
    if view == "skill":
        for f in graph.get("fails") or []:
            lines.append("%s  %s %s %s  %s" % (indent, icon["FAIL"], f["material"],
                                              f["step"], f["diag"]))
    if legend:
        lines.append(indent + "图例 " + " ".join(
            "%s%s" % (icon[k], _WORD[k]) for k in _ORDER if k != "OTHER"))
    return "\n".join(lines)


def _mermaid_label(*parts):
    """各段分别转义后用 <br/> 换行（先拼再转义会把 <br/> 也转掉）。"""
    esc = [str(x).replace('"', "#quot;").replace("<", "&lt;").replace(">", "&gt;")
           for x in parts if x not in (None, "")]
    return "<br/>".join(esc)


def render_mermaid(graph, marks=None):
    """任务图 → Mermaid flowchart（能渲染 Markdown 的 agent 客户端里就是一张图）。"""
    marks = marks or {}
    nodes = graph.get("nodes") or []
    ids = {n["label"]: "n%d" % i for i, n in enumerate(nodes)}
    out = ["flowchart TD"]
    for n in nodes:
        if graph.get("view") == "skill":
            cnt = n.get("counts") or {}
            parts = [n["label"], " ".join("%s%d" % (_ICON[k], cnt[k])
                                          for k in _ORDER if cnt.get(k))]
        else:
            parts = ["%s %s" % (_ICON.get(n.get("kind"), "?"), n["label"]),
                     n.get("status") or ""]
        parts.append(marks.get(n["label"]))
        out.append('  %s["%s"]' % (ids[n["label"]], _mermaid_label(*parts)))
    for a, b in graph.get("edges") or []:
        if a in ids and b in ids:
            out.append("  %s --> %s" % (ids[a], ids[b]))
    out += ["  classDef ok fill:#dcfce7,stroke:#15803d,color:#14532d",
            "  classDef run fill:#dbeafe,stroke:#1d4ed8,color:#1e3a8a",
            "  classDef pd fill:#fef9c3,stroke:#a16207,color:#713f12",
            "  classDef fail fill:#fee2e2,stroke:#b91c1c,color:#7f1d1d",
            "  classDef ready fill:#f3f4f6,stroke:#4b5563,color:#111827",
            "  classDef wait fill:#ffffff,stroke:#9ca3af,color:#6b7280,stroke-dasharray:3 3",
            "  classDef cancel fill:#f5f5f4,stroke:#57534e,color:#44403c",
            "  classDef warn fill:#ffedd5,stroke:#c2410c,color:#7c2d12",
            "  classDef target stroke:#dc2626,stroke-width:3px"]
    if graph.get("view") != "skill":
        for n in nodes:
            out.append("  class %s %s" % (ids[n["label"]],
                                          _MERMAID_CLASS.get(n.get("kind"), "wait")))
    for lb in marks:
        if lb in ids and str(marks[lb]).startswith("◀"):
            out.append("  class %s target" % ids[lb])
    return "\n".join(out)


def with_renderings(graph, marks=None):
    g = dict(graph)
    g["text"] = render_text(graph, marks=marks, ascii_mode=False)
    g["mermaid"] = render_mermaid(graph, marks=marks)
    if marks:
        g["marks"] = dict(marks)
    return g


# ===== 命令入口 =====
def cmd_graph(data, projs=None, json_out=False, mermaid=False):
    """autozt [-tt T] [-p M] graph [--json|--mermaid]：带 -p 画单材料，不带画整技能。"""
    from autozt.report import find_material
    graphs = []
    if projs:
        for pj in projs:
            t, m = find_material(data, pj)
            graphs.append(material_graph(t, m))
    else:
        for t in data.get("types") or []:
            if t.get("materials"):
                graphs.append(skill_graph(t))
    if not graphs:
        print("没有可画的材料（先 register/init 材料，或检查 -tt/--project）。")
        return 1
    if json_out:
        out = [with_renderings(g) for g in graphs]
        print(json.dumps(out[0] if len(out) == 1 else out, ensure_ascii=False, indent=1))
        return 0
    if mermaid:
        print("\n\n".join("```mermaid\n%s\n```" % render_mermaid(g) for g in graphs))
        return 0
    asc = not _unicode_ok()
    print("\n\n".join(render_text(g, ascii_mode=asc, legend=(i == len(graphs) - 1))
                      for i, g in enumerate(graphs)))
    return 0


def fetch_graphs(base_argv, tt=None, materials=None, timeout=90, env=None):
    """在子进程里跑只读的 `autozt graph --json`（审批卡片、MCP 用）。

    base_argv = 启动 autozt 的命令前缀（含 -c 配置）。失败返回 (None, 原因)。"""
    argv = list(base_argv)
    if tt:
        argv += ["-tt", str(tt)]
    if materials:
        argv += ["-p", ",".join(str(x) for x in materials)]
    argv += ["graph", "--json"]
    e = dict(os.environ if env is None else env)
    e["AUTOZT_AGENT_GATEWAY"] = "act"       # 只读子调用：不重复记审计
    try:
        p = subprocess.run(argv, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, env=e)
    except subprocess.TimeoutExpired:
        return None, "采集超时（%ss）" % timeout
    except OSError as ex:
        return None, str(ex)
    if p.returncode != 0:
        msg = ((p.stdout or "") + (p.stderr or "")).strip().splitlines()
        return None, (msg[-1] if msg else "rc=%s" % p.returncode)[:300]
    try:
        got = json.loads(p.stdout)
    except ValueError:
        return None, "graph 输出不是 JSON"
    return (got if isinstance(got, list) else [got]), None
