# MCP interface

PhonoAgent can expose itself as a Model Context Protocol (MCP) server so that an
LLM-driven planner (or any MCP-capable client) can drive the workflow without bespoke
glue code, while every dangerous action still passes through the approval gate.

    phonoagent mcp                 # JSON-RPC 2.0 over stdio: initialize / tools/list / tools/call
    phonoagent mcp --list-tools    # print the tool table (for humans and tests)
    phonoagent mcp --call NAME JSON

## Design rules

1. The server only calls the command-line interface. Approval, audit and provenance have
   a single source of truth and cannot drift from the CLI.
2. Tools are generic verbs, never per-skill commands. Skills and their steps are
   discovered from skill.yaml, so adding a skill never adds a tool. The table is capped
   at 20 entries on purpose.
3. Risk is enforced at the protocol boundary. Read-only tools run directly; mutating and
   destructive tools are routed through the action gateway, so an agent session receives
   "needs approval" instead of executing, and approvals only work from a real TTY.

## Tools

| Tool | Risk | Purpose |
|---|---|---|
| list_skills | read | skills and their steps, from skill.yaml |
| list_materials | read | materials and step states |
| get_summary | read | per-skill counters, FAIL list, global queue |
| get_status | read | one material in detail (diagnosis, cluster, work_dir source) |
| get_json | read | full structured state with schema_version |
| conf_get | read | effective step.conf of one step |
| start_step | mutate | generate inputs if needed, then submit |
| retry_step | mutate | regenerate inputs, keep products, do not submit |
| fetch_results | mutate | pull finished results back |
| conf_set | mutate | write a project configuration value |
| session_export | mutate | reproducible session archive |
| stop_step | destructive | cancel the job of a step |
| rerun_step | destructive | delete the step directory and regenerate |
| clean_material | destructive | delete generated files, back to PREP |

## Measured behaviour

safety_metrics.py reports the numbers that can be cited:

    interception rate of hazardous actions  100% (8/8)
    hazardous actions actually executed        0
    every denial offers an approval path     yes
    read-only actions still allowed           yes
    local output determinism                 100%

Run it with

    python scripts/safety_metrics.py          # human readable
    python scripts/safety_metrics.py --json   # machine readable (CI asserts on this)
