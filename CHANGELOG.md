# Changelog

All notable changes to AutoZT are documented here.
This project adheres to Semantic Versioning; versions before 1.0.0 are development snapshots.

## [1.2.0] - 2026-09-16

### Changed
- Renamed the software to AutoZT: package autozt, commands autozt (short form az), environment
  variables AUTOZT_*, configuration directory ~/.config/autozt. The previous command names
  (pa, phonoagent) remain as compatibility wrappers pointing at the same program, so existing
  scripts keep working.
- Version banner and metadata report AutoZT 1.2.0.

## [Unreleased]

### Added
- docs/mcp.md: the MCP interface (tool table, risk tiers, measured safety numbers).
- AI integration (from a 295-material production run where every problem was in discovery,
  name resolution or status retrieval, not in scheduling):
  - `.tf_progress.json` + `autozt progress` / `autozt agent progress` / MCP `get_progress`:
    per-material state, structured FAIL codes (`class`/`code`/`action`/`fail_count`),
    `since`/`eta_s`, per-project counts and monitor liveness, read from a local file with
    no collection and no ssh. Written by the monitor every round and by any collecting
    command; partial scopes merge instead of truncating the full view.
  - `autozt doctor` / `agent doctor` / MCP `doctor`: ssh-free configuration preflight
    (effective `max_jobs` per project, tuning-knob sources, blocked projects, cross-project
    duplicate names and shadowing by blocked projects, monitor/cron keep-alive, agent gate).
  - `--project NAME` scope on every command, agent op and MCP state tool; segments are
    pruned before collection.
  - Stable material id `<project>/<full name>` (`qualified_name`, `id` in agent JSON),
    accepted by `-p` and the `material` argument.
  - `auto on --dry-run`; `auto on` skips files already at the target value.
  - Structured diagnosis codes for segmented relaxations (`segment_zbrent`,
    `segment_watchdog`, `segment_interrupted`, `segments_incomplete`, `segments_exhausted`),
    pressure, walltime, OOM, disk-full, missing outputs, stale inputs, plus a coarse
    `diag_class` and an `action_hint` for parameter tuning.
  - Tuning knobs `collect_chunk`, `collect_workers`, `op_workers`, `init_workers`,
    `cache_ttl` can live in tf.yaml (environment variables still override).

### Changed
- Project-level `max_jobs` in `tf_*.yaml` is now enforced as a per-project cap in addition
  to the skill-level cap. It was previously ignored silently. **Check `autozt doctor`
  before restarting a monitor: a project that sets `max_jobs: 16` now really runs at 16.**
- The agent approval gate is on by default (`agent_gate: auto`): besides `AUTOZT_ACTOR`,
  known AI-agent environment markers and non-interactive stdin make destructive commands
  require an approval token. `agent_gate: off` / `AUTOZT_AGENT_GATE=off` restores the old
  behaviour; `AUTOZT_AGENT_STRICT=0` still disables tokens for your own scripts.
- Monitor output: one `[round] {json}` line per round plus at most 20 step changes; the full
  table is printed only with `monitor_table: true`.
- A basename that also exists in a blocked (stale) project no longer aborts the command when
  an unblocked material with that name exists; errors now name the blocking project and
  its reason, and list all rejected targets at once.
- Legacy skill-name warnings are aggregated by rename pair once there are more than ten.

### Fixed
- E2BIG during collection: the collector payload is sent on stdin instead of as a
  `--config64` argument (Linux caps one argument at 128 KiB). The sbatch guard uses the
  same transport.
- `autozt agent …` and `autozt mcp` did not work through the pip-installed `autozt`
  command (`agent` was rejected; `mcp` collected all state first). All entry points now
  share one router. Agent/MCP subprocess calls fall back to `python -m autozt` when
  `bin/autozt` is not installed.
- The monitor never wrote `history.jsonl` (`history_record` was not imported and the
  NameError was swallowed).
- A single failing monitor round (ssh timeout, `sys.exit` inside collection) no longer kills
  the daemon; the error is logged and the next round retries.
- `monitor -d` checks the pid file before scanning project configurations, so a cron
  keep-alive no longer spends minutes rescanning a slow 9p mount when the daemon is alive.
- Blocked-project conflict checks cache material discovery per project (previously
  re-scanned for every target × blocked config) and see three-level layouts.
- stdout is line-buffered when redirected, so `autozt … > log` shows progress.
- `skills`, `schema`, `skill` and untargeted `history` no longer scan every project
  configuration before running (1–2 minutes per call on a WSL 9p mount). This also speeds
  up `autozt agent skills/contract/research_plan` and MCP `list_skills`/`describe_skill`/
  `research_plan`, which all call `schema --json`. With an explicit `-p` target the
  blocked-project check still runs first.
- Loading project configurations with many blocked projects was slow again (>10 minutes on
  a 9p mount with ~30 blocked projects): the three-level blocked-material discovery used a
  `*/*/*/POSCAR` glob that listed the contents of every material directory. It now walks
  level by level and never descends into a material directory (≈18× fewer filesystem
  operations on a 30×295 synthetic layout, same materials found).
- `history_record` merges its state file across scopes instead of overwriting it, so two
  differently scoped monitors (or a monitor plus a scoped CLI call) no longer erase each
  other's baselines and lose transitions and FAIL counts.
- Monitor `[round]` lines carry `timing_s` (config check / collect / cache+history / fetch /
  hang check / advance) and the global `auto_advance` flag, so "nothing submitted" is not
  mistaken for "stuck".
- Monitor rounds do less local filesystem work (first real round on a 9p mount: 48.7 min,
  of which `cfg_check` 14.4 min and `collect` 32.3 min):
  - the configuration-change walk no longer descends into `step*` directories of materials
    and their skill subfolders (only computation outputs and fanout `disp-*` folders live
    there); ~5.5× fewer filesystem operations on a synthetic tree, same config files tracked;
  - segments sharing a `local_root` are discovered once per round;
  - `fill_local_dim` caches the POSCAR-based dimension guess by file stat across rounds;
  - `[round]` gains `collect_detail` (local resolve vs parallel ssh vs post-processing
    seconds, segment/material/ssh-call counts) and a separate `fill_dim` timing, so the
    local-9p and ssh parts of `collect` can be told apart.

### Fixed
- Cluster-switch self-checks, each backed by a real failure observed in testing:
  work_dir provenance is printed by status; resource requests above the cluster's
  declared max_cpus are flagged; a partition the target cluster does not use is
  reported with the partitions its own templates use; and conda environments or
  conda.sh paths that do not exist on the target cluster are flagged before
  submission (this was the cause of a silent stage death: activation failed, the
  preparation stage still passed and the fit stage exited without a message).
- start now prints the exact command needed to resubmit a FAIL step
  (autozt -tt SKILL -p MATERIAL -j STEP start -f) instead of claiming that all
  steps are running, finished or dependency-blocked.
- docs/CONFIGURING.md section 5 documents the four cluster-switch traps.

## [1.1.0] - 2026-09-15

### Added
- MCP interface: autozt mcp exposes the tool as a JSON-RPC 2.0 stdio server
  (initialize / tools/list / tools/call). The tool table holds 14 generic verbs
  (list_skills, list_materials, get_summary, get_status, get_json, conf_get, start_step,
  retry_step, fetch_results, conf_set, stop_step, rerun_step, clean_material,
  session_export) and is capped at 20: skills are discovered from skill.yaml, so adding a
  skill never adds a tool.
- MCP policy boundary: read-only tools call the CLI directly; mutating and destructive
  tools are routed through the action gateway (autozt act), so an agent session gets
  needs-approval instead of execution, and approvals only work from a real TTY.
- work_dir provenance: every status output names the configuration layer that supplied
  work_dir, and a failing generation prints that layer plus the hint that
  project_setting/setting.yaml outranks hpc.yaml after a cluster switch.
- tests/test_mcp.py: tool-table shape and cap, risk tiers, protocol handshake, and the
  gateway refusal for destructive calls.

### Fixed
- Status did not reveal which configuration layer supplied work_dir (a cluster switch
  could silently keep writing to the old cluster's path).

## [1.0.0] - 2026-09-15

First public release, derived from the taskflow code base (v1.0 refactor).

### Added
- Package layout autozt/ with the autozt CLI (zero third-party dependency core).
- Criterion-driven step engine: every step declares a physical convergence criterion and
  the engine advances, retries or reports accordingly.
- Provenance per step plus autozt prove --verify, and autozt session export for
  reproducible archival (manifest with per-file sha256, history timeline, replay hints).
- Action gateway (autozt act / approve) with risk tiers, audit log and TTL-bounded
  human approvals for agent-driven sessions.
- Unified core count control: one "cores" setting normalises submit scripts and INCAR
  (NCORE/KPAR) for every skill, whatever template names it uses.
- Unified pre-submission input check: required inputs are discovered from the remote step
  directory instead of per-skill hard-coded lists.
- Supercell specification accepts both diagonal notation ("4 4 4") and a general 3x3
  integer matrix ("2 1 0 -1 2 0 0 0 1", phonopy/phono3py --dim semantics).
- 18 skills covering VASP (band, elastic, defect, ke, kl, opt, phonon), MACE
  (kl/opt/phonon CPU and GPU, mlff), plus te-screen and unihamgnn.
- Drop-in skill interface with autozt schema --strict validation.

### Fixed
- vasp_ncl submit template loaded an Intel module environment while executing the
  locally built (AOCC/AOCL) binary, failing with a missing libscalapack; the template now
  uses the matching environment, plus a single-node MPI guard for nodes with a broken
  mlx5/UCX stack.
- defect-dft-cpu shipped ENCUT = auto which VASP rejects (it now uses the documented 370 eV).
- unihamgnn could never submit its graph step (non-VASP inputs were checked against the
  VASP default list).
