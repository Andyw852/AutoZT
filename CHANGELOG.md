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
- Supercell symmetry guard (`SUPERCELL_SYMMETRY = strict | warn | off`, default strict):
  the displacement generators (kl-dft-cpu S4, phonon-dft-cpu S2, phonon-mlff S2, kl-mlff S2,
  zt-dft-cpu SK4), the fit-fc-thermal S1 fit and S2_kappa stop when the supercell does not
  keep every point-group operation of the crystal, and suggest the smallest symmetric
  diagonal supercell.  Mn2In2Se5 (R-3m) on a 3x3x1 supercell kept 4 of 12 operations; every
  fitting engine then produced a symmetry-broken fc2 (spurious zone-interior imaginary
  modes, kappa jumping by 20-30 % between q meshes).  `symmetry_audit.supercell_symmetry_gate`.
- fit-fc-thermal S1 gate: fc2 symmetry check (`FC2_SYM_TOL`, frequency spread between
  symmetry-equivalent q points; new status `fc2_symmetry_broken`, shown by autozt as
  "fc2 breaks crystal symmetry") and a denser full mesh (`IMAG_MESH_LENGTH = 100`, no
  symmetry reduction) instead of `run_mesh(60.0)` with reduction.
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
- `_converge_mesh` compares each mesh with the latest earlier mesh that is coarser along every
  axis (axes frozen below `MESH_CONV_MAX_LENGTH` exempt), not with the previous mesh: with
  `L *= 1.25` the long axis of an anisotropic cell only steps every few meshes, and two meshes
  with the same k_z division were declared converged without k_z ever being tested.
  `mesh_convergence.json` records `axis_counts`, `compared_with` and `frozen_axes`.
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
- fit-fc-thermal S1 gate: the fc2 (dataset SPOSCAR order) is reordered into the supercell order
  phonopy rebuilds from the wrapped POSCAR; an atom at fractional coordinate 1.0 had its
  periodic images permuted, so the gate evaluated a mislabelled fc2.
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
- Agent calls no longer collect the whole skill for every read or action (on the 328-material
  production set each `inspect`/`cycle` or `act start -p X` took tens of minutes on 9p):
  - `agent snapshot/inspect/plan/propose/cycle/run` and the MCP `get_snapshot`/`inspect`/
    `cycle`/`propose_actions`/`apply_actions` tools take `source` (`auto` default,
    `progress`, `live`). `auto` reads `.tf_progress.json` when it is fresh (live monitor, or
    not older than `narrow_max_age`) and covers the requested scope, and falls back to live
    collection otherwise. Every result says which source it used and how old it is.
  - Execution groups actions by (action, skill, step) into one `act -p A,B,C` call: one
    collection and one shared max_jobs gate per group instead of one full collection per
    action. `start`/`retry` carry `--expect-state`, so each target is re-checked against the
    freshly collected state and skipped if it changed (a plan from the progress file cannot
    resubmit a finished step or `retry` away a running job). Each action reports an
    `outcome`: `submitted`, `deferred` (max_jobs reached), `skipped_stale`, `ok`, `failed`.
  - Commands with explicit `-p` targets (`start`, `retry`, `fetch`, `advance`, `list`, `json`,
    `diagnose`, `conf`, `dir`, `prove`) resolve and probe only the matching materials of a
    shared-layout segment when the progress file recorded a whole-skill collection within
    `narrow_max_age` (default 7200 s; 0 disables). The matcher is a superset of `-p`
    resolution, so ambiguity errors are unchanged. The max_jobs gate adds the running/queued
    jobs of the materials that were not collected, from the progress file plus a submission
    ledger (`.tf_submit_ledger.jsonl`) that covers jobs submitted after the last collection.
    `autozt doctor` shows per skill whether this fast path is available.
- `autozt agent setup` prints a ready-to-paste MCP server entry (absolute interpreter and
  entry paths, `AUTOZT_CONFIG`, `workflow` profile) and a short rules card for the model;
  `--rules` prints only the card.
- Progress views carry `fail_summary` (FAILs bucketed by suggested action with codes and
  ids), and jobs submitted after the last collection show as `PD`
  (`submitted_after_collect`) instead of being proposed again.
- max_jobs was not enforced for explicit `start -p X -j S` (the path every agent
  `start_step` takes), and `start -p A,B,C` built a fresh gate per material from the same
  pre-submission counts, so a batch could overshoot the cap. Both now share one gate; `-f`
  still lets a human exceed the cap on purpose. The gate also counts jobs of materials
  outside the collected scope (`--project`, skill-subdir `-p` pruning) when the progress
  file is fresh.
- `act --project P …` was classified as an unknown (destructive) command because
  `--project`'s value was taken for the command word.

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
