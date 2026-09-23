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
