# Changelog

All notable changes to PhonoAgent are documented here.
This project adheres to Semantic Versioning; versions before 1.0.0 are development snapshots.

## [1.1.0] - 2026-09-15

### Added
- MCP interface: phonoagent mcp exposes the tool as a JSON-RPC 2.0 stdio server
  (initialize / tools/list / tools/call). The tool table holds 14 generic verbs
  (list_skills, list_materials, get_summary, get_status, get_json, conf_get, start_step,
  retry_step, fetch_results, conf_set, stop_step, rerun_step, clean_material,
  session_export) and is capped at 20: skills are discovered from skill.yaml, so adding a
  skill never adds a tool.
- MCP policy boundary: read-only tools call the CLI directly; mutating and destructive
  tools are routed through the action gateway (phonoagent act), so an agent session gets
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
- Package layout phonoagent/ with the phonoagent CLI (zero third-party dependency core).
- Criterion-driven step engine: every step declares a physical convergence criterion and
  the engine advances, retries or reports accordingly.
- Provenance per step plus phonoagent prove --verify, and phonoagent session export for
  reproducible archival (manifest with per-file sha256, history timeline, replay hints).
- Action gateway (phonoagent act / approve) with risk tiers, audit log and TTL-bounded
  human approvals for agent-driven sessions.
- Unified core count control: one "cores" setting normalises submit scripts and INCAR
  (NCORE/KPAR) for every skill, whatever template names it uses.
- Unified pre-submission input check: required inputs are discovered from the remote step
  directory instead of per-skill hard-coded lists.
- Supercell specification accepts both diagonal notation ("4 4 4") and a general 3x3
  integer matrix ("2 1 0 -1 2 0 0 0 1", phonopy/phono3py --dim semantics).
- 18 skills covering VASP (band, elastic, defect, ke, kl, opt, phonon), MACE
  (kl/opt/phonon CPU and GPU, mlff-mace), plus te-screen and unihamgnn.
- Drop-in skill interface with phonoagent schema --strict validation.

### Fixed
- vasp_ncl submit template loaded an Intel module environment while executing the
  locally built (AOCC/AOCL) binary, failing with a missing libscalapack; the template now
  uses the matching environment, plus a single-node MPI guard for nodes with a broken
  mlx5/UCX stack.
- defect-dft-cpu shipped ENCUT = auto which VASP rejects (it now uses the documented 370 eV).
- unihamgnn could never submit its graph step (non-VASP inputs were checked against the
  VASP default list).
