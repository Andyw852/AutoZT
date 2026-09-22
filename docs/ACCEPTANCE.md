# Acceptance checklist (AutoZT 1.2.0)

Everything below was produced and verified inside this repository. Each row names the
artifact that carries the evidence, so a reviewer can re-run it.

## 1. Interfaces

| Item | Evidence |
|---|---|
| CLI works under both names (autozt, pa) | bin/autozt, bin/pa; pyproject scripts |
| English help | AUTOZT_LANG=en pa --help (50 lines, 29 commands, no Chinese left) |
| MCP server (JSON-RPC over stdio) | autozt mcp; tests/test_mcp.py (32 checks, including profiles, capabilities, bounded proposals, configured skill resources, failure data, and cursor snapshots) |
| 24 generic MCP tools in the explicit full profile; compact exposes 7, workflow 12, monitor 2 | autozt mcp --list-tools; MCP tests assert skill-agnostic names and profile output |

## 2. Safety

| Metric | Value | Evidence |
|---|---|---|
| interception rate of hazardous actions | 100% (9/9) | scripts/safety_metrics.py |
| hazardous actions actually executed | 0 | scripts/safety_metrics.py |
| every denial offers an approval path | yes | scripts/safety_metrics.py |
| read-only actions still allowed for agents | yes | scripts/safety_metrics.py |
| approval only from a real TTY | enforced by the gateway | tests/suite_agentgate.py |

## 3. Reproducibility

| Item | Value | Evidence |
|---|---|---|
| provenance verification on a real cluster | 17/17 comparable inputs byte-identical | pa ... prove --verify on a completed material |
| local output determinism | 100% | scripts/safety_metrics.py |
| reproducible archive | session export (manifest + per-file sha256 + replay hints) | autozt/session.py; tests/suite_session.py |

## 4. Skill coverage (real runs)

| Scope | Result |
|---|---|
| input generation, all 18 skills | 16/18 (cohp-cogito and fc-fit need upstream products, as designed) |
| real submissions | 15/15 finished OK (VASP 4-core projections converged, MACE relaxations finished) |
| force-constant fitting (fc-fit, hiphive, Si, 10 frames x 250 atoms) | fc2/fc3 + ShengBTE export written; phonon verdict stable (min_freq ~ 0 THz) |

## 5. Cluster-switch traps (all fixed and self-checked)

1. work_dir provenance is printed by status (project setting.yaml outranks hpc.yaml).
2. resource requests above the cluster's declared max_cpus are flagged.
3. a partition the target cluster does not use is reported with the partitions its own
   templates use.
4. conda environments or conda.sh paths missing on the target cluster are flagged before
   submission (this was the cause of a stage that died without a single log line).

## 6. Test inventory (all green locally; the same set runs in CI)

    tests/test_smoke.py            5 assertions
    tests/test_mcp.py              31
    tests/test_safety_metrics.py   4
    tests/test_preflight_checks.py 6
    tests/test_i18n.py             6
    tests/test_suites.py           run via `python3 -m pytest tests/test_suites.py`; drives the
                                  local end-to-end suites (cluster ones tagged). Running the
                                  file directly with `python3` is a NO-OP.

## 7. Measured numbers

See docs/EVALUATION.md for the scripts, the table of measured values and the
controlled ablation, plus the threats-to-validity section.

## 8. Known gaps

- Runtime messages are only partly translated: help is 100% English, but the body of
  status/ops messages remains Chinese; scripts/i18n_report.py reports the lower bound.
- ShengBTE's CONTROL file cannot express a non-diagonal supercell; that path reports an
  explicit error instead of approximating.
- The repository is private and has no DOI (a deliberate choice); a GitHub release is
  optional and the notes are ready in docs/RELEASE-NOTES-v1.1.0.md.
