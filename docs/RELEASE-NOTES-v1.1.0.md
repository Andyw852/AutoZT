# AutoZT v1.1.0 — release notes (ready to paste into the GitHub Release form)

**Title**: v1.1.0 — MCP interface, safety metrics, cluster-switch self-checks

**Tag**: v1.1.0   (superseded by v1.2.0: see docs/RELEASE-NOTES-v1.2.0.md)

## Highlights

- **MCP interface.** AutoZT now speaks the Model Context Protocol over stdio
  (initialize / tools/list / tools/call) so an LLM planner can drive it directly. The
  tool table holds 14 generic verbs and is capped at 20: skills are discovered from
  skill.yaml, so adding a skill never adds a tool. Read-only tools call the CLI
  directly; mutating and destructive tools are routed through the action gateway, so an
  agent session receives "needs approval" instead of executing, and approvals only work
  from a real TTY.
- **Measured safety and reproducibility.** scripts/safety_metrics.py reports the numbers
  this release is verified against: interception rate of hazardous actions 100% (8/8),
  hazardous actions actually executed 0, every denial offers an approval path, read-only
  actions still allowed, local output determinism 100%. On a real cluster,
  prove --verify reports byte-identical inputs and provenance (17/17 comparable files).
- **Cluster-switch self-checks.** Four traps that used to fail silently are now reported
  before submission: work_dir provenance (status prints which configuration layer
  supplied it), resource requests above the cluster's declared max_cpus, a partition the
  target cluster does not use, and conda environments or conda.sh paths that do not exist
  on the target cluster (that one used to kill a stage without writing any log line).
- **Short command.** Both autozt and pa are installed; they are the same program.
- **English help.** AUTOZT_LANG=en (or --lang en) prints the complete command
  reference in English; the i18n mechanism and an inventory script
  (scripts/i18n_report.py) make the remaining runtime messages incrementally translatable.

## Verification in this release

- Unit and smoke suites: smoke 5/5, MCP 6/6, safety metrics 4/4, preflight checks 6/6,
  i18n 4/4; ten end-to-end self-test suites pass locally and in CI.
- Cluster runs: all 18 skills generated inputs (16/18 without upstream products, as
  designed), 15/15 real submissions finished OK, and the force-constant fitting skill
  (hiphive engine, Si, 10 frames, 250 atoms) produced fc2/fc3 and a stable phonon
  verdict (min_freq ≈ 0 THz).

## Notes

- AutoZT does not redistribute pseudopotentials, MACE model weights or third-party
  training data; obtain them from their licensed sources.
- Configuration for your own clusters: docs/CONFIGURING.md.
