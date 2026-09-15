# Delivery inventory

Where everything lives, so a collaborator can pick the project up without asking.

## Code

| Path | What |
|---|---|
| bin/autozt, bin/pa | entry points (pa is the two-letter alias; same program) |
| autozt/ | the package: bootstrap (discovery/config), cli, collect, data, ops, report, workflow (engine), preflight (pre-submission checks), i18n, mcp (MCP server), agentgate (approval gate + audit), session (reproducible archives), prov (provenance), history, corrections, skillspec, yamlmini |
| skill/ | 18 skills, each a directory with skill.yaml plus generation scripts and templates |
| setting/ | cluster templates (setting/template-cluster.yaml) and per-cluster configuration |
| scripts/ | safety_metrics.py, i18n_report.py, check_site_specific.py and helpers |
| tests/ | pytest suites (smoke, mcp, safety metrics, preflight checks, i18n) plus suite_*.py end-to-end scripts driven by test_suites.py |

## How to run it

    pip install -e .            # CLI only, no third-party dependency
    pa skills                   # list the skills
    pa -c <config> -tt <skill> -p <material> status
    AUTOZT_LANG=en pa --help

## Documentation

| File | Audience |
|---|---|
| README.md | full manual (Chinese) |
| README.en.md | English overview and quick start |
| docs/index.md, installation.md, quickstart.md, skills.md, cli.md | users |
| docs/CONFIGURING.md | site configuration, including the four cluster-switch traps |
| docs/mcp.md | the MCP interface (tools, resources, prompts, risk tiers) |
| docs/I18N.md | how runtime messages are translated, and the safe migration recipe |
| docs/ACCEPTANCE.md | acceptance checklist with evidence per item |
| docs/EVALUATION.md | measured numbers, controlled ablation, threats to validity |
| docs/RELEASING.md, RELEASE-NOTES-v1.1.0.md | release procedure and ready-to-paste notes |
| docs/RELATED-WORK.md, RELATED-WORK-EXTENDED.md | positioning against comparable systems |
| docs/dev/ | internal notes (code map, isolation rules, roadmap, ops memo) |
| CHANGELOG.md, CITATION.cff, LICENSE, .zenodo.json | project metadata |

## Verified state

- 33 assertions across five pytest suites pass locally; the same set runs in CI
  (.github/workflows/ci.yml), with suites that need a cluster tagged and skipped.
- Real-cluster runs: input generation for 18 skills (16 without upstream products, as
  designed), 15/15 real submissions finished, force-constant fitting verified end to end.
- Safety: interception rate of hazardous actions 100% (8/8), hazardous actions executed 0,
  read-only profile exposes 6 of 14 tools while read-only tasks still complete.
- Reproducibility: prove --verify reported 17/17 comparable inputs byte-identical on a
  cluster; local output determinism 100%.

## Repository

- Remote: github.com/Andyw852/AutoZT (private), branch master, tags v1.0.0 and v1.1.0.
- No DOI by choice; a GitHub release is optional and its notes are ready.

## Not done (honest list)

- Runtime message translation is partial (help is complete; see docs/I18N.md for progress).
- ShengBTE CONTROL cannot express a non-diagonal supercell; that path errors explicitly.
- Throughput was not benchmarked; the cluster evidence demonstrates correctness, not speed.
