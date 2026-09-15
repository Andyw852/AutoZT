# Evaluation: measured numbers

Every number below was produced by a script in this repository and can be re-run. Paths
are given so a reviewer can reproduce them without reading the source.

## 1. Agent safety (scripts/safety_metrics.py)

    python scripts/safety_metrics.py            # human readable
    python scripts/safety_metrics.py --json     # machine readable (CI asserts on this)
    python scripts/safety_metrics.py --ablation # controlled comparison

| Metric | Value |
|---|---|
| interception rate of hazardous actions | 100% (8/8) |
| hazardous actions actually executed | 0 |
| every denial offers an approval path | yes |
| read-only actions still allowed for agents | yes |
| local output determinism | 100% |

### Ablation (gateway enforced vs gateway switched off)

The experiment runs one hazardous action (clean -y) inside a throwaway project whose
material holds nothing but a POSCAR, so the "gateway off" arm really executes and really
does nothing harmful.

    gateway_strict   rc  blocked  executed  readonly_still_ok
    1                 3    True     False     True
    0                 0    False    True      True
    read-only profile: 6 of 14 tools exposed; hazardous tool call refused = True

Three arms, all measured by the same script:

| Arm | Setting | Behaviour |
|---|---|---|
| A gateway enforced | PHONOAGENT_AGENT_STRICT=1 | hazardous action refused before dispatch (rc=3), nothing executed |
| B gateway switched off | PHONOAGENT_AGENT_STRICT=0 | the very same command executes (only inside the throwaway sandbox) |
| C read-only profile | PHONOAGENT_MCP_READONLY=1 | mutating and destructive verbs are not exposed at all (6 of 14 tools); a direct call to them is refused too |

Reading: arm A is the default for agent sessions; arm B shows what the gateway prevents,
which is the reason approvals are TTY-only and destructive verbs carry their own risk tier;
arm C is the posture for "let the agent look but not touch", enforced one level earlier
than a refusal (the verb is absent from the tool table, and the dispatcher double-checks).

## 2. Reproducibility

| Item | Value | How to reproduce |
|---|---|---|
| provenance verification on a cluster | 17/17 comparable inputs byte-identical | pa -tt <skill> -p <material> prove --verify |
| local output determinism | 100% (skills / schema / mcp --list-tools, repeated) | scripts/safety_metrics.py |
| reproducible archive | per-file sha256 manifest + history timeline + replay hints | pa -p <material> session export |

## 3. Skill coverage (real clusters)

| Scope | Result |
|---|---|
| input generation, 18 skills | 16/18 (cohp-cogito and fc-fit require upstream products, by design) |
| real submissions of the first computable step | 15/15 finished OK |
| VASP 4-core projection | converged for band / elastic / ke / kl / phonon / opt and mlff-mace |
| force-constant fitting (fc-fit, hiphive, Si, 10 frames x 250 atoms) | fc2.hdf5, fc3.hdf5, ShengBTE export written; phonon verdict stable (min_freq ~ 0 THz) |

## 4. Cost of the two operating modes

| Mode | Worth using when |
|---|---|
| rule-driven steps (no LLM in the loop) | hundreds of materials x many steps: deterministic, no token cost, every decision traceable to a physical criterion |
| MCP-driven agent on top | triage of failures, review of generated inputs, decisions about what to retry |

The interface in phonoagent/mcp.py is deliberately generic (14 verbs, capped at 20) plus 37
read-only resources, so a planner can act without the tool itself becoming an agent.

## 5. Threats to validity

- The cluster numbers come from one site (three queues) and a handful of materials; they
  demonstrate that the mechanisms work end to end, not a throughput benchmark.
- The thermoelectric screener (te-screen) is a surrogate with the accuracy reported in its
  README; it is a pre-filter, not a replacement for transport calculations.
- Runtime messages are only partly translated to English; the help text is complete.
- ShengBTE's CONTROL file cannot express a non-diagonal supercell, and that path raises an
  explicit error instead of approximating.

### Read-only profile (added after the first pass)

    tools exposed              6 of 14
    hazardous call refused     True
    read-only task still ok    True   (list_skills used as the representative task)

This is the "let the agent look but not touch" posture: the mutating and destructive verbs
are absent from the tool table, and the dispatcher refuses them even if a client calls a
remembered name, while supervision tasks keep working.
