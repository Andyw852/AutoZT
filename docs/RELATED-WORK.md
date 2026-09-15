# Related work: PhonoAgent against agentic and workflow systems

Sources read locally (extracted text): AtomAgents (PNAS 2025, 10.1073/pnas.2414074122,
12 pages) and DREAMS (arXiv:2507.14267v2, 91 pages). Line numbers below refer to the
extracted text used for this check.

## 0. Corrections to earlier drafts of this file (important)

My first draft of this comparison asserted five things about these systems that the
primary sources do NOT support. They are retracted here:

| Earlier claim | Verdict from the papers |
|---|---|
| "DREAMS runs on a single cluster / locally" | WRONG. DREAMS has an "HPC agent [that] schedules, monitors, and retrieves simulations" on an "HPC cluster" (p. ~289, 543). |
| "DREAMS has no safety guard / approval before acting" | WRONG. DREAMS states it is "built around a multi-tier safety guard" applying "deterministic" validation (abstract, line 19; section "Safety Guards", line 253). |
| "DREAMS has no provenance tracking" | WRONG. DREAMS reports that "a judge audits the full provenance graph behind every claim" (line 24) and discusses "circumvention of provenance checks" as a failure mode (line 121). |
| "Deterministic replay is not a concern for these systems" | WRONG for DREAMS: it "applies deterministic validation when explicit ..." (line 149), "deterministically verify values and their provenance" (line 284), and argues that autonomous systems "must ... reproduce both the execution and ..." (line 50). |
| "AtomAgents: LLM holds tools directly, no gate" | PARTLY WRONG. AtomAgents routes plans through a "critic" agent: the plan "is then evaluated and approved by the accompanying critic agent" (line 429). That is an agent-internal approval, not a human authorisation, but it is not "no gate". |

Also note: the AtomAgents paper (12 pages) contains no occurrence of cluster, HPC,
reproducibility, determinism, token/API cost or provenance. The correct statement is
"not discussed in that paper", NOT "does not exist".

## 1. What the two systems actually are

- AtomAgents: physics-aware multimodal multi-agent AI for alloy design and discovery
  (planner + critic agents, LLM-driven, atomistic simulation tools).
- DREAMS: a "Density Functional Theory Based Research Engine for Agentic Materials
  Simulation" with a multi-tier safety guard, a judge that audits a provenance graph,
  an HPC agent that schedules simulations, and an explicit determinism argument.

## 2. Defensible differences (rewritten after reading the sources)

| Dimension | DREAMS / AtomAgents | PhonoAgent | Basis |
|---|---|---|---|
| Where the loop decision lives | LLM agents decide and a critic/judge validates the decision | Physical criteria decide (converged / imaginary frequency / marker); the steady-state loop contains no LLM call | own repository (verifiable) vs DREAMS/AtomAgents architecture |
| Nature of the safety guard | Guard lives inside the agent stack (validation of values, provenance checks, judge) | Guard lives at the tool boundary: risk tiers, TTY-only time-bounded human approval, audit log; agent never receives an ungated destructive command | DREAMS line 19/253 vs PhonoAgent gateway |
| Provenance granularity | Claim-level provenance graph audited by a judge | Artifact-level: per-file sha256 per step, verifiable archive, replayable session export | DREAMS line 24 vs PhonoAgent prov/session |
| Determinism target | Deterministic validation of values and parameters inside the agent's reasoning | Deterministic execution: same input, same artifacts, hash-checkable reruns | DREAMS line 149/284 vs PhonoAgent prove --verify |
| Infrastructure model | One HPC cluster, scheduler managed by an HPC agent | Heterogeneous: multiple real SLURM clusters plus a scheduler-less GPU host, per-material cluster choice, mixed DFT and MLIP steps | DREAMS line 289 vs PhonoAgent config/status |
| Cost model | Token cost per agent decision (reported by DREAMS) | Zero LLM calls in steady state; LLM only on failure paths | DREAMS line token counts vs PhonoAgent rule engine |
| Scope | Open-ended research and discovery campaigns | Production throughput on known pipelines, with explicit refusal when a capability is unsupported | both sources |

## 3. Honest positioning sentence

PhonoAgent does not claim to be the first system with a safety guard, provenance or HPC
awareness; DREAMS already combines those inside an agentic research engine. PhonoAgent
differs in WHERE those guarantees live and WHAT they cover: its guarantees are attached
to the execution substrate (criteria, per-file artifacts, action gate, multi-cluster
scheduler) rather than to the reasoning stack, they are measurable (interception rate,
hash agreement), and they remain in force when the planning layer is replaced.

## 4. Experiments that would make this comparison empirical

1. Scale-cost: LLM tokens and wall clock versus number of materials (flat vs growing).
2. Safety ablation: destructive-action interception rate and operator errors, gated vs
   direct tool access. Compare against DREAMS' multi-tier guard if it can be run.
3. Reproducibility: artifact hash agreement over N reruns of identical input.
4. Cross-code consistency: phonon spectra and kappa for Si/MgO/PbTe.

## 5. Open items before submission

- Read the DREAMS safety-guard section in full (what each tier gates) and quote it.
- Read the AtomAgents critic/planner workflow in full and state its approval granularity.
- Check DREAMS reproducibility claims for bit-exactness of artifacts (vs claim-level).
- Decide whether the comparison table cites page numbers or section names in the paper.
