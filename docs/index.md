# PhonoAgent

Multi-material, multi-step, multi-cluster orchestration for VASP and MACE phonon and
force-constant workflows.

PhonoAgent takes a directory of structures and drives them through multi-step pipelines
(relaxation, static, supercell and displacement generation, force-constant fitting,
phonon and thermal-conductivity post-processing) across several HPC clusters, while
keeping every step inspectable and every input traceable.

## Why it exists

- One command drives hundreds of materials; the engine decides per step whether the
  physical criterion is met, and retries or reports instead of silently continuing.
- VASP (CPU) and MACE (CPU/GPU) pipelines coexist; a material can move between clusters
  with a single configuration change.
- Every step stores provenance (per-file SHA-256) that can be verified later, and a whole
  material can be exported as a reproducible, replayable session archive.
- Agent-driven operation is gated: risk tiers, TTL-bounded human approval and an audit log.
- Mechanisms are uniform across skills: one core-count setting, one pre-submission input
  check, one supercell specification (diagonal or general 3x3 matrix).

## Next steps

- [Installation](installation.md)
- [Quick start](quickstart.md)
- [Configuring clusters](CONFIGURING.md)
- [CLI reference](cli.md)
