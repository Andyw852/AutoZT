# Skills

A skill is a directory under skill/ with a skill.yaml declaring its steps, generation
scripts, convergence criteria and default cluster. AutoZT ships 20 skills:

| Family | Skills |
|---|---|
| VASP (CPU) | band-dft-cpu, defect-dft-cpu, elastic-dft-cpu, ke-dft-cpu, kl-dft-cpu, opt-dft-cpu, phonon-dft-cpu, zt-dft-cpu |
| MACE | kl-mlff-cpu, kl-mlff-gpu, opt-mlff-cpu, opt-mlff-gpu, phonon-mlff-cpu, phonon-mlff-gpu, mlff |
| Auxiliary | cohp-cogito (COHP/ICOHP bonding analysis), eph-qe-cpu (Quantum ESPRESSO electron-phonon), te-screen (thermoelectric surrogate screening), unihamgnn (graph data generation) |
| Fitting | fc-fit (phono3py / pheasy / hiphive force-constant fitting) |

`zt-dft-cpu` is a composite skill: it declares the electron-transport steps of `ke-dft-cpu`
and the lattice-thermal-conductivity steps of `kl-dft-cpu` by **symlinking** their step and
template directories (no generator is copied), and adds one new generator, `S20_zt`, which
computes ZT = S²σT/(κ_e + κ_L) over the full temperature × doping grid. Either segment can be
switched off (`electronic: false` / `lattice: false`) so the summary step can reuse results
produced by the standalone skills. Details: `skill/zt-dft-cpu/README.md`.

Skills that fit the existing generator/checker/step lifecycle can be added as a directory
without changing the MCP tool table. `autozt schema --strict` validates declaration
structure and some references; it does not prove complete scientific I/O coverage.

## Developing skills for LLM access

The maintained development specification is **TASKFLOW.md, chapter 7** in the software
package. Sections **7.15–7.18** cover model-facing contracts, when core/interface changes
are required, failure reporting, and validation against actual generators and outputs.
Chapter **9.1** explains CLI/MCP entrypoints, process lifetime and generated records.

- Declare `schema: 2`, `io_schema`, `flow` and applicable `corrections`; keep them aligned
  with generators, templates, `CONF_SPEC`, checkers and the skill README.
- New file formats and parameters normally stay inside the skill. New execution semantics
  need core support first, followed by shared Agent CLI/MCP changes when necessary.
- State reporting is not automatic scientific-result parsing. Output declarations do not
  automatically fetch, convert or chain artifacts.
- Test missing inputs, invalid parameters, failed/partial results and recovery behavior;
  distinguish declaration checks from authorized real-cluster verification.

See [Agent CLI](agent-cli.md) and [MCP](mcp.md) for transport details. `agent serve` is an
optional JSONL integration endpoint, not a prerequisite for ordinary use or monitoring.
