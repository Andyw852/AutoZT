# PhonoAgent

Multi-material, multi-step, multi-cluster orchestration for VASP and MACE phonon /
force-constant workflows.

PhonoAgent runs hundreds of materials through multi-step pipelines (relaxation, static,
supercell generation, force-constant fitting, phonon and thermal-conductivity
post-processing) on several HPC clusters, and keeps the whole run inspectable: every
step declares a physical convergence criterion, every input carries a provenance record,
and every failure mode is diagnosed instead of hidden.

## Highlights

- **Criterion-driven engine.** A step advances only when its own physical criterion is
  met; otherwise PhonoAgent retries, escalates or reports, never silently continues.
- **Multi-cluster.** VASP (CPU) and MACE (CPU/GPU) pipelines across SLURM clusters and a
  bare GPU server; switching a material between clusters is a one-word configuration change.
- **Reproducibility.** Per-step provenance (per-file SHA-256) with phonoagent prove --verify,
  and phonoagent session export producing a self-contained, replayable archive.
- **Agent safety.** phonoagent act / approve routes agent-driven actions through risk tiers,
  a TTL-bounded human approval and an audit log.
- **Uniform mechanisms.** One "cores" setting normalises queue scripts and INCAR
  (NCORE/KPAR) for every skill, and required inputs are discovered from the remote step
  directory instead of per-skill hard-coded lists.
- **Supercells.** Diagonal notation ("4 4 4") or a general 3x3 integer matrix
  ("2 1 0 -1 2 0 0 0 1", phonopy/phono3py --dim semantics) for every generation path.
- **Drop-in skills.** A new skill is a directory with a skill.yaml; phonoagent schema
  --strict validates it.

## Skills

18 skills ship today: VASP (band, elastic, defect, electronic thermal conductivity,
lattice thermal conductivity, structure optimisation, phonons), MACE (lattice thermal
conductivity, optimisation, phonons, CPU and GPU variants, MLFF training), plus
te-screen (thermoelectric surrogate screening) and unihamgnn.

## Install

    git clone https://github.com/Andyw852/PhonoAgent
    cd PhonoAgent
    pip install -e .              # the CLI core has no third-party dependency
    pip install -e ".[yaml]"      # + PyYAML, to read yaml configuration
    pip install -e ".[vasp]"      # + the VASP-side scientific stack
    pip install -e ".[mace]"      # + the MACE-side scientific stack

## Quick start (no cluster needed)

    phonoagent -c my.yaml -tt te-screen -p Si_demo start

This runs a two-step, pure-python thermoelectric surrogate screen on a POSCAR and writes
te_features.json plus te_screen_summary.json. See docs/CONFIGURING.md for connecting a
real cluster; setting/template-cluster.yaml is the configuration template.

## Documentation

- README.md - the full manual (Chinese, includes the operational playbook)
- CONTEXT.md - code map for contributors
- docs/CONFIGURING.md - cluster and configuration guide
- CHANGELOG.md - release history

## Citation

See CITATION.cff. If you use PhonoAgent in published work, please cite the software and
the underlying methods (phonopy, phono3py, hiPhive, MACE, spglib).

## License

MIT (see LICENSE). PhonoAgent does not redistribute pseudopotentials, MACE model weights
or third-party training data; users obtain those from their own licensed sources.
